"""Refresh TESS sector coverage columns for the CPV concat table."""

from __future__ import annotations

import argparse
import logging
import shutil
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from astropy.time import Time

ROOT_DIR = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = ROOT_DIR / (
    "results/tables/jan2026_compilation/"
    "concat_R16_S17_S18_B20_S21_Z19_G22_P23_B24_qlp_0to100pc.csv"
)
DEFAULT_SEPARATOR = "|"
DEFAULT_TIC_COL = "ticid"
DEFAULT_RA_COL = "tic8_ra"
DEFAULT_DEC_COL = "tic8_dec"
DEFAULT_LATEST_NAME = "concat_R16_S17_S18_B20_S21_Z19_G22_P23_B24_qlp_0to100pc_latest.csv"

LOGGER = logging.getLogger(__name__)


def _setup_logging(verbose: bool) -> None:
    """Configure logging for the script.

    Args:
        verbose: When True, enable DEBUG logging.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")


def _current_sector_number() -> int | None:
    """Return the current TESS sector number based on midtimes.

    Returns:
        The latest sector whose midpoint is not in the future, or None if
        unavailable.
    """
    try:
        from tess_stars2px import TESS_Spacecraft_Pointing_Data as tspd
    except ImportError as exc:
        raise ImportError(
            "Missing dependency 'tess_stars2px'. Install it in the cpv environment."
        ) from exc

    sectors = np.array(tspd.sectors, dtype=int)
    times = Time(tspd.midtimes, format="jd")
    now = Time.now()

    in_past = times <= now
    if not np.any(in_past):
        return None

    return int(np.max(sectors[in_past]))


def _normalize_sector_value(value: object, max_sector: int | None = None) -> tuple[str, int]:
    """Normalize a tess-point sector field into a string and count.

    Args:
        value: Raw sector value returned by tess-point.
        max_sector: Optional maximum sector number to include.

    Returns:
        Tuple of (comma-separated sector string, number of sectors).
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "", 0

    if isinstance(value, str):
        tokens = [token.strip() for token in value.split(",") if token.strip()]
    elif isinstance(value, (list, tuple, np.ndarray, pd.Series)):
        tokens = [str(token).strip() for token in value if str(token).strip()]
    else:
        tokens = [str(value).strip()]

    numeric_tokens: list[int] = []
    other_tokens: list[str] = []
    for token in tokens:
        try:
            numeric_tokens.append(int(float(token)))
        except (TypeError, ValueError):
            other_tokens.append(token)

    if max_sector is not None:
        numeric_tokens = [token for token in numeric_tokens if token <= max_sector]

    numeric_unique = sorted(set(numeric_tokens))
    other_unique = sorted(set(other_tokens))

    combined = [str(token) for token in numeric_unique] + other_unique
    return ",".join(combined), len(combined)


def _resolve_column(df: pd.DataFrame, preferred: str, fallback: str | None = None) -> str:
    """Resolve a column name with an optional fallback.

    Args:
        df: DataFrame containing columns.
        preferred: Preferred column name.
        fallback: Optional fallback column name.

    Returns:
        The resolved column name.

    Raises:
        KeyError: If neither column exists.
    """
    if preferred in df.columns:
        return preferred
    if fallback and fallback in df.columns:
        return fallback
    raise KeyError(f"Missing required column: {preferred}")


def _fetch_sectors(
    ra: float, dec: float, ticid: object, max_sector: int | None
) -> tuple[str, int]:
    """Query tess-point for sector coverage.

    Args:
        ra: Right ascension in degrees.
        dec: Declination in degrees.
        ticid: TIC identifier.

    Returns:
        Tuple of (sectors string, number of sectors).

    Raises:
        RuntimeError: If tess-point results are missing the sector column.
    """
    try:
        from tars.tess import check_tesspoint
    except ImportError as exc:
        raise ImportError(
            "Missing dependency 'tars'. Install it in the cpv environment before running."
        ) from exc

    result = check_tesspoint(ra, dec, ticid, return_tuple=False)
    if "sector" not in result.columns:
        raise RuntimeError("tess-point result missing 'sector' column")
    return _normalize_sector_value(result["sector"].iloc[0], max_sector=max_sector)


def _update_sectors(
    df: pd.DataFrame,
    ra_col: str,
    dec_col: str,
    tic_col: str,
    fallback_sectors: Iterable[object],
    fallback_counts: Iterable[object],
    current_sector: int | None,
) -> tuple[list[str], list[int]]:
    """Update sector values using tess-point.

    Args:
        df: Input DataFrame.
        ra_col: Column containing right ascension values.
        dec_col: Column containing declination values.
        tic_col: Column containing TIC IDs.
        fallback_sectors: Existing sector values.
        fallback_counts: Existing N_sectors values.

    Returns:
        Lists of updated sector strings and counts.
    """
    cache: dict[str, tuple[str, int]] = {}
    sectors_list: list[str] = []
    counts_list: list[int] = []

    for idx, row in df.reset_index(drop=True).iterrows():
        ra = row.get(ra_col)
        dec = row.get(dec_col)
        ticid = row.get(tic_col)
        fallback_sector = fallback_sectors[idx]
        fallback_count = fallback_counts[idx]

        if pd.isna(ra) or pd.isna(dec) or pd.isna(ticid):
            LOGGER.warning(
                f"Row {idx} missing RA/Dec/TICID; keeping existing sectors value."
            )
            sectors_list.append(str(fallback_sector) if pd.notna(fallback_sector) else "")
            counts_list.append(int(fallback_count) if pd.notna(fallback_count) else 0)
            continue

        cache_key = str(ticid)
        if cache_key in cache:
            sectors, n_sectors = cache[cache_key]
        else:
            try:
                sectors, n_sectors = _fetch_sectors(
                    float(ra), float(dec), ticid, max_sector=current_sector
                )
            except Exception as exc:  # pylint: disable=broad-except
                LOGGER.warning(
                    f"Row {idx} TIC {ticid}: {exc}; keeping existing sectors value."
                )
                sectors = str(fallback_sector) if pd.notna(fallback_sector) else ""
                n_sectors = int(fallback_count) if pd.notna(fallback_count) else 0
            cache[cache_key] = (sectors, n_sectors)

        sectors_list.append(sectors)
        counts_list.append(n_sectors)

        if (idx + 1) % 25 == 0:
            LOGGER.info(f"Processed {idx + 1} rows")

    return sectors_list, counts_list


def refresh_table(
    input_path: Path,
    output_path: Path,
    sep: str,
    ra_col: str,
    dec_col: str,
    tic_col: str,
    backup: bool,
    latest_name: str,
) -> Path:
    """Refresh the sectors and N_sectors columns in the input table.

    Args:
        input_path: Path to the input CSV table.
        output_path: Path to write the updated CSV table.
        sep: CSV delimiter.
        ra_col: Right ascension column name.
        dec_col: Declination column name.
        tic_col: TIC ID column name.
        backup: When True, write a timestamped backup of the input table.
    """
    df = pd.read_csv(input_path, sep=sep)

    resolved_ra = _resolve_column(df, ra_col, fallback="tic8_RA_orig")
    resolved_dec = _resolve_column(df, dec_col, fallback="tic8_Dec_orig")
    resolved_tic = _resolve_column(df, tic_col, fallback="ticid_x")

    fallback_sectors = df.get("sectors", pd.Series([""] * len(df))).tolist()
    fallback_counts = df.get("N_sectors", pd.Series([0] * len(df))).tolist()

    current_sector = _current_sector_number()
    if current_sector is None:
        LOGGER.warning("Could not determine current sector; keeping all sectors.")
    else:
        LOGGER.info(f"Masking future sectors beyond S{current_sector}.")

    LOGGER.info(
        f"Refreshing sectors for {len(df)} rows using RA='{resolved_ra}', "
        f"Dec='{resolved_dec}', TICID='{resolved_tic}'."
    )

    updated_sectors, updated_counts = _update_sectors(
        df,
        resolved_ra,
        resolved_dec,
        resolved_tic,
        fallback_sectors,
        fallback_counts,
        current_sector,
    )

    if backup:
        timestamp = datetime.utcnow().strftime("%Y%m%d")
        backup_path = input_path.with_suffix(f".bak_{timestamp}.csv")
        shutil.copy2(input_path, backup_path)
        LOGGER.info(f"Wrote backup to {backup_path}")

    df["sectors"] = updated_sectors
    df["N_sectors"] = updated_counts

    df.to_csv(output_path, sep=sep, index=False)
    LOGGER.info(f"Wrote updated table to {output_path}")

    timestamp = datetime.utcnow().strftime("%Y%m%d")
    dated_path = output_path.with_suffix(f".{timestamp}.csv")
    shutil.copy2(output_path, dated_path)
    LOGGER.info(f"Wrote dated snapshot to {dated_path}")

    latest_path = output_path.with_name(latest_name)
    shutil.copy2(output_path, latest_path)
    LOGGER.info(f"Wrote latest snapshot to {latest_path}")

    return output_path


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Refresh the sectors and N_sectors columns using tess-point via tars.tess."
        )
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--sep", default=DEFAULT_SEPARATOR)
    parser.add_argument("--ra-col", default=DEFAULT_RA_COL)
    parser.add_argument("--dec-col", default=DEFAULT_DEC_COL)
    parser.add_argument("--tic-col", default=DEFAULT_TIC_COL)
    parser.add_argument("--latest-name", default=DEFAULT_LATEST_NAME)
    parser.add_argument("--backup", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Entry point for the script."""
    args = _parse_args()
    _setup_logging(args.verbose)

    input_path = args.input
    output_path = args.output or args.input

    if not input_path.exists():
        raise FileNotFoundError(f"Input table not found: {input_path}")

    refresh_table(
        input_path=input_path,
        output_path=output_path,
        sep=args.sep,
        ra_col=args.ra_col,
        dec_col=args.dec_col,
        tic_col=args.tic_col,
        backup=args.backup,
        latest_name=args.latest_name,
    )


if __name__ == "__main__":
    main()
