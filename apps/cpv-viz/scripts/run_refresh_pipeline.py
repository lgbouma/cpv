"""Refresh TESS sectors and rebuild the cpv-viz data asset."""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
from datetime import datetime, timezone
from pathlib import Path


LOGGER = logging.getLogger(__name__)
ROOT_DIR = Path(__file__).resolve().parents[3]
DEFAULT_STATUS_PATH = ROOT_DIR / "apps/cpv-viz/data/refresh_status.json"
DEFAULT_LIT_INPUT = ROOT_DIR / (
    "results/tables/jan2026_compilation/"
    "20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24_TIC8_obs_supplemented.csv"
)
DEFAULT_LIT_LATEST_NAME = (
    "20240304_CPV_lit_compilation_R16_S17_S18_B20_S21_Z19_G22_P23_B24_TIC8_obs_supplemented_latest.csv"
)


def _utc_now() -> str:
    """Return a UTC timestamp in ISO-8601 format."""
    return datetime.now(timezone.utc).isoformat()


def _write_status(status_path: Path, payload: dict) -> None:
    """Write a JSON status payload to disk.

    Args:
        status_path: Output path for the status file.
        payload: Dictionary payload to serialize.
    """
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _setup_logging(verbose: bool) -> None:
    """Configure logging for the script.

    Args:
        verbose: When True, enable DEBUG logging.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")


def _run_command(command: list[str]) -> None:
    """Run a subprocess command and raise on failure.

    Args:
        command: Command list suitable for subprocess.

    Raises:
        subprocess.CalledProcessError: If the command fails.
    """
    LOGGER.info("Running: %s", " ".join(command))
    subprocess.run(command, check=True)


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run the weekly cpv-viz refresh pipeline.",
    )
    parser.add_argument(
        "--python",
        default="/Users/luke/local/miniconda3/envs/cpv/bin/python",
        help="Python executable to use for pipeline steps.",
    )
    parser.add_argument(
        "--status-path",
        default=str(DEFAULT_STATUS_PATH),
        help="Path to write the pipeline status JSON.",
    )
    parser.add_argument(
        "--skip-lit-update",
        action="store_true",
        help="Skip updating the literature compilation table.",
    )
    parser.add_argument(
        "--lit-input",
        default=str(DEFAULT_LIT_INPUT),
        help="Literature compilation CSV to update.",
    )
    parser.add_argument(
        "--lit-sep",
        default=",",
        help="Delimiter for the literature compilation CSV.",
    )
    parser.add_argument(
        "--lit-latest-name",
        default=DEFAULT_LIT_LATEST_NAME,
        help="Filename for the literature compilation latest snapshot.",
    )
    parser.add_argument(
        "--lit-ra-col",
        default="tic8_ra",
        help="RA column for the literature compilation table.",
    )
    parser.add_argument(
        "--lit-dec-col",
        default="tic8_dec",
        help="Dec column for the literature compilation table.",
    )
    parser.add_argument(
        "--lit-tic-col",
        default="ticid",
        help="TIC column for the literature compilation table.",
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        help="Create a backup of the input concat table.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def _resolve_status_path(raw_path: str) -> Path:
    """Resolve the status path relative to the repo root if needed.

    Args:
        raw_path: Raw path string supplied on the CLI.

    Returns:
        Absolute path anchored at the repository root.
    """
    status_path = Path(raw_path)
    if status_path.is_absolute():
        return status_path
    return ROOT_DIR / status_path


def main() -> None:
    """Run the refresh pipeline."""
    args = _parse_args()
    _setup_logging(args.verbose)

    python_exe = Path(args.python)
    if not python_exe.exists():
        raise FileNotFoundError(f"Python executable not found: {python_exe}")

    status_path = _resolve_status_path(args.status_path)
    status_payload = {
        "started_at": _utc_now(),
        "finished_at": None,
        "status": "running",
        "steps": [],
    }
    _write_status(status_path, status_payload)

    try:
        update_cmd = [
            str(python_exe),
            str(ROOT_DIR / "apps/cpv-viz/scripts/update_tess_sectors.py"),
        ]
        if args.backup:
            update_cmd.append("--backup")

        step_update = {
            "name": "update_tess_sectors",
            "command": " ".join(update_cmd),
            "started_at": _utc_now(),
            "finished_at": None,
            "status": "running",
            "error": None,
        }
        status_payload["steps"].append(step_update)
        _write_status(status_path, status_payload)
        _run_command(update_cmd)
        step_update["finished_at"] = _utc_now()
        step_update["status"] = "success"
        _write_status(status_path, status_payload)

        if not args.skip_lit_update:
            lit_input = Path(args.lit_input)
            lit_cmd = [
                str(python_exe),
                str(ROOT_DIR / "apps/cpv-viz/scripts/update_tess_sectors.py"),
                "--input",
                str(lit_input),
                "--sep",
                args.lit_sep,
                "--latest-name",
                args.lit_latest_name,
                "--ra-col",
                args.lit_ra_col,
                "--dec-col",
                args.lit_dec_col,
                "--tic-col",
                args.lit_tic_col,
            ]
            if args.backup:
                lit_cmd.append("--backup")

            step_lit = {
                "name": "update_lit_compilation",
                "command": " ".join(lit_cmd),
                "started_at": _utc_now(),
                "finished_at": None,
                "status": "running",
                "error": None,
            }
            status_payload["steps"].append(step_lit)
            _write_status(status_path, status_payload)
            _run_command(lit_cmd)
            step_lit["finished_at"] = _utc_now()
            step_lit["status"] = "success"
            _write_status(status_path, status_payload)

        build_cmd = [
            str(python_exe),
            str(ROOT_DIR / "apps/cpv-viz/scripts/build_data.py"),
        ]
        step_build = {
            "name": "build_data",
            "command": " ".join(build_cmd),
            "started_at": _utc_now(),
            "finished_at": None,
            "status": "running",
            "error": None,
        }
        status_payload["steps"].append(step_build)
        _write_status(status_path, status_payload)
        _run_command(build_cmd)
        step_build["finished_at"] = _utc_now()
        step_build["status"] = "success"
        status_payload["status"] = "success"
        status_payload["finished_at"] = _utc_now()
        _write_status(status_path, status_payload)

    except subprocess.CalledProcessError as exc:
        status_payload["status"] = "failed"
        status_payload["finished_at"] = _utc_now()
        if status_payload["steps"]:
            status_payload["steps"][-1]["finished_at"] = _utc_now()
            status_payload["steps"][-1]["status"] = "failed"
            status_payload["steps"][-1]["error"] = str(exc)
        _write_status(status_path, status_payload)
        raise
    except Exception as exc:  # pylint: disable=broad-except
        status_payload["status"] = "failed"
        status_payload["finished_at"] = _utc_now()
        if status_payload["steps"]:
            status_payload["steps"][-1]["finished_at"] = _utc_now()
            status_payload["steps"][-1]["status"] = "failed"
            status_payload["steps"][-1]["error"] = str(exc)
        _write_status(status_path, status_payload)
        raise


if __name__ == "__main__":
    main()
