"""Build a manifest of CPV vetter PDFs for cpv-viz."""

from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

ROOT_DIR = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT = ROOT_DIR / "apps/cpv-viz/data/cpv_pdf_manifest.json"
DEFAULT_BASE_URL = "./pdfs"
DEFAULT_ROOT = Path.home() / "local" / "complexrotators"
PRIMARY_PDF_DIR = DEFAULT_ROOT / "tars_jan2026_knownCPVs"
FALLBACK_PDF_DIR = DEFAULT_ROOT / "cpv_finding" / "tars_jan2026_knownCPVs"
DEFAULT_PATTERN = "*_cpvvetter.pdf"


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments.

    Returns:
        Parsed arguments namespace.
    """
    parser = argparse.ArgumentParser(
        description="Build cpv-viz PDF manifest from CPV vetter outputs.",
    )
    parser.add_argument(
        "--pdf-dir",
        default=None,
        help=(
            "Directory containing the *_cpvvetter.pdf files. Defaults to "
            f"{PRIMARY_PDF_DIR} (fallback {FALLBACK_PDF_DIR})."
        ),
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT),
        help="Output JSON path for the PDF manifest.",
    )
    parser.add_argument(
        "--base-url",
        default=DEFAULT_BASE_URL,
        help=(
            "Base URL prefix used by the web app when linking PDFs. "
            "Use a relative path like ./pdfs after creating a symlink."
        ),
    )
    parser.add_argument(
        "--pattern",
        default=DEFAULT_PATTERN,
        help="Glob pattern to match PDF files.",
    )
    return parser.parse_args()


def _resolve_pdf_dir(raw_path: str | None) -> Path:
    """Resolve the PDF directory, applying fallback defaults.

    Args:
        raw_path: Optional directory supplied on the CLI.

    Returns:
        Resolved PDF directory path.

    Raises:
        FileNotFoundError: If no valid PDF directory is found.
    """
    if raw_path:
        candidate = Path(raw_path).expanduser()
        if candidate.exists():
            return candidate
        raise FileNotFoundError(f"PDF directory not found: {candidate}")

    if PRIMARY_PDF_DIR.exists():
        return PRIMARY_PDF_DIR
    if FALLBACK_PDF_DIR.exists():
        return FALLBACK_PDF_DIR
    raise FileNotFoundError(
        "No default PDF directory found. "
        f"Tried {PRIMARY_PDF_DIR} and {FALLBACK_PDF_DIR}."
    )


def _extract_ticid(filename: str) -> str | None:
    """Extract the TIC ID from a PDF filename.

    Args:
        filename: PDF filename to parse.

    Returns:
        TIC ID string, or None if parsing fails.
    """
    token = filename.split("_")[0]
    digits = "".join(ch for ch in token if ch.isdigit())
    return digits if digits else None


def _extract_sector_id(filename: str) -> int | None:
    """Extract the sector ID from a PDF filename.

    Args:
        filename: PDF filename to parse.

    Returns:
        Sector number, or None if no sector is found.
    """
    parts = filename.split("_")
    if len(parts) < 2:
        return None
    digits = "".join(ch for ch in parts[1] if ch.isdigit())
    if not digits:
        return None
    return int(digits)


def _sort_key(filename: str) -> Tuple[int, int, str]:
    """Build a sort key for PDF filenames.

    Args:
        filename: PDF filename to sort.

    Returns:
        Tuple suitable for sorting by sector then filename.
    """
    sector = _extract_sector_id(filename)
    if sector is None:
        return (1, 0, filename)
    return (0, sector, filename)


def _collect_pdfs(pdf_dir: Path, pattern: str) -> Dict[str, List[str]]:
    """Collect PDF filenames grouped by TIC ID.

    Args:
        pdf_dir: Directory containing the PDF files.
        pattern: Glob pattern for matching PDFs.

    Returns:
        Mapping of TIC ID to sorted list of PDF filenames.
    """
    grouped: Dict[str, List[str]] = {}
    for pdf_path in pdf_dir.glob(pattern):
        if not pdf_path.is_file():
            continue
        ticid = _extract_ticid(pdf_path.name)
        if not ticid:
            continue
        grouped.setdefault(ticid, []).append(pdf_path.name)

    for ticid, pdfs in grouped.items():
        pdfs.sort(key=_sort_key)
        grouped[ticid] = pdfs

    return grouped


def build_manifest(pdf_dir: Path, output_path: Path, base_url: str, pattern: str) -> Path:
    """Build and write the PDF manifest JSON.

    Args:
        pdf_dir: Directory containing PDF files.
        output_path: Path to write the manifest JSON.
        base_url: Base URL prefix for linking to PDFs.
        pattern: Glob pattern for matching PDF files.

    Returns:
        Path to the written manifest file.
    """
    pdfs_by_tic = _collect_pdfs(pdf_dir, pattern)
    manifest = {
        "generated_date": datetime.now().strftime("%Y-%m-%d"),
        "pdf_base_url": base_url,
        "source_dir": str(pdf_dir),
        "tic_count": len(pdfs_by_tic),
        "pdf_count": sum(len(pdfs) for pdfs in pdfs_by_tic.values()),
        "pdfs": pdfs_by_tic,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return output_path


def main() -> None:
    """Run the manifest build CLI."""
    args = _parse_args()
    pdf_dir = _resolve_pdf_dir(args.pdf_dir)
    output_path = Path(args.output)
    manifest_path = build_manifest(pdf_dir, output_path, args.base_url, args.pattern)
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
