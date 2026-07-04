#!/usr/bin/env python
"""Interactively label coincident dips ('|') and flares ('v') in the DBSP
Hα-EW + photometry mosaics.

Opens one interactive matplotlib window per figure.  Left-click adds a dip,
right-click adds a flare, clicking an existing marker removes it; labels are
auto-saved to output/coincidence_labels.csv, which the plotting drivers read.

Usage:
    python run_label.py                 # label both figures in turn
    python run_label.py --figure mosaic # main mosaic only
    python run_label.py --figure supp   # supplement only

env: cpv  (/Users/luke/local/miniconda3/envs/cpv/bin/python)
"""
from __future__ import annotations

import argparse

from label_coincidences.interactive import LabelSession
from label_coincidences.labels import LABELS_PATH


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figure", choices=["mosaic", "supp", "both"],
                    default="both", help="which figure(s) to label")
    args = ap.parse_args()

    figs = ["mosaic", "supp"] if args.figure == "both" else [args.figure]
    for fig in figs:
        print(f"Labeling '{fig}' — close the window (or press 'q') to continue.")
        LabelSession(fig).run()

    print(f"Labels written to {LABELS_PATH}")


if __name__ == "__main__":
    main()
