#!/usr/bin/env bash
#
# Wed Oct 15 11:58:59 2025
# USAGE: From within /cpv/data/spectra/MIKE/RAW/2025XXYY/, launch CarPy
# installation and update the science indices (and target info) below.
#

set -euo pipefail

DB="mike_20250102_RDXMIKE.db"
DB_TEMPLATE="mike_20250102_RDXMIKE_TEMPLATE.db"

# Cleanup function
cleanup() {
  rm targ* || true
  rm -f -- *.fits || true
  rm -f -- Makefile* || true
  rm -rf -- flatblue slitblue lampblue flatred slitred lampred || true
  rm -rf -- Final-Products || true
}

mkdir -p Final-Science-Products

# Loop over science indices: 0079, 0081, ..., 0117
for i in $(seq 79 2 117); do
  SCIENCE_INDEX=$(printf "%04d" "$i")
  THAR_INDEX=$(printf "%04d" $((i-1)))

  echo "$(date '+%Y-%m-%d %H:%M:%S') Starting ${SCIENCE_INDEX}..."

  # Build a fresh DB from template, uncommenting ONLY matching FITS lines
  patterns="b${SCIENCE_INDEX}\.fits|r${SCIENCE_INDEX}\.fits|b${THAR_INDEX}\.fits|r${THAR_INDEX}\.fits"
  if [[ ! -f "$DB_TEMPLATE" ]]; then
    echo "Template file '$DB_TEMPLATE' not found" >&2
    exit 1
  fi
  perl -pe "if (/$patterns/) { s/^\\s*[#;!]+\\s*// }" "$DB_TEMPLATE" > "$DB"

  # Run cleanup for this iteration
  cleanup

  # Setup and build
  mikesetup -db "$DB" -red -blue -all -mk Makefile
  make

  # Save final products with the science index
  if [[ -f "tic3006blue/tic3006blue_multi.fits" ]]; then
    cp "tic3006blue/tic3006blue_multi.fits" "Final-Science-Products/tic3006blue_multi_${SCIENCE_INDEX}.fits"
  else
    echo "Missing tic3006blue/tic3006blue_multi.fits" >&2
  fi

  if [[ -f "tic3006red/tic3006red_multi.fits" ]]; then
    cp "tic3006red/tic3006red_multi.fits" "Final-Science-Products/tic3006red_multi_${SCIENCE_INDEX}.fits"
  else
    echo "Missing tic3006red/tic3006red_multi.fits" >&2
  fi

  echo "$(date '+%Y-%m-%d %H:%M:%S') Completed ${SCIENCE_INDEX}"
done
