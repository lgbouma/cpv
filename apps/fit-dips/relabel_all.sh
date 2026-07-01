#!/usr/bin/env bash
#
# Relabel + refit every dip from scratch, EXCEPT the Onitsuka 2016-02-09 group
# (gp/rp/zs), which is already done and is left untouched.
#
# Two subcommands:
#   clear    -- full reset (drop BOTH labels and fits + diagnostics) for every
#               dataset except the gp/rp/zs group. Run this ONCE. Destructive.
#   relabel  -- (default) launch the interactive labeler once per epoch group,
#               skipping any group already re-done. You draw the dip window
#               (shared across the group's bands), 'y' to fit (the cross-band
#               shared-duration joint fit runs automatically for >=2 bands),
#               'q' to close; the next group then opens. Re-runnable: it skips
#               groups already 'done', so you can relabel over several sittings.
#
# TIC 300651846 keeps only NB2.09 + TESS_S88 (per request): the raw TESS window
# and the TESS_pm2cycle fold are cleared and NOT relabeled.
#
# Usage:
#   cd apps/fit-dips
#   bash relabel_all.sh clear      # once, wipes labels+fits (except gp/rp/zs)
#   bash relabel_all.sh            # walk the interactive labeling, resumable
set -uo pipefail

PY="${CPV_PYTHON:-/Users/luke/local/miniconda3/envs/cpv/bin/python}"
cd "$(dirname "$0")"

# one id per epoch group -> launching it auto-groups the night's other bands.
GROUP_PRIMARIES=(
  LP12-502_20231208_g
  LP12-502_20231215_g
  PTFO8-8695_20140119_I+z
  PTFO8-8695_20140223_I_C
  PTFO8-8695_20141227_I_C
  PTFO8-8695_20150110_I_C
  PTFO8-8695_20150123_I_C
  PTFO8-8695_20150210_I_C
  PTFO8-8695_20150214_I_C
  PTFO8-8695_20150223_I_C
  PTFO8-8695_20161125_I_C
  PTFO8-8695_20161209_I_C
  TIC20178925_20191107_rprime
  TIC262400835_20201213_g
  TIC262400835_20201215_g
  TIC262400835_20201216_g
  TIC435899024_20210525_B
)

# never touched (already jointly refit)
KEEP=(
  PTFO8-8695_20160209_gp
  PTFO8-8695_20160209_rp
  PTFO8-8695_20160209_zs
)

is_keep() { local id="$1"; for k in "${KEEP[@]}"; do [[ "$k" == "$id" ]] && return 0; done; return 1; }

status_of() {  # echo a dataset's status, or "missing"
  "$PY" - "$1" <<'PY' 2>/dev/null
import sys
from fit_dips import registry
try:
    print(registry.load(sys.argv[1])["status"])
except Exception:
    print("missing")
PY
}

launch() {  # id [extra app flags...]; skip if already done
  local id="$1"; shift
  local st; st="$(status_of "$id")"
  if [[ "$st" == "done" ]]; then echo "skip (done): $id"; return; fi
  echo ">>> $id $*   --  label dip(s), 'y' to fit, 'q' to close"
  "$PY" -m fit_dips.app --dataset "$id" "$@"
}

cmd="${1:-relabel}"

if [[ "$cmd" == "clear" ]]; then
  echo "============ CLEAR: drop labels+fits for all except gp/rp/zs ==========="
  for f in datasets/*.json; do
    id="$(basename "$f" .json)"
    if is_keep "$id"; then echo "keep  $id"; continue; fi
    "$PY" -m fit_dips.app --reset "$id"      # no --keep-labels: labels dropped too
  done
  echo "done. Now run:  bash relabel_all.sh"
  exit 0
fi

if [[ "$cmd" == "relabel" ]]; then
  echo "============ RELABEL: interactive, one group at a time ================="
  echo "(close each window with 'q' after 'y' to move to the next)"
  for primary in "${GROUP_PRIMARIES[@]}"; do
    launch "$primary"
  done
  launch TIC300651846_20250210_NB2.09 --solo      # not tied to the raw TESS window
  launch TIC300651846_20250210_TESS_S88           # full-sector folded variant
  echo
  echo "============ remaining unfit ==========================================="
  "$PY" -m fit_dips.app --list | grep -E "not_started|started_failed" \
    && echo "(^ expected leftovers: TIC300651846 raw TESS window + pm2cycle)" \
    || echo "everything fit."
  # Optional next step:  $PY -m fit_dips.make_table
  exit 0
fi

echo "unknown subcommand: $cmd  (use 'clear' or 'relabel')"
exit 2
