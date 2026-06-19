"""
Per-dataset JSON registry and processing-status tracking.

One JSON file per dataset lives in apps/fit-dips/datasets/<id>.json and holds
metadata, the data pointer, the user's labels, and the results of *all* fitted
models plus the preferred model.
"""
import json
import os
from datetime import datetime, timezone

DATASETS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "datasets")

STATUS_NOT_STARTED = "not_started"
STATUS_FAILED = "started_failed"
STATUS_DONE = "done"
STATUSES = (STATUS_NOT_STARTED, STATUS_FAILED, STATUS_DONE)


def _path(dataset_id):
    return os.path.join(DATASETS_DIR, f"{dataset_id}.json")


def new_dataset(dataset_id, star, date_ut, band, instrument, ref,
                data=None, period_day=None, note=""):
    """Construct a fresh dataset record dict (not yet saved)."""
    return {
        "id": dataset_id,
        "star": star,
        "date_ut": date_ut,
        "band": band,
        "instrument": instrument,
        "ref": ref,
        "period_day": period_day,
        "note": note,
        "data": data,  # None until a data file is resolved
        "status": STATUS_NOT_STARTED,
        "sigma0_p2p_rms": None,
        "labels": {"dips": [], "flares": [], "confirmed": False},
        "fits": {},
        "preferred_model": None,
        "updated": None,
    }


def save(record):
    os.makedirs(DATASETS_DIR, exist_ok=True)
    record["updated"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with open(_path(record["id"]), "w") as f:
        json.dump(record, f, indent=2)
    return _path(record["id"])


def load(dataset_id):
    with open(_path(dataset_id)) as f:
        return json.load(f)


def exists(dataset_id):
    return os.path.exists(_path(dataset_id))


def list_datasets():
    """Return all dataset records sorted by id."""
    if not os.path.isdir(DATASETS_DIR):
        return []
    out = []
    for fn in sorted(os.listdir(DATASETS_DIR)):
        if fn.endswith(".json"):
            with open(os.path.join(DATASETS_DIR, fn)) as f:
                out.append(json.load(f))
    return out


def set_status(record, status):
    assert status in STATUSES, status
    record["status"] = status
    return record


def reset(dataset_id, keep_labels=False):
    """Clear a dataset's fit results and mark it not_started.

    Drops fits, preferred model, and the derived fit metadata; by default also
    clears the dip/flare labels. The data pointer, period, and other dataset
    metadata are preserved. Returns the saved record.
    """
    rec = load(dataset_id)
    rec["fits"] = {}
    rec["preferred_model"] = None
    rec["sigma0_p2p_rms"] = None
    for key in ("t_ref", "error_rescale_f", "sigma_eff", "diagnostics_dir"):
        rec.pop(key, None)
    if not keep_labels:
        rec["labels"] = {"dips": [], "flares": [], "confirmed": False}
    set_status(rec, STATUS_NOT_STARTED)
    save(rec)
    return rec
