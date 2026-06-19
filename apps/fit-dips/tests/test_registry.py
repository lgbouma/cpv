import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fit_dips import registry  # noqa: E402


def _make_done_record(monkeypatch, tmp_path):
    monkeypatch.setattr(registry, "DATASETS_DIR", str(tmp_path))
    rec = registry.new_dataset("X_1_g", "X", "2020-01-01", "g", "Inst", 1,
                               data={"loader": "muscat_csv", "path": "p"},
                               period_day=0.5)
    rec["labels"] = {"dips": [[1.0, 2.0]], "flares": [[3.0, 3.1]],
                     "confirmed": True}
    rec["fits"] = {"poly1": {"BIC": 1.0}}
    rec["preferred_model"] = "poly1"
    rec["sigma_eff"] = 0.01
    rec["error_rescale_f"] = 1.2
    registry.set_status(rec, registry.STATUS_DONE)
    registry.save(rec)
    return rec


def test_reset_clears_fit_and_labels(monkeypatch, tmp_path):
    _make_done_record(monkeypatch, tmp_path)
    registry.reset("X_1_g")
    r = registry.load("X_1_g")
    assert r["status"] == registry.STATUS_NOT_STARTED
    assert r["fits"] == {}
    assert r["preferred_model"] is None
    assert r["labels"]["dips"] == [] and r["labels"]["flares"] == []
    assert "sigma_eff" not in r and "error_rescale_f" not in r
    # metadata preserved
    assert r["period_day"] == 0.5 and r["data"]["loader"] == "muscat_csv"


def test_reset_keep_labels(monkeypatch, tmp_path):
    _make_done_record(monkeypatch, tmp_path)
    registry.reset("X_1_g", keep_labels=True)
    r = registry.load("X_1_g")
    assert r["status"] == registry.STATUS_NOT_STARTED
    assert r["fits"] == {}
    assert r["labels"]["dips"] == [[1.0, 2.0]]   # labels kept
