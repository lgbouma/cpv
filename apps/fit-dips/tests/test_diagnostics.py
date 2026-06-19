import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fit_dips import models, fitting, diagnostics  # noqa: E402


def test_write_diagnostics_creates_figures(tmp_path, monkeypatch):
    monkeypatch.setattr(diagnostics, "OUTBASE", str(tmp_path))

    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0, 0.25, 300))
    t_ref = float(np.median(t))
    flux = models.full_model("poly1", [1.0, 0.05], [0.03, 0.12, 0.01],
                             t, t_ref) + rng.normal(0, 0.001, 300)
    results, preferred = fitting.fit_all_models(
        t, flux, [(0.08, 0.16)], sigma0=0.001, t_ref=t_ref)

    outdir = diagnostics.write_diagnostics(
        {"id": "synthtest"}, t, flux, [(0.08, 0.16)], [], t, flux, t_ref,
        None, results, preferred)

    pngs = glob.glob(os.path.join(outdir, "*.png"))
    # one per model + overview + comparison + preferred residual hist
    assert len(pngs) == len(results) + 3
    assert any("__overview.png" in p for p in pngs)
    assert any("__model_comparison.png" in p for p in pngs)
    assert any(f"__model_{preferred}.png" in p for p in pngs)
