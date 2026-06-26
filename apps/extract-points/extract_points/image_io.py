"""Image loading and basic color/grayscale helpers."""
from __future__ import annotations

import os
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FIGURES = {
    "f12": os.path.join(HERE, "pasj_72_2_23_f12.jpeg"),
    "f13": os.path.join(HERE, "pasj_72_2_23_f13.jpeg"),
}


def load_rgb(figure: str) -> np.ndarray:
    """Return an (H, W, 3) uint8 RGB array for a named figure."""
    im = Image.open(FIGURES[figure]).convert("RGB")
    return np.asarray(im)


def gray(rgb: np.ndarray) -> np.ndarray:
    """Luminance, float32 in [0, 255]."""
    return rgb.astype(np.float32) @ np.array([0.299, 0.587, 0.114], np.float32)


def dark_mask(rgb: np.ndarray, thresh: float = 110.0) -> np.ndarray:
    """Anything dark (frame lines, markers, text): luminance below thresh."""
    return gray(rgb) < thresh
