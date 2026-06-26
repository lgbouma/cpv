"""Separate the black data markers from the red model curve.

Markers are black (low R,G,B); the model is red (high R, low G,B). JPEG
compression leaves a reddish halo around the model, so the red mask is dilated
a little before being subtracted from the black mask.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def red_mask(rgb: np.ndarray, dilate: int = 2) -> np.ndarray:
    R = rgb[..., 0].astype(np.int16)
    G = rgb[..., 1].astype(np.int16)
    B = rgb[..., 2].astype(np.int16)
    m = (R > 120) & (R - G > 40) & (R - B > 40)
    if dilate > 0:
        m = ndi.binary_dilation(m, iterations=dilate)
    return m


def black_mask(rgb: np.ndarray, thresh: int = 100) -> np.ndarray:
    R = rgb[..., 0].astype(np.int16)
    G = rgb[..., 1].astype(np.int16)
    B = rgb[..., 2].astype(np.int16)
    return (R < thresh) & (G < thresh) & (B < thresh)


def marker_mask(rgb: np.ndarray) -> np.ndarray:
    """Black markers with the (dilated) red model removed."""
    return black_mask(rgb) & ~red_mask(rgb)
