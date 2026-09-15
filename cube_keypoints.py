"""Sticker localisation with the trained CNN (no vision-language model needed).

Drop-in replacement for the model call in cube_locate: returns the same
(image array, [points]) so cube_locate.measure() does the face split, grid
orientation and pixel colours.

    from cube_keypoints import locate_points
    arr, groups = locate_points("photo.jpg")
"""
import os

import numpy as np
import torch
from PIL import Image

from synth.model import KeypointNet, extract_peaks
from synth.train import STRIDE, letterbox
import cube_locate

ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CKPT = os.path.join(ROOT, "synth", "ckpt", "best.pt")
ALT_CKPT = os.path.join(ROOT, "synth", "ckpt", "alt.pt")       # second opinion when the split fails
_MODEL = {}


def load_model(ckpt=DEFAULT_CKPT):
    if ckpt not in _MODEL:
        state = torch.load(ckpt, map_location="cpu")
        model = KeypointNet(state.get("width", 32))
        model.load_state_dict(state["model"])
        model.eval()
        _MODEL[ckpt] = model
    return _MODEL[ckpt]


@torch.no_grad()
def detect(arr, ckpt=DEFAULT_CKPT, threshold=0.25):
    """Sticker centres (n,3): x, y in the coordinates of the given RGB array and
    the peak score, sorted by score (highest first)."""
    model = load_model(ckpt)
    x, (scale, dx, dy) = letterbox(Image.fromarray(arr))
    logits = model(x[None])
    peaks = extract_peaks(logits[0, 0], threshold=threshold)
    out = np.array([[(px * STRIDE) / scale - dx, (py * STRIDE) / scale - dy, sc] for px, py, sc in peaks], dtype=float)
    return out.reshape(-1, 3)


def locate_points(path, max_side=1024, ckpt=DEFAULT_CKPT, keep=27):
    """(image array, [points]) like cube_locate.locate_points, from the CNN.
    Only the `keep` highest-scoring peaks are used (extras are usually spurious)."""
    arr = cube_locate.load_image(path, max_side)
    det = detect(arr, ckpt)[:keep]
    pts = det[:, :2]
    if len(pts) >= 2:
        d = np.sort(np.linalg.norm(pts[:, None] - pts[None], axis=2), 1)[:, 1]
        size = float(np.median(d)) * 0.9
        pts = np.array([cube_locate.snap(arr, x, y, size) for x, y in pts])
    return arr, [pts]


def measure_photo(path, max_side=1024, ckpt=DEFAULT_CKPT):
    """cube_locate.measure() on the CNN points; retries with a few extra peaks
    when the top-27 do not split into three faces."""
    last = None
    ckpts = [ckpt] + ([ALT_CKPT] if ckpt == DEFAULT_CKPT and os.path.exists(ALT_CKPT) else [])
    for ck in ckpts:
        for keep in (27, 28, 30, 33):
            try:
                arr, groups = locate_points(path, max_side, ck, keep)
                return cube_locate.measure(arr, groups)
            except cube_locate.LocateError as e:
                last = e
    raise cube_locate.LocateError(f"{path}: {last}")


def read_photos(paths, max_side=1024, ckpt=DEFAULT_CKPT):
    """Two photos -> list of view dicts, CNN localisation + pixel colours."""
    return cube_locate.label_views([measure_photo(p, max_side, ckpt) for p in paths])
