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
def detect(arr, ckpt=DEFAULT_CKPT, threshold=0.3):
    """Sticker centres (n,2) in the coordinates of the given RGB array."""
    model = load_model(ckpt)
    x, (scale, dx, dy) = letterbox(Image.fromarray(arr))
    logits = model(x[None])
    peaks = extract_peaks(logits[0, 0], threshold=threshold)
    return np.array([[(px * STRIDE) / scale - dx, (py * STRIDE) / scale - dy] for px, py, _ in peaks], dtype=float)


def locate_points(path, max_side=1024, ckpt=DEFAULT_CKPT):
    """(image array, [points]) like cube_locate.locate_points, from the CNN."""
    arr = cube_locate.load_image(path, max_side)
    pts = detect(arr, ckpt)
    if len(pts) >= 2:
        d = np.sort(np.linalg.norm(pts[:, None] - pts[None], axis=2), 1)[:, 1]
        size = float(np.median(d)) * 0.9
        pts = np.array([cube_locate.snap(arr, x, y, size) for x, y in pts])
    return arr, [pts]


def read_photos(paths, max_side=1024, ckpt=DEFAULT_CKPT):
    """Two photos -> list of view dicts, CNN localisation + pixel colours."""
    per_photo = []
    for p in paths:
        arr, groups = locate_points(p, max_side, ckpt)
        per_photo.append(cube_locate.measure(arr, groups))
    all_samples = np.concatenate([s for _, s in per_photo])
    labels = cube_locate.classify(all_samples, per_color=9 if len(all_samples) == 54 else len(all_samples))
    views, k = [], 0
    for cells, samples in per_photo:
        view = {g: [[None] * 3 for _ in range(3)] for g in cube_locate.GRIDS}
        for (g, i, j) in cells:
            view[g][i][j] = labels[k]
            k += 1
        views.append(view)
    return views
