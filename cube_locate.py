"""Locate-then-measure perception: the model only finds the stickers.

The vision-language model is asked for the bounding boxes of the 27 visible
stickers, grouped by face. Everything else is geometry and pixels:

  1. snap      each box centre is moved to the centroid of the sticker blob
               under it (stickers are bright patches separated by black gaps)
  2. lattice   the nine points of a face are sorted into a 3x3 grid from the
               middle point and its four nearest neighbours (two basis vectors)
  3. roles     which face is top / lower-left / lower-right, and which lattice
               direction runs along which shared edge, from positions alone
  4. colours   the median colour of a patch around every point, classified into
               the six sticker colours with the constraint that each colour
               appears exactly nine times across the two photos

The result is the same {top, lower_left, lower_right} grid dict that
cube_vision.state_from_views consumes, so validation, repair and solving are
unchanged. Position slips and misread centre stickers, the failure modes of
asking the model for colours, cannot happen here; the remaining risks are a
missed or duplicated sticker box (detected: every face needs nine points on a
clean lattice) and colour confusion between red and orange or white and
yellow (bounded by the nine-per-colour constraint).
"""
import base64
import io
import json
import math

import numpy as np
from PIL import Image, ImageOps

GRIDS = ("top", "lower_left", "lower_right")
COLORS = ("blue", "green", "orange", "red", "white", "yellow")

LOCATE_PROMPT = """This photo shows a 3x3x3 Rubik's cube with three visible faces, 27 coloured stickers in
total. Detect EVERY visible sticker and output its bounding box. Group the boxes by face:
one group per visible face, each with exactly 9 stickers.
Output ONLY JSON: {"faces": [{"stickers": [[x1,y1,x2,y2], ...]}, {"stickers": [...]},
{"stickers": [...]}]} with coordinates normalised to 0-1000 of the image width and height."""

LOCATE_SCHEMA = {"type": "object", "properties": {"faces": {"type": "array", "items": {
    "type": "object",
    "properties": {"stickers": {"type": "array", "items": {
        "type": "array", "minItems": 4, "maxItems": 4, "items": {"type": "number"}}}},
    "required": ["stickers"], "additionalProperties": False}}},
    "required": ["faces"], "additionalProperties": False}


class LocateError(Exception):
    """The boxes do not form three clean 3x3 faces."""


# ----------------------------------------------------------------------------
# image helpers
# ----------------------------------------------------------------------------
def load_image(path, max_side=1024):
    """EXIF-corrected RGB image downscaled to max_side, as a numpy uint8 array."""
    img = ImageOps.exif_transpose(Image.open(path)).convert("RGB")
    scale = max_side / max(img.size)
    if scale < 1:
        img = img.resize((round(img.width * scale), round(img.height * scale)), Image.LANCZOS)
    return np.asarray(img)


def jpeg_b64(arr):
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="JPEG", quality=92)
    return base64.b64encode(buf.getvalue()).decode()


def hsv(rgb):
    """HSV in [0,1] for an (..., 3) float array in [0,1]."""
    r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    mx, mn = rgb.max(-1), rgb.min(-1)
    d = mx - mn + 1e-9
    h = np.where(mx == r, ((g - b) / d) % 6, np.where(mx == g, (b - r) / d + 2, (r - g) / d + 4)) / 6.0
    s = d / (mx + 1e-9)
    return np.stack([h, s, mx], -1)


# ----------------------------------------------------------------------------
# 1. snap boxes to sticker blobs
# ----------------------------------------------------------------------------
def snap(arr, cx, cy, size):
    """Refine a sticker centre: among positions within +-25 % of the box size,
    pick the one whose small patch is most uniform in colour and not dark.
    Sticker interiors are flat; gaps, edges, logos and the table are not."""
    h, w = arr.shape[:2]
    r = max(2, int(size * 0.18))
    best, best_score = (cx, cy), None
    for dy in np.linspace(-0.25, 0.25, 7) * size:
        for dx in np.linspace(-0.25, 0.25, 7) * size:
            x, y = cx + dx, cy + dy
            x0, x1 = int(max(0, x - r)), int(min(w, x + r + 1))
            y0, y1 = int(max(0, y - r)), int(min(h, y + r + 1))
            if x1 - x0 < 3 or y1 - y0 < 3:
                continue
            patch = arr[y0:y1, x0:x1].reshape(-1, 3).astype(np.float32) / 255
            v = patch.max(1).mean()
            score = patch.std(0).sum() + (0.5 if v < 0.3 else 0.0) + 0.02 * math.hypot(dx, dy) / size
            if best_score is None or score < best_score:
                best, best_score = (x, y), score
    return best


# ----------------------------------------------------------------------------
# 2. lattice ordering of nine points
# ----------------------------------------------------------------------------
def lattice(points, max_resid=0.35, allow_missing=0):
    """Fit a 3x3 lattice to N points. Returns (pts9, used, leftover): pts9 is a (9, 2)
    array ordered by cell (row-major, cell k = (k // 3, k % 3)), used maps cell index
    -> input index (or None for a cell predicted from the affine fit), leftover lists
    the input indices not used.

    Every point is tried as the middle cell and every pair of the others as a
    basis; the candidate whose best-fitting points cover {-1,0,1}^2 (up to
    allow_missing cells absent) with the smallest least-squares affine residual
    wins. This copes with perspective shear (a corner cell can be nearer to the
    middle than an edge cell), with extra points from another face, and with a
    sticker the model missed."""
    pts = np.asarray(points, dtype=float)
    n = len(pts)
    if n < 9 - allow_missing:
        raise LocateError(f"a face needs at least {9 - allow_missing} stickers, got {n}")
    full = {(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)}
    best = None
    for ci in range(n):
        rel = pts - pts[ci]
        for a in range(n):
            for b in range(a + 1, n):
                if ci in (a, b):
                    continue
                basis = np.stack([rel[a], rel[b]], 1)
                if abs(np.linalg.det(basis)) < 1e-6:
                    continue
                ab = rel @ np.linalg.inv(basis).T
                rounded = np.round(ab)
                ok = (np.abs(ab - rounded).max(1) <= 0.45) & (np.abs(rounded).max(1) <= 1)
                chosen = {}
                for k in np.flatnonzero(ok):                       # one point per cell, best residual
                    cell = (int(rounded[k][0]), int(rounded[k][1]))
                    r = np.abs(ab[k] - rounded[k]).max()
                    if cell not in chosen or r < chosen[cell][1]:
                        chosen[cell] = (k, r)
                if len(full - set(chosen)) > allow_missing:
                    continue
                cells = sorted(chosen)
                sel = np.array([chosen[c][0] for c in cells])
                grid = np.array(cells, dtype=float)
                A = np.hstack([grid, np.ones((len(cells), 1))])
                coef, *_ = np.linalg.lstsq(A, pts[sel], rcond=None)
                spacing = min(np.linalg.norm(coef[0]), np.linalg.norm(coef[1]))
                resid = np.linalg.norm(A @ coef - pts[sel], axis=1).max() / spacing
                score = resid + 0.1 * (9 - len(cells))             # prefer complete faces
                if resid < max_resid and (best is None or score < best[0]):
                    best = (score, chosen, coef)
    if best is None:
        raise LocateError("sticker points do not sit on a 3x3 lattice")
    _, chosen, coef = best
    pts9, used = [], []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if (i, j) in chosen:
                k = chosen[(i, j)][0]
                pts9.append(pts[k])
                used.append(int(k))
            else:                                                  # predicted position
                pts9.append(coef[0] * i + coef[1] * j + coef[2])
                used.append(None)
    leftover = [k for k in range(n) if k not in used]
    return np.asarray(pts9), used, leftover


def dedupe(points, min_dist):
    """Merge points closer than min_dist (keeps the first of each cluster)."""
    kept = []
    for p in points:
        if all(np.linalg.norm(p - q) >= min_dist for q in kept):
            kept.append(np.asarray(p, dtype=float))
    return kept


def regroup(groups, strict=0.28):
    """Turn the model's rough face groups into three (9, 2) arrays ordered by cell.

    The model's grouping is only a hint (it is wrong on tilted photos), so the
    de-duplicated points are split geometrically: pull out the nine points that
    best form a lattice, then the best nine of the rest, then the last nine. A
    face may be completed from eight points when the model missed a sticker.
    If that fails, fall back to the model's groups with leftovers passed on."""
    all_pts = [np.asarray(p, dtype=float) for g in groups for p in g]
    if len(all_pts) < 24:
        raise LocateError(f"only {len(all_pts)} sticker boxes for three faces")
    d = np.sort(np.linalg.norm(np.array(all_pts)[:, None] - np.array(all_pts)[None], axis=2), 1)[:, 1]
    spacing = float(np.median(d))
    pool = dedupe(all_pts, 0.4 * spacing)
    remaining, faces = list(pool), []
    try:
        for step in range(3):
            faces_left = 3 - step
            allow = 1 if len(remaining) < 9 * faces_left else 0
            pts9, used, left = lattice(remaining, max_resid=strict, allow_missing=allow)
            faces.append(pts9)
            remaining = [remaining[i] for i in left]
        return faces
    except LocateError:
        pass
    clean, seen = [], []
    for g in groups:
        keep = []
        for p in g:
            p = np.asarray(p, dtype=float)
            if all(np.linalg.norm(p - q) >= 0.4 * spacing for q in seen):
                seen.append(p)
                keep.append(p)
        clean.append(keep)
    if len(clean) != 3:
        raise LocateError(f"need 3 faces, got {len(clean)}")
    faces, leftovers = [None] * 3, []
    for k in sorted(range(3), key=lambda k: -len(clean[k])):
        cand = clean[k] + leftovers
        pts9, used, left = lattice(cand, allow_missing=1)
        faces[k] = pts9
        leftovers = [cand[i] for i in left]
    return faces


# ----------------------------------------------------------------------------
# 3. face roles and grid orientation
# ----------------------------------------------------------------------------
def assign_faces(face_points):
    """face_points: list of three (9,2) arrays. Returns {grid_name: (points, cell_of_index)}
    where cell_of_index maps point index -> (i, j) in the image-relative grid convention
    of cube_vision (top[i][j], lower_left[i][j], lower_right[i][j])."""
    if len(face_points) != 3:
        raise LocateError(f"need 3 faces, got {len(face_points)}")
    # each face is a (9, 2) array ordered by cell, as returned by lattice()/regroup()
    lat = [{k: (k // 3, k % 3) for k in range(9)} for _ in face_points]
    cents = [np.asarray(p).mean(0) for p in face_points]
    C = np.mean(np.concatenate([np.asarray(p) for p in face_points]), 0)     # cube centre
    # clockwise angle of each face centroid around C, measured from straight up
    ang = [(math.degrees(math.atan2(c[0] - C[0], -(c[1] - C[1]))) % 360) for c in cents]
    top = int(np.argmin([min(a, 360 - a) for a in ang]))
    rest = sorted([k for k in range(3) if k != top], key=lambda k: (ang[k] - ang[top]) % 360)
    lr, ll = rest[0], rest[1]                          # clockwise from top: lower-right, then lower-left
    roles = {"top": top, "lower_right": lr, "lower_left": ll}
    others = {"top": ("lower_left", "lower_right"), "lower_left": ("top", "lower_right"),
              "lower_right": ("top", "lower_left")}
    out = {}
    for name, k in roles.items():
        pts = np.asarray(face_points[k], dtype=float)
        idx = lat[k]
        cell_pt = {ab: pts[i] for i, ab in idx.items()}
        near = min(cell_pt, key=lambda ab: np.linalg.norm(cell_pt[ab] - C))      # near-corner cell
        if near not in {(0, 0), (0, 2), (2, 0), (2, 2)}:
            raise LocateError(f"{name}: the cell nearest the cube centre is not a corner cell")
        da = (1 if near[0] == 0 else -1, 0)           # lattice step away from the corner along a
        db = (0, 1 if near[1] == 0 else -1)           # ... along b
        # which step runs along the edge shared with which neighbouring face:
        # the first cell along a step is closer to the centroid of that neighbour
        nA, nB = others[name]
        cA, cB = cents[roles[nA]], cents[roles[nB]]
        pa = cell_pt[(near[0] + da[0], near[1] + da[1])]
        pb = cell_pt[(near[0] + db[0], near[1] + db[1])]
        a_is_A = (np.linalg.norm(pa - cA) - np.linalg.norm(pa - cB)) < (np.linalg.norm(pb - cA) - np.linalg.norm(pb - cB))
        dA, dB = (da, db) if a_is_A else (db, da)
        cell_of = {}
        for i_pt, ab in idx.items():
            sa = (ab[0] - near[0]) * dA[0] + (ab[1] - near[1]) * dA[1]     # steps along the A edge
            sb = (ab[0] - near[0]) * dB[0] + (ab[1] - near[1]) * dB[1]     # steps along the B edge
            if name == "top":            # A = lower_left, B = lower_right
                ij = (2 - sb, 2 - sa)
            elif name == "lower_left":   # A = top, B = lower_right
                ij = (sb, 2 - sa)
            else:                        # A = top, B = lower_left
                ij = (sb, sa)
            cell_of[i_pt] = ij
        out[name] = (pts, cell_of)
    return out


# ----------------------------------------------------------------------------
# 4. colour measurement and balanced classification
# ----------------------------------------------------------------------------
def patch_color(arr, x, y, r):
    h, w = arr.shape[:2]
    x0, x1 = max(0, int(x - r)), min(w, int(x + r) + 1)
    y0, y1 = max(0, int(y - r)), min(h, int(y + r) + 1)
    return np.median(arr[y0:y1, x0:x1].reshape(-1, 3), 0) / 255.0


def features(rgb):
    """(s*cos h, s*sin h, v): hue direction scaled by saturation, plus brightness."""
    x = hsv(np.asarray(rgb, dtype=float))
    ang = x[..., 0] * 2 * np.pi
    return np.stack([x[..., 1] * np.cos(ang), x[..., 1] * np.sin(ang), x[..., 2]], -1)


REFERENCE = {                 # typical sticker colours under daylight, as (hue deg, sat, val)
    "blue": (212, 0.7, 0.9), "green": (105, 0.65, 0.8), "orange": (25, 0.75, 0.9),
    "red": (4, 0.72, 0.85), "white": (0, 0.1, 0.95), "yellow": (56, 0.65, 0.95),
}


def _ref_features():
    out = []
    for c in COLORS:
        hh, ss, vv = REFERENCE[c]
        a = math.radians(hh)
        out.append([ss * math.cos(a), ss * math.sin(a), vv])
    return np.array(out)


def _balanced(feat, centres, per_color):
    """Greedy balanced assignment of samples to colour centres (per_color each)."""
    d = np.linalg.norm(feat[:, None] - centres[None], axis=2)
    n = len(feat)
    order = np.dstack(np.unravel_index(np.argsort(d, axis=None), d.shape))[0]
    label, count = [None] * n, [0] * len(COLORS)
    for i, c in order:
        if label[i] is None and count[c] < per_color:
            label[i], count[c] = int(c), count[c] + 1
        if all(l is not None for l in label):
            break
    return np.array(label)


def classify(rgb, per_color=9, iterations=4):
    """Classify sticker colours with exactly per_color samples per colour.

    Starts from reference colours, then alternates between recomputing each
    colour's mean feature from its members and re-assigning, so the photo's
    white balance and exposure are absorbed."""
    feat = features(rgb)
    centres = _ref_features()
    label = _balanced(feat, centres, per_color)
    for _ in range(iterations):
        new = np.array([feat[label == c].mean(0) if (label == c).any() else centres[c] for c in range(len(COLORS))])
        nl = _balanced(feat, new, per_color)
        centres = new
        if (nl == label).all():
            break
        label = nl
    return [COLORS[c] for c in label]


# ----------------------------------------------------------------------------
# glue: locate with the model, then measure
# ----------------------------------------------------------------------------
def locate_points(reader, path, max_side=1024, seed=None):
    """Ask the model for sticker boxes; return (image array, list of point arrays per face)."""
    arr = load_image(path, max_side)
    h, w = arr.shape[:2]
    out = reader.chat_raw(jpeg_b64(arr), LOCATE_PROMPT, LOCATE_SCHEMA, think=False, seed=seed)
    faces = out.get("faces", [])
    groups = []
    for f in faces:
        pts = []
        for b in f.get("stickers", []):
            x1, y1, x2, y2 = [float(t) for t in b]
            cx, cy = (x1 + x2) / 2 / 1000 * w, (y1 + y2) / 2 / 1000 * h
            size = max(8.0, ((x2 - x1) / 1000 * w + (y2 - y1) / 1000 * h) / 2)
            pts.append(snap(arr, cx, cy, size))
        groups.append(np.asarray(pts, dtype=float))
    return arr, groups


def measure(arr, groups):
    """Geometry + colour samples for one photo: returns (cells, samples) where cells is a
    list of (grid, i, j) and samples the matching (n,3) colours."""
    roles = assign_faces(regroup(groups))
    cells, samples = [], []
    for name, (pts, cell_of) in roles.items():
        # patch radius: a quarter of the smallest lattice spacing
        d = np.sort(np.linalg.norm(pts[:, None] - pts[None], axis=2), 1)[:, 1].min()
        for k, ij in cell_of.items():
            cells.append((name, ij[0], ij[1]))
            samples.append(patch_color(arr, pts[k][0], pts[k][1], max(2, d * 0.22)))
    return cells, np.asarray(samples)


def read_photos(reader, paths, max_side=1024):
    """Two photos -> list of view dicts (one per photo) using locate-then-measure."""
    per_photo = []
    for p in paths:
        last = None
        for attempt in range(3):                                   # new seed each time
            try:
                arr, groups = locate_points(reader, p, max_side, seed=reader.seed + attempt)
                per_photo.append(measure(arr, groups))
                break
            except LocateError as e:
                last = e
        else:
            raise LocateError(f"{p}: {last}")
    all_samples = np.concatenate([s for _, s in per_photo])
    labels = classify(all_samples, per_color=9 if len(all_samples) == 54 else len(all_samples))
    views, k = [], 0
    for cells, samples in per_photo:
        view = {g: [[None] * 3 for _ in range(3)] for g in GRIDS}
        for (g, i, j) in cells:
            view[g][i][j] = labels[k]
            k += 1
        views.append(view)
    return views
