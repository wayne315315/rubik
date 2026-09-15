"""Synthetic corner-view photos of a 3x3x3 cube with free sticker-centre labels.

A small numpy/PIL rasteriser: 26 cubies drawn as black blocks with inset
coloured stickers, pinhole perspective camera placed near one of the two body
diagonals, Lambert shading with a coloured light, random background, then
photo-like degradation (blur, noise, JPEG, exposure, white balance). Every
sample returns the image and the projected centres of the 27 visible stickers
with their face id, which is all the keypoint detector needs.

    python -m synth.render            # writes synth/preview.jpg (4x4 samples)
"""
import io
import math
import os
import random

import numpy as np
from PIL import Image, ImageDraw, ImageFilter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# sticker colours (RGB 0-1), same scheme as visual.py, jittered per sample
BASE_COLORS = {
    (0, 1): (0.25, 0.32, 0.71), (0, -1): (0.30, 0.69, 0.31),
    (1, 1): (0.71, 0.12, 0.12), (1, -1): (1.00, 0.44, 0.00),
    (2, 1): (1.00, 0.92, 0.23), (2, -1): (0.96, 0.96, 0.96),
}
BODY = (0.06, 0.06, 0.07)


# ----------------------------------------------------------------------------
# cube geometry
# ----------------------------------------------------------------------------
def cube_quads(state_colors, gap=0.06, inset=0.10):
    """All quads of the cube: (points (4,3), colour, normal, sticker_id or None).

    state_colors: {(cell(x,y,z), axis, sign): colour} for the 54 stickers."""
    quads = []
    half = 0.5 - gap / 2
    sid = 0
    for x in (-1, 0, 1):
        for y in (-1, 0, 1):
            for z in (-1, 0, 1):
                if x == y == z == 0:
                    continue
                c = np.array([x, y, z], dtype=float)
                for axis in range(3):
                    for sign in (1, -1):
                        n = np.zeros(3)
                        n[axis] = sign
                        i, j = [k for k in range(3) if k != axis]
                        def corners(h, lift):
                            pts = []
                            for di, dj in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                                p = c.copy()
                                p[axis] += sign * (half + lift)
                                p[i] += di * h
                                p[j] += dj * h
                                pts.append(p)
                            return np.array(pts)
                        quads.append((corners(half, 0.0), BODY, n, None))
                        key = ((x, y, z), axis, sign)
                        if key in state_colors:                    # outer face: sticker
                            quads.append((corners(half - inset, 0.012), state_colors[key], n, sid))
                            sid += 1
    return quads


def random_state_colors(rng):
    """A random assignment of face colours to the 54 stickers, 9 per colour (a
    permutation; whether it is a reachable cube state is irrelevant for the
    detector, which only needs to find sticker centres). Returns
    {sticker key: face key (axis, sign)} so colours can be jittered per face."""
    keys = [((x, y, z), a, s) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)
            if not (x == y == z == 0) for a in range(3) for s in (1, -1) if (x, y, z)[a] == s]
    faces = [f for f in BASE_COLORS for _ in range(9)]
    rng.shuffle(faces)
    return dict(zip(keys, faces))


# ----------------------------------------------------------------------------
# camera and shading
# ----------------------------------------------------------------------------
def look_at(eye, target, roll):
    f = target - eye
    f /= np.linalg.norm(f)
    up = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(f, up)) > 0.99:
        up = np.array([0.0, 1.0, 0.0])
    r = np.cross(f, up)
    r /= np.linalg.norm(r)
    u = np.cross(r, f)
    cr, sr = math.cos(roll), math.sin(roll)
    r, u = cr * r + sr * u, -sr * r + cr * u
    return np.stack([r, u, f])                       # rows: right, up, forward


def project(points, eye, R, focal, size):
    """Pinhole projection of (n,3) world points to pixel coords (n,2) and depths."""
    rel = (points - eye) @ R.T
    depth = rel[:, 2]
    x = focal * rel[:, 0] / depth + size / 2
    y = -focal * rel[:, 1] / depth + size / 2
    return np.stack([x, y], 1), depth


def sample_camera(rng, size):
    """Camera near a body diagonal (either sign), random distance, fov and roll."""
    corner = np.array([rng.choice([-1, 1]) for _ in range(3)], dtype=float)
    d = corner / np.linalg.norm(corner)
    # random tilt away from the diagonal (up to ~35 degrees)
    perturb = rng.normal(size=3)
    perturb -= np.dot(perturb, d) * d
    perturb /= np.linalg.norm(perturb) + 1e-9
    ang = math.radians(rng.uniform(0, 35))
    d = math.cos(ang) * d + math.sin(ang) * perturb
    dist = rng.uniform(4.5, 8.0)
    eye = d * dist
    target = rng.normal(scale=0.25, size=3)
    R = look_at(eye, target, rng.uniform(0, 2 * math.pi))
    fov = math.radians(rng.uniform(35, 60))
    focal = (size / 2) / math.tan(fov / 2)
    return eye, R, focal


def shade(color, normal, light_dir, light_color, ambient, view_dir):
    lam = max(0.0, float(np.dot(normal, light_dir)))
    h = light_dir + view_dir
    h /= np.linalg.norm(h) + 1e-9
    spec = max(0.0, float(np.dot(normal, h))) ** 40 * 0.35
    c = np.array(color) * (ambient + (1 - ambient) * lam * np.array(light_color)) + spec
    return tuple(int(255 * min(1.0, max(0.0, v))) for v in c)


# ----------------------------------------------------------------------------
# backgrounds and degradation
# ----------------------------------------------------------------------------
_BG_CACHE = None


def background_images():
    """Photos to crop backgrounds from: everything under examples/ and uploads/."""
    global _BG_CACHE
    if _BG_CACHE is None:
        paths = []
        for d in ("examples", "uploads"):
            base = os.path.join(ROOT, d)
            for dirpath, _, files in os.walk(base):
                paths += [os.path.join(dirpath, f) for f in files if f.lower().endswith((".jpg", ".jpeg", ".png"))]
        _BG_CACHE = []
        for p in paths[:40]:
            try:
                im = Image.open(p).convert("RGB")
                im.thumbnail((1200, 1200))
                _BG_CACHE.append(im)
            except OSError:
                pass
    return _BG_CACHE


def random_background(rng, size):
    kind = rng.random()
    bgs = background_images()
    if kind < 0.6 and bgs:
        im = bgs[rng.integers(len(bgs))]
        w, h = im.size
        s = int(rng.uniform(0.25, 0.7) * min(w, h))
        # corners of photos are more likely to be cube-free; pick an edge strip
        x0 = int(rng.uniform(0, w - s))
        y0 = int(rng.choice([0, h - s])) if rng.random() < 0.7 else int(rng.uniform(0, h - s))
        bg = im.crop((x0, y0, x0 + s, y0 + s)).resize((size, size), Image.BILINEAR)
        bg = bg.rotate(rng.uniform(0, 360), resample=Image.BILINEAR, expand=False)
        return bg
    if kind < 0.8:                                   # procedural wood-like stripes
        y, x = np.mgrid[0:size, 0:size].astype(np.float32)
        a = rng.uniform(0, math.pi)
        t = x * math.cos(a) + y * math.sin(a)
        stripes = 0.5 + 0.5 * np.sin(t / rng.uniform(4, 20) + 2 * np.sin(t / rng.uniform(40, 120)))
        base = np.array([rng.uniform(0.35, 0.75), rng.uniform(0.22, 0.5), rng.uniform(0.1, 0.3)])
        img = (base[None, None] * (0.7 + 0.3 * stripes[..., None]) * 255).astype(np.uint8)
        return Image.fromarray(img)
    if kind < 0.9:                                   # gradient
        y, x = np.mgrid[0:size, 0:size].astype(np.float32) / size
        c1, c2 = rng.uniform(0.1, 0.9, 3), rng.uniform(0.1, 0.9, 3)
        t = (x * rng.uniform(-1, 1) + y * rng.uniform(-1, 1) + 1) / 2
        img = ((c1[None, None] * (1 - t[..., None]) + c2[None, None] * t[..., None]) * 255).astype(np.uint8)
        return Image.fromarray(img)
    return Image.new("RGB", (size, size), tuple(int(v) for v in rng.uniform(20, 235, 3)))


def degrade(img, rng):
    if rng.random() < 0.5:
        img = img.filter(ImageFilter.GaussianBlur(rng.uniform(0.3, 1.6)))
    arr = np.asarray(img).astype(np.float32)
    # white balance / exposure
    arr *= np.array([rng.uniform(0.85, 1.15), rng.uniform(0.9, 1.1), rng.uniform(0.85, 1.15)])[None, None]
    arr = arr * rng.uniform(0.7, 1.3) + rng.uniform(-20, 20)
    arr += rng.normal(scale=rng.uniform(0, 8), size=arr.shape)
    img = Image.fromarray(np.clip(arr, 0, 255).astype(np.uint8))
    if rng.random() < 0.7:                           # JPEG round trip
        buf = io.BytesIO()
        img.save(buf, format="JPEG", quality=int(rng.uniform(40, 95)))
        img = Image.open(io.BytesIO(buf.getvalue())).convert("RGB")
    return img


# ----------------------------------------------------------------------------
# sample
# ----------------------------------------------------------------------------
def render(rng, size=256, supersample=2):
    """One synthetic sample: (PIL image size x size, keypoints (27,2), face ids (27,))."""
    S = size * supersample
    assignment = random_state_colors(rng)
    jitter = {k: tuple(min(1, max(0, c + rng.normal(scale=0.04))) for c in v) for k, v in BASE_COLORS.items()}
    colors = {k: jitter[f] for k, f in assignment.items()}
    quads = cube_quads(colors, gap=rng.uniform(0.03, 0.09), inset=rng.uniform(0.06, 0.14))
    for _ in range(20):                                # keep the whole cube inside the frame
        eye, R, focal = sample_camera(rng, S)
        corners = np.array([[x, y, z] for x in (-1.5, 1.5) for y in (-1.5, 1.5) for z in (-1.5, 1.5)])
        pix, _ = project(corners, eye, R, focal, S)
        if pix.min() > 0.04 * S and pix.max() < 0.96 * S:
            break
    view_dir = -eye / np.linalg.norm(eye)
    light_dir = rng.normal(size=3)
    light_dir[2] = abs(light_dir[2]) + 0.5
    light_dir /= np.linalg.norm(light_dir)
    if np.dot(light_dir, -view_dir) < 0.2:           # keep the light roughly on the camera side
        light_dir = light_dir * 0.5 + (-view_dir) * 0.8
        light_dir /= np.linalg.norm(light_dir)
    warm = rng.uniform(-0.15, 0.15)
    light_color = (1 + warm, 1.0, 1 - warm)
    ambient = rng.uniform(0.25, 0.55)

    bg = random_background(rng, S)
    draw = ImageDraw.Draw(bg)
    drawn = []
    for pts, col, n, sid in quads:
        if np.dot(n, eye - pts.mean(0)) <= 0:          # back-facing
            continue
        pix, depth = project(pts, eye, R, focal, S)
        drawn.append((depth.mean(), pix, col, n, sid))
    drawn.sort(key=lambda t: -t[0])                    # far to near
    keypoints, faces = [], []
    for _, pix, col, n, sid in drawn:
        draw.polygon([tuple(p) for p in pix], fill=shade(col, n, light_dir, light_color, ambient, -view_dir))
        if sid is not None:
            keypoints.append(pix.mean(0))
            faces.append(int(np.argmax(np.abs(n))) * 2 + (0 if n[np.argmax(np.abs(n))] > 0 else 1))
    img = bg.resize((size, size), Image.LANCZOS) if supersample > 1 else bg
    img = degrade(img, rng)
    kp = np.asarray(keypoints) / supersample
    return img, kp, np.asarray(faces)


def preview(path=os.path.join(ROOT, "synth", "preview.jpg"), n=4, size=256, seed=0):
    rng = np.random.default_rng(seed)
    sheet = Image.new("RGB", (n * size, n * size))
    for k in range(n * n):
        img, kp, _ = render(rng, size)
        d = ImageDraw.Draw(img)
        for x, y in kp:
            d.ellipse((x - 2, y - 2, x + 2, y + 2), fill=(255, 0, 255))
        sheet.paste(img, ((k % n) * size, (k // n) * size))
    sheet.save(path, quality=90)
    return path


if __name__ == "__main__":
    print(preview())
