"""Fixed-camera 3D renderer: every move is animated as a real slice turn.

The camera never moves. Two views are shown side by side, the (+x, +y, +z)
corner on the left and the opposite (-x, -y, -z) corner on the right, so every
face is visible. Each move rotates the moving slice smoothly through 90 degrees
(right-hand rule about the positive axis, exactly as r3 defines it), followed
by a short hold.

Note on r3 middle-slice rotations (x0p, ...): r3 implements them as the two
outer slices turning the opposite way, so that is what is animated here.

Usage::

    from visual import export_video
    export_video(coords, seq, "out.mp4")

    python visual.py          # scramble, solve with thistlethwaite, export
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from r3 import coords, b2i

# prerequisites : sudo apt install ffmpeg

# face colors
# (x+, x-, y+, y-, z+, z-) : (blue, green, red, orange, yellow, white)
hex = {
    (0, 1): "#3F51B5",
    (0, -1): "#4CAF50",
    (1, 1): "#B51F1F",
    (1, -1): "#FF6F00",
    (2, 1): "#FFEB3B",
    (2, -1): "#FFFFFF"
}

# colour of every sticker, by its home position (never changes when the cube turns)
colors = np.array([hex[(coords[i][-1], coords[i][coords[i][-1]])] for i in range(len(coords))])

CELLS = sorted(b2i)                                   # the 26 cubie positions
FACES = [(axis, sign) for axis in range(3) for sign in (1, -1)]
BODY_COLOR = "#111111"
STICKER_LIFT = 0.012                                  # keep stickers above the body
STICKER_INSET = 0.08                                  # black border around stickers


def draw_axes(ax, arrow_offset=3, text_offset=3.5):
    """Draw the x, y, z arrows with labels."""
    x0, y0, z0 = np.zeros((3, 3))
    x1, y1, z1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) * arrow_offset
    ax.quiver(x0, y0, z0, x1, y1, z1, arrow_length_ratio=0.1, color="black")
    for i, a in enumerate(["x", "y", "z"]):
        args = [0, 0, 0, 0]
        args[i] = text_offset
        args[-1] = a
        ax.text(*args, color="black")


def _quad(cell, axis, sign, half=0.5, lift=0.0):
    """Four corners of the (axis, sign) face of the unit cube centred on cell."""
    c = np.asarray(cell, dtype=float)
    i, j = [k for k in range(3) if k != axis]
    pts = []
    for di, dj in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        p = c.copy()
        p[axis] += sign * (0.5 + lift)
        p[i] += di * half
        p[j] += dj * half
        pts.append(p)
    return np.array(pts)


# Fixed geometry: 26 x 6 black body faces, then 54 sticker slots on the outside.
BODY_VERTS = np.array([_quad(cell, axis, sign) for cell in CELLS for axis, sign in FACES])
BODY_CELLS = [cell for cell in CELLS for _ in FACES]
SLOTS = [(cell, axis, sign) for cell in CELLS for axis, sign in FACES if cell[axis] == sign]
SLOT_VERTS = np.array([_quad(cell, axis, sign, half=0.5 - STICKER_INSET, lift=STICKER_LIFT)
                       for cell, axis, sign in SLOTS])
SLOT_INDEX = {slot: k for k, slot in enumerate(SLOTS)}
ALL_VERTS = np.concatenate([BODY_VERTS, SLOT_VERTS])
ALL_CELLS = np.array(BODY_CELLS + [cell for cell, _, _ in SLOTS])   # (210, 3)
N_BODY = len(BODY_VERTS)


def slot_colors(state):
    """Sticker colour per slot for an r3 coordinate array (54, 4)."""
    fc = [BODY_COLOR] * len(SLOTS)
    for i, (x, y, z, n) in enumerate(state):
        fc[SLOT_INDEX[((x, y, z), n, state[i][n])]] = colors[i]
    return fc


def rotation_matrix(axis, theta):
    """Right-hand rotation by theta about +axis (cyclic order x->y->z->x)."""
    c, s = np.cos(theta), np.sin(theta)
    i, j = (axis + 1) % 3, (axis + 2) % 3
    m = np.eye(3)
    m[i, i], m[i, j], m[j, i], m[j, j] = c, -s, s, c
    return m


def turning(r):
    """(moving-cell mask, signed angle) for an r3 Rotation."""
    if r.level == 0:                 # r3: both outer slices, opposite direction
        mask = ALL_CELLS[:, r.axis] != 0
        sign = -r.orient
    else:
        mask = ALL_CELLS[:, r.axis] == r.level
        sign = r.orient
    return mask, sign * np.pi / 2


def frame_verts(r, t):
    """Vertices with the slice of rotation r turned by fraction t in [0, 1]."""
    verts = ALL_VERTS.copy()
    if r is None or t == 0:
        return verts
    mask, angle = turning(r)
    m = rotation_matrix(r.axis, angle * t)
    verts[mask] = verts[mask] @ m.T
    return verts


def smoothstep(t):
    return t * t * (3 - 2 * t)


def _setup_axes(ax, elev, azim, limit=2.2):
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    draw_axes(ax)
    coll = Poly3DCollection(ALL_VERTS, edgecolor="black", linewidth=0.4, shade=False)
    ax.add_collection3d(coll)
    return coll


def export_video(coords_init, seq, filename, frames_per_turn=18, hold=6, fps=24,
                 views=((30, 45), (-30, 225)), figsize=(10, 5)):
    """Render seq applied to coords_init as a fixed-camera slice-turn video."""
    seq = list(seq)
    states = [np.asarray(coords_init)]
    for r in seq:
        states.append(r(states[-1]))

    # timeline: hold, then for each move (turn, hold)
    timeline = [(None, 0, 0.0)] * hold
    for k, r in enumerate(seq):
        timeline += [(r, k, smoothstep((t + 1) / frames_per_turn)) for t in range(frames_per_turn)]
        timeline += [(None, k + 1, 0.0)] * hold

    fig = plt.figure(figsize=figsize)
    fig.patch.set_facecolor("white")
    colls = [_setup_axes(fig.add_subplot(1, len(views), i + 1, projection="3d"), *v)
             for i, v in enumerate(views)]
    label = fig.text(0.02, 0.95, "", fontsize="large", family="monospace", va="top")

    def animate(f):
        r, k, t = timeline[f]
        verts = frame_verts(r, t)
        state = states[k]
        fc = [BODY_COLOR] * N_BODY + slot_colors(state)
        for coll in colls:
            coll.set_verts(verts)
            coll.set_facecolor(fc)
        prior = seq[k - 1] if k > 0 else None
        nxt = r if r is not None else (seq[k] if k < len(seq) else None)
        label.set_text(f"move {min(k + (r is not None), len(seq))}/{len(seq)}   "
                       f"prior: {prior}   next: {nxt}")
        return colls

    anim = animation.FuncAnimation(fig, animate, frames=len(timeline), interval=1000 / fps, blit=False)
    anim.save(filename, writer=animation.FFMpegWriter(fps=fps))
    plt.close(fig)


if __name__ == "__main__":
    import random
    from r3 import index, rs, RotationSequence
    from thistlethwaite import solve

    question = RotationSequence(random.choices(rs, k=10))
    coords_q = question(coords)
    answer = RotationSequence(solve(question(index)))
    print("Question: %s" % question)
    print("Answer: %s" % answer)
    export_video(coords, question.seq, "question.mp4")
    export_video(coords_q, answer.seq, "answer.mp4")
    print("Exported question.mp4 and answer.mp4")
