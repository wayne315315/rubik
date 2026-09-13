# rubik

A small, NumPy-based Rubik's cube toolkit: a coordinate model of the cube, a
permutation-based rotation system, a brute-force solver, and a 3D animation
exporter built on matplotlib.

## Features

- **Coordinate model** – every one of the 54 stickers is a 4-vector
  `(x, y, z, n)` where `x, y, z ∈ {-1, 0, 1}` is the cubie position and
  `n ∈ {0, 1, 2}` is the axis of the sticker's normal vector.
- **18 rotations** – every slice (`+1`, `0`, `-1`) on every axis (`x`, `y`, `z`)
  in both directions, each precomputed as a permutation of sticker indices.
- **Composable sequences** – `RotationSequence` composes rotations, compares
  them by their resulting permutation, and supports `+` for concatenation.
- **Brute-force solver** – a multiprocess, breadth-first search over all
  move sequences of increasing length.
- **Video export** – render a scramble or a solution as a rotating 3D
  animation and save it as `.mp4`.

## Requirements

- Python 3.10+
- [ffmpeg](https://ffmpeg.org/) (only needed for video export)

```bash
sudo apt install ffmpeg          # Debian / Ubuntu
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Quick start

### Rotations

Rotations are named `<axis><level><orient>`:

| Part     | Values                                                  |
| -------- | ------------------------------------------------------- |
| axis     | `x`, `y`, `z`                                           |
| level    | `p` (slice at +1), `0` (middle slice), `n` (slice at -1) |
| orient   | `p` (positive, right-hand rule), `n` (negative)          |

So `xpp` turns the `x = +1` slice in the positive direction, and `y0n` turns
the middle `y` slice in the negative direction. All 18 are available as module
attributes in `r3.py`, together with the list `rs` that holds them all.

```python
import numpy as np
from r3 import index, coords, rs, xpp, y0p, zpn, RotationSequence

# A rotation is callable on either an index array or a coordinate array.
state = xpp(index)            # permuted sticker indices
state = y0p(state)

# Or compose a sequence.
scramble = RotationSequence([xpp, y0p, zpn])
state = scramble(index)

# Sequences are compared by the permutation they produce.
assert RotationSequence([xpp, xpp]) == RotationSequence([xpp] * 2)
assert RotationSequence([xpp]) + RotationSequence([y0p]) == RotationSequence([xpp, y0p])

# Coordinates work the same way, useful for drawing.
moved_coords = scramble(coords)
```

The solved state is `index`, the identity permutation `0..53`. A state is
solved when `np.all(state == index)`.

> **Note on middle slices.** A middle-slice turn such as `x0p` is implemented
> as turning both outer slices the opposite way. The resulting cube is the
> same as a true middle-slice turn up to a whole-cube reorientation, which
> keeps the face centers fixed in absolute coordinates.

### Solving

```python
from r3 import index, xpp, y0p, zpn, y0n
from solver import brute_force_multi

scramble = (xpp, y0p, zpn, y0n)
state = index
for r in scramble:
    state = r(state)

answers = brute_force_multi(state)      # list of solving sequences
best = min(answers, key=len)
print(" -> ".join(map(str, best)))
```

`brute_force_multi` enumerates every sequence of length 0, 1, 2, … and splits
the work round-robin across all CPU cores. It returns as soon as any worker
finds a solution, collecting every solution found by that moment. Because the
search space grows as 18ⁿ, this is practical for short scrambles (roughly six
moves or fewer) and is meant as a reference implementation rather than a
competitive solver.

A plain single-process depth-first search, `dfs`, is also provided.

Running the solver directly generates a random scramble and solves it:

```bash
python solver.py
```

### Video export

```bash
python visual.py
```

This scrambles the cube with a fixed four-move sequence, solves it, and writes
two files:

- `question.mp4` – the scramble applied to a solved cube.
- `answer.mp4` – the solution applied to the scrambled cube.

Each move is shown as one full camera orbit, with a red arc indicating the
slice and direction of the upcoming turn and the previous and next move
labelled in the corner.

To render your own sequence:

```python
from r3 import coords, xpp, ynn, z0p
from visual import export_video

export_video(coords, (xpp, ynn, z0p), "my_sequence.mp4")
```

Face colors follow the standard scheme: `x+` blue, `x-` green, `y+` red,
`y-` orange, `z+` yellow, `z-` white.

## Project layout

| File              | Purpose                                                        |
| ----------------- | -------------------------------------------------------------- |
| `r3.py`           | Coordinate system, index mappings, `Rotation`, `RotationSequence` |
| `solver.py`       | `dfs` and the multiprocess `brute_force_multi` solver           |
| `visual.py`       | matplotlib 3D drawing and `export_video`                        |
| `requirements.txt`| Pinned Python dependencies                                     |

## How rotations work

A 90° turn about an axis is computed with a cross product between the axis
unit vector and each sticker position. Stickers outside the rotating slice
are left untouched. Stickers inside the slice move to the cross-product
position, and their normal axis is remapped: a normal parallel to the rotation
axis stays put, while a normal perpendicular to it swaps to the other
perpendicular axis.

Each `Rotation` runs this once for a positive turn and three times for a
negative turn, then records the result as a permutation of the 54 sticker
indices. Applying a rotation afterwards is a single NumPy gather, which is
what makes the brute-force search feasible.
