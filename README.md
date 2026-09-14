# rubik

A small, NumPy-based Rubik's cube toolkit: a coordinate model of the cube, a
permutation-based rotation system, a fast Thistlethwaite solver, an optimal
IDA* solver, and 3D animation exporters built on matplotlib.

## Features

- **Coordinate model** – every one of the 54 stickers is a 4-vector
  `(x, y, z, n)` where `x, y, z ∈ {-1, 0, 1}` is the cubie position and
  `n ∈ {0, 1, 2}` is the axis of the sticker's normal vector.
- **18 rotations** – every slice (`+1`, `0`, `-1`) on every axis (`x`, `y`, `z`)
  in both directions, each precomputed as a permutation of sticker indices.
- **Composable sequences** – `RotationSequence` composes rotations, compares
  them by their resulting permutation, and supports `+` for concatenation.
- **Fast solver** – a Thistlethwaite four-phase solver with complete
  distance tables per phase. Solves any scramble in about 2 ms with roughly
  50 rotations, after a one-off 2 s table build that is cached to disk.
- **Optimal solver** – IDA* with four pattern databases, vectorised with
  numpy. Guarantees the shortest solution in the r3, quarter-turn, or
  half-turn metric.
- **God's number bounds** – a script that bounds God's number for the r3
  move set by counting canonical sequences and by solving hard positions.
- **Photos to moves** – `cube_vision.py` reads two corner-view photos, taken at
  any angle, with an open-source vision-language model (Qwen on a local Ollama
  server), validates and repairs the reading, and returns the solving and
  scrambling `RotationSequence`s. `upload_server.py` does the same from a phone
  browser. See [COOKBOOK.md](COOKBOOK.md).
- **Brute-force solver** – a multiprocess, breadth-first search over all
  move sequences of increasing length, kept as a reference implementation.
- **Video export** – render a scramble or a solution as an `.mp4` with a fixed
  camera and smoothly turning slices (`visual.py`).

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
from r3 import index, rs, RotationSequence
from thistlethwaite import solve
import random

question = RotationSequence(random.choices(rs, k=50))
state = question(index)

answer = RotationSequence(solve(state))       # list of r3 Rotation objects
assert answer(state).tolist() == index.tolist()
print(answer)
```

`solve` accepts either an index vector or a coordinate array and returns a
list of `Rotation` objects, so the result plugs straight into
`RotationSequence` and `visual.export_video`. It raises `ValueError` for a
state that no sequence of rotations can produce.

The first call builds four breadth-first tables (about 2 s) and caches them
in `thistlethwaite_tables.pickle`. Later runs load the cache. Run the module
directly to solve and verify a batch of random scrambles:

```bash
python thistlethwaite.py 500
```

#### How the fast solver works

The solver walks the cube group down a chain of nested subgroups, using the
`z` axis as up/down, `x` as front/back, and `y` as right/left:

| Phase | Goal                                                 | Moves allowed                   | States    | Max depth |
| ----- | ---------------------------------------------------- | ------------------------------- | --------- | --------- |
| 1     | Orient all 12 edges                                  | any                             | 2,048     | 6         |
| 2     | Orient all 8 corners, E-slice edges into the E slice | `z`, `y` quarter turns, half turns | 1,082,565 | 10        |
| 3     | M-slice edges into the M slice, corners into tetrads | `z` quarter turns, half turns   | 2,822,400 | 13        |
| 4     | Finish                                               | half turns only                 | 663,552   | 15        |

Each phase tracks a small coordinate (for example the positions of the
twelve edge reference stickers) that fully captures its sub-problem. A
complete distance table is built over that coordinate space, so each phase is
solved by greedily stepping to any neighbour one move closer to the goal.
Every move's effect on a coordinate is derived by applying the r3 `Rotation`
objects themselves, so no cube mechanics are hard-coded in the solver.

#### Optimal solver

```python
from optimal import solve

answer = solve(state)                 # fewest r3 rotations, guaranteed
answer = solve(state, metric="htm")   # fewest half-turn-metric moves
```

`optimal.solve` runs IDA* (iterative deepening A*) with an admissible
heuristic taken from four pattern databases, so the first solution it finds
is provably the shortest for the chosen metric:

| Metric | Moves, each costing 1                             |
| ------ | ------------------------------------------------- |
| `r3`   | the 18 r3 rotations (12 quarter turns + 6 middle-slice rotations) |
| `qtm`  | 12 quarter turns                                  |
| `htm`  | 12 quarter turns + 6 half turns                   |

The state is a 6-tuple of small coordinates (corner orientation and
permutation, edge orientation, and the ordered positions of the E, M, and S
slice edges), each with a tiny move table. The pattern databases are exact
distances over pairs of those coordinates (88 M + 3 × 24 M entries, about
160 MB), built once in about 35 s and cached in `optimal_tables_<metric>.pickle`.
The search itself is depth-first over numpy arrays of nodes, examining tens of
millions of nodes per second per core, and prunes redundant same-axis move
runs with a small automaton. Once an iteration grows past a million nodes,
the tree is split a few levels below the root into thousands of independent
subtrees that are searched on every CPU core in parallel (`workers` argument,
default all cores). Any solution found in an iteration is optimal because all
smaller bounds were exhausted first.

Running time still grows exponentially with the optimal length, roughly a
factor of 10 per extra move. On a 192-core machine a fully random state
(optimal 18 r3 rotations, about 5×10¹⁰ nodes) takes about 35 s. Run the module
directly to compare against the Thistlethwaite answer on a random scramble:

```bash
python optimal.py r3 12          # metric, scramble length, optional seed
```

#### God's number for the r3 move set

```bash
python gods_number.py              # rigorous counting lower bound, instant
python gods_number.py --random 5   # optimal lengths of 5 random states
python gods_number.py --hard       # four-spot, superflip, superflip + four-spot
```

The exact value is out of reach on a single machine (the half-turn-metric
result of 20 took about 35 CPU-years). The script prints what can be
established: a lower bound from counting canonical sequences (18 for the r3
move set), optionally raised by any position it solves optimally, and the
upper bound 26 inherited from the quarter-turn metric, since every
quarter-turn sequence is also an r3 sequence.

#### Brute-force reference solver

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

### Demo

```bash
python demo.py                    # 30 random moves, then solve and render
python demo.py -n 50 --seed 7     # reproducible 50-move scramble
python demo.py --no-video         # solve and verify only
python demo.py -n 12 --optimal    # shortest possible answer
```

The demo scrambles the cube, solves it with the fast solver, verifies the
answer, and writes two videos (under a minute of rendering for a 30-move
scramble):

- `question.mp4` – the scramble applied to a solved cube.
- `answer.mp4` – the solution applied to the scrambled cube.

The videos come from `visual.py`: the camera stays fixed, two views show
opposite corners of the cube so every face is visible, and each move is
animated as a smooth 90° slice turn followed by a short hold. Use
`--frames-per-turn` to change the turn speed. Middle-slice rotations are
shown the way r3 defines them, as the two outer slices turning the opposite
way.

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
| `thistlethwaite.py` | Fast four-phase solver, `solve`                              |
| `optimal.py`      | Optimal IDA* solver with pattern databases, `solve`             |
| `gods_number.py`  | Bounds on God's number for the r3 move set                     |
| `cube_vision.py`  | Two photos → cube state → `RotationSequence` via an open-source VLM |
| `upload_server.py`| Phone-friendly upload page that runs the photo pipeline         |
| `COOKBOOK.md`     | Walkthrough for the photo pipeline                              |
| `examples/`       | Example photos and hand-read `reading.json` / `reading2.json` fixtures |
| `solver.py`       | `dfs` and the multiprocess `brute_force_multi` reference solver |
| `visual.py`       | Fixed-camera renderer with animated slice turns, `export_video` |
| `demo.py`         | Scramble, solve, verify, and export both videos                 |
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
