# Cookbook: from two photos to a `RotationSequence` with an open-source VLM

This guide turns two photographs of a physical cube into an r3 `RotationSequence`,
using an open-source vision-language model for the part that needs eyes and the
project's solvers for the part that needs planning. Everything runs on your own
hardware: the model on a DGX Spark on the LAN (Ollama at `192.168.1.254:11434`),
the rest on any machine with this repo. The whole pipeline is `cube_vision.py`.

```
photo ──► VLM reads 27 stickers (per photo) ──► JSON grids ──► geometry ──► colours at 54 positions
                 ▲                                                               │
                 │  re-read per face / re-read with the error / majority vote     ▼
                 └────────────────────────────────── reject ◄── validation ──► r3 state
                                                       │ (unique fix)              │
                                                       └── bounded repair ─────────┘
r3 state ──► thistlethwaite.solve / optimal.solve ──► RotationSequence (answer)
                                                  └─► inverse = scramble
```

Measured on the two example photos (54 stickers, hand-read ground truth in
`examples/reading.json`), one request per photo unless stated:

| Model / setting                           | Correct stickers | Time   |
| ----------------------------------------- | ---------------- | ------ |
| `qwen3.8:27b`, thinking on (**default**)  | **54 / 54**      | 100 s  |
| `qwen3.8:27b`, thinking off               | 50 / 54          | 10 s   |
| `qwen3.8:27b`, thinking off, per face     | 52 / 54          | 108 s  |
| `qwen3.8:27b`, thinking off, 5-sample vote| 51 / 54          | 42 s   |
| `qwen2.5vl:72b`                           | 43 / 54          | 101 s  |
| `qwen2.5vl:7b`                            | 30 / 54          | 13 s   |
| `qwen3.8:27b`, both photos in one request | fails (1 view)   |        |
| `qwen3.8:27b`, images at 2048 px          | 43 / 54          | 24 s   |

Two facts shaped the design: the model must see **one photo per request** (given
two, it returns one grid set), and its mistakes are **position slips within a
grid**, not colour confusion, which is why thinking mode and per-face requests
help and higher resolution does not.

A second, harder pair of photos (`examples/scramble2_*.jpg`, ground truth in
`examples/reading2.json`, one of them shot at a 45° tilt with the white face at
the bottom) needed the full ladder: round 1 misread 3 stickers, the per-face round
and one feedback re-read plus a single-sticker repair produced the correct state
after 13 model calls. The recovered answer matched the ground-truth answer exactly.

## 1. Why this shape

A vision-language-action model maps pixels to actions. For a cube, asking a model to
emit moves directly is fragile: the move sequence depends on all 54 stickers at once,
and a single misread sticker makes the whole sequence wrong with no way to notice.
Splitting the job fixes that:

- **Perception (the VLM):** read colours into a fixed JSON layout. Easy to check.
- **Geometry (code):** map each grid cell to an r3 sticker position. Deterministic.
- **Validation (code):** every colour nine times, every piece a real piece, no piece
  twice, the state physically possible. A failure produces a sentence the model can
  act on.
- **Recovery (code):** re-read per face, re-read with the error as feedback, majority
  vote over all readings, and a bounded repair that only accepts a *unique* legal fix.
- **Planning (code):** `thistlethwaite.solve` for an instant answer, `optimal.solve`
  for the shortest one. The inverse of the answer is the scramble that recreates the
  photographed cube from a solved one.

Together the pipeline behaves like a VLA: images in, actions out.

## 2. Serve the model

The defaults point at the Ollama server on the DGX Spark and the model
`qwen3.8:27b`. Ollama's native API is used because it exposes the `think` switch
and constrained JSON output (`format` = JSON schema).

```bash
ollama pull qwen3.8:27b
ollama serve                              # http://<host>:11434
curl http://192.168.1.254:11434/api/tags  # smoke test from the client machine
```

Thinking mode matters: without it the 27B model slips a row or column a few times
per photo; with it the example photos read perfectly. The first request after
`ollama serve` also loads the model, so it is slower.

## 3. Take the photos

Each photo is a **corner view**: the camera sits on a body diagonal of the cube and
looks at the centre, so exactly three faces are visible and they meet at the corner
nearest the camera, in the middle of the picture.

- **Photo 1:** camera at (2, 2, 2) looking toward the origin. With the project's
  colour scheme that shows yellow on top, blue lower-left, red lower-right.
- **Photo 2:** camera at (−2, −2, −2) looking toward the origin. In practice you turn
  the cube over. Recommended: white on top, orange lower-left, green lower-right.

The in-plane orientation does not matter at all: before reading, the pipeline asks
the model for the clock position of each face's centre sticker (a cheap question
without thinking) and rotates the photo so that one face sits at the top with the
other two near 4 and 8 o'clock. Only the clock numbers are used, so a wrong colour
in that answer is harmless; if the numbers do not form a corner-view pattern the
question is repeated with thinking, and a refinement pass runs on the rotated
image. The code then identifies each face from its centre sticker, and a
handedness check rejects impossible (mirror-image) views. The hard requirements are
that the two photos show complementary sets of faces and that all nine stickers of
every face are visible. Diffuse light, no flash, and a plain
background help; white stickers photograph grey, which the prompt anticipates.
Images are downscaled to 1024 px on the long side before sending; that read better
than 2048 px.

The example photos in `examples/` show both views.

## 4. What the model is asked for

Per photo, the model returns three 3×3 grids, indexed relative to the image, never
to the cube's colours:

```
                     [T00]
                [T10]     [T01]
           [T20]     [T11]     [T02]
                [T21]     [T12]
      [L00]          [T22]          [R02]
           [L01]  [L02] | [R00]  [R01]
      [L10]  [L11]  [L12] | [R10]  [R11]  [R12]
      [L20]  [L21]  [L22] | [R20]  [R21]  [R22]
```

| Grid          | Row index `i`            | Column index `j`                          | Touches the near corner |
| ------------- | ------------------------ | ----------------------------------------- | ----------------------- |
| `top`         | grows moving down-left   | grows moving down-right                   | `top[2][2]`             |
| `lower_left`  | top → bottom             | left → right, column 2 at the centre edge | `lower_left[0][2]`      |
| `lower_right` | top → bottom             | left → right, column 0 at the centre edge | `lower_right[0][0]`     |

The prompts are `cube_vision.VIEW_PROMPT` (all three grids) and
`cube_vision.FACE_PROMPT` (one grid); `VIEW_SCHEMA` / `FACE_SCHEMA` are sent as the
JSON schema for constrained decoding, so colour names are restricted to the six real
ones and the shape is always 3×3.

## 5. From grids to r3 coordinates

`visual.py` fixes the colour scheme: blue x+, green x−, red y+, orange y−,
yellow z+, white z−. The centre sticker of each grid therefore tells the code which
physical face it is: axis `a` and sign `s`.

Let the top, lower-left and lower-right faces be `(at, st)`, `(al, sl)`, `(ar, sr)`.
The view is a real (not mirrored) view exactly when

```
(at, al, ar) is an even permutation of (x, y, z)   ⇔   st · sl · sr = +1
```

Yellow/blue/red is `(z, x, y)`, even, with sign product +1: valid. Green on top with
white lower-left and orange lower-right, as in the second example photo, is
`(x, z, y)`, odd, with product −1: also valid. Swapping the two side faces breaks
the equivalence and the reading is rejected.

Each cell then maps to a cube cell `(x, y, z)` and a normal axis:

| Grid cell            | Along `at`    | Along `al`    | Along `ar`    | Normal |
| -------------------- | ------------- | ------------- | ------------- | ------ |
| `top[i][j]`          | `st`          | `sl · (i − 1)`| `sr · (j − 1)`| `at`   |
| `lower_left[i][j]`   | `st · (1 − i)`| `sl`          | `sr · (j − 1)`| `al`   |
| `lower_right[i][j]`  | `st · (1 − i)`| `sl · (1 − j)`| `sr`          | `ar`   |

and `r3.c2i[(x, y, z, normal)]` gives the sticker position index. This is
`cube_vision.view_positions`.

Positions only say which colour sits where. r3 wants a permutation: where each
*home* sticker went. `state_from_views` groups the 54 positions by cubie, looks up
the home cubie with that colour set, and records `state[home_sticker] = position`.

## 6. Validation, recovery, repair

Every check raises `CubeReadError` with a sentence meant for the model:

| Check                                   | Typical cause                                   |
| --------------------------------------- | ----------------------------------------------- |
| three different centre colours          | a side face was read twice                      |
| handedness                              | side centres swapped, grid filled mirrored      |
| all six faces seen, none twice          | both photos show the same corner                |
| each colour exactly nine times          | one sticker slipped to a neighbour's colour     |
| every piece is a real corner or edge    | a sticker from the wrong cell                   |
| no piece appears twice                  | two slips that happen to balance the counts     |
| state physically possible               | a piece twisted, flipped, or two pieces swapped |

The last check is `thistlethwaite.solve` itself, which raises for any state outside
the cube group.

`read_state` escalates through rounds until a candidate passes (every read uses
the rotation found by the orientation step):

1. **Round 1:** one request per photo, thinking on. On easy photos this is enough.
2. **Round 2:** one request per face (three per photo) without thinking: fast extra
   votes, and the model tracks one grid at a time.
3. **Round 3+:** re-read each photo with the rejection sentence appended, thinking
   on, at temperature 0.7 for diversity.

After each round two candidates are tried: the latest reading, and from round 2 on
the per-sticker **majority vote** over every reading so far. When a candidate is
rejected, `repair` looks for the smallest edit that makes it legal:

- **One sticker:** if colour counts are off by one, move an over-represented sticker
  to an under-represented colour, trying colours other readings proposed first.
- **One piece:** if counts are fine but the state is impossible, flip an edge or
  twist a corner in place. Legality alone cannot choose the piece (flipping *any*
  other edge also restores parity), so only pieces the readings disagree on are
  candidates.

A repair is accepted only when exactly one legal edit exists; it is printed as
`repaired one sticker (photo 1 top[2][2]: yellow -> red); please confirm against the
photo` and stored in the saved JSON under `log.repairs`. Anything ambiguous falls
through to the next round rather than guessing.

## 7. Solve and play back

```python
from cube_vision import Reader, read_state, solve_state

reader = Reader("http://192.168.1.254:11434", "qwen3.8:27b")   # think=True by default
state, views, log = read_state(["photo1.jpg", "photo2.jpg"], reader)
answer, scramble = solve_state(state)              # or optimal=True
```

- `answer` is a `RotationSequence` that solves the photographed cube. Apply it to
  the physical cube in order.
- `scramble` is its inverse: from a solved cube it reproduces exactly what the
  photos show, so `scramble(index)` equals `state`.
- `optimal=True` gives the shortest answer in r3 rotations; see the README for its
  running time.
- To watch the answer, `visual.export_video(coords[state], answer.seq, "answer.mp4")`.

## 8. From the phone

`upload_server.py` is a one-file web page (standard library only) that runs the
whole pipeline on wayne-kv:

```bash
python upload_server.py            # http://0.0.0.0:8080
```

On the phone, open `http://<wayne-kv address>:8080`, pick or shoot the two photos,
tap **Upload and solve**. The photos land in `uploads/<timestamp>/`, the page
refreshes every few seconds while the model reads, then shows the answer and
scramble sequences and plays `answer.mp4` inline. wayne-kv is reachable from an
iPhone with NordVPN Meshnet at its meshnet address (`nordlynx` interface,
`100.83.197.164` at the time of writing) or on the LAN at `192.168.1.20`.

## 9. Running it from the shell

Offline check without a model, using the example photos read by hand
(`reading.json` for the first pair, `reading2.json` for the tilted second pair):

```bash
python cube_vision.py --from-json examples/reading.json --video
python cube_vision.py --from-json examples/reading2.json
```

With the model on the DGX Spark (defaults: `http://192.168.1.254:11434`,
`qwen3.8:27b`, thinking on, four rounds):

```bash
python cube_vision.py examples/photo_yellow_blue_red.jpg examples/photo_green_white_orange.jpg \
    --save-json reading.json --video
```

Expected output on the examples: round 1 accepted after two model calls in about
100 s, a 48-move answer, `verified: True`, and `reading.json` identical to
`examples/reading.json`. Useful switches:

- `--no-think`: each call is ten times faster, but expect several rounds. On the
  examples it converged in round 4 after 12 calls (58 s) with one sticker repaired,
  and produced the same answer as the thinking run.
- `--model qwen2.5vl:72b` (or any model on the server) for a second opinion.
- `--optimal` for the shortest answer, `--attempts N` for more rounds.

Troubleshooting:

- **Handedness error on every attempt:** the photo is not a corner view, or the
  cube's colour scheme differs from the project's. Check `visual.hex`; a cube with
  a different scheme needs those six colours changed and nothing else.
- **Colour counts off in every round:** lighting. Retake with diffuse light and no
  flash; make sure every sticker of all three faces is inside the frame.
- **Impossible piece / physically impossible state that never resolves:** the model
  keeps slipping the same row. Crop the photo so the cube fills the frame, or try
  the 72B model; the majority vote across rounds usually settles it.
- **Empty or truncated replies:** the thinking budget ran out; raise `num_predict`
  in `Reader` (default 16000).
