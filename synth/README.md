# Sticker keypoint detector (branch `exp`)

A local replacement for the vision-language model in round 0 of the photo
pipeline. The network only *locates* the 27 stickers of a corner-view photo;
`cube_locate.py` still does the face split, grid orientation and pixel colour
classification, and `cube_vision.py` the validation, repair and solving. Nothing
here needs the DGX Spark: rendering, training and inference run on the CPU of
wayne-kv.

```
photo ──► KeypointNet (0.53 M params, 256 px, stride-4 heatmap) ──► 27 peaks
      ──► snap to sticker patch ──► lattice split into 3 faces ──► roles, colours ──► state
```

## Files

| File                 | Purpose                                                                   |
| -------------------- | ------------------------------------------------------------------------- |
| `synth/render.py`    | Synthetic corner-view photos with labelled sticker centres (domain randomisation) |
| `synth/model.py`     | The heatmap network, focal loss, peak extraction                          |
| `synth/train.py`     | CPU training on rendered data, evaluation on real photos, checkpoints     |
| `synth/real/*.json`  | Sticker-centre labels of real photos (normalised x, y), used for evaluation |
| `cube_keypoints.py`  | Inference: photo → points → `cube_locate.measure()`                       |
| `synth/ckpt/best.pt`, `alt.pt` | Trained weights (committed, 2 MB each; `synth/train.py` regenerates them) |

## Synthetic data

`render.py` is a small numpy/PIL rasteriser: 26 cubies as black blocks with
inset stickers, a pinhole camera placed near either body diagonal with up to
35° of tilt, random distance, field of view and roll, Lambert shading with a
coloured light plus a specular highlight, and a background that is a crop of the
real photos, a procedural wood texture, a gradient or a flat colour. The image is
then blurred, exposure- and white-balance-shifted, noised and JPEG-compressed.
Sticker colours are a random 9-per-colour permutation with per-face jitter. The
labels are the projected sticker centres, so they cost nothing. A 256 px sample
takes about 65 ms; `python -m synth.render` writes a 4×4 preview sheet.

## Training

```bash
.venv/bin/python -m synth.train --steps 8000 --batch 64 --workers 48 --threads 96
```

Data is rendered on the fly by DataLoader workers. Every 250 steps the model is
scored on the real photos in `synth/real/` (PCK: fraction of labelled centres
with a peak within 30 % of the sticker spacing) and `synth/ckpt/best.pt` keeps
the best checkpoint.

Throughput on the 96-core EPYC 9655P (192 logical CPUs):

| Setting                                  | samples/s |
| ---------------------------------------- | --------- |
| fp32, 96 threads                         | 46        |
| fp32, 192 threads                        | 30        |
| bf16 autocast + channels_last, 96 threads| 361       |

Threads beyond the 96 physical cores do not help (the workload is bandwidth
bound, hyper-threads only compete), so a 50 % reading in `top` is the machine
fully used. bf16 uses the CPU's avx512_bf16 units and is the real speed-up.

Real-photo PCK reaches 0.95 after 250 steps and plateaus at 0.98-0.99 from
step 500 on (8000 steps, batch 64, about 50 minutes). `best.pt` is the step-6250
checkpoint (PCK 0.989); `alt.pt` is the step-500 checkpoint, used as a second
opinion when the three-face split fails with the first, since the two make
different localisation slips on the same hard photo. Both are committed (2 MB
each).

## Using it

```bash
.venv/bin/python cube_vision.py A.jpg B.jpg --detector cnn        # no Ollama call
```

or in Python:

```python
from cube_keypoints import read_photos
from cube_vision import state_from_views, solve_state
state = state_from_views(read_photos(["A.jpg", "B.jpg"]))
answer, scramble = solve_state(state)
```

Two photos take 2-7 s in total, most of it the geometry and colour code, not the
network. Results on the three real pairs (`best.pt` with the `alt.pt` fallback):

| Pair                          | CNN result                                |
| ----------------------------- | ----------------------------------------- |
| first example                 | state equals ground truth                 |
| tilted second example         | state equals ground truth                 |
| phone pair                    | legal state, identical to the VLM result  |

## Dependencies

Only `torch` (CPU wheel) on top of the main requirements, installed into `.venv`:

```bash
.venv/bin/pip install torch --index-url https://download.pytorch.org/whl/cpu
```

## What could still go wrong

- A sticker completely missed by the network: one missing sticker per face is
  reconstructed from the lattice; two on one face are not.
- A spurious peak with a high score: the top-27 peaks are tried first and a
  few extras are admitted only if the three-face split fails.
- New cube colour schemes or unusual lighting: the colour classifier adapts
  per photo, the renderer can be widened if a real failure shows up.
