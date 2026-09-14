"""Two corner-view photos -> cube state -> RotationSequence.

An open-source vision-language model (default qwen3.8:27b on the Ollama server at
192.168.1.254) does the *perception*: it reads the 54 sticker colours into JSON.
Deterministic code does the *geometry* (image grid -> r3 coordinates), validates
the reading, repairs or re-reads when it cannot be a real cube, and the existing
solvers do the *planning* (state -> RotationSequence). See COOKBOOK.md.

Photo convention ("corner view"): the camera looks along a body diagonal, e.g. from
(2,2,2) toward the origin and then from (-2,-2,-2) toward the origin. Three faces are
visible and meet at the corner nearest the camera, in the image centre:

                     [T00]                    T = top face (a rhombus)
                [T10]     [T01]               L = lower-left face
           [T20]     [T11]     [T02]          R = lower-right face
                [T21]     [T12]
      [L00]          [T22]          [R02]     T22, L02, R00 touch the near corner
           [L01]  [L02] | [R00]  [R01]        '|' is the centre vertical edge
      [L10]  [L11]  [L12] | [R10]  [R11]  [R12]
      [L20]  [L21]  [L22] | [R20]  [R21]  [R22]

    top[i][j]         i grows moving down-left, j grows moving down-right
    lower_left[i][j]  i = row top->bottom, j = column left->right (col 2 at the edge)
    lower_right[i][j] i = row top->bottom, j = column left->right (col 0 at the edge)

Which physical face is which is detected from the three centre stickers, so the
cube may be turned in-plane however you like; a mirror-image (impossible) view is
rejected by a handedness check.

Reading strategy (measured on the example photos, see COOKBOOK.md):
  round 1  one request per photo, thinking on                 -> 54/54 stickers
  round 2  one request per face                               -> 52/54 without thinking
  round 3+ re-read with the validation error as feedback and majority-vote all
           readings so far
  after every round a bounded repair tries single-sticker / single-piece fixes
  that make the reading a legal cube, and reports them.

Usage::

    python cube_vision.py A.jpg B.jpg                 # Ollama on the DGX Spark by default
    python cube_vision.py --from-json examples/reading.json --video
"""
import argparse
import base64
import io
import json
import re
import time
import urllib.error
import urllib.request
from collections import Counter

import numpy as np
from PIL import Image

from r3 import index, coords, c2i, b2i, rs, RotationSequence
from visual import hex as FACE_HEX
import thistlethwaite

# ----------------------------------------------------------------------------
# Colour scheme: taken from visual.py so the whole project agrees on it
# ----------------------------------------------------------------------------
_NAMES = {"#3F51B5": "blue", "#4CAF50": "green", "#B51F1F": "red",
          "#FF6F00": "orange", "#FFEB3B": "yellow", "#FFFFFF": "white"}
FACE_OF_COLOR = {_NAMES[h]: face for face, h in FACE_HEX.items()}     # name -> (axis, sign)
COLOR_OF_FACE = {face: name for name, face in FACE_OF_COLOR.items()}  # (axis, sign) -> name
COLORS = sorted(FACE_OF_COLOR)
AXIS_NAME = "xyz"
GRIDS = ("top", "lower_left", "lower_right")
CELLS = [(g, i, j) for g in GRIDS for i in range(3) for j in range(3)]

# home colour of every sticker, and colour-set -> home cubie
HOME_COLOR = {f: COLOR_OF_FACE[(int(c[3]), int(c[c[3]]))] for f, c in enumerate(coords)}
HOME_BLOCK = {frozenset(HOME_COLOR[f] for f in stickers): block
              for block, stickers in b2i.items() if len(stickers) > 1}
ROT = {str(r): r for r in rs}


def face_name(face):
    axis, sign = face
    return f"{COLOR_OF_FACE[face]} ({AXIS_NAME[axis]}{'+' if sign > 0 else '-'})"


# ----------------------------------------------------------------------------
# Prompts and JSON schemas for the vision-language model
# ----------------------------------------------------------------------------
GRID_SCHEMA = {"type": "array", "minItems": 3, "maxItems": 3,
               "items": {"type": "array", "minItems": 3, "maxItems": 3,
                         "items": {"type": "string", "enum": COLORS}}}
VIEW_SCHEMA = {"type": "object", "properties": {g: GRID_SCHEMA for g in GRIDS},
               "required": list(GRIDS), "additionalProperties": False}
FACE_SCHEMA = {"type": "object", "properties": {"grid": GRID_SCHEMA},
               "required": ["grid"], "additionalProperties": False}

GEOMETRY = """You are reading the stickers of a 3x3x3 Rubik's cube from a photograph.

The photo is a "corner view": the camera looks at the cube along a body diagonal, so
exactly three faces are visible and they meet at the corner closest to the camera,
near the centre of the image. On screen the three faces are:
- TOP: the face on top, drawn as a rhombus (diamond shape).
- LOWER_LEFT: the face on the lower left, a parallelogram.
- LOWER_RIGHT: the face on the lower right, a parallelogram.

Each face is a 3x3 grid of colour names, indexed relative to the image:
- top[i][j]: top[0][0] is the sticker at the TOPMOST vertex of the rhombus. The index
  i increases as you move DOWN-LEFT (toward the lower-left face); the index j increases
  as you move DOWN-RIGHT (toward the lower-right face). So top[1][1] is the centre,
  top[2][2] is the sticker touching the near corner, top[2][0] is the left vertex and
  top[0][2] is the right vertex.
- lower_left[i][j]: i is the row from top (touching the top face) to bottom; j is the
  column from left to right, so column 2 touches the vertical edge in the middle of
  the picture and lower_left[0][2] touches the near corner.
- lower_right[i][j]: i is the row from top to bottom; j is the column from left to
  right, so column 0 touches the vertical edge in the middle of the picture and
  lower_right[0][0] touches the near corner.

Sketch (each label is one sticker, T = top, L = lower_left, R = lower_right):

                     [T00]
                [T10]     [T01]
           [T20]     [T11]     [T02]
                [T21]     [T12]
      [L00]          [T22]          [R02]
           [L01]  [L02] | [R00]  [R01]
      [L10]  [L11]  [L12] | [R10]  [R11]  [R12]
      [L20]  [L21]  [L22] | [R20]  [R21]  [R22]

Colours: exactly these six names: white, yellow, red, orange, blue, green.
- White stickers often look grey in photos; report them as white.
- Orange is lighter and more yellowish than red; look at several stickers to calibrate.
- The centre sticker of a face (index [1][1]) may carry a logo; report its colour.
"""

VIEW_PROMPT = GEOMETRY + """
Work face by face and row by row: first identify the three centre colours, then fill
every grid, checking each sticker's neighbours so that no row or column is shifted.

Output ONLY a JSON object of the form
{"top": [[...],[...],[...]], "lower_left": [[...],[...],[...]], "lower_right": [[...],[...],[...]]}"""

FACE_PROMPT = GEOMETRY + """
Report ONLY the {face} face of this photo, as a JSON object {{"grid": [[...],[...],[...]]}}
using the {face}[i][j] indexing defined above. Use the other two faces only to orient
yourself."""

FEEDBACK = """

A previous reading of this photo was rejected for this reason:
{error}
Look at the photo again, sticker by sticker, and output the corrected JSON."""


# ----------------------------------------------------------------------------
# Talking to the model
# ----------------------------------------------------------------------------
def encode_image(path, max_side=1024):
    """Downscale to max_side px and return JPEG bytes; 1024 read best in tests."""
    img = Image.open(path).convert("RGB")
    scale = max_side / max(img.size)
    if scale < 1:
        img = img.resize((round(img.width * scale), round(img.height * scale)), Image.LANCZOS)
    buf = io.BytesIO()
    img.save(buf, format="JPEG", quality=92)
    return buf.getvalue()


def _post(url, payload, api_key=None, timeout=1800):
    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def extract_json(text):
    """Strip <think> blocks and parse the first {...} object in the reply."""
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.S)
    start, end = text.find("{"), text.rfind("}")
    if start < 0 or end < 0:
        raise CubeReadError("the reply contained no JSON object")
    try:
        return json.loads(text[start:end + 1])
    except json.JSONDecodeError as e:
        raise CubeReadError(f"the reply was not valid JSON: {e}")


class Reader:
    """One model behind either Ollama's native API or an OpenAI-compatible endpoint."""

    def __init__(self, base_url="http://192.168.1.254:11434", model="qwen3.8:27b", backend="auto",
                 think=True, api_key=None, max_side=1024, num_ctx=16384, num_predict=16000, verbose=True):
        self.base_url = base_url.rstrip("/")
        if backend == "auto":
            backend = "ollama" if ":11434" in self.base_url or not self.base_url.endswith("/v1") else "openai"
        self.backend, self.model, self.think, self.api_key = backend, model, think, api_key
        self.max_side, self.num_ctx, self.num_predict, self.verbose = max_side, num_ctx, num_predict, verbose
        self.calls = 0

    def chat(self, image, prompt, schema, temperature=0.0):
        """Send one image + prompt, return the parsed JSON object."""
        self.calls += 1
        t0 = time.time()
        data = encode_image(image, self.max_side)
        if self.backend == "ollama":
            url = (self.base_url[:-3] if self.base_url.endswith("/v1") else self.base_url) + "/api/chat"
            payload = {"model": self.model, "stream": False, "think": self.think, "format": schema,
                       "options": {"temperature": temperature, "num_ctx": self.num_ctx,
                                   "num_predict": self.num_predict},
                       "messages": [{"role": "user", "content": prompt,
                                     "images": [base64.b64encode(data).decode()]}]}
            reply = _post(url, payload)["message"]["content"]
        else:
            url = self.base_url + "/chat/completions"
            content = [{"type": "image_url", "image_url": {
                "url": "data:image/jpeg;base64," + base64.b64encode(data).decode()}},
                {"type": "text", "text": prompt}]
            payload = {"model": self.model, "temperature": temperature, "max_tokens": self.num_predict,
                       "messages": [{"role": "user", "content": content}]}
            fmt = {"type": "json_schema", "json_schema": {"name": "cube", "schema": schema}}
            try:
                resp = _post(url, dict(payload, response_format=fmt), self.api_key)
            except urllib.error.HTTPError as e:
                if e.code != 400:
                    raise
                resp = _post(url, payload, self.api_key)
            reply = resp["choices"][0]["message"]["content"] or ""
            if isinstance(reply, list):
                reply = "".join(part.get("text", "") for part in reply)
        if self.verbose:
            print(f"    {self.model} replied in {time.time() - t0:.0f}s")
        return extract_json(reply)

    def read_view(self, image, feedback=None, temperature=0.0):
        """All three grids of one photo in one request."""
        view = self.chat(image, VIEW_PROMPT + (FEEDBACK.format(error=feedback) if feedback else ""),
                         VIEW_SCHEMA, temperature)
        return {g: view.get(g) for g in GRIDS}

    def read_view_by_faces(self, image, temperature=0.0):
        """One request per face; the model tracks a single grid at a time."""
        return {g: self.chat(image, FACE_PROMPT.format(face=g.upper()), FACE_SCHEMA, temperature).get("grid")
                for g in GRIDS}


# ----------------------------------------------------------------------------
# Geometry: image-relative grids -> r3 sticker positions -> r3 state
# ----------------------------------------------------------------------------
class CubeReadError(Exception):
    """A reading that cannot be a real cube; the message is written for the model."""


def _check_grid(view, grid, k):
    g = view.get(grid) if isinstance(view, dict) else None
    ok = (isinstance(g, list) and len(g) == 3
          and all(isinstance(r, list) and len(r) == 3 for r in g)
          and all(c in FACE_OF_COLOR for r in g for c in r))
    if not ok:
        raise CubeReadError(f"photo {k}: '{grid}' must be a 3x3 grid of {COLORS}")
    return g


def view_faces(view, k=1):
    """Physical face of each grid, from the centre stickers, with a handedness check."""
    faces = {g: FACE_OF_COLOR[_check_grid(view, g, k)[1][1]] for g in GRIDS}
    if len(set(faces.values())) < 3:
        raise CubeReadError(f"photo {k}: the three centre stickers must be three different colours, "
                            f"got {[view[g][1][1] for g in GRIDS]}")
    (at, st), (al, sl), (ar, sr) = (faces[g] for g in GRIDS)
    even = (at, al, ar) in {(0, 1, 2), (1, 2, 0), (2, 0, 1)}
    if even != (st * sl * sr == 1):
        raise CubeReadError(
            f"photo {k}: {view['top'][1][1]} on top, {view['lower_left'][1][1]} lower-left and "
            f"{view['lower_right'][1][1]} lower-right is a mirror image of a real cube. The two side "
            f"centres are probably swapped, or the grids were filled in the wrong orientation.")
    return faces


def view_positions(faces):
    """{(grid, i, j): sticker position index} for a corner view with the given faces."""
    (at, st), (al, sl), (ar, sr) = (faces[g] for g in GRIDS)
    pos = {}
    for i in range(3):
        for j in range(3):
            c = [0, 0, 0]
            c[at], c[al], c[ar] = st, sl * (i - 1), sr * (j - 1)
            pos[("top", i, j)] = c2i[(c[0], c[1], c[2], at)]
            c = [0, 0, 0]
            c[al], c[at], c[ar] = sl, st * (1 - i), sr * (j - 1)
            pos[("lower_left", i, j)] = c2i[(c[0], c[1], c[2], al)]
            c = [0, 0, 0]
            c[ar], c[at], c[al] = sr, st * (1 - i), sl * (1 - j)
            pos[("lower_right", i, j)] = c2i[(c[0], c[1], c[2], ar)]
    return pos


def state_from_views(views, verbose=False):
    """Turn the grids into an r3 index vector, or raise CubeReadError."""
    if not isinstance(views, list) or not views:
        raise CubeReadError("'views' must be a non-empty list, one entry per photo")
    color_at = {}
    for k, view in enumerate(views, 1):
        faces = view_faces(view, k)
        if verbose:
            print(f"  photo {k}: top {face_name(faces['top'])}, lower-left {face_name(faces['lower_left'])}, "
                  f"lower-right {face_name(faces['lower_right'])}")
        for (grid, i, j), f in view_positions(faces).items():
            color = view[grid][i][j]
            if color_at.get(f, color) != color:
                raise CubeReadError(f"photo {k}: the {view[grid][1][1]} face was already read from another "
                                    f"photo with different colours; each face must appear in one photo only")
            color_at[f] = color
    if len(color_at) < len(coords):
        seen = {(int(coords[f][3]), int(coords[f][coords[f][3]])) for f in color_at}
        missing = [face_name(face) for face in COLOR_OF_FACE if face not in seen]
        raise CubeReadError(f"these faces were not seen in any photo: {missing}")
    counts = Counter(color_at.values())
    if any(counts[c] != 9 for c in COLORS):
        raise CubeReadError("each colour must appear exactly 9 times, got " +
                            ", ".join(f"{c}: {counts[c]}" for c in COLORS))
    state = index.copy()
    used = {}
    for block, stickers in b2i.items():
        if len(stickers) == 1:                      # centres never move
            continue
        colors = [color_at[f] for f in stickers]
        key = frozenset(colors)
        if len(key) != len(colors) or key not in HOME_BLOCK:
            raise CubeReadError(f"the piece at cell {block} shows colours {colors}, but no real "
                                f"{'corner' if len(colors) == 3 else 'edge'} piece has that combination")
        if key in used:
            raise CubeReadError(f"the cells {used[key]} and {block} both show the piece with colours "
                                f"{sorted(key)}; a real cube has only one such piece")
        used[key] = block
        for f, color in zip(stickers, colors):
            home = next(h for h in b2i[HOME_BLOCK[key]] if HOME_COLOR[h] == color)
            state[home] = f
    assert sorted(state.tolist()) == list(range(len(coords)))
    try:
        thistlethwaite.solve(state, verbose=False)
    except (ValueError, RuntimeError):
        raise CubeReadError("every piece exists but the arrangement is physically impossible (a piece is "
                            "flipped or twisted, or two pieces are swapped). At least one sticker is "
                            "misread; re-check the edge and corner pieces sticker by sticker")
    return state


# ----------------------------------------------------------------------------
# Robustness: majority vote over readings, and bounded legality-guided repair
# ----------------------------------------------------------------------------
def vote(pool):
    """Per-sticker majority over several readings of one photo (ties -> latest)."""
    out = {}
    for g in GRIDS:
        out[g] = [[Counter(v[g][i][j] for v in reversed(pool)).most_common(1)[0][0]
                   for j in range(3)] for i in range(3)]
    return out


def _legal(views):
    try:
        return state_from_views(views)
    except CubeReadError:
        return None


def _with(views, k, cell, color):
    new = [{g: [row[:] for row in v[g]] for g in GRIDS} for v in views]
    g, i, j = cell
    new[k][g][i][j] = color
    return new


def repair(views, pools=None, verbose=True):
    """Try the smallest edits that make the reading a legal cube.

    1. one sticker: if colour counts are off, move one over-represented sticker
       to an under-represented colour (alternatives seen in the vote pool first);
    2. one piece: flip an edge or twist a corner in place, only among pieces the
       vote pool disagrees on (needs pools).
    Accepts an edit only if it is the unique legal one, and reports it.
    Returns (state, views, edits) or (None, views, [])."""
    counts = Counter(v[g][i][j] for v in views for (g, i, j) in CELLS)
    over = [c for c in COLORS if counts[c] > 9]
    under = [c for c in COLORS if counts[c] < 9]
    if over and under:
        found = []
        for k, v in enumerate(views):
            for cell in CELLS:
                g, i, j = cell
                if v[g][i][j] not in over or (i, j) == (1, 1):
                    continue
                seen = {p[g][i][j] for p in (pools[k] if pools else [])}
                for color in sorted(under, key=lambda c: c not in seen):
                    cand = _with(views, k, cell, color)
                    if _legal(cand) is not None:
                        found.append((k, cell, v[g][i][j], color, cand))
        if len(found) == 1:
            k, (g, i, j), old, new, cand = found[0]
            edit = f"photo {k + 1} {g}[{i}][{j}]: {old} -> {new}"
            if verbose:
                print(f"  repaired one sticker ({edit}); please confirm against the photo")
            return _legal(cand), cand, [edit]
        if verbose and found:
            print(f"  {len(found)} different single-sticker fixes would be legal; not guessing")
        return None, views, []
    if not over and not under and pools:
        # Every colour 9 times but illegal: flip/twist one piece in place. Legality
        # alone cannot pick the piece (flipping any other edge also restores parity),
        # so only pieces whose stickers the readings disagree on are candidates.
        found = []
        pos = [(k, view_positions(view_faces(v, k + 1))) for k, v in enumerate(views)]
        cell_of = {f: (k, cell) for k, p in pos for cell, f in p.items()}
        for block, stickers in b2i.items():
            if len(stickers) == 1:
                continue
            cells = [cell_of[f] for f in stickers]
            if all(len({p[g][i][j] for p in pools[k]}) == 1 for k, (g, i, j) in cells):
                continue
            colors = [views[k][g][i][j] for k, (g, i, j) in cells]
            for shift in range(1, len(stickers)):
                cand = views
                for (k, cell), color in zip(cells, colors[shift:] + colors[:shift]):
                    cand = _with(cand, k, cell, color)
                if _legal(cand) is not None:
                    found.append((block, cand, colors, shift))
        if len(found) == 1:
            block, cand, colors, shift = found[0]
            edit = f"piece at {block}: rotated its colours {colors} by {shift}"
            if verbose:
                print(f"  repaired one piece ({edit}); please confirm against the photo")
            return _legal(cand), cand, [edit]
        if verbose and found:
            print(f"  {len(found)} different single-piece fixes would be legal; not guessing")
    return None, views, []


def read_state(image_paths, reader, attempts=4, verbose=True):
    """Read photos with escalating strategies until the reading is a legal cube.

    Returns (state, views, log) where log records rounds, errors and repairs."""
    pools = [[] for _ in image_paths]
    log = {"rounds": [], "repairs": []}
    error = None
    for rnd in range(1, attempts + 1):
        if verbose:
            how = {1: "one request per photo", 2: "one request per face"}.get(rnd, "re-read with feedback")
            print(f"round {rnd}: {how}")
        for k, path in enumerate(image_paths):
            if rnd == 1:
                view = reader.read_view(path)
            elif rnd == 2:
                view = reader.read_view_by_faces(path)
            else:
                view = reader.read_view(path, feedback=error, temperature=0.7)
            pools[k].append(view)
        candidates = [("latest reading", [p[-1] for p in pools])]
        if rnd > 1:
            candidates.append(("majority vote", [vote(p) for p in pools]))
        for name, views in candidates:
            try:
                state = state_from_views(views, verbose=verbose)
                log["rounds"].append({"round": rnd, "accepted": name})
                return state, views, log
            except CubeReadError as e:
                error = str(e)
                log["rounds"].append({"round": rnd, "candidate": name, "error": error})
                if verbose:
                    print(f"  {name} rejected: {error}")
            state, fixed, edits = repair(views, pools, verbose=verbose)
            if state is not None:
                log["repairs"] = edits
                log["rounds"].append({"round": rnd, "accepted": f"{name} + repair"})
                return state, fixed, log
    raise CubeReadError(f"no legal reading after {attempts} rounds; last error: {error}")


# ----------------------------------------------------------------------------
# Planning: state -> RotationSequence
# ----------------------------------------------------------------------------
def inverse(seq):
    """Rotations that undo seq: reversed order, each orientation flipped."""
    return [ROT[str(r)[:2] + {"p": "n", "n": "p"}[str(r)[2]]] for r in reversed(list(seq))]


def solve_state(state, optimal=False):
    """(answer, scramble): answer solves the photographed cube, scramble recreates it."""
    if optimal:
        import optimal as opt
        ans = opt.solve(state)
    else:
        ans = thistlethwaite.solve(state)
    answer = RotationSequence(ans)
    assert np.all(answer(state) == index)
    return answer, RotationSequence(inverse(ans))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("images", nargs="*", help="corner-view photos, all six faces between them")
    parser.add_argument("--base-url", default="http://192.168.1.254:11434",
                        help="Ollama server, or an OpenAI-compatible endpoint ending in /v1 (vLLM)")
    parser.add_argument("--model", default="qwen3.8:27b", help="model name on the server")
    parser.add_argument("--backend", choices=["auto", "ollama", "openai"], default="auto")
    parser.add_argument("--no-think", action="store_true", help="disable thinking (faster, less accurate)")
    parser.add_argument("--api-key", default=None)
    parser.add_argument("--attempts", type=int, default=4, help="reading rounds before giving up")
    parser.add_argument("--save-json", help="write the accepted reading (and log) here")
    parser.add_argument("--from-json", help="skip the model and use a saved reading")
    parser.add_argument("--optimal", action="store_true", help="shortest answer (optimal.py) instead of Thistlethwaite")
    parser.add_argument("--video", action="store_true", help="render answer.mp4 with visual2")
    args = parser.parse_args()

    log = {}
    if args.from_json:
        with open(args.from_json) as fh:
            views = json.load(fh)["views"]
        try:
            state = state_from_views(views, verbose=True)
        except CubeReadError as e:
            print(f"reading rejected: {e}")
            state, views, edits = repair(views)
            if state is None:
                raise SystemExit("could not repair the reading")
            log["repairs"] = edits
    elif args.images:
        reader = Reader(args.base_url, args.model, args.backend, think=not args.no_think, api_key=args.api_key)
        t0 = time.time()
        state, views, log = read_state(args.images, reader, args.attempts)
        print(f"reading accepted after {reader.calls} model calls, {time.time() - t0:.0f}s")
    else:
        parser.error("give photo paths or --from-json")
    if args.save_json:
        with open(args.save_json, "w") as fh:
            json.dump({"views": views, "log": log}, fh, indent=1)

    answer, scramble = solve_state(state, optimal=args.optimal)
    print(f"Answer   ({len(answer)} moves): {answer}")
    print(f"Scramble ({len(scramble)} moves): {scramble}")
    print("verified:", bool(np.all(answer(state) == index)))
    if args.video:
        from visual2 import export_video
        export_video(coords[state], answer.seq, "answer.mp4")
        print("Exported answer.mp4")


if __name__ == "__main__":
    main()
