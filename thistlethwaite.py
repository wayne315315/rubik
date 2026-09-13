"""Thistlethwaite four-phase solver built on the r3 coordinate model.

Axis roles (fixed by r3's coordinate system):
    UD = z axis (2), FB = x axis (0), RL = y axis (1)

The cube group is walked down through a chain of nested subgroups:

    phase 1  G0 -> G1   orient all 12 edges                 any move
    phase 2  G1 -> G2   orient all 8 corners and bring the   z, y quarter turns
                        E-slice edges into the E slice       + all half turns
    phase 3  G2 -> G3   M-slice edges into the M slice and   z quarter turns
                        corners into a half-turn-reachable   + all half turns
                        arrangement
    phase 4  G3 -> I    finish                               half turns only

(z0 / y0 middle-slice rotations from r3 are allowed wherever their two
outer-slice equivalents are.)

Every phase works on a small coordinate that fully captures its sub-problem,
and a complete breadth-first distance table is built over that coordinate
space.  Solving a phase is then a greedy descent along decreasing distances,
so no search is needed at solve time.  All move effects are derived from the
r3 ``Rotation`` objects themselves; nothing about the cube is hard-coded here
beyond which axis plays which role.

Usage::

    from thistlethwaite import solve
    ans = solve(index_q)        # index_q: a 54-vector produced by r3 rotations
    # ans is a list of r3 Rotation objects; RotationSequence(ans)(index_q) == index
"""
import os
import pickle
import time

import numpy as np

from r3 import coords, index, b2i, i2b, c2i, rs

UD, FB, RL = 2, 0, 1

# ----------------------------------------------------------------------------
# Cubies derived from r3's 54 stickers
# ----------------------------------------------------------------------------
CORNERS = sorted(b for b in b2i if all(v != 0 for v in b))
EDGES = sorted(b for b in b2i if sum(v == 0 for v in b) == 1)
CUBIES = CORNERS + EDGES                        # 8 corners then 12 edges
cubie_id = {b: i for i, b in enumerate(CUBIES)}
cubie_of_facelet = np.array([cubie_id.get(i2b[f], -1) for f in range(len(coords))])
facelet_of_cubie = np.array([b2i[b][0] for b in CUBIES])   # one sticker per cubie


def _facelet(block, normal):
    return next(f for f in b2i[block] if coords[f][3] == normal)


# Corner orientation is tracked through the UD-coloured sticker of each corner.
CORNER_UD_FACELETS = np.array([_facelet(b, UD) for b in CORNERS])
# Edge orientation is tracked through a reference sticker per edge:
# the UD sticker for UD-layer edges, the FB sticker for E-slice edges.
EDGE_REF_FACELETS = np.array([_facelet(b, UD if b[UD] != 0 else FB) for b in EDGES])
CORNER_IDS = np.array([cubie_id[b] for b in CORNERS])
EDGE_IDS = np.array([cubie_id[b] for b in EDGES])
E_SLICE_IDS = np.array([cubie_id[b] for b in EDGES if b[UD] == 0])
M_SLICE_IDS = np.array([cubie_id[b] for b in EDGES if b[UD] != 0 and b[RL] == 0])
S_SLICE_IDS = np.array([cubie_id[b] for b in EDGES if b[UD] != 0 and b[FB] == 0])


# ----------------------------------------------------------------------------
# Moves: one or two r3 rotations, with their action on stickers and cubies
# ----------------------------------------------------------------------------
class Move:
    def __init__(self, name, rots):
        self.name = name
        self.rots = tuple(rots)
        perm = index.copy()
        for r in self.rots:
            perm = r(perm)
        self.fperm = perm                                       # sticker -> sticker
        self.cperm = cubie_of_facelet[perm[facelet_of_cubie]]   # cubie  -> cubie

    def __repr__(self):
        return self.name


ROT = {str(r): r for r in rs}
MOVES = {}
for _axis in "xyz":
    for _level in "pn":
        for _orient in "pn":                       # quarter turns, e.g. xpp
            _n = _axis + _level + _orient
            MOVES[_n] = Move(_n, (ROT[_n],))
        _n = _axis + _level + "2"                  # half turns, e.g. xp2
        MOVES[_n] = Move(_n, (ROT[_axis + _level + "p"],) * 2)
    for _orient in "pn":                           # r3 middle-slice moves, e.g. x0p
        _n = _axis + "0" + _orient
        MOVES[_n] = Move(_n, (ROT[_n],))

def _quarter(axis):
    return [n for n in MOVES if n[0] == axis and n[1] in "pn" and n[2] in "pn"]


HALF = [n for n in MOVES if n.endswith("2")]
ALL_MOVES = list(MOVES)                                            # G0 generators
G1_MOVES = _quarter("z") + _quarter("y") + ["z0p", "z0n", "y0p", "y0n"] + HALF   # U D R L F2 B2
G2_MOVES = _quarter("z") + ["z0p", "z0n"] + HALF                   # U D R2 L2 F2 B2
G3_MOVES = HALF                                                    # U2 D2 R2 L2 F2 B2


# ----------------------------------------------------------------------------
# Coordinates: each is (extract(state) -> value, apply(value, move) -> value)
# A state is r3's 54-vector: state[f] = current sticker position of sticker f.
# ----------------------------------------------------------------------------
def facelet_set(facelets):
    """Unordered set of positions currently occupied by the given stickers."""
    def extract(state):
        return tuple(sorted(state[facelets]))

    def apply(value, move):
        return tuple(sorted(move.fperm[list(value)]))
    return extract, apply


def cubie_set(cubies):
    """Unordered set of cubie positions currently occupied by the given cubies."""
    facelets = facelet_of_cubie[cubies]

    def extract(state):
        return tuple(sorted(cubie_of_facelet[state[facelets]]))

    def apply(value, move):
        return tuple(sorted(move.cperm[list(value)]))
    return extract, apply


def cubie_perm(cubies):
    """Ordered tuple: where each of the given cubies currently is."""
    facelets = facelet_of_cubie[cubies]

    def extract(state):
        return tuple(cubie_of_facelet[state[facelets]])

    def apply(value, move):
        return tuple(move.cperm[list(value)])
    return extract, apply


def enumerate_coord(start, moves, apply):
    """BFS over one coordinate. Returns (value -> idx, move table[idx, move])."""
    states = [start]
    idx = {start: 0}
    rows = []
    i = 0
    while i < len(states):
        row = []
        for m in moves:
            t = apply(states[i], m)
            if t not in idx:
                idx[t] = len(states)
                states.append(t)
            row.append(idx[t])
        rows.append(row)
        i += 1
    return idx, np.array(rows, dtype=np.int32)


# ----------------------------------------------------------------------------
# Product of coordinates: mixed-radix index, vectorised neighbours and BFS
# ----------------------------------------------------------------------------
def product_strides(sizes):
    strides, s = [], 1
    for n in sizes:
        strides.append(s)
        s *= n
    return strides, s


def product_next(tables, sizes, strides, idxs, m):
    """Combined index of every state in idxs after move column m."""
    out = np.zeros_like(idxs)
    for table, size, stride in zip(tables, sizes, strides):
        out += table[(idxs // stride) % size, m] * stride
    return out


def product_bfs(tables, sizes, strides, N, sources, cols):
    """Exact distance (int8, -1 = unreachable) from every state to sources."""
    dist = np.full(N, -1, dtype=np.int8)
    dist[sources] = 0
    frontier = np.asarray(sources)
    d = 0
    while len(frontier):
        d += 1
        for m in cols:
            n = product_next(tables, sizes, strides, frontier, m)
            n = n[dist[n] < 0]
            dist[n] = d
        frontier = np.flatnonzero(dist == d)
    return dist


# ----------------------------------------------------------------------------
# A phase = product of coordinates + full BFS distance table to its goal set
# ----------------------------------------------------------------------------
class Phase:
    def __init__(self, name, move_names, specs, goal_move_names=()):
        self.name = name
        self.move_names = list(move_names)
        self.moves = [MOVES[n] for n in self.move_names]
        self.specs = specs
        self.goal_cols = [i for i, n in enumerate(self.move_names) if n in goal_move_names]

    # -- building -------------------------------------------------------
    def build(self):
        self.index = []
        self.tables = []
        for extract, apply in self.specs:
            idx, table = enumerate_coord(extract(index), self.moves, apply)
            self.index.append(idx)
            self.tables.append(table)
        self._finish()
        goal = self._reach(np.array([0]), self.goal_cols) if self.goal_cols else np.array([0])
        self.dist = self._bfs(goal, range(len(self.moves)))
        return self

    def _finish(self):
        self.sizes = [len(i) for i in self.index]
        self.strides, self.N = product_strides(self.sizes)

    def _next(self, idxs, m):
        return product_next(self.tables, self.sizes, self.strides, idxs, m)

    def _reach(self, start, cols):
        dist = product_bfs(self.tables, self.sizes, self.strides, self.N, start, cols)
        return np.flatnonzero(dist >= 0)

    def _bfs(self, goal, cols):
        return product_bfs(self.tables, self.sizes, self.strides, self.N, goal, cols)

    # -- (de)serialisation ------------------------------------------------
    def data(self):
        return {"index": self.index, "tables": self.tables, "dist": self.dist}

    def load(self, data):
        self.index, self.tables, self.dist = data["index"], data["tables"], data["dist"]
        self._finish()
        return self

    # -- solving ---------------------------------------------------------
    def encode(self, state):
        c = 0
        for (extract, _), idx, stride in zip(self.specs, self.index, self.strides):
            value = extract(state)
            if value not in idx:
                raise ValueError(f"state is not in the expected group for {self.name}")
            c += idx[value] * stride
        return c

    def solve(self, state):
        c = self.encode(state)
        d = int(self.dist[c])
        if d < 0:
            raise ValueError(f"state is unreachable in {self.name}")
        moves = []
        one = np.array([c])
        while d > 0:
            for m, move in enumerate(self.moves):
                n = int(self._next(one, m)[0])
                if self.dist[n] == d - 1:
                    moves.append(move)
                    one[0] = n
                    d -= 1
                    break
        return moves


def make_phases():
    return [
        Phase("phase 1 (edge orientation)", ALL_MOVES,
              [facelet_set(EDGE_REF_FACELETS)]),
        Phase("phase 2 (corner orientation, E slice)", G1_MOVES,
              [facelet_set(CORNER_UD_FACELETS), cubie_set(E_SLICE_IDS)]),
        Phase("phase 3 (M slice, corner tetrads)", G2_MOVES,
              [cubie_perm(CORNER_IDS), cubie_set(M_SLICE_IDS)],
              goal_move_names=G3_MOVES),
        Phase("phase 4 (half turns)", G3_MOVES,
              [cubie_perm(CORNER_IDS), cubie_perm(EDGE_IDS)]),
    ]


CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "thistlethwaite_tables.pickle")
CACHE_VERSION = 1
_PHASES = None


def get_phases(verbose=True):
    """Build (or load from cache) the four phase tables. Done once per process."""
    global _PHASES
    if _PHASES is not None:
        return _PHASES
    phases = make_phases()
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as fh:
            data = pickle.load(fh)
        if data.get("version") == CACHE_VERSION and len(data["phases"]) == len(phases):
            for p, d in zip(phases, data["phases"]):
                p.load(d)
            _PHASES = phases
            return phases
    if verbose:
        print("Building Thistlethwaite tables (first run only)...")
    t0 = time.time()
    for p in phases:
        t1 = time.time()
        p.build()
        if verbose:
            print(f"  {p.name}: {p.N:,} states, max depth {int(p.dist.max())}, {time.time() - t1:.1f}s")
    if verbose:
        print(f"  total {time.time() - t0:.1f}s, cached to {os.path.basename(CACHE)}")
    with open(CACHE, "wb") as fh:
        pickle.dump({"version": CACHE_VERSION, "phases": [p.data() for p in phases]}, fh)
    _PHASES = phases
    return phases


def to_index(state):
    """Accept either an r3 index vector (54,) or coordinate array (54, 4)."""
    state = np.asarray(state)
    if state.shape == coords.shape:
        return np.array([c2i[tuple(c)] for c in state])
    if state.shape != index.shape:
        raise ValueError("state must be an index vector (54,) or coordinates (54, 4)")
    return state.copy()


def solve(state, verbose=False):
    """Return a list of r3 Rotation objects that solve ``state``."""
    state = to_index(state)
    ans = []
    for phase in get_phases(verbose=verbose):
        moves = phase.solve(state)
        for m in moves:
            ans.extend(m.rots)
            for r in m.rots:
                state = r(state)
        if verbose:
            print(f"{phase.name}: {' '.join(map(str, moves)) or '-'}")
    if not np.all(state == index):
        raise RuntimeError("solver produced a wrong answer")
    return ans


if __name__ == "__main__":
    import random
    import sys
    from r3 import RotationSequence

    get_phases()
    n_tests = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    lengths, times = [], []
    for _ in range(n_tests):
        seq = RotationSequence(random.choices(rs, k=random.randrange(20, 100)))
        state = seq(index)
        t0 = time.time()
        ans = solve(state)
        times.append(time.time() - t0)
        assert np.all(RotationSequence(ans)(state) == index)
        lengths.append(len(ans))
    print(f"{n_tests} random scrambles solved and verified")
    print(f"solution length (r3 rotations): mean {np.mean(lengths):.1f}, max {max(lengths)}")
    print(f"solve time: mean {np.mean(times) * 1000:.2f} ms, max {max(times) * 1000:.2f} ms")
