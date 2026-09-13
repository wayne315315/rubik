"""Optimal solver: IDA* with pattern databases, vectorised with numpy.

Guarantees the shortest solution for the chosen move metric:

    "r3"   the 18 r3 rotations, each costing 1 (12 quarter turns + 6 middle
           slice rotations)                                   -- default
    "qtm"  quarter-turn metric: 12 quarter turns only
    "htm"  half-turn metric: 12 quarter turns + 6 half turns, each costing 1

The cube state is a 6-tuple of small coordinates, each with a tiny move table:

    co  corner orientation           2,187      cp  corner permutation   40,320
    eo  edge orientation             2,048      eE/eM/eS  ordered positions of
                                                the 4 E/M/S-slice edges  11,880

Four pattern databases (exact distances in the chosen metric, built once by
BFS and cached) give an admissible, consistent heuristic:

    h = max( PDB[co, cp], PDB[eo, eE], PDB[eo, eM], PDB[eo, eS] )

IDA* then runs as a depth-first search over numpy arrays of nodes, so tens of
millions of nodes per second are examined per core.  Once an iteration gets
large, the search tree is split a few levels below the root into thousands of
independent subtrees that are searched in parallel on every CPU core (fork,
tables shared copy-on-write).  Redundant move sequences (same-axis runs that
are not the cheapest canonical word) are pruned by a small automaton derived
from the move set.

Running time grows exponentially with the optimal length.  On a 192-core
machine a fully random state (optimal 18 moves, ~5e10 nodes) takes about half
a minute; each extra move costs roughly a factor of 10.

Usage::

    from optimal import solve
    ans = solve(index_q)                 # list of r3 Rotation objects
    ans = solve(index_q, metric="htm")   # shortest in half-turn metric
"""
import itertools
import multiprocessing as mp
import os
import pickle
import time

import numpy as np

from r3 import index
from thistlethwaite import (MOVES, facelet_set, cubie_perm, enumerate_coord, to_index,
                            product_strides, product_bfs,
                            CORNER_UD_FACELETS, EDGE_REF_FACELETS, CORNER_IDS,
                            E_SLICE_IDS, M_SLICE_IDS, S_SLICE_IDS)

QUARTER = [n for n in MOVES if n[1] in "pn" and n[2] in "pn"]
HALF = [n for n in MOVES if n[2] == "2"]
SLICE = [n for n in MOVES if n[1] == "0"]
METRICS = {
    "r3": QUARTER + SLICE,
    "qtm": QUARTER,
    "htm": QUARTER + HALF,
}


# ----------------------------------------------------------------------------
# Canonical same-axis runs: all moves about one axis commute, so a run of
# same-axis moves is only allowed if it is the cheapest (then lexicographically
# first) word for the group element it produces.
# ----------------------------------------------------------------------------
def build_automaton(move_names):
    moves = [MOVES[n] for n in move_names]
    axis = [n[0] for n in move_names]
    canonical = {()}
    for a in "xyz":
        ids = [i for i, ax in enumerate(axis) if ax == a]
        best = {}
        for length in range(1, 5):
            for word in itertools.product(ids, repeat=length):
                perm = index
                for i in word:
                    perm = moves[i].fperm[perm]
                key = perm.tobytes()
                if key not in best or (len(word), word) < (len(best[key]), best[key]):
                    best[key] = word
        canonical |= set(best.values())
    words = sorted(canonical, key=lambda w: (len(w), w))
    wid = {w: i for i, w in enumerate(words)}
    nxt = np.full((len(words), len(moves)), -1, dtype=np.int16)
    for w, i in wid.items():
        for m in range(len(moves)):
            if w and axis[w[-1]] == axis[m]:
                w2 = w + (m,)
                if w2 in wid:
                    nxt[i, m] = wid[w2]
            else:
                nxt[i, m] = wid[(m,)]
    return nxt


# ----------------------------------------------------------------------------
# Tables per metric: coordinate move tables, pattern databases, automaton
# ----------------------------------------------------------------------------
CACHE_VERSION = 1
_TABLES = {}


class Tables:
    def __init__(self, metric):
        if metric not in METRICS:
            raise ValueError(f"metric must be one of {list(METRICS)}")
        self.metric = metric
        self.move_names = METRICS[metric]
        self.moves = [MOVES[n] for n in self.move_names]
        self.specs = [facelet_set(CORNER_UD_FACELETS), cubie_perm(CORNER_IDS),
                      facelet_set(EDGE_REF_FACELETS), cubie_perm(E_SLICE_IDS),
                      cubie_perm(M_SLICE_IDS), cubie_perm(S_SLICE_IDS)]
        self.pdb_pairs = [(0, 1), (2, 3), (2, 4), (2, 5)]

    @property
    def cache_file(self):
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), f"optimal_tables_{self.metric}.pickle")

    def build(self, verbose=True):
        t0 = time.time()
        self.index, self.tables = [], []
        for extract, apply in self.specs:
            idx, table = enumerate_coord(extract(index), self.moves, apply)
            self.index.append(idx)
            self.tables.append(table)
        self.sizes = [len(i) for i in self.index]
        self.pdbs = []
        for a, b in self.pdb_pairs:
            t1 = time.time()
            sizes = [self.sizes[a], self.sizes[b]]
            strides, N = product_strides(sizes)
            dist = product_bfs([self.tables[a], self.tables[b]], sizes, strides, N,
                               np.array([0]), range(len(self.moves)))
            self.pdbs.append(dist)
            if verbose:
                print(f"  pattern database {(a, b)}: {N:,} states, max {int(dist.max())}, "
                      f"mean {dist.mean():.2f}, {time.time() - t1:.0f}s")
        self.automaton = build_automaton(self.move_names)
        if verbose:
            print(f"  total {time.time() - t0:.0f}s, cached to {os.path.basename(self.cache_file)}")
        with open(self.cache_file, "wb") as fh:
            pickle.dump({"version": CACHE_VERSION, "index": self.index, "tables": self.tables,
                         "pdbs": self.pdbs, "automaton": self.automaton}, fh, protocol=4)
        return self

    def load(self):
        with open(self.cache_file, "rb") as fh:
            data = pickle.load(fh)
        if data.get("version") != CACHE_VERSION:
            raise ValueError("stale cache")
        self.index, self.tables = data["index"], data["tables"]
        self.pdbs, self.automaton = data["pdbs"], data["automaton"]
        self.sizes = [len(i) for i in self.index]
        return self

    def encode(self, state):
        coords = []
        for (extract, _), idx in zip(self.specs, self.index):
            value = extract(state)
            if value not in idx:
                raise ValueError("state is not a legal cube state")
            coords.append(idx[value])
        return np.array(coords, dtype=np.int64)

    def heuristic(self, coords):
        """Admissible lower bound for a (6,) coordinate vector."""
        return max(int(pdb[coords[a] + coords[b] * self.sizes[a]]) for pdb, (a, b) in zip(self.pdbs, self.pdb_pairs))


def get_tables(metric="r3", verbose=True):
    if metric in _TABLES:
        return _TABLES[metric]
    t = Tables(metric)
    try:
        t.load()
    except (FileNotFoundError, ValueError, KeyError):
        if verbose:
            print(f"Building optimal-solver tables for metric '{metric}' (first run only, ~2 min)...")
        t.build(verbose=verbose)
    _TABLES[metric] = t
    return t


# ----------------------------------------------------------------------------
# IDA*: depth-first search over arrays of nodes, one chunk of parents at a time
# ----------------------------------------------------------------------------
class _Frame:
    __slots__ = ("g", "coords", "run", "parent", "move", "pos")

    def __init__(self, g, coords, run, parent, move):
        self.g, self.coords, self.run, self.parent, self.move, self.pos = g, coords, run, parent, move, 0


def _search(T, root, bound, chunk, run0=0):
    """One IDA* iteration. Returns (list of move ids or None, nodes generated)."""
    tables, sizes, pdbs, pairs, nxt = T.tables, T.sizes, T.pdbs, T.pdb_pairs, T.automaton
    nm = len(T.moves)
    nodes = 0
    stack = [_Frame(0, root.reshape(6, 1), np.array([run0], dtype=np.int16), None, None)]
    while stack:
        fr = stack[-1]
        if fr.pos >= fr.coords.shape[1]:
            stack.pop()
            continue
        base, stop = fr.pos, min(fr.pos + chunk, fr.coords.shape[1])
        fr.pos = stop
        g = fr.g + 1
        slack = bound - g                         # h must be <= slack
        P, R = fr.coords[:, base:stop], fr.run[base:stop]
        kids, runs, parents, mids = [], [], [], []
        for m in range(nm):
            nr = nxt[R, m]
            idx = np.flatnonzero(nr >= 0)
            if not len(idx):
                continue
            nodes += len(idx)
            c = [tables[k][P[k, idx], m] for k in range(6)]
            for pdb, (a, b) in zip(pdbs, pairs):       # prune coordinate by coordinate
                keep = pdb[c[a] + c[b] * sizes[a]] <= slack
                if not keep.all():
                    idx = idx[keep]
                    c = [ck[keep] for ck in c]
                if not len(idx):
                    break
            if not len(idx):
                continue
            c = np.stack(c)
            solved = np.flatnonzero(~c.any(axis=0))
            if len(solved):                             # goal: walk the frames back up
                path = [m]
                i = base + idx[solved[0]]
                for f in reversed(stack):
                    if f.move is None:
                        break
                    path.append(int(f.move[i]))
                    i = int(f.parent[i])
                return path[::-1], nodes
            if g < bound:
                kids.append(c)
                runs.append(nr[idx])
                parents.append(base + idx)
                mids.append(np.full(len(idx), m, dtype=np.int8))
        if kids:
            stack.append(_Frame(g, np.concatenate(kids, axis=1), np.concatenate(runs),
                                np.concatenate(parents), np.concatenate(mids)))
    return None, nodes


# ----------------------------------------------------------------------------
# Parallel iteration: split the tree a few levels below the root, search the
# subtrees on all cores. Any solution found at this bound is optimal because
# every smaller bound was already exhausted.
# ----------------------------------------------------------------------------
def _expand(T, root, bound, depth):
    """All (path, coords, run) at the given depth whose subtree can still
    contain a solution of length bound. Returns (solution path or None, list)."""
    tables, sizes, pdbs, pairs, nxt = T.tables, T.sizes, T.pdbs, T.pdb_pairs, T.automaton
    level = [((), root, 0)]
    for g in range(1, depth + 1):
        new = []
        for path, c, run in level:
            for m in range(len(T.moves)):
                nr = nxt[run, m]
                if nr < 0:
                    continue
                c2 = np.array([tables[k][c[k], m] for k in range(6)], dtype=np.int64)
                h = max(int(pdb[c2[a] + c2[b] * sizes[a]]) for pdb, (a, b) in zip(pdbs, pairs))
                if h == 0:
                    return list(path) + [m], []
                if g + h <= bound:
                    new.append((path + (m,), c2, int(nr)))
        level = new
    return None, level


def _worker(args):
    metric, path, coords, run, bound, chunk = args
    path2, nodes = _search(_TABLES[metric], coords, bound, chunk, run0=run)
    return (None if path2 is None else list(path) + path2), nodes


def _search_parallel(T, root, bound, chunk, workers):
    depth = 1
    while True:
        path, tasks = _expand(T, root, bound, depth)
        if path is not None:
            return path, 0
        if not tasks:
            return None, 0
        if len(tasks) >= 8 * workers or depth == 4:
            break
        depth += 1
    args = [(T.metric, p, c, r, bound - depth, chunk) for p, c, r in tasks]
    nodes = 0
    ctx = mp.get_context("fork")                       # tables inherited, copy-on-write
    with ctx.Pool(min(workers, len(args))) as pool:
        for path, n in pool.imap_unordered(_worker, args, chunksize=1):
            nodes += n
            if path is not None:
                pool.terminate()
                return path, nodes
    return None, nodes


def solve(state, metric="r3", max_depth=40, chunk=16384, verbose=True, workers=None,
          parallel_threshold=1_000_000):
    """Return an optimal solution (list of r3 Rotation objects) for ``state``.

    An iteration runs on all cores (``workers``, default every CPU) once the
    previous iteration generated more than ``parallel_threshold`` nodes.
    """
    T = get_tables(metric, verbose=verbose)
    root = T.encode(to_index(state))
    workers = workers or os.cpu_count()
    total = 0
    last = 0
    t0 = time.time()
    for bound in range(T.heuristic(root), max_depth + 1):
        t1 = time.time()
        if workers > 1 and last > parallel_threshold:
            path, nodes = _search_parallel(T, root, bound, chunk, workers)
            how = f"{workers} cores"
        else:
            path, nodes = _search(T, root, bound, chunk)
            how = "1 core"
        total += nodes
        last = nodes
        if verbose:
            print(f"  depth {bound:2d}: {nodes:>14,} nodes, {time.time() - t1:8.1f}s ({how})"
                  + ("  found" if path is not None else ""))
        if path is not None:
            if verbose:
                print(f"optimal length {len(path)} ({metric}), {total:,} nodes, {time.time() - t0:.1f}s")
            ans = []
            for m in path:
                ans.extend(T.moves[m].rots)
            return ans
    raise RuntimeError(f"no solution within max_depth={max_depth}")


if __name__ == "__main__":
    import random
    import sys
    from r3 import rs, RotationSequence
    from thistlethwaite import solve as solve_fast

    metric = sys.argv[1] if len(sys.argv) > 1 else "r3"
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else None
    random.seed(seed)
    T = get_tables(metric)
    question = RotationSequence(random.choices([MOVES[n].rots[0] for n in QUARTER + SLICE], k=k))
    state = question(index)
    print(f"Question ({len(question)}): {question}")
    ans = RotationSequence(solve(state, metric=metric))
    assert np.all(ans(state) == index)
    print(f"Optimal ({len(ans)} rotations): {ans}")
    fast = RotationSequence(solve_fast(state))
    print(f"Thistlethwaite for comparison: {len(fast)} rotations")
