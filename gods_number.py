"""Bounds on God's number for the r3 move set (18 rotations, each costing 1).

Computing God's number exactly is far beyond a laptop: the half-turn-metric
result (20) took about 35 CPU-years with heavy symmetry and coset methods.
What this script does establish rigorously:

  lower bound  (a) counting: if there are fewer distinct canonical sequences
                   of length <= d-1 than cube states, some state needs >= d;
               (b) any position we solve optimally with optimal.py whose
                   solution has length L proves God's number >= L.
  upper bound  every quarter-turn sequence is an r3 sequence, and God's number
               in the quarter-turn metric is 26 (Rokicki & Davidson, 2014),
               so the r3 God's number is <= 26.

    python gods_number.py                 # counting bound only (instant)
    python gods_number.py --random 20     # + optimal lengths of 20 random states
    python gods_number.py --hard          # + four-spot, superflip, superflip+four-spot
                                          #   (superflip ones can take a long time)
"""
import argparse
import random
import time

import numpy as np

from r3 import index, rs, b2i, c2i, coords, RotationSequence
from optimal import METRICS, build_automaton, solve, get_tables
from thistlethwaite import MOVES

N_STATES = 43_252_003_274_489_856_000            # 8!*3^7*12!*2^10

# Standard face notation -> r3 rotation (clockwise as seen from that face).
FACE = {"U": "zpn", "D": "znp", "F": "xpn", "B": "xnp", "R": "ypn", "L": "ynp"}
ROT = {str(r): r for r in rs}


def from_notation(s):
    """'F2 B2 U D'' -> list of r3 rotations."""
    seq = []
    for tok in s.split():
        name = FACE[tok[0]]
        if tok.endswith("'"):
            name = name[:2] + {"p": "n", "n": "p"}[name[2]]
        seq.extend([ROT[name]] * (2 if tok.endswith("2") else 1))
    return seq


def superflip():
    """Every edge flipped in place: swap the two stickers of each edge."""
    state = index.copy()
    for block, stickers in b2i.items():
        if sum(v == 0 for v in block) == 1:
            i, j = stickers
            state[i], state[j] = j, i
    return state


def four_spot():
    """Four dots: opposite side centres swapped. With r3's fixed centres this is
    every corner and edge rotated 180 degrees about z while centres stay."""
    state = index.copy()
    for i, (x, y, z, n) in enumerate(coords):
        if sum(v == 0 for v in (x, y, z)) < 2:
            state[i] = c2i[(-x, -y, z, n)]
    return state


# ----------------------------------------------------------------------------
# (a) counting lower bound: number of canonical sequences per length
# ----------------------------------------------------------------------------
def counting_bound(metric="r3"):
    nxt = build_automaton(METRICS[metric])
    n = len(nxt)
    # transfer matrix with exact integers: A[s][t] = #moves leading s -> t
    A = [[0] * n for _ in range(n)]
    for s in range(n):
        for t in nxt[s]:
            if t >= 0:
                A[s][t] += 1
    v = [0] * n
    v[0] = 1                                  # start: empty run
    total, d = 1, 0
    per_length = [1]
    while total < N_STATES:
        v = [sum(v[s] * A[s][t] for s in range(n)) for t in range(n)]
        d += 1
        per_length.append(sum(v))
        total += per_length[-1]
    return d, per_length, total


def solve_and_report(name, state, metric):
    t0 = time.time()
    print(f"  {name}:")
    ans = solve(state, metric=metric, verbose=True)
    assert np.all(RotationSequence(ans)(state) == index)
    print(f"  {name}: optimal length {len(ans)}  ({time.time() - t0:.0f}s)  {RotationSequence(ans)}", flush=True)
    return len(ans)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metric", default="r3", choices=list(METRICS))
    parser.add_argument("--random", type=int, default=0, help="optimally solve this many random states")
    parser.add_argument("--hard", action="store_true", help="optimally solve superflip and superflip+four-spot")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    print(f"Move set '{args.metric}': {len(METRICS[args.metric])} moves, {N_STATES:,} states")

    d, per_length, total = counting_bound(args.metric)
    print("\nCounting argument (canonical sequences per length):")
    for i, c in enumerate(per_length):
        print(f"  length {i:2d}: {c:>26,}")
    print(f"  sequences of length <= {d - 1}: {total - per_length[-1]:,} < {N_STATES:,} states")
    print(f"  => God's number >= {d}")
    lower = d

    if args.random or args.hard:
        get_tables(args.metric)
    if args.random:
        print(f"\nOptimal lengths of {args.random} random states (seed {args.seed}):")
        random.seed(args.seed)
        lengths = []
        for i in range(args.random):
            state = RotationSequence(random.choices(rs, k=200))(index)
            lengths.append(solve_and_report(f"random {i}", state, args.metric))
        print(f"  mean {np.mean(lengths):.2f}, max {max(lengths)}")
        lower = max(lower, max(lengths))

    if args.hard:
        print("\nKnown hard positions:")
        sf, fs = superflip(), four_spot()
        hard = {
            "four-spot": fs,
            "superflip": sf,
            "superflip + four-spot": sf[fs],         # superflip is central, order is irrelevant
        }
        for name, state in hard.items():
            lower = max(lower, solve_and_report(name, state, args.metric))

    upper = {"r3": 26, "qtm": 26, "htm": 20}[args.metric]
    print(f"\nGod's number for '{args.metric}': {lower} <= N <= {upper}")
    if args.metric == "r3":
        print("  upper bound: quarter-turn metric God's number is 26 and every quarter-turn\n"
              "  sequence is an r3 sequence. The exact value is out of reach on one machine.")


if __name__ == "__main__":
    main()
