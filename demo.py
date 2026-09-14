"""Scramble a cube, solve it with the Thistlethwaite solver, and render both
the scramble and the solution as videos.

    python demo.py                    # 30 random moves, question.mp4 + answer.mp4
    python demo.py -n 50 --seed 7     # reproducible 50-move scramble
    python demo.py --no-video         # just solve and verify
    python demo.py -n 12 --optimal    # shortest possible answer (slow for long scrambles)

The renderer is visual.py: fixed camera, each move animated as a smooth slice turn.
"""
import argparse
import random
import time

import numpy as np

from r3 import coords, index, rs, RotationSequence
import thistlethwaite
import optimal
import visual


def main():
    class Help(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=Help)
    parser.add_argument("-n", "--scramble", type=int, default=30, help="number of random rotations")
    parser.add_argument("--seed", type=int, default=None, help="random seed for a reproducible scramble")
    parser.add_argument("--no-video", action="store_true", help="skip rendering the mp4 files")
    parser.add_argument("--optimal", action="store_true",
                        help="use the optimal solver (guaranteed shortest, exponential time)")
    parser.add_argument("--frames-per-turn", type=int, default=18, help="frames per 90-degree turn")
    parser.add_argument("--fps", type=int, default=24, help="video frame rate")
    args = parser.parse_args()

    random.seed(args.seed)
    if args.optimal:
        optimal.get_tables()                       # build/load tables up front
    else:
        thistlethwaite.get_phases()

    # question
    question = RotationSequence(random.choices(rs, k=args.scramble))
    index_q = question(index)
    coords_q = question(coords)

    # answer
    t0 = time.time()
    answer = RotationSequence(optimal.solve(index_q) if args.optimal else thistlethwaite.solve(index_q))
    dt = time.time() - t0

    solved = np.all(answer(index_q) == index)
    print(f"Question ({len(question)} moves): {question}")
    print(f"Answer   ({len(answer)} moves): {answer}")
    print(f"Solved in {dt * 1000:.1f} ms, verified: {solved}")
    assert solved

    if args.no_video:
        return
    kwargs = dict(frames_per_turn=args.frames_per_turn, fps=args.fps)
    t0 = time.time()
    visual.export_video(coords, question.seq, "question.mp4", **kwargs)
    visual.export_video(coords_q, answer.seq, "answer.mp4", **kwargs)
    print(f"Exported question.mp4 and answer.mp4 in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
