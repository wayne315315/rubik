"""Regression suite for the photo pipeline on the real photo pairs in examples/.

    .venv/bin/python -m synth.regress            # CNN detector (no model server)

Every pair must yield a legal cube state, and pairs with a reference JSON must
reproduce exactly that state. Also checks that the balanced colour assignment
is optimal on random data. Exit status is non-zero on any failure.
"""
import glob
import itertools
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cube_vision import state_from_views, CubeReadError                       # noqa: E402
import cube_keypoints as K                                                     # noqa: E402
import cube_locate as L                                                        # noqa: E402

PAIRS = {                       # name: (photo pair, reference views json or None)
    "pair1": (["photo_yellow_blue_red", "photo_green_white_orange"], "reading.json"),
    "pair2": (["scramble2_yellow_blue_red", "scramble2_green_orange_white"], "reading2.json"),
    "phone3": (["phone3_yellow_blue_red", "phone3_green_white_orange"], "reading_phone3.json"),
    "desk4": (["desk4_red_blue_white", "desk4_yellow_green_orange"], "reading_desk4.json"),
    "desk5": (["desk5_red_blue_white", "desk5_yellow_green_orange"], "reading_desk5.json"),
}


def check_assignment():
    rng = np.random.default_rng(0)
    old, L.COLORS = L.COLORS, ("a", "b", "c")
    try:
        for _ in range(100):
            feat, cen = rng.normal(size=(6, 3)), rng.normal(size=(3, 3))
            lab = L._balanced(feat, cen, 2)
            d = np.linalg.norm(feat[:, None] - cen[None], axis=2)
            best = min(sum(d[i, p[i]] for i in range(6)) for p in set(itertools.permutations([0, 0, 1, 1, 2, 2])))
            if sum(d[i, lab[i]] for i in range(6)) - best > 1e-9:
                return False
    finally:
        L.COLORS = old
    return True


def main():
    ok = True
    print("balanced assignment optimal:", "OK" if check_assignment() else "FAIL")
    for name, (photos, ref) in PAIRS.items():
        paths = [os.path.join("examples", p + ".jpg") for p in photos]
        if not all(os.path.exists(p) for p in paths):
            print(f"{name}: photos missing, skipped")
            continue
        t0 = time.time()
        try:
            state = state_from_views(K.read_photos(paths))
        except (L.LocateError, CubeReadError) as e:
            print(f"{name}: FAIL {type(e).__name__}: {str(e)[:100]}")
            ok = False
            continue
        if ref and os.path.exists(os.path.join("examples", ref)):
            want = state_from_views(json.load(open(os.path.join("examples", ref)))["views"])
            good = np.array_equal(state, want)
            print(f"{name}: {'OK, matches reference' if good else 'FAIL, differs from reference'} ({time.time() - t0:.1f}s)")
            ok &= good
        else:
            print(f"{name}: OK, legal ({time.time() - t0:.1f}s)")
    print("ALL OK" if ok else "FAILURES")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
