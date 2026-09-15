"""Train the sticker keypoint detector on synthetic renders, CPU only.

    python -m synth.train --steps 6000 --workers 48 --threads 96

Data is rendered on the fly by DataLoader workers (synth.render); every
`--eval-every` steps the model is scored on the real photos in synth/real/
(PCK: fraction of labelled sticker centres with a detected peak within 30 % of
the sticker spacing, and the mean number of peaks). Checkpoints go to
synth/ckpt/ (best real-photo PCK kept as best.pt).
"""
import argparse
import glob
import json
import math
import os
import time

import numpy as np
import torch
from PIL import Image
from torch.utils.data import DataLoader, IterableDataset

from synth.model import KeypointNet, extract_peaks, focal_loss, gaussian_heatmap
from synth.render import render

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIZE, STRIDE, SIGMA = 256, 4, 1.3


class Synthetic(IterableDataset):
    def __init__(self, seed=0):
        self.seed = seed

    def __iter__(self):
        info = torch.utils.data.get_worker_info()
        wid = info.id if info else 0
        rng = np.random.default_rng(self.seed * 1000 + wid + int(time.time()) % 100000)
        while True:
            img, kp, _ = render(rng, SIZE)
            x = torch.from_numpy(np.asarray(img).astype(np.float32) / 255).permute(2, 0, 1)
            hm = gaussian_heatmap(kp / STRIDE, SIZE // STRIDE, SIGMA)
            yield x, hm[None]


def letterbox(img, size=SIZE):
    """Pad an RGB PIL image to a square and resize; returns tensor and (scale, dx, dy)
    mapping model-input coords back to original: orig = (p / scale) - (dx, dy)."""
    w, h = img.size
    side = max(w, h)
    canvas = Image.new("RGB", (side, side), (96, 96, 96))
    dx, dy = (side - w) // 2, (side - h) // 2
    canvas.paste(img, (dx, dy))
    scale = size / side
    x = torch.from_numpy(np.asarray(canvas.resize((size, size), Image.BILINEAR)).astype(np.float32) / 255).permute(2, 0, 1)
    return x, (scale, dx, dy)


def load_real():
    items = []
    for jp in sorted(glob.glob(os.path.join(ROOT, "synth", "real", "*.json"))):
        meta = json.load(open(jp))
        ip = jp[:-5] + ".jpg"
        if not os.path.exists(ip):
            continue
        img = Image.open(ip).convert("RGB")
        pts = np.array(meta["points"]) * np.array(img.size)
        items.append((os.path.basename(ip), img, pts))
    return items


@torch.no_grad()
def evaluate(model, real):
    model.eval()
    pcks, counts = [], []
    for name, img, pts in real:
        x, (scale, dx, dy) = letterbox(img)
        logits = model(x[None])
        peaks = extract_peaks(logits[0, 0])
        pred = np.array([[(px * STRIDE) / scale - dx, (py * STRIDE) / scale - dy] for px, py, _ in peaks]) if peaks else np.zeros((0, 2))
        d = np.sort(np.linalg.norm(pts[:, None] - pts[None], axis=2), 1)[:, 1]
        spacing = float(np.median(d))
        hit = 0
        for p in pts:
            if len(pred) and np.linalg.norm(pred - p, axis=1).min() < 0.3 * spacing:
                hit += 1
        pcks.append(hit / len(pts))
        counts.append(len(pred))
    model.train()
    return float(np.mean(pcks)), float(np.mean(counts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=6000)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--workers", type=int, default=48)
    ap.add_argument("--threads", type=int, default=96)
    ap.add_argument("--eval-every", type=int, default=250)
    ap.add_argument("--width", type=int, default=32)
    ap.add_argument("--resume", default=None)
    args = ap.parse_args()
    torch.set_num_threads(args.threads)
    ckdir = os.path.join(ROOT, "synth", "ckpt")
    os.makedirs(ckdir, exist_ok=True)
    model = KeypointNet(args.width)
    if args.resume:
        model.load_state_dict(torch.load(args.resume, map_location="cpu")["model"])
    print(f"params: {sum(p.numel() for p in model.parameters()) / 1e6:.2f} M", flush=True)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=args.lr, total_steps=args.steps, pct_start=0.1)
    loader = DataLoader(Synthetic(), batch_size=args.batch, num_workers=args.workers, prefetch_factor=4, persistent_workers=True)
    real = load_real()
    print(f"real validation photos: {len(real)}", flush=True)
    best, t0, run = -1.0, time.time(), 0.0
    for step, (x, hm) in enumerate(loader, 1):
        logits = model(x)
        loss = focal_loss(logits, hm)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        opt.step()
        sched.step()
        run = 0.98 * run + 0.02 * loss.item() if step > 1 else loss.item()
        if step % 25 == 0:
            print(f"step {step:5d}  loss {run:.4f}  lr {sched.get_last_lr()[0]:.2e}  {(time.time() - t0) / step:.2f}s/step", flush=True)
        if step % args.eval_every == 0 or step == args.steps:
            pck, cnt = evaluate(model, real) if real else (-1, -1)
            print(f"eval step {step}: real PCK {pck:.3f}  mean peaks {cnt:.1f}", flush=True)
            torch.save({"model": model.state_dict(), "step": step, "pck": pck, "width": args.width},
                       os.path.join(ckdir, "last.pt"))
            if pck >= best:
                best = pck
                torch.save({"model": model.state_dict(), "step": step, "pck": pck, "width": args.width},
                           os.path.join(ckdir, "best.pt"))
        if step >= args.steps:
            break
    print(f"done, best real PCK {best:.3f}", flush=True)


if __name__ == "__main__":
    main()
