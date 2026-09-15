"""Small fully-convolutional sticker-centre heatmap network (CPU friendly).

Input: RGB 256x256 in [0,1]. Output: one heatmap at stride 4 (64x64) whose
peaks are sticker centres. About 1.2 M parameters; a forward pass on a few
CPU cores takes tens of milliseconds.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


def conv_bn(cin, cout, k=3, s=1, d=1):
    p = d * (k // 2)
    return nn.Sequential(nn.Conv2d(cin, cout, k, s, p, dilation=d, bias=False),
                         nn.BatchNorm2d(cout), nn.ReLU(inplace=True))


class Residual(nn.Module):
    def __init__(self, c, d=1):
        super().__init__()
        self.a = conv_bn(c, c, d=d)
        self.b = nn.Sequential(nn.Conv2d(c, c, 3, 1, d, dilation=d, bias=False), nn.BatchNorm2d(c))

    def forward(self, x):
        return F.relu(x + self.b(self.a(x)))


class KeypointNet(nn.Module):
    def __init__(self, width=32):
        super().__init__()
        w = width
        self.stem = nn.Sequential(conv_bn(3, w, s=2), conv_bn(w, w))                 # /2
        self.down = nn.Sequential(conv_bn(w, 2 * w, s=2), conv_bn(2 * w, 2 * w))     # /4
        self.body = nn.Sequential(Residual(2 * w), Residual(2 * w, d=2), Residual(2 * w, d=4),
                                  Residual(2 * w, d=8), Residual(2 * w, d=2), Residual(2 * w))
        self.head = nn.Sequential(conv_bn(2 * w, w), nn.Conv2d(w, 1, 1))
        nn.init.constant_(self.head[-1].bias, -3.0)                                  # sparse peaks prior

    def forward(self, x):
        return self.head(self.body(self.down(self.stem(x))))


def focal_loss(logits, target, alpha=2, beta=4):
    """CenterNet penalty-reduced focal loss; target is a Gaussian heatmap in [0,1]."""
    p = torch.sigmoid(logits).clamp(1e-4, 1 - 1e-4)
    pos = target.eq(1).float()
    neg = 1 - pos
    pos_loss = torch.log(p) * (1 - p) ** alpha * pos
    neg_loss = torch.log(1 - p) * p ** alpha * (1 - target) ** beta * neg
    n = pos.sum().clamp(min=1)
    return -(pos_loss.sum() + neg_loss.sum()) / n


def gaussian_heatmap(points, size, sigma):
    """(n,2) points in output-pixel coords -> (size,size) heatmap with peak 1 at each."""
    hm = torch.zeros(size, size)
    ys, xs = torch.meshgrid(torch.arange(size, dtype=torch.float32), torch.arange(size, dtype=torch.float32), indexing="ij")
    for x, y in points:
        if not (0 <= x < size and 0 <= y < size):
            continue
        g = torch.exp(-((xs - x) ** 2 + (ys - y) ** 2) / (2 * sigma ** 2))
        hm = torch.maximum(hm, g)
    # exact 1 at the nearest cell so the focal loss sees a positive
    for x, y in points:
        xi, yi = int(round(float(x))), int(round(float(y)))
        if 0 <= xi < size and 0 <= yi < size:
            hm[yi, xi] = 1.0
    return hm


def extract_peaks(logits, threshold=0.3, max_peaks=40):
    """Peaks of a (1,1,H,W) or (H,W) logit map: list of (x, y, score) with
    sub-pixel refinement, in heatmap pixel coordinates."""
    hm = torch.sigmoid(logits.detach()).reshape(1, 1, *logits.shape[-2:])
    pooled = F.max_pool2d(hm, 3, 1, 1)
    keep = (hm == pooled) & (hm > threshold)
    ys, xs = torch.nonzero(keep[0, 0], as_tuple=True)
    scores = hm[0, 0, ys, xs]
    order = torch.argsort(scores, descending=True)[:max_peaks]
    H, W = hm.shape[-2:]
    out = []
    h = hm[0, 0]
    for k in order:
        x, y = int(xs[k]), int(ys[k])
        dx = dy = 0.0
        if 0 < x < W - 1:
            l, c, r = float(h[y, x - 1]), float(h[y, x]), float(h[y, x + 1])
            den = l - 2 * c + r
            dx = 0.5 * (l - r) / den if den < 0 else 0.0
        if 0 < y < H - 1:
            u, c, d = float(h[y - 1, x]), float(h[y, x]), float(h[y + 1, x])
            den = u - 2 * c + d
            dy = 0.5 * (u - d) / den if den < 0 else 0.0
        out.append((x + max(-0.5, min(0.5, dx)), y + max(-0.5, min(0.5, dy)), float(scores[k])))
    return out
