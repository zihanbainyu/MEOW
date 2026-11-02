from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image

D = Path(__file__).resolve().parent
P = D / "stim_favila"
plt.rcParams.update({"font.family": "Helvetica", "font.weight": "light", "svg.fonttype": "none"})
INK, GRAY = "#26262f", "#9a9aa2"

def im(n):
    img = Image.open(P / n).convert("RGB")
    s = min(img.size)
    return img.crop(((img.width - s) // 2, (img.height - s) // 2, (img.width + s) // 2, (img.height + s) // 2))

A, C, B = im("barn1.jpg"), im("barn2.jpg"), im("face09.jpg")
nA, nC = im("boardwalk1.jpg"), im("boardwalk2.jpg")
mA, mC, Dm, Em = im("busstop1.jpg"), im("busstop2.jpg"), im("face01.jpg"), im("face11.jpg")

W, H = 7.95, 4.2
sz = 0.62
fig = plt.figure(figsize=(W, H))

def box(x, y, w, h, img=None, ec=None, alpha=1.0, ls="-"):
    ax = fig.add_axes([x / W, y / H, w / W, h / H])
    if img is not None:
        ax.imshow(img, alpha=alpha)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if ec:
        ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, fill=False, ec=ec, lw=1.2, ls=ls))
    return ax

def lab(x, y, t, **k):
    fig.text(x / W, y / H, t, ha="center", va="center", color=k.pop("color", INK), fontsize=k.pop("fs", 8), **k)

def link(xa, xb, y):
    fig.add_artist(plt.Line2D([xa / W, xb / W], [y / H, y / H], color=INK, lw=0.8))

fig.text(1.05 / W, (H - 0.22) / H, "Encoding", ha="left", va="center", color=INK, fontsize=11)
fig.text(4.45 / W, (H - 0.22) / H, "Test", ha="left", va="center", color=INK, fontsize=11)

rows = [("Inferred", "#e0913a", (A, C), (B, "B", B, "B")), ("Matched", "#5b8fc9", (mA, mC), (Dm, "D", Em, "E")), ("Novel", "#6aa84f", (nA, nC), None)]
for r, (name, col, (a, c), enc) in enumerate(rows):
    y = H - 0.55 - sz - r * 1.25
    fig.text(0.12 / W, (y + sz / 2) / H, name, ha="left", va="center", color=col, fontsize=10)
    if enc:
        for x, img, t, e in ((1.05, a, "A", col), (1.75, enc[0], enc[1], GRAY), (2.70, enc[2], enc[3], GRAY), (3.40, c, "C", col)):
            box(x, y, sz, sz, img, e)
            lab(x + sz / 2, y - 0.13, t)
        link(1.05 + sz, 1.75, y + sz / 2)
        link(2.70 + sz, 3.40, y + sz / 2)
    box(4.45, y, sz, sz, a, col)
    lab(4.45 + sz / 2, y - 0.13, "cue A")
    dx, dw = 5.30, 1.55
    box(dx, y, dw, sz, None, GRAY, ls=(0, (3, 2)))
    lab(dx + dw / 2, y - 0.13, "delay")
    ghosts = (a, enc[0]) if enc else (a,)
    gs = 0.40
    gx0 = dx + dw / 2 - (len(ghosts) * gs + (len(ghosts) - 1) * 0.1) / 2
    for k, g in enumerate(ghosts):
        box(gx0 + k * (gs + 0.1), y + (sz - gs) / 2, gs, gs, g, None, alpha=0.35)
    box(7.10, y, sz, sz, c, col)
    lab(7.10 + sz / 2, y - 0.13, "probe C")

for ext in ("png", "svg"):
    fig.savefig(D / f"paradigm_inference.{ext}", dpi=220, facecolor="white")
