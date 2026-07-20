#!/usr/bin/env python3
"""
Predicted 16x16 representational dissimilarity matrix (RDM) for the MST-Back
fMRI RSA. Rendered in the lab heatmap style (Purples, pixel cells, bottom
gradient bar).

Design space (16 conditions):
    item  in {A1, B1, A2, B2}      A/B = pair exemplar, 1/2 = encoding/retrieval
    cond  in {compared, isolated}
    acc   in {correct, incorrect}

Mechanisms live in B1's row of each 4x4 block:
    A1 . B1  -> pattern separation   (2a)
    A2 . B1  -> predictive recall    (2c)
    B2 . B1  -> pattern completion   (2b)

The matrix is generated from latent representation vectors so it is a valid,
symmetric RDM; parameters are tuned to the predicted directions, not fit data.
"""
import os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap

mpl.rcParams.update({
    "svg.fonttype": "none",
    "font.family": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 11,
    "text.color": "#26262F",
})

OUT = "/Users/bai/Documents/GitHub/MEOW/Experiment_2/expected_results"
os.makedirs(OUT, exist_ok=True)

INK   = "#26262F"
MUTED = "#6B6B78"
COMP  = "#8E7CC3"   # compared / purple
ISO   = "#5B8FC9"   # isolated / blue
# mechanism highlight colors
C_SEP = "#C7603F"   # separation
C_COM = "#1F8A8A"   # completion
C_REC = "#C9A227"   # predictive recall

# heatmap colormap matched to the lab gaze heatmaps (white -> deep purple)
PURPLES = LinearSegmentedColormap.from_list(
    "labpurple", ["#FFFFFF", "#EDE9F5", "#CFC3E6", "#A996D2", "#7B62B3", "#4B3B86"])


# --------------------------------------------------------------------------
# build predicted representations -> RDM
# --------------------------------------------------------------------------
def build_rdm():
    rng = np.random.default_rng(0)
    D = 120

    s  = rng.normal(size=D)                 # shared component (pair are similar)
    uA = rng.normal(size=D)
    uB = rng.normal(size=D)
    baseA = 1.0 * s + 0.95 * uA
    baseB = 1.0 * s + 0.95 * uB             # baseline corr(A,B) ~ 0.5
    dirAB = baseA - baseB

    def enc_pattern(item, cond):
        """stored encoding representation for exemplar under a condition."""
        base = baseA if item == "A" else baseB
        sep = 0.95 if cond == "comp" else 0.12          # comparison differentiates
        sign = +1 if item == "A" else -1
        return base + sign * sep * 0.5 * dirAB

    def rep(item, phase, cond, acc):
        e = enc_pattern(item, cond)
        if phase == 1:                                   # encoding
            return e + 0.10 * rng.normal(size=D)
        # phase 2: retrieval  -> reinstatement of encoding pattern (acc-scaled)
        reinstate = 0.70 if acc == "cor" else 0.30
        r = reinstate * e + (1 - reinstate) * 1.7 * rng.normal(size=D)
        if item == "A" and cond == "comp":               # predictive recall of B
            recall = 0.65 if acc == "cor" else 0.15
            r = r + recall * enc_pattern("B", cond)
        return r + 0.10 * rng.normal(size=D)

    # ordering: acc (outer) -> cond -> item
    items = [("A1", "A", 1), ("B1", "B", 1), ("A2", "A", 2), ("B2", "B", 2)]
    order, labels = [], []
    for acc in ["cor", "inc"]:
        for cond in ["comp", "iso"]:
            for name, it, ph in items:
                order.append(rep(it, ph, cond, acc))
                labels.append(name)
    X = np.vstack(order)
    R = np.corrcoef(X)
    RDM = 1.0 - R
    np.fill_diagonal(RDM, 0.0)
    return RDM, labels


# --------------------------------------------------------------------------
# plot in lab heatmap style
# --------------------------------------------------------------------------
def plot_rdm(RDM, labels):
    n = 16
    fig = plt.figure(figsize=(11.4, 8.7))
    ax = fig.add_axes([0.085, 0.17, 0.537, 0.70])   # square heatmap

    vmax = 1.15
    im = ax.imshow(RDM, cmap=PURPLES, vmin=0, vmax=vmax, interpolation="nearest")

    # pixel gridlines (fine, white)
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", length=0)

    # block dividers: every 4 (cond) thin, at 8 (accuracy) thick
    for b in [4, 8, 12]:
        lw = 2.6 if b == 8 else 1.4
        ax.axhline(b - 0.5, color=INK, lw=lw)
        ax.axvline(b - 0.5, color=INK, lw=lw)
    for sp in ax.spines.values():
        sp.set_edgecolor(INK); sp.set_linewidth(1.6)

    # item tick labels
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(labels, fontsize=8.5, rotation=90)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.tick_params(length=0)

    # ---- group brackets (condition = inner, accuracy = outer) ----
    def bracket(axis, start, stop, offset, text, color, pad, fs, lw=1.6):
        """draw a group bracket outside the axis in data coords."""
        mid = (start + stop) / 2
        if axis == "x":
            y = n - 0.5 + offset
            ax.plot([start - 0.5, stop - 0.5], [y, y], color=color, lw=lw,
                    clip_on=False)
            ax.text(mid - 0.5, y + pad, text, ha="center", va="top",
                    color=color, fontsize=fs, fontweight="bold", clip_on=False)
        else:
            x = -0.5 - offset
            ax.plot([x, x], [start - 0.5, stop - 0.5], color=color, lw=lw,
                    clip_on=False)
            ax.text(x - pad, mid - 0.5, text, ha="right", va="center",
                    color=color, fontsize=fs, fontweight="bold",
                    rotation=90, clip_on=False)

    # condition brackets (x)
    for x0, cond, col in [(0, "Compared", COMP), (4, "Isolated", ISO),
                          (8, "Compared", COMP), (12, "Isolated", ISO)]:
        bracket("x", x0, x0 + 4, 0.9, cond, col, 0.55, 9)
    # accuracy brackets (x, outer)
    bracket("x", 0, 8, 2.5, "Correct", INK, 0.6, 11, lw=2.2)
    bracket("x", 8, 16, 2.5, "Incorrect", MUTED, 0.6, 11, lw=2.2)
    # condition brackets (y)
    for y0, cond, col in [(0, "Compared", COMP), (4, "Isolated", ISO),
                          (8, "Compared", COMP), (12, "Isolated", ISO)]:
        bracket("y", y0, y0 + 4, 0.9, cond, col, 0.55, 9)
    bracket("y", 0, 8, 2.5, "Correct", INK, 0.6, 11, lw=2.2)
    bracket("y", 8, 16, 2.5, "Incorrect", MUTED, 0.6, 11, lw=2.2)

    # ---- highlight the three mechanism cells in Correct-Compared block ----
    # block indices: A1=0, B1=1, A2=2, B2=3  (B1 is the reference row/col)
    mech = [((0, 1), C_SEP, "Pattern separation", "A1 · B1"),
            ((2, 1), C_REC, "Predictive recall",  "A2 · B1"),
            ((3, 1), C_COM, "Pattern completion", "B2 · B1")]
    for (c, r), col, _, _ in mech:
        ax.add_patch(Rectangle((c - 0.5, r - 0.5), 1, 1, fill=False,
                               edgecolor=col, lw=2.8, zorder=6))
        ax.add_patch(Rectangle((r - 0.5, c - 0.5), 1, 1, fill=False,
                               edgecolor=col, lw=2.8, zorder=6))

    # ---- mechanism key panel (own axes on the right) ----
    kax = fig.add_axes([0.66, 0.17, 0.325, 0.70]); kax.axis("off")
    kax.set_xlim(0, 1); kax.set_ylim(0, 1)
    kax.text(0.0, 0.97, "Mechanisms", fontsize=13, fontweight="bold", color=INK)
    kax.text(0.0, 0.925, "read from B1's row within each 4×4 block",
             fontsize=9.5, color=MUTED)
    for i, ((_, _), col, name, pair) in enumerate(mech):
        yy = 0.85 - i * 0.10
        kax.add_patch(Rectangle((0.0, yy - 0.028), 0.05, 0.05, fill=False,
                                edgecolor=col, lw=2.8))
        kax.text(0.09, yy + 0.006, name, fontsize=11, color=INK, va="center")
        kax.text(0.09, yy - 0.035, pair, fontsize=9.5, color=col, va="center",
                 fontweight="bold")

    kax.text(0.0, 0.50, "Predicted  ·  compared, correct", fontsize=10.5,
             fontweight="bold", color=INK)
    reads = [(C_SEP, "separation", "high dissimilarity  (dark)"),
             (C_COM, "completion", "low dissimilarity  (light)"),
             (C_REC, "predictive recall", "low dissimilarity  (light)")]
    for i, (col, a, b) in enumerate(reads):
        yy = 0.44 - i * 0.055
        kax.text(0.0, yy, f"{a}", fontsize=10, color=col, va="center",
                 fontweight="bold")
        kax.text(0.42, yy, f"→  {b}", fontsize=10, color=MUTED, va="center")
    kax.text(0.0, 0.24,
             "Accuracy contrast: reinstatement (completion\n"
             "& recall) is predicted to weaken on incorrect\n"
             "trials — those cells darken relative to correct.",
             fontsize=9.2, color=MUTED, va="top", linespacing=1.5)

    # ---- bottom gradient colorbar in lab style ----
    cax = fig.add_axes([0.19, 0.075, 0.33, 0.020])
    grad = np.linspace(0, 1, 256)[None, :]
    cax.imshow(grad, aspect="auto", cmap=PURPLES)
    cax.set_xticks([]); cax.set_yticks([])
    for sp in cax.spines.values():
        sp.set_edgecolor(INK); sp.set_linewidth(1.0)
    cax.text(-0.02, 0.5, "low", transform=cax.transAxes, ha="right",
             va="center", fontsize=10, color=INK)
    cax.text(1.02, 0.5, "high", transform=cax.transAxes, ha="left",
             va="center", fontsize=10, color=INK)
    cax.text(0.5, -1.7, "neural dissimilarity  (1 − r)", transform=cax.transAxes,
             ha="center", va="top", fontsize=10.5, color=INK)

    # titles
    fig.text(0.085, 0.955, "Predicted RDM — MST-Back RSA",
             fontsize=15, fontweight="bold", color=INK, ha="left")
    fig.text(0.085, 0.925,
             "16 conditions:  {A1, B1, A2, B2} × {compared, isolated} × {correct, incorrect}",
             fontsize=10.5, color=MUTED, ha="left")
    fig.text(0.985, 0.955, "PREDICTED", ha="right", va="center", fontsize=8.5,
             fontweight="bold", color="#8A6D3B",
             bbox=dict(boxstyle="round,pad=0.32", fc="#FBF3DD",
                       ec="#E3CE93", lw=0.8))

    fig.savefig(os.path.join(OUT, "mechanisms_RDM_16x16.svg"), format="svg",
                transparent=True)
    if os.environ.get("QA_PNG"):
        qa = "/tmp/qa_svg"; os.makedirs(qa, exist_ok=True)
        fig.savefig(os.path.join(qa, "mechanisms_RDM_16x16.png"), dpi=120,
                    facecolor="white")
    plt.close(fig)
    print("wrote mechanisms_RDM_16x16.svg")


if __name__ == "__main__":
    RDM, labels = build_rdm()
    # sanity print of the two correct blocks (rows/cols A1,B1,A2,B2)
    np.set_printoptions(precision=2, suppress=True)
    print("Correct-Compared block (A1,B1,A2,B2):\n", RDM[0:4, 0:4])
    print("Correct-Isolated block (A1,B1,A2,B2):\n", RDM[4:8, 4:8])
    plot_rdm(RDM, labels)
