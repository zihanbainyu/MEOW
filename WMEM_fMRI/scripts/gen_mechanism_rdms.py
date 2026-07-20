#!/usr/bin/env python3
"""
Desired (predicted) RDM per mechanism for the MST-Back RSA.

Each mechanism is the correct-vs-incorrect contrast on one B1 pairing:
    pattern separation   sim(A1·B1)   correct < incorrect   (more separated when correct)
    pattern completion   sim(B2·B1)   correct > incorrect   (reinstated when correct)
    predictive recall    sim(A2·B1)   correct > incorrect   (pair-mate recalled when correct)

Each RDM spans {B1, X} x {compared, isolated} x {correct, incorrect} = 8 conditions,
ordered within each condition block as [B1_cor, X_cor, B1_inc, X_inc] so the diagnostic
correct-B1X and incorrect-B1X cells sit adjacent to the diagonal for direct reading.

Style: lab minimal — Purples, continuous cells, black block dividers, labels top+left,
Helvetica Light, bottom gradient bar.
"""
import os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

mpl.rcParams.update({
    "svg.fonttype": "none",
    "font.family": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.weight": "light",
    "text.color": "#26262f",
})

OUT = "/Users/bai/Documents/GitHub/MEOW/Experiment_2/expected_results"
os.makedirs(OUT, exist_ok=True)

INK  = "#26262f"
COMP = "#8e7cc3"   # compared
ISO  = "#5b8fc9"   # isolated
PURPLES = LinearSegmentedColormap.from_list(
    "labpurple", ["#FFFFFF", "#EDE9F5", "#CFC3E6", "#A996D2", "#7B62B3", "#4B3B86"])

# diagnostic B1.X neural SIMILARITY  [condition][accuracy]   (dark = high similarity)
DIAG = {
    "separation": {"comp": {"cor": 0.12, "inc": 0.70},   # correct LESS similar (separated)
                   "iso":  {"cor": 0.45, "inc": 0.50}},   # flat (comparison-driven)
    "completion": {"comp": {"cor": 0.82, "inc": 0.18},   # correct MORE similar
                   "iso":  {"cor": 0.82, "inc": 0.18}},   # both conditions
    "recall":     {"comp": {"cor": 0.82, "inc": 0.18},   # correct MORE similar
                   "iso":  {"cor": 0.45, "inc": 0.45}},   # needs a studied pair-mate
}
MECHS = [
    ("separation", "Pattern separation", "A1", "A1 · B1"),
    ("completion", "Pattern completion", "B2", "B2 · B1"),
    ("recall",     "Predictive recall",  "A2", "A2 · B1"),
]

# condition order: cond -> acc -> item ; within block [B1_cor, X_cor, B1_inc, X_inc]
CONDS = []
for cond in ["comp", "iso"]:
    for acc in ["cor", "inc"]:
        for item in ["B1", "X"]:
            CONDS.append((item, cond, acc))


def build(mech):
    """Return an 8x8 neural SIMILARITY matrix (diagonal = 1, dark = high r)."""
    n = 8
    S = np.ones((n, n))                            # diagonal = self-similarity = 1
    dg = DIAG[mech]
    for i, (it_i, cd_i, ac_i) in enumerate(CONDS):
        for j, (it_j, cd_j, ac_j) in enumerate(CONDS):
            if i == j:
                continue
            if it_i == it_j:                       # same item (B1-B1 or X-X): high similarity
                s = 0.85
                if ac_i != ac_j:
                    s -= 0.12
                if cd_i != cd_j:
                    s -= 0.18
            else:                                  # B1 vs X
                if cd_i == cd_j and ac_i == ac_j:  # diagnostic cell
                    s = dg[cd_i][ac_i]
                else:                              # non-diagnostic cross
                    s = 0.38 - (0.06 if cd_i != cd_j else 0.0)
            S[i, j] = min(max(s, 0.0), 1.0)
    return S


def draw(mech_key, title, xitem, pair):
    n = 8
    M = build(mech_key)
    items = ["B1", xitem, "B1", xitem, "B1", xitem, "B1", xitem]

    fig = plt.figure(figsize=(4.4, 5.0))
    ax = fig.add_axes([0.20, 0.17, 0.67, 0.58])
    ax.imshow(M, cmap=PURPLES, vmin=0, vmax=1.0, interpolation="nearest")

    # black block dividers: thin at accuracy (2,6), thick at condition (4)
    for b, lw in [(2, 1.1), (4, 2.4), (6, 1.1)]:
        ax.axhline(b - 0.5, color=INK, lw=lw)
        ax.axvline(b - 0.5, color=INK, lw=lw)
    for sp in ax.spines.values():
        sp.set_edgecolor(INK); sp.set_linewidth(1.4)

    # item ticks (top + left)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(items, fontsize=9); ax.set_yticklabels(items, fontsize=9)
    ax.xaxis.set_ticks_position("top")
    ax.tick_params(length=0)

    tf_x = ax.get_xaxis_transform()   # x=data, y=axes fraction
    tf_y = ax.get_yaxis_transform()   # y=data, x=axes fraction

    # accuracy labels (black) centered over each 2-block
    for c, lab in [(0.5, "Correct"), (2.5, "Incorrect"),
                   (4.5, "Correct"), (6.5, "Incorrect")]:
        ax.text(c, 1.09, lab, transform=tf_x, ha="center", va="bottom",
                fontsize=8.5, color=INK, clip_on=False)
        ax.text(-0.11, c, lab, transform=tf_y, ha="right", va="center",
                fontsize=8.5, color=INK, rotation=90, clip_on=False)
    # condition labels (colored) centered over each 4-block
    for c, lab, col in [(1.5, "Compared", COMP), (5.5, "Isolated", ISO)]:
        ax.text(c, 1.20, lab, transform=tf_x, ha="center", va="bottom",
                fontsize=11, color=col, clip_on=False)
        ax.text(-0.22, c, lab, transform=tf_y, ha="right", va="center",
                fontsize=11, color=col, rotation=90, clip_on=False)

    # minimal mechanism label (top-left, above the group labels)
    fig.text(0.20, 0.965, title, fontsize=12.5, color=INK, ha="left", va="top")
    fig.text(0.20, 0.928, pair, fontsize=10, color="#6b6b78", ha="left", va="top")

    # bottom gradient colorbar
    cax = fig.add_axes([0.32, 0.075, 0.42, 0.026])
    cax.imshow(np.linspace(0, 1, 256)[None, :], aspect="auto", cmap=PURPLES)
    cax.set_xticks([]); cax.set_yticks([])
    for sp in cax.spines.values():
        sp.set_edgecolor(INK); sp.set_linewidth(1.0)
    cax.text(-0.03, 0.5, "low", transform=cax.transAxes, ha="right", va="center",
             fontsize=9.5, color=INK)
    cax.text(1.03, 0.5, "high", transform=cax.transAxes, ha="left", va="center",
             fontsize=9.5, color=INK)
    cax.text(0.5, -1.9, "neural similarity  (r)", transform=cax.transAxes,
             ha="center", va="top", fontsize=9.5, color=INK)

    name = f"mech_rdm_{mech_key}.svg"
    fig.savefig(os.path.join(OUT, name), format="svg", transparent=True)
    if os.environ.get("QA_PNG"):
        qa = "/tmp/qa_svg"; os.makedirs(qa, exist_ok=True)
        fig.savefig(os.path.join(qa, name.replace(".svg", ".png")), dpi=150,
                    facecolor="white")
    plt.close(fig)
    print("wrote", name, " diag(comp cor,inc)=", DIAG[mech_key]["comp"])


if __name__ == "__main__":
    for k, t, x, p in MECHS:
        draw(k, t, x, p)
