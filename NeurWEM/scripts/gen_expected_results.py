#!/usr/bin/env python3
"""
Generate expected-results (predicted) figures for the MST-Back fMRI study.
One editable vector SVG per research question. Styled for a talk.

Palette follows the behavioral task legend:
    isolation (A)  -> blue      comparison (B) -> purple     novel -> green

Every panel carries a "Predicted" tag so the audience reads them as
hypotheses, not data.
"""
import os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from matplotlib.lines import Line2D

# ---- global style ---------------------------------------------------------
mpl.rcParams.update({
    "svg.fonttype": "none",          # keep text as editable <text> in SVG
    "font.family": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 11,
    "axes.edgecolor": "#3A3A46",
    "axes.linewidth": 1.1,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.color": "#3A3A46",
    "ytick.color": "#3A3A46",
    "text.color": "#26262F",
    "axes.labelcolor": "#26262F",
    "figure.dpi": 100,
})

INK      = "#26262F"
MUTED    = "#6B6B78"
GRID     = "#E7E7EE"
ISO      = "#5B8FC9"   # isolation / A
ISO_L    = "#AEC9E6"
COMP     = "#8E7CC3"   # comparison / B
COMP_L   = "#CBBFE4"
NOVEL    = "#7FA97F"   # novel
NOVEL_L  = "#C0D6C0"
WARM     = "#C7603F"   # predicted-effect accent (lines / fits)
GREY     = "#A7ABB3"
GREY_L   = "#D6D8DD"

OUT = "/Users/bai/Documents/GitHub/MEOW/Experiment_2/expected_results"
os.makedirs(OUT, exist_ok=True)


def predicted_tag(ax):
    """Small 'Predicted' pill, top-right of the axes."""
    ax.annotate(
        "PREDICTED", xy=(1.0, 1.045), xycoords="axes fraction",
        ha="right", va="bottom", fontsize=8, fontweight="bold",
        color="#8A6D3B",
        bbox=dict(boxstyle="round,pad=0.32", fc="#FBF3DD", ec="#E3CE93", lw=0.8),
    )


def fig_tag(fig):
    """Figure-level 'Predicted' pill (for multi-panel figures)."""
    fig.text(0.985, 0.905, "PREDICTED", ha="right", va="top", fontsize=8,
             fontweight="bold", color="#8A6D3B",
             bbox=dict(boxstyle="round,pad=0.32", fc="#FBF3DD",
                       ec="#E3CE93", lw=0.8))


def titles(fig, main, sub):
    fig.text(0.015, 0.965, main, fontsize=13.5, fontweight="bold",
             color=INK, ha="left", va="top")
    fig.text(0.015, 0.905, sub, fontsize=10, color=MUTED,
             ha="left", va="top")


def ygrid(ax):
    ax.yaxis.grid(True, color=GRID, lw=1.0, zorder=0)
    ax.set_axisbelow(True)


def sig(ax, x1, x2, y, text, lw=1.2, h=None):
    """significance bracket"""
    if h is None:
        h = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.03
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color=INK, lw=lw)
    ax.text((x1 + x2) / 2, y + h * 1.05, text, ha="center", va="bottom",
            fontsize=12, color=INK)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), format="svg",
                transparent=True, bbox_inches="tight")
    if os.environ.get("QA_PNG"):
        qa = "/tmp/qa_svg"; os.makedirs(qa, exist_ok=True)
        fig.savefig(os.path.join(qa, name.replace(".svg", ".png")),
                    dpi=120, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ===========================================================================
# Q1 (primary) - regional BOLD: EM available vs unavailable during
# successful 'similar' discrimination
# ===========================================================================
def fig1():
    fig, ax = plt.subplots(figsize=(6.1, 4.3))
    fig.subplots_adjust(top=0.80, left=0.12, right=0.97, bottom=0.20)
    rois = ["Hippocampus", "Perirhinal", "Parahipp.", "mPFC", "Early visual\n(control)"]
    prev  = np.array([0.42, 0.31, 0.23, 0.26, 0.05])
    novel = np.array([0.18, 0.15, 0.16, 0.13, 0.06])
    e = np.array([0.05, 0.05, 0.05, 0.05, 0.04])
    x = np.arange(len(rois)); w = 0.38
    ygrid(ax)
    ax.bar(x - w/2, prev,  w, yerr=e, color=COMP,  edgecolor="white",
           label="Previously encountered  (EM available)",
           error_kw=dict(ecolor=INK, lw=1.1, capsize=3), zorder=3)
    ax.bar(x + w/2, novel, w, yerr=e, color=NOVEL, edgecolor="white",
           label="Novel  (EM unavailable)",
           error_kw=dict(ecolor=INK, lw=1.1, capsize=3), zorder=3)
    for i, s in enumerate(["***", "**", "*", "*", "n.s."]):
        yy = max(prev[i], novel[i]) + e[i] + 0.02
        ax.text(x[i], yy, s, ha="center", va="bottom", fontsize=11,
                color=INK if s != "n.s." else MUTED)
    ax.set_xticks(x); ax.set_xticklabels(rois, fontsize=9.5)
    ax.set_ylabel("BOLD response during correct\n'similar' discrimination (a.u.)")
    ax.set_ylim(0, 0.56)
    ax.legend(frameon=False, fontsize=9, loc="upper right",
              bbox_to_anchor=(1.0, 1.02))
    predicted_tag(ax)
    titles(fig,
           "Episodic availability elevates MTL activity during lure discrimination",
           "Q1  •  Which regions favor 'similar' discrimination when EM can support WM?")
    save(fig, "Q1_regions_EM_available_vs_novel.svg")


# ===========================================================================
# Q2a (primary) - pattern separation: comparison differentiates similar items,
# and differentiation predicts later discrimination
# ===========================================================================
def fig2a():
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.6, 4.3),
                                   gridspec_kw=dict(wspace=0.38))
    fig.subplots_adjust(top=0.70, left=0.085, right=0.975, bottom=0.17)

    # Panel A: within-pair neural similarity, isolated vs compared
    ygrid(axA)
    conds = ["Isolated\n(A1)", "Compared\n(B1)"]
    vals = [0.55, 0.31]; e = [0.05, 0.05]
    axA.bar([0, 1], vals, 0.55, yerr=e, color=[ISO, COMP], edgecolor="white",
            error_kw=dict(ecolor=INK, lw=1.1, capsize=3), zorder=3)
    axA.set_xticks([0, 1]); axA.set_xticklabels(conds, fontsize=10)
    axA.set_ylabel("Neural similarity between\nsimilar pair-mates (r)")
    axA.set_ylim(0, 0.72)
    axA.set_xlim(-0.6, 1.6)
    sig(axA, 0, 1, 0.63, "**")
    axA.set_title("Comparison differentiates\nsimilar representations",
                  fontsize=10.5, color=INK, pad=8)

    # Panel B: differentiation vs discrimination (LDI)
    ygrid(axB)
    rng = np.random.default_rng(7)
    n = 30
    diff = rng.normal(0.24, 0.10, n)
    ldi = 0.28 + 1.25 * (diff - diff.mean()) + rng.normal(0, 0.10, n)
    axB.scatter(diff, ldi, s=42, color=COMP, edgecolor="white", lw=0.8,
                alpha=0.9, zorder=3)
    b, a = np.polyfit(diff, ldi, 1)
    xs = np.linspace(diff.min(), diff.max(), 100)
    axB.plot(xs, a + b * xs, color=WARM, lw=2.2, zorder=4)
    axB.fill_between(xs, a + b*xs - 0.08, a + b*xs + 0.08, color=WARM, alpha=0.12)
    axB.set_xlabel("Neural differentiation of pair\n(isolated − compared similarity)")
    axB.set_ylabel("Later discrimination (LDI)")
    axB.text(0.05, 0.92, "r ≈ 0.55", transform=axB.transAxes,
             fontsize=11, color=WARM, fontweight="bold", va="top")
    axB.set_title("Differentiation supports\nlater discrimination",
                  fontsize=10.5, color=INK, pad=8)

    fig_tag(fig)
    titles(fig,
           "Pattern separation: encoding-phase comparison sharpens similar representations",
           "Q2a  •  Does comparing similar items differentiate them, and does that aid discrimination?")
    save(fig, "Q2a_pattern_separation.svg")


# ===========================================================================
# Q2b (primary) - pattern completion: studied items reinstate their encoding
# pattern during WM discrimination (stronger when correct)
# ===========================================================================
def fig2b():
    fig, ax = plt.subplots(figsize=(6.2, 4.3))
    fig.subplots_adjust(top=0.80, left=0.12, right=0.97, bottom=0.15)
    conds = ["Prev. isolated\n(A2)", "Prev. compared\n(B2)", "Novel"]
    correct   = np.array([0.27, 0.34, 0.04])
    incorrect = np.array([0.11, 0.13, 0.03])
    e = np.array([0.04, 0.04, 0.03])
    x = np.arange(3); w = 0.38
    ygrid(ax)
    ax.axhline(0, color="#B9BCC3", lw=1.0, zorder=1)
    cols = [ISO, COMP, NOVEL]
    ax.bar(x - w/2, correct, w, yerr=e, color=cols, edgecolor="white",
           error_kw=dict(ecolor=INK, lw=1.1, capsize=3), zorder=3)
    ax.bar(x + w/2, incorrect, w, yerr=e, color=cols, edgecolor="white",
           hatch="////", alpha=0.55,
           error_kw=dict(ecolor=INK, lw=1.1, capsize=3), zorder=3)
    for i, s in enumerate(["**", "***", "n.s."]):
        yy = correct[i] + e[i] + 0.015
        ax.text(x[i] - w/2, yy, s, ha="center", va="bottom", fontsize=11,
                color=INK if s != "n.s." else MUTED)
    ax.set_xticks(x); ax.set_xticklabels(conds, fontsize=9.5)
    ax.set_ylabel("Encoding→retrieval reinstatement (r)")
    ax.set_ylim(-0.03, 0.44)
    leg = [Line2D([0],[0], marker="s", color="w", markerfacecolor=GREY,
                  markersize=11, label="Correct discrimination"),
           Line2D([0],[0], marker="s", color="w", markerfacecolor="w",
                  markeredgecolor=GREY, markersize=11, label="Incorrect")]
    # represent hatch in legend
    leg[1] = Line2D([0],[0], marker="s", color="w", markerfacecolor=GREY_L,
                    markersize=11, label="Incorrect")
    ax.legend(handles=leg, frameon=False, fontsize=9, loc="upper right")
    predicted_tag(ax)
    titles(fig,
           "Pattern completion: studied items reinstate their encoding pattern",
           "Q2b  •  Is the encoding representation restored from a partial cue during WM discrimination?")
    save(fig, "Q2b_pattern_completion.svg")


# ===========================================================================
# Q2c (primary) - predictive recall: the similar pair-mate is reinstated
# BEFORE it reappears (anticipatory), and this precedes correct discrimination
# ===========================================================================
def fig2c():
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    fig.subplots_adjust(top=0.80, left=0.12, right=0.97, bottom=0.15)
    t = np.linspace(-3.0, 2.5, 200)

    def bump(center, amp, width, base=0.02):
        return base + amp * np.exp(-0.5 * ((t - center) / width) ** 2)

    correct = bump(-0.5, 0.30, 0.9)          # anticipatory rise before onset
    incorrect = bump(0.4, 0.10, 1.0, 0.02)   # weak, only after onset
    ygrid(ax)
    # anticipatory window shading
    ax.axvspan(-2.4, 0, color=COMP, alpha=0.07, zorder=0)
    ax.text(-1.2, 0.345, "anticipatory\nwindow", ha="center", va="top",
            fontsize=9, color=COMP)
    ax.axvline(0, color="#7C7F87", ls="--", lw=1.3, zorder=2)
    ax.text(0.06, 0.36, "pair-mate\nre-enters WM", fontsize=9, color=MUTED,
            va="top", ha="left")
    ax.plot(t, correct, color=WARM, lw=2.6, zorder=4, label="Correct discrimination")
    ax.plot(t, incorrect, color=GREY, lw=2.2, ls=(0, (4, 2)), zorder=3,
            label="Incorrect")
    ax.axhline(0, color="#B9BCC3", lw=1.0)
    ax.set_xlabel("Time relative to pair-mate onset (s)")
    ax.set_ylabel("Pair-mate reinstatement evidence")
    ax.set_ylim(-0.02, 0.40)
    ax.set_xlim(-3.0, 2.5)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    predicted_tag(ax)
    titles(fig,
           "Predictive recall: the pair-mate is reinstated before it reappears",
           "Q2c  •  Is the similar pair-mate predictively recalled from EM to guide upcoming discrimination?")
    save(fig, "Q2c_predictive_recall.svg")


# ===========================================================================
# Q3 (secondary) - HPC-PFC coupling tracks the individual EM->WM benefit
# ===========================================================================
def fig3():
    fig, ax = plt.subplots(figsize=(6.0, 4.3))
    fig.subplots_adjust(top=0.80, left=0.13, right=0.96, bottom=0.15)
    rng = np.random.default_rng(3)
    n = 28
    benefit = rng.normal(0.15, 0.07, n)
    coupling = 0.20 + 1.6 * (benefit - benefit.mean()) + rng.normal(0, 0.07, n)
    ygrid(ax)
    ax.scatter(benefit, coupling, s=46, color=COMP, edgecolor="white",
               lw=0.8, alpha=0.9, zorder=3)
    b, a = np.polyfit(benefit, coupling, 1)
    xs = np.linspace(benefit.min(), benefit.max(), 100)
    ax.plot(xs, a + b*xs, color=WARM, lw=2.3, zorder=4)
    ax.fill_between(xs, a+b*xs-0.06, a+b*xs+0.06, color=WARM, alpha=0.12)
    ax.set_xlabel("EM→WM behavioral benefit\n(Δ discrimination: prev. encountered − novel)")
    ax.set_ylabel("Hippocampal–prefrontal coupling\nduring WM discrimination (PPI β)")
    ax.text(0.05, 0.93, "r ≈ 0.52", transform=ax.transAxes, fontsize=11,
            color=WARM, fontweight="bold", va="top")
    predicted_tag(ax)
    titles(fig,
           "Hippocampal–prefrontal coupling scales with the episodic benefit",
           "Q3  •  Does HPC–PFC coupling track how much performance benefits from available EM?")
    save(fig, "Q3_hpc_pfc_coupling.svg")


# ===========================================================================
# Q4 (exploratory, CPM) - distributed connectivity predicts the EM->WM benefit
# ===========================================================================
def fig4():
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(8.8, 4.4),
                                   gridspec_kw=dict(width_ratios=[1.05, 1]))
    fig.subplots_adjust(top=0.70, left=0.085, right=0.975, bottom=0.16, wspace=0.30)

    # Panel A: predicted vs observed (leave-one-out CPM)
    rng = np.random.default_rng(11)
    n = 28
    obs = rng.normal(0.15, 0.08, n)
    pred = 0.15 + 0.55 * (obs - obs.mean()) + rng.normal(0, 0.055, n)
    ygrid(axA)
    axA.scatter(obs, pred, s=44, color=COMP, edgecolor="white", lw=0.8,
                alpha=0.9, zorder=3)
    lim = [min(obs.min(), pred.min()) - 0.02, max(obs.max(), pred.max()) + 0.02]
    axA.plot(lim, lim, ls="--", color=GREY, lw=1.2, zorder=2)
    b, a = np.polyfit(obs, pred, 1)
    xs = np.linspace(obs.min(), obs.max(), 100)
    axA.plot(xs, a + b*xs, color=WARM, lw=2.2, zorder=4)
    axA.set_xlim(lim); axA.set_ylim(lim)
    axA.set_xlabel("Observed EM→WM benefit")
    axA.set_ylabel("CPM-predicted benefit\n(cross-validated)")
    axA.text(0.05, 0.93, "r ≈ 0.43\np < .01", transform=axA.transAxes,
             fontsize=10.5, color=WARM, fontweight="bold", va="top")
    axA.set_title("Connectome predicts the benefit", fontsize=10.5, pad=8)

    # Panel B: schematic positive/negative networks (circular node layout)
    axB.set_aspect("equal"); axB.axis("off")
    axB.set_xlim(-1.35, 1.35); axB.set_ylim(-1.35, 1.45)
    nodes = ["HPC", "PFC\n(FPN)", "Parietal", "DMN", "Visual"]
    ang = np.linspace(90, 90 - 360, len(nodes) + 1)[:-1] * np.pi / 180
    pos = {nm: (np.cos(a), np.sin(a)) for nm, a in zip(nodes, ang)}
    pos_edges = [("HPC", "PFC\n(FPN)"), ("HPC", "Parietal"), ("PFC\n(FPN)", "Parietal")]
    neg_edges = [("DMN", "Visual"), ("DMN", "HPC")]
    for u, v in pos_edges:
        axB.plot(*zip(pos[u], pos[v]), color=WARM, lw=3.0, alpha=0.85, zorder=2)
    for u, v in neg_edges:
        axB.plot(*zip(pos[u], pos[v]), color=ISO, lw=3.0, alpha=0.85, zorder=2)
    for nm, (xx, yy) in pos.items():
        axB.add_patch(Circle((xx, yy), 0.24, fc="white", ec=INK, lw=1.4, zorder=3))
        axB.text(xx, yy, nm, ha="center", va="center", fontsize=8.5,
                 color=INK, zorder=4)
    axB.plot([], [], color=WARM, lw=3, label="Positive network (+)")
    axB.plot([], [], color=ISO, lw=3, label="Negative network (−)")
    axB.legend(frameon=False, fontsize=9, loc="lower center",
               bbox_to_anchor=(0.5, -0.12), ncol=2, handlelength=1.4)
    axB.set_title("Edges predicting the benefit", fontsize=10.5, pad=2)

    fig_tag(fig)
    titles(fig,
           "A distributed connectivity pattern predicts the episodic benefit",
           "Q4  •  Exploratory connectome-based predictive modeling (CPM; Finn et al., 2015)")
    save(fig, "Q4_cpm_connectome.svg")


# ===========================================================================
# Q5 (confirmatory, Kirwan & Stark 2007) - MTL/parietal old > new at recognition
# ===========================================================================
def fig5():
    fig, ax = plt.subplots(figsize=(6.2, 4.3))
    fig.subplots_adjust(top=0.80, left=0.12, right=0.97, bottom=0.18)
    rois = ["Hippocampus", "Perirhinal", "Lateral\nparietal", "Precuneus", "mPFC"]
    hits = np.array([0.36, 0.30, 0.40, 0.34, 0.22])
    cr   = np.array([0.16, 0.16, 0.18, 0.16, 0.14])
    e = np.array([0.05]*5)
    x = np.arange(5); w = 0.38
    ygrid(ax)
    ax.bar(x - w/2, hits, w, yerr=e, color=COMP, edgecolor="white",
           label="Hits (old)", error_kw=dict(ecolor=INK, lw=1.1, capsize=3), zorder=3)
    ax.bar(x + w/2, cr, w, yerr=e, color=GREY, edgecolor="white",
           label="Correct rejections (new)",
           error_kw=dict(ecolor=INK, lw=1.1, capsize=3), zorder=3)
    for i, s in enumerate(["***", "**", "***", "**", "*"]):
        yy = max(hits[i], cr[i]) + e[i] + 0.015
        ax.text(x[i], yy, s, ha="center", va="bottom", fontsize=11, color=INK)
    ax.set_xticks(x); ax.set_xticklabels(rois, fontsize=9.5)
    ax.set_ylabel("BOLD response at recognition (a.u.)")
    ax.set_ylim(0, 0.52)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    predicted_tag(ax)
    titles(fig,
           "MTL and parietal cortex signal old items above novel foils",
           "Q5  •  Confirmatory old > new recognition contrast (Kirwan & Stark, 2007)")
    save(fig, "Q5_recognition_old_vs_new.svg")


# ===========================================================================
# Q6 (exploratory) - item recognition strength predicts in-task reinstatement
# ===========================================================================
def fig6():
    fig, ax = plt.subplots(figsize=(6.0, 4.3))
    fig.subplots_adjust(top=0.80, left=0.13, right=0.96, bottom=0.15)
    rng = np.random.default_rng(21)
    n = 40
    strength = rng.normal(0.0, 1.0, n)
    reinst = 0.20 + 0.09 * strength + rng.normal(0, 0.09, n)
    ygrid(ax)
    ax.scatter(strength, reinst, s=42, color=COMP, edgecolor="white", lw=0.8,
               alpha=0.9, zorder=3)
    b, a = np.polyfit(strength, reinst, 1)
    xs = np.linspace(strength.min(), strength.max(), 100)
    ax.plot(xs, a + b*xs, color=WARM, lw=2.3, zorder=4)
    ax.fill_between(xs, a+b*xs-0.07, a+b*xs+0.07, color=WARM, alpha=0.12)
    ax.set_xlabel("Item recognition strength (memory-strength / confidence, z)")
    ax.set_ylabel("In-task episodic reinstatement (r)")
    ax.text(0.05, 0.93, "r ≈ 0.38", transform=ax.transAxes, fontsize=11,
            color=WARM, fontweight="bold", va="top")
    predicted_tag(ax)
    titles(fig,
           "Recognition strength tracks in-task episodic reinstatement",
           "Q6  •  Does item-level recognition strength predict that item's reinstatement magnitude?")
    save(fig, "Q6_recognition_predicts_reinstatement.svg")


if __name__ == "__main__":
    fig1(); fig2a(); fig2b(); fig2c(); fig3(); fig4(); fig5(); fig6()
    print("\nAll figures written to:", OUT)
