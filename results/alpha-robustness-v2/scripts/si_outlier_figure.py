"""Supplementary figures for outlier runs (carrying-capacity model, alpha = 0.1): full-window and zoomed trajectories,
relevant processes and process weights, full vs transient (Boolean) model. Style matches figure.py; colours follow the
manuscript (blue = bacteria / bacterial processes, red = phage / phage processes, dashed black = transient model,
grey band = irrelevant zone w < 1).
usage: si_outlier_figure.py <model> <sample> <stem>  ->  figures/<stem>.{png,pdf,svg}"""
import os, sys, numpy as np, pandas as pd, matplotlib, matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models, figure                                  # figure: shared rcParams, pow10(), save()
from prerelevance import trajectories
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BLUE, RED, INK, MID, PALE = "#3a53a4", "#ed2224", "#231f20", "#7f8080", "#e6e6e6"
PCOL = {"growth": BLUE, "predation": BLUE, "burst": RED, "decay": RED, "carrying capacity": BLUE}
ALPHA, N_SWITCH_ZOOM = 0.1, 10
TITLES = {"a": "Population dynamics (full window)", "b": "Persistent delay (last cycles)",
          "c": "Relevant processes (last cycles)", "d": "Process weights (last cycles)"}

def title(fig, ax, s, gap_pt=6):
    y = ax.get_position().y1 + gap_pt / (fig.get_figheight() * 72); x0 = ax.get_position().x0 - 0.075
    t = fig.text(x0, y, f"{s})", fontsize=8.5, fontweight="bold", va="baseline", color=INK)
    fig.canvas.draw(); bb = t.get_window_extent().transformed(fig.transFigure.inverted())
    fig.text(bb.x1 + 0.004, y, TITLES[s], fontsize=8.5, va="baseline", color=INK)

def populations(ax, t, yf, yb):
    for y_f, y_b, col in ((yf[0], yb[0], BLUE), (yf[1], yb[1], RED)):
        ax.semilogy(t, y_f, color=col, lw=0.8); ax.semilogy(t, y_b, color=INK, lw=0.6, ls=(0, (3, 2)))
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_ylabel("Concentration (ml$^{-1}$)"); ax.set_xlabel("Time (h)")

if __name__ == "__main__":
    model, sample, stem = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    sel = pd.read_csv(os.path.join(HERE, "tables", f"{model}_analysis_sets.csv.gz"))
    row = sel[(sel["sample"] == sample) & (sel.alpha == ALPHA) & (sel.set == "completed")].iloc[0]
    p = {k: row[k] for k in ["r", "a", "c", "m", "tobs", "b0", "p0"] + (["K"] if model == "cc" else [])}
    t, yf, yb = trajectories(model, p)
    wf = models.weights(model, p, yf[0], yf[1], ALPHA); wb = models.weights(model, p, yb[0], yb[1], ALPHA)
    sw = np.where(np.any((wb >= 1)[:, 1:] != (wb >= 1)[:, :-1], axis=0))[0]
    i0 = sw[-N_SWITCH_ZOOM] if len(sw) >= N_SWITCH_ZOOM else 0
    procs = models.PROCESSES[model]; P = len(procs)

    fig = plt.figure(figsize=(7.2, 5.0))
    gs = GridSpec(2, 2, figure=fig, left=0.09, right=0.85, bottom=0.08, top=0.93, wspace=0.32, hspace=0.42)
    ax_a, ax_b, ax_c = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0])
    populations(ax_a, t, yf, yb)
    ax_a.plot([], [], color=BLUE, lw=0.8, label="Bacteria"); ax_a.plot([], [], color=RED, lw=0.8, label="Phage")
    ax_a.plot([], [], color=INK, lw=0.6, ls=(0, (3, 2)), label="Transient model")
    ax_a.legend(fontsize=6, frameon=False, ncol=3, loc="lower left", bbox_to_anchor=(0, -0.02), handlelength=1.6, columnspacing=0.9)
    lo, hi = ax_a.get_ylim(); ax_a.set_ylim(lo / 1e3, hi)                       # room for the legend below the data
    ax_a.axvspan(t[i0], t[-1], color=PALE, lw=0, zorder=0)                     # zoom window shown in b-d
    populations(ax_b, t[i0:], yf[:, i0:], yb[:, i0:])
    ax_c.step(t[i0:], (wf[:, i0:] >= 1).sum(0), where="post", color=MID, lw=1.2, label="Full model")
    ax_c.step(t[i0:], (wb[:, i0:] >= 1).sum(0), where="post", color=INK, lw=0.7, ls=(0, (3, 2)), label="Transient model")
    ax_c.set_ylim(-0.3, P + 0.6); ax_c.set_yticks(range(P + 1)); ax_c.set_ylabel("Relevant processes (#)"); ax_c.set_xlabel("Time (h)")
    ax_c.legend(fontsize=6, frameon=False, ncol=2, loc="lower left")
    gd = gs[1, 1].subgridspec(P, 1, hspace=0.15)
    ax_d = []
    for j, name in enumerate(procs):
        ax = fig.add_subplot(gd[j]); ax_d.append(ax)
        ax.axhspan(1e-4, 1, color=PALE, lw=0); ax.axhline(1, color=MID, lw=0.5)
        ax.semilogy(t[i0:], np.clip(wf[j, i0:], 1e-4, 1e4), color=PCOL[name], lw=0.8)
        ax.semilogy(t[i0:], np.clip(wb[j, i0:], 1e-4, 1e4), color=INK, lw=0.55, ls=(0, (3, 2)))
        ax.set_ylim(1e-4, 1e4); ax.set_yticks([1e-3, 1, 1e3]); ax.set_yticklabels([figure.pow10(-3), "1", figure.pow10(3)], fontsize=5.5)
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.text(1.02, 0.5, name.replace("carrying capacity", "carrying\ncapacity"), transform=ax.transAxes, fontsize=6,
                color=PCOL[name], va="center", ha="left")
        if j < P - 1: ax.tick_params(labelbottom=False)
        else: ax.set_xlabel("Time (h)")
    fig.text(ax_d[0].get_position().x0 - 0.06, (ax_d[0].get_position().y1 + ax_d[-1].get_position().y0) / 2, "Weight (w)",
             rotation=90, va="center", ha="center", fontsize=7)
    for ax, l in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c"), (ax_d[0], "d")): title(fig, ax, l)
    figure.save(fig, os.path.join(HERE, "figures", stem))
    print(f"saved {stem}: sample {sample}, mean error {row.mean_err:.3g}%, zoom from t = {t[i0]:.0f} h of {p['tobs']:.0f} h,",
          ", ".join(f"{k} = {v:.3g}" for k, v in p.items()))
