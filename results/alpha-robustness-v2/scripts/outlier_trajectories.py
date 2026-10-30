"""Diagnostic figures for the completed runs above tolerance at alpha = 0.1 (one per run), plus a summary of the
pre-relevance test. Colours follow the manuscript: blue = bacteria and bacterial processes, red = phage and phage
processes, dashed black = transient (Boolean) model.
usage: outlier_trajectories.py  ->  figures/outlier_diagnostics/*.{png,pdf}"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models, figure
from prerelevance import trajectories, ALPHA, LOW
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "figures", "outlier_diagnostics")
BLUE, RED, INK, MID, PALE, ACCENT = "#3a53a4", "#ed2224", "#231f20", "#7f8080", "#e6e6e6", "#e5a41c"
PCOL = {"growth": BLUE, "predation": BLUE, "burst": RED, "decay": RED, "carrying capacity": BLUE}

def run_figure(model, row, zoom=False):
    p = {k: row[k] for k in ["r", "a", "c", "m", "tobs", "b0", "p0"] + (["K"] if model == "cc" else [])}
    t, yf, yb = trajectories(model, p)
    if zoom:  # last ~10 process switches of the transient model: where the cycle-to-cycle drift is largest
        wb_all = models.weights(model, p, yb[0], yb[1], ALPHA) >= 1
        sw = np.where(np.any(wb_all[:, 1:] != wb_all[:, :-1], axis=0))[0]
        i0 = sw[-10] if len(sw) >= 10 else 0
        t, yf, yb = t[i0:], yf[:, i0:], yb[:, i0:]
    wf = models.weights(model, p, yf[0], yf[1], ALPHA); wb = models.weights(model, p, yb[0], yb[1], ALPHA)
    com = (np.abs(yf[0] - yb[0]) / yf[0] + np.abs(yf[1] - yb[1]) / yf[1]) / 2 * 100
    procs = models.PROCESSES[model]; P = len(procs)
    fig = plt.figure(figsize=(7.2, 1.15 * P + 1.6))
    gs = fig.add_gridspec(P, 2, width_ratios=[1.15, 1], hspace=0.12, wspace=0.42, left=0.09, right=0.86, top=0.84, bottom=0.09)
    nleft = 3; gl = gs[:, 0].subgridspec(nleft, 1, hspace=0.12, height_ratios=[1.4, 1, 1])
    ax1, ax2, ax3 = [fig.add_subplot(gl[i]) for i in range(nleft)]
    for y_f, y_b, col, lab in ((yf[0], yb[0], BLUE, "Bacteria"), (yf[1], yb[1], RED, "Phage")):
        ax1.semilogy(t, y_f, color=col, lw=0.9, label=f"{lab} (full)")
        ax1.semilogy(t, y_b, color=INK, lw=0.7, ls=(0, (3, 2)))
    ax1.plot([], [], color=INK, lw=0.7, ls=(0, (3, 2)), label="transient model")
    ax1.set_ylabel("Concentration (ml$^{-1}$)")
    ax1.legend(fontsize=6, frameon=False, ncol=3, loc="lower left", bbox_to_anchor=(0, 1.0), borderaxespad=0.2)
    ax2.semilogy(t, np.maximum(com, 1e-6), color=INK, lw=0.7); ax2.axhline(100 * ALPHA, color=MID, lw=0.8, ls=(0, (3, 2)))
    ax2.text(t[0], 100 * ALPHA * 1.3, "tolerance 10%", fontsize=6, color=MID, ha="left", va="bottom")
    ax2.set_ylabel("Combined error (%)")
    ax3.step(t, (wf >= 1).sum(0), where="post", color=MID, lw=1.1, label="full model")
    ax3.step(t, (wb >= 1).sum(0), where="post", color=INK, lw=0.7, ls=(0, (3, 2)), label="transient model")
    ax3.set_ylim(-0.3, P + 0.5); ax3.set_yticks(range(P + 1)); ax3.set_ylabel("Relevant\nprocesses (#)")
    ax3.legend(fontsize=6, frameon=False, ncol=2, loc="lower left"); ax3.set_xlabel("Time (h)")
    for ax in (ax1, ax2): ax.tick_params(labelbottom=False)
    for j, name in enumerate(procs):
        ax = fig.add_subplot(gs[j, 1])
        ax.axhspan(LOW, 1, color=PALE, lw=0); ax.axhline(1, color=MID, lw=0.6)
        ax.semilogy(t, np.clip(wf[j], 1e-4, 1e4), color=PCOL[name], lw=0.9)
        ax.semilogy(t, np.clip(wb[j], 1e-4, 1e4), color=INK, lw=0.6, ls=(0, (3, 2)))
        pre = ((wb[j] >= LOW) & (wb[j] < 1)).mean()
        ax.set_ylim(1e-4, 1e4); ax.set_yticks([1e-4, 1, 1e4]); ax.set_yticklabels([figure.pow10(-4), "1", figure.pow10(4)])
        ax.text(1.02, 0.5, f"{name}\n{100 * pre:.0f}% pre-relevant", transform=ax.transAxes, fontsize=6, va="center", ha="left",
                color=PCOL[name])
        if j < P - 1: ax.tick_params(labelbottom=False)
        else: ax.set_xlabel("Time (h)")
        if j == 0: ax.set_title("Weights (shaded: pre-relevant 0.1 ≤ w < 1)", fontsize=7, loc="left", pad=3)
    act = wb >= 1; nsw = int(np.sum(np.any(act[:, 1:] != act[:, :-1], axis=0)))
    pars = ", ".join(f"{k} = {v:.3g}" for k, v in p.items())
    what = f"zoom on the last {nsw} switches" if zoom else f"{nsw} switches"
    fig.suptitle(f"{'Carrying-capacity' if model == 'cc' else 'Standard'} model, run {int(row['sample'])}: mean error "
                 f"{row['mean_err']:.3g}% at α = 0.1 ({what})\n{pars}", x=0.02, ha="left", fontsize=7.5)
    return fig

def summary_figure(df):
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4), gridspec_kw=dict(wspace=0.45, left=0.08, right=0.98, bottom=0.2, top=0.86))
    for ax, (c, lab) in zip(axes, [("frac_time_any_pre", "Time with a pre-relevant\nprocess (fraction)"),
                                   ("longest_pre_stay", "Longest pre-relevant stay\n(fraction of window)"),
                                   ("n_switches", "Process switches (#)")]):
        for above, col, z, s in ((False, MID, 2, 2), (True, ACCENT, 3, 10)):
            g = df[df.above_tolerance == above]
            x = g[c] + (1 if c == "n_switches" else 0)
            ax.scatter(x, np.maximum(g.mean_err, 1e-8), s=s, color=col, lw=0, alpha=0.7 if not above else 1, zorder=z,
                       label="above tolerance" if above else "within tolerance")
        ax.set_yscale("log"); ax.axhline(10, color=INK, lw=0.7, ls=(0, (3, 2)))
        if c == "n_switches": ax.set_xscale("log"); ax.set_xlabel("Process switches + 1")
        else: ax.set_xlabel(lab)
        ax.set_ylabel("Mean error (%)")
    axes[0].legend(fontsize=6, frameon=False, loc="lower right")
    fig.suptitle("Pre-relevance test, completed runs at α = 0.1 (both models): dashed line = 10% tolerance", x=0.02, ha="left", fontsize=7.5)
    return fig

def _switches(model, row):
    p = {k: row[k] for k in ["r", "a", "c", "m", "tobs", "b0", "p0"] + (["K"] if model == "cc" else [])}
    t, yf, yb = trajectories(model, p); act = models.weights(model, p, yb[0], yb[1], ALPHA) >= 1
    return int(np.sum(np.any(act[:, 1:] != act[:, :-1], axis=0)))

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    pre = pd.read_csv(os.path.join(HERE, "tables", "prerelevance_alpha0.1.csv"))
    f = summary_figure(pre)
    for ext in ("png", "pdf"): f.savefig(os.path.join(OUT, f"summary_prerelevance_test.{ext}"), dpi=300 if ext == "png" else None)
    plt.close(f)
    with PdfPages(os.path.join(OUT, "all_runs_above_tolerance.pdf")) as pdf:
        for model in ("standard", "cc"):
            sel = pd.read_csv(os.path.join(HERE, "tables", f"{model}_analysis_sets.csv.gz"))
            sel = sel[(sel.alpha == ALPHA) & (sel.set == "completed") & (sel.mean_err > 100 * ALPHA)].sort_values("mean_err", ascending=False)
            for _, row in sel.iterrows():
                fig = run_figure(model, row)
                fig.savefig(os.path.join(OUT, f"{model}_run{int(row['sample'])}.png"), dpi=200); pdf.savefig(fig); plt.close(fig)
                if row["mean_err"] > 0 and _switches(model, row) > 20:
                    fig = run_figure(model, row, zoom=True)
                    fig.savefig(os.path.join(OUT, f"{model}_run{int(row['sample'])}_zoom.png"), dpi=200); pdf.savefig(fig); plt.close(fig)
                print("saved", model, int(row["sample"]), f"{row['mean_err']:.3g}%")
