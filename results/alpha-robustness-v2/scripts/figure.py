"""Manuscript figure (7.2 in wide, three panels) for one model and one run set.
a) mean combined relative error at alpha = 0.1; b) fraction of time with 0..P relevant processes (alpha = 0.1);
c) alpha sweep: error distribution (top) and fraction of time with all processes relevant (bottom).
usage: figure.py <standard|cc> <completed|stopped>   ->  figures/<model>_<set>_{figure,panel_a,panel_b,panel_c}.{png,pdf,svg}"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
                     "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "legend.fontsize": 6, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False,
                     "axes.spines.right": False, "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans",
                     "mathtext.it": "Liberation Sans:italic", "mathtext.bf": "Liberation Sans:bold", "mathtext.default": "regular"})
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
BLUE, ORANGE, LIGHT = "#2a78d6", "#eb6834", "#cfe0f6"
FLOOR = 1e-8          # errors below this (including exact zeros) are drawn at the floor, labelled "0"
ALPHA_A = 0.1

def _err_for_log(e): return np.maximum(e, FLOOR)

def pow10(k):
    """10^k typeset with the body font (mathtext set to Liberation Sans, an Arial metric clone)."""
    return "1" if k == 0 else rf"$10^{{{int(k)}}}$"

def panel_a(ax, s, P):
    e = s.mean_err.to_numpy(); pos = e[e > 0]; nz = int((e == 0).sum())
    lo = int(np.floor(np.log10(max(pos.min(), FLOOR)))) if len(pos) else -8
    bins = np.arange(lo, 2.25, 0.25)
    ax.hist(np.log10(np.clip(pos, 10.0 ** lo, 10 ** 1.99)), bins=bins, color=BLUE, edgecolor="white", linewidth=0.4)
    ticks = list(range(lo, 3, 2 if 2 - lo > 6 else 1))
    if nz:  # exact zeros: separate bar left of the axis range
        ax.bar(lo - 1, nz, width=0.6, color=MUTED, edgecolor="white", linewidth=0.4)
        ax.set_xticks([lo - 1] + ticks); ax.set_xticklabels(["0"] + [pow10(k) for k in ticks])
        ax.set_xlim(lo - 1.6, 2.2)
    else:
        ax.set_xticks(ticks); ax.set_xticklabels([pow10(k) for k in ticks]); ax.set_xlim(lo - 0.2, 2.2)
    ax.axvline(1, color=ORANGE, lw=0.9)                      # 10 % = observer tolerance alpha = 0.1
    med = np.median(e); ax.axvline(np.log10(max(med, FLOOR)), color=INK2, lw=0.7, ls=(0, (3, 2)))
    top = ax.get_ylim()[1] * 1.45; ax.set_ylim(0, top)
    ax.text(1.06, top * 0.98, "tolerance\n(10%)", color=ORANGE, fontsize=6, va="top", ha="left")
    ax.text(np.log10(max(med, FLOOR)) - 0.08, top * 0.98, f"median\n{med:.2g}%", color=INK2, fontsize=6, va="top", ha="right")
    ax.text(0.02, 0.78, f"{int(np.sum(e <= 10)):,} of {len(e):,} runs\nwithin 10%", transform=ax.transAxes, fontsize=6, color=INK2, va="top")
    ax.set_xlabel("Mean combined relative error (%)"); ax.set_ylabel("Number of runs")

def panel_b(ax, s, P):
    fr = [s[f"frac_{j}_processes"].to_numpy() for j in range(P + 1)]
    parts = ax.violinplot(fr, positions=range(P + 1), widths=0.8, showextrema=False)
    for b in parts["bodies"]: b.set_facecolor(LIGHT); b.set_edgecolor("none"); b.set_alpha(1)
    ax.boxplot(fr, positions=range(P + 1), widths=0.18, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="white", edgecolor=BLUE, lw=0.7), medianprops=dict(color=BLUE, lw=1),
               whiskerprops=dict(color=BLUE, lw=0.6), capprops=dict(color=BLUE, lw=0.6))
    rng = np.random.default_rng(2026); idx = rng.choice(len(s), size=min(200, len(s)), replace=False)
    for j in range(P + 1):
        ax.scatter(j + rng.uniform(-0.28, 0.28, len(idx)), fr[j][idx], s=1.5, color=INK2, alpha=0.45, lw=0, zorder=3)
    ax.set_xticks(range(P + 1)); ax.set_xlim(-0.6, P + 0.6); ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("Number of relevant processes"); ax.set_ylabel("Fraction of time")

def panel_c(ax_top, ax_bot, by_alpha, P):
    alphas = sorted(by_alpha); x = np.log10(alphas)
    data = [_err_for_log(by_alpha[a].mean_err.to_numpy()) for a in alphas]
    ax_top.boxplot([np.log10(d) for d in data], positions=x, widths=0.45, showfliers=True, patch_artist=True,
                   boxprops=dict(facecolor=LIGHT, edgecolor=BLUE, lw=0.6), medianprops=dict(color=BLUE, lw=1),
                   whiskerprops=dict(color=BLUE, lw=0.6), capprops=dict(color=BLUE, lw=0.6),
                   flierprops=dict(marker="o", ms=1.2, mfc=MUTED, mec="none", alpha=0.5))
    means = [by_alpha[a].mean_err.mean() for a in alphas]
    ax_top.plot(x, np.log10(_err_for_log(np.array(means))), "D", ms=2.6, color=INK, zorder=4, label="mean")
    ax_top.plot(x, np.log10(100 * np.array(alphas)), color=ORANGE, lw=0.9, ls=(0, (3, 2)), label="tolerance (100·α %)")
    lo = int(np.log10(FLOOR))
    hi = max(2.0, max(np.log10(d.max()) for d in data), float(np.log10(_err_for_log(np.array(means))).max()))
    step = 2 if hi - lo <= 12 else 4
    yt = list(range(lo, int(np.ceil(hi)) + 1, step))
    ax_top.set_yticks(yt); ax_top.set_yticklabels(["0"] + [pow10(k) for k in yt[1:]])
    ax_top.set_ylim(lo - 0.5, hi + 0.2 * (hi - lo) + 1)  # every value in range, plus room for the legend
    ax_top.set_ylabel("Mean error (%)")
    ax_top.legend(loc="upper left", frameon=False, handlelength=1.4, borderaxespad=0.1, ncol=2, columnspacing=0.8)
    ax_top.tick_params(labelbottom=False)
    m = np.array([by_alpha[a][f"frac_{P}_processes"].mean() for a in alphas])
    sd = np.array([by_alpha[a][f"frac_{P}_processes"].std() for a in alphas])
    ax_bot.errorbar(x, m, yerr=[m - np.clip(m - sd, 0, 1), np.clip(m + sd, 0, 1) - m], fmt="o", ms=2.8, color=ORANGE,
                    ecolor=ORANGE, elinewidth=0.7, capsize=1.5, capthick=0.6)
    ax_bot.set_ylim(-0.05, 1.05); ax_bot.set_yticks([0, 0.5, 1])
    ax_bot.set_ylabel(f"Time with all {P}\nrelevant (fraction)")
    for a_ in (ax_top, ax_bot): a_.set_xlim(x[0] - 0.5, x[-1] + 0.5)
    ax_bot.set_xticks(x); ax_bot.set_xticklabels([pow10(v) for v in x])
    ax_bot.set_xlabel("Observer tolerance, α")

def letter(fig, ax, s, y=0.97):  # panel letters on one baseline, left of each panel's y-axis label
    fig.text(ax.get_position().x0 - 0.055, y, s, fontsize=9, fontweight="bold", va="top", ha="left")

def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{stem}.{ext}", dpi=300 if ext == "png" else None)
    plt.close(fig)
    svg = open(f"{stem}.svg").read()  # name Arial first everywhere, so editors render the text in Arial
    open(f"{stem}.svg", "w").write(svg.replace("font-family: 'Liberation Sans'", "font-family: 'Arial', 'Liberation Sans', sans-serif"))

if __name__ == "__main__":
    model, which = sys.argv[1], sys.argv[2]
    P = len(models.PROCESSES[model])
    sel = pd.read_csv(os.path.join(HERE, "tables", f"{model}_analysis_sets.csv.gz"))
    sel = sel[sel.set == which]
    by_alpha = {a: g for a, g in sel.groupby("alpha")}
    s01 = by_alpha[ALPHA_A]
    stem = os.path.join(HERE, "figures", f"{model}_{which}")
    # combined figure
    fig = plt.figure(figsize=(7.2, 2.5))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 1.1, 1], wspace=0.5, left=0.07, right=0.99, bottom=0.17, top=0.9)
    ax_a, ax_b = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    gc = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[2], hspace=0.12, height_ratios=[1.35, 1])
    ax_c1, ax_c2 = fig.add_subplot(gc[0]), fig.add_subplot(gc[1])
    panel_a(ax_a, s01, P); panel_b(ax_b, s01, P); panel_c(ax_c1, ax_c2, by_alpha, P)
    for ax, l in ((ax_a, "a"), (ax_b, "b"), (ax_c1, "c")): letter(fig, ax, l)
    save(fig, stem + "_figure")
    # individual panels
    for name, w, fn in (("panel_a", 2.4, lambda f: panel_a(f.add_subplot(111), s01, P)),
                        ("panel_b", 2.6, lambda f: panel_b(f.add_subplot(111), s01, P))):
        f = plt.figure(figsize=(w, 2.5)); fn(f); f.subplots_adjust(left=0.2, right=0.97, bottom=0.17, top=0.95); save(f, f"{stem}_{name}")
    f = plt.figure(figsize=(2.4, 2.5)); g2 = GridSpec(2, 1, figure=f, hspace=0.12, height_ratios=[1.35, 1], left=0.27, right=0.97, bottom=0.17, top=0.95)
    panel_c(f.add_subplot(g2[0]), f.add_subplot(g2[1]), by_alpha, P); save(f, f"{stem}_panel_c")
    print("saved", stem + "_*")
