import pandas as pd, numpy as np, matplotlib, matplotlib.ticker, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
alphas = ["3.7e-14", "0.0001", "0.001", "0.01", "0.1", "1", "10", "12", "13"]  # 3.7e-14: just below alpha* = min_t,process (weight * alpha)
rows = []
for a in alphas:
    d = pd.read_excel(f"{S}/alpha_sweep/alpha_{a}/generated-data/standard_model/results_standard/array_outputs.xlsx")
    rows.append(dict(alpha=float(a), mean_com=d.rel_error_com_per.mean(), median_com=d.rel_error_com_per.median(),
                     max_com=d.rel_error_com_per.max(), mean_bac=d.rel_error_bac_per.mean(), mean_pha=d.rel_error_phage_per.mean()))
t = pd.DataFrame(rows); t.to_csv(f"{S}/figures/scenario1_alpha_sweep_errors.csv", index=False); print(t.to_string(index=False))

INK, INK2, MUTED, GRID, SURF, S1 = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb", "#2a78d6"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK2, "text.color": INK})
YMIN = 3e-5
zero, t = t[t.mean_com == 0].reset_index(drop=True), t[t.mean_com > 0].reset_index(drop=True)  # log axis cannot show 0
# broken x-axis: narrow left panel for alpha* (exact zero error), main panel for the sweep
fig, (ax0, ax) = plt.subplots(1, 2, sharey=True, figsize=(6.8, 4.3), dpi=200, facecolor=SURF,
                              gridspec_kw=dict(width_ratios=[1, 6], wspace=0.06))
for x in zero.alpha:  # drawn hollow on the bottom edge
    ax0.set_facecolor(SURF); ax0.set_xscale("log"); ax0.set_xlim(x / 3, x * 3)
    ax0.plot([x], [YMIN], marker="o", ms=8, mfc=SURF, mec=S1, mew=2, clip_on=False, zorder=4)
    m, e = f"{x:.1e}".split("e")
    ax0.set_xticks([x]); ax0.set_xticklabels([rf"${m}\times10^{{{int(e)}}}$"], fontsize=8)
    ax0.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax0.annotate("α*: error = 0\nall processes\nrelevant at\nall times", (x, YMIN), textcoords="offset points",
                 xytext=(0, 14), ha="center", va="bottom", color=INK2, fontsize=8)
ax0.spines[["top", "right"]].set_visible(False)
ax0.grid(True, axis="y", which="major", color=GRID, lw=0.6); ax0.set_axisbelow(True)
ax.set_facecolor(SURF)
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", which="both", left=False)
kw = dict(color=MUTED, lw=1, clip_on=False, transform=fig.transFigure)
ax.plot(t.alpha, t.mean_com, color=S1, lw=2, marker="o", ms=8, mec=SURF, mew=2, zorder=3)
def fmt(y):
    if 1e-2 <= y < 1e4:
        return f"{y:.3g}%"
    m, e = f"{y:.1e}".split("e")
    return rf"${m}\times10^{{{int(e)}}}$%"
ys = list(t.mean_com)
for i, (x, y) in enumerate(zip(t.alpha, ys)):
    if i > 0 and y == ys[i - 1]:
        continue  # identical value (collapsed model): label it once
    collapsed_next = i + 1 < len(ys) and ys[i + 1] == y
    jump_next = i + 1 < len(ys) and ys[i + 1] > 100 * y
    near_next = i + 1 < len(ys) and t.alpha[i + 1] / x < 1.5  # crowded neighbour to the right
    off, ha, va = (((10, 0), "left", "center") if i <= 1 else
                   ((-2, 8), "right", "bottom") if near_next and not jump_next else
                   ((10, -6), "left", "top") if jump_next else
                   ((-10, 0), "right", "center") if collapsed_next or i == len(ys) - 1 else
                   ((0, 10), "center", "bottom"))
    decade = np.isclose(np.log10(x), round(np.log10(x)))
    label = (fmt(y) + f" (α ≥ {x:g})") if collapsed_next else fmt(y) if decade else f"α = {x:g}: {fmt(y)}"
    if y >= 1e4:  # model collapse (every process weight below 1): only then can the error exceed 10^4 %
        label = f"α = {x:g}: {fmt(y)}\nno process observable"
    ax.annotate(label, (x, y), textcoords="offset points", xytext=off, ha=ha, va=va, color=INK2, fontsize=9)
ax.set_xscale("log"); ax.set_yscale("log")
decades = [x for x in t.alpha if np.isclose(np.log10(x), round(np.log10(x)))]
ax.set_xticks(decades); ax.set_xticklabels([f"{x:g}" for x in decades])
ax.set_xticks([x for x in t.alpha if x not in decades], minor=True)  # tested non-decade alphas: unlabeled ticks
ax.tick_params(axis="x", which="minor", length=4)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_ylim(YMIN, 1e13)
ax.set_xlabel("Observer tolerance, α (reference level)")
ax0.set_ylabel("Mean combined relative error (%)")
fig.suptitle("Scenario 1 (E. coli–T4): Boolean model error vs. α", x=0.02, ha="left", fontsize=11, color=INK)
ax.grid(True, which="major", color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.spines[["top", "right"]].set_visible(False)
ax0.text(0, -0.3, "Mean over 10,000 time points (0–14 h) of the average of bacteria and phage\nrelative errors between the full and Boolean models.",
        transform=ax.transAxes, fontsize=7.5, color=MUTED)
fig.subplots_adjust(left=0.11, right=0.86, top=0.9, bottom=0.25)
for axx, side in ((ax0, 1), (ax, 0)):  # diagonal break marks where the x-axis is cut
    pos = axx.get_position(); xb = pos.x0 + side * pos.width
    fig.add_artist(plt.Line2D([xb - 0.006, xb + 0.006], [pos.y0 - 0.015, pos.y0 + 0.015], **kw))
fig.savefig(f"{S}/figures/scenario1_mean_error_vs_alpha.png", facecolor=SURF)
