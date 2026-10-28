import pandas as pd, numpy as np, matplotlib, os, sys
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lhs_tally_stop import classify
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 2026
N = sys.argv[2] if len(sys.argv) > 2 else "1000"
TITLE_N = "500 completed + 500 stopped runs" if N.startswith("1000balanced") else f"{int(N):,} Latin hypercube samples"
SUB_N = "Random draw from 2,000 Latin hypercube samples; " if N.startswith("1000balanced") else ""
df = classify(pd.read_csv(f"{S}/figures/lhs_b0p0min1_stop_alpha0.1_n{N}_seed{seed}.csv"))
df.to_csv(f"{S}/figures/lhs_b0p0min1_stop_alpha0.1_n{N}_seed{seed}_classified.csv", index=False)
INK, INK2, MUTED, GRID, SURF, BLUE, AQUA = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb", "#2a78d6", "#1baf7a"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK2, "text.color": INK})
plot = df[(df.mean_com > 0) & df.stop_group.isin(["completed full window", "stopped: below 10⁻³/ml"])]
groups = [("Completed the full window", BLUE, plot[plot.stop_group == "completed full window"]),
          ("Stopped when a population fell below 10⁻³/ml", AQUA, plot[plot.stop_group == "stopped: below 10⁻³/ml"])]
v = plot.mean_com
LO, HI = int(np.floor(np.log10(v.min()))), 2
bins = np.arange(LO, HI + 0.5, 0.5)
fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=200, facecolor=SURF); ax.set_facecolor(SURF)
ax.hist([np.log10(g.mean_com) for _, _, g in groups], bins=bins, stacked=True, color=[c for _, c, _ in groups],
        edgecolor=SURF, linewidth=1.5, label=[f"{n} (n = {len(g)})" for n, _, g in groups])
ax.set_ylim(0, ax.get_ylim()[1] * 1.3)
for q, lab in ((0.5, "median"), (0.9, "90th percentile")):
    x = np.log10(v.quantile(q)); left = q == 0.5
    ax.axvline(x, color=INK2, lw=1, ls=(0, (4, 3)))
    ax.text(x - 0.1 if left else x + 0.1, ax.get_ylim()[1] * 0.97, f"{lab}*\n{v.quantile(q):.2g}%", va="top",
            ha="right" if left else "left", fontsize=8, color=INK2)
step = 2 if HI - LO > 8 else 1
ax.set_xlim(LO, HI); ax.set_xticks(range(LO, HI + 1, step)); ax.set_xticklabels([rf"$10^{{{k}}}$" for k in range(LO, HI + 1, step)])
ax.set_xlabel("Mean combined relative error over the simulated time (%)"); ax.set_ylabel("Number of runs")
ax.grid(True, axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.spines[["top", "right"]].set_visible(False)
ax.legend(loc="upper left", fontsize=8, frameon=False)
c = df.stop_group.value_counts()
n_early = int((classify(pd.read_csv(f"{S}/figures/lhs_b0p0min1_stop_alpha0.1_n2000_seed{seed}.csv")).stop_group == "stopped before 2nd time point").sum()) if N == "1000balanced2" else 0
nz = df[df.zero_error & df.stop_group.isin(["completed full window", "stopped: below 10⁻³/ml"])]
fig.suptitle(f"Robustness of the Boolean model at α = 0.1: {TITLE_N} (seed {seed})" if not SUB_N else f"Robustness of the Boolean model at α = 0.1: {TITLE_N}", x=0.02, ha="left", fontsize=11)
fig.text(0.02, 0.905, SUB_N + ("runs stop" if SUB_N else "Simulations stop") + " when a population in either model leaves [10⁻³, 10¹⁵] per ml", fontsize=9, color=INK2)
fig.text(0.02, 0.02,
         "Log-uniform sampling of r, a, c, m, B₀, P₀ and τ over the Table 1 ranges, with B₀ and P₀ from 1 per ml. Both models are compared up to the\n"
         "earlier stop time. " + (f"Stopped runs drawn only from those past the 1st time point ({n_early} of the 2,000 collapsed sooner). Not plotted: "
         if N == "1000balanced2" else f"Not plotted: {c.get('stopped before 2nd time point', 0)} runs stopped before the 2nd of 10,000 time points (nothing to compare), ") +
         f"{len(nz)} runs with error exactly 0,\n{c.get('stopped: above 10¹⁵/ml', 0)} runs stopped above 10¹⁵/ml, "
         f"{c.get('solver failed', 0)} solver failures. *Median and 90th percentile of the plotted runs.",
         fontsize=7, color=MUTED)
fig.subplots_adjust(left=0.09, right=0.97, top=0.88, bottom=0.25)
fig.savefig(f"{S}/figures/lhs_b0p0min1_stop_alpha0.1_n{N}_error_distribution_seed{seed}.png", facecolor=SURF)
print(seed, {n: len(g) for n, _, g in groups}, "| median %.3g%%, 90th %.3g%%" % (v.median(), v.quantile(.9)))
