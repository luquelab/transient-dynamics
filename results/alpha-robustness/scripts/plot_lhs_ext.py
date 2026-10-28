import pandas as pd, numpy as np, matplotlib, os, sys
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lhs_classify_ext import classify
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 2026
df = classify(pd.read_csv(f"{S}/figures/lhs_b0p0min1_ext_alpha0.1_n1000_seed{seed}.csv"))
df.to_csv(f"{S}/figures/lhs_b0p0min1_ext_alpha0.1_n1000_seed{seed}_classified.csv", index=False)
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK2, "text.color": INK})
p = df[df.outcome.isin(["valid", "no process observable"])]
groups = [("No extinction, some processes observable", BLUE, p[(p.outcome == "valid") & ~p.extinction_any]),
          ("No process observable", ORANGE, p[p.outcome == "no process observable"]),
          ("Extinction in either model (<10⁻³/ml)", AQUA, p[(p.outcome == "valid") & p.extinction_any])]
allv = p.mean_com
LO, HI = int(np.floor(np.log10(allv.min()))), 2
bins = np.arange(LO, HI + 0.5, 0.5)
fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=200, facecolor=SURF); ax.set_facecolor(SURF)
ax.hist([np.log10(g.mean_com) for _, _, g in groups], bins=bins, stacked=True, color=[c for _, c, _ in groups],
        edgecolor=SURF, linewidth=1.5, label=[f"{n} (n = {len(g)})" for n, _, g in groups])
ax.set_ylim(0, ax.get_ylim()[1] * 1.3)
for q, lab in ((0.5, "median"), (0.9, "90th percentile")):
    x = np.log10(allv.quantile(q)); left = q == 0.5
    ax.axvline(x, color=INK2, lw=1, ls=(0, (4, 3)))
    ax.text(x - 0.1 if left else x + 0.1, ax.get_ylim()[1] * 0.97, f"{lab}*\n{allv.quantile(q):.2g}%", va="top",
            ha="right" if left else "left", fontsize=8, color=INK2)
step = 2 if HI - LO > 8 else 1
ax.set_xlim(LO, HI); ax.set_xticks(range(LO, HI + 1, step)); ax.set_xticklabels([rf"$10^{{{k}}}$" for k in range(LO, HI + 1, step)])
ax.set_xlabel("Mean combined relative error, full vs. Boolean model (%)"); ax.set_ylabel("Number of runs")
ax.grid(True, axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.spines[["top", "right"]].set_visible(False)
ax.legend(loc="upper left", fontsize=8, frameon=False)
c = df.outcome.value_counts()
fig.suptitle(f"Robustness of the Boolean model at α = 0.1: 1,000 Latin hypercube samples (seed {seed})", x=0.02, ha="left", fontsize=11)
fig.text(0.02, 0.905, "Extinction threshold 10⁻³ per ml (absorbing, both models); ceiling 10¹³ per ml", fontsize=9, color=INK2)
fig.text(0.02, 0.02,
         "Log-uniform sampling of r, a, c, m, B₀, P₀ and τ over the Table 1 ranges, with B₀ and P₀ from 1 per ml. After an extinction,\n"
         "both models at 0 count as agreement (0%) and one model at 0 as 100%. Not plotted: "
         f"{c.get('excluded: above 1e13/ml', 0)} runs above 10¹³/ml, {c.get('zero error (all processes observable)', 0)} runs with error exactly 0, "
         f"{c.get('solver failed', 0)} solver failures.\n*Median and 90th percentile of all plotted runs.",
         fontsize=7, color=MUTED)
fig.subplots_adjust(left=0.09, right=0.97, top=0.88, bottom=0.25)
fig.savefig(f"{S}/figures/lhs_b0p0min1_ext_alpha0.1_error_distribution_seed{seed}.png", facecolor=SURF)
print(seed, "plotted", len(p), {n: len(g) for n, _, g in groups})
