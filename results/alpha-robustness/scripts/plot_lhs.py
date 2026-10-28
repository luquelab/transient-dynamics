import pandas as pd, numpy as np, matplotlib, os, sys
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lhs_classify import classify
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 2026
df = classify(pd.read_csv(f"{S}/figures/lhs_b0p0min1_alpha0.1_n1000_seed{seed}.csv"))
df.to_csv(f"{S}/figures/lhs_b0p0min1_alpha0.1_n1000_seed{seed}_bounded_classified.csv", index=False)
INK, INK2, MUTED, GRID, SURF, BLUE, ORANGE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb", "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK2, "text.color": INK})
LO, HI = -8, 2  # log10 % range; values >= 10^HI are grouped into the last bin
bins = np.arange(LO, HI + 0.5, 0.5)
def lg(s): return np.clip(np.log10(s), LO, HI - 1e-9)
v = df[df.outcome == "valid"].mean_com; nv = df[df.outcome == "no process observable"].mean_com
fig, ax = plt.subplots(figsize=(7.5, 4.6), dpi=200, facecolor=SURF); ax.set_facecolor(SURF)
ax.hist([lg(v), lg(nv)], bins=bins, stacked=True, color=[BLUE, ORANGE], edgecolor=SURF, linewidth=1.5,
        label=[f"Some processes observable (n = {len(v)})", f"No process observable (n = {len(nv)})"])
allp = pd.concat([v, nv])
ax.set_ylim(0, ax.get_ylim()[1] * 1.3)  # headroom for the percentile labels
for q, lab in ((0.5, "median"), (0.9, "90th percentile")):
    x = np.log10(allp.quantile(q))
    ax.axvline(x, color=INK2, lw=1, ls=(0, (4, 3)))
    left = q == 0.5  # median label left of its line, 90th percentile right of its line
    ax.text(x - 0.08 if left else x + 0.08, ax.get_ylim()[1] * 0.97, f"{lab}*\n{allp.quantile(q):.2g}%", va="top",
            ha="right" if left else "left", fontsize=8, color=INK2)
n_hi = int((v >= 10 ** HI).sum() + (nv >= 10 ** HI).sum())
if n_hi:
    ax.annotate(f"{n_hi} runs ≥ $10^{{{HI}}}$%\n(phage → $10^{{-100}}$/ml)", (HI - 0.25, n_hi), textcoords="offset points",
                xytext=(-8, 18), ha="right", fontsize=8, color=INK2, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
ax.set_xlim(LO, HI); ax.set_xticks(range(LO, HI + 1, 1))
ax.set_xticklabels([rf"$10^{{{k}}}$" for k in range(LO, HI + 1, 1)])
ax.set_xlabel("Mean combined relative error, full vs. Boolean model (%)"); ax.set_ylabel("Number of runs")
ax.grid(True, axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.spines[["top", "right"]].set_visible(False)
ax.legend(loc="upper left", fontsize=8, frameon=False)
c = df.outcome.value_counts()
fig.suptitle(f"Robustness of the Boolean model at α = 0.1: 1,000 Latin hypercube samples (seed {seed})", x=0.02, ha="left", fontsize=11)
fig.text(0.02, 0.905, "Populations bounded to [10⁻³, 10¹³] per ml", fontsize=9, color=INK2)
ex = df.outcome.str.startswith("excluded")
fig.text(0.02, 0.02,
         "Log-uniform sampling of r, a, c, m, B₀, P₀ and τ over the Table 1 ranges, with B₀ and P₀ from 1 per ml.\n"
         f"Excluded (not plotted): {int(ex.sum())} runs where a population in either model left [10⁻³, 10¹³] per ml: "
         f"{c.get('excluded: below 1e-3/ml', 0)} below 10⁻³/ml" + (f" ({nsb} of all excluded started below it)" if (nsb := int((ex & df.start_below_floor).sum())) else "") + ",\n"
         f"{c.get('excluded: above 1e13/ml', 0)} above 10¹³/ml, {c.get('excluded: below floor and above ceiling', 0)} both. "
         f"Also not plotted: {c.get('zero error (all processes observable)', 0)} runs with error exactly 0; "
         f"{c.get('solver failed', 0)} solver failures.\n*Median and 90th percentile of all plotted runs.",
         fontsize=7, color=MUTED)
fig.subplots_adjust(left=0.09, right=0.97, top=0.9, bottom=0.29)
fig.savefig(f"{S}/figures/lhs_b0p0min1_alpha0.1_error_distribution_seed{seed}.png", facecolor=SURF)
print(c.to_string())
