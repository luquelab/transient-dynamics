import pandas as pd, numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
d = pd.read_excel(f"{S}/alpha_sweep/alpha_0.5/generated-data/standard_model/results_standard/array_outputs.xlsx")
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
BLUE, ORANGE = "#2a78d6", "#eb6834"
STYLE = {"growth": "-", "predation": (0, (6, 3)), "burst": (0, (1.5, 2)), "decay": (0, (6, 2, 1.5, 2))}
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK2, "text.color": INK})

def crossings(t, w):
    s = w > 1
    i = np.where(s[1:] != s[:-1])[0] + 1
    return [(t[k], "above 1" if s[k] else "below 1") for k in i]

for kind, b, p, suf, title in [("full", "bacteria_full", "phage_full", "_full", "Full Lotka–Volterra model"),
                               ("boolean", "bacteria_adaptive", "phage_adaptive", "", "Boolean (approximated) model")]:
    t = d.time.to_numpy()
    fig, ax = plt.subplots(figsize=(7.5, 4.6), dpi=200, facecolor=SURF); ax.set_facecolor(SURF)
    ax.plot(t, d[b], color=BLUE, lw=2, zorder=3)
    ax.plot(t, d[p], color=ORANGE, lw=2, zorder=3)
    ax.set_yscale("log"); ax.set_ylim(1e-5, 1e13); ax.set_xlim(0, 14)
    ax.text(12.3, d[p].iloc[-1] / 12, "Phage", color=INK2, fontsize=9, ha="left")
    ax.text(12.3, d[b].iloc[int(12.3/14*len(t))] * 8, "Bacteria", color=INK2, fontsize=9, ha="left")
    always = []
    events = []
    for proc in ["growth", "predation", "burst", "decay"]:
        w = d[f"weight_{proc}{suf}"].to_numpy()
        cr = crossings(t, w)
        if not cr:
            always.append(f"{proc.capitalize()}: always {'above' if w[0] > 1 else 'below'} 1 (w = {w[0]:.3g})")
        for tc, direction in cr:
            events.append((tc, proc, direction))
    for j, (tc, proc, direction) in enumerate(sorted(events)):
        ax.axvline(tc, color=INK2, lw=1.2, ls=STYLE[proc], zorder=2)
        arrow = "↑" if direction == "above 1" else "↓"
        ax.text(tc + 0.12, 2e12 if j % 2 == 0 else 1.5e11, f"{proc.capitalize()} {arrow} {tc:.2f} h",
                color=INK, fontsize=8.5, ha="left", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc=SURF, ec="none"))
    ax.set_xlabel("Time (h)"); ax.set_ylabel("Population (per ml)")
    ax.set_title(f"Scenario 1, α = 0.5 — {title}", loc="left", fontsize=11)
    ax.grid(True, which="major", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s_ in ["top", "right"]: ax.spines[s_].set_visible(False)
    handles = [Line2D([], [], color=BLUE, lw=2, label="Bacteria"), Line2D([], [], color=ORANGE, lw=2, label="Phage")]
    handles += [Line2D([], [], color=INK2, lw=1.2, ls=STYLE[pr], label=f"{pr.capitalize()} weight crosses 1")
                for pr in sorted({e[1] for e in events}, key=list(STYLE).index)]
    ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=True, facecolor=SURF, edgecolor=GRID)
    ax.text(0, -0.2, "↑ weight rises above 1 (process becomes relevant)   ↓ weight falls below 1 (process becomes irrelevant)\n"
            + "   ·   ".join(always), transform=ax.transAxes, fontsize=7.5, color=MUTED, va="top")
    fig.tight_layout(); fig.savefig(f"{S}/figures/scenario1_alpha0.5_populations_{kind}.png", facecolor=SURF)
    print(kind, sorted(events), always)
