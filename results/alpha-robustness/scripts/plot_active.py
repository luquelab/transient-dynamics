import pandas as pd, numpy as np, matplotlib, os
matplotlib.use("Agg")
import matplotlib.pyplot as plt
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INK, INK2, MUTED, GRID, SURF, BLUE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb", "#2a78d6"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK2, "text.color": INK})
PROCS = ["growth", "predation", "burst", "decay"]
alphas = ["3.7e-14", "0.0001", "0.001", "0.01", "0.1", "1", "10", "12", "13"]
# load each run; consecutive alphas with an identical switching pattern share one panel
groups = []
for a in alphas:
    d = pd.read_excel(f"{S}/alpha_sweep/alpha_{a}/generated-data/standard_model/results_standard/array_outputs.xlsx")
    act = np.column_stack([d[f"weight_{p}"].to_numpy() >= 1 for p in PROCS])  # same test as the model's theta switches
    if groups and np.array_equal(groups[-1][2], act):
        groups[-1][0].append(a)
    else:
        groups.append([[a], d.time.to_numpy(), act])
fig, axes = plt.subplots(len(groups), 1, sharex=True, figsize=(7.5, 1.6 * len(groups) + 1), dpi=200, facecolor=SURF)
for ax, (names_a, t, act) in zip(axes, groups):
    def sci(x):
        m, e = f"{float(x):.1e}".split("e")
        return rf"${m}\times10^{{{int(e)}}}$"
    a = ", ".join(sci(x) if "e" in x else x for x in names_a)
    if names_a == ["3.7e-14"]: a += "  (α*: all processes relevant at all times)"
    n = act.sum(axis=1)
    ax.set_facecolor(SURF)
    ax.fill_between(t, 0, n, step="post", color=BLUE, alpha=0.15, lw=0)
    ax.step(t, n, where="post", color=BLUE, lw=2)
    ch = np.r_[0, np.where((act[1:] != act[:-1]).any(axis=1))[0] + 1, len(t)]
    for s, e in zip(ch[:-1], ch[1:]):
        names = ", ".join(p.capitalize() for p, on in zip(PROCS, act[s]) if on) or "None (populations frozen at initial values)"
        mid = (t[s] + t[e - 1]) / 2
        if t[e - 1] - t[s] < 3.3 and names.count(",") >= 2:  # short segment: wrap onto two lines
            parts = names.split(", "); names = ", ".join(parts[:2]) + ",\n" + ", ".join(parts[2:])
        ax.text(mid, n[s] + 0.25, names, ha="center", va="bottom", fontsize=7.5, color=INK2)
        if s > 0:
            close_next = e < len(t) and t[e] - t[s] < 1.0  # next switch too close: label on the left side of the line
            ax.text(t[s] - 0.1 if close_next else t[s] + 0.1, 0.3, f"{t[s]:.2f} h", ha="right" if close_next else "left",
                    va="bottom", fontsize=7.5, color=INK2)
    ax.set_ylim(0, 5.2); ax.set_yticks(range(5)); ax.set_xlim(0, 14)
    ax.grid(True, axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s_ in ["top", "right"]: ax.spines[s_].set_visible(False)
    ax.set_title(f"α = {a}", loc="left", fontsize=10, color=INK, pad=2)
axes[-1].set_xlabel("Time (h)")
fig.supylabel("Active processes (weight ≥ 1)", color=INK2, fontsize=10)
fig.suptitle("Scenario 1 (E. coli–T4): active processes in the Boolean model", x=0.02, ha="left", fontsize=11, color=INK)
fig.tight_layout(h_pad=1.6, rect=(0, 0, 1, 0.985))
fig.savefig(f"{S}/figures/scenario1_active_processes_by_alpha.png", facecolor=SURF)
