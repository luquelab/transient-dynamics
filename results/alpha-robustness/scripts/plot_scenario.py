"""Error-vs-alpha and active-process figures for one scenario sweep.
usage: plot_scenario.py <scenario> "<title>" <alpha_all> <alpha_none> <alpha> [<alpha> ...]   (alphas as folder strings)"""
import pandas as pd, numpy as np, matplotlib, matplotlib.ticker, os, sys, re
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))  # repository root
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = ROOT + "/data/transient-dynamics_output-current/example_scenarios/standard_model"
scen, title, a_all, a_none, alphas = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:]
alphas = sorted(alphas, key=float)
tobs = float(eval(re.search(r"^tobs\s*=\s*([^#\n]+)", open(f"{REPO}/{scen}/input.txt").read(), re.M).group(1)))
INK, INK2, MUTED, GRID, SURF, BLUE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb", "#2a78d6"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK2, "text.color": INK})
PROCS = ["growth", "predation", "burst", "decay"]
runs = {a: pd.read_excel(f"{S}/sweeps/{scen}/alpha_{a}/generated-data/standard_model/results_standard/array_outputs.xlsx") for a in alphas}

def sci(x, pct=False):
    x = float(x)
    if x == 0: return "0"
    if 1e-2 <= x < 1e3: return f"{x:.3g}" + ("%" if pct else "")
    m, e = f"{x:.1e}".split("e")
    return rf"${m}\times10^{{{int(e)}}}$" + ("%" if pct else "")
def alabel(a):
    return sci(a) if float(a) < 1e-4 else f"{float(a):g}"

# ---------- error vs alpha ----------
err = pd.DataFrame([dict(a=a, alpha=float(a), y=runs[a].rel_error_com_per.mean()) for a in alphas])
YMIN = 10 ** np.floor(np.log10(err.y[err.y > 0].min()) - 1)
YMAX = 10 ** np.ceil(np.log10(err.y.max()) + 1)
far = err[(err.y == 0) & (err.alpha < 1e-5)]          # zero-error alphas far left: broken-axis panel
near = err[~err.index.isin(far.index)]
if len(far):
    fig, (ax0, ax) = plt.subplots(1, 2, sharey=True, figsize=(6.8, 4.3), dpi=200, facecolor=SURF,
                                  gridspec_kw=dict(width_ratios=[1, 6], wspace=0.06))
else:
    fig, ax = plt.subplots(figsize=(6.8, 4.3), dpi=200, facecolor=SURF); ax0 = None
ax.set_facecolor(SURF); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(YMIN, YMAX)
ALL_TXT = "error = 0\nall processes\nobservable"
if ax0 is not None:
    x = far.alpha.iloc[0]
    ax0.set_facecolor(SURF); ax0.set_xscale("log"); ax0.set_xlim(x / 3, x * 3)
    ax0.plot([x], [YMIN], marker="o", ms=8, mfc=SURF, mec=BLUE, mew=2, clip_on=False, zorder=4)
    ax0.set_xticks([x]); ax0.set_xticklabels([sci(x)], fontsize=8); ax0.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax0.annotate("α_all:\n" + ALL_TXT, (x, YMIN), textcoords="offset points", xytext=(0, 14), ha="center", va="bottom", color=INK2, fontsize=8)
    ax0.spines[["top", "right"]].set_visible(False); ax0.grid(True, axis="y", color=GRID, lw=0.6); ax0.set_axisbelow(True)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", which="both", left=False)
    ax0.set_ylabel("Mean combined relative error (%)")
else:
    ax.set_ylabel("Mean combined relative error (%)")
right_label = False
pos = near[near.y > 0]
ax.plot(pos.alpha, pos.y, color=BLUE, lw=2, marker="o", ms=8, mec=SURF, mew=2, zorder=3)
zer = near[near.y == 0]
ax.plot(zer.alpha, [YMIN] * len(zer), ls="none", marker="o", ms=8, mfc=SURF, mec=BLUE, mew=2, clip_on=False, zorder=4)
if len(zer) and len(pos):  # dashed link from the last zero point to the first positive one
    xz = zer.alpha.max(); nxt = pos[pos.alpha > xz]
    if len(nxt): ax.plot([xz, nxt.alpha.iloc[0]], [YMIN, nxt.y.iloc[0]], color=BLUE, lw=1.2, ls=(0, (3, 2)), zorder=2, clip_on=False)
for _, r in near.iterrows():
    if r.a == a_all and r.y == 0:
        ax.annotate("α_all = " + alabel(r.a) + "\n" + ALL_TXT.replace("\n", " ", 1), (r.alpha, YMIN), textcoords="offset points",
                    xytext=(-4, 12), ha="right", va="bottom", color=INK2, fontsize=8)
    elif r.y > 0:
        lab = sci(r.y, pct=True)
        off, ha, va = (0, 10), "center", "bottom"
        nxt = near[near.alpha > r.alpha]
        if len(nxt) and nxt.a.iloc[0] == a_none and nxt.y.iloc[0] > 10 * r.y:  # steep rise into the collapse: label above-left
            off, ha, va = (-4, 8), "right", "bottom"
        if r.a == a_none:
            lab = f"α_none = {alabel(r.a)}: {lab}\nno process observable"
            prev = near[(near.alpha < r.alpha) & (near.y > 0)]
            below = len(prev) and prev.y.iloc[-1] >= r.y  # line arrives from above or level: label below-left
            off, ha, va = ((10, 0), "left", "center") if below else ((-10, 0), "right", "center")
            right_label = bool(below)
        ax.annotate(lab, (r.alpha, r.y), textcoords="offset points", xytext=off, ha=ha, va=va, color=INK2, fontsize=8.5)
dec = [x for x in near.alpha if np.isclose(np.log10(x), round(np.log10(x)))]
ax.set_xticks(dec); ax.set_xticklabels([alabel(a) for a, x in zip(near.a, near.alpha) if x in dec])
ax.set_xticks([x for x in near.alpha if x not in dec], minor=True)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.tick_params(axis="x", which="minor", length=4)
lo, hi = near.alpha.min(), near.alpha.max()
ax.set_xlim(lo / (12 if len(zer) else 2.5), hi * (150 if right_label else 2.5))  # room for the alpha_all / alpha_none labels
ax.set_xlabel("Observer tolerance, α (reference level)")
ax.grid(True, which="major", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"{title}: Boolean model error vs. α", x=0.02, ha="left", fontsize=11, color=INK)
fig.text(0.02, 0.02, f"Mean over 10,000 time points (0–{tobs:g} h) of the average of bacteria and phage relative errors between the full and\n"
         "Boolean models. Open circles: error exactly 0 (drawn on the bottom edge; a log axis cannot show 0).",
         fontsize=7, color=MUTED)
fig.subplots_adjust(left=0.11, right=0.97, top=0.9, bottom=0.24)
if ax0 is not None:
    kw = dict(color=MUTED, lw=1, clip_on=False, transform=fig.transFigure)
    for axx, side in ((ax0, 1), (ax, 0)):
        p = axx.get_position(); xb = p.x0 + side * p.width
        fig.add_artist(plt.Line2D([xb - 0.006, xb + 0.006], [p.y0 - 0.015, p.y0 + 0.015], **kw))
fig.savefig(f"{S}/figures/{scen}_mean_error_vs_alpha.png", facecolor=SURF); plt.close(fig)

# ---------- active processes over time ----------
groups = []
for a in alphas:
    d = runs[a]; act = np.column_stack([d[f"weight_{p}"].to_numpy() >= 1 for p in PROCS])  # same test as the model's theta switches
    if groups and np.array_equal(groups[-1][2], act): groups[-1][0].append(a)
    else: groups.append([[a], d.time.to_numpy(), act])
fig, axes = plt.subplots(len(groups), 1, sharex=True, figsize=(7.5, 1.6 * len(groups) + 1), dpi=200, facecolor=SURF, squeeze=False)
axes = axes[:, 0]
for ax, (names_a, t, act) in zip(axes, groups):
    n = act.sum(axis=1)
    ax.set_facecolor(SURF)
    ax.fill_between(t, 0, n, step="post", color=BLUE, alpha=0.15, lw=0); ax.step(t, n, where="post", color=BLUE, lw=2)
    ch = np.r_[0, np.where((act[1:] != act[:-1]).any(axis=1))[0] + 1, len(t)]
    span = t[-1]; nsw = len(ch) - 2
    dense = nsw > 8  # too many switches to label individually: summarize in the title instead
    last_lab = -np.inf
    for s_, e_ in zip(ch[:-1], ch[1:]):
        names = ", ".join(p.capitalize() for p, on in zip(PROCS, act[s_]) if on) or "None (populations frozen at initial values)"
        seg = t[e_ - 1] - t[s_]
        if not dense and seg >= 0.12 * span:
            if seg < 0.24 * span and names.count(",") >= 2:
                parts = names.split(", "); names = ", ".join(parts[:2]) + ",\n" + ", ".join(parts[2:])
            ax.text((t[s_] + t[e_ - 1]) / 2, n[s_] + 0.25, names, ha="center", va="bottom", fontsize=7.5, color=INK2)
        if s_ > 0 and not dense and t[s_] - last_lab >= 0.07 * span:
            ax.text(t[s_] + 0.007 * span, 0.3, f"{t[s_]:.3g} h", ha="left", va="bottom", fontsize=7.5, color=INK2)
            last_lab = t[s_]
    extra = ""
    if dense:
        used = [p.capitalize() for j, p in enumerate(PROCS) if act[:, j].any() and not act[:, j].all()]
        extra = f" — {nsw} switches ({', '.join(used)} toggling)"
    ax.set_ylim(0, 5.2); ax.set_yticks(range(5)); ax.set_xlim(0, span)
    ax.grid(True, axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.spines[["top", "right"]].set_visible(False)
    tag = "  (α_all: all processes observable at all times)" if a_all in names_a else "  (α_none: no process observable)" if a_none in names_a else ""
    ax.set_title("α = " + ", ".join(alabel(a) for a in names_a) + tag + extra, loc="left", fontsize=10, color=INK, pad=2)
axes[-1].set_xlabel("Time (h)")
fig.supylabel("Active processes (weight ≥ 1)", color=INK2, fontsize=10)
fig.suptitle(f"{title}: active processes in the Boolean model", x=0.02, ha="left", fontsize=11, color=INK)
fig.tight_layout(h_pad=1.6, rect=(0, 0, 1, 0.985))
fig.savefig(f"{S}/figures/{scen}_active_processes_by_alpha.png", facecolor=SURF); plt.close(fig)
print("ok", scen)
