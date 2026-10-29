"""Test of the 'pre-relevant process' hypothesis at alpha = 0.1, completed runs.
For every run: time each process spends with weight in the pre-relevant zone [LOW, 1) in the transient (Boolean) model,
its longest continuous stay there, process switching, and the neglected-effect index sum_p mean_t(w_p * 1[w_p < 1])
(the time-integrated per-capita change of the neglected processes, in units of alpha).
usage: prerelevance.py  ->  tables/prerelevance_alpha0.1.csv"""
import os, sys, numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ALPHA, LOW = 0.1, 0.1

def trajectories(model, p, alpha=ALPHA):
    t = np.linspace(0, p["tobs"], models.N_POINTS)
    if model == "standard":
        ff, fb = models.lotka_volterra_standard, models.theta_approach_model
        args = (p["r"], p["a"], p["c"], p["m"], p["tobs"], alpha, 1.0)
    else:
        ff, fb = models.lotka_volterra_standard_cc, models.theta_approach_model_carrying_capacity
        args = (p["r"], p["a"], p["c"], p["m"], p["K"], p["tobs"], alpha, 1.0)
    sf, tf, _ = models._solve(ff, [p["b0"], p["p0"]], p["tobs"], args, t)
    sb, tb, _ = models._solve(fb, [p["b0"], p["p0"]], p["tobs"], args, t)
    k = min(int(np.sum(t <= min(tf, tb))), sf.y.shape[1], sb.y.shape[1])
    return t[:k], sf.y[:, :k], sb.y[:, :k]

def longest_run(mask):
    best = cur = 0
    for v in mask:
        cur = cur + 1 if v else 0; best = max(best, cur)
    return best

def metrics(arg):
    model, row = arg
    p = {k: row[k] for k in ["r", "a", "c", "m", "tobs", "b0", "p0"] + (["K"] if model == "cc" else [])}
    t, yf, yb = trajectories(model, p)
    w = models.weights(model, p, yb[0], yb[1], ALPHA)          # transient-model weights, as used by its switches
    pre = (w >= LOW) & (w < 1); act = w >= 1; n = w.shape[1]
    frac_pre = pre.mean(axis=1)
    return dict(model=model, sample=int(row["sample"]), mean_err=row["mean_err"], above_tolerance=row["mean_err"] > 100 * ALPHA,
                tobs=p["tobs"], frac_time_any_pre=float(pre.any(axis=0).mean()), frac_time_pre_max_process=float(frac_pre.max()),
                longest_pre_stay=float(max(longest_run(pre[j]) for j in range(len(w))) / n),
                neglected_index=float(np.sum(np.where(w < 1, w, 0).mean(axis=1))),
                n_switches=int(np.sum(np.any(act[:, 1:] != act[:, :-1], axis=0))),
                **{f"frac_pre_{name}": float(f) for name, f in zip(models.PROCESSES[model], frac_pre)})

if __name__ == "__main__":
    jobs = []
    for model in ("standard", "cc"):
        sel = pd.read_csv(os.path.join(HERE, "tables", f"{model}_analysis_sets.csv.gz"))
        sel = sel[(sel.alpha == ALPHA) & (sel.set == "completed")]
        jobs += [(model, r) for _, r in sel.iterrows()]
    with Pool(4) as pool:
        df = pd.DataFrame(pool.map(metrics, jobs, chunksize=4))
    df.to_csv(os.path.join(HERE, "tables", "prerelevance_alpha0.1.csv"), index=False)
    cols = ["tobs", "frac_time_any_pre", "frac_time_pre_max_process", "longest_pre_stay", "neglected_index", "n_switches"]
    pd.set_option("display.width", 220)
    print(df.groupby(["model", "above_tolerance"]).size().to_string())
    print(df.groupby(["model", "above_tolerance"])[cols].median().to_string(float_format=lambda x: f"{x:.3g}"))
