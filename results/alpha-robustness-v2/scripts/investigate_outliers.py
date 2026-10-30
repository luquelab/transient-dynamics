"""Characterise runs whose mean error exceeds the observer tolerance (100*alpha %) and compare them with the rest.
Re-simulates each run and records trajectory features: closest approach to the floor and ceilings, number of
process switches, oscillation, and the share of the error coming from near-floor time points.
usage: investigate_outliers.py  ->  tables/outlier_diagnostics.csv"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def features(model, p, alpha):
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
    bf, pf, b, ph = sf.y[0, :k], sf.y[1, :k], sb.y[0, :k], sb.y[1, :k]
    com = (np.abs(bf - b) / bf + np.abs(pf - ph) / pf) / 2 * 100
    act = models.weights(model, p, b, ph, alpha) >= 1
    lo = np.minimum.reduce([bf, pf, b, ph])
    near = lo < 1.0                                   # any population below 1 per ml at that time point
    lp = np.log10(pf)                                 # full-model phage oscillation: count local minima
    n_troughs = int(np.sum((lp[1:-1] < lp[:-2]) & (lp[1:-1] < lp[2:])))
    return dict(min_pop_log10=float(np.log10(lo.min())), max_B_log10=float(np.log10(max(bf.max(), b.max()))),
                max_P_log10=float(np.log10(max(pf.max(), ph.max()))),
                frac_time_below_1_per_ml=float(near.mean()),
                share_of_error_from_below_1=float(com[near].sum() / com.sum()) if com.sum() > 0 else 0.0,
                mean_err_excluding_below_1=float(com[~near].mean()) if (~near).any() else np.nan,
                n_switches=int(np.sum(np.any(act[:, 1:] != act[:, :-1], axis=0))),
                n_troughs_phage=n_troughs, avg_relevant=float(act.sum(0).mean()))

if __name__ == "__main__":
    rows = []
    for model in ("standard", "cc"):
        sel = pd.read_csv(os.path.join(HERE, "tables", f"{model}_analysis_sets.csv.gz"))
        for (alpha, which), g in sel.groupby(["alpha", "set"]):
            if alpha not in (0.1, 1.0): continue
            above = g[g.mean_err > 100 * alpha]
            ref = g[g.mean_err <= 100 * alpha].sample(min(150, (g.mean_err <= 100 * alpha).sum()), random_state=2026)
            for label, part in (("above tolerance", above), ("within tolerance (sample)", ref)):
                for _, s in part.iterrows():
                    p = {k: s[k] for k in ["r", "a", "c", "m", "tobs", "b0", "p0"] + (["K"] if model == "cc" else [])}
                    rows.append(dict(model=model, alpha=alpha, set=which, group=label, sample=int(s["sample"]),
                                     mean_err=s.mean_err, **p, **features(model, p, alpha)))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(HERE, "tables", "outlier_diagnostics.csv"), index=False)
    cols = ["mean_err", "min_pop_log10", "frac_time_below_1_per_ml", "share_of_error_from_below_1", "mean_err_excluding_below_1",
            "n_switches", "n_troughs_phage", "avg_relevant", "max_B_log10", "max_P_log10", "tobs", "r", "m"]
    pd.set_option("display.width", 250)
    print(df.groupby(["model", "alpha", "set", "group"]).size().to_string())
    print(df.groupby(["model", "alpha", "set", "group"])[cols].median().to_string(float_format=lambda x: f"{x:.3g}"))
