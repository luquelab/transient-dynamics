"""Summaries for one model's sweep: per alpha, the completed set (exactly TARGET runs drawn at random from all
completed runs) and the stopped set (all stopped runs with at least two compared time points).
usage: analyze.py <standard|cc> <tables/lhs_*.csv.gz>  ->  tables/<model>_*.csv"""
import os, sys, numpy as np, pandas as pd
from scipy.stats import binomtest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARGET, SEED = 1000, 2026

def sets(df):
    """dict alpha -> (completed selection, stopped comparable runs)"""
    out = {}
    for a, g in df.groupby("alpha"):
        done = g[g.status == "completed"]
        sel = done.sample(TARGET, random_state=SEED) if len(done) >= TARGET else done
        stop = g[(g.status == "stopped") & g.comparable.astype(bool)]
        out[a] = tuple(x.assign(far_out=flag_outliers(x)) for x in (sel, stop))
    return out

def flag_outliers(s):
    """Tukey far-out test on log10(mean error): flag runs above Q3 + 3*IQR. Exact zeros are never outliers."""
    pos = s.mean_err[s.mean_err > 0]
    if len(pos) < 4: return pd.Series(False, index=s.index)
    q1, q3 = np.log10(pos).quantile([.25, .75])
    return s.mean_err > 10 ** (q3 + 3 * (q3 - q1))

def ci(k, n):
    if n == 0: return (np.nan, np.nan)
    c = binomtest(int(k), int(n)).proportion_ci(method="wilson"); return (c.low, c.high)

def summarize(model, df):
    P = len(models.PROCESSES[model]); fr = [f"frac_{j}_processes" for j in range(P + 1)]
    rows, counts = [], []
    for a, (done, stop) in sets(df).items():
        g = df[df.alpha == a]; n_all = len(g)
        counts.append(dict(alpha=a, samples=n_all, completed=int((g.status == "completed").sum()),
                           stopped=int((g.status == "stopped").sum()),
                           stopped_before_2nd_point=int(((g.status == "stopped") & ~g.comparable.astype(bool)).sum()),
                           solver_failed=int((g.status == "solver failed").sum()),
                           **{f"stop: {r}": int((g.stop_reason == r).sum()) for r in sorted(g.stop_reason.dropna().unique()) if r}))
        for name, s_all in (("completed", done), ("stopped", stop)):
            s = s_all[~s_all.far_out]                   # statistics exclude the flagged extreme outliers
            e = s.mean_err; allon = s[fr[-1]]
            k_sus = int(s.all_on_sustained.sum()); lo, hi = ci(k_sus, len(s))
            # whole-window cases in the sampled space: completed runs (all of them, not only the selection) with every
            # process relevant at every time point, out of all samples for this alpha
            comp = g[g.status == "completed"]
            k_all = int(comp.all_on_sustained.astype(bool).sum()); lo2, hi2 = ci(k_all, n_all)
            rows.append(dict(alpha=a, set=name, n=len(s), n_far_out_excluded=int(s_all.far_out.sum()),
                             far_out_errors_pct="; ".join(f"{v:.3g}" for v in sorted(s_all.mean_err[s_all.far_out])),
                             err_mean_incl_outliers=s_all.mean_err.mean(), err_mean=e.mean(), err_sd=e.std(), err_var=e.var(), err_median=e.median(),
                             err_q25=e.quantile(.25), err_q75=e.quantile(.75), err_p90=e.quantile(.9), err_max=e.max(),
                             n_err_zero=int((e == 0).sum()), n_err_above_10pct=int((e > 10).sum()),
                             frac_err_within_alpha=float((e <= 100 * a).mean()),
                             time_all_relevant_mean=allon.mean(), time_all_relevant_sd=allon.std(), time_all_relevant_var=allon.var(),
                             all_relevant_whole_run=k_sus, all_relevant_whole_run_pct=100 * k_sus / max(len(s), 1),
                             all_relevant_whole_run_ci95=f"{100*lo:.1f}-{100*hi:.1f}%",
                             all_relevant_whole_window_completed_runs=k_all, completed_runs=len(comp),
                             all_relevant_whole_window_pct_of_completed=100 * k_all / max(len(comp), 1),
                             all_relevant_whole_window_pct_of_all_samples=100 * k_all / n_all,
                             all_relevant_whole_window_pct_of_all_samples_ci95=f"{100*lo2:.2f}-{100*hi2:.2f}%",
                             **{f"{c}_mean": s[c].mean() for c in fr}))
    return pd.DataFrame(counts), pd.DataFrame(rows)

if __name__ == "__main__":
    model, path = sys.argv[1], sys.argv[2]
    df = pd.read_csv(path)
    counts, summ = summarize(model, df)
    counts.to_csv(os.path.join(HERE, "tables", f"{model}_outcome_counts.csv"), index=False)
    summ.to_csv(os.path.join(HERE, "tables", f"{model}_summary_by_alpha.csv"), index=False)
    sel = pd.concat([s.assign(set=name) for a, (d, st) in sets(df).items() for name, s in (("completed", d), ("stopped", st))])
    sel.to_csv(os.path.join(HERE, "tables", f"{model}_analysis_sets.csv.gz"), index=False)
    pd.set_option("display.width", 250)
    print(counts.to_string(index=False))
    print(summ[["alpha", "set", "n", "n_far_out_excluded", "err_mean", "err_sd", "err_median", "err_p90", "err_max", "n_err_above_10pct",
                "frac_err_within_alpha", "time_all_relevant_mean", "all_relevant_whole_run", "all_relevant_whole_run_pct",
                "all_relevant_whole_run_ci95", "all_relevant_whole_window_completed_runs", "completed_runs",
                "all_relevant_whole_window_pct_of_all_samples", "all_relevant_whole_window_pct_of_all_samples_ci95"]].to_string(index=False, float_format=lambda x: f"{x:.3g}"))
