"""Tally the stop-at-threshold sweep: why each run ended, and the error over the simulated time."""
import pandas as pd, numpy as np, sys
def classify(df):
    reason = df.stop.fillna("solver failed")
    group = np.where(reason == "completed", "completed full window",
            np.where(reason.str.contains("below floor"), "stopped: below 10⁻³/ml",
            np.where(reason.str.contains("above ceiling"), "stopped: above 10¹⁵/ml", "solver failed")))
    group = np.where((group == "stopped: below 10⁻³/ml") & (df.n_points < 2), "stopped before 2nd time point", group)
    return df.assign(stop=reason, stop_group=group, zero_error=df.mean_com == 0)
if __name__ == "__main__":
    df = classify(pd.read_csv(sys.argv[1]))
    print(df.stop_group.value_counts().to_string()); print(); print(df.stop.value_counts().to_string())
    print(f"\nexactly zero error: {int(df.zero_error.sum())} | no process observable: {int(df.none_ever.sum())} | non-finite error: {int((~df.finite.astype(bool)).sum())}")
    for g, d in df.groupby("stop_group"):
        e = d.mean_com[d.mean_com > 0]
        print(f"{g:28s} n={len(d):4d} | time simulated: median {d.frac_run.median()*100:.2g}% of window | "
              f"error median {e.median():.3g}% 90th {e.quantile(.9):.3g}% max {e.max():.3g}% | >10%: {int((e > 10).sum())}")
