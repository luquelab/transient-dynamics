"""Mutually exclusive outcome classes for an LHS result file (applied in this order).
Bounds apply to both models at every time point: below FLOOR (1 per microliter) or above CEIL is excluded."""
import pandas as pd, numpy as np, sys
FLOOR, CEIL = 1e-3, 1e13
def classify(df):
    lo = np.fmin(df.bf_min, df.pf_min); hi = np.fmax(df.bf_max, df.pf_max)  # NaN-safe
    below, above = lo < FLOOR, (hi > CEIL) | ~np.isfinite(hi)
    cls = pd.Series("valid", index=df.index)
    cls[~df.ok] = "solver failed"
    cls[(cls == "valid") & below & above] = "excluded: below floor and above ceiling"
    cls[(cls == "valid") & below] = "excluded: below 1e-3/ml"
    cls[(cls == "valid") & above] = "excluded: above 1e13/ml"
    cls[(cls == "valid") & ~df.finite.astype(bool)] = "non-finite error within bounds"
    cls[(cls == "valid") & (df.mean_com == 0)] = "zero error (all processes observable)"
    cls[(cls == "valid") & df.none_ever.astype(bool)] = "no process observable"
    return df.assign(outcome=cls, start_below_floor=(df.b0 < FLOOR) | (df.p0 < FLOOR))
if __name__ == "__main__":
    df = classify(pd.read_csv(sys.argv[1]))
    print(df.outcome.value_counts().to_string())
    ex = df.outcome.str.startswith("excluded")
    print(f"\nexcluded total {int(ex.sum())}; of those starting below the floor (B0 or P0 < 1e-3): {int((ex & df.start_below_floor).sum())}")
    for name in ["valid", "no process observable"]:
        v = df[df.outcome == name].mean_com
        print(f"{name}: n={len(v)}  quantiles 10/50/90/max: " + ", ".join(f"{x:.3g}" for x in v.quantile([.1, .5, .9, 1])) + f"  | >10%: {int((v > 10).sum())}")
