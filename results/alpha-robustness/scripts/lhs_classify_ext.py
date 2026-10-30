"""Outcome classes for the extinction-threshold sweep (populations below 1e-3/ml are extinct, not excluded).
Applied in order: solver failed; above ceiling (either model, any time); zero error; no process observable; valid."""
import pandas as pd, numpy as np, sys
CEIL = 1e13
def classify(df):
    hi = np.fmax(df.bf_max, df.pf_max)
    cls = pd.Series("valid", index=df.index)
    cls[~df.ok] = "solver failed"
    cls[(cls == "valid") & ((hi > CEIL) | ~np.isfinite(hi))] = "excluded: above 1e13/ml"
    cls[(cls == "valid") & (df.mean_com == 0)] = "zero error (all processes observable)"
    cls[(cls == "valid") & df.none_ever.astype(bool)] = "no process observable"
    ext_f = df.ext_B_full.notna() | df.ext_P_full.notna(); ext_b = df.ext_B_bool.notna() | df.ext_P_bool.notna()
    same = (df.ext_B_full.notna() == df.ext_B_bool.notna()) & (df.ext_P_full.notna() == df.ext_P_bool.notna())
    return df.assign(outcome=cls, extinction_any=ext_f | ext_b, extinction_same_species=same)
if __name__ == "__main__":
    df = classify(pd.read_csv(sys.argv[1]))
    print(df.outcome.value_counts().to_string())
    p = df[df.outcome.isin(["valid", "no process observable"])]
    print(f"\nplotted: {len(p)} | with an extinction in either model: {int(p.extinction_any.sum())} "
          f"| full and Boolean agree on which species go extinct: {int((p.extinction_any & p.extinction_same_species).sum())}/{int(p.extinction_any.sum())}")
    both = p[p.ext_B_full.notna() & p.ext_B_bool.notna()]
    if len(both):
        d = (both.ext_B_bool - both.ext_B_full).abs() / both.tobs * 100
        print(f"bacterial extinction in both models (n={len(both)}): timing difference as % of window, median {d.median():.2g}, 90th {d.quantile(.9):.2g}, max {d.max():.2g}")
    bothp = p[p.ext_P_full.notna() & p.ext_P_bool.notna()]
    if len(bothp):
        d = (bothp.ext_P_bool - bothp.ext_P_full).abs() / bothp.tobs * 100
        print(f"phage extinction in both models (n={len(bothp)}): timing difference as % of window, median {d.median():.2g}, 90th {d.quantile(.9):.2g}, max {d.max():.2g}")
    for name, g in [("all plotted", p), ("with extinction", p[p.extinction_any]), ("without extinction", p[~p.extinction_any])]:
        print(f"{name}: n={len(g)} median {g.mean_com.median():.3g}% 90th {g.mean_com.quantile(.9):.3g}% max {g.mean_com.max():.3g}% | >10%: {int((g.mean_com > 10).sum())}")
