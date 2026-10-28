"""Are floor exclusions real model behaviour or solver artefacts?
Re-solve the full Lotka-Volterra model in log variables (x = ln B, y = ln P):
  dx/dt = r - a e^y,  dy/dt = c a e^x - m
Positivity is exact by construction, so any crossing below the floor is the model's, not the solver's."""
import pandas as pd, numpy as np, os, sys
from scipy.integrate import solve_ivp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lhs_classify import classify, FLOOR
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
df = classify(pd.read_csv(f"{S}/figures/lhs_b0p0min1_alpha0.1_n1000_seed2026.csv"))
ex = df[df.outcome == "excluded: below 1e-3/ml"]
neg = (ex.bf_min < 0) | (ex.pf_min < 0)
print(f"floor-excluded runs: {len(ex)}; with a negative value in the original (linear) solve: {int(neg.sum())}")
rows = []
for _, s in ex.iterrows():
    f = lambda t, z: [s.r - s.a * np.exp(z[1]), s.c * s.a * np.exp(z[0]) - s.m]
    t = np.linspace(0, s.tobs, 10000)
    sol = solve_ivp(f, (0, s.tobs), [np.log(s.b0), np.log(s.p0)], method="LSODA", t_eval=t, rtol=1e-10, atol=1e-12)
    if not sol.success or sol.y.shape[1] < len(t):
        rows.append(dict(ok=False)); continue
    lb, lp = sol.y; lf = np.log(FLOOR)
    below = (lb < lf) | (lp < lf)
    first = t[np.argmax(below)] / s.tobs if below.any() else np.nan
    rows.append(dict(ok=True, crosses=bool(below.any()), first_cross_frac=first, frac_time_below=below.mean(),
                     log10_min=min(lb.min(), lp.min()) / np.log(10)))
r = pd.DataFrame(rows)
ok = r[r.ok]
print(f"log-space re-solve succeeded: {len(ok)}/{len(r)}")
print(f"full model goes below 1e-3/ml in the exact-positive solve: {int(ok.crosses.sum())}/{len(ok)}")
c = ok[ok.crosses]
print("first crossing, as a fraction of the window: quartiles", c.first_cross_frac.quantile([.25, .5, .75]).round(2).tolist())
print("fraction of the window spent below the floor: quartiles", c.frac_time_below.quantile([.25, .5, .75]).round(2).tolist())
print("lowest population reached (log10 per ml): quartiles", c.log10_min.quantile([.25, .5, .75]).round(1).tolist(), " min", round(c.log10_min.min(), 1))
