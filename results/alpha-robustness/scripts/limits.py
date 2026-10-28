"""alpha limits per scenario. weight_p = W_p(t)/alpha with W independent of alpha:
   growth tobs*r, predation tobs*a*P, burst tobs*c*a*B, decay tobs*m.
   alpha_all  = min_t,p W_p(t) on the full trajectory (all on => Boolean == full)
   alpha_none = max_p W_p(t=0): above it every weight starts < 1 and the Boolean model never moves."""
import pandas as pd, numpy as np, re, os, json
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))  # repository root
REPO = ROOT + "/data/transient-dynamics_output-current/example_scenarios/standard_model"
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = {}
for scen in ["scenario_1", "scenario_2", "scenario_3_quasi_equilibrium", "scenario_3_transient", "scenario_4_quasi_static", "scenario_4_transient"]:
    txt = open(f"{REPO}/{scen}/input.txt").read()
    v = {k: float(eval(re.search(rf"^{k}\s*=\s*([^#\n]+)", txt, re.M).group(1))) for k in ["r","a","c","m","tobs"]}
    d = pd.read_excel(f"{REPO}/{scen}/results_standard/array_outputs.xlsx")  # Lucas's alpha = 0.1 run
    B, P = d.bacteria_full.to_numpy(), d.phage_full.to_numpy()
    W = {"growth": np.full_like(B, v["tobs"]*v["r"]), "predation": v["tobs"]*v["a"]*P,
         "burst": v["tobs"]*v["c"]*v["a"]*B, "decay": np.full_like(B, v["tobs"]*v["m"])}
    ok = all(np.allclose(W[p]/0.1, d[f"weight_{p}_full"], rtol=1e-6) for p in W)
    amin = {p: (w.min(), d.time[w.argmin()]) for p, w in W.items()}
    w0 = {p: w[0] for p, w in W.items()}
    pa = min(amin, key=lambda p: amin[p][0]); pn = max(w0, key=w0.get)
    out[scen] = dict(alpha_all=amin[pa][0], all_by=pa, all_at=amin[pa][1], alpha_none=w0[pn], none_by=pn, weights_ok=ok,
                     bmin=B.min(), pmin=P.min())
    print(f"{scen:30s} weights match: {ok} | alpha_all = {amin[pa][0]:.4g} ({pa} at t={amin[pa][1]:.1f} h; Bmin {B.min():.3g}, Pmin {P.min():.3g})"
          f" | alpha_none > {w0[pn]:.4g} ({pn} at t=0)  all W0: " + ", ".join(f"{p} {x:.3g}" for p, x in w0.items()))
json.dump(out, open(f"{S}/figures/alpha_limits.json", "w"), indent=1, default=float)
