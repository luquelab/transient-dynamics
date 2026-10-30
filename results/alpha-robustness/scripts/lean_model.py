import os
"""Stripped-down notebook pipeline: same models (notebook cell 14), solver, tolerances, time grid,
weight test and error definitions as transient_dynamics.ipynb, without figures or Excel output."""
import json, numpy as np
from scipy.integrate import solve_ivp
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))  # repository root
_cells = json.load(open(ROOT + "/notebooks/transient_dynamics.ipynb"))["cells"]
exec("".join(_cells[14]["source"]))  # defines theta_approach_model, lotka_volterra_standard

def simulate(r, a, c, m, tobs, b0, p0, alpha=0.1, threshold=1.0, n=10000):
    para = (r, a, c, m, tobs, alpha, threshold)
    t_eval = np.linspace(0, tobs, n)
    kw = dict(args=para, method="LSODA", t_eval=t_eval, rtol=1e-13, atol=1e-16)
    ad = solve_ivp(theta_approach_model, (0, tobs), [b0, p0], **kw)
    fu = solve_ivp(lotka_volterra_standard, (0, tobs), [b0, p0], **kw)
    out = dict(ok=bool(ad.success and fu.success and ad.y.shape[1] == n and fu.y.shape[1] == n),
               msg=(ad.message if not ad.success else fu.message if not fu.success else ""))
    if not out["ok"]:
        return out
    b, p, bf, pf = ad.y[0], ad.y[1], fu.y[0], fu.y[1]
    with np.errstate(all="ignore"):
        # Boolean-model weights (notebook cell 26), weight >= threshold means the process is on
        w = np.vstack([tobs * r * b / (alpha * b), tobs * a * p * b / (alpha * b),
                       tobs * c * a * p * b / (alpha * p), tobs * m * p / (alpha * p)])
        bac = np.abs(bf - b) / bf * 100; pha = np.abs(pf - p) / pf * 100
        com = (bac + pha) / 2
    act = w >= threshold
    out.update(mean_com=float(np.mean(com)), median_com=float(np.median(com)), max_com=float(np.max(com)),
               mean_bac=float(np.mean(bac)), mean_pha=float(np.mean(pha)),
               all_on=bool(act.all()), none_ever=bool(not act.any()), avg_active=float(act.sum(0).mean()),
               bf_max=float(max(bf.max(), b.max())), pf_max=float(max(pf.max(), p.max())),
               bf_min=float(min(bf.min(), b.min())), pf_min=float(min(pf.min(), p.min())),
               finite=bool(np.isfinite(com).all()))
    return out
