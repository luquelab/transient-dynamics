"""lean_model with stopping rules: integration stops when either population, in either model,
falls below FLOOR or rises above CEIL. Both models are compared up to the earlier stop time.
Same equations, solver, tolerances and time grid as the notebook."""
import numpy as np
from scipy.integrate import solve_ivp
from lean_model import theta_approach_model, lotka_volterra_standard
FLOOR, CEIL = 1e-3, 1e15
NAMES = ["B below floor", "P below floor", "B above ceiling", "P above ceiling"]

def _events():
    evs = []
    for k, (lim, d) in enumerate([(FLOOR, -1), (FLOOR, -1), (CEIL, 1), (CEIL, 1)]):
        ev = (lambda i, lim: (lambda t, z, *a: z[i] - lim))(k % 2, lim)
        ev.terminal, ev.direction = True, d
        evs.append(ev)
    return evs

def _solve(f, y0, tobs, para, t_eval):
    with np.errstate(all="ignore"):
        sol = solve_ivp(f, (0, tobs), y0, args=para, method="LSODA", t_eval=t_eval, events=_events(), rtol=1e-13, atol=1e-16)
    hit = [(sol.t_events[k][0], k) for k in range(4) if len(sol.t_events[k])]
    t_stop, why = min(hit) if hit else (tobs, None)
    return sol, t_stop, why

def simulate(r, a, c, m, tobs, b0, p0, alpha=0.1, threshold=1.0, n=10000):
    para = (r, a, c, m, tobs, alpha, threshold)
    t_eval = np.linspace(0, tobs, n)
    sa, ta, wa = _solve(theta_approach_model, [b0, p0], tobs, para, t_eval)
    sf, tf, wf = _solve(lotka_volterra_standard, [b0, p0], tobs, para, t_eval)
    out = dict(ok=bool(sa.status in (0, 1) and sf.status in (0, 1)), msg=sa.message if sa.status < 0 else sf.message if sf.status < 0 else "")
    if not out["ok"]: return out
    t_stop = min(ta, tf)
    stop = "completed" if t_stop >= tobs else \
        (("full: " if tf <= ta else "Boolean: ") + NAMES[wf if tf <= ta else wa])
    k = int(np.sum(t_eval <= t_stop)); k = min(k, sa.y.shape[1], sf.y.shape[1])  # time points both models reached
    b, p = sa.y[0, :k], sa.y[1, :k]; bf, pf = sf.y[0, :k], sf.y[1, :k]
    with np.errstate(all="ignore"):
        bac = np.abs(bf - b) / bf * 100; pha = np.abs(pf - p) / pf * 100; com = (bac + pha) / 2
        w = np.vstack([tobs * r * b / (alpha * b), tobs * a * p * b / (alpha * b),
                       tobs * c * a * p * b / (alpha * p), tobs * m * p / (alpha * p)])
    act = w >= threshold
    out.update(stop=stop, t_stop=float(t_stop), frac_run=float(t_stop / tobs), n_points=k,
               mean_com=float(np.mean(com)), median_com=float(np.median(com)), max_com=float(np.max(com)),
               mean_bac=float(np.mean(bac)), mean_pha=float(np.mean(pha)),
               all_on=bool(act.all()), none_ever=bool(not act.any()), avg_active=float(act.sum(0).mean()),
               finite=bool(np.isfinite(com).all()))
    return out
