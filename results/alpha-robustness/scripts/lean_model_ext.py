"""lean_model with an extinction threshold: a population that falls to FLOOR is set to 0 and stays at 0
(absorbing), in both the full and the Boolean model. Same equations, solver, tolerances and grid as the notebook."""
import numpy as np
from scipy.integrate import solve_ivp
from lean_model import theta_approach_model, lotka_volterra_standard
FLOOR = 1e-3

def _solve(f, y0, tobs, para, t_eval):
    """Integrate f with absorbing extinction at FLOOR; returns y on t_eval, extinction times, success, message."""
    y = np.full((2, len(t_eval)), np.nan); y0 = np.array(y0, float); t0 = 0.0
    dead = y0 < FLOOR; y0[dead] = 0.0; t_ext = [0.0 if d else np.nan for d in dead]
    while True:
        evs = []
        for k in range(2):
            ev = (lambda k: (lambda t, z, *a: z[k] - FLOOR))(k)
            ev.terminal, ev.direction = (not dead[k]), -1
            evs.append(ev)
        idx = np.where(t_eval >= t0)[0]
        if len(idx) == 0: break
        rhs = lambda t, z, *a: _masked(f, t, z, a, dead)  # extinct populations stay at 0
        with np.errstate(all="ignore"):
            sol = solve_ivp(rhs, (t0, tobs), y0, args=para, method="LSODA", t_eval=t_eval[idx],
                            events=evs, rtol=1e-13, atol=1e-16)
        if not sol.success: return y, t_ext, False, sol.message
        y[:, idx[:sol.y.shape[1]]] = sol.y
        hit = [k for k in range(2) if not dead[k] and len(sol.t_events[k])]
        if not hit or sol.status != 1: break
        te = min(sol.t_events[k][0] for k in hit); k = min(hit, key=lambda k: sol.t_events[k][0])
        dead[k] = True; t_ext[k] = te
        y0 = sol.y_events[k][0].copy(); y0[k] = 0.0; t0 = te
        y[k, t_eval >= te] = 0.0
    return y, t_ext, True, ""

def _masked(f, t, z, a, dead):
    d = f(t, z, *a)
    return [0.0 if dead[0] else d[0], 0.0 if dead[1] else d[1]]

def simulate(r, a, c, m, tobs, b0, p0, alpha=0.1, threshold=1.0, n=10000):
    para = (r, a, c, m, tobs, alpha, threshold)
    t_eval = np.linspace(0, tobs, n)
    ya, ea, oka, ma = _solve(theta_approach_model, [b0, p0], tobs, para, t_eval)
    yf, ef, okf, mf = _solve(lotka_volterra_standard, [b0, p0], tobs, para, t_eval)
    out = dict(ok=bool(oka and okf and np.isfinite(ya).all() and np.isfinite(yf).all()), msg=ma or mf,
               ext_B_full=ef[0], ext_P_full=ef[1], ext_B_bool=ea[0], ext_P_bool=ea[1])
    if not (oka and okf): return out
    b, p = ya; bf, pf = yf
    with np.errstate(all="ignore"):
        def rel(full, approx):  # % relative error; both extinct = agree (0), only full extinct = 100
            e = np.abs(full - approx) / full * 100
            return np.where(full == 0, np.where(approx == 0, 0.0, 100.0), e)
        bac, pha = rel(bf, b), rel(pf, p); com = (bac + pha) / 2
        w = np.vstack([tobs * r * b / (alpha * b), tobs * a * p * b / (alpha * b),
                       tobs * c * a * p * b / (alpha * p), tobs * m * p / (alpha * p)])
    act = np.nan_to_num(w, nan=0.0) >= threshold  # processes of an extinct population count as off
    t_first = np.nanmin(ef + ea) if np.isfinite(ef + ea).any() else np.inf  # first extinction in either model
    pre = t_eval < t_first
    out.update(mean_com_pre=float(np.mean(com[pre])) if pre.any() else np.nan, frac_pre=float(pre.mean()),
               mean_com=float(np.mean(com)), median_com=float(np.median(com)), max_com=float(np.max(com)),
               mean_bac=float(np.mean(bac)), mean_pha=float(np.mean(pha)),
               all_on=bool(act.all()), none_ever=bool(not act.any()), avg_active=float(act.sum(0).mean()),
               bf_max=float(np.nanmax([bf.max(), b.max()])), pf_max=float(np.nanmax([pf.max(), p.max()])),
               bf_min=float(min(bf.min(), b.min())), pf_min=float(min(pf.min(), p.min())),
               finite=bool(np.isfinite(com).all()))
    return out
