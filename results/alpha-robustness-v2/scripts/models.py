"""Full vs. Boolean (adaptive) phage-bacteria models with stopping bounds.

Model functions come from notebooks/transient_dynamics.ipynb (code cell 14: standard Lotka-Volterra;
code cell 60: with carrying capacity), with the notebook's solver settings and 10,000-point time grid.
A run stops when either population, in either model, leaves its bounds; both models are then compared
up to the earlier stop time. The full model does not depend on alpha, so it is solved once per sample."""
import json, os
import numpy as np
from scipy.integrate import solve_ivp

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
_cells = json.load(open(os.path.join(ROOT, "notebooks", "transient_dynamics.ipynb")))["cells"]
exec("".join(_cells[14]["source"]))  # theta_approach_model, lotka_volterra_standard
exec("".join(_cells[60]["source"]))  # theta_approach_model_carrying_capacity, lotka_volterra_standard_cc

FLOOR = 1e-3                     # per ml, both agents
CEIL = {"B": 1e15, "P": 1e18}    # per ml: smallest bacterium (~1e-3 um^3) and phage (~1e-6 um^3) packed in 1 ml
N_POINTS = 10000
RTOL, ATOL = 1e-13, 1e-16
PROCESSES = {"standard": ["growth", "predation", "burst", "decay"],
             "cc": ["growth", "predation", "burst", "decay", "carrying capacity"]}
_REASONS = ["B below floor", "P below floor", "B above ceiling", "P above ceiling"]

def _events():
    evs = []
    for i, lim, d in [(0, FLOOR, -1), (1, FLOOR, -1), (0, CEIL["B"], 1), (1, CEIL["P"], 1)]:
        ev = (lambda i, lim: (lambda t, z, *a: z[i] - lim))(i, lim)
        ev.terminal, ev.direction = True, d
        evs.append(ev)
    return evs

def _solve(f, y0, tobs, args, t_eval):
    with np.errstate(all="ignore"):
        sol = solve_ivp(f, (0, tobs), y0, args=args, method="LSODA", t_eval=t_eval, events=_events(), rtol=RTOL, atol=ATOL)
    hits = [(sol.t_events[k][0], k) for k in range(4) if len(sol.t_events[k])]
    t_stop, why = min(hits) if hits else (tobs, None)
    return sol, t_stop, (None if why is None else _REASONS[why])

def weights(model, p, b, ph, alpha):
    """Process weights (notebook definitions) for trajectories b, ph; rows follow PROCESSES[model]."""
    tobs = p["tobs"]
    w = [tobs * p["r"] * b / (alpha * b), tobs * p["a"] * ph * b / (alpha * b),
         tobs * p["c"] * p["a"] * ph * b / (alpha * ph), tobs * p["m"] * ph / (alpha * ph)]
    if model == "cc":
        w.append(tobs * p["r"] * b ** 2 / (alpha * b * p["K"]))
    return np.vstack(w)

def simulate(model, p, alphas, threshold=1.0):
    """One parameter sample p (dict: r, a, c, m, tobs, b0, p0[, K]) at several alphas -> list of result dicts."""
    t_eval = np.linspace(0, p["tobs"], N_POINTS)
    if model == "standard":
        full_f, bool_f = lotka_volterra_standard, theta_approach_model
        args = lambda al: (p["r"], p["a"], p["c"], p["m"], p["tobs"], al, threshold)
    else:
        full_f, bool_f = lotka_volterra_standard_cc, theta_approach_model_carrying_capacity
        args = lambda al: (p["r"], p["a"], p["c"], p["m"], p["K"], p["tobs"], al, threshold)
    nproc = len(PROCESSES[model])
    sf, tf, wf = _solve(full_f, [p["b0"], p["p0"]], p["tobs"], args(alphas[0]), t_eval)  # alpha-independent
    rows = []
    for al in alphas:
        row = dict(alpha=al)
        sb, tb, wb = _solve(bool_f, [p["b0"], p["p0"]], p["tobs"], args(al), t_eval)
        if sf.status < 0 or sb.status < 0:
            row.update(status="solver failed", msg=(sf.message if sf.status < 0 else sb.message)); rows.append(row); continue
        t_stop = min(tf, tb)
        if t_stop >= p["tobs"]:
            status, reason = "completed", ""
        else:
            status = "stopped"
            reason = ("full: " + wf) if tf <= tb else ("Boolean: " + wb)
        k = min(int(np.sum(t_eval <= t_stop)), sf.y.shape[1], sb.y.shape[1])  # points both models reached
        bf, pf, b, ph = sf.y[0, :k], sf.y[1, :k], sb.y[0, :k], sb.y[1, :k]
        with np.errstate(all="ignore"):
            com = (np.abs(bf - b) / bf + np.abs(pf - ph) / pf) / 2 * 100
            act = weights(model, p, b, ph, al) >= threshold      # Boolean-model activity, as in the notebook
        n_on = act.sum(axis=0)
        row.update(status=status, stop_reason=reason, t_stop=float(t_stop), frac_window=float(t_stop / p["tobs"]),
                   n_points=k, comparable=k >= 2, mean_err=float(np.mean(com)), median_err=float(np.median(com)),
                   max_err=float(np.max(com)), all_on_sustained=bool(act.all()), none_ever=bool(not act.any()),
                   **{f"frac_{j}_processes": float(np.mean(n_on == j)) for j in range(nproc + 1)})
        rows.append(row)
    return rows
