"""Log-uniform Latin hypercube sweep of the full vs. Boolean models over alpha = 1e-5 ... 1 (decades).
The design size doubles (2,000, 4,000, ...) until every alpha has at least TARGET completed runs.
usage: run_lhs.py <standard|cc> [seed]"""
import os, sys, signal, time, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import qmc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import models

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ALPHAS = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
TARGET, SIZES, TIMEOUT = 1000, [2000, 4000, 8000, 16000], 900
# Table 1 ranges; initial populations and K as revised (see README)
RANGES = {"r": (3e-5, 6), "a": (7e-10, 1e-6), "c": (2, 5000), "m": (1e-3, 0.8), "tobs": (0.8, 16800),
          "b0": (1e3, 1e12), "p0": (1e3, 1e14)}
K_RANGE = (1e4, 1e13)

def design(model, n, seed):
    names = list(RANGES) + (["K"] if model == "cc" else [])
    u = qmc.LatinHypercube(d=len(names), seed=seed).random(n)
    lo = {k: np.log10(v[0]) for k, v in {**RANGES, "K": K_RANGE}.items()}
    hi = {k: np.log10(v[1]) for k, v in {**RANGES, "K": K_RANGE}.items()}
    X = {k: 10 ** (lo[k] + u[:, j] * (hi[k] - lo[k])) for j, k in enumerate(names)}
    if model == "cc":  # constraint B0 <= K: B0 log-uniform on [1e3, min(K, 1e12)] given K
        top = np.log10(np.minimum(X["K"], RANGES["b0"][1]))
        X["b0"] = 10 ** (lo["b0"] + u[:, names.index("b0")] * (top - lo["b0"]))
    return [dict(sample=i, **{k: float(X[k][i]) for k in names}) for i in range(n)]

def _alarm(*_): raise TimeoutError
def run_one(arg):
    model, p = arg
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(TIMEOUT); t0 = time.time()
    try:
        rows = models.simulate(model, p, ALPHAS)
    except TimeoutError:
        rows = [dict(alpha=a, status="solver failed", msg=f"timeout > {TIMEOUT} s") for a in ALPHAS]
    except Exception as e:
        rows = [dict(alpha=a, status="solver failed", msg=f"{type(e).__name__}: {e}") for a in ALPHAS]
    finally:
        signal.alarm(0)
    dt = time.time() - t0
    return [{**p, **r, "seconds": dt} for r in rows]

if __name__ == "__main__":
    model = sys.argv[1]; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 2026
    for n in SIZES:
        t0 = time.time()
        with Pool(4) as pool:
            res = pool.map(run_one, [(model, p) for p in design(model, n, seed)], chunksize=4)
        df = pd.DataFrame([r for rows in res for r in rows])
        done = df[df.status == "completed"].groupby("alpha").size().reindex(ALPHAS, fill_value=0)
        print(f"{model} n={n}: {time.time() - t0:.0f} s | completed per alpha: " + ", ".join(f"{a:g}:{c}" for a, c in done.items()), flush=True)
        if done.min() >= TARGET:
            break
    out = os.path.join(HERE, "tables", f"lhs_{model}_n{n}_seed{seed}.csv.gz")
    df.to_csv(out, index=False)
    print("saved", out, "| solver failures:", int((df.status == "solver failed").sum()), flush=True)
