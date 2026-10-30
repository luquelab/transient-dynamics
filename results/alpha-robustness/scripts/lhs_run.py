"""Latin hypercube (log-uniform) sweep of the standard model at alpha = 0.1.
Ranges: Table 1 of the manuscript (2026-10-26). usage: lhs_run.py <n> <seed>"""
import sys, os, signal, time, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import qmc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
MODEL = os.environ.get('LEAN_MODEL', 'lean_model')  # lean_model_ext: extinction threshold at 1e-3/ml
simulate = importlib.import_module(MODEL).simulate
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RANGES = {"r": (3e-5, 6), "a": (7e-10, 1e-6), "c": (2, 5000), "m": (1e-3, 0.8),
          "b0": (1, 8e9), "p0": (1, 2e10), "tobs": (0.8, 16800)}  # B0, P0 lower bound 1/ml (Table 1's 1e-4 is a typo)
TIMEOUT = 300  # s per sample; a stalled solve is flagged, not waited on

def _alarm(*_): raise TimeoutError
def run(row):
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(TIMEOUT); t0 = time.time()
    try:
        out = simulate(row["r"], row["a"], row["c"], row["m"], row["tobs"], row["b0"], row["p0"], alpha=0.1)
    except TimeoutError:
        out = dict(ok=False, msg=f"timeout > {TIMEOUT} s")
    except Exception as e:
        out = dict(ok=False, msg=f"{type(e).__name__}: {e}")
    finally:
        signal.alarm(0)
    return {**row, **out, "seconds": time.time() - t0}

if __name__ == "__main__":
    n, seed = int(sys.argv[1]), int(sys.argv[2])
    lo = np.log10([v[0] for v in RANGES.values()]); hi = np.log10([v[1] for v in RANGES.values()])
    u = qmc.LatinHypercube(d=len(RANGES), seed=seed).random(n)
    X = 10 ** qmc.scale(u, lo, hi)
    rows = [dict(sample=i, **dict(zip(RANGES, x))) for i, x in enumerate(X)]
    t0 = time.time()
    with Pool(4) as pool:
        res = pool.map(run, rows, chunksize=1)
    df = pd.DataFrame(res)
    df.to_csv(f"{S}/figures/lhs_b0p0min1{ {'lean_model_ext': '_ext', 'lean_model_stop': '_stop'}.get(MODEL, '') }_alpha0.1_n{n}_seed{seed}.csv", index=False)
    print(f"n={n} seed={seed}: {time.time()-t0:.0f} s wall; per-sample median {df.seconds.median():.2f} s, max {df.seconds.max():.1f} s")
    print("solver ok:", int(df.ok.sum()), "| failed:", int((~df.ok).sum()))
    if (~df.ok).any(): print(df.loc[~df.ok, "msg"].str[:80].value_counts().to_string())
