import pandas as pd, numpy as np, re, glob, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lean_model import simulate
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))  # repository root
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = ROOT + "/data/transient-dynamics_output-current/example_scenarios/standard_model"
def inputs(scen):
    txt = open(f"{REPO}/{scen}/input.txt").read()
    return {k: float(eval(re.search(rf"^{k}\s*=\s*([^#\n]+)", txt, re.M).group(1))) for k in ["r","a","c","m","tobs","bi","pi"]}
cases = [(s, 0.1, f"{REPO}/{s}/results_standard/array_outputs.xlsx") for s in sorted(os.listdir(REPO)) if os.path.isdir(f"{REPO}/{s}")]
cases += [(f.split("/sweeps/")[1].split("/")[0], float(f.split("/alpha_")[1].split("/")[0]), f)
          for f in sorted(glob.glob(f"{S}/sweeps/*/alpha_*/generated-data/standard_model/results_standard/array_outputs.xlsx"))]
worst = 0; t0 = time.time()
for scen, alpha, f in cases:
    v = inputs(scen); ref = pd.read_excel(f).rel_error_com_per.mean()
    got = simulate(v["r"], v["a"], v["c"], v["m"], v["tobs"], v["bi"], v["pi"], alpha=alpha)["mean_com"]
    rel = abs(got - ref) / ref if ref else abs(got)
    worst = max(worst, rel)
    if rel > 1e-9: print(f"MISMATCH {scen} alpha={alpha}: lean {got:.6g} vs notebook {ref:.6g}")
print(f"{len(cases)} cases checked (6 from Lucas + {len(cases)-6} from our sweeps); worst relative difference {worst:.2e}; {time.time()-t0:.0f} s incl. Excel reads")
