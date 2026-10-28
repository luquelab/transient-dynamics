import json, re, os, subprocess, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))  # repository root
REPO = ROOT + ""
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cells = json.load(open(f"{REPO}/notebooks/transient_dynamics.ipynb"))["cells"]
txt = open(f"{REPO}/data/transient-dynamics_output-current/example_scenarios/standard_model/scenario_1/input.txt").read()
base = {k: re.search(rf"^{k}\s*=\s*([^#\n]+)", txt, re.M).group(1).strip()
        for k in ["r","a","c","m","tobs","bi","pi","reference_level","thresholdlevel"]}
for alpha in (sys.argv[1:] or ["0.001", "0.01", "0.1", "0.5"]):
    p = dict(base, reference_level=alpha)
    out = f"{S}/alpha_sweep/alpha_{alpha}"
    os.makedirs(out, exist_ok=True)
    code = ["import matplotlib\nmatplotlib.use('Agg')\n"]
    for i, c in enumerate(cells[:50]):  # standard model through output storage
        if c["cell_type"] != "code": continue
        src = "".join(c["source"])
        if i == 9:
            src = src.replace('os.path.join("..", "data", "transient-dynamics_output-current", "data_output")', repr(out))
        if i == 17:
            for k, v in p.items():
                src, n = re.subn(rf"^{k}\s*=.*$", f"{k} = {v}", src, count=1, flags=re.M); assert n == 1
        # plot-only guard: axis limits are hard-coded for alpha in {0.1, 1}; fall back to a wide range otherwise
        src = src.replace('raise ValueError("Unsupported reference_level")', 'y_lim = (1e-6, 1e6); y_ticks = [1e-4, 1, 1e4]')
        # print-only helper assumes both models have the same number of tipping points
        src = src.replace("for i in range(0, len(tp_full)):", "for i in range(0, min(len(tp_full), len(tp))):")
        code.append(f"# ---- cell {i}\n{src}\n")
    open(f"{out}/run.py", "w").write("\n".join(code))
    r = subprocess.run([sys.executable, "run.py"], cwd=out, capture_output=True, text=True)
    open(f"{out}/log.txt", "w").write(r.stdout + "\n--- stderr\n" + r.stderr)
    print(alpha, "exit", r.returncode, r.stderr.strip().splitlines()[-1:] if r.returncode else "")
