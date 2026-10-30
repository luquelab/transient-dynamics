import json, re, sys, os, subprocess, csv
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))  # repository root
REPO = ROOT + ""
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
nb = json.load(open(f"{REPO}/notebooks/transient_dynamics.ipynb"))
cells = nb["cells"]
rows = list(csv.DictReader(open(f"{REPO}/data/transient-dynamics_output-current/example_scenarios/standard_model/scenario_inputs.csv")))
for row in rows:
    name = row["folder"]
    # parse the scenario's own input.txt (source of truth) for exact expressions
    txt = open(f"{REPO}/data/transient-dynamics_output-current/example_scenarios/standard_model/{name}/input.txt").read()
    params = {}
    for k in ["r","a","c","m","tobs","bi","pi","reference_level","thresholdlevel"]:
        mt = re.search(rf"^{k}\s*=\s*([^#\n]+)", txt, re.M)
        params[k] = mt.group(1).strip()
    out = f"{S}/runs/{name}"
    os.makedirs(out, exist_ok=True)
    code = ["import matplotlib\nmatplotlib.use('Agg')\n"]
    for i, c in enumerate(cells[:57]):
        if c["cell_type"] != "code":
            continue
        src = "".join(c["source"])
        if i == 9:
            src = src.replace('parent_dir = os.path.join("..", "data", "transient-dynamics_output-current", "data_output")',
                              f'parent_dir = {out!r}')
            assert repr(out) in src
        if i == 17:
            for k, v in params.items():
                src, n = re.subn(rf"^{k}\s*=.*$", f"{k} = {v}", src, count=1, flags=re.M)
                assert n == 1, (name, k)
        code.append(f"# ---- cell {i}\n{src}\n")
    open(f"{out}/run.py", "w").write("\n".join(code))
    r = subprocess.run([sys.executable, "run.py"], cwd=out, capture_output=True, text=True)
    open(f"{out}/log.txt", "w").write(r.stdout + "\n--- stderr\n" + r.stderr)
    print(name, "exit", r.returncode, params, r.stderr.strip().splitlines()[-1:] if r.returncode else "")
