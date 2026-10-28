import os
# solve-only cost: same models, solver and tolerances as the notebook, no figures or Excel
import time, numpy as np, json
from scipy.integrate import solve_ivp
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))  # repository root
cells = json.load(open(ROOT + "/notebooks/transient_dynamics.ipynb"))["cells"]
exec("import numpy as np\n" + "".join(cells[14]["source"]))  # theta_approach_model, lotka_volterra_standard
for name, para, y0 in [("scenario 1", (0.9, 3e-8, 50, 1.5e-3, 14, 0.1, 1.0), [1e3, 1e4]),
                       ("scenario 3b", (0.1, 3e-9, 50, 0.1, 300, 0.1, 1.0), [2.65e3, 4.93e4])]:
    t_eval = np.linspace(0, para[4], 10000)
    t0 = time.perf_counter()
    for f in (theta_approach_model, lotka_volterra_standard):
        solve_ivp(f, (0, para[4]), y0, args=para, method="LSODA", t_eval=t_eval, rtol=1e-13, atol=1e-16)
    print(f"{name}: full + Boolean solve = {time.perf_counter() - t0:.2f} s")
