import pandas as pd, numpy as np, os, glob
S = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = ["growth", "predation", "burst", "decay"]
rows = []
for f in sorted(glob.glob(f"{S}/sweeps/*/alpha_*/generated-data/standard_model/results_standard/array_outputs.xlsx")):
    scen = f.split("/sweeps/")[1].split("/")[0]; a = f.split("/alpha_")[1].split("/")[0]
    d = pd.read_excel(f)
    act = np.column_stack([d[f"weight_{p}"].to_numpy() >= 1 for p in P])
    rows.append(dict(scenario=scen, alpha_str=a, alpha=float(a), mean_com=d.rel_error_com_per.mean(), median_com=d.rel_error_com_per.median(),
                     max_com=d.rel_error_com_per.max(), mean_bac=d.rel_error_bac_per.mean(), mean_pha=d.rel_error_phage_per.mean(),
                     all_on_always=bool(act.all()), none_ever=bool(not act.any()), avg_active=act.sum(1).mean(),
                     n_nan=int(d.rel_error_com_per.isna().sum())))
t = pd.DataFrame(rows).sort_values(["scenario", "alpha"])
t.to_csv(f"{S}/figures/all_scenarios_alpha_sweep_errors.csv", index=False)
pd.set_option("display.width", 200)
print(t.drop(columns="alpha_str").to_string(index=False, float_format=lambda x: f"{x:.3g}"))
