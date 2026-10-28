# WHAT
Robustness analysis of the adaptive Boolean model against the full Lotka-Volterra (standard) model of `notebooks/transient_dynamics.ipynb`. It covers:

1. Reproducing the six `standard_model` example scenarios in `data/transient-dynamics_output-current/example_scenarios/`.
2. Sweeping the observer tolerance α for each scenario, including the two limit values of α (all processes observable; no process observable).
3. Latin hypercube sampling (LHS) of the parameter space in Table 1 of the manuscript at α = 0.1, with successive criteria for populations that collapse or explode.

All simulations use the notebook's own model functions (code cell 14, loaded from the notebook at run time), the same solver (`solve_ivp`, LSODA, `rtol = 1e-13`, `atol = 1e-16`), the same 10,000-point time grid, and the same weight and error definitions:
- Relative error: |full − Boolean| / full × 100, for bacteria and for phage.
- Combined error: the average of the bacteria and phage errors at each time point.
- Mean combined relative error: the average of the combined error over the time points compared.

# WHEN
Started: 2026-10-28

# ENVIRONMENT
Python 3.11 with the package versions in `requirements.txt`. Run every script from the `scripts/` folder:

```bash
pip install -r requirements.txt
```

```bash
cd results/alpha-robustness/scripts
```

# FILES & FOLDERS
FOLDER: /scripts
--> Analysis code. Outputs go to `../figures`. The full notebook runs (per-run figures and Excel files) go to `../sweeps`, `../alpha_sweep` and `../runs`. Those folders are git-ignored: they are large and can be regenerated.

FOLDER: /figures
--> Figures (PNG), per-run tables (CSV) and the α limits (JSON) described below.

# 1. Reproducing the example scenarios
- `build_and_run.py`: runs the notebook's standard-model cells with each scenario's `input.txt` (outputs to `../runs`).
- `compare.py`: compares those runs with the stored outputs.
- `validate_lean.py`: checks that the stripped-down simulator in `lean_model.py` (no per-run figures or Excel files) reproduces the notebook. The worst relative difference in the mean combined error is 3×10⁻¹¹.

Result: all six scenarios are reproduced. Five match exactly. Scenario 1 differs only at floating-point noise level (relative differences of about 10⁻¹¹ in the full-model columns).

# 2. α sweeps per scenario
Weights scale as W/α, where W does not depend on α:
- growth: τ·r
- decay: τ·m
- predation: τ·a·P(t)
- burst: τ·c·a·B(t)

This defines two limits (`limits.py`, saved to `figures/alpha_limits.json`):
- **α_all** = min over time and processes of W along the full trajectory. Below it, every process stays observable and the Boolean model equals the full model (error 0).
- **α_none** = max over processes of W at t = 0. Above it, no process is ever observable and the Boolean model stays at its starting populations.

| Scenario | α_all | α_none |
|---|---|---|
| 1 | 3.8×10⁻¹⁴ | 12.6 |
| 2 | 1.1×10⁻⁸ | 392 |
| 3 quasi-equilibrium | 5.37 | 30 |
| 3 transient | 4.8×10⁻⁴ | 30 |
| 4 quasi-static | 4.2×10⁻⁴ | 0.0195 |
| 4 transient | 4.2×10⁻⁴ | 16.6 |

Scripts:
- Scenario 1: `alpha_sweep.py <alpha> ...` runs the notebook; then `plot_alpha.py`, `plot_active.py` and `plot_pops.py` make `scenario1_*`.
- Other scenarios: `sweep_scenario.py <scenario_folder> <alpha> ...` runs the notebook; then `summarize.py` and `plot_scenario.py` make `scenario_*` and `all_scenarios_alpha_sweep_errors.csv`. `limit_runs.txt` lists the limit values that were run.

To run the notebook for α values other than 0.1 and 1, the scripts make two changes to their copy of the code, leaving the notebook untouched:
- the weight plot gets fallback axis limits;
- the tipping-point print helper tolerates different numbers of tipping points in the two models.

Neither change affects the model or the outputs.

# 3. Latin hypercube robustness at α = 0.1
The design is log-uniform over r, a, c, m, B₀, P₀ and τ (`lhs_run.py <n> <seed>`), using the ranges in Table 1 of the manuscript. The B₀ and P₀ ranges start at 1 per ml; the 10⁻⁴ lower bound in Table 1 is a typo. Each step below was a decision made during the analysis:

| Step | Treatment of extreme populations | Files (`figures/`) |
|---|---|---|
| a | None. B₀ and P₀ from 10⁻⁴; runs with negative populations or overflow flagged | `lhs_alpha0.1_*` (not bounded) |
| b | Exclude runs leaving [10⁻³, 10¹³] per ml; B₀ and P₀ from 10⁻⁴ | `lhs_alpha0.1_*bounded*` |
| c | As b, B₀ and P₀ from 1 per ml | `lhs_b0p0min1_alpha0.1_*` |
| d | Absorbing extinction at 10⁻³ per ml (set to 0), ceiling 10¹³ | `lhs_b0p0min1_ext_*` (`lean_model_ext.py`) |
| **e** | **Stop the run when a population in either model leaves [10⁻³, 10¹⁵] per ml; compare both models up to the earlier stop** | `lhs_b0p0min1_stop_*` (`lean_model_stop.py`) |

The scripts are in their final state (step e). Steps a to c were produced with earlier settings of `lhs_run.py` (B₀ and P₀ lower bound), `lhs_classify.py` (bounds) and `plot_lhs.py` (file names). `check_extinction.py` re-solves the excluded runs in log variables. It confirms that falling below the floor is real model behaviour (636 of 638 runs), not a solver artefact.

## Main results (step e)
- **2,000 samples, seed 2026** (`*_n2000_*`):
  - Outcomes: 538 runs completed the window, 1,145 stopped at the floor, and 317 stopped before the 2nd time point (nothing to compare). None hit the ceiling.
  - Errors across plotted runs: median 0.042%, 90th percentile 1.2%, maximum 8.9%.
- **Balanced sample, 500 completed + 500 comparable stopped runs** (`*_n1000balanced2_*`):
  - Construction: drawn at random (seed 2026) from the 2,000 runs. Stopped runs come only from those that got past the 1st time point.
  - Errors: median 0.10%, 90th percentile 1.6%, maximum 8.9%. No run exceeds 10%.
  - This is a design choice that over-represents completed runs. Report it alongside the 2,000-run figures.
- The earlier balanced draw that included early-collapsing runs is kept as `*_n1000balanced_*`.

Classified CSVs have one row per sample, including:
- the parameters,
- `stop`: why the run ended,
- `t_stop`: when it ended,
- `n_points`: the number of time points compared,
- the error summaries,
- the process-activity flags.
