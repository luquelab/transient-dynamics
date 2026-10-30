# WHAT
Revised robustness analysis of the adaptive Boolean model against the full model, for the standard phage-bacteria (Lotka-Volterra) model and the model with bacterial carrying capacity. It supersedes the exploratory analysis in `../alpha-robustness/`.

# WHEN
2026-10-29

# DESIGN
**Models.** Both models are taken unchanged from `notebooks/transient_dynamics.ipynb`:
- standard model: code cell 14, with 4 processes (growth, predation, burst, decay);
- carrying-capacity model: code cell 60, with 5 processes (the four above plus the carrying-capacity term r·B²/K).

The solver (LSODA, rtol 10⁻¹³, atol 10⁻¹⁶), the 10,000-point time grid and the weights (a process is relevant when its weight is ≥ 1) are the notebook's.

**Error.** The mean combined relative error is the time average of ½(|B_full − B_Boolean|/B_full + |P_full − P_Boolean|/P_full) × 100.

**Parameter ranges.** Sampling is log-uniform for r, a, c, m and τ over the Table 1 ranges. The other inputs use these ranges:

| Input | Range |
|---|---|
| B₀ | 10³–10¹² per ml |
| P₀ | 10³–10¹⁴ per ml |
| K (carrying-capacity model only) | 10⁴–10¹³ per ml |

B₀ ≤ K is imposed in the carrying-capacity model: B₀ is sampled log-uniformly on [10³, min(K, 10¹²)] given K.

**Stopping rules.** A run stops when either population, in either model, leaves these bounds:
- lower bound: 10⁻³ per ml for both populations;
- upper bound: 10¹⁵ per ml for bacteria and 10¹⁸ per ml for phage. These are the smallest bacterium (about 10⁻³ µm³) and the smallest phage (about 10⁻⁶ µm³) packed into 1 ml.

The two models are compared up to the earlier stop time. A run is "completed" if it lasts the whole observation window τ. "Stopped" runs keep the portion up to the stop. Stopped runs with fewer than 2 compared time points (collapse within τ/10,000) are counted but not analysed.

**Sampling.** α takes the values 10⁻⁵, 10⁻⁴, …, 1 (decades), and the same parameter samples are used at every α. The Latin hypercube size doubles (2,000 → 4,000 → …) until every α has at least 1,000 completed runs; this took 16,000 samples for the standard model and 8,000 for the carrying-capacity model. For each α, exactly 1,000 completed runs are then drawn at random (seed 2026). All analysable stopped runs are used for the supplementary analysis.

**Fraction of sampled space.** The share of the sampled space where all processes stay relevant for the whole window is: completed runs with every process relevant at every time point, divided by all samples. It is reported with a Wilson 95% confidence interval.

# RESULTS (α = 0.1)
| | Standard, completed | Standard, stopped | Carrying capacity, completed | Carrying capacity, stopped |
|---|---|---|---|---|
| Runs analysed | 1,000 | 6,290 | 999 (+1 extreme outlier excluded) | 3,513 |
| Median error | 0.68% | 0.0035% | 0.77% | 0.0065% |
| 90th percentile | 2.2% | 0.28% | 2.9% | 0.38% |
| Maximum | 13% | 12% | 276% | 12% |
| Runs within 10% | 997 | 6,289 | 988 | 3,511 |
| Mean fraction of time with all processes relevant | 0.19 | 0.27 | 0.13 | 0.08 |
| All processes relevant for the whole window | 9.5% (7.8–11.5%) | — | 6.3% (5.0–8.0%) | — |
| … as a share of all samples | 0.96% (0.82–1.13%) | | 1.29% (1.06–1.56%) | |

**Extreme outliers.** Within each model, α and set, a run is flagged by Tukey's far-out test on log₁₀(mean error): it is an outlier if the value is above Q3 + 3·IQR. Flagged runs are excluded from all statistics and shown in ochre in the figures.

Only 5 runs are flagged, all in the carrying-capacity model and all completed: 1 at α = 0.1 (8.6×10⁵%) and 4 at α = 1. They are the runs that pushed the means above the tolerance line:
- at α = 0.1, the carrying-capacity mean is 1.65% without the outlier and 865% with it;
- at α = 1, it is 13.6% without the outliers and 9.3×10⁵% with them.

**What the runs above tolerance have in common** (`scripts/investigate_outliers.py`, `tables/outlier_diagnostics.csv`). Each was re-simulated and compared with 150 random runs within tolerance.
- Runs that completed the window are long-window oscillations. The median observation time is 5,000–10,000 h, against about 10–20 h for runs within tolerance. Over that window the phage oscillate 15–70 times and the set of relevant processes switches 30–190 times, against 0–1 switches for runs within tolerance. A small timing shift between the Boolean and full models builds up from cycle to cycle.
- Their troughs bring populations down to a few individuals per ml, still above the 10⁻³ floor. Relative errors at these troughs are very large (the worst run reaches 0.7 vs 3.4×10⁶ phage per ml).
- None comes close to the upper bounds: maxima are about 10⁵–10⁷ bacteria and 10⁸–10⁹ phage per ml.
- The exception is the standard model at α = 1. Its 3 runs above tolerance are short (about 3 h) with fast growth and decay, and the coarse tolerance simply ignores sizeable changes.
- Stopped runs above tolerance get most of their error (55–75%) from time points below 1 per ml, just before they hit the floor.

**Across all α:** the mean error stays at or below the observer tolerance (100·α %) in at least 98.5% of completed runs, for every α and both models (column `frac_err_within_alpha` in the summary tables). As α decreases, the fraction of time with all processes relevant rises towards 1, and the Boolean model converges to the full model.

**Outcome counts:**

| | Standard model | Carrying-capacity model |
|---|---|---|
| Samples | 16,000 | 8,000 |
| Completed (α ≤ 0.1) | 1,689 (10.6%) | about 1,505 (18.8%) |
| Stopped before the 2nd time point | 8,021 | 2,980 |
| Solver failures | 0 | 4 samples (solver timeout > 900 s) |
| Upper bound reached | 0 | 0 |

All stopped runs stop at the lower bound.

# PRE-RELEVANCE DIAGNOSTICS (α = 0.1)
**Hypothesis tested.** Runs above tolerance spend a long time with a process just below the threshold (pre-relevant, 0.1 ≤ w < 1). Such a process is ignored, yet its effect accumulates.

**Approach.** `scripts/prerelevance.py` computes pre-relevance metrics from the transient-model weights for all 2,000 completed runs; results are in `tables/prerelevance_alpha0.1.csv`. `scripts/outlier_trajectories.py` draws, for each of the 15 runs above tolerance:
- the trajectories (full model vs transient model),
- the combined error over time,
- the number of relevant processes in each model,
- each process's weight over time, with the pre-relevant zone shaded.

Runs with more than 20 switches get an extra zoomed view of their last 10 switches. These are in `figures/outlier_diagnostics/`, with all runs in `all_runs_above_tolerance.pdf` and the summary in `summary_prerelevance_test.png`.

**Findings.**
- **Pre-relevance sets the general error level.** Across all runs, the error correlates with the neglected-effect index (Spearman ρ = 0.70–0.71) and with the time spent pre-relevant (ρ = 0.61–0.64).
- **Pre-relevance alone does not single out the runs above tolerance.** Many runs within tolerance also spend their whole window with a process pre-relevant (AUC 0.43–0.55). The runs above tolerance are instead marked by many process switches (AUC 0.87) and long windows (AUC 0.88–0.90).
- **The trajectories show that a pre-relevant process causes trouble when the dynamics amplify its neglected effect.** Two mechanisms appear:
  - *Sustained pre-relevance with a long decline* (e.g. standard run 7773). Phage decay stays at w ≈ 0.55 for the whole window and is ignored. The transient model keeps slightly more phage, the extra predation makes its bacteria decline slightly faster, and the gap compounds over ten orders of magnitude of decline.
  - *Oscillation phase drift* (e.g. carrying-capacity run 682, 351 switches). The carrying-capacity term reaches the pre-relevant zone (w ≈ 0.4) at every bacterial peak, and predation crosses the threshold every cycle. The transient model lags slightly in each cycle, and the lag builds up to about 10 h after about 170 cycles.
- **A lower tolerance fixes these runs.** With α = 0.01, 14 of the 15 runs are within 10%; with α = 0.001, all 15 are. This supports the proposed remedy: monitor pre-relevant processes and recommend a lower α. The trigger should combine pre-relevance with amplification (many switches or cycles, or a long monotonic decline), rather than use time spent pre-relevant alone.

# FILES
**`scripts/`**
- `models.py`: the two models, the stopping bounds and the per-run metrics.
- `run_lhs.py <standard|cc>`: the progressive Latin hypercube sweep.
- `analyze.py`: the per-α statistics and the analysis sets.
- `figure.py <standard|cc> <completed|stopped>`: the manuscript figures.
- `investigate_outliers.py`: re-simulates the runs above tolerance and compares their trajectories with runs within tolerance.
- `prerelevance.py`, `outlier_trajectories.py`: the pre-relevance test and the per-run diagnostic figures.
- `si_outlier_figure.py <model> <sample> <stem>`: the supplementary outlier figures (S6: carrying-capacity run 178; S7: run 682).
- `make_docx.js`: builds the Word documents in `manuscript/` (main text: Results, Methods, Figure 6; supplementary: Figures S3–S7). Run with `node make_docx.js`.

**`tables/`**
- `lhs_*_seed2026.csv.gz`: every sample × α, the raw sweep.
- `*_analysis_sets.csv.gz`: the runs used in each analysis set.
- `*_summary_by_alpha.csv`: per-α statistics for the error (mean, SD, variance, median, IQR, 90th percentile, max), the fraction of time with k = 0…P relevant processes, and the counts of runs with all processes relevant (with 95% CIs).
- `*_outcome_counts.csv`: completed, stopped (by reason), early-stopped and failed runs.

**`figures/`** contains `<model>_<set>_figure.*` (three panels, 7.2 in wide) and each panel separately (`_panel_a/b/c`), in PNG (300 dpi), PDF and SVG. The panels are:
- a: error distribution at α = 0.1;
- b: fraction of time with each number of relevant processes at α = 0.1 (violins and 200 random runs; n at the top is the number of runs that spend any time in that state);
- c: the α sweep, with errors at the top (box plots, mean ◆, tolerance line, extreme outliers in ochre) and the mean ± SD fraction of time with all processes relevant at the bottom.

Each panel letter is followed by a short description. Colours are neutral greys matched to the manuscript figures (`#7f8080`, `#505151`, `#231f20`), with ochre `#e5a41c` (also from the manuscript palette) for extreme outliers; red and blue are left for phage and bacteria. Fonts are Arial at 6–9 pt. The SVG text names Arial and stays editable. PNG and PDF are rendered and embedded with Liberation Sans, which has the same letter widths as Arial (Arial is not available in the cloud environment).

**Rerun** (from `scripts/`; package versions as in `../alpha-robustness/requirements.txt`):

```bash
python run_lhs.py standard 2026
```

```bash
python run_lhs.py cc 2026
```

```bash
python analyze.py standard ../tables/lhs_standard_n16000_seed2026.csv.gz
```

```bash
python analyze.py cc ../tables/lhs_cc_n8000_seed2026.csv.gz
```

```bash
python figure.py standard completed
```

Run `figure.py` likewise for the other three combinations of model and set.

# NOTE ON THE NOTEBOOK
Before 2026-10-29, the notebook's carrying-capacity combined error (`rel_error_com_per` in `array_outputs_cc`) was wrong. Code cell 71 computed it as `relative_error_percentage_finder(p_rel_error_percent_cc, ...)`, which is the phage error multiplied by 100 again. The bacteria and phage error columns were correct.

Fixed on 2026-10-29. The carrying-capacity section now calls `combined_error_finder` (the average of the bacteria and phage errors, as in the standard model). `combined_error_finder` itself also used to append to the global `p_com_relative_error_percent` instead of a fresh list, so calling it a second time would have doubled the list; it now builds its own list.

The fix was verified by running the notebook through the carrying-capacity section (τ = 300 h): the combined error is 0.9750869526%, matching this analysis. The standard-model output is unchanged and identical, column by column, to the stored scenario 3 (transient) results.

This analysis computes the error directly and was never affected.
