# Progress log

Current state of the transient-dynamics project. Keep this file short: update it at the end of each work session. Details live in the folder READMEs. This repository is public, so do not record reviewer comments or unpublished material here.

Last updated: 2026-10-30

## Status

| Area | State |
|---|---|
| Manuscript | Revision in progress (Cobo-López et al.; preprint: https://doi.org/10.1101/2023.12.27.573307) |
| Robustness analysis | Done: Figure 6 and Supplementary Figures S3–S6 (`results/alpha-robustness-v2/`) |
| Notebook | Carrying-capacity combined-error bug fixed and example outputs regenerated (PR #27) |
| Repository | Legacy material removed from the working tree (PR #28); preserved in the `legacy` branch and the `v0.1.0-alpha` pre-release |
| Releases | [`v0.2.2-alpha`](https://github.com/luquelab/transient-dynamics/releases/tag/v0.2.2-alpha) (2026-10-30): **version to cite for the PLOS Computational Biology revision**, DOI [10.5281/zenodo.23810547](https://doi.org/10.5281/zenodo.23810547); cites bioRxiv version 3 and states that it accompanies the revision (PRs #33–#37). Same code and analyses as `v0.2.1-alpha` ([10.5281/zenodo.23809339](https://doi.org/10.5281/zenodo.23809339)) and `v0.2.0-alpha` (PRs #27–#30). Concept DOI, always the latest release: [10.5281/zenodo.23809338](https://doi.org/10.5281/zenodo.23809338). Previous: `v0.1.0-alpha` (before the robustness analysis) |
| Refactor | Not started; planned as `v0.3.0-alpha` (see Next steps) |

**Versioning plan:** alpha while the code structure and interface change; beta once the package supports other models through a stable, tested interface; `v1.0.0` for a stable public interface.

## Key results (robustness analysis, α = 0.1)

Results are for the runs that completed the observational window.

| | Standard model | With carrying capacity |
|---|---|---|
| Runs within the 10% tolerance | 997 of 1,000 | 988 of 999 (+1 extreme outlier excluded) |
| Median mean combined relative error | 0.68% | 0.77% |
| Mean fraction of time with all processes relevant | 0.19 | 0.13 |
| All processes relevant for the whole window (share of sampled space) | 0.96% | 1.29% |

- **Across α = 10⁻⁵ to 1:** the error stays within the tolerance (100·α %) in at least 98.5% of completed runs, and it decreases as α decreases.
- **Runs that exceed the tolerance:** these are long oscillatory trajectories in which the Boolean model accumulates a persistent delay; their dynamics remain qualitatively correct (Figure S6).

## Decisions and conventions

- **Language:** American English (see `CLAUDE.md`).
- **Default observer tolerance:** α = 0.1 (10%).
- **Physical bounds:** a run stops when a population falls below 10⁻³ per ml or rises above 10¹⁵ bacteria or 10¹⁸ phage per ml. Both models are compared up to the stop time.
- **Sampling:** log-uniform Latin hypercube over the Table 1 ranges.
  - Initial populations: 10³–10¹² bacteria and 10³–10¹⁴ phage per ml.
  - Carrying capacity: K = 10⁴–10¹³ per ml, with B₀ ≤ K.
  - 1,000 completed runs per α.
- **Outliers:** Tukey's far-out test (above Q3 + 3·IQR) on log₁₀ error; outliers are excluded from the statistics and shown in ochre.
- **Figures:** Arial at 6–9 pt, 7.2 in wide, neutral grays. Red and blue are reserved for phage and bacteria. Panel titles use the style "**a)** Description".
- **Supplementary numbering:**
  - S1: Scenario 2, decay relevant.
  - S2: Scenario 4, transient.
  - S3–S5: robustness analysis (with carrying capacity; stopped runs; both).
  - S6: extreme outlier.
- **Stable path:** `results/alpha-robustness-v2` is cited in the Methods. Keep the path until acceptance, or update the citation.

## Next steps

1. **Finish the manuscript and rebuttal** (in progress). Cite the Zenodo DOI of `v0.2.2-alpha` (10.5281/zenodo.23810547) for the code. After the revision, release again with the new bioRxiv version (or the published article) and remove the "under review" wording. Until the refactor:
   - do not restructure the notebook: the analysis scripts load its model functions by cell number;
   - if the notebook must change, rerun one analysis script to confirm it still reproduces its reference value;
   - make revision fixes on the current structure and release them as a patch (`v0.2.3-alpha`).
2. **Refactor the code:**
   - move the model code into `src/transient_dynamics/` so analyses no longer run notebook cells by number;
   - add regression tests for the example scenarios;
   - switch the notebook to the package;
   - add `pyproject.toml` (`CITATION.cff` already exists).
3. **Explain the results theoretically:** give a heuristic argument for why the error tracks α, with oscillations as the amplifying case.
4. **Follow-up projects:**
   - test the approach on a larger system to show computational savings;
   - position it against tropical and dominant-balance model reduction;
   - consider a timing-aware error metric for oscillations.

## Where to find things

- **Notebook:** `notebooks/transient_dynamics.ipynb`
- **Example scenarios and data dictionary:** `data/transient-dynamics_output-current/`
- **Manuscript figures 1–5, S1 and S2:** `results/transient-dynamics-generated/`
- **Robustness analysis:** `results/alpha-robustness-v2/` (see its README)
- **Tutorial and changelog:** `docs/transient-dynamics/`
