# Progress log

Current state of the transient-dynamics project. Keep this file short: update it at the end of each work session. Details live in the folder READMEs. This repository is public, so do not record reviewer comments or unpublished material here.

Last updated: 2026-11-06

## Status

| Area | State |
|---|---|
| Manuscript | Revision in progress (Cobo-López et al.; preprint: https://doi.org/10.1101/2023.12.27.573307) |
| Robustness analysis | Done: Figure 6 and Supplementary Figures S3–S6 (`results/alpha-robustness-v2/`); labels say "Boolean model" (PR #41); PNG, PDF, SVG and EPS (PR #42); panel b gives n and n_all, and the supplementary figures follow the revised order (2026-11-06) |
| Notebook | Carrying-capacity combined-error bug fixed and example outputs regenerated (PR #27) |
| Repository | Legacy material removed from the working tree (PR #28); preserved in the `legacy` branch and the `v0.1.0-alpha` pre-release |
| Releases | [`v0.2.3-alpha`](https://github.com/luquelab/transient-dynamics/releases/tag/v0.2.3-alpha) (2026-11-06): **version to cite for the resubmission to PLOS Computational Biology**, DOI [10.5281/zenodo.24014219](https://doi.org/10.5281/zenodo.24014219); relabeled and EPS figures, n_all counts and the revised supplementary order (PRs #41, #42, #44 and #46). Previous: `v0.2.2-alpha` ([10.5281/zenodo.23810547](https://doi.org/10.5281/zenodo.23810547)), `v0.2.1-alpha` ([10.5281/zenodo.23809339](https://doi.org/10.5281/zenodo.23809339)), `v0.2.0-alpha` (PRs #27–#30) and `v0.1.0-alpha` (before the robustness analysis). Concept DOI, always the latest release: [10.5281/zenodo.23809338](https://doi.org/10.5281/zenodo.23809338) |
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
  - S3–S5: robustness analysis (S3, stopped runs; S4, with carrying capacity; S5, both). Reordered on 2026-11-06; S3 and S4 were swapped before.
  - S6: extreme outlier.
- **Releases:** create tags and releases in the GitHub release UI; Zenodo archives each release and assigns it a version DOI. Update `version` in `CITATION.cff` before each release. The README badge shows the concept DOI, which always resolves to the latest release, so it needs no update after a release.
- **Stable path:** `results/alpha-robustness-v2` is cited in the Methods. Keep the path until acceptance, or update the citation.

## Next steps

1. **Finish the manuscript and rebuttal** (in progress). Cite the Zenodo DOI of `v0.2.3-alpha` (10.5281/zenodo.24014219) for the code, in the reference list, the Data Availability Statement and Table S6. After the revision, release again with the new bioRxiv version (or the published article) and remove the "under review" wording. Until the refactor:
   - do not restructure the notebook: the analysis scripts load its model functions by cell number;
   - if the notebook must change, rerun one analysis script to confirm it still reproduces its reference value;
   - make revision fixes on the current structure and release them as a patch (`v0.2.4-alpha`);
   - find the parameter-ranges spreadsheet (S1 Dataset), probably with Lucas or in the original submission files; add it to `data/` in a patch release if it should be archived with the code.
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
- **Citation metadata:** `CITATION.cff` (used by Zenodo and GitHub's "Cite this repository" button)
- **Project conventions:** `CLAUDE.md`

## Session log

Newest first. One short entry per work session.

- **2026-11-06:**
  - computed the values for the revised Results paragraph on the robustness analysis;
  - added the number of runs that spend all of their time with each number of relevant processes (n_all) to panel b of Figures 6 and S3–S5;
  - renumbered the supplementary figures to the revised order (S3, stopped runs; S4, with carrying capacity; S5, both) in `make_docx.js`, the README and this file (PR #44);
  - updated the supplementary figure files and the summary table for the revised manuscript to the new numbering;
  - released `v0.2.3-alpha` (PR #46), archived on Zenodo as 10.5281/zenodo.24014219, and recorded it here (PR #47);
  - kept the README badge on the concept DOI (it resolves to `v0.2.3-alpha`);
  - searched the repository, the `legacy` branch and the history for the two supplementary datasets: the simulation outputs are `data/transient-dynamics_output-current/` (packaged as `supplementary_file_simulations_outputs.zip` from `v0.2.3-alpha`); the parameter-ranges spreadsheet (`supplementary_file_parameter_ranges.xlsx`) is not in the repository.
- **2026-11-04:**
  - computed values and a supplementary summary table of the robustness analysis for the revised manuscript;
  - relabeled Figures 6 and S3–S6 from "transient model" to "Boolean model" (PR #41);
  - added EPS versions of the robustness figures, as the journal requires (PR #42);
  - planned `v0.2.3-alpha` for before resubmission.
- **2026-10-30:**
  - finished the manuscript Word documents for the robustness analysis, aligned with the final supplementary figures (S7 dropped), and added `CLAUDE.md` (PR #27);
  - removed legacy and superseded material (PR #28), added `PROGRESS.md` (PR #29) and the default GitHub assignee in `CLAUDE.md` (PR #30);
  - released `v0.2.0-alpha` and recorded it here (PR #31);
  - added `CITATION.cff` and the AI-assistance acknowledgment (PR #32), connected the repository to Zenodo, and released `v0.2.1-alpha`; added the DOI badge to the README (PR #34);
  - updated the preprint reference to bioRxiv version 3 (PR #35), added the PLOS Computational Biology revision note (PR #37), and released `v0.2.2-alpha`;
  - recorded the release DOIs and revision-period guidelines here (PRs #33, #38 and #39).
