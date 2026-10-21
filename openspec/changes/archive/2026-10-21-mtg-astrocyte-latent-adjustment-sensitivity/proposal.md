# Proposal

## Why

In both cohorts the CV-selected rank-1 component of the primary analysis is a donor-level library-quality axis. `axis-diagnostics` gives it Spearman −0.84 with `ribosomal_mean`, 0.73–0.79 with the mitochondrial fractions and +0.49 with stage (`reports/motco-astrocyte-technical-axis-2026-10-19.md`). The 0.8.0 covariates adjust the stage means within the basis, but they never touch which basis PLS-DA chooses. The upstream fixes were label-blind and failed: the `residualized` variant removed the stage signal along with the axis. The trend-normalized candidates did not meet `ADOPTION_RULE` in both cohorts (`reports/motco-astrocyte-trend-normalization-2026-10-20.md` §4).

MOTCO 0.9.0 (`covariate-aware-latent-construction`) adds `covariate_latent_adjustment`, which fits the basis on features residualized on the declared covariates conditional on stage. Its Type I check returned `adopt_stage`, and its power comparison returned `candidate_fails` for both conditionings. MOTCO's binding conditions (`../MOTCO/docs/roadmap.md` Phase 6; `docs/api/analysis.md` "Latent adjustment") are:

- use `stage` conditioning only;
- use it only as a sensitivity analysis beside the 0.8.0 covariate-adjusted analysis;
- justify the covariate set in the case study's own proposal, which carries the mediator limitation.

This change is that proposal.

## What Changes

- **MOTCO pin moves from `v0.8.0` to `v0.9.0`** (`uv add "motco @ git+https://github.com/tomszar/MOTCO@v0.9.0"`).
  - MOTCO guarantees that every output without the new setting equals 0.8.0's except the runtime. The reference checks (`_check_reference`) compare the primary's recorded MOTCO version with the installed one, so the primary `analyze` is rerun under 0.9.0.
  - Before the old directories are replaced, the rerun is checked against a rerun under MOTCO 0.8.0 on the same machine (the pre-change code in its own locked environment). `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and `analysis_metadata.csv` must be byte-identical. `result.json` must be equal apart from the runtime and the MOTCO version.
  - The rerun is also compared with the kept 0.8.0 directories, which may have been made on another machine. That comparison is recorded in the report, not a gate: floating-point results are not byte-reproducible across machines.
  - A same-machine difference stops the change, and the 0.8.0 directories are restored. Otherwise the 0.9.0 rerun becomes `results/motco_analysis/`, the reference of the arms (its `result.json` hashes are new and reported); the kept directories stay as `results/motco_analysis_0.8.0/` until the report is committed; and the cluster's copy is replaced, with its old one kept as `results/motco_analysis_0.8.0/`.
  - Kept outputs of `rank-sensitivity`, `preprocessing-sensitivity` and `axis-diagnostics` made under 0.8.0 stay as they are. They are rerun only if they are needed again.
- **New `latent-sensitivity` step** with a predeclared registry `LATENT_SENSITIVITY_ARMS`. Per cohort (`cohort_30` primary, `cohort_50` sensitivity) it runs two arms on the **unchanged primary export**:
  - `extended_covariates`: the six `ANALYSIS_COVARIATES` plus `ribosomal_pseudobulk` and `mitochondrial_pseudobulk` from `results/pseudobulk/donor_technical.csv`, with the covariate-blind basis (no latent adjustment).
  - `latent_adjusted`: the same eight covariates with `covariate_latent_adjustment = "stage"`.
  - Everything else equals `analyze`: blocks, stage order, orientation, CV rank selection, `ANALYSIS_SEED`, `ANALYSIS_PERMUTATIONS`, `ANALYSIS_MISSING_COVARIATES = error`, and `n_jobs = 1` for reported results.
  - The primary analysis (six covariates, covariate-blind basis) is the reference. It is read and checked as `rank-sensitivity` and `preprocessing-sensitivity` do, never rerun or modified.
  - The three-way layout separates two effects. Primary → `extended_covariates` is the effect of the larger covariate set. `extended_covariates` → `latent_adjusted` is the effect of the basis choice alone, since MOTCO guarantees that the trajectory model, LS-means and RRPP are identical between those two and only the basis and the CV rank change.
- **The covariate set is a constant**, `LATENT_SENSITIVITY_COVARIATES = ANALYSIS_COVARIATES + ("ribosomal_pseudobulk", "mitochondrial_pseudobulk")`. The two fractions enter on the raw fraction scale used by the `residualized` variant.
  - Justification: they are the two fractions the technical-axis report classified as residualization candidates (between-donor share dominant, within-donor supertype share small). `hemoglobin` stays excluded as negligible (`RESIDUALIZATION_CELL_STATE_EXCLUSIONS`).
  - The export's `log10_lib_<omic>` stay out. The nuclei-count covariates already enter, and MOTCO's Type I check used eight covariates.
  - The set is never a flag and never chosen from results.
- **Predeclared mediator check, written before the runs and reported, never a gate.** Per cohort, `mediator_check.csv` gives each added covariate's Spearman correlation with stage (ordinal), overall and within each sex. Cells at `|ρ| ≥ DIAGNOSTIC_CORRELATION_THRESHOLD` are flagged as mediator candidates.
  - No simulation can tell a confounder from a mediator. So a flag removes nothing, and it does not change which arm runs.
  - A flag binds the report to say that the `latent_adjusted` basis may have removed stage signal carried by that fraction.
  - No p-value is computed.
- **Outputs** go to `results/motco_latent_sensitivity/`:
  - `<arm>/motco_analysis/<cohort>/` per arm, the same files `analyze` writes. MOTCO's `latent_adjustment` block, including `feature_r2`, is in `result.json`, and the technical table's hash is in the manifest. This mirrors the primary layout, so `axis-diagnostics --analysis-dir` and `rank-sensitivity --analysis-dir` can read an arm unchanged.
  - `comparison/<cohort>/` with `latent_sensitivity.csv` (reference row plus one row per arm: statistics, p-values, rejections, null summaries, selected rank, CV AUROC, covariates, conditioning), `latent_sensitivity.md`, `mediator_check.csv` and `latent_sensitivity_manifest.json`.
- **A committed report**, `reports/motco-astrocyte-latent-adjustment-<date>.md`. Its numbers are pasted from the generated files and checked against every `result.json`. It includes `axis-diagnostics` on both arms (the first latent axis's technical correlations, a measurement-space finding), and it binds the reading:
  - The primary analysis stays primary.
  - A rejection status that differs between an arm and the primary is covariate-set or basis dependence, not a result of the arm.
  - A different CV rank is reported beside `angle`.
  - The non-rejection rules hold in every arm.
  - Both MOTCO reports' conditions are carried as limitations: `shape` power loss with a nuisance axis, small-magnitude power loss without one, `candidate_fails`, the mediator limitation, and Phase 6 design-only validation.
  - The eight-covariate adjustment has no Type I check of its own beyond MOTCO's.
- README and AGENTS.md document the command, the constants and the reading rules.

Out of scope, for later changes:

- a revised primary analysis, since MOTCO's `candidate_fails` rules latent adjustment out as a replacement;
- `group_stage` conditioning (not admissible);
- other covariate sets;
- latent adjustment on a preprocessing variant's export;
- block weighting and the rank-1 `angle` degeneracy (MOTCO roadmap items 4–5).

## Capabilities

### New Capabilities
- `motco-latent-adjustment-sensitivity`: the three-way covariate-set and latent-basis sensitivity analysis of the astrocyte trajectory comparison through MOTCO's `stage`-conditioned latent adjustment. It covers the predeclared covariate set and arms, the primary as a read-only reference, the predeclared mediator check, the generated comparison table, and a committed report that reads the arms as sensitivity only.

### Modified Capabilities
<!-- None. `motco-trajectory-analysis` requires a pinned MOTCO ≥ 0.8.0, which 0.9.0 satisfies; its outputs are reproduced unchanged under the new pin (a task, not a requirement change). `latent-axis-diagnostic` already reads any analysis root given as `--analysis-dir`. -->

## Impact

- `src/sead_sexdiff/data.py`:
  - new constants `LATENT_SENSITIVITY_COVARIATES`, `LATENT_SENSITIVITY_ARMS`, `LATENT_SENSITIVITY_CONDITIONING` (`"stage"`), `LATENT_MOTCO_MIN_VERSION` (`0.9.0`) and `LATENT_SENSITIVITY_*` output names;
  - `_analyze_cohort` and `load_analysis_inputs` gain the covariate set and the latent-adjustment setting as keyword arguments. Their defaults reproduce the current behaviour byte-identically;
  - a `run_latent_sensitivity` function and a `render_latent_sensitivity` renderer;
  - `render_analysis_summary` prints the covariate set and, when present, the latent-adjustment block, without changing its output for the primary.
- `src/sead_sexdiff/cli.py`: new `latent-sensitivity` subcommand (`--cohort`, `--input-dir`, `--analysis-dir`, `--pseudobulk-dir`, `--out-dir`, `--n-jobs`, `--force`). The arms and the covariate set are constants, not flags.
- `pyproject.toml` / `uv.lock`: MOTCO `v0.8.0` → `v0.9.0`. MOTCO pins scikit-learn exactly, so if the lock needs to move it, it moves with MOTCO and is not loosened.
- `reports/`: one new committed report. `results/` stays gitignored.
- Compute: two CV analyses per cohort, about a minute each like `analyze`. The latent-adjusted CV adds one OLS per training set. Runs on the workstation holding `results/motco_input/`, `results/motco_analysis/` and `results/pseudobulk/` (incl. the technical table). This machine has none of them.
- Scientific choices introduced here (flag and ask before changing): the covariate set, the arms, `stage` conditioning, the mediator check and its threshold reuse, the raw fraction scale, `n_jobs = 1` for reported results, and the rule that the primary analysis stays primary.
