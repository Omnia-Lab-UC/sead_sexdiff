# Proposal

## Why

The MOTCO inputs exist (`results/motco_input/cohort_{30,50}/`, archived change `2026-10-16-mtg-astrocyte-motco-export`) and the MOTCO side is unblocked: MOTCO 0.8.0 (tag `v0.8.0`, 2026-10-16) adds covariates to `motco analyze` and its Type I check returned `covariates_validated` (`../MOTCO/docs/reports/covariate-type-i-2026-10-16.md` §5). Nothing in this repo runs the trajectory analysis yet, so the case study has no result. This change runs it, in the way the MOTCO roadmap binds (`../MOTCO/docs/roadmap.md`, Next changes item 1): through `motco analyze` / `motco.analysis.analyze_trajectories` on a pinned release, never through `motco plsr` + `motco de`.

## What Changes

- **Pin MOTCO.** Add `motco` (≥ 0.8.0, from the public GitHub tag) as a dependency with `uv add`, and check that it installs on this project's Python ≥ 3.12 lock. MOTCO pins `scikit-learn == 1.8.0`; the lock currently holds 1.7.2, so the resolution has to be verified.
- **New `analyze` step.** For each cohort (`cohort_30` primary, `cohort_50` sensitivity) read only `results/motco_input/<cohort>/` and run the covariate-adjusted trajectory analysis once:
  - blocks `rna` and `atac` in identity units; group `sex`; stage `stage` with the explicit order `0, 1, 2`;
  - the MOTCO defaults the operating studies ran (net-displacement orientation, double-CV rank with `cv2_splits=4`, `cv1_splits=3`, `n_repeats=5`, `max_components=20`), 999 RRPP permutations and a fixed seed;
  - covariates, additive and in this order: `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`. The two nuclei counts enter on the log10 scale (derived at analysis time and recorded). `education_level` is not entered (the user's decision: it is categorical and would spend one column per level). `age_capped` is not entered because it is constant in this release (0 capped ages) and MOTCO rejects constant covariates.
  - missing-covariate policy `error`: no value is missing in this release, and nothing is imputed.
  - There is no unadjusted run. The adjusted run is the analysis.
- **Analysis outputs per cohort** (`results/motco_analysis/<cohort>/`): MOTCO's artifact (`result.json`, `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv`), the merged metadata table the model saw (`analysis_metadata.csv`, so the equivalent `motco analyze` command can be rerun by hand), an `analysis_manifest.json` (input file hashes, the export manifest hash, the covariate set and derivation, the preprocessing summary copied from the export manifest, MOTCO and software versions, git revision, findings) and a `summary.md` table generated from the record.
- **A committed report** (`reports/motco-astrocyte-<date>.md`) with, for each cohort: `delta`, `angle` and `shape` with their p-values and null summaries; the CV-selected rank beside `angle`; the covariate block from the record (design columns, missing policy, dropped donors); the design cells; the feature filtering beside the results; and the binding interpretation statements and limitations carried from the MOTCO reports. `results/` stays gitignored.
- README and AGENTS.md document the command, the constants and the interpretation rules.

Out of scope, as follow-up changes: the forced-rank (dimensionality) sensitivity sweep; preprocessing sensitivity exports (other N, keeping chrX); pathway/module enrichment; attribution (MOTCO refuses it with covariates); block weighting (MOTCO roadmap item 3).

## Capabilities

### New Capabilities
- `motco-trajectory-analysis`: the covariate-adjusted sex × ADNC-stage trajectory comparison of the astrocyte cohorts through MOTCO's analysis entry point, with a fixed, recorded covariate set, reproducible outputs, and a committed report that follows the binding interpretation rules.

### Modified Capabilities
<!-- None. motco-input-export keeps "covariates exported, not applied"; the log10 derivation happens in the analysis step and is recorded there. -->

## Impact

- `pyproject.toml`, `uv.lock`: new dependency `motco` from `https://github.com/tomszar/MOTCO` at `v0.8.0` (public, HTTPS, so the cluster login node can also sync). Likely upgrades `scikit-learn` to 1.8.0; `scipy`, `matplotlib`, `tqdm` are already in the lock.
- `src/sead_sexdiff/data.py`: new constants (`ANALYSIS_COVARIATES`, `ANALYSIS_LOG10_COLUMNS`, `ANALYSIS_PERMUTATIONS`, `ANALYSIS_SEED`, `ANALYSIS_ALPHA`) and `run_analysis`. Export, pseudobulk, covariates, subset and DE functions are untouched.
- `src/sead_sexdiff/cli.py`: new `analyze` subcommand.
- `reports/`: new directory, committed.
- Compute: runs on any workstation holding `results/motco_input/` (n = 80, 10,000 features, 999 permutations; minutes, not hours). No cluster job.
- Scientific choices introduced here (flag and ask before changing): the covariate set and order, the log10 scale for nuclei counts, 999 permutations, the seed, α = 0.05.
