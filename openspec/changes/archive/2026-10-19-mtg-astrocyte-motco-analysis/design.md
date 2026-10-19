# Design

## Context

See `proposal.md` for the motivation. The state that shapes the approach:

- **Export outputs** (`results/motco_input/<cohort>/`, laptop): `rna.csv` and `atac.csv` (donors × 5,000 log2-CPM features each, first column `donor_id`), `metadata.csv` (`donor_id`, `sex`, `stage`, `stage_label`), `covariates.csv` (`donor_id`, `age_at_death`, `age_capped`, `pmi`, `rin`, `years_education`, `education_level`, `n_rna`, `n_atac`) and `manifest.json` (thresholds, feature counts per step, TMM factors, annotation, findings). `cohort_30` has 80 donors, `cohort_50` 74; `k = 9` in both. Ages are uncapped and no covariate is missing in this release.
- **MOTCO 0.8.0** (`../MOTCO`, public at `github.com/tomszar/MOTCO`, tag `v0.8.0`, Python ≥ 3.11, pins `scikit-learn == 1.8.0`). `motco analyze` takes blocks as `--block name=path` CSVs and one metadata CSV that must already hold the covariate columns (`--covariate COL`, repeatable). The Python entry point `motco.analysis.analyze_trajectories(TrajectoryDataset, TrajectoryAnalysisParams)` takes the same inputs as DataFrames; `write_analysis_outputs(result, out_dir, overwrite=False)` writes `result.json`, `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv` and refuses an existing `result.json` unless `overwrite`. The record's `covariates` block has `columns`, `kinds`, `design_columns`, `reference_row`, `missing_policy`, `dropped_sample_ids`; `latent` has `selected_lv`, `component_selection`, `cv_mean_auroc`; `measurement` has the estimator, permutations and seed. `TrajectoryAnalysisParams` defaults: `orientation_estimator="net_displacement"`, `cv1_splits=3`, `cv2_splits=4`, `n_repeats=5`, `max_components=20`, `random_state=1203`, `permutations=0`, `seed=None`, `n_jobs=1`.
- **MOTCO guards that matter here:** a constant covariate (`age_capped`, all `False`) is rejected; a non-finite covariate is rejected; attribution with covariates is rejected.
- **Binding interpretation** (`../MOTCO/docs/roadmap.md` Phase 6 and Next changes item 1; `reports/net-displacement-power-2026-10-14.md` §7; `reports/phase6-small-n-2026-10-09.md` §8; `reports/covariate-type-i-2026-10-16.md` §4–5): CV-selected rank beside `angle`; `angle`/`shape` non-rejections uninformative; `delta` rejection is evidence of a difference, not of size; feature filtering reported beside results; covariate block reported; additive-linear, covariate-blind-basis and structural-analogue limitations carried; no power claim for the adjusted analysis.
- **Existing helpers in `data.py`:** `_git_revision`, `_software_versions`, `EXPORT_COHORTS`, `STAGE_LABELS`, the sha256 helper used for the annotation and input-manifest hashes, and the `.tmp/` + rename pattern from `run_export`.

## Goals / Non-Goals

**Goals:**
- One command per cohort, from the export outputs, to a MOTCO artifact plus a manifest that makes the run reproducible by hand (`motco analyze` on `analysis_metadata.csv`).
- Every scientific setting is a named constant in `data.py`, recorded in the manifest, and never a CLI flag.
- The numbers in the committed report come from a generated `summary.md`, not from transcription.

**Non-Goals:**
- Any change to `export` or its manifest. The log10 columns are derived here.
- Forced-rank reruns, unadjusted runs, attribution, block weighting, enrichment.
- A pytest suite; verification is synthetic checks plus the real run, as before.

## Decisions

**D1. Call the MOTCO Python API, and also write the CLI-equivalent inputs.**
`run_analysis` builds `TrajectoryDataset` and calls `analyze_trajectories` + `write_analysis_outputs` directly. The roadmap allows either entry point; the API avoids a subprocess, a temporary merged CSV and parsing MOTCO's stdout, and lets this repo hash the inputs it actually passed. For reproduction by hand, the merged table is still written as `analysis_metadata.csv` and the manifest records the equivalent `motco analyze` command line (`cli_equivalent`), including the six `--covariate` flags and `--stage-order 0,1,2`.
- *Alternative: subprocess `motco analyze`.* Rejected: it needs the merged CSV anyway, and the record would be identical.

**D2. Covariate table assembled at analysis time; constants own the set.**
`ANALYSIS_COVARIATES = ("age_at_death", "pmi", "rin", "years_education", "log10_n_rna", "log10_n_atac")` and `ANALYSIS_LOG10_COLUMNS = {"log10_n_rna": "n_rna", "log10_n_atac": "n_atac"}`. The merge is `metadata.csv ⋈ covariates.csv` on `donor_id` (inner join, then an equality check on the donor sets, error on any difference), the log10 columns are added, and only `donor_id`, `sex`, `stage`, `stage_label` and `ANALYSIS_COVARIATES` are kept in `analysis_metadata.csv`. `education_level` and `age_capped` never reach MOTCO. The manifest lists the exported columns that were not entered and why (`education_level`: by decision; `age_capped`: constant in this release).
- *Alternative: have `export` write the log10 columns.* Rejected: it would change an archived spec for a derivation that belongs to the model.
- *Alternative: enter `age_capped` when it varies.* Rejected for this change: the set is fixed; a future release with capped ages should prompt a decision, which the manifest's `n_age_capped` (copied from the export manifest) makes visible.

**D3. Measurement settings.**
`ANALYSIS_PERMUTATIONS = 999`, `ANALYSIS_SEED = 20261019`, `ANALYSIS_ALPHA = 0.05` with rejection when p ≤ α (`ANALYSIS_REJECTION_RULE`; decided during apply, 2026-10-19: it is the rule of MOTCO's covariate Type I check, which validated the adjusted analysis, and gives exact size α for (1 + b)/(1 + B) p-values; MOTCO's earlier operating studies count p < 0.05, which differs only at p = 0.050), MOTCO CV defaults untouched, `forced_components=None`, attribution off, `missing_covariates="error"`. Both cohorts use the same seed. `--n-jobs` is the only runtime flag (default 1); it must not change the result, and task 3.3 checks that `n_jobs=1` and `n_jobs=4` give identical statistics on `cohort_50`. *Outcome (task 3.3, 2026-10-19):* the statistics, rank, latent scores, loadings and stage means are identical, but the p-values and null summaries are not. MOTCO 0.8.0's parallel RRPP seeds one stream per worker (`SeedSequence(seed).spawn(n_workers)`), whereas the serial path uses one stream. Results are therefore reproducible at a fixed `n_jobs` only. The manifest records `n_jobs` with a note and adds a finding when `n_jobs != 1`, and the report uses `n_jobs = 1`.

**D4. Version gate before reading inputs.**
`run_analysis` compares `motco.__version__` with `ANALYSIS_MOTCO_MIN_VERSION = "0.8.0"` (`packaging.version` if available, else a tuple compare) and raises before touching the files. The dependency is pinned in `pyproject.toml` as `motco @ git+https://github.com/tomszar/MOTCO@v0.8.0` through `uv add`, which locks the commit hash in `uv.lock`.
- *Alternative: PyPI.* MOTCO is not published there.
- *Alternative: a path dependency on `../MOTCO`.* Rejected: machines differ, and the lock would not pin a revision.

**D5. Output layout and write discipline.**
`results/motco_analysis/<cohort>/` is written to `<cohort>.tmp/` and renamed when complete, as `run_export` does; an existing cohort directory is an error without `--force`. `analysis_manifest.json` is written last and holds: `cohort`, `inputs` (path and sha256 per file, plus the export manifest's own hash), `preprocessing` (copied from the export manifest's `thresholds` and `feature_counts`, plus `n_age_capped` and its findings), `covariates` (`entered`, `derived`, `excluded`, `missing_policy`), `measurement` (estimator, `component_selection`, permutations, seed, alpha, `n_jobs`), `cli_equivalent`, `software` (existing `_software_versions` plus `motco` and `scikit-learn`), `git`, `findings`. `summary.md` is rendered from `result.json` with a small formatter (`_render_summary`) that emits the statistics table (value, p, null mean/sd/q95, with the selected rank in the `angle` row), the design cells, the latent block and the covariate block.

**D6. The committed report is hand-written from the two `summary.md` files, with the binding statements.**
`reports/motco-astrocyte-<date>.md` follows the layout of MOTCO's reports: run provenance (MOTCO version and commit, seed, git revision of this repo), one results section per cohort pasting the generated tables, a "Feature filtering" section copied from the export manifests, then "Interpretation" that applies the binding rules verbatim and "Limitations". The tasks require a check that every number in the report equals the record (a short script comparing the Markdown numbers with `result.json`, run once before commit).

## Risks / Trade-offs

- [`uv add` cannot resolve `scikit-learn == 1.8.0` against scanpy/decoupler/pydeseq2] → the lock currently shows no upper bound on scikit-learn from any dependent, so 1.8.0 should resolve; if it does not, stop and report rather than loosening MOTCO's pin.
- [Six covariate columns at cells of 9–28 donors make the model near-singular or the null unstable] → MOTCO rejects a rank-deficient design with an error; the report states `design_columns` and the null sd per statistic. The covariate Type I check covered three design columns at this n and does not bound six; the report says so.
- [PLS spends a component on a covariate aligned with stage] → recorded limitation (covariate-blind basis); the residualization sensitivity analysis is a follow-up.
- [A `delta` rejection gets read as a size difference, or a non-rejection as agreement] → the spec fixes the wording; the report's interpretation section is a required artifact.
- [Runtime] → n = 80 and 10,000 features: the CV is 60 PLS fits of ≤ 20 components and RRPP is 999 OLS refits in a ≤ 20-dimensional latent space; minutes on the workstation.
- [`n_jobs` changes the result] → checked once (D3); if it does, the manifest records `n_jobs` and the report uses `n_jobs=1`.
