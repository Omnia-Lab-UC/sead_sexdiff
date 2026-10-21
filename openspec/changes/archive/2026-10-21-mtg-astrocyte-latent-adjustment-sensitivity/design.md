# Design

## Context

See `proposal.md` for the motivation. The state that shapes the approach:

- **`_analyze_cohort`** (shared by `analyze`, `rank-sensitivity`, `preprocessing-sensitivity`) builds `TrajectoryDataset(..., covariate_cols=ANALYSIS_COVARIATES, missing_covariates=ANALYSIS_MISSING_COVARIATES)` with the covariate list **hardcoded**, prints it, and writes the manifest's `covariates` block and `cli_equivalent` from the same constant. `load_analysis_inputs` builds the metadata table as `_ANALYSIS_DESIGN_COLUMNS + ANALYSIS_COVARIATES` from the export's `covariates.csv`, which holds no technical column.
- **The technical table** (`donor_technical.csv`) has `ribosomal_pseudobulk` and `mitochondrial_pseudobulk` per donor. `_check_technical_reference(pseudobulk_dir, donors=, columns=)` already checks the pseudobulk manifest hash, `TECHNICAL_GENE_SETS` and completeness. The `residualized` variant uses these columns on the raw fraction scale.
- **`_check_reference(cohort_dir, analysis_dir, inputs, motco_version)`** checks the five export hashes, MOTCO version, seed and permutations of the primary analysis. With the pin moved to 0.9.0, it fails against the 0.8.0 primary until `analyze` is rerun.
- **MOTCO 0.9.0** adds `TrajectoryAnalysisParams.covariate_latent_adjustment` (`None` | `"stage"` | `"group_stage"`). The setting without covariates is an error. The residualization covariates are exactly `covariate_cols`. The record gains a `latent_adjustment` block (`conditioning`, `covariates`, `design_columns`, `n_samples`, `residual_df`, `cv_min_training_df`, `feature_scale`, `feature_r2`) and `params.covariate_latent_adjustment`. Without the setting, every output equals 0.8.0's except the runtime (`tests/test_analysis_identity.py`). Under `stage` the CV folds equal the unadjusted analysis's folds. The rare-level fold failure needs a categorical covariate, and all eight covariates here are continuous.
- **`axis-diagnostics`** takes `--input-dir`, `--analysis-dir`, `--rank-dir`, `--pseudobulk-dir`, `--out-dir`. `_check_analysis_inputs` checks every entry of `manifest["inputs"]["files"]` against `<input-dir>/<name>`. Its partial correlations condition on the six `ANALYSIS_COVARIATES`, which stays the diagnostic's definition for the arms too.
- **Data location**: the real inputs are on the workstation that ran `analyze` and the sensitivity steps. This machine has none of them.

## Goals / Non-Goals

**Goals:**
- One command that produces both arms and the comparison per cohort, with the provenance discipline of the other sensitivity steps.
- An arm directory is a primary directory with exactly the covariate set (and, for `latent_adjusted`, the latent setting) changed.
- `analyze`, `rank-sensitivity` and `preprocessing-sensitivity` keep byte-identical outputs across the refactor.
- The mediator check is on disk before any arm's result exists.

**Non-Goals:**
- Any adoption rule or revised-primary path. MOTCO's `candidate_fails` excludes it, so `ADOPTION_RULE` is not evaluated here.
- Latent adjustment of preprocessing-variant exports, `group_stage`, forced ranks with latent adjustment.
- Rerunning the 0.8.0 sensitivity outputs.

## Decisions

**D1. Parameterize the covariate set and the latent setting, with defaults equal to today's behaviour.**
- `load_analysis_inputs(cohort_dir, *, extra_covariates=None)` merges `extra_covariates` (a donors × columns DataFrame) by `donor_id` after the log10 derivation, and returns `metadata` with `_ANALYSIS_DESIGN_COLUMNS + covariates`. `None` gives the current table.
- `_analyze_cohort(..., covariates=ANALYSIS_COVARIATES, ...)` uses `covariates` for the dataset, the progress line, the manifest's covariate block and `cli_equivalent`. The latent setting travels in `params` (`params.covariate_latent_adjustment`). `cli_equivalent` then adds `--latent-adjustment stage`, and the manifest's `measurement` gains `latent_adjustment`. Both are absent, not `null`, when the setting is `None`, as MOTCO leaves `covariate_latent_adjustment` out of its record: a `null` key in the primary's manifest would break the JSON equality of the byte-identity check below and of D7. Two further keywords carry the arm's provenance and are likewise absent by default: `covariate_sources` (recorded as `covariates.sources`) and `extra_inputs` (merged into `inputs` beside `files`, for `inputs.technical`).
- `render_analysis_summary` prints the covariate list from the record, which it may already do, and a latent-adjustment line when `record` has the block.
- Verification: on the synthetic dataset, outputs of `analyze` and `rank-sensitivity` before and after the refactor are byte-identical (CSVs and `summary.md`), and the JSON files are equal apart from `runtime_seconds`, `created_at` and `git`. This repeats the rank-sensitivity change's D1 check.
- *Alternative: a separate analysis function for the arms.* Rejected: it would duplicate the manifest and summary code, and the comparison rests on nothing else differing.

**D2. Registry and constants.**
```python
LATENT_SENSITIVITY_ADDED_COVARIATES = ("ribosomal_pseudobulk", "mitochondrial_pseudobulk")
LATENT_SENSITIVITY_COVARIATES = ANALYSIS_COVARIATES + LATENT_SENSITIVITY_ADDED_COVARIATES
LATENT_SENSITIVITY_CONDITIONING = "stage"
LATENT_SENSITIVITY_ARMS = {
    "extended_covariates": {"description": ..., "covariates": LATENT_SENSITIVITY_COVARIATES, "latent_adjustment": None},
    "latent_adjusted": {"description": ..., "covariates": LATENT_SENSITIVITY_COVARIATES, "latent_adjustment": LATENT_SENSITIVITY_CONDITIONING},
}
LATENT_MOTCO_MIN_VERSION = "0.9.0"
```
The added fractions are exactly `RESIDUALIZATION_COVARIATES` minus the export's `log10_lib_<omic>`. An assertion at import ties them, so a later change to the residualization set cannot silently diverge from this one. `check_motco_version` gains a `required=` keyword. `analyze` keeps its 0.8.0 minimum, and `latent-sensitivity` requires 0.9.0.
- *Why not `log10_lib_<omic>`:* the nuclei-count covariates are already in. Ten covariates would go beyond the eight that MOTCO's Type I check used, at cells of 9–28 donors. And the library sizes are export-derived, which would make the arm depend on the export's step 2.
- *Why raw fractions:* the `residualized` variant used this scale, which keeps the two sensitivity analyses comparable. A logit would be a new choice made with the results already partly seen.

**D3. Mediator check first.** `mediator_check(technical, metadata, columns)` returns rows `cohort, covariate, subset (all | F | M), n, spearman, flagged` with `scipy.stats.spearmanr` on stage coded 0/1/2, where `flagged = |ρ| ≥ DIAGNOSTIC_CORRELATION_THRESHOLD`. It is written to `<cohort>.tmp/mediator_check.csv` before the first arm runs, and the flags are copied into the manifest's `mediator_flags`. The flags change nothing downstream.
- *Alternative: partial on the six covariates.* Rejected: the question is whether the fraction moves with stage at all. Conditioning on RIN/PMI could hide exactly the pathway a mediator would take.
- *Alternative: a gate (drop flagged fractions).* Rejected: MOTCO's limitation says no data can decide confounder vs mediator, and a gate would choose the covariate set from a stage association.

**D4. Output layout and atomicity.**
```
results/motco_latent_sensitivity/
  extended_covariates/motco_analysis/<cohort>/   # full analysis directory
  latent_adjusted/motco_analysis/<cohort>/
  comparison/<cohort>/{latent_sensitivity.csv,.md,_manifest.json, mediator_check.csv}
```
Per cohort, both arms are written via `_analyze_cohort(out_dir=<arm>/motco_analysis, final_name=f"{cohort}.tmp")` into temporary directories. Once the comparison is written, all three are renamed. If anything fails, the temporaries are removed and the finished directories of other cohorts are untouched. The arm manifest's `inputs.files` holds only the five export files, and the technical table goes under `inputs.technical`. This lets `axis-diagnostics --input-dir results/motco_input --analysis-dir results/motco_latent_sensitivity/<arm>/motco_analysis` pass its hash check without changes. Existing outputs are an error without `--force`.
- *Alternative: one directory per cohort containing both arms.* Rejected: `axis-diagnostics` and `rank-sensitivity` expect `<root>/<cohort>/`.

**D5. Comparison table.** Reuse `_sensitivity_row(record, cohort=, component_selection=, directory=)` and extend each row with `arm`, `covariates` (`;`-joined), `conditioning`, `residual_df`, `cv_min_training_df` and `feature_r2_{q05,q25,q50,q75,q95,mean,frac_gt_0_5}` (empty for the reference and `extended_covariates`). Floats are written with `repr`, and `runtime_seconds` stays out of the CSV (as in `rank_sensitivity.csv`). `render_latent_sensitivity` writes one table per statistic (reference / arm × observed, p, rejected, null mean / SD / q95), a rank and AUROC table, the latent-adjustment block, and the mediator check with flagged rows first.

**D6. CLI.** `sead-sexdiff latent-sensitivity [--cohort {cohort_30,cohort_50,all}] [--input-dir results/motco_input] [--analysis-dir results/motco_analysis] [--pseudobulk-dir results/pseudobulk] [--out-dir results/motco_latent_sensitivity] [--n-jobs 1] [--force]`, printing `[latent-sensitivity]` progress per arm and the comparison at the end.

**D7. Pin bump procedure.** `uv add "motco @ git+https://github.com/tomszar/MOTCO@v0.9.0"`. If `uv.lock` moves scikit-learn, that is MOTCO's pin and is accepted, never loosened.
- On the workstation, copy `results/motco_analysis/` to `results/motco_analysis_0.8.0/`, then run `analyze --force` for both cohorts under 0.9.0. On the same machine, run the pre-change code (a `git archive` of the commit before the change, `uv sync --locked`, MOTCO 0.8.0) on the same export into a scratch directory.
- A scratch script compares the two reruns: the four CSVs byte-for-byte, and `result.json` / `analysis_manifest.json` after dropping `runtime_seconds`, `created_at`, `git`, `software.motco` and the `cli_equivalent` output path. This is the gate: if it fails, the copy is restored and the change stops.
- The same script compares the 0.9.0 rerun with the kept directories, and the report records the outcome. It is not a gate, because floating-point results are not byte-reproducible across machines. Decision of 2026-10-21 (user), on this evidence from `nixos-office`: against the kept primary, made on the laptop, the rerun differed only at rounding level (max abs diff 2.6e-13 on latent scores up to 87, 1e-16 on loadings, statistics about 1e-15 relative; identical p-values, ranks and `analysis_metadata.csv`), plus five `preprocessing.*` manifest keys that later archived changes added. MOTCO 0.8.0 against 0.9.0 on the same machine was byte-identical.
- If the gate passes, the 0.9.0 rerun stays as `results/motco_analysis/`, the reference of the arms; its `result.json` hashes are new and reported. The kept copy stays as `results/motco_analysis_0.8.0/` until the report is committed. The cluster's `results/motco_analysis/` is replaced by the rerun, and its old one is kept as `results/motco_analysis_0.8.0/`, so the hashes the earlier reports cite stay resolvable.

**D8. Report.** `reports/motco-astrocyte-latent-adjustment-<date>.md` is assembled by a scratch script. It pastes each cohort's `latent_sensitivity.md` and `mediator_check.csv`, and the `pls_0` rows of `axis-diagnostics` on the primary and on both arms. The arm diagnostics go to `results/motco_latent_sensitivity/diagnostics/<arm>/`, and the rank directory is absent there, which that command records as a finding. Provenance and the reading rules are hand-written. A check script finds every statistic, p-value and rank from the six records in the text.

## Risks / Trade-offs

- [The fractions are mediators: pathology lowers astrocyte RNA quality] → the mediator check is written before the arms, the report carries it, and `extended_covariates` shows how much of any change comes from the covariates themselves, before the basis changes.
- [Eight covariates at cells of 9–28 donors make the model near-singular] → MOTCO errors on a rank-deficient design and records `residual_df` and `cv_min_training_df`, and these are in the table. If an arm fails, the cohort fails without partial output and the error is reported, not worked around.
- [The `latent_adjusted` arm rejects something and is read as the result] → the spec fixes the report wording, the primary is the reference row everywhere, and MOTCO's `candidate_fails` is quoted.
- [The refactor changes the primary or the other steps' outputs] → D1's synthetic byte-identity check, then D7's real-data check.
- [0.9.0 changes the primary's numbers] → D7's same-machine comparison detects it, restores, stops. A cross-machine rounding difference is recorded, not mistaken for one.
- [Collinearity of the fractions with `log10_n_*` or RIN] → harmless for the statistics (adjustment is joint). The `feature_r2` summary and the arm-vs-arm comparison show how much the set removes.

## Migration Plan

1. Pin bump and refactor (code, synthetic checks) on any machine.
2. On the workstation, sync the branch, `uv sync --locked`, D7 check (same-machine reruns under 0.9.0 and 0.8.0, and the comparison with the kept directories), then `latent-sensitivity --cohort all`, then `axis-diagnostics` on both arms.
3. Report, check, commit.

Rollback: revert the pin, and restore `results/motco_analysis_0.8.0/` (on the workstation and on the cluster).
