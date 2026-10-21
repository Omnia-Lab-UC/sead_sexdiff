# Tasks

## 1. Constants and synthetic baseline

- [ ] 1.1 Add `STAGE_DE_MODELS`, `STAGE_DE_COVARIATES`, `STAGE_DE_ALPHA`, `STAGE_DE_REFERENCE_LEVELS`, `STAGE_DE_FEATURE_STATUSES`, `STAGE_DE_CONCORDANCE_RULE` and the `STAGE_DE_*` output names to `data.py` (design D4–D7). Verify by importing `sead_sexdiff.data` and printing the registry: four models, the tested coefficients named per design D3, `rna` and `atac` covariate tuples of five.
- [ ] 1.2 Rebuild the synthetic pseudobulk + export + analysis dataset in the scratchpad (the builder of the latent-adjustment change, with a planted stage effect on 50 genes / 50 peaks and a planted sex × stage interaction on 10 genes / 10 peaks, all in the abundance-passing set) and run `export` and `analyze` on it. Verify that `features_<omic>.csv` carries all four statuses and that the planted features have status `kept` or `variance`.

## 2. Fitting

- [ ] 2.1 Implement the per-fit helper: build the design table (sex and stage as categoricals with the reference levels, `stage_num`, the omic's covariates as continuous), fit `DeseqDataSet` on the export-bound feature set, resolve each tested coefficient's design column by name, run `DeseqStats` with a unit contrast vector, return the result table with `exported`. Verify on the synthetic data that the resolved column names are the ones design D3 expects and that an unresolvable name raises an error naming it.
- [ ] 2.2 Verify the fits on the synthetic data: the planted stage genes dominate the smallest `padj` of `factor_additive` and `trend_additive`, the planted interaction genes dominate `factor_interaction` and `trend_interaction`, and unplanted features show no enrichment at `padj ≤ 0.05` beyond the nominal rate.
- [ ] 2.3 Implement `run_stage_de(cohort, omic, pseudobulk_dir, input_dir, analysis_dir, out_dir, n_cpus, force)`: input existence checks with the producing command named, the pseudobulk ↔ export and export ↔ analysis reference checks (design D8), the cohort donors in the export's order, one table per model × coefficient, resumable per table, the summary with both feature sets and the union rows, the concordance cells, the manifest with hashes and versions. Verify on the synthetic data that a modified `rna_counts.h5ad` stops the step naming the file and that nothing is written on failure.
- [ ] 2.4 Implement `render_stage_de` (`stage_de.md` from summary + manifest) and the `stage-de` CLI subcommand with `--cohort`, `--omic`, `--pseudobulk-dir`, `--input-dir`, `--analysis-dir`, `--out-dir`, `--n-cpus`, `--force`. Verify with `uv run sead-sexdiff stage-de --help` and a synthetic run that `stage_de.md` lists the models, covariates, counts and concordance table.

## 3. Real data

- [ ] 3.1 Run `uv run sead-sexdiff stage-de --cohort all --omic both` on the workstation holding `results/pseudobulk/`, `results/motco_input/` and `results/motco_analysis/` (the 0.9.0 primary). Verify that both cohorts' directories hold twelve result tables, the summary and the manifest, that every exported `n_tested` is 5,000 and every abundance `n_tested` equals the export manifest's post-abundance count, and that the sha256 of every input file is unchanged.
- [ ] 3.2 Rerun `cohort_30` with `--force` into a scratch directory. Verify that `stage_de_summary.csv` is byte-identical and every table's `padj` equal.
- [ ] 3.3 Rerun one `cohort_30` RNA model with a different `--n-cpus` into the scratchpad. Verify that `pvalue` and `padj` are equal to 3.1's; if they are not, record the difference in the manifest's findings and in the report, and report results at `n_cpus = 1`.
- [ ] 3.4 Cross-check one table with R DESeq2 (the same design, counts and reference levels; add `DESeq2` to `flake.nix` if absent): the `factor_interaction` RNA fit of `cohort_30`. Verify that log2 fold changes and Wald statistics agree to 1e-3 relative for the tested coefficients and that the set at `padj ≤ 0.05` is the same; record any difference in the report.

## 4. Report and docs

- [ ] 4.1 Write `reports/motco-astrocyte-stage-de-<date>.md` from the generated files: provenance (inputs with hashes, pydeseq2 version, git), the count tables per cohort × omic on both feature sets, the concordance table, the top features per coefficient (descriptive), the reading rules and limitations the spec requires. Verify with a check script that every count and concordance cell equals `stage_de_summary.csv` and that the `delta` statuses equal the primary `result.json` records.
- [ ] 4.2 Update README (command, flags, models, covariates, feature sets, outputs, reading rules) and AGENTS.md (pipeline step 11, the constants under "flag and ask", the CLI name). Verify that every `STAGE_DE_*` name in both files exists in `data.py`.
- [ ] 4.3 Run `openspec validate mtg-astrocyte-stage-de` and verify that it passes.
