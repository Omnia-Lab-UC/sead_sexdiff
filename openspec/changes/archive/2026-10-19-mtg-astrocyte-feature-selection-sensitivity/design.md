# Design

## Context

See `proposal.md` for the motivation. The state that shapes the approach:

- **The export's per-block pipeline** is `preprocess_block(counts, names, chromosomes, *, k, min_cpm, prior_count, top_n, excluded_chromosomes)`: (2) chromosome removal with status `sex_chromosome`, (3) abundance with status `abundance`, (4) TMM, (5) log2-CPM, (6) top `N` by population variance with ties by name and status `variance`; it returns `values`, `kept`, `status`, `variance`, `steps`, `findings` and the filtered intermediates. `_export_cohort` calls it per omic, writes the block CSV with `EXPORT_FLOAT_FORMAT`, `_features_frame` (feature, gene id, chromosome, status, variance), `metadata.csv`, `covariates.csv` and the manifest (`preprocessing_order`, `thresholds`, `feature_counts`, `tmm`, an optional `variant` block). `run_export` output must stay byte-identical (the archived sensitivity change checks this on synthetic data and reproduced the kept export on `nixos-home`).
- **The sensitivity runner** (`run_preprocessing_sensitivity`) checks, per variant and cohort, `_check_export_reference` and `_check_reference`, then calls `_export_cohort(**variant["export"])` and `_analyze_cohort` with the primary parameters, and writes `comparison/<cohort>/` through `_variant_row`, `_feature_accounting(primary_export_dir, variant_export_dir)` (chrX-specific: kept chrX counts and names, shared / only counts, a `sex_chromosome` check) and `render_preprocessing_sensitivity`. The CLI's `--variant` choices are the registry keys plus `all`.
- **The technical table** (`donor_technical.csv`, change `mtg-astrocyte-technical-axis-diagnostic`) holds per donor the gene-set pseudobulk fractions (`ribosomal_pseudobulk`, `mitochondrial_pseudobulk`, `hemoglobin_pseudobulk`), nucleus QC means and ATAC library metrics, with `technical_manifest.json` recording the pseudobulk manifest hash. `TECHNICAL_GENE_SETS` and the gene-set matcher, the Spearman / partial-correlation functions and `DIAGNOSTIC_CORRELATION_THRESHOLD` are that change's. The technical-axis report classifies each gene-set fraction as preparation property or cell state under a fixed rule.
- **What residualization does to the primary's checks.** `metadata.csv` and `covariates.csv` of a variant must stay byte-identical to the primary's (reference check), so the design matrix needs its own file. `load_analysis_inputs` reads only the five analysis inputs, so an extra `residualization.csv` in the variant export directory is ignored by the analysis and must be excluded from the five-file hash set but recorded in the export manifest.
- **The 2026-10-19 scratch results that fix expectations.** Appending the technical variables to X leaves the PLS axes unchanged; residualizing the exported blocks on the TMM factors, library sizes, nuclei counts, RIN and PMI makes the axes orthogonal to them, equalizes component spreads (17.6 and 15.8 against 35 and 17) and raises the first axis's stage correlation to 0.59, while a −0.50 correlation with the ribosomal score remains because the TMM factor is only a proxy; hence the gene-set fractions in the set and the TMM factors kept out of it as check variables.

## Goals / Non-Goals

**Goals:**
- Each new variant is one more registry entry whose `export` dict names the setting it changes; the export code grows three keyword parameters with primary-preserving defaults and no new code path for the primary.
- The adoption check is a generated file computed from constants, so the decision to propose a revised primary is mechanical and auditable.
- The report can attribute the change in the first axis to each setting because the single-setting variants and the candidate share every other setting.

**Non-Goals:**
- Choosing `N`, the trend method's parameters or the technical set from results.
- Changing the variant analysis: it stays the primary's (`ANALYSIS_COVARIATES`, CV rank selection, seed, permutations).
- Running the forced-rank sweep or the axis diagnostics on the variants as part of this step (the layout serves both; the report may cite an `axis-diagnostics` run on a variant as a descriptive check).

## Decisions

**D1. Three keyword parameters on `preprocess_block`, threaded from the registry.**
`preprocess_block(..., excluded_gene_sets=(), selection_rule="variance", residualize=None)`. `excluded_gene_sets` names sets of `TECHNICAL_GENE_SETS` matched on `names` (RNA only; ATAC names never match), removed at step 2 with status `excluded_gene_set`, counted in `removed_by_gene_set`; the step-3 library sizes exclude them, as they exclude chromosome-removed features. `selection_rule` is `variance`, `trend_ratio` or `mean`; the returned `variance` array is kept and a `selection_score` array is added (equal to variance under `variance`), with `trend` and `ratio` under `trend_ratio`; dropped status is `variance` under `variance` (byte-identity) and `selection` otherwise. `residualize` is `None` or a donors × covariates design `DataFrame` aligned to the counts' rows; when given, after step 5 each column of the log2-CPM matrix is replaced by `mean + residual` from an OLS fit on `[1, design]` (`np.linalg.lstsq` once for all features), the per-feature R² is kept, and step 6 ranks the residualized values. `_export_cohort` gets the same three keywords, builds the design from `donor_technical.csv` and the export manifest's library sizes for the cohort's donors, writes `residualization.csv`, and records `excluded_gene_sets`, `removed_by_gene_set`, `selection_rule` (per omic, with `trend_bins` edges and medians) and `residualization` (covariates, design sha256, R² quantiles per omic) in the manifest. With defaults the manifest records `excluded_gene_sets: []`, `selection_rule: {rna: variance, atac: variance}` and `residualization: null`.
- *Alternative: separate functions per variant.* Rejected: the claim "one setting changed" needs one code path.
- *Primary byte-identity:* the synthetic export B of the archived change (or a rebuilt one) before and after; the kept `results/motco_input/<cohort>/` reproduced on `nixos-home`. The manifest gains three keys, so manifest equality is checked apart from them; the task records this as the one accepted difference.

**D2. Per-omic selection rules as a dict in the registry.**
`"rna_trend_ratio": {"export": {"selection_rule": {"rna": "trend_ratio", "atac": "variance"}}}`, `"atac_mean_accessibility": {"export": {"selection_rule": {"rna": "variance", "atac": "mean"}}}`. `trend_ratio`: sort abundance-passing features by mean log2-CPM, cut into twenty equal-count bins (`np.array_split` of the sorted index), take each bin's median variance at its median mean, interpolate linearly (`np.interp`, flat beyond the end bins), rank by `variance / trend`, ties by name. `mean`: rank by mean log2-CPM, ties by name.
- *Alternative: loess (Seurat `vst`, scran `modelGeneVar`).* Rejected: a smoothing span is a tuning parameter and would need a new dependency; twenty equal-count bins with medians is deterministic, robust to the low-count tail and recordable in the manifest.
- *Alternative: a DESeq2 VST before raw-variance ranking.* Rejected here: it changes the transform the edgeR cross-check covered and MOTCO z-scores features anyway; it can be a later registry entry.

**D3. The residualization set and its cell-state exclusion rule.**
`RESIDUALIZATION_COVARIATES = ("log10_lib_rna", "log10_lib_atac", "ribosomal_pseudobulk", "mitochondrial_pseudobulk", "hemoglobin_pseudobulk")`, minus any fraction the technical-axis report classifies as cell state (its within-donor supertype share exceeding its between-donor share); `RESIDUALIZATION_CELL_STATE_EXCLUSIONS` records the removed names with the report's section as the reason, filled when this change is applied and empty if nothing is excluded. Library sizes are the export's own (post-chromosome totals, log10), so the design needs no analysis covariate. Fractions enter untransformed (they are bounded and small; the report shows their distribution). Nuclei counts, RIN and PMI stay out: they are analysis covariates already in the RRPP model and leaving them out keeps the double adjustment to the two library-size proxies.
- *Alternative: the TMM factors as regressors.* Rejected: they are normalization by-products and the proxy of the composition effect, not its measurement; they serve as check variables instead.
- *Alternative: also residualize against the six analysis covariates.* Rejected: the variant analysis must stay the primary's, and adjusting the same covariates twice in a different space makes the comparison harder to read; stated as a limitation for the two library-size proxies only.

**D4. The analysis is unchanged, residualization included.**
`_analyze_cohort` runs with the primary parameters on the variant blocks; `covariates.csv` is the primary's. The residualized blocks are orthogonal to the design by construction, so the regression on latent scores sees no library-size effect from them; `log10_n_rna` / `log10_n_atac` still enter the RRPP model as in the primary.

**D5. Reference checks gain the technical table.**
`_check_technical_reference(pseudobulk_dir)` reads `technical_manifest.json`, compares its recorded pseudobulk manifest sha256 with `results/pseudobulk/manifest.json`, and checks completeness of `RESIDUALIZATION_COVARIATES` (after exclusions) for the cohort's donors; run only for variants whose `export` dict declares `residualize` or `excluded_gene_sets`. The five-file hash set of the export stays as is; `residualization.csv`'s sha256 goes into the export manifest and the comparison manifest.

**D6. Variant-aware accounting and comparison.**
`_feature_accounting` keeps the shared / only counts and the chrX block (emitted when the variant's excluded chromosomes differ from the primary's), and adds: `removed_by_gene_set` when sets were excluded; `selection_rule` per omic with the kept-mean quantiles of both exports when a rule differs from `variance`; the residualization covariates and R² quantiles when present. `render_preprocessing_sensitivity` renders whichever blocks exist. `sensitivity.csv` keeps its columns (the archived `keep_chrx` tables must stay regenerable byte-identically; the task checks this on the synthetic run).

**D7. The adoption check is a constant rule over generated numbers.**
`ADOPTION_CHECK_VARIABLES = ("tmm_rna", "tmm_atac", "fraction_mitochondrial_umis_mean", "doublet_score_mean", "gex_reads_mapped_confidently_to_intronic_regions_mean", "gex_fraction_of_transcriptomic_reads_in_cells_mean", "rin", "pmi")` (the technical-table column names as the diagnostic change writes them; none is a regressor). For `revised_candidate` the runner computes, from the candidate's `latent_scores.csv`, the primary export's manifest (`tmm`) and `covariates.csv` (RIN, PMI) and the technical table, the Spearman correlation of every component with every variable (the diagnostic change's function), flags `|rho| ≥ DIAGNOSTIC_CORRELATION_THRESHOLD`, and writes `adoption_check.csv` with a final `holds` row. `ADOPTION_RULE` is the string the manifest and report quote: "holds when no candidate component reaches the threshold against any check variable, in both cohorts; p-values, rejections and CV AUROC are not read". A variable missing for a cohort donor is an error, not a skipped cell.
- *Alternative: include a stage-signal floor (CV AUROC not below the primary's).* Rejected as a rule: it is a label-based performance number; it is reported beside the check instead.

**D8. CLI and registry order.**
`--variant` gains the five names; `all` runs the registry in order (`keep_chrx` is skipped with a message when its comparison directory exists and `--force` is absent, so `all` is resumable). Output layout, `.tmp` renames and `--force` semantics are the archived change's.

**D9. A separate dated report.**
`reports/motco-astrocyte-feature-selection-<date>.md`, assembled by a scratch script from the ten `sensitivity.md` files, the two `adoption_check.csv`, the comparison manifests and, as descriptive checks, an `axis-diagnostics` run on each variant (`--input-dir <variant>/motco_input --analysis-dir <variant>/motco_analysis`) for the first axis's technical correlations, so the attribution table is generated. Hand-written: provenance, the reading rules the spec binds, the attribution reading, the adoption outcome, limitations. A check script finds every number in the generated files.

## Risks / Trade-offs

- [A variant rejection gets read as the result] → registry constant, primary row first, spec wording; the adoption rule reads no p-value.
- [Residualization removes stage or sex signal carried by the technical set] → stated limitation; the attribution table shows the first axis's stage correlation per variant beside the technical correlations.
- [The composition axis is partly biology and the candidate deletes it] → the technical-axis report's cell-state rule removes such fractions from the set before this change runs; what remains is stated.
- [The trend-ratio rule favours the high-mean tail if the trend is misfit at the ends] → flat extrapolation beyond the end bins; the kept-mean quantiles are reported for both exports.
- [Byte-identity of `export` breaks through the new keywords] → defaults reproduce the primary path exactly; synthetic before/after and the `nixos-home` reproduction; manifest differences limited to the three recorded defaults.
- [Technical table made from a different pseudobulk than the export reads] → D5's hash check.
- [Donors missing a technical value] → error naming the donor; no complete-case subset, since the cohorts are preregistered.
- [Ten CV analyses of two minutes each plus five diagnostics] → one workstation session; `all` is resumable.
