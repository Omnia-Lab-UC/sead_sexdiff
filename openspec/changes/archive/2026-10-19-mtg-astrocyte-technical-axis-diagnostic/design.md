# Design

## Context

See `proposal.md` for the motivation. The state that shapes the approach:

- **The pseudobulk's per-donor read** (`pseudobulk_donor`, archived change `2026-10-16-mtg-astrocyte-donor-pseudobulks`) opens one donor object with h5py, reads the `obs` columns `Donor ID`, `Sex`, the ADNC column, `Subclass` and `Supertype`, maps the release's taxonomy, excludes `Reference` donors, reads the astrocyte rows of the `UMIs` layer with `_read_rows`, drops nuclei with fewer than `RNA_MIN_GENES` detected genes, sums the rest and caches `n_astro`, `n_nuclei`, the summed `row` and `supertype_counts` in `results/pseudobulk/per_donor/rna/<donor>.npz`. The cache exists for all 84 RNA donors on the cluster and is the ground truth for "the nuclei the pseudobulk summed".
- **The covariate table is frozen by hashes.** `_check_export_reference` compares `input.covariates_manifest_sha256` in the primary export manifest with `results/pseudobulk/covariates_manifest.json`, and `_check_reference` binds the five analysis inputs, `covariates.csv` among them, to the primary analysis. Adding a column to `donor_covariates.csv` would invalidate the reference checks of the rank sweep and of `preprocessing-sensitivity`. Technical variables therefore need their own files.
- **What the donor objects hold.** RNA `obs` has per-nucleus QC (`Fraction mitochondrial UMIs`, `Number of UMIs`, `Genes detected`, `Doublet score`) and library-level 10x metrics (`GEX_*`, `load_name`, `library_prep`, `pcr_cycles`) that take up to four distinct values within one donor object (one donor spans up to four libraries), so library metrics have to be summarized over the astrocyte nuclei, not read as a single value. ATAC `obs` carries the `ATAC_*` library metrics, constant within a donor object (one library), plus multiome GEX columns that this change does not read. The RNA `UMIs` layer has 36,601 genes; the pseudobulk `rna_counts.h5ad` keeps them all, including the 13 `MT-` genes and 103 `RPL`/`RPS`-prefixed names (pseudogenes included, hence the name pattern).
- **What the analyses write.** `latent_scores.csv` (`donor_id`, `pls_k`), `pls_loadings.csv` (`block`, `feature`, `pls_k`), `stage_means.csv` (`group`, `stage`, `pls_k`), `result.json` and `analysis_manifest.json` with the input hashes; the forced-rank analyses under `results/motco_rank_sensitivity/<cohort>/rank_kk/` have the same files. The export manifest holds `tmm` per omic and donor (`lib_size_filter`, `lib_size`, `factor`) and `feature_counts`; `features_<omic>.csv` holds `status` and the variance of every abundance-passing feature; the block CSVs hold the kept features' log2-CPM. Everything the 2026-10-19 scratch checks used is therefore in kept files, except the gene-set fractions.
- **Where the data are.** The cluster holds the raw objects and `results/pseudobulk/`; `nixos-home` holds `results/motco_input/`, `results/motco_analysis/`, `results/motco_rank_sensitivity/` and `results/motco_sensitivity/` and pulls `results/pseudobulk/` with `scripts/sync_pseudobulk.sh`, which excludes `per_donor/`.

## Goals / Non-Goals

**Goals:**
- One cluster command that writes the technical table from the same nuclei as the pseudobulk, resumable, without touching the covariate files.
- One workstation command that writes the axis, geometry and feature-selection diagnostics from kept outputs, deterministic, with input hashes.
- A report whose every number is in a generated file, read under fixed rules, that gives the next change its residualization candidates and baseline values.

**Non-Goals:**
- Reading per-nucleus ATAC counts: no per-nucleus fragments-in-peaks metric exists in the objects and a peak-set fraction would be a new quantity without an established technical meaning.
- The multiome GEX columns of the ATAC objects: the RNA objects are the RNA source.
- Any adoption rule for a preprocessing variant: it is written in the next change's proposal, before that change runs.
- A pytest suite; verification stays synthetic checks plus the real runs.

## Decisions

**D1. Extract the nucleus selection from `pseudobulk_donor` into `select_astrocyte_nuclei`.**
`select_astrocyte_nuclei(f, path, omic, *, min_genes, reference_features)` returns the donor header (`donor_id`, `sex`, `adnc`, `excluded`, release), the selected nuclei's `Supertype` values and the CSR matrix of their raw counts after the nucleus filter (RNA: `UMIs` rows with ≥ `min_genes` genes; ATAC: `X` rows, no filter), or an excluded record. `pseudobulk_donor` becomes that call plus the sum and the supertype counts; `technical_donor` becomes that call plus the per-nucleus fractions and the `obs` QC read for the same rows.
- *Alternative: duplicate the selection in `technical_donor`.* Rejected: "the same nuclei" would be a claim, not a property of the code.
- *Verification:* a synthetic donor object (two releases' taxonomy labels, a `Reference` donor, nuclei below and above `min_genes`) summed before and after the refactor gives identical `row`, `n_astro`, `n_nuclei` and `supertype_counts`; on the cluster, `pseudobulk_donor` on two real donors is compared with their cached `.npz` (`row` equal, counts equal). `pseudobulk` output is not regenerated.

**D2. Gene sets are name patterns, recorded as matched lists.**
`TECHNICAL_GENE_SETS = {"ribosomal": (r"^RP[LS]\d+[A-Z]?$", r"^RPLP[0-2]$", r"^RPSA$"), "mitochondrial": (r"^MT-",), "hemoglobin": (r"^HBB$", r"^HBA1$", r"^HBA2$")}` applied to `var_names`; the manifest lists every matched gene per set. The ribosomal pattern excludes pseudogenes (`RPL13AP5`) and kinases (`RPS6KA1`).
- *Alternative: HGNC gene-family identifiers.* Rejected: the GTF carries no family annotation and a second download would add provenance for no gain at this size.

**D3. Per-nucleus fractions, two per-donor summaries, one biology-versus-preparation statistic.**
Per nucleus: UMIs in the set over total UMIs of the nucleus (all genes, the raw layer). Per donor and set: the nucleus mean and SD, the UMI-weighted pseudobulk fraction (sum over nuclei of set UMIs over sum of total UMIs), and the within-donor supertype share, the between-supertype sum of squares over the total sum of squares of the per-nucleus fraction (η²) over the supertypes present in that donor, undefined (empty) for a donor with one supertype or fewer than two nuclei in every other supertype. QC columns: nucleus mean and SD per donor; library-level columns are summarized the same way (nucleus-weighted), with `n_libraries = ` the number of distinct `load_name` values among the selected nuclei. ATAC columns: nucleus mean, with a per-donor `constant` flag in the manifest.
- *Alternative: one library-level value per donor.* Rejected: a quarter of the donors span several libraries; the astrocyte-weighted mean is the quantity a donor-level model can use.
- *Alternative: a mixed model for the supertype share.* Rejected for this change: η² is descriptive, needs no fitting choices, and the report reads it beside the between-donor share of the same fraction (variance of donor means over total variance of the per-nucleus fractions, pooled over donors).

**D4. The pseudobulk cross-check is the "same nuclei" check.**
For every donor, the UMI-weighted fraction from the per-nucleus rows is compared with the fraction computed from the donor's row of `rna_counts.h5ad` (`np.isclose`, `rtol = COVARIATE_RTOL`), and the selected nucleus count with the cache's `n_nuclei` when the cache is present; either mismatch stops the run naming the donor. The pseudobulk manifest's sha256 is recorded, so the next change can check that the technical table was made from the pseudobulk its variant exports read.

**D5. Separate files, separate cache.**
`results/pseudobulk/donor_technical.csv` (columns: `donor_id`, `n_nuclei_rna`, `n_libraries_rna`, then per set `<set>_mean`, `<set>_sd`, `<set>_pseudobulk`, `<set>_supertype_share`, then the RNA QC means and SDs with snake_case names, then the ATAC means), `nucleus_technical_summary.csv` (one row per donor and set: `donor_id`, `gene_set`, `n_nuclei`, `mean`, `sd`, `pseudobulk`, `supertype_share`, `n_supertypes`), `technical_manifest.json`; per-donor cache `results/pseudobulk/per_donor_technical/<omic>/<donor>.json` holding the summaries and the per-supertype nucleus counts and means (no per-nucleus arrays), skipped on rerun unless `--force`. `donor_covariates.csv` and `covariates_manifest.json` are never opened for writing; the task verifies their hashes before and after.
- *Alternative: add the columns to `donor_covariates.csv`.* Rejected: see Context; it would invalidate the reference checks of two archived steps.

**D6. The diagnostic's statistics are the scratch checks, defined once.**
`DIAGNOSTIC_VARIABLES`: every numeric donor column of `donor_technical.csv`; `tmm_<omic>` and `log10_lib_<omic>` from the export manifest's `tmm`; the six analysis covariates as the analysis derives them (`log10_n_rna`, `log10_n_atac` from `covariates.csv`); `stage` (0, 1, 2) and `male` (1 for Male). Spearman correlation via scipy; the partial correlation is the Spearman correlation of the OLS residuals of both variables on an intercept and the six analysis covariates (the definition is written in the manifest). Block PCs: SVD of the z-scored block (population SD), three components, variance share recorded. Geometry shares from `stage_means.csv`: for a vector `v` in component coordinates the share along component `j` is `v_j² / |v|²`, computed for each sex's net displacement (last stage minus first), for the three per-stage sex-difference vectors (summed squares) and for the stage-centred sex differences. Feature accounting from `features_<omic>.csv` (`status`, `variance`) and the block CSV (means). `DIAGNOSTIC_CORRELATION_THRESHOLD = 0.3` flags table rows; no p-value is attached to any correlation and the renderer says so.
- *Alternative: Pearson correlations.* Rejected: the technical variables are skewed (library sizes, fractions) and the scratch checks that motivate the change were Spearman; the two must agree.

**D7. Forced-rank analyses are read when present, checked the same way.**
`_check_analysis_inputs(export_dir, analysis_dir)` compares every hash in an analysis manifest's `inputs` with the export files, and is applied to the primary analysis and to every `rank_*/` directory found; a mismatch is an error. The nestedness finding compares each forced rank's `pls_0` with the primary's (`np.allclose`).

**D8. Output layout and determinism.**
`results/motco_diagnostics/<cohort>/` written as `<cohort>.tmp/` and renamed; existing output is an error without `--force`. CSV floats use `_csv_float`; rows are sorted by analysis, component and variable name; the manifest records every input path and sha256, the variable list with sources, the threshold, the analyses read, findings, software versions and git revision. `diagnostics.md` lists flagged correlations first, then the full tables.

**D9. CLI and machines.**
`sead-sexdiff technical-covariates [--omic {rna,both}] [--data-root …] [--pseudobulk-dir results/pseudobulk] [--out-dir results/pseudobulk] [--force]` runs on the cluster: `scripts/technical.sbatch` mirrors `pseudobulk.sbatch` (partition `512x1024`, sequential, resumable), and a login-node `nohup` run is acceptable since the per-donor read is the pseudobulk's. `sead-sexdiff axis-diagnostics [--cohort {cohort_30,cohort_50,all}] [--input-dir results/motco_input] [--analysis-dir results/motco_analysis] [--rank-dir results/motco_rank_sensitivity] [--pseudobulk-dir results/pseudobulk] [--out-dir results/motco_diagnostics] [--force]` runs on a workstation after `scripts/sync_pseudobulk.sh`, whose excludes gain `per_donor_technical/`.

**D10. A separate dated report.**
`reports/motco-astrocyte-technical-axis-<date>.md`, assembled by a scratch script from both cohorts' `diagnostics.md`, the donor-level technical summary (distribution of each fraction across donors, between-donor share, median supertype share) and the two manifests, with hand-written provenance and the reading rules the spec binds. A check script finds every number of the tables in the CSVs.

## Risks / Trade-offs

- [Eighty donors, dozens of variables, correlations read as results] → the threshold is a reporting flag; no p-values; the report states that the tables describe the measurement space and support no inferential claim.
- [The diagnostic reads sex and stage and could steer a later choice] → the spec forbids any downstream read of its outputs and documents that variants are predeclared in a proposal; the next change writes its adoption rule before running.
- [The supertype share confounds cell state with supertype-specific nucleus quality] → stated as a limitation; the report reads it beside the between-donor share and the nucleus QC columns, not alone.
- [Donors without an ATAC object, or with missing QC values] → empty cells listed in the manifest; completeness within a cohort is the next change's requirement for the residualized variant.
- [The refactor changes `pseudobulk` output] → D1's synthetic and cached-donor checks; the pseudobulk outputs are not regenerated.
- [Cluster time for 84 UMIs reads] → the per-donor cache; `technical.sbatch`; a stopped run resumes.
- [Cross-machine floating point] → the diagnostics are read from files and written with `_csv_float`; byte-identity is checked between two runs on the same machine only.
