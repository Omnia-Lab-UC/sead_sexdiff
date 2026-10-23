# Design

## Context

See `proposal.md` — Why. The state that shapes the approach:

- **The flag is in `obs`.** Read on the cluster on 2026-10-22 (`h5py`, obs keys only): the RNA 2026-06-22 objects hold `Severely Affected Donor` (`Y` / `N`), `Cognitive Status`, `Braak`, `CERAD score`, `Thal`, `LATE`, `APOE Genotype`, `Highest Lewy Body Disease`, `Overall CAA Score`; the ATAC 2024-12-06 objects hold the same plus `Continuous Pseudo-progression Score`. Reference donors carry `Reference` in the staging columns and an empty flag. `read_donor_covariates` / `_read_obs_columns` / `_single_value` already read one-value-per-object columns without loading matrices; `run_covariates` is the pattern (per omic, per object, cross-check, manifest).
- **Cohorts and `k`.** `donor_table.csv` holds `donor_id`, `sex`, `stage`, `n_rna`, `n_atac`, `cohort_30`, `cohort_50`; `_cohort_rows(table, covariates, cohort, pb)` gives the cohort's rows; `_export_cohort(cohort, sub, blocks, covariates, annotation, pb_info, ...)` computes `k, k_cell, cells = smallest_cell(sub)` from the rows it is given and records them. So an exclusion export is `_export_cohort` on `sub` minus the flagged donors, and `k` follows by the existing rule. `preprocessing-sensitivity` requires the variant's `metadata.csv` and `covariates.csv` byte-identical to the primary's; that check is specific to export-setting variants and is not reused here.
- **Stage DE binds donors to the donor table.** `_stage_de_prepare` reads `donor_table.csv`, takes the cohort's donors and requires them to equal the export's `metadata.csv` donors (error naming the differing donors). The limma step reuses `_stage_de_prepare` unchanged. The fit-record fingerprints carry the input sha256s, not the donor list.
- **Readers of a variant.** `axis-diagnostics`, `stage-de` and `stage-de-limma` take `--input-dir` / `--analysis-dir` / `--out-dir` and check the analysis manifest's `inputs.files` hashes against `<input-dir>/<cohort>/`, so a variant laid out as `<variant>/motco_input/<cohort>/` + `<variant>/motco_analysis/<cohort>/` is readable, as for the preprocessing variants.
- **E3's technical variables** are `POSITIVE_CONTROL_TECHNICAL_SOURCES`: the technical table's numeric columns and the export's TMM factors and log10 library sizes (`DIAGNOSTIC_VARIABLES["technical"]`, `["export"]`).

## Goals / Non-Goals

**Goals:**
- One obs-only step for the flags, one step for the variant, both per cell type, with the provenance discipline of the other sensitivity steps.
- The variant export is the primary export with the donor rows removed and nothing else; the variant analysis is the primary analysis on it.
- Every primary output and every existing stage-DE output byte-identical.

**Non-Goals:**
- Any adoption rule or revised-primary path.
- Reading the flag anywhere else (covariate, diagnostic variable, filter).
- ATAC stage DE on the variant; sensitivity steps on the variant.

## Decisions

**D1. `donor-flags` modelled on `covariates`.**
```python
DONOR_FLAG_COLUMNS: dict[str, str] = {
    "Severely Affected Donor": "severely_affected", "Cognitive Status": "cognitive_status", "Braak": "braak",
    "CERAD score": "cerad_score", "Thal": "thal", "LATE": "late",
}
DONOR_FLAG_ATAC_ONLY_COLUMNS: dict[str, str] = {"Continuous Pseudo-progression Score": "cps"}
DONOR_FLAG_VALUES: dict[str, bool] = {"Y": True, "N": False}
DONOR_FLAG_OUTPUTS = {"table": "donor_flags.csv", "manifest": "donor_flags_manifest.json"}
```
`read_donor_flags(path, omic)` reads `Donor ID` + the columns (the ATAC-only ones when present, recorded otherwise), one value per object through `_single_value`. `run_donor_flags(*, omic="both", cell_type, data_root, out_dir)` reads every object of each omic, drops donors whose staging columns read `Reference` (listed in the manifest), codes the flag (`DONOR_FLAG_VALUES`, anything else an error naming the donor), cross-checks RNA against ATAC (flag disagreement an error; other columns findings, RNA value used), takes `cps` from ATAC, and writes the two files into `_resolve_dir(out_dir, cell_type, "pseudobulk")`. It writes nothing else. The donor-level table does not depend on the cell type; the directory does, as for `covariates`.
- *Alternative: add the columns to `donor_covariates.csv`.* Rejected: that file's hash is recorded by every export and analysis manifest; rewriting it would make every reference check fail.
- *Alternative: the flag only.* Rejected in favour of recording the staging labels and CPS beside it: the same obs read, no cost, and a later staging proposal needs them; recorded-only, read by nothing.

**D2. Exclusion registry and rule.**
```python
DONOR_EXCLUSION_VARIANTS: dict[str, dict] = {
    "excluding_severely_affected": {
        "description": "the cohort donors minus those SEA-AD flagged as severely affected (Gabitto et al. 2024: QC-outlying high-ADNC donors excluded from the gene-expression tests)",
        "flag": "severely_affected", "exclude_when": True,
    },
}
DONOR_EXCLUSION_CHECK_VARIABLES = ("technical table numeric columns", "primary export tmm_<omic>, log10_lib_<omic>", "primary pls_0")
DONOR_EXCLUSION_AXIS_RULE = (
    "along the variant analysis's first latent component, |Spearman| with stage exceeds the largest |Spearman| with a "
    "technical variable of POSITIVE_CONTROL_TECHNICAL_SOURCES (technical table columns for the variant donors; the variant "
    "export's TMM factors and log10 library sizes); the control cell type's value is the expectation, the reference's descriptive"
)
```
Output names `DONOR_EXCLUSION_*` (`flag_check.csv`, `axis_check.csv`, `exclusion_sensitivity.csv`, `.md`, `_manifest.json`), subdirectories `flag_check/`, `motco_input/`, `motco_analysis/`, `comparison/`; `CELL_TYPE_SUBDIRS["motco_exclusion_sensitivity"]`.

**D3. `run_exclusion_sensitivity(variant="all", cohort="all", *, cell_type, pseudobulk_dir=None, input_dir=None, analysis_dir=None, out_dir=None, data_root=None, n_jobs=1, force=False)`.** Per cohort, before any write: `_load_export_inputs`, `_cohort_rows`, `smallest_cell`, `_check_export_reference` (pseudobulk, covariates, annotation hashes, thresholds, `k`, cells), `load_analysis_inputs` + `_check_reference` (export hashes, MOTCO version, seed, permutations, CV selection), `_check_technical_reference(pb, donors=cohort donors)`, the flag table (`check_manifest_cell_type` on its manifest, a flag value for every cohort donor, sha256 recorded). Then per cohort: `flag_check` → `flag_check/<cohort>.tmp/flag_check.csv`; `sub_v = sub[~flagged]`; `_export_cohort(c, sub_v, ..., out_dir=<variant>/motco_input, final_name=c, donor_exclusion={...}, excluded_chromosomes=EXPORT_EXCLUDED_CHROMOSOMES)` (primary settings; the new keyword is recorded as the manifest's `donor_exclusion` block and is absent by default, so the primary export's manifest is byte-identical); `load_analysis_inputs(exp["dir"])` → `_analyze_cohort(..., params=TrajectoryAnalysisParams(permutations=ANALYSIS_PERMUTATIONS, seed=ANALYSIS_SEED, n_jobs=n_jobs), cell_type=cell_type, findings=[...primary stays primary...])`; `axis_check` from the variant `latent_scores.csv`, the technical table rows of the variant donors and the variant export manifest's `tmm` block; the comparison rows via `_sensitivity_row` plus `n_donors`, `cells`, `k`, `excluded_donors`; rename the four temporaries. Existing output without `--force` is an error; `--variant all` skips complete variant × cohort pairs (resumable, as `preprocessing-sensitivity`).
- *Alternative: a cohort variant inside `PREPROCESSING_VARIANTS`.* Rejected: that registry changes export settings on the primary donors and its spec requires the donors unchanged; a donor change is a different kind of variant and gets its own registry and layout.

**D4. Flag check statistics.** For a variable `x` and the binary flag over the cohort donors (and over the stage-2 donors), the rank-biserial AUC is `P(x_flagged > x_unflagged) + ½ P(=)` (Mann–Whitney U / (n₁ n₀)); the Spearman correlations use `scipy.stats.spearmanr`. Variables: every numeric column of the technical table, `tmm_<omic>` and `log10_lib_<omic>` from the primary export manifest, and the primary `pls_0` from `latent_scores.csv`. Rows: `cohort, variable, source, n_flagged, n_unflagged, auc_all, auc_stage2, spearman_unflagged_stage, spearman_unflagged_pls0` (the last two columns only for the technical variables and stage). Plus the counts of flagged donors per sex × stage cell in the manifest and the Markdown. No p-value anywhere.

**D5. Stage DE on the variant.** `_stage_de_prepare` reads `export_manifest.get("donor_exclusion")`; when present it removes the recorded excluded donors from the donor table's cohort before the agreement check, records the block in the fit manifests (not in the fit-record fingerprint, which already carries the export hashes), and prints it. Without the block the behaviour is unchanged and every existing fit record stays complete. The two steps are then run by hand per cell type and cohort with `--omic rna --input-dir <variant>/motco_input --analysis-dir <variant>/motco_analysis --out-dir <variant>/stage_de` (and `--deseq2-dir <variant>/stage_de --out-dir <variant>/stage_de_limma` for limma). `stage_de_limma.R` is not edited.
- *Why RNA only:* E1 reads RNA; ATAC DESeq2 on 130,000 peaks is the costly part of the DE steps and the question here is the RNA stage signal E1 counted.

**D6. Cross-cell-type comparisons in the report script.** The step is per cell type; the comparisons that need both roots (microglia vs astrocyte variant CV AUROC; each cell type's variant RNA stage counts against its primary counts) are computed by the report's check script from the generated `exclusion_sensitivity.csv`, `stage_de_summary.csv` and `stage_de_limma_summary.csv` files, with the rule written in the proposal. They are compared, not tested.
- *Alternative: a positive-control basis for the variant.* Deferred: `mtg-microglia-latent-sensitivity` introduces a basis registry on `positive-control`; adding an exclusion-cohort basis would make this change depend on that one. If both are merged, a later change may register the basis.

**D7. Byte-identity.** The synthetic dataset builder (scratchpad) runs `export`, `analyze`, `stage-de` and `stage-de-limma` before and after the refactor (CSVs and Markdown byte-identical, JSON equal apart from runtime, timestamps and git; the `donor_exclusion` key absent). A synthetic flag table (two flagged donors) exercises `exclusion-sensitivity` and the two DE steps on its variant.

**D8. Report.** `reports/motco-donor-exclusion-sensitivity-<date>.md`, assembled by a scratch script: the flag table summary (flagged donors per cell type, cohort, sex and stage, from `donor_flags.csv` and the two `donor_table.csv`), both cell types' flag checks, comparison tables, axis checks, the stage-DE counts and concordance cells, the `pls_0` rows of `axis-diagnostics` on the variant (`<variant>/diagnostics/`), the cross-cell-type comparisons (D6); provenance and reading hand-written. A check script finds every statistic, p-value, rank, AUROC, flag-check value, axis-check value and DE count in the text.

## Risks / Trade-offs

- [The flagged donors are mostly female (29 F / 13 M high-ADNC donors), so the exclusion changes the female stage-2 cell more than the male one] → the sex split of the flag is in the flag check and the report; the variant is read as cohort dependence, and the imbalance is a stated limitation.
- [Excluding the flagged donors removes the late epoch: a stage effect concentrated in them disappears with them] → predeclared limitation; the axis check and the DE counts say what remains, and nothing is adopted.
- [The flag's sex × stage cells shrink `k` and the tested feature set changes] → recorded (cells and `k` before and after; the feature accounting of the variant export); the variant export is a different cohort by design, which the reading rules state.
- [The donor-table check in `_stage_de_prepare` is relaxed] → only when the export manifest records an exclusion, which only `exclusion-sensitivity` writes; the recorded donors must still match exactly; D7's byte-identity guards the default path.
- [A flag value other than `Y` / `N` appears in a release] → an error naming the donor; the proposal is amended before any run, never worked around.
- [The cluster is unreachable for `donor-flags`] → the step reads `obs` only; as for the microglia spike, the public bucket over HTTPS range requests with size and ETag checks against the pseudobulk manifest is the recorded fallback.

## Migration Plan

1. Code (D1–D5), synthetic byte-identity and variant checks (D7), on any machine.
2. On the cluster login node: `donor-flags --omic both` and `donor-flags --cell-type microglia --omic both`; sync both pseudobulk directories.
3. On the workstation: `exclusion-sensitivity --cohort all` for both cell types, `axis-diagnostics` on each variant, `stage-de --omic rna` and `stage-de-limma --omic rna` on each variant (both cohorts); copy the variant roots to the cluster.
4. Report (D8), README and AGENTS.md, check, commit.

Rollback: revert the commit; no primary output is touched, and the new directories and the two flag files are removed.
