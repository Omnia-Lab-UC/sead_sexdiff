# Design

## Context

See proposal.md — Why. The pieces already exist: raw donor × feature counts in `results/pseudobulk/{rna,atac}_counts.h5ad` (`obs`: `sex`, `stage`, `n_rna`, `n_atac`, cohort flags), the export's per-feature status tables (`features_<omic>.csv`: `sex_chromosome` / `abundance` / `variance` / `kept`), the export's `covariates.csv` and `metadata.csv`, the primary analysis's `result.json`, the reference checks `_check_reference` (analysis ↔ export hashes, MOTCO version, seed, permutations) and `check_export_inputs` (export ↔ pseudobulk hashes). pydeseq2 0.5.3 is pinned (`DeseqDataSet(design=<formula>, continuous_factors=...)`, `DeseqStats(contrast=<vector>)`, Wald only). The legacy PFC `differential_expression` uses pydeseq2 on a cell × supertype pseudobulk with a hard-coded contrast vector; it is left alone.

## Goals / Non-Goals

**Goals:**
- One function, `run_stage_de`, that fits the four registered models per cohort × omic on the export-bound feature set, writes one table per tested coefficient, a summary and a manifest, evaluates the concordance rule, and reads the astrocyte directories by default but any cell-type root through its directory arguments.
- Every number in the report traceable to a generated file.

**Non-Goals:**
- A joint test of the two factor-interaction coefficients (no LRT in pydeseq2 0.5.3).
- Effect-size shrinkage, enrichment, peak-to-gene links, plots.
- Using the DE results anywhere upstream.

## Decisions

- **D1 — Feature set from `features_<omic>.csv`, not recomputed.** The step reads the export's status table and keeps `variance` + `kept`, in the export's order. Alternative: redo the chromosome and abundance filter from the annotation. Rejected: it would duplicate the export's code path and need the GTF; the status table is hash-bound to the pseudobulk through the export manifest, so reading it binds the DE to the exact primary filtering.
- **D2 — DESeq2's own normalization on the abundance-passing set.** Size factors (median of ratios) and dispersions are fitted on the ~20k genes / ~125k peaks that passed abundance, not on MOTCO's 5,000 or on all 36,601 / 218,882 features. Alternative: all raw features. Rejected: thousands of all-zero or near-zero peaks and the sex-chromosome genes would enter the dispersion trend. The TMM factors of the export are not reused: DESeq2 is its own method and the manifest says so.
- **D3 — Design formulas through pydeseq2's formula interface.** `design="~ sex + stage + sex:stage + age_at_death + pmi + rin + years_education + log10_n_rna"` with `continuous_factors` listing the five covariates (and `stage_num` in the trend models); `sex` and `stage` as string categoricals with the reference levels set by `ref_level` / ordered categories so the design columns are `sex[T.Male]`, `stage[T.1]`, `stage[T.2]`, `sex[T.Male]:stage[T.1]`, `sex[T.Male]:stage[T.2]`. Each tested coefficient is a unit contrast vector on `dds.obsm["design_matrix"].columns`, located by column name, never by position. If the installed pydeseq2 names the columns differently, the fit helper resolves the names once and records them in the manifest; a name that cannot be resolved is an error. Verified at task 2.2 on a synthetic dataset.
- **D4 — Registry `STAGE_DE_MODELS`** maps model name → `{"formula": ..., "stage_coding": "factor" | "numeric", "tested": [(label, column), ...], "interaction": bool}`. `STAGE_DE_COVARIATES = {"rna": (..., "log10_n_rna"), "atac": (..., "log10_n_atac")}`. The union per factor model is computed in the summary, not as a fifth test.
- **D5 — Two reporting sets.** Every table carries `exported`; the summary counts on `exported` (primary) and `abundance` (all tested). `padj` is DESeq2's over all tested features, not recomputed on the exported subset, so the exported-set count is "features MOTCO measured that are significant in the full test", which is what the concordance rule needs.
- **D6 — Concordance rule** `STAGE_DE_CONCORDANCE_RULE` evaluated from `stage_de_summary.csv` rows (`factor_interaction`, union, `exported`) and the primary `result.json` (`p_values.delta ≤ ANALYSIS_ALPHA` under `ANALYSIS_REJECTION_RULE`). Written as a row per cohort × omic into the summary and into the manifest. The rule names only the factor model, so the trend model's interaction count is reported beside it but does not enter the rule.
- **D7 — Layout and flags.** `results/stage_de/<cohort>/<omic>/<model>__<coefficient>.csv`, `stage_de_summary.csv`, `stage_de.md`, `stage_de_manifest.json`. CLI: `--cohort`, `--omic {rna,atac,both}`, `--pseudobulk-dir`, `--input-dir`, `--analysis-dir`, `--out-dir`, `--n-cpus`, `--force`. `--n-cpus` affects runtime only; task 3.3 checks it does not change `padj`.
- **D8 — Reference checks reused.** `check_export_inputs`-style hash comparison of the pseudobulk files against the export manifest's `input` block and `_check_reference` for the analysis, run before any fit. The step writes nothing on failure.
- **D9 — Rendering.** `render_stage_de` writes `stage_de.md` from the summary and manifest only (counts table per omic, concordance table, covariates, models, findings), so the report can paste it.

## Risks / Trade-offs

- [Two Wald tests per factor model inflate the union count] → the union is labelled as such in every output and the per-coefficient counts are always beside it; the concordance rule is predeclared on the union and the report says what that means.
- [ATAC: ~125k peaks × 4 fits is slow and memory-heavy in pydeseq2] → `n_cpus`; fits run one model at a time and free the `DeseqDataSet`; the step is resumable per cohort × omic × model (an existing table is skipped without `--force`).
- [DESeq2 normalization differs from the export's TMM] → recorded in the manifest and stated in the report; the feature set is the same, the normalization is DESeq2's by design (D2).
- [Interaction tests at n = 80 with cells of 9–28 have little power] → the report reads a zero interaction count as "none detected at this n", never as "no interaction"; the stage main-effect counts are the informative description.
- [pydeseq2 column naming or `ref_level` semantics differ from the assumption in D3] → resolved by name with an error on mismatch; checked on synthetic data before real data.

## Open Questions

None that change the specs or tasks. Whether `independent_filter` removes many low-count peaks is read from the output (`padj` NaN count, recorded in the manifest).
