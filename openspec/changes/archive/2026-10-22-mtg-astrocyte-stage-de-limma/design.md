# Design

## Context

See proposal.md — Why. What exists:

- **`stage-de`** (`run_stage_de`, data.py §"Pseudobulk stage DE") already does everything before the fit, inside `_stage_de_prepare(cohort, ..., log=...)`: input checks naming the producing command, hashing of every file, the pseudobulk-manifest hash against the export, the reproduction of the export from the counts (`_stage_de_reproduce`, which reruns `preprocess_block`), the donors / sex / stage / k / cells cross-checks, the feature statuses, `_check_reference`, and the `delta` status. It returns the tested counts (donors × tested features, export order), the `exported` masks and the per-omic design tables (`stage_de_design`), and it pre-checks each model's design matrix (`stage_de_design_check`, built with `FormulaicContrasts` exactly as pydeseq2 builds it).
- **`stage-de`'s outputs** follow a fixed layout: the per-coefficient tables `<omic>/<model>__<label>.csv`, the fit records (fingerprint plus table sha256), `stage_de_summary.csv` (`STAGE_DE_SUMMARY_COLUMNS`) and `stage_de_manifest.json` (`inputs.files` with every sha256, `tables`, `concordance`, `complete`).
- **The export's TMM** is computed in `preprocess_block` step 4 on the abundance-passing counts, with library sizes equal to those counts' row sums (`lib = c3.sum(axis=1)`). So TMM on `stage-de`'s tested set reproduces it. `scripts/check_tmm_edger.R` showed that edgeR matches it to ≤ 1e-8.
- **R in the Nix shell:** R 4.6.1 with edgeR 4.10.1 and limma 3.68.4. limma is installed only as an edgeR dependency today. In limma 3.68.4 the relevant defaults are:
  - `voomWithQualityWeights(span = 0.5, adaptive.span = TRUE, var.design = NULL, method = "genebygene", maxiter = 50, tol = 1e-5)`;
  - `arrayWeights(prior.n = 10)`;
  - `eBayes(trend = FALSE, robust = FALSE, legacy = NULL)`.
- **rpy2 3.6.4** is declared in `pyproject.toml` and used nowhere. In the Nix shell it imports only in ABI mode: API mode fails with "No module named `_rinterface_cffi_api`" and falls back, printing a message on every import.
- **limma's joint F with weights** (checked in limma 3.68.4's source): `lm.series` stores `cov.coefficients` from the unweighted design's QR, and `classifyTestsF` combines the gene-specific moderated t's through `cov2cor(cov.coefficients)`. With voom's observation weights the coefficients' correlation differs by feature, so limma's F is an approximation.

## Goals / Non-Goals

**Goals:**
- One function, `run_stage_de_limma`, with `stage-de`'s structure (checks for every cohort before any write, resumable per omic × model, summary, Markdown, manifest), so the two steps read alike.
- The same donors, counts, design matrices and feature set as `stage-de`, proven by reusing its code path rather than copying it.
- Every number of the comparison report traceable to a generated file.

**Non-Goals:**
- Any change to `stage-de`'s code path, outputs or fingerprints.
- A plain-voom output, a robust-eBayes output, or any second limma arm.
- Plots (mean–variance trend, weights), enrichment, effect-size shrinkage.

## Decisions

- **D1 — A separate step, not `stage-de --method`.** New subcommand `stage-de-limma` with its own root `results/stage_de_limma/`.
  - Alternative: a `--method` flag on `stage-de` with method subdirectories. Rejected: it would change `stage-de`'s fingerprints, manifest and summary layout. It would also touch a merged capability whose outputs the committed report hashes.
  - A separate step leaves `stage-de` byte-identical (verified on the synthetic data, task 2.5) and makes the method comparison an explicit read of checked outputs.

- **D2 — Reuse, not copy, of the DESeq2 step's preparation and designs.** `run_stage_de_limma` calls `_stage_de_prepare(..., log="stage-de-limma")` unchanged, so every reference check and the reproduction are the same code.
  - The design matrix passed to R is `FormulaicContrasts(design, stage_de_formula(model, omic)).design_matrix`, the matrix `stage_de_design_check` already builds and checks for rank. It goes to R as numbers with its column names (`Intercept`, `sex[T.Male]`, `stage[T.1]`, …), and the tested and joint coefficients are resolved by those names. So the limma design is the DESeq2 design by construction.
  - R does not parse the formula for the fit. As a check, the R script rebuilds each design with `model.matrix` from the design table (factor levels in `STAGE_DE_REFERENCE_LEVELS` order) and the same formula. It requires exact equality under the name map `sex[T.Male]` ↔ `sexMale`, `stage[T.k]` ↔ `stagek`, `:` kept, `Intercept` ↔ `(Intercept)`. This proves that formulaic and R agree on treatment coding; a mismatch stops the model.
  - Private helpers of `stage-de` may be split for reuse (e.g. the table and record writing), but only behind a byte-identity check of `stage-de` on the synthetic data.

- **D3 — Normalization: the export's TMM, recomputed by edgeR and checked.**
  - The R script runs `calcNormFactors(DGEList(counts), method = "TMM", logratioTrim = TMM_LOGRATIO_TRIM, sumTrim = TMM_SUM_TRIM)` on the tested counts, with `lib.size` equal to their column sums and edgeR's other defaults (`doWeighting = TRUE`, `Acutoff = -1e10`, automatic reference column).
  - The factors are compared with the export manifest's `tmm.<omic>.<donor>.factor`. The tolerance is absolute `STAGE_DE_LIMMA_TMM_ATOL = 1e-8`, the tolerance of `check_tmm_edger.R`. Beyond it the step stops before any model of that cohort is fitted.
  - So limma-voom runs on the export's own normalization, the one MOTCO's blocks were built with. DESeq2 used median of ratios.
  - voom's log-CPM is `log2((y + 0.5) / (lib.size × factor + 1) × 1e6)`, not the export's prior-2 log2-CPM. The manifest and report state the difference.
  - Alternative: pass the export's recorded factors to voom. Rejected: recomputing in R and checking binds the normalization twice, the second time across implementations.

- **D4 — `voomWithQualityWeights`, decided yes.** Why, given this data:
  - **The heterogeneity is real.** Library quality varies strongly between donors (proposal Why): sizes span 45× / 23×, TMM factors run 0.73–1.24, and the primary's first latent axis tracks the ribosomal fraction and the TMM factor.
  - **The covariates cover the mean, not the noise.** `rin`, `pmi` and the log10 nuclei count model part of the mean shift. Nothing in either method models a donor whose counts are noisier than others' at matched depth:
    - voom's precision weights follow the fitted count size, so they handle depth and not a donor-level excess variance;
    - DESeq2's Cook's handling was inactive by construction.
  - **What the weights do.** `arrayWeights` gives one weight per donor, estimated from the residuals of the same model and shrunk toward equality (`prior.n = 10`). Noisy donors count less instead of being removed or left at full weight. When donors are of equal quality the weights stay near 1 and cost little (Liu et al. 2015, NAR 43:e97).
  - **Settings:** limma's defaults as listed in Context. One weight per donor (`var.design = NULL`, no `var.group`). `normalize.method = "none"`, because TMM is in the DGEList. Recorded in `STAGE_DE_LIMMA_METHOD` and in every fit record.
  - **Per model:** each model gets its own weights, because they are estimated with its design. Per donor and model, the record holds the weight, plus the minimum, median and maximum and the mean per sex × stage cell. The Markdown shows them.
  - **Alternatives rejected:**
    - plain voom as the method: it ignores the documented heterogeneity. It is run once, in the synthetic calibration only (task 2.3), as the reference for this decision.
    - a second limma arm with plain voom: one more method-dependence axis, twice the outputs, and no rule to choose between the arms.
    - removing low-quality donors: changes the preregistered cohorts.
    - `var.group` weights: no predeclared quality grouping exists.

- **D5 — Fit and moderation with limma's defaults.** `lmFit(v, design)` then `eBayes(fit)`:
  - `trend = FALSE`: voom's weights carry the mean–variance trend.
  - `robust = FALSE`.
  - `legacy = NULL`: resolves to the legacy hyperparameter estimate, because the residual df are equal across features. The record keeps `df.prior`, `s2.prior` and the resolved choice.
  - Alternatives: `robust = TRUE` (Phipson et al. 2016) protects the prior from outlier variances but changes how every other feature is moderated. `treat` tests a different null. Both are rejected so that the method stays limma's canonical voom pipeline, predeclared without tuning.
  - `duplicateCorrelation` / `dream` do not apply: there is one pseudobulk per donor.

- **D6 — Tests, registry and BH family.**
  - **Per tested coefficient** (`STAGE_DE_MODELS[model]["tested"]`, same labels as DESeq2): `topTable(fit, coef = <column>, number = Inf, sort.by = "none", adjust.method = "BH")`. Columns: `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`, and `SE = sqrt(s2.post) × stdev.unscaled`.
  - **Joint tests**, registry `STAGE_DE_LIMMA_JOINT_TESTS = {"factor_additive": ("stage_joint", ("stage_1_vs_0", "stage_2_vs_0")), "factor_interaction": ("interaction_joint", ("male_x_stage_1", "male_x_stage_2"))}`, tied at import to the factor models' tested labels. Run as `topTable(fit, coef = c(a, b), number = Inf, sort.by = "none", adjust.method = "BH")`: a moderated F with df1 = 2, df2 = `df.prior + df.residual`.
  - **Trend models** have one coefficient each, so the t is the test and no joint row exists (F = t²).
  - **Union rows** (either moderated t at adj.P ≤ α) are computed in the summary, as `stage-de` does, for the like-for-like comparison only.
  - **BH family:** all tested features. limma has no independent filter, and the export's abundance filter is the only filter. pydeseq2's per-coefficient baseMean cutoff makes DESeq2's BH families smaller in ATAC. The report states the difference.
  - Labels are suffixed `__<label>.csv` as in `stage-de`, so a joint table is `factor_interaction__interaction_joint.csv`.

- **D7 — Exact-F diagnostic for the joint tests.** For each joint test the R script also computes, per feature g:
  - `F_exact = b' V⁻¹ b / (2 · s2.post_g)`, where `b` is the two coefficients and `V` the matching 2 × 2 block of `(XᵀW_gX)⁻¹`, from the feature's own weights `W_g` (QR of `√w ⊙ X` per feature).
  - Its p-value from F(2, df.total), then BH-adjusted over all features.

  The fit record keeps the largest relative difference from limma's F and the number and names of features whose call at α changes. Each changed call is a finding.
  - limma's F stays the predeclared test. The diagnostic measures the approximation that Context names; it never replaces the F.
  - Cost: one small QR per feature, seconds for RNA and about a minute expected for ~129k peaks.

- **D8 — Limma concordance rule, and how it sits beside the DESeq2 cell.** `STAGE_DE_LIMMA_CONCORDANCE_RULE = {"model": "factor_interaction", "test": "interaction_joint", "feature_set": "exported", ...}`, with texts parallel to `STAGE_DE_CONCORDANCE_RULE`, tied at import to the joint-test registry. `delta_rejected` comes from `_stage_de_prepare` (same derivation).
  - **Per cohort × omic**, the summary row for `factor_interaction` / `interaction_joint` / `exported` carries `interaction_detected`, `delta_rejected` and `concordance`.
  - **Cohort level:** `method_comparison.csv` and the Markdown carry a concordance table: omic | DESeq2 cell (union of two Wald tests, recomputed from the checked DESeq2 tables and required to equal `stage_de_summary.csv`) | limma cell (joint F) | `methods_agree`.
  - **Why the joint F and not limma's union:** the F is the joint test `stage-de` could not run, and it is the only one of the two counts with a per-feature FDR claim on the joint null. Limma's union is reported beside it as a row, not as a cell.
  - **Neither cell replaces the other.** The DESeq2 cell stays the reported cell of `mtg-astrocyte-stage-de`. A disagreement is a method-dependence finding, and a discordance with `delta` in either method is a finding as in `stage-de`.

- **D9 — Method comparison.** `method_comparison.csv`, one row per cohort × omic × model × test × feature set:
  - tests: every tested coefficient (Wald ↔ moderated t), every union (union ↔ union), and every joint F against DESeq2's union of the same model (`kind = joint_f_vs_union`, labelled as different tests);
  - columns: `n_tested`, `n_significant_deseq2`, `n_significant_limma`, `n_both`, `n_only_deseq2`, `n_only_limma`, `detected_deseq2`, `detected_limma`, `detection_agrees`;
  - per-coefficient rows also carry `spearman_log2fc` (DESeq2 `log2FoldChange` ↔ limma `logFC`), `spearman_stat` (Wald `stat` ↔ moderated `t`) and `n_compared` (features with both values finite).

  DESeq2 padj NaN (its independent filter) counts as not significant, as in `stage-de`. The table is descriptive: it carries no p-value of its own and gates nothing. Each `detection_agrees = false` row is a finding.

- **D10 — DESeq2 outputs as a checked reference.** `--deseq2-dir` (default `results/stage_de`); the step reads `<deseq2-dir>/<cohort>/stage_de_manifest.json`. It requires:
  - `complete = true`;
  - `inputs.files` sha256 equal, file for file, to the hashes `_stage_de_prepare` just computed;
  - each table under `tables` hashing as recorded;
  - the summary's `n_significant`, `n_tested` and concordance cells equal to a recomputation from those tables with `_stage_de_summary_rows`.

  Any failure stops the step before any fit, naming the file and `stage-de` (`--force` for a stale run). The DESeq2 directory is hashed before and after, never written, and recorded under `inputs.deseq2` in the manifest.

- **D11 — R bridge: `Rscript` subprocess, not rpy2.** The script is `src/sead_sexdiff/stage_de_limma.R`, shipped in the wheel (hatch includes package files, as it ships `DEG.R`) and found via `Path(__file__).with_name(...)`.
  - **Per fit,** Python writes the exchange files into a temporary directory (or `debug_dir`, Python API only, kept so the script can be rerun by hand):
    - `counts.csv`: integers, features × donors in export and donor order;
    - `design_matrix.csv` and `design_table.csv`: `%.17g`;
    - `settings.json`: formula, tested and joint columns, α, TMM trims, the export's TMM factors, the voom / eBayes settings.
  - It runs `Rscript --vanilla stage_de_limma.R <dir>` with `OMP_NUM_THREADS = OPENBLAS_NUM_THREADS = 1`, so reruns are byte-identical.
  - R writes back `tmm.csv`, one CSV per test (`formatC(x, digits = 17, format = "g")`), `sample_weights.csv` and `fit.json` (versions, df.prior, s2.prior, the adaptive span voom used, the design check, the TMM check, the F diagnostic, runtime).
  - Python checks the exit status. A non-zero status raises with the last lines of stderr, and nothing of that model is written. It parses with round-trip precision, writes the final tables with `_csv_float` (`.tmp` then rename, as `stage-de`) and removes the temporary directory.
  - **Before any write:** `Rscript` must be on `PATH`, otherwise the error names the Nix shell. `Rscript stage_de_limma.R --versions` reports R, limma and edgeR, which enter the fingerprints.
  - **Why not rpy2:**
    - it imports only in ABI mode here;
    - it embeds R's allocator and BLAS in the process that holds the 129k-peak matrices;
    - R errors would surface as Python exceptions from inside R;
    - the repository's three R checks already use `Rscript` on files, so one script with a recorded sha256, runnable by hand, keeps one convention.

    The cost is file exchange, seconds per ATAC fit. rpy2 stays declared and unused (removing it is out of scope).

- **D12 — Layout, resume, CLI.** `results/stage_de_limma/<cohort>/` holds:
  - `<omic>/<model>__<label>.csv` per tested coefficient and joint test;
  - `<omic>/<model>__fit.json`;
  - `stage_de_limma_summary.csv` (`stage-de`'s summary columns; `kind` ∈ {`moderated_t`, `union`, `joint_f`});
  - `method_comparison.csv`, `stage_de_limma.md`, `stage_de_limma_manifest.json`.

  The fingerprint contains everything `stage-de`'s does except pydeseq2, plus the joint tests, `STAGE_DE_LIMMA_METHOD`, the TMM trims, the R script's sha256 and the R / limma / edgeR versions. Resume semantics are `stage-de`'s D7: `complete` / `missing` / `stale`; a stale record is an error without `--force`; `--force` deletes only after every check.
  - The comparison and the summary are rebuilt each run from the records and tables, cheaply. They are not fingerprinted, but the manifest records the DESeq2 hashes they read.
  - CLI: `--cohort`, `--omic {rna,atac,both}`, `--pseudobulk-dir`, `--input-dir`, `--analysis-dir`, `--deseq2-dir`, `--out-dir`, `--force`. There is no `--n-cpus`, because limma runs single-threaded here.

- **D13 — Rendering.** `render_stage_de_limma` writes `stage_de_limma.md` from the summary, the comparison and the manifest only: the reading rules, the side-by-side concordance table, the count tables per omic with the DESeq2 counts beside them, the comparison rows, the sample-weight summaries, the TMM check, the F diagnostic, the inputs and the findings. The report can paste it.

## Risks / Trade-offs

- [Written after the DESeq2 null: a second method could become method shopping] → Every limma setting is fixed here before any limma run. The DESeq2 results stay the reported ones. A call present in one method only is method dependence by rule, and no later change may promote whichever method finds more.
- [Quality weights read the design and may down-weight a sex × stage cell, e.g. if low quality concentrates at stage 2 (stage ρ 0.49 with the primary's first axis)] → Per-donor weights and their per-cell means are recorded and shown, without a gate. The report states each cell's effective down-weighting.
- [limma's joint F is approximate with observation weights] → the exact-F diagnostic (D7) measures it, and changed calls are findings.
- [The calibration simulations follow DESeq2's own NB model, which favours DESeq2] → stated beside the numbers. The heterogeneous-quality condition (task 2.3) is the one designed for the quality-weights rationale. Neither is a gate on real data. The one stop point is limma's realized FDR exceeding DESeq2's on the homogeneous conditions, which is decided before real data is read.
- [voom on low-count ATAC peaks] → the export's abundance filter (CPM ≥ 1 in ≥ k donors) is the filterByExpr-type filter voom expects. The trend span used is recorded.
- [R or limma versions differ between machines] → the versions are in every fingerprint, so a record from another version is stale.
- [Same donors and counts in both methods] → agreement between methods is not replication. The report says so.
- [Power is unchanged: n = 80 / 74, cells 9–28] → a zero count means "none detected at this n"; the joint F does not add donors.

## Migration Plan

None needed. The step is additive, and no existing output, constant or fingerprint changes. Rollback is deleting `results/stage_de_limma/`.

## Open Questions

None that change the specs or tasks. Whether the ATAC exact-F diagnostic fits in about a minute is read from the first synthetic and real runs. If it does not, it moves into vectorized R, with the same definition.
