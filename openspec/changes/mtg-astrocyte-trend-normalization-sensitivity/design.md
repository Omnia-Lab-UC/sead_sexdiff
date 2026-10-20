# Design

## Context

See `proposal.md` for the motivation. The state that shapes the approach:

- **The export's per-block pipeline** is `preprocess_block(counts, names, chromosomes, *, k, min_cpm, prior_count, top_n, excluded_chromosomes, excluded_gene_sets, selection_rule, residualize)`: (2) chromosome and gene-set removal, (3) abundance, (4) TMM, (5) log2-CPM, (5b) residualization on a design matrix when given, (6) top `N` by the selection rule. `_export_cohort` resolves the per-omic settings, builds the residualization design once per cohort, calls `preprocess_block` per omic, writes the seven files plus `residualization.csv`, and records `preprocessing_order`, `thresholds`, the settings and `removed_by_gene_set` / `trend_bins` / `residualization` blocks in the manifest. `_write_debug` writes the filtered counts, TMM factors and log2-CPM per omic for `scripts/check_tmm_edger.R`.
- **The sensitivity runner** reads a registry entry's `export` dict, passes it to `_export_cohort`, runs `_analyze_cohort` with the `analyze` parameters, and writes the comparison through `_variant_row`, `_feature_accounting`, `render_preprocessing_sensitivity` and, for an entry with `adoption_check: True`, `_adoption_variables` + `adoption_check` into `adoption_check.csv`. `_needs_technical` decides whether the technical table is read. The CLI's `--variant` choices come from the registry.
- **The check of 2026-10-20** (`scratchpad/ma_trend.py`, run on this laptop against the kept exports, which it reproduced exactly through `preprocess_block`) defines the tilt statistics this change records: per donor, the OLS slope of M on A over the abundance-passing features, and the median M of the lowest and highest of twenty equal-count abundance bins. Both are label-free.
- **limma's `normalizeCyclicLoess(x, span = 0.7, adaptive.span = TRUE, iterations = 3, method = "fast")`** (limma 3.68.4 is in the Nix shell through `rPackages.edgeR`): with `method = "fast"`, each iteration computes the per-feature mean `a` over samples, and for each sample fits `loessFit(x[, i] − a, a, span)` and subtracts the fit; `adaptive.span` sets the span from the number of features by `chooseLowessSpan(n, small.n = 50, min.span = 0.3, power = 1/3)` (≈ 0.40 for 20,000 genes, ≈ 0.35 for 129,000 peaks). `loessFit` without weights calls R's `lowess` (Cleveland's algorithm, 3 robustifying iterations, `delta = 0.01 · range(x)`). statsmodels' `lowess(endog, exog, frac, it, delta)` is the same algorithm; statsmodels 0.14.5 is already in `uv.lock` as a transitive dependency of scanpy.
- **AGENTS.md rules**: the export uses no sex, stage or covariate; variants are registry entries written in a proposal before any run; the primary export's output must stay byte-identical; TMM / log2-CPM changes are cross-checked in R; `n_jobs = 1` for reported results.

## Goals / Non-Goals

**Goals:**
- A `normalization` keyword per omic on `preprocess_block` / `_export_cohort`, default `none`, with the primary's output byte-identical.
- A numpy cyclic loess that matches limma's fast method to a recorded tolerance on the real blocks, with the recorded span rule, iterations and reference.
- Two registry entries, the adoption check on both, the precedence rule as a constant, and the comparison / report additions, all reusing the archived change's machinery.

**Non-Goals:**
- A per-gene-class normalization (the lncRNA-versus-coding offset at matched abundance). Recorded as a finding; a later proposal.
- Changing the primary export, the adoption rule, its variables or the threshold.
- Any change to the analysis settings, the six analysis covariates or MOTCO.
- A per-library (rather than per-donor) normalization: the export works on donor pseudobulks and has no per-library counts.

## Decisions

### D1. Position: after log2-CPM, before residualization and selection, on all abundance-passing features

The tilt is a property of the normalized log-ratios, so it is corrected where TMM's single factor stops: on the step-5 log2-CPM, over the step-3 feature set (the same set TMM is estimated on), not on the selected set. Residualization (5b) and selection (6) then see the corrected values, which is what the candidate needs and keeps the primary's order. *Alternative:* an abundance-dependent offset inside the count model (csaw's `offset` path). It would need a count-level loess and a different log-CPM; the log-scale correction of the same values is what limma's cyclic loess does and is what the R cross-check can compare.

### D2. Method: limma's fast cyclic loess, iterated three times, against the per-feature mean, adaptive span

The reference is the per-feature mean over donors (the "average donor" of the check), updated at each of three iterations; the span follows limma's default adaptive rule from the number of abundance-passing features, so the RNA and ATAC blocks get their own spans (recorded). Choosing limma's defaults means the cross-check compares like with like and no span is tuned on these data. *Alternatives:* the full pairwise cyclic loess (every pair of donors; 3,160 pairs × 129,000 peaks × 3 iterations is minutes-to-hours for no gain over the fast method at n = 80); the 20-bin median trend already in `trend_ratio` (coarsest where the tilt is steepest, and not an established normalization); a fixed span (no reason to depart from the default). The loess is `statsmodels.nonparametric.lowess(M, A, frac = span, it = 3, delta = 0.01 · (A.max() − A.min()), return_sorted = False)`, matching R's `lowess` arguments; identical ties in A are handled by `lowess` itself.

### D3. The `normalization` keyword and constants

`normalization: dict[str, str] | str | None` on `_export_cohort` (per omic, one rule for both, or `None` → `none`), resolved like `selection_rule`; `normalization: str = "none"` on `preprocess_block`. Constants: `NORMALIZATION_RULES = ("none", "cyclic_loess")`, `LOESS_ITERATIONS = 3`, `LOESS_SPAN_RULE = {"small_n": 50, "min_span": 0.3, "power": 1/3}` (limma's `chooseLowessSpan` defaults, applied as `min_span + (1 − min_span) · (small_n / n) ** power`, capped at 1). The manifest records `normalization: {rna: ..., atac: ...}` always (the primary records `none`, as it records `selection_rule: variance`), and a `normalization_fit` block (span and iterations per omic, reference `feature_mean`) only when a rule is not `none`; `preprocessing_order` gains the line between log2-CPM and residualization only then. `PREPROCESSING_PRIMARY_SETTINGS` gains `normalization: {rna: none, atac: none}`.

### D4. Tilt accounting: `normalization_trend.csv` and the comparison

`preprocess_block` returns, when it normalizes, per donor the four tilt statistics before and after (slope of M on A; median M of the lowest and highest of twenty equal-count bins of A; A and M defined against the per-feature mean of the values at that stage) and the largest absolute correction. `_export_cohort` writes them as `normalization_trend.csv` (columns `donor_id`, `omic`, `slope_before`, `slope_after`, `low_bin_before`, `low_bin_after`, `high_bin_before`, `high_bin_after`, `max_abs_correction`) and records the quantiles (5, 25, 50, 75, 95) over donors in the manifest; `_feature_accounting` carries them into the comparison and `render_preprocessing_sensitivity` renders a "Normalization" block, as it renders "Residualization". The statistics are label-free, so the export may compute them; their correlation with `pls_0` is the report's business (`axis-diagnostics` on the variant root, as in the archived change's task 4.4).

### D5. Two entries and the check on both

`trend_normalized = {"export": {"normalization": "cyclic_loess"}, "adoption_check": True}`; `revised_candidate_trend` = the archived candidate's `export` dict plus `"normalization": "cyclic_loess"`, `adoption_check: True`. `_needs_technical` is already true for the candidate (gene sets, residualize, adoption check) and becomes true for `trend_normalized` through `adoption_check` (the check variables come from the technical table). Running the check on the single-setting variant costs nothing and answers the question the precedence rule asks; it is declared here, before any run. *Alternative considered:* a third entry, loess plus the archived candidate's settings without the residualization. Not added: the two entries bracket it (loess alone; loess with everything), and a third candidate would need its own place in the precedence rule, which should stay a two-case rule.

### D6. `ADOPTION_PRECEDENCE`

A constant string beside `ADOPTION_RULE`, applied by the step from the `holds` values of the adoption tables of both cohorts and written into each candidate's comparison manifest as `adoption_precedence: {candidate: <name> | null, reason: ...}` once both cohorts and both candidates have run (otherwise `null` with the reason "incomplete"). `trend_normalized` wins ties because it changes the least (no gene removed, no degrees of freedom spent per feature, no technical covariate read); the rule never reads a p-value, a rank or a CV AUROC.

### D7. R cross-check

`scripts/check_cyclic_loess_limma.R <debug_dir> [cohort] [omic ...]` reads `<cohort>_<omic>_logcpm.csv` and the new `<cohort>_<omic>_normalized.csv` written by `_write_debug`, runs `limma::normalizeCyclicLoess(t(logcpm), span = <recorded span>, adaptive.span = FALSE, iterations = 3, method = "fast")` and reports the largest absolute difference per omic. Tolerance: 1e-6 on log2-CPM; statsmodels and R share the algorithm but not the floating-point order, and R's `lowess` interpolates between fitted points within `delta` while statsmodels does the same, so differences beyond rounding would point to a span, delta or iteration mismatch. If the tolerance is not met, the task records the observed difference and the cause, and the proposal's tolerance is revised by decision, not silently. The script runs inside the Nix shell (limma present); `flake.nix` is unchanged.

### D8. Reuse of the comparison and report machinery

`_variant_settings_record` serializes `normalization` like `selection_rule`; `render_preprocessing_sensitivity` adds the block of D4; the report script of the archived change (paste of `sensitivity.md`, `adoption_check.csv`, the attribution table from `axis-diagnostics` on the variant roots) is reused with the two new variants and the precedence outcome. The report's descriptive section reproduces the per-donor tilt check of the proposal from `normalization_trend.csv` beside `pls_0`, so the proposal's scratch numbers become step outputs.

## Risks / Trade-offs

- [The normalization removes a real stage or sex effect shaped like a tilt] → Stated as the variant's main limitation; the variant is sensitivity with the adoption rule, never the primary; a candidate becomes primary only by a later change, with the current primary kept as the reference.
- [The curve is least determined at the rare end, where the selection draws most features] → limma's default span and robustifying iterations; the tilt statistics before and after are recorded per donor so the report can show what the curve did at the low end; `max_abs_correction` is recorded.
- [statsmodels and limma differ beyond the tolerance] → The cross-check is a task with a recorded outcome; the cause is diagnosed (span, delta, iteration count, tie handling) and either fixed or the tolerance revised by decision.
- [ATAC runtime: 129,000 peaks × 80 donors × 3 iterations] → statsmodels' `lowess` with `delta` is O(n) per fit after sorting; about 240 fits of ~130,000 points per cohort, minutes. Recorded in the manifest's runtime.
- [The lncRNA-versus-coding offset survives and the RNA check fails on it] → Expected and stated in the proposal; the failure is informative (it names the next variant) and does not change this change's scope.
- [Double adjustment in `revised_candidate_trend`: loess, then residualization on library sizes and fractions, then the six analysis covariates] → Inherited from the archived candidate, stated; `trend_normalized` has none of it, which is one reason for the precedence rule.
- [Uncommitted implementation at run time, as in the earlier changes] → Deterministic outputs; a rerun after commit reproduces them (task).

## Open Questions

None that change the specs or the tasks. The tolerance of D7 may be revised after the cross-check runs, by decision recorded in the task.
