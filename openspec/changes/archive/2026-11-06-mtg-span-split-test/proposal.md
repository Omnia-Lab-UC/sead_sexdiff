# Proposal

## Why

The case study's non-rejections can only see a sex × stage difference inside the fitted stage-PLS span. MOTCO's plasmode on the astrocyte `cohort_30` export (2026-11-02, `attribution_leg_failed`) missed a 20-feature plant at 3 within-stage SDs. The plant lay almost entirely outside that span. The latent-space comparison (2026-11-04, `document_tradeoff`) then found no measurement space that sees such a difference without losing in-span power.

MOTCO 0.18.0 (`span-split-interaction-test`, 2026-11-05, decision **`split_validated`**) asks the out-of-span question with a test of its own. `TrajectoryAnalysisParams.span_split_test` takes the group × stage interaction energy of the analysed features and splits it exactly into three parts: the fitted span (`span`), its orthogonal complement (`complement`) and all of it (`total`). Each part is tested by a trace pseudo-F under Freedman–Lane permutation with RRPP's draws at the analysis's seed, and the record gains a `span_split` block. MOTCO states that every other output is unchanged. Its validation found:
- **Size holds** on three simulated designs (18 control rates 0.032–0.058 against 0.0695). It also holds on real noise at the case study's configuration (six covariates, rank 5, `df`): total 0.040, span 0.050, complement 0.012. That real-noise check ran on the astrocyte `cohort_30` export whose five files hash exactly as ours (`plasmode-fidelity/fidelity.json`).
- **The complement sees a dense out-of-span rotation that `angle` misses:** 1.00 / 0.68 against 0.08 / 0.10, at n = 200 and on the SEA-AD table at a strong stage signal.
- **It has no power for the plasmode's sparse plant** (0.02–0.05), and at SEA-AD's own stage signal no part has power for any plant (0.03–0.08).

MOTCO's roadmap (item 16) assigns this repository the follow-up: apply the test to SEA-AD, pin 0.18.0, and fix a multiplicity policy over the three parts before any real-label analysis. 0.18.0 is tagged `v0.18.0` (commit `bceead1`) and published on PyPI (wheel sha256 `b142a095c9fd6833…`). No real-label split test has been run on these exports. MOTCO's validation used only within-stage-permuted labels, and its fidelity gate reads no statistic.

## What Changes

- **Pin `motco==0.18.0` from PyPI** with the existing rule:
  - the wheel's `motco/` tree must be byte-identical to `src/motco/` at `v0.18.0`;
  - `uv export` may change only the motco line (scikit-learn stays at `==1.8.0`, as in 0.13.0);
  - `ANALYSIS_MOTCO_MIN_VERSION` becomes `0.18.0`.

  Between `v0.13.0` and `v0.18.0` the analysis path gained covariate-adjusted attribution (0.14.0, off here as before), `latent_artifacts` (in memory only) and the opt-in split test. The core statistics modules are unchanged. The pin-move gate checks that nothing else moved.
- **The split test in every analysis.** A new constant `ANALYSIS_SPAN_SPLIT_TEST = True`:
  - **Where it is passed:** explicitly in every `TrajectoryAnalysisParams` the package builds: the primaries, the forced ranks, the preprocessing variants, the `extended_covariates` arm and the exclusion variants.
  - **The one exception:** the `latent_adjusted` arm passes `False`, because MOTCO rejects the test with a latent adjustment (it is defined on the covariate-blind basis only). That arm's manifest records the test as not available, with MOTCO's reason.
  - **Checked after each run, before anything is written:** the block must be present exactly when requested, with the three parts, the analysis's permutations, seed and `df` correction, a span dimension equal to the selected rank, and `draws_match_rrpp` at `n_jobs = 1`. Otherwise nothing is written.
  - **Recorded:** in the manifest's `measurement.span_split`, as `--span-split-test` in `cli_equivalent`, and as a parts table in `summary.md` (pseudo-F p-values only; MOTCO's interaction-SS p-value stays in `result.json` as its diagnostic and is never read).
  - **Checked by readers:** `_check_reference` requires the block on the primary.
  - **In the comparison tables:** the shared row builder gives every table (`rank_sensitivity.csv`, `sensitivity.csv`, `latent_sensitivity.csv`, `exclusion_sensitivity.csv`) the parts' F and p, the complement's status, the span concordance cell and the in-span share. These cells are empty for `latent_adjusted`.
- **Figure 11 `11_span_split`** in `FIGURES`:
  - the three parts' p-values per primary analysis, each marked by its role;
  - the complement p-value of every analysis that carries the test;
  - the in-span share against the fitted rank.

  The consistency pre-pass requires the block in every analysis drawn except `latent_adjusted`. A new reading-rule text restates the policy below.
- **Re-issue under `motco-release-migration`.** Every reader compares MOTCO versions, and the figures need one version throughout, so the pin move re-issues every MOTCO-dependent output of both cell types:
  - the same-machine pin-move gate (0.18.0 without the test against 0.13.0);
  - the contract gate (0.18.0 with and without the test: everything but the block and its parameter equal);
  - the move to `results/superseded/mtg-span-split-test/`;
  - every step rerun in dependency order, followed by the version and measurement consistency check;
  - `figures --report` last.
- **Predeclared reading**, written here before any real-label split test:
  - **Multiplicity policy over the three parts.**
    - The **`complement`** is the one confirmatory part. Its pseudo-F p-value is read at α = 0.05 (p ≤ α), unadjusted. It asks the one question `delta`, `angle` and `shape` cannot.
    - The **`span`** part is a concordance cell. Its status is compared with "any of `delta`, `angle`, `shape` rejected" in the same analysis, and a discordance in either direction is a finding. It is never a result.
    - **`total`**, the in-span share (`span_share`) and every interaction-SS p-value are reported and not read.
    - No adjustment is made across `delta`, `angle`, `shape` and the complement, or across cohorts and cell types, as for the production statistics today. The report states that each primary carries four tests read at α.
  - **Roles.** Each cell type's `cohort_30` primary gives the result, and `cohort_50` is sensitivity. In every other analysis a complement status that differs from the primary's is rank, preprocessing, covariate-set or cohort dependence, not a result of that analysis.
  - **A complement rejection** means that the sexes' stage trajectories differ with interaction energy outside the fitted stage-PLS span. It does not mean a different molecular program. It can also contain:
    - leakage from an in-span difference: a 45° in-span rotation is rejected by the complement in 0.32 of replicates on the SEA-AD table;
    - what the fitted span missed of the true program: at SEA-AD's stage signal the fitted span keeps about 0.10 of the population program.

    It carries no geometry, direction or feature list, and it starts no attribution.
  - **A complement non-rejection is uninformative**, as `angle` and `shape` non-rejections are. No power was shown at SEA-AD's stage signal for any plant, nor for sparse differences at any signal. The dense-rotation power (0.68) was measured at a strong stage signal.
  - **Size** was measured on real noise only for the astrocyte `cohort_30` export at this configuration (complement 0.012, conservative). For `cohort_50`, microglia, the eight-covariate arm, the exclusion variants, the forced ranks and the preprocessing variants it is inferred, not measured. Each is carried as a limitation.
  - **Known when this change was written.** No current analysis rejects any statistic (`delta`, `angle` or `shape`) in any of the 44 analyses, and every stage-DE concordance cell is concordant.

    The test only adds a block, so every other value is expected to reproduce. The primaries are checked by the contract gate and every other analysis by the comparison script. A failed reproduction is a finding and changes no reading.

    Nothing else is re-read: the stage-DE concordance rules (which read `delta`), the positive control (which reads no p-value), the adoption, mediator, flag and axis checks, and the `latent_adjusted` reading.
- **One committed report**, `reports/motco-span-split-<date>.md`, covering both cell types. It holds:
  - provenance (pin, gates, move);
  - MOTCO's validation with its numbers;
  - the split results of every analysis, with the primaries in full;
  - per re-issued step, its outcome table (superseded → new);
  - the reading in the outcome that occurred;
  - findings and limitations;
  - the frozen eleven-figure set.

  Numbers are pasted from generated files and checked by a script.
- **README and AGENTS.md:** the pin, the constant and the policy, figure 11, and the superseded entry.

Out of scope:
- MOTCO's threads 2 (a group-aware space refitted per permutation) and 3 (a sparse per-feature interaction scan);
- attribution;
- MOTCO's orthonormal span coordinates (item 15);
- a stage-DE concordance rule for the complement;
- a null in-span share computed here;
- any change to a registry, expectation, adoption rule, threshold or the production statistics' reading;
- rerunning the MOTCO-independent steps;
- edits to earlier reports.

## Capabilities

### New Capabilities
- `motco-span-split-test`: the split interaction test carried by every analysis that MOTCO allows it in. This covers the registered constant, the record check, the manifest, summary and table columns, the reference check, the predeclared multiplicity policy (complement confirmatory, span concordance cell, total and share descriptive) and the reading and limitations the report carries.

### Modified Capabilities
- `motco-trajectory-analysis`: MOTCO ≥ 0.18.0. The split test joins the fixed measurement settings, and `summary.md` gains the parts table.
- `motco-rank-sensitivity`: the cross-rank table gains the split columns, and the forced runs carry the test.
- `motco-preprocessing-sensitivity`: the comparison table gains the split columns.
- `motco-latent-adjustment-sensitivity`: `extended_covariates` carries the test and `latent_adjusted` cannot (recorded as not available). The comparison table gains the split columns, empty for that arm.
- `donor-exclusion-sensitivity`: the comparison table gains the split columns.
- `motco-release-migration`: the contract gate covers a setting that only adds a record block.
- `case-study-figures`: eleven figures, with `11_span_split` added. The consistency pre-pass checks the split setting, and the binding context includes the split policy.

## Impact

- `pyproject.toml`, `uv.lock`: `motco==0.18.0` from PyPI.
- `src/sead_sexdiff/data.py`:
  - `ANALYSIS_MOTCO_MIN_VERSION` and `ANALYSIS_SPAN_SPLIT_TEST`, plus the policy constants;
  - `span_split_test` in the five `TrajectoryAnalysisParams` constructions;
  - the record check and `measurement.span_split` in `_analyze_cohort`;
  - `_cli_equivalent`, `render_analysis_summary`, `_check_reference`, `_sensitivity_row` and the latent arm registry.

  No CLI flag is added.
- `src/sead_sexdiff/figures.py`: figure 11, its caption, `_measurement` and the consistency check, plus a reading rule. The module's sha256 changes, so every figure is redrawn (it is anyway in the re-issue).
- `README.md`, `AGENTS.md`.
- `results/` (gitignored), on `nixos-home` and the cluster:
  - the superseded outputs are moved, and every MOTCO-dependent output is re-made;
  - the compute is that of the 2026-10-31 re-issue (a few hours), plus 999 O(n²) draws per analysis for the split test.
- `reports/`: one new report with its frozen figure set.
- Scientific choices introduced here (flag and ask before changing any of them):
  - the split test in every analysis that MOTCO allows it in;
  - the complement as the one confirmatory part, at α = 0.05 unadjusted;
  - the span part as a concordance cell against the production statistics;
  - `total`, the in-span share and the interaction-SS p-values reported and not read;
  - no adjustment across the four tests per analysis;
  - the reading and limitations above;
  - re-issuing every MOTCO-dependent output under 0.18.0.
