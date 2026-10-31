# Proposal

## Why

The current primaries are the floored analyses of 2026-10-26 (MOTCO 0.10.0, `reports/motco-rank-floor-primary-2026-10-26.md`), and their RRPP null is too narrow at this configuration.

MOTCO ran three checks after that report:
- **Size check (2026-10-29, `docs/reports/case-study-type-i-2026-10-29.md`).** It measured `motco analyze` where the case study runs: rank floor 2, six covariate columns, a weak stage signal. It returned `diagnose_configuration_excess`. The covariate-adjusted `delta` read 0.060–0.072 in every reading (0.066 pooled) against the 0.0695 bound. The `stage` latent adjustment exceeded the bound at the anchor on `angle` (0.084–0.090) and `shape` (0.078–0.092).
- **Diagnosis (2026-10-30, `case-study-diagnosis-2026-10-30.md`).** It traced the `delta` excess to `residual_df_shrinkage`: RRPP permutes the reduced model's residuals, whose variance is RSS/n rather than RSS/(n − p_red), with p_red = 10 at n = 80. It traced the `stage` adjustment's excess to `residualized_basis`. MOTCO's roadmap then **withheld** this repository's primary until a corrected permutation was released and validated.
- **Corrected-permutation check (MOTCO 0.13.0, change `residual-df-corrected-permutation`, `corrected-permutation-2026-10-30.md`).** 0.13.0 adds the opt-in `residual_correction`. Under `df` every residual is multiplied by κ = √(n / (n − r)). The check ran at this configuration on both case-study `cohort_30` tables, at 500 replicates × 999 permutations, and returned **`adopt_df`**:
  - the corrected null reads 0.029–0.058 pooled against 0.0597, at every reading and every pooled rate;
  - the uncorrected null was liberal on pooled `shape` (0.0685 / 0.0675) and on the strong-signal anchor's `angle` (0.086);
  - the correction costs 0.014–0.054 of power on the target statistics.

  The decision file states the consequence: this repository re-issues both cell types' primaries with `--residual-correction df`, pinned to 0.13.0, in its own change.

MOTCO 0.13.0 is tagged `v0.13.0` and published on PyPI (2026-10-31). This change makes the corrected analysis the primary of both cell types and re-issues every output that depends on a MOTCO analysis, so no output mixes the two nulls. The correction was chosen by MOTCO's predeclared check before any corrected run on these data, not from these results.

## What Changes

- **Pin `motco==0.13.0` from PyPI** (the user's decision), replacing the git dependency on tag `v0.10.0`.
  - The `[tool.uv.sources]` entry is removed, and `uv.lock` records the wheel's sha256.
  - The wheel's `motco/` package is checked byte-identical to `src/motco/` at tag `v0.13.0`, and the result is recorded.
  - No other locked package may move. scikit-learn stays at MOTCO's exact pin.
  - `ANALYSIS_MOTCO_MIN_VERSION` becomes `0.13.0`.
  - AGENTS.md's pin rule becomes: an exact PyPI version, moved with `uv add "motco==X.Y.Z"`, checked against the release tag of `Omnia-Lab-UC/MOTCO`; never a git or path dependency.
- **A recorded residual correction.** A new constant `ANALYSIS_RESIDUAL_CORRECTION = "df"`:
  - **Where it is passed:** explicitly in every `TrajectoryAnalysisParams` the package builds (`analyze`, `rank-sensitivity`, `preprocessing-sensitivity`, `latent-sensitivity`, `exclusion-sensitivity`).
  - **A MOTCO default it deliberately overrides:** MOTCO keeps the correction opt-in (default off). This is the one MOTCO default the package overrides, and it does so by MOTCO's own verdict. It is never a flag.
  - **Checked after each run:** the record must state the correction it applied, with κ. Otherwise nothing is written.
  - **Recorded:** in the analysis manifest's `measurement` (kind, κ, n, reduced rank), in `cli_equivalent` (`--residual-correction df`) and in `summary.md`.
  - **Checked by readers:** `_check_reference` requires it, so no sensitivity, diagnostic or stage-DE step can compare against an uncorrected reference.
- **Same-machine pin-move gate before any corrected run** (spec `motco-release-migration`). The four primaries run under 0.13.0 without the correction and under 0.10.0 from the pre-change commit. Their CSVs must be byte-identical and their records equal apart from versions, times and revisions. Otherwise the change stops.
- **Contract gate before the move.** MOTCO states that the correction moves only the null. On the same machine, the four primaries run under 0.13.0 with and without it:
  - the four CSVs must be byte-identical;
  - `result.json` must be equal apart from the p-values, null summaries and the correction's own record;
  - each recorded κ must equal √(n / (n − 10)).

  Otherwise the change stops.
- **Superseded outputs are moved, not deleted.** Every MOTCO-dependent directory of both cell types goes to `results/superseded/mtg-residual-correction-primary/` with `scripts/move_superseded.py`, on the workstation and on the cluster. The MOTCO-independent outputs (pseudobulks, the covariate, technical and flag tables, the exports, the annotation) are untouched, and their hashes are checked unchanged at the end.
- **Every MOTCO-dependent step is rerun** under 0.13.0 with the correction, from scratch, in dependency order. The steps and commands are those of the 2026-10-26 re-issue (its section 1.5), for both cell types, with `n_jobs = 1`.
- **A version-consistency check at the end.** Every manifest outside `results/superseded/` must record MOTCO 0.13.0. Every analysis record must carry `residual_correction = df`, and either the floor 2 or a forced rank.
- **Predeclared reading**, written here before any corrected run:
  - **Primary.** The corrected analysis is the primary of each cell type. The uncorrected floored analyses are superseded records, not a sensitivity arm, and are not cited as results. `cohort_30` stays primary and `cohort_50` sensitivity.
  - **Binding rules.** They hold unchanged: the CV-selected rank, eigengap and null summary are given beside `angle`; `angle` and `shape` non-rejections are uninformative; a `delta` rejection is a difference, not specifically of size. MOTCO's statement is added: under `df` the p-values are size-controlled at this configuration. `delta` is conservative there (0.029 / 0.031 pooled), so a `delta` non-rejection stays weak evidence. The `angle` null is too wide to inform a non-rejection.
  - **What the correction cannot change.** It moves only the null. Every analysis keeps its latent space, rank, CV AUROC, observed statistics, stage means and eigengap. Every reading that uses no p-value reproduces: `pls_0` and the axis diagnostics, the adoption check and precedence, the mediator, flag and axis checks, E1–E3 and both bases' control outcomes, and the stage-DE tables, counts and calls. The primaries are checked by the contract gate; every other analysis is checked by the comparison script. A failed reproduction is a finding with both values and changes no reading.
  - **Known when this change was written.** No superseded floored analysis rejects any statistic at α = 0.05: no primary, forced rank, preprocessing variant, latent arm or exclusion variant, in either cell type or cohort. Every stage-DE concordance cell is concordant. The correction widens the residual variance of the null (κ about 1.07–1.12), so a new rejection is not expected. A rejection status that differs in either direction is reported as dependence on the null correction. The corrected status is then the result, read under the binding rules, and nothing is re-tuned. The concordance cells re-read the corrected `delta` status under their unchanged rules.
  - **The `latent_adjusted` arm.** MOTCO's diagnosis classes the `stage` latent adjustment at this configuration as `residualized_basis`. That excess is a property of the basis, which the correction was not designed or validated for. From this change on, the arm's `angle` and `shape` p-values are reported but not read: neither a rejection nor a non-rejection is interpretable. The arm stays in `LATENT_SENSITIVITY_ARMS` as a sensitivity analysis. The restated positive control on that basis reads only the CV AUROC and first-component correlations, which no null affects, so its predeclared reading stands.
  - **Size inferred, not measured.** MOTCO measured the corrected null only at six covariates on the two `cohort_30` tables. Its size is inferred there from the mechanism (it depends only on n and the reduced rank) for:
    - the `cohort_50` tables;
    - the eight-covariate `extended_covariates` and `latent_adjusted` arms (r = 12);
    - the smaller-n exclusion variants;
    - the forced ranks;
    - the preprocessing variants.

    Each is carried as a limitation.
  - **Other readings.** Every sensitivity reading is re-read against the corrected primary under its own unchanged rules. The positive control's E1–E3 and all-three rule are unchanged. Its outcome on the corrected primary is the primary-basis outcome, with the superseded outcome beside it.
- **One committed report**, `reports/motco-residual-correction-primary-<date>.md`, covering both cell types. It contains:
  - provenance: the PyPI pin and wheel check, both gates, and the move;
  - MOTCO's three reports as the reason, stating that the 2026-10-26 primaries' p-values were withheld;
  - the corrected primaries in full, with κ;
  - per re-issued step, its regenerated tables and an outcome table "superseded → corrected" (p-values and rejections, with the reproduced values stated equal);
  - the reading above in the outcome that occurred;
  - findings and limitations.

  Numbers are pasted from generated files and checked by a script. Earlier reports stay unchanged.
- **README and AGENTS.md**: the pin source and rule, the constant, the superseded layout, and the `latent_adjusted` reading.

Out of scope:
- MOTCO's `leverage` correction;
- any change to a registry, expectation, adoption rule, threshold or `ANALYSIS_FORCED_RANKS`;
- removing the `latent_adjusted` arm;
- making `df` MOTCO's default (a MOTCO change);
- MOTCO's block weighting;
- `rank-sensitivity` or `preprocessing-sensitivity` on microglia;
- rerunning the MOTCO-independent steps;
- edits to earlier reports.

## Capabilities

### New Capabilities

None.

### Modified Capabilities
- `motco-trajectory-analysis`: MOTCO ≥ 0.13.0. `ANALYSIS_RESIDUAL_CORRECTION = "df"` joins the fixed measurement settings, passed explicitly, checked in the record, recorded in the manifest, `cli_equivalent` and `summary.md`, and checked by every reader of the primary. The report carries MOTCO's statement for the corrected null.
- `motco-rank-sensitivity`: the reference check and the sweep manifest include the residual correction. The forced runs carry it.
- `motco-release-migration`: the pin is an exact PyPI version whose wheel matches the release tag, with no other locked package moving. A measurement change that MOTCO documents as moving only named outputs is contract-checked on the primaries before the move.
- `motco-latent-adjustment-sensitivity`: under MOTCO's `residualized_basis` diagnosis, the `latent_adjusted` arm's `angle` and `shape` p-values are reported but not read. The limitations add that diagnosis and the corrected null's unmeasured size at eight covariates.
- `cell-type-positive-control`: the control is re-evaluated on the corrected primary on both bases. Its outcome is expected to reproduce because no expectation reads a p-value, and it is reported with the superseded outcome beside it.

## Impact

- `pyproject.toml`, `uv.lock`: `motco==0.13.0` from PyPI; the git source is removed.
- `src/sead_sexdiff/data.py`:
  - `ANALYSIS_MOTCO_MIN_VERSION` and `ANALYSIS_RESIDUAL_CORRECTION`;
  - `residual_correction` in the five `TrajectoryAnalysisParams` constructions;
  - the record check and the manifest's `measurement` in `_analyze_cohort`;
  - `_cli_equivalent`, `render_analysis_summary`, `_check_reference` and the rank-sweep manifest.

  No CLI flag is added.
- `README.md`, `AGENTS.md`.
- `results/` (gitignored), on the workstation that runs the re-issue (`nixos-home`, where the current outputs were made) and on the cluster: superseded outputs moved, and every MOTCO-dependent output re-made. The compute equals the 2026-10-26 re-issue: about 50 CV analyses of 1–2 min, 12 forced ones, sixteen variant exports, and eight stage-DE plus eight limma runs; a few hours.
- `reports/`: one new report.
- Scientific choices introduced here (flag and ask before changing any of them):
  - the `df`-corrected null as the primary's for both cell types;
  - `ANALYSIS_RESIDUAL_CORRECTION = "df"`, set away from MOTCO's default by its verdict;
  - the `latent_adjusted` arm's `angle` and `shape` not read;
  - re-issuing every MOTCO-dependent output rather than keeping mixed nulls;
  - the superseded outputs as records, not an arm;
  - the predeclared reading above.
