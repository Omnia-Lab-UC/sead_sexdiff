# MOTCO residual-correction primary: the case study re-issued under MOTCO 0.13.0 with the `df`-corrected RRPP null, astrocytes and microglia — 2026-10-31

> The re-issue of OpenSpec change `mtg-residual-correction-primary` (specs `motco-trajectory-analysis`, `motco-rank-sensitivity`, `motco-release-migration`, `motco-latent-adjustment-sensitivity`, `cell-type-positive-control`). MOTCO is pinned to `motco==0.13.0` from PyPI (the wheel checked byte-identical to tag `v0.13.0`), whose RRPP can rescale the reduced model's residuals before permuting them. `ANALYSIS_RESIDUAL_CORRECTION = "df"` passes MOTCO's validated `df` correction (κ = √(n / (n − r)), r the reduced model's rank) explicitly in every analysis the package runs, checks it in every record and in every reader of a reference analysis, and records it in the manifests, the `cli_equivalent` and `summary.md`. The `df`-corrected analysis is the primary analysis of both cell types. Two same-machine gates came first: without the correction the pin move changes no artifact, and with it only the null moves. Every MOTCO-dependent output of both cell types was then moved to `results/superseded/mtg-residual-correction-primary/` and re-made from scratch: the primary analyses, the rank sweep, the axis diagnostics, the eight preprocessing variants with their adoption check and precedence, the latent-adjustment arms, both stage-DE steps, the donor-exclusion variant with its diagnostics and stage DE, and the positive control on both bases. The superseded values are given beside the corrected ones throughout; the earlier reports stay as dated records.
>

> **Result.** The corrected primaries keep their latent spaces, ranks and observed statistics: rank 5 (astrocyte `cohort_30`, CV mean AUROC 0.6356), rank 2 (astrocyte `cohort_50`, 0.6965), rank 2 (microglia `cohort_30`, 0.6852) and rank 4 (microglia `cohort_50`, 0.6565), with κ 1.0690, 1.0753, 1.0731 and 1.0903. Their p-values rise as the null widens: `delta` 0.468 → 0.501, 0.596 → 0.628, 0.689 → 0.705, 0.819 → 0.834; `angle` 0.396 → 0.449, 0.222 → 0.242, 0.424 → 0.456, 0.227 → 0.291; `shape` 0.284 → 0.331, 0.484 → 0.514, 0.99 → 0.991, 0.598 → 0.644. **No statistic is rejected in any analysis**: not in the four primaries, the twelve forced ranks, the sixteen preprocessing-variant analyses, the eight latent-adjustment arms or the four donor-exclusion analyses. No rejection status differs from its superseded value, so the re-issue shows no dependence on the null correction at α = 0.05. Across the 44 analyses every p-value moved up or stayed (`delta` up in 43, unchanged at 0.99 in one; `angle` and `shape` up in all 44), as a correction with κ > 1 and the same permutation indices predicts. Everything that reads no p-value reproduces: the ranks, CV AUROCs, observed statistics and stage means of all 44 analyses, the first latent components and axis diagnostics, the adoption check (no candidate holds in both cohorts; the precedence names none), the mediator, flag and axis checks, every stage-DE and limma concordance cell (all concordant) and the positive control on both bases: on the corrected primary it **does not hold** in `cohort_30` (E1 holds, E2 holds at 0.6852 against 0.6356, E3 fails at 0.6882 against `mitochondrial_pseudobulk` 0.8438), and the restated control on `latent_adjusted` **holds** in `cohort_30` (E2 0.6575 against 0.6433, E3 0.5279 against `tmm_rna` 0.4153). The re-issue ran on `nixos-dell-laptop`, not on `nixos-home` where the superseded outputs were made: the analyses reproduce at rounding level, the per-feature DESeq2 tables differ beyond rounding with every exported call equal, and one degenerate pydeseq2 fit of a non-exported microglia ATAC peak moves a recorded abundance-set count from 3 to 4 (section 12).


**Reading rules (binding, predeclared in the proposal before any corrected run).** The corrected analysis is the primary of each cell type; the uncorrected floored analyses are superseded records, not a sensitivity arm, and are not cited as results. `cohort_30` stays primary and `cohort_50` sensitivity. The binding rules hold unchanged: the CV-selected rank, the eigengap and the null summary beside `angle`; `angle` and `shape` non-rejections uninformative; a `delta` rejection a difference, not specifically of size. MOTCO's statement is added: under `df` the p-values are size-controlled at the configuration MOTCO measured (six covariates, both `cohort_30` tables), with `delta` conservative there (0.029 / 0.031 pooled), so a `delta` non-rejection stays weak evidence, and the `angle` null is too wide to inform a non-rejection. The correction moves only the null, so every reading that uses no p-value was expected to reproduce; a failed reproduction is a finding with both values and changes no reading. A rejection status that differed in either direction would be dependence on the null correction, with the corrected status the result under the binding rules; nothing is re-tuned. The `latent_adjusted` arm's `angle` and `shape` p-values are reported but not read (MOTCO's `residualized_basis` diagnosis); the arm stays a registered sensitivity analysis. The corrected null's size is inferred, not measured, for `cohort_50`, the eight-covariate arms, the exclusion variants, the forced ranks and the preprocessing variants. Every sensitivity reading is re-read against the corrected primary under its own unchanged rules; the positive control's E1–E3 and all-three rule are unchanged, and its outcome on the corrected primary is the primary-basis outcome, with the superseded outcome beside it.


## 1. Provenance

- Software: python 3.12.15, anndata 0.12.7, numpy 2.3.5, scipy 1.15.3, pandas 2.3.3, h5py 3.15.1, boto3 1.42.12, motco 0.13.0, scikit-learn 1.8.0 (MOTCO from PyPI, `motco==0.13.0`, wheel sha256 `614c2e9cc7fccb5579e4e46c04e5356514ab5bb268b6a821ab89b1a029df83ba` in `uv.lock`); pydeseq2 0.5.3; R 4.6.1 with limma 3.68.4, edgeR 4.10.1 and jsonlite 2.0.0 (the Nix shell; `stage_de_limma.R` unchanged). sead-sexdiff revision `725d5ca27b959b9dd3b33415f45746c1d41e8f7f` for every re-issued output (recorded `dirty: True`: the uncommitted tracked files at run time were the regenerated OpenSpec skill files under `.claude/` and `.junie/`, the IDE file `.idea/sead-sexdiff.iml` and this change's `tasks.md` as its boxes were ticked; nothing under `src/`). `n_jobs = 1` throughout, `--n-cpus 10` for `stage-de`.
- Machine: `nixos-dell-laptop` (Nix shell). The design's re-issue machine, `nixos-home`, which made the superseded outputs, was unreachable on 2026-10-31 (SSH timeout), so the design's fallback applied: the cluster's current directories of both cell types (everything under `results/` except `superseded/` and the per-donor caches, 1,204 files) were copied here with rsync and every file's sha256 verified equal to the cluster's. The gates are unaffected (each compares runs made on this machine); the comparisons with the superseded outputs are across machines (section 12).

### 1.1 The PyPI pin and the wheel check

- `pyproject.toml` lists `motco==0.13.0`; the `[tool.uv.sources]` git entry (`tomszar/MOTCO`, tag `v0.10.0`) is removed. `uv.lock` records motco 0.13.0 from the PyPI registry with the wheel `motco-0.13.0-py3-none-any.whl` (sha256 `614c2e9cc7fccb5579e4e46c04e5356514ab5bb268b6a821ab89b1a029df83ba`).
- The downloaded wheel hashes as the lock records: True. Its `motco/` tree (51 files) equals `git -C ../MOTCO archive v0.13.0 src/motco` (51 files; tag commit `5758c29b9138860509314834b0d57a2f91b81db1`) file by file: True (only in the wheel: none; only in the tag: none; differing: none).
- `uv export --no-hashes` before and after the lock change differs in one line only:

```
100c100
< motco @ git+https://github.com/tomszar/MOTCO@098b170e618959591c7a3d2277487765dcd0459a
---
> motco==0.13.0
```

  scikit-learn stays at 1.8.0, MOTCO's exact pin (`Requires-Dist: scikit-learn==1.8.0`); no other locked package moved. `import motco` prints 0.13.0.

### 1.2 Same-machine pin-move gate (before any corrected run)

The four primary analyses were run on this machine on the real exports twice: from the pre-change commit `fddf295` (a scratch worktree, `uv sync --locked`, MOTCO 0.10.0) with `analyze --cell-type {astrocyte,microglia} --cohort all`, and from this change's code (MOTCO 0.13.0) through `_analyze_cohort` with the `analyze` parameters and `residual_correction = None`. Compared: `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv`, `analysis_metadata.csv` byte for byte; `result.json` and `analysis_manifest.json` after dropping `runtime_seconds`, `created_at`, `git`, `software.motco`, `cli_equivalent` and `measurement.residual_correction`. Every 0.10.0 manifest records MOTCO 0.10.0 and the real exports' sha256; every uncorrected 0.13.0 record carries no `residual_correction` and MOTCO 0.13.0 under `software`.


| cell type / cohort | latent_scores.csv | pls_loadings.csv | stage_means.csv | analysis_metadata.csv | result.json | analysis_manifest.json | gate |
|---|---|---|---|---|---|---|---|
| astrocyte/cohort_30 | identical | identical | identical | identical | equal | equal | pass |
| astrocyte/cohort_50 | identical | identical | identical | identical | equal | equal | pass |
| microglia/cohort_30 | identical | identical | identical | identical | equal | equal | pass |
| microglia/cohort_50 | identical | identical | identical | identical | equal | equal | pass |

Gate: **pass** — without the correction the pin move changes no artifact of the primary analyses.

### 1.3 Contract gate: the correction moves only the null

This change's `analyze --cell-type {astrocyte,microglia} --cohort all` (MOTCO 0.13.0, `df`) was compared on this machine with the uncorrected 0.13.0 runs of 1.2: the four CSVs byte for byte; `result.json` after dropping `p_values`, `null_summary`, `config_spectrum.permutation_pooled_eigengap` (drawn from the permutations), `params.residual_correction`, `measurement.residual_correction` and `measurement.runtime_seconds`; and each κ against √(n / (n − 10)) to 1e-12.


| cell type / cohort | four CSVs | result.json outside the null keys | n | reduced rank | κ recorded | √(n / (n − 10)) | p `delta` none → df | p `angle` none → df | p `shape` none → df | gate |
|---|---|---|---:|---:|---:|---:|---|---|---|---|
| astrocyte/cohort_30 | identical | equal | 80 | 10 | 1.0690449676 | 1.0690449676 | 0.468 → 0.501 | 0.396 → 0.449 | 0.284 → 0.331 | pass |
| astrocyte/cohort_50 | identical | equal | 74 | 10 | 1.0752906584 | 1.0752906584 | 0.596 → 0.628 | 0.222 → 0.242 | 0.484 → 0.514 | pass |
| microglia/cohort_30 | identical | equal | 76 | 10 | 1.0730867400 | 1.0730867400 | 0.689 → 0.705 | 0.424 → 0.456 | 0.99 → 0.991 | pass |
| microglia/cohort_50 | identical | equal | 63 | 10 | 1.0902656765 | 1.0902656765 | 0.819 → 0.834 | 0.227 → 0.291 | 0.598 → 0.644 | pass |

Gate: **pass** — the correction changes no latent score, loading, stage mean, rank, CV AUROC, observed statistic or observed spectrum; it moves only the null, as MOTCO documents. These `df` runs were also the reproducibility check of the re-issued primaries (section 1.6).

### 1.4 The uncorrected rerun against the kept primaries (recorded, not a gate)

- astrocyte/cohort_30 (kept primary made on `nixos-home` at `4eff306`, 2026-10-26T14:46:24+00:00): `latent_scores.csv` max |diff| 2.8e-14; `pls_loadings.csv` max |diff| 3.1e-16; `stage_means.csv` max |diff| 9e-13; `analysis_metadata.csv` differs in the last bit of 3 `log10_n_atac` values; statistics equal to ≤ 4e-13 relative; p-values equal: True, CV AUROC equal: True, CV-selected rank 5 on both.
- astrocyte/cohort_50 (kept primary made on `nixos-home` at `4eff306`, 2026-10-26T14:48:18+00:00): `latent_scores.csv` max |diff| 2.8e-14; `pls_loadings.csv` max |diff| 2e-16; `stage_means.csv` max |diff| 6.1e-13; `analysis_metadata.csv` differs in the last bit of 3 `log10_n_atac` values; statistics equal to ≤ 4e-13 relative; p-values equal: True, CV AUROC equal: True, CV-selected rank 2 on both.
- microglia/cohort_30 (kept primary made on `nixos-home` at `4eff306`, 2026-10-26T14:46:20+00:00): `latent_scores.csv` max |diff| 2.1e-14; `pls_loadings.csv` max |diff| 2e-16; `stage_means.csv` max |diff| 5.1e-13; `analysis_metadata.csv` differs in the last bit of 2 `log10_n_rna`, 6 `log10_n_atac` values; statistics equal to ≤ 4e-13 relative; p-values equal: True, CV AUROC equal: True, CV-selected rank 2 on both.
- microglia/cohort_50 (kept primary made on `nixos-home` at `4eff306`, 2026-10-26T14:48:08+00:00): `latent_scores.csv` max |diff| 1.8e-14; `pls_loadings.csv` max |diff| 2e-16; `stage_means.csv` max |diff| 5.6e-13; `analysis_metadata.csv` differs in the last bit of 3 `log10_n_atac` values; statistics equal to ≤ 4e-13 relative; p-values equal: True, CV AUROC equal: True, CV-selected rank 4 on both.

The kept primaries were made on `nixos-home`; this machine reproduces them at rounding level. Excluding run times, every record value agrees to ≤ 2.6e-11 relative, except the third eigenvalue of each three-point configuration spectrum, which is numerically zero (≤ 1.6e-32) on both machines. `analysis_metadata.csv` differs only in the last bit of a few `log10` values (`np.log10` between machines), which carries through to the scores at rounding level, as in the 2026-10-26 re-issue.

### 1.5 Superseded outputs moved

Before the reruns, `scripts/move_superseded.py` moved 15 directories (1130 files) to `results/superseded/mtg-residual-correction-primary/` with their paths relative to `results/` kept (`motco_analysis` (14), `motco_rank_sensitivity` (90), `motco_diagnostics` (12), `motco_sensitivity` (384), `motco_latent_sensitivity` (72), `stage_de` (46), `stage_de_limma` (56), `motco_exclusion_sensitivity` (128), `microglia/motco_analysis` (14), `microglia/motco_diagnostics` (12), `microglia/motco_latent_sensitivity` (60), `microglia/stage_de` (46), `microglia/stage_de_limma` (56), `microglia/motco_exclusion_sensitivity` (128), `microglia/positive_control` (12); absent here: `microglia/motco_rank_sensitivity`, `microglia/motco_sensitivity`), recorded in `MOVED.json`. A snapshot of every file under `results/` and `data/annotation/` (1,211 files) taken before the move was re-hashed after it: the 1,130 moved files hash as before at their new paths, the 81 MOTCO-independent files are unchanged in place, and no MOTCO-dependent file is left outside `results/superseded/`. The 44 `result.json` sha256 recorded in the moved manifests, including the 12 recorded by the four positive-control manifests whose `positive_control.md` the 2026-10-26 report pasted, resolve under `results/superseded/mtg-residual-correction-primary/`.

### 1.6 Checks after the re-issue

- Version consistency (design D8): outside `results/superseded/`, 156 manifests (86 of which record a MOTCO version) and 44 `result.json` (12 forced, 32 CV-selected) were read; every recorded MOTCO version is 0.13.0, every record carries `params.residual_correction = "df"` and a `df` correction block, and every CV-selected record carries `params.min_components = 2`. No stale manifest.
- MOTCO-independent outputs: the 81 files under `results/pseudobulk/` (33, including six local per-donor cache files), `results/motco_input/` (14), `results/microglia/pseudobulk/` (19), `results/microglia/motco_input/` (14) and `data/annotation/` (1), hashed before the move, hash unchanged after the last rerun.
- Reproducibility: the four re-issued primaries are byte-identical to the contract gate's `df` runs (1.3) in all four CSVs, with equal statistics, p-values, null summaries, ranks and configuration spectra.
- Section-4 checks: every analysis record carries `df` with κ = √(n / (n − r)) (r = 10 at six covariates, 12 at the eight of the latent arms), floor 2, seed 20261019 and 999 permutations; the rank sweep's manifest records the floor and `df` and its reference's `df`; `pls_0` is identical across the primary and the six forced ranks; all eight variants are complete in both cohorts with `adoption_check.csv` for the three candidates and the precedence set in both trend candidates' manifests; `mediator_check.csv` and `flag_check.csv` were written before the arms and the variant exports; every stage-DE and limma manifest binds the re-issued primary's `result.json` sha256 and `delta` status, and the exclusion roots' record the `donor_exclusion` block; every positive-control input resolves to the re-issued file now on disk.
- Cluster (`cluster.ing.uc.cl`, the canonical copy): section 1.7.

### 1.7 Cluster

`cluster.ing.uc.cl` holds the canonical copy. There the branch was checked out at `ecf077e` and `uv sync --locked` replaced the git-pinned MOTCO 0.10.0 with `motco==0.13.0` from PyPI on the login node (`import motco` prints 0.13.0, scikit-learn 1.8.0). The same move (`scripts/move_superseded.py --change mtg-residual-correction-primary`, the directory list of 1.8) moved 15 directories (1,130 files; absent there too: `microglia/motco_rank_sensitivity`, `microglia/motco_sensitivity`), every file hashing as before; the moved tree is identical in files and sha256 to this machine's superseded tree and to the cluster's own hashes taken before its outputs were copied to this machine (section 1), and the 758 files of `pseudobulk/` and `motco_input/` of both cell types hash unchanged. Then the 15 re-issued directories (1,130 files) were copied from `nixos-dell-laptop` with rsync, every file's sha256 equal on both machines; every MOTCO version recorded in a JSON file outside `results/superseded/` on the cluster (172 fields) is 0.13.0.

### 1.8 Commands

```
# pin (design D1)
uv lock                                   # after pyproject: motco==0.13.0, git source removed
# gates (scratch worktree of fddf295, MOTCO 0.10.0; then this change, MOTCO 0.13.0)
sead-sexdiff analyze --cell-type {astrocyte,microglia} --cohort all --input-dir <real export> --out-dir <scratch>/v0100/<cell type>
_analyze_cohort(..., TrajectoryAnalysisParams(..., residual_correction=None))       # <scratch>/v0130_none/<cell type>
sead-sexdiff analyze --cell-type {astrocyte,microglia} --cohort all --input-dir <real export> --out-dir <scratch>/v0130_df/<cell type>
# move
python3 scripts/move_superseded.py --change mtg-residual-correction-primary motco_analysis motco_rank_sensitivity motco_diagnostics \
    motco_sensitivity motco_latent_sensitivity stage_de stage_de_limma motco_exclusion_sensitivity microglia/motco_analysis \
    microglia/motco_rank_sensitivity microglia/motco_diagnostics microglia/motco_sensitivity microglia/motco_latent_sensitivity \
    microglia/stage_de microglia/stage_de_limma microglia/motco_exclusion_sensitivity microglia/positive_control
# re-issue (nix develop -c uv run ...; astrocyte and microglia chains side by side)
sead-sexdiff analyze --cohort all --n-jobs 1
sead-sexdiff analyze --cell-type microglia --cohort all --n-jobs 1
sead-sexdiff rank-sensitivity --cohort all --n-jobs 1
sead-sexdiff axis-diagnostics --cohort all
sead-sexdiff axis-diagnostics --cell-type microglia --cohort all
sead-sexdiff preprocessing-sensitivity --variant all --cohort all --n-jobs 1
V=results/motco_sensitivity/<variant>   # each of the eight
sead-sexdiff axis-diagnostics --cohort all --input-dir $V/motco_input --analysis-dir $V/motco_analysis --rank-dir $V/motco_rank_sensitivity --out-dir $V/diagnostics
sead-sexdiff latent-sensitivity --cohort all --n-jobs 1                       # and --cell-type microglia
sead-sexdiff axis-diagnostics --cohort all --analysis-dir <primary or arm root> --rank-dir results/motco_latent_sensitivity/motco_rank_sensitivity --out-dir results/motco_latent_sensitivity/diagnostics/<name>
sead-sexdiff axis-diagnostics --cell-type microglia --cohort all --analysis-dir <arm root> --out-dir results/microglia/motco_latent_sensitivity/diagnostics/<arm>
sead-sexdiff stage-de --cohort all --omic both --n-cpus 10                    # and --cell-type microglia
sead-sexdiff stage-de-limma --cohort all --omic both                          # and --cell-type microglia
sead-sexdiff exclusion-sensitivity --cohort all --n-jobs 1                    # and --cell-type microglia
V=results/motco_exclusion_sensitivity/excluding_severely_affected            # results/microglia/... with --cell-type microglia
sead-sexdiff axis-diagnostics --cohort all --input-dir $V/motco_input --analysis-dir $V/motco_analysis --rank-dir $V/motco_rank_sensitivity --out-dir $V/diagnostics
sead-sexdiff stage-de --cohort all --omic both --input-dir $V/motco_input --analysis-dir $V/motco_analysis --out-dir $V/stage_de --n-cpus 10
sead-sexdiff stage-de-limma --cohort all --omic rna --input-dir $V/motco_input --analysis-dir $V/motco_analysis --deseq2-dir $V/stage_de --out-dir $V/stage_de_limma
sead-sexdiff positive-control --cohort all
sead-sexdiff positive-control --basis latent_adjusted --cohort all
```

## 2. Why the primaries were re-issued: MOTCO's three reports

The 2026-10-26 re-issue made the floored analyses (MOTCO 0.10.0) the primaries of both cell types. MOTCO then measured `motco analyze` at this case study's configuration and found its RRPP null too narrow there. Three MOTCO reports (`../MOTCO/docs/reports/`) and its roadmap (Phase 6; Next changes items 1 and 8) are the reason for this change; their numbers are quoted here, not re-derived.

1. **Case-study Type I check (2026-10-29, `case-study-type-i-2026-10-29.md`), decision `diagnose_configuration_excess`.** At the astrocyte `cohort_30` design (n = 80, six covariate columns, floor 2, net displacement, 500 replicates × 999 permutations at the anchor and the translation control, stage signal δ = 0.1 and 0.05), against the bound 0.0695 (≤ 34 / 500): the covariate-adjusted analysis held every statistic at δ = 0.05 and five of six at δ = 0.1, exceeding it on the translation control's `delta` (0.072, 36 / 500); pooled over the four readings, `delta` 0.066, `shape` 0.062, `angle` 0.046. The `stage` latent adjustment exceeded the bound at the anchor on `angle` (0.090 / 0.084) and `shape` (0.092 / 0.078). MOTCO's consequence: the primary analysis of this repository was "**withheld pending diagnosis**"; its `delta` p-values were to be read against a realized size of about 0.066 rather than 0.05, and the `stage` arm's `angle` and `shape` rejections were not interpretable.
2. **Diagnosis (2026-10-30, `case-study-diagnosis-2026-10-30.md`), decision `diagnosed`.** The covariate-adjusted `delta` excess is **`residual_df_shrinkage`**: RRPP permutes the reduced model's residuals, whose mean square is RSS/n rather than RSS/(n − p_red), with p_red = 10 at n = 80. Pooled `delta` rose 0.0545 → 0.0595 → 0.066 with one, three and six covariate columns (column contrast one-sided p = 0.015), the residuals rescaled by κ = √(n / (n − p_red)) read 0.0475–0.0545 and held every reading, and the covariate's confounding was not the cause (p = 0.39). The `stage` latent adjustment is **`residualized_basis`**: at a forced rank of 2, with no cross-validation, it still exceeded on anchor `angle` (0.084 / 0.082). MOTCO: "The primary stays withheld until that change is released and validated; then the primary is re-issued under it."
3. **Corrected-permutation check (MOTCO 0.13.0, change `residual-df-corrected-permutation`, `corrected-permutation-2026-10-30.md`), verdict **`adopt_df`**.** 0.13.0 adds the opt-in `residual_correction`; under `df` every residual is multiplied by κ = √(n / (n − r)). At this configuration on both case-study `cohort_30` tables (astrocyte n = 80, κ 1.069; microglia n = 76, κ 1.073), 500 replicates × 999 permutations on a fresh seed family, against 0.0597 pooled (≤ 119 / 2,000) and 0.0792 per reading (≤ 39 / 500): `df` read 0.029–0.058 pooled per table and statistic and at most 32 / 500 per reading (`delta` 0.031 / 0.029, `angle` 0.0425 / 0.05, `shape` 0.058 / 0.052, astrocyte / microglia). The uncorrected null was liberal on pooled `shape` (0.0685 / 0.0675) and on the strong-signal anchor's `angle` (43 / 500 = 0.086), not on `delta` there (0.046 / 0.045). The correction costs 0.014–0.054 of power on the target statistics (pooled magnitude `delta` 0.449 → 0.4135, orientation `angle` 0.41 → 0.382, `shape` 0.188 → 0.1595). `leverage` held too and showed no paired advantage, so the rule chose `df`. MOTCO's consequence: this repository "re-issues both cell types' primaries with `--residual-correction df`, pinned to the release that contains it (0.13.0), in its own change"; "`delta` is conservative under `df` here (0.029 and 0.031 pooled). A `delta` non-rejection therefore stays weak evidence"; the `angle` null is too wide to inform a non-rejection; the `stage` latent-adjusted arm "stays inadmissible at this configuration (diagnosis, `residualized_basis`)". MOTCO keeps `df` opt-in, so it is the one MOTCO default this package overrides, by MOTCO's verdict.

**The 2026-10-26 primaries' p-values were withheld.** From the Type I check (2026-10-29) until this re-issue, MOTCO's roadmap withheld this repository's primary analyses ("The primary therefore stayed withheld until the corrected permutation (item 8) was released and validated"; until the re-issue "the uncorrected primary reads against a liberal null: `delta` about 0.066 on the diagnosis's seed family, and pooled `shape` 0.0675–0.0685 on the check's fresh family"). No earlier report of this repository records the withholding; this report does, and the uncorrected analyses are superseded records from here on. The correction was chosen by MOTCO's predeclared check before any corrected run on these data.


## 3. The corrected primary analyses

Each primary is pasted from its generated `summary.md`, which now states the correction the null was drawn under. The eigengap table gives each analysis's `config_spectrum` relative eigengap and the permutation-pooled eigengap (drawn from the null, so it moves with the correction) beside the `angle` null summary, as the binding rules require. In every primary the observed relative eigengap sits inside its permutation distribution: the floor guarantees two dimensions, not a second axis with signal.

### 3.1 Microglia

#### cohort_30 (≥ 30 microglia nuclei in both omics, n = 76) — primary cohort

Pasted from `results/microglia/motco_analysis/cohort_30/summary.md`:

##### microglia, MOTCO trajectory analysis: cohort_30

MOTCO 0.13.0; n = 76; groups Female vs Male; stages 0 → 1 → 2; orientation `net_displacement`; 999 RRPP permutations, seed 20261019; α = 0.05 (rejected when p ≤ α).

RRPP null: residual correction `df` (κ = 1.07309; n = 76, reduced rank 10).

###### Statistics

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 5.10924 | 0.705 | no | 10.8174 | 8.13722 | 26.9672 |  |
| angle | 12.5272 | 0.456 | no | 13.6061 | 10.1863 | 32.1219 | degrees; CV-selected rank 2 (rank floor 2) |
| shape | 0.0242177 | 0.991 | no | 0.272507 | 0.148066 | 0.544 |  |

###### Cell sizes

| group | stage 0 | stage 1 | stage 2 | total |
|---|---:|---:|---:|---:|
| Female | 11 | 10 | 24 | 45 |
| Male | 9 | 11 | 11 | 31 |

###### Latent space

- PLS-DA rank: 2 (`cv`)
- cv_mean_auroc: 0.6852
- CV: cv2_splits 4, cv1_splits 3, n_repeats 5 (5 completed), max_components 20, rank floor 2, random_state 1203
- blocks: rna (5000 features, identity), atac (5000 features, identity)

###### Covariates

- entered (in order): `age_at_death` (numeric), `pmi` (numeric), `rin` (numeric), `years_education` (numeric), `log10_n_rna` (numeric), `log10_n_atac` (numeric)
- design columns: 6 (age_at_death, pmi, rin, years_education, log10_n_rna, log10_n_atac)
- missing policy: `error`; dropped donors: 0
- LS-means evaluated at the reference row (`pooled_mean`):

| design column | reference value |
|---|---:|
| age_at_death | 88.8684 |
| pmi | 6.95281 |
| rin | 8.335 |
| years_education | 16.2632 |
| log10_n_rna | 2.5772 |
| log10_n_atac | 1.98488 |

#### cohort_50 (≥ 50 microglia nuclei in both omics, n = 63) — sensitivity cohort

Pasted from `results/microglia/motco_analysis/cohort_50/summary.md`:

##### microglia, MOTCO trajectory analysis: cohort_50

MOTCO 0.13.0; n = 63; groups Female vs Male; stages 0 → 1 → 2; orientation `net_displacement`; 999 RRPP permutations, seed 20261019; α = 0.05 (rejected when p ≤ α).

RRPP null: residual correction `df` (κ = 1.09027; n = 63, reduced rank 10).

###### Statistics

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 1.6686 | 0.834 | no | 6.37019 | 4.88809 | 15.3753 |  |
| angle | 22.2604 | 0.291 | no | 18.0464 | 8.60842 | 34.1958 | degrees; CV-selected rank 4 (rank floor 2) |
| shape | 0.0797433 | 0.644 | no | 0.106788 | 0.0560953 | 0.21359 |  |

###### Cell sizes

| group | stage 0 | stage 1 | stage 2 | total |
|---|---:|---:|---:|---:|
| Female | 8 | 10 | 18 | 36 |
| Male | 8 | 10 | 9 | 27 |

###### Latent space

- PLS-DA rank: 4 (`cv`)
- cv_mean_auroc: 0.6565
- CV: cv2_splits 4, cv1_splits 3, n_repeats 5 (5 completed), max_components 20, rank floor 2, random_state 1203
- blocks: rna (5000 features, identity), atac (5000 features, identity)

###### Covariates

- entered (in order): `age_at_death` (numeric), `pmi` (numeric), `rin` (numeric), `years_education` (numeric), `log10_n_rna` (numeric), `log10_n_atac` (numeric)
- design columns: 6 (age_at_death, pmi, rin, years_education, log10_n_rna, log10_n_atac)
- missing policy: `error`; dropped donors: 0
- LS-means evaluated at the reference row (`pooled_mean`):

| design column | reference value |
|---|---:|
| age_at_death | 89 |
| pmi | 7.05048 |
| rin | 8.3927 |
| years_education | 16.2063 |
| log10_n_rna | 2.58304 |
| log10_n_atac | 2.06888 |

#### Eigengap and `angle` null (microglia, from `result.json` `config_spectrum` and `null_summary`)

| cohort | rank | κ | relative eigengap (pooled) | permutation-pooled eigengap superseded: mean, q05, q95 | corrected: mean, q05, q95 | `angle` | `angle` null mean, sd, q95 superseded | corrected |
|---|---:|---:|---:|---|---|---:|---|---|
| cohort_30 | 2 | 1.0731 | 0.7753 | 0.7720, 0.5003, 0.9770 | 0.7706, 0.4795, 0.9824 | 12.5272 | 12.6966, 9.5119, 30.1659 | 13.6061, 10.1863, 32.1219 |
| cohort_50 | 4 | 1.0903 | 0.2248 | 0.2328, 0.1237, 0.3448 | 0.2337, 0.1146, 0.3544 | 22.2604 | 16.6216, 7.9428, 31.3913 | 18.0464, 8.6084, 34.1958 |

### 3.2 Astrocyte

#### cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

Pasted from `results/motco_analysis/cohort_30/summary.md`:

##### astrocyte, MOTCO trajectory analysis: cohort_30

MOTCO 0.13.0; n = 80; groups Female vs Male; stages 0 → 1 → 2; orientation `net_displacement`; 999 RRPP permutations, seed 20261019; α = 0.05 (rejected when p ≤ α).

RRPP null: residual correction `df` (κ = 1.06904; n = 80, reduced rank 10).

###### Statistics

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |  |
| angle | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 | degrees; CV-selected rank 5 (rank floor 2) |
| shape | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |  |

###### Cell sizes

| group | stage 0 | stage 1 | stage 2 | total |
|---|---:|---:|---:|---:|
| Female | 11 | 10 | 28 | 49 |
| Male | 9 | 10 | 12 | 31 |

###### Latent space

- PLS-DA rank: 5 (`cv`)
- cv_mean_auroc: 0.6356
- CV: cv2_splits 4, cv1_splits 3, n_repeats 5 (5 completed), max_components 20, rank floor 2, random_state 1203
- blocks: rna (5000 features, identity), atac (5000 features, identity)

###### Covariates

- entered (in order): `age_at_death` (numeric), `pmi` (numeric), `rin` (numeric), `years_education` (numeric), `log10_n_rna` (numeric), `log10_n_atac` (numeric)
- design columns: 6 (age_at_death, pmi, rin, years_education, log10_n_rna, log10_n_atac)
- missing policy: `error`; dropped donors: 0
- LS-means evaluated at the reference row (`pooled_mean`):

| design column | reference value |
|---|---:|
| age_at_death | 89.025 |
| pmi | 7.0935 |
| rin | 8.15087 |
| years_education | 16.1375 |
| log10_n_rna | 2.82437 |
| log10_n_atac | 2.22423 |

#### cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74) — sensitivity cohort

Pasted from `results/motco_analysis/cohort_50/summary.md`:

##### astrocyte, MOTCO trajectory analysis: cohort_50

MOTCO 0.13.0; n = 74; groups Female vs Male; stages 0 → 1 → 2; orientation `net_displacement`; 999 RRPP permutations, seed 20261019; α = 0.05 (rejected when p ≤ α).

RRPP null: residual correction `df` (κ = 1.07529; n = 74, reduced rank 10).

###### Statistics

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |  |
| angle | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 | degrees; CV-selected rank 2 (rank floor 2) |
| shape | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |  |

###### Cell sizes

| group | stage 0 | stage 1 | stage 2 | total |
|---|---:|---:|---:|---:|
| Female | 10 | 9 | 25 | 44 |
| Male | 9 | 9 | 12 | 30 |

###### Latent space

- PLS-DA rank: 2 (`cv`)
- cv_mean_auroc: 0.6965
- CV: cv2_splits 4, cv1_splits 3, n_repeats 5 (5 completed), max_components 20, rank floor 2, random_state 1203
- blocks: rna (5000 features, identity), atac (5000 features, identity)

###### Covariates

- entered (in order): `age_at_death` (numeric), `pmi` (numeric), `rin` (numeric), `years_education` (numeric), `log10_n_rna` (numeric), `log10_n_atac` (numeric)
- design columns: 6 (age_at_death, pmi, rin, years_education, log10_n_rna, log10_n_atac)
- missing policy: `error`; dropped donors: 0
- LS-means evaluated at the reference row (`pooled_mean`):

| design column | reference value |
|---|---:|
| age_at_death | 89.3108 |
| pmi | 7.23847 |
| rin | 8.14784 |
| years_education | 16 |
| log10_n_rna | 2.832 |
| log10_n_atac | 2.2729 |

#### Eigengap and `angle` null (astrocyte, from `result.json` `config_spectrum` and `null_summary`)

| cohort | rank | κ | relative eigengap (pooled) | permutation-pooled eigengap superseded: mean, q05, q95 | corrected: mean, q05, q95 | `angle` | `angle` null mean, sd, q95 superseded | corrected |
|---|---:|---:|---:|---|---|---:|---|---|
| cohort_30 | 5 | 1.0690 | 0.1709 | 0.2143, 0.0792, 0.3615 | 0.2164, 0.0779, 0.3668 | 21.9941 | 20.7155, 9.2780, 37.7882 | 22.0577, 9.8721, 40.3393 |
| cohort_50 | 2 | 1.0753 | 0.9700 | 0.9275, 0.7444, 0.9997 | 0.9195, 0.7162, 0.9996 | 23.1596 | 15.1566, 12.4608, 40.0109 | 16.3255, 13.4558, 43.0401 |

### 3.3 Superseded → corrected (primaries)

| cell type | cohort | analysis | rank | CV AUROC | κ (n, reduced rank) | statistic | observed | p superseded → corrected | change | rejected superseded → corrected |
|---|---|---|---:|---:|---|---|---:|---|---|---|
| astrocyte | cohort_30 | primary | 5 | 0.6356 | 1.0690 (80, 10) | delta | 5.95773 | 0.468 → 0.501 | up | no → no |
| astrocyte | cohort_30 | primary |  |  |  | angle | 21.9941 | 0.396 → 0.449 | up | no → no |
| astrocyte | cohort_30 | primary |  |  |  | shape | 0.15645 | 0.284 → 0.331 | up | no → no |
| astrocyte | cohort_50 | primary | 2 | 0.6965 | 1.0753 (74, 10) | delta | 6.01569 | 0.596 → 0.628 | up | no → no |
| astrocyte | cohort_50 | primary |  |  |  | angle | 23.1596 | 0.222 → 0.242 | up | no → no |
| astrocyte | cohort_50 | primary |  |  |  | shape | 0.284932 | 0.484 → 0.514 | up | no → no |
| microglia | cohort_30 | primary | 2 | 0.6852 | 1.0731 (76, 10) | delta | 5.10924 | 0.689 → 0.705 | up | no → no |
| microglia | cohort_30 | primary |  |  |  | angle | 12.5272 | 0.424 → 0.456 | up | no → no |
| microglia | cohort_30 | primary |  |  |  | shape | 0.0242177 | 0.99 → 0.991 | up | no → no |
| microglia | cohort_50 | primary | 4 | 0.6565 | 1.0903 (63, 10) | delta | 1.6686 | 0.819 → 0.834 | up | no → no |
| microglia | cohort_50 | primary |  |  |  | angle | 22.2604 | 0.227 → 0.291 | up | no → no |
| microglia | cohort_50 | primary |  |  |  | shape | 0.0797433 | 0.598 → 0.644 | up | no → no |

Superseded against corrected, generated by the comparison script over 4 analyses: of the 16 analysis CSVs (four per analysis), 16 equal at rounding level (largest |difference| 9e-13); of the 4 records outside the null keys, 4 equal at rounding level (largest relative difference 2.6e-11, apart from the numerically-zero third eigenvalues of the configuration spectra, ≤ 1.6e-32); rank equal in 4 of 4, CV AUROC equal in 4 of 4; κ = √(n / (n − r)) to 1e-12 in 4 of 4; p-values up in 12, equal in 0, down in 0 of 12; rejection statuses changed: 0.

## 4. Rank sensitivity (astrocytes)

The sweep was rerun against the corrected primaries with the correction at every forced rank (r = 10, the same κ as its cohort's primary). No forced rank rejects any statistic in either cohort, and no rejection status differs from the superseded sweep: no rank dependence and no dependence on the null correction. The CV-selected rank 5 of astrocyte `cohort_30` is outside the ladder and the reference row stands in for it; in `cohort_50` the forced rank-2 run equals the CV-selected reference in statistics and p-values. MOTCO's forced-rank operating characteristics stay limitations, and the corrected null's size at a forced rank is inferred, not measured.

### 4.1 cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

#### MOTCO rank sensitivity: cohort_30

Forced-rank sweep of the covariate-adjusted trajectory analysis. The reference row (`cv`) is the CV-selected analysis at rank 5, which stays the primary result; the forced ranks 2, 3, 4, 6, 9, 12 are read as sensitivity to it, and no rank is chosen from the results. 999 RRPP permutations with seed 20261019 at every rank (the same permutation indices throughout); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative at every rank; a `delta` rejection is a difference, not specifically of size.

##### delta

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| forced | 2 | 3.91847 | 0.764 | no | 10.323 | 7.75591 | 24.9881 |
| forced | 3 | 10.2517 | 0.379 | no | 9.15613 | 7.3164 | 23.2727 |
| forced | 4 | 19.247 | 0.079 | no | 8.42977 | 6.54143 | 21.1146 |
| forced | 6 | 6.08738 | 0.41 | no | 6.2179 | 5.05874 | 16.2431 |
| forced | 9 | 1.79797 | 0.792 | no | 5.77816 | 4.75609 | 14.9617 |
| forced | 12 | 1.15591 | 0.879 | no | 5.85988 | 4.74306 | 15.0591 |

##### angle (degrees)

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| forced | 2 | 22.832 | 0.271 | no | 16.8725 | 13.578 | 43.6636 |
| forced | 3 | 23.0297 | 0.357 | no | 20.0799 | 11.2524 | 41.674 |
| forced | 4 | 22.1423 | 0.403 | no | 20.9569 | 10.3038 | 39.8287 |
| forced | 6 | 21.9664 | 0.518 | no | 23.4375 | 9.52376 | 41.4501 |
| forced | 9 | 31.5192 | 0.28 | no | 27.0358 | 8.58543 | 42.8006 |
| forced | 12 | 35.8438 | 0.211 | no | 30.15 | 7.86997 | 44.1979 |

##### shape

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| forced | 2 | 0.295865 | 0.471 | no | 0.305979 | 0.165745 | 0.615801 |
| forced | 3 | 0.121943 | 0.812 | no | 0.232822 | 0.12249 | 0.448958 |
| forced | 4 | 0.0296112 | 0.977 | no | 0.172843 | 0.0873044 | 0.327695 |
| forced | 6 | 0.153228 | 0.266 | no | 0.118237 | 0.0617083 | 0.228477 |
| forced | 9 | 0.147921 | 0.201 | no | 0.101946 | 0.0542056 | 0.199609 |
| forced | 12 | 0.131817 | 0.246 | no | 0.0983606 | 0.052514 | 0.192181 |

##### Runs

| selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|
| cv | 5 | 0.6356 | 113.8 | `results/motco_analysis/cohort_30` |
| forced | 2 | n/a | 0.4 | `rank_02` |
| forced | 3 | n/a | 0.4 | `rank_03` |
| forced | 4 | n/a | 0.8 | `rank_04` |
| forced | 6 | n/a | 0.8 | `rank_06` |
| forced | 9 | n/a | 1.7 | `rank_09` |
| forced | 12 | n/a | 2.3 | `rank_12` |

### 4.2 cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74) — sensitivity cohort

#### MOTCO rank sensitivity: cohort_50

Forced-rank sweep of the covariate-adjusted trajectory analysis. The reference row (`cv`) is the CV-selected analysis at rank 2, which stays the primary result; the forced ranks 2, 3, 4, 6, 9, 12 are read as sensitivity to it, and no rank is chosen from the results. 999 RRPP permutations with seed 20261019 at every rank (the same permutation indices throughout); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative at every rank; a `delta` rejection is a difference, not specifically of size.

##### delta

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| forced | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| forced | 3 | 9.36583 | 0.527 | no | 11.9926 | 9.26716 | 29.3778 |
| forced | 4 | 14.1822 | 0.248 | no | 9.77049 | 7.1315 | 23.309 |
| forced | 6 | 10.7046 | 0.245 | no | 7.46128 | 5.73876 | 18.4253 |
| forced | 9 | 6.51697 | 0.42 | no | 6.80668 | 5.46153 | 17.6321 |
| forced | 12 | 5.88663 | 0.473 | no | 6.83839 | 5.48645 | 17.4329 |

##### angle (degrees)

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| forced | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| forced | 3 | 23.4732 | 0.337 | no | 20.6913 | 12.448 | 45.5826 |
| forced | 4 | 24.5502 | 0.348 | no | 21.9741 | 11.4489 | 44.8719 |
| forced | 6 | 25.2911 | 0.403 | no | 24.1728 | 10.1042 | 43.9033 |
| forced | 9 | 30.5012 | 0.319 | no | 27.5517 | 9.16436 | 45.6858 |
| forced | 12 | 32.5938 | 0.349 | no | 30.7019 | 8.52074 | 46.9294 |

##### shape

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| forced | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| forced | 3 | 0.16487 | 0.754 | no | 0.266607 | 0.133073 | 0.507623 |
| forced | 4 | 0.110409 | 0.75 | no | 0.17881 | 0.0945905 | 0.360001 |
| forced | 6 | 0.0845199 | 0.685 | no | 0.124553 | 0.0671707 | 0.246623 |
| forced | 9 | 0.0826203 | 0.628 | no | 0.110816 | 0.0605131 | 0.220582 |
| forced | 12 | 0.0785119 | 0.641 | no | 0.106864 | 0.0590265 | 0.213553 |

##### Runs

| selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|
| cv | 2 | 0.6965 | 102.3 | `results/motco_analysis/cohort_50` |
| forced | 2 | n/a | 0.4 | `rank_02` |
| forced | 3 | n/a | 0.7 | `rank_03` |
| forced | 4 | n/a | 1.0 | `rank_04` |
| forced | 6 | n/a | 1.0 | `rank_06` |
| forced | 9 | n/a | 0.5 | `rank_09` |
| forced | 12 | n/a | 1.6 | `rank_12` |

### 4.3 Superseded → corrected (forced ranks)

| cell type | cohort | analysis | rank | CV AUROC | κ (n, reduced rank) | statistic | observed | p superseded → corrected | change | rejected superseded → corrected |
|---|---|---|---:|---:|---|---|---:|---|---|---|
| astrocyte | cohort_30 | rank_02 | 2 | — (forced) | 1.0690 (80, 10) | delta | 3.91847 | 0.744 → 0.764 | up | no → no |
| astrocyte | cohort_30 | rank_02 |  |  |  | angle | 22.832 | 0.241 → 0.271 | up | no → no |
| astrocyte | cohort_30 | rank_02 |  |  |  | shape | 0.295865 | 0.445 → 0.471 | up | no → no |
| astrocyte | cohort_30 | rank_03 | 3 | — (forced) | 1.0690 (80, 10) | delta | 10.2517 | 0.338 → 0.379 | up | no → no |
| astrocyte | cohort_30 | rank_03 |  |  |  | angle | 23.0297 | 0.316 → 0.357 | up | no → no |
| astrocyte | cohort_30 | rank_03 |  |  |  | shape | 0.121943 | 0.789 → 0.812 | up | no → no |
| astrocyte | cohort_30 | rank_04 | 4 | — (forced) | 1.0690 (80, 10) | delta | 19.247 | 0.052 → 0.079 | up | no → no |
| astrocyte | cohort_30 | rank_04 |  |  |  | angle | 22.1423 | 0.356 → 0.403 | up | no → no |
| astrocyte | cohort_30 | rank_04 |  |  |  | shape | 0.0296112 | 0.975 → 0.977 | up | no → no |
| astrocyte | cohort_30 | rank_06 | 6 | — (forced) | 1.0690 (80, 10) | delta | 6.08738 | 0.369 → 0.41 | up | no → no |
| astrocyte | cohort_30 | rank_06 |  |  |  | angle | 21.9664 | 0.455 → 0.518 | up | no → no |
| astrocyte | cohort_30 | rank_06 |  |  |  | shape | 0.153228 | 0.224 → 0.266 | up | no → no |
| astrocyte | cohort_30 | rank_09 | 9 | — (forced) | 1.0690 (80, 10) | delta | 1.79797 | 0.778 → 0.792 | up | no → no |
| astrocyte | cohort_30 | rank_09 |  |  |  | angle | 31.5192 | 0.21 → 0.28 | up | no → no |
| astrocyte | cohort_30 | rank_09 |  |  |  | shape | 0.147921 | 0.161 → 0.201 | up | no → no |
| astrocyte | cohort_30 | rank_12 | 12 | — (forced) | 1.0690 (80, 10) | delta | 1.15591 | 0.866 → 0.879 | up | no → no |
| astrocyte | cohort_30 | rank_12 |  |  |  | angle | 35.8438 | 0.152 → 0.211 | up | no → no |
| astrocyte | cohort_30 | rank_12 |  |  |  | shape | 0.131817 | 0.208 → 0.246 | up | no → no |
| astrocyte | cohort_50 | rank_02 | 2 | — (forced) | 1.0753 (74, 10) | delta | 6.01569 | 0.596 → 0.628 | up | no → no |
| astrocyte | cohort_50 | rank_02 |  |  |  | angle | 23.1596 | 0.222 → 0.242 | up | no → no |
| astrocyte | cohort_50 | rank_02 |  |  |  | shape | 0.284932 | 0.484 → 0.514 | up | no → no |
| astrocyte | cohort_50 | rank_03 | 3 | — (forced) | 1.0753 (74, 10) | delta | 9.36583 | 0.496 → 0.527 | up | no → no |
| astrocyte | cohort_50 | rank_03 |  |  |  | angle | 23.4732 | 0.299 → 0.337 | up | no → no |
| astrocyte | cohort_50 | rank_03 |  |  |  | shape | 0.16487 | 0.739 → 0.754 | up | no → no |
| astrocyte | cohort_50 | rank_04 | 4 | — (forced) | 1.0753 (74, 10) | delta | 14.1822 | 0.217 → 0.248 | up | no → no |
| astrocyte | cohort_50 | rank_04 |  |  |  | angle | 24.5502 | 0.292 → 0.348 | up | no → no |
| astrocyte | cohort_50 | rank_04 |  |  |  | shape | 0.110409 | 0.719 → 0.75 | up | no → no |
| astrocyte | cohort_50 | rank_06 | 6 | — (forced) | 1.0753 (74, 10) | delta | 10.7046 | 0.208 → 0.245 | up | no → no |
| astrocyte | cohort_50 | rank_06 |  |  |  | angle | 25.2911 | 0.33 → 0.403 | up | no → no |
| astrocyte | cohort_50 | rank_06 |  |  |  | shape | 0.0845199 | 0.652 → 0.685 | up | no → no |
| astrocyte | cohort_50 | rank_09 | 9 | — (forced) | 1.0753 (74, 10) | delta | 6.51697 | 0.386 → 0.42 | up | no → no |
| astrocyte | cohort_50 | rank_09 |  |  |  | angle | 30.5012 | 0.251 → 0.319 | up | no → no |
| astrocyte | cohort_50 | rank_09 |  |  |  | shape | 0.0826203 | 0.593 → 0.628 | up | no → no |
| astrocyte | cohort_50 | rank_12 | 12 | — (forced) | 1.0753 (74, 10) | delta | 5.88663 | 0.436 → 0.473 | up | no → no |
| astrocyte | cohort_50 | rank_12 |  |  |  | angle | 32.5938 | 0.279 → 0.349 | up | no → no |
| astrocyte | cohort_50 | rank_12 |  |  |  | shape | 0.0785119 | 0.598 → 0.641 | up | no → no |

Superseded against corrected, generated by the comparison script over 12 analyses: of the 48 analysis CSVs (four per analysis), 48 equal at rounding level (largest |difference| 9e-13); of the 12 records outside the null keys, 12 equal at rounding level (largest relative difference 4.3e-11, apart from the numerically-zero third eigenvalues of the configuration spectra, ≤ 3.6e-32); rank equal in 12 of 12, CV AUROC equal in 12 of 12; κ = √(n / (n − r)) to 1e-12 in 12 of 12; p-values up in 36, equal in 0, down in 0 of 36; rejection statuses changed: 0.

## 5. Axis diagnostics: the first latent component

`axis-diagnostics` was rerun on every re-issued analysis root (primary with the sweep, latent arms, preprocessing variants, exclusion variant). Its input-hash checks passed everywhere, and the astrocyte diagnostics record `pls_0` as identical across the primary and the six forced ranks. The correction moves no latent score, so every correlation reproduces: the diagnostics' CSVs are byte-identical or equal at rounding level to the superseded ones (largest |difference| 6.0e-10, relative 4.2e-11, across machines). The first-axis readings are those of the earlier reports: in microglia `pls_0` tracks `mitochondrial_pseudobulk` more strongly than stage, in astrocytes the ribosomal/library-composition axis. An axis correlated with a technical variable stays a measurement-space finding, not a reinterpretation of the statistics.

Generated from each `axis_correlations.csv` (rows `analysis = primary`, `component = pls_0`; technical variables: the technical table and the export's TMM factors and log10 library sizes), corrected; the last column compares with the superseded diagnostic (to 1e-9):

| cell type | cohort | analysis | rank | stage ρ | largest technical |ρ| | technical variable | n flagged (pls_0, raw or partial) | same as superseded |
|---|---|---|---|---|---|---|---|---|
| microglia | cohort_30 | primary | 2 | 0.688175 | 0.843828 | mitochondrial_pseudobulk | 28 | yes |
| microglia | cohort_30 | extended_covariates | 2 | 0.688175 | 0.843828 | mitochondrial_pseudobulk | 28 | yes |
| microglia | cohort_30 | latent_adjusted | 3 | -0.527877 | 0.415256 | tmm_rna | 10 | yes |
| microglia | cohort_30 | excluding_severely_affected | 4 | 0.820134 | 0.689149 | fraction_mitochondrial_umis_sd | 18 | yes |
| microglia | cohort_50 | primary | 4 | 0.730773 | 0.784562 | mitochondrial_pseudobulk | 25 | yes |
| microglia | cohort_50 | extended_covariates | 4 | 0.730773 | 0.784562 | mitochondrial_pseudobulk | 25 | yes |
| microglia | cohort_50 | latent_adjusted | 3 | -0.532681 | 0.406702 | gex_reads_mapped_confidently_to_intronic_regions_sd | 8 | yes |
| microglia | cohort_50 | excluding_severely_affected | 4 | 0.824801 | 0.674613 | fraction_mitochondrial_umis_sd | 22 | yes |
| astrocyte | cohort_30 | primary | 5 | 0.491212 | 0.840858 | ribosomal_mean | 20 | yes |
| astrocyte | cohort_30 | extended_covariates | 5 | 0.491212 | 0.840858 | ribosomal_mean | 20 | yes |
| astrocyte | cohort_30 | latent_adjusted | 4 | -0.636917 | 0.382068 | log10_lib_atac | 3 | yes |
| astrocyte | cohort_30 | excluding_severely_affected | 2 | 0.578399 | 0.728598 | ribosomal_mean | 18 | yes |
| astrocyte | cohort_50 | primary | 2 | 0.503605 | 0.839822 | ribosomal_mean | 21 | yes |
| astrocyte | cohort_50 | extended_covariates | 2 | 0.503605 | 0.839822 | ribosomal_mean | 21 | yes |
| astrocyte | cohort_50 | latent_adjusted | 3 | -0.752715 | 0.326442 | log10_lib_atac | 2 | yes |
| astrocyte | cohort_50 | excluding_severely_affected | 2 | 0.595467 | 0.671612 | ribosomal_mean | 22 | yes |

## 6. Preprocessing sensitivity (astrocytes)

All eight registered variants were rerun in both cohorts (the export at the primary's thresholds and N, then the analysis at the primary's settings, the floor and the correction). The variant exports reproduce the superseded ones: four variants and both exclusion roots byte for byte, and the four that residualize or trend-normalize (`residualized`, `revised_candidate`, `trend_normalized`, `revised_candidate_trend`) to the last printed digit (|difference| ≤ 1e-8, the same features selected), across machines. No variant rejects any statistic in either cohort, and no rejection status differs from the corrected primary or from the superseded variant: no preprocessing dependence and no dependence on the null correction at α = 0.05. The variants' CV-selected ranks (2 to 7) differ from the primary's in several rows; each is reported beside `angle`, not read as better or worse.

The adoption rule reads no p-value, so it reproduces exactly: `adoption_check.csv` is byte-identical in all six candidate cells. `revised_candidate` fails in `cohort_30` (`pls_1` × `tmm_rna` 0.377, `pls_1` × `tmm_atac` −0.400, `pls_2` × `tmm_atac` 0.366) and in `cohort_50`; `trend_normalized` fails in both cohorts on `pls_0`; `revised_candidate_trend` fails in `cohort_30` (`pls_1` × `tmm_rna` −0.323) and holds in `cohort_50`. The precedence names no candidate, as before: no revised-primary change follows and the variants stay sensitivity analyses. MOTCO's operating characteristics were measured on none of the variant blocks; the corrected null's size there is inferred.

### 6.1 cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

#### Variant `keep_chrx`

##### MOTCO preprocessing sensitivity: cohort_30, variant `keep_chrx`

chrX kept; chrY and chrM removed; every other export setting equal to the primary. Export setting(s) changed: `excluded_chromosomes` = ['chrY', 'chrM'] (primary: `excluded_chromosomes` = ['chrX', 'chrY', 'chrM']); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `keep_chrx` row (CV-selected rank 4) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| keep_chrx | cv | 4 | 19.0474 | 0.083 | no | 8.46697 | 6.57309 | 21.188 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| keep_chrx | cv | 4 | 22.3677 | 0.396 | no | 21.0629 | 10.3721 | 40.2918 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| keep_chrx | cv | 4 | 0.0237388 | 0.991 | no | 0.174101 | 0.087966 | 0.33274 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports; `kept on chrX` is per export.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept on chrX | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20387 | 5000 | 0 | 4875 | 125 |
| rna | keep_chrx | 36601 | 36477 | 20995 | 5000 | 125 | 4875 | 125 |
| atac | primary | 218882 | 212816 | 128916 | 5000 | 0 | 4898 | 102 |
| atac | keep_chrx | 218882 | 218702 | 130697 | 5000 | 102 | 4898 | 102 |

- `XIST` (chrX, rna): status `sex_chromosome` in the primary export, `kept` in the keep_chrx export.

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 5 | 0.6356 | 113.8 | `results/motco_analysis/cohort_30` |
| keep_chrx | cv | 4 | 0.6341 | 110.0 | `motco_analysis/cohort_30` |

###### Findings

- CV-selected rank differs: 5 in the primary analysis, 4 in the keep_chrx analysis; reported beside angle, not read as better or worse
- keep_chrx export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `exclude_families`

##### MOTCO preprocessing sensitivity: cohort_30, variant `exclude_families`

RNA genes of the ribosomal, mitochondrial and hemoglobin sets (TECHNICAL_GENE_SETS) removed before the abundance filter; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'] (primary: `excluded_gene_sets` = []); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `exclude_families` row (CV-selected rank 5) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| exclude_families | cv | 5 | 6.05273 | 0.493 | no | 7.38995 | 5.76405 | 18.6499 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| exclude_families | cv | 5 | 21.9731 | 0.45 | no | 22.0566 | 9.86138 | 40.2619 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| exclude_families | cv | 5 | 0.156521 | 0.328 | no | 0.13198 | 0.0683819 | 0.260052 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20387 | 5000 | 4990 | 10 |
| rna | exclude_families | 36601 | 35329 | 35249 | 20317 | 5000 | 4990 | 10 |
| atac | primary | 218882 | 212816 | 212816 | 128916 | 5000 | 5000 | 0 |
| atac | exclude_families | 218882 | 212816 | 212816 | 128916 | 5000 | 5000 | 0 |

###### Gene sets removed

RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps `sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.

| omic | set | matched | removed here | already on an excluded chromosome |
|---|---|---:|---:|---:|
| rna | ribosomal | 81 | 77 | 4 |
| rna | mitochondrial | 13 | 0 | 13 |
| rna | hemoglobin | 3 | 3 | 0 |
| atac | ribosomal | 0 | 0 | 0 |
| atac | mitochondrial | 0 | 0 | 0 |
| atac | hemoglobin | 0 | 0 | 0 |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 5 | 0.6356 | 113.8 | `results/motco_analysis/cohort_30` |
| exclude_families | cv | 5 | 0.6369 | 109.3 | `motco_analysis/cohort_30` |

###### Findings

- exclude_families export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `rna_trend_ratio`

##### MOTCO preprocessing sensitivity: cohort_30, variant `rna_trend_ratio`

RNA features ranked at the selection step by the ratio of their log2-CPM variance to a binned-median mean–variance trend instead of by raw variance; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'trend_ratio', 'atac': 'variance'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `rna_trend_ratio` row (CV-selected rank 6) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| rna_trend_ratio | cv | 6 | 5.89665 | 0.48 | no | 6.89413 | 5.51164 | 17.8636 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| rna_trend_ratio | cv | 6 | 22.2081 | 0.53 | no | 23.755 | 9.20883 | 40.2582 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| rna_trend_ratio | cv | 6 | 0.144967 | 0.316 | no | 0.12014 | 0.0625566 | 0.23536 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20387 | 5000 | 2536 | 2464 |
| rna | rna_trend_ratio | 36601 | 35329 | 20387 | 5000 | 2536 | 2464 |
| atac | primary | 218882 | 212816 | 128916 | 5000 | 5000 | 0 |
| atac | rna_trend_ratio | 218882 | 212816 | 128916 | 5000 | 5000 | 0 |

###### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | primary | variance | -0.721 | -0.137 | 0.774 | 2.160 | 4.968 |
| rna | rna_trend_ratio | trend_ratio | 0.084 | 2.086 | 3.819 | 5.362 | 7.521 |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 5 | 0.6356 | 113.8 | `results/motco_analysis/cohort_30` |
| rna_trend_ratio | cv | 6 | 0.6351 | 138.9 | `motco_analysis/cohort_30` |

###### Findings

- CV-selected rank differs: 5 in the primary analysis, 6 in the rna_trend_ratio analysis; reported beside angle, not read as better or worse
- rna_trend_ratio export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `atac_mean_accessibility`

##### MOTCO preprocessing sensitivity: cohort_30, variant `atac_mean_accessibility`

ATAC peaks ranked at the selection step by mean log2-CPM (total accessibility) instead of by variance; the RNA block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'variance', 'atac': 'mean'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `atac_mean_accessibility` row (CV-selected rank 7) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| atac_mean_accessibility | cv | 7 | 0.18785 | 0.99 | no | 11.5291 | 9.63552 | 31.2613 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| atac_mean_accessibility | cv | 7 | 18.6469 | 0.814 | no | 30.3014 | 13.5315 | 57.3551 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| atac_mean_accessibility | cv | 7 | 0.102546 | 0.819 | no | 0.200624 | 0.10309 | 0.396866 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20387 | 5000 | 5000 | 0 |
| rna | atac_mean_accessibility | 36601 | 35329 | 20387 | 5000 | 5000 | 0 |
| atac | primary | 218882 | 212816 | 128916 | 5000 | 0 | 5000 |
| atac | atac_mean_accessibility | 218882 | 212816 | 128916 | 5000 | 0 | 5000 |

###### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| atac | primary | variance | -0.321 | 0.046 | 0.371 | 0.778 | 1.531 |
| atac | atac_mean_accessibility | mean | 5.436 | 5.613 | 5.876 | 6.274 | 7.156 |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 5 | 0.6356 | 113.8 | `results/motco_analysis/cohort_30` |
| atac_mean_accessibility | cv | 7 | 0.6706 | 124.7 | `motco_analysis/cohort_30` |

###### Findings

- CV-selected rank differs: 5 in the primary analysis, 7 in the atac_mean_accessibility analysis; reported beside angle, not read as better or worse
- atac_mean_accessibility export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `residualized`

##### MOTCO preprocessing sensitivity: cohort_30, variant `residualized`

after log2-CPM and before the selection step, each abundance-passing feature replaced by its mean plus the residual of its regression on an intercept and RESIDUALIZATION_COVARIATES (the export's log10 library sizes and the technical table's gene-set fractions); the selection ranks the residual variance; every other export setting equal to the primary. Export setting(s) changed: `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'] (primary: `residualize` = None); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `residualized` row (CV-selected rank 2) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| residualized | cv | 2 | 10.6959 | 0.446 | no | 11.2855 | 8.40232 | 27.4098 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| residualized | cv | 2 | 18.1967 | 0.336 | no | 15.0518 | 11.2436 | 36.3423 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| residualized | cv | 2 | 0.380087 | 0.265 | no | 0.282339 | 0.152251 | 0.570231 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20387 | 5000 | 4252 | 748 |
| rna | residualized | 36601 | 35329 | 20387 | 5000 | 4252 | 748 |
| atac | primary | 218882 | 212816 | 128916 | 5000 | 3447 | 1553 |
| atac | residualized | 218882 | 212816 | 128916 | 5000 | 3447 | 1553 |

###### Residualization

After log2-CPM and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `4aad61a2afe8…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.039 | 0.105 | 0.186 | 0.304 | 0.528 |
| atac | 0.034 | 0.081 | 0.131 | 0.198 | 0.365 |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 5 | 0.6356 | 113.8 | `results/motco_analysis/cohort_30` |
| residualized | cv | 2 | 0.3991 | 130.5 | `motco_analysis/cohort_30` |

###### Findings

- CV-selected rank differs: 5 in the primary analysis, 2 in the residualized analysis; reported beside angle, not read as better or worse
- residualized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `revised_candidate`

##### MOTCO preprocessing sensitivity: cohort_30, variant `revised_candidate`

the four settings together in the primary's order: gene sets removed before the abundance filter, residualization after log2-CPM, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `revised_candidate` row (CV-selected rank 3) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| revised_candidate | cv | 3 | 15.1154 | 0.42 | no | 15.3025 | 12.2357 | 39.8494 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| revised_candidate | cv | 3 | 29.1593 | 0.234 | no | 21.4986 | 11.8591 | 42.7198 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| revised_candidate | cv | 3 | 0.185726 | 0.75 | no | 0.314441 | 0.166225 | 0.61389 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20387 | 5000 | 2282 | 2718 |
| rna | revised_candidate | 36601 | 35329 | 35249 | 20317 | 5000 | 2282 | 2718 |
| atac | primary | 218882 | 212816 | 212816 | 128916 | 5000 | 0 | 5000 |
| atac | revised_candidate | 218882 | 212816 | 212816 | 128916 | 5000 | 0 | 5000 |

###### Gene sets removed

RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps `sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.

| omic | set | matched | removed here | already on an excluded chromosome |
|---|---|---:|---:|---:|
| rna | ribosomal | 81 | 77 | 4 |
| rna | mitochondrial | 13 | 0 | 13 |
| rna | hemoglobin | 3 | 3 | 0 |
| atac | ribosomal | 0 | 0 | 0 |
| atac | mitochondrial | 0 | 0 | 0 |
| atac | hemoglobin | 0 | 0 | 0 |

###### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | primary | variance | -0.721 | -0.137 | 0.774 | 2.160 | 4.968 |
| rna | revised_candidate | trend_ratio | 0.127 | 2.387 | 3.969 | 5.478 | 7.641 |
| atac | primary | variance | -0.321 | 0.046 | 0.371 | 0.778 | 1.531 |
| atac | revised_candidate | mean | 5.436 | 5.613 | 5.876 | 6.274 | 7.156 |

###### Residualization

After log2-CPM and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `e95605f9069f…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.039 | 0.105 | 0.185 | 0.302 | 0.522 |
| atac | 0.034 | 0.081 | 0.131 | 0.198 | 0.365 |

###### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule does not hold: flagged pls_1 × `tmm_rna` (0.38), pls_1 × `tmm_atac` (-0.40), pls_2 × `tmm_atac` (0.37). The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.183 | no |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.165 | no |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.139 | no |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 80 | 0.001 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.219 | no |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.122 | no |
| pls_0 | `rin` | covariates.csv (primary export) | 80 | -0.033 | no |
| pls_0 | `pmi` | covariates.csv (primary export) | 80 | -0.020 | no |
| pls_1 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | 0.377 | yes |
| pls_1 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | -0.400 | yes |
| pls_1 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | -0.087 | no |
| pls_1 | `doublet_score_mean` | donor_technical.csv | 80 | -0.099 | no |
| pls_1 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | 0.191 | no |
| pls_1 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | 0.160 | no |
| pls_1 | `rin` | covariates.csv (primary export) | 80 | 0.130 | no |
| pls_1 | `pmi` | covariates.csv (primary export) | 80 | 0.018 | no |
| pls_2 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.021 | no |
| pls_2 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.366 | yes |
| pls_2 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.117 | no |
| pls_2 | `doublet_score_mean` | donor_technical.csv | 80 | -0.018 | no |
| pls_2 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.120 | no |
| pls_2 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.182 | no |
| pls_2 | `rin` | covariates.csv (primary export) | 80 | -0.064 | no |
| pls_2 | `pmi` | covariates.csv (primary export) | 80 | -0.152 | no |
| all | all | | 80 | | yes → holds = false |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 5 | 0.6356 | 113.8 | `results/motco_analysis/cohort_30` |
| revised_candidate | cv | 3 | 0.4945 | 80.2 | `motco_analysis/cohort_30` |

###### Findings

- adoption check does not hold in cohort_30: pls_1 x tmm_rna (0.38), pls_1 x tmm_atac (-0.40), pls_2 x tmm_atac (0.37); the variants stay sensitivity analyses
- CV-selected rank differs: 5 in the primary analysis, 3 in the revised_candidate analysis; reported beside angle, not read as better or worse
- revised_candidate export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `trend_normalized`

##### MOTCO preprocessing sensitivity: cohort_30, variant `trend_normalized`

after log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features corrected by limma's fast cyclic loess (a lowess curve of its deviation from the per-feature mean over donors against that mean, three iterations, adaptive span) in both omics; every other export setting equal to the primary; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'} (primary: `normalization` = {'rna': 'none', 'atac': 'none'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `trend_normalized` row (CV-selected rank 4) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| trend_normalized | cv | 4 | 6.32863 | 0.655 | no | 11.5817 | 9.10815 | 28.5388 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| trend_normalized | cv | 4 | 24.277 | 0.457 | no | 24.0838 | 11.937 | 45.1637 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| trend_normalized | cv | 4 | 0.0847488 | 0.869 | no | 0.194297 | 0.0937329 | 0.351891 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20387 | 5000 | 4329 | 671 |
| rna | trend_normalized | 36601 | 35329 | 20387 | 5000 | 4329 | 671 |
| atac | primary | 218882 | 212816 | 128916 | 5000 | 3360 | 1640 |
| atac | trend_normalized | 218882 | 212816 | 128916 | 5000 | 3360 | 1640 |

###### Normalization

After log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features was corrected by limma's fast cyclic loess (`normalizeCyclicLoess(method = "fast")`): at each of 3 iterations the reference is the per-feature mean over donors (A), and each donor's deviation from it (M) is fitted by a lowess curve on A (3 robustifying iterations, delta = 0.01 of the range of A) and subtracted. The span is limma's adaptive default from the number of features (`LOESS_SPAN_RULE`: 0.3 + (1 − 0.3) · (50 / n)^(1/3)). Only the block's values are read (no label, covariate or technical variable). The tilt statistics are per donor, with A and M taken against the per-feature mean of the values before or after the normalization: the OLS slope of M on A and the median M of the lowest and highest of 20 equal-count bins of A; the quantiles are over donors and the per-donor values are in `normalization_trend.csv` beside the variant blocks (sha256 `f7192b22689e…`). Limitations: the normalization assumes that most features do not change at any abundance level, so a stage or sex effect shaped like a tilt is removed with the artefact; the curve is least determined at the rare end, where the selection draws most features; the lncRNA-versus-coding offset at matched abundance is not corrected; and no Type I check exists for the normalized blocks.

| omic | primary | variant | features | span | iterations | largest absolute correction |
|---|---|---|---:|---:|---:|---:|
| rna | none | cyclic_loess | 20387 | 0.394400 | 3 | 1.2155 |
| atac | none | cyclic_loess | 128916 | 0.351049 | 3 | 1.1286 |

| omic | statistic | values | q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | slope of M on A | before | -0.0639 | -0.0342 | 0.0012 | 0.0294 | 0.0655 |
| rna | slope of M on A | after | -0.0140 | -0.0066 | -0.0026 | 0.0043 | 0.0188 |
| rna | median M, lowest bin | before | -0.3761 | -0.1549 | 0.0267 | 0.1926 | 0.4588 |
| rna | median M, lowest bin | after | -0.1090 | -0.0214 | 0.0102 | 0.0419 | 0.0651 |
| rna | median M, highest bin | before | -0.2846 | -0.1030 | 0.0048 | 0.1100 | 0.2457 |
| rna | median M, highest bin | after | -0.0145 | -0.0086 | -0.0031 | 0.0022 | 0.0164 |
| rna | largest absolute correction | per donor | 0.1311 | 0.2392 | 0.3396 | 0.5065 | 0.9425 |
| atac | slope of M on A | before | -0.1610 | -0.0483 | 0.0200 | 0.0550 | 0.0826 |
| atac | slope of M on A | after | -0.0247 | -0.0043 | 0.0035 | 0.0072 | 0.0154 |
| atac | median M, lowest bin | before | -0.9403 | -0.8906 | -0.0687 | 0.1940 | 0.6540 |
| atac | median M, lowest bin | after | -0.6540 | -0.4934 | -0.0842 | -0.0499 | 0.1628 |
| atac | median M, highest bin | before | -0.4164 | -0.1617 | 0.0252 | 0.1658 | 0.3549 |
| atac | median M, highest bin | after | -0.0118 | -0.0007 | 0.0085 | 0.0160 | 0.0297 |
| atac | largest absolute correction | per donor | 0.3428 | 0.4949 | 0.6157 | 0.7275 | 0.8875 |

###### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule does not hold: flagged pls_0 × `tmm_rna` (-0.67), pls_0 × `tmm_atac` (0.51), pls_0 × `fraction_mitochondrial_umis_mean` (0.73), pls_0 × `gex_reads_mapped_confidently_to_intronic_regions_mean` (-0.56), pls_0 × `gex_fraction_of_transcriptomic_reads_in_cells_mean` (-0.63), pls_0 × `rin` (-0.33). The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result. Between `trend_normalized` and `revised_candidate_trend`, `ADOPTION_PRECEDENCE` decides from both cohorts' `holds` values: if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read; its outcome is recorded as `adoption_precedence` in the comparison manifests.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.673 | yes |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.505 | yes |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.730 | yes |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 80 | 0.018 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.565 | yes |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.632 | yes |
| pls_0 | `rin` | covariates.csv (primary export) | 80 | -0.335 | yes |
| pls_0 | `pmi` | covariates.csv (primary export) | 80 | 0.032 | no |
| pls_1 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | 0.148 | no |
| pls_1 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.040 | no |
| pls_1 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | -0.106 | no |
| pls_1 | `doublet_score_mean` | donor_technical.csv | 80 | -0.025 | no |
| pls_1 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.102 | no |
| pls_1 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.026 | no |
| pls_1 | `rin` | covariates.csv (primary export) | 80 | -0.039 | no |
| pls_1 | `pmi` | covariates.csv (primary export) | 80 | -0.170 | no |
| pls_2 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.166 | no |
| pls_2 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | -0.073 | no |
| pls_2 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.148 | no |
| pls_2 | `doublet_score_mean` | donor_technical.csv | 80 | -0.141 | no |
| pls_2 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | 0.096 | no |
| pls_2 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | 0.174 | no |
| pls_2 | `rin` | covariates.csv (primary export) | 80 | 0.212 | no |
| pls_2 | `pmi` | covariates.csv (primary export) | 80 | 0.119 | no |
| pls_3 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | 0.238 | no |
| pls_3 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | -0.001 | no |
| pls_3 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | -0.110 | no |
| pls_3 | `doublet_score_mean` | donor_technical.csv | 80 | -0.137 | no |
| pls_3 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | 0.126 | no |
| pls_3 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.079 | no |
| pls_3 | `rin` | covariates.csv (primary export) | 80 | -0.132 | no |
| pls_3 | `pmi` | covariates.csv (primary export) | 80 | -0.116 | no |
| all | all | | 80 | | yes → holds = false |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 5 | 0.6356 | 113.8 | `results/motco_analysis/cohort_30` |
| trend_normalized | cv | 4 | 0.5991 | 152.0 | `motco_analysis/cohort_30` |

###### Findings

- adoption check does not hold in cohort_30: pls_0 x tmm_rna (-0.67), pls_0 x tmm_atac (0.51), pls_0 x fraction_mitochondrial_umis_mean (0.73), pls_0 x gex_reads_mapped_confidently_to_intronic_regions_mean (-0.56), pls_0 x gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63), pls_0 x rin (-0.33); the variants stay sensitivity analyses
- CV-selected rank differs: 5 in the primary analysis, 4 in the trend_normalized analysis; reported beside angle, not read as better or worse
- trend_normalized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `revised_candidate_trend`

##### MOTCO preprocessing sensitivity: cohort_30, variant `revised_candidate_trend`

the revised candidate's four settings plus the cyclic-loess trend normalization of both omics, in the primary's order: gene sets removed before the abundance filter, cyclic loess after log2-CPM, residualization after the normalization, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'}, `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `normalization` = {'rna': 'none', 'atac': 'none'}, `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `revised_candidate_trend` row (CV-selected rank 3) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| revised_candidate_trend | cv | 3 | 14.3836 | 0.197 | no | 9.09659 | 7.00138 | 23.7116 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| revised_candidate_trend | cv | 3 | 16.7741 | 0.472 | no | 17.5615 | 9.94492 | 36.3862 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| revised_candidate_trend | cv | 3 | 0.0792388 | 0.87 | no | 0.200861 | 0.110567 | 0.406137 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20387 | 5000 | 2238 | 2762 |
| rna | revised_candidate_trend | 36601 | 35329 | 35249 | 20317 | 5000 | 2238 | 2762 |
| atac | primary | 218882 | 212816 | 212816 | 128916 | 5000 | 0 | 5000 |
| atac | revised_candidate_trend | 218882 | 212816 | 212816 | 128916 | 5000 | 0 | 5000 |

###### Gene sets removed

RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps `sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.

| omic | set | matched | removed here | already on an excluded chromosome |
|---|---|---:|---:|---:|
| rna | ribosomal | 81 | 77 | 4 |
| rna | mitochondrial | 13 | 0 | 13 |
| rna | hemoglobin | 3 | 3 | 0 |
| atac | ribosomal | 0 | 0 | 0 |
| atac | mitochondrial | 0 | 0 | 0 |
| atac | hemoglobin | 0 | 0 | 0 |

###### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | primary | variance | -0.721 | -0.137 | 0.774 | 2.160 | 4.968 |
| rna | revised_candidate_trend | trend_ratio | -0.031 | 2.348 | 4.000 | 5.532 | 7.735 |
| atac | primary | variance | -0.321 | 0.046 | 0.371 | 0.778 | 1.531 |
| atac | revised_candidate_trend | mean | 5.433 | 5.614 | 5.881 | 6.285 | 7.179 |

###### Normalization

After log2-CPM and before the residualization and the selection step, each donor's log2-CPM over the abundance-passing features was corrected by limma's fast cyclic loess (`normalizeCyclicLoess(method = "fast")`): at each of 3 iterations the reference is the per-feature mean over donors (A), and each donor's deviation from it (M) is fitted by a lowess curve on A (3 robustifying iterations, delta = 0.01 of the range of A) and subtracted. The span is limma's adaptive default from the number of features (`LOESS_SPAN_RULE`: 0.3 + (1 − 0.3) · (50 / n)^(1/3)). Only the block's values are read (no label, covariate or technical variable). The tilt statistics are per donor, with A and M taken against the per-feature mean of the values before or after the normalization: the OLS slope of M on A and the median M of the lowest and highest of 20 equal-count bins of A; the quantiles are over donors and the per-donor values are in `normalization_trend.csv` beside the variant blocks (sha256 `e60d1a63f939…`). Limitations: the normalization assumes that most features do not change at any abundance level, so a stage or sex effect shaped like a tilt is removed with the artefact; the curve is least determined at the rare end, where the selection draws most features; the lncRNA-versus-coding offset at matched abundance is not corrected; and no Type I check exists for the normalized blocks.

| omic | primary | variant | features | span | iterations | largest absolute correction |
|---|---|---|---:|---:|---:|---:|
| rna | none | cyclic_loess | 20317 | 0.394508 | 3 | 1.2130 |
| atac | none | cyclic_loess | 128916 | 0.351049 | 3 | 1.1286 |

| omic | statistic | values | q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | slope of M on A | before | -0.0639 | -0.0342 | 0.0013 | 0.0296 | 0.0660 |
| rna | slope of M on A | after | -0.0139 | -0.0065 | -0.0025 | 0.0043 | 0.0189 |
| rna | median M, lowest bin | before | -0.3774 | -0.1554 | 0.0277 | 0.1933 | 0.4569 |
| rna | median M, lowest bin | after | -0.0996 | -0.0237 | 0.0094 | 0.0385 | 0.0657 |
| rna | median M, highest bin | before | -0.2835 | -0.1029 | 0.0028 | 0.1104 | 0.2436 |
| rna | median M, highest bin | after | -0.0141 | -0.0087 | -0.0033 | 0.0022 | 0.0161 |
| rna | largest absolute correction | per donor | 0.1312 | 0.2384 | 0.3403 | 0.5049 | 0.9387 |
| atac | slope of M on A | before | -0.1610 | -0.0483 | 0.0200 | 0.0550 | 0.0826 |
| atac | slope of M on A | after | -0.0247 | -0.0043 | 0.0035 | 0.0072 | 0.0154 |
| atac | median M, lowest bin | before | -0.9403 | -0.8906 | -0.0687 | 0.1940 | 0.6540 |
| atac | median M, lowest bin | after | -0.6540 | -0.4934 | -0.0842 | -0.0499 | 0.1628 |
| atac | median M, highest bin | before | -0.4164 | -0.1617 | 0.0252 | 0.1658 | 0.3549 |
| atac | median M, highest bin | after | -0.0118 | -0.0007 | 0.0085 | 0.0160 | 0.0297 |
| atac | largest absolute correction | per donor | 0.3428 | 0.4949 | 0.6157 | 0.7275 | 0.8875 |

###### Residualization

After log2-CPM, the trend normalization and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `e95605f9069f…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.037 | 0.096 | 0.169 | 0.267 | 0.478 |
| atac | 0.015 | 0.041 | 0.074 | 0.124 | 0.251 |

###### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule does not hold: flagged pls_1 × `tmm_rna` (-0.32). The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result. Between `trend_normalized` and `revised_candidate_trend`, `ADOPTION_PRECEDENCE` decides from both cohorts' `holds` values: if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read; its outcome is recorded as `adoption_precedence` in the comparison manifests.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.159 | no |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.142 | no |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.159 | no |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 80 | 0.016 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.213 | no |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.121 | no |
| pls_0 | `rin` | covariates.csv (primary export) | 80 | -0.020 | no |
| pls_0 | `pmi` | covariates.csv (primary export) | 80 | -0.002 | no |
| pls_1 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.323 | yes |
| pls_1 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.146 | no |
| pls_1 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.089 | no |
| pls_1 | `doublet_score_mean` | donor_technical.csv | 80 | 0.129 | no |
| pls_1 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.115 | no |
| pls_1 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.039 | no |
| pls_1 | `rin` | covariates.csv (primary export) | 80 | -0.048 | no |
| pls_1 | `pmi` | covariates.csv (primary export) | 80 | 0.084 | no |
| pls_2 | `tmm_rna` | manifest.json (primary export, TMM factor) | 80 | -0.234 | no |
| pls_2 | `tmm_atac` | manifest.json (primary export, TMM factor) | 80 | 0.172 | no |
| pls_2 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.167 | no |
| pls_2 | `doublet_score_mean` | donor_technical.csv | 80 | 0.017 | no |
| pls_2 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.063 | no |
| pls_2 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | 0.046 | no |
| pls_2 | `rin` | covariates.csv (primary export) | 80 | 0.056 | no |
| pls_2 | `pmi` | covariates.csv (primary export) | 80 | 0.035 | no |
| all | all | | 80 | | yes → holds = false |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 5 | 0.6356 | 113.8 | `results/motco_analysis/cohort_30` |
| revised_candidate_trend | cv | 3 | 0.4907 | 151.5 | `motco_analysis/cohort_30` |

###### Findings

- adoption check does not hold in cohort_30: pls_1 x tmm_rna (-0.32); the variants stay sensitivity analyses
- CV-selected rank differs: 5 in the primary analysis, 3 in the revised_candidate_trend analysis; reported beside angle, not read as better or worse
- revised_candidate_trend export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### 6.2 cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74) — sensitivity cohort

#### Variant `keep_chrx`

##### MOTCO preprocessing sensitivity: cohort_50, variant `keep_chrx`

chrX kept; chrY and chrM removed; every other export setting equal to the primary. Export setting(s) changed: `excluded_chromosomes` = ['chrY', 'chrM'] (primary: `excluded_chromosomes` = ['chrX', 'chrY', 'chrM']); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `keep_chrx` row (CV-selected rank 5) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| keep_chrx | cv | 5 | 13.6639 | 0.242 | no | 9.29499 | 6.85111 | 22.9002 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| keep_chrx | cv | 5 | 24.9388 | 0.37 | no | 22.8531 | 10.3079 | 43.3069 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| keep_chrx | cv | 5 | 0.0743536 | 0.804 | no | 0.144712 | 0.0786753 | 0.286751 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports; `kept on chrX` is per export.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept on chrX | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20238 | 5000 | 0 | 4870 | 130 |
| rna | keep_chrx | 36601 | 36477 | 20820 | 5000 | 129 | 4870 | 130 |
| atac | primary | 218882 | 212816 | 124818 | 5000 | 0 | 4885 | 115 |
| atac | keep_chrx | 218882 | 218702 | 126408 | 5000 | 110 | 4885 | 115 |

- `XIST` (chrX, rna): status `sex_chromosome` in the primary export, `kept` in the keep_chrx export.

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 2 | 0.6965 | 102.3 | `results/motco_analysis/cohort_50` |
| keep_chrx | cv | 5 | 0.6917 | 111.7 | `motco_analysis/cohort_50` |

###### Findings

- CV-selected rank differs: 2 in the primary analysis, 5 in the keep_chrx analysis; reported beside angle, not read as better or worse
- keep_chrx export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `exclude_families`

##### MOTCO preprocessing sensitivity: cohort_50, variant `exclude_families`

RNA genes of the ribosomal, mitochondrial and hemoglobin sets (TECHNICAL_GENE_SETS) removed before the abundance filter; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'] (primary: `excluded_gene_sets` = []); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `exclude_families` row (CV-selected rank 5) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| exclude_families | cv | 5 | 13.4939 | 0.247 | no | 9.2985 | 6.85123 | 22.4268 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| exclude_families | cv | 5 | 24.9552 | 0.369 | no | 22.7937 | 10.2553 | 43.4392 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| exclude_families | cv | 5 | 0.0745555 | 0.807 | no | 0.144366 | 0.07841 | 0.286294 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20238 | 5000 | 4977 | 23 |
| rna | exclude_families | 36601 | 35329 | 35249 | 20166 | 5000 | 4977 | 23 |
| atac | primary | 218882 | 212816 | 212816 | 124818 | 5000 | 5000 | 0 |
| atac | exclude_families | 218882 | 212816 | 212816 | 124818 | 5000 | 5000 | 0 |

###### Gene sets removed

RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps `sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.

| omic | set | matched | removed here | already on an excluded chromosome |
|---|---|---:|---:|---:|
| rna | ribosomal | 81 | 77 | 4 |
| rna | mitochondrial | 13 | 0 | 13 |
| rna | hemoglobin | 3 | 3 | 0 |
| atac | ribosomal | 0 | 0 | 0 |
| atac | mitochondrial | 0 | 0 | 0 |
| atac | hemoglobin | 0 | 0 | 0 |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 2 | 0.6965 | 102.3 | `results/motco_analysis/cohort_50` |
| exclude_families | cv | 5 | 0.7028 | 121.6 | `motco_analysis/cohort_50` |

###### Findings

- CV-selected rank differs: 2 in the primary analysis, 5 in the exclude_families analysis; reported beside angle, not read as better or worse
- exclude_families export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `rna_trend_ratio`

##### MOTCO preprocessing sensitivity: cohort_50, variant `rna_trend_ratio`

RNA features ranked at the selection step by the ratio of their log2-CPM variance to a binned-median mean–variance trend instead of by raw variance; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'trend_ratio', 'atac': 'variance'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `rna_trend_ratio` row (CV-selected rank 4) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| rna_trend_ratio | cv | 4 | 7.5918 | 0.593 | no | 11.0692 | 8.26933 | 27.3417 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| rna_trend_ratio | cv | 4 | 25.2598 | 0.392 | no | 23.8522 | 11.7136 | 45.58 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| rna_trend_ratio | cv | 4 | 0.156475 | 0.649 | no | 0.204294 | 0.106282 | 0.404664 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20238 | 5000 | 2572 | 2428 |
| rna | rna_trend_ratio | 36601 | 35329 | 20238 | 5000 | 2572 | 2428 |
| atac | primary | 218882 | 212816 | 124818 | 5000 | 5000 | 0 |
| atac | rna_trend_ratio | 218882 | 212816 | 124818 | 5000 | 5000 | 0 |

###### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | primary | variance | -0.662 | -0.096 | 0.818 | 2.229 | 4.957 |
| rna | rna_trend_ratio | trend_ratio | 0.100 | 2.105 | 3.813 | 5.342 | 7.526 |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 2 | 0.6965 | 102.3 | `results/motco_analysis/cohort_50` |
| rna_trend_ratio | cv | 4 | 0.6918 | 136.4 | `motco_analysis/cohort_50` |

###### Findings

- CV-selected rank differs: 2 in the primary analysis, 4 in the rna_trend_ratio analysis; reported beside angle, not read as better or worse
- rna_trend_ratio export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `atac_mean_accessibility`

##### MOTCO preprocessing sensitivity: cohort_50, variant `atac_mean_accessibility`

ATAC peaks ranked at the selection step by mean log2-CPM (total accessibility) instead of by variance; the RNA block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'variance', 'atac': 'mean'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `atac_mean_accessibility` row (CV-selected rank 7) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| atac_mean_accessibility | cv | 7 | 9.05688 | 0.491 | no | 11.6188 | 9.9148 | 30.3005 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| atac_mean_accessibility | cv | 7 | 23.1367 | 0.696 | no | 32.4082 | 15.1021 | 62.5449 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| atac_mean_accessibility | cv | 7 | 0.0993176 | 0.853 | no | 0.221136 | 0.113437 | 0.420583 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20238 | 5000 | 5000 | 0 |
| rna | atac_mean_accessibility | 36601 | 35329 | 20238 | 5000 | 5000 | 0 |
| atac | primary | 218882 | 212816 | 124818 | 5000 | 0 | 5000 |
| atac | atac_mean_accessibility | 218882 | 212816 | 124818 | 5000 | 0 | 5000 |

###### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| atac | primary | variance | -0.432 | -0.057 | 0.271 | 0.698 | 1.479 |
| atac | atac_mean_accessibility | mean | 5.443 | 5.621 | 5.883 | 6.281 | 7.163 |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 2 | 0.6965 | 102.3 | `results/motco_analysis/cohort_50` |
| atac_mean_accessibility | cv | 7 | 0.7271 | 126.3 | `motco_analysis/cohort_50` |

###### Findings

- CV-selected rank differs: 2 in the primary analysis, 7 in the atac_mean_accessibility analysis; reported beside angle, not read as better or worse
- atac_mean_accessibility export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `residualized`

##### MOTCO preprocessing sensitivity: cohort_50, variant `residualized`

after log2-CPM and before the selection step, each abundance-passing feature replaced by its mean plus the residual of its regression on an intercept and RESIDUALIZATION_COVARIATES (the export's log10 library sizes and the technical table's gene-set fractions); the selection ranks the residual variance; every other export setting equal to the primary. Export setting(s) changed: `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'] (primary: `residualize` = None); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `residualized` row (CV-selected rank 2) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| residualized | cv | 2 | 8.94403 | 0.486 | no | 10.2432 | 8.11982 | 26.0801 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| residualized | cv | 2 | 30.0936 | 0.134 | no | 15.5567 | 12.1783 | 38.6128 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| residualized | cv | 2 | 0.317151 | 0.376 | no | 0.281117 | 0.150903 | 0.54998 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20238 | 5000 | 4218 | 782 |
| rna | residualized | 36601 | 35329 | 20238 | 5000 | 4218 | 782 |
| atac | primary | 218882 | 212816 | 124818 | 5000 | 3423 | 1577 |
| atac | residualized | 218882 | 212816 | 124818 | 5000 | 3423 | 1577 |

###### Residualization

After log2-CPM and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `4773eec7a84e…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.043 | 0.115 | 0.201 | 0.324 | 0.550 |
| atac | 0.032 | 0.080 | 0.131 | 0.202 | 0.378 |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 2 | 0.6965 | 102.3 | `results/motco_analysis/cohort_50` |
| residualized | cv | 2 | 0.4942 | 112.4 | `motco_analysis/cohort_50` |

###### Findings

- residualized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `revised_candidate`

##### MOTCO preprocessing sensitivity: cohort_50, variant `revised_candidate`

the four settings together in the primary's order: gene sets removed before the abundance filter, residualization after log2-CPM, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `revised_candidate` row (CV-selected rank 2) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| revised_candidate | cv | 2 | 6.33891 | 0.603 | no | 10.5045 | 8.2407 | 27.5045 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| revised_candidate | cv | 2 | 2.15631 | 0.906 | no | 13.3699 | 10.5336 | 32.7749 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| revised_candidate | cv | 2 | 0.230511 | 0.602 | no | 0.293981 | 0.165487 | 0.60537 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20238 | 5000 | 2300 | 2700 |
| rna | revised_candidate | 36601 | 35329 | 35249 | 20166 | 5000 | 2300 | 2700 |
| atac | primary | 218882 | 212816 | 212816 | 124818 | 5000 | 0 | 5000 |
| atac | revised_candidate | 218882 | 212816 | 212816 | 124818 | 5000 | 0 | 5000 |

###### Gene sets removed

RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps `sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.

| omic | set | matched | removed here | already on an excluded chromosome |
|---|---|---:|---:|---:|
| rna | ribosomal | 81 | 77 | 4 |
| rna | mitochondrial | 13 | 0 | 13 |
| rna | hemoglobin | 3 | 3 | 0 |
| atac | ribosomal | 0 | 0 | 0 |
| atac | mitochondrial | 0 | 0 | 0 |
| atac | hemoglobin | 0 | 0 | 0 |

###### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | primary | variance | -0.662 | -0.096 | 0.818 | 2.229 | 4.957 |
| rna | revised_candidate | trend_ratio | 0.148 | 2.413 | 3.980 | 5.501 | 7.651 |
| atac | primary | variance | -0.432 | -0.057 | 0.271 | 0.698 | 1.479 |
| atac | revised_candidate | mean | 5.443 | 5.621 | 5.883 | 6.281 | 7.163 |

###### Residualization

After log2-CPM and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `443e90258051…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.043 | 0.114 | 0.199 | 0.321 | 0.544 |
| atac | 0.032 | 0.080 | 0.131 | 0.202 | 0.378 |

###### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule does not hold: flagged pls_0 × `tmm_atac` (0.37), pls_0 × `gex_reads_mapped_confidently_to_intronic_regions_mean` (-0.33), pls_1 × `tmm_atac` (-0.37). The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | -0.287 | no |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | 0.375 | yes |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | 0.204 | no |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 74 | 0.034 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | -0.325 | yes |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | -0.197 | no |
| pls_0 | `rin` | covariates.csv (primary export) | 74 | -0.113 | no |
| pls_0 | `pmi` | covariates.csv (primary export) | 74 | -0.015 | no |
| pls_1 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | 0.209 | no |
| pls_1 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | -0.371 | yes |
| pls_1 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | -0.030 | no |
| pls_1 | `doublet_score_mean` | donor_technical.csv | 74 | -0.115 | no |
| pls_1 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | 0.134 | no |
| pls_1 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | 0.161 | no |
| pls_1 | `rin` | covariates.csv (primary export) | 74 | 0.202 | no |
| pls_1 | `pmi` | covariates.csv (primary export) | 74 | 0.065 | no |
| all | all | | 74 | | yes → holds = false |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 2 | 0.6965 | 102.3 | `results/motco_analysis/cohort_50` |
| revised_candidate | cv | 2 | 0.5514 | 116.5 | `motco_analysis/cohort_50` |

###### Findings

- adoption check does not hold in cohort_50: pls_0 x tmm_atac (0.37), pls_0 x gex_reads_mapped_confidently_to_intronic_regions_mean (-0.33), pls_1 x tmm_atac (-0.37); the variants stay sensitivity analyses
- revised_candidate export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `trend_normalized`

##### MOTCO preprocessing sensitivity: cohort_50, variant `trend_normalized`

after log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features corrected by limma's fast cyclic loess (a lowess curve of its deviation from the per-feature mean over donors against that mean, three iterations, adaptive span) in both omics; every other export setting equal to the primary; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'} (primary: `normalization` = {'rna': 'none', 'atac': 'none'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `trend_normalized` row (CV-selected rank 4) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| trend_normalized | cv | 4 | 9.24514 | 0.554 | no | 12.5908 | 9.63945 | 30.7129 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| trend_normalized | cv | 4 | 25.992 | 0.356 | no | 23.0856 | 11.5907 | 44.6216 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| trend_normalized | cv | 4 | 0.116966 | 0.737 | no | 0.18593 | 0.0953054 | 0.361841 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20238 | 5000 | 4342 | 658 |
| rna | trend_normalized | 36601 | 35329 | 20238 | 5000 | 4342 | 658 |
| atac | primary | 218882 | 212816 | 124818 | 5000 | 3257 | 1743 |
| atac | trend_normalized | 218882 | 212816 | 124818 | 5000 | 3257 | 1743 |

###### Normalization

After log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features was corrected by limma's fast cyclic loess (`normalizeCyclicLoess(method = "fast")`): at each of 3 iterations the reference is the per-feature mean over donors (A), and each donor's deviation from it (M) is fitted by a lowess curve on A (3 robustifying iterations, delta = 0.01 of the range of A) and subtracted. The span is limma's adaptive default from the number of features (`LOESS_SPAN_RULE`: 0.3 + (1 − 0.3) · (50 / n)^(1/3)). Only the block's values are read (no label, covariate or technical variable). The tilt statistics are per donor, with A and M taken against the per-feature mean of the values before or after the normalization: the OLS slope of M on A and the median M of the lowest and highest of 20 equal-count bins of A; the quantiles are over donors and the per-donor values are in `normalization_trend.csv` beside the variant blocks (sha256 `8a60e38aab48…`). Limitations: the normalization assumes that most features do not change at any abundance level, so a stage or sex effect shaped like a tilt is removed with the artefact; the curve is least determined at the rare end, where the selection draws most features; the lncRNA-versus-coding offset at matched abundance is not corrected; and no Type I check exists for the normalized blocks.

| omic | primary | variant | features | span | iterations | largest absolute correction |
|---|---|---|---:|---:|---:|---:|
| rna | none | cyclic_loess | 20238 | 0.394631 | 3 | 1.2402 |
| atac | none | cyclic_loess | 124818 | 0.351601 | 3 | 1.1331 |

| omic | statistic | values | q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | slope of M on A | before | -0.0626 | -0.0337 | 0.0001 | 0.0299 | 0.0673 |
| rna | slope of M on A | after | -0.0138 | -0.0061 | -0.0029 | 0.0049 | 0.0176 |
| rna | median M, lowest bin | before | -0.3766 | -0.1769 | 0.0506 | 0.2020 | 0.4238 |
| rna | median M, lowest bin | after | -0.1472 | -0.0196 | 0.0179 | 0.0407 | 0.0767 |
| rna | median M, highest bin | before | -0.2576 | -0.1153 | 0.0109 | 0.1169 | 0.2688 |
| rna | median M, highest bin | after | -0.0146 | -0.0080 | -0.0035 | 0.0025 | 0.0143 |
| rna | largest absolute correction | per donor | 0.1478 | 0.2433 | 0.3612 | 0.5583 | 0.9439 |
| atac | slope of M on A | before | -0.1629 | -0.0510 | 0.0186 | 0.0605 | 0.0913 |
| atac | slope of M on A | after | -0.0235 | -0.0040 | 0.0019 | 0.0069 | 0.0169 |
| atac | median M, lowest bin | before | -1.0531 | -0.8624 | -0.0226 | 0.2940 | 0.6626 |
| atac | median M, lowest bin | after | -0.6301 | -0.3593 | -0.0674 | -0.0205 | 0.2357 |
| atac | median M, highest bin | before | -0.4462 | -0.1349 | 0.0489 | 0.1695 | 0.3353 |
| atac | median M, highest bin | after | -0.0134 | -0.0016 | 0.0081 | 0.0161 | 0.0298 |
| atac | largest absolute correction | per donor | 0.3368 | 0.4902 | 0.6114 | 0.8138 | 0.9805 |

###### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule does not hold: flagged pls_0 × `tmm_rna` (-0.66), pls_0 × `tmm_atac` (0.58), pls_0 × `fraction_mitochondrial_umis_mean` (0.78), pls_0 × `gex_reads_mapped_confidently_to_intronic_regions_mean` (-0.60), pls_0 × `gex_fraction_of_transcriptomic_reads_in_cells_mean` (-0.63), pls_0 × `rin` (-0.34). The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result. Between `trend_normalized` and `revised_candidate_trend`, `ADOPTION_PRECEDENCE` decides from both cohorts' `holds` values: if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read; its outcome is recorded as `adoption_precedence` in the comparison manifests.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | -0.661 | yes |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | 0.578 | yes |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | 0.783 | yes |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 74 | 0.011 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | -0.596 | yes |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | -0.634 | yes |
| pls_0 | `rin` | covariates.csv (primary export) | 74 | -0.341 | yes |
| pls_0 | `pmi` | covariates.csv (primary export) | 74 | 0.072 | no |
| pls_1 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | 0.123 | no |
| pls_1 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | 0.033 | no |
| pls_1 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | -0.082 | no |
| pls_1 | `doublet_score_mean` | donor_technical.csv | 74 | -0.023 | no |
| pls_1 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | -0.139 | no |
| pls_1 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | 0.011 | no |
| pls_1 | `rin` | covariates.csv (primary export) | 74 | -0.069 | no |
| pls_1 | `pmi` | covariates.csv (primary export) | 74 | -0.055 | no |
| pls_2 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | -0.145 | no |
| pls_2 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | -0.171 | no |
| pls_2 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | 0.047 | no |
| pls_2 | `doublet_score_mean` | donor_technical.csv | 74 | -0.198 | no |
| pls_2 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | 0.164 | no |
| pls_2 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | 0.236 | no |
| pls_2 | `rin` | covariates.csv (primary export) | 74 | 0.217 | no |
| pls_2 | `pmi` | covariates.csv (primary export) | 74 | 0.077 | no |
| pls_3 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | -0.193 | no |
| pls_3 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | -0.072 | no |
| pls_3 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | 0.058 | no |
| pls_3 | `doublet_score_mean` | donor_technical.csv | 74 | 0.079 | no |
| pls_3 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | 0.001 | no |
| pls_3 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | 0.170 | no |
| pls_3 | `rin` | covariates.csv (primary export) | 74 | 0.165 | no |
| pls_3 | `pmi` | covariates.csv (primary export) | 74 | 0.095 | no |
| all | all | | 74 | | yes → holds = false |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 2 | 0.6965 | 102.3 | `results/motco_analysis/cohort_50` |
| trend_normalized | cv | 4 | 0.6722 | 159.2 | `motco_analysis/cohort_50` |

###### Findings

- adoption check does not hold in cohort_50: pls_0 x tmm_rna (-0.66), pls_0 x tmm_atac (0.58), pls_0 x fraction_mitochondrial_umis_mean (0.78), pls_0 x gex_reads_mapped_confidently_to_intronic_regions_mean (-0.60), pls_0 x gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63), pls_0 x rin (-0.34); the variants stay sensitivity analyses
- CV-selected rank differs: 2 in the primary analysis, 4 in the trend_normalized analysis; reported beside angle, not read as better or worse
- trend_normalized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

#### Variant `revised_candidate_trend`

##### MOTCO preprocessing sensitivity: cohort_50, variant `revised_candidate_trend`

the revised candidate's four settings plus the cyclic-loess trend normalization of both omics, in the primary's order: gene sets removed before the abundance filter, cyclic loess after log2-CPM, residualization after the normalization, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'}, `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `normalization` = {'rna': 'none', 'atac': 'none'}, `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `revised_candidate_trend` row (CV-selected rank 2) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| revised_candidate_trend | cv | 2 | 5.41504 | 0.618 | no | 8.81492 | 6.83489 | 22.4536 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| revised_candidate_trend | cv | 2 | 19.9809 | 0.249 | no | 13.7096 | 10.5005 | 33.6364 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| revised_candidate_trend | cv | 2 | 0.259816 | 0.506 | no | 0.27208 | 0.138222 | 0.52147 |

###### Features

Counts after each export step per omic; `kept only here` is per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after gene-set removal | after abundance | kept (after selection) | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 35329 | 20238 | 5000 | 2272 | 2728 |
| rna | revised_candidate_trend | 36601 | 35329 | 35249 | 20166 | 5000 | 2272 | 2728 |
| atac | primary | 218882 | 212816 | 212816 | 124818 | 5000 | 0 | 5000 |
| atac | revised_candidate_trend | 218882 | 212816 | 212816 | 124818 | 5000 | 0 | 5000 |

###### Gene sets removed

RNA genes matched by the excluded sets of `TECHNICAL_GENE_SETS` are removed at step 2 after the chromosome removal (status `excluded_gene_set` in `features_rna.csv`); a matched gene on an excluded chromosome keeps `sex_chromosome` and is counted under `already on an excluded chromosome`. ATAC names match no set.

| omic | set | matched | removed here | already on an excluded chromosome |
|---|---|---:|---:|---:|
| rna | ribosomal | 81 | 77 | 4 |
| rna | mitochondrial | 13 | 0 | 13 |
| rna | hemoglobin | 3 | 3 | 0 |
| atac | ribosomal | 0 | 0 | 0 |
| atac | mitochondrial | 0 | 0 | 0 |
| atac | hemoglobin | 0 | 0 | 0 |

###### Selection rule

`variance`: population variance (ddof 0) of log2-CPM, the primary's rule. `trend_ratio`: variance over a mean–variance trend, the median variance of 20 equal-count bins of mean log2-CPM placed at each bin's median mean and linearly interpolated (flat beyond the end bins; the bins are in the variant export's manifest). `mean`: mean log2-CPM. Ties by name under every rule; the quantiles are of the kept features' mean log2-CPM.

| omic | export | rule | kept mean q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | primary | variance | -0.662 | -0.096 | 0.818 | 2.229 | 4.957 |
| rna | revised_candidate_trend | trend_ratio | -0.027 | 2.354 | 4.000 | 5.524 | 7.718 |
| atac | primary | variance | -0.432 | -0.057 | 0.271 | 0.698 | 1.479 |
| atac | revised_candidate_trend | mean | 5.440 | 5.620 | 5.886 | 6.290 | 7.183 |

###### Normalization

After log2-CPM and before the residualization and the selection step, each donor's log2-CPM over the abundance-passing features was corrected by limma's fast cyclic loess (`normalizeCyclicLoess(method = "fast")`): at each of 3 iterations the reference is the per-feature mean over donors (A), and each donor's deviation from it (M) is fitted by a lowess curve on A (3 robustifying iterations, delta = 0.01 of the range of A) and subtracted. The span is limma's adaptive default from the number of features (`LOESS_SPAN_RULE`: 0.3 + (1 − 0.3) · (50 / n)^(1/3)). Only the block's values are read (no label, covariate or technical variable). The tilt statistics are per donor, with A and M taken against the per-feature mean of the values before or after the normalization: the OLS slope of M on A and the median M of the lowest and highest of 20 equal-count bins of A; the quantiles are over donors and the per-donor values are in `normalization_trend.csv` beside the variant blocks (sha256 `a88681aba998…`). Limitations: the normalization assumes that most features do not change at any abundance level, so a stage or sex effect shaped like a tilt is removed with the artefact; the curve is least determined at the rare end, where the selection draws most features; the lncRNA-versus-coding offset at matched abundance is not corrected; and no Type I check exists for the normalized blocks.

| omic | primary | variant | features | span | iterations | largest absolute correction |
|---|---|---|---:|---:|---:|---:|
| rna | none | cyclic_loess | 20166 | 0.394743 | 3 | 1.2390 |
| atac | none | cyclic_loess | 124818 | 0.351601 | 3 | 1.1331 |

| omic | statistic | values | q05 | q25 | q50 | q75 | q95 |
|---|---|---|---:|---:|---:|---:|---:|
| rna | slope of M on A | before | -0.0627 | -0.0337 | 0.0009 | 0.0301 | 0.0674 |
| rna | slope of M on A | after | -0.0135 | -0.0061 | -0.0028 | 0.0049 | 0.0175 |
| rna | median M, lowest bin | before | -0.3886 | -0.1829 | 0.0441 | 0.2033 | 0.4315 |
| rna | median M, lowest bin | after | -0.1469 | -0.0235 | 0.0148 | 0.0415 | 0.0769 |
| rna | median M, highest bin | before | -0.2553 | -0.1044 | 0.0040 | 0.1160 | 0.2627 |
| rna | median M, highest bin | after | -0.0149 | -0.0085 | -0.0033 | 0.0027 | 0.0142 |
| rna | largest absolute correction | per donor | 0.1561 | 0.2340 | 0.3503 | 0.5367 | 0.9357 |
| atac | slope of M on A | before | -0.1629 | -0.0510 | 0.0186 | 0.0605 | 0.0913 |
| atac | slope of M on A | after | -0.0235 | -0.0040 | 0.0019 | 0.0069 | 0.0169 |
| atac | median M, lowest bin | before | -1.0531 | -0.8624 | -0.0226 | 0.2940 | 0.6626 |
| atac | median M, lowest bin | after | -0.6301 | -0.3593 | -0.0674 | -0.0205 | 0.2357 |
| atac | median M, highest bin | before | -0.4462 | -0.1349 | 0.0489 | 0.1695 | 0.3353 |
| atac | median M, highest bin | after | -0.0134 | -0.0016 | 0.0081 | 0.0161 | 0.0298 |
| atac | largest absolute correction | per donor | 0.3368 | 0.4902 | 0.6114 | 0.8138 | 0.9805 |

###### Residualization

After log2-CPM, the trend normalization and before the selection step, each abundance-passing feature was replaced by its mean plus the residual of its OLS fit on an intercept and 4 covariates: `log10_lib_rna`, `log10_lib_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (`RESIDUALIZATION_COVARIATES`: the export's own log10 library sizes after step 2 and the technical table's gene-set pseudobulk fractions; excluded as cell state under the technical-axis report's rule: `hemoglobin_pseudobulk` (reports/motco-astrocyte-technical-axis-2026-10-19.md, section 2, 'Classification under the spec's rule': median within-donor supertype share 0.00394 above the between-donor share 0.00285 (both negligible; median pseudobulk fraction 6.7e-07 of UMIs); not a residualization candidate)). The design (`residualization.csv` beside the variant blocks, sha256 `443e90258051…`) holds no sex, stage or analysis covariate. Limitations: the fit spends 5 degrees of freedom per feature (intercept plus covariates) before the permutation null, which does not see them; it removes whatever part of the stage and sex signal the covariates carry; the analysis covariates were left unchanged, so the library-size proxies are adjusted twice (log10 library sizes here, log10 nuclei counts in the RRPP model); and no Type I check exists for the residualized blocks.

| omic | R² q05 | q25 | q50 | q75 | q95 |
|---|---:|---:|---:|---:|---:|
| rna | 0.039 | 0.104 | 0.182 | 0.286 | 0.501 |
| atac | 0.016 | 0.044 | 0.079 | 0.135 | 0.277 |

###### Adoption check

`ADOPTION_RULE`: holds when no latent component of the candidate's CV-selected analysis reaches |Spearman| >= DIAGNOSTIC_CORRELATION_THRESHOLD against any variable of ADOPTION_CHECK_VARIABLES, in both cohorts; p-values, rejection statuses and CV AUROC are not read. If it holds, a later change may propose the candidate as a revised primary with the current primary kept as the reference; if not, the variants stay sensitivity analyses. Threshold 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`). In this cohort the rule holds. The rule needs both cohorts; whether the candidate becomes a revised primary is a later change with the current primary kept as the reference, and the candidate's statistics above are sensitivity, not the result. Between `trend_normalized` and `revised_candidate_trend`, `ADOPTION_PRECEDENCE` decides from both cohorts' `holds` values: if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read; its outcome is recorded as `adoption_precedence` in the comparison manifests.

| component | variable | source | n | Spearman | flagged |
|---|---|---|---:|---:|:---:|
| pls_0 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | -0.209 | no |
| pls_0 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | 0.233 | no |
| pls_0 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | 0.206 | no |
| pls_0 | `doublet_score_mean` | donor_technical.csv | 74 | 0.026 | no |
| pls_0 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | -0.275 | no |
| pls_0 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | -0.157 | no |
| pls_0 | `rin` | covariates.csv (primary export) | 74 | -0.037 | no |
| pls_0 | `pmi` | covariates.csv (primary export) | 74 | 0.035 | no |
| pls_1 | `tmm_rna` | manifest.json (primary export, TMM factor) | 74 | 0.220 | no |
| pls_1 | `tmm_atac` | manifest.json (primary export, TMM factor) | 74 | -0.158 | no |
| pls_1 | `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | -0.070 | no |
| pls_1 | `doublet_score_mean` | donor_technical.csv | 74 | -0.132 | no |
| pls_1 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | 0.119 | no |
| pls_1 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | 0.099 | no |
| pls_1 | `rin` | covariates.csv (primary export) | 74 | 0.133 | no |
| pls_1 | `pmi` | covariates.csv (primary export) | 74 | -0.063 | no |
| all | all | | 74 | | no → holds = true |

###### Runs

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 2 | 0.6965 | 102.3 | `results/motco_analysis/cohort_50` |
| revised_candidate_trend | cv | 2 | 0.5588 | 156.0 | `motco_analysis/cohort_50` |

###### Findings

- adoption check holds in cohort_50: no candidate component reaches |Spearman| >= 0.3 against any of 8 check variables (the rule needs both cohorts; a revised primary is a later change)
- revised_candidate_trend export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### 6.3 Adoption precedence

From the comparison manifests of `trend_normalized` and `revised_candidate_trend` (identical in all four): holds `trend_normalized` cohort_30 no, cohort_50 no; `revised_candidate_trend` cohort_30 no, cohort_50 yes; precedence: none — "neither trend_normalized nor revised_candidate_trend meets ADOPTION_RULE in both cohorts; no revised-primary change follows and the variants stay sensitivity analyses". Equal to the superseded precedence: True.

### 6.4 Superseded → corrected (variants)

| cell type | cohort | analysis | rank | CV AUROC | κ (n, reduced rank) | statistic | observed | p superseded → corrected | change | rejected superseded → corrected |
|---|---|---|---:|---:|---|---|---:|---|---|---|
| astrocyte | cohort_30 | keep_chrx | 4 | 0.6341 | 1.0690 (80, 10) | delta | 19.0474 | 0.055 → 0.083 | up | no → no |
| astrocyte | cohort_30 | keep_chrx |  |  |  | angle | 22.3677 | 0.353 → 0.396 | up | no → no |
| astrocyte | cohort_30 | keep_chrx |  |  |  | shape | 0.0237388 | 0.988 → 0.991 | up | no → no |
| astrocyte | cohort_50 | keep_chrx | 5 | 0.6917 | 1.0753 (74, 10) | delta | 13.6639 | 0.207 → 0.242 | up | no → no |
| astrocyte | cohort_50 | keep_chrx |  |  |  | angle | 24.9388 | 0.297 → 0.37 | up | no → no |
| astrocyte | cohort_50 | keep_chrx |  |  |  | shape | 0.0743536 | 0.78 → 0.804 | up | no → no |
| astrocyte | cohort_30 | exclude_families | 5 | 0.6369 | 1.0690 (80, 10) | delta | 6.05273 | 0.455 → 0.493 | up | no → no |
| astrocyte | cohort_30 | exclude_families |  |  |  | angle | 21.9731 | 0.396 → 0.45 | up | no → no |
| astrocyte | cohort_30 | exclude_families |  |  |  | shape | 0.156521 | 0.282 → 0.328 | up | no → no |
| astrocyte | cohort_50 | exclude_families | 5 | 0.7028 | 1.0753 (74, 10) | delta | 13.4939 | 0.211 → 0.247 | up | no → no |
| astrocyte | cohort_50 | exclude_families |  |  |  | angle | 24.9552 | 0.295 → 0.369 | up | no → no |
| astrocyte | cohort_50 | exclude_families |  |  |  | shape | 0.0745555 | 0.78 → 0.807 | up | no → no |
| astrocyte | cohort_30 | rna_trend_ratio | 6 | 0.6351 | 1.0690 (80, 10) | delta | 5.89665 | 0.446 → 0.48 | up | no → no |
| astrocyte | cohort_30 | rna_trend_ratio |  |  |  | angle | 22.2081 | 0.479 → 0.53 | up | no → no |
| astrocyte | cohort_30 | rna_trend_ratio |  |  |  | shape | 0.144967 | 0.261 → 0.316 | up | no → no |
| astrocyte | cohort_50 | rna_trend_ratio | 4 | 0.6918 | 1.0753 (74, 10) | delta | 7.5918 | 0.562 → 0.593 | up | no → no |
| astrocyte | cohort_50 | rna_trend_ratio |  |  |  | angle | 25.2598 | 0.346 → 0.392 | up | no → no |
| astrocyte | cohort_50 | rna_trend_ratio |  |  |  | shape | 0.156475 | 0.606 → 0.649 | up | no → no |
| astrocyte | cohort_30 | atac_mean_accessibility | 7 | 0.6706 | 1.0690 (80, 10) | delta | 0.18785 | 0.99 → 0.99 | equal | no → no |
| astrocyte | cohort_30 | atac_mean_accessibility |  |  |  | angle | 18.6469 | 0.772 → 0.814 | up | no → no |
| astrocyte | cohort_30 | atac_mean_accessibility |  |  |  | shape | 0.102546 | 0.794 → 0.819 | up | no → no |
| astrocyte | cohort_50 | atac_mean_accessibility | 7 | 0.7271 | 1.0753 (74, 10) | delta | 9.05688 | 0.455 → 0.491 | up | no → no |
| astrocyte | cohort_50 | atac_mean_accessibility |  |  |  | angle | 23.1367 | 0.641 → 0.696 | up | no → no |
| astrocyte | cohort_50 | atac_mean_accessibility |  |  |  | shape | 0.0993176 | 0.84 → 0.853 | up | no → no |
| astrocyte | cohort_30 | residualized | 2 | 0.3991 | 1.0690 (80, 10) | delta | 10.6959 | 0.419 → 0.446 | up | no → no |
| astrocyte | cohort_30 | residualized |  |  |  | angle | 18.1967 | 0.308 → 0.336 | up | no → no |
| astrocyte | cohort_30 | residualized |  |  |  | shape | 0.380087 | 0.236 → 0.265 | up | no → no |
| astrocyte | cohort_50 | residualized | 2 | 0.4942 | 1.0753 (74, 10) | delta | 8.94403 | 0.444 → 0.486 | up | no → no |
| astrocyte | cohort_50 | residualized |  |  |  | angle | 30.0936 | 0.104 → 0.134 | up | no → no |
| astrocyte | cohort_50 | residualized |  |  |  | shape | 0.317151 | 0.338 → 0.376 | up | no → no |
| astrocyte | cohort_30 | revised_candidate | 3 | 0.4945 | 1.0690 (80, 10) | delta | 15.1154 | 0.391 → 0.42 | up | no → no |
| astrocyte | cohort_30 | revised_candidate |  |  |  | angle | 29.1593 | 0.19 → 0.234 | up | no → no |
| astrocyte | cohort_30 | revised_candidate |  |  |  | shape | 0.185726 | 0.735 → 0.75 | up | no → no |
| astrocyte | cohort_50 | revised_candidate | 2 | 0.5514 | 1.0753 (74, 10) | delta | 6.33891 | 0.585 → 0.603 | up | no → no |
| astrocyte | cohort_50 | revised_candidate |  |  |  | angle | 2.15631 | 0.902 → 0.906 | up | no → no |
| astrocyte | cohort_50 | revised_candidate |  |  |  | shape | 0.230511 | 0.564 → 0.602 | up | no → no |
| astrocyte | cohort_30 | trend_normalized | 4 | 0.5991 | 1.0690 (80, 10) | delta | 6.32863 | 0.63 → 0.655 | up | no → no |
| astrocyte | cohort_30 | trend_normalized |  |  |  | angle | 24.277 | 0.417 → 0.457 | up | no → no |
| astrocyte | cohort_30 | trend_normalized |  |  |  | shape | 0.0847488 | 0.856 → 0.869 | up | no → no |
| astrocyte | cohort_50 | trend_normalized | 4 | 0.6722 | 1.0753 (74, 10) | delta | 9.24514 | 0.529 → 0.554 | up | no → no |
| astrocyte | cohort_50 | trend_normalized |  |  |  | angle | 25.992 | 0.296 → 0.356 | up | no → no |
| astrocyte | cohort_50 | trend_normalized |  |  |  | shape | 0.116966 | 0.709 → 0.737 | up | no → no |
| astrocyte | cohort_30 | revised_candidate_trend | 3 | 0.4907 | 1.0690 (80, 10) | delta | 14.3836 | 0.166 → 0.197 | up | no → no |
| astrocyte | cohort_30 | revised_candidate_trend |  |  |  | angle | 16.7741 | 0.435 → 0.472 | up | no → no |
| astrocyte | cohort_30 | revised_candidate_trend |  |  |  | shape | 0.0792388 | 0.852 → 0.87 | up | no → no |
| astrocyte | cohort_50 | revised_candidate_trend | 2 | 0.5588 | 1.0753 (74, 10) | delta | 5.41504 | 0.585 → 0.618 | up | no → no |
| astrocyte | cohort_50 | revised_candidate_trend |  |  |  | angle | 19.9809 | 0.22 → 0.249 | up | no → no |
| astrocyte | cohort_50 | revised_candidate_trend |  |  |  | shape | 0.259816 | 0.459 → 0.506 | up | no → no |

Superseded against corrected, generated by the comparison script over 16 analyses: of the 64 analysis CSVs (four per analysis), 64 equal at rounding level (largest |difference| 1.3e-10); of the 16 records outside the null keys, 16 equal at rounding level (largest relative difference 5.3e-11, apart from the numerically-zero third eigenvalues of the configuration spectra, ≤ 2.5e-32); rank equal in 16 of 16, CV AUROC equal in 16 of 16; κ = √(n / (n − r)) to 1e-12 in 16 of 16; p-values up in 47, equal in 1, down in 0 of 48; rejection statuses changed: 0.

## 7. Latent-adjustment sensitivity

Both arms were rerun for both cell types on the unchanged primary exports with the eight covariates (reduced rank 12, κ 1.0847 / 1.0925 for astrocyte `cohort_30` / `cohort_50`, 1.0897 / 1.1114 for microglia) and the correction. No arm rejects any statistic, and no rejection status differs from the corrected primary or from the superseded arms. The `extended_covariates` arms keep their primaries' rank and CV AUROC in all four cohorts, as MOTCO's contract states. The mediator check reads no p-value and is byte-identical to the superseded one; its flags stand (astrocyte: `ribosomal_pseudobulk` −0.3993 and `mitochondrial_pseudobulk` 0.3118 in `cohort_30`, −0.3518 and 0.356 in `cohort_50`, all donors; microglia: `mitochondrial_pseudobulk` 0.4978 in `cohort_30` and 0.4397 in `cohort_50`, all donors; per-sex flags in the pasted tables), so the `latent_adjusted` basis may have removed stage signal carried by those fractions.

**The `latent_adjusted` arm's `angle` and `shape` are reported but not read.** MOTCO's diagnosis classes the `stage` latent adjustment at this configuration as `residualized_basis`: liberal on anchor `angle` (0.084–0.090) and `shape` (0.078–0.092) against 0.0695, and on `angle` still at a forced rank of 2 (0.084 / 0.082), an excess of the basis that the residual correction was neither designed nor validated for. Its `angle` p-values (astrocyte 0.714 / 0.202, microglia 0.868 / 0.659 in `cohort_30` / `cohort_50`) and `shape` p-values (0.321 / 0.395, 0.355 / 0.17) are therefore neither a rejection nor an informative non-rejection; no reading is drawn from them. Its `delta` is read under the binding rules: not rejected in any cohort of either cell type (astrocyte 0.197 / 0.535, microglia 0.735 / 0.884). The arm stays in `LATENT_SENSITIVITY_ARMS` as a sensitivity analysis. The corrected null's size at eight covariates (r = 12) was not measured by MOTCO; it is inferred from the mechanism.

### 7.1 microglia, cohort_30 (≥ 30 microglia nuclei in both omics, n = 76) — primary cohort

#### MOTCO latent-adjustment sensitivity: cohort_30

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 2, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

##### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 5.10924 | 0.705 | no | 10.8174 | 8.13722 | 26.9672 |
| extended_covariates | 8 | none | 2 | 6.89041 | 0.574 | no | 9.51245 | 7.01675 | 22.5416 |
| latent_adjusted | 8 | stage | 3 | 2.43939 | 0.735 | no | 5.50034 | 4.1478 | 13.5441 |

##### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 12.5272 | 0.456 | no | 13.6061 | 10.1863 | 32.1219 |
| extended_covariates | 8 | none | 2 | 5.24863 | 0.7 | no | 11.2635 | 8.53848 | 26.9175 |
| latent_adjusted | 8 | stage | 3 | 5.99532 | 0.868 | no | 16.8715 | 10.5818 | 36.6757 |

##### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 0.0242177 | 0.991 | no | 0.272507 | 0.148066 | 0.544 |
| extended_covariates | 8 | none | 2 | 0.192379 | 0.567 | no | 0.225934 | 0.116701 | 0.445358 |
| latent_adjusted | 8 | stage | 3 | 0.130018 | 0.355 | no | 0.113833 | 0.0604416 | 0.216083 |

##### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 2 | 0.6852 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 107.7 | `results/microglia/motco_analysis/cohort_30` |
| extended_covariates | cv | 2 | 0.6852 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 104.3 | `extended_covariates/motco_analysis/cohort_30` |
| latent_adjusted | cv | 3 | 0.6575 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 111.9 | `latent_adjusted/motco_analysis/cohort_30` |

##### Latent adjustment

MOTCO's `latent_adjustment` block: the residualization's degrees of freedom on the full sample (`n − 1 − k − q`) and the smallest over the CV training sets, and the quantiles over features of the partial R² of the covariates given the intercept and the conditioning indicators (pooled-standardized features).

| arm | conditioning | residual df | smallest CV training df | R² q05 | R² q25 | R² q50 | R² q75 | R² q95 | R² mean | R² fraction > 0.5 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| latent_adjusted | stage | 65 | 27 | 0.0772 | 0.1414 | 0.2005 | 0.2697 | 0.4172 | 0.2158 | 0.0173 |

##### Mediator check

Spearman correlation of each added covariate with stage (ordinal 0, 1, 2) over all cohort donors and within each sex, written before the arms ran; |ρ| ≥ 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`) flags a mediator candidate (flagged rows first). No p-value is computed, and a flag removes no covariate and changes no arm: no simulation can tell a confounder from a mediator, so a flagged fraction may carry stage signal that the `latent_adjusted` basis removed (and that the covariate adjustment of both arms absorbed).

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `mitochondrial_pseudobulk` | all | 76 | 0.498 | yes |
| `mitochondrial_pseudobulk` | F | 45 | 0.577 | yes |
| `mitochondrial_pseudobulk` | M | 31 | 0.329 | yes |
| `ribosomal_pseudobulk` | all | 76 | 0.053 | no |
| `ribosomal_pseudobulk` | F | 45 | -0.030 | no |
| `ribosomal_pseudobulk` | M | 31 | 0.127 | no |

##### Findings

- CV-selected rank differs: 2 in extended_covariates, 3 in latent_adjusted; reported beside angle, not read as better or worse
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.498 with stage (all, n = 76); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.577 with stage (F, n = 45); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.329 with stage (M, n = 31); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)

### 7.2 microglia, cohort_50 (≥ 50 microglia nuclei in both omics, n = 63) — sensitivity cohort

#### MOTCO latent-adjustment sensitivity: cohort_50

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 4, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

##### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 4 | 1.6686 | 0.834 | no | 6.37019 | 4.88809 | 15.3753 |
| extended_covariates | 8 | none | 4 | 1.55793 | 0.839 | no | 6.0961 | 4.63177 | 14.9806 |
| latent_adjusted | 8 | stage | 3 | 1.11008 | 0.884 | no | 6.58876 | 5.17931 | 16.6025 |

##### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 4 | 22.2604 | 0.291 | no | 18.0464 | 8.60842 | 34.1958 |
| extended_covariates | 8 | none | 4 | 10.8156 | 0.7 | no | 14.8633 | 6.72838 | 27.037 |
| latent_adjusted | 8 | stage | 3 | 12.1741 | 0.659 | no | 18.5678 | 11.5285 | 40.6364 |

##### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 4 | 0.0797433 | 0.644 | no | 0.106788 | 0.0560953 | 0.21359 |
| extended_covariates | 8 | none | 4 | 0.0750599 | 0.606 | no | 0.0937548 | 0.0495105 | 0.183592 |
| latent_adjusted | 8 | stage | 3 | 0.176571 | 0.17 | no | 0.121459 | 0.0626905 | 0.236358 |

##### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 4 | 0.6565 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 96.5 | `results/microglia/motco_analysis/cohort_50` |
| extended_covariates | cv | 4 | 0.6565 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 97.9 | `extended_covariates/motco_analysis/cohort_50` |
| latent_adjusted | cv | 3 | 0.6261 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 112.0 | `latent_adjusted/motco_analysis/cohort_50` |

##### Latent adjustment

MOTCO's `latent_adjustment` block: the residualization's degrees of freedom on the full sample (`n − 1 − k − q`) and the smallest over the CV training sets, and the quantiles over features of the partial R² of the covariates given the intercept and the conditioning indicators (pooled-standardized features).

| arm | conditioning | residual df | smallest CV training df | R² q05 | R² q25 | R² q50 | R² q75 | R² q95 | R² mean | R² fraction > 0.5 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| latent_adjusted | stage | 52 | 20 | 0.0813 | 0.1480 | 0.2078 | 0.2794 | 0.4109 | 0.2217 | 0.0141 |

##### Mediator check

Spearman correlation of each added covariate with stage (ordinal 0, 1, 2) over all cohort donors and within each sex, written before the arms ran; |ρ| ≥ 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`) flags a mediator candidate (flagged rows first). No p-value is computed, and a flag removes no covariate and changes no arm: no simulation can tell a confounder from a mediator, so a flagged fraction may carry stage signal that the `latent_adjusted` basis removed (and that the covariate adjustment of both arms absorbed).

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `mitochondrial_pseudobulk` | all | 63 | 0.440 | yes |
| `mitochondrial_pseudobulk` | F | 36 | 0.541 | yes |
| `ribosomal_pseudobulk` | all | 63 | -0.029 | no |
| `ribosomal_pseudobulk` | F | 36 | -0.226 | no |
| `ribosomal_pseudobulk` | M | 27 | 0.158 | no |
| `mitochondrial_pseudobulk` | M | 27 | 0.254 | no |

##### Findings

- CV-selected rank differs: 4 in extended_covariates, 3 in latent_adjusted; reported beside angle, not read as better or worse
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.440 with stage (all, n = 63); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.541 with stage (F, n = 36); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)

### 7.3 astrocyte, cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

#### MOTCO latent-adjustment sensitivity: cohort_30

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 5, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

##### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| extended_covariates | 8 | none | 5 | 0.492704 | 0.954 | no | 6.48472 | 5.05692 | 16.4666 |
| latent_adjusted | 8 | stage | 4 | 6.88452 | 0.197 | no | 4.34258 | 3.35709 | 10.7608 |

##### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| extended_covariates | 8 | none | 5 | 21.5237 | 0.368 | no | 19.2029 | 7.54841 | 32.5527 |
| latent_adjusted | 8 | stage | 4 | 14.0475 | 0.714 | no | 20.5718 | 9.63653 | 38.1246 |

##### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| extended_covariates | 8 | none | 5 | 0.156527 | 0.18 | no | 0.107141 | 0.0562527 | 0.208107 |
| latent_adjusted | 8 | stage | 4 | 0.105552 | 0.321 | no | 0.0869865 | 0.045953 | 0.16805 |

##### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 5 | 0.6356 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 113.8 | `results/motco_analysis/cohort_30` |
| extended_covariates | cv | 5 | 0.6356 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 169.1 | `extended_covariates/motco_analysis/cohort_30` |
| latent_adjusted | cv | 4 | 0.6433 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 175.1 | `latent_adjusted/motco_analysis/cohort_30` |

##### Latent adjustment

MOTCO's `latent_adjustment` block: the residualization's degrees of freedom on the full sample (`n − 1 − k − q`) and the smallest over the CV training sets, and the quantiles over features of the partial R² of the covariates given the intercept and the conditioning indicators (pooled-standardized features).

| arm | conditioning | residual df | smallest CV training df | R² q05 | R² q25 | R² q50 | R² q75 | R² q95 | R² mean | R² fraction > 0.5 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| latent_adjusted | stage | 69 | 29 | 0.0862 | 0.1547 | 0.2230 | 0.3192 | 0.5743 | 0.2578 | 0.0907 |

##### Mediator check

Spearman correlation of each added covariate with stage (ordinal 0, 1, 2) over all cohort donors and within each sex, written before the arms ran; |ρ| ≥ 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`) flags a mediator candidate (flagged rows first). No p-value is computed, and a flag removes no covariate and changes no arm: no simulation can tell a confounder from a mediator, so a flagged fraction may carry stage signal that the `latent_adjusted` basis removed (and that the covariate adjustment of both arms absorbed).

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `ribosomal_pseudobulk` | all | 80 | -0.399 | yes |
| `ribosomal_pseudobulk` | F | 49 | -0.432 | yes |
| `mitochondrial_pseudobulk` | all | 80 | 0.312 | yes |
| `mitochondrial_pseudobulk` | F | 49 | 0.342 | yes |
| `ribosomal_pseudobulk` | M | 31 | -0.285 | no |
| `mitochondrial_pseudobulk` | M | 31 | 0.189 | no |

##### Findings

- CV-selected rank differs: 5 in extended_covariates, 4 in latent_adjusted; reported beside angle, not read as better or worse
- mediator candidate: ribosomal_pseudobulk has Spearman -0.399 with stage (all, n = 80); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: ribosomal_pseudobulk has Spearman -0.432 with stage (F, n = 49); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.312 with stage (all, n = 80); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.342 with stage (F, n = 49); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)

### 7.4 astrocyte, cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74) — sensitivity cohort

#### MOTCO latent-adjustment sensitivity: cohort_50

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 2, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

##### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| extended_covariates | 8 | none | 2 | 0.416021 | 0.967 | no | 7.9949 | 6.44495 | 20.5372 |
| latent_adjusted | 8 | stage | 3 | 4.96136 | 0.535 | no | 6.59857 | 4.9744 | 16.0365 |

##### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| extended_covariates | 8 | none | 2 | 23.6532 | 0.109 | no | 11.8711 | 8.5424 | 28.1463 |
| latent_adjusted | 8 | stage | 3 | 24.7286 | 0.202 | no | 16.9351 | 9.81851 | 35.058 |

##### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| extended_covariates | 8 | none | 2 | 0.151547 | 0.794 | no | 0.281622 | 0.147201 | 0.558166 |
| latent_adjusted | 8 | stage | 3 | 0.132168 | 0.395 | no | 0.121001 | 0.0634925 | 0.239543 |

##### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 2 | 0.6965 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 102.3 | `results/motco_analysis/cohort_50` |
| extended_covariates | cv | 2 | 0.6965 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 200.2 | `extended_covariates/motco_analysis/cohort_50` |
| latent_adjusted | cv | 3 | 0.6645 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 203.9 | `latent_adjusted/motco_analysis/cohort_50` |

##### Latent adjustment

MOTCO's `latent_adjustment` block: the residualization's degrees of freedom on the full sample (`n − 1 − k − q`) and the smallest over the CV training sets, and the quantiles over features of the partial R² of the covariates given the intercept and the conditioning indicators (pooled-standardized features).

| arm | conditioning | residual df | smallest CV training df | R² q05 | R² q25 | R² q50 | R² q75 | R² q95 | R² mean | R² fraction > 0.5 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| latent_adjusted | stage | 63 | 25 | 0.0906 | 0.1613 | 0.2350 | 0.3406 | 0.6031 | 0.2721 | 0.1082 |

##### Mediator check

Spearman correlation of each added covariate with stage (ordinal 0, 1, 2) over all cohort donors and within each sex, written before the arms ran; |ρ| ≥ 0.3 (`DIAGNOSTIC_CORRELATION_THRESHOLD`) flags a mediator candidate (flagged rows first). No p-value is computed, and a flag removes no covariate and changes no arm: no simulation can tell a confounder from a mediator, so a flagged fraction may carry stage signal that the `latent_adjusted` basis removed (and that the covariate adjustment of both arms absorbed).

| covariate | subset | n | Spearman | flagged |
|---|---|---:|---:|:---:|
| `ribosomal_pseudobulk` | all | 74 | -0.352 | yes |
| `ribosomal_pseudobulk` | F | 44 | -0.360 | yes |
| `mitochondrial_pseudobulk` | all | 74 | 0.356 | yes |
| `mitochondrial_pseudobulk` | F | 44 | 0.422 | yes |
| `ribosomal_pseudobulk` | M | 30 | -0.286 | no |
| `mitochondrial_pseudobulk` | M | 30 | 0.189 | no |

##### Findings

- CV-selected rank differs: 2 in extended_covariates, 3 in latent_adjusted; reported beside angle, not read as better or worse
- mediator candidate: ribosomal_pseudobulk has Spearman -0.352 with stage (all, n = 74); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: ribosomal_pseudobulk has Spearman -0.360 with stage (F, n = 44); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.356 with stage (all, n = 74); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)
- mediator candidate: mitochondrial_pseudobulk has Spearman 0.422 with stage (F, n = 44); the latent_adjusted basis may have removed stage signal carried by it (no simulation can tell a confounder from a mediator; the flag changes no arm)

### 7.5 Superseded → corrected (arms)

| cell type | cohort | analysis | rank | CV AUROC | κ (n, reduced rank) | statistic | observed | p superseded → corrected | change | rejected superseded → corrected |
|---|---|---|---:|---:|---|---|---:|---|---|---|
| astrocyte | cohort_30 | extended_covariates | 5 | 0.6356 | 1.0847 (80, 12) | delta | 0.492704 | 0.942 → 0.954 | up | no → no |
| astrocyte | cohort_30 | extended_covariates |  |  |  | angle | 21.5237 | 0.286 → 0.368 | up | no → no |
| astrocyte | cohort_30 | extended_covariates |  |  |  | shape | 0.156527 | 0.137 → 0.18 | up | no → no |
| astrocyte | cohort_50 | extended_covariates | 2 | 0.6965 | 1.0925 (74, 12) | delta | 0.416021 | 0.961 → 0.967 | up | no → no |
| astrocyte | cohort_50 | extended_covariates |  |  |  | angle | 23.6532 | 0.075 → 0.109 | up | no → no |
| astrocyte | cohort_50 | extended_covariates |  |  |  | shape | 0.151547 | 0.764 → 0.794 | up | no → no |
| astrocyte | cohort_30 | latent_adjusted | 4 | 0.6433 | 1.0847 (80, 12) | delta | 6.88452 | 0.162 → 0.197 | up | no → no |
| astrocyte | cohort_30 | latent_adjusted |  |  |  | angle | 14.0475 | 0.673 → 0.714 | up | no → no |
| astrocyte | cohort_30 | latent_adjusted |  |  |  | shape | 0.105552 | 0.256 → 0.321 | up | no → no |
| astrocyte | cohort_50 | latent_adjusted | 3 | 0.6645 | 1.0925 (74, 12) | delta | 4.96136 | 0.494 → 0.535 | up | no → no |
| astrocyte | cohort_50 | latent_adjusted |  |  |  | angle | 24.7286 | 0.154 → 0.202 | up | no → no |
| astrocyte | cohort_50 | latent_adjusted |  |  |  | shape | 0.132168 | 0.337 → 0.395 | up | no → no |
| microglia | cohort_30 | extended_covariates | 2 | 0.6852 | 1.0897 (76, 12) | delta | 6.89041 | 0.546 → 0.574 | up | no → no |
| microglia | cohort_30 | extended_covariates |  |  |  | angle | 5.24863 | 0.677 → 0.7 | up | no → no |
| microglia | cohort_30 | extended_covariates |  |  |  | shape | 0.192379 | 0.509 → 0.567 | up | no → no |
| microglia | cohort_50 | extended_covariates | 4 | 0.6565 | 1.1114 (63, 12) | delta | 1.55793 | 0.825 → 0.839 | up | no → no |
| microglia | cohort_50 | extended_covariates |  |  |  | angle | 10.8156 | 0.616 → 0.7 | up | no → no |
| microglia | cohort_50 | extended_covariates |  |  |  | shape | 0.0750599 | 0.537 → 0.606 | up | no → no |
| microglia | cohort_30 | latent_adjusted | 3 | 0.6575 | 1.0897 (76, 12) | delta | 2.43939 | 0.715 → 0.735 | up | no → no |
| microglia | cohort_30 | latent_adjusted |  |  |  | angle | 5.99532 | 0.849 → 0.868 | up | no → no |
| microglia | cohort_30 | latent_adjusted |  |  |  | shape | 0.130018 | 0.294 → 0.355 | up | no → no |
| microglia | cohort_50 | latent_adjusted | 3 | 0.6261 | 1.1114 (63, 12) | delta | 1.11008 | 0.875 → 0.884 | up | no → no |
| microglia | cohort_50 | latent_adjusted |  |  |  | angle | 12.1741 | 0.611 → 0.659 | up | no → no |
| microglia | cohort_50 | latent_adjusted |  |  |  | shape | 0.176571 | 0.114 → 0.17 | up | no → no |

Superseded against corrected, generated by the comparison script over 8 analyses: of the 32 analysis CSVs (four per analysis), 32 equal at rounding level (largest |difference| 9.7e-13); of the 8 records outside the null keys, 8 equal at rounding level (largest relative difference 5.3e-11, apart from the numerically-zero third eigenvalues of the configuration spectra, ≤ 1.4e-31); rank equal in 8 of 8, CV AUROC equal in 8 of 8; κ = √(n / (n − r)) to 1e-12 in 8 of 8; p-values up in 24, equal in 0, down in 0 of 24; rejection statuses changed: 0.

## 8. Stage DE (DESeq2 and limma-voom)

The stage-DE fits do not depend on MOTCO; they were refitted so that every output rests on one pin, and their concordance cells now read the corrected primaries' `delta` status (not rejected in any). Every concordance cell, DESeq2 (`factor_interaction` union, exported) and limma (`interaction_joint`, exported), in both omics, both cohorts and both cell types, and on the exclusion roots, is concordant with 0 exported interaction features, and the methods agree. Every exported stage count is unchanged, including the counts E1 reads (microglia 61 against astrocyte 1 under DESeq2, 2 against 1 under limma, in `cohort_30`).

### 8.1 Concordance cells (from `stage_de_summary.csv`, `stage_de_limma_summary.csv` and `method_comparison.csv`), corrected

| cell type | cohort | root | method | omic | cell | n significant | delta rejected | concordance | superseded |
|---|---|---|---|---|---|---|---|---|---|
| astrocyte | cohort_30 | primary | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant | same |
| astrocyte | cohort_30 | primary | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant | same |
| astrocyte | cohort_30 | primary | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| astrocyte | cohort_30 | primary | limma | atac | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| astrocyte | cohort_50 | primary | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant | same |
| astrocyte | cohort_50 | primary | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant | same |
| astrocyte | cohort_50 | primary | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| astrocyte | cohort_50 | primary | limma | atac | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant | same |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant | same |
| astrocyte | cohort_30 | excluding_severely_affected | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant | same |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant | same |
| astrocyte | cohort_50 | excluding_severely_affected | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| microglia | cohort_30 | primary | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant | same |
| microglia | cohort_30 | primary | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant | same |
| microglia | cohort_30 | primary | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| microglia | cohort_30 | primary | limma | atac | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| microglia | cohort_50 | primary | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant | same |
| microglia | cohort_50 | primary | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant | same |
| microglia | cohort_50 | primary | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| microglia | cohort_50 | primary | limma | atac | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant | same |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant | same |
| microglia | cohort_30 | excluding_severely_affected | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant | same |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant | same |
| microglia | cohort_50 | excluding_severely_affected | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (methods agree True) | same |

### 8.2 Exported-RNA stage counts at padj ≤ 0.05 (from the summaries), corrected / superseded

| cell type | cohort | root | method | model | coefficient | n significant (exported RNA) | superseded |
|---|---|---|---|---|---|---:|---:|
| microglia | cohort_30 | primary | DESeq2 | factor_additive | stage_1_vs_0 | 0 | 0 |
| microglia | cohort_30 | primary | DESeq2 | factor_additive | stage_2_vs_0 | 61 | 61 |
| microglia | cohort_30 | primary | DESeq2 | factor_additive | union | 61 | 61 |
| microglia | cohort_30 | primary | DESeq2 | trend_additive | stage_num | 57 | 57 |
| microglia | cohort_30 | primary | limma | factor_additive | stage_1_vs_0 | 0 | 0 |
| microglia | cohort_30 | primary | limma | factor_additive | stage_2_vs_0 | 28 | 28 |
| microglia | cohort_30 | primary | limma | factor_additive | stage_joint | 2 | 2 |
| microglia | cohort_30 | primary | limma | factor_additive | union | 28 | 28 |
| microglia | cohort_30 | primary | limma | trend_additive | stage_num | 31 | 31 |
| microglia | cohort_50 | primary | DESeq2 | factor_additive | stage_1_vs_0 | 1 | 1 |
| microglia | cohort_50 | primary | DESeq2 | factor_additive | stage_2_vs_0 | 34 | 34 |
| microglia | cohort_50 | primary | DESeq2 | factor_additive | union | 35 | 35 |
| microglia | cohort_50 | primary | DESeq2 | trend_additive | stage_num | 34 | 34 |
| microglia | cohort_50 | primary | limma | factor_additive | stage_1_vs_0 | 0 | 0 |
| microglia | cohort_50 | primary | limma | factor_additive | stage_2_vs_0 | 36 | 36 |
| microglia | cohort_50 | primary | limma | factor_additive | stage_joint | 0 | 0 |
| microglia | cohort_50 | primary | limma | factor_additive | union | 36 | 36 |
| microglia | cohort_50 | primary | limma | trend_additive | stage_num | 35 | 35 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | stage_1_vs_0 | 0 | 0 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | stage_2_vs_0 | 20 | 20 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | union | 20 | 20 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | trend_additive | stage_num | 21 | 21 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_1_vs_0 | 0 | 0 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_2_vs_0 | 2 | 2 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_joint | 2 | 2 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_additive | union | 2 | 2 |
| microglia | cohort_30 | excluding_severely_affected | limma | trend_additive | stage_num | 2 | 2 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | stage_1_vs_0 | 1 | 1 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | stage_2_vs_0 | 2 | 2 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | union | 3 | 3 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | trend_additive | stage_num | 0 | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_1_vs_0 | 0 | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_2_vs_0 | 0 | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_joint | 0 | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_additive | union | 0 | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | trend_additive | stage_num | 0 | 0 |
| astrocyte | cohort_30 | primary | DESeq2 | factor_additive | stage_1_vs_0 | 0 | 0 |
| astrocyte | cohort_30 | primary | DESeq2 | factor_additive | stage_2_vs_0 | 1 | 1 |
| astrocyte | cohort_30 | primary | DESeq2 | factor_additive | union | 1 | 1 |
| astrocyte | cohort_30 | primary | DESeq2 | trend_additive | stage_num | 2 | 2 |
| astrocyte | cohort_30 | primary | limma | factor_additive | stage_1_vs_0 | 0 | 0 |
| astrocyte | cohort_30 | primary | limma | factor_additive | stage_2_vs_0 | 1 | 1 |
| astrocyte | cohort_30 | primary | limma | factor_additive | stage_joint | 1 | 1 |
| astrocyte | cohort_30 | primary | limma | factor_additive | union | 1 | 1 |
| astrocyte | cohort_30 | primary | limma | trend_additive | stage_num | 1 | 1 |
| astrocyte | cohort_50 | primary | DESeq2 | factor_additive | stage_1_vs_0 | 0 | 0 |
| astrocyte | cohort_50 | primary | DESeq2 | factor_additive | stage_2_vs_0 | 2 | 2 |
| astrocyte | cohort_50 | primary | DESeq2 | factor_additive | union | 2 | 2 |
| astrocyte | cohort_50 | primary | DESeq2 | trend_additive | stage_num | 2 | 2 |
| astrocyte | cohort_50 | primary | limma | factor_additive | stage_1_vs_0 | 0 | 0 |
| astrocyte | cohort_50 | primary | limma | factor_additive | stage_2_vs_0 | 1 | 1 |
| astrocyte | cohort_50 | primary | limma | factor_additive | stage_joint | 1 | 1 |
| astrocyte | cohort_50 | primary | limma | factor_additive | union | 1 | 1 |
| astrocyte | cohort_50 | primary | limma | trend_additive | stage_num | 1 | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | stage_1_vs_0 | 0 | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | stage_2_vs_0 | 1 | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | union | 1 | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | trend_additive | stage_num | 1 | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_1_vs_0 | 0 | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_2_vs_0 | 1 | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_joint | 1 | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_additive | union | 1 | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | trend_additive | stage_num | 1 | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | stage_1_vs_0 | 0 | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | stage_2_vs_0 | 1 | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | union | 1 | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | trend_additive | stage_num | 1 | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_1_vs_0 | 0 | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_2_vs_0 | 1 | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_joint | 1 | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_additive | union | 1 | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | trend_additive | stage_num | 1 | 1 |

### 8.3 Stage-DE tables against their superseded copies (generated)

Generated by the comparison script; the superseded tables were fitted on `nixos-home`, these on `nixos-dell-laptop`. The 96 limma tables are byte-identical (48) or equal at rounding level (48; largest |difference| 2.9e-12). The 96 DESeq2 tables differ beyond rounding, as across machines in the 2026-10-26 re-issue: excluding degenerate fits (|Wald statistic| > 1e6 in either copy), the largest |Δ log2FoldChange| is 0.0872, |Δ stat| 0.527 and |Δ padj| 0.0636. Every exported padj ≤ 0.05 call is equal in all 96 tables; all calls are equal in 94 (the other two in section 12). Per table:

| table | max |Δ log2FoldChange| | max |Δ lfcSE| | max |Δ stat| | max |Δ pvalue| | max |Δ padj| | calls equal (exported) | n calls (superseded, corrected) |
|---|---|---|---|---|---|---|---|
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/factor_additive__stage_1_vs_0.csv` (degenerate fits superseded 2, corrected 2) | 0.00159 | 0.00481 | 0.0718 | 0.0375 | 1.79e-08 | yes (yes) | 3, 3 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/factor_additive__stage_2_vs_0.csv` (degenerate fits superseded 2, corrected 2) | 0.00145 | 0.00434 | 0.0767 | 0.0375 | 0.0202 | yes (yes) | 8, 8 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/factor_interaction__male_x_stage_1.csv` (degenerate fits superseded 1, corrected 1) | 0.00632 | 0.00951 | 0.0529 | 0.0259 | 5.21e-10 | yes (yes) | 1, 1 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/factor_interaction__male_x_stage_2.csv` (degenerate fits superseded 1, corrected 1) | 0.0059 | 0.00882 | 0.0711 | 0.0299 | 1.59e-09 | yes (yes) | 1, 1 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/trend_additive__stage_num.csv` (degenerate fits superseded 1, corrected 1) | 0.000428 | 0.0021 | 0.0616 | 0.0354 | 0.0165 | yes (yes) | 11, 11 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/trend_interaction__male_x_stage_num.csv` | 0.00184 | 0.00374 | 0.025 | 0.0155 | 1.21e-10 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/factor_additive__stage_1_vs_0.csv` | 0.000246 | 0.00237 | 0.0739 | 0.0152 | 0.00969 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/factor_additive__stage_2_vs_0.csv` | 0.000577 | 0.00218 | 0.0601 | 0.0138 | 0.0153 | yes (yes) | 45, 45 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/factor_interaction__male_x_stage_1.csv` | 0.000517 | 0.00541 | 0.0333 | 0.0194 | 1.26e-08 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/factor_interaction__male_x_stage_2.csv` | 0.00168 | 0.00508 | 0.0644 | 0.0197 | 0.00292 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/trend_additive__stage_num.csv` | 0.000234 | 0.00125 | 0.104 | 0.0124 | 0.0261 | yes (yes) | 46, 46 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/trend_interaction__male_x_stage_num.csv` | 0.000749 | 0.00745 | 0.202 | 0.0448 | 0.0225 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/factor_additive__stage_1_vs_0.csv` | 0.00412 | 0.0062 | 0.0462 | 0.0183 | 5.2e-08 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/factor_additive__stage_2_vs_0.csv` | 0.00106 | 0.00576 | 0.0719 | 0.0402 | 0.0034 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/factor_interaction__male_x_stage_1.csv` | 0.003 | 0.00918 | 0.0715 | 0.0178 | 1.25e-07 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/factor_interaction__male_x_stage_2.csv` | 0.0024 | 0.00866 | 0.0532 | 0.0395 | 2.15e-09 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/trend_additive__stage_num.csv` | 0.00123 | 0.00299 | 0.0897 | 0.0527 | 0.0262 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/trend_interaction__male_x_stage_num.csv` (degenerate fits superseded 1, corrected 1) | 0.00164 | 0.00652 | 0.0832 | 0.0224 | 3.37e-09 | yes (yes) | 1, 1 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/factor_additive__stage_1_vs_0.csv` | 0.000226 | 0.00319 | 0.0827 | 0.0222 | 0.014 | yes (yes) | 1, 1 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/factor_additive__stage_2_vs_0.csv` | 0.000346 | 0.00307 | 0.0589 | 0.0204 | 0.0228 | yes (yes) | 10, 10 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/factor_interaction__male_x_stage_1.csv` | 0.00169 | 0.00706 | 0.00492 | 0.00369 | 2.37e-09 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/factor_interaction__male_x_stage_2.csv` | 0.000275 | 0.00674 | 0.0234 | 0.0161 | 1.9e-10 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/trend_additive__stage_num.csv` | 0.000445 | 0.00163 | 0.109 | 0.0221 | 0.0317 | yes (yes) | 0, 0 |
| `microglia/motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/trend_interaction__male_x_stage_num.csv` | 0.0005 | 0.00316 | 0.0766 | 0.0115 | 2.31e-08 | yes (yes) | 0, 0 |
| `microglia/stage_de/cohort_30/atac/factor_additive__stage_1_vs_0.csv` (degenerate fits superseded 1, corrected 1) | 0.00769 | 0.225 | 0.126 | 0.0346 | 3.51e-05 | yes (yes) | 1, 1 |
| `microglia/stage_de/cohort_30/atac/factor_additive__stage_2_vs_0.csv` (degenerate fits superseded 1, corrected 1) | 0.0218 | 0.2 | 0.133 | 0.0697 | 0.0384 | yes (yes) | 113, 113 |
| `microglia/stage_de/cohort_30/atac/factor_interaction__male_x_stage_1.csv` (degenerate fits superseded 3, corrected 3) | 0.00504 | 0.0111 | 0.129 | 0.0566 | 5.31e-08 | yes (yes) | 3, 3 |
| `microglia/stage_de/cohort_30/atac/factor_interaction__male_x_stage_2.csv` (degenerate fits superseded 3, corrected 3) | 0.00475 | 0.00993 | 0.0528 | 0.0241 | 8.66e-06 | yes (yes) | 3, 3 |
| `microglia/stage_de/cohort_30/atac/trend_additive__stage_num.csv` | 0.00181 | 0.00181 | 0.048 | 0.0255 | 0.0176 | yes (yes) | 116, 116 |
| `microglia/stage_de/cohort_30/atac/trend_interaction__male_x_stage_num.csv` (degenerate fits superseded 1, corrected 1) | 0.0872 | 0.00388 | 0.118 | 0.065 | 5.97e-06 | yes (yes) | 1, 1 |
| `microglia/stage_de/cohort_30/rna/factor_additive__stage_1_vs_0.csv` | 0.00171 | 0.00286 | 0.0934 | 0.0215 | 0.0107 | yes (yes) | 0, 0 |
| `microglia/stage_de/cohort_30/rna/factor_additive__stage_2_vs_0.csv` | 0.00108 | 0.00257 | 0.175 | 0.016 | 0.023 | yes (yes) | 156, 156 |
| `microglia/stage_de/cohort_30/rna/factor_interaction__male_x_stage_1.csv` | 0.00158 | 0.00629 | 0.0927 | 0.0239 | 3.18e-08 | yes (yes) | 0, 0 |
| `microglia/stage_de/cohort_30/rna/factor_interaction__male_x_stage_2.csv` | 0.0029 | 0.00564 | 0.047 | 0.0184 | 0.00803 | yes (yes) | 0, 0 |
| `microglia/stage_de/cohort_30/rna/trend_additive__stage_num.csv` | 0.00082 | 0.00201 | 0.108 | 0.0302 | 0.041 | yes (yes) | 145, 145 |
| `microglia/stage_de/cohort_30/rna/trend_interaction__male_x_stage_num.csv` | 0.00085 | 0.00269 | 0.0842 | 0.0243 | 0.0309 | yes (yes) | 0, 0 |
| `microglia/stage_de/cohort_50/atac/factor_additive__stage_1_vs_0.csv` (degenerate fits superseded 2, corrected 2) | 0.00293 | 0.0074 | 0.146 | 0.0204 | 3.87e-09 | yes (yes) | 3, 3 |
| `microglia/stage_de/cohort_50/atac/factor_additive__stage_2_vs_0.csv` (degenerate fits superseded 2, corrected 2) | 0.00425 | 0.00665 | 0.0247 | 0.0171 | 0.00743 | yes (yes) | 121, 121 |
| `microglia/stage_de/cohort_50/atac/factor_interaction__male_x_stage_1.csv` (degenerate fits superseded 2, corrected 3) | 0.00811 | 0.0136 | 0.214 | 0.0649 | 1.19e-07 | no (yes) | 3, 4 |
| `microglia/stage_de/cohort_50/atac/factor_interaction__male_x_stage_2.csv` (degenerate fits superseded 2, corrected 3) | 0.00779 | 0.0128 | 0.146 | 0.0627 | 3.9e-06 | no (yes) | 3, 4 |
| `microglia/stage_de/cohort_50/atac/trend_additive__stage_num.csv` (degenerate fits superseded 2, corrected 2) | 0.00174 | 0.00264 | 0.148 | 0.0199 | 0.0468 | yes (yes) | 130, 130 |
| `microglia/stage_de/cohort_50/atac/trend_interaction__male_x_stage_num.csv` (degenerate fits superseded 3, corrected 3) | 0.00259 | 0.00936 | 0.0749 | 0.046 | 1.29e-10 | yes (yes) | 3, 3 |
| `microglia/stage_de/cohort_50/rna/factor_additive__stage_1_vs_0.csv` | 0.00195 | 0.00405 | 0.0662 | 0.0258 | 0.0239 | yes (yes) | 1, 1 |
| `microglia/stage_de/cohort_50/rna/factor_additive__stage_2_vs_0.csv` | 0.00212 | 0.00379 | 0.0736 | 0.0176 | 0.021 | yes (yes) | 64, 64 |
| `microglia/stage_de/cohort_50/rna/factor_interaction__male_x_stage_1.csv` | 0.00325 | 0.00873 | 0.0895 | 0.0247 | 4.81e-08 | yes (yes) | 0, 0 |
| `microglia/stage_de/cohort_50/rna/factor_interaction__male_x_stage_2.csv` | 0.00567 | 0.00833 | 0.0731 | 0.0339 | 0.00193 | yes (yes) | 0, 0 |
| `microglia/stage_de/cohort_50/rna/trend_additive__stage_num.csv` | 0.000509 | 0.00183 | 0.0774 | 0.0404 | 0.0264 | yes (yes) | 59, 59 |
| `microglia/stage_de/cohort_50/rna/trend_interaction__male_x_stage_num.csv` | 0.00103 | 0.00354 | 0.0926 | 0.0128 | 0.00507 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/factor_additive__stage_1_vs_0.csv` | 0.00546 | 0.0119 | 0.527 | 0.0835 | 0.0636 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/factor_additive__stage_2_vs_0.csv` | 0.00548 | 0.0103 | 0.31 | 0.112 | 0.046 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/factor_interaction__male_x_stage_1.csv` | 0.00685 | 0.0176 | 0.112 | 0.0707 | 0.00129 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/factor_interaction__male_x_stage_2.csv` | 0.00746 | 0.0161 | 0.182 | 0.0706 | 0.00794 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/trend_additive__stage_num.csv` | 0.00219 | 0.00363 | 0.131 | 0.0895 | 0.0208 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/atac/trend_interaction__male_x_stage_num.csv` | 0.00416 | 0.00654 | 0.25 | 0.0673 | 0.00564 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/factor_additive__stage_1_vs_0.csv` | 0.000496 | 0.00357 | 0.0711 | 0.0186 | 0.000111 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/factor_additive__stage_2_vs_0.csv` | 0.00118 | 0.00318 | 0.0717 | 0.0332 | 0.0145 | yes (yes) | 1, 1 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/factor_interaction__male_x_stage_1.csv` | 0.00981 | 0.0115 | 0.164 | 0.079 | 2.04e-07 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/factor_interaction__male_x_stage_2.csv` | 0.00331 | 0.0105 | 0.124 | 0.0298 | 1.98e-07 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/trend_additive__stage_num.csv` | 0.000273 | 0.00141 | 0.124 | 0.0178 | 0.0152 | yes (yes) | 1, 1 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_30/rna/trend_interaction__male_x_stage_num.csv` | 0.00135 | 0.00271 | 0.114 | 0.0195 | 0.000202 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/factor_additive__stage_1_vs_0.csv` | 0.00589 | 0.00858 | 0.198 | 0.0648 | 0.0127 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/factor_additive__stage_2_vs_0.csv` | 0.00388 | 0.00721 | 0.235 | 0.0542 | 0.0497 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/factor_interaction__male_x_stage_1.csv` | 0.0055 | 0.015 | 0.151 | 0.0719 | 4.21e-08 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/factor_interaction__male_x_stage_2.csv` | 0.00474 | 0.0136 | 0.216 | 0.0535 | 1.32e-07 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/trend_additive__stage_num.csv` | 0.00429 | 0.00397 | 0.16 | 0.119 | 0.0328 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/atac/trend_interaction__male_x_stage_num.csv` | 0.00418 | 0.00864 | 0.166 | 0.0311 | 0.000624 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/factor_additive__stage_1_vs_0.csv` | 0.000974 | 0.00299 | 0.118 | 0.0139 | 0.0156 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/factor_additive__stage_2_vs_0.csv` | 0.00114 | 0.00265 | 0.0984 | 0.0189 | 0.0181 | yes (yes) | 1, 1 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/factor_interaction__male_x_stage_1.csv` | 0.00216 | 0.00601 | 0.0571 | 0.0204 | 6.39e-08 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/factor_interaction__male_x_stage_2.csv` | 0.0019 | 0.00547 | 0.0828 | 0.0226 | 7e-09 | yes (yes) | 0, 0 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/trend_additive__stage_num.csv` | 0.000564 | 0.00155 | 0.0552 | 0.0243 | 0.0177 | yes (yes) | 1, 1 |
| `motco_exclusion_sensitivity/excluding_severely_affected/stage_de/cohort_50/rna/trend_interaction__male_x_stage_num.csv` | 0.0018 | 0.00361 | 0.0916 | 0.032 | 0.000271 | yes (yes) | 0, 0 |
| `stage_de/cohort_30/atac/factor_additive__stage_1_vs_0.csv` | 0.00312 | 0.00585 | 0.117 | 0.0716 | 0.00398 | yes (yes) | 0, 0 |
| `stage_de/cohort_30/atac/factor_additive__stage_2_vs_0.csv` | 0.00333 | 0.00501 | 0.16 | 0.0733 | 0.0419 | yes (yes) | 22, 22 |
| `stage_de/cohort_30/atac/factor_interaction__male_x_stage_1.csv` | 0.00493 | 0.0148 | 0.188 | 0.0509 | 8.69e-08 | yes (yes) | 0, 0 |
| `stage_de/cohort_30/atac/factor_interaction__male_x_stage_2.csv` | 0.00822 | 0.013 | 0.289 | 0.0616 | 0.0242 | yes (yes) | 0, 0 |
| `stage_de/cohort_30/atac/trend_additive__stage_num.csv` | 0.00231 | 0.00348 | 0.163 | 0.0664 | 0.0449 | yes (yes) | 23, 23 |
| `stage_de/cohort_30/atac/trend_interaction__male_x_stage_num.csv` | 0.00317 | 0.00918 | 0.182 | 0.0686 | 0.0101 | yes (yes) | 0, 0 |
| `stage_de/cohort_30/rna/factor_additive__stage_1_vs_0.csv` | 0.000887 | 0.00236 | 0.0415 | 0.0153 | 0.000655 | yes (yes) | 0, 0 |
| `stage_de/cohort_30/rna/factor_additive__stage_2_vs_0.csv` | 0.0011 | 0.00206 | 0.092 | 0.0183 | 0.0161 | yes (yes) | 1, 1 |
| `stage_de/cohort_30/rna/factor_interaction__male_x_stage_1.csv` | 0.000937 | 0.00546 | 0.0343 | 0.025 | 4.48e-05 | yes (yes) | 0, 0 |
| `stage_de/cohort_30/rna/factor_interaction__male_x_stage_2.csv` | 0.00195 | 0.0048 | 0.0506 | 0.0302 | 0.00992 | yes (yes) | 0, 0 |
| `stage_de/cohort_30/rna/trend_additive__stage_num.csv` | 0.000434 | 0.291 | 0.231 | 0.161 | 0.0582 | yes (yes) | 2, 2 |
| `stage_de/cohort_30/rna/trend_interaction__male_x_stage_num.csv` | 0.00154 | 0.00314 | 0.0711 | 0.0262 | 0.0217 | yes (yes) | 0, 0 |
| `stage_de/cohort_50/atac/factor_additive__stage_1_vs_0.csv` | 0.0029 | 0.00675 | 0.125 | 0.0426 | 7.75e-08 | yes (yes) | 0, 0 |
| `stage_de/cohort_50/atac/factor_additive__stage_2_vs_0.csv` | 0.00428 | 0.0056 | 0.222 | 0.032 | 0.0544 | yes (yes) | 36, 36 |
| `stage_de/cohort_50/atac/factor_interaction__male_x_stage_1.csv` | 0.00865 | 0.0146 | 0.309 | 0.038 | 5e-08 | yes (yes) | 0, 0 |
| `stage_de/cohort_50/atac/factor_interaction__male_x_stage_2.csv` | 0.0119 | 0.0127 | 0.316 | 0.0272 | 0.0174 | yes (yes) | 0, 0 |
| `stage_de/cohort_50/atac/trend_additive__stage_num.csv` | 0.00272 | 0.00328 | 0.143 | 0.063 | 0.039 | yes (yes) | 30, 30 |
| `stage_de/cohort_50/atac/trend_interaction__male_x_stage_num.csv` | 0.00187 | 0.00602 | 0.0896 | 0.0589 | 0.0202 | yes (yes) | 0, 0 |
| `stage_de/cohort_50/rna/factor_additive__stage_1_vs_0.csv` | 0.000621 | 0.00272 | 0.0537 | 0.0221 | 0.00566 | yes (yes) | 0, 0 |
| `stage_de/cohort_50/rna/factor_additive__stage_2_vs_0.csv` | 0.000769 | 0.00233 | 0.0437 | 0.0135 | 0.014 | yes (yes) | 2, 2 |
| `stage_de/cohort_50/rna/factor_interaction__male_x_stage_1.csv` | 0.00751 | 0.013 | 0.172 | 0.0652 | 3.64e-05 | yes (yes) | 0, 0 |
| `stage_de/cohort_50/rna/factor_interaction__male_x_stage_2.csv` | 0.00961 | 0.0111 | 0.128 | 0.02 | 0.0136 | yes (yes) | 0, 0 |
| `stage_de/cohort_50/rna/trend_additive__stage_num.csv` | 0.000194 | 0.00126 | 0.108 | 0.0303 | 0.0278 | yes (yes) | 2, 2 |
| `stage_de/cohort_50/rna/trend_interaction__male_x_stage_num.csv` | 0.00103 | 0.00249 | 0.0593 | 0.0238 | 0.00843 | yes (yes) | 0, 0 |


## 9. Donor-exclusion sensitivity

The variant (`excluding_severely_affected`) was rerun for both cell types with the correction (κ 1.0801 / 1.0887 for astrocyte `cohort_30` / `cohort_50` at n = 70 / 64, 1.0814 / 1.0992 for microglia at n = 69 / 58). The flag check and the variant exports read no p-value and are byte-identical to the superseded ones; `flag_check.csv` was written before each variant export. No variant analysis rejects any statistic, and no rejection status differs from the corrected primary or the superseded variant: no cohort dependence and no dependence on the null correction. The axis check is byte-identical: it holds in microglia, the expectation (`cohort_30` |ρ stage| 0.8201 against `fraction_mitochondrial_umis_sd` 0.6891; `cohort_50` 0.8248 against 0.6746), and fails in astrocytes, descriptive (0.5784 against `ribosomal_mean` 0.7286; 0.5955 against 0.6716). The exclusion roots' stage DE is concordant in every cell. The corrected null's size at the variants' smaller n is inferred, not measured.

### 9.1 microglia, cohort_30 (≥ 30 microglia nuclei in both omics, n = 76) — primary cohort

#### Donor-exclusion sensitivity: microglia, cohort_30, excluding_severely_affected

The `primary` row is the primary analysis (CV-selected rank 2, 76 donors), read from its `result.json` and never rerun; it stays the primary result. The `excluding_severely_affected` row reruns the whole export and the primary analysis on the cohort donors minus those SEA-AD flagged as severely affected (Gabitto et al. 2024: QC-outlying high-ADNC donors excluded from the gene-expression tests): 7 cohort donor(s) excluded, the cohort thresholds unchanged and k recomputed from the remaining cells by the export's rule (k = 9 → 9). The variant is a sensitivity analysis, never a revised primary: the excluded donors carry the most severe pathology, so the exclusion truncates the late epoch along with the quality tail, and the flag was made from QC metrics, so the variant conditions on quality by donor selection. A rejection status that differs between the rows is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in both rows. 999 RRPP permutations with seed 20261019 in both rows; α = 0.05 (rejected when p ≤ α). MOTCO's operating characteristics were measured on neither cohort.

##### Flagged donors per cell

| cell | flagged | unflagged |
|---|---:|---:|
| Female/stage 0 | 0 | 11 |
| Female/stage 1 | 0 | 10 |
| Female/stage 2 | 5 | 19 |
| Male/stage 0 | 0 | 9 |
| Male/stage 1 | 0 | 11 |
| Male/stage 2 | 2 | 9 |
| all | 7 | 69 |

##### Flag check

Written before the variant export ran. The rank-biserial AUC is P(flagged > unflagged) + ½ P(=) over the cohort donors (`auc_all`) and within stage 2 (`auc_stage2`); 0.5 means no separation. Among the unflagged donors, each variable's Spearman correlation with stage and with the primary `pls_0`. Descriptive, no p-value; the variant runs whatever it says.

| variable | source | n flagged | n unflagged | AUC all | AUC stage 2 | ρ stage (unflagged) | ρ pls_0 (unflagged) |
|---|---|---:|---:|---:|---:|---:|---:|
| `pls_0` | latent_scores.csv (primary analysis) | 7 | 69 | 0.986 | 0.964 | 0.643 | n/a |
| `stage` | metadata.csv | 7 | 69 | 0.797 | n/a | n/a | 0.643 |
| `n_nuclei_rna` | donor_technical.csv | 7 | 69 | 0.334 | 0.332 | -0.009 | 0.064 |
| `n_libraries_rna` | donor_technical.csv | 7 | 69 | 0.539 | 0.533 | 0.076 | 0.090 |
| `ribosomal_mean` | donor_technical.csv | 7 | 69 | 0.337 | 0.321 | 0.076 | -0.055 |
| `ribosomal_sd` | donor_technical.csv | 7 | 69 | 0.377 | 0.357 | 0.071 | 0.069 |
| `ribosomal_pseudobulk` | donor_technical.csv | 7 | 69 | 0.447 | 0.418 | 0.079 | 0.034 |
| `ribosomal_supertype_share` | donor_technical.csv | 7 | 69 | 0.832 | 0.770 | 0.292 | 0.311 |
| `mitochondrial_mean` | donor_technical.csv | 7 | 69 | 0.967 | 0.934 | 0.387 | 0.790 |
| `mitochondrial_sd` | donor_technical.csv | 7 | 69 | 0.950 | 0.903 | 0.447 | 0.781 |
| `mitochondrial_pseudobulk` | donor_technical.csv | 7 | 69 | 0.961 | 0.918 | 0.415 | 0.799 |
| `mitochondrial_supertype_share` | donor_technical.csv | 7 | 69 | 0.870 | 0.786 | 0.327 | 0.513 |
| `hemoglobin_mean` | donor_technical.csv | 7 | 69 | 0.284 | 0.301 | -0.037 | -0.254 |
| `hemoglobin_sd` | donor_technical.csv | 7 | 69 | 0.286 | 0.306 | -0.012 | -0.252 |
| `hemoglobin_pseudobulk` | donor_technical.csv | 7 | 69 | 0.271 | 0.291 | -0.045 | -0.275 |
| `hemoglobin_supertype_share` | donor_technical.csv | 3 | 61 | 0.366 | 0.431 | -0.178 | -0.203 |
| `fraction_mitochondrial_umis_mean` | donor_technical.csv | 7 | 69 | 0.967 | 0.934 | 0.387 | 0.790 |
| `fraction_mitochondrial_umis_sd` | donor_technical.csv | 7 | 69 | 0.950 | 0.903 | 0.447 | 0.781 |
| `number_of_umis_mean` | donor_technical.csv | 7 | 69 | 0.464 | 0.423 | 0.135 | 0.535 |
| `number_of_umis_sd` | donor_technical.csv | 7 | 69 | 0.522 | 0.454 | 0.137 | 0.532 |
| `genes_detected_mean` | donor_technical.csv | 7 | 69 | 0.348 | 0.286 | 0.151 | 0.500 |
| `genes_detected_sd` | donor_technical.csv | 7 | 69 | 0.553 | 0.485 | 0.104 | 0.413 |
| `doublet_score_mean` | donor_technical.csv | 7 | 69 | 0.754 | 0.750 | -0.051 | 0.092 |
| `doublet_score_sd` | donor_technical.csv | 7 | 69 | 0.903 | 0.918 | 0.014 | 0.175 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 7 | 69 | 0.012 | 0.020 | -0.275 | -0.417 |
| `gex_reads_mapped_confidently_to_intronic_regions_sd` | donor_technical.csv | 7 | 69 | 0.716 | 0.694 | 0.086 | 0.045 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 7 | 69 | 0.027 | 0.041 | -0.189 | -0.319 |
| `gex_fraction_of_transcriptomic_reads_in_cells_sd` | donor_technical.csv | 7 | 69 | 0.824 | 0.791 | 0.013 | 0.007 |
| `gex_median_umi_counts_per_cell_mean` | donor_technical.csv | 7 | 69 | 0.029 | 0.041 | -0.130 | -0.256 |
| `gex_median_umi_counts_per_cell_sd` | donor_technical.csv | 7 | 69 | 0.327 | 0.342 | -0.067 | 0.052 |
| `gex_sequencing_saturation_mean` | donor_technical.csv | 7 | 69 | 0.938 | 0.893 | 0.313 | 0.641 |
| `gex_sequencing_saturation_sd` | donor_technical.csv | 7 | 69 | 0.392 | 0.426 | -0.043 | -0.074 |
| `atac_tss_enrichment_score_mean` | donor_technical.csv | 7 | 69 | 0.017 | 0.010 | -0.125 | -0.633 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | donor_technical.csv | 7 | 69 | 0.014 | 0.036 | -0.169 | -0.642 |
| `atac_fraction_of_high_quality_fragments_in_cells_mean` | donor_technical.csv | 7 | 69 | 0.002 | 0.005 | -0.057 | -0.359 |
| `atac_percent_duplicates_mean` | donor_technical.csv | 7 | 69 | 0.406 | 0.321 | 0.191 | -0.111 |
| `atac_median_high_quality_fragments_per_cell_mean` | donor_technical.csv | 7 | 69 | 0.207 | 0.235 | -0.139 | -0.014 |
| `tmm_rna` | manifest.json (export, tmm) | 7 | 69 | 0.404 | 0.372 | 0.037 | -0.372 |
| `log10_lib_rna` | manifest.json (export, tmm) | 7 | 69 | 0.360 | 0.321 | 0.024 | 0.237 |
| `tmm_atac` | manifest.json (export, tmm) | 7 | 69 | 0.948 | 0.939 | 0.204 | 0.300 |
| `log10_lib_atac` | manifest.json (export, tmm) | 7 | 69 | 0.079 | 0.107 | -0.159 | -0.175 |

##### delta

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 76 | 2 | 5.10924 | 0.705 | no | 10.8174 | 8.13722 | 26.9672 |
| excluding_severely_affected | 69 | 4 | 2.4993 | 0.798 | no | 7.5077 | 5.55404 | 18.2186 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 76 | 2 | 12.5272 | 0.456 | no | 13.6061 | 10.1863 | 32.1219 |
| excluding_severely_affected | 69 | 4 | 22.861 | 0.277 | no | 18.4684 | 8.424 | 33.627 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 76 | 2 | 0.0242177 | 0.991 | no | 0.272507 | 0.148066 | 0.544 |
| excluding_severely_affected | 69 | 4 | 0.0853244 | 0.61 | no | 0.106695 | 0.0558344 | 0.201529 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 76 | 0 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=24; Male/stage 0=9; Male/stage 1=11; Male/stage 2=11 | cv | 2 | 0.6852 | 107.7 | `results/microglia/motco_analysis/cohort_30` |
| excluding_severely_affected | 69 | 7 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=19; Male/stage 0=9; Male/stage 1=11; Male/stage 2=9 | cv | 4 | 0.6649 | 163.8 | `motco_analysis/cohort_30` |

Excluded donors: `H20.33.028`, `H20.33.033`, `H20.33.037`, `H20.33.038`, `H20.33.046`, `H21.33.009`, `H21.33.034`.

##### Axis check

Predeclared rule (`DONOR_EXCLUSION_AXIS_RULE`): along the variant analysis's first latent component, |Spearman| with stage exceeds the largest |Spearman| with a technical variable of POSITIVE_CONTROL_TECHNICAL_SOURCES (technical table columns for the variant donors; the variant export's TMM factors and log10 library sizes); the control cell type's value is the expectation, the reference's descriptive. A comparison, not a test; no p-value.

| cell type | role | reading | component | donors | ρ stage | largest technical | source | ρ technical | technical compared | holds |
|---|---|---|---|---:|---:|---|---|---:|---:|:---:|
| microglia | control | expectation | `pls_0` | 69 | 0.8201 | `fraction_mitochondrial_umis_sd` | donor_technical.csv | 0.6891 | 39 | yes |

##### Reading

Predeclared (`DONOR_EXCLUSION_READING`): the flag check and the axis check are read together: if the flagged donors carry the quality tail (high AUC for the mitochondrial fraction and the primary pls_0) and the axis check holds in the control cell type, the quality gradient of the primary analyses sits in the donors SEA-AD excluded from expression tests, and excluding them is a donor-selection route to the same end as the latent adjustment; if the axis check fails although the flagged donors carry the tail, the gradient extends beyond them and donor selection does not resolve it; if the flagged donors do not carry the tail, the flag is not the gradient's carrier and the variant is an ordinary cohort sensitivity. In every case the primary analyses and the primary-basis positive control stay as they are; a rejection status that differs between the variant and the primary is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside angle; the non-rejection rules hold in the variant.

##### Findings

- CV-selected rank differs: 2 in the primary analysis, 4 in the excluding_severely_affected analysis; reported beside angle, not read as better or worse
- axis check holds in cohort_30 (expectation, role control): |Spearman| of pls_0 with stage 0.8201 vs 0.6891 with fraction_mitochondrial_umis_sd (donor_technical.csv)
- excluding_severely_affected export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### 9.2 microglia, cohort_50 (≥ 50 microglia nuclei in both omics, n = 63) — sensitivity cohort

#### Donor-exclusion sensitivity: microglia, cohort_50, excluding_severely_affected

The `primary` row is the primary analysis (CV-selected rank 4, 63 donors), read from its `result.json` and never rerun; it stays the primary result. The `excluding_severely_affected` row reruns the whole export and the primary analysis on the cohort donors minus those SEA-AD flagged as severely affected (Gabitto et al. 2024: QC-outlying high-ADNC donors excluded from the gene-expression tests): 5 cohort donor(s) excluded, the cohort thresholds unchanged and k recomputed from the remaining cells by the export's rule (k = 8 → 8). The variant is a sensitivity analysis, never a revised primary: the excluded donors carry the most severe pathology, so the exclusion truncates the late epoch along with the quality tail, and the flag was made from QC metrics, so the variant conditions on quality by donor selection. A rejection status that differs between the rows is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in both rows. 999 RRPP permutations with seed 20261019 in both rows; α = 0.05 (rejected when p ≤ α). MOTCO's operating characteristics were measured on neither cohort.

##### Flagged donors per cell

| cell | flagged | unflagged |
|---|---:|---:|
| Female/stage 0 | 0 | 8 |
| Female/stage 1 | 0 | 10 |
| Female/stage 2 | 4 | 14 |
| Male/stage 0 | 0 | 8 |
| Male/stage 1 | 0 | 10 |
| Male/stage 2 | 1 | 8 |
| all | 5 | 58 |

##### Flag check

Written before the variant export ran. The rank-biserial AUC is P(flagged > unflagged) + ½ P(=) over the cohort donors (`auc_all`) and within stage 2 (`auc_stage2`); 0.5 means no separation. Among the unflagged donors, each variable's Spearman correlation with stage and with the primary `pls_0`. Descriptive, no p-value; the variant runs whatever it says.

| variable | source | n flagged | n unflagged | AUC all | AUC stage 2 | ρ stage (unflagged) | ρ pls_0 (unflagged) |
|---|---|---:|---:|---:|---:|---:|---:|
| `pls_0` | latent_scores.csv (primary analysis) | 5 | 58 | 0.983 | 0.955 | 0.699 | n/a |
| `stage` | metadata.csv | 5 | 58 | 0.810 | n/a | n/a | 0.699 |
| `n_nuclei_rna` | donor_technical.csv | 5 | 58 | 0.436 | 0.455 | -0.147 | -0.002 |
| `n_libraries_rna` | donor_technical.csv | 5 | 58 | 0.522 | 0.541 | -0.015 | 0.092 |
| `ribosomal_mean` | donor_technical.csv | 5 | 58 | 0.269 | 0.273 | -0.012 | -0.157 |
| `ribosomal_sd` | donor_technical.csv | 5 | 58 | 0.303 | 0.300 | 0.048 | 0.016 |
| `ribosomal_pseudobulk` | donor_technical.csv | 5 | 58 | 0.403 | 0.409 | 0.002 | -0.060 |
| `ribosomal_supertype_share` | donor_technical.csv | 5 | 58 | 0.762 | 0.673 | 0.279 | 0.282 |
| `mitochondrial_mean` | donor_technical.csv | 5 | 58 | 0.955 | 0.909 | 0.334 | 0.730 |
| `mitochondrial_sd` | donor_technical.csv | 5 | 58 | 0.934 | 0.873 | 0.391 | 0.722 |
| `mitochondrial_pseudobulk` | donor_technical.csv | 5 | 58 | 0.945 | 0.882 | 0.355 | 0.733 |
| `mitochondrial_supertype_share` | donor_technical.csv | 5 | 58 | 0.848 | 0.773 | 0.257 | 0.425 |
| `hemoglobin_mean` | donor_technical.csv | 5 | 58 | 0.379 | 0.400 | 0.004 | -0.229 |
| `hemoglobin_sd` | donor_technical.csv | 5 | 58 | 0.386 | 0.400 | 0.033 | -0.218 |
| `hemoglobin_pseudobulk` | donor_technical.csv | 5 | 58 | 0.366 | 0.382 | 0.013 | -0.219 |
| `hemoglobin_supertype_share` | donor_technical.csv | 3 | 50 | 0.373 | 0.426 | -0.127 | -0.116 |
| `fraction_mitochondrial_umis_mean` | donor_technical.csv | 5 | 58 | 0.955 | 0.909 | 0.334 | 0.730 |
| `fraction_mitochondrial_umis_sd` | donor_technical.csv | 5 | 58 | 0.934 | 0.873 | 0.391 | 0.722 |
| `number_of_umis_mean` | donor_technical.csv | 5 | 58 | 0.583 | 0.500 | 0.156 | 0.513 |
| `number_of_umis_sd` | donor_technical.csv | 5 | 58 | 0.610 | 0.509 | 0.148 | 0.511 |
| `genes_detected_mean` | donor_technical.csv | 5 | 58 | 0.428 | 0.309 | 0.174 | 0.475 |
| `genes_detected_sd` | donor_technical.csv | 5 | 58 | 0.645 | 0.555 | 0.134 | 0.447 |
| `doublet_score_mean` | donor_technical.csv | 5 | 58 | 0.714 | 0.709 | 0.005 | 0.078 |
| `doublet_score_sd` | donor_technical.csv | 5 | 58 | 0.886 | 0.900 | -0.014 | 0.150 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 5 | 58 | 0.021 | 0.036 | -0.281 | -0.387 |
| `gex_reads_mapped_confidently_to_intronic_regions_sd` | donor_technical.csv | 5 | 58 | 0.652 | 0.636 | 0.046 | 0.079 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 5 | 58 | 0.038 | 0.055 | -0.288 | -0.412 |
| `gex_fraction_of_transcriptomic_reads_in_cells_sd` | donor_technical.csv | 5 | 58 | 0.855 | 0.818 | 0.025 | 0.061 |
| `gex_median_umi_counts_per_cell_mean` | donor_technical.csv | 5 | 58 | 0.048 | 0.073 | -0.200 | -0.323 |
| `gex_median_umi_counts_per_cell_sd` | donor_technical.csv | 5 | 58 | 0.383 | 0.445 | -0.127 | -0.011 |
| `gex_sequencing_saturation_mean` | donor_technical.csv | 5 | 58 | 0.903 | 0.827 | 0.329 | 0.669 |
| `gex_sequencing_saturation_sd` | donor_technical.csv | 5 | 58 | 0.426 | 0.486 | -0.078 | -0.105 |
| `atac_tss_enrichment_score_mean` | donor_technical.csv | 5 | 58 | 0.028 | 0.018 | -0.056 | -0.529 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | donor_technical.csv | 5 | 58 | 0.024 | 0.064 | -0.154 | -0.615 |
| `atac_fraction_of_high_quality_fragments_in_cells_mean` | donor_technical.csv | 5 | 58 | 0.003 | 0.009 | -0.104 | -0.390 |
| `atac_percent_duplicates_mean` | donor_technical.csv | 5 | 58 | 0.376 | 0.309 | 0.117 | -0.193 |
| `atac_median_high_quality_fragments_per_cell_mean` | donor_technical.csv | 5 | 58 | 0.238 | 0.309 | -0.227 | -0.043 |
| `tmm_rna` | manifest.json (export, tmm) | 5 | 58 | 0.321 | 0.291 | -0.008 | -0.418 |
| `log10_lib_rna` | manifest.json (export, tmm) | 5 | 58 | 0.483 | 0.464 | -0.093 | 0.171 |
| `tmm_atac` | manifest.json (export, tmm) | 5 | 58 | 0.983 | 0.964 | 0.266 | 0.325 |
| `log10_lib_atac` | manifest.json (export, tmm) | 5 | 58 | 0.086 | 0.136 | -0.193 | -0.160 |

##### delta

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 63 | 4 | 1.6686 | 0.834 | no | 6.37019 | 4.88809 | 15.3753 |
| excluding_severely_affected | 58 | 4 | 2.55257 | 0.779 | no | 7.40384 | 5.56414 | 18.4315 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 63 | 4 | 22.2604 | 0.291 | no | 18.0464 | 8.60842 | 34.1958 |
| excluding_severely_affected | 58 | 4 | 19.1427 | 0.516 | no | 20.4516 | 9.34034 | 37.3416 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 63 | 4 | 0.0797433 | 0.644 | no | 0.106788 | 0.0560953 | 0.21359 |
| excluding_severely_affected | 58 | 4 | 0.0880231 | 0.583 | no | 0.106181 | 0.0546538 | 0.204244 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 63 | 0 | 8 | Female/stage 0=8; Female/stage 1=10; Female/stage 2=18; Male/stage 0=8; Male/stage 1=10; Male/stage 2=9 | cv | 4 | 0.6565 | 96.5 | `results/microglia/motco_analysis/cohort_50` |
| excluding_severely_affected | 58 | 5 | 8 | Female/stage 0=8; Female/stage 1=10; Female/stage 2=14; Male/stage 0=8; Male/stage 1=10; Male/stage 2=8 | cv | 4 | 0.6450 | 156.7 | `motco_analysis/cohort_50` |

Excluded donors: `H20.33.028`, `H20.33.037`, `H20.33.038`, `H20.33.046`, `H21.33.009`.

##### Axis check

Predeclared rule (`DONOR_EXCLUSION_AXIS_RULE`): along the variant analysis's first latent component, |Spearman| with stage exceeds the largest |Spearman| with a technical variable of POSITIVE_CONTROL_TECHNICAL_SOURCES (technical table columns for the variant donors; the variant export's TMM factors and log10 library sizes); the control cell type's value is the expectation, the reference's descriptive. A comparison, not a test; no p-value.

| cell type | role | reading | component | donors | ρ stage | largest technical | source | ρ technical | technical compared | holds |
|---|---|---|---|---:|---:|---|---|---:|---:|:---:|
| microglia | control | expectation | `pls_0` | 58 | 0.8248 | `fraction_mitochondrial_umis_sd` | donor_technical.csv | 0.6746 | 39 | yes |

##### Reading

Predeclared (`DONOR_EXCLUSION_READING`): the flag check and the axis check are read together: if the flagged donors carry the quality tail (high AUC for the mitochondrial fraction and the primary pls_0) and the axis check holds in the control cell type, the quality gradient of the primary analyses sits in the donors SEA-AD excluded from expression tests, and excluding them is a donor-selection route to the same end as the latent adjustment; if the axis check fails although the flagged donors carry the tail, the gradient extends beyond them and donor selection does not resolve it; if the flagged donors do not carry the tail, the flag is not the gradient's carrier and the variant is an ordinary cohort sensitivity. In every case the primary analyses and the primary-basis positive control stay as they are; a rejection status that differs between the variant and the primary is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside angle; the non-rejection rules hold in the variant.

##### Findings

- axis check holds in cohort_50 (expectation, role control): |Spearman| of pls_0 with stage 0.8248 vs 0.6746 with fraction_mitochondrial_umis_sd (donor_technical.csv)
- excluding_severely_affected export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### 9.3 astrocyte, cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

#### Donor-exclusion sensitivity: astrocyte, cohort_30, excluding_severely_affected

The `primary` row is the primary analysis (CV-selected rank 5, 80 donors), read from its `result.json` and never rerun; it stays the primary result. The `excluding_severely_affected` row reruns the whole export and the primary analysis on the cohort donors minus those SEA-AD flagged as severely affected (Gabitto et al. 2024: QC-outlying high-ADNC donors excluded from the gene-expression tests): 10 cohort donor(s) excluded, the cohort thresholds unchanged and k recomputed from the remaining cells by the export's rule (k = 9 → 9). The variant is a sensitivity analysis, never a revised primary: the excluded donors carry the most severe pathology, so the exclusion truncates the late epoch along with the quality tail, and the flag was made from QC metrics, so the variant conditions on quality by donor selection. A rejection status that differs between the rows is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in both rows. 999 RRPP permutations with seed 20261019 in both rows; α = 0.05 (rejected when p ≤ α). MOTCO's operating characteristics were measured on neither cohort.

##### Flagged donors per cell

| cell | flagged | unflagged |
|---|---:|---:|
| Female/stage 0 | 0 | 11 |
| Female/stage 1 | 0 | 10 |
| Female/stage 2 | 8 | 20 |
| Male/stage 0 | 0 | 9 |
| Male/stage 1 | 0 | 10 |
| Male/stage 2 | 2 | 10 |
| all | 10 | 70 |

##### Flag check

Written before the variant export ran. The rank-biserial AUC is P(flagged > unflagged) + ½ P(=) over the cohort donors (`auc_all`) and within stage 2 (`auc_stage2`); 0.5 means no separation. Among the unflagged donors, each variable's Spearman correlation with stage and with the primary `pls_0`. Descriptive, no p-value; the variant runs whatever it says.

| variable | source | n flagged | n unflagged | AUC all | AUC stage 2 | ρ stage (unflagged) | ρ pls_0 (unflagged) |
|---|---|---:|---:|---:|---:|---:|---:|
| `pls_0` | latent_scores.csv (primary analysis) | 10 | 70 | 0.961 | 0.937 | 0.384 | n/a |
| `stage` | metadata.csv | 10 | 70 | 0.786 | n/a | n/a | 0.384 |
| `n_nuclei_rna` | donor_technical.csv | 10 | 70 | 0.566 | 0.588 | -0.013 | 0.300 |
| `n_libraries_rna` | donor_technical.csv | 10 | 70 | 0.416 | 0.438 | -0.004 | 0.211 |
| `ribosomal_mean` | donor_technical.csv | 10 | 70 | 0.076 | 0.113 | -0.272 | -0.791 |
| `ribosomal_sd` | donor_technical.csv | 10 | 70 | 0.083 | 0.070 | -0.159 | -0.607 |
| `ribosomal_pseudobulk` | donor_technical.csv | 10 | 70 | 0.081 | 0.123 | -0.276 | -0.753 |
| `ribosomal_supertype_share` | donor_technical.csv | 10 | 70 | 0.326 | 0.310 | 0.094 | 0.029 |
| `mitochondrial_mean` | donor_technical.csv | 10 | 70 | 0.964 | 0.953 | 0.181 | 0.647 |
| `mitochondrial_sd` | donor_technical.csv | 10 | 70 | 0.923 | 0.910 | 0.092 | 0.523 |
| `mitochondrial_pseudobulk` | donor_technical.csv | 10 | 70 | 0.966 | 0.963 | 0.160 | 0.617 |
| `mitochondrial_supertype_share` | donor_technical.csv | 10 | 70 | 0.683 | 0.687 | -0.072 | 0.094 |
| `hemoglobin_mean` | donor_technical.csv | 10 | 70 | 0.499 | 0.480 | 0.146 | -0.153 |
| `hemoglobin_sd` | donor_technical.csv | 10 | 70 | 0.486 | 0.440 | 0.179 | -0.085 |
| `hemoglobin_pseudobulk` | donor_technical.csv | 10 | 70 | 0.497 | 0.483 | 0.100 | -0.231 |
| `hemoglobin_supertype_share` | donor_technical.csv | 8 | 67 | 0.584 | 0.612 | -0.048 | -0.321 |
| `fraction_mitochondrial_umis_mean` | donor_technical.csv | 10 | 70 | 0.964 | 0.953 | 0.181 | 0.647 |
| `fraction_mitochondrial_umis_sd` | donor_technical.csv | 10 | 70 | 0.923 | 0.910 | 0.092 | 0.523 |
| `number_of_umis_mean` | donor_technical.csv | 10 | 70 | 0.417 | 0.427 | -0.009 | 0.245 |
| `number_of_umis_sd` | donor_technical.csv | 10 | 70 | 0.423 | 0.417 | -0.053 | 0.290 |
| `genes_detected_mean` | donor_technical.csv | 10 | 70 | 0.354 | 0.353 | -0.012 | 0.097 |
| `genes_detected_sd` | donor_technical.csv | 10 | 70 | 0.396 | 0.390 | -0.074 | 0.238 |
| `doublet_score_mean` | donor_technical.csv | 10 | 70 | 0.557 | 0.553 | -0.006 | -0.000 |
| `doublet_score_sd` | donor_technical.csv | 10 | 70 | 0.723 | 0.663 | 0.069 | 0.181 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 10 | 70 | 0.021 | 0.033 | -0.231 | -0.329 |
| `gex_reads_mapped_confidently_to_intronic_regions_sd` | donor_technical.csv | 10 | 70 | 0.681 | 0.673 | 0.032 | 0.224 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 10 | 70 | 0.046 | 0.067 | -0.211 | -0.449 |
| `gex_fraction_of_transcriptomic_reads_in_cells_sd` | donor_technical.csv | 10 | 70 | 0.663 | 0.643 | -0.012 | 0.150 |
| `gex_median_umi_counts_per_cell_mean` | donor_technical.csv | 10 | 70 | 0.034 | 0.047 | -0.112 | -0.233 |
| `gex_median_umi_counts_per_cell_sd` | donor_technical.csv | 10 | 70 | 0.247 | 0.273 | -0.070 | 0.131 |
| `gex_sequencing_saturation_mean` | donor_technical.csv | 10 | 70 | 0.923 | 0.883 | 0.281 | 0.553 |
| `gex_sequencing_saturation_sd` | donor_technical.csv | 10 | 70 | 0.369 | 0.390 | -0.047 | -0.008 |
| `atac_tss_enrichment_score_mean` | donor_technical.csv | 10 | 70 | 0.017 | 0.013 | -0.050 | -0.677 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | donor_technical.csv | 10 | 70 | 0.013 | 0.030 | -0.121 | -0.674 |
| `atac_fraction_of_high_quality_fragments_in_cells_mean` | donor_technical.csv | 10 | 70 | 0.003 | 0.007 | -0.002 | -0.350 |
| `atac_percent_duplicates_mean` | donor_technical.csv | 10 | 70 | 0.389 | 0.323 | 0.129 | -0.187 |
| `atac_median_high_quality_fragments_per_cell_mean` | donor_technical.csv | 10 | 70 | 0.209 | 0.223 | -0.121 | 0.163 |
| `tmm_rna` | manifest.json (export, tmm) | 10 | 70 | 0.153 | 0.150 | -0.076 | -0.696 |
| `log10_lib_rna` | manifest.json (export, tmm) | 10 | 70 | 0.497 | 0.507 | -0.018 | 0.386 |
| `tmm_atac` | manifest.json (export, tmm) | 10 | 70 | 0.967 | 0.973 | 0.110 | 0.280 |
| `log10_lib_atac` | manifest.json (export, tmm) | 10 | 70 | 0.406 | 0.417 | -0.002 | 0.318 |

##### delta

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 80 | 5 | 5.95773 | 0.501 | no | 7.38883 | 5.76462 | 18.6479 |
| excluding_severely_affected | 70 | 2 | 4.67327 | 0.636 | no | 7.86545 | 6.27135 | 19.8853 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 80 | 5 | 21.9941 | 0.449 | no | 22.0577 | 9.8721 | 40.3393 |
| excluding_severely_affected | 70 | 2 | 10.6202 | 0.545 | no | 13.7715 | 10.3493 | 33.3216 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 80 | 5 | 0.15645 | 0.331 | no | 0.132161 | 0.0684683 | 0.26132 |
| excluding_severely_affected | 70 | 2 | 0.248331 | 0.407 | no | 0.227417 | 0.120507 | 0.439706 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 80 | 0 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=28; Male/stage 0=9; Male/stage 1=10; Male/stage 2=12 | cv | 5 | 0.6356 | 113.8 | `results/motco_analysis/cohort_30` |
| excluding_severely_affected | 70 | 10 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=20; Male/stage 0=9; Male/stage 1=10; Male/stage 2=10 | cv | 2 | 0.5892 | 139.3 | `motco_analysis/cohort_30` |

Excluded donors: `H20.33.028`, `H20.33.037`, `H20.33.038`, `H20.33.046`, `H21.33.009`, `H21.33.020`, `H21.33.033`, `H21.33.034`, `H21.33.039`, `H21.33.045`.

##### Axis check

Predeclared rule (`DONOR_EXCLUSION_AXIS_RULE`): along the variant analysis's first latent component, |Spearman| with stage exceeds the largest |Spearman| with a technical variable of POSITIVE_CONTROL_TECHNICAL_SOURCES (technical table columns for the variant donors; the variant export's TMM factors and log10 library sizes); the control cell type's value is the expectation, the reference's descriptive. A comparison, not a test; no p-value.

| cell type | role | reading | component | donors | ρ stage | largest technical | source | ρ technical | technical compared | holds |
|---|---|---|---|---:|---:|---|---|---:|---:|:---:|
| astrocyte | reference | descriptive | `pls_0` | 70 | 0.5784 | `ribosomal_mean` | donor_technical.csv | -0.7286 | 39 | no |

##### Reading

Predeclared (`DONOR_EXCLUSION_READING`): the flag check and the axis check are read together: if the flagged donors carry the quality tail (high AUC for the mitochondrial fraction and the primary pls_0) and the axis check holds in the control cell type, the quality gradient of the primary analyses sits in the donors SEA-AD excluded from expression tests, and excluding them is a donor-selection route to the same end as the latent adjustment; if the axis check fails although the flagged donors carry the tail, the gradient extends beyond them and donor selection does not resolve it; if the flagged donors do not carry the tail, the flag is not the gradient's carrier and the variant is an ordinary cohort sensitivity. In every case the primary analyses and the primary-basis positive control stay as they are; a rejection status that differs between the variant and the primary is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside angle; the non-rejection rules hold in the variant.

##### Findings

- CV-selected rank differs: 5 in the primary analysis, 2 in the excluding_severely_affected analysis; reported beside angle, not read as better or worse
- axis check does not hold in cohort_30 (descriptive, role reference): |Spearman| of pls_0 with stage 0.5784 vs 0.7286 with ribosomal_mean (donor_technical.csv)
- excluding_severely_affected export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### 9.4 astrocyte, cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74) — sensitivity cohort

#### Donor-exclusion sensitivity: astrocyte, cohort_50, excluding_severely_affected

The `primary` row is the primary analysis (CV-selected rank 2, 74 donors), read from its `result.json` and never rerun; it stays the primary result. The `excluding_severely_affected` row reruns the whole export and the primary analysis on the cohort donors minus those SEA-AD flagged as severely affected (Gabitto et al. 2024: QC-outlying high-ADNC donors excluded from the gene-expression tests): 10 cohort donor(s) excluded, the cohort thresholds unchanged and k recomputed from the remaining cells by the export's rule (k = 9 → 9). The variant is a sensitivity analysis, never a revised primary: the excluded donors carry the most severe pathology, so the exclusion truncates the late epoch along with the quality tail, and the flag was made from QC metrics, so the variant conditions on quality by donor selection. A rejection status that differs between the rows is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in both rows. 999 RRPP permutations with seed 20261019 in both rows; α = 0.05 (rejected when p ≤ α). MOTCO's operating characteristics were measured on neither cohort.

##### Flagged donors per cell

| cell | flagged | unflagged |
|---|---:|---:|
| Female/stage 0 | 0 | 10 |
| Female/stage 1 | 0 | 9 |
| Female/stage 2 | 8 | 17 |
| Male/stage 0 | 0 | 9 |
| Male/stage 1 | 0 | 9 |
| Male/stage 2 | 2 | 10 |
| all | 10 | 64 |

##### Flag check

Written before the variant export ran. The rank-biserial AUC is P(flagged > unflagged) + ½ P(=) over the cohort donors (`auc_all`) and within stage 2 (`auc_stage2`); 0.5 means no separation. Among the unflagged donors, each variable's Spearman correlation with stage and with the primary `pls_0`. Descriptive, no p-value; the variant runs whatever it says.

| variable | source | n flagged | n unflagged | AUC all | AUC stage 2 | ρ stage (unflagged) | ρ pls_0 (unflagged) |
|---|---|---:|---:|---:|---:|---:|---:|
| `pls_0` | latent_scores.csv (primary analysis) | 10 | 64 | 0.959 | 0.926 | 0.392 | n/a |
| `stage` | metadata.csv | 10 | 64 | 0.789 | n/a | n/a | 0.392 |
| `n_nuclei_rna` | donor_technical.csv | 10 | 64 | 0.560 | 0.565 | 0.062 | 0.339 |
| `n_libraries_rna` | donor_technical.csv | 10 | 64 | 0.408 | 0.420 | 0.065 | 0.237 |
| `ribosomal_mean` | donor_technical.csv | 10 | 64 | 0.069 | 0.093 | -0.213 | -0.780 |
| `ribosomal_sd` | donor_technical.csv | 10 | 64 | 0.069 | 0.026 | -0.084 | -0.577 |
| `ribosomal_pseudobulk` | donor_technical.csv | 10 | 64 | 0.073 | 0.100 | -0.202 | -0.746 |
| `ribosomal_supertype_share` | donor_technical.csv | 10 | 64 | 0.336 | 0.311 | 0.123 | 0.037 |
| `mitochondrial_mean` | donor_technical.csv | 10 | 64 | 0.963 | 0.948 | 0.213 | 0.695 |
| `mitochondrial_sd` | donor_technical.csv | 10 | 64 | 0.920 | 0.904 | 0.120 | 0.567 |
| `mitochondrial_pseudobulk` | donor_technical.csv | 10 | 64 | 0.966 | 0.959 | 0.201 | 0.662 |
| `mitochondrial_supertype_share` | donor_technical.csv | 10 | 64 | 0.675 | 0.674 | -0.040 | 0.082 |
| `hemoglobin_mean` | donor_technical.csv | 10 | 64 | 0.502 | 0.467 | 0.149 | -0.080 |
| `hemoglobin_sd` | donor_technical.csv | 10 | 64 | 0.487 | 0.433 | 0.177 | -0.028 |
| `hemoglobin_pseudobulk` | donor_technical.csv | 10 | 64 | 0.505 | 0.474 | 0.105 | -0.162 |
| `hemoglobin_supertype_share` | donor_technical.csv | 8 | 61 | 0.590 | 0.625 | -0.092 | -0.309 |
| `fraction_mitochondrial_umis_mean` | donor_technical.csv | 10 | 64 | 0.963 | 0.948 | 0.213 | 0.695 |
| `fraction_mitochondrial_umis_sd` | donor_technical.csv | 10 | 64 | 0.920 | 0.904 | 0.120 | 0.567 |
| `number_of_umis_mean` | donor_technical.csv | 10 | 64 | 0.417 | 0.430 | -0.053 | 0.196 |
| `number_of_umis_sd` | donor_technical.csv | 10 | 64 | 0.427 | 0.426 | -0.079 | 0.246 |
| `genes_detected_mean` | donor_technical.csv | 10 | 64 | 0.350 | 0.348 | -0.057 | 0.052 |
| `genes_detected_sd` | donor_technical.csv | 10 | 64 | 0.400 | 0.400 | -0.076 | 0.196 |
| `doublet_score_mean` | donor_technical.csv | 10 | 64 | 0.552 | 0.552 | -0.031 | -0.033 |
| `doublet_score_sd` | donor_technical.csv | 10 | 64 | 0.714 | 0.656 | 0.054 | 0.192 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 10 | 64 | 0.019 | 0.026 | -0.263 | -0.344 |
| `gex_reads_mapped_confidently_to_intronic_regions_sd` | donor_technical.csv | 10 | 64 | 0.683 | 0.663 | 0.087 | 0.266 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 10 | 64 | 0.045 | 0.067 | -0.213 | -0.444 |
| `gex_fraction_of_transcriptomic_reads_in_cells_sd` | donor_technical.csv | 10 | 64 | 0.669 | 0.641 | 0.026 | 0.181 |
| `gex_median_umi_counts_per_cell_mean` | donor_technical.csv | 10 | 64 | 0.025 | 0.026 | -0.113 | -0.268 |
| `gex_median_umi_counts_per_cell_sd` | donor_technical.csv | 10 | 64 | 0.233 | 0.267 | -0.055 | 0.118 |
| `gex_sequencing_saturation_mean` | donor_technical.csv | 10 | 64 | 0.923 | 0.881 | 0.292 | 0.628 |
| `gex_sequencing_saturation_sd` | donor_technical.csv | 10 | 64 | 0.356 | 0.385 | -0.019 | -0.061 |
| `atac_tss_enrichment_score_mean` | donor_technical.csv | 10 | 64 | 0.016 | 0.015 | -0.067 | -0.678 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | donor_technical.csv | 10 | 64 | 0.014 | 0.033 | -0.152 | -0.695 |
| `atac_fraction_of_high_quality_fragments_in_cells_mean` | donor_technical.csv | 10 | 64 | 0.003 | 0.007 | -0.059 | -0.403 |
| `atac_percent_duplicates_mean` | donor_technical.csv | 10 | 64 | 0.386 | 0.330 | 0.101 | -0.263 |
| `atac_median_high_quality_fragments_per_cell_mean` | donor_technical.csv | 10 | 64 | 0.222 | 0.241 | -0.128 | 0.166 |
| `tmm_rna` | manifest.json (export, tmm) | 10 | 64 | 0.142 | 0.148 | -0.078 | -0.670 |
| `log10_lib_rna` | manifest.json (export, tmm) | 10 | 64 | 0.494 | 0.489 | 0.051 | 0.421 |
| `tmm_atac` | manifest.json (export, tmm) | 10 | 64 | 0.992 | 0.989 | 0.145 | 0.342 |
| `log10_lib_atac` | manifest.json (export, tmm) | 10 | 64 | 0.361 | 0.370 | 0.019 | 0.402 |

##### delta

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 74 | 2 | 6.01569 | 0.628 | no | 10.5899 | 8.56419 | 27.3174 |
| excluding_severely_affected | 64 | 2 | 3.27569 | 0.753 | no | 8.5507 | 7.0053 | 22.2778 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 74 | 2 | 23.1596 | 0.242 | no | 16.3255 | 13.4558 | 43.0401 |
| excluding_severely_affected | 64 | 2 | 9.0544 | 0.611 | no | 13.3587 | 9.43745 | 31.1722 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 74 | 2 | 0.284932 | 0.514 | no | 0.313708 | 0.17022 | 0.630147 |
| excluding_severely_affected | 64 | 2 | 0.255222 | 0.384 | no | 0.229393 | 0.126116 | 0.468445 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 74 | 0 | 9 | Female/stage 0=10; Female/stage 1=9; Female/stage 2=25; Male/stage 0=9; Male/stage 1=9; Male/stage 2=12 | cv | 2 | 0.6965 | 102.3 | `results/motco_analysis/cohort_50` |
| excluding_severely_affected | 64 | 10 | 9 | Female/stage 0=10; Female/stage 1=9; Female/stage 2=17; Male/stage 0=9; Male/stage 1=9; Male/stage 2=10 | cv | 2 | 0.6012 | 129.0 | `motco_analysis/cohort_50` |

Excluded donors: `H20.33.028`, `H20.33.037`, `H20.33.038`, `H20.33.046`, `H21.33.009`, `H21.33.020`, `H21.33.033`, `H21.33.034`, `H21.33.039`, `H21.33.045`.

##### Axis check

Predeclared rule (`DONOR_EXCLUSION_AXIS_RULE`): along the variant analysis's first latent component, |Spearman| with stage exceeds the largest |Spearman| with a technical variable of POSITIVE_CONTROL_TECHNICAL_SOURCES (technical table columns for the variant donors; the variant export's TMM factors and log10 library sizes); the control cell type's value is the expectation, the reference's descriptive. A comparison, not a test; no p-value.

| cell type | role | reading | component | donors | ρ stage | largest technical | source | ρ technical | technical compared | holds |
|---|---|---|---|---:|---:|---|---|---:|---:|:---:|
| astrocyte | reference | descriptive | `pls_0` | 64 | 0.5955 | `ribosomal_mean` | donor_technical.csv | -0.6716 | 39 | no |

##### Reading

Predeclared (`DONOR_EXCLUSION_READING`): the flag check and the axis check are read together: if the flagged donors carry the quality tail (high AUC for the mitochondrial fraction and the primary pls_0) and the axis check holds in the control cell type, the quality gradient of the primary analyses sits in the donors SEA-AD excluded from expression tests, and excluding them is a donor-selection route to the same end as the latent adjustment; if the axis check fails although the flagged donors carry the tail, the gradient extends beyond them and donor selection does not resolve it; if the flagged donors do not carry the tail, the flag is not the gradient's carrier and the variant is an ordinary cohort sensitivity. In every case the primary analyses and the primary-basis positive control stay as they are; a rejection status that differs between the variant and the primary is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside angle; the non-rejection rules hold in the variant.

##### Findings

- axis check does not hold in cohort_50 (descriptive, role reference): |Spearman| of pls_0 with stage 0.5955 vs 0.6716 with ribosomal_mean (donor_technical.csv)
- excluding_severely_affected export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### 9.5 Superseded → corrected (exclusion variants)

| cell type | cohort | analysis | rank | CV AUROC | κ (n, reduced rank) | statistic | observed | p superseded → corrected | change | rejected superseded → corrected |
|---|---|---|---:|---:|---|---|---:|---|---|---|
| astrocyte | cohort_30 | excluding_severely_affected | 2 | 0.5892 | 1.0801 (70, 10) | delta | 4.67327 | 0.599 → 0.636 | up | no → no |
| astrocyte | cohort_30 | excluding_severely_affected |  |  |  | angle | 10.6202 | 0.503 → 0.545 | up | no → no |
| astrocyte | cohort_30 | excluding_severely_affected |  |  |  | shape | 0.248331 | 0.356 → 0.407 | up | no → no |
| astrocyte | cohort_50 | excluding_severely_affected | 2 | 0.6012 | 1.0887 (64, 10) | delta | 3.27569 | 0.734 → 0.753 | up | no → no |
| astrocyte | cohort_50 | excluding_severely_affected |  |  |  | angle | 9.0544 | 0.576 → 0.611 | up | no → no |
| astrocyte | cohort_50 | excluding_severely_affected |  |  |  | shape | 0.255222 | 0.336 → 0.384 | up | no → no |
| microglia | cohort_30 | excluding_severely_affected | 4 | 0.6649 | 1.0814 (69, 10) | delta | 2.4993 | 0.781 → 0.798 | up | no → no |
| microglia | cohort_30 | excluding_severely_affected |  |  |  | angle | 22.861 | 0.218 → 0.277 | up | no → no |
| microglia | cohort_30 | excluding_severely_affected |  |  |  | shape | 0.0853244 | 0.569 → 0.61 | up | no → no |
| microglia | cohort_50 | excluding_severely_affected | 4 | 0.6450 | 1.0992 (58, 10) | delta | 2.55257 | 0.758 → 0.779 | up | no → no |
| microglia | cohort_50 | excluding_severely_affected |  |  |  | angle | 19.1427 | 0.45 → 0.516 | up | no → no |
| microglia | cohort_50 | excluding_severely_affected |  |  |  | shape | 0.0880231 | 0.52 → 0.583 | up | no → no |

Superseded against corrected, generated by the comparison script over 4 analyses: of the 16 analysis CSVs (four per analysis), 16 equal at rounding level (largest |difference| 8.7e-13); of the 4 records outside the null keys, 4 equal at rounding level (largest relative difference 2.9e-12, apart from the numerically-zero third eigenvalues of the configuration spectra, ≤ 4.6e-33); rank equal in 4 of 4, CV AUROC equal in 4 of 4; κ = √(n / (n − r)) to 1e-12 in 4 of 4; p-values up in 12, equal in 0, down in 0 of 12; rejection statuses changed: 0.

## 10. Positive control (microglia against astrocytes)

No expectation reads a p-value, so E1–E3 were expected to reproduce on both bases, and they do: the four `positive_control.csv` are byte-identical to the superseded ones. On the corrected primary the control's outcome in `cohort_30` is unchanged: E1 holds (DESeq2 61 against 1; limma 2 against 1), E2 holds (0.6852 against 0.6356), E3 fails (0.6882 against `mitochondrial_pseudobulk` 0.8438), so the control **does not hold**, through E3 alone; `cohort_50` stays descriptive (E1 `method_dependent`, E2 0.6565 against 0.6965 and E3 0.7308 against 0.7846 fail). The restated control on `latent_adjusted` keeps its predeclared reading and its outcome: in `cohort_30` it holds (E1 holds, E2 0.6575 against 0.6433, E3 0.5279 against `tmm_rna` 0.4153); `cohort_50` stays descriptive (E1 `method_dependent`, E2 0.6261 against 0.6645 fails, E3 0.5327 against 0.4067 holds). E2 and E3 on that basis read the arm's CV AUROC and first-component correlations only, which MOTCO's `residualized_basis` diagnosis (a size excess of the arm's `angle` and `shape` tests) does not bear on.

### 10.1 basis `primary`, cohort_30

Pasted from `results/microglia/positive_control/cohort_30/positive_control.md` (inputs with sha256 in the file):

#### Positive control: microglia against astrocyte (cohort_30)

Predeclared expectations (POSITIVE_CONTROL_EXPECTATIONS), evaluated from the generated files of both cell types. The control holds when E1, E2 and E3 all hold in the primary cohort (POSITIVE_CONTROL_PRIMARY_COHORT); if it holds, the pipeline resolves stage structure on real data and the reference's null is not a pipeline failure; if not, the latent construction is in question for any cell type. No expectation reads a p-value, a rejection status or a trajectory statistic; the values are compared, not tested.

This is a positive control for the stage signal, not for the sex × stage hypothesis: the microglia trajectory result is read under the same binding rules as the astrocyte one, and nothing here changes any astrocyte output or reading.

**Primary cohort.** Overall: **does not hold** (E3 not holding).

| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |
|---|---|---|---|---|---|---|---|---|---|
| E1 | deseq2 | exported RNA features at FDR 0.05: factor_additive union | microglia | 61 | astrocyte | 1 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | limma | exported RNA features at FDR 0.05: factor_additive stage_joint | microglia | 2 | astrocyte | 1 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | both | E1 under both stage-DE methods |  |  |  |  | yes | holds | holds under both methods |
| E2 |  | CV mean AUROC of the primary analysis (latent.cv_mean_auroc) | microglia | 0.685213 | astrocyte | 0.635556 | yes | holds |  |
| E3 |  | \|Spearman\| of pls_0 (primary) with stage vs the largest with a technical variable | stage | 0.688175 | mitochondrial_pseudobulk | 0.843828 | no | fails | 39 technical variables compared (donor_technical.csv, manifest.json (export, tmm)) |
| overall |  | E1, E2 and E3 all hold |  |  |  |  | no | does_not_hold | the primary cohort: this row is the control's outcome; not holding: E3 |

##### Expectations

- **E1 (stage signal)**: under each stage-DE method, the number of exported RNA features with a stage effect at FDR 0.05 is larger for the control than for the reference, compared within the method. Holds only when it holds under both methods; when it holds under one only it is method_dependent, which counts as not holding, and the method under which it held is named.
- **E2 (latent separation)**: the control's primary analysis has a larger CV mean AUROC than the reference's.
- **E3 (axis identity)**: along the first latent component of the control's primary analysis, |Spearman| with stage exceeds |Spearman| with every technical variable of the diagnostic's variable set.

##### Inputs

- control (microglia): `results/microglia`
  - `motco_analysis/cohort_30/result.json` sha256 `953a1c453e16…`
  - `motco_analysis/cohort_30/analysis_manifest.json` sha256 `827b16619b1c…`
  - `stage_de/cohort_30/stage_de_summary.csv` sha256 `d688b767a7bc…`
  - `stage_de/cohort_30/stage_de_manifest.json` sha256 `84207c2cbc0b…`
  - `stage_de_limma/cohort_30/stage_de_limma_summary.csv` sha256 `cca81af3a00d…`
  - `stage_de_limma/cohort_30/stage_de_limma_manifest.json` sha256 `2a28d21edfcf…`
  - `motco_diagnostics/cohort_30/axis_correlations.csv` sha256 `80e59859092c…`
  - `motco_diagnostics/cohort_30/diagnostics_manifest.json` sha256 `72400b6c7b43…`
- reference (astrocyte): `results`
  - `motco_analysis/cohort_30/result.json` sha256 `cbeeba6e1666…`
  - `motco_analysis/cohort_30/analysis_manifest.json` sha256 `2cf659eb1b66…`
  - `stage_de/cohort_30/stage_de_summary.csv` sha256 `6b4947e8db81…`
  - `stage_de/cohort_30/stage_de_manifest.json` sha256 `84d43e6c031d…`
  - `stage_de_limma/cohort_30/stage_de_limma_summary.csv` sha256 `f3a91146cafc…`
  - `stage_de_limma/cohort_30/stage_de_limma_manifest.json` sha256 `810d4fb89a56…`

### 10.2 basis `primary`, cohort_50

Pasted from `results/microglia/positive_control/cohort_50/positive_control.md` (inputs with sha256 in the file):

#### Positive control: microglia against astrocyte (cohort_50)

Predeclared expectations (POSITIVE_CONTROL_EXPECTATIONS), evaluated from the generated files of both cell types. The control holds when E1, E2 and E3 all hold in the primary cohort (POSITIVE_CONTROL_PRIMARY_COHORT); if it holds, the pipeline resolves stage structure on real data and the reference's null is not a pipeline failure; if not, the latent construction is in question for any cell type. No expectation reads a p-value, a rejection status or a trajectory statistic; the values are compared, not tested.

This is a positive control for the stage signal, not for the sex × stage hypothesis: the microglia trajectory result is read under the same binding rules as the astrocyte one, and nothing here changes any astrocyte output or reading.

**Not the primary cohort** (cohort_30 is): descriptive. Overall: **does not hold** (E1, E2, E3 not holding).

| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |
|---|---|---|---|---|---|---|---|---|---|
| E1 | deseq2 | exported RNA features at FDR 0.05: factor_additive union | microglia | 35 | astrocyte | 2 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | limma | exported RNA features at FDR 0.05: factor_additive stage_joint | microglia | 0 | astrocyte | 1 | no | fails | tested 5000 (control), 5000 (reference) |
| E1 | both | E1 under both stage-DE methods |  |  |  |  | no | method_dependent | held under deseq2 only; counts as not holding |
| E2 |  | CV mean AUROC of the primary analysis (latent.cv_mean_auroc) | microglia | 0.656488 | astrocyte | 0.696504 | no | fails |  |
| E3 |  | \|Spearman\| of pls_0 (primary) with stage vs the largest with a technical variable | stage | 0.730773 | mitochondrial_pseudobulk | 0.784562 | no | fails | 39 technical variables compared (donor_technical.csv, manifest.json (export, tmm)) |
| overall |  | E1, E2 and E3 all hold |  |  |  |  | no | does_not_hold | not the primary cohort (cohort_30); descriptive; not holding: E1, E2, E3 |

##### Expectations

- **E1 (stage signal)**: under each stage-DE method, the number of exported RNA features with a stage effect at FDR 0.05 is larger for the control than for the reference, compared within the method. Holds only when it holds under both methods; when it holds under one only it is method_dependent, which counts as not holding, and the method under which it held is named.
- **E2 (latent separation)**: the control's primary analysis has a larger CV mean AUROC than the reference's.
- **E3 (axis identity)**: along the first latent component of the control's primary analysis, |Spearman| with stage exceeds |Spearman| with every technical variable of the diagnostic's variable set.

##### Inputs

- control (microglia): `results/microglia`
  - `motco_analysis/cohort_50/result.json` sha256 `5adb341b3c5e…`
  - `motco_analysis/cohort_50/analysis_manifest.json` sha256 `ffd44822b8fa…`
  - `stage_de/cohort_50/stage_de_summary.csv` sha256 `c21538e81570…`
  - `stage_de/cohort_50/stage_de_manifest.json` sha256 `37bb4b3848fc…`
  - `stage_de_limma/cohort_50/stage_de_limma_summary.csv` sha256 `3966f52e7798…`
  - `stage_de_limma/cohort_50/stage_de_limma_manifest.json` sha256 `16d833a4b3dd…`
  - `motco_diagnostics/cohort_50/axis_correlations.csv` sha256 `21277a529547…`
  - `motco_diagnostics/cohort_50/diagnostics_manifest.json` sha256 `e9e68afaab27…`
- reference (astrocyte): `results`
  - `motco_analysis/cohort_50/result.json` sha256 `3f72fa6c2a59…`
  - `motco_analysis/cohort_50/analysis_manifest.json` sha256 `95075d4aa6be…`
  - `stage_de/cohort_50/stage_de_summary.csv` sha256 `9f33acf5b59e…`
  - `stage_de/cohort_50/stage_de_manifest.json` sha256 `09d2b15bd6dc…`
  - `stage_de_limma/cohort_50/stage_de_limma_summary.csv` sha256 `37e1bb49243c…`
  - `stage_de_limma/cohort_50/stage_de_limma_manifest.json` sha256 `bb960fb52053…`

### 10.3 basis `latent_adjusted`, cohort_30

Pasted from `results/microglia/positive_control/latent_adjusted/cohort_30/positive_control.md` (inputs with sha256 in the file):

#### Positive control: microglia against astrocyte (cohort_30), basis `latent_adjusted`

Predeclared expectations (POSITIVE_CONTROL_EXPECTATIONS), evaluated from the generated files of both cell types. The control holds when E1, E2 and E3 all hold in the primary cohort (POSITIVE_CONTROL_PRIMARY_COHORT); if it holds, the pipeline resolves stage structure on real data and the reference's null is not a pipeline failure; if not, the latent construction is in question for any cell type. No expectation reads a p-value, a rejection status or a trajectory statistic; the values are compared, not tested.

This is a positive control for the stage signal, not for the sex × stage hypothesis: the microglia trajectory result is read under the same binding rules as the astrocyte one, and nothing here changes any astrocyte output or reading.

**Basis `latent_adjusted`** (POSITIVE_CONTROL_BASES): E2 and E3 read the latent_adjusted arm of latent-sensitivity of each cell type and the diagnostic made from it; E1 is basis-independent and equals the primary-basis rows. The expectations, their comparisons, E3's technical variables and the all-three rule are the primary basis's. Read as a diagnosis of the primary-basis control, which stays the control's outcome: if the restated control holds in the primary cohort, the measured quality gradient explains the primary control's failure (once the library-quality direction is conditioned out of the basis, the pipeline resolves stage in the cell type with the stronger stage signal), without promoting the arm (MOTCO's power comparison, candidate_fails, excludes the latent adjustment as a replacement of the primary), changing the primary-basis outcome or changing the reference cell type's reading; if it does not hold, the latent construction is in question beyond the quality gradient (the composition of the control population, the rank-1 selection and the control's ATAC nuclei per donor are the stated candidates) and the diagnosis 'technical gradient' is not supported. E3 on this basis carries the mediator caveat in both directions: the adjustment removes the quality direction and, with it, any stage signal the added fractions carry, so the stage |rho| may rise or fall; a failure for that reason is still a failure and the expectation is not re-tuned. The sensitivity cohort is descriptive, as on the primary basis.

**Primary cohort.** Overall: **holds**.

| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |
|---|---|---|---|---|---|---|---|---|---|
| E1 | deseq2 | exported RNA features at FDR 0.05: factor_additive union | microglia | 61 | astrocyte | 1 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | limma | exported RNA features at FDR 0.05: factor_additive stage_joint | microglia | 2 | astrocyte | 1 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | both | E1 under both stage-DE methods |  |  |  |  | yes | holds | holds under both methods |
| E2 |  | CV mean AUROC of the latent_adjusted arm analysis (basis latent_adjusted) (latent.cv_mean_auroc) | microglia | 0.65748 | astrocyte | 0.643333 | yes | holds |  |
| E3 |  | \|Spearman\| of pls_0 (latent_adjusted arm, labelled primary by its diagnostic) with stage vs the largest with a technical variable | stage | 0.527877 | tmm_rna | 0.415256 | yes | holds | 39 technical variables compared (donor_technical.csv, manifest.json (export, tmm)) |
| overall |  | E1, E2 and E3 all hold |  |  |  |  | yes | holds | the primary cohort: this row is the control's outcome |

##### Expectations

- **E1 (stage signal)**: under each stage-DE method, the number of exported RNA features with a stage effect at FDR 0.05 is larger for the control than for the reference, compared within the method. Holds only when it holds under both methods; when it holds under one only it is method_dependent, which counts as not holding, and the method under which it held is named.
- **E2 (latent separation)**: the control's primary analysis has a larger CV mean AUROC than the reference's.
- **E3 (axis identity)**: along the first latent component of the control's primary analysis, |Spearman| with stage exceeds |Spearman| with every technical variable of the diagnostic's variable set.

##### Inputs

- control (microglia): `results/microglia`
  - `motco_analysis/cohort_30/result.json` sha256 `953a1c453e16…`
  - `motco_analysis/cohort_30/analysis_manifest.json` sha256 `827b16619b1c…`
  - `motco_latent_sensitivity/latent_adjusted/motco_analysis/cohort_30/result.json` sha256 `92a6e44e3c39…`
  - `motco_latent_sensitivity/latent_adjusted/motco_analysis/cohort_30/analysis_manifest.json` sha256 `1a7374f41565…`
  - `motco_input/cohort_30/rna.csv` sha256 `66a1f0d6abcc…`
  - `motco_input/cohort_30/atac.csv` sha256 `215f0bee3009…`
  - `motco_input/cohort_30/metadata.csv` sha256 `1a52aa11f13a…`
  - `motco_input/cohort_30/covariates.csv` sha256 `53b50e34ae5c…`
  - `motco_input/cohort_30/manifest.json` sha256 `fb6afd9c24d8…`
  - `motco_latent_sensitivity/comparison/cohort_30/latent_sensitivity_manifest.json` sha256 `dd18fc0944f7…`
  - `stage_de/cohort_30/stage_de_summary.csv` sha256 `d688b767a7bc…`
  - `stage_de/cohort_30/stage_de_manifest.json` sha256 `84207c2cbc0b…`
  - `stage_de_limma/cohort_30/stage_de_limma_summary.csv` sha256 `cca81af3a00d…`
  - `stage_de_limma/cohort_30/stage_de_limma_manifest.json` sha256 `2a28d21edfcf…`
  - `motco_latent_sensitivity/diagnostics/latent_adjusted/cohort_30/axis_correlations.csv` sha256 `53f980afd0f5…`
  - `motco_latent_sensitivity/diagnostics/latent_adjusted/cohort_30/diagnostics_manifest.json` sha256 `7812bbeff4ea…`
- reference (astrocyte): `results`
  - `motco_analysis/cohort_30/result.json` sha256 `cbeeba6e1666…`
  - `motco_analysis/cohort_30/analysis_manifest.json` sha256 `2cf659eb1b66…`
  - `motco_latent_sensitivity/latent_adjusted/motco_analysis/cohort_30/result.json` sha256 `a7466ac719bf…`
  - `motco_latent_sensitivity/latent_adjusted/motco_analysis/cohort_30/analysis_manifest.json` sha256 `4d172e0646a6…`
  - `motco_input/cohort_30/rna.csv` sha256 `186198b9c52d…`
  - `motco_input/cohort_30/atac.csv` sha256 `32ebd055a8f5…`
  - `motco_input/cohort_30/metadata.csv` sha256 `e2d8a4127eee…`
  - `motco_input/cohort_30/covariates.csv` sha256 `7ed4b1bffd19…`
  - `motco_input/cohort_30/manifest.json` sha256 `18f2289e0d4c…`
  - `motco_latent_sensitivity/comparison/cohort_30/latent_sensitivity_manifest.json` sha256 `2bc1a6430a69…`
  - `stage_de/cohort_30/stage_de_summary.csv` sha256 `6b4947e8db81…`
  - `stage_de/cohort_30/stage_de_manifest.json` sha256 `84d43e6c031d…`
  - `stage_de_limma/cohort_30/stage_de_limma_summary.csv` sha256 `f3a91146cafc…`
  - `stage_de_limma/cohort_30/stage_de_limma_manifest.json` sha256 `810d4fb89a56…`

### 10.4 basis `latent_adjusted`, cohort_50

Pasted from `results/microglia/positive_control/latent_adjusted/cohort_50/positive_control.md` (inputs with sha256 in the file):

#### Positive control: microglia against astrocyte (cohort_50), basis `latent_adjusted`

Predeclared expectations (POSITIVE_CONTROL_EXPECTATIONS), evaluated from the generated files of both cell types. The control holds when E1, E2 and E3 all hold in the primary cohort (POSITIVE_CONTROL_PRIMARY_COHORT); if it holds, the pipeline resolves stage structure on real data and the reference's null is not a pipeline failure; if not, the latent construction is in question for any cell type. No expectation reads a p-value, a rejection status or a trajectory statistic; the values are compared, not tested.

This is a positive control for the stage signal, not for the sex × stage hypothesis: the microglia trajectory result is read under the same binding rules as the astrocyte one, and nothing here changes any astrocyte output or reading.

**Basis `latent_adjusted`** (POSITIVE_CONTROL_BASES): E2 and E3 read the latent_adjusted arm of latent-sensitivity of each cell type and the diagnostic made from it; E1 is basis-independent and equals the primary-basis rows. The expectations, their comparisons, E3's technical variables and the all-three rule are the primary basis's. Read as a diagnosis of the primary-basis control, which stays the control's outcome: if the restated control holds in the primary cohort, the measured quality gradient explains the primary control's failure (once the library-quality direction is conditioned out of the basis, the pipeline resolves stage in the cell type with the stronger stage signal), without promoting the arm (MOTCO's power comparison, candidate_fails, excludes the latent adjustment as a replacement of the primary), changing the primary-basis outcome or changing the reference cell type's reading; if it does not hold, the latent construction is in question beyond the quality gradient (the composition of the control population, the rank-1 selection and the control's ATAC nuclei per donor are the stated candidates) and the diagnosis 'technical gradient' is not supported. E3 on this basis carries the mediator caveat in both directions: the adjustment removes the quality direction and, with it, any stage signal the added fractions carry, so the stage |rho| may rise or fall; a failure for that reason is still a failure and the expectation is not re-tuned. The sensitivity cohort is descriptive, as on the primary basis.

**Not the primary cohort** (cohort_30 is): descriptive. Overall: **does not hold** (E1, E2 not holding).

| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |
|---|---|---|---|---|---|---|---|---|---|
| E1 | deseq2 | exported RNA features at FDR 0.05: factor_additive union | microglia | 35 | astrocyte | 2 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | limma | exported RNA features at FDR 0.05: factor_additive stage_joint | microglia | 0 | astrocyte | 1 | no | fails | tested 5000 (control), 5000 (reference) |
| E1 | both | E1 under both stage-DE methods |  |  |  |  | no | method_dependent | held under deseq2 only; counts as not holding |
| E2 |  | CV mean AUROC of the latent_adjusted arm analysis (basis latent_adjusted) (latent.cv_mean_auroc) | microglia | 0.626085 | astrocyte | 0.664479 | no | fails |  |
| E3 |  | \|Spearman\| of pls_0 (latent_adjusted arm, labelled primary by its diagnostic) with stage vs the largest with a technical variable | stage | 0.532681 | gex_reads_mapped_confidently_to_intronic_regions_sd | 0.406702 | yes | holds | 39 technical variables compared (donor_technical.csv, manifest.json (export, tmm)) |
| overall |  | E1, E2 and E3 all hold |  |  |  |  | no | does_not_hold | not the primary cohort (cohort_30); descriptive; not holding: E1, E2 |

##### Expectations

- **E1 (stage signal)**: under each stage-DE method, the number of exported RNA features with a stage effect at FDR 0.05 is larger for the control than for the reference, compared within the method. Holds only when it holds under both methods; when it holds under one only it is method_dependent, which counts as not holding, and the method under which it held is named.
- **E2 (latent separation)**: the control's primary analysis has a larger CV mean AUROC than the reference's.
- **E3 (axis identity)**: along the first latent component of the control's primary analysis, |Spearman| with stage exceeds |Spearman| with every technical variable of the diagnostic's variable set.

##### Inputs

- control (microglia): `results/microglia`
  - `motco_analysis/cohort_50/result.json` sha256 `5adb341b3c5e…`
  - `motco_analysis/cohort_50/analysis_manifest.json` sha256 `ffd44822b8fa…`
  - `motco_latent_sensitivity/latent_adjusted/motco_analysis/cohort_50/result.json` sha256 `9c73ccb23478…`
  - `motco_latent_sensitivity/latent_adjusted/motco_analysis/cohort_50/analysis_manifest.json` sha256 `06238c238722…`
  - `motco_input/cohort_50/rna.csv` sha256 `5b8d6cc3f83d…`
  - `motco_input/cohort_50/atac.csv` sha256 `4d813e5b3aea…`
  - `motco_input/cohort_50/metadata.csv` sha256 `1383ea3cb132…`
  - `motco_input/cohort_50/covariates.csv` sha256 `99c171578215…`
  - `motco_input/cohort_50/manifest.json` sha256 `8f0fd84ab291…`
  - `motco_latent_sensitivity/comparison/cohort_50/latent_sensitivity_manifest.json` sha256 `d58a518d5744…`
  - `stage_de/cohort_50/stage_de_summary.csv` sha256 `c21538e81570…`
  - `stage_de/cohort_50/stage_de_manifest.json` sha256 `37bb4b3848fc…`
  - `stage_de_limma/cohort_50/stage_de_limma_summary.csv` sha256 `3966f52e7798…`
  - `stage_de_limma/cohort_50/stage_de_limma_manifest.json` sha256 `16d833a4b3dd…`
  - `motco_latent_sensitivity/diagnostics/latent_adjusted/cohort_50/axis_correlations.csv` sha256 `f7509dde6422…`
  - `motco_latent_sensitivity/diagnostics/latent_adjusted/cohort_50/diagnostics_manifest.json` sha256 `f14ad111b576…`
- reference (astrocyte): `results`
  - `motco_analysis/cohort_50/result.json` sha256 `3f72fa6c2a59…`
  - `motco_analysis/cohort_50/analysis_manifest.json` sha256 `95075d4aa6be…`
  - `motco_latent_sensitivity/latent_adjusted/motco_analysis/cohort_50/result.json` sha256 `a7076726f035…`
  - `motco_latent_sensitivity/latent_adjusted/motco_analysis/cohort_50/analysis_manifest.json` sha256 `efecf60dac27…`
  - `motco_input/cohort_50/rna.csv` sha256 `270afbcdaca1…`
  - `motco_input/cohort_50/atac.csv` sha256 `a41ab8a3d2ab…`
  - `motco_input/cohort_50/metadata.csv` sha256 `fcb6639c2114…`
  - `motco_input/cohort_50/covariates.csv` sha256 `c1ba106dceec…`
  - `motco_input/cohort_50/manifest.json` sha256 `dd14a6e7d558…`
  - `motco_latent_sensitivity/comparison/cohort_50/latent_sensitivity_manifest.json` sha256 `1b9ff216c0ce…`
  - `stage_de/cohort_50/stage_de_summary.csv` sha256 `9f33acf5b59e…`
  - `stage_de/cohort_50/stage_de_manifest.json` sha256 `09d2b15bd6dc…`
  - `stage_de_limma/cohort_50/stage_de_limma_summary.csv` sha256 `37e1bb49243c…`
  - `stage_de_limma/cohort_50/stage_de_limma_manifest.json` sha256 `bb960fb52053…`

### 10.5 Superseded → corrected (E1–E3 values and outcomes)

| basis | cohort | expectation | method | control value superseded → corrected | compared value superseded → corrected | compared with | holds | outcome |
|---|---|---|---|---|---|---|---|---|
| primary | cohort_30 | E1 | deseq2 | 61 → 61 | 1 → 1 | astrocyte | True → True | holds → holds |
| primary | cohort_30 | E1 | limma | 2 → 2 | 1 → 1 | astrocyte | True → True | holds → holds |
| primary | cohort_30 | E1 | both |  |  |  | True → True | holds → holds |
| primary | cohort_30 | E2 |  | 0.685213 → 0.685213 | 0.635556 → 0.635556 | astrocyte | True → True | holds → holds |
| primary | cohort_30 | E3 |  | 0.688175 → 0.688175 | 0.843828 → 0.843828 | mitochondrial_pseudobulk | False → False | fails → fails |
| primary | cohort_30 | overall |  |  |  |  | False → False | does_not_hold → does_not_hold |
| primary | cohort_50 | E1 | deseq2 | 35 → 35 | 2 → 2 | astrocyte | True → True | holds → holds |
| primary | cohort_50 | E1 | limma | 0 → 0 | 1 → 1 | astrocyte | False → False | fails → fails |
| primary | cohort_50 | E1 | both |  |  |  | False → False | method_dependent → method_dependent |
| primary | cohort_50 | E2 |  | 0.656488 → 0.656488 | 0.696504 → 0.696504 | astrocyte | False → False | fails → fails |
| primary | cohort_50 | E3 |  | 0.730773 → 0.730773 | 0.784562 → 0.784562 | mitochondrial_pseudobulk | False → False | fails → fails |
| primary | cohort_50 | overall |  |  |  |  | False → False | does_not_hold → does_not_hold |
| latent_adjusted | cohort_30 | E1 | deseq2 | 61 → 61 | 1 → 1 | astrocyte | True → True | holds → holds |
| latent_adjusted | cohort_30 | E1 | limma | 2 → 2 | 1 → 1 | astrocyte | True → True | holds → holds |
| latent_adjusted | cohort_30 | E1 | both |  |  |  | True → True | holds → holds |
| latent_adjusted | cohort_30 | E2 |  | 0.65748 → 0.65748 | 0.643333 → 0.643333 | astrocyte | True → True | holds → holds |
| latent_adjusted | cohort_30 | E3 |  | 0.527877 → 0.527877 | 0.415256 → 0.415256 | tmm_rna | True → True | holds → holds |
| latent_adjusted | cohort_30 | overall |  |  |  |  | True → True | holds → holds |
| latent_adjusted | cohort_50 | E1 | deseq2 | 35 → 35 | 2 → 2 | astrocyte | True → True | holds → holds |
| latent_adjusted | cohort_50 | E1 | limma | 0 → 0 | 1 → 1 | astrocyte | False → False | fails → fails |
| latent_adjusted | cohort_50 | E1 | both |  |  |  | False → False | method_dependent → method_dependent |
| latent_adjusted | cohort_50 | E2 |  | 0.626085 → 0.626085 | 0.664479 → 0.664479 | astrocyte | False → False | fails → fails |
| latent_adjusted | cohort_50 | E3 |  | 0.532681 → 0.532681 | 0.406702 → 0.406702 | gex_reads_mapped_confidently_to_intronic_regions_sd | True → True | holds → holds |
| latent_adjusted | cohort_50 | overall |  |  |  |  | False → False | does_not_hold → does_not_hold |

## 11. Reading

- **The corrected primaries are the primary results.** In both cell types and both cohorts the corrected analysis rejects none of `delta`, `angle` and `shape` at α = 0.05. Under the binding rules: the `delta` non-rejections are weak evidence against a size difference: MOTCO measured the corrected `delta` as conservative at this configuration (0.029 / 0.031 pooled), and the adjusted analysis carries no power claim (the correction costs 0.014–0.054 of power on MOTCO's targets); the `angle` and `shape` non-rejections are uninformative about orientation and shape at this n, and the `angle` null is too wide to inform a non-rejection. The ranks (5, 2, 2, 4) are given beside `angle` with the eigengaps and null summaries (section 3); the floor guarantees two dimensions, not a second axis with signal. In `cohort_30` the corrected null's size is the one MOTCO measured at this configuration (its `cohort_30` design table, six covariates); in `cohort_50` it is inferred from the mechanism.
- **No dependence on the null correction.** No rejection status differs between a corrected analysis and its superseded one, in any step: the predeclared case "status differs" did not occur. Every p-value moved up or stayed, as κ > 1 with the same permutation indices predicts.
- **What the correction cannot change, reproduced.** The latent spaces, ranks, CV AUROCs, observed statistics, stage means and observed spectra of all 44 analyses; `pls_0` and the axis diagnostics; the adoption check and precedence; the mediator, flag and axis checks; E1–E3 and both bases' control outcomes; every stage-DE concordance cell and every exported count and call. The only reproduction beyond rounding concerns the machine-dependent per-feature DESeq2 statistics, with one recorded abundance-set ATAC count moved by a degenerate fit (section 12); it changes no reading.
- **Sensitivity readings, re-read against the corrected primary under their own rules.** Rank: no forced rank rejects; no rank dependence. Preprocessing: no variant rejects; the adoption rule holds for no candidate in both cohorts and the precedence names none, so the variants stay sensitivity analyses. Latent adjustment: no arm rejects; the `extended_covariates` arms keep their primaries' rank and AUROC; the `latent_adjusted` arm's `angle` and `shape` are reported and not read (`residualized_basis`), its `delta` is not rejected; the mediator limitation stands with its flags. Donor exclusion: no variant rejects; the microglia axis check holds and the astrocyte one does not, as before. Stage DE: every concordance cell concordant, methods agree.
- **Positive control.** On the corrected primary the control **does not hold** in `cohort_30`, through E3 alone, as on the superseded primary; the restated control on `latent_adjusted` **holds** in `cohort_30`, as before, so the measured quality gradient still explains the primary-basis failure, without promoting the arm and with the mediator caveat; `cohort_50` stays descriptive on both bases. The control is about the stage signal, not the sex × stage hypothesis.

## 12. Findings

- **Re-issue machine.** `nixos-home` was unreachable, so the re-issue ran on `nixos-dell-laptop` after a verified copy of the cluster's current outputs (design D10's fallback). The gates and the reproducibility check compare runs made here and are byte-identical. Every comparison with the superseded outputs is across machines: the 44 analyses reproduce at rounding level (the four CSVs within 1.3e-10; the records outside the null keys within 5.3e-11 relative, apart from the numerically-zero third eigenvalues of the configuration spectra, ≤ 1.4e-31), with every rank and CV AUROC equal.
- **DESeq2 across machines.** The 96 per-feature DESeq2 tables differ beyond rounding (section 8.3), as the astrocyte tables did in the 2026-10-26 re-issue; every exported call, every exported count and every concordance cell is equal. The limma tables reproduce at rounding level.
- **Degenerate pydeseq2 fits in the microglia ATAC tables.** A few non-exported, near-zero-count microglia ATAC peaks (baseMean 0.4–0.7) have collapsed standard errors (lfcSE < 1e-10, |Wald statistic| > 1e6) and padj 0 in both the superseded and the corrected tables, in the primary and the exclusion roots; they enter only the abundance-set ATAC counts, which are recorded, not read. This machine adds one: `chr13:40895349-40895968` in microglia `cohort_50` `factor_interaction` (superseded lfcSE 2.95 and 2.89, padj ≈ 1; corrected lfcSE 2.3e-71 and 4.5e-57, padj 0), which moves the abundance-set ATAC counts of `male_x_stage_1`, `male_x_stage_2`, their union and limma's method comparison (DESeq2 side) from 3 to 4. No exported count, concordance cell, E1 value or reading changes. Pre-existing, recorded here, not acted on.
- **Variant exports to the last digit.** The exports of the four variants that residualize or trend-normalize differ from their superseded copies in the last printed digit only (|difference| ≤ 1e-8; `residualized` `cohort_50` `rna.csv` 1.0e-9, the comparison script's one "differs" above its 1e-9 rounding cut), with the same features selected; the other variant exports and both exclusion roots' exports are byte-identical.
- **The p-values moved as predicted.** Up in 131 of 132 statistic × analysis cells and unchanged in one (`atac_mean_accessibility` `cohort_30` `delta`, 0.99); none moved down.

Generated by the comparison script (superseded against corrected; every difference listed, none a gate):

- outputs not byte-identical to their superseded copies: 337; rejection-status changes: 0. By output kind and level:

| output | byte-identical | equal at rounding level (≤ 1e-9) | different beyond rounding |
|---|---:|---:|---:|
| adoption check | 6 | 0 | 0 |
| analysis CSVs (latent scores, loadings, stage means, metadata) | 0 | 176 | 0 |
| analysis records outside the null keys | 0 | 44 | 0 |
| axis check | 4 | 0 | 0 |
| axis diagnostics | 5 | 11 | 0 |
| axis diagnostics (exclusion) | 5 | 11 | 0 |
| axis diagnostics (latent) | 13 | 27 | 0 |
| axis diagnostics (variant) | 15 | 49 | 0 |
| exclusion export | 24 | 0 | 0 |
| flag check | 4 | 0 | 0 |
| mediator check | 4 | 0 | 0 |
| positive_control.csv | 4 | 0 | 0 |
| stage-DE per-feature tables, DESeq2 | 0 | 0 | 96 |
| stage-DE per-feature tables, limma | 48 | 48 | 0 |
| stage-DE summaries and method comparisons | 3 | 5 | 16 |
| variant export | 78 | 27 | 1 |

The outputs different beyond rounding are the DESeq2 tables and the summaries that carry their `min_padj` (section 8.3), the microglia `cohort_50` summaries whose abundance-set ATAC interaction count moved through the degenerate fit above, and `residualized` `cohort_50` `rna.csv` (last printed digit). Every adoption check, axis check, flag check, mediator check, exclusion export and positive-control table is byte-identical.



## 13. Limitations

- The `df`-corrected null's size was measured by MOTCO at six covariates on the two `cohort_30` tables only (astrocyte n = 80, microglia n = 76). It is inferred from the mechanism, which depends only on n and the reduced rank, for `cohort_50` (n = 74 and 63, κ 1.0753 and 1.0903), the eight-covariate latent arms (r = 12), the exclusion variants (n = 58–70), the forced ranks and the preprocessing variants. No Type I claim is made beyond MOTCO's measurement.
- `delta` is conservative under `df` at the measured configuration (0.029 / 0.031 pooled), and the correction costs 0.014–0.054 of power on MOTCO's targets; a non-rejection is weak evidence, and the adjusted analysis carries no power claim at these blocks.
- MOTCO's `residualized_basis` diagnosis makes the `latent_adjusted` arm inadmissible for `angle` and `shape` at this configuration; their p-values are reported and not read, and the residual correction does not address the excess.
- The floor's limitations carry over: it defines `angle` and `shape` in at least two dimensions without making them powerful, and the forced-rank operating characteristics (Type I at forced 4 only, unadjusted) stay limitations of the sweep.
- Every limitation of the earlier reports carries over unchanged: the additive-linear adjustment on a covariate-blind basis; the structural-analogue operating characteristics (InterSIM blocks, 498 features); the mediator limitation of the latent adjustment; the exclusion variant's truncation of the late epoch and quality conditioning; the stage-DE methods' approximate FDR and no joint error-rate claim for the unions.
- The per-feature DESeq2 statistics depend on the machine beyond rounding, and pydeseq2 fits of near-zero-count features can degenerate (collapsed standard errors; section 12); exported counts and calls did not change in this comparison.
