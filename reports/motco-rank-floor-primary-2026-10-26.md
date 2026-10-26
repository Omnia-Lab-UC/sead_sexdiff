# MOTCO rank-floor primary: the case study re-issued under MOTCO 0.10.0 (latent-rank floor 2), astrocytes and microglia — 2026-10-26

> The re-issue of OpenSpec change `mtg-rank-floor-primary` (specs `motco-trajectory-analysis`, `motco-rank-sensitivity`, `cell-type-positive-control`, new `motco-release-migration`): MOTCO pinned to `v0.10.0`, whose cross-validated latent space selects the rank from candidates at or above a floor of 2 (`min_components = 2`, MOTCO's validated default; `ANALYSIS_RANK_FLOOR = 2` here, passed explicitly, checked against MOTCO's default and by every reader of a reference analysis). The floored analysis is the primary analysis of both cell types. A same-machine gate showed that at floor 1 the pin move changes no artifact; every MOTCO-dependent output of both cell types was then moved to `results/superseded/mtg-rank-floor-primary/` and re-made from scratch: the primary analyses, the rank sweep, the axis diagnostics, the eight preprocessing variants with their adoption check and precedence, the latent-adjustment arms, both stage-DE steps, the donor-exclusion variant with its diagnostics and stage DE, and the positive control on both bases. The superseded values are given beside the floored ones throughout; the earlier reports stay as dated records.
>
> **Result.** The floored primaries select rank 5 (astrocyte `cohort_30`, CV mean AUROC 0.6356), rank 2 (astrocyte `cohort_50`, 0.6965), rank 2 (microglia `cohort_30`, 0.6852) and rank 4 (microglia `cohort_50`, 0.6565), where every superseded primary selected rank 1. `angle` is now measured in at least two dimensions (21.9941°, 23.1596°, 12.5272°, 22.2604°, p 0.396, 0.222, 0.424, 0.227) and, like `delta` and `shape`, **is rejected in no analysis**: not in the four primaries, the twelve forced ranks, the sixteen preprocessing-variant analyses, the eight latent-adjustment arms or the four donor-exclusion analyses. No rejection status differs from its superseded value anywhere, so the floor shows no dependence on the rank rule at α = 0.05. The first latent component reproduces in all 40 analyses that have a superseded copy, so every first-axis reading (E3, the axis check, the flag check's AUCs, the diagnostics' first-axis correlations) is unchanged. The positive control's outcome is unchanged on both bases: on the floored primary it does **not hold** in `cohort_30` (E1 holds, E2 holds at 0.6852 against 0.6356, E3 fails at 0.6882 against `mitochondrial_pseudobulk` 0.8438), and the restated control on `latent_adjusted` **holds** in `cohort_30` (E2 0.6575 against 0.6433, E3 0.5279 against `tmm_rna` 0.4153). No preprocessing candidate meets the adoption rule in both cohorts and the precedence names none, as before; `revised_candidate`'s `cohort_30` cell now fails because its floored analysis (rank 3) exposes two more components to the unchanged rule. Every stage-DE and limma concordance cell is concordant. One expected reproduction failed: the 24 per-feature astrocyte DESeq2 tables, superseded copies made on another workstation, differ beyond rounding (largest |Δ stat| 0.316), with every padj ≤ 0.05 call and count equal.

**Reading rules (binding, predeclared in the proposal before any floored run).** The floored analysis is the primary of each cell type; the unfloored analyses are superseded records, not a sensitivity arm, and are not cited as results. `cohort_30` stays primary and `cohort_50` sensitivity. The binding rules hold unchanged: the CV-selected rank (now ≥ 2) beside `angle`; `angle` and `shape` non-rejections uninformative; a `delta` rejection a difference, not specifically of size. The floor defines `angle` and `shape` in at least two dimensions; it does not make them powerful (MOTCO's floor check found power at the null at its strengths), so each analysis's `config_spectrum` eigengap and null summary are given beside `angle`. A rejection status that differs between a floored analysis and its superseded one would be dependence on the rank rule, with the floored status the primary result; nothing is re-tuned. Every sensitivity reading is re-read against the floored primary under its own unchanged rules; no registry, expectation, adoption rule or threshold changes. The positive control's floored outcome is the primary-basis outcome from this change on, with the superseded outcome beside it; the restated control on `latent_adjusted` keeps its predeclared reading.

## 1. Provenance

- Software: MOTCO 0.10.0 (git dependency pinned to `v0.10.0`, commit `098b170e618959591c7a3d2277487765dcd0459a` in `uv.lock`), scikit-learn 1.8.0, numpy 2.3.5, scipy 1.15.3, pandas 2.3.3, pydeseq2 0.5.3, Python 3.12.15; R 4.6.1 with limma 3.68.4, edgeR 4.10.1 and jsonlite 2.0.0 (the Nix shell; `stage_de_limma.R` sha256 `0543a17dc2d12ba77ef00d27b53b08a68095b8b1ce00c7b98dc63e31eb5d70a8`). sead-sexdiff revision `4eff306c1b9d7557906cede718deaf8192cca63f` for every re-issued output (recorded `dirty: True`: the uncommitted tracked files at run time were the IDE file `.idea/sead-sexdiff.iml` and, from the rank sweep on, this change's `tasks.md` as its boxes were ticked; nothing under `src/`). Machine `nixos-home` (Nix shell), `n_jobs = 1` throughout, `--n-cpus 10` for `stage-de`.

### 1.1 Same-machine pin-move gate (before any floored run)

The four primary analyses were run on `nixos-home` on the real exports twice: from the pre-change commit `b1a67a7` (a scratch worktree, `uv sync --locked`, MOTCO 0.9.0) with `analyze --cell-type {astrocyte,microglia} --cohort all`, and from this change's code (MOTCO 0.10.0) through `_analyze_cohort` with `min_components = 1`. Compared: `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv`, `analysis_metadata.csv` byte for byte; `result.json` and `analysis_manifest.json` after dropping `runtime_seconds`, `created_at`, `git`, `software.motco`, `cli_equivalent` and `measurement.min_components`. Every 0.9.0 manifest records MOTCO 0.9.0 and the real exports' sha256; every floor-1 record carries no `min_components` in `params` and MOTCO 0.10.0 under `software`.

| cell type / cohort | latent_scores.csv | pls_loadings.csv | stage_means.csv | analysis_metadata.csv | result.json | analysis_manifest.json | gate |
|---|---|---|---|---|---|---|---|
| astrocyte/cohort_30 | identical | identical | identical | identical | equal | equal | pass |
| astrocyte/cohort_50 | identical | identical | identical | identical | equal | equal | pass |
| microglia/cohort_30 | identical | identical | identical | identical | equal | equal | pass |
| microglia/cohort_50 | identical | identical | identical | identical | equal | equal | pass |

Gate: **pass** — at floor 1 the pin move changes no artifact of the primary analyses, so every difference below comes from the floor.

### 1.2 Floor-1 rerun against the kept primaries (recorded, not a gate)

- microglia cohort_30 (kept primary made on this machine): the four CSVs byte-identical; `result.json` equal after the dropped keys: True; rank 1, p-values equal: True.
- microglia cohort_50 (kept primary made on this machine): the four CSVs byte-identical; `result.json` equal after the dropped keys: True; rank 1, p-values equal: True.
- astrocyte cohort_30 (kept primary made on the cluster at `ef388f7`, 2026-10-21T15:45:08+00:00): `latent_scores.csv` max |diff| 2.63e-13 (largest |value| 87.2; column(s) pls_0); `pls_loadings.csv` max |diff| 1.01e-16 (largest |value| 0.0327; column(s) pls_0); `stage_means.csv` max |diff| 9.66e-13 (largest |value| 24.1; column(s) pls_0); `analysis_metadata.csv` max |diff| 4.44e-16 (largest |value| 102; column(s) log10_n_atac); statistics equal to 1.3e-13 relative, null summaries to 2e-13 relative; p-values equal: True; CV-selected rank 1 on both; CV AUROC equal: True. `analysis_metadata.csv` differs only in the last bit of three `log10_n_atac` values (`np.log10` between machines), which carries through to the scores at rounding level.
- astrocyte cohort_50 (kept primary made on the cluster at `ef388f7`, 2026-10-21T15:46:11+00:00): `latent_scores.csv` max |diff| 7.82e-14 (largest |value| 84.6; column(s) pls_0); `pls_loadings.csv` max |diff| 1.01e-16 (largest |value| 0.0324; column(s) pls_0); `stage_means.csv` max |diff| 5.97e-13 (largest |value| 23; column(s) pls_0); `analysis_metadata.csv` max |diff| 4.44e-16 (largest |value| 102; column(s) log10_n_atac); statistics equal to 9.4e-14 relative, null summaries to 7.4e-14 relative; p-values equal: True; CV-selected rank 1 on both; CV AUROC equal: True. `analysis_metadata.csv` differs only in the last bit of three `log10_n_atac` values (`np.log10` between machines), which carries through to the scores at rounding level.

### 1.3 Superseded outputs moved

Before the reruns, `scripts/move_superseded.py` moved 15 directories (956 files) to `results/superseded/mtg-rank-floor-primary/` with their paths relative to `results/` kept (`motco_analysis` (14), `motco_rank_sensitivity` (90), `motco_diagnostics` (12), `motco_sensitivity` (210), `motco_latent_sensitivity` (72), `stage_de` (46), `stage_de_limma` (56), `motco_exclusion_sensitivity` (128), `microglia/motco_analysis` (14), `microglia/motco_diagnostics` (12), `microglia/motco_latent_sensitivity` (60), `microglia/stage_de` (46), `microglia/stage_de_limma` (56), `microglia/motco_exclusion_sensitivity` (128), `microglia/positive_control` (12); absent here: `microglia/motco_rank_sensitivity`, `microglia/motco_sensitivity`), and `results/motco_analysis_0.8.0/` (14 files) to `results/superseded/mtg-astrocyte-latent-adjustment-sensitivity/motco_analysis/`, each recorded in `MOVED.json`. A snapshot of every file under `results/` and `data/` taken before the move (1045 files) was re-hashed after it: every file hashes as before at its new path, no MOTCO-dependent file is left outside `results/superseded/`, and all 42 sha256 values cited by the earlier reports that resolved to a file before the move still do (the 43rd is `src/sead_sexdiff/stage_de_limma.R`), including the 8 cited `result.json` hashes, now under `results/superseded/`.

### 1.4 Checks after the re-issue

- Version consistency (design D7): outside `results/superseded/`, 156 manifests (86 of which record a MOTCO version) and 44 `result.json` (12 forced, 32 CV-selected) were read; every recorded MOTCO version is 0.10.0 and every CV-selected record carries `params.min_components = 2`. No stale manifest.
- MOTCO-independent outputs: the 75 files under `results/pseudobulk/` (27), `results/motco_input/` (14), `results/microglia/pseudobulk/` (19), `results/microglia/motco_input/` (14) and `data/annotation/` (1), hashed before the move, hash unchanged after the last rerun.
- Variant exports: the 100 export CSVs of the variant roots that have a superseded copy (six preprocessing variants, both donor-exclusion roots) are byte-identical to it; the 30 CSVs of `trend_normalized` and `revised_candidate_trend` had no superseded copy on this machine.
- Reproducibility (`analyze --cohort cohort_30 --force` into a scratch directory, per cell type): astrocyte: the four CSVs byte-identical: True, statistics, p-values and null summaries equal: True, rank 5 on both, `config_spectrum` equal: True; microglia: the four CSVs byte-identical: True, statistics, p-values and null summaries equal: True, rank 2 on both, `config_spectrum` equal: True.

### 1.5 Commands

```
# gate (scratch worktree of b1a67a7, MOTCO 0.9.0; then this change, MOTCO 0.10.0, floor 1 through _analyze_cohort)
sead-sexdiff analyze --cell-type {astrocyte,microglia} --cohort all --input-dir <real export> --out-dir <scratch>/v090/<cell type>
# move
python3 scripts/move_superseded.py --change mtg-rank-floor-primary motco_analysis motco_rank_sensitivity motco_diagnostics \
    motco_sensitivity motco_latent_sensitivity stage_de stage_de_limma motco_exclusion_sensitivity microglia/motco_analysis \
    microglia/motco_rank_sensitivity microglia/motco_diagnostics microglia/motco_sensitivity microglia/motco_latent_sensitivity \
    microglia/stage_de microglia/stage_de_limma microglia/motco_exclusion_sensitivity microglia/positive_control
python3 scripts/move_superseded.py --change mtg-astrocyte-latent-adjustment-sensitivity motco_analysis_0.8.0=motco_analysis
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
sead-sexdiff axis-diagnostics --cohort all --analysis-dir <primary or arm root> --rank-dir <absent> --out-dir results/motco_latent_sensitivity/diagnostics/<name>
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
# reproducibility
sead-sexdiff analyze --cohort cohort_30 --n-jobs 1 --force --out-dir <scratch>/repro/astrocyte    # and --cell-type microglia
```

## 2. The floored primary analyses

The four primary analyses, run with every setting of the superseded ones (blocks, covariates, orientation, CV settings, 999 permutations, seed 20261019, `n_jobs = 1`) and the rank floor 2. Each selects a rank of at least 2; none rejects `delta`, `angle` or `shape` at α = 0.05. The `angle` rows are now orientation statistics in two to five dimensions; their non-rejections are uninformative about orientation at this n, not evidence that the sexes' directions agree. The relative eigengap ((λ1 − λ2) / Σλ of the stage-mean configuration, MOTCO's `config_spectrum`) says how strongly one axis dominates the configuration and so how well a principal orientation could be resolved; it is given beside the RRPP permutation distribution of the same quantity, as MOTCO advises, and is read as geometry, not as a test. A superseded rank-1 record has eigengap 1 and an `angle` of 0° with p = 1 by construction.

### 2.1 Microglia

#### cohort_30 (≥ 30 microglia nuclei in both omics, n = 76) — primary cohort

Pasted from `results/microglia/motco_analysis/cohort_30/summary.md`:

##### microglia, cohort_30

MOTCO 0.10.0; n = 76; groups Female vs Male; stages 0 → 1 → 2; orientation `net_displacement`; 999 RRPP permutations, seed 20261019; α = 0.05 (rejected when p ≤ α).

###### Statistics

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 5.10924 | 0.689 | no | 10.2655 | 7.68786 | 25.6788 |  |
| angle | 12.5272 | 0.424 | no | 12.6966 | 9.51188 | 30.1659 | degrees; CV-selected rank 2 (rank floor 2) |
| shape | 0.0242177 | 0.99 | no | 0.259273 | 0.142347 | 0.515887 |  |

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

##### microglia, cohort_50

MOTCO 0.10.0; n = 63; groups Female vs Male; stages 0 → 1 → 2; orientation `net_displacement`; 999 RRPP permutations, seed 20261019; α = 0.05 (rejected when p ≤ α).

###### Statistics

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 1.6686 | 0.819 | no | 5.78012 | 4.4338 | 14.0167 |  |
| angle | 22.2604 | 0.227 | no | 16.6216 | 7.94278 | 31.3913 | degrees; CV-selected rank 4 (rank floor 2) |
| shape | 0.0797433 | 0.598 | no | 0.0979335 | 0.0515193 | 0.196376 |  |

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

| analysis | rank | pooled spectrum | eigengap (pooled) | eigengap F / M | permutation eigengap mean (q05–q95) | angle | angle p | angle null mean / sd / q95 |
|---|---|---|---|---|---|---|---|---|
| microglia cohort_30 | 2 | 0.8876, 0.1124 | 0.775273 | 0.7892 / 0.7587 | 0.7720 (0.5003–0.9770) | 12.5272 | 0.424 | 12.6966 / 9.51188 / 30.1659 |
| microglia cohort_50 | 4 | 0.6124, 0.3876, 0.0000 | 0.224779 | 0.2882 / 0.1919 | 0.2328 (0.1237–0.3448) | 22.2604 | 0.227 | 16.6216 / 7.94278 / 31.3913 |


#### Superseded → floored (microglia)

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| microglia cohort_30 | 1 → 2 | 0.6921 → 0.6852 | 8.0642 → 5.10924 | 0.391 → 0.689 | no → no | 0 → 12.5272 | 1 → 0.424 | no → no | 0.154684 → 0.0242177 | 0.637 → 0.99 | no → no |
| microglia cohort_50 | 1 → 4 | 0.6414 → 0.6565 | 12.4641 → 1.6686 | 0.227 → 0.819 | no → no | 0 → 22.2604 | 1 → 0.227 | no → no | 0.269568 → 0.0797433 | 0.399 → 0.598 | no → no |


### 2.2 Astrocyte

#### cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

Pasted from `results/motco_analysis/cohort_30/summary.md`:

##### astrocyte, cohort_30

MOTCO 0.10.0; n = 80; groups Female vs Male; stages 0 → 1 → 2; orientation `net_displacement`; 999 RRPP permutations, seed 20261019; α = 0.05 (rejected when p ≤ α).

###### Statistics

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |  |
| angle | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 | degrees; CV-selected rank 5 (rank floor 2) |
| shape | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |  |

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

##### astrocyte, cohort_50

MOTCO 0.10.0; n = 74; groups Female vs Male; stages 0 → 1 → 2; orientation `net_displacement`; 999 RRPP permutations, seed 20261019; α = 0.05 (rejected when p ≤ α).

###### Statistics

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |  |
| angle | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 | degrees; CV-selected rank 2 (rank floor 2) |
| shape | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |  |

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

| analysis | rank | pooled spectrum | eigengap (pooled) | eigengap F / M | permutation eigengap mean (q05–q95) | angle | angle p | angle null mean / sd / q95 |
|---|---|---|---|---|---|---|---|---|
| astrocyte cohort_30 | 5 | 0.5854, 0.4146, 0.0000 | 0.170868 | 0.2818 / 0.0461 | 0.2143 (0.0792–0.3615) | 21.9941 | 0.396 | 20.7155 / 9.27804 / 37.7882 |
| astrocyte cohort_50 | 2 | 0.9850, 0.0150 | 0.969992 | 0.9989 / 0.8476 | 0.9275 (0.7444–0.9997) | 23.1596 | 0.222 | 15.1566 / 12.4608 / 40.0109 |


#### Superseded → floored (astrocyte)

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| astrocyte cohort_30 | 1 → 5 | 0.6313 → 0.6356 | 14.941 → 5.95773 | 0.3 → 0.468 | no → no | 0 → 21.9941 | 1 → 0.396 | no → no | 0.608832 → 0.15645 | 0.333 → 0.284 | no → no |
| astrocyte cohort_50 | 1 → 2 | 0.6616 → 0.6965 | 18.7091 → 6.01569 | 0.213 → 0.596 | no → no | 0 → 23.1596 | 1 → 0.222 | no → no | 0.56471 → 0.284932 | 0.377 → 0.484 | no → no |


## 3. Rank sensitivity (astrocytes)

The forced-rank sweep reruns the analysis at the fixed ladder 2, 3, 4, 6, 9, 12 with every other setting of the primary; a forced rank skips cross-validation, so the floor does not act on it, and every forced record reproduces its superseded one exactly (statistics, p-values, null summaries; `pls_0` max |diff| 0). The reference row is now the floored primary. In `cohort_30` the CV-selected rank 5 is not in the ladder, so the reference row stands in for it (the manifest's finding); in `cohort_50` the forced rank-2 run equals the reference (statistics and p-values). No forced rank rejects any statistic and no rejection status differs from the reference: no rank dependence. MOTCO's forced-rank operating characteristics (Type I at forced 4 only, unadjusted; `shape` power loss and `angle` leakage at forced 2(K − 1); a wider `angle` null above `n_stages − 1`) remain limitations, and the adjusted analysis has no Type I check at any forced rank.

### 3.1 cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

Pasted from `results/motco_rank_sensitivity/cohort_30/rank_sensitivity.md`:

#### MOTCO rank sensitivity: cohort_30

Forced-rank sweep of the covariate-adjusted trajectory analysis. The reference row (`cv`) is the CV-selected analysis at rank 5, which stays the primary result; the forced ranks 2, 3, 4, 6, 9, 12 are read as sensitivity to it, and no rank is chosen from the results. 999 RRPP permutations with seed 20261019 at every rank (the same permutation indices throughout); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative at every rank; a `delta` rejection is a difference, not specifically of size.

##### delta

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| forced | 2 | 3.91847 | 0.744 | no | 9.63089 | 7.18465 | 23.2427 |
| forced | 3 | 10.2517 | 0.338 | no | 8.48407 | 6.74633 | 21.5709 |
| forced | 4 | 19.247 | 0.052 | no | 7.78873 | 6.01757 | 19.3458 |
| forced | 6 | 6.08738 | 0.369 | no | 5.69703 | 4.5972 | 14.8882 |
| forced | 9 | 1.79797 | 0.778 | no | 5.26197 | 4.29954 | 13.432 |
| forced | 12 | 1.15591 | 0.866 | no | 5.32251 | 4.29312 | 13.705 |

##### angle (degrees)

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| forced | 2 | 22.832 | 0.241 | no | 15.7454 | 12.6226 | 40.5824 |
| forced | 3 | 23.0297 | 0.316 | no | 18.8414 | 10.5607 | 39.0162 |
| forced | 4 | 22.1423 | 0.356 | no | 19.6796 | 9.68581 | 37.4196 |
| forced | 6 | 21.9664 | 0.455 | no | 22.0223 | 8.9536 | 38.8732 |
| forced | 9 | 31.5192 | 0.21 | no | 25.4414 | 8.08815 | 40.3043 |
| forced | 12 | 35.8438 | 0.152 | no | 28.4101 | 7.42662 | 41.5921 |

##### shape

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| forced | 2 | 0.295865 | 0.445 | no | 0.290846 | 0.158807 | 0.589015 |
| forced | 3 | 0.121943 | 0.789 | no | 0.22084 | 0.117102 | 0.427894 |
| forced | 4 | 0.0296112 | 0.975 | no | 0.162616 | 0.08214 | 0.307316 |
| forced | 6 | 0.153228 | 0.224 | no | 0.111023 | 0.0578495 | 0.215647 |
| forced | 9 | 0.147921 | 0.161 | no | 0.0957448 | 0.0508231 | 0.187395 |
| forced | 12 | 0.131817 | 0.208 | no | 0.0925515 | 0.0492176 | 0.181653 |

##### Runs

| selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|
| cv | 5 | 0.6356 | 115.2 | `results/motco_analysis/cohort_30` |
| forced | 2 | n/a | 0.8 | `rank_02` |
| forced | 3 | n/a | 1.1 | `rank_03` |
| forced | 4 | n/a | 0.8 | `rank_04` |
| forced | 6 | n/a | 1.0 | `rank_06` |
| forced | 9 | n/a | 1.9 | `rank_09` |
| forced | 12 | n/a | 1.4 | `rank_12` |


- Reference: CV-selected rank 5 at rank floor 2, in the ladder: False; forced run at the CV-selected rank: none (not in the ladder).
- Manifest findings: CV-selected rank 5 is not in the ladder (2, 3, 4, 6, 9, 12); the reference row stands in for it.
- `pls_0` across the primary and the forced ranks (`axis-diagnostics` nestedness): rank_02 max |diff| 0, rank_03 max |diff| 0, rank_04 max |diff| 0, rank_06 max |diff| 0, rank_09 max |diff| 0, rank_12 max |diff| 0.

Superseded → floored (cohort_30; a forced rank's record is unchanged by the floor, so its row reproduces):

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CV-selected | 1 → 5 | 0.6313 → 0.6356 | 14.941 → 5.95773 | 0.3 → 0.468 | no → no | 0 → 21.9941 | 1 → 0.396 | no → no | 0.608832 → 0.15645 | 0.333 → 0.284 | no → no |
| forced 2 | 2 → 2 | n/a → n/a | 3.91847 → 3.91847 | 0.744 → 0.744 | no → no | 22.832 → 22.832 | 0.241 → 0.241 | no → no | 0.295865 → 0.295865 | 0.445 → 0.445 | no → no |
| forced 3 | 3 → 3 | n/a → n/a | 10.2517 → 10.2517 | 0.338 → 0.338 | no → no | 23.0297 → 23.0297 | 0.316 → 0.316 | no → no | 0.121943 → 0.121943 | 0.789 → 0.789 | no → no |
| forced 4 | 4 → 4 | n/a → n/a | 19.247 → 19.247 | 0.052 → 0.052 | no → no | 22.1423 → 22.1423 | 0.356 → 0.356 | no → no | 0.0296112 → 0.0296112 | 0.975 → 0.975 | no → no |
| forced 6 | 6 → 6 | n/a → n/a | 6.08738 → 6.08738 | 0.369 → 0.369 | no → no | 21.9664 → 21.9664 | 0.455 → 0.455 | no → no | 0.153228 → 0.153228 | 0.224 → 0.224 | no → no |
| forced 9 | 9 → 9 | n/a → n/a | 1.79797 → 1.79797 | 0.778 → 0.778 | no → no | 31.5192 → 31.5192 | 0.21 → 0.21 | no → no | 0.147921 → 0.147921 | 0.161 → 0.161 | no → no |
| forced 12 | 12 → 12 | n/a → n/a | 1.15591 → 1.15591 | 0.866 → 0.866 | no → no | 35.8438 → 35.8438 | 0.152 → 0.152 | no → no | 0.131817 → 0.131817 | 0.208 → 0.208 | no → no |


### 3.2 cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74) — sensitivity cohort

Pasted from `results/motco_rank_sensitivity/cohort_50/rank_sensitivity.md`:

#### MOTCO rank sensitivity: cohort_50

Forced-rank sweep of the covariate-adjusted trajectory analysis. The reference row (`cv`) is the CV-selected analysis at rank 2, which stays the primary result; the forced ranks 2, 3, 4, 6, 9, 12 are read as sensitivity to it, and no rank is chosen from the results. 999 RRPP permutations with seed 20261019 at every rank (the same permutation indices throughout); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative at every rank; a `delta` rejection is a difference, not specifically of size.

##### delta

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| forced | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| forced | 3 | 9.36583 | 0.496 | no | 11.1486 | 8.59682 | 27.3245 |
| forced | 4 | 14.1822 | 0.217 | no | 9.01492 | 6.56113 | 21.4611 |
| forced | 6 | 10.7046 | 0.208 | no | 6.81048 | 5.21583 | 16.7231 |
| forced | 9 | 6.51697 | 0.386 | no | 6.1817 | 4.92541 | 15.8524 |
| forced | 12 | 5.88663 | 0.436 | no | 6.20365 | 4.94692 | 15.8513 |

##### angle (degrees)

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| forced | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| forced | 3 | 23.4732 | 0.299 | no | 19.2555 | 11.5529 | 42.4812 |
| forced | 4 | 24.5502 | 0.292 | no | 20.4756 | 10.6534 | 42.0583 |
| forced | 6 | 25.2911 | 0.33 | no | 22.5825 | 9.44256 | 40.8348 |
| forced | 9 | 30.5012 | 0.251 | no | 25.7816 | 8.57882 | 42.7693 |
| forced | 12 | 32.5938 | 0.279 | no | 28.7713 | 7.98762 | 43.9592 |

##### shape

| selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|:---:|---:|---:|---:|
| cv | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| forced | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| forced | 3 | 0.16487 | 0.739 | no | 0.251509 | 0.126233 | 0.483667 |
| forced | 4 | 0.110409 | 0.719 | no | 0.168021 | 0.0888066 | 0.33709 |
| forced | 6 | 0.0845199 | 0.652 | no | 0.116727 | 0.062907 | 0.231105 |
| forced | 9 | 0.0826203 | 0.593 | no | 0.103863 | 0.0566925 | 0.208077 |
| forced | 12 | 0.0785119 | 0.598 | no | 0.100357 | 0.0553857 | 0.200242 |

##### Runs

| selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|
| cv | 2 | 0.6965 | 114.4 | `results/motco_analysis/cohort_50` |
| forced | 2 | n/a | 1.1 | `rank_02` |
| forced | 3 | n/a | 0.9 | `rank_03` |
| forced | 4 | n/a | 1.1 | `rank_04` |
| forced | 6 | n/a | 1.1 | `rank_06` |
| forced | 9 | n/a | 2.8 | `rank_09` |
| forced | 12 | n/a | 2.1 | `rank_12` |


- Reference: CV-selected rank 2 at rank floor 2, in the ladder: True; forced run at the CV-selected rank: {"rank": 2, "statistics_equal": true, "p_values_equal": true}.
- Manifest findings: none.
- `pls_0` across the primary and the forced ranks (`axis-diagnostics` nestedness): rank_02 max |diff| 0, rank_03 max |diff| 0, rank_04 max |diff| 0, rank_06 max |diff| 0, rank_09 max |diff| 0, rank_12 max |diff| 0.

Superseded → floored (cohort_50; a forced rank's record is unchanged by the floor, so its row reproduces):

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CV-selected | 1 → 2 | 0.6616 → 0.6965 | 18.7091 → 6.01569 | 0.213 → 0.596 | no → no | 0 → 23.1596 | 1 → 0.222 | no → no | 0.56471 → 0.284932 | 0.377 → 0.484 | no → no |
| forced 2 | 2 → 2 | n/a → n/a | 6.01569 → 6.01569 | 0.596 → 0.596 | no → no | 23.1596 → 23.1596 | 0.222 → 0.222 | no → no | 0.284932 → 0.284932 | 0.484 → 0.484 | no → no |
| forced 3 | 3 → 3 | n/a → n/a | 9.36583 → 9.36583 | 0.496 → 0.496 | no → no | 23.4732 → 23.4732 | 0.299 → 0.299 | no → no | 0.16487 → 0.16487 | 0.739 → 0.739 | no → no |
| forced 4 | 4 → 4 | n/a → n/a | 14.1822 → 14.1822 | 0.217 → 0.217 | no → no | 24.5502 → 24.5502 | 0.292 → 0.292 | no → no | 0.110409 → 0.110409 | 0.719 → 0.719 | no → no |
| forced 6 | 6 → 6 | n/a → n/a | 10.7046 → 10.7046 | 0.208 → 0.208 | no → no | 25.2911 → 25.2911 | 0.33 → 0.33 | no → no | 0.0845199 → 0.0845199 | 0.652 → 0.652 | no → no |
| forced 9 | 9 → 9 | n/a → n/a | 6.51697 → 6.51697 | 0.386 → 0.386 | no → no | 30.5012 → 30.5012 | 0.251 → 0.251 | no → no | 0.0826203 → 0.0826203 | 0.593 → 0.593 | no → no |
| forced 12 | 12 → 12 | n/a → n/a | 5.88663 → 5.88663 | 0.436 → 0.436 | no → no | 32.5938 → 32.5938 | 0.279 → 0.279 | no → no | 0.0785119 → 0.0785119 | 0.598 → 0.598 | no → no |


## 4. Axis diagnostics: the first latent component

`axis-diagnostics` was rerun on every re-issued analysis root (primary, latent arms, preprocessing variants, exclusion variant). Its input-hash checks passed everywhere. Because PLS components nest, the first latent component `pls_0` of every floored analysis equals the superseded one (40 comparisons, max |diff| ≤ 5.8e-13 against scores up to 148; the comparison table is in section 11), and the astrocyte diagnostics record `pls_0` as identical across the primary and the six forced ranks. The first-axis readings are therefore those of the earlier reports: in microglia `pls_0` tracks `mitochondrial_pseudobulk` more strongly than stage, in astrocytes the ribosomal/library-composition axis; the floor adds components but does not move the first one. Correlations along `pls_1` and higher, which the floored analyses now carry, are in each `axis_correlations.csv`; as before, an axis correlated with a technical variable is a measurement-space finding, not a reinterpretation of the statistics.

Generated from each `axis_correlations.csv` (rows `analysis = primary`, `component = pls_0`; technical variables: the technical table and the export's TMM factors and log10 library sizes):

| cell type | cohort | analysis | rank | stage ρ | largest technical |ρ| | technical variable | n flagged (pls_0, raw or partial) |
|---|---|---|---|---|---|---|---|
| microglia | cohort_30 | primary | 2 | 0.688175 | 0.843828 | mitochondrial_pseudobulk | 28 |
| microglia | cohort_30 | extended_covariates | 2 | 0.688175 | 0.843828 | mitochondrial_pseudobulk | 28 |
| microglia | cohort_30 | latent_adjusted | 3 | -0.527877 | 0.415256 | tmm_rna | 10 |
| microglia | cohort_30 | excluding_severely_affected | 4 | 0.820134 | 0.689149 | fraction_mitochondrial_umis_sd | 18 |
| microglia | cohort_50 | primary | 4 | 0.730773 | 0.784562 | mitochondrial_pseudobulk | 25 |
| microglia | cohort_50 | extended_covariates | 4 | 0.730773 | 0.784562 | mitochondrial_pseudobulk | 25 |
| microglia | cohort_50 | latent_adjusted | 3 | -0.532681 | 0.406702 | gex_reads_mapped_confidently_to_intronic_regions_sd | 8 |
| microglia | cohort_50 | excluding_severely_affected | 4 | 0.824801 | 0.674613 | fraction_mitochondrial_umis_sd | 22 |
| astrocyte | cohort_30 | primary | 5 | 0.491212 | 0.840858 | ribosomal_mean | 20 |
| astrocyte | cohort_30 | extended_covariates | 5 | 0.491212 | 0.840858 | ribosomal_mean | 20 |
| astrocyte | cohort_30 | latent_adjusted | 4 | -0.636917 | 0.382068 | log10_lib_atac | 3 |
| astrocyte | cohort_30 | excluding_severely_affected | 2 | 0.578399 | 0.728598 | ribosomal_mean | 18 |
| astrocyte | cohort_50 | primary | 2 | 0.503605 | 0.839822 | ribosomal_mean | 21 |
| astrocyte | cohort_50 | extended_covariates | 2 | 0.503605 | 0.839822 | ribosomal_mean | 21 |
| astrocyte | cohort_50 | latent_adjusted | 3 | -0.752715 | 0.326442 | log10_lib_atac | 2 |
| astrocyte | cohort_50 | excluding_severely_affected | 2 | 0.595467 | 0.671612 | ribosomal_mean | 22 |


## 5. Preprocessing sensitivity (astrocytes)

All eight registered variants were rerun in both cohorts (the export at the primary's thresholds and N, then the analysis at the primary's settings and the floor). The variant exports reproduce the superseded ones byte for byte where a superseded copy exists, so every difference below is the floor's. No variant rejects any statistic in either cohort, and no rejection status differs from the floored primary or from the superseded variant: no preprocessing dependence at α = 0.05. The variants' CV-selected ranks (2 to 7) differ from the primary's in several rows; each is reported beside `angle`, not read as better or worse.

The adoption rule (no latent component of the candidate's analysis with |Spearman| ≥ 0.3 against a check variable, in both cohorts) is unchanged and reads every component, so a floored candidate exposes more components to it. `revised_candidate` now fails in `cohort_30` (rank 3: `pls_1` × `tmm_rna` 0.3772, `pls_1` × `tmm_atac` −0.3995, `pls_2` × `tmm_atac` 0.3664), where its superseded rank-1 analysis held; it fails in `cohort_50` exactly as before (rank 2 in both). `trend_normalized` fails in both cohorts on `pls_0` itself (six check variables), as in the trend-normalization report. `revised_candidate_trend` fails in `cohort_30` (`pls_1` × `tmm_rna` −0.3225) and holds in `cohort_50`. `adoption_precedence` names no candidate ("neither trend_normalized nor revised_candidate_trend meets ADOPTION_RULE in both cohorts"), as before: no revised-primary change follows and the variants stay sensitivity analyses. The adoption outcome of every candidate over both cohorts is unchanged from the superseded records.

### 5.1 cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

Superseded → floored, per variant (the superseded variants were made against the MOTCO 0.8.0 analysis; `trend_normalized` and `revised_candidate_trend` had no superseded copy on this machine):

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| primary | 1 → 5 | 0.6313 → 0.6356 | 14.941 → 5.95773 | 0.3 → 0.468 | no → no | 0 → 21.9941 | 1 → 0.396 | no → no | 0.608832 → 0.15645 | 0.333 → 0.284 | no → no |
| keep_chrx | 1 → 4 | 0.6294 → 0.6341 | 14.9283 → 19.0474 | 0.301 → 0.055 | no → no | 0 → 22.3677 | 1 → 0.353 | no → no | 0.607073 → 0.0237388 | 0.337 → 0.988 | no → no |
| exclude_families | 1 → 5 | 0.6313 → 0.6369 | 14.9196 → 6.05273 | 0.301 → 0.455 | no → no | 0 → 21.9731 | 1 → 0.396 | no → no | 0.609158 → 0.156521 | 0.334 → 0.282 | no → no |
| rna_trend_ratio | 1 → 6 | 0.6326 → 0.6351 | 13.101 → 5.89665 | 0.354 → 0.446 | no → no | 0 → 22.2081 | 1 → 0.479 | no → no | 0.610539 → 0.144967 | 0.315 → 0.261 | no → no |
| atac_mean_accessibility | 1 → 7 | 0.6346 → 0.6706 | 2.2945 → 0.18785 | 0.933 → 0.99 | no → no | 0 → 18.6469 | 1 → 0.772 | no → no | 0.885871 → 0.102546 | 0.31 → 0.794 | no → no |
| residualized | 1 → 2 | 0.4078 → 0.3991 | 4.04454 → 10.6959 | 0.39 → 0.419 | no → no | 0 → 18.1967 | 1 → 0.308 | no → no | 0.0442514 → 0.380087 | 0.795 → 0.236 | no → no |
| revised_candidate | 1 → 3 | 0.4822 → 0.4945 | 1.45283 → 15.1154 | 0.829 → 0.391 | no → no | 0 → 29.1593 | 1 → 0.19 | no → no | 0.224695 → 0.185726 | 0.514 → 0.735 | no → no |
| trend_normalized | n/a → 4 | n/a → 0.5991 | n/a → 6.32863 | n/a → 0.63 | n/a → no | n/a → 24.277 | n/a → 0.417 | n/a → no | n/a → 0.0847488 | n/a → 0.856 | n/a → no |
| revised_candidate_trend | n/a → 3 | n/a → 0.4907 | n/a → 14.3836 | n/a → 0.166 | n/a → no | n/a → 16.7741 | n/a → 0.435 | n/a → no | n/a → 0.0792388 | n/a → 0.852 | n/a → no |


#### Variant `keep_chrx`

Pasted from `results/motco_sensitivity/keep_chrx/comparison/cohort_30/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_30, variant `keep_chrx`

chrX kept; chrY and chrM removed; every other export setting equal to the primary. Export setting(s) changed: `excluded_chromosomes` = ['chrY', 'chrM'] (primary: `excluded_chromosomes` = ['chrX', 'chrY', 'chrM']); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `keep_chrx` row (CV-selected rank 4) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| keep_chrx | cv | 4 | 19.0474 | 0.055 | no | 7.82005 | 6.04517 | 19.3555 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| keep_chrx | cv | 4 | 22.3677 | 0.353 | no | 19.7796 | 9.75014 | 37.6798 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| keep_chrx | cv | 4 | 0.0237388 | 0.988 | no | 0.163805 | 0.0827693 | 0.311155 |

###### Findings

- CV-selected rank differs: 5 in the primary analysis, 4 in the keep_chrx analysis; reported beside angle, not read as better or worse
- keep_chrx export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


#### Variant `exclude_families`

Pasted from `results/motco_sensitivity/exclude_families/comparison/cohort_30/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_30, variant `exclude_families`

RNA genes of the ribosomal, mitochondrial and hemoglobin sets (TECHNICAL_GENE_SETS) removed before the abundance filter; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'] (primary: `excluded_gene_sets` = []); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `exclude_families` row (CV-selected rank 5) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| exclude_families | cv | 5 | 6.05273 | 0.455 | no | 6.8158 | 5.28424 | 16.936 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| exclude_families | cv | 5 | 21.9731 | 0.396 | no | 20.7146 | 9.26825 | 37.723 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| exclude_families | cv | 5 | 0.156521 | 0.282 | no | 0.123928 | 0.0641013 | 0.246375 |

###### Findings

- exclude_families export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


#### Variant `rna_trend_ratio`

Pasted from `results/motco_sensitivity/rna_trend_ratio/comparison/cohort_30/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_30, variant `rna_trend_ratio`

RNA features ranked at the selection step by the ratio of their log2-CPM variance to a binned-median mean–variance trend instead of by raw variance; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'trend_ratio', 'atac': 'variance'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `rna_trend_ratio` row (CV-selected rank 6) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| rna_trend_ratio | cv | 6 | 5.89665 | 0.446 | no | 6.3606 | 5.04761 | 16.4416 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| rna_trend_ratio | cv | 6 | 22.2081 | 0.479 | no | 22.3244 | 8.66347 | 37.8822 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| rna_trend_ratio | cv | 6 | 0.144967 | 0.261 | no | 0.112787 | 0.0586212 | 0.220731 |

###### Findings

- CV-selected rank differs: 5 in the primary analysis, 6 in the rna_trend_ratio analysis; reported beside angle, not read as better or worse
- rna_trend_ratio export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


#### Variant `atac_mean_accessibility`

Pasted from `results/motco_sensitivity/atac_mean_accessibility/comparison/cohort_30/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_30, variant `atac_mean_accessibility`

ATAC peaks ranked at the selection step by mean log2-CPM (total accessibility) instead of by variance; the RNA block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'variance', 'atac': 'mean'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `atac_mean_accessibility` row (CV-selected rank 7) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| atac_mean_accessibility | cv | 7 | 0.18785 | 0.99 | no | 10.5688 | 8.7947 | 28.2406 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| atac_mean_accessibility | cv | 7 | 18.6469 | 0.772 | no | 28.5283 | 12.7374 | 53.9926 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| atac_mean_accessibility | cv | 7 | 0.102546 | 0.794 | no | 0.190476 | 0.0981652 | 0.378288 |

###### Findings

- CV-selected rank differs: 5 in the primary analysis, 7 in the atac_mean_accessibility analysis; reported beside angle, not read as better or worse
- atac_mean_accessibility export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


#### Variant `residualized`

Pasted from `results/motco_sensitivity/residualized/comparison/cohort_30/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_30, variant `residualized`

after log2-CPM and before the selection step, each abundance-passing feature replaced by its mean plus the residual of its regression on an intercept and RESIDUALIZATION_COVARIATES (the export's log10 library sizes and the technical table's gene-set fractions); the selection ranks the residual variance; every other export setting equal to the primary. Export setting(s) changed: `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'] (primary: `residualize` = None); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `residualized` row (CV-selected rank 2) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| residualized | cv | 2 | 10.6959 | 0.419 | no | 10.6588 | 7.88308 | 25.7971 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| residualized | cv | 2 | 18.1967 | 0.308 | no | 14.1035 | 10.5223 | 34.1378 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| residualized | cv | 2 | 0.380087 | 0.236 | no | 0.26892 | 0.145072 | 0.540248 |

###### Findings

- CV-selected rank differs: 5 in the primary analysis, 2 in the residualized analysis; reported beside angle, not read as better or worse
- residualized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


#### Variant `revised_candidate`

Pasted from `results/motco_sensitivity/revised_candidate/comparison/cohort_30/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_30, variant `revised_candidate`

the four settings together in the primary's order: gene sets removed before the abundance filter, residualization after log2-CPM, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `revised_candidate` row (CV-selected rank 3) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| revised_candidate | cv | 3 | 15.1154 | 0.391 | no | 14.4112 | 11.4564 | 37.3465 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| revised_candidate | cv | 3 | 29.1593 | 0.19 | no | 20.139 | 11.0734 | 40.0134 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| revised_candidate | cv | 3 | 0.185726 | 0.735 | no | 0.301173 | 0.159827 | 0.59153 |

###### Findings

- adoption check does not hold in cohort_30: pls_1 x tmm_rna (0.38), pls_1 x tmm_atac (-0.40), pls_2 x tmm_atac (0.37); the variants stay sensitivity analyses
- CV-selected rank differs: 5 in the primary analysis, 3 in the revised_candidate analysis; reported beside angle, not read as better or worse
- revised_candidate export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


`adoption_check.csv` (revised_candidate, cohort_30):

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---|---|---|---|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.18335677449601503 | False |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.16530707923112986 | False |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.1391701828410689 | False |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 80 | 0.0012658227848101266 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.21857840236194812 | False |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.12182372245663385 | False |  |
| pls_0 | rin | covariates.csv (primary export) | 80 | -0.033373971799641314 | False |  |
| pls_0 | pmi | covariates.csv (primary export) | 80 | -0.020000469609411455 | False |  |
| pls_1 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | 0.3772386310360994 | True |  |
| pls_1 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | -0.3995311767463667 | True |  |
| pls_1 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | -0.08661509610876701 | False |  |
| pls_1 | doublet_score_mean | donor_technical.csv | 80 | -0.09861697140178155 | False |  |
| pls_1 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | 0.1905660487963384 | False |  |
| pls_1 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | 0.15982184716361933 | False |  |
| pls_1 | rin | covariates.csv (primary export) | 80 | 0.13019014078216243 | False |  |
| pls_1 | pmi | covariates.csv (primary export) | 80 | 0.01766747227513661 | False |  |
| pls_2 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.02086263478668542 | False |  |
| pls_2 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.36640881387716834 | True |  |
| pls_2 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.11664322550398501 | False |  |
| pls_2 | doublet_score_mean | donor_technical.csv | 80 | -0.017510548523206753 | False |  |
| pls_2 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.1204531203321219 | False |  |
| pls_2 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.18223159868729488 | False |  |
| pls_2 | rin | covariates.csv (primary export) | 80 | -0.06405176744406046 | False |  |
| pls_2 | pmi | covariates.csv (primary export) | 80 | -0.15216066533443803 | False |  |
| all | all |  | 80 |  | True | False |


Adoption check holds in cohort_30: **False** (components pls_0, pls_1, pls_2; flagged: [{'component': 'pls_1', 'variable': 'tmm_rna', 'spearman': 0.3772386310360994}, {'component': 'pls_1', 'variable': 'tmm_atac', 'spearman': -0.3995311767463667}, {'component': 'pls_2', 'variable': 'tmm_atac', 'spearman': 0.36640881387716834}]).

#### Variant `trend_normalized`

Pasted from `results/motco_sensitivity/trend_normalized/comparison/cohort_30/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_30, variant `trend_normalized`

after log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features corrected by limma's fast cyclic loess (a lowess curve of its deviation from the per-feature mean over donors against that mean, three iterations, adaptive span) in both omics; every other export setting equal to the primary; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'} (primary: `normalization` = {'rna': 'none', 'atac': 'none'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `trend_normalized` row (CV-selected rank 4) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| trend_normalized | cv | 4 | 6.32863 | 0.63 | no | 10.8406 | 8.49241 | 26.7835 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| trend_normalized | cv | 4 | 24.277 | 0.417 | no | 22.6035 | 11.1794 | 42.4136 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| trend_normalized | cv | 4 | 0.0847488 | 0.856 | no | 0.184064 | 0.088751 | 0.333665 |

###### Findings

- adoption check does not hold in cohort_30: pls_0 x tmm_rna (-0.67), pls_0 x tmm_atac (0.51), pls_0 x fraction_mitochondrial_umis_mean (0.73), pls_0 x gex_reads_mapped_confidently_to_intronic_regions_mean (-0.56), pls_0 x gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63), pls_0 x rin (-0.33); the variants stay sensitivity analyses
- CV-selected rank differs: 5 in the primary analysis, 4 in the trend_normalized analysis; reported beside angle, not read as better or worse
- trend_normalized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


`adoption_check.csv` (trend_normalized, cohort_30):

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---|---|---|---|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.673230192217534 | True |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.5054149085794655 | True |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.7298406000937647 | True |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 80 | 0.01809657759024848 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.5649118900231207 | True |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.6321143928738866 | True |  |
| pls_0 | rin | covariates.csv (primary export) | 80 | -0.334618905873116 | True |  |
| pls_0 | pmi | covariates.csv (primary export) | 80 | 0.031630285466701116 | False |  |
| pls_1 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | 0.14753867791842476 | False |  |
| pls_1 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.040412564463197374 | False |  |
| pls_1 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | -0.10628223159868729 | False |  |
| pls_1 | doublet_score_mean | donor_technical.csv | 80 | -0.024706985466479137 | False |  |
| pls_1 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.10186416938858336 | False |  |
| pls_1 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.026207219878105953 | False |  |
| pls_1 | rin | covariates.csv (primary export) | 80 | -0.038637376554835894 | False |  |
| pls_1 | pmi | covariates.csv (primary export) | 80 | -0.1697929667954901 | False |  |
| pls_2 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.16570557899671823 | False |  |
| pls_2 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | -0.07273792780121895 | False |  |
| pls_2 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.1475855602437881 | False |  |
| pls_2 | doublet_score_mean | donor_technical.csv | 80 | -0.1409282700421941 | False |  |
| pls_2 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | 0.0964609497259281 | False |  |
| pls_2 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | 0.17419127988748243 | False |  |
| pls_2 | rin | covariates.csv (primary export) | 80 | 0.21174360822512153 | False |  |
| pls_2 | pmi | covariates.csv (primary export) | 80 | 0.11865460311656115 | False |  |
| pls_3 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | 0.2376933895921238 | False |  |
| pls_3 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | -0.0005860290670417253 | False |  |
| pls_3 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | -0.10989217065166432 | False |  |
| pls_3 | doublet_score_mean | donor_technical.csv | 80 | -0.13745897796530707 | False |  |
| pls_3 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | 0.12600870844513404 | False |  |
| pls_3 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.0792076887013596 | False |  |
| pls_3 | rin | covariates.csv (primary export) | 80 | -0.13151478384972812 | False |  |
| pls_3 | pmi | covariates.csv (primary export) | 80 | -0.11627471136350692 | False |  |
| all | all |  | 80 |  | True | False |


Adoption check holds in cohort_30: **False** (components pls_0, pls_1, pls_2, pls_3; flagged: [{'component': 'pls_0', 'variable': 'tmm_rna', 'spearman': -0.673230192217534}, {'component': 'pls_0', 'variable': 'tmm_atac', 'spearman': 0.5054149085794655}, {'component': 'pls_0', 'variable': 'fraction_mitochondrial_umis_mean', 'spearman': 0.7298406000937647}, {'component': 'pls_0', 'variable': 'gex_reads_mapped_confidently_to_intronic_regions_mean', 'spearman': -0.5649118900231207}, {'component': 'pls_0', 'variable': 'gex_fraction_of_transcriptomic_reads_in_cells_mean', 'spearman': -0.6321143928738866}, {'component': 'pls_0', 'variable': 'rin', 'spearman': -0.334618905873116}]).

`adoption_precedence` (comparison manifest): 

```json
{
  "rule": "if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read",
  "candidates": [
    "trend_normalized",
    "revised_candidate_trend"
  ],
  "holds": {
    "trend_normalized": {
      "cohort_30": false,
      "cohort_50": false
    },
    "revised_candidate_trend": {
      "cohort_30": false,
      "cohort_50": true
    }
  },
  "flagged": {
    "trend_normalized": {
      "cohort_30": [
        "pls_0 \u00d7 tmm_rna (-0.67)",
        "pls_0 \u00d7 tmm_atac (0.51)",
        "pls_0 \u00d7 fraction_mitochondrial_umis_mean (0.73)",
        "pls_0 \u00d7 gex_reads_mapped_confidently_to_intronic_regions_mean (-0.56)",
        "pls_0 \u00d7 gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63)",
        "pls_0 \u00d7 rin (-0.33)"
      ],
      "cohort_50": [
        "pls_0 \u00d7 tmm_rna (-0.66)",
        "pls_0 \u00d7 tmm_atac (0.58)",
        "pls_0 \u00d7 fraction_mitochondrial_umis_mean (0.78)",
        "pls_0 \u00d7 gex_reads_mapped_confidently_to_intronic_regions_mean (-0.60)",
        "pls_0 \u00d7 gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63)",
        "pls_0 \u00d7 rin (-0.34)"
      ]
    },
    "revised_candidate_trend": {
      "cohort_30": [
        "pls_1 \u00d7 tmm_rna (-0.32)"
      ],
      "cohort_50": []
    }
  },
  "candidate": null,
  "reason": "neither trend_normalized nor revised_candidate_trend meets ADOPTION_RULE in both cohorts; no revised-primary change follows and the variants stay sensitivity analyses"
}
```

#### Variant `revised_candidate_trend`

Pasted from `results/motco_sensitivity/revised_candidate_trend/comparison/cohort_30/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_30, variant `revised_candidate_trend`

the revised candidate's four settings plus the cyclic-loess trend normalization of both omics, in the primary's order: gene sets removed before the abundance filter, cyclic loess after log2-CPM, residualization after the normalization, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'}, `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `normalization` = {'rna': 'none', 'atac': 'none'}, `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 5, which stays the primary result; the `revised_candidate_trend` row (CV-selected rank 3) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| revised_candidate_trend | cv | 3 | 14.3836 | 0.166 | no | 8.49356 | 6.54021 | 22.1517 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| revised_candidate_trend | cv | 3 | 16.7741 | 0.435 | no | 16.4397 | 9.28967 | 34.0151 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| revised_candidate_trend | cv | 3 | 0.0792388 | 0.852 | no | 0.188868 | 0.104185 | 0.3839 |

###### Findings

- adoption check does not hold in cohort_30: pls_1 x tmm_rna (-0.32); the variants stay sensitivity analyses
- CV-selected rank differs: 5 in the primary analysis, 3 in the revised_candidate_trend analysis; reported beside angle, not read as better or worse
- revised_candidate_trend export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


`adoption_check.csv` (revised_candidate_trend, cohort_30):

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---|---|---|---|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.15850914205344585 | False |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.14174871073605252 | False |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.15860290670417254 | False |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 80 | 0.015822784810126583 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.2125188262977556 | False |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.12116736990154713 | False |  |
| pls_0 | rin | covariates.csv (primary export) | 80 | -0.020385436234484106 | False |  |
| pls_0 | pmi | covariates.csv (primary export) | 80 | -0.0018875003558706002 | False |  |
| pls_1 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.3225035161744023 | True |  |
| pls_1 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.1456399437412096 | False |  |
| pls_1 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.0890295358649789 | False |  |
| pls_1 | doublet_score_mean | donor_technical.csv | 80 | 0.12894983591186124 | False |  |
| pls_1 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.11486237026902739 | False |  |
| pls_1 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | -0.03900609470229723 | False |  |
| pls_1 | rin | covariates.csv (primary export) | 80 | -0.04762853790725067 | False |  |
| pls_1 | pmi | covariates.csv (primary export) | 80 | 0.08449201903577276 | False |  |
| pls_2 | tmm_rna | manifest.json (primary export, TMM factor) | 80 | -0.2342240975152368 | False |  |
| pls_2 | tmm_atac | manifest.json (primary export, TMM factor) | 80 | 0.1721050164088139 | False |  |
| pls_2 | fraction_mitochondrial_umis_mean | donor_technical.csv | 80 | 0.16727613689639007 | False |  |
| pls_2 | doublet_score_mean | donor_technical.csv | 80 | 0.017229254571026725 | False |  |
| pls_2 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 80 | -0.06265859504675718 | False |  |
| pls_2 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 80 | 0.04631973745897796 | False |  |
| pls_2 | rin | covariates.csv (primary export) | 80 | 0.056408694169256764 | False |  |
| pls_2 | pmi | covariates.csv (primary export) | 80 | 0.03531149734088353 | False |  |
| all | all |  | 80 |  | True | False |


Adoption check holds in cohort_30: **False** (components pls_0, pls_1, pls_2; flagged: [{'component': 'pls_1', 'variable': 'tmm_rna', 'spearman': -0.3225035161744023}]).

`adoption_precedence` (comparison manifest): 

```json
{
  "rule": "if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read",
  "candidates": [
    "trend_normalized",
    "revised_candidate_trend"
  ],
  "holds": {
    "trend_normalized": {
      "cohort_30": false,
      "cohort_50": false
    },
    "revised_candidate_trend": {
      "cohort_30": false,
      "cohort_50": true
    }
  },
  "flagged": {
    "trend_normalized": {
      "cohort_30": [
        "pls_0 \u00d7 tmm_rna (-0.67)",
        "pls_0 \u00d7 tmm_atac (0.51)",
        "pls_0 \u00d7 fraction_mitochondrial_umis_mean (0.73)",
        "pls_0 \u00d7 gex_reads_mapped_confidently_to_intronic_regions_mean (-0.56)",
        "pls_0 \u00d7 gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63)",
        "pls_0 \u00d7 rin (-0.33)"
      ],
      "cohort_50": [
        "pls_0 \u00d7 tmm_rna (-0.66)",
        "pls_0 \u00d7 tmm_atac (0.58)",
        "pls_0 \u00d7 fraction_mitochondrial_umis_mean (0.78)",
        "pls_0 \u00d7 gex_reads_mapped_confidently_to_intronic_regions_mean (-0.60)",
        "pls_0 \u00d7 gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63)",
        "pls_0 \u00d7 rin (-0.34)"
      ]
    },
    "revised_candidate_trend": {
      "cohort_30": [
        "pls_1 \u00d7 tmm_rna (-0.32)"
      ],
      "cohort_50": []
    }
  },
  "candidate": null,
  "reason": "neither trend_normalized nor revised_candidate_trend meets ADOPTION_RULE in both cohorts; no revised-primary change follows and the variants stay sensitivity analyses"
}
```

### 5.2 cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74) — sensitivity cohort

Superseded → floored, per variant (the superseded variants were made against the MOTCO 0.8.0 analysis; `trend_normalized` and `revised_candidate_trend` had no superseded copy on this machine):

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| primary | 1 → 2 | 0.6616 → 0.6965 | 18.7091 → 6.01569 | 0.213 → 0.596 | no → no | 0 → 23.1596 | 1 → 0.222 | no → no | 0.56471 → 0.284932 | 0.377 → 0.484 | no → no |
| keep_chrx | 1 → 5 | 0.6642 → 0.6917 | 18.9249 → 13.6639 | 0.212 → 0.207 | no → no | 0 → 24.9388 | 1 → 0.297 | no → no | 0.556064 → 0.0743536 | 0.391 → 0.78 | no → no |
| exclude_families | 1 → 5 | 0.6674 → 0.7028 | 19.0327 → 13.4939 | 0.202 → 0.211 | no → no | 0 → 24.9552 | 1 → 0.295 | no → no | 0.563235 → 0.0745555 | 0.377 → 0.78 | no → no |
| rna_trend_ratio | 1 → 4 | 0.6716 → 0.6918 | 18.0822 → 7.5918 | 0.211 → 0.562 | no → no | 0 → 25.2598 | 1 → 0.346 | no → no | 0.559086 → 0.156475 | 0.353 → 0.606 | no → no |
| atac_mean_accessibility | 1 → 7 | 0.6820 → 0.7271 | 9.1382 → 9.05688 | 0.665 → 0.455 | no → no | 0 → 23.1367 | 1 → 0.641 | no → no | 0.645472 → 0.0993176 | 0.465 → 0.84 | no → no |
| residualized | 2 → 2 | 0.4791 → 0.4942 | 8.94403 → 8.94403 | 0.444 → 0.444 | no → no | 30.0936 → 30.0936 | 0.104 → 0.104 | no → no | 0.317151 → 0.317151 | 0.338 → 0.338 | no → no |
| revised_candidate | 2 → 2 | 0.5514 → 0.5514 | 6.33891 → 6.33891 | 0.585 → 0.585 | no → no | 2.15631 → 2.15631 | 0.902 → 0.902 | no → no | 0.230511 → 0.230511 | 0.564 → 0.564 | no → no |
| trend_normalized | n/a → 4 | n/a → 0.6722 | n/a → 9.24514 | n/a → 0.529 | n/a → no | n/a → 25.992 | n/a → 0.296 | n/a → no | n/a → 0.116966 | n/a → 0.709 | n/a → no |
| revised_candidate_trend | n/a → 2 | n/a → 0.5588 | n/a → 5.41504 | n/a → 0.585 | n/a → no | n/a → 19.9809 | n/a → 0.22 | n/a → no | n/a → 0.259816 | n/a → 0.459 | n/a → no |


#### Variant `keep_chrx`

Pasted from `results/motco_sensitivity/keep_chrx/comparison/cohort_50/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_50, variant `keep_chrx`

chrX kept; chrY and chrM removed; every other export setting equal to the primary. Export setting(s) changed: `excluded_chromosomes` = ['chrY', 'chrM'] (primary: `excluded_chromosomes` = ['chrX', 'chrY', 'chrM']); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `keep_chrx` row (CV-selected rank 5) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| keep_chrx | cv | 5 | 13.6639 | 0.207 | no | 8.55689 | 6.28015 | 20.7816 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| keep_chrx | cv | 5 | 24.9388 | 0.297 | no | 21.3378 | 9.62679 | 40.4569 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| keep_chrx | cv | 5 | 0.0743536 | 0.78 | no | 0.135753 | 0.0737208 | 0.271685 |

###### Findings

- CV-selected rank differs: 2 in the primary analysis, 5 in the keep_chrx analysis; reported beside angle, not read as better or worse
- keep_chrx export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


#### Variant `exclude_families`

Pasted from `results/motco_sensitivity/exclude_families/comparison/cohort_50/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_50, variant `exclude_families`

RNA genes of the ribosomal, mitochondrial and hemoglobin sets (TECHNICAL_GENE_SETS) removed before the abundance filter; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'] (primary: `excluded_gene_sets` = []); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `exclude_families` row (CV-selected rank 5) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| exclude_families | cv | 5 | 13.4939 | 0.211 | no | 8.5631 | 6.27947 | 20.453 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| exclude_families | cv | 5 | 24.9552 | 0.295 | no | 21.2822 | 9.57813 | 40.4704 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| exclude_families | cv | 5 | 0.0745555 | 0.78 | no | 0.135411 | 0.0734519 | 0.270234 |

###### Findings

- CV-selected rank differs: 2 in the primary analysis, 5 in the exclude_families analysis; reported beside angle, not read as better or worse
- exclude_families export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


#### Variant `rna_trend_ratio`

Pasted from `results/motco_sensitivity/rna_trend_ratio/comparison/cohort_50/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_50, variant `rna_trend_ratio`

RNA features ranked at the selection step by the ratio of their log2-CPM variance to a binned-median mean–variance trend instead of by raw variance; the ATAC block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'trend_ratio', 'atac': 'variance'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `rna_trend_ratio` row (CV-selected rank 4) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| rna_trend_ratio | cv | 4 | 7.5918 | 0.562 | no | 10.2648 | 7.66266 | 25.1727 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| rna_trend_ratio | cv | 4 | 25.2598 | 0.346 | no | 22.2226 | 10.8818 | 42.3442 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| rna_trend_ratio | cv | 4 | 0.156475 | 0.606 | no | 0.192497 | 0.100341 | 0.382166 |

###### Findings

- CV-selected rank differs: 2 in the primary analysis, 4 in the rna_trend_ratio analysis; reported beside angle, not read as better or worse
- rna_trend_ratio export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


#### Variant `atac_mean_accessibility`

Pasted from `results/motco_sensitivity/atac_mean_accessibility/comparison/cohort_50/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_50, variant `atac_mean_accessibility`

ATAC peaks ranked at the selection step by mean log2-CPM (total accessibility) instead of by variance; the RNA block and every other export setting equal to the primary. Export setting(s) changed: `selection_rule` = {'rna': 'variance', 'atac': 'mean'} (primary: `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `atac_mean_accessibility` row (CV-selected rank 7) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| atac_mean_accessibility | cv | 7 | 9.05688 | 0.455 | no | 10.5116 | 8.9503 | 27.7114 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| atac_mean_accessibility | cv | 7 | 23.1367 | 0.641 | no | 30.3072 | 14.1015 | 58.6036 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| atac_mean_accessibility | cv | 7 | 0.0993176 | 0.84 | no | 0.210347 | 0.108098 | 0.401405 |

###### Findings

- CV-selected rank differs: 2 in the primary analysis, 7 in the atac_mean_accessibility analysis; reported beside angle, not read as better or worse
- atac_mean_accessibility export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


#### Variant `residualized`

Pasted from `results/motco_sensitivity/residualized/comparison/cohort_50/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_50, variant `residualized`

after log2-CPM and before the selection step, each abundance-passing feature replaced by its mean plus the residual of its regression on an intercept and RESIDUALIZATION_COVARIATES (the export's log10 library sizes and the technical table's gene-set fractions); the selection ranks the residual variance; every other export setting equal to the primary. Export setting(s) changed: `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'] (primary: `residualize` = None); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `residualized` row (CV-selected rank 2) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| residualized | cv | 2 | 8.94403 | 0.444 | no | 9.36668 | 7.40074 | 23.7553 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| residualized | cv | 2 | 30.0936 | 0.104 | no | 14.4779 | 11.3133 | 36.2531 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| residualized | cv | 2 | 0.317151 | 0.338 | no | 0.264388 | 0.142759 | 0.517872 |

###### Findings

- residualized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


#### Variant `revised_candidate`

Pasted from `results/motco_sensitivity/revised_candidate/comparison/cohort_50/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_50, variant `revised_candidate`

the four settings together in the primary's order: gene sets removed before the abundance filter, residualization after log2-CPM, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `revised_candidate` row (CV-selected rank 2) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected nor residualized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| revised_candidate | cv | 2 | 6.33891 | 0.585 | no | 9.67538 | 7.57586 | 25.2956 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| revised_candidate | cv | 2 | 2.15631 | 0.902 | no | 12.4302 | 9.75934 | 30.5331 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| revised_candidate | cv | 2 | 0.230511 | 0.564 | no | 0.276132 | 0.155841 | 0.565357 |

###### Findings

- adoption check does not hold in cohort_50: pls_0 x tmm_atac (0.37), pls_0 x gex_reads_mapped_confidently_to_intronic_regions_mean (-0.33), pls_1 x tmm_atac (-0.37); the variants stay sensitivity analyses
- revised_candidate export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


`adoption_check.csv` (revised_candidate, cohort_50):

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---|---|---|---|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | -0.2872269529803776 | False |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | 0.37469085523880036 | True |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | 0.20393928174750087 | False |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 74 | 0.03386893743058126 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | -0.32545224326877603 | True |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | -0.19686042206590149 | False |  |
| pls_0 | rin | covariates.csv (primary export) | 74 | -0.11312147045596636 | False |  |
| pls_0 | pmi | covariates.csv (primary export) | 74 | -0.0145176589153464 | False |  |
| pls_1 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | 0.20873750462791554 | False |  |
| pls_1 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | -0.37131432802665676 | True |  |
| pls_1 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | -0.029692706405035168 | False |  |
| pls_1 | doublet_score_mean | donor_technical.csv | 74 | -0.1146094039244724 | False |  |
| pls_1 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | 0.13405504669043325 | False |  |
| pls_1 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | 0.1605479452054794 | False |  |
| pls_1 | rin | covariates.csv (primary export) | 74 | 0.202069281126135 | False |  |
| pls_1 | pmi | covariates.csv (primary export) | 74 | 0.06529983724372136 | False |  |
| all | all |  | 74 |  | True | False |


Adoption check holds in cohort_50: **False** (components pls_0, pls_1; flagged: [{'component': 'pls_0', 'variable': 'tmm_atac', 'spearman': 0.37469085523880036}, {'component': 'pls_0', 'variable': 'gex_reads_mapped_confidently_to_intronic_regions_mean', 'spearman': -0.32545224326877603}, {'component': 'pls_1', 'variable': 'tmm_atac', 'spearman': -0.37131432802665676}]).

#### Variant `trend_normalized`

Pasted from `results/motco_sensitivity/trend_normalized/comparison/cohort_50/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_50, variant `trend_normalized`

after log2-CPM and before the selection step, each donor's log2-CPM over the abundance-passing features corrected by limma's fast cyclic loess (a lowess curve of its deviation from the per-feature mean over donors against that mean, three iterations, adaptive span) in both omics; every other export setting equal to the primary; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'} (primary: `normalization` = {'rna': 'none', 'atac': 'none'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `trend_normalized` row (CV-selected rank 4) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| trend_normalized | cv | 4 | 9.24514 | 0.529 | no | 11.7307 | 8.93463 | 28.563 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| trend_normalized | cv | 4 | 25.992 | 0.296 | no | 21.5415 | 10.8058 | 41.6538 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| trend_normalized | cv | 4 | 0.116966 | 0.709 | no | 0.174987 | 0.0897564 | 0.339846 |

###### Findings

- adoption check does not hold in cohort_50: pls_0 x tmm_rna (-0.66), pls_0 x tmm_atac (0.58), pls_0 x fraction_mitochondrial_umis_mean (0.78), pls_0 x gex_reads_mapped_confidently_to_intronic_regions_mean (-0.60), pls_0 x gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63), pls_0 x rin (-0.34); the variants stay sensitivity analyses
- CV-selected rank differs: 2 in the primary analysis, 4 in the trend_normalized analysis; reported beside angle, not read as better or worse
- trend_normalized export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


`adoption_check.csv` (trend_normalized, cohort_50):

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---|---|---|---|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | -0.6610144390966308 | True |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | 0.5777267678637541 | True |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | 0.7825101814142909 | True |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 74 | 0.010914476119955569 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | -0.5957837525802175 | True |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | -0.6339726027397259 | True |  |
| pls_0 | rin | covariates.csv (primary export) | 74 | -0.3407567667722281 | True |  |
| pls_0 | pmi | covariates.csv (primary export) | 74 | 0.07235127157403248 | False |  |
| pls_1 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | 0.12293224731580893 | False |  |
| pls_1 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | 0.03306923361717881 | False |  |
| pls_1 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | -0.08244353942984077 | False |  |
| pls_1 | doublet_score_mean | donor_technical.csv | 74 | -0.02305812661977045 | False |  |
| pls_1 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | -0.13949011100057343 | False |  |
| pls_1 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | 0.011447611995557199 | False |  |
| pls_1 | rin | covariates.csv (primary export) | 74 | -0.06908452772117678 | False |  |
| pls_1 | pmi | covariates.csv (primary export) | 74 | -0.05546338263169074 | False |  |
| pls_2 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | -0.14458348759718617 | False |  |
| pls_2 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | -0.1708848574601999 | False |  |
| pls_2 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | 0.04678267308404294 | False |  |
| pls_2 | doublet_score_mean | donor_technical.csv | 74 | -0.1978970751573491 | False |  |
| pls_2 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | 0.1641182089508817 | False |  |
| pls_2 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | 0.23631247686042203 | False |  |
| pls_2 | rin | covariates.csv (primary export) | 74 | 0.21655570224777113 | False |  |
| pls_2 | pmi | covariates.csv (primary export) | 74 | 0.07683989468765488 | False |  |
| pls_3 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | -0.19330618289522394 | False |  |
| pls_3 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | -0.07160310995927432 | False |  |
| pls_3 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | 0.05806738245094408 | False |  |
| pls_3 | doublet_score_mean | donor_technical.csv | 74 | 0.07915586819696406 | False |  |
| pls_3 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | 0.0008589474931556703 | False |  |
| pls_3 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | 0.16972972972972972 | False |  |
| pls_3 | rin | covariates.csv (primary export) | 74 | 0.16500892770452602 | False |  |
| pls_3 | pmi | covariates.csv (primary export) | 74 | 0.09482401501748194 | False |  |
| all | all |  | 74 |  | True | False |


Adoption check holds in cohort_50: **False** (components pls_0, pls_1, pls_2, pls_3; flagged: [{'component': 'pls_0', 'variable': 'tmm_rna', 'spearman': -0.6610144390966308}, {'component': 'pls_0', 'variable': 'tmm_atac', 'spearman': 0.5777267678637541}, {'component': 'pls_0', 'variable': 'fraction_mitochondrial_umis_mean', 'spearman': 0.7825101814142909}, {'component': 'pls_0', 'variable': 'gex_reads_mapped_confidently_to_intronic_regions_mean', 'spearman': -0.5957837525802175}, {'component': 'pls_0', 'variable': 'gex_fraction_of_transcriptomic_reads_in_cells_mean', 'spearman': -0.6339726027397259}, {'component': 'pls_0', 'variable': 'rin', 'spearman': -0.3407567667722281}]).

`adoption_precedence` (comparison manifest): 

```json
{
  "rule": "if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read",
  "candidates": [
    "trend_normalized",
    "revised_candidate_trend"
  ],
  "holds": {
    "trend_normalized": {
      "cohort_30": false,
      "cohort_50": false
    },
    "revised_candidate_trend": {
      "cohort_30": false,
      "cohort_50": true
    }
  },
  "flagged": {
    "trend_normalized": {
      "cohort_30": [
        "pls_0 \u00d7 tmm_rna (-0.67)",
        "pls_0 \u00d7 tmm_atac (0.51)",
        "pls_0 \u00d7 fraction_mitochondrial_umis_mean (0.73)",
        "pls_0 \u00d7 gex_reads_mapped_confidently_to_intronic_regions_mean (-0.56)",
        "pls_0 \u00d7 gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63)",
        "pls_0 \u00d7 rin (-0.33)"
      ],
      "cohort_50": [
        "pls_0 \u00d7 tmm_rna (-0.66)",
        "pls_0 \u00d7 tmm_atac (0.58)",
        "pls_0 \u00d7 fraction_mitochondrial_umis_mean (0.78)",
        "pls_0 \u00d7 gex_reads_mapped_confidently_to_intronic_regions_mean (-0.60)",
        "pls_0 \u00d7 gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63)",
        "pls_0 \u00d7 rin (-0.34)"
      ]
    },
    "revised_candidate_trend": {
      "cohort_30": [
        "pls_1 \u00d7 tmm_rna (-0.32)"
      ],
      "cohort_50": []
    }
  },
  "candidate": null,
  "reason": "neither trend_normalized nor revised_candidate_trend meets ADOPTION_RULE in both cohorts; no revised-primary change follows and the variants stay sensitivity analyses"
}
```

#### Variant `revised_candidate_trend`

Pasted from `results/motco_sensitivity/revised_candidate_trend/comparison/cohort_50/sensitivity.md` (statistics and findings; the export and adoption-check sections are in the file, the adoption table is pasted below from its CSV):

##### MOTCO preprocessing sensitivity: cohort_50, variant `revised_candidate_trend`

the revised candidate's four settings plus the cyclic-loess trend normalization of both omics, in the primary's order: gene sets removed before the abundance filter, cyclic loess after log2-CPM, residualization after the normalization, trend-ratio ranking for RNA and mean-accessibility ranking for ATAC at the selection step; the adoption rule ADOPTION_RULE is evaluated for this variant. Export setting(s) changed: `excluded_gene_sets` = ['ribosomal', 'mitochondrial', 'hemoglobin'], `normalization` = {'rna': 'cyclic_loess', 'atac': 'cyclic_loess'}, `residualize` = ['log10_lib_rna', 'log10_lib_atac', 'ribosomal_pseudobulk', 'mitochondrial_pseudobulk'], `selection_rule` = {'rna': 'trend_ratio', 'atac': 'mean'} (primary: `excluded_gene_sets` = [], `normalization` = {'rna': 'none', 'atac': 'none'}, `residualize` = None, `selection_rule` = {'rna': 'variance', 'atac': 'variance'}); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 2, which stays the primary result; the `revised_candidate_trend` row (CV-selected rank 2) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse. MOTCO's operating characteristics were measured on none of the variant blocks (neither chrX-kept, family-excluded, trend- or mean-selected, residualized nor trend-normalized), so the variant carries no Type I or power claim of its own.

###### delta

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| revised_candidate_trend | cv | 2 | 5.41504 | 0.585 | no | 8.1077 | 6.27099 | 20.5446 |

###### angle (degrees)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| revised_candidate_trend | cv | 2 | 19.9809 | 0.22 | no | 12.7642 | 9.78384 | 31.2805 |

###### shape

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| revised_candidate_trend | cv | 2 | 0.259816 | 0.459 | no | 0.256019 | 0.130477 | 0.493029 |

###### Findings

- adoption check holds in cohort_50: no candidate component reaches |Spearman| >= 0.3 against any of 8 check variables (the rule needs both cohorts; a revised primary is a later change)
- revised_candidate_trend export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')


`adoption_check.csv` (revised_candidate_trend, cohort_50):

| component | variable | source | n | spearman | flag | holds |
|---|---|---|---|---|---|---|
| pls_0 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | -0.2094187338022954 | False |  |
| pls_0 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | 0.2332913735653461 | False |  |
| pls_0 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | 0.20592373195112915 | False |  |
| pls_0 | doublet_score_mean | donor_technical.csv | 74 | 0.026019992595335054 | False |  |
| pls_0 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | -0.2753667187540782 | False |  |
| pls_0 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | -0.1572602739726027 | False |  |
| pls_0 | rin | covariates.csv (primary export) | 74 | -0.03738622383538812 | False |  |
| pls_0 | pmi | covariates.csv (primary export) | 74 | 0.035153474087874496 | False |  |
| pls_1 | tmm_rna | manifest.json (primary export, TMM factor) | 74 | 0.21969640873750457 | False |  |
| pls_1 | tmm_atac | manifest.json (primary export, TMM factor) | 74 | -0.1584154017030729 | False |  |
| pls_1 | fraction_mitochondrial_umis_mean | donor_technical.csv | 74 | -0.07041836356904849 | False |  |
| pls_1 | doublet_score_mean | donor_technical.csv | 74 | -0.1317882265827471 | False |  |
| pls_1 | gex_reads_mapped_confidently_to_intronic_regions_mean | donor_technical.csv | 74 | 0.11865322957177983 | False |  |
| pls_1 | gex_fraction_of_transcriptomic_reads_in_cells_mean | donor_technical.csv | 74 | 0.09870418363569047 | False |  |
| pls_1 | rin | covariates.csv (primary export) | 74 | 0.13326618694413112 | False |  |
| pls_1 | pmi | covariates.csv (primary export) | 74 | -0.06279628177770755 | False |  |
| all | all |  | 74 |  | False | True |


Adoption check holds in cohort_50: **True** (components pls_0, pls_1; flagged: none).

`adoption_precedence` (comparison manifest): 

```json
{
  "rule": "if trend_normalized meets ADOPTION_RULE in both cohorts it is the candidate a later revised-primary change proposes, whatever revised_candidate_trend gives (no gene removed, no degrees of freedom spent, no technical covariate read); if only revised_candidate_trend meets it in both cohorts, that one is; if neither does, the variants stay sensitivity analyses. Only the holds values of the adoption tables are read",
  "candidates": [
    "trend_normalized",
    "revised_candidate_trend"
  ],
  "holds": {
    "trend_normalized": {
      "cohort_30": false,
      "cohort_50": false
    },
    "revised_candidate_trend": {
      "cohort_30": false,
      "cohort_50": true
    }
  },
  "flagged": {
    "trend_normalized": {
      "cohort_30": [
        "pls_0 \u00d7 tmm_rna (-0.67)",
        "pls_0 \u00d7 tmm_atac (0.51)",
        "pls_0 \u00d7 fraction_mitochondrial_umis_mean (0.73)",
        "pls_0 \u00d7 gex_reads_mapped_confidently_to_intronic_regions_mean (-0.56)",
        "pls_0 \u00d7 gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63)",
        "pls_0 \u00d7 rin (-0.33)"
      ],
      "cohort_50": [
        "pls_0 \u00d7 tmm_rna (-0.66)",
        "pls_0 \u00d7 tmm_atac (0.58)",
        "pls_0 \u00d7 fraction_mitochondrial_umis_mean (0.78)",
        "pls_0 \u00d7 gex_reads_mapped_confidently_to_intronic_regions_mean (-0.60)",
        "pls_0 \u00d7 gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.63)",
        "pls_0 \u00d7 rin (-0.34)"
      ]
    },
    "revised_candidate_trend": {
      "cohort_30": [
        "pls_1 \u00d7 tmm_rna (-0.32)"
      ],
      "cohort_50": []
    }
  },
  "candidate": null,
  "reason": "neither trend_normalized nor revised_candidate_trend meets ADOPTION_RULE in both cohorts; no revised-primary change follows and the variants stay sensitivity analyses"
}
```

## 6. Latent-adjustment sensitivity

Both arms were rerun for both cell types on the unchanged primary exports with the eight covariates and the floor; `mediator_check.csv` was written before the arms of each cohort and is byte-identical to the superseded one. Each `extended_covariates` arm reproduces its primary's CV-selected rank and CV mean AUROC (astrocytes 5 / 0.6356 and 2 / 0.6965; microglia 2 / 0.6852 and 4 / 0.6565), as MOTCO's contract says (the covariate-blind basis does not see the covariates). The `latent_adjusted` arms select ranks 4 and 3 (astrocytes) and 3 and 3 (microglia); where the superseded arm already selected rank 3 (astrocyte and microglia `cohort_50`) the floored statistics are identical and only the microglia CV AUROC moves (0.6205 → 0.6261), because the outer folds now select at ≥ 2. No arm rejects any statistic; a CV-selected rank that differs between arms is reported beside `angle`. MOTCO's power comparison (`candidate_fails`) still excludes the latent adjustment as a replacement, so the arms are sensitivity only.

### 6.1 microglia, cohort_30 (≥ 30 microglia nuclei in both omics, n = 76) — primary cohort

Pasted from `results/microglia/motco_latent_sensitivity/comparison/cohort_30/latent_sensitivity.md`:

#### MOTCO latent-adjustment sensitivity: cohort_30

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 2, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

##### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 5.10924 | 0.689 | no | 10.2655 | 7.68786 | 25.6788 |
| extended_covariates | 8 | none | 2 | 6.89041 | 0.546 | no | 8.81982 | 6.49261 | 20.9763 |
| latent_adjusted | 8 | stage | 3 | 2.43939 | 0.715 | no | 5.00415 | 3.76187 | 12.0786 |

##### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 12.5272 | 0.424 | no | 12.6966 | 9.51188 | 30.1659 |
| extended_covariates | 8 | none | 2 | 5.24863 | 0.677 | no | 10.3465 | 7.85108 | 24.7464 |
| latent_adjusted | 8 | stage | 3 | 5.99532 | 0.849 | no | 15.5632 | 9.78946 | 33.8817 |

##### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 0.0242177 | 0.99 | no | 0.259273 | 0.142347 | 0.515887 |
| extended_covariates | 8 | none | 2 | 0.192379 | 0.509 | no | 0.209276 | 0.108795 | 0.419364 |
| latent_adjusted | 8 | stage | 3 | 0.130018 | 0.294 | no | 0.104228 | 0.0554984 | 0.198081 |

##### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 2 | 0.6852 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 111.7 | `results/microglia/motco_analysis/cohort_30` |
| extended_covariates | cv | 2 | 0.6852 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 164.6 | `extended_covariates/motco_analysis/cohort_30` |
| latent_adjusted | cv | 3 | 0.6575 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 189.5 | `latent_adjusted/motco_analysis/cohort_30` |

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


`mediator_check.csv`:

| cohort | covariate | subset | n | spearman | flagged |
|---|---|---|---|---|---|
| cohort_30 | ribosomal_pseudobulk | all | 76 | 0.053354167238899325 | False |
| cohort_30 | ribosomal_pseudobulk | F | 45 | -0.03006050866078459 | False |
| cohort_30 | ribosomal_pseudobulk | M | 31 | 0.12715085021216893 | False |
| cohort_30 | mitochondrial_pseudobulk | all | 76 | 0.4978299890219748 | True |
| cohort_30 | mitochondrial_pseudobulk | F | 45 | 0.5768858485747186 | True |
| cohort_30 | mitochondrial_pseudobulk | M | 31 | 0.3292222350611378 | True |


Superseded → floored:

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| microglia cohort_30 primary | 1 → 2 | 0.6921 → 0.6852 | 8.0642 → 5.10924 | 0.391 → 0.689 | no → no | 0 → 12.5272 | 1 → 0.424 | no → no | 0.154684 → 0.0242177 | 0.637 → 0.99 | no → no |
| microglia cohort_30 extended_covariates | 1 → 2 | 0.6921 → 0.6852 | 3.00017 → 6.89041 | 0.587 → 0.546 | no → no | 0 → 5.24863 | 1 → 0.677 | no → no | 0.190112 → 0.192379 | 0.611 → 0.509 | no → no |
| microglia cohort_30 latent_adjusted | 1 → 3 | 0.6415 → 0.6575 | 4.19879 → 2.43939 | 0.35 → 0.715 | no → no | 0 → 5.99532 | 1 → 0.849 | no → no | 0.170999 → 0.130018 | 0.274 → 0.294 | no → no |


### 6.2 microglia, cohort_50 (≥ 50 microglia nuclei in both omics, n = 63) — sensitivity cohort

Pasted from `results/microglia/motco_latent_sensitivity/comparison/cohort_50/latent_sensitivity.md`:

#### MOTCO latent-adjustment sensitivity: cohort_50

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 4, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

##### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 4 | 1.6686 | 0.819 | no | 5.78012 | 4.4338 | 14.0167 |
| extended_covariates | 8 | none | 4 | 1.55793 | 0.825 | no | 5.48154 | 4.15429 | 13.439 |
| latent_adjusted | 8 | stage | 3 | 1.11008 | 0.875 | no | 5.86056 | 4.59823 | 14.6532 |

##### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 4 | 22.2604 | 0.227 | no | 16.6216 | 7.94278 | 31.3913 |
| extended_covariates | 8 | none | 4 | 10.8156 | 0.616 | no | 13.4235 | 6.09246 | 24.5714 |
| latent_adjusted | 8 | stage | 3 | 12.1741 | 0.611 | no | 16.8349 | 10.5021 | 36.7317 |

##### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 4 | 0.0797433 | 0.598 | no | 0.0979335 | 0.0515193 | 0.196376 |
| extended_covariates | 8 | none | 4 | 0.0750599 | 0.537 | no | 0.08453 | 0.04461 | 0.165572 |
| latent_adjusted | 8 | stage | 3 | 0.176571 | 0.114 | no | 0.108847 | 0.0564956 | 0.216514 |

##### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 4 | 0.6565 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 107.4 | `results/microglia/motco_analysis/cohort_50` |
| extended_covariates | cv | 4 | 0.6565 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 149.3 | `extended_covariates/motco_analysis/cohort_50` |
| latent_adjusted | cv | 3 | 0.6261 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 165.4 | `latent_adjusted/motco_analysis/cohort_50` |

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


`mediator_check.csv`:

| cohort | covariate | subset | n | spearman | flagged |
|---|---|---|---|---|---|
| cohort_50 | ribosomal_pseudobulk | all | 63 | -0.02905110889935568 | False |
| cohort_50 | ribosomal_pseudobulk | F | 36 | -0.22648732511857167 | False |
| cohort_50 | ribosomal_pseudobulk | M | 27 | 0.1584580826700018 | False |
| cohort_50 | mitochondrial_pseudobulk | all | 63 | 0.43969662434939916 | True |
| cohort_50 | mitochondrial_pseudobulk | F | 36 | 0.5412710703019331 | True |
| cohort_50 | mitochondrial_pseudobulk | M | 27 | 0.25405140452613784 | False |


Superseded → floored:

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| microglia cohort_50 primary | 1 → 4 | 0.6414 → 0.6565 | 12.4641 → 1.6686 | 0.227 → 0.819 | no → no | 0 → 22.2604 | 1 → 0.227 | no → no | 0.269568 → 0.0797433 | 0.399 → 0.598 | no → no |
| microglia cohort_50 extended_covariates | 1 → 4 | 0.6414 → 0.6565 | 2.46906 → 1.55793 | 0.699 → 0.825 | no → no | 0 → 10.8156 | 1 → 0.616 | no → no | 0.270584 → 0.0750599 | 0.351 → 0.537 | no → no |
| microglia cohort_50 latent_adjusted | 3 → 3 | 0.6205 → 0.6261 | 1.11008 → 1.11008 | 0.875 → 0.875 | no → no | 12.1741 → 12.1741 | 0.611 → 0.611 | no → no | 0.176571 → 0.176571 | 0.114 → 0.114 | no → no |


### 6.3 astrocyte, cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

Pasted from `results/motco_latent_sensitivity/comparison/cohort_30/latent_sensitivity.md`:

#### MOTCO latent-adjustment sensitivity: cohort_30

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 5, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

##### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| extended_covariates | 8 | none | 5 | 0.492704 | 0.942 | no | 5.93102 | 4.6327 | 15.0437 |
| latent_adjusted | 8 | stage | 4 | 6.88452 | 0.162 | no | 3.95376 | 3.02878 | 9.75925 |

##### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| extended_covariates | 8 | none | 5 | 21.5237 | 0.286 | no | 17.7736 | 6.99856 | 30.0125 |
| latent_adjusted | 8 | stage | 4 | 14.0475 | 0.673 | no | 19.0842 | 8.9644 | 35.5561 |

##### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| extended_covariates | 8 | none | 5 | 0.156527 | 0.137 | no | 0.0990181 | 0.0518795 | 0.191642 |
| latent_adjusted | 8 | stage | 4 | 0.105552 | 0.256 | no | 0.0794449 | 0.0419201 | 0.152948 |

##### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 5 | 0.6356 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 115.2 | `results/motco_analysis/cohort_30` |
| extended_covariates | cv | 5 | 0.6356 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 176.5 | `extended_covariates/motco_analysis/cohort_30` |
| latent_adjusted | cv | 4 | 0.6433 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 198.0 | `latent_adjusted/motco_analysis/cohort_30` |

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


`mediator_check.csv`:

| cohort | covariate | subset | n | spearman | flagged |
|---|---|---|---|---|---|
| cohort_30 | ribosomal_pseudobulk | all | 80 | -0.39934905947404925 | True |
| cohort_30 | ribosomal_pseudobulk | F | 49 | -0.4317995185367991 | True |
| cohort_30 | ribosomal_pseudobulk | M | 31 | -0.2850545728135508 | False |
| cohort_30 | mitochondrial_pseudobulk | all | 80 | 0.31182399404299566 | True |
| cohort_30 | mitochondrial_pseudobulk | F | 49 | 0.34218723283064734 | True |
| cohort_30 | mitochondrial_pseudobulk | M | 31 | 0.18925051713862057 | False |


Superseded → floored:

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| astrocyte cohort_30 primary | 1 → 5 | 0.6313 → 0.6356 | 14.941 → 5.95773 | 0.3 → 0.468 | no → no | 0 → 21.9941 | 1 → 0.396 | no → no | 0.608832 → 0.15645 | 0.333 → 0.284 | no → no |
| astrocyte cohort_30 extended_covariates | 1 → 5 | 0.6313 → 0.6356 | 3.93908 → 0.492704 | 0.641 → 0.942 | no → no | 0 → 21.5237 | 1 → 0.286 | no → no | 1.01735 → 0.156527 | 0.127 → 0.137 | no → no |
| astrocyte cohort_30 latent_adjusted | 1 → 4 | 0.6299 → 0.6433 | 1.3709 → 6.88452 | 0.717 → 0.162 | no → no | 0 → 14.0475 | 1 → 0.673 | no → no | 0.16446 → 0.105552 | 0.327 → 0.256 | no → no |


### 6.4 astrocyte, cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74) — sensitivity cohort

Pasted from `results/motco_latent_sensitivity/comparison/cohort_50/latent_sensitivity.md`:

#### MOTCO latent-adjustment sensitivity: cohort_50

Covariate-set and latent-basis sensitivity of the covariate-adjusted trajectory analysis. The `primary` row is the primary analysis (CV-selected rank 2, the six analysis covariates, covariate-blind basis), read from its `result.json` and never rerun; it stays the primary result, and MOTCO's power comparison (`candidate_fails`) excludes the latent adjustment as its replacement. `extended_covariates` adds `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` (raw fractions from `donor_technical.csv`) with the covariate-blind basis; `latent_adjusted` enters the same eight covariates and fits the basis on features residualized on them conditional on `stage`, with the trajectory model, LS-means and RRPP of `extended_covariates`, so those two rows differ through the basis and the CV-selected rank only. Both arms are read as sensitivity, never as the result: a rejection status that differs between `primary` and `extended_covariates` is covariate-set dependence, one that differs between `extended_covariates` and `latent_adjusted` is basis dependence; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in every row. 999 RRPP permutations with seed 20261019 in every row (the same permutation indices); α = 0.05 (rejected when p ≤ α). MOTCO measured the latent adjustment's operating characteristics at its Phase 6 design only (Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one), on neither these blocks nor this eight-covariate set, so neither arm carries a Type I or power claim of its own.

##### delta

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| extended_covariates | 8 | none | 2 | 0.416021 | 0.961 | no | 7.3652 | 5.91377 | 18.8894 |
| latent_adjusted | 8 | stage | 3 | 4.96136 | 0.494 | no | 5.99914 | 4.52218 | 14.382 |

##### angle (degrees)

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| extended_covariates | 8 | none | 2 | 23.6532 | 0.075 | no | 10.8714 | 7.82222 | 25.8055 |
| latent_adjusted | 8 | stage | 3 | 24.7286 | 0.154 | no | 15.5719 | 9.05433 | 32.1693 |

##### shape

| arm | covariates | conditioning | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | 6 | none | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| extended_covariates | 8 | none | 2 | 0.151547 | 0.764 | no | 0.264643 | 0.13837 | 0.518351 |
| latent_adjusted | 8 | stage | 3 | 0.132168 | 0.337 | no | 0.110548 | 0.0582071 | 0.219488 |

##### Runs

| arm | selection | rank | cv_mean_auroc | covariates (in order) | conditioning | runtime (s) | directory |
|---|---|---:|---:|---|---|---:|---|
| primary | cv | 2 | 0.6965 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` | none | 114.4 | `results/motco_analysis/cohort_50` |
| extended_covariates | cv | 2 | 0.6965 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | none | 158.4 | `extended_covariates/motco_analysis/cohort_50` |
| latent_adjusted | cv | 3 | 0.6645 | `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`, `ribosomal_pseudobulk`, `mitochondrial_pseudobulk` | stage | 168.5 | `latent_adjusted/motco_analysis/cohort_50` |

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


`mediator_check.csv`:

| cohort | covariate | subset | n | spearman | flagged |
|---|---|---|---|---|---|
| cohort_50 | ribosomal_pseudobulk | all | 74 | -0.351804210204838 | True |
| cohort_50 | ribosomal_pseudobulk | F | 44 | -0.36026658902869707 | True |
| cohort_50 | ribosomal_pseudobulk | M | 30 | -0.28623379151947786 | False |
| cohort_50 | mitochondrial_pseudobulk | all | 74 | 0.3560280532565137 | True |
| cohort_50 | mitochondrial_pseudobulk | F | 44 | 0.4224032985543566 | True |
| cohort_50 | mitochondrial_pseudobulk | M | 30 | 0.18892850755379928 | False |


Superseded → floored:

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| astrocyte cohort_50 primary | 1 → 2 | 0.6616 → 0.6965 | 18.7091 → 6.01569 | 0.213 → 0.596 | no → no | 0 → 23.1596 | 1 → 0.222 | no → no | 0.56471 → 0.284932 | 0.377 → 0.484 | no → no |
| astrocyte cohort_50 extended_covariates | 1 → 2 | 0.6616 → 0.6965 | 5.61916 → 0.416021 | 0.585 → 0.961 | no → no | 0 → 23.6532 | 1 → 0.075 | no → no | 0.697892 → 0.151547 | 0.201 → 0.764 | no → no |
| astrocyte cohort_50 latent_adjusted | 3 → 3 | 0.6645 → 0.6645 | 4.96136 → 4.96136 | 0.494 → 0.494 | no → no | 24.7286 → 24.7286 | 0.154 → 0.154 | no → no | 0.132168 → 0.132168 | 0.337 → 0.337 | no → no |


## 7. Stage DE (DESeq2 and limma-voom)

The stage-DE fits do not depend on MOTCO; they were refitted so that every output rests on one pin, and their concordance cells now read the floored primaries' `delta` status (not rejected in any). Every concordance cell, DESeq2 (`factor_interaction` union, exported) and limma (`interaction_joint`, exported), in both omics, both cohorts and both cell types, and on the exclusion roots, is concordant with 0 exported interaction features, and the methods agree. The microglia DESeq2 tables, all limma tables and the exclusion roots' tables are byte-identical to their superseded copies (made on this machine); the astrocyte DESeq2 tables, whose superseded copies were fitted on `nixos-dell-laptop` on 2026-10-21, differ beyond rounding with every call equal (section 11). The stage counts that E1 reads are unchanged.

### 7.1 Concordance cells (from `stage_de_summary.csv` and `method_comparison.csv`)

| cell type | cohort | cohorts' root | method | omic | cell | n significant | delta rejected | concordance |
|---|---|---|---|---|---|---|---|---|
| microglia | cohort_30 | primary | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant |
| microglia | cohort_30 | primary | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant |
| microglia | cohort_30 | primary | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| microglia | cohort_30 | primary | limma | atac | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant |
| microglia | cohort_30 | excluding_severely_affected | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| microglia | cohort_50 | primary | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant |
| microglia | cohort_50 | primary | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant |
| microglia | cohort_50 | primary | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| microglia | cohort_50 | primary | limma | atac | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant |
| microglia | cohort_50 | excluding_severely_affected | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| astrocyte | cohort_30 | primary | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant |
| astrocyte | cohort_30 | primary | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant |
| astrocyte | cohort_30 | primary | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| astrocyte | cohort_30 | primary | limma | atac | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant |
| astrocyte | cohort_30 | excluding_severely_affected | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| astrocyte | cohort_50 | primary | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant |
| astrocyte | cohort_50 | primary | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant |
| astrocyte | cohort_50 | primary | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| astrocyte | cohort_50 | primary | limma | atac | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | rna | factor_interaction union (exported) | 0 | False | concordant |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | atac | factor_interaction union (exported) | 0 | False | concordant |
| astrocyte | cohort_50 | excluding_severely_affected | limma | rna | factor_interaction interaction_joint (exported) | 0 | False | concordant (DESeq2 concordant; methods agree True) |


### 7.2 Exported-RNA counts at padj ≤ 0.05 (from the summaries)

| cell type | cohort | root | method | model | coefficient | n significant (exported RNA) |
|---|---|---|---|---|---|---|
| microglia | cohort_30 | primary | DESeq2 | factor_additive | stage_1_vs_0 | 0 |
| microglia | cohort_30 | primary | DESeq2 | factor_additive | stage_2_vs_0 | 61 |
| microglia | cohort_30 | primary | DESeq2 | factor_additive | union | 61 |
| microglia | cohort_30 | primary | DESeq2 | factor_interaction | male_x_stage_1 | 0 |
| microglia | cohort_30 | primary | DESeq2 | factor_interaction | male_x_stage_2 | 0 |
| microglia | cohort_30 | primary | DESeq2 | factor_interaction | union | 0 |
| microglia | cohort_30 | primary | DESeq2 | trend_additive | stage_num | 57 |
| microglia | cohort_30 | primary | DESeq2 | trend_interaction | male_x_stage_num | 0 |
| microglia | cohort_30 | primary | limma | factor_additive | stage_1_vs_0 | 0 |
| microglia | cohort_30 | primary | limma | factor_additive | stage_2_vs_0 | 28 |
| microglia | cohort_30 | primary | limma | factor_additive | stage_joint | 2 |
| microglia | cohort_30 | primary | limma | factor_additive | union | 28 |
| microglia | cohort_30 | primary | limma | factor_interaction | male_x_stage_1 | 0 |
| microglia | cohort_30 | primary | limma | factor_interaction | male_x_stage_2 | 0 |
| microglia | cohort_30 | primary | limma | factor_interaction | interaction_joint | 0 |
| microglia | cohort_30 | primary | limma | factor_interaction | union | 0 |
| microglia | cohort_30 | primary | limma | trend_additive | stage_num | 31 |
| microglia | cohort_30 | primary | limma | trend_interaction | male_x_stage_num | 0 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | stage_1_vs_0 | 0 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | stage_2_vs_0 | 20 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | union | 20 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | factor_interaction | male_x_stage_1 | 0 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | factor_interaction | male_x_stage_2 | 0 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | factor_interaction | union | 0 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | trend_additive | stage_num | 21 |
| microglia | cohort_30 | excluding_severely_affected | DESeq2 | trend_interaction | male_x_stage_num | 0 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_1_vs_0 | 0 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_2_vs_0 | 2 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_joint | 2 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_additive | union | 2 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_interaction | male_x_stage_1 | 0 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_interaction | male_x_stage_2 | 0 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_interaction | interaction_joint | 0 |
| microglia | cohort_30 | excluding_severely_affected | limma | factor_interaction | union | 0 |
| microglia | cohort_30 | excluding_severely_affected | limma | trend_additive | stage_num | 2 |
| microglia | cohort_30 | excluding_severely_affected | limma | trend_interaction | male_x_stage_num | 0 |
| microglia | cohort_50 | primary | DESeq2 | factor_additive | stage_1_vs_0 | 1 |
| microglia | cohort_50 | primary | DESeq2 | factor_additive | stage_2_vs_0 | 34 |
| microglia | cohort_50 | primary | DESeq2 | factor_additive | union | 35 |
| microglia | cohort_50 | primary | DESeq2 | factor_interaction | male_x_stage_1 | 0 |
| microglia | cohort_50 | primary | DESeq2 | factor_interaction | male_x_stage_2 | 0 |
| microglia | cohort_50 | primary | DESeq2 | factor_interaction | union | 0 |
| microglia | cohort_50 | primary | DESeq2 | trend_additive | stage_num | 34 |
| microglia | cohort_50 | primary | DESeq2 | trend_interaction | male_x_stage_num | 0 |
| microglia | cohort_50 | primary | limma | factor_additive | stage_1_vs_0 | 0 |
| microglia | cohort_50 | primary | limma | factor_additive | stage_2_vs_0 | 36 |
| microglia | cohort_50 | primary | limma | factor_additive | stage_joint | 0 |
| microglia | cohort_50 | primary | limma | factor_additive | union | 36 |
| microglia | cohort_50 | primary | limma | factor_interaction | male_x_stage_1 | 0 |
| microglia | cohort_50 | primary | limma | factor_interaction | male_x_stage_2 | 0 |
| microglia | cohort_50 | primary | limma | factor_interaction | interaction_joint | 0 |
| microglia | cohort_50 | primary | limma | factor_interaction | union | 0 |
| microglia | cohort_50 | primary | limma | trend_additive | stage_num | 35 |
| microglia | cohort_50 | primary | limma | trend_interaction | male_x_stage_num | 0 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | stage_1_vs_0 | 1 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | stage_2_vs_0 | 2 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | union | 3 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | factor_interaction | male_x_stage_1 | 0 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | factor_interaction | male_x_stage_2 | 0 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | factor_interaction | union | 0 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | trend_additive | stage_num | 0 |
| microglia | cohort_50 | excluding_severely_affected | DESeq2 | trend_interaction | male_x_stage_num | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_1_vs_0 | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_2_vs_0 | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_joint | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_additive | union | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_interaction | male_x_stage_1 | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_interaction | male_x_stage_2 | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_interaction | interaction_joint | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | factor_interaction | union | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | trend_additive | stage_num | 0 |
| microglia | cohort_50 | excluding_severely_affected | limma | trend_interaction | male_x_stage_num | 0 |
| astrocyte | cohort_30 | primary | DESeq2 | factor_additive | stage_1_vs_0 | 0 |
| astrocyte | cohort_30 | primary | DESeq2 | factor_additive | stage_2_vs_0 | 1 |
| astrocyte | cohort_30 | primary | DESeq2 | factor_additive | union | 1 |
| astrocyte | cohort_30 | primary | DESeq2 | factor_interaction | male_x_stage_1 | 0 |
| astrocyte | cohort_30 | primary | DESeq2 | factor_interaction | male_x_stage_2 | 0 |
| astrocyte | cohort_30 | primary | DESeq2 | factor_interaction | union | 0 |
| astrocyte | cohort_30 | primary | DESeq2 | trend_additive | stage_num | 2 |
| astrocyte | cohort_30 | primary | DESeq2 | trend_interaction | male_x_stage_num | 0 |
| astrocyte | cohort_30 | primary | limma | factor_additive | stage_1_vs_0 | 0 |
| astrocyte | cohort_30 | primary | limma | factor_additive | stage_2_vs_0 | 1 |
| astrocyte | cohort_30 | primary | limma | factor_additive | stage_joint | 1 |
| astrocyte | cohort_30 | primary | limma | factor_additive | union | 1 |
| astrocyte | cohort_30 | primary | limma | factor_interaction | male_x_stage_1 | 0 |
| astrocyte | cohort_30 | primary | limma | factor_interaction | male_x_stage_2 | 0 |
| astrocyte | cohort_30 | primary | limma | factor_interaction | interaction_joint | 0 |
| astrocyte | cohort_30 | primary | limma | factor_interaction | union | 0 |
| astrocyte | cohort_30 | primary | limma | trend_additive | stage_num | 1 |
| astrocyte | cohort_30 | primary | limma | trend_interaction | male_x_stage_num | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | stage_1_vs_0 | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | stage_2_vs_0 | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | factor_additive | union | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | factor_interaction | male_x_stage_1 | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | factor_interaction | male_x_stage_2 | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | factor_interaction | union | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | trend_additive | stage_num | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | DESeq2 | trend_interaction | male_x_stage_num | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_1_vs_0 | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_2_vs_0 | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_additive | stage_joint | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_additive | union | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_interaction | male_x_stage_1 | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_interaction | male_x_stage_2 | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_interaction | interaction_joint | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | factor_interaction | union | 0 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | trend_additive | stage_num | 1 |
| astrocyte | cohort_30 | excluding_severely_affected | limma | trend_interaction | male_x_stage_num | 0 |
| astrocyte | cohort_50 | primary | DESeq2 | factor_additive | stage_1_vs_0 | 0 |
| astrocyte | cohort_50 | primary | DESeq2 | factor_additive | stage_2_vs_0 | 2 |
| astrocyte | cohort_50 | primary | DESeq2 | factor_additive | union | 2 |
| astrocyte | cohort_50 | primary | DESeq2 | factor_interaction | male_x_stage_1 | 0 |
| astrocyte | cohort_50 | primary | DESeq2 | factor_interaction | male_x_stage_2 | 0 |
| astrocyte | cohort_50 | primary | DESeq2 | factor_interaction | union | 0 |
| astrocyte | cohort_50 | primary | DESeq2 | trend_additive | stage_num | 2 |
| astrocyte | cohort_50 | primary | DESeq2 | trend_interaction | male_x_stage_num | 0 |
| astrocyte | cohort_50 | primary | limma | factor_additive | stage_1_vs_0 | 0 |
| astrocyte | cohort_50 | primary | limma | factor_additive | stage_2_vs_0 | 1 |
| astrocyte | cohort_50 | primary | limma | factor_additive | stage_joint | 1 |
| astrocyte | cohort_50 | primary | limma | factor_additive | union | 1 |
| astrocyte | cohort_50 | primary | limma | factor_interaction | male_x_stage_1 | 0 |
| astrocyte | cohort_50 | primary | limma | factor_interaction | male_x_stage_2 | 0 |
| astrocyte | cohort_50 | primary | limma | factor_interaction | interaction_joint | 0 |
| astrocyte | cohort_50 | primary | limma | factor_interaction | union | 0 |
| astrocyte | cohort_50 | primary | limma | trend_additive | stage_num | 1 |
| astrocyte | cohort_50 | primary | limma | trend_interaction | male_x_stage_num | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | stage_1_vs_0 | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | stage_2_vs_0 | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | factor_additive | union | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | factor_interaction | male_x_stage_1 | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | factor_interaction | male_x_stage_2 | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | factor_interaction | union | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | trend_additive | stage_num | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | DESeq2 | trend_interaction | male_x_stage_num | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_1_vs_0 | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_2_vs_0 | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_additive | stage_joint | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_additive | union | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_interaction | male_x_stage_1 | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_interaction | male_x_stage_2 | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_interaction | interaction_joint | 0 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | factor_interaction | union | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | trend_additive | stage_num | 1 |
| astrocyte | cohort_50 | excluding_severely_affected | limma | trend_interaction | male_x_stage_num | 1 |


## 8. Donor-exclusion sensitivity

The variant (cohort donors minus SEA-AD's severely affected donors) was rerun for both cell types: the flag check was written before the variant export (about 2 s earlier in each cohort) and is byte-identical to the superseded one, the variant export is byte-identical, and the variant analysis carries the floor. The variant analyses select ranks 2 and 2 (astrocytes) and 4 and 4 (microglia; `cohort_50` already selected 4, and its statistics are identical with the CV AUROC 0.6280 → 0.6450); none rejects any statistic. The predeclared axis check reads `pls_0`, which nests, so its values reproduce: in microglia it holds in both cohorts (stage 0.8201 and 0.8248 against `fraction_mitochondrial_umis_sd` 0.6891 and 0.6746), in astrocytes (descriptive) it does not (stage 0.5784 and 0.5955 against `ribosomal_mean` 0.7286 and 0.6716).

### 8.1 microglia, cohort_30 (≥ 30 microglia nuclei in both omics, n = 76) — primary cohort

Pasted from `results/microglia/motco_exclusion_sensitivity/excluding_severely_affected/comparison/cohort_30/exclusion_sensitivity.md` (the flag check's per-variable rows other than `pls_0` and stage are in the file):

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

##### delta

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 76 | 2 | 5.10924 | 0.689 | no | 10.2655 | 7.68786 | 25.6788 |
| excluding_severely_affected | 69 | 4 | 2.4993 | 0.781 | no | 6.9452 | 5.14197 | 16.822 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 76 | 2 | 12.5272 | 0.424 | no | 12.6966 | 9.51188 | 30.1659 |
| excluding_severely_affected | 69 | 4 | 22.861 | 0.218 | no | 17.1426 | 7.83356 | 31.246 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 76 | 2 | 0.0242177 | 0.99 | no | 0.259273 | 0.142347 | 0.515887 |
| excluding_severely_affected | 69 | 4 | 0.0853244 | 0.569 | no | 0.0989314 | 0.0517103 | 0.186701 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 76 | 0 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=24; Male/stage 0=9; Male/stage 1=11; Male/stage 2=11 | cv | 2 | 0.6852 | 111.7 | `results/microglia/motco_analysis/cohort_30` |
| excluding_severely_affected | 69 | 7 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=19; Male/stage 0=9; Male/stage 1=11; Male/stage 2=9 | cv | 4 | 0.6649 | 108.1 | `motco_analysis/cohort_30` |

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


Flag check rows for `pls_0` and stage (`flag_check.csv`):

| cohort | variable | source | n_flagged | n_unflagged | auc_all | auc_stage2 | spearman_unflagged_stage | spearman_unflagged_pls0 |
|---|---|---|---|---|---|---|---|---|
| cohort_30 | pls_0 | latent_scores.csv (primary analysis) | 7 | 69 | 0.9855072463768116 | 0.9642857142857143 | 0.6434314934106844 |  |
| cohort_30 | stage | metadata.csv | 7 | 69 | 0.7971014492753623 |  |  | 0.6434314934106845 |

`axis_check.csv`:

| cohort | cell_type | role | reading | component | n_donors | stage_spearman | stage_abs_spearman | technical_variable | technical_source | technical_spearman | technical_abs_spearman | n_technical | holds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cohort_30 | microglia | control | expectation | pls_0 | 69 | 0.8201337894290228 | 0.8201337894290228 | fraction_mitochondrial_umis_sd | donor_technical.csv | 0.6891487029594446 | 0.6891487029594446 | 39 | True |


Superseded → floored:

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| microglia cohort_30 primary | 1 → 2 | 0.6921 → 0.6852 | 8.0642 → 5.10924 | 0.391 → 0.689 | no → no | 0 → 12.5272 | 1 → 0.424 | no → no | 0.154684 → 0.0242177 | 0.637 → 0.99 | no → no |
| microglia cohort_30 excluding_severely_affected | 1 → 4 | 0.6388 → 0.6649 | 5.34849 → 2.4993 | 0.379 → 0.781 | no → no | 0 → 22.861 | 1 → 0.218 | no → no | 0.187608 → 0.0853244 | 0.366 → 0.569 | no → no |


### 8.2 microglia, cohort_50 (≥ 50 microglia nuclei in both omics, n = 63) — sensitivity cohort

Pasted from `results/microglia/motco_exclusion_sensitivity/excluding_severely_affected/comparison/cohort_50/exclusion_sensitivity.md` (the flag check's per-variable rows other than `pls_0` and stage are in the file):

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

##### delta

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 63 | 4 | 1.6686 | 0.819 | no | 5.78012 | 4.4338 | 14.0167 |
| excluding_severely_affected | 58 | 4 | 2.55257 | 0.758 | no | 6.67341 | 5.02936 | 16.7839 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 63 | 4 | 22.2604 | 0.227 | no | 16.6216 | 7.94278 | 31.3913 |
| excluding_severely_affected | 58 | 4 | 19.1427 | 0.45 | no | 18.7309 | 8.58051 | 34.5261 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 63 | 4 | 0.0797433 | 0.598 | no | 0.0979335 | 0.0515193 | 0.196376 |
| excluding_severely_affected | 58 | 4 | 0.0880231 | 0.52 | no | 0.0962721 | 0.0495251 | 0.183762 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 63 | 0 | 8 | Female/stage 0=8; Female/stage 1=10; Female/stage 2=18; Male/stage 0=8; Male/stage 1=10; Male/stage 2=9 | cv | 4 | 0.6565 | 107.4 | `results/microglia/motco_analysis/cohort_50` |
| excluding_severely_affected | 58 | 5 | 8 | Female/stage 0=8; Female/stage 1=10; Female/stage 2=14; Male/stage 0=8; Male/stage 1=10; Male/stage 2=8 | cv | 4 | 0.6450 | 97.3 | `motco_analysis/cohort_50` |

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


Flag check rows for `pls_0` and stage (`flag_check.csv`):

| cohort | variable | source | n_flagged | n_unflagged | auc_all | auc_stage2 | spearman_unflagged_stage | spearman_unflagged_pls0 |
|---|---|---|---|---|---|---|---|---|
| cohort_50 | pls_0 | latent_scores.csv (primary analysis) | 5 | 58 | 0.9827586206896551 | 0.9545454545454546 | 0.6989510265184621 |  |
| cohort_50 | stage | metadata.csv | 5 | 58 | 0.8103448275862069 |  |  | 0.6989510265184621 |

`axis_check.csv`:

| cohort | cell_type | role | reading | component | n_donors | stage_spearman | stage_abs_spearman | technical_variable | technical_source | technical_spearman | technical_abs_spearman | n_technical | holds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cohort_50 | microglia | control | expectation | pls_0 | 58 | 0.8248014782033875 | 0.8248014782033875 | fraction_mitochondrial_umis_sd | donor_technical.csv | 0.6746131840413425 | 0.6746131840413425 | 39 | True |


Superseded → floored:

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| microglia cohort_50 primary | 1 → 4 | 0.6414 → 0.6565 | 12.4641 → 1.6686 | 0.227 → 0.819 | no → no | 0 → 22.2604 | 1 → 0.227 | no → no | 0.269568 → 0.0797433 | 0.399 → 0.598 | no → no |
| microglia cohort_50 excluding_severely_affected | 4 → 4 | 0.6280 → 0.6450 | 2.55257 → 2.55257 | 0.758 → 0.758 | no → no | 19.1427 → 19.1427 | 0.45 → 0.45 | no → no | 0.0880231 → 0.0880231 | 0.52 → 0.52 | no → no |


### 8.3 astrocyte, cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

Pasted from `results/motco_exclusion_sensitivity/excluding_severely_affected/comparison/cohort_30/exclusion_sensitivity.md` (the flag check's per-variable rows other than `pls_0` and stage are in the file):

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

##### delta

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 80 | 5 | 5.95773 | 0.468 | no | 6.81334 | 5.28618 | 16.9444 |
| excluding_severely_affected | 70 | 2 | 4.67327 | 0.599 | no | 7.23162 | 5.75846 | 18.3082 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 80 | 5 | 21.9941 | 0.396 | no | 20.7155 | 9.27804 | 37.7882 |
| excluding_severely_affected | 70 | 2 | 10.6202 | 0.503 | no | 12.7808 | 9.619 | 31.1427 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 80 | 5 | 0.15645 | 0.284 | no | 0.124103 | 0.0641842 | 0.24556 |
| excluding_severely_affected | 70 | 2 | 0.248331 | 0.356 | no | 0.214012 | 0.113554 | 0.412346 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 80 | 0 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=28; Male/stage 0=9; Male/stage 1=10; Male/stage 2=12 | cv | 5 | 0.6356 | 115.2 | `results/motco_analysis/cohort_30` |
| excluding_severely_affected | 70 | 10 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=20; Male/stage 0=9; Male/stage 1=10; Male/stage 2=10 | cv | 2 | 0.5892 | 113.6 | `motco_analysis/cohort_30` |

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


Flag check rows for `pls_0` and stage (`flag_check.csv`):

| cohort | variable | source | n_flagged | n_unflagged | auc_all | auc_stage2 | spearman_unflagged_stage | spearman_unflagged_pls0 |
|---|---|---|---|---|---|---|---|---|
| cohort_30 | pls_0 | latent_scores.csv (primary analysis) | 10 | 70 | 0.9614285714285714 | 0.9366666666666666 | 0.3836662026450225 |  |
| cohort_30 | stage | metadata.csv | 10 | 70 | 0.7857142857142857 |  |  | 0.38366620264502244 |

`axis_check.csv`:

| cohort | cell_type | role | reading | component | n_donors | stage_spearman | stage_abs_spearman | technical_variable | technical_source | technical_spearman | technical_abs_spearman | n_technical | holds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cohort_30 | astrocyte | reference | descriptive | pls_0 | 70 | 0.5783987803892782 | 0.5783987803892782 | ribosomal_mean | donor_technical.csv | -0.7285976729944886 | 0.7285976729944886 | 39 | False |


Superseded → floored:

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| astrocyte cohort_30 primary | 1 → 5 | 0.6313 → 0.6356 | 14.941 → 5.95773 | 0.3 → 0.468 | no → no | 0 → 21.9941 | 1 → 0.396 | no → no | 0.608832 → 0.15645 | 0.333 → 0.284 | no → no |
| astrocyte cohort_30 excluding_severely_affected | 1 → 2 | 0.5559 → 0.5892 | 6.35558 → 4.67327 | 0.534 → 0.599 | no → no | 0 → 10.6202 | 1 → 0.503 | no → no | 0.495893 → 0.248331 | 0.173 → 0.356 | no → no |


### 8.4 astrocyte, cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74) — sensitivity cohort

Pasted from `results/motco_exclusion_sensitivity/excluding_severely_affected/comparison/cohort_50/exclusion_sensitivity.md` (the flag check's per-variable rows other than `pls_0` and stage are in the file):

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

##### delta

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 74 | 2 | 6.01569 | 0.596 | no | 9.78355 | 7.88773 | 25.288 |
| excluding_severely_affected | 64 | 2 | 3.27569 | 0.734 | no | 7.77224 | 6.34976 | 20.1871 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 74 | 2 | 23.1596 | 0.222 | no | 15.1566 | 12.4608 | 40.0109 |
| excluding_severely_affected | 64 | 2 | 9.0544 | 0.576 | no | 12.2999 | 8.69314 | 28.6555 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 74 | 2 | 0.284932 | 0.484 | no | 0.296547 | 0.162566 | 0.599664 |
| excluding_severely_affected | 64 | 2 | 0.255222 | 0.336 | no | 0.215239 | 0.118741 | 0.43538 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 74 | 0 | 9 | Female/stage 0=10; Female/stage 1=9; Female/stage 2=25; Male/stage 0=9; Male/stage 1=9; Male/stage 2=12 | cv | 2 | 0.6965 | 114.4 | `results/motco_analysis/cohort_50` |
| excluding_severely_affected | 64 | 10 | 9 | Female/stage 0=10; Female/stage 1=9; Female/stage 2=17; Male/stage 0=9; Male/stage 1=9; Male/stage 2=10 | cv | 2 | 0.6012 | 104.4 | `motco_analysis/cohort_50` |

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


Flag check rows for `pls_0` and stage (`flag_check.csv`):

| cohort | variable | source | n_flagged | n_unflagged | auc_all | auc_stage2 | spearman_unflagged_stage | spearman_unflagged_pls0 |
|---|---|---|---|---|---|---|---|---|
| cohort_50 | pls_0 | latent_scores.csv (primary analysis) | 10 | 64 | 0.959375 | 0.9259259259259259 | 0.39169692303289294 |  |
| cohort_50 | stage | metadata.csv | 10 | 64 | 0.7890625 |  |  | 0.39169692303289294 |

`axis_check.csv`:

| cohort | cell_type | role | reading | component | n_donors | stage_spearman | stage_abs_spearman | technical_variable | technical_source | technical_spearman | technical_abs_spearman | n_technical | holds |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cohort_50 | astrocyte | reference | descriptive | pls_0 | 64 | 0.595467344790454 | 0.595467344790454 | ribosomal_mean | donor_technical.csv | -0.6716117216117216 | 0.6716117216117216 | 39 | False |


Superseded → floored:

| analysis | rank | CV AUROC | delta | delta p | delta rejected | angle | angle p | angle rejected | shape | shape p | shape rejected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| astrocyte cohort_50 primary | 1 → 2 | 0.6616 → 0.6965 | 18.7091 → 6.01569 | 0.213 → 0.596 | no → no | 0 → 23.1596 | 1 → 0.222 | no → no | 0.56471 → 0.284932 | 0.377 → 0.484 | no → no |
| astrocyte cohort_50 excluding_severely_affected | 1 → 2 | 0.5661 → 0.6012 | 7.48789 → 3.27569 | 0.476 → 0.734 | no → no | 0 → 9.0544 | 1 → 0.576 | no → no | 0.384785 → 0.255222 | 0.296 → 0.336 | no → no |


## 9. Positive control (microglia against astrocytes)

E1–E3, their comparisons, technical variables and the all-three rule are unchanged. E1 reads the stage-DE counts (unchanged), E3 the first component (unchanged by nesting), so only E2, the CV mean AUROC, moves with the floor.

### 9.1 basis `primary`, cohort_30

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


Superseded → floored (`positive_control.csv`):

| expectation | method | control value | reference | reference value | outcome |
|---|---|---|---|---|---|
| E1 | deseq2 | 61 → 61 | astrocyte → astrocyte | 1 → 1 | holds → holds |
| E1 | limma | 2 → 2 | astrocyte → astrocyte | 1 → 1 | holds → holds |
| E1 | both | – → – | – → – | – → – | holds → holds |
| E2 |  | 0.6921091408591409 → 0.6852132127132127 | astrocyte → astrocyte | 0.6312777777777778 → 0.6355555555555557 | holds → holds |
| E3 |  | 0.6881745856580481 → 0.6881745856580481 | mitochondrial_pseudobulk → mitochondrial_pseudobulk | 0.8438277511961723 → 0.8438277511961723 | fails → fails |
| overall |  | – → – | – → – | – → – | does_not_hold → does_not_hold |


### 9.2 basis `primary`, cohort_50

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


Superseded → floored (`positive_control.csv`):

| expectation | method | control value | reference | reference value | outcome |
|---|---|---|---|---|---|
| E1 | deseq2 | 35 → 35 | astrocyte → astrocyte | 2 → 2 | holds → holds |
| E1 | limma | 0 → 0 | astrocyte → astrocyte | 1 → 1 | fails → fails |
| E1 | both | – → – | – → – | – → – | method_dependent → method_dependent |
| E2 |  | 0.6414312369729036 → 0.6564882756132756 | astrocyte → astrocyte | 0.6616446773391218 → 0.6965036403925293 | fails → fails |
| E3 |  | 0.7307728100677887 → 0.7307728100677887 | mitochondrial_pseudobulk → mitochondrial_pseudobulk | 0.784562211981567 → 0.784562211981567 | fails → fails |
| overall |  | – → – | – → – | – → – | does_not_hold → does_not_hold |


### 9.3 basis `latent_adjusted`, cohort_30

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


Superseded → floored (`positive_control.csv`):

| expectation | method | control value | reference | reference value | outcome |
|---|---|---|---|---|---|
| E1 | deseq2 | 61 → 61 | astrocyte → astrocyte | 1 → 1 | holds → holds |
| E1 | limma | 2 → 2 | astrocyte → astrocyte | 1 → 1 | holds → holds |
| E1 | both | – → – | – → – | – → – | holds → holds |
| E2 |  | 0.6415025252525253 → 0.6574797887297887 | astrocyte → astrocyte | 0.6299444444444445 → 0.6433333333333333 | holds → holds |
| E3 |  | 0.5278766546319672 → 0.5278766546319672 | tmm_rna → tmm_rna | 0.4152563226247436 → 0.4152563226247436 | holds → holds |
| overall |  | – → – | – → – | – → – | holds → holds |


### 9.4 basis `latent_adjusted`, cohort_50

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


Superseded → floored (`positive_control.csv`):

| expectation | method | control value | reference | reference value | outcome |
|---|---|---|---|---|---|
| E1 | deseq2 | 35 → 35 | astrocyte → astrocyte | 2 → 2 | holds → holds |
| E1 | limma | 0 → 0 | astrocyte → astrocyte | 1 → 1 | fails → fails |
| E1 | both | – → – | – → – | – → – | method_dependent → method_dependent |
| E2 |  | 0.6204978354978354 → 0.6260852573352572 | astrocyte → astrocyte | 0.6644790677429566 → 0.6644790677429566 | fails → fails |
| E3 |  | 0.5326807217991494 → 0.5326807217991494 | gex_reads_mapped_confidently_to_intronic_regions_sd → gex_reads_mapped_confidently_to_intronic_regions_sd | 0.4067015510422894 → 0.4067015510422894 | holds → holds |
| overall |  | – → – | – → – | – → – | does_not_hold → does_not_hold |


## 10. Reading

- **The floored primaries are the primary results.** In both cell types and both cohorts the floored analysis selects a rank of at least 2 and rejects none of `delta`, `angle` and `shape` at α = 0.05. Under the binding rules: a `delta` non-rejection is not evidence against a size difference (the adjusted analysis carries no power claim); the `angle` and `shape` non-rejections are uninformative about orientation and shape at this n. The floor guarantees that `angle` and `shape` are measured in at least two dimensions, not that a second axis carries signal: the observed eigengaps sit inside their permutation distributions (section 2), and MOTCO's floor check found power at the null at its strengths.
- **No dependence on the rank rule.** No rejection status differs between a floored analysis and its superseded one, in any step. The predeclared case "status differs after the move" did not occur.
- **Sensitivity readings, re-read against the floored primary under their own rules.** Rank: no forced rank rejects, no rank dependence; the CV-selected rank 5 of astrocyte `cohort_30` is outside the ladder and the reference row stands in for it. Preprocessing: no variant rejects; the adoption rule holds for no candidate in both cohorts and the precedence names none, so the variants stay sensitivity analyses (unchanged outcome; one candidate cell flipped as the design anticipated). Latent adjustment: no arm rejects; the `extended_covariates` arms keep their primaries' rank and AUROC; the mediator limitation stands with its flags. Donor exclusion: no variant rejects; the microglia axis check holds and the astrocyte one does not, as before. Stage DE: every concordance cell concordant, methods agree.
- **Positive control.** On the floored primary the control's outcome in `cohort_30` is unchanged: E1 holds (DESeq2 61 against 1; limma 2 against 1), E2 holds (0.6852 against 0.6356; superseded 0.6921 against 0.6313), E3 fails (0.6882 against 0.8438; superseded the same values), so the control **does not hold**, through E3 alone; `cohort_50` stays descriptive (E1 `method_dependent`, E2 and E3 fail). The restated control on `latent_adjusted` keeps its predeclared reading and its outcome: in `cohort_30` it holds (E1 holds, E2 0.6575 against 0.6433, E3 0.5279 against 0.4153), so the measured quality gradient explains the primary-basis failure, without promoting the arm and with the mediator caveat; `cohort_50` stays descriptive (E2 fails, E3 holds, E1 `method_dependent`). The reading of the control is the earlier reports' reading: about the stage signal, not the sex × stage hypothesis.
- **Expected reproductions.** `pls_0` (40 of 40), E3, the axis checks, the flag checks, the mediator checks, the variant exports, the forced-rank records and the stage-DE counts and calls reproduce. The per-feature astrocyte DESeq2 statistics do not reproduce beyond rounding across machines; this changes no count, call, concordance cell or reading.

## 11. Findings

- **Astrocyte DESeq2 across machines.** The 24 per-feature astrocyte DESeq2 tables re-fitted on `nixos-home` differ from the superseded ones fitted on `nixos-dell-laptop` (2026-10-21, same pydeseq2 0.5.3, numpy 2.3.5, Python 3.12.15, `--n-cpus 10`) beyond rounding: `baseMean` agrees to within 6e-11, but the dispersion-dependent columns differ (largest |Δ log2FoldChange| 0.0119, |Δ stat| 0.316, |Δ padj| 0.0582). Every padj ≤ 0.05 call, every count and every concordance cell is equal. Same-machine refits are byte-identical (the microglia tables and the astrocyte exclusion-root tables, both previously fitted here), so the differences come from the machine, through pydeseq2's iterative dispersion fits. Recorded, no reading changed.
- **Kept astrocyte primary across machines.** The floor-1 rerun on this machine differs from the kept cluster-made astrocyte primary at rounding level only (section 1.2); the gate itself is byte-identical.
- **CV AUROC under an unchanged rank.** Where the floored CV selects the same rank as the superseded one (astrocyte `cohort_50` `latent_adjusted` and `residualized`; microglia `cohort_50` `latent_adjusted` and donor exclusion; astrocyte `cohort_50` `revised_candidate`), the statistics, p-values and null summaries are identical, and the CV mean AUROC changes in three of them, because it averages outer folds whose selected ranks the floor now constrains.

Generated by the comparison script (superseded against re-issued):

- 24 per-feature DESeq2 tables differ from their superseded copies beyond rounding (largest |Δ stat| 0.316, largest |Δ padj| 0.0582); padj <= 0.05 calls equal in 24 of 24 (table in de_differences.md)

### 11.1 Per-feature DESeq2 tables against their superseded copies (astrocytes; generated)

| table | max |Δ log2FoldChange| | max |Δ lfcSE| | max |Δ stat| | max |Δ pvalue| | max |Δ padj| | calls equal | n calls (old, new) |
|---|---|---|---|---|---|---|---|
| `stage_de/cohort_30/atac/factor_additive__stage_1_vs_0.csv` | 0.00312 | 0.00585 | 0.117 | 0.0716 | 0.00398 | yes | 0, 0 |
| `stage_de/cohort_30/atac/factor_additive__stage_2_vs_0.csv` | 0.00333 | 0.00501 | 0.16 | 0.0733 | 0.0419 | yes | 22, 22 |
| `stage_de/cohort_30/atac/factor_interaction__male_x_stage_1.csv` | 0.00493 | 0.0148 | 0.188 | 0.0509 | 8.69e-08 | yes | 0, 0 |
| `stage_de/cohort_30/atac/factor_interaction__male_x_stage_2.csv` | 0.00822 | 0.013 | 0.289 | 0.0616 | 0.0242 | yes | 0, 0 |
| `stage_de/cohort_30/atac/trend_additive__stage_num.csv` | 0.00231 | 0.00348 | 0.163 | 0.0664 | 0.0449 | yes | 23, 23 |
| `stage_de/cohort_30/atac/trend_interaction__male_x_stage_num.csv` | 0.00317 | 0.00918 | 0.182 | 0.0686 | 0.0101 | yes | 0, 0 |
| `stage_de/cohort_30/rna/factor_additive__stage_1_vs_0.csv` | 0.000887 | 0.00236 | 0.0415 | 0.0153 | 0.000655 | yes | 0, 0 |
| `stage_de/cohort_30/rna/factor_additive__stage_2_vs_0.csv` | 0.0011 | 0.00206 | 0.092 | 0.0183 | 0.0161 | yes | 1, 1 |
| `stage_de/cohort_30/rna/factor_interaction__male_x_stage_1.csv` | 0.000937 | 0.00546 | 0.0343 | 0.025 | 4.48e-05 | yes | 0, 0 |
| `stage_de/cohort_30/rna/factor_interaction__male_x_stage_2.csv` | 0.00195 | 0.0048 | 0.0506 | 0.0302 | 0.00992 | yes | 0, 0 |
| `stage_de/cohort_30/rna/trend_additive__stage_num.csv` | 0.000434 | 0.291 | 0.231 | 0.161 | 0.0582 | yes | 2, 2 |
| `stage_de/cohort_30/rna/trend_interaction__male_x_stage_num.csv` | 0.00154 | 0.00314 | 0.0711 | 0.0262 | 0.0217 | yes | 0, 0 |
| `stage_de/cohort_50/atac/factor_additive__stage_1_vs_0.csv` | 0.0029 | 0.00675 | 0.125 | 0.0426 | 7.75e-08 | yes | 0, 0 |
| `stage_de/cohort_50/atac/factor_additive__stage_2_vs_0.csv` | 0.00428 | 0.0056 | 0.222 | 0.032 | 0.0544 | yes | 36, 36 |
| `stage_de/cohort_50/atac/factor_interaction__male_x_stage_1.csv` | 0.00865 | 0.0146 | 0.309 | 0.038 | 5e-08 | yes | 0, 0 |
| `stage_de/cohort_50/atac/factor_interaction__male_x_stage_2.csv` | 0.0119 | 0.0127 | 0.316 | 0.0272 | 0.0174 | yes | 0, 0 |
| `stage_de/cohort_50/atac/trend_additive__stage_num.csv` | 0.00272 | 0.00328 | 0.143 | 0.063 | 0.039 | yes | 30, 30 |
| `stage_de/cohort_50/atac/trend_interaction__male_x_stage_num.csv` | 0.00187 | 0.00602 | 0.0896 | 0.0589 | 0.0202 | yes | 0, 0 |
| `stage_de/cohort_50/rna/factor_additive__stage_1_vs_0.csv` | 0.000621 | 0.00272 | 0.0537 | 0.0221 | 0.00566 | yes | 0, 0 |
| `stage_de/cohort_50/rna/factor_additive__stage_2_vs_0.csv` | 0.000769 | 0.00233 | 0.0437 | 0.0135 | 0.014 | yes | 2, 2 |
| `stage_de/cohort_50/rna/factor_interaction__male_x_stage_1.csv` | 0.00751 | 0.013 | 0.172 | 0.0652 | 3.64e-05 | yes | 0, 0 |
| `stage_de/cohort_50/rna/factor_interaction__male_x_stage_2.csv` | 0.00961 | 0.0111 | 0.128 | 0.02 | 0.0136 | yes | 0, 0 |
| `stage_de/cohort_50/rna/trend_additive__stage_num.csv` | 0.000194 | 0.00126 | 0.108 | 0.0303 | 0.0278 | yes | 2, 2 |
| `stage_de/cohort_50/rna/trend_interaction__male_x_stage_num.csv` | 0.00103 | 0.00249 | 0.0593 | 0.0238 | 0.00843 | yes | 0, 0 |

### 11.2 First latent component against the superseded analyses (generated)

| analysis | allclose | max_abs_diff | max_abs_value | sign_flipped_allclose |
|---|---|---|---|---|
| astrocyte primary cohort_30 | yes | 2.629e-13 | 87.21 | no |
| astrocyte extended_covariates cohort_30 | yes | 2.629e-13 | 87.21 | no |
| astrocyte latent_adjusted cohort_30 | yes | 5.755e-13 | 91.81 | no |
| astrocyte excluding_severely_affected cohort_30 | yes | 0 | 50.52 | no |
| astrocyte forced rank 2 cohort_30 | yes | 0 | 87.21 | no |
| astrocyte forced rank 3 cohort_30 | yes | 0 | 87.21 | no |
| astrocyte forced rank 4 cohort_30 | yes | 0 | 87.21 | no |
| astrocyte forced rank 6 cohort_30 | yes | 0 | 87.21 | no |
| astrocyte forced rank 9 cohort_30 | yes | 0 | 87.21 | no |
| astrocyte forced rank 12 cohort_30 | yes | 0 | 87.21 | no |
| astrocyte variant keep_chrx cohort_30 | yes | 0 | 87.45 | no |
| astrocyte variant exclude_families cohort_30 | yes | 0 | 86.94 | no |
| astrocyte variant rna_trend_ratio cohort_30 | yes | 0 | 85.39 | no |
| astrocyte variant atac_mean_accessibility cohort_30 | yes | 0 | 139.4 | no |
| astrocyte variant residualized cohort_30 | yes | 0 | 34.17 | no |
| astrocyte variant revised_candidate cohort_30 | yes | 0 | 44.83 | no |
| astrocyte primary cohort_50 | yes | 7.816e-14 | 84.63 | no |
| astrocyte extended_covariates cohort_50 | yes | 7.816e-14 | 84.63 | no |
| astrocyte latent_adjusted cohort_50 | yes | 5.471e-13 | 64.27 | no |
| astrocyte excluding_severely_affected cohort_50 | yes | 0 | 47.85 | no |
| astrocyte forced rank 2 cohort_50 | yes | 0 | 84.63 | no |
| astrocyte forced rank 3 cohort_50 | yes | 0 | 84.63 | no |
| astrocyte forced rank 4 cohort_50 | yes | 0 | 84.63 | no |
| astrocyte forced rank 6 cohort_50 | yes | 0 | 84.63 | no |
| astrocyte forced rank 9 cohort_50 | yes | 0 | 84.63 | no |
| astrocyte forced rank 12 cohort_50 | yes | 0 | 84.63 | no |
| astrocyte variant keep_chrx cohort_50 | yes | 0 | 85.22 | no |
| astrocyte variant exclude_families cohort_50 | yes | 0 | 84.61 | no |
| astrocyte variant rna_trend_ratio cohort_50 | yes | 0 | 85.1 | no |
| astrocyte variant atac_mean_accessibility cohort_50 | yes | 0 | 148.2 | no |
| astrocyte variant residualized cohort_50 | yes | 0 | 40.13 | no |
| astrocyte variant revised_candidate cohort_50 | yes | 0 | 56.84 | no |
| microglia primary cohort_30 | yes | 0 | 69.84 | no |
| microglia extended_covariates cohort_30 | yes | 0 | 69.84 | no |
| microglia latent_adjusted cohort_30 | yes | 2.789e-13 | 58.36 | no |
| microglia excluding_severely_affected cohort_30 | yes | 0 | 40.09 | no |
| microglia primary cohort_50 | yes | 0 | 75.67 | no |
| microglia extended_covariates cohort_50 | yes | 0 | 75.67 | no |
| microglia latent_adjusted cohort_50 | yes | 0 | 62.15 | no |
| microglia excluding_severely_affected cohort_50 | yes | 0 | 40.49 | no |

Stage-DE / limmatables compared: 216; byte-identical: 188.

Variant export CSVs compared: 100; byte-identical: 100; without a superseded copy: 30.

## 12. Limitations

- The floor's interaction with the covariate adjustment and with the latent adjustment was not measured by MOTCO: its floor check is at the Phase 6 SEA-AD design weakened to this case study's AUROC regime, without covariates and without `covariate_latent_adjustment`. No Type I claim is made for the floored covariate-adjusted or latent-adjusted analyses beyond what MOTCO measured.
- MOTCO's floor check found the floored `angle` and `shape` at the null at its strengths (power at the null); the floor defines the statistics, it does not make them powerful. No operating characteristic was measured on these blocks, on these cohorts or at six or eight covariates.
- The CV-selected ranks (2 to 7) are now data-dependent above the floor; a rank differing between analyses is reported beside `angle`, never read as better or worse, and the ladder's forced ranks are sensitivity only.
- Every limitation of the earlier reports carries over unchanged: the additive-linear adjustment on a covariate-blind basis; the structural-analogue operating characteristics (InterSIM blocks, 498 features); the mediator limitation of the latent adjustment; the exclusion variant's truncation of the late epoch and quality conditioning; the stage-DE methods' approximate FDR and no joint error-rate claim for the unions.
- The astrocyte DESeq2 per-feature statistics depend on the machine beyond rounding (section 11); counts and calls do not, on the evidence of this comparison.
