# MOTCO astrocyte preprocessing sensitivity (chrX kept): sex × ADNC stage — 2026-10-19

> Preprocessing sensitivity analysis of the covariate-adjusted comparison of female and male astrocyte trajectories
> across ADNC stages in the SEA-AD MTG (OpenSpec change `mtg-astrocyte-chrx-sensitivity`, spec
> `motco-preprocessing-sensitivity`). The primary analysis is the CV-selected one in `reports/motco-astrocyte-2026-10-19.md`,
> whose export removed chrX, chrY and chrM features before any filtering; this report reruns the export with chrX kept
> (variant `keep_chrx` of the registry `PREPROCESSING_VARIANTS`, chrY and chrM still removed), reruns the analysis unchanged
> on the variant blocks, and reads the result as sensitivity to the primary analysis, never as its replacement. Every
> statistic, p-value, null summary, rank and feature count below is pasted by script from the generated
> `sensitivity.md` of each cohort (`results/motco_sensitivity/keep_chrx/comparison/<cohort>/`), which is rendered from the
> `result.json` of the primary and of the variant analysis and from both exports' `features_<omic>.csv`, and was checked
> against those four records and the two comparison manifests before commit. The stage means come from `stage_means.csv`
> and the latent-score correlation from `latent_scores.csv`.

**Run:** `sead-sexdiff preprocessing-sensitivity` (all variants, all cohorts) on the NixOS workstation `nixos-home`,
2026-10-19, at commit `50bbb5a` with the change's implementation uncommitted (the manifests record `dirty: true`),
`n_jobs = 1`: variant analyses of 104.4 s and 94.6 s (CV rank selection
included) beside primary analyses of 148.7 s and 151.6 s. The primary
export and analyses were produced on the laptop workstation and pulled from the cluster; before the run the primary export
was reproduced on this machine byte for byte (section 5) and the step checked the pseudobulk, covariates and annotation
hashes, the thresholds, `k` and the cells against the primary export manifest, and the five export hashes, MOTCO version,
seed, permutations and CV selection against the primary analysis manifest. **Result:** nothing is rejected at
α = 0.05 in the variant or the primary analysis of either cohort, CV selects rank 1 in both as before, and the
statistics, p-values and null summaries of the two rows are nearly equal. Keeping chrX gives no evidence of a difference
either, and that is not evidence of agreement in size, direction or shape. The primary analysis remains the primary result.

**In one paragraph.** For each cohort, the export was rerun from the synced pseudobulk outputs with `excluded_chromosomes = (chrY, chrM)` and nothing else changed, so library sizes, the CPM ≥ 1 in ≥ k filter, TMM factors, log2-CPM and the top-5,000 variance ranking were all recomputed with chrX present, and MOTCO 0.8.0's `analyze_trajectories` was then called on the variant blocks with the primary analysis's settings (six covariates, net-displacement orientation, CV rank selection, 999 RRPP permutations, seed 20261019). In the primary cohort 125 chrX genes (`XIST` and `TSIX` first) and 102 chrX peaks entered the blocks, displacing 125 genes and 102 peaks; 4,875 genes and 4,898 peaks are shared with the primary export. CV again selected rank 1 in both cohorts, the single component correlates 0.99997 / 0.99993 with the primary's, and every statistic and p-value moves by a hair: `delta` 14.941 → 14.9283 (p 0.3 → 0.301) and `shape` 0.608832 → 0.607073 (p 0.333 → 0.337) in the primary cohort, `delta` 18.7091 → 18.9249 (p 0.213 → 0.212) and `shape` 0.56471 → 0.556064 (p 0.377 → 0.391) in the sensitivity cohort, with `angle` 0° (p = 1) in every row because the latent space has one dimension. No rejection status changes, so there is no preprocessing dependence to report; that is a statement about this one setting, not evidence that the trajectories agree, and the primary analysis stays the primary result.

## 1. Provenance

| item | value |
|---|---|
| MOTCO | 0.8.0, tag `v0.8.0`, commit `979657c3d0e6c9ad4d4bff0dc595a18df854ca61` (pinned in `uv.lock`) |
| entry point | `sead_sexdiff.data.run_preprocessing_sensitivity`: `_export_cohort` (the `export` step's per-cohort helper) with `excluded_chromosomes = ("chrY", "chrM")`, then `motco.analysis.analyze_trajectories` + `write_analysis_outputs` through `_analyze_cohort` (the `analyze` step's helper) with `TrajectoryAnalysisParams(permutations=999, seed=20261019)` |
| sead_sexdiff | `50bbb5a` (the change's implementation uncommitted at run time; both comparison manifests record `dirty: true`) |
| machine | `nixos-home` (NixOS workstation); the primary export and analyses were produced on the laptop workstation and pulled from the cluster, the pseudobulk outputs were synced with `scripts/sync_pseudobulk.sh` |
| pseudobulk inputs | `results/pseudobulk/` (`manifest.json` sha256 `2e93cc0b5268…`, created 2026-10-16T18:46:13+00:00; `covariates_manifest.json` sha256 `b8bf9091b507…`), equal to the hashes the primary export recorded |
| annotation | GENCODE v32 basic, `data/annotation/gencode.v32.basic.annotation.gtf.gz` sha256 `a61ab50fc8b5…`, equal to the primary export's |
| primary export | `results/motco_input/<cohort>/` (`reports/motco-astrocyte-2026-10-19.md` §4): manifest sha256 `18f2289e0d4c…` (cohort_30), `dd14a6e7d558…` (cohort_50); `k` = 9 / 9; `excluded_chromosomes` ['chrX', 'chrY', 'chrM']; read and checked, never rerun or modified |
| primary analyses | `results/motco_analysis/<cohort>/` (sead_sexdiff `c2f25e3`, `dirty: false`): CV-selected rank 1 in both cohorts; `result.json` sha256 `ca92876a1a61…` (cohort_30), `1e7eefe3e63e…` (cohort_50); `analysis_manifest.json` sha256 `17b2150ecba0…`, `43fe406690fa…`; read and checked, never rerun or modified |
| variant | `keep_chrx`: chrX kept; chrY and chrM removed; every other export setting equal to the primary; export setting changed: `excluded_chromosomes = ['chrY', 'chrM']`; every other export setting (`k`, CPM ≥ 1 in ≥ k donors, TMM, prior count 2, N = 5,000 per omic, the annotation) equal to the primary's; the registry is a constant and no variant was added, dropped or chosen from the results |
| variant export | `results/motco_sensitivity/keep_chrx/motco_input/<cohort>/`: manifest sha256 `a55d9c224193…` (cohort_30), `e4cba3cd44cf…` (cohort_50); `metadata.csv` and `covariates.csv` byte-identical to the primary export's |
| variant analyses | `results/motco_sensitivity/keep_chrx/motco_analysis/<cohort>/`: `result.json` sha256 `a2634732d478…` (cohort_30), `26fc46233a91…` (cohort_50); `component_selection = cv` |
| measurement | as the primary analysis: blocks `rna`, `atac` in identity units (log2-CPM); group `sex`; stages 0, 1, 2; covariates `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` (additive, `missing_covariates = error`); net-displacement orientation; double-CV rank selection at MOTCO's defaults; no attribution |
| permutations, seed | 999 RRPP permutations, seed 20261019, equal to the primary analysis's, so both rows draw the same permutation indices; `n_jobs` 1 (both cohorts) |
| α | 0.05, rejected when p ≤ α (`p <= alpha`) |
| software | python 3.12.15, numpy 2.3.5, scipy 1.15.3, scikit-learn 1.8.0, pandas 2.3.3 |

## 2. Primary cohort: cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80)

Primary: `results/motco_analysis/cohort_30` (CV-selected rank 1, CV mean AUROC 0.6313;
`result.json` sha256 `ca92876a1a61…`, `analysis_manifest.json` sha256
`17b2150ecba0…`; export manifest sha256 `18f2289e0d4c…`), read and
checked, never rerun. Variant: `results/motco_sensitivity/keep_chrx/motco_analysis/cohort_30` (CV-selected rank 1, CV mean
AUROC 0.6294; `result.json` sha256 `a2634732d478…`; variant export
manifest sha256 `a55d9c224193…`). Nothing is rejected at α = 0.05 in either row; the smallest p-value
outside the degenerate `angle` row is `delta` in the primary analysis (p = 0.3). The two rows select
the same rank (1), so `angle` is 0° by construction in both and carries no orientation information. Findings
recorded by the step: keep_chrx export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped').

chrX kept; chrY and chrM removed; every other export setting equal to the primary. Export setting(s) changed: `excluded_chromosomes` = ['chrY', 'chrM'] (primary: `excluded_chromosomes` = ['chrX', 'chrY', 'chrM']); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `keep_chrx` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse.

### delta (pasted from `sensitivity.md`)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| keep_chrx | cv | 1 | 14.9283 | 0.301 | no | 11.3597 | 8.31813 | 26.8093 |

### angle (degrees) (pasted from `sensitivity.md`)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| keep_chrx | cv | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |

### shape (pasted from `sensitivity.md`)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| keep_chrx | cv | 1 | 0.607073 | 0.337 | no | 0.492838 | 0.357284 | 1.1929 |

### Features (pasted from `sensitivity.md`)

Counts after each export step per omic; `kept on chrX` and `kept only here` are per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after variance) | kept on chrX | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20387 | 5000 | 0 | 4875 | 125 |
| rna | keep_chrx | 36601 | 36477 | 20995 | 5000 | 125 | 4875 | 125 |
| atac | primary | 218882 | 212816 | 128916 | 5000 | 0 | 4898 | 102 |
| atac | keep_chrx | 218882 | 218702 | 130697 | 5000 | 102 | 4898 | 102 |

- `XIST` (chrX, rna): status `sex_chromosome` in the primary export, `kept` in the keep_chrx export.

### Runs (pasted from `sensitivity.md`)

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6313 | 148.7 | `results/motco_analysis/cohort_30` |
| keep_chrx | cv | 1 | 0.6294 | 104.4 | `motco_analysis/cohort_30` |

### Findings (pasted from `sensitivity.md`)

- keep_chrx export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Descriptive checks

These checks describe the geometry behind the tables and carry no inference.

- **The single latent component is the same axis.** The Pearson correlation of `pls_0` scores between the primary and
  the variant analysis over the 80 donors is 0.99997 (same sign). The latent space is
  fitted on 5,000 + 5,000 features of which 4,875 RNA and 4,898 ATAC are
  shared by the two exports, so the axis moves very little when 125 chrX genes and 102 chrX
  peaks replace the 125 genes and 102 peaks that only the primary
  export keeps.
- **Stage means along `pls_0`** (LS-means at the pooled covariate means, `stage_means.csv`) shift by at most a few
  tenths of a unit, and the female path stays the longer one in both analyses:

| group | analysis | stage 0 | stage 1 | stage 2 | net displacement (2 − 0) along `pls_0` |
|---|---|---:|---:|---:|---:|
| Female | primary | -11.75 | -5.33 | 20.89 | 32.63 |
| Female | keep_chrx | -11.58 | -5.14 | 21.02 | 32.61 |
| Male | primary | -24.11 | -10.83 | -6.42 | 17.69 |
| Male | keep_chrx | -24.34 | -11.08 | -6.66 | 17.68 |

- **What chrX contributed.** Of the 1,148 chrX genes, 464 fail
  the variant's own abundance filter (CPM ≥ 1 in ≥ 9 donors), 559
  pass it but fall below the variance cut and 125 are kept; of the 5,886 chrX peaks,
  2,779 fail abundance, 3,005 fall
  below the variance cut and 102 are kept. `XIST` has status `sex_chromosome` in the primary export and
  `kept` in the variant export, and its log2-CPM variance (23.20) is the
  largest of all 5,000 kept RNA features (rank 1; the kept RNA
  variances run from 0.65 upward, the kept ATAC variances from
  1.18 to 7.78). The ten kept chrX genes
  with the highest variance are `XIST` (23.20), `TSIX` (5.45), `SRPX2` (3.03), `GABRE` (3.02), `ALG13-AS1` (2.44), `AL035427.1` (2.25), `MTCP1` (2.11), `AL807742.1` (2.02), `DIPK2B` (1.95), `GK-AS1` (1.90); the ten kept chrX peaks are `chrX:73849830-73853086` (5.21), `chrX:131791829-131793778` (4.12), `chrX:131740917-131742318` (3.81), `chrX:131704058-131705198` (3.53), `chrX:131728490-131729969` (3.39), `chrX:131752284-131753684` (3.37), `chrX:131821916-131822742` (2.99), `chrX:131746492-131747458` (2.86), `chrX:19127-20649` (2.79), `chrX:131762647-131764200` (2.79), the first of them at the `XIST` locus and most of
  the others in one cluster near chrX:131.7–131.8 Mb.

## 3. Sensitivity cohort: cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74)

Primary: `results/motco_analysis/cohort_50` (CV-selected rank 1, CV mean AUROC 0.6616;
`result.json` sha256 `1e7eefe3e63e…`, `analysis_manifest.json` sha256
`43fe406690fa…`; export manifest sha256 `dd14a6e7d558…`), read and
checked, never rerun. Variant: `results/motco_sensitivity/keep_chrx/motco_analysis/cohort_50` (CV-selected rank 1, CV mean
AUROC 0.6642; `result.json` sha256 `26fc46233a91…`; variant export
manifest sha256 `e4cba3cd44cf…`). Nothing is rejected at α = 0.05 in either row; the smallest p-value
outside the degenerate `angle` row is `delta` in the keep_chrx analysis (p = 0.212). The two rows select
the same rank (1), so `angle` is 0° by construction in both and carries no orientation information. Findings
recorded by the step: keep_chrx export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped').

chrX kept; chrY and chrM removed; every other export setting equal to the primary. Export setting(s) changed: `excluded_chromosomes` = ['chrY', 'chrM'] (primary: `excluded_chromosomes` = ['chrX', 'chrY', 'chrM']); every other export setting and every analysis setting equal the primary's. The `primary` row is the CV-selected primary analysis at rank 1, which stays the primary result; the `keep_chrx` row (CV-selected rank 1) is read as sensitivity to it, never as its replacement, and no variant is chosen from the results. 999 RRPP permutations with seed 20261019 in both rows (the same permutation indices); α = 0.05 (rejected when p ≤ α). `angle` and `shape` non-rejections are uninformative in both rows; a `delta` rejection is a difference, not specifically of size; a rejection status that differs between the rows is preprocessing dependence of the sensitivity analysis, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse.

### delta (pasted from `sensitivity.md`)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| keep_chrx | cv | 1 | 18.9249 | 0.212 | no | 11.7631 | 9.29566 | 30.2908 |

### angle (degrees) (pasted from `sensitivity.md`)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| keep_chrx | cv | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |

### shape (pasted from `sensitivity.md`)

| variant | selection | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---|---:|---:|---:|:---:|---:|---:|---:|
| primary | cv | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| keep_chrx | cv | 1 | 0.556064 | 0.391 | no | 0.498063 | 0.344727 | 1.16125 |

### Features (pasted from `sensitivity.md`)

Counts after each export step per omic; `kept on chrX` and `kept only here` are per export, `kept shared` is the number of kept features common to both exports.

| omic | export | input | after chromosome removal | after abundance | kept (after variance) | kept on chrX | kept shared | kept only here |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| rna | primary | 36601 | 35329 | 20238 | 5000 | 0 | 4870 | 130 |
| rna | keep_chrx | 36601 | 36477 | 20820 | 5000 | 129 | 4870 | 130 |
| atac | primary | 218882 | 212816 | 124818 | 5000 | 0 | 4885 | 115 |
| atac | keep_chrx | 218882 | 218702 | 126408 | 5000 | 110 | 4885 | 115 |

- `XIST` (chrX, rna): status `sex_chromosome` in the primary export, `kept` in the keep_chrx export.

### Runs (pasted from `sensitivity.md`)

| variant | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---|---:|---:|---:|---|
| primary | cv | 1 | 0.6616 | 151.6 | `results/motco_analysis/cohort_50` |
| keep_chrx | cv | 1 | 0.6642 | 94.6 | `motco_analysis/cohort_50` |

### Findings (pasted from `sensitivity.md`)

- keep_chrx export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### Descriptive checks

These checks describe the geometry behind the tables and carry no inference.

- **The single latent component is the same axis.** The Pearson correlation of `pls_0` scores between the primary and
  the variant analysis over the 74 donors is 0.99993 (same sign). The latent space is
  fitted on 5,000 + 5,000 features of which 4,870 RNA and 4,885 ATAC are
  shared by the two exports, so the axis moves very little when 129 chrX genes and 110 chrX
  peaks replace the 130 genes and 115 peaks that only the primary
  export keeps.
- **Stage means along `pls_0`** (LS-means at the pooled covariate means, `stage_means.csv`) shift by at most a few
  tenths of a unit, and the female path stays the longer one in both analyses:

| group | analysis | stage 0 | stage 1 | stage 2 | net displacement (2 − 0) along `pls_0` |
|---|---|---:|---:|---:|---:|
| Female | primary | -16.43 | -3.39 | 20.72 | 37.15 |
| Female | keep_chrx | -16.41 | -3.14 | 20.89 | 37.31 |
| Male | primary | -22.95 | -6.94 | -4.51 | 18.44 |
| Male | keep_chrx | -23.12 | -7.22 | -4.74 | 18.38 |

- **What chrX contributed.** Of the 1,148 chrX genes, 465 fail
  the variant's own abundance filter (CPM ≥ 1 in ≥ 9 donors), 554
  pass it but fall below the variance cut and 129 are kept; of the 5,886 chrX peaks,
  2,893 fail abundance, 2,883 fall
  below the variance cut and 110 are kept. `XIST` has status `sex_chromosome` in the primary export and
  `kept` in the variant export, and its log2-CPM variance (23.60) is the
  largest of all 5,000 kept RNA features (rank 1; the kept RNA
  variances run from 0.64 upward, the kept ATAC variances from
  1.16 to 7.85). The ten kept chrX genes
  with the highest variance are `XIST` (23.60), `TSIX` (5.66), `SRPX2` (3.04), `GABRE` (3.02), `ALG13-AS1` (2.49), `AL035427.1` (2.31), `MTCP1` (2.19), `AL807742.1` (1.98), `GK-AS1` (1.90), `DGAT2L6` (1.72); the ten kept chrX peaks are `chrX:73849830-73853086` (5.32), `chrX:131791829-131793778` (4.33), `chrX:131740917-131742318` (4.15), `chrX:131704058-131705198` (3.88), `chrX:131752284-131753684` (3.67), `chrX:131728490-131729969` (3.47), `chrX:131821916-131822742` (3.13), `chrX:131746492-131747458` (3.05), `chrX:19127-20649` (2.90), `chrX:131762647-131764200` (2.87), the first of them at the `XIST` locus and most of
  the others in one cluster near chrX:131.7–131.8 Mb.

## 4. What keeping chrX changed in the export

The variant is the primary pipeline with chrX present, not the primary blocks with chrX appended. With the
1,148 chrX genes and 5,886 chrX peaks
left in, each donor's library size for the CPM filter includes their counts, so every feature's CPM shifts slightly and the
abundance filter moves at its margin: 20,995 genes pass in the variant
export of the primary cohort against 20,387 in the primary export (the
difference, 608, is
not simply the 684 chrX genes that pass: the larger
library sizes push 76 autosomal genes that passed in the primary export below CPM ≥ 1 in ≥ k donors), and 130,697 peaks against
128,916. TMM factors and log2-CPM are then recomputed on the new filtered
set, and the variance ranking is run with chrX features competing for the 5,000 slots: 125 chrX genes and
102 chrX peaks win a slot in the primary cohort (129 and 110 in
the sensitivity cohort), so exactly that many primary-only features drop out of the top 5,000 and the shared sets are
4,875 / 4,898 (cohort_30) and 4,870 / 4,885
(cohort_50). The chrX genes that enter are dominated by the X-inactivation locus: `XIST` has the largest log2-CPM variance of all
5,000 kept RNA features in both cohorts (about 23.2, female-only expression), `TSIX`
is the next chrX gene, and the top kept chrX peak is at the same locus. These are features whose values differ between the sexes by X dosage, not by stage;
a constant sex offset shifts the two sexes' centroids but not the within-sex displacement between stages, and the PLS basis is
stage-supervised, which is consistent with the near-identical `pls_0` axis and stage means in sections 2 and 3. The variant's
`features_<omic>.csv` labels every chrX feature `abundance`, `variance` or `kept`, never `sex_chromosome`.

## 5. Reproducibility checks

- **Primary export reproduced on this machine.** Before the run, `run_export(cohort="all")` into a scratch directory from the
  synced pseudobulk outputs and the freshly downloaded annotation gave, in both cohorts, all six CSVs byte-identical to the kept
  `results/motco_input/<cohort>/` and equal `feature_counts`, `thresholds` and `removed_by_chromosome`, so the refactored export
  (the variant and the primary now share one code path) and the machine both reproduce the laptop-made export exactly. The kept
  export was not replaced.
- **Second run of the variant.** `preprocessing-sensitivity --cohort cohort_50 --out-dir <scratch> --force` on the same machine, same `n_jobs`, gave the six variant export CSVs, `latent_scores.csv`, `pls_loadings.csv`, `stage_means.csv`, `analysis_metadata.csv` and `sensitivity.csv` byte-identical to the kept run, and equal `statistics`, `p_values`, `null_summary` and `selected_lv`.
- **MOTCO's own CLI.** The `cli_equivalent` of the `cohort_30` variant manifest (`motco analyze` on the variant blocks and `analysis_metadata.csv`, with the six `--covariate` flags, the CV flags, 999 permutations, seed 20261019 and `--n-jobs 1`), run from the repository root into a scratch directory, gave equal `statistics`, `p_values`, `null_summary` and `selected_lv`.
- **Synthetic checks** (before the real run, in the change's tasks): the refactor left the `export` output byte-identical on a
  synthetic pseudobulk set; the reference checks stop the step, writing nothing, when the pseudobulk manifest, the annotation,
  a threshold constant or the donor table's cells differ from what the primary export recorded; the mirrored layout lets
  `rank-sensitivity` run on the variant without code changes; and the rank sweep's own table is unchanged by this change.

## 6. Interpretation

- **The primary analysis stays primary.** The variant is sensitivity to it and not its replacement; it was predeclared in the
  registry and no variant was chosen after seeing the results. `cohort_30` is primary and `cohort_50` sensitivity, as before.
- **No preprocessing dependence at this setting.** No statistic changes rejection status between the primary and the variant
  in either cohort, and no statistic is rejected in either row. Keeping chrX therefore gives no evidence of a difference either,
  and this is not evidence of agreement in size, direction or shape. Had a status differed, it would have been preprocessing
  dependence of the sensitivity analysis, not a result of the variant.
- **`delta`.** Not rejected in either row (0.3 → 0.301 in the primary
  cohort, 0.213 → 0.212 in the sensitivity cohort). A `delta` rejection
  would have been a difference, not specifically of size; the adjusted analysis carries no power claim, so these non-rejections are
  not evidence against a size difference.
- **`angle`.** CV selected rank 1 in the variant as in the primary, so `angle` is 0° (p = 1) by construction in both rows and is
  uninformative about orientation; the ranks agree, so there is no rank difference to report beside it. This non-rejection is not
  evidence that the sexes' directions agree.
- **`shape`.** Not rejected in either row (0.333 → 0.337;
  0.377 → 0.391); uninformative, as in the primary analysis.
- **Descriptively**, the single latent component and the stage means along it are nearly unchanged (section 2 and 3), which says
  that the primary result does not hinge on the chrX removal; it says nothing about whether the result is right.
- The covariate block of the primary analysis (reference row, design columns, missing policy) is unchanged and is reported in
  `reports/motco-astrocyte-2026-10-19.md` §2–3; the variant's `result.json` carries the same six covariate columns in the same order.

## 7. Limitations

- **X-dosage features entered the blocks.** With chrX kept, features whose expression or accessibility differs between the sexes
  by dosage (`XIST`, `TSIX` and the XIST-locus peak first, and possibly X-inactivation escapees among the other kept chrX features) enter
  the latent space. The primary export removed them by design because sex is the group. MOTCO's operating characteristics (Type I
  error, power, the orientation probe) were measured on InterSIM blocks without such features, so the variant carries no Type I or
  power claim of its own, and its p-values are a sensitivity description, not calibrated tests.
- **Rank 1 again.** The variant's latent space has one dimension, like the primary's, so `angle` is degenerate and the rank-1
  limitations of the primary report apply unchanged (rank 1 is outside MOTCO's simulated range).
- **One variant, one axis.** Only the chrX setting was varied, at the primary's `k`, thresholds and N. Other N values and
  residualized features are deferred by decision; a feature-selection variant (trend-corrected variance or a ribosomal /
  mitochondrial exclusion list) was suggested by a diagnostic of the primary export but is not part of this change.
- **No forced-rank sweep on the variant.** The layout allows `rank-sensitivity` to read the variant, but it was not run here.
- **Monte Carlo error and `n_jobs`.** p-values from 999 permutations vary by about 0.015 (SD) between seeds near
  p = 0.2–0.4, so differences of a few thousandths between the rows are within the draw's resolution even though both rows use the
  same permutation indices; p-values depend on `n_jobs` through MOTCO's per-worker seeding, and the reported runs use `n_jobs = 1`.
- **Uncommitted implementation at run time.** The comparison manifests record commit `50bbb5a` with `dirty: true`; the variant
  outputs are deterministic (section 5), so a rerun after the change is committed reproduces them.

## 8. Next steps

- Other preprocessing variants as registry entries: a selection-rule variant (trend-corrected variance or an RP/MT exclusion
  list), other N values and residualized features (the last two deferred by decision on 2026-10-19).
- The forced-rank sweep on the `keep_chrx` variant (`rank-sensitivity --input-dir results/motco_sensitivity/keep_chrx/motco_input
  --analysis-dir results/motco_sensitivity/keep_chrx/motco_analysis --out-dir …`), read as sensitivity.
- MOTCO side: whether CV may select rank 1, where `angle` is degenerate; whether parallel RRPP should reproduce the serial draws.
- Pathway/module enrichment remains exploratory without an orientation result.

## 9. Reproduction

```
scripts/sync_pseudobulk.sh                                   # results/pseudobulk/ from the cluster (no per-donor cache)
uv run sead-sexdiff preprocessing-sensitivity                 # all variants, both cohorts; downloads the GTF once
uv run sead-sexdiff preprocessing-sensitivity --cohort cohort_50 --out-dir /tmp/rerun --force   # reproducibility check
```

The primary export (`results/motco_input/`) and analyses (`results/motco_analysis/`) must be the ones the comparison manifests
hash; the step refuses to run otherwise. Each variant `analysis_manifest.json` carries a `cli_equivalent` that reproduces that
analysis with `motco analyze` on the variant blocks and `analysis_metadata.csv`.
