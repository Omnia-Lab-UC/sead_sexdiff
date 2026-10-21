# Pseudobulk stage DE of the MOTCO astrocyte blocks: ADNC stage and sex × stage — 2026-10-21

> Per-feature DESeq2 of ADNC stage and of the sex × stage interaction on the astrocyte donor pseudobulks of the SEA-AD MTG, on the features of the primary MOTCO export, as a description of the stage signal and a predeclared cross-check of the primary MOTCO analysis (OpenSpec change `mtg-astrocyte-stage-de`, spec `pseudobulk-stage-de`). Nothing is rejected in the MOTCO analysis or in any of its sensitivity analyses, and the latent-adjustment arm left two readings open (`reports/motco-astrocyte-latent-adjustment-2026-10-21.md`): little stage signal in the pseudobulks, or a stage signal that is the same in both sexes. Per cohort and omic, four predeclared covariate-adjusted models (`factor_additive`, `factor_interaction`, `trend_additive`, `trend_interaction`) were fitted with pydeseq2 on the raw counts of the cohort donors, and the predeclared concordance rule compared the exported interaction count of `factor_interaction` with the primary analysis's `delta` rejection status. Every number below is pasted by script from the generated files (`stage_de.md`, `stage_de_summary.csv`, `stage_de_manifest.json`, the fit records and per-coefficient tables, the primary `result.json`, the export manifests, the run log and the logs of the reproducibility and R checks) or, where named, from the verification outputs of the change (synthetic data and simulations), and the tables were checked against `stage_de_summary.csv` and the `result.json` records before commit. The hand-written text carries only the reading, under the rules below.

**Reading rules (binding).** The stage counts describe the stage signal in the pseudobulks. The interaction counts and the concordance cells cross-check the primary analysis; they are not a second test of its hypothesis. A discordance in either direction would be reported as such, since per-feature FDR over thousands of features and one latent-space permutation test ask related, not identical, questions. A zero interaction count is "none detected at this n", never "no interaction". The `union` counts carry no joint error-rate claim. The primary analysis stays primary with its non-rejection rules: `angle` and `shape` non-rejections are uninformative, and a `delta` non-rejection is not evidence that the trajectories agree. The step read sex and stage, so it is a check beside the analysis and never a filter: no feature set, export or analysis is derived from it. `cohort_30` is primary and `cohort_50` sensitivity (a donor subset of `cohort_30`, not an independent replication).

## 1. Provenance

- Software: pydeseq2 0.5.3, formulaic 1.2.1, formulaic-contrasts 1.0.0, statsmodels 0.14.5, numpy 2.3.5, scipy 1.15.3, pandas 2.3.3, anndata 0.12.7, Python 3.12.15 (the Nix shell); MOTCO 0.9.0 installed for the reference check; R DESeq2 1.52.0 (Nix) for the cross-check. sead-sexdiff revision `f8a1683e7288c6fa86939d2c5d41c7c95a1bf378` (dirty: True; the implementation of this change was uncommitted at run time; the reruns of section 4 reproduced the outputs byte for byte).
- Machine: `nixos-dell-laptop`, 2026-10-21. The primary analysis under MOTCO 0.9.0 was pulled from the cluster (`/home1/tgonza/sead_sexdiff/results/motco_analysis/`) after the laptop's kept MOTCO 0.8.0 directories were moved to `results/motco_analysis_0.8.0/`, the cluster's layout; the pulled `result.json` and `analysis_manifest.json` hash as `reports/motco-astrocyte-latent-adjustment-2026-10-21.md` records. All 32 input files hashed identically before and after every run.
- Run: `sead-sexdiff stage-de --cohort all --omic both --n-cpus 10`; fits 496 s (cohort_30) and 517 s (cohort_50). The tables do not depend on `n_cpus` (section 4); the warning records do (below).
- Pseudobulk inputs: `rna_counts.h5ad` sha256 `1253bed9e5f4a762393ecfbb74ae00c550c38cecd648c8a98ea817170885de40`, `atac_counts.h5ad` sha256 `2cfbe937d8e94cb1f3e0024ba0dcc1df0fd731d658b55be62f4e867fb46cff0d`, `donor_table.csv` sha256 `1b34658b2d65132472a6af2d7be7a26f239df305e1f454f289983c89f0a45e71`, `manifest.json` sha256 `2e93cc0b526856fc043f908a424a11266187146636bc75afed7edd187989a25c`.
- cohort_30: 80 donors, cells Female/stage 0 11, Female/stage 1 10, Female/stage 2 28, Male/stage 0 9, Male/stage 1 10, Male/stage 2 12, k = 9; primary export manifest sha256 `18f2289e0d4c84ed18f4a437ae569052544e9d150ed2726385273694b4a3b0b6`; primary analysis `result.json` sha256 `1d6154b02d7704ef9313ce752035a124cf6cb3fbd5350f0428439a542d5190cd` (CV-selected rank 1, MOTCO 0.9.0, seed 20261019, 999 permutations; `delta` rejected: no).
- cohort_50: 74 donors, cells Female/stage 0 10, Female/stage 1 9, Female/stage 2 25, Male/stage 0 9, Male/stage 1 9, Male/stage 2 12, k = 9; primary export manifest sha256 `dd14a6e7d5583fb527e8a99c7db43d12cf237f5451522e8875c5a0a3401c5cb3`; primary analysis `result.json` sha256 `4c6e7e58e39d8ce82c9b93907b70b2110fc2c415e32da40fef8668b74faff9fb` (CV-selected rank 1, MOTCO 0.9.0, seed 20261019, 999 permutations; `delta` rejected: no).
- Reference checks (before any fit): the export manifest records the sha256 of the pseudobulk `manifest.json` on disk; the export records no count-file hash, so the counts were bound to it by reproduction — the export's `preprocess_block` rerun on the cohort donors' counts gave the recorded statuses and per-donor library sizes exactly, and:
  - cohort_30 rna: TMM factors to 0 (relative), every tested feature's log2-CPM variance to 4.95e-10 (relative), exported values to 4.99e-09 (absolute); 20387 features tested, 5000 exported.
  - cohort_30 atac: TMM factors to 0 (relative), every tested feature's log2-CPM variance to 4.98e-10 (relative), exported values to 5e-10 (absolute); 128916 features tested, 5000 exported.
  - cohort_50 rna: TMM factors to 0 (relative), every tested feature's log2-CPM variance to 4.96e-10 (relative), exported values to 4.99e-09 (absolute); 20238 features tested, 5000 exported.
  - cohort_50 atac: TMM factors to 0 (relative), every tested feature's log2-CPM variance to 4.99e-10 (relative), exported values to 5e-10 (absolute); 124818 features tested, 5000 exported.
  The primary analysis matched the export (`_check_reference`: the five export files' sha256, MOTCO version, seed, permutations, CV selection); donors, sex and stage agreed between the export, the donor table and the count files; every model's design matrix was full rank and every tested column resolved by name.
- Models (`STAGE_DE_MODELS`), each followed by the omic's covariates `age_at_death`, `pmi`, `rin`, `years_education` and `log10_n_rna` (RNA) or `log10_n_atac` (ATAC), continuous, on their recorded scale; sex reference `Female`, stage reference `0`:
  - `factor_additive`: `~ sex + stage`; tested `stage_1_vs_0` (`stage[T.1]`), `stage_2_vs_0` (`stage[T.2]`).
  - `factor_interaction`: `~ sex + stage + sex:stage`; tested `male_x_stage_1` (`sex[T.Male]:stage[T.1]`), `male_x_stage_2` (`sex[T.Male]:stage[T.2]`).
  - `trend_additive`: `~ sex + stage_num`; tested `stage_num` (`stage_num`).
  - `trend_interaction`: `~ sex + stage_num + sex:stage_num`; tested `male_x_stage_num` (`sex[T.Male]:stage_num`).
- Significance: padj ≤ 0.05, pydeseq2's Benjamini–Hochberg over the tested features that pass its independent filter (a baseMean cutoff chosen per coefficient; a filtered feature has no padj and counts as not significant), no permutation null, no LFC shrinkage. Tested features: those of `features_<omic>.csv` with status `variance` or `kept` (the abundance-passing set; size factors and dispersions fitted on it); `exported` = status `kept`, MOTCO's 5,000 features per block, the primary reading.
- Fits:
  - cohort_30 rna: dispersion trend parametric; features not converged per fit gene-wise 5–12, MAP 7–14, LFC 18–28; Cook's replaceable donors, refitted and outlier features summed over the four fits 0 / 0 / 0; warnings captured in the main process 0 (n_cpus 10).
  - cohort_30 atac: dispersion trend parametric; features not converged per fit gene-wise 10–18, MAP 21–29, LFC 2–6; Cook's replaceable donors, refitted and outlier features summed over the four fits 0 / 0 / 0; warnings captured in the main process 0 (n_cpus 10).
  - cohort_50 rna: dispersion trend parametric; features not converged per fit gene-wise 7–13, MAP 10–18, LFC 17–27; Cook's replaceable donors, refitted and outlier features summed over the four fits 0 / 0 / 0; warnings captured in the main process 0 (n_cpus 10).
  - cohort_50 atac: dispersion trend parametric; features not converged per fit gene-wise 7–15, MAP 14–24, LFC 2–6; Cook's replaceable donors, refitted and outlier features summed over the four fits 0 / 0 / 0; warnings captured in the main process 0 (n_cpus 10).
  With n_cpus 10 the fit records capture no warning raised in pydeseq2's worker processes. The run log holds 50 such RuntimeWarnings (10 × 'invalid value encountered in divide', 10 × 'invalid value encountered in multiply', 10 × 'invalid value encountered in subtract', 20 × 'overflow encountered in exp'), raised in pydeseq2's IRLS fallback (`utils.py`); in the n_cpus = 1 rerun of task 3.3 every cohort_30 RNA fit recorded 4 (invalid value encountered in divide; invalid value encountered in multiply; invalid value encountered in subtract; overflow encountered in exp). They go with the 2–28 features per fit whose LFC fit did not converge; the tables are byte-identical at n_cpus 1, 8 and 10.

## 2. Primary cohort: cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80)

### Concordance (pasted from `results/stage_de/cohort_30/stage_de.md`)

Rule (`STAGE_DE_CONCORDANCE_RULE`, predeclared): `interaction_detected` = at least one exported feature (status kept) with padj <= STAGE_DE_ALPHA on either interaction coefficient of factor_interaction; `delta_rejected` = the primary analysis rejects delta (p <= ANALYSIS_ALPHA, ANALYSIS_REJECTION_RULE; result.json); concordant when both are true or both are false; discordant otherwise, written as a finding. The primary analysis (`results/motco_analysis/cohort_30`, CV-selected rank 1, MOTCO 0.9.0) is read, never rerun or modified; a DE result changes no MOTCO reading.

| omic | interaction detected | exported features with interaction padj ≤ α | delta rejected | cell |
|---|:---:|---:|:---:|---|
| rna | no | 0 | no | concordant |
| atac | no | 0 | no | concordant |

### Counts: rna (pasted from `results/stage_de/cohort_30/stage_de.md`)

5000 exported and 20387 abundance-passing features tested; significant / tested at padj ≤ α, missing padj (removed by the independent filter; for `union`, missing on both coefficients) and the smallest padj per set.

| model | coefficient | design column | exported: significant / tested | abundance: significant / tested | padj missing (exported / abundance) | min padj (exported) | min padj (abundance) |
|---|---|---|---:|---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | `stage[T.1]` | 0 / 5000 | 0 / 20387 | 0 / 0 | 0.526 | 0.526 |
| factor_additive | stage_2_vs_0 | `stage[T.2]` | 1 / 5000 | 1 / 20387 | 0 / 0 | 7.17e-06 | 7.17e-06 |
| factor_additive | union (no joint error rate) | `stage[T.1];stage[T.2]` | 1 / 5000 | 1 / 20387 | 0 / 0 | 7.17e-06 | 7.17e-06 |
| factor_interaction | male_x_stage_1 | `sex[T.Male]:stage[T.1]` | 0 / 5000 | 0 / 20387 | 0 / 0 | 0.095 | 0.095 |
| factor_interaction | male_x_stage_2 | `sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 20387 | 0 / 0 | 0.131 | 0.131 |
| factor_interaction | union (no joint error rate) | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 20387 | 0 / 0 | 0.095 | 0.095 |
| trend_additive | stage_num | `stage_num` | 2 / 5000 | 2 / 20387 | 0 / 0 | 4.34e-06 | 4.34e-06 |
| trend_interaction | male_x_stage_num | `sex[T.Male]:stage_num` | 0 / 5000 | 0 / 20387 | 0 / 0 | 0.108 | 0.108 |

### Counts: atac (pasted from `results/stage_de/cohort_30/stage_de.md`)

5000 exported and 128916 abundance-passing features tested; significant / tested at padj ≤ α, missing padj (removed by the independent filter; for `union`, missing on both coefficients) and the smallest padj per set.

| model | coefficient | design column | exported: significant / tested | abundance: significant / tested | padj missing (exported / abundance) | min padj (exported) | min padj (abundance) |
|---|---|---|---:|---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | `stage[T.1]` | 0 / 5000 | 0 / 128916 | 0 / 0 | 0.602 | 0.283 |
| factor_additive | stage_2_vs_0 | `stage[T.2]` | 0 / 5000 | 22 / 128916 | 202 / 29993 | 0.0581 | 0.0257 |
| factor_additive | union (no joint error rate) | `stage[T.1];stage[T.2]` | 0 / 5000 | 22 / 128916 | 0 / 0 | 0.0581 | 0.0257 |
| factor_interaction | male_x_stage_1 | `sex[T.Male]:stage[T.1]` | 0 / 5000 | 0 / 128916 | 0 / 0 | 0.463 | 0.463 |
| factor_interaction | male_x_stage_2 | `sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 128916 | 0 / 0 | 0.421 | 0.421 |
| factor_interaction | union (no joint error rate) | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 128916 | 0 / 0 | 0.421 | 0.421 |
| trend_additive | stage_num | `stage_num` | 0 / 5000 | 23 / 128916 | 318 / 32492 | 0.0571 | 0.0185 |
| trend_interaction | male_x_stage_num | `sex[T.Male]:stage_num` | 0 / 5000 | 0 / 128916 | 0 / 0 | 0.845 | 0.595 |

### Features at padj ≤ 0.05 (generated from the per-coefficient tables; descriptive)

Every feature at padj ≤ 0.05 per tested coefficient, by padj (ties by p-value, then name); a coefficient without one lists its five smallest padj (BH ties many features at the same padj, so these are the five smallest p-values among them, not top interaction features). `exported` = one of MOTCO's 5,000 features.

- rna `factor_additive` `stage_1_vs_0`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `LINC02150` | 10.3 | -1.202 | 0.286 | 0.526 | yes |
  | `STAC3` | 20.7 | -0.729 | 0.186 | 0.613 | no |
  | `SELENOV` | 5.5 | -1.264 | 0.324 | 0.613 | yes |
  | `TMEM176A` | 30.6 | 0.639 | 0.166 | 0.613 | no |
  | `MTRNR2L8` | 52.6 | 1.851 | 0.494 | 0.727 | yes |

- rna `factor_additive` `stage_2_vs_0`: 1 at padj ≤ 0.05 (1 exported):

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `CREB5` | 256.3 | 1.210 | 0.193 | 7.17e-06 | yes |

- rna `factor_interaction` `male_x_stage_1`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `ROBO2` | 1000.2 | 2.104 | 0.459 | 0.095 | yes |
  | `AC007326.5` | 15.4 | 2.635 | 0.666 | 0.707 | yes |
  | `HVCN1` | 302.3 | 1.446 | 0.373 | 0.707 | no |
  | `EGF` | 248.3 | 2.392 | 0.654 | 1 | yes |
  | `LINC02475` | 2.4 | 4.811 | 1.356 | 1 | yes |

- rna `factor_interaction` `male_x_stage_2`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `AC006305.1` | 28.3 | 3.564 | 0.790 | 0.131 | yes |
  | `AC133785.1` | 16.3 | 1.686 | 0.387 | 0.135 | yes |
  | `AC104083.1` | 9.4 | 1.952 | 0.495 | 0.461 | yes |
  | `AC091074.2` | 8.4 | 1.961 | 0.515 | 0.461 | yes |
  | `ITPR3` | 7.6 | -1.480 | 0.393 | 0.461 | no |

- rna `trend_additive` `stage_num`: 2 at padj ≤ 0.05 (2 exported):

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `CREB5` | 256.3 | 0.606 | 0.095 | 4.34e-06 | yes |
  | `ADAM33` | 172.8 | 0.618 | 0.135 | 0.0484 | yes |

- rna `trend_interaction` `male_x_stage_num`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `AC006305.1` | 28.3 | 1.803 | 0.401 | 0.108 | yes |
  | `AC133785.1` | 16.3 | 0.829 | 0.188 | 0.108 | yes |
  | `AC091074.2` | 8.4 | 1.041 | 0.257 | 0.349 | yes |
  | `AC104083.1` | 9.4 | 0.929 | 0.243 | 0.569 | yes |
  | `CDH12` | 292.2 | 0.574 | 0.154 | 0.569 | no |

- atac `factor_additive` `stage_1_vs_0`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr5:10592071-10592905` | 2.4 | 2.584 | 0.561 | 0.283 | no |
  | `chr1:26809972-26810485` | 13.8 | 0.862 | 0.188 | 0.283 | no |
  | `chr8:119152756-119153569` | 4.4 | 1.759 | 0.407 | 0.602 | no |
  | `chr22:35154787-35155489` | 6.3 | 1.297 | 0.308 | 0.602 | no |
  | `chr3:9848511-9849427` | 22.5 | 0.600 | 0.143 | 0.602 | no |

- atac `factor_additive` `stage_2_vs_0`: 22 at padj ≤ 0.05 (0 exported):

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr19:12647443-12648424` | 18.5 | 1.098 | 0.213 | 0.0257 | no |
  | `chr8:119152756-119153569` | 4.4 | 1.826 | 0.378 | 0.0402 | no |
  | `chr17:76786356-76788339` | 6.8 | 1.567 | 0.327 | 0.0402 | no |
  | `chr1:172321371-172323154` | 23.6 | 0.895 | 0.189 | 0.0402 | no |
  | `chr19:46755222-46758121` | 50.7 | 0.439 | 0.093 | 0.0402 | no |
  | `chr7:44037666-44039504` | 36.0 | 0.594 | 0.126 | 0.0402 | no |
  | `chr1:156158129-156164411` | 71.4 | 0.501 | 0.109 | 0.0448 | no |
  | `chr15:74791714-74793922` | 35.0 | 0.604 | 0.132 | 0.0448 | no |
  | `chr1:26809972-26810485` | 13.8 | 0.787 | 0.173 | 0.0448 | no |
  | `chr2:131671819-131674704` | 58.5 | 0.467 | 0.103 | 0.0448 | no |
  | `chr1:12166811-12167277` | 9.9 | 1.033 | 0.228 | 0.0448 | no |
  | `chr17:78427475-78428426` | 6.0 | 1.265 | 0.279 | 0.0448 | no |
  | `chr3:194880040-194881125` | 85.7 | -0.471 | 0.104 | 0.0448 | no |
  | `chr12:106328829-106330471` | 68.6 | -0.412 | 0.092 | 0.0448 | no |
  | `chr4:173360734-173361293` | 11.7 | 0.939 | 0.210 | 0.0448 | no |
  | `chr14:49791015-49792973` | 8.2 | 1.133 | 0.254 | 0.0448 | no |
  | `chr8:116089472-116089647` | 12.1 | -0.775 | 0.173 | 0.0448 | no |
  | `chr1:24201835-24202481` | 33.3 | -0.726 | 0.163 | 0.0448 | no |
  | `chr7:44143794-44148504` | 120.7 | 0.439 | 0.099 | 0.0473 | no |
  | `chr8:134480391-134482291` | 58.1 | 0.553 | 0.125 | 0.0489 | no |
  | `chr3:143134550-143135514` | 38.9 | -0.674 | 0.153 | 0.0489 | no |
  | `chr16:69772503-69773933` | 13.6 | 0.904 | 0.206 | 0.0489 | no |

- atac `factor_interaction` `male_x_stage_1`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr7:30924156-30925045` | 8.0 | -3.550 | 0.766 | 0.463 | yes |
  | `chr5:116738352-116738721` | 13.5 | -2.051 | 0.485 | 0.938 | no |
  | `chr8:144802240-144811386` | 78.2 | 1.121 | 0.267 | 0.938 | no |
  | `chr12:48753307-48754450` | 14.4 | 1.571 | 0.379 | 0.938 | no |
  | `chr9:132977315-132978071` | 31.2 | 1.212 | 0.294 | 0.938 | no |

- atac `factor_interaction` `male_x_stage_2`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr7:30924156-30925045` | 8.0 | -3.040 | 0.684 | 0.421 | yes |
  | `chr4:153212080-153213297` | 22.3 | -1.418 | 0.321 | 0.421 | no |
  | `chr2:232307884-232311263` | 39.2 | -1.608 | 0.366 | 0.421 | no |
  | `chr3:160808131-160812952` | 138.2 | -0.657 | 0.151 | 0.421 | no |
  | `chr10:29589886-29591981` | 13.6 | -1.767 | 0.418 | 0.517 | no |

- atac `trend_additive` `stage_num`: 23 at padj ≤ 0.05 (0 exported):

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr19:12647443-12648424` | 18.5 | 0.543 | 0.104 | 0.0185 | no |
  | `chr17:76786356-76788339` | 6.8 | 0.772 | 0.158 | 0.0361 | no |
  | `chr7:44037666-44039504` | 36.0 | 0.301 | 0.062 | 0.0361 | no |
  | `chr19:46755222-46758121` | 50.7 | 0.219 | 0.046 | 0.0365 | no |
  | `chr1:12166811-12167277` | 9.9 | 0.510 | 0.110 | 0.0498 | no |
  | `chr7:44143794-44148504` | 120.7 | 0.224 | 0.049 | 0.0498 | no |
  | `chr1:172321371-172323154` | 23.6 | 0.427 | 0.093 | 0.0498 | no |
  | `chr17:78427475-78428426` | 6.0 | 0.611 | 0.134 | 0.0498 | no |
  | `chr3:194880040-194881125` | 85.7 | -0.234 | 0.052 | 0.0498 | no |
  | `chr1:207890323-207890716` | 14.8 | 0.460 | 0.102 | 0.0498 | no |
  | `chr8:119152756-119153569` | 4.4 | 0.793 | 0.177 | 0.0498 | no |
  | `chr12:106328829-106330471` | 68.6 | -0.204 | 0.046 | 0.0498 | no |
  | `chr8:134480391-134482291` | 58.1 | 0.276 | 0.062 | 0.0498 | no |
  | `chr8:139650153-139659350` | 42.3 | 0.298 | 0.067 | 0.0498 | no |
  | `chr5:127996542-127998934` | 23.2 | 0.391 | 0.088 | 0.0498 | no |
  | `chr2:131671819-131674704` | 58.5 | 0.225 | 0.051 | 0.0498 | no |
  | `chr14:49791015-49792973` | 8.2 | 0.543 | 0.122 | 0.0498 | no |
  | `chr1:156158129-156164411` | 71.4 | 0.238 | 0.054 | 0.0498 | no |
  | `chr21:45833370-45834589` | 69.0 | -0.283 | 0.064 | 0.0498 | no |
  | `chr15:74791714-74793922` | 35.0 | 0.284 | 0.065 | 0.0498 | no |
  | `chr12:75378102-75379687` | 8.0 | 0.597 | 0.136 | 0.0498 | no |
  | `chr8:116089472-116089647` | 12.1 | -0.381 | 0.087 | 0.0498 | no |
  | `chr1:24201835-24202481` | 33.3 | -0.354 | 0.081 | 0.0498 | no |

- atac `trend_interaction` `male_x_stage_num`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr5:142674265-142675471` | 51.9 | -0.568 | 0.129 | 0.595 | no |
  | `chr9:87563421-87565541` | 29.7 | -0.624 | 0.145 | 0.595 | no |
  | `chr14:64701693-64702058` | 8.1 | -1.169 | 0.271 | 0.595 | no |
  | `chr2:232307884-232311263` | 39.2 | -0.754 | 0.177 | 0.595 | no |
  | `chr3:160808131-160812952` | 138.2 | -0.310 | 0.073 | 0.595 | no |

### Findings (pasted from `results/stage_de/cohort_30/stage_de.md`)

- fits made with n_cpus != 1: rna/factor_additive (10), rna/factor_interaction (10), rna/trend_additive (10), rna/trend_interaction (10), atac/factor_additive (10), atac/factor_interaction (10), atac/trend_additive (10), atac/trend_interaction (10); pydeseq2's per-feature fits are expected not to depend on n_cpus (the report states whether this was checked), but warnings raised in pydeseq2's worker processes were not captured for these fits, so their warning records and findings are incomplete

## 3. Sensitivity cohort: cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74)

### Concordance (pasted from `results/stage_de/cohort_50/stage_de.md`)

Rule (`STAGE_DE_CONCORDANCE_RULE`, predeclared): `interaction_detected` = at least one exported feature (status kept) with padj <= STAGE_DE_ALPHA on either interaction coefficient of factor_interaction; `delta_rejected` = the primary analysis rejects delta (p <= ANALYSIS_ALPHA, ANALYSIS_REJECTION_RULE; result.json); concordant when both are true or both are false; discordant otherwise, written as a finding. The primary analysis (`results/motco_analysis/cohort_50`, CV-selected rank 1, MOTCO 0.9.0) is read, never rerun or modified; a DE result changes no MOTCO reading.

| omic | interaction detected | exported features with interaction padj ≤ α | delta rejected | cell |
|---|:---:|---:|:---:|---|
| rna | no | 0 | no | concordant |
| atac | no | 0 | no | concordant |

### Counts: rna (pasted from `results/stage_de/cohort_50/stage_de.md`)

5000 exported and 20238 abundance-passing features tested; significant / tested at padj ≤ α, missing padj (removed by the independent filter; for `union`, missing on both coefficients) and the smallest padj per set.

| model | coefficient | design column | exported: significant / tested | abundance: significant / tested | padj missing (exported / abundance) | min padj (exported) | min padj (abundance) |
|---|---|---|---:|---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | `stage[T.1]` | 0 / 5000 | 0 / 20238 | 0 / 0 | 0.229 | 0.12 |
| factor_additive | stage_2_vs_0 | `stage[T.2]` | 2 / 5000 | 2 / 20238 | 0 / 0 | 3.21e-05 | 3.21e-05 |
| factor_additive | union (no joint error rate) | `stage[T.1];stage[T.2]` | 2 / 5000 | 2 / 20238 | 0 / 0 | 3.21e-05 | 3.21e-05 |
| factor_interaction | male_x_stage_1 | `sex[T.Male]:stage[T.1]` | 0 / 5000 | 0 / 20238 | 0 / 0 | 0.414 | 0.414 |
| factor_interaction | male_x_stage_2 | `sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 20238 | 0 / 0 | 0.0681 | 0.0681 |
| factor_interaction | union (no joint error rate) | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 20238 | 0 / 0 | 0.0681 | 0.0681 |
| trend_additive | stage_num | `stage_num` | 2 / 5000 | 2 / 20238 | 0 / 0 | 2.01e-05 | 2.01e-05 |
| trend_interaction | male_x_stage_num | `sex[T.Male]:stage_num` | 0 / 5000 | 0 / 20238 | 0 / 0 | 0.0827 | 0.0827 |

### Counts: atac (pasted from `results/stage_de/cohort_50/stage_de.md`)

5000 exported and 124818 abundance-passing features tested; significant / tested at padj ≤ α, missing padj (removed by the independent filter; for `union`, missing on both coefficients) and the smallest padj per set.

| model | coefficient | design column | exported: significant / tested | abundance: significant / tested | padj missing (exported / abundance) | min padj (exported) | min padj (abundance) |
|---|---|---|---:|---:|---:|---:|---:|
| factor_additive | stage_1_vs_0 | `stage[T.1]` | 0 / 5000 | 0 / 124818 | 0 / 0 | 0.964 | 0.385 |
| factor_additive | stage_2_vs_0 | `stage[T.2]` | 4 / 5000 | 36 / 124818 | 1471 / 38719 | 0.0295 | 0.0139 |
| factor_additive | union (no joint error rate) | `stage[T.1];stage[T.2]` | 4 / 5000 | 36 / 124818 | 0 / 0 | 0.0295 | 0.0139 |
| factor_interaction | male_x_stage_1 | `sex[T.Male]:stage[T.1]` | 0 / 5000 | 0 / 124818 | 0 / 0 | 0.714 | 0.714 |
| factor_interaction | male_x_stage_2 | `sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 124818 | 0 / 0 | 0.46 | 0.392 |
| factor_interaction | union (no joint error rate) | `sex[T.Male]:stage[T.1];sex[T.Male]:stage[T.2]` | 0 / 5000 | 0 / 124818 | 0 / 0 | 0.46 | 0.392 |
| trend_additive | stage_num | `stage_num` | 4 / 5000 | 30 / 124818 | 0 / 0 | 0.0264 | 0.0112 |
| trend_interaction | male_x_stage_num | `sex[T.Male]:stage_num` | 0 / 5000 | 0 / 124818 | 0 / 0 | 0.628 | 0.554 |

### Features at padj ≤ 0.05 (generated from the per-coefficient tables; descriptive)

Every feature at padj ≤ 0.05 per tested coefficient, by padj (ties by p-value, then name); a coefficient without one lists its five smallest padj (BH ties many features at the same padj, so these are the five smallest p-values among them, not top interaction features). `exported` = one of MOTCO's 5,000 features.

- rna `factor_additive` `stage_1_vs_0`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `STAC3` | 20.9 | -0.824 | 0.182 | 0.12 | no |
  | `GINS2` | 13.1 | -0.719 | 0.168 | 0.197 | no |
  | `SELENOV` | 5.5 | -1.430 | 0.347 | 0.229 | yes |
  | `LINC02150` | 10.6 | -1.245 | 0.305 | 0.229 | yes |
  | `MTRNR2L8` | 56.2 | 1.972 | 0.514 | 0.437 | yes |

- rna `factor_additive` `stage_2_vs_0`: 2 at padj ≤ 0.05 (2 exported):

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `CREB5` | 265.7 | 1.225 | 0.203 | 3.21e-05 | yes |
  | `ADAM33` | 172.2 | 1.336 | 0.273 | 0.00973 | yes |

- rna `factor_interaction` `male_x_stage_1`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `ROBO2` | 1000.1 | 2.111 | 0.496 | 0.414 | yes |
  | `AC068408.2` | 4.7 | 3.979 | 1.034 | 0.999 | yes |
  | `AC007326.5` | 15.5 | 2.623 | 0.721 | 0.999 | yes |
  | `HVCN1` | 308.8 | 1.348 | 0.404 | 0.999 | no |
  | `EGF` | 251.6 | 2.348 | 0.708 | 0.999 | yes |

- rna `factor_interaction` `male_x_stage_2`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `AC006305.1` | 29.9 | 3.857 | 0.830 | 0.0681 | yes |
  | `AC133785.1` | 16.2 | 1.710 | 0.405 | 0.24 | yes |
  | `DNMT1` | 142.5 | -0.610 | 0.155 | 0.569 | no |
  | `CCDC175` | 41.2 | 1.799 | 0.483 | 0.653 | yes |
  | `AC104083.1` | 9.1 | 1.951 | 0.524 | 0.653 | yes |

- rna `trend_additive` `stage_num`: 2 at padj ≤ 0.05 (2 exported):

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `CREB5` | 265.7 | 0.613 | 0.100 | 2.01e-05 | yes |
  | `ADAM33` | 172.2 | 0.672 | 0.135 | 0.00651 | yes |

- rna `trend_interaction` `male_x_stage_num`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `AC006305.1` | 29.9 | 1.934 | 0.420 | 0.0827 | yes |
  | `AC133785.1` | 16.2 | 0.857 | 0.195 | 0.117 | yes |
  | `LINC01605` | 39.6 | 1.507 | 0.380 | 0.447 | yes |
  | `AC091074.2` | 8.4 | 1.079 | 0.275 | 0.447 | yes |
  | `DNMT1` | 142.5 | -0.286 | 0.075 | 0.567 | no |

- atac `factor_additive` `stage_1_vs_0`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr1:26809972-26810485` | 14.7 | 0.904 | 0.194 | 0.385 | no |
  | `chr14:59376737-59378608` | 51.4 | -0.827 | 0.188 | 0.63 | no |
  | `chr1:245589383-245590226` | 27.3 | -0.812 | 0.189 | 0.63 | no |
  | `chr1:178241119-178243592` | 3.2 | 2.525 | 0.605 | 0.63 | no |
  | `chr4:173360734-173361293` | 12.8 | 0.973 | 0.234 | 0.63 | no |

- atac `factor_additive` `stage_2_vs_0`: 36 at padj ≤ 0.05 (4 exported):

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr19:12647443-12648424` | 19.9 | 1.125 | 0.217 | 0.0139 | no |
  | `chr7:44037666-44039504` | 39.4 | 0.648 | 0.127 | 0.0139 | no |
  | `chr21:45833370-45834589` | 77.2 | -0.601 | 0.125 | 0.0295 | no |
  | `chr14:49791015-49792973` | 9.0 | 1.236 | 0.258 | 0.0295 | no |
  | `chr5:127996542-127998934` | 25.1 | 0.842 | 0.178 | 0.0295 | no |
  | `chr17:76786356-76788339` | 7.2 | 1.595 | 0.339 | 0.0295 | yes |
  | `chr1:156158129-156164411` | 76.9 | 0.508 | 0.108 | 0.0295 | no |
  | `chr8:134480391-134482291` | 63.4 | 0.595 | 0.127 | 0.0295 | no |
  | `chr2:131671819-131674704` | 63.8 | 0.485 | 0.104 | 0.031 | no |
  | `chr19:46755222-46758121` | 55.8 | 0.428 | 0.093 | 0.031 | no |
  | `chr1:172321371-172323154` | 25.2 | 0.865 | 0.188 | 0.031 | no |
  | `chr1:26809972-26810485` | 14.7 | 0.813 | 0.177 | 0.031 | no |
  | `chr12:106328829-106330471` | 74.8 | -0.429 | 0.094 | 0.033 | no |
  | `chr1:200749095-200750613` | 23.4 | -0.739 | 0.163 | 0.0341 | no |
  | `chr1:4652099-4677045` | 133.7 | 0.583 | 0.129 | 0.0348 | no |
  | `chr3:143134550-143135514` | 42.9 | -0.701 | 0.156 | 0.0354 | no |
  | `chr16:69772503-69773933` | 15.2 | 0.908 | 0.203 | 0.0354 | no |
  | `chr1:12166811-12167277` | 10.4 | 1.038 | 0.233 | 0.0354 | no |
  | `chr4:173360734-173361293` | 12.8 | 0.934 | 0.210 | 0.0354 | no |
  | `chr14:35336521-35337063` | 43.8 | 0.679 | 0.153 | 0.0354 | no |
  | `chr1:207890323-207890716` | 16.1 | 0.954 | 0.215 | 0.0354 | no |
  | `chr10:86865478-86866802` | 52.5 | -0.428 | 0.097 | 0.0354 | no |
  | `chr12:75378102-75379687` | 8.6 | 1.239 | 0.280 | 0.0354 | no |
  | `chr6:11605152-11608407` | 54.3 | -0.594 | 0.135 | 0.0383 | no |
  | `chr5:14780301-14782149` | 12.2 | 1.581 | 0.361 | 0.0412 | yes |
  | `chr8:116089472-116089647` | 13.3 | -0.765 | 0.176 | 0.0455 | no |
  | `chr9:109881578-109891357` | 38.0 | -0.661 | 0.153 | 0.0473 | no |
  | `chr2:58427385-58433920` | 23.9 | 0.926 | 0.215 | 0.0486 | no |
  | `chr17:51174146-51176074` | 6.9 | 1.494 | 0.348 | 0.0486 | yes |
  | `chr1:153539635-153542451` | 14.5 | 0.992 | 0.232 | 0.0486 | no |
  | `chr8:139650153-139659350` | 46.6 | 0.599 | 0.140 | 0.0486 | no |
  | `chr3:194880040-194881125` | 92.3 | -0.454 | 0.106 | 0.0486 | no |
  | `chr7:44143794-44148504` | 131.7 | 0.435 | 0.102 | 0.0486 | no |
  | `chr17:78427475-78428426` | 6.3 | 1.226 | 0.288 | 0.0486 | no |
  | `chr6:166488142-166489979` | 10.1 | 1.152 | 0.270 | 0.0486 | yes |
  | `chr1:231053444-231054140` | 24.5 | -0.668 | 0.157 | 0.0486 | no |

- atac `factor_interaction` `male_x_stage_1`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr14:95544578-95546242` | 5.5 | 3.689 | 0.836 | 0.714 | no |
  | `chr7:90932085-90932904` | 3.9 | -3.504 | 0.823 | 0.714 | no |
  | `chr7:30924156-30925045` | 8.8 | -3.308 | 0.778 | 0.714 | yes |
  | `chr5:16664584-16665195` | 6.1 | -2.691 | 0.641 | 0.714 | no |
  | `chr11:19972885-19973951` | 5.9 | -2.605 | 0.626 | 0.714 | no |

- atac `factor_interaction` `male_x_stage_2`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr4:153212080-153213297` | 24.5 | -1.438 | 0.323 | 0.392 | no |
  | `chr5:16664584-16665195` | 6.1 | -2.601 | 0.590 | 0.392 | no |
  | `chr1:145954696-145956095` | 31.1 | -1.334 | 0.304 | 0.392 | no |
  | `chr2:54634115-54635891` | 20.3 | -1.669 | 0.382 | 0.392 | no |
  | `chr7:30924156-30925045` | 8.8 | -2.957 | 0.690 | 0.46 | yes |

- atac `trend_additive` `stage_num`: 30 at padj ≤ 0.05 (4 exported):

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr7:44037666-44039504` | 39.4 | 0.326 | 0.062 | 0.0112 | no |
  | `chr19:12647443-12648424` | 19.9 | 0.554 | 0.106 | 0.0112 | no |
  | `chr21:45833370-45834589` | 77.2 | -0.300 | 0.061 | 0.0264 | no |
  | `chr1:178241119-178243592` | 3.2 | 1.172 | 0.240 | 0.0264 | no |
  | `chr5:127996542-127998934` | 25.1 | 0.425 | 0.087 | 0.0264 | no |
  | `chr17:76786356-76788339` | 7.2 | 0.794 | 0.164 | 0.0264 | yes |
  | `chr14:49791015-49792973` | 9.0 | 0.596 | 0.124 | 0.0267 | no |
  | `chr8:119152756-119153569` | 4.7 | 0.866 | 0.183 | 0.0329 | no |
  | `chr8:134480391-134482291` | 63.4 | 0.293 | 0.062 | 0.0329 | no |
  | `chr19:46755222-46758121` | 55.8 | 0.213 | 0.045 | 0.0329 | no |
  | `chr12:106328829-106330471` | 74.8 | -0.214 | 0.046 | 0.0407 | no |
  | `chr1:207890323-207890716` | 16.1 | 0.481 | 0.105 | 0.0407 | no |
  | `chr1:12166811-12167277` | 10.4 | 0.515 | 0.113 | 0.0407 | no |
  | `chr14:35336521-35337063` | 43.8 | 0.342 | 0.075 | 0.0407 | no |
  | `chr2:131671819-131674704` | 63.8 | 0.234 | 0.051 | 0.0407 | no |
  | `chr1:156158129-156164411` | 76.9 | 0.242 | 0.053 | 0.0407 | no |
  | `chr2:165142365-165143845` | 4.8 | 0.836 | 0.184 | 0.0409 | yes |
  | `chr1:4652099-4677045` | 133.7 | 0.287 | 0.063 | 0.0415 | no |
  | `chr12:75378102-75379687` | 8.6 | 0.607 | 0.135 | 0.0481 | no |
  | `chr1:172321371-172323154` | 25.2 | 0.414 | 0.093 | 0.0481 | no |
  | `chr19:39150443-39151190` | 2.7 | 1.007 | 0.226 | 0.0481 | no |
  | `chr8:139650153-139659350` | 46.6 | 0.305 | 0.069 | 0.0481 | no |
  | `chr5:14780301-14782149` | 12.2 | 0.781 | 0.176 | 0.0481 | yes |
  | `chr9:109881578-109891357` | 38.0 | -0.333 | 0.075 | 0.0481 | no |
  | `chr6:11605152-11608407` | 54.3 | -0.294 | 0.067 | 0.0481 | no |
  | `chr7:44143794-44148504` | 131.7 | 0.220 | 0.050 | 0.0481 | no |
  | `chr10:45547979-45548812` | 5.9 | 0.843 | 0.191 | 0.0481 | yes |
  | `chr1:108984064-108985617` | 5.7 | 0.759 | 0.172 | 0.0481 | no |
  | `chr2:58427385-58433920` | 23.9 | 0.461 | 0.105 | 0.0481 | no |
  | `chr1:153539635-153542451` | 14.5 | 0.495 | 0.113 | 0.0481 | no |

- atac `trend_interaction` `male_x_stage_num`: none at padj ≤ 0.05; smallest five:

  | feature | baseMean | log2FC | lfcSE | padj | exported |
  |---|---:|---:|---:|---:|:---:|
  | `chr2:54634115-54635891` | 20.3 | -0.839 | 0.183 | 0.554 | no |
  | `chr1:145954696-145956095` | 31.1 | -0.640 | 0.145 | 0.575 | no |
  | `chr14:64701693-64702058` | 8.8 | -1.193 | 0.275 | 0.575 | no |
  | `chr9:87563421-87565541` | 32.5 | -0.634 | 0.148 | 0.575 | no |
  | `chr2:114535868-114539024` | 26.5 | -0.760 | 0.183 | 0.575 | no |

### Findings (pasted from `results/stage_de/cohort_50/stage_de.md`)

- fits made with n_cpus != 1: rna/factor_additive (10), rna/factor_interaction (10), rna/trend_additive (10), rna/trend_interaction (10), atac/factor_additive (10), atac/factor_interaction (10), atac/trend_additive (10), atac/trend_interaction (10); pydeseq2's per-feature fits are expected not to depend on n_cpus (the report states whether this was checked), but warnings raised in pydeseq2's worker processes were not captured for these fits, so their warning records and findings are incomplete

## 4. Reproducibility and the R DESeq2 cross-check

- Rerun (task 3.2): `run_stage_de("cohort_30", omic="both", force=True, n_cpus=8)` into a scratch directory (418 s) gave a `stage_de_summary.csv` and 12 tables byte-identical to the reported run's (`n_cpus = 10`).
- `n_cpus` (task 3.3): `stage-de --cohort cohort_30 --omic rna --n-cpus 1` (185 s) gave 6 tables and 16 summary rows byte-identical to the reported run's; on the synthetic data `n_cpus` 1, 2 and 4 gave byte-identical tables for all 16 fits. The manifests' finding on fits made with `n_cpus ≠ 1` stays as a record.
- R DESeq2 1.52.0 (task 3.4), `scripts/check_stage_de_deseq2.R` on the debug files of the 3.2 rerun, all four RNA models of cohort_30. Two comparisons per coefficient: `full`, a `DESeq()` fit with the dispersion trend type pydeseq2 used (parametric here), and `plug-in`, R's Wald test on pydeseq2's size factors and dispersions. The 1e-3 relative (+ 1e-6 absolute) tolerance is evaluated on the features whose fitted means are all ≥ 0.5 (DESeq2 floors mu at 0.5 in the standard errors) and whose pydeseq2 LFC fit converged; the others are counted apart, and every feature enters the comparison of the significant sets. Output lines, verbatim:

```
cohort_30 rna factor_additive: ~sex + stage + age_at_death + pmi + rin + years_education + log10_n_rna; 20387 features x 80 donors; DESeq2 1.52.0; 14.2 s; dispersion trend pydeseq2 parametric, DESeq2 parametric; rows not converged in beta: 1; size factors max rel diff 1.65e-15; final dispersions rel diff median 0.000188, 99% 0.00277, max 2.29
  trend coefficients a0 0.0955378 / 0.0955129, a1 1.61526 / 1.61703 (pydeseq2 / DESeq2)
  stage_1_vs_0 (stage_1_vs_0)
    full:    log2FC max abs 0.651 (234 beyond tolerance); stat max abs 0.733 (328 beyond tolerance); 2002 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 1.47, stat 1.54); padj <= 0.05: 0 pydeseq2, 0 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
    plug-in: log2FC max abs 0.0108 (5 beyond tolerance); stat max abs 0.0185 (5 beyond tolerance); 2002 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 2.22, stat 2.12); padj <= 0.05: 0 pydeseq2, 0 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
  stage_2_vs_0 (stage_2_vs_0)
    full:    log2FC max abs 1.02 (214 beyond tolerance); stat max abs 0.868 (301 beyond tolerance); 2002 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 0.885, stat 0.959); padj <= 0.05: 1 pydeseq2, 1 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
    plug-in: log2FC max abs 0.00635 (7 beyond tolerance); stat max abs 0.0128 (9 beyond tolerance); 2002 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 2.04, stat 2.27); padj <= 0.05: 1 pydeseq2, 1 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
tolerance 0.001 relative + 1e-06 absolute on features with fitted mu >= 0.5 and a converged pydeseq2 LFC fit (18 of 20387 with a non-converged one); full FAIL, plug-in FAIL
exit 4
cohort_30 rna factor_interaction: ~sex + stage + sex:stage + age_at_death + pmi + rin + years_education + log10_n_rna; 20387 features x 80 donors; DESeq2 1.52.0; 16.0 s; dispersion trend pydeseq2 parametric, DESeq2 parametric; rows not converged in beta: 1; size factors max rel diff 1.65e-15; final dispersions rel diff median 8.14e-05, 99% 0.00332, max 2.48
  trend coefficients a0 0.0954903 / 0.0954881, a1 1.61784 / 1.61126 (pydeseq2 / DESeq2)
  male_x_stage_1 (sexMale.stage1)
    full:    log2FC max abs 0.0173 (99 beyond tolerance); stat max abs 0.0343 (165 beyond tolerance); 2069 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 2.44, stat 2.46); padj <= 0.05: 0 pydeseq2, 0 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
    plug-in: log2FC max abs 0.0568 (2 beyond tolerance); stat max abs 0.0474 (2 beyond tolerance); 2069 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 2.8, stat 1.51); padj <= 0.05: 0 pydeseq2, 0 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
  male_x_stage_2 (sexMale.stage2)
    full:    log2FC max abs 0.0983 (88 beyond tolerance); stat max abs 0.277 (135 beyond tolerance); 2069 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 1.44, stat 1.76); padj <= 0.05: 0 pydeseq2, 0 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
    plug-in: log2FC max abs 0.105 (3 beyond tolerance); stat max abs 0.0989 (4 beyond tolerance); 2069 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 1.6, stat 1.24); padj <= 0.05: 0 pydeseq2, 0 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
tolerance 0.001 relative + 1e-06 absolute on features with fitted mu >= 0.5 and a converged pydeseq2 LFC fit (24 of 20387 with a non-converged one); full FAIL, plug-in FAIL
exit 4
cohort_30 rna trend_additive: ~sex + stage_num + age_at_death + pmi + rin + years_education + log10_n_rna; 20387 features x 80 donors; DESeq2 1.52.0; 13.2 s; dispersion trend pydeseq2 parametric, DESeq2 parametric; rows not converged in beta: 0; size factors max rel diff 1.65e-15; final dispersions rel diff median 0.000261, 99% 0.00284, max 2.28
  trend coefficients a0 0.0951036 / 0.0952415, a1 1.61154 / 1.60422 (pydeseq2 / DESeq2)
  stage_num (stage_num)
    full:    log2FC max abs 0.00333 (203 beyond tolerance); stat max abs 0.317 (357 beyond tolerance); 1958 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 0.719, stat 2.68); padj <= 0.05: 2 pydeseq2, 2 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
    plug-in: log2FC max abs 0.00333 (5 beyond tolerance); stat max abs 0.0142 (10 beyond tolerance); 1958 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 0.379, stat 1.04); padj <= 0.05: 2 pydeseq2, 2 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
tolerance 0.001 relative + 1e-06 absolute on features with fitted mu >= 0.5 and a converged pydeseq2 LFC fit (22 of 20387 with a non-converged one); full FAIL, plug-in FAIL
exit 4
cohort_30 rna trend_interaction: ~sex + stage_num + sex:stage_num + age_at_death + pmi + rin + years_education + log10_n_rna; 20387 features x 80 donors; DESeq2 1.52.0; 14.3 s; dispersion trend pydeseq2 parametric, DESeq2 parametric; rows not converged in beta: 3; size factors max rel diff 1.65e-15; final dispersions rel diff median 0.000191, 99% 0.00277, max 2.28
  trend coefficients a0 0.0948007 / 0.0949024, a1 1.6014 / 1.60229 (pydeseq2 / DESeq2)
  male_x_stage_num (sexMale.stage_num)
    full:    log2FC max abs 0.0308 (222 beyond tolerance); stat max abs 0.357 (318 beyond tolerance); 1982 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 0.805, stat 1.49); padj <= 0.05: 0 pydeseq2, 0 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
    plug-in: log2FC max abs 0.000798 (1 beyond tolerance); stat max abs 0.00203 (1 beyond tolerance); 1982 features with fitted mu < 0.5 or a non-converged pydeseq2 LFC fit apart (max abs log2FC 1.16, stat 1.42); padj <= 0.05: 0 pydeseq2, 0 DESeq2, 0 only pydeseq2, 0 only DESeq2; padj NA: 0 pydeseq2, 0 DESeq2 (same features) -> FAIL
tolerance 0.001 relative + 1e-06 absolute on features with fitted mu >= 0.5 and a converged pydeseq2 LFC fit (28 of 20387 with a non-converged one); full FAIL, plug-in FAIL
exit 4
```

  The criterion of task 3.4 (log2 fold changes and Wald statistics within 1e-3 on every feature) is not met; the difference is recorded here as the task requires. The significant sets at padj ≤ 0.05 are identical in every coefficient under both comparisons, and so are the features without padj. The plug-in comparison fails its tolerance too (exit 4 in every model), on 1–7 log2 fold changes and 1–10 Wald statistics of the comparable features per coefficient (up to 0.105 in log2FC and 0.0989 in stat): given the same dispersions R's GLM, Wald test and independent filter reproduce pydeseq2's on all but these few features, and that small remainder is not explained. The full comparison's further differences come from the dispersion estimation (pydeseq2 0.5.3 lacks DESeq2's reset of non-improving gene-wise fits, found on the synthetic data). The features set apart differ by construction: every feature whose pydeseq2 LFC fit did not converge differs from R's; in cohort_30 the smallest padj among them is 0.712 (RNA) and 0.761 (ATAC). Only the cohort_30 RNA fits were cross-checked; the ATAC fits, whose stage peaks sit closer to α, and cohort_50 were not.

## 5. Reading

- **Concordant in every cell, both negative.** No exported feature has padj ≤ 0.05 on an interaction coefficient of `factor_interaction` in either cohort or omic, and the primary analysis rejected `delta` in neither cohort, so all four cells are `concordant` (cohort_30 rna: 0 exported features, `delta` rejected no; cohort_30 atac: 0 exported features, `delta` rejected no; cohort_50 rna: 0 exported features, `delta` rejected no; cohort_50 atac: 0 exported features, `delta` rejected no). The cross-check finds nothing that contradicts the primary analysis; the `delta` non-rejection stays uninformative, and so does the absence of DE interactions.
- **No sex × stage interaction detected at this n, in any model.** On the abundance-passing set too, no feature reaches padj ≤ 0.05 on any interaction coefficient of `factor_interaction` or `trend_interaction`; the smallest interaction padj is 0.095 (cohort_30 rna), 0.421 (cohort_30 atac), 0.0681 (cohort_50 rna), 0.392 (cohort_50 atac). With cells of 9–28 donors this is "none detected at this n". On the synthetic data (60 / 54 donors, 9–10 per cell, dispersions 0.05–0.10) a female-only effect of 0.8 log2 units at stage 1 was found for 0–1 of 10 planted features and one of 1.6 at stage 2 for 9–10; the real fits' median dispersions are higher (0.10 RNA, 0.21 ATAC, cohort_30 `factor_interaction`), so the real power is lower, and these recalls do not bound the size of an undetected interaction.
- **Little per-feature stage signal.** In RNA cohort_30: 1 gene at stage 2 vs 0 (`CREB5`) and 2 on the stage trend (`CREB5`, `ADAM33`); cohort_50: 2 genes at stage 2 vs 0 (`CREB5`, `ADAM33`) and 2 on the stage trend (`CREB5`, `ADAM33`), none at stage 1 vs 0, all of them exported. In ATAC, cohort_30: 22 peaks at stage 2 vs 0 (0 exported) and 23 on the trend (0 exported); cohort_50: 36 peaks at stage 2 vs 0 (4 exported) and 30 on the trend (4 exported), none at stage 1 vs 0. Of the proposal's two readings this description favours little per-feature stage signal: the stage effect that reaches FDR is confined to a handful of features at the High stage. It does not settle the multivariate reading — the latent-adjusted arm's stage-tracking component can summarise many small shifts that no single feature carries at FDR 0.05 — and it cannot show that a stage signal is the same in both sexes, because no interaction test here has the power to show an absence.
- **MOTCO's ATAC block holds few of the stage peaks.** The variance-ranked 5,000 exported peaks contain 0 of the 22 stage-2 peaks (cohort_30), 4 of the 36 stage-2 peaks (cohort_50); both RNA genes are exported. This is a description of the measurement space only. The export's selection is label-free by design, and a feature set chosen from these results would read sex and stage and be circular, so it is no reason to change the export.

## 6. Limitations

- Power: n = 80 / 74 with sex × stage cells of 9–28 donors; the interaction tests have little power (section 5; the synthetic recalls overstate it), and a zero count is "none detected at this n".
- The FDR is approximate: DESeq2's Wald test plugs in estimated dispersions. In simulations on the synthetic 60-donor design under DESeq2's own dispersion model (`scratchpad/syn/round2_cal/calibration_pooled.csv` of the change), the realized FDR at padj ≤ 0.05 in the additive tables was 0.06 on average with 200 of 4,000 features non-null and 0.09 with 60 of 5,830, against 0.03–0.05 with the true dispersions. A count of a few features at padj ≤ 0.05 can hold a false discovery.
- Independent filtering chooses a baseMean cutoff per coefficient, so the BH families differ between coefficients: in ATAC `stage_2_vs_0` it removed 29993 (cohort_30), 38719 (cohort_50) and `stage_num` 32492 (cohort_30), 0 (cohort_50) of the tested peaks from the padj family.
- No outlier protection: with continuous covariates no two donors share a design row, so Cook's refit and filtering are inactive by construction (as in R's DESeq2). A feature at padj ≤ 0.05 can be driven by few donors; no outlier diagnostic was predeclared, and none was run.
- Covariates enter on their recorded scale; pydeseq2's IRLS did not converge for 2–28 features per fit (fitted by its fallback optimizer, which also raised the uncaptured worker RuntimeWarnings of section 1); in cohort_30 the smallest padj among them is 0.712 (RNA) and 0.761 (ATAC); cohort_50 was not checked. The union counts carry no joint error-rate claim (pydeseq2 has no likelihood-ratio test).
- DESeq2's own normalization (median of ratios on the abundance-passing set) is not the export's TMM; the feature set is the export's, the normalization is not.
- pydeseq2 and R's DESeq2 differ in the dispersion estimation, and slightly in the fits themselves (section 4). In the cohort_30 RNA fits checked the significant sets agree; the ATAC fits, whose stage peaks sit at padj close to α, and all of cohort_50 were not cross-checked with R.
- The counts are bound to the export by reproduction, not by a hash (the export records none): a change that left every status, library size, TMM factor, tested-feature variance and exported value unchanged would pass.
- `cohort_50` is a donor subset of `cohort_30`; agreement between them is not an independent replication.
- The `delta` side of the concordance cells is MOTCO's joint two-block latent-space test with the six analysis covariates; its operating characteristics were measured in MOTCO's own studies, on none of these blocks, and every limitation of the primary report (`reports/motco-astrocyte-2026-10-19.md` and its sensitivity reports) applies to it.

## 7. Next steps

- The proposed `mtg-microglia-positive-control` change runs this step (through its directory flags) and MOTCO on microglia, where a stage signal is expected, as a positive control of the pipeline's ability to find one. Nothing here changes the primary MOTCO analysis, which stays primary.

## 8. Reproduction

```
scripts/sync_pseudobulk.sh                                    # results/pseudobulk/ from the cluster
rsync -a <cluster>:/home1/tgonza/sead_sexdiff/results/{motco_input,motco_analysis} results/
uv run sead-sexdiff stage-de --cohort all --omic both --n-cpus 10
# R cross-check: debug files from the Python API, then the script per model
uv run python -c 'from sead_sexdiff.data import run_stage_de; run_stage_de("cohort_30", out_dir="/tmp/sd", debug_dir="/tmp/sd_debug", n_cpus=8, force=True)'
Rscript scripts/check_stage_de_deseq2.R /tmp/sd_debug /tmp/sd cohort_30 rna factor_interaction
```
