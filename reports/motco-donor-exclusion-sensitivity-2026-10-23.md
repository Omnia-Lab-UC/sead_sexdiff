# MOTCO donor-exclusion sensitivity: SEA-AD's severely affected donors, microglia and astrocytes — 2026-10-23

> The predeclared donor-exclusion sensitivity analysis of OpenSpec change `mtg-severely-affected-donors` (specs `donor-flag-table`, `donor-exclusion-sensitivity`, `pseudobulk-stage-de`): SEA-AD's donor-level `Severely Affected Donor` flag read from the donor objects (`donor-flags`), the flag check written before any run, the whole export and the primary analysis rerun on the cohort donors without the flagged ones (`exclusion-sensitivity`, variant `excluding_severely_affected` of `DONOR_EXCLUSION_VARIANTS`), the predeclared axis check, the first-axis correlations of the variant (`axis-diagnostics`) and both stage-DE methods on the exclusion cohort, for the control cell type (microglia) and the reference cell type (astrocytes). Gabitto et al. (2024) flagged 11 of the 42 high-ADNC donors as severely affected from a PCA of library QC metrics and excluded them from every gene-expression test; our pipeline kept them everywhere, and the primary analyses' first latent component of both cell types follows a donor-level quality axis on which stage and the mitochondrial fraction are entangled (`reports/motco-microglia-positive-control-2026-10-22.md`, `reports/motco-astrocyte-technical-axis-2026-10-19.md`).
>
> **Result.** The flagged donors carry the quality tail in both cell types: in `cohort_30` their rank-biserial AUC against the unflagged donors is 0.986 (microglia) and 0.961 (astrocytes) on the primary `pls_0` and 0.961 / 0.966 on `mitochondrial_pseudobulk`. **In microglia the predeclared axis check holds in both cohorts**: along the variant's first component |ρ| with stage is 0.8201 (`cohort_30`) and 0.8248 (`cohort_50`), against 0.6891 and 0.6746 for the largest technical variable, `fraction_mitochondrial_umis_sd`. In astrocytes (descriptive) it does not hold: the variant's first component follows `ribosomal_mean` (0.7286 / 0.6716) more than stage (0.5784 / 0.5955). No statistic is rejected in any of the eight rows; the microglia `cohort_50` variant selects rank 4 (primary 1). The exported RNA stage counts of microglia fall from 61 to 20 (`cohort_30`) and from 35 to 3 (`cohort_50`) under DESeq2's `factor_additive` union, while limma's joint F gives 2 and 2 (0 and 0); astrocyte counts stay at 1 to 2. Every concordance cell is concordant. Under the predeclared reading the quality gradient of the microglia primary analysis sits in the donors SEA-AD excluded, and excluding them is a donor-selection route to the same end as the latent adjustment; the primary analyses, the primary cohorts and the primary-basis positive control are unchanged.

**Reading rules (binding).** The variant is a sensitivity analysis, never a revised primary: the microglia and astrocyte primary analyses stay the primary results, `cohort_30` is primary and `cohort_50` sensitivity, and the primary-basis positive control (`reports/motco-microglia-positive-control-2026-10-22.md`) keeps its outcome. A rejection status that differs between the variant and the primary is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections stay uninformative and a `delta` rejection is a difference, not specifically of size, in every row. The flag check is descriptive (no p-value, never a gate); the axis check is a comparison, not a test, whose microglia outcome is the change's expectation and whose astrocyte outcome is descriptive. The stage-DE counts describe how much of the stage signal survives the exclusion, and the concordance cells cross-check the variant analysis as they cross-check the primary; neither is a second test of the sex × stage hypothesis. The cross-cell-type comparisons are compared, not tested. The flag check and the axis check are read together under the reading written in the proposal before any run (`DONOR_EXCLUSION_READING`, §6).

## 1. Provenance

- Software: MOTCO 0.9.0 (git dependency pinned to `v0.9.0`, commit `0150adc27b48528f1999a31a146a3df17ae3047f` in `uv.lock`), scikit-learn 1.8.0, numpy 2.3.5, scipy 1.15.3, pandas 2.3.3, pydeseq2 0.5.3, Python 3.12.15, R 4.6.1 with limma 3.68.4 and edgeR 4.10.1 (the Nix shell). sead-sexdiff revision `404be4bb2e1f643ae9df36e1b4d932db84d22497` for every workstation run (recorded `dirty: True`: at run time the uncommitted files were `README.md`, `AGENTS.md`, the change's `tasks.md` and the IDE file `.idea/sead-sexdiff.iml`, none under `src/`); `donor-flags` at `82cd6b1f3573e63275a16317a2b775b33cfc84f4` (clean).
- Flags: `donor-flags --omic both` and `donor-flags --cell-type microglia --omic both` on the cluster login node (`cluster.ing.uc.cl`) on 2026-10-23 over the 84 RNA (2026-06-22) and 86 ATAC (2024-12-06) donor objects, `obs` only (under 0.1 s per object). Both runs list 84 non-Reference donors and 11 flagged (SEA-AD's count), leave out the Reference donors `H19.30.002` and `H19.30.004`, find no flag and no label disagreement between omics, and record one finding: the RNA objects lack `Continuous Pseudo-progression Score` (the score comes from ATAC). The two `donor_flags.csv` are identical (sha256 `58d38a9269378da844a71dcbc10fa71fddf0b38668e6e2c24789bea8d4653136`); manifests `cc64c5ee8bc41f219b6074bb4b4a5e536321f1758f4fd9079f9361ec1b0e7076` (astrocyte) and `ddd14570b10e3a7d0308967f77ecb317d6aa1a7f4c5a7f2f8475846c922dedd9` (microglia). Each pseudobulk directory gained exactly the two files; the 726 files there before (per-donor caches and logs included) hashed unchanged on the cluster, and the synced directories hash as the cluster's.
- Machine: `nixos-home` (NixOS workstation), 2026-10-23, `n_jobs = 1`, both cell types' `exclusion-sensitivity` runs in parallel. Every manifest read was checked for the cell type; the astrocyte pseudobulk, technical, export and analysis manifests predate the cell-type registry and read as `astrocyte` with a finding; the microglia ones record `microglia`.
- Microglia inputs: pseudobulk manifest `49b937223a485a5233fcc7953fc51f9cd343ba65ccbc5c2dc2c3800ce8d8cb1c`, technical table `ac54eef6a253135db7260eec0365b837c103495f1d3e1041c9b9e06ff288b5cd`; `cohort_30` primary export manifest `fb6afd9c24d8e6eee88bad4d8cb5a3ead2ddbf282b5ecff4b196bb893b8a9936` and `result.json` `c72967f920b2c766b31cd1a6957761f5e1d2c94c64f3cb25bda0c1e60ec2ceb4`; `cohort_50` `8f0fd84ab291f93c60ea40f8f864ce8aad830005331918a8f2c65d97f2389bad` and `bf56d2c959dbffbf31fb8bd44698b191b52ceb2b3c2eeb7065df0e24eb56d77e`. Variants: `cohort_30` export manifest `10ebd0e43fe0cc09c59253487dad23c573dc345b6cce688cb90304c6f9b5b356`, `result.json` `ca205407495c25ed58cc98a2a8ce9dd8b7ee57e74b0f28e5a7edc67d4148ffe3`; `cohort_50` `089ba00e19e498d6422f7edcbbbe8cb77783d32bff40d069ca45a4ca6416d0d9`, `67246cad7ae14c24eaf9e02748f5b8b51dcb798a0c3babf1d225bb9b96e7e498`.
- Astrocyte inputs: pseudobulk manifest `2e93cc0b526856fc043f908a424a11266187146636bc75afed7edd187989a25c`, technical table `eac1b9a15977dfab683738a4bd0b0ffc207b896a102032409b60fccbfb96247b`; `cohort_30` primary export manifest `18f2289e0d4c84ed18f4a437ae569052544e9d150ed2726385273694b4a3b0b6` and `result.json` `1d6154b02d7704ef9313ce752035a124cf6cb3fbd5350f0428439a542d5190cd`; `cohort_50` `dd14a6e7d5583fb527e8a99c7db43d12cf237f5451522e8875c5a0a3401c5cb3` and `4c6e7e58e39d8ce82c9b93907b70b2110fc2c415e32da40fef8668b74faff9fb`. Variants: `cohort_30` export manifest `dab3c61d6090cfbfc180870d05c02b65dfa4d5e678086b248717a9b9d3b95fd0`, `result.json` `3137eec201c48d3458f63bb4489dba3a580a02236cc0880e1db7a86adb5ae35f`; `cohort_50` `dd36a8a6aaa8ca6f95feb2a6f616ebbda740dafda002c9f93ed1431603671760`, `29cadd9c1bbfee592976301edf5e3e380607e51321d526b37071179a83b0b8fb`.
- Checks before any write (every cohort, both cell types): the primary export against the pseudobulk, covariates and annotation hashes, the export constants, `k`, the cells and the donor count; the primary analysis against the export hashes, MOTCO 0.9.0, seed 20261019, 999 permutations, CV selection, orientation, CV settings and the six covariates; the technical table against the pseudobulk manifest and `TECHNICAL_GENE_SETS`, complete for the cohort donors; a flag for every cohort donor. Each `flag_check.csv` was put in place before its variant export ran (recorded `written_at`, file times). The variant exports' `metadata.csv` and `covariates.csv` equal the primary's minus the excluded rows. The 306 files under both cell types' `pseudobulk/`, `motco_input/`, `motco_analysis/`, `stage_de/` and `stage_de_limma/` hashed unchanged after every run of this change.
- Settings shared by the eight rows: seed 20261019, 999 RRPP permutations, α = 0.05 (p ≤ α), `n_jobs = 1`; blocks `rna` and `atac` in identity units, stage order 0, 1, 2, net-displacement orientation, CV rank selection at MOTCO's defaults, the six analysis covariates (`age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`), missing-covariate policy `error`. The variant exports keep every primary export setting (chrX/Y/M removed, CPM ≥ 1 in ≥ k donors, TMM, prior count 2, top 5,000 by variance per block); the cohort thresholds are unchanged and `k` is recomputed from the remaining cells by the export's rule. No flagged donor sits in a smallest cell, so `k` is unchanged in every cohort (astrocyte 9 / 9, microglia 9 / 8).
- Diagnostics: `axis-diagnostics [--cell-type microglia] --input-dir <variant>/motco_input --analysis-dir <variant>/motco_analysis --rank-dir <variant>/motco_rank_sensitivity --out-dir <variant>/diagnostics` per root (`<variant>` = `<root>/motco_exclusion_sensitivity/excluding_severely_affected`); the input-hash and cell-type checks passed; `--rank-dir` names a directory that does not exist (no sweep on a variant; the astrocyte default sweep was made from the primary export), recorded as a finding. The diagnostic labels the analysis it reads `primary`; the tables below name it the variant. Its partial correlations condition on the six analysis covariates.
- Stage DE on the variants: `stage-de --omic both --n-cpus 10` and then `stage-de-limma --omic rna` per root and cohort, with the variant's `motco_input/` and `motco_analysis/`. **Deviation from the proposal** (decided with the user on 2026-10-23): DESeq2 ran on both omics, not RNA only, because `stage-de-limma` requires a DESeq2 run complete for both omics (spec `pseudobulk-stage-de-limma`), so the limma step and its spec stay unchanged; the ATAC DESeq2 tables are recorded under the variant roots and are not read here. Both step manifests record the export's `donor_exclusion` block; the tested RNA sets are the variant exports' abundance-passing features (astrocyte 20,191 / 20,037, microglia 20,421 / 20,154). `--n-cpus 10` leaves pydeseq2's worker warnings uncaptured (a finding of the step, as for the primary runs).
- Reproducibility and copies: `exclusion-sensitivity --cell-type microglia --cohort cohort_30 --force` into a scratch directory gave byte-identical variant export CSVs (six), `flag_check.csv`, `axis_check.csv`, `exclusion_sensitivity.csv` and analysis CSVs (13 files) and equal statistics, p-values, null summaries and rank. Before the real runs, the refactored package reproduced a reduced fixture of the astrocyte pseudobulk through `export`, `analyze`, `stage-de` and `stage-de-limma` byte for byte (82 CSVs, both `summary.md`; manifests apart from times, runtimes and the hashes of files carrying them; no `donor_exclusion` key on the primary path). `results/motco_exclusion_sensitivity/` and `results/microglia/motco_exclusion_sensitivity/` were copied to the cluster (`/home1/tgonza/sead_sexdiff/results/`) and all 256 files hashed identically there.
- Numbers below are generated from the output files by a scratch script and pasted; the comparison sections are the generated `exclusion_sensitivity.md` files verbatim (their runtimes are the only values not in a CSV). A check script found every statistic, p-value, null summary, rank and CV AUROC of the eight `result.json` records, every flag-check and axis-check value and every stage-DE count in the text.

## 2. The flag table

All 11 flagged donors are stage 2 (high ADNC), 8 female and 3 male; 10 of them have dementia, all are Braak V or VI, and their continuous pseudo-progression scores (recorded, read by nothing) range from 0.729331972 to 0.924064541. Ten are in both astrocyte cohorts, but only 7 enter microglia `cohort_30` and 5 microglia `cohort_50`: four flagged donors (H21.33.020, H21.33.033, H21.33.039, H21.33.045) have fewer than 30 microglia ATAC nuclei (11 to 24; SEA-AD reports many fewer ATAC peaks for the flagged donors), and one astrocyte one (H20.33.033, 20 astrocyte ATAC nuclei) misses the astrocyte cohorts, so the flag removes 9 % of the microglia `cohort_30` donors and 12.5 % of the astrocyte one. Every flagged donor sits in a stage-2 cell, so the exclusion shrinks only the stage-2 cells, the female one most (astrocyte `cohort_30` 28 → 20, microglia 24 → 19), and changes no smallest cell.

Flagged donors per sex × stage cell (flagged / donors; generated from `donor_flags.csv` and each cell type's `donor_table.csv`):

| cell type | donors | flagged | Female/stage 0 | Female/stage 1 | Female/stage 2 | Male/stage 0 | Male/stage 1 | Male/stage 2 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| microglia | donor table (84) | 11 | 0 / 12 | 0 / 10 | 8 / 29 | 0 / 9 | 0 / 11 | 3 / 13 |
| microglia | cohort_30 (76) | 7 | 0 / 11 | 0 / 10 | 5 / 24 | 0 / 9 | 0 / 11 | 2 / 11 |
| microglia | cohort_50 (63) | 5 | 0 / 8 | 0 / 10 | 4 / 18 | 0 / 8 | 0 / 10 | 1 / 9 |
| astrocyte | donor table (84) | 11 | 0 / 12 | 0 / 10 | 8 / 29 | 0 / 9 | 0 / 11 | 3 / 13 |
| astrocyte | cohort_30 (80) | 10 | 0 / 11 | 0 / 10 | 8 / 28 | 0 / 9 | 0 / 10 | 2 / 12 |
| astrocyte | cohort_50 (74) | 10 | 0 / 10 | 0 / 9 | 8 / 25 | 0 / 9 | 0 / 9 | 2 / 12 |

The flagged donors (labels as recorded, read by nothing):

| donor | sex | stage | cognitive status | Braak | CERAD | Thal | LATE | CPS | in astrocyte cohort_30 / 50 | in microglia cohort_30 / 50 |
|---|---|---:|---|---|---|---|---|---:|---|---|
| `H20.33.028` | Female | 2 | Dementia | Braak V | Moderate | Thal 4 | LATE Stage 2 | 0.814933639 | yes / yes | yes / yes |
| `H20.33.033` | Male | 2 | Dementia | Braak VI | Frequent | Thal 5 | Not Identified | 0.874937976 | no / no | yes / no |
| `H20.33.037` | Female | 2 | Dementia | Braak V | Frequent | Thal 5 | LATE Stage 2 | 0.891215423 | yes / yes | yes / yes |
| `H20.33.038` | Female | 2 | Dementia | Braak V | Frequent | Thal 4 | LATE Stage 2 | 0.748876479 | yes / yes | yes / yes |
| `H20.33.046` | Male | 2 | Dementia | Braak VI | Frequent | Thal 5 | LATE Stage 2 | 0.868382665 | yes / yes | yes / yes |
| `H21.33.009` | Female | 2 | Dementia | Braak VI | Frequent | Thal 5 | Not Identified | 0.911019291 | yes / yes | yes / yes |
| `H21.33.020` | Male | 2 | Dementia | Braak V | Frequent | Thal 4 | LATE Stage 2 | 0.911076764 | yes / yes | no / no |
| `H21.33.033` | Female | 2 | No dementia | Braak V | Moderate | Thal 5 | LATE Stage 1 | 0.899911521 | yes / yes | no / no |
| `H21.33.034` | Female | 2 | Dementia | Braak VI | Frequent | Thal 5 | LATE Stage 2 | 0.869521183 | yes / yes | yes / no |
| `H21.33.039` | Female | 2 | Dementia | Braak V | Moderate | Thal 4 | LATE Stage 3 | 0.729331972 | yes / yes | no / no |
| `H21.33.045` | Female | 2 | Dementia | Braak VI | Frequent | Thal 4 | LATE Stage 3 | 0.924064541 | yes / yes | no / no |

## 3. Microglia, the control cell type

Microglia is the control cell type: its axis-check outcome is the change's expectation. The primary-basis positive control failed on E3 here (`cohort_30`: |ρ| of the primary `pls_0` with stage 0.688 against 0.844 with the mitochondrial fraction), and the restated control on MOTCO's `latent_adjusted` basis held (stage 0.5279 against `tmm_rna` 0.4153; `reports/motco-microglia-latent-adjustment-2026-10-23.md`); both reports' numbers are cited, not regenerated.

### 3.1 cohort_30 (≥ 30 microglia nuclei in both omics, n = 76) — primary cohort

#### Overview (generated from `exclusion_sensitivity.csv`)

| analysis | donors | excluded | k | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| primary | 76 | 0 | 9 | 1 | 0.6921 | 8.0642 (0.391) | 0 (1) | 0.154684 (0.637) | none |
| excluding_severely_affected | 69 | 7 | 9 | 1 | 0.6388 | 5.34849 (0.379) | 0 (1) | 0.187608 (0.366) | none |

#### Axis check (from `axis_check.csv`)

Along the variant's `pls_0` (n = 69): |ρ| with stage 0.8201 (ρ = 0.8201) against the largest of 39 technical variables, `fraction_mitochondrial_umis_sd` (donor_technical.csv) |ρ| 0.6891 (ρ = 0.6891): **holds** (expectation, role control).

#### Stage DE on the exclusion cohort, RNA (from `stage_de_summary.csv` and `stage_de_limma_summary.csv`; counts at FDR 0.05)

| method | model | coefficient / test | meaning | set | primary (n tested) | variant (n tested) |
|---|---|---|---|---|---:|---:|
| DESeq2 | `factor_additive` | `union` | E1's DESeq2 quantity: stage, either coefficient | exported | 61 (5000) | 20 (5000) |
| DESeq2 | `factor_additive` | `union` | E1's DESeq2 quantity: stage, either coefficient | abundance | 156 (20660) | 45 (20421) |
| DESeq2 | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | exported | 61 (5000) | 20 (5000) |
| DESeq2 | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | abundance | 156 (20660) | 45 (20421) |
| DESeq2 | `trend_additive` | `stage_num` | stage trend | exported | 57 (5000) | 21 (5000) |
| DESeq2 | `trend_additive` | `stage_num` | stage trend | abundance | 145 (20660) | 46 (20421) |
| DESeq2 | `factor_interaction` | `union` | sex × stage, either coefficient (concordance quantity) | exported | 0 (5000) | 0 (5000) |
| DESeq2 | `factor_interaction` | `union` | sex × stage, either coefficient (concordance quantity) | abundance | 0 (20660) | 0 (20421) |
| limma | `factor_additive` | `stage_joint` | E1's limma quantity: joint F of stage | exported | 2 (5000) | 2 (5000) |
| limma | `factor_additive` | `stage_joint` | E1's limma quantity: joint F of stage | abundance | 2 (20660) | 2 (20421) |
| limma | `factor_additive` | `union` | stage, either moderated t | exported | 28 (5000) | 2 (5000) |
| limma | `factor_additive` | `union` | stage, either moderated t | abundance | 64 (20660) | 2 (20421) |
| limma | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | exported | 28 (5000) | 2 (5000) |
| limma | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | abundance | 64 (20660) | 2 (20421) |
| limma | `trend_additive` | `stage_num` | stage trend | exported | 31 (5000) | 2 (5000) |
| limma | `trend_additive` | `stage_num` | stage trend | abundance | 67 (20660) | 2 (20421) |
| limma | `factor_interaction` | `interaction_joint` | joint F of sex × stage (concordance quantity) | exported | 0 (5000) | 0 (5000) |
| limma | `factor_interaction` | `interaction_joint` | joint F of sex × stage (concordance quantity) | abundance | 0 (20660) | 0 (20421) |
| limma | `factor_interaction` | `union` | sex × stage, either moderated t | exported | 0 (5000) | 0 (5000) |
| limma | `factor_interaction` | `union` | sex × stage, either moderated t | abundance | 0 (20660) | 0 (20421) |

| analysis | DESeq2 interaction count (exported) | DESeq2 cell | limma interaction count (exported) | limma cell | delta rejected | methods agree |
|---|---:|---|---:|---|:---:|:---:|
| primary | 0 | concordant | 0 | concordant | no | yes |
| excluding_severely_affected | 0 | concordant | 0 | concordant | no | yes |

#### First latent axis of the variant (generated from `diagnostics/cohort_30/axis_correlations.csv`, `pls_0`, the 12 largest |ρ| and stage)

| variable | source | Spearman | partial Spearman |
|---|---|---:|---:|
| `stage` | metadata.csv | 0.820 | 0.800 |
| `fraction_mitochondrial_umis_sd` | donor_technical.csv | 0.689 | 0.708 |
| `mitochondrial_sd` | donor_technical.csv | 0.689 | 0.708 |
| `mitochondrial_pseudobulk` | donor_technical.csv | 0.688 | 0.687 |
| `fraction_mitochondrial_umis_mean` | donor_technical.csv | 0.671 | 0.678 |
| `mitochondrial_mean` | donor_technical.csv | 0.671 | 0.678 |
| `gex_sequencing_saturation_mean` | donor_technical.csv | 0.528 | 0.505 |
| `atac_tss_enrichment_score_mean` | donor_technical.csv | -0.475 | -0.463 |
| `mitochondrial_supertype_share` | donor_technical.csv | 0.458 | 0.500 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | donor_technical.csv | -0.445 | -0.419 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | -0.399 | -0.413 |
| `genes_detected_mean` | donor_technical.csv | 0.390 | 0.331 |

#### Comparison (pasted from `results/microglia/motco_exclusion_sensitivity/excluding_severely_affected/comparison/cohort_30/exclusion_sensitivity.md`)

#### Donor-exclusion sensitivity: microglia, cohort_30, excluding_severely_affected

The `primary` row is the primary analysis (CV-selected rank 1, 76 donors), read from its `result.json` and never rerun; it stays the primary result. The `excluding_severely_affected` row reruns the whole export and the primary analysis on the cohort donors minus those SEA-AD flagged as severely affected (Gabitto et al. 2024: QC-outlying high-ADNC donors excluded from the gene-expression tests): 7 cohort donor(s) excluded, the cohort thresholds unchanged and k recomputed from the remaining cells by the export's rule (k = 9 → 9). The variant is a sensitivity analysis, never a revised primary: the excluded donors carry the most severe pathology, so the exclusion truncates the late epoch along with the quality tail, and the flag was made from QC metrics, so the variant conditions on quality by donor selection. A rejection status that differs between the rows is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in both rows. 999 RRPP permutations with seed 20261019 in both rows; α = 0.05 (rejected when p ≤ α). MOTCO's operating characteristics were measured on neither cohort.

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
| primary | 76 | 1 | 8.0642 | 0.391 | no | 7.23504 | 5.27064 | 17.4616 |
| excluding_severely_affected | 69 | 1 | 5.34849 | 0.379 | no | 4.90174 | 3.75421 | 12.1065 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 76 | 1 | 0 | 1 | no | 0 | 0 | 0 |
| excluding_severely_affected | 69 | 1 | 0 | 1 | no | 0 | 0 | 0 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 76 | 1 | 0.154684 | 0.637 | no | 0.25893 | 0.192607 | 0.638734 |
| excluding_severely_affected | 69 | 1 | 0.187608 | 0.366 | no | 0.161436 | 0.122576 | 0.397789 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 76 | 0 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=24; Male/stage 0=9; Male/stage 1=11; Male/stage 2=11 | cv | 1 | 0.6921 | 88.9 | `results/microglia/motco_analysis/cohort_30` |
| excluding_severely_affected | 69 | 7 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=19; Male/stage 0=9; Male/stage 1=11; Male/stage 2=9 | cv | 1 | 0.6388 | 108.8 | `motco_analysis/cohort_30` |

Excluded donors: `H20.33.028`, `H20.33.033`, `H20.33.037`, `H20.33.038`, `H20.33.046`, `H21.33.009`, `H21.33.034`.

##### Axis check

Predeclared rule (`DONOR_EXCLUSION_AXIS_RULE`): along the variant analysis's first latent component, |Spearman| with stage exceeds the largest |Spearman| with a technical variable of POSITIVE_CONTROL_TECHNICAL_SOURCES (technical table columns for the variant donors; the variant export's TMM factors and log10 library sizes); the control cell type's value is the expectation, the reference's descriptive. A comparison, not a test; no p-value.

| cell type | role | reading | component | donors | ρ stage | largest technical | source | ρ technical | technical compared | holds |
|---|---|---|---|---:|---:|---|---|---:|---:|:---:|
| microglia | control | expectation | `pls_0` | 69 | 0.8201 | `fraction_mitochondrial_umis_sd` | donor_technical.csv | 0.6891 | 39 | yes |

##### Reading

Predeclared (`DONOR_EXCLUSION_READING`): the flag check and the axis check are read together: if the flagged donors carry the quality tail (high AUC for the mitochondrial fraction and the primary pls_0) and the axis check holds in the control cell type, the quality gradient of the primary analyses sits in the donors SEA-AD excluded from expression tests, and excluding them is a donor-selection route to the same end as the latent adjustment; if the axis check fails although the flagged donors carry the tail, the gradient extends beyond them and donor selection does not resolve it; if the flagged donors do not carry the tail, the flag is not the gradient's carrier and the variant is an ordinary cohort sensitivity. In every case the primary analyses and the primary-basis positive control stay as they are; a rejection status that differs between the variant and the primary is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside angle; the non-rejection rules hold in the variant.

##### Findings

- axis check holds in cohort_30 (expectation, role control): |Spearman| of pls_0 with stage 0.8201 vs 0.6891 with fraction_mitochondrial_umis_sd (donor_technical.csv)
- excluding_severely_affected export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### 3.2 cohort_50 (≥ 50 microglia nuclei in both omics, n = 63) — sensitivity cohort

#### Overview (generated from `exclusion_sensitivity.csv`)

| analysis | donors | excluded | k | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| primary | 63 | 0 | 8 | 1 | 0.6414 | 12.4641 (0.227) | 0 (1) | 0.269568 (0.399) | none |
| excluding_severely_affected | 58 | 5 | 8 | 4 | 0.6280 | 2.55257 (0.758) | 19.1427 (0.45) | 0.0880231 (0.52) | none |

#### Axis check (from `axis_check.csv`)

Along the variant's `pls_0` (n = 58): |ρ| with stage 0.8248 (ρ = 0.8248) against the largest of 39 technical variables, `fraction_mitochondrial_umis_sd` (donor_technical.csv) |ρ| 0.6746 (ρ = 0.6746): **holds** (expectation, role control).

#### Stage DE on the exclusion cohort, RNA (from `stage_de_summary.csv` and `stage_de_limma_summary.csv`; counts at FDR 0.05)

| method | model | coefficient / test | meaning | set | primary (n tested) | variant (n tested) |
|---|---|---|---|---|---:|---:|
| DESeq2 | `factor_additive` | `union` | E1's DESeq2 quantity: stage, either coefficient | exported | 35 (5000) | 3 (5000) |
| DESeq2 | `factor_additive` | `union` | E1's DESeq2 quantity: stage, either coefficient | abundance | 65 (20272) | 11 (20154) |
| DESeq2 | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | exported | 34 (5000) | 2 (5000) |
| DESeq2 | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | abundance | 64 (20272) | 10 (20154) |
| DESeq2 | `trend_additive` | `stage_num` | stage trend | exported | 34 (5000) | 0 (5000) |
| DESeq2 | `trend_additive` | `stage_num` | stage trend | abundance | 59 (20272) | 0 (20154) |
| DESeq2 | `factor_interaction` | `union` | sex × stage, either coefficient (concordance quantity) | exported | 0 (5000) | 0 (5000) |
| DESeq2 | `factor_interaction` | `union` | sex × stage, either coefficient (concordance quantity) | abundance | 0 (20272) | 0 (20154) |
| limma | `factor_additive` | `stage_joint` | E1's limma quantity: joint F of stage | exported | 0 (5000) | 0 (5000) |
| limma | `factor_additive` | `stage_joint` | E1's limma quantity: joint F of stage | abundance | 0 (20272) | 0 (20154) |
| limma | `factor_additive` | `union` | stage, either moderated t | exported | 36 (5000) | 0 (5000) |
| limma | `factor_additive` | `union` | stage, either moderated t | abundance | 70 (20272) | 0 (20154) |
| limma | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | exported | 36 (5000) | 0 (5000) |
| limma | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | abundance | 70 (20272) | 0 (20154) |
| limma | `trend_additive` | `stage_num` | stage trend | exported | 35 (5000) | 0 (5000) |
| limma | `trend_additive` | `stage_num` | stage trend | abundance | 67 (20272) | 0 (20154) |
| limma | `factor_interaction` | `interaction_joint` | joint F of sex × stage (concordance quantity) | exported | 0 (5000) | 0 (5000) |
| limma | `factor_interaction` | `interaction_joint` | joint F of sex × stage (concordance quantity) | abundance | 0 (20272) | 0 (20154) |
| limma | `factor_interaction` | `union` | sex × stage, either moderated t | exported | 0 (5000) | 0 (5000) |
| limma | `factor_interaction` | `union` | sex × stage, either moderated t | abundance | 0 (20272) | 0 (20154) |

| analysis | DESeq2 interaction count (exported) | DESeq2 cell | limma interaction count (exported) | limma cell | delta rejected | methods agree |
|---|---:|---|---:|---|:---:|:---:|
| primary | 0 | concordant | 0 | concordant | no | yes |
| excluding_severely_affected | 0 | concordant | 0 | concordant | no | yes |

#### First latent axis of the variant (generated from `diagnostics/cohort_50/axis_correlations.csv`, `pls_0`, the 12 largest |ρ| and stage)

| variable | source | Spearman | partial Spearman |
|---|---|---:|---:|
| `stage` | metadata.csv | 0.825 | 0.785 |
| `fraction_mitochondrial_umis_sd` | donor_technical.csv | 0.675 | 0.743 |
| `mitochondrial_sd` | donor_technical.csv | 0.675 | 0.743 |
| `mitochondrial_pseudobulk` | donor_technical.csv | 0.666 | 0.724 |
| `fraction_mitochondrial_umis_mean` | donor_technical.csv | 0.659 | 0.734 |
| `mitochondrial_mean` | donor_technical.csv | 0.659 | 0.734 |
| `gex_sequencing_saturation_mean` | donor_technical.csv | 0.604 | 0.654 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | donor_technical.csv | -0.493 | -0.509 |
| `number_of_umis_mean` | donor_technical.csv | 0.434 | 0.445 |
| `number_of_umis_sd` | donor_technical.csv | 0.428 | 0.389 |
| `genes_detected_mean` | donor_technical.csv | 0.418 | 0.379 |
| `atac_tss_enrichment_score_mean` | donor_technical.csv | -0.417 | -0.488 |

#### Comparison (pasted from `results/microglia/motco_exclusion_sensitivity/excluding_severely_affected/comparison/cohort_50/exclusion_sensitivity.md`)

#### Donor-exclusion sensitivity: microglia, cohort_50, excluding_severely_affected

The `primary` row is the primary analysis (CV-selected rank 1, 63 donors), read from its `result.json` and never rerun; it stays the primary result. The `excluding_severely_affected` row reruns the whole export and the primary analysis on the cohort donors minus those SEA-AD flagged as severely affected (Gabitto et al. 2024: QC-outlying high-ADNC donors excluded from the gene-expression tests): 5 cohort donor(s) excluded, the cohort thresholds unchanged and k recomputed from the remaining cells by the export's rule (k = 8 → 8). The variant is a sensitivity analysis, never a revised primary: the excluded donors carry the most severe pathology, so the exclusion truncates the late epoch along with the quality tail, and the flag was made from QC metrics, so the variant conditions on quality by donor selection. A rejection status that differs between the rows is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in both rows. 999 RRPP permutations with seed 20261019 in both rows; α = 0.05 (rejected when p ≤ α). MOTCO's operating characteristics were measured on neither cohort.

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
| primary | 63 | 1 | 12.4641 | 0.227 | no | 8.45222 | 5.99525 | 19.5118 |
| excluding_severely_affected | 58 | 4 | 2.55257 | 0.758 | no | 6.67341 | 5.02936 | 16.7839 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 63 | 1 | 0 | 1 | no | 0 | 0 | 0 |
| excluding_severely_affected | 58 | 4 | 19.1427 | 0.45 | no | 18.7309 | 8.58051 | 34.5261 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 63 | 1 | 0.269568 | 0.399 | no | 0.255498 | 0.194219 | 0.624989 |
| excluding_severely_affected | 58 | 4 | 0.0880231 | 0.52 | no | 0.0962721 | 0.0495251 | 0.183762 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 63 | 0 | 8 | Female/stage 0=8; Female/stage 1=10; Female/stage 2=18; Male/stage 0=8; Male/stage 1=10; Male/stage 2=9 | cv | 1 | 0.6414 | 84.1 | `results/microglia/motco_analysis/cohort_50` |
| excluding_severely_affected | 58 | 5 | 8 | Female/stage 0=8; Female/stage 1=10; Female/stage 2=14; Male/stage 0=8; Male/stage 1=10; Male/stage 2=8 | cv | 4 | 0.6280 | 99.4 | `motco_analysis/cohort_50` |

Excluded donors: `H20.33.028`, `H20.33.037`, `H20.33.038`, `H20.33.046`, `H21.33.009`.

##### Axis check

Predeclared rule (`DONOR_EXCLUSION_AXIS_RULE`): along the variant analysis's first latent component, |Spearman| with stage exceeds the largest |Spearman| with a technical variable of POSITIVE_CONTROL_TECHNICAL_SOURCES (technical table columns for the variant donors; the variant export's TMM factors and log10 library sizes); the control cell type's value is the expectation, the reference's descriptive. A comparison, not a test; no p-value.

| cell type | role | reading | component | donors | ρ stage | largest technical | source | ρ technical | technical compared | holds |
|---|---|---|---|---:|---:|---|---|---:|---:|:---:|
| microglia | control | expectation | `pls_0` | 58 | 0.8248 | `fraction_mitochondrial_umis_sd` | donor_technical.csv | 0.6746 | 39 | yes |

##### Reading

Predeclared (`DONOR_EXCLUSION_READING`): the flag check and the axis check are read together: if the flagged donors carry the quality tail (high AUC for the mitochondrial fraction and the primary pls_0) and the axis check holds in the control cell type, the quality gradient of the primary analyses sits in the donors SEA-AD excluded from expression tests, and excluding them is a donor-selection route to the same end as the latent adjustment; if the axis check fails although the flagged donors carry the tail, the gradient extends beyond them and donor selection does not resolve it; if the flagged donors do not carry the tail, the flag is not the gradient's carrier and the variant is an ordinary cohort sensitivity. In every case the primary analyses and the primary-basis positive control stay as they are; a rejection status that differs between the variant and the primary is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside angle; the non-rejection rules hold in the variant.

##### Findings

- CV-selected rank differs: 1 in the primary analysis, 4 in the excluding_severely_affected analysis; reported beside angle, not read as better or worse
- axis check holds in cohort_50 (expectation, role control): |Spearman| of pls_0 with stage 0.8248 vs 0.6746 with fraction_mitochondrial_umis_sd (donor_technical.csv)
- excluding_severely_affected export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

## 4. Astrocytes, the reference cell type

Astrocytes are the reference cell type: their axis-check outcome is descriptive. The astrocyte primary analyses select rank 1 with a first axis that the technical-axis report identified as a donor-level library-quality axis (`reports/motco-astrocyte-technical-axis-2026-10-19.md`).

### 4.1 cohort_30 (≥ 30 astrocyte nuclei in both omics, n = 80) — primary cohort

#### Overview (generated from `exclusion_sensitivity.csv`)

| analysis | donors | excluded | k | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| primary | 80 | 0 | 9 | 1 | 0.6313 | 14.941 (0.3) | 0 (1) | 0.608832 (0.333) | none |
| excluding_severely_affected | 70 | 10 | 9 | 1 | 0.5559 | 6.35558 (0.534) | 0 (1) | 0.495893 (0.173) | none |

#### Axis check (from `axis_check.csv`)

Along the variant's `pls_0` (n = 70): |ρ| with stage 0.5784 (ρ = 0.5784) against the largest of 39 technical variables, `ribosomal_mean` (donor_technical.csv) |ρ| 0.7286 (ρ = -0.7286): **does not hold** (descriptive, role reference).

#### Stage DE on the exclusion cohort, RNA (from `stage_de_summary.csv` and `stage_de_limma_summary.csv`; counts at FDR 0.05)

| method | model | coefficient / test | meaning | set | primary (n tested) | variant (n tested) |
|---|---|---|---|---|---:|---:|
| DESeq2 | `factor_additive` | `union` | E1's DESeq2 quantity: stage, either coefficient | exported | 1 (5000) | 1 (5000) |
| DESeq2 | `factor_additive` | `union` | E1's DESeq2 quantity: stage, either coefficient | abundance | 1 (20387) | 1 (20191) |
| DESeq2 | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | exported | 1 (5000) | 1 (5000) |
| DESeq2 | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | abundance | 1 (20387) | 1 (20191) |
| DESeq2 | `trend_additive` | `stage_num` | stage trend | exported | 2 (5000) | 1 (5000) |
| DESeq2 | `trend_additive` | `stage_num` | stage trend | abundance | 2 (20387) | 1 (20191) |
| DESeq2 | `factor_interaction` | `union` | sex × stage, either coefficient (concordance quantity) | exported | 0 (5000) | 0 (5000) |
| DESeq2 | `factor_interaction` | `union` | sex × stage, either coefficient (concordance quantity) | abundance | 0 (20387) | 0 (20191) |
| limma | `factor_additive` | `stage_joint` | E1's limma quantity: joint F of stage | exported | 1 (5000) | 1 (5000) |
| limma | `factor_additive` | `stage_joint` | E1's limma quantity: joint F of stage | abundance | 1 (20387) | 1 (20191) |
| limma | `factor_additive` | `union` | stage, either moderated t | exported | 1 (5000) | 1 (5000) |
| limma | `factor_additive` | `union` | stage, either moderated t | abundance | 1 (20387) | 1 (20191) |
| limma | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | exported | 1 (5000) | 1 (5000) |
| limma | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | abundance | 1 (20387) | 1 (20191) |
| limma | `trend_additive` | `stage_num` | stage trend | exported | 1 (5000) | 1 (5000) |
| limma | `trend_additive` | `stage_num` | stage trend | abundance | 1 (20387) | 1 (20191) |
| limma | `factor_interaction` | `interaction_joint` | joint F of sex × stage (concordance quantity) | exported | 0 (5000) | 0 (5000) |
| limma | `factor_interaction` | `interaction_joint` | joint F of sex × stage (concordance quantity) | abundance | 0 (20387) | 0 (20191) |
| limma | `factor_interaction` | `union` | sex × stage, either moderated t | exported | 0 (5000) | 0 (5000) |
| limma | `factor_interaction` | `union` | sex × stage, either moderated t | abundance | 0 (20387) | 0 (20191) |

| analysis | DESeq2 interaction count (exported) | DESeq2 cell | limma interaction count (exported) | limma cell | delta rejected | methods agree |
|---|---:|---|---:|---|:---:|:---:|
| primary | 0 | concordant | 0 | concordant | no | yes |
| excluding_severely_affected | 0 | concordant | 0 | concordant | no | yes |

#### First latent axis of the variant (generated from `diagnostics/cohort_30/axis_correlations.csv`, `pls_0`, the 12 largest |ρ| and stage)

| variable | source | Spearman | partial Spearman |
|---|---|---:|---:|
| `ribosomal_mean` | donor_technical.csv | -0.729 | -0.723 |
| `ribosomal_pseudobulk` | donor_technical.csv | -0.700 | -0.692 |
| `tmm_rna` | manifest.json (export, tmm) | -0.627 | -0.596 |
| `atac_tss_enrichment_score_mean` | donor_technical.csv | -0.606 | -0.515 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | donor_technical.csv | -0.599 | -0.515 |
| `stage` | metadata.csv | 0.578 | 0.597 |
| `fraction_mitochondrial_umis_mean` | donor_technical.csv | 0.557 | 0.513 |
| `mitochondrial_mean` | donor_technical.csv | 0.557 | 0.513 |
| `ribosomal_sd` | donor_technical.csv | -0.552 | -0.542 |
| `mitochondrial_pseudobulk` | donor_technical.csv | 0.523 | 0.496 |
| `gex_sequencing_saturation_mean` | donor_technical.csv | 0.520 | 0.501 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | -0.454 | -0.427 |

#### Comparison (pasted from `results/motco_exclusion_sensitivity/excluding_severely_affected/comparison/cohort_30/exclusion_sensitivity.md`)

#### Donor-exclusion sensitivity: astrocyte, cohort_30, excluding_severely_affected

The `primary` row is the primary analysis (CV-selected rank 1, 80 donors), read from its `result.json` and never rerun; it stays the primary result. The `excluding_severely_affected` row reruns the whole export and the primary analysis on the cohort donors minus those SEA-AD flagged as severely affected (Gabitto et al. 2024: QC-outlying high-ADNC donors excluded from the gene-expression tests): 10 cohort donor(s) excluded, the cohort thresholds unchanged and k recomputed from the remaining cells by the export's rule (k = 9 → 9). The variant is a sensitivity analysis, never a revised primary: the excluded donors carry the most severe pathology, so the exclusion truncates the late epoch along with the quality tail, and the flag was made from QC metrics, so the variant conditions on quality by donor selection. A rejection status that differs between the rows is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in both rows. 999 RRPP permutations with seed 20261019 in both rows; α = 0.05 (rejected when p ≤ α). MOTCO's operating characteristics were measured on neither cohort.

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
| primary | 80 | 1 | 14.941 | 0.3 | no | 11.2868 | 8.27982 | 26.5775 |
| excluding_severely_affected | 70 | 1 | 6.35558 | 0.534 | no | 7.89107 | 5.76103 | 19.0545 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 80 | 1 | 0 | 1 | no | 3.06306 | 23.2919 | 0 |
| excluding_severely_affected | 70 | 1 | 0 | 1 | no | 0 | 0 | 0 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 80 | 1 | 0.608832 | 0.333 | no | 0.489849 | 0.355498 | 1.18287 |
| excluding_severely_affected | 70 | 1 | 0.495893 | 0.173 | no | 0.296309 | 0.217351 | 0.736906 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 80 | 0 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=28; Male/stage 0=9; Male/stage 1=10; Male/stage 2=12 | cv | 1 | 0.6313 | 61.9 | `results/motco_analysis/cohort_30` |
| excluding_severely_affected | 70 | 10 | 9 | Female/stage 0=11; Female/stage 1=10; Female/stage 2=20; Male/stage 0=9; Male/stage 1=10; Male/stage 2=10 | cv | 1 | 0.5559 | 108.3 | `motco_analysis/cohort_30` |

Excluded donors: `H20.33.028`, `H20.33.037`, `H20.33.038`, `H20.33.046`, `H21.33.009`, `H21.33.020`, `H21.33.033`, `H21.33.034`, `H21.33.039`, `H21.33.045`.

##### Axis check

Predeclared rule (`DONOR_EXCLUSION_AXIS_RULE`): along the variant analysis's first latent component, |Spearman| with stage exceeds the largest |Spearman| with a technical variable of POSITIVE_CONTROL_TECHNICAL_SOURCES (technical table columns for the variant donors; the variant export's TMM factors and log10 library sizes); the control cell type's value is the expectation, the reference's descriptive. A comparison, not a test; no p-value.

| cell type | role | reading | component | donors | ρ stage | largest technical | source | ρ technical | technical compared | holds |
|---|---|---|---|---:|---:|---|---|---:|---:|:---:|
| astrocyte | reference | descriptive | `pls_0` | 70 | 0.5784 | `ribosomal_mean` | donor_technical.csv | -0.7286 | 39 | no |

##### Reading

Predeclared (`DONOR_EXCLUSION_READING`): the flag check and the axis check are read together: if the flagged donors carry the quality tail (high AUC for the mitochondrial fraction and the primary pls_0) and the axis check holds in the control cell type, the quality gradient of the primary analyses sits in the donors SEA-AD excluded from expression tests, and excluding them is a donor-selection route to the same end as the latent adjustment; if the axis check fails although the flagged donors carry the tail, the gradient extends beyond them and donor selection does not resolve it; if the flagged donors do not carry the tail, the flag is not the gradient's carrier and the variant is an ordinary cohort sensitivity. In every case the primary analyses and the primary-basis positive control stay as they are; a rejection status that differs between the variant and the primary is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside angle; the non-rejection rules hold in the variant.

##### Findings

- axis check does not hold in cohort_30 (descriptive, role reference): |Spearman| of pls_0 with stage 0.5784 vs 0.7286 with ribosomal_mean (donor_technical.csv)
- excluding_severely_affected export: rna: 29 of 36601 genes are absent from GENCODE v32 basic; they are kept (chromosome 'unmapped')

### 4.2 cohort_50 (≥ 50 astrocyte nuclei in both omics, n = 74) — sensitivity cohort

#### Overview (generated from `exclusion_sensitivity.csv`)

| analysis | donors | excluded | k | rank | cv_mean_auroc | delta (p) | angle (p) | shape (p) | rejected |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| primary | 74 | 0 | 9 | 1 | 0.6616 | 18.7091 (0.213) | 0 (1) | 0.56471 (0.377) | none |
| excluding_severely_affected | 64 | 10 | 9 | 1 | 0.5661 | 7.48789 (0.476) | 0 (1) | 0.384785 (0.296) | none |

#### Axis check (from `axis_check.csv`)

Along the variant's `pls_0` (n = 64): |ρ| with stage 0.5955 (ρ = 0.5955) against the largest of 39 technical variables, `ribosomal_mean` (donor_technical.csv) |ρ| 0.6716 (ρ = -0.6716): **does not hold** (descriptive, role reference).

#### Stage DE on the exclusion cohort, RNA (from `stage_de_summary.csv` and `stage_de_limma_summary.csv`; counts at FDR 0.05)

| method | model | coefficient / test | meaning | set | primary (n tested) | variant (n tested) |
|---|---|---|---|---|---:|---:|
| DESeq2 | `factor_additive` | `union` | E1's DESeq2 quantity: stage, either coefficient | exported | 2 (5000) | 1 (5000) |
| DESeq2 | `factor_additive` | `union` | E1's DESeq2 quantity: stage, either coefficient | abundance | 2 (20238) | 1 (20037) |
| DESeq2 | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | exported | 2 (5000) | 1 (5000) |
| DESeq2 | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | abundance | 2 (20238) | 1 (20037) |
| DESeq2 | `trend_additive` | `stage_num` | stage trend | exported | 2 (5000) | 1 (5000) |
| DESeq2 | `trend_additive` | `stage_num` | stage trend | abundance | 2 (20238) | 1 (20037) |
| DESeq2 | `factor_interaction` | `union` | sex × stage, either coefficient (concordance quantity) | exported | 0 (5000) | 0 (5000) |
| DESeq2 | `factor_interaction` | `union` | sex × stage, either coefficient (concordance quantity) | abundance | 0 (20238) | 0 (20037) |
| limma | `factor_additive` | `stage_joint` | E1's limma quantity: joint F of stage | exported | 1 (5000) | 1 (5000) |
| limma | `factor_additive` | `stage_joint` | E1's limma quantity: joint F of stage | abundance | 2 (20238) | 1 (20037) |
| limma | `factor_additive` | `union` | stage, either moderated t | exported | 1 (5000) | 1 (5000) |
| limma | `factor_additive` | `union` | stage, either moderated t | abundance | 1 (20238) | 1 (20037) |
| limma | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | exported | 1 (5000) | 1 (5000) |
| limma | `factor_additive` | `stage_2_vs_0` | stage 2 vs 0 | abundance | 1 (20238) | 1 (20037) |
| limma | `trend_additive` | `stage_num` | stage trend | exported | 1 (5000) | 1 (5000) |
| limma | `trend_additive` | `stage_num` | stage trend | abundance | 1 (20238) | 1 (20037) |
| limma | `factor_interaction` | `interaction_joint` | joint F of sex × stage (concordance quantity) | exported | 0 (5000) | 0 (5000) |
| limma | `factor_interaction` | `interaction_joint` | joint F of sex × stage (concordance quantity) | abundance | 0 (20238) | 0 (20037) |
| limma | `factor_interaction` | `union` | sex × stage, either moderated t | exported | 0 (5000) | 1 (5000) |
| limma | `factor_interaction` | `union` | sex × stage, either moderated t | abundance | 0 (20238) | 1 (20037) |

| analysis | DESeq2 interaction count (exported) | DESeq2 cell | limma interaction count (exported) | limma cell | delta rejected | methods agree |
|---|---:|---|---:|---|:---:|:---:|
| primary | 0 | concordant | 0 | concordant | no | yes |
| excluding_severely_affected | 0 | concordant | 0 | concordant | no | yes |

#### First latent axis of the variant (generated from `diagnostics/cohort_50/axis_correlations.csv`, `pls_0`, the 12 largest |ρ| and stage)

| variable | source | Spearman | partial Spearman |
|---|---|---:|---:|
| `ribosomal_mean` | donor_technical.csv | -0.672 | -0.638 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | donor_technical.csv | -0.652 | -0.615 |
| `ribosomal_pseudobulk` | donor_technical.csv | -0.639 | -0.614 |
| `fraction_mitochondrial_umis_mean` | donor_technical.csv | 0.629 | 0.553 |
| `mitochondrial_mean` | donor_technical.csv | 0.629 | 0.553 |
| `atac_tss_enrichment_score_mean` | donor_technical.csv | -0.627 | -0.589 |
| `gex_sequencing_saturation_mean` | donor_technical.csv | 0.620 | 0.586 |
| `tmm_rna` | manifest.json (export, tmm) | -0.613 | -0.566 |
| `mitochondrial_pseudobulk` | donor_technical.csv | 0.599 | 0.532 |
| `stage` | metadata.csv | 0.595 | 0.595 |
| `fraction_mitochondrial_umis_sd` | donor_technical.csv | 0.495 | 0.465 |
| `mitochondrial_sd` | donor_technical.csv | 0.495 | 0.465 |

#### Comparison (pasted from `results/motco_exclusion_sensitivity/excluding_severely_affected/comparison/cohort_50/exclusion_sensitivity.md`)

#### Donor-exclusion sensitivity: astrocyte, cohort_50, excluding_severely_affected

The `primary` row is the primary analysis (CV-selected rank 1, 74 donors), read from its `result.json` and never rerun; it stays the primary result. The `excluding_severely_affected` row reruns the whole export and the primary analysis on the cohort donors minus those SEA-AD flagged as severely affected (Gabitto et al. 2024: QC-outlying high-ADNC donors excluded from the gene-expression tests): 10 cohort donor(s) excluded, the cohort thresholds unchanged and k recomputed from the remaining cells by the export's rule (k = 9 → 9). The variant is a sensitivity analysis, never a revised primary: the excluded donors carry the most severe pathology, so the exclusion truncates the late epoch along with the quality tail, and the flag was made from QC metrics, so the variant conditions on quality by donor selection. A rejection status that differs between the rows is cohort dependence, not a result of the variant; a CV-selected rank that differs is reported beside `angle`, not read as better or worse; `angle` and `shape` non-rejections are uninformative and a `delta` rejection is a difference, not specifically of size, in both rows. 999 RRPP permutations with seed 20261019 in both rows; α = 0.05 (rejected when p ≤ α). MOTCO's operating characteristics were measured on neither cohort.

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
| primary | 74 | 1 | 18.7091 | 0.213 | no | 11.6918 | 9.218 | 30.0679 |
| excluding_severely_affected | 64 | 1 | 7.48789 | 0.476 | no | 8.65938 | 6.61897 | 21.2566 |

##### angle (degrees)

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 74 | 1 | 0 | 1 | no | 1.98198 | 18.7931 | 0 |
| excluding_severely_affected | 64 | 1 | 0 | 1 | no | 0 | 0 | 0 |

##### shape

| analysis | donors | rank | observed | p (RRPP) | rejected | null mean | null sd | null q95 |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| primary | 74 | 1 | 0.56471 | 0.377 | no | 0.495877 | 0.343673 | 1.159 |
| excluding_severely_affected | 64 | 1 | 0.384785 | 0.296 | no | 0.291485 | 0.225515 | 0.738017 |

##### Runs

| analysis | donors | excluded | k | cells | selection | rank | cv_mean_auroc | runtime (s) | directory |
|---|---:|---:|---:|---|---|---:|---:|---:|---|
| primary | 74 | 0 | 9 | Female/stage 0=10; Female/stage 1=9; Female/stage 2=25; Male/stage 0=9; Male/stage 1=9; Male/stage 2=12 | cv | 1 | 0.6616 | 62.4 | `results/motco_analysis/cohort_50` |
| excluding_severely_affected | 64 | 10 | 9 | Female/stage 0=10; Female/stage 1=9; Female/stage 2=17; Male/stage 0=9; Male/stage 1=9; Male/stage 2=10 | cv | 1 | 0.5661 | 100.2 | `motco_analysis/cohort_50` |

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

## 5. Cross-cell-type comparisons (compared, not tested)

Predeclared in the proposal and computed from the generated tables (D6): the variant CV mean AUROC of microglia against astrocytes, and each cell type's exported RNA stage counts on the exclusion cohort against its primary counts (DESeq2 `factor_additive` union, limma `stage_joint`, E1's two quantities). The microglia variant AUROC exceeds the astrocyte variant AUROC in both cohorts (0.6388 against 0.5559; 0.6280 against 0.5661), as the primary AUROCs did in `cohort_30` (0.6921 against 0.6313) but not in `cohort_50` (0.6414 against 0.6616). The exclusion lowers every AUROC. The microglia DESeq2 counts keep a third (`cohort_30`) and under a tenth (`cohort_50`) of the primary's; limma's joint F, which already counted few microglia stage features on the primary (2 and 0), is unchanged; the astrocyte counts are 1 or 2 throughout.

| cohort | microglia variant cv_mean_auroc | astrocyte variant cv_mean_auroc | microglia > astrocyte | microglia primary | astrocyte primary |
|---|---:|---:|:---:|---:|---:|
| cohort_30 | 0.6388 | 0.5559 | yes | 0.6921 | 0.6313 |
| cohort_50 | 0.6280 | 0.5661 | yes | 0.6414 | 0.6616 |

| cell type | cohort | method | quantity (exported RNA) | primary | variant | variant / primary |
|---|---|---|---|---:|---:|---:|
| microglia | cohort_30 | DESeq2 | `factor_additive` `union` | 61 | 20 | 0.33 |
| microglia | cohort_30 | limma | `factor_additive` `stage_joint` | 2 | 2 | 1.00 |
| microglia | cohort_50 | DESeq2 | `factor_additive` `union` | 35 | 3 | 0.09 |
| microglia | cohort_50 | limma | `factor_additive` `stage_joint` | 0 | 0 | n/a |
| astrocyte | cohort_30 | DESeq2 | `factor_additive` `union` | 1 | 1 | 1.00 |
| astrocyte | cohort_30 | limma | `factor_additive` `stage_joint` | 1 | 1 | 1.00 |
| astrocyte | cohort_50 | DESeq2 | `factor_additive` `union` | 2 | 1 | 0.50 |
| astrocyte | cohort_50 | limma | `factor_additive` `stage_joint` | 1 | 1 | 1.00 |

## 6. Reading (predeclared, in the outcome that occurred)

The reading was written in the proposal before any run (`DONOR_EXCLUSION_READING`): the flag check and the axis check are read together; if the flagged donors carry the quality tail (high AUC for the mitochondrial fraction and the primary `pls_0`) and the axis check holds in microglia, the quality gradient of the primary analyses sits in the donors SEA-AD excluded from expression tests, and excluding them is a donor-selection route to the same end as the latent adjustment; if the axis check fails although the flagged donors carry the tail, the gradient extends beyond them and donor selection does not resolve it; if the flagged donors do not carry the tail, the flag is not the gradient's carrier and the variant is an ordinary cohort sensitivity.

**The first branch occurred in microglia, in both cohorts.** The flagged donors carry the tail (AUC on the primary `pls_0` 0.986 / 0.983, within stage 2 0.964 / 0.955; on `mitochondrial_pseudobulk` 0.961 / 0.945, within stage 2 0.918 / 0.882), and on the cohort without them the first latent component orders donors by stage first (|ρ| 0.8201 / 0.8248 against at most 0.6891 / 0.6746 for a technical variable). So the quality gradient that made E3 fail sits in the donors SEA-AD excluded from its expression tests, and removing them reaches what the latent adjustment reached by conditioning the basis: an axis led by stage. The two routes are not the same operation and their stage |ρ| differ (0.8201 here, 0.5279 on the `latent_adjusted` basis, 0.688 on the primary); neither is promoted. Three qualifications bind this reading:

- The gradient is not confined to the flagged donors along the primary axis: among the unflagged donors the primary `pls_0` still correlates more with `mitochondrial_pseudobulk` (0.799 / 0.733) than with stage (0.643 / 0.699), and the mitochondrial fraction still rises with stage among them (0.415 / 0.355). What changes is the refit: without the tail, PLS-DA's first component on the variant blocks is led by stage, and the mitochondrial measures follow at 0.66 to 0.69, above the diagnostic's 0.3 flag, which is a measurement-space finding about the variant.
- The stage signal shrinks with the excluded donors: the exported RNA stage features under DESeq2 fall from 61 to 20 and from 35 to 3, and limma's stage 2 vs 0 t from 28 to 2 and from 36 to 0. All flagged donors are stage 2, so their removal truncates the late epoch with the quality tail, and this analysis cannot tell how much of the lost per-feature stage signal was quality and how much was late pathology (the predeclared limitation). The latent separation also falls (CV AUROC 0.6921 → 0.6388, 0.6414 → 0.6280).
- The microglia `cohort_50` variant selects rank 4, so its `angle` (19.14°, p 0.45) carries orientation information there and none in the rank-1 rows; this is reported beside `angle` and is cohort dependence, not a result.

In astrocytes (descriptive) the flagged donors carry the tail as well (primary `pls_0` AUC 0.961 / 0.959, `mitochondrial_pseudobulk` 0.966 / 0.966, `ribosomal_mean` 0.076 / 0.069), but the axis check does not hold: the variant's first component follows `ribosomal_mean` (0.7286 / 0.6716) more than stage (0.5784 / 0.5955). For the reference cell type that is the second branch read descriptively: its quality axis extends beyond the flagged donors (among the unflagged donors the primary `pls_0` keeps |ρ| 0.791 / 0.780 with `ribosomal_mean` against 0.384 / 0.392 with stage), and donor selection does not resolve it. The astrocyte stage counts stay at 1 to 2 features, as on the primary. One interaction count moves: on the astrocyte `cohort_50` variant limma's per-coefficient union of `factor_interaction` holds one exported feature (`AL160408.1`, `male_x_stage_2` moderated t, adj. p 0.040; 0 on the primary), while limma's joint F (`interaction_joint`, the concordance quantity) and DESeq2 (smallest padj 0.14) detect none; it is method dependence on a union count that carries no joint error-rate claim, and every concordance cell stays concordant.

In every case, as predeclared: the primary analyses of both cell types, the primary cohorts and the primary-basis positive control stay as they are; no rejection status differs between a variant and its primary (nothing is rejected in any row), so there is no cohort dependence of a rejection to report; the non-rejections of `angle` and `shape` stay uninformative; every concordance cell is concordant under both methods, on the variants as on the primaries.

## 7. Limitations

- **Late-epoch truncation.** The excluded donors carry the most severe pathology (all stage 2, Braak V or VI, 10 of 11 with dementia), so the exclusion truncates the late epoch along with the quality tail, and a stage effect concentrated in them is removed with them; the fall of the microglia stage counts is read in that light.
- **The flag was made from QC metrics.** SEA-AD derived it from a PCA of library QC metrics, so the variant conditions on quality by donor selection rather than by regression; that the axis check then holds is partly built into the selection.
- **Smaller cells, less power.** The stage-2 cells shrink (female astrocyte `cohort_30` 28 → 20, microglia 24 → 19; male 12 → 10 and 11 → 9), so the variants have less power than the primaries, and the absent rejections say less.
- **No operating characteristics.** MOTCO's operating characteristics were measured on neither cohort's blocks, so no variant carries a Type I or power claim of its own; the microglia `cohort_50` variant selects rank 4, above `n_stages − 1`, where MOTCO's forced-rank studies report a wider `angle` null; a CV-selected rank 4 carries no calibration of its own.
- **The flag's sex split.** 8 of the 11 flagged donors are female (5 of 7 in microglia `cohort_30`), so the exclusion changes the female stage-2 cell more than the male one, which bears on the sex × stage comparison itself.
- **One flag, one release.** The flag is SEA-AD's, taken as recorded (11 donors, as published); the recorded staging labels and the pseudo-progression score are read by nothing here (CPS-based staging is a later proposal).

## 8. Reproduction

```
# cluster login node (obs only)
sead-sexdiff donor-flags --omic both
sead-sexdiff donor-flags --cell-type microglia --omic both
# workstation, after scripts/sync_pseudobulk.sh [--cell-type microglia]
sead-sexdiff exclusion-sensitivity --cohort all --n-jobs 1
sead-sexdiff exclusion-sensitivity --cell-type microglia --cohort all --n-jobs 1
V=results/motco_exclusion_sensitivity/excluding_severely_affected   # results/microglia/... with --cell-type microglia
sead-sexdiff axis-diagnostics --input-dir $V/motco_input --analysis-dir $V/motco_analysis --rank-dir $V/motco_rank_sensitivity --out-dir $V/diagnostics
sead-sexdiff stage-de --omic both --input-dir $V/motco_input --analysis-dir $V/motco_analysis --out-dir $V/stage_de --n-cpus 10
sead-sexdiff stage-de-limma --omic rna --input-dir $V/motco_input --analysis-dir $V/motco_analysis --deseq2-dir $V/stage_de --out-dir $V/stage_de_limma
```
