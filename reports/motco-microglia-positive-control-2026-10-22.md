# MOTCO microglia positive control: does the pipeline resolve ADNC stage on real pseudobulks? — 2026-10-22

> OpenSpec change `mtg-microglia-positive-control` (spec `cell-type-positive-control`). The astrocyte MOTCO analyses are non-rejections with a weak, rank-1 latent space, and nothing had shown that the pipeline, run on real SEA-AD donor pseudobulks, builds a latent space that resolves ADNC stage at all. Microglia carry the strongest reported stage response in the SEA-AD MTG, so the unchanged pipeline was run on microglia pseudobulks and three expectations, written into the code before any microglia run (`POSITIVE_CONTROL_EXPECTATIONS`), were evaluated from the generated files by the `positive-control` step. Every number below is pasted by script from the generated files (`positive_control.csv` / `.md` / manifest, the stage-DE summaries of both cell types and both methods, `summary.md` and `result.json`, `diagnostics.md` and `axis_correlations.csv`, the cohort summaries, the manifests and the metadata-only spike's tables) and was checked against them before commit. The hand-written text carries only the reading, under the rules below.

**Reading rules (binding).** This is a positive control for the stage signal, not for the sex × stage hypothesis. If it holds, the pipeline resolves stage structure on real data and the astrocyte null is not a pipeline failure; if it does not, the latent construction is in question for any cell type, and that is reported. The expectations read no p-value, rejection status or trajectory statistic; their values are compared, not tested, and the expectations are not changed from results. The microglia sex × stage analysis is a second case-study result read under the existing rules: a `delta` rejection is a difference, not specifically of size; `angle` and `shape` non-rejections are uninformative; the CV-selected rank beside `angle`; MOTCO's operating characteristics were measured on neither cell type's blocks. A stage-DE call that differs between the two methods is method dependence, a result of neither. Nothing found in microglia changes any astrocyte output or reading.

**Result.** The predeclared positive control does **not** hold: in the primary cohort (`cohort_30`) E1 and E2 hold but E3 fails (not holding: E3). Microglia carry more per-feature stage signal than astrocytes under both stage-DE methods (exported RNA genes with a stage effect at FDR 0.05: DESeq2 61 against 1; limma 2 against 1), and the microglia latent space separates the stages somewhat better (CV mean AUROC 0.6921 against 0.6313). But the first microglia latent axis follows the technical variables more closely than stage: its |Spearman| with stage is 0.688, below the 0.844 of `mitochondrial_pseudobulk`. In `cohort_50` (descriptive) nothing holds: E1 is `method_dependent` (DESeq2 35 against 2; limma 0 against 1), E2 fails (0.6414 against 0.6616) and E3 fails (0.731 against 0.785, `mitochondrial_pseudobulk`). The microglia sex × stage analysis rejects nothing (`delta` p 0.391 and 0.227; CV-selected rank 1 and 1), and every stage-DE concordance cell of both methods is concordant (section 4).

## 1. Provenance

- Code: commit `ea6bc4f` (`ea6bc4f99430b5ca0983a884d03056410c4c31f7`) for every step. Cluster steps (`pseudobulk`, `covariates`, `technical-covariates`) ran on the `cluster.ing.uc.cl` login node from a clean tree (dirty = False); the workstation steps ran on `nixos-home` at the same commit (`ea6bc4f`, dirty = True: the tracked changes at run time were the task checkboxes in `tasks.md` and an IDE file (`.idea/sead-sexdiff.iml`); `git diff ea6bc4f -- src scripts` was empty).
- Software (workstation): MOTCO 0.9.0, scikit-learn 1.8.0, pydeseq2 0.5.3, Python 3.12.15, numpy 2.3.5, pandas 2.3.3; R 4.6.1, limma 3.68.4, edgeR 4.10.1. Cluster: Python 3.13.15, anndata 0.12.7.
- Runs: `covariates --cell-type microglia --omic both` on the cluster login node, 2026-10-22 18:32Z; `pseudobulk --cell-type microglia --omic both` 2026-10-22T18:59:38Z → `technical-covariates --cell-type microglia --omic both` 2026-10-22T19:05:00Z → end 2026-10-22T19:07:38Z (login node, nohup); on `nixos-home` `export`, `analyze` (`n_jobs = 1`; cohort_30 88.9 s, cohort_50 84.1 s), `axis-diagnostics`, `stage-de --n-cpus 10`, `stage-de-limma` and `positive-control`, all with `--cell-type microglia` defaults, 2026-10-22T19:10:23Z → 2026-10-22T19:38:55Z.
- Microglia pseudobulk (`results/microglia/pseudobulk/`): `manifest.json` sha256 `49b937223a485a5233fcc7953fc51f9cd343ba65ccbc5c2dc2c3800ce8d8cb1c`, `rna_counts.h5ad` `450b72285cf35525d5670d3297d5129fa6149002b22dc64551accc5333ba3360`, `atac_counts.h5ad` `1fd7a67bf676b32f707a6662d079a73e141171a3f15b0f8d103a420b0ef00f9d`, `donor_covariates.csv` `b14398d88cd8c60c0a64c37d1229e2281ff26e847301039fcc228ebd46c453d7`, `donor_technical.csv` `ac54eef6a253135db7260eec0365b837c103495f1d3e1041c9b9e06ff288b5cd`; synced from the cluster and checked by sha256 (11 files).
- Astrocyte reference (`results/`, read by `positive-control`, never rerun): cohort_30: `motco_analysis/cohort_30/result.json` `1d6154b02d77…`, `stage_de/cohort_30/stage_de_summary.csv` `6b4947e8db81…`, `stage_de_limma/cohort_30/stage_de_limma_summary.csv` `fdd010a57299…`; cohort_50: `motco_analysis/cohort_50/result.json` `4c6e7e58e39d…`, `stage_de/cohort_50/stage_de_summary.csv` `9f33acf5b59e…`, `stage_de_limma/cohort_50/stage_de_limma_summary.csv` `c63ab8f92d8e…`.
- Deviations from the tasks as written: (1) the metadata-only spike read the 170 objects from the public bucket over HTTPS range requests while the cluster was unreachable, each object's size and ETag checked against the astrocyte pseudobulk manifest (all matched; the spike's astrocyte tallies reproduce every astrocyte `donor_table.csv` count); (2) the microglia `pseudobulk` and `technical-covariates` ran with nohup on the cluster login node instead of `scripts/pseudobulk.sbatch` / `scripts/technical.sbatch`: Slurm job 925058 spent 28 min in NFS waits on n11 (the login node had rebooted at 15:24 local time) without finishing a donor, and was cancelled with its chained technical job. The same code and arguments ran; only the host differs.
- Astrocyte outputs untouched: on the cluster, the 367 files under `results/pseudobulk/` hash as before the microglia runs (`sha256sum -c`); on `nixos-home`, the 167 astrocyte files under `results/{pseudobulk,motco_input,motco_analysis,motco_diagnostics,stage_de,stage_de_limma}/` hash as before the microglia chain (`sha256sum -c`).

## 2. The microglia population and cohorts

Selection (registry `CELL_TYPES["microglia"]`, from the metadata-only spike in `openspec/changes/archive/2026-10-22-mtg-microglia-positive-control/spike.md`): RNA `Subclass` `Immune`, ATAC `Subclass` `Microglia-PVM`, the six `Micro-PVM_*` supertypes of both releases (no release has a PVM-only supertype), `Lymphocyte` and `Monocyte` not selected; `Reference` donors excluded; RNA nuclei with ≥ 200 genes. Nuclei per supertype over all donors (spike, Reference included):

| supertype | RNA (Immune) | ATAC (Microglia-PVM) | selected |
|---|---:|---:|:---:|
| `Lymphocyte` | 767 | 219 | no |
| `Micro-PVM_1` | 1383 | 447 | yes |
| `Micro-PVM_2` | 21090 | 6130 | yes |
| `Micro-PVM_2_1-SEAAD` | 282 | 35 | yes |
| `Micro-PVM_2_3-SEAAD` | 8278 | 1796 | yes |
| `Micro-PVM_3-SEAAD` | 4468 | 910 | yes |
| `Micro-PVM_4-SEAAD` | 406 | 69 | yes |
| `Monocyte` | 12 | 57 | no |

Nuclei per donor (non-Reference donors; pasted from the two `donor_table.csv`):

| cell type | omic | donors | min | median | max | below 30 | below 50 |
|---|---|---:|---:|---:|---:|---:|---:|
| microglia | rna | 84 | 63 | 376 | 1034 | 0 | 0 |
| microglia | atac | 84 | 11 | 88 | 478 | 8 | 21 |
| astrocyte | rna | 84 | 76 | 667.5 | 2266 | 0 | 0 |
| astrocyte | atac | 84 | 20 | 159 | 1348 | 4 | 10 |

Cohorts (`cohort_summary.csv` of each cell type; donors per sex × stage, stages 0 / 1 / 2):

| cohort | sex | microglia | astrocyte |
|---|---|---|---|
| cohort_30 | Female | 11 / 10 / 24 (n = 45) | 11 / 10 / 28 (n = 49) |
| cohort_30 | Male | 9 / 11 / 11 (n = 31) | 9 / 10 / 12 (n = 31) |
| cohort_30 | all | 20 / 21 / 35 (n = 76) | 20 / 20 / 40 (n = 80) |
| cohort_50 | Female | 8 / 10 / 18 (n = 36) | 10 / 9 / 25 (n = 44) |
| cohort_50 | Male | 8 / 10 / 9 (n = 27) | 9 / 9 / 12 (n = 30) |
| cohort_50 | all | 16 / 20 / 27 (n = 63) | 19 / 18 / 37 (n = 74) |

The microglia cohorts equal the spike's preregistered cells (`manifest.json` findings: none); no cell is empty, so the astrocyte thresholds were kept. `k` (smallest cell): cohort_30 9, cohort_50 8. Donors in both cell types' cohorts: cohort_30 73, cohort_50 60.

## 3. Positive control: E1–E3

### cohort_30 (primary cohort; pasted from `results/microglia/positive_control/cohort_30/positive_control.md`)

| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |
|---|---|---|---|---|---|---|---|---|---|
| E1 | deseq2 | exported RNA features at FDR 0.05: factor_additive union | microglia | 61 | astrocyte | 1 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | limma | exported RNA features at FDR 0.05: factor_additive stage_joint | microglia | 2 | astrocyte | 1 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | both | E1 under both stage-DE methods |  |  |  |  | yes | holds | holds under both methods |
| E2 |  | CV mean AUROC of the primary analysis (latent.cv_mean_auroc) | microglia | 0.692109 | astrocyte | 0.631278 | yes | holds |  |
| E3 |  | \|Spearman\| of pls_0 (primary) with stage vs the largest with a technical variable | stage | 0.688175 | mitochondrial_pseudobulk | 0.843828 | no | fails | 39 technical variables compared (donor_technical.csv, manifest.json (export, tmm)) |
| overall |  | E1, E2 and E3 all hold |  |  |  |  | no | does_not_hold | the primary cohort: this row is the control's outcome; not holding: E3 |

### cohort_50 (descriptive; pasted from `results/microglia/positive_control/cohort_50/positive_control.md`)

| Expectation | Method | Quantity | Control | Value | Compared with | Value | Holds | Outcome | Detail |
|---|---|---|---|---|---|---|---|---|---|
| E1 | deseq2 | exported RNA features at FDR 0.05: factor_additive union | microglia | 35 | astrocyte | 2 | yes | holds | tested 5000 (control), 5000 (reference) |
| E1 | limma | exported RNA features at FDR 0.05: factor_additive stage_joint | microglia | 0 | astrocyte | 1 | no | fails | tested 5000 (control), 5000 (reference) |
| E1 | both | E1 under both stage-DE methods |  |  |  |  | no | method_dependent | held under deseq2 only; counts as not holding |
| E2 |  | CV mean AUROC of the primary analysis (latent.cv_mean_auroc) | microglia | 0.641431 | astrocyte | 0.661645 | no | fails |  |
| E3 |  | \|Spearman\| of pls_0 (primary) with stage vs the largest with a technical variable | stage | 0.730773 | mitochondrial_pseudobulk | 0.784562 | no | fails | 39 technical variables compared (donor_technical.csv, manifest.json (export, tmm)) |
| overall |  | E1, E2 and E3 all hold |  |  |  |  | no | does_not_hold | not the primary cohort (cohort_30); descriptive; not holding: E1, E2, E3 |

Findings of the step (both cohorts): 6 astrocyte manifests predate the registry and read as astrocyte; other findings: none.

## 4. Stage DE of both cell types under both methods

Counts of features at FDR 0.05 (DESeq2 `padj`, limma adjusted p) on the exported set (MOTCO's 5,000 features) and, after the slash, on the abundance-passing set; pasted from `stage_de_summary.csv` and `stage_de_limma_summary.csv` of `results/` (astrocyte) and `results/microglia/`. Stage rows first; E1 reads the RNA `factor_additive` `union` (DESeq2) and `stage_joint` (limma) rows, exported set.

### cohort_30

| omic | model | DESeq2 test | limma test | DESeq2 microglia | DESeq2 astrocyte | limma microglia | limma astrocyte |
|---|---|---|---|---:|---:|---:|---:|
| rna | `factor_additive` | `stage_1_vs_0` | `stage_1_vs_0` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| rna | `factor_additive` | `stage_2_vs_0` | `stage_2_vs_0` | 61 / 156 | 1 / 1 | 28 / 64 | 1 / 1 |
| rna | `factor_additive` | `union` | `stage_joint` | 61 / 156 | 1 / 1 | 2 / 2 | 1 / 1 |
| rna | `trend_additive` | `stage_num` | `stage_num` | 57 / 145 | 2 / 2 | 31 / 67 | 1 / 1 |
| rna | `factor_interaction` | `union` | `interaction_joint` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| rna | `trend_interaction` | `male_x_stage_num` | `male_x_stage_num` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| atac | `factor_additive` | `stage_1_vs_0` | `stage_1_vs_0` | 0 / 1 | 0 / 0 | 0 / 1 | 0 / 0 |
| atac | `factor_additive` | `stage_2_vs_0` | `stage_2_vs_0` | 20 / 113 | 0 / 22 | 14 / 267 | 0 / 9 |
| atac | `factor_additive` | `union` | `stage_joint` | 20 / 114 | 0 / 22 | 3 / 49 | 0 / 0 |
| atac | `trend_additive` | `stage_num` | `stage_num` | 23 / 116 | 0 / 23 | 14 / 218 | 0 / 0 |
| atac | `factor_interaction` | `union` | `interaction_joint` | 0 / 3 | 0 / 0 | 0 / 0 | 0 / 0 |
| atac | `trend_interaction` | `male_x_stage_num` | `male_x_stage_num` | 0 / 1 | 0 / 0 | 0 / 1 | 0 / 0 |

### cohort_50

| omic | model | DESeq2 test | limma test | DESeq2 microglia | DESeq2 astrocyte | limma microglia | limma astrocyte |
|---|---|---|---|---:|---:|---:|---:|
| rna | `factor_additive` | `stage_1_vs_0` | `stage_1_vs_0` | 1 / 1 | 0 / 0 | 0 / 0 | 0 / 0 |
| rna | `factor_additive` | `stage_2_vs_0` | `stage_2_vs_0` | 34 / 64 | 2 / 2 | 36 / 70 | 1 / 1 |
| rna | `factor_additive` | `union` | `stage_joint` | 35 / 65 | 2 / 2 | 0 / 0 | 1 / 2 |
| rna | `trend_additive` | `stage_num` | `stage_num` | 34 / 59 | 2 / 2 | 35 / 67 | 1 / 1 |
| rna | `factor_interaction` | `union` | `interaction_joint` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| rna | `trend_interaction` | `male_x_stage_num` | `male_x_stage_num` | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| atac | `factor_additive` | `stage_1_vs_0` | `stage_1_vs_0` | 0 / 3 | 0 / 0 | 0 / 0 | 0 / 0 |
| atac | `factor_additive` | `stage_2_vs_0` | `stage_2_vs_0` | 13 / 121 | 4 / 36 | 31 / 582 | 0 / 12 |
| atac | `factor_additive` | `union` | `stage_joint` | 13 / 124 | 4 / 36 | 10 / 136 | 0 / 0 |
| atac | `trend_additive` | `stage_num` | `stage_num` | 16 / 130 | 4 / 30 | 33 / 614 | 0 / 7 |
| atac | `factor_interaction` | `union` | `interaction_joint` | 0 / 3 | 0 / 0 | 0 / 1 | 0 / 0 |
| atac | `trend_interaction` | `male_x_stage_num` | `male_x_stage_num` | 0 / 3 | 0 / 0 | 0 / 2 | 0 / 0 |

### Concordance cells (`STAGE_DE_CONCORDANCE_RULE`, DESeq2, the reported one; `STAGE_DE_LIMMA_CONCORDANCE_RULE` beside it)

| cohort | omic | cell type | DESeq2 exported interaction count | DESeq2 cell | limma exported interaction count | limma cell | delta rejected | methods agree |
|---|---|---|---:|---|---:|---|:---:|:---:|
| cohort_30 | rna | microglia | 0 | concordant | 0 | concordant | no | yes |
| cohort_30 | rna | astrocyte | 0 | concordant | 0 | concordant | no | yes |
| cohort_30 | atac | microglia | 0 | concordant | 0 | concordant | no | yes |
| cohort_30 | atac | astrocyte | 0 | concordant | 0 | concordant | no | yes |
| cohort_50 | rna | microglia | 0 | concordant | 0 | concordant | no | yes |
| cohort_50 | rna | astrocyte | 0 | concordant | 0 | concordant | no | yes |
| cohort_50 | atac | microglia | 0 | concordant | 0 | concordant | no | yes |
| cohort_50 | atac | astrocyte | 0 | concordant | 0 | concordant | no | yes |

## 5. The microglia sex × stage analysis

### cohort_30 (pasted from `results/microglia/motco_analysis/cohort_30/summary.md`)

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 8.0642 | 0.391 | no | 7.23504 | 5.27064 | 17.4616 |  |
| angle | 0 | 1 | no | 0 | 0 | 0 | degrees; CV-selected rank 1 (at rank 1 the angle is 0° or 180° by construction) |
| shape | 0.154684 | 0.637 | no | 0.25893 | 0.192607 | 0.638734 |  |

- PLS-DA rank: 1 (`cv`)
- cv_mean_auroc: 0.6921

### cohort_50 (pasted from `results/microglia/motco_analysis/cohort_50/summary.md`)

| statistic | observed | p (RRPP) | rejected | null mean | null sd | null q95 | note |
|---|---:|---:|:---:|---:|---:|---:|---|
| delta | 12.4641 | 0.227 | no | 8.45222 | 5.99525 | 19.5118 |  |
| angle | 0 | 1 | no | 0 | 0 | 0 | degrees; CV-selected rank 1 (at rank 1 the angle is 0° or 180° by construction) |
| shape | 0.269568 | 0.399 | no | 0.255498 | 0.194219 | 0.624989 |  |

- PLS-DA rank: 1 (`cv`)
- cv_mean_auroc: 0.6414

Astrocyte primary analysis beside it (`results/motco_analysis/<cohort>/result.json`): cohort_30 rank 1, CV mean AUROC 0.6313, delta p 0.3; cohort_50 rank 1, CV mean AUROC 0.6616, delta p 0.213.

Covariates (both cell types): `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac` (the cell type's own nuclei counts), additive, missing values an error (none). Feature filtering per cohort (export manifest `feature_counts`, input → after chrX/Y/M → after CPM ≥ 1 in ≥ k → kept):

| cohort | omic | microglia | astrocyte |
|---|---|---|---|
| cohort_30 | rna | 36601 → 35329 → 20660 → 5000 | 36601 → 35329 → 20387 → 5000 |
| cohort_30 | atac | 218882 → 212816 → 130566 → 5000 | 218882 → 212816 → 128916 → 5000 |
| cohort_50 | rna | 36601 → 35329 → 20272 → 5000 | 36601 → 35329 → 20238 → 5000 |
| cohort_50 | atac | 218882 → 212816 → 128618 → 5000 | 218882 → 212816 → 124818 → 5000 |

## 6. What the first microglia latent axis tracks

Raw Spearman correlations of the primary `pls_0` with the diagnostic variables, and partial correlations given the six analysis covariates (pasted from `results/microglia/motco_diagnostics/<cohort>/diagnostics.md`, flagged rows at |ρ| ≥ 0.3, the first 12 by |ρ|). The diagnostic reads sex and stage, so it is a measurement-space finding and feeds no analysis.

### cohort_30

| analysis | rank | component | sd | variable | n | rho | rho \| covariates | flag |
|---|---:|---|---:|---|---:|---:|---:|---|
| primary | 1 | pls_0 | 24.24 | `mitochondrial_pseudobulk` | 76 | 0.844 | 0.801 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `fraction_mitochondrial_umis_mean` | 76 | 0.837 | 0.804 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `mitochondrial_mean` | 76 | 0.837 | 0.804 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `fraction_mitochondrial_umis_sd` | 76 | 0.828 | 0.817 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `mitochondrial_sd` | 76 | 0.828 | 0.817 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 76 | -0.725 | -0.685 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `gex_sequencing_saturation_mean` | 76 | 0.718 | 0.661 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `atac_tss_enrichment_score_mean` | 76 | -0.717 | -0.695 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `stage` | 76 | 0.688 | 0.654 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `mitochondrial_supertype_share` | 76 | 0.595 | 0.630 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 76 | -0.561 | -0.631 | raw+partial |
| primary | 1 | pls_0 | 24.24 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 76 | -0.514 | -0.512 | raw+partial |

### cohort_50

| analysis | rank | component | sd | variable | n | rho | rho \| covariates | flag |
|---|---:|---|---:|---|---:|---:|---:|---|
| primary | 1 | pls_0 | 23.64 | `mitochondrial_pseudobulk` | 63 | 0.785 | 0.816 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `fraction_mitochondrial_umis_mean` | 63 | 0.783 | 0.811 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `mitochondrial_mean` | 63 | 0.783 | 0.811 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `fraction_mitochondrial_umis_sd` | 63 | 0.772 | 0.798 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `mitochondrial_sd` | 63 | 0.772 | 0.798 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `stage` | 63 | 0.731 | 0.675 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `gex_sequencing_saturation_mean` | 63 | 0.725 | 0.725 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 63 | -0.692 | -0.700 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `atac_tss_enrichment_score_mean` | 63 | -0.621 | -0.653 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 63 | -0.532 | -0.540 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 63 | -0.520 | -0.716 | raw+partial |
| primary | 1 | pls_0 | 23.64 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 63 | -0.518 | -0.573 | raw+partial |

Descriptive, not an expectation (from `donor_technical.csv` and the export's `metadata.csv`): Spearman correlation of the donor's mitochondrial pseudobulk fraction with stage over the cohort donors: cohort_30 microglia 0.498, astrocyte 0.312; cohort_50 microglia 0.440, astrocyte 0.356.

## 7. Reading

**The control does not hold, and by the predeclared rule that puts the latent construction in question for any cell type.** The three expectations separate where the pipeline succeeds from where it does not. The stage signal is there to be found: by DESeq2, microglia show tens of exported RNA genes with a stage effect where astrocytes show one or two. The latent space built from the same blocks separates the stages only modestly better than the astrocytes' in the primary cohort, and not at all better in `cohort_50`. Its first and only component (rank 1) is most closely tied to the donor's mitochondrial UMI fraction; the ATAC library-quality metrics and the sequencing saturation sit close to stage, above it in `cohort_30` (section 6). Stage is strongly present on that axis (|ρ| 0.688), so the axis is not stage-blind. But it is not a stage axis in the sense E3 asked for, and the partial correlations given the six analysis covariates order the two the same way (0.654 with stage against 0.801 with `mitochondrial_pseudobulk`). This repeats the astrocyte finding that the CV-selected rank-1 axis is dominated by a technical, library-quality direction (`reports/motco-astrocyte-technical-axis-2026-10-19.md`), now in a cell type whose pseudobulks hold far more per-feature stage signal.

**What this says about the astrocyte null.** The control was built so that a holding result would rule out "the pipeline cannot see stage structure in this data" as the explanation of the astrocyte non-rejections. It did not hold, so that explanation is not ruled out. The failure does not show that the pipeline cannot resolve stage: E1 and E2 held in the primary cohort. It shows that, on these blocks, the stage-supervised latent space with CV-selected rank picks a single direction that technical variation shares with stage. The astrocyte reading is unchanged (its non-rejections stay uninformative under the binding rules); this result only removes the support the control would have given it.

**A caution on E3 (descriptive, not an expectation).** Among the microglia cohort donors the mitochondrial pseudobulk fraction itself correlates with stage (Spearman 0.498 in `cohort_30`, 0.440 in `cohort_50`; astrocytes 0.312 and 0.356). So the technical variable that beats stage is partly stage-associated, and the diagnostic cannot tell a preparation artefact from a biological correlate of pathology (a rising mitochondrial transcript share in diseased microglia would look the same). The design stated this risk for E3 before any run ("a strongly technical cell type could fail it even with real stage signal"); the expectation is not re-tuned, and the reading above is limited accordingly.

**E1's method difference is mostly a difference of test.** As predeclared, E1 reads DESeq2's union of the two Wald tests of `factor_additive` and limma's joint moderated F over the same two coefficients. On the single coefficient stage 2 vs 0 the methods are of the same order in microglia: 61 exported RNA genes by DESeq2 and 28 by limma in `cohort_30` (34 and 36 in `cohort_50`). The joint F spends two degrees of freedom on a stage 1 vs 0 contrast that carries almost no signal in either method (section 4) and finds 2 and 0. So E1 holds under both methods in `cohort_30` by the predeclared comparison, under limma by a margin of one gene, and is method dependent in `cohort_50`. Under the rules of the limma change this is method dependence, a result of neither method. No method or test is promoted, and E1 is not re-read with another test.

**The microglia sex × stage result** is a second case-study result read under the binding rules. Nothing is rejected in either cohort. The `angle` and `shape` non-rejections are uninformative, and at rank 1 `angle` is 0° or 180° by construction. The `delta` non-rejection is not evidence that the sexes' trajectories are the same. Both concordance cells are concordant in every cohort and omic: no exported feature shows a sex × stage interaction under either method, and the primary `delta` is not rejected. This cross-check agrees with the primary analysis; it is not a second test of its hypothesis.

**Next.** Nothing here changes an astrocyte output or reading. A follow-up, for example the latent-adjustment arm on the microglia blocks, or a control on a cell type with a different technical profile, would need its own proposal, with expectations written before its runs.

## 8. Limitations

- E1–E3 compare values; they test nothing and carry no uncertainty. Ties do not hold, and the limma E1 margin in `cohort_30` is a single gene.
- The two cell types' cohorts are not the same donors (sizes in section 2): 73 donors are in both `cohort_30`s and 60 in both `cohort_50`s. The comparison is between the cell types' own primary analyses, as predeclared, not a paired comparison.
- Microglia ATAC nuclei per donor are few (section 2), so the ATAC pseudobulks are noisier than the astrocytes', and the ATAC counts set the cohorts.
- The `Micro-PVM_*` population bundles microglia and perivascular macrophages and includes the disease-associated supertypes, whose share may change with stage. Part of the microglia stage signal may therefore be compositional; this change does not separate composition from expression.
- MOTCO's operating characteristics were measured on neither cell type's blocks. Rank 1 is below the ranks at which they were measured, and the adjusted analysis has no Type I check.
- `stage-de` ran with `--n-cpus 10`, as for the astrocytes: the tables do not depend on it, but warnings raised in pydeseq2's worker processes are not captured (recorded in the manifests).
- The two deviations (the spike over S3, the cluster steps on the login node) change where the code ran, not what it read or computed. Object sizes and ETags were checked against the manifest, the code was the same commit, and the astrocyte files hash unchanged.

## 9. Reproduction

```
# cluster (login node; the sbatch scripts with CELL_TYPE=microglia are the documented route)
sead-sexdiff covariates --cell-type microglia --omic both
sead-sexdiff pseudobulk --cell-type microglia --omic both
sead-sexdiff technical-covariates --cell-type microglia --omic both
# workstation, after scripts/sync_pseudobulk.sh --cell-type microglia (in the Nix shell)
sead-sexdiff export --cell-type microglia
sead-sexdiff analyze --cell-type microglia --n-jobs 1
sead-sexdiff axis-diagnostics --cell-type microglia
sead-sexdiff stage-de --cell-type microglia --omic both --n-cpus 10
sead-sexdiff stage-de-limma --cell-type microglia --omic both
sead-sexdiff positive-control
```
