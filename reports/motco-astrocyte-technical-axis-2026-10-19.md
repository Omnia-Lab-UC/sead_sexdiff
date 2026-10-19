# MOTCO astrocyte technical-axis diagnostic: sex × ADNC stage — 2026-10-19

> Measurement-space diagnostic of the covariate-adjusted comparison of female and male astrocyte trajectories across ADNC
> stages in the SEA-AD MTG (OpenSpec change `mtg-astrocyte-technical-axis-diagnostic`, specs `donor-technical-table` and
> `latent-axis-diagnostic`). It asks what the latent axis the primary analysis (`reports/motco-astrocyte-2026-10-19.md`), the
> forced-rank sweep and the chrX variant all measured on correlates with, and whether the technical quantities it tracks are
> donor-level properties of the preparation or vary with astrocyte state within donors. Every number below is pasted by script
> from the generated files: both cohorts' `diagnostics.md` (rendered from the four CSVs of `results/motco_diagnostics/<cohort>/`),
> the donor-level technical summary (`results/pseudobulk/donor_technical.csv`, `technical_manifest.json`) and the manifests,
> and was checked against those files before commit. The hand-written text carries only the reading, under the rules the spec
> binds (section 5). Nothing here changes the primary analysis, its statistics or their reading, and no preprocessing variant
> is chosen in this report.

**Runs.** `sead-sexdiff technical-covariates --omic both` on the cluster login node (data root `/home1/tgonza/sead_sexdiff/data`; a `nohup` run after
the sbatch job 922572 on node n11 produced no Python output in nine minutes and was cancelled with nothing cached), 2026-10-19,
about 14 min for the 84 RNA and 86 ATAC objects, at commit `3ae7553` with
the change's implementation uncommitted (`dirty: true`); then `scripts/sync_pseudobulk.sh` and
`sead-sexdiff axis-diagnostics --cohort all` on the NixOS workstation `nixos-home`, seconds per cohort, same commit (`dirty:
true`). Every cross-check passed: each donor's gene-set fractions from the per-nucleus `UMIs` rows equal the
fractions from its `rna_counts.h5ad` row, and every nucleus count equals the pseudobulk cache's and `donor_table.csv`'s; every
analysis manifest's export hashes equal the export files on disk. **Result.** In both cohorts the single CV-selected latent
component is, in the measurement space, a donor-level composition and library-quality axis: in the primary cohort it correlates
-0.84 with the astrocyte nuclei's mean ribosomal-protein UMI fraction (-0.82 given the six analysis covariates),
0.75 with their mitochondrial fraction, -0.77 and -0.77 with the ATAC library's TSS enrichment and fraction
of fragments in peaks and -0.73 with the RNA TMM factor, against 0.49 with stage and -0.23 with sex; the
sensitivity cohort gives the same picture (-0.84 with the ribosomal fraction). PLS components nest, so this axis is the first
component of every forced rank of the sweep and carries 98 % of the squared sex-difference geometry at rank 2,
65 % at rank 4 and 49 % at rank 12. The ribosomal and mitochondrial fractions are donor-level properties
(between-donor shares 0.23 and 0.42 of the per-nucleus variance against median within-donor supertype
shares of 0.04 and 0.03), so under the spec's rule both are residualization candidates for the next
change; the hemoglobin fraction is negligible (median 6.7e-07 of UMIs) and is not. These are findings about the
measurement space of the primary analysis, not a reinterpretation of its statistics: the primary analysis stays primary and
nothing was rejected in it.

## 1. Provenance

| item | value |
|---|---|
| MOTCO | 0.8.0, tag `v0.8.0`, commit `979657c3d0e6c9ad4d4bff0dc595a18df854ca61` (pinned in `uv.lock`); the analyses read were made with MOTCO 0.8.0 |
| entry points | `sead_sexdiff.data.run_technical_covariates` (`select_astrocyte_nuclei`, the pseudobulk's nucleus selection, then `technical_donor` per object) and `sead_sexdiff.data.run_axis_diagnostics` |
| sead_sexdiff | `3ae7553` (the change's implementation uncommitted at run time; both manifests record `dirty: true`) |
| machines | technical table: cluster login node (`/home1/tgonza/sead_sexdiff/data`); diagnostics: `nixos-home` (NixOS workstation), after `scripts/sync_pseudobulk.sh` |
| pseudobulk inputs | `results/pseudobulk/` (`manifest.json` sha256 `2e93cc0b5268…`, created 2026-10-16T18:46:13+00:00; `rna_counts.h5ad` sha256 `1253bed9e5f4…`), the pseudobulk the primary export recorded |
| covariate files | untouched: `donor_covariates.csv` sha256 `b14398d88cd8…` and `covariates_manifest.json` sha256 `b8bf9091b507…` before and after the technical run; the latter is the hash the primary export recorded |
| technical table | `results/pseudobulk/donor_technical.csv` sha256 `eac1b9a15977…`; `technical_manifest.json` sha256 `c2bdb3449964…`; 84 donors with a row; gene sets `ribosomal` (81 genes), `mitochondrial` (13), `hemoglobin` (3) |
| primary exports | `results/motco_input/<cohort>/`: manifest sha256 `18f2289e0d4c…` (cohort_30), `dd14a6e7d558…` (cohort_50); read, hashed, never written |
| primary analyses | `results/motco_analysis/<cohort>/`: CV-selected rank 1 in both; `result.json` sha256 `ca92876a1a61…` (cohort_30), `1e7eefe3e63e…` (cohort_50); read, checked against the export hashes, never rerun |
| forced-rank analyses | `results/motco_rank_sensitivity/<cohort>/rank_{02,03,04,06,09,12}/` (`reports/motco-astrocyte-rank-sensitivity-2026-10-19.md`), each checked against the export hashes |
| variables | 47: the 35 numeric donor columns of the technical table, `tmm_rna`, `log10_lib_rna`, `tmm_atac`, `log10_lib_atac` from the export manifest, the six analysis covariates (`age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`), `stage` (0, 1, 2) and `male` |
| statistics | Spearman correlations over the donors with the variable; partial correlation = Spearman correlation of the OLS residuals of both variables on an intercept and the six analysis covariates; block PCs from the SVD of each z-scored block (population SD); geometry shares from `stage_means.csv`; threshold |rho| ≥ 0.3 is a reporting flag, not an adoption rule; no p-values |
| software | python 3.12.15, numpy 2.3.5, scipy 1.15.3, pandas 2.3.3, anndata 0.12.7, h5py 3.15.1 |

## 2. The donor technical table

Per RNA donor object the astrocyte nuclei were selected by the pseudobulk's own code path (taxonomy mapping, `Reference`
donors excluded, ≥ 200 detected genes in the raw `UMIs` layer), and for each nucleus the UMI fraction in each predeclared gene
set was computed from the raw layer (all genes) beside the `obs` QC columns; the normalized RNA `X` was never read. The raw
layer reproduces SEA-AD's own per-nucleus `Number of UMIs` and `Genes detected` exactly (maximum absolute difference 0 in every
donor), and the computed mitochondrial fraction is the same quantity as SEA-AD's `Fraction mitochondrial UMIs`, which is why the
two appear with identical correlations in sections 3 and 4. 84 donors have a technical row (the two ATAC-only `Reference` donors
`H19.30.002` and `H19.30.004` are excluded, as in the pseudobulk); every donor also has an ATAC object. Astrocyte nuclei per donor:
median 667.5, range 76–2,266. 80 of the 84 donors' astrocytes span more than one
10x library (`load_name`: 1 library for 4 donors, 2 for 46, 3 for 29, 4 for 3, 5 for 2),
so the library-level `GEX_*` metrics are astrocyte-weighted means per donor, and 29 ATAC objects also carry two or three
distinct values of each `ATAC_*` metric (several libraries per donor there too; 145 donor × column cells flagged
non-constant in the manifest), so the ATAC metrics are astrocyte-weighted means as well. 28 donors have nuclei with no
`GEX_sequencing_saturation` value (dropped from that mean; every donor keeps a mean). The hemoglobin fraction has no variance
among the astrocytes of 6 donors, so its supertype share is undefined there.

Per set and donor the table holds the nucleus-weighted mean and SD of the per-nucleus fraction, the UMI-weighted pseudobulk
fraction (equal, to 1e-9, to the fraction from the donor's `rna_counts.h5ad` row) and the within-donor supertype share (η²: the
share of the per-nucleus fraction's variance explained by astrocyte supertype within that donor, over the supertypes present).
The manifest adds the between-donor share: the share of the pooled per-nucleus variance, over all nuclei of all donors, that lies
between donors (donor η²). Across the 84 donors:

| set | genes | pseudobulk fraction: median [IQR] (range) | nucleus-mean fraction: median | between-donor share | supertype share: median (q75, max; n defined) | donors with share < 0.1 | donors with share > between-donor share |
|---|---:|---|---:|---:|---|---:|---:|
| `ribosomal` | 81 | 0.00179 [0.00127, 0.00237] (0.000792–0.00398) | 0.00194 | 0.227 | 0.0395 (0.08, 0.271; 84) | 67 | 1 |
| `mitochondrial` | 13 | 0.000892 [0.000499, 0.00186] (8.62e-05–0.0127) | 0.00109 | 0.422 | 0.0311 (0.0424, 0.262; 84) | 81 | 0 |
| `hemoglobin` | 3 | 6.68e-07 [2.7e-07, 1.4e-06] (0–5.37e-05) | 9.06e-07 | 0.00285 | 0.00394 (0.0104, 0.0684; 78) | 78 | 48 |

| QC column (astrocyte-weighted mean per donor) | median | range | between-donor share |
|---|---:|---|---:|
| `fraction_mitochondrial_umis` | 0.001088 | 8.998e-05–0.01395 | 0.422 |
| `number_of_umis` | 1.12e+04 | 4127–1.838e+04 | 0.149 |
| `genes_detected` | 3856 | 2069–4761 | 0.158 |
| `doublet_score` | 0.02021 | 0.0085–0.04628 | 0.0898 |
| `gex_reads_mapped_confidently_to_intronic_regions` | 0.6258 | 0.144–0.6906 | 0.965 |
| `gex_fraction_of_transcriptomic_reads_in_cells` | 0.8781 | 0.4658–0.9356 | 0.783 |
| `gex_median_umi_counts_per_cell` | 1.625e+04 | 2656–2.731e+04 | 0.811 |
| `gex_sequencing_saturation` | 0.6932 | 0.5087–0.919 | 0.74 |

**Classification under the spec's rule** (a fraction whose between-donor share dominates and whose within-donor supertype share
is small in most donors is a preparation property and a residualization candidate; one whose within-donor supertype share exceeds
its between-donor share is cell state and not a candidate):

- `ribosomal`: between-donor share 0.227; supertype share median 0.039, below 0.1 in 67 of 84 donors
  and above the between-donor share in 1. A donor-level preparation property → **residualization candidate**.
- `mitochondrial`: between-donor share 0.422; supertype share median 0.031, below 0.1 in 81 of 84 donors and
  above the between-donor share in 0. A donor-level preparation property → **residualization candidate** (it is the same
  quantity as SEA-AD's `Fraction mitochondrial UMIs`).
- `hemoglobin`: median pseudobulk fraction 6.7e-07 (at most 5.4e-05) of UMIs, between-donor share 0.00285 and supertype share
  median 0.00394, both negligible (the 48 donors whose supertype share exceeds the between-donor share do so
  because both are near zero). There is no blood-contamination signal to residualize → **not a candidate**.
- The library-level `GEX_*` metrics have between-donor shares of 0.74–0.96 by construction (they vary only across
  libraries), and the ATAC library metrics are donor-level by construction; whether any of them enters a residualized variant is
  the next change's decision, written in its proposal before it runs.

## 3. Primary cohort: cohort_30 (n = 80)

Primary analysis `results/motco_analysis/cohort_30` (CV-selected rank 1; `result.json` sha256 `ca92876a1a61…`;
`delta` 14.941, p 0.3; `angle` 0°, p 1; `shape`
0.608832, p 0.333; nothing rejected at α = 0.05), read and checked against the export hashes, never
rerun. Forced-rank analyses read: rank_02 (rank 2), rank_03 (rank 3), rank_04 (rank 4), rank_06 (rank 6), rank_09 (rank 9), rank_12 (rank 12), each checked the same way.
47 variables; threshold |rho| ≥ 0.3 (a reporting flag); no p-values.

**What `pls_0` is, in the measurement space.** The single CV-selected component (score SD 34.99) correlates
-0.841 with the astrocyte nuclei's mean ribosomal-protein UMI fraction (`ribosomal_mean`; -0.821 given the
six analysis covariates) and -0.813 with the UMI-weighted pseudobulk ribosomal fraction, 0.751 with the mean
mitochondrial fraction (`mitochondrial_mean`, identical to SEA-AD's own `Fraction mitochondrial UMIs`; partial 0.593),
-0.767 and -0.767 with the two ATAC library metrics (TSS enrichment and the fraction of fragments in peaks,
partials -0.673 and -0.642), and -0.733 with the RNA TMM factor (partial -0.670), against
0.491 with stage (partial 0.406) and -0.229 with sex (`male`; partial -0.325).
20 of 47 variables are flagged at |rho| ≥ 0.3 on the raw or the partial correlation; the six analysis covariates
leave the technical correlations largely in place (the partial correlations in the table), so the covariate adjustment of the primary
analysis does not remove this axis from the latent space, which MOTCO fits without the covariates (covariate-blind basis).

**Geometry.** The sexes' net displacements along `pls_0` are 32.63 (Female) and 17.69 (Male) latent units; at rank 1 every share
is 1 by construction. `pls_0` is identical (`np.allclose`) across the primary and all six forced ranks, so the forced-rank
sweep inherited this axis as its first and largest component: the share of the squared per-stage sex-difference geometry along `pls_0`
is 0.976 at rank 2, 0.650 at rank 4 and 0.494 at rank 12 (stage-centred 0.923, 0.506, 0.280); of the
female net displacement 0.671, 0.391, 0.375 and of the male 0.283, 0.202, 0.155. `delta` in the primary analysis
measured the difference between those two displacements along this axis.

**Block axes.** The first principal component of the z-scored RNA block holds 0.277 of its variance and correlates
0.803 with `ribosomal_mean`, 0.772 with `ribosomal_pseudobulk` and 0.749 with `tmm_rna`
(19 variables flagged); the ATAC block's first component holds 0.145 of its variance and correlates
0.637 with `log10_n_atac`, 0.608 with `log10_lib_atac` and -0.561 with `tmm_rna`
(17 flagged). The composition axis is therefore the dominant axis of the exported RNA block itself, not an artefact of
the PLS fit.

**Feature selection.** Of 20,387 abundance-passing genes, the 5,000 kept by the top-N variance rule have a median mean
log2-CPM of 0.77 (5–95 % quantiles -0.72 to 4.97), 29.3 % of them below 0 and 54.7 % below 1,
and hold 56.7 % of the abundance-passing variance (cut at 0.646; dropped median 0.212);
62.3 % (3,116) of the kept genes have |rho| ≥ 0.3 with `pls_0`. Of 128,916 abundance-passing
peaks, the 5,000 kept have a median mean of 0.37, 84.2 % below 1, and hold 7.7 % of the variance (cut
1.173 against a dropped median of 0.732: a near-flat cut); 40.2 % (2,009) correlate with `pls_0`
at that level.

Findings recorded by the step: (1) variables with missing donors (pairwise-complete correlations): {'hemoglobin_supertype_share': 5} (2) pls_0 is identical (np.allclose) across the primary and every forced rank ['rank_02', 'rank_03', 'rank_04', 'rank_06', 'rank_09', 'rank_12']: PLS components nest, so every forced rank keeps the primary's first axis (3) primary pls_0: 20 of 47 variables flagged at |rho| ≥ 0.3 (raw or partial): ribosomal_mean (-0.84; partial -0.82), ribosomal_pseudobulk (-0.81; partial -0.79), atac_fraction_of_high_quality_fragments_overlapping_peaks_mean (-0.77; partial -0.64), atac_tss_enrichment_score_mean (-0.77; partial -0.67), fraction_mitochondrial_umis_mean (0.75; partial 0.59), mitochondrial_mean (0.75; partial 0.59), tmm_rna (-0.73; partial -0.67), mitochondrial_pseudobulk (0.73; partial 0.56), ribosomal_sd (-0.70; partial -0.68), gex_sequencing_saturation_mean (0.67; partial 0.52), fraction_mitochondrial_umis_sd (0.65; partial 0.51), mitochondrial_sd (0.65; partial 0.51), gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.61; partial -0.41), atac_fraction_of_high_quality_fragments_in_cells_mean (-0.55; partial -0.37), gex_reads_mapped_confidently_to_intronic_regions_mean (-0.54; partial -0.52), stage (0.49; partial 0.41), tmm_atac (0.49; partial 0.50), gex_median_umi_counts_per_cell_mean (-0.47; partial -0.38), rin (-0.31; partial n/a), male (-0.23; partial -0.33)

Descriptive measurement-space diagnostics of the covariate-adjusted trajectory analysis: Spearman correlations of every latent component of the primary analysis (CV-selected rank 1) and of the forced-rank analyses (2, 3, 4, 6, 9, 12) with 47 predeclared variables (technical table, export TMM factors and library sizes, the six analysis covariates, stage and sex), the partial correlations given the six analysis covariates (Spearman correlation of the OLS residuals of both variables on an intercept and the covariates), the same for the first 3 principal components of each z-scored block, the share of the trajectory geometry along each component, and the feature-selection accounting. Rows with |rho| ≥ 0.3 are flagged and listed first; the flag is a reporting threshold, not an adoption rule. No p-value is attached to any correlation, the tables support no inferential claim, and nothing here chooses a preprocessing variant: an axis correlated with a technical variable is a finding about the measurement space of the primary analysis, not a reinterpretation of its statistics, and the primary analysis stays primary with its non-rejection reading rules.

### Flagged correlations (pasted from `diagnostics.md`)

| analysis | rank | component | sd | variable | n | rho | rho \| covariates | flag |
|---|---:|---|---:|---|---:|---:|---:|---|
| primary | 1 | pls_0 | 34.99 | `ribosomal_mean` | 80 | -0.841 | -0.821 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `ribosomal_pseudobulk` | 80 | -0.813 | -0.794 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | -0.767 | -0.642 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `atac_tss_enrichment_score_mean` | 80 | -0.767 | -0.673 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `fraction_mitochondrial_umis_mean` | 80 | 0.751 | 0.593 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `mitochondrial_mean` | 80 | 0.751 | 0.593 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `tmm_rna` | 80 | -0.733 | -0.670 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `mitochondrial_pseudobulk` | 80 | 0.729 | 0.558 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `ribosomal_sd` | 80 | -0.703 | -0.680 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `gex_sequencing_saturation_mean` | 80 | 0.672 | 0.522 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `fraction_mitochondrial_umis_sd` | 80 | 0.648 | 0.505 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `mitochondrial_sd` | 80 | 0.648 | 0.505 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | -0.607 | -0.412 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | -0.552 | -0.369 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | -0.540 | -0.520 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `stage` | 80 | 0.491 | 0.406 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `tmm_atac` | 80 | 0.486 | 0.505 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `gex_median_umi_counts_per_cell_mean` | 80 | -0.470 | -0.376 | raw+partial |
| primary | 1 | pls_0 | 34.99 | `rin` | 80 | -0.307 | n/a | raw |
| primary | 1 | pls_0 | 34.99 | `male` | 80 | -0.229 | -0.325 | partial |
| rank_02 | 2 | pls_0 | 34.99 | `ribosomal_mean` | 80 | -0.841 | -0.821 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `ribosomal_pseudobulk` | 80 | -0.813 | -0.794 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | -0.767 | -0.642 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `atac_tss_enrichment_score_mean` | 80 | -0.767 | -0.673 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `fraction_mitochondrial_umis_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `mitochondrial_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `tmm_rna` | 80 | -0.733 | -0.670 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `mitochondrial_pseudobulk` | 80 | 0.729 | 0.558 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `ribosomal_sd` | 80 | -0.703 | -0.680 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `gex_sequencing_saturation_mean` | 80 | 0.672 | 0.522 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `fraction_mitochondrial_umis_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `mitochondrial_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | -0.607 | -0.412 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | -0.552 | -0.369 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | -0.540 | -0.520 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `stage` | 80 | 0.491 | 0.406 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `tmm_atac` | 80 | 0.486 | 0.505 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `gex_median_umi_counts_per_cell_mean` | 80 | -0.470 | -0.376 | raw+partial |
| rank_02 | 2 | pls_0 | 34.99 | `rin` | 80 | -0.307 | n/a | raw |
| rank_02 | 2 | pls_0 | 34.99 | `male` | 80 | -0.229 | -0.325 | partial |
| rank_02 | 2 | pls_1 | 17.09 | `stage` | 80 | -0.616 | -0.742 | raw+partial |
| rank_02 | 2 | pls_1 | 17.09 | `log10_lib_atac` | 80 | 0.524 | 0.287 | raw |
| rank_02 | 2 | pls_1 | 17.09 | `log10_n_atac` | 80 | 0.440 | n/a | raw |
| rank_03 | 3 | pls_0 | 34.99 | `ribosomal_mean` | 80 | -0.841 | -0.821 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `ribosomal_pseudobulk` | 80 | -0.813 | -0.794 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | -0.767 | -0.642 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `atac_tss_enrichment_score_mean` | 80 | -0.767 | -0.673 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `fraction_mitochondrial_umis_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `mitochondrial_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `tmm_rna` | 80 | -0.733 | -0.670 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `mitochondrial_pseudobulk` | 80 | 0.729 | 0.558 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `ribosomal_sd` | 80 | -0.703 | -0.680 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `gex_sequencing_saturation_mean` | 80 | 0.672 | 0.522 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `fraction_mitochondrial_umis_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `mitochondrial_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | -0.607 | -0.412 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | -0.552 | -0.369 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | -0.540 | -0.520 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `stage` | 80 | 0.491 | 0.406 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `tmm_atac` | 80 | 0.486 | 0.505 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `gex_median_umi_counts_per_cell_mean` | 80 | -0.470 | -0.376 | raw+partial |
| rank_03 | 3 | pls_0 | 34.99 | `rin` | 80 | -0.307 | n/a | raw |
| rank_03 | 3 | pls_0 | 34.99 | `male` | 80 | -0.229 | -0.325 | partial |
| rank_03 | 3 | pls_1 | 17.09 | `stage` | 80 | -0.616 | -0.742 | raw+partial |
| rank_03 | 3 | pls_1 | 17.09 | `log10_lib_atac` | 80 | 0.524 | 0.287 | raw |
| rank_03 | 3 | pls_1 | 17.09 | `log10_n_atac` | 80 | 0.440 | n/a | raw |
| rank_03 | 3 | pls_2 | 20.1 | `stage` | 80 | 0.480 | 0.535 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `ribosomal_mean` | 80 | -0.841 | -0.821 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `ribosomal_pseudobulk` | 80 | -0.813 | -0.794 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | -0.767 | -0.642 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `atac_tss_enrichment_score_mean` | 80 | -0.767 | -0.673 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `fraction_mitochondrial_umis_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `mitochondrial_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `tmm_rna` | 80 | -0.733 | -0.670 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `mitochondrial_pseudobulk` | 80 | 0.729 | 0.558 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `ribosomal_sd` | 80 | -0.703 | -0.680 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `gex_sequencing_saturation_mean` | 80 | 0.672 | 0.522 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `fraction_mitochondrial_umis_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `mitochondrial_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | -0.607 | -0.412 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | -0.552 | -0.369 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | -0.540 | -0.520 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `stage` | 80 | 0.491 | 0.406 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `tmm_atac` | 80 | 0.486 | 0.505 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `gex_median_umi_counts_per_cell_mean` | 80 | -0.470 | -0.376 | raw+partial |
| rank_04 | 4 | pls_0 | 34.99 | `rin` | 80 | -0.307 | n/a | raw |
| rank_04 | 4 | pls_0 | 34.99 | `male` | 80 | -0.229 | -0.325 | partial |
| rank_04 | 4 | pls_1 | 17.09 | `stage` | 80 | -0.616 | -0.742 | raw+partial |
| rank_04 | 4 | pls_1 | 17.09 | `log10_lib_atac` | 80 | 0.524 | 0.287 | raw |
| rank_04 | 4 | pls_1 | 17.09 | `log10_n_atac` | 80 | 0.440 | n/a | raw |
| rank_04 | 4 | pls_2 | 20.1 | `stage` | 80 | 0.480 | 0.535 | raw+partial |
| rank_04 | 4 | pls_3 | 14.47 | `log10_n_atac` | 80 | -0.442 | n/a | raw |
| rank_04 | 4 | pls_3 | 14.47 | `log10_lib_atac` | 80 | -0.424 | -0.077 | raw |
| rank_04 | 4 | pls_3 | 14.47 | `stage` | 80 | -0.236 | -0.303 | partial |
| rank_06 | 6 | pls_0 | 34.99 | `ribosomal_mean` | 80 | -0.841 | -0.821 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `ribosomal_pseudobulk` | 80 | -0.813 | -0.794 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | -0.767 | -0.642 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `atac_tss_enrichment_score_mean` | 80 | -0.767 | -0.673 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `fraction_mitochondrial_umis_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `mitochondrial_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `tmm_rna` | 80 | -0.733 | -0.670 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `mitochondrial_pseudobulk` | 80 | 0.729 | 0.558 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `ribosomal_sd` | 80 | -0.703 | -0.680 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `gex_sequencing_saturation_mean` | 80 | 0.672 | 0.522 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `fraction_mitochondrial_umis_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `mitochondrial_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | -0.607 | -0.412 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | -0.552 | -0.369 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | -0.540 | -0.520 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `stage` | 80 | 0.491 | 0.406 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `tmm_atac` | 80 | 0.486 | 0.505 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `gex_median_umi_counts_per_cell_mean` | 80 | -0.470 | -0.376 | raw+partial |
| rank_06 | 6 | pls_0 | 34.99 | `rin` | 80 | -0.307 | n/a | raw |
| rank_06 | 6 | pls_0 | 34.99 | `male` | 80 | -0.229 | -0.325 | partial |
| rank_06 | 6 | pls_1 | 17.09 | `stage` | 80 | -0.616 | -0.742 | raw+partial |
| rank_06 | 6 | pls_1 | 17.09 | `log10_lib_atac` | 80 | 0.524 | 0.287 | raw |
| rank_06 | 6 | pls_1 | 17.09 | `log10_n_atac` | 80 | 0.440 | n/a | raw |
| rank_06 | 6 | pls_2 | 20.1 | `stage` | 80 | 0.480 | 0.535 | raw+partial |
| rank_06 | 6 | pls_3 | 14.47 | `log10_n_atac` | 80 | -0.442 | n/a | raw |
| rank_06 | 6 | pls_3 | 14.47 | `log10_lib_atac` | 80 | -0.424 | -0.077 | raw |
| rank_06 | 6 | pls_3 | 14.47 | `stage` | 80 | -0.236 | -0.303 | partial |
| rank_06 | 6 | pls_4 | 15.82 | `log10_n_atac` | 80 | -0.309 | n/a | raw |
| rank_06 | 6 | pls_5 | 12.99 | `log10_n_rna` | 80 | -0.351 | n/a | raw |
| rank_06 | 6 | pls_5 | 12.99 | `n_nuclei_rna` | 80 | -0.351 | -0.003 | raw |
| rank_06 | 6 | pls_5 | 12.99 | `doublet_score_mean` | 80 | 0.347 | 0.197 | raw |
| rank_06 | 6 | pls_5 | 12.99 | `log10_lib_rna` | 80 | -0.338 | -0.020 | raw |
| rank_09 | 9 | pls_0 | 34.99 | `ribosomal_mean` | 80 | -0.841 | -0.821 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `ribosomal_pseudobulk` | 80 | -0.813 | -0.794 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | -0.767 | -0.642 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `atac_tss_enrichment_score_mean` | 80 | -0.767 | -0.673 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `fraction_mitochondrial_umis_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `mitochondrial_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `tmm_rna` | 80 | -0.733 | -0.670 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `mitochondrial_pseudobulk` | 80 | 0.729 | 0.558 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `ribosomal_sd` | 80 | -0.703 | -0.680 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `gex_sequencing_saturation_mean` | 80 | 0.672 | 0.522 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `fraction_mitochondrial_umis_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `mitochondrial_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | -0.607 | -0.412 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | -0.552 | -0.369 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | -0.540 | -0.520 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `stage` | 80 | 0.491 | 0.406 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `tmm_atac` | 80 | 0.486 | 0.505 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `gex_median_umi_counts_per_cell_mean` | 80 | -0.470 | -0.376 | raw+partial |
| rank_09 | 9 | pls_0 | 34.99 | `rin` | 80 | -0.307 | n/a | raw |
| rank_09 | 9 | pls_0 | 34.99 | `male` | 80 | -0.229 | -0.325 | partial |
| rank_09 | 9 | pls_1 | 17.09 | `stage` | 80 | -0.616 | -0.742 | raw+partial |
| rank_09 | 9 | pls_1 | 17.09 | `log10_lib_atac` | 80 | 0.524 | 0.287 | raw |
| rank_09 | 9 | pls_1 | 17.09 | `log10_n_atac` | 80 | 0.440 | n/a | raw |
| rank_09 | 9 | pls_2 | 20.1 | `stage` | 80 | 0.480 | 0.535 | raw+partial |
| rank_09 | 9 | pls_3 | 14.47 | `log10_n_atac` | 80 | -0.442 | n/a | raw |
| rank_09 | 9 | pls_3 | 14.47 | `log10_lib_atac` | 80 | -0.424 | -0.077 | raw |
| rank_09 | 9 | pls_3 | 14.47 | `stage` | 80 | -0.236 | -0.303 | partial |
| rank_09 | 9 | pls_4 | 15.82 | `log10_n_atac` | 80 | -0.309 | n/a | raw |
| rank_09 | 9 | pls_5 | 12.99 | `log10_n_rna` | 80 | -0.351 | n/a | raw |
| rank_09 | 9 | pls_5 | 12.99 | `n_nuclei_rna` | 80 | -0.351 | -0.003 | raw |
| rank_09 | 9 | pls_5 | 12.99 | `doublet_score_mean` | 80 | 0.347 | 0.197 | raw |
| rank_09 | 9 | pls_5 | 12.99 | `log10_lib_rna` | 80 | -0.338 | -0.020 | raw |
| rank_09 | 9 | pls_8 | 12.94 | `n_nuclei_rna` | 80 | 0.037 | -0.324 | partial |
| rank_12 | 12 | pls_0 | 34.99 | `ribosomal_mean` | 80 | -0.841 | -0.821 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `ribosomal_pseudobulk` | 80 | -0.813 | -0.794 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | -0.767 | -0.642 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `atac_tss_enrichment_score_mean` | 80 | -0.767 | -0.673 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `fraction_mitochondrial_umis_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `mitochondrial_mean` | 80 | 0.751 | 0.593 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `tmm_rna` | 80 | -0.733 | -0.670 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `mitochondrial_pseudobulk` | 80 | 0.729 | 0.558 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `ribosomal_sd` | 80 | -0.703 | -0.680 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `gex_sequencing_saturation_mean` | 80 | 0.672 | 0.522 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `fraction_mitochondrial_umis_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `mitochondrial_sd` | 80 | 0.648 | 0.505 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | -0.607 | -0.412 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | -0.552 | -0.369 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | -0.540 | -0.520 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `stage` | 80 | 0.491 | 0.406 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `tmm_atac` | 80 | 0.486 | 0.505 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `gex_median_umi_counts_per_cell_mean` | 80 | -0.470 | -0.376 | raw+partial |
| rank_12 | 12 | pls_0 | 34.99 | `rin` | 80 | -0.307 | n/a | raw |
| rank_12 | 12 | pls_0 | 34.99 | `male` | 80 | -0.229 | -0.325 | partial |
| rank_12 | 12 | pls_1 | 17.09 | `stage` | 80 | -0.616 | -0.742 | raw+partial |
| rank_12 | 12 | pls_1 | 17.09 | `log10_lib_atac` | 80 | 0.524 | 0.287 | raw |
| rank_12 | 12 | pls_1 | 17.09 | `log10_n_atac` | 80 | 0.440 | n/a | raw |
| rank_12 | 12 | pls_10 | 12.87 | `atac_median_high_quality_fragments_per_cell_mean` | 80 | 0.315 | 0.140 | raw |
| rank_12 | 12 | pls_2 | 20.1 | `stage` | 80 | 0.480 | 0.535 | raw+partial |
| rank_12 | 12 | pls_3 | 14.47 | `log10_n_atac` | 80 | -0.442 | n/a | raw |
| rank_12 | 12 | pls_3 | 14.47 | `log10_lib_atac` | 80 | -0.424 | -0.077 | raw |
| rank_12 | 12 | pls_3 | 14.47 | `stage` | 80 | -0.236 | -0.303 | partial |
| rank_12 | 12 | pls_4 | 15.82 | `log10_n_atac` | 80 | -0.309 | n/a | raw |
| rank_12 | 12 | pls_5 | 12.99 | `log10_n_rna` | 80 | -0.351 | n/a | raw |
| rank_12 | 12 | pls_5 | 12.99 | `n_nuclei_rna` | 80 | -0.351 | -0.003 | raw |
| rank_12 | 12 | pls_5 | 12.99 | `doublet_score_mean` | 80 | 0.347 | 0.197 | raw |
| rank_12 | 12 | pls_5 | 12.99 | `log10_lib_rna` | 80 | -0.338 | -0.020 | raw |
| rank_12 | 12 | pls_8 | 12.94 | `n_nuclei_rna` | 80 | 0.037 | -0.324 | partial |
| rank_12 | 12 | pls_9 | 11.18 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | -0.170 | -0.328 | partial |

Block principal components:

| block | component | variance share | variable | n | rho | rho \| covariates | flag |
|---|---|---:|---|---:|---:|---:|---|
| atac | pc1 | 0.145 | `log10_n_atac` | 80 | 0.637 | n/a | raw |
| atac | pc1 | 0.145 | `log10_lib_atac` | 80 | 0.608 | 0.162 | raw |
| atac | pc1 | 0.145 | `tmm_rna` | 80 | -0.561 | -0.443 | raw+partial |
| atac | pc1 | 0.145 | `gex_sequencing_saturation_mean` | 80 | 0.454 | 0.310 | raw+partial |
| atac | pc1 | 0.145 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | -0.449 | -0.323 | raw+partial |
| atac | pc1 | 0.145 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 80 | 0.441 | 0.154 | raw |
| atac | pc1 | 0.145 | `fraction_mitochondrial_umis_sd` | 80 | 0.426 | 0.217 | raw |
| atac | pc1 | 0.145 | `mitochondrial_sd` | 80 | 0.426 | 0.217 | raw |
| atac | pc1 | 0.145 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | -0.423 | -0.278 | raw |
| atac | pc1 | 0.145 | `fraction_mitochondrial_umis_mean` | 80 | 0.404 | 0.214 | raw |
| atac | pc1 | 0.145 | `mitochondrial_mean` | 80 | 0.404 | 0.214 | raw |
| atac | pc1 | 0.145 | `mitochondrial_pseudobulk` | 80 | 0.390 | 0.225 | raw |
| atac | pc1 | 0.145 | `atac_tss_enrichment_score_mean` | 80 | -0.377 | -0.400 | raw+partial |
| atac | pc1 | 0.145 | `ribosomal_mean` | 80 | -0.370 | -0.308 | raw+partial |
| atac | pc1 | 0.145 | `ribosomal_pseudobulk` | 80 | -0.354 | -0.299 | raw |
| atac | pc1 | 0.145 | `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 80 | 0.349 | 0.103 | raw |
| atac | pc1 | 0.145 | `n_libraries_rna` | 80 | 0.300 | 0.100 | raw |
| atac | pc2 | 0.110 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | 0.613 | 0.512 | raw+partial |
| atac | pc2 | 0.110 | `atac_tss_enrichment_score_mean` | 80 | 0.578 | 0.498 | raw+partial |
| atac | pc2 | 0.110 | `ribosomal_mean` | 80 | 0.571 | 0.449 | raw+partial |
| atac | pc2 | 0.110 | `ribosomal_pseudobulk` | 80 | 0.555 | 0.450 | raw+partial |
| atac | pc2 | 0.110 | `mitochondrial_pseudobulk` | 80 | -0.548 | -0.452 | raw+partial |
| atac | pc2 | 0.110 | `fraction_mitochondrial_umis_mean` | 80 | -0.547 | -0.466 | raw+partial |
| atac | pc2 | 0.110 | `mitochondrial_mean` | 80 | -0.547 | -0.466 | raw+partial |
| atac | pc2 | 0.110 | `fraction_mitochondrial_umis_sd` | 80 | -0.503 | -0.450 | raw+partial |
| atac | pc2 | 0.110 | `mitochondrial_sd` | 80 | -0.503 | -0.450 | raw+partial |
| atac | pc2 | 0.110 | `gex_sequencing_saturation_mean` | 80 | -0.461 | -0.285 | raw |
| atac | pc2 | 0.110 | `ribosomal_sd` | 80 | 0.447 | 0.328 | raw+partial |
| atac | pc2 | 0.110 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | 0.438 | 0.355 | raw+partial |
| atac | pc2 | 0.110 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | 0.416 | 0.351 | raw+partial |
| atac | pc2 | 0.110 | `tmm_rna` | 80 | 0.400 | 0.233 | raw |
| atac | pc2 | 0.110 | `gex_median_umi_counts_per_cell_mean` | 80 | 0.392 | 0.312 | raw+partial |
| atac | pc2 | 0.110 | `tmm_atac` | 80 | -0.339 | -0.403 | raw+partial |
| atac | pc2 | 0.110 | `log10_n_atac` | 80 | -0.335 | n/a | raw |
| atac | pc2 | 0.110 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 80 | -0.321 | -0.106 | raw |
| atac | pc2 | 0.110 | `atac_percent_duplicates_mean` | 80 | 0.314 | 0.362 | raw+partial |
| atac | pc2 | 0.110 | `log10_lib_rna` | 80 | -0.302 | -0.054 | raw |
| atac | pc2 | 0.110 | `doublet_score_sd` | 80 | -0.302 | -0.247 | raw |
| atac | pc2 | 0.110 | `male` | 80 | 0.197 | 0.310 | partial |
| atac | pc3 | 0.036 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | 0.445 | 0.275 | raw |
| atac | pc3 | 0.036 | `atac_tss_enrichment_score_mean` | 80 | 0.399 | 0.369 | raw+partial |
| atac | pc3 | 0.036 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | 0.398 | 0.330 | raw+partial |
| atac | pc3 | 0.036 | `tmm_atac` | 80 | -0.366 | -0.208 | raw |
| atac | pc3 | 0.036 | `log10_lib_atac` | 80 | 0.350 | 0.267 | raw |
| atac | pc3 | 0.036 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | 0.343 | 0.366 | raw+partial |
| atac | pc3 | 0.036 | `rin` | 80 | 0.301 | n/a | raw |
| atac | pc3 | 0.036 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | 0.300 | 0.119 | raw |
| rna | pc1 | 0.277 | `ribosomal_mean` | 80 | 0.803 | 0.785 | raw+partial |
| rna | pc1 | 0.277 | `ribosomal_pseudobulk` | 80 | 0.772 | 0.757 | raw+partial |
| rna | pc1 | 0.277 | `tmm_rna` | 80 | 0.749 | 0.690 | raw+partial |
| rna | pc1 | 0.277 | `fraction_mitochondrial_umis_mean` | 80 | -0.734 | -0.552 | raw+partial |
| rna | pc1 | 0.277 | `mitochondrial_mean` | 80 | -0.734 | -0.552 | raw+partial |
| rna | pc1 | 0.277 | `atac_tss_enrichment_score_mean` | 80 | 0.722 | 0.644 | raw+partial |
| rna | pc1 | 0.277 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | 0.708 | 0.583 | raw+partial |
| rna | pc1 | 0.277 | `mitochondrial_pseudobulk` | 80 | -0.706 | -0.514 | raw+partial |
| rna | pc1 | 0.277 | `ribosomal_sd` | 80 | 0.691 | 0.661 | raw+partial |
| rna | pc1 | 0.277 | `fraction_mitochondrial_umis_sd` | 80 | -0.651 | -0.488 | raw+partial |
| rna | pc1 | 0.277 | `mitochondrial_sd` | 80 | -0.651 | -0.488 | raw+partial |
| rna | pc1 | 0.277 | `gex_sequencing_saturation_mean` | 80 | -0.641 | -0.489 | raw+partial |
| rna | pc1 | 0.277 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | 0.593 | 0.396 | raw+partial |
| rna | pc1 | 0.277 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | 0.547 | 0.339 | raw+partial |
| rna | pc1 | 0.277 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | 0.534 | 0.513 | raw+partial |
| rna | pc1 | 0.277 | `tmm_atac` | 80 | -0.480 | -0.494 | raw+partial |
| rna | pc1 | 0.277 | `gex_median_umi_counts_per_cell_mean` | 80 | 0.457 | 0.365 | raw+partial |
| rna | pc1 | 0.277 | `doublet_score_sd` | 80 | -0.327 | -0.247 | raw |
| rna | pc1 | 0.277 | `stage` | 80 | -0.316 | -0.237 | raw |
| rna | pc2 | 0.080 | `ribosomal_supertype_share` | 80 | -0.383 | -0.249 | raw |
| rna | pc2 | 0.080 | `ribosomal_sd` | 80 | -0.365 | -0.221 | raw |
| rna | pc2 | 0.080 | `ribosomal_pseudobulk` | 80 | -0.359 | -0.243 | raw |
| rna | pc2 | 0.080 | `ribosomal_mean` | 80 | -0.357 | -0.225 | raw |
| rna | pc3 | 0.059 | `log10_n_rna` | 80 | -0.353 | n/a | raw |
| rna | pc3 | 0.059 | `n_nuclei_rna` | 80 | -0.353 | 0.167 | raw |

### Primary analysis (`results/motco_analysis/cohort_30`, cv rank 1) (pasted from `diagnostics.md`)

Component score SD: pls_0 34.99. Columns per component: Spearman rho, then rho given the six analysis covariates (n/a for a covariate itself). Variables sorted by name; `n` is the number of donors with the variable.

| variable | source | n | pls_0 rho | pls_0 rho \| cov |
|---|---|---:|---:|---:|
| `age_at_death` | covariates.csv | 80 | -0.121 | n/a |
| `atac_fraction_of_high_quality_fragments_in_cells_mean` | donor_technical.csv | 80 | -0.552* | -0.369* |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | donor_technical.csv | 80 | -0.767* | -0.642* |
| `atac_median_high_quality_fragments_per_cell_mean` | donor_technical.csv | 80 | -0.047 | 0.211 |
| `atac_percent_duplicates_mean` | donor_technical.csv | 80 | -0.234 | -0.264 |
| `atac_tss_enrichment_score_mean` | donor_technical.csv | 80 | -0.767* | -0.673* |
| `doublet_score_mean` | donor_technical.csv | 80 | 0.041 | 0.146 |
| `doublet_score_sd` | donor_technical.csv | 80 | 0.289 | 0.217 |
| `fraction_mitochondrial_umis_mean` | donor_technical.csv | 80 | 0.751* | 0.593* |
| `fraction_mitochondrial_umis_sd` | donor_technical.csv | 80 | 0.648* | 0.505* |
| `genes_detected_mean` | donor_technical.csv | 80 | -0.069 | -0.030 |
| `genes_detected_sd` | donor_technical.csv | 80 | 0.065 | 0.056 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 80 | -0.607* | -0.412* |
| `gex_fraction_of_transcriptomic_reads_in_cells_sd` | donor_technical.csv | 80 | 0.223 | 0.113 |
| `gex_median_umi_counts_per_cell_mean` | donor_technical.csv | 80 | -0.470* | -0.376* |
| `gex_median_umi_counts_per_cell_sd` | donor_technical.csv | 80 | -0.060 | -0.015 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 80 | -0.540* | -0.520* |
| `gex_reads_mapped_confidently_to_intronic_regions_sd` | donor_technical.csv | 80 | 0.288 | 0.086 |
| `gex_sequencing_saturation_mean` | donor_technical.csv | 80 | 0.672* | 0.522* |
| `gex_sequencing_saturation_sd` | donor_technical.csv | 80 | -0.070 | -0.088 |
| `hemoglobin_mean` | donor_technical.csv | 80 | -0.136 | -0.116 |
| `hemoglobin_pseudobulk` | donor_technical.csv | 80 | -0.190 | -0.084 |
| `hemoglobin_sd` | donor_technical.csv | 80 | -0.084 | -0.007 |
| `hemoglobin_supertype_share` | donor_technical.csv | 75 | -0.243 | -0.095 |
| `log10_lib_atac` | manifest.json (export, tmm) | 80 | 0.208 | -0.001 |
| `log10_lib_rna` | manifest.json (export, tmm) | 80 | 0.265 | 0.076 |
| `log10_n_atac` | covariates.csv | 80 | 0.282 | n/a |
| `log10_n_rna` | covariates.csv | 80 | 0.264 | n/a |
| `male` | metadata.csv | 80 | -0.229 | -0.325* |
| `mitochondrial_mean` | donor_technical.csv | 80 | 0.751* | 0.593* |
| `mitochondrial_pseudobulk` | donor_technical.csv | 80 | 0.729* | 0.558* |
| `mitochondrial_sd` | donor_technical.csv | 80 | 0.648* | 0.505* |
| `mitochondrial_supertype_share` | donor_technical.csv | 80 | 0.170 | 0.230 |
| `n_libraries_rna` | donor_technical.csv | 80 | 0.112 | 0.119 |
| `n_nuclei_rna` | donor_technical.csv | 80 | 0.264 | 0.167 |
| `number_of_umis_mean` | donor_technical.csv | 80 | 0.079 | 0.115 |
| `number_of_umis_sd` | donor_technical.csv | 80 | 0.128 | 0.191 |
| `pmi` | covariates.csv | 80 | 0.054 | n/a |
| `ribosomal_mean` | donor_technical.csv | 80 | -0.841* | -0.821* |
| `ribosomal_pseudobulk` | donor_technical.csv | 80 | -0.813* | -0.794* |
| `ribosomal_sd` | donor_technical.csv | 80 | -0.703* | -0.680* |
| `ribosomal_supertype_share` | donor_technical.csv | 80 | -0.077 | 0.187 |
| `rin` | covariates.csv | 80 | -0.307* | n/a |
| `stage` | metadata.csv | 80 | 0.491* | 0.406* |
| `tmm_atac` | manifest.json (export, tmm) | 80 | 0.486* | 0.505* |
| `tmm_rna` | manifest.json (export, tmm) | 80 | -0.733* | -0.670* |
| `years_education` | covariates.csv | 80 | -0.051 | n/a |

`*` marks |rho| ≥ 0.3.

### Forced-rank analyses (pasted from `diagnostics.md`)

Per component: score SD, the number of variables flagged on the raw correlation, the variable with the largest |rho| and that rho, the share of each sex's net displacement and of the per-stage sex-difference vectors along the component. The full rows are in `axis_correlations.csv` and `geometry_shares.csv`.

| analysis | rank | component | sd | flagged | top variable | rho | share net Female | share net Male | share sex diff | share sex diff centred |
|---|---:|---|---:|---:|---|---:|---:|---:|---:|---:|
| rank_02 | 2 | pls_0 | 34.99 | 19 | `ribosomal_mean` | -0.841 | 0.671 | 0.283 | 0.976 | 0.923 |
| rank_02 | 2 | pls_1 | 17.09 | 3 | `stage` | -0.616 | 0.329 | 0.717 | 0.024 | 0.077 |
| rank_03 | 3 | pls_0 | 34.99 | 19 | `ribosomal_mean` | -0.841 | 0.428 | 0.234 | 0.674 | 0.562 |
| rank_03 | 3 | pls_1 | 17.09 | 3 | `stage` | -0.616 | 0.210 | 0.591 | 0.017 | 0.047 |
| rank_03 | 3 | pls_2 | 20.1 | 1 | `stage` | 0.480 | 0.362 | 0.175 | 0.310 | 0.391 |
| rank_04 | 4 | pls_0 | 34.99 | 19 | `ribosomal_mean` | -0.841 | 0.391 | 0.202 | 0.650 | 0.506 |
| rank_04 | 4 | pls_1 | 17.09 | 3 | `stage` | -0.616 | 0.192 | 0.512 | 0.016 | 0.042 |
| rank_04 | 4 | pls_2 | 20.1 | 1 | `stage` | 0.480 | 0.330 | 0.151 | 0.299 | 0.352 |
| rank_04 | 4 | pls_3 | 14.47 | 2 | `log10_n_atac` | -0.442 | 0.087 | 0.134 | 0.035 | 0.100 |
| rank_06 | 6 | pls_0 | 34.99 | 19 | `ribosomal_mean` | -0.841 | 0.384 | 0.197 | 0.625 | 0.473 |
| rank_06 | 6 | pls_1 | 17.09 | 3 | `stage` | -0.616 | 0.188 | 0.499 | 0.016 | 0.039 |
| rank_06 | 6 | pls_2 | 20.1 | 1 | `stage` | 0.480 | 0.324 | 0.147 | 0.287 | 0.329 |
| rank_06 | 6 | pls_3 | 14.47 | 2 | `log10_n_atac` | -0.442 | 0.085 | 0.131 | 0.034 | 0.093 |
| rank_06 | 6 | pls_4 | 15.82 | 1 | `log10_n_atac` | -0.309 | 0.013 | 0.017 | 0.023 | 0.065 |
| rank_06 | 6 | pls_5 | 12.99 | 4 | `log10_n_rna` | -0.351 | 0.005 | 0.009 | 0.016 | 0.000 |
| rank_09 | 9 | pls_0 | 34.99 | 19 | `ribosomal_mean` | -0.841 | 0.379 | 0.164 | 0.536 | 0.326 |
| rank_09 | 9 | pls_1 | 17.09 | 3 | `stage` | -0.616 | 0.186 | 0.415 | 0.013 | 0.027 |
| rank_09 | 9 | pls_2 | 20.1 | 1 | `stage` | 0.480 | 0.320 | 0.123 | 0.246 | 0.227 |
| rank_09 | 9 | pls_3 | 14.47 | 2 | `log10_n_atac` | -0.442 | 0.084 | 0.109 | 0.029 | 0.064 |
| rank_09 | 9 | pls_4 | 15.82 | 1 | `log10_n_atac` | -0.309 | 0.013 | 0.014 | 0.020 | 0.045 |
| rank_09 | 9 | pls_5 | 12.99 | 4 | `log10_n_rna` | -0.351 | 0.005 | 0.007 | 0.014 | 0.000 |
| rank_09 | 9 | pls_6 | 12.23 | 0 | `stage` | -0.263 | 0.008 | 0.028 | 0.006 | 0.011 |
| rank_09 | 9 | pls_7 | 10.8 | 0 | `log10_n_rna` | 0.266 | 0.002 | 0.063 | 0.030 | 0.064 |
| rank_09 | 9 | pls_8 | 12.94 | 0 | `n_libraries_rna` | -0.217 | 0.003 | 0.076 | 0.106 | 0.236 |
| rank_12 | 12 | pls_0 | 34.99 | 19 | `ribosomal_mean` | -0.841 | 0.375 | 0.155 | 0.494 | 0.280 |
| rank_12 | 12 | pls_1 | 17.09 | 3 | `stage` | -0.616 | 0.184 | 0.392 | 0.012 | 0.023 |
| rank_12 | 12 | pls_2 | 20.1 | 1 | `stage` | 0.480 | 0.317 | 0.116 | 0.227 | 0.195 |
| rank_12 | 12 | pls_3 | 14.47 | 2 | `log10_n_atac` | -0.442 | 0.083 | 0.103 | 0.027 | 0.055 |
| rank_12 | 12 | pls_4 | 15.82 | 1 | `log10_n_atac` | -0.309 | 0.013 | 0.014 | 0.018 | 0.038 |
| rank_12 | 12 | pls_5 | 12.99 | 4 | `log10_n_rna` | -0.351 | 0.005 | 0.007 | 0.013 | 0.000 |
| rank_12 | 12 | pls_6 | 12.23 | 0 | `stage` | -0.263 | 0.008 | 0.027 | 0.006 | 0.010 |
| rank_12 | 12 | pls_7 | 10.8 | 0 | `log10_n_rna` | 0.266 | 0.002 | 0.059 | 0.028 | 0.055 |
| rank_12 | 12 | pls_8 | 12.94 | 0 | `n_libraries_rna` | -0.217 | 0.003 | 0.071 | 0.098 | 0.203 |
| rank_12 | 12 | pls_9 | 11.18 | 0 | `doublet_score_sd` | 0.274 | 0.003 | 0.051 | 0.061 | 0.106 |
| rank_12 | 12 | pls_10 | 12.87 | 1 | `atac_median_high_quality_fragments_per_cell_mean` | 0.315 | 0.000 | 0.005 | 0.008 | 0.016 |
| rank_12 | 12 | pls_11 | 10.2 | 0 | `genes_detected_sd` | -0.227 | 0.008 | 0.000 | 0.010 | 0.019 |

### Block principal components (pasted from `diagnostics.md`)

#### atac (variance share: pc1 0.145, pc2 0.110, pc3 0.036)

| variable | n | pc1 rho | pc1 rho \| cov | pc2 rho | pc2 rho \| cov | pc3 rho | pc3 rho \| cov |
|---|---:|---:|---:|---:|---:|---:|---:|
| `age_at_death` | 80 | -0.231 | n/a | -0.128 | n/a | 0.186 | n/a |
| `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | -0.423* | -0.278 | 0.438* | 0.355* | 0.398* | 0.330* |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | -0.449* | -0.323* | 0.613* | 0.512* | 0.343* | 0.366* |
| `atac_median_high_quality_fragments_per_cell_mean` | 80 | -0.232 | 0.161 | 0.134 | -0.049 | -0.058 | -0.003 |
| `atac_percent_duplicates_mean` | 80 | -0.150 | -0.239 | 0.314* | 0.362* | 0.226 | 0.092 |
| `atac_tss_enrichment_score_mean` | 80 | -0.377* | -0.400* | 0.578* | 0.498* | 0.399* | 0.369* |
| `doublet_score_mean` | 80 | 0.126 | 0.107 | 0.028 | -0.052 | 0.024 | 0.020 |
| `doublet_score_sd` | 80 | 0.121 | 0.037 | -0.302* | -0.247 | -0.142 | -0.107 |
| `fraction_mitochondrial_umis_mean` | 80 | 0.404* | 0.214 | -0.547* | -0.466* | -0.233 | -0.169 |
| `fraction_mitochondrial_umis_sd` | 80 | 0.426* | 0.217 | -0.503* | -0.450* | -0.229 | -0.111 |
| `genes_detected_mean` | 80 | -0.145 | -0.087 | 0.004 | -0.033 | 0.264 | 0.217 |
| `genes_detected_sd` | 80 | -0.136 | -0.117 | -0.046 | -0.014 | 0.191 | 0.140 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | -0.161 | -0.118 | 0.416* | 0.351* | 0.300* | 0.119 |
| `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 80 | 0.349* | 0.103 | -0.160 | -0.082 | 0.020 | -0.108 |
| `gex_median_umi_counts_per_cell_mean` | 80 | -0.158 | -0.164 | 0.392* | 0.312* | 0.184 | 0.052 |
| `gex_median_umi_counts_per_cell_sd` | 80 | 0.157 | 0.080 | 0.055 | 0.167 | 0.117 | -0.035 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | -0.180 | -0.291 | 0.270 | 0.226 | 0.445* | 0.275 |
| `gex_reads_mapped_confidently_to_intronic_regions_sd` | 80 | 0.441* | 0.154 | -0.321* | -0.106 | 0.015 | -0.096 |
| `gex_sequencing_saturation_mean` | 80 | 0.454* | 0.310* | -0.461* | -0.285 | -0.176 | -0.136 |
| `gex_sequencing_saturation_sd` | 80 | 0.015 | -0.074 | 0.030 | 0.083 | 0.045 | -0.140 |
| `hemoglobin_mean` | 80 | -0.053 | -0.084 | 0.053 | 0.027 | -0.027 | 0.022 |
| `hemoglobin_pseudobulk` | 80 | -0.103 | -0.094 | 0.089 | 0.034 | -0.076 | -0.009 |
| `hemoglobin_sd` | 80 | -0.037 | -0.014 | 0.031 | 0.090 | -0.019 | -0.055 |
| `hemoglobin_supertype_share` | 75 | -0.184 | 0.012 | 0.097 | 0.123 | -0.048 | 0.008 |
| `log10_lib_atac` | 80 | 0.608* | 0.162 | -0.247 | 0.069 | 0.350* | 0.267 |
| `log10_lib_rna` | 80 | 0.190 | -0.030 | -0.302* | -0.054 | 0.088 | 0.155 |
| `log10_n_atac` | 80 | 0.637* | n/a | -0.335* | n/a | 0.273 | n/a |
| `log10_n_rna` | 80 | 0.243 | n/a | -0.295 | n/a | 0.005 | n/a |
| `male` | 80 | 0.124 | -0.023 | 0.197 | 0.310* | 0.215 | 0.232 |
| `mitochondrial_mean` | 80 | 0.404* | 0.214 | -0.547* | -0.466* | -0.233 | -0.169 |
| `mitochondrial_pseudobulk` | 80 | 0.390* | 0.225 | -0.548* | -0.452* | -0.240 | -0.152 |
| `mitochondrial_sd` | 80 | 0.426* | 0.217 | -0.503* | -0.450* | -0.229 | -0.111 |
| `mitochondrial_supertype_share` | 80 | 0.131 | 0.194 | 0.077 | 0.053 | -0.150 | -0.186 |
| `n_libraries_rna` | 80 | 0.300* | 0.100 | -0.180 | -0.006 | 0.159 | -0.093 |
| `n_nuclei_rna` | 80 | 0.243 | -0.024 | -0.295 | -0.174 | 0.005 | -0.082 |
| `number_of_umis_mean` | 80 | 0.000 | -0.010 | -0.094 | -0.081 | 0.242 | 0.147 |
| `number_of_umis_sd` | 80 | -0.013 | 0.009 | -0.095 | -0.058 | 0.220 | 0.102 |
| `pmi` | 80 | 0.181 | n/a | -0.080 | n/a | -0.150 | n/a |
| `ribosomal_mean` | 80 | -0.370* | -0.308* | 0.571* | 0.449* | 0.135 | 0.194 |
| `ribosomal_pseudobulk` | 80 | -0.354* | -0.299 | 0.555* | 0.450* | 0.125 | 0.165 |
| `ribosomal_sd` | 80 | -0.286 | -0.248 | 0.447* | 0.328* | 0.139 | 0.158 |
| `ribosomal_supertype_share` | 80 | -0.016 | 0.226 | 0.281 | -0.008 | -0.130 | -0.195 |
| `rin` | 80 | 0.038 | n/a | 0.202 | n/a | 0.301* | n/a |
| `stage` | 80 | 0.096 | 0.063 | -0.231 | -0.146 | -0.260 | -0.265 |
| `tmm_atac` | 80 | 0.034 | 0.268 | -0.339* | -0.403* | -0.366* | -0.208 |
| `tmm_rna` | 80 | -0.561* | -0.443* | 0.400* | 0.233 | 0.175 | 0.238 |
| `years_education` | 80 | -0.023 | n/a | 0.129 | n/a | -0.219 | n/a |

#### rna (variance share: pc1 0.277, pc2 0.080, pc3 0.059)

| variable | n | pc1 rho | pc1 rho \| cov | pc2 rho | pc2 rho \| cov | pc3 rho | pc3 rho \| cov |
|---|---:|---:|---:|---:|---:|---:|---:|
| `age_at_death` | 80 | 0.131 | n/a | 0.249 | n/a | 0.105 | n/a |
| `atac_fraction_of_high_quality_fragments_in_cells_mean` | 80 | 0.547* | 0.339* | -0.084 | 0.021 | 0.151 | 0.066 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 80 | 0.708* | 0.583* | -0.095 | -0.003 | 0.208 | 0.096 |
| `atac_median_high_quality_fragments_per_cell_mean` | 80 | 0.021 | -0.197 | -0.087 | 0.036 | -0.118 | 0.096 |
| `atac_percent_duplicates_mean` | 80 | 0.239 | 0.244 | -0.038 | -0.072 | 0.221 | 0.069 |
| `atac_tss_enrichment_score_mean` | 80 | 0.722* | 0.644* | -0.132 | -0.051 | 0.278 | 0.109 |
| `doublet_score_mean` | 80 | -0.111 | -0.198 | -0.039 | -0.022 | 0.246 | 0.127 |
| `doublet_score_sd` | 80 | -0.327* | -0.247 | 0.149 | 0.098 | 0.064 | 0.054 |
| `fraction_mitochondrial_umis_mean` | 80 | -0.734* | -0.552* | 0.265 | 0.173 | -0.143 | 0.001 |
| `fraction_mitochondrial_umis_sd` | 80 | -0.651* | -0.488* | 0.247 | 0.213 | -0.138 | 0.028 |
| `genes_detected_mean` | 80 | 0.081 | 0.027 | 0.248 | 0.272 | 0.098 | 0.105 |
| `genes_detected_sd` | 80 | -0.101 | -0.132 | 0.283 | 0.271 | 0.053 | 0.125 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 80 | 0.593* | 0.396* | -0.124 | -0.098 | 0.137 | -0.083 |
| `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 80 | -0.221 | -0.110 | -0.034 | -0.065 | 0.149 | 0.096 |
| `gex_median_umi_counts_per_cell_mean` | 80 | 0.457* | 0.365* | -0.128 | -0.059 | 0.070 | -0.063 |
| `gex_median_umi_counts_per_cell_sd` | 80 | 0.032 | -0.027 | -0.011 | -0.131 | 0.007 | -0.057 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | 80 | 0.534* | 0.513* | -0.013 | 0.029 | 0.289 | 0.076 |
| `gex_reads_mapped_confidently_to_intronic_regions_sd` | 80 | -0.292 | -0.103 | 0.034 | -0.046 | -0.118 | -0.021 |
| `gex_sequencing_saturation_mean` | 80 | -0.641* | -0.489* | 0.137 | 0.060 | -0.059 | 0.029 |
| `gex_sequencing_saturation_sd` | 80 | 0.036 | 0.060 | 0.093 | 0.030 | -0.052 | -0.026 |
| `hemoglobin_mean` | 80 | 0.124 | 0.122 | -0.006 | -0.049 | 0.040 | -0.013 |
| `hemoglobin_pseudobulk` | 80 | 0.188 | 0.094 | -0.034 | -0.040 | -0.010 | -0.053 |
| `hemoglobin_sd` | 80 | 0.069 | -0.008 | -0.004 | -0.126 | 0.027 | -0.033 |
| `hemoglobin_supertype_share` | 75 | 0.243 | 0.067 | -0.065 | -0.101 | 0.121 | 0.059 |
| `log10_lib_atac` | 80 | -0.203 | -0.021 | 0.172 | 0.055 | 0.172 | 0.203 |
| `log10_lib_rna` | 80 | -0.229 | -0.093 | 0.180 | 0.247 | -0.277 | 0.123 |
| `log10_n_atac` | 80 | -0.259 | n/a | 0.212 | n/a | 0.171 | n/a |
| `log10_n_rna` | 80 | -0.227 | n/a | 0.104 | n/a | -0.353* | n/a |
| `male` | 80 | 0.204 | 0.294 | -0.082 | -0.128 | -0.093 | -0.053 |
| `mitochondrial_mean` | 80 | -0.734* | -0.552* | 0.265 | 0.173 | -0.143 | 0.001 |
| `mitochondrial_pseudobulk` | 80 | -0.706* | -0.514* | 0.262 | 0.138 | -0.143 | -0.009 |
| `mitochondrial_sd` | 80 | -0.651* | -0.488* | 0.247 | 0.213 | -0.138 | 0.028 |
| `mitochondrial_supertype_share` | 80 | -0.201 | -0.274 | -0.061 | -0.077 | 0.033 | -0.096 |
| `n_libraries_rna` | 80 | -0.129 | -0.141 | 0.014 | -0.134 | 0.042 | 0.088 |
| `n_nuclei_rna` | 80 | -0.227 | -0.167 | 0.104 | 0.010 | -0.353* | 0.167 |
| `number_of_umis_mean` | 80 | -0.079 | -0.127 | 0.250 | 0.237 | 0.089 | 0.134 |
| `number_of_umis_sd` | 80 | -0.158 | -0.241 | 0.259 | 0.207 | 0.104 | 0.171 |
| `pmi` | 80 | -0.072 | n/a | 0.037 | n/a | -0.221 | n/a |
| `ribosomal_mean` | 80 | 0.803* | 0.785* | -0.357* | -0.225 | -0.039 | -0.091 |
| `ribosomal_pseudobulk` | 80 | 0.772* | 0.757* | -0.359* | -0.243 | -0.055 | -0.111 |
| `ribosomal_sd` | 80 | 0.691* | 0.661* | -0.365* | -0.221 | -0.014 | -0.056 |
| `ribosomal_supertype_share` | 80 | 0.051 | -0.176 | -0.383* | -0.249 | 0.046 | -0.024 |
| `rin` | 80 | 0.266 | n/a | -0.150 | n/a | 0.291 | n/a |
| `stage` | 80 | -0.316* | -0.237 | -0.020 | -0.098 | -0.141 | -0.158 |
| `tmm_atac` | 80 | -0.480* | -0.494* | 0.021 | -0.002 | -0.074 | 0.089 |
| `tmm_rna` | 80 | 0.749* | 0.690* | -0.042 | 0.057 | 0.128 | -0.057 |
| `years_education` | 80 | 0.028 | n/a | -0.156 | n/a | -0.097 | n/a |

### Geometry shares (pasted from `diagnostics.md`)

Share of the squared length of each vector of the trajectory geometry (from `stage_means.csv`) that lies along each component; `value` is the signed component of a single vector (the net displacements), empty for the summed quantities. At rank 1 every share is 1 by construction.

| analysis | rank | quantity | component | value | along | total | share |
|---|---:|---|---|---:|---:|---:|---:|
| primary | 1 | net_displacement_Female | pls_0 | 32.63 | 1065 | 1065 | 1.000 |
| primary | 1 | net_displacement_Male | pls_0 | 17.69 | 313 | 313 | 1.000 |
| primary | 1 | sex_difference_by_stage | pls_0 | n/a | 929.1 | 929.1 | 1.000 |
| primary | 1 | sex_difference_centred | pls_0 | n/a | 248.6 | 248.6 | 1.000 |
| rank_02 | 2 | net_displacement_Female | pls_0 | 32.63 | 1065 | 1587 | 0.671 |
| rank_02 | 2 | net_displacement_Female | pls_1 | -22.86 | 522.4 | 1587 | 0.329 |
| rank_02 | 2 | net_displacement_Male | pls_0 | 17.69 | 313 | 1105 | 0.283 |
| rank_02 | 2 | net_displacement_Male | pls_1 | -28.14 | 791.8 | 1105 | 0.717 |
| rank_02 | 2 | sex_difference_by_stage | pls_0 | n/a | 929.1 | 952.3 | 0.976 |
| rank_02 | 2 | sex_difference_by_stage | pls_1 | n/a | 23.23 | 952.3 | 0.024 |
| rank_02 | 2 | sex_difference_centred | pls_0 | n/a | 248.6 | 269.2 | 0.923 |
| rank_02 | 2 | sex_difference_centred | pls_1 | n/a | 20.63 | 269.2 | 0.077 |
| rank_03 | 3 | net_displacement_Female | pls_0 | 32.63 | 1065 | 2487 | 0.428 |
| rank_03 | 3 | net_displacement_Female | pls_1 | -22.86 | 522.4 | 2487 | 0.210 |
| rank_03 | 3 | net_displacement_Female | pls_2 | 30 | 899.7 | 2487 | 0.362 |
| rank_03 | 3 | net_displacement_Male | pls_0 | 17.69 | 313 | 1339 | 0.234 |
| rank_03 | 3 | net_displacement_Male | pls_1 | -28.14 | 791.8 | 1339 | 0.591 |
| rank_03 | 3 | net_displacement_Male | pls_2 | 15.29 | 233.9 | 1339 | 0.175 |
| rank_03 | 3 | sex_difference_by_stage | pls_0 | n/a | 929.1 | 1379 | 0.674 |
| rank_03 | 3 | sex_difference_by_stage | pls_1 | n/a | 23.23 | 1379 | 0.017 |
| rank_03 | 3 | sex_difference_by_stage | pls_2 | n/a | 426.9 | 1379 | 0.310 |
| rank_03 | 3 | sex_difference_centred | pls_0 | n/a | 248.6 | 442 | 0.562 |
| rank_03 | 3 | sex_difference_centred | pls_1 | n/a | 20.63 | 442 | 0.047 |
| rank_03 | 3 | sex_difference_centred | pls_2 | n/a | 172.8 | 442 | 0.391 |
| rank_04 | 4 | net_displacement_Female | pls_0 | 32.63 | 1065 | 2724 | 0.391 |
| rank_04 | 4 | net_displacement_Female | pls_1 | -22.86 | 522.4 | 2724 | 0.192 |
| rank_04 | 4 | net_displacement_Female | pls_2 | 30 | 899.7 | 2724 | 0.330 |
| rank_04 | 4 | net_displacement_Female | pls_3 | -15.39 | 236.8 | 2724 | 0.087 |
| rank_04 | 4 | net_displacement_Male | pls_0 | 17.69 | 313 | 1547 | 0.202 |
| rank_04 | 4 | net_displacement_Male | pls_1 | -28.14 | 791.8 | 1547 | 0.512 |
| rank_04 | 4 | net_displacement_Male | pls_2 | 15.29 | 233.9 | 1547 | 0.151 |
| rank_04 | 4 | net_displacement_Male | pls_3 | -14.42 | 208 | 1547 | 0.134 |
| rank_04 | 4 | sex_difference_by_stage | pls_0 | n/a | 929.1 | 1429 | 0.650 |
| rank_04 | 4 | sex_difference_by_stage | pls_1 | n/a | 23.23 | 1429 | 0.016 |
| rank_04 | 4 | sex_difference_by_stage | pls_2 | n/a | 426.9 | 1429 | 0.299 |
| rank_04 | 4 | sex_difference_by_stage | pls_3 | n/a | 50.06 | 1429 | 0.035 |
| rank_04 | 4 | sex_difference_centred | pls_0 | n/a | 248.6 | 490.9 | 0.506 |
| rank_04 | 4 | sex_difference_centred | pls_1 | n/a | 20.63 | 490.9 | 0.042 |
| rank_04 | 4 | sex_difference_centred | pls_2 | n/a | 172.8 | 490.9 | 0.352 |
| rank_04 | 4 | sex_difference_centred | pls_3 | n/a | 48.87 | 490.9 | 0.100 |
| rank_06 | 6 | net_displacement_Female | pls_0 | 32.63 | 1065 | 2773 | 0.384 |
| rank_06 | 6 | net_displacement_Female | pls_1 | -22.86 | 522.4 | 2773 | 0.188 |
| rank_06 | 6 | net_displacement_Female | pls_2 | 30 | 899.7 | 2773 | 0.324 |
| rank_06 | 6 | net_displacement_Female | pls_3 | -15.39 | 236.8 | 2773 | 0.085 |
| rank_06 | 6 | net_displacement_Female | pls_4 | 6.013 | 36.15 | 2773 | 0.013 |
| rank_06 | 6 | net_displacement_Female | pls_5 | 3.584 | 12.85 | 2773 | 0.005 |
| rank_06 | 6 | net_displacement_Male | pls_0 | 17.69 | 313 | 1588 | 0.197 |
| rank_06 | 6 | net_displacement_Male | pls_1 | -28.14 | 791.8 | 1588 | 0.499 |
| rank_06 | 6 | net_displacement_Male | pls_2 | 15.29 | 233.9 | 1588 | 0.147 |
| rank_06 | 6 | net_displacement_Male | pls_3 | -14.42 | 208 | 1588 | 0.131 |
| rank_06 | 6 | net_displacement_Male | pls_4 | 5.226 | 27.31 | 1588 | 0.017 |
| rank_06 | 6 | net_displacement_Male | pls_5 | 3.683 | 13.56 | 1588 | 0.009 |
| rank_06 | 6 | sex_difference_by_stage | pls_0 | n/a | 929.1 | 1487 | 0.625 |
| rank_06 | 6 | sex_difference_by_stage | pls_1 | n/a | 23.23 | 1487 | 0.016 |
| rank_06 | 6 | sex_difference_by_stage | pls_2 | n/a | 426.9 | 1487 | 0.287 |
| rank_06 | 6 | sex_difference_by_stage | pls_3 | n/a | 50.06 | 1487 | 0.034 |
| rank_06 | 6 | sex_difference_by_stage | pls_4 | n/a | 34.1 | 1487 | 0.023 |
| rank_06 | 6 | sex_difference_by_stage | pls_5 | n/a | 23.51 | 1487 | 0.016 |
| rank_06 | 6 | sex_difference_centred | pls_0 | n/a | 248.6 | 525.1 | 0.473 |
| rank_06 | 6 | sex_difference_centred | pls_1 | n/a | 20.63 | 525.1 | 0.039 |
| rank_06 | 6 | sex_difference_centred | pls_2 | n/a | 172.8 | 525.1 | 0.329 |
| rank_06 | 6 | sex_difference_centred | pls_3 | n/a | 48.87 | 525.1 | 0.093 |
| rank_06 | 6 | sex_difference_centred | pls_4 | n/a | 34.04 | 525.1 | 0.065 |
| rank_06 | 6 | sex_difference_centred | pls_5 | n/a | 0.23 | 525.1 | 0.000 |
| rank_09 | 9 | net_displacement_Female | pls_0 | 32.63 | 1065 | 2807 | 0.379 |
| rank_09 | 9 | net_displacement_Female | pls_1 | -22.86 | 522.4 | 2807 | 0.186 |
| rank_09 | 9 | net_displacement_Female | pls_2 | 30 | 899.7 | 2807 | 0.320 |
| rank_09 | 9 | net_displacement_Female | pls_3 | -15.39 | 236.8 | 2807 | 0.084 |
| rank_09 | 9 | net_displacement_Female | pls_4 | 6.013 | 36.15 | 2807 | 0.013 |
| rank_09 | 9 | net_displacement_Female | pls_5 | 3.584 | 12.85 | 2807 | 0.005 |
| rank_09 | 9 | net_displacement_Female | pls_6 | -4.697 | 22.06 | 2807 | 0.008 |
| rank_09 | 9 | net_displacement_Female | pls_7 | 2.226 | 4.953 | 2807 | 0.002 |
| rank_09 | 9 | net_displacement_Female | pls_8 | 2.737 | 7.491 | 2807 | 0.003 |
| rank_09 | 9 | net_displacement_Male | pls_0 | 17.69 | 313 | 1906 | 0.164 |
| rank_09 | 9 | net_displacement_Male | pls_1 | -28.14 | 791.8 | 1906 | 0.415 |
| rank_09 | 9 | net_displacement_Male | pls_2 | 15.29 | 233.9 | 1906 | 0.123 |
| rank_09 | 9 | net_displacement_Male | pls_3 | -14.42 | 208 | 1906 | 0.109 |
| rank_09 | 9 | net_displacement_Male | pls_4 | 5.226 | 27.31 | 1906 | 0.014 |
| rank_09 | 9 | net_displacement_Male | pls_5 | 3.683 | 13.56 | 1906 | 0.007 |
| rank_09 | 9 | net_displacement_Male | pls_6 | -7.337 | 53.84 | 1906 | 0.028 |
| rank_09 | 9 | net_displacement_Male | pls_7 | 10.96 | 120.1 | 1906 | 0.063 |
| rank_09 | 9 | net_displacement_Male | pls_8 | -12.01 | 144.2 | 1906 | 0.076 |
| rank_09 | 9 | sex_difference_by_stage | pls_0 | n/a | 929.1 | 1733 | 0.536 |
| rank_09 | 9 | sex_difference_by_stage | pls_1 | n/a | 23.23 | 1733 | 0.013 |
| rank_09 | 9 | sex_difference_by_stage | pls_2 | n/a | 426.9 | 1733 | 0.246 |
| rank_09 | 9 | sex_difference_by_stage | pls_3 | n/a | 50.06 | 1733 | 0.029 |
| rank_09 | 9 | sex_difference_by_stage | pls_4 | n/a | 34.1 | 1733 | 0.020 |
| rank_09 | 9 | sex_difference_by_stage | pls_5 | n/a | 23.51 | 1733 | 0.014 |
| rank_09 | 9 | sex_difference_by_stage | pls_6 | n/a | 10.75 | 1733 | 0.006 |
| rank_09 | 9 | sex_difference_by_stage | pls_7 | n/a | 52.27 | 1733 | 0.030 |
| rank_09 | 9 | sex_difference_by_stage | pls_8 | n/a | 183.5 | 1733 | 0.106 |
| rank_09 | 9 | sex_difference_centred | pls_0 | n/a | 248.6 | 762.4 | 0.326 |
| rank_09 | 9 | sex_difference_centred | pls_1 | n/a | 20.63 | 762.4 | 0.027 |
| rank_09 | 9 | sex_difference_centred | pls_2 | n/a | 172.8 | 762.4 | 0.227 |
| rank_09 | 9 | sex_difference_centred | pls_3 | n/a | 48.87 | 762.4 | 0.064 |
| rank_09 | 9 | sex_difference_centred | pls_4 | n/a | 34.04 | 762.4 | 0.045 |
| rank_09 | 9 | sex_difference_centred | pls_5 | n/a | 0.23 | 762.4 | 0.000 |
| rank_09 | 9 | sex_difference_centred | pls_6 | n/a | 8.494 | 762.4 | 0.011 |
| rank_09 | 9 | sex_difference_centred | pls_7 | n/a | 49.13 | 762.4 | 0.064 |
| rank_09 | 9 | sex_difference_centred | pls_8 | n/a | 179.7 | 762.4 | 0.236 |
| rank_12 | 12 | net_displacement_Female | pls_0 | 32.63 | 1065 | 2838 | 0.375 |
| rank_12 | 12 | net_displacement_Female | pls_1 | -22.86 | 522.4 | 2838 | 0.184 |
| rank_12 | 12 | net_displacement_Female | pls_2 | 30 | 899.7 | 2838 | 0.317 |
| rank_12 | 12 | net_displacement_Female | pls_3 | -15.39 | 236.8 | 2838 | 0.083 |
| rank_12 | 12 | net_displacement_Female | pls_4 | 6.013 | 36.15 | 2838 | 0.013 |
| rank_12 | 12 | net_displacement_Female | pls_5 | 3.584 | 12.85 | 2838 | 0.005 |
| rank_12 | 12 | net_displacement_Female | pls_6 | -4.697 | 22.06 | 2838 | 0.008 |
| rank_12 | 12 | net_displacement_Female | pls_7 | 2.226 | 4.953 | 2838 | 0.002 |
| rank_12 | 12 | net_displacement_Female | pls_8 | 2.737 | 7.491 | 2838 | 0.003 |
| rank_12 | 12 | net_displacement_Female | pls_9 | -2.668 | 7.116 | 2838 | 0.003 |
| rank_12 | 12 | net_displacement_Female | pls_10 | 0.01603 | 0.0002569 | 2838 | 0.000 |
| rank_12 | 12 | net_displacement_Female | pls_11 | -4.878 | 23.8 | 2838 | 0.008 |
| rank_12 | 12 | net_displacement_Male | pls_0 | 17.69 | 313 | 2020 | 0.155 |
| rank_12 | 12 | net_displacement_Male | pls_1 | -28.14 | 791.8 | 2020 | 0.392 |
| rank_12 | 12 | net_displacement_Male | pls_2 | 15.29 | 233.9 | 2020 | 0.116 |
| rank_12 | 12 | net_displacement_Male | pls_3 | -14.42 | 208 | 2020 | 0.103 |
| rank_12 | 12 | net_displacement_Male | pls_4 | 5.226 | 27.31 | 2020 | 0.014 |
| rank_12 | 12 | net_displacement_Male | pls_5 | 3.683 | 13.56 | 2020 | 0.007 |
| rank_12 | 12 | net_displacement_Male | pls_6 | -7.337 | 53.84 | 2020 | 0.027 |
| rank_12 | 12 | net_displacement_Male | pls_7 | 10.96 | 120.1 | 2020 | 0.059 |
| rank_12 | 12 | net_displacement_Male | pls_8 | -12.01 | 144.2 | 2020 | 0.071 |
| rank_12 | 12 | net_displacement_Male | pls_9 | 10.15 | 103 | 2020 | 0.051 |
| rank_12 | 12 | net_displacement_Male | pls_10 | 3.223 | 10.39 | 2020 | 0.005 |
| rank_12 | 12 | net_displacement_Male | pls_11 | 0.8011 | 0.6418 | 2020 | 0.000 |
| rank_12 | 12 | sex_difference_by_stage | pls_0 | n/a | 929.1 | 1881 | 0.494 |
| rank_12 | 12 | sex_difference_by_stage | pls_1 | n/a | 23.23 | 1881 | 0.012 |
| rank_12 | 12 | sex_difference_by_stage | pls_2 | n/a | 426.9 | 1881 | 0.227 |
| rank_12 | 12 | sex_difference_by_stage | pls_3 | n/a | 50.06 | 1881 | 0.027 |
| rank_12 | 12 | sex_difference_by_stage | pls_4 | n/a | 34.1 | 1881 | 0.018 |
| rank_12 | 12 | sex_difference_by_stage | pls_5 | n/a | 23.51 | 1881 | 0.013 |
| rank_12 | 12 | sex_difference_by_stage | pls_6 | n/a | 10.75 | 1881 | 0.006 |
| rank_12 | 12 | sex_difference_by_stage | pls_7 | n/a | 52.27 | 1881 | 0.028 |
| rank_12 | 12 | sex_difference_by_stage | pls_8 | n/a | 183.5 | 1881 | 0.098 |
| rank_12 | 12 | sex_difference_by_stage | pls_9 | n/a | 114.1 | 1881 | 0.061 |
| rank_12 | 12 | sex_difference_by_stage | pls_10 | n/a | 14.11 | 1881 | 0.008 |
| rank_12 | 12 | sex_difference_by_stage | pls_11 | n/a | 19.16 | 1881 | 0.010 |
| rank_12 | 12 | sex_difference_centred | pls_0 | n/a | 248.6 | 887 | 0.280 |
| rank_12 | 12 | sex_difference_centred | pls_1 | n/a | 20.63 | 887 | 0.023 |
| rank_12 | 12 | sex_difference_centred | pls_2 | n/a | 172.8 | 887 | 0.195 |
| rank_12 | 12 | sex_difference_centred | pls_3 | n/a | 48.87 | 887 | 0.055 |
| rank_12 | 12 | sex_difference_centred | pls_4 | n/a | 34.04 | 887 | 0.038 |
| rank_12 | 12 | sex_difference_centred | pls_5 | n/a | 0.23 | 887 | 0.000 |
| rank_12 | 12 | sex_difference_centred | pls_6 | n/a | 8.494 | 887 | 0.010 |
| rank_12 | 12 | sex_difference_centred | pls_7 | n/a | 49.13 | 887 | 0.055 |
| rank_12 | 12 | sex_difference_centred | pls_8 | n/a | 179.7 | 887 | 0.203 |
| rank_12 | 12 | sex_difference_centred | pls_9 | n/a | 93.86 | 887 | 0.106 |
| rank_12 | 12 | sex_difference_centred | pls_10 | n/a | 13.81 | 887 | 0.016 |
| rank_12 | 12 | sex_difference_centred | pls_11 | n/a | 16.92 | 887 | 0.019 |

### Feature selection (pasted from `diagnostics.md`)

Per omic, from the export's `features_<omic>.csv` and block: input features, features passing the abundance filter, kept features; quantiles of the kept features' mean log2-CPM; the share of kept features with mean below 0 and below 1; the share of the abundance-passing variance held by the kept features; the variance at the cut (kept minimum, dropped median); the share of kept features whose |Spearman| with the primary analysis's `pls_0` is ≥ 0.3.

| omic | input | abundance-passing | kept | mean q05 | mean q25 | mean q50 | mean q75 | mean q95 | share mean < 0 | share mean < 1 | variance share kept | kept min variance | dropped median variance | share \|rho(pls_0)\| ≥ thr |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rna | 36601 | 20387 | 5000 | -0.721 | -0.137 | 0.774 | 2.160 | 4.968 | 0.293 | 0.547 | 0.567 | 0.6459 | 0.2125 | 0.623 (3116) |
| atac | 218882 | 128916 | 5000 | -0.321 | 0.046 | 0.371 | 0.778 | 1.531 | 0.214 | 0.842 | 0.077 | 1.173 | 0.7316 | 0.402 (2009) |

### Analyses read (pasted from `diagnostics.md`)

| analysis | selection | rank | directory |
|---|---|---:|---|
| primary | cv | 1 | `results/motco_analysis/cohort_30` |
| rank_02 | forced | 2 | `results/motco_rank_sensitivity/cohort_30/rank_02` |
| rank_03 | forced | 3 | `results/motco_rank_sensitivity/cohort_30/rank_03` |
| rank_04 | forced | 4 | `results/motco_rank_sensitivity/cohort_30/rank_04` |
| rank_06 | forced | 6 | `results/motco_rank_sensitivity/cohort_30/rank_06` |
| rank_09 | forced | 9 | `results/motco_rank_sensitivity/cohort_30/rank_09` |
| rank_12 | forced | 12 | `results/motco_rank_sensitivity/cohort_30/rank_12` |

### Findings (pasted from `diagnostics.md`)

- variables with missing donors (pairwise-complete correlations): {'hemoglobin_supertype_share': 5}
- pls_0 is identical (np.allclose) across the primary and every forced rank ['rank_02', 'rank_03', 'rank_04', 'rank_06', 'rank_09', 'rank_12']: PLS components nest, so every forced rank keeps the primary's first axis
- primary pls_0: 20 of 47 variables flagged at |rho| ≥ 0.3 (raw or partial): ribosomal_mean (-0.84; partial -0.82), ribosomal_pseudobulk (-0.81; partial -0.79), atac_fraction_of_high_quality_fragments_overlapping_peaks_mean (-0.77; partial -0.64), atac_tss_enrichment_score_mean (-0.77; partial -0.67), fraction_mitochondrial_umis_mean (0.75; partial 0.59), mitochondrial_mean (0.75; partial 0.59), tmm_rna (-0.73; partial -0.67), mitochondrial_pseudobulk (0.73; partial 0.56), ribosomal_sd (-0.70; partial -0.68), gex_sequencing_saturation_mean (0.67; partial 0.52), fraction_mitochondrial_umis_sd (0.65; partial 0.51), mitochondrial_sd (0.65; partial 0.51), gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.61; partial -0.41), atac_fraction_of_high_quality_fragments_in_cells_mean (-0.55; partial -0.37), gex_reads_mapped_confidently_to_intronic_regions_mean (-0.54; partial -0.52), stage (0.49; partial 0.41), tmm_atac (0.49; partial 0.50), gex_median_umi_counts_per_cell_mean (-0.47; partial -0.38), rin (-0.31; partial n/a), male (-0.23; partial -0.33)


## 4. Sensitivity cohort: cohort_50 (n = 74)

Primary analysis `results/motco_analysis/cohort_50` (CV-selected rank 1; `result.json` sha256 `1e7eefe3e63e…`;
`delta` 18.7091, p 0.213; `angle` 0°, p 1; `shape`
0.56471, p 0.377; nothing rejected at α = 0.05), read and checked against the export hashes, never
rerun. Forced-rank analyses read: rank_02 (rank 2), rank_03 (rank 3), rank_04 (rank 4), rank_06 (rank 6), rank_09 (rank 9), rank_12 (rank 12), each checked the same way.
47 variables; threshold |rho| ≥ 0.3 (a reporting flag); no p-values.

**What `pls_0` is, in the measurement space.** The single CV-selected component (score SD 35.76) correlates
-0.840 with the astrocyte nuclei's mean ribosomal-protein UMI fraction (`ribosomal_mean`; -0.797 given the
six analysis covariates) and -0.815 with the UMI-weighted pseudobulk ribosomal fraction, 0.787 with the mean
mitochondrial fraction (`mitochondrial_mean`, identical to SEA-AD's own `Fraction mitochondrial UMIs`; partial 0.605),
-0.774 and -0.785 with the two ATAC library metrics (TSS enrichment and the fraction of fragments in peaks,
partials -0.686 and -0.682), and -0.724 with the RNA TMM factor (partial -0.616), against
0.504 with stage (partial 0.429) and -0.191 with sex (`male`; partial -0.291).
21 of 47 variables are flagged at |rho| ≥ 0.3 on the raw or the partial correlation; the six analysis covariates
leave the technical correlations largely in place (the partial correlations in the table), so the covariate adjustment of the primary
analysis does not remove this axis from the latent space, which MOTCO fits without the covariates (covariate-blind basis).

**Geometry.** The sexes' net displacements along `pls_0` are 37.15 (Female) and 18.44 (Male) latent units; at rank 1 every share
is 1 by construction. `pls_0` is identical (`np.allclose`) across the primary and all six forced ranks, so the forced-rank
sweep inherited this axis as its first and largest component: the share of the squared per-stage sex-difference geometry along `pls_0`
is 0.987 at rank 2, 0.700 at rank 4 and 0.459 at rank 12 (stage-centred 0.973, 0.719, 0.369); of the
female net displacement 0.662, 0.534, 0.469 and of the male 0.270, 0.253, 0.162. `delta` in the primary analysis
measured the difference between those two displacements along this axis.

**Block axes.** The first principal component of the z-scored RNA block holds 0.282 of its variance and correlates
0.816 with `ribosomal_mean`, 0.791 with `ribosomal_pseudobulk` and -0.762 with `fraction_mitochondrial_umis_mean`
(19 variables flagged); the ATAC block's first component holds 0.150 of its variance and correlates
-0.611 with `tmm_rna`, 0.607 with `log10_n_atac` and 0.546 with `log10_lib_atac`
(18 flagged). The composition axis is therefore the dominant axis of the exported RNA block itself, not an artefact of
the PLS fit.

**Feature selection.** Of 20,238 abundance-passing genes, the 5,000 kept by the top-N variance rule have a median mean
log2-CPM of 0.82 (5–95 % quantiles -0.66 to 4.96), 28.2 % of them below 0 and 53.8 % below 1,
and hold 57.4 % of the abundance-passing variance (cut at 0.634; dropped median 0.207);
62.9 % (3,145) of the kept genes have |rho| ≥ 0.3 with `pls_0`. Of 124,818 abundance-passing
peaks, the 5,000 kept have a median mean of 0.27, 85.8 % below 1, and hold 8.2 % of the variance (cut
1.156 against a dropped median of 0.689: a near-flat cut); 46.2 % (2,308) correlate with `pls_0`
at that level.

Findings recorded by the step: (1) variables with missing donors (pairwise-complete correlations): {'hemoglobin_supertype_share': 5} (2) pls_0 is identical (np.allclose) across the primary and every forced rank ['rank_02', 'rank_03', 'rank_04', 'rank_06', 'rank_09', 'rank_12']: PLS components nest, so every forced rank keeps the primary's first axis (3) primary pls_0: 21 of 47 variables flagged at |rho| ≥ 0.3 (raw or partial): ribosomal_mean (-0.84; partial -0.80), ribosomal_pseudobulk (-0.82; partial -0.77), fraction_mitochondrial_umis_mean (0.79; partial 0.60), mitochondrial_mean (0.79; partial 0.60), atac_fraction_of_high_quality_fragments_overlapping_peaks_mean (-0.78; partial -0.68), atac_tss_enrichment_score_mean (-0.77; partial -0.69), mitochondrial_pseudobulk (0.77; partial 0.57), gex_sequencing_saturation_mean (0.73; partial 0.56), tmm_rna (-0.72; partial -0.62), ribosomal_sd (-0.69; partial -0.66), fraction_mitochondrial_umis_sd (0.68; partial 0.49), mitochondrial_sd (0.68; partial 0.49), gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.61; partial -0.40), atac_fraction_of_high_quality_fragments_in_cells_mean (-0.60; partial -0.44), gex_reads_mapped_confidently_to_intronic_regions_mean (-0.56; partial -0.54), tmm_atac (0.56; partial 0.58), gex_median_umi_counts_per_cell_mean (-0.51; partial -0.42), stage (0.50; partial 0.43), gex_reads_mapped_confidently_to_intronic_regions_sd (0.33; partial 0.09), log10_n_atac (0.32; partial n/a), atac_percent_duplicates_mean (-0.30; partial -0.33)

Descriptive measurement-space diagnostics of the covariate-adjusted trajectory analysis: Spearman correlations of every latent component of the primary analysis (CV-selected rank 1) and of the forced-rank analyses (2, 3, 4, 6, 9, 12) with 47 predeclared variables (technical table, export TMM factors and library sizes, the six analysis covariates, stage and sex), the partial correlations given the six analysis covariates (Spearman correlation of the OLS residuals of both variables on an intercept and the covariates), the same for the first 3 principal components of each z-scored block, the share of the trajectory geometry along each component, and the feature-selection accounting. Rows with |rho| ≥ 0.3 are flagged and listed first; the flag is a reporting threshold, not an adoption rule. No p-value is attached to any correlation, the tables support no inferential claim, and nothing here chooses a preprocessing variant: an axis correlated with a technical variable is a finding about the measurement space of the primary analysis, not a reinterpretation of its statistics, and the primary analysis stays primary with its non-rejection reading rules.

### Flagged correlations (pasted from `diagnostics.md`)

| analysis | rank | component | sd | variable | n | rho | rho \| covariates | flag |
|---|---:|---|---:|---|---:|---:|---:|---|
| primary | 1 | pls_0 | 35.76 | `ribosomal_mean` | 74 | -0.840 | -0.797 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `ribosomal_pseudobulk` | 74 | -0.815 | -0.769 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `fraction_mitochondrial_umis_mean` | 74 | 0.787 | 0.605 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `mitochondrial_mean` | 74 | 0.787 | 0.605 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | -0.785 | -0.682 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `atac_tss_enrichment_score_mean` | 74 | -0.774 | -0.686 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `mitochondrial_pseudobulk` | 74 | 0.766 | 0.573 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `gex_sequencing_saturation_mean` | 74 | 0.727 | 0.560 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `tmm_rna` | 74 | -0.724 | -0.616 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `ribosomal_sd` | 74 | -0.694 | -0.663 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `fraction_mitochondrial_umis_sd` | 74 | 0.683 | 0.491 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `mitochondrial_sd` | 74 | 0.683 | 0.491 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | -0.614 | -0.401 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | -0.598 | -0.444 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | -0.563 | -0.541 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `tmm_atac` | 74 | 0.557 | 0.581 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `gex_median_umi_counts_per_cell_mean` | 74 | -0.511 | -0.417 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `stage` | 74 | 0.504 | 0.429 | raw+partial |
| primary | 1 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | 0.326 | 0.090 | raw |
| primary | 1 | pls_0 | 35.76 | `log10_n_atac` | 74 | 0.324 | n/a | raw |
| primary | 1 | pls_0 | 35.76 | `atac_percent_duplicates_mean` | 74 | -0.296 | -0.334 | partial |
| rank_02 | 2 | pls_0 | 35.76 | `ribosomal_mean` | 74 | -0.840 | -0.797 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `ribosomal_pseudobulk` | 74 | -0.815 | -0.769 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `fraction_mitochondrial_umis_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `mitochondrial_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | -0.785 | -0.682 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `atac_tss_enrichment_score_mean` | 74 | -0.774 | -0.686 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `mitochondrial_pseudobulk` | 74 | 0.766 | 0.573 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `gex_sequencing_saturation_mean` | 74 | 0.727 | 0.560 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `tmm_rna` | 74 | -0.724 | -0.616 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `ribosomal_sd` | 74 | -0.694 | -0.663 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `fraction_mitochondrial_umis_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `mitochondrial_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | -0.614 | -0.401 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | -0.598 | -0.444 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | -0.563 | -0.541 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `tmm_atac` | 74 | 0.557 | 0.581 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `gex_median_umi_counts_per_cell_mean` | 74 | -0.511 | -0.417 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `stage` | 74 | 0.504 | 0.429 | raw+partial |
| rank_02 | 2 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | 0.326 | 0.090 | raw |
| rank_02 | 2 | pls_0 | 35.76 | `log10_n_atac` | 74 | 0.324 | n/a | raw |
| rank_02 | 2 | pls_0 | 35.76 | `atac_percent_duplicates_mean` | 74 | -0.296 | -0.334 | partial |
| rank_02 | 2 | pls_1 | 18.17 | `stage` | 74 | -0.616 | -0.761 | raw+partial |
| rank_02 | 2 | pls_1 | 18.17 | `log10_lib_atac` | 74 | 0.496 | 0.363 | raw+partial |
| rank_02 | 2 | pls_1 | 18.17 | `log10_n_atac` | 74 | 0.403 | n/a | raw |
| rank_03 | 3 | pls_0 | 35.76 | `ribosomal_mean` | 74 | -0.840 | -0.797 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `ribosomal_pseudobulk` | 74 | -0.815 | -0.769 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `fraction_mitochondrial_umis_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `mitochondrial_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | -0.785 | -0.682 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `atac_tss_enrichment_score_mean` | 74 | -0.774 | -0.686 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `mitochondrial_pseudobulk` | 74 | 0.766 | 0.573 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `gex_sequencing_saturation_mean` | 74 | 0.727 | 0.560 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `tmm_rna` | 74 | -0.724 | -0.616 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `ribosomal_sd` | 74 | -0.694 | -0.663 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `fraction_mitochondrial_umis_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `mitochondrial_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | -0.614 | -0.401 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | -0.598 | -0.444 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | -0.563 | -0.541 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `tmm_atac` | 74 | 0.557 | 0.581 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `gex_median_umi_counts_per_cell_mean` | 74 | -0.511 | -0.417 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `stage` | 74 | 0.504 | 0.429 | raw+partial |
| rank_03 | 3 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | 0.326 | 0.090 | raw |
| rank_03 | 3 | pls_0 | 35.76 | `log10_n_atac` | 74 | 0.324 | n/a | raw |
| rank_03 | 3 | pls_0 | 35.76 | `atac_percent_duplicates_mean` | 74 | -0.296 | -0.334 | partial |
| rank_03 | 3 | pls_1 | 18.17 | `stage` | 74 | -0.616 | -0.761 | raw+partial |
| rank_03 | 3 | pls_1 | 18.17 | `log10_lib_atac` | 74 | 0.496 | 0.363 | raw+partial |
| rank_03 | 3 | pls_1 | 18.17 | `log10_n_atac` | 74 | 0.403 | n/a | raw |
| rank_03 | 3 | pls_2 | 17.53 | `atac_percent_duplicates_mean` | 74 | 0.345 | 0.312 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `ribosomal_mean` | 74 | -0.840 | -0.797 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `ribosomal_pseudobulk` | 74 | -0.815 | -0.769 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `fraction_mitochondrial_umis_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `mitochondrial_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | -0.785 | -0.682 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `atac_tss_enrichment_score_mean` | 74 | -0.774 | -0.686 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `mitochondrial_pseudobulk` | 74 | 0.766 | 0.573 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `gex_sequencing_saturation_mean` | 74 | 0.727 | 0.560 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `tmm_rna` | 74 | -0.724 | -0.616 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `ribosomal_sd` | 74 | -0.694 | -0.663 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `fraction_mitochondrial_umis_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `mitochondrial_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | -0.614 | -0.401 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | -0.598 | -0.444 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | -0.563 | -0.541 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `tmm_atac` | 74 | 0.557 | 0.581 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `gex_median_umi_counts_per_cell_mean` | 74 | -0.511 | -0.417 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `stage` | 74 | 0.504 | 0.429 | raw+partial |
| rank_04 | 4 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | 0.326 | 0.090 | raw |
| rank_04 | 4 | pls_0 | 35.76 | `log10_n_atac` | 74 | 0.324 | n/a | raw |
| rank_04 | 4 | pls_0 | 35.76 | `atac_percent_duplicates_mean` | 74 | -0.296 | -0.334 | partial |
| rank_04 | 4 | pls_1 | 18.17 | `stage` | 74 | -0.616 | -0.761 | raw+partial |
| rank_04 | 4 | pls_1 | 18.17 | `log10_lib_atac` | 74 | 0.496 | 0.363 | raw+partial |
| rank_04 | 4 | pls_1 | 18.17 | `log10_n_atac` | 74 | 0.403 | n/a | raw |
| rank_04 | 4 | pls_2 | 17.53 | `atac_percent_duplicates_mean` | 74 | 0.345 | 0.312 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `ribosomal_mean` | 74 | -0.840 | -0.797 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `ribosomal_pseudobulk` | 74 | -0.815 | -0.769 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `fraction_mitochondrial_umis_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `mitochondrial_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | -0.785 | -0.682 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `atac_tss_enrichment_score_mean` | 74 | -0.774 | -0.686 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `mitochondrial_pseudobulk` | 74 | 0.766 | 0.573 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `gex_sequencing_saturation_mean` | 74 | 0.727 | 0.560 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `tmm_rna` | 74 | -0.724 | -0.616 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `ribosomal_sd` | 74 | -0.694 | -0.663 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `fraction_mitochondrial_umis_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `mitochondrial_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | -0.614 | -0.401 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | -0.598 | -0.444 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | -0.563 | -0.541 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `tmm_atac` | 74 | 0.557 | 0.581 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `gex_median_umi_counts_per_cell_mean` | 74 | -0.511 | -0.417 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `stage` | 74 | 0.504 | 0.429 | raw+partial |
| rank_06 | 6 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | 0.326 | 0.090 | raw |
| rank_06 | 6 | pls_0 | 35.76 | `log10_n_atac` | 74 | 0.324 | n/a | raw |
| rank_06 | 6 | pls_0 | 35.76 | `atac_percent_duplicates_mean` | 74 | -0.296 | -0.334 | partial |
| rank_06 | 6 | pls_1 | 18.17 | `stage` | 74 | -0.616 | -0.761 | raw+partial |
| rank_06 | 6 | pls_1 | 18.17 | `log10_lib_atac` | 74 | 0.496 | 0.363 | raw+partial |
| rank_06 | 6 | pls_1 | 18.17 | `log10_n_atac` | 74 | 0.403 | n/a | raw |
| rank_06 | 6 | pls_2 | 17.53 | `atac_percent_duplicates_mean` | 74 | 0.345 | 0.312 | raw+partial |
| rank_06 | 6 | pls_4 | 15.74 | `log10_n_atac` | 74 | 0.525 | n/a | raw |
| rank_06 | 6 | pls_4 | 15.74 | `stage` | 74 | 0.483 | 0.559 | raw+partial |
| rank_06 | 6 | pls_4 | 15.74 | `log10_lib_atac` | 74 | 0.404 | -0.034 | raw |
| rank_06 | 6 | pls_4 | 15.74 | `atac_median_high_quality_fragments_per_cell_mean` | 74 | -0.333 | -0.144 | raw |
| rank_06 | 6 | pls_5 | 13.3 | `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 74 | -0.327 | -0.145 | raw |
| rank_06 | 6 | pls_5 | 13.3 | `log10_lib_rna` | 74 | -0.314 | -0.158 | raw |
| rank_06 | 6 | pls_5 | 13.3 | `doublet_score_mean` | 74 | 0.310 | 0.266 | raw |
| rank_09 | 9 | pls_0 | 35.76 | `ribosomal_mean` | 74 | -0.840 | -0.797 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `ribosomal_pseudobulk` | 74 | -0.815 | -0.769 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `fraction_mitochondrial_umis_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `mitochondrial_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | -0.785 | -0.682 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `atac_tss_enrichment_score_mean` | 74 | -0.774 | -0.686 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `mitochondrial_pseudobulk` | 74 | 0.766 | 0.573 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `gex_sequencing_saturation_mean` | 74 | 0.727 | 0.560 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `tmm_rna` | 74 | -0.724 | -0.616 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `ribosomal_sd` | 74 | -0.694 | -0.663 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `fraction_mitochondrial_umis_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `mitochondrial_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | -0.614 | -0.401 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | -0.598 | -0.444 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | -0.563 | -0.541 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `tmm_atac` | 74 | 0.557 | 0.581 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `gex_median_umi_counts_per_cell_mean` | 74 | -0.511 | -0.417 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `stage` | 74 | 0.504 | 0.429 | raw+partial |
| rank_09 | 9 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | 0.326 | 0.090 | raw |
| rank_09 | 9 | pls_0 | 35.76 | `log10_n_atac` | 74 | 0.324 | n/a | raw |
| rank_09 | 9 | pls_0 | 35.76 | `atac_percent_duplicates_mean` | 74 | -0.296 | -0.334 | partial |
| rank_09 | 9 | pls_1 | 18.17 | `stage` | 74 | -0.616 | -0.761 | raw+partial |
| rank_09 | 9 | pls_1 | 18.17 | `log10_lib_atac` | 74 | 0.496 | 0.363 | raw+partial |
| rank_09 | 9 | pls_1 | 18.17 | `log10_n_atac` | 74 | 0.403 | n/a | raw |
| rank_09 | 9 | pls_2 | 17.53 | `atac_percent_duplicates_mean` | 74 | 0.345 | 0.312 | raw+partial |
| rank_09 | 9 | pls_4 | 15.74 | `log10_n_atac` | 74 | 0.525 | n/a | raw |
| rank_09 | 9 | pls_4 | 15.74 | `stage` | 74 | 0.483 | 0.559 | raw+partial |
| rank_09 | 9 | pls_4 | 15.74 | `log10_lib_atac` | 74 | 0.404 | -0.034 | raw |
| rank_09 | 9 | pls_4 | 15.74 | `atac_median_high_quality_fragments_per_cell_mean` | 74 | -0.333 | -0.144 | raw |
| rank_09 | 9 | pls_5 | 13.3 | `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 74 | -0.327 | -0.145 | raw |
| rank_09 | 9 | pls_5 | 13.3 | `log10_lib_rna` | 74 | -0.314 | -0.158 | raw |
| rank_09 | 9 | pls_5 | 13.3 | `doublet_score_mean` | 74 | 0.310 | 0.266 | raw |
| rank_12 | 12 | pls_0 | 35.76 | `ribosomal_mean` | 74 | -0.840 | -0.797 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `ribosomal_pseudobulk` | 74 | -0.815 | -0.769 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `fraction_mitochondrial_umis_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `mitochondrial_mean` | 74 | 0.787 | 0.605 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | -0.785 | -0.682 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `atac_tss_enrichment_score_mean` | 74 | -0.774 | -0.686 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `mitochondrial_pseudobulk` | 74 | 0.766 | 0.573 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `gex_sequencing_saturation_mean` | 74 | 0.727 | 0.560 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `tmm_rna` | 74 | -0.724 | -0.616 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `ribosomal_sd` | 74 | -0.694 | -0.663 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `fraction_mitochondrial_umis_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `mitochondrial_sd` | 74 | 0.683 | 0.491 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | -0.614 | -0.401 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | -0.598 | -0.444 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | -0.563 | -0.541 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `tmm_atac` | 74 | 0.557 | 0.581 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `gex_median_umi_counts_per_cell_mean` | 74 | -0.511 | -0.417 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `stage` | 74 | 0.504 | 0.429 | raw+partial |
| rank_12 | 12 | pls_0 | 35.76 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | 0.326 | 0.090 | raw |
| rank_12 | 12 | pls_0 | 35.76 | `log10_n_atac` | 74 | 0.324 | n/a | raw |
| rank_12 | 12 | pls_0 | 35.76 | `atac_percent_duplicates_mean` | 74 | -0.296 | -0.334 | partial |
| rank_12 | 12 | pls_1 | 18.17 | `stage` | 74 | -0.616 | -0.761 | raw+partial |
| rank_12 | 12 | pls_1 | 18.17 | `log10_lib_atac` | 74 | 0.496 | 0.363 | raw+partial |
| rank_12 | 12 | pls_1 | 18.17 | `log10_n_atac` | 74 | 0.403 | n/a | raw |
| rank_12 | 12 | pls_2 | 17.53 | `atac_percent_duplicates_mean` | 74 | 0.345 | 0.312 | raw+partial |
| rank_12 | 12 | pls_4 | 15.74 | `log10_n_atac` | 74 | 0.525 | n/a | raw |
| rank_12 | 12 | pls_4 | 15.74 | `stage` | 74 | 0.483 | 0.559 | raw+partial |
| rank_12 | 12 | pls_4 | 15.74 | `log10_lib_atac` | 74 | 0.404 | -0.034 | raw |
| rank_12 | 12 | pls_4 | 15.74 | `atac_median_high_quality_fragments_per_cell_mean` | 74 | -0.333 | -0.144 | raw |
| rank_12 | 12 | pls_5 | 13.3 | `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 74 | -0.327 | -0.145 | raw |
| rank_12 | 12 | pls_5 | 13.3 | `log10_lib_rna` | 74 | -0.314 | -0.158 | raw |
| rank_12 | 12 | pls_5 | 13.3 | `doublet_score_mean` | 74 | 0.310 | 0.266 | raw |
| rank_12 | 12 | pls_9 | 10.88 | `atac_percent_duplicates_mean` | 74 | -0.242 | -0.308 | partial |
| rank_12 | 12 | pls_9 | 10.88 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | -0.166 | -0.359 | partial |
| rank_12 | 12 | pls_9 | 10.88 | `tmm_atac` | 74 | 0.164 | 0.320 | partial |

Block principal components:

| block | component | variance share | variable | n | rho | rho \| covariates | flag |
|---|---|---:|---|---:|---:|---:|---|
| atac | pc1 | 0.150 | `tmm_rna` | 74 | -0.611 | -0.442 | raw+partial |
| atac | pc1 | 0.150 | `log10_n_atac` | 74 | 0.607 | n/a | raw |
| atac | pc1 | 0.150 | `log10_lib_atac` | 74 | 0.546 | 0.113 | raw |
| atac | pc1 | 0.150 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | 0.512 | 0.227 | raw |
| atac | pc1 | 0.150 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | -0.503 | -0.348 | raw+partial |
| atac | pc1 | 0.150 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | -0.498 | -0.393 | raw+partial |
| atac | pc1 | 0.150 | `gex_sequencing_saturation_mean` | 74 | 0.441 | 0.302 | raw+partial |
| atac | pc1 | 0.150 | `atac_tss_enrichment_score_mean` | 74 | -0.436 | -0.441 | raw+partial |
| atac | pc1 | 0.150 | `fraction_mitochondrial_umis_sd` | 74 | 0.433 | 0.228 | raw |
| atac | pc1 | 0.150 | `mitochondrial_sd` | 74 | 0.433 | 0.228 | raw |
| atac | pc1 | 0.150 | `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 74 | 0.433 | 0.167 | raw |
| atac | pc1 | 0.150 | `fraction_mitochondrial_umis_mean` | 74 | 0.420 | 0.246 | raw |
| atac | pc1 | 0.150 | `mitochondrial_mean` | 74 | 0.420 | 0.246 | raw |
| atac | pc1 | 0.150 | `ribosomal_mean` | 74 | -0.416 | -0.321 | raw+partial |
| atac | pc1 | 0.150 | `mitochondrial_pseudobulk` | 74 | 0.398 | 0.247 | raw |
| atac | pc1 | 0.150 | `ribosomal_pseudobulk` | 74 | -0.395 | -0.305 | raw+partial |
| atac | pc1 | 0.150 | `ribosomal_sd` | 74 | -0.327 | -0.269 | raw |
| atac | pc1 | 0.150 | `n_libraries_rna` | 74 | 0.324 | 0.192 | raw |
| atac | pc1 | 0.150 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | -0.236 | -0.325 | partial |
| atac | pc1 | 0.150 | `tmm_atac` | 74 | 0.195 | 0.335 | partial |
| atac | pc2 | 0.123 | `ribosomal_mean` | 74 | 0.627 | 0.514 | raw+partial |
| atac | pc2 | 0.123 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | 0.620 | 0.518 | raw+partial |
| atac | pc2 | 0.123 | `ribosomal_pseudobulk` | 74 | 0.616 | 0.518 | raw+partial |
| atac | pc2 | 0.123 | `fraction_mitochondrial_umis_mean` | 74 | -0.598 | -0.502 | raw+partial |
| atac | pc2 | 0.123 | `mitochondrial_mean` | 74 | -0.598 | -0.502 | raw+partial |
| atac | pc2 | 0.123 | `atac_tss_enrichment_score_mean` | 74 | 0.598 | 0.505 | raw+partial |
| atac | pc2 | 0.123 | `mitochondrial_pseudobulk` | 74 | -0.593 | -0.473 | raw+partial |
| atac | pc2 | 0.123 | `fraction_mitochondrial_umis_sd` | 74 | -0.556 | -0.500 | raw+partial |
| atac | pc2 | 0.123 | `mitochondrial_sd` | 74 | -0.556 | -0.500 | raw+partial |
| atac | pc2 | 0.123 | `gex_sequencing_saturation_mean` | 74 | -0.540 | -0.392 | raw+partial |
| atac | pc2 | 0.123 | `ribosomal_sd` | 74 | 0.504 | 0.390 | raw+partial |
| atac | pc2 | 0.123 | `tmm_rna` | 74 | 0.487 | 0.338 | raw+partial |
| atac | pc2 | 0.123 | `gex_median_umi_counts_per_cell_mean` | 74 | 0.475 | 0.391 | raw+partial |
| atac | pc2 | 0.123 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | 0.449 | 0.333 | raw+partial |
| atac | pc2 | 0.123 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | 0.448 | 0.370 | raw+partial |
| atac | pc2 | 0.123 | `tmm_atac` | 74 | -0.425 | -0.443 | raw+partial |
| atac | pc2 | 0.123 | `atac_percent_duplicates_mean` | 74 | 0.334 | 0.386 | raw+partial |
| atac | pc2 | 0.123 | `doublet_score_sd` | 74 | -0.331 | -0.290 | raw |
| atac | pc2 | 0.123 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | 0.316 | 0.254 | raw |
| atac | pc2 | 0.123 | `log10_n_atac` | 74 | -0.313 | n/a | raw |
| atac | pc2 | 0.123 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | -0.306 | -0.114 | raw |
| atac | pc2 | 0.123 | `male` | 74 | 0.226 | 0.307 | partial |
| atac | pc3 | 0.034 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | 0.504 | 0.257 | raw |
| atac | pc3 | 0.034 | `log10_lib_atac` | 74 | 0.415 | 0.298 | raw |
| atac | pc3 | 0.034 | `tmm_atac` | 74 | -0.406 | -0.172 | raw |
| atac | pc3 | 0.034 | `rin` | 74 | 0.403 | n/a | raw |
| atac | pc3 | 0.034 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | 0.374 | 0.258 | raw |
| atac | pc3 | 0.034 | `stage` | 74 | -0.360 | -0.337 | raw+partial |
| atac | pc3 | 0.034 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | 0.354 | 0.093 | raw |
| atac | pc3 | 0.034 | `atac_tss_enrichment_score_mean` | 74 | 0.354 | 0.284 | raw |
| atac | pc3 | 0.034 | `log10_n_atac` | 74 | 0.301 | n/a | raw |
| rna | pc1 | 0.282 | `ribosomal_mean` | 74 | 0.816 | 0.775 | raw+partial |
| rna | pc1 | 0.282 | `ribosomal_pseudobulk` | 74 | 0.791 | 0.747 | raw+partial |
| rna | pc1 | 0.282 | `fraction_mitochondrial_umis_mean` | 74 | -0.762 | -0.549 | raw+partial |
| rna | pc1 | 0.282 | `mitochondrial_mean` | 74 | -0.762 | -0.549 | raw+partial |
| rna | pc1 | 0.282 | `mitochondrial_pseudobulk` | 74 | -0.736 | -0.517 | raw+partial |
| rna | pc1 | 0.282 | `tmm_rna` | 74 | 0.731 | 0.616 | raw+partial |
| rna | pc1 | 0.282 | `atac_tss_enrichment_score_mean` | 74 | 0.729 | 0.667 | raw+partial |
| rna | pc1 | 0.282 | `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | 0.720 | 0.631 | raw+partial |
| rna | pc1 | 0.282 | `ribosomal_sd` | 74 | 0.698 | 0.653 | raw+partial |
| rna | pc1 | 0.282 | `fraction_mitochondrial_umis_sd` | 74 | -0.686 | -0.466 | raw+partial |
| rna | pc1 | 0.282 | `gex_sequencing_saturation_mean` | 74 | -0.686 | -0.515 | raw+partial |
| rna | pc1 | 0.282 | `mitochondrial_sd` | 74 | -0.686 | -0.466 | raw+partial |
| rna | pc1 | 0.282 | `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | 0.602 | 0.390 | raw+partial |
| rna | pc1 | 0.282 | `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | 0.591 | 0.397 | raw+partial |
| rna | pc1 | 0.282 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | 0.550 | 0.520 | raw+partial |
| rna | pc1 | 0.282 | `tmm_atac` | 74 | -0.545 | -0.550 | raw+partial |
| rna | pc1 | 0.282 | `gex_median_umi_counts_per_cell_mean` | 74 | 0.485 | 0.369 | raw+partial |
| rna | pc1 | 0.282 | `doublet_score_sd` | 74 | -0.345 | -0.252 | raw |
| rna | pc1 | 0.282 | `stage` | 74 | -0.309 | -0.215 | raw |
| rna | pc1 | 0.282 | `atac_percent_duplicates_mean` | 74 | 0.297 | 0.306 | partial |
| rna | pc1 | 0.282 | `mitochondrial_supertype_share` | 74 | -0.207 | -0.310 | partial |
| rna | pc2 | 0.085 | `ribosomal_sd` | 74 | -0.411 | -0.216 | raw |
| rna | pc2 | 0.085 | `ribosomal_pseudobulk` | 74 | -0.410 | -0.281 | raw |
| rna | pc2 | 0.085 | `ribosomal_mean` | 74 | -0.405 | -0.259 | raw |
| rna | pc2 | 0.085 | `mitochondrial_pseudobulk` | 74 | 0.352 | 0.235 | raw |
| rna | pc2 | 0.085 | `fraction_mitochondrial_umis_mean` | 74 | 0.348 | 0.271 | raw |
| rna | pc2 | 0.085 | `mitochondrial_mean` | 74 | 0.348 | 0.271 | raw |
| rna | pc2 | 0.085 | `ribosomal_supertype_share` | 74 | -0.346 | -0.177 | raw |
| rna | pc2 | 0.085 | `fraction_mitochondrial_umis_sd` | 74 | 0.342 | 0.317 | raw+partial |
| rna | pc2 | 0.085 | `mitochondrial_sd` | 74 | 0.342 | 0.317 | raw+partial |
| rna | pc3 | 0.060 | `rin` | 74 | 0.361 | n/a | raw |
| rna | pc3 | 0.060 | `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | 0.338 | 0.081 | raw |

### Primary analysis (`results/motco_analysis/cohort_50`, cv rank 1) (pasted from `diagnostics.md`)

Component score SD: pls_0 35.76. Columns per component: Spearman rho, then rho given the six analysis covariates (n/a for a covariate itself). Variables sorted by name; `n` is the number of donors with the variable.

| variable | source | n | pls_0 rho | pls_0 rho \| cov |
|---|---|---:|---:|---:|
| `age_at_death` | covariates.csv | 74 | -0.122 | n/a |
| `atac_fraction_of_high_quality_fragments_in_cells_mean` | donor_technical.csv | 74 | -0.598* | -0.444* |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | donor_technical.csv | 74 | -0.785* | -0.682* |
| `atac_median_high_quality_fragments_per_cell_mean` | donor_technical.csv | 74 | -0.050 | 0.189 |
| `atac_percent_duplicates_mean` | donor_technical.csv | 74 | -0.296 | -0.334* |
| `atac_tss_enrichment_score_mean` | donor_technical.csv | 74 | -0.774* | -0.686* |
| `doublet_score_mean` | donor_technical.csv | 74 | 0.014 | 0.097 |
| `doublet_score_sd` | donor_technical.csv | 74 | 0.295 | 0.198 |
| `fraction_mitochondrial_umis_mean` | donor_technical.csv | 74 | 0.787* | 0.605* |
| `fraction_mitochondrial_umis_sd` | donor_technical.csv | 74 | 0.683* | 0.491* |
| `genes_detected_mean` | donor_technical.csv | 74 | -0.116 | -0.077 |
| `genes_detected_sd` | donor_technical.csv | 74 | 0.023 | 0.001 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | donor_technical.csv | 74 | -0.614* | -0.401* |
| `gex_fraction_of_transcriptomic_reads_in_cells_sd` | donor_technical.csv | 74 | 0.256 | 0.125 |
| `gex_median_umi_counts_per_cell_mean` | donor_technical.csv | 74 | -0.511* | -0.417* |
| `gex_median_umi_counts_per_cell_sd` | donor_technical.csv | 74 | -0.090 | -0.053 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | donor_technical.csv | 74 | -0.563* | -0.541* |
| `gex_reads_mapped_confidently_to_intronic_regions_sd` | donor_technical.csv | 74 | 0.326* | 0.090 |
| `gex_sequencing_saturation_mean` | donor_technical.csv | 74 | 0.727* | 0.560* |
| `gex_sequencing_saturation_sd` | donor_technical.csv | 74 | -0.121 | -0.085 |
| `hemoglobin_mean` | donor_technical.csv | 74 | -0.082 | -0.046 |
| `hemoglobin_pseudobulk` | donor_technical.csv | 74 | -0.132 | -0.008 |
| `hemoglobin_sd` | donor_technical.csv | 74 | -0.043 | 0.008 |
| `hemoglobin_supertype_share` | donor_technical.csv | 69 | -0.225 | -0.098 |
| `log10_lib_atac` | manifest.json (export, tmm) | 74 | 0.232 | -0.040 |
| `log10_lib_rna` | manifest.json (export, tmm) | 74 | 0.281 | 0.018 |
| `log10_n_atac` | covariates.csv | 74 | 0.324* | n/a |
| `log10_n_rna` | covariates.csv | 74 | 0.288 | n/a |
| `male` | metadata.csv | 74 | -0.191 | -0.291 |
| `mitochondrial_mean` | donor_technical.csv | 74 | 0.787* | 0.605* |
| `mitochondrial_pseudobulk` | donor_technical.csv | 74 | 0.766* | 0.573* |
| `mitochondrial_sd` | donor_technical.csv | 74 | 0.683* | 0.491* |
| `mitochondrial_supertype_share` | donor_technical.csv | 74 | 0.157 | 0.262 |
| `n_libraries_rna` | donor_technical.csv | 74 | 0.120 | 0.112 |
| `n_nuclei_rna` | donor_technical.csv | 74 | 0.288 | 0.235 |
| `number_of_umis_mean` | donor_technical.csv | 74 | 0.033 | 0.054 |
| `number_of_umis_sd` | donor_technical.csv | 74 | 0.086 | 0.118 |
| `pmi` | covariates.csv | 74 | 0.101 | n/a |
| `ribosomal_mean` | donor_technical.csv | 74 | -0.840* | -0.797* |
| `ribosomal_pseudobulk` | donor_technical.csv | 74 | -0.815* | -0.769* |
| `ribosomal_sd` | donor_technical.csv | 74 | -0.694* | -0.663* |
| `ribosomal_supertype_share` | donor_technical.csv | 74 | -0.070 | 0.245 |
| `rin` | covariates.csv | 74 | -0.299 | n/a |
| `stage` | metadata.csv | 74 | 0.504* | 0.429* |
| `tmm_atac` | manifest.json (export, tmm) | 74 | 0.557* | 0.581* |
| `tmm_rna` | manifest.json (export, tmm) | 74 | -0.724* | -0.616* |
| `years_education` | covariates.csv | 74 | -0.013 | n/a |

`*` marks |rho| ≥ 0.3.

### Forced-rank analyses (pasted from `diagnostics.md`)

Per component: score SD, the number of variables flagged on the raw correlation, the variable with the largest |rho| and that rho, the share of each sex's net displacement and of the per-stage sex-difference vectors along the component. The full rows are in `axis_correlations.csv` and `geometry_shares.csv`.

| analysis | rank | component | sd | flagged | top variable | rho | share net Female | share net Male | share sex diff | share sex diff centred |
|---|---:|---|---:|---:|---|---:|---:|---:|---:|---:|
| rank_02 | 2 | pls_0 | 35.76 | 20 | `ribosomal_mean` | -0.840 | 0.662 | 0.270 | 0.987 | 0.973 |
| rank_02 | 2 | pls_1 | 18.17 | 3 | `stage` | -0.616 | 0.338 | 0.730 | 0.013 | 0.027 |
| rank_03 | 3 | pls_0 | 35.76 | 20 | `ribosomal_mean` | -0.840 | 0.646 | 0.269 | 0.870 | 0.915 |
| rank_03 | 3 | pls_1 | 18.17 | 3 | `stage` | -0.616 | 0.330 | 0.726 | 0.011 | 0.025 |
| rank_03 | 3 | pls_2 | 17.53 | 1 | `atac_percent_duplicates_mean` | 0.345 | 0.024 | 0.005 | 0.119 | 0.060 |
| rank_04 | 4 | pls_0 | 35.76 | 20 | `ribosomal_mean` | -0.840 | 0.534 | 0.253 | 0.700 | 0.719 |
| rank_04 | 4 | pls_1 | 18.17 | 3 | `stage` | -0.616 | 0.273 | 0.684 | 0.009 | 0.020 |
| rank_04 | 4 | pls_2 | 17.53 | 1 | `atac_percent_duplicates_mean` | 0.345 | 0.020 | 0.005 | 0.096 | 0.047 |
| rank_04 | 4 | pls_3 | 18 | 0 | `gex_median_umi_counts_per_cell_sd` | 0.283 | 0.173 | 0.058 | 0.195 | 0.214 |
| rank_06 | 6 | pls_0 | 35.76 | 20 | `ribosomal_mean` | -0.840 | 0.478 | 0.190 | 0.604 | 0.599 |
| rank_06 | 6 | pls_1 | 18.17 | 3 | `stage` | -0.616 | 0.244 | 0.513 | 0.008 | 0.016 |
| rank_06 | 6 | pls_2 | 17.53 | 1 | `atac_percent_duplicates_mean` | 0.345 | 0.018 | 0.004 | 0.083 | 0.039 |
| rank_06 | 6 | pls_3 | 18 | 0 | `gex_median_umi_counts_per_cell_sd` | 0.283 | 0.155 | 0.043 | 0.169 | 0.178 |
| rank_06 | 6 | pls_4 | 15.74 | 4 | `log10_n_atac` | 0.525 | 0.106 | 0.245 | 0.055 | 0.135 |
| rank_06 | 6 | pls_5 | 13.3 | 3 | `gex_fraction_of_transcriptomic_reads_in_cells_sd` | -0.327 | 0.000 | 0.005 | 0.082 | 0.032 |
| rank_09 | 9 | pls_0 | 35.76 | 20 | `ribosomal_mean` | -0.840 | 0.472 | 0.168 | 0.507 | 0.418 |
| rank_09 | 9 | pls_1 | 18.17 | 3 | `stage` | -0.616 | 0.241 | 0.454 | 0.007 | 0.011 |
| rank_09 | 9 | pls_2 | 17.53 | 1 | `atac_percent_duplicates_mean` | 0.345 | 0.018 | 0.003 | 0.069 | 0.027 |
| rank_09 | 9 | pls_3 | 18 | 0 | `gex_median_umi_counts_per_cell_sd` | 0.283 | 0.153 | 0.038 | 0.142 | 0.124 |
| rank_09 | 9 | pls_4 | 15.74 | 4 | `log10_n_atac` | 0.525 | 0.104 | 0.217 | 0.046 | 0.094 |
| rank_09 | 9 | pls_5 | 13.3 | 3 | `gex_fraction_of_transcriptomic_reads_in_cells_sd` | -0.327 | 0.000 | 0.004 | 0.069 | 0.022 |
| rank_09 | 9 | pls_6 | 10.57 | 0 | `stage` | 0.246 | 0.005 | 0.051 | 0.026 | 0.048 |
| rank_09 | 9 | pls_7 | 12.33 | 0 | `log10_lib_atac` | 0.288 | 0.007 | 0.001 | 0.034 | 0.049 |
| rank_09 | 9 | pls_8 | 13.06 | 0 | `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 0.210 | 0.000 | 0.062 | 0.100 | 0.205 |
| rank_12 | 12 | pls_0 | 35.76 | 20 | `ribosomal_mean` | -0.840 | 0.469 | 0.162 | 0.459 | 0.369 |
| rank_12 | 12 | pls_1 | 18.17 | 3 | `stage` | -0.616 | 0.239 | 0.439 | 0.006 | 0.010 |
| rank_12 | 12 | pls_2 | 17.53 | 1 | `atac_percent_duplicates_mean` | 0.345 | 0.017 | 0.003 | 0.063 | 0.024 |
| rank_12 | 12 | pls_3 | 18 | 0 | `gex_median_umi_counts_per_cell_sd` | 0.283 | 0.152 | 0.037 | 0.128 | 0.110 |
| rank_12 | 12 | pls_4 | 15.74 | 4 | `log10_n_atac` | 0.525 | 0.104 | 0.210 | 0.042 | 0.083 |
| rank_12 | 12 | pls_5 | 13.3 | 3 | `gex_fraction_of_transcriptomic_reads_in_cells_sd` | -0.327 | 0.000 | 0.004 | 0.062 | 0.020 |
| rank_12 | 12 | pls_6 | 10.57 | 0 | `stage` | 0.246 | 0.005 | 0.050 | 0.024 | 0.042 |
| rank_12 | 12 | pls_7 | 12.33 | 0 | `log10_lib_atac` | 0.288 | 0.007 | 0.001 | 0.031 | 0.044 |
| rank_12 | 12 | pls_8 | 13.06 | 0 | `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 0.210 | 0.000 | 0.060 | 0.090 | 0.181 |
| rank_12 | 12 | pls_9 | 10.88 | 0 | `atac_percent_duplicates_mean` | -0.242 | 0.000 | 0.021 | 0.022 | 0.037 |
| rank_12 | 12 | pls_10 | 12.63 | 0 | `gex_reads_mapped_confidently_to_intronic_regions_sd` | 0.234 | 0.005 | 0.002 | 0.026 | 0.041 |
| rank_12 | 12 | pls_11 | 11.36 | 0 | `atac_median_high_quality_fragments_per_cell_mean` | 0.214 | 0.002 | 0.011 | 0.047 | 0.040 |

### Block principal components (pasted from `diagnostics.md`)

#### atac (variance share: pc1 0.150, pc2 0.123, pc3 0.034)

| variable | n | pc1 rho | pc1 rho \| cov | pc2 rho | pc2 rho \| cov | pc3 rho | pc3 rho \| cov |
|---|---:|---:|---:|---:|---:|---:|---:|
| `age_at_death` | 74 | -0.293 | n/a | -0.137 | n/a | 0.092 | n/a |
| `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | -0.503* | -0.348* | 0.449* | 0.333* | 0.374* | 0.258 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | -0.498* | -0.393* | 0.620* | 0.518* | 0.287 | 0.292 |
| `atac_median_high_quality_fragments_per_cell_mean` | 74 | -0.214 | 0.121 | 0.128 | -0.029 | 0.022 | 0.046 |
| `atac_percent_duplicates_mean` | 74 | -0.164 | -0.247 | 0.334* | 0.386* | 0.206 | 0.063 |
| `atac_tss_enrichment_score_mean` | 74 | -0.436* | -0.441* | 0.598* | 0.505* | 0.354* | 0.284 |
| `doublet_score_mean` | 74 | 0.096 | 0.111 | 0.003 | -0.115 | -0.045 | -0.046 |
| `doublet_score_sd` | 74 | 0.142 | 0.061 | -0.331* | -0.290 | -0.233 | -0.185 |
| `fraction_mitochondrial_umis_mean` | 74 | 0.420* | 0.246 | -0.598* | -0.502* | -0.279 | -0.147 |
| `fraction_mitochondrial_umis_sd` | 74 | 0.433* | 0.228 | -0.556* | -0.500* | -0.266 | -0.080 |
| `genes_detected_mean` | 74 | -0.192 | -0.151 | -0.002 | -0.065 | 0.254 | 0.177 |
| `genes_detected_sd` | 74 | -0.135 | -0.124 | -0.051 | -0.045 | 0.189 | 0.124 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | -0.243 | -0.167 | 0.448* | 0.370* | 0.354* | 0.093 |
| `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 74 | 0.433* | 0.167 | -0.143 | -0.094 | 0.068 | -0.045 |
| `gex_median_umi_counts_per_cell_mean` | 74 | -0.214 | -0.211 | 0.475* | 0.391* | 0.256 | 0.102 |
| `gex_median_umi_counts_per_cell_sd` | 74 | 0.179 | 0.119 | 0.112 | 0.183 | 0.153 | 0.003 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | -0.236 | -0.325* | 0.316* | 0.254 | 0.504* | 0.257 |
| `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | 0.512* | 0.227 | -0.306* | -0.114 | 0.082 | -0.010 |
| `gex_sequencing_saturation_mean` | 74 | 0.441* | 0.302* | -0.540* | -0.392* | -0.240 | -0.179 |
| `gex_sequencing_saturation_sd` | 74 | 0.021 | -0.016 | 0.099 | 0.172 | 0.018 | -0.114 |
| `hemoglobin_mean` | 74 | -0.035 | -0.080 | -0.001 | 0.072 | -0.092 | -0.032 |
| `hemoglobin_pseudobulk` | 74 | -0.104 | -0.108 | 0.041 | 0.047 | -0.130 | -0.087 |
| `hemoglobin_sd` | 74 | 0.005 | -0.050 | -0.018 | 0.116 | -0.085 | -0.128 |
| `hemoglobin_supertype_share` | 69 | -0.206 | 0.010 | 0.063 | 0.165 | -0.088 | -0.020 |
| `log10_lib_atac` | 74 | 0.546* | 0.113 | -0.218 | 0.081 | 0.415* | 0.298 |
| `log10_lib_rna` | 74 | 0.170 | -0.068 | -0.288 | -0.087 | 0.177 | 0.114 |
| `log10_n_atac` | 74 | 0.607* | n/a | -0.313* | n/a | 0.301* | n/a |
| `log10_n_rna` | 74 | 0.218 | n/a | -0.271 | n/a | 0.108 | n/a |
| `male` | 74 | 0.112 | -0.016 | 0.226 | 0.307* | 0.268 | 0.228 |
| `mitochondrial_mean` | 74 | 0.420* | 0.246 | -0.598* | -0.502* | -0.279 | -0.147 |
| `mitochondrial_pseudobulk` | 74 | 0.398* | 0.247 | -0.593* | -0.473* | -0.277 | -0.143 |
| `mitochondrial_sd` | 74 | 0.433* | 0.228 | -0.556* | -0.500* | -0.266 | -0.080 |
| `mitochondrial_supertype_share` | 74 | 0.138 | 0.226 | 0.115 | 0.071 | -0.153 | -0.160 |
| `n_libraries_rna` | 74 | 0.324* | 0.192 | -0.141 | -0.004 | 0.203 | -0.086 |
| `n_nuclei_rna` | 74 | 0.218 | 0.084 | -0.271 | -0.186 | 0.108 | -0.078 |
| `number_of_umis_mean` | 74 | -0.019 | -0.054 | -0.108 | -0.110 | 0.247 | 0.115 |
| `number_of_umis_sd` | 74 | -0.010 | -0.014 | -0.108 | -0.085 | 0.235 | 0.115 |
| `pmi` | 74 | 0.113 | n/a | -0.055 | n/a | -0.102 | n/a |
| `ribosomal_mean` | 74 | -0.416* | -0.321* | 0.627* | 0.514* | 0.214 | 0.273 |
| `ribosomal_pseudobulk` | 74 | -0.395* | -0.305* | 0.616* | 0.518* | 0.208 | 0.239 |
| `ribosomal_sd` | 74 | -0.327* | -0.269 | 0.504* | 0.390* | 0.237 | 0.246 |
| `ribosomal_supertype_share` | 74 | 0.036 | 0.265 | 0.284 | 0.016 | -0.038 | -0.132 |
| `rin` | 74 | -0.009 | n/a | 0.221 | n/a | 0.403* | n/a |
| `stage` | 74 | 0.134 | 0.048 | -0.270 | -0.225 | -0.360* | -0.337* |
| `tmm_atac` | 74 | 0.195 | 0.335* | -0.425* | -0.443* | -0.406* | -0.172 |
| `tmm_rna` | 74 | -0.611* | -0.442* | 0.487* | 0.338* | 0.086 | 0.192 |
| `years_education` | 74 | -0.006 | n/a | 0.104 | n/a | -0.181 | n/a |

#### rna (variance share: pc1 0.282, pc2 0.085, pc3 0.060)

| variable | n | pc1 rho | pc1 rho \| cov | pc2 rho | pc2 rho \| cov | pc3 rho | pc3 rho \| cov |
|---|---:|---:|---:|---:|---:|---:|---:|
| `age_at_death` | 74 | 0.144 | n/a | 0.174 | n/a | 0.082 | n/a |
| `atac_fraction_of_high_quality_fragments_in_cells_mean` | 74 | 0.591* | 0.397* | -0.185 | -0.074 | 0.135 | 0.070 |
| `atac_fraction_of_high_quality_fragments_overlapping_peaks_mean` | 74 | 0.720* | 0.631* | -0.193 | -0.115 | 0.143 | 0.072 |
| `atac_median_high_quality_fragments_per_cell_mean` | 74 | 0.022 | -0.191 | -0.046 | 0.032 | -0.078 | 0.122 |
| `atac_percent_duplicates_mean` | 74 | 0.297 | 0.306* | -0.106 | -0.145 | 0.214 | 0.040 |
| `atac_tss_enrichment_score_mean` | 74 | 0.729* | 0.667* | -0.245 | -0.159 | 0.223 | 0.085 |
| `doublet_score_mean` | 74 | -0.096 | -0.157 | -0.081 | -0.033 | 0.230 | 0.133 |
| `doublet_score_sd` | 74 | -0.345* | -0.252 | 0.146 | 0.129 | 0.018 | 0.040 |
| `fraction_mitochondrial_umis_mean` | 74 | -0.762* | -0.549* | 0.348* | 0.271 | -0.129 | 0.017 |
| `fraction_mitochondrial_umis_sd` | 74 | -0.686* | -0.466* | 0.342* | 0.317* | -0.129 | 0.037 |
| `genes_detected_mean` | 74 | 0.117 | 0.061 | 0.216 | 0.245 | 0.116 | 0.126 |
| `genes_detected_sd` | 74 | -0.064 | -0.094 | 0.272 | 0.254 | 0.083 | 0.131 |
| `gex_fraction_of_transcriptomic_reads_in_cells_mean` | 74 | 0.602* | 0.390* | -0.149 | -0.150 | 0.165 | -0.070 |
| `gex_fraction_of_transcriptomic_reads_in_cells_sd` | 74 | -0.248 | -0.123 | 0.009 | -0.078 | 0.199 | 0.098 |
| `gex_median_umi_counts_per_cell_mean` | 74 | 0.485* | 0.369* | -0.153 | -0.101 | 0.147 | -0.004 |
| `gex_median_umi_counts_per_cell_sd` | 74 | 0.075 | 0.031 | -0.018 | -0.185 | 0.069 | -0.037 |
| `gex_reads_mapped_confidently_to_intronic_regions_mean` | 74 | 0.550* | 0.520* | -0.054 | -0.053 | 0.338* | 0.081 |
| `gex_reads_mapped_confidently_to_intronic_regions_sd` | 74 | -0.299 | -0.090 | 0.141 | -0.013 | -0.052 | -0.024 |
| `gex_sequencing_saturation_mean` | 74 | -0.686* | -0.515* | 0.189 | 0.134 | -0.084 | -0.002 |
| `gex_sequencing_saturation_sd` | 74 | 0.096 | 0.071 | 0.069 | -0.063 | -0.013 | 0.002 |
| `hemoglobin_mean` | 74 | 0.087 | 0.075 | -0.019 | -0.017 | -0.046 | -0.050 |
| `hemoglobin_pseudobulk` | 74 | 0.143 | 0.054 | -0.025 | 0.002 | -0.085 | -0.089 |
| `hemoglobin_sd` | 74 | 0.052 | -0.011 | -0.033 | -0.098 | -0.053 | -0.052 |
| `hemoglobin_supertype_share` | 69 | 0.211 | 0.073 | -0.077 | -0.126 | 0.055 | 0.019 |
| `log10_lib_atac` | 74 | -0.226 | -0.006 | 0.189 | 0.016 | 0.260 | 0.221 |
| `log10_lib_rna` | 74 | -0.232 | -0.042 | 0.286 | 0.210 | -0.176 | 0.127 |
| `log10_n_atac` | 74 | -0.297 | n/a | 0.225 | n/a | 0.252 | n/a |
| `log10_n_rna` | 74 | -0.232 | n/a | 0.227 | n/a | -0.262 | n/a |
| `male` | 74 | 0.177 | 0.237 | -0.090 | -0.133 | -0.115 | -0.058 |
| `mitochondrial_mean` | 74 | -0.762* | -0.549* | 0.348* | 0.271 | -0.129 | 0.017 |
| `mitochondrial_pseudobulk` | 74 | -0.736* | -0.517* | 0.352* | 0.235 | -0.123 | 0.002 |
| `mitochondrial_sd` | 74 | -0.686* | -0.466* | 0.342* | 0.317* | -0.129 | 0.037 |
| `mitochondrial_supertype_share` | 74 | -0.207 | -0.310* | -0.052 | -0.037 | 0.065 | -0.058 |
| `n_libraries_rna` | 74 | -0.110 | -0.112 | 0.051 | -0.193 | 0.127 | 0.083 |
| `n_nuclei_rna` | 74 | -0.232 | -0.187 | 0.227 | -0.065 | -0.262 | 0.152 |
| `number_of_umis_mean` | 74 | -0.042 | -0.077 | 0.238 | 0.205 | 0.119 | 0.152 |
| `number_of_umis_sd` | 74 | -0.124 | -0.190 | 0.255 | 0.188 | 0.147 | 0.196 |
| `pmi` | 74 | -0.106 | n/a | 0.088 | n/a | -0.169 | n/a |
| `ribosomal_mean` | 74 | 0.816* | 0.775* | -0.405* | -0.259 | -0.029 | -0.064 |
| `ribosomal_pseudobulk` | 74 | 0.791* | 0.747* | -0.410* | -0.281 | -0.039 | -0.103 |
| `ribosomal_sd` | 74 | 0.698* | 0.653* | -0.411* | -0.216 | 0.018 | -0.031 |
| `ribosomal_supertype_share` | 74 | 0.035 | -0.246 | -0.346* | -0.177 | 0.044 | -0.051 |
| `rin` | 74 | 0.247 | n/a | -0.133 | n/a | 0.361* | n/a |
| `stage` | 74 | -0.309* | -0.215 | -0.026 | -0.069 | -0.191 | -0.216 |
| `tmm_atac` | 74 | -0.545* | -0.550* | 0.077 | 0.097 | -0.095 | 0.058 |
| `tmm_rna` | 74 | 0.731* | 0.616* | -0.159 | -0.020 | 0.019 | -0.048 |
| `years_education` | 74 | -0.028 | n/a | -0.098 | n/a | -0.114 | n/a |

### Geometry shares (pasted from `diagnostics.md`)

Share of the squared length of each vector of the trajectory geometry (from `stage_means.csv`) that lies along each component; `value` is the signed component of a single vector (the net displacements), empty for the summed quantities. At rank 1 every share is 1 by construction.

| analysis | rank | quantity | component | value | along | total | share |
|---|---:|---|---|---:|---:|---:|---:|
| primary | 1 | net_displacement_Female | pls_0 | 37.15 | 1380 | 1380 | 1.000 |
| primary | 1 | net_displacement_Male | pls_0 | 18.44 | 340.1 | 340.1 | 1.000 |
| primary | 1 | sex_difference_by_stage | pls_0 | n/a | 691.6 | 691.6 | 1.000 |
| primary | 1 | sex_difference_centred | pls_0 | n/a | 276.3 | 276.3 | 1.000 |
| rank_02 | 2 | net_displacement_Female | pls_0 | 37.15 | 1380 | 2084 | 0.662 |
| rank_02 | 2 | net_displacement_Female | pls_1 | -26.53 | 703.9 | 2084 | 0.338 |
| rank_02 | 2 | net_displacement_Male | pls_0 | 18.44 | 340.1 | 1260 | 0.270 |
| rank_02 | 2 | net_displacement_Male | pls_1 | -30.32 | 919.5 | 1260 | 0.730 |
| rank_02 | 2 | sex_difference_by_stage | pls_0 | n/a | 691.6 | 700.5 | 0.987 |
| rank_02 | 2 | sex_difference_by_stage | pls_1 | n/a | 8.935 | 700.5 | 0.013 |
| rank_02 | 2 | sex_difference_centred | pls_0 | n/a | 276.3 | 283.9 | 0.973 |
| rank_02 | 2 | sex_difference_centred | pls_1 | n/a | 7.542 | 283.9 | 0.027 |
| rank_03 | 3 | net_displacement_Female | pls_0 | 37.15 | 1380 | 2135 | 0.646 |
| rank_03 | 3 | net_displacement_Female | pls_1 | -26.53 | 703.9 | 2135 | 0.330 |
| rank_03 | 3 | net_displacement_Female | pls_2 | 7.165 | 51.33 | 2135 | 0.024 |
| rank_03 | 3 | net_displacement_Male | pls_0 | 18.44 | 340.1 | 1266 | 0.269 |
| rank_03 | 3 | net_displacement_Male | pls_1 | -30.32 | 919.5 | 1266 | 0.726 |
| rank_03 | 3 | net_displacement_Male | pls_2 | 2.621 | 6.867 | 1266 | 0.005 |
| rank_03 | 3 | sex_difference_by_stage | pls_0 | n/a | 691.6 | 795.3 | 0.870 |
| rank_03 | 3 | sex_difference_by_stage | pls_1 | n/a | 8.935 | 795.3 | 0.011 |
| rank_03 | 3 | sex_difference_by_stage | pls_2 | n/a | 94.76 | 795.3 | 0.119 |
| rank_03 | 3 | sex_difference_centred | pls_0 | n/a | 276.3 | 301.9 | 0.915 |
| rank_03 | 3 | sex_difference_centred | pls_1 | n/a | 7.542 | 301.9 | 0.025 |
| rank_03 | 3 | sex_difference_centred | pls_2 | n/a | 18.08 | 301.9 | 0.060 |
| rank_04 | 4 | net_displacement_Female | pls_0 | 37.15 | 1380 | 2582 | 0.534 |
| rank_04 | 4 | net_displacement_Female | pls_1 | -26.53 | 703.9 | 2582 | 0.273 |
| rank_04 | 4 | net_displacement_Female | pls_2 | 7.165 | 51.33 | 2582 | 0.020 |
| rank_04 | 4 | net_displacement_Female | pls_3 | 21.14 | 447 | 2582 | 0.173 |
| rank_04 | 4 | net_displacement_Male | pls_0 | 18.44 | 340.1 | 1344 | 0.253 |
| rank_04 | 4 | net_displacement_Male | pls_1 | -30.32 | 919.5 | 1344 | 0.684 |
| rank_04 | 4 | net_displacement_Male | pls_2 | 2.621 | 6.867 | 1344 | 0.005 |
| rank_04 | 4 | net_displacement_Male | pls_3 | 8.804 | 77.51 | 1344 | 0.058 |
| rank_04 | 4 | sex_difference_by_stage | pls_0 | n/a | 691.6 | 988.5 | 0.700 |
| rank_04 | 4 | sex_difference_by_stage | pls_1 | n/a | 8.935 | 988.5 | 0.009 |
| rank_04 | 4 | sex_difference_by_stage | pls_2 | n/a | 94.76 | 988.5 | 0.096 |
| rank_04 | 4 | sex_difference_by_stage | pls_3 | n/a | 193.2 | 988.5 | 0.195 |
| rank_04 | 4 | sex_difference_centred | pls_0 | n/a | 276.3 | 384.1 | 0.719 |
| rank_04 | 4 | sex_difference_centred | pls_1 | n/a | 7.542 | 384.1 | 0.020 |
| rank_04 | 4 | sex_difference_centred | pls_2 | n/a | 18.08 | 384.1 | 0.047 |
| rank_04 | 4 | sex_difference_centred | pls_3 | n/a | 82.13 | 384.1 | 0.214 |
| rank_06 | 6 | net_displacement_Female | pls_0 | 37.15 | 1380 | 2889 | 0.478 |
| rank_06 | 6 | net_displacement_Female | pls_1 | -26.53 | 703.9 | 2889 | 0.244 |
| rank_06 | 6 | net_displacement_Female | pls_2 | 7.165 | 51.33 | 2889 | 0.018 |
| rank_06 | 6 | net_displacement_Female | pls_3 | 21.14 | 447 | 2889 | 0.155 |
| rank_06 | 6 | net_displacement_Female | pls_4 | 17.46 | 305 | 2889 | 0.106 |
| rank_06 | 6 | net_displacement_Female | pls_5 | -1.165 | 1.358 | 2889 | 0.000 |
| rank_06 | 6 | net_displacement_Male | pls_0 | 18.44 | 340.1 | 1792 | 0.190 |
| rank_06 | 6 | net_displacement_Male | pls_1 | -30.32 | 919.5 | 1792 | 0.513 |
| rank_06 | 6 | net_displacement_Male | pls_2 | 2.621 | 6.867 | 1792 | 0.004 |
| rank_06 | 6 | net_displacement_Male | pls_3 | 8.804 | 77.51 | 1792 | 0.043 |
| rank_06 | 6 | net_displacement_Male | pls_4 | 20.95 | 439 | 1792 | 0.245 |
| rank_06 | 6 | net_displacement_Male | pls_5 | 2.967 | 8.801 | 1792 | 0.005 |
| rank_06 | 6 | sex_difference_by_stage | pls_0 | n/a | 691.6 | 1145 | 0.604 |
| rank_06 | 6 | sex_difference_by_stage | pls_1 | n/a | 8.935 | 1145 | 0.008 |
| rank_06 | 6 | sex_difference_by_stage | pls_2 | n/a | 94.76 | 1145 | 0.083 |
| rank_06 | 6 | sex_difference_by_stage | pls_3 | n/a | 193.2 | 1145 | 0.169 |
| rank_06 | 6 | sex_difference_by_stage | pls_4 | n/a | 63.26 | 1145 | 0.055 |
| rank_06 | 6 | sex_difference_by_stage | pls_5 | n/a | 93.62 | 1145 | 0.082 |
| rank_06 | 6 | sex_difference_centred | pls_0 | n/a | 276.3 | 461.4 | 0.599 |
| rank_06 | 6 | sex_difference_centred | pls_1 | n/a | 7.542 | 461.4 | 0.016 |
| rank_06 | 6 | sex_difference_centred | pls_2 | n/a | 18.08 | 461.4 | 0.039 |
| rank_06 | 6 | sex_difference_centred | pls_3 | n/a | 82.13 | 461.4 | 0.178 |
| rank_06 | 6 | sex_difference_centred | pls_4 | n/a | 62.42 | 461.4 | 0.135 |
| rank_06 | 6 | sex_difference_centred | pls_5 | n/a | 14.86 | 461.4 | 0.032 |
| rank_09 | 9 | net_displacement_Female | pls_0 | 37.15 | 1380 | 2922 | 0.472 |
| rank_09 | 9 | net_displacement_Female | pls_1 | -26.53 | 703.9 | 2922 | 0.241 |
| rank_09 | 9 | net_displacement_Female | pls_2 | 7.165 | 51.33 | 2922 | 0.018 |
| rank_09 | 9 | net_displacement_Female | pls_3 | 21.14 | 447 | 2922 | 0.153 |
| rank_09 | 9 | net_displacement_Female | pls_4 | 17.46 | 305 | 2922 | 0.104 |
| rank_09 | 9 | net_displacement_Female | pls_5 | -1.165 | 1.358 | 2922 | 0.000 |
| rank_09 | 9 | net_displacement_Female | pls_6 | 3.641 | 13.25 | 2922 | 0.005 |
| rank_09 | 9 | net_displacement_Female | pls_7 | -4.394 | 19.31 | 2922 | 0.007 |
| rank_09 | 9 | net_displacement_Female | pls_8 | -0.4361 | 0.1902 | 2922 | 0.000 |
| rank_09 | 9 | net_displacement_Male | pls_0 | 18.44 | 340.1 | 2024 | 0.168 |
| rank_09 | 9 | net_displacement_Male | pls_1 | -30.32 | 919.5 | 2024 | 0.454 |
| rank_09 | 9 | net_displacement_Male | pls_2 | 2.621 | 6.867 | 2024 | 0.003 |
| rank_09 | 9 | net_displacement_Male | pls_3 | 8.804 | 77.51 | 2024 | 0.038 |
| rank_09 | 9 | net_displacement_Male | pls_4 | 20.95 | 439 | 2024 | 0.217 |
| rank_09 | 9 | net_displacement_Male | pls_5 | 2.967 | 8.801 | 2024 | 0.004 |
| rank_09 | 9 | net_displacement_Male | pls_6 | 10.2 | 104.1 | 2024 | 0.051 |
| rank_09 | 9 | net_displacement_Male | pls_7 | -1.324 | 1.753 | 2024 | 0.001 |
| rank_09 | 9 | net_displacement_Male | pls_8 | 11.23 | 126.2 | 2024 | 0.062 |
| rank_09 | 9 | sex_difference_by_stage | pls_0 | n/a | 691.6 | 1364 | 0.507 |
| rank_09 | 9 | sex_difference_by_stage | pls_1 | n/a | 8.935 | 1364 | 0.007 |
| rank_09 | 9 | sex_difference_by_stage | pls_2 | n/a | 94.76 | 1364 | 0.069 |
| rank_09 | 9 | sex_difference_by_stage | pls_3 | n/a | 193.2 | 1364 | 0.142 |
| rank_09 | 9 | sex_difference_by_stage | pls_4 | n/a | 63.26 | 1364 | 0.046 |
| rank_09 | 9 | sex_difference_by_stage | pls_5 | n/a | 93.62 | 1364 | 0.069 |
| rank_09 | 9 | sex_difference_by_stage | pls_6 | n/a | 36.08 | 1364 | 0.026 |
| rank_09 | 9 | sex_difference_by_stage | pls_7 | n/a | 46.14 | 1364 | 0.034 |
| rank_09 | 9 | sex_difference_by_stage | pls_8 | n/a | 136.1 | 1364 | 0.100 |
| rank_09 | 9 | sex_difference_centred | pls_0 | n/a | 276.3 | 661 | 0.418 |
| rank_09 | 9 | sex_difference_centred | pls_1 | n/a | 7.542 | 661 | 0.011 |
| rank_09 | 9 | sex_difference_centred | pls_2 | n/a | 18.08 | 661 | 0.027 |
| rank_09 | 9 | sex_difference_centred | pls_3 | n/a | 82.13 | 661 | 0.124 |
| rank_09 | 9 | sex_difference_centred | pls_4 | n/a | 62.42 | 661 | 0.094 |
| rank_09 | 9 | sex_difference_centred | pls_5 | n/a | 14.86 | 661 | 0.022 |
| rank_09 | 9 | sex_difference_centred | pls_6 | n/a | 31.63 | 661 | 0.048 |
| rank_09 | 9 | sex_difference_centred | pls_7 | n/a | 32.64 | 661 | 0.049 |
| rank_09 | 9 | sex_difference_centred | pls_8 | n/a | 135.4 | 661 | 0.205 |
| rank_12 | 12 | net_displacement_Female | pls_0 | 37.15 | 1380 | 2944 | 0.469 |
| rank_12 | 12 | net_displacement_Female | pls_1 | -26.53 | 703.9 | 2944 | 0.239 |
| rank_12 | 12 | net_displacement_Female | pls_2 | 7.165 | 51.33 | 2944 | 0.017 |
| rank_12 | 12 | net_displacement_Female | pls_3 | 21.14 | 447 | 2944 | 0.152 |
| rank_12 | 12 | net_displacement_Female | pls_4 | 17.46 | 305 | 2944 | 0.104 |
| rank_12 | 12 | net_displacement_Female | pls_5 | -1.165 | 1.358 | 2944 | 0.000 |
| rank_12 | 12 | net_displacement_Female | pls_6 | 3.641 | 13.25 | 2944 | 0.005 |
| rank_12 | 12 | net_displacement_Female | pls_7 | -4.394 | 19.31 | 2944 | 0.007 |
| rank_12 | 12 | net_displacement_Female | pls_8 | -0.4361 | 0.1902 | 2944 | 0.000 |
| rank_12 | 12 | net_displacement_Female | pls_9 | -0.8295 | 0.6881 | 2944 | 0.000 |
| rank_12 | 12 | net_displacement_Female | pls_10 | 4.012 | 16.09 | 2944 | 0.005 |
| rank_12 | 12 | net_displacement_Female | pls_11 | 2.46 | 6.052 | 2944 | 0.002 |
| rank_12 | 12 | net_displacement_Male | pls_0 | 18.44 | 340.1 | 2095 | 0.162 |
| rank_12 | 12 | net_displacement_Male | pls_1 | -30.32 | 919.5 | 2095 | 0.439 |
| rank_12 | 12 | net_displacement_Male | pls_2 | 2.621 | 6.867 | 2095 | 0.003 |
| rank_12 | 12 | net_displacement_Male | pls_3 | 8.804 | 77.51 | 2095 | 0.037 |
| rank_12 | 12 | net_displacement_Male | pls_4 | 20.95 | 439 | 2095 | 0.210 |
| rank_12 | 12 | net_displacement_Male | pls_5 | 2.967 | 8.801 | 2095 | 0.004 |
| rank_12 | 12 | net_displacement_Male | pls_6 | 10.2 | 104.1 | 2095 | 0.050 |
| rank_12 | 12 | net_displacement_Male | pls_7 | -1.324 | 1.753 | 2095 | 0.001 |
| rank_12 | 12 | net_displacement_Male | pls_8 | 11.23 | 126.2 | 2095 | 0.060 |
| rank_12 | 12 | net_displacement_Male | pls_9 | 6.577 | 43.26 | 2095 | 0.021 |
| rank_12 | 12 | net_displacement_Male | pls_10 | -2.237 | 5.003 | 2095 | 0.002 |
| rank_12 | 12 | net_displacement_Male | pls_11 | 4.771 | 22.76 | 2095 | 0.011 |
| rank_12 | 12 | sex_difference_by_stage | pls_0 | n/a | 691.6 | 1506 | 0.459 |
| rank_12 | 12 | sex_difference_by_stage | pls_1 | n/a | 8.935 | 1506 | 0.006 |
| rank_12 | 12 | sex_difference_by_stage | pls_2 | n/a | 94.76 | 1506 | 0.063 |
| rank_12 | 12 | sex_difference_by_stage | pls_3 | n/a | 193.2 | 1506 | 0.128 |
| rank_12 | 12 | sex_difference_by_stage | pls_4 | n/a | 63.26 | 1506 | 0.042 |
| rank_12 | 12 | sex_difference_by_stage | pls_5 | n/a | 93.62 | 1506 | 0.062 |
| rank_12 | 12 | sex_difference_by_stage | pls_6 | n/a | 36.08 | 1506 | 0.024 |
| rank_12 | 12 | sex_difference_by_stage | pls_7 | n/a | 46.14 | 1506 | 0.031 |
| rank_12 | 12 | sex_difference_by_stage | pls_8 | n/a | 136.1 | 1506 | 0.090 |
| rank_12 | 12 | sex_difference_by_stage | pls_9 | n/a | 33.25 | 1506 | 0.022 |
| rank_12 | 12 | sex_difference_by_stage | pls_10 | n/a | 38.96 | 1506 | 0.026 |
| rank_12 | 12 | sex_difference_by_stage | pls_11 | n/a | 70.23 | 1506 | 0.047 |
| rank_12 | 12 | sex_difference_centred | pls_0 | n/a | 276.3 | 749 | 0.369 |
| rank_12 | 12 | sex_difference_centred | pls_1 | n/a | 7.542 | 749 | 0.010 |
| rank_12 | 12 | sex_difference_centred | pls_2 | n/a | 18.08 | 749 | 0.024 |
| rank_12 | 12 | sex_difference_centred | pls_3 | n/a | 82.13 | 749 | 0.110 |
| rank_12 | 12 | sex_difference_centred | pls_4 | n/a | 62.42 | 749 | 0.083 |
| rank_12 | 12 | sex_difference_centred | pls_5 | n/a | 14.86 | 749 | 0.020 |
| rank_12 | 12 | sex_difference_centred | pls_6 | n/a | 31.63 | 749 | 0.042 |
| rank_12 | 12 | sex_difference_centred | pls_7 | n/a | 32.64 | 749 | 0.044 |
| rank_12 | 12 | sex_difference_centred | pls_8 | n/a | 135.4 | 749 | 0.181 |
| rank_12 | 12 | sex_difference_centred | pls_9 | n/a | 27.55 | 749 | 0.037 |
| rank_12 | 12 | sex_difference_centred | pls_10 | n/a | 30.54 | 749 | 0.041 |
| rank_12 | 12 | sex_difference_centred | pls_11 | n/a | 29.87 | 749 | 0.040 |

### Feature selection (pasted from `diagnostics.md`)

Per omic, from the export's `features_<omic>.csv` and block: input features, features passing the abundance filter, kept features; quantiles of the kept features' mean log2-CPM; the share of kept features with mean below 0 and below 1; the share of the abundance-passing variance held by the kept features; the variance at the cut (kept minimum, dropped median); the share of kept features whose |Spearman| with the primary analysis's `pls_0` is ≥ 0.3.

| omic | input | abundance-passing | kept | mean q05 | mean q25 | mean q50 | mean q75 | mean q95 | share mean < 0 | share mean < 1 | variance share kept | kept min variance | dropped median variance | share \|rho(pls_0)\| ≥ thr |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rna | 36601 | 20238 | 5000 | -0.662 | -0.096 | 0.818 | 2.229 | 4.957 | 0.282 | 0.538 | 0.574 | 0.6341 | 0.2066 | 0.629 (3145) |
| atac | 218882 | 124818 | 5000 | -0.432 | -0.057 | 0.271 | 0.698 | 1.479 | 0.294 | 0.858 | 0.082 | 1.156 | 0.6894 | 0.462 (2308) |

### Analyses read (pasted from `diagnostics.md`)

| analysis | selection | rank | directory |
|---|---|---:|---|
| primary | cv | 1 | `results/motco_analysis/cohort_50` |
| rank_02 | forced | 2 | `results/motco_rank_sensitivity/cohort_50/rank_02` |
| rank_03 | forced | 3 | `results/motco_rank_sensitivity/cohort_50/rank_03` |
| rank_04 | forced | 4 | `results/motco_rank_sensitivity/cohort_50/rank_04` |
| rank_06 | forced | 6 | `results/motco_rank_sensitivity/cohort_50/rank_06` |
| rank_09 | forced | 9 | `results/motco_rank_sensitivity/cohort_50/rank_09` |
| rank_12 | forced | 12 | `results/motco_rank_sensitivity/cohort_50/rank_12` |

### Findings (pasted from `diagnostics.md`)

- variables with missing donors (pairwise-complete correlations): {'hemoglobin_supertype_share': 5}
- pls_0 is identical (np.allclose) across the primary and every forced rank ['rank_02', 'rank_03', 'rank_04', 'rank_06', 'rank_09', 'rank_12']: PLS components nest, so every forced rank keeps the primary's first axis
- primary pls_0: 21 of 47 variables flagged at |rho| ≥ 0.3 (raw or partial): ribosomal_mean (-0.84; partial -0.80), ribosomal_pseudobulk (-0.82; partial -0.77), fraction_mitochondrial_umis_mean (0.79; partial 0.60), mitochondrial_mean (0.79; partial 0.60), atac_fraction_of_high_quality_fragments_overlapping_peaks_mean (-0.78; partial -0.68), atac_tss_enrichment_score_mean (-0.77; partial -0.69), mitochondrial_pseudobulk (0.77; partial 0.57), gex_sequencing_saturation_mean (0.73; partial 0.56), tmm_rna (-0.72; partial -0.62), ribosomal_sd (-0.69; partial -0.66), fraction_mitochondrial_umis_sd (0.68; partial 0.49), mitochondrial_sd (0.68; partial 0.49), gex_fraction_of_transcriptomic_reads_in_cells_mean (-0.61; partial -0.40), atac_fraction_of_high_quality_fragments_in_cells_mean (-0.60; partial -0.44), gex_reads_mapped_confidently_to_intronic_regions_mean (-0.56; partial -0.54), tmm_atac (0.56; partial 0.58), gex_median_umi_counts_per_cell_mean (-0.51; partial -0.42), stage (0.50; partial 0.43), gex_reads_mapped_confidently_to_intronic_regions_sd (0.33; partial 0.09), log10_n_atac (0.32; partial n/a), atac_percent_duplicates_mean (-0.30; partial -0.33)


## 5. Reading

- **A measurement-space finding, not a new result.** The tables say what the latent axis of the primary analysis correlates
  with across donors; they do not re-read its statistics. The primary analysis stays primary: `delta` was not rejected, `angle`
  is 0° at rank 1 by construction and uninformative, and the `shape` non-rejection is uninformative, exactly as in
  `reports/motco-astrocyte-2026-10-19.md`. A `delta` rejection would have been a difference, not specifically of size.
- **What the axis is.** In both cohorts `pls_0` is a donor-level composition and library-quality axis: the ribosomal-protein
  fraction of the astrocyte nuclei (-0.84 / -0.84), the mitochondrial fraction, the ATAC library's TSS
  enrichment and fraction of fragments in peaks and the RNA TMM factor all load on it more strongly than stage (0.49 /
  0.50), and the six analysis covariates leave those correlations largely in place. The first principal
  component of the exported RNA block tracks the same variables, so the axis is a property of the exported blocks that the
  stage-supervised PLS fit then takes first; the forced-rank sweep inherited it at every rank and therefore tested
  dimensionality, not feature content.
- **Stage is on the same axis.** `pls_0` also correlates 0.49 with stage (partial 0.41) in the primary cohort, so the
  technical variables and stage are not independent in these donors. Residualizing a technical variable from the blocks can
  therefore remove stage-related signal as well as preparation noise, and this report cannot say how much of each. That is a reason to run a
  residualized variant as predeclared sensitivity, beside the primary analysis, not as a correction.
- **Residualization candidates for the next change** (`mtg-astrocyte-feature-selection-sensitivity`), under the rule of
  section 2: the ribosomal-protein UMI fraction and the mitochondrial UMI fraction of the astrocyte nuclei (donor-level
  preparation properties with small within-donor supertype shares). The hemoglobin fraction is not a candidate. The RNA TMM
  factor and the ATAC library metrics are donor-level technical variables correlated with the same axis and are in the
  technical table for that change to consider; its adoption rule is written in its proposal before any run, and no variant is
  chosen here.
- **No downstream read.** The diagnostic reads stage and sex. Its outputs are an input of no `export`, `analyze`,
  `rank-sensitivity` or `preprocessing-sensitivity` run.

## 6. Limitations

- **Descriptive correlations at n = 80 and 74.** No p-value is attached to any correlation; the flag is a reporting
  threshold. With 47 variables, several of them near-duplicates (`mitochondrial_mean` and `fraction_mitochondrial_umis_mean`
  are the same quantity; mean, SD and pseudobulk versions of each fraction move together), the tables describe one measurement
  space, not 47 independent tests. The two cohorts share 74 donors and are not independent replicates.
- **The supertype share confounds cell state with supertype-specific nucleus quality.** A small η² says the fraction does not
  differ much between astrocyte supertypes within a donor; it cannot separate a preparation effect from a state effect that is
  shared across supertypes, and a large η² could reflect quality differences between supertypes rather than biology. It is read
  beside the between-donor share and the nucleus QC columns, never alone.
- **Library structure.** 80 donors span several RNA libraries and 29 ATAC objects several ATAC libraries, so the
  library-level metrics are astrocyte-weighted means per donor; a library-level model would need the per-library values, which
  the step summarizes but does not export.
- **Missing values.** The hemoglobin supertype share is undefined for 6 donors (no variance), and 28 donors have nuclei
  without a `GEX_sequencing_saturation` value; correlations are pairwise-complete and the `n` column records the donors used. No
  donor lacks an ATAC object in this release.
- **Partial correlations.** They condition on the six analysis covariates by linear OLS, as the analysis does; a covariate itself
  has no partial correlation (n/a). The latent space is fitted without the covariates (MOTCO's covariate-blind basis), so a
  partial correlation that stays large says the adjustment did not touch the axis, not that it should have.
- **Geometry shares on unscaled scores.** MOTCO measures on the raw PLS scores, so the share of the geometry along `pls_0`
  reflects its larger score SD as well as its direction; the shares are a description of what the forced-rank statistics saw,
  not a decomposition of a test.
- **Uncommitted implementation at run time.** Both manifests record commit `3ae7553` with `dirty: true`; the diagnostics are
  deterministic (the cohort_50 rerun is byte-identical), so a rerun after the change is committed reproduces them.

## 7. Next steps

- `mtg-astrocyte-feature-selection-sensitivity`: predeclared export variants (a selection rule that is not raw log2-CPM variance,
  a ribosomal / mitochondrial family exclusion, residualization of the candidates named in section 5), each run as sensitivity
  beside the primary analysis with an adoption rule written in its proposal before any run.
- MOTCO side (its own repository): PLS takes the highest-variance axis first and CV may select rank 1 on it; unscaled scores
  weight components by variance in the geometry; an opt-in covariate residualization before the latent fit.
- Pathway/module enrichment remains exploratory without an orientation result.

## 8. Reproduction

```
# cluster login node (raw donor objects + results/pseudobulk/), after scripts/pseudobulk.sbatch
sbatch scripts/technical.sbatch            # or: nohup .venv/bin/sead-sexdiff technical-covariates --omic both &
# workstation
scripts/sync_pseudobulk.sh                 # pulls donor_technical.csv, nucleus_technical_summary.csv, technical_manifest.json
uv run sead-sexdiff axis-diagnostics       # both cohorts; reads the kept export, analysis and rank-sweep outputs
uv run sead-sexdiff axis-diagnostics --cohort cohort_50 --out-dir /tmp/rerun   # reproducibility check
```

The primary export, analysis and rank-sweep directories must be the ones whose hashes the analysis manifests record; the
diagnostic refuses to run otherwise and writes nothing under them.
