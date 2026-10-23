# Spike: microglia labels and cohort cells (2026-10-22)

Metadata only (task 1.2, design D6): the `obs` columns `Donor ID`, `Sex`, `Overall AD neuropathological Change`,
`Subclass`, `Supertype` (RNA also `Genes detected`) of all 84 RNA (release 2026-06-22) and 86 ATAC (release
2024-12-06) donor objects. No count matrix, no `var`, no pseudobulk.

## Where it ran

The cluster was unreachable on 2026-10-22 (SSH banner-exchange timeouts, then TCP timeouts), so the spike read the
same release objects from the public bucket instead of the cluster's copies: each object was opened with h5py over
HTTPS range requests (1 MiB blocks, only the blocks the `obs` columns occupy), and its size and ETag had to equal
the ones the astrocyte pseudobulk manifest records for it (`results/pseudobulk/manifest.json`
`releases.<omic>.files`, the cluster's download manifest), so these are the objects the cluster's pseudobulk
read. Every one of the 170 matched.

## Objects read

| omic | objects | bytes_object | requests | bytes_read |
|---|---|---|---|---|
| atac | 86 | 17933512911 | 396 | 370765381 |
| rna | 84 | 32788099111 | 459 | 436375234 |

(bytes_object: the objects' total size; bytes_read: what the range requests transferred.)

Check of the obs-only nucleus count: the same script tallied the astrocyte nuclei of every donor (RNA:
`Genes detected` >= `RNA_MIN_GENES` = 200, which the astrocyte spike found equal to the per-nucleus nonzero count
of `UMIs`; ATAC: no nucleus QC) and passed them through the package's `build_donor_table`. All 84 non-Reference
donors' `n_rna` and `n_atac` equal the astrocyte `results/pseudobulk/donor_table.csv` (0 mismatches), so the
microglia counts below are the counts the pseudobulk step will sum, up to the label decision.

## Subclass inventory against TAXONOMY

- atac 2024-12-06: 24 Subclass labels; not in TAXONOMY: none; in TAXONOMY but absent: none
- rna 2026-06-22: 24 Subclass labels; not in TAXONOMY: none; in TAXONOMY but absent: none

## Supertype labels of the immune / microglia subclass (all donors, Reference included)

| omic | release | subclass | supertype | nuclei | nuclei_qc | donors |
|---|---|---|---|---|---|---|
| atac | 2024-12-06 | Microglia-PVM | Lymphocyte | 219 |  | 54 |
| atac | 2024-12-06 | Microglia-PVM | Micro-PVM_1 | 447 |  | 74 |
| atac | 2024-12-06 | Microglia-PVM | Micro-PVM_2 | 6130 |  | 86 |
| atac | 2024-12-06 | Microglia-PVM | Micro-PVM_2_1-SEAAD | 35 |  | 13 |
| atac | 2024-12-06 | Microglia-PVM | Micro-PVM_2_3-SEAAD | 1796 |  | 82 |
| atac | 2024-12-06 | Microglia-PVM | Micro-PVM_3-SEAAD | 910 |  | 81 |
| atac | 2024-12-06 | Microglia-PVM | Micro-PVM_4-SEAAD | 69 |  | 19 |
| atac | 2024-12-06 | Microglia-PVM | Monocyte | 57 |  | 27 |
| rna | 2026-06-22 | Immune | Lymphocyte | 767 | 767 | 79 |
| rna | 2026-06-22 | Immune | Micro-PVM_1 | 1383 | 1383 | 83 |
| rna | 2026-06-22 | Immune | Micro-PVM_2 | 21090 | 21090 | 84 |
| rna | 2026-06-22 | Immune | Micro-PVM_2_1-SEAAD | 282 | 282 | 44 |
| rna | 2026-06-22 | Immune | Micro-PVM_2_3-SEAAD | 8278 | 8278 | 84 |
| rna | 2026-06-22 | Immune | Micro-PVM_3-SEAAD | 4468 | 4468 | 84 |
| rna | 2026-06-22 | Immune | Micro-PVM_4-SEAAD | 406 | 406 | 48 |
| rna | 2026-06-22 | Immune | Monocyte | 12 | 12 | 9 |

Microglia / PVM-named supertypes under any other Subclass: none

The RNA `Immune` and the ATAC `Microglia-PVM` subclasses carry the same eight `Supertype` labels. Every RNA
immune nucleus passes the 200-gene filter.

## Label decision (design Open Question)

- Neither release has a PVM-only supertype: the six `Micro-PVM_*` labels bundle microglia and perivascular
  macrophages in both releases. Following the design's default, every one of them is kept, so the population is
  the same set of labels in both omics.
- Population (`CELL_TYPES["microglia"]`, identity raw → canonical mapping in both releases): `Micro-PVM_1`,
  `Micro-PVM_2`, `Micro-PVM_2_1-SEAAD`, `Micro-PVM_2_3-SEAAD`, `Micro-PVM_3-SEAAD`, `Micro-PVM_4-SEAAD`.
- Recorded non-population labels of the subclass (`other_supertypes`, never selected): `Lymphocyte`,
  `Monocyte`. Any other label of the subclass stops the run.
- Subclass: `Immune` (RNA 2026-06-22), `Microglia-PVM` (ATAC 2024-12-06).

## Microglia / PVM nuclei per donor (RNA after Genes detected >= 200; Reference donors excluded)

| donor_id | sex | adnc_raw | stage | n_rna | n_atac | cohort_30 | cohort_50 |
|---|---|---|---|---|---|---|---|
| H19.33.004 | Female | Not AD | 0 | 711 | 19 | False | False |
| H20.33.001 | Male | Low | 0 | 782 | 216 | True | True |
| H20.33.002 | Female | Not AD | 0 | 526 | 176 | True | True |
| H20.33.004 | Male | High | 2 | 462 | 249 | True | True |
| H20.33.005 | Female | Intermediate | 1 | 777 | 192 | True | True |
| H20.33.008 | Female | High | 2 | 184 | 226 | True | True |
| H20.33.011 | Female | High | 2 | 386 | 49 | True | False |
| H20.33.012 | Female | Low | 0 | 472 | 162 | True | True |
| H20.33.013 | Male | Intermediate | 1 | 345 | 139 | True | True |
| H20.33.014 | Female | Intermediate | 1 | 262 | 156 | True | True |
| H20.33.015 | Male | Intermediate | 1 | 348 | 167 | True | True |
| H20.33.016 | Female | Intermediate | 1 | 260 | 118 | True | True |
| H20.33.017 | Male | High | 2 | 504 | 274 | True | True |
| H20.33.018 | Female | High | 2 | 519 | 89 | True | True |
| H20.33.019 | Female | Intermediate | 1 | 331 | 116 | True | True |
| H20.33.020 | Male | High | 2 | 569 | 175 | True | True |
| H20.33.024 | Male | Intermediate | 1 | 117 | 59 | True | True |
| H20.33.025 | Male | High | 2 | 727 | 222 | True | True |
| H20.33.026 | Female | High | 2 | 149 | 164 | True | True |
| H20.33.027 | Female | Intermediate | 1 | 214 | 121 | True | True |
| H20.33.028 | Female | High | 2 | 984 | 337 | True | True |
| H20.33.029 | Female | High | 2 | 347 | 157 | True | True |
| H20.33.030 | Female | High | 2 | 595 | 73 | True | True |
| H20.33.031 | Female | High | 2 | 343 | 101 | True | True |
| H20.33.032 | Male | High | 2 | 208 | 148 | True | True |
| H20.33.033 | Male | High | 2 | 158 | 33 | True | False |
| H20.33.034 | Female | Intermediate | 1 | 340 | 59 | True | True |
| H20.33.035 | Female | Not AD | 0 | 191 | 30 | True | False |
| H20.33.036 | Female | High | 2 | 138 | 51 | True | True |
| H20.33.037 | Female | High | 2 | 63 | 62 | True | True |
| H20.33.038 | Female | High | 2 | 288 | 149 | True | True |
| H20.33.039 | Female | High | 2 | 384 | 235 | True | True |
| H20.33.040 | Male | Intermediate | 1 | 509 | 133 | True | True |
| H20.33.041 | Female | High | 2 | 259 | 114 | True | True |
| H20.33.043 | Male | Intermediate | 1 | 731 | 478 | True | True |
| H20.33.044 | Male | Not AD | 0 | 315 | 64 | True | True |
| H20.33.045 | Female | High | 2 | 238 | 77 | True | True |
| H20.33.046 | Male | High | 2 | 483 | 51 | True | True |
| H21.33.001 | Male | Low | 0 | 814 | 116 | True | True |
| H21.33.002 | Female | High | 2 | 419 | 38 | True | False |
| H21.33.003 | Male | Not AD | 0 | 442 | 192 | True | True |
| H21.33.004 | Male | Not AD | 0 | 317 | 191 | True | True |
| H21.33.005 | Male | Intermediate | 1 | 536 | 113 | True | True |
| H21.33.006 | Male | Intermediate | 1 | 782 | 257 | True | True |
| H21.33.007 | Female | High | 2 | 344 | 32 | True | False |
| H21.33.008 | Female | High | 2 | 694 | 425 | True | True |
| H21.33.009 | Female | High | 2 | 245 | 61 | True | True |
| H21.33.010 | Female | High | 2 | 895 | 43 | True | False |
| H21.33.011 | Female | Not AD | 0 | 460 | 43 | True | False |
| H21.33.012 | Female | Intermediate | 1 | 207 | 62 | True | True |
| H21.33.013 | Female | High | 2 | 471 | 115 | True | True |
| H21.33.014 | Male | Intermediate | 1 | 232 | 112 | True | True |
| H21.33.015 | Male | Low | 0 | 425 | 103 | True | True |
| H21.33.016 | Female | Low | 0 | 245 | 105 | True | True |
| H21.33.017 | Female | High | 2 | 863 | 53 | True | True |
| H21.33.018 | Female | Low | 0 | 746 | 74 | True | True |
| H21.33.019 | Male | Low | 0 | 263 | 269 | True | True |
| H21.33.020 | Male | High | 2 | 79 | 24 | False | False |
| H21.33.021 | Male | Intermediate | 1 | 241 | 115 | True | True |
| H21.33.022 | Female | Intermediate | 1 | 883 | 79 | True | True |
| H21.33.023 | Male | Not AD | 0 | 382 | 47 | True | False |
| H21.33.025 | Female | Intermediate | 1 | 866 | 104 | True | True |
| H21.33.026 | Female | High | 2 | 285 | 13 | False | False |
| H21.33.027 | Male | High | 2 | 1034 | 41 | True | False |
| H21.33.028 | Male | Low | 0 | 748 | 54 | True | True |
| H21.33.029 | Male | High | 2 | 183 | 25 | False | False |
| H21.33.030 | Male | Intermediate | 1 | 240 | 101 | True | True |
| H21.33.031 | Male | High | 2 | 401 | 87 | True | True |
| H21.33.032 | Female | Low | 0 | 320 | 89 | True | True |
| H21.33.033 | Female | High | 2 | 584 | 15 | False | False |
| H21.33.034 | Female | High | 2 | 209 | 32 | True | False |
| H21.33.035 | Female | High | 2 | 404 | 43 | True | False |
| H21.33.036 | Female | High | 2 | 307 | 29 | False | False |
| H21.33.037 | Female | Low | 0 | 263 | 55 | True | True |
| H21.33.038 | Female | Low | 0 | 156 | 34 | True | False |
| H21.33.039 | Female | High | 2 | 197 | 11 | False | False |
| H21.33.040 | Male | High | 2 | 372 | 68 | True | True |
| H21.33.041 | Female | Not AD | 0 | 444 | 69 | True | True |
| H21.33.042 | Female | High | 2 | 433 | 61 | True | True |
| H21.33.043 | Female | Low | 0 | 962 | 96 | True | True |
| H21.33.044 | Female | Intermediate | 1 | 386 | 81 | True | True |
| H21.33.045 | Female | High | 2 | 287 | 18 | False | False |
| H21.33.046 | Male | High | 2 | 242 | 61 | True | True |
| H21.33.047 | Male | Intermediate | 1 | 380 | 33 | True | False |

RNA microglia / PVM nuclei per donor range from 63 to 1,034 (median 376); ATAC from 11 to 478 (median 88), with
8 donors below 30 and 21 below 50. The ATAC side sets the cohorts.

## Cohort cells (the preregistered microglia cells)

| cohort | sex | stage 0 | stage 1 | stage 2 | total |
|---|---|---|---|---|---|
| cohort_30 | Female | 11 | 10 | 24 | 45 |
| cohort_30 | Male | 9 | 11 | 11 | 31 |
| cohort_30 | all |  |  |  | 76 |
| cohort_50 | Female | 8 | 10 | 18 | 36 |
| cohort_50 | Male | 8 | 10 | 9 | 27 |
| cohort_50 | all |  |  |  | 63 |

empty cells: {'cohort_30': [], 'cohort_50': []}

`k` (smallest sex × stage cell): 9 in `cohort_30`, 8 in `cohort_50`.

## Thresholds (task 1.3)

No sex × stage cell is empty at >= 30 (nor at >= 50), so the astrocyte thresholds are kept for microglia:
`COHORT_THRESHOLDS` (>= 30 / >= 50 nuclei in both omics), `ADNC_STAGES`, `k` per cohort as for the astrocytes.
The proposal is unchanged. The registry's preregistered microglia cells are the cells above:
`cohort_30` n = 76, Female 11/10/24, Male 9/11/11; `cohort_50` n = 63, Female 8/10/18, Male 8/10/9.
