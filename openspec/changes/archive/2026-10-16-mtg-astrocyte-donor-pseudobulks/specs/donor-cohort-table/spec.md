# Spec Delta

## Purpose

Defines the donor-level metadata table of the case study: sex as the group, ADNC merged into three ordered stages, astrocyte nucleus counts per omic, and the ≥30 and ≥50 nuclei cohort flags that select the primary and sensitivity cohorts.

## ADDED Requirements

### Requirement: One row per donor across both omics
The system SHALL write `results/pseudobulk/donor_table.csv` with one row per donor present in either omic, with columns `donor_id`, `sex`, `adnc_raw`, `stage`, `n_rna`, `n_atac`, `cohort_30`, `cohort_50`. A donor absent from one omic SHALL have a nucleus count of 0 for that omic.

#### Scenario: Donor only in ATAC
- **WHEN** a donor exists among the 86 ATAC objects but not among the 84 RNA objects
- **THEN** its row has `n_rna = 0` and both cohort flags false

### Requirement: ADNC merged into three ordered stages
The system SHALL recode `Overall AD neuropathological Change` into the string stage labels `0` (`Not AD` and `Low`), `1` (`Intermediate`) and `2` (`High`). Any other value SHALL be an error naming the donor and the value. The recoding SHALL be recorded in the manifest.

#### Scenario: Low merges with Not AD
- **WHEN** a donor's ADNC is `Low`
- **THEN** its `stage` is `0`

#### Scenario: Unexpected ADNC value
- **WHEN** a donor's ADNC is a value outside the four expected labels
- **THEN** the run stops with an error naming the donor and the value

### Requirement: Sex has exactly two labels and agrees across omics
The system SHALL require that `sex` takes exactly two distinct values across the table and that, for donors present in both omics, sex and ADNC agree between the RNA and ATAC objects. Disagreement SHALL be an error naming the donor.

#### Scenario: Metadata disagreement
- **WHEN** a donor is labelled `Female` in the RNA object and `Male` in the ATAC object
- **THEN** the run stops with an error naming the donor and both values

#### Scenario: Single sex label in a partial run
- **WHEN** a partial run (donor limit, donor list, a single omic, or fewer files than the full release) yields only one sex label
- **THEN** the run records a warning in the manifest instead of failing; more than two labels, or fewer than two in a full run, still stop the run

### Requirement: Cohort flags require both omics
The system SHALL set `cohort_30` true when `n_rna >= 30` and `n_atac >= 30`, and `cohort_50` true when `n_rna >= 50` and `n_atac >= 50`. Thresholds SHALL be recorded in the manifest.

#### Scenario: Borderline donor
- **WHEN** a donor has `n_rna = 45` and `n_atac = 60`
- **THEN** `cohort_30` is true and `cohort_50` is false

### Requirement: Cohort summary
The system SHALL write `results/pseudobulk/cohort_summary.csv` with, for each cohort, the count of donors per sex × stage cell and the cohort total, and SHALL print the same table at the end of the run. The manifest SHALL record whether every sex × stage cell is non-empty in each cohort.

#### Scenario: Summary matches the preregistered cohort
- **WHEN** the table is built from the complete MTG releases
- **THEN** the ≥30 cohort summary reports n = 80 with F 11/10/28 and M 9/10/12 across stages 0/1/2, and the ≥50 cohort reports n = 74; any deviation is written to the manifest as a finding rather than silently accepted

#### Scenario: Empty cell is flagged
- **WHEN** a sex × stage cell has zero donors in a cohort
- **THEN** the manifest marks that cohort as not analysable and the summary prints a warning
