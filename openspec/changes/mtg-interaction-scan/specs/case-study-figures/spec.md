## MODIFIED Requirements

### Requirement: Registered figure set
The `figures` command SHALL draw exactly the figures of the registry `FIGURES`, in registry order, each written as `<nn>_<name>` (`01_cohorts` … `12_interaction_scan`). No flag SHALL select, skip or add a figure. Entries SHALL be added or changed by a change proposal only.

#### Scenario: Default run
- **WHEN** `sead-sexdiff figures` runs with every input present
- **THEN** twelve figures `01_cohorts` to `12_interaction_scan` are written, each as PNG, PDF and SVG, with one caption each in `figures.md`

### Requirement: One measurement across the analyses drawn
Every MOTCO analysis record drawn SHALL carry the same MOTCO version and the registered measurement. That is:
- the residual correction `ANALYSIS_RESIDUAL_CORRECTION`;
- the rank floor `ANALYSIS_RANK_FLOOR`, or a forced rank of `ANALYSIS_FORCED_RANKS` in `rank-sensitivity`;
- `ANALYSIS_PERMUTATIONS` and `ANALYSIS_SEED`;
- the split test `ANALYSIS_SPAN_SPLIT_TEST`, except in the `latent_adjusted` arms, whose records SHALL NOT carry it;
- the interaction scan `ANALYSIS_INTERACTION_SCAN`, in every record.

A record that differs SHALL stop the run naming it and the differing setting.

#### Scenario: Re-issue half finished
- **WHEN** the microglia primary records MOTCO 0.18.0 and the astrocyte primary MOTCO 0.13.0
- **THEN** the run stops naming both records and versions, and no figure mixes the two

#### Scenario: Uncorrected analysis
- **WHEN** an arm's `result.json` has no `params.residual_correction`
- **THEN** the run stops naming the file and the correction expected (`df`)

#### Scenario: Analysis without the split test
- **WHEN** a preprocessing variant's `result.json` has no `span_split` block
- **THEN** the run stops naming the file and the split test expected

#### Scenario: Analysis without the scan
- **WHEN** a `latent_adjusted` arm's `result.json` has no `interaction_scan` block
- **THEN** the run stops naming the file and the scan expected

### Requirement: Binding reporting context drawn with MOTCO results
Every figure that draws a MOTCO statistic SHALL show, in the figure itself, the covariate block, the feature filtering, the RRPP permutations and residual correction, and the non-rejection rules. Every `angle` SHALL be shown with its CV-selected rank, and in figure 04 with its eigengap and null summary. The `latent_adjusted` arm's `angle` and `shape` p-values SHALL be marked as reported, not read. Figure 11 SHALL show the split test's multiplicity policy and state that a complement non-rejection is uninformative. Figure 12 SHALL show the scan's multiplicity policy and the reading of a non-rejection.

#### Scenario: Overview figure
- **WHEN** figure 05 is drawn
- **THEN** each row shows the selected or forced rank per cohort, the `latent_adjusted` `angle` and `shape` markers are drawn as not read, and the footnote states that sensitivity analyses are read as dependence and that `angle` / `shape` non-rejections are uninformative

#### Scenario: Split figure
- **WHEN** figure 11 is drawn
- **THEN** its footnote states that only the complement is confirmatory, that the span part is a concordance cell, that `total` and the in-span share are not read, that sensitivity analyses are read as dependence, and that a complement non-rejection is uninformative

#### Scenario: Scan figure
- **WHEN** figure 12 is drawn
- **THEN** its footnote states that the global p-value is confirmatory with five tests per primary and no adjustment across them, that hits are familywise at 0.05, and that the per-stage differences and null-maximum shares are not read. It also states that a non-rejection carries the power statement for the astrocyte `cohort_30` primary only and is uninformative elsewhere

## ADDED Requirements

### Requirement: Interaction-scan figure
Figure 12 (`12_interaction_scan`) SHALL draw, from the analysis records and scan tables only:
- per primary analysis, the smallest adjusted p-values per block on a log axis with α, hits marked;
- the global scan p-value of every registered analysis, enumerated from the registries, with the forced ranks and the `latent_adjusted` arms marked as identical to their reference by construction;
- per primary, the share of null maxima per block beside its share of features, with no threshold.

#### Scenario: Analyses drawn
- **WHEN** figure 12 is drawn from the re-issued outputs
- **THEN** it holds 44 global p-values, 28 of them distinct analyses. Its caption gives, per primary, the global p-value, the maximum F with its feature and block, the hits per block and the null-maximum shares, read from `result.json`, and each cohort's scan concordance cells read from the stage-DE summaries

#### Scenario: Table bound to its record
- **WHEN** a scan table read for figure 12 has a different number of rows or hits than its record states
- **THEN** the run stops naming the table and the record, and nothing is written
