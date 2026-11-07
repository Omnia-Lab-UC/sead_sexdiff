# Spec Delta

## ADDED Requirements

### Requirement: Split-test figure
Figure 11 (`11_span_split`) SHALL draw, from the analysis records only:
- per primary analysis (both cell types and cohorts), the pseudo-F p-value of each part, marked by its role (`complement` confirmatory, `span` the concordance cell with its status beside the production statistics', `total` descriptive), with α;
- the complement p-value of every registered analysis that carries the test, enumerated from the registries, with each `latent_adjusted` arm listed as not available;
- the in-span share of every such analysis against its fitted rank, with no reference line.

#### Scenario: Analyses drawn
- **WHEN** figure 11 is drawn from the re-issued outputs
- **THEN** it holds 40 complement p-values (44 analyses minus the four `latent_adjusted` arms), and its caption gives the four primaries' three p-values, ranks and in-span shares read from `result.json`

#### Scenario: Share not read
- **WHEN** figure 11 is drawn
- **THEN** its caption states that the in-span share is descriptive, grows with the fitted rank, and is read against no null here

## MODIFIED Requirements

### Requirement: Registered figure set
The `figures` command SHALL draw exactly the figures of the registry `FIGURES`, in registry order, each written as `<nn>_<name>` (`01_cohorts` … `11_span_split`). No flag SHALL select, skip or add a figure. Entries SHALL be added or changed by a change proposal only.

#### Scenario: Default run
- **WHEN** `sead-sexdiff figures` runs with every input present
- **THEN** eleven figures `01_cohorts` to `11_span_split` are written, each as PNG, PDF and SVG, with one caption each in `figures.md`

### Requirement: One measurement across the analyses drawn
Every MOTCO analysis record drawn SHALL carry the same MOTCO version and the registered measurement. That is:
- the residual correction `ANALYSIS_RESIDUAL_CORRECTION`;
- the rank floor `ANALYSIS_RANK_FLOOR`, or a forced rank of `ANALYSIS_FORCED_RANKS` in `rank-sensitivity`;
- `ANALYSIS_PERMUTATIONS` and `ANALYSIS_SEED`;
- the split test `ANALYSIS_SPAN_SPLIT_TEST`, except in the `latent_adjusted` arms, whose records SHALL NOT carry it.

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

### Requirement: Binding reporting context drawn with MOTCO results
Every figure that draws a MOTCO statistic SHALL show, in the figure itself, the covariate block, the feature filtering, the RRPP permutations and residual correction, and the non-rejection rules. Every `angle` SHALL be shown with its CV-selected rank, and in figure 04 with its eigengap and null summary. The `latent_adjusted` arm's `angle` and `shape` p-values SHALL be marked as reported, not read. Figure 11 SHALL show the split test's multiplicity policy and state that a complement non-rejection is uninformative.

#### Scenario: Overview figure
- **WHEN** figure 05 is drawn
- **THEN** each row shows the selected or forced rank per cohort, the `latent_adjusted` `angle` and `shape` markers are drawn as not read, and the footnote states that sensitivity analyses are read as dependence and that `angle` / `shape` non-rejections are uninformative

#### Scenario: Split figure
- **WHEN** figure 11 is drawn
- **THEN** its footnote states that only the complement is confirmatory, that the span part is a concordance cell, that `total` and the in-span share are not read, that sensitivity analyses are read as dependence, and that a complement non-rejection is uninformative
