## MODIFIED Requirements

### Requirement: Predeclared arms
The step SHALL run exactly two arms per cohort, each once:
- `extended_covariates`, with the covariate-blind basis (no latent adjustment);
- `latent_adjusted`, with MOTCO's covariate latent adjustment conditioned on `stage`.

The arms SHALL be a fixed registry, not a runtime option. `group_stage` conditioning SHALL NOT be run. Apart from the covariate set and the latent-adjustment setting, both arms SHALL use the primary analysis's settings: the same export, blocks and units, group and stage columns, stage order, orientation estimator, CV rank selection, residual correction, seed and permutation count.

#### Scenario: Arms complete
- **WHEN** the step finishes for `cohort_30`
- **THEN** two arm records exist: `extended_covariates` without a `latent_adjustment` block, and `latent_adjusted` with one whose `conditioning` is `stage`; both have `component_selection = cv`, `residual_correction = df` with a reduced rank of 12 (eight covariates), the primary's seed and 999 permutations

#### Scenario: Latent adjustment changes only the basis
- **WHEN** both arms of a cohort have finished
- **THEN** the two records carry the same covariate block, and any difference between their statistics arises through the latent basis and the CV-selected rank only, as MOTCO's latent-adjustment contract guarantees

### Requirement: Report keeps the primary analysis primary
The report SHALL state that:
- the primary analysis stays primary, and that MOTCO's power comparison (`candidate_fails`) excludes latent adjustment as its replacement;
- a rejection status that differs between an arm and the primary is covariate-set dependence (primary vs `extended_covariates`) or basis dependence (`extended_covariates` vs `latent_adjusted`), not a result of the arm;
- a CV-selected rank that differs is reported beside `angle`, not as better or worse;
- `angle` and `shape` non-rejections stay uninformative and a `delta` rejection is a difference, not specifically of size, in every arm;
- the `latent_adjusted` arm's `angle` and `shape` p-values are reported but not read, because MOTCO's case-study diagnosis (`residualized_basis`) found the `stage` latent adjustment liberal on both at this configuration, an excess of the basis that the residual correction does not address: neither a rejection nor a non-rejection of either is interpretable.

A technical correlation of an arm's first axis SHALL be described as a measurement-space finding only, never as a reinterpretation of any arm's statistics, and never as grounds to adopt an arm.

#### Scenario: Arm rejects where the primary does not
- **WHEN** `extended_covariates` rejects `delta` in `cohort_30` and the primary does not
- **THEN** the report calls this covariate-set dependence, reports the arm's selected rank beside `angle`, and keeps the primary's non-rejection as the result

#### Scenario: Latent-adjusted angle or shape
- **WHEN** `latent_adjusted` rejects `angle` or `shape` in a cohort, or does not
- **THEN** the report gives the p-value, states that it is not interpretable under MOTCO's `residualized_basis` diagnosis, and draws no reading from it, while the arm's `delta` is read under the rules above

### Requirement: Report carries the limitations
The report SHALL carry:
- the mediator limitation, with the mediator check's flags: a flagged fraction may carry stage signal that the `latent_adjusted` basis removed, and no simulation can tell a confounder from a mediator;
- MOTCO's latent-adjustment operating characteristics, measured at the Phase 6 design only: Type I `adopt_stage`; power `candidate_fails`, with `shape` power lost under a nuisance axis and small-magnitude power lost without one;
- MOTCO's case-study diagnosis of the `stage` latent adjustment at this configuration (`residualized_basis`: anchor `angle` 0.084–0.090 and `shape` 0.078–0.092 against 0.0695, `angle` still 0.084 / 0.082 at a forced rank of 2), which makes that arm inadmissible for `angle` and `shape` here;
- the facts that no operating characteristic was measured on the SEA-AD blocks or at this eight-covariate set, that the residual-corrected null's size was measured at six covariates only and is inferred at eight (reduced rank 12), that the two arms carry no Type I or power claim of their own, and that the nuclei counts and the fractions both proxy library quality, which the report states.

#### Scenario: Mediator flag reported
- **WHEN** a fraction is flagged in `mediator_check.csv`
- **THEN** the report names the fraction, the cohort and ρ, and states that the `latent_adjusted` arm may have removed stage signal carried by it

#### Scenario: Diagnosis carried
- **WHEN** the re-issue report presents the latent-adjustment sensitivity
- **THEN** it cites MOTCO's diagnosis report with the anchor rates above and states that the `latent_adjusted` arm's `angle` and `shape` are not read
