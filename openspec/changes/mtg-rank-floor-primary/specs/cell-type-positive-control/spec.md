# Spec Delta

## ADDED Requirements

### Requirement: Control outcome re-issued on the floored primary
When the primary analyses are re-issued under the rank floor, the control SHALL be evaluated again on both bases with E1–E3, their comparisons, technical variables and the all-three rule unchanged. The floored outcome SHALL be the primary-basis outcome from then on; the report SHALL give the superseded outcome beside it, per expectation and cohort, with the values compared.

#### Scenario: Outcome unchanged
- **WHEN** E1–E3 give the same outcomes on the floored primary as on the superseded one in `cohort_30`
- **THEN** the report states that the floor leaves the control's outcome as it was, with both sets of values

#### Scenario: Outcome changes through E2
- **WHEN** the floored microglia CV mean AUROC no longer exceeds the astrocytes' and the superseded one did
- **THEN** E2 is recorded as not holding on the floored primary, the overall outcome follows the all-three rule, the superseded outcome is reported beside it as dependence on the rank rule, and no expectation is changed

#### Scenario: First component reproduces
- **WHEN** E3 is evaluated on the floored primary
- **THEN** the report states whether the first-component correlations equal the superseded ones (PLS components nest across ranks), and a difference beyond rounding is recorded as a finding
