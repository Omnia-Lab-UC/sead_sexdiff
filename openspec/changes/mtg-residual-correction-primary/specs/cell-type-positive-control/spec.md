## ADDED Requirements

### Requirement: Control outcome re-issued on the corrected primary
When the primary analyses are re-issued under the residual correction, the control SHALL be evaluated again on both bases with E1–E3, their comparisons, technical variables and the all-three rule unchanged. The corrected outcome SHALL be the primary-basis outcome from then on, with the superseded outcome beside it per expectation and cohort. Because no expectation reads a p-value, E1–E3's values SHALL be compared with the superseded ones.

#### Scenario: Values reproduce
- **WHEN** E1–E3 are evaluated on the corrected primary and on the `latent_adjusted` arms re-issued with the correction
- **THEN** the report states, per basis, expectation and cohort, that each value equals the superseded one, and that the control's outcome is unchanged on both bases

#### Scenario: Value differs
- **WHEN** an E2 CV mean AUROC or an E3 correlation differs from the superseded value beyond rounding
- **THEN** the report records it as a finding with both values, the outcome follows the all-three rule on the new values, and no expectation is changed

#### Scenario: Restated control and the basis diagnosis
- **WHEN** the report presents the restated control on `latent_adjusted`
- **THEN** it states that E2 and E3 read the arm's CV AUROC and first-component correlations only, which MOTCO's `residualized_basis` diagnosis (a size excess of the arm's `angle` and `shape` tests) does not bear on, so the restated control keeps its predeclared reading
