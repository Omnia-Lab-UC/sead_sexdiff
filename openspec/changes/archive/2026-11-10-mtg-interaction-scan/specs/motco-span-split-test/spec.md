## MODIFIED Requirements

### Requirement: Complement is the one confirmatory part
The `complement` part's pseudo-F p-value SHALL be the split test's only confirmatory reading. It is read at `ANALYSIS_ALPHA` under `ANALYSIS_REJECTION_RULE` (p ≤ 0.05), unadjusted. No adjustment SHALL be made across `delta`, `angle`, `shape`, the complement and the interaction scan's global p-value (spec `motco-interaction-scan`), or across cohorts and cell types. The report SHALL state that each primary carries five tests read at α. This policy is predeclared and SHALL NOT be changed from results.

#### Scenario: Complement rejected in a primary
- **WHEN** the `cohort_30` primary of a cell type has a complement pseudo-F p-value ≤ 0.05
- **THEN** the report gives it as that cell type's split-test result under the reading requirement, whatever `total` and `span` read

#### Scenario: Five tests stated
- **WHEN** the report gives a primary's results
- **THEN** it states that the primary carries five tests read at α with no adjustment across them: `delta`, `angle`, `shape`, the complement and the scan
