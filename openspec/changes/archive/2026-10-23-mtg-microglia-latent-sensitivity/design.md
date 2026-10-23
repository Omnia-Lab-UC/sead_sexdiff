# Design

## Context

See `proposal.md` — Why. The state that shapes the approach:

- **`run_latent_sensitivity`** has plain directory defaults (`results/motco_input`, `results/motco_analysis`, `results/pseudobulk`, `results/motco_latent_sensitivity`) and no `cell_type`; it calls `_check_reference`, `_check_technical_reference` and `_analyze_cohort` (which already takes `cell_type` and otherwise reads it from the export manifest), and writes the comparison manifest itself. The eight other MTG commands resolve their defaults with `_resolve_dir(value, cell_type, key)` over `cell_type_dirs` / `CELL_TYPE_SUBDIRS`, and check every manifest they read with `check_manifest_cell_type` (a mismatch is an error before anything is written; an absent key reads as `astrocyte` with a finding).
- **`run_positive_control`** resolves the control and reference roots from `POSITIVE_CONTROL_CELL_TYPES`, reads per root through `_positive_control_read(root, cohort, role, diagnostics=)` with the analysis and diagnostics subdirectories taken from `POSITIVE_CONTROL_EXPECTATIONS["E2"]["subdir"]` and `["E3"]["subdir"]` (`CELL_TYPE_SUBDIRS["motco_analysis"]`, `["motco_diagnostics"]`), checks the chain (stage-DE manifests and the diagnostic made from the analysis read for E2; the limma step from the DESeq2 summary), and writes `<out_dir>/<cohort>/` with `out_dir = <control_root>/positive_control`.
- **The astrocyte arms** exist at `results/motco_latent_sensitivity/` on `nixos-office` (report of 2026-10-21); their diagnostics at `results/motco_latent_sensitivity/diagnostics/<arm>/`. `nixos-home` holds the microglia chain and the astrocyte primary outputs, not the astrocyte arms. The arm manifests predate the registry and have no `cell_type` key.
- **MOTCO 0.9.0 contract**: with the covariate-blind basis the latent space, CV folds, selected rank and CV AUROC do not depend on the covariate set; under `stage` conditioning the folds equal the unadjusted analysis's folds. So `extended_covariates` repeats the primary's `pls_0` and `cv_mean_auroc`, and only `latent_adjusted` has a basis of its own.
- **`axis-diagnostics`** labels whatever analysis it reads as `primary` in its tables, so E3's row lookup (`analysis == DIAGNOSTIC_PRIMARY_LABEL`, `component == pls_0`) works on an arm's diagnostic unchanged.

## Goals / Non-Goals

**Goals:**
- `latent-sensitivity --cell-type microglia` is the astrocyte command with every default directory moved under `results/microglia/`, and nothing else different.
- The restated control is the primary-basis control with E2's and E3's source directories moved to the arm; the expectation constants are shared, not copied.
- Astrocyte outputs of `latent-sensitivity` and `positive-control` byte-identical under the defaults.

**Non-Goals:**
- A cell-type flag on `rank-sensitivity` or `preprocessing-sensitivity`.
- Any new expectation, threshold or technical-variable set.
- Rerunning the astrocyte arms or any primary analysis.

## Decisions

**D1. `run_latent_sensitivity(cohort, *, cell_type=DEFAULT_CELL_TYPE, input_dir=None, analysis_dir=None, pseudobulk_dir=None, out_dir=None, n_jobs=1, force=False)`.** Directory keywords default to `None` and resolve through `_resolve_dir` (as `run_analysis` does); `CELL_TYPE_SUBDIRS` gains `"motco_latent_sensitivity": "motco_latent_sensitivity"` if it is not there. The three manifests read (export, analysis, technical) go through `check_manifest_cell_type` before any write, their findings merged into the comparison manifest's `findings`; `_analyze_cohort(..., cell_type=cell_type)` so the arm manifests record it; the comparison manifest gains `"cell_type"`. The CLI mirrors it (`--cell-type` first, directory flags `default=None`, help text showing the astrocyte defaults), as for the eight commands.
- *Alternative: a separate microglia command.* Rejected: the arms, covariate set and checks are the same code, and the registry rule is one flag, one path.

**D2. Basis registry.**
```python
POSITIVE_CONTROL_BASES: dict[str, dict] = {
    "primary": {
        "description": "the primary analysis of each cell type and the diagnostic made from it",
        "analysis_subdir": CELL_TYPE_SUBDIRS["motco_analysis"],
        "diagnostics_subdir": CELL_TYPE_SUBDIRS["motco_diagnostics"],
        "arm": None,
    },
    "latent_adjusted": {
        "description": "the latent_adjusted arm of latent-sensitivity of each cell type and the diagnostic made from it",
        "analysis_subdir": f"{CELL_TYPE_SUBDIRS['motco_latent_sensitivity']}/latent_adjusted/{LATENT_SENSITIVITY_ANALYSIS_DIR}",
        "diagnostics_subdir": f"{CELL_TYPE_SUBDIRS['motco_latent_sensitivity']}/diagnostics/latent_adjusted",
        "arm": "latent_adjusted",
    },
}
POSITIVE_CONTROL_DEFAULT_BASIS = "primary"
```
An assertion at import ties `"latent_adjusted"` to a key of `LATENT_SENSITIVITY_ARMS` and requires `LATENT_SENSITIVITY_ARMS["latent_adjusted"]["latent_adjustment"]` to be set, so `extended_covariates` cannot be registered as a basis by mistake. `POSITIVE_CONTROL_EXPECTATIONS["E2"]["subdir"]` and `["E3"]["subdir"]` keep their values (they describe the primary basis); `_positive_control_read` takes a `basis` argument and uses the registry's subdirectories in their place. Nothing else in the expectations changes.
- *Alternative: a `--analysis-dir` / `--diagnostics-dir` pair on `positive-control`.* Rejected: a free directory could point the control at any analysis, and the reading rules are written for a named basis.

**D3. Output path.** `basis == "primary"` writes `<out_dir>/<cohort>/` as today; any other basis writes `<out_dir>/<basis>/<cohort>/`. The manifest records `basis` and the registry entry; the CSV gains no column (the basis is a directory level and a manifest key), so the primary-basis CSV stays byte-identical. `render_positive_control` prints the basis in its title line for a non-primary basis only.

**D4. Chain on the adjusted basis.** `_positive_control_read(..., basis="latent_adjusted")` reads, per root, the arm's `analysis_manifest.json` and the arm's comparison manifest (`<root>/motco_latent_sensitivity/comparison/<cohort>/latent_sensitivity_manifest.json`) and checks: cell type of both (an absent key reads as `astrocyte` with a finding, as everywhere; the astrocyte arm manifests predate the registry); the five export hashes recorded by the arm equal the files on disk under `<root>/motco_input/<cohort>/` (the same check the diagnostics run); the comparison manifest's reference analysis `result.json` sha256 equals the primary `result.json` on disk, which is the analysis E1's stage-DE manifests are checked against — so E1, E2 and E3 hang off one primary analysis per cell type; the diagnostic's manifest records the arm's `result.json` sha256 (the existing diagnostic chain check, pointed at the arm). E1 is read from the stage-DE summaries exactly as on the primary basis and its rows are identical by construction; the manifest records `e1_basis_independent: true`.

**D5. `extended_covariates` equality check in the report script, not in the step.** The report's check script reads, per cell type and cohort, `result.json` of the primary and of `extended_covariates` and compares `latent.n_components` (or the record's selected-rank key) and `latent.cv_mean_auroc`; equal → stated; different → a finding in the report (MOTCO's contract would be in question, which is a MOTCO issue, not resolved here). The step does not evaluate a restated control for that arm (D2's assertion).

**D6. Arm diagnostics layout.** `axis-diagnostics --cell-type microglia --input-dir results/microglia/motco_input --analysis-dir results/microglia/motco_latent_sensitivity/<arm>/motco_analysis --out-dir results/microglia/motco_latent_sensitivity/diagnostics/<arm>`, both arms, both cohorts; no rank directory (recorded as a finding, as for the astrocyte arms). The basis registry's `diagnostics_subdir` is this path.

**D7. Astrocyte arms on the evaluating machine.** Before the restated control runs, `results/motco_latent_sensitivity/` (arms, comparison and `diagnostics/`) is copied from `nixos-office` (or the cluster if it holds a verified copy) to the machine that runs it, and every file the astrocyte report cites is checked by sha256 against the report's provenance. If the arms cannot be retrieved, they are regenerated with `latent-sensitivity --force` at the recorded code and MOTCO versions, and the regenerated `result.json` records are compared with the astrocyte report's numbers (statistics, p-values, ranks, AUROC); the comparison is recorded in the microglia report's provenance, not a gate, as cross-machine floating-point differences are expected at rounding level. The primary analyses are never regenerated.

**D8. Byte-identity.** The synthetic dataset builder (scratchpad, rebuilt as in the previous changes) runs `latent-sensitivity` and `positive-control` before and after the refactor; CSVs and Markdown byte-identical, JSON equal apart from runtime, timestamps, git and the added `cell_type` / `basis` keys.

**D9. Report.** `reports/motco-microglia-latent-adjustment-<date>.md`, assembled by a scratch script from `latent_sensitivity.md`, `mediator_check.csv`, the `pls_0` rows of the diagnostics on the primary and both arms, `positive_control.md` of both bases, and the equality check; provenance and reading hand-written. A check script finds every statistic, p-value, rank and AUROC of the six microglia records and every value of both positive-control tables in the text.

## Risks / Trade-offs

- [The mitochondrial fraction is a mediator in microglia (ρ 0.50 with stage): the adjusted basis removes stage signal with the quality direction, E3 fails on the adjusted basis too, and the diagnosis is inconclusive] → predeclared in the proposal as a reading: the failure counts, the mediator caveat is stated, and the `extended_covariates` → `latent_adjusted` comparison shows how much of the trajectory statistics the basis alone changed.
- [Eight covariates in cells of 9–24 microglia donors near-singular] → MOTCO errors on a rank-deficient design; an arm failure fails the cohort without partial output and is reported, as in the astrocyte change.
- [The astrocyte arms cannot be retrieved byte-identically] → D7: regenerate at the recorded versions, record the comparison as provenance, never regenerate a primary.
- [The restated control holds and is read as adopting the arm] → the spec binds the wording: MOTCO's `candidate_fails` is quoted, the primary-basis outcome and the astrocyte reading are unchanged.
- [`CELL_TYPE_SUBDIRS` or the comparison manifest gain keys that change astrocyte bytes] → D8's synthetic check; the astrocyte arm outputs on disk are not rewritten.

## Migration Plan

1. Code (D1–D4), synthetic byte-identity check (D8), on any machine.
2. On the workstation holding the microglia chain: `latent-sensitivity --cell-type microglia --cohort all`, `axis-diagnostics` on both arms (D6), copy and check the astrocyte arms (D7), `positive-control --basis latent_adjusted --cohort all`.
3. Report (D9), README and AGENTS.md, check, commit; copy `results/microglia/motco_latent_sensitivity/` and `results/microglia/positive_control/latent_adjusted/` to the cluster and verify.

Rollback: revert the commit; no primary output is touched, and the new directories are removed.
