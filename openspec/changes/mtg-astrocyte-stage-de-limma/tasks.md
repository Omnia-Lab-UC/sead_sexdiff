# Tasks

## 1. Environment, constants and synthetic baseline

- [ ] 1.1 Add `rPackages.limma` to `flake.nix` beside `rPackages.edgeR`. Verify: `nix develop -c Rscript -e 'packageVersion("limma"); packageVersion("edgeR")'` prints 3.68.4 and 4.10.1 (R 4.6.1), and the versions are recorded here. If they differ, record them and re-read the defaults that design Context lists.
- [ ] 1.2 Add the constants to `data.py`: `STAGE_DE_LIMMA_METHOD` (design D3–D5), `STAGE_DE_LIMMA_JOINT_TESTS` (D6), `STAGE_DE_LIMMA_CONCORDANCE_RULE` (D8), `STAGE_DE_LIMMA_TMM_ATOL = 1e-8` and the `STAGE_DE_LIMMA_*` output names (D12). Add import-time assertions:
  - the joint tests name exactly the two tested labels of each factor model of `STAGE_DE_MODELS`;
  - the concordance rule names `factor_interaction` / `interaction_joint`;
  - the method's TMM trims are `TMM_LOGRATIO_TRIM` / `TMM_SUM_TRIM`.

  Verify by importing `sead_sexdiff.data` and printing the registries, and check that a patched joint-test label fails the import.
- [ ] 1.3 Synthetic baseline. Copy `make_synthetic.py`, `run_export_analyze.py`, `cal_sim.py` and `cal_sim_syn.py` from the `mtg-astrocyte-stage-de` session scratchpad (`/tmp/claude-1000/-home-tomas-Documents-Repos-sead-sexdiff/a4d4b57f-2de0-45ae-99ba-1e3c848dcc49/scratchpad/scripts/`) if it still exists. Otherwise rebuild them to the spec in that change's task 1.2 (seed 20261021, 60 donors, 10 per sex × stage cell, the planted stage and female-only interaction features). Then run `export`, `analyze` (19 permutations) and `stage-de --cohort all --omic both` with the unmodified package, and record the sha256 of every `stage-de` output as the byte-identity baseline. Verify that the statuses and planted features match that change's task 1.2 record (RNA 6,509 → 6,156 → 5,706 → 5,000; ATAC 6,500 → 6,280 → 5,830 → 5,000; all 120 planted features `kept`).

## 2. Implementation on the synthetic data

- [ ] 2.1 Write `src/sead_sexdiff/stage_de_limma.R` (design D2, D3–D7, D11). It reads the exchange files and performs, in order:
  - checks that the counts' donor columns, the design rows and the export TMM donors are the same ids in the same order;
  - TMM and the check against the export factors (`STAGE_DE_LIMMA_TMM_ATOL`);
  - the `model.matrix` rebuild and its equality under the name map;
  - `voomWithQualityWeights`, `lmFit`, `eBayes`;
  - `topTable` per tested coefficient and per joint test;
  - the exact-F diagnostic.

  It writes the outputs with 17 significant digits, plus `fit.json`. A `--versions` mode reports R, limma and edgeR. Verify by hand, on exchange files written from the synthetic cohort_30 RNA `factor_interaction` design:
  - the TMM factors equal the export manifest's to ≤ 1e-8;
  - the design check passes, and fails naming the model when the `sex` levels are reversed in `design_table.csv`;
  - a factor perturbed by 3e-6 stops the script naming the omic;
  - every double in the outputs round-trips through Python's `float()` unchanged.
- [ ] 2.2 Implement `stage_de_limma_fit` (Python). It writes the exchange files, runs `Rscript --vanilla` with `OMP_NUM_THREADS = OPENBLAS_NUM_THREADS = 1`, checks the exit status, parses the results with round-trip precision, and returns the tables (with `exported`) and the fit (TMM check, design check, sample weights with per-cell means, df.prior, s2.prior, span, F diagnostic, versions, script sha256, runtime). Verify on both synthetic cohorts and omics, all four models:
  - the resolved design columns equal `stage-de`'s fit records;
  - the planted stage features dominate the smallest adj.P of the stage tests and of `stage_joint`;
  - the planted interaction features dominate the interaction tests and `interaction_joint`;
  - recall at adj.P ≤ 0.05 is recorded per test, with `interaction_joint` beside the union;
  - null features have p ≤ 0.05 at about the nominal rate;
  - sample weights stay near 1, since the synthetic quality is homogeneous;
  - the exact-F diagnostic's largest relative difference and changed-call count are recorded.

  Error paths: with `PATH` lacking `Rscript`, the error names the Nix shell; a corrupted `counts.csv` in a debug run gives an error ending with R's stderr; neither writes a table.
- [ ] 2.3 Calibration on synthetic data (verification only; design D4 and Risks). Rerun the earlier change's simulations: every condition × seed of `calibration_runs.csv`, regenerated with `cal_sim.py` / `cal_sim_syn.py` from their seeds. Fit each simulated dataset with limma-voom with quality weights (through `stage_de_limma_fit`) and with plain voom (a scratchpad R script, the packaged script with `voom` in place of `voomWithQualityWeights`; the packaged script gets no switch). Add the predeclared heterogeneous-quality condition: the `trend_rna` and `logunif_atacsyn` settings, the same seeds, and one donor per sex × stage cell whose dispersions are multiplied by 4. Fit it with pydeseq2 too.
  - Record per method and table: the realized FDR at α (pooled over seeds) and recall. Compare with DESeq2's figures from `calibration_pooled.csv` (0.06 with 200 / 4,000 non-null; 0.09 with 60 / 5,830).
  - Verify that the comparison table exists for every condition.
  - **Stop point:** if limma-voom with quality weights has a higher pooled realized FDR than DESeq2 on the homogeneous conditions' additive tables, stop and ask the user before group 3 (no real data read).
- [ ] 2.4 Implement `run_stage_de_limma`, in this order:
  - `_stage_de_prepare(..., log="stage-de-limma")` for every requested cohort;
  - the DESeq2 reference check (D10);
  - the fingerprints and the resume states (D12);
  - the fits (tables, then the record, through `.tmp` files);
  - the summary with union and joint rows;
  - both concordance cells and `methods_agree` (D8);
  - `method_comparison.csv` (D9);
  - the findings, the hash re-check of the inputs and of the DESeq2 outputs, and the manifest.

  Verify on copies of the synthetic inputs. Each of these stops naming the file and the producing command, and writes nothing:
  - missing `stage_de_manifest.json` (`stage-de`);
  - a DESeq2 manifest recording another `rna_counts.h5ad` hash (`stage-de --force`);
  - an edited DESeq2 table;
  - a DESeq2 summary whose count disagrees with its tables;
  - modified counts (`export`);
  - an out-dir inside an input or inside the DESeq2 directory.

  Resume: deleting one table refits only that model with a byte-identical summary and comparison; a rerun without `--force` skips all 16 fits; a record whose script sha256 differs is stale and stops without `--force`.

  Inputs and DESeq2 outputs hash unchanged after every run. After rerunning `stage-de --cohort all --omic both --force` with the changed package, every `stage-de` table and `stage_de_summary.csv` equals the 1.3 baseline byte for byte, and the fit records and manifest are equal apart from timestamps, runtimes and `git`. This shows that the DESeq2 path is unchanged.
- [ ] 2.5 Implement `render_stage_de_limma` (D13) and the `stage-de-limma` CLI subcommand (`--cohort`, `--omic`, `--pseudobulk-dir`, `--input-dir`, `--analysis-dir`, `--deseq2-dir`, `--out-dir`, `--force`). Document the command in README with a section "Pseudobulk stage DE, limma-voom (second method)": command, flags, method chain and settings, tests and labels, BH family, TMM check, concordance cells side by side, method comparison, outputs, `debug_dir` files and rerunning the R script by hand, reading rules. Document it in AGENTS.md:
  - pipeline step 12;
  - a "Stage-DE limma constants" bullet under flag-and-ask with the reading rules;
  - the CLI name;
  - `stage_de_limma.R` in Layout;
  - limma in the flake sentence.

  Verify:
  - `uv run sead-sexdiff stage-de-limma --help` lists the eight flags;
  - a synthetic `stage_de_limma.md` lists the reading rules, the concordance table, the count tables with the DESeq2 counts beside them, the comparison, the weights, the TMM check, the F diagnostic and the findings, and every count equals the summary (parser check);
  - every `STAGE_DE_LIMMA_*` name in README and AGENTS.md exists in `data.py`.

## 3. Real data

- [ ] 3.1 Run `uv run sead-sexdiff stage-de-limma --cohort all --omic both` on nixos-dell-laptop, which holds `results/pseudobulk/`, `results/motco_input/`, `results/motco_analysis/` (the MOTCO 0.9.0 primary) and `results/stage_de/`. Verify:
  - each cohort holds 16 tables, 8 fit records, the summary, `method_comparison.csv`, `stage_de_limma.md` and the manifest (`complete = true`);
  - every exported `n_tested` is 5,000, and every abundance `n_tested` equals `stage-de`'s (rna 20,387 / atac 128,916 in cohort_30; 20,238 / 124,818 in cohort_50);
  - the TMM check's largest difference is ≤ 1e-8, and every design check passed;
  - every input file and every `results/stage_de/` file hashes as before (`sha256sum -c`);
  - the exact-F changed-call counts and the per-cell sample-weight means are recorded here.
- [ ] 3.2 Rerun `cohort_30` with `--force` into a scratch directory. Verify that `stage_de_limma_summary.csv`, `method_comparison.csv` and all 16 tables are byte-identical to 3.1's.
- [ ] 3.3 Copy `results/stage_de_limma/` to the cluster (`/home1/tgonza/sead_sexdiff/results/`). Verify with `sha256sum -c` on the cluster.

## 4. Report and validation

- [ ] 4.1 Write `reports/motco-astrocyte-stage-de-limma-<date>.md` from the generated files with a script, with these sections:
  - provenance: inputs with hashes, R / limma / edgeR, the script sha256, git, machine;
  - the statement that the proposal was written after the DESeq2 results;
  - the method chain;
  - per cohort, the side-by-side concordance table, the count tables of both methods on both feature sets, the comparison rows, the features at adj.P ≤ 0.05 per test (or the five smallest), the sample weights per cell, the TMM check and the F diagnostic;
  - the synthetic calibration (2.3), labelled as simulations under DESeq2's model;
  - the reading rules and limitations the spec requires;
  - next steps, naming the microglia amendment.

  Verify with a check script that every count, concordance cell and comparison value equals `stage_de_limma_summary.csv` or `method_comparison.csv`, every DESeq2 number equals `stage_de_summary.csv`, and the `delta` statuses equal the primary `result.json` records.
- [ ] 4.2 Run `openspec validate mtg-astrocyte-stage-de-limma` and verify that it passes.

## Workflow follow-up

- Commit in the project's pattern (proposal; implementation with the R script; report), open the PR. The user merges.
- Archive the change and sync the `pseudobulk-stage-de-limma` spec.
- After the merge, amend `mtg-microglia-positive-control` (branch `origin/mtg-microglia-positive-control`, rebased onto main, via `/opsx:update`) before any microglia run:
  - `--cell-type` and the `cell_type` manifest check on `stage-de-limma` (eight MTG commands);
  - `stage-de-limma` in its D10 byte-identity baseline;
  - the limma concordance rule for microglia beside the DESeq2 one;
  - a predeclared choice for E1 (DESeq2 `factor_additive` union as written, limma `stage_joint`, or both) and how a method disagreement enters the all-three rule.
