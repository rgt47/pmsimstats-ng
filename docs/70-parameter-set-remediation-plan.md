---
title: 'Restoring a single response-parameter set across the compendium'
author: 'pmsimstats team'
---

*2026-11-02 17:07 PST*

## 1. Purpose

Referee comment M1 on paper 01 (`analysis/report/01-dgp-mean-
moderation-vs-mvn/referee-report-2026-11-01.md`) found that paper
01's G10 results were simulated under response parameters different
from those the paper specifies. Tracing the cause
(`analysis/report/01-dgp-mean-moderation-vs-mvn/m1-parameter-
provenance-2026-11-02.md`) showed that the compendium uses three
parameter sets, two of which descend from a synthetic unit-test
fixture and have no scientific basis. This document plans the return
of every reported result to one documented parameter set.

## 2. The three parameter sets

| Set | TV curve (max, rate) | SD TV, PB, BR | Biomarker; baseline | SD of change from baseline |
|---|---|---|---|---|
| Published (Hendrickson 58b32a9, `data/extracted_rp.rda`, `extracted_bp.rda`) | 6.506, 0.35 | 10, $10\eta_t$, 8 | 124.3/15.4; 83.1/18.5 | 16.1 to 19.1 (simulated) |
| Sim-core (`simulation-core.R` `default_*()`) | 10.986, 0.42 | 6, $2\eta_t$, 8 | 0/1; 70/10 | about 11.5 (by hand) |
| Raw fixture (test fixture of 2026-03-25) | 10.986, 0.42 | 5, $2\eta_t$, 5 | 0/1; 70/10 | about 8.5 (by hand) |
| CSP 563 observed (`docs/69`) | not identifiable | not identifiable | 14.2; 16.0 | 17 to 24, increasing |

With the published values passed in, the `simulation-core.R`
generator reproduces the closed form exactly (E9 power 0.280, MCSE
0.014, against 0.285; SD of change 19.1 and 16.1 against 19.2 and
16.0). Parameters are the only difference between the generators.

## 3. Decision required before any rerun

**Which set is the target?** Every rerun below is expensive, so the
target should be fixed once.

- **Option A (recommended): the published set.** It is the set
  paper 01 specifies and the closed form uses, it is traceable to a
  public source, and it is within or just below the CSP 563 range on
  total variance. It becomes the parameter table of the refined
  process.
- **Option B: a CSP 563-calibrated set.** It matches the observed
  total variance and its growth, but the TV/PB split, the expectancy
  weight and the carryover half-life are not identifiable from CSP
  563 (`docs/69`), so the calibration would still need assumptions.
  The working protocol `analysis/scripts/csp563-calibrated/00-plan.md`
  already treats it as a sensitivity arm.

The plan below assumes Option A, with the CSP 563 set retained as a
sensitivity arm for paper 01 and paper 02 only.

## 4. Inventory of affected results

Driver parameter sources were read from code on 2026-11-02
(inspected). Whether each paper's claims are relative (comparisons
between analyses or designs on common data) or absolute (power,
sample size, bias magnitudes) is a first reading and is confirmed in
Phase 0.

| Paper | Drivers behind reported results | Set | Claims (first reading) | Action |
|---|---|---|---|---|
| 01 DGP | `hendrickson-problems/04-power-simulation.R` (closed form, study 1), `11-strawman-stress-test.R`, `figure4/` | Published | Absolute | None |
| 01 DGP | `19-g9-planning-and-fig4b.R` (G10 planning, Fig 4B G10) | Sim-core | Absolute | Rerun (Phase 1) |
| 01 DGP | `26-co-dbc-fix-paired.R` (Appendix, CO guard rerun) | Sim-core | Relative | Rerun (Phase 1, small) |
| 02 Carryover | `01-run-factorial.R`, `06`, `27`, `32`, `36`, `38`, `41`, `cr2-olbdc-check` | Sim-core | Mixed: rankings relative, size and power absolute | Rerun (Phase 2) |
| 03 Latent class | `10-study5-production.R`, `09-calibration-substudy.R` | Raw fixture | Mixed | Rerun (Phase 3) |
| 03 Latent class | `12-referee-rerun.R` | Sim-core | Mixed | Rerun (Phase 3) |
| 04 Main effect | `treatment-main-effect/production-run.R` (`nof1power`) | Own main-effect parameters | Absolute | Out of scope (Section 5) |
| 05 Design sensitivity | `nof1-design-sensitivity/` (`baseline_resp_param()`, unit scale) | Own, deliberate | Relative | Out of scope (Section 5) |
| 06 Component decomposition | `component-decomposition/run-study-*.R`, `00-common.R` | Raw fixture | Mixed | Rerun (Phase 3) |
| 07 Gompertz evaluation | `gompertz-evaluation/0[1-3]-*.R` | Published | | None |
| 08 Test procedure | `01-study1-test-procedure.R` | Raw fixture | Mixed | Rerun (Phase 3) |
| 09 Dropout | `20-paired-dropout.R` (BR SD varied as `SIGMA_BR`), `00-common.R`, `cr2-olbdc-check` | Raw fixture, sim-core | Mixed | Rerun (Phase 3) |
| 10 Test calibration | `paper10-common.R`, `run-paper10-refined.sh` | Sim-core | Absolute (size) | Rerun (Phase 2, with paper 02) |
| 11 Combined DGP | `dgp-combined/01-run-combined-factorial.R` | Published | | None |
| 12 Design efficiency | `dgp-combined` summary, `01-closed-form-validation.R` | Published | | None (confirm the closed-form validation in Phase 0) |
| 13, 14 Two-component | `two-component-dgp/01-run-grid.R` | Published | | None |

**Results that do not depend on the response SDs.** The feasibility
bounds (Schur complement, Kronecker factor) and all ceilings depend
only on the correlation matrix, so the compendium's ceiling results,
`docs/34`, `docs/35`, `docs/65` and the `compendium-mathematics`
ceiling checks are unaffected. The checks in `compendium-mathematics/`
that compute power already use the published values.

## 5. Out of scope, with reasons

- **Paper 04** studies the treatment main effect with the vendored
  `nof1power` engine and its own parameterization. It does not claim
  to use the Hendrickson process. Confirm in Phase 0 that its text
  makes no such claim.
- **Paper 05** uses a deliberately unit-scaled parameter set and
  varies the BR SD as a design axis (`00-designs.R` labels its design
  'not the published Hybrid'). Its claims are about relative design
  performance. Add one sentence stating its parameter table.
- **Unit and parity tests** (`inst/tinytest/`, `implementations/*/
  tests/`, `analysis/scripts/parity/`) may keep the fixture: an
  equality test needs only valid inputs. Rename the fixture
  `test_fixture_resp_param` in a comment so it is not mistaken for a
  default again.
- **Prototype and quick-sim drivers** (`quick-sim/0*-prototype*.R`,
  `quick-sim-driver*.R`, smoke tests) that feed no reported number
  are left as they are.

## 6. Phases

### Phase 0: infrastructure and triage (no simulation)

1. Add `published_resp_param()` and `published_baseline_param()`,
   built from `data/extracted_rp.rda` and `data/extracted_bp.rda` (not
   retyped), to the package (`R/utilities.R`) and export them. Have
   `simulation-core.R` call them.
2. Rename `default_resp_param()` and `default_baseline_param()` to
   `simcore2026_resp_param()` and `simcore2026_baseline_param()`, and
   keep them only so stored runs reproduce; correct the comment that
   calls them Hendrickson's defaults. Make the published functions the
   default in every driver through a `--params published|simcore2026|
   fixture` flag, so each output records the set it used.
3. Add a moment check (`analysis/scripts/quick-sim/carryover-closed-
   form/24-check-generator-moments.R`) asserting that each generator
   (package `generateData`, `implementations/tidyverse`,
   `simulation-core.R`) reproduces the closed-form E9 power within
   2 MCSE and the SD of change from baseline within 2 percent at three
   cells. Run it in the tinytest suite at a small replicate count.
4. Record the parameter table in NOTATION.md, in CLAUDE.md under
   'Named processes', and in paper 01 Table `tab:spec`, so the refined
   process is defined by its correlation structure and its response
   parameters together.
5. For each paper marked Rerun, list every reported number, table and
   figure and classify it as relative or absolute. Estimate the
   compute cost of each rerun from the stored job grids and per-fit
   timings (most stored jobs do not record elapsed time). Record
   decisions in an appendix to this document.
6. Confirm the 'None' and 'Out of scope' verdicts for papers 04, 05,
   07 and 11 to 14 by reading their Methods text for any claim of
   matching paper 01 or the Hendrickson parameters.

### Phase 1: paper 01 (one simrun job and one small local run)

1. Rerun driver 19, blocks P and F, under the published set, fitting
   G10 and G9 on the same data so the paired comparison
   (`pg_vs9_max`) survives. Add a G10 $N_{80}$ bisection for the
   Hybrid cells (referee comment M5).
2. Rerun `26-co-dbc-fix-paired.R` under the published set.
3. Rewrite the planning results, Table `tab:planningg10`, the G10
   Figure 4B paragraph, the appendix paragraph on the crossover guard,
   the Limitations bullet and the Conclusions sentence. Expect the
   claim that the closed form is conservative for G10 to weaken or
   reverse.
4. Address the other referee comments that need no simulation (M2 to
   M4, M6 to M8, minor comments) in the same revision.

### Phase 2: papers 02 and 10 (the bulk of the compute)

Paper 02's ladder rests on many paired jobs (31 to 42 and 56 and the
G10 backfill). Rerun in dependency order so the pairing guards hold:

1. The base factorial (`01-run-factorial.R`) and Tier 1 G-vs-F-E0
   (`27`).
2. Blocks A to F (`32-run-g9-gap-fill.R`) and the referee variants
   (`36`), with G9 and G10 fitted on the same data in each job, which
   removes the need for a separate G10 backfill (`41`).
3. The E0 isolation (`38`) and the CR2 calibration checks
   (`cr2-olbdc-check`, `diagnose-null-type1-cr2.R`).
4. Paper 10 (`paper10-common.R`, `run-paper10-refined.sh`), which
   shares the generator.
5. Regenerate the ladder outputs (`37-ladder-tables-figures.R`), the
   summary (`build-summary.R`) and every typed number in both
   reports; replace typed numbers with values read from stored
   output wherever practical.

All arms that are paired trial by trial must run on the same host
(container and local results are not draw-for-draw identical).

### Phase 3: papers 03, 06, 08 and 09 (raw-fixture papers)

These papers used BR SD 5 rather than 8, so under mean moderation the
effect per biomarker SD ($c_{bm}\sigma_{BR}$) changes as well as the
noise. For each:

1. Switch the drivers to `published_*()` through the flag.
2. Check that the $c_{bm}$ grid stays feasible and meaningful at the
   new scale (the feasibility bound does not change; effect sizes on
   the outcome scale do).
3. Rerun the production and calibration studies; paper 09's
   `SIGMA_BR` axis should be re-centered on 8.
4. Re-narrate and re-render.

The order within Phase 3 follows the readiness of each paper for
submission, since none blocks another.

### Phase 4: secondary documents

1. Update the whitepaper series (`docs/series/`) and numbered docs
   that quote numbers from the rerun papers; `grep` for each changed
   number.
2. Update `docs/25` (Gompertz table) and the CSP 563 working protocol
   (`csp563-calibrated/00-plan.md`), whose Section 2 already notes the
   mismatch.
3. Update the provenance note and this plan with the outcome.

### Phase 5: verification

1. The moment check passes for every generator.
2. A search of `analysis/scripts/` for hard-coded `sd = c(` blocks in
   drivers that feed reported results returns none outside the
   fixture and the deliberate paper 05 set.
3. Every reported result records its parameter set in its output
   metadata.
4. All affected manuscripts re-render through `tools/render.sh`, and
   the notation lint passes.

## 7. Order and dependencies

Phase 0 precedes everything. Phase 1 is independent of Phases 2 and 3
and is the critical path for paper 01. Phase 2 is the largest and
should start as soon as Phase 0 items 1 to 3 are done, running on the
simulation hosts while Phase 1 text work proceeds. Phase 3 can share
the second host. Phase 4 follows each paper's rerun.

## 8. Risks

- **Conclusions may change.** Under 40 to 70 percent more noise SD,
  power falls everywhere. Rankings among analyses will probably
  survive (inferred), but size results under misspecification,
  convergence rates and the boundary between adequate and inadequate
  power can move. Each paper's rewrite must be driven by the new
  numbers, not by fitting them to existing text.
- **Feasibility interacts with the biomarker scale only through the
  correlation matrix**, so no grid cell becomes infeasible; but mean-
  moderation cells calibrated as $c_{bm}\sigma_{BR}$ change effect
  size in the raw-fixture papers.
- **Compute.** Paper 02 alone was produced by more than a dozen
  simrun jobs. Phase 0 item 5 must give a credible estimate before
  Phase 2 is committed; if it is prohibitive, paper 02 falls back to
  disclosure (stating its parameter table and dropping the claim that
  it uses paper 01's process unchanged).
- **History before 2026-05-01** is reconstructed from file
  modification times; the lineage of the fixture is inferred.

## 9. Evidence status

- Verified: the published values (downloaded from 58b32a9 and
  loaded); generator equivalence under the published values (Monte
  Carlo, 1,000 replicates); the parameter source of each driver named
  in Section 4 (read from code).
- Computed by hand: SDs of change for the sim-core and fixture sets.
- Inferred: the relative or absolute character of each paper's
  claims; that rankings survive the change; the fixture lineage.
- Not checked: compute cost of each rerun; papers 04 and 05 Methods
  text; `nof1-design-efficiency/01-closed-form-validation.R`
  parameters.
