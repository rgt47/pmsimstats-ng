---
title: 'Provenance of the response and baseline parameters (referee comment M1)'
author: 'pmsimstats team'
---

*2026-11-02 10:44 PST; updated 2026-11-02 16:57 PST with the origin
of the `simulation-core.R` values*

## Question

Referee comment M1 (`referee-report-2026-11-01.md`) found that the
G10 planning results and the G10 Figure 4B panel in paper 01 were
generated under response parameters different from those paper 01
specifies. This note establishes how and when the two parameter sets
arose, and whether the departure is a mistake or a deliberate fix.

## Finding

The repository has carried two parameter sets in parallel since its
first commit. The set in `simulation-core.R` is the departure from
the published process. Its values are specific because they descend
from a synthetic unit-test fixture written in March 2026 to check a
code port, not from any calibration. The fixture was reused as the
default parameter set of the early carryover drivers. Two of its
standard deviations were then edited and the result was labeled as
Hendrickson's defaults. It is a provenance error, not a fix.

## Which set is Hendrickson's

The published package data (`data/extracted_rp.rda` and
`data/extracted_bp.rda` at commit 58b32a9 of
`rchendrickson/pmsimstats`) were downloaded and loaded in R on
2026-11-02 (verified).

| Parameter | Hendrickson 58b32a9 | Paper 01 `build()` | `simulation-core.R` |
|---|---|---|---|
| TV curve (max, disp, rate) | 6.506, 5, 0.35 | same | 10.986, 5, 0.42 |
| TV SD | 10 | 10 | 6 |
| PB SD | $10\eta_t$ | $10\eta_t$ | $2\eta_t$ |
| BR SD | 8 | 8 | 8 |
| Biomarker mean, SD | 124.3, 15.36 | same | 0, 1 |
| Baseline mean, SD | 83.07, 18.48 | same | 70, 10 |

Paper 01's closed-form DGP
(`hendrickson-problems/04-power-simulation.R`) reproduces the
published values exactly. `simulation-core.R`
(`default_resp_param()`, `default_baseline_param()`) is used by paper
02 and by driver 19 (`19-g9-planning-and-fig4b.R`), which produced
the G10 results in paper 01.

## How and when

- **In the repository.** `git log -L` shows that both default
  functions entered the repository in the initial commit
  (`bbed702e9`, 2026-05-01, a bulk import of 201 files) and have not
  been edited since. `simulation-grid-plan.md`, committed with them,
  is dated 2026-04-15 and already describes the driver as built, so
  the file predates the repository.
- **Earlier history.** The pre-repository material is in the outer
  workspace at `~/prj/res/36-pmsimstats-ng/archive/` (not the Dropbox
  path named in `CLAUDE.md`). The archived 2025 repository has a
  single commit dated 2026-05-01, so the dates below are file
  modification times.

| Stage | File | TV curve | SDs TV, PB, BR | Biomarker; baseline |
|---|---|---|---|---|
| Published | `extracted_rp.rda`, `extracted_bp.rda` at 58b32a9 | 6.506, 0.35 | 10, $10\eta_t$, 8 | 124.3/15.4; 83.1/18.5 |
| 1. Test fixture (2026-03-25) | `archive/2025/02-test-alignment.R` | 10.986, 0.42 | 5, 2, 5 | 0/1; 70/10 |
| 2. Early carryover drivers (2026-04-22) | `archive/carryover_factorial/01-run-factorial-orig.R` | 10.986, 0.42 | 5, 2, 5 | 0/1; 70/10 |
| 3. Paper 02 driver (by 2026-04-15) | `simulation-core.R` | 10.986, 0.42 | 6, 2, 8 | 0/1; 70/10 |

The stage 2 file date is a later modification; the drivers must
have existed before stage 3 for the lineage to hold (inferred).

- **Stage 1, the fixture.** `02-test-alignment.R` checked that the
  2025 tidyverse port reproduces the original package's Gompertz
  values, covariance matrix and generated data exactly. An equality
  test needs only some valid inputs, so the fixture used simple
  ones: a standardized biomarker, baseline 70/10, TV on the BR curve
  and SDs 5, 2, 5. The same fixture remains in `inst/tinytest/` and
  `implementations/original/tests/`.
- **Stage 2, reuse as defaults.** The `carryover_factorial` drivers,
  the precursors of paper 02, copied the fixture verbatim as their
  'default' parameter set.
- **Stage 3, two edits.** `simulation-core.R` raised the BR SD from
  5 to 8 and the TV SD from 5 to 6. Its comment reads: 'br sd = 8
  matches the docs/19 calibration target ... tv / pb values follow
  the Hendrickson (2020) defaults; tv sd is bumped to 6 so the TV
  trajectory is not dominated by BR noise.' The BR change has a
  stated reason: `docs/19` calibrated mean moderation against the
  published parameters, whose BR SD is 8. 'Bumped to 6' is accurate
  as an edit from 5. The claim that the TV and PB values follow
  Hendrickson is false: TV uses the BR curve, and PB has SD 2, not
  10.
- **Strength of the inference.** The Gompertz values agree to six
  significant figures, PB SD 2 and the baseline values carry over
  unchanged, and 'bumped to 6' is coherent only as a change from 5.
  The lineage is therefore well supported, but it is inferred: no
  commit records the copy.
- **The published values were at hand.** The 2025 codebase that the
  fixture tested defined them explicitly
  (`pmsimstats2025/analysis/scripts/pm_functions.R`, lines
  1376-1385, `hendrickson_resp_params` and `hendrickson_bl_params`),
  with SDs 10, 10, 8, biomarker 124.3/15.4 and baseline 83.1/18.5.
  The drivers took the fixture instead, presumably because it was
  the parameter block that already ran end to end (inferred).
- **Detection.** `docs/25` recorded the departure on 2026-10-07 and
  noted that an earlier version of that document had listed the
  `simulation-core.R` values as the Hendrickson defaults.
  `analysis/scripts/csp563-calibrated/00-plan.md` noted the mismatch
  again and left open which set is Hendrickson's. The download above
  settles that question.
- **Why it stayed hidden.** Paper 02's data-generating section says
  it uses 'the extended DGP of the companion study' and does not
  state its response parameters.
- **Not a CSP 563 calibration.** The defaults predate the CSP 563
  reanalysis (`docs/69`, 2026-10-30) by six months. Of the values
  that differ from the published set, only the BR SD of 8 has a
  rationale, and it equals the published value. The TV curve, TV SD
  6 and PB SD 2 are inherited from a test fixture and have no
  scientific basis.

## Which differences matter

- **Biomarker mean and SD (0, 1).** Harmless (inferred). The test of
  the biomarker coefficient is invariant to linear rescaling of the
  biomarker.
- **Baseline mean (70).** Harmless (inferred). It cancels from change
  scores and is absorbed by the intercept.
- **Baseline SD (18.5 to 10).** Probably minor (inferred). Baseline
  is constant within a patient, so the random intercept absorbs most
  of it, and it cancels exactly in the paired-difference statistic
  E9.
- **TV SD (10 to 6) and PB SD ($10\eta_t$ to $2\eta_t$).** Material.
  The SD of change from baseline, computed by hand with the
  within-occasion cross-component correlation $c_1 = 0.2$:

| Source | SD of change from baseline |
|---|---|
| Hendrickson parameters | 16.0 ($\eta_t = 0.5$) to 19.2 ($\eta_t = 1$) |
| `simulation-core.R` parameters | about 11.2 to 11.6 |
| CSP 563 observed (`docs/69`) | 17 to 24, increasing over time |

The 16.0 agrees with the value `docs/69` reports for the DGP, which
checks the arithmetic. The empirical data therefore favor the
published parameters, which if anything are slightly too quiet; the
`simulation-core.R` set is further from the data. The lower noise is
what inflates the G10 advantage over E9 in paper 01: under the
`simulation-core.R` generator, E9 power at Hybrid, $N = 70$,
$c_{bm} = 0.25$, $t_{1/2} = 0.2$ week is 0.549 (Monte Carlo, 1,000
replicates), against 0.285 in closed form under the published
parameters.

## Consequences and options

- **Paper 01.** Rerun driver 19 with the `build()` parameters. This
  resolves M1 directly and is a single driver.
- **Paper 02.** All results were produced under the
  `simulation-core.R` set. Switching to the published values would
  lower power throughout and require rerunning the ladder. The
  cheaper alternative is to keep the runs, state the parameter table
  explicitly in paper 02, and remove the claim that it uses paper
  01's DGP unchanged.
- **Code.** Correct the comment in `simulation-core.R` in either
  case, and add a test that the generator reproduces the closed-form
  mean and variance of E9 at a few cells.

## Evidence status

- Verified: published parameter values (downloaded and loaded); the
  parameters in both generators (read from code); the git history of
  the default functions; the E9 Monte Carlo value; the fixture values
  in `archive/2025/02-test-alignment.R`, the archived
  `carryover_factorial` drivers and the package tests; the
  `hendrickson_resp_params` values in the 2025 codebase.
- Computed by hand: SDs of change from baseline.
- Inferred: the lineage from fixture to `simulation-core.R` (no
  commit records it); the irrelevance of biomarker and baseline
  scaling to power; the reason the drivers used the fixture.
- Not checked: exact creation dates before 2026-05-01; file
  modification times are the only record.
