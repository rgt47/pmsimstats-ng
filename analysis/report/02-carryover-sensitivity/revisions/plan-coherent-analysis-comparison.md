# Plan: A Coherent Experiment Comparing Analysis Strategies (Paper 02)

*2026-09-08 12:04 PDT*

## Purpose

Paper 02 (`analysis/report/02-carryover-sensitivity/`) compares nine
analysis specifications (G1-G9) for the biomarker-treatment
interaction test under Architecture B. This document records an audit
of the manuscript against its own simulation code, and a plan to
restructure the study so that each reported claim is supported by an
experiment designed to answer the question that claim makes.

The audit was conducted 2026-09-08 by reading `report.Rmd` and
`supplement.Rmd` in full, inspecting the heatmap and simulation
drivers, and querying the saved `.rds` outputs directly.

## Audit findings

Findings are ordered by severity. Each is labeled with its epistemic
status: `verified` (executed and observed), `inspected` (read in
source and confirmed), or `inferred`.

### 1. The decay-shape axis does not measure mis-specification

**Status: inspected, with a verified supporting calculation.**

Neither `23-run-decay-shape-sensitivity.R` nor
`25-run-decay-shape-sensitivity-g9.R` sets `analysis_form`,
`analysis_shape`, or `analysis_t1half` in its grid. Confirmed by
reading the saved `$grid` object from both `.rds` files: the columns
are absent.

`simulation-core.R:815-817` (in `simulate_cell_s7`, and identically
at 417-419 in `simulate_cell`) defaults each analysis-side parameter
to the corresponding DGP value:

```r
analysis_form <- default_val(cell[['analysis_form']],
                             cell$carryover_form)
```

`prepare_long_data()` then builds the `Dbc` predictor from
`analysis_form` and `analysis_shape`. The exposure-weighted predictor
was therefore constructed with the *correct* Weibull shape in all 30
cells of both runs.

The two predictors differ materially. At a true half-life of one week
and shape `k = 0.25` (verified by calling `carryover_decay()`
directly):

| `tsd` (weeks) | 0.5 | 1.0 | 2.5 | 5.0 |
|---|---|---|---|---|
| matched, `k = 0.25` | 0.558 | 0.500 | 0.418 | 0.355 |
| assumed exponential | 0.707 | 0.500 | 0.177 | 0.031 |

They agree only at `tsd = t_half`, by construction, and diverge by an
order of magnitude in the tail.

The manuscript describes this axis as mis-specification in at least
five places: the Abstract ('widens under a light-tailed decay-shape
mis-specification'), Section 3.5 ('the cost of assuming exponential
decay when the truth is Weibull'), Section 4.3 ('Under incorrect
decay-shape assumptions'), the Section 4.5 decision rule, and
Conclusion 2. The figure filename
(`02xs-heatmap-matched-vs-mismatched.pdf`) carries the same reading.

The numbers themselves are sound. What is measured is how power
varies with the true decay shape when the analyst specifies it
correctly, which is a different and less interesting finding than the
one reported. A heavy-tailed DGP leaves less recoverable signal
regardless of what the analyst assumes.

Section 3.5 additionally refers to 'the remaining Exposure-weighted
row' when Exposure-weighted is a column in that figure, suggesting the
paragraph was written against an earlier figure orientation.

### 2. The headline result is adjudicated against the wrong yardstick

**Status: inspected.**

Conclusion 8 calls G9 beating G8 'the single best result in this
manuscript'. The margin is 0.007 (0.870 against 0.863). G8 was not
rerun for Block S10; it is carried over from Block S6 on a different
seed. That single contrast is therefore unpaired, while the
manuscript's stated justification for reading small gaps is that all
specifications share common random numbers.

The supplement discloses this twice, including in the Table S10
caption. The main text's Section 3.7 and Conclusion 8 do not.

### 3. A Monte Carlo noise floor is derived and then ignored

**Status: inspected.**

Supplement Section S2 establishes a 'roughly four- to five-point noise
floor' from the scatter of specifications that are invariant to the
axis being varied, and uses it to explain away an apparent anomaly.
Section S8 then reports 2.2 to 2.6 point reversals as 'a real, if
modest, reversal rather than pure noise', and Section S10 treats 0.7
points as decisive. The three passages cannot all be right.

### 4. Cell count is internally inconsistent

**Status: inspected.**

The Abstract and Sections 3.4 and 4.1 state 216 cells. Section 2.6's
own arithmetic gives `3 x 3 x 3 x 2 x 3 = 162`. Section 2.9 describes
filtering a 540-cell grid 'down to the 216-cell slice'. The two
figures are used interchangeably without reconciliation.

### 5. The oracle caveat is not propagated

**Status: inspected.**

Section 4.3 correctly notes that outside Block S2 the analysis-side
half-life equals the data-generating one, granting Exposure-weighted
an oracle no analyst possesses, and that Tier 1 margins should
therefore be read as upper bounds. The Abstract and Conclusion 2 both
quote the ten-point advantage without that qualification.

### 6. Lower-severity items

**Status: verified unless noted.**

- Supplement cross-references into main-text Section 3 are uniformly
  off by one (S1 and S2 cite 3.6 for material in 3.7; S10 cites 3.5
  for Block S6, which is 3.6; S8 cites 3.3 for the null in 3.4).
  Inspected.
- The supplement states that 'Block, specification, and G-code
  notation match the main manuscript exactly', but Sections S1, S2,
  and all of S8 are written in E-codes, which the main text reserves
  for storage. Inspected.
- `dgp_arch` is absent from the summary and `meta` of
  `02-grid-summary-hendrickson-g9.rds` and
  `02-decay-shape-sensitivity-g9.rds`. The architecture is fixed only
  in the driver source. The Architecture A counterparts write
  identically shaped files, so the two cannot be distinguished from
  the artifact alone.
- `meta$elapsed_secs` for the 2026-09-01 run records 886 seconds,
  against 57 minutes spanned by the checkpoint mtimes.
- `13-hendrickson-heatmap-tier2.R` describes the G1-G9 extension
  panels as 'Preliminary n_sim = 100 (medium smoke test)' at two
  places. All three G9 files are `n = 500` in every cell.
- Six heatmap PDFs in `analysis/figures/` are built but referenced by
  neither document: `hendrickson-d`, `sens-S6-cr2`, `sens-S6-model`,
  `sens-S6-g9-null`, `sens-S6-g9-power`, and a stale `sens-S8.pdf`.
- `report.pdf` (2026-09-01) is older than `report.Rmd` (2026-09-07).

## The underlying problem

The manuscript conflates three questions that require different
experimental structures.

- **Q1, oracle ranking.** Which specification performs best when the
  analyst knows the truth? Requires the analysis side set to truth.
- **Q2, cost of error.** What does an incorrect assumption cost?
  Requires the analysis side crossed with truth as independent
  factors.
- **Q3, estimation against assumption.** Does AIC selection recover
  what a correct assumption would give? Requires both.

Tier 1 is constructed for Q1 and is reported as answering Q2 on the
decay-shape axis. Block S2 answers Q2 for half-life alone, with three
specifications. Q3 is answered at a single cell. This is why the same
result reads as robust in one section and fragile in another: the
sections are not measuring the same quantity.

A second structural problem is that Tier 2 accreted rather than being
designed. Nine blocks across seven drivers and three seeds, with
inconsistent specification coverage, means no two blocks are strictly
comparable, and the manuscript's strongest claim rests on the least
comparable of them.

## Plan

### Phase 0: Decide the decay-shape fork

This decision gates Phase 2.

- **Re-run with the analysis side pinned to exponential.** Makes the
  existing prose true. Approximately 15 to 60 minutes given per-cell
  checkpointing.
- **Rewrite the prose to describe sensitivity under correct
  specification.** Free, touches five passages, but leaves the
  manuscript with no decay-shape mis-specification evidence at all.

**Recommendation: do both, as separate factors.** They answer
different questions and the paper needs both. This is absorbed into
Phase 2 below.

### Phase 1: Consolidate the simulation

To be completed before any new simulation is run. One driver, one
seed, one output schema.

1. A single `simulate_cell()` entry point fitting all nine
   specifications to every generated dataset, so that every contrast
   in the paper is paired, G8 included.
2. Make `analysis_t1half`, `analysis_form`, and `analysis_shape`
   required grid columns. The current silent fallback to DGP truth is
   the direct cause of finding 1 and should raise an error rather
   than default.
3. Record `dgp_arch`, `n_reps`, the seed, and both DGP and
   analysis-side parameters in `meta`, and retain them as summary
   columns even when constant.
4. Compute `elapsed_secs` with `difftime(units = 'secs')`.

This phase costs roughly a day and structurally removes findings 1,
2, and part of 6.

### Phase 2: Two clean studies

**Study A, oracle ranking (Q1).** Analysis set to truth throughout,
and stated as such. Design (3) x `N` (2) x `c_bm` (3) x `t_half` (3)
x decay form (5) x 9 specifications, approximately 270 cells.
Absorbs the current Panels A, B, and C into one run under one seed.

**Study B, mis-specification (Q2 and Q3).** The genuinely new work.
DGP truth crossed against analyst assumption as independent factors:

- true shape in {exponential, 0.25, 0.5, 2.0, 4.0} at `t_half = 1.0`
- assumed shape in {exponential, 0.5, 2.0}
- assumed `t_half` in {0.5, 1.0, 2.0}
- designs {CO, Hybrid}, `N = 70`

This gives 90 cells and subsumes Block S2 as the half-life slice of a
single coherent factorial rather than a separate block. G4 and G9
belong here, since AIC selection is precisely the strategy that gets
to see the data, so Q3 is answered against a real baseline rather
than an oracle one.

The existing S1, S3, and S4 robustness runs are retained. They are
already at 500 replicates with all nine specifications.

### Phase 3: Fix the inferential reporting

Independent of the simulation work, and startable immediately. This
phase uses existing data only.

1. **Report paired contrasts with paired standard errors.** The
   machinery already exists: Section 3.8 applies McNemar correctly to
   the 24-cell rerun. Every headline gap should be reported as a
   paired difference with its own standard error, rather than as two
   point estimates compared against a marginal MCSE of 0.022. This
   will settle whether the 0.007, 2.4, and 10 point gaps are each
   real, without any new simulation.
2. Adopt one MCSE convention and apply it throughout, including the
   supplement's noise-floor paragraph (finding 3).
3. Reconcile 162 against 216 and state the counting rule once
   (finding 4).

### Phase 4: Manuscript restructure

Gated on Phase 2.

1. Retitle the decay-shape material to match what it measures, and
   add Study B's results as a distinct mis-specification section.
2. Propagate the oracle caveat to the Abstract and Conclusion 2
   (finding 5).
3. Demote Conclusion 8 to match its evidence, or promote the evidence
   via Study B.
4. Correct the Section 3.5 row and column orientation error.

### Phase 5: Housekeeping

1. Correct the supplement's Section 3 cross-references.
2. Sweep E-codes to G-codes in supplement Sections S1, S2, and S8, or
   withdraw the claim that notation matches.
3. Delete or regenerate the six orphaned figure PDFs.
4. Update the two stale `n_sim = 100` script comments.
5. Rebuild `report.pdf` and `supplement.pdf` through
   `bash tools/render.sh`.

## Sequencing

Phases 3 and 5 are independent of the simulation work and can begin
immediately. Phase 1 gates Phase 2, which gates Phase 4. Phase 0 is a
decision, not work, and should be taken before Phase 2 is scoped.

The highest-value immediate action is Phase 3, item 1. It is free, it
uses data already on disk, and it determines which of the
manuscript's disputed gaps survive before any compute is spent
deciding what to re-run.

## Provenance

Audit conducted 2026-09-08 against `report.Rmd` at its 2026-09-07
state and `supplement.Rmd` at its 2026-08-19 state. Simulation
outputs inspected: `02-grid-summary.rds`,
`02-grid-summary-hendrickson-g9.rds`,
`02-decay-shape-sensitivity.rds`,
`02-decay-shape-sensitivity-g9.rds`,
`02-sensitivity-summary-g9.rds`. Drivers inspected:
`19-run-tier1-hendrickson-g9.R`, `23-run-decay-shape-sensitivity.R`,
`25-run-decay-shape-sensitivity-g9.R`, `simulation-core.R`,
`12-hendrickson-heatmap.R`, `13-hendrickson-heatmap-tier2.R`,
`03-render-figures-extra-slim.R`.

Not verified: no cell was re-run under an alternative analysis-side
specification, so the magnitude of finding 1's effect on reported
power is unquantified. Only its existence is established.
