# Plan: Architecture C (Combined DGP) for Paper 01
*2026-06-10 11:19 PDT*

## Motivation

Paper 01 compares Architecture A (mean moderation) and Architecture B
(MVN differential correlation) as mutually exclusive DGPs for
encoding biomarker-treatment interactions. The paper notes that under
carryover, Architecture B is attacked via two channels (mean blurring
and correlation erosion) while Architecture A faces only one. A natural
extension is a third DGP that activates both channels simultaneously.

## What Architecture C is

Architecture C applies both mechanisms in sequence:

1. The MVN covariance matrix encodes BM-BR differential correlation
   with weight `c.bm_b` (the Architecture B channel).
2. After the MVN draw, BR is additively shifted on-drug by
   `c.bm_a * bm_z * br_sd` (the Architecture A channel).

Special cases: `c.bm_b = 0` recovers pure Architecture A; `c.bm_a = 0`
recovers pure Architecture B.

## Parameterization

Two free parameters (`c.bm_a`, `c.bm_b`), not a mixing weight.
A single weight would impose false commensurability because the two
mechanisms act on different moments of the joint distribution.

Primary factorial grid: `c.bm_a` x `c.bm_b` in {0, 0.22, 0.45}^2
(nine cells, 1000 reps each). The (0.45, 0) and (0, 0.45) cells
reproduce existing Architecture A and B baselines exactly and serve
as the validation gate.

Note on calibration: with both parameters at 0.45 the total
interaction effect is approximately double that of either pure
architecture alone. Section 2.3 should acknowledge this and present
an optional 'matched total' parameterization (hold `c.bm_a + c.bm_b =
0.45` with mixing weight `alpha`) for readers requiring equal effect
sizes across all three architectures. Both parameterizations should
be defined; the factorial grid is primary.

## Phase 1: R package changes

### `R/generateData.R`

- Extend `match.arg` to accept `'combined'` as a third level
  alongside `'mvn'` and `'mean_moderation'`.
- Add a `combined` branch after the existing `mean_moderation` block
  (lines 111-128). The branch applies the mean shift using
  `modelparam$c.bm_a`; `buildSigma` uses `modelparam$c.bm_b` for
  the covariance encoding.
- Backward compatibility: when `dgp_architecture` is `'mvn'` or
  `'mean_moderation'`, the existing single `c.bm` entry is used
  unchanged.

### `buildSigma` (inside `generateData.R`, lines 262-272)

- Thread a conditional lookup: use `modelparam$c.bm_b` for the
  BM-BR correlation block when `dgp_architecture == 'combined'`,
  fall back to `modelparam$c.bm` otherwise. No structural change.

### `R/generateSimulatedResults.R`

- Verify that `modelparams` rows with both `c.bm_a` and `c.bm_b`
  columns are passed through to `generateData` without being dropped
  or coerced. No deep change expected; the function passes
  `modelparam` as a list row.

## Phase 2: Simulation driver

New scripts under `analysis/scripts/dgp-combined/`:

- `01-run-combined-factorial.R`: same design matrix as Section 3
  (CO, Hybrid, OL+BDC; `carryover_t1half` in {0, 0.5, 1.0};
  `N = 35` primary, `N = 70` robustness). Outer grid:
  `c.bm_a` x `c.bm_b` in {0, 0.22, 0.45}^2. Save chunks to
  `analysis/data/combined-factorial-*.rds`.
- `02-summarise-combined.R`: tidy summary table matching the
  structure of existing Section 3 tables.

Validation gate before accepting results: confirm the (0.45, 0) and
(0, 0.45) cells reproduce the existing Architecture A and B power
numbers.

## Phase 3: Manuscript additions (paper 01)

### New subsection 2.3 (current 2.3 Mathematical Comparison renumbers to 2.4)

`2.3 Architecture C: Combined Mean Moderation and Differential
Correlation`

- Defines the DGP formally.
- States the boundary cases.
- Motivates separate calibration: the mean-moderation pathway
  (pharmacodynamic scaling by biomarker) and the correlation pathway
  (differential co-regulation) are mechanistically distinct and need
  not scale proportionally.
- Includes the calibration note on matched-total parameterization.

### Section 3 extensions

- Third results block (after current Architecture A and B tables):
  the 3x3 `c.bm_a` x `c.bm_b` power grid at `N = 35`.
- New figure: 3x3 heatmap panel where rows = `c.bm_a`, columns =
  `c.bm_b`, cells = power at `carryover_t1half = 1.0`. Edges of
  the grid reproduce existing A and B heatmaps.
- Extend Section 3.2 to address Architecture C: the two-channel
  attack intensifies monotonically with `c.bm_b`; the A-channel
  signal degrades in the same way as pure Architecture A.
- Robustness check at `N = 70` (new block in Section 3.3).

### New subsection 4.4

`4.4 When Architecture C is appropriate`: both channels are active
when (a) the drug affects mean BR response via a biomarker-scaled
pharmacodynamic mechanism (Architecture A) and (b) on-drug and
off-drug periods have genuinely different BM-BR covariance structure
(Architecture B), for instance when the drug alters signal
transduction pathways that couple biomarker expression to treatment
response.

### Section 5 update

Update 5.1 ('Choosing the appropriate architecture') to present
Architecture C as a third case with guidance: if both mechanisms are
plausible but their relative weights are unknown, run sensitivity
analyses spanning the `(c.bm_a, c.bm_b)` grid rather than committing
to one architecture.

## Ordering

1. Phase 1 (package changes) -- required before simulation can run.
2. Validation gate: (0.45, 0) and (0, 0.45) cells reproduce existing
   results.
3. Phase 2 (simulation driver) -- can run overnight.
4. Phase 3 (manuscript) -- draft Section 2.3 and 4.4 before results
   arrive; fill tables and figure once simulation completes.
