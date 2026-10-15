# Note: patching 58b32a9

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** The eight-line graded patch, the anchored fix, and their
joint effect on power.  
**Prerequisites:** ch26, ch27

---

## Contents

1. [Notation](#notation)
2. [Summary](#summary)
3. [Construct box](#construct-box)
4. [Setting and argument](#setting-and-argument)
5. [Derivations](#derivations)
6. [Numerical evidence](#numerical-evidence)
7. [Relation to the manuscripts](#relation-to-the-manuscripts)
8. [Limitations and open items](#limitations-and-open-items)
9. [Reproduce](#reproduce)
10. [References](#references)
11. [Change log](#change-log)

---

## Notation

*Local symbols used in this chapter only. For project-wide notation see
`analysis/report/NOTATION.md`.*

<!-- notation entries -->

---

## Summary

<!-- evidence-labeled results go here -->

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 37 Appendix A; 53 Section 9; 55 Sections 9 to 10 (Table
  8, runs A to C); 36 patch material (not found in docs/36 -- see
  note below); the patch README at analysis/patches/orig-to-refined/.

- New to write: the single power table combining runs A to C.

- Note: docs/36 Section 6.8 covers Type I error under adjustments,
  not patch content. No patch material was found in docs/36.
-->

### Two defects in commit 58b32a9

<!-- From docs/55 Section 10 -->

Commit `58b32a9` (the published Hendrickson code, vendored at
`implementations/original/`) contains two independent defects that
together collapse power at any nonzero carryover half-life.

**Defect 1: step coupling** (ch26). The sigma matrix construction
uses a step rule: $\mathrm{Cor}(BM, BR_t) = c_{bm}$ whenever the
mean drug response is positive, and 0 otherwise. Under the
pharmacokinetic model, the drug response decays continuously after
discontinuation, so the mean is positive throughout the trial. This
sets $\mathrm{Cor}(BM, BR_t) = c_{bm}$ at off-drug visits for any
$t_{1/2} > 0$, removing the between-path covariance that identifies
the interaction. Power collapses to approximately 0.09 at
$t_{1/2} \ge 0.1$ week (inspected; `docs/55` Table 8, column
"Published").

**Defect 2: carryover recursion** (ch27). The DGP adds the carryover
contribution to the current mean at each off-drug step, rather than
computing it directly from the time since last discontinuation. This
compounds the decay across off-drug visits and produces carryover
effects far larger than the specified $t_{1/2}$ implies. At
$t_{1/2} = 1$ week the compounding increases the effective half-life
by a factor of approximately 3 to 4 (`docs/53` Section 9; inspected).

The two defects interact: step coupling drives most of the power
collapse at short half-lives; compounding carryover inflates the
off-drug response at longer half-lives, masking the interaction by
increasing off-drug noise. The defects are independent in the code
and can be patched separately (`docs/55` Section 10; inspected).

### The graded coupling patch

<!-- From docs/37 Appendix A -->

The 8-line graded coupling patch replaces the step rule with an
exponential decay of the biomarker-response correlation during
off-drug periods:

**Published code (step rule):**

```r
# Coupling active whenever mean drug response > 0
sigma_bm_br[t] <- ifelse(mean_br[t] > 0, c.bm, 0) * sigma_bm * sigma_br
```

**Replacement (graded coupling):**

```r
# Coupling decays with time since discontinuation
coupling_t <- ifelse(
  on_drug[t],
  c.bm,
  ifelse(tsd[t] > 0, c.bm * exp(-lambda_cor * tsd[t]), 0)
)
sigma_bm_br[t] <- coupling_t * sigma_bm * sigma_br
```

The graded rule sets coupling to $c_{bm} e^{-\lambda t_{sd}}$ during
off-drug periods after first exposure and 0 before first exposure.
This matches the pharmacokinetic interpretation of the biomarker-
response association as mediated by drug plasma concentration.

**Power at check run** (inspected from `docs/37` Appendix A Table):

| $t_{1/2}$ (weeks) | Published power | Graded power |
|---|---|---|
| 0 | 0.668 | 0.668 |
| 0.1 | 0.126 | 0.638 |
| 0.2 | 0.092 | 0.656 |
| 1.0 | 0.094 | 0.462 |

The graded patch restores power to near the no-carryover baseline at
short half-lives. At $t_{1/2} = 1$ week the power loss is real
(pharmacologically, carryover of that length contaminates the off-drug
observations) and is not an artifact.

**What the graded patch does not fix:**

- The positive-definiteness ceiling (still 0.256 under compound
  symmetry; fixed separately by separable AR(1)).
- Visits before first exposure incorrectly coded as exposed in some
  paths (`docs/37`, Appendix A; inspected).
- Carryover recursion (compounding carryover; patched separately).

### The anchored carryover patch

<!-- From docs/53 Section 9 -->

The anchored carryover fix changes one line in the DGP's off-drug
mean computation. The published code carries the previous visit's mean
forward, compounding the decay:

```r
# Published (compounding):
mu_off[t] <- mu_off[t-1] * exp(-lambda_cor * dt)
```

The patch anchors the decay to the last on-drug mean:

```r
# Anchored:
mu_off[t] <- mu_last_on * exp(-lambda_cor * tsd[t])
```

where `mu_last_on` is the last on-drug visit's mean response and
`tsd[t]` is the cumulative time since discontinuation. This is the
pharmacokinetically correct form: the off-drug response is a fraction
of the last on-drug response, not a fraction of a fraction.

The upstream companion change (required for the anchored DGP to
pass the analysis model's carryover term) is to replace the
recursive `Dbc` computation in the analysis data-preparation function
with the direct formula `exp(-lambda * tsd)` (verified;
`docs/53` Section 9).

**Re-exposure.** After a second on-drug period, the anchored formula
takes `mu_last_on` from the second on-drug period's final visit.
This is the correct pharmacokinetic interpretation: each on-drug
exposure resets the reference point. The published code also resets
the reference point (because it multiplies from the current mean),
but the path to that reset is through compounded decay and may
accumulate error at long half-lives.

### Runs A, B, C under both defects and individual patches

<!-- From docs/55 Sections 9 to 10 (Table 8) -->

Three simulation runs characterize the patches:

- **Run A:** the published Hendrickson allocation (three paths at 46%,
  31%, and 23%; inspected from `docs/55`).
- **Run B:** equal allocation (25% per path).
- **Run C:** the allocation from the published paper's Figure 4 (1%,
  66%, 33%, traced from the figure; inspected).

Power by run and coupling rule ($N = 70$, Hybrid, $n_{rep} = 1000$):

| $t_{1/2}$ | Pub A | Graded A | Pub B | Graded B | Pub C | Graded C |
|---|---|---|---|---|---|---|
| 0 | 0.668 | 0.668 | 0.628 | 0.628 | 0.670 | 0.670 |
| 0.1 | 0.126 | 0.638 | 0.078 | 0.658 | 0.100 | 0.656 |
| 0.2 | 0.092 | 0.656 | 0.082 | 0.670 | -- | -- |
| 1.0 | 0.094 | 0.462 | 0.088 | 0.486 | 0.076 | 0.448 |

(Inspected from `docs/55` Table 8; "--" indicates cell not run in
Table 8.)

For the 2024 upstream package (vendored at
`implementations/original-extended/`), graded coupling at $t_{1/2} = 1$
week gives 0.368 (run A), compared to 0.462 in run A above. The 2024
version incorporates a partial coupling correction that recovers
some power but does not reach the separable AR(1) baseline
(`docs/55`, Section 10; inspected).

### Applying the patches

<!-- From analysis/patches/orig-to-refined/README.md lines 1-199 -->

The `analysis/patches/orig-to-refined/` directory contains two patch
files and a README:

- **`graded-coupling.patch`:** applies the 8-line graded coupling
  change to `implementations/original/` and
  `implementations/original-extended/`.
- **`anchored-carryover.patch`:** applies the 1-line anchored
  carryover change to the same two collections.

Apply with:

```bash
cd implementations/original
patch -p1 < ../../analysis/patches/orig-to-refined/graded-coupling.patch
patch -p1 < ../../analysis/patches/orig-to-refined/anchored-carryover.patch
```

Both patches pass the parity tests at zero carryover (`make parity`;
verified from the README). Non-zero carryover parity is expected to
diverge from the published code, because the published code has the
defects the patches correct.

The `options` patch table in the README documents default changes:
`carryover_t1half` default changed from `NA` (no carryover) to 0.5,
and `dgp_architecture` default changed from `'mvn'` to `'mvn'` with
anchored carryover active. The `analysis` patch table documents the
companion analysis-side change to the `Dbc` computation.

<!-- TODO: single combined power table summarizing runs A, B, C
across t_half values for both published and graded coupling, plus
anchored carryover vs compounding carryover, in one 6-column table
for easy comparison. This requires re-running the simulations or
assembling from existing output files in analysis/data/. -->

---

## Derivations

<!-- one row per result: statement | CM label | verification script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the numbers
agree? -->

---

## Limitations and open items

<!-- link each to the open-items register at
docs/series/registry/open-items.csv -->

---

## Reproduce

<!-- one command per output -->

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | Setting and argument | Initial content move from docs/37 App A, docs/53 §9, docs/55 §§9-10, patch README | Phase 2 content moves |
