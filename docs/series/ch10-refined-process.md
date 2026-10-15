# The refined process: specification, the adjustments evaluated,
and the roadmap

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Authoritative definition; difference table from
published; operating characteristics under the refined process.
**Prerequisites:** ch04, ch09

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

*Local symbols used in this chapter only. For project-wide notation
see `analysis/report/NOTATION.md`.*

Named processes follow NOTATION.md "Named data-generating processes":
- **the published process:** Hendrickson's `generateData()` at commit
  58b32a9.
- **the refined process:** the recommended DGP (Paper 01).
- **configuration A-F:** the covariance and coupling constructs
  compared in `docs/36`.

---

## Summary

<!-- evidence-labeled results go here -->

The refined process is the recommended data-generating process for
the compendium's power simulations. It is defined in Paper 01, Table
`tab:spec`. Its five defining features are:

1. Graded coupling ($g_t = 0$ before first exposure; ch08).
2. Separable AR(1) with $\rho = 0.7$ and $c_1 = 0.2$ (configuration C;
   ch09).
3. Anchored carryover on the drug-response mean (ch06).
4. A stop on matrices that are not positive definite (ch09).
5. Randomized path allocation.

It differs from the published process in four respects: the coupling
rule (defect, not preference), the within-factor temporal correlation
form (modeling choice), the cross-factor temporal correlation form
(consequence of AR(1)), and the carryover recursion (repair; ch06).

Eight roadmap steps (`docs/44`) take the published process to the
refined process; steps 1 to 6 are implemented (steps 5 and 6 in the
quick-sim drivers; packages have steps 2, 3 and 4); steps 7 and 8
are proposals.

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3.3 for this chapter:
- Sources: 36 Sections 6.1 to 6.6, 6.8; 44 Sections 1, 3, 5, 6;
  32 Section 2 (repair, choice, consequence); NOTATION "Named
  data-generating processes"; paper 01 tab:spec; 50 Section 1;
  51 Section 20; 23 Section 2 (history only).

- New to write: the difference table; the refined-process operating
  characteristics from job 14 and job 15; the unreported common-grid
  contrast behind 36-fig4-contrast-reps100.png.
-->

### Named processes

<!-- From NOTATION.md "Named data-generating processes" -->

| Name | Meaning | Defined in |
|---|---|---|
| the published process | Hendrickson's `generateData()` at 58b32a9: CS ($\rho_{cs} = 0.8$), constant cross-component, step coupling, compounding recursion, zero restart | Paper 01 |
| the refined process | graded coupling (0 before first exposure), separable AR(1) with $\rho = 0.7$ and $c_1 = 0.2$, anchored carryover, stop on non-PD matrices, randomized allocation | Paper 01, Table `tab:spec` |
| configuration A | published compound-symmetry construct with whatever coupling the context states | `docs/36` |
| configuration C | separable AR(1), graded coupling | `docs/36` |
| configuration F | AR(1) with `covar` cross block, graded coupling | `docs/36` |
| `orig` | the published process's covariance construct as implemented | `docs/32`, `docs/36` |
| `covar` | the package's covariance construct: AR(1) within component with decaying cross-component block (configuration F at package parameters) | `docs/32` |

### The four differences from the published process

<!-- From docs/32 Section 2 -->

| Component | Published process (`orig`) | Refined process |
|---|---|---|
| Within-factor correlation | CS, $\rho_{cs} = 0.8$ at every lag | AR(1), $\rho^{\lvert w_t - w_s \rvert}$ with $\rho = 0.7$ |
| Cross-factor correlation ($t \ne s$) | constant $c_\times = 0.1$ | separable: $c_1 \rho^{\lvert w_t - w_s \rvert}$ |
| Coupling off drug | step: $c_{bm}$ wherever $\mu_{BR,t} \ne 0$ | graded: $c_{bm} e^{-\lambda t_{sd,t}}$ off drug, 0 before first exposure |
| Carryover on $\mu_{BR}$ | compounding recursion | anchored decay from last on-drug value |

The coupling change is a defect repair (the step rule removes the
interaction signal at any positive half-life; `docs/32` Section 3,
`docs/37`). The within-factor form is a modeling choice (`docs/32`
Section 4). The cross-factor form follows from AR(1) when separability
is adopted (`docs/32` Section 5). The carryover recursion is a repair
(`docs/53`; ch06).

<!-- TODO: full difference table per Paper 01 tab:spec: list every
parameter and whether it is the same as or different from the
published process, with the rationale column from docs/32 -->

### The eight roadmap steps

<!-- From docs/44 Sections 1 and 3 -->

Eight steps take the published process to a fully specified generative
model (`docs/44`):

| Step | Change | Status in refined process |
|---|---|---|
| 1 | Fail loudly instead of repairing | yes: cells with non-PD matrix are stopped |
| 2 | Graded coupling | yes: $g_t = 2^{-t_{sd}/t_{1/2}}$ off drug |
| 3 | Anchored carryover | yes: $\mu_{BR,t} = \mu_{BR,\text{last}} \cdot 2^{-t_{sd}/t_{1/2}}$ |
| 4 | AR(1) within-component | yes: $\rho = 0.7$ per week |
| 5 | Separable cross-component and factored generator | yes (quick-sim drivers); not yet in package |
| 6 | $\rho$ and $c_{bm}$ inside the ceiling | yes: $\rho = 0.7$, $c_{bm} = 0.45$ is feasible |
| 7 | Restart from the residual on re-exposure | no (zero restart) |
| 8 | Pharmacodynamic sensitivity model | no (proposal; ch14) |

Steps 1 to 6 define the refined process. Steps 7 and 8 are proposals
awaiting implementation. The package has steps 2, 3 and 4 but not 1
(silent Cholesky fallback), 5 (uses decaying cross-factor form) or 6
in full.

### Power under the adjustments

<!-- From docs/36 Section 6.8 -->

Nominal Hybrid power at $c_{bm} = 0.25$, `corCAR1` analysis:

| Configuration | $t_{1/2} = 0$ | 0.1 | 0.2 | 0.5 | 1 |
|---|---|---|---|---|---|
| A Published (step, CS) | 0.59 | 0.08 | 0.08 | 0.09 | 0.10 |
| F AR(1) + graded, `covar` | 0.43 | 0.41 | 0.44 | 0.36 | 0.23 |
| C AR(1) + graded, sep. ($c_{bm} = 0.25$) | 0.51 | 0.48 | 0.52 | 0.41 | 0.28 |
| C AR(1) + graded, sep. ($c_{bm} = 0.30$) | 0.69 | 0.70 | 0.72 | 0.55 | 0.36 |
| D Mean moderation ($c_{bm} = 0.30$) | 0.78 | 0.75 | 0.74 | 0.77 | 0.73 |

(Verified, `power-summary-reps500-car1-anchored.csv`; `docs/36`
Table 9.)

Five findings (`docs/36` Section 6.8):

1. AR(1) alone does not help: configuration B collapses as A does
   (step rule retained).
2. Graded coupling removes the collapse: at $t_{1/2} = 0.1$ and 0.2,
   power under C equals its no-carryover value.
3. Mean moderation is flat: no power loss at any half-life (the
   plain profile stops when the drug stops, which the binary
   indicator codes correctly).
4. The response structure matters for CO: D gives 0.85, E gives 0.44
   at $c_{bm} = 0.3$ (compound symmetry versus separable AR(1) at
   2.5-week spacing).
5. Separability changes feasibility, not power: F and C have similar
   shaped curves; C's case rests on structure (valid for every
   schedule).

### The ceiling and feasibility

At the refined process's $\rho = 0.7$, the Hybrid coupling ceiling
for configuration C is 0.480 (path A; 0.491 on path C; `docs/35`,
Section 6). The compendium's standard effect size $c_{bm} = 0.45$ is
feasible with a margin of at least 0.02 in every design and half-life.

### Analysis implications

From `docs/44` Section 5:

| DGP step | Analysis implication |
|---|---|
| Graded coupling, anchored carryover | exposure-weighted $D_{bc}$ matches the coupling decay |
| AR(1) | working covariance must allow serial correlation: `corCAR1` with CR2 |
| Feasible effect sizes | larger effects become testable; compare at matched feasible strengths |

The selected analysis is F1: random intercept with `corCAR1` residuals
and CR2 inference, baseline row, linear time, $D_{bc}$ at a fixed
$t_{1/2}^a = 0.2$ week (or F1A0B with AIC-chosen half-life), and the
expectancy term `De` with its biomarker interaction (`docs/49`;
`docs/46` Section 7).

### New-to-write items

<!-- TODO: the difference table per Paper 01 tab:spec: every
parameter, orig vs refined, rationale column -->

<!-- TODO: the refined-process operating characteristics from job 14
and job 15 (the stress-test reruns under the refined process with
anchored carryover); summarize size, maximin power and key survivor
analyses; see docs/47 "Rerun under the refined process" section -->

<!-- TODO: the unreported common-grid contrast behind
36-fig4-contrast-reps100.png: both orig and the fixes run on a common
grid ($c_{bm} \in \{0.25, 0.30\}$, $t_{1/2} \in \{0, 0.1, 0.2, 0.5,
1\}$, $N = 70$) for a side-by-side power comparison not in docs/36's
main figures; report it here -->

---

## Derivations

| Statement | Source | Script |
|---|---|---|
| Ceiling at $\rho = 0.7$, config. C, Hybrid: 0.480 | `docs/35` Section 6 | `check-landmarks-hybrid-ceilings.R` |
| AR(1) analysis needed to avoid inflated Type I error | `docs/36` Section 6.8 | `power-summary-reps500-car1-anchored.csv` |

---

## Numerical evidence

See power table above (docs/36 Section 6.8) and ceiling table in
ch09.

<!-- TODO: refined-process stress-test results (jobs 14 and 15) -->

---

## Relation to the manuscripts

- **Paper 01:** defines the refined process; Table `tab:spec`; reports
  power contrast against the published process.
- **`docs/36`:** the three adjustments and their power consequences.
- **`docs/44`:** the eight roadmap steps and package status.
- **`docs/32`:** the four differences from the published process,
  classified as repair / choice / consequence.

---

## Limitations and open items

- The package does not yet implement the separable cross-component
  form (configuration C) or the factored generator (step 5).
- The silent Cholesky fallback in the package is still present (step 1
  partial); inspected 2026-10-07, not yet removed.
- Stored covariance-study and stress-test results used the published
  recursive carryover; the anchored reruns (jobs 14 and 15) are the
  source of the refined-process operating characteristics.

---

## Reproduce

```bash
# Adjustment power simulation (docs/36 Section 6.8)
Rscript analysis/scripts/quick-sim/hendrickson-problems/\
  04-power-simulation.R
# Common-grid contrast
Rscript analysis/scripts/quick-sim/hendrickson-problems/\
  12-common-grid-contrast.R
```

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | All | Phase 2 content move from docs/36 sec 6.1-6.8, docs/44 sec 1,3,5,6, docs/32 sec 2, NOTATION.md, docs/51 sec 20 | Whitepaper refactor |
