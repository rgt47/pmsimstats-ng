---
title: "Simulation grid plan for manuscript 02 (carryover sensitivity)"
author: "pmsimstats team"
date: "2026-04-15"
---

# Simulation grid plan: manuscript 02

This document defines the parameter grid for the manuscript-02
sensitivity simulation and records the rationale for each
parameter choice. Sections labelled **Rationale (for manuscript)**
are drafted for lift into the Methods section.

## Design principle

A full-factorial over all parameters of interest is
combinatorially infeasible (adding one three-level axis multiplies
the grid by three). We therefore adopt a common two-tier
sensitivity-analysis design:

- **Tier 1: principal factorial.** A compact factorial over the
  axes that together define the headline scientific question.
  Every combination is evaluated at production replicate counts.
- **Tier 2: marginal sensitivity blocks.** One-factor-at-a-time
  (OFAT) sweeps anchored at a *reference configuration* drawn
  from Tier 1. Each block answers a single sensitivity question
  without multiplying the principal grid.

This pattern is standard in simulation methodology studies where
the main claim is a ranking or contrast rather than a precise
power estimate. It is used in Schork (2022) for N-of-1 serial
correlation sensitivity and in Sturdevant & Lumley (2021) for
carryover-analysis-model sensitivity; the current manuscript
sits at the intersection of those two literatures.

## Tier 1: principal factorial

**Scientific question.** Does the ranking of analysis-model
carryover specifications (A1 binary, A2 exposure-weighted, A3
lagged-nuisance) depend on (a) the true DGP decay form, (b) the
DGP architecture, and (c) the carryover half-life, across the
trial-design and sample-size settings used in precision-medicine
N-of-1 work?

Tier 1 crosses every axis whose variation is load-bearing for
that question. Cell count: $5 \times 3 \times 2 \times 3 \times
3 \times 2 \times 3 = 540$.

| Axis | Values | Levels |
|---|---|---|
| DGP carryover form | linear, exponential, Weibull($k \in \{0.7, 1.0, 1.5\}$) | 5 |
| Analysis specification | A1 binary, A2 exposure-weighted, A3 lagged | 3 |
| DGP architecture | mean moderation, MVN differential correlation | 2 |
| Carryover half-life $t_{1/2}$ (weeks) | 0, 0.5, 1.0 | 3 |
| Trial design | CO, Hybrid, OL+BDC | 3 |
| Sample size $N$ | 35, 70 | 2 |
| Biomarker moderation $c_{bm}$ | 0, 0.30, 0.45 | 3 |

### Rationale (for manuscript)

**DGP carryover form.** Three forms are evaluated: linear,
exponential, and Weibull with three shape parameters. Exponential
is the first-order-elimination default in pharmacokinetics and is
the form used throughout the Hendrickson (2020) framework; it is
the natural reference cell. Linear is biologically implausible
but represents the minimum-information assumption a designer can
make from "carryover has cleared by week $X$"; quantifying how
much the linear assumption costs answers a practical question
trialists face when only a bound is available. Weibull
generalises exponential by adding a shape parameter $k$ that
captures accelerated washout ($k > 1$; drugs with active-transport
or saturable clearance) and slow elimination ($k < 1$; drugs with
deep-compartment sequestration). Shape values $k \in \{0.7, 1.5\}$
bracket a realistic deviation from exponential while remaining
identifiable in a half-life parameterisation; $k = 1$ recovers
exponential and serves as an internal consistency check.

**Analysis specification.** Three specifications span the practical
analyst choices:

- A1 (binary-treatment) tests the damage from ignoring carryover
  entirely and provides a reference for what the other
  specifications are buying.
- A2 (exposure-weighted) is the current pmsimstats-ng default and
  the specification used by Hendrickson (2020). It assumes the
  analyst posits a specific decay form and half-life.
- A3 (lagged-nuisance) follows the Jones & Kenward (2014)
  crossover convention of including a lagged-treatment indicator
  as a nuisance covariate rather than parameterising the decay.
  It requires no assumption about functional form but is
  typically less efficient when that form is correctly specified.

**DGP architecture.** Both architectures from the companion
manuscript are retained. Architecture A encodes the
biomarker-treatment interaction in the mean structure via a
standardised-biomarker-scaled shift on on-drug BR; Architecture B
encodes it in treatment-state-dependent correlation within a
joint multivariate normal. The companion manuscript establishes
that the two produce qualitatively different power behaviour
under exponential carryover; the question here is whether that
qualitative gap is robust to decay-form variation.

**Carryover half-life.** Values $\{0, 0.5, 1.0\}$ weeks span no
carryover to a half-life comparable with the dosing interval in
the motivating prazosin-PTSD application (Raskind et al. 2013).
Longer half-lives produce exposure-decay profiles that are
nearly indistinguishable from sustained on-drug exposure for the
analysis-model predictors evaluated here; shorter half-lives
collapse onto the no-carryover case.

**Trial design.** Three designs are retained: traditional
crossover (CO), Hybrid N-of-1, and open-label with blinded
discontinuation (OL+BDC). Open-label (OL) without any off-drug
phase is omitted because it lacks the treatment contrast
necessary to identify A1 and A3; OL appears only as a reference
point in the companion manuscript.

**Sample size.** $N \in \{35, 70\}$ spans small single-site
pilots to pooled multi-site samples typical of the aggregated
N-of-1 literature. Larger samples ($N > 100$) push most cells
above 0.9 power and reduce the comparative sensitivity of the
ranking we are evaluating.

**Biomarker moderation $c_{bm}$.** Three levels serve three
distinct roles. $c_{bm} = 0$ is the null cell for type I error
verification. $c_{bm} = 0.30$ is a moderate effect typical of
biomarkers with partial predictive value (polygenic risk scores,
continuous severity indices). $c_{bm} = 0.45$ matches the
companion manuscript's reference case, chosen in Hendrickson
(2020) as the largest value that passes positive-definiteness
under the Architecture B covariance structure.

## Tier 2: marginal sensitivity blocks

Each block varies one additional axis against the **reference
configuration** below, drawn from the most demanding but still
representative cell in Tier 1:

```
reference = {
  carryover_form = 'exponential',
  weibull_shape  = 1.0,
  dgp_arch       = 'mvn',           # the sensitive architecture
  t1half         = 1.0,              # the hardest carryover setting
  design         = 'Hybrid',         # the motivating N-of-1 design
  N              = 70,
  c_bm           = 0.45
}
```

All Tier-2 blocks cross the varied axis with the three analysis
specifications (A1 / A2 / A3), so every block yields a 3-column
comparison of how the spec ranking changes.

### Block S1: autocorrelation sensitivity

**Axis.** AR(1) within-factor autocorrelation $\rho \in \{0.5,
0.7, 0.8, 0.9\}$.

**Cell count.** $4 \times 3 = 12$.

**Rationale.** Hendrickson (2020) uses $\rho = 0.8$; the revised
analysis in `docs/02` uses $\rho = 0.7$ after positive-definiteness
issues at $\rho = 0.8$. Within-factor autocorrelation is the
single largest non-varied determinant of effective sample size in
longitudinal N-of-1 analysis (Schork 2022), and its value is
typically unknown at the design stage. Block S1 quantifies how
the analysis-specification ranking moves as $\rho$ varies across
plausible operating points. We predict that A2 (which uses the
continuous-time corCAR1 structure most heavily) is most sensitive
to $\rho$, and that A3 (which partially sidesteps the
autocorrelation via the lagged-treatment indicator) is least
sensitive.

### Block S2: analyst-truth half-life mismatch

**Axis.** Cartesian product of DGP half-life $t_{1/2}^{\text{true}}
\in \{0.5, 1.0\}$ and analyst-assumed half-life in A2 and DGP
decay form $t_{1/2}^{\text{analyst}} \in \{0.25, 0.5, 1.0, 2.0\}$.

**Cell count.** $2 \times 4 \times 3 = 24$.

**Rationale.** The Tier 1 grid holds the analyst's assumed
half-life equal to the DGP's, corresponding to perfect prospective
knowledge of the drug's pharmacokinetic half-life. In practice
the analyst works from a point estimate that may be off by a
factor of two or more. Block S2 implements the classical Jones &
Kenward (2014) "model misspecification" comparison: how does
power (and convergence) degrade when the assumed half-life
differs from truth, and does A3 (which does not require the
analyst to specify a half-life) become preferred as the mismatch
grows? This block directly answers a question trialists face when
pharmacokinetic information is limited.

### Block S3: dropout sensitivity

**Axis.** Dropout rate $d \in \{0, 0.1, 0.2, 0.3\}$ under two
mechanisms (MCAR, MAR biased by baseline severity), using the
`R/censordata.R` facility. Crossed with analysis specification.

**Cell count.** $4 \times 2 \times 3 = 24$.

**Rationale.** Biomarker-moderated trials in clinical populations
typically experience 10-30% dropout. The A2 and A3 specifications
rely on different portions of the data for identification: A2 on
smoothly exposure-weighted observations across the entire
trajectory, A3 on the presence of a lagged-on timepoint. MAR
dropout biased by baseline severity will differentially cost
each specification, and quantifying this cost is important for
trialists deciding between them at the design stage.

### Block S4: effect-size curve

**Axis.** $c_{bm} \in \{0, 0.10, 0.20, 0.30, 0.45, 0.60\}$
crossed with analysis specification.

**Cell count.** $6 \times 3 = 18$.

**Rationale.** Tier 1's three levels establish type-I control at
$c_{bm} = 0$ and two points on the power curve. Block S4 adds
three levels to characterise the full power curve's slope, which
is what trialists actually need for sample-size planning. The
$c_{bm} = 0.60$ value is only feasible under Architecture A (the
MVN covariance fails positive-definiteness above approximately
$c_{bm} = 0.45$); we report $c_{bm} \in \{0.60\}$ only for
Architecture A with a note.

### Block S5: rho-by-carryover interaction (exploratory)

**Axis.** $\rho \in \{0.5, 0.8\}$ crossed with $t_{1/2} \in
\{0, 0.5, 1.0\}$ crossed with analysis specification.

**Cell count.** $2 \times 3 \times 3 = 18$.

**Rationale.** Block S1 varies $\rho$ at a single carryover
level; Block S5 lets us check whether the $\rho$-sensitivity of
the spec ranking interacts with carryover strength. This is
flagged exploratory because we have no strong prior and the
conclusion may be "no interaction," in which case the block
becomes a methodological reassurance rather than a finding.

## Aggregate computational cost

| Block | Cells | Replicates | Cell-fits | Wall-clock (8-core) |
|---|---|---|---|---|
| Tier 1 | 540 | 500 | 810,000 | ~85 min |
| S1 | 12 | 500 | 18,000 | ~2 min |
| S2 | 24 | 500 | 36,000 | ~4 min |
| S3 | 24 | 500 | 36,000 | ~4 min |
| S4 | 18 | 500 | 27,000 | ~3 min |
| S5 | 18 | 500 | 27,000 | ~3 min |
| **Total** | **636** | | **954,000** | **~100 min** |

Marginal blocks add approximately 15 minutes to the Tier 1
production run. The total fits in one overnight execution.

## Implementation status

- **Tier 1:** implemented in `01-run-factorial.R`; production run
  currently in flight.
- **S1 (rho):** requires plumbing `rho` through `model_param`
  in the tidyverse implementation. Currently fixed at 0.7.
- **S2 (half-life mismatch):** requires splitting
  `carryover_t1half` into DGP and analysis values in
  `prepare_long_data()`. Single-line change.
- **S3 (dropout):** requires wiring `censor_data()` into the
  per-replicate loop of `simulate_cell()`.
- **S4 (effect-size curve):** requires no code change; extend
  the `c_bm` vector in a separate driver script to avoid bloating
  Tier 1.
- **S5 (rho-by-carryover):** depends on S1 plumbing.

Script `04-run-sensitivity-blocks.R` will orchestrate the five
blocks after Tier 1 is archived.

## What is *not* in the plan

The following axes are deliberately excluded from the manuscript-02
grid:

- **Washout period length / timepoint density.** Design-level
  mitigation is orthogonal to the analysis-specification question
  and is covered in `docs/02` §6.6. A dedicated paper on design-
  level mitigation is tracked in the research portfolio.
- **Gompertz response trajectory parameters (max, rate, disp).**
  Already explored in `analysis/figure5`; secondary to the
  central question here.
- **Random slopes vs intercept-only random effects.** Orthogonal
  to carryover specification and would double the grid; a brief
  note in the Discussion will point to the Hendrickson (2020)
  sensitivity analysis that addressed this.
- **Multiple-testing correction.** Out of scope; each simulation
  tests a single prespecified biomarker.
- **Alpha level.** Type-I error is evaluated at the standard
  $\alpha = 0.05$; a supplementary check at $\alpha = 0.10$ is
  straightforward post hoc if a reviewer requests it.

---

*Rendered on 2026-04-15.*
*Source: ~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/02-carryover-sensitivity/simulation-grid-plan.md*
