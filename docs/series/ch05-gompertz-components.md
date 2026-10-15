# Response components and the modified Gompertz trajectory

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** BR, PB, TV, BL_i; their clocks; modified Gompertz
G(t_od; m, d, r); alternative families.
**Prerequisites:** ch02

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

| Symbol | Meaning |
|---|---|
| $m$, $d$, $r$ | Gompertz ceiling, displacement, rate (code: `maxr`, `disp`, `rate`) |
| $t_{od}$ | cumulative time on drug (clock for BR) |
| $t_{pb}$ | expectancy-weighted time (clock for PB) |
| $t$ | calendar week (clock for TV) |
| $t^\ast$ | inflection time $\ln(d)/r$ |
| $\eta_t$ | design expectancy weight at visit $t$ |
| $G(t; m, d, r)$ | rescaled modified Gompertz curve |

---

## Summary

<!-- evidence-labeled results go here -->

The symptom score decomposes as
$Y_{it} = \mathrm{BL}_i - [BR_{it} + PB_{it} + TV_{it}] +
\varepsilon_{it}$. Each latent component follows a modified Gompertz
trajectory in its own time argument. The rescaled curve passes through
zero and converges to $m$; the displacement $d$ sets the lag length
(not the steepness) of the rise; the inflection is at $t^\ast =
\ln(d)/r$.

At the shipped drug-response parameters ($m = 10.986$, $d = 5$,
$r = 0.42$ per week) the inflection is at 3.83 weeks, half the
ceiling is reached at 4.73 weeks, and 90% at 9.21 weeks (verified,
`check-landmarks-hybrid-ceilings.R`). Participant heterogeneity is
level-only: $(m, d, r)$ are population constants; individuals differ
through MVN deviations around the mean curve.

Paper 07 finds that power differences across four response-curve
families (Gompertz, logistic, hyperbolic-tangent,
piecewise-linear) are small under both architectures (no carryover,
$N = 70$, OL+BDC design). The covariance channel is shape-insensitive
by construction; the mean channel is mildly shape-sensitive.

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3.3 for this chapter:
- Sources: 24 Components 1 to 3; 25 Sections 2 to 3, 5 to 8
  (literature, alternatives, ranked impact); 40 Sections 3 to 4;
  CM "The response trajectory" (4 subsections); paper 07;
  18 via Appendix I.

- New to write: a family-overlay figure; 25 Section 6 and
  Recommendation 4 rewritten to the corrected account of Section 3;
  the link from 25 Recommendation 5 to P04's v_i*Psi_it.
-->

### Three components and their clocks

<!-- From docs/24 Components 1 to 3 -->

The BR-PB-TV decomposition writes the observed change in outcome as a
sum of three biologically interpretable components plus measurement
noise:

$$
Y_{it} = \mathrm{BL}_i - [\,BR_{it} + PB_{it} + TV_{it}\,]
            + \varepsilon_{it}.
$$

Here $\mathrm{BL}_i$ is participant $i$'s baseline symptom score,
$BR$, $PB$, and $TV$ are the three response components (each a
non-negative reduction in symptoms), and $\varepsilon$ is
observation-level noise.

**Component 1: Biological Response ($BR$).** BR is the physiological
response to the drug itself. It is what happens because of the
medication's chemical action, independent of belief, suggestion, or
natural history. BR rises with cumulative time on drug $t_{od}$,
reaches a ceiling $m_{BR}$, and collapses to zero (through
carryover) when the drug is stopped. Its clock resets to zero at each
new on-drug period (zero restart; the turnover model of ch06 and ch14
would restart from the residual).

The prazosin pharmacokinetic profile is short (plasma half-life
approximately 2.3 hours; Hobbs et al. 1978). The carryover that
matters at weekly visits is pharmacodynamic: the response half-life,
not the plasma half-life.

**Component 2: Placebo-Belief Response ($PB$).** PB is the change in
symptom score arising from the patient's belief that they are
receiving treatment. It is scaled at each visit by the design
expectancy weight $\eta_t$ (1.0 open-label, approximately 0.5
blinded discontinuation, less at open-label placebo), which acts on
both the mean and the SD of $PB$. Its clock is $t_{pb}$,
expectancy-weighted time. The reset of $t_{pb}$ at visits with
$\eta_t = 0$ means blinded-phase PB builds from the belief that
started accumulating in the most recent non-zero phase
(inspected: `R/buildtrialdesign.R`).

**Component 3: Natural Course ($TV$).** TV is the improvement (or
worsening) that would have occurred regardless of treatment.
It runs in calendar time $t$ and is not modulated by drug exposure
or expectancy.

### The modified Gompertz curve

<!-- From CM "The response trajectory", lines 572-878 -->

Every latent factor ($BR$, $PB$, $TV$) is generated as a trajectory
that rises from zero toward an asymptote. The default is the modified
Gompertz curve. The package computes (`modgompertz()` in
`R/utilities.R`, inspected):

$$
G(t; m, d, r) = m\cdot
  \frac{\exp(-d\cdot\exp(-r\cdot t)) -
        \exp(-d)}{1 - \exp(-d)},
  \qquad G(0) = 0, \quad G(\infty) = m.
$$

The three parameters:

| Parameter | Symbol | Meaning |
|---|---|---|
| Ceiling | $m$ | Asymptotic maximum response |
| Rate | $r$ | Onset speed; multiplying $r$ by $k$ divides every landmark by $k$ |
| Displacement | $d$ | Length of initial lag; inflection at $t^\ast = \ln(d)/r$; large $d$ delays the inflection but makes the steepest rise slightly less steep, not steeper (derived; corrected 2026-10-08) |

The kernel $y(t) = m e^{-d e^{-r t}}$ does not satisfy $y(0) = 0$;
the rescaling subtracts $m e^{-d}$ and restores the asymptote. As
$d \to 0$ the curve tends to $m(1 - e^{-rt})$ (the exponential
approach). For $d < 1$ no inflection occurs on $t \ge 0$.

**Inflection fraction.** The rescaling is affine, so the inflection
stays at $t^\ast = \ln(d)/r$. The share of the ceiling reached there
is $(e^{-1} - e^{-d})/(1 - e^{-d})$, which is 0.364 at $d = 5$
(derived; verified, `check-gompertz-rmanova-cycles.R`). A logistic
curve inflects at 0.5, so the Gompertz completes less of its rise
before the inflection.

**The displacement delays, it does not steepen.** The maximum slope
on $t \ge 0$ is $mr/[e(1 - e^{-d})]$ for $d \ge 1$, which is
strictly decreasing in $d$ (derived; verified numerically). At the
shipped drug-response values the steepest rise is 1.71 points per
week at week 3.83. An earlier version of this document said a large
$d$ gives "a long slow start and a steep rise"; only the lag
statement is correct (corrected 2026-10-08).

### Shipped parameters and landmarks

<!-- From CM Table 1 and Section "Shipped parameters" -->

The package dataset `extracted_rp` (`data/extracted_rp.rda`) gives
the Hendrickson et al. calibration:

| Component | $m$ | $d$ | $r$ (per week) | SD |
|---|---|---|---|---|
| $BR$ | 10.986 | 5 | 0.42 | 8 |
| $TV$ | 6.506 | 5 | 0.35 | 10 |
| $PB$ | 6.506 | 5 | 0.35 | 10 |

$TV$ and $PB$ share parameters because the source trial fitted only
their sum; separating them requires an open-label-placebo arm or a
direct expectancy measurement (source: `docs/40`, Section 3).

Landmarks of the shipped curves (both at $d = 5$; verified,
`check-landmarks-hybrid-ceilings.R`):

| Curve | Inflection (weeks) | Half-ceiling (weeks) | 90% ceiling (weeks) | Share at week 8 |
|---|---|---|---|---|
| $BR$ ($r = 0.42$) | 3.83 | 4.73 | 9.21 | 0.839 |
| $TV$, $PB$ ($r = 0.35$) | 4.60 | 5.67 | 11.05 | 0.736 |

A four-week on-drug crossover period recovers 0.390 of the BR
ceiling; one week on drug recovers 0.031. Half the drug-response
ceiling is reached at 33 days (verified).

### Participant heterogeneity is level-only

<!-- From CM "Heterogeneity is level-only" -->

In every current data-generating process the parameters $(m, d, r)$
of each component are population constants. Participants differ only
through the multivariate-normal deviation of each component around
its mean curve, with a constant per-visit standard deviation (scaled
by $\eta_t$ for $PB$) and the serial correlation, plus the
architecture's biomarker term. No participant has a different ceiling
or onset rate (`buildSigma()` and `generateData()`, inspected).

### Alternative response-curve families

<!-- From CM "Alternative response-curve families" and docs/25 -->

Paper 07 evaluates three alternatives, each calibrated to match the
Gompertz reference at its asymptote and at $t = 5$ weeks:

- **Logistic:** $y(t) = m / [1 + \exp(-r(t - t_0))] - \text{offset}$
- **Hyperbolic tangent:** $y(t) = m \cdot [\tanh(r(t - t_0)) + 1]/2
  - \text{offset}$
- **Piecewise-linear breakpoint:** $y(t) = m \cdot \min(1,
  t/t_\text{bp}) + \text{post\_slope} \cdot \max(0, t - t_\text{bp})$

The hyperbolic-tangent form is a rescaled logistic: $(\tanh(x) +
1)/2 = 1/(1 + e^{-2x})$, so it is equivalent to the logistic at
rate $2r$ (derived).

The pharmacodynamic Hill (Emax) family is:

$$
y(t) = E_\text{max} \cdot \frac{t^{n_H}}{EC_{50}^{n_H} + t^{n_H}}
$$

For $n_H > 1$ it has an inflection at $t^\ast = EC_{50}((n_H - 1)/
(n_H + 1))^{1/n_H}$, where it has completed $(n_H - 1)/(2n_H)$ of
its rise (derived; verified, `check-gompertz-rmanova-cycles.R`). The
Hill family is not among the package's `br_family` options and was
not evaluated by Paper 07.

### Trajectory shape and the architectures (Paper 07)

<!-- From CM "How trajectory shape interacts with carryover", lines
1527-1583 -->

Under mean moderation Paper 07 uses a trajectory-scaled form whose
biomarker-dependent component rides the response curve
($\text{moderation\_scaling = 'trajectory'}$). A family with a
slower early rise therefore places less of the interaction signal on
the earliest on-drug visits, while the marginal interaction effect
size is held fixed.

Under covariance moderation the interaction lives in the
biomarker-response correlation $\mathrm{Cor}(B_i, BR_{it}) =
c_{bm} g_t$, whose coupling profile $g_t$ is set by exposure state
and $t_{sd}$, not by the mean curve. Trajectory-shape
misspecification cannot move the covariance-channel interaction test;
Paper 07 confirms this with a matched-seed check. Paper 07's
16-cell factorial (four families, two architectures, two biomarker
levels, $N = 70$, OL+BDC, $t_{1/2} = 0$) finds power differences of
only a few percentage points across families under both architectures.

### New-to-write items

<!-- TODO: a family-overlay figure showing Gompertz, logistic, tanh
and piecewise-linear calibrated at the same asymptote and t=5 anchor
-->

<!-- TODO: 25 Section 6 and Recommendation 4 rewritten to the
corrected account: the displacement delays the rise but does not
steepen it (Eq. invc-maxslope); the logistic inflects at 50% of its
ceiling, the Gompertz at ~36%, so the Gompertz completes less of its
rise before the inflection -->

<!-- TODO: the link from 25 Recommendation 5 to P04's v_i*Psi_it:
the trajectory-scaled mean-moderation shift is
beta_bm * b_i * sigma_BR * mu_BR(t) / mean(mu_BR, on-drug), and
paper 04's sensitivity-model notation writes this as v_i * Psi_it -->

---

## Derivations

<!-- one row per result: statement | CM label | verification script -->

| Statement | CM label / source | Script |
|---|---|---|
| Inflection at $t^\ast = \ln(d)/r$ | `eq:invc-share-inflection` | `check-gompertz-rmanova-cycles.R` |
| Inflection fraction $(e^{-1} - e^{-d})/(1 - e^{-d})$ | `eq:invc-share-inflection` | `check-gompertz-rmanova-cycles.R` |
| Max slope decreasing in $d$ | `eq:invc-maxslope` | `check-gompertz-rmanova-cycles.R` block 2 |
| Landmarks of shipped curves | Table `tab:invc-landmarks` | `check-landmarks-hybrid-ceilings.R` block 1 |
| Rate rescaling: $G(t; m, d, kr) = G(kt; m, d, r)$ | CM text | verified at $k=2$ |

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

See landmark table in "Setting and argument". Paper 07 simulation
results (power differences $< 5\%$ across families) are in the
paper 07 report (`analysis/report/07-gompertz-evaluation/report.Rmd`).

---

## Relation to the manuscripts

- **Paper 01:** uses $G(t; m, d, r)$ for BR, TV, PB with shipped
  parameters; landmarks stated in Methods.
- **Paper 05:** uses the same parameterization; the modified form is
  stated as Eq. (1).
- **Paper 06:** writes the asymptote-unrescaled form
  $C(t) = m_C \exp(-d_C \exp(-r_C t)) - m_C \exp(-d_C)$; asymptote
  is $m_C(1 - e^{-d_C})$, not $m_C$ (a 0.67% difference at $d = 5$).
- **Paper 07:** evaluates four families; see "Trajectory shape"
  subsection above.

---

## Limitations and open items

- Participant heterogeneity in Gompertz parameters (individual onset
  rates, ceilings) is not implemented; all heterogeneity is
  level-only. This is a known structural limitation of every current
  DGP.
- The numerical instability of the rescaling for $d \to 0$
  (`expm1` form needed for $d < 10^{-6}$; see CM
  `eq:invb-gompertz-stable`).
- Paper 06's biological-response section still states that a high
  $d_{BR}$ gives an early-onset curve (incorrect; corrected in the
  compendium but not yet in the paper source).

---

## Reproduce

```bash
Rscript analysis/scripts/compendium-mathematics/\
  check-gompertz-rmanova-cycles.R
Rscript analysis/scripts/compendium-mathematics/\
  check-landmarks-hybrid-ceilings.R
```

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | All | Phase 2 content move from docs/24, docs/25, docs/40, CM sec:gompertz | Whitepaper refactor |
