# Component decomposition, identification and reduction

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Identification of BR, PB, TV by Hybrid design; rank of
the phase-augmented model; reductions.
**Prerequisites:** ch05, ch07

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
| $\Delta_{PB}$ | $E[PB \mid \text{on}] - E[PB \mid \text{off}]$ (Paper 06 notation) |
| $\Delta_{TV}$ | natural-history contribution that fails to cancel |
| $w_{PB}$, $w_{TV}$ | design-determined loadings on placebo/natural-history slope |
| $\beta_{bm}^C$ | $\mathrm{Cov}(m_C, B) / \mathrm{Var}(B)$ for component $C$ |
| $\gamma$ | $BR \times PB$ subadditivity parameter |
| $c_{bm,PB}$ | correlation of biomarker with individual PB level |

---

## Summary

<!-- evidence-labeled results go here -->

Paper 06 asks what a one-component ("lumped") analysis of the treatment
effect and its biomarker interaction actually estimates when the
generative process has three distinguishable components. The
probability limit of the lumped treatment estimate is:

$$
\operatorname*{plim}\hat\beta_D = m_{BR} + \Delta_{PB} + \Delta_{TV},
$$

a systematic upward bias from placebo and natural-history contributions.
For the biomarker slope:

$$
\beta_{bm}^{\mathrm{lumped}} = \beta_{bm}^{BR} + w_{PB}\,\beta_{bm}^{PB}
  + w_{TV}\,\beta_{bm}^{TV},
$$

an identity that cannot be eliminated by increasing sample size. A
biomarker correlated with placebo responsiveness registers a nonzero
lumped slope even when its true pharmacological association is exactly
zero.

Recovery of $\beta_{bm}^{BR}$ is design-dependent: under the
compendium's standard Hybrid design, decomposition analyses fail to
recover the true pharmacological slope because $D_{bc}$ and the
expectancy multiplier $\eta$ are positively correlated ($\approx +0.6$).
A balanced-placebo design (Paper 06) reduces that correlation to
$\approx -0.66$ and restores recovery (verified, study B production run).

The phase-augmented analysis formula is rank-deficient by one whenever
the open-label phase is entirely on drug; `nlme::lme` cannot fit the
full formula, and the driver refits without `Dbc:phase` (inspected).

<!-- TODO: sections 18.1, 18.3 and 22.2 of docs/31 brought up to date:
the prognostic couplings of the stress-test driver exist; F1 carries
De; these change the contamination story from the original docs/31
account -->

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3.3 for this chapter:
- Sources: 24 identification and "Why the analysis does not match";
  31 Parts II, III (paper 06 studies, totals), IV, V; 28 Sections 3
  to 4 (2026-10-08 corrections); CM "Partitioning the treatment
  effect" (4 subsections); papers 06, 13, 14.

- New to write: 31 Sections 18.1, 18.3 and 22.2 brought up to date
  (the prognostic couplings of the stress-test driver exist; F1
  carries De).
-->

### The outcome decomposition

<!-- From CM sec:partition, lines 4076-4127 -->

The outcome decomposes as:

$$
Y_{it} = \mathrm{BL}_i - \big[BR_{it} + PB_{it} + TV_{it}\big] +
  \varepsilon_{it}.
$$

The placebo-belief component $PB$ is scaled at each visit by the
expectancy multiplier $\eta(\mathrm{phase})$:

$$
PB(t \mid \mathrm{phase}) = \eta(\mathrm{phase})\cdot m_{PB}
  \exp\!\big(-d_{PB}\exp(-r_{PB}t)\big) + \text{(offset)},
$$

$$
\eta = \begin{cases} 1.0 & \text{open-label, on drug} \\
  \approx 0.5 & \text{blinded discontinuation} \\
  \in(0,0.5) & \text{open-label placebo.}
\end{cases}
$$

Because $\eta$ scales both the mean and SD of $PB$ by the same factor,
only the product $\eta\cdot m_{PB}$ is identified from a two-phase
design; recovering $m_{PB}$ itself requires an open-label-placebo phase
or a direct expectancy measurement.

A generalized non-additive form allows a $BR \times PB$ interaction:

$$
Y_{it} = \mathrm{BL}_i - \big[BR_{it} + PB_{it} + TV_{it} +
  \gamma\cdot BR_{it}\cdot PB_{it}\big] + \varepsilon_{it},
$$

with subadditivity defined by $\partial Y/\partial(BR\cdot PB) > 0$ on
the reduction scale. Paper 06 distinguishes true $\gamma \ne 0$
subadditivity from apparent subadditivity that a latent responder-class
gating function could generate through confounding.

### What a one-component analysis actually estimates

<!-- From CM sec:partition, lines 4128-4221 -->

Writing within-participant symptom reduction as $R_{it} = \mathrm{BL}_i
- Y_{it}$ and fitting a one-component model, Paper 06 derives:

$$
\operatorname*{plim}\hat\beta_D = m_{BR} + \Delta_{PB} + \Delta_{TV},
$$

where $\Delta_{PB} = E[PB \mid \text{on}] - E[PB \mid \text{off}]$ and
$\Delta_{TV}$ is the natural-history contribution that fails to cancel
under the contrast used (zero for a clean on/off within-participant
contrast; $E[TV]$ itself for a baseline-anchored pre-post analysis).
Both are generally non-negative for an improving outcome, so the
one-component estimate systematically overstates $m_{BR}$. This is a
property of the estimator's probability limit, not a finite-sample
artifact.

For the biomarker slope:

$$
\beta_{bm}^{\mathrm{lumped}} = \beta_{bm}^{BR} + w_{PB}\,\beta_{bm}^{PB}
  + w_{TV}\,\beta_{bm}^{TV}, \qquad
\beta_{bm}^{C} = \frac{\operatorname{Cov}(m_C, B)}{\operatorname{Var}(B)},
$$

with $w_{PB}, w_{TV} \ge 0$ design-determined loadings. A biomarker
correlated with placebo responsiveness therefore registers a nonzero
lumped slope even when its true pharmacological association
$\beta_{bm}^{BR}$ is exactly zero, and because this is an identity, no
additional sample size can separate a genuinely pharmacological signal
from a confounded one whenever $w_{PB}\beta_{bm}^{PB} +
w_{TV}\beta_{bm}^{TV} \ne 0$.

The operative confound is a biomarker's correlation with
individual-level variance in placebo responsiveness, not with its mean:

$$
\text{contamination driver} = c_{bm,PB}\cdot\sigma_{PB}, \qquad
  \text{not } c_{bm,PB}\cdot m_{PB}.
$$

Paper 06's calibration target is $\beta_{bm:D} = -c_{bm}\sigma_{BR}/
\sigma_{bm} = -2.25$. Study A (uncontaminated biomarker, $c_{bm,PB} = 0$)
shows bias within one Monte Carlo standard error of zero at every
$N \in \{35, 70, 200, 300\}$ (totals; verified from stored summary).
A contaminated pilot at $c_{bm,PB} \in \{0, 0.15, 0.30, 0.45\}$ shows
bias at $N = 300$ of $-0.06$, $-0.12$, $+0.17$, $+0.51$ (the last
about 30 MCSEs from zero; verified). Paper 06 explains the positive sign
at high contamination as attenuation toward zero: the elevated placebo
response of high-$c_{bm,PB}$ participants is partially absorbed by the
biomarker main effect.

The bias lives in the location of the estimand, not in the test: a
significance test of $\hat\beta_{bm:D} = 0$ cannot detect a
displacement away from the null. Detecting it requires comparing the
estimate with $\beta_{bm}^{BR}$, which only a decomposition or a phase
contrast makes available.

### Recovery, and its dependence on trial design

<!-- From CM sec:partition, lines 4223-4261 -->

A phase-augmented analysis model:

```
Sx ~ bm + t + Dbc + phase + Dbc:phase + bm:Dbc + bm:Dbc:phase
```

with random intercept and `corCAR1` residual correlation, where
`phase` (open-label / blinded / post-crossover) absorbs systematic
expectancy differences.

Whether this, or full component-matched Gompertz-basis fitting,
recovers $\beta_{bm}^{BR}$ is design-dependent. Under the compendium's
standard Hybrid design, belief-covariate decomposition, blinded-stratum
contrasts, and full Gompertz-basis decomposition all fail to recover
$\beta_{bm}^{BR}$, performing no better (and sometimes worse) than the
lumped analysis, because $D_{bc}$ and the expectancy multiplier $\eta$
are positively correlated ($\approx +0.6$) under that design.

A balanced-placebo design (Paper 06: three open-label, six covert
on-drug with $\eta = 0$, and six open-placebo with $\eta = 1$
timepoints) reduces that correlation to $\approx -0.66$ and restores
recovery. In the production run (1,000 replicates per cell, single-path
design, $N \in \{70, 150\}$ totals), decomposition-analysis bias at
$N = 150$ is $-0.022$, $+0.008$, $+0.021$, $-0.026$ across
$c_{bm,PB} \in \{0, 0.1, 0.2, 0.3\}$ (MCSE about 0.02, so within 1.4
MCSE of zero; verified from `study-b-balanced-summary.rds`), while the
lumped estimator's bias grows monotonically: $-0.014$, $+0.161$,
$+0.340$, $+0.476$.

### Rank of the phase-augmented model

<!-- From CM sec:invb-collinearity, lines 4263-4358 -->

The phase-augmented formula is rank-deficient whenever a phase is
entirely on drug. If every visit of the open-label level is on drug,
$D_{bc,it} = 1$ there, so $D_{bc,it}$ is an exact linear combination
of the intercept, the phase main effects, and the `Dbc:phase` columns.
The fixed-effects design matrix therefore loses one rank, whatever the
coding of the factor (derived).

Consequences:

- With three phases (open label, blinded discontinuation, crossover),
  the extension adds six columns but only five identified parameters,
  not the six that `docs/24` counts. In the two-level coding of the
  Paper 06 driver (phase set by $\eta_t \ge 0.99$), three columns add
  two identified parameters. (Verified by the rank of the model matrix
  for the 8-visit Hybrid paths and for the 16-visit two-path design
  of the Paper 06 recovery study, `check-gompertz-rmanova-cycles.R`.)
- `nlme::lme` cannot fit the full formula: on the 16-visit design it
  stops with "Singularity in backsolve" (verified). The driver catches
  the error and refits without `Dbc:phase`, so every phase-augmented
  result is from the reduced formula and is recorded with
  `formula_dropped = 'Dbc:phase'` (inspected).
- In the reduced formula `Sx ~ ... + bm:Dbc + bm:Dbc:phase`, the
  biomarker slope is $\beta_B + \beta_{bm:D}$ at open-label visits and
  $\beta_B + (\beta_{bm:D} + \beta_{bm:D:\mathrm{BDC}})D_{bc,it}$ at
  blinded visits. The reported coefficient $\beta_{bm:D}$ is therefore
  the open-label slope minus the blinded slope extrapolated to
  $D_{bc} = 0$: an open-label versus blinded contrast in which drug
  state and expectancy change together, while the blinded-window
  interaction is the sum $\beta_{bm:D} + \beta_{bm:D:\mathrm{BDC}}$.

The drug-expectancy correlation by design, for a single-path design
with $n_1$ on-drug full-belief, $n_2$ on-drug covert ($\eta = 0$), and
$n_3$ off-drug full-belief visits:

$$
\mathrm{Cor}(D, \eta) = \frac{p_{11} - p_D p_\eta}
  {\sqrt{p_D(1 - p_D)\,p_\eta(1 - p_\eta)}}, \qquad
p_D = \frac{n_1 + n_2}{n}, \quad p_\eta = \frac{n_1 + n_3}{n},
\quad p_{11} = \frac{n_1}{n}.
$$

For the balanced-placebo design of Paper 06 ($n_1 = 3$,
$n_2 = n_3 = 6$) this is exactly $-2/3$; with decayed $D_{bc}$ at the
driver's $t_{1/2} = 1$ week and 1.5-week spacing it is $-0.656$, which
matches Paper 06's "$\approx -0.66$" (derived; verified,
`check-gompertz-rmanova-cycles.R`).

### New-to-write items

<!-- TODO: docs/31 Sections 18.1 and 18.3 brought up to date: the
stress-test driver includes prognostic couplings of the biomarker with
PB and TV (c_{bm,PB} and c_{bm,TV} in the DGP); state whether the
current F1 analysis (which carries De and bm:De) is adequate to detect
contamination from prognostic couplings or whether the decomposition is
still needed -->

<!-- TODO: docs/31 Section 22.2 brought up to date: F1 carries De and
bm:De (the expectancy term and its biomarker interaction); explain how
this changes the identification problem and whether the lumped-slope
identity of eq:lumped-slope still applies unchanged -->

<!-- TODO: the two-component reductions (BR+PB and BR+TV) and their
power implications: which reduced models remain adequate for detecting
bm:Dbc; sourced from docs/31 Parts IV-V -->

---

## Derivations

| Statement | Source | Script |
|---|---|---|
| $\operatorname*{plim}\hat\beta_D = m_{BR} + \Delta_{PB} + \Delta_{TV}$ | CM `eq:plim-betaD`, Paper 06 | `check-gompertz-rmanova-cycles.R` |
| Lumped-slope identity | CM `eq:lumped-slope`, Paper 06 | derived |
| Phase-augmented formula loses one rank when OL is entirely on drug | CM `eq:invb-alias` | `check-gompertz-rmanova-cycles.R` |
| $\mathrm{Cor}(D, \eta) = -2/3$ for balanced-placebo | CM `eq:invb-cor-de` | `check-gompertz-rmanova-cycles.R` |

---

## Numerical evidence

From the Paper 06 production run (verified from stored summaries):

| Analysis | $c_{bm,PB}$ | Bias at $N = 150$ | MCSE |
|---|---|---|---|
| Decomposition (balanced-placebo) | 0 | $-0.022$ | 0.02 |
| Decomposition (balanced-placebo) | 0.1 | $+0.008$ | 0.02 |
| Decomposition (balanced-placebo) | 0.2 | $+0.021$ | 0.02 |
| Decomposition (balanced-placebo) | 0.3 | $-0.026$ | 0.02 |
| Lumped (Hybrid) | 0 | $-0.014$ | |
| Lumped (Hybrid) | 0.1 | $+0.161$ | |
| Lumped (Hybrid) | 0.2 | $+0.340$ | |
| Lumped (Hybrid) | 0.3 | $+0.476$ | |

---

## Relation to the manuscripts

- **Paper 06:** defines the phase-augmented model; studies A and B;
  demonstrates design-dependent recovery.
- **Paper 13:** tests the lumped-slope identity in a within-participant
  decomposition.
- **Paper 14:** two-component reduction studies; TV-only identification.

---

## Limitations and open items

- The stress-test driver includes prognostic couplings ($c_{bm,PB}$,
  $c_{bm,TV}$); their interaction with F1's expectancy term `De` has
  not been evaluated.
- The lumped-slope identity under F1's additional `De` and `bm:De`
  terms (compared with the original analysis model) has not been
  derived.
- Sections 18.1, 18.3 and 22.2 of `docs/31` are not yet updated to
  account for these two changes (open item; see TODO above).

---

## Reproduce

```bash
Rscript analysis/scripts/component-decomposition/\
  01-study-a-run.R
Rscript analysis/scripts/component-decomposition/\
  02-study-b-run.R
Rscript analysis/scripts/compendium-mathematics/\
  check-gompertz-rmanova-cycles.R
```

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | All | Phase 2 content move from CM sec:partition, docs/24, docs/28, docs/31 Parts II-V | Whitepaper refactor |
