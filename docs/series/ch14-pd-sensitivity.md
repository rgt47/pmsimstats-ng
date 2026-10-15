# Beyond the refined process: a pharmacodynamic sensitivity model

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** Pharmacodynamic turnover model; sensitivity binding; consequences for the test.  
**Prerequisites:** ch09, ch10

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

<!-- table: process, carryover, coupling, allocation, design, N, analysis ID,
t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 38 (spine, renamed $f(t) \to h(t)$, $F \to \Delta_h$); 44 Steps 7
  to 8; CM "A pharmacodynamic sensitivity model" (4 subsections); 36 Section 6.7.

- New to write: $\sigma_v^2$ rewritten as $\tau^2$ (DQ-6, decided).
-->

<!-- Spine: docs/38-pharmacodynamic-carryover-model.md (copied verbatim; Phase 2) -->

# Modeling Carryover from Pharmacodynamic Principles: A Generative Proposal for the Mean, Moderation and Covariance of the Drug Response {.unlisted .unnumbered}
*2026-10-02 17:31 PDT*

**Author.** pmsimstats team

**Purpose.** This paper proposes, from pharmacodynamic principles, a
generative model of the drug response in which carryover is specified
once and the mean, the biomarker moderation and the covariance of the
response all follow from it. It sets out the model, derives its implied
moments and the consequences for the interaction test, and compares it
with the constructs in use. It is a proposal: its results are derived,
not implemented or simulated.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

The notation follows `analysis/report/NOTATION.md`. Symbols marked
*canonical* carry the meaning given there; symbols marked *local* are
introduced here and chosen not to collide with canonical symbols. The
proposal's sensitivity, moderation-strength and noise symbols were
renamed on 2026-10-07 ($\theta$, $\kappa$, $u$, $\varepsilon$, $\gamma$
in the first draft) because those letters carry other canonical
meanings. On 2026-10-09 the contrast-biomarker correlation was renamed
from $\rho_{\Delta B}$ to $r_{\Delta B}$, as in `docs/37` and
`analysis/report/compendium-mathematics.tex` ($\rho$ is the canonical
serial correlation), and the normal law written $\mathcal{N}$.

| Symbol | Meaning | Status |
|----------|----------------------------------------|------|
| $i$, $t$, $s$ | patient; visit (two visits $t$, $s$) | canonical |
| $B_i$, $b_i$ | biomarker; standardized biomarker $(B_i - \mu_{bm})/\sigma_{bm}$ | canonical |
| $\mu_{bm}$, $\sigma_{bm}$ | biomarker mean and standard deviation | canonical |
| $BR_{it}$ | drug (biological) response component | canonical |
| $Y$ | symptom score, the outcome | canonical |
| $c_{bm}$ | covariance-moderation strength of the published construct | canonical |
| $D_{bc}$ | exposure-decayed continuous drug indicator | canonical |
| $\beta_{bm:D}$ | interaction slope: population slope of the contrast $\Delta$ on $B$ | canonical |
| $t_{sd}$, $t_{1/2}$ | time since discontinuation; carryover half-life as used by the published construct | canonical |
| $t_{1/2}^{R}$ | response (pharmacodynamic offset) half-life | local |
| $k_{\text{out}}$ | first-order turnover rate of the response, $\ln 2 / t_{1/2}^{R}$ | local |
| $f(t)$ | population drug-response time course, scaled to 1 at full on-drug effect | local |
| $f^{\ast}(t)$ | drug-driven target of the turnover model: 1 on drug, 0 off drug | local |
| $t_{\text{stop}}$ | time of discontinuation | local |
| $\omega_i$, $\bar\omega$ | patient $i$'s drug sensitivity and its population mean | local |
| $\zeta$ | strength of multiplicative moderation of sensitivity by $b_i$ | local |
| $v_i$, $\sigma_v$ | sensitivity variation not explained by the biomarker, and its SD | local |
| $\xi_{it}$, $\sigma_\xi$, $r_\xi(t,s)$ | serial within-patient noise of $BR$, its SD and its correlation | local |
| $c_{bm}^{\text{peak}}$ | biomarker-response correlation at full effect, a derived value | local |
| $a_t$ | E9 contrast weights, $1/n_{\text{on}}$ on drug, $-1/n_{\text{off}}$ off drug | canonical |
| $\Delta$, $\Delta_i$ | patient's on-minus-off contrast (E9) | canonical |
| $F$ | contrast-weighted exposure, $\bar f_{\text{on}} - \bar f_{\text{off}}$ | local |
| $\bar f_{\text{on}}$, $\bar f_{\text{off}}$ | mean of $f(t)$ over on-drug and off-drug visits | local |
| $V_0$ | contrast variance from serial noise and the other response factors | local |
| $r_{\Delta B}$ | correlation between contrast and biomarker | canonical |
| $\mathrm{BL}$, $TV$, $PB$ | baseline level; natural-course and placebo-belief components | canonical |
| $n_{\text{on}}$, $n_{\text{off}}$ | on-drug and off-drug visits of a path | canonical |
| $G$, $H$ | shorthand $\bar\omega\zeta$ and $\bar\omega^2\zeta^2 + \sigma_v^2$ in the proof of Section 6 | local |

Table: Notation index for the pharmacodynamic carryover model

**Glossary.**

- **Pharmacokinetic carryover.** Drug remaining in the body after the
  last dose.
- **Pharmacodynamic (disease-state) carryover.** Improvement that
  persists after the drug has been eliminated.
- **Turnover (indirect-response) model.** A response that relaxes
  toward a drug-driven target at a first-order rate.
- **Anchored decay.** Off-drug decay measured from the level reached at
  discontinuation.
- **Moderation profile.** The fraction of the full biomarker moderation
  present at a visit (`docs/37`).
- **Step rule, graded coupling, decayed mean moderation.** The coupling
  rules of `docs/37`, Section 2.
- **E9.** The paired-difference interaction statistic of `docs/37`,
  Section 3.
- **`orig`, `covar`, configurations C and D.** The published
  (`58b32a9`) construct, the package's AR(1) construct, and the
  configurations of `docs/36`, Appendix A.6.

## 1. Summary

`docs/36` and `docs/37` traced both problems of the Hendrickson et al.
(2020) simulation framework, invalid correlation matrices and the
collapse of power under small carryover, to the way the framework
specifies its correlation matrix: entry by entry, with the carryover
half-life acting on the off-drug mean through one formula and on the
biomarker correlation through a separate rule. This paper argues, from
pharmacodynamic principles, that carryover should not be specified
separately for the mean and for each element of the correlation
matrix. It should be specified once, in a generative model of each
patient's drug response, from which the mean, the variance, the
biomarker correlation and the cross-visit correlations all follow.

The proposal:

1. **Carryover is pharmacodynamic.** For the prazosin and
   posttraumatic stress disorder (PTSD) context, the drug leaves plasma
   with a half-life of approximately 2.3 hours (Hobbs et al., 1978), so
   pharmacokinetic carryover is essentially complete within a day.
   Persistence of benefit beyond that is a property of the response,
   and the half-life parameter should be the response (offset)
   half-life (Section 2).
2. **The mean follows a turnover model.** The response relaxes toward a
   drug-driven target at a first-order rate (Dayneka et al., 1993).
   Off drug this gives decay anchored at the level reached at
   discontinuation, and on restart the response builds from its
   residual level (Section 3).
3. **Moderation scales with the response.** A predictive biomarker
   acts on drug sensitivity, so the biomarker-dependent part of the
   response rises and decays with the response itself (Section 4).
4. **The covariance is implied, not specified.** Writing each patient's
   drug response as a random sensitivity times the population time
   course, plus serial noise, makes the covariance a sum of positive
   semidefinite terms. It is valid for every design and half-life by
   construction, the biomarker correlation follows the response, and
   the drug-response variance shrinks as the effect washes out
   (Section 5).
5. **The closed-form results of `docs/37` extend.** The interaction
   slope is proportional to the contrast-weighted exposure
   $F = \bar f_{\text{on}} - \bar f_{\text{off}}$, and power is strictly
   increasing in $F$; unexplained variation in drug sensitivity caps
   the attainable power (Section 6).

The paper is a proposal. Sections 3 to 6 are derivations, not
simulation results; nothing in it has been implemented or simulated.

## 2. Two kinds of carryover

**Pharmacokinetic carryover** is drug remaining in the body after the
last dose. For prazosin it is short: in 24 normotensive subjects, drug left
plasma with a half-life of approximately 2.3 hours (Hobbs et al.,
1978). Five half-lives, about 12 hours, remove 97% of it. Hendrickson
et al.'s half-lives of 0.1 and 0.2 weeks (17 and 34 hours) are on this
scale or slightly above it.

**Pharmacodynamic, or disease-state, carryover** is improvement that
persists after the drug has gone: consolidated sleep, fewer
nightmares, behavioral learning, recovery of a physiological set
point. Its time course is set by the turnover of the response, not by
drug elimination, and when the response is slow to turn over it can
last days to weeks. Effect-compartment models (Sheiner et al., 1979)
and indirect-response models (Dayneka et al., 1993) are the standard
quantitative treatments of this delay between concentration and
effect; in indirect-response models in particular, the response
returns to baseline slowly after the drug is gone.

For a symptom outcome measured weekly or less often, pharmacokinetic
carryover is negligible and the carryover that matters is
pharmacodynamic. The half-life parameter of the simulation should
therefore be read, and chosen, as the response half-life
$t_{1/2}^{R}$. The appropriate value for prazosin in PTSD is not known
to us; it is an empirical quantity that the simulation should treat as
uncertain.

## 3. The mean: a turnover model

Let $f(t) \in [0, 1]$ denote the population drug-response time course,
scaled so that $f = 1$ at the full on-drug effect. A minimal turnover
model is

$$
\frac{df}{dt} = k_{\text{out}}\,\big(f^{\ast}(t) - f(t)\big), \qquad
k_{\text{out}} = \frac{\ln 2}{t_{1/2}^{R}},
$$

where $f^{\ast}(t)$ is the drug-driven target: the full on-drug effect
(1 on this scale) while on drug, 0 off drug. For a target that is constant within each on-drug
or off-drug interval, the solution is exact and piecewise: off drug,

$$
f(t) = f(t_{\text{stop}})\; 2^{-t_{sd}/t_{1/2}^{R}},
$$

decaying from the level reached at discontinuation; on drug, $f$
approaches the target from wherever it started. Three properties
follow, each of which differs from the published construct:

- **Anchored decay.** Off-drug decay is measured from the level at
  discontinuation. The published code multiplies the previous visit's
  already-decayed mean by $2^{-t_{sd}/t_{1/2}}$ with $t_{sd}$
  cumulative, which decays faster than the stated half-life (`docs/36`;
  for off-drug visits one and two weeks after stopping, the second is
  at $2^{-3/t_{1/2}}$ rather than $2^{-2/t_{1/2}}$).
- **Restart from the residual.** When the drug is restarted the
  response builds from its residual level. The published code restarts
  the Gompertz onset curve from zero, as if the earlier course had not
  occurred.
- **One mechanism.** Onset and offset come from the same process. If a
  Gompertz onset is retained for fidelity to the published trajectories,
  the offset should at least be anchored and the restart should begin
  from the residual.

The simplest change compatible with the published model is to keep
the Gompertz onset, anchor the offset, and start each new on-drug
interval from the residual level. The full turnover model is the
principled version.

## 4. Moderation scales with the response

A predictive biomarker modifies the drug's effect: in pharmacological
terms it acts on sensitivity, the maximal effect or the potency of the
drug in that patient. Write patient $i$'s drug response as

$$
BR_{it} = \omega_i\, f(t) + \xi_{it}, \qquad
\omega_i = \bar\omega\,(1 + \zeta b_i) + v_i ,
$$

with $b_i$ the standardized biomarker, $\zeta$ the strength of
moderation, $v_i \sim \mathcal{N}(0, \sigma_v^2)$ the variation in sensitivity
that the biomarker does not explain, and $\xi_{it}$ serial
within-patient noise with variance $\sigma_\xi^2$, independent
of $b_i$ and $v_i$. The mean response is $\bar\omega f(t)$, and the
biomarker-dependent part, $\bar\omega\zeta b_i f(t)$, has the same
time course as the response.

The *moderation profile*, the fraction of the full moderation present
at visit $t$, is therefore $f(t)$ itself:

- **On drug**, moderation grows with the response. It is not at full
  strength from the first on-drug visit; it is small early in
  titration, when the response is small.
- **Off drug**, moderation decays with the response. It neither holds
  at full strength while the drug effect washes out, as under the
  published step rule (`docs/37`, Section 5.1), nor stops the moment
  the drug does, as under the current mean-moderation construct
  (configuration D of `docs/36`), which shifts the response at on-drug
  visits only.

The decayed mean moderation of `docs/37` (Section 5.6), with an
anchored profile equal to 1 on drug and $2^{-t_{sd}/t_{1/2}}$ off drug,
is an approximation to this; the proportional profile $f(t)$ is the
form the pharmacology implies.

## 5. The implied covariance

Under the model of Section 4 every second moment of the drug response
follows from four quantities, $\bar\omega\zeta$, $\sigma_v$,
$\sigma_\xi$ and the serial correlation $r_\xi$ of $\xi$
(derived):

$$
\begin{aligned}
\mathrm{Var}(BR_t) &= (\bar\omega^2\zeta^2 + \sigma_v^2)\, f(t)^2 + \sigma_\xi^2, \\
\mathrm{Cov}(BR_t, BR_s) &= (\bar\omega^2\zeta^2 + \sigma_v^2)\, f(t) f(s)
  + \sigma_\xi^2\, r_\xi(t, s), \\
\mathrm{Cov}(B, BR_t) &= \sigma_{bm}\,\bar\omega\zeta\, f(t), \qquad
\mathrm{Cor}(B, BR_t) = \frac{\bar\omega\zeta\, f(t)}
  {\sqrt{(\bar\omega^2\zeta^2 + \sigma_v^2) f(t)^2 + \sigma_\xi^2}} .
\end{aligned}
$$

Four consequences:

- **The biomarker correlation follows the response.** It is zero where
  $f = 0$, rises and decays with $f$, and is never a step. The
  quantity the published framework calls $c_{bm}$ becomes a derived
  value, most naturally the correlation at full effect,
  $c_{bm}^{\text{peak}} = \bar\omega\zeta /
  \sqrt{\bar\omega^2\zeta^2 + \sigma_v^2 + \sigma_\xi^2}$.
- **The drug-response variance shrinks off drug.** As the effect washes
  out, the sensitivity term vanishes and only the serial noise remains.
  The published construct gives the drug-response factor its full
  standard deviation (8) at every visit, including off-drug visits with
  zero mean, as if patients had a full-sized random drug response with
  no drug. Of the construct's features this is, in our judgment, the
  least defensible, and it adds noise to every on-off contrast.
- **Cross-visit correlation has structure.** It is a rank-one
  sensitivity term, $f(t)f(s)$, plus the serial correlation of the
  noise. For symptom data a decaying serial correlation (AR(1) or
  continuous-time AR(1)) is more plausible than compound symmetry.
- **Validity by construction.** The covariance of the drug responses is
  a sum of a positive semidefinite rank-one term and a positive
  definite noise covariance, so it is positive definite for every
  design, every schedule and every half-life. There is no ceiling on
  $c_{bm}$ to exceed and nothing to repair: an impossible combination
  of correlations cannot be written down, because correlations are
  outputs rather than inputs.

The other response factors follow the same logic. The time-related
(natural history) and expectancy (placebo) factors are separate
processes; correlation between them and the drug response should
arise from shared patient-level effects, such as a patient's general
responsiveness entering both $\omega_i$ and the expectancy effect,
rather than from constants placed in matrix cells. The resulting
covariance is again a sum of positive semidefinite terms, in the form
of the linear model of coregionalization of `docs/35`.

Data from this model need no correlation matrix at all: draw $b_i$,
$v_i$ and the serial noise, compute $f(t)$ for the patient's schedule,
and assemble $BR_{it}$. The construction is the factored generator of
`docs/36`, Section 6.7, without the matrix.

## 6. Consequences for the interaction test

The closed-form approach of `docs/37` extends directly (derived; not
verified by simulation). With the paired-difference statistic E9 and
contrast weights $a_t$ ($1/n_{\text{on}}$ on drug, $-1/n_{\text{off}}$
off drug), define the contrast-weighted exposure

$$
F = \sum_t a_t f(t) = \bar f_{\text{on}} - \bar f_{\text{off}} .
$$

The drug-response part of the sum $\sum_t a_t BR_{it}$ is
$\omega_i F + \sum_t a_t\xi_{it}$; because $Y = \mathrm{BL} - (TV + PB
+ BR)$, it enters the contrast $\Delta_i = \sum_t a_t Y_{it}$ with a
negative sign. Writing $V_0$ for the variance contributed by the serial
noise and the other response factors, assumed independent of
$\omega_i$,

$$
\beta_{bm:D} = -\frac{\bar\omega\zeta\,F}{\sigma_{bm}}, \qquad
\mathrm{Var}(\Delta) = (\bar\omega^2\zeta^2 + \sigma_v^2)\,F^2 + V_0, \qquad
r_{\Delta B} = -\frac{\bar\omega\zeta\, F}
  {\sqrt{(\bar\omega^2\zeta^2 + \sigma_v^2) F^2 + V_0}} .
$$

(The first draft gave $r_{\Delta B}$ without the minus sign; the
sign follows from that of $\beta_{bm:D}$ and does not affect power,
which depends on $\lvert r_{\Delta B}\rvert$. Corrected 2026-10-07.)

**Proposition.** *$\lvert r_{\Delta B}\rvert$ is strictly increasing
in $\lvert F\rvert$.* With $G = \bar\omega\zeta > 0$ and
$H = \bar\omega^2\zeta^2 + \sigma_v^2$,
$\lvert r_{\Delta B}\rvert = G\lvert F\rvert/\sqrt{HF^2 + V_0}$, whose
derivative in $\lvert F\rvert$ is $GV_0/(HF^2 + V_0)^{3/2} > 0$ for
$V_0 > 0$. Since power is increasing in $\lvert r_{\Delta B}\rvert$
(`docs/37`, Section 4.4), power is strictly increasing in the magnitude
of the contrast-weighted exposure (derived).

The half-life enters through $F$. Under anchored offset with the
on-drug course held fixed, each off-drug $f(t)$ increases with
$t_{1/2}^{R}$, so $F$ decreases and power falls monotonically, as in
`docs/37`, Proposition 3. Under the full turnover model a longer
response half-life also slows onset and changes the restart level, so
$F$ is not guaranteed to be monotone in $t_{1/2}^{R}$ for every
schedule; it must be checked per design. The monotone relation
between power and $F$ holds regardless.

Two features differ from the published construct:

- **Variance depends on carryover.** In `docs/37`, Proposition 1, the
  contrast variance did not depend on the half-life. Here the term
  $(\bar\omega^2\zeta^2 + \sigma_v^2)F^2$ does, because sensitivity
  variation contributes noise in proportion to the exposure contrast.
- **Power is capped by what the biomarker explains.** The bound
  $\lvert r_{\Delta B}\rvert < \bar\omega\zeta/\sqrt{\bar\omega^2\zeta^2
  + \sigma_v^2}$ holds for every $F$ and is approached as $HF^2$
  dominates $V_0$; it is the square root of the share of sensitivity
  variation the biomarker accounts for. Since $f \in [0, 1]$ gives
  $\lvert F\rvert \le 1$, the attainable value is lower still. No
  design and no amount of exposure contrast can exceed the bound. This is a property of the biology the model
  encodes, not of the simulation's construction, and it gives the
  effect size a direct interpretation.

## 7. Comparison with the constructs in use

| Feature | `orig` (58b32a9) | `covar` (package) | Graded + separable (`docs/36` C) | Mean moderation (D) | Proposed |
|---|---|---|---|---|---|
| Off-drug mean | recursive decay | anchored decay | recursive (as published) | recursive (as published) | anchored, from turnover |
| Restart after off-drug interval | from zero | from zero | from zero | from zero | from residual level |
| Moderation on drug | full at every visit | full | full | full | proportional to response |
| Moderation off drug | full if any residual mean (step) | decays | decays | none | proportional to response |
| Drug-response variance off drug | full | full | full | full | shrinks with response |
| Valid for every schedule | no (ceiling, repair) | no (ceiling, Cholesky fallback) | response block yes; ceiling remains | yes | yes, by construction |
| Ceiling on effect size | yes | yes | yes | none | none (correlations are outputs) |

Table: Features of orig, covar, graded separable, mean moderation and proposed constructs

The mean and restart entries for `covar` were inspected in the package
source: `R/generateData.R` decays the off-drug mean from the last
on-drug value, and `R/buildtrialdesign.R` resets time on drug to zero
after an off-drug interval. The remaining `covar` entries summarize
`docs/32` and `docs/36`.

(Updated 2026-10-09: the "recursive (as published)" entries for
configurations C and D describe the drivers as they were when this
paper was written and the results stored from them. Since 2026-10-08
the quick-sim drivers (`hendrickson-problems/04-power-simulation.R`,
`carryover-closed-form/08-covariance-study.R` and
`11-strawman-stress-test.R`) default to anchored carryover, with
`--carryover recursive` reproducing the stored runs (`docs/53`;
`analysis/scripts/quick-sim/carryover-closed-form/RERUN-anchored.md`).
The power simulation reported in `docs/36`, the covariance study and
the stress test were produced with the recursion. An anchored rerun of
the `docs/36` power simulation with the published analysis is stored
(`power-summary-reps500-anchored.csv`, used in `docs/37`, Section 6);
anchored reruns of the covariance study and the stress test are
pending on the simulation server.)

## 8. Implications for analysis and design

- **Analysis.** The analysis should use an exposure variable built from
  the same time course: the continuous indicator $D_{bc}$ of the
  compendium, with the response half-life rather than the drug
  half-life, weights each visit by its expected exposure. A nonlinear
  mixed model of the turnover process, fitted directly, is the
  principled alternative and would estimate the response half-life.
- **Design.** Off-drug visits early after discontinuation carry the
  least contrast when the response half-life is long. Washout or
  off-drug assessment at four to five response half-lives restores
  $\bar f_{\text{off}} \approx 0$; if the response half-life is a week
  or more, an assessment one week after stopping contributes little to
  the interaction whatever the analysis.
- **Sensitivity.** Because the response half-life is uncertain, power
  should be reported across a range of plausible values rather than at
  a single assumed one.

## 9. What adopting the model would involve

- **Package.** A new data-generating option alongside
  `dgp_architecture = 'mvn'` and `'mean_moderation'` in
  `R/generateData.R`, with parameters for sensitivity
  ($\bar\omega$, matching the Gompertz maximum), moderation ($\zeta$
  or $c_{bm}^{\text{peak}}$), unexplained sensitivity variation
  ($\sigma_v$), serial noise ($\sigma_\xi$ and its correlation)
  and the response half-life. No correlation matrix is formed.
- **Checks.** With $\zeta = 0$ the interaction test must hold its size;
  with $t_{1/2}^{R} \to 0$ the results must reduce to the no-carryover
  case; the simulated moments must match Section 5.
- **Compendium.** Every paper that simulates through `generateData()`
  is affected in its data-generating assumptions, and results would
  need to be compared with the current architectures before any are
  replaced. The closed-form results of `docs/37` carry over with the
  coupling gap replaced by the contrast-weighted exposure $F$.

This is a substantial change to the data-generating process and is
offered as a proposal for decision, not as a completed revision.

## 10. Limitations and evidence status

- **Derived, not tested.** Sections 3 to 6 are derivations. No part of
  the proposed model has been implemented or simulated.
- **Pharmacological premises.** The prazosin plasma half-life is from
  one study in 24 normotensive volunteers (Hobbs et al., 1978). The claim that
  persistent benefit in PTSD is pharmacodynamic rather than
  pharmacokinetic follows from that half-life and from general
  pharmacodynamic principles; the response half-life itself is unknown
  to us.
- **Multiplicative moderation is a modeling choice.** It assumes the
  biomarker acts on sensitivity. A biomarker that acted on the onset or
  offset rate would require moderation of $k_{\text{out}}$ instead, and would give a
  different moderation profile.
- **Normality.** The derivations of Section 6 assume joint normality of
  $b_i$, $v_i$ and the noise; with a multiplicative random effect
  ($\omega_i f(t)$, with $f$ a fixed function of the schedule), the
  responses are normal for each schedule, but
  extensions with a nonlinear link would not be.
- **Simplifications.** The other response factors are assumed
  independent of $\omega_i$ in Section 6; shared patient-level effects
  would add covariance terms to $V_0$ and to $\beta_{bm:D}$.

## 11. References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020; 2:13. doi:10.3389/fdgth.2020.00013. Code:
   `github.com/rchendrickson/pmsimstats`,
   commit `58b32a9`.
2. Hobbs DC, Twomey TM, Palmer RF. Pharmacokinetics of prazosin in man.
   *Journal of Clinical Pharmacology* 1978; 18(8-9):402-406.
   doi:10.1002/j.1552-4604.1978.tb02456.x (PMID 690251).
3. Dayneka NL, Garg V, Jusko WJ. Comparison of four basic models of
   indirect pharmacodynamic responses. *Journal of Pharmacokinetics and
   Biopharmaceutics* 1993; 21(4):457-478. doi:10.1007/BF01061691
   (PMID 8133465).
4. Sheiner LB, Stanski DR, Vozeh S, Miller RD, Ham J. Simultaneous
   modeling of pharmacokinetics and pharmacodynamics: application to
   d-tubocurarine. *Clinical Pharmacology and Therapeutics* 1979;
   25(3):358-371. doi:10.1002/cpt1979253358 (PMID 761446).
5. pmsimstats team. `docs/35-separable-response-covariance.md`,
   `docs/36-hendrickson-pd-and-carryover.md`,
   `docs/37-carryover-power-closed-form.md`.

---

## Derivations

<!-- one row per result: statement | CM label | verification script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the numbers agree? -->

---

## Limitations and open items

<!-- link each to the open-items register at docs/series/registry/open-items.csv -->

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
