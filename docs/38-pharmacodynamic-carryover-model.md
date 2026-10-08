# Modeling Carryover from Pharmacodynamic Principles: A Generative Proposal for the Mean, Moderation and Covariance of the Drug Response {.unlisted .unnumbered}
*2026-10-02 17:31 PDT*

**Author.** pmsimstats team

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

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
   with a half-life of about 2.3 hours (Hobbs et al., 1978), so
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
last dose. For prazosin it is short: in 24 healthy subjects, drug left
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
\frac{df}{dt} = k\,\big(E(t) - f(t)\big), \qquad
k = \frac{\ln 2}{t_{1/2}^{R}},
$$

where $E(t)$ is the drug-driven target: the on-drug effect while on
drug, 0 off drug. For a target that is constant within each on-drug
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
BR_{it} = \theta_i\, f(t) + \varepsilon_{it}, \qquad
\theta_i = \bar\theta\,(1 + \kappa z_i) + u_i ,
$$

with $z_i$ the standardized biomarker, $\kappa$ the strength of
moderation, $u_i \sim N(0, \sigma_u^2)$ the variation in sensitivity
that the biomarker does not explain, and $\varepsilon_{it}$ serial
within-patient noise with variance $\sigma_\varepsilon^2$, independent
of $z_i$ and $u_i$. The mean response is $\bar\theta f(t)$, and the
biomarker-dependent part, $\bar\theta\kappa z_i f(t)$, has the same
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

The decayed mean moderation of the `docs/37` discussion, with an
anchored profile equal to 1 on drug and $2^{-t_{sd}/t_{1/2}}$ off drug,
is an approximation to this; the proportional profile $f(t)$ is the
form the pharmacology implies.

## 5. The implied covariance

Under the model of Section 4 every second moment of the drug response
follows from four quantities, $\bar\theta\kappa$, $\sigma_u$,
$\sigma_\varepsilon$ and the serial correlation of $\varepsilon$
(derived):

$$
\begin{aligned}
\mathrm{Var}(BR_t) &= (\bar\theta^2\kappa^2 + \sigma_u^2)\, f(t)^2 + \sigma_\varepsilon^2, \\
\mathrm{Cov}(BR_t, BR_s) &= (\bar\theta^2\kappa^2 + \sigma_u^2)\, f(t) f(s)
  + \sigma_\varepsilon^2\, r_\varepsilon(t, s), \\
\mathrm{Cov}(B, BR_t) &= \sigma_B\,\bar\theta\kappa\, f(t), \qquad
\mathrm{Cor}(B, BR_t) = \frac{\bar\theta\kappa\, f(t)}
  {\sqrt{(\bar\theta^2\kappa^2 + \sigma_u^2) f(t)^2 + \sigma_\varepsilon^2}} .
\end{aligned}
$$

Four consequences:

- **The biomarker correlation follows the response.** It is zero where
  $f = 0$, rises and decays with $f$, and is never a step. The
  quantity the published framework calls $c_{bm}$ becomes a derived
  value, most naturally the correlation at full effect,
  $c_{\text{peak}} = \bar\theta\kappa /
  \sqrt{\bar\theta^2\kappa^2 + \sigma_u^2 + \sigma_\varepsilon^2}$.
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
responsiveness entering both $\theta_i$ and the expectancy effect,
rather than from constants placed in matrix cells. The resulting
covariance is again a sum of positive semidefinite terms, in the form
of the linear model of coregionalization of `docs/35`.

Data from this model need no correlation matrix at all: draw $z_i$,
$u_i$ and the serial noise, compute $f(t)$ for the patient's schedule,
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

The drug-response part of each patient's contrast is
$\theta_i F + \sum_t a_t\varepsilon_{it}$. Writing $V_0$ for the
variance contributed by the serial noise and the other response
factors, assumed independent of $\theta_i$,

$$
\gamma = -\frac{\bar\theta\kappa\,F}{\sigma_B}, \qquad
\mathrm{Var}(\Delta) = (\bar\theta^2\kappa^2 + \sigma_u^2)\,F^2 + V_0, \qquad
\rho_{\Delta B} = \frac{\bar\theta\kappa\, F}
  {\sqrt{(\bar\theta^2\kappa^2 + \sigma_u^2) F^2 + V_0}} .
$$

**Proposition.** *$|\rho_{\Delta B}|$ is strictly increasing in $|F|$.*
With $A = \bar\theta\kappa$, $C = \bar\theta^2\kappa^2 + \sigma_u^2$,
$\rho = AF/\sqrt{CF^2 + V_0}$ has derivative
$AV_0/(CF^2 + V_0)^{3/2} > 0$. Since power is increasing in
$|\rho_{\Delta B}|$ (`docs/37`, Section 4.4), power is strictly
increasing in the contrast-weighted exposure.

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
  $(\bar\theta^2\kappa^2 + \sigma_u^2)F^2$ does, because sensitivity
  variation contributes noise in proportion to the exposure contrast.
- **Power is capped by what the biomarker explains.** As $F$ grows,
  $\rho_{\Delta B} \to \bar\theta\kappa/\sqrt{\bar\theta^2\kappa^2 +
  \sigma_u^2}$, the square root of the share of sensitivity variation
  the biomarker accounts for. No design and no amount of exposure
  contrast can exceed it. This is a property of the biology the model
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
  ($\bar\theta$, matching the Gompertz maximum), moderation ($\kappa$
  or $c_{\text{peak}}$), unexplained sensitivity variation
  ($\sigma_u$), serial noise ($\sigma_\varepsilon$ and its correlation)
  and the response half-life. No correlation matrix is formed.
- **Checks.** With $\kappa = 0$ the interaction test must hold its size;
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
  one study in healthy volunteers (Hobbs et al., 1978). The claim that
  persistent benefit in PTSD is pharmacodynamic rather than
  pharmacokinetic follows from that half-life and from general
  pharmacodynamic principles; the response half-life itself is unknown
  to us.
- **Multiplicative moderation is a modeling choice.** It assumes the
  biomarker acts on sensitivity. A biomarker that acted on the onset or
  offset rate would require moderation of $k$ instead, and would give a
  different moderation profile.
- **Normality.** The derivations of Section 6 assume joint normality of
  $z_i$, $u_i$ and the noise; with a multiplicative random effect
  ($\theta_i f(t)$), the responses are normal for each schedule, but
  extensions with a nonlinear link would not be.
- **Simplifications.** The other response factors are assumed
  independent of $\theta_i$ in Section 6; shared patient-level effects
  would add covariance terms to $V_0$ and to $\gamma$.

## 11. References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020; 2:13. Code: `github.com/rchendrickson/pmsimstats`,
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
