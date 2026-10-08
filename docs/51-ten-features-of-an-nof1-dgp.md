# Seventeen Features of a Data-Generating Process for Simulating Aggregated N-of-1 Trials: A Framework for Comparison {.unlisted .unnumbered}
*2026-10-07 15:56 PDT*

**Author.** pmsimstats team

**Revision 2026-10-07.** Expanded from ten to seventeen features.
Features 11 to 17 (Sections 13 to 19) were added after review. They
meet the second and third selection criteria of Section 2, but not
yet the first: their effects are derived or taken from the literature,
not simulated in this compendium.

**Purpose.** A Monte Carlo study of an N-of-1 design can only show
what its data-generating process (DGP) lets the data contain. This
paper identifies seventeen features of a DGP that determine what a
simulation of an aggregated N-of-1 trial can and cannot reveal. They
can be used to describe, compare and choose among DGPs.

For each feature we give the intuition, the mathematics, the relevant
literature, the choices made in this compendium and the evidence for
them, and a best-practice recommendation. Section 20 compares five
DGPs on all seventeen features:

- the published process of Hendrickson et al. (2020; "RH");
- the package's current construct;
- the refined process proposed in paper 01;
- mean moderation;
- the pharmacodynamic sensitivity model of `docs/38`.

Section 21 is a summary written for the introduction of paper 01.
Evidence labels follow the compendium convention: derived, verified
(computed or simulated), inspected, argued.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. Symbols marked *local*
are defined here.

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $Y_{it}$ | symptom score of participant $i$ at visit $t$; $t = 0$ is baseline | canonical |
| $\mathrm{BL}_i$ | baseline level | canonical |
| $BR_{it}$, $PB_{it}$, $TV_{it}$ | drug, expectancy (placebo) and natural-course responses | canonical |
| $B_i$, $b_i$ | biomarker; standardized biomarker | canonical |
| $t_{od}$, $t_{sd}$ | time on drug in the current on-drug period; time since discontinuation | canonical |
| $t_{1/2}$ | carryover (response) half-life | canonical |
| $e_t$ | expectancy: 1 open label, 0.5 blinded, 0 at baseline | canonical |
| $f(t)$ | population drug-response time course, scaled to 1 at full effect | local (`docs/38`) |
| $g_t$ | coupling profile: fraction of the biomarker's moderation present at visit $t$ | local (`docs/37`) |
| $c_{bm}$ | covariance-moderation strength (correlation) | canonical |
| $\beta_{bm}$ | mean-moderation strength (multiplier) | canonical |
| $\omega_i$, $\bar\omega$, $\zeta$, $u_i$ | participant drug sensitivity, its mean, the biomarker's effect on it, unexplained sensitivity | local (`docs/38`, renamed in `docs/44`) |
| $\rho$ | AR(1) serial correlation per week | canonical |
| $K$, $A$ | component correlation; temporal correlation kernel | local (`docs/35`) |
| $\beta_{bm:D}$ | biomarker-by-treatment interaction, the estimand | canonical |

Table: Notation index: symbols, meanings and status

**Glossary.**

- **DGP.** The probability model from which simulated trials are drawn.
- **Predictive / prognostic biomarker.** Associated with the size of the
  treatment effect, or with the outcome irrespective of treatment.
- **Coupling.** How strongly the biomarker is tied to the drug response
  at a given visit.
- **Ceiling.** The largest biomarker effect for which a covariance-based
  DGP has a valid joint distribution.

## 1. Summary

| # | Feature | Question it answers | Recommended choice |
|--|--------------|--------------------|-----------------------------------|
| 1 | Outcome decomposition | What is the score made of? | Baseline minus drug, expectancy and natural-course responses, each with its own clock |
| 2 | Onset of the drug response | How does benefit build on drug? | A saturating time course in time on drug, with clinically read landmarks; an effect-compartment or turnover model as the principled form |
| 3 | Offset and carryover | How does benefit decay off drug? | Pharmacodynamic half-life, anchored at discontinuation, zero before first exposure |
| 4 | Re-exposure | What happens when the drug restarts? | Drug-specific: restart from zero for re-titrated drugs, from the residual otherwise |
| 5 | Encoding of the predictive biomarker | How is the interaction placed in the model? | A moderation profile that fades with the drug effect, in the mean or the covariance; never a step |
| 6 | Prognostic pathways | Can the biomarker act without the drug? | Optional couplings to natural course, placebo response and baseline |
| 7 | Within-person dependence | How do deviations persist over visits and across components? | Separable AR(1) plus a person-level term, valid for every schedule |
| 8 | Between-person heterogeneity of response | Do participants differ in response beyond the biomarker? | Explicit unexplained sensitivity; drug-response variance that shrinks off drug |
| 9 | Observation process | How are the latent responses observed? | Phase-specific measurement error, the visit schedule as data, specified missingness mechanisms |
| 10 | Validity, estimand and calibration | Is the DGP coherent, and what does the simulation estimate? | Valid by construction (stop, never repair); an explicit estimand; parameters calibrated and reported |
| 11 | Biomarker measurement and distribution | Is the biomarker observed exactly, and how is it distributed? | A reliability parameter below 1; a realistic distribution where the encoding allows it |
| 12 | Functional form of moderation | Is the drug response linear in the biomarker? | Linear as the reference, with a threshold or responder-mixture scenario |
| 13 | Exposure versus assignment | Does the participant receive the drug as assigned? | Exposure separated from assignment: titration ramp and adherence below 1 as stress scenarios |
| 14 | Integrity of blinding | Does expectancy depend on the actual drug state? | An unblinding parameter, zero in the reference and positive in a stress scenario |
| 15 | Between-person variation in carryover | Do participants share one carryover half-life? | Log-normal person-level half-lives as a stress scenario |
| 16 | Selection and regression to the mean | How do enrollment criteria shape the sample? | A noisy screening score with an eligibility threshold, kept separate from the baseline visit |
| 17 | Outcome scale | Can simulated scores leave the instrument's range? | Check the range; use a bounded model when out-of-range scores are not rare |

Table: The seventeen DGP features, the question each answers, and the recommended choice

Features 1 to 10 are the ones on which existing DGPs differ. Features
11 to 17 are shared gaps: all five DGPs compared in Section 20 make
the same simplifying choice on each.

Three cross-cutting findings motivate the list:

1. **Defects in features 3 and 5 together produced the published power
   collapse under carryover** (`docs/37`).
2. **Feature 7 bounds what feature 5 can represent** (`docs/34` to
   `docs/36`).
3. **Feature 6 decides which analysis is valid** (`docs/48`, `docs/49`).

A DGP comparison that varies one feature while holding the others fixed
is the only way to attribute a result to a modeling choice.

## 2. Framing: what a DGP must do

A simulation study defines aims, data-generating mechanisms, estimands,
methods and performance measures (the ADEMP structure of Morris, White
and Crowther 2019 [1]; see also Burton et al. 2006 [2]). For N-of-1
designs, the DGP must do three things.

- **Represent the clinical phenomena the design is meant to resolve:**
  within-person contrasts over time, carryover, expectancy, and
  between-person differences in response.
- **Encode a well-defined estimand,** here the biomarker-by-treatment
  interaction $\beta_{bm:D}$.
- **Be valid and verifiable.** It must be a proper probability model,
  with every parameter interpretable and every result reproducible.

Published simulation studies of N-of-1 designs make these choices
implicitly and differently [31, 32, 33, 34], and the CENT reporting
statement asks for several of them to be reported for real trials
[35]. The first ten features below are the places where DGPs for
these designs actually differ. We selected them by three criteria:

- the choice changes at least one simulation result materially
  (verified in this compendium);
- it has a pharmacological or clinical interpretation;
- it can be stated and compared across DGPs.

Features 11 to 17 meet the second and third criteria. For the first,
the evidence is a derivation or the literature, not a simulation in
this compendium; each section says which.

## 3. Feature 1: outcome decomposition

**What it is.** The symptom score is built from separately generated
components, each with its own driver:

$$
Y_{it} = \mathrm{BL}_i - \big(TV_{it} + PB_{it} + BR_{it}\big).
$$

| Component | Driver |
|---|---|
| $TV$ (natural course) | time in the trial |
| $PB$ (expectancy) | time in the trial, scaled by the expectancy $e_t$ |
| $BR$ (drug response) | time on drug |

Table: Components of the symptom score and the clock that drives each

**Intuition.** An N-of-1 contrast compares a participant on drug with
the same participant off drug. What makes that contrast informative,
or misleading, are the processes that change over the same weeks: the
condition improves on its own, belief in treatment produces benefit,
and blinding changes that belief. A DGP that generates only "drug
effect plus noise" cannot show how a design copes with these.

**Mathematics.** Each component is observed only through the sum.
Identification therefore comes from design variation:

- expectancy varies between the open-label and blinded phases;
- drug state varies within the blinded phase;
- time varies throughout.

In the Hybrid design, the open-label visits are all on drug at
$e_t = 1$ and every off-drug visit is blinded at $e_t = 0.5$. Expectancy
and drug are therefore partly confounded, and only the blinded visits
separate them cleanly (`docs/40`, Section 7; `docs/49`, Section 2).

**Literature.** Decomposing response into drug action and disease
progression is standard in pharmacometrics [3, 4], and placebo response
has its own longitudinal models [28]. That expectancy adds benefit when
a treatment is known to be given is shown by the open-versus-hidden
paradigm [29]. Its dependence on trial design (open label against
placebo-controlled) is modeled by [30]. The three-component
construction for N-of-1 designs is RH's [5].

**This compendium.** All DGPs use RH's decomposition and Gompertz
parameters. The component decomposition paper (06) and `docs/31`
examine whether the components can be recovered.

**Best practice.**

- Generate each component on its own clock.
- State which design contrasts identify which component.
- Report component-level parameters on the outcome scale.

## 4. Feature 2: onset of the drug response

**What it is.** The time course of benefit after the drug is started.

**Intuition.** Benefit builds over weeks, through titration and the
accumulation of effect. Early on-drug visits therefore show a partial
effect, and a design whose on-drug periods are short sees mostly
partial effects.

**Mathematics.** RH use a modified Gompertz curve in time on drug,

$$
G(t_{od}; m, d, r) = m\,\frac{\exp(-d e^{-r t_{od}}) - e^{-d}}{1 - e^{-d}},
$$

with $m = 10.99$, $d = 5$ and $r = 0.42$ per week. Its landmarks:

- the inflection is at $t = \ln d / r = 3.8$ weeks;
- half of the maximum is reached at 4.7 weeks;
- 90% is reached at 9.2 weeks (`docs/40`, Section 4).

The pharmacometric alternative is a turnover (effect-compartment) model,

$$
\frac{df}{dt} = k\big(E(t) - f(t)\big), \qquad k = \frac{\ln 2}{t_{1/2}^{R}},
$$

in which onset and offset share one rate. On drug from level $f_0$,
$f(t) = 1 - (1 - f_0)\,2^{-t/t_{1/2}^{R}}$, so onset is exponential
rather than sigmoid. A sigmoid onset follows if the target $E(t)$
itself rises with titration.

**Literature.**

- **Effect-compartment models** describe the lag between concentration
  and effect [6].
- **Indirect-response (turnover) models** describe responses that
  build and decay through production and loss [7].
- **Emax and sigmoid-Emax models** relate exposure to effect [8].

**This compendium.** RH's Gompertz onset is kept in every current DGP.
`docs/38` proposes the turnover form, which is not yet implemented.

**Best practice.**

- **Read the onset clinically** before simulating: how long until half
  of the effect is reached. Check it against the design's on-drug
  period lengths.
- **Prefer an onset that is coherent with the offset** (Feature 3).

## 5. Feature 3: offset and carryover

**What it is.** How the drug response decays after the drug is stopped.

**Intuition.** Two kinds of carryover exist, and they act on very
different time scales.

- **Pharmacokinetic carryover** is drug still in the body. For prazosin
  it is about 2.3 hours of plasma half-life [9], gone within a day.
- **Pharmacodynamic carryover** is benefit that persists after the drug
  has gone, set by how fast the response turns over: days to weeks.

At weekly visits, only the pharmacodynamic kind matters, so the
simulation's half-life is a response half-life.

**Mathematics.** Anchored decay from the level at discontinuation:

$$
\mu_{BR,t} = \mu_{BR,\text{stop}}\; 2^{-t_{sd,t}/t_{1/2}} .
$$

RH's recursion instead decays the previous visit's already-decayed mean
by the total time since discontinuation. At off-drug visits one and two
weeks after stopping, the second is then at $2^{-(1+2)/t_{1/2}}$ instead
of $2^{-2/t_{1/2}}$, so the decay compounds (derived; `docs/37`,
Appendix A). Before first exposure $t_{sd} = 0$, but there is no drug
effect, so the response is zero. A decay function evaluated at
$t_{sd} = 0$ returns 1 and must not be used there.

That error, made in the analysis-side indicator of manuscript 02, cost
the exposure-weighted analysis 0.26 to 0.30 in power in the crossover
design (verified, `docs/45`, Section 8).

**Literature.**

- **Carryover is the classical problem of crossover designs** [10, 11].
- **Behavioral carryover** that outlasts the pharmacological effect has
  been modeled explicitly [12].
- **Washout and analytical treatments** in N-of-1 designs [13].
- **Causal estimation under carryover** in aggregated N-of-1 trials,
  with a simulated DGP [27].

**This compendium.**

- **Decay forms:** exponential, with Weibull as a sensitivity option
  (manuscript 02, the stress test).
- **Half-life grid:** 0, 0.5 and 1 week.
- **Pharmacological range:** 0.1 to 0.2 week for prazosin (paper 01).

**Best practice.**

- **Specify carryover as a pharmacodynamic half-life,** with its
  justification.
- **Anchor the decay at discontinuation.**
- **Code pre-exposure visits as unexposed.**
- **Vary the decay form as a sensitivity analysis.**

## 6. Feature 4: re-exposure

**What it is.** The drug response when the drug is restarted after an
off-drug period.

**Intuition.** Either the response rebuilds from zero, which is right
if the drug is re-titrated, or it continues from whatever remains of
the earlier effect, which is right if the residual carries forward.

**Mathematics.**

- **Restart from zero:** $\mu_{BR,t} = G(t_{od})$ with $t_{od}$ reset.
  In the Hybrid crossover, a participant at near-full response at
  week 10 is credited with $G(4) = 4.28$ of 10.99 at week 16, 39% of
  the maximum.
- **Turnover model:** restart from the residual is automatic,
  $f(t) = 1 - (1 - f_{\text{res}})\,2^{-t/t_{1/2}^{R}}$.
- **Combined rules:** anchored offset with a zero restart can step the
  mean down at re-exposure, whenever the residual exceeds the restarted
  curve.

**Literature.** Re-titration after an interruption is standard
prescribing practice for prazosin, because of first-dose hypotension.
Turnover models give restart from the residual by construction [7].

**This compendium.** RH, the package and the refined process all
restart from zero. Restart from the residual is roadmap step 7 of
`docs/44` and is not implemented.

**Best practice.**

- **Make the rule an explicit, drug-specific choice.**
- **For re-titrated drugs, restart from zero.**
- **Otherwise, restart from the residual,** preferably through a
  turnover model.

## 7. Feature 5: encoding of the predictive biomarker

**What it is.** How the biomarker-by-treatment interaction enters the
model. There are three forms:

- **Covariance binding (RH):** a correlation $c_{bm}g_t$ between the
  biomarker and $BR_{it}$.
- **Mean binding:** a shift of $\beta_{bm} b_i \sigma_{BR}\, g_t$ in
  $BR_{it}$.
- **Sensitivity binding** (`docs/38`): $BR_{it} = \omega_i f(t) +
  \varepsilon_{it}$ with $\omega_i = \bar\omega(1 + \zeta b_i) + u_i$.

**Intuition.** A predictive biomarker tells you how much drug effect to
expect. Wherever the drug effect is absent, the biomarker should
therefore predict nothing. What matters most is not *where* the
interaction is placed but its profile over visits, $g_t$, and above all
whether that profile fades with the drug effect.

**Mathematics.**

**(a) Conditional equivalence** (derived). Under joint normality,

$$
E(BR_{it} \mid B_i) = \mu_{BR,t} + c_{bm} g_t \sigma_{BR} b_i,
\qquad
\mathrm{Cov}(BR_{i\cdot} \mid B_i) = \Sigma_{BR} - c_{bm}^2\sigma_{BR}^2\, g g^\top .
$$

So covariance binding is mean binding with exposure profile $g_t$ and
$\beta_{bm} = c_{bm}$, and a rank-one smaller residual covariance.

**(b) Power depends on the profile through the coupling gap** (derived,
`docs/37`). For the paired-difference statistic,

$$
\beta_{bm:D} = -c_{bm}\frac{\sigma_{BR}}{\sigma_{bm}}\big(1 - \bar g_{\text{off}}\big),
$$

where $\bar g_{\text{off}}$ is the mean off-drug coupling. Power falls
with the gap $1 - \bar g_{\text{off}}$, and the required sample size
scales roughly as $1/\text{gap}^2$.

| Profile | Off-drug coupling | Gap under carryover |
|---|---|---|
| RH's step rule | 1 wherever any drug effect remains | 0: the interaction vanishes at any positive half-life |
| Graded coupling | $2^{-t_{sd}/t_{1/2}}$ | falls smoothly with the half-life |
| Sensitivity model | $f(t)$, proportional to the response | follows the response |

Table: Coupling profiles, their off-drug coupling, and the coupling gap under carryover

On RH's own code, graded coupling took power from 0.13 back to 0.64
at a half-life of 0.1 week (verified, `docs/37`, Appendix A).

**(c) The ceiling** (derived; `docs/34`, `docs/35`). Covariance binding
requires

$$
c_{bm}^2\, u^\top M^{-1} u < 1,
$$

by the Schur complement, where $M$ is the response correlation block
and $u$ the coupling pattern. Under RH's compound symmetry the ceiling
is 0.256, below RH's effect sizes of 0.3 and 0.6. RH's code repaired
the infeasible matrices silently: 20 of 54 at those effect sizes. Mean
and sensitivity binding have no ceiling: under sensitivity binding the
correlation is an output,

$$
\mathrm{Cor}(B, BR_t) = \frac{\bar\omega\zeta f(t)}
  {\sqrt{(\bar\omega^2\zeta^2 + \sigma_u^2) f(t)^2 + \sigma_\varepsilon^2}} .
$$

**Literature.**

- **Mean moderation** is the standard representation of
  treatment-effect heterogeneity by a covariate [14, 15].
- **Random participant-by-treatment variation** is the standard
  description of individual response in replicated crossovers
  [16, 17].
- **RH were, to our knowledge, the first** to parameterize the
  interaction as a covariance in the N-of-1 literature [5].
- **Report 03** shows that covariance binding matches the second
  moments of a latent responder-class mixture, with
  $c_{bm}^2 = f_B f_{BR}$, the product of the between-class variance
  fractions.

**This compendium.**

- **The refined process** uses graded coupling.
- **Paper 01** compares the three encodings.
- **The sensitivity model** is proposed and not yet implemented.

**Best practice.**

- **Choose the profile first, and make it fade with the drug effect.**
- **Choose the placement second.** Covariance binding is acceptable
  below its ceiling. Use mean or sensitivity binding when larger
  effects are needed.
- **Report the strength on two scales:** as the correlation or
  multiplier, and as the implied slope in outcome points per SD of
  biomarker.

## 8. Feature 6: prognostic pathways

**What it is.** Whether the biomarker can be associated with the
outcome through routes that do not involve the drug.

**Intuition.** Blood pressure might predict who improves on placebo, or
who improves anyway. A DGP in which the biomarker touches only the drug
response cannot reveal analyses that mistake such an association for a
drug interaction.

**Mathematics.** Each pathway implies a slope of $Y$ on $b$ at baseline
($s_0$) and afterwards ($s_t$):

| Pathway | $s_0$ | $s_t$ |
|---|---|---|
| natural course, constant | 0 | $-c_{bm,TV}\sigma_{TV}$ |
| natural course, growing | 0 | $-\beta^{TV}_{bm}\sigma_{TV}\,G_{TV}(t)/\max G_{TV}$ |
| placebo response | 0 | $-c_{bm,PB}\, 10\, e_t$ |
| baseline severity | $c_{bm,BL}\sigma_{BL}$ | $c_{bm,BL}\sigma_{BL}$ |

Table: Baseline and post-baseline biomarker slopes implied by each prognostic pathway

An analysis that shares one biomarker slope between baseline and the
off-drug visits fits $\hat\beta_B \approx f\,s_{\text{prog}}$, for some
$0 < f < 1$ set by the weight of the baseline row. It then reports a
spurious interaction of about $(1 - f)\,s_{\text{prog}}$ (heuristic,
`docs/48`, Section 2).

At a prognostic coupling of 0.3, the null rejection rate was 0.10 for
a random-intercept-plus-AR(1) analysis and 0.42 for an unstructured
covariance. Adding the expectancy terms restored it to 0.040 to 0.047
(verified, `docs/49`).

Under separability the joint validity condition generalizes to
$c_{bl}^2 + \sum_{cc'}[K^{-1}]_{cc'} v_c^\top A^{-1} v_{c'} < 1$
(derived and verified, `docs/35`, Section 4.8).

**Literature.**

- **The predictive/prognostic distinction** is central to the PATH
  statement [14, 15].
- **Treating baseline as a response** imposes a constancy constraint
  [18].

**This compendium.** Prognostic couplings to the natural course and to
the placebo response exist in the stress-test driver. A
baseline-severity coupling is being run. Only the placebo coupling is
in the package (`c.bm.pb`).

**Best practice.** Include prognostic pathways as standard stress
scenarios whenever analyses are compared. A DGP without them can
confirm an analysis but never refute it.

## 9. Feature 7: within-person dependence and cross-component structure

**What it is.** How a participant's deviations persist over visits,
and how the components co-vary.

**Intuition.** The structure can carry permanent traits (compound
symmetry), fading memory (AR(1)), or both. Separability means one
clock for all three components: the link between components fades at
the same rate as each component's own memory (`docs/35`, "In brief").

**Mathematics.**

- **Separable form.** $M = K \otimes A$, with $K$ the component
  correlation and $A_{ts} = \rho^{\lvert w_t - w_s\rvert}$. Its
  eigenvalues are products of those of $K$ and $A$, so it is valid for
  every visit schedule.
- **The non-separable alternative** fails whenever
  $\lambda_{\min}(A) \le (c_1 - c_\times)/(1 - c_\times)$, which
  happens on dense schedules.
- **RH's construct** is a sum of a person-level and an occasion-level
  separable term, each valid.
- **The ceiling factors.** Its schedule term is a sum of per-switch
  costs, for example $1/(1 - \phi^2)$ for off to on, with
  $\phi = \rho^{\text{gap}}$.
- **Implied outcome covariance.** With the person-level baseline,
  $$
  \mathrm{Cov}(Y_t, Y_s) = \sigma_{BL}^2 + (s^\top K s)\,\rho^{\lvert t-s\rvert},
  $$
  a random intercept plus AR(1). That is the structure assumed by the
  random-intercept-plus-`corCAR1` analysis.

**Consequences** (verified).

- **Power.** At the same $c_{bm}$, AR(1) data give much less power
  than compound symmetry: 0.335 against 0.721 in Hybrid at 0.25.
- **Analysis validity.** A compound-symmetry working model applied to
  AR(1) data rejected a true null 0.15 to 0.16 of the time in the
  crossover design (`docs/45`).

**Literature.**

- **Random effects plus serial correlation plus measurement error** as
  the general longitudinal model [19, 20].
- **AR(1) with unequally spaced visits** [21].
- **Kronecker (separable) covariances** for multivariate repeated
  measures [22].

**This compendium.**

- **RH:** compound symmetry, $\rho = 0.8$.
- **The package:** non-separable AR(1).
- **The refined process:** separable AR(1) at $\rho = 0.7$.

**Best practice.**

- **Specify the structure generatively** (person-level term plus
  serial process), so that validity holds by construction.
- **State $\rho$ and its basis.**
- **Compare effect sizes across structures on the slope scale.**
- **Check the analysis model against the structure the DGP implies.**

## 10. Feature 8: between-person heterogeneity of response

**What it is.** How much participants differ in their drug response
beyond what the biomarker explains, and how the drug-response variance
behaves over time.

**Intuition.** Two participants with the same blood pressure need not
respond equally. Under RH's construct, the drug-response component has
its full SD of 8 even off drug, where its mean is zero. That is as if
participants had a random drug response with no drug present, which
adds noise to every on-off contrast.

**Mathematics.** Under the sensitivity model,

$$
\mathrm{Var}(BR_t) = (\bar\omega^2\zeta^2 + \sigma_u^2)\, f(t)^2 + \sigma_\varepsilon^2,
$$

which shrinks to the serial-noise variance as the drug effect washes
out. The participant-by-treatment variance $\sigma_u^2$ is the quantity
that replicated designs estimate [16, 17]. It bounds the attainable
biomarker correlation by
$\bar\omega\zeta / \sqrt{\bar\omega^2\zeta^2 + \sigma_u^2}$ (derived,
`docs/38`, Section 5). A latent-class (bimodal) response is a different
heterogeneity, with nonlinear conditional mean (report 03).

**Literature.** Individual response and its variance components in
N-of-1 and replicated crossover designs [16, 17, 23].

**This compendium.** RH, the package and the refined process all use a
constant drug-response variance. The sensitivity model makes it
proportional to the response. Report 03 studies latent classes.

**Best practice.**

- **Separate explained from unexplained heterogeneity,** that is, $\zeta$
  from $\sigma_u$.
- **Let the drug-response variance follow the drug effect.**
- **Treat the distribution of response as a stated assumption,** normal
  or mixture.

## 11. Feature 9: the observation process

**What it is.** How the latent responses become recorded data:
measurement error, its dependence on the study phase, the visit
schedule, the scale, and missingness.

**Intuition.**

- **Ratings can be noisier in some phases.** Open-label and blinded
  phases can differ in rating noise, as well as in mean.
- **Visits are not equally spaced.** For example, the Hybrid design
  measures at weeks 4, 8, 9, 10, 11, 12, 16 and 20.
- **Participants drop out,** sometimes because they are doing badly.

**Mathematics.**

- **Phase-specific measurement error** adds a variance
  $\sigma^2_{\varepsilon,\text{phase}}$ to each visit, which a
  homoscedastic working model omits.
- **Missingness mechanisms.** Under MCAR or MAR, a likelihood analysis
  with all observed data is valid; under MNAR it is not [24, 25].
- **A MAR hazard** depending on the previous score,
  $\text{logit}\,h_t = a + \gamma\, z(Y_{t-1})$, can be calibrated to a
  target dropout rate per simulated trial (`docs/46`).

**Literature.** Missing-data theory [24, 25]. Informative dropout in
N-of-1 designs is the subject of paper 09.

**This compendium.** RH's censoring patterns; the stress test's phase
measurement error and its MCAR and MAR dropout; paper 09's informative
dropout. Under 20% MCAR or MAR dropout, the random-intercept-plus-AR(1)
analyses with CR2 standard errors held their size (at most 0.068). The
published analysis rejected 0.079 under MCAR dropout and 0.071 with
phase-specific measurement error. The unstructured model with AIC
selection rejected 0.075 under MAR (verified, stress test).

**Best practice.**

- **Treat the schedule as data.** Use calendar time in every
  correlation.
- **Include phase-specific error and at least MCAR and MAR dropout** as
  stress scenarios.
- **State the mechanism of any missingness simulated.**

## 12. Feature 10: validity, estimand and calibration

**What it is.** The properties that make a DGP trustworthy as the
truth against which analyses are judged.

**Intuition.** A simulation tests analyses against a stated truth. If
the truth silently changes (repaired matrices), is undefined (no
estimand), or is unrealistic (uncalibrated parameters), the results
cannot be interpreted.

**Mathematics and practice.**

- **Validity.** Either the DGP is generative (valid by construction), or
  every matrix is checked and infeasible cells stop the run. RH's
  silent repair changed correlations by up to 0.10 at $c_{bm} = 0.6$
  (verified, `docs/36`).
- **Estimand.** State the target, here the on-drug moderation slope,
  and whether an analysis with a binary drug indicator estimates it or
  an attenuated contrast.
- **Calibration.** Choose parameters from data where possible, and
  report effect sizes on the clinical scale. For prazosin, published
  trials give the plasma half-life [9] and a biomarker association on
  the Clinician-Administered PTSD Scale [26]. The simulated scale
  differs, so the mapping must be explicit.
- **Verifiability:**
  - a faithful reference arm reproducing the published process;
  - one change at a time;
  - common random numbers across arms;
  - a Monte Carlo standard error for every reported rate;
  - pre-registration of comparisons [1].

**This compendium.**

- **Reference arm:** a reproduction of RH's code at commit 58b32a9.
- **Paired designs** throughout.
- **Pre-registered stress test:** docs/46, with one amendment
  disclosed.
- **Pending:** calibration of effect sizes to trial data.

**Best practice.**

- **Stop, never repair.**
- **Define the estimand before the DGP.**
- **Keep a reference arm.**
- **Report ADEMP** [1].

## 13. Feature 11: biomarker measurement and distribution

**What it is.** Whether the biomarker that enters the analysis is the
biomarker that moderates the response, and what distribution it has.

**Intuition.** A blood pressure reading, a genotype score or an imaging
measure is observed with error. The drug response follows the true
value; the analysis sees the measured one. A noisy biomarker spreads
participants who are alike in truth across the measured range, which
flattens the estimated interaction.

**Mathematics.**

- **Attenuation.** Let the observed biomarker be $B^*_i = B_i + \delta_i$
  with independent error and reliability
  $\lambda = \operatorname{Var}(B)/\operatorname{Var}(B^*)$. For a normal
  biomarker, $E[B_i \mid B^*_i]$ is linear in $B^*_i$ with slope
  $\lambda$, so the interaction coefficient on $B^*$ is attenuated by
  $\lambda$ [36].
- **Power.** With the biomarker standardized before analysis, the
  coefficient is $\beta_{bm:D}\sqrt{\lambda}$ and its standard error is
  essentially unchanged when the interaction explains a small share of
  the residual variance. The noncentrality therefore falls by
  $\sqrt{\lambda}$, and the sample size needed for a given power rises
  by about $1/\lambda$ (derived). At $\lambda = 0.8$ a trial needs about
  25% more participants.
- **Covariance encoding.** Under covariance binding the correlation the
  analysis can see is $c_{bm}\sqrt{\lambda}$.
- **Distribution.** Covariance binding draws $B_i$ inside the
  multivariate normal, so the biomarker is necessarily normal. A binary
  or skewed biomarker (a genotype, a thresholded laboratory value) can
  be represented only under mean binding. Precision for the interaction
  is proportional to the biomarker's variance in the enrolled sample
  (see also Feature 16).

**Literature.** Random measurement error and regression dilution [36].
The prazosin motivating biomarker is pretreatment blood pressure [26],
which is measured with error.

**This compendium.** Every DGP measures the biomarker once, at
baseline, without error, from a normal distribution. The published
power figures are therefore upper bounds with respect to biomarker
reliability (argued).

**Best practice.**

- **Add a reliability parameter** $\lambda$ and report power at
  plausible values, not only at $\lambda = 1$.
- **Draw the biomarker from a realistic distribution** when the
  encoding allows it.
- **State whether power is quoted for the true or the measured
  biomarker.**

## 14. Feature 12: functional form of moderation

**What it is.** The shape of the relation between the biomarker and the
size of the drug response.

**Intuition.** A linear interaction says each unit of biomarker adds the
same benefit. Biology often suggests otherwise: a drug may work only
above a receptor or blood-pressure threshold, or participants may be
responders or non-responders, with the biomarker shifting the odds of
being a responder. When the effect changes sign across the biomarker
range, the interaction is qualitative [37].

**Mathematics.** Write the moderation as a function $m(b_i)$, so the
on-drug response is scaled by $1 + m(b_i)$.

- **Linear:** $m(b) = \beta b$, the form used throughout.
- **Threshold:** $m(b) = \Delta\,\mathbb{1}(b > \tau)$.
- **Responder mixture:** $\omega_i = \omega_R R_i$ with
  $R_i \sim \text{Bernoulli}\{\pi(b_i)\}$ and
  $\operatorname{logit}\pi(b) = a + \gamma b$.

A linear working model estimates the linear projection
$\beta_{\text{lin}} = \operatorname{Cov}\{m(b), b\}/\operatorname{Var}(b)$.
For a threshold at the median of a standard normal biomarker,
$\beta_{\text{lin}} = \Delta\,\phi(0) \approx 0.40\,\Delta$, and the
linear term captures $\phi(0)^2/(1/4) \approx 0.64$ of the moderation
variance. The linear test then has about 80% of the noncentrality of a
test that uses the true form (derived).

**Literature.** Qualitative and non-crossover interactions [37]; the
PATH statement on modeling heterogeneous treatment effects [14].

**This compendium.** Every DGP compared here is linear in the biomarker,
including the sensitivity model, where $\omega_i$ is linear in $b_i$.
Paper 03 (latent class mixture application) works with responder
classes in the analysis, not in the DGP.

**Best practice.**

- **Keep linear moderation as the reference.**
- **Add one threshold and one responder-mixture scenario** to show what
  the linear working model loses.
- **Report the moderation on the clinical scale** at several biomarker
  values, not only as a coefficient.

## 15. Feature 13: exposure versus assignment

**What it is.** The difference between the drug a participant is
assigned and the drug they receive: dose titration and adherence.

**Intuition.** Prazosin is not given at full dose from the first day.
In the largest veterans trial it was escalated in divided doses over
five weeks to a daily maximum of 20 mg in men and 12 mg in women [38].
Participants also miss doses. In both cases the response is driven by
exposure, while the design and the analysis clocks ($t_{od}$, $t_{sd}$)
run on assignment.

**Mathematics.**

- **Exposure.** Let $E_{it} = a_{it}\, d(t_{od})$, with adherence
  $a_{it} \in [0, 1]$ and a titration ramp
  $d(t_{od}) = \min(1, t_{od}/T)$, for example $T = 5$ weeks [38]. The
  drug response is driven by $E_{it}$, for example as the input of the
  turnover model of `docs/38`.
- **Attenuation.** If the response is linear in exposure and adherence
  is independent of the biomarker, the assigned-drug contrast and the
  interaction are both scaled by the mean adherence $\bar a$. Power
  falls accordingly, and the sample size rises by about $1/\bar a^{2}$
  (derived). Variable adherence also adds between-person variance
  (Feature 8).
- **Dependence on the biomarker.** If adherence depends on the
  biomarker, for example through side effects, the interaction is
  biased, in a direction set by that dependence (argued).

**Literature.** Undetected nonadherence in randomized trials leads to
underestimated efficacy and inflated outcome variability [39].
Titration schedule for prazosin in PTSD [38].

**This compendium.** All DGPs assume full adherence. Titration is not
modeled separately; the Gompertz onset of Feature 2 may absorb part of
it during the first on-drug period. Feature 4 (restart from zero on
re-exposure) is the drug-specific answer for re-titration.

**Best practice.**

- **Separate exposure from assignment** in the DGP.
- **Simulate adherence below 1** as a stress scenario, and one scenario
  in which adherence depends on the biomarker.
- **Model titration explicitly** when its duration is comparable with a
  treatment period.

## 16. Feature 14: integrity of blinding

**What it is.** Whether participants and raters can tell the drug state
during the blinded phases, so that expectancy follows the actual drug.

**Intuition.** The expectancy $e_t$ is fixed by design: 1 in open label
and 0.5 throughout the blinded phases. If a drug has a perceptible
effect, blinded participants may guess when they are on it. Expectancy
then rises on drug and falls off drug, and part of the placebo response
is mistaken for drug response. Prazosin lowered supine systolic blood
pressure by 6.7 mm Hg relative to placebo at 10 weeks [38]; whether
such an effect unblinds participants is unknown (argued).

**Mathematics.**

- **Unblinding parameter.** In blinded phases let
  $e_{it} = \tfrac12 + \eta\,(D_{it} - \tfrac12)$, with $\eta = 0$ for
  intact blinding and $\eta = 1$ for full unblinding.
- **Bias in the drug effect.** With $PB_{it} = e_{it} P_i(t)$, the
  blinded on-off contrast picks up $\eta P_i(t)$, which biases the drug
  main effect by $\eta \bar P$.
- **Bias in the interaction.** If the placebo response depends on the
  biomarker, $P_i = \bar p + s_{pb}\, b_i$ (a prognostic placebo
  pathway, Feature 6), the interaction picks up $\eta\, s_{pb}$
  (derived). Unblinding thus turns a prognostic pathway into a spurious
  predictive one. An expectancy term coded from the design (`docs/49`)
  is constant within the blinded phase and cannot absorb it.

**Literature.** The blinding index and the assessment of blinding
success [40].

**This compendium.** All DGPs fix $e_t$ by design; $\eta = 0$
throughout.

**Best practice.**

- **Add an unblinding parameter** $\eta$, zero in the reference
  scenario, and include $\eta > 0$ combined with a prognostic placebo
  pathway as a stress scenario.
- **Recommend that real trials report a blinding index** [40].

## 17. Feature 15: between-person variation in carryover

**What it is.** Whether participants share one carryover half-life.

**Intuition.** Pharmacokinetic and pharmacodynamic parameters vary
between people, and clinical trial simulation in pharmacometrics
usually represents that variation explicitly [41]. If carryover lasts
longer in some participants than in others, a single half-life in the
DGP understates how heterogeneous the off-drug visits are, and a single
half-life in the analysis is partly wrong for everyone.

**Mathematics.**

- **Person-level half-lives.** Let
  $t_{1/2,i} = t_{1/2} \exp(\sigma_h z_i)$ with $z_i \sim N(0, 1)$, and
  coupling $g_{it} = 2^{-t_{sd}/t_{1/2,i}}$ off drug.
- **Mean profile.** The population-mean profile
  $\bar g(t_{sd}) = E\{2^{-t_{sd}/t_{1/2,i}}\}$ is a mixture of
  exponentials. Its decay rate falls with time since discontinuation,
  so at long $t_{sd}$ it is dominated by the longest half-lives and
  decays more slowly than the profile at the median half-life (derived).
- **Consequences.** An analysis with a common half-life in $D_{bc}$
  fits the mean profile poorly at both short and long $t_{sd}$. The
  person-level coupling gaps also vary, which adds variance. To first
  order, power under a binary drug indicator still follows the mean
  coupling gap of `docs/37` (argued).

**Literature.** Between-subject variability in clinical trial
simulation [41]; prazosin's plasma half-life [9]. We found no estimate
of between-person variation in the pharmacodynamic carryover of
prazosin in PTSD.

**This compendium.** All DGPs, including the sensitivity model, use one
half-life (or one turnover rate) for all participants.

**Best practice.**

- **Simulate log-normal person-level half-lives** as a stress scenario.
- **Check that the AIC-selected $D_{bc}$** remains well sized when the
  half-life varies between participants.

## 18. Feature 16: selection and regression to the mean

**What it is.** The enrollment mechanism: who enters the trial, and on
what measured value.

**Intuition.** PTSD trials enroll people whose screening score exceeds
a threshold. Part of a high screening score is transient, so scores
fall afterward whatever the treatment. This regression to the mean can
look like natural course or placebo response. Selection also narrows
the distribution of anything correlated with the score.

**Mathematics.**

- **Regression to the mean.** Let the screening score be
  $Y_{i,\text{scr}} = \mathrm{BL}_i + \varepsilon_i$ with between- and
  within-person variances $\sigma_b^2$ and $\sigma_w^2$, and enroll when
  $Y_{i,\text{scr}} > \tau$. The expected fall after selection is
  $\sigma_w^2 / \sqrt{\sigma_b^2 + \sigma_w^2}\; C(z)$, with
  $C(z) = \phi(z)/\{1 - \Phi(z)\}$ and $z$ the standardized threshold
  [42].
- **Range restriction.** If the biomarker is prognostic for baseline
  severity (Feature 6), selection on the screening score reduces its
  variance in the enrolled sample to
  $1 - r^2 C(z)\{C(z) - z\}$, where $r$ is the biomarker's correlation
  with the screening score. Interaction precision is proportional to
  that variance, so power falls (derived).
- **Within-person contrasts.** The blinded on-off contrast is largely
  protected, because regression to the mean affects all post-screening
  visits alike. Contrasts that lean on the baseline visit, including
  analyses that share the biomarker effect with baseline (`docs/48`),
  are not.

**Literature.** Regression to the mean and how to address it in design
and analysis [42].

**This compendium.** All DGPs draw the baseline from an unrestricted
distribution; there is no screening step. The stress-test cells with a
prognostic baseline (cells 40 and 41) have no selection.

**Best practice.**

- **Simulate a noisy screening score** with an eligibility threshold.
- **Keep the screening value separate from the baseline visit** of the
  analysis.

## 19. Feature 17: outcome scale

**What it is.** The range and distribution of the recorded score.

**Intuition.** Symptom scales are bounded sums of item ratings. A
Gaussian DGP can produce scores below the floor, particularly for
strong responders late in an on-drug period. Under positive moderation
those are the high-biomarker participants who carry the interaction.

**Mathematics.**

- **Bounded models.** Either generate on a latent scale and map to the
  instrument, $Y = L + (U - L)\operatorname{logit}^{-1}(\eta)$, or use a
  beta model, which keeps simulated scores within the instrument's range
  even when residuals are added [43].
- **Censoring at the bounds** (clipping) is simple but distorts the mean
  and variance near the floor.
- **Consequence.** Near the floor, responses are compressed. The
  interaction is attenuated and the residual variance becomes
  heteroscedastic, falling where it matters most (argued).

**Literature.** Beta regression for simulating a bounded clinical scale
[43].

**This compendium.** All DGPs generate an unbounded Gaussian score. The
fraction of simulated scores outside the instrument's range has not
been computed.

**Best practice.**

- **Report the fraction of simulated scores outside the range.**
- **Use a bounded model** when that fraction is not negligible.

## 20. The five DGPs compared

| Feature | RH (58b32a9) | Package `covar` | Refined (paper 01) | Mean moderation | Sensitivity model (`docs/38`) |
|-----------|-----------|-----------|-------------|-----------|-------------|
| 1 Decomposition | BL, TV, PB, BR | same | same | same | same |
| 2 Onset | Gompertz in $t_{od}$ | Gompertz | Gompertz | Gompertz | turnover $f(t)$ |
| 3 Offset | compounding recursion | compounding | anchored | anchored or compounding | anchored (turnover) |
| 4 Re-exposure | restart at 0 | restart at 0 | restart at 0 (option: residual) | restart at 0 | from residual |
| 5 Biomarker encoding | covariance, step profile | covariance, graded | covariance, graded | mean, on-drug only (or decayed) | mean through sensitivity, proportional to $f$ |
| 6 Prognostic pathways | none | placebo only | natural course, placebo, baseline (driver) | possible | possible |
| 7 Dependence | CS, $\rho = 0.8$ | non-separable AR(1) | separable AR(1), $\rho = 0.7$ | as chosen | rank-one plus serial |
| 8 Response heterogeneity | constant variance, full off drug | same | same | same | shrinks with the response |
| 9 Observation | censoring sets | censoring sets | phase error, MCAR/MAR (driver) | as chosen | as chosen |
| 10 Validity | silent repair | silent repair in sampler | stop on infeasible | always valid | valid by construction |
| 11 Biomarker measurement | exact, normal | exact, normal | exact, normal | exact; any distribution possible | exact; any distribution possible |
| 12 Form of moderation | linear | linear | linear | linear | linear in $\omega_i$ |
| 13 Exposure | full adherence; no titration | same | same | same | same (exposure could enter the turnover input) |
| 14 Blinding | $e_t$ fixed by design | same | same | same | same |
| 15 Carryover variation | one half-life | same | same | same | one turnover rate |
| 16 Selection | none | none | none | none | none |
| 17 Outcome scale | unbounded Gaussian | same | same | same | same |

Table: The five DGPs compared on all seventeen features

*Driver: implemented in the simulation drivers, not yet in the package.*

The five DGPs agree on every row from 11 to 17. Those features do not
discriminate among the existing processes; they mark limitations that
all of them share.

## 21. Summary for the introduction of paper 01

Monte Carlo studies of aggregated N-of-1 trials can only show what
their data-generating process lets the data contain. We identify
seventeen features of such processes. The first ten are the features
on which existing processes differ:

1. the decomposition of the score into baseline, natural course,
   expectancy and drug responses;
2. the onset of the drug response;
3. its offset and carryover;
4. the response on re-exposure;
5. the encoding of the predictive biomarker;
6. prognostic pathways for that biomarker;
7. the within-person dependence over visits and across components;
8. between-person heterogeneity of response;
9. the observation process, including measurement error and
   missingness;
10. the validity, estimand and calibration of the process.

The remaining seven are simplifications that every published process
we examined shares:

11. a biomarker measured without error;
12. moderation linear in the biomarker;
13. exposure equal to assignment;
14. intact blinding;
15. one carryover half-life for all participants;
16. no selection on a screening score;
17. an unbounded outcome scale.

Three findings motivate the list. First, the power collapse under small
carryover reported by Hendrickson et al. is produced jointly by
Features 3 and 5. Their biomarker coupling stays at full strength
wherever any drug effect remains, so the simulated interaction vanishes
at any positive half-life. A coupling that fades with the drug effect
removes the collapse. Second, Feature 7 bounds Feature 5. Under
compound symmetry the covariance encoding cannot represent biomarker
correlations above 0.256, and the published effect sizes relied on
silently repaired matrices; a separable AR(1) structure is valid for
every schedule and raises the bound. Third, Feature 6 decides which
analysis is valid. If the biomarker is also prognostic, analyses that
share its effect with the baseline visit report spurious interactions.
We use the first ten features to specify a refined process and to
compare it with the published one, and treat the remaining seven as
stated limitations and directions for stress testing.

## 22. Limitations

- **The features are drawn from this compendium's experience** with one
  disease area (PTSD) and one drug (prazosin). Other settings may add
  features, for example binary or count outcomes, or multiple
  biomarkers.
- **Some recommendations rest on derivations** and have not been
  simulated: the turnover onset, restart from the residual, and the
  sensitivity model.
- **Features 11 to 17 have not been simulated in this compendium.**
  Their stated effects are first-order derivations (attenuation by
  reliability, the linear projection of a threshold, scaling by
  adherence, the unblinding bias, range restriction) or arguments; the
  derivations assume linear models and normal distributions.
- **"Every published process we examined"** (Section 21) refers to the
  five DGPs of Section 20, not to a systematic review of the
  literature.
- **Parameter calibration to trial data is outstanding.**

## 23. References

Each reference was verified against PubMed (PMID given) or the
publisher's page, on 2026-10-07. Exceptions are noted.

1. Morris TP, White IR, Crowther MJ. Using simulation studies to
   evaluate statistical methods. *Stat Med* 2019;38(11):2074-2102.
   doi:10.1002/sim.8086. PMID 30652356.
2. Burton A, Altman DG, Royston P, Holder RL. The design of simulation
   studies in medical statistics. *Stat Med* 2006;25(24):4279-4292.
   doi:10.1002/sim.2673. PMID 16947139.
3. Holford N. Clinical pharmacology = disease progression + drug
   action. *Br J Clin Pharmacol* 2015;79(1):18-27. doi:10.1111/bcp.12170.
   PMID 23713816.
4. Chan PL, Holford NHG. Drug treatment effects on disease progression.
   *Annu Rev Pharmacol Toxicol* 2001;41:625-659.
   doi:10.1146/annurev.pharmtox.41.1.625. PMID 11264471.
5. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Front Digit Health*
   2020;2:13. doi:10.3389/fdgth.2020.00013. PMID 34713026.
6. Sheiner LB, Stanski DR, Vozeh S, Miller RD, Ham J. Simultaneous
   modeling of pharmacokinetics and pharmacodynamics: application to
   d-tubocurarine. *Clin Pharmacol Ther* 1979;25(3):358-371.
   doi:10.1002/cpt1979253358. PMID 761446.
7. Dayneka NL, Garg V, Jusko WJ. Comparison of four basic models of
   indirect pharmacodynamic responses. *J Pharmacokinet Biopharm*
   1993;21(4):457-478. doi:10.1007/BF01061691. PMID 8133465.
8. Holford NHG, Sheiner LB. Understanding the dose-effect relationship:
   clinical application of pharmacokinetic-pharmacodynamic models.
   *Clin Pharmacokinet* 1981;6(6):429-453.
   doi:10.2165/00003088-198106060-00002. PMID 7032803.
9. Hobbs DC, Twomey TM, Palmer RF. Pharmacokinetics of prazosin in man.
   *J Clin Pharmacol* 1978;18(8-9):402-406.
   doi:10.1002/j.1552-4604.1978.tb02456.x. PMID 690251.
10. Senn S. *Cross-over Trials in Clinical Research*, 2nd ed.
    Chichester: Wiley; 2002. (Publisher page; DOI not confirmed.)
11. Jones B, Kenward MG. *Design and Analysis of Cross-Over Trials*, 3rd
    ed. Boca Raton: Chapman and Hall/CRC; 2014. doi:10.1201/b17537.
12. Shi D, Ye T. Behavioral carry-over effect and power consideration
    in crossover trials. *Biometrics* 2024;80(2):ujae023.
    doi:10.1093/biomtc/ujae023. PMID 38563531.
13. Duan N, Kravitz RL, Schmid CH. Single-patient (n-of-1) trials: a
    pragmatic clinical decision methodology for patient-centered
    comparative effectiveness research. *J Clin Epidemiol*
    2013;66(8 Suppl):S21-S28. doi:10.1016/j.jclinepi.2013.04.006.
    PMID 23849149.
14. Kent DM, Paulus JK, van Klaveren D, et al. The Predictive Approaches
    to Treatment effect Heterogeneity (PATH) statement. *Ann Intern Med*
    2020;172(1):35-45. doi:10.7326/M18-3667. PMID 31711134.
15. Kent DM, Steyerberg E, van Klaveren D. Personalized evidence based
    medicine: predictive approaches to heterogeneous treatment effects.
    *BMJ* 2018;363:k4245. doi:10.1136/bmj.k4245. PMID 30530757.
16. Senn S. Mastering variation: variance components and personalised
    medicine. *Stat Med* 2016;35(7):966-977. doi:10.1002/sim.6739.
    PMID 26415869.
17. Araujo A, Julious S, Senn S. Understanding variation in sets of
    N-of-1 trials. *PLoS ONE* 2016;11(12):e0167167.
    doi:10.1371/journal.pone.0167167. PMID 27907056.
18. Liu GF, Lu K, Mogg R, Mallick M, Mehrotra DV. Should baseline be a
    covariate or dependent variable in analyses of change from baseline
    in clinical trials? *Stat Med* 2009;28(20):2509-2530.
    doi:10.1002/sim.3639. PMID 19610129.
19. Diggle PJ. An approach to the analysis of repeated measurements.
    *Biometrics* 1988;44(4):959-971. PMID 3233259. (DOI 10.2307/2531727
    from a secondary source.)
20. Verbeke G, Molenberghs G. *Linear Mixed Models for Longitudinal
    Data*. New York: Springer; 2000.
21. Jones RH, Boadi-Boateng F. Unequally spaced longitudinal data with
    AR(1) serial correlation. *Biometrics* 1991;47(1):161-175.
    PMID 2049497. (DOI 10.2307/2532504 from a secondary source.)
22. Galecki AT. General class of covariance structures for two or more
    repeated factors in longitudinal data analysis. *Commun Stat Theory
    Methods* 1994;23(11):3105-3119. doi:10.1080/03610929408831436.
23. Senn S. Sample size considerations for n-of-1 trials. *Stat Methods
    Med Res* 2019;28(2):372-383. doi:10.1177/0962280217726801.
    PMID 28882093.
24. Little RJA, Rubin DB. *Statistical Analysis with Missing Data*, 3rd
    ed. Hoboken: Wiley; 2019. doi:10.1002/9781119482260.
25. Molenberghs G, Kenward MG. *Missing Data in Clinical Studies*.
    Chichester: Wiley; 2007. doi:10.1002/9780470510445.
26. Raskind MA, Millard SP, Petrie EC, et al. Higher pretreatment blood
    pressure is associated with greater posttraumatic stress disorder
    symptom reduction in soldiers treated with prazosin. *Biol
    Psychiatry* 2016;80(10):736-742. doi:10.1016/j.biopsych.2016.03.2108.
    PMID 27320368.
27. Gartner T, Schneider J, Arnrich B, Konigorski S. Comparison of
    Bayesian networks, G-estimation and linear models to estimate causal
    treatment effects in aggregated N-of-1 trials with carry-over
    effects. *BMC Med Res Methodol* 2023;23:191.
    doi:10.1186/s12874-023-02012-5. PMID 37605171.
28. Gomeni R, Lavergne A, Merlo-Pich E. Modelling placebo response in
    depression trials using a longitudinal model with informative
    dropout. *Eur J Pharm Sci* 2009;36(1):4-10.
    doi:10.1016/j.ejps.2008.10.025. PMID 19041717.
29. Colloca L, Lopiano L, Lanotte M, Benedetti F. Overt versus covert
    treatment for pain, anxiety, and Parkinson's disease. *Lancet
    Neurol* 2004;3(11):679-684. doi:10.1016/S1474-4422(04)00908-1.
    PMID 15488461.
30. Rutherford BR, Roose SP. A model of placebo response in
    antidepressant clinical trials. *Am J Psychiatry*
    2013;170(7):723-733. doi:10.1176/appi.ajp.2012.12040474.
    PMID 23318413.
31. Wang Y, Schork NJ. Power and design issues in crossover-based N-of-1
    clinical trials with fixed data collection periods. *Healthcare
    (Basel)* 2019;7(3):84. doi:10.3390/healthcare7030084. PMID 31269712.
32. Blackston JW, Chapple AG, McGree JM, McDonald S, Nikles J. Comparison
    of aggregated N-of-1 trials with parallel and crossover randomized
    controlled trials using simulation studies. *Healthcare (Basel)*
    2019;7(4):137. doi:10.3390/healthcare7040137. PMID 31698799.
33. Tang J, Landes RD. Some t-tests for N-of-1 trials with serial
    correlation. *PLoS ONE* 2020;15(2):e0228077.
    doi:10.1371/journal.pone.0228077. PMID 32017772.
34. Chen X, Chen P. A comparison of four methods for the analysis of
    N-of-1 trials. *PLoS ONE* 2014;9(2):e87752.
    doi:10.1371/journal.pone.0087752. PMID 24503561.
35. Vohra S, Shamseer L, Sampson M, et al. CONSORT extension for
    reporting N-of-1 trials (CENT) 2015 Statement. *BMJ* 2015;350:h1738.
    doi:10.1136/bmj.h1738. PMID 25976398.
36. Hutcheon JA, Chiolero A, Hanley JA. Random measurement error and
    regression dilution bias. *BMJ* 2010;340:c2289.
    doi:10.1136/bmj.c2289. PMID 20573762.
37. Gail M, Simon R. Testing for qualitative interactions between
    treatment effects and patient subsets. *Biometrics*
    1985;41(2):361-372. PMID 4027319. (JSTOR stable 2530862; DOI not
    confirmed.)
38. Raskind MA, Peskind ER, Chow B, et al. Trial of prazosin for
    post-traumatic stress disorder in military veterans. *N Engl J Med*
    2018;378(6):507-517. doi:10.1056/NEJMoa1507598. PMID 29414272.
39. Le Flohic E, Vrijens B, Hiligsmann M. The impacts of undetected
    nonadherence in phase II, III and post-marketing clinical trials:
    an overview. *Br J Clin Pharmacol* 2024;90(8):1984-2003.
    doi:10.1111/bcp.16089. PMID 38752447.
40. Bang H, Ni L, Davis CE. Assessment of blinding in clinical trials.
    *Control Clin Trials* 2004;25(2):143-156.
    doi:10.1016/j.cct.2003.10.016. PMID 15020033.
41. Holford N, Ma SC, Ploeger BA. Clinical trial simulation: a review.
    *Clin Pharmacol Ther* 2010;88(2):166-182.
    doi:10.1038/clpt.2010.114. PMID 20613720. (Cited for general
    practice; its treatment of between-subject variability was not
    checked against the full text.)
42. Barnett AG, van der Pols JC, Dobson AJ. Regression to the mean: what
    it is and how to deal with it. *Int J Epidemiol* 2005;34(1):215-220.
    doi:10.1093/ije/dyh299. PMID 15333621.
43. Rogers JA, Polhamus D, Gillespie WR, et al. Combining patient-level
    and summary-level data for Alzheimer's disease modeling and
    simulation: a beta regression meta-analysis. *J Pharmacokinet
    Pharmacodyn* 2012;39(5):479-498. doi:10.1007/s10928-012-9263-3.
    PMID 22821139.
