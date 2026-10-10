# Understanding the BR-PB-TV Decomposition {.unlisted .unnumbered}

## An Intuitive Guide to Treatment Response in N-of-1 Trials {.unlisted .unnumbered}

*2026-05-07 08:51 PDT*

**Author.** pmsimstats team

**Purpose.** An introductory guide to the three-component
decomposition of the symptom score into a drug response ($BR$), a
placebo-belief response ($PB$) and a natural-course response ($TV$),
as used by the package's data-generating process (DGP). It explains
what each component represents, how the Hybrid design identifies it,
and why the analysis model deliberately does not mirror the DGP. The
worked examples are illustrative constructions, not simulation
output.

**Revision 2026-10-07.** Reviewed against
`analysis/report/NOTATION.md` and the current package source
(`R/generateData.R`, `R/buildtrialdesign.R`, `R/utilities.R`). The
Gompertz formula, its parameter interpretations and landmarks, the
Hybrid design schedule and expectancies, the size of the joint
covariance matrix, the description of participant heterogeneity, and
the variance-of-a-sum identity were corrected; symbols were renamed to
the canonical ones. Evidence labels used: *inspected* (read in the
source), *derived*, *argued*, and *illustrative*.

**Revision 2026-10-08.** Two statements were corrected against
`analysis/report/compendium-mathematics.tex`: the reset of $t_{pb}$
at visits with $\eta_t = 0$ (Component 2), and the count of
identified parameters in the phase-indicator extension, which is
five, not six, because the full formula is exactly aliased (section
"A small extension that usually pays off").

---

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. Symbols marked *local*
are defined in this guide only.

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $i$, $t$ | participant; visit (week) | canonical |
| $Y_{it}$ | symptom score (code column `Sx`) | canonical |
| $\mathrm{BL}_i$ | baseline symptom level | canonical |
| $BR_{it}$, $PB_{it}$, $TV_{it}$ | drug, placebo-belief and natural-course responses (reductions in symptoms) | canonical |
| $\varepsilon_{it}$ | observation-level residual | canonical |
| $u_i$, $v_i$ | participant random intercept; participant random slope | canonical; local |
| $\beta_0$, $\beta_t$ | fixed intercept; fixed linear time slope | local |
| $D_{bc,it}$ | exposure-decayed drug indicator, 1 on drug and $e^{-\lambda t_{sd}}$ off drug (code `Dbc`) | canonical |
| $t_{od}$, $t_{sd}$, $t_{pb}$ | time on drug; time since discontinuation; expectancy-weighted time (code `tod`, `tsd`, `tpb`) | canonical |
| $t_{1/2}$, $\lambda$ | carryover half-life; decay rate $\ln 2/t_{1/2}$ | canonical |
| $\eta$ | expectancy weight of a visit (code `e`, the `expectancies` argument): 1 open label, 0.5 blinded | canonical (paper 06, $\eta(\text{phase})$) |
| $m$, $r$, $d$ | Gompertz maximum, rate and displacement (code `maxr`, `rate`, `disp`); subscripts $BR$, $PB$, $TV$ name the component | local |
| $\sigma_{PB}$ | open-label SD of the placebo-belief response | local |
| $c_{cf1t}$, $c_{cfct}$ | cross-component correlations at the same visit and at different visits (code `c.cf1t`, `c.cfct`) | local |
| $\rho$, $w_t$ | AR(1) serial correlation per week; calendar week of visit $t$ | canonical; local |
| $K$ | number of visits | local |
| $I_\ell$ | indicator of phase level $\ell$ (open label OL is the reference) | local |
| $\beta_{bm:D}$ | biomarker-by-treatment interaction, the estimand (`bm:Dbc`) | canonical |
| $\beta_{bm}^{\text{lumped}}$ | biomarker slope of a one-component (lumped) regression | canonical |
| $\beta_{bm}^{BR}$, $\beta_{bm}^{PB}$, $\beta_{bm}^{TV}$ | component-specific biomarker slopes | canonical |
| $N$ | total number of participants across randomization paths | canonical |

Table: Notation index: symbols, meanings and status

**Glossary.**

- **Component.** One of the three latent responses whose sum, subtracted
  from baseline, gives the symptom score.
- **Expectancy.** The design-fixed weight of the participant's belief
  that they are receiving active drug at a visit.
- **Open label (OL), blinded discontinuation (BDC), crossover (CO).**
  The three phases of the Hybrid design.
- **Identifiability.** Whether the data, under a given design, contain
  enough contrasts to estimate a component separately from the others.
- **Predictive / prognostic biomarker.** Associated with the size of
  the drug effect, or with the outcome irrespective of treatment.

## The Big Picture: Why Split Responses?

Imagine you are a PTSD patient trying prazosin (a blood-pressure
medication that also reduces nightmares). After eight weeks of
treatment, your symptom score has dropped by twelve points. The
clinical question is straightforward: 'is the drug working?' The
statistical question is harder: 'how much of that twelve-point drop
should I attribute to prazosin, and how much to other things?'

Three plausible mechanisms produce a drop in symptom score on a
clinical trial:

1. **Drug effect.** The medication chemically modulates the
   physiology of nightmares. Prazosin is an alpha-1-adrenergic
   antagonist that reduces noradrenergic tone during REM sleep,
   and its pharmacological action is independent of whether you
   know you are taking it.
2. **Expectation effect.** You believe the treatment will work, and
   your belief activates well-characterized neurobiological
   pathways: dopaminergic reward signaling, descending pain
   modulation, anticipatory cortisol regulation. The placebo
   literature is now thirty years deep on this; it is not 'all in
   your head' in any dismissive sense.
3. **Natural history.** Some symptoms wax and wane on their own.
   PTSD severity drifts with anniversaries, life events, sleep
   hygiene, and the slow erosion that any chronic condition shows
   when measured repeatedly over months.

If we lump these three causes together as a single 'treatment
response', we cannot distinguish them. The patient's twelve-point
drop might be eight points of drug, three of expectation, and one
of natural drift; or it might be two of drug, eight of expectation,
and two of drift. The total is the same. The clinical implications
are completely different.

The BR-PB-TV decomposition addresses this directly. It writes the
observed change in outcome as a sum of three biologically
interpretable components plus measurement noise:

$$
Y_{it} \;=\; \mathrm{BL}_i - [\,BR_{it} + PB_{it} + TV_{it}\,]
            + \varepsilon_{it}.
$$

Here $\mathrm{BL}_i$ is participant $i$'s baseline symptom score,
$BR$, $PB$, and $TV$ are the three response components defined
below (each a non-negative reduction in symptoms), and
$\varepsilon$ is observation-level noise. By decomposing the
response into these three causes, we can ask cause-specific
questions that the lumped form cannot answer: how much of the
improvement is pharmacological, how much is expectation-driven,
and how much would have occurred without any treatment at all.

### The problem with one-component models

A naive analysis writes:

```
improvement = treatment_effect + error
```

The mathematical mistake is subtle but consequential. The 'error'
term in such a model becomes a wastebasket for everything the
analyst did not measure. If improvement is, in truth,

```
improvement = 4 (drug) + 6 (expectation) + 1 (natural history) + 1 (noise) = 12
```

the naive model fits this as

```
improvement = 12 (treatment_effect) + 0 (error)
```

at the population mean and infers that the drug is three times more
powerful than it actually is. The structure is missed, and the
inference is biased toward the lumped quantity rather than the
quantity of clinical interest. Worse, the lumped quantity does not
generalize: a population with stronger placebo expectations would
appear to have a more powerful drug, and a population with milder
natural-history drift would appear to have a less powerful one.
Neither inference would be correct.

Decomposing the response separates these failure modes.

---

## Component 1. Biological Response (BR)

### What it is

BR is the **physiological response to the drug itself.** It is what
happens because of the medication's chemical action on the body,
independent of belief, suggestion, or natural history. A
participant who somehow received prazosin without knowing it would
still develop a measurable BR; a participant who took an inert
sugar pill while believing it was prazosin would have a BR of zero
no matter how strong their expectation.

Methodologically, BR is the quantity that drug regulators and
prescribers actually want to know. Approval decisions, dose
selection, and biomarker-guided patient selection are all framed
around BR. The other two components are real and inferentially
relevant, but they are not the substance of pharmacological efficacy.

### Real example: PTSD and prazosin

Prazosin is an alpha-1-adrenergic antagonist originally developed
for hypertension. It crosses the blood-brain barrier, and at the
relevant doses it antagonizes noradrenergic activity in the locus
coeruleus and prefrontal cortex during REM sleep. The mechanism for
nightmare reduction is reasonably specific: prazosin damps the
adrenergic surge that accompanies trauma-related dream content,
and the resulting suppression of REM-stage arousal reduces
nightmare frequency and severity.

The pharmacokinetic profile is short. Prazosin has a plasma
half-life of about 2.3 hours (Hobbs et al. 1978); brain concentrations equilibrate quickly with
plasma. This produces a characteristic on-drug / off-drug pattern:

- **Week 1 on drug.** Brain prazosin reaches steady-state within
  the first few doses. Nightmare scores begin to fall. The first-
  week reduction is typically small because the patient has not
  yet had many drug-on REM cycles to integrate.
- **Weeks 2-3 on drug.** Effect accumulates as the patient
  experiences several weeks of drug-on REM sleep. The reduction
  approaches its plateau.
- **Steady state.** Once the brain has fully adapted, additional
  weeks on drug produce essentially no further reduction. This is
  the saturating-curve property: the body has a ceiling for how
  much benefit any given drug can produce, and that ceiling is
  reached on a timescale set by the drug's own pharmacology and
  the disease's response time.
- **Discontinuation.** Within a few half-lives, plasma and brain
  concentrations return to negligible. BR drops to zero. The patient
  no longer has any pharmacological effect.

This pattern (saturating rise, plateau, sharp washout) is what BR
captures. It is **not** linear in time; it is **not** symmetric on
and off drug; and it is **not** modulated by what the patient
believes. Those properties are the empirical signatures that let us
identify BR separately from PB and TV.

### Mathematical pattern: modified Gompertz

BR follows a saturating curve known as the modified Gompertz
function, which has the form

$$
BR(t_{od}) \;=\; m\,\frac{\exp\bigl(-d\,e^{-r t_{od}}\bigr) - e^{-d}}
                        {1 - e^{-d}},
$$

a Gompertz curve in time on drug $t_{od}$, shifted so that
$BR(0) = 0$ and rescaled so that $BR(t_{od}) \to m$ as
$t_{od} \to \infty$ (inspected: `modgompertz()` in `R/utilities.R`).
The three parameters have direct interpretations:

| Parameter | Symbol | Meaning |
|---|---|---|
| Maximum response | $m$ | The biological ceiling. The largest reduction the drug can produce, even with infinite time on treatment. |
| Onset rate | $r$ | How fast the patient approaches the ceiling, given $d$. Higher $r$ compresses the whole curve in time. |
| Displacement | $d$ | The length of the initial lag. Large $d$ gives a sigmoid curve with a slow start and an inflection at $t_{od} = \ln d / r$; small $d$ (below 1) gives a curve that rises fastest at the start, with no lag. |

Table: Modified Gompertz parameters for the BR component and their meanings

For PTSD prazosin, the values used by Hendrickson et al. (2020) and
throughout this compendium are $m = 10.99$ (about eleven points of
symptom reduction at saturation), $d = 5$ and $r = 0.42$ per week.
Their landmarks (derived from the formula; `docs/40`, Section 4):
the inflection is at 3.8 weeks, half of the maximum is reached at
4.7 weeks (about 33 days), and 90% at 9.2 weeks. This is a slow,
sigmoid onset.

The Gompertz family was chosen because it captures the
'asymptote-from-below' shape that physiological dose-response
curves typically display. Two alternatives that do **not** capture
this shape:

- A linear function in $t$ ($BR = \beta_t \cdot t$) climbs without
  bound and is biologically implausible past a few weeks.
- An exponential approach to the ceiling ($BR = m(1 - e^{-rt})$)
  is acceptable for fast-onset drugs but lacks the slow-start phase
  that characterizes most chronic-condition responses, where
  receptor adaptation, feedback regulation, and downstream
  gene-expression changes produce a measurable lag.

The modified Gompertz contains the exponential approach as its
small-$d$ limit (for small $d$, $\exp(-d e^{-rt}) \approx
1 - d e^{-rt}$, so the rescaled curve tends to $m(1 - e^{-rt})$;
derived), and becomes an increasingly delayed sigmoid as $d$ grows.
Its three-parameter form is flexible enough to describe a range of
pharmacodynamic profiles (argued; no fit to individual-level data
is reported in this guide).

### Why this matters for design

If we knew the participant's BR curve in isolation, we could read
off the drug effect directly: $BR$ at any timepoint is the answer
to 'how much pharmacological reduction has accumulated by now?'.
The trick of N-of-1 trial design is to construct contrasts between
on-drug and off-drug timepoints that isolate $BR$ from $PB$ and
$TV$. The blinded discontinuation phase, discussed in the trial-
design section below, is precisely such a contrast: if the
participant is taken off drug while still believing they are on
drug, then $BR$ collapses to zero while $PB$ remains roughly intact.
The difference between on-drug and off-drug response in that
window is therefore a clean estimate of $BR$ that does not require
the analyst to model $PB$ or $TV$ explicitly.

---

## Component 2. Placebo-Belief Response (PB)

### What it is

PB is the change in symptom score that arises from the patient's
**belief** that they are receiving treatment, independent of what
they are actually receiving. This is the placebo response in the
strict sense: the part of the improvement that would persist if
the active drug were silently replaced with an inert sugar pill,
provided the patient continued to believe they were on the active
drug.

It is essential to be precise about what PB is and is not. PB is
not 'fake' improvement, and it is not 'all in the patient's head'
in any dismissive sense. The placebo response is mediated by well-
characterized neurobiological pathways. Belief-driven release of
endogenous opioids and dopamine produces measurable analgesia and
mood elevation. Anticipatory regulation of the
hypothalamic-pituitary-adrenal axis modulates cortisol and
inflammatory cytokines. Top-down attentional shifts change how
symptoms are perceived and reported. Each of these mechanisms is a
real biological response; what makes them 'placebo' rather than
'drug' is that they are triggered by belief rather than by
pharmacology.

For trial design, the relevant point is that PB is **separately
identifiable** from BR if and only if the trial design includes a
phase in which the patient's belief about treatment differs from
the actual treatment they are receiving. The blinded discontinuation
phase of the Hendrickson Hybrid design is the canonical such phase.

### Real example: the blinded discontinuation contrast

Consider the Hendrickson Hybrid design as implemented in this
compendium (inspected: `analysis/scripts/carryover-sensitivity/
simulation-core.R`), a twenty-week N-of-1 trial with visits at weeks
4, 8, 9, 10, 11, 12, 16 and 20, structured in three phases:

- **Weeks 1-8 (open-label on drug).** The patient knows they are
  receiving prazosin. Belief in treatment is at full strength
  (expectancy weight $\eta = 1$, in the model parameterization
  introduced below).
- **Weeks 9-12 (blinded discontinuation).** The patient is told
  that, at some unannounced point during this window, they may be
  silently switched to placebo. They cannot reliably tell whether
  they are still receiving the active drug. Belief in treatment
  drops to a partial level ($\eta = 0.5$) because the patient's
  expectation is now hedged. Everyone is on drug at week 9; the drug
  stops after week 9 in two of the four randomization paths and
  after week 10 in the other two.
- **Weeks 13-20 (blinded crossover).** Four weeks on drug and four
  on placebo, in randomized order. The crossover is also blinded,
  so $\eta = 0.5$ throughout it (`docs/40`, Section 11).

At week 10, a patient whose drug stopped after week 9 will have:

- $BR = 0$, or only its carryover residual (no pharmacology).
- $PB = 0.5 \cdot PB_{\text{open-label}}$, reflecting hedged belief.
- $TV$ continuing on its slow trajectory.

At the same visit, a patient still on active drug will have:

- $BR$ continuing at its near-plateau value.
- The same hedged $PB = 0.5 \cdot PB_{\text{open-label}}$ as the
  placebo-switched patient.
- The same $TV$ trajectory.

The contrast between the two groups of paths at that visit is
therefore a nearly clean estimate of $BR$, with $PB$ matched (because both
groups have the same belief) and $TV$ matched (because $TV$ is by
construction independent of treatment). This is the inferential
move that makes the BR-PB-TV decomposition identifiable from the
data: the trial design supplies phases in which different
combinations of components are present, and the analyst recovers
the components by differencing those phases.

### Mathematical pattern

PB also follows a modified Gompertz curve, with its own parameters,
in the expectancy-weighted time $t_{pb}$ (time accumulated over
consecutive visits with positive expectancy; it resets to zero at any
visit with $\eta_t = 0$), and its amplitude is scaled by the
expectancy weight $\eta_t$ of the current visit:

$$
PB_t \;=\; \eta_t \, m_{PB}\,
\frac{\exp\bigl(-d_{PB}\,e^{-r_{PB} t_{pb}}\bigr) - e^{-d_{PB}}}
     {1 - e^{-d_{PB}}}
$$

(inspected: `buildSigma()` in `R/generateData.R` multiplies the
Gompertz mean by the design's `e`). (Corrected 2026-10-08: this
sentence previously described $t_{pb}$ as time accumulated while the
expectancy is positive. `buildtrialdesign()` in `R/buildtrialdesign.R`
resets `tpb` to zero at every visit with `e = 0` and restarts the
accumulation afterward (inspected), which matters only for designs
with covert visits, such as a balanced-placebo design. See
`analysis/report/compendium-mathematics.tex`, Section "Time arguments
of the components, the shared clock, and the phase contrasts"
(`sec:invb-clocks`).) Values of the
expectancy weight:

| Phase | $\eta$ | Reasoning |
|---|---|---|
| Open-label on drug | 1.0 | Patient knows they are on the active drug. |
| Blinded discontinuation and blinded crossover | 0.5 | Patient believes they may or may not be on drug; expectation is hedged. |
| Open-label placebo | 0.0-0.5 | Patient told they are on placebo retains some residual expectation. Not used by any design in this compendium. |

Table: Expectancy weight $\eta$ by trial phase with reasoning

The choice of $\eta = 0.5$ rather than $\eta = 0$ for the blinded
phases is itself a modeling assumption, inherited from Hendrickson
et al. (2020). The empirical question is how much expectation a
patient carries when uncertain about allocation. This guide does not
cite a source for a specific range, and none was located in the
project bibliography (flagged 2026-10-07; an earlier version quoted
'30-60% of the open-label level' without a source). Setting
$\eta = 0.5$ in the blinded phases is a convenient default, but
sensitivity analyses should vary it to confirm that conclusions are
not driven by this single number.

### Why expectancy variance differs by phase

A subtle property of $PB$ is that **its variance also scales with
$\eta$**, not just its mean. When the expectancy factor is high
(open-label on drug), there is substantial variability across
patients in how strongly the placebo response is expressed: high-
suggestibility patients may show large $PB$, low-suggestibility
patients may show essentially none, and the population variance is
correspondingly large. When $\eta$ is low (blinded), the placebo
response is uniformly attenuated, and the patient-to-patient
variance shrinks accordingly. The model parameterization reflects
this with a visit-dependent standard deviation (inspected:
`buildSigma()` scales the PB SD by `e`):

$$
\mathrm{SD}(PB_t) \;=\; \eta_t \, \sigma_{PB}.
$$

This heteroscedastic structure matters for the analysis (argued): a
working model with a common residual variance misstates the
variability of the open-label visits (where $PB$ contributes most
variance) relative to the blinded ones, with consequences for both
confidence intervals on the drug effect and the calibration of any
subsequent biomarker test. The stress test of `docs/46` includes
analyses with phase-specific residual variances for this reason.

---

## Component 3. Time-Varying Response (TV)

### What it is

TV is the change in symptom score that would have happened **even
if the patient had received no treatment at all.** It is the
participant's natural-history trajectory, integrated over the
duration of the trial. For some participants TV is positive
(symptoms naturally improve over the course of the trial because
of life events, therapy, time-in-condition, or regression to the
mean); for others TV is negative (symptoms naturally worsen because
of trauma anniversaries, accumulating stress, or progressive
underlying disease).

Like PB, TV is not 'noise'. It is a structured, participant-specific
signal that the model can in principle estimate from data, given
sufficient repeated measurements and a trial design that includes
a long enough off-drug baseline or follow-up window. What
distinguishes TV from PB is its origin: TV is independent of
treatment status and of belief about treatment, whereas PB is
modulated by both.

### When does natural history actually exist?

A common simplifying assumption in trial analysis is 'subjects
will improve over time, with or without treatment'. This is true
for some conditions and dangerously wrong for others. The shape
of TV depends on the underlying disease class, the trial timescale,
and the selection process by which participants entered the trial.
None of these can be assumed without thought, and the empirical
shape of TV is one of the things the decomposition is designed to
recover.

**Acute conditions tend to TV-positive.** A cold, a sprained ankle,
a post-operative pain syndrome, an episode of acute back pain:
these conditions resolve, on average, on a timescale of days to
weeks even without intervention. A trial of ibuprofen for tension
headache or of a topical analgesic for ankle sprain will see
substantial TV-improvement in the placebo arm, and the apparent
drug effect must be measured *against* this baseline trajectory.
Failing to model TV in an acute-condition trial is a common
source of inflated effect estimates.

**Chronic stable conditions tend toward TV-zero on the population
mean, with substantial individual variance.** Type 2 diabetes
under stable management, treated hypertension, controlled chronic
asthma, and long-standing PTSD without active triggers are all
reasonably steady on the population mean across a sixteen-week
trial. But individual participants show substantial variance: one
patient is doing better this month because they have started
sleeping more regularly, another is doing worse because their job
has become more stressful. The population mean of TV may be near
zero, but the participant-level TV variance is not, and the
decomposition recovers both.

**Chronic progressive conditions tend toward TV-negative.**
Alzheimer's disease and other neurodegenerative conditions, ALS,
late-stage congestive heart failure, untreated progressive multiple
sclerosis, and many forms of cancer all worsen on the timescale of
months. A trial of a putatively neuroprotective intervention in
early Alzheimer's must explicitly model the negative TV
trajectory: the apparent treatment effect is the slowing of
progression relative to the worsening that would otherwise occur.
Treating TV as zero in such a trial would systematically
under-estimate efficacy, because the placebo arm continues to
worsen and the active arm's stability looks like 'no effect' rather
than 'effect is preventing the natural decline'.

**Cyclical and triggered conditions show structured TV that is
neither monotone nor stationary.** PTSD severity often spikes near
trauma anniversaries; major depressive disorder has a measurable
seasonal component; relapsing-remitting multiple sclerosis
oscillates on a timescale of months. In these conditions, TV is
not well-described by a monotonically increasing or decreasing
Gompertz curve, and the model needs either a more flexible TV
specification or an explicit covariate (e.g., season, days from
nearest trigger) that absorbs the cyclic component before TV
estimates the residual trend.

**Selection-induced TV is real and rarely zero.** Even a chronic
stable condition can show TV-positive average trajectory in a
trial population because of how participants are recruited.
Patients enroll at a moment of clinical engagement, which often
coincides with a recent worsening that triggered the
help-seeking behavior. Regression to the mean alone (the
statistical fact that an extreme observation tends to be followed
by less extreme observations on average) can produce several
points of apparent TV-improvement in the early weeks of a trial,
even in a population whose underlying condition is stable. This
is a recognized failure mode of single-arm trials in chronic-pain
conditions: large 'treatment' effects can shrink when a control arm
is added, because much of the apparent effect was regression to the
mean (argued; see Barnett et al. 2005 for the mechanism).

**Trial participation itself induces TV.** The Hawthorne effect
(participants improve simply because they are being observed),
the structured-visit effect (regular visits with clinicians may
improve symptom management directly), and the
diary-completion effect (the act of monitoring symptoms changes
how they are perceived and reported) are all real, all
participant-mediated, and all confounded with calendar time. They
appear in the data as TV-improvement that is independent of any
pharmacological or placebo response, and they apply equally to
all arms of the trial. Modeling them as part of TV is the
appropriate handling.

### Can we assume no TV effect?

The short answer is: rarely safely, and never without empirical
support. The longer answer:

- **TV may be assumed zero only if** the disease is in a stable
  chronic phase, the trial is short relative to the disease
  timescale, the recruitment process does not favor participants
  at peak severity, and the trial protocol does not introduce
  ancillary care that itself improves symptoms. For most chronic
  trials in psychiatric or pain conditions, at least one of these
  conditions fails. The conservative default is to estimate TV from
  the data, even at the cost of some statistical efficiency, rather
  than to assume it is zero.
- **The decomposition makes the assumption testable.** A model
  that includes a TV component will return $\hat{m}_{TV}$ and a
  variance estimate; if both are small relative to $BR$ and $PB$,
  the analyst has empirical evidence for the simplifying
  assumption. If $\hat{m}_{TV}$ is meaningfully large, the
  assumption was wrong, and the trial needs the TV component to
  produce unbiased inference about the others.
- **Setting TV to zero and treating it as part of $\varepsilon$
  is worse than estimating it.** The participant-level TV is
  structured (each participant has their own slow trajectory) and
  not independent of the on-drug versus off-drug schedule. Lumping
  TV into noise produces residuals that are autocorrelated within
  participants and heteroscedastic across them, and the resulting
  standard errors on the BR and PB estimates are wrong in both
  directions: too small for participants whose TV is small, too
  large for participants whose TV is large. Estimating TV
  explicitly fixes both problems simultaneously.

The general principle: TV is a *structured* component of the
response, not random noise, and the decomposition treats it
accordingly.

### Why this matters

If TV is ignored, every TV-induced point of symptom change gets
attributed to either BR or PB or to noise. This is the most
common sources of bias in chronic-condition trial inference
(argued). A typical pattern: an uncontrolled trial that enrolls
patients at peak symptom severity shows a large improvement that a
concurrent natural-history arm would reveal to be largely regression
to the mean; the patients would have improved similarly with no
intervention (illustrative; an earlier version of this guide
attributed this to 'early trials of chronic-fatigue therapies'
without a source, and that attribution has been removed). TV in the
model is the formal way to handle this risk.

The attribution problem is not symmetric. TV is positively
correlated with the *opportunity* to detect a drug effect (chronic-
condition patients are typically at peak severity at trial entry,
and natural-history regression brings their scores down on
average), and it is *negatively* correlated with the *credibility*
of the detected effect (a regulator cannot tell, from a
single-arm trial without a TV-control phase, how much of the
observed improvement was treatment versus regression). Both
problems are addressed by decomposing TV explicitly and by
including trial phases in which BR and PB are absent.

### Mathematical pattern

TV also follows a modified Gompertz curve, in time since trial entry
and parameterized independently of BR:

$$
TV_t \;=\; m_{TV}\,
\frac{\exp\bigl(-d_{TV}\,e^{-r_{TV} t}\bigr) - e^{-d_{TV}}}
     {1 - e^{-d_{TV}}} .
$$

The maximum $m_{TV}$ can be either positive (natural improvement)
or negative (natural worsening). In the package's DGP the Gompertz
parameters are population constants: the curve gives the mean of
$TV_t$, and participants differ through a multivariate-normal
deviation around that mean, with a common SD and AR(1) correlation
across visits (inspected: `buildSigma()`). Hendrickson et al. gave
$TV$ the same parameters as $PB$ ($m = 6.51$, $d = 5$, $r = 0.35$,
SD 10; `docs/40`, Section 11), so the split between the two is an
assumption rather than an estimate. A trial design that aggregates
over many participants can identify the population mean trajectory
and the participant-level variance, but only if there are enough
visits to constrain the curve.

### Why not just include 'week' as a covariate?

A natural objection is that a model with `week` as a covariate
already controls for time. Three reasons this is not enough:

1. **Linearity assumption.** A `week` covariate enforces a linear
   trajectory: every additional week contributes the same fixed
   amount of natural-history change. Real disease trajectories are
   not linear; they curve, plateau, or accelerate. The Gompertz
   form is flexible enough to capture S-shaped or saturating
   patterns that linear time cannot.
2. **Participant-specific TV.** A single `week` coefficient
   estimates a population-mean trajectory, but participant
   heterogeneity in natural history is large and clinically
   relevant. The full BR-PB-TV decomposition uses participant-
   specific TV components, drawn from a population distribution,
   which is the correct way to express the heterogeneity.
3. **Confounding with BR onset.** Both BR and TV evolve smoothly
   over the same trial timescale. A `week` covariate cannot tell
   them apart from observational data alone; the design contrasts
   (on-drug vs. blinded vs. off-drug phases) are what supply
   identifiability. Without the explicit decomposition, the
   `week` coefficient absorbs both, and the analyst loses any
   ability to attribute the resulting time-trend.

For trials with rich within-person time structure (and N-of-1
trials are precisely such trials), the explicit Gompertz TV is the
right tool, and the simpler `week` covariate is a coarse
approximation that loses information.

---

## How the trial design lets us identify the components

The mathematical decomposition above is meaningless unless the
data contain enough information to pin down each component
separately. Trial design is what supplies that information. The
Hendrickson Hybrid design, used as the worked example throughout
the rest of this guide, has three phases that each contribute a
different combination of components.

| Phase | Drug? | Belief? | Components present |
|---|---|---|---|
| Open-label on drug | Yes | Full | $BR + PB(\eta = 1) + TV$ |
| Blinded discontinuation | Hidden | Hedged | $BR \cdot \mathbf{1}_{\text{on drug}} + PB(\eta = 0.5) + TV$ |
| Blinded crossover | Either | Hedged | $BR \cdot \mathbf{1}_{\text{on drug}} + PB(\eta = 0.5) + TV$ |

Table: Components present in each phase of the Hendrickson Hybrid design

Three contrasts then bear on the three components:

- **Phase 1 vs. phase 2 (within drug):** the change from
  $PB(\eta = 1)$ to $PB(\eta = 0.5)$, about half of $PB$. This
  isolates the belief component only to the extent that $BR$ and
  $TV$ change little between week 8 and week 9.
- **Phase 2 or 3 on drug vs. placebo at the same visit:** $BR$ alone
  (because $PB$ and $TV$ match across paths at a visit).
- **Off-drug visits:** $TV + PB(\eta = 0.5)$, plus any carryover
  residual of $BR$. Every off-drug visit is blinded, so $TV$ is never
  observed without $PB$; and since the expectancy is positive at
  every visit, $t_{pb}$ equals calendar time and the two components
  share one clock. $TV$ and $PB$ are separated only through the drop
  in expectancy at week 9 and their assumed curve shapes (inspected;
  `docs/40`, Section 3).

Each contrast targets a different component, and the analyst
combines them through the linear-mixed-effects model fitted to the
full trajectory. The model is the formal mechanism, but the
inference is driven by the design contrasts. In the Hybrid design
the drug-versus-placebo contrast is the clean one; the $TV$ versus
$PB$ split is weakly identified.

A trial that lacks any of these phases loses the corresponding
identifiability. A pure open-label design (no blinded
discontinuation) cannot separate $BR$ from $PB$, because both are
present at full strength throughout. A pure blinded design (no
open-label phase) cannot estimate $\eta = 1$ $PB$, because no
phase has belief at full strength. The decomposition machinery is
useful only to the extent that the design supplies identifiability,
and trial designers should evaluate any proposed design by checking
that each component has at least one phase or contrast in which
it is uniquely present.

---

## Worked example: what's lost when you don't decompose

Consider two participants, A and B, enrolled in a twenty-week
N-of-1 trial of prazosin for PTSD. The example is illustrative: the
numbers are chosen by hand to make the arithmetic transparent, and
it simplifies the Hybrid schedule. Both are otherwise similar:
same age, same baseline severity, same trial design. They differ
in two latent traits that the trial does not measure directly.

- **Participant A** is a high placebo responder and a moderate
  pharmacological responder. Their true component values:
  $m_{BR}^A = 4$ (drug ceiling 4 points), $m_{PB}^A = 6$ (belief
  ceiling 6 points), $m_{TV}^A = 1$ (mild natural improvement).
- **Participant B** is a low placebo responder and a strong
  pharmacological responder. Their true component values:
  $m_{BR}^B = 8$, $m_{PB}^B = 1$, $m_{TV}^B = 1$.

In the open-label phase (weeks 1-8), both BR and PB run at full
strength, and TV accumulates slowly. By week 8, both participants
have reached close to their saturation values:

|  | $BR$ | $PB(\eta = 1)$ | $TV$ | Total improvement |
|---|---|---|---|---|
| Participant A | 4.0 | 6.0 | 1.0 | **11.0** |
| Participant B | 8.0 | 1.0 | 1.0 | **10.0** |

Table: Worked example components for participants A and B, open-label phase

A clinician examining only the totals would conclude that the
two participants are responding similarly: both improved by about
ten points, both look like clinical successes. A simple t-test on
the open-label change would estimate a population drug effect of
about ten points and would not distinguish the two participants.

Now consider the blinded discontinuation phase (weeks 9-12). Both
participants are silently switched to placebo early in the phase
(in the implemented design, after week 9 or week 10), and the
expectancy weight drops from 1.0 to 0.5. By the end of
week 12, the components are:

|  | $BR$ | $PB(\eta = 0.5)$ | $TV$ | Total improvement |
|---|---|---|---|---|
| Participant A | 0.0 | 3.0 | 1.0 | **4.0** |
| Participant B | 0.0 | 0.5 | 1.0 | **1.5** |

Table: Worked example components for participants A and B, blinded discontinuation phase

The two participants now look quite different. Participant A still
shows a sizeable improvement (four points) under blinded placebo;
participant B has nearly returned to baseline (one and a half
points). This is the diagnostic contrast that the decomposition
exploits.

What does each modeling strategy infer?

**One-component model.** A fixed-effects regression of total
symptom change on phase-by-allocation indicators sees:

- Open-label means: $11.0$ (A), $10.0$ (B).
- Blinded-placebo means: $4.0$ (A), $1.5$ (B).
- Estimated drug effect: $11.0 - 4.0 = 7.0$ (A), $10.0 - 1.5 =
  8.5$ (B). Population mean: $7.75$.
- Estimated placebo effect: not identifiable from this contrast.
- Estimated natural history: not identifiable from this contrast.

The one-component model has produced a single number (the drug
effect, or rather the mean of (BR + the drop in PB induced by
discontinuation)), conflated with the placebo response that washes
out under blinding. It cannot say whether the drug effect is the
same for participants A and B. It estimates 'drug' at 7.75 points
on average, which is an average of the true 4.0 (A) and 8.0 (B)
mixed with the placebo wash-out, and it has no principled way to
decompose this further.

**Three-component model.** A linear-mixed-effects model with the
full BR-PB-TV decomposition, fitted to the same data, can recover
the individual components by exploiting the differential expectancy
between phases and the on-drug-versus-blinded-placebo contrast.
In the idealized case of noise-free data and a correctly specified
model, the fit returns (illustrative, not a simulation result; with
real data each estimate carries sampling error, and the $TV$ versus
$PB$ split is weakly identified, as noted above):

- $\hat{m}_{BR}^A = 4.0$ (true 4.0). $\hat{m}_{BR}^B = 8.0$ (true
  8.0). The drug effect is identified separately for each
  participant.
- $\hat{m}_{PB}^A = 6.0$ (true 6.0). $\hat{m}_{PB}^B = 1.0$ (true
  1.0). The placebo response is identified, and the individual
  difference between A and B is preserved.
- $\hat{m}_{TV}^A = 1.0$, $\hat{m}_{TV}^B = 1.0$. Natural history
  is identified and matched, as it should be for participants who
  truly have similar trajectories.

### What is at stake

Compare the inferences:

| Question | One-component answer | Three-component answer |
|---|---|---|
| 'Does the drug work for A?' | 'Some response, magnitude unclear, ~7.75 points overall.' | 'Yes: $m_{BR}^A = 4.0$ points, distinct from her larger placebo response.' |
| 'Does the drug work for B?' | 'Some response, magnitude unclear, ~7.75 points overall.' | 'Yes: $m_{BR}^B = 8.0$ points, twice the magnitude of A.' |
| 'Are A and B different responders?' | 'No discernible difference; both improved by ~10 points.' | 'Yes: B has 2x A's BR but 1/6 of A's PB.' |
| 'If a biomarker predicts BR, is it useful here?' | Cannot answer (no $BR$ estimate). | Yes if the biomarker correlates with $\hat{m}_{BR}$. |
| 'Should A be continued on this drug long-term?' | Cannot say (the response could be largely placebo). | Probably yes (BR = 4 is real and pharmacological), but the larger placebo component will not persist if expectations decay. |
| 'What is the population mean drug effect?' | 7.75 points (biased; mixes BR and the discontinuation-induced drop in PB). | 6.0 points = $(4.0 + 8.0)/2$ (correctly the mean of $BR$). |

Table: One-component versus three-component answers to clinical questions

The one-component model has produced one biased average where the
three-component model has produced six clinically actionable
quantities. The information loss from lumping the components is not
abstract: it makes the difference between concluding 'these patients
respond similarly' and 'B has twice the pharmacological response
that A does, and a biomarker predicting BR would distinguish them
prospectively.' Predictive-biomarker validation, which is the whole
point of running an N-of-1 trial in this design family, becomes
impossible without the decomposition.

### A second example: the placebo-cohort fallacy

A subtler example shows the same issue at the population level.
Consider two cohorts of N-of-1 participants enrolled in
methodologically identical trials of the same drug.

- **Cohort 1 (early trial, naive patients).** Mean $m_{BR} = 5$,
  mean $m_{PB} = 7$, mean $m_{TV} = 1$.
- **Cohort 2 (later trial, patients who have had multiple prior
  failed treatments and are skeptical).** Mean $m_{BR} = 5$, mean
  $m_{PB} = 2$, mean $m_{TV} = 1$.

By construction the drug works equally well in both cohorts ($BR =
5$ points). But the open-label totals differ substantially: cohort
1 reports thirteen-point improvements on average, cohort 2 reports
eight-point improvements. A one-component analysis comparing the
two cohorts would conclude that the drug 'works less well' in
cohort 2 and might trigger a search for confounders, sub-population
effects, or methodological differences between the trials. None
of this is real; the difference is entirely in $PB$, and
specifically in the population baseline level of placebo
responsiveness.

The three-component decomposition, applied to either cohort,
returns $\hat{m}_{BR} = 5$ in both. The population difference is
correctly attributed to $PB$, and the apparent drug-effect
discrepancy disappears. This is the kind of reproducibility crisis
that pharmacology has spent the last decade learning to navigate;
the decomposition is one of the formal tools that makes it
tractable.

---

## Why the analysis model doesn't match the DGP

A reasonable reaction to the worked example is the following.
'If I am simulating data with three response components ($BR$,
$PB$, $TV$) plus noise, why is the analysis model written as
`Sx ~ bm + t + Dbc + bm:Dbc` with a single random intercept
and AR(1) residual correlation? Should I not include three
parametric Gompertz components in the analysis, mirroring the
DGP, so that the analyst sees the same structure that generated
the data?' The question is reasonable, and the standard answer
is that four constraints make a component-matched analysis
impractical for most aggregated N-of-1 trial sizes, and that a
modest extension (a phase-by-treatment indicator) captures the
practically important part of the gap at minimal cost.

### Four reasons not to match

**1. The trial design must supply identifiability, and many
designs supply only part of it.** The BR-PB-TV decomposition is
identified by phase contrasts: open-label versus blinded
identifies the expectancy modulation of $PB$, on-drug versus
blinded-placebo identifies $BR$, and off-drug timepoints
identify $TV$. A pure open-label trial has none of these
contrasts; in such a trial $BR$ and $PB$ both rise from zero on
treatment initiation and are statistically indistinguishable
from each other, regardless of how the analyst writes the
model. Forcing an analysis to fit three components when the
design supplies contrasts for only one or two of them does not
recover the missing components; it produces unidentified or
weakly identified parameter estimates with implausibly tight
standard errors that misrepresent the analyst's actual
inferential precision.

**2. Each component term costs degrees of freedom that come
out of the moderation estimand.** The biomarker-moderation
question is asked through the `bm:Dbc` term. A component-matched
analysis adds a Gompertz BR (three parameters: maximum, onset
rate, displacement), a Gompertz PB (three parameters plus a
phase-specific expectancy multiplier), and a Gompertz TV (three
parameters), each typically with participant-specific random
effects on the maxima. That is on the order of a dozen new
parameters, half of them on a non-linear scale. At trial
sample sizes of $N \in [30, 150]$, the variance of
$\hat{\beta}_{bm:D}$ is dominated by participant and timepoint
counts; adding parametric structure to absorb residual
confounding usually inflates the moderation estimand's
variance by more than it reduces bias. The linear-mixed model
`Sx ~ bm + t + Dbc + bm:Dbc` with `corCAR1` residual
correlation is, in this regime, the more powerful test for
the question the analyst actually asks.

**3. Component-matched analyses are non-linear in the
parameters and convergence is fragile.** Fitting Gompertz BR,
PB, TV simultaneously requires either a non-linear mixed
model (`nlme::nlme` with self-starting Gompertz forms), a
Bayesian hierarchical fit, or a structural-equation
approximation. Each is more sensitive to starting values, more
prone to local optima, and more likely to fail to converge on
a single trial dataset than the linear-mixed analogue. At
trial-relevant $N$, the non-linear component-matched fit can be
expected to fail to converge far more often than the linear
analogue, and the non-converging replicates contaminate any
aggregated analysis (argued; the convergence rate has not been
measured in this compendium, and an earlier version's figure of
'below 90%' had no source).

**4. The moderation estimand is robust to the BR-versus-PB
split for the most common scientific question.** For
predictive-biomarker validation, the regulator wants to know
whether the biomarker predicts treatment response under the
trial conditions. They typically do not require a clean
pharmacological-versus-placebo decomposition of that
prediction. The simple `bm:Dbc` interaction is the right
marginal estimand for the validation question whether the
moderation is biological or expectation-mediated. A
component-matched analysis would be required only for a
secondary question of the form 'is this biomarker predicting
the *pharmacological* response specifically, and not the
expectation response?', and that secondary question typically
admits a simpler answer through phase indicators rather than
through full component fitting.

### A small extension that usually pays off: phase indicators

The smallest useful step toward component-aware analysis,
without paying the costs above, is to add a phase indicator
and its interaction with the drug indicator. Schematically,

```{.r}
Sx ~ bm + t + Dbc + phase + Dbc:phase
       + bm:Dbc + bm:Dbc:phase
       + (1 | ptID),      correlation = corCAR1(...)
```

`phase` is a categorical covariate (open-label, blinded, post-
crossover) that absorbs the systematic differences in PB
expectancy across phases. `Dbc:phase` tests whether the
apparent drug effect shrinks under blinding, a signature
that the open-label drug effect was partly PB-mediated.
`bm:Dbc:phase` tests whether the moderation itself shrinks
under blinding, a signature that the biomarker-by-treatment
interaction was partly biomarker-by-PB rather than
biomarker-by-BR. (In the Hybrid design every open-label visit is
on drug, so `Dbc:phase` and `bm:Dbc:phase` are identified only
from the blinded phases and the open-label versus blinded on-drug
difference; argued.)

With three phase levels, this extension adds six fixed-effect
columns (two each for `phase`, `Dbc:phase` and `bm:Dbc:phase`) but
only five identified parameters, within budget at realistic $N$.
(Corrected 2026-10-08: this sentence previously counted six
fixed-effect parameters. When every open-label visit is on drug, as
in the Hybrid design, $D_{bc} = 1$ throughout the reference phase, so
$D_{bc} = 1 - \sum_{\ell \ne \mathrm{OL}} I_\ell +
\sum_{\ell \ne \mathrm{OL}} D_{bc} I_\ell$ exactly: the `Dbc` column
is an exact linear combination of the intercept, the `phase` columns
and the `Dbc:phase` columns, at every $N$ and under any coding of the
factor. The `bm:Dbc:phase` columns are not aliased. `nlme::lme` stops
on the full formula with "Singularity in backsolve" (verified), and
the Paper 06 driver, `phase_augmented_fit()` in
`analysis/scripts/component-decomposition/00-common.R`, catches the
error and refits without `Dbc:phase` (and skips the full formula for
$N < 70$), so every reported phase-augmented fit is the reduced
formula, recorded as `formula_dropped = 'Dbc:phase'` (inspected).
Verified by the rank of the model matrix in
`analysis/scripts/compendium-mathematics/check-gompertz-rmanova-cycles.R`
(three-phase Hybrid: six columns added, rank gain five). See
`analysis/report/compendium-mathematics.tex`, Section "Rank of the
phase-augmented model, and the drug-expectancy correlation by design"
(`sec:invb-collinearity`).) It does not separately
estimate $BR$ and $PB$ as components, but it does identify
the *change* in apparent drug effect that occurs when belief
is hedged. For most analyses, this is the practically
relevant question, and it can be answered without committing
to Gompertz parametric form.

### When to actually fit the full decomposition

A component-matched analysis is worth the effort when:

- Sample size is in the several hundreds (a phase 3 trial
  fielded as aggregated N-of-1, or a meta-analysis pooling
  multiple smaller trials).
- The trial design provides timepoints in every phase
  combination (open-label on-drug, open-label off-drug,
  blinded on-drug, blinded off-drug, with enough density per
  cell to constrain the Gompertz curves).
- The scientific question requires the BR-PB split (a
  pharmacology paper claiming a specific mechanism, where
  ruling out a placebo contribution is the central point).
- A Bayesian hierarchical fit with informative priors on the
  Gompertz parameters is feasible and the analyst is prepared
  to defend the priors.

In that regime, `lcmm` (for component-specific class-based
extensions), `nlme::nlme` with self-starting Gompertz forms,
or a `brms` hierarchical formulation can return
participant-specific BR, PB, and TV component estimates that
mirror the DGP. Outside that regime, the linear-mixed
analysis with phase indicators is the better tool, and the
DGP and the analysis model deliberately do not match.

The asymmetry is not a flaw of the framework; it is a
deliberate consequence of the asymmetry between simulation
(where the analyst controls the DGP and can specify any
structure, however elaborate) and inference (where the
analyst must fit a model to a finite dataset and pay the
identifiability cost of every additional parameter).
Simulation studies use the rich DGP to characterize how a
simpler analysis behaves under realistic data generation;
they do not require the analysis to match the DGP in order
to produce useful inference.

---

## The math behind it

### The core model

For participant $i$ at trial timepoint $t$:

$$
Y_{it} \;=\; \mathrm{BL}_i \;-\; [\,BR_{it} + PB_{it} + TV_{it}\,]
            \;+\; \varepsilon_{it},
$$

where the symptom score $Y_{it}$ decreases as the components grow.
Each component is a Gompertz function of its own time-on-state
variable:

- $BR_{it}$ depends on time-on-drug (cumulative exposure since
  drug initiation).
- $PB_{it}$ depends on time-on-belief, modulated by the
  phase-specific expectancy factor $\eta(\text{phase})$.
- $TV_{it}$ depends on time-since-trial-entry (calendar time).

Each component is a participant-specific random trajectory. In the
package's DGP the Gompertz parameters $(m, r, d)$ of each component
are fixed population values that set the mean trajectory; each
participant's trajectory is that mean plus a multivariate-normal
deviation with AR(1) correlation over visits (inspected:
`buildSigma()` in `R/generateData.R`). The parameters themselves
are not drawn per participant.

### Joint covariance structure

The components are not independent of each other. A patient with
strong pharmacological response is, on average, also somewhat more
likely to develop a strong placebo expectation ('this drug is
clearly working, my belief in it strengthens'), and a patient with
favorable natural history may be more likely to attribute that
improvement to whatever treatment they happen to be on. The
Hendrickson framework parameterizes this with two cross-component
correlations:

- $c_{cf1t}$: correlation between different factors at a single
  timepoint. Captures the within-time covariance of $BR$, $PB$, and
  $TV$.
- $c_{cfct}$: correlation between different factors at different
  timepoints. Captures the across-time spillover (in the package it
  is further multiplied by $\rho^{\lvert w_t - w_s\rvert}$).

Together with the within-component AR(1) correlations, these define
one joint covariance matrix per randomization path, of dimension
$(3K + 2) \times (3K + 2)$: the biomarker, the baseline level, and
$K$ visits of each of $TV$, $PB$ and $BR$, where $K$ is the number
of visits (26 variables in the eight-visit Hybrid design;
inspected: the `labels` vector of `buildSigma()`). The simulation
draws each participant's full trajectory from a single multivariate
normal of that structure.

### A common confusion: variance decomposition

A frequent slip is to write

$$
\mathrm{Var}(Y) \;=\; \mathrm{Var}(BR) + \mathrm{Var}(PB) +
\mathrm{Var}(TV) + \mathrm{Var}(\varepsilon)
$$

as if this were an alternative to working with the full covariance
matrix. By the bilinearity of covariance (the variance of a sum),
the correct identity for any joint distribution is

$$
\mathrm{Var}(BR + PB + TV + \varepsilon) \;=\;
\sum_{c} \mathrm{Var}(c) \;+\; 2 \sum_{c < c'}
\mathrm{Cov}(c, c'),
$$

with the second sum over unordered pairs of distinct components
(derived). So the simple sum-of-variances form holds when the
components are pairwise uncorrelated, and more generally whenever
the covariances sum to zero. Under the Hendrickson parameterization
with $c_{cf1t}$ and $c_{cfct}$ both non-zero, the components are
correlated, and the cross-covariance terms are non-negligible.
Working with the full covariance matrix is therefore mandatory: the
sum-of-variances shorthand is a special case that does not apply
to the parameter regimes used in this project.

The practical consequence is that any variance-component analysis
performed on simulation output should be done after the components
are extracted (using the participant-specific posterior means or
variance-component estimates from the fitted mixed model), not by
reading off the marginal variance of $Y$. The marginal variance of
$Y$ is a single number that contains no information about the
internal structure.

---

## Common questions

### Q: Could I just use a random intercept and random slope per participant?

A random intercept-and-slope mixed model is

$$
Y_{it} \;=\; \beta_0 + u_i + (\beta_t + v_i)\,t + \varepsilon_{it},
$$

with participant random intercept $u_i$ and random slope $v_i$.
This captures heterogeneity in baseline severity and in linear
time trends, but it cannot decompose the time trend into BR-style
saturating-drug, PB-style expectancy-modulated, and TV-style
natural-history pieces. It also forces each participant's trend
to be linear in $t$, which is biologically wrong for chronic-condition responses
on the timescale of weeks. The BR-PB-TV decomposition is more
specific about *why* the trajectory has the shape it does, and
that specificity translates directly into identifiability of the
clinically distinct components.

### Q: Isn't TV just disease natural history? Why not measure it independently?

In principle yes. In practice, three obstacles:

1. **Ethical.** Withholding active treatment from patients with
   active symptoms to measure natural history is rarely
   acceptable.
2. **Practical.** Natural-history studies require long observation
   windows in the absence of intervention, and they typically lack
   the within-person matching that an N-of-1 design naturally
   provides.
3. **Statistical efficiency.** When the same trial collects on-
   drug and off-drug data on the same participants, the model
   extracts $TV$ from the off-drug timepoints with high efficiency,
   and a separate natural-history study is redundant.

The decomposition is the way to recover TV from the trial itself,
without paying the cost of a parallel natural-history study.

### Q: Why is the biomarker only correlated with BR, not with PB?

This is a modeling assumption, not a fact about reality. The
Hendrickson framework correlates the baseline biomarker (pretreatment
standing systolic blood pressure in the motivating prazosin data;
Raskind et al. 2016) with the BR factor only, and treats the
biomarker as conditionally independent of PB and TV. The
justification is that the biomarker is intended to be a
*pharmacological* predictor: it measures something about the
participant's likelihood of responding to the drug's mechanism,
not their suggestibility or their natural-history trajectory.

If a candidate biomarker turned out to be correlated with PB,
that would be a substantive finding (the biomarker is partly
psychological in its action), and the model could be re-specified
to allow that correlation. The package already offers it as an
option (`c.bm.pb` in `buildSigma()`; inspected), and the stress test
of `docs/46` adds prognostic couplings to the natural course and the
baseline level as well (`docs/49`). The default specification reflects the
prior that a useful predictive biomarker should be a clean
pharmacological predictor; deviations from that prior are testable.

### Q: What if TV is negative (the patient naturally worsens)?

In the standard parameterization with positive $m_{TV}$, the
Gompertz form is monotonically increasing and TV is non-negative
by construction. To allow naturally-worsening trajectories, the
maximum $m_{TV}$ can be drawn from a population distribution that
includes negative values; the resulting Gompertz curve is then
monotonically decreasing in symptom-reduction (equivalently,
monotonically increasing in symptom severity). The model handles
this case without modification, and detecting $\hat{m}_{TV} < 0$
is itself informative: it implies the underlying condition is
deteriorating, and the apparent treatment effect is partly
camouflaging deterioration that would otherwise be visible.

### Q: How many timepoints do I need to identify all three components?

Roughly, at least one timepoint per phase contrast that you want
to exploit. A pure open-label design with eight weekly timepoints
can identify $BR + PB$ jointly but cannot separate them. Adding a
blinded discontinuation phase with at least three or four
timepoints lets the model identify $BR$ separately. Adding a
longer off-drug follow-up extends identifiability of $TV$,
particularly if the participant's natural-history trajectory is
still evolving. The design contrasts table earlier in this
document is the formal version of the answer: each component
needs at least one phase in which it is uniquely present (or
uniquely absent in a comparable arm), and the more such phases the
better.

### Q: What if my drug has a fast onset and a long washout?

The Gompertz parameters absorb this asymmetry. A drug with a fast
$BR$ onset has a high $r$; a drug with a long pharmacological
washout has a long off-drug tail captured by the carryover
half-life $t_{1/2}$ used in the analysis-side exposure-decayed drug
indicator, $D_{bc,it} = 1$ on drug and $e^{-\lambda t_{sd}}$ off
drug after exposure (0 before the first exposure; `docs/45`,
Section 8). The decomposition is
agnostic to the specific pharmacokinetics so long as the trial
includes enough timepoints in the relevant transition windows
(on-to-off and off-to-on) to constrain the curves.

---

## Why this matters for biomarker validation

The motivating clinical question for the trial designs in this
project is: 'does the biomarker predict treatment response?', and
specifically 'does the biomarker predict the *pharmacological*
component of treatment response?'. This is the precision-medicine
question, and it depends critically on the BR-PB-TV decomposition.

### The wrong way

A direct regression of total response on baseline biomarker

```
total_response = beta_0 + beta_bm * biomarker + error
```

estimates the biomarker's correlation with whatever the lumped
total contains. If $PB$ varies systematically across the biomarker
distribution (because, say, the biomarker is correlated with
optimism, suggestibility, or treatment-seeking history), then
$\hat{\beta}_{bm}^{\text{lumped}}$ (the coefficient `beta_bm` above)
contains a placebo-prediction component that is
inseparable from the pharmacological-prediction component. The
resulting biomarker, applied prospectively, would mis-stratify
patients: those flagged as 'high responders' by the biomarker
would include both true high-BR responders and high-PB responders,
and the latter would not benefit from prescription decisions made
on a pharmacological basis.

### The right way

Regress *each component* on the biomarker:

```
BR    = beta_0_BR  + beta_bm_BR  * biomarker + error_BR
PB    = beta_0_PB  + beta_bm_PB  * biomarker + error_PB
TV    = beta_0_TV  + beta_bm_TV  * biomarker + error_TV
```

A biomarker is a useful pharmacological predictor if and only if
$\hat{\beta}_{bm}^{BR}$ is meaningfully non-zero. Significant
$\hat{\beta}_{bm}^{PB}$ or $\hat{\beta}_{bm}^{TV}$ are diagnostic
findings about the biomarker (it is partly psychological, or it
correlates with natural-history characteristics), but they are
not the validation of pharmacological prediction. A biomarker that
correlates only with $PB$ would be valuable for selecting placebo
responders in run-in designs but worthless for predicting
pharmacological efficacy; a biomarker that correlates only with
$TV$ would be a prognostic marker, not a predictive one.

The decomposition is therefore not optional for the precision-
medicine question: it is the formal mechanism by which the question
is asked at all.

---

## Summary table

| Component | Cause | Timescale | Identified by | Inferential role |
|---|---|---|---|---|
| **BR** | Drug pharmacology | Pharmacokinetic (hours-weeks) | Phase contrasts (on-drug vs. blinded placebo) | Pharmacological efficacy, biomarker-prediction target |
| **PB** | Belief in treatment | Expectancy-modulated (hours-weeks) | Open-label vs. blinded contrasts | Placebo response, accounted for, not the target |
| **TV** | Natural history | Disease-trajectory (weeks-months) | Off-drug timepoints | Confounding control |
| **$\varepsilon$** | Measurement noise | Per-observation | Residuals | Fit diagnostics |

Table: Summary of the BR, PB, TV, and noise components and their roles

---

## Key takeaway

The BR-PB-TV decomposition is not technical bookkeeping. It is the
mechanism by which the trial design and the model together
identify clinically distinct causes of symptom change from a single
within-person trajectory. Without the decomposition, the analyst
has only the lumped total, and the lumped total mixes
pharmacological efficacy with placebo response and natural
history. Pharmacological efficacy is the regulatory and clinical
quantity of interest; the other two are biases that must be
controlled for before efficacy can be inferred.

The trial design (open-label phase, blinded discontinuation,
crossover) supplies the identifiability. The model extracts the
components from the design contrasts. The biomarker validation
question, if asked at the level of total response rather than the
BR component, is unanswerable.

---

## Further reading

- Hendrickson, Thomas, Schork, and Raskind (2020). *Optimizing
  aggregated N-of-1 trial designs for predictive biomarker
  validation: statistical methods and theoretical findings.*
  Frontiers in Digital Health 2020;2:13.
  doi:10.3389/fdgth.2020.00013. PMID 34713026. The original paper
  that introduced this decomposition framework for N-of-1 trials
  targeting predictive biomarker validation.
- Senn (2018). *Statistical pitfalls of personalized medicine.*
  Nature 563(7733):619-621. doi:10.1038/d41586-018-07535-2.
  PMID 30482931. On variance decomposition and per-participant
  inference in personalized-medicine settings.
- Benedetti (2014). *Placebo effects: from the neurobiological
  paradigm to translational implications.* Neuron 84(3):623-637.
  doi:10.1016/j.neuron.2014.10.023. PMID 25442940. On the
  neurobiology of the PB component.
- Hobbs DC, Twomey TM, Palmer RF (1978). *Pharmacokinetics of
  prazosin in man.* J Clin Pharmacol 18(8-9):402-406.
  doi:10.1002/j.1552-4604.1978.tb02456.x. PMID 690251. Source of
  the 2.3-hour plasma half-life.
- Raskind MA, Millard SP, Petrie EC, et al. (2016). *Higher
  pretreatment blood pressure is associated with greater
  posttraumatic stress disorder symptom reduction in soldiers
  treated with prazosin.* Biol Psychiatry 80(10):736-742.
  doi:10.1016/j.biopsych.2016.03.2108. PMID 27320368. The motivating
  biomarker.
- Barnett AG, van der Pols JC, Dobson AJ (2005). *Regression to the
  mean: what it is and how to deal with it.* Int J Epidemiol
  34(1):215-220. doi:10.1093/ije/dyh299. PMID 15333621.

Each entry above was checked against PubMed on 2026-10-07.

- Project documentation:
  `docs/06-ar1-residual-correlation.tex` for the within-factor
  serial-correlation structure;
  `docs/08-biomarker-correlation-decay.tex` for biomarker-BR
  correlation decay during off-drug periods;
  `docs/09a-carryover-bias-variance.md` for the bias-variance
  trade in carryover specification;
  `docs/19-mean-moderation-implementation-notes.tex` for the
  mean-moderation architecture's implementation;
  `docs/40-hendrickson-hybrid-design-primer.md` for the Hybrid
  design, the published parameters, and the corrections applied to
  this guide in 2026-10;
  `analysis/report/01-dgp-mean-moderation-vs-mvn/` for the manuscript
  treatment of the BR moderation question.

---

*Audience: graduate students in biostatistics, clinical research
methodologists, and trialists encountering the BR-PB-TV
decomposition for the first time. Familiarity with linear mixed
models and Gompertz curves is assumed; deep knowledge of
psychometric latent-variable methods is not.*
