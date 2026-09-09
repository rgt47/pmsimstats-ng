# The BR-PB-TV Decomposition, Its Identification, and Its Reductions

*2026-09-09 09:26 PDT*

Author: pmsimstats team

## 0. What this document is

This white paper consolidates two earlier documents and extends them
in two directions.

Document 24, the component decomposition pedagogy, supplies Parts I
through III. It develops the three-component BR-PB-TV decomposition
with worked examples and an intuition narrative. Document 24 remains
in place as a standalone introduction and is not superseded. A reader
meeting the decomposition for the first time may prefer it, since it
is shorter and written to be read straight through. Parts I through
III reproduce its content so that the present document is complete on
its own.

Document 28, on the Hendrickson approach and the hybrid design, is
folded into Parts II and III and is superseded.

Two parts are new. Part IV places the decomposition in the context of
the clinical trials and N-of-1 literatures, and asks which of its
three components are conventional and which are innovations of this
program. Part V reports papers 13 and 14, which evaluate the two
available two-component reductions.

We should state the organizing observation at the outset, since it
runs through both new parts. The three components are not equally
well established. A trend term is standard in the crossover and
N-of-1 literature. A modeled placebo-belief trajectory is not, and
the empirical placebo literature regards belief and natural history
as separable only under designs that most trials do not field. That
asymmetry turns out to govern which reductions of the decomposition
are defensible.

### Reading order

| Part | Content | Source |
|---|---|---|
| I | The three components and why they are separated | doc 24 |
| II | Identification from trial design | docs 24, 28 |
| III | Analysis models, and why they need not match the DGP | docs 24, 28 |
| IV | **The decomposition in the wider literature** | **new** |
| V | **The two-component reductions** | **papers 13, 14 (new)** |
| VI | Reference: questions, summary table, bibliography | doc 24 |

A reader familiar with the decomposition may begin at Part IV. A
reader who wants only the practical guidance may read section 21,
which supersedes the earlier guidance of sections 12 and 13.

### A note on maintenance

Parts I through III duplicate document 24, which is maintained
separately. Document 24 is canonical for that material. A correction
to the components pedagogy should be made there first and mirrored
here.

# Part I. The three components

## 1. Why split the response?

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
   your belief activates well-characterised neurobiological
   pathways: dopaminergic reward signalling, descending pain
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
generalise: a population with stronger placebo expectations would
appear to have a more powerful drug, and a population with milder
natural-history drift would appear to have a less powerful one.
Neither inference would be correct.

Decomposing the response separates these failure modes.

---


## 2. Biological response (BR)

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
relevant doses it antagonises noradrenergic activity in the locus
coeruleus and prefrontal cortex during REM sleep. The mechanism for
nightmare reduction is reasonably specific: prazosin damps the
adrenergic surge that accompanies trauma-related dream content,
and the resulting suppression of REM-stage arousal reduces
nightmare frequency and severity.

The pharmacokinetic profile is short. Prazosin has a half-life of
two to three hours; brain concentrations equilibrate quickly with
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
BR(t) \;=\; m \,\exp\bigl(-d \,\exp(-r t)\bigr) - m\exp(-d),
$$

rescaled so that $BR(0) = 0$ and $BR(t) \to m$ as $t \to \infty$.
The three parameters have direct interpretations:

| Parameter | Symbol | Meaning |
|---|---|---|
| Maximum response | $m$ | The biological ceiling. The largest reduction the drug can produce, even with infinite time on treatment. |
| Onset rate | $r$ | How fast the patient approaches the ceiling. High $r$ means a fast onset (effect plateaus in one to two weeks); low $r$ means a slow climb (six to eight weeks to plateau). |
| Displacement | $d$ | The shape of the climb. High $d$ is an early-onset curve that rises sharply at first; low $d$ is a more gradual climb. |

For PTSD prazosin, calibrated values from the Hendrickson et al.
(2020) reference dataset are roughly $m \approx 11$ (eleven points
of nightmare score reduction at saturation), $r \approx 0.42$ per
week (half-maximum at about eighteen days), $d \approx 5$
(moderately steep early climb).

The Gompertz family was chosen because it captures the
'asymptote-from-below' shape that physiological dose-response
curves typically display. Two alternatives that do **not** capture
this shape:

- A linear function in $t$ ($BR = \beta_t \cdot t$) climbs without
  bound and is biologically implausible past a few weeks.
- An exponential decay ($BR = m(1 - \exp(-rt))$) is acceptable for
  fast-onset drugs but lacks the slow-start phase that characterises
  most chronic-condition responses, where receptor adaptation,
  feedback regulation, and downstream gene-expression changes
  produce a measurable lag.

The Gompertz interpolates between exponential decay (in the
high-$d$ limit) and a linear increase (in the low-$d$ limit), and
its three-parameter form is flexible enough to fit a wide range of
real pharmacodynamic profiles without overfitting.

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


## 3. Placebo-belief response (PB)

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
characterised neurobiological pathways. Belief-driven release of
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
phase of the Hendrickson hybrid design is the canonical such phase.

### Real example: the blinded discontinuation contrast

Consider a sixteen-week N-of-1 trial structured in three phases:

- **Weeks 1-8 (open-label on drug).** The patient knows they are
  receiving prazosin. Belief in treatment is at full strength
  (expectancy multiplier 1.0, in the model parameterisation
  introduced below).
- **Weeks 9-12 (blinded discontinuation).** The patient is told
  that, at some unannounced point during this window, they may be
  silently switched to placebo. They cannot reliably tell whether
  they are still receiving the active drug. Belief in treatment
  drops to a partial level (multiplier roughly 0.5) because the
  patient's expectation is now hedged.
- **Weeks 13-16 (open-label crossover, possibly back on drug).**
  The patient is told whether they are now on active drug or
  placebo, and belief returns to a level commensurate with that
  knowledge.

A patient who is silently switched to placebo at the start of week
9 will have:

- $BR = 0$ from week 9 onward (no pharmacology).
- $PB \approx 0.5 \cdot PB_{\text{open-label}}$ during weeks 9-12,
  reflecting hedged belief.
- $TV$ continuing on its slow trajectory.

A patient who remains on active drug during weeks 9-12 will have:

- $BR$ continuing at its plateau value.
- The same hedged $PB \approx 0.5 \cdot PB_{\text{open-label}}$ as
  the placebo-switched patient.
- The same $TV$ trajectory.

The contrast between the two arms during weeks 9-12 is therefore a
nearly-clean estimate of $BR$, with $PB$ matched (because both
groups have the same belief) and $TV$ matched (because $TV$ is by
construction independent of treatment). This is the inferential
move that makes the BR-PB-TV decomposition identifiable from the
data: the trial design supplies phases in which different
combinations of components are present, and the analyst recovers
the components by differencing those phases.

### Mathematical pattern

PB also follows a Gompertz curve, but scaled by an
expectancy factor $\eta(\text{phase})$ that depends on the trial
phase the patient is currently in:

$$
PB(t \mid \text{phase}) \;=\; \eta(\text{phase}) \,\cdot\,
m_{PB} \,\exp\bigl(-d_{PB} \exp(-r_{PB} t)\bigr) \,+\,
\text{constant adjustment to start at zero.}
$$

Typical values of the expectancy factor:

| Phase | $\eta$ | Reasoning |
|---|---|---|
| Open-label on drug | 1.0 | Patient knows they are on the active drug. |
| Blinded discontinuation | 0.5 | Patient believes they may or may not still be on drug; expectation is hedged. |
| Open-label placebo | 0.0-0.5 | Patient told they are on placebo retains some residual expectation, especially if the previous active phase was beneficial. |

The choice of $\eta = 0.5$ rather than $\eta = 0$ for the blinded
phase is itself a modelling assumption. The empirical question is
how much expectation a patient carries when uncertain about
allocation, and the published placebo literature suggests it is
substantial: residual expectations of 30-60% of the open-label
level are reported across pain, depression, and PTSD trials. Setting
$\eta = 0.5$ in the blinded phase is a reasonable default, but
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
variance shrinks accordingly. The model parameterisation reflects
this with a phase-dependent standard deviation:

$$
\mathrm{SD}(PB \mid \text{phase}) \;=\;
\eta(\text{phase}) \cdot \sigma_{PB}^{\text{baseline}}.
$$

This heteroscedastic structure is not optional: ignoring it leads
to inflated standard errors in the open-label phase (where $PB$
contributes most variance) and deflated standard errors in the
blinded phase, with consequences for both confidence intervals on
the drug effect and the calibration of any subsequent biomarker
test.

---


## 4. Natural history (TV)

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
a post-operative pain syndrome, an episode of acute back pain --
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
Patients enrol at a moment of clinical engagement, which often
coincides with a recent worsening that triggered the
help-seeking behaviour. Regression to the mean alone -- the
statistical fact that an extreme observation tends to be followed
by less extreme observations on average -- can produce several
points of apparent TV-improvement in the early weeks of a trial,
even in a population whose underlying condition is stable. This
is the famous failure mode of single-arm trials in chronic-pain
conditions: large 'treatment' effects evaporate when a control arm
is added, because the apparent effect was almost entirely
regression to the mean.

**Trial participation itself induces TV.** The Hawthorne effect
(participants improve simply because they are being observed),
the structured-visit effect (regular visits with clinicians may
improve symptom management directly), and the
diary-completion effect (the act of monitoring symptoms changes
how they are perceived and reported) are all real, all
participant-mediated, and all confounded with calendar time. They
appear in the data as TV-improvement that is independent of any
pharmacological or placebo response, and they apply equally to
all arms of the trial. Modelling them as part of TV is the
appropriate handling.

### Can we assume no TV effect?

The short answer is: rarely safely, and never without empirical
support. The longer answer:

- **TV may be assumed zero only if** the disease is in a stable
  chronic phase, the trial is short relative to the disease
  timescale, the recruitment process does not favour participants
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
common single source of bias in chronic-condition trial inference.
A famous example: many early trials of chronic-fatigue therapies
appeared to show large effects, until a careful natural-history
arm revealed that the apparent effect was almost entirely a
regression-to-the-mean phenomenon driven by the selection of
patients at peak symptom severity. The patients would have improved
similarly with no intervention. TV in the model is the formal way
to handle this risk.

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

TV also follows a Gompertz curve but parameterised independently of
BR and PB:

$$
TV(t) \;=\; m_{TV} \,\exp\bigl(-d_{TV} \exp(-r_{TV} t)\bigr) \,+\,
\text{offset to start at zero.}
$$

The maximum $m_{TV}$ can be either positive (natural improvement)
or negative (natural worsening), and is a participant-level random
effect drawn from the population distribution of natural-history
trajectories. A trial design that aggregates over many participants
can identify both the population mean and the population variance
of $m_{TV}$, but only if there are enough off-drug timepoints to
constrain the curve.

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


---

# Part II. Identification from trial design

## 5. How the design identifies the components

The mathematical decomposition above is meaningless unless the
data contain enough information to pin down each component
separately. Trial design is what supplies that information. The
Hendrickson hybrid design, used as the worked example throughout
the rest of this guide, has three phases that each contribute a
different combination of components.

| Phase | Drug? | Belief? | Components present |
|---|---|---|---|
| Open-label on drug | Yes | Full | $BR + PB(\eta = 1) + TV$ |
| Blinded discontinuation | Hidden | Hedged | $BR \cdot \mathbf{1}_{\text{on drug}} + PB(\eta = 0.5) + TV$ |
| Open-label crossover | Either | Knows | $BR \cdot \mathbf{1}_{\text{on drug}} + PB(\eta = 1) + TV$ |

Three phase-by-allocation contrasts then identify the three
components:

- **Phase 1 vs. phase 2 (within drug):** $PB(\eta = 1) - PB(\eta =
  0.5) \approx 0.5 \cdot PB$. Isolates the belief component.
- **Phase 2 on drug vs. phase 2 placebo:** $BR$ alone (because $PB$
  and $TV$ match across allocation arms within phase 2).
- **All phases off drug:** $TV$ alone (no $BR$, no $PB$ at the
  off-drug expectation level).

Each contrast targets a different component, and the analyst
combines them through the linear-mixed-effects model fitted to the
full trajectory. The model is the formal mechanism, but the
inference is driven by the design contrasts.

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


## 6. The Hendrickson approach

The founding pmsimstats question is the detection of a
biomarker-by-treatment interaction in an aggregated N-of-1 trial. The
Hendrickson reference implementation
(`implementations/original/`) is the historical anchor for the whole
programme. Its defining feature, and the starting point for paper 06,
is an asymmetry between how it generates data and how it analyses it.

### 2.1 The data-generating process is fully decomposed

The DGP in `implementations/original/R/generateData.R` builds the
observed symptom score as the sum of three causally distinct
components, each a modified Gompertz trajectory with
participant-specific random effects:

- `BR` (biological response, the pharmacological component), driven
  by cumulative time on drug, decaying off drug with a carryover
  half-life;
- `PB` (placebo-belief response), driven by cumulative time on
  belief and scaled by the phase expectancy multiplier eta;
- `TV` (natural-history component), driven by time
  since trial entry.

The components are summed and subtracted from baseline (lines near
99 and 104 of `generateData.R`):

```
D_t  = rowSums(t.br + t.pb + t.tv)   # total response at timepoint t
Sx_t = BL - D_t                       # symptom = baseline minus total
```

The raw simulated data retain the individual component values
(`t1.br`, `t1.pb`, `t1.tv`, and so on for each timepoint), so the
ground truth is three-component by construction. The prazosin-PTSD
calibration used throughout the programme draws its Gompertz
parameters from the Hendrickson et al. (2020) dataset
(`m_BR` approximately 11, `r_BR` approximately 0.42 per week,
`d_BR` approximately 5).

### 2.2 The analysis model is a single lumped drug indicator

The analysis in `implementations/original/R/lme_analysis.R` does not
mirror the DGP. It fits one continuous drug indicator, not three
components:

```
Sx ~ bm + t + Dbc + bm:Dbc
random      = ~ 1 | ptID
correlation = corCAR1(form = ~ t | ptID)
```

where `Dbc` is the continuous drug state (equal to 1 on drug and to
`exp(-lambda * tsd)` off drug, with `lambda = ln 2 / t_half`), and
`bm:Dbc` is the biomarker-by-drug interaction. Variant forms in the
same file add an expectancy indicator (`De`) or a carryover term
(`tsd`), but they remain single-drug-indicator models. The model
estimates one lumped drug main effect and one lumped biomarker
slope; it does not recover `BR`, `PB`, and `TV` separately.

### 2.3 The asymmetry that paper 06 studies

The original therefore applies decomposition only as a
data-generating device. It knows the response is `BR + PB + TV` and
simulates accordingly, but its inferential model collapses the three
into a single treatment coefficient. This lumped estimator is exactly
the baseline that paper 06 formalises, and the central question of
that paper is the consequence of the asymmetry: when does the lumped
analysis of three-component data give a biased or uninformative
answer, and when does adding analysis-side decomposition fix it.


## 7. The hybrid design

Component decomposition is meaningless unless the data contain enough
information to pin down each component separately. That information is
supplied by trial design, not by the analysis model. The hybrid
open-label-plus-blinded-discontinuation design is the worked example
throughout paper 06 and the primary identifiability vehicle.

### 7.1 Phase structure

The hybrid design strings together three phases, each contributing a
different mixture of the three components:

| Phase | Drug | Belief | Components present |
|---|---|---|---|
| Open-label on drug | Yes | Full (eta = 1) | BR + PB(eta=1) + TV |
| Blinded discontinuation | Hidden | Hedged (eta = 0.5) | BR (on drug) + PB(eta=0.5) + TV |
| Open-label crossover | Either | Knows (eta = 1) | BR (on drug) + PB(eta=1) + TV |

The `on drug` qualifier on `BR` is shorthand for the continuous-decay
drug state `Dbc`, not a binary switch.

### 7.2 The three identifying contrasts

Three phase-by-allocation contrasts identify the three components:

- **Open-label versus blinded (within drug)** yields
  `PB(eta=1) - PB(eta=0.5)` which is approximately `0.5 * PB`, and so
  isolates the belief component, conditional on a known blinded-phase
  expectancy.
- **Blinded on-drug versus blinded placebo** yields `BR` once drug
  carryover has elapsed, because `PB` and `TV` match across allocation
  arms within the blinded phase.
- **Off-drug phases** yield `TV` plus any residual off-drug `PB`, with
  no active `BR` contribution beyond carryover.

The lumped LME combines these into one model, but the inference is
driven by the design. The one-component analysis already exploits the
drug-on/off contrast through its `Dbc` term; that is why, as section 11
shows, it is not as impoverished as its single coefficient suggests.

### 7.3 The expectancy identification limit

The hybrid design has one principled limitation that deserves
emphasis. With only an open-label and a blinded phase, the data
identify the product `eta * m_PB`, not the two factors separately. The
working value $\eta \approx 0.5$ in the blinded phase is a
modelling assumption, not an estimate, and absent a belief measurement
it is not falsifiable from outcome data alone. The identified set for
$m_{PB}$ is the interval obtained by dividing the estimated product by
$\eta$ over its plausible range, collapsing to a point only when $\eta$
is pinned down. Two design features close the gap:

- an open-label placebo phase, which supplies a second equation that
  separates $\eta$ from $m_{PB}$;
- direct measurement of treatment expectancy, as in the open-hidden
  and balanced-placebo designs that manipulate or record belief.

Any trial relying on the decomposition for a PB-versus-BR attribution
should include such a phase or collect an expectancy measure, and
should report sensitivity of the attribution to $\eta$ across its
plausible range.

### 7.4 Identifiability is necessary but not sufficient

Design-contrast identifiability does not guarantee stable estimation
at finite sample sizes. The three Gompertz trajectories share an onset
shape and overlap in calendar time, so their parameters can be
near-collinear even when the identifying contrasts are present. Joint
recovery additionally requires a phase in which each component varies
while the others are held fixed, and an off-drug window long enough,
relative to the carryover half-life, for `TV` to separate from a still
decaying `BR`. Where these conditions are weak the components are
identified only up to a near-collinear direction; the symptom is the
rank deficiency the Study B pilot encountered when the full
phase-augmented formula was fitted at `N = 35`. Recoverability is to be
mapped by identifiability diagnostics (convergence rates,
variance-inflation factors on the component-specific parameters, and
sensitivity of the `BR` estimate to the assumed $\eta$), not asserted.


---

# Part III. Analysis models

## 8. Why the analysis model does not match the DGP

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
trial-relevant N, the convergence rate of the non-linear
component-matched fit is typically below 90%, and the
non-converging cells contaminate any aggregated analysis.

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

```r
Sx ~ bm + t + Dbc + phase + Dbc:phase
       + bm:Dbc + bm:Dbc:phase
       + (1 | ptID),      correlation = corCAR1(...)
```

`phase` is a categorical covariate (open-label, blinded, post-
crossover) that absorbs the systematic differences in PB
expectancy across phases. `Dbc:phase` tests whether the
apparent drug effect shrinks under blinding -- a signature
that the open-label drug effect was partly PB-mediated.
`bm:Dbc:phase` tests whether the moderation itself shrinks
under blinding -- a signature that the biomarker-by-treatment
interaction was partly biomarker-by-PB rather than
biomarker-by-BR.

This extension adds three to four parameters at most, well
within budget at any realistic N. It does not separately
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
Simulation studies use the rich DGP to characterise how a
simpler analysis behaves under realistic data generation;
they do not require the analysis to match the DGP in order
to produce useful inference.

---


## 9. The mathematics of the decomposition

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

Each component is a participant-specific random function: the
Gompertz parameters $(m, r, d)$ for BR, PB, TV are drawn from
population distributions, and the participant-level deviations are
the random effects in the mixed model.

### Joint covariance structure

The components are not independent of each other. A patient with
strong pharmacological response is, on average, also somewhat more
likely to develop a strong placebo expectation ('this drug is
clearly working, my belief in it strengthens'), and a patient with
favourable natural history may be more likely to attribute that
improvement to whatever treatment they happen to be on. The
Hendrickson framework parameterises this with two cross-component
correlations:

- $c_{cf1t}$: correlation between different factors at a single
  timepoint. Captures the within-time covariance of $BR$, $PB$, and
  $TV$.
- $c_{cfct}$: correlation between different factors at different
  timepoints. Captures the across-time spillover.

Together these define a $4 \times K \times K$ block covariance
matrix (one $K \times K$ block per pair of factors, where $K$ is
the number of timepoints), and the simulation draws each
participant's full trajectory from a single multivariate normal of
that structure.

### A common confusion: variance decomposition

A frequent slip is to write

$$
\mathrm{Var}(Y) \;=\; \mathrm{Var}(BR) + \mathrm{Var}(PB) +
\mathrm{Var}(TV) + \mathrm{Var}(\varepsilon)
$$

as if this were an alternative to working with the full covariance
matrix. By the law of total variance for sums, the correct
identity for any joint distribution is

$$
\mathrm{Var}(BR + PB + TV + \varepsilon) \;=\;
\sum_{c} \mathrm{Var}(c) \;+\; 2 \sum_{c \neq c'}
\mathrm{Cov}(c, c').
$$

So the simple sum-of-variances form holds **if and only if** the
components are uncorrelated. Under the Hendrickson parameterisation
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


## 10. Worked example: what is lost without decomposition

Consider two participants, A and B, enrolled in a sixteen-week
N-of-1 trial of prazosin for PTSD. Both are otherwise similar:
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
strength, and TV accumulates linearly. By week 8, both participants
have reached close to their saturation values:

|  | $BR$ | $PB(\eta = 1)$ | $TV$ | Total improvement |
|---|---|---|---|---|
| Participant A | 4.0 | 6.0 | 1.0 | **11.0** |
| Participant B | 8.0 | 1.0 | 1.0 | **10.0** |

A clinician examining only the totals would conclude that the
two participants are responding similarly: both improved by about
ten points, both look like clinical successes. A simple t-test on
the open-label change would estimate a population drug effect of
about ten points and would not distinguish the two participants.

Now consider the blinded discontinuation phase (weeks 9-12). Both
participants are silently switched to placebo at the start of week
9. The expectancy multiplier drops from 1.0 to 0.5. By the end of
week 12, the components are:

|  | $BR$ | $PB(\eta = 0.5)$ | $TV$ | Total improvement |
|---|---|---|---|---|
| Participant A | 0.0 | 3.0 | 1.0 | **4.0** |
| Participant B | 0.0 | 0.5 | 1.0 | **1.5** |

The two participants now look quite different. Participant A still
shows a sizeable improvement (four points) under blinded placebo;
participant B has nearly returned to baseline (one and a half
points). This is the diagnostic contrast that the decomposition
exploits.

What does each modelling strategy infer?

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
After fitting, the model returns:

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
| 'Does the drug work for A?' | 'Some response, magnitude unclear, ~7.75 points overall.' | 'Yes: $BR_A = 4.0$ points, distinct from her larger placebo response.' |
| 'Does the drug work for B?' | 'Some response, magnitude unclear, ~7.75 points overall.' | 'Yes: $BR_B = 8.0$ points, twice the magnitude of A.' |
| 'Are A and B different responders?' | 'No discernible difference; both improved by ~10 points.' | 'Yes: B has 2x A's BR but 1/6 of A's PB.' |
| 'If a biomarker predicts BR, is it useful here?' | Cannot answer (no $BR$ estimate). | Yes if the biomarker correlates with $\hat{m}_{BR}$. |
| 'Should A be continued on this drug long-term?' | Cannot say (the response could be largely placebo). | Probably yes (BR = 4 is real and pharmacological), but the larger placebo component will not persist if expectations decay. |
| 'What is the population mean drug effect?' | 7.75 points (biased; mixes BR and the discontinuation-induced drop in PB). | 6.0 points = $(4.0 + 8.0)/2$ (correctly the mean of $BR$). |

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


## 11. The paper 06 analysis

### 11.1 The omitted-variable-bias identity

The analytic core of paper 06 is an exact identity for what the
lumped estimator estimates. By linearity of covariance, the lumped
biomarker-by-treatment slope decomposes as

```
beta_bm^lumped = beta_bm^BR + w_PB * beta_bm^PB + w_TV * beta_bm^TV
```

so the lumped slope is displaced from the pharmacological slope
`beta_bm^BR` by a term proportional to the biomarker's covariance with
`PB` and `TV`, weighted by design-dependent weights `w_PB` and `w_TV`.
The displacement is in the probability limit, not the variance: a
larger one-component trial converges on the same biased value. The
bias is zero if and only if the biomarker is orthogonal to `PB` and
`TV`, regardless of the magnitude of those components. The companion
identity for the lumped treatment effect is displaced by the `PB` and
`TV` contributions themselves, so a lumped pre-post treatment effect
is biased whenever `PB` or `TV` is non-zero.

### 11.2 Three analysis strategies compared

Paper 06 compares three analysis strategies against three-component
data:

1. **One-component (lumped).** The Hendrickson model,
   `Sx ~ bm + t + Dbc + bm:Dbc`.
2. **Phase-augmented LME.** A modest extension,
   `Sx ~ bm + t + Dbc + phase + Dbc:phase + bm:Dbc + bm:Dbc:phase`,
   which tests whether the drug and biomarker-by-drug effects attenuate
   under blinding, without the full Gompertz parametric form.
3. **Full component-matched decomposition.** A parametric model that
   mirrors the DGP and recovers component-specific slopes, feasible only
   inside the recoverable region of section 7.4.

### 11.3 What the simulations show

**Study A (orthogonal-biomarker arm, hybrid design; 1,000 replicates
per alternative cell, 5,000 per null).** With the biomarker coupled to
`BR` only (`c_bm,PB = 0`), the one-component estimate of the
biomarker-by-treatment slope is essentially unbiased across the entire
`(m_PB, m_TV)` grid: mean absolute bias 0.016 at `N = 35` and 0.006 at
`N = 150` against a true value of -2.25, both within one Monte Carlo
standard error of zero, with near-nominal type I error (0.031 to 0.056
for the one-component analysis). The phase-augmented analysis reduces
no bias (0.028 versus 0.016 at `N = 35`) and is markedly less powerful
(mean power 0.59 versus 0.35 at `N = 35`; 0.90 versus 0.62 at
`N = 70`; both reach approximately 1.0 by `N = 100`). Under an
uncontaminated biomarker, decomposition is therefore strictly
dominated.

**Study B (coupling-arm pilot; 100 replicates per cell).** With the
biomarker coupled to `PB`, the one-component bias scales with the
coupling strength $c_{bm,PB}$, not with the magnitude of `PB`, exactly as
the identity predicts (the one-component bias is approximately 0.51 at
the contamination level studied). The phase-augmented analysis offers
no protection: its bias is statistically indistinguishable from the
one-component bias at every contamination level.

**Study C (recovery under a belief-decoupling design; 1,000 replicates
per cell, $N$ in {70, 150}, $c_{bm,PB}$ in {0, 0.1, 0.2, 0.3}).** This is
the decisive design comparison. On the standard hybrid design the
drug-state and belief contrasts are collinear (analysis-stage
correlation approximately +0.6), and the belief-covariate,
blinded-stratum, and Gompertz-basis decompositions all fail to recover
`beta_bm^BR`, no better than the one-component analysis. A
balanced-placebo (open-hidden) design adds a covert on-drug phase
(`eta = 0`) and an open-placebo phase (`eta = 1`), reducing the
drug-state-by-belief correlation to approximately -0.66 and breaking
the collinearity. On that design the belief-covariate decomposition
returns the unbiased pharmacological slope across the feasible
contamination range (bias -0.022, +0.008, +0.021, -0.026 at `N = 150`,
MCSE approximately 0.02; coverage 0.94 to 0.96), whereas the
one-component estimator degrades monotonically (bias +0.16, +0.34,
+0.48; coverage down to 0.86). The grid is capped at `c_bm,PB = 0.30`
because beyond approximately 0.45 the implied covariance is
non-positive-definite.

The power comparison between the two estimators on the balanced-placebo
design is more nuanced than a uniform decomposition penalty (rejection
rate for the interaction at alpha = 0.05, computed from the
1,000-replicate per-cell p-values; paired standard error of the gap at
most 0.014):

| c_bm,PB | N | One-component | Decomposition | Gap |
|---|---|---|---|---|
| 0.0 | 70 | 0.85 | 0.66 | +0.19 |
| 0.1 | 70 | 0.78 | 0.67 | +0.11 |
| 0.2 | 70 | 0.74 | 0.69 | +0.05 |
| 0.3 | 70 | 0.64 | 0.67 | -0.03 |
| 0.0 | 150 | 0.99 | 0.93 | +0.07 |
| 0.1 | 150 | 0.98 | 0.93 | +0.05 |
| 0.2 | 150 | 0.97 | 0.93 | +0.04 |
| 0.3 | 150 | 0.93 | 0.93 | -0.00 |

At zero contamination the one-component analysis is the more powerful of
the two (by 0.19 at N = 70, 0.07 at N = 150), reflecting the degrees of
freedom the decomposition spends. That advantage erodes monotonically as
contamination grows and reverses by `c_bm,PB = 0.3` (0.64 versus 0.67 at
N = 70), because contamination in this direction drags the biased
one-component estimate toward the null and so costs it power, while the
decomposition holds the estimate at the pharmacological target. The
one-component power advantage is therefore confined to the
low-contamination regime in which decomposition is unnecessary. (The sign
of the gap depends on the direction of the coupling: a coupling that
reinforced the BR sign would inflate the one-component estimate and give
it spuriously higher power, which is power to detect a contaminated
quantity rather than a genuine advantage.)

On the hybrid design the same comparison is uninformative because power
saturates near 1.0 for both estimators (the BR effect of -2.25 is large
relative to its standard error at N = 70 to 150), and both are biased
under contamination in any case; there the operative axis is bias and
coverage, not power.


## 12. When decomposition pays off

*The guidance of this section and of section 13 was written before
the reductions of Part V were evaluated. Section 21 revises it. Where
the two differ, section 21 governs.*

The inferential value of decomposition is conditional, not universal,
and the condition is a joint function of the biomarker's correlation
structure and the trial design, not of either alone. Three points
summarise the result.

1. **The design contrast matters, not the parametric model.** The
   trial must supply the drug-on/off and belief contrasts that identify
   the components. The lumped one-component analysis already exploits
   the drug-on/off contrast through `Dbc`. Adding parametric component
   structure on top is not free: it spends degrees of freedom and, under
   an uncontaminated biomarker, reduces no bias while losing substantial
   power.

2. **For the biomarker-by-treatment slope, the bias is governed by
   covariance, not magnitude.** A biomarker orthogonal to `PB` and `TV`
   gives an unbiased lumped slope no matter how large `PB` and `TV` are.
   A biomarker coupled to `PB` or `TV` gives a contaminated slope that no
   amount of data fixes, because the failure is identification, not
   precision.

3. **Recovery under contamination is a design property.** When the
   biomarker is coupled to `PB`, the bias is recoverable only under a
   design that varies belief independently of drug exposure. On the
   standard hybrid design no analysis recovers it; on a balanced-placebo
   design the belief-covariate decomposition does.

### The bottom line specific to the hybrid design

For the biomarker-by-treatment interaction estimand on the hybrid
design, component decomposition never wins:

- if the biomarker is orthogonal to `PB` and `TV`, the lumped analysis
  is already unbiased and is the more powerful of the two, so
  decomposition only costs power;
- if the biomarker is coupled to `PB`, the lumped analysis is biased
  but decomposition cannot rescue it on the hybrid design, because the
  identifying contrasts are collinear.

Decomposition earns its place on the hybrid design only for a
different estimand: the **attribution** of response to pharmacology
versus belief versus natural history. For attribution the hybrid
design is the primary identification vehicle, subject to the $\eta$
limit of section 7.3 and the collinearity caveat of section 7.4.

### Was the decomposition in Hendrickson's DGP unnecessary?

A natural reading of the paper-06 result is that, because the lumped
analysis suffices for the biomarker-by-treatment estimand, the
three-component structure Hendrickson built into the simulator was
itself superfluous. That reading conflates two distinct choices. The
decomposition was unnecessary in the *analysis* model; it was
load-bearing in the *data-generating process*. The two claims are not
interchangeable.

The DGP decomposition was necessary, or at least not removable, for two
reasons.

- It defines the ground truth the whole programme tests against. The
  omitted-variable-bias identity of section 11.1, and every simulation
  that confirms it, can be exhibited only if the simulator generates
  `BR`, `PB`, and `TV` as separately controllable channels. A DGP that
  produced a single lumped treatment effect would offer no $c_{bm,PB}$
  knob to turn, no way to construct the contamination, and therefore no
  way to show that the lumped analysis is unbiased under orthogonality
  but biased under coupling. The DGP decomposition is what makes the
  bias question askable.
- It encodes genuinely distinct mechanisms with different time-drivers:
  `BR` by cumulative time on drug with off-drug carryover decay, `PB` by
  time on belief scaled by the phase expectancy $\eta$, and `TV` by time
  since entry. These produce different trajectories across trial phases,
  and the design-contrast structure of section 7.2 (open-label versus
  blinded yielding approximately `0.5 * PB`, blinded on-drug versus
  blinded placebo yielding `BR`) cannot be reproduced from a single
  lumped curve. The components are not cosmetic; they drive the phase
  behaviour the whole identification argument relies on.

What paper 06 shows to be unnecessary is only the symmetric move of
mirroring that structure in the analysis. Fitting the components back
out is unnecessary, and under an orthogonal biomarker counterproductive,
for the biomarker-by-treatment estimand, and unrecoverable on the hybrid
design when the biomarker is coupled. Hendrickson's choice to analyse
with a single lumped `Dbc` was, in hindsight, the appropriate one for
that estimand.

The qualifier is the estimand. The analysis-side lumping is validated as
sufficient only for the interaction-slope target under an orthogonal
biomarker. Had the scientific target been attribution, or a `PB`-coupled
biomarker, the lumped analysis would be insufficient, and the DGP's
decomposition would then require both a matching analysis and a
belief-decoupling design. The accurate summary is therefore that the DGP
decomposition was necessary to pose and answer the question, whereas the
analysis-side decomposition was unnecessary for the specific estimand
Hendrickson actually pursued.


## 13. Practical advice for the hybrid design

The hybrid design is the right default for aggregated N-of-1 trials
whose primary estimand is the pharmacological treatment effect or a
biomarker-by-treatment interaction with a mechanistically
pharmacological biomarker. It is not the right design when the
biomarker may couple to placebo responsiveness and the pharmacological
slope must be recovered. The following guidance operationalises that.

**Use the hybrid design, with the lumped one-component analysis, when:**

- the estimand is the lumped treatment effect or a biomarker-by-drug
  interaction, and
- the candidate biomarker is mechanistically orthogonal to belief and
  to natural history. Pharmacokinetic or receptor-binding biomarkers
  (for example a CYP metabolizer phenotype, or an adrenergic-receptor
  variant affecting prazosin binding) qualify: enzymatic metabolism and
  receptor affinity cannot plausibly couple to belief state or untreated
  disease trajectory. Here the lumped slope is already the target, and a
  parametric component model is counterproductive.

**Use the hybrid design, but add an open-label placebo phase or an
expectancy measure, when:**

- the estimand is attribution (how much of the response is
  pharmacological versus belief versus natural history). The base
  two-phase hybrid identifies only the product `eta * m_PB`; the
  open-label placebo phase or a recorded expectancy separates them. Report
  attribution sensitivity to $\eta$ across its plausible range, and verify
  the off-drug window is long enough, relative to the carryover
  half-life, for `TV` to separate from a decaying `BR`.

**Do not rely on the hybrid design; field a balanced-placebo
(open-hidden) design instead, when:**

- the estimand is the pharmacological biomarker slope and the
  candidate biomarker may correlate with placebo responsiveness or with
  natural-history drift. Psychological or symptom-severity biomarkers
  fall here: baseline PCL-5 severity, for example, couples to both
  expectancy and regression-to-the-mean `TV`, and a COMT genotype is a
  documented predictor of placebo-response magnitude. On the hybrid
  design the contamination is unrecoverable; only the covert-drug and
  open-placebo phases of the balanced-placebo design decouple belief from
  drug exposure enough to recover the unbiased slope.

**Design-stage checklist before committing to the hybrid design.**

- Identify the estimand first: lumped effect, attribution, or
  pharmacological biomarker slope. The estimand, not convenience,
  selects the design.
- Classify the biomarker as orthogonal to or potentially coupled with
  `PB` and `TV`, on mechanistic grounds, before enrolment. If coupling
  cannot be ruled out, the hybrid design is the wrong instrument.
- Verify that every component the estimand requires has at least one
  phase or contrast in which it is uniquely present. A pure open-label
  design cannot separate `BR` from `PB`; a pure blinded design cannot
  reach full-strength `PB`.
- If attribution is in scope, plan the open-label placebo phase or
  expectancy measurement at design time, and size the off-drug window
  against the carryover half-life.
- Do not plan to buy bias protection with the analysis model. On the
  hybrid design the phase-augmented and full-decomposition analyses do
  not reduce biomarker-slope bias; design is the only lever.

# Part IV. The decomposition in the wider literature

Parts I through III have developed the decomposition on its own
terms. We turn now to the question of how it sits in relation to the
clinical trials and N-of-1 literatures, since a reader deciding
whether to adopt it, or which of its components to carry, will want
to know what is conventional and what is not.

We should say at the outset that the three components are not equally
well established. One of them is standard, one is contested, and one
is an innovation of this program. That asymmetry turns out to matter
for the reductions considered in Part V, and it is the reason this
part precedes them.

## 14. Placebo response and the three-arm tradition

The separation of placebo response from natural history has a long
and contentious history in clinical trials methodology, and it is
worth recalling why.

A two-arm trial compares active treatment against placebo. The
difference between the arms estimates the pharmacological effect,
which is what regulators want. But what does the placebo arm itself
measure? It is tempting to read improvement in the placebo arm as the
placebo response, that is, as the benefit conferred by the belief
that one is being treated. This reading is mistaken, or at least
unsupported, because the placebo arm also contains everything that
would have happened anyway. Regression to the mean, spontaneous
remission, seasonal variation, and the ancillary care that trial
participation itself supplies are all present in the placebo arm and
none of them is a placebo effect.

Separating the two requires a third arm in which participants receive
no treatment at all. Hrobjartsson and Gotzsche assembled the trials
that had done this, and their conclusion was deflationary. Across
their reviews, the last covering 202 trials in 60 conditions, they
found little evidence of a placebo effect on objective outcomes, and
what remained on subjective outcomes was modest. Much of what the
field had called placebo response was, on their reading, natural
history and regression to the mean wearing a different name.

We do not need to adjudicate that debate here. The relevant point for
present purposes is narrower. The empirical literature regards the
belief channel and the natural-history channel as separable only
under a design that most trials do not field. A simulation that
instantiates both as distinct latent trajectories is therefore
asserting an identification which is, at best, difficult to obtain in
practice.

This bears directly on the hybrid design of section 7. That design
supplies an open-label and a blinded phase, and the contrast between
them identifies the product of the expectancy weight and the
placebo-belief maximum. It does not supply a no-treatment arm. As
section 7 sets out, the working value of the expectancy weight in the
blinded phase is a modeling assumption rather than an estimate.

## 15. Trend terms in the N-of-1 and crossover literature

The natural-history component stands on quite different ground.

A time trend is a standard element of the crossover and N-of-1
model. Period effects, drift, and serial correlation are routine in
that literature, and the reason is straightforward. A within-person
design measures the same participant repeatedly over weeks or
months, and any process that moves the outcome over that span will
be confounded with the treatment schedule unless it is modeled. Senn
treats period effects as a first-order concern throughout his work on
crossover trials. Schmid and Yang's Bayesian hierarchical formulation
for N-of-1 data carries an individual-level trend parameter alongside
the treatment effect and the autocorrelation.

The natural-history component of this decomposition is the
formalization of that trend. It is more elaborate than a linear drift
term, since it takes a Gompertz form and carries participant-specific
random effects, but it occupies the same position in the model and
answers the same concern.

We note one consequence which is easy to overlook. The analysis model
fitted throughout this program includes a linear term in time. That
term is the analysis-side counterpart of the natural-history
component. There is no corresponding term for placebo-belief, since
the expectancy regressor exists in the software but is disabled in
the simulation drivers. The decomposition is therefore matched to the
analysis model in one component and not in the other.

## 16. What this program adds, and what it assumes

A modeled placebo-belief trajectory, with its own Gompertz parameters
and its own participant-level random effects, is not inherited from
either literature. It is an innovation of the Hendrickson framework
and of this program, and it is worth being explicit about what it
buys and what it costs.

On the plus side, it makes the attribution question askable. A trial
that can separate pharmacological response from expectancy response
can ask whether a candidate biomarker predicts the former or merely
the latter, and that is a question with real regulatory and clinical
consequence. Part III's worked example shows what is lost when the
question cannot be posed.

On the minus side, the identification rests on the expectancy weight,
which the design does not identify. It also asserts a separation
which, as section 14 records, the empirical placebo literature
regards as attainable only under designs this program does not
simulate. The assumption is not unreasonable, but it is an assumption
and it should be labeled as one.

The reader may reasonably ask how much of the framework depends on
it. That is precisely the question Part V takes up, by removing one
component at a time and observing what changes.

---

# Part V. The two-component reductions

*This part is new. It reports and interprets papers 13 and 14, and
where it disagrees with Parts I to III the disagreement is flagged
explicitly rather than smoothed over.*

## 17. The question

The three-component decomposition is expensive. At eight measurement
occasions it makes the correlation matrix `(2 + 3n) = 26` square, with
three within-component blocks and three cross-component blocks, and it
imposes positive-definiteness constraints that bound the
biomarker-moderation parameter.

For a study whose object is the *comparison* between data-generating
architectures, or between analysis specifications, it is fair to ask
whether all three components earn their place. Two reductions to two
components are available:

| paper | components retained | dropped |
|---|---|---|
| 13 | BR + TV | placebo-belief (PB) |
| 14 | BR + PB | natural history (TV) |

Both take the matrix to `(2 + 2n) = 18` square with one cross-component
block. Both reproduce paper 01's architecture comparison: covariance
moderation loses power as carryover grows, mean moderation does not,
and the crossover design is insensitive to carryover under either.
**Neither reduction changes the substantive conclusion**, so the choice
between them is not a choice about results. It is a choice about which
assumption to carry.

## 18. The four axes

The reductions differ on four axes, and they do not all point the same
way.

| axis | TV-drop (paper 14) | PB-drop (paper 13) | favours |
|---|---|---|---|
| **Safety** | unconditional | conditional on `c.bm.pb = 0` | 14 |
| **Design structure** | nothing becomes inert | expectancy weight inert | 14 |
| **Analysis-model match** | `t` term fits nothing; PB unmodelled | DGP matches the fitted model | 13 |
| **Convention** | deletes the conventional trend term | matches N-of-1 practice | 13 |

### 18.1 Safety

The biomarker cannot couple to TV. The implementation provides an
optional biomarker-PB correlation (`c.bm.pb`), used in paper 06's
contaminated-biomarker study, but **there is no biomarker-TV
counterpart anywhere** in the package, the drivers, or any
implementation collection. The biomarker-TV entry of the correlation
matrix is identically zero under every parameter setting the software
can express.

By paper 06's omitted-variable identity, the one-component biomarker
slope is

    beta_bm^BR  +  w_PB * beta_bm^PB  +  w_TV * beta_bm^TV

so a component contributes bias only through the biomarker's
covariance with it. Deleting a block whose coupling is *structurally*
zero cannot change the estimand. Deleting a block whose coupling is
merely *set* to zero changes the estimand as soon as someone sets it
otherwise.

Paper 14's reduction is of the first kind, paper 13's of the second.
A study setting `c.bm.pb` nonzero cannot use paper 13's DGP at all.

**Caveat that matters.** Paper 14's safety rests on a *gap in the
software*, not on a property of natural history. Under the same
identity, a biomarker correlated with natural-history drift would
displace the slope exactly as a PB-correlated biomarker does. Closing
that gap, by adding the missing parameter, would both test the
untested half of paper 06's central claim and tell us how much
paper 14's advantage is really worth.

### 18.2 Design structure

The design's expectancy weight $\eta$, which distinguishes open-label
from blinded occasions, enters the data-generating process at exactly
two points, both inside the PB component: it scales the PB mean
trajectory and the PB standard deviation. It appears nowhere else.

With PB removed, $\eta$ has no effect on the generated data. The
open-label versus blinded contrast becomes invisible to the
simulation, and the three trial designs are distinguished only by
their on-drug and off-drug patterns.

The mirror is not true. Dropping TV renders **nothing** inert: the
calendar-time argument that drives the TV mean also drives the AR(1)
week gaps of every within-component block, so unequal visit spacing
survives, and $\eta$ continues to act through PB.

This asymmetry is one-sided and it favours paper 14.

**Mitigation, from Part II.** Section 7.3 above establishes that the
hybrid design identifies only the product `eta * m_PB`, and that the
working value $\eta \approx 0.5$ is a modelling assumption rather
than an estimate, not falsifiable from outcome data without an
expectancy measurement. If $\eta$ is not identified in the first place,
paper 13's cost in rendering it inert is smaller than it first
appears. The cost is real for any study that *manipulates* expectancy,
as the balanced-placebo designs do, and largely notional for a study
that merely carries a blinded phase.

### 18.3 Analysis-model match

This axis was not visible until the analysis model was read alongside
the DGP, and it reverses the balance of the first two.

The fitted model is assembled as

    modelbase <- "Sx ~ bm + t"          # always
    ... + Dbc + bm:Dbc                  # when drug status varies
    ... + De                            # only when useDE = TRUE

Both papers' drivers set `useDE = FALSE`. So:

- **TV has a fitted analysis-side counterpart**, the `t` term, present
  in every fit.
- **PB has none** in these runs. The expectancy regressor `De` exists
  but is disabled.

The consequence is that the two reductions differ in how well the DGP
matches what the analysis actually estimates:

- **Paper 13 (BR + TV).** The DGP generates TV, and the analysis fits
  a trend term that partially absorbs it. The DGP generates no PB, and
  the analysis models no PB. The two are coherent.
- **Paper 14 (BR + PB).** The DGP generates PB, and the analysis
  models it nowhere, so PB is lumped into the residual. The DGP
  generates no TV, and the analysis fits a `t` term with nothing to
  fit.

Part I warns specifically against the first of those: lumping a
structured participant-level component into noise produces residuals
that are autocorrelated within participants and heteroscedastic across
them, with standard errors wrong in both directions. Paper 14 does
this to PB.

The effect did not show up as a Type I error problem in either paper's
results, which is reassuring but not decisive: both papers' null cells
show the same design-dependent calibration pattern, so any additional
distortion from unmodelled PB is not separable from it in the runs
performed.

### 18.4 Convention

A time trend is a standard element of crossover and N-of-1 models.
Period effects, drift, and serial correlation are routine, and the TV
component is the formalization of that trend. A modelled
placebo-belief trajectory with its own Gompertz parameters is an
innovation of this program rather than an inherited convention;
the wider literature handles placebo response by blinding.

The placebo literature reinforces the point. Hrobjartsson and
Gotzsche's systematic reviews, the last covering 202 trials across 60
conditions, find that much of what is attributed to placebo response
is not separable from regression to the mean and natural-history drift
without a no-treatment arm. Their three-arm methodology exists because
the two channels are hard to distinguish empirically. A simulation
instantiating both as separate latent trajectories asserts an
identification the empirical literature regards as attainable only
under designs most trials do not field.

On this axis paper 13's model is the more familiar and the less
committed, and paper 14's is neither.

## 19. Conflict with section 4

Section 4 contains a passage headed *"Can we assume no TV effect?"* whose
answer is "rarely safely, and never without empirical support." It
lists four conditions under which TV may be set to zero and observes
that "for most chronic trials in psychiatric or pain conditions, at
least one of these conditions fails." It then states that setting TV
to zero and treating it as part of the residual "is worse than
estimating it."

**Paper 14 does what that section warns against, in a psychiatric
indication.** The conflict is real and should not be smoothed over.

Two things reduce but do not eliminate it.

First, Part I argues about an *analysis* that omits TV from data which
contains it. Paper 14 removes TV from the *data-generating process*,
so there is nothing to omit and no residual structure to mis-model.
The specific harm Part I names, autocorrelated and heteroscedastic
residuals, does not arise.

Second, paper 14 is explicitly a methodological device rather than a
clinical model, and its stated scope is the comparison between
architectures. Part I's warning is addressed to analysts of real
trials.

What survives the two qualifications is this: **paper 14's DGP is not
a defensible model of a chronic psychiatric indication**, and should
not be described as one. On this program's own taxonomy PTSD is a
cyclical and triggered condition whose natural-history trajectory is
"structured and non-stationary," with selection at moments of clinical
engagement contributing regression to the mean on top. A simulation
that deletes TV cannot represent any of that.

Paper 14's limitations section states this. Readers should not carry
its power figures into any claim about what a real prazosin trial
would achieve.

## 20. What the reductions do not buy

Both papers expected the reduction to relax the positive-definiteness
constraint on the moderation parameter. It does not, and the reason is
analytic rather than empirical.

The constraint is a Schur complement on the biomarker row. Writing `M`
for the component block and `v` for the unit-scaled coupling pattern,
the exact ceiling is

    c_bm* = ( v' M^-1 v )^(-1/2)

Under compound symmetry `M` is a Kronecker sum and its eigenvalues
collapse to those of two small reductions, giving a smallest
eigenvalue of

    (1 - rho) - (c_1 - c_x)

which contains **neither the number of occasions nor the number of
components**. Deleting a third of the rows and columns does not touch
the binding structural constraint.

Numerically the reduction buys between 0.008 and 0.013 of feasible
$c_{bm}$ range across designs and half-lives. By contrast, halving the
number of measurement occasions moves the ceiling by roughly 0.2.

**The case for either reduction is transparency and cost, not
numerical headroom.** Anyone adopting one for the latter reason has
misdiagnosed the constraint.

## 21. Revised guidance

Superseding the guidance of Parts I and III where they differ.

### Use the full three-component decomposition when

- the estimand is the attribution of response to pharmacology versus
  belief versus natural history;
- a candidate biomarker may correlate with a non-pharmacological
  component;
- the study manipulates or measures expectancy, including any
  balanced-placebo or open-hidden design;
- absolute power is reported as a design input rather than as a
  comparison;
- the simulation is offered as a model of a specific clinical
  indication.

### A two-component reduction is defensible when

- the object is a comparison among architectures, analysis
  specifications, designs, or inference procedures;
- absolute power figures are used only relative to one another;
- the reduction and its costs are stated rather than assumed.

### Choosing between the two

Neither dominates; the axes are split two and two.

- **Prefer dropping TV (paper 14)** for work internal to this program
  where the contaminated-biomarker case is live, or where the blinding
  contrast is used, since it is safe under every parameter setting and
  costs no design structure.
- **Prefer dropping PB (paper 13)** for work intended for the wider
  N-of-1 methodological community, since it presents a familiar model,
  avoids asserting a contested decomposition, and is the better matched
  to the analysis model actually fitted.
- **Prefer neither for clinical modelling** of a chronic psychiatric
  indication.

What should not happen is that either reduction be adopted without
noticing that a choice was made.

### Flattening versus removing

A third option exists and is distinct from both. Setting a component's
Gompertz maximum to zero flattens its mean trajectory while leaving
its variance contribution and its covariance blocks in place, because
the implementation sets component standard deviations independently of
the maximum. For PB, flattening also preserves $\eta$'s effect on the
PB standard deviation, so a flattened model **retains the blinding
contrast that a removed model loses**.

A study wanting the simpler mental picture without the change in
residual variance, or without losing the design structure, should
flatten rather than remove.

## 22. Open problems

Two sets of open items are collected here. The first is inherited
from the hybrid-design analysis of Part III, and the second arises
from the reductions of the present part.

### 22.1 Inherited from the hybrid-design analysis

The Study C recovery result is established at 1,000 replicates over
$c_{bm,PB}$ in {0, 0.1, 0.2, 0.3} and $N$ in {70, 150}; recovery beyond
the positive-definite range ($c_{bm,PB}$ above approximately $0.45$),
participant-level component recovery, and replication under an
Architecture A (mean-moderation) coupling remain open. The Study B
coupling-arm evidence is a 100-replicate pilot; its production run is
pre-registered. The $\eta$ working value of 0.5 in the blinded phase
remains a modelling assumption absent a belief measurement. The
framework's participant-level claims are conditional on falling inside
the recoverable region of section 7.4.

### 22.2 Arising from the reductions


1. **Add `c.bm.tv`.** The absence of a biomarker-TV coupling parameter
   is why paper 06's central claim has been demonstrated for PB and
   only asserted for TV, and it is why paper 14's safety argument
   rests on a software gap. Adding it would close both.
2. **Run a matched three-component grid** under the corrected
   carryover recursion, so that the power difference between two and
   three components can be decomposed into its variance and recursion
   parts. Neither paper can currently do this.
3. **Diagnose the design-dependent Type I error.** Both papers find
   pooled rejection rates running from roughly 0.06 under CO to 0.02
   under OL+BDC against a nominal 0.05, architecture-independent and
   monotone in how much off-drug information the design supplies.
   Neither diagnosed it. Whether the cluster-robust correction that
   removes a milder conservatism in paper 02 also flattens this
   pattern is untested.
4. **Test whether unmodelled PB degrades inference** in paper 14's
   configuration, which Section 18.3 identifies as a theoretical
   concern that the existing runs cannot separate from the Type I
   pattern of item 3.


---

---

# Part VI. Reference

## 23. Common questions

### Q: Could I just use a random intercept and random slope per participant?

A random intercept-and-slope mixed model is

$$
Y_{it} \;=\; \alpha_i + \beta_i t + \varepsilon_{it},
$$

with participant-specific intercept $\alpha_i$ and slope $\beta_i$.
This captures heterogeneity in baseline severity and in linear
time trends, but it cannot decompose the time trend into BR-style
saturating-drug, PB-style expectancy-modulated, and TV-style
natural-history pieces. It also forces $\beta_i$ to be linear in
$t$, which is biologically wrong for chronic-condition responses
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

This is a modelling assumption, not a fact about reality. The
Hendrickson framework correlates the baseline biomarker (e.g., a
blood pressure summary) with the BR factor only, and treats the
biomarker as conditionally independent of PB and TV. The
justification is that the biomarker is intended to be a
*pharmacological* predictor: it measures something about the
participant's likelihood of responding to the drug's mechanism,
not their suggestibility or their natural-history trajectory.

If a candidate biomarker turned out to be correlated with PB,
that would be a substantive finding (the biomarker is partly
psychological in its action), and the model could be re-specified
to allow that correlation. The default specification reflects the
prior that a useful predictive biomarker should be a clean
pharmacological predictor; deviations from that prior are testable.

We must add a qualification which is not a matter of modeling
preference. The implementation supplies an optional biomarker to
placebo-belief correlation, and paper 06 exercises it. It supplies no
biomarker to natural-history counterpart at all. The independence of
biomarker and PB is therefore an assumption that can be relaxed,
whereas the independence of biomarker and TV is at present a property
of the software. Section 18.1 of Part V sets out why the difference
matters.

### Q: What if TV is negative (the patient naturally worsens)?

In the standard parameterisation with positive $m_{TV}$, the
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
still evolving. We note separately that adding occasions tightens the
positive-definiteness bound on the moderation parameter, and section
20 gives the closed form. Identifiability and feasibility therefore
pull in opposite directions as the schedule grows. The design contrasts table earlier in this
document is the formal version of the answer: each component
needs at least one phase in which it is uniquely present (or
uniquely absent in a comparable arm), and the more such phases the
better.

### Q: What if my drug has a fast onset and a long washout?

The Gompertz parameters absorb this asymmetry. A drug with a fast
$BR$ onset has a high $r$; a drug with a long pharmacological
washout has a long off-drug tail captured by the carryover
half-life $t_{1/2}$ used in the analysis-side drug indicator
$D_{it} = \exp(-\lambda \cdot t_{sd})$. The decomposition is
agnostic to the specific pharmacokinetics so long as the trial
includes enough timepoints in the relevant transition windows
(on-to-off and off-to-on) to constrain the curves.

---


## 24. Why this matters for biomarker validation

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
$\hat{\beta}_{bm}$ contains a placebo-prediction component that is
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


## 25. Summary table

| Component | Cause | Timescale | Identified by | Inferential role |
|---|---|---|---|---|
| **BR** | Drug pharmacology | Pharmacokinetic (hours-weeks) | Phase contrasts (on-drug vs. blinded placebo) | Pharmacological efficacy, biomarker-prediction target |
| **PB** | Belief in treatment | Expectancy-modulated (hours-weeks) | Open-label vs. blinded contrasts | Placebo response, accounted for, not the target |
| **TV** | Natural history | Disease-trajectory (weeks-months) | Off-drug timepoints | Confounding control |
| **$\varepsilon$** | Measurement noise | Per-observation | Residuals | Fit diagnostics |

---

### Key takeaway

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

---

## 26. References

- Hendrickson, R. C., et al. (2020). Optimizing aggregated N-of-1
  trial designs for the detection of biomarker-treatment
  interactions. Source of the prazosin-PTSD Gompertz calibration used
  throughout this document.
- Hrobjartsson, A. and Gotzsche, P. C. (2001, 2004, 2010). Placebo
  interventions for all clinical conditions. The 2010 review covers
  202 trials across 60 conditions and is the principal source for
  section 14.
- Kirsch, I. (2008) and Locher, C., et al. (2015). Antidepressant
  trial analyses bearing on the separation of placebo response from
  natural-history drift.
- Senn, S. Cross-over Trials in Clinical Research, and related work
  on period effects and on variance components in N-of-1 designs.
  Principal source for section 15.
- Schmid, C. H. and Yang, J. (2022). Bayesian hierarchical methods
  for aggregated N-of-1 trials, carrying an individual-level trend
  parameter alongside treatment effect and autocorrelation.
- Colloca, L., et al. (2004) and Rohsenow, D. J. (1981).
  Balanced-placebo and open-hidden designs, cited in section 7 as the
  established means of separating drug exposure from belief.
- pmsimstats team. Paper 06, Three-component decomposition of
  treatment response in aggregated N-of-1 trials
  (`analysis/report/06-component-decomposition/`). Sections 4
  through 7 supply the identifiability conditions, the
  omitted-variable-bias identity, and Studies A, B and C.
- pmsimstats team. Paper 13, dropping the placebo-belief component
  (`analysis/report/13-two-component-drop-pb/`), and Paper 14,
  dropping natural history
  (`analysis/report/14-two-component-drop-tv/`). The evidence base
  for Part V.
- pmsimstats team. Component-decomposition pedagogy
  (`docs/24-component-decomposition-pedagogy.md`). The standalone
  introduction from which Parts I through III are drawn.
