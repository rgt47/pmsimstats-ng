# The question: a predictive biomarker in aggregated N-of-1 trials

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** The estimand and its biological rationale; history of the
program; why carryover matters.  
**Prerequisites:** fm2

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
- Sources: 24 "Big picture"; 40 Sections 1 to 2; 31 Part IV
  (Sections 14 to 16); 02 Sections 1 and 7; 21 "Biological
  Assumptions"; papers 01 and 06 introductions.

- New to write: the program history; a framing that answers the
  co-author's comments on 02 pages 2 and 3.
-->

<!-- TODO: program history -->

<!-- TODO: framing that answers the co-author's comments on
     docs/02 pages 2 and 3 -->

<!-- docs/02 and docs/27 not found in current or legacy docs;
     skipped. -->

### The clinical question

<!-- From docs/24 "The Big Picture: Why Split Responses?" -->

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
   antagonist that reduces noradrenergic tone during REM sleep, and
   its pharmacological action is independent of whether you know you
   are taking it.
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
drop might be eight points of drug, three of expectation, and one of
natural drift; or it might be two of drug, eight of expectation, and
two of drift. The total is the same. The clinical implications are
completely different.

The BR-PB-TV decomposition addresses this directly. It writes the
observed change in outcome as a sum of three biologically
interpretable components plus measurement noise:

$$
Y_{it} \;=\; \mathrm{BL}_i - [\,BR_{it} + PB_{it} + TV_{it}\,]
            + \varepsilon_{it}.
$$

Here $\mathrm{BL}_i$ is participant $i$'s baseline symptom score,
$BR$, $PB$, and $TV$ are the three response components (each a
non-negative reduction in symptoms), and $\varepsilon$ is
observation-level noise.

### The estimand

<!-- From docs/40 Section 1 (Summary) and Section 2 -->

The design answers a question about individuals, not arms. The
estimand is whether baseline standing blood pressure predicts the
size of a participant's response to prazosin. Every participant is
on drug and off drug at some point, so each participant contributes
a within-person contrast, and the interaction becomes a question
about how those contrasts vary with the biomarker (docs/40,
Section 1).

Hendrickson et al. designed the simulation to plan a trial of
prazosin for posttraumatic stress disorder (PTSD), registered as
NCT03539614. The hypothesis was that standing systolic blood
pressure at baseline predicts the reduction in PTSD symptoms that
prazosin produces. That relationship had been seen in a post hoc
analysis of an earlier parallel-group trial. The estimand is
therefore not a treatment effect but a treatment-by-biomarker
interaction: whether participants with higher blood pressure respond
more to the drug (docs/40, Section 2).

A parallel-group trial answers this question poorly, because a
participant on placebo carries no information about their own drug
response. A design in which every participant is observed both on
and off drug supplies a within-person contrast for each participant.
Hendrickson et al. compared four designs of equal length (20 weeks)
and equal number of assessments (8 after baseline): open label
(OL), open label followed by blinded discontinuation (OL+BDC), a
traditional crossover (CO), and the Hybrid design (docs/40,
Section 2).

### Predictive versus prognostic biomarkers

<!-- From docs/40 Section 1, item 1 -->

The estimand $\beta_{bm:D}$ is the biomarker-by-treatment
interaction. A biomarker is *predictive* if it correlates with the
response to the drug specifically, and *prognostic* if it correlates
with the outcome irrespective of treatment. Both may share the
baseline biomarker value, but their inference roles are distinct:
a predictive biomarker guides patient selection for a specific
therapy; a prognostic biomarker guides prognosis regardless of
therapy. Under the sign convention of `NOTATION.md`, the three
response components are reductions in symptom severity, so they
enter $Y$ with a negative sign. Drug effects and interaction
coefficients are therefore negative, and the moderation parameters
$c_{bm}$ and $\beta_{bm}$ are positive.

### The Hendrickson framework

<!-- From analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd
     Background sections -->

In an aggregated N-of-1 trial each participant moves between drug
and placebo, or between drug and no drug, several times, and the
individual series are analyzed together. Because every participant
contributes both an on-drug and an off-drug response, the design
can estimate how the drug effect varies between people, and so test
whether a baseline measurement predicts who benefits. The motivating
example is prazosin for post-traumatic stress disorder (PTSD).
Higher pretreatment standing systolic blood pressure was associated
with greater symptom reduction on prazosin but not on placebo, which
makes blood pressure a candidate predictive biomarker.

Hendrickson et al. (2020) introduced the simulation framework this
series builds on. Its contributions include:

- **A component decomposition of the symptom score.** The score is
  the baseline level minus three responses, each following a
  modified Gompertz curve: a biological drug response (BR), driven
  by time since the current on-drug period began; an expectancy
  response (PB), driven by time in the trial and scaled by the
  participant's expectancy of receiving drug; a time-varying natural
  course (TV).
- **Design machinery.** It covers open-label, crossover, blinded
  discontinuation and a Hybrid design combining all three. It tracks
  time on drug, time since discontinuation and expectancy for each
  visit.
- **Open software** (the `pmsimstats` R package).

Hendrickson et al. also encoded the predictive biomarker in a
distinctive way. They placed the biomarker in the same multivariate
normal distribution as the response components, and gave it a
nonzero correlation $c_{bm}$ with the drug response at on-drug
visits and zero correlation elsewhere. A predictive biomarker is
then one whose association with the outcome differs between on-drug
and off-drug visits.

### Why carryover matters

<!-- From analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd
     "Why a follow-up is needed" -->

Reusing the published process for a series of simulation studies
exposed two problems. Both trace to the implementation, not to the
idea.

1. **The joint distribution often does not exist.** Under the
   compound-symmetry structure used, the correlation matrix is
   positive definite only for $c_{bm}$ below a ceiling, 0.256 in
   the Hybrid design without carryover. The published nonzero effect
   sizes were 0.3 and 0.6. The code repairs any non-positive-
   definite matrix silently.
2. **Power collapses under any carryover.** In the Hybrid design at
   $N = 70$ and $c_{bm} = 0.3$ it fell from 0.74 to 0.12 at a
   half-life of 0.1 week, when only 0.1% of the drug effect remains
   a week after stopping. This is an artifact of the coupling rule,
   not a property of the designs.

### Placebo response and the three-arm tradition

<!-- From docs/31 Section 14 -->

The separation of placebo response from natural history has a long
and contentious history in clinical trials methodology. A two-arm
trial compares active treatment against placebo. The difference
between the arms estimates the pharmacological effect. But what does
the placebo arm itself measure? It is tempting to read improvement
in the placebo arm as the placebo response, that is, as the benefit
conferred by the belief that one is being treated. This reading is
mistaken, or at least unsupported, because the placebo arm also
contains everything that would have happened anyway. Regression to
the mean, spontaneous remission, seasonal variation, and the
ancillary care that trial participation itself supplies are all
present in the placebo arm and none of them is a placebo effect.

Separating the two requires a third arm in which participants
receive no treatment at all. Hróbjartsson and Gøtzsche assembled
the trials that had done this; across their reviews they found
little evidence of a placebo effect on objective outcomes, and
what remained on subjective outcomes was modest.

The relevant point for present purposes is narrower. The empirical
literature regards the belief channel and the natural-history
channel as separable only under a design that most trials do not
field. A simulation that instantiates both as distinct latent
trajectories is therefore asserting an identification which is, at
best, difficult to obtain in practice (docs/31, Section 14).

### Trend terms in the N-of-1 and crossover literature

<!-- From docs/31 Section 15 -->

The natural-history component stands on quite different ground. A
time trend is a standard element of the crossover and N-of-1 model.
Period effects, drift, and serial correlation are routine in that
literature. A within-person design measures the same participant
repeatedly over weeks or months, and any process that moves the
outcome over that span will be confounded with the treatment
schedule unless it is modeled. Senn treats period effects as a
first-order concern throughout his work on crossover trials. The
natural-history component of this decomposition is the
formalization of that trend (docs/31, Section 15).

### What this program adds, and what it assumes

<!-- From docs/31 Section 16 -->

A modeled placebo-belief trajectory, with its own Gompertz
parameters and its own participant-level random effects, is not
inherited from either literature. It is an innovation of the
Hendrickson framework and of this program.

On the plus side, it makes the attribution question askable. A
trial that can separate pharmacological response from expectancy
response can ask whether a candidate biomarker predicts the former
or merely the latter, and that is a question with real regulatory
and clinical consequence.

On the minus side, the identification rests on the expectancy
weight, which the design does not identify. It also asserts a
separation which the empirical placebo literature regards as
attainable only under designs with a no-treatment arm, which this
program does not simulate. The assumption is not unreasonable, but
it is an assumption and it should be labeled as one (docs/31,
Section 16).

### Biological assumptions underlying the architectures

<!-- From docs/legacy/21 "Biological Assumptions" -->

The choice between interaction architectures reflects fundamentally
different assumptions about the mechanism by which the biomarker
moderates treatment response.

**When mean moderation is appropriate.** Mean moderation is
appropriate when the biomarker governs the magnitude of the drug's
biological effect for each individual participant. Under this
mechanism, the drug effect for a participant is proportional to
their biomarker value: higher biomarker predicts larger treatment
response.

**When covariance moderation is appropriate.** Covariance
moderation is appropriate when the biomarker's predictive value
emerges from shared variance with a drug-responsive physiological
component, and this shared variance decays as the drug is
eliminated. Under covariance moderation the biomarker and drug
response are jointly determined by the same underlying pathway
state, and their correlation weakens off drug. The interaction is a
property of the joint response to the drug, not an inherent
property of the individual.

For simulation studies of predictive biomarker-moderated treatment
effects, the choice should be driven by the clinical hypothesis. A
baseline biomarker that predicts intrinsic drug responsiveness
corresponds to mean moderation; one whose predictive value depends
on ongoing drug exposure corresponds to covariance moderation
(docs/legacy/21, Biological Assumptions).

---

## Derivations

<!-- one row per result: statement | CM label | verification
script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the
numbers agree? -->

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
| 2026-10-15 | Setting and argument | Initial content moves from docs/24, docs/40, docs/31, docs/legacy/21, paper 01 | Phase 2 content move |
