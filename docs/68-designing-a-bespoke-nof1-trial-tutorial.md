---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# Designing an N-of-1 Trial from First Principles: Matching the Schedule to the Drug and the Patients, and How the Hendrickson Hybrid Compares with a Crossover {.unlisted .unnumbered}
*2026-10-30 10:23 PDT*

**Author.** pmsimstats team

**Audience.** Undergraduate statistics students who have taken a
course in linear regression and know what variance, correlation,
power and Type I error are. No clinical background is assumed. The
companion tutorial `docs/67` builds the *analysis* of these trials.
This one builds the *design*. The two can be read in either order.

**Purpose.** A trial design is not chosen from a catalog. It is
derived from three sets of facts: what the drug does and how fast, who
the patients are and what they will accept, and how the outcome
behaves statistically over time. This tutorial shows how each fact
turns into a design decision. It then uses the trial that motivated
this project, prazosin for posttraumatic stress disorder (PTSD), to
show how those decisions produce the Hybrid design of Hendrickson et
al. (2020). Section 9 compares the Hybrid design with a standard
one-cycle crossover and states where each is better.

**Evidence.** Every number below comes from an earlier project
document, which is cited beside it. Most come from `docs/39`, which
evaluates designs with an exact (closed-form) power formula and checks
it by simulation, and from paper 01, which simulates the project's
full data-generating process. Section 11 lists what the evidence does
not cover.

```{=latex}
\clearpage
\tableofcontents
\clearpage
```

## Notation

Symbols follow `analysis/report/NOTATION.md`. The last column gives
the section where each symbol first appears.

| Symbol | Meaning | Section |
|------------|----------------------------------------------------------------|---------|
| $N$ | total number of participants | 2.1 |
| $i$, $t$ | participant and visit | 2.1 |
| $w_t$ | calendar week of visit $t$ | 2.1 |
| $Y_{it}$ | symptom score; lower is better | 2.1 |
| $B_i$ | biomarker measured once before treatment | 2.1 |
| $D_{it}$ | drug state: 1 on drug, 0 off drug | 2.1 |
| $\beta_{bm:D}$ | biomarker-by-treatment interaction, the target of the trial | 2.1 |
| $\Delta_i$ | participant's on-drug mean minus off-drug mean | 2.2 |
| $a_t$ | contrast weights: $1/n_{\text{on}}$ on drug, $-1/n_{\text{off}}$ off drug | 6.1 |
| $n_{\text{on}}$, $n_{\text{off}}$ | on-drug and off-drug visits per participant | 6.1 |
| $t_{od}$ | time on drug since the last start, weeks | 4.1 |
| $t_{sd}$ | time since the drug was stopped, weeks | 4.2 |
| $t_{1/2}$ | response half-life: time for the remaining drug effect to halve | 4.2 |
| $\eta_t$ | design expectancy: 1 open label, 0.5 blinded | 5.1 |
| $\rho_{cs}$ | compound-symmetry correlation, the same for every pair of visits | 6.2 |
| $\rho$ | AR(1) correlation per week, $\rho^{\lvert w_t - w_s \rvert}$ | 6.2 |
| $A$ | matrix of AR(1) correlations between visits | 6.2 |
| $a^\top A a$ | time part of the contrast variance under AR(1) | 6.2 |
| $P$ | number of randomization paths | 6.5 |

Table: Notation

## Glossary

Terms are listed alphabetically. The number in brackets is the
section where the term is introduced.

- **Aggregated N-of-1 trial.** Many single-patient on/off trials,
  analyzed together [2.2].
- **AR(1) correlation.** Correlation between two visits that falls
  geometrically with the time between them, $\rho^{\lvert w_t - w_s
  \rvert}$ [6.2].
- **Balance.** Equal numbers of on-drug and off-drug visits per
  participant [6.1].
- **Blinded.** The participant does not know whether they are taking
  the drug or a placebo [5.1].
- **Blinded discontinuation.** A phase in which every participant
  starts on drug and is switched to placebo, without being told when
  [7.2].
- **Carryover.** Drug effect that persists after the drug is stopped
  [4.2].
- **Compound symmetry (CS).** The same correlation between every pair
  of visits, regardless of the time between them [6.2].
- **Contrast.** A weighted sum of a participant's scores whose weights
  sum to zero, such as the on-minus-off difference $\Delta_i$ [2.2].
- **Counterbalancing.** Arranging sequences so that each treatment
  appears equally often in each period [6.5].
- **Crossover (one-cycle, AB/BA).** Two periods, drug and placebo, in
  randomized order [3].
- **Expectancy.** The participant's belief about being on drug. In the
  project's model it scales the placebo response [5.1].
- **Interaction.** Here, the change in drug effect per unit of
  biomarker [2.1].
- **Interleaving.** Placing on-drug and off-drug visits close together
  in time [6.2].
- **Measurement error.** Visit-to-visit noise in a rating that is not
  shared with neighboring visits [6.2].
- **Natural history.** Change in symptoms over time that would occur
  without treatment [5.3].
- **Onset.** How quickly the drug effect builds after starting [4.1].
- **Open label.** The participant knows they are taking the drug [5.1].
- **Path (sequence).** One planned order of drug states through the
  trial. Participants are randomized to paths [6.5].
- **Period effect.** A systematic difference between trial periods
  unrelated to treatment, such as natural improvement [3].
- **Placebo response.** Improvement caused by the belief in treatment
  rather than by the drug [5.1].
- **Power.** The probability that the trial detects an interaction that
  is present [2.1].
- **Rebound.** Worsening beyond baseline after a drug is stopped [4.4].
- **Response half-life.** The time for the remaining drug *effect* to
  fall by half after the drug is stopped. It can differ from the
  drug's half-life in blood [4.2].
- **Switch.** A change from on drug to off drug, or back, between
  consecutive visits [6.3].
- **Titration.** Raising the dose in steps to the target dose [4.1].
- **Washout.** An unmeasured gap after stopping a drug, so that the
  next measured period starts free of carryover [4.2].

## 1. What a design is

A design for this kind of trial is a set of answers to five
questions.

1. **Which drug states, in which order?** The on/off sequence for each
   group of participants.
2. **Who knows?** Which phases are open label and which are blinded.
3. **When are the measurements?** The visit weeks.
4. **What is randomized?** The number of paths and how participants
   are allocated to them.
5. **How long, and how many visits?** Total duration and burden.

A good design answers them so that the trial estimates the quantity of
interest precisely, without bias, with participants who will enroll
and stay. Precision, bias and feasibility often pull in different
directions. A design is a compromise among them, and the compromise
depends on facts about the drug and the patients. Sections 4 to 6
collect those facts. Section 7 turns them into a design.

## 2. The question decides the design family

### 2.1 The estimand

The trial asks whether a biomarker $B_i$, measured once before
treatment, predicts how much a participant benefits from the drug. In
regression terms it asks whether the interaction $\beta_{bm:D}$ in

$$
Y_{it} = \beta_0 + \beta_B B_i + \beta_T w_t + \beta_D D_{it}
       + \beta_{bm:D} B_i D_{it} + \ldots
$$

is different from zero. The design's job is to make that coefficient
estimable with as little variance as possible.

### 2.2 Why each participant must be seen on and off drug

In a parallel-group trial each participant is on drug or on placebo,
never both. The drug effect of one person is never observed. The
interaction must be estimated by comparing the biomarker slope in the
drug arm with the slope in the placebo arm. Two slopes from two
groups of different people carry all the between-person variation in
symptom level.

If instead every participant is observed both ways, each has their own
on-minus-off contrast

$$
\Delta_i = \overline{Y}_{i,\text{on}} - \overline{Y}_{i,\text{off}}.
$$

The person's own baseline severity cancels out of $\Delta_i$. The
interaction is then the slope of $\Delta_i$ on $B_i$. Because the
largest source of noise has been removed, far fewer participants are
needed. This is the rationale for an **aggregated N-of-1 design**. An
interaction test needs about four times the sample of a main-effect
test of the same size (Brookes et al., 2004), so this saving matters
more for interactions than for ordinary treatment effects.

**Design consequence.** Every participant should spend time both on
and off drug. The remaining questions are how much, in what order,
and when to measure.

## 3. The baseline candidate: the one-cycle crossover

The simplest design in which every participant is seen both ways is
the **two-period crossover**, often written AB/BA. Half the
participants receive drug then placebo; the other half receive
placebo then drug. The project's version, taken from Hendrickson et
al. (2020), has two 10-week blocks, eight visits every 2.5 weeks,
blinding throughout and no washout (`docs/series/ch02`).

| Path | Drug state at weeks 2.5, 5, ..., 20 |
|---|---|
| A | 1 1 1 1 0 0 0 0 |
| B | 0 0 0 0 1 1 1 1 |

Table: The one-cycle crossover (CO)

Its strengths are well known (Senn, 2002).

- **Balance.** Four on-drug and four off-drug visits per participant.
- **Constant blinding.** Expectancy is the same at every visit, so the
  on/off contrast is not mixed up with a change in what the participant
  believes.
- **Period effects are separable.** Because the two sequences are
  mirror images, a systematic change between the first and second
  half of the trial can be estimated and removed.
- **A carryover-free path.** In the placebo-first path, every off-drug
  visit comes before any drug exposure.

Its weaknesses are also well known.

- **Carryover in the drug-first path.** The first placebo visits
  follow the drug directly. Classical analyses of the two-period
  crossover struggle to separate carryover from the treatment effect
  (Grizzle, 1965; Freeman, 1989). A washout gap is the usual remedy.
- **Half the participants start on placebo.** Some patients will not
  enroll on those terms (Section 5.2).
- **One switch per participant.** Each participant's on and off visits
  are far apart in time. Section 6.2 shows why that matters when
  correlation fades with time.

The rest of this tutorial asks how each property of the drug and the
population should change this starting point.

## 4. Properties of the drug

### 4.1 Onset: how fast the effect builds

Some drugs act within hours. Others take weeks. Prazosin must be
titrated, and Hendrickson et al. assumed about 2.5 weeks of titration
each time it is started. The project's model of the drug response
follows a slow-start S-shaped (Gompertz) curve in time on drug,
$t_{od}$. At the published parameters the response reaches 3 percent
of its ceiling after one week on drug, 39 percent after four weeks
and 84 percent after eight (`docs/40`, Section 4).

**Design consequences.**

- An on-drug period must be long enough for the effect to appear. A
  one-week on-drug period would measure almost nothing.
- Each restart costs time. In the model, time on drug resets at every
  restart, so a four-week re-exposure recovers only 39 percent of the
  ceiling. (This reset is a modeling choice. `docs/38` argues the
  response may rebuild from its residual level instead.)
- Slow onset therefore favors **few, long on-drug periods**, and it
  argues for putting a long on-drug phase early, where it also serves
  as titration.

### 4.2 Offset: carryover and the response half-life

When the drug stops, its effect decays. Write $t_{1/2}$ for the
**response half-life**: the time for the remaining effect on symptoms
to halve. It is not the same as the drug's half-life in blood. A drug
can leave the blood in hours while its effect on symptoms lasts
longer (`docs/38`). One half-life after stopping, half the effect
remains. Two half-lives after, a quarter.

An off-drug visit taken before the effect has decayed is
contaminated. It is coded "off" but still carries drug effect. This
shrinks the measured on-minus-off contrast and costs power.

**Design consequences.**

- Off-drug visits should fall several half-lives after a
  discontinuation. `docs/39` (Section 5) finds the best spacing
  between visits to be about two to six half-lives, depending on how
  noisy the measurement is.
- Each switch from on to off creates one contaminated window. A long
  half-life therefore favors **fewer switches**. With a two-week
  half-life under compound symmetry, the best one-switch design has
  power 0.731 and the best seven-switch design 0.559 (`docs/39`,
  Section 4).
- The alternatives are a **washout**, an unmeasured gap that lengthens
  the trial, or an analysis that models the decay. The second option
  is the subject of `docs/67`, Section 8.

### 4.3 Reversibility

An N-of-1 design assumes that the condition returns toward its
untreated state when the drug stops. A drug that cures, or that
changes the course of the disease permanently, cannot be studied this
way. Every off-drug period after the first exposure would still show
the benefit. Reversibility is a precondition, not a design choice.

### 4.4 Safety of stopping

Some drugs cause rebound or withdrawal symptoms when stopped abruptly.
A design that stops the drug must then stop it in a controlled way
and monitor closely afterwards. Frequent visits right after a
discontinuation serve both safety and measurement.

## 5. Properties of the population

### 5.1 Expectancy and the placebo response

Patients who know they are on drug improve more than patients who are
unsure. The project models this as a placebo-response component
scaled by the design expectancy $\eta_t$: 1 when the participant knows
they are on drug (open label), 0.5 when blinded (`docs/series/ch02`).

If the open-label visits are all on drug and the blinded visits are a
mixture, part of the on-minus-off contrast is really an
open-label-versus-blinded contrast. The design has mixed the drug
effect with an expectancy effect. In the Hybrid design the
correlation between the drug indicator and the expectancy variable is
0.59 (`docs/40`, Section 1). The analysis must then model expectancy,
and any error in that model adds noise to the contrast.

**Design consequence.** Keep expectancy the same at the on-drug and
off-drug visits that form the contrast. Blind the contrast.

### 5.2 Willingness to enroll

Patients who are acutely ill may refuse a trial in which they could
start on placebo. Hendrickson et al. included an open-label phase so
that such patients could enroll; they reported that recruitment then
outpaced the budget (`docs/40`, Section 2). A trial that cannot recruit
has no power at all, so this consideration can outweigh any
statistical one.

**Design consequence.** For a symptomatic population, start everyone
on drug, open label. Section 5.1 says this has a statistical cost.
`docs/39` (Section 7) measures it. With eight visits in total, each
open-label visit that *replaces* a blinded visit lowers power. Two
open-label visits cost 0.09 to 0.15 of power under compound symmetry
and 0.09 to 0.13 under AR(1) with 20 percent measurement error. The
same document notes a cheaper option that was not evaluated: *add*
the open-label phase before the blinded visits rather than in place of
them, and leave it out of the on/off contrast.

### 5.3 Natural history

Many conditions improve or worsen on their own over months. If
participants are on drug early and off drug late, a natural
improvement over time looks like a drug effect, or hides one. The
analysis adjusts for time, but only as well as its model of time is
correct.

**Design consequence.** Randomize *when* the drug stops, so that at
some visits some participants are on drug and others are off at the
same calendar time. A comparison between paths at the same visit is
protected by randomization, not by a model of time (`docs/40`,
Section 7).

### 5.4 Dropout and burden

Longer trials and more visits lose more participants. If participants
drop out because they feel worse off drug, dropout depends on
treatment, and the missing visits are not missing at random. Paper 09
studies this in detail.

**Design consequences.** Keep the trial no longer than the onset and
offset require. Collect informative on-drug visits early. Do not make
off-drug periods longer than the carryover requires.

## 6. Statistical properties of the outcome

### 6.1 Balance

The variance of $\Delta_i$ contains a term proportional to
$1/n_{\text{on}} + 1/n_{\text{off}}$. With eight visits, that sum is
0.50 for four and four, 0.53 for five and three, and 0.67 for six and
two. The balance term is therefore 7 percent larger at five and three
and 33 percent larger at six and two (`docs/39`, Section 3).

**Design consequence.** Aim for equal numbers of on-drug and off-drug
visits in the contrast.

### 6.2 How correlation behaves over time

The contrast $\Delta_i = \sum_t a_t Y_{it}$, with weights
$a_t = 1/n_{\text{on}}$ on drug and $-1/n_{\text{off}}$ off drug, has
variance $a^\top C a$, where $C$ is the covariance of a participant's
eight scores. What the design can do about this variance depends on
the shape of $C$. Two shapes bracket the possibilities.

**Compound symmetry (CS).** Every pair of visits has the same
correlation $\rho_{cs}$, whether one week or sixteen weeks apart. Then
the order of the on and off visits is irrelevant. Only the balance
and the blinding matter. Every blinded four-and-four design has the
same contrast variance (`docs/39`, Section 3).

**AR(1).** Correlation fades with the time between visits,
$\rho^{\lvert w_t - w_s \rvert}$. A participant's slowly drifting
"current state" is shared by nearby visits and cancels from a
difference of nearby visits. It does not cancel from a difference of
distant visits. Under AR(1) the variance factors into a part from the
components and a time part $a^\top A a$, where
$A_{ts} = \rho^{\lvert w_t - w_s \rvert}$. Only the time part depends
on the arrangement of visits.

The following R code computes the time part for several designs at
$\rho = 0.8$ per week.

```r
time_part <- function(weeks, on, rho = 0.8) {
  a <- ifelse(on == 1, 1 / sum(on), -1 / sum(1 - on))
  A <- rho^abs(outer(weeks, weeks, '-'))
  drop(t(a) %*% A %*% a)
}
co <- seq(2.5, 20, by = 2.5)
hybrid <- c(4, 8, 9, 10, 11, 12, 16, 20)
time_part(co, c(1, 1, 1, 1, 0, 0, 0, 0))      # crossover
time_part(hybrid, c(1, 1, 1, 0, 0, 0, 1, 0))  # Hybrid, path 3
time_part(co, c(1, 0, 1, 0, 1, 0, 1, 0))      # alternating
```

| Design | Switches per path | $a^\top A a$ |
|---|---|---|
| Crossover (CO) | 1 | 0.828 |
| Open label + blinded discontinuation | 1 | 0.649 to 0.664 |
| Hybrid | 2 to 3 | 0.361 to 0.449 |
| Alternating every visit, every 2.5 weeks | 7 | 0.165 |

Table: Time part of the contrast variance under AR(1), $\rho = 0.8$.
The crossover, Hybrid and alternating values were recomputed with the
code above for this tutorial; all agree with `docs/39`, Section 3.

The alternating design has about a fifth of the crossover's time
part. The Hybrid sits between them, because its weekly
blinded-discontinuation visits put off-drug visits right next to
on-drug ones.

**Design consequence.** If correlation fades with time, interleave on
and off visits. If it does not, interleaving buys nothing and only
adds switches, which carryover penalizes (Section 4.2). Which shape
real symptom data follow is an empirical question. Pilot data or
earlier trials should be used to judge it.

### 6.3 Switches against carryover

Sections 4.2 and 6.2 pull in opposite directions. More switches put
on and off visits closer together, which helps under AR(1). Each
switch also creates a contaminated off-drug window, which hurts when
the half-life is long. `docs/39` (Section 4) gives the best power
attainable with a given number of switches. Under AR(1) with 20
percent measurement error and a one-week half-life, it rises from
0.221 with one switch to 0.435 with seven. Under compound symmetry
with a two-week half-life, it falls from 0.731 to 0.559.

### 6.4 Measurement error

Rating scales add noise that is not shared with the next visit. Under
pure AR(1), visits placed very close together would be almost
perfectly correlated, and the contrast noise would vanish. Measurement
error puts a floor under that. The more of the variance is
measurement error, the less close spacing helps, and the further
apart the visits should be (`docs/39`, Section 5).

### 6.5 Randomization and paths

Randomizing participants to two or more paths does two things.

- It allows **between-path comparisons at the same visit**, which
  identify the drug effect without relying on a model of time or
  expectancy (Section 5.3).
- It allows **counterbalancing**, so that period effects can be
  estimated and removed.

For the interaction specifically, the mixture of paths matters less
than one might expect. In the Hybrid design, putting all 70
participants on one path changed interaction power by at most 0.04
(`docs/40`, Section 7). The interaction is estimated from each
participant's own contrast. The randomization protects against a
biomarker that is related to the natural history or to expectancy.

## 7. Worked example: deriving a design for prazosin in PTSD

This section applies Sections 4 to 6 to the trial that motivated the
project (NCT03539614). The question was whether standing systolic
blood pressure at baseline predicts the reduction in PTSD symptoms
that prazosin produces (`docs/40`, Section 2).

### 7.1 The facts

| Fact | Value | Source |
|---|---|---|
| Onset | slow; titration about 2.5 weeks; 84 percent of effect by 8 weeks | 4.1 |
| Offset | plausibly short pharmacologically (0.1 to 0.2 week); longer behavioral persistence possible | 4.2 |
| Reversibility | symptomatic treatment; benefit expected to fade off drug | 4.3 |
| Population | acutely symptomatic; reluctance to accept placebo at entry | 5.2 |
| Natural history | symptoms may improve over months | 5.3 |
| Correlation | unknown; plausibly fades with time | 6.2 |

Table: Inputs for the prazosin design

### 7.2 The decisions

1. **Start everyone on drug, open label, for eight weeks.** This
   answers the recruitment problem (5.2). It doubles as titration
   (4.1). Eight weeks is long enough for the slow-onset response to
   reach most of its ceiling. Visits at weeks 4 and 8 record the
   on-drug response. The cost is expectancy mixing (5.1), which the
   analysis must model.
2. **Then blind, while everyone is still on drug.** From week 8,
   participants are told they may be switched to placebo. Everyone
   remains on drug at week 9, so the drop in expectancy can be seen
   without a change in drug.
3. **Randomize when the drug stops.** Half the participants stop after
   week 9, half after week 10. At week 10 the two groups differ only
   in drug state, at the same expectancy and calendar time. This is a
   randomized comparison (5.3, 6.5). Weekly visits at weeks 9 to 12
   place off-drug visits next to on-drug ones, which helps if
   correlation fades with time (6.2). They also monitor safety after
   stopping (4.4).
4. **End with a short blinded crossover.** Two four-week periods, drug
   and placebo in randomized order, with one visit at the end of each
   (weeks 16 and 20). This re-exposure adds a second randomized on/off
   comparison at the same calendar time. Four weeks also lets a short
   carryover decay before the off-drug visit (4.2).
5. **Two randomizations make four paths.** Two stop times crossed with
   two crossover orders.

The result is the Hybrid design (Table 5, from `docs/40`, Section 2).

| Path | Drug state at weeks 4, 8, 9, 10, 11, 12, 16, 20 | Weeks since stopping |
|---|---|---|
| 1 | 1 1 1 1 0 0 1 0 | 0 0 0 0 1 2 0 4 |
| 2 | 1 1 1 1 0 0 0 1 | 0 0 0 0 1 2 6 0 |
| 3 | 1 1 1 0 0 0 1 0 | 0 0 0 1 2 3 0 4 |
| 4 | 1 1 1 0 0 0 0 1 | 0 0 0 1 2 3 7 0 |

Table: The Hybrid design. Weeks 4 and 8 are open label; weeks 9 to 20
are blinded.

### 7.3 Different facts, different designs

The same reasoning gives other designs for other drugs and
populations. These are illustrations of the principles, not evaluated
recommendations.

- **Fast onset, short carryover, a stable chronic condition, patients
  who accept placebo.** No open label is needed. If correlation fades
  with time, a blinded design that alternates every two to four weeks
  is close to the best design in every condition `docs/39` examined
  (Section 6 there).
- **Long carryover, compound-symmetry-like correlation.** Few, long
  blocks: a blinded crossover with wide spacing, or with a washout
  before the off-drug visits.
- **Slow onset and patients who accept placebo.** A long blinded
  on-drug period first, then a blinded discontinuation, without the
  open-label phase. Under compound symmetry, running the Hybrid
  schedule fully blinded raised its power from 0.663 to 0.739
  (`docs/39`, Section 3, decayed mean moderation, no carryover).
- **Patients who refuse placebo at entry, and a need for maximum
  power.** Add the open-label phase *before* the eight blinded visits
  instead of replacing two of them, and leave it out of the contrast
  (`docs/39`, Section 7, not evaluated there).

## 8. A design procedure

The steps below summarize Sections 4 to 7. They follow the procedure
in `docs/39`, Section 10.

1. **State the estimand.** Here, the biomarker-by-treatment
   interaction.
2. **Check reversibility.** If the effect does not fade off drug, do
   not use an N-of-1 design.
3. **Estimate onset and the response half-life,** from pharmacology,
   earlier trials or pilot data. Give a plausible range, not a point.
4. **Judge the correlation structure and the measurement-error
   share,** from earlier repeated-measures data on the same outcome.
5. **List the population constraints:** acceptance of placebo at
   entry, tolerance of withdrawal, expected dropout, maximum length.
6. **Blind and balance the contrast.** Equal on and off visits at
   equal expectancy. If an open-label phase is required, add it
   before the blinded visits rather than in place of them.
7. **Choose the number of switches and the spacing** against the
   half-life and the correlation structure (Sections 4.2, 6.2, 6.3).
8. **Randomize the timing of at least one switch,** so that some
   visits compare paths at the same calendar time.
9. **Evaluate candidates by simulation or closed form across the
   plausible ranges,** and prefer a design that stays near the best
   in all of them over one that is best in only one.
10. **Fix the analysis together with the design.** The analysis of
    `docs/67` and the design interact: for example, estimating the
    half-life from the data needs more than one switch per path
    (Section 9).

## 9. The Hybrid design against the one-cycle crossover

Both designs have eight visits over 20 weeks. Table 6 compares them
on each criterion of Sections 4 to 6.

| Criterion | Hybrid | Crossover | Better |
|-------------------------|----------------------------------|----------------------------------|-----------|
| Recruitment of symptomatic patients | everyone starts on drug, open label | half start on placebo | Hybrid |
| Titration | built into the open-label phase, before any blinded visit | falls inside a blinded measured period, at week 10 for the placebo-first path | Hybrid |
| Blinding of the contrast | open-label visits are all on drug; expectancy mixed with drug state | constant expectancy | Crossover |
| Balance of on and off visits | 5 and 3 (paths 1, 2) or 4 and 4 (paths 3, 4) | 4 and 4 | Crossover, slightly |
| Contrast noise under compound symmetry | larger | smaller | Crossover |
| Contrast noise under AR(1) | time part 0.36 to 0.45 | time part 0.83 | Hybrid |
| Robustness to carryover with a binary drug indicator | off-drug visits one to three weeks after stopping | placebo-first path free of carryover; long off-drug blocks | Crossover |
| Estimating the half-life from the data | two or three switches per path | one switch per path | Hybrid |
| Dependence on the model of time | drug state mixed with early trial time | period effect separable by design | Crossover |
| Randomized same-visit comparisons | weeks 10, 16, 20 | every visit | Crossover |

Table: Hybrid and one-cycle crossover compared

The evidence behind the rows:

- **Compound symmetry favors the crossover.** Under CS the contrast
  variance is 35.9 for the blinded crossover against 44.1 and 45.0 for
  the Hybrid paths (`docs/39`, Section 3). With decayed mean
  moderation and no carryover the power of the analysis studied there
  was 0.765 for the crossover and 0.663 for the Hybrid. Most of the
  difference comes from the open-label phase: the same Hybrid schedule
  run blinded reached 0.739. The rest comes from the five-and-three
  imbalance of paths 1 and 2, whose blinded contrast variance is 38.3
  against 35.9 for paths 3 and 4 and the crossover. Under CS the arrangement of visits does
  not matter, so the Hybrid's interleaving brings no gain there.
- **AR(1) favors the Hybrid.** Its time part is about half the
  crossover's (Section 6.2). The gain comes from the weekly
  blinded-discontinuation visits and the extra switches, not from the
  open-label phase.
- **Under the project's full simulation the Hybrid has more power.**
  Paper 01 simulated both designs under the extended DGP, which uses
  AR(1) correlation, and analyzed them with G9 (`docs/67`). At
  $N = 70$, a moderate biomarker effect and a 0.2-week half-life, G9
  power was 0.650 (MCSE 0.015) in the Hybrid design and 0.413 (MCSE
  0.016) in the crossover. Paper 01 attributes part of the crossover's
  deficit to the analysis: the crossover has constant expectancy and
  only one switch per path, which gives G9 little information to
  choose a half-life. Paper 01 did not separate that from the
  correlation structure.
- **A simple analysis loses more to carryover in the Hybrid.** At a
  one-week half-life, an analysis with a binary drug indicator lost a
  fifth to two thirds of its power in the Hybrid and OL+BDC designs,
  and a tenth to a sixth in the crossover (paper 02, Introduction).
  The Hybrid's off-drug visits come soon after stopping. The
  crossover's placebo-first path has none of that contamination.
- **The crossover depends less on the model of time.** In the Hybrid
  design, a flexible time adjustment (a spline or visit as a factor)
  raised interaction power but also raised the Type I error to 0.075
  and 0.080 under one-week carryover, against 0.056 for a linear time
  term. In the crossover all three time adjustments held their size
  (`docs/40`, Section 1).

**Verdict.** Neither design dominates. The Hybrid's statistical
advantage is its interleaving. That advantage exists only if
correlation fades with time, which the project's main simulation
assumes but which must be checked for a real outcome. Its open-label
phase is a statistical cost, paid for recruitment and titration. The
crossover is better when correlation is roughly constant over time,
when the analysis is simple, and when separating period effects
cleanly matters most. If patients will accept placebo at entry and
correlation fades with time, a blinded design that alternates more
often than either is likely to beat both (`docs/39`, Section 6). If
patients will not accept placebo at entry, the Hybrid's open-label
start is the reason the trial can run at all.

## 10. Summary of principles

1. The estimand decides the family: an interaction needs each
   participant on and off drug.
2. The drug decides period lengths: slow onset needs long on-drug
   periods; long carryover needs few switches and wide spacing.
3. The population decides what is feasible: an open-label start
   recruits symptomatic patients, at a measured cost in power.
4. The correlation structure decides the arrangement: interleave if
   correlation fades with time; otherwise only balance and blinding
   matter.
5. Randomize the timing of a switch, to protect the drug-effect
   comparison from the model of time.
6. Evaluate designs across the plausible range of the uncertain
   inputs, and choose one that is robust.

## 11. What the evidence does not cover

- **All numbers are from models.** The closed-form results of
  `docs/39` assume a mean-shift form of the interaction, $N = 70$, two
  sequences of equal size and equal visit spacing. Paper 01 uses a
  correlation form of the interaction and AR(1) correlation. Neither
  was fitted to real PTSD symptom data.
- **The correlation structure is unknown.** The ranking of the
  Hybrid and the crossover reverses between compound symmetry and
  AR(1). This is the most important unknown for a real design.
- **The restart rule is an assumption.** The model restarts the drug
  response from zero at every re-exposure. If the response rebuilds
  from its residual level, short re-exposure periods are less costly.
- **Expectancy enters only the noise.** In the models used here the
  biomarker is unrelated to the placebo response. If it were related,
  an open-label phase could bias the interaction, not only add noise.
- **Dropout and recruitment were not modeled jointly with the
  design.** The recruitment benefit of the open-label phase is
  reported by Hendrickson et al.; it is not a simulated quantity.

## 12. Exercises

1. Use `time_part()` from Section 6.2 to compute $a^\top A a$ for all
   four Hybrid paths and for an ABAB design at weeks 2.5 to 20. Repeat
   at $\rho = 0.5$. Does the ranking of the designs change?
2. Show that for a blinded design under compound symmetry the
   contrast variance depends only on $n_{\text{on}}$ and
   $n_{\text{off}}$, not on their order.
3. With a response half-life of one week, what fraction of the drug
   effect remains at each off-drug visit of each Hybrid path? Which
   path's off-drug visits are least contaminated?
4. A drug reaches 90 percent of its effect within three days and has a
   response half-life of one day. Patients accept placebo. Propose an
   eight-visit, eight-week design and justify each choice with a
   section of this tutorial.
5. Redesign the Hybrid trial so that the open-label phase is added
   before eight blinded visits instead of replacing two of them. How
   long is the trial? What is gained, and what does it cost the
   participants?
6. Explain why randomizing the stop time (week 9 or week 10) protects
   the estimate of the drug effect against natural improvement over
   time, and why it does less for the interaction.
7. A colleague claims that the Hybrid design is more powerful than the
   crossover. Write a one-paragraph reply stating under what
   conditions the claim holds and under what conditions it fails.

## References

- Brookes ST, Whitely E, Egger M, Smith GD, Mulheran PA, Peters TJ
  (2004). Subgroup analyses in randomized trials: risks of subgroup-
  specific analyses; power and sample size for the interaction test.
  *Journal of Clinical Epidemiology* 57(3):229-236.
- Freeman PR (1989). The performance of the two-stage analysis of
  two-treatment, two-period crossover trials. *Statistics in Medicine*
  8(12):1421-1432.
- Grizzle JE (1965). The two-period change-over design and its use in
  clinical trials. *Biometrics* 21(2):467-480.
- Hendrickson RC, Thomas RG, Schork NJ, Raskind MA (2020). Optimizing
  aggregated N-of-1 trial designs for predictive biomarker validation:
  statistical methods and theoretical findings. *Frontiers in Digital
  Health* 2:13.
- Senn S (2002). *Cross-over Trials in Clinical Research*, 2nd ed.
  Chichester: Wiley.
- Project documents: `docs/38` (response half-life), `docs/39` (design
  optimization by closed form), `docs/40` (Hybrid design primer),
  `docs/series/ch02` (design tables), `docs/67` (companion analysis
  tutorial), paper 01 (`analysis/report/01-dgp-mean-moderation-vs-mvn/`,
  rendered 2026-10-29), paper 02
  (`analysis/report/02-carryover-sensitivity/`, rendered 2026-10-29),
  paper 09 (`analysis/report/09-informative-dropout-by-design/`).
