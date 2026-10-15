# Trial designs, randomization paths and exposure clocks

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** OL, CO, OL+BDC, Hybrid; paths P; clocks t_od, t_sd,
t_pb; matched-N correction.  
**Prerequisites:** ch01

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
- Sources: 13 Section 2 (after a check against `58b32a9`); 33
  Section 2.5 (design table); 40 Sections 2 to 3; 24 Component 2;
  37 Section 2; 27 Sections 2a and 2c (cycle designs, Hybrid round
  trip corrected); NOTATION Part 2; CM "Trial designs and treatment
  exposure", "Time arguments of the components".

- New to write: one table per design; a Hybrid schedule figure.
-->

### Design summary tables

The four designs share 8 post-baseline visits over 20 weeks but
differ in blinding, number of paths, and on-drug fractions. The
tables below give the structural parameters for each design; visit
weeks and expectancy weights are from the `buildtrialdesign()` calls
in the published code (inspected, `implementations/original/`).

**Open Label (OL)**

| Parameter | Value |
|---|---|
| Periods | 1 (drug only) |
| Paths | 1 |
| Sequences | 1 |
| Visit weeks | 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20 |
| Expectancy $\eta$ | 1.0 at all visits |
| On-drug visits | 8 of 8 |
| Off-drug visits | 0 of 8 |
| Balance | unbalanced (no within-subject contrast) |

Table: OL design parameters

**Traditional Crossover (CO)**

| Parameter | Value |
|---|---|
| Periods | 2 (one on drug, one on placebo, order randomized) |
| Paths | 2 |
| Sequences | drug-then-placebo (path A), placebo-then-drug (path B) |
| Visit weeks | 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20 (uniform 2.5-wk intervals) |
| Expectancy $\eta$ | 0.5 at all visits |
| On-drug visits | 4 of 8 per path |
| Off-drug visits | 4 of 8 per path |
| Balance | balanced (equal on/off within path) |

Table: CO design parameters

**Open Label + Blinded Discontinuation (OL+BDC)**

| Parameter | Value |
|---|---|
| Periods | 2 (16 wk open label, 4 wk blinded discontinuation) |
| Paths | 2 |
| Sequences | discontinuation at week 17 (path A), at week 18 (path B) |
| Visit weeks | 4, 8, 12, 16, 17, 18, 19, 20 |
| Expectancy $\eta$ | 1.0 at visits 1--4; 0.5 at visits 5--8 |
| On-drug visits | 6 of 8 (path A), 5 of 8 (path B) |
| Off-drug visits | 2 of 8 (path A), 3 of 8 (path B) |
| Balance | unbalanced (majority on drug; short off-drug window) |

Table: OL+BDC design parameters

**Hybrid (N-of-1)**

| Parameter | Value |
|---|---|
| Periods | 3 (open label, blinded discontinuation, crossover) |
| Paths | 4 |
| Sequences | 2 discontinuation times $\times$ 2 crossover orders |
| Visit weeks | 4, 8, 9, 10, 11, 12, 16, 20 (clustered around wk 9--12) |
| Expectancy $\eta$ | 1.0 at visits 1--2; 0.5 at visits 3--8 |
| On-drug visits | 5 of 8 (paths A, C) or variable (paths B, D) |
| Off-drug visits | 3 of 8 (paths A, C) or variable (paths B, D) |
| Balance | moderately balanced; all four paths identical at weeks 4, 8 and 9 |

Table: Hybrid design parameters

### Hybrid design schedule

The table below gives the complete visit timeline for the Hybrid
design. $D_{bc}$ values at $t_{1/2}^a = 0.2$ wk are computed as
$D_{bc} = \exp(-\ln 2 \cdot t_{sd}/0.2)$ off drug after first
exposure, 0 before first exposure, and 1 on drug (derived;
docs/legacy/13, Section 2).

| Visit | Week | Phase | Path A drug | Path B drug | Path C drug | Path D drug | $D_{bc}$, $t_{1/2}^a = 0.2$ wk |
|---|---|---|---|---|---|---|---|
| 1 | 4 | Open label | on | on | on | on | 1.000 |
| 2 | 8 | Open label | on | on | on | on | 1.000 |
| 3 | 9 | Blinded discontinuation | on | on | on | on | 1.000 |
| 4 | 10 | Blinded discontinuation | on | on | off | off | 1.000 (A/B); 0.031 (C/D) |
| 5 | 11 | Blinded discontinuation | off | off | off | off | 0.031 (A/B); 0.001 (C/D) |
| 6 | 12 | Blinded discontinuation | off | off | off | off | 0.001 (A/B); $<$0.001 (C/D) |
| 7 | 16 | Crossover | on | off | on | off | 1.000 (A/C); $<$0.001 (B/D) |
| 8 | 20 | Crossover | off | on | off | on | $<$0.001 (A/C); 1.000 (B/D) |

Table: Hybrid design visit schedule with drug states and $D_{bc}$
at analysis half-life $t_{1/2}^a = 0.2$ wk (inspected,
`buildtrialdesign()` source and docs/40, Section 2)

Notes: (1) Path labels A--D correspond to `pathA`--`pathD` in the
published `buildtrialdesign()` call (ondrug vectors from ch02,
"The four trial designs"). (2) $D_{bc}$ during the crossover
off-drug visits (visits 7--8 for paths B and D, paths A and C
respectively) is essentially zero because the four-week off-drug
gap at $t_{1/2}^a = 0.2$ wk implies $\exp(-\ln 2 \cdot 4/0.2)
\approx 1.2 \times 10^{-6}$. (3) At re-exposure in the crossover
the time-on-drug clock resets to zero; $D_{bc} = 1$ at that visit.

<!-- docs/33 not found in current or legacy docs; skipped. -->
<!-- docs/27 not found in current or legacy docs; skipped. -->

### The four trial designs

<!-- From docs/legacy/13 Section 2 "The Four Trial Designs" -->

Each design spans 20 weeks with 8 assessment timepoints (excluding
baseline). They differ in blinding, randomization structure, and
the number of distinct paths a participant can follow through the
trial.

**Open-Label (OL) -- 1 path.**
All participants receive active drug for 20 weeks, unblinded.
Expectancy is 1.0 at all timepoints. There is no within-subject
drug variation, so the analysis tests whether the biomarker predicts
the rate of improvement over time (`bm:t` interaction).

```r
tdOL <- buildtrialdesign(
  name_longform = "open label",
  name_shortform = "OL",
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = paste("OL", 1:8, sep = ""),
  expectancies = rep(1, 8),
  ondrug = list(pathA = rep(1, 8))
)
```

**Open-Label + Blinded Discontinuation (OL+BDC) -- 2 paths.**
16 weeks open-label on drug, then a 4-week blinded discontinuation
block. In the blinded block, participants are randomized to receive
active drug for the first 1 or 2 weeks before switching to placebo.
Expectancy drops to 0.5 during the blinded phase.

```r
tdOLBDC <- buildtrialdesign(
  name_longform = "open label+blinded discontinuation",
  name_shortform = "OL+BDC",
  timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptname = c("OL1","OL2","OL3","OL4",
                 "BD1","BD2","BD3","BD4"),
  expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)
```

**Traditional Crossover (CO) -- 2 paths.**
Two 10-week blocks: one on drug, one on placebo, order randomized.
All blinded (expectancy = 0.5). No washout period.

```r
tdCO <- buildtrialdesign(
  name_longform = "traditional crossover",
  name_shortform = "CO",
  timepoints = cumulative(rep(2.5, 8)),
  timeptname = c(paste("COa", 1:4, sep = ""),
                 paste("COb", 1:4, sep = "")),
  expectancies = rep(.5, 8),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
)
```

**N-of-1 (Hybrid) -- 4 paths.**
8-week open-label phase followed by a 4-week blinded
discontinuation and two 4-week crossover blocks. Two independent
randomization points create 4 possible paths. The timepoints are
non-uniform: c(4, 8, 9, 10, 11, 12, 16, 20), with weekly
assessments clustered around the blinded-discontinuation
transition.

```r
tdNof1 <- buildtrialdesign(
  name_longform = "primary N-of-1 design",
  name_shortform = "N-of-1",
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptname = c("OL1","OL2","BD1","BD2",
                 "BD3","BD4","COd","COp"),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
  )
)
```

The `buildtrialdesign()` function computes three derived timing
variables from the `ondrug` vector:

- **tod** (time on drug): cumulative weeks on active treatment,
  resetting to 0 when off drug
- **tsd** (time since discontinuation): cumulative weeks off drug
  since last on-drug period, zeroed before the participant is ever
  on drug
- **tpb** (time in positive belief): cumulative weeks where
  expectancy > 0

These drive the Gompertz trajectory calculations in the DGP
(docs/legacy/13, Section 2).

### The Hybrid design in detail

<!-- From docs/40 Section 2 "The trial and its question" -->

**Why each phase is there.**

- **Open label, weeks 0 to 8.** Every participant starts on active
  drug, knowing it. Hendrickson et al. included this phase so that
  acutely symptomatic patients, who would decline possible placebo
  assignment, could enroll. Statistically, the open-label phase
  contributes on-drug visits at full expectancy ($\eta_t = 1$).
  Drug and expectancy response are both present at full strength
  there and cannot be separated within the phase.
- **Blinded discontinuation, weeks 8 to 12.** Participants are
  told they may be switched to placebo during the block. Every
  path is on drug at week 9. Randomization 1 then decides whether
  the drug stops after week 9 or after week 10. Expectancy is 0.5
  throughout the block.
- **Crossover, weeks 12 to 20.** Two 4-week blinded periods, one
  on drug and one on placebo, in an order set by randomization 2.
  Each period has one assessment at its end (weeks 16 and 20).
  Expectancy is 0.5.

**The four paths** are the two discontinuation times crossed with
the two crossover orders. The published code allocates patients in
fixed numbers (18, 18, 17, 17 at $N = 70$) and generates each
path's data separately. Visit by visit:

| Path | Drug state | $t_{od}$ (weeks) | $t_{sd}$ (weeks) |
|------|------------|-----------------|-----------------|
| 1 | 1 1 1 1 0 0 1 0 | 4 8 9 10 0 0 4 0 | 0 0 0 0 1 2 0 4 |
| 2 | 1 1 1 1 0 0 0 1 | 4 8 9 10 0 0 0 4 | 0 0 0 0 1 2 6 0 |
| 3 | 1 1 1 0 0 0 1 0 | 4 8 9 0 0 0 4 0 | 0 0 0 1 2 3 0 4 |
| 4 | 1 1 1 0 0 0 0 1 | 4 8 9 0 0 0 0 4 | 0 0 0 1 2 3 7 0 |

Table: Drug state, $t_{od}$ and $t_{sd}$ by visit for the four
Hybrid paths

Two features of this table matter. First, the paths differ in drug
state only at weeks 10, 16 and 20. At weeks 4, 8, 9, 11 and 12
every path is in the same state. Second, on re-exposure in the
crossover, time on drug restarts from zero rather than continuing
from the earlier run, so it is 4 weeks at the re-exposure visit
(week 16 or 20) (docs/40, Section 2).

### The three response components and their clocks

<!-- From docs/40 Section 3 "The three response components" -->

Hendrickson et al. model each participant's symptom score as
baseline minus the sum of three responses:

$$
Y_{it} = \mathrm{BL}_i - \bigl[\,TV_{it} + PB_{it} + BR_{it}\,\bigr].
$$

Each component has its own mean trajectory and its own clock.

| Component | Clock | Mean at visit $t$ | SD |
|-----------|-------|-------------------|-----|
| $BR$ | time on drug $t_{od}$ | $G(t_{od}; 10.99, 5, 0.42)$, plus carryover | 8 |
| $PB$ | time under expectancy $t_{pb}$ | $\eta_t \, G(t_{pb}; 6.51, 5, 0.35)$ | $10\,\eta_t$ |
| $TV$ | calendar time $w_t$ | $G(w_t; 6.51, 5, 0.35)$ | 10 |

Table: Clock, mean trajectory and SD of each response component

**Expectancy acts twice.** The design expectancy $\eta_t$ scales
both the mean and the standard deviation of $PB$. Blinding at
week 9 halves the expectancy response, and it also halves that
response's variability.

**What differs between paths.** In the Hybrid design every path
has the same visit weeks, and expectancy stays positive
throughout, so $t_{pb} = w_t$ in every path. Hence $TV$ and $PB$
have identical means in all four paths. Only $BR$ differs, through
the drug-state pattern (docs/40, Section 3).

### The placebo-belief response and expectancy

<!-- From docs/24 Component 2 "Placebo-Belief Response (PB)" -->

PB is the change in symptom score that arises from the patient's
belief that they are receiving treatment, independent of what they
are actually receiving. It is the placebo response in the strict
sense: the part of the improvement that would persist if the active
drug were silently replaced with an inert sugar pill, provided the
patient continued to believe they were on the active drug.

PB also follows a modified Gompertz curve, with its own parameters,
in the expectancy-weighted time $t_{pb}$ (time accumulated over
consecutive visits with positive expectancy; it resets to zero at
any visit with $\eta_t = 0$), and its amplitude is scaled by the
expectancy weight $\eta_t$ of the current visit:

$$
PB_t \;=\; \eta_t \, m_{PB}\,
\frac{\exp\bigl(-d_{PB}\,e^{-r_{PB} t_{pb}}\bigr) - e^{-d_{PB}}}
     {1 - e^{-d_{PB}}}.
$$

Values of the expectancy weight:

| Phase | $\eta$ | Reasoning |
|-------|--------|-----------|
| Open-label on drug | 1.0 | Patient knows they are on the active drug. |
| Blinded discontinuation and blinded crossover | 0.5 | Patient believes they may or may not be on drug; expectation is hedged. |
| Open-label placebo | 0.0--0.5 | Patient told they are on placebo retains some residual expectation. Not used by any design in this compendium. |

Table: Expectancy weight $\eta$ by trial phase with reasoning
(docs/24, Component 2)

### The setting of docs/37

<!-- From docs/37 Section 2 "Setting" -->

The Hybrid design has eight post-baseline visits at weeks 4, 8, 9,
10, 11, 12, 16 and 20, with expectancy 1 at the two open-label
visits and 0.5 afterwards, and four paths. $N = 70$ patients are
allocated 18, 18, 17, 17 across the paths, as in the published
code. Every off-drug visit follows an exposure, so every off-drug
visit has $t_{sd} \ge 1$ week.

The biomarker-DR correlation at visit $t$ is:

$$
\mathrm{Cor}(B, BR_t) = c_{bm}\, g_t,
$$

where $g_t$ is the *coupling* at visit $t$: the fraction of the
biomarker effect $c_{bm}$ that is present in the biomarker's
correlation with the drug response at that visit. The interaction
the analysis is meant to detect is the difference between on-drug
and off-drug coupling: if the biomarker predicts the response as
strongly off drug as on drug, there is no biomarker-by-treatment
interaction to find (docs/37, Section 2).

### Sample size and randomization paths

<!-- From analysis/report/NOTATION.md Part 2 "Sample size and
     randomization paths" -->

Rule 4 fixes $N$ as the total. The complication is that
`modelparam$N` does not mean the same thing at every level of the
call stack, and two families of driver read it differently.

`generateData` generates one path: it takes a single `trialdesign`
and treats `N` as the count for that path alone.
`generateSimulatedResults` generates a whole design: it takes `N`
as the total and allocates it, `Ns <- N %/% nP` with the remainder
spread one per path.

Path counts by design:

| Design | $P$ | Paths |
|--------|-----|-------|
| OL | 1 | on drug throughout |
| CO | 2 | drug-then-placebo, placebo-then-drug |
| OL+BDC | 2 | discontinuation at either of two visits |
| Hybrid | 4 | two discontinuation points crossed with two crossover orders |

The **matched-N correction** (applied 2026-08-08): fixing the
per-path count across designs with different $P$ leaves a
comparison that is not matched on $N$. Papers 01 and 09 were
re-run with every design at a common total of $N = 70$. In paper 01
the Hybrid column moved substantially: Hybrid power fell from about
0.99 to 0.84, which roughly tripled its measured carryover loss
(from 1.4 to 6.9 percent under mean moderation; from 5.1 to 15.4
percent under covariance moderation). The paper's headline
architecture contrast survives, but the claim that mean-moderation
losses are uniformly 1 to 3 percent did not.

The lesson: an unmatched design comparison does not merely inflate
one design, it can manufacture a robustness finding out of a
ceiling effect (NOTATION.md, Part 2).

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
| 2026-10-15 | Setting and argument | Initial content moves from docs/legacy/13, docs/40, docs/24, docs/37, NOTATION.md | Phase 2 content move |
