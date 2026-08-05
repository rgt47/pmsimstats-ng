---
geometry: margin=0.85in
fontsize: 10pt
---

# White paper: the three carryover analysis specifications, and why these three

*2026-08-05 10:02 PDT*

**pmsimstats team.** A design-justification note for the
carryover-sensitivity manuscript
(`analysis/report/02-carryover-sensitivity/report.Rmd`). It names the
three analysis-model specifications the study compares, states them in
the canonical notation of `analysis/report/NOTATION.md`, and sets out
why these three, rather than some other three, are the right set to
evaluate. It closes with the one specification the current set omits,
and with three gaps in the canonical notation that the comparison
brings to light.

**Intended reader.** A statistically literate physician: someone
comfortable with regression and with the idea of an interaction term,
but who should not be expected to supply the simulation-methods
background. Statistical terms are glossed at first use.

## 1. The clinical problem, and the choice it forces

In an aggregated N-of-1 trial, each patient is crossed repeatedly
between active drug and comparator, and many such single-patient
series are pooled into one mixed model. The question these trials
exist to answer is whether a baseline measurement predicts who
benefits, that is, whether a **predictive biomarker** identifies the
patients in whom the drug works. Statistically this is a
biomarker-by-treatment **interaction**: the treatment effect is
allowed to differ according to the patient's baseline biomarker value.

The complication is that drugs do not stop working the moment they are
stopped. The effect decays over days to weeks, so measurements taken
during a nominally off-drug period still reflect partial drug
exposure. This residual effect is called **carryover**. In the
prazosin and post-traumatic stress disorder setting this program is
calibrated to, the relevant half-life is on the order of half a week
to a week, while measurements are weekly, so the first off-drug
measurement typically still carries roughly half the drug effect.

This forces a decision on the analyst that has no default answer. The
model needs a variable representing drug exposure at each occasion.
Coding that variable as simply on or off is one choice among several,
and the manuscript exists to determine which choice performs best. The
present note explains why the field of candidates was narrowed to
three.

### 1.1 The three trial designs

The comparison is run under three protocols, which differ in how the
on-drug and off-drug measurement occasions are arranged. Because the
whole problem concerns what happens *after* the drug stops, a design
that generates more post-discontinuation observations gives the
analyst more carryover to handle, and the ranking of the three
specifications turns out to depend on which design is used. All three
use eight measurement occasions per patient. Patients are randomized
between two or four **paths**, meaning fixed schedules of on-drug and
off-drug periods.

- **CO (crossover).** The traditional two-period within-patient
  crossover. One path takes drug for the first four visits and stops;
  the other takes comparator first and starts drug at visit five. Half
  the patients therefore never discontinue at all.
- **OL+BDC (open-label plus blinded discontinuation).** Patients are
  titrated openly on active drug, then blindly switched to placebo.
  Both paths discontinue, but late, so most visits are on drug and
  only two or three follow discontinuation.
- **Hybrid.** The Hendrickson design and this program's reference:
  open-label titration, then blinded discontinuation, then a brief
  crossover at the end. Four paths, and the richest in on-off
  transitions, which is why it is used as the reference design
  throughout.

## 2. Notation

Symbols follow `analysis/report/NOTATION.md`, the compendium's single
source of truth for notation. The plain-language column is added here
for readability and is not part of that file.

| Symbol | `NOTATION.md` meaning | In plain terms |
|---|---|---|
| $Y_{it}$ | observed symptom score for patient $i$ at timepoint $t$ (code column `Sx`) | the outcome being measured |
| $B_i$ | pre-treatment biomarker value (code `bm`) | the baseline test result |
| $b_i$ | standardized biomarker, $(B_i - \bar{B})/s_B$ | the same, rescaled so one unit is one standard deviation |
| $D_{it}$ | binary drug state: 1 on drug, 0 off drug | the simple on/off switch |
| $D_{bc,it}$ | continuous exposure-decayed indicator: 1 on drug, $e^{-\lambda t_{sd}}$ off drug (code `Dbc`) | a dimmer switch that fades out after the drug stops |
| $t_{sd}$ | time since discontinuation (code `tsd`) | weeks since the drug was stopped |
| $t_{1/2}$ | carryover half-life (code `carryover_t1half`) | time for the residual effect to halve |
| $\lambda$ | decay rate, $\ln 2 / t_{1/2}$ | how fast the residual effect fades |
| $u_i$ | patient random intercept | each patient's own baseline level |
| $\varepsilon_{it}$ | observation-level residual | unexplained measurement-to-measurement variation |
| $c_{bm}$ | moderation parameter of the covariance-moderation architecture | how strongly the biomarker predicts benefit |
| $\sigma_{BR}$ | standard deviation of the response component | the scale the effect is expressed on |
| $\theta_{\text{true}}$ | calibrated true value, $-c_{bm}\sigma_{BR}$ | the correct answer the simulation is checked against |
| $N$ | total patients, across all randomization paths | trial size |

A hat marks an analyst-side quantity as against a true one, so
$\hat{t}_{1/2}$ is the half-life the analyst assumes while $t_{1/2}$
is the one the data were generated with. Per rule 5 of
`NOTATION.md`, the response components are non-negative *reductions*
in symptom severity, so a beneficial treatment effect lowers $Y_{it}$
and interaction coefficients are negative.

Three symbols this note requires are not in `NOTATION.md`. They are
used here provisionally and proposed formally in Section 9:
$L_{it}$ for the lagged just-off-drug indicator, $\beta_{bm:D}$
disambiguated from $\beta_{bm:D_{bc}}$ according to which exposure
variable carries the interaction, and the main-effect coefficients
$\beta_b$, $\beta_t$, $\beta_D$, $\beta_L$. Note that $\beta_b$, the
biomarker main effect, is a different object from $\beta_{bm}$, which
`NOTATION.md` reserves for the moderation parameter of the
mean-moderation architecture.

Abbreviations used below, per the `NOTATION.md` glossary: **CO**
(crossover), **BDC** (blinded discontinuation), **OL+BDC**
(open-label titration followed by blinded discontinuation), **Hybrid**
(the Hendrickson design: open-label, then blinded discontinuation,
then a brief crossover), **LME** (linear mixed-effects model),
**MCSE** (Monte Carlo standard error), and **ADEMP** (the Morris,
White and Crowther reporting framework for simulation studies).

All three specifications share one model,

$$
Y_{it} = \beta_0 + \beta_b\, b_i + \beta_t\, t
       + \beta_D\, X_{it} + \beta_{bm:X}\, b_i X_{it}
       + u_i + \varepsilon_{it},
$$

fitted as an LME with a patient random intercept and `corCAR1`
residual correlation, which allows measurements taken closer together
in time to be more similar. In words: the outcome depends on the
biomarker, on time, on drug exposure, and on the product of biomarker
and exposure. That last product term, $\beta_{bm:X}$, is the estimand,
meaning the specific quantity the trial is trying to estimate. It says
how much the treatment effect changes per standard deviation of the
biomarker.

The three specifications differ only in how the exposure variable
$X_{it}$ is built, and in whether an extra term is added.

## 3. The three specifications

### 3.1 Unadjusted

$$X_{it} = D_{it}, \qquad D_{it} \in \{0, 1\}$$

In code, `Sx ~ bm + t + Db + bm:Db`. Carryover is not represented at
all. A measurement taken one week after stopping, still carrying about
half the drug effect, is coded exactly like one taken after the drug
has fully cleared. The estimand is $\beta_{bm:D}$, the coefficient on
$b_i D_{it}$.

This is what an analyst gets by default from a standard mixed model
without thinking about carryover, which is precisely why it belongs in
the comparison: it is the benchmark any remedy must beat to justify
its extra complexity.

### 3.2 Lag-adjusted

$$X_{it} = D_{it}, \qquad L_{it} = D_{i,t-1}\,(1 - D_{it})$$

with $L_{it}$ added as an extra term, in code
`Sx ~ bm + t + Db + bm:Db + L`. The indicator $L_{it}$ equals 1 at the
first off-drug measurement following an on-drug one, and 0 everywhere
else. It flags the occasions most likely to be contaminated by
residual drug and lets the model estimate how much those occasions
differ, rather than requiring the analyst to say in advance how fast
the drug washes out. This is the standard remedy in the
crossover-trial literature, following Jones and Kenward.

Two structural features determine how its results must be read. First,
$L_{it}$ is added only as a main effect, with no $b_i L_{it}$ product
term. It can therefore adjust the average outcome at contaminated
occasions, but it cannot adjust how the *biomarker effect* behaves
there. Second, it involves no half-life, so it is completely unaffected
by what the analyst assumes about the pharmacokinetics. The estimand is
$\beta_{bm:D}$, the same quantity Unadjusted targets.

### 3.3 Exposure-weighted

$$X_{it} = D_{bc,it} =
\begin{cases}
1 & \text{if } D_{it} = 1,\\[2pt]
\exp(-\hat{\lambda}\, t_{sd,it}) & \text{if } D_{it} = 0,
\end{cases}
\qquad \hat{\lambda} = \frac{\ln 2}{\hat{t}_{1/2}}
$$

In code, `Sx ~ bm + t + Dbc + bm:Dbc`. Instead of an on/off switch,
exposure is a dimmer that starts at 1 and halves every $\hat{t}_{1/2}$
weeks after the drug stops. A measurement one week after stopping, at
an assumed half-life of one week, enters as 0.5 rather than 0. This is
the current default in the pmsimstats-ng package.

It is the only one of the three requiring the analyst to commit to
both a decay shape and a half-life. Because those assumptions enter
through the values of the variable rather than through the model
formula, it is effectively a different model for every assumed
half-life even though it is written identically each time. Its
estimand is $\beta_{bm:D_{bc}}$, the coefficient on $b_i D_{bc,it}$,
which is *not* the same quantity as $\beta_{bm:D}$.

### 3.4 The three codings on one patient

The differences are easiest to see on a single schedule. The table
below takes one OL+BDC path and one Hybrid path, and shows what each
specification writes into the model at each visit, assuming a
half-life of one week. Verified by direct computation from the design
definitions in `simulation-core.R`.

**OL+BDC, path A.** Drug is stopped after week 18.

| Week | On drug | $t_{sd}$ | Unadjusted $D_{it}$ | Lag flag $L_{it}$ | Exposure-weighted $D_{bc,it}$ |
|---|---|---|---|---|---|
| 4 | yes | 0 | 1 | 0 | 1 |
| 8 | yes | 0 | 1 | 0 | 1 |
| 12 | yes | 0 | 1 | 0 | 1 |
| 16 | yes | 0 | 1 | 0 | 1 |
| 17 | yes | 0 | 1 | 0 | 1 |
| 18 | yes | 0 | 1 | 0 | 1 |
| 19 | no | 1 | 0 | 1 | 0.50 |
| 20 | no | 2 | 0 | 0 | 0.25 |

Read the last two rows. At week 19 the patient still has about half
the drug effect on board. Unadjusted records a 0, exactly as it does
for a drug-free patient. Lag-adjusted records a 0 as well, but flags
the visit so the model can allow its average outcome to differ.
Exposure-weighted records 0.50, which is the value that actually
reflects the residual exposure. At week 20 the difference persists:
Unadjusted and Lag-adjusted both record 0 for a visit still carrying
about a quarter of the effect, and only Exposure-weighted distinguishes
it.

**Hybrid, path A.** Drug is stopped after week 10, restarted at week
16, and stopped again after week 16.

| Week | On drug | $t_{sd}$ | Unadjusted $D_{it}$ | Lag flag $L_{it}$ | Exposure-weighted $D_{bc,it}$ |
|---|---|---|---|---|---|
| 4 | yes | 0 | 1 | 0 | 1 |
| 8 | yes | 0 | 1 | 0 | 1 |
| 9 | yes | 0 | 1 | 0 | 1 |
| 10 | yes | 0 | 1 | 0 | 1 |
| 11 | no | 1 | 0 | 1 | 0.50 |
| 12 | no | 2 | 0 | 0 | 0.25 |
| 16 | yes | 0 | 1 | 0 | 1 |
| 20 | no | 4 | 0 | 1 | 0.0625 |

This path exposes the weakness of the lag flag directly. Weeks 11 and
20 both receive $L_{it} = 1$, because each is the first visit after a
discontinuation, and the model therefore assigns them a single shared
coefficient. But week 11 is one week after stopping and retains about
50% of the drug effect, while week 20 is four weeks after stopping and
retains about 6%, an eightfold difference. Exposure-weighted separates
them, 0.50 against 0.0625, because it is built from elapsed time,
whereas the lag flag is built from position in the visit sequence and
cannot tell the two apart. This is the mechanism behind the Hybrid
caveat in Section 5.2.

## 4. Why these three

Four arguments. The second is the substantive one.

### 4.1 They span the three possible positions on washout

Write $\phi(t_{sd})$ for the fraction of drug effect still present
$t_{sd}$ weeks after stopping. Every analyst takes a position on
$\phi$, whether or not they state one, and only three positions are
qualitatively distinct:

- **Deny it.** Assume $\phi \equiv 0$, so exposure is strictly on or
  off. This is Unadjusted.
- **Estimate it.** Leave the shape of $\phi$ unspecified but give it
  one free parameter, applied at the first off-drug occasion. This is
  Lag-adjusted.
- **Assume it.** Commit to a specific decay curve and rate,
  $\phi(t_{sd}) = \exp(-\hat\lambda t_{sd})$. This is
  Exposure-weighted.

The three are therefore not an arbitrary sample from a large space of
codings. They are one representative of each available stance toward
the unknown washout, ordered by how much the analyst is willing to
assume.

There is also a formal relationship that makes this ladder continuous
rather than merely rhetorical. As the assumed half-life goes to zero,
the dimmer switch collapses to the on/off switch, since
$\hat{t}_{1/2} \to 0$ gives $\hat{\lambda} \to \infty$ and hence
$D_{bc,it} \to D_{it}$. **Unadjusted is therefore the limiting case of
Exposure-weighted at an assumed half-life of zero**, not a separate
method sitting outside it. At the other extreme, as
$\hat{t}_{1/2} \to \infty$, the dimmer never dims,
$D_{bc,it} \to 1$ everywhere, the on-off contrast disappears and the
model can no longer estimate anything. The manuscript's
half-life-mis-specification sweep is consequently a traverse along a
path whose left endpoint is Unadjusted itself, and Unadjusted should
be read as the anchor of that curve.

### 4.2 They separate the two ways carryover does damage

This is the strongest reason for this particular set, because it makes
the comparison a controlled experiment rather than a horse race.

Carryover degrades the interaction test through two distinct
mechanisms:

- **Dilution.** The exposure variable is a mis-measured version of the
  exposure the patient actually received. Contaminated occasions are
  entered as though they were drug-free, so the estimated interaction
  is pulled toward zero. This is the same phenomenon clinicians know
  as regression dilution, where a single noisy blood-pressure reading
  understates the true association with outcome. The signal shrinks.
- **Added noise.** Those same occasions differ from one another in
  ways the model does not describe, since some are one week off drug
  and others four. That unexplained variation is absorbed into the
  residual $\varepsilon_{it}$, which inflates the standard error of
  the interaction estimate. The noise grows.

Power falls on both counts at once, the signal shrinking while the
noise grows. The three specifications occupy three different cells of
the resulting classification:

| Specification | Dilution | Added noise |
|---|---|---|
| Unadjusted | not addressed | not addressed |
| Lag-adjusted | not addressed | addressed at the first off-drug occasion |
| Exposure-weighted | addressed | addressed |

The placements follow from where each term sits. Lag-adjusted adds
$L_{it}$ to the average-outcome part of the model only, so it can
absorb noise from contaminated occasions, but because the interaction
is still carried by the on/off variable $D_{it}$, the dilution is
untouched. Exposure-weighted instead replaces the variable *inside*
the interaction, so it addresses dilution directly and mops up the
associated noise as a by-product.

Two things follow, and together they are the reason this set is well
chosen:

1. **Comparing Lag-adjusted with Unadjusted isolates the noise
   mechanism.** The two estimate the same quantity, are fitted to the
   same simulated patients, and differ by exactly one term. Any
   difference between them is therefore attributable to the noise
   mechanism alone, with nothing else varying. Controlled contrasts
   this clean are rare in a simulation study.
2. **What is left over is interpretable.** Whatever advantage
   Exposure-weighted shows beyond that first comparison is the value
   of fixing dilution. The design therefore reveals not only which
   method wins but which of the two mechanisms actually matters in
   this setting, which is the more portable finding.

A set of three methods that all attacked the same mechanism, or that
differed on several dimensions at once, would support neither
conclusion.

### 4.3 Each is the working default of a literature

Each specification is the standard approach of a distinct research
tradition bearing on this problem, so the comparison speaks to all
three audiences rather than one:

- Unadjusted is ordinary applied practice in aggregated N-of-1
  analysis, and the default of any off-the-shelf mixed model.
- Lag-adjusted is the crossover-trial tradition, where a lagged
  nuisance term is the orthodox remedy for carryover.
- Exposure-weighted is the pharmacokinetically motivated approach,
  used by Hendrickson and implemented as the pmsimstats-ng default.

Omitting any one would expose the study to the objection that it never
tested the method its readers actually use.

### 4.4 They are equally estimable, so no cell is decided by failure

All three can be fitted in all three trial designs (CO, OL+BDC and
Hybrid) at both sample sizes, need no data beyond what the trial
already collects, and add at most one term to the model. Verified by
direct recomputation from the stored cell-level summaries: across the
full 540-cell production grid, none of the three produced a missing
power value or a missing estimate. This matters because a comparison
in which one method quietly fails to converge in some conditions is
not a comparison of methods but of numerical robustness.

## 5. What the set excludes, and the one genuine gap

### 5.1 Candidates set aside

Other specifications were considered and not carried forward:

- **Washout deletion**, discarding off-drug measurements taken within
  some window of stopping. This is exposure weighting with crude 0/1
  weights and an arbitrary cutoff, so it is a coarser version of
  Exposure-weighted, and it discards observations the other methods
  use.
- **Estimating the half-life from the data** rather than fixing
  $\hat{t}_{1/2}$ in advance. This is the principled answer to
  Exposure-weighted's main weakness and is the natural next step, but
  the half-life is poorly determined at $N = 35$ from a handful of
  off-drug measurements per patient, so it introduces a
  convergence-failure mode the other three do not have.
- **A full set of lag indicators**, one per off-drug occasion. The
  off-drug window here is only three to four measurements long, so a
  full set is nearly indistinguishable from the time term $t$ and
  exhausts the available information.

### 5.2 The genuine gap: the $b_i L_{it}$ specification

The classification in Section 4.2 has an empty cell. No specification
addresses dilution *within the lag family*. That specification would
be

$$
Y_{it} = \beta_0 + \beta_b\, b_i + \beta_t\, t
       + \beta_D D_{it} + \beta_L L_{it}
       + \beta_{bm:D}\, b_i D_{it}
       + \beta_{bm:L}\, b_i L_{it} + u_i + \varepsilon_{it},
$$

in code `Sx ~ bm + t + Db + L + bm:Db + bm:L`. The current
Lag-adjusted specification has $\beta_L$ but not $\beta_{bm:L}$, and
that one missing term is the entire difference.

**Why the extra term reduces dilution.** When the interaction is
carried only by $b_i D_{it}$, the estimate $\beta_{bm:D}$ compares
on-drug occasions against *all* off-drug occasions lumped together.
Because that comparison group still contains contaminated occasions,
which display part of the biomarker effect themselves, the contrast is
blurred. Adding $b_i L_{it}$ takes the most contaminated occasions out
of the comparison group and gives them their own coefficient. What is
left behind is a cleaner comparison, so $\beta_{bm:D}$ recovers more
of the true effect.

This makes the specification the assumption-free twin of
Exposure-weighted. Exposure-weighted asserts that the biomarker effect
at an off-drug occasion is the on-drug effect multiplied by
$\exp(-\hat\lambda t_{sd})$, a fixed curve chosen in advance. The
$b_i L_{it}$ specification instead estimates a free two-step pattern,
full effect on drug and $\beta_{bm:L}$ at the first off-drug occasion,
committing to no decay shape at all. Same goal, opposite stance on
whether the washout is assumed or estimated.

**Would it work in these designs?** The table below reconstructs the
exposure pattern of every randomization path, applying the same
lag rule the pipeline uses, and reports how many weeks after stopping
the flagged occasions sit, alongside those left in the comparison
group. Verified by direct computation from the design definitions in
`simulation-core.R`; $t_{sd}$ is in weeks.

| Design | Path | $D_{it}$ pattern | Flagged ($L=1$) | $t_{sd}$ where $L=1$ | $t_{sd}$ in comparison group |
|---|---|---|---|---|---|
| CO | A | 1 1 1 1 0 0 0 0 | 1 | 1 | 2, 3, 4 |
| CO | B | 0 0 0 0 1 1 1 1 | 0 | none | never on drug |
| OL+BDC | A | 1 1 1 1 1 1 0 0 | 1 | 1 | 2 |
| OL+BDC | B | 1 1 1 1 1 0 0 0 | 1 | 1 | 2, 3 |
| Hybrid | A | 1 1 1 1 0 0 1 0 | 2 | 1, 4 | 2 |
| Hybrid | B | 1 1 1 1 0 0 0 1 | 1 | 1 | 2, 6 |
| Hybrid | C | 1 1 1 0 0 0 1 0 | 2 | 1, 4 | 2, 3 |
| Hybrid | D | 1 1 1 0 0 0 0 1 | 1 | 1 | 2, 3, 7 |

Three consequences follow, and each tempers the expected benefit.

- **The improvement would be partial.** The comparison group still
  contains occasions two or more weeks after stopping, which at a
  one-week half-life retain about 25% of the drug effect and less. The
  extra term removes the worst contamination and leaves the rest, so
  it should reduce dilution rather than eliminate it.
- **Half the crossover sample contributes nothing.** Path B of the CO
  design starts off drug, goes on, and never stops, so it contains no
  flagged occasions at all. The new coefficient would be estimated
  from path A alone.
- **The flag is crude in the Hybrid design.** Paths A and C each carry
  two flagged occasions, one at one week after stopping and one at
  four weeks, holding roughly 50% and 6% of the drug effect
  respectively. A single indicator pools them under one coefficient.
  This is exactly the case that favors using elapsed time directly, as
  Exposure-weighted does, because the flag is defined by position in
  the measurement sequence and is blind to the uneven spacing of the
  Hybrid visits.

**Expected effect, marked as inference.** This specification has not
been run. Reasoning from the residual-exposure arithmetic above, it
should move the estimate at the Hybrid reference cell from about
$-2.51$, where Unadjusted and Lag-adjusted both sit, partway toward
the correct value $\theta_{\text{true}} = -3.6$ without reaching it,
and should recover some power. Whether it closes the roughly
ten-percentage-point gap to Exposure-weighted at that cell is
genuinely open, and the third bullet argues it will not close it under
the Hybrid design specifically.

**Cost.** The code change is small: one branch in `fit_spec()` and one
more level in the grid, with no change to the data-generating side. A
four-arm paper must be consistent across every table it prints,
however, so it requires re-running the Tier 1 factorial (roughly 85
minutes on eight cores), the Tier 2 sensitivity blocks (about four
minutes) and the S6 recalibration, then re-summarizing and
regenerating all figures. The `quick-sim` rerun would not pick up the
fourth arm automatically, since it runs through the package rather
than `simulation-core.R`. The larger cost is editorial: the
manuscript's title, abstract and discussion are all organized around a
two-way contrast, and a fourth arm means restructuring that spine a
second time.

**Recommendation: a staged approach rather than a blind fourth arm.**
Run $b_i L_{it}$ first as a narrow diagnostic, at the reference cells
only, which is a few minutes of compute rather than a rebuild of the
study. Two outcomes are possible, and both are useful:

- **It moves the estimate materially toward
  $\theta_{\text{true}} = -3.6$ and recovers power.** Then it is a
  genuine finding, it may weaken the case for Exposure-weighted (since
  comparable power without committing to a half-life is the more
  attractive option for a trialist), and it justifies the full
  four-arm rebuild.
- **It lands on top of Lag-adjusted.** Then one paragraph and a
  footnote citing the diagnostic close the gap completely, at
  negligible cost, and the study is defended against the objection
  without restructuring anything.

Either way the diagnostic is worth running before the decision is
taken. Note that there are now two candidate fourth arms, not one:
$b_i L_{it}$ as above, and the linear $t_{sd}$ term described in
Section 5.3, which the reference implementation already provides and
which existing output suggests does recover power. If only one
diagnostic is run, the linear $t_{sd}$ term is the better first
choice, because it is already implemented, because it is the
alternative the original authors offered, and because the `quick-sim`
outputs give a prior expectation that it works.

What should not happen is publishing the three-arm comparison
with no acknowledgement of the empty cell: as things stand, a finding
that Lag-adjusted fails to recover power is partly a consequence of
how the lag device was implemented rather than of the idea itself, and
a reader from the crossover tradition is entitled to press on it. If
the grid stays at three, the Methods section must state plainly that
the lag family is evaluated only in its main-effect-only form.

### 5.3 A second omission: the original implementation's own alternative

The reference implementation this program derives from, vendored at
`implementations/original/`, supports **three** carryover
representations, not two, and none of them is the binary lag
indicator. Reading `R/lme_analysis.R`, the fixed-effects formula is
assembled at lines 142 to 153 and the exposure variable at lines 108
to 118, giving:

| Setting | Fitted exposure terms | Equivalent to |
|---|---|---|
| `carryover_t1half = 0`, `simplecarryover = FALSE` | $D_{bc} = D_{it}$, so $D_{it} + b_i D_{it}$ | Unadjusted (Section 3.1) |
| `carryover_t1half = 0`, `simplecarryover = TRUE` | $D_{it} + b_i D_{it} + \beta_{sd} t_{sd,it}$ | no counterpart in this study |
| `carryover_t1half > 0` | $D_{bc,it} + b_i D_{bc,it}$, $D_{bc} = (1/2)^{t_{sd}/\hat{t}_{1/2}}$ | Exposure-weighted (Section 3.3) |

The first and third are mutually exclusive with the second by
construction: supplying both a non-zero half-life and
`simplecarryover = TRUE` raises an error (lines 67 to 68).

Two consequences follow. First, the binary lag indicator $L_{it}$
evaluated here as Lag-adjusted is *not* inherited from the reference
implementation. It was introduced by this project in
`simulation-core.R`, following Jones and Kenward. The present study
therefore adopts two of the original three specifications and
substitutes a device of its own for the third, without evaluating the
third itself.

Second, that unevaluated third option, the linear $t_{sd}$ term, is a
serious candidate. It sits between the other two on the assumption
ladder of Section 4.1, conceding that residual exposure declines with
elapsed time, as Exposure-weighted does, while declining to specify a
half-life or a decay shape, as Lag-adjusted does. It is also the only
alternative the original authors themselves thought worth
implementing, and the `quick-sim` outputs discussed in Section 6 give
a prior expectation that it recovers power.

**A related divergence in the analysis model.** Hendrickson's default
`useDE = TRUE` appends the expectancy covariate $De$ whenever
expectancy varies across occasions, which it does in the OL+BDC and
Hybrid designs though not in CO. The specifications compared here omit
it: `simulation-core.R` passes `expectancies` to the design builder,
so expectancy shapes the simulated data, but $De$ never enters any of
the three fitted formulas. This is orthogonal to the carryover
question and applies equally to all three arms, so it should not
disturb their ranking, but it does mean the fitted model is not
Hendrickson's, and the close agreement between this study's 0.860 and
Hendrickson's reported 0.82 at the reference cell has not been checked
against that difference.

## 6. What the choice has already revealed

Recorded because it bears on whether the design succeeded, not as a
results summary. Epistemic status: recomputed directly from the stored
cell-level summaries under `analysis/data/`, not re-simulated.

**How the comparison is scored.** The study follows ADEMP, which
prescribes a fixed set of performance measures. For a clinical reader
the relevant ones are:

- **Power.** How often the trial correctly detects a biomarker effect
  that is genuinely present. Higher is better; 0.80 is the
  conventional target.
- **Type I error.** How often it declares an effect when none exists.
  This should sit at the nominal 0.05, and a method that achieves high
  power by exceeding 0.05 has not earned it.
- **Bias.** The average gap between the estimate and the true value.
- **Coverage.** How often the 95% confidence interval contains the
  true value. Should be 0.95.
- **Mean squared error.** A single number combining bias and spread.
- **MCSE.** The simulation's own margin of error, since each condition
  is run a finite number of times. Differences smaller than the MCSE
  should not be interpreted.

Each condition, or **cell**, is one combination of design, sample size,
biomarker strength, decay shape and half-life, and is run 500 times.

Within the manuscript's 216-cell scope, Lag-adjusted and Unadjusted
are statistically indistinguishable. Mean power differs by 0.003
against an MCSE of about 0.018, only 3 of 216 conditions differ by
more than 0.02, and mean bias (0.373 versus 0.372), mean squared error
(1.885 versus 1.876) and confidence-interval coverage (0.936 versus
0.938) agree to three decimal places. The same equality holds in every
sensitivity block and in the cluster-robust recalibration cells.

Read through Section 4.2 this is informative rather than merely
negative. It says the added-noise mechanism is negligible in this
setting, and therefore that the whole of the Exposure-weighted
advantage comes from fixing dilution. A two-method comparison could
not have reached that conclusion.

One qualification is required, and it turns out to be a naming
collision rather than a contradiction. A separate high-precision rerun
under `analysis/data/quick-sim/` does show its arm labelled
`A3_lagged` ahead of its `A1_binary` arm by 1.3 to 6.5 percentage
points on paired tests. That rerun does not fit the specification
described in Section 3.2. It calls the package's `lme_analysis()` with
`simplecarryover = TRUE`, which, verified by reading
`implementations/original/R/lme_analysis.R:151-153`, appends
$t_{sd}$ itself as a linear main effect rather than the binary
indicator $L_{it}$. Lag-adjusted adds $\beta_L L_{it}$; the rerun's
arm adds $\beta_{sd} t_{sd,it}$. A continuous $t_{sd}$ term can track a graded decline across every
off-drug occasion; a single binary flag marks one occasion and is
constant thereafter. The first carries substantially more information
than the second, which is a sufficient explanation for why one gains
power and the other does not, and it supersedes the reading that the
gap is merely a matter of test calibration. Test size does differ in
the same direction, the rerun's arm carrying the highest
false-positive rate of its three, up to 0.053 against 0.037 where 0.05
is nominal, so both effects are probably present. Disentangling them
would require fitting both models to the same datasets, which has not
been done.

The practical consequence is that the two sources were never measuring
the same thing, and the label `A3_lagged` in the `quick-sim` outputs
should not be read as this manuscript's Lag-adjusted specification.

## 7. Which comparisons this design supports

A choice of methods determines which conclusions the study is entitled
to draw. Four rules follow from Sections 3 and 4, and they should
govern how the manuscript's tables are read.

**Power and Type I error are comparable across all three.** These are
properties of the test, not of the coefficient being tested, so they
can be compared freely even though the three specifications do not all
estimate the same quantity.

**Bias, mean squared error and coverage are comparable only within the
Unadjusted and Lag-adjusted pair.** Those two estimate the coefficient
on $b_i D_{it}$; Exposure-weighted estimates the coefficient on
$b_i D_{bc,it}$. All three are scored in the manuscript against the
single target $\theta_{\text{true}} = -3.6$, which is the correct
target for Exposure-weighted but not for the other two. The
consequence is easy to misread: Exposure-weighted shows near-zero bias
and nominal coverage while the other two show substantial bias and
under-coverage, and that gap is largely a statement about which
quantity each one targets rather than about which estimator is better
behaved. Those columns should not be used to rank the three.

**The Lag-adjusted versus Unadjusted contrast is the controlled one.**
Same estimand, same simulated patients, one term of difference.
Whatever separates them is the added-noise mechanism, cleanly
identified. This is the single most trustworthy comparison in the
study, and it is the reason the third arm was worth adding.

**Any conclusion about the lag device is bounded by its main-effect
form.** Because the specification lacks $b_i L_{it}$, a finding that
Lag-adjusted does not recover power says that a lag *main effect* does
not recover power. It does not establish that the crossover
tradition's approach fails, only that this implementation of it does.
Section 5.2 sets out what would be needed to close that gap.

## 8. Summary

The three specifications are Unadjusted, which codes exposure as a
simple on/off switch; Lag-adjusted, which adds a flag for the first
off-drug measurement; and Exposure-weighted, which replaces the switch
with a dimmer that fades at an assumed half-life. They are the right
three because they occupy the three distinct positions available on
the unknown washout curve, because two of them estimate the same
quantity and so isolate one damage mechanism under otherwise identical
conditions while the third addresses the other, and because each is
the working default of a literature that will read the result. Their
principal weakness is the empty cell at $b_i L_{it}$, which leaves the
crossover tradition's remedy represented only in a form that cannot
address dilution, and which should either be filled or explicitly
disclaimed.

## 9. Proposed additions to NOTATION.md

Three gaps in the canonical notation are exposed by this comparison.

**9.1 Add $L_{it}$.** The lagged just-off-drug indicator has no symbol
in Part 1. Proposed entry for the treatment-exposure table:

| Symbol | Meaning |
|---|---|
| $L_{it}$ | lagged just-off-drug indicator, $D_{i,t-1}(1 - D_{it})$; 1 at the first off-drug occasion following an on-drug one (code column `L`) |

**9.2 Split the estimand symbol.** Part 1 currently carries a single
entry, $\beta_{bm:D}$, glossed as the biomarker-by-treatment
interaction coefficient and annotated `bm:Dbc`. That conflates two
different coefficients. Unadjusted and Lag-adjusted estimate the
coefficient on $b_i D_{it}$; Exposure-weighted estimates the
coefficient on $b_i D_{bc,it}$. The distinction is what licenses the
controlled comparison of Section 4.2, and it is also why the
manuscript's bias, mean-squared-error and coverage columns are
comparable within the first pair but not across to the third.
Proposed replacement:

| Symbol | Meaning |
|---|---|
| $\beta_{bm:D}$ | interaction coefficient on the binary regressor, $b_i D_{it}$ (`bm:Db`); the estimand of Unadjusted and Lag-adjusted |
| $\beta_{bm:D_{bc}}$ | interaction coefficient on the exposure-decayed regressor, $b_i D_{bc,it}$ (`bm:Dbc`); the estimand of Exposure-weighted |

**9.3 Add the main-effect coefficients.** Part 1 has no symbols for
the intercept or for the non-interaction terms of the analysis model,
so any paper writing the model out must invent them. Proposed:
$\beta_0$ intercept, $\beta_b$ biomarker main effect, $\beta_t$ time,
$\beta_D$ drug-exposure main effect, and $\beta_L$ lagged-indicator
main effect. The entry for $\beta_b$ should note explicitly that it is
distinct from $\beta_{bm}$, which the file reserves for the moderation
parameter of the mean-moderation architecture; the two are easily
confused and mean quite different things.

**9.4 A policy tension to resolve before renaming.** Part 2 of
`NOTATION.md` requires that a stored data value doubling as a
reporting label be kept identical to it, and records `spec` as stored
`E1/E2/E3` and reported unmapped. That rule was enforced on
2026-07-31, when twenty-one `.rds` files were migrated from
`A1/A2/A3` precisely to retire a display-time mapping. Renaming the
specifications reintroduces the layer that migration removed.

The same file supplies the governing precedent, however: the
architectures were renamed to mechanism names while their stored
values (`mean_moderation`, `mvn`, `combined`) were left untouched and
the mapping recorded. That is the operation proposed here, for the
same reason. Follow that precedent rather than migrating twice. Keep
`E1/E2/E3` stored permanently, apply display names through one shared
helper, and mark the Part 2 row as mapped. A second migration would
invalidate the archived drafts, the pre-migration backups, and every
summary quoted in papers 01 and 10, all to satisfy a rule the file
already suspends for the architectures.

---
*Rendered on 2026-08-05 at 10:26 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/02-carryover-sensitivity/whitepaper-three-analysis-specifications.md*
