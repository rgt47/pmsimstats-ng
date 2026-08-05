---
geometry: margin=0.8in
fontsize: 10pt
---

# White paper: the three carryover analysis specifications, and why these three

*2026-08-05 10:02 PDT; restructured 2026-08-05 10:46 PDT*

**pmsimstats team.** A design-justification note for the
carryover-sensitivity manuscript
(`analysis/report/02-carryover-sensitivity/report.Rmd`). It names the
three analysis-model specifications the study compares, states them in
the canonical notation of `analysis/report/NOTATION.md`, and sets out
why these three are the right set to evaluate. It also records why a
fourth candidate is retained as a reported side-check rather than an
arm, identifies the one cell the set still leaves empty, and lists
four additions the comparison requires of the canonical notation.

**Revision note.** This document originally presented the middle arm
as the Jones and Kenward lagged-treatment indicator. That arm has been
replaced by the linear washout term supplied by the reference
implementation, for the reasons in Section 5. The lagged indicator is
retained as a side-check (Section 3.4).

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
model needs a variable representing drug exposure at each occasion,
and coding it as simply on or off is one choice among several. The
manuscript exists to determine which choice performs best; this note
explains why the field was narrowed to three.

### 1.1 The three trial designs

The comparison runs under three protocols, which differ in how the
on-drug and off-drug occasions are arranged. Because the whole problem
concerns what happens *after* the drug stops, a design generating more
post-discontinuation observations gives the analyst more carryover to
handle, and the ranking of the specifications depends on which design
is used. All three use eight measurement occasions per patient.
Patients are randomized between two or four **paths**, meaning fixed
schedules of on-drug and off-drug periods.

- **CO (crossover).** The traditional two-period within-patient
  crossover. One path takes drug for the first four visits and stops;
  the other takes comparator first and starts drug at visit five. Half
  the patients therefore never discontinue at all.
- **OL+BDC (open-label plus blinded discontinuation).** Open titration
  on active drug, then a blinded switch to placebo. Both paths
  discontinue, but late, so only two or three visits follow.
- **Hybrid.** The Hendrickson design and this program's reference:
  open-label titration, blinded discontinuation, then a brief
  crossover. Four paths, and the richest in on-off transitions.

## 2. Notation

Symbols follow `analysis/report/NOTATION.md`, the compendium's single
source of truth. The plain-language column is added here and is not
part of that file.

| Symbol | `NOTATION.md` meaning | In plain terms |
|---|---|---|
| $Y_{it}$ | symptom score for patient $i$ at timepoint $t$ (code `Sx`) | the outcome measured |
| $B_i$ | pre-treatment biomarker (code `bm`) | the baseline test result |
| $b_i$ | standardized biomarker, $(B_i - \bar{B})/s_B$ | the same, in standard-deviation units |
| $D_{it}$ | binary drug state: 1 on drug, 0 off | the simple on/off switch |
| $D_{bc,it}$ | exposure-decayed indicator: 1 on drug, $e^{-\lambda t_{sd}}$ off (code `Dbc`) | a dimmer fading after the drug stops |
| $t_{sd}$ | time since discontinuation (code `tsd`) | weeks since the drug was stopped |
| $t_{1/2}$ | carryover half-life (code `carryover_t1half`) | time for residual effect to halve |
| $\lambda$ | decay rate, $\ln 2 / t_{1/2}$ | how fast the residual effect fades |
| $u_i$ | patient random intercept | each patient's own baseline level |
| $\varepsilon_{it}$ | observation-level residual | unexplained variation |
| $c_{bm}$ | moderation parameter, covariance-moderation architecture | how strongly the biomarker predicts benefit |
| $\sigma_{BR}$ | standard deviation of the response component | the scale the effect is expressed on |
| $\theta_{\text{true}}$ | calibrated true value, $-c_{bm}\sigma_{BR}$ | the correct answer to check against |
| $N$ | total patients, across all paths | trial size |

A hat marks an analyst-side quantity, so $\hat{t}_{1/2}$ is the
half-life the analyst assumes while $t_{1/2}$ is the one the data were
generated with. Per rule 5 of `NOTATION.md`, response components are
non-negative *reductions* in severity, so a beneficial treatment
effect lowers $Y_{it}$ and interaction coefficients are negative.

Four symbols used here are not in `NOTATION.md`; Section 9 proposes
them formally. They are $L_{it}$, the lagged just-off-drug indicator;
$\beta_{bm:D}$ disambiguated from $\beta_{bm:D_{bc}}$ according to
which exposure variable carries the interaction; the main-effect
coefficients $\beta_b$, $\beta_t$, $\beta_D$, $\beta_L$,
$\beta_{sd}$; and the stored specification code `E4`. Note that
$\beta_b$, the biomarker main effect, is a different object from
$\beta_{bm}$, which `NOTATION.md` reserves for the mean-moderation
architecture's moderation parameter.

Abbreviations, per the `NOTATION.md` glossary: **CO** (crossover),
**BDC** (blinded discontinuation), **OL+BDC**, **Hybrid**, **LME**
(linear mixed-effects model), **MCSE** (Monte Carlo standard error),
**ADEMP** (the Morris, White and Crowther reporting framework).

All three specifications share one model,

$$
Y_{it} = \beta_0 + \beta_b\, b_i + \beta_t\, t
       + \beta_D\, X_{it} + \beta_{bm:X}\, b_i X_{it}
       + u_i + \varepsilon_{it},
$$

fitted as an LME with a patient random intercept and `corCAR1`
residual correlation, which lets measurements closer in time be more
similar. In words: the outcome depends on the biomarker, on time, on
drug exposure, and on the product of biomarker and exposure. That
product term is the **estimand**, the specific quantity the trial
is trying to estimate; it says how much the treatment effect changes
per standard deviation of the biomarker. The specifications differ
only in how $X_{it}$ is built and in whether a nuisance term is added.

## 3. The three specifications

### 3.1 Unadjusted (stored code `E1`)

$$X_{it} = D_{it}, \qquad D_{it} \in \{0, 1\}$$

In code, `Sx ~ bm + t + Db + bm:Db`. Carryover is not represented. A
measurement one week after stopping, still carrying about half the
drug effect, is coded exactly like one taken after full clearance. The
estimand is $\beta_{bm:D}$, the coefficient on $b_i D_{it}$.

This is what an analyst gets by default from a standard mixed model
without thinking about carryover, which is why it belongs in the
comparison: it is the benchmark any remedy must beat.

### 3.2 Washout-adjusted (stored code `E4`)

$$X_{it} = D_{it}, \qquad
\text{plus } \beta_{sd}\, t_{sd,it} \text{ as a nuisance main effect}$$

In code, `Sx ~ bm + t + Db + bm:Db + tsd`. Exposure stays binary, but
the model is given the elapsed time since discontinuation as a
covariate, so it can estimate how the average outcome drifts across
the washout without being told how fast the drug clears. This is the
alternative the reference implementation provides
(`implementations/original/R/lme_analysis.R`, `simplecarryover =
TRUE`).

Two properties matter. It consumes no half-life, so it is completely
unaffected by what the analyst assumes about the pharmacokinetics. And
because $t_{sd}$ enters the average-outcome part of the model only,
with no $b_i t_{sd}$ product, it can absorb mean-level carryover but
cannot repair distortion of the *biomarker* effect. The estimand is
$\beta_{bm:D}$, the same quantity Unadjusted targets.

### 3.3 Exposure-weighted (stored code `E2`)

$$X_{it} = D_{bc,it} =
\begin{cases}
1 & \text{if } D_{it} = 1,\\[2pt]
\exp(-\hat{\lambda}\, t_{sd,it}) & \text{if } D_{it} = 0,
\end{cases}
\qquad \hat{\lambda} = \frac{\ln 2}{\hat{t}_{1/2}}
$$

In code, `Sx ~ bm + t + Dbc + bm:Dbc`. Instead of an on/off switch,
exposure is a dimmer starting at 1 and halving every $\hat{t}_{1/2}$
weeks. A measurement one week after stopping, at an assumed half-life
of one week, enters as 0.5 rather than 0. This is the current
pmsimstats-ng default.

It is the only one of the three requiring a commitment to both a decay
shape and a half-life. Because those enter through the variable's
values rather than the formula, it is effectively a different model
for every assumed half-life though written identically each time. Its
estimand is $\beta_{bm:D_{bc}}$, the coefficient on $b_i D_{bc,it}$,
which is *not* the same quantity as $\beta_{bm:D}$.

### 3.4 Retained as a side-check: Lag-adjusted (stored code `E3`)

$$L_{it} = D_{i,t-1}\,(1 - D_{it})$$

In code, `Sx ~ bm + t + Db + bm:Db + L`. The indicator equals 1 at the
first off-drug measurement following an on-drug one, and 0 elsewhere.
This is the classical crossover remedy of Jones and Kenward, and it
occupies the same conceptual slot as Washout-adjusted: a nuisance main
effect, no half-life, estimand $\beta_{bm:D}$. It differs in being
indexed by position in the visit sequence rather than by elapsed time,
so it marks one occasion and is constant thereafter.

It is not an arm of the study, for the reasons in Section 5.1, but its
results are already present in every stored summary and are reported
as a side-check. Section 6 gives them.

### 3.5 The three codings on one patient

The differences are clearest on a single schedule. Both tables assume
a half-life of one week. Verified by direct computation from the
design definitions in `simulation-core.R`.

**OL+BDC, path A.** Drug stopped after week 18.

| Week | On drug | $t_{sd}$ | Unadjusted $D_{it}$ | Washout $t_{sd}$ | Lag flag $L_{it}$ | Exposure-weighted $D_{bc,it}$ |
|---|---|---|---|---|---|---|
| 4 to 18 | yes | 0 | 1 | 0 | 0 | 1 |
| 19 | no | 1 | 0 | 1 | 1 | 0.50 |
| 20 | no | 2 | 0 | 2 | 0 | 0.25 |

At week 19 the patient still has about half the drug effect on board.
Unadjusted records a 0, exactly as for a drug-free patient.
Washout-adjusted also records 0 for exposure, but hands the model the
elapsed time so the average outcome may drift. Only Exposure-weighted
records a value reflecting the residual exposure, 0.50.

**Hybrid, path A.** Stopped after week 10, restarted at 16, stopped
after 16.

| Week | On drug | $t_{sd}$ | Unadjusted $D_{it}$ | Washout $t_{sd}$ | Lag flag $L_{it}$ | Exposure-weighted $D_{bc,it}$ |
|---|---|---|---|---|---|---|
| 4 to 10 | yes | 0 | 1 | 0 | 0 | 1 |
| 11 | no | 1 | 0 | 1 | 1 | 0.50 |
| 12 | no | 2 | 0 | 2 | 0 | 0.25 |
| 16 | yes | 0 | 1 | 0 | 0 | 1 |
| 20 | no | 4 | 0 | 4 | 1 | 0.0625 |

This path shows why the lag flag was replaced. Weeks 11 and 20 both
receive $L_{it} = 1$, each being the first visit after a
discontinuation, so the model gives them one shared coefficient. Yet
week 11 retains about 50% of the drug effect and week 20 about 6%, an
eightfold difference. Washout-adjusted distinguishes them, 1 against
4, and Exposure-weighted distinguishes them further, 0.50 against
0.0625. The lag flag alone cannot tell them apart.

## 4. Why these three

### 4.1 They span the three positions on washout

Write $\phi(t_{sd})$ for the fraction of drug effect still present
$t_{sd}$ weeks after stopping. Every analyst takes a position on
$\phi$, and only three are qualitatively distinct:

- **Deny it.** Assume $\phi \equiv 0$; exposure is strictly on or off.
  This is Unadjusted.
- **Estimate it.** Leave the shape unspecified but let the data say
  how the outcome drifts across the washout. This is
  Washout-adjusted.
- **Assume it.** Commit to a decay curve and rate,
  $\phi = \exp(-\hat\lambda t_{sd})$. This is Exposure-weighted.

The three are therefore one representative of each available stance
toward the unknown washout, ordered by how much is assumed.

A formal relationship makes the ladder continuous. As the assumed
half-life goes to zero the dimmer collapses to the on/off switch,
since $\hat{t}_{1/2} \to 0$ gives $\hat{\lambda} \to \infty$ and hence
$D_{bc,it} \to D_{it}$. **Unadjusted is therefore the limiting case of
Exposure-weighted at an assumed half-life of zero**, not a separate
method outside it. As $\hat{t}_{1/2} \to \infty$ the dimmer never
dims, the on-off contrast vanishes, and nothing is estimable. The
manuscript's half-life-mis-specification sweep is thus a traverse
along a path anchored at Unadjusted.

### 4.2 They separate the two ways carryover does damage

This is the strongest reason for the particular set, because it makes
the comparison a controlled experiment rather than a horse race.

Carryover degrades the interaction test through two mechanisms:

- **Dilution.** The exposure variable is a mis-measured version of
  what the patient received. Contaminated occasions are entered as
  drug-free, so the estimated interaction is pulled toward zero, the
  same phenomenon clinicians know as regression dilution, where a
  single noisy blood-pressure reading understates a true association.
  The signal shrinks.
- **Added noise.** Those occasions also differ from one another in
  ways the model does not describe, some being one week off drug and
  others four. That variation is absorbed into $\varepsilon_{it}$,
  inflating the standard error of the interaction estimate. The noise
  grows.

Power falls on both counts at once. The three arms occupy three cells:

| Specification | Dilution | Added noise |
|---|---|---|
| Unadjusted | not addressed | not addressed |
| Washout-adjusted | not addressed | addressed, graded across the whole washout |
| Exposure-weighted | addressed | addressed |
| *(Lag-adjusted, side-check)* | *not addressed* | *addressed at one occasion only* |

The placements follow from where each term sits. Washout-adjusted adds
$t_{sd}$ to the average-outcome part only, so it absorbs noise while
the interaction stays on $D_{it}$ and dilution is untouched.
Exposure-weighted replaces the variable *inside* the interaction, so
it addresses dilution directly and mops up the noise as a by-product.

Two consequences follow:

1. **Comparing Washout-adjusted with Unadjusted isolates the noise
   mechanism.** The two estimate the same quantity, are fitted to the
   same simulated patients, and differ by exactly one term. Any
   difference is attributable to the noise mechanism alone. Controlled
   contrasts this clean are rare in a simulation study.
2. **What is left over is interpretable.** Whatever advantage
   Exposure-weighted shows beyond that is the value of fixing
   dilution. The design reveals not only which method wins but which
   mechanism matters, which is the more portable finding.

### 4.3 Each is a working default, and the set is now complete

Unadjusted is ordinary applied practice and the default of any
off-the-shelf mixed model. Washout-adjusted and Exposure-weighted are
the two options the reference implementation actually provides
(Section 5.2). The set therefore evaluates that implementation
completely, rather than two-thirds of it plus a substitute.

This involves a deliberate trade, which should be stated plainly. The
crossover tradition's own remedy, the lagged indicator, moves from a
headline arm to a side-check. What is bought is that the study now
tests every option the original authors offered; what is given up is
that the crossover literature is represented by a reported result
rather than by a full arm. Section 6 reports that result, so the
tradition is answered, but with less prominence.

### 4.4 They are equally estimable

All three fit in all three designs at both sample sizes, need no data
beyond what the trial collects, and add at most one term. Verified for
Unadjusted and Exposure-weighted by recomputation from the stored
summaries: across the full 540-cell grid neither produced a missing
power value or estimate, and the same holds for the retained
side-check. Washout-adjusted requires at least one post-discontinuation
occasion, which every path of every design except CO path B provides;
that path never discontinues, so it contributes to the other terms but
not to $\beta_{sd}$. This has not yet been verified by a run.

## 5. Why the middle arm was changed, and what is still missing

### 5.1 Why Lag-adjusted was replaced

Three reasons, in increasing order of weight.

**Provenance.** The lagged indicator is not inherited from the
reference implementation. It was introduced by this project in
`simulation-core.R` following Jones and Kenward. Reporting it as the
comparator while leaving the original authors' own alternative
untested is hard to defend.

**Resolution.** As Section 3.5 shows, the flag marks one occasion and
cannot distinguish a visit one week after stopping from one four weeks
after. In the Hybrid design it assigns a single coefficient to
occasions differing eightfold in residual exposure.

**It may have produced a false negative.** Recomputed from the stored
summaries, Lag-adjusted and Unadjusted are statistically
indistinguishable across all 216 cells in scope (Section 6). The
manuscript currently reads that as evidence that the noise mechanism
is negligible, and therefore that all of Exposure-weighted's advantage
comes from fixing dilution. That inference is only as good as the
arm it rests on. The `quick-sim` outputs show the reference
implementation's $t_{sd}$ arm reaching an *identical point estimate*
to its unadjusted arm while achieving materially higher power, which
is the exact signature of a gain through the noise channel with
dilution untouched. If that reproduces here, the noise mechanism is
not negligible; a binary flag was simply too crude to exploit it, and
the manuscript's mechanistic conclusion would change.

### 5.2 The reference implementation's three options

`implementations/original/R/lme_analysis.R` supports three carryover
representations, assembled at lines 142 to 153 with the exposure
variable built at 108 to 118. None is a lagged indicator.

| Setting | Fitted exposure terms | This study |
|---|---|---|
| `carryover_t1half = 0`, `simplecarryover = FALSE` | $D_{it} + b_i D_{it}$ | Unadjusted |
| `carryover_t1half = 0`, `simplecarryover = TRUE` | $D_{it} + b_i D_{it} + \beta_{sd} t_{sd,it}$ | Washout-adjusted |
| `carryover_t1half > 0` | $D_{bc,it} + b_i D_{bc,it}$ | Exposure-weighted |

The first and third are mutually exclusive with the second by
construction: supplying both a non-zero half-life and
`simplecarryover = TRUE` raises an error (lines 67 to 68).

**A divergence to disclose.** Hendrickson's default `useDE = TRUE`
appends the expectancy covariate $De$ wherever expectancy varies,
which it does in OL+BDC and Hybrid though not CO. The specifications
here omit it: `simulation-core.R` passes `expectancies` to the design
builder, so expectancy shapes the simulated data, but $De$ never
enters any fitted formula. It applies equally to all arms so should
not disturb their ranking, but the fitted model is not Hendrickson's,
and the agreement between this study's 0.860 and Hendrickson's
reported 0.82 has not been checked against that difference.

### 5.3 The cell that remains empty

Replacing the middle arm improves it but does not fill the gap in
Section 4.2. Washout-adjusted still carries the interaction on
$D_{it}$, with $t_{sd}$ entering as a main effect only, so dilution
remains unaddressed by any assumption-free method. The missing
specification is now

$$
Y_{it} = \ldots + \beta_D D_{it} + \beta_{sd} t_{sd,it}
       + \beta_{bm:D}\, b_i D_{it}
       + \beta_{bm:sd}\, b_i t_{sd,it} + u_i + \varepsilon_{it},
$$

in code `Sx ~ bm + t + Db + tsd + bm:Db + bm:tsd`. This is a more
interesting candidate than the $b_i L_{it}$ term it replaces, because
it is essentially a linear approximation to what Exposure-weighted
does parametrically: it lets the biomarker effect decline across the
washout without specifying a rate. Adding it would complete the ladder
of Section 4.1.

It is not proposed for this manuscript. If the grid stays at three,
the Methods section should state that every assumption-free
specification evaluated here enters carryover as a main effect only,
and that the study therefore bounds what such methods can achieve
against dilution rather than settling it.

### 5.4 Candidates set aside

- **Washout deletion**, discarding off-drug measurements within a
  window of stopping. Exposure weighting with crude 0/1 weights and an
  arbitrary cutoff; discards observations the others use.
- **Estimating the half-life from the data.** The principled answer to
  Exposure-weighted's main weakness and the natural next step, but the
  half-life is poorly determined at $N = 35$ from a handful of
  off-drug measurements, introducing a convergence-failure mode the
  others do not have.
- **A full set of lag indicators**, one per off-drug occasion. The
  window is three to four measurements long, so a full set is nearly
  indistinguishable from the time term $t$.

## 6. What is known, and what must still be checked

Epistemic status: the Lag-adjusted and Unadjusted figures below are
recomputed directly from the stored cell-level summaries under
`analysis/data/`. **No Washout-adjusted results exist yet**; the
specification has been added to `fit_spec()` but not run.

**How the comparison is scored.** Following ADEMP: **power**, how
often a real effect is detected, target 0.80; **Type I error**, how
often an effect is declared when none exists, which should sit at
0.05; **bias**, the average gap between estimate and truth;
**coverage**, how often the 95% interval contains the truth, target
0.95; **MSE**, combining bias and spread; and **MCSE**, the
simulation's own margin of error, below which differences should not
be interpreted. A **cell** is one combination of design, sample size,
biomarker strength, decay shape and half-life, run 500 times.

**The side-check result.** Within the 216-cell scope, Lag-adjusted and
Unadjusted are statistically indistinguishable. Mean power differs by
0.003 against an MCSE of about 0.018, only 3 of 216 conditions differ
by more than 0.02, and mean bias (0.373 versus 0.372), MSE (1.885
versus 1.876) and coverage (0.936 versus 0.938) agree to three decimal
places. The same equality holds in every sensitivity block and in the
cluster-robust cells. The classical lagged device therefore provides
no measurable benefit over ignoring carryover altogether, which is a
reportable finding in its own right.

**What the reference implementation's output suggests.** Its $t_{sd}$
arm leads its unadjusted arm by 1.3 to 6.5 percentage points on paired
tests, with an identical point estimate (-1.61 in both) and nearly
identical empirical spread. Identical estimate with higher power means
the gain came through the noise channel, consistent with the
classification in Section 4.2.

**Two checks that must pass before that is believed.** First, the
`quick-sim` run uses a different parameterization and a different code
path, so it establishes only that the term can work, not that it will
here. Second, and more seriously, its $t_{sd}$ arm carried the highest
false-positive rate of its three, up to 0.053 against 0.037 where 0.05
is nominal. If Washout-adjusted proves anticonservative under this
study's parameterization, part or all of its power advantage is not
real. Type I error must therefore be reported alongside power, and no
power gain should be claimed until that check passes.

## 7. Which comparisons this design supports

**Power and Type I error are comparable across all three.** These are
properties of the test, not of the coefficient tested, so they can be
compared freely even though the three do not all estimate the same
quantity.

**Bias, MSE and coverage are comparable only within the Unadjusted and
Washout-adjusted pair.** Those two estimate the coefficient on
$b_i D_{it}$; Exposure-weighted estimates the coefficient on
$b_i D_{bc,it}$. All are scored against the single target
$\theta_{\text{true}} = -3.6$, correct for Exposure-weighted but not
for the others. Exposure-weighted will therefore show near-zero bias
and nominal coverage while the others show substantial bias and
under-coverage, and that gap is largely a statement about which
quantity each targets. Those columns must not be used to rank the
three.

**The Washout-adjusted versus Unadjusted contrast is the controlled
one.** Same estimand, same simulated patients, one term of difference.
Whatever separates them is the noise mechanism, cleanly identified.
This is the most trustworthy comparison in the study.

**Conclusions about assumption-free methods are bounded by the
main-effect form.** Neither Washout-adjusted nor the retained
side-check includes a biomarker interaction with the carryover term,
so a null result for either says that a carryover *main effect* does
not recover the diluted signal. It does not establish that no
assumption-free method can (Section 5.3).

## 8. Summary

The three specifications are Unadjusted, which codes exposure as an
on/off switch; Washout-adjusted, which adds elapsed time since
stopping as a covariate; and Exposure-weighted, which replaces the
switch with a dimmer fading at an assumed half-life. They are the
right three because they occupy the three distinct positions on the
unknown washout curve, because two share an estimand and so isolate
one damage mechanism under otherwise identical conditions while the
third addresses the other, and because the set now covers every option
the reference implementation provides.

The middle arm was changed from the Jones and Kenward lagged indicator
because that indicator is a project invention rather than the
published alternative, because it cannot distinguish occasions
differing eightfold in residual exposure, and because its null result
may have been an artifact of that crudeness rather than evidence about
the mechanism. It is retained as a reported side-check. The set still
leaves one cell empty, at $b_i t_{sd}$, which should be disclaimed in
the Methods if it is not filled.

## 9. Proposed additions to NOTATION.md

**9.1 Add $L_{it}$ and $t_{sd}$ as model terms.** The lagged
just-off-drug indicator has no symbol in Part 1. Proposed:
$L_{it} = D_{i,t-1}(1 - D_{it})$, equal to 1 at the first off-drug
occasion following an on-drug one (code column `L`). $t_{sd}$ is
already listed as a data-generating quantity; the entry should note
that it also serves as an analysis-model covariate.

**9.2 Split the estimand symbol.** Part 1 carries a single entry,
$\beta_{bm:D}$, glossed as the interaction coefficient and annotated
`bm:Dbc`. That conflates two coefficients. Proposed:
$\beta_{bm:D}$ for the coefficient on $b_i D_{it}$ (`bm:Db`), the
estimand of Unadjusted, Washout-adjusted and the side-check; and
$\beta_{bm:D_{bc}}$ for the coefficient on $b_i D_{bc,it}$
(`bm:Dbc`), the estimand of Exposure-weighted. This distinction is
what licenses the controlled comparison of Section 4.2 and why the
bias and coverage columns are comparable within one pair but not
across to the third.

**9.3 Add the main-effect coefficients.** Part 1 has no symbols for
the intercept or non-interaction terms. Proposed: $\beta_0$
intercept, $\beta_b$ biomarker, $\beta_t$ time, $\beta_D$ drug
exposure, $\beta_L$ lagged indicator, $\beta_{sd}$ linear time since
discontinuation. The $\beta_b$ entry should note that it is distinct
from $\beta_{bm}$, reserved for the mean-moderation architecture's
moderation parameter.

**9.4 Register the `E4` stored value.** The Part 2 table lists `spec`
as taking `E1`, `E2`, `E3`. Washout-adjusted adds `E4`. It is a new
code rather than a redefinition of `E3` deliberately: redefining would
silently change the meaning of `E3` in the twenty-one migrated `.rds`
files and in whatever papers 01 and 10 quote, which is the precise
failure mode Part 2 exists to prevent. The Part 3 glossary entry for
analysis specifications needs the fourth line accordingly.

**9.5 A policy tension to resolve before renaming.** Part 2 requires
that a stored value doubling as a reporting label be kept identical to
it, and records `spec` as unmapped. That rule was enforced on
2026-07-31, when twenty-one `.rds` files were migrated from
`A1/A2/A3` to retire a display-time mapping, so renaming reintroduces
the layer that migration removed. The same file supplies the governing
precedent, however: the architectures were renamed to mechanism names
while their stored values were left untouched and the mapping
recorded. Follow that precedent rather than migrating twice. Keep
`E1` to `E4` stored permanently, apply display names through one
shared helper, and mark the Part 2 row as mapped.

---
*Rendered on 2026-08-05 at 10:46 PDT. Source:
`~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/02-carryover-sensitivity/whitepaper-three-analysis-specifications.md`*
