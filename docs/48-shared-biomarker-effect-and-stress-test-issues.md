# Shared Biomarker Effects, Prognostic Biomarkers and the Symmetry Screen: Five Issues from the Hybrid Stress Test {.unlisted .unnumbered}
*2026-10-07 08:54 PDT*

**Author.** pmsimstats team

**Purpose.** The stress test pre-registered in `docs/46` was run on
2026-10-06 and 07: 39 cells, 26,500 simulated trials and 11 analyses,
with 27 failed fits in 291,500. Its results turn on five issues. This
paper explains each from first principles and gives the evidence for
it:

1. what it means for an analysis to *share* the biomarker effect with
   baseline;
2. why a prognostic biomarker breaks every analysis that does so;
3. what *freeing* the baseline biomarker effect does;
4. the pre-registered symmetry screen, and why it failed;
5. the untested combination of the recommended analysis with
   phase-specific variances.

The full results report, against the pre-registration, is reserved for
`docs/47`.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`, `docs/45` and
`docs/46`. Symbols marked *local* are defined here.

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $Y_{it}$ | symptom score (code `Sx`); $t = 0$ is baseline | canonical |
| $B_i$ | biomarker (blood pressure); centered in every analysis (code `bmc`) | canonical |
| $D_{bc,it}$ | exposure-decayed drug indicator (code `Dbc`), 0 at baseline | canonical |
| $\beta_B$ | biomarker main effect in the analysis model | local (`docs/41`) |
| $\beta_{bm:D}$ | biomarker-by-drug interaction, the estimand | canonical |
| $s_0$, $s_{\text{off}}$, $s_{\text{on}}$ | true slope of $Y$ on $B$ at baseline, at post-baseline off-drug visits, at on-drug visits | local |
| $s_{\text{prog}}$ | prognostic slope: the slope of $Y$ on $B$ that does not depend on drug | local |
| $f$ | the share of $s_{\text{prog}}$ that a shared $\beta_B$ absorbs | local |
| $c_{bm}$ | biomarker-BR coupling (the predictive effect) | canonical |
| $\pi$ | power | canonical |

Table: Notation index and glossary of symbols, meanings and status

**Analyses** (`docs/46`, Section 4). All of them use a random intercept,
continuous-time AR(1) residuals (`corCAR1`) and CR2 standard errors
unless stated, with baseline as a response row, linear time, and
$D_{bc,it}$ at an assumed half-life of 0.5 week.

- **S** (the strawman): `Sx ~ bmc + t + Dbc + bmc:Dbc`.
- **ES**: S plus `bmc:bl0`, where `bl0` is 1 at baseline.
- **F1**: S plus `De + bmc:De`, where `De` is the expectancy: 0 at
  baseline, 1 open label, 0.5 blinded.
- **F2**: S plus `bmc:t`.
- **C1**: S with residual variance by phase, `varIdent(~ 1 | phase)`,
  for the phases baseline, open label and blinded. **C2** is C1 with
  model-based inference.
- **A**: unstructured covariance with Kenward-Roger. **EA** is A plus
  `bmc:bl0`.
- **DS** and **DA**: S and A with the half-life chosen by AIC.
- **R**: the published analysis.

**Glossary.**

- **Predictive biomarker.** Changes the effect of the drug.
- **Prognostic biomarker.** Associated with the outcome whatever the
  treatment, for example through the natural course or the placebo
  response.
- **Size.** The false-positive rate of the interaction test; 0.05
  nominal.
- **Size-adjusted power.** Power when each analysis rejects at its own
  5% null quantile, so analyses are compared at equal size.

## 1. The shared biomarker effect

How strongly does the symptom score move with blood pressure? Asked
separately in three groups of rows, the question gives three slopes:

| Rows | True slope of $Y$ on $B$ |
|---|---|
| baseline (week 0, before any treatment) | $s_0$ |
| post-baseline visits off drug | $s_{\text{off}}$ |
| post-baseline visits on drug | $s_{\text{on}}$ |

Table: True slope of $Y$ on $B$ in three groups of rows

The estimand is the extra biomarker effect on drug,
$\beta_{bm:D} = s_{\text{on}} - s_{\text{off}}$.

The analyses that keep baseline as an outcome row coded off drug
(S, A, C1, C2, DS, DA, F2, R) have only two slope parameters:

- $\beta_B$ (`bmc`) is the slope in *every* row coded off drug: the
  baseline row and the post-baseline off-drug visits together;
- $\beta_{bm:D}$ (`bmc:Dbc`) is the additional slope on drug.

So the model imposes $s_0 = s_{\text{off}}$: a single $\beta_B$
describes baseline and the later off-drug visits. This is what *sharing
the biomarker effect with baseline* means. The condition that makes it
correct, $s_0 = s_{\text{off}}$, is the constancy assumption.

**When the assumption holds,** for a purely predictive biomarker,
$s_0 = s_{\text{off}} = 0$ and only $s_{\text{on}}$ differs from zero.
The constraint is then correct, and the baseline row adds a clean
off-drug observation that sharpens $\beta_B$, and through it
$\beta_{bm:D}$. That free precision is why the shared-slope analyses
were the most powerful under AR(1) data (verified):

| Size-adjusted power, AR(1) construct, $c_{bm} = 0.30$ | |
|---|---|
| S (shared) | 0.514 |
| ES (freed) | 0.364 |

Table: Size-adjusted power of shared and freed analyses, AR(1) construct

The mechanism, the extra information the baseline row contributes to
$\beta_B$, was isolated in `docs/45`, Section 6, by freeing that one
term.

## 2. Why a prognostic biomarker breaks it

Suppose blood pressure predicts the natural course of symptoms. For
example, participants with higher pressure improve more over time,
whether or not they take the drug. With no true interaction:

- $s_0 = 0$: no course has happened yet at baseline;
- $s_{\text{off}} = s_{\text{on}} = s_{\text{prog}}$;
- the true interaction is $s_{\text{on}} - s_{\text{off}} = 0$.

The shared-slope model must fit one $\beta_B$ to the baseline rows
(true slope 0) and to the off-drug visits (true slope
$s_{\text{prog}}$). The fit settles between them,
$\hat\beta_B \approx f\, s_{\text{prog}}$ with $0 < f < 1$. The on-drug
rows have slope $s_{\text{prog}}$, and the only parameter left to
account for the gap is the interaction:

$$
\hat\beta_{bm:D} \approx s_{\text{prog}} - f\, s_{\text{prog}}
  = (1 - f)\, s_{\text{prog}} \neq 0 .
$$

The model therefore reports an interaction that does not exist. The
error is a fixed bias, not noise, so the false-positive rate grows with
the sample size.

How much leaks depends on how much weight the baseline row receives
when $\beta_B$ is estimated. The more weight it gets, the smaller $f$
and the larger the leak.

| Analysis | Weight on baseline | Null rejection: TV, constant | TV, growing | PB |
|---|---|---|---|---|
| S | moderate | 0.099 | 0.049 | 0.100 |
| DS (S with AIC) | moderate | 0.115 | 0.058 | 0.112 |
| C1 (S, phase variances) | larger | 0.382 | 0.173 | 0.205 |
| A (unstructured, KR) | heavy | 0.416 | 0.195 | 0.222 |
| R (published) | moderate | 0.107 | 0.053 | 0.159 |

Table: Baseline weight and null rejection under prognostic biomarker cells

*1,000 simulated trials per cell (verified). TV: blood pressure
correlated 0.3 with the natural-course component, at every visit or
growing along its curve. PB: correlated 0.3 with the placebo
component.*

The rates are verified; the weighting account is inferred.

- **The unstructured model fails worst.** It weights the baseline row
  most heavily, because baseline is the off-drug observation least
  contaminated by treatment.
- **Phase-specific variances make the leak worse.** The baseline
  variance is the smallest, so the phase-variance model up-weights the
  baseline row.
- **A placebo-prognostic biomarker** breaks the constraint the same way:
  the placebo response is 0 at baseline and nonzero afterwards, so again
  $s_0 \ne s_{\text{off}}$.
- **A natural-course association that grows over time** leaks least into
  S, because it is small early, which is also where most on-drug visits
  are.

## 3. Freeing the baseline biomarker effect

The direct remedy gives baseline its own slope:

`Sx ~ bmc + t + Dbc + bmc:Dbc + bmc:bl0`

| Rows | Slope in the model |
|---|---|
| baseline | $\beta_B$ + coefficient of `bmc:bl0` (whatever the data say) |
| post-baseline off drug | $\beta_B$ |
| on drug | $\beta_B + \beta_{bm:D}$ |

Table: Slope of the biomarker in the model for each group of rows

The interaction is now identified from post-baseline visits alone, and
nothing that happens at baseline can leak into it. This is ES. It is
equivalent to analyzing baseline as a covariate (`docs/41`, Section 7,
confirmed in `docs/45`, Section 6).

**F1 frees baseline implicitly.** It adds the expectancy and its
interaction with the biomarker, `De + bmc:De`. Because the expectancy
is 0 at baseline and positive afterwards, the biomarker slope at
baseline is no longer tied to the post-baseline slope. In addition, the
slope may scale with the expectancy, which matches the placebo-prognostic
case directly.

| Analysis | Pooled null, 11 standard cells | TV, constant | TV, growing | PB |
|---|---|---|---|---|
| ES | 0.048 | 0.049 | 0.063 | 0.049 |
| EA | 0.052 | 0.049 | 0.057 | 0.050 |
| F1 | 0.051 | 0.047 | 0.046 | 0.040 |
| F2 (`bmc:t` only) | 0.049 | 0.090 | 0.048 | 0.098 |

Table: Null rejection of analyses that free the baseline biomarker effect

*Verified.*

F2 shows that modeling the time course alone is not enough. It repairs
only the growing natural-course case, because it does not free the
baseline slope.

**The cost.** When the biomarker is purely predictive, freeing the slope
throws away the information the baseline row contributed. Under AR(1)
data this costs about a third of the power (0.364 against 0.514 above).
Under compound symmetry the participant effect already supplies that
information, and the cost largely disappears (verified):

| Size-adjusted power | S | ES | F1 |
|---|---|---|---|
| CS construct, $c_{bm} = 0.25$ | 0.554 | 0.662 | 0.522 |

Table: Size-adjusted power of S, ES and F1 in the CS construct

**Choosing between ES and F1.** Both pass every size screen. They
differ in power by condition (verified):

| Size-adjusted power, alternative cells | ES | F1 |
|---|---|---|
| maximin shortfall, non-prognostic cells | 0.046 | 0.140 |
| TV, constant prognosis, $c_{bm} = 0.30$ | 0.370 | 0.482 |
| TV, growing prognosis, $c_{bm} = 0.30$ | 0.192 | 0.378 |
| PB prognosis, $c_{bm} = 0.30$ | 0.498 | 0.446 |
| maximin shortfall, all cells | 0.186 | 0.140 |

Table: Size-adjusted power and maximin shortfall of ES and F1 by cell

- When the biomarker is purely predictive, ES is slightly better.
- When it is prognostic, F1 keeps its power across all three forms;
  ES loses most of it when the association grows over time.

The ranking therefore depends on whether the prognostic alternative
cells count. The pre-registration does not settle this; `docs/47` will
report both versions. `docs/46` recommends F1 because it is robust in
size and in power.

## 4. The symmetry screen

**What was pre-registered.** One component of sensitivity (`docs/46`,
Section 3) required that power at $-c_{bm}$ be within 2 Monte Carlo SE
of power at $+c_{bm}$. It was tested at $\pm 0.30$ in the AR(1)
construct and $\pm 0.20$ in the compound-symmetry construct. An analysis
failing either test could not be selected.

**Why power is exactly symmetric** (derived). In these cells the
biomarker is correlated only with the drug-response component. Let
$B^* = 2\mu_B - B$, the biomarker reflected about its mean. $B^*$ has
the same mean and variance as $B$, and its covariance with every other
variable has the opposite sign. So the joint distribution of
$(B^*, \text{everything else})$ under $+c_{bm}$ is exactly that of
$(B, \text{everything else})$ under $-c_{bm}$.

Every analysis uses the centered biomarker linearly, and negating that
column negates $\hat\beta_B$ and $\hat\beta_{bm:D}$ while leaving every
fitted value, variance estimate, standard error, degree of freedom and
p-value unchanged. The two-sided rejection rate at $-c_{bm}$ therefore
equals the rate at $+c_{bm}$. The size-adjusted rates are equal too,
because both cells use the same null quantile.

Any observed difference is Monte Carlo error.

**What happened.** The $-c_{bm}$ and $+c_{bm}$ cells are separate,
independently seeded simulations. The analyses within a cell are fitted
to the same trials and are therefore strongly correlated.

- In both constructs, power came out higher at $+c_{bm}$ for every
  analysis, by 0.004 to 0.084 (verified).
- At the threshold, about 0.06, 9 of the 22 tests failed: S, DS and F2
  in the AR(1) construct; A, C1, EA, ES and R in the
  compound-symmetry construct.
- The screen thereby removed the strawman and ES for no real reason.

**Rerun with fresh seeds** (`13-symmetry-check.R`, 1,000 trials per
cell, verified):

| Construct | Analysis | Power at $-c_{bm}$ | Power at $+c_{bm}$ | Difference / SE |
|---|---|---|---|---|
| CS, $\pm 0.20$ | S | 0.363 | 0.367 | 0.2 |
| CS, $\pm 0.20$ | ES | 0.418 | 0.416 | 0.1 |
| AR(1), $\pm 0.30$ | S | 0.447 | 0.481 | 1.5 |
| AR(1), $\pm 0.30$ | ES | 0.368 | 0.391 | 1.1 |

Table: Power at $-c_{bm}$ and $+c_{bm}$ in the symmetry rerun

The estimates are mirror images (for S in the AR(1) construct, $+0.155$
and $-0.158$). The compound-symmetry cells agree exactly. The AR(1)
cells again lean toward $+c_{bm}$, but within chance. Since the
identity is exact, no amount of leaning can make it real.

**Why the screen was a design error.** It tested, with noisy and
correlated comparisons, a property that holds by construction.

- **The false-failure rate was high.** Each test had about a 5% chance
  of false failure on its own. Because one simulated data set served
  all analyses, a single unlucky draw failed several analyses together.
- **The screen could only remove good analyses by chance.** It could not
  detect any real insensitivity, because none can exist in these cells.

There are two correct designs. Either state the symmetry as a property
and test only one sign, or generate the $-c_{bm}$ trials from the same
random draws as the $+c_{bm}$ trials with the biomarker reflected. The
second makes any asymmetry in the code visible and asymmetry from
sampling impossible.

**Effect on the selection.**

| Rule | Survivors | Selected |
|---|---|---|
| As pre-registered | C2, F1 | C2 by maximin, but flagged by the prognostic screen; only F1 passes every screen |
| Symmetry screen dropped (a deviation) | S, C1, C2, DS, ES, EA, F1, F2 | C1 by maximin (0.034), but flagged; among the unflagged, ES on non-prognostic cells, F1 on all cells |

Table: Survivors and selected analysis with and without the symmetry screen

Either way, the analyses that are both nominal and robust to a
prognostic biomarker are ES, EA and F1, and the recommendation of F1
does not depend on the faulty screen. `docs/47` will present the
as-registered outcome first and the deviation second.

## 5. Combining F1 with phase-specific variances (not yet tested)

**The idea.** C1 and C2 were the most powerful analyses in nearly every
non-prognostic cell. For example, at $c_{bm} = 0.30$ in the AR(1)
construct, C2 reached 0.606 and C1 0.578, against 0.514 for S and 0.364
for F1. F1 is the most robust. The combination is C1 with F1's terms:

```{.r}
nlme::lme(Sx ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De,
          random = ~ 1 | ptID,
          correlation = nlme::corCAR1(form = ~ t | ptID),
          weights = nlme::varIdent(form = ~ 1 | phase))
# CR2 standard errors as in S
```

**Why it might combine the two strengths** (inferred):

- **The power of C1** probably comes from weighting. The implied variance
  differs by phase (`docs/40`: baseline 342, open label 710, blinded 599
  in the published construct). A model that lets the variance differ
  weights each visit by its actual precision.
- **The robustness of F1** comes from its mean model. `bmc:De` frees the
  baseline slope and absorbs an expectancy-scaled prognostic slope.
- **The two operate on different parts of the model,** variance and
  mean, so in principle neither undoes the other.

**Why it might not:**

- **The weighting that gives C1 its power is the weighting that made it
  fail worst under a prognostic biomarker** (0.382 against 0.099 for S).
  It up-weights the baseline row. If F1's terms free the baseline slope,
  that extra weight should become harmless, because it now feeds only
  the baseline-specific slope, not the interaction. But this is the
  untested crux.
- **Freeing the baseline slope removes the baseline information** that
  helped the shared-slope analyses. The power gain from phase variances
  may shrink once baseline no longer informs $\beta_B$. In the
  non-prognostic cells, the phase variances might then add little over
  F1 alone.
- **More variance parameters** (two more) may cost something at
  $N = 35$, where C1 was already the highest of the CR2 analyses
  (0.069, within the cell limit).

**Status.** Untested. It was not among the pre-registered analyses, so
any result is post hoc. It can motivate a recommendation, but it cannot
replace the pre-registered one without a confirmatory run.

**The test.** One arm on the same 39 cells and the same trials, run as a
post hoc extension of `11-strawman-stress-test.R`, takes about 3 to 4
hours. It succeeds if two conditions hold:

- it is nominal in all 11 standard null cells and in all three
  prognostic cells;
- its size-adjusted power exceeds F1's by a margin worth its
  complexity, say 0.05 in the AR(1) cells.

## 6. Summary

- The shared-slope analyses, the strawman, the unstructured model and
  the published analysis among them, gain power by assuming the
  biomarker's off-drug association is the same at baseline and
  afterwards.
- A prognostic biomarker violates that assumption and produces false
  interactions. These reach 0.10 for the strawman and 0.42 for the
  unstructured model.
- Freeing the baseline slope (ES, or F1 through the expectancy term)
  removes the problem. Under AR(1) data it costs about a third of the
  power; under compound symmetry it costs little.
- The symmetry screen was a design error: it tested an exact identity
  with noisy comparisons. Dropping it does not change which analyses
  are both nominal and robust.
- Phase-specific variances combined with F1's terms is the most
  promising untested option.

## 7. Reproducibility

```bash
Rscript analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R \
  --cores 8 --chunk 25
Rscript analysis/scripts/quick-sim/carryover-closed-form/12-strawman-stress-decision.R
Rscript analysis/scripts/quick-sim/carryover-closed-form/13-symmetry-check.R \
  --reps 1000 --cores 8
```

Outputs are in
`analysis/data/quick-sim/carryover-closed-form/strawman-stress/`:

- `summary-full.csv` and `replicates-full.rds`;
- the `decision-*.csv` tables;
- `symmetry/summary-reps1000.csv`.

## 8. References

1. pmsimstats team. `docs/40-hendrickson-hybrid-design-primer.md`,
   `docs/41-baseline-as-covariate.md`,
   `docs/45-correctly-sized-hybrid-analysis.md`,
   `docs/46-strawman-stress-test-preregistration.md`.
2. Liu GF, Lu K, Mogg R, Mallick M, Mehrotra DV. Should baseline be a
   covariate or dependent variable in analyses of change from baseline
   in clinical trials? *Statistics in Medicine* 2009; 28(20):2509-2530.
   doi:10.1002/sim.3639.
3. Bell RM, McCaffrey DF. Bias reduction in standard errors for linear
   regression with multi-stage samples. *Survey Methodology* 2002;
   28(2):169-181.
