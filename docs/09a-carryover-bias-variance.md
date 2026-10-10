# Carryover Effects and the Bias-Variance Tradeoff in Model Misspecification {.unlisted .unnumbered}
*2026-02-18*

## A White Paper on Non-Monotonic Power in Precision Medicine Trials {.unlisted .unnumbered}

**Author.** pmsimstats team

**Purpose.** This paper, written in the pmsimstats2025 era (source
`white_paper_carryover_bias_variance.md` of that project, now archived
under `archive/pmsimstats2025-2026-04-23/`), examines how the power of
a biomarker-treatment interaction test in a simplified open-label then
discontinuation design varies with an unmodeled carryover half-life,
and proposes a bias-variance explanation for the pattern.

**Status note (2026-10-07).** This is a historical document; a review
on 2026-10-07 corrected its description of the simulation against the
archived scripts, corrected its references, and found that several of
its central claims are not supported by its own results. In brief: (i)
the simulation is a single-path mean-moderation model with linear
response growth, not the Hendrickson et al. (2020) model, and its
"without carryover" analysis is misspecified even at $t_{1/2} = 0$; (ii)
the simulated power does not show a clear drop-then-recovery pattern
(Section 5); (iii) the derivation of Section 3.5 does not produce an
interior maximum as originally stated; and (iv) the "precipitous
decline" in the published Hendrickson results has since been traced to
the step rule for the biomarker-response coupling under covariance
moderation, not to a bias-variance tradeoff (`docs/36`, Sections 5.3
to 5.5; `docs/37`). Corrections are marked in place with the same
date.

---

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Symbols follow `analysis/report/NOTATION.md` where it defines one.
Code names (`weeks_on_drug`, `bm_centered`) refer to the archived
scripts and appear only in code font.

| Symbol | Meaning | Status |
|---|---|---|
| $Y_{it}$ | observed symptom score of participant $i$ at week $t$ | canonical |
| $BR_{it}$ | biological response component | canonical |
| $\varepsilon_{it}$ | observation-level residual, $N(0, 2^2)$ in the simulation | canonical |
| $b_i$ | standardized biomarker (code `bm_centered`) | canonical |
| $\beta_{bm}$ | mean-moderation multiplier, 0.3 in the simulation | canonical |
| $\beta_{bm:D}$ | biomarker-by-exposure interaction coefficient of the analysis model; true value $-a\beta_{bm} = -0.15$ | canonical |
| $N$ | total number of participants (70; one path) | canonical |
| $t_{1/2}$ | carryover half-life (weeks) | canonical |
| $t_{od}$, $t_{sd}$ | cumulative weeks on drug (code `weeks_on_drug`); time since discontinuation | canonical |
| $t_{od}^{\text{disc}}$ | weeks on drug at discontinuation (8 in the simulation) | local |
| $g(t_{sd})$ | carryover retention, $2^{-t_{sd}/t_{1/2}}$ for $t_{1/2} > 0$, 0 for $t_{1/2} = 0$ | local |
| $a$ | response growth per week on drug (code `BR_rate`, 0.5) | local |
| $m$, $T$ | constant off-drug prediction of the analysis model; length of the off-drug span | local |
| $V(t_{1/2})$ | integrated squared misspecification over the off-drug span | local |
| $\tau_1$, $\tilde\tau_2$, $\omega_0$, $\omega_1$ | period treatment effects and carryover effects in the crossover estimator of Section 2.3 | local |
| $c_{bm}$ | covariance-moderation correlation of the Hendrickson model (Section 5.5 only) | canonical |
| $\pi$ | power, $\Pr(p < \alpha)$ | canonical |
| $\alpha$ | nominal significance level, 0.05 | canonical |
| $n_{\text{reps}}$ | Monte Carlo replicates per cell (100) | canonical |

Table: Notation index for this paper

**Glossary.**

- *OL+BDC*: open-label treatment followed by blinded discontinuation.
  In this paper's simulation every participant is on drug in weeks 1
  to 8 and off drug thereafter; there is no randomized second path.
- *Carryover*: persistence of the drug effect after discontinuation,
  here exponential with half-life $t_{1/2}$.
- *Without / with carryover modeling*: analysis with exposure
  `weeks_on_drug` (constant after discontinuation) against analysis
  with an exponentially decayed exposure at an assumed half-life fixed
  at 0.2 weeks.
- *Bias %*: $100\,(\bar{\hat\beta}_{bm:D} - \beta_{bm:D})/\lvert\beta_{bm:D}\rvert$;
  since $\beta_{bm:D} < 0$, a positive value denotes attenuation toward
  zero.
- *MCSE*: Monte Carlo standard error, $\sqrt{\pi(1 - \pi)/n_{\text{reps}}}$
  for power.

## Executive Summary

This white paper examines a counterintuitive phenomenon observed in simulations of
Open-Label plus Blinded Discontinuation (OL+BDC) trial designs: when carryover
effects are not modeled in the analysis, statistical power exhibits non-monotonic
behavior as carryover half-life increases. Specifically, power drops sharply from
zero carryover to short carryover ($t_{1/2}$ = 0.1-0.5 weeks), then partially recovers
as carryover extends further ($t_{1/2}$ = 1-2 weeks).

This phenomenon reflects a fundamental bias-variance tradeoff arising from model
misspecification and has important implications for trial design and analysis.

*Correction (2026-10-07).* The paper's own simulation (Section 5.1)
does not reproduce this pattern clearly. Without carryover modeling,
power is 0.06 at $t_{1/2} = 0$, 0.00-0.02 at 0.05-0.5 weeks and rises
to 0.36 at 2 weeks (simulated, 100 replicates per cell; MCSE up to
0.05). The drop from 0.06 to 0.01 is about two MCSEs of a difference,
and at $t_{1/2} = 0$ the mean estimate has the wrong sign, so the
evidence is for a power that increases with half-life beyond about 0.2
weeks rather than for a drop followed by recovery. The pattern quoted
in this summary comes from the unsourced table of Section 1.1.

---

## 1. Introduction

### 1.1 The Observed Phenomenon

In simulations comparing trial designs for precision medicine validation
(Hendrickson et al., 2020), we observed the following pattern for OL+BDC designs
analyzed WITHOUT explicit carryover modeling:

| Carryover $t_{1/2}$ | Statistical Power |
|--------------|-------------------|
| 0 weeks      | ~85% (highest)    |
| 0.1 weeks    | ~60% (sharp drop) |
| 0.2 weeks    | ~55% (near minimum) |
| 0.5 weeks    | ~58% (slight recovery) |
| 1.0 weeks    | ~62% (continued recovery) |
| 2.0 weeks    | ~68% (partial recovery) |

Table: OL+BDC power by carryover half-life without explicit carryover modeling (approximate values; source not recorded)

*Note (2026-10-07).* The source of these approximate values is not
recorded, and a search of the archived pmsimstats2025 analysis scripts
for this half-life grid found none; they are unverified. The half-life grid
matches neither the Hendrickson et al. (2020) grid (0, 0.1, 0.2 weeks)
nor exactly the grid of Section 5, and the baseline of about 85%
differs from the 6% simulated at $t_{1/2} = 0$ in Section 5.1.

This non-monotonic pattern contradicts the naive expectation that increasing
carryover should monotonically degrade power when unmodeled.

### 1.2 Significance

Hendrickson et al. (2020) noted that "the presence of even a short (0.1 weeks)
carryover component resulted in a **precipitous decline in power**" in both N-of-1
and OL+BDC designs (the full sentence reads "in both the N-of-1 and, to a
lesser but still very significant degree, the open label + blinded
discontinuation design"; quotation verified against the PubMed Central
full text). The magnitude of this drop underscores the importance of
carryover as a design parameter. The published grid stops at 0.2 weeks,
so any recovery at longer half-lives lies outside it.

Understanding this phenomenon is essential for:

1. Appropriate trial design selection
2. Deciding whether to model carryover in the analysis
3. Setting realistic power expectations across carryover scenarios

---

## 2. Background: Carryover in Crossover and N-of-1 Designs

### 2.1 Definition of Carryover

Carryover (or residual) effects occur when the pharmacological or biological
effect of a treatment persists into subsequent periods after treatment
discontinuation. In the context of our simulations, carryover follows an
exponential decay model:

$$
g(t_{sd}) = \left(\tfrac{1}{2}\right)^{t_{sd}/t_{1/2}},
$$

where $t_{sd}$ is time since discontinuation and $t_{1/2}$ is the
half-life; $g = 0$ off drug when $t_{1/2} = 0$.

### 2.2 Literature Guidance

The crossover trial literature provides extensive guidance on carryover
(Wellek & Blettner, 2012; Lim & In, 2021):

- **Washout periods**: Standard practice sets washout at 3-5× the
  pharmacokinetic half-life (attributed; not checked against the full
  texts of the two reviews)
- **Statistical tests**: Grizzle's (1965) two-stage procedure tests for
  carryover before analyzing treatment effects, though Freeman (1989)
  showed that it inflates the Type I error and judged it too misleading
  for practical use
- **Model-based approaches**: Mixed effects models can incorporate carryover
  terms, but are sensitive to misspecification

Recent work on behavioral carryover (Shi & Ye, 2024) distinguishes between:

- **Biological carryover**: Pharmacokinetic persistence (addressable via washout)
- **Behavioral carryover**: Treatment-induced behavior changes that persist
  regardless of washout

### 2.3 Mathematical Framework in the Literature

Shi and Ye (2024) formalize carryover in potential outcomes notation
(their symbols are renamed here to avoid collision with the project's
canonical $\lambda$ and $\theta$):

- $\omega_0$: carryover effect when the previous treatment was control
- $\omega_1$: carryover effect when the previous treatment was active

With $\tau_1$ and $\tilde\tau_2$ the treatment effects in the two
periods, the basic crossover estimator $\hat\tau$ has expectation

$$
E(\hat\tau) = \tfrac{1}{2}\bigl(\tau_1 + \tilde\tau_2 - \omega_1 + \omega_0\bigr).
$$

This shows the estimator is biased by half the difference of the two
carryover effects, $-(\omega_1 - \omega_0)/2$. The direction and
magnitude of bias depend on the sign and size of this difference; with
no carryover from control ($\omega_0 = 0$) the bias is $-\omega_1/2$,
toward zero when the carryover has the sign of the treatment effect.
(Unverified in detail: the published abstract confirms that under a
sign condition on the carryover effect the basic estimator
underestimates the treatment effect, but this expression was not
checked against the full text.)

(Corrected 2026-10-08: this section previously wrote the expectation as
$\tfrac12(\tau_1 + \tilde\tau_2 - \omega_0 - \omega_1)$, with the
carryover terms summed, and described the bias as driven by the sum
$\omega_0 + \omega_1$. Under the additive two-period model, with
$\hat\tau = \tfrac12[(\bar Y_{AB,1} - \bar Y_{AB,2}) - (\bar Y_{BA,1} -
\bar Y_{BA,2})]$, the carryover from active treatment enters the
period-2 control mean of the AB sequence and the carryover from control
enters the period-2 active mean of the BA sequence, so the two terms
enter with opposite signs (derived; checked numerically). The summed
form agrees only if $\omega_0$ is defined with the opposite sign, and
both forms agree when $\omega_0 = 0$. The Shi and Ye (2024) convention
could not be checked because the full text was not available. See
`analysis/report/compendium-mathematics.tex`, subsection "Carryover
bias of the classical two-period crossover", Equation
`eq:inva-crossover`.)

---

## 3. Mathematical Framework for the Non-Monotonic Phenomenon

### 3.1 The OL+BDC Design Structure

In an OL+BDC design:

- **Open-label phase** (weeks 1-8): All participants receive active treatment
- **Blinded discontinuation** (weeks 9+): Participants randomized to continue
  or discontinue

*Correction (2026-10-07, inspected `simulation_carryover_spectrum.R`).*
The simulation does not randomize at discontinuation: every
participant is on drug in weeks 1 to 8 and off drug from week 9, with
measurements at weeks 1 to 12, 16 and 20. It is therefore a single-path
design with no concurrent on-drug comparison during the discontinuation
phase.

The analysis predictor in a simple model (without carryover) is `weeks_on_drug`,
which accumulates during treatment but remains constant after discontinuation.

### 3.2 True Data-Generating Process

Let $BR_{it}$ denote the true biological response of participant $i$ at
week $t$. *Correction (2026-10-07).* The original text called this
"the Hendrickson model", used a Gompertz growth function, and applied
the moderation factor twice off drug. The archived script
(`simulation_carryover_spectrum.R`, inspected) instead uses a
mean-moderation model with linear growth, which differs from the
covariance-moderation model of Hendrickson et al. (2020):

**During treatment:**

$$
BR_{it} = a\,t_{od,it}\,(1 + \beta_{bm} b_i),
$$

**After discontinuation ($t_{1/2} > 0$):**

$$
BR_{it} = a\,t_{od}^{\text{disc}}\,(1 + \beta_{bm} b_i)\,g(t_{sd,it}),
$$

and $BR_{it} = 0$ off drug when $t_{1/2} = 0$. The outcome is
$Y_{it} = 50 - BR_{it} + \varepsilon_{it}$, with no random intercept
or serial correlation in the data-generating process, where:

- $a = 0.5$ is the response growth per week on drug
- $\beta_{bm} = 0.3$ is the biomarker moderation multiplier (the
  parameter of interest), applied to the standardized biomarker $b_i$
- $t_{od}^{\text{disc}} = 8$ is the number of weeks on drug at
  discontinuation
- $t_{sd}$ is time since discontinuation, counted from the first
  off-drug visit (week 9, where $t_{sd} = 0$ and $g = 1$ for any
  $t_{1/2} > 0$)
- $t_{1/2}$ is the carryover half-life
- $\varepsilon_{it} \sim N(0, 2^2)$ independently

### 3.3 Analysis Model (Without Carryover)

The analysis model assumes:

```
symptoms ~ biomarker + weeks_on_drug + biomarker:weeks_on_drug + (1|participant)
```

The key interaction term `biomarker:weeks_on_drug` is used to detect biomarker
moderation of treatment response. In the script the model is fitted
with `lme4::lmer` by maximum likelihood, with the standardized
biomarker `bm_centered` in place of `biomarker`, and the interaction
is tested with a $t$ reference distribution on $N - 4 = 66$ degrees of
freedom. Its target is $\beta_{bm:D} = -a\beta_{bm} = -0.15$, the
on-drug slope. The "with carryover" analysis of Section 5 replaces
`weeks_on_drug` off drug by $t_{od}^{\text{disc}}\,2^{-t_{sd}/0.2}$,
that is, it assumes a half-life of 0.2 weeks whatever the true value.

### 3.4 Source of the Bias-Variance Tradeoff

The non-monotonic power arises from competing effects:

#### Component 1: Signal Attenuation (Monotonically Increasing with $t_{1/2}$)

As $t_{1/2}$ increases, the true effect during off-drug periods increasingly differs
from zero. But the analysis predictor (`weeks_on_drug`) remains constant at its
discontinuation value. This creates:

- **At $t_{1/2}$ = 0**: Perfect on/off separation. Off-drug has zero true effect and
  predictor is irrelevant for those periods.
- **At $t_{1/2}$ > 0**: The "off-drug" periods have residual effect that the
  predictor cannot capture.

*Correction (2026-10-07, derived from the model of Section 3.2).* The
first bullet does not hold for the predictor actually used. Because
`weeks_on_drug` stays at 8 after discontinuation, at $t_{1/2} = 0$ the
model predicts the full on-drug interaction at every off-drug visit
while the true effect is zero: this is the largest mismatch between
predictor and truth, not perfect separation. The mismatch shrinks
monotonically as $t_{1/2}$ grows and the true off-drug effect
approaches the constant predictor. This is consistent with the
simulated estimates (Section 5.1), which are attenuated by 92-96% at short
half-lives and less at long ones, and with the wrong-signed mean
estimate at $t_{1/2} = 0$.

#### Component 2: Unexplained Variance (Non-Monotonic with $t_{1/2}$)

The original argument turns on how unexplained variance in off-drug observations
varies with $t_{1/2}$:

**At $t_{1/2} \approx$ 0.1-0.5 weeks (short carryover):**
- True effect decays rapidly during off-drug (at $t_{1/2}$ = 0.1, from 100% at
  the first off-drug visit to 0.1% one week later; corrected 2026-10-07
  from "12% in 3 weeks", which is the $t_{1/2} = 1$ value)
- This creates high WITHIN-PERIOD variance in response
- The predictor is constant, but the outcome varies substantially
- Maximum mismatch → Maximum residual variance → **Minimum power**

**At $t_{1/2} \approx$ 2+ weeks (long carryover):**
- True effect decays slowly (at $t_{1/2}$ = 2, to 71% after one week and 35%
  after three; corrected 2026-10-07 from "70% in 3 weeks")
- Off-drug response remains relatively stable
- The predictor is constant, AND the outcome is roughly constant
- Reduced mismatch → Reduced residual variance → **Power recovery**

### 3.5 Formal Expression

The unexplained variance in the model residuals can be described
heuristically as the sum of within-participant noise, between-participant
variation, and a misspecification term, the mean squared difference
between the true response and the model's prediction. For off-drug
observations the true response is $A\,g(t_{sd})$, with
$A = a\,t_{od}^{\text{disc}}(1 + \beta_{bm} b_i)$ the response at
discontinuation, while the analysis model, whose predictor is constant
off drug, predicts a constant $m$. Over an off-drug span of length $T$
the integrated squared misspecification is

$$
V(t_{1/2}) = \int_0^{T} \bigl[A\,2^{-s/t_{1/2}} - m\bigr]^2 \, ds .
$$

*Correction (2026-10-07, derived).* The original text stated that this
integral has an interior maximum in $t_{1/2}$. Whether it does depends
on $m$:

- **$m$ fixed at the on-drug level, $m = A$.** This is what a constant
  `weeks_on_drug` predictor with a common slope implies. For every $s$,
  $2^{-s/t_{1/2}}$ increases with $t_{1/2}$, so
  $(A\,2^{-s/t_{1/2}} - A)^2$ decreases, and $V$ falls monotonically
  from $A^2 T$ as $t_{1/2} \to 0$ to 0 as $t_{1/2} \to \infty$. There
  is no interior maximum.
- **$m$ the least-squares constant over the span.** Then
  $V = T\,\mathrm{Var}_s\bigl(A\,2^{-s/t_{1/2}}\bigr)$, which is 0 in
  both limits (the trajectory is constant at 0 or at $A$) and positive
  in between, so an interior maximum exists.

(The original limits were also misstated: at $t_{1/2} \to 0$ the
integral is $m^2 T$, not $m^2$.) The fitted model is neither case
exactly, because one interaction slope serves both phases, but the
simulated estimates of Section 5.1, whose attenuation is largest at
$t_{1/2} = 0$ (where the estimate changes sign), roughly flat at 95%
up to 0.2 weeks and shrinking thereafter, behave like the first case. The bias-variance account
of a minimum at intermediate half-life is therefore not supported by
this derivation.

---

## 4. Intuitive Explanation

Consider an analogy: imagine trying to predict a person's location using only
their home address as a predictor.

### Scenario A: Person stays at home ($t_{1/2}$ = 0)
- Predictor: Home address
- True location: Always at home
- Result: **Perfect prediction**

### Scenario B: Person commutes daily ($t_{1/2}$ = 0.5)
- Predictor: Home address
- True location: Varies between home, transit, and work
- Result: **Poor prediction**: high unexplained variance

### Scenario C: Person works from home ($t_{1/2}$ = ∞)
- Predictor: Home address
- True location: Almost always at home
- Result: **Good prediction**: effect persists, matching the constant predictor

The intermediate case (Scenario B) creates the most unexplained variance
because the true outcome varies while the predictor remains constant.

*Note (2026-10-07).* The analogy follows the original hypothesis. For
the constant `weeks_on_drug` predictor actually simulated, $t_{1/2} = 0$
is the case of largest mismatch (Section 3.4), so Scenario A does not
describe the simulated analysis.

---

## 5. Simulation Evidence

### 5.1 Fine-Grained Carryover Spectrum Simulation

We conducted a simulation study across 11 carryover half-life values
(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0 weeks) with 100
iterations each (simulated; $N = 70$, seed 42). The results show a
large effect of carryover modeling. All values below were re-read from
the archived `carryover_spectrum_results.RData` on 2026-10-07 and agree
with the table (verified); no fit failed.

| $t_{1/2}$ (weeks) | Power (Without) | Power (With) | Bias % (Without) | Bias % (With) |
|------------|-----------------|--------------|------------------|---------------|
| 0.00       | 6%              | 100%         | +113%            | +28%          |
| 0.05       | 1%              | 100%         | +95%             | -2%           |
| 0.10       | 1%              | 100%         | +96%             | -1%           |
| 0.15       | 2%              | 100%         | +95%             | -1%           |
| 0.20       | 0%              | 100%         | +96%             | -1%           |
| 0.30       | 2%              | 100%         | +92%             | +4%           |
| 0.50       | 2%              | 100%         | +88%             | +3%           |
| 0.75       | 8%              | 100%         | +84%             | +12%          |
| 1.00       | 12%             | 100%         | +80%             | +16%          |
| 1.50       | 17%             | 100%         | +73%             | +24%          |
| 2.00       | 36%             | 100%         | +64%             | +32%          |

Table: Power and bias by carryover half-life with and without carryover modeling (source: simulation_carryover_spectrum.R, N = 70, 100 iterations; the with-carryover analysis assumes a half-life of 0.2 weeks; positive bias denotes attenuation toward zero)

Power is a proportion of 100 replicates, so its MCSE is at most 0.05
(0.048 at 0.36, 0.024 at 0.06, 0 at 100%). The bias is relative to the
on-drug target $\beta_{bm:D} = -0.15$; a bias above +100%, as at
$t_{1/2} = 0$ without carryover modeling (mean estimate $+0.020$),
means the estimate has the wrong sign.

### 5.2 Key Findings from Simulation

**Without Carryover Modeling:**

1. **Severe bias at short carryover**: At $t_{1/2}$ = 0.05-0.3 weeks, estimates are
   biased toward zero by 92-96%, making the biomarker moderation effect
   nearly undetectable.

2. **Power increases with longer carryover**: Counter to naive expectation,
   power INCREASES from 0-2% at short carryover to 36% at $t_{1/2}$ = 2 weeks.

3. **Mechanism**: Longer carryover means the off-drug periods retain more
   of the treatment effect. Since the analysis predictor (`weeks_on_drug`)
   remains constant after discontinuation, it becomes a better proxy for
   the true cumulative effect when carryover persists.

**With Carryover Modeling:**

1. **Uniform high power**: The carryover-adjusted analysis rejects in
   100 of 100 replicates in every carryover scenario.

2. **Moderate bias**: Bias ranges from -2% to +32%. It is smallest
   where the true half-life is near the assumed 0.2 weeks (-2% to -1%
   at 0.05-0.2) and grows as the true half-life departs from it, in
   either direction (+28% at 0, +32% at 2 weeks). It is therefore a
   half-life misspecification bias of the adjusted analysis (corrected
   2026-10-07 from "a separate issue related to model
   parameterization").

### 5.3 Predictor-Outcome Correlation Analysis

Analysis from `diagnose_nonmonotonic_power.R` shows how the correlation between
the simple predictor (`weeks_on_drug × biomarker`) and the true biological
response varies with carryover:

| Carryover $t_{1/2}$ | Cor(predictor, true_BR) | Squared correlation |
|--------------|-------------------------|-----|
| 0.0 weeks    | 0.214                   | 0.046 |
| 0.5 weeks    | 0.342                   | 0.117 |
| 1.0 weeks    | 0.421                   | 0.177 |
| 2.0 weeks    | 0.525                   | 0.276 |

Table: Correlation between the simple predictor (weeks on drug times standardized biomarker) and true BR, and its square, by half-life (one simulated data set per half-life)

*Correction (2026-10-07, verified by rerunning
`diagnose_nonmonotonic_power.R`).* The correlations reproduce exactly.
The third column was labeled $R^2$ but is the square of the second
column; the script's own $R^2$ of `weeks_on_drug` alone against true
$BR$ is 0.084, 0.013, 0.004 and 0.000, which falls with half-life.
Each row is a single data set of 70 participants, the diagnostic uses
$\beta_{bm} = 0.5$ rather than 0.3 and measurement weeks 4, 8, 9-12, 16
and 20, so it is a descriptive illustration, not a replicated result.

Counterintuitively, the correlation INCREASES with longer carryover. This is
because:

1. At $t_{1/2}$ = 0, off-drug observations have zero true effect, contributing zero
   covariance but positive predictor variance
2. At $t_{1/2}$ > 0, off-drug observations have positive true effect that correlates
   with the (constant) predictor value

### 5.4 Off-Drug Retention Profiles

The carryover retention at off-drug timepoints (the last dose is in
week 8; $t_{sd}$ is counted from the first off-drug visit, week 9):

| Week | $t_{sd}$ | $t_{1/2}$ = 0.1 | $t_{1/2}$ = 0.2 | $t_{1/2}$ = 0.5 | $t_{1/2}$ = 1.0 | $t_{1/2}$ = 2.0 |
|------|-----|--------|--------|--------|--------|--------|
| 9    | 0   | 1.000  | 1.000  | 1.000  | 1.000  | 1.000  |
| 10   | 1   | 0.001  | 0.031  | 0.250  | 0.500  | 0.707  |
| 12   | 3   | 0.000  | 0.000  | 0.016  | 0.125  | 0.354  |
| 16   | 7   | 0.000  | 0.000  | 0.000  | 0.008  | 0.088  |
| 20   | 11  | 0.000  | 0.000  | 0.000  | 0.000  | 0.022  |

Table: Carryover retention at off-drug weeks by carryover half-life (values of $2^{-t_{sd}/t_{1/2}}$, checked 2026-10-07)

At $t_{1/2}$ = 0.1-0.2 weeks, the effect has essentially vanished by week 10, creating
a mismatch between the constant predictor and the rapidly decaying true effect.

At $t_{1/2}$ = 2 weeks, substantial effect (35-70%) persists through week 12, making
the off-drug period appear more similar to the on-drug period, reducing the
model misspecification.

### 5.5 Reconciling with Hendrickson et al. Findings

Hendrickson et al. (2020) reported that "even a short (0.1 weeks) carryover
component resulted in a precipitous decline in power." Our simulation confirms
this finding:

- Power drops from ~6% at $t_{1/2}$ = 0 to 0-2% at $t_{1/2}$ = 0.1-0.3 weeks
- This represents a collapse in detection ability

*Correction (2026-10-07).* This simulation does not confirm the
published finding. The drop from 6% to 0-2% is within about two MCSEs,
starts from a baseline that is already near the nominal level, and at
$t_{1/2} = 0$ the mean estimate has the wrong sign; there is no
collapse from high power to reproduce. The published collapse (Hybrid
power at $N = 70$, $c_{bm} = 0.3$ falling from 0.74 to 0.12 at
$t_{1/2} = 0.1$) arises in a covariance-moderation model and has since
been explained by the published step rule for the biomarker-response
coupling: with any positive half-life the coupling is constant across
visits, the coupling gap between on-drug and off-drug visits is zero,
and the interaction vanishes; graded coupling ($2^{-t_{sd}/t_{1/2}}$
off drug, 0 before first exposure) removes the collapse (`docs/36`,
Sections 5.3 to 5.5 and 6.8; `docs/37`). A mean-moderation model with
linear growth, as here, cannot exhibit that mechanism.

The apparent non-monotonicity observed in the original heatmaps (power dropping
then recovering) depends on the specific design parameters and the baseline
power at $t_{1/2}$ = 0. When baseline power is already low (as in our simplified
simulation), the pattern appears as monotonic recovery with longer carryover.
When baseline power is high (as in the full Hendrickson simulation with more
complex data generation), the drop-then-recovery pattern is more pronounced.
(This explanation was not tested: no high-baseline run of the
simplified model was made, and the published Hendrickson grid does not
extend beyond 0.2 weeks.)

---

## 6. Implications and Recommendations

### 6.1 Trial Design Implications

1. **Very short carryover ($t_{1/2}$ < 0.2 weeks)** represents the "worst case" for
   analysis without carryover modeling, not zero carryover.

2. **Zero carryover** is actually the **best case** for simple analysis models,
   providing clean on/off separation.

*Note (2026-10-07).* Items 1 and 2 restate the original hypothesis.
In this paper's simulation power without carryover modeling is 0.00-0.06
over 0 to 0.5 weeks with no reliable difference between zero and very
short carryover (Section 5.1), and for the constant `weeks_on_drug`
predictor zero carryover is the case of largest mismatch (Section 3.4).

3. **Long carryover ($t_{1/2}$ > 1 week)** partially recovers power even without
   explicit carryover modeling, because the effect remains relatively stable.

### 6.2 Analysis Recommendations

1. **When carryover is suspected but unknown:**
   - Always include carryover modeling in the analysis
   - The exponential decay model is robust across different true decay forms
     (exponential, linear, Weibull) per our robustness simulations (not
     reported or cited in this paper; unverified here)

2. **Power calculations:**
   - Do not assume power degrades monotonically with carryover
   - The minimum power scenario is at intermediate carryover ($t_{1/2} \approx$ 0.1-0.5
     weeks for typical trial durations; in this paper's simulation the
     minimum is not distinguishable from the $t_{1/2} = 0$ value)

3. **Washout periods:**
   - If carryover modeling is not feasible, longer washout periods should
     target at least 5× the expected half-life
   - For behavioral or unknown carryover, analytical approaches (modeling or
     down-weighting early post-discontinuation observations) are preferable

### 6.3 Future Research Directions

1. Optimal down-weighting schemes for early post-discontinuation observations
2. Robust estimation methods that are insensitive to carryover misspecification
3. Adaptive designs that can estimate carryover parameters during the trial

---

## 7. Conclusions

The relationship between carryover half-life and statistical power in simple
analysis models reflects a fundamental bias-variance tradeoff arising from
model misspecification. Our key findings are:

### 7.1 Primary Conclusions

1. **Carryover modeling is essential**: Without explicit carryover modeling,
   power ranges from 0-36% across the carryover spectrum. With carryover
   modeling, power is uniformly 100% in this simplified simulation. The
   difference shows the importance of appropriate model specification.

2. **Short carryover is the worst case**: At $t_{1/2}$ = 0.1-0.3 weeks, the simple
   model exhibits 92-96% bias toward zero, making effect detection nearly
   impossible (0-2% power).

3. **Longer carryover paradoxically helps simple models**: As carryover
   increases beyond 0.5 weeks, power gradually recovers (up to 36% at
   $t_{1/2}$ = 2 weeks) because the persistent effect better matches the constant
   predictor assumption.

4. **The "precipitous decline" is real**: Hendrickson et al.'s observation
   that even 0.1 weeks of carryover causes dramatic power loss is confirmed.
   This represents the critical transition from a well-specified to a
   misspecified model.

*Correction (2026-10-07).* Conclusion 4 is withdrawn. The simulation
starts from a power of 0.06, not from a well-specified high-power
model, so it cannot confirm a collapse; and the published decline has
been traced to the step coupling rule of the covariance-moderation
model (`docs/36`, Section 5; `docs/37`), not to the misspecification
studied here. Conclusion 2 should also be read with Section 3.4: the
worst case for the constant predictor is $t_{1/2} = 0$ in terms of
bias, and power at 0 and at 0.1-0.3 weeks is not reliably different.

### 7.2 Practical Recommendations

1. **Always model carryover** when any carryover is suspected, regardless
   of the expected half-life.

2. **Exponential decay is robust**: Our robustness analysis shows the
   exponential carryover model performs well even when the true decay
   follows linear or Weibull patterns.

3. **Power calculations must account for model choice**: Power under
   carryover misspecification cannot be extrapolated from no-carryover
   scenarios.

### 7.3 Implications for Trial Design

These findings have direct implications for precision medicine trial design:

- **OL+BDC designs** require carryover modeling to achieve adequate power
- **N-of-1 designs** with multiple crossovers are particularly vulnerable
  to carryover effects (not simulated in this paper)
- **Washout periods** of 5× the expected half-life may be insufficient
  for behavioral carryover, which a washout does not remove (Section
  2.2); analytical approaches are preferred

Understanding this tradeoff enables more informed trial design decisions and
more accurate power calculations for precision medicine validation studies.

---

## Appendix: Generated Figures

The following figures were generated by `simulation_carryover_spectrum.R` and
are available in the archived project at
`archive/pmsimstats2025-2026-04-23/pmsimstats2025/analysis/output/`
(not reproduced here):

1. `carryover_spectrum_power.png`: power curves across the carryover spectrum
2. `carryover_spectrum_bias.png`: bias curves across the carryover spectrum
3. `carryover_spectrum_tradeoff.png`: combined bias-power visualization

---

## References

*Note (2026-10-07).* All entries were checked in PubMed. Three were
wrongly attributed in the original and are corrected: the Hendrickson
et al. (2020) entry gave a nonexistent *Journal of Neurotrauma* article,
the crossover-design review is by Lim and In (2021), not Lee and Ahn
(2022), and the behavioral carryover paper is by Shi and Ye (2024) in
*Biometrics*, not Zhou et al. in *Biostatistics*.

Freeman, P. R. (1989). The performance of the two-stage analysis of two-treatment,
two-period crossover trials. *Statistics in Medicine*, 8(12), 1421-1432.
doi:10.1002/sim.4780081202. (PMID 2616932.)

Grizzle, J. E. (1965). The two-period change-over design and its use in clinical
trials. *Biometrics*, 21(2), 467-480. (PMID 14338679.)

Hendrickson, R. C., Thomas, R. G., Schork, N. J., & Raskind, M. A.
(2020). Optimizing aggregated N-of-1 trial designs for predictive
biomarker validation: Statistical methods and theoretical findings.
*Frontiers in Digital Health*, 2, 13. doi:10.3389/fdgth.2020.00013.
(PMID 34713026.)

Lim, C.-Y., & In, J. (2021). Considerations for crossover design in
clinical study. *Korean Journal of Anesthesiology*, 74(4), 293-299.
doi:10.4097/kja.21165. (PMID 34344139.)

Shi, D., & Ye, T. (2024). Behavioral carry-over effect and power
consideration in crossover trials. *Biometrics*, 80(2), ujae023.
doi:10.1093/biomtc/ujae023. (PMID 38563531.)

Wellek, S., & Blettner, M. (2012). On the proper use of the crossover design in
clinical trials: Part 18 of a series on evaluation of scientific publications.
*Deutsches Ärzteblatt International*, 109(15), 276-281.
doi:10.3238/arztebl.2012.0276. (PMID 22567063.)
