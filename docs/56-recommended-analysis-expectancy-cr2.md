---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# The Recommended Analysis of the Biomarker-by-Drug Interaction in Hybrid N-of-1 Trials: Expectancy Terms, a Working Covariance and CR2 Inference, with Worked Examples {.unlisted .unnumbered}
*2026-10-09 09:41 PDT*

**Author.** pmsimstats team

**Purpose.** The exploratory stress test of `docs/46` and `docs/47`
compared fourteen analyses of the biomarker-by-drug interaction
$\beta_{bm:D}$ in the Hybrid design and selected one of them, called
F1. This note explains that analysis term by term: what each part of
the model is for, what goes wrong without it, and what it costs. Two
worked examples carry the explanation. The first fits F1 and the
simpler analysis S to one simulated trial, with every coefficient,
standard error and degree of freedom shown. The second reads the
stored stress-test results to show the same contrast over thousands
of trials. A short section explains why choosing the carryover
half-life by AIC, as the variant F1A did, fails the size criterion,
and what the project intends instead; the full treatment of a
data-chosen half-life is in the companion note `docs/57`. The note
ends with a recipe in R. It is written for trial statisticians and
methodologists who will use the analysis, including the authors of
the original pmsimstats code [1, 2], and can be read without the other
project documents.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\listoffigures
\clearpage
```

## Notation index and glossary

Notation follows the project's canonical table
(`analysis/report/NOTATION.md`). Symbols marked *local* are defined
here for this note only. (Notation, 2026-10-09: the rate in the MCSE
formula is now $\hat\pi$, since $q$ is the canonical path index; the
normal law is $\mathcal N$; the time term is $\beta_T w_t$, with $w_t$
the calendar week; symbols previously marked local or proposed
that `NOTATION.md` now lists are marked canonical. "The published
construct" now means the compound-symmetry covariance of the published
process only, so the AR(1) variants of the `04-power-simulation.R`
construct quoted from the compendium are named as such.)

| Symbol | Meaning | Status |
|------------|--------------------------------------------------|----------|
| $Y_{it}$ | symptom score of participant $i$ at visit $t$ (code `Sx`); $t = 0$ is baseline; higher is worse | canonical |
| $\mathrm{BL}_i$, $TV_{it}$, $PB_{it}$, $BR_{it}$ | baseline level; natural-course, placebo (expectancy) and drug responses; $Y_{it} = \mathrm{BL}_i - (TV_{it} + PB_{it} + BR_{it})$ | canonical |
| $B_i$ | biomarker of participant $i$ (code `bm`) | canonical |
| $B_i^c$ | biomarker centered at the trial's baseline mean, $B_i - \bar B$, not scaled (code `bmc`) | canonical |
| $D_{it}$ | binary drug state at visit $t$ | canonical |
| $D_{bc,it}$ | exposure-decayed drug indicator: 1 on drug, $2^{-t_{sd}/t_{1/2}^a}$ off drug after first exposure, 0 at baseline and before first exposure (code `Dbc`; the stress-test driver calls it `X`) | canonical |
| $t_{sd}$ | time since discontinuation, in weeks (code `tsd`); 0 on drug and before first exposure | canonical |
| $t_{1/2}$ | true carryover half-life of the data-generating process | canonical |
| $t_{1/2}^a$ | half-life assumed by the analysis in $D_{bc,it}$, in weeks (code `t_half_a`) | canonical |
| $\eta_t$ | expectancy: 0 at baseline, 1 in the open-label phase, 0.5 when blinded (code `De`) | canonical |
| $\beta_0$, $\beta_T$ | intercept and coefficient of linear time $w_t$ (code `t`) | canonical |
| $w_t$ | calendar week of visit $t$ (code column `t`) | canonical |
| $\beta_B$ | biomarker main effect (coefficient of `bmc`) | canonical |
| $\beta_D$ | drug main effect (coefficient of `Dbc`) | canonical |
| $\beta_{bm:D}$ | biomarker-by-drug interaction, the estimand (coefficient of `bmc:Dbc`) | canonical |
| $\beta_\eta$, $\beta_{bm:\eta}$ | coefficients of `De` and `bmc:De` | canonical |
| $u_i$ | participant random intercept, variance $\sigma_u^2$ | canonical |
| $\varepsilon_{it}$ | observation-level residual, variance $\sigma_\varepsilon^2$; $\hat\varepsilon_i$ the residual vector of participant $i$ | canonical |
| $\rho_a$ | working CAR(1) correlation per week | canonical |
| $c_{bm}$ | biomarker-drug-response coupling (the predictive effect) | canonical |
| $c_{bm,TV}$ | biomarker-natural-course coupling (a prognostic effect) | canonical |
| $\sigma_{BR}$, $\sigma_{bm}$, $\sigma_{TV}$ | standard deviations of the drug response, the biomarker and the natural-course component | canonical |
| $s_{\text{prog}}$ | prognostic slope of $Y$ on $B$ that does not depend on drug | canonical |
| $f$ | share of $s_{\text{prog}}$ absorbed by a shared $\beta_B$ | canonical |
| $n_0, n_1, n_2, n_3$ | total weights (row counts) of baseline, open-label on-drug, blinded on-drug and blinded off-drug rows | local |
| $X_i$, $\hat V_i$, $H_{ii}$ | participant $i$'s design matrix, fitted working covariance and hat-matrix block | local |
| $\ell_g$ | terms of the Satterthwaite degrees-of-freedom formula for CR2 (compendium, eq:cr2-df) | canonical |
| $\kappa$ | mean standard error over the empirical SD of the estimate | canonical |
| $\nu$ | denominator degrees of freedom | canonical |
| $N$ | total sample size across paths | canonical |
| $\alpha$ | nominal significance level, 0.05 | canonical |
| $n_{\text{sim}}$ | Monte Carlo replicates (trials) in a cell | canonical |
| $\hat\pi$ | an estimated rejection rate, in the MCSE formula (formerly $q$, which is the canonical path index) | canonical |

Table: Notation index for this note

**Glossary.**

- **Hybrid design.** Eight visits at weeks 4, 8, 9, 10, 11, 12, 16 and
  20 after a baseline at week 0: open-label treatment at weeks 4 and 8,
  a blinded discontinuation, and a brief crossover, over four
  randomization paths [1]. Every open-label visit is on drug.
- **Predictive biomarker.** Changes the size of the drug effect; this is
  what $\beta_{bm:D}$ measures.
- **Prognostic biomarker.** Associated with the symptom score whatever
  the treatment, for example through the natural course or the placebo
  response.
- **Working covariance.** The covariance structure an analysis assumes
  for a participant's repeated scores in order to weight them. It need
  not be the true covariance.
- **CR2.** The bias-reduced cluster-robust (sandwich) variance estimator
  of Bell and McCaffrey [3], clustered on participant, with the
  Satterthwaite-type degrees of freedom of Pustejovsky and Tipton [4].
- **Size.** The false-positive rate of the two-sided test of
  $\beta_{bm:D} = 0$ at $\alpha = 0.05$.
- **Size-adjusted power.** Power when each analysis rejects at the 5%
  quantile of its own $p$-values in the matching null cell, so that
  analyses are compared at equal size.
- **MCSE.** Monte Carlo standard error, $\sqrt{\hat\pi(1 - \hat\pi)/n_{\text{sim}}}$ for a
  rate $\hat\pi$ over $n_{\text{sim}}$ trials [10]: 0.0069 at $\hat\pi = 0.05$ with
  $n_{\text{sim}} = 1{,}000$, at most 0.022 with $n_{\text{sim}} = 500$.
- **The analyses.** All are fitted to the long data with baseline as a
  response row at $t = 0$, linear time in weeks, and $D_{bc,it}$ at
  $t_{1/2}^a = 0.5$ week unless stated.

| ID | Mean model | Working covariance | Inference |
|----|------------------------------|------------------|----------------|
| S | `Sx ~ bmc + t + Dbc + bmc:Dbc` | $u_i$ + CAR(1) | CR2 |
| F1 | S + `De + bmc:De` | $u_i$ + CAR(1) | CR2 |
| F1A | F1, $t_{1/2}^a$ chosen by AIC over {0, 0.25, 0.5, 1, 2} | $u_i$ + CAR(1) | CR2 |
| ES | S + `bmc:bl0` (`bl0` = 1 at baseline) | $u_i$ + CAR(1) | CR2 |
| G1, G2 | F1 | $u_i$ + CAR(1), residual SD by phase | CR2 (G1), model-based (G2) |
| DS | S, $t_{1/2}^a$ chosen by AIC | $u_i$ + CAR(1) | CR2 |
| R | the published-type reference: S with binary $D_{it}$ | compound symmetry | model-based, Satterthwaite |

Table: The analyses discussed in this note (`docs/47`, notation table; inspected in the driver)

**Evidence labels.** Each claim carries one of six labels: *quoted*
(reproduced from a source), *inspected* (read from source code or a
stored file), *derived* (follows algebraically), *verified* (computed
by running code; the script is named in Section 10), *argued*
(reasoned but not computed), or *proposal* (a procedure not yet run).

**Status of the stored results.** Every Monte Carlo number in this note
is from the stress test as run on 2026-10-07. That run used the
published recursive carryover on the drug-response mean and a fixed
allocation (18, 18, 17, 17 at $N = 70$), in configuration C at
$\rho = 0.7$ (`docs/47`, Section 2). The driver now defaults to
anchored carryover and randomized allocation; a rerun under those
defaults is pending. The worked examples of Sections 4 and 6 use the
same recursive, fixed-allocation construct so that they match the
stored results (inspected).

## 1. Summary for the busy reader

**The model.** For participant $i$ at visit $t$ (baseline included),

$$
\begin{aligned}
Y_{it} = {} & \beta_0 + \beta_B B_i^c + \beta_T w_t
  + \beta_D D_{bc,it} + \beta_{bm:D}\, B_i^c D_{bc,it}
  + \beta_\eta \eta_t + \beta_{bm:\eta}\, B_i^c \eta_t
  + u_i + \varepsilon_{it}, \\
& u_i \sim \mathcal N(0, \sigma_u^2), \qquad
  \mathrm{Cor}(\varepsilon_{it}, \varepsilon_{is}) = \rho_a^{|t - s|} \ \text{(weeks)},
\end{aligned}
$$

fitted by REML with `nlme::lme`, and the test of $\beta_{bm:D} = 0$
made with the CR2 cluster-robust standard error and its
Satterthwaite degrees of freedom (`clubSandwich::coef_test`). Here
$D_{bc,it}$ uses a half-life $t_{1/2}^a$ fixed in advance from
pharmacology (0.5 week in the stress test) and is 0 before any
exposure.

**When to use it.** For a Hybrid or similar N-of-1 design whose
open-label phase is all on drug, whenever the biomarker might also be
prognostic, which for a physiological marker such as blood pressure
cannot usually be excluded (`docs/49`, Section 6).

**What it protects against** (verified, stored stress test).

- *A prognostic biomarker.* With the biomarker correlated 0.3 with the
  natural course or with the placebo response and no interaction, S
  rejects 0.099 and 0.100 and the published-type analysis R 0.107 and
  0.159; F1 rejects 0.047 and 0.040.
- *A misspecified working covariance.* Over 11 standard null cells F1
  rejects 0.051 pooled, and at most 0.062 in any cell, including the
  cells where the true half-life (0 or 1 week) differs from the assumed
  0.5 week. R, with a compound-symmetry working covariance and
  model-based inference, rejects 0.069 pooled and 0.081 at $N = 35$.

**What it costs.** Power. At $c_{bm} = 0.30$, $N = 70$, under AR(1)
data, F1's size-adjusted power is 0.364 against 0.514 for S, a
relative loss of 29%; over ten AR(1) alternative cells the loss ranges
from 9% to 30%. Under compound-symmetry data it is small (0.522
against 0.554) (verified). Plan the sample size for F1, not for S.

**The one thing not to do.** Do not choose the half-life by AIC and then
test as if it had been fixed. That is F1A, and it fails the size
criterion: 0.057 pooled and 0.073 when the true half-life is 1 week
(verified; Section 6).

**Status.** The validated analysis is F1 with a half-life fixed in
advance. The project's intended analysis lets the data choose the
half-life, by a procedure that does not reuse the interaction's noise
(the proposal F1A0) or that calibrates the test for the selection
(the proposal F1AB). Neither has been evaluated yet; both are set out
in `docs/57`. Until they are, report F1 at the fixed half-life, with a
sensitivity analysis over a grid of half-lives and the data's own
estimate of the half-life reported alongside (Section 8).

## 2. Scope

This note covers one analysis family for one estimand in one design
family: the biomarker-by-drug interaction in aggregated Hybrid N-of-1
trials with a Gaussian outcome. It does not reopen the selection; the
rule, the cells and the amendments are in `docs/46` and `docs/47`. It
explains the selected analysis and gives worked examples. Where an
argument is approximate, it says so.

## 3. The model, term by term

### 3.1 The centered biomarker and the time adjustment

**Centering.** The biomarker enters as $B_i^c = B_i - \bar B$, centered
at the trial's baseline mean and not scaled (inspected, driver
`make_trial()`). Centering makes $\beta_D$ the drug effect at an
average biomarker and keeps the biomarker block separate from the rest
of the model: because $B_i^c$ has mean zero and is independent of the
randomized path, its cross-products with the intercept, time and
$D_{bc,it}$ vanish in the limit, and the interaction's probability limit
depends only on the biomarker columns (derived;
`analysis/report/compendium-mathematics.tex`, eq:bp-plim). Scaling
changes nothing that is tested: $t$-statistics and $p$-values are
invariant to rescaling $B_i^c$, and $\hat\beta_{bm:D}$ is then per unit
of the raw biomarker (inspected; `docs/47`, note below the notation
table). In this note's construct $\sigma_{bm} = 15.36$, so a slope of
$-0.195$ per unit is $-3.0$ points per SD.

**Time.** Linear time in weeks absorbs the average drift of the natural
course and the placebo response. It enters only the non-biomarker
block, so it does not change the interaction's probability limit, and
`docs/42` found estimates that agree to within 0.0009 across time
adjustments (quoted there; derived in the compendium, eq:bp-plim). A
richer time adjustment is a matter of efficiency, not of bias, for
$\beta_{bm:D}$.

### 3.2 $D_{bc,it}$ and the pre-exposure guard

**What it does.** The drug effect does not vanish on the day the drug is
stopped. $D_{bc,it}$ codes exposure as 1 on drug and
$2^{-t_{sd}/t_{1/2}^a}$ afterwards, so the first blinded off-drug visit
still counts as partly exposed; carryover in within-participant
designs, and the cost of ignoring it, are long-standing concerns [9].
When $t_{1/2}^a$ equals the true
half-life and the biomarker coupling decays the same way, the true
biomarker slopes lie exactly in the span of $(1, D_{bc,it})$ and the
interaction is recovered in full under any working covariance (derived;
compendium, eq:invf-dbc-exact). With $t_{1/2}^a = 0$ it reduces to the
binary $D_{it}$.

Table 3 shows the coding on the four Hybrid paths at
$t_{1/2}^a = 0.5$ week (verified).

| Week | $\eta_t$ | Path 1 | Path 2 | Path 3 | Path 4 |
|------|------|------|------|------|------|
| 0 (baseline) | 0 | 0 | 0 | 0 | 0 |
| 4 | 1 | 1 | 1 | 1 | 1 |
| 8 | 1 | 1 | 1 | 1 | 1 |
| 9 | 0.5 | 1 | 1 | 1 | 1 |
| 10 | 0.5 | 1 | 1 | 0.25 | 0.25 |
| 11 | 0.5 | 0.25 | 0.25 | 0.0625 | 0.0625 |
| 12 | 0.5 | 0.0625 | 0.0625 | 0.0156 | 0.0156 |
| 16 | 0.5 | 1 | 0.0002 | 1 | 0.0001 |
| 20 | 0.5 | 0.0039 | 1 | 0.0039 | 1 |

Table: $D_{bc,it}$ at $t_{1/2}^a = 0.5$ week on the four Hybrid paths (verified)

**The guard.** The time since discontinuation is 0 both on drug and
before the first exposure, and $2^{0} = 1$. A decay function written as
`ifelse(on_drug, 1, 0.5^(tsd / t_half))` therefore codes a participant
who starts on placebo as fully exposed at every placebo visit. In the
Hybrid design every path starts on drug, so the guard never acts there
(inspected, Table 3); it matters in any design with a placebo-first
sequence. On the placebo-first path of the two-period crossover (eight
visits 2.5 weeks apart, drug from week 12.5) at $t_{1/2}^a = 1$ week,
the guarded indicator is 0, 0, 0, 0, 1, 1, 1, 1 and the unguarded one
is 1 at all eight visits (verified). The unguarded participant then
has no within-participant variation in exposure and contributes nothing
to $\beta_{bm:D}$. The project found this in its own carryover
manuscript: in a paired rerun of its crossover cells ($N = 70$,
$c_{bm} = 0.45$, $t_{1/2} = 1$ week, 500 replicates), restoring the
guard raised the power of the exposure-weighted analysis from 0.512 to
0.814, the level of the binary indicator, 0.810 (quoted; report 01,
Appendix B.4). The guarded form, as in the driver and in the patched
`exposuredecay()` (inspected):

```{.r}
exposure_decay <- function(on_drug, tsd, t_half) {
  if (t_half <= 0) return(as.numeric(on_drug))
  ifelse(on_drug, 1, ifelse(tsd > 0, 0.5^(tsd / t_half), 0))
}
```

### 3.3 The expectancy terms $\eta_t$ and $B_i^c\eta_t$

This is the part of F1 that the other analyses lack, and the reason it
was selected.

**The confounding the Hybrid design builds in.** Every open-label visit
is on drug and carries full expectancy, $\eta_t = 1$; every blinded
visit has $\eta_t = 0.5$; baseline has neither. Drug and expectancy
therefore move together. Over the nine rows at the stress-test
allocation, the allocation-weighted correlation of $\eta_t$ with the
binary $D_{it}$ is 0.588 (0.611 with $D_{bc,it}$; 0.507 over the eight
post-baseline rows) (verified; it reproduces the compendium's 0.588).
Anything that scales with expectancy looks, on average, like something
that scales with drug.

**What a prognostic biomarker does to a shared slope.** Suppose the
biomarker predicts the natural course but not the drug effect, so the
true slope of $Y$ on $B$ is 0 at baseline and $s_{\text{prog}}$ at
every later visit, on drug or off. S has two biomarker slopes:
$\beta_B$ for every row coded off drug, baseline included, and
$\beta_B + \beta_{bm:D}$ on drug. It must fit one $\beta_B$ to the
baseline rows (true slope 0) and the off-drug rows (true slope
$s_{\text{prog}}$), so $\hat\beta_B$ settles between them,
$\hat\beta_B \approx f s_{\text{prog}}$. The on-drug rows have slope
$s_{\text{prog}}$, and the only parameter left to close the gap is the
interaction:

$$
\hat\beta_{bm:D} \approx (1 - f)\, s_{\text{prog}} \neq 0
$$

(derived; `docs/48`, Section 2; compendium, eq:bp-leak). Under working
independence with binary $D_{it}$, $1 - f = n_0/(n_0 + n_3)$, the weight
share of the baseline rows among the rows coded off drug. In the
Hybrid design at the stress-test allocation the row counts are
$n_0 = 70$, $n_1 = 140$, $n_2 = 176$ and $n_3 = 244$ (verified), so
$1 - f = 0.223$. A correlated working covariance changes the factor:
the compendium computes 0.276 for S's random-intercept CAR(1) fit in
the `04-power-simulation.R` construct with AR(1) components at
$\rho = 0.7$, and 0.629 in configuration C under the true
covariance (quoted; compendium, tab:bp-leak and tab:invd-leakC). The
bias is fixed, not random, so the false-positive rate grows with $N$.
That a baseline response row imposes such a constancy constraint is
known from the literature on baseline as a covariate or a response
[8].

**What the expectancy terms do.** F1 adds $\eta_t$ and $B_i^c \eta_t$.
The biomarker slope the model can now fit is different in each of the
four groups of rows:

| Rows | S | F1 |
|---|---|---|
| baseline ($\eta_t = 0$) | $\beta_B$ | $\beta_B$ |
| open label, on drug ($\eta_t = 1$) | $\beta_B + \beta_{bm:D}$ | $\beta_B + \beta_{bm:\eta} + \beta_{bm:D}$ |
| blinded, on drug ($\eta_t = 0.5$) | $\beta_B + \beta_{bm:D}$ | $\beta_B + 0.5\beta_{bm:\eta} + \beta_{bm:D}$ |
| blinded, off drug ($\eta_t = 0.5$) | $\beta_B$ | $\beta_B + 0.5\beta_{bm:\eta}$ |

Table: Biomarker slope fitted by S and F1 in each group of rows, with $D_{bc,it}$ near 1 on drug and near 0 off drug (derived; `docs/49`, Section 3)

Two things follow (derived; compendium, Section "Freeing the baseline
slope, and the expectancy terms").

1. *Baseline is freed.* Because $\eta_0 = 0$, the baseline slope is no
   longer tied to the blinded off-drug slope; $\beta_{bm:\eta}$ absorbs
   the difference.
2. *An expectancy-scaled prognostic slope is absorbed exactly.* If the
   biomarker acts through the placebo response, its slope is
   proportional to $\eta_t$, which lies in the span of F1's biomarker
   columns, so the leak is zero under any working covariance.

The interaction is then identified mainly from the comparison the
Hybrid design was built to make clean: blinded on drug against blinded
off drug, at the same expectancy.

**Candidly: F1 is approximately, not exactly, unbiased here.** For a
prognostic slope that is constant after baseline (the natural-course
pathway), F1's fitted slopes obey one linear restriction, and the
weighted projection under working independence with binary $D_{it}$ gives

$$
\hat\beta_{bm:D} \to s_{\text{prog}}\,
  \frac{1/n_3 - 1/n_2}{1/n_0 + 1/n_1 + 1/n_2 + 1/n_3}
$$

(derived; compendium, eq:bp-f1-leak). The leak vanishes only when the
blinded on-drug and off-drug rows carry equal weight, $n_2 = n_3$. In
the Hybrid design $n_2 = 176$ and $n_3 = 244$, so the factor is
$-0.051$ under independence with binary $D_{it}$, $-0.020$ with $D_{bc,it}$,
0.089 under a random-intercept CAR(1) working covariance and 0.171
under the true AR(1) covariance ($\rho = 0.8$) of the
`04-power-simulation.R` construct, and 0.162
under the true covariance of configuration C (quoted; compendium,
tab:bp-leak and tab:invd-leakC). These factors are a fifth to a third
of S's under the same covariances, which is why F1's size stayed at
0.047 in that cell; but the protection is quantitative, not exact.

**ES has the mirror-image weakness.** ES frees baseline directly with
`bmc:bl0` and is exact for any pathway constant after baseline, but an
expectancy-scaled slope is larger at the open-label visits, which are
all on drug, and leaks into ES's interaction by 0.14 to 0.24 of
$s_{\text{prog}}$ depending on the working covariance (quoted;
compendium, tab:bp-leak). Its stress-test size in the placebo cell,
0.049, shows the leak was small at a coupling of 0.3 (verified), and
its mean estimate there, $-0.027$, against F1's $-0.001$, shows it was
not zero (verified, `summary-full.csv`).

### 3.4 Random intercept plus CAR(1): a working covariance

The random intercept $u_i$ and the continuous-time AR(1) residual,
`corCAR1(form = ~ t | ptID)` on the week scale, are a *working*
covariance: a parsimonious device, three parameters
($\sigma_u^2$, $\sigma_\varepsilon^2$, $\rho_a$), for weighting a participant's
nine scores sensibly. It is not a claim that the data have this
structure, and in the stress-test construct they do not. The true
outcome covariance has three phase-specific variance levels (341.6 at
baseline, 709.6 open label, 598.6 blinded), no noise of its own at
baseline, and cross-phase correlations the structure cannot reproduce
exactly (quoted; compendium, Section "The implied outcome covariance
and phase-specific variances"). A random intercept absorbs only the
long-lag part of AR(1) correlation (compendium, eq:inva-ri). The
working model need only give reasonable weights; the correctness of
the standard error is delegated to CR2.

Richer working covariances bought no robustness in the stress test, and
some made the prognostic leak worse. Residual variances by phase (C1,
C2) up-weight the baseline row, whose variance is smallest, and S's
leak factor rises toward 0.63; C1 rejected 0.382 in the
natural-course prognostic cell (quoted; compendium, tab:invd-leakC;
verified, `decision-screen-ds-nosym-all14.csv`). Once F1's expectancy
terms free the baseline slope, the same phase variances are harmless
(G1 0.055, G2 0.062 in that cell) but add little power (Section 7).

### 3.5 CR2 standard errors with Satterthwaite degrees of freedom

**The intuition.** A model-based standard error trusts the working
covariance: it computes the variance of $\hat\beta$ as if the data had
exactly that structure. The sandwich estimator [5] does not. It keeps
the working model's weights (the bread) but measures the variability
of each participant's contribution to the estimating equation (the
meat) from the observed residuals, participant by participant,

$$
\widehat{\mathrm{Var}}(\hat\beta) = \Big(\sum_i X_i' \hat V_i^{-1} X_i\Big)^{-1}
\Big(\sum_i X_i' \hat V_i^{-1} \hat\varepsilon_i \hat\varepsilon_i' \hat V_i^{-1} X_i\Big)
\Big(\sum_i X_i' \hat V_i^{-1} X_i\Big)^{-1}
$$

(quoted; compendium, eq:sandwich). Because participants are
independent, the observed spread of their contributions is a consistent
estimate of the true variance whatever the working covariance. A wrong
working covariance makes the estimate less efficient; it does not make
the standard error wrong.

**Why CR2 and not the plain sandwich.** Residuals are smaller than the
errors they estimate, because each participant's own data pulled the
fit toward them; with 35 to 70 participants the plain sandwich is
noticeably too small. CR2 inflates each participant's residuals by
$(I - H_{ii})^{-1/2}$, with $H_{ii}$ that participant's block of the hat
matrix, which makes the estimator exactly unbiased when the working
model is right and nearly so when it is not [3] (compendium, eq:cr2).

**Why the degrees of freedom matter.** The precision of a sandwich
standard error is governed by the number of participants and by how
unevenly they inform the coefficient, not by the number of rows. For
$\beta_{bm:D}$, a participant-level biomarker times a within-participant
exposure contrast, a few participants with extreme biomarker values
carry much of the information. The Satterthwaite approximation of
Pustejovsky and Tipton [4] matches the first two moments of the CR2
variance to a scaled $\chi^2$ and returns $\nu = (\sum_g \ell_g)^2 /
\sum_g \ell_g^2$ (compendium, eq:cr2-df). In the worked example of
Section 4 it is 25.8 to 30.5 at $N = 70$, against 555 to 557 for the
`nlme` containment degrees of freedom of the model-based test
(verified); over five verification data sets the mean CR2 degrees of
freedom were 21.2 for S and 20.5 for F1 (quoted; patch README,
"Verification results"). At $\nu \approx 25$ the two-sided critical
value is 2.06 rather than 1.96, so the correction matters at the
sample sizes these trials use.

**The evidence.** Three pieces bear on it.

- *The published-type analysis is liberal.* R (compound symmetry,
  model-based standard errors, Satterthwaite degrees of freedom, binary
  $D_{it}$) rejected 0.069 pooled over the 11 standard null cells, 0.081 at
  $N = 35$, and its $\kappa$ was 0.88 to 0.98; the CR2 analyses S and
  F1 rejected 0.048 and 0.051 (verified, Table 8 and
  `56-kappa.csv`). R differs from S in working covariance and exposure
  coding as well as in inference, so this comparison does not isolate
  the inference.
- *The clean comparison is `docs/45`.* There the same compound-symmetry
  working covariance was fitted with model-based and with CR2
  inference: mean size 0.079 against 0.051 over eight null cells, and
  0.158 against 0.060 in the crossover under AR(1) data, with $\kappa$
  0.73 against 0.95 (quoted; `docs/45`, Table 4.1). Kenward-Roger
  inference [6] did not repair it (0.079), because it corrects for
  estimating a covariance assumed correct, not for a wrong one
  (quoted; `docs/45`, Section 4.2).
- *Where the working covariance is close, model-based inference is
  also fine.* With phase-specific residual variances, which approach
  the stress-test construct's covariance, model-based and CR2 inference
  on the same fits gave pooled sizes of 0.052 and 0.053 (C2, C1) and
  0.047 and 0.049 (G2, G1) (verified, `decision-screen-ds-nosym-all14.csv`).
  CR2 is recommended because a trial's analyst cannot know in advance
  that the working covariance is close.

## 4. Worked example 1: one trial, step by step

### 4.1 The trial

One Hybrid trial of $N = 70$ (18, 18, 17, 17 per path) is drawn from the
stress-test cell `prog-TV-const` with the driver's own `build()`,
`draw()` and `make_trial()`: configuration C at $\rho = 0.7$, $t_{1/2} =
0.5$ week, recursive carryover, and the biomarker correlated
$c_{bm,TV} = 0.3$ with the natural-course component at every visit
(inspected). The biomarker is therefore prognostic. The true
prognostic slope is

$$
s_{\text{prog}} = -\frac{c_{bm,TV}\,\sigma_{TV}}{\sigma_{bm}}
  = -\frac{0.3 \times 10}{15.36} = -0.195
$$

points per unit of biomarker after baseline and 0 at baseline (derived;
the minus sign because $Y = \mathrm{BL} - (TV + \cdots)$). The first
version has no interaction, $c_{bm} = 0$. The second uses the same seed,
20261009, with $c_{bm} = 0.30$, for which the true interaction is
$\beta_{bm:D} = -c_{bm}\sigma_{BR}/\sigma_{bm} = -0.30 \times 8/15.36 =
-0.156$ (derived). The two draws share the seed but not the covariance,
so they are not the same participants shifted by an effect; they are
two trials from the same random stream.

Loading the driver's functions without running the driver:

```{.r}
source('analysis/scripts/recommended-analysis/56-load-driver.R')
env <- load_driver(carryover = 'recursive', allocation = 'fixed')
cl <- env$cell('prog-TV-const', 0, c_tv = 0.3)    # c_bm = 0
st <- env$cell_setup(cl)
set.seed(20261009L)
d <- env$make_trial(cl, st$bs, env$alloc(cl, st), 1)
dd <- env$with_x(d, 0.5)                           # X = D_bc at 0.5 week
fit_S  <- env$fit_lme(env$F_S, dd)                 # Sx ~ bmc + t + X + bmc:X
fit_F1 <- env$fit_lme(env$F_F1, dd)                # + De + bmc:De
clubSandwich::coef_test(fit_F1, vcov = 'CR2', test = 'Satterthwaite')
```

The fits are the driver's `fit_lme()` with its `lmeControl`, so they
are exactly the stress-test analyses (inspected; the script checks
that the driver's `cr2()` returns the same standard error). Both fits
together took 0.4 to 1.2 seconds (verified).

### 4.2 The coefficients

| Term | S estimate | S SE, model / CR2 | F1 estimate | F1 SE, model / CR2 |
|--------------------------|------|----------|------|----------|
| intercept | 78.42 | 2.98 / 2.77 | 79.11 | 3.06 / 2.59 |
| `bmc` ($\beta_B$) | $-0.159$ | 0.193 / 0.183 | $-0.171$ | 0.200 / 0.169 |
| `t` | $-0.528$ | 0.119 / 0.090 | $-0.514$ | 0.119 / 0.094 |
| `Dbc` ($\beta_D$) | $-7.50$ | 1.22 / 1.37 | $-6.83$ | 1.38 / 1.42 |
| `De` ($\beta_\eta$) | | | $-2.21$ | 2.11 / 1.93 |
| `bmc:De` ($\beta_{bm:\eta}$) | | | 0.030 | 0.147 / 0.116 |
| `bmc:Dbc` ($\beta_{bm:D}$) | $-0.040$ | 0.085 / 0.085 | $-0.050$ | 0.097 / 0.100 |

Table: Worked example 1, no interaction ($c_{bm} = 0$, prognostic natural course): coefficients of S and F1 with model-based and CR2 standard errors (verified)

The drug lowers the score by about 7 points on drug, and the time slope
is about half a point per week of improvement. The biomarker terms are
all small relative to their standard errors.

### 4.3 The test of the interaction

| Case | Analysis | $\hat\beta_{bm:D}$ | SE (model) | $\nu$ (model) | $p$ (model) | SE (CR2) | $\nu$ (CR2) | $p$ (CR2) |
|------------|------|------|------|------|------|------|------|------|
| $c_{bm} = 0$ | S | $-0.040$ | 0.085 | 557 | 0.636 | 0.085 | 30.5 | 0.637 |
| $c_{bm} = 0$ | F1 | $-0.050$ | 0.097 | 555 | 0.607 | 0.100 | 30.5 | 0.623 |
| $c_{bm} = 0.30$ | S | $-0.296$ | 0.097 | 557 | 0.0024 | 0.072 | 25.8 | 0.0004 |
| $c_{bm} = 0.30$ | F1 | $-0.224$ | 0.108 | 555 | 0.039 | 0.088 | 25.8 | 0.017 |

Table: Worked example 1: the interaction test under S and F1, model-based and CR2 inference, seed 20261009 (verified)

**Reading.**

- *With no interaction, neither analysis rejects*, and in this trial S
  shows no visible spurious interaction: its estimate, $-0.040$, is half
  a standard error from zero, and F1's is similar. This is the typical
  case. S's bias, $(1 - f)s_{\text{prog}}$, is a fraction of a standard
  error at $N = 70$, so a single trial cannot reveal it; it shows up as
  a doubled false-positive rate over many trials (Section 5) and as a
  clear offset in a large trial (Section 4.4).
- *With the interaction, both reject at 0.05.* F1's estimate is closer
  to the true $-0.156$ and has the larger standard error (0.088 against
  0.072 under CR2); that standard error is the price of not trusting
  the baseline row (Section 3.3).
- *Model-based and CR2 standard errors differ in both directions*, and
  the degrees of freedom differ twentyfold. In the trial with the
  interaction CR2 gives the smaller standard error; over many trials the two
  are close on average for these analyses ($\kappa$ near 0.94 to 1.03;
  Section 5), and CR2's value is that it remains right when the working
  covariance is wrong.
- *$\hat\beta_{bm:\eta}$*, the coefficient of `bmc:De`, is near zero
  without the interaction (0.030) and $-0.241$ (CR2 SE 0.172) with it:
  in a single trial it also picks up part of the predictive signal,
  because the open-label visits are on drug (verified,
  `56-example-coefficients.csv`).

### 4.4 A large trial: the bias made visible

The same construct at $N = 2000$ (500 per path, seed 20261010) shows
what the $N = 70$ trial hides (verified; fits took 170 seconds).

| Case | Analysis | $\hat\beta_{bm:D}$ | CR2 SE | $\hat\beta_{bm:D}/s_{\text{prog}}$ | Target |
|------------|------|------|------|------|------|
| $c_{bm} = 0$ | S | $-0.0565$ | 0.0150 | 0.289 | 0 |
| $c_{bm} = 0$ | F1 | $-0.0184$ | 0.0170 | 0.094 | 0 |
| $c_{bm} = 0.30$ | S | $-0.1908$ | 0.0140 | | $-0.156$ |
| $c_{bm} = 0.30$ | F1 | $-0.1464$ | 0.0156 | | $-0.156$ |

Table: Worked example 1 at $N = 2000$: the interaction estimate against its target (verified)

With no interaction, S reports an interaction of $-0.057$, 3.8 CR2
standard errors from zero ($p = 0.0002$), and its ratio to the
prognostic slope, 0.289, is close to the leak factor computed for S's
working covariance in the `04-power-simulation.R` construct with AR(1)
components, 0.276 (quoted;
compendium, tab:bp-leak). F1's estimate, $-0.018$, is 1.1 standard
errors from zero, and its ratio, 0.094, is close to that computed for
F1, 0.089: the small residual bias of Section 3.3, visible at this
size. With the interaction, S overshoots the target by $-0.035$, about
the size of its null leak, while F1 is within one standard error. The
1,000 stored null trials of this cell at $N = 70$ agree: mean
estimates $-0.051$ for S, $-0.017$ for F1 and $0.002$ for ES (verified,
`summary-full.csv`, cell 25).

### 4.5 The fitted slopes by visit type

Figure 1 plots, for each group of rows, the biomarker slope implied by
each fit, averaged over the rows of the group, against the true slope
of the construct and the slope observed in the trial.

![Worked example 1: slope of $Y$ on the biomarker by visit type, implied by the S fit (solid, triangles) and the F1 fit (dashed, circles), with the true slope of the construct (dotted, squares) and the slope observed within each visit group of the simulated trial (crosses). Top: the $N = 70$ trial; bottom: the $N = 2000$ trial. Left: no interaction; right: $c_{bm} = 0.30$. The biomarker is prognostic through the natural course, so its true slope is 0 at baseline and $-0.195$ at every later visit, plus the predictive part on drug. S must use one slope for baseline and the blinded off-drug visits; it compromises between them, and the on-drug visits absorb the rest as an interaction. F1 can bend at baseline and between open-label and blinded visits.](figures/56-fig1-fitted-slopes.png)

The bottom-left panel is the mechanism in one picture (verified). The
true slope is 0 at baseline and $-0.195$ after. S fits $-0.129$ at
baseline and $-0.134$ at the blinded off-drug visits (one shared slope,
moved slightly by the carryover in $D_{bc,it}$) and $-0.185$ on drug, so the difference, its interaction, is
$-0.051$ per row-average (the coefficient, $-0.057$, applies to the
full on-drug contrast). F1 fits $-0.082$ at baseline, $-0.224$ open
label, $-0.162$ blinded on drug and $-0.146$ blinded off drug: it still
misfits each group, but the blinded on-drug and off-drug slopes differ
by only $-0.017$.

## 5. Worked example 2: the stored Monte Carlo results

### 5.1 Size

Table 8 gives the null rejection rates of six analyses in each of the
15 null cells, 1,000 trials each, read from the stored summaries
(verified; `summary-full.csv`, `summary-full-extra.csv`,
`summary-full-extra2.csv`; pooled rates agree with
`decision-screen-ds-nosym-all14.csv` to $3 \times 10^{-5}$).

| Null cell | S | F1 | ES | F1A | DS | R |
|------------------------------------|-----|-----|-----|-----|-----|-----|
| AR(1) curve | 0.049 | 0.062 | 0.059 | 0.068 | 0.054 | 0.076 |
| CS curve | 0.048 | 0.057 | 0.048 | 0.069 | 0.050 | 0.050 |
| measurement error | 0.051 | 0.047 | 0.041 | 0.057 | 0.053 | 0.071 |
| $N = 35$ | 0.058 | 0.053 | 0.050 | 0.059 | 0.058 | 0.081 |
| $N = 150$ | 0.047 | 0.047 | 0.050 | 0.047 | 0.050 | 0.063 |
| $t_{1/2} = 0$ | 0.036 | 0.049 | 0.042 | 0.043 | 0.040 | 0.072 |
| $t_{1/2} = 1$ week | 0.043 | 0.061 | 0.050 | 0.073 | 0.055 | 0.061 |
| Weibull, shape 0.5 | 0.047 | 0.046 | 0.047 | 0.058 | 0.054 | 0.073 |
| Weibull, shape 2 | 0.037 | 0.045 | 0.044 | 0.040 | 0.040 | 0.066 |
| 20% MCAR dropout | 0.059 | 0.051 | 0.049 | 0.064 | 0.068 | 0.079 |
| 20% MAR dropout | 0.056 | 0.044 | 0.046 | 0.051 | 0.059 | 0.064 |
| **pooled, 11 standard cells** | **0.048** | **0.051** | **0.048** | **0.057** | **0.053** | **0.069** |
| prognostic: natural course, constant | 0.099 | 0.047 | 0.049 | 0.059 | 0.115 | 0.107 |
| prognostic: natural course, growing | 0.049 | 0.046 | 0.063 | 0.052 | 0.058 | 0.053 |
| prognostic: placebo response | 0.100 | 0.040 | 0.049 | 0.047 | 0.112 | 0.159 |
| prognostic: baseline severity | 0.050 | 0.040 | 0.050 | 0.045 | 0.055 | 0.075 |

Table: Null rejection rates by cell and analysis, stress test as run on 2026-10-07 (verified; MCSE 0.0069 per cell, 0.0021 pooled; criteria: pooled within 0.046 to 0.054, every cell at most 0.075)

**Reading.**

- *Only the analyses that free the baseline slope survive the
  prognostic cells.* S, DS and R double their false-positive rate when
  the biomarker tracks the natural course or the placebo response; F1
  and ES do not. A baseline-severity association enters every visit
  equally and harms no analysis materially (R reaches the cap, 0.075).
- *F1's largest cell is 0.062*, the AR(1) curve's null point, within
  the cap and about 1.7 MCSE above 0.05. Its $\kappa$ there is 0.94
  (verified, `56-kappa.csv`): CR2 is slightly small at $N = 70$ in this
  cell, as the small-sample theory allows.
- *F1 is right-sized when its half-life is wrong.* With the true
  half-life 0 or 1 week and $t_{1/2}^a$ fixed at 0.5, F1 rejects 0.049
  and 0.061; under Weibull decay 0.046 and 0.045 (Section 6.2).

### 5.2 Size-adjusted power

| Alternative cell | S | F1 | ES | F1A | G1 | G2 | DS | R |
|------------------------------------|-----|-----|-----|-----|-----|-----|-----|-----|
| AR(1) curve, $c_{bm} = 0.30$ | 0.514 | 0.364 | 0.364 | 0.338 | 0.354 | 0.380 | 0.486 | 0.398 |
| CS curve, $c_{bm} = 0.25$ | 0.554 | 0.522 | 0.662 | 0.446 | 0.546 | 0.528 | 0.480 | 0.596 |
| $N = 35$ | 0.222 | 0.192 | 0.206 | 0.218 | 0.214 | 0.218 | 0.232 | 0.208 |
| $N = 150$ | 0.828 | 0.748 | 0.704 | 0.748 | 0.752 | 0.754 | 0.822 | 0.734 |
| $t_{1/2} = 0$ | 0.594 | 0.496 | 0.504 | 0.522 | 0.510 | 0.532 | 0.588 | 0.496 |
| $t_{1/2} = 1$ week | 0.388 | 0.272 | 0.288 | 0.282 | 0.282 | 0.306 | 0.412 | 0.334 |
| measurement error | 0.400 | 0.318 | 0.336 | 0.296 | 0.320 | 0.330 | 0.402 | 0.382 |
| 20% MAR dropout | 0.418 | 0.382 | 0.356 | 0.372 | 0.404 | 0.378 | 0.424 | 0.378 |
| prognostic: natural course, constant | (0.606) | 0.482 | 0.370 | 0.448 | 0.520 | 0.524 | (0.590) | (0.522) |
| prognostic: natural course, growing | 0.516 | 0.378 | 0.192 | 0.386 | 0.484 | 0.490 | 0.522 | 0.448 |
| prognostic: placebo response | (0.652) | 0.446 | 0.498 | 0.434 | 0.428 | 0.424 | (0.610) | (0.590) |

Table: Size-adjusted power at $c_{bm} = 0.30$ (0.25 in the CS curve), $N = 70$ unless stated, 500 trials per cell (verified, `decision-power-ds-nosym-all14.csv`; MCSE at most 0.022). Values in parentheses belong to analyses whose test is biased in that cell; size adjustment corrects the level, not the bias

![Worked example 2. Left: null rejection rates in the 15 null cells for S, F1, ES, F1A and R; the shaded band is $0.05 \pm 2$ MCSE for one cell, the dashed line the cap of 0.075, and the dotted vertical line separates the 11 standard cells from the four prognostic ones. R's rate in the placebo cell, 0.159, is drawn at the top edge. Right: size-adjusted power along the AR(1) curve ($N = 70$; the stored curve has no null point, so lines join $c_{bm} > 0$ only, and $c_{bm} = -0.30$ is shown as an isolated point).](figures/56-fig2-size-power.png)

**Reading.**

- *The cost of F1 is real and depends on the serial correlation.* Under
  AR(1) data F1 loses 29% of S's power at the reference cell, and
  between 9% (MAR dropout) and 30% ($t_{1/2} = 1$ week)
  over ten AR(1) alternative cells (verified, `56-cost-F1-vs-S.csv`).
  Under compound symmetry the loss is 6% (0.522 against 0.554). The
  loss comes from no longer using the baseline row as a clean
  observation of the off-drug biomarker slope (`docs/45`, Section 6;
  `docs/48`, Section 1).
- *F1 keeps its power under every prognostic form*; ES loses most of
  it when the natural-course association grows over time (0.192).
- *G1 and G2 are effectively equivalent to F1* in the standard cells
  and somewhat better in the natural-course prognostic cells (0.520,
  0.524 against 0.482; 0.484, 0.490 against 0.378).
- *None of the robust analyses reaches power 0.80 on the grid at
  $N = 70$*; at $N = 150$, F1 reaches 0.748.

## 6. The analysis half-life

The full treatment of choosing $t_{1/2}^a$ from the data is in
`docs/57`. This section states only why the variant that was run
failed, why a fixed half-life does not, and what the project intends.

### 6.1 What F1A did, and why it failed size

F1A fitted F1 by maximum likelihood at each $t_{1/2}^a \in \{0, 0.25,
0.5, 1, 2\}$ weeks, took the value with the smallest AIC, refitted that
model by REML and tested $\beta_{bm:D}$ with CR2 as if the half-life had
been fixed (inspected, driver `fit_extra2()`). Because the candidates
differ only in the values of one regressor, they have the same number
of parameters, and the AIC rule is the maximum-likelihood estimate of
the half-life over the grid (derived; compendium, eq:invd-aic).

It is tempting because carryover is real and the data carry some
information about it. It fails because the selection and the test use
the same noise.

- *The likelihood includes the interaction term.* The half-life is
  chosen with `bmc:Dbc` in the model, so a value of $t_{1/2}^a$ that
  happens to line up the exposure regressor with chance biomarker
  structure in the drug contrast raises the likelihood and the
  $t$-statistic together (argued).
- *The test then ignores the choice.* The reported standard error is
  that of a model whose regressor was fixed, not estimated.

The stored paired results show the consequence (verified;
`replicates-full.rds` and `replicates-full-extra2.rds`, the same
trials, matched by cell and replicate; S agrees exactly between the
two runs).

| Selected $t_{1/2}^a$ | Trials | F1 size (fixed 0.5) | F1A size | Mean $\lvert t\rvert$, F1 | Mean $\lvert t\rvert$, F1A |
|------|------|------|------|------|------|
| 0 | 3,148 | 0.051 | 0.053 | 0.812 | 0.848 |
| 0.25 | 1,655 | 0.048 | 0.051 | 0.823 | 0.830 |
| 0.5 | 2,934 | 0.051 | 0.051 | 0.839 | 0.839 |
| 1 | 1,967 | 0.056 | 0.064 | 0.842 | 0.858 |
| 2 | 1,295 | 0.046 | 0.079 | 0.807 | 0.915 |
| all | 10,999 | 0.051 | 0.057 | 0.826 | 0.853 |

Table: F1A against F1 on the same null trials of the 11 standard cells, by the half-life F1A selected; complete pairs only (verified)

(Corrected 2026-10-09: the 0.5-week row previously gave 0.052 for
both sizes; the stored rate in `56-f1a-by-selected-h.csv` is
151/2,934 = 0.05147, which rounds to 0.051, so 0.052 was most likely
a second rounding of the three-digit printout 0.0515. The
table is restricted to complete pairs: of the 11,000 null trials, one
F1A refit failed after its selection (cell 2, replicate 641, selected
0.25 week) and is dropped from both columns, which gives 10,999 trials
and 1,655 at 0.25 week; F1 converged in all 11,000. `docs/57`, Table
2, now uses the same rule.)

When AIC moves the half-life away from 0.5 week, and especially to 2
weeks, the test rejects more often than F1 on the same trials. Over
the 10,999 trials, 135 are rejected by F1A only and 68 by F1 only. In
the null cell with true half-life 1 week, 77% of selections were 1 or
2 weeks and F1A rejected 0.073 (0.083 when it chose 2 weeks), against
F1's 0.061 (verified, `56-f1a-thalf1-by-h.csv`).

**How well AIC identifies the half-life.** Coarsely. Over all cells the
selections were spread 24%, 16%, 29%, 19% and 12% over 0, 0.25, 0.5, 1
and 2 weeks (verified; as in `docs/47`, Section 8), but that pools cells
with different true half-lives. By cell, AIC picked 0 in 87% of the
null trials when the true half-life was 0, and 1 or 2 weeks in 77% when
it was 1 week, but the true 0.5 week in only 38% of the null trials of
the base cell, and 1 or 2 weeks in 85% of the compound-symmetry curve's
null trials, whose true half-life is also 0.5 week (verified,
`56-stress-summary.R` output). The likelihood distinguishes "little
carryover" from "a week or more", not 0.25 from 0.5 or 1, and its
choice depends on the correlation structure.

**A small illustration.** Figure 3 shows the ML AIC profile over a fine
grid for three null trials of the base cell and three of the
$t_{1/2} = 1$ week cell, with the CR2 $t$-statistic of the interaction
across the same grid (verified; 40 seconds for the six profiles). The
AIC profiles are shallow, often within 2 units of their minimum over a
factor of four in the half-life, and the $t$-statistic changes
smoothly with $t_{1/2}^a$, by up to 1.1 over the grid in one trial. The
dashed profiles are those of F1 without the interaction term, the
proposed F1A0 rule (Section 6.3); in these trials they are nearly
identical to F1's, and on 50 further null trials per cell the two
rules chose the same grid value in 86% and 88% of trials (verified,
`56-aic-selection-summary.csv`; 164 seconds). This is an illustration
on a few trials, not a size estimate.

![AIC selection of the analysis half-life on six null trials. Top: ML AIC minus its minimum over the half-life grid, for F1 with the interaction (solid; the F1A rule) and without it (dashed; the proposed F1A0 rule), true half-life 0.5 week (left) and 1 week (right). The dotted line marks 2 AIC units. Bottom: the CR2 $t$-statistic of the interaction in the same trials across the grid; filled circles mark the half-life chosen by the F1 AIC, open circles that chosen by the F1A0 AIC. Colors identify trials.](figures/56-fig3-aic-profiles.png)

### 6.2 Why a fixed half-life is right-sized even when it is wrong

Under the null with a purely predictive biomarker, the biomarker has no
association with $Y$ that depends on exposure, so the true slope vector
has no component along $B_i^c D_{bc,it}$ for *any* regressor $D_{bc,it}$
chosen without looking at the outcome, and the interaction's limit is
zero whatever $t_{1/2}^a$ is (derived; compendium, eq:bp-plim). A wrong
half-life mis-scales the alternative, so it costs power; it does not
create false positives. The stress test agrees: F1 at the fixed 0.5
week rejected 0.049 when the true half-life was 0, 0.061 when it was 1
week, 0.046 and 0.045 under Weibull decay, and at most 0.062 in any
cell (verified, Table 8). With a prognostic biomarker the leak factors
of Section 3.3 depend mildly on $t_{1/2}^a$, through $D_{bc,it}$ in the
fitted slopes (derived), but no stress-test cell combined a prognostic
pathway with a wrong half-life.

### 6.3 Status and intention

The project's position, set out in `docs/57`:

- **Validated now:** F1 with $t_{1/2}^a$ fixed in advance from
  pharmacology, with a sensitivity analysis over a grid of half-lives.
- **Intended, pending evaluation:** F1 with the half-life chosen by the
  data, because a reader can fairly ask why the analysis should use a
  half-life the data contradict. Two procedures are proposed and not
  yet run. **F1A0** chooses $t_{1/2}^a$ by AIC on F1 *without* the
  interaction term, where carryover is identified from the drug-response
  trajectory, then refits with the interaction and tests with CR2; it
  is approximately, not exactly, independent of the interaction test
  under the null (argued). **F1AB** calibrates the selection-then-test
  procedure by a parametric bootstrap from the fitted null model,
  repeating the selection in every resample; its operating
  characteristics would be evaluated with the warp-speed method of
  Giacomini, Politis and White [7] (proposal).
- **Report the half-life estimate in either case,** with a
  profile-likelihood interval, and state that a Hybrid trial carries
  little information about it (Section 6.1).

## 7. Alternatives and when to prefer them

| Analysis | Prefer it when | Size, pooled / worst prognostic cell | Power, AR(1) / CS reference cell | Caveat |
|------|----------------------|----------|----------|------------------------------|
| F1 | the biomarker may be prognostic (the default) | 0.051 / 0.047 | 0.364 / 0.522 | approximately unbiased under a constant prognostic pathway; about a third less power than S under AR(1) data |
| G1, G2 | as F1, with phase-specific residual variances; effectively equivalent | 0.049 / 0.055; 0.047 / 0.062 | 0.354 / 0.546; 0.380 / 0.528 | two more covariance parameters; G2's model-based inference relies on the working covariance |
| ES | the prognostic pathway, if any, is constant after baseline; under CS-like data | 0.048 / 0.063 | 0.364 / 0.662 | leaks an expectancy-scaled slope; loses most power when the prognostic slope grows (0.192) |
| S | independent evidence shows the biomarker is not prognostic (or only through baseline severity) | 0.048 / 0.100 | 0.514 / 0.554 | doubles its false-positive rate under a natural-course or placebo pathway |
| R (published type) | not recommended | 0.069 / 0.159 | 0.398 / 0.596 | liberal in standard cells and prognostic cells |
| F1A | not recommended | 0.057 / 0.059 | 0.338 / 0.446 | fails size through post-selection inference (Section 6) |

Table: Alternatives to F1 (verified from the stored stress-test outputs; prognostic column: the largest of the four prognostic null cells)

Three remarks. First, the selection of F1 over G2 was a tie broken by
the smaller number of covariance parameters (`docs/47`, Section 7); if
all alternative cells, prognostic ones included, govern the selection,
G1 and G2 lead F1 by about 0.05 in largest shortfall (quoted, `docs/47`,
Section 7). Either is a defensible primary analysis. Second, report S
as a sensitivity analysis: if F1 and S agree, the conclusion does not
depend on the constancy assumption; if they disagree materially, the
disagreement is itself evidence of a prognostic association
(`docs/49`, Section 6). Third, in a fully blinded crossover the
expectancy is constant and the expectancy terms are unnecessary and
unidentified (`docs/49`, Section 8).

## 8. Practical recipe

The steps below fit F1 to a long data set `dat` with one row per
participant and visit, baseline included, and columns `ptID`
(participant), `t` (weeks from baseline, 0 at baseline), `Sx` (score),
`bm` (biomarker), `on_drug` (logical), `tsd` (weeks since
discontinuation: 0 on drug and before first exposure, cumulative
within an off-drug run) and `De` (expectancy: 0 at baseline, 1 open
label, 0.5 blinded).

1. **Fix $t_{1/2}^a$ before unblinding, from pharmacology.** For prazosin
   the pharmacological range is 0.1 to 0.2 week (report 01, Methods);
   the stress test used 0.5 week as its base setting, which is longer
   than the pharmacological range. Write down the grid for the
   sensitivity analysis (for example 0, 0.25, 0.5, 1 and 2 weeks).
2. **Build the exposure indicator with the pre-exposure guard**, and
   check that no participant who starts on placebo is coded exposed.
3. **Center the biomarker** at its baseline mean. Do not dichotomize it.
4. **Fit the mean model with the expectancy terms** and the working
   covariance, by REML.
5. **Test with CR2 and Satterthwaite degrees of freedom**, clustered on
   participant, two-sided at 0.05.
6. **Report** the estimate per unit (and per SD) of the biomarker, its
   CR2 standard error, degrees of freedom, $p$-value and 95% interval;
   the coefficient of `bmc:De` as an exploratory measure of a
   prognostic or expectancy-linked association; S and ES as sensitivity
   analyses; F1 over the half-life grid; and the data's own estimate of
   the half-life (the ML AIC choice of F1 without the interaction term)
   with the statement that it is weakly identified. State the
   serial-correlation assumption behind the sample size.

```{.r}
library(nlme)
library(clubSandwich)

exposure_decay <- function(on_drug, tsd, t_half) {
  if (t_half <= 0) return(as.numeric(on_drug))
  ## tsd is 0 before first exposure too; without the guard a
  ## placebo-first participant is coded fully exposed
  ifelse(on_drug, 1, ifelse(tsd > 0, 0.5^(tsd / t_half), 0))
}

t_half_a <- 0.5
dat <- dat |>
  transform(bmc = bm - mean(bm[t == 0]),
            Dbc = exposure_decay(on_drug, tsd, t_half_a))
stopifnot(all(dat$Dbc[!dat$on_drug & dat$tsd == 0] == 0))

fit <- lme(Sx ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De,
           random = ~ 1 | ptID,
           correlation = corCAR1(form = ~ t | ptID),
           data = dat, method = 'REML',
           control = lmeControl(opt = 'optim', maxIter = 200,
                                msMaxIter = 200))
coef_test(fit, vcov = 'CR2', test = 'Satterthwaite') |>
  subset(Coef == 'bmc:Dbc')
conf_int(fit, vcov = 'CR2', level = 0.95) |>
  subset(Coef == 'bmc:Dbc')
```

Applied to the $c_{bm} = 0.30$ trial of Section 4, this returns
$\hat\beta_{bm:D} = -0.224$, SE 0.0877, $\nu = 25.8$, $p = 0.0167$ and
the interval $(-0.405, -0.044)$, identical to the driver's F1 (verified;
`56-recipe` check described in Section 10).

**With the patched pmsimstats.** On `rchendrickson/pmsimstats` at
`58b32a9` with `analysis/patches/orig-to-refined/` applied, the same
analysis is the default of `lme_analysis()`: `op$analysis = 'F1'`,
`op$t_half_a = 0.5`, and the result carries the estimate, the CR2
standard error, the $p$-value and the degrees of freedom (inspected;
patch README, analysis table; verified there identical to the driver
on five data sets).

**Do not:** choose $t_{1/2}^a$ by AIC with the interaction in the model
and test as if fixed (F1A); use model-based standard errors with a
compound-symmetry working covariance (R); drop the guard in
$D_{bc,it}$; or plan the sample size for S and analyze with F1.

## 9. Limitations

- **Approximate unbiasedness.** F1 removes an expectancy-scaled
  prognostic slope exactly but a constant one only approximately in the
  Hybrid design, with residual leak factors of $-0.02$ to 0.17 of the
  prognostic slope depending on the working covariance (Section 3.3).
  The stress test simulated prognostic couplings of 0.3 only, one
  pathway at a time; larger couplings, combined pathways and a
  prognostic pathway combined with a wrong half-life were not run.
- **The stored results predate the refined process.** They used the
  published recursive carryover and a fixed allocation; the rerun
  under anchored carryover and randomized allocation is pending.
  Anchored and recursive means differ by up to 0.48 points at
  $t_{1/2} = 0.5$ week in the Hybrid design (quoted; patch README), and
  the published-code analysis was insensitive to that change
  (`docs/53`, Section 8.3), but the stress-test rates have not been
  recomputed.
- **One design family, one effect size.** Only the Hybrid design was
  stressed, with the alternative at $c_{bm} = 0.30$ (0.25 under
  compound symmetry), and $N$ from 35 to 150.
- **Gaussian outcome.** The construct is multivariate normal; skewed
  or bounded scores, and a biomarker acting on the rate rather than the
  size of response, are outside it.
- **Exploratory selection.** The study's criteria and analyses were
  revised as results came in (`docs/46`, Amendments 1 and 2), so F1 is
  a recommendation to be confirmed in an independent simulation.
- **The worked examples are single trials.** Section 4's $N = 70$
  numbers illustrate the computation; the claims about bias rest on the
  large trial and the stored Monte Carlo results.
- **The data-chosen half-life is not yet evaluated** (Section 6.3,
  `docs/57`).

## 10. Reproducibility

Every number labeled verified in Sections 3 to 8 was computed by the
scripts below, in `analysis/scripts/` `recommended-analysis/`, run from
the repository root with R 4.6, `nlme` 3.1-169 and `clubSandwich` 0.7.0
(verified, host session; the container was not used).

- `56-load-driver.R`: loads the data-generating process and analyses of
  the stress-test driver (`analysis/scripts/quick-sim/`
  `carryover-closed-form/` `11-strawman-stress-test.R`) by evaluating its top-level definitions only, with
  `--carryover recursive --allocation fixed` semantics, so that every
  fit is the driver's own code.
- `56-example-trial.R`: Table 3 and the crossover guard, the design
  weights and the drug-expectancy correlation, Tables 5 to 7 and
  Figure 1 (about 3 minutes, most of it the $N = 2000$ fits).
- `56-stress-summary.R`: Tables 8 to 10, the CR2 versus model-based
  comparison, $\kappa$, the F1A selection counts and Figure 2, from the
  stored files in
  `analysis/data/quick-sim/` `carryover-closed-form/strawman-stress/`
  (`summary-full*.csv`, `replicates-full.rds`,
  `replicates-full-extra2.rds`, `decision-*-ds-nosym-all14.csv`); no
  new simulation (a few seconds).
- `56-aic-illustration.R`: Figure 3 and the F1A versus F1A0 selection
  comparison on 2 by 50 null trials (about 3.5 minutes).

For example, the command

```bash
Rscript analysis/scripts/recommended-analysis/56-example-trial.R
```

writes its tables to `analysis/data/recommended-analysis/` and the
figure to `docs/figures/`. The recipe of Section 8 was run on the
$c_{bm} = 0.30$ trial of Section 4 and matched the driver's F1 to all
printed digits (verified, a scratch check not kept in the repository).

## 11. References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker
   validation: statistical methods and theoretical findings.
   *Frontiers in Digital Health*. 2020;2:13.
   doi:10.3389/fdgth.2020.00013.
2. Hendrickson RC. pmsimstats: R package for simulation of aggregated
   N-of-1 trials. GitHub repository rchendrickson/pmsimstats, commit
   58b32a92b99a8beb145de00c7bd8a0880f06f0f0 (2020-11-19).
   https://github.com/rchendrickson/pmsimstats.
3. Bell RM, McCaffrey DF. Bias reduction in standard errors for linear
   regression with multi-stage samples. *Survey Methodology*.
   2002;28(2):169-181. (No DOI; published by Statistics Canada.)
4. Pustejovsky JE, Tipton E. Small-sample methods for cluster-robust
   variance estimation and hypothesis testing in fixed effects models.
   *Journal of Business and Economic Statistics*. 2018;36(4):672-683.
   doi:10.1080/07350015.2016.1247004.
5. Liang KY, Zeger SL. Longitudinal data analysis using generalized
   linear models. *Biometrika*. 1986;73(1):13-22.
   doi:10.1093/biomet/73.1.13.
6. Kenward MG, Roger JH. Small sample inference for fixed effects from
   restricted maximum likelihood. *Biometrics*. 1997;53(3):983-997.
   doi:10.2307/2533558.
7. Giacomini R, Politis DN, White H. A warp-speed method for conducting
   Monte Carlo experiments involving bootstrap estimators.
   *Econometric Theory*. 2013;29(3):567-589.
   doi:10.1017/S0266466612000655.
8. Liu GF, Lu K, Mogg R, Mallick M, Mehrotra DV. Should baseline be a
   covariate or dependent variable in analyses of change from baseline
   in clinical trials? *Statistics in Medicine*. 2009;28(20):2509-2530.
   doi:10.1002/sim.3639.
9. Senn S. *Cross-over Trials in Clinical Research*. 2nd ed.
   Chichester: Wiley; 2002. doi:10.1002/0470854596.
10. Morris TP, White IR, Crowther MJ. Using simulation studies to
    evaluate statistical methods. *Statistics in Medicine*.
    2019;38(11):2074-2102. doi:10.1002/sim.8086.
11. pmsimstats team. `docs/41`, `docs/42`, `docs/43`, `docs/45`,
    `docs/46`, `docs/47`, `docs/48`, `docs/49`, `docs/53`, `docs/57`;
    `analysis/report/` `compendium-mathematics.tex`; report 01
    (`analysis/report/` `01-dgp-mean-moderation-vs-mvn/report.Rmd`).
