# Choosing the analysis half-life from the data

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** Why AIC with interaction inflates size; eight methods; F1A0B evidence; extension results.  
**Prerequisites:** ch22

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

<!-- table: process, carryover, coupling, allocation, design, N, analysis ID,
t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 57 (spine; Sections 6 and 7.3 moved to an annex); 47 Section 8; 56
  Section 6; CM "AIC selection of the analysis half-life", "The assumed
  half-life of $D_{bc}$"; Appendix I Section 9 (`characterize_carryover`).

- New to write: F1A0, the bootstrap, max-$\lvert t\rvert$ and the
  orthogonality result moved into CM; the full-bootstrap check (GP-15) when
  run.
-->

<!-- Spine: docs/57-data-chosen-carryover-half-life.md (copied verbatim; Phase 2) -->

---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# Testing the Biomarker-by-Drug Interaction When the Carryover Half-Life Is Chosen From the Data {.unlisted .unnumbered}
*2026-10-09 09:54 PDT*

**Author.** pmsimstats team

**Purpose.** The recommended analysis of the biomarker-by-drug
interaction in an aggregated N-of-1 (Hybrid) trial, F1 (`docs/56`),
weights each off-drug visit by an exposure-decayed drug indicator
$D_{bc,it}$ whose half-life $t_{1/2}^a$ is set by the analyst. With
$t_{1/2}^a$ fixed in advance the test holds its size, even when the
value is wrong. Face validity, however, asks that the analysis not use
a half-life the trial's own data contradict, and the project decided
on 2026-10-09 that the half-life should be chosen from the data. This
note sets out how to do that without losing control of the Type I
error. It explains why the one data-driven variant already run (F1A,
AIC over a grid with the interaction in the model) is liberal, derives
the mechanism, and presents eight procedures and a permutation
approach, each with its idea, an exact procedure, a validity argument
and its conditions, its cost, what it reports and its status. Two
worked examples, one trial and a small Monte Carlo illustration on
stored stress-test trials, show what each procedure does. A decision
table and a provisional plan close the note. It is written for trial
statisticians and methodologists, including the authors of the
original `pmsimstats` work [1], and assumes no other pmsimstats-ng
document beyond the model statement of `docs/56`.

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
here for this note only. (Notation, 2026-10-09: the centered,
unscaled biomarker of the code column `bmc`, formerly written $b_i$
and, in Section 4.2, $\tilde b_i$, is now the canonical $B^c_i$; $b_i$
is kept for the standardized biomarker. The time term is
$\beta_T w_t$, the normal law $\mathcal N$, and the compound-symmetry
correlation $\rho_{cs}$; probability is $\Pr$ ($P$ is the number of
paths). Symbols previously marked local or proposed
that `NOTATION.md` now lists are marked canonical.)

| Symbol | Meaning | Status |
|------------|--------------------------------------------------|----------|
| $Y_{it}$ | symptom score of participant $i$ at visit $t$; $t = 0$ is the baseline row (code `Sx`) | canonical |
| $B_i$, $B^c_i$ | biomarker; biomarker centered at the trial's baseline mean, $B_i - \bar B$, not scaled (code `bmc`). Test statistics are invariant to the scaling, so the results hold equally for the standardized $b_i$ | canonical |
| $D_{it}$ | binary drug state: 1 on drug, 0 off | canonical |
| $t_{sd}$ | time since discontinuation, cumulative within an off-drug run, 0 on drug and before first exposure (code `tsd`) | canonical |
| $t_{1/2}$ | true carryover half-life of the data-generating process, in weeks | canonical |
| $t_{1/2}^a$ | half-life assumed by the analysis in $D_{bc,it}$, in weeks | canonical |
| $h$ | a candidate value of $t_{1/2}^a$, used as a function argument | local |
| $D_{bc,it}(h)$ | exposure-decayed drug indicator at analysis half-life $h$: 1 on drug, $2^{-t_{sd}/h}$ off drug after first exposure, 0 before it; $h = 0$ gives $D_{it}$ (code `Dbc`, `X` in the driver) | canonical |
| $\eta_t$ | expectancy: 1 open label, 0.5 blinded, 0 at baseline (code `De`) | canonical |
| $\beta_{bm:D}$ | the biomarker-by-drug interaction, the estimand | canonical |
| $\beta_D$ | coefficient of $D_{bc,it}$, the drug main effect | canonical |
| $\beta_0$, $\beta_B$, $\beta_T$ | intercept, biomarker main effect and linear-time coefficient | canonical |
| $\beta_\eta$, $\beta_{bm:\eta}$ | coefficients of $\eta_t$ and $B^c_i\eta_t$ | canonical |
| $u_i$, $\sigma_u^2$ | participant random intercept and its variance | canonical |
| $\varepsilon_{it}$, $\sigma_\varepsilon$ | observation-level residual and its standard deviation | canonical |
| $\rho_a$ | working CAR(1) correlation per week of the analysis model | canonical |
| $\mathcal T_a$ | the five-point grid $\{0, 0.25, 0.5, 1, 2\}$ weeks (compendium, eq:invd-aic) | local |
| $\ell_1(h)$, $\ell_0(h)$ | maximized ML log-likelihood of F1 at $h$, with and without the interaction term | local |
| $\hat h_1$, $\hat h_0$ | $\arg\max_h \ell_1(h)$ (the F1A choice) and $\arg\max_h \ell_0(h)$ (the F1A0 choice) | local |
| $T(h)$ | CR2 $t$-statistic for $\beta_{bm:D}$ in F1 refitted by REML at $h$ | local |
| $\mathrm{LR}(h)$ | $2\{\ell_1(h) - \ell_0(h)\}$, the likelihood-ratio statistic for the interaction at $h$ | local |
| $\mathrm{AIC}_1(h)$, $\mathrm{AIC}_0(h)$ | AIC of F1 at $h$ with and without the interaction term | local |
| $a_i(h, Y)$ | participant $i$'s contribution to the exposure contrast after partialling (Section 4.2) | local |
| $g(Y)$ | a half-life selection rule that uses the outcome and the design only | local |
| $h_A$, $h_B$ | half-lives chosen on the two halves of a split sample | local |
| $\omega$ | $\log t_{1/2}^a$, the half-life as a continuous parameter (Section 4.6) | local |
| $D_{bc,i}$, $V_i$ | participant $i$'s vector of $D_{bc,it}$ and working covariance | local |
| $I_{\beta\omega}$, $U_\beta$ | blocks of the information matrix and of the score | local |
| $q$, $c_q$ | randomization-path index (canonical), and a path-level constant (Section 4.6, local) | canonical, local |
| $W_i$, $\hat\xi$ | participant $i$'s null-model design matrix and the fitted null-model fixed effects (Section 4.5) | local |
| $c_{bm}$ | biomarker-BR coupling in the data-generating process (the predictive effect) | canonical |
| $\rho$ | AR(1) serial correlation per week | canonical |
| $N$ | total sample size across paths | canonical |
| $M$, $m$ | number of bootstrap resamples per trial, and the resample index | local |
| $T^{*}_m$, $p^{*}$ | bootstrap statistic and bootstrap $p$-value | local |
| $n_{\text{sim}}$ | number of Monte Carlo replicates in a cell | canonical |
| $c_{bm,TV}$ | biomarker-natural-course coupling of a prognostic cell | canonical |
| $c_{bm,BL}$ | biomarker-baseline-level coupling of the prognostic-baseline cell (code `c_bl`) | canonical |
| $\mathrm{BL}_i$, $\sigma_{BL}$, $\sigma_{bm}$ | baseline level of participant $i$, its standard deviation, and that of the biomarker | canonical |
| $\beta_{B,0}$ | extra biomarker slope at baseline (code `bmc:bl0`) | canonical |
| $\mu_{BL}$, $\mu_B$, $\gamma$, $e_i$ | population means of $\mathrm{BL}_i$ and $B_i$, the regression slope of $\mathrm{BL}_i$ on $B_i$, and its residual (Section 8.5) | local |
| $\alpha$ | nominal level, 0.05 | canonical |

Table: Notation index for this note

**Glossary.**

- **Hybrid design.** Eight visits at weeks 4, 8, 9, 10, 11, 12, 16 and
  20 plus baseline, open-label titration, blinded discontinuation and
  a brief crossover, four randomization paths [1].
- **F1.** The recommended analysis (`docs/56`): linear time, $D_{bc,it}$
  at $t_{1/2}^a$, the expectancy terms $\eta_t$ and $B^c_i\eta_t$, a
  random intercept, continuous-time AR(1) residuals (`corCAR1`) and
  CR2 cluster-robust inference with Satterthwaite degrees of freedom
  [3, 4].
- **F1A.** F1 with $t_{1/2}^a$ chosen by AIC over $\mathcal T_a$, ML
  for the choice, REML refit, then tested as if fixed (the stress-test
  variant of `docs/47`).
- **F1A0.** F1 with $t_{1/2}^a$ chosen by AIC over $\mathcal T_a$ on F1
  *without* the interaction term, then refitted with it (proposed in
  `docs/56`; illustrated here).
- **F1AB.** F1A, or F1A0, calibrated by a parametric bootstrap from the
  fitted null model that repeats the selection in every resample.
- **Size.** The null rejection rate at nominal $\alpha$; the stress
  test caps any cell at 0.075 and asks the pooled rate to lie in the
  2-MCSE band around 0.05 (`docs/47`).
- **Base cell.** The stress-test cell with $N = 70$, AR(1) at
  $\rho = 0.7$, $t_{1/2} = 0.5$ week, no dropout and no prognostic
  pathway (`docs/46`, `docs/47`).
- **Working parameter.** A parameter of the analysis model that need
  not correspond to any parameter of the data-generating process; here
  $t_{1/2}^a$ (Section 3.2).

**Evidence labels.** Each claim carries one of five labels: *quoted*
(text reproduced from a source), *inspected* (read from source code or
stored output), *derived* (follows algebraically, possibly to first
order, as stated), *verified* (computed by running code, with the
script named in Section 11) or *argued* (a reasoned claim that has not
been derived in full or checked by computation).

## 1. Summary for the busy reader

**The question.** F1 tests $\beta_{bm:D}$ with $D_{bc,it}$ built at a
half-life $t_{1/2}^a$. Fixed in advance at 0.5 week, F1 rejected 0.051
of 11,000 null trials pooled over the 11 non-prognostic null cells,
and at most 0.062 in any cell, including the cells whose true
half-life was 0 or 1 week (quoted from `docs/47`, Section 3; stored
setting: configuration C at $\rho = 0.7$, the published recursive
carryover and fixed allocation; an anchored rerun is pending). The
project now wants the half-life chosen from the data.

**Why the obvious way fails.** F1A chose $t_{1/2}^a$ by AIC with the
interaction in the model and rejected 0.057 pooled (outside the band
0.046 to 0.054) and 0.073 at true $t_{1/2} = 1$ week (quoted,
`docs/47`). Because the candidates have equal parameter counts,
$\mathrm{AIC}_1(h) = \mathrm{AIC}_0(h) - \mathrm{LR}(h)$: the rule
rewards the half-life at which the interaction looks strongest, and
the test then ignores that it was chosen (derived, Section 3.4). On
the stored trials, when F1A chose 2 weeks it rejected 0.079 against
F1's 0.046 on the same trials (verified).

**What the data say about the half-life.** Coarsely something. AIC
chose 0 in 87% of null trials when the truth was 0, and 1 or 2 weeks
in 77% when it was 1 week, but the true 0.5 week in only 38% of
base-cell trials, and 1 or 2 weeks in 85% of compound-symmetry trials
whose truth was also 0.5 (verified from stored replicates). In the
example trial the profile likelihood cannot separate 0 from 0.5 week
but rejects a week or more decisively (verified, Section 5). The
analysis half-life is also a working parameter, not the pharmacological
one (Section 3.2).

**What to do instead.** Choose the half-life by a criterion that does
not see the interaction (F1A0), or on data the test does not use
(sample splitting and cross-fitting), or calibrate the whole
select-then-test procedure (F1AB); report the half-life with a
profile-likelihood interval. Under the null with the biomarker
independent of the outcome, a selection that uses only the outcome is
harmless (derived to first order, Section 4.2), and the half-life is
asymptotically information-orthogonal to $\beta_{bm:D}$ when paths are
randomized independently of the biomarker (derived, Section 4.6).

**What the illustration shows.** On stored stress-test trials
(Section 6; illustration, not a size study; 1,000 pooled null trials,
MCSE about 0.007), F1 at the fixed 0.5 week rejected 0.061, F1A 0.074
and F1A0 0.067 on the same trials: F1A0 removes about half of F1A's
excess, not all. The procedures that do not commit to one selected
value held size: AIC-weight averaging 0.053 (0.049 with null-model
weights), cross-fitting with a Cauchy combination 0.050, and the
maximum $|t|$ over the grid with bootstrap calibration 0.048. A rule
that keeps 0.5 week unless the null-model profile rejects it (F1G,
post hoc) behaved as F1 (0.063) and switched in 1% to 15% of trials.
Power in the base alternative cell was 0.400 for F1 and F1G, 0.39 for
averaging and the bootstrap max-$|t|$, 0.29 for the Cauchy cross-fit
and 0.21 for a single split (verified).

**What the stress-test extension shows** (2026-10-10, Section 8; 58
cells, 40,000 fresh trials, refined process; verified). Pooled over the
12 non-prognostic null cells (band 0.046 to 0.054), F1 at the fixed
0.5 week rejects 0.053, F1 at 0.2 week 0.050, F1A 0.060 (fails again),
F1A0 0.053, AIC-weight averaging (MA) 0.0455 (just below the band) and
the grid maximum 0.018 with Bonferroni, 0.072 with a multivariate
normal reference (liberal) and 0.053 with a multivariate $t$. Every
warp-speed bootstrap arm rejects 0.048 to 0.051. The illustration's
F1A0 excess does not recur. Mean size-adjusted power over the 26
non-prognostic alternative cells is 0.38 to 0.40 for every rule that
passes; the gains are where 0.5 week is far from the truth, and they
belong to the bootstrap arms: at true $t_{1/2} = 2$ weeks and $c_{bm} =
0.45$, F1 has 0.352, F1AB 0.456 and F1A0B 0.500. Maximin shortfall is
0.076 for F1A0B, 0.078 for F1AB and 0.148 for F1 (0.062 for the post
hoc MA0B). The bootstrap arms reach 0.067 to 0.077 in the
prognostic-baseline null cell; that is Monte Carlo error, not a
failure of the null model (derived and verified, Section 8.5).

**Recommendation** (revised 2026-10-10, superseding the provisional
recommendation of 2026-10-09; see the dated note in Section 7.2). The
primary data-driven candidate is **F1A0B**: choose the half-life by AIC
without the interaction term (F1A0) and calibrate the test by the
parametric null-model bootstrap that repeats the selection in every
resample (Section 4.5). The alternative is **F1AB**, the same bootstrap
of the selection with the interaction in the model. **F1 at a fixed
half-life** stays the validated simple fallback. MA is slightly
conservative (0.0455) and its bootstrap MAB passes; MaxTMVN is liberal
and MaxTBonf too conservative; the pharmacological fixed 0.2 week
gives no power gain over 0.5. The bootstrap evidence is warp-speed
(one resample per trial; the binomial MCSE understates its Monte Carlo
error), and a real analysis needs about $M = 999$ resamples, about 17
worker-minutes for F1A0B at $N = 70$ (derived, Section 8.4). Report the
half-life estimate and its profile-likelihood interval in every case.
For the main effect the conclusions differ (`docs/58`): AIC selection
with the tested term in the model is exactly a max-$|t|$ test judged
against a single-test critical value, and a fixed half-life does not
by itself guarantee size unless time is modeled with visit effects.

(Superseded 2026-10-10. The provisional recommendation of 2026-10-09
made AIC-weight averaging with $\ell_1$ weights the primary candidate,
the bootstrap max-$|t|$ test the calibration reference and F1A0 a
partial fix, on the illustration of Section 6. At scale, MA falls just
below the size band, F1A0 holds size, and the bootstrap calibration
carries the power gain.)

## 2. Scope and relation to `docs/56`

`docs/56` states the recommended model F1, its expectancy terms, the
pre-exposure guard on $D_{bc,it}$, its covariance and its CR2 inference,
and summarizes the half-life issue in its Section 6. This note is the
full treatment of that issue and does not repeat the model's
justification. It concerns only how $t_{1/2}^a$ is chosen and how the
test of $\beta_{bm:D}$ accounts for the choice. It does not revisit the
expectancy terms, the covariance structure or the robustness of F1 to
a prognostic biomarker, except where a half-life procedure interacts
with them. All numerical work is on the Hybrid design at $N = 70$.

## 3. The problem

### 3.1 Where the half-life enters

F1 (`docs/56`, Section 2) is the linear mixed model

$$
Y_{it} = \beta_0 + \beta_B B^c_i + \beta_T w_t + \beta_D D_{bc,it}(h)
  + \beta_{bm:D}\, B^c_i D_{bc,it}(h) + \beta_\eta \eta_t
  + \beta_{bm:\eta}\, B^c_i \eta_t + u_i + \varepsilon_{it},
$$

with $u_i \sim \mathcal N(0, \sigma_u^2)$, $\varepsilon_{it}$ continuous-time AR(1)
in weeks, the baseline as a response row at $t = 0$ with
$D_{bc,it} = \eta_t = 0$, and $H_0: \beta_{bm:D} = 0$ tested with the CR2
standard error and Satterthwaite degrees of freedom [3, 4]. The
half-life enters through one regressor, $D_{bc,it}(h)$, in two places:
the drug main effect and the interaction. In R, as in the stress-test
driver (`11-strawman-stress-test.R`, lines 203 to 206 and 239;
inspected):

```{.r}
dbc <- function(Db, tsd, h) {
  if (h <= 0) return(as.numeric(Db))
  ifelse(Db, 1, ifelse(tsd > 0, 0.5^(tsd / h), 0))
}
F_F1 <- Sx ~ bmc + t + X + bmc:X + De + bmc:De      # X = dbc(Db, tsd, h)
```

At $h = 0$ the regressor is the binary $D_{it}$; as $h \to \infty$ it
becomes 1 at every visit after first exposure. In the Hybrid design
the half-life changes $D_{bc,it}$ only at the off-drug visits, all of
which follow exposure: weeks 11, 12 and 20 on path 1; 11, 12 and 16 on
path 2; 10, 11, 12 and 20 on path 3; and 10, 11, 12 and 16 on path 4
(inspected, the `designs` list of
`hendrickson-problems/04-power-simulation.R`, which the driver
reads). The fits use `nlme::lme` with `corCAR1` [2].

### 3.2 The analysis half-life is a working parameter

$t_{1/2}^a$ and $t_{1/2}$ are different objects (argued). In the
data-generating process the drug-response mean rises along a modified
Gompertz curve while on drug and decays from its last on-drug value
off drug [1]; the analysis replaces the whole on-drug curve with a
constant $\beta_D$ and the off-drug decay with $\beta_D 2^{-t_{sd}/h}$.
The value of $h$ that best fits a given trial therefore reflects the
on-drug shape, the visit schedule and the correlation structure as
well as the pharmacology. In the stored stress test, which used the
published recursive carryover that decays faster than its stated
half-life (`docs/53`), AIC chose 1 or 2 weeks in 85% of null trials of
the compound-symmetry cell although the stated half-life was 0.5 week
(verified, Table 3). A data-chosen $t_{1/2}^a$ is an estimate of the
best-fitting exposure regressor, not of the drug's carryover
half-life, and should be reported as such.

### 3.3 Why a fixed half-life keeps its size even when it is wrong

Under the null with a purely predictive biomarker ($c_{bm} = 0$), the
biomarker is unrelated to how the outcome depends on exposure. For any
regressor $D_{bc,it}(h)$ fixed before the outcome is seen, the true
mean has no component along $B^c_i D_{bc,it}(h)$ that is not already
absorbed by the other terms, so the probability limit of
$\hat\beta_{bm:D}$ is 0 for every $h$ (derived; compendium,
eq:bp-plim, and `docs/56`, Section 6.2). A wrong $h$ mis-scales the
alternative: it moves the estimand and costs power, but it does not
create false positives. The stress test agrees (quoted, `docs/47`,
Section 3, and `docs/56`, Section 6.2): F1 at 0.5 week rejected 0.049
when the true half-life was 0 and 0.061 when it was 1 week. With a
prognostic biomarker the leakage terms of `docs/56` depend mildly on
$h$; no stored cell combined a prognostic pathway with a wrong
half-life.

### 3.4 Why choosing the half-life by AIC inflates size

**The identity.** Every candidate model has the same parameters, so
AIC differences are likelihood differences (compendium, eq:invd-aic).
Writing $\ell_1(h)$ and $\ell_0(h)$ for the ML log-likelihoods with and
without the interaction,

$$
-\tfrac12\mathrm{AIC}_1(h) + \text{const}
  = \ell_1(h) = \ell_0(h) + \tfrac12 \mathrm{LR}(h),
\qquad \mathrm{LR}(h) \approx T(h)^2
$$

(derived; the approximation is the usual equivalence of the likelihood
ratio and the squared Wald statistic, here with a CR2 rather than a
model-based standard error). In the example trial of Section 5,
$\mathrm{LR}(0) = 0.678$ against $T(0)^2 = 0.641$, and
$\mathrm{LR}(2) = 0.0016$ against $T(2)^2 = 0.0016$ (verified). F1A
therefore maximizes

$$
\ell_0(h) + \tfrac12 T(h)^2 \quad\text{approximately:}
$$

the drug-trajectory information in $\ell_0$ plus a bonus for a large
interaction statistic. It then compares $T(\hat h_1)$ with the
reference distribution of a test at a fixed $h$, the general
post-selection problem [16]. Two limiting cases
bracket its behavior (derived from the display):

- **$\ell_0$ steep.** If the drug main effect pins the half-life, so
  that differences in $\ell_0$ across $\mathcal T_a$ are large compared
  with the variation of $T(h)^2/2$, then $\hat h_1 = \hat h_0$ in most
  trials and F1A behaves like F1A0.
- **$\ell_0$ flat.** If the candidate regressors are nearly
  indistinguishable through the drug main effect, then
  $\hat h_1 \approx \arg\max_h |T(h)|$ and F1A becomes the max-$|t|$
  test of Section 4.7 referred to a single-$t$ critical value. Its
  size is then $\Pr(\max_h |T(h)| > t_{0.975})$, which exceeds
  $\alpha$ by an amount set by how correlated the five $T(h)$ are.

The inflation is mild because the $T(h)$ are strongly correlated
across $h$ (in the example trial they move from $-0.80$ to $-0.04$
over the grid, smoothly; Figure 1), and worst where $\ell_0$ is flat,
that is, where the half-life is weakly identified (argued).

**The stored evidence.** Table 2 compares F1A with F1 on the same
null trials, by the half-life F1A selected.

| Selected $t_{1/2}^a$ (weeks) | Trials | F1 at 0.5 | F1A |
|------|------|------|------|
| 0 | 3,148 | 0.051 | 0.053 |
| 0.25 | 1,655 | 0.048 | 0.051 |
| 0.5 | 2,934 | 0.051 | 0.051 |
| 1 | 1,967 | 0.056 | 0.064 |
| 2 | 1,295 | 0.046 | 0.079 |
| all | 10,999 | 0.051 | 0.057 |

Table: Null rejection rates of F1 (fixed 0.5 week) and F1A on the same trials of the 11 non-prognostic null cells, by the half-life F1A selected; complete pairs only (verified from `replicates-full.rds` and `replicates-full-extra2.rds`; stored setting)

(Corrected 2026-10-09: the table previously gave 11,000 trials and
1,656 at 0.25 week. Of the 11,000 null trials, one F1A fit failed
after its selection (cell 2, replicate 641, selected 0.25 week; F1A
$p$-value missing), while F1 converged in all 11,000. The earlier
version kept that trial in the count and in F1's rate but dropped it
from F1A's, so the two columns did not refer to the same trials. The
table now uses the 10,999 complete pairs, as `docs/56` does; no
rounded rate changes. F1's pooled rate over all 11,000 trials, quoted
in Section 1, is 0.051 either way.)

When F1A keeps 0.5 week it is F1. When it moves away, and most of all
to 2 weeks, it rejects more often than F1 on the same trials; F1's own
rate in those trials stays near 0.05, so the excess comes from the
selection, not from the trials (verified). `docs/56`, Section 6.1,
reports the same comparison with the mean $|t|$.

### 3.5 How much the data say about the half-life

| Null cell | True $t_{1/2}$ | 0 | 0.25 | 0.5 | 1 | 2 |
|------|------|------|------|------|------|------|
| base (AR(1), $\rho = 0.7$) | 0.5 | 0.232 | 0.230 | 0.379 | 0.121 | 0.038 |
| $t_{1/2} = 0$ | 0 | 0.870 | 0.082 | 0.044 | 0.003 | 0.001 |
| $t_{1/2} = 1$ week | 1 | 0.012 | 0.027 | 0.194 | 0.417 | 0.350 |
| compound symmetry, $\rho_{cs} = 0.8$ | 0.5 | 0.019 | 0.025 | 0.110 | 0.382 | 0.464 |
| $N = 35$ | 0.5 | 0.317 | 0.147 | 0.286 | 0.146 | 0.104 |
| $N = 150$ | 0.5 | 0.138 | 0.286 | 0.479 | 0.092 | 0.005 |
| Weibull decay, shape 2 | 0.5 | 0.710 | 0.167 | 0.107 | 0.012 | 0.004 |
| all 41 cells | | 0.236 | 0.163 | 0.292 | 0.190 | 0.120 |

Table: Share of null trials in which F1A's AIC chose each analysis half-life (weeks), selected cells, 1,000 trials each, and pooled over all cells and hypotheses (verified; stored setting)

Three readings (verified unless stated):

- **The pooled split** of `docs/47` (24, 16, 29, 19 and 12%) mixes
  cells with different truths and should not be read as "AIC cannot
  find the half-life".
- **The likelihood separates "little carryover" from "a week or
  more"**, not 0.25 from 0.5 or 1: 87% at 0 when the truth is 0, 77%
  at 1 or 2 weeks when it is 1 week, but the true 0.5 only 38% of the
  time in the base cell. At $N = 150$ the share at 0.5 rises to 0.48.
- **The choice depends on things other than the half-life**: the
  correlation structure (compound symmetry pushes it to 1 or 2 weeks)
  and the decay shape (a Weibull shape-2 decay with the same half-life
  is chosen as 0 in 71% of trials). This is the working-parameter
  point of Section 3.2 (argued).

## 4. Methods

Each method is described by its idea, the intuition, an exact
procedure, a validity argument with its conditions, its cost, what it
reports, and its status. *Evaluated* means run at the stress-test
scale (`docs/47`); *illustrated* means run only in the small Monte
Carlo of Section 6; *proposed* means not run. Costs are on the
machine used for `docs/47` (a 4 performance plus 4 efficiency core
machine): one F1 fit with CR2 takes a fraction of a second, and the
five-point grid with and without the interaction (15 fits and 5 CR2
computations) about 1.1 seconds on one worker (verified).

The code below uses the helpers of
`analysis/scripts/data-chosen-half-life/57-common.R`: `with_x(d, h)`
adds `X` $= D_{bc,it}(h)$, `fit_lme(f, d, method)` fits F1-type models
with a random intercept and `corCAR1`, and `cr2t(fit)` returns the
CR2 estimate, standard error, $t$, Satterthwaite df and $p$ for
`bmc:X`.

```{.r}
FM1 <- Sx ~ bmc + t + X + bmc:X + De + bmc:De   # F1
FM0 <- Sx ~ bmc + t + X + De + bmc:De           # F1 without the interaction
G5  <- c(0, 0.25, 0.5, 1, 2)                    # the grid T_a, weeks
ll  <- function(f, d, h) as.numeric(logLik(fit_lme(f, with_x(d, h), 'ML')))
```

### 4.1 Method 1: fixed in advance, with a sensitivity grid

**Idea.** Set $t_{1/2}^a$ before the data are seen, test at that
value, and report F1 over the grid as a sensitivity analysis. In this
note the fixed value is 0.5 week, the base setting of the stress test
(`docs/46`). It is not prazosin's pharmacological half-life, which
report 01 (Methods) puts at 0.1 to 0.2 week at the symptom level.
(Revised 2026-10-09: this paragraph formerly called 0.5 week the
value for the project's PTSD and prazosin context.)

**Intuition.** A regressor chosen without the outcome cannot borrow
the outcome's noise.

**Procedure.**

```{.r}
fit <- fit_lme(FM1, with_x(d, 0.5))            # REML
cr2t(fit)                                       # primary test
sapply(G5, function(h) cr2t(fit_lme(FM1, with_x(d, h)))['t'])   # sensitivity
```

**Validity.** Exact in the sense that matters here: no selection
occurs, so the size is that of the CR2 test at a fixed regressor,
whatever the true half-life (Section 3.3). Conditions: those of F1
itself (`docs/56`). The sensitivity grid is descriptive; it must not
be used to pick the reported value.

**Cost.** One fit per trial; five for the sensitivity grid.

**Reports.** $\hat\beta_{bm:D}$ at the stated $t_{1/2}^a$ with its
CR2 interval, and the five sensitivity $t$-statistics.

**Why it is not enough.** If the trial's profile likelihood excludes
the stated value (as it would exclude 1 week in the example trial of
Section 5), a reader can fairly ask why the analysis uses a
regressor the data contradict. A fixed value may also lose power when
it is far from the best-fitting one: with $t_{1/2}^a = 0.5$, F1's
size-adjusted power at $c_{bm} = 0.30$ was 0.496 when the true
half-life was 0 and 0.272 when it was 1 week (quoted, `docs/47`,
Section 6). How much of that fall is due to the mismatch, rather than
to the longer off-drug fading of the coupling itself, has not been
separated.

**Status.** Evaluated (`docs/47`): size 0.051 pooled, at most 0.062 in
a cell, stored setting. Extension (Section 8, refined process): 0.053
pooled over 12 cells, at most 0.062; at 0.2 week 0.050, with no power
gain over 0.5 week. The validated simple fallback.

### 4.2 Method 2: selection blind to the interaction (F1A0)

**Idea.** Choose $t_{1/2}^a$ by AIC (equivalently ML, Section 3.4) on
F1 *without* the interaction term; refit F1 with the interaction at
the chosen value by REML; test with CR2.

**Intuition.** The half-life is identified, as far as it is
identified at all, by the drug-response trajectory: how far the
outcome falls back at the off-drug visits. That information sits in
the drug main effect. Removing the interaction from the selection
model removes the $\tfrac12 T(h)^2$ bonus of Section 3.4.

**Procedure.**

```{.r}
ll0 <- sapply(G5, function(h) ll(FM0, d, h))   # ML, no interaction
h0  <- G5[which.max(ll0)]
cr2t(fit_lme(FM1, with_x(d, h0)))              # REML refit, CR2 test
```

**Validity (derived to first order, with conditions).** Two arguments,
for the two regimes of Section 3.4.

1. *Weak identification: selection on the outcome alone.* Suppose,
   as in every non-prognostic null cell, that the biomarker is
   independent of the whole outcome vector $Y$ and of the path
   assignment. Consider first a selection rule that uses $Y$ and the
   design but no biomarker term, $\hat h = g(Y)$. Conditional on $Y$,
   $\hat h$ is fixed and the biomarker values keep their marginal
   distribution. To first order the CR2 statistic at any fixed $h$ is
   a self-normalized sum,
   $$
   T(h) \approx \frac{\sum_i B^c_i\, a_i(h, Y)}
     {\big\{\sum_i (B^c_i)^2\, a_i(h, Y)^2\big\}^{1/2}},
   $$
   where $B^c_i$ is the centered biomarker and $a_i(h, Y)$ is
   participant $i$'s contribution to the exposure contrast after the
   other terms are partialled out (treating, to first order, the
   partialling of the biomarker terms $B^c_i$ and $B^c_i\eta_t$ as fixed).
   Given $Y$, the $B^c_i$ are
   exchangeable and mean zero, so $T(h)$ is approximately standard
   normal for every $h$ and every $Y$ for which no single participant
   dominates the weights. The rejection probability given $Y$ is then
   approximately $\alpha$ whatever $g(Y)$ selects, and so is the
   unconditional size (derived to first order). The argument does not
   use identification of the half-life at all; it covers the flat
   case, where F1A fails.
2. *Strong identification: orthogonality.* When the drug main effect
   identifies the half-life, $\hat h_0$ is a consistent estimator of
   the best-fitting value, and the half-life is asymptotically
   information-orthogonal to $\beta_{bm:D}$ under the null when paths
   are randomized independently of the biomarker (derived, Section
   4.6). Estimating an orthogonal nuisance parameter does not change
   the first-order distribution of the test.

F1A0 as specified keeps $B^c_i$ and $B^c_i\eta_t$ in the selection model,
so $\hat h_0$ depends on the biomarker through those terms. Their
contribution to $\ell_0(h)$ changes with $h$ only through the
partialling of $D_{bc,it}(h)$, a second-order effect (argued). A stricter
variant drops every biomarker term from the selection model and
satisfies argument 1 as stated.

**What would make it fail** (argued):

- *A prognostic biomarker.* Then $B^c_i$ is not independent of $Y$, and
  argument 1 no longer applies. F1's protection against prognostic
  leakage is approximate and depends mildly on $h$ (`docs/56`), so a
  selection that moves $h$ with the outcome could interact with the
  leakage. The illustration includes a prognostic null cell
  (Section 6).
- *Informative dropout related to the biomarker,* or allocation that
  depends on the biomarker, which break the independence and the
  orthogonality.
- *Very small $N$,* where the normal approximation to the
  self-normalized sum is poor; this affects the fixed-$h$ test too,
  which is why CR2 with Satterthwaite degrees of freedom is used.
- *Under the alternative,* the selection model omits a true term and
  can choose $h$ to absorb part of the interaction through the main
  effect. That costs power, not size.

**Cost.** Five ML fits and one REML fit with CR2 per trial, about 0.5
seconds.

**Reports.** $\hat h_0$, and $\hat\beta_{bm:D}$ with its CR2 interval
at $\hat h_0$.

**A guarded variant (F1G; post hoc).** Keep the value fixed in
advance, $t_{1/2}^a = 0.5$, unless the null-model profile rejects it at
level 0.05, $2\{\max_h \ell_0(h) - \ell_0(0.5)\} > \chi^2_{1,0.95}$, and
only then use $\hat h_0$:

```{.r}
dev05 <- 2 * (max(ll0) - ll0[G5 == 0.5])
h_use <- if (dev05 > qchisq(0.95, 1)) G5[which.max(ll0)] else 0.5
```

It answers the face-validity objection directly (the analysis never
uses a half-life that its own data reject) while selecting in only
the few trials where the default is contradicted, so whatever
selection effect F1A0 has is confined to those trials (argued). It
was defined after the illustration's results were seen.

**Status.** Illustrated (Section 6): 0.067 pooled null against F1's
0.061 on the same trials; F1G 0.063. Evaluated in the extension
(Section 8): 0.053 pooled, F1's rate on the same trials, at most 0.062
in a cell; F1G 0.053. As nominal tests neither gains power over F1;
with bootstrap calibration (F1A0B) F1A0 is the primary data-driven
candidate (Section 7.2).
`docs/56`, Section 6.1, found that on 2 × 50 null trials F1A0 picked
the same grid value as F1A in 86% and 88% of trials (quoted); here it
did in 87% and 89% (verified, Table 8).

### 4.3 Method 3: sample splitting and cross-fitting

**Idea.** Choose the half-life on one random half of the participants
and test on the other half, so that the test never sees the data that
chose its regressor [8]. Cross-fitting swaps the halves and combines
the two tests, recovering most of the lost information [9].

**Intuition.** Participants are independent. Once the half-life is a
function of half A only, it is a constant as far as half B is
concerned.

**Procedure.** Split within path, so that both halves keep the
design's balance; choose with the F1A0 criterion (any criterion is
valid here, including F1A's); test with F1 at the other half's choice.

```{.r}
ids <- unique(data.table(ptID = d$ptID, path = path_of(d)))
ids[, half := sample(rep(1:2, length.out = .N)), by = path]
dA <- d[ptID %in% ids[half == 1, ptID]]; dB <- d[ptID %in% ids[half == 2, ptID]]
choose <- function(dd) G5[which.max(sapply(G5, function(h) ll(FM0, dd, h)))]
hA <- choose(dA); hB <- choose(dB)
tB <- cr2t(fit_lme(FM1, with_x(dB, hA)))   # single split: test on B at hA
tA <- cr2t(fit_lme(FM1, with_x(dA, hB)))   # swapped
## cross-fit by averaging, with a Welch-Satterthwaite df
est <- (tA['est'] + tB['est']) / 2
se  <- sqrt(tA['se']^2 + tB['se']^2) / 2
df  <- (tA['se']^2 + tB['se']^2)^2 / (tA['se']^4 / tA['df'] + tB['se']^4 / tB['df'])
p_avg <- 2 * pt(-abs(est / se), df)
## or by the Cauchy combination of the two p-values [10]
p_cct <- 0.5 - atan(0.5 * tan((0.5 - tA['p']) * pi) + 0.5 * tan((0.5 - tB['p']) * pi)) / pi
```

**Validity.**

- *Single split: exact with respect to the selection* (derived).
  Conditional on half A, $h_A$ is fixed and half B is independent of
  half A, so the test on B is a fixed-$h$ test. Its size is the
  fixed-$h$ size of F1 at $N/2$; no further condition is needed.
- *Cross-fit average: approximate* (derived to first order). Under the
  null, $E[\hat\beta_A \mid \text{half B}] \approx 0$ because given
  half B the regressor for half A is fixed, so the two half estimates
  are approximately uncorrelated and the variance of their average is
  approximately $(\mathrm{SE}_A^2 + \mathrm{SE}_B^2)/4$. They are not
  independent, since each half's data choose the other's regressor.
- *Cauchy combination: approximate.* Valid in the upper tail under
  arbitrary dependence as $\alpha \to 0$ [10]; at $\alpha = 0.05$ it is
  close to, and slightly below, nominal for two moderately dependent
  $p$-values (argued). Bonferroni, $2\min(p_A, p_B)$, is exact and
  conservative.

**Cost.** Ten ML fits on half the data and two REML fits with CR2:
about the cost of F1A0. The cost that matters is power. A single split
tests on $N/2$, so the noncentrality falls by a factor of about
$\sqrt 2$; the cross-fit average uses all participants but each at a
regressor chosen on half the data (argued; Section 6 measures it).

**Reports.** The two choices $h_A$ and $h_B$, the cross-fit estimate
and its interval or the combined $p$-value. The result depends on the
random split; the split seed should be fixed in advance and reported,
or the analysis repeated over many splits and summarized (for
example by the median $p$-value, which needs its own calibration).

**Status.** Illustrated (Section 6): single split 0.044 pooled null,
power 0.208; cross-fit average 0.061 and 0.372; Cauchy 0.050 and
0.292.

### 4.4 Method 4: reporting the half-life

**Idea.** Whatever the test uses, report what the trial says about the
half-life: the ML estimate under the null model, the profile
likelihood over a fine grid or a continuous range, and a
profile-likelihood interval.

**Intuition.** The profile shows directly which half-lives the data
contradict, which is the face-validity question.

**Procedure.**

```{.r}
fine <- c(0, 0.02, 0.05, seq(0.1, 1, 0.05), seq(1.1, 3, 0.1), 3.5, 4, 5, 6, 8, 10, 15, 20, 50)
prof <- sapply(fine, function(h) ll(FM0, d, h))
opt  <- optimize(function(lh) ll(FM0, d, exp(lh)), log(c(0.01, 100)), maximum = TRUE)
lmax <- max(prof, opt$objective)
ci   <- range(fine[2 * (lmax - prof) <= qchisq(0.95, 1)])
```

**Validity.** The interval is the set of half-lives not rejected by a
likelihood-ratio test at level 0.05 with a $\chi^2_1$ reference. That
reference is not exact here (derived from standard theory): the
parameter has a boundary at 0, where the binary indicator lies, and a
second one at $\infty$, where $D_{bc,it}$ becomes an ever-exposed
indicator [11]; and when the profile is flat, the half-life is weakly
identified and the $\chi^2$ approximation fails. The interval should
be read as a description of the likelihood, not as an interval with
known coverage. Its coverage of the best-fitting working value has
not been evaluated. It is an interval for the working parameter
(Section 3.2), not for the pharmacological half-life.

**Cost.** On a 52-point grid, two ML fits per point, plus a REML refit
and CR2 for the $t$-profile: 11 seconds per trial (verified, Section
5); once per trial, not per resample.

**Reports.** $\hat t_{1/2}^a$, the interval, and the deviance at 0 and
at long half-lives. The honest statement for a Hybrid trial is that
its data usually separate a short half-life from one of a week or
more, but not 0 from 0.25 or 0.5 week (Sections 3.5 and 5).

**Status.** Illustrated on one trial (Section 5); coverage not
evaluated. Recommended for reporting in every case.

### 4.5 Method 5: parametric bootstrap of the whole procedure (F1AB)

**Idea.** Calibrate the select-then-test statistic by simulating it
under a fitted null model that reproduces the trial's design, its
estimated drug trajectory and its estimated covariance, rerunning the
selection in every resample.

**Intuition.** F1A's flaw is that its reference distribution ignores
the selection. A bootstrap that repeats the selection includes it.

**Procedure.**

1. Fit the null model FM0 with its own selected half-life $\hat h_0$
   (ML choice, REML refit), giving $\hat\xi$, $\hat\sigma_u$,
   $\hat\sigma_\varepsilon$ and the CAR(1) parameter $\hat\rho_a$.
2. For $m = 1, \ldots, M$: simulate
   $Y^{*}_i = W_i(\hat h_0)\hat\xi + u^{*}_i + \varepsilon^{*}_i$ for
   every participant, with the observed visits, biomarker values and
   paths; rerun the whole procedure (selection and test) on $Y^*$ and
   record $|T^{*}_m|$.
3. Report $p^{*} = (1 + \#\{m : |T^{*}_m| \ge |T_{\text{obs}}|\})/(M + 1)$.

```{.r}
fit0 <- fit_lme(FM0, with_x(d, h0))                    # REML null fit at h0
sim_null <- function(fit, d, h0) {
  dd  <- with_x(d, h0)
  mu  <- drop(model.matrix(FM0[-2], dd) %*% nlme::fixef(fit))
  tau <- as.numeric(nlme::VarCorr(fit)['(Intercept)', 'StdDev'])
  phi <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  y <- mu
  for (idx in split(seq_len(nrow(dd)), dd$ptID)) {
    tt <- dd$t[idx]
    S  <- tau^2 + fit$sigma^2 * phi^abs(outer(tt, tt, '-'))
    y[idx] <- mu[idx] + drop(crossprod(chol(S), rnorm(length(idx))))
  }
  ds <- copy(d); ds[, Sx := y]; ds
}
tstar <- replicate(M, abs(stat_F1A(sim_null(fit0, d, h0))))   # stat_F1A: select, refit, CR2 t
p_boot <- (1 + sum(tstar >= abs(t_obs))) / (M + 1)
```

**Validity.** By construction, up to the fit of the null model
(argued, with the usual parametric-bootstrap conditions). The
bootstrap distribution of the selected statistic is computed under a
null model with the estimated parameters; if that model describes the
null data well, $p^{*}$ is approximately uniform. It can fail when the
null model is misspecified in ways that matter for the statistic's
distribution: a residual covariance other than random intercept plus
CAR(1), on-drug means that a constant $\beta_D$ fits poorly, or a
prognostic biomarker whose association with $Y$ the null model
reproduces only through $B^c_i$ and $B^c_i\eta_t$. Because the resamples
keep the observed biomarker values, paths and visits, the bootstrap
conditions on the design. The bootstrap is calibrated against the
selection under the fitted null model; it does not repair a
misspecified analysis model.

**Choice of $M$.** For a test at 0.05, $M = 199$ or $M = 999$ (so that
$(M + 1)\alpha$ is an integer). The Monte Carlo standard error of
$p^*$ near 0.05 is about $\sqrt{0.05 \times 0.95 / M}$, 0.015 at
$M = 199$ and 0.007 at $M = 999$ (derived).

**Cost.** Each resample costs one run of the procedure: about 1.1
seconds for the grid with and without the interaction. For one trial,
$M = 199$ took 128 seconds on 7 workers (verified, Section 5); $M =
999$ would take about 11 minutes. In a simulation study, a full
bootstrap per replicate multiplies the cost by $M$. The warp-speed
method of Giacomini, Politis and White [5] draws one resample per Monte
Carlo replicate and uses the pooled resampled statistics of a cell as
the bootstrap distribution, which estimates the rejection rate of the
bootstrap test at about twice the cost of the plain procedure. Its
justification requires the bootstrap distribution to vary little
across replicates of a cell [5] (assumed here); its rejection rate
has a Monte Carlo standard error somewhat larger than the binomial
one reported below.

**Reports.** $p^{*}$ for the chosen statistic (F1A, F1A0 or max-$|t|$),
the bootstrap critical value, and $\hat h_0$.

**Status.** Illustrated by the warp-speed method (Section 6: F1AB
0.062 pooled null, F1's own rate on the same trials being 0.061; power
0.364) and on one trial with $M = 199$ (Section 5). Evaluated by the
warp-speed method in the extension (Section 8): F1A0B 0.050 and F1AB
0.050 pooled, mean size-adjusted power 0.395 and 0.400, maximin
shortfall 0.076 and 0.078. The prognostic-baseline cell, where the B
arms reached 0.067 to 0.077, does not show the null-model failure
described under Validity: there the null model reproduces the
association exactly and the excess is Monte Carlo error (Section 8.5).

### 4.6 Method 6: joint estimation of the half-life as a continuous parameter

**Idea.** Treat $\omega = \log t_{1/2}^a$ as a parameter of a
nonlinear mixed model in which $D_{bc,it}(e^\omega)$ enters both the
main effect and the interaction, estimate it with the other
parameters, and test $\beta_{bm:D}$ with a standard error that
includes the uncertainty in $\hat\omega$.

**Intuition.** This is the textbook answer: estimate every unknown,
and let the information matrix decide how much the estimation of
$\omega$ costs the test.

**Procedure** (sketch; not implemented).

- *Estimation:* either `nlme::nlme()` with the mean function
  $\beta_0 + \beta_B B^c_i + \beta_T w_t + (\beta_D + \beta_{bm:D}
  B^c_i)\,D_{bc,it}(e^\omega) + \beta_\eta\eta_t + \beta_{bm:\eta}
  B^c_i\eta_t$, a random intercept and `corCAR1`, or profile
  $\hat\omega = \arg\max_\omega \ell_1(e^\omega)$ by `optimize()` and
  fit the linear model at $\hat\omega$.
- *Inference:* a Wald test with the sandwich variance of the stacked
  estimating equations for $(\beta, \omega)$, in which the derivative
  $\partial D_{bc,it}/\partial\omega = D_{bc,it}\, \ln 2 \;
  t_{sd}\, e^{-\omega}$ (off drug after exposure; 0 otherwise)
  enters the bread, or a score test of $\beta_{bm:D} = 0$ at the null
  fit $(\hat\omega_0, \hat\xi_0)$ with the efficient score
  $U_{\beta} - I_{\beta\omega} I_{\omega\omega}^{-1} U_\omega$.

**Validity and the orthogonality result.** Under the null, the
cross-information between $\beta_{bm:D}$ and $\omega$ is

$$
I_{\beta\omega} = \beta_D \sum_i B^c_i\,
  D_{bc,i}(e^\omega)^{\top} V_i^{-1}\,
  \frac{\partial D_{bc,i}}{\partial \omega},
$$

plus terms with mean-zero residuals (derived). The vectors
$D_{bc,i}$ and $\partial D_{bc,i}/\partial\omega$ are the same for all
participants on the same path, so the sum is
$\beta_D \sum_q c_q \sum_{i \in q} B^c_i$ over paths $q$, and is
$O_p(N^{1/2})$ when paths are assigned independently of the
biomarker, against $O(N)$ for $I_{\beta\beta}$ and $I_{\omega\omega}$
(derived). The half-life is therefore asymptotically orthogonal to the
interaction under the null: its estimation does not change the
first-order null distribution of the score test, and F1A0 is that
score test with $\omega$ restricted to a grid. The orthogonality needs
three conditions: randomization independent of the biomarker,
$\beta_D$ and $I_{\omega\omega}$ bounded away from zero (a drug effect
that identifies the half-life), and $\omega$ in the interior. The last
two fail in practice: $I_{\omega\omega}$ is small in a Hybrid trial
(Section 3.5), the ML estimate is often on the boundary at 0 (it is in
the example trial, Section 5), and then the asymptotic argument says
nothing [11]. The Wald test at the joint estimate also inherits F1A's
problem when $\omega$ is weakly identified, since the joint estimate
of $\omega$ then responds to the interaction's noise.

**Small samples.** No CR2 analog exists for a nonlinear mean
parameter; a CR2-type correction would have to be developed and
checked (argued).

**Cost.** One `nlme` fit per trial, typically several times an `lme`
fit, with convergence failures near the boundaries; or a profile over
$\omega$ (about 25 one-dimensional evaluations) plus a sandwich.

**Reports.** $\hat t_{1/2}^a$ with a standard error, and the
interaction test adjusted for its estimation.

**Status.** Proposed. Its main contribution here is the orthogonality
argument, which supports F1A0.

### 4.7 Method 7: the maximum statistic over the grid

**Idea.** Compute $T(h)$ at every $h \in \mathcal T_a$ and test with
$\max_h |T(h)|$, calibrated by its joint null distribution.

**Intuition.** Instead of guessing the half-life, ask whether there is
an interaction at *some* half-life on the grid, and pay for the
multiplicity explicitly. This is the problem of a nuisance parameter
present only under the alternative studied by Davies [6, 7] and
Hansen [12]: if $\beta_D$ were 0, $h$ would appear only in the
interaction term and be unidentified under the null.

**Procedure.**

```{.r}
Th <- sapply(G5, function(h) cr2t(fit_lme(FM1, with_x(d, h)))['t'])
Tmax <- max(abs(Th))
p_bonf <- min(1, length(G5) * min(2 * pt(-abs(Th), df)))   # Bonferroni: exact, conservative
## calibrated: compare Tmax with the max over the grid in null-model resamples (Section 4.5)
## or with the max of a multivariate normal with the correlation of the T(h)
```

The correlation of the $T(h)$ can be estimated from the stacked
participant-level influence contributions of the five fits (a
cluster-robust covariance of the five interaction estimates), and the
critical value read from the maximum of a multivariate normal or $t$
with that correlation (`mvtnorm`), or taken from the null-model
bootstrap.

**Validity.** Bonferroni is exact with respect to the multiplicity if
each $T(h)$ is a valid fixed-$h$ test, and conservative because the
$T(h)$ are strongly correlated. The bootstrap calibration is valid up
to the fit of the null model (Section 4.5). The multivariate normal
calibration is approximate (first order). Davies' bounds [6, 7] cover
the continuous-$h$ version.

**What it tests.** The same null hypothesis as F1 (no interaction,
for any exposure regressor), but against a different alternative: an
interaction whose best-fitting half-life is unknown. The reported
quantity is a $p$-value and the half-life at which $|T|$ is largest;
there is no single estimand, which weakens its face validity for a
reader who wants an effect size.

**Cost.** Five REML fits with CR2 (about 1 second); bootstrap
calibration adds $M$ times that per trial.

**Status.** Illustrated (Section 6): Bonferroni 0.024 pooled null,
power 0.220; warp-speed bootstrap 0.048 and 0.384. Evaluated in the
extension (Section 8): Bonferroni 0.018 (too conservative), the
multivariate normal reference 0.072 (liberal), the multivariate $t$
0.053 (flagged at 0.077 in one prognostic cell), the warp-speed
bootstrap 0.051 with mean size-adjusted power 0.388.

### 4.8 Method 8: marginalizing over the half-life

**Idea.** Put a prior on $t_{1/2}^a$ informed by pharmacology, and base
the decision on the posterior of $\beta_{bm:D}$ integrated over the
half-life.

**Intuition.** The analyst does not know the half-life; averaging
over the plausible values with weights given by prior and data is the
coherent way to carry that uncertainty.

**Procedure.** Full Bayesian: a Stan model with the F1 mean function,
$D_{bc,it}(e^\omega)$ computed inside the model, a random intercept and
CAR(1) residuals (not available directly in `brms`, which has no
continuous-time AR(1) for irregular visits, so a custom Stan program
is needed), a prior on $\omega$ (for example normal on the log scale
centered at $\log 0.5$ with a standard deviation of 1), and weakly
informative priors elsewhere. The decision rule might be: reject if
the 95% central posterior interval of $\beta_{bm:D}$ excludes 0. A
cheaper frequentist stand-in is model averaging over the grid with
smoothed AIC weights (Buckland, Burnham and Augustin [13]):

```{.r}
w   <- exp(ll1 - max(ll1)); w <- w / sum(w)          # ll1: ML of F1 on the grid
est <- sum(w * est_h)                                 # est_h, se_h: REML + CR2 at each h
se  <- sum(w * sqrt(se_h^2 + (est_h - est)^2))
```

**Validity.** Unknown. A posterior-interval decision rule has no
guaranteed frequentist size; it must be calibrated by simulation, and
the calibration depends on the prior. The weights use $\ell_1$, which
includes the interaction, so averaging inherits a damped version of
F1A's preference for half-lives at which $|T|$ is large (derived from
Section 3.4). Bayesian model averaging is reviewed by Hoeting et al.
[14].

**Cost.** An MCMC fit per trial, perhaps 1 to 5 minutes with a custom
Stan program (assumed, not measured), which makes a 28,000-trial
calibration expensive. The AIC-weight stand-in costs the grid fits
only.

**Reports.** The posterior of $\beta_{bm:D}$ and of $t_{1/2}^a$; the
stand-in reports a model-averaged estimate and interval.

**Status.** Proposed (Bayesian); the AIC-weight stand-in is
illustrated (Section 6: 0.053 pooled null with $\ell_1$ weights and
0.049 with $\ell_0$ weights; power 0.388 and 0.356). Evaluated in the
extension (Section 8): 0.0455 with $\ell_1$ weights and 0.040 with
$\ell_0$ weights, both below the band; with bootstrap calibration, MAB
0.048 and MA0B 0.049, mean size-adjusted power 0.401 and 0.399.

### 4.9 Permutation tests

**Naive permutation of the biomarker is not valid for this null.**
Permuting $B^c_i$ across participants generates the distribution of the
statistic under the *strong* null that the biomarker is unrelated to
the outcome in every respect. The hypothesis of interest is weaker:
$\beta_{bm:D} = 0$ while $B^c_i$ may be prognostic, with a main effect
or an association that varies with expectancy. When the biomarker is
prognostic, permutation destroys that association in the resamples
but not in the observed data, so the observed statistic is compared
with the wrong distribution (derived). For a studentized statistic
with the prognostic terms in the model the damage may be small, but
nothing guarantees it, and the test then answers a question nobody
asked. When the biomarker is not prognostic the strong and weak nulls
coincide and permutation within path is exact, but that is the case in
which no help is needed.

**What a valid version would require.** A residual permutation in the
manner of Freedman and Lane [15]: fit the null model FM0 (with its
selected half-life), permute the participant-level residual vectors
(all visits of a participant together, to keep the within-participant
correlation) among participants *of the same path* (whose design rows
are identical), add them to the fitted null values, rerun the
selection and the test, and refer the observed statistic to the
permutation distribution. This requires the residual vectors to be
exchangeable within path, which holds only approximately (estimated
residuals, dropout, phase-dependent variances), and it is valid only
approximately, at about the cost of the parametric bootstrap.

**Status.** Proposed with conditions; not recommended over the
parametric bootstrap, which needs fewer exchangeability assumptions
and handles unequal path sizes and dropout directly.

## 5. Worked example 1: one null trial

### 5.1 Set-up

One trial was drawn from the base null cell of the stress test with
the *refined* process of the driver (anchored carryover, randomized
allocation; `--carryover anchored --allocation random`): $N = 70$,
path sizes 14, 23, 13 and 20, AR(1) at $\rho = 0.7$, true
$t_{1/2} = 0.5$ week, exponential decay, $c_{bm} = 0$, no dropout and
no prognostic pathway (inspected, `57-example-trial.R`). The seed
(5701) was fixed before the trial was seen, and the trial was not
selected for any property.

### 5.2 The profile and the test across the half-life

| $h$ (weeks) | $\ell_1(h)$ | $\ell_0(h)$ | $\mathrm{LR}(h)$ | $\hat\beta_{bm:D}$ | CR2 SE | $T(h)$ | df | $p$ |
|------|------|------|------|------|------|------|------|------|
| 0 | $-2591.658$ | $-2591.997$ | 0.678 | $-0.0602$ | 0.0752 | $-0.801$ | 26.0 | 0.430 |
| 0.25 | $-2591.899$ | $-2592.201$ | 0.605 | $-0.0589$ | 0.0781 | $-0.753$ | 25.9 | 0.458 |
| 0.5 | $-2593.667$ | $-2593.852$ | 0.371 | $-0.0510$ | 0.0874 | $-0.584$ | 26.0 | 0.564 |
| 1 | $-2597.846$ | $-2597.903$ | 0.113 | $-0.0330$ | 0.1022 | $-0.323$ | 26.5 | 0.749 |
| 2 | $-2601.938$ | $-2601.939$ | 0.002 | $-0.0048$ | 0.1202 | $-0.040$ | 27.3 | 0.968 |

Table: Example 1: ML log-likelihoods with and without the interaction, the likelihood-ratio statistic, and the REML + CR2 test of the interaction at each grid half-life (verified)

![Example 1, one null Hybrid trial ($N = 70$, true $t_{1/2} = 0.5$ week, refined process). A: ML deviance $2\{\ell_{\max} - \ell(h)\}$ of F1 with the interaction (solid) and without it (dashed) over a fine grid of analysis half-lives; the horizontal line is $\chi^2_{1,0.95} = 3.84$, and points mark the grid $\mathcal T_a$. B: the CR2 $t$-statistic of the interaction across the same half-lives, with $\pm 1.96$ marked.](figures/57-fig1-example-trial.png)

Four observations (verified):

- **The profile is flat where it matters and steep elsewhere.** The
  null-model deviance is 0 at $h = 0$, 0.41 at 0.25 week and 3.71 at
  0.5 week, then 11.8 at 1 week and 19.9 at 2 weeks, rising to 30.8 as
  $h \to \infty$. The trial cannot tell 0 from 0.5 week but excludes a
  week or more.
- **The ML estimate is on the boundary.** The continuous maximum over
  $[0.01, 100]$ weeks is at 0.024 week, the lower end of the range in
  practice, and $\ell_0(0)$ is at least as large: the binary indicator
  fits best, although the true half-life is 0.5 week. The
  profile-likelihood interval is $[0, 0.50]$ weeks for the null model
  and $[0, 0.45]$ for F1, so the true value lies at the upper edge of
  one and just outside the other; with a working parameter and a
  boundary estimate, neither interval has a known coverage (Section
  4.4).
- **The two profiles are almost parallel.** $\mathrm{LR}(h)$ is below
  0.7 at every $h$, so in this trial F1A and F1A0 make the same choice.
- **The interaction statistic changes smoothly and modestly with
  $h$:** from $-0.80$ at 0 to $-0.04$ at 2 weeks and $+1.1$ near 20
  weeks.

### 5.3 What each method reports

| Method | Half-life used (weeks) | $\hat\beta_{bm:D}$ | Statistic | $p$ |
|------|------|------|------|------|
| 1. F1, fixed 0.5 | 0.5 | $-0.051$ | $t = -0.58$ | 0.564 |
| F1A (AIC with the interaction) | 0 | $-0.060$ | $t = -0.80$ | 0.430 |
| 2. F1A0 (AIC without it) | 0 | $-0.060$ | $t = -0.80$ | 0.430 |
| 3. single split (choose on A, test on B) | 0 | 0.048 | $t = 0.42$ (df 8.8) | 0.684 |
| 3. cross-fit, average | 0 and 0.25 | $-0.022$ | $t = -0.28$ (df 21.7) | 0.781 |
| 3. cross-fit, Cauchy | 0 and 0.25 | | | 0.544 |
| 4. profile interval | ML at 0; interval $[0, 0.50]$ | | | |
| 5. F1AB, bootstrap of F1A ($M = 199$) | 0 | | $\lvert t\rvert = 0.80$; critical value 1.93 | 0.435 |
| 5. bootstrap of F1A0 | 0 | | critical value 1.92 | 0.430 |
| 5. bootstrap of F1 at 0.5 | 0.5 | | critical value 1.92 | 0.550 |
| 7. max-$\lvert t\rvert$, Bonferroni | 0 | | $\max\lvert t\rvert = 0.80$ | 1.000 |
| 7. max-$\lvert t\rvert$, bootstrap ($M = 199$) | 0 | | critical value 2.19 | 0.560 |
| 8. AIC-weight averaging | weighted mean 0.14 | $-0.059$ | $t = -0.76$ | 0.452 |

Table: Example 1: what each procedure reports for the same null trial (verified)

In this trial every method agrees that there is no evidence of an
interaction, as it should. Three details are worth noting (verified):

- **The bootstrap critical values.** For a single fixed-$h$ test the
  bootstrap gives 1.92, close to the $t_{26}$ value 2.06 and below it
  (199 resamples, so the quantile has a standard error of roughly
  0.13, derived from the binomial error of a 95% quantile).
  For F1A it gives 1.93, barely more: in this trial the selection adds
  little, because $\ell_0$ is steep enough to pin the choice near 0.
  For the maximum over the grid it gives 2.19, the price of
  multiplicity across five strongly correlated statistics.
- **The single split** tests on 35 participants with only 8.8
  Satterthwaite degrees of freedom; its estimate has the opposite sign
  and a standard error 1.5 times F1's. The cross-fit average recovers
  a standard error close to F1A0's (0.076 against 0.075).
- **Timing.** The fine profile (52 half-lives, three fits and one CR2
  each) took 11 seconds; the five-point grid with the split 1.5
  seconds; the bootstrap with $M = 199$ 128 seconds on 7 workers; the
  whole example 146 seconds.

## 6. Worked example 2: a small Monte Carlo illustration

### 6.1 Set-up

(2026-10-10: the stress-test extension of Section 8 supersedes this
illustration as evidence on size and power. The illustration is kept
for its paired comparison with `docs/47` and for the procedures the
extension did not run, cross-fitting and the single split.)

This is an illustration, not a replacement for the stress-test
extension of Section 7.3. It uses the trials of the stored stress test
exactly: the driver's data-generating functions, cells and seeds
(chunks of 25 replicates, seed $20261007 + 1000 \times \text{cell} +
\text{chunk}$) with `--carryover recursive --allocation fixed`, the
setting of the stored results, so that replicate $r$ here is replicate
$r$ of `docs/47` (inspected, `57-null-illustration.R`). The pairing
was checked: F1 and F1A reproduce the stored estimates exactly in
every replicate (maximum absolute difference 0 in all 1,500). The split and the single bootstrap resample
per replicate use a separate random stream, so they do not disturb the
trials.

Four cells were run: the base null cell (cell 2, $n_{\text{sim}} = 500$), the null
cell with true $t_{1/2} = 1$ week where F1A failed (cell 20, $n_{\text{sim}} =
500$), the null cell with a constant prognostic natural-course pathway
$c_{bm,TV} = 0.3$ (cell 25, $n_{\text{sim}} = 250$), and the base alternative cell
$c_{bm} = 0.30$ (cell 6, $n_{\text{sim}} = 250$). The replicate counts were cut from
1,000 and 500 to keep the run within the time budget (Section 6.4).
Each replicate fits F1 at every grid value, the selection models, the
split and cross-fit, and one null-model bootstrap resample on which the
F1, F1A, F1A0 and max-$|t|$ procedures are rerun (warp-speed method
[5]).

### 6.2 Results

| Method | Base null | Null, $t_{1/2} = 1$ | Null, prognostic TV | Pooled null (base and $t_{1/2} = 1$) | Power, $c_{bm} = 0.30$ |
|------|------|------|------|------|------|
| $n_{\text{sim}}$ | 500 | 500 | 250 | 1,000 | 250 |
| F1, fixed 0.5 | 0.064 (0.011) | 0.058 (0.010) | 0.032 (0.011) | 0.061 (0.008) | 0.400 (0.031) |
| F1A | 0.070 (0.011) | 0.078 (0.012) | 0.040 (0.012) | 0.074 (0.008) | 0.428 (0.031) |
| F1A0 | 0.064 (0.011) | 0.070 (0.011) | 0.040 (0.012) | 0.067 (0.008) | 0.396 (0.031) |
| F1G, guarded F1A0 (post hoc) | 0.064 (0.011) | 0.062 (0.011) | 0.032 (0.011) | 0.063 (0.008) | 0.400 (0.031) |
| AIC-weight averaging, $\ell_1$ weights | 0.054 (0.010) | 0.052 (0.010) | 0.032 (0.011) | 0.053 (0.007) | 0.388 (0.031) |
| AIC-weight averaging, $\ell_0$ weights (post hoc) | 0.050 (0.010) | 0.048 (0.010) | 0.032 (0.011) | 0.049 (0.007) | 0.356 (0.030) |
| single split | 0.046 (0.009) | 0.042 (0.009) | 0.044 (0.013) | 0.044 (0.006) | 0.208 (0.026) |
| cross-fit, average | 0.056 (0.010) | 0.066 (0.011) | 0.032 (0.011) | 0.061 (0.008) | 0.372 (0.031) |
| cross-fit, Cauchy | 0.050 (0.010) | 0.050 (0.010) | 0.032 (0.011) | 0.050 (0.007) | 0.292 (0.029) |
| max-$\lvert t\rvert$, Bonferroni | 0.024 (0.007) | 0.024 (0.007) | 0.016 (0.008) | 0.024 (0.005) | 0.220 (0.026) |
| F1, warp-speed bootstrap | 0.062 (0.011) | 0.060 (0.011) | 0.028 (0.010) | 0.061 | 0.388 (0.031) |
| F1AB, warp-speed bootstrap of F1A | 0.066 (0.011) | 0.058 (0.010) | 0.036 (0.012) | 0.062 | 0.364 (0.030) |
| F1A0, warp-speed bootstrap | 0.070 (0.011) | 0.058 (0.010) | 0.032 (0.011) | 0.064 | 0.348 (0.030) |
| max-$\lvert t\rvert$, warp-speed bootstrap | 0.052 (0.010) | 0.044 (0.009) | 0.040 (0.012) | 0.048 | 0.384 (0.031) |

Table: Example 2: rejection rates at nominal 0.05 (binomial MCSE in parentheses) on stored stress-test trials, stored setting (verified; the warp-speed MCSEs understate the error, and their pooled column is the mean of the two cell rates)

| Procedure | Base cell | $t_{1/2} = 1$ | Prognostic TV | $c_{bm} = 0.30$ |
|------|------|------|------|------|
| F1, bootstrap | 2.09 | 2.03 | 2.10 | 2.13 |
| F1AB | 2.11 | 2.24 | 2.13 | 2.21 |
| F1A0, bootstrap | 2.01 | 2.15 | 2.13 | 2.19 |
| max-$\lvert t\rvert$, bootstrap | 2.39 | 2.43 | 2.27 | 2.30 |

Table: Example 2: warp-speed bootstrap critical values for $\lvert T\rvert$, the 95% quantile of the pooled resampled statistics of each cell (verified)

![Example 2: rejection rates at nominal 0.05 with 2-MCSE bars, by procedure, in the four stress-test cells (stored setting). Circles: nominal reference distribution; triangles: warp-speed bootstrap calibration. The vertical line marks 0.05 in the null cells; the alternative cell (top right) shows power.](figures/57-fig2-null-illustration.png)

The post hoc rows were defined after the other results were seen and
are computed from the same grid fits (`57-posthoc-rules.R`). F1G keeps
$t_{1/2}^a = 0.5$ week unless the null-model profile rejects it
($2\{\max_h \ell_0(h) - \ell_0(0.5)\} > 3.84$ on the grid), and only
then uses $\hat h_0$; it switched in 1.0%, 15.0%, 2.0% and 3.6% of the
trials of the four cells (verified).

### 6.3 Reading

All rates below are from Table 6 and the paired counts of
`57-mc-discordance.csv` (verified); with 1,000 pooled null trials the
MCSE of a rate near 0.05 is 0.007, so only differences of about 0.015
or more are resolved.

- **The pairing holds.** F1 and F1A reproduce the stored estimates
  and decisions in every replicate (maximum absolute difference 0), so
  every comparison is paired with `docs/47`. In these first 500 trials
  of the base cell F1 itself rejects 0.064 (0.062 over all 1,000
  stored trials): the excess there is F1's, not the selection's.
- **F1A0 removes about half of F1A's excess, not all of it.** Pooled,
  F1A rejects 0.074, F1A0 0.067 and F1 0.061 on the same trials. In
  the $t_{1/2} = 1$ cell F1A0 rejected 10 trials that F1 did not, and
  F1 4 that F1A0 did not; F1A 12 and 2. F1A0's choices differ from
  F1A's in only 11 to 13% of null trials (Table 8), which is where the
  difference comes from. The first-order argument of Section 4.2 is
  therefore not borne out cleanly at $N = 70$, though the excess over
  F1 (0.006, paired) is within the Monte Carlo error.
- **The guard costs nothing measurable.** F1G is F1 except in the few
  trials whose data contradict 0.5 week; it rejected 2 null trials
  that F1 did not, out of 1,250, and had the same power (0.400).
- **Averaging over the grid holds size.** With $\ell_1$ weights, 0.053
  pooled; with $\ell_0$ weights, 0.049. Both are slightly more
  conservative than F1 on the same trials (in the base cell F1
  rejected 6 and 7 trials that they did not, against 1 and 0 the other
  way), because the averaged
  standard error includes the spread of the estimates across the grid.
  Their power cost against F1 is 0.012 and 0.044, within one and a
  half MCSEs.
- **Splitting is valid and expensive.** The single split holds size
  (0.044) and loses half the power (0.208 against 0.400). The
  cross-fit average recovers most of the power (0.372) but inherits a
  little of the selection effect (0.066 in the $t_{1/2} = 1$ cell);
  the Cauchy combination is at 0.050 with power 0.292.
- **The bootstrap calibrates the selection but not F1's own excess.**
  The warp-speed critical value of F1A exceeds that of F1 by 0.02 to
  0.21, and F1AB's pooled rate (0.062) is F1's (0.061). The bootstrap
  of F1 at a fixed half-life does not bring F1 to 0.05 in the base
  cell (0.062), which is consistent with a null model (random intercept
  plus CAR(1), constant drug effect) that is lighter-tailed than the
  data-generating process (argued). The max-$|t|$ statistic with
  bootstrap calibration is at 0.048 pooled with power 0.384; its
  Bonferroni version is conservative (0.024) and loses half the power.
- **The prognostic cell** ($c_{bm,TV} = 0.3$) is below 0.05 for every
  procedure (0.016 to 0.044, $n_{\text{sim}} = 250$); no selection procedure
  reintroduced the prognostic leakage that F1 removes. This cell has
  the true half-life 0.5, so it does not test a prognostic biomarker
  combined with a wrong half-life.
- **Power is that of a cell where 0.5 week is right.** In the
  alternative cell the fixed value matches the truth, so no
  data-chosen procedure can gain; F1A's 0.428 is inflated by its size.
  The cells where a data-chosen half-life could gain power ($t_{1/2} =
  0$ or 1 week under the alternative) were not run.

| Cell | Procedure | 0 | 0.25 | 0.5 | 1 | 2 |
|------|------|------|------|------|------|------|
| base null | F1A | 0.234 | 0.234 | 0.394 | 0.110 | 0.028 |
| base null | F1A0 | 0.218 | 0.242 | 0.416 | 0.096 | 0.028 |
| base null | split (half A) | 0.282 | 0.180 | 0.346 | 0.126 | 0.066 |
| $t_{1/2} = 1$ null | F1A | 0.008 | 0.036 | 0.192 | 0.416 | 0.348 |
| $t_{1/2} = 1$ null | F1A0 | 0.004 | 0.028 | 0.218 | 0.386 | 0.364 |
| $t_{1/2} = 1$ null | split (half A) | 0.070 | 0.054 | 0.202 | 0.294 | 0.380 |

Table: Example 2: share of trials choosing each analysis half-life (weeks); F1A and F1A0 agreed in 86.6% and 89.4% of the trials of these two cells (verified)

### 6.4 Run time

The 1,500 replicates took 10,590 worker-seconds, 7.1 worker-seconds
per replicate (one grid of 15 fits and 5 CR2 computations on the
data, the split, and the same grid on one bootstrap resample),
measured on 7 workers on the machine of `docs/47`; the final run took
23 minutes of wall time, after a first launch (with twice the
replicates) that was stopped to stay within the time budget and whose
completed chunks were reused (verified, `57-mc-timing.csv`).
Example 1 took 146 seconds. The whole illustration therefore used
about 30 minutes of wall time on 7 workers.

## 7. Recommendation

### 7.1 Decision table

(Revised 2026-10-10: the table now leads with the stress-test
extension of Section 8 and keeps the illustration's figures for the
procedures the extension did not run.) "Size" is the extension's
pooled null rate over the 12 non-prognostic null cells (band 0.046 to
0.054), and "power" its mean size-adjusted power over the 26
non-prognostic alternative cells, with the value at true $t_{1/2} = 2$
weeks and $c_{bm} = 0.45$ after the slash; "max." is the maximin
shortfall. Figures marked "ill." are from the illustration of Section
6 (stored setting). "Face validity" asks whether the analysis avoids a
half-life the trial's data contradict and reports what the data say
about it.

| Method | Validity | Cost per trial | Size | Power | Face validity | Status |
|------|------|------|------|------|------|------|
| 1. F1, fixed 0.5 (stress-test base setting) | exact with respect to selection | 1 fit | 0.053 | 0.388 / 0.352; max. 0.148 | weak | evaluated; validated simple fallback |
| 1. F1, fixed 0.2 (pharmacological) | exact with respect to selection | 1 fit | 0.050 | 0.376 / 0.348; max. 0.152 | weak | evaluated; no gain over 0.5 |
| F1A (AIC with the interaction) | none: post-selection | 6 fits | 0.060 (fails) | not ranked | good | evaluated; not to be used |
| 2. F1A0 | approximate, first order | 6 fits | 0.053 | 0.381 / 0.428; max. 0.126 | good | evaluated; holds size, no power gain as a nominal test |
| F1G, guarded F1A0 | approximate; equals F1 unless 0.5 is rejected | 6 fits | 0.053 | 0.381 / 0.376; max. 0.142 | good | evaluated, post hoc; no gain over F1 |
| 3. single split | exact with respect to selection | 12 half-size fits | ill. 0.044 | ill. 0.208 | good | illustrated; not run in the extension |
| 3. cross-fit, average | approximate | 12 half-size fits | ill. 0.061 | ill. 0.372 | good | illustrated; not run in the extension |
| 3. cross-fit, Cauchy | approximate | 12 half-size fits | ill. 0.050 | ill. 0.292 | good | illustrated; not run in the extension |
| 4. profile interval | descriptive; coverage unknown | about 150 fits, 11 s | not a test | not a test | essential | illustrated, one trial; report in every case |
| 5. F1A0B, bootstrap of F1A0 | by construction, up to the null model | $M$ times F1A0; about 17 worker-minutes at $M = 999$ | 0.050 | 0.395 / 0.500; max. 0.076 | good | evaluated (warp-speed); **primary data-driven candidate** |
| 5. F1AB, bootstrap of F1A | by construction, up to the null model | $M$ times F1A | 0.050 | 0.400 / 0.456; max. 0.078 | good | evaluated (warp-speed); alternative |
| 5. F1B, bootstrap of fixed F1 | by construction, up to the null model | $M$ times F1 | 0.050 | 0.387 / 0.374; max. 0.126 | weak | evaluated (warp-speed); reference |
| 6. joint estimation | approximate; fails near boundaries | 1 `nlme` fit and a sandwich | not run | not run | good | proposed |
| 7. max-$\lvert t\rvert$, bootstrap (MaxTB) | by construction, up to the null model | $M$ times 5 fits | 0.051 | 0.388 / 0.420; max. 0.086 | moderate: no single half-life | evaluated (warp-speed) |
| 7. max-$\lvert t\rvert$, Bonferroni | exact, conservative | 5 fits | 0.018 (fails) | not ranked | moderate | evaluated; too conservative |
| 7. max-$\lvert t\rvert$, multivariate normal | asymptotic | 5 fits and influence vectors | 0.072 (fails) | not ranked | moderate | evaluated; liberal |
| 7. max-$\lvert t\rvert$, multivariate $t$ | approximate | 5 fits and influence vectors | 0.053 | 0.375 / 0.408; max. 0.160 | moderate | evaluated; prognostic flag (0.077) |
| 8. AIC-weight averaging, $\ell_1$ (MA) | unknown | 15 fits | 0.0455 (fails, conservative) | not ranked | good | evaluated; slightly conservative |
| 8. MA with bootstrap (MAB) | by construction, up to the null model | $M$ times 15 fits | 0.048 | 0.401 / 0.454; max. 0.088 | good | evaluated (warp-speed) |
| 8. MA0B, $\ell_0$ weights with bootstrap | by construction, up to the null model | $M$ times 15 fits | 0.049 | 0.399 / 0.440; max. 0.062 | good | evaluated (warp-speed), post hoc |
| 8. Bayesian marginalization | unknown; calibration needed | MCMC, minutes | not run | not run | good | proposed |
| Permutation of $B^c_i$ | invalid for a prognostic biomarker | | not run | not run | | not recommended |
| Residual permutation within path | approximate | $M$ times the procedure | not run | not run | good | proposed |

Table: Decision table for choosing the analysis half-life from the data (extension of Section 8, refined process, unless marked "ill."; revised 2026-10-10)

### 7.2 Plan

(Revised 2026-10-10, on the stress-test extension of Section 8. This
supersedes the provisional plan of 2026-10-09, which made AIC-weight
averaging with $\ell_1$ weights the primary data-driven candidate, the
bootstrap-calibrated max-$|t|$ test the calibration reference, F1A0 a
partial fix and F1G a post hoc candidate, on the illustration of
Section 6 alone. That plan had in turn replaced an earlier one that
made F1G primary and $\ell_0$-weight averaging the fallback, because
both had been defined after the illustration's results were seen. At
scale, on fresh trials under the refined process, AIC-weight averaging
falls just below the size band, F1A0 holds size but gains no power as
a nominal test, F1G behaves as F1, and the bootstrap calibration of a
data-chosen rule holds size and carries the power gain where 0.5 week
is wrong.)

The study is exploratory: the plan rests on one design (Hybrid), mostly
$N = 70$, and a warp-speed estimate of the bootstrap arms.

1. **Primary data-driven candidate: F1A0B.** Choose $t_{1/2}^a$ over
   $\mathcal T_a$ by ML without the interaction term ($\hat h_0$,
   Section 4.2), refit F1 by REML at $\hat h_0$ and compute the CR2
   $|T|$; refer it to $M = 999$ resamples from the null model FM0 at
   $\hat h_0$, each rerunning the selection and the test (Section 4.5).
   In the extension it rejected 0.050 pooled over the 12 non-prognostic
   null cells, at most 0.063 in one of them; its mean size-adjusted
   power was 0.395 against F1's 0.388, 0.500 against 0.352 at true
   $t_{1/2} = 2$ weeks and $c_{bm} = 0.45$, and its maximin shortfall
   0.076 against 0.148 (Section 8). Its selection never sees the
   interaction and its calibration includes the selection.
2. **Alternative: F1AB.** The same bootstrap applied to F1A (selection
   with the interaction in the model): 0.050 pooled, mean power 0.400,
   0.456 at $t_{1/2} = 2$ weeks and $c_{bm} = 0.45$, maximin 0.078. It is
   at the 0.075 cap in the prognostic-baseline cell, which is Monte
   Carlo error (Section 8.5); the reported half-life is the one that
   maximizes the full model's likelihood.
3. **Validated simple fallback: F1 at a fixed half-life.** When a
   bootstrap is not feasible, F1 with $t_{1/2}^a$ fixed before the data
   are seen holds size (0.053 pooled; 0.051 in `docs/47`) and loses
   power only when the fixed value is far from the truth. The
   pharmacological 0.2 week gave no power gain over 0.5 week (mean
   0.376 against 0.388; behind F1 in 8 of 9 significant paired
   comparisons), so a fixed value need not be pushed toward the drug's
   plasma half-life for power; 0.5 week is the stress test's base
   setting, not prazosin's pharmacological half-life, which report 01
   (Methods) puts at 0.1 to 0.2 week.
4. **Not recommended, for stated reasons.** F1A tested as if the
   half-life were fixed (0.060, fails); AIC-weight averaging without the
   bootstrap (MA 0.0455, slightly conservative; its bootstrap version
   MAB passes, 0.048, with the best mean power, 0.401, and is a
   further alternative to F1A0B); the max-$|t|$ test with a
   multivariate normal reference (0.072, liberal) or Bonferroni (0.018,
   too conservative); the multivariate $t$ reference (passes pooled, but
   flagged at 0.077 in one prognostic cell and lowest in mean power);
   F1G and MA0B, which are post hoc (MA0B has the best maximin
   shortfall, 0.062, but was defined after the illustration and should
   not be primary on that ground).
5. **The warp-speed caveat.** Every bootstrap figure is a warp-speed
   estimate: one resample per trial, pooled within a cell. Its Monte
   Carlo error is about 0.010 per cell rather than the binomial 0.0069
   (Section 8.1), and it assumes the bootstrap distribution varies
   little within a cell, which is less secure under the alternative.
   The power gain of F1A0B over F1A0 (same statistic, different
   critical value) is the part most exposed to that assumption. A
   check with $M = 999$ resamples per trial in a few cells
   (the base null cell, $t_{1/2} = 2$ at $c_{bm} = 0.45$ and its null
   twin) would cost about 17 worker-minutes per trial for F1A0B
   (Section 8.4) and is the next step before the recommendation is
   written into a paper.
6. **The prognostic caveat.** In the prognostic-baseline null cell the
   bootstrap arms reached 0.067 to 0.077 against F1's 0.057. That cell
   is, for every statistic, a non-prognostic null cell with a 4.6%
   smaller baseline SD (derived, Section 8.5), and the excess is Monte
   Carlo error. For prognostic pathways that vary over visits the null
   model reproduces the biomarker association only through $B^c_i$ and
   $B^c_i\eta_t$ (Section 4.5); in the nine other prognostic null cells
   the bootstrap arms stayed within the cap.
7. **The main effect is different** (`docs/58`). AIC selection with the
   tested term in the model is, for the main effect exactly and for the
   interaction approximately (Section 3.4, flat-$\ell_0$ case), a
   max-$|t|$ test judged against a single-test critical value, and is
   therefore liberal. For the main effect there is no analog of F1A0,
   and a fixed half-life does not by itself guarantee size unless time
   is modeled with visit effects (`docs/58`, Section 1, results 2 and
   3). Nothing in this plan transfers to the main effect without that
   qualification.
8. **Reporting, in every case.** The null-model ML estimate
   $\hat h_0$, the profile likelihood over a fine grid with the
   deviance at 0, 0.5, 1 and 2 weeks, the profile-likelihood interval
   (described as such, with unknown coverage), the half-life actually
   used and why, the bootstrap $p^*$ and critical value when F1A0B or
   F1AB is used, and F1 over the grid as a sensitivity analysis. State
   that the analysis half-life is a working parameter and that a
   Hybrid trial distinguishes a short half-life from one of a week or
   more but not values within $[0, 0.5]$ week.
9. **Not to be used:** F1A with its nominal reference, the max-$|t|$
   statistic with a single-$t$ reference, and naive permutation of the
   biomarker.

### 7.3 The stress-test extension: design

(Revised 2026-10-10: the extension has run; its results are in
Section 8. As implemented in
`analysis/scripts/quick-sim/carryover-closed-form/RERUN-ext3.md`, it
added five non-prognostic cells to the 12 below, the F1 arm at 0.2
week and the multivariate normal and $t$ calibrations of the grid
maximum, and it omitted cross-fitting and the profile-interval
coverage. The outline below is the design as written on 2026-10-09.)

- **Process.** The refined process (anchored carryover, randomized
  allocation, the driver's defaults), together with the pending
  anchored rerun of `docs/47`, so that F1 and the data-chosen rules
  are judged on the same trials.
- **Cells.** All 41 cells of `docs/46` (15 null at 1,000 replicates,
  26 alternative at 500; 28,000 trials), plus 12 new cells crossing a
  wrong half-life with a prognostic biomarker: the constant and
  growing natural-course pathways and the placebo pathway, each at
  true $t_{1/2} = 0$ and 1 week, null (1,000) and $c_{bm} = 0.30$
  (500); 9,000 more trials, 37,000 in all. The alternative cells with
  true $t_{1/2} = 0$ and 1 week are where a data-chosen half-life
  could recover power and must be included.
- **Arms.** F1, F1G, F1A0, AIC-weight averaging with $\ell_0$ and
  $\ell_1$ weights, cross-fit (average and Cauchy), F1AB and the
  max-$|t|$ test by warp-speed bootstrap, and the profile-interval
  coverage of the best-fitting grid value (estimated from a large
  noise-free fit per cell).
- **Cost.** The illustration's procedure took 7.1 worker-seconds per
  replicate with one bootstrap resample (Section 6.4), against about
  2.1 for the driver's `--extra2` arm (`RERUN-anchored.md`, Section 3;
  quoted). At 7.1 seconds, 37,000 replicates need about 263,000
  worker-seconds, 73 worker-hours: about 9 hours on 8 workers or 24
  hours on the 3 cores of a default `simrun` job (derived; assumes the
  server's cores match this machine's, which is not measured).
- **Criteria.** Those of `docs/46` as amended: pooled null rate within
  the 2-MCSE band of 0.05 over the non-prognostic null cells, no cell
  above 0.075, prognostic cells judged separately, and power compared
  on paired trials.

## 8. Stress-test extension results (2026-10-10)

### 8.1 What was run

The extension of Section 7.3 ran as `simrun` job 14 on the committed
tree of commit e077526: `11-strawman-stress-test.R --extra3`, specified
in `analysis/scripts/quick-sim/carryover-closed-form/RERUN-ext3.md`
(inspected). It has 58 cells and 40,000 trials: the 41 cells of
`docs/46` and 17 new ones (the 12 prognostic-by-wrong-half-life cells
of Section 7.3, the null and alternative cells at true $t_{1/2} = 2$
weeks, and the alternatives at $c_{bm} = 0.45$ with true $t_{1/2} = 0$,
1 and 2 weeks). There are 22 null cells (12 non-prognostic, 10
prognostic counting `prog-BL`) at $n_{\text{sim}} = 1{,}000$ and 36
alternative cells at 500. The trials use the refined process
(anchored carryover, randomized allocation) and fresh seeds (base
57261009), so they share no trial with `docs/47`, the anchored rerun or
Section 6. Every analysis is fitted to the same trials, so comparisons
between analyses are paired. The 18 analyses are those of
`RERUN-ext3.md`, Section 2: the eight nominal rules of Sections 4.1,
4.2, 4.7 and 4.8 plus F1 at the pharmacological 0.2 week (F1h02), the
max-$|T|$ statistic calibrated by a multivariate normal (MaxTMVN) or
$t$ (MaxTMVT) with the correlation of the $T(h)$ from the participants'
model-based influence contributions, and the warp-speed bootstrap of
eight of them (suffix B). F1G, MA0 and their bootstrap versions are
post hoc (Section 6.2); their rates here are on trials their
definition did not see.

The tables below are read from the decision files of
`12-strawman-stress-decision.R --ext3 --no-symmetry --tag -anchored`
(`decision-{screen,screen41,power,maximin,sensitivity}-f1-nosym-ext3-anchored.csv`)
and from `replicates-full-extra3-anchored.rds`, all in
`analysis/data/quick-sim/carryover-closed-form/strawman-stress/`
(verified). The screens are those of `docs/46` as amended: the pooled
null rate over the 12 non-prognostic null cells must lie in the band
0.046 to 0.054 (2 binomial MCSEs at 12,000 trials, 0.0040), no
non-prognostic cell may exceed 0.075, and a prognostic cell above 0.075
raises a flag. The symmetry screen was not applied (`--no-symmetry`).

**The bootstrap rates are warp-speed estimates.** Each trial drew one
null-model resample, and the B arms refer each trial's $|T|$ to the
pooled resampled statistics of its cell [5]. That estimates the
rejection rate of the bootstrap test, but its Monte Carlo error is
larger than the binomial one the screens use: resampling trials (pairs
of observed and resampled statistics) gives a per-cell MCSE of 0.0097
to 0.0112 for a rate near 0.05 at $n_{\text{sim}} = 1{,}000$, against
the binomial 0.0069 (verified, `57-progbl-check.R`, part 1), so about
0.003 rather than 0.002 for the 12-cell pooled rate (derived, assuming
independent cells). The warp-speed method also assumes that the
bootstrap distribution varies little between the trials of a cell [5]
(assumed, not checked). A real analysis would draw $M = 999$ resamples
for its own trial (Section 8.4).

### 8.2 Size

| Analysis | Pooled, 12 cells | Pooled, 11 cells of `docs/47` | Worst non-prognostic cell | `prog-BL` | Worst prognostic cell | Screens |
|------|------|------|------|------|------|------|
| F1, fixed 0.5 | 0.053 | 0.052 | 0.062 (Weibull $k = 2$) | 0.057 | 0.067 (`prog-PB`) | pass |
| F1h02, fixed 0.2 | 0.050 | 0.051 | 0.058 (Weibull $k = 2$) | 0.061 | 0.064 (`prog-PB`) | pass |
| F1A | 0.060 | 0.059 | 0.077 (MAR dropout) | 0.063 | 0.072 (`prog-PB`) | fails band and cap |
| F1A0 | 0.053 | 0.052 | 0.062 (MAR dropout) | 0.055 | 0.065 (`prog-PB`) | pass |
| F1G (post hoc) | 0.053 | 0.053 | 0.065 (Weibull $k = 2$) | 0.057 | 0.067 (`prog-PB`) | pass |
| MA, $\ell_1$ weights | 0.0455 | 0.045 | 0.059 (Weibull $k = 2$) | 0.042 | 0.063 (`prog-PB`) | fails band, conservative |
| MA0, $\ell_0$ weights (post hoc) | 0.040 | 0.040 | 0.058 (Weibull $k = 2$) | 0.039 | 0.059 (`prog-PB-thalf0`) | fails band, conservative |
| MaxTBonf | 0.018 | 0.018 | 0.026 (Weibull $k = 2$) | 0.024 | 0.029 (`prog-PB`) | fails band, conservative |
| MaxTMVN | 0.072 | 0.071 | 0.088 ($N = 35$) | 0.066 | 0.095 (`prog-TV-const-thalf0`) | fails band and cap, liberal |
| MaxTMVT | 0.053 | 0.052 | 0.063 (Weibull $k = 2$) | 0.050 | 0.077 (`prog-TV-const-thalf0`) | pass; prognostic flag |
| F1B | 0.050 | 0.050 | 0.059 ($N = 150$) | 0.077 | 0.077 (`prog-BL`) | pass; prognostic flag |
| F1h02B | 0.050 | 0.051 | 0.061 (Weibull $k = 0.5$) | 0.072 | 0.072 (`prog-BL`) | pass |
| F1AB | 0.050 | 0.050 | 0.064 (MAR dropout) | 0.075 | 0.075 (`prog-BL`) | pass, at the cap |
| F1A0B | 0.050 | 0.049 | 0.063 (MAR dropout) | 0.071 | 0.071 (`prog-BL`) | pass |
| F1GB (post hoc) | 0.051 | 0.050 | 0.065 (Weibull $k = 2$) | 0.076 | 0.076 (`prog-BL`) | pass; prognostic flag |
| MAB | 0.048 | 0.048 | 0.061 (Weibull $k = 2$) | 0.076 | 0.076 (`prog-BL`) | pass; prognostic flag |
| MA0B (post hoc) | 0.049 | 0.049 | 0.062 (Weibull $k = 2$) | 0.072 | 0.072 (`prog-BL`) | pass |
| MaxTB | 0.051 | 0.050 | 0.063 ($N = 150$) | 0.067 | 0.073 (`prog-PB`) | pass |

Table: Extension: null rejection rates at nominal 0.05 (verified; binomial MCSE 0.0020 pooled over 12,000 trials and 0.0069 in a cell; the B rows carry the larger warp-speed error of Section 8.1)

- **F1 at the fixed 0.5 week holds size under the refined process**
  (0.053 over 12 cells, 0.052 over the 11 of `docs/47`, against 0.051
  there), on fresh trials.
- **F1A fails again** (0.060 pooled, 0.077 under MAR dropout), as in
  `docs/47`; the failure replicates on independent trials.
- **F1A0 is at F1's rate at scale** (0.053; 0.052 over the 11 cells).
  The excess seen in the illustration (0.067 against F1's 0.061,
  Section 6.2) does not recur, which is consistent with its having been
  Monte Carlo error in a 1,000-trial illustration (argued).
- **Averaging is slightly conservative.** MA rejects 0.0455, just below
  the band's lower edge 0.0460; MA0 0.040. Their bootstrap versions MAB
  (0.048) and MA0B (0.049) pass.
- **The grid maximum needs the bootstrap.** Bonferroni is very
  conservative (0.018); the multivariate normal reference is liberal
  (0.072, 0.088 at $N = 35$), because it ignores the small-sample $t$
  tails of each $T(h)$ (argued); the multivariate $t$ on the smallest
  Satterthwaite df passes pooled (0.053) but is flagged in one
  prognostic cell (0.077).
- **Every bootstrap arm passes the pooled screen** (0.048 to 0.051), and
  each reaches its highest rate, 0.067 to 0.077, in the `prog-BL` cell,
  where F1 itself rejects 0.057. Three exceed the 0.075 flag there (F1B,
  F1GB, MAB) and F1AB sits at it. Section 8.5 shows that this is Monte
  Carlo error rather than a failure of the null model.

### 8.3 Power where the fixed value is wrong

| Cell | F1 | F1h02 | F1A0 | MA | F1AB | F1A0B | MAB | MA0B | MaxTB |
|------|------|------|------|------|------|------|------|------|------|
| base, $c_{bm} = 0.30$ ($t_{1/2} = 0.5$) | 0.438 | 0.428 | 0.424 | 0.438 | 0.478 | 0.464 | 0.478 | 0.470 | 0.448 |
| $t_{1/2} = 0$, $c_{bm} = 0.30$ | 0.488 | 0.494 | 0.500 | 0.508 | 0.546 | 0.550 | 0.526 | 0.528 | 0.504 |
| $t_{1/2} = 0$, $c_{bm} = 0.45$ | 0.830 | 0.834 | 0.828 | 0.842 | 0.860 | 0.868 | 0.866 | 0.868 | 0.834 |
| $t_{1/2} = 1$, $c_{bm} = 0.30$ | 0.282 | 0.252 | 0.272 | 0.302 | 0.304 | 0.264 | 0.314 | 0.312 | 0.314 |
| $t_{1/2} = 1$, $c_{bm} = 0.45$ | 0.564 | 0.516 | 0.568 | 0.574 | 0.570 | 0.584 | 0.574 | 0.580 | 0.556 |
| $t_{1/2} = 2$, $c_{bm} = 0.30$ | 0.188 | 0.176 | 0.206 | 0.194 | 0.226 | 0.218 | 0.236 | 0.236 | 0.214 |
| $t_{1/2} = 2$, $c_{bm} = 0.45$ | 0.352 | 0.348 | 0.428 | 0.412 | 0.456 | 0.500 | 0.454 | 0.440 | 0.420 |
| mean over 26 non-prognostic alternative cells | 0.388 | 0.376 | 0.381 | 0.387 | 0.400 | 0.395 | 0.401 | 0.399 | 0.388 |

Table: Extension: size-adjusted power, each analysis rejecting at the 5% quantile of its own $p$-values in the matching null cell ($n_{\text{sim}} = 500$; MCSE at most 0.022; verified from `decision-power-f1-nosym-ext3-anchored.csv`)

| Analysis | Maximin shortfall, 26 non-prognostic cells | Worst cell | Mean size-adjusted power | Maximin, all 36 cells | McNemar against F1: significant cells, ahead / behind |
|------|------|------|------|------|------|
| MA0B (post hoc) | 0.062 | $N = 150$, 0.30 | 0.399 | 0.090 | 18: 12 / 6 |
| F1A0B | 0.076 | MAR dropout, 0.30 | 0.395 | 0.092 | 12: 9 / 3 |
| F1AB | 0.078 | base curve, 0.375 | 0.400 | 0.092 | 13: 10 / 3 |
| MaxTB | 0.086 | $N = 150$, 0.30 | 0.388 | 0.132 | 11: 6 / 5 |
| MAB | 0.088 | MAR dropout, 0.30 | 0.401 | 0.102 | 17: 12 / 5 |
| F1GB (post hoc) | 0.110 | $t_{1/2} = 2$, 0.45 | 0.391 | 0.132 | 15: 10 / 5 |
| F1B | 0.126 | $t_{1/2} = 2$, 0.45 | 0.387 | 0.126 | 22: 13 / 9 |
| F1A0 | 0.126 | Weibull $k = 2$, 0.30 | 0.381 | 0.126 | 7: 2 / 5 |
| F1G (post hoc) | 0.142 | Weibull $k = 2$, 0.30 | 0.381 | 0.142 | 7: 2 / 5 |
| F1 | 0.148 | $t_{1/2} = 2$, 0.45 | 0.388 | 0.148 | (reference) |
| F1h02 | 0.152 | $t_{1/2} = 2$, 0.45 | 0.376 | 0.152 | 9: 1 / 8 |
| F1h02B | 0.158 | $t_{1/2} = 2$, 0.45 | 0.385 | 0.158 | 16: 7 / 9 |
| MaxTMVT | 0.160 | Weibull $k = 2$, 0.30 | 0.375 | 0.160 | 13: 2 / 11 |

Table: Extension: maximin shortfall in size-adjusted power among the analyses that passed the nominal screen, and paired McNemar comparisons with F1 over the 26 non-prognostic alternative cells (verified; MA, MA0, MaxTBonf, MaxTMVN and F1A failed the screen and are not ranked)

- **No rule gains much on average.** Mean size-adjusted power lies
  between 0.375 and 0.401 for every survivor; F1 is at 0.388.
- **The gains are where the fixed value is far from the truth.** At
  true $t_{1/2} = 2$ weeks and $c_{bm} = 0.45$, F1 has 0.352, F1A0
  0.428, F1AB 0.456 and F1A0B 0.500; at $t_{1/2} = 0$ and $c_{bm} =
  0.30$, F1 0.488 and F1A0B 0.550. At $t_{1/2} = 1$ week the rules are
  within the Monte Carlo error of one another.
- **F1A0 needs the bootstrap to gain.** The nominal F1A0 is behind F1 in
  more significant paired comparisons (5) than it is ahead (2), and its
  mean power (0.381) is below F1's; with the bootstrap calibration, F1A0B
  is ahead in 9 of 12 and F1AB in 10 of 13. The statistic of F1A0B is
  that of F1A0, so the difference comes from the critical value, which
  in the warp-speed estimate is the cell's pooled resampled quantile
  and in a real analysis the trial's own (inspected). How much of this
  gain survives with a per-trial bootstrap of $M = 999$ resamples, where
  the critical value varies from trial to trial, is not known (the
  warp-speed assumption of Section 8.1).
- **The pharmacological 0.2 week gives no gain over 0.5.** F1h02 has a
  mean power of 0.376 against 0.388 and is behind F1 in 8 of the 9
  cells where they differ significantly (ahead in 1).
- **The maximin ranking resolves groups, not places.** Each power has
  an MCSE of up to 0.022, and the maximin shortfall is a maximum over
  26 cells of differences between noisy estimates, so it is biased
  upward and differences below about 0.04 between analyses are not
  resolved (argued). MA0B, F1A0B, F1AB, MaxTB and MAB (0.062 to 0.088)
  form one group, ahead of F1 and the other nominal rules (0.126 to
  0.160).

**Selection.** Over all 58 cells, F1A chose 0, 0.25, 0.5, 1 and 2
weeks in 33%, 13%, 19%, 15% and 20% of trials, and F1A0 in 34%, 13%,
19%, 14% and 20% (verified). In the non-prognostic null cells F1A0
chose 0 in 85% of trials when the truth was 0, 1 or 2 weeks in 70%
when it was 1 week and 2 weeks in 84% when it was 2, but spread its
choices 33%, 17%, 25%, 14% and 12% when it was 0.5 (verified). The
guard F1G kept 0.5 week in 82% of trials and switched in 18%.

### 8.4 Cost

The `secs` column of the replicate file records the wall time of each
trial on one `simrun` worker: the full set of analyses on the data plus
one resample on which every bootstrapped procedure was rerun. Its
median is 25.7 seconds at $N = 70$ (10th and 90th percentiles over all
58 cells 22.4 and 29.1 seconds), 14.8 seconds at $N = 35$ and 49.0 seconds at $N =
150$ (verified). The resample is roughly half of that work (derived
from the operation counts of `RERUN-ext3.md`, Section 5), so a real
analysis with $M = 999$ resamples of all the bootstrapped procedures
would take about $999 \times 13 \approx 13{,}000$ worker-seconds, 3.6
worker-hours or about 30 minutes on 8 workers, at $N = 70$ (derived).
For F1A0B alone each resample needs five ML fits without the
interaction, one REML fit with it and one CR2 computation. On the
machine of `docs/47`, the check of Section 8.5 took a median of 4.6
seconds per trial under heavy load (verified) for about four such units
plus three null-model fits, roughly 1 second per unit, so $M = 999$
would take about 17 worker-minutes (derived). Inside a simulation study, a full bootstrap
per trial multiplies the cost by $M$, which is why the extension used
the warp-speed method.

### 8.5 The prognostic-baseline cell

The B arms reach 0.067 to 0.077 in the `prog-BL` null cell (cell 40;
biomarker correlated $c_{bm,BL} = 0.3$ with the baseline level
$\mathrm{BL}_i$, no interaction), against 0.057 for F1 on the same
trials. The hypothesis to be checked was that the null model, which
reproduces a biomarker association with $Y$ only through $\beta_B B^c_i$
and $\beta_{bm:\eta} B^c_i\eta_t$, does not reproduce the association
with $\mathrm{BL}_i$, so that the resamples lack the prognostic
structure and the bootstrap reference is too narrow. Full details are
in `57-progbl-check.R`; the findings are these.

- **The null model does reproduce it, exactly** (derived; inspected in
  the driver's `build`, `draw` and `make_trial`). In the driver,
  $Y_{i0} = \mathrm{BL}_i$ and $Y_{it} = \mathrm{BL}_i - (TV_{it} +
  PB_{it} + BR_{it})$, so $\mathrm{BL}_i$ enters every row. In a null
  `prog-BL` trial the biomarker is correlated with $\mathrm{BL}_i$ only,
  and joint normality gives $\mathrm{BL}_i = \mu_{BL} + \gamma(B_i -
  \mu_B) + e_i$ with $\gamma = c_{bm,BL}\sigma_{BL}/\sigma_{bm}$ and $e_i$
  normal with variance $\sigma_{BL}^2(1 - c_{bm,BL}^2)$, independent of
  $B_i$, $TV$, $PB$ and $BR$. The prognostic term $\gamma B_i$ is the same
  at every visit and lies in the span of the intercept and $B^c_i$. REML
  and ML likelihoods, the residuals, the CR2 standard errors and the
  Satterthwaite df are invariant to adding a vector in the span of the
  fixed effects, and so is the resampled statistic, whose mean
  $W_i\hat\xi$ is in that span. Every statistic of every procedure,
  observed or resampled, therefore has the same distribution in
  `prog-BL` as in a non-prognostic null cell whose baseline SD is
  $\sigma_{BL}\sqrt{1 - 0.09} = 0.954\,\sigma_{BL}$. The baseline-specific
  slope $\beta_{B,0}$ (`bmc:bl0`) is zero in this process.
- **The stored rates are consistent with Monte Carlo error**
  (verified, `57-progbl-check.R`, part 1). With the trial-resampling
  MCSE of about 0.010, F1B's 0.077 is 2.5 MCSEs above the 22-cell mean
  of 0.053; the largest of 22 such deviations is expected near 2
  (derived, normal approximation), and the spread of the B rates over
  the 22 null cells (SD 0.0103 to 0.0119) is no larger than their MCSE,
  so there is no sign of real between-cell differences. The base null
  cell, which differs from `prog-BL` only by the 4.6% in $\sigma_{BL}$,
  has the opposite extreme: F1B 0.032, F1AB 0.033. In `prog-BL` the
  resampled 95% quantile is the lowest of the 22 cells for seven of the
  eight statistics and the second lowest for F1A0 (F1: 1.89 against a mean of 2.07, SD 0.077, over the
  other 21 cells), while the observed $|T|$ quantile is ordinary (2.09,
  F1's nominal rate 0.057); the resampled draws of the eight procedures
  share each trial's resample, so their low quantiles are one event,
  not eight.
- **Fresh trials** (verified, `57-progbl-check.R`, part 2; 200 new
  trials in each of the base null cell and `prog-BL`, refined process,
  seeds $57990000 + 10000 \times \text{cell} + r$). Each trial was
  scored by F1 and F1A0 and given one resample from each of three null
  models, rerunning F1A0's selection on the resample: the stored
  parametric model (P0), the same with the baseline-specific slope
  `bmc:bl0` (P1), and a wild cluster bootstrap that multiplies each
  participant's null-model residual vector by a Rademacher sign (W).
  The run was cut to 200 trials per cell because of machine load, so
  the warp-speed sizes carry trial-resampling MCSEs of 0.016 to 0.043
  and resolve only large effects.

| Null model | Statistic | Base null cell | `prog-BL` |
|------|------|------|------|
| P0, stored parametric | F1 | 0.050 (0.031) | 0.065 (0.027) |
| P0, stored parametric | F1A0 | 0.055 (0.035) | 0.055 (0.022) |
| P1, with `bmc:bl0` | F1 | 0.025 (0.025) | 0.060 (0.027) |
| P1, with `bmc:bl0` | F1A0 | 0.025 (0.016) | 0.065 (0.034) |
| W, wild cluster | F1 | 0.115 (0.036) | 0.045 (0.035) |
| W, wild cluster | F1A0 | 0.095 (0.043) | 0.065 (0.036) |
| nominal CR2 (no bootstrap) | F1 | 0.035 (0.015) | 0.045 (0.015) |

Table: Fresh-trial check of the prognostic-baseline cell: warp-speed size of each null-model bootstrap at nominal 0.05, $n_{\text{sim}} = 200$ per cell (trial-resampling MCSE in parentheses; binomial for the nominal row; verified)

Under the stored parametric model the two cells do not differ
detectably (0.050 and 0.065 for F1, 0.055 and 0.055 for F1A0), and
the 95% quantile of the resampled F1 statistic is 1.95 in the base
cell and 1.91 in `prog-BL`, a difference well inside the sampling
error of a quantile from 200 draws, about 0.13 (derived). Adding
`bmc:bl0` to the null model (P1) does not lower the `prog-BL` rate, as
the derivation predicts. The wild cluster bootstrap is not a remedy:
in the base null cell its resampled statistics are too narrow (95%
quantile 1.67 for F1 against 2.00 for the observed $|T|$) and it
rejects 0.115, consistent with the known shrinkage of unadjusted
residuals (argued).

**Conclusion.** The `prog-BL` flag on F1B, F1GB and MAB, and F1AB's
0.075, are Monte Carlo error of the warp-speed estimate, not a failure
of the null model to reproduce a baseline association: no remedy is
needed for a biomarker associated with baseline severity alone. The
equivalence is derived; the fresh-trial check supports it but, at 200
trials per cell, cannot exclude a difference of about 0.03 between the
cells. The
caveat of Section 4.5 still stands for prognostic pathways that change
across visits, which the null model reproduces only through $B^c_i$ and
$B^c_i\eta_t$; in the six prognostic null cells with a wrong half-life
the B arms reject 0.037 to 0.073, within the cap.

### 8.6 Reading

- The extension confirms that a fixed half-life holds size and that
  AIC selection with the interaction in the model (F1A) does not.
- Choosing the half-life without the interaction (F1A0) holds size at
  scale, but as a nominal test it gains no power over F1.
- The bootstrap calibration of a data-chosen rule holds size (pooled
  0.048 to 0.051) and is where the power gain lies, in the cells whose
  true half-life is far from 0.5 week. This rests on the warp-speed
  estimate and its assumption.
- Among the bootstrap arms, F1A0B and F1AB are specified procedures
  (Sections 4.2 and 4.5) with the best maximin values after the post hoc
  MA0B; MAB and MaxTB are close; their differences are not resolved.
- Averaging without the bootstrap is slightly conservative; the
  analytic calibrations of the grid maximum are either liberal (MVN) or
  conservative (Bonferroni); and the pharmacological 0.2 week brings no
  gain over 0.5.

## 9. Limitations

(Revised 2026-10-10 for the extension of Section 8.)

- **Two settings.** Sections 1 to 7 quote the stored stress test and
  the illustration of Section 6, which use configuration C at
  $\rho = 0.7$ with the published recursive carryover and fixed
  allocation. The extension of Section 8 uses the refined process
  (anchored carryover, randomized allocation) on fresh trials, so its
  figures are not paired with those of `docs/47` or Section 6.
- **Warp-speed bootstrap.** Every bootstrap rate of Section 8 rests on
  one resample per trial, pooled within a cell. Its Monte Carlo error is
  about 0.010 per cell (Section 8.1), and the method assumes that the
  bootstrap distribution varies little between the trials of a cell,
  which was not checked and is less secure under the alternative. No
  per-trial bootstrap with $M = 999$ has been run except on the single
  trial of Section 5 ($M = 199$).
- **Power at 500 trials per cell.** The size-adjusted powers carry an
  MCSE of up to 0.022, and maximin differences below about 0.04 are not
  resolved (Section 8.3).
- **Arguments, not proofs.** The validity arguments for F1A0
  (Section 4.2) and the orthogonality result (Section 4.6) are first
  order. They do not cover a prognostic biomarker, informative
  dropout, or the boundary at $t_{1/2}^a = 0$, where the example
  trial's estimate lay. The extension found F1A0 at F1's size in every
  null cell run, including MAR dropout (0.062) and six prognostic
  cells with a wrong half-life, which is evidence for, not a proof of,
  the argument.
- **Post hoc rules.** F1G, averaging with $\ell_0$ weights and their
  bootstrap versions were defined after the illustration's results
  were seen. The extension evaluated them on fresh trials, but MA0B's
  leading maximin value is still the rate of a rule chosen after
  seeing related results, and it is not made primary.
- **Prognostic pathways.** The prognostic-baseline cell is shown in
  Section 8.5 to be equivalent to a non-prognostic cell. The other nine
  prognostic null cells held the cap for every bootstrap arm, but the
  null model reproduces a time-varying prognostic association only
  through $B^c_i$ and $B^c_i\eta_t$, and pathways other than those of
  `docs/46` were not run.
- **One design, mostly one sample size.** Hybrid at $N = 70$, with
  $N = 35$ and 150 in one null and one alternative cell each.
- **Not implemented:** joint estimation with a sandwich that includes
  $\partial D_{bc,it}/\partial\omega$ (Section 4.6), the Bayesian
  analysis (Section 4.8), the residual permutation test (Section 4.9),
  the profile interval's coverage (Section 4.4), and cross-fitting at
  scale (Section 4.3).

## 10. How to reproduce

From the repository root:

```{.bash}
Rscript analysis/scripts/data-chosen-half-life/57-example-trial.R --cores 7 --boot 199
Rscript analysis/scripts/data-chosen-half-life/57-null-illustration.R --cores 7
Rscript analysis/scripts/data-chosen-half-life/57-posthoc-rules.R
Rscript analysis/scripts/data-chosen-half-life/57-stored-tables.R
```

The extension of Section 8 (`RERUN-ext3.md`, Section 6), its decision
files, and the check of Section 8.5:

```{.bash}
S=analysis/scripts/quick-sim/carryover-closed-form
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript $S/11-strawman-stress-test.R --extra3 --cores 16
Rscript $S/12-strawman-stress-decision.R --ext3 --no-symmetry --tag -anchored
Rscript analysis/scripts/data-chosen-half-life/57-progbl-check.R --reps 1000 --cores 6
```

Both source `57-common.R`, which evaluates the top-level definitions
of `analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R`
(everything before its pre-run checks, skipping the lines that create
output directories), so the data-generating process and the fitting
helpers are the driver's own.

## 11. Reproducibility

Every number labeled verified in Sections 1 to 8 was computed by the
scripts in `analysis/scripts/data-chosen-half-life/`:

- `57-common.R`: the procedures (grid fits, F1, F1A, F1A0, AIC-weight
  averaging, max-$|t|$, split and cross-fit, null-model bootstrap).
- `57-example-trial.R` and `57-example-figure.R`: Section 5, Tables 4
  and 5 and Figure 1; outputs `57-example-*.csv` in
  `analysis/data/quick-sim/data-chosen-half-life/`.
- `57-null-illustration.R`: Section 6, Tables 6 to 8 and Figure 2; outputs
  `57-mc-*.csv`, `57-mc-replicates.rds` and the checkpointed chunks in
  `57-mc-chunks/`.
- `57-posthoc-rules.R`: the post hoc rows of Table 6 (F1G, averaging
  with $\ell_0$ weights), the pooled null rates and the paired
  discordance counts; outputs `57-mc-posthoc-rules.csv`,
  `57-mc-pooled-null.csv` and `57-mc-discordance.csv`.

- `57-progbl-check.R`: Section 8.5 (part 1 reads the extension's
  replicate file; part 2 fits fresh trials of cells 2 and 40); output
  `progbl-check.rds` in
  `analysis/data/quick-sim/carryover-closed-form/strawman-stress/`.

The tables of Section 8 are read from the decision files
`decision-*-f1-nosym-ext3-anchored.csv` and from
`replicates-full-extra3-anchored.rds` in the same directory, written by
`11-strawman-stress-test.R --extra3` (`simrun` job 14) and
`12-strawman-stress-decision.R --ext3 --no-symmetry --tag -anchored`;
the McNemar counts and selection shares are those the decision script
prints (recomputed here from the replicate file with its rules), and
the `secs` quantiles were read directly from the replicate file.

Tables 2 and 3 were computed by `57-stored-tables.R` from the stored
stress-test replicates
(`analysis/data/quick-sim/carryover-closed-form/strawman-stress/replicates-full.rds`
and `replicates-full-extra2.rds`) by matching cell and replicate.

## 12. References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker
   validation: statistical methods and theoretical findings.
   *Frontiers in Digital Health*. 2020;2:13.
   doi:10.3389/fdgth.2020.00013.
2. Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS*. New
   York: Springer; 2000. doi:10.1007/b98882.
3. Bell RM, McCaffrey DF. Bias reduction in standard errors for linear
   regression with multi-stage samples. *Survey Methodology*.
   2002;28(2):169-181.
4. Pustejovsky JE, Tipton E. Small-sample methods for cluster-robust
   variance estimation and hypothesis testing in fixed effects models.
   *Journal of Business and Economic Statistics*. 2018;36(4):672-683.
   doi:10.1080/07350015.2016.1247004.
5. Giacomini R, Politis DN, White H. A warp-speed method for
   conducting Monte Carlo experiments involving bootstrap estimators.
   *Econometric Theory*. 2013;29(3):567-589.
   doi:10.1017/S0266466612000655.
6. Davies RB. Hypothesis testing when a nuisance parameter is present
   only under the alternative. *Biometrika*. 1977;64(2):247-254.
   doi:10.2307/2335690.
7. Davies RB. Hypothesis testing when a nuisance parameter is present
   only under the alternative. *Biometrika*. 1987;74(1):33-43.
   doi:10.1093/biomet/74.1.33.
8. Cox DR. A note on data-splitting for the evaluation of significance
   levels. *Biometrika*. 1975;62(2):441-444.
   doi:10.1093/biomet/62.2.441.
9. Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C, Newey
   W, Robins J. Double/debiased machine learning for treatment and
   structural parameters. *The Econometrics Journal*.
   2018;21(1):C1-C68. doi:10.1111/ectj.12097.
10. Liu Y, Xie J. Cauchy combination test: a powerful test with
    analytic p-value calculation under arbitrary dependency structures.
    *Journal of the American Statistical Association*.
    2020;115(529):393-402. doi:10.1080/01621459.2018.1554485.
11. Self SG, Liang KY. Asymptotic properties of maximum likelihood
    estimators and likelihood ratio tests under nonstandard
    conditions. *Journal of the American Statistical Association*.
    1987;82(398):605-610. doi:10.1080/01621459.1987.10478472.
12. Hansen BE. Inference when a nuisance parameter is not identified
    under the null hypothesis. *Econometrica*. 1996;64(2):413-430.
    doi:10.2307/2171789.
13. Buckland ST, Burnham KP, Augustin NH. Model selection: an integral
    part of inference. *Biometrics*. 1997;53(2):603-618.
    doi:10.2307/2533961.
14. Hoeting JA, Madigan D, Raftery AE, Volinsky CT. Bayesian model
    averaging: a tutorial. *Statistical Science*. 1999;14(4):382-417.
    doi:10.1214/ss/1009212519.
15. Freedman D, Lane D. A nonstochastic interpretation of reported
    significance levels. *Journal of Business and Economic
    Statistics*. 1983;1(4):292-298. doi:10.1080/07350015.1983.10509354.
16. Leeb H, Pötscher BM. Model selection and inference: facts and
    fiction. *Econometric Theory*. 2005;21(1):21-59.
    doi:10.1017/S0266466605050036.
17. Morris TP, White IR, Crowther MJ. Using simulation studies to
    evaluate statistical methods. *Statistics in Medicine*.
    2019;38(11):2074-2102. doi:10.1002/sim.8086.


---

## Derivations

<!-- one row per result: statement | CM label | verification script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the numbers agree? -->

---

## Limitations and open items

<!-- link each to the open-items register at docs/series/registry/open-items.csv -->

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
