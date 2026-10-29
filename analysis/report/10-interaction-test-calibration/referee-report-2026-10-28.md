# Referee report: Type I calibration of fixed-effect interaction tests in linear mixed models (paper 10)
*2026-10-28 22:35 PDT*

**Author.** pmsimstats team (referee pass)

**Manuscript.** `analysis/report/10-interaction-test-calibration/report.Rmd`
(last modified 2026-09-07), with `references.bib`, `title-page.md` and
`type1-examination-plan.md`. `bullets.Rmd` exists in the same directory
and was not reviewed.

**Target journal.** *Statistics in Medicine* (methods contribution). The
paper is a poor fit for *Trials*, which would require an applied trial
question rather than a methodological one.

**Evidence labels.** VERIFIED: I executed code and observed the result.
INSPECTED: I read the source or artifact and confirmed the point.
SUSPECTED: the point is consistent with what I read but was not
confirmed. UNVERIFIED: not checked.

**Scratch scripts** (all under the session scratchpad
`/private/tmp/claude-501/-Users-zenn-prj-res-36-pmsimstats-ng-pmsimstats-ng-analysis-report-01-dgp-mean-moderation-vs-mvn/513155aa-ec85-4a1d-889a-5d57e8607d18/scratchpad/`):
`p10-verify.R` (reproduces every number from the stored `.rds`),
`p10-kappa-analytic.R` and `p10-kappa-lib.R` (exact population, or
pseudo-true, $\kappa$ from the DGP covariance and mean),
`p10-kappa-meanfit-run.R` (the same with the mean model made correct),
`p10-meanmisfit-sim.R` (a paired 1,500-replicate simulation check),
`p10-misc.R`.

## 1. Summary

The manuscript asks why the model-based Wald test of a biomarker by
drug interaction in an aggregated N-of-1 simulation rejects at about
0.03 under the null. The analysis is `nlme::lme` with a random
intercept and a `corCAR1` residual (called M3). The authors set out a
"four-channel decomposition" of the realized Type I error: bias,
standard-error scale $\kappa$, degrees of freedom, and reference shape.
They add a working-covariance ladder (M0 ordinary least squares, M1
random intercept, M2 `corCAR1` only, M3 both) and a sample-size
gradient ($N$ = 35 to 280). They conclude four things:

- $\kappa$ is 1.06 to 1.10.
- The `corCAR1` layer causes the excess.
- The excess does not shrink with $N$, so it is asymptotic.
- A CR2 cluster-robust standard error with Satterthwaite degrees of
  freedom restores the size to about 0.05.

They also report three side checks. Kenward-Roger is inert on an
autoregressive marginal fit. A GEE with the Mancl-DeRouen correction
also approximately calibrates. At one alternative cell, CR2 raises
power by 3 to 4 points without reordering the specifications.

The arithmetic of the paper is sound. Every reported number that I
checked reproduces from the stored artifacts, apart from one misquoted
range. The design is the four-path Hendrickson Hybrid, as claimed. No
stored cell needed a positive-definiteness repair. The CR2 remedy works
as reported.

The central mechanistic claim is wrong. The excess $\kappa$ of M3 does
not come mainly from a working covariance that "cannot match a mixture
of autoregressive processes". It comes mainly from the misspecified
**mean** model. The fixed effects $\beta_0 + \beta_1 b + \beta_2 t +
\beta_3 X + \theta\, b X$ cannot follow the Gompertz natural-history,
placebo and drug-response mean trajectories. The residual mean pattern
is common to every participant on a path, and REML absorbs it into the
variance parameters. I show this analytically and confirm it with a
paired simulation (Major comment 1). This changes the diagnosis, the
explanation of the ladder, and the paper's main practical lesson.

The stored runs also predate two corrections to the project's DGP: the
cross-factor correlation constant (`docs/61`) and anchored carryover.
The manuscript never mentions either (Major comment 2).

## 2. Overall assessment and recommendation

**Recommendation: reject in the current form, with encouragement to
resubmit as a substantially revised paper.** If the editor prefers a
graded decision, this is a major revision. The required revision is
not cosmetic. The causal story in Sections 4.2, 4.3 and 5.1 to 5.2 has
to be rewritten around mean-model misspecification. The headline runs
have to be regenerated under the current DGP. The novelty claims have
to be scaled back to what the literature leaves open.

The CR2 conclusion survives. It is the most useful result in the
paper, and it agrees with paper 02's recommended analysis (G9) and
with paper 08.

Strengths, briefly:

- The stored artifacts are complete at the replicate level.
- The pairing of seeds across analyses is genuine, so paired
  comparisons are available.
- The ladder, CR2 and Kenward-Roger computations are correctly coded.

## 3. Significance and novelty

**What is known.** A misspecified working covariance leaves the
model-based fixed-effect standard error inconsistent, and the sandwich
repairs it. This is classical: Liang and Zeger (1986), Royall (1986),
White (1982). Gurka, Edwards and Muller (2011, *Stat Med*) proved that
underspecifying the covariance biases mixed-model fixed-effect
inference in large as well as small samples. They state explicitly
that "recruiting more participants will not alleviate" it. Lu and
Mehrotra (2010, *Stat Med*) and Gosho et al. (2017, *SMMR*) study the
sandwich and its bias corrections against Kenward-Roger under
misspecified MMRM covariances. None of these three papers is cited.

Given this literature, the "asymptotic" finding (Section 4.3) is not
new. The manuscript's own Section 2.3 also states that $\kappa$ is a
per-cluster functional and constant in $N$ by construction. A flat
gradient therefore confirms a mathematical identity. It is not an
empirical discovery, yet the paper presents it as one of "two features
that are not obvious in advance".

**The four-channel decomposition.** The decomposition restates the
standard performance measures of Morris, White and Crowther (2019):
bias, the ratio of model-based to empirical SE, and rejection rate. The
manuscript concedes this (lines 349 to 382). The claims of
exhaustivity and independent estimability are also not correct as
stated (Major comment 4).

**What could be new.** Two things:

- A quantified demonstration that mean-model misspecification, through
  REML absorption, makes a mixed-model interaction test conservative.
- A demonstration that the conservatism and the apparent calibration
  of the random-intercept-only model are both artifacts of that
  absorption.

The current manuscript does not make either point. With them, the
paper would be a useful short methods note. Without them, its
contribution is a single case study of a known phenomenon with a
partly wrong explanation.

**Contemporariness.** The citation base stops at 2018 for the method
itself, apart from McNeish (2017). It omits Luke (2017), Gosho et al.
(2017), Wang et al. (2016), Li and Redden (2015), and the N-of-1
methodology literature (Zucker et al. 2010; Araujo, Julious and Senn
2016).

## 4. Major comments

### 4.1 The conservatism comes mainly from mean-model misspecification, not from the working covariance (correctness)

**Location.** Sections 3.2, 4.2, 4.3 (lines 526 to 535, 692 to 705,
810 to 822). Section 5.1 to 5.2. Abstract ("localizes the
miscalibration to the corCAR1 residual layer"). Conclusions, points 3
and 4.

**Problem.** In the manuscript's sandwich formula (eq. vsand), $V$ must
be the second moment of $y$ about the fitted working mean. Under a
misspecified mean, that moment is $V_p + r_p r_p^\top$, where $r_p$ is
the path-specific residual mean pattern. The analysis model has a
linear trend and a drug indicator. The DGP mean is a sum of modified
Gompertz curves with an expectancy-scaled placebo response and
carryover. The root-mean-square of $r_p$ is about 4.1 to 4.4 points
per path. The sum is $\sum_p r_p r_p^\top$. REML treats this pattern
as covariance.

I computed the exact pseudo-true M3 parameters and the population
$\kappa$ from the DGP's own covariance and mean for every Hybrid path
(`p10-kappa-analytic.R`).

The analytic machinery is validated against the stored simulation:

| Quantity (stored DGP, E1) | Analytic | Stored simulation |
|---|---|---|
| $\kappa$, M3 | 1.071 | 1.068 |
| $\kappa$, M3 (E2) | 1.082 | 1.097 |
| `corCAR1` $\rho_a$ | 0.659 | median 0.664 |
| $\sigma_u$ | 9.85 | median 9.69 |
| residual SD | 12.10 | median 12.15 |

A further check: with a homoscedastic placebo component, matched
cross-factor constants and a correct mean model, M3 is exactly
correctly specified. In that case the computation returns $\kappa =
1.000$, $\sigma_u = 10 = \sigma_{BL}$ and $\rho_a = 0.70$. This is
VERIFIED.

The decomposition of the excess:

| Scenario | $\kappa$ M3 (E1, E2, E3) | $\kappa$ M1 (E1, E2, E3) |
|---|---|---|
| Stored DGP, misspecified mean (as in paper) | 1.071, 1.082, 1.064 | 0.972, 0.949, 0.967 |
| Stored DGP, mean model correct | 1.017, 1.023, 1.017 | 0.901, 0.887, 0.901 |
| Refined DGP ($c_\times = c_1 = 0.2$), anchored, misspecified mean | 1.058, 1.062, 1.055 | 0.959, 0.935, 0.957 |
| Refined DGP, mean model correct | 1.006, 1.006, 1.006 | 0.888, 0.872, 0.888 |

A paired simulation confirms this. It used 1,487 replicates of the
primary null cell, analyzed with M3 E1 twice on the same data: once
on $Y$ and once on $Y$ minus its true path mean (`p10-meanmisfit-sim.R`).

- $\kappa$ fell from 1.072 to 1.014. The paired difference was 0.057,
  with a bootstrap SE of 0.0026.
- The Type I error rose from 0.033 to 0.047.

This is VERIFIED.

**Why it matters.**

1. The covariance structure accounts for at most 1 to 2 percent of the
   6 to 10 percent excess. That part comes from the placebo-variance
   heteroscedasticity and the `docs/61` lag-0 discontinuity. The mean
   misfit accounts for the rest.
2. The ladder's central reading is an artifact of offsetting errors.
   This applies to "M1, a random intercept alone, is well calibrated"
   and to "the less elaborate model was the calibrated one"
   (lines 1060 to 1063). Under a correct mean, M1 has $\kappa$ of 0.87
   to 0.90 and would be anti-conservative, as the plan originally
   expected (plan Section 6). It appears calibrated only because the
   absorbed mean misfit pushes it upward by a similar amount.
3. The practical lesson changes. A richer mean model reduces the
   conservatism, for instance with expectancy terms (the F1 model of
   `docs/56`) or with visit as a factor. Removing the `corCAR1` term
   does not.

**Remedy.**

- Rewrite Sections 3.2, 4.2, 4.3 and 5.1 to 5.2 around mean-model
  misspecification.
- Add the analytic pseudo-true $\kappa$ computation, or an equivalent
  large-$N$ fit, as a fourth instrument. It decomposes the excess into
  its mean and covariance parts exactly and without Monte Carlo error.
- Add ladder rungs that vary the mean model while holding the
  covariance fixed:
  - visit as a factor crossed with path, or a saturated mean;
  - the F1 mean, `bmc * (Dbc + De) + t`.
- State the generalized condition in Section 2.2: $V$ is replaced by
  $V + r r^\top$ when the mean is wrong.

**Evidence.** VERIFIED. Analytic computation, plus a paired simulation
of 1,487 replicates.

### 4.2 The stored runs predate two DGP corrections, and the manuscript does not say so (correctness, reproducibility)

**Location.** Section 3.2. Reproducibility section. Every table.

**Problem.** All seven drivers (`diagnose-null-type1*.R`,
`diagnose-alt-cell-cr2.R`) hard-code `c.cf1t = 0.2, c.cfct = 0.1`. This
is the cross-factor mismatch documented in `docs/61` (2026-10-21). The
refined DGP specifies a single $c_1$, so `c.cfct = 0.2`.

The `.rds` files were produced on 2026-05-20 and relabeled on
2026-07-31. In May, the tidyverse `apply_carryover_to_component()` used
the compounding recursion. The anchored rule replaced it later (git
diff from commit `d337bc90a` to HEAD). Rerunning the drivers today
would therefore not reproduce the stored numbers. The table captions
do not identify which DGP produced them.

`docs/61` (Section 8) asserts that these diagnostic scripts do not
"feed the manuscript's reported tables or figures directly". That is
true of paper 02 and false of paper 10.

**Why it matters.** Under the refined DGP the population $\kappa$ of M3
falls from 1.064 to 1.082 down to 1.055 to 1.062 (table in 4.1). The
predicted size moves from 0.031 to 0.036 up to about 0.037 to 0.038.
The direction of the finding holds. The quoted magnitudes do not
describe the project's current DGP.

**Remedy.**

- Set `c.cfct = 0.2` in all seven drivers.
- Rerun under the current package. The runs are minutes of compute
  per the plan's Section 12.
- State in Methods the DGP version, the commit, and the carryover
  rule.
- Correct the sentence in `docs/61` Section 8.

**Evidence.** INSPECTED for the driver code, the git history and the
`docs/61` text. VERIFIED for the effect on population $\kappa$, by
analytic computation.

### 4.3 The null point-estimate "bias" is Monte Carlo noise, not a property of the estimator (correctness)

**Location.** Lines 577 to 583. Limitations, lines 1244 to 1248.
Plan Section 13.1.

**Problem.** At $c_{bm} = 0$ the biomarker is $\mathcal N(0, 1)$ and
independent of the outcome. The REML fit is equivariant under
$b \to -b$: the interaction estimate changes sign and the variance
parameters do not change. Hence $\mathbb E(\hat\theta) = 0$ exactly.

The observed mean of 0.032 to 0.039 is 2.3 to 2.6 Monte Carlo SE from
zero (VERIFIED). It recurs across E1 to E3, M0 to M3 and the GEE
(0.026 to 0.039) only because every analysis was fitted to the same
5,000 datasets (VERIFIED: identical estimates for d13 and the M3 rung
of d4, and identical estimates for d13 and the CR2 run). This is one
chance fluctuation seen many times, not replicated evidence.

**Why it matters.** The Limitations section calls it "a genuine
finite-sample feature of the estimator that a separate examination
might characterize". That statement is false, and it invites wasted
follow-up work.

**Remedy.** Replace the statement with the symmetry argument. Report
the mean with its MCSE.

**Evidence.** VERIFIED for the numbers and the shared datasets.
INSPECTED for the derivation.

### 4.4 The four-channel decomposition is internally inconsistent and overclaimed (correctness)

**Location.** Section 2.1, equations (eq:sdt), (eq:alpha) and lines
383 to 393.

**Problems.**

1. $Z^\star = (\hat\theta - \mu)/\sigma$ has unit variance by
   construction (line 257). The text nevertheless reads
   $\mathrm{SD}(T) < 1$ as possibly "a sub-normal reference shape"
   (line 303). Under the stated approximation,
   $\mathrm{SD}(Z^\star) = 1$ always. Any gap between
   $\mathrm{SD}(T)$ and $1/\kappa$ comes from the variability of
   $\widehat{\mathrm{SE}}$ and its dependence on $\hat\theta$, which
   the approximation has discarded.
2. The text says the channels are "exhaustive", that "a discrepancy
   absent from all four cannot exist", and that they "are estimated
   independently". It also says (line 309) that they are "only
   approximately separable". The reference-shape channel is a
   residual catch-all, so exhaustivity holds by definition. The joint
   law of $(\hat\theta, \widehat{\mathrm{SE}})$ is not a channel at
   all.
3. This matters in practice. The GEE naive sandwich has
   $\kappa = 0.955$ to 0.963. The closed form (eq:closed) predicts a
   size of about 0.058 to 0.061, but 0.067 to 0.069 is observed
   (VERIFIED). CR2 has $\kappa$ of 0.976 to 1.000 but a size of
   0.045 to 0.049, because the coefficient of variation of the CR2 SE
   is 0.17 (VERIFIED). The promised return to "its limits where the
   cluster-robust standard error ... is taken up" (line 311) never
   occurs.

**Remedy.**

- Present the decomposition as a first-order heuristic.
- Remove the exhaustivity language.
- Add $\mathrm{CV}(\widehat{\mathrm{SE}})$ as an explicit fifth
  diagnostic.
- Show where the closed form fails: GEE naive, CR2.

**Evidence.** INSPECTED for the derivation. VERIFIED for the
numerical counterexamples.

### 4.5 The sample-size gradient cannot carry the "decisive" weight placed on it (correctness, design)

**Location.** Section 4.3. Abstract. Conclusions, point 4.

**Problem.** At 1,500 replicates each $\kappa$ has an MC SE of about
$\kappa/\sqrt{2n} \approx 0.02$. The fitted slopes of $\kappa$ on
$\log_2 N$ are -0.003, -0.008 and -0.002 (VERIFIED). These are far
below what the design can resolve: a slope of -0.01 per doubling would
go undetected.

The size estimates at $N = 70$ for E1 and E3 (0.043 and 0.045) are
within 1.3 MC SE of 0.05 (VERIFIED). So "no trend toward 0.05" rests
mainly on E2. In addition, Section 2.3 already argues that $\kappa$ is
constant in $N$ by construction.

**Remedy.**

- Report the slope with a confidence interval.
- State that the gradient checks for finite-sample REML terms only.
- Use the analytic pseudo-true $\kappa$ (Major comment 1) as the
  decisive evidence that the excess is asymptotic.

**Evidence.** VERIFIED.

### 4.6 The Kenward-Roger check tests a different model from the one it is used to exonerate (correctness)

**Location.** Section 5.2, lines 1082 to 1106. Table tab-kr.
Conclusions, point 4.

**Problems.**

1. The `mmrm` `sp_exp(t | ptID)` fit has no random intercept. It is
   the M2 rung, not the production M3. The manuscript calls it "the
   marginal counterpart of the implicated autoregressive layer" and
   draws conclusions about M3.
2. The `sp_exp` $\kappa$ values (E1 1.050, E2 1.137, E3 1.040) agree
   with my analytic M2 values (1.064, 1.132, 1.054). They do not agree
   with the ladder's nlme `gls` M2 rung (1.023, 1.070). On the same
   first 2,000 datasets, the mean `gls` SE for E2 is 1.347 against
   1.373 for `mmrm` (VERIFIED). The two fits of nominally the same
   model differ by 2 percent in mean SE. That points to an
   optimization or convergence problem in the ladder's `gls` rung
   (SUSPECTED).
3. The manuscript says that Kenward-Roger "is not implemented for the
   production corCAR1 mixed model in standard software". This is not
   correct. SAS PROC MIXED supports `DDFM=KR` with a random intercept
   and `REPEATED / TYPE=SP(POW)(t)` (INSPECTED, general knowledge;
   not run here).
4. Kenward-Roger on the unstructured covariance still rejects at
   0.0605 and 0.0620 for E1 and E3, which is 2.2 and 2.5 MC SE above
   0.05. The text says "toward 0.059" and does not flag the residual
   anti-conservatism (VERIFIED).
5. The Limitations call fitting an unstructured analysis a "natural
   extension" (lines 1221 to 1223). The paper already fits one, in
   Table tab-kr.

**Why it matters.** Paper 08 makes the more general and correct
argument: a df correction cannot fix a conservative test, because KR
inflates the variance and Satterthwaite lowers the df. Paper 10's
empirical leg does not test M3 at all.

**Remedy.**

- Relabel the KR rows as M2.
- Either run KR on M3 in SAS, or state the logical argument of
  paper 08 and drop the claim about M3.
- Resolve the `gls` versus `mmrm` discrepancy before using the M2 rung.
- Delete the contradictory limitation.

**Evidence.** VERIFIED for points 2 and 4. INSPECTED for points 1, 3
and 5.

### 4.7 Comparisons across procedures lack paired MCSEs and formal size tests (reporting)

**Location.** Section 4 throughout. Section 4.5.

**Problem.** The analyses share datasets: ladder rungs, model versus
CR2, and E1 to E3. Paired comparisons are therefore available, but
none is reported. MCSEs are quoted only generically (line 545), and
STYLE.md Section 1 requires them with every number. My paired values
(VERIFIED):

- **CR2 minus model-based size:** +0.0128 (paired MC SE 0.0027),
  +0.0161 (0.0029) and +0.0104 (0.0028) for E1 to E3.
- **CR2 against nominal:**
  - E1: 0.0483 (z = -0.42)
  - E2: 0.0450 (z = -1.26)
  - E3: 0.0491 (z = -0.24)
- **M1 minus M3 size, E1:** +0.0165 (0.0028).
- **Alternative cell, E2 minus E1 power:** +0.083 (0.0060)
  model-based and +0.081 (0.0058) CR2. The difference-in-differences
  is -0.0024 (0.0059).
- **Alternative cell, E2 minus E3:** +0.086 and +0.089.

The text quotes +0.086 and +0.082, and +0.088 and +0.089. These come
from unpaired marginal summaries with different numbers of usable
replicates (2,966, 2,981 and 2,972).

- **Mancl-DeRouen GEE, E2:** 0.058, which is 2.6 MC SE above 0.05
  (binomial p = 0.011). The text describes this as "approximate
  calibration".

**Remedy.**

- Add a size table with MCSE, binomial p-value and 2-MCSE flag.
- Report paired differences for every between-procedure claim.
- Compute the ranking gaps on the common replicate set.

**Evidence.** VERIFIED.

### 4.8 The GEE comparison omits its efficiency cost and mixes reference distributions (correctness)

**Location.** Section 5.4, lines 1157 to 1211.

**Problems.**

- The GEE empirical SD of $\hat\theta$ is 1.045 against 0.895 for the
  M3 `lme` on the same data, for E1 (VERIFIED). That is a relative
  efficiency of about 0.73. The claim that the working correlation
  "governs only the efficiency" understates this, and "the choice
  between them is practical rather than principled" is untenable when
  27 percent of the information is lost.
- The `geepack` `ar1` working correlation is indexed by visit order.
  It ignores the irregular Hybrid spacing (weeks 4, 8, 9, 10, 11, 12,
  16, 20) (INSPECTED).
- The naive test uses a $z$ reference. The Mancl-DeRouen test uses
  $t_{65}$. The SE correction and the reference change are therefore
  confounded (INSPECTED).
- The MD standard error comes from a separate `geesmv` refit. The
  estimate comes from `geeglm` (INSPECTED).

**Remedy.**

- Report efficiency.
- Use a calendar-time working correlation, or state the limitation.
- Hold the reference fixed when comparing SE estimators.
- Take estimate and SE from one fit.

**Evidence.** As labeled per item above.

### 4.9 The analysis model is not the project's recommended analysis, and the relation to G9 and paper 08 is not stated (scope, framing)

**Location.** Sections 3.1 and 4.4. Discussion.

**Problem.** Paper 02 recommends G9. In the notation of paper 10, G9
is E2 with the half-life chosen by AIC over {0.25, 0.5, 1, 2} weeks,
a REML refit, and CR2 Satterthwaite
(`02-carryover-sensitivity/whitepaper-g9.Rmd`). The CR2 remedy in
paper 10 is G8 (E2 with CR2 at the true half-life) or G6 (E1 with
CR2). The paper does not say this. It does not evaluate the selection
step, and it does not cite paper 02's evidence that selection adds a
measurable size increment. `docs/56` and `docs/57` recommend F1 with
expectancy terms, which bears directly on Major comment 1. Paper 08
reports $\kappa$ of 1.10 to 1.16 for the same model class and makes
the df-correction argument. Neither is cited.

The analysis model also uses the uncentered `bm`, a binary or
exposure-weighted drug term and no expectancy term. This is the
historical E-family model (`docs/43`), not the current recommendation.

**Remedy.**

- Position paper 10 explicitly as the diagnosis behind G6, G8 and G9.
- Cite paper 02 and paper 08.
- Add at least one null cell analyzed with G9 and with F1.

**Evidence.** INSPECTED.

### 4.10 The generality claims exceed the evidence (exposition)

**Location.** Abstract "Conclusions". Section 5.1. Conclusions, sixth
point. Title.

**Problem.** The paper studies one DGP, one cell, one estimand and one
design. It then recommends routine practice for "any mixed-model
interaction test". The title promises a "staged diagnosis of
working-covariance misspecification". By Major comment 1, the
diagnosis in the worked example is mostly of mean misspecification.
The claim that the bias may run "in either direction" is supported
only by the unstructured fit and the OLS and GEE rows. Those are other
models, not the production model.

**Remedy.**

- Retitle.
- Confine the claims to the example.
- Present generality as a hypothesis, supported by the
  $V + r r^\top$ argument.

**Evidence.** INSPECTED.

## 5. Minor comments

1. **Misquoted range.** Line 629 says the observed Type I error is
   0.029 to 0.033. Table tab-decomp gives 0.029 to 0.035 (E3 0.0354).
   The closed-form predictions (0.036, 0.031, 0.037) all exceed the
   observed values (0.033, 0.029, 0.035) by 1 to 1.2 MC SE. They err
   in the same direction because the datasets are shared. VERIFIED.
2. **Fallback contamination.** `diagnose-null-type1.R` refits without
   `corCAR1` on error and still marks the fit as converged. 32, 40 and
   40 replicates of E1, E2 and E3 in Table tab-decomp are therefore M1
   fits counted as M3. The effect is small, since $\kappa$ moves by
   about 0.001, but the protocol is wrong. The ladder correctly drops
   them. VERIFIED.
3. **Boundary estimates.** The random-intercept SD collapses below 1
   in 8.7 percent of M3 replicates, against a pseudo-true value of
   9.8. These fits sit at the boundary and are not discussed.
   VERIFIED.
4. **CR2 failures.** CR2 failed on 19 to 24 of 3,000 replicates per
   specification. These replicates were dropped silently. Report them.
   VERIFIED.
5. **CR2 degrees of freedom.** The Satterthwaite df of about 24 at
   $N = 70$ (5th to 95th percentile 17 to 30) deserves a sentence on
   which design features drive it. Paths of unequal on-drug fractions
   are the likely cause. VERIFIED (values).
6. **"Within-cluster contrast" is approximate.** $b_i X_{it}$ has a
   between-participant component, because the mean on-drug fraction
   differs by path (5/8 against 4/8). This qualifies the OLS argument
   (lines 729 to 745). INSPECTED.
7. **Retired label.** Line 951 says "Architecture-B". This label is
   retired (NOTATION.md Part 3), and the adjacent caption already says
   "covariance-moderation architecture". INSPECTED.
8. **Notation departures from NOTATION.md.** The linter is clean
   because these are not mechanized. INSPECTED.
   - $S_{x,it}$ should be $Y_{it}$.
   - $\mathrm{bm}_i$ should be $b_i$ or $B^c_i$.
   - $b_i$ and $\sigma_b^2$ for the random intercept collide with the
     canonical $b_i$ (standardized biomarker). Use $u_i$ and
     $\sigma_u^2$.
   - $X_{it}$ is the exposure and $X$ is also the design matrix in
     Section 2.2. Use $D_{it}$ or $D_{bc,it}$.
   - $\beta_2 t$ should be $\beta_T w_t$.
   - $\phi$ for the normal density should be $\varphi$, since $\phi$
     is reserved for the decay family.
   - $\mathrm{N}(0,1)$ should be $\mathcal N(0,1)$.
   - $\mathrm{E}$ should be $\mathbb E$.
   - The eigenvalues $\lambda_g$ should be $\ell_g$, since $\lambda$
     is the decay rate.
   - The working CAR(1) parameter should be named $\rho_a$.
   - The sign convention (rule 5) is not stated.
9. **STYLE.md compliance.** INSPECTED/VERIFIED.
   - No Glossary or Abbreviations section, although Section 4 of
     STYLE.md requires one, including for $\kappa$, M0 to M3, E1 to
     E3, Hybrid and corCAR1.
   - Three abstract sentences exceed 30 words. The longest has 57.
   - About 35 prose lines join clauses with semicolons.
   - Rhetorical questions appear at lines 184 to 188.
   - Metaphors: "spends precision to buy validity ... that price would
     rise" (lines 1231 to 1233).
   - Hyperbole: "by far the most common" (line 149).
   - "Endeavor" appears twice.
   - British forms: "decentred" (318), "idealised" (734), "studentised"
     (375). The title page has "organisation" and the hyphenated
     "Type-I" and "Mis-specification".
   - The STYLE grep has one hit, "graded" (line 436). It is used in
     its ordinary sense and is harmless, but rephrase it to keep the
     check clean.
10. **ADEMP.** "The reporting follows the ADEMP structure" (line 487)
    is asserted, but no ADEMP table appears. Add one.
11. **Build instructions.** The build comment (lines 62 to 63) tells
    the reader to call `rmarkdown::render()` directly. That
    contradicts the project rule to render through `tools/render.sh`.
    `in_header: claudecode.tex` is an inert leftover. INSPECTED.
12. **Title page.** The title page carries a placeholder
    corresponding-author identity and a placeholder repository URL.
    Its title differs in hyphenation from the Rmd title. INSPECTED.
13. **Linearization.** Equation (eq:deriv) gives
    $\alpha \approx 0.05 - 0.23(\kappa - 1)$. At $\kappa = 1.097$ this
    gives 0.028, against an exact closed-form value of 0.031. Say that
    the linearization is used only for the elasticity. VERIFIED (the
    derivative of -0.229 and the elasticity of -4.58 are correct).
14. **Plan mismatch.** The examination plan (Section 6) states that M3
    is "the well-specified member of the ladder". The manuscript
    states the opposite. The plan's stage-4b bootstrap and stage-5
    design-dependence analyses were never run. Say so in Limitations.
    INSPECTED.
15. **Positive definiteness.** The positive-definiteness statements
    are correct. `pd_corrected` is FALSE on all four paths. Rebuilt
    under the current package, which refuses repair, every stored cell
    builds:
    - $c_{bm}$ of 0 and 0.45;
    - `c.cfct` of 0.1 and 0.2;
    - minimum eigenvalue 0.10 to 0.18.

    No stored run depended on a repaired matrix. VERIFIED.
16. **Design.** The simulated design is the four-path Hendrickson
    Hybrid: weeks 4, 8, 9, 10, 11, 12, 16, 20, with on-drug patterns
    11110010, 11110001, 11100010 and 11100001. $N = 70$ is a total
    allocated 18, 18, 17, 17. The sibling-paper mislabeling problem
    does not arise here. VERIFIED.

## 6. Missing and questionable references

| Reference | Status | Note |
|---|---|---|
| Gurka, Edwards, Muller 2011, *Stat Med* 30:2696 ([doi](https://doi.org/10.1002/sim.4293)) | Missing, essential | Proves that covariance underspecification biases fixed-effect inference in large samples. Anticipates the "asymptotic" claim. |
| Lu, Mehrotra 2010, *Stat Med* 29:474 ([doi](https://doi.org/10.1002/sim.3820)) | Missing | Covariance misspecification and the sandwich in longitudinal RCTs. |
| Gosho et al. 2017, *SMMR* 26:2389 ([doi](https://doi.org/10.1177/0962280215597938)) | Missing | Mancl-DeRouen and Kauermann-Carroll against KR with UN in MMRM. Directly comparable to Tables tab-gee and tab-kr. |
| Luke 2017, *Behav Res Methods* 49:1494 ([doi](https://doi.org/10.3758/s13428-016-0809-y)) | Missing | Type I error of KR and Satterthwaite in LMMs. Requested by the editor's checklist. |
| Wang, Kong, Li, Zhang 2016, *Stat Med* 35:1706 ([doi](https://doi.org/10.1002/sim.6817)) | Missing | Source of `geesmv`, which the GEE driver uses. Software must be cited. |
| Li, Redden 2015, *Stat Med* 34:281 ([doi](https://doi.org/10.1002/sim.6344)) | Missing | Bias-corrected sandwich with a $t$ reference. |
| Zucker, Ruthazer, Schmid 2010, *J Clin Epidemiol* 63:1312 ([doi](https://doi.org/10.1016/j.jclinepi.2010.04.020)); Araujo, Julious, Senn 2016, *PLoS One* 11:e0167167 ([doi](https://doi.org/10.1371/journal.pone.0167167)) | Missing | Analysis of combined N-of-1 trials, including covariance choice. |
| Fai, Cornelius 1996; Schaalje, McBride, Fellingham 2002 | Missing (not located in PubMed; verify) | Satterthwaite df for mixed models, and KR under covariance misspecification. |
| Companion papers 02 and 08 (`@pmsimstats-paper02`, `@pmsimstats-paper08`) | Missing | G6, G8, G9 and the $\kappa$ findings of paper 08. |
| `clubSandwich` (Pustejovsky), `nlme`, `geepack` software citations | Missing | Only `mmrm` is cited. Versions installed here: clubSandwich 0.7.0, nlme 3.1.169, geesmv 1.3. |
| verbeke1997 | Questionable support | It concerns a misspecified random-effects *distribution*, not covariance misspecification. Line 1037 cites it for the latter. |
| sturdevant2021carryover | Uncited, bad metadata | It is in the bib but never cited. The first name is "S. Gwynn" in PubMed ([doi](https://doi.org/10.1016/j.conctc.2021.100711)), not "Sarah G." Remove it or correct it. |
| hendrickson2020 | Verified | PubMed metadata match ([doi](https://doi.org/10.3389/fdgth.2020.00013)). |
| mancl2001, kenward1997, panwall2002, fay2001, laird1982 | Verified | PubMed matches (PMIDs 11252587, 9333350, 12185894, 11764261, 7168798). |
| white1982, kauermann2001, verbeke1997, neuhaus1998, scott1982, watson1955, carter2017, fay2001, panwall2002, zeger1988, neuhaus1991, patterson1971, harville1977 | Bib hygiene | DOIs are missing, and the bib comment says "to be added after verification". Add them. |
| pmsimstats_carryover2026, pmsimstats_dgp_architectures2026 | Bib hygiene | The titles do not match the current titles of papers 02 and 01. Use the project-wide keys. |
| mmrm2026 | Minor | Version 0.3.17 matches the installed version. The year "2026" should be the release year of that version. |

Citation-to-claim support for the classical sources is broadly
appropriate. The exception is the Box (1954) papers, which concern
ANOVA F tests. They are cited as the origin of the scale channel in a
mixed-model Wald test. That is a loose analogy, not a source.

## 7. Suggestions for strengthening

1. Make the analytic pseudo-true $\kappa$ the centerpiece. It is
   exact, cheap, separates the mean and covariance contributions, and
   replaces the underpowered gradient.
2. Build a two-factor ladder. One factor is the mean model: as fitted,
   F1 with expectancy terms, saturated visit by path. The other is the
   covariance: M0 to M3, unstructured. Report population $\kappa$ and
   simulated size for every combination.
3. Rerun the null and alternative cells under the current refined DGP:
   `c.cfct = 0.2`, anchored carryover, stop on non-PD. Add G9 and F1
   as procedures.
4. Report the parametric-bootstrap reference promised in plan stage
   4b on a subsample. It is the natural gold standard against which
   CR2 is judged.
5. Add at least one cell with $N = 35$ to show the price of CR2 at
   small cluster counts. A cell with heavier carryover or a different
   design (OL+BDC) would test the ranking-invariance claim beyond one
   cell.
6. Consider reframing the paper as a short communication. Its message
   would be that in mixed-model interaction tests a misspecified mean
   absorbed by REML can make the model-based test conservative, and
   that CR2 repairs the size while a richer mean model repairs the
   cause. That message is novel and practically useful.

## 8. Scope of this review

**Done.**

- Read the full manuscript, including code chunks, the bib, the title
  page, the examination plan, STYLE.md, NOTATION.md, the seven
  drivers, the relevant parts of `simulation-core.R` and the tidyverse
  `functions.R`, `docs/61`, paper 02's `whitepaper-g9.Rmd` and paper
  08's conservatism section.
- Reproduced every table value and most prose numbers from the stored
  `.rds` files.
- Computed binomial tests and paired MCSEs.
- Rebuilt every stored cell's covariance under the current package.
- Computed population $\kappa$ analytically for four DGP variants and
  the full ladder.
- Ran one paired simulation of 1,487 replicates at about two minutes
  of compute.
- Ran the STYLE grep and `tools/notation-lint.pl`, which was clean.
- Checked 13 citations against PubMed.

**Not done.**

- I did not rerun the paper's drivers, the full ladder, CR2, GEE or KR
  under the corrected DGP. The magnitudes in Major comment 2 for the
  refined DGP are analytic population values, not simulated sizes.
- I did not run KR on M3 in SAS.
- I did not render the manuscript.
- I did not review `bullets.Rmd`.
- I did not verify pre-1960 page numbers or the econometrics citations
  outside PubMed.
- I did not check `check-null-architecture-equivalence.R` beyond
  reading the plan's description of it.
- The `gls` versus `mmrm` M2 discrepancy is SUSPECTED and its cause is
  undiagnosed.
- The analytic $\kappa$ treats the biomarker as independent and
  mean-zero at the null, and uses equal path allocation. Finite-sample
  terms are excluded by construction.
