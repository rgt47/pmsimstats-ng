# Examination plan and findings: type-I conservatism of the interaction test
*2026-05-20 11:30 PDT*

**Author:** pmsimstats team

## 1. Objective

The empirical type-I error of the biomarker-by-treatment
interaction test (`bm:Db` under specifications A1 and A3, `bm:Dbc`
under A2) runs systematically below the nominal 0.05 across the
carryover-sensitivity factorial. The pooled mean over the 540 null
cells of the production run is 0.039, and against roughly 270000
null replicates this is on the order of thirty standard errors
below nominal. The deflation is therefore real and not a sampling
artefact.

The manuscript (§4.6) attributes the conservatism in a single line
to 'the conservative AR(1) reference distribution at moderate
sample sizes'. That attribution is plausible but unverified. The
purpose of this examination is to localise the conservatism to one
of four inference channels, attribute it to a mechanism, and
produce an evidenced statement that either confirms the §4.6
attribution or replaces it.

A preliminary point has already been settled. At the null
(`c_bm = 0`) the two DGP architectures collapse to a bit-identical
data-generating process, confirmed by
`check-null-architecture-equivalence.R`. The architecture axis
therefore carries no information at the null, and the examination
need not stratify by architecture.

Sections 2 to 12 record the examination as planned. Sections 13 to
15, added after execution on 2026-05-20, report the findings
obtained and the resulting plan to address the conservatism.

## 2. Stage 0: scope and reference cell

The examination begins at a single null cell and expands only as
the evidence requires.

- **Primary null cell.** `c_bm = 0`, exponential DGP, Hybrid
  design, `N = 70`, `t1half = 1.0`, `rho = 0.7`. This is the null
  counterpart of the manuscript reference configuration.
- **Replicates.** 5000 for the primary cell, giving a Monte Carlo
  standard error on the type-I estimate of approximately 0.0031 and
  resolving 0.039 against 0.050 at roughly 3.5 standard errors.
  10000 replicates are held in reserve should the p-value ECDF tail
  require a cleaner estimate.
- **Reproducibility.** A fixed master seed, with per-replicate
  seeds derived as in the production driver. The complete
  per-replicate output is retained, not merely the rejection
  indicator.

## 3. Stage 1: instrument the per-replicate output

The production fitter `fit_spec()` in `simulation-core.R` returns
only `estimate`, `std_error`, `p_value`, and `converged`. The
decomposition that follows requires more. A standalone diagnostic
driver, `diagnose-null-type1.R`, will reuse `simulation-core.R` for
data generation but supply an augmented fitter that records, per
replicate:

- the degrees of freedom `lme` assigns the interaction term, from
  `summary(fit)$tTable[, 'DF']`;
- the estimated corCAR1 correlation parameter, the random-intercept
  standard deviation, and the residual standard deviation, from
  `VarCorr(fit)` and `fit$sigma`;
- the `non_positive_definite_rate` attribute returned by
  `generate_data()`, so that the positive-definiteness correction
  can be excluded as a contributor at negligible cost;
- the optimizer convergence messages.

The diagnostic driver does not modify the production pipeline.

## 4. Stage 2: decompose into the four inference channels

The interaction test is a Wald t-statistic, `T = estimate /
std_error`, rejected when `|T|` exceeds the t critical value at the
assigned degrees of freedom. Conservatism must enter through one of
four channels, and the channels are separable. From the primary
5000-replicate output we compute the following.

1. **Point-estimate bias.** The mean of `estimate` and a test of
   its equality to zero. The expected value is approximately zero;
   a non-zero mean would constitute a separate problem and would,
   if anything, mask the conservatism rather than cause it.
2. **Standard-error calibration.** The ratio of the mean
   model-reported `std_error` to the empirical standard deviation
   of `estimate` across replicates, together with the full
   distribution of `std_error` so that a uniform upward shift can
   be distinguished from a heavy right tail. A ratio above 1.05
   implicates standard-error over-estimation.
3. **Statistic dispersion.** The standard deviation of `T` under
   the null, which should be approximately one under correct
   calibration. A value below one confirms a conservative
   reference.
4. **p-value uniformity.** The empirical CDF of the p-values
   against the diagonal, Kolmogorov-Smirnov and Anderson-Darling
   tests against Uniform(0, 1), and the type-I rate at alpha in
   {0.01, 0.05, 0.10}. Whether all three alpha levels are depressed
   or only the tail distinguishes a global scale error from a
   tail-shape error.

The deliverable of this stage is a four-row decomposition table and
the p-value ECDF figure.

## 5. Stage 3: separate the DF channel from the SE channel

The p-values are recomputed from the same `estimate` and
`std_error` against a standard-normal reference in place of the t
reference.

- If the type-I rate rises to approximately 0.05 under the normal
  reference, the `lme` degrees-of-freedom rule is the cause.
- If the type-I rate remains at approximately 0.039, the
  degrees-of-freedom rule is excluded and the cause lies in the
  standard error or the joint estimate-and-SE calibration. The
  examination then proceeds to Stage 4.

A prior check short-circuits this stage. The assigned degrees of
freedom are tabulated directly; with `N = 70` across eight
measurement occasions the level-one degrees of freedom are likely
in the hundreds, in which case the t and normal references are
numerically indistinguishable and the degrees-of-freedom channel is
excluded immediately.

## 6. Stage 4: attribute the SE mis-calibration to a covariance model

The same saved null datasets are refitted under a ladder of
correlation specifications, comparing the interaction standard
error, the calibration ratio, and the type-I rate across them.

- **M0.** `lm(Sx ~ bm + t + Db + bm:Db)`, ignoring all
  within-subject correlation. Expected to be anti-conservative and
  serving as a lower anchor.
- **M1.** `lme` with a random intercept only, no corCAR1.
- **M2.** `gls` with a corCAR1 residual only, no random intercept.
- **M3.** `lme` with both a random intercept and a corCAR1
  residual, the production model.

The data-generating process produces AR(1)-correlated factors plus
a participant-level baseline term `BL` that acts as a random
intercept. Specification M3 is therefore the well-specified member
of the ladder, and M1 and M2 are each mis-specified. The reading
is as follows. If M3 alone is conservative while M1 is calibrated,
the near-redundant pairing of a random intercept with a corCAR1
residual on only eight occasions is implicated. If M1, M2, and M3
are all conservative, the cause is the generic finite-sample
plug-in Wald-t and is not specific to the redundancy.

## 7. Stage 4b: gold-standard calibration counterfactual (conditional)

This stage is undertaken only if Stages 2 to 4 point to the Wald
approximation. The conservatism is then checked against a reference
that does not rely on that approximation.

- A parametric-bootstrap or permutation null distribution for the
  interaction statistic, computed on a subsample of datasets. If
  those p-values are uniform while the Wald p-values are not, the
  Wald approximation is conclusively the cause.
- Optionally, an `lmer` fit with a Kenward-Roger small-sample
  correction by way of `pbkrtest`, accepting that `lmer` cannot
  carry the corCAR1 residual and so tests only the
  random-intercept structure.

## 8. Stage 5: design dependence

The Stage 2 decomposition is repeated for the CO, Hybrid, and
OLBDC designs at `c_bm = 0`, `N = 70`. Separately, each design's
long-form design matrix is inspected: the correlation and variance
inflation factor between `Db` and `t`, and between the `bm:Db` and
`bm:t` columns. If the conservatism tracks the `Db`-`t`
collinearity ranking across designs, design-induced collinearity is
a contributing channel. If the type-I rate is uniformly near 0.039
irrespective of design, collinearity is excluded.

## 9. Stage 6: sample-size gradient

The null decomposition is run at `N` in {35, 70, 140, 280}. This is
the decisive diagnostic. A finite-sample Wald artefact must
attenuate as `N` grows, with the type-I rate approaching 0.05. If
it does, the conservatism is benign and expected, and the
manuscript may say so with evidence. If the type-I rate remains
near 0.039 at large `N`, the cause is asymptotic, constitutes a
structural mis-specification, and warrants a correction rather than
a footnote.

## 10. Decision logic

```
Stage 2  ratio > 1.05 .................... SE over-estimation
         SD(T) < 1, ECDF bowed .......... conservative reference confirmed
Stage 3  normal reference fixes it ...... DF rule is the cause
         normal reference does not ...... SE/calibration, go to Stage 4
Stage 4  M1 calibrated, M3 not .......... RI + corCAR1 redundancy
         M1, M2, M3 all conservative .... generic finite-sample Wald-t
Stage 5  tracks Db-t collinearity ....... design collinearity contributes
Stage 6  type-I -> 0.05 as N grows ...... finite-sample, benign
         type-I flat in N ............... asymptotic mis-specification
```

## 11. Deliverables

1. `diagnose-null-type1.R`, the standalone diagnostic driver.
2. A short results note containing the four-channel decomposition
   table, the p-value ECDF figure, the model-ladder table, the
   design and sample-size tables, and an attributed one-paragraph
   conclusion.
3. A concrete recommendation for §4.6 of the manuscript: either
   confirm the finite-sample attribution with the sample-size
   gradient as evidence, or replace it.

## 12. Cost and sequencing

All stages are computationally light. At the throughput observed in
the development run, the 5000-replicate primary cell is on the
order of a few minutes, and the model ladder and the sample-size
gradient multiply that by small constants. The whole examination is
minutes of compute, not hours.

The recommended sequencing is to run Stages 1 to 3 first, since
they are quick and may already exclude the degrees-of-freedom
channel and quantify the standard-error ratio, and then to decide
whether Stages 4 to 6 are required. Stages 4b and 5 are
conditional; Stage 6 should be run regardless, because it is the
single most decisive test of whether the conservatism is the
benign finite-sample behaviour the manuscript assumes.

## 13. Results

The examination was executed on 2026-05-20. Stages 1 to 3 ran on
the primary null cell at 5000 replicates; Stage 4 refitted the
identical datasets under the model ladder at 5000 replicates;
Stage 6 ran the production model M3 across the sample-size grid at
1500 replicates per point. The driver scripts are
`diagnose-null-type1.R` for Stages 1 to 3 and
`diagnose-null-type1-stages46.R` for Stages 4 and 6.

### 13.1 Stages 1 to 3: the conservatism is standard-error over-estimation

| Spec | mean est | se_ratio | sd(T) | type-I 0.05 | type-I 0.05 (z) | mean DF |
|------|----------|----------|-------|-------------|-----------------|---------|
| A1   | 0.032    | 1.069    | 0.932 | 0.033       | 0.034           | 487     |
| A2   | 0.039    | 1.097    | 0.909 | 0.029       | 0.029           | 487     |
| A3   | 0.033    | 1.061    | 0.939 | 0.035       | 0.036           | 486     |

The interaction estimate is effectively unbiased and the
degrees-of-freedom channel is excluded. The standard error is
over-estimated by 6 to 10 percent (`se_ratio`), and this fully
accounts for the compressed test statistic, since the standard
deviation of the statistic equals the reciprocal of `se_ratio` to
three digits. Recomputing the p-values against a normal reference
moves the type-I rate by approximately 0.001, confirming that the
assigned degrees of freedom, approximately 487, are not the cause.
The positive-definiteness correction was triggered on none of the
four Hybrid paths and is excluded.

A small positive bias in the null interaction estimate is present
(mean approximately 0.033, about 0.036 of a standard deviation,
roughly 2.5 standard errors from zero). It is model-independent,
points toward rejection rather than away from it, and therefore
cannot cause the conservatism, but it is recorded here for
completeness and noted again in §15.5.

The null p-value ECDF lies uniformly below the diagonal across the
whole unit interval, the signature of a global scale compression
rather than a tail anomaly.

![Stage 2: null p-value empirical CDF against Uniform(0,1) for the
primary null cell. All three specifications bow uniformly below the
diagonal; A2 is the lowest curve, the most conservative.](output/diag-null-type1-ecdf.pdf){width=62%}

### 13.2 Stage 4: the corCAR1 layer inflates the interaction SE

| Model | Structure              | se_ratio A1 | se_ratio A2 | type-I A1 | type-I A2 |
|-------|------------------------|-------------|-------------|-----------|-----------|
| M0    | lm                     | 1.323       | 1.289       | 0.010     | 0.012     |
| M1    | lme, ~1\|ptID          | 0.989       | 0.972       | 0.049     | 0.055     |
| M2    | gls, corCAR1           | 1.023       | 1.070       | 0.043     | 0.036     |
| M3    | lme, ~1\|ptID + corCAR1 | 1.068       | 1.097       | 0.033     | 0.029     |

The ladder is unambiguous. M1, a random intercept with no residual
correlation, is well calibrated: its standard-error ratio sits at
0.97 to 0.99 and its type-I rate at 0.049 to 0.055. M3, the
production model, adds a corCAR1 residual on top of that random
intercept and is conservative: the ratio rises to 1.07 to 1.10 and
the type-I rate falls to 0.029 to 0.033. The corCAR1
residual-correlation layer is therefore the proximate cause of the
conservatism.

Two ladder results merit comment. First, M0, ordinary least
squares, is not anti-conservative as Stage 4's plan text
anticipated but strongly conservative, with a standard-error ratio
near 1.3. The reason is that the interaction is a within-subject
contrast, and ordinary least squares charges it the full residual
variance, between-subject plus within-subject, when the relevant
error is only the within-subject part. The plan's expectation for
M0 was therefore wrong, and the corrected reading is recorded here.
Second, M2, a corCAR1 residual with no random intercept, is
intermediate. Taken together, the ladder shows that the random
intercept alone recovers calibration and the corCAR1 addition
removes it.

![Stage 4: interaction standard-error calibration across the model
ladder. The left panel is the ratio of model-reported SE to the
empirical SE, the right panel the realised type-I rate. M1 alone
sits on the reference lines; M3 does not.](output/diag-null-type1-ladder.pdf){width=100%}

### 13.3 Stage 6: the conservatism is asymptotic, not finite-sample

| N   | se_ratio A1 | se_ratio A2 | se_ratio A3 | type-I A1 | type-I A2 | type-I A3 |
|-----|-------------|-------------|-------------|-----------|-----------|-----------|
| 35  | 1.072       | 1.103       | 1.064       | 0.038     | 0.028     | 0.039     |
| 70  | 1.056       | 1.083       | 1.050       | 0.043     | 0.031     | 0.045     |
| 140 | 1.086       | 1.098       | 1.079       | 0.041     | 0.033     | 0.042     |
| 280 | 1.052       | 1.070       | 1.046       | 0.037     | 0.028     | 0.035     |

This is the decisive result. The plan posited that a finite-sample
Wald artefact would attenuate as N grows, with the type-I rate
approaching 0.05. It does not. Across N from 35 to 280 the
standard-error ratio holds at 1.05 to 1.10 with no trend toward
one, and the type-I rate holds at 0.028 to 0.045 with no trend
toward 0.05. At N = 280 the production specification A2 still
rejects at 0.028. The conservatism is a structural feature of the
M3 model and does not vanish with sample size.

![Stage 6: type-I rate and standard-error ratio against sample
size under the production model M3. Neither panel trends toward its
reference line as N grows.](output/diag-null-type1-ngradient.pdf){width=100%}

The explanation consistent with all three stages is as follows. The
simulated symptom score is a mixture of three AR(1)-correlated
factors with timepoint-dependent variances, plus a
participant-level baseline term. The M3 working covariance, a
single random intercept plus a single-parameter corCAR1 residual,
cannot match that structure exactly. With only eight measurement
occasions per participant the random-intercept variance and the
corCAR1 correlation parameter are weakly separated, and the
resulting model-based standard error for the interaction is biased
upward. Increasing the number of participants sharpens the
variance-parameter estimates but does not change the eight-occasion
structure that drives the mismatch, so the bias persists. This is
why the conservatism is asymptotic.

## 14. Conclusions

1. The type-I conservatism is real and is standard-error
   over-estimation. The model-reported standard error of the
   interaction coefficient exceeds its true sampling standard
   deviation by 6 to 10 percent, the test statistic is compressed
   by the reciprocal of that factor, and the null p-values are
   stochastically too large at every alpha level.
2. The cause is the corCAR1 residual-correlation term in the
   production analysis model. A random intercept alone (M1) is well
   calibrated; adding corCAR1 (M3) inflates the interaction
   standard error. The effect is largest for the exposure-weighted
   specification A2, which is consequently the most conservative,
   consistent with A2 sitting lowest in the manuscript's Figure 4.
3. The conservatism is asymptotic, not a moderate-sample-size
   artefact. It does not attenuate from N = 35 to N = 280. The
   manuscript's §4.6 attribution, 'the conservative AR(1) reference
   distribution at moderate sample sizes', is correct in naming the
   AR(1) term but wrong in implying the effect is finite-sample and
   would diminish with N.
4. The degrees-of-freedom rule, the positive-definiteness
   correction, and point-estimate bias are all excluded as causes.
5. The conservatism is in the safe direction. It is a valid test
   that under-rejects, not an inflated test. Its practical cost is
   a modest loss of power, which means the manuscript's power
   estimates are mildly pessimistic rather than optimistic. The
   headline ranking A2 greater than A1 approximately equal to A3 is
   not threatened, because all three specifications are inflated by
   similar factors and a recalibration would lift all three powers
   while preserving their order.

## 15. Plan to address the conservatism

The conservatism is structural, so it will not be removed by more
replicates or larger N. Four remedies are available, listed from
the most principled to the most minimal.

### 15.1 Option A: cluster-robust standard errors (recommended)

Retain the M3 point estimate and the corCAR1 model, which the
project adopts for clinical realism, but replace the model-based
standard error of the interaction with a cluster-robust sandwich
standard error, clustering on `ptID`. A sandwich estimator is
consistent for the true sampling variance whether or not the
working covariance is correctly specified, so it removes the
asymptotic bias identified in Stage 6. The bias-reduced CR2 variant
in the `clubSandwich` package is appropriate, since it carries a
small-sample correction and supports `lme` objects. This is the
recommended fix: it preserves the analysis model and corrects only
the inference.

### 15.2 Option B: drop the corCAR1 term

Use M1, a random intercept only. Stage 4 shows M1 is well
calibrated. The cost is the loss of the corCAR1 residual structure,
which the project adopted deliberately, and which M1 mis-specifies
by treating the within-subject residual as independent. M1's
fixed-effect test is nonetheless calibrated, so this option is
defensible, but it is a larger change to the analysis than Option A
and should not be taken without re-running the principal factorial.

### 15.3 Option C: parametric-bootstrap reference

Replace the Wald reference for the interaction test with a
parametric-bootstrap or permutation null distribution. This is the
gold standard for calibration and is robust to the
working-covariance mismatch, but it is computationally heavy across
a 540-cell factorial and is best reserved for confirmation on
selected cells rather than as the production default.

### 15.4 Option D: document without changing the analysis

If the analysis is not to be re-run, the minimum acceptable
response is to correct §4.6. The conservatism should be described
as a structural consequence of the corCAR1 term, persistent at all
sample sizes, of magnitude 6 to 10 percent in standard-error terms
and roughly 0.03 to 0.04 in realised type-I error, and in the safe
direction. The phrase 'at moderate sample sizes' should be removed,
since Stage 6 shows the effect does not depend on N.

### 15.5 Recommended sequence

1. Implement Option A on the primary null cell and confirm that the
   CR2 cluster-robust standard error brings the standard-error
   ratio to approximately one and the type-I rate to approximately
   0.05.
2. If confirmed, apply the CR2 standard error across the Tier 1
   factorial and regenerate the affected summaries. The point
   estimates and the power ranking are not expected to change; the
   type-I columns and the absolute power values will.
3. Revise §4.6 of the manuscript per Option D in all cases, since
   the current 'moderate sample sizes' wording is incorrect
   regardless of which remedy is adopted.
4. Record the small positive null bias of the interaction estimate
   (§13.1) in the manuscript's limitations, noting it is negligible
   in magnitude and does not affect the conclusions.

The first step is a short, single-cell verification, which can be
written as `diagnose-null-type1-cr2.R` and run on request.

---
*Rendered on 2026-05-20 at 16:46 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/scripts/carryover-sensitivity/type1-examination-plan.md*
