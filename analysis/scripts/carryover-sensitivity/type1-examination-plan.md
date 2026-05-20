# Examination plan: type-I conservatism of the interaction test
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

---
*Rendered on 2026-05-20 at 11:30 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/scripts/carryover-sensitivity/type1-examination-plan.md*
