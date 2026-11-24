---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# Normalizing the Biomarker in Papers 01 and 02 {.unlisted .unnumbered}
*2026-11-23 09:32 PST*

**Author.** pmsimstats team

**Purpose.** Papers 01 and 02 put the biomarker into the analysis model
on different scales. Some analyses use it raw (mean 124.3, SD 15.4 in
parameter set 1). Others center it or standardize it. This note asks
whether that matters, and sets out the case for and against
normalizing the biomarker throughout.

**Status.** Exploratory. The checks in Sections 2 and 3 use synthetic
data, not the papers' data-generating process (DGP). No manuscript
text or result has been changed.

## Summary

1. For every analysis whose model contains the biomarker, the drug
   term and their interaction, the scale of the biomarker does not
   change the test. Size, power and the AIC choice of half-life are
   the same raw or standardized. Only the units of the coefficient
   change. We confirmed this numerically to within $10^{-9}$.
2. One reported analysis is not of this kind. The free-baseline model
   (ES) of paper 01's stress test has a biomarker-by-baseline term
   without a baseline main effect. Its result depends on where the
   biomarker is centered. The stress test centers it at the sample
   baseline mean, which fits the intended model. A raw biomarker
   would fit a different one.
3. The two papers report coefficients on different scales. Paper 01
   reports per raw biomarker unit; paper 02 reports per biomarker SD.
   They differ by a factor of $\sigma_{bm} = 15.4$, and neither paper
   says so.
4. We recommend no refits. We recommend one sentence in each paper
   stating the scale and the conversion, and a stated centering rule
   for ES.

## 1. What the code does now

The DGP draws the biomarker as one component of a multivariate normal
vector, so it is normal with the mean and SD of the parameter set:
$\mathcal{N}(124.3, 15.4^2)$ in set 1 (Hendrickson et al.'s values) and
$\mathcal{N}(0, 1)$ in set 2. What the analysis model receives depends on the
driver.

| Analysis | Script | Biomarker in the model |
|---|---|---|
| Reruns of Hendrickson et al.'s code (paper 01) | `quick-sim/hendrickson-problems/14-`, `15-`, `17-` | Raw |
| Package `lme_analysis()` | `R/lme_analysis.R` | Raw |
| A4 on the Figure 4B grid, Appendix B.4 rerun (paper 01) | `19-g9-planning-and-fig4b.R`, `26-co-dbc-fix-paired.R` | Standardized by the DGP mean and SD |
| All paper 02 analyses, docs/77 | `carryover-sensitivity/simulation-core.R` | Standardized by the DGP mean and SD |
| Stress test, study 2 (paper 01) | `11-strawman-stress-test.R` | Centered at the sample baseline mean, raw SD units |

The standardization in `simulation-core.R` (`prepare_long_data()`,
line 348) divides by the DGP's own mean and SD. These are known
constants in a simulation, so the map is a fixed linear one. In a real
trial they would be estimated (Section 5).

Paper 01 states its estimand per raw unit,
$\beta_{bm:D} = -c_{bm}\sigma_{BR}/\sigma_{bm}$. At $c_{bm} = 0.25$
this is $-0.25 \times 8 / 15.36 \approx -0.130$, and the reported mean
estimate of $-0.126$ matches it. From the standardized A4 runs paper
01 reports power only. Paper 02 states its target per biomarker SD,
$\theta = -c_{bm}\sigma_{BR}$, which is $-3.6$ at $c_{bm} = 0.45$.

## 2. When the scale does not matter

Suppose the model contains an intercept, the biomarker $B$, the drug
term $X$ and the product $B X$. Replace $B$ by $(B - a)/s$. The new
columns are linear combinations of the old ones, and the old columns
of the new. The two models therefore span the same mean space. They
have the same fitted values, the same likelihood and the same
residuals. The interaction coefficient is multiplied by $s$, and its
SE is multiplied by $s$ as well. In other words, the $t$ statistic is
unchanged.

The same argument covers the parts of the analysis that go beyond
least squares. The REML likelihood changes only by a constant, so the
variance estimates do not move. The CR2 adjustment depends on the hat
matrix, which depends only on the column space. The Satterthwaite
degrees of freedom follow. The AIC comparison of candidate half-lives
compares models that all shift by the same constant, so the chosen
half-life does not change.

We checked this numerically
(`analysis/scripts/quick-sim/biomarker-scaling/01-invariance-check.R`).
Twenty synthetic data sets of 70 participants with eight visits each
were fitted twice, raw and standardized. The ratio of the interaction
estimates was 15.4 in every data set. The largest differences in $t$,
degrees of freedom and $p$ were below $10^{-9}$. This held for the
published analysis (`lmer`, random intercept, Satterthwaite) and for
the paper 02 core (`nlme` with `corCAR1`, CR2 SEs).

The argument requires only that every term containing $B$ comes with
its lower-order terms. A0 to A4 in paper 02, A4 with the baseline row
in docs/77, the expectancy-term model F1 (with $B D_e$ and $D_e$), and
Hendrickson et al.'s analysis all satisfy it.

## 3. When the scale does matter

The free-baseline model of paper 01's stress test,
`Sx ~ bmc + t + X + bmc:X + bmc:bl0`, does not satisfy it. Here
`bl0` is 1 at the baseline visit and 0 elsewhere, and there is no
`bl0` main effect. Shifting the biomarker by $a$ adds $a \cdot$`bl0`
to the span of the `bmc:bl0` column. With the biomarker centered, the
term lets only the slope on the biomarker differ at baseline. With the
raw biomarker, it also moves the baseline mean, by an amount tied to
that slope. That is a different and unintended model.

We checked the size of the difference on synthetic data with a
baseline row
(`analysis/scripts/quick-sim/biomarker-scaling/02-free-baseline-check.R`,
300 replicates per setting).

| Setting | Rejection, raw | Rejection, centered | Mean estimate, raw | Mean estimate, centered |
|---|--:|--:|--:|--:|
| Null, baseline on trend | 0.067 | 0.067 | 0.0042 | 0.0053 |
| Null, baseline 5 points off trend | 0.060 | 0.067 | 0.0126 | 0.0053 |
| $\beta = -0.13$, baseline 5 points off trend | 0.947 | 0.923 | $-0.1174$ | $-0.1247$ |

Estimates are per raw biomarker unit. The rejection rates differ by
less than two Monte Carlo SEs (about 0.015 here), so we do not claim a
difference in size or power. The individual fits do differ. The
$p$-values of the raw and centered fits differed by up to 0.66 within
one data set. When the baseline lies off the linear trend, the raw fit
moves the mean estimate away from the true value ($-0.13$), and the
centered fit does not. Standardizing a centered biomarker changed
nothing, as Section 2 predicts.

The stress test as run centers at the sample baseline mean, so it fits
the intended model. The point is that the centering is part of the
model's definition and should be stated with it.

## 4. The case for normalizing

**A common scale across parameter sets and papers.** The target of
inference is $c_{bm}$, which has no units. Per biomarker SD, the
coefficient is $-c_{bm}\sigma_{BR}$ in both parameter sets. Per raw
unit, it is 15.4 times smaller in set 1 than in set 2 for the same
$c_{bm}$. Bias, empirical SE and coverage tables are then comparable
across sets only after standardization. Paper 02 reports both sets,
and this is presumably why its driver standardizes.

**A drug coefficient that means something.** With a raw biomarker, the
coefficient of $X$ is the drug effect at a biomarker value of 0. In
set 1 that is about eight SDs below the mean, far outside the data.
Centering makes it the drug effect at the average biomarker. Neither
paper reports this coefficient, so the gain is for readers and for
real-data use rather than for the current results.

**Numerical stability.** A raw biomarker near 124 makes the biomarker
and product columns nearly collinear with the intercept and the drug
term. In our check, `lmer` gave no scale warnings in 900 raw-scale
fits. We did not compare convergence rates, and we found no practical
problem here. It remains a reason to center in general use.

**Robustness to the model form.** Section 3 shows that a centered
biomarker keeps models such as ES meaning what they are meant to mean.
Centering at a stated value removes that source of ambiguity.

## 5. The case against normalizing

**Fidelity to the published analysis.** Paper 01 reruns Hendrickson et
al.'s code to reproduce their results. That code fits the raw
biomarker. Its test would not change on another scale, but its
reported estimates would. Keeping the raw scale keeps the rerun a
rerun of their analysis.

**Natural units for clinical readers.** A clinician reads an effect
per mmHg of systolic blood pressure more easily than per SD. A per-SD
effect also depends on the spread of the biomarker in the enrolled
population. An enrichment trial with a narrow biomarker range would
report a smaller SD and a different per-SD effect for the same biology.
The per-unit coefficient does not have that dependence.

**Estimated, not known, constants in real data.** The simulation
standardizes by the DGP's mean and SD, which are known. A real trial
would standardize by the sample mean and SD. The per-SD effect is then
an effect per sample SD, which is itself random. At $N = 70$ the
sample SD varies by roughly 8.5 percent ($1/\sqrt{2(N-1)}$). The test
is unaffected (Section 2). The SE of the per-population-SD effect,
however, would ignore that variation, and coverage would be slightly
below nominal. We have not checked the size of this effect.

**Centering becomes part of the model.** For models such as ES, the
choice of center (sample mean, population mean, baseline mean) changes
the result. Normalizing does not remove that choice. It only makes it
visible.

## 6. Recommendations

1. **No refits.** Every reported test is either scale-invariant
   (Section 2) or already centered as intended (Section 3).
2. **Paper 01.** State in Methods that the coefficient is per raw
   biomarker unit. Note that multiplying by $\sigma_{bm} = 15.4$ gives
   the effect per biomarker SD, the scale of paper 02. Keep the raw
   scale for the reruns of Hendrickson et al.'s code.
3. **Paper 02.** State in the estimand section that `bm` in the model
   is standardized by the DGP mean and SD. The model line
   `Sx ~ bm * X + t` currently reads as if it used the raw biomarker,
   which in set 1 has mean 124.3.
4. **Paper 01, ES.** State that the biomarker is centered at the
   sample baseline mean, and that the model depends on this choice
   because it has no baseline main effect.
5. **Real-data guidance.** In the vignette and in
   `analyze_trial_extended()`, center the biomarker at the sample mean.
   Report the interaction per natural unit, with the per-SD value
   alongside it. Say that the per-SD value uses the sample SD.

## 7. Open questions

- The checks use synthetic data with a simple design, not the papers'
  DGP. We expect the invariance in Section 2 to carry over exactly,
  since it is algebraic. The size of the ES difference in Section 3
  under the stress test's own DGP is not known.
- We have not rerun the stress test with a raw biomarker in ES.
- The coverage loss from standardizing by a sample SD (Section 5) has
  not been measured.
- We have not traced the scale of the 'own implementation' estimates
  quoted in paper 01 ($-0.126$ falling to $-0.026$ to $-0.032$). They
  match the per-raw-unit target, so we infer a raw or centered
  biomarker.
