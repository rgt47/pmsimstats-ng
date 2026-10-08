# Baseline as a Covariate in the Analysis of the Hybrid N-of-1 Design: Advantages and Disadvantages {.unlisted .unnumbered}
*2026-10-03 17:17 PDT*

**Author.** pmsimstats team

**Purpose.** The published analysis of Hendrickson et al. [1], and every
analysis in this project so far, treats the baseline symptom score as
an observation of the response, a ninth row coded off drug. This paper
weighs the alternative used in most confirmatory trials: baseline as a
covariate in a model for the post-baseline visits. It also covers a
third option, constrained longitudinal data analysis, which keeps
baseline in the response but changes how it is modeled. Section 4 gives
a formal account of what adjusting for baseline does to each quantity
the analysis estimates, with each result checked numerically against
the published construct. No new Monte Carlo simulation was run;
Section 9 lists what the paused covariance study
(`08-covariance-study.R`) would settle.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md` and `docs/40`. Symbols
marked *local* are defined here.

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $Y_{it}$ | symptom score of participant $i$ at visit $t$ (code `Sx`); $t = 0$ is baseline | canonical |
| $\mathrm{BL}_i$ | baseline symptom score, $Y_{i0} = \mathrm{BL}_i$ | canonical |
| $BR_{it}$, $PB_{it}$, $TV_{it}$ | drug, expectancy and natural-course responses | canonical |
| $S_{it}$ | sum of the three responses, $TV_{it} + PB_{it} + BR_{it}$ | local |
| $B_i$ | biomarker | canonical |
| $D_{it}$ | binary drug state (code `Db`); 0 at baseline | canonical |
| $\beta_B$ | biomarker main effect | local |
| $\beta_D$, $\beta_{bm:D}$ | drug main effect; biomarker-by-drug interaction | canonical |
| $\mathbf{Y}_{ip}$ | vector of the 8 post-baseline outcomes | local |
| $\sigma_{00}$, $\boldsymbol\sigma_{p0}$, $\Sigma_{pp}$ | baseline variance, baseline-by-visit covariances, post-baseline covariance | local |
| $\boldsymbol\gamma$ | baseline regression slopes, $\boldsymbol\gamma = \boldsymbol\sigma_{p0} / \sigma_{00}$ | local |
| $\mathbf{a}$ | contrast weights on the post-baseline visits, $\sum_t a_t = 0$ | local |
| $n_{\text{on}}$, $n_{\text{off}}$ | numbers of on-drug and off-drug post-baseline visits in a path | local |

Table: Notation index and glossary of symbols used in this document

- **Response coding.** Baseline enters as a row of the outcome, at
  $t = 0$ with $D = 0$ (the published analysis).
- **Covariate coding (ANCOVA).** Only post-baseline visits are
  outcomes; $\mathrm{BL}_i$ is a fixed covariate.
- **cLDA.** Constrained longitudinal data analysis [5, 6]: baseline
  stays in the response vector, its mean is constrained to be common
  across randomization groups, and the covariance of all visits is
  modeled, usually unstructured.
- **E9 and E9b.** The paired-difference statistics of `docs/37`: each
  participant's on-drug mean minus off-drug mean, regressed on the
  biomarker. E9 uses post-baseline visits; E9b also counts baseline as
  an off-drug visit, which mirrors the response coding.
- **Unstructured covariance (`us`).** Every visit has its own variance
  and every pair its own covariance; the MMRM default.

## 1. Summary

1. **Baseline adjustment changes the analysis, not the data.** The same
   replicates can be fitted both ways (Section 2).
2. **In this construct the baseline slopes are all exactly 1.**
   Post-baseline scores are $\mathrm{BL}_i - S_{it}$, and baseline is
   uncorrelated with every response component and with the biomarker.
   Conditioning on baseline is therefore change from baseline
   (Section 3, Proposition 1).
3. **A covariate gives a within-participant contrast nothing when the
   baseline slopes are equal.** The variance reduction from conditioning
   a contrast $\mathbf{a}'\mathbf{Y}_{ip}$ on baseline is
   $\sigma_{00}(\mathbf{a}'\boldsymbol\gamma)^2$. That is zero whenever
   every visit has the same slope, because the weights sum to zero. The
   interaction is such a contrast. The construct's E9 variance is 44.13
   with or without conditioning (Proposition 2).
4. **Between participants the gain is large.** Conditioning reduces the
   variance of a blinded post-baseline value from 599 to 257. That
   benefits the biomarker main effect and comparisons of path means,
   which are not the target (Proposition 3).
5. **What matters for the interaction is whether the baseline row is
   counted as an off-drug visit.** Counting it leaves part of each
   participant's shared response level in the contrast. That raises the
   contrast variance by 15% to 16%. It also creates a residual
   interaction of $1/(n_{\text{off}} + 1)$ of the full value under the
   published step rule with carryover: the power floor of the published
   analysis (Proposition 4).
6. **The baseline row is incompatible with a random-intercept
   covariance.** It has no noise of its own, which compound symmetry
   cannot represent. Dropping it, or modeling an unstructured
   covariance, removes that misspecification of the model-based
   standard errors (Proposition 5).
7. **Recommendation.** Analyze the post-baseline visits with baseline as
   a covariate in an MMRM with an unstructured covariance, or use cLDA,
   which coincides with it for complete data. Avoid the published
   combination of a baseline row, a random intercept and linear time.
   The covariate's value for the interaction is robustness, not
   efficiency: it pays when baseline predicts on-drug and off-drug
   visits differently, which the construct rules out but a real trial
   may not (Sections 4 and 8).

## 2. The two codings

**Response coding (published).** For each participant the outcome
vector has nine entries, $(\mathrm{BL}_i, Y_{i1}, \ldots, Y_{i8})$.
The baseline row has $t = 0$, $D = 0$ and expectancy 0 (`lme_analysis.R`
at commit `58b32a9`, lines 62 to 87). The published model,
`lmer(Sx ~ bm + Db + t + bm * Db + (1 | ptID), data = long)`, treats it
as one more off-drug visit.

**Covariate coding.** The outcome vector has the eight post-baseline
entries, and baseline is a covariate:

$$
Y_{it} = \beta_0 + \gamma_{\mathrm{BL}}\,\mathrm{BL}_i + \beta_B B_i
+ \beta_D D_{it} + f(t) + \beta_{bm:D}\,B_i D_{it} + \varepsilon_{it},
\qquad t = 1, \ldots, 8,
$$

with $f(t)$ the time adjustment and $\varepsilon_i$ given a working
covariance, for example `us(visit | ptID)` in `mmrm`.

The simulated data are unchanged: the same replicates, with the
baseline value moved from a row of the response to a column of the
design matrix.

## 3. The data-generating process

In the published construct $Y_{it} = \mathrm{BL}_i - S_{it}$ for
$t \geq 1$, and $\mathrm{BL}_i$ is drawn independently of the
biomarker and of all response components. Hence

$$
\sigma_{t0} = \operatorname{Cov}(Y_{it}, \mathrm{BL}_i) = \operatorname{Var}(\mathrm{BL}_i) = \sigma_{00}
\quad \text{for every } t \geq 1,
\qquad
\Sigma_{pp} = \sigma_{00} J + \operatorname{Cov}(\mathbf{S}_{ip}),
$$

with $J$ the matrix of ones. Computed exactly from the construct (path
1, no interaction; `docs/40`, Section 6):

| | Variance | Covariance within phase | Covariance with baseline |
|---|---|---|---|
| Baseline | 342 | | |
| Open label (weeks 4, 8) | 710 | 605 | 342 |
| Blinded (weeks 9 to 20) | 599 | 527 | 342 |

Table: Exact variances and covariances by phase, path 1, no interaction

with 556 between open-label and blinded visits. Two features matter
below. Baseline has no noise of its own: its covariance with every
later visit equals its own variance. And every later visit has the same
covariance with baseline.

## 4. A formal assessment of adjusting for baseline

Throughout, the outcomes are multivariate normal, the mean model is
correctly specified, and the data are complete. Write
$\mathbf{Y}_i = (Y_{i0}, \mathbf{Y}_{ip}')'$ with

$$
E(Y_{i0}) = \mu_0 + \beta_B B_i, \qquad
E(\mathbf{Y}_{ip}) = X_{ip}\boldsymbol\beta, \qquad
\Sigma = \begin{pmatrix} \sigma_{00} & \boldsymbol\sigma_{p0}' \\ \boldsymbol\sigma_{p0} & \Sigma_{pp} \end{pmatrix}.
$$

Baseline is off drug, so $D_{i0} = 0$ and neither $\beta_D$ nor
$\beta_{bm:D}$ enters its mean.

### Proposition 1: the covariate model is the conditional model

Conditional on baseline,

$$
\mathbf{Y}_{ip} \mid Y_{i0} \;\sim\; N\!\Bigl( X_{ip}\boldsymbol\beta + \boldsymbol\gamma\,(Y_{i0} - \mu_0 - \beta_B B_i),\;\; \Sigma_{pp\cdot 0} \Bigr),
\qquad
\boldsymbol\gamma = \frac{\boldsymbol\sigma_{p0}}{\sigma_{00}}, \quad
\Sigma_{pp\cdot 0} = \Sigma_{pp} - \frac{\boldsymbol\sigma_{p0}\boldsymbol\sigma_{p0}'}{\sigma_{00}}.
$$

*Proof.* The standard conditional distribution of a partitioned
multivariate normal. $\square$

ANCOVA with visit-specific baseline slopes $\gamma_t$ and residual
covariance $\Sigma_{pp\cdot 0}$ is exactly this model. A common slope
$\gamma_{\mathrm{BL}}$ is the special case $\gamma_t \equiv
\gamma_{\mathrm{BL}}$. Conditioning changes the intercept and the
biomarker main effect, which absorbs $-\gamma_t \beta_B$. It leaves
$\beta_D$ and $\beta_{bm:D}$ unchanged, because they do not appear in the
baseline mean. The estimand of the interaction is the same under both
codings.

**In the construct,** $\sigma_{t0} = \sigma_{00}$ for every $t$, so
$\boldsymbol\gamma = \mathbf{1}$ (verified: every slope is 1 to six
decimals) and $\Sigma_{pp\cdot 0} = \operatorname{Cov}(\mathbf{S}_{ip})$:
conditioning on baseline removes exactly its variance from every entry.
The covariate analysis then coincides with a change-from-baseline
analysis, and the classic argument between the two [2, 3] does not
arise.

### Proposition 2: a covariate adds nothing to a within-participant contrast when the slopes are equal

Let $\mathbf{a}$ be contrast weights on the post-baseline visits with
$\sum_t a_t = 0$; the E9 weights $1/n_{\text{on}}$ and
$-1/n_{\text{off}}$ are an example. Then

$$
\operatorname{Var}(\mathbf{a}'\mathbf{Y}_{ip} \mid Y_{i0})
= \mathbf{a}'\Sigma_{pp}\mathbf{a} - \frac{(\mathbf{a}'\boldsymbol\sigma_{p0})^2}{\sigma_{00}}
= \mathbf{a}'\Sigma_{pp}\mathbf{a} - \sigma_{00}\,(\mathbf{a}'\boldsymbol\gamma)^2 .
$$

The reduction $G = \sigma_{00}(\mathbf{a}'\boldsymbol\gamma)^2$ is zero
whenever the slopes are equal, $\gamma_t \equiv \bar\gamma$, because
then $\mathbf{a}'\boldsymbol\gamma = \bar\gamma \sum_t a_t = 0$.

*Proof.* Proposition 1 and $\operatorname{Var}(\mathbf{a}'\mathbf{Z}) =
\mathbf{a}'\operatorname{Cov}(\mathbf{Z})\mathbf{a}$. $\square$

**Consequence for the interaction.** The E9 slope regresses
$\Delta_i = \mathbf{a}'\mathbf{Y}_{ip}$ on $B_i$. Adding baseline to that
regression reduces the residual variance by the factor
$1 - \rho^2_{\Delta,\mathrm{BL}}$, where
$\rho^2_{\Delta,\mathrm{BL}} = (\mathbf{a}'\boldsymbol\sigma_{p0})^2 /
(\sigma_{00}\,\mathbf{a}'\Sigma_{pp}\mathbf{a})$. With equal slopes this
is zero, and the extra parameter costs one degree of freedom for
nothing.

**In the construct** $\mathbf{a}'\boldsymbol\sigma_{p0} = 0$ exactly, and
the variance of the E9 contrast is the same with and without
conditioning: 44.13 for paths 1 and 2, 45.03 for paths 3 and 4
(verified).

**When the covariate does help a within-participant contrast.**
$\mathbf{a}'\boldsymbol\gamma \neq 0$ requires baseline to predict
on-drug visits differently from off-drug visits. Two mechanisms give
that:

- **Baseline severity modifies the drug response.** Then $\gamma_t$ is
  larger at on-drug visits, and the appropriate model also includes
  $\mathrm{BL}_i \times D_{it}$.
- **Baseline predicts a time-varying natural course.** An example is
  regression to the mean that is strongest early in the trial. Then
  $\gamma_t$ changes with time, and in the Hybrid design drug state is
  confounded with time: on drug early, off drug late. So
  $\mathbf{a}'\boldsymbol\gamma \neq 0$.

Measurement error in baseline alone does not qualify. It attenuates
every slope by the same reliability factor, so the slopes stay equal
and $G$ stays zero.

### Proposition 3: between participants, the classic ANCOVA gain

For a comparison of post-baseline means between participants, such as
the biomarker main effect or a difference between paths, the relevant
variance is that of a post-baseline value. For a single visit $t$, the
three classic analyses [2, 3] have

$$
\operatorname{Var}(Y_{it}) = \sigma_{tt}, \qquad
\operatorname{Var}(Y_{it} - Y_{i0}) = \sigma_{tt} + \sigma_{00} - 2\sigma_{t0}, \qquad
\operatorname{Var}(Y_{it} \mid Y_{i0}) = \sigma_{tt} - \frac{\sigma_{t0}^2}{\sigma_{00}},
$$

and the conditional variance is never larger than the other two. It
equals the change-score variance exactly when $\gamma_t = 1$.

*Proof.* $\operatorname{Var}(Y_{it} - c\,Y_{i0})$ is a quadratic in $c$,
minimized at $c = \gamma_t$, with minimum $\sigma_{tt} -
\sigma_{t0}^2/\sigma_{00}$. $c = 0$ gives the first expression and
$c = 1$ the second. $\square$

**In the construct,** for a blinded visit: 599 unadjusted, 257 as change
or conditional on baseline (verified). That is a 57% reduction, the
largest virtue of baseline adjustment. It accrues to between-participant
terms, which in this design are the biomarker main effect and the path
means. Neither is the target.

### Proposition 4: the baseline row as an off-drug visit

The response coding counts baseline as an off-drug observation. For a
path's paired-difference contrast, that is the E9b statistic, with
off-drug mean $(Y_{i0} + \sum_{\text{off}} Y_{it})/(n_{\text{off}} + 1)$.
Substituting $Y_{i0} = \mathrm{BL}_i$ and $Y_{it} = \mathrm{BL}_i -
S_{it}$:

$$
\Delta_i^{\mathrm{E9}} = -\bar S_{i,\text{on}} + \bar S_{i,\text{off}},
\qquad
\Delta_i^{\mathrm{E9b}} = -\bar S_{i,\text{on}} + \frac{n_{\text{off}}}{n_{\text{off}} + 1}\,\bar S_{i,\text{off}} .
$$

Baseline cancels from both. The only difference is that E9b shrinks the
off-drug mean of the responses toward zero, the value they have at
baseline. Two consequences follow.

**(a) The interaction slope.** With $g_t$ the coupling of the biomarker
to the drug response at visit $t$ (`docs/40`, Section 5),

$$
\beta^{\mathrm{E9}}_{bm:D} = -c_{bm}\frac{\sigma_{BR}}{\sigma_{bm}}\bigl(\bar g_{\text{on}} - \bar g_{\text{off}}\bigr),
\qquad
\beta^{\mathrm{E9b}}_{bm:D} = -c_{bm}\frac{\sigma_{BR}}{\sigma_{bm}}\Bigl(\bar g_{\text{on}} - \frac{n_{\text{off}}}{n_{\text{off}} + 1}\,\bar g_{\text{off}}\Bigr).
$$

- **Without carryover** ($\bar g_{\text{off}} = 0$) the two coincide.
- **Under the published step rule with any carryover**
  ($\bar g_{\text{on}} = \bar g_{\text{off}} = 1$), the E9 slope is zero.
  The E9b slope is $1/(n_{\text{off}} + 1)$ of the no-carryover value,
  produced entirely by the uncoupled baseline row. That is the residual
  power floor of the published analysis (0.081 in the closed form; 0.08
  to 0.10 for `lmer`; `docs/37`, Section 6). It measures a
  pre-treatment comparison, not an interaction during the trial.

**(b) The precision.** Write $S_{it} = s_i + d_{it}$, where $s_i$ is the
participant's shared response level, with variance $\sigma_s$ equal to
the common covariance of the responses across visits. E9 removes $s_i$
exactly, since its weights on the responses sum to zero. E9b leaves
$-s_i/(n_{\text{off}} + 1)$ in the contrast, which adds
$\sigma_s/(n_{\text{off}} + 1)^2$ to its variance, apart from smaller
changes in the visit-specific terms. In the construct $\sigma_s \approx
200$, and the contrast variance rises from 44.13 to 51.04 (paths 1 and
2, $n_{\text{off}} = 3$) and from 45.03 to 51.77 (paths 3 and 4,
$n_{\text{off}} = 4$): increases of 16% and 15% (verified). In the
closed form this costs 0.06 to 0.07 in power (Section 6).

### Proposition 5: the working covariance and the validity of the standard errors

For a working covariance $V$, the GLS estimator
$\hat{\boldsymbol\beta}_V = \bigl(\sum_i X_i'V^{-1}X_i\bigr)^{-1}\sum_i
X_i'V^{-1}\mathbf{Y}_i$ is unbiased under a correct mean model, with

$$
\operatorname{Var}(\hat{\boldsymbol\beta}_V) = A^{-1} M A^{-1},
\qquad A = \sum_i X_i' V^{-1} X_i, \quad M = \sum_i X_i' V^{-1} \Sigma V^{-1} X_i .
$$

The model-based variance reported by the software is $A^{-1}$. It is
correct when $V \propto \Sigma$, and otherwise can be too large or too
small. The sandwich estimator of $A^{-1}MA^{-1}$ is consistent for any
$V$.

*Proof.* Standard; $\hat{\boldsymbol\beta}_V$ is linear in
$\mathbf{Y}_i$ with $E(\mathbf{Y}_i) = X_i\boldsymbol\beta$. $\square$

**In the construct,** the published random intercept implies $V =
\sigma_u^2 J + \sigma_e^2 I$ over the nine rows. Every row of $V$ then
carries the same residual variance $\sigma_e^2 > 0$, whereas the baseline
row of $\Sigma$ has none ($\sigma_{00} = \sigma_{t0}$). No choice of
$(\sigma_u^2, \sigma_e^2)$ matches it, so the model-based standard errors
are not guaranteed correct. The residual mismatch between open-label
and blinded visits remains. Dropping the baseline row (covariate coding)
or modeling all nine rows with an unstructured covariance (cLDA) removes
the baseline part of the mismatch. Only the unstructured covariance, or
a sandwich estimator, removes all of it.

### What this means for each quantity

| Quantity | Effect of baseline as a covariate | Result |
|-----------------------------|------------------------------------------|--------|
| Interaction $\beta_{bm:D}$ | no efficiency gain from the covariate when slopes are equal; gains from no longer counting baseline as off drug | Propositions 2 and 4 |
| Drug main effect $\beta_D$ | no gain for its within-participant part; loses baseline's off-drug information | Propositions 2 and 4 |
| Biomarker main effect $\beta_B$ | large gain (57% variance reduction per visit) | Proposition 3 |
| Model-based standard errors | removes the baseline row's incompatibility with compound symmetry | Proposition 5 |
| Robustness | protects the interaction when baseline predicts on and off visits differently | Proposition 2 |

Table: Effect of baseline as a covariate on each estimated quantity

## 5. Advantages of baseline as a covariate

1. **It removes the row the published covariance cannot represent**
   (Proposition 5). The baseline row has the smallest variance (342), no
   noise of its own, and lower correlation with later visits (0.69 to
   0.76, against 0.85 to 0.88 among later visits). That is the largest
   departure from compound symmetry in the outcome covariance.
2. **It stops a pre-treatment value from diluting the on-off contrast**
   (Proposition 4). Without carryover this improves the contrast's
   precision by 15% to 16% in variance. Under the step rule it removes a
   residual interaction that does not describe the trial.
3. **It protects the interaction when baseline is informative about the
   course of response** (Proposition 2). A covariate removes the part of
   the contrast predictable from baseline when $\mathbf{a}'\boldsymbol\gamma
   \neq 0$: for example, when severity moderates the drug response or
   predicts a time-varying natural course. The construct excludes both;
   real trials may not.
4. **It estimates the baseline coefficient instead of fixing it.** With
   measurement error or regression to the mean, the slopes fall below
   1. A covariate estimates them; the response coding under a random
   intercept implicitly ties baseline to the same participant level as
   every later visit. ANCOVA is the standard choice for this reason
   [2, 3].
5. **It is the convention reviewers expect.** Baseline-adjusted analysis
   of post-baseline visits, usually as an MMRM, is the default in
   confirmatory trials, and regulatory guidance on covariate adjustment
   assumes it [4].
6. **It opens a direct test of severity as a moderator.** A
   $\mathrm{BL}_i \times D_{it}$ term is the natural extension. In the
   construct it is exactly zero.

## 6. Disadvantages of baseline as a covariate

1. **It gives the interaction no efficiency gain in this construct**
   (Proposition 2). Any power change comes from dropping the baseline
   row, not from the covariate, and costs one degree of freedom.
2. **It discards baseline's off-drug information for the drug main
   effect.** Baseline is the only off-drug value with no expectancy,
   carryover or natural-course change. How much that matters depends on
   the time adjustment: with visit effects, baseline has its own mean and
   carries no drug information anyway (`docs/42`).
3. **It changes the published estimand.** Results are no longer directly
   comparable with Hendrickson et al.'s, including the power floor under
   the step rule (Proposition 4a).
4. **A participant without a baseline value is lost** unless it is
   imputed. cLDA is more efficient than ANCOVA when baselines are
   missing [5, 6]. In this design baseline precedes randomization and
   should rarely be missing.
5. **The covariate does not fix the post-baseline covariance**
   (Proposition 5). Open-label and blinded visits still differ (variance
   368 against 257 after conditioning, covariance 263 against 185). A
   random intercept on the post-baseline visits remains misspecified; the
   covariate coding is fully effective only with `csh` or `us`.

**Closed-form evidence on the baseline row** (E9 against E9b; Hybrid,
compound symmetry, $N = 70$, strength 0.25; verified, `docs/37`,
Sections 4 and 6):

| Rule, $t_{1/2}$ | E9 power | E9b power | E9 SD | E9b SD |
|------------------------|----------|----------|----------|----------|
| Step, 0 | 0.721 | 0.654 | 0.0507 | 0.055 |
| Step, 0.1 or more | 0.050 | 0.081 | 0.0531 | 0.057 |
| Mean moderation, 0 to 1 | 0.681 | 0.619 | 0.0531 | 0.057 |

Table: Power and SD of E9 and E9b by carryover rule and $t_{1/2}$

These are results for the paired-difference statistic. They show the
direction and size of Proposition 4's effects. A mixed model weights the
visits by its working covariance, so its numbers depend on that
covariance (Section 7).

## 7. Constrained longitudinal data analysis: when the choice matters

Liang and Zeger proposed keeping baseline in the response vector,
constraining its mean to be common across randomization groups, and
modeling the covariance of all visits [5].

**With complete data and an unstructured covariance, cLDA and ANCOVA
with visit-specific slopes give the same estimates of the post-baseline
parameters.** The joint likelihood factorizes as $f(Y_{i0}) \times
f(\mathbf{Y}_{ip} \mid Y_{i0})$. Under an unstructured $\Sigma$, the
parameters of the conditional factor ($\boldsymbol\beta$,
$\boldsymbol\gamma$, $\Sigma_{pp\cdot 0}$) vary independently of those of
the marginal factor ($\mu_0$, $\sigma_{00}$), apart from any parameter
shared by the two means (Proposition 1). In the Hybrid design every path
is off drug at baseline with the same time and expectancy, so
randomization's constraint holds automatically. The only shared mean
parameter is the biomarker main effect $\beta_B$, and allowing it to
differ at baseline restores exact equivalence. The equivalence is a
property of the likelihood; it has not been simulated here. Lu [6]
shows that cLDA is more efficient when baselines are missing.

**The choice matters most under the published model.** With a random
intercept and linear time, the baseline row is not modeled by its own
mean or covariance. It enters as an off-drug visit at $t = 0$ under a
covariance it does not have (Proposition 5), it dilutes the on-off
contrast (Proposition 4), and it anchors the time trend.

| Analysis | Effect of the baseline coding |
|------------------------------------------|----------------------------------------|
| Random intercept, linear `t` (published) | large: the row is misfit, dilutes the contrast and has leverage |
| Random intercept, visit factor | moderate: the row has its own mean, but the covariance is still wrong |
| Unstructured covariance, visit factor | none for complete data: cLDA and ANCOVA coincide |

Table: Effect of the baseline coding under each analysis model

## 8. Recommendation

**For a Hybrid trial analysis.** Analyze the eight post-baseline visits
with baseline as a covariate, in an MMRM with an unstructured covariance
and visit effects. Allow visit-specific baseline slopes if their
variation is plausible, and consider $\mathrm{BL}_i \times D_{it}$ as a
secondary moderator. With complete baselines this coincides with cLDA,
which is the more efficient form when baselines may be missing. The
covariate's value for the interaction is robustness against baseline
predicting the course of response (Proposition 2), not efficiency. Avoid
the published combination of a baseline row, a random intercept and
linear time.

**For simulations that reproduce Hendrickson et al.** Keep the response
coding, since it defines the published estimand. Report the covariate
coding alongside it wherever the comparison bears on conclusions.

## 9. Open questions

The paused covariance study would answer these if the covariate coding
were added: the post-baseline `cs` and `us` analyses with each time
adjustment, about 40% more run time.

1. **Size.** Does the covariate coding remove the anticonservatism of
   the flexible time adjustments in Hybrid (`docs/40`, Section 9)?
   Proposition 5 predicts that it removes part of the cause.
2. **Power.** Does the mixed-model gain from dropping the baseline row
   match Proposition 4's 15% to 16% variance reduction?
3. **Equivalence.** Do cLDA and ANCOVA coincide under an unstructured
   covariance, as Section 7 derives?
4. **Robustness.** Under a construct in which baseline predicts the
   natural course or moderates the drug response, how large is
   Proposition 2's gain? That requires extending the construct; the
   software has no baseline-to-component link (`docs/31`, Section 18.1).

## 10. Limitations and evidence status

- **Derivations and numeric checks, no Monte Carlo.** Propositions 1 to
  5 are standard multivariate-normal and GLS results applied to this
  design (derived). Their numeric consequences were computed exactly
  from the construct's implied covariance (verified,
  `analysis/scripts/quick-sim/hybrid-design-primer/02-baseline-propositions.R`).
  The power figures are closed-form results for E9 and
  E9b, validated by Monte Carlo in `docs/37` (verified). Statements
  about regression to the mean and measurement error come from the
  literature [2, 3].
- **One construct.** The published compound-symmetry construct has no
  baseline measurement error, no regression to the mean and no link
  between baseline and the response. Each would strengthen the case for
  the covariate (Proposition 2).
- **One design.** Hybrid at $N = 70$. In CO the baseline row is equally
  uncoupled, but the paths differ in drug state from the first visit.

## 11. References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020; 2:13. doi:10.3389/fdgth.2020.00013.
2. Vickers AJ, Altman DG. Statistics notes: analysing controlled trials
   with baseline and follow up measurements. *BMJ* 2001;
   323(7321):1123-1124. doi:10.1136/bmj.323.7321.1123.
3. Senn S. Change from baseline and analysis of covariance revisited.
   *Statistics in Medicine* 2006; 25(24):4334-4344. doi:10.1002/sim.2682.
4. European Medicines Agency. Guideline on adjustment for baseline
   covariates in clinical trials. EMA/CHMP/295050/2013, 2015. U.S. Food
   and Drug Administration. Adjusting for covariates in randomized
   clinical trials for drugs and biological products: guidance for
   industry, 2023.
5. Liu GF, Lu K, Mogg R, Mallick M, Mehrotra DV. Should baseline be a
   covariate or dependent variable in analyses of change from baseline
   in clinical trials? *Statistics in Medicine* 2009; 28(20):2509-2530.
   doi:10.1002/sim.3639. (Describes the cLDA model of Liang KY and
   Zeger SL, *Sankhya Series B* 2000; 62:134-148.)
6. Lu K. On efficiency of constrained longitudinal data analysis versus
   longitudinal analysis of covariance. *Biometrics* 2010; 66(3):891-896.
   doi:10.1111/j.1541-0420.2009.01332.x.
