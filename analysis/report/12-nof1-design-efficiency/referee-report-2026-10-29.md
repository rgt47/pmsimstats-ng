# Referee report: paper 12, closed-form design efficiency
*2026-10-29 12:26 PDT*

**Manuscript.** 'A closed-form approximation to the biomarker-treatment
interaction estimator in aggregated N-of-1 trials, and what it implies
for design efficiency under a mandatory open-label lead-in'
(`analysis/report/12-nof1-design-efficiency/report.Rmd`, last edited
2026-10-08; rendered `report.pdf` of the same date).

**Submitted as.** Short communication, *Statistics in Medicine*
(cover letter).

**Reviewer.** pmsimstats team, internal referee pass.

**Evidence labels.** VERIFIED: I executed code and observed the result.
INSPECTED: I read the source and confirmed it. SUSPECTED: consistent
with what I read but not confirmed. UNVERIFIED: not checked.

**Scratch scripts** (not part of the repository):
`p12-check1.R` (reproduction of every inline number and the null
cells), `p12-check2.R` (exact variances under the stored DGP,
positive-definiteness check, attenuation MCSE), `p12-check3.R`
(covariance rebuilt under the current DGP with repair refused), all in
`/private/tmp/claude-501/-Users-zenn-prj-res-36-pmsimstats-ng-pmsimstats-ng-analysis-report-01-dgp-mean-moderation-vs-mvn/513155aa-ec85-4a1d-889a-5d57e8607d18/scratchpad/`.

## 1. Summary

The manuscript proposes a closed-form screen for the precision and
power of the biomarker-by-drug interaction coefficient in aggregated
N-of-1 designs. It reduces the analysis model to a per-participant
on-minus-off contrast, approximates the contrast variance by the sum of
the AR(1) variances of the on-drug and off-drug visit means ($v_1$),
and multiplies the effect by a carryover attenuation factor
$A = 1 - \overline{D_{bc}}^{\,\text{off}}$. One scale parameter
$\sigma^2$ is fitted by least squares through the origin to the
empirical $\mathrm{sd}(\hat\beta)^2 N$ of 18 stored Monte Carlo cells
from the dual-channel driver (mean channel only, specification E1,
three designs, $N \in \{35, 70\}$, $t_{1/2} \in \{0, 0.5, 1\}$ week,
1,000 replicates). The manuscript reports $R^2 = 0.948$ for the
variance fit, a power mean absolute error of 0.049, a correct design
ranking by $v_1$, and an attenuation factor that ranks designs
correctly but overstates the magnitude. It closes with design
guidance for trials whose first phase must be open label.

The question is useful and the idea of a visit-schedule screen is
sound. However, the central quantitative claims do not survive
checking against the stored artifacts. The reported $R^2$ is an
uncentered statistic from a no-intercept fit. The centered $R^2$ is
$-0.49$, so the closed form predicts the 18 variances worse than their
mean. The $v_1$ ranking of the designs is inverted for the crossover.
The attenuation factor describes a carryover of the interaction that
the validating DGP does not contain. The 'order-of-magnitude
calibration question' is a biomarker-units issue that resolves
exactly. A parameter-free closed form exists: the GLS sandwich
variance under the DGP's own covariance. It predicts the empirical
standard deviations to within 0 to 10 percent with no fitted
parameter. The compendium already derives it (Section
'Separation of the biomarker block'), and paper 01 now contains its
own exact E9 power result. The validation data are also stale: they
were archived on 2026-10-29, they predate the `c.cfct` correction, and
they use the compounding carryover rule.

## 2. Overall assessment and recommendation

**Recommendation: reject in current form. Resubmission is encouraged
after the closed form is rebuilt on the GLS variance and revalidated
against current-DGP output.**

The defects are not presentational. Three of the four headline
results are incorrect as stated: the variance fit, the design ranking
and the attenuation result. The fourth, the power agreement, is
largely in-sample, because both the scale and the effect size are
taken from the cells being predicted. The design guidance in Section 5
depends on these results and on paper-01 numbers that exist only in
an archived draft. A correct version of the paper is within reach, and
most of the required algebra is already in the compendium. That
version would be a different paper: a validated GLS screen for design
choice under a known covariance, with the expectancy confound of an
open-label lead-in treated explicitly.

## 3. Significance and novelty

- **Within the field.** Closed-form power for within-subject designs
  under serial correlation is established. Wang and Schork (2019)
  derived analytic power for crossover-based N-of-1 trials under AR(1)
  errors and carryover. Optimal crossover design under autocorrelated
  errors goes back to Kunert (1985) and Matthews (1987), with
  orthogonal-array optimality results in Kunert and Martin (2000).
  Senn (2019) gives the main-effect and patient-by-treatment variance
  formulas for N-of-1 series. Closed-form power for interaction terms
  in mixed models exists for longitudinal designs (Heo and Leon). None
  of these is cited, except Senn. The specific combination here (a
  continuous biomarker by drug interaction, irregular Hendrickson
  schedules, an open-label lead-in) is new in its application. The
  method, GLS information for a covariate-by-treatment term under a
  known covariance, is standard.
- **Within the compendium.** The claim to be 'the first to attempt a
  closed-form treatment' of this estimand (README; Introduction) is no
  longer accurate. Paper 01 (current `report.Rmd`, lines 405 to 450
  and Table at line 709) derives exact E9 and E9b slopes and
  noncentral-$t$ power and validates them against 5,000-replicate
  Monte Carlo cells. The compendium derives the GLS separation of the
  biomarker block (`eq:bp-plim`) and the most efficient contrast
  (`eq:invf-optcontrast`). Both are more general and more accurate
  than $v_1$. The compendium already annotates paper 12 as needing
  correction (`compendium-mathematics.tex`, subsection 'A closed-form
  screen from visit schedules'). INSPECTED.
- **Contemporariness.** The validation uses specification E1 with
  model-based standard errors. The compendium's current analysis
  recommendation is F1 (biomarker centered, expectancy terms, CR2 with
  Satterthwaite df), and paper 02 now recommends G9. A screen
  calibrated to E1 with model-based inference does not reflect the
  analysis a reader would be advised to run. INSPECTED.

## 4. Major comments

### 4.1 The reported $R^2$ is uncentered; the fit is worse than a constant

- **Location.** Section 4.2, first paragraph; Conclusion ('reproduces
  95% of the variation'); cover letter ('explains 95% of the
  variance'); README.
- **Problem.** The fit `lm(I(sd_beta^2 * N) ~ 0 + v_sum)` has no
  intercept. For such a fit, `summary.lm` reports
  $1 - \mathrm{RSS}/\sum y^2$, not $1 - \mathrm{RSS}/\sum(y - \bar y)^2$.
  The outcome varies over a narrow band (0.49 to 0.88, coefficient of
  variation 0.20), so almost any positive regressor gives an
  uncentered $R^2$ near 1. The centered $R^2$ is $-0.493$. The
  residual standard deviation of the closed form, 0.157, exceeds that
  of the constant-only model, 0.132. VERIFIED (`p12-check1.R`).
  The 18 cells also carry only three distinct $v_1$ values, one per
  design, because $v_1$ does not depend on $t_{1/2}$. The fit is
  therefore a three-point fit replicated six times, not an 18-point
  validation.
- **Why it matters.** This is the paper's main quantitative claim. As
  stated, the variance component 'performs well' is false. The closed
  form explains none of the between-design variation in precision.
- **Remedy.** Report centered $R^2$ or, better, the ratio
  $\mathrm{sd}(\hat\beta)/\widehat{\mathrm{sd}}$ with its Monte Carlo
  interval, per design. Replace $v_1$ with the exact contrast variance
  or the GLS variance (comment 4.3). Remove the 95 percent claim from
  the Conclusion, the cover letter and the README.

### 4.2 The $v_1$ design ranking is inverted for the crossover (correctness)

- **Location.** Section 4.2 ('the design ranking by $v_1$ ... is
  reproduced in direction'); Section 4.3 (the CO versus OL+BDC
  'close call'); Section 5, steps 2 and 4.
- **Problem.** $v_1$ ranks the crossover as the most precise design
  (0.909) and OL+BDC as the least precise (1.142). The empirical
  standard deviations rank them the other way in all six
  design-by-half-life comparisons. At $N = 70$, $t_{1/2} = 0$:
  Hybrid 0.084, OL+BDC 0.096, CO 0.110. VERIFIED. The text states the
  ranking is 'reproduced in direction'. Only the Hybrid versus OL+BDC
  order is correct. The CO versus OL+BDC disagreement is described as
  a close call between designs with 'similar' $v_1$. The two values
  differ by 26 percent, in the wrong direction, and the empirical
  variances differ by a ratio of 1.32. This is a sign error in the
  screen, not a close call.
- **Why it matters.** Section 5 recommends using the $v_1$ ranking to
  eliminate dominated designs. Applied to these three designs, it
  would retain the least efficient one and discard the second most
  efficient one.
- **Remedy.** Retract the ranking claim. The exact contrast variance
  ($v_1$ minus the on/off cross term, already given in the
  2026-10-08 correction note) ranks the designs correctly: Hybrid
  0.492, OL+BDC 0.714, CO 0.770. VERIFIED. Rebuild Section 5 on the
  exact form or on the GLS form.

### 4.3 A parameter-free closed form exists and fits the data

- **Location.** Sections 3.2 and 4.2; Limitations.
- **Problem.** Under the validating DGP, the participant covariance of
  the outcome is fully known: $Y_{it} = \mathrm{BL}_i - TV_{it} -
  PB_{it} - BR_{it}$, with the component correlation matrix built by
  the driver. Because the biomarker is independent of design and
  centered, the biomarker block of the GLS normal equations separates
  (compendium `eq:bp-plim`). The variance of the $B_i D_{it}$
  coefficient then follows from $\sum_p n_p Z_p' V_p^{-1} Z_p$ with
  $Z_p = (\mathbf 1, \mathbf D_p)$, or from the corresponding sandwich
  when the working covariance is a random intercept plus AR(1). I
  computed both from the stored DGP settings, with no fitted
  parameter (`p12-check2.R`). VERIFIED results:

  | Quantity (through-origin fit) | Uncentered $R^2$ | Centered $R^2$ |
  |---|---|---|
  | $v_1$ (manuscript) | 0.948 | $-0.493$ |
  | exact AR(1) contrast variance | 0.994 | 0.821 |
  | contrast variance under true $V$ | 0.984 | 0.537 |
  | AR(1) working-covariance sandwich | 0.997 | 0.905 |
  | GLS under true $V$ | 0.998 | 0.932 |

  Without any fitting, the AR(1) sandwich predicts the empirical
  standard deviation with ratios
  $\mathrm{sd}(\hat\beta)/\widehat{\mathrm{sd}}$ of 0.996 to 1.036 for
  Hybrid, 0.998 to 1.034 for OL+BDC and 1.044 to 1.095 for CO.
  The random intercept cancels exactly from the biomarker-by-drug
  information under any working covariance that contains it.
- **Why it matters.** The paper's justification for a fitted
  $\sigma^2$ ('an empirical quantity like any other nuisance variance
  component') does not hold. The variance components are inputs that a
  trialist must specify anyway, and the exact information costs no
  more to compute than $v_1$. The fitted $\sigma^2$ also has no stated
  interpretation or units (points squared per squared biomarker unit,
  absorbing $\sigma_{bm}^2$).
- **Remedy.** Present the GLS information (or its sandwich under the
  analysis model's working covariance) as the closed form. Validate it
  out of sample, without fitting anything to the cells it predicts.
  Explain the residual 5 to 9 percent excess for CO. A plausible
  source, SUSPECTED, is the uniform blinded expectancy of the CO
  design interacting with the linear-time mean model.

### 4.4 The attenuation factor models an effect absent from the DGP

- **Location.** Equations `delta-decomp` and `e9-like`; Section 4.4
  and Table `tab:attenuation`; Section 5, step 3; Conclusion.
- **Problem.** The derivation assumes the biomarker moderation is
  carried into off-drug visits in proportion to $D_{bc}$:
  $\beta_1(1 + \beta_{bm}B_i)\bar D_{bc,i}^{\,\text{off}}$. In the
  driver that produced the validation data, the mean-channel shift
  $\beta_{bm} b_i \sigma_{BR}$ is added to $BR$ at on-drug visits only
  (`implementations/tidyverse/R/functions.R`, Jul 31 version, lines
  426 to 441; unchanged in the current version, line 531). Carryover
  acts only on the $BR$ mean, which is not a function of the
  biomarker. INSPECTED. The interaction therefore has no carryover in
  this DGP, and E1's binary indicator is correctly specified for it.
  The stored estimates confirm this. The mean coefficient is $-0.235$
  across all 18 cells, against $-0.234$ implied by the DGP. Its
  change from $t_{1/2} = 0$ to $1$ at $N = 70$ is $-0.0072$ (CO),
  $+0.0041$ (Hybrid) and $+0.0059$ (OL+BDC), with unpaired MCSE 0.0049,
  0.0038 and 0.0043 ($|z| \le 1.5$). VERIFIED.

  Three statements in Section 4.4 follow from this and are incorrect.
  First, 'the direction is right in every case' is false: CO's
  magnitude rose (ratio 1.031). Second, the 'correct ranking' is a
  ranking of differences within Monte Carlo noise. Third, the
  structural explanation (the time covariate 'recovers gradient
  information within the off-drug period') is not the cause. There
  was no attenuation to recover. The cells at different half-lives
  also use different random-number streams
  (`furrr_options(seed = seed)` over grid rows), so the comparison is
  unpaired. INSPECTED.
- **Why it matters.** The paper presents a 'genuine, useful negative
  result' that is an artifact of applying an attenuation model to a
  DGP without the attenuated quantity. Section 5 then recommends $A$
  as a screen.
- **Remedy.** Either validate $A$ against a DGP in which the
  moderation decays with the drug effect (the covariance channel with
  the decaying correlation, or a mean channel whose shift scales with
  the carried-over $BR$ mean), or remove the attenuation component.
  Paper 01's Proposition 2 already gives the exact E9 slope under the
  decaying correlation and should be cited for this purpose. Report
  paired differences with paired MCSE if any cross-half-life
  comparison is retained.

### 4.5 The calibration question is a biomarker-units error

- **Location.** Limitations ('3.6 nominal against approximately 0.24
  observed ... an open calibration question for the companion
  architecture manuscripts').
- **Problem.** The analysis model uses the raw biomarker. The
  calibration in `data/extracted_bp.rda` has $\sigma_{bm} = 15.36$. The
  coefficient per raw biomarker unit is therefore
  $-c_{bm,A}\,\sigma_{BR}/\sigma_{bm} = -0.45 \times 8 / 15.36 =
  -0.2344$. The stored mean estimate is $-0.2353$. VERIFIED. There is
  no discrepancy. NOTATION.md lists the true value as
  $-c_{bm}\sigma_{BR}$ only for a standardized biomarker, and paper
  01's Proposition 2 carries the $1/\sigma_{bm}$ factor explicitly.
- **Why it matters.** The manuscript flags a nonexistent defect in the
  companion papers. It also uses the empirical mean estimate as the
  'effective true effect' to avoid a problem that does not exist.
- **Remedy.** Use the DGP's true value in the power formula. Delete
  the open-question paragraph. State the sign convention (NOTATION
  rule 5): the coefficients are negative.

### 4.6 The power validation is in-sample and conflates calibration

- **Location.** Section 4.3 and Table `tab:power-fit`; Conclusion.
- **Problem.** Both inputs to the predicted power are taken from the
  cells being predicted: $\hat\sigma^2$ is fitted to their standard
  deviations, and the effect is each cell's own mean estimate. The
  simulated power, however, uses model-based `lme` standard errors,
  whose calibration differs by design. In the null cells of the same
  run ($c_{bm,A} = c_{bm,B} = 0$), E1's rejection rate is 0.021 to
  0.033 for OL+BDC, 0.031 to 0.045 for Hybrid and 0.043 to 0.060 for
  CO (MCSE about 0.006). VERIFIED. Power computed from the empirical
  standard deviation and the true effect therefore exceeds simulated
  power for OL+BDC by 0.04 to 0.11, and falls below it for CO at
  several cells. VERIFIED. Part of the reported error pattern is test
  conservatism ($\kappa > 1$), the subject of paper 10, not design
  efficiency.
- **Why it matters.** A screen for design efficiency should predict
  the precision of the estimator, and separately the operating
  characteristics of a calibrated test. Mixing the two makes neither
  claim checkable.
- **Remedy.** Validate the variance against $\mathrm{sd}(\hat\beta)$
  and power against a calibrated test (CR2 with Satterthwaite df, the
  compendium's recommended inference). Report $\kappa$ by design.
  Give MCSE for every power value (STYLE section 1).

### 4.7 Stale validation data and a broken reproducibility chain

- **Location.** Section 4.1; Reproducibility; title-page data
  statement.
- **Problem.** The validation reads
  `analysis/data/dgp-combined/01-combined-summary.rds`. That file was
  moved to `archive-2026-07-31/` in the working tree, while paper 11
  is being rerun. `01-closed-form-validation.R` now fails at its
  `readRDS` path. INSPECTED. The report itself reads only the derived
  CSV of 2026-08-20, so it renders, but the stated chain does not run.
  The archived run used the July DGP (`simulation-core.R` at commit
  26b8c16a3): `c.cfct = 0.1` rather than 0.2, the compounding
  carryover rule rather than carryover from the last on-drug mean,
  and `make_positive_definite = TRUE`. INSPECTED.

  Impact checks, by item:

  - **Positive definiteness.** All nine design-by-half-life
    covariance matrices are positive definite under the July settings
    (minimum eigenvalue 2.36). They are also positive definite under
    the current settings with repair refused (minimum eigenvalue
    4.49). VERIFIED. No stored cell depended on a repair.
  - **Cross-component correlation.** Under `c.cfct = 0.2`, the
    predicted standard deviation at $N = 70$ falls by about 6 percent
    for Hybrid (0.084 to 0.079) and 7 percent for OL+BDC (0.096 to
    0.089). CO is unchanged (0.102 to 0.103). This is INFERRED from
    the sandwich formula, not re-simulated. The stored numbers are
    therefore stale in magnitude but probably not in direction.
  - **Carryover rule.** It affects only the main-effect mean under
    this DGP, so its effect on the interaction is expected to be nil.
    SUSPECTED.
  - **Pre-exposure guard.** The July `prepare_long_data` sets
    $D_{bc} = 1$ for pre-exposure visits at $t_{1/2} > 0$. That error
    affects E2, not the E1 results used here. INSPECTED.
- **Why it matters.** A short communication's whole evidence base is
  one stored run. That run must be current and reproducible.
- **Remedy.** Point the script at a current-DGP run, ideally paper
  11's rerun or a dedicated paired run. Use the same seeds across
  designs and half-lives. Record the commit and DGP settings in the
  Reproducibility section.

### 4.8 E1 with model-based errors is not the recommended analysis

- **Location.** Sections 3.1 and 4.1.
- **Problem.** The fitted model is `Sx ~ bm + t + Db + bm:Db`, with a
  random intercept, `corCAR1` and model-based standard errors. If the
  `corCAR1` fit fails, the driver silently refits without the
  correlation structure and counts the fit as converged
  (`fit_spec`, July version). INSPECTED. How often this fallback
  occurred is UNVERIFIED, because the summary records `conv_rate = 1`
  in every cell. The estimand is the same across designs (the raw-unit
  $\beta_{bm:D}$, VERIFIED unbiased in all cells), which is correct.
  The mean model omits the expectancy term. For an open-label lead-in,
  that term is the central threat. On-drug status coincides with
  expectancy 1 in the open-label visits. A biomarker correlated with
  the placebo response then loads onto `bm:Db`. The compendium's
  recommended F1 model includes `bmc:De` for exactly this reason
  (CLAUDE.md, 'Analysis model'; `docs/47` to `docs/49`, `docs/56`).
- **Why it matters.** The paper's motivating question concerns the
  open-label-first constraint. The guidance never addresses the bias
  that the constraint creates. It only notes, in the last sentence of
  the Limitations, that a placebo-correlated biomarker is outside its
  scope.
- **Remedy.** Base the screen on the F1 mean model (whose biomarker
  block is $Z = (\mathbf 1, \mathbf D, \boldsymbol\eta)$). Report the
  precision cost of adding `bmc:De` for each open-label-first design.
  That cost is computable in the same GLS framework and is the
  quantity a trialist facing this constraint actually needs.

### 4.9 Design guidance rests on archived numbers and an infeasible design

- **Location.** Section 2 ('0.842 versus 0.751 ... 0.784 versus
  0.730'); Section 5.
- **Problem.** These powers appear in no current manuscript. They come
  from the archived pre-restatement draft of paper 01
  (`archive/report-2026-10-07-pre-restatement.Rmd`, lines 760 to
  766). Paper 01's current text contains no mean-moderation design
  comparison. INSPECTED. They also disagree with this paper's own
  stored cells for the same architecture and $N = 70$ (Hybrid 0.738,
  OL+BDC 0.612, CO 0.595 at $t_{1/2} = 0$), and the difference is not
  explained. Separately, the crossover has no open-label phase and is
  therefore not a candidate under the motivating constraint. The step
  'adding visits after the mandatory open-label phase is the largest
  lever in $v_1$' is not shown anywhere in the paper.
- **Remedy.** Remove or re-source the paper-01 figures. Restrict the
  design comparison to designs that begin open label (for example,
  OL+BDC and Hybrid with varying discontinuation and crossover
  schedules). Use the GLS screen to rank them, and verify the ranking
  with a paired simulation.

### 4.10 Mis-described companion results

- **Location.** Introduction and Section 2.
- **Problem.**
  - The dichotomized-biomarker RM-ANOVA closed form is attributed to
    paper 10. It is derived in paper 08 (paper 08 title and Methods,
    line 53). The README says paper 08.
  - Section 2, 'Fourth', attributes to paper 10 the finding that a
    misspecified residual correlation is the larger
    anti-conservatism risk. Paper 10's current Conclusion says the
    opposite: the conservatism came mainly from the mean model, not
    the residual correlation.
  - The 'general four-channel decomposition ... its $\kappa = 1$,
    $\nu = \infty$ special case' does not appear in paper 10's current
    text (searched for 'four-channel').

  All three findings are INSPECTED.
- **Remedy.** Recheck every companion citation against the current
  text of the cited paper.

## 5. Minor comments

1. **No abstract.** The manuscript has no abstract, even though it is
   a short communication. STYLE section 5 requires an abstract check.
2. **No glossary.** STYLE section 4 requires a glossary of analysis
   codes (E1, E9), design names and project quantities. None exists.
3. **Notation conflicts with NOTATION.md** (not caught by the linter,
   which reported the file clean):
   - $\eta_i$ is used for the residual difference. $\eta$ is reserved
     for the design expectancy.
   - $A$ is used for the attenuation factor. $A$ is the temporal
     correlation matrix.
   - $k_{\text{on}}$ and $k_{\text{off}}$ should be $n_{\text{on}}$
     and $n_{\text{off}}$, since $k$ is the cycle count and the
     Weibull shape.
   - $\beta_1$ should be $\beta_D$.
   - $\gamma(h)$ and $\gamma_0$ collide with paper 06's $\gamma$.
   - Calendar time in Equation `var-mean-ar1` should be $w_t$, not
     $t$.
   - $\mathrm{Var}(b_i) = 1$ is assumed although the fitted model uses
     the raw $B_i$.
   - The multiplicative form $\beta_1(1 + \beta_{bm}B_i)$ is not the
     canonical additive shift $\beta_{bm}\sigma_{BR}b_i$ (rule 2). The
     compendium already flags this.
4. **Embedded correction notes.** The two 'Corrected 2026-10-08'
   notes contradict the surrounding prose they correct. A submitted
   paper cannot carry both. The note's per-path values (0.926/0.464,
   0.910/0.451) and the CO values (0.909/0.770) reproduce. VERIFIED.
5. **Hard-coded section reference.** '§4.4' (Section 5) should be a
   `\ref`.
6. **Prose against STYLE.** Many sentences exceed 30 words, and
   semicolons join clauses. Evaluative words recur: 'honestly',
   'genuine' (three times), 'crude', 'badly', 'simply', 'reported
   rather than concealed', 'not ... smoothed over'. The STYLE grep
   flags 'ceiling' at lines 176 and 179 (outside code). 'Artefact'
   (line 470) is British spelling.
7. **Missing MCSE.** No table reports MCSE, although the summary
   holds `mcse_power`. The ratio column in Table `tab:variance-fit`
   needs an interval.
8. **Rounding.** The $1/N$ check ('1.439 against 1.414') is an
   expected property of any estimator. It does not confirm the closed
   form.
9. **Unused bibliography entries.** `hendricksonOptimizing2020` and
   `morris2019simulation` are never cited, so pandoc drops them.
   Hendrickson et al. must be cited as the source of the Hybrid design
   (STYLE section 2). Jones and Kenward is not cited.
10. **Stale companion titles.** Six of the seven `@unpublished` titles
    no longer match the current manuscripts (papers 01, 02, 05, 09,
    10 and 11).
11. **Front matter.**
    - The title page names a placeholder corresponding author with a
      personal email address. This conflicts with the project's
      'pmsimstats team' convention.
    - The repository URL is a placeholder ('[organisation]', a British
      spelling).
    - The cover letter repeats the false 95 percent claim.
    - The README says '`rgt` paragraphs are placeholders', which is no
      longer true. The `in_header: claudecode.tex` line is an inert
      leftover.
12. **Units.** $t_{1/2}$ is in weeks, which is consistent with
    NOTATION. No `**Units.**` paragraph is needed.
13. **Rendered artifact.** `report.pdf` (2026-10-08) and `report.tex`
    sit in the directory. Rendering should go through
    `tools/render.sh`, which stages a stamped copy in `share/`. No
    stamped copy of paper 12 exists in `share/MANIFEST.md`.
    SUSPECTED; only the file name was checked.

## 6. Missing and questionable references

| Reference | Status | Relevance |
|---|---|---|
| Wang Y, Schork NJ. Power and design issues in crossover-based N-of-1 clinical trials with fixed data collection periods. *Healthcare* 2019;7(3):84. doi:10.3390/healthcare7030084 | Missing (PubMed 31269712) | Analytic AR(1) power for N-of-1 designs, with carryover; the closest prior work |
| Percha B, Baskerville EB, Johnson M, et al. Designing robust N-of-1 studies for precision medicine. *J Med Internet Res* 2019;21(4):e12641. doi:10.2196/12641 | Missing (PubMed 30932871) | Design recommendations under autocorrelation and carryover |
| Blackston JW, Chapple AG, McGree JM, et al. Comparison of aggregated N-of-1 trials with parallel and crossover RCTs using simulation. *Healthcare* 2019;7(4):137. doi:10.3390/healthcare7040137 | Missing (PubMed 31698799) | Aggregated N-of-1 versus crossover efficiency |
| Araujo A, Julious S, Senn S. Understanding variation in sets of N-of-1 trials. *PLoS One* 2016;11(12):e0167167. doi:10.1371/journal.pone.0167167 | Missing (PubMed 27907056) | Mixed-model structure and purpose of analysis for N-of-1 series |
| Matthews JNS. Optimal crossover designs for the comparison of two treatments in the presence of carryover effects and autocorrelated errors. *Biometrika* 1987;74:311-320 | Missing (web search) | Optimal designs under AR(1) errors |
| Kunert J. Optimal repeated measurements designs for correlated observations and analysis by weighted least squares. *Biometrika* 1985;72:375-389 | Missing; content not confirmed | GLS design optimality under correlation |
| Jones B, Kenward MG. *Design and Analysis of Cross-Over Trials*, 3rd ed. CRC Press, 2014 | Missing | Standard reference for crossover variance and carryover |
| Heo M, Leon AC. Sample sizes required to detect two-way and three-way interactions involving slope differences in mixed-effects linear models. *J Biopharm Stat* 2010 | Missing; details UNVERIFIED | Closed-form interaction power in mixed models |
| Hendrickson RC, et al. 2020, *Front Digit Health* | In bib, never cited | Source of the Hybrid design; must be cited |
| Senn S. Sample size considerations for n-of-1 trials. *SMMR* 2019;28(2):372-383 | Cited; supports main-effect formulas | Bib correct; epub 2017 |
| Zucker DR, Ruthazer R, Schmid CH. *J Clin Epidemiol* 2010;63(12):1312-1323 | Cited; questionable | Compares combining models; it does not supply a 'classical variance formula' |
| pmsimstats papers 01, 08, 10 | Questionable as cited | See comments 4.9 and 4.10 |
| Morris, White, Crowther 2019 | In bib, never cited | Should be cited if an ADEMP-style validation is added |

PubMed was searched for Percha, Blackston, Wang and Schork, Araujo,
and Zucker. Details for Matthews, Kunert and Martin, and Heo and Leon
come from a web search and have not been checked against the
originals.

## 7. Suggestions for strengthening

1. Rebuild the screen as the GLS (or working-covariance sandwich)
   variance of the biomarker-by-drug coefficient, using the separated
   biomarker block. Give it as one displayed formula. Report it beside
   the exact E9 contrast variance, to show what the paired-difference
   reduction loses (the compendium's Table `tab:invf-opt` gives 22.5
   percent for Hybrid under AR(1)).
2. Validate out of sample. Fix every input from the DGP, and compare
   with a new paired run under the current DGP. Use common seeds
   across designs and half-lives, report paired MCSE, and test the
   screen at several $\rho$ values and sample sizes.
3. Add a decision-relevant open-label-first comparison:
   - designs that begin open label, with varying discontinuation
     timing and crossover length, at matched $N$;
   - the precision cost of the expectancy terms;
   - power under CR2 inference.
4. If carryover attenuation is kept, take it from paper 01's
   Proposition 2 under the decaying correlation. Validate it against
   a DGP that contains it.
5. Reframe the paper as a tutorial-style short communication. Its
   contribution would be a validated, closed-form design screen for
   the interaction under irregular schedules, positioned against Wang
   and Schork (2019) and the optimal-crossover literature.

## 8. Scope of this review

- **Read in full.**
  - `report.Rmd` with all code chunks, `references.bib`, `README.md`,
    `cover-letter.md` and `title-page.md`.
  - `analysis/scripts/nof1-design-efficiency/01-closed-form-validation.R`
    and the stored CSV.
  - The July versions (commit 26b8c16a3) of the dgp-combined driver,
    `simulation-core.R` and `implementations/tidyverse/R/functions.R`.
  - STYLE.md and the relevant parts of NOTATION.md.
- **Read in part.** Paper 01 (closed-form section and archived
  table), paper 08 and paper 10 (title and abstract), and the
  compendium subsections cited above.
- **No prior referee report exists for paper 12.** None was found in
  the manuscript directory, in `docs/reviews/` or in the August
  `pub_review` files, so there are no earlier comments to track. The
  compendium's annotations of paper 12 (2026-10-08) were addressed
  only by inserted notes, not by rerunning the validation.
- **Executed.**
  - Reproduction of every inline statistic from the stored CSV and
    archived summary. All matched to machine precision.
  - Centered and uncentered $R^2$.
  - Null-cell rejection rates.
  - Exact contrast, GLS and sandwich variances under the July DGP.
  - Positive-definiteness of every path covariance under the July
    settings and under the current settings with repair refused.
  - MCSE of the attenuation differences.
  - Verification of the correction-note values.
- **Not done.**
  - No new Monte Carlo simulation.
  - The current-DGP standard deviations are a closed-form inference,
    not simulated.
  - The frequency of the silent no-correlation fallback in `fit_spec`
    was not measured.
  - The replicate-level file was not opened.
  - The PDF was not re-rendered.
  - Paper 02's claim about decay shape was not checked against its
    current text.
  - Heo and Leon, Kunert (1985), and Kunert and Martin were not read
    in the original.
- **Not applicable.**
  - The design is correctly labeled: the Hybrid simulated is weeks 4,
    8, 9, 10, 11, 12, 16 and 20 with four paths, and OL+BDC is the
    two-path schedule. INSPECTED.
  - No dropout was applied, so there is no dependence on
    `R/censordata.R`. INSPECTED.
- **Alternative short versions.** None exist in this directory. No
  `bullets.Rmd` is present.
