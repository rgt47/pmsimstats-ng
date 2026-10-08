# Stress-Testing a Provisional Analysis for the Biomarker-by-Drug Interaction in the Hybrid Design: Pre-Registration {.unlisted .unnumbered}
*2026-10-06 19:34 PDT*

**Author.** pmsimstats team

**Purpose.** This document fixes, before any of the simulations it
describes are run, how the project will choose its analysis for the
biomarker-by-drug interaction in the Hybrid design.

It names a *strawman*: a provisional analysis, put forward to be
attacked. It defines four criteria the analysis must meet: nominal,
robust, powerful and sensitive. It lists the challengers, and a set of
stress conditions designed to break the strawman. It then states the
rule by which the strawman survives or is replaced.

The strawman is chosen on the evidence in `docs/45`, read under the
strict size criterion of Section 3. It is not the primary analysis
recommended there. Section 8 explains the change.

The study follows the ADEMP structure: aims, data-generating
mechanisms, estimands, methods, performance measures (Morris, White
and Crowther [1]).

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md` and `docs/45`. Symbols
marked *local* are defined here.

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $Y_{it}$ | symptom score (code `Sx`); $t = 0$ is baseline | canonical |
| $BR_{it}$, $PB_{it}$, $TV_{it}$ | drug, expectancy and natural-course responses | canonical |
| $B_i$, $b_i$ | biomarker; standardized biomarker | canonical |
| $D_{it}$, $D_{bc,it}$ | binary drug state; exposure-decayed indicator, 0 before any exposure | canonical |
| $t_{1/2}$, $t_{1/2}^{a}$ | true carryover half-life; half-life assumed by the analysis | canonical; local (`docs/45`) |
| $c_{bm}$ | biomarker-BR coupling (the predictive effect) | canonical |
| $c_{bm,TV}$, $c_{bm,PB}$ | biomarker-TV and biomarker-PB correlations (prognostic couplings) | local; $c_{bm,PB}$ is the package's `c.bm.pb` |
| $\beta_{bm}^{TV}$ | mean moderation of the natural course by $b_i$ | canonical (component slope) |
| $\beta_B$, $\beta_{bm:D}$ | biomarker main effect; the interaction, the estimand | local (`docs/41`); canonical |
| $\rho$ | AR(1) serial correlation per week | canonical |
| $\sigma_{\varepsilon,\text{OL}}$, $\sigma_{\varepsilon,\text{bl}}$ | SD of added measurement error, open-label and blinded visits | local |
| $\kappa$ | mean model SE over empirical SD | canonical |
| $\pi$ | power | canonical |
| $N$ | total sample size | canonical |

Table: Notation index and glossary of symbols used in this document

**Glossary.**

- **Strawman.** A provisional analysis put forward to be criticized and
  tested, not a deliberately weak one. In `docs/36` and `docs/44` the
  word named the published-code reference arm; that arm is called the
  *published reference arm* here.
- **Challenger.** An analysis that could replace the strawman.
- **Stress cell.** A condition chosen because it is expected to break at
  least one candidate.
- **Power curve.** Size-adjusted power as a function of the true coupling
  $c_{bm}$, the basis of the sensitivity criterion.
- **Base cell.** The condition from which each stress varies one factor.
- **Prognostic biomarker.** A biomarker associated with the outcome
  regardless of treatment, here through TV or PB.

## 1. Aims

1. Test whether the strawman (Section 4) holds its nominal size in every
   stress cell.
2. Compare its power and sensitivity with every challenger that also
   holds its size.
3. Replace the strawman if, and only if, the rule of Section 7 says so.
4. Report the published reference arm alongside, for comparability with
   Hendrickson et al. [2].

## 2. Scope

- **Design:** Hybrid only. Visits at weeks 4, 8, 9, 10, 11, 12, 16 and 20,
  plus baseline. Four paths with the published 18/18/17/17 allocation at
  $N = 70$, scaled proportionally at other $N$.
- **Estimand:** $\beta_{bm:D}$.
- **Test:** two-sided at $\alpha = 0.05$.

## 3. Criteria (fixed before the run)

**Nominal.** Two conditions must both hold:

- the pooled rejection rate over all null cells lies within 2 Monte
  Carlo SE of 0.05;
- no null cell exceeds 0.075.

Each null cell has 1,000 replicates, so the MCSE is 0.0069 per cell.
There is no "slightly liberal" class. An analysis that fails either
condition fails the criterion.

**Robust.** Nominal in every stress cell of Section 5, each assessed
separately with the cell condition above. A failure in a single stress
cell is recorded with the condition that caused it.

**Powerful.** Power is size-adjusted: each analysis rejects at the 5%
quantile of its own p-values in the matching null cell. The summary is
the maximin shortfall. In each alternative cell, compute each
analysis's shortfall from the best nominal analysis in that cell. The
summary for an analysis is its largest shortfall over the cells.
Unadjusted power is also reported.

**Sensitive.** The test should respond to the strength of the true
interaction: power should rise steadily as the true coupling $c_{bm}$
moves away from zero in either direction, and reach $\alpha$ at zero.
It is assessed on the power curve, size-adjusted, over the grid of
Section 5. Its components are:

- **monotonicity**: no decrease in power between adjacent grid points
  larger than 2 MCSE, on either side of zero;
- **symmetry**: power at $-c_{bm}$ within 2 MCSE of power at $+c_{bm}$;
- **steepness**: the average slope of power against $|c_{bm}|$ between
  0.15 and 0.30;
- **detectable coupling**: the smallest $|c_{bm}|$ at which power
  reaches 0.80, by linear interpolation on the grid. Smaller is more
  sensitive.

The curve is estimated in two constructs. In the base construct (AR(1))
$|c_{bm}|$ runs to 0.45. In the compound-symmetry construct it runs to
0.25, the largest value below that construct's ceiling of 0.256.

Estimand fidelity is also reported, though it is not a selection
criterion: the relative bias of $\hat\beta_{bm:D}$ against the value
the same analysis gives at $t_{1/2} = 0$, when $t_{1/2}^{a}$ is close
to $t_{1/2}$.

## 4. The strawman and the challengers

**Strawman S.** Random intercept plus `corCAR1` residuals, with CR2
standard errors and Satterthwaite degrees of freedom (Bell and
McCaffrey [3]). Baseline is a response row, coded off drug. The model
uses linear time and $D_{bc,it}$ at $t_{1/2}^{a} = 0.5$ week:

```{.r}
nlme::lme(Sx ~ bmc + t + Dbc + bmc:Dbc, random = ~ 1 | ptID,
          correlation = nlme::corCAR1(form = ~ t | ptID))
# then clubSandwich::coef_test(fit, vcov = 'CR2',
#                              test = 'Satterthwaite')
```

It is chosen because it is the most powerful analysis that was nominal,
under the strict criterion, in all eight null cells of `docs/45`. Its
sizes there were 0.040 to 0.064, with a mean of 0.052.

**Challengers.** All are fitted to the same trials.

| ID | Analysis | Reason |
|---|---|---|
| A | unstructured (`us`), Kenward-Roger, baseline row, linear $t$, $D_{bc,it}$ at 0.5 | primary of `docs/45`; most powerful under AR(1) there |
| C1 | S with `varIdent(form = ~ 1 \| phase)`, phase = baseline / open label / blinded; CR2 | gives the parsimonious model the phase heteroscedasticity that likely drives A's advantage |
| C2 | C1 with model-based inference (the `nlme` Wald t-test; `nlme` has no Satterthwaite option for this model) | whether CR2 costs power once the variance model is closer |
| DS | S with $t_{1/2}^{a}$ chosen by AIC over {0, 0.25, 0.5, 1, 2} (0 meaning $D_{it}$); ML for selection, REML refit | data-chosen half-life |
| DA | A with the same AIC rule | the same, for the unstructured model |
| ES, EA | S and A with $\beta_B$ freed at baseline (`+ bmc:bl0`); both are fitted, so the comparison is not chosen after the results | robustness to failure of the constancy assumption |
| F1 | S plus `De + bmc:De`, where `De` is the expectancy (1 open label, 0.5 blinded, 0 at baseline) | remedy for a placebo-prognostic biomarker |
| F2 | S plus `bmc:t` | remedy for a natural-course-prognostic biomarker |
| R | published reference arm: random intercept, model-based, linear $t$, binary $D_{it}$ (`cs-Satt` in `docs/45`) | comparability |

Table: Challenger analyses with identifiers and reasons for inclusion

ES, EA, F1 and F2 change the model, not only the inference. They are
evaluated in every cell, so their cost when nothing is wrong is
measured as well as their benefit when something is.

## 5. Data-generating mechanisms

**Base construct.** Configuration C of `docs/36`:

- separable AR(1) within each component at $\rho = 0.7$ per week;
- graded coupling of the biomarker to BR;
- published Gompertz means and SDs;
- the biomarker coupled only to BR.

This is the construct `docs/44` identifies as the best implemented. The
base cell is:

- $N = 70$, $t_{1/2} = 0.5$ week, no dropout, exponential decay.

**Power-curve cells** (the sensitivity criterion):

| Construct | Grid of $c_{bm}$ |
|---|---|
| base (AR(1), $\rho = 0.7$) | $-0.30$, $0$, $0.075$, $0.15$, $0.225$, $0.30$, $0.375$, $0.45$ |
| compound symmetry, 0.8 at every lag | $-0.20$, $0$, $0.05$, $0.10$, $0.15$, $0.20$, $0.25$ |

Table: Power-curve cells by construct with their $c_{bm}$ grids

The zero cells have 1,000 replicates, all others 500. A negative
$c_{bm}$ means the biomarker is coupled to a smaller drug response.
Under graded coupling the joint distribution exists for negative values
on the same terms as for positive ones, which is checked before the
run. The base-construct cells at 0 and 0.30 double as the null and
alternative of the base cell.

**Stress cells.** Each varies one factor from the base cell. Each has a
null twin ($c_{bm} = 0$, 1,000 replicates) and an alternative at
$c_{bm} = 0.30$ (500 replicates), unless the table says otherwise.

| # | Factor | Level | Expected to break |
|---|---|---|---|
| 1 | correlation | compound symmetry, 0.8 at every lag (`docs/45` CS): covered by the compound-symmetry power curve, whose largest value is 0.25 because 0.30 exceeds that construct's ceiling | pure-AR(1) assumptions; tests S under its weakest structure |
| 2 | correlation | base plus measurement error, $\sigma_{\varepsilon,\text{OL}} = 3$, $\sigma_{\varepsilon,\text{bl}} = 6$ | homoscedastic working models (S, F1, F2) |
| 3 | sample size | $N = 35$ | `us` with Kenward-Roger (A, DA); CR2 degrees of freedom |
| 4 | sample size | $N = 150$ | none expected; shows whether any excess persists with $N$ |
| 5 | carryover | $t_{1/2} = 0$ | $D_{bc,it}$ when there is no carryover |
| 6 | carryover | $t_{1/2} = 1$ week | binary $D_{it}$ (R); short $t_{1/2}^{a}$ |
| 7 | decay form | Weibull, $k = 0.5$ (heavy tail) | the exponential $D_{bc,it}$ |
| 8 | decay form | Weibull, $k = 2$ (light tail) | the exponential $D_{bc,it}$ |
| 9 | dropout | 20% monotone, completely at random | covariance estimation (A) |
| 10 | dropout | 20% monotone, at random given the last observed $Y$ | all; MAR is assumed by every likelihood analysis |
| 11 | prognostic | $c_{bm,TV} = 0.3$ (time-constant natural-course association) | the shared $\beta_B$ of the row coding (S, A, C) |
| 12 | prognostic | $\beta_{bm}^{TV} = 0.3$ (association growing along the TV curve) | every analysis without `bmc:t` |
| 13 | prognostic | $c_{bm,PB} = 0.3$ (expectancy-scaled placebo association) | every analysis without `bmc:De` |

Table: Stress-test factors, levels and the analyses each is expected to break

For the prognostic cells 11 to 13, the alternative twin keeps the
prognostic coupling as well. The question is whether the remedies
recover the interaction, not only whether the null holds.

**Implementation of the new mechanisms** in the `04-power-simulation.R`
construct. The biomarker is index 1 of the joint vector.

- **$c_{bm,TV}$, $c_{bm,PB}$.** Set the biomarker's correlation with every
  post-baseline TV (or PB) element to the stated value. The PB SD is
  already scaled by the expectancy $e_t$, so the implied slope is
  expectancy-scaled. If the matrix is not positive definite the cell
  stops with an error; no repair is permitted.
- **$\beta_{bm}^{TV}$.** After the draw, add
  $\beta_{bm}^{TV}\, b_i\, \sigma_{TV}\, G_{TV}(t)/\max_t G_{TV}(t)$ to $TV_{it}$.
  This is mean moderation that grows along the natural-course curve.
- **Measurement error.** Add independent normal error to $Y_{it}$, with
  SD 3 at open-label visits and 6 at blinded visits; none at baseline.
- **Weibull decay.** Applied to both the BR carryover and the graded
  coupling, as in manuscript 02's decay family. The analysis keeps the
  exponential $D_{bc,it}$.
- **Dropout.** Monotone, starting at a post-baseline visit:
  - MCAR: drawn uniformly over visits 3 to 8 for 20% of participants;
  - MAR: the hazard at each visit increases with the previous observed
    $Y_{i,t-1}$ (logistic, calibrated to 20% overall).

## 6. Methods and performance measures

**Replicates:** 1,000 per null cell and 500 per alternative cell, with
per-chunk seeds as in `08-covariance-study.R`. All analyses are fitted
to the same trials.

**Reported for each analysis and cell:**

- rejection rate with its MCSE;
- size-adjusted power;
- mean estimate, empirical SD, mean model SE and $\kappa$;
- relative bias where defined;
- the distribution of the selected half-life (DS, DA);
- failed fits.

Failed fits count as non-rejections, and their rate is reported. An
analysis with more than 2% failed fits in any cell is flagged.

**Paired comparisons.** Differences in power between two analyses in a
cell are tested with McNemar's test on the paired rejections.

## 7. Decision rule

1. **Screen.** Discard every analysis that fails the nominal criterion in
   any null cell, base or stress, except the prognostic stress cells
   11 to 13.
2. **Prognostic screen.** Cells 11 to 13 are judged separately. An
   analysis that fails one of them is not discarded. It is flagged, and
   recommended only together with the remedy (ES, EA, F1 or F2) that is
   nominal in that cell.
3. **Sensitivity screen.** Flag any survivor whose power curve fails
   the monotonicity or symmetry condition in either construct; a flagged
   analysis may not be selected.
4. **Select.** Among the remaining survivors, select the analysis with
   the smallest maximin shortfall over the alternative cells, power-curve
   cells included. Report the detectable coupling and steepness of
   every survivor, and break ties within 0.02 of maximin shortfall by
   the smaller detectable coupling in the base construct.
5. **Strawman rule.** S is retained unless a survivor beats it on the
   maximin shortfall by more than 0.05, with the McNemar difference
   significant at 0.05 in at least half the alternative cells. Ties go
   to the analysis with fewer covariance parameters.
6. **Replacement.** If S is replaced, the winner is reported as the
   recommended analysis. It is not stress-tested again in this study,
   because doing so would make the selection adaptive.

These rules are fixed by this document. Any change after results are
seen will be reported as a deviation, with its reason.

## 8. Relation to earlier documents

- **`docs/45` recommended A as primary.** Under the strict criterion of
  Section 3, A's mean size of 0.058 fails the nominal criterion (the
  2-MCSE band is 0.043 to 0.057). The "slightly liberal" class of
  `docs/45` was defined after the results were seen. This document
  removes it, and S becomes the strawman.
- **The parsimonious challengers are given the feature A exploited.** In
  `docs/45` they were not. C1 adds the phase-specific variance that was
  the likely source of A's power advantage.
- **The constancy assumption is now tested.** `docs/45` could not test
  it, because every construct satisfied it. Cells 11 to 13 do.

## 9. Compute and execution

**Size.** 15 power-curve cells (8 base, 7 compound symmetry) plus 13
stress cells × 2, about 40 cells once the curve cells shared with the
stress set are counted once. About 12 analyses per replicate, with
about 24,000 replicates in total. A, DA and
EA use unstructured Kenward-Roger fits and dominate the cost.

**Time.** The pilot (20 replicates per cell, 780 in all, 16.7 minutes
on 8 cores, no failed fits) puts the full run at about 9 to 10 hours.
An earlier estimate of four to five days, scaled from the
exposure-coding pass of `docs/45`, was too high. The run is local and
detached, with each chunk of 25 replicates checkpointed.

**Order of work.**

1. Implement the mechanisms of Section 5 and the analyses of Section 4
   in a new driver, `11-strawman-stress-test.R`. Add unit checks:
   - the positive-definiteness check on every new matrix;
   - the prognostic couplings reproduce their target correlations in a
     large draw;
   - the dropout rates match their targets;
   - the AIC rule reproduces the fixed $t_{1/2}^{a}$ analysis when it
     selects it.
2. A pilot of 20 replicates per cell, to check every analysis fits and to
   time the run.
3. The full run.
4. A results document (`docs/47`), reporting against this
   pre-registration, deviations included.

## 10. Limitations known in advance

- One design (Hybrid) and one construct family (Gompertz components with
  separable AR(1) or CS).
- The prognostic couplings and the measurement-error SDs are chosen
  values, not estimates from trial data.
- Stress factors are varied one at a time, so interactions between
  stresses (for example small $N$ with dropout) are not tested.
- The strawman's own analysis half-life (0.5 week) and the AIC grid are
  fixed here. Other values are not explored.

## Amendment 1 (2026-10-07, after the results were seen)

Two changes to the decision rule, made after the full run. They are
deviations and are reported as such. The as-registered rule remains
reproducible with the defaults of `12-strawman-stress-decision.R`.

1. **The symmetry screen is dropped.** Without prognostic coupling,
   power is exactly symmetric in the sign of $c_{bm}$ by construction
   (`docs/48`, Section 4). Its 9 failures in 22 tests were Monte Carlo
   error, confirmed by a fresh-seed rerun. The monotonicity screen is
   kept.
2. **The strawman selects its half-life by AIC.** It is analysis DS:
   S with $t_{1/2}^{a}$ chosen over {0, 0.25, 0.5, 1, 2}, ML for
   selection and REML refit. All other rules are unchanged.

Run with `--strawman DS --no-symmetry`. Outputs carry the suffix
`-ds-nosym`.

**Outcome.**

- **Survivors** of the nominal and monotonicity screens: S, C1, C2,
  DS, ES, EA, F1 and F2.
- **Prognostic screen:** S, C1, C2, DS and F2 are flagged. ES, EA and
  F1 are not.
- **The strawman is replaced.**
  - C1 beats DS on maximin shortfall by 0.148 (0.034 against 0.182),
    and is significantly ahead by McNemar's test in 19 of 22
    alternative cells.
  - C2 beats it by 0.142, ahead in 20 of 22.
  - S itself beats DS by only 0.038. AIC selection costs a little
    power against the fixed half-life of 0.5 week.
- **The maximin winner is C1, but it is flagged** (null rejection 0.382
  under a constant natural-course association). By rule 2 it may be
  recommended only together with a remedy that is nominal in that
  cell, and no such C1 variant was pre-registered.
- **F1 is the best of the unflagged analyses** (maximin shortfall
  0.242, against 0.278 for ES and 0.318 for EA).

The combination of C1 with F1's remedy is being run as a post hoc arm
(`11-strawman-stress-test.R --extra`). It can motivate a
recommendation but cannot replace this outcome without a confirmatory
run.

## 11. References

1. Morris TP, White IR, Crowther MJ. Using simulation studies to
   evaluate statistical methods. *Statistics in Medicine* 2019;
   38(11):2074-2102. doi:10.1002/sim.8086.
2. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020; 2:13. doi:10.3389/fdgth.2020.00013.
3. Bell RM, McCaffrey DF. Bias reduction in standard errors for linear
   regression with multi-stage samples. *Survey Methodology* 2002;
   28(2):169-181.
4. pmsimstats team. `docs/36`, `docs/37`, `docs/41`, `docs/44`,
   `docs/45`.
