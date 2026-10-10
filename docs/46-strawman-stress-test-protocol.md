# Stress Test of Candidate Analyses for the Biomarker-by-Drug Interaction in the Hybrid Design: Working Protocol {.unlisted .unnumbered}
*2026-10-06 19:34 PDT*

**Author.** pmsimstats team

**Purpose.** This working protocol sets out how the project chooses its
analysis for the biomarker-by-drug interaction in the Hybrid design.
The study is exploratory: the protocol was revised as results came in
(Amendments 1 and 2), and its findings should be confirmed in an
independent study.

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
| $t_{1/2}$, $t_{1/2}^{a}$ | true carryover half-life; half-life assumed by the analysis | canonical |
| $c_{bm}$ | biomarker-BR coupling (the predictive effect) | canonical |
| $c_{bm,TV}$, $c_{bm,PB}$ | biomarker-TV and biomarker-PB correlations (prognostic couplings) | canonical; $c_{bm,PB}$ is the package's `c.bm.pb` |
| $\beta_{bm}^{TV}$ | mean moderation of the natural course by $b_i$ | canonical (component slope) |
| $\beta_B$, $\beta_{bm:D}$ | biomarker main effect; the interaction, the estimand | canonical |
| $\rho$ | AR(1) serial correlation per week | canonical |
| $\rho_{cs}$ | common correlation of the compound-symmetry construct, the same at every lag | canonical; row added 2026-10-09 |
| $\sigma_{\varepsilon,\text{OL}}$, $\sigma_{\varepsilon,\text{blind}}$ | SD of added measurement error, open-label and blinded visits | local |
| $\kappa$ | mean model SE over empirical SD | canonical |
| $\pi$ | power | canonical |
| $N$ | total sample size | canonical |
| $\alpha$ | nominal level, 0.05 | canonical; row added 2026-10-07 |
| $\eta_t$ | expectancy at visit $t$ (code `De`): 0 baseline, 1 open label, 0.5 blinded | canonical (NOTATION.md) |
| $\sigma_{TV}$, $G_{TV}(t)$ | SD and Gompertz mean curve of the natural-course component | canonical; row added 2026-10-07 |
| $k$ | Weibull shape of the decay family | canonical; row added 2026-10-07 |
| $c_{bm,BL}$ | biomarker-baseline correlation (baseline-severity prognostic coupling; Amendment 2) | canonical; row added 2026-10-07 |

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

## 3. Criteria

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

(Clarified 2026-10-08: the term "maximin shortfall" is kept because the
decision outputs use it, but the rule is a minimax-regret rule: each
analysis is scored by its largest shortfall, and the analysis with the
smallest such score is preferred. The best analysis in each cell is
taken over every analysis passing the size and sensitivity screens,
including those flagged by the prognostic screen, such as C1 and C2,
as `12-strawman-stress-decision.R` implements it; the shortfalls of the
unflagged analyses are therefore measured against the phase-variance
models. See `analysis/report/compendium-mathematics.tex`, Section
"Size, power screens and the selection criteria of the stress test".)

**Sensitive.** The test should respond to the strength of the true
interaction: power should rise steadily as the true coupling $c_{bm}$
moves away from zero in either direction, and reach $\alpha$ at zero.
It is assessed on the power curve, size-adjusted, over the grid of
Section 5. Its components are:

- **monotonicity**: no decrease in power between adjacent grid points
  larger than 2 MCSE, on either side of zero (corrected 2026-10-08:
  `12-strawman-stress-decision.R` checks the positive side of the grid
  only, from zero upward; the negative side has the single point
  $-0.30$, or $-0.20$ in the compound-symmetry construct, so there is
  no adjacent pair to check there);
- **symmetry**: power at $-c_{bm}$ within 2 MCSE of power at $+c_{bm}$;
- **steepness**: the average slope of power against $|c_{bm}|$ between
  0.15 and 0.30;
- **detectable coupling**: the smallest $|c_{bm}|$ at which power
  reaches 0.80, by linear interpolation on the grid. Smaller is more
  sensitive.

The curve is estimated in two constructs. In the base construct
(configuration C, separable AR(1) at $\rho = 0.7$; Section 5)
$|c_{bm}|$ runs to 0.45. In the compound-symmetry construct it runs to
0.25, the largest value below that construct's ceiling of 0.256.
(Corrected 2026-10-08: 0.256 is the compound-symmetry Hybrid ceiling at
$t_{1/2} = 0$; at the curve's half-life of 0.5 week the graded ceiling
is 0.274. The grid, and the exclusion of 0.30, are unaffected. See
`analysis/report/compendium-mathematics.tex`, Section "The stress-test
construct and its fourteen analyses".)

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
| C1 | S with `nlme::varIdent()` residual variances by `phase`, phase = baseline / open label / blinded; CR2 | gives the parsimonious model the phase heteroscedasticity that likely drives A's advantage |
| C2 | C1 with model-based inference (the `nlme` Wald t-test; `nlme` has no Satterthwaite option for this model) | whether CR2 costs power once the variance model is closer |
| DS | S with $t_{1/2}^{a}$ chosen by AIC over {0, 0.25, 0.5, 1, 2} (0 meaning $D_{it}$); ML for selection, REML refit | data-chosen half-life |
| DA | A with the same AIC rule | the same, for the unstructured model |
| ES, EA | S and A with $\beta_B$ freed at baseline (`+ bmc:bl0`); both are fitted, so the comparison is not chosen after the results | robustness to failure of the constancy assumption |
| F1 | S plus `De + bmc:De`, where `De` is the expectancy (1 open label, 0.5 blinded, 0 at baseline) | remedy for a placebo-prognostic biomarker |
| F2 | S plus `bmc:t` | remedy for a natural-course-prognostic biomarker |
| R | published reference arm: random intercept, model-based, linear $t$, binary $D_{it}$ (`cs-Satt` in `docs/45`) | comparability |

Table: Challenger analyses with identifiers and reasons for inclusion

(Corrected 2026-10-08: R is described above as a random intercept with
model-based inference. The driver `11-strawman-stress-test.R` fits it
as the marginal compound-symmetry structure of `mmrm`,
`cs(visit | ptID)`, with model-based standard errors and Satterthwaite
degrees of freedom, which is the `cs-Satt` arm of `docs/45`. A
compound-symmetry covariance equals a random intercept with
independent errors when the common correlation is nonnegative, so the
model is the same; the inference is Satterthwaite, not the `nlme`
containment degrees of freedom. See
`analysis/report/compendium-mathematics.tex`, Section "The stress-test
construct and its fourteen analyses".)

ES, EA, F1 and F2 change the model, not only the inference. They are
evaluated in every cell, so their cost when nothing is wrong is
measured as well as their benefit when something is.

## 5. Data-generating mechanisms

**Base construct.** Configuration C of `docs/36`:

- separable AR(1) within each component at $\rho = 0.7$ per week;
- graded coupling of the biomarker to BR;
- published Gompertz means and SDs;
- the biomarker coupled only to BR.

This is the covariance and coupling construct of the refined process
(paper 01; `docs/44`, whose near-term target adds anchored carryover
to it). It is not the
refined process in full: the driver keeps the published recursive
carryover on the BR mean and the fixed published allocation (inspected
in `11-strawman-stress-test.R`, 2026-10-08), and it stops on a matrix
that is not positive definite, as the refined process does. The base
cell is:

- $N = 70$, $t_{1/2} = 0.5$ week, no dropout, exponential decay.

**Power-curve cells** (the sensitivity criterion):

| Construct | Grid of $c_{bm}$ |
|---|---|
| base (AR(1), $\rho = 0.7$) | $-0.30$, $0$, $0.075$, $0.15$, $0.225$, $0.30$, $0.375$, $0.45$ |
| compound symmetry, $\rho_{cs} = 0.8$ | $-0.20$, $0$, $0.05$, $0.10$, $0.15$, $0.20$, $0.25$ |

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
| 1 | correlation | compound symmetry, $\rho_{cs} = 0.8$ (`docs/45` CS): covered by the compound-symmetry power curve, whose largest value is 0.25 because 0.30 exceeds that construct's ceiling | pure-AR(1) assumptions; tests S under its weakest structure |
| 2 | correlation | base plus measurement error, $\sigma_{\varepsilon,\text{OL}} = 3$, $\sigma_{\varepsilon,\text{blind}} = 6$ | homoscedastic working models (S, F1, F2) |
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
  already scaled by the expectancy $\eta_t$, so the implied slope is
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

Changes to these rules made after results were seen are recorded in
the amendments below as revisions of the working protocol, each with
its reason.

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
4. A results document (`docs/47`), reporting against this working
   protocol, its revisions included.

## 10. Limitations

- One design (Hybrid) and one construct family (Gompertz components with
  separable AR(1) or CS).
- The prognostic couplings and the measurement-error SDs are chosen
  values, not estimates from trial data.
- Stress factors are varied one at a time, so interactions between
  stresses (for example small $N$ with dropout) are not tested.
- The strawman's own analysis half-life (0.5 week) and the AIC grid are
  fixed here. Other values are not explored.

## Amendment 1 (2026-10-07, after the results were seen)

Two revisions of the working protocol's decision rule, made after the
full run. The original rule remains reproducible with the defaults of `12-strawman-stress-decision.R`.

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
  cell, and no such C1 variant was in the original protocol.
- **F1 is the best of the unflagged analyses** (maximin shortfall
  0.242, against 0.278 for ES and 0.318 for EA).

The combination of C1 with F1's remedy is being run as a post hoc arm
(`11-strawman-stress-test.R --extra`). It can motivate a
recommendation but cannot replace this outcome without an independent
run.

## Amendment 2 (2026-10-07, after the results were seen)

This amendment records four further revisions of the working protocol,
all post hoc, and the outcome of the decision rule rerun with every
analysis. The original protocol text above, and Amendment 1, are unchanged in substance.
Two edits to them are typographical: the C1 row of the challenger table
showed a stray backslash in its code span and was reworded to name
`nlme::varIdent()` by `phase`; and the expectancy symbol $e_t$ was
renamed $\eta_t$, the canonical symbol of `analysis/report/NOTATION.md`,
with no change of meaning. (Notation pass 2026-10-09: the subscript of
the blinded-visit measurement-error SD was renamed from bl to blind,
$\sigma_{\varepsilon,\text{blind}}$, so that it cannot be read as BL,
the baseline level of $c_{bm,BL}$; also typographical. Notation,
2026-10-09: the compound-symmetry correlation, previously written "0.8
at every lag", is written $\rho_{cs} = 0.8$, the symbol of
`analysis/report/NOTATION.md`; also typographical.)

**1. Two baseline-severity cells (40 and 41).** The biomarker is
correlated $c_{bm,BL} = 0.3$ with the baseline level $\mathrm{BL}_i$,
and with no response component: a null twin (cell 40, 1,000 replicates)
and an alternative at $c_{bm} = 0.30$ (cell 41, 500). Such an
association is the same at every visit, so the constancy assumption
holds. The cells are judged as prognostic cells under rule 2, alongside
cells 11 to 13. The study now has 41 cells, 15 null and 26 alternative,
and 28,000 simulated trials.

**2. Three analyses (G1, G2, F1A).**

- **G1:** the F1 mean model with phase-specific residual variances
  (`nlme::varIdent()` by phase, as in C1), CR2 inference;
- **G2:** the same model with model-based inference, as in C2;
- **F1A:** F1 with $t_{1/2}^{a}$ chosen by AIC over {0, 0.25, 0.5, 1,
  2}, ML for selection and REML refit, as in DS.

They are fitted to the same trials and judged by the same rule as the
original analyses. Being post hoc, none can replace the outcome
of Amendment 1 without an independent run.

**3. Failed fits.** Section 6 states that failed fits count as
non-rejections. The summaries of `11-strawman-stress-test.R` and the
decision script instead exclude them from the denominator. Failed fits
number 0 to 12 per analysis in 28,000 trials, at most 6 of 1,000 (0.6%)
in any cell (G1 and G2), below the 2% flag. Counting them as
non-rejections would lower a cell's rejection rate by at most 0.0025
(verified, `summary-full*.csv`). This is far smaller than the margin
of any size screen outcome, so no outcome depends on the convention
(inferred from the margins).

**4. Outcome with all 14 analyses** (verified; the added analyses come
from `11-strawman-stress-test.R --extra` and `--extra2`, and
`12-strawman-stress-decision.R --strawman DS --no-symmetry
--with-extras` writes `decision-*-ds-nosym-all14.csv`).

- **Size screen.** A (0.057 pooled; 0.087 at $N = 35$), DA (0.065),
  R (0.069) and F1A fail. F1A rejects 0.057 pooled, outside the
  2-MCSE band of 0.046 to 0.054, and 0.073 in the $t_{1/2} = 1$ cell.
  AIC selection therefore does not hold its size once the expectancy
  terms are in the model.
- **Survivors** of the size and sensitivity screens: S, C1, C2, DS, ES,
  EA, F1, F2, G1 and G2.
- **Baseline-severity cells.** In the null cell every analysis is
  within the cell limit, 0.040 (F1) to 0.075 (R), as the constancy
  argument predicts; S rejects 0.050 and A 0.051.
- **Prognostic screen.** The analyses nominal in all four prognostic
  null cells are ES, EA, F1, G1 and G2 (G1 0.042 to 0.055, G2 0.043 to
  0.062). S, C1, C2, DS and F2 remain flagged.
- **Maximin shortfall** over the non-prognostic alternative cells,
  among the prognostic-robust survivors: G2 0.232, F1 0.242, G1 0.252,
  ES 0.278, EA 0.318. C1 has the smallest value of all survivors
  (0.034) but is flagged.
- **Selection.** G2 and F1 tie within 0.02. The tie-break of rule 4,
  the detectable coupling in the base construct, is undefined, because
  neither reaches power 0.80 on the grid. By the tie rule of rule 5
  (fewer covariance parameters), F1 is selected. The outcome of
  Amendment 1 is unchanged.
- **Over all 26 alternative cells,** prognostic cells included, the
  order among the robust survivors is G1 0.308, G2 0.312, F1 0.358,
  EA 0.528, ES 0.544. On this reading the phase-variance versions lead
  F1 by about 0.05, mainly through the natural-course prognostic
  alternatives. Because the protocol does not say whether the
  prognostic alternative cells count (`docs/48`, Section 3), and G1 and
  G2 are post hoc, this is recorded and does not change the selection.

Size-adjusted power of the added analyses in the two reference cells:
at $c_{bm} = 0.30$ in the AR(1) construct, F1 0.364, G1 0.354, G2 0.380
and F1A 0.338; at $c_{bm} = 0.25$ under compound symmetry, F1 0.522,
G1 0.546, G2 0.528 and F1A 0.446 (verified,
`decision-power-ds-nosym-all14.csv`).

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
