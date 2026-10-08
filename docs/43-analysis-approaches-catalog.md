# Analysis Approaches for the Biomarker-by-Drug Interaction in Aggregated N-of-1 Trials: A Catalog {.unlisted .unnumbered}
*2026-10-04 16:41 PDT*

**Author.** pmsimstats team

**Purpose.** This paper lists every analysis approach for the
biomarker-by-drug interaction that has been used, evaluated or proposed
in this project. It covers:

- the published analysis of Hendrickson et al. (`orig`, commit
  `58b32a9`) and its options;
- the nine specifications of manuscript 02;
- the test procedures of manuscripts 08 and 10, and the analysis
  strategies of manuscript 06;
- the summary-measures statistics of `docs/37`;
- the alternatives discussed in `docs/40`, `docs/41` and `docs/42`;
- the analyses in the covariance study (`08-covariance-study.R`),
  including its follow-up run.

The approaches are organized along the dimensions on which an analysis
can vary (Section 2). The named approaches are then listed with their
source, their status, and what is known about the size of their test
(Sections 3 and 4). Section 5 summarizes which approaches are known to
hold their size, and Section 6 lists the gaps.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation and glossary

Notation follows `analysis/report/NOTATION.md`.

- **$\beta_{bm:D}$.** The biomarker-by-drug interaction, the estimand.
- **`Db`, `Dbc`.** Binary drug state, $D_{it}$; exposure-decayed drug
  indicator, $D_{bc,it}$ (1 on drug, $e^{-\lambda t_{sd}}$ off drug).
- **Size.** Type I error of the interaction test; nominal 0.05. An
  approach is *correctly sized* in a cell if its null rejection rate is
  within about 2 Monte Carlo standard errors of 0.05.
- **$\kappa$.** Mean model-based standard error over the empirical
  standard deviation of the estimate. Below 1 means anticonservative.
- **CS, AR(1).** Compound symmetry; first-order autoregressive
  correlation.
- **CR2.** Bias-reduced cluster-robust (sandwich) standard error, with
  Satterthwaite (Bell-McCaffrey) degrees of freedom.
- **KR.** Kenward-Roger standard-error and degrees-of-freedom
  correction.
- **MD.** Mancl-DeRouen bias-corrected sandwich for GEE.
- **cLDA.** Constrained longitudinal data analysis: baseline kept in the
  response with a common mean.
- **Construct.** The data-generating process a result was obtained
  under. *Published CS* is the 58b32a9 construct (`docs/36`);
  *configuration C* is separable AR(1) with graded coupling; *package*
  is the current package construct used by manuscripts 02, 06, 08 and
  10.

## 1. Summary

1. **An analysis is a combination of about a dozen choices.** They are:
   the estimator family, the exposure coding, the time adjustment, the
   expectancy term, the baseline handling, the random effects, the
   working covariance, the standard error and degrees of freedom, the
   path term, the biomarker coding, the target coefficient, and the
   design-specific branches (Section 2).
2. **About 60 distinct approaches** have been used, evaluated or
   proposed (Section 3). Most vary one or two choices against a
   reference.
3. **Correct size depends mainly on the working covariance and the
   small-sample correction**, not on the exposure coding, the time
   adjustment or the baseline coding. Each of those moves the Type I
   error by at most about 0.02. The covariance and its correction move
   it by up to 0.10 (Section 4).
4. **Known to hold their size in every construct tested:**
   - an unstructured covariance with Kenward-Roger, baseline as a row
     (covariance study: 0.046 to 0.072);
   - CR2 standard errors with a random intercept and `corCAR1` (papers
     02 and 10: 0.045 to 0.049);
   - GEE with the Mancl-DeRouen correction (papers 08 and 10).

   **Known to fail somewhere:**
   - the published random intercept under AR(1) data (Type I error
     0.14 to 0.16 in CO);
   - any unstructured covariance with Satterthwaite only (0.072 to
     0.104);
   - naive GEE (1.5 to 1.9 times nominal);
   - flexible time adjustments in Hybrid under carryover (about 0.07,
     under every covariance tested).
5. **Open:** the Kenward-Roger and CR2 versions of the baseline-covariate
   coding, CR2 with structured covariances under AR(1) data, and the
   mechanism behind the row coding's power under AR(1). All are in the
   follow-up run now in progress (Section 6).

## 2. Dimensions of an analysis

Each dimension lists its options, where each was used or proposed, and
its current status. "Evaluated" means a simulation result exists; the
source column gives the paper.

### 2.1 Estimator family

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| Linear mixed model, conditional (`lmer`, `nlme::lme`) | orig, 02, 06, 08, 10, covariance study | evaluated |
| Marginal model with structured covariance (MMRM, `mmrm`, `gls`) | covariance study; `docs/41`, `docs/42` | evaluated |
| Summary measures: per-participant on-minus-off contrast, regressed on the biomarker (E9) | `docs/37`; 02 (G5); 08 (continuous-biomarker contrast) | evaluated |
| Strict repeated-measures ANOVA, dichotomized biomarker | `docs/26`; 08 | evaluated (08) |
| GEE with sandwich variance | 08; 10 | evaluated |
| Two-stage approaches other than E9 | 08 (literature review) | not evaluated |
| Nonlinear component-matched model (Gompertz BR, PB, TV) | `docs/31`; 06 | partly evaluated (06, recovery study) |
| Bayesian hierarchical N-of-1 model | 02 (named as out of scope) | not evaluated |
| Randomization or permutation test | 08 (literature) | not evaluated |
| Latent-class mixture model | manuscript 03 | evaluated (different estimand) |

Table: Estimator family options with source and evaluation status

### 2.2 Drug-exposure coding

| Option | Code | Source | Status |
|---------------------------------|------|------------------------|----------------|
| Binary on-drug indicator `Db` | E1 | orig; 02 G1 | evaluated |
| Exposure-decayed `Dbc`, assumed half-life | E2 | 02 G3; implementations/original | evaluated |
| Binary plus lagged just-off-drug term $L_{it}$ | E3 | 02 G2 (Jones-Kenward) | evaluated: no gain over E1 |
| Binary plus linear time-since-discontinuation | E4 | implementations/original (`simplecarryover`) | retired in 02 |
| Lag term crossed with biomarker ($L$, $\text{bm}{:}L$) | E5 | 02, earlier S7 | retired: worse than E1 |
| $t_{sd}$ crossed with biomarker | E6 | 02, earlier S7 | retired: worse than E1 |
| `Dbc` with half-life chosen by AIC | E7 | 02 G4 | evaluated (one cell) |
| `Dbc` with alternative decay forms (Weibull, linear, power) | | 02 sensitivity blocks; NOTATION | evaluated as misspecification |

Table: Drug-exposure coding options E1 to E7 with source and status

### 2.3 Time adjustment

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| None | E9 (`docs/37`) | evaluated |
| Linear `t` | orig; 02; all mixed models above | evaluated |
| Linear `t` plus expectancy `De` | orig option (`useDE`); `docs/40` | population-mean limit only |
| Natural spline in `t`, 3 df | `docs/40`, `docs/42`; covariance study | evaluated |
| Visit as factor | `docs/40`, `docs/42`; covariance study | evaluated |
| Phase indicators, phase by drug (attribution) | 06 (phase-augmented) | evaluated (06): no bias reduction, less power |
| Random slope on `t` | orig option (`t_random_slope`) | not evaluated |
| Path-specific intercepts (E9 stage 2) | `docs/37` (stratified E9) | closed form only |

Table: Time adjustment options with source and evaluation status

### 2.4 Expectancy

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| Not modeled | orig as published | evaluated |
| Design expectancy `De` as a fixed effect | orig option | population-mean limit only (`docs/40`) |
| Absorbed by visit effects | `docs/40`, `docs/42` | evaluated |
| Measured belief covariate $\eta$ | 06 (balanced-placebo design) | evaluated (06) |

Table: Expectancy modeling options with source and evaluation status

### 2.5 Baseline handling

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| Ninth response row, coded off drug | orig | evaluated |
| Excluded, no covariate | 02 (all nine) | evaluated (02) |
| Covariate (ANCOVA) | `docs/41`; covariance study | evaluated with `cs` and `us` (Satterthwaite) |
| cLDA: row with common baseline mean, unstructured covariance | `docs/41` | equivalent to row plus `us` with visit effects |
| Change from baseline | `docs/41` | identical to ANCOVA in the construct |
| Baseline-by-drug moderator | `docs/41` | not evaluated (zero in the construct) |
| Row, with a separate biomarker effect at baseline | follow-up run | running |

Table: Baseline handling options with source and evaluation status

### 2.6 Random effects

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| None (marginal covariance only) | MMRM; 10 (M0, M2) | evaluated |
| Random intercept | orig; 02; 10 (M1, M3) | evaluated |
| Random intercept and slope on `t` | orig option | not evaluated |
| Random intercept, post-baseline shift and expectancy coefficient, phase-specific residual variance | this discussion (implied by the published construct) | not evaluated |

Table: Random-effects structures with source and evaluation status

### 2.7 Working covariance of the residuals

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| Independent | orig; 10 (M1) | evaluated |
| `corCAR1` (continuous-time AR(1)) | 02; 06; 08; 10 (M2, M3) | evaluated |
| Compound symmetry `cs` | covariance study (equals orig's random intercept) | evaluated |
| Heterogeneous CS `csh` | covariance study | evaluated |
| Unstructured `us` | covariance study; 10 | evaluated |
| Spatial exponential `sp_exp` | covariance study | evaluated |
| AR(1), Toeplitz, ante-dependence | `mmrm` options | not evaluated (assume equal spacing) |
| Phase-specific variance (`varIdent`) | this discussion | not evaluated |
| GEE working correlation: independence, exchangeable, AR(1) | 08; 10 (AR(1)) | evaluated |

Table: Working residual covariance options with source and evaluation status

### 2.8 Standard error and degrees of freedom

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| Model-based SE, Satterthwaite df (`lmerTest`) | orig | evaluated |
| Model-based SE, containment df (`nlme`) | 02; 06; 08 | evaluated |
| Model-based SE, Satterthwaite df (`mmrm`) | covariance study | evaluated |
| Kenward-Roger | covariance study (`us`); 10 | evaluated for `us`; follow-up adds `cs`, `csh` and baseline covariate |
| CR2 sandwich, Satterthwaite df (`clubSandwich`) | 02 (G6 to G9); 10 | evaluated |
| CR2 sandwich, Satterthwaite df (`mmrm`) | follow-up run | running |
| CR0, CR1, CR3 sandwiches | 10 (discussed) | not evaluated |
| Naive GEE sandwich | 08; 10 | evaluated |
| Mancl-DeRouen GEE sandwich | 08; 10 | evaluated |
| OLS t, N − 2 df | E9 | evaluated |

Table: Standard error and degrees-of-freedom options with source and status

### 2.9 Path or sequence term

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| None | orig and most | evaluated |
| Fixed path effect | `docs/37`, single-path study | evaluated: no effect on `lmer` power |
| Path-specific intercepts in E9 | `docs/37` | closed form only |

Table: Path or sequence term options with source and evaluation status

### 2.10 Biomarker coding and target

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| Uncentered biomarker | orig; 02 | evaluated |
| Centered biomarker | `docs/40`; covariance study | evaluated (target unchanged) |
| Dichotomized biomarker | 08 (strict RM-ANOVA) | evaluated: costs 0.08 to 0.21 power |
| Target `bm:Db` | orig; 02 G1, G2 | evaluated |
| Target `bm:Dbc` | 02 G3, G4 | evaluated |
| Component-specific slopes | 06 | evaluated |
| Joint 2-df test (E5, E6) | 02 | retired |

Table: Biomarker coding and target coefficient options with source and status

### 2.11 Design-specific branches

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| Biomarker-by-time test when no participant varies in drug state (OL design): `Sx ~ bm + t + bm*t` | orig | evaluated in the published results; confounded with natural course and expectancy |

Table: Design-specific analysis branch for the OL design

## 3. Named approaches

Each row is a complete specification that has been fitted or proposed.
Status is *evaluated*, *running* or *proposed*; the evidence column
gives the construct and the key size result where one exists.

### 3.1 The published analysis and its options (orig, 58b32a9)

| ID | Specification | Status | Evidence on size |
|------|------------------------------------------|----------|----------------------------------------|
| O1 | `lmer(Sx ~ bm + Db + t + bm*Db + (1|ptID))`, baseline row, Satterthwaite (published default) | evaluated | published CS: 0.049 to 0.058 (Hybrid), 0.04 to 0.05 (CO). Published results: 0.01 to 0.09. AR(1) data: 0.14 to 0.23 in CO (`docs/36`) |
| O2 | O1 plus `De` (`useDE = TRUE`) | proposed | dropped by Hendrickson et al. for collinearity; correlation with `Db` 0.59 (`docs/40`) |
| O3 | O1 with `(1 + t | ptID)` (`t_random_slope = TRUE`) | proposed | |
| O4 | OL branch: `Sx ~ bm + t + bm*t + (1|ptID)` | evaluated | confounded test; not a drug-by-biomarker interaction |
| O5 | Later vendored implementation: `nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc)`, random intercept, `corCAR1` | evaluated as 02 G3 | see 02 G3 |

Table: Published analysis O1 and its options O2 to O5, with evidence on size

### 3.2 Manuscript 02 (carryover specifications; package construct)

All use `nlme::lme` with a random intercept and `corCAR1`, no baseline
row, linear `t`.

| ID | Specification | Status | Evidence on size |
|------|------------------------------------------|----------|----------------------------------------|
| G1 | `Db`, target `bm:Db` | evaluated | mildly conservative (mean null 0.036 to 0.039) |
| G2 | G1 plus lag $L$ | evaluated | as G1; no power gain over G1 |
| G3 | `Dbc`, target `bm:Dbc` | evaluated | as G1 |
| G4 | G3, half-life chosen by AIC | evaluated (one cell) | nominal despite selection |
| G5 | Paired difference (E9) | evaluated | |
| G6 to G9 | G1, G2, G3, G4 with CR2 | evaluated | CR2 restores nominal size (0.045 to 0.049 in 10) |

Table: Manuscript 02 carryover specifications G1 to G9 with evidence on size

### 3.3 Manuscripts 08 and 10 (test procedures and calibration; package construct)

| ID | Specification | Status | Evidence on size |
|------|------------------------------------------|----------|----------------------------------------|
| T1 | Strict RM-ANOVA, dichotomized biomarker | evaluated (08) | dichotomization costs 0.08 to 0.21 power |
| T2 | Continuous-biomarker within-participant contrast | evaluated (08) | calibrated |
| T3 | `lme`, random intercept, `corCAR1` (production) | evaluated (08, 10 M3) | conservative, 0.029 to 0.033 (10) |
| T4 | GEE, naive sandwich | evaluated (08, 10) | anticonservative, 1.5 to 1.9 times nominal (08); 0.067 to 0.069 (10) |
| T5 | GEE, Mancl-DeRouen | evaluated (08, 10) | restores nominal size |
| M0 | OLS, no within-participant correlation | evaluated (10) | not anticonservative |
| M1 | Random intercept only | evaluated (10) | well calibrated in 10's construct |
| M2 | `corCAR1` only | evaluated (10) | intermediate |
| M3 | Random intercept plus `corCAR1` | evaluated (10) | conservative; CR2 fixes it; KR inert on it |
| M-us | Unstructured, KR against Satterthwaite | evaluated (10) | Satterthwaite 0.076, KR 0.059 |

Table: Manuscript 08 and 10 test procedures with evidence on size

### 3.4 Manuscript 06 (component decomposition)

| ID | Specification | Status | Evidence |
|------|------------------------------------------|----------|----------------------------------------|
| C1 | One-component (`bm:Dbc`) | evaluated | unbiased if the biomarker is unrelated to PB and TV; size 0.031 to 0.056 |
| C2 | Phase-augmented: phase, phase by drug, phase by biomarker by drug | evaluated | no bias reduction; less power; size 0.038 to 0.049 |
| C3 | Belief-covariate decomposition | evaluated | recovers the pharmacological slope only on a balanced-placebo design |
| C4 | Blinded-stratum contrast; Gompertz-basis decomposition | evaluated | fail on the Hybrid design |

Table: Manuscript 06 component-decomposition strategies C1 to C4 with evidence

### 3.5 Summary measures (`docs/37`)

| ID | Specification | Status | Evidence |
|------|------------------------------------------|----------|----------------------------------------|
| S1 | E9: post-baseline on-minus-off contrast, OLS on biomarker | evaluated (closed form and Monte Carlo) | exact closed form; validated over 566 cell-variants |
| S2 | E9b: baseline counted as off drug | evaluated | mirrors O1's baseline coding; 15% to 16% noisier than E9 without carryover (`docs/41`) |
| S3 | Stratified E9: path-specific intercepts | closed form only | removes the CO period-effect loss |

Table: Summary-measures statistics S1 to S3 with status and evidence

### 3.6 This discussion (`docs/40` to `docs/42`; published CS construct)

| ID | Specification | Status | Evidence on size |
|------|------------------------------------------|----------|----------------------------------------|
| D1 | O1 with `ns(t, 3)` | evaluated | Hybrid t½ 1: 0.073 (O1 0.058) |
| D2 | O1 with visit factor | evaluated | Hybrid t½ 1: 0.080 |
| D3 | O1 plus fixed path effect | evaluated | power changes by at most 0.001 |
| D4 | O1 on a single path against the mixture | evaluated | no material difference in Hybrid |
| D5 | Randomized (multinomial or permuted-block) path allocation | proposed | |
| D6 | Random intercept, post-baseline shift and expectancy coefficient, `varIdent` by phase | proposed | |
| D7 | Baseline-by-drug moderator | proposed | |

Table: Approaches D1 to D7 from docs/40 to docs/42 with evidence on size

### 3.7 The covariance study: main run (published CS and configuration C)

All `mmrm`, centered biomarker, binary `Db`. Fitted to the same trials:
16 cells, 250 alternative and 500 null replicates.

| ID | Coding | Covariance | Time | Test | Size range |
|------|----------|------------|------------------|----------------|----------------|
| V1 | row | `cs` | linear, spline, visit | Satterthwaite | CS 0.036 to 0.072; configuration C, CO 0.140 to 0.158 |
| V2 | row | `csh` | linear, spline, visit | Satterthwaite | as V1 (CO, configuration C: 0.138 to 0.164) |
| V3 | row | `us` | linear, spline, visit | Satterthwaite | 0.072 to 0.104 |
| V4 | row | `sp_exp` | linear, spline, visit | Satterthwaite | 0.026 to 0.062; low power under CS |
| V5 | row | `us` | linear, visit | Kenward-Roger | 0.046 to 0.072 |
| V6 | covariate | `cs` | linear, spline, visit | Satterthwaite | CS 0.034 to 0.064; configuration C, CO 0.146 to 0.164 |
| V7 | covariate | `us` | linear, spline, visit | Satterthwaite | 0.072 to 0.096 |

Table: Covariance study main-run approaches V1 to V7 with size ranges

### 3.8 The covariance study: follow-up run (running; same trials)

| ID | Coding | Covariance | Time | Test | Question |
|------|----------|------------|------------------|----------------|----------------------------------------|
| F1 | covariate | `us` | linear, spline, visit | Kenward-Roger | is the covariate coding correctly sized? |
| F2 | covariate | `cs` | linear, visit | Kenward-Roger | does KR change a parsimonious model? |
| F3 | row | `us` | spline | Kenward-Roger | completes V5 |
| F4 | row | `cs`, `csh` | linear, visit | Kenward-Roger | can KR repair a misspecified covariance? |
| F5 | row | `cs`, `csh`, `us`, `sp_exp` | linear, visit | CR2 | does CR2 hold size under every covariance? |
| F6 | covariate | `cs`, `us` | linear, visit | CR2 | as F5, covariate coding |
| F7 | row, separate baseline biomarker effect | `us` | visit | Kenward-Roger | does the shared biomarker effect explain the row coding's power? |

Table: Covariance study follow-up approaches F1 to F7 and their questions

## 4. What is known about size

**Main conclusions:**

1. **The covariance and its small-sample correction dominate.** Changing
   the exposure coding (02), the time adjustment (`docs/40`) or the
   baseline coding (covariance study) moves the Type I error by at most
   about 0.02. Changing the working covariance or the correction moves
   it by up to 0.10.
2. **A parsimonious covariance that is wrong fails in either
   direction.** The published random intercept is anticonservative
   under AR(1) data in CO (0.14 to 0.16). The random intercept with
   `corCAR1` is conservative in the package construct (0.029 to 0.039).
   CR2 restores the latter. Kenward-Roger does not, being inert on
   parsimonious structures (10).
3. **An unstructured covariance needs Kenward-Roger.** With Satterthwaite
   alone it is anticonservative in every construct tested (0.072 to
   0.104 in the covariance study, 0.076 in 10). Kenward-Roger brings it
   to nominal.
4. **Sandwich estimators need small-sample corrections.** Naive GEE is
   anticonservative. Mancl-DeRouen and CR2 restore the size.
5. **A residual size problem remains for flexible time adjustments in
   Hybrid under carryover** (about 0.07). It persists under every
   covariance tested, including `us` with Kenward-Roger, and its cause
   is unexplained.

**Construct dependence.** Results come from three constructs: published
CS, configuration C, and the package construct. Approaches have not all
been tested in all three. The parsimonious covariances are the most
construct-dependent: the random intercept is well calibrated under CS
data (10 M1; covariance study V1 under CS) and fails under AR(1) data.

## 5. Approaches correctly sized in every cell tested

| Approach | Constructs | Size range | Power relative to the published analysis |
|----------------------------------------|----------------|------------|----------------------------------------|
| Row, `us`, Kenward-Roger, linear or visit (V5) | published CS, configuration C | 0.046 to 0.072 | similar under CS; +0.10 to +0.12 in Hybrid under configuration C |
| Random intercept plus `corCAR1`, CR2 (02 G6 to G9; 10) | package | 0.045 to 0.049 | slightly higher than model-based (recovers its conservatism) |
| GEE, Mancl-DeRouen (T5) | package | nominal | modest power cost against the mixed model (08) |
| `sp_exp`, Satterthwaite (V4) | published CS, configuration C | 0.026 to 0.062 | much lower under CS (0.31 to 0.52 against about 0.70) |

Table: Approaches correctly sized in every cell tested, with relative power

V5's upper value of 0.072 occurs with the visit-factor time adjustment
in Hybrid at t½ = 1 (published CS), the residual problem of Section 4,
point 5. With linear time V5 stays within 0.052 to 0.066.

## 6. Gaps and the evaluation in progress

**Being answered by the follow-up run (F1 to F7):** Kenward-Roger and
CR2 for the baseline-covariate coding; whether CR2 rescues the
misspecified parsimonious covariances under AR(1) data; whether
Kenward-Roger does (expected not); and the mechanism behind the row
coding's power under AR(1).

**Not yet scheduled:**

1. The cause of the flexible-time size excess in Hybrid under carryover
   (Section 4, point 5).
2. The random-effects structure implied by the published construct, with
   phase-specific variance (D6).
3. A random slope on time (O3) and the expectancy term (O2) by
   simulation, not only in the population-mean limit.
4. GEE with Mancl-DeRouen and the summary-measures statistic in the
   published CS and configuration C constructs, for a direct comparison
   with V5.
5. Exposure-weighted coding (`Dbc`) crossed with the covariance choices.
   Manuscript 02 evaluated it only under `corCAR1`.
6. Randomized path allocation (D5).
7. Power at matched size, for every approach found correctly sized.
8. Bayesian hierarchical models and randomization tests, named in the
   literature reviews of 02 and 08 and not evaluated.

**Suggested core comparison.** Once the follow-up run is in, a single
paired study in the three constructs could cover:

- the approaches correctly sized so far (V5, CR2 variants, GEE with
  Mancl-DeRouen);
- the published analysis as reference (O1);
- crossed with the exposure coding (`Db` against `Dbc`) and the baseline
  coding.

That would give one table of size and power from which a pre-specified
primary analysis for a Hybrid trial could be chosen.

## 7. Sources

- `orig`: `github.com/rchendrickson/pmsimstats`, commit `58b32a9`,
  `R/lme_analysis.R` and vignette 1 (inspected).
- Manuscripts: `analysis/report/02-carryover-sensitivity/`,
  `06-component-decomposition/`, `08-test-procedure-design-sensitivity/`,
  `10-interaction-test-calibration/` (read for their specifications and
  headline results; the numbers quoted are those reported in the
  manuscripts).
- Project documents: `docs/26`, `docs/31`, `docs/36`, `docs/37`,
  `docs/40`, `docs/41`, `docs/42`.
- Covariance study: `analysis/scripts/quick-sim/carryover-closed-form/08-covariance-study.R`;
  results in
  `analysis/data/quick-sim/carryover-closed-form/covariance-study/summary-alt250-null500.csv`
  (verified); follow-up running.
