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
  including its follow-up run;
- the 14 analyses of the stress test (`docs/46`,
  `11-strawman-stress-test.R`), added on 2026-10-07 (Section 3.9).

The approaches are organized along the dimensions on which an analysis
can vary (Section 2). The named approaches are then listed with their
source, their status, and what is known about the size of their test
(Sections 3 and 4). Section 5 summarizes which approaches are known to
hold their size, and Section 6 lists the gaps.

**Status note (2026-10-07).** The catalog was written on 2026-10-04,
while the covariance-study follow-up was running. Its status entries
have been updated with dated notes: the follow-up is reported in
`docs/45`, and the stress test of `docs/46` has since evaluated 14
analyses in 41 cells. Two statements of the original are corrected
where they occur: the unstructured covariance with Kenward-Roger is not
correctly sized everywhere, and the flexible-time size excess is
confined to model-based and Kenward-Roger inference.

(Updated 2026-10-09: results at $t_{1/2} > 0$ now use anchored
carryover; see `docs/53`. Values changed by at most 0.002 (the
`docs/40` time-adjustment sizes of O1 and D1); conclusions unchanged.
The covariance-study and stress-test entries still report the runs with
the published recursive carryover, whose anchored reruns are pending.)

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. Symbols marked *local*
are defined here. (Notation, 2026-10-09: the measured belief covariate,
formerly $\eta$, is now $\tilde\eta_{it}$, because $\eta_t$ is the design
expectancy. Manuscript 06's $\eta(\text{phase})$ is the design
expectancy, not a measured covariate; see Section 2.4.)

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $\beta_{bm:D}$ | biomarker-by-drug interaction, the estimand | canonical |
| $D_{it}$, $D_{bc,it}$ | binary drug state (`Db`); exposure-decayed drug indicator (`Dbc`) | canonical |
| $\lambda$, $t_{sd}$, $t_{1/2}$ | decay rate; time since discontinuation; carryover half-life | canonical |
| $t_{1/2}^{a}$ | half-life assumed by the analysis (`docs/45`) | canonical |
| $L_{it}$ | lagged just-off-drug indicator (E3) | canonical |
| $\tilde\eta_{it}$ | measured placebo-belief covariate, reserved; not used by any project analysis (written $\eta$ before 2026-10-09; $\eta_t$ is the design expectancy, code `De`, and manuscript 06's $\eta(\text{phase})$) | canonical |
| $c_{bm}$ | biomarker-BR coupling in the stress-test constructs | canonical |
| $B_i$ | biomarker (code `bm`; centered, `bmc`, in the stress test) | canonical |
| $\beta_B$ | biomarker main effect (`docs/41`) | canonical |
| $\rho$ | AR(1) serial correlation per week | canonical |
| $N$ | total sample size across paths | canonical |

Table: Notation index and glossary of symbols used in this document

**Identifiers.** The IDs of Section 3 are local to this catalog. The
stress test of `docs/46` uses its own IDs (S, A, C1, C2, DS, DA, ES, EA,
F1, F2, R, G1, G2, F1A), several of which collide with IDs here (C1,
C2, F1, F2, G1, G2). Following the NOTATION convention for paper-local
labels, they are cited here with a prefix: 46-F1 is the stress test's
F1, and F1 alone is the covariance-study follow-up arm of Section 3.8.

- **$\beta_{bm:D}$.** The biomarker-by-drug interaction, the estimand.
- **`Db`, `Dbc`.** Binary drug state, $D_{it}$; exposure-decayed drug
  indicator, $D_{bc,it}$ (1 on drug, $e^{-\lambda t_{sd}}$ off drug).
- **Size.** Type I error of the interaction test; nominal 0.05. An
  approach is *correctly sized* in a cell if its null rejection rate is
  within about 2 Monte Carlo standard errors of 0.05.
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
  *configuration C* is separable AR(1) with graded coupling, at
  $\rho = 0.7$ in the covariance study and the stress test (the
  covariance construct of the refined process of paper 01; both studies
  kept the published recursive carryover, so neither ran the refined
  process in full); *package*
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

   *Corrected (2026-10-07).* The unstructured covariance with
   Kenward-Roger is slightly liberal over the eight covariance-study
   null cells (mean 0.058; `docs/45`, Section 4.1) and fails the
   stress test of `docs/46`, where it is 46-A: pooled null rejection
   0.057, 0.087 at $N = 35$, and 0.195 to 0.416 in the three
   original prognostic cells. The flexible-time excess is real
   under model-based and Kenward-Roger inference but absent under CR2
   (`docs/45`, Section 4.4; `docs/42`, update after Section 1).
5. **Open:** the Kenward-Roger and CR2 versions of the baseline-covariate
   coding, CR2 with structured covariances under AR(1) data, and the
   mechanism behind the row coding's power under AR(1). All are in the
   follow-up run now in progress (Section 6). *Answered (2026-10-07)* by
   `docs/45`: the covariate coding is slightly liberal with
   Kenward-Roger (0.059) and nominal with `cs` and CR2 (0.048); CR2 with
   a parsimonious covariance is nominal in every construct (0.044 to
   0.054); the row coding's power under AR(1) comes from the shared
   biomarker effect $\beta_B$ (Section 3.8).
6. **The stress test (2026-10-07; `docs/46`, Section 3.9).** 14 analyses,
   41 cells, 28,000 trials. Ten survive the size and sensitivity
   screens: 46-S, 46-C1, 46-C2, 46-DS, 46-ES, 46-EA, 46-F1, 46-F2,
   46-G1 and 46-G2. Five of them are also nominal in every prognostic
   cell: 46-ES, 46-EA, 46-F1, 46-G1 and 46-G2. Among these, 46-G2 and
   46-F1 tie on maximin shortfall (0.232 and 0.242), and the rule
   selects 46-F1, which has fewer covariance parameters: a random
   intercept with `corCAR1` residuals, CR2 inference, and the
   expectancy term `De` with `bmc:De` (`docs/49`).

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
| Lag term crossed with biomarker ($L_{it}$ and $B_i L_{it}$) | E5 | 02, earlier S7 | retired: worse than E1 |
| $t_{sd}$ crossed with biomarker | E6 | 02, earlier S7 | retired: worse than E1 |
| `Dbc` with half-life chosen by AIC | E7 | 02 G4; stress test (46-DS, 46-DA, 46-F1A) | evaluated: nominal in one cell (02); in the stress test 46-DS passes, 46-DA and 46-F1A fail size (pooled 0.065 and 0.057) |
| `Dbc` with alternative decay forms (Weibull, linear, power) | | 02 sensitivity blocks; NOTATION | evaluated as misspecification |

Table: Drug-exposure coding options E1 to E7 with source and status

### 2.3 Time adjustment

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| None | E9 (`docs/37`) | evaluated |
| Linear `t` | orig; 02; all mixed models above | evaluated |
| Linear `t` plus expectancy `De` | orig option (`useDE`); `docs/40`; stress test (46-F1, with `bmc:De`) | evaluated (2026-10-07): 46-F1 nominal and prognostic-robust; the selected analysis |
| Natural spline in `t`, 3 df | `docs/40`, `docs/42`; covariance study | evaluated: size excess in the Hybrid carryover cell under model-based and Kenward-Roger inference, none under CR2 |
| Visit as factor | `docs/40`, `docs/42`; covariance study | evaluated: as the spline |
| Biomarker by linear time (`bmc:t`) | stress test (46-F2) | evaluated (2026-10-07): nominal, but fails the expectancy-prognostic cell |
| Phase indicators, phase by drug (attribution) | 06 (phase-augmented) | evaluated (06): no bias reduction, less power |
| Random slope on `t` | orig option (`t_random_slope`) | not evaluated |
| Path-specific intercepts (E9 stage 2) | `docs/37` (stratified E9) | closed form only |

Table: Time adjustment options with source and evaluation status

### 2.4 Expectancy

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| Not modeled | orig as published | evaluated |
| Design expectancy `De` as a fixed effect | orig option | population-mean limit (`docs/40`); evaluated with `bmc:De` in the stress test (46-F1, 46-G1, 46-G2, 46-F1A; 2026-10-07) |
| Absorbed by visit effects | `docs/40`, `docs/42` | evaluated |
| Design expectancy assigned by phase, $\eta(\text{phase})$, with phases that cross belief and exposure | 06 (balanced-placebo design) | evaluated (06) |
| Measured belief covariate $\tilde\eta_{it}$ | none; recommended by 06 as a design addition | not evaluated |

Table: Expectancy modeling options with source and evaluation status

(Corrected 2026-10-09: this table formerly listed a *measured* belief
covariate as evaluated in manuscript 06. Manuscript 06's
$\eta(\text{phase})$ is the design expectancy, assigned by phase (1
open label, about 0.5 blinded, 0 in the covert on-drug phase and 1 in
the open-placebo phase of its balanced-placebo design). The paper
calls the blinded value "a modeling assumption, not an estimate" and
recommends that trials "collect an expectancy measure", so no measured
belief covariate was used. In canonical notation the object is
$\eta_t$; $\tilde\eta_{it}$ stays reserved for a measured covariate,
which no project analysis has yet used.)

### 2.5 Baseline handling

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| Ninth response row, coded off drug | orig | evaluated |
| Excluded, no covariate | 02 (all nine) | evaluated (02) |
| Covariate (ANCOVA) | `docs/41`; covariance study | evaluated with `cs` and `us`: Satterthwaite (main run), Kenward-Roger and CR2 (follow-up, `docs/45`) |
| cLDA: row with common baseline mean, unstructured covariance | `docs/41` | equivalent to the covariate coding only with $\beta_B$ freed at baseline (`docs/45`, Section 6) |
| Change from baseline | `docs/41` | identical to ANCOVA in the construct |
| Baseline-by-drug moderator | `docs/41` | not evaluated (zero in the construct) |
| Row, with a separate biomarker effect at baseline | follow-up run (`us-KR-sepbm`); stress test (46-ES, 46-EA) | evaluated (2026-10-07): 46-ES and 46-EA nominal and prognostic-robust |

Table: Baseline handling options with source and evaluation status

### 2.6 Random effects

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| None (marginal covariance only) | MMRM; 10 (M0, M2) | evaluated |
| Random intercept | orig; 02; 10 (M1, M3) | evaluated |
| Random intercept and slope on `t` | orig option | not evaluated |
| Random intercept, post-baseline shift and expectancy coefficient, phase-specific residual variance | this discussion (implied by the published construct) | partly evaluated (2026-10-07): 46-G1 and 46-G2 have a fixed expectancy term and phase-specific residual variances, without the random shift or random expectancy coefficient |

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
| Phase-specific variance (`varIdent`) | this discussion; stress test | evaluated (2026-10-07): 46-C1, 46-C2, 46-G1, 46-G2 |
| GEE working correlation: independence, exchangeable, AR(1) | 08; 10 (AR(1)) | evaluated |

Table: Working residual covariance options with source and evaluation status

### 2.8 Standard error and degrees of freedom

| Option | Source | Status |
|---------------------------------|------------------------|----------------|
| Model-based SE, Satterthwaite df (`lmerTest`) | orig | evaluated |
| Model-based SE, containment df (`nlme`) | 02; 06; 08 | evaluated |
| Model-based SE, Satterthwaite df (`mmrm`) | covariance study | evaluated |
| Kenward-Roger | covariance study (`us`); 10 | evaluated for `us`, `cs`, `csh` and the baseline covariate (`docs/45`): repairs `us` only |
| CR2 sandwich, Satterthwaite df (`clubSandwich`) | 02 (G6 to G9); 10; stress test | evaluated; the inference of 46-S, 46-C1, 46-DS, 46-ES, 46-F1, 46-F2, 46-G1 and 46-F1A |
| CR2 sandwich, Satterthwaite df (`mmrm`) | follow-up run | evaluated (`docs/45`): nominal with parsimonious covariances, liberal with `us` (0.086) |
| CR0, CR1, CR3 sandwiches | 10 (discussed) | not evaluated |
| Naive GEE sandwich | 08; 10 | evaluated |
| Mancl-DeRouen GEE sandwich | 08; 10 | evaluated |
| OLS $t$, $N - 2$ df (paths pooled, common intercept) | E9 | evaluated |
| OLS $t$, $N - P - 1$ df (path-specific intercepts) | stratified E9 (compendium, eq:invf-strat; `04-design-search.R`) | closed form only (design search) |

Table: Standard error and degrees-of-freedom options with source and status

(Checked 2026-10-09: E9 regresses $\Delta_i$ on $B_i$ by OLS over all
$N$ participants with one common intercept, pooling the paths, so
$N - 2$ degrees of freedom is correct (`docs/37`, Sections 3 and 4.4;
compendium, Section "The paired-difference statistics E9 and E9b, and
their moments"; inspected in `01-e9-closed-form.R` and
`06-single-path-lme.R`). The path differences in the mean contrast
enter its residual variance, not its degrees of freedom. A regression
on a single path of $n$ participants has $n - 2$; the stratified
variant with path-specific intercepts has $N - P - 1$.)

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
| O1 | `lmer`, fixed part `Sx ~ bm + Db + t + bm*Db`, random intercept by `ptID`, baseline row, Satterthwaite (published default) | evaluated | published CS: 0.049 to 0.056 (Hybrid), 0.04 to 0.05 (CO). Published results: 0.01 to 0.09. AR(1) data: 0.14 to 0.23 in CO (`docs/36`) |
| O2 | O1 plus `De` (`useDE = TRUE`) | proposed | dropped by Hendrickson et al. for collinearity; correlation with `Db` 0.59 (`docs/40`). Not simulated as specified; a CR2 analogue with `bmc:De` is 46-F1 |
| O3 | O1 with a random intercept and random slope on `t` by `ptID` (`t_random_slope = TRUE`) | proposed | |
| O4 | OL branch: fixed part `Sx ~ bm + t + bm*t`, random intercept by `ptID` | evaluated | confounded test; not a drug-by-biomarker interaction |
| O5 | Later vendored implementation: `nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc)`, random intercept, `corCAR1` | evaluated as 02 G3 | see 02 G3 |

Table: Published analysis O1 and its options O2 to O5, with evidence on size

### 3.2 Manuscript 02 (carryover specifications; package construct)

All use `nlme::lme` with a random intercept and `corCAR1`, no baseline
row, linear `t`.

| ID | Specification | Status | Evidence on size |
|------|------------------------------------------|----------|----------------------------------------|
| G1 | `Db`, target `bm:Db` | evaluated | mildly conservative (mean null 0.036 to 0.039) |
| G2 | G1 plus lag $L_{it}$ | evaluated | as G1; no power gain over G1 |
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
| D1 | O1 with `ns(t, 3)` | evaluated | Hybrid, $t_{1/2} = 1$: 0.075 (O1 0.056); recurs in `docs/45` (0.070), absent under CR2 |
| D2 | O1 with visit factor | evaluated | Hybrid, $t_{1/2} = 1$: 0.080; recurs in `docs/45` (0.072), absent under CR2 |
| D3 | O1 plus fixed path effect | evaluated | power changes by at most 0.001 |
| D4 | O1 on a single path against the mixture | evaluated | no material difference in Hybrid |
| D5 | Randomized (multinomial or permuted-block) path allocation | proposed | |
| D6 | Random intercept, post-baseline shift and expectancy coefficient, `varIdent` by phase | proposed; partly evaluated as 46-G1, 46-G2 (fixed, not random, expectancy coefficient) | 46-G1 pooled 0.049, 46-G2 0.047 |
| D7 | Baseline-by-drug moderator | proposed | |

Table: Approaches D1 to D7 from docs/40 to docs/42 with evidence on size

### 3.7 The covariance study: main run (published CS and configuration C at $\rho = 0.7$)

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

### 3.8 The covariance study: follow-up run (completed; same trials)

The run was in progress when this catalog was written. The last column
was added on 2026-10-07 from `docs/45` (Section 4.1 mean sizes over the
eight null cells, linear time; Section 6 for F7).

| ID | Coding | Covariance | Time | Test | Question | Answer (`docs/45`) |
|------|----------|------------|------------------|----------------|------------------------|------------------------|
| F1 | covariate | `us` | linear, spline, visit | Kenward-Roger | is the covariate coding correctly sized? | slightly liberal, 0.059 |
| F2 | covariate | `cs` | linear, visit | Kenward-Roger | does KR change a parsimonious model? | no; liberal, 0.077 |
| F3 | row | `us` | spline | Kenward-Roger | completes V5 | same class as linear time, within 0.005 |
| F4 | row | `cs`, `csh` | linear, visit | Kenward-Roger | can KR repair a misspecified covariance? | no; 0.079 and 0.073 |
| F5 | row | `cs`, `csh`, `us`, `sp_exp` | linear, visit | CR2 | does CR2 hold size under every covariance? | nominal for `cs`, `csh`, `sp_exp` (0.048 to 0.054); liberal for `us` (0.086) |
| F6 | covariate | `cs`, `us` | linear, visit | CR2 | as F5, covariate coding | `cs` nominal (0.048); `us` liberal (0.080) |
| F7 | row, separate baseline biomarker effect | `us` | visit | Kenward-Roger | does the shared biomarker effect explain the row coding's power? | yes: freed, power matches the covariate coding (0.780 against 0.792) |

Table: Covariance study follow-up approaches F1 to F7, their questions and answers

### 3.9 The stress test (`docs/46`; added 2026-10-07)

All 14 analyses were fitted to the same simulated Hybrid trials: 41
cells (39 in the original working protocol, plus two prognostic-baseline cells, 40 and
41), 15 null cells of 1,000 replicates and 26 alternative cells of 500,
28,000 trials in all. The base cell is configuration C at
$\rho = 0.7$, $N = 70$, $t_{1/2} = 0.5$ week, with graded coupling, the
positive-definite stop and the published recursive carryover (inspected
in `11-strawman-stress-test.R`). Every analysis has
linear time and a centered biomarker (`bmc`), and uses $D_{bc,it}$ at
$t_{1/2}^{a} = 0.5$ week, except 46-R (binary $D_{it}$) and the AIC
arms 46-DS, 46-DA and 46-F1A (half-life chosen over
$\{0, 0.25, 0.5, 1, 2\}$). IDs carry the 46- prefix of the Identifiers
note above.

| ID | Specification | Pooled size | Worst non-prognostic cell | Prognostic cells (BL, PB, TV constant, TV growing) | Outcome |
|------------|----------------------------|--------|------------|------------------------|------------------|
| 46-S | random intercept, `corCAR1`, CR2, baseline row | 0.048 | 0.059 | 0.050, 0.100, 0.099, 0.049 | survives; flagged |
| 46-A | `us`, Kenward-Roger, baseline row (`docs/45` primary) | 0.057 | 0.087 ($N = 35$) | 0.051, 0.222, 0.416, 0.195 | fails size |
| 46-C1 | 46-S with `varIdent` by phase, CR2 | 0.053 | 0.069 | 0.054, 0.205, 0.382, 0.173 | survives; flagged |
| 46-C2 | 46-C1, model-based | 0.052 | 0.066 | 0.050, 0.193, 0.393, 0.172 | survives; flagged |
| 46-DS | 46-S, $t_{1/2}^{a}$ chosen by AIC | 0.053 | 0.068 | 0.055, 0.112, 0.115, 0.058 | survives; flagged |
| 46-DA | 46-A, $t_{1/2}^{a}$ chosen by AIC | 0.065 | 0.103 ($N = 35$) | 0.064, 0.256, 0.509, 0.245 | fails size |
| 46-ES | 46-S with $\beta_B$ freed at baseline | 0.048 | 0.059 | 0.050, 0.049, 0.049, 0.063 | survives; robust |
| 46-EA | 46-A with $\beta_B$ freed at baseline | 0.052 | 0.070 | 0.066, 0.050, 0.049, 0.057 | survives; robust |
| 46-F1 | 46-S plus `De + bmc:De` | 0.051 | 0.062 | 0.040, 0.040, 0.047, 0.046 | survives; robust; **selected** |
| 46-F2 | 46-S plus `bmc:t` | 0.049 | 0.059 | 0.048, 0.098, 0.090, 0.048 | survives; flagged |
| 46-R | published reference: compound symmetry (`mmrm` `cs`, equivalent to a random intercept), model-based SE with Satterthwaite df, binary $D_{it}$ | 0.069 | 0.081 | 0.075, 0.159, 0.107, 0.053 | fails size |
| 46-G1 | 46-F1 mean model, phase-specific residual variances, CR2 | 0.049 | 0.060 | 0.042, 0.042, 0.055, 0.044 | survives; robust |
| 46-G2 | 46-G1, model-based | 0.047 | 0.066 | 0.046, 0.045, 0.062, 0.043 | survives; robust |
| 46-F1A | 46-F1, $t_{1/2}^{a}$ chosen by AIC | 0.057 | 0.073 ($t_{1/2} = 1$) | 0.045, 0.047, 0.059, 0.052 | fails size |

Table: Stress-test analyses with null rejection rates and screen outcomes

Pooled size is over the 11 non-prognostic null cells (11,000 trials);
the nominal band is 0.0458 to 0.0542 (2 MCSE) with no cell above 0.075.
*Flagged* means the analysis fails at least one prognostic cell;
*robust* means it is nominal in all four. All 14 pass the monotonicity
screen; the symmetry screen was dropped in Amendment 1 of `docs/46`,
because power is symmetric in the sign of $c_{bm}$ by construction.
Maximin shortfall of size-adjusted power over the non-prognostic
alternative cells: 46-C1 0.034, 46-C2 0.040, 46-S 0.144, 46-F2 0.160,
46-DS 0.182, 46-G2 0.232, 46-F1 0.242, 46-G1 0.252, 46-ES 0.278,
46-EA 0.318. 46-C1 has the smallest maximin but fails the prognostic
cells. Among the robust survivors 46-G2 and 46-F1 tie within 0.02;
neither reaches power 0.80 anywhere on the power-curve grid, so the
tie goes to the analysis with fewer covariance parameters, 46-F1.
(Clarified 2026-10-08: "maximin shortfall" is the term used by the
decision outputs; selecting the smallest largest shortfall is a
minimax-regret rule, and the best analysis in each cell is taken over
all ten survivors, the flagged 46-C1 and 46-C2 included, so the robust
analyses' shortfalls are measured against the phase-variance models.
The 46-R row was also reworded on 2026-10-08: the driver fits it as the
`mmrm` compound-symmetry structure with Satterthwaite degrees of
freedom, not as an `lme` random intercept with containment degrees of
freedom. Neither change alters any rate or the selection. See
`analysis/report/compendium-mathematics.tex`, Sections "Size, power
screens and the selection criteria of the stress test" and "The
stress-test construct and its fourteen analyses".)
Failed fits were 0 to 12 per analysis out of 28,000 (at most 0.6% in a
cell; verified from the `fails` column of `summary-full*.csv`). The working protocol counted failed fits as non-rejections; the
summaries exclude them, which changes no rate by more than 0.0025.
(Verified against
`analysis/data/quick-sim/carryover-closed-form/strawman-stress/decision-*-ds-nosym-all14.csv`;
decision rule rerun with all 14 analyses, 2026-10-07.)

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
   to nominal. *Corrected (2026-10-07):* Kenward-Roger brings it close
   to nominal but not to it: slightly liberal over the covariance
   study's null cells (mean 0.058, `docs/45`), and failing the stress
   test's size criterion (46-A: pooled 0.057, 0.087 at $N = 35$).
4. **Sandwich estimators need small-sample corrections.** Naive GEE is
   anticonservative. Mancl-DeRouen and CR2 restore the size.
5. **A residual size problem remains for flexible time adjustments in
   Hybrid under carryover** (about 0.07). It persists under every
   covariance tested, including `us` with Kenward-Roger, and its cause
   is unexplained. *Corrected (2026-10-07):* the follow-up run tested
   CR2, under which the excess is absent (0.058 to 0.060 in that cell);
   under model-based and Kenward-Roger inference it is real, recurring
   in a second run at 0.068 to 0.072 (`docs/45`, Section 4.4, whose
   description of this evidence as weak is itself corrected in
   `docs/42`). Its cause remains unexplained.
6. **The stress test separates the prognostic-robust analyses (added
   2026-10-07).** Every analysis that neither frees $\beta_B$ at
   baseline nor adds the expectancy term fails at least one cell in
   which the biomarker is prognostic for $TV$ or $PB$ (null rejection up
   to 0.509). Freeing $\beta_B$ at baseline (46-ES, 46-EA) or adding
   `De + bmc:De` (46-F1, 46-G1, 46-G2) gives nominal size in every
   prognostic cell; `bmc:t` (46-F2) repairs only the cell with a growing
   natural-course association (Section 3.9).

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
in Hybrid at $t_{1/2} = 1$ (published CS), the residual problem of Section 4,
point 5. With linear time V5 stays within 0.052 to 0.066.

*Corrected (2026-10-07).* Two rows of this table do not meet the
catalog's own criterion. V5 is slightly liberal (mean 0.058, `docs/45`)
and its stress-test counterpart 46-A, the same model with $D_{bc,it}$
at 0.5 week, fails the size criterion (pooled 0.057; 0.087 at
$N = 35$) and all three original prognostic cells (0.195 to
0.416). It is therefore not recommended, and `docs/45`'s recommendation
of it as the primary analysis is superseded by `docs/46` and `docs/49`.
V4 (`sp_exp`, Satterthwaite) is conservative (mean 0.042, minimum
0.028; `docs/45`). The analyses now known to be correctly sized and
robust in every stress-test cell, prognostic cells included, are 46-ES,
46-EA, 46-F1, 46-G1 and 46-G2 (Section 3.9), all in the configuration C
and compound-symmetry constructs of `docs/46`; 46-F1 is the selected
analysis.

## 6. Gaps and the evaluation in progress

**Being answered by the follow-up run (F1 to F7):** Kenward-Roger and
CR2 for the baseline-covariate coding; whether CR2 rescues the
misspecified parsimonious covariances under AR(1) data; whether
Kenward-Roger does (expected not); and the mechanism behind the row
coding's power under AR(1). *Answered (2026-10-07)* in `docs/45`; see
the last column of Section 3.8.

**Not yet scheduled** (status of each item added 2026-10-07):

1. The cause of the flexible-time size excess in Hybrid under carryover
   (Section 4, point 5). *Open;* CR2 is known to remove it.
2. The random-effects structure implied by the published construct, with
   phase-specific variance (D6). *Partly done:* 46-G1 and 46-G2 (fixed
   expectancy term, phase-specific residual variances).
3. A random slope on time (O3) and the expectancy term (O2) by
   simulation, not only in the population-mean limit. *Partly done:* the
   expectancy term with `bmc:De` under CR2 is 46-F1; O3 is open.
4. GEE with Mancl-DeRouen and the summary-measures statistic in the
   published CS and configuration C constructs, for a direct comparison
   with V5. *Open.*
5. Exposure-weighted coding (`Dbc`) crossed with the covariance choices.
   Manuscript 02 evaluated it only under `corCAR1`. *Done:* `docs/45`
   exposure pass (`us-KR`, `us-CR2`, `RI+CAR1-CR2`; no change in size)
   and every stress-test analysis except 46-R.
6. Randomized path allocation (D5). *Open.*
7. Power at matched size, for every approach found correctly sized.
   *Done* for the analyses of `docs/45` and `docs/46`.
8. Bayesian hierarchical models and randomization tests, named in the
   literature reviews of 02 and 08 and not evaluated. *Open.*

**Suggested core comparison.** Once the follow-up run is in, a single
paired study in the three constructs could cover:

- the approaches correctly sized so far (V5, CR2 variants, GEE with
  Mancl-DeRouen);
- the published analysis as reference (O1);
- crossed with the exposure coding (`Db` against `Dbc`) and the baseline
  coding.

That would give one table of size and power from which a primary
analysis for a Hybrid trial could be chosen. *Status
(2026-10-07):* the exploratory stress test of `docs/46` is that study
for the Hybrid design, in configuration C and compound symmetry but not
the package construct, and without GEE or the summary-measures
statistic. It selected 46-F1 (Section 3.9; `docs/49`).

## 7. Sources

- `orig`: `github.com/rchendrickson/pmsimstats`, commit `58b32a9`,
  `R/lme_analysis.R` and vignette 1 (inspected).
- Manuscripts: `analysis/report/02-carryover-sensitivity/`,
  `06-component-decomposition/`, `08-test-procedure-design-sensitivity/`,
  `10-interaction-test-calibration/` (read for their specifications and
  headline results; the numbers quoted are those reported in the
  manuscripts).
- Project documents: `docs/26`, `docs/31`, `docs/36`, `docs/37`,
  `docs/40`, `docs/41`, `docs/42`; added 2026-10-07: `docs/45`,
  `docs/46`, `docs/48`, `docs/49`.
- Covariance study: `analysis/scripts/quick-sim/carryover-closed-form/08-covariance-study.R`;
  results in
  `analysis/data/quick-sim/carryover-closed-form/covariance-study/summary-alt250-null500.csv`
  (verified); follow-up completed and reported in `docs/45`.
- Stress test: `analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R`
  and `12-strawman-stress-decision.R`; results in
  `analysis/data/quick-sim/carryover-closed-form/strawman-stress/`
  (`decision-*-ds-nosym-all14.csv`, `summary-full*.csv`; verified).
