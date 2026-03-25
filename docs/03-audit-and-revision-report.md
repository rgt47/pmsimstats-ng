# pmsimstats: Status Report on N-of-1 Trial Simulation Framework

**Ronald G. Thomas, Ph.D.**
Department of Family Medicine and Public Health, UC San Diego

**Collaborator:** Rebecca C. Hendrickson, M.D., Ph.D.
VA Puget Sound Health Care System

---

## 1. Background and Objective

The pmsimstats framework provides Monte Carlo simulation-based
power analysis for aggregated N-of-1 clinical trial designs. It
compares four trial designs -- Open-Label (OL), Crossover (CO),
Hybrid N-of-1, and Open-Label with Blinded Drug Cessation
(OL+BDC) -- on statistical power to detect a biomarker-treatment
interaction. The motivating clinical application is the use of
resting blood pressure as a predictive biomarker for prazosin
response in PTSD, based on the Raskind et al. (2013) RCT analyzed
in Murray et al. (in preparation).

The original simulation code was published in Hendrickson et al.
(2020), *Frontiers in Digital Health* 2:13. The current work
represents a comprehensive audit, revision, and consolidation of
the codebase into a single reproducible research repository
(`pmsimstats-ng`) housed as a zzcollab Docker workspace.

---

## 2. Audit and DGP Revisions

A technical audit of the original publication code identified
seven substantive corrections to the data-generating process
(DGP), each justified by either statistical theory or clinical
plausibility.

### 2.1 AR(1) autocorrelation (was compound symmetry)

The original code used compound symmetry (constant correlation
$\rho$ between all timepoint pairs). The revised code implements
AR(1) decay: $\text{Cor}(X_{t_i}, X_{t_j}) = \rho^{|t_i - t_j|}$.
This change is both more clinically realistic (nearby measurements
are more correlated than distant ones) and numerically beneficial:
it expands the positive-definite-feasible region of the biomarker
correlation parameter $c_{bm}$ from a maximum of 0.25 to 0.49+.

### 2.2 Exponential BM-BR correlation decay

During off-drug periods, the correlation between the biomarker
and the biological response should decay as the drug effect wanes.
The original code used a step function (correlation drops to zero
instantly at discontinuation). The revised code implements
exponential decay tied to the pharmacokinetic half-life:
$\text{Cor}(bm, br_t) = c_{bm} \cdot e^{-\lambda_{cor} \cdot
t_{sd}}$, where $\lambda_{cor} = \ln(2) / t_{1/2}$. This
derivation is presented in the correlation decay derivation white
paper.

### 2.3 nlme::lme with corCAR1 analysis model

The original analysis model used `lme4::lmer` with a random
intercept only. Under the AR(1) DGP, this misspecified model
absorbs residual autocorrelation into the error term, inflating
Type I error to 13-17% for CO and OL designs. The revised code
uses `nlme::lme` with `corCAR1(form = ~t|ptID)`, which models
the continuous-time AR(1) residual structure directly. This
corrects Type I error to the nominal 5% across all designs. The
analysis model white paper provides the full derivation and
validation.

### 2.4 Carryover scale factor removal

The original code applied a scale factor of 2 to the carryover
half-life parameter, effectively halving the half-life without
documentation. The parameter `carryover_t1half` should mean what
it says: the time for the carryover effect to decay to half its
initial value. The scale factor has been removed.

### 2.5 Timepoint-1 correlation fix

The original code guarded the BM-BR correlation assignment with
`if (p > 1)`, skipping the first timepoint. This was a coding
artifact to avoid an index error in the proportional scaling
formula. Participants are on drug for 2.5-4 weeks at timepoint 1
and should have a biomarker-response correlation. The guard has
been removed.

### 2.6 Hendrickson-aligned Gompertz function

The Gompertz growth function `f(t) = M \cdot e^{-d \cdot
e^{-r \cdot t}}` does not pass through the origin ($f(0) \neq 0$).
The original publication code applied an offset-rescaling to
guarantee $f(0) = 0$ and preserve the asymptote at $M$. The 2025
tidyverse rewrite inadvertently used a different formula that does
not pass through the origin. This has been corrected: both
pipelines now use the Hendrickson-aligned formula.

### 2.7 Updated default parameters

Autocorrelation base rate reduced from 0.8 to 0.7 (consistent
with AR(1) decay). Biomarker correlation $c_{bm}$ constrained to
$\{0, 0.25, 0.50\}$ (maximum PD-feasible value under AR(1) is
0.49 for OL). Carryover half-life range updated from $\{0, 0.1,
0.2\}$ weeks (too short to produce genuine carryover effects) to
$\{0, 0.5, 1.0\}$ weeks.

---

## 3. Resolving the Positive Definiteness Problem

### 3.1 The problem

The DGP constructs a $26 \times 26$ covariance matrix $\Sigma$
from user-specified correlation parameters and passes it to an
MVN draw. When the resulting matrix is not positive definite (PD),
the code silently projects it to the nearest PD matrix via
`corpcor::make.positive.definite()`. Instrumentation of the
original publication code revealed that this is not a rare edge
case: 24.7% of covariance matrices required silent correction at
the publication commit, rising to 63.0% after subsequent edits.
All 162 matrices in the standard parameter sweep were
ill-conditioned ($\kappa > 100$).

### 3.2 Root causes

The failures are concentrated in the OL+BDC design (83-92%
failure rate) due to an abrupt expectancy-driven change in the
ER (expectancy response) variance mid-trial. The OL+BDC design
transitions from full expectancy ($e = 1$) during the open-label
phase to half expectancy ($e = 0.5$) during the blinded drug
cessation phase. Because the ER standard deviation is scaled by
$e$, this creates a factor-of-two variance discontinuity within
the same correlation block, violating the implicit assumption
that the specified correlations are consistent with the marginal
variances. The compound symmetry autocorrelation structure
compounds the problem: it requires all same-factor timepoint
pairs to have the same correlation regardless of their temporal
separation, which becomes untenable when variances change
mid-trial.

The OL design never fails because all timepoints have the same
expectancy ($e = 1$), so there is no variance heterogeneity.

### 3.3 Resolution

The revised code resolves the PD problem through two mechanisms:

1. **AR(1) autocorrelation** replaces compound symmetry. Because
   AR(1) correlations decay with temporal distance ($\rho^{|t_i -
   t_j|}$), the correlation matrix has smaller off-diagonal entries
   for distant timepoints, making it substantially easier to
   satisfy the PD constraint. The maximum feasible $c_{bm}$
   increases from 0.25 to 0.49+.

2. **Pre-validation with constrained parameters.** The
   `validateParameterGrid()` function tests all parameter
   combinations for PD before simulation and reports failures
   upfront. The default parameter grid ($c_{bm} \leq 0.45$,
   $\rho = 0.7$) produces zero PD failures across all 162
   combinations for all four trial designs.

The net result: the revised code runs the complete parameter
sweep with zero silent corrections, eliminating the uncontrolled
perturbation that confounded cross-design power comparisons in
the original publication.

---

## 4. Resolving the Power Dropoff Problem

### 4.1 The problem

When carryover is introduced into the DGP (non-zero $t_{1/2}$),
the statistical power to detect the biomarker-treatment
interaction should decline because the carryover effect blurs the
distinction between on-drug and off-drug periods. However, the
magnitude and pattern of the power decline depends on the entire
variable chain from trial design specification through the
analysis model interaction term.

Five derived variables form this chain:

1. `ondrug` -- binary path assignment (design-determined)
2. `tod` -- cumulative time on drug (Gompertz input)
3. `tsd` -- time since discontinuation (carryover input)
4. `Dbc` -- continuous drug indicator in the analysis model
5. `bm:Dbc` or `bm:t` -- the interaction coefficient tested

Carryover enters at step 3 and propagates through steps 4 and 5
to reduce the variance of the interaction term, inflating its
standard error and reducing power.

### 4.2 The mechanism

For designs with both on-drug and off-drug periods (CO, Hybrid,
OL+BDC), the analysis model uses `Dbc` as the treatment
indicator. On drug, $Dbc = 1$. Off drug, $Dbc = (1/2)^{s \cdot
tsd / t_{1/2}}$, where $s$ is the scale factor and $t_{1/2}$
is the carryover half-life.

Without carryover ($t_{1/2} = 0$), $Dbc$ is effectively binary
(1 vs 0), maximizing the contrast between treatment conditions.
With carryover, $Dbc$ during off-drug periods is strictly positive
(decaying from 1 toward 0), compressing the treatment contrast.
The interaction coefficient `bm:Dbc` depends on the variance of
$Dbc$, which decreases as carryover increases, inflating the
standard error.

For the open-label design (no off-drug periods), the analysis
model uses `bm:t` instead of `bm:Dbc`, so carryover has no
effect on power -- there are no off-drug periods for carryover
to operate on.

### 4.3 Resolution

The revised code addresses the power dropoff problem at two
levels:

1. **Correct DGP-analysis alignment.** The AR(1) autocorrelation
   in the DGP is matched by `corCAR1` in the analysis model.
   This prevents residual autocorrelation from being absorbed
   into the error term, which would artificially inflate or
   deflate power estimates.

2. **Clinically realistic carryover parameters.** The original
   half-life values (0.1-0.2 weeks) produce negligible carryover
   effects that decay to essentially zero within one inter-visit
   interval. The revised values (0.5-1.0 weeks) produce genuine
   carryover that persists across multiple off-drug timepoints,
   enabling meaningful power comparisons between designs that
   differ in their susceptibility to carryover.

3. **Conditional analysis formula.** The analysis model adapts
   to the design: `bm:Dbc` when there is within-subject treatment
   variation (CO, Hybrid, OL+BDC), and `bm:t` when there is not
   (OL). This prevents the open-label power estimate from being
   confounded by a degenerate treatment indicator.

The 1,000-replicate validation confirms that Type I error rates
are within 2.5-5.8% across all designs and carryover conditions,
and power estimates decrease monotonically with increasing
carryover half-life as expected.

---

## 5. Consolidated Repository

The pmsimstats-ng repository consolidates the original
publication code (`orig`, data.table) and the 2025 tidyverse
rewrite (`2025`) into a single zzcollab Docker workspace with
shared data, documentation, and test infrastructure.

### 5.1 Architecture

The two simulation pipelines coexist independently:

- **orig** (`R/`, 1973 lines, 10 files): data.table-based
  package functions loaded via `devtools::load_all(".")`. Uses
  `buildSigma()`, `generateData()`, `lme_analysis()`,
  `generateSimulatedResults()`. This is the gold standard.

- **2025** (`analysis/2025/01-pm-functions.R`, 1023 lines):
  tidyverse-based monolithic file loaded via `source()`. Uses
  `build_sigma_matrix()`, `generate_data()`, `lme_analysis()`,
  `generate_simulated_results()`. Validated bit-identical to
  orig for the same inputs.

- **simple** (`simple/simulation.R`, 577 lines): self-contained
  pedagogical simulation using direct mean moderation (not MVN
  correlation) and OLS on phase means (not LME). Serves a
  different purpose and is architecturally independent.

### 5.2 Validation

- 67 unit tests covering `cumulative()`, `modgompertz()`,
  `buildtrialdesign()` (OL, CO, Hybrid, tpb), `buildSigma()`
  (PD, dimensions, labels, AR(1), lambda_cor, Cholesky),
  `generateData()` (output shape, caching, outcome identity),
  and `lme_analysis()` (OL, CO, p-value bounds).
- 3-level alignment test: Gompertz, sigma matrix, and generated
  data are bit-identical between orig and 2025 for the same seed.
- A bug in `lme_analysis()` (uninitialized `op` list when called
  without arguments) was discovered and fixed through the test
  suite.

### 5.3 Reproducibility infrastructure

- Docker image via zzcollab `analysis` profile
  (`rocker/tidyverse:4.5` + renv)
- `renv.lock` capturing all dependencies (nlme, data.table,
  furrr, corpcor, MASS, ggplot2)
- `Makefile` with targets: `make r` (container shell),
  `make rstudio`, `make check-renv`, `make docker-build`,
  `make test`
- GitHub repository: `rgt47/pmsimstats-ng`

---

## 6. Collaborator Deliverables

### 6.1 Goal 1: Did prazosin improve PTSD symptoms?

Implemented as `analyze_trial_extended()` in
`R/carryover_analysis.R`. Fits the LME model to actual trial data
and extracts both the main drug effect (`Dbc` coefficient) and
the biomarker-treatment interaction (`bm:Dbc`), with 95%
confidence intervals, variance components, predicted response
curves per participant, and a clinical decision threshold
(the blood pressure value at which predicted benefit exceeds a
specified CAPS improvement). Status: complete.

### 6.2 Goal 2a: Characterize carryover in actual data

Implemented as `characterize_carryover()` in
`R/carryover_analysis.R`. Sweeps a range of candidate
carryover half-lives (default 0.25 to 8 weeks), fits the LME
model at each half-life with exponentially decaying `Dbc`,
compares model fit via AIC/BIC and log-likelihood, and
identifies the best-fitting half-life. Validated on simulated
data: correctly recovers the true $t_{1/2} = 1.0$ from data
generated with that half-life. Status: complete.

### 6.3 Goal 2b: Identify best sequence of response measures

Conceptual design only. Would require per-timepoint residual
analysis to identify where the carryover signal is strongest
in the observed data. The existing model infrastructure
(`full_model_out = TRUE`) provides all necessary outputs
(fitted values, residuals, random effects per participant).
Status: not yet implemented.

### 6.4 Goal 3a: Simulation power calculations

Complete. The revised simulation produces 1,000-replicate
power estimates with proper Type I error control (2.5-5.8%
across all designs), zero PD failures across all 162 parameter
combinations, sigma matrix caching for performance (8 minutes
for 200 reps, 47 minutes for 1,000 reps), and parallel
execution via furrr. Status: complete.

### 6.5 Goal 3b: Actual analysis output structure

Implemented in `analyze_trial_extended()`. Returns the full
model object, predicted treatment effect as a continuous
function of biomarker value, per-participant predicted response
curves, and clinical decision thresholds. The function can be
applied directly to the bundled `CTdata.rda` dataset from the
Raskind et al. trial. Status: complete.

---

## 7. White Paper Documentation

The consolidated repository includes seven technical white
papers documenting the audit findings, theoretical derivations,
and validation results.

### 7.1 Revised Power Analysis (04-revised-power-analysis.pdf)

A 13-page comprehensive summary of all changes from v0.1.0
to v0.2.0, written for the original authors. Documents each
DGP revision with before/after code, parameter tables, and
heatmap comparisons. Includes the validation results showing
Type I error correction and the complete parameter grid with
power estimates at 200 and 1,000 replicates. Serves as the
primary reference for understanding what changed and why.

### 7.2 AR(1) Residual Correlation (06-ar1-residual-correlation.pdf)

Derives why the original `lme4::lmer` random-intercept model
produces inflated Type I error (13-17%) under the AR(1) DGP
and validates that `nlme::lme` with `corCAR1` corrects this
to the nominal 5%. Presents a simulation study comparing the
two analysis models across 1,000 replicates under the null
($c_{bm} = 0$) for all four trial designs. Shows that the
random-intercept model absorbs the AR(1) residual
autocorrelation into a biased variance estimate, shrinking
standard errors and inflating test statistics.

### 7.3 Biomarker Correlation Decay (08-biomarker-correlation-decay.pdf)

Derives the exponential decay formula for the biomarker-response
correlation during off-drug periods from pharmacokinetic first
principles. When each participant's drug effect decays by the
same fraction per unit time ($e^{-\lambda t}$, where $\lambda =
\ln 2 / t_{1/2}$), the between-participant spread of drug
effects decays proportionally, preserving the coefficient of
variation. This implies the correlation between biomarker and
drug-mediated response decays at the same rate:
$\text{Cor}(bm, br_t) = c_{bm} \cdot e^{-\lambda t_{sd}}$.
The derivation shows this is the unique decay function that
maintains the proportional relationship between individual drug
effects and the biomarker.

### 7.4 Positive Definiteness Failures (07-positive-definiteness-failures.pdf)

Presents the instrumented measurement of positive definiteness
failures across four code versions and 162 parameter
combinations. Documents that 24.7% of covariance matrices fail
PD at the publication commit, rising to 63.0% after subsequent
edits, with all matrices ill-conditioned ($\kappa > 100$).
Traces the root cause to the OL+BDC design's abrupt
expectancy-driven variance change (83-92% failure rate) and the
compound symmetry autocorrelation structure. Evaluates the
perturbation introduced by `make.positive.definite()` and its
consequences for cross-design power comparisons. Recommends
AR(1) autocorrelation, pre-validation, and parameter
constraints as solutions.

### 7.5 Carryover Correlation Artifact (09-carryover-correlation-artifact.pdf)

Analyzes the BM-BR correlation mechanism in the Hendrickson
DGP. The biomarker-treatment interaction emerges from
differential correlation in the MVN draw: when on drug,
$\text{Cor}(bm, br) = c_{bm}$; when off drug, the correlation
decays or drops to zero. This means participants with high
biomarker values tend to have large drug effects (high BR
means) and participants with low biomarker values tend to have
small drug effects, creating the interaction detectable by the
analysis model. The paper traces how this correlation
propagates through the sigma matrix to produce the observed
power pattern across designs and parameter values.

### 7.6 Response Parameter Sensitivity (18-response-parameter-sensitivity.pdf)

Compares Figure 5 (response parameter sensitivity analysis)
between the original publication code and the revised code.
Figure 5 sweeps the Gompertz response parameters (max, disp,
rate) and standard deviations across a grid to assess which
parameters most strongly influence power. Documents the
differences introduced by the DGP revisions and validates
that the qualitative conclusions about parameter sensitivity
are robust to the corrections.

### 7.7 Revision Summary (05-revision-summary.pdf)

Traces the impact of each individual commit (from the
publication code through the Ron Thomas edits to the current
HEAD) on the power estimates. Reproduces the Figure 4 heatmap
at each commit to show how each change affected the power
surface. Documents which changes improved power estimates
(AR(1) expanding the feasible region), which degraded them
(proportional BM-BR scaling reducing effective correlation),
and which were neutral (code cleanup without DGP changes).

---

## 8. Presentation Outline (7 Slides)

### Slide 1: Problem and Motivation

- N-of-1 trial designs for biomarker-treatment interactions
- Four designs compared: OL, CO, Hybrid, OL+BDC
- Clinical context: prazosin for PTSD, blood pressure as
  predictive biomarker (Raskind et al. 2013, Hendrickson et al.
  2020)
- Question: which trial design has the best power to detect
  who responds?

### Slide 2: Audit Findings -- Seven DGP Corrections

- AR(1) autocorrelation (was compound symmetry)
- Exponential BM-BR correlation decay (was step function)
- nlme + corCAR1 analysis model (was lmer, Type I error
  13-17% corrected to 5%)
- Scale factor removal, timepoint-1 fix, Gompertz alignment
- Updated defaults: rho=0.7, c.bm <= 0.45, t_half 0-1.0 wks

### Slide 3: Resolving the PD Problem

- Original: 24.7-63.0% of covariance matrices silently corrected
- Root cause: OL+BDC expectancy-driven variance discontinuity
  under compound symmetry
- Solution: AR(1) + parameter constraints = zero PD failures
  across all 162 combinations
- Consequence: cross-design power comparisons are now
  unconfounded

### Slide 4: Resolving the Power Dropoff Problem

- Carryover erodes on/off-drug contrast through the variable
  chain: ondrug -> tod -> tsd -> Dbc -> bm:Dbc
- Original half-lives (0.1-0.2 wks) too short for genuine
  carryover
- Revised half-lives (0.5-1.0 wks) produce clinically
  meaningful power differences across designs
- Conditional formula: bm:Dbc (crossover) vs bm:t (open-label)

### Slide 5: Consolidated Framework

- Two validated pipelines: data.table (orig, gold standard) +
  tidyverse (2025, bit-identical)
- 67 unit tests, Docker reproducibility, zzcollab workspace
- 1,023-line tidyverse pipeline (38% reduction from cleanup)
- GitHub: rgt47/pmsimstats-ng

### Slide 6: Collaborator Deliverables

- analyze_trial_extended(): drug effect + interaction + predicted
  curves + clinical thresholds
- characterize_carryover(): sweep half-lives, compare AIC, find
  best fit
- 1,000-rep validated power with Type I error 2.5-5.8%
- 7 white papers documenting methodology

### Slide 7: Next Steps

- Goal 2b: per-timepoint residual analysis for carryover signal
  localization
- Final figure reproduction with consolidated code and revised
  defaults
- Manuscript describing the revised methodology and consolidated
  results
- Sensitivity analysis: carryover decay shape (exponential vs
  Weibull), separate correlation decay rate

---
*Rendered on 2026-03-25 at 08:31 PDT.*
*Source: ~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/03-audit-and-revision-report.md*
