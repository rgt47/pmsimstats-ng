# pmsimstats-ng Documentation Compendium
*2026-04-21 06:35 PDT*

**Author:** pmsimstats team

This compendium consolidates the novel technical content from all
source documents in `~/Dropbox/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/`.
It unifies mathematical derivations, empirical results, audit findings,
and implementation guidance that is currently scattered across twenty-six
numbered .md and .tex files. Each section cites the source documents that
contributed to it. Sections or claims that have been superseded by
subsequent development are flagged explicitly.

The documentation corpus was last substantially updated on 2026-04-16
(audit report) and 2026-04-15 (parity workflow, revision response to
reviewer on document 02). The current repository state reflects the
2026-04-08 merge that reorganized code into `implementations/`, added
Architecture A (mean moderation) alongside Architecture B (multivariate
normal differential correlation), and moved the pre-merge tidyverse
pipeline (`analysis/2025/`) to an external archive.

---

## 0. Source Document Inventory and Antiquation Map

The source files under `docs/` are listed below with their role and
their standing as of 2026-04-21. Files marked 'Current' are faithful
to the code. Files marked 'Historical' document decisions or states
that have been implemented; their technical content is preserved here
because the reasoning remains instructive. Files marked 'Superseded'
describe methodology that the revised code no longer uses; readers
should consult this compendium or the currently-labeled source rather
than those files.

| File | Role | Standing |
|---|---|---|
| 00-documentation-index.tex | Navigation index | Current |
| 01-codebase-overview.tex | Package structure and parameter table | Current |
| 02-revision-response-latent-class-expanded.md | Psychometric connection | Current |
| 02-revision-response-to-reviewer.md | Editorial revision log | Current |
| 03-audit-and-revision-report.tex | Seven DGP corrections | Current (foundational) |
| 04-revised-power-analysis.tex | Validated power estimates | Current |
| 06-ar1-residual-correlation.tex | nlme plus corCAR1 transition | Current |
| 07-positive-definiteness-failures.tex | PD failure analysis | Current |
| 08-biomarker-correlation-decay.tex | lambda_cor derivation | Current |
| 09-carryover-correlation-artifact.tex | Step-function artifact | Current (diagnostic) |
| 10-carryover-analysis-model-assessment.tex | 2x2 factorial plan | Current (proposed study) |
| 13-figure4-walkthrough.tex | Publication methodology | Superseded |
| 14-trial-analysis-implementation-guide.tex | Collaborator deliverables | Current |
| 15-tidyverse-alignment-plan.tex | Pre-merge port plan | Historical |
| 16-tidyverse-code-review.tex | Pre-merge code review | Historical |
| 17-simplification-plan.tex | Alternative simplification | Historical (not adopted) |
| 18-response-parameter-sensitivity.tex | Figure 5 validation | Current |
| 19-mean-moderation-implementation-notes.tex | Architecture A derivation | Current |
| 19-mean-moderation-implementation-notes.md | Quick-reference summary | Current |
| 20-annotated-bibliography.tex | Literature synthesis | Current reference |
| 21-architecture-comparison-summary.md | Empirical comparison | Current |
| 22-journal-target-recommendation.md | Submission strategy | Current |
| 22-tidyverse-development-parity-workflow.md | Parity workflow | Current |
| audit-2026-04-16.md | Technical audit | Current (most recent) |
| drafting-02-rgt/working.tex | Working draft of doc 02 | Working file (no material differences) |
| 07a-pd-constraints-theory.tex | PD constraints derivation (Sylvester, Gershgorin, correlation hierarchy) | Current (migrated 2026-04-23 from pmsimstats2025) |
| 09a-carryover-bias-variance.md | Non-monotonic power under model misspecification (bias-variance tradeoff) | Current (migrated 2026-04-23 from pmsimstats2025) |
| 24-component-decomposition-pedagogy.md | BR/ER/TV conceptual guide with prazosin-PTSD narrative | Current (migrated 2026-04-23 from pmsimstats2025) |
| archive/correlation-structure-alternatives.md | Alternative correlation structures (Matern, rational-quadratic, power-law) | Historical (migrated 2026-04-23; exploratory) |
| archive/hendrickson-pseudocode-comparison.md | Side-by-side pseudocode audit trail vs. Hendrickson original | Historical (migrated 2026-04-23) |

Three key antiquation points to highlight:

1. **Doc 13 (Figure 4 Walkthrough)** describes the publication
   methodology (compound symmetry, `lme4::lmer`, step-function BM-BR
   correlation, carryover half-lives of 0.1 and 0.2 weeks). The revised
   code no longer uses any of these elements. The document is retained
   as a record of the baseline from which the seven DGP corrections
   (Section 2) were derived.

2. **Doc 15 (Tidyverse Alignment Plan)** describes a pre-merge roadmap
   for porting `analysis/2025/pm_functions.R` to parity with the
   original data.table code. The merge on 2026-04-08 executed most of
   this plan, moved the pre-merge tidyverse file to an external
   archive, and reorganized the code into `implementations/`. Function
   names and paths referenced in doc 15 no longer resolve.

3. **Doc 16 (Tidyverse Code Review)** reviews the same pre-merge
   file. Its findings informed the 2026-04-08 rewrite. Line numbers
   and function names in doc 16 refer to the archived file, not the
   current `implementations/tidyverse/R/functions.R`.

4. **Doc 17 (Simplification Plan)** proposed a substantially simpler
   framework (random intercept with independent errors, ANCOVA on
   phase-mean change scores, reduced parameter space). The merge took
   a different path: it added Architecture A rather than simplifying.
   The document is retained as a record of considered alternatives
   and contains the most thorough formal specification of six
   alternative analysis methods in the project.

One execution-blocking issue was carried forward at the time of writing
(audit-2026-04-16, finding A1):
`analysis/carryover_factorial/02-run-factorial-2025.R:20` sources
`analysis/2025/01-pm-functions.R`, which has been moved to an external
archive. The script fails at load time until the source target is
updated to `implementations/tidyverse/R/functions.R` or a copy from the
archive is restored.

---

## 1. The Dual Data-Generating Process Architectures

Source documents: 02-dgp-mean-moderation-vs-mvn.tex;
19-mean-moderation-implementation-notes.tex; 21-architecture-comparison-summary.md;
02-revision-response-latent-class-expanded.md.

### 1.1 Architecture A: Direct Mean Moderation

The biomarker enters the mean structure of the outcome as a
multiplicative scaling of the treatment effect:

$$Y_{it} = \mu_0 + \beta_1 \cdot D_{it} \cdot (1 + \beta_{bm} \cdot B_i)
  + \epsilon_{it}$$

The additive form is equivalent under centering of the biomarker:

$$Y_{it} = \beta_0 + \beta_1 D_{it} + \beta_2 B_i + \beta_3 D_{it} B_i
  + \epsilon_{it}$$

Properties:

- The interaction is deterministic given biomarker and treatment
  status. This is a population-level mean-structure claim, not a
  statement about individual outcomes.
- `beta_bm` is a population parameter, uniform across participants.
- The signal lives in the first moment of the outcome distribution.
- Under carryover, the ratio of drug effect between high- and
  low-biomarker participants is preserved. Power loss is modest.

Implementation in the pmsimstats-ng framework: after an MVN draw from
a covariance matrix that omits the BM-BR block, each participant's BR
factor at on-drug timepoints is shifted by a standardized biomarker
deviation:

```
BR_{it}^{mod} = BR_{it}^{base} + c.bm * (B_i - mu_B)/sigma_B * 1[on drug] * sigma_BR
```

The shift per one standard deviation of biomarker is `c.bm * sigma_BR`.
For `c.bm = 0.45` and `sigma_BR = 8` this is 3.6 BR units, calibrated
to match the conditional expectation under Architecture B.

### 1.2 Architecture B: Multivariate Normal Differential Correlation

The biomarker and the biologic response factor are jointly multivariate
normal, with a treatment-state-dependent correlation:

$$\begin{pmatrix} B_i \\ BR_{it} \end{pmatrix} \sim
  \mathrm{MVN}\left( \begin{pmatrix} \mu_B \\ \mu_{BR}(t, D_{it}) \end{pmatrix},
  \Sigma(D_{it}) \right)$$

$$\mathrm{Cor}(B_i, BR_{it}) = \begin{cases}
  c_{bm} & \text{if } D_{it} = 1 \\
  c_{bm} \cdot \exp(-\lambda \cdot t_{sd}) & \text{if } D_{it} = 0
\end{cases}$$

The conditional expectation on drug is:

$$E[BR_{it} \mid B_i = b, D_{it} = 1] = \mu_{BR}(t)
  + c_{bm} \cdot \frac{\sigma_{BR}}{\sigma_B} \cdot (b - \mu_B)$$

Off-drug, with carryover half-life `t_half` and
`lambda = ln(2) / t_half`:

$$E[BR \mid b, D=0, t_{sd}] = \mu_{BR,\text{off}}(t_{sd})
  + c_{bm} \cdot e^{-\lambda t_{sd}} \cdot \frac{\sigma_{BR}}{\sigma_B}
  \cdot (b - \mu_B)$$

Properties:

- The interaction is probabilistic: a tendency, not a deterministic
  scaling.
- `c_bm` is a correlation coefficient with `|c_bm| <= 1`.
- The signal lives in the second moment (the covariance structure).
- Carryover erodes the differential correlation. Power loss is
  substantial.

### 1.3 Mathematical Comparison

Architecture A yields an expected difference that scales cleanly with
drug exposure:

$$E[Y \mid b_1, D] - E[Y \mid b_2, D]
  = \beta_1 \cdot \beta_{bm} \cdot (b_1 - b_2) \cdot D$$

The ratio of effects between high- and low-biomarker participants is
preserved as `D` varies. If `D` is replaced by a carryover-decayed
continuous indicator, the difference scales by the same factor, but
the high/low ratio is unchanged.

Architecture B yields an on-drug difference:

$$E[BR \mid b_1, D=1] - E[BR \mid b_2, D=1]
  = c_{bm} \cdot \frac{\sigma_{BR}}{\sigma_B} \cdot (b_1 - b_2)$$

Off drug, the difference is multiplied by `exp(-lambda * t_sd)`. The
biomarker-specific signal erodes exponentially during washout.

### 1.4 Empirical Power Comparison

Setup (doc 21): N = 70, c.bm = 0.45, nlme::lme plus corCAR1 with
`bm:Dbc` interaction, 50 replicates. Half-lives in weeks.

Architecture B (MVN differential correlation):

| Design | t1/2 = 0 | t1/2 = 0.5 | t1/2 = 1.0 |
|---|---|---|---|
| N-of-1 | 0.82 | 0.64 | 0.50 |
| CO | 0.40 | 0.40 | 0.32 |
| OL+BDC | 0.62 | 0.38 | 0.24 |

Relative power loss at t1/2 = 1.0: N-of-1 39 percent, CO 20 percent,
OL+BDC 61 percent.

Architecture A (mean moderation):

| Design | t1/2 = 0 | t1/2 = 0.5 | t1/2 = 1.0 |
|---|---|---|---|
| N-of-1 | 0.74 | 0.72 | 0.68 |
| CO | 0.38 | 0.36 | 0.34 |
| OL+BDC | 0.62 | 0.58 | 0.54 |

Relative power loss at t1/2 = 1.0: N-of-1 8 percent, CO 11 percent,
OL+BDC 13 percent.

Monte Carlo standard error at power near 0.5 with 50 replicates is
approximately `sqrt(0.5 * 0.5 / 50) = 0.07`.

### 1.5 Why the Architectures Diverge

Two mechanisms operate simultaneously under Architecture B:

1. **Mean blurring** (shared with Architecture A). Off-drug BR means
   are inflated by residual carryover, compressing the on-drug versus
   off-drug treatment contrast. Both architectures experience this.

2. **Correlation erosion** (unique to Architecture B). The BM-BR
   correlation decays exponentially with `t_sd`. Because the
   interaction signal lives in the second moment, this directly
   weakens the signal. This mechanism has no analogue in Architecture
   A, where the interaction lives in the first moment.

Architecture A's interaction coefficient sees reduced `Dbc` contrast
under carryover, but the signal-to-noise ratio degrades only modestly
because the biomarker-by-drug signal scales proportionally with the
reduced contrast. Architecture B's interaction coefficient is attacked
from both sides: the treatment contrast in `Dbc` is compressed, and the
biomarker-BR correlation generating the covariance signal is
simultaneously decaying.

### 1.6 Biological Rationale and Architecture Selection

The choice between architectures is a claim about biology, not
statistical convenience.

Architecture A is appropriate when the biomarker governs the
**magnitude of the drug's biological effect** for each individual.
Examples:

- Renal clearance (eGFR determines drug elimination).
- Receptor density moderation (target availability determines
  pharmacodynamics).
- CYP2D6 metabolizer status.
- Dose-response with biomarker-dependent effective dose.

Architecture B is appropriate when the biomarker predicts which
participants are likely to respond, through shared variance with a
drug-responsive physiological component that decays as the drug is
eliminated. Examples:

- Blood pressure as PTSD subtype marker (elevated resting SBP indexes
  heightened central noradrenergic tone, a subtype preferentially
  responsive to prazosin as alpha-adrenergic antagonist).
- Baseline disease severity as treatment predictor.
- Inflammatory biomarkers in psychiatry.
- Genetic risk scores.

Decision table (doc 02):

| Criterion | Architecture A | Architecture B |
|---|---|---|
| Biomarker role | Causal mediator | Statistical predictor |
| Mechanism | Deterministic scaling | Probabilistic association |
| Signal location | Mean structure | Covariance structure |
| Carryover impact | Modest (10-15 pct) | Substantial (40-60 pct) |
| Appropriate when | Biomarker determines effective dose or PK | Biomarker indexes latent subtype |

Reporting recommendation: simulation studies should explicitly state
whether the biomarker-treatment interaction is mean moderation or
differential correlation, provide the biological rationale, confirm
that the carryover model operates consistently with the chosen
architecture, and include sensitivity analyses under the alternative
architecture when the biological mechanism is uncertain. When
uncertain, Architecture B produces more conservative (lower) power
estimates under carryover and is therefore the safer default for
trial design.

### 1.7 Prevalence of Each Architecture in the Literature

Architecture A is the near-universal standard in clinical trial
simulation, appearing in biomarker-stratified designs (Simon 2010),
adaptive enrichment trials (Freidlin and Korn 2014), basket and
umbrella trials (Renfro and Sargent 2017), platform trials, and
biomarker-by-treatment interaction studies (Haller and Ulm 2018).
It is standard in N-of-1 simulation (Zucker et al. 1997; Araujo et
al. 2016; Duan et al. 2013) and in N-of-1 methodology on power and
design (Wang and Schork 2019; Schork 2022).

Architecture B (MVN differential correlation) appears primarily in
latent variable and structural equation modeling, Bayesian joint
modeling, and copula-based simulation. The Hendrickson et al. (2020)
framework appears to be unique in the N-of-1 trial literature in
using differential correlation as the interaction mechanism.

Implication: published power estimates for biomarker-moderated
crossover and N-of-1 designs may be systematically optimistic for
biomarkers that operate through an Architecture B mechanism, because
those studies do not account for the additional power loss that
carryover produces when the interaction signal resides in the
covariance structure rather than the mean structure.

---

## 2. The Seven DGP Corrections

Source document: 03-audit-and-revision-report.tex (dated 2026-03-24).

The revised codebase differs from the Hendrickson et al. (2020)
publication code in seven substantive ways. Each is documented below
with the mathematical form, the reason for the change, and the effect
on inference.

### 2.1 Correction 1: AR(1) Autocorrelation

Replaces compound symmetry within each response factor.

$$\mathrm{Cor}(X_{t_i}, X_{t_j}) = \rho^{|t_i - t_j|}$$

Under compound symmetry the correlation is constant regardless of
time gap. Under AR(1) it decays exponentially with temporal
separation, which is clinically more realistic. The change expands
the PD-feasible maximum for the BM-BR correlation from 0.25 under
compound symmetry with `rho = 0.8` to 0.49 or higher under AR(1)
with `rho = 0.7`.

### 2.2 Correction 2: Exponential BM-BR Correlation Decay

Replaces the step-function rule under which the BM-BR correlation
dropped to zero the moment a participant was off drug regardless of
carryover residual.

$$\mathrm{Cor}(\mathrm{bm}, \mathrm{br}_t) = c_{bm} \cdot
  e^{-\lambda_{\mathrm{cor}} \cdot t_{\mathrm{sd}}}$$

with `lambda_cor = ln(2) / t_half` by default. The correlation decays
in proportion to the drug's pharmacokinetic elimination rate. The
derivation from first principles appears in Section 5 below.

### 2.3 Correction 3: nlme::lme with corCAR1

Replaces `lme4::lmer` with random intercept only.

```
nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc,
          random = ~1|ptID,
          correlation = nlme::corCAR1(form = ~t|ptID),
          data = datamerged,
          control = nlme::lmeControl(opt = 'optim',
                                     maxIter = 200,
                                     msMaxIter = 200))
```

Under the AR(1) DGP, `lmer` with random intercept alone inflates
Type I error to 13-17 percent for CO and OL designs. The `corCAR1`
residual structure restores nominal Type I control (Section 4
below).

### 2.4 Correction 4: Carryover Scale Factor Removal

Removes the undocumented scale factor of 2 that commit 8609f12
added to the carryover decay exponent:

- Original (commit 42ac030): `(1/2)^(t_sd / t_half)`.
- Ron Thomas edit (commit 8609f12): `(1/2)^(2 * t_sd / t_half)`,
  halving the effective half-life.
- Revised HEAD: `(1/2)^(t_sd / t_half)`. The parameter `t_half`
  means what it says.

### 2.5 Correction 5: Timepoint-1 Guard Removal

Commit 8609f12 introduced a guard `if (p > 1)` in the BM-BR
correlation assignment loop that was intended to avoid an index
error in the proportional-scaling formula. The guard inadvertently
forced the BM-BR correlation at timepoint 1 to zero even when
participants were on drug at that timepoint. The revised code
iterates from `p = 1` without the guard.

### 2.6 Correction 6: Gompertz Origin-Passing

Replaces the naive parameterization `f(t) = M * exp(-d * exp(-r*t))`,
which does not pass through the origin, with an offset-rescaled form
that guarantees `f(0) = 0` and preserves the asymptote at `M`:

```
y <- maxr * exp(-disp * exp(-rate * t))
vert_offset <- maxr * exp(-disp)
y <- y - vert_offset
y <- y * (maxr / (maxr - vert_offset))
```

Doc 15 (tidyverse alignment plan) identifies this correction as the
single largest source of divergence between the original and the
pre-merge tidyverse implementation. The naive form shifts every
component mean in the MVN distribution.

### 2.7 Correction 7: Updated Default Parameters

- Autocorrelation: 0.8 (compound symmetry) to 0.7 (AR(1)).
- Biomarker correlation: publication grid `{0, 0.3, 0.6}` reduced
  to `{0, 0.25, 0.45}`. PD-feasible maximum under AR(1) is 0.49
  for OL; `c.bm = 0.6` fails PD checks for 7 of 9 trial paths at
  zero carryover and all 9 at higher carryover under Architecture B.
- Carryover half-life: publication grid `{0, 0.1, 0.2}` weeks
  replaced with `{0, 0.5, 1.0}` weeks. The publication values are
  too short to be clinically meaningful and were used precisely
  where the step-function correlation artifact was largest.

### 2.8 Correction Status Across Implementations

From the 2026-04-16 audit:

| Correction | original | original-extended | tidyverse | simple |
|---|---|---|---|---|
| 1. AR(1) autocorrelation | yes | yes | yes | n/a |
| 2. Exponential BM-BR decay | yes | yes | yes | n/a |
| 3. nlme::lme with corCAR1 | yes | yes | yes | no (OLS) |
| 4. Scale factor removal | yes | yes | yes | n/a |
| 5. Timepoint-1 guard removal | yes | yes | yes | n/a |
| 6. Gompertz origin-passing | yes | yes | yes | partial |
| 7. Updated defaults | yes | yes | yes | partial |

All seven corrections are embodied in the three production
implementations. The simple sandbox omits some by design (it uses
ANCOVA on phase-mean change scores rather than LME).

---

## 3. Positive Definiteness of the Covariance Matrix

Source documents: 07-positive-definiteness-failures.tex;
03-audit-and-revision-report.tex; 04-revised-power-analysis.tex.

### 3.1 The Problem

The DGP builds a 26-by-26 covariance matrix per trial path from
user-specified correlations and standard deviations. When the result
is not positive definite, `generateData()` projects to the nearest PD
matrix via `corpcor::make.positive.definite(sigma, tol = 1e-3)`.
In the publication code the correction was silent.

Failure rates by commit:

| Commit | Description | Failures | Rate |
|---|---|---|---|
| 42ac030 | Initial publication | 40 of 162 | 24.7 pct |
| 8609f12 | Ron Thomas edits | 102 of 162 | 63.0 pct |
| f6ee86d | Autocorr typo fix | 102 of 162 | 63.0 pct |
| f70f86d | Carryover added | 102 of 162 | 63.0 pct |

All 162 matrices in the standard parameter sweep were ill-conditioned
(`kappa > 100`) under compound symmetry.

Condition number statistics (all 162 matrices):

| Commit | Median kappa | kappa > 100 | kappa > 1e3 | kappa > 1e6 |
|---|---|---|---|---|
| 42ac030 | 229 | 162/162 | 40/162 | 40/162 |
| 8609f12 | 3.1 times 10^18 | 162/162 | 102/162 | 102/162 |

Minimum eigenvalues of failed matrices ranged from -33.26 (worst) to
-0.52 (least severe), on a dominant eigenvalue of approximately 700.
The most extreme corrections displaced about 5 percent of the matrix
structure, concentrated in eigenvectors affecting the BM-BR and
cross-phase ER correlations disproportionately.

### 3.2 Root Cause: OL+BDC Expectancy Discontinuity

Failures were concentrated in the OL+BDC design (83-92 percent
failure rate, 32 of 40 failures under the initial commit). The OL
design never fails. The mechanism:

- Expectancy changes abruptly: `e = 1.0` in the open-label phase
  (4 timepoints) and `e = 0.5` in the blinded phase (4 timepoints).
- The ER standard deviation scales with expectancy:
  `sigma_ER,t = sigma_ER_base * e_t`, producing the pattern
  `[10, 10, 10, 10, 5, 5, 5, 5]`.
- Cross-phase covariance under compound symmetry with `rho = 0.8`:
  `0.8 * 10 * 5 = 40`, yielding a pairwise correlation of
  `40 / sqrt(100 * 25) = 0.8`, consistent at the pair level.
- However, compound symmetry requires all same-factor timepoint
  pairs to have identical correlation regardless of temporal
  separation. With heterogeneous variances across phases this
  becomes self-contradictory: a single correlation coefficient
  cannot simultaneously satisfy the within-phase structure in both
  phases and the cross-phase structure.

OL never fails because expectancy is constant (`e_t = 1`) and all
variances are homogeneous within each factor block.

### 3.3 Resolution

Two mechanisms restore PD:

1. **AR(1) autocorrelation replaces compound symmetry.** Under AR(1),
   nearby timepoints are more correlated than distant ones, which
   accommodates variance heterogeneity without forcing a single
   correlation to cover all pairs. Maximum feasible `c.bm` increases
   from 0.25 (CS with `rho = 0.8`) to 0.49 or higher (AR(1) with
   `rho = 0.7`) depending on design.

2. **Pre-validation via `validateParameterGrid()`.** The revised
   code tests every parameter combination for PD before simulation
   begins and reports which cells would require correction. Under
   the revised parameter grid (`c.bm <= 0.45`, `rho = 0.7`,
   clinically realistic half-lives), zero of 162 matrices require
   correction.

Maximum feasible `c.bm` by design:

| Design | Original (CS, rho = 0.8) | Revised (AR(1), rho = 0.7) |
|---|---|---|
| OL | 0.89 | 0.49 |
| OL+BDC | 0.26 | 0.45 |
| CO | 0.25 | 0.61 |
| N-of-1 | 0.25 | 0.53 |

The publication value of `c.bm = 0.6` exceeds the PD-feasible range
for every design except OL. The revised value of `c.bm = 0.45` is
within range for all designs.

### 3.4 Architecture A Removes the Constraint

Under Architecture A, the BM-BR correlation block is not entered
into the covariance matrix (it is applied post-draw as an additive
shift). The `c.bm` parameter therefore has no effect on the
eigenvalue structure. Testing `c.bm` from 0 to 0.95 in steps of 0.05
across 4 designs and 3 half-lives:

| | Architecture B | Architecture A |
|---|---|---|
| First PD failure | c.bm = 0.50 | none through 0.95 |
| Feasible ceiling | approximately 0.45 | unconstrained |
| Paths failing at 0.60 | 7 of 9 | 0 of 9 |
| Paths failing at 0.95 | 9 of 9 | 0 of 9 |

### 3.5 Residual PD Instrumentation Gap (Finding C1)

The audit-2026-04-16 flagged that all three production implementations
still call `make.positive.definite()` silently when a matrix is not PD
(root `generateData.R:284-299`; corresponding lines in original,
original-extended, tidyverse). No flag is exposed to downstream
consumers. The current PD correction rate therefore cannot be verified
without instrumentation. The recommended fix is to return a
`pd_corrected` flag from `buildSigma()` and aggregate it into
per-parameter-set summaries.

The doc 19 (mean moderation implementation notes) also flags that
`validateParameterGrid()` builds sigma matrices with the Architecture B
BM-BR correlation block regardless of the simulation's actual
`dgp_architecture` setting. Under Architecture A with high `c.bm`, the
validator may report spurious PD failures. The data generation itself
is unaffected because the simulation uses the Architecture A-specific
sigma. The warning can be safely ignored when running Architecture A.

---

## 4. Analysis Model: nlme::lme with corCAR1

Source documents: 06-ar1-residual-correlation.tex;
03-audit-and-revision-report.tex; 04-revised-power-analysis.tex.

### 4.1 Formula Construction

The analysis model adapts to the design. For designs with
within-subject drug variation (CO, N-of-1, OL+BDC):

$$S_{i,t} = \beta_0 + \beta_1 \mathrm{bm}_i + \beta_2 t
  + \beta_3 D_{bc,t} + \beta_4 \mathrm{bm}_i \cdot D_{bc,t}
  + u_i + \varepsilon_{i,t}$$

with `Dbc` the continuous drug indicator (1 on drug, exponential
decay off drug) and `beta_4` the interaction of interest.

For open-label (no within-subject drug variation):

$$S_{i,t} = \beta_0 + \beta_1 \mathrm{bm}_i + \beta_2 t
  + \beta_3 \mathrm{bm}_i \cdot t + u_i + \varepsilon_{i,t}$$

with `beta_3` the biomarker-by-time interaction. Random intercept
only (`random = ~1|ptID`). Residual correlation
`corCAR1(form = ~t|ptID)`.

### 4.2 Continuous-Time AR(1) Rationale

Trial timepoints are unequally spaced. The N-of-1 design uses weeks
`c(4, 8, 9, 10, 11, 12, 16, 20)` with gaps of 1 to 4 weeks. Discrete
AR(1) (`corAR1`) assumes equal spacing and is incorrect for these
designs. Continuous-time AR(1) (`corCAR1`) handles unequal spacing
via `phi^|t_1 - t_2|` where `t` is in weeks.

### 4.3 Why the Random Intercept Alone Is Insufficient

Under compound symmetry the within-person correlation is constant,
and a random intercept absorbs it perfectly (residuals independent).

Under AR(1):

$$\mathrm{Cor}(\varepsilon_{i,t_1}, \varepsilon_{i,t_2} \mid u_i)
  = \rho^{|t_1 - t_2|}
  - \frac{\sigma_u^2}{\sigma_u^2 + \sigma_\varepsilon^2}$$

The random intercept absorbs only the average correlation. Residual
correlation remains positive for nearby timepoints (where
`rho^|dt|` exceeds the average) and negative for distant timepoints.
`lme4::lmer` treats residuals as independent and therefore
underestimates standard errors.

### 4.4 Type I Error Correction

Under AR(1) DGP with `lme4::lmer` analysis (200 replicates, `c.bm = 0`,
no censoring, no carryover):

| Design | N=35 | N=70 |
|---|---|---|
| OL | 0.13 | 0.09 |
| OL+BDC | 0.09 | 0.05 |
| CO | 0.17 | 0.13 |
| N-of-1 | 0.06 | 0.08 |

The CO design reaches 17 percent false positives at N=35, 3.4 times
nominal. CO and OL are most affected because their drug indicators
are most strongly confounded with time-dependent residual structure.

Under the revised `nlme::lme` plus `corCAR1` model:

| Design | N=35 | N=70 |
|---|---|---|
| OL | 0.03 | 0.02 |
| OL+BDC | 0.03 | 0.04 |
| CO | 0.06 | 0.08 |
| N-of-1 | 0.04 | 0.06 |

All within the 95 percent binomial CI for `p = 0.05` at n = 200
(approximately [0.02, 0.09]). At 1,000 replicates (doc 04), Type I
error is 2.5 to 5.8 percent across all designs.

### 4.5 Error Handling

`lme_analysis()` wraps the fit in a two-level fallback:

1. If `lme` with `corCAR1` fails to converge: fall back to `lme`
   without correlation structure.
2. If the fallback fails: return NA for all outputs.

Rows with NA outcomes are removed before fitting (nlme uses
`na.fail` by default). Audit finding C2 notes that the returned
`issingular` column is currently hard-coded FALSE
(`lme_analysis.R:178, 221`); downstream consumers cannot
distinguish clean fits from singular ones. The recommended fix is
to inspect `fit$apVar` for the 'Non-positive-definite approximate
variance-covariance' string and check random-effects variance for
zero eigenvalues.

### 4.6 Coefficient Extraction

nlme returns `summary(fit)$tTable` with columns `Value`, `Std.Error`,
`p-value` (contrasted with `lme4::lmer`, which uses `Estimate`,
`Std. Error`, `Pr(>|t|)`). The extraction searches the row names for
`'bm:Dbc'` or `'Dbc:bm'` and selects the first match:

```r
target <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(c))
if (length(target) == 0) {
  beta <- betaSE <- p <- as.numeric(NA)
} else {
  target <- target[1]
  p <- c[target, 'p-value']
  beta <- c[target, 'Value']
  betaSE <- c[target, 'Std.Error']
}
```

Audit finding C4 flags this as fragile; `carryover_analysis.R`
uses a more flexible pattern with a longer alias list.

---

## 5. BM-BR Correlation Decay: Derivation from First Principles

Source documents: 08-biomarker-correlation-decay.tex;
09-carryover-correlation-artifact.tex.

### 5.1 Three Historical Rules

The BM-BR correlation at an off-drug timepoint has been computed
three different ways across commits:

1. **Step function** (commit 42ac030, publication): full correlation
   `c_bm` whenever the adjusted BR mean is nonzero, zero otherwise.
   `rho_t = c_bm * 1[mu_BR,t != 0]`.

2. **Proportional scaling** (commit 8609f12, Ron Thomas):
   `rho_t = (mu_t / mu_{t-1}) * c_bm`. Correlation decays
   proportionally with the ratio of adjacent-timepoint BR means.

3. **Independent exponential** (current HEAD):
   `rho_t = c_bm * exp(-lambda_cor * t_sd)`. Correlation decays
   exponentially with a parameterized rate.

### 5.2 Individual-Level Argument

Consider two participants during an on-drug period when the
population BR mean is `mu`:

- Participant A with biomarker one SD above mean.
- Participant B with biomarker one SD below mean.

From the bivariate normal conditional expectation:

$$E[\mathrm{BR}_A] = \mu + c_{bm} \cdot \sigma_{BR}$$
$$E[\mathrm{BR}_B] = \mu - c_{bm} \cdot \sigma_{BR}$$

Biomarker-specific spread:
`Delta_on = 2 * c_bm * sigma_BR`.

### 5.3 Behavior at Drug Discontinuation

If carryover is pharmacokinetic, drug molecules clear at a rate
determined by the half-life, independent of any participant's
biomarker. Every participant's drug effect decays by the same factor
`alpha = (1/2)^(t_sd / t_half)`:

$$E[\mathrm{BR}_A^{\mathrm{off}}] = \alpha (\mu + c_{bm} \sigma_{BR}),
  \quad E[\mathrm{BR}_B^{\mathrm{off}}] = \alpha (\mu - c_{bm}
  \sigma_{BR})$$

Population mean off-drug: `mu_off = alpha * mu`.

Biomarker-specific spread off-drug: `Delta_off = alpha * Delta_on`.
The spread decays by the same factor as the mean.

### 5.4 Implied Correlation

Setting the implied spread equal to the MVN-derived form
`Delta_off = 2 * rho_off * sigma_BR`:

$$\boxed{\; \rho_{\mathrm{off}} = \alpha \cdot c_{bm}
  = (1/2)^{t_{sd}/t_{1/2}} \cdot c_{bm} \;}$$

Equivalently in exponential form, since
`alpha = exp(-(ln 2 / t_half) * t_sd)`:

$$\boxed{\; \lambda_{\mathrm{cor}} = \lambda_{\mathrm{drug}}
  = \frac{\ln 2}{t_{1/2}} \;}$$

The correlation decay rate should equal the drug elimination rate
constant.

### 5.5 Numerical Validation

N-of-1 Path A, `c.bm = 0.5`, `t_half = 1.0` week:

| Metric | On drug (BD2, week 10) | Off drug (BD3, week 11) |
|---|---|---|
| mu_BR | 10.19 | 5.09 |
| rho | 0.50 | 0.25 |
| E[BR | +1 SD] | 14.19 | 7.09 |
| E[BR | -1 SD] | 6.19 | 3.09 |
| Spread | 8.0 | 4.0 |
| Coefficient of variation | 0.785 | 0.786 |

The spread halves as the mean halves. The coefficient of variation
of the biomarker-specific component is preserved at 0.79. No
participant has negative expected BR.

### 5.6 Comparison with Alternative Rules

Off-drug conditional expectations at BD3 (`mu_BR = 5.09`,
`c.bm = 0.5`):

| Rule | rho | E[BR | +1 SD] | E[BR | -1 SD] | Issue |
|---|---|---|---|---|
| Step function (lambda_cor = 0) | 0 | 5.09 | 5.09 | No biomarker signal |
| Proportional (lambda_cor = lambda_drug) | 0.25 | 7.09 | 3.09 | Consistent |
| Slow decay (lambda_cor = 1) | 0.18 | 6.57 | 3.61 | Under-predicts spread |
| Full correlation (rho = c_bm) | 0.50 | 9.09 | 1.09 | Over-predicts spread |

Only proportional scaling produces a spread consistent with
pharmacokinetic decay.

### 5.7 The Publication Artifact

The publication used the step-function rule at `t_half = 0.1` and
`t_half = 0.2` weeks. At these values the carryover residual in the
mean is less than 0.1 percent of peak. The step-function rule
assigned the full `c_bm` correlation to these timepoints, producing
large biomarker-specific variation around a near-zero mean.
Numerical example at BD3 with 0.1-week half-life, `c.bm = 0.6`,
`sigma_BR = 8`, `sigma_BM = 15.36`:

| Biomarker | E[BR] original | E[BR] revised |
|---|---|---|
| 109 (-1 SD) | -4.8 | +0.01 |
| 124 (mean) | +0.01 | +0.01 |
| 140 (+1 SD) | +4.9 | +0.01 |

Under the original rule, a low-biomarker participant's conditional
expected BR is -4.8 at a timepoint where the drug effect has decayed
to 0.01 percent of peak. This is clinically incoherent: it implies
that discontinuing the drug worsens symptoms beyond the natural
trajectory for low-biomarker participants, specifically because the
drug is gone. No rebound effect was hypothesized.

The analysis model sees this spurious biomarker-specific variation
as structured noise in the `bm:Dbc` term. It inflates the residual
variance of the interaction estimate, widens the standard error of
`beta_4`, and reduces the probability of rejecting the null. The
effect is largest for designs with many off-drug timepoints (N-of-1
has 4 of 8) and smallest for OL (always on drug).

Under the revised rule, the off-drug correlation at BD3 with a
0.1-week half-life is approximately `0.001 * c_bm`, effectively
zero. The Hendrickson et al. (2020) finding of a power drop from
0.95 to 0.18 at `t_half = 0.1` weeks is therefore primarily a
modeling artifact. At clinically realistic half-lives (0.5 to 1.0
weeks), the revised correlation is substantial and carryover does
reduce power, but moderately rather than catastrophically.

### 5.8 The Remaining Variance Limitation

The derivation assumes `sigma_BR = 8` is constant. This is the
current implementation. The conditional variance during off-drug
periods is nearly unchanged (`sigma_BR^2 (1 - rho_t^2)` with small
`rho_t`), so participants have large random fluctuations around a
decaying mean. This conflates two sources of variation that should
be separated:

- **Trait variability**: between-person differences in latent drug
  responsiveness, existing before drug and persisting regardless of
  drug status. Should remain constant.
- **State-dependent variability**: variation in acute drug effect,
  scaling with effect magnitude. Should decay during carryover.

Addressing the limitation would require time-varying `sigma_BR,t`
or a variance decomposition, which would be substantial restructuring.
For the current application of comparing relative power across designs,
the fixed-variance assumption affects all designs similarly and is
acceptable as a conservative approximation.

---

## 6. Validated Results

Source documents: 04-revised-power-analysis.tex;
06-ar1-residual-correlation.tex; 14-trial-analysis-implementation-guide.tex.

### 6.1 Revised Parameter Grid

| Parameter | Value | Rationale |
|---|---|---|
| dgp_architecture | 'mvn' (default); 'mean_moderation' optional | Dual support |
| Autocorrelation | AR(1), rho = 0.7 | PD-feasible, clinically realistic |
| c.bm (Arch B) | {0, 0.25, 0.45} | Max PD-feasible approximately 0.45 for OL+BDC |
| c.bm (Arch A) | up to 0.95 | PD constraint absent |
| Carryover t_half | {0, 0.5, 1.0} weeks | Clinically realistic |
| lambda_cor | ln(2) / t_half (auto) | Pharmacokinetic derivation |
| Cross-factor (same t) | c.cf1t = 0.2 | Unchanged from publication |
| Cross-factor (diff t) | c.cfct times rho^|t_i - t_j| | AR(1) decay |
| Analysis model | nlme::lme + corCAR1 | Required for Type I control |
| PD corrections | 0 of 162 | All matrices valid |

### 6.2 Type I Error at 1,000 Replicates (c.bm = 0)

| Design | N=35 | N=70 |
|---|---|---|
| OL | 0.029 | 0.025 |
| OL+BDC | 0.030 | 0.031 |
| CO | 0.053 | 0.058 |
| N-of-1 | 0.040 | 0.030 |

### 6.3 Power at c.bm = 0.45, No Censoring, 1,000 Replicates

| Design | N | t_half = 0 | t_half = 0.5 | t_half = 1.0 | Interpretation |
|---|---|---|---|---|---|
| OL | 35 | 0.06 | 0.06 | 0.06 | Unaffected |
| OL | 70 | 0.09 | 0.10 | 0.09 | Unaffected |
| CO | 35 | 0.19 | 0.21 | 0.18 | Robust |
| CO | 70 | 0.36 | 0.35 | 0.30 | Robust |
| OL+BDC | 35 | 0.31 | 0.21 | 0.12 | Moderate drop |
| OL+BDC | 70 | 0.62 | 0.44 | 0.21 | Large drop |
| N-of-1 | 35 | 0.45 | 0.33 | 0.23 | Large drop |
| N-of-1 | 70 | 0.75 | 0.63 | 0.43 | Large drop |

The qualitative design ranking from the publication is preserved
(N-of-1 greater than CO greater than OL+BDC greater than OL at
low N, with ordering adjustments at high N). Carryover sensitivity
is moderate and gradual, not the catastrophic drop from 0.95 to
0.18 at `t_half = 0.1` reported in the publication.

### 6.4 Runtime

Current implementation runtime on modern hardware (from doc 14):

- 200 replicates: approximately 8 minutes.
- 1,000 replicates: approximately 47 minutes.

Contributing optimizations (from doc 04 appendix): sigma caching
(pre-build 162 covariance matrices once), `furrr::future_map`
parallel execution across available cores, pre-computed Cholesky
factor, vectorized `data.table` outcome construction, pre-allocated
output collection, PD pre-validation, automatic core detection.
Together these yield roughly a 12x speedup over the publication code.

---

## 7. Repository Structure and Implementations

Source documents: 01-codebase-overview.tex; audit-2026-04-16.md;
and implementations/README.md.

### 7.1 Four Parallel Implementations

| Collection | Style | Architectures | Location | Role |
|---|---|---|---|---|
| original | data.table | B (MVN) only | `implementations/original/` | Historical baseline |
| original-extended | data.table | A + B | `implementations/original-extended/` | Production |
| tidyverse | tidyverse | A + B | `implementations/tidyverse/` | Modern alternative |
| simple | base R + tidyverse | A only | `implementations/simple/` | Pedagogical sandbox |

The installed package at root `R/` is effectively a copy of
`original-extended` with additional plotting, carryover analysis, and
documentation functions. The seven core R files are byte-identical
between root and original-extended except for the `dgp_architecture`
threading in `generateData.R` and `generateSimulatedResults.R`.

### 7.2 Installed Package File Inventory (root R/)

| File | Lines | Exports |
|---|---|---|
| buildtrialdesign.R | 119 | `buildtrialdesign` |
| censordata.R | 89 | `censordata` |
| datadocumentation.R | 61 | (none; roxygen stubs) |
| generateData.R | 355 | `generateData`, `buildSigma`, `validateParameterGrid` |
| generateSimulatedResults.R | 330 | `generateSimulatedResults` |
| lme_analysis.R | 237 | `lme_analysis` |
| packagedocumentation.R | 19 | (package NULL) |
| plottingfunctions.R | 341 | `PlotModelingResults`, `plotfactortrajectories` |
| utilities.R | 108 | `cumulative`, `modgompertz`, `reknitsimresults` |
| carryover_analysis.R | 377 | `characterize_carryover`, `analyze_trial_extended`, `print_carryover_summary`, `print_trial_summary` |

### 7.3 Where Each Implementation Stands

**original**: the seven-file data.table code as corrected in the
revised repository. Architecture B only; no `dgp_architecture`
parameter. Supports OL, CO, OL+BDC, and N-of-1 designs. Has an
independent test suite (6 test files, 458 lines) using tinytest
naming conventions, driven by a custom `run_tests.R` rather than
standard package test infrastructure. Retained as the
Architecture B-only reference implementation.

**original-extended**: adds the `dgp_architecture` argument with
`match.arg` enforcement, the Architecture A additive shift block in
`generateData()`, and the BM-BR guard in `buildSigma()`. Threads the
argument through `generateSimulatedResults()`. Fifty-two additional
lines in two files relative to `original`. The `tests/` directory is
empty; coverage relies on the root `inst/tinytest/` suite and the
parity tests.

**tidyverse**: a complete reimplementation as a single 991-line
`functions.R` with snake_case naming, tidyverse selection helpers,
and three forms of carryover decay (exponential, linear, Weibull)
exposed via `carryover_decay()`. The `tests/` directory is empty.
A non-automated `test-alignment.R` lives alongside `functions.R`
and is sourced manually to validate sigma, Gompertz, and data
generation against the `original` collection.

**simple**: a single 577-line `simulation.R` that demonstrates
Architecture A mean moderation with a reduced-complexity model
(ANCOVA on phase-level means rather than LME plus corCAR1; no sigma
caching; no carryover-aware residual structure). Not package-
structured. By design, not used for publication-grade power
calculations.

### 7.4 Test Inventory

From the 2026-04-16 audit:

| Location | Files | Coverage | Assertions | Gaps |
|---|---|---|---|---|
| `inst/tinytest/` | 5 | cumulative, modgompertz, buildtrialdesign, buildSigma, generateData, lme_analysis | 57 | No Type I error check; no mean_moderation; no corCAR1 correctness; no PD recovery; no convergence/singular fit test |
| `implementations/original/tests/` | 6 | independent | 458 lines | Tinytest-style naming but sourced manually |
| `implementations/original-extended/tests/` | 0 | none | 0 | Empty |
| `implementations/tidyverse/R/test-alignment.R` | 1 | sigma, Gompertz, data parity | 5 stopifnot | Manual sourcing; no formal framework |
| `implementations/test-parity-extended-tidyverse.R` | 1 | cross-implementation | approximately 144 cells | Parity-only |

The root `inst/tinytest/` suite covers core functions but has
known coverage gaps flagged as audit finding C7: no test exercises
`dgp_architecture = 'mean_moderation'`, null-effect Type I error,
`corCAR1` correctness under known `phi` ground truth, or the silent
PD correction path.

Note: the project CLAUDE.md refers to `testthat` framework. The
installed package actually uses `tinytest` (configured as
`Suggests: tinytest (>= 1.4.0)` in DESCRIPTION, driven by
`tests/tinytest.R`). This discrepancy is flagged in the audit.

### 7.5 Bundled Data

- `results_core.rda`: publication Figure 4 results (initial commit,
  500 replicates).
- `extracted_bp.rda`: baseline parameters (BP mean/SD, CAPS mean/SD
  from pilot data).
- `extracted_rp.rda`: response parameters (Gompertz max, disp, rate,
  sd for tv, pb, br).
- `CTdata.rda`: actual clinical trial data for Vignette 3.
- `results_maxes.rda`, `results_rates.rda`, `results_trajectories.rda`:
  Figure 5 results.

### 7.6 Reproducibility Infrastructure

- **Docker**: `rocker/tidyverse:4.5.3` (ARG `R_VERSION=4.5.3`);
  non-root `analyst` user; `languageserver` for IDE support. Working
  directory `/home/analyst/project`.
- **renv.lock**: R 4.5.2 (4.4.2 in the CI workflow), approximately
  200 packages pinned from Posit PPM.
- **Makefile targets**: `help`, `r`, `rstudio`, `document`, `build`,
  `check`, `install`, `test`, `vignettes`, `manuscripts`, `deps`,
  `check-renv`, `check-system-deps`, `docker-build`, `docker-test`,
  `manuscript-01`, `manuscript-02`.
- **GitHub Actions**: `.github/workflows/r-package.yml` builds in
  `rocker/tidyverse:4.4.2`, restores from renv.lock, and runs
  `R CMD check --no-manual --ignore-vignettes` with `rcmdcheck` and
  `tinytest` installed to a separate CI tool library.
- **Pre-commit hook**: `scripts/pre-commit-parity.sh` runs
  `implementations/test-parity-original-tidyverse.R` to validate
  cross-implementation parity on each commit.

---

## 8. Architecture A: Implementation History

Source documents: 19-mean-moderation-implementation-notes.tex;
19-mean-moderation-implementation-notes.md;
02-revision-response-to-reviewer.md.

### 8.1 Four Attempts

The current Architecture A implementation was reached after three
incorrect formulations. Each is recorded because the failure modes
illuminate the calibration constraints the correct formulation must
satisfy.

**Attempt 1: Multiplicative scaling (incorrect).**

```
BR_{it}^{mod} = BR_{it}^{base} * (1 + beta_bm * B_i^*)
```

Problem: during off-drug timepoints, `BR_{it}^{base}` approximately
zero (the Gompertz evaluates to approximately 0 when `tod = 0`).
Multiplying near-zero by `(1 + small)` produces near-zero. The
signal is proportional to the magnitude of the base BR, not to the
drug exposure. Result: power dropped approximately 60 percent under
carryover, indistinguishable from Architecture B. Lesson:
multiplicative moderation conflates the interaction signal with the
response magnitude.

**Attempt 2: Additive with drug exposure decay (incorrect).**

```
BR_{it}^{mod} = BR_{it}^{base} + beta_bm * B_i^* * phi(t) * sigma_BR
with phi(t) = 1 on drug,
     (1/2)^(t_sd / t_half) off drug,
     0 if never treated
```

Problem: drug exposure decays during off-drug. At 2.5-week intervals
with `t_half = 1.0`, exposure values are {0.18, 0.031, 0.006, 0.001}.
Moderation signal at the first off-drug timepoint is
`0.45 * 1 * 0.18 * 8 = 0.65` BR units, falling to 0.004. Against
per-timepoint noise `sigma_BR = 8`, the signal is undetectable.
Result: power under carryover remained low, comparable to
Architecture B. Lesson: tracking drug exposure decay in the
moderation signal reintroduces the very phenomenon Architecture A is
intended to avoid.

**Attempt 3: Additive with time-on-drug and uncentered biomarker
(incorrect scaling).**

```
BR_{it}^{mod} = BR_{it}^{base} + B_i * tod_t * beta_bm
where B_i = raw (uncentered) biomarker
```

This was the form in the pre-merge 2025 tidyverse codebase (commit
957c8bd, `pm_functions.R` lines approximately 650-670). Problem: the
biomarker mean is approximately 124 mmHg (systolic BP) with SD
approximately 15; `tod` reaches 10 to 20 weeks. Moderation shift at
the last on-drug timepoint is `124 * 10 * 0.45 = 558` BR units
against a BR noise SD of 8: a signal-to-noise ratio of about 70 to 1.
Result: power = 1.00 everywhere. Lesson: dimensional analysis would
have caught this immediately. The product `B_i * tod_t` has units
of mmHg-weeks, which is not commensurate with BR (unitless symptom-
score change).

**Attempt 4: Centered, on-drug-binary, calibrated (correct).**

```
BR_{it}^{mod} = BR_{it}^{base}
              + c_bm * (B_i - mu_B) / sigma_B * 1[on drug at t] * sigma_BR
```

Derivation: this matches Architecture B's conditional expectation
from the MVN model,

$$E[BR \mid \mathrm{bm}, \text{on drug}]
  = \mu_{BR} + c_{bm} \cdot \frac{\sigma_{BR}}{\sigma_B}
  \cdot (\mathrm{bm} - \mu_B)$$

Why it works:

1. **Centering**: the standardized biomarker `B_i^*` has mean zero.
   The population-average moderation is zero. Only deviations
   produce differential BR. Correct for an interaction, not a main
   effect.

2. **Binary indicator**: `1[on drug]` is 1 during treatment, 0
   during off-drug. Unlike `tod` (grows with exposure) or `phi(t)`
   (decays), it produces a clean constant signal during all
   on-drug timepoints and zero during off-drug. Carryover affects
   the BR population means but not the moderation term.

3. **Scaling by sigma_BR**: a one-SD biomarker shift produces a
   `c_bm * sigma_BR` shift in BR. For `c_bm = 0.45` and
   `sigma_BR = 8` this is 3.6 BR units, calibrated to match the
   Architecture B conditional expectation.

### 8.2 Empirical Verification

Drawing 500 participants with identical parameters
(`c.bm = 0.45`, crossover design, no carryover):

| | Architecture A | Architecture B |
|---|---|---|
| BR shift per SD of biomarker | 3.36 | 3.48 |
| Residual SD | 5.23 | 3.89 |

The shift magnitudes are comparable. Architecture A has slightly
higher residual variance because the BM-BR correlation in the
covariance matrix is absent (the noise components are independent),
whereas Architecture B's MVN draw produces correlated noise that
reduces the residual variance in the regression.

### 8.3 Hendrickson Publication Parameters Under Architecture A

With the original publication grid (`c.bm` up to 0.6, `t_half` in
{0, 0.1, 0.2}), Architecture A makes previously infeasible
combinations evaluable. Power at `c.bm = 0.6`, 50 replicates:

| Design | N | t_half = 0 | t_half = 0.1 | t_half = 0.2 |
|---|---|---|---|---|
| N-of-1 | 70 | 0.98 | 0.96 | 0.96 |
| OL+BDC | 70 | 0.92 | 0.92 | 0.92 |
| CO | 70 | 0.66 | 0.62 | 0.70 |
| OL | 70 | 0.14 | 0.16 | 0.20 |
| N-of-1 | 35 | 0.68 | 0.78 | 0.70 |
| OL+BDC | 35 | 0.58 | 0.66 | 0.68 |
| CO | 35 | 0.40 | 0.28 | 0.26 |

Power at `c.bm = 0.3`, 50 replicates:

| Design | N | t_half = 0 |
|---|---|---|
| N-of-1 | 70 | 0.48 |
| OL+BDC | 70 | 0.36 |
| CO | 70 | 0.26 |

Type I error at `c.bm = 0` is 0.00 to 0.08 across all designs
(nominal 5 percent). Observations:

1. Carryover resilience is confirmed: power is stable across all
   three carryover half-lives for every design under Architecture A.
2. The design ranking is preserved:
   `N-of-1 > OL+BDC > CO > OL`.
3. `c.bm = 0.6` is now evaluable. N-of-1 at `N = 70` achieves 0.96 to
   0.98 power, previously unobtainable because the parameter failed
   PD validation under Architecture B.

---

## 9. Clinical Application and Trial Analysis

Source document: 14-trial-analysis-implementation-guide.tex.

### 9.1 Collaborator Deliverables

Three research goals mapped to implemented functions:

**Goal 1: Did prazosin improve PTSD symptoms?**

```r
result <- analyze_trial_extended(td, dat, op, threshold = 5)
print_trial_summary(result)
```

Defined in `R/carryover_analysis.R`. Outputs:

- Main drug effect (`Dbc` coefficient).
- Biomarker-treatment interaction (`bm:Dbc`).
- 95 percent confidence intervals.
- Variance components (random intercept SD, residual SD, AR(1) phi).
- Predicted response curves per participant via `predict(fit, newdata)`.
- Clinical decision threshold:
  `bm_threshold = mean_bm + (threshold - beta_Dbc) / beta_interaction`.

Status: complete.

**Goal 2a: Characterize carryover in actual data.**

```r
cc <- characterize_carryover(td, dat,
                             half_lives = c(0.25, 0.5, 1.0, 2.0, 4.0))
print_carryover_summary(cc)
# cc$best_t_half, cc$best_model
```

Sweeps candidate half-lives, fits `lme_analysis()` at each,
compares AIC/BIC/log-likelihood, returns the best-fitting half-life
and full model object. Validated to recover true `t_half = 1.0`
from simulated data.

Status: complete.

**Goal 2b: Identify best sequence of response measures.**

Status: not yet implemented; conceptual only.

**Goal 3a: Simulation power calculations.**

Status: complete. Validated 1,000-replicate power estimates with
Type I error 2.5 to 5.8 percent across designs. Zero PD failures
across 162 parameter combinations. Runtime 8 minutes (200 reps),
47 minutes (1,000 reps).

**Goal 3b: Actual analysis output structure.**

`lme_analysis(..., full_model_out = TRUE)` returns:

- `form`: the model formula.
- `fit`: the full nlme::lme model object.
- `data`: analysis-ready data with derived variables (`Dbc`, `tsd`,
  `De`).
- `summary`: beta, SE, p-value.

From the full model object, downstream analyses access
`summary(fit)$tTable`, `fitted(fit)`, `residuals(fit)`, `ranef(fit)`,
`intervals(fit)`, and `VarCorr(fit)`.

### 9.2 Prazosin-PTSD Context

From docs 02 and 20:

- Resting systolic blood pressure is a proxy for central noradrenergic
  tone, a PTSD subtype preferentially responsive to prazosin (an
  alpha-adrenergic antagonist).
- The biomarker is a statistical predictor of response class, not a
  causal pharmacokinetic mediator. Architecture B is the appropriate
  DGP.
- Clinical trial history: Raskind et al. 2003, 2007, 2013 (smaller
  trials, favorable results on nightmares and PTSD symptoms);
  Raskind et al. 2018 (PACT trial, N = 304, null result on primary
  endpoint); Raskind et al. 2023 (posttraumatic headache pilot).
- The PACT null result motivates the search for a predictive
  biomarker and is the central clinical puzzle the pmsimstats
  framework addresses.

### 9.3 Implication for Trial Design

From doc 02, Section 5 (precision medicine trial design): the
architectural choice matters most in designs that deliberately
expose each participant to both treatment states to estimate
within-subject responsiveness (crossover, basket, platform, N-of-1
hybrid). Parallel-group enrichment designs escape Architecture B's
carryover sensitivity because they have no off-drug phase.
Recommendation: when the biological mechanism is uncertain, run the
power calculation under both architectures and treat the
Architecture B estimate as the conservative anchor for sample size
planning. This is consistent with the PATH Statement (Kent et al.
2020).

---

## 10. Alternative Analysis Strategies Under Architecture B

Source document: 02-dgp-mean-moderation-vs-mvn.tex, Section 6.

Six candidate strategies for recovering power under Architecture B,
evaluated comparatively.

**6.1 Restrict analysis to on-drug observations only.**

Eliminates correlation erosion (all retained observations have full
`c_bm`). Costs: reduced effective sample size and loss of the
off-drug reference, which confounds biological, placebo, and
time-varying components. Best for designs with many on-drug
timepoints (OL has 8, Hybrid has 5 to 6). Least useful for CO, which
would lose half its observations.

**6.2 Weighted analysis.**

On-drug observations weight 1. Off-drug observations weight
`w_t = exp(-lambda * t_sd)`. Rationale: information-theoretic.
Observations with higher drug exposure carry more information about
the biomarker-treatment interaction because the correlation signal
is stronger.

**6.3 Within-subject contrast.**

For designs with both on-drug and off-drug periods (CO, Hybrid,
OL+BDC):

$$\bar{Y}_{i,\mathrm{on}} - \bar{Y}_{i,\mathrm{off}}
  = \alpha + \gamma \cdot B_i + \epsilon_i$$

Sidesteps `Dbc` parameterization and the repeated-measures model.
Within-subject averaging reduces noise.

**6.4 Two-stage random slopes.**

Stage 1: fit a random-slopes LME.

$$Y_{it} = \beta_0 + \beta_1 t + \beta_2 D_{it}
  + u_{0i} + u_{1i} D_{it} + \epsilon_{it}$$

Stage 2: regress the estimated random deviations on the biomarker.

$$\hat{u}_{1i} = \delta_0 + \delta_1 B_i + \eta_i$$

Separates within-subject treatment effect estimation (affected by
carryover) from between-subject biomarker association (which uses
cross-sectional relationships).

**6.5 Exclusion of contaminated observations.**

Exclude the first one or two off-drug timepoints, which carry the
most carryover contamination. Later off-drug timepoints provide a
cleaner contrast with the washed-out state.

**6.6 Design-level solutions.**

- Longer washout periods between drug phases.
- More on-drug timepoints relative to off-drug.
- Front-loaded drug exposure (OL+BDC design): a long initial on-drug
  block establishes the BM-BR correlation before off-drug
  observations begin.

**6.7 Comparative evaluation.**

| Strategy | Discards data | Complexity | Expected benefit |
|---|---|---|---|
| On-drug only | Yes | Low | High for OL/Hybrid; poor for CO |
| Weighted | No | Low | Moderate; depends on weight calibration |
| Within-subject contrast | Aggregates | Low | Moderate; robust to misspecification |
| Two-stage random slopes | No | Moderate | Unknown; separates estimation stages |
| Exclude early off-drug | Some | Low | Moderate; improves contrast clarity |
| Design modification | N/A | N/A | High; addresses root cause |

---

## 11. Latent-Class and Psychometric Connections

Source document: 02-revision-response-latent-class-expanded.md.

### 11.1 The Faithful Generative Model

If the biomarker indexes an unobserved responder/non-responder
dichotomy, the mechanistically correct DGP is a finite mixture:

$$Z_i \mid B_i \sim \mathrm{Bernoulli}(\pi(B_i)),\quad
  BR_{it} \mid Z_i = z \sim f_z(t, D_{it})$$

with `Z_i` the latent class, `pi(B_i) = Pr(Z_i = 1 | B_i)` typically
a logistic function of the biomarker, and class-specific response
distributions `f_0, f_1` differing in drug-response parameters.

Key property: carryover under this DGP operates only on the
class-specific means and does not erode any correlation signal,
because class membership `Z_i` is drug-state-invariant.

### 11.2 MVN Differential Correlation as Second-Moment Approximation

A two-component mixture with class-dependent means and common
within-class covariance produces a marginal joint distribution of
`(B_i, BR_{it})` whose first two moments match, to leading order,
those of a single MVN with `Cor(B, BR) = c_bm` chosen so that
`c_bm^2` equals the between-class variance fraction. Under mild
regularity conditions (class overlap, moderate `pi`), the mixture's
second-moment structure is well-approximated by MVN with
differential correlation in treatment states where class-specific
means separate. Architecture B as implemented captures the
covariance implications of the latent-class mechanism without
committing to the discrete generative form.

### 11.3 Where the Approximation Breaks Down

Three regimes:

1. **Bimodal responses**: when `f_0` and `f_1` are well-separated
   (non-responders show near-zero drug effect, responders show a
   large effect), the marginal BR distribution is bimodal and no
   MVN matches its tails. Power under the true mixture DGP is
   **higher** than Architecture B predicts because mixture analysis
   can exploit bimodality.

2. **Strong class-membership gradient**: when `pi(B_i)` is steep
   (approaching a step function at a threshold), the biomarker
   behaves almost deterministically as a class label and
   Architecture A becomes the better approximation.

3. **Class-varying covariance**: if autocorrelation, residual
   variance, or placebo-response parameters themselves differ
   between classes, no single-component MVN can represent the
   generative structure. Analysis models assuming constant
   covariance will be misspecified.

### 11.4 Taxonomy of Finite Mixture Models

From the psychometric tradition:

- **Latent class and latent profile analysis** (Lazarsfeld and
  Henry 1968; Goodman 1974; Clogg 1995; McLachlan and Peel 2000):
  small number of discrete subpopulations. LCA for categorical
  indicators; LPA for continuous. R: `poLCA`.
- **Growth mixture models and group-based trajectory models**
  (Muthen and Shedden 1999; Nagin 1999, 2005): class-specific
  growth curves with within-class random-effect variation. Closest
  psychometric analogue of the mixture DGP discussed above.
  Identifiability hazards (Bauer and Curran 2003; Bauer 2007)
  apply with full force at N-of-1 sample sizes.
- **Factor mixture models** (Yung 1997; Arminger, Stein, and
  Wittenberg 1999; Lubke and Muthen 2005): continuous latent
  factors combined with discrete latent classes. Subsumes both
  pure continuous (single-class factor model) and pure discrete
  (LCA with no within-class structure) as limits.
- **Regression mixtures and mixtures of experts** (DeSarbo and
  Cron 1988; Wedel and DeSarbo 1995; Jacobs et al. 1991): class-
  specific regression coefficients with covariate-dependent class
  membership probabilities.

### 11.5 The Architecture B Spectrum

Architecture B as implemented is the **most restrictive** operational
case of a broader family that ranges from mean-only moderation
(Architecture A) through combined mean-plus-covariance moderation
(Arminger et al. 1999; the canonical Psychometrika treatment) to
heterogeneous-random-effects models (Verbeke and Lesaffre 1996;
Proust-Lima et al. 2017).

Implication: the substantial power loss under carryover (40 to 60
percent) documented in Section 1.4 is a specific property of the
restrictive covariance-only case and may not generalize to
mean-plus-covariance or heterogeneous-random-effects variants.

Recommendation for a complete sensitivity analysis: report power
under at least three points on the spectrum: Architecture A
(mean-only), Architecture B as implemented (covariance-only), and a
combined mean-plus-covariance specification of the Arminger form.

### 11.6 Practical Constraint and Next Extension

Identifiability concerns for mixture-based analysis at N-of-1 trial
sample sizes (typically 30 to 100) are acute. Published FMM
applications involve hundreds to thousands of participants. The
open question is not whether mixture modeling is more faithful (it
is), but whether the fidelity gain is recoverable at realistic
sample sizes. Analysis-strategy mitigation (Section 10 above) is
more likely than mixture-model analysis to deliver recoverable
power gains at realistic trial sizes.

If the pmsimstats program pursues a richer DGP, the natural target
is the **heterogeneous-random-slopes form** rather than a full FMM.
This preserves the analysis model's compatibility with existing
N-of-1 estimation machinery (`nlme::lme` with participant-specific
random treatment effects), admits a biologically faithful latent-
class generative structure at the DGP level, preserves the audit
chain of the current framework, and permits direct comparison with
the covariance-only baseline under carryover.

R packages relevant in approximate order of fit: `lcmm` (closest
tooling for the heterogeneous-random-slopes form), `flexmix`
(regression-mixture framework), `mclust`, `gamlss`, `poLCA`,
`OpenMx` and `lavaan` with mixture extensions. Reference
implementation outside R: Mplus.

---

## 12. The 2x2 Factorial: DGP Carryover Versus Analysis-Model Carryover

Source document: 10-carryover-analysis-model-assessment.tex.

### 12.1 The Question and the Framework

Standard methodology (Senn 2002; Jones and Kenward 2014) for
assessing carryover in the analysis model is a 2x2 factorial with
DGP carryover on one axis and analysis-model carryover on the
other, reporting Type I error, power, bias, and confidence interval
coverage in each cell. Senn's modern consensus:

- The two-stage procedure (test for carryover first, then condition)
  is invalid; carryover tests have low power.
- If carryover is scientifically plausible, either design the trial
  to eliminate it via adequate washout, or always include carryover
  in the model and accept the associated power loss. Do not
  test-then-decide.
- For simulation studies: generate data under a range of carryover
  magnitudes (including zero), analyze under models that include
  and exclude the carryover term, and report operating
  characteristics for each combination.

### 12.2 Three Analysis Approaches

**Approach A: Continuous drug indicator (Dbc).**

Formula: `Sx ~ bm + t + Dbc + bm:Dbc`. `Dbc` equals 1 on drug and
`(1/2)^(s * t_sd / t_half)` off drug. Conflates treatment and
carryover; captures both treatment-biomarker and carryover-
biomarker interactions simultaneously.

**Approach B: Separate carryover covariate (tsd).**

Formula: `Sx ~ bm + t + Db + bm:Db + tsd`. The linear time-since-
discontinuation term absorbs carryover. Closer to standard
crossover methodology but uses a linear carryover term, which is
misspecified relative to the exponential DGP.

**Approach C: Ignore carryover.**

Formula: `Sx ~ bm + t + Db + bm:Db`. Off-drug observations treated
as if no drug were present. Biases the treatment contrast toward
zero, reducing power, and may inflate Type I error if the bias is
systematic.

### 12.3 Proposed Simulation Study

Design:
- DGP factor: `t_half` in `{0, 0.25, 0.5, 1.0, 2.0}` weeks.
- Analysis factor: three models (Ignore, Separate, Continuous).
- Fixed: N = 70, rho = 0.7, all four trial designs.
- Biomarker correlation: `c_bm` in `{0, 0.3, 0.45}`.
- 1,000 replicates per cell.
- Total cells: 5 * 3 * 3 * 4 = 180.

Metrics: Type I error (at `c_bm = 0`), power, bias, RMSE, coverage.

Expected findings:
- Approach 1 (Ignore): correct Type I error at `t_half = 0`, may
  inflate at higher carryover; power declines with increasing
  carryover.
- Approach 2 (Separate): correct Type I error across all carryover
  values; power slightly lower at `t_half = 0` (cost of extra
  parameter), higher at `t_half > 0`; mild bias from linear versus
  exponential misspecification.
- Approach 3 (Continuous): correct Type I error if decay rate
  matches DGP; power depends on whether the conflation helps or
  hurts.

Infrastructure: models supported via existing codebase through
`op$carryover_t1half` and `op$simplecarryover` parameters.

Status: proposed simulation study; not yet complete.

---

## 13. Simplification Alternatives (Not Adopted)

Source document: 17-simplification-plan.tex.

Doc 17 proposed reducing complexity of the pmsimstats framework
substantially. The 2026-04-08 merge took a different path (adding
Architecture A rather than simplifying). The document's novel
contribution is the most thorough formal specification of six
analysis methods available in the project, including precision
equations and design implications. It is preserved here because the
analytical framework remains useful for cross-method comparisons.

### 13.1 Six Interaction Test Methods

**1. ANCOVA on phase-mean change scores (recommended in doc 17).**

Collapse to phase means: `Y_bar_on = (1/n_on) sum(Y_it | D=1)` and
similarly for off. Compute `Delta_i = Y_bar_on - Y_bar_off`. Regress
`Delta_i = gamma_0 + gamma_1 * bm_i + e_i`. Test `H_0: gamma_1 = 0`
via t-test on `gamma_1_hat` with df = N - 2. Precision:
`Var(gamma_1_hat) = sigma_Delta^2 / (N * Var(bm))`.

**2. Paired t-test on change scores.**

Tests `H_0: E[Delta_i] = 0` (main effect, not interaction).
Requires dichotomization of the biomarker (e.g., median split) for
the interaction. Efficiency loss from dichotomization approximately
36 percent (`2/pi` approximately 0.64 relative efficiency).

**3. Generalized Estimating Equations (GEE).**

Population-averaged model with sandwich estimator. Key precision
parameter: `Var_within(D_i) = D_bar_i (1 - D_bar_i)`, maximized at
`D_bar = 0.5` (equal on/off time). For the OL design,
`Var_within(D) = 0`, so the interaction is not estimable.

**4. Random-intercept LMM.**

Subject-specific model. Precision depends on
`SS_within(D_i) = sum_t (D_it - D_bar_i)^2`. Design-specific
values: OL has 0, OL+BDC approximately 1.5, CO 2.0, N-of-1
approximately 2.9. Short blocks (high K, low L) suffer less from
autocorrelation than long blocks because fewer consecutive
observations are highly correlated.

**5. Simple linear regression (ignoring clustering).**

Downward-biased standard errors; Type I error inflation with
ICC = 0.3 and T = 8 is approximately 15 to 25 percent at nominal
alpha = 0.05. Not recommended for inference.

**6. Repeated-measures ANOVA (biomarker dichotomized).**

Requires dichotomization (36 percent information loss). Sphericity
assumption violation requires Greenhouse-Geisser or Huynh-Feldt
correction, further reducing power. RM-ANOVA retains approximately
0.64 times ANCOVA efficiency. Design-specific behavior: OL uses
time as the within-subject factor; CO is a clean two-level factor
with balanced phases; N-of-1 with multiple blocks raises
sphericity concerns.

### 13.2 Cross-Design Summary (doc 17, Table 8.1)

| | OL | OL+BDC | CO | N-of-1 |
|---|---|---|---|---|
| Interaction term | bm times t | bm times D | bm times D | bm times D |
| Target | Total improvement rate | Drug-specific | Drug-specific | Drug-specific |
| Confounded with | ER, TV slopes | -- | -- | -- |
| LMM info | sum bm^2 SS(t) | sum bm^2 (1.5) | sum bm^2 (2.0) | sum bm^2 (2.9) |
| Transitions | 0 | 1 | 1 | 3 |
| Autocorr vulnerability | Low | Medium | Medium | Low per block |
| Carryover vulnerability | None | Low | Low | High |

### 13.3 The OL Design Is a Fundamentally Different Test

For the open-label design, all six methods face the same problem:
there is no within-subject drug variation. The interaction must be
formulated as biomarker-by-time. This is a weaker test than the
biomarker-by-drug interaction for two reasons: confounding (the
slope `dY/dt` includes BR, ER, and TV; the biomarker-by-time
interaction conflates all three) and precision (the within-subject
information for the time slope depends on
`SS_within(t) = sum(t - t_bar)^2`, which is typically smaller than
`SS_within(D)`).

### 13.4 Adoption Status

Not adopted. The 2026-04-08 merge preserved the LME with corCAR1
analysis model and added Architecture A rather than simplifying to
ANCOVA on phase-mean change scores. The doc is retained because
the formal specification of analysis method precision across
designs is the most detailed in the project and remains useful for
power planning in related contexts.

---

## 14. Code Audit Findings (2026-04-16)

Source document: audit-2026-04-16.md.

The audit identified nine issues across three severity tiers.

### 14.1 P0 (Execution or Correctness Blockers)

**Finding T1**: Unqualified dplyr calls in the tidyverse Architecture
A path.

Location: `implementations/tidyverse/R/functions.R:421-423`.

```r
bm_mean <- baseline_param |> filter(cat == 'bm') |> pull(m)
bm_sd   <- baseline_param |> filter(cat == 'bm') |> pull(sd)
br_sd   <- resp_param     |> filter(cat == 'br') |> pull(sd)
```

These lines execute only when `dgp_architecture == 'mean_moderation'`.
If called without `library(dplyr)` in scope, they will either fail
at runtime or resolve `filter()` to `stats::filter`. The parity test
masks the issue because it runs in a tidyverse-loaded session. Fix:
qualify both calls with `dplyr::` and add a regression test
invoking `generate_data(..., dgp_architecture = 'mean_moderation')`
with only declared Imports loaded.

**Finding A1**: Broken source path.

Location: `analysis/carryover_factorial/02-run-factorial-2025.R:20`.

Sources `analysis/2025/01-pm-functions.R`, which has been moved to
an external archive. The script fails at load time. Fix: change the
source target to `implementations/tidyverse/R/functions.R` and
update downstream function-name references, or restore a copy from
the external archive if historical comparison is needed.

### 14.2 P1 (Publication and CRAN Blockers)

**Finding C1**: Silent positive-definiteness correction.

Locations: `generateData.R:284-299` (root), `orig/generateData.R:242-253`,
`ext/generateData.R:284-299`, `tid/functions.R:337-356`.

When `is.positive.definite(sigma)` returns FALSE, all
implementations call `make.positive.definite(sigma, tol = 1e-3)`
without recording the correction. Doc 04 reports 24.7 percent of
matrices in the publication required correction. The current audit
cannot verify the present rate without instrumentation. Fix: return
a `pd_corrected` flag from `buildSigma()` and aggregate into
per-parameter-set summaries.

**Finding C2**: Hard-coded `issingular = FALSE`.

Locations: `lme_analysis.R:178, 221` (copies in orig, ext);
`tid/functions.R:654`.

The `issingular` column is always FALSE even when the underlying
fit is singular. nlme exposes `fit$apVar` and
`fit$modelStruct`; singularity is inferable. Fix: detect
singularity from `fit$apVar` (character string 'Non-positive-
definite approximate variance-covariance') and check random-effects
variance for zero eigenvalues; populate `issingular` accordingly.

**Finding C5**: Unvalidated `carryover_t1half`.

Location: `generateData.R:172-177` (copies).

A negative `carryover_t1half` produces a negative `lambda_cor`,
turning `exp(-lambda_cor * tsd)` into an exponentially growing
correlation. No validation. Fix: add `stopifnot(carryover_t1half >= 0)`
at the top of `generateData()`.

**Finding C6**: Stale imports.

The root DESCRIPTION Imports list includes `lme4`, `lmerTest`,
`ggpubr`, and `gridExtra`, none of which are called anywhere in
`R/`, `vignettes/`, or the three package-style implementations.
These are residues from the pre-nlme analysis model and earlier
plotting code. Fix: remove them from DESCRIPTION; if any surfaces in
a vignette, promote to Suggests. The audit notes these will generate
NOTEs on CRAN.

**Finding C7**: Missing test coverage.

No test exercises `dgp_architecture = 'mean_moderation'`, null-effect
Type I error at nominal alpha, `corCAR1` correctness under known
`phi` ground truth, or the silent PD correction path. Fix: add four
`tinytest` cases (or migrate to `testthat` 3e) covering these gaps.

**Finding A2**: Missing READMEs in `analysis/` subfolders.

Locations: `analysis/02-carryover-sensitivity/`,
`analysis/architecture_comparison/`, `analysis/carryover_factorial/`,
`analysis/figure5/`. Each contains numbered driver scripts without a
README. Fix: add a one-paragraph README per subfolder describing
research question, inputs, outputs, and runtime.

### 14.3 P2 (Hygiene and Maintainability)

**Finding C3**: `eval(parse(text = ...))` for column selection.

Locations: `censordata.R:50, 84`; `plottingfunctions.R:136`;
`lme_analysis.R:93` (partial). Fix: replace with `get()`,
`dat[[name]]`, or `dplyr::all_of(name)`. No user-facing behavior
change.

**Finding C4**: Positional coefficient extraction.

Locations: `lme_analysis.R:193-217` (copies); `tid/functions.R:637-650`;
`carryover_analysis.R:71-93`. Relies on `intersect(c('bm:Dbc', 'Dbc:bm'),
rownames(...))` and `target[1]`. `carryover_analysis` uses a more
flexible pattern. Fix: harmonize around the `carryover_analysis`
pattern.

### 14.4 Parity Status

From the audit: `parity_baseline.rds` and the recent
`parity_report.rds` differ by two bytes, indicating stable parity
between `original-extended` and `tidyverse` at tolerances of 1e-10
(sigma, data) and 1e-6 (beta, SE, p-value). No expected divergences
are recorded in `parity_divergences.csv` (zero data rows).

### 14.5 Publication Readiness Summary

- CRAN blockers: T1 will fail `R CMD check --as-cran`; C6 will
  generate NOTEs; C7 weakens the submission narrative.
- Manuscript blockers: A1 will break reviewer reproduction; C1
  (silent PD correction) should surface in methods or appendix but
  the current rate is not instrumented; C2 (singularity flag always
  FALSE) may affect reported convergence statistics.

With P0 and P1 items addressed, the package and manuscript are
well positioned for submission.

---

## 15. Parity Development Workflow

Source document: 22-tidyverse-development-parity-workflow.md.

### 15.1 Baseline and Inner Loop

Before starting a tidyverse refactor:

```bash
git checkout -b tidyverse-refactor
Rscript implementations/test-parity-extended-tidyverse.R
cp implementations/parity_report.rds implementations/parity_baseline.rds
```

Commit `parity_baseline.rds` with the branch. At the end of
development, diff the current `parity_report.rds` against it to
confirm only intended differences.

Two inner-loop cadences:

- `--quick` flag: approximately 2 seconds, 16 cells. Use for
  iteration on a specific function.
- Full sweep: approximately 10 seconds. Run before every commit.

An optional pre-commit hook blocks commits on failure:

```bash
Rscript implementations/test-parity-extended-tidyverse.R --quick || \
  { echo 'parity broken; run --full for details'; exit 1; }
```

### 15.2 Edit Classification

- **Pure refactors** (rename, split function, swap dplyr for purrr,
  extract helper): parity must hold. If it breaks, a bug was
  introduced.
- **Feature additions with default values**: parity must still
  hold for default call paths.
- **Behavior changes** (different DGP parameterization, different
  MVN ordering, different sigma caching): parity will break by
  design; decide whether to apply to both implementations in
  lockstep.

Most development falls into the first two categories. The third
should be treated as a separate project with explicit scope.

### 15.3 Breakage Diagnosis

When parity breaks mid-refactor, the test prints the first failing
cell. Diagnose in this order:

1. **Sigma failure**: the change likely touched
   `build_sigma_matrix` or `build_correlation_matrix`. The BM-BR
   correlation branch is the most common culprit.
2. **Data failure**: touched MVN sampling, post-draw mean
   moderation, or participant/path assembly.
3. **Analysis failure**: touched formula construction, `Dbc`
   computation, or `nlme::lme` invocation.

The `--quick` run covers 16 cells in about 2 seconds. Re-running
after every approximately 20 lines of changes localizes bugs to
the edit that introduced them.

### 15.4 RNG Guardrails

The parity test calls `set.seed()` at the per-cell level. If a
refactor consumes random draws in a different order (for example
swapping `MASS::mvrnorm` for a Cholesky-based draw), parity will
break even though the code is semantically correct. Defenses:

- Do not change RNG-consuming loops unless `original-extended`
  is updated in the same commit.
- Keep `n_cores = 1` for parity-relevant comparison. Introduce
  parallelism only after functional parity is locked in.

### 15.5 Intentional Divergence

When a refactor intends to change tidyverse behavior in a way that
`original-extended` should not follow (for example fixing a bug
that lives only in data.table code, or using tidyverse as a
testbed for a new feature):

- Note in the commit message or `DIVERGENCES.md` which cells are
  allowed to fail and why.
- Mark those cells in the script with an `expected_divergent` flag.
- Re-baseline after the intentional divergence lands:
  `cp parity_report.rds parity_baseline.rds`.

The trap to avoid is letting an intentional divergence mask a
subsequent unintentional one.

### 15.6 Pre-Merge Checklist

- [ ] cp parity_report.rds parity_baseline.rds at start of branch.
- [ ] --quick after every focused edit (approximately 2 s).
- [ ] Full sweep before every commit (approximately 10 s).
- [ ] If parity breaks, read the first failing cell; fix before
  proceeding.
- [ ] Never change RNG consumption order in a pure refactor.
- [ ] Add edge cells when new features are added.
- [ ] Document any intentional divergence explicitly.

One-liner for diff before merging:

```r
a <- readRDS('implementations/parity_baseline.rds')
b <- readRDS('implementations/parity_report.rds')
all.equal(a$grid, b$grid)  # TRUE, or character vector of diffs
```

---

## 16. Publication Strategy

Source document: 22-journal-target-recommendation.md.

### 16.1 Manuscript Summary

Manuscript: 'Two Architectures for Simulating Biomarker-Treatment
Interactions: Implications for Statistical Power Under Carryover'.
The manuscript formalizes the distinction between direct mean
moderation (Architecture A) and differential correlation via MVN
(Architecture B), with the central finding that architecture choice
determines carryover impact on power: 10 to 15 percent loss under
Architecture A, 40 to 60 percent under Architecture B.

Methodological positioning: intersection of statistical simulation
methodology, N-of-1 and crossover trial design, and biomarker-
treatment interaction modeling.

### 16.2 Ranked Journal Targets

1. **Statistics in Medicine** (primary). Premier journal for
   statistical methodology in clinical trial design. Primary
   audience (biostatisticians designing trials) would act on the
   finding. Collection precedent: Senn (2016).

2. **Biometrics**. Yap, Wang, and Harhay (2024) is the closest
   comparator (carryover impact on power in crossover designs).
   The manuscript extends this work by showing that the
   mechanism determining the interaction affects the magnitude of
   power loss.

3. **Trials**. Highest concentration of N-of-1 methodology papers
   in the related literature (4 of 50 in the annotated
   bibliography). Open access, practical orientation.
   Collection precedents: Haller and Ulm (2018), Chen et al.
   (2024), Senn and Julious (2024).

4. **Statistical Methods in Medical Research**. Publishes
   simulation-based methodological work for trial design.
   Collection precedent: Arends et al. (2019).

5. **Contemporary Clinical Trials Communications**. Practical
   framing and clinical application fit well. Collection
   precedents: Sturdevant and Lumley (2021), Liu et al. (2021).

6. **Harvard Data Science Review**. Collection precedent:
   Schork (2022). More of a perspective venue; heavy simulation
   results may not fit format.

7. **Healthcare (MDPI)**. Collection precedents: Blackston et al.
   (2019), Wang and Schork (2019). Fast review; lower stigma
   risk if strengthened before submission.

Primary recommendation: target **Statistics in Medicine**. Backup:
target **Trials** if Statistics in Medicine considers the scope
too specialized.

---

## 17. Literature Synthesis

Source document: 20-annotated-bibliography.tex.

Fifty references in the annotated bibliography organized by topic.
The summaries below focus on the claims most directly relevant to
the pmsimstats framework.

### 17.1 Core N-of-1 Simulation Methodology

- **Thomas et al. (2020)**, Frontiers in Digital Health 2:13.
  Derives analytic power formulae for aggregated N-of-1 trials
  designed to detect predictive biomarker effects, validates via
  Monte Carlo. Finding: aggregated N-of-1 achieves adequate power
  with fewer participants than parallel-group designs. The
  framework pmsimstats-ng revises.
- **Blackston et al. (2019)**, Healthcare 7(4):137. 5,000
  simulated trials per scenario across four variance structures
  comparing aggregated N-of-1, parallel RCT, crossover.
  Finding: aggregated N-of-1 outperforms traditional designs but
  has elevated Type I error under moderate-to-strong unmodeled
  carryover.
- **Wang and Schork (2019)**, Healthcare 7(3):84. Analytical
  power formulae under AR(1) serial correlation. Finding: serial
  correlation substantially reduces power; increasing evaluation
  periods partially compensates.

### 17.2 Carryover and Temporal Correlation

- **Yap, Wang, and Harhay (2024)**, Biometrics 80(2):ujae023.
  Formal framework for behavioral carryover in crossover trials
  (psychological mechanisms, not pharmacokinetic). Closest
  comparator for the pmsimstats manuscript. Finding: behavioral
  carryover can bias treatment effect estimates even when
  pharmacokinetic carryover is negligible. Directly relevant to
  the prazosin-PTSD application.
- **Sturdevant and Lumley (2021)**, Contemporary Clinical Trials
  Communications 22:100711. Mixed-effects model approach for
  testing carryover. Finding: 'carryover has been largely ignored
  except as a nuisance parameter'.
- **Schork (2022)**, Harvard Data Science Review SI3. Serial
  correlation in personalized and aggregated studies.
  Finding: ignoring serial correlation inflates false positives.
- **Tang and Landes (2020)**, PLoS One 15(2):e0228077. Modified
  t-tests accounting for serial correlation in N-of-1 data.
  Finding: proposed tests maintain nominal Type I where standard
  t-tests fail.

### 17.3 Biomarker-Treatment Interactions

- **Haller and Ulm (2018)**, Trials 19:128. Simulation study of
  biomarker-treatment interaction estimation in randomized trials
  with time-to-event outcomes. Same genre as the pmsimstats
  manuscript (Architecture A paradigm).
- **Grenet et al. (2020)**, Basic and Clinical Pharmacology and
  Toxicology 126:59-64. Simulation comparing crossover versus
  parallel-group for binary predictive biomarker identification.
  Finding: crossover more efficient when within-subject
  correlation is high (supports N-of-1 use in pmsimstats).

### 17.4 Analysis Methods

- **Senn, Julious, and Araujo (2014)**, PLoS One 9(2):e87752.
  Four methods for N-of-1 analysis (fixed-effects, random-
  effects, Bayesian hierarchical, paired). Finding: all similar
  for population inference; differ substantially for individual
  treatment effect estimation.
- **Peterlin, Kejzar, and Blagus (2023)**, TEST 32:184-210.
  Diagnostic tests for correct specification of fixed- and
  random-effects design matrices in LMM.

### 17.5 Sample Size

- **Senn (2019)**, Statistical Methods in Medical Research
  28(2):372-383. Five distinct criteria for planning sample size
  in N-of-1 trials. Finding: criterion depends on goal
  (individual-level, population-level, heterogeneity detection).
- **Schmid, Duan, and Kravitz (2014)**. AHRQ technical report.
  Comprehensive methods for N-of-1 trial design.
- **Senn, Julious, and Chen (2024)**, Trials 25:264. Methodological
  review of published randomized N-of-1 trials. Finding:
  widespread inconsistencies in reporting and analytic practice.

### 17.6 Bayesian Approaches

- **Schmid and Staudenmayer (2021)**, Harvard Data Science Review.
  Review of Bayesian modeling for N-of-1 trials.
- **Daza, Wac, and Oppezzo (2018)**, AMIA Annual Symposium
  Proceedings 2018:532-541. Counterfactual framework for causal
  inference in N-of-1 trials via self-tracked time series.

### 17.7 Prazosin-PTSD Clinical Context

- **Raskind et al. (2003)**, American Journal of Psychiatry
  160(2):371-373. Initial proof-of-concept (N=10), double-blind
  crossover; significant reduction in nightmares.
- **Raskind et al. (2013)**, American Journal of Psychiatry
  170(9):1003-1010. Parallel-group RCT (N=67) in active-duty
  soldiers; significant effects; observed effect sizes and
  heterogeneity calibrate simulations.
- **Raskind et al. (2018)**, New England Journal of Medicine
  378(6):507-517. PACT trial (N=304); no significant difference
  versus placebo on primary outcome. The null result motivates
  the pmsimstats framework search for a predictive biomarker.
- **Raskind et al. (2023)**, Headache 63(6):751-762. Pilot for
  posttraumatic headache prophylaxis (N=48).

---

## 18. Core Contribution Statement

The pmsimstats-ng project establishes that the choice between direct
mean moderation (Architecture A) and differential correlation
(Architecture B) as the DGP for biomarker-treatment interaction in
aggregated N-of-1 designs determines whether carryover effects
reduce statistical power modestly (approximately 10 to 15 percent,
Architecture A) or substantially (approximately 40 to 60 percent,
Architecture B). The choice is a claim about biology: Architecture
A for biomarkers that determine the magnitude of the drug's
biological effect (renal clearance, receptor density, CYP2D6 status);
Architecture B for biomarkers that predict response class through
shared variance with a drug-responsive physiological component
(blood pressure as PTSD subtype marker, baseline severity,
inflammatory markers). Most published simulation studies use
Architecture A. The Hendrickson et al. (2020) framework, which
pmsimstats-ng revises and extends, is unique in the N-of-1 literature
in using Architecture B. Published power estimates for biomarker-
moderated crossover and N-of-1 designs may be systematically
optimistic for biomarkers operating through an Architecture B
mechanism.

The revised implementation embodies seven DGP corrections relative
to the publication (AR(1) autocorrelation, exponential BM-BR decay,
nlme with corCAR1, scale-factor removal, timepoint-1 guard removal,
origin-passing Gompertz, updated defaults), supports both
architectures via the `dgp_architecture` parameter, validates Type
I error at nominal 5 percent across all designs, and executes 1,000
replicates in 47 minutes with zero silent positive-definiteness
corrections.

---

*Rendered on 2026-04-21 at 15:07 PDT.*<br>
*Source: ~/Dropbox/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/compendium-2026-04-21.md*
