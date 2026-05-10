# ADEMP pre-registration: paper 02 (carryover-sensitivity)
*2026-05-08 07:35 PDT*

**Author.** pmsimstats team

**Manuscript.** `analysis/report/02-carryover-sensitivity/report.Rmd`

**Pre-registration date.** 2026-04-15 (original grid-plan;
reformatted into ADEMP structure 2026-05-08).

This document pre-registers the simulation programme described in
the carryover-sensitivity manuscript (paper 02) following the
ADEMP reporting framework of @morris2019simulation. It restates
the parameter grid, performance measures, and replicate budget
specified in `simulation-grid-plan.md` (the original design
document) in the formal Aims / Data-generating-mechanism /
Estimands / Methods / Performance-measures structure SIM editors
expect. Items in this document are committed prior to the
production runs reported in the manuscript Results section; any
post-hoc deviations are documented at the foot of this
pre-registration.

## A. Aims

**A1 (principal aim).** Determine whether the ranking of
analysis-model carryover specifications (A1 binary, A2
exposure-weighted, A3 lagged-nuisance) for detecting the
biomarker-treatment interaction in aggregated N-of-1 trials
depends on (i) the true data-generating carryover decay form,
(ii) the data-generating-process architecture (mean moderation
versus multivariate-normal differential correlation), and
(iii) the carryover half-life, across the trial-design and
sample-size settings used in precision-medicine N-of-1 work.

**A2 (sensitivity aims).** Evaluate the robustness of the A1
ranking to (S1) within-factor autocorrelation strength,
(S2) analyst-versus-truth half-life mismatch, (S3) informative
and non-informative dropout, (S4) the full effect-size curve
in the biomarker moderation parameter, and (S5) interactions
between autocorrelation strength and carryover half-life.

**A3 (secondary aim).** Quantify the absolute power penalty of
analysis-model carryover mis-specification under each DGP
architecture. The A3 outputs feed sample-size guidance in the
Discussion.

## D. Data-generating mechanisms

The simulation generates aggregated N-of-1 trial data using the
pmsimstats-ng `tidyverse` simulation collection
(`implementations/tidyverse/`). The data-generating model has
four conceptual components:

- **Trial design.** One of `CO` (traditional crossover), `Hybrid`
  (Hendrickson hybrid), and `OL+BDC` (open-label with blinded
  discontinuation), built by `buildtrialdesign()`.
- **Response trajectory.** Three modGompertz components
  ($BR$ pharmacological, $TV$ time-variant, $PB$ placebo-belief),
  parameterised at the prazosin/PTSD reference cell.
- **Biomarker moderation.** Either Architecture A (direct
  mean-moderation: post-draw shift on on-drug $BR$ scaled by
  the centred biomarker) or Architecture B (MVN differential
  correlation: $\mathrm{Cor}(BM, BR)$ on-drug equals $c_{bm}$;
  off-drug decays as $c_{bm} \exp(-\lambda t_{sd})$).
- **Carryover decay.** One of linear, exponential, or Weibull
  with shape parameter $k \in \{0.7, 1.0, 1.5\}$.

### Tier 1 principal factorial (540 cells)

Every combination of the seven axes below is evaluated:

| Axis | Values | Levels |
|---|---|---|
| DGP carryover form | linear, exponential, Weibull($k \in \{0.7, 1.0, 1.5\}$) | 5 |
| Analysis specification | A1 binary, A2 exposure-weighted, A3 lagged | 3 |
| DGP architecture | mean moderation, MVN differential correlation | 2 |
| Carryover half-life $t_{1/2}$ (weeks) | 0, 0.5, 1.0 | 3 |
| Trial design | CO, Hybrid, OL+BDC | 3 |
| Sample size $N$ | 35, 70 | 2 |
| Biomarker moderation $c_{bm}$ | 0, 0.30, 0.45 | 3 |

Cell count: $5 \times 3 \times 2 \times 3 \times 3 \times 2
\times 3 = 540$.

### Tier 2 sensitivity blocks (5 blocks, 31 cells)

Each block varies one additional axis at the **reference cell**
$(\text{exponential}, k=1.0, \text{MVN}, t_{1/2}=1.0,
\text{Hybrid}, N=70, c_{bm}=0.45)$ crossed with the three
analysis specifications.

- **S1 autocorrelation.** $\rho \in \{0.5, 0.7, 0.8, 0.9\}$.
  4 cells.
- **S2 half-life mismatch.** $t_{1/2}^{\text{true}} \in \{0.5,
  1.0\}$ × $t_{1/2}^{\text{analyst}} \in \{0.25, 0.5, 1.0,
  2.0\}$. 8 cells.
- **S3 dropout.** Rate $d \in \{0, 0.1, 0.2, 0.3\}$ × mechanism
  $\{$MCAR, MAR-by-baseline-severity$\}$, with the
  $(d=0, \text{MAR})$ cell collapsed to MCAR. 7 cells.
- **S4 effect-size curve.** $c_{bm} \in \{0, 0.10, 0.20, 0.30,
  0.45, 0.60\}$. 6 cells. ($c_{bm} = 0.60$ feasible only under
  Architecture A; reported with note.)
- **S5 rho × carryover (exploratory).** $\rho \in \{0.5, 0.8\}$
  × $t_{1/2} \in \{0, 0.5, 1.0\}$. 6 cells.

Tier 1 and Tier 2 share an anchor at the reference cell; results
are reported separately and not pooled. The reasoning for the
two-tier separation is documented in `simulation-grid-plan.md`
§"Relationship between Tier 1 and Tier 2".

## E. Estimands

The principal target of inference is the
biomarker-by-treatment interaction coefficient
$\beta_{bm{:}D_{bc}}$ in the linear-mixed analysis model
$Sx \sim bm + t + D_{bc} + bm{:}D_{bc}$, with random intercept
per participant and continuous-time AR(1) residual correlation
(`corCAR1(form = ~ t | ptID)`).

Three estimands are evaluated for each cell:

- **E1: Power.** The rejection probability for the $bm{:}D_{bc}$
  Wald $t$-test at $\alpha = 0.05$ under the alternative
  ($c_{bm} > 0$).
- **E2: Type I error.** The rejection probability under the null
  ($c_{bm} = 0$). Reported on the `c_bm = 0` cells in Tier 1
  and Tier 2 blocks where applicable.
- **E3: Estimator bias.** The mean of
  $\hat\beta_{bm{:}D_{bc}}$ minus the implied true value of the
  interaction coefficient under each architecture. Architecture
  A induces a true value of $c_{bm} \cdot \sigma_{BR} /
  \sigma_{bm}$; Architecture B induces a true value derivable
  from the Sigma matrix's
  $\partial \mathrm{Cor}(BM, BR)/\partial D_{bc}$.

Secondary estimands include the empirical standard deviation
of $\hat\beta_{bm{:}D_{bc}}$ across replicates, the median
absolute residual of the model fit, and the convergence rate
of the `nlme::lme` optimiser.

## M. Methods

Three estimators are compared at every cell:

- **A1 binary.** Fit the linear-mixed model with `Db` (binary
  on-drug indicator: 1 on-drug, 0 off-drug). Strict-RM-ANOVA
  analogue.
- **A2 matched.** Fit the linear-mixed model with `Dbc`
  (continuous-decay drug indicator: 1 on-drug,
  $\exp(-\lambda^{\text{analyst}} t_{sd})$ off-drug, where
  $\lambda^{\text{analyst}}$ matches the DGP except in S2).
  pmsimstats-ng default.
- **A3 lagged.** Fit the linear-mixed model with two indicators:
  the binary on-drug `Db` and a lagged-on indicator that equals
  1 at the first off-drug timepoint after each on-drug episode
  and 0 elsewhere. The lagged indicator absorbs short-range
  carryover without committing to a parametric decay form
  [@joneskenward2014].

All three fits use `corCAR1(form = ~ t | ptID)` for the residual
correlation structure, random intercept per participant, REML
estimation, and the `nlme::lme` `optim` optimiser with default
tolerances. Convergence is recorded for each replicate; failed
fits are excluded from the performance summaries (with the
counts reported).

## P. Performance measures

For every cell, the following are reported with their Monte
Carlo standard errors (MCSE):

| Performance measure | Formula | MCSE formula |
|---|---|---|
| Power | $\hat\pi = \frac{1}{n}\sum I(p_i < 0.05)$ | $\sqrt{\hat\pi(1-\hat\pi)/n}$ |
| Type I error (null cells) | as above with $c_{bm} = 0$ | as above |
| Mean estimate | $\bar\beta = \frac{1}{n}\sum \hat\beta_i$ | $s_\beta / \sqrt{n}$ |
| Empirical SE | $s_\beta = \mathrm{sd}(\hat\beta_i)$ | $s_\beta / \sqrt{2(n-1)}$ |
| Bias | $\bar\beta - \beta_{\text{true}}$ | $s_\beta / \sqrt{n}$ |
| Coverage | $\frac{1}{n}\sum I(\beta_{\text{true}} \in \mathrm{CI}_i)$ | $\sqrt{\hat c(1-\hat c)/n}$ |
| Convergence rate | $\frac{1}{n}\sum I(\text{conv}_i = 1)$ | $\sqrt{\hat c(1-\hat c)/n}$ |

For paired cross-specification comparisons within a cell
(A2 versus A1, A2 versus A3), McNemar $\chi^2$ tests on the
per-replicate rejection indicators are reported with their
$p$-values; this leverages the variance reduction from applying
all three specifications to the same simulated datasets.

## Replication budget and justification

The principal factorial uses **500 replicates per cell**, giving
MCSE on a power estimate near 0.5 of approximately $0.022$.
At the manuscript's headline cells (power approximately 0.86)
the MCSE is approximately $0.015$, so the A2-vs-A1 separation
of approximately 10 percentage points is approximately five MCSE
of separation, sufficient for unambiguous discrimination. The
paired McNemar comparisons leverage the within-replicate
correlation among A1, A2, and A3 fits, reducing the effective
SE of a within-cell ranking comparison below the marginal
binomial MCSE.

A high-precision rerun at 600-1000 replicates per cell is
specified for a 24-cell slice (3 specs × 2 $c_{bm}$ × 2 $N$ ×
2 $t_{1/2}^{\text{DGP}}$) to confirm the headline ranking with
MCSE near 0.018, reported as a Prototype Monte Carlo
confirmation in the manuscript Discussion section.

## Convergence and failure reporting

The `nlme::lme` fit at each replicate produces a `conv` flag.
Replicates with $\text{conv} = 0$ are counted but excluded from
the rejection-rate and bias estimators. Per-cell convergence
rates are reported alongside the performance measures.
Cells with convergence rate below 90% are flagged in the
manuscript Results section and discussed in the Discussion.

In the production run as of 2026-05-08, convergence is
$\ge 99.5\%$ across the 540 Tier-1 cells.

## Software, versions, and seeds

- R 4.4.x inside the project's zzcollab Docker container
  (image hash recorded in `Dockerfile`).
- Package versions captured in `renv.lock`. Principal packages:
  `nlme`, `data.table`, `furrr`, `purrr`, `parallel`,
  `rmarkdown`, `bookdown`.
- Master RNG seed: `set.seed(20260415)` at the top of the
  driver `01-run-factorial.R`. Per-replicate seeds derive from
  the master seed by `+ rep_idx`. Tier 2 sensitivity blocks use
  the same convention.

## Pre-specified outputs

| Output | Path |
|---|---|
| Tier 1 per-replicate | `analysis/scripts/carryover-sensitivity/output/01-factorial.rds` |
| Tier 2 per-replicate | `analysis/scripts/carryover-sensitivity/output/04-sensitivity.rds` |
| Tier 1 cell summary | `analysis/data/02-grid-summary.rds` |
| Tier 2 cell summary | `analysis/data/02-sensitivity-summary.rds` |
| High-precision rerun | `analysis/data/quick-sim/02-carryover-replicates.rds` |
| Prototype summary | `analysis/data/quick-sim/02-carryover-summary.txt` |
| Tier 1 figures | `analysis/figures/02-power-by-spec.pdf`, `02-heatmap-matched-vs-mismatched.pdf`, `02-type1-boxplots.pdf` |
| Tier 2 figures | `analysis/figures/02-sens-S{1..5}.pdf` |

## Scope exclusions

The following axes are deliberately excluded from this
pre-registration and are noted in the Discussion section
of the manuscript:

- **Washout period length and timepoint density.** Design-level
  mitigation is orthogonal to the analysis-specification
  question evaluated here. A dedicated paper on design-level
  mitigation is in the research portfolio.
- **Gompertz response trajectory parameters
  (max, rate, disp).** Out of scope; the prazosin/PTSD
  reference parameters are held fixed.
- **Random slopes versus intercept-only random effects.**
  Orthogonal to carryover specification and would double the
  grid; a brief note in the manuscript Discussion points to the
  Hendrickson 2020 sensitivity analysis that addressed this.
- **Multiple-testing correction.** Each simulation tests a
  single pre-specified biomarker.
- **Alpha level.** Type I error is evaluated at the standard
  $\alpha = 0.05$; a supplementary check at $\alpha = 0.10$ is
  available post hoc.

## Post-hoc deviations from this pre-registration

(None to report at the pre-registration date. Any subsequent
deviations should be appended below with date, scope, and
justification.)

---

*Reformatted from `simulation-grid-plan.md` on 2026-05-08
07:35 PDT.*
*Source:
~/Dropbox/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/scripts/carryover-sensitivity/00-ademp-pre-registration.md.*
