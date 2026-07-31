# Canonical notation, identifiers, and glossary

*Originally written 2026-07-30 11:57 PDT; expanded to a full reference
2026-07-30 20:20 PDT*

This file is the single source of truth for mathematical notation,
code identifiers, and controlled vocabulary across the manuscripts
under `analysis/report/`. New manuscripts should adopt it rather than
re-deriving symbols, and any deviation should be stated explicitly in
the deviating paper's Methods section.

The audit behind it is `analysis/report/whitepaper-notation-audit.md`,
whose PDF binds this file in as its appendix so that the two cannot
drift apart.

Symbols and identifiers listed here were extracted from the manuscript
sources, from `R/`, from the drivers under `analysis/scripts/`, and
from the archived `.rds` summaries.

# Part 1. Mathematical symbols

## Indices and structure

| Symbol | Meaning |
|---|---|
| $i$ | participant index |
| $t$ | timepoint index; weekly measurement occasions throughout |
| $j$, $ij$ | period or occasion index within participant (papers 04, 06) |
| $N$ | number of participants. **State whether per randomization path or in total**; both conventions are in use |
| $n_{\text{reps}}$, $n_{\text{sim}}$ | Monte Carlo replicates per cell |
| $k$ | cycles per participant in design sweeps; also the Weibull shape parameter in the decay family |

## Outcome and latent components

| Symbol | Meaning |
|---|---|
| $Y_{it}$ | observed symptom score; the modelled outcome (code column `Sx`) |
| $\mathrm{BL}_i$ | participant-specific baseline level |
| $BR_{it}$ | biological (pharmacological) response component |
| $PB_{it}$ | placebo-belief response component |
| $TV_{it}$ | time-variant natural-history component |
| $\varepsilon_{it}$ | observation-level residual |
| $u_i$ | participant random intercept |
| $\eta$ | expectancy covariate (measured belief), paper 06 |

## Treatment exposure and carryover

| Symbol | Meaning |
|---|---|
| $D_{it}$ | binary drug state: 1 on drug, 0 off drug |
| $D_{bc,it}$ | continuous exposure-decayed drug indicator: 1 on drug, $e^{-\lambda t_{sd}}$ off drug (code column `Dbc`) |
| $t_{sd}$ | time since discontinuation; input to the carryover decay (code `tsd`) |
| $t_{od}$ | cumulative time on drug; input to the Gompertz response (code `tod`) |
| $t_{pb}$ | time-by-placebo-belief accumulation (code `tpb`) |
| $t_{1/2}$ | carryover half-life (code `carryover_t1half`) |
| $\lambda$ | carryover decay rate, $\ln 2 / t_{1/2}$ (code `lambda_cor`) |
| $\lambda_w$ | Weibull rate, $(\ln 2)^{1/k} / t_{1/2}$ |

Decay families implemented on the data-generating side: exponential
$e^{-\lambda t}$; linear $\max(0, 1 - t/(2t_{1/2}))$; Weibull
$e^{-(\lambda_w t)^k}$; power $\max(0, (1 - t/(3t_{1/2}))^p)$.

## Biomarker and moderation

| Symbol | Meaning |
|---|---|
| $B_i$ | participant-level pre-treatment biomarker value (code `bm`) |
| $b_i$ | standardized biomarker, $(B_i - \bar{B})/s_B$ (code `bm_z`) |
| $c_{bm}$ | moderation parameter of the **covariance-moderation** architecture: on-drug correlation between biomarker and $BR$. A correlation, bounded in $[-1,1]$ |
| $\beta_{bm}$ | moderation parameter of the **mean-moderation** architecture: dimensionless multiplier scaling the treatment effect by $b_i$ |
| $c_{bm,A}$, $c_{bm,B}$ | independent mean-channel and covariance-channel weights of the dual-channel architecture |
| $\sigma_{BR}$, $\sigma_{bm}$ | standard deviations of the response component and of the biomarker |

## Estimand, inference, and performance

| Symbol | Meaning |
|---|---|
| $\beta_{bm:D}$ | the target estimand: the biomarker-by-treatment interaction coefficient (`bm:Dbc`) |
| $\beta_{bm}^{BR}$, $\beta_{bm}^{PB}$, $\beta_{bm}^{TV}$ | component-specific biomarker slopes (paper 06) |
| $\beta_{bm}^{\text{lumped}}$ | the one-component (single-indicator) biomarker slope |
| $\beta_D$ | treatment main effect (papers 04, 05) |
| $\theta$, $\hat\theta$ | a generic interaction coefficient; used only where the argument is deliberately not N-of-1 specific (paper 10) |
| $\theta_{\text{true}}$ | calibrated true value, $-c_{bm}\,\sigma_{BR}$ |
| $\rho$ | AR(1) within-factor serial correlation |
| $\kappa$ | ratio of mean model-based standard error to empirical standard deviation; $\kappa > 1$ means a conservative test |
| $\delta$ | standardized bias of the test statistic |
| $\nu$ | denominator degrees of freedom |
| $\alpha$ | nominal significance level; 0.05 throughout |
| $\pi$ | power, $\Pr(p < \alpha)$ |
| $\gamma$ | latent-class gating slope (paper 03); a $BR \times PB$ interaction coefficient (paper 06) |
| $Z_i$ | unobserved latent class label (paper 03) |

## The five rules the symbols encode

1. **$D_{it}$ is binary; $D_{bc,it}$ is continuous.** These are
   different objects, and the distinction is the subject of the
   carryover-specification manuscript. Never use one for the other.
2. **$c_{bm}$ is a correlation; $\beta_{bm}$ is a multiplier.** They
   are calibrated to a common numeric scale, because the mean-channel
   shift is applied as $\beta_{bm}\,b_i\,\sigma_{BR}$ on the
   standardized biomarker, so the reference value $0.45$ denotes a
   matched moderation strength under either architecture. The symbols
   are nonetheless not interchangeable. A paper reporting both
   architectures at a matched value may use one shared label, but must
   say so.
3. **Mathematics in model definitions, code identifiers in code
   contexts.** Write $Y_{ij} = \beta_0 + \beta_D D_{bc,ij} + \ldots$
   for the model; write `Sx ~ bm + t + Dbc + bm:Dbc` in `\texttt{}`
   when quoting the R formula. Never mix the two inside one
   expression.
4. **State the sample-size convention.** $N$ per randomization path
   and $N$ in total are both in use and are not interchangeable; every
   Methods section must say which it means.
5. **State the sign convention.** The three components are
   non-negative *reductions* in symptom severity, so an increase in
   any component lowers $Y$. Treatment effects and interaction
   coefficients are consequently negative; moderation parameters are
   positive.

## Units

Carryover half-life $t_{1/2}$ is quoted in **weeks** in the
interaction-focused manuscripts, on the canonical grid
$\{0, 0.5, 1.0\}$, and in **days** in the main-effect and
test-procedure manuscripts, on a pharmacokinetic scale with a 3-day
baseline. Both conventions are retained because each paper's figures
are drawn on its own scale. Every paper using days carries an explicit
units paragraph giving the conversion (1 week = 7 days). Timepoints
$t$ are weekly in both conventions.

# Part 2. Code identifiers and stored data values

These are the names as they appear in `R/`, in the drivers, and in the
archived `.rds` summaries.

Where a stored value is also a reporting label, the two are kept
**identical**, so that no display-time mapping stands between the data
and the page. A mapping layer is a standing hazard: it has to be
applied at every display site, it fails silently when a rename reaches
one layer and not the other, and it makes a stored table unreadable
without the manuscript beside it. Where a stored value is a package
API argument (`dgp_architecture = 'mvn'`), it stays as the API defines
it and the manuscript maps it to a reporting name, which is the
ordinary factor-level-to-label relationship.

## Exported package functions

`buildtrialdesign`, `buildSigma`, `generateData`, `censordata`,
`lme_analysis`, `generateSimulatedResults`, `validateParameterGrid`,
`cumulative`, `modgompertz`, `trajectoryShape`, `shape_logistic`,
`shape_hyperbolic_tangent`, `shape_piecewise_linear_breakpoint`,
`characterize_carryover`, `analyze_trial_extended`,
`print_carryover_summary`, `print_trial_summary`,
`plotfactortrajectories`, `PlotModelingResults`, `reknitsimresults`.

## Principal arguments

| Identifier | Role |
|---|---|
| `modelparam` | list: sample size `N`, moderation `c.bm` |
| `respparam` | Gompertz parameters `maxr`, `rate`, `disp` per component |
| `blparam` | biomarker and nuisance means and standard deviations |
| `trialdesign` | output of `buildtrialdesign`: `timepoints`, `timeptnames`, `expectancies`, `ondrug` |
| `dgp_architecture` | `'mvn'`, `'mean_moderation'`, or `'combined'` |
| `carryover_t1half` | carryover half-life |
| `lambda_cor` | decay rate; derived from the half-life when `NA` |
| `br_family` | response-curve family: `'gompertz'` and alternatives |
| `br_p2`, `br_p3` | shape parameters of the alternative families |
| `moderation_scaling` | `'constant'` or trajectory-scaled |
| `cached_sigma` | pre-built covariance matrix, reused across parallel workers |
| `makePositiveDefinite` | positive-definiteness repair switch |
| `n_cores` | worker count; `-1` auto-detects |
| `save_chunks` | progressive checkpointing for long runs |

## Data columns

`ptID` (participant), `Sx` (outcome, the $Y_{it}$ above), `bm`
(biomarker), `bm_z` (standardized biomarker), `Dbc` (continuous
exposure indicator, $D_{bc,it}$), `tsd`, `tod`, `tpb`, `t_wk` (time in
weeks), `path` (randomization path).

## Stored factor levels and their reporting labels

| Column | Stored values | Reported as | Mapped? |
|---|---|---|---|
| `spec` | `E1`, `E2`, `E3` | E1, E2, E3 | no, identical |
| `carryover_form` | `exponential`, `linear`, `weibull` | as stored | no, identical |
| `design` | `CO`, `Hybrid`, `OLBDC` | CO, Hybrid, OL+BDC | only `OLBDC` to `OL+BDC` |
| `dgp_arch` | `mean_moderation`, `mvn`, `combined` | Mean, Covariance, Dual | yes; these are package API values |

The `spec` column stored `A1`, `A2`, `A3` until 2026-07-31. Twenty-one
`.rds` files were migrated to `E1/E2/E3` and the producing drivers
updated to emit them, so the display-time mapping that previously
bridged the two has been retired. Pre-migration copies are in
`analysis/.spec-migration-backup/`. Output archived before that date,
including anything under a paper's own `archive/`, still uses the old
values.

## Summary-table columns

`n_reps`, `n_converged`, `converged_frac`, `non_convergence_rate`,
`power`, `bias`, `mean_estimate`, `true_value`, `empirical_se`, `mse`,
`coverage`, and an `mc_se_`-prefixed Monte Carlo standard error for
each (`mc_se_power`, `mc_se_bias`, `mc_se_mse`, `mc_se_coverage`,
`mc_se_non_convergence`).

# Part 3. Controlled vocabulary and glossary

The project lexicon at `docs/29-nof1-precision-medicine-lexicon.md`
holds the full terminology set (241 terms across 13 thematic sections,
with cross-field synonyms). The entries below are the subset carrying
a fixed, enforced meaning in this compendium.

## Trial designs

**OL** (open-label). Unblinded titration; no off-drug contrast, so it
supplies no within-subject interaction information.
**BDC** (blinded discontinuation). A blinded switch to placebo after
open-label response.
**OL+BDC**. Open-label titration followed by blinded discontinuation.
**CO** (crossover). Traditional two-period within-subject crossover.
**Hybrid**. The Hendrickson (2020) design: open-label titration,
blinded discontinuation, then a brief crossover. The programme's
reference design. Capitalize it.
**Aggregated N-of-1 trial**. Many single-patient on/off series pooled
into one mixed model.

Define each abbreviation at first use in every paper.

## Data-generating architectures

Named by the channel that carries the interaction. Single letters are
no longer used.

| Name in prose | Legends and tables | Data value | Former label |
|---|---|---|---|
| the mean-moderation architecture | Mean | `mean_moderation` | Architecture A |
| the covariance-moderation architecture | Covariance | `mvn` | Architecture B |
| the dual-channel architecture | Dual | `combined` | Architecture C |

**Mean-moderation.** The biomarker scales the treatment effect in the
conditional mean; the interaction lives in the first moment. The
near-universal convention in the trial-simulation literature.
**Covariance-moderation.** The biomarker-response correlation is
treatment-state dependent; the interaction lives in the second moment.
The Hendrickson formulation, and the one whose power is most eroded by
carryover.
**Dual-channel.** Both channels active with independent weights
$c_{bm,A}$ and $c_{bm,B}$; the other two are its single-channel
limits. A data-generating construct only: the two weights are not
separately identifiable from a fitted `bm:Dbc` coefficient.

## Analysis specifications for carryover

**E1.** Binary on-drug indicator; carryover ignored.
**E2.** Exposure-weighted continuous predictor $D_{bc}$, committing to
an assumed decay form and half-life.
**E3.** Binary indicator plus a lagged just-off-drug nuisance term;
the classical Jones-Kenward crossover device, assuming no decay form.

The letter E stands for the exposure regressor that distinguishes
them. These were labeled A1/A2/A3 until 2026-07-31; the stored data
and the drivers were migrated to E1/E2/E3, so the labels now agree
everywhere and nothing is mapped at display time.

**Why not A/B/C and M1/M2/M3.** The original scheme used A/B/C for
architectures and A1/A2/A3 for specifications, colliding on the letter
A within a single manuscript. Relabeling the specifications M1/M2/M3
removed that collision but created a second one: the calibration
manuscript already uses M0-M3 for its working-covariance model ladder,
and those labels are keys into archived data. E1/E2/E3 is free across
the whole corpus, and naming the architectures by mechanism removes
the class of collision entirely.

## Response-curve families

**Modified Gompertz** (default), $y = \mathrm{maxr}\cdot
\exp(-\mathrm{disp}\cdot\exp(-\mathrm{rate}\cdot t))$ with a vertical
offset so the curve passes through the origin. Alternatives:
**logistic** (symmetric sigmoidal), **hyperbolic tangent**
(fast-saturating), **piecewise-linear breakpoint** (no smoothness
assumption).

## Missing data and dropout

**MCAR**, **MAR**, **MNAR** in the standard Rubin sense.
**Informative dropout.** Hazard depending on cumulative symptom
worsening since baseline.
**Happy accident.** The randomization-path selection effect by which
informative dropout preferentially preserves the most informative
crossover blocks.

## Inference

**LME.** Linear mixed-effects model; the programme's default analysis.
**corCAR1.** Continuous-time AR(1) residual correlation structure.
**GEE.** Generalized estimating equations; the marginal comparator.
**Mancl-DeRouen**, **CR2.** Small-sample bias-corrected sandwich
variance estimators.
**RM-ANOVA.** Repeated-measures analysis of variance; the strict form
requires a dichotomized biomarker.
**Anti-conservative.** Realized Type I error above nominal.
**Working covariance.** The assumed residual covariance of the LME;
reserve **working correlation structure** for the GEE object.

## Estimands and reporting

**Biomarker-treatment interaction.** The programme's core estimand;
equivalently the patient-by-treatment interaction, or heterogeneity of
treatment effects (HTE) in the precision-medicine literature.
**Predictive biomarker.** Predicts differential benefit, as opposed to
a **prognostic** biomarker, which predicts outcome regardless of
treatment.
**ADEMP.** The Morris, White and Crowther (2019) reporting framework:
aims, data-generating mechanisms, estimands, methods, performance
measures.
**MCSE.** Monte Carlo standard error; reported with every performance
measure.

## Enumeration conventions

**Simulation studies.** Number them (Study 1, Study 2, ...) rather
than lettering them.
**Sensitivity sweeps.** Labels S1, S2, ... are paper-local and are not
comparable across papers. Cite them with the paper prefix when
crossing manuscripts (05-S2, not S2).

## Prose terminology

Symbol conventions are settled here; prose terminology is settled by
`docs/30-terminology-consistency-plan.md`, executed 2026-07-30 with
its D1 spelling decision locked to **US English**, overriding that
document's original frequency-based recommendation of British forms.
The two are complementary: this file governs mathematics and
identifiers, that one governs surface forms and acronym expansion.

