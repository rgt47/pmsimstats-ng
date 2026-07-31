---
geometry: margin=0.85in
fontsize: 10pt
header-includes:
  - \linespread{0.96}
  - \setlength{\parskip}{0.4em}
---

# White paper: notation and terminology across the compendium

*2026-07-30 20:10 PDT (second edition; supersedes the 2026-07-29
diagnosis)*

**A re-audit of mathematical notation and controlled vocabulary in the
manuscripts under `analysis/report/`, after two repair passes. Records
what was fixed, what the repairs themselves broke, and what remains.**

## Method

Sixteen manuscript sources were searched for the symbols carrying the
programme's substantive quantities and for the labels naming its
designs, architectures, specifications, and studies. The checks are now
mechanized as a linter (`notation-lint.pl`) that skips fenced R chunks
and `verbatim` blocks, strips `\texttt{}` and backtick spans, and
reports violations of the rules in `analysis/report/NOTATION.md`. Every
manuscript was re-rendered after each change; counts below are from the
sources and from `pdftotext` on the rendered output.

The corpus is smaller than in the first edition. Papers 01 and 02 have
since been reorganized: their narrative variants were promoted to
`report.Rmd` and the earlier full and slim versions moved to
per-paper `archive/` directories. This audit reflects the promoted
files.

## State after the second pass

The linter reports **clean** for all sixteen manuscripts except one
known false positive (a bare `Sx` inside a backtick span that crosses a
line break). All sixteen render. The nine inconsistencies of the first
edition now stand as follows.

| # | Finding (first edition) | State |
|---|---|---|
| N1 | Drug indicator: $D_{it}$ used for both the binary state and the continuous exposure | Fixed |
| N2 | $c_{bm}$ used for both a correlation and a regression multiplier | Fixed |
| N3 | Half-life quoted in weeks and in days with no conversion | Fixed by units paragraphs |
| N4 | Four symbols for one estimand | Standardized on $\beta_{bm:D}$ |
| N5 | Outcome alternates between $Y$ and the code identifier `Sx` | Fixed in model definitions |
| N6 | Label collisions across architectures, specifications, sweeps, studies | Fixed, after one false start |
| N7 | Two sample-size conventions (per path, total) | Documented, not resolved |
| N8 | Design labels alternate abbreviation and prose | Open (editorial) |
| N9 | Sign convention stated only in one paper | Open (editorial) |

## What the second pass found that the first pass missed

**The specification relabel was incomplete.** The first pass renamed
the three carryover specifications A1/A2/A3 to M1/M2/M3 in the
carryover manuscript and its variants, but four other papers cite those
specifications and were not touched. For a day the compendium said M2
in one paper and A2 in four others, which is worse than the original
state, because a reader cross-referencing them would conclude they were
different objects.

**The chosen replacement collided with an existing scheme.** The
calibration manuscript already uses M0 through M3 for its
working-covariance model ladder, and those labels are keys into
archived simulation output, not free text. M1/M2/M3 for specifications
therefore reintroduced exactly the class of collision the pass existed
to remove. The specifications are now **E1/E2/E3**, after E for the
exposure regressor that distinguishes them; E is unused anywhere else
in the corpus.

**The architectures were the deeper problem.** A, B, and C are opaque
(a reader must memorize which is which) and they are the reason the
specification labels collided in the first place. They are now named by
the channel that carries the interaction:

| Prose | Legends and tables | Data value |
|---|---|---|
| the mean-moderation architecture | Mean | `mean_moderation` |
| the covariance-moderation architecture | Covariance | `mvn` |
| the dual-channel architecture | Dual | `combined` |

This touched 581 label occurrences across 16 files plus the figure
drivers. Data values are unchanged, so no simulation was re-run; the
figures were regenerated from existing summaries.

**A second British-spelling family survived the first pass**, because
the first pass enumerated word stems rather than patterns: `flavour`,
`rigour`, `endeavour`, `labour`, `tumour`. Fixed in four files.

## What the repairs themselves broke

Reporting this is the point of a second audit. Three defects were
introduced by the first pass and caught here.

**Data-side labels were renamed along with display labels.** The
specification rename reached inside R chunks, changing
`method_levels <- c('A1', ...)` and three `filter(spec %in% c('A2',
'A3'))` calls. Factor levels and filters index archived data values;
renaming them silently produced empty selections and `NA` rows rather
than an error. Fixed by keeping data values at A1/A2/A3 and mapping to
E1/E2/E3 at display time through a `spec_display()` helper, with the
mapping documented in each paper's reproducibility section.

**A display helper returned `NA` for unmapped inputs.** The first
version of `spec_display()` indexed a named vector, so the two
robust-inference methods (`lme+CR2`, `GEE+MD`) that are not
specifications came back `NA` and printed as `NA` in the
development-review tables. Now passes unmapped values through.

**Three manuscript files were truncated to zero bytes.** A shell loop
copied each source to a temporary file, transformed it, and moved the
result back. When the cloud file provider returned a spurious
`ENOENT` for a source that had just been moved by a concurrent
reorganization, the `cp` failed but the shell had already created an
empty output file, which `mv` then placed over the source. The files
were restored from scratchpad backups. Every subsequent script
verifies that the source is non-empty before producing output and that
the result has not shrunk before writing back. This is a sharper
version of the hazard the project guardrails already warn about: the
danger is not only in-place editors, but any pipeline whose failure
mode is an empty file.

## What remains open

**N7, the sample-size convention.** One paper reports $N$ per
randomization path and others report it in total. The combined-DGP
manuscript documents its own instance and defers to a pending
boundary re-run. This needs a simulation, not an edit.

**N8 and N9** are editorial: design labels still alternate between
abbreviation and spelled-out form, and the sign convention is stated
in only one paper. Both are best folded into the `rgt` writing pass
rather than done mechanically now.

**The `rgt` prose remains unwritten** in every paper that still carries
the three-part scaffold, so a third terminology pass will be needed
once the author's text exists. That is the argument for keeping the
linter: it converts these two one-off passes into a check that can run
on every render.

**A caution on scope.** This audit verifies internal consistency of
symbols and labels. It does not verify that any numerical result is
correct, and the label changes were deliberately confined to display
so that no stored result was touched. Two papers were promoted to
narrative variants during the pass; their archived predecessors were
not re-audited and will not match this scheme if they are ever
revived.

\newpage

# Appendix A. Mathematical symbols

Extracted from the manuscript sources and reconciled against
`analysis/report/NOTATION.md`, which is the canonical reference. A
symbol is listed here only if it appears in at least one manuscript.

## Indices and structure

| Symbol | Meaning |
|---|---|
| $i$ | participant index |
| $t$ | timepoint index; weekly measurement occasions throughout |
| $j$, $ij$ | period or occasion index within participant (papers 04, 06) |
| $N$ | number of participants. **State whether per randomization path or in total**; both conventions are in use |
| $n_{\text{reps}}$, $n_{\text{sim}}$ | Monte Carlo replicates per cell |
| $k$ | cycles per participant (design sweeps); also the Weibull shape parameter in the decay family |

## Outcome and latent components

| Symbol | Meaning |
|---|---|
| $Y_{it}$ | observed symptom score; the modelled outcome |
| $\mathrm{BL}_i$ | participant-specific baseline level |
| $BR_{it}$ | biological (pharmacological) response component |
| $PB_{it}$ | placebo-belief response component |
| $TV_{it}$ | time-variant natural-history component |
| $\varepsilon_{it}$ | observation-level residual |
| $u_i$ | participant random intercept |
| $\eta$ | expectancy covariate (belief measure), paper 06 |

Sign convention: the three components are non-negative *reductions* in
symptom severity, so an increase in any component lowers $Y$.
Treatment effects and interaction coefficients are consequently
negative; moderation parameters are positive.

## Treatment exposure and carryover

| Symbol | Meaning |
|---|---|
| $D_{it}$ | binary drug state: 1 on drug, 0 off drug |
| $D_{bc,it}$ | continuous exposure-decayed drug indicator: 1 on drug, $e^{-\lambda t_{sd}}$ off drug |
| $t_{sd}$ | time since discontinuation (input to the carryover decay) |
| $t_{od}$ | time on drug, cumulative (input to the Gompertz response) |
| $t_{pb}$ | time-by-placebo-belief accumulation |
| $t_{1/2}$ | carryover half-life. **Weeks** in the interaction manuscripts, **days** in the main-effect manuscripts; each paper carries a units paragraph |
| $\lambda$ | carryover decay rate, $\lambda = \ln 2 / t_{1/2}$ |
| $\lambda_w$ | Weibull rate, $(\ln 2)^{1/k} / t_{1/2}$ |

Decay families implemented on the data-generating side: exponential
$e^{-\lambda t}$, linear $\max(0, 1 - t/(2t_{1/2}))$, Weibull
$e^{-(\lambda_w t)^k}$, and power $\max(0, (1 - t/(3t_{1/2}))^p)$.

## Biomarker and moderation

| Symbol | Meaning |
|---|---|
| $B_i$ | participant-level pre-treatment biomarker value |
| $b_i$ | standardized biomarker, $(B_i - \bar{B})/s_B$ |
| $c_{bm}$ | moderation parameter of the **covariance-moderation** architecture: the on-drug correlation between biomarker and $BR$. A correlation, bounded in $[-1,1]$ |
| $\beta_{bm}$ | moderation parameter of the **mean-moderation** architecture: the dimensionless multiplier scaling the treatment effect by $b_i$ |
| $c_{bm,A}$, $c_{bm,B}$ | independent mean-channel and covariance-channel weights of the dual-channel architecture |
| $\sigma_{BR}$, $\sigma_{bm}$ | standard deviations of the response component and the biomarker |

$c_{bm}$ and $\beta_{bm}$ are calibrated to a common numeric scale
(the mean-channel shift is applied as
$\beta_{bm}\,b_i\,\sigma_{BR}$), so the reference value $0.45$ denotes
a matched moderation strength under either architecture. The symbols
are not interchangeable.

## Estimand, inference, and performance

| Symbol | Meaning |
|---|---|
| $\beta_{bm:D}$ | the target estimand: the biomarker-by-treatment interaction coefficient (`bm:Dbc`) |
| $\beta_{bm}^{BR}$, $\beta_{bm}^{PB}$, $\beta_{bm}^{TV}$ | component-specific biomarker slopes (paper 06) |
| $\beta_{bm}^{\text{lumped}}$ | the one-component (single-indicator) biomarker slope |
| $\beta_D$ | treatment main effect (papers 04, 05) |
| $\theta$, $\hat\theta$ | a generic interaction coefficient; used where the argument is deliberately not N-of-1 specific (paper 10) |
| $\theta_{\text{true}}$ | calibrated true value, $-c_{bm}\sigma_{BR}$ |
| $\rho$ | AR(1) within-factor serial correlation |
| $\kappa$ | ratio of mean model-based standard error to empirical standard deviation; $\kappa > 1$ indicates a conservative test |
| $\delta$ | standardized bias of the test statistic |
| $\nu$ | denominator degrees of freedom |
| $\alpha$ | nominal significance level, 0.05 throughout |
| $\pi$ | power, $\Pr(p < \alpha)$ |
| $\gamma$ | latent-class gating slope (paper 03); also a $BR \times PB$ interaction coefficient (paper 06) |
| $Z_i$ | unobserved latent class label (paper 03) |

# Appendix B. Code identifiers and stored data values

These are the names as they appear in `R/`, in the drivers under
`analysis/scripts/`, and in the archived `.rds` summaries. They are
deliberately **not** renamed when reporting labels change, so that
stored results stay readable; manuscripts map them at display time.

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
| `blparam` | biomarker and nuisance means and SDs |
| `trialdesign` | output of `buildtrialdesign`: `timepoints`, `timeptnames`, `expectancies`, `ondrug` |
| `dgp_architecture` | `'mvn'`, `'mean_moderation'`, or `'combined'` |
| `carryover_t1half` | carryover half-life |
| `lambda_cor` | decay rate; derived from the half-life when `NA` |
| `br_family` | response-curve family: `'gompertz'` and alternatives |
| `br_p2`, `br_p3` | shape parameters of the alternative families |
| `moderation_scaling` | `'constant'` or trajectory-scaled |
| `cached_sigma` | pre-built covariance matrix, reused across workers |
| `makePositiveDefinite` | positive-definiteness repair switch |
| `n_cores` | worker count; `-1` auto-detects |

## Data columns

`ptID` (participant), `Sx` (outcome, the $Y_{it}$ of Appendix A), `bm`
(biomarker), `bm_z` (standardized biomarker), `Dbc` (continuous
exposure indicator, $D_{bc,it}$), `tsd`, `tod`, `tpb`, `t_wk` (time in
weeks), `path` (randomization path).

## Stored factor levels

| Column | Stored values | Reported as |
|---|---|---|
| `design` | `CO`, `Hybrid`, `OLBDC` | CO, Hybrid, OL+BDC |
| `dgp_arch` | `mean_moderation`, `mvn`, `combined` | Mean, Covariance, Dual |
| `spec` | `A1`, `A2`, `A3` | E1, E2, E3 |
| `carryover_form` | `exponential`, `linear`, `weibull` | as stored |

## Summary-table columns

`n_reps`, `n_converged`, `converged_frac`, `non_convergence_rate`,
`power`, `bias`, `mean_estimate`, `true_value`, `empirical_se`, `mse`,
`coverage`, and an `mc_se_` prefixed Monte Carlo standard error for
each (`mc_se_power`, `mc_se_bias`, `mc_se_mse`, `mc_se_coverage`,
`mc_se_non_convergence`).

# Appendix C. Glossary of controlled vocabulary

The project lexicon at
`docs/29-nof1-precision-medicine-lexicon.md` holds the full
terminology set (241 terms across 13 thematic sections, with
cross-field synonyms). The entries below are the subset that carries a
fixed, enforced meaning in this compendium.

## Trial designs

**OL** (open-label). Unblinded titration; no off-drug contrast, so it
supplies no within-subject interaction information.
**BDC** (blinded discontinuation). A blinded switch to placebo after
open-label response.
**OL+BDC**. Open-label titration followed by blinded discontinuation.
**CO** (crossover). Traditional two-period within-subject crossover.
**Hybrid**. The Hendrickson (2020) design: open-label titration,
blinded discontinuation, then a brief crossover. The programme's
reference design.
**Aggregated N-of-1 trial**. Many single-patient on/off series pooled
into one mixed model.

## Data-generating architectures

**Mean-moderation architecture** (`mean_moderation`). The biomarker
scales the treatment effect in the conditional mean; the interaction
lives in the first moment.
**Covariance-moderation architecture** (`mvn`). The biomarker-response
correlation is treatment-state dependent; the interaction lives in the
second moment. The Hendrickson formulation.
**Dual-channel architecture** (`combined`). Both channels active with
independent weights; the other two are its single-channel limits.

## Analysis specifications for carryover

**E1**. Binary on-drug indicator; carryover ignored.
**E2**. Exposure-weighted continuous predictor $D_{bc}$, committing to
an assumed decay form and half-life.
**E3**. Binary indicator plus a lagged just-off-drug nuisance term;
the classical Jones-Kenward crossover device. Assumes no decay form.

## Response-curve families

**Modified Gompertz** (default), $y = \mathrm{maxr}\cdot
\exp(-\mathrm{disp}\cdot\exp(-\mathrm{rate}\cdot t))$ with a vertical
offset so the curve passes through the origin. Alternatives:
**logistic** (symmetric sigmoidal), **hyperbolic tangent**
(fast-saturating), **piecewise-linear breakpoint** (no smoothness
assumption).

## Missing data and dropout

**MCAR**, **MAR**, **MNAR** in the standard Rubin sense.
**Informative dropout**. Hazard depending on cumulative symptom
worsening since baseline.
**Happy accident**. The randomization-path selection effect by which
informative dropout preferentially preserves the most informative
crossover blocks.

## Inference

**LME**. Linear mixed-effects model; the programme's default analysis.
**corCAR1**. Continuous-time AR(1) residual correlation structure.
**GEE**. Generalized estimating equations; the marginal comparator.
**Mancl-DeRouen**, **CR2**. Small-sample bias-corrected sandwich
variance estimators.
**RM-ANOVA**. Repeated-measures analysis of variance; requires a
dichotomized biomarker in the strict form.
**Anti-conservative**. Realized Type I error above nominal.

## Estimands and reporting

**Biomarker-treatment interaction**. The programme's core estimand;
equivalently the patient-by-treatment interaction, or heterogeneity of
treatment effects (HTE) in the precision-medicine literature.
**Predictive biomarker**. Predicts differential benefit, as opposed to
a **prognostic** biomarker, which predicts outcome regardless of
treatment.
**ADEMP**. The Morris, White and Crowther (2019) reporting framework:
aims, data-generating mechanisms, estimands, methods, performance
measures.
**MCSE**. Monte Carlo standard error, reported with every performance
measure.

---
*Rendered on 2026-07-30 at 20:10 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/whitepaper-notation-audit.md*
