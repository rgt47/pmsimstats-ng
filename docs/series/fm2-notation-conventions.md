# Notation, vocabulary and conventions

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Pointer to NOTATION; named processes; analysis IDs;
evidence labels; Monte Carlo conventions; construct box.
**Prerequisites:** None

---

## Contents

1. [Pointer to NOTATION](#pointer-to-notation)
2. [Evidence labels](#evidence-labels)
3. [Monte Carlo conventions](#monte-carlo-conventions)
4. [Analysis-ID scheme](#analysis-id-scheme)
5. [Canonical names](#canonical-names)
6. [Mathematical symbols](#mathematical-symbols)
7. [Code identifiers and stored data values](#code-identifiers-and-stored-data-values)
8. [Controlled vocabulary and glossary](#controlled-vocabulary-and-glossary)
9. [Change log](#change-log)

---

## Pointer to NOTATION

<!-- From analysis/report/NOTATION.md preamble -->

The single source of truth for mathematical notation, code
identifiers, and controlled vocabulary across the manuscripts under
`analysis/report/` is:

`analysis/report/NOTATION.md`

New manuscripts should adopt it rather than re-deriving symbols, and
any deviation should be stated explicitly in the deviating paper's
Methods section. Every canonical symbol used in a chapter must map to
NOTATION by its table name. A local symbol may not use a letter on
NOTATION's reserved list. This chapter holds the project vocabulary
that NOTATION Part 3 holds for the manuscripts (named processes,
configurations, analysis IDs, evidence labels), so that the lexicon
(Appendix E) can stay a cross-field glossary.

The notation of the whitepapers `docs/34` to `docs/58` and of
`analysis/report/compendium-mathematics.tex` was reviewed in five
passes (2026-10-09) and the proposed symbols were merged into
NOTATION.md. Seven user decisions of 2026-10-09 settle the collisions
the passes found:

1. $\tau^2$ is the between-participant variance of the drug effect
   (treatment heterogeneity, as in Senn, paper 04, the paper 04
   working protocol and `docs/58`). The residual variance of the E9
   contrast given the biomarker, previously written $\tau^2$ and
   $\tau_p^2$ in `docs/37` to `docs/39`, the compendium and paper 01,
   is now $\sigma_\Delta^2$ and $\sigma_{\Delta,p}^2$.
2. $\eta_t$ is the design expectancy (1 open label, 0.5 blinded, 0 at
   baseline; code `e`, regressor `De`). A measured placebo-belief
   covariate, if one is used, is $\tilde\eta_{it}$.
3. $\rho$ is the AR(1) correlation per week; the common correlation of
   a compound-symmetry construct is $\rho_{cs}$ (the published process
   has $\rho_{cs} = 0.8$).
4. The linear-time coefficient is $\beta_T$, multiplying the calendar
   week $w_t$.
5. The modified Gompertz curve is written with the mathematical
   symbols $m$ (ceiling), $d$ (displacement) and $r$ (rate); the code
   names are `maxr`, `disp`, `rate`.
6. The share of the prognostic slope absorbed by a shared biomarker
   main effect is $f$. $w_t$ is reserved for the calendar week, so
   row-group counts are $n_{\text{on}}$, $n_{\text{off}}$, $n_v$.
7. The prognostic couplings are $c_{bm,TV}$, $c_{bm,PB}$,
   $c_{bm,BL}$ ($c_{bl}$ is retired).

---

## Evidence labels

<!-- From docs/59-whitepaper-refactoring-plan.md Section 2.5 -->

The series uses seven labels, defined once here and used in every
chapter:

| Label | Use |
|------------|-----------------------------------------------------|
| verified | computed or run, and the output observed; cite the script |
| derived | proved analytically; cite the CM label or the chapter section |
| inspected | read in code or in a source document; cite file and line or commit |
| quoted | taken from another document; cite document and section |
| argued | a reasoned position, not a proof |
| inferred | consistent with the surrounding pattern, not checked |
| unverified | not checked at all |

Table: The evidence-label vocabulary of the series (merging the sets
of `docs/51`, `docs/53`, `docs/54`, `docs/55`, `docs/56`, `docs/57`
and `docs/58`)

---

## Monte Carlo conventions

<!-- From docs/59-whitepaper-refactoring-plan.md Section 2.6 -->

Every chapter applies the following conventions, collected from
`docs/45` Section 2, `docs/46` Section 3, `docs/57` Section 8.1 and
CM "Shared performance measures":

- **Reporting frame.** ADEMP for every study (aims, data-generating
  mechanisms, estimands, methods, performance measures), with the
  construct box of Section 2.10.
- **Size band.** For a pooled null rate over $n$ null trials the band
  is $0.05 \pm 2\sqrt{0.05 \times 0.95 / n}$, stated in the caption
  with $n$; for 12,000 trials, 0.046 to 0.054. Per-cell cap 0.075;
  a prognostic cell above 0.075 raises a flag. These are the screens
  of `docs/46` as amended (quoted, `docs/57` Section 8.1).
- **Distances from nominal.** "$k$ MCSE above 0.05" uses the MCSE at
  0.05, $\sqrt{0.05 \times 0.95 / n_{\text{reps}}}$, as `docs/42`
  Section 6 does; the convention of `docs/40` Section 9 (MCSE at the
  observed rate) is retired (resolves CR-26).
- **MCSE with every number.** Every rejection rate carries its
  binomial MCSE $\sqrt{\hat\pi(1-\hat\pi)/n_{\text{reps}}}$; every
  pooled rate its pooled MCSE; every difference between independent
  runs $\sqrt{\text{MCSE}_1^2 + \text{MCSE}_2^2}$
  (`RERUN-anchored.md` Section 7); paired comparisons within a run
  use McNemar.
- **Warp-speed bootstrap.** Rates from one resample per trial
  (Giacomini, Politis and White, 2013) are labeled warp-speed; their
  MCSE is computed by resampling trials and is larger than the
  binomial one (0.0097 to 0.0112 against 0.0069 per cell at
  $n_{\text{sim}} = 1{,}000$; quoted, `docs/57` Section 8.1). The
  binomial MCSE is never quoted alone for a warp-speed rate.
- **Size-adjusted power** follows CM "Size-adjusted power"; maximin
  shortfall states its cell set in the caption (non-prognostic
  alternative cells, or all 26 alternative cells), so the two
  readings of `docs/43` and `docs/47` Section 7 stay apart (CR-29).
- **Replicate bands.** Null cells $n_{\text{reps}} \ge 1{,}000$ and
  alternative cells $\ge 500$ for any number that enters a decision;
  smaller runs are labeled "illustration" or "pilot" in the caption.
- **Sample size.** $N$ is always the total (NOTATION rule 4), and
  every study table states the allocation (fixed or randomized).

---

## Analysis-ID scheme

<!-- From docs/59-whitepaper-refactoring-plan.md Section 2.7 -->

The stored interaction-track labels (S, C1, C2, A, ES, EA, F1, F2,
R, DS, DA, G1, G2, F1A, F1A0, F1A0B, F1AB, MA, MaxT and their
suffixes) are keys in the stored stress-test files, and NOTATION
Part 2 forbids a display-time mapping layer. The scheme (decided
2026-10-11, DQ-5) therefore keeps them and changes what is new or
local. A document names the track at first use ("interaction analysis
F1", "main-effect analysis `M-F1`").

| Class | Form | Examples | Change |
|----------------------|----------|--------------------|----------------------|
| Interaction track (stress test, covariance study selection, data-chosen half-life) | bare label | S, F1, F1A0B, G2 | none; the default track of Parts III and IV |
| Main-effect track (P04, `docs/58`) | `M-` prefix | `M-S`, `M-F1`, `M-VF`, `M-RV`, `M-F1A` | rename in P04, `docs/58` content and `ch24` before paper 04's stored data exist (DQ-5, decided) |
| Covariance study analyses (`docs/45`, `docs/43` Sections 3.7 and 3.8) | `CV-` prefix | `CV-V1` to `CV-V7`, `CV-F1` to `CV-F7` | removes the collision of catalog F1 with stress-test F1 and configuration F |
| Paper-local catalog IDs | paper prefix, as NOTATION's "05-S2" rule | `02-G1` to `02-G9`, `06-C1` to `06-C4`, `08-T1`, `10-M0` | removes collisions with G1, G2, C1, C2 and configuration C |
| Published analyses (`docs/43` Section 3.1) | `PUB-` prefix | `PUB-O1` to `PUB-O5` | none in substance |
| Summary measures and time adjustments (`docs/43` S1 to S3, D1 to D7) | `SM-`, `TA-` prefix | `SM-1`, `TA-3` | removes collisions with sensitivity sweeps S1, S2 and with $D_{it}$ |
| Exposure codings | E1 to E9, E9b, E`n`cr2 | unchanged | NOTATION extended E series |
| Constructs | "configuration A" to "configuration F", always spelled out, with $\rho$ | "configuration C at $\rho = 0.7$" | never a bare letter |

Table: The analysis-ID scheme for the series

The registry of all IDs, with formula, track, source and stored-data
key, is a machine-readable table (`docs/series/registry/analyses.csv`)
from which Appendix A is generated.

---

## Canonical names

<!-- From docs/59-whitepaper-refactoring-plan.md Section 2.8 -->

- **The published process**: Hendrickson's `generateData()` and
  `buildtrialdesign()` at commit `58b32a9` (NOTATION Part 3). Never
  "the original code", "orig" or `42ac030` for the published process;
  `42ac030` is the initial commit, whose `generateData.R` equals
  `58b32a9`'s, but analysis defaults may differ.
- **The 2024 revision**: commits `8609f12`, `f6ee86d`, `f70f86d`,
  `325314f` (scale factor $s_f$, ratio coupling, timepoint-1 guard,
  `Dbc` in the analysis). **The vendored copy**: the 2024 state
  copied as `06dac83`, which `docs/54` calls `orig`.
- **The refined process**: as NOTATION Part 3; its validated
  fixed-half-life analysis uses $t_{1/2}^a = 0.2$ week.
- **Configurations A to F**: as NOTATION Part 3, always with $\rho$.
  `covar` is replaced in prose by "configuration F at the package's
  parameters"; `orig` by "the published construct" or "the vendored
  copy", whichever is meant.
- **The published construct**: configuration A only (compound
  symmetry with $\rho_{cs} = 0.8$, constant cross block). AR(1)
  variants of the `04-power-simulation.R` construct are named by
  configuration (resolves CR-30).
- **Architectures**: mean-moderation, covariance-moderation,
  dual-channel. **Designs**: OL, CO, OL+BDC, Hybrid (capitalized).
- **The estimand**: the biomarker-treatment interaction $\beta_{bm:D}$
  (`docs/30` D3); "biomarker-by-drug" in the titles of `docs/43` to
  `docs/57` becomes "biomarker-treatment" (DQ-7, decided).
- **Participant**, never patient; **visit**, never occasion;
  **off-drug period**, never washout for the span after
  discontinuation.

---

## Mathematical symbols

<!-- From analysis/report/NOTATION.md Part 1 -->

Each table gives the symbol, its meaning, its units ("none" for a
dimensionless quantity, "points" for outcome-scale points), the code
identifier where one exists, and where it is used. "CM" is
`analysis/report/compendium-mathematics.tex`; "docs/NN" is a
whitepaper under `docs/`; "all" means the manuscripts and whitepapers
generally. Entries marked (2026-10-09) were added or changed in the
merge described above.

### Indices and structure

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $i$ | participant index | none | `ptID` | all |
| $t$ | visit (timepoint) index; weekly measurement occasions throughout | none | row order within `ptID` | all |
| $j$, $ij$ | period or occasion index within participant | none | | papers 04, 06 |
| $w_t$, $w_s$ | calendar week of visit $t$ (Hybrid: 4, 8, 9, 10, 11, 12, 16, 20); reserved for calendar time, never a weight (2026-10-09) | weeks | column `t` in drivers, `t_wk` in stored data | docs/34 to docs/58, CM |
| $d_t$ | gap between visits, $w_t - w_{t-1}$ (2026-10-09) | weeks | | docs/35, docs/36, CM |
| $N$ | **total** number of participants in the trial, across all randomization paths; never a per-path count. See Part 2, sample size and randomization paths | participants | `N` (caller-dependent, Part 2) | all |
| $P$ | number of randomization paths in a design; the expected allocation to each is $N/P$ | none | `length(td$trialpaths)` | all |
| $q$, $p$ | randomization-path index, $q = 1, \ldots, P$ (2026-10-09) | none | `path` | docs/40, docs/57; docs/37, docs/58, CM |
| $\varpi_p$ | share of participants allocated to path $p$, $n_p/N$ (2026-10-09; formerly $q_p$) | none | | docs/37, docs/58, paper 04 protocol, CM |
| $n_{\text{on}}$, $n_{\text{off}}$, $n_v$ | numbers of on-drug and off-drug post-baseline visits (rows) of a path, and of rows of visit kind $v$ (2026-10-09) | rows | | docs/34, docs/37 to docs/39, docs/49, docs/55, docs/56 |
| $n_c$ | number of response components (3; 2 in papers 13 and 14) (2026-10-09) | none | | CM |
| $n_{\text{reps}}$, $n_{\text{sim}}$ | Monte Carlo replicates per cell | replicates | `n_reps` | all |
| $k$ | cycles per participant in design sweeps; also the Weibull shape parameter in the decay family | none | | papers 05, 12; CM |

### Outcome and latent components

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $Y_{it}$ | observed symptom score; the modeled outcome | points | column `Sx` | all |
| $\mathrm{BL}_i$ | participant-specific baseline level | points | | all |
| $BR_{it}$ | biological (pharmacological) response component | points | | all |
| $PB_{it}$ | placebo-belief response component | points | | all |
| $TV_{it}$ | time-variant natural-history component | points | | all |
| $\varepsilon_{it}$ | observation-level residual | points | | all |
| $u_i$ | participant random intercept | points | `random = ~1` by `ptID` | all |
| $\sigma_u^2$, $\sigma_\varepsilon^2$ | variances of the random intercept and of the residual of the analysis model (2026-10-09) | points$^2$ | | docs/41, docs/56, docs/57, CM |
| $\eta_t$ | design expectancy at visit $t$: 1 open label, 0.5 blinded, 0 at baseline (2026-10-09; formerly $e_t$) | none | column `e`; regressor `De` | papers 01, 06; docs/36 to docs/58, CM |
| $\tilde\eta_{it}$ | measured (participant-reported) placebo-belief covariate, if an analysis uses one; reserved (2026-10-09) | none | | docs/43 |
| $\sigma_{TV}$, $\sigma_{PB}$, $\sigma_{BL}$ | standard deviations of the natural-course and placebo-belief components and of the baseline level (2026-10-09) | points | `blparam`, `respparam` SDs | docs/35 to docs/38, docs/46, docs/49, docs/51, CM |

### Treatment exposure and carryover

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $D_{it}$ | binary drug state: 1 on drug, 0 off drug | none | `ondrug` | all |
| $D_{bc,it}$ | continuous exposure-decayed drug indicator: 1 on drug, $e^{-\lambda t_{sd}}$ off drug after exposure, 0 at off-drug visits before the first exposure (2026-10-09) | none | column `Dbc` | all |
| $t_{sd}$ | time since discontinuation; input to the carryover decay | weeks | `tsd` | all |
| $t_{od}$ | cumulative time on drug; input to the Gompertz response | weeks | `tod` | all |
| $t_{pb}$ | time-by-placebo-belief accumulation | weeks | `tpb` | all |
| $t_{1/2}$ | carryover half-life of the data-generating process | weeks | `carryover_t1half` | all |
| $t_{1/2}^a$ | half-life assumed by the analysis when it builds $D_{bc,it}$ (2026-10-09) | weeks | `t_half_a` | docs/43, docs/45 to docs/47, docs/56 to docs/58, CM |
| $\lambda$ | carryover decay rate, $\ln 2 / t_{1/2}$ | per week | `lambda_cor` | all |
| $\lambda_w$ | Weibull rate, $(\ln 2)^{1/k} / t_{1/2}$ | per week | | CM |
| $\phi(\cdot)$, $\phi_t$ | decay-family weight and its value at visit $t$ (2026-10-09) | none | `carryover_decay()` | docs/34, docs/36, docs/37, CM |
| $G(t_{od}; m, d, r)$ | modified Gompertz mean of a component, $G(0) = 0$ (2026-10-09) | points | `modgompertz` | docs/40, docs/46, docs/49, docs/51, docs/53, docs/55, CM |
| $\mu_{BR,t}$, $\mu_{\text{last}}$ | mean drug response at visit $t$ after carryover; its value at the last on-drug visit, the anchor of anchored carryover (2026-10-09) | points | `brmeans` | docs/36, docs/37, docs/53 to docs/55 |
| $s_f$ | scale factor on the carryover exponent in the 2024 upstream code, default 2 (2026-10-09) | none | `scalefactor` | docs/53 to docs/55 |

Decay families implemented on the data-generating side: exponential
$e^{-\lambda t}$; linear $\max(0, 1 - t/(2t_{1/2}))$; Weibull
$e^{-(\lambda_w t)^k}$; power $\max(0, (1 - t/(3t_{1/2}))^p)$.

### Biomarker and moderation

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $B_i$ | participant-level pre-treatment biomarker value | biomarker units | `bm` | all |
| $b_i$ | standardized biomarker, $(B_i - \bar{B})/s_B$ | none | `bm_z` | all |
| $B^c_i$ | biomarker centered at the trial's baseline sample mean, $B_i - \bar B$, not scaled (2026-10-09) | biomarker units | `bmc` | docs/47, docs/56, docs/57, CM |
| $\mu_{bm}$, $\sigma_{bm}$ | population mean and standard deviation of the biomarker | biomarker units | `blparam` | docs/36 to docs/38, docs/48, CM |
| $\sigma_{BR}$ | standard deviation of the drug-response component | points | `respparam` | all |
| $c_{bm}$ | moderation parameter of the **covariance-moderation** architecture: on-drug correlation between biomarker and $BR$. A correlation, bounded in $[-1,1]$ | none | `c.bm` | all |
| $\beta_{bm}$ | moderation parameter of the **mean-moderation** architecture: dimensionless multiplier scaling the treatment effect by $b_i$ | none | `c.bm` with `dgp_architecture = 'mean_moderation'` | all |
| $c_{bm,A}$, $c_{bm,B}$ | independent mean-channel and covariance-channel weights of the dual-channel architecture | none | `c.bm_a`, `c.bm_b` | paper 11 |
| $g_t$, $\mathbf g$ | coupling profile: the fraction of $c_{bm}$ present at visit $t$ (2026-10-09) | none | built in `buildSigma()` | papers 01, 11; docs/34 to docs/37, docs/40, docs/41, docs/44, docs/49, docs/51, docs/55, CM |
| $c_{bm}^{\ast}$ | coupling ceiling: the largest $c_{bm}$ for which every path's joint correlation matrix is positive definite (2026-10-09) | none | | docs/34 to docs/37, docs/39, docs/44, docs/54, CM |
| $c_{bm,TV}$, $c_{bm,PB}$, $c_{bm,BL}$ | prognostic couplings: correlation of the biomarker with $TV$, $PB$, and the baseline level (2026-10-09; $c_{bl}$ retired) | none | `c_tv`, `c_pb`, `c_bl` | docs/35, docs/44, docs/46 to docs/49, docs/51, docs/56, docs/57, CM |
| $s_{\text{prog}}$, $f$ | prognostic slope and the share absorbed by a shared biomarker main effect (2026-10-09) | points per biomarker unit; none | | docs/48, docs/49, docs/51, docs/56, CM |

### Estimand, inference, and performance

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $\beta_{bm:D}$ | the target estimand: the biomarker-by-treatment interaction coefficient | points per biomarker unit | `bm:Dbc` | all |
| $\beta_D$ | treatment main effect | points | coefficient of `Dbc` or `D` | papers 04, 05; docs/56 to docs/58, CM |
| $\tau^2$ | between-participant variance of the drug effect: treatment heterogeneity (Senn); in the paper 04 protocol the variance of $v_i$ (2026-10-09) | points$^2$ | `tau_squared`, `tau` | paper 04, its protocol, docs/58, CM |
| $\kappa$ | ratio of mean model-based standard error to empirical standard deviation; $\kappa > 1$ means a conservative test | none | | papers 08, 10 |
| $\nu$ | denominator degrees of freedom | none | | all |
| $\alpha$ | nominal significance level; 0.05 throughout | none | | all |
| $\pi$, $\hat\pi$ | power and an estimated rejection rate, with $\mathrm{MCSE}(\hat\pi) = \sqrt{\hat\pi(1-\hat\pi)/n_{\text{reps}}}$ (2026-10-09) | none | `power`, `mc_se_power` | all |

### Covariance construct

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $\rho$ | AR(1) within-component correlation per week (refined process 0.7; configuration C in `docs/36` 0.8) | per week | `rho` | all |
| $\rho_{cs}$ | common within-component correlation of a compound-symmetry construct (published process 0.8) (2026-10-09) | none | `rho` in the published code | docs/34 to docs/36, docs/44, docs/48, docs/51, docs/54, CM |
| $c_1$, $c_\times$ | cross-component correlation at the same visit and at different visits | none | `c.cf1t`, `c.cfct` | docs/34 to docs/37, docs/40, docs/44, docs/50, docs/51, docs/54, CM |
| $K$, $A$, $M$ | $n_c \times n_c$ cross-component correlation, temporal correlation ($A_{ts} = \rho^{|w_t - w_s|}$ under AR(1)), and the response block | none | | docs/34 to docs/37, docs/39, docs/44, docs/51, CM |

**Reserved letters (2026-10-09).** Besides the letters of the tables
above, these may not be reused for another meaning in a document that
also uses the canonical one: $w$ (calendar week), $f$ (absorbed
share), $\tau$ (heterogeneity variance), $\eta$ (design expectancy),
$\rho$ (AR(1) correlation), $\phi$ (decay family) against $\varphi$
(normal density), $q$ (path index), $g$ (coupling profile),
$\Delta$ (E9 contrast), and $m$, $d$, $r$ when they appear
undecorated beside a Gompertz curve.

### The five rules the symbols encode

1. **$D_{it}$ is binary; $D_{bc,it}$ is continuous.** These are
   different objects, and the distinction is the subject of the
   carryover-specification manuscript. Never use one for the other.
2. **$c_{bm}$ is a correlation; $\beta_{bm}$ is a multiplier.** They
   are calibrated to a common numeric scale, because the mean-channel
   shift is applied as $\beta_{bm}\,b_i\,\sigma_{BR}$ on the
   standardized biomarker, so the reference value $0.45$ denotes a
   matched moderation strength under either architecture. The symbols
   are nonetheless not interchangeable.
3. **Mathematics in model definitions, code identifiers in code
   contexts.** Write $Y_{ij} = \beta_0 + \beta_D D_{bc,ij} + \ldots$
   for the model; write `Sx ~ bm + t + Dbc + bm:Dbc` in `\texttt{}`
   when quoting the R formula. Never mix the two inside one expression.
4. **$N$ is the total sample size.** $N$ is the total number of
   participants across all randomization paths, and the expected
   allocation to each of $P$ paths is $N/P$: with $P = 4$ and
   $N = 140$, each path receives about 35. Never use $N$ for a
   per-path count.
5. **State the sign convention.** The three components are
   non-negative *reductions* in symptom severity, so an increase in
   any component lowers $Y$. Treatment effects and interaction
   coefficients are consequently negative; moderation parameters are
   positive.

---

## Code identifiers and stored data values

<!-- From analysis/report/NOTATION.md Part 2 -->

Where a stored value is also a reporting label, the two are kept
**identical**, so that no display-time mapping stands between the data
and the page.

### Exported package functions

`buildtrialdesign`, `buildSigma`, `generateData`, `censordata`,
`lme_analysis`, `generateSimulatedResults`, `validateParameterGrid`,
`cumulative`, `modgompertz`, `trajectoryShape`, `shape_logistic`,
`shape_hyperbolic_tangent`, `shape_piecewise_linear_breakpoint`,
`characterize_carryover`, `analyze_trial_extended`,
`print_carryover_summary`, `print_trial_summary`,
`plotfactortrajectories`, `PlotModelingResults`, `reknitsimresults`.

### Principal arguments

| Identifier | Role |
|---|---|
| `modelparam` | list: sample size `N`, moderation `c.bm` |
| `respparam` | Gompertz parameters `maxr`, `rate`, `disp` per component |
| `blparam` | biomarker and nuisance means and standard deviations |
| `trialdesign` | output of `buildtrialdesign`: `timepoints`, `timeptnames`, `expectancies`, `ondrug` |
| `dgp_architecture` | `'mvn'`, `'mean_moderation'`, or `'combined'` |
| `carryover_t1half` | carryover half-life |
| `lambda_cor` | decay rate; derived from the half-life when `NA` |
| `cached_sigma` | pre-built covariance matrix, reused across parallel workers |
| `makePositiveDefinite` | positive-definiteness repair switch |
| `n_cores` | worker count; `-1` auto-detects |
| `save_chunks` | progressive checkpointing for long runs |

### Randomization-path counts by design

| Design | $P$ | Paths |
|---|---|---|
| OL | 1 | on drug throughout |
| CO | 2 | drug-then-placebo, placebo-then-drug |
| OL+BDC | 2 | discontinuation at either of two visits |
| Hybrid | 4 | two discontinuation points crossed with two crossover orders |

### Data columns

`ptID` (participant), `Sx` (outcome, the $Y_{it}$ above), `bm`
(biomarker), `bm_z` (standardized biomarker), `Dbc` (continuous
exposure indicator, $D_{bc,it}$), `tsd`, `tod`, `tpb`, `t_wk` (time
in weeks), `path` (randomization path).

### Stored factor levels and their reporting labels

| Column | Stored values | Reported as | Mapped? |
|---|---|---|---|
| `spec` | `E1`, `E2`, `E3` | E1, E2, E3 | no, identical |
| `carryover_form` | `exponential`, `linear`, `weibull` | as stored | no, identical |
| `design` | `CO`, `Hybrid`, `OLBDC` | CO, Hybrid, OL+BDC | only `OLBDC` to `OL+BDC` |
| `dgp_arch` | `mean_moderation`, `mvn`, `combined` | Mean, Covariance, Dual | yes; these are package API values |

### Summary-table columns

`n_reps`, `n_converged`, `converged_frac`, `non_convergence_rate`,
`power`, `bias`, `mean_estimate`, `true_value`, `empirical_se`,
`mse`, `coverage`, and an `mc_se_`-prefixed Monte Carlo standard
error for each.

---

## Controlled vocabulary and glossary

<!-- From analysis/report/NOTATION.md Part 3 -->

The project lexicon at `docs/29-nof1-precision-medicine-lexicon.md`
holds the full terminology set (241 terms across 13 thematic
sections, with cross-field synonyms). The entries below are the
subset carrying a fixed, enforced meaning in this compendium.

### Trial designs

**OL** (open-label). Unblinded titration; no off-drug contrast, so it
supplies no within-subject interaction information.
**BDC** (blinded discontinuation). A blinded switch to placebo after
open-label response.
**OL+BDC**. Open-label titration followed by blinded discontinuation.
**CO** (crossover). Traditional two-period within-subject crossover.
**Hybrid**. The Hendrickson (2020) design: open-label titration,
blinded discontinuation, then a brief crossover. The program's
reference design. Capitalize it.
**Aggregated N-of-1 trial**. Many single-patient on/off series pooled
into one mixed model.

### Data-generating architectures

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

### Named data-generating processes (added 2026-10-08)

The current best DGP is called **refined** (user decision,
2026-10-08): "the refined process" in prose, "Refined" in legends and
tables, against "the published process" / "Published" for
Hendrickson's code.

| Name in prose | Meaning | Defined in |
|---|---|---|
| the published process | Hendrickson's `generateData()` and `buildtrialdesign()` at commit `58b32a9`: compound symmetry within components ($\rho_{cs} = 0.8$), constant cross-component correlation, step coupling, compounding carryover recursion, time on drug restarting at zero on re-exposure | Paper 01 |
| the refined process | The recommended DGP: graded coupling (0 before first exposure), separable AR(1) with $\rho = 0.7$ and $c_1 = 0.2$, anchored carryover, a stop on matrices that are not positive definite, randomized path allocation; re-exposure as in the published process; optional prognostic couplings, phase-specific measurement error, Weibull decay and dropout | Paper 01, Table `tab:spec` |
| the published construct | The covariance construct of the published process alone (configuration A: compound symmetry within components, $\rho_{cs} = 0.8$, constant cross-component correlation), with whatever coupling the context states; not the whole published process (added 2026-10-09) | CM, glossary |
| configuration A to F | The covariance and coupling constructs compared in `docs/36` (A published; B AR(1) with the `covar` cross block, step coupling; C separable AR(1) with graded coupling; D the published compound-symmetry construct with mean moderation; E separable AR(1) with mean moderation and no covariance coupling; F AR(1) with the `covar` cross block, graded coupling) | `docs/36`, Section 6 and Appendix A.6 |

Use **the refined process** for the DGP of the stress test and the
recommendation. Use **configuration C** only for the covariance and
coupling construct itself, as in the `docs/36` comparison, and state
$\rho$ with it (`docs/36` evaluates it at $\rho = 0.8$; the refined
process uses $\rho = 0.7$). Configuration C is not the refined
process: it does not by itself include anchored carryover, the
positive-definite stop or the options.

### Analysis specifications for carryover

**E1.** Binary on-drug indicator; carryover ignored.
**E2.** Exposure-weighted continuous predictor $D_{bc}$, committing to
an assumed decay form and half-life.
**E3.** Binary indicator plus a lagged just-off-drug nuisance term;
the classical Jones-Kenward crossover device, assuming no decay form.

The letter E stands for the exposure regressor that distinguishes
them. These were labeled `A1`/`A2`/`A3` until 2026-07-31; the stored
data and the drivers were migrated to E1/E2/E3, so the labels now
agree everywhere and nothing is mapped at display time.

**The extended E series (2026-10-09):**

| Label | Meaning | Defined in |
|---|---|---|
| E4 | binary drug state plus linear time since discontinuation (retired) | `docs/43` |
| E5 | lag term crossed with the biomarker, $L_{it}$ and $B_i L_{it}$ (retired) | `docs/43`, paper 02 |
| E6 | $t_{sd}$ crossed with the biomarker (retired) | `docs/43`, paper 02 |
| E7 | E2 form with the half-life chosen by AIC | paper 02, `docs/43`, CM |
| E9 | paired-difference statistic: each participant's on-drug mean minus off-drug mean, $\Delta_i$, regressed on the biomarker by OLS | paper 01, paper 02, `docs/37`, CM |
| E9b | E9 with baseline counted as one more off-drug visit | paper 01, `docs/37`, CM |
| E1cr2, E2cr2, E3cr2, E7cr2 | the same formulas with the CR2 cluster-robust standard error | paper 02, CM |

### Response-curve families

**Modified Gompertz** (default), $y = m\,\exp(-d\,\exp(-r\,t))$ with
a vertical offset so the curve passes through the origin; written
$G(t_{od}; m, d, r)$ for a component's mean, with ceiling $m$
(points), displacement $d$ (none) and rate $r$ (per week), and
inflection time $t^\ast = \ln(d)/r$ (weeks). The code names are
`maxr`, `disp` and `rate` (`respparam`); write the mathematical
symbols in equations and the code names only in code (2026-10-09).
Alternatives: **logistic** (symmetric sigmoidal), **hyperbolic
tangent** (fast-saturating), **piecewise-linear breakpoint** (no
smoothness assumption).

### Missing data and dropout

**MCAR**, **MAR**, **MNAR** in the standard Rubin sense.
**Informative dropout.** Hazard depending on cumulative symptom
worsening since baseline.
**Happy accident.** The randomization-path selection effect by which
informative dropout preferentially preserves the most informative
crossover blocks.

### Inference

**LME.** Linear mixed-effects model; the program's default analysis.
**corCAR1.** Continuous-time AR(1) residual correlation structure.
**GEE.** Generalized estimating equations; the marginal comparator.
**Mancl-DeRouen**, **CR2.** Small-sample bias-corrected sandwich
variance estimators.
**RM-ANOVA.** Repeated-measures analysis of variance.
**Anti-conservative.** Realized Type I error above nominal.
**Working covariance.** The assumed residual covariance of the LME.

### Estimands and reporting

**Biomarker-treatment interaction.** The program's core estimand;
equivalently the patient-by-treatment interaction, or heterogeneity
of treatment effects (HTE) in the precision-medicine literature.
**Predictive biomarker.** Predicts differential benefit, as opposed
to a **prognostic** biomarker, which predicts outcome regardless of
treatment.
**ADEMP.** The Morris, White and Crowther (2019) reporting framework:
aims, data-generating mechanisms, estimands, methods, performance
measures.
**MCSE.** Monte Carlo standard error; reported with every performance
measure.

### Enumeration conventions

**Simulation studies.** Number them (Study 1, Study 2, ...) rather
than lettering them.
**Sensitivity sweeps.** Labels S1, S2, ... are paper-local and are
not comparable across papers. Cite them with the paper prefix when
crossing manuscripts (05-S2, not S2).

### Prose terminology

Symbol conventions are settled above; prose terminology is settled by
`docs/30-terminology-consistency-plan.md`, executed 2026-07-30 with
its D1 spelling decision locked to **US English**. The two are
complementary: this chapter governs mathematics and identifiers, that
document governs surface forms and acronym expansion.

## Monte Carlo simulation methodology

Each simulation study is a Monte Carlo experiment that generates
$n_{\text{reps}}$ independent trial replications per parameter cell
under a fully specified data-generating process (DGP) and analysis
pipeline. A cell is defined by a unique combination of design
parameters (trial design, sample size $N$, moderation strength
$c_{bm}$, carryover half-life $t_{1/2}$, DGP architecture) and
analysis method. Cells are independent: no data are shared across
cells.

Within each replication, the pipeline proceeds in three steps.
First, one synthetic trial's data are generated by `generateData()`,
which draws participant outcomes from the multivariate normal DGP
and applies the specified DGP architecture (covariance moderation or
mean moderation). Second, the analysis model is fitted to those data
by `lme_analysis()`, which estimates the interaction term
$\hat\beta_{bm:D}$ and its standard error via the specified analysis
method. Third, the result is recorded: whether
$\hat\beta_{bm:D} / \widehat{\mathrm{SE}}$ exceeds the critical
value at $\alpha = 0.05$ (two-sided), and the point estimate,
model-based SE, and any requested ancillary statistics.

Power is the empirical rejection rate
$\hat\pi = R / n_{\text{reps}}$, where $R$ is the count of
replications that reject the null hypothesis
$H_0 : \beta_{bm:D} = 0$ under the alternative DGP. The empirical
size is $\hat\pi$ when the DGP is generated under the null
($c_{bm} = 0$). The Monte Carlo standard error (MCSE) of $\hat\pi$
is the binomial standard error
$\sqrt{\hat\pi(1 - \hat\pi) / n_{\text{reps}}}$.

Parameter sweeps dispatch cells to parallel workers via
`generateSimulatedResults()`, which uses `furrr::future_map` with
independent random-number streams per cell. All cells within a run
are computed in the same session to ensure consistent
platform-level numerical behaviour. Cell independence means the
order of computation does not affect results.

Reproducibility relies on per-cell seeds recorded in the parameter
grid alongside each cell's output. The same cell, DGP, and analysis
method will produce draw-for-draw identical results when run from
the same platform with the same seed. Results generated on different
platforms (different BLAS/LAPACK libraries, different operating
systems) are not draw-for-draw identical even under the same seed,
because floating-point operations in the linear-algebra back-end
may differ. Arms that must be paired trial-by-trial -- for example,
paired power comparisons that use McNemar's test -- must therefore
be run on the same platform.

The target $n_{\text{reps}}$ for main production runs is 500 per
cell, yielding an MCSE of approximately 0.010 at $\hat\pi = 0.5$
and 0.022 at the nominal size $\alpha = 0.05$; the 2-MCSE band at
the nominal size is 0.006 to 0.094 for 500 replications. Rerun
cells -- cells where a methodological change requires an updated
estimate -- use 1,000 replications, halving the MCSE to 0.007 at
$\hat\pi = 0.5$.

<!-- TODO: [project glossary: cross-field synonyms table from docs/29;
cross-reference to appE] -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | all | Phase 2 content move: NOTATION.md Parts 1-3; evidence labels, MC conventions, analysis-ID scheme and canonical names from docs/59 Sections 2.5-2.8 | Phase 2 skeleton fill |
