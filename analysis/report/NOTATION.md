# Canonical notation, identifiers, and glossary

*Originally written 2026-07-30 11:57 PDT; expanded to a full reference
2026-07-30 20:20 PDT; revised 2026-08-06 to settle $N$ as the total
across randomization paths; revised 2026-10-09 to merge the symbols
of the whitepapers `docs/34` to `docs/58` and the mathematics
compendium (see "Revision 2026-10-09" below)*

This file is the single source of truth for mathematical notation,
code identifiers, and controlled vocabulary across the manuscripts
under `analysis/report/`. New manuscripts should adopt it rather than
re-deriving symbols, and any deviation should be stated explicitly in
the deviating paper's Methods section.

`analysis/report/whitepaper-notation-audit.md` was intended to hold a
cross-reference audit of this symbol table against each manuscript's
actual usage, but as of 2026-08-20 that file contains only YAML front
matter and no audit content (pub_review whitepaper, Major Issue 2.5).
Until that audit is written, this file's consistency with the eleven
manuscripts rests on the compendium-wide `tools/notation-lint.pl`
mechanized check (Minor Issue 3.9) rather than on a written audit
trail.

Symbols and identifiers listed here were extracted from the manuscript
sources, from `R/`, from the drivers under `analysis/scripts/`, and
from the archived `.rds` summaries.

## Revision 2026-10-09

The notation of the whitepapers `docs/34` to `docs/58` and of
`analysis/report/compendium-mathematics.tex` was reviewed in five
passes (2026-10-09) and the proposed symbols were merged here. Seven
user decisions of 2026-10-09 settle the collisions the passes found:

1. $\tau^2$ is the between-participant variance of the drug effect
   (treatment heterogeneity, as in Senn, paper 04, the paper 04
   working protocol and `docs/58`). The residual variance of the E9
   contrast given the biomarker, previously written $\tau^2$ and
   $\tau_p^2$ in `docs/37` to `docs/39`, the compendium and paper 01,
   is now $\sigma_\Delta^2$ and $\sigma_{\Delta,p}^2$.
2. $\eta_t$ is the design expectancy (1 open label, 0.5 blinded, 0 at
   baseline; code `e`, regressor `De`). A measured placebo-belief
   covariate, if one is used, is $\tilde\eta_{it}$. Paper 06's
   $\eta(\text{phase})$ is the design expectancy, not a measured
   quantity (paper 06 assigns it by phase, calls the blinded value "a
   modeling assumption, not an estimate", and recommends that trials
   collect an expectancy measure), so the former reading of
   $\eta$ as "measured belief, paper 06" was withdrawn.
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
7. The prognostic couplings are $c_{bm,TV}$, $c_{bm,PB}$, $c_{bm,BL}$
   ($c_{bl}$ is retired).

Adopted without a decision, by consensus of the passes: the normal law
is $\mathcal N(\cdot,\cdot)$, never $N(\cdot,\cdot)$; $\varphi$ and
$\Phi$ are the standard normal density and distribution function;
$\mathbb E$ is the expectation operator; the remaining entries marked
"(2026-10-09)" in the tables below. Where two documents still
disagreed, the canonical meaning here won, and otherwise the
compendium's choice won.

# Part 1. Mathematical symbols

Each table gives the symbol, its meaning, its units ("none" for a
dimensionless quantity, "points" for outcome-scale points), the code
identifier where one exists, and where it is used. "CM" is
`analysis/report/compendium-mathematics.tex`; "docs/NN" is a
whitepaper under `docs/`; "all" means the manuscripts and whitepapers
generally. Entries marked (2026-10-09) were added or changed in the
merge described above.

## Indices and structure

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $i$ | participant index | none | `ptID` | all |
| $t$ | visit (timepoint) index; weekly measurement occasions throughout | none | row order within `ptID` | all |
| $j$, $ij$ | period or occasion index within participant | none | | papers 04, 06 |
| $w_t$, $w_s$ | calendar week of visit $t$ (Hybrid: 4, 8, 9, 10, 11, 12, 16, 20); reserved for calendar time, never a weight (2026-10-09) | weeks | column `t` in drivers, `t_wk` in stored data | docs/34 to docs/58, CM |
| $d_t$ | gap between visits, $w_t - w_{t-1}$ (2026-10-09) | weeks | | docs/35, docs/36, CM |
| $N$ | **total** number of participants in the trial, across all randomization paths; never a per-path count. See Part 2, sample size and randomization paths | participants | `N` (caller-dependent, Part 2) | all |
| $P$ | number of randomization paths in a design; the expected allocation to each is $N/P$ | none | `length(td$trialpaths)` | all |
| $q$, $p$ | randomization-path index, $q = 1, \ldots, P$, where a document also uses $p$ for a $p$-value or an exponent (2026-10-09). The subscript $p$ on path-specific quantities fixed by the 2026-10-09 decisions and the compendium ($\sigma_{\Delta,p}^2$, $\tilde\beta_p$, the path share $\varpi_p$) is retained; a document uses one of the two letters for paths throughout | none | `path` | $q$: docs/40, docs/57; $p$: docs/37, docs/58, CM |
| $\varpi_p$ | share of participants allocated to path $p$, $n_p/N$ ($1/P$ under equal allocation) (2026-10-09; formerly $q_p$ in `docs/37` and the compendium) | none | | docs/37, docs/58, paper 04 protocol, CM |
| $n_{\text{on}}$, $n_{\text{off}}$, $n_v$ | numbers of on-drug and off-drug post-baseline visits (rows) of a path, and of rows of visit kind $v$ (2026-10-09) | rows | | docs/34, docs/37 to docs/39, docs/49, docs/55, docs/56 |
| $n_c$ | number of response components (3; 2 in papers 13 and 14) (2026-10-09) | none | | CM |
| $n_{\text{reps}}$, $n_{\text{sim}}$ | Monte Carlo replicates per cell | replicates | `n_reps` | all |
| $k$ | cycles per participant in design sweeps; also the Weibull shape parameter in the decay family | none | | papers 05, 12; CM |

## Outcome and latent components

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $Y_{it}$ | observed symptom score; the modeled outcome | points | column `Sx` | all |
| $\mathrm{BL}_i$ | participant-specific baseline level | points | | all |
| $BR_{it}$ | biological (pharmacological) response component | points | | all |
| $PB_{it}$ | placebo-belief response component | points | | all |
| $TV_{it}$ | time-variant natural-history component | points | | all |
| $\varepsilon_{it}$ | observation-level residual | points | | all |
| $u_i$ | participant random intercept | points | `random = ~1` by `ptID` | all |
| $\sigma_u^2$, $\sigma_\varepsilon^2$ | variances of the random intercept and of the residual of the analysis model (2026-10-09; not $\sigma_b^2$, $\sigma_e^2$) | points$^2$ | | docs/41, docs/56, docs/57, CM |
| $\eta_t$ | design expectancy at visit $t$: 1 open label, 0.5 blinded, 0 at baseline; paper 06 writes $\eta(\text{phase})$ for the same object (2026-10-09; formerly $e_t$) | none | column `e`; regressor `De` | papers 01, 06; docs/36 to docs/58, CM |
| $\tilde\eta_{it}$ | measured (participant-reported) placebo-belief covariate, if an analysis uses one; reserved (2026-10-09) | none | | docs/43 |
| $\sigma_{TV}$, $\sigma_{PB}$, $\sigma_{BL}$ | standard deviations of the natural-course and placebo-belief components and of the baseline level; $\sigma_{PB,t} = \sigma_{PB}\,\eta_t$ (2026-10-09) | points | `blparam`, `respparam` SDs | docs/35 to docs/38, docs/46, docs/49, docs/51, CM |

## Treatment exposure and carryover

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $D_{it}$ | binary drug state: 1 on drug, 0 off drug | none | `ondrug` | all |
| $D_{bc,it}$ | continuous exposure-decayed drug indicator: 1 on drug, $e^{-\lambda t_{sd}}$ off drug after exposure, 0 at off-drug visits before the first exposure; in an analysis it is built at the assumed half-life $t_{1/2}^a$ (2026-10-09: pre-exposure value stated) | none | column `Dbc` | all |
| $t_{sd}$ | time since discontinuation; input to the carryover decay | weeks | `tsd` | all |
| $t_{od}$ | cumulative time on drug; input to the Gompertz response | weeks | `tod` | all |
| $t_{pb}$ | time-by-placebo-belief accumulation | weeks | `tpb` | all |
| $t_{1/2}$ | carryover half-life of the data-generating process | weeks (see Units) | `carryover_t1half` | all |
| $t_{1/2}^a$ | half-life assumed by the analysis when it builds $D_{bc,it}$ (2026-10-09) | weeks | `t_half_a`; suffixes `Dbc0.5`, `Dbc1`, `Dbc2` | docs/43, docs/45 to docs/47, docs/56 to docs/58, CM |
| $\lambda$ | carryover decay rate, $\ln 2 / t_{1/2}$ | per week | `lambda_cor` | all |
| $\lambda_w$ | Weibull rate, $(\ln 2)^{1/k} / t_{1/2}$ | per week | | CM |
| $\phi(\cdot)$, $\phi_t$ | decay-family weight (exponential, linear, Weibull, power) and its value at visit $t$, the residual fraction $e^{-\lambda t_{sd}}$ off drug after exposure and 0 before (2026-10-09) | none | `carryover_decay()`; `exp(-lambda_cor * tsd)` | docs/34, docs/36, docs/37, CM |
| $G(t_{od}; m, d, r)$ | modified Gompertz mean of a component (Part 3), $G(0) = 0$ (2026-10-09) | points | `modgompertz` | docs/40, docs/46, docs/49, docs/51, docs/53, docs/55, CM |
| $\mu_{BR,t}$, $\mu_{\text{last}}$ | mean drug response at visit $t$ after carryover; its value at the last on-drug visit, the anchor of anchored carryover (2026-10-09) | points | `brmeans` | docs/36, docs/37, docs/53 to docs/55 |
| $s_f$ | scale factor on the carryover exponent in the 2024 upstream code (commit 8609f12), default 2 (2026-10-09) | none | `scalefactor` | docs/53 to docs/55 |

Decay families implemented on the data-generating side: exponential
$e^{-\lambda t}$; linear $\max(0, 1 - t/(2t_{1/2}))$; Weibull
$e^{-(\lambda_w t)^k}$; power $\max(0, (1 - t/(3t_{1/2}))^p)$.

## Biomarker and moderation

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $B_i$ | participant-level pre-treatment biomarker value | biomarker units | `bm` | all |
| $b_i$ | standardized biomarker, $(B_i - \bar{B})/s_B$ | none | `bm_z` | all |
| $B^c_i$ | biomarker centered at the trial's baseline sample mean, $B_i - \bar B$, not scaled (2026-10-09) | biomarker units | `bmc` | docs/47, docs/56, docs/57, CM |
| $\mu_{bm}$, $\sigma_{bm}$ | population mean and standard deviation of the biomarker ($\mu_{bm}$ 2026-10-09) | biomarker units | `blparam` | docs/36 to docs/38, docs/48, CM |
| $\sigma_{BR}$ | standard deviation of the drug-response component | points | `respparam` | all |
| $c_{bm}$ | moderation parameter of the **covariance-moderation** architecture: on-drug correlation between biomarker and $BR$. A correlation, bounded in $[-1,1]$ | none | `c.bm` | all |
| $\beta_{bm}$ | moderation parameter of the **mean-moderation** architecture: dimensionless multiplier scaling the treatment effect by $b_i$ | none | `c.bm` with `dgp_architecture = 'mean_moderation'` | all |
| $c_{bm,A}$, $c_{bm,B}$ | independent mean-channel and covariance-channel weights of the dual-channel architecture | none | `c.bm_a`, `c.bm_b` | paper 11 |
| $g_t$, $\mathbf g$ | coupling profile: the fraction of $c_{bm}$ present at visit $t$, $\mathrm{Cor}(B_i, BR_{it}) = c_{bm}\,g_t$; 1 on drug; off drug after exposure $\phi_t$ (graded) or $\mathbf 1[\mu_{BR,t} \ne 0]$ (step); 0 before first exposure (graded); $\mathbf g$ the vector over visits (2026-10-09) | none | built in `buildSigma()` | papers 01, 11; docs/34 to docs/37, docs/40, docs/41, docs/44, docs/49, docs/51, docs/55, CM |
| $\tilde g$ | the coupling profile padded with zeros at the $TV$ and $PB$ entries (2026-10-09) | none | | docs/36, CM |
| $\bar g_{\text{on}}$, $\bar g_{\text{off}}$, $\mathcal G$ | mean coupling over a path's on-drug and off-drug visits, and the coupling gap $\mathcal G = \bar g_{\text{on}} - \bar g_{\text{off}}$ (2026-10-09) | none | | docs/40, docs/49, docs/55, CM |
| $c_{bm}^{\ast}$ | coupling ceiling: the largest $c_{bm}$ for which every path's joint correlation matrix is positive definite, $(\tilde g^\top M^{-1}\tilde g)^{-1/2}$ (2026-10-09) | none | | docs/34 to docs/37, docs/39, docs/44, docs/54, CM |
| $c_{bm,TV}$, $c_{bm,PB}$, $c_{bm,BL}$ | prognostic couplings: correlation of the biomarker with every post-baseline $TV$ element, every post-baseline $PB$ element, and the baseline level (2026-10-09; $c_{bl}$ retired) | none | `c_tv`, `c_pb`, `c_bl` (stress-test driver); `c.bm.pb` (package) | docs/35, docs/44, docs/46 to docs/49, docs/51, docs/56, docs/57, CM |
| $s_{\text{prog}}$, $f$ | prognostic slope of $Y$ on the biomarker that does not depend on drug, and the share of it absorbed by a shared biomarker main effect $\beta_B$ (2026-10-09; not $q$) | points per biomarker unit; none | | docs/48, docs/49, docs/51, docs/56, CM |

## Estimand, inference, and performance

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $\beta_{bm:D}$ | the target estimand: the biomarker-by-treatment interaction coefficient | points per biomarker unit | `bm:Dbc` | all |
| $\beta_{bm}^{BR}$, $\beta_{bm}^{PB}$, $\beta_{bm}^{TV}$ | component-specific biomarker slopes (paper 06); $\beta_{bm}^{TV}$ is also the stress-test DGP parameter for a mean moderation of $TV$ (docs/44) | points per biomarker unit | | paper 06, docs/44, CM |
| $\beta_{bm}^{\text{lumped}}$ | the one-component (single-indicator) biomarker slope | points per biomarker unit | | paper 06 |
| $\beta_D$ | treatment main effect | points | coefficient of `Dbc` or `D` | papers 04, 05; docs/56 to docs/58, CM |
| $\tau^2$ | between-participant variance of the drug effect: treatment heterogeneity (Senn); in the paper 04 protocol the variance of the participant deviation $v_i$ of the drug-response ceiling (2026-10-09) | points$^2$ | `tau_squared`, `tau` | paper 04, its protocol, docs/58, CM (main-effect track) |
| $\theta$, $\hat\theta$ | a generic interaction coefficient; used only where the argument is deliberately not N-of-1 specific | points per biomarker unit | | paper 10, CM |
| $\theta_{\text{true}}$ | calibrated true value, $-c_{bm}\,\sigma_{BR}$ | points per biomarker unit | `true_value` | paper 10 |
| $\kappa$ | ratio of mean model-based standard error to empirical standard deviation; $\kappa > 1$ means a conservative test | none | | papers 08, 10 |
| $\delta$ | standardized bias of the test statistic | none | | paper 10 |
| $\nu$ | denominator degrees of freedom | none | | all |
| $\alpha$ | nominal significance level; 0.05 throughout | none | | all |
| $\pi$, $\hat\pi$ | power, $\Pr(p < \alpha)$, and an estimated rejection rate, with $\mathrm{MCSE}(\hat\pi) = \sqrt{\hat\pi(1-\hat\pi)/n_{\text{reps}}}$ ($\hat\pi$ 2026-10-09) | none | `power`, `mc_se_power` | all |
| $\gamma$ | a $BR \times PB$ interaction coefficient, the subadditivity/superadditivity term | points | | paper 06 |
| $\psi_g$ | latent-class gating slope: steepness of $\pi(B_i)$'s logistic dependence on the standardized biomarker (called `gs` in code and earlier prose) | none | `gs` | paper 03 |
| $Z_i$ | unobserved latent class label; also the responder label of a responder mixture (docs/51) | none | | paper 03, docs/51 |

## Analysis-model, covariance-construct and stress-test symbols

Added 2026-10-09. Coefficients of the analysis model (written in R as
`Sx ~ bmc * (Dbc + De) + t`, with variants), of the paired-difference
statistic E9, and generic conventions.

### Covariance construct

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $\rho$ | AR(1) within-component correlation per week, $\mathrm{Cor} = \rho^{\lvert w_t - w_s\rvert}$ (refined process 0.7; configuration C in `docs/36` 0.8) | per week | `rho` | all |
| $\rho_{cs}$ | common within-component correlation of a compound-symmetry construct, the same at every lag (published process 0.8) (2026-10-09) | none | `rho` in the published code | docs/34 to docs/36, docs/44, docs/48, docs/51, docs/54, CM |
| $\rho_t$ | AR(1) step correlation between visits $t-1$ and $t$, $\rho^{d_t}$ | none | | docs/35, docs/36, CM |
| $c_1$, $c_\times$ | cross-component correlation at the same visit and at different visits (the latter times $\rho^{\lvert w_t - w_s\rvert}$ under the separable form) | none | `c.cf1t`, `c.cfct` | docs/34 to docs/37, docs/40, docs/44, docs/50, docs/51, docs/54, CM |
| $K$, $A$, $M$ | $n_c \times n_c$ cross-component correlation, temporal correlation ($A_{ts} = \rho^{\lvert w_t - w_s\rvert}$ under AR(1)), and the response block, $M = K \otimes A$ when separable | none | | docs/34 to docs/37, docs/39, docs/44, docs/51, CM |

### Analysis model, contrasts and conventions

| Symbol | Meaning | Units | Code | Used in |
|---|---|---|---|---|
| $\beta_0$ | intercept of the analysis model | points | `(Intercept)` | docs/56, docs/57 |
| $\beta_B$ | biomarker main effect (prognostic coefficient of the raw or centered biomarker); distinct from the mean-moderation multiplier $\beta_{bm}$ | points per biomarker unit | coefficient of `bm` or `bmc` | docs/40, docs/41, docs/43, docs/45, docs/46, docs/48, docs/49, docs/51, docs/56, docs/57, CM |
| $\beta_T$ | coefficient of linear time, $\beta_T\,w_t$ (not $\beta_t$) | points per week | coefficient of `t` | docs/40, docs/42, docs/45, docs/56, docs/57, CM |
| $\mu_t$ | visit effect, the free mean of visit $t$ when time enters as a factor | points | `factor(t)` | docs/40, docs/42 |
| $\beta_\eta$, $\beta_{bm:\eta}$ | coefficients of the design expectancy and of the biomarker-by-expectancy term | points; points per biomarker unit | `De`, `bmc:De` | docs/40, docs/42, docs/49, docs/56, docs/57, CM |
| $\beta_{B,0}$ | extra biomarker slope at baseline | points per biomarker unit | `bmc:bl0` | docs/49, CM |
| $\beta_{BL}$, $\beta_{BL,t}$ | coefficient on baseline in the covariate coding, and the slope of visit $t$ on baseline | none | | docs/41, CM |
| $\rho_a$ | working CAR(1) correlation per week of the analysis model, distinct from the data-generating $\rho$ | per week | `corCAR1` parameter `Phi` | docs/56 to docs/58, CM |
| $\Delta_i$, $a_t$ | participant $i$'s on-minus-off contrast $\sum_t a_t Y_{it}$ (E9) and its visit weights $1/n_{\text{on}}$, $-1/n_{\text{off}}$ | points; none | | papers 01, 02; docs/37 to docs/41, docs/43, docs/55, CM |
| $\sigma_\Delta^2$, $\sigma_{\Delta,p}^2$ | residual variance of $\Delta_i$ given the biomarker, overall and in path $p$ (2026-10-09; formerly $\tau^2$, $\tau_p^2$) | points$^2$ | | paper 01; docs/37 to docs/39, CM |
| $r_{\Delta B}$ | correlation between the contrast $\Delta_i$ and the biomarker (not $\rho_{\Delta B}$) | none | | paper 01; docs/37 to docs/39, CM |
| $\Lambda$ | noncentrality parameter of a $t$ or $F$ test | none | | paper 01; docs/37, docs/58, CM |
| $L_{it}$ | lagged just-off-drug indicator of specification E3, $D_{i,t-1}(1 - D_{it})$ | none | | docs/43, CM |
| $\ell_g$ | eigenvalue weights of the CR2 Satterthwaite degrees of freedom | none | | docs/56, CM |
| $\Psi_{it}$, $\boldsymbol\Psi_p$, $\bar\Psi_d$ | unit drug-response profile: the mean of $BR_{it}$ on participant $i$'s path divided by its Gompertz ceiling $m_{BR}$ (onset and anchored carryover included, 0 at baseline and before exposure); its vector over the visits of path $p$; the design's on-minus-off mean of it. Participant-level heterogeneity enters as $v_i\Psi_{it}$ with $\mathrm{Var}(v_i) = \tau^2$ (2026-10-09; formerly $q_{it}$, $\mathbf q_p$, $\bar q_d$, which collided with the path index $q$; capital $\Psi$, distinct from the gating slope $\psi_g$) | none | | paper 04 protocol, docs/58 |
| $\mathcal N(\mu, \sigma^2)$ | the normal law; never $N(\cdot)$, since $N$ is the total sample size | | `rnorm` | all |
| $\varphi(\cdot)$, $\Phi(\cdot)$ | standard normal density and distribution function ($\phi$ is the decay family) | none | `dnorm`, `pnorm` | docs/51, docs/58, CM |
| $\mathbb E$ | expectation operator (distinct from specification labels E1 to E9) | | | docs/51, CM |

**Reserved letters (2026-10-09).** Besides the letters of the tables
above, these may not be reused for another meaning in a document that
also uses the canonical one: $w$ (calendar week; use $n$ for row
counts and another letter for weights), $f$ (absorbed share; a
drug-response time course is $h(t)$, as in `docs/44`), $\tau$
(heterogeneity variance), $\eta$ (design expectancy), $\rho$ (AR(1)
correlation; $\rho_{cs}$, $\rho_t$, $\rho_a$ are its decorated
relatives), $\phi$ (decay family) against $\varphi$ (normal density),
$q$ (path index, or $p$ as in the $q$, $p$ row), $g$ (coupling profile), $\Delta$ (E9 contrast), and
$m$, $d$, $r$ when they appear undecorated beside a Gompertz curve.
(Added 2026-10-09: $\Psi$ is the unit drug-response profile of the
main-effect track; the paper 04 working protocol and `docs/58`, which
wrote it $q$, were renamed accordingly.)

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
4. **$N$ is the total sample size.** $N$ is the total number of
   participants across all randomization paths, and the expected
   allocation to each of $P$ paths is $N/P$: with $P = 4$ and
   $N = 140$, each path receives about 35. Never use $N$ for a per-path
   count. A paper may name a per-path count alongside the total, in
   words or as $N/P$, but must give $N$ itself as the total. Because
   two of the drivers fix the per-path count instead, the reported $N$
   cannot be taken on trust: check the driver against the table in
   Part 2 before quoting a sample size.
5. **State the sign convention.** The three components are
   non-negative *reductions* in symptom severity, so an increase in
   any component lowers $Y$. Treatment effects and interaction
   coefficients are consequently negative; moderation parameters are
   positive.

## What the linter checks

`tools/notation-lint.pl` mechanizes the part of this file that can be
checked from the manuscript sources. Run it over the manuscripts:

```
perl tools/notation-lint.pl analysis/report/*/report.Rmd
```

It skips fenced R chunks and `verbatim` blocks and strips `\texttt{}`
and backtick spans (and, in `.tex` files, the compendium's `\code{}`
macro), including a `\texttt{}` or `\code{}` span that continues across
a line break, so it sees reader-facing prose and mathematics only
(2026-10-09: `\code{}` and multi-line spans added). Four checks are
implemented.

| Check | Fires on | Rule |
|---|---|---|
| `bare-Dbc` | `Dbc` as a symbol outside a code span | 1 |
| `bare-Sx` | `Sx` as a symbol outside a code span | 3 |
| `legacy-spec-label` | `A1`, `A2`, `A3` | Part 3, specifications |
| `N-per-path` | a value of $N$ directly qualified as per-path | 4 |

`N-per-path` fires on `$N = 35$ per randomization path` but not on
`$N = 140$ (70 per path)`, which is the compliant way to give both
numbers. A day-denominated $t_{1/2}$ is exempted when the file carries
a `**Units.**` paragraph (or, in LaTeX, a `\textbf{Units.}`
paragraph), which is how rule 3's units requirement is enforced.

Rule 2 and rule 5 are not mechanized: whether a shared $c_{bm}$ label
is declared, and whether the sign convention is stated, both need a
reader. One known false positive stands, a `bare-Sx` inside a backtick
span that crosses a line break in paper 02.

## Units

Carryover half-life $t_{1/2}$ is quoted in **weeks** in the
interaction-focused manuscripts, on the canonical grid
$\{0, 0.5, 1.0\}$, and in the main-effect manuscript (paper 04); it
is quoted in **days** in the test-procedure manuscripts, on a pharmacokinetic
scale with a 3-day baseline. Both conventions are retained because
each paper's figures are drawn on its own scale. Every paper using
days carries an explicit units paragraph giving the conversion
(1 week = 7 days). Timepoints $t$ are weekly in both conventions.

*Revised 2026-10-09 (user decision):* paper 04 is being rebuilt under
the refined process and uses weeks, as the interaction manuscripts do
(`analysis/scripts/treatment-main-effect/00-ademp-protocol.md`,
"Units"). The day scale with a 3-day baseline belonged to the version
1 study of paper 04 (`00-ademp-protocol-v1.md`), which this paragraph
previously covered.

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
| `modelparam` | list: sample size `N`, moderation `c.bm`. The meaning of `N` depends on the caller; see sample size and randomization paths below |
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

## Sample size and randomization paths

Rule 4 fixes $N$ as the total. The complication is that `modelparam$N`
does not mean the same thing at every level of the call stack, and two
families of driver read it differently.

**The two levels.** `generateData` generates one path: it takes a
single `trialdesign` and treats `N` as the count for that path alone.
`generateSimulatedResults` generates a whole design: it takes `N` as
the total and allocates it, `Ns <- N %/% nP` with the remainder spread
one per path. A driver that loops over `td$trialpaths` itself is
therefore working at the lower level and must divide the total before
the loop. A driver that passes the same full `N` to every iteration
generates $N \times P$ participants, and its reported $N$ is a
per-path count wearing the wrong name.

**Path counts.** $P$ is a property of the design, not a free
parameter.

| Design | $P$ | Paths |
|---|---|---|
| OL | 1 | on drug throughout |
| CO | 2 | drug-then-placebo, placebo-then-drug |
| OL+BDC | 2 | discontinuation at either of two visits |
| Hybrid | 4 | two discontinuation points crossed with two crossover orders |

Paper 06 also uses a reduced single-path OL+BDC variant for its
smaller cells and a two-path sixteen-visit variant for its larger
ones, so $P$ changes within that paper's sample-size axis.

**Which drivers allocate and which duplicate.** Established by
inspection of every driver reachable from a current manuscript, on
2026-08-05.

| Driver | Reads `N` as | Mechanism |
|---|---|---|
| `R/generateSimulatedResults.R` | total | `Ns <- N %/% nP` plus remainder |
| `analysis/scripts/carryover-sensitivity/simulation-core.R` | total | `generate_data_multi_path`, `floor(N / n_paths)` |
| `implementations/nof1power/R/simulation.R` | total | same integer split |
| `analysis/scripts/test-procedure-design-sensitivity/01-study1-test-procedure.R` | total | explicit `split_N` |
| `analysis/scripts/quick-sim/*.R` | per path | hand-rolled `for (g in seq_along(paths))` passing one `mp` |
| `analysis/scripts/component-decomposition/*.R` | per path | same shape |
| `analysis/scripts/gompertz-evaluation/02-faithful-trajectory-sweep.R` | per path | same shape |

**Status by paper.** The manuscripts report totals throughout, and as
of 2026-08-08 the drivers and stored output for papers 01, 06, 07 and
09 were brought into line with them. One entry is weaker than the
others: the driver behind paper 01's Section 3.3 robustness run was
not located, so its per-path reading is taken from that paper's own
statement of its convention rather than from inspected code.

| Paper | Driver was | Now |
|---|---|---|
| 01 | per path | driver allocates the total; re-run at a matched $N = 70$ for all three designs |
| 02, 05, 08, 10, 11 | total | already correct; untouched |
| 03, 04 | coincident (single-path or path-free DGP) | unaffected |
| 06 | per path | all drivers allocate the total; stored `N` migrated to totals; see the caveat below |
| 07 | per path | driver allocates the total; re-run; stored `N` added as 70 |
| 09 | per path | driver allocates the total; re-run at a matched $N = 70$ for all four designs |

**How the corrections were validated.** The work ran in two stages,
and it is worth keeping them apart. The first stage corrected the
label only. Re-expressing a grid in totals and allocating across paths
leaves the per-path counts unchanged, so the generated data should be
identical and nothing but the recorded `N` should move. That was
checked before each re-run by generating one replicate both ways and
comparing: for every driver converted, the drawn data and the fitted
coefficients were identical. The re-runs then reproduced the published
cell power values exactly, in all 36 cells of paper 01, all 16 of
paper 07 and all 16 of paper 09, with per-replicate coefficients
moving by at most $7 \times 10^{-9}$ (paper 01), $0$ (paper 07) and
$4 \times 10^{-5}$ (paper 09), which is `lme` convergence-tolerance
noise on identical data.

The second stage, for papers 01 and 09 only, changed the sample sizes
themselves to match the designs, so it was expected to change results
and did. Its internal check is that the cells whose totals were
already 70 must still reproduce exactly while the others move: in
paper 01 the CO and OL+BDC cells came back bit-identical and only the
Hybrid column changed, by up to 0.233 in power.

**Paper 06 was migrated, not re-run.** Study A alone holds 200,000
stored fits, so re-running it to change a column that the equivalence
check proves is mislabeled would cost hours to reproduce the same
numbers. Instead `analysis/scripts/component-decomposition/
migrate-n-to-total.R` recoded the stored `N` in place, 109 files
across five output directories, with pre-migration copies in
`analysis/.n-migration-backup/` and read-back verification of row
counts and of the absence of legacy values. A column-by-column
comparison against the backups confirms that nothing but `N` changed.
The mapping is per directory because the studies differ in path count:
the study-A and contamination grids switch to the two-path design
above a threshold, so 100 and 150 became 200 and 300 while 35 and 70
were already single-path totals; study-B recovery is two-path
throughout, so 70 and 150 became 140 and 300; study-B balanced is
single-path and was left alone.

**All four paper-06 drivers are now converted**: `run-study-a-prod.R`,
`run-study-contam.R` and `run-study-b-recovery.R` read `N` as a total
and allocate it, and `run-study-b-balanced.R` needed nothing because
its balanced-placebo design has a single path. Each conversion was
checked two ways. An equivalence gate regenerated one replicate per
design type under both parameterizations and got identical
coefficients, standard errors, p-values and formula-dropped status.
Then every stored checkpoint was matched against the new grid's `N`
for its `cell_id`: 18 of 18 for the contamination study and 8 of 8 for
study-B recovery, no mismatches, so each driver agrees with the data
already migrated on disk. Where a fit function takes an `N` argument
it is passed the per-path count, which is what it carried before, so
the formula ladders keyed on it are unaffected. Study-B recovery's
three fits declare the argument but never read it; it is kept per-path
there for consistency with the other drivers.

**The substantive problem was fixed on 2026-08-08.** Fixing the
per-path count across designs with different $P$ leaves a comparison
that is not matched on $N$, so papers 01 and 09 were re-run with every
design at a common total of $N = 70$. This changes numbers rather than
labels, and it changed conclusions in both papers.

In paper 01 only the Hybrid column moved, since CO and OL+BDC were
already at 70 and reproduced exactly. Hybrid power fell from about
0.99 to 0.84, which took it off the ceiling and roughly tripled its
measured carryover loss: 1.4 to 6.9 percent under mean moderation and
5.1 to 15.4 percent under covariance moderation. The paper's headline
architecture contrast survives, with the difference-of-differences
$z$ at Hybrid moving from 3.96 to 2.80 and OL+BDC unchanged at 7.64,
but its claim that mean-moderation losses are uniformly 1 to 3 percent
did not, and the Hybrid design is no longer saturated anywhere in the
grid.

In paper 09 the OL and hybrid rows moved. The happy-accident ordering
held, hybrid 0.968 above OL+BDC 0.950 above crossover 0.854 above OL
0.088 at every dropout pattern, but the hybrid margin over OL+BDC
collapsed to within one MCSE, and the paper's second conclusion, that
the most powerful design is also the most dropout-robust, was
withdrawn: hybrid loss under 40 percent biased dropout rose from 0.4
to 5.8 percentage points, level with OL+BDC's 5.6, so the earlier
near-immunity was a ceiling effect. The OL row also stopped
coinciding with $\alpha$, moving from 0.022-0.046 to 0.086-0.112,
which showed that reading it as a null-calibration check had been a
numerical coincidence at one sample size rather than a property of the
design.

The lesson worth carrying: an unmatched design comparison does not
merely inflate one design, it can manufacture a robustness finding out
of a ceiling effect.

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
blinded discontinuation, then a brief crossover. The program's
reference design. Capitalize it.
**Aggregated N-of-1 trial**. Many single-patient on/off series pooled
into one mixed model.

Define each abbreviation at first use in every paper. Each design's
randomization-path count is given in Part 2, under sample size and
randomization paths, and is not repeated here.

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

## Named data-generating processes (added 2026-10-08)

The current best DGP is called **refined** (user decision,
2026-10-08): "the refined process" in prose, "Refined" in legends and
tables, against "the published process" / "Published" for
Hendrickson's code.

| Name in prose | Meaning | Defined in |
|---|---|---|
| the published process | Hendrickson's `generateData()` and `buildtrialdesign()` at commit 58b32a9: compound symmetry within components ($\rho_{cs} = 0.8$), constant cross-component correlation, step coupling, compounding carryover recursion, time on drug restarting at zero on re-exposure | Paper 01 |
| the refined process | The recommended DGP: graded coupling (0 before first exposure), separable AR(1) with $\rho = 0.7$ and $c_1 = 0.2$, anchored carryover, a stop on matrices that are not positive definite, randomized path allocation; re-exposure as in the published process; optional prognostic couplings, phase-specific measurement error, Weibull decay and dropout | Paper 01, Table `tab:spec` |
| the published construct | The covariance construct of the published process alone (configuration A: compound symmetry within components, $\rho_{cs} = 0.8$, constant cross-component correlation), with whatever coupling the context states; not the whole published process (added 2026-10-09) | CM, glossary; used in CM, `docs/43`, `docs/56` |
| configuration A to F | The covariance and coupling constructs compared in `docs/36` (A published; B AR(1) with the `covar` cross block, step coupling; C separable AR(1) with graded coupling; D the published compound-symmetry construct with mean moderation; E separable AR(1) with mean moderation and no covariance coupling; F AR(1) with the `covar` cross block, graded coupling). (Corrected 2026-10-09: E was listed as "all three" adjustments; graded coupling has no effect in E, whose coupling vector is zero, so E is adjustment 1 with the separable cross block $c_1\rho^{\lvert w_t - w_s\rvert}$ plus adjustment 3, as arm E of `hendrickson-problems/04-power-simulation.R` builds it) | `docs/36`, Section 6 and Appendix A.6 |
| `orig`, `covar` | Construct labels of the whitepapers. `orig`: the published process's construct as implemented in a given copy of the code (in `docs/54` the 2024 vendored copy, not commit 58b32a9, so the document states which); `covar`: the package's covariance-moderation construct, AR(1) within component with the decaying cross-component block (configuration F at the package's parameters) (added 2026-10-09) | `docs/32`, `docs/34` to `docs/38`, `docs/54` |

Use **the refined process** for the DGP of the stress test and the
recommendation. Use **configuration C** only for the covariance and
coupling construct itself, as in the `docs/36` comparison, and state
$\rho$ with it (`docs/36` evaluates it at $\rho = 0.8$; the refined
process uses $\rho = 0.7$). Configuration C is not the refined process:
it does not by itself include anchored carryover, the positive-definite
stop or the options. Do not confuse configuration C with the former
label "Architecture C" (the dual-channel architecture).

## Analysis specifications for carryover

**E1.** Binary on-drug indicator; carryover ignored.
**E2.** Exposure-weighted continuous predictor $D_{bc}$, committing to
an assumed decay form and half-life.
**E3.** Binary indicator plus a lagged just-off-drug nuisance term;
the classical Jones-Kenward crossover device, assuming no decay form.

The letter E stands for the exposure regressor that distinguishes
them. These were labeled `A1`/`A2`/`A3` until 2026-07-31; the stored data
and the drivers were migrated to E1/E2/E3, so the labels now agree
everywhere and nothing is mapped at display time.

**Why not A/B/C and M1/M2/M3.** The original scheme used A/B/C for
architectures and `A1`/`A2`/`A3` for specifications, colliding on the letter
A within a single manuscript. Relabeling the specifications M1/M2/M3
removed that collision but created a second one: the calibration
manuscript already uses M0-M3 for its working-covariance model ladder,
and those labels are keys into archived data. E1/E2/E3 was free across
the whole corpus when it was chosen, and naming the architectures by
mechanism removes the class of collision entirely.

**The extended E series (2026-10-09).** The E labels are no longer
free beyond E3: later documents extend the series, and every extension
is listed here so that no label is reused with a second meaning.

| Label | Meaning | Defined in |
|---|---|---|
| E4 | binary drug state plus linear time since discontinuation (retired) | `docs/43` (from `implementations/original`) |
| E5 | lag term crossed with the biomarker, $L_{it}$ and $B_i L_{it}$ (retired) | `docs/43`, paper 02 |
| E6 | $t_{sd}$ crossed with the biomarker (retired) | `docs/43`, paper 02 |
| E7 | E2 form with the half-life chosen by AIC | paper 02 (stored code), `docs/43`, CM |
| E9 | paired-difference statistic: each participant's on-drug mean minus off-drug mean, $\Delta_i$, regressed on the biomarker by OLS | paper 01, paper 02 (stored code), `docs/37`, CM |
| E9b | E9 with baseline counted as one more off-drug visit | paper 01, `docs/37`, CM |
| E1cr2, E2cr2, E3cr2, E7cr2 | the same formulas with the CR2 cluster-robust standard error | paper 02 (stored codes), CM |

E8 is unassigned. E9 is a summary statistic rather than an exposure
coding; it keeps its label because the label is a key in paper 02's
stored data.

## Response-curve families

**Modified Gompertz** (default), $y = m\,\exp(-d\,\exp(-r\,t))$ with a
vertical offset so the curve passes through the origin; written
$G(t_{od}; m, d, r)$ for a component's mean, with ceiling $m$
(points), displacement $d$ (none) and rate $r$ (per week), and
inflection time $t^\ast = \ln(d)/r$ (weeks). The code names are
`maxr`, `disp` and `rate` (`respparam`); write the mathematical
symbols in equations and the code names only in code (2026-10-09; the
equation was previously written with the code names). A component
subscript names the curve where needed ($m_{BR}$, $G_{TV}$). Alternatives:
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

**LME.** Linear mixed-effects model; the program's default analysis.
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

**Biomarker-treatment interaction.** The program's core estimand;
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

---
*Rendered on 2026-08-08 at 21:18 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/NOTATION.md*

