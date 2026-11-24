# Notation (lite)
*2026-10-28 10:22 PDT*

Symbol and meaning columns only. Units, code identifiers and
usage columns are in `NOTATION.md`, which governs.

## Indices and structure

| Symbol | Meaning |
|------------------|----------------------------------------------------------------------------------|
| $i$ | participant index |
| $t$ | visit (timepoint) index; weekly measurement occasions throughout |
| $j$, $ij$ | period or occasion index within participant |
| $w_t$, $w_s$ | calendar week of visit $t$ (Hybrid: 4, 8, 9, 10, 11, 12, 16, 20); reserved for calendar time, never a weight (2026-10-09) |
| $d_t$ | gap between visits, $w_t - w_{t-1}$ (2026-10-09) |
| $N$ | **total** number of participants in the trial, across all randomization paths; never a per-path count. See Part 2, sample size and randomization paths |
| $P$ | number of randomization paths in a design; the expected allocation to each is $N/P$ |
| $q$, $p$ | randomization-path index, $q = 1, \ldots, P$, where a document also uses $p$ for a $p$-value or an exponent (2026-10-09). The subscript $p$ on path-specific quantities fixed by the 2026-10-09 decisions and the compendium ($\sigma_{\Delta,p}^2$, $\tilde\beta_p$, the path share $\varpi_p$) is retained; a document uses one of the two letters for paths throughout |
| $\varpi_p$ | share of participants allocated to path $p$, $n_p/N$ ($1/P$ under equal allocation) (2026-10-09; formerly $q_p$ in `docs/37` and the compendium) |
| $n_{\text{on}}$, $n_{\text{off}}$, $n_v$ | numbers of on-drug and off-drug post-baseline visits (rows) of a path, and of rows of visit kind $v$ (2026-10-09) |
| $n_c$ | number of response components (3; 2 in papers 13 and 14) (2026-10-09) |
| $n_{\text{reps}}$, $n_{\text{sim}}$ | Monte Carlo replicates per cell |
| $k$ | cycles per participant in design sweeps; also the Weibull shape parameter in the decay family |

## Outcome and latent components

| Symbol | Meaning |
|------------------|----------------------------------------------------------------------------------|
| $Y_{it}$ | observed symptom score; the modeled outcome |
| $\mathrm{BL}_i$ | participant-specific baseline level |
| $BR_{it}$ | biological (pharmacological) response component |
| $PB_{it}$ | placebo-belief response component |
| $TV_{it}$ | time-variant natural-history component |
| $\varepsilon_{it}$ | observation-level residual |
| $u_i$ | participant random intercept |
| $\sigma_u^2$, $\sigma_\varepsilon^2$ | variances of the random intercept and of the residual of the analysis model (2026-10-09; not $\sigma_b^2$, $\sigma_e^2$) |
| $\eta_t$ | design expectancy at visit $t$: 1 open label, 0.5 blinded, 0 at baseline; paper 06 writes $\eta(\text{phase})$ for the same object (2026-10-09; formerly $e_t$) |
| $\tilde\eta_{it}$ | measured (participant-reported) placebo-belief covariate, if an analysis uses one; reserved (2026-10-09) |
| $\sigma_{TV}$, $\sigma_{PB}$, $\sigma_{BL}$ | standard deviations of the natural-course and placebo-belief components and of the baseline level; $\sigma_{PB,t} = \sigma_{PB}\,\eta_t$ (2026-10-09) |

## Treatment exposure and carryover

| Symbol | Meaning |
|------------------|----------------------------------------------------------------------------------|
| $D_{it}$ | binary drug state: 1 on drug, 0 off drug |
| $D_{bc,it}$ | continuous exposure-decayed drug indicator: 1 on drug, $e^{-\lambda t_{sd}}$ off drug after exposure, 0 at off-drug visits before the first exposure; in an analysis it is built at the assumed half-life $t_{1/2}^a$ (2026-10-09: pre-exposure value stated) |
| $t_{sd}$ | time since discontinuation; input to the carryover decay |
| $t_{od}$ | cumulative time on drug; input to the Gompertz response |
| $t_{pb}$ | time-by-placebo-belief accumulation |
| $t_{1/2}$ | carryover half-life of the data-generating process |
| $t_{1/2}^a$ | half-life assumed by the analysis when it builds $D_{bc,it}$ (2026-10-09) |
| $\lambda$ | carryover decay rate, $\ln 2 / t_{1/2}$ |
| $\lambda_w$ | Weibull rate, $(\ln 2)^{1/k} / t_{1/2}$ |
| $\phi(\cdot)$, $\phi_t$ | decay-family weight (exponential, linear, Weibull, power) and its value at visit $t$, the residual fraction $e^{-\lambda t_{sd}}$ off drug after exposure and 0 before (2026-10-09) |
| $G(t_{od}; m, d, r)$ | modified Gompertz mean of a component (Part 3), $G(0) = 0$ (2026-10-09) |
| $\mu_{BR,t}$, $\mu_{\text{last}}$ | mean drug response at visit $t$ after carryover; its value at the last on-drug visit, the anchor of anchored carryover (2026-10-09) |
| $s_f$ | scale factor on the carryover exponent in the 2024 upstream code (commit 8609f12), default 2 (2026-10-09) |

## Biomarker and moderation

| Symbol | Meaning |
|------------------|----------------------------------------------------------------------------------|
| $B_i$ | participant-level pre-treatment biomarker value |
| $b_i$ | standardized biomarker, $(B_i - \bar{B})/s_B$ |
| $B^c_i$ | biomarker centered at the trial's baseline sample mean, $B_i - \bar B$, not scaled (2026-10-09) |
| $\mu_{bm}$, $\sigma_{bm}$ | population mean and standard deviation of the biomarker ($\mu_{bm}$ 2026-10-09) |
| $\sigma_{BR}$ | standard deviation of the drug-response component |
| $c_{bm}$ | moderation parameter of the **covariance-moderation** architecture: on-drug correlation between biomarker and $BR$. A correlation, bounded in $[-1,1]$ |
| $\beta_{bm}$ | moderation parameter of the **mean-moderation** architecture: dimensionless multiplier scaling the treatment effect by $b_i$ |
| $c_{bm,A}$, $c_{bm,B}$ | independent mean-channel and covariance-channel weights of the dual-channel architecture |
| $g_t$, $\mathbf g$ | coupling profile: the fraction of $c_{bm}$ present at visit $t$, $\mathrm{Cor}(B_i, BR_{it}) = c_{bm}\,g_t$; 1 on drug; off drug after exposure $\phi_t$ (graded) or $\mathbf 1[\mu_{BR,t} \ne 0]$ (step); 0 before first exposure (graded); $\mathbf g$ the vector over visits (2026-10-09) |
| $\tilde g$ | the coupling profile padded with zeros at the $TV$ and $PB$ entries (2026-10-09) |
| $\bar g_{\text{on}}$, $\bar g_{\text{off}}$, $\mathcal G$ | mean coupling over a path's on-drug and off-drug visits, and the coupling gap $\mathcal G = \bar g_{\text{on}} - \bar g_{\text{off}}$ (2026-10-09) |
| $c_{bm}^{\ast}$ | coupling ceiling: the largest $c_{bm}$ for which every path's joint correlation matrix is positive definite, $(\tilde g^\top M^{-1}\tilde g)^{-1/2}$ (2026-10-09) |
| $c_{bm,TV}$, $c_{bm,PB}$, $c_{bm,BL}$ | prognostic couplings: correlation of the biomarker with every post-baseline $TV$ element, every post-baseline $PB$ element, and the baseline level (2026-10-09; $c_{bl}$ retired) |
| $s_{\text{prog}}$, $f$ | prognostic slope of $Y$ on the biomarker that does not depend on drug, and the share of it absorbed by a shared biomarker main effect $\beta_B$ (2026-10-09; not $q$) |

## Estimand, inference, and performance

| Symbol | Meaning |
|------------------|----------------------------------------------------------------------------------|
| $\beta_{bm:D}$ | the target estimand: the biomarker-by-treatment interaction coefficient |
| $\beta_{bm}^{BR}$, $\beta_{bm}^{PB}$, $\beta_{bm}^{TV}$ | component-specific biomarker slopes (paper 06); $\beta_{bm}^{TV}$ is also the stress-test DGP parameter for a mean moderation of $TV$ (docs/44) |
| $\beta_{bm}^{\text{lumped}}$ | the one-component (single-indicator) biomarker slope |
| $\beta_D$ | treatment main effect |
| $\tau^2$ | between-participant variance of the drug effect: treatment heterogeneity (Senn); in the paper 04 protocol the variance of the participant deviation $v_i$ of the drug-response ceiling (2026-10-09) |
| $\theta$, $\hat\theta$ | a generic interaction coefficient; used only where the argument is deliberately not N-of-1 specific |
| $\theta_{\text{true}}$ | calibrated true value, $-c_{bm}\,\sigma_{BR}$ |
| $\kappa$ | ratio of mean model-based standard error to empirical standard deviation; $\kappa > 1$ means a conservative test |
| $\delta$ | standardized bias of the test statistic |
| $\nu$ | denominator degrees of freedom |
| $\alpha$ | nominal significance level; 0.05 throughout |
| $\pi$, $\hat\pi$ | power, $\Pr(p < \alpha)$, and an estimated rejection rate, with $\mathrm{MCSE}(\hat\pi) = \sqrt{\hat\pi(1-\hat\pi)/n_{\text{reps}}}$ ($\hat\pi$ 2026-10-09) |
| $\gamma$ | a $BR \times PB$ interaction coefficient, the subadditivity/superadditivity term |
| $\psi_g$ | latent-class gating slope: steepness of $\pi(B_i)$'s logistic dependence on the standardized biomarker (called `gs` in code and earlier prose) |
| $Z_i$ | unobserved latent class label; also the responder label of a responder mixture (docs/51) |

## Analysis-model, covariance-construct and stress-test symbols

### Covariance construct

| Symbol | Meaning |
|------------------|----------------------------------------------------------------------------------|
| $\rho$ | AR(1) within-component correlation per week, $\mathrm{Cor} = \rho^{\lvert w_t - w_s\rvert}$ (refined process 0.7; configuration C in `docs/36` 0.8) |
| $\rho_{cs}$ | common within-component correlation of a compound-symmetry construct, the same at every lag (published process 0.8) (2026-10-09) |
| $\rho_t$ | AR(1) step correlation between visits $t-1$ and $t$, $\rho^{d_t}$ |
| $c_1$, $c_\times$ | cross-component correlation at the same visit and at different visits (the latter times $\rho^{\lvert w_t - w_s\rvert}$ under the separable form) |
| $K$, $A$, $M$ | $n_c \times n_c$ cross-component correlation, temporal correlation ($A_{ts} = \rho^{\lvert w_t - w_s\rvert}$ under AR(1)), and the response block, $M = K \otimes A$ when separable |

### Analysis model, contrasts and conventions

| Symbol | Meaning |
|------------------|----------------------------------------------------------------------------------|
| $\beta_0$ | intercept of the analysis model |
| $\beta_B$ | biomarker main effect (prognostic coefficient of the raw or centered biomarker); distinct from the mean-moderation multiplier $\beta_{bm}$ |
| $\beta_T$ | coefficient of linear time, $\beta_T\,w_t$ (not $\beta_t$) |
| $\mu_t$ | visit effect, the free mean of visit $t$ when time enters as a factor |
| $\beta_\eta$, $\beta_{bm:\eta}$ | coefficients of the design expectancy and of the biomarker-by-expectancy term |
| $\beta_{B,0}$ | extra biomarker slope at baseline |
| $\beta_{BL}$, $\beta_{BL,t}$ | coefficient on baseline in the covariate coding, and the slope of visit $t$ on baseline |
| $\rho_a$ | working CAR(1) correlation per week of the analysis model, distinct from the data-generating $\rho$ |
| $\Delta_i$, $a_t$ | participant $i$'s on-minus-off contrast $\sum_t a_t Y_{it}$ (E9) and its visit weights $1/n_{\text{on}}$, $-1/n_{\text{off}}$ |
| $\sigma_\Delta^2$, $\sigma_{\Delta,p}^2$ | residual variance of $\Delta_i$ given the biomarker, overall and in path $p$ (2026-10-09; formerly $\tau^2$, $\tau_p^2$) |
| $r_{\Delta B}$ | correlation between the contrast $\Delta_i$ and the biomarker (not $\rho_{\Delta B}$) |
| $\Lambda$ | noncentrality parameter of a $t$ or $F$ test |
| $L_{it}$ | lagged just-off-drug indicator of specification E3, $D_{i,t-1}(1 - D_{it})$ |
| $\ell_g$ | eigenvalue weights of the CR2 Satterthwaite degrees of freedom |
| $\Psi_{it}$, $\boldsymbol\Psi_p$, $\bar\Psi_d$ | unit drug-response profile: the mean of $BR_{it}$ on participant $i$'s path divided by its Gompertz ceiling $m_{BR}$ (onset and anchored carryover included, 0 at baseline and before exposure); its vector over the visits of path $p$; the design's on-minus-off mean of it. Participant-level heterogeneity enters as $v_i\Psi_{it}$ with $\mathrm{Var}(v_i) = \tau^2$ (2026-10-09; formerly $q_{it}$, $\mathbf q_p$, $\bar q_d$, which collided with the path index $q$; capital $\Psi$, distinct from the gating slope $\psi_g$) |
| $\mathcal N(\mu, \sigma^2)$ | the normal law; never $N(\cdot)$, since $N$ is the total sample size |
| $\varphi(\cdot)$, $\Phi(\cdot)$ | standard normal density and distribution function ($\phi$ is the decay family) |
| $\mathbb E$ | expectation operator (distinct from specification labels E1 to E9) |
