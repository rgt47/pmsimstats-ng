---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# Testing the Drug Main Effect in Aggregated N-of-1 Trials When the Carryover Half-Life Is Unknown: Estimands, Probability Limits, Randomization Tests and Data-Chosen Half-Lives {.unlisted .unnumbered}
*2026-10-09 12:06 PDT*

**Author.** pmsimstats team

**Purpose.** The rebuilt paper 04 will test the drug main effect
$\beta_D$ in the Hybrid design and its comparators under the refined
process, with analyses that weight off-drug visits by an
exposure-decayed indicator $D_{bc,it}$ built at an analysis half-life
$t_{1/2}^a$ (the working protocol,
`analysis/scripts/treatment-main-effect/00-ademp-protocol.md`, here
"the protocol"). `docs/57` asked how $t_{1/2}^a$ may be chosen from the
data when the target is the biomarker-by-drug interaction. This note
asks the same question for the main effect, where the answer is
different, and develops the mathematics the protocol's closed-form
step needs: what the coefficient of $D_{bc,it}(h)$ estimates under
Gompertz onset, expectancy and carryover; its probability limit under a
working covariance when the true half-life is $t_{1/2}$ and the analysis
uses $h$; why the half-life cannot be estimated under the null; how
size can be controlled when it is chosen anyway (maximum statistics,
randomization tests, the bootstrap); what model averaging and joint
estimation do to the estimand; and a closed-form power development for
the Hybrid design. A worked example on one simulated trial and a modest
Monte Carlo check accompany the derivations. A comparison table, a
provisional recommendation for paper 04 and a list of proposed edits to
the protocol close the note. It is written for trial statisticians and
methodologists and assumes the model statements of `docs/56` and the
protocol.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\listoffigures
\clearpage
```

## Notation index and glossary

Notation follows the project's canonical table
(`analysis/report/NOTATION.md`). Symbols marked *local* are defined
here or in the protocol for the main-effect work; Appendix A lists the
source of every symbol.

| Symbol | Meaning | Status |
|------------|--------------------------------------------------|----------|
| $i$, $t$ | participant; visit, with $t = 0$ the baseline row | canonical |
| $w_t$ | calendar week of visit $t$ | canonical |
| $N$, $P$ | total participants across paths; number of randomization paths | canonical |
| $p$, $p(i)$ | a path; the path of participant $i$ (the protocol's path letter, allowed by `NOTATION.md`'s $q$, $p$ row) | local (protocol) |
| $\varpi_p$ | expected share of path $p$, $1/P$ | local (protocol) |
| $Y_{it}$, $\mathbf Y_i$ | symptom score (code `Sx`); participant $i$'s nine scores | canonical |
| $\mathrm{BL}_i$, $TV_{it}$, $PB_{it}$, $BR_{it}$ | baseline level; natural-course, placebo-belief and drug responses | canonical |
| $\eta_t$ | design expectancy: 1 open label, 0.5 blinded, 0 at baseline (code `De`) | canonical |
| $D_{it}$, $D_{bc,it}$ | binary drug state; exposure-decayed indicator with the pre-exposure guard (code `Dbc`, `X` in the scripts) | canonical |
| $t_{od}$, $t_{sd}$, $t_{pb}$ | time on drug, time since discontinuation, time under expectancy (weeks) | canonical |
| $t_{1/2}$, $t_{1/2}^a$ | true carryover half-life; half-life assumed by an analysis (weeks) | canonical |
| $h$ | a candidate value of $t_{1/2}^a$; $D_{bc,it}(h)$ is the indicator built with it | local (`docs/57`) |
| $\mathcal T$ | the analysis grid $\{0, 0.2, 0.5, 1, 2\}$ weeks | local |
| $r_1$ | $2^{-1/h}$, the share of the drug effect that $D_{bc}(h)$ keeps one week after stopping | local |
| $m_{BR}$, $f_m$ | Gompertz ceiling of $BR$ (published 10.986); drug-effect multiplier $m_{BR}/10.986$ | canonical ($m_{BR}$); local, protocol ($f_m$) |
| $\Psi_{pt}$, $\boldsymbol\Psi_p$ | unit drug-response profile of path $p$ at visit $t$, and its vector over visits: the $BR$ mean divided by $m_{BR}$, anchored carryover included, 0 at baseline (renamed 2026-10-09 from $q_{pt}$, $\mathbf q_p$, which collided with the canonical path index $q$) | canonical |
| $v_i$, $\tau^2$ | participant-specific deviation of the ceiling; its variance (participant-by-treatment variance) | local, protocol ($v_i$); canonical ($\tau^2$) |
| $\beta_D$ | drug main effect | canonical |
| $\beta_D^{\text{full}}$, $\beta_D^{\text{exp}}$, $\beta_D^{\text{bl}}$ | the three candidate estimands of the protocol, Section 4.1 | local (protocol) |
| $\boldsymbol\mu_p$, $\boldsymbol\mu^0$ | mean of $\mathbf Y_i$ on path $p$; its non-drug part, common to all paths | local |
| $\Sigma_p$ | covariance of $\mathbf Y_i$ on path $p$ | local |
| $X_p(h)$, $\mathbf x_p(h)$ | design matrix of an analysis on path $p$; its drug column $D_{bc}(h)$ | local |
| $V$ | working covariance of an analysis | local |
| $\boldsymbol\omega_p(h)$ | weights of the drug coefficient on path $p$ (Section 5.1) | local |
| $b_0(h)$ | null bias of the drug coefficient: its limit when $f_m = 0$ | local |
| $\bar\Psi_\omega(h; t_{1/2})$ | exposure-weighted drug profile, $\sum_p \varpi_p \boldsymbol\omega_p(h)'\boldsymbol\Psi_p$ (formerly $\bar q_\omega$) | local |
| $\bar d(h)$ | on-minus-off mean of $D_{bc}(h)$ over post-baseline visits | local |
| $\hat\beta_D^{\text{exp}}(h)$ | model-implied exposure contrast, $\hat\beta_D(h)\,\bar d(h)$ | local |
| $T(h)$, $\mathrm{LR}(h)$ | CR2 $t$-statistic of the drug coefficient at $h$; likelihood-ratio statistic for it | local (`docs/57`) |
| $\ell_1(h)$, $\ell_0$ | maximized ML log-likelihood with $D_{bc}(h)$; without any drug term | local |
| $c_{\max}$ | 95% critical value of $\max_{h \in \mathcal T}\lvert T(h)\rvert$ | local |
| $\hat C$, $\hat K$ | the contrast part and the relabeling-invariant part of a linear statistic (Section 9.3) | local |
| $\Lambda$ | noncentrality of a test statistic | canonical |
| $\sigma_u^2$, $\rho_a$ | random-intercept variance and CAR(1) parameter of the working covariance | canonical |
| $\rho$, $c_1$ | AR(1) serial correlation per week (0.7); cross-component correlation (0.2) | canonical |
| $\alpha$, $\pi$ | nominal level (0.05); power | canonical |
| $M$ | relabelings or resamples per trial | local |
| $n_{\text{sim}}$ | Monte Carlo replicates per cell (formerly $R$) | canonical |
| $U_i(h)$ | participant $i$'s contribution to the CR2 score of the drug coefficient at $h$ (formerly $u_i(h)$, which collides with the canonical random intercept $u_i$) | local |
| $\Pr(\cdot)$, $\mathbb E(\cdot)$ | probability and expectation ($P$ is the number of paths) | canonical |

Table: Notation index for this note

Two collisions are avoided deliberately. `docs/57` wrote $\tau^2$ for
the random-intercept variance and $\theta$ for $\log t_{1/2}^a$ before
its 2026-10-09 notation pass; here
$\tau^2$ is the participant-by-treatment variance, as in the protocol,
the random-intercept variance is $\sigma_u^2$, and $\theta$ (reserved
for a generic interaction coefficient) is not used. The reserved
symbols $\gamma$, $\delta$, $\kappa$, $\nu$, $\theta$ and $\lambda$ do
not appear; $\pi$ is power only.

**Sign convention** (NOTATION.md, rule 5). The components are
non-negative reductions in symptom severity, so $\beta_D$ is negative
when the drug helps, and a positive bias is a bias toward no effect.

**Units.** Carryover half-lives are in weeks throughout, as in the
protocol and, since 2026-10-09, in `NOTATION.md` for paper 04 (1 week
= 7 days). Effects are in CAPS points.

**Glossary.**

- **Designs.** *Hybrid*: open label at weeks 4 and 8, blinded
  discontinuation and a brief blinded crossover, visits at weeks 4, 8,
  9, 10, 11, 12, 16 and 20, four paths [1]. *OL+BDC*: open label to
  week 16, then blinded discontinuation at week 17 or 18, two paths.
  *CO*: two-period crossover, visits every 2.5 weeks, two paths
  (protocol, Table "Designs").
- **The refined process.** The recommended data-generating process of
  paper 01 (Table `tab:spec`): separable AR(1) with $\rho = 0.7$ and
  $c_1 = 0.2$, anchored carryover, a stop on matrices that are not
  positive definite, randomized path allocation. **The published
  process**: Hendrickson's code at commit 58b32a9 (`NOTATION.md`).
- **Anchored and compounding carryover.** Decay of the off-drug drug
  response from its last on-drug value by the elapsed time, against the
  published recursion that decays an already decayed value (`docs/53`).
- **Analyses.** All have a random intercept and CAR(1) residuals on
  weeks (`nlme::lme`, REML) and CR2 inference with Satterthwaite
  degrees of freedom [3, 4], with the baseline as a response row.
  **F1**: $Y \sim t + D_{bc}(h) + \eta_t$, the protocol's proposed
  reference. **S**: $Y \sim t + D_{bc}(h)$. **VF**: visit effects plus
  $D_{bc}(h)$ (protocol, Section 5.2). **F1A**, **VFA**: the same with
  $h$ chosen by AIC over $\mathcal T$ (ML for the choice, REML refit).
- **Sharp null.** No drug effect for anyone: $f_m = 0$ and $\tau = 0$,
  so each participant's outcome distribution does not depend on the
  path. **Weak null.** Zero mean effect with heterogeneity: $f_m = 0$,
  $\tau > 0$ (the protocol's "heterogeneous null").
- **Working parameter.** A parameter of the analysis model that need
  not equal any parameter of the data-generating process; here $h$
  (`docs/57`, Section 3.2).

**Evidence labels.** *Quoted* (reproduced from a source), *inspected*
(read from code or stored output), *derived* (follows algebraically,
to first order where stated), *verified* (computed by running code;
scripts in Section 15), *argued* (reasoned, not derived in full or
computed). Results labeled verified come from a closed-form
computation unless the text says Monte Carlo.

## 1. Summary for the busy reader

**The question.** Paper 04 tests $\beta_D$ with $D_{bc,it}(h)$ in the
model. The project wants the half-life chosen from the data
(`docs/57`). For the interaction, `docs/57` found that a selection
blind to the tested term (F1A0) is approximately harmless. For the main
effect no such selection exists: $h$ enters the model only through the
tested column.

**Five results** (Sections 4 to 9).

1. **The coefficient of $D_{bc}(h)$ is a weighted average of the drug
   profile plus a null bias** (derived): its probability limit is
   $b_0(h) + \beta_D^{\text{full}}\,\bar\Psi_\omega(h; t_{1/2})$, where
   the weights come from the design, the time adjustment and the
   working covariance. Its scale changes with $h$; the model-implied
   exposure contrast $\hat\beta_D(h)\,\bar d(h)$ is the quantity
   comparable across $h$ and with $\beta_D^{\text{exp}}$.
2. **A fixed half-life does not by itself guarantee size for the main
   effect.** The null limit $b_0(h)$ is exactly 0 for every $h$ and any
   working covariance when the time adjustment contains visit effects
   (VF), or in the CO design (derived). Under F1 it is $+0.33$ points
   in the Hybrid design (SE 1.22 at $N = 70$), which gives size 0.058
   at $N = 70$ and 0.068 at $N = 150$ by the normal approximation;
   under S it is $-0.61$ (size 0.083) in the Hybrid design and $-2.23$
   (size 0.47) in OL+BDC (verified). This corrects the expectation,
   carried over from the interaction, that any fixed $h$ is right-sized.
3. **AIC choice of $h$ is the maximum-$|t|$ test with a one-test
   critical value** (derived): because the model without the drug term
   does not depend on $h$, $\ell_1(h) = \ell_0 + \mathrm{LR}(h)/2$, so
   AIC picks the $h$ with the largest likelihood-ratio statistic. Its
   size is $\Pr(\max_h |T(h)| > 1.96)$: 0.088 for F1 and 0.077 for VF in
   the Hybrid design, against a calibrated critical value
   $c_{\max} \approx 2.15$ (verified, closed form); the Monte Carlo
   check gives 0.093 and 0.087 (MCSE 0.012; Section 11).
4. **A randomization test is exact under the sharp null and can contain
   the selection,** and it is cheap if the working covariance is taken
   from the null model, whose design does not depend on the path labels
   (derived). It uses only the randomized between-path contrasts: at
   $f_m = 0.4$ its power is 0.39 whatever statistic it is applied to,
   against 0.63 for the model-based F1 test, whose extra power rests on
   within-path contrasts protected only by the time and expectancy
   model (verified, closed form). Under the weak null an unstudentized
   permutation test is not exact; in the Hybrid design the variance
   ratio is 1.027 at $\tau = 5.49$ (size 0.053), and with two equally
   allocated paths it is exactly 1 (derived, verified).
5. **The half-life is barely identified for the main effect at
   realistic effect sizes.** At $f_m = 0.4$ and $N = 70$ the expected
   ML deviance of F1 between $h = 0$ and $h = 2$ weeks is below 2 for
   every true half-life examined; at the published ceiling it reaches
   10 to 15 against a short truth (verified). The best-fitting $h$ is
   a working value that depends on the effect size (Section 9.6).

**Anchored against compounding carryover** matters more for the main
effect than for the interaction: at a true half-life of 1 week the
compounding recursion inflates $\beta_D^{\text{exp}}$ by 11% in the
Hybrid design and 14% in OL+BDC and F1's power by 0.047 and 0.068
(verified), where the interaction's closed-form power moved by at most
0.002 (quoted, paper 01, Appendix B.3).

**The worked example** (one trial, $N = 70$, $f_m = 0.4$, true
$t_{1/2} = 0.2$, truth $\beta_D^{\text{exp}} = -2.85$; Section 10):
F1 at $h = 0.2$ estimated an exposure contrast of $-2.61$ ($p = 0.031$);
paying for a data-chosen $h$ gave $p$ = 0.047 (multivariate normal),
0.049 (Davies) and 0.056 (exact randomization with the maximum inside),
while every test resting on the randomized contrast alone gave 0.11 to
0.19. The profile could not exclude any half-life. The Monte Carlo
check (Section 11) reproduced the closed-form moments, sizes and
powers within its error.

**Provisional recommendation** (Section 12). Keep a half-life fixed in
advance (0.2 week) in the primary analysis and add a randomization test
of a visit-effects statistic with the grid maximum inside it, as the
size-exact companion; report the profile of $h$ and the model-implied
exposure contrasts across the grid; do not use F1A or AIC averaging as
tests. Everything here is provisional.

## 2. Scope and relation to other documents

The protocol (Sections 3 to 5) fixes the data-generating process, the
designs and the analyses; its Appendix A reviews 38 Monte Carlo and
closed-form studies of main-effect tests in N-of-1 and single-case
designs, and this note does not repeat that review (Section 13
cross-references it). `docs/56` states F1 for the interaction and its
pre-exposure guard; `docs/57` treats the data-chosen half-life for the
interaction. The closed form here extends the compendium's separation
result (`compendium-mathematics.tex`, eq:bp-plim) and its probability
limit of the drug coefficient (eq:bp-fwl; `docs/42`, Section 4) from
working independence to the random-intercept CAR(1) working covariance,
and adds the decomposition by half-life. All numbers are for the
refined process with $c_{bm} = 0$ (no biomarker term), $\rho = 0.7$,
anchored carryover and equal expected path shares, unless a table says
otherwise. Nothing here is a stress test.

## 3. Setting

### 3.1 The data-generating process

For participant $i$ on path $p$, under the refined process with the
protocol's random drug slope (protocol, Sections 3.1 and 3.2),

$$
Y_{i0} = \mathrm{BL}_i, \qquad
Y_{it} = \mathrm{BL}_i - (TV_{it} + PB_{it} + BR_{it}) - v_i\,\Psi_{pt},
\quad t \ge 1,
$$

with $\mathbb E(BR_{it}) = f_m m_{BR}\,\Psi_{pt}$, so that

$$
\boldsymbol\mu_p = \boldsymbol\mu^0 + \beta_D^{\text{full}}\,\boldsymbol\Psi_p,
\qquad \beta_D^{\text{full}} = -f_m m_{BR},
\qquad \Sigma_p = \Sigma + \tau^2\,\boldsymbol\Psi_p\boldsymbol\Psi_p' .
$$

$\boldsymbol\mu^0 = \mathbb E(\mathrm{BL}) - \mathbb E(TV) - \mathbb E(PB)$ is the non-drug
mean. In all three designs it is the same on every path, because the
paths share the visit schedule and the expectancy pattern (inspected,
the `designs` list of `04-power-simulation.R`). $\Sigma$, the covariance
of the components without heterogeneity, is also common to all paths,
because the $BR$ standard deviation is 8 on and off drug (protocol,
Section 3.3). The profile $\boldsymbol\Psi_p$ is the Gompertz curve in time on
drug, $G(t_{od}; 1, 5, 0.42)$ rescaled to the unit ceiling, restarting
at re-exposure, with anchored decay $\Psi_{\text{last}}\,2^{-t_{sd}/t_{1/2}}$
off drug. On Hybrid path 1 at $t_{1/2} = 0.5$ week it is 0, 0.390,
0.839, 0.891, 0.927, 0.232, 0.058, 0.390, 0.002 at weeks 0 to 20
(verified, `58-common.R`).

**The sharp null holds when $f_m = 0$ and $\tau = 0$** (derived). Then
$\boldsymbol\mu_p = \boldsymbol\mu^0$ and $\Sigma_p = \Sigma$ for every
path, so the distribution of $\mathbf Y_i$ does not depend on the path:
the outcome vectors are exchangeable across path labels. This is the
property a randomization test needs (Section 9.3).

### 3.2 The analyses and the half-life

An analysis fits $\mathbf Y_i = X_{p(i)}(h)\boldsymbol\xi + \mathbf e_i$
by GLS with a working covariance $V$, where $X_p(h)$ contains the
non-drug columns (intercept, linear time, $\eta_t$; or visit
indicators) and the drug column $\mathbf x_p(h)$:

$$
x_{pt}(h) = D_{bc,pt}(h) =
\begin{cases} 1 & \text{on drug} \\
2^{-t_{sd,pt}/h} = r_1^{\,t_{sd,pt}} & \text{off drug after first exposure} \\
0 & \text{at baseline and before first exposure,}\end{cases}
$$

with $h = 0$ giving $D_{pt}$ (the guard of `docs/56`, Section 3.2). In
the Hybrid design the off-drug visits all follow exposure, at $t_{sd}$
of 1, 2 and 4 weeks on path 1; 1, 2 and 6 on path 2; 1, 2, 3 and 4 on
path 3; 1, 2, 3 and 7 on path 4 (inspected). $D_{bc}(h)$ is therefore
a polynomial in $r_1 = 2^{-1/h}$, the share of the effect the indicator
keeps one week after stopping: $r_1 = 0.031$ at $h = 0.2$, 0.25 at
$h = 0.5$, 0.5 at $h = 1$ and 0.71 at $h = 2$ (derived). Any
$h \le 0.15$ gives $r_1 < 0.01$, so all half-lives in the
pharmacological range of the protocol (0.1 to 0.2 week) give nearly
the binary indicator at weekly visits.

The analysis grid is $\mathcal T = \{0, 0.2, 0.5, 1, 2\}$ weeks: the
grid of `docs/57` with 0.25 replaced by the protocol's default
$t_{1/2}^a = 0.2$.

## 4. Estimands

### 4.1 Definitions and values

The protocol (Section 4.1) defines three candidates on the
$\beta_D$ scale:

$$
\beta_D^{\text{full}} = -f_m m_{BR}, \qquad
\beta_D^{\text{exp}} = \beta_D^{\text{full}}\;
  \frac1P\sum_p\Big(\bar\Psi_p^{\,\text{on}}
  - \bar\Psi_p^{\,\text{off}}\Big), \qquad
\beta_D^{\text{bl}} = \text{the same over blinded visits only,}
$$

where $\bar\Psi_p^{\,\text{on}}$ and $\bar\Psi_p^{\,\text{off}}$ are means of
$\boldsymbol\Psi_p$ over the post-baseline on-drug and off-drug visits. The
script reproduces the protocol's values at $t_{1/2} = 0$: Hybrid
$-7.22$, OL+BDC $-9.35$, CO $-6.70$ for $\beta_D^{\text{exp}}$, and
Hybrid $-7.56$ for $\beta_D^{\text{bl}}$ (verified against the
protocol's `estimand-check.R` values, quoted). Table 2 adds the
dependence on the true half-life and on the carryover form.

| Design | $t_{1/2}$ (weeks) | $\beta_D^{\text{exp}}$, anchored | $\beta_D^{\text{exp}}$, compounding | $\beta_D^{\text{bl}}$, anchored |
|------|------|------|------|------|
| Hybrid | 0 | $-7.22$ | $-7.22$ | $-7.56$ |
| Hybrid | 0.2 | $-7.13$ | $-7.13$ | $-7.47$ |
| Hybrid | 0.5 | $-6.29$ | $-6.44$ | $-6.63$ |
| Hybrid | 1 | $-4.82$ | $-5.34$ | $-5.16$ |
| Hybrid | 2 | $-2.95$ | $-3.77$ | $-3.29$ |
| OL+BDC | 0.2 | $-9.20$ | $-9.20$ | $-10.80$ |
| OL+BDC | 1 | $-5.70$ | $-6.47$ | $-7.30$ |
| CO | 0.2 | $-6.70$ | $-6.70$ | $-6.70$ |
| CO | 1 | $-6.43$ | $-6.47$ | $-6.43$ |

Table: Estimand values at the published ceiling ($f_m = 1$; every value scales with $f_m$), equal path weights (verified, `58-cf-estimands.csv`)

Two readings (derived from the table). In the pharmacological range
the estimands barely move ($-7.22$ to $-7.13$ in the Hybrid design),
because the effect is almost gone one week after stopping. At one week
or more carryover fills the off-drug visits and $\beta_D^{\text{exp}}$
shrinks by a third or more, which is a property of the estimand, not
a bias.

### 4.2 What the coefficient of $D_{bc}(h)$ estimates

**Proposition 1 (decomposition; derived).** Let
$B = \sum_p \varpi_p X_p'V^{-1}X_p$ and let
$\boldsymbol\omega_p(h)' = \mathbf e_D'B^{-1}X_p(h)'V^{-1}$ be the row
of the GLS projection that gives the drug coefficient on path $p$.
Then $\hat\beta_D(h) = N^{-1}\sum_i
\boldsymbol\omega_{p(i)}(h)'\mathbf Y_i + o_p(N^{-1/2})$ with the
weights at the limit of the working covariance, and

$$
\operatorname{plim}\hat\beta_D(h)
= \sum_p \varpi_p\,\boldsymbol\omega_p(h)'\boldsymbol\mu_p
= \underbrace{\sum_p \varpi_p\,
  \boldsymbol\omega_p(h)'\boldsymbol\mu^0}_{b_0(h)}
\;+\; \beta_D^{\text{full}}\,
  \underbrace{\sum_p \varpi_p\,\boldsymbol\omega_p(h)'\boldsymbol\Psi_p}
  _{\bar\Psi_\omega(h;\,t_{1/2})} .
$$

The weights reproduce the columns: $\sum_p \varpi_p
\boldsymbol\omega_p'\mathbf x_p(h) = 1$ and $\sum_p \varpi_p
\boldsymbol\omega_p'\mathbf z_p = 0$ for every other column
$\mathbf z_p$ (derived, from $B^{-1}B = I$). This is the compendium's
eq:bp-fwl generalized to any working covariance and any half-life; at
working independence with visit effects it reduces to that equation.

Three consequences (derived):

- **The coefficient is a weighted average of the drug profile, on the
  scale of $D_{bc}(h)$.** If the true profile lay in the span of
  $\mathbf x_p(h)$ and the non-drug columns, $\bar\Psi_\omega = 1$ and the
  coefficient would be $\beta_D^{\text{full}}$. Under Gompertz onset
  the on-drug profile varies from 0.39 to 0.93 across the Hybrid
  visits, so it never does, and $\bar\Psi_\omega$ depends on which visits
  the analysis weights.
- **The scale changes with $h$.** A larger $h$ codes the off-drug
  visits as partly exposed, so the coefficient is the effect "per unit
  of $D_{bc}(h)$" and grows in magnitude as $h$ grows, even when the
  data do not change (Table 4, column "plim").
- **The comparable quantity is the model-implied exposure contrast.**
  The fitted drug part of the mean is $\hat\beta_D(h)\,x_{pt}(h)$; its
  on-minus-off mean over the design's visits, the definition of
  $\beta_D^{\text{exp}}$ applied to the working model, is
  $$
  \hat\beta_D^{\text{exp}}(h) = \hat\beta_D(h)\,\bar d(h), \qquad
  \bar d(h) = 1 - \frac1P\sum_p \bar x_p^{\,\text{off}}(h),
  $$
  because $x = 1$ on drug. $\bar d(h)$ is 1 at $h = 0$, 0.991 at 0.2,
  0.906 at 0.5, 0.755 at 1 and 0.551 at 2 weeks in the Hybrid design
  (verified). Any comparison or average across $h$ should be made on
  this scale (Section 9.4).

### 4.3 Sharp and weak nulls

Every estimand is zero when $f_m = 0$, so size refers to one hypothesis
in every design (protocol, Section 4.1). Two nulls must nevertheless be
distinguished, as in the literature on series of N-of-1 trials [43,
44] and on randomization inference [28, 29]:

- **Sharp null**, $f_m = 0$ and $\tau = 0$: no participant has a drug
  effect. The outcome vectors are exchangeable across paths (Section
  3.1), so a randomization test is exact (Section 9.3).
- **Weak null**, $f_m = 0$ and $\tau > 0$: individual effects
  $-v_i\boldsymbol\Psi_p$ average to zero. The means are still common to all
  paths, but the covariances $\Sigma_p = \Sigma + \tau^2\boldsymbol\Psi_p
  \boldsymbol\Psi_p'$ differ by path, so exchangeability fails. A model-based
  or CR2 test of the mean is still asymptotically valid if its standard
  error includes the $\tau^2$ term; CR2 does, because it is computed
  from participant residuals (derived from the sandwich form, `docs/56`,
  Section 3.5). A randomization test is then valid only asymptotically
  and only for a suitably studentized statistic (Section 9.3).

## 5. The probability limit of the drug coefficient

### 5.1 The limit under a working covariance

**Proposition 2 (limit and variance; derived).** Under randomized
allocation with shares $\varpi_p$, participants' contributions
$X_{p(i)}'V^{-1}(\mathbf Y_i - X_{p(i)}\boldsymbol\xi^*)$ are
independent with mean zero at the limit $\boldsymbol\xi^*$, so

$$
\operatorname{Var}\hat\beta_D(h) \approx \frac1N\sum_p \varpi_p\,
\boldsymbol\omega_p(h)'\big(\Sigma_p + \mathbf r_p\mathbf r_p'\big)
\boldsymbol\omega_p(h),
\qquad \mathbf r_p = \boldsymbol\mu_p - X_p(h)\boldsymbol\xi^* ,
$$

the sandwich with the mean misfit $\mathbf r_p$ included; under fixed
path counts the $\mathbf r_p\mathbf r_p'$ term drops. CR2 estimates the
first form. The working covariance itself converges to the
random-intercept CAR(1) member closest to the truth in Kullback-Leibler
distance, which the script finds by minimizing
$\sum_p \varpi_p[\log|V| + \operatorname{tr}\{V^{-1}(\Sigma_p +
\mathbf r_p\mathbf r_p')\}]$, the expected ML criterion per participant,
with $\boldsymbol\xi^*$ recomputed for each candidate $V$ (derived;
REML has the same limit). In the Hybrid design under the null the
limit has $\rho_a = 0.646$ per week for F1 (verified).

The limit of F1 under the true covariance is close to that under the
working covariance, and the working-independence limit of `docs/42` is
not: F1's null bias is $+0.29$ under independence, $+0.18$ under the
true covariance and $+0.33$ under the random-intercept CAR(1) limit,
with standard errors 1.63, 1.21 and 1.22 at $N = 70$ (verified,
`58-cf-null.csv`). The working covariance, not only the time
adjustment, decides the weights.

### 5.2 Under the null: when a fixed half-life keeps its size

**Proposition 3 (null limit; derived).** If the non-drug mean
$\boldsymbol\mu^0$ lies in the span of the non-drug columns, and those
columns are the same on every path, then $b_0(h) = 0$ for every $h$ and
every working covariance. Proof: $\boldsymbol\mu^0 = Z\mathbf a$ for the
common non-drug block $Z$, and $\sum_p \varpi_p\boldsymbol\omega_p'Z = 0$
by Proposition 1. Visit effects satisfy the condition in every design
(they span all of $\mathbb R^9$ on the common schedule). Linear time
with or without $\eta_t$ does not, because the non-drug mean is a
sigmoid natural course plus an expectancy curve that halves at week 9
(`docs/42`, Section 2).

So, unlike the interaction, whose null limit is zero for any fixed $h$
because the biomarker block separates (compendium, eq:bp-plim), the
main effect is protected for any fixed $h$ only by the time adjustment
or by the design's symmetry. Table 3 gives the size this implies.

| Design | Analysis | $b_0$, $h = 0$ | $b_0$, $h = 0.2$ | $b_0$, $h = 1$ | $b_0$, $h = 2$ | SE, $h = 0.2$ | Size, $h = 0.2$ |
|------|------|------|------|------|------|------|------|
| Hybrid | S | $-0.58$ | $-0.61$ | $-1.12$ | $-1.63$ | 1.14 | 0.083 |
| Hybrid | F1 | $+0.32$ | $+0.33$ | $+0.25$ | $+0.13$ | 1.22 | 0.058 |
| Hybrid | VF | 0 | 0 | 0 | 0 | 1.66 | 0.050 |
| OL+BDC | S | $-2.16$ | $-2.23$ | $-3.29$ | $-3.92$ | 1.18 | 0.473 |
| OL+BDC | F1 | $-0.22$ | $-0.22$ | $-0.09$ | $+0.60$ | 1.42 | 0.053 |
| OL+BDC | VF | 0 | 0 | 0 | 0 | 2.35 | 0.050 |
| CO | S | $-0.26$ | $-0.26$ | $-0.33$ | $-0.43$ | 1.48 | 0.053 |
| CO | F1 | 0 | 0 | $-0.06$ | $-0.15$ | 1.59 | 0.050 |
| CO | VF | 0 | 0 | 0 | 0 | 1.59 | 0.050 |

Table: Null bias $b_0(h)$ (points), standard error and normal-approximation size at $N = 70$, random-intercept CAR(1) working covariance at its limit, $t_{1/2} = 0.2$ (verified, `58-cf-null.csv`)

Readings (verified unless stated):

- **F1 is nearly, not exactly, right-sized in the Hybrid design.** Its
  null bias of $+0.33$ points is a quarter of a standard error at
  $N = 70$; because the bias is fixed and the standard error falls as
  $N^{-1/2}$, the size grows with $N$: 0.054 at $N = 35$, 0.058 at 70,
  0.068 at 150 and 0.086 at 300 (derived from Table 3 by scaling the
  standard error). The protocol's size screen (pooled rate within
  2 MCSE of 0.05, no cell above 0.075) would pass F1 at $N \le 150$
  and fail it at large $N$.
- **S, linear time without $\eta_t$, is not usable for the main effect
  in OL+BDC** and is liberal in the Hybrid design, as the protocol
  predicted from `docs/40` (protocol, Section 5.1).
- **The null bias depends on $h$** under F1 and S, so a data-chosen
  $h$ selects among biased statistics (Section 9.1).
- **VF and CO are protected exactly.** In CO the two paths are mirror
  images and F1's columns absorb what remains; the small $b_0$ at
  $h \ge 1$ comes from the asymmetric guard on the placebo-first path.

### 5.3 Under the alternative: bias and power as $h$ departs from $t_{1/2}$

Table 4 gives the limit, the model-implied exposure contrast and the
power of F1 and VF in the Hybrid design at $f_m = 0.4$ ($m_{BR} = 4.39$
points), for true and analysis half-lives on the grid.

| True $t_{1/2}$ | $h$ | F1 plim | F1 $\hat\beta_D^{\text{exp}}$ limit | $\beta_D^{\text{exp}}$ | F1 SE | F1 power | VF $\hat\beta_D^{\text{exp}}$ limit | VF SE | VF power |
|------|------|------|------|------|------|------|------|------|------|
| 0 | 0 | $-2.81$ | $-2.81$ | $-2.89$ | 1.20 | 0.651 | $-2.88$ | 1.62 | 0.428 |
| 0 | 1 | $-3.26$ | $-2.46$ | $-2.89$ | 1.60 | 0.533 | $-2.18$ | 2.21 | 0.257 |
| 0.2 | 0 | $-2.74$ | $-2.74$ | $-2.85$ | 1.20 | 0.630 | $-2.80$ | 1.62 | 0.407 |
| 0.2 | 0.2 | $-2.78$ | $-2.76$ | $-2.85$ | 1.22 | 0.627 | $-2.81$ | 1.66 | 0.400 |
| 0.2 | 1 | $-3.23$ | $-2.44$ | $-2.85$ | 1.60 | 0.526 | $-2.14$ | 2.21 | 0.250 |
| 0.2 | 2 | $-3.62$ | $-2.00$ | $-2.85$ | 1.96 | 0.455 | $-1.62$ | 2.71 | 0.192 |
| 0.5 | 0 | $-2.24$ | $-2.24$ | $-2.52$ | 1.20 | 0.466 | $-2.22$ | 1.62 | 0.278 |
| 0.5 | 0.5 | $-2.59$ | $-2.35$ | $-2.52$ | 1.37 | 0.473 | $-2.24$ | 1.92 | 0.251 |
| 1 | 0 | $-1.59$ | $-1.59$ | $-1.93$ | 1.20 | 0.263 | $-1.63$ | 1.62 | 0.171 |
| 1 | 1 | $-2.45$ | $-1.85$ | $-1.93$ | 1.59 | 0.336 | $-1.59$ | 2.21 | 0.158 |
| 1 | 2 | $-3.10$ | $-1.71$ | $-1.93$ | 1.96 | 0.353 | $-1.32$ | 2.71 | 0.142 |
| 2 | 0 | $-0.91$ | $-0.91$ | $-1.18$ | 1.20 | 0.118 | $-1.08$ | 1.62 | 0.102 |
| 2 | 2 | $-2.41$ | $-1.33$ | $-1.18$ | 1.96 | 0.234 | $-0.98$ | 2.71 | 0.101 |

Table: Hybrid design, $f_m = 0.4$, $N = 70$: probability limit of the drug coefficient and of the model-implied exposure contrast, standard error and normal-approximation power for F1 and VF (verified, `58-cf-alt.csv`; all grid values in the file)

Readings (verified unless stated):

- **In the pharmacological range the half-life hardly matters.** At
  true $t_{1/2} = 0.2$, every $h \le 0.5$ gives F1 power between 0.590
  and 0.630; a wrong $h = 1$ costs 0.10 and $h = 2$ costs 0.17,
  because a long $h$ codes the off-drug visits as exposed and
  inflates the standard error (1.22 to 1.96).
- **Under long carryover the best $h$ is longer than the truth.** At
  true $t_{1/2} = 1$, F1's power is 0.263 at $h = 0$, 0.336 at $h = 1$
  and 0.353 at $h = 2$. Power peaks near the working value that best
  fits the misspecified mean, not at the pharmacological half-life
  (the working-parameter point of `docs/57`, Section 3.2; Section 9.6).
- **F1 targets $\beta_D^{\text{exp}}$ closely when $h$ is near the
  truth**: its exposure contrast limit is within 0.09 point (3%) of
  $\beta_D^{\text{exp}}$ at $t_{1/2} = h = 0.2$, including the null
  bias, and within 0.17 at $t_{1/2} = h = 0.5$. A wrong long $h$
  attenuates it by 14% to 30%.
- **VF is unbiased under the null but pays in precision** (SE 1.66
  against 1.22, power 0.40 against 0.63) and targets the randomized
  contrast, whose value is close to $\beta_D^{\text{exp}}$ in the
  Hybrid design and equal to $\beta_D^{\text{bl}}$ in OL+BDC
  ($-4.38$ at $f_m = 0.4$; verified).

In OL+BDC F1's exposure contrast overstates $\beta_D^{\text{exp}}$ by
21% at $t_{1/2} = h = 0.2$ ($-4.40$ against $-3.68$), because the four
open-label visits carry the rising Gompertz response at full expectancy
and a linear time term with $\eta_t$ cannot separate the two; in CO
the limit is within 10% of the estimand at every $h$ and true
half-life (verified,
`58-cf-alt.csv`).

## 6. Expectancy confounding and the half-life

The open-label visits are all on drug at $\eta_t = 1$, the blinded
visits have $\eta_t = 0.5$ and baseline has neither. Over the nine
rows with equal path weights the correlation of $\eta_t$ with $D_{it}$
is 0.588 (0.591 with $D_{bc}(0.2)$; 0.509 over the eight post-baseline
rows), reproducing the compendium's 0.588 (verified; compendium,
Section "The published Hybrid construct").

**How F1 identifies $\beta_D$** (derived from Proposition 1). F1's
weights are orthogonal to the intercept, linear time and $\eta_t$
columns. The open-label rows carry the only $\eta_t = 1$ values, so
their contribution to the drug coefficient is limited to what the
linear time and $\eta_t$ restrictions leave over. Table 5 splits the
weights by row group.

| Row group | S: share of $\sum\lvert\omega\rvert$ | F1: share | VF: share | F1: contribution to $b_0$ | F1: contribution to $\bar\Psi_\omega$ times $m_{BR}$ |
|------|------|------|------|------|------|
| baseline | 0.104 | 0.030 | 0.005 | $-0.52$ | 0 |
| open label, on drug | 0.149 | 0.097 | 0.007 | $+0.61$ | $+1.03$ |
| blinded, on drug | 0.351 | 0.466 | 0.488 | $-1.21$ | $-8.96$ |
| blinded, off drug | 0.396 | 0.406 | 0.501 | $+1.44$ | $+0.17$ |
| total | 1 | 1 | 1 | $+0.33$ | $-7.76$ |

Table: Weights of the drug coefficient by row group (Hybrid, $h = t_{1/2} = 0.2$, random-intercept CAR(1) limit), and F1's null bias and drug signal by row group; the null mean is centered, which changes no total (verified, `58-cf-rowgroups.csv`)

Readings (derived from the table, verified numbers):

- **Adding $\eta_t$ moves the weight to the blinded rows**, from 75%
  under S to 87% under F1; VF puts 99% there. The blinded on-off
  contrast at matched expectancy is what identifies $\beta_D$ once the
  expectancy step is modeled, as the protocol's estimand discussion
  intends (protocol, Section 4.1).
- **The open-label rows enter F1 with net negative weight** and reduce
  the drug signal by 1.03 points per unit ceiling: F1 uses them mainly
  to estimate the time and expectancy terms, and the residual Gompertz
  rise between weeks 4 and 8 then feeds back into the drug coefficient.
- **The null bias is a sum of four non-small parts that nearly cancel**
  ($-0.52$, $+0.61$, $-1.21$, $+1.44$). Its size is therefore fragile
  to the shape of the natural course and the placebo curve, which a
  real trial does not know (argued).

**The between-path information.** At weeks 10, 16 and 20 the paths
differ in drug state; under $h > 0$ they also differ at weeks 11 and 12,
where paths 1 and 2 are one week from stopping and paths 3 and 4 two
weeks (inspected). With visit effects only these between-path
differences identify $\beta_D$, and the plim depends only on the
between-path differences of the drug profile:
$\operatorname{plim}\hat\beta_D = \beta_D^{\text{full}}
\sum_p\varpi_p\boldsymbol\omega_p'(\boldsymbol\Psi_p - \bar{\boldsymbol\Psi})$,
because $\sum_p\varpi_p\boldsymbol\omega_p = 0$ under visit effects
(derived). A positive $h$ therefore lets VF borrow weeks 11 and 12, at
the price of assuming the carryover shape there; this is how $h$
interacts with the expectancy protection (argued). The randomization
test of Section 9.3 uses exactly this between-path information whatever
statistic it is applied to.

## 7. Anchored against compounding carryover

**Why the form matters more for the main effect** (derived). For the
main effect the off-drug means are part of the signal: by Proposition
1 the limit is $\sum_p\varpi_p\boldsymbol\omega_p'\boldsymbol\mu_p$, and
the compounding recursion of the published process lowers the off-drug
drug response from the second off-drug visit of each run onward
(`docs/53`, Section 1), which raises the on-minus-off contrast directly.
For the interaction the slope depends on the coupling $g_t$, not on the
mean, and the carryover of the mean enters only through the path-mean
term of the pooled contrast, under 0.7% of the residual variance in the
Hybrid design (quoted, `docs/37`, Section 4.4); paper 01 found the
closed-form power changed by at most 0.002 (quoted, Appendix B.3).

| Design | $t_{1/2}$ | $\beta_D^{\text{exp}}$, anchored | compounding | F1 power, anchored | compounding | VF power, anchored | compounding |
|------|------|------|------|------|------|------|------|
| Hybrid | 0.2 | $-2.85$ | $-2.85$ | 0.627 | 0.627 | 0.400 | 0.400 |
| Hybrid | 0.5 | $-2.52$ | $-2.58$ | 0.473 | 0.487 | 0.251 | 0.249 |
| Hybrid | 1 | $-1.93$ | $-2.14$ | 0.336 | 0.383 | 0.158 | 0.158 |
| Hybrid | 2 | $-1.18$ | $-1.51$ | 0.234 | 0.302 | 0.101 | 0.112 |
| OL+BDC | 1 | $-2.28$ | $-2.59$ | 0.490 | 0.558 | 0.129 | 0.125 |
| OL+BDC | 2 | $-1.28$ | $-1.71$ | 0.180 | 0.258 | 0.075 | 0.075 |
| CO | 1 | $-2.57$ | $-2.59$ | 0.329 | 0.332 | 0.335 | 0.338 |

Table: Anchored against compounding carryover at $f_m = 0.4$, $N = 70$, $h = t_{1/2}$ (verified, `58-cf-anchored-recursive.csv`)

Readings (verified): in the pharmacological range the two forms agree
to three decimals. At $t_{1/2} = 1$ week the compounding form inflates
$\beta_D^{\text{exp}}$ by 11% (Hybrid) and 14% (OL+BDC) and F1's power
by 0.047 and 0.068; at 2 weeks by 28% and 34%, and 0.068 and 0.078.
VF changes little, because the second off-drug visits are at weeks 11
and 12, where VF's weights are smaller. CO has a single off-drug run of
four visits 2.5 weeks apart and is barely affected. Any main-effect
result from the stored stress test, which used the compounding form,
should not be carried into paper 04 at $t_{1/2} \ge 1$ week (argued).

## 8. Closed-form power for the main effect

**Proposition 4 (power; derived).** For a linear analysis with
weights $\boldsymbol\omega_p(h)$, by Propositions 1 and 2 the
statistic $\hat\beta_D(h)/\mathrm{SE}$ is approximately normal with
noncentrality

$$
\Lambda(h) = \sqrt N\;
\frac{b_0(h) + \beta_D^{\text{full}}\,\bar\Psi_\omega(h; t_{1/2})}
{\Big[\sum_p\varpi_p\boldsymbol\omega_p'(\Sigma + \mathbf r_p
\mathbf r_p')\boldsymbol\omega_p + \tau^2\sum_p\varpi_p
(\boldsymbol\omega_p'\boldsymbol\Psi_p)^2\Big]^{1/2}},
\qquad
\pi \approx \Phi(\Lambda - z_{0.975}) + \Phi(-\Lambda - z_{0.975}).
$$

This is the main-effect analog of `docs/37`'s E9 development: E9 is a
fixed per-participant contrast $\mathbf a$ of the outcome regressed on
the biomarker; here the contrast is $\boldsymbol\omega_p$, its mean is
the exposure-weighted drug profile instead of the coupling gap, and
heterogeneity enters as $\tau^2(\boldsymbol\omega_p'\boldsymbol\Psi_p)^2$,
which is the protocol's $\tau^2(\mathbf a'\boldsymbol\Psi)^2$ (protocol,
Section 5.1). For VF under working independence, $\boldsymbol\omega_p$
is the FWL contrast $D_{pt} - \bar D_t$ of the compendium's eq:bp-fwl
(the protocol's summary measure PR), and the moments are exact given
the path counts. The normal reference is adequate here: the CR2
Satterthwaite degrees of freedom for the drug coefficient were 64 to 69
at $N = 70$ in the worked example (Section 10), against 21 to 30 for
the interaction (`docs/56`, Section 3.5), because every participant
contributes an exposure contrast, not only those with extreme biomarker
values (verified for the example; argued in general).

| Design | Analysis | $\rho$ | $\tau$ = 0, $N$ = 35 | 70 | 150 | $\tau$ = 5.49, $N$ = 70 | SE ($\tau = 0$, $N = 70$) |
|------|------|------|------|------|------|------|------|
| Hybrid | F1 | 0.5 | 0.248 | 0.439 | 0.753 | 0.408 | 1.48 |
| Hybrid | F1 | 0.7 | 0.365 | 0.627 | 0.917 | 0.568 | 1.22 |
| Hybrid | F1 | 0.8 | 0.507 | 0.798 | 0.983 | 0.717 | 1.02 |
| Hybrid | VF | 0.7 | 0.226 | 0.400 | 0.704 | 0.379 | 1.66 |
| OL+BDC | F1 | 0.7 | 0.602 | 0.881 | 0.996 | 0.817 | 1.42 |
| OL+BDC | VF | 0.7 | 0.261 | 0.462 | 0.779 | 0.435 | 2.35 |
| CO | F1 | 0.7 | 0.204 | 0.358 | 0.646 | 0.340 | 1.59 |
| CO | VF | 0.7 | 0.205 | 0.360 | 0.650 | 0.342 | 1.59 |

Table: Normal-approximation power at $f_m = 0.4$, $h = t_{1/2} = 0.2$, by design, analysis, serial correlation, heterogeneity and $N$ (verified, `58-cf-power-grid.csv`, which also has $\tau = 2.75$ and every $\rho$)

Readings (verified):

- **Power rises with $\rho$** in the N-of-1 designs (Hybrid F1 0.44,
  0.63, 0.80 at $\rho$ = 0.5, 0.7, 0.8), because a stronger serial
  correlation makes nearby on and off visits more alike apart from the
  drug, so within-participant contrasts are less noisy.
- **Heterogeneity costs little at these values**: $\tau = 5.49$
  points (half the published ceiling) lowers Hybrid F1 power from
  0.627 to 0.568 at $N = 70$, because $\boldsymbol\omega_p'\boldsymbol\Psi_p$
  is about 0.7, so it adds $0.5\,\tau^2$ to a contrast variance of
  about 100.
- **OL+BDC F1's high power is partly its bias** away from zero
  (Section 5.3), and CO's F1 and VF coincide because its time terms
  are almost aliased with the visit effects.
- **The half-life enters through $\bar\Psi_\omega$ and the SE only**
  (Table 4). With $h$ fixed the closed form is a function of $h$,
  $t_{1/2}$, $\rho$, $\tau^2$ and $N$, all tabulated in the scripts'
  outputs.

Figure 1 summarizes Sections 5 to 8 and Section 9.3.

![Closed-form results, refined process, $N = 70$. A: limit of the model-implied exposure contrast relative to $\beta_D^{\text{exp}}$ as the analysis half-life $h$ varies, for five true half-lives, F1 and VF, Hybrid, $f_m = 0.4$. B: normal-approximation power in the same cells. C: size under the sharp null implied by the null bias $b_0(h)$, by analysis, Hybrid and OL+BDC. D: power of the model-based test and of the randomization test (Section 9.3) as the drug effect grows, by design, F1 and VF at $h = 0.2$.](figures/58-fig1-closed-form.png)

## 9. Choosing the half-life from the data

### 9.1 The half-life is not identified under the null

When $\beta_D = 0$ the model $\mathbf Y_i = X_p(h)\boldsymbol\xi +
\mathbf e_i$ does not depend on $h$, because $h$ enters only through
the column whose coefficient is zero. The half-life is a nuisance
parameter present only under the alternative, the problem of Davies
[19, 20], Andrews and Ploberger [21] and Hansen [22]: the score and
likelihood-ratio statistics for $\beta_D$ are functions of $h$, and
their null distributions are those of a process indexed by $h$, not of
a single $\chi^2_1$ variable.

**Proposition 5 (AIC selection with the tested term; derived).** Every
candidate model on $\mathcal T$ has the same parameters, so AIC
differences are log-likelihood differences. With $\ell_0$ the maximized
ML log-likelihood of the model without any drug term,

$$
\ell_1(h) = \ell_0 + \tfrac12\,\mathrm{LR}(h), \qquad
\hat h_{\text{AIC}} = \arg\max_{h\in\mathcal T}\mathrm{LR}(h),
$$

and $\ell_0$ does not depend on $h$. With the working covariance held
fixed, $\mathrm{LR}(h) = \hat\beta_D(h)^2/\widehat{\mathrm{Var}}_{\text{mb}}
\{\hat\beta_D(h)\}$ exactly, the square of the model-based Wald
statistic; with the covariance re-estimated it holds to first order.
F1A therefore selects the $h$ with the largest model-based $|t|$ and
refers the CR2 $t$ at that $h$ to a single-test critical value: to
first order it is the max-$|t|$ test of Section 9.2 with critical value
$t_{0.975}$, and its size is $\Pr(\max_{h}|T(h)| > t_{0.975}) > \alpha$.

This is the "flat $\ell_0$" limiting case of `docs/57`, Section 3.4,
holding exactly rather than approximately. The remedy of `docs/57`,
choosing $h$ on the model without the tested term (F1A0), has no analog:
removing $D_{bc}$ removes every trace of $h$, so any likelihood-based
choice of $h$ must look at the tested term. A data-chosen half-life for
the main effect can be valid only if the test is calibrated for the
choice (Sections 9.2, 9.3, 9.5), if $h$ is chosen on data the test does
not use (sample splitting, `docs/57`, Section 4.3), or if the choice
does not use the outcome (pharmacology, fixed in advance).

### 9.2 The maximum statistic over the grid

**Joint null distribution** (derived). With the working covariance at
its limit for each $h$, the estimates are asymptotically jointly normal
with

$$
\operatorname{Cov}\{\hat\beta_D(h), \hat\beta_D(h')\} \approx
\frac1N\sum_p\varpi_p\,\boldsymbol\omega_p(h)'\big(\Sigma_p +
\mathbf r_p(h)\mathbf r_p(h')'\big)\boldsymbol\omega_p(h'),
$$

so the correlation of $T(h)$ and $T(h')$ is determined by the overlap
of the weight vectors, that is by how similar the regressors
$D_{bc}(h)$ and $D_{bc}(h')$ are after the other columns are partialled
out, measured in the metric of the true covariance. In a trial the same
quantity is estimated from the participants' influence contributions
$U_i(h) = \mathbf e_D'\hat B_h^{-1}X_{p(i)}(h)'\hat V_h^{-1}\hat{\mathbf
e}_i(h)$, as $\sum_i U_i(h)U_i(h')$, the cluster-robust covariance of
the stacked estimates (the multiplier construction of Hansen [22]).
The level-$\alpha$ critical value is the $c_{\max}$ with
$\Pr(\max_h|Z_h| \le c_{\max}) = 1 - \alpha$ for $Z$ multivariate normal
with that correlation (`mvtnorm::qmvnorm`); a multivariate $t$ with the
smallest CR2 degrees of freedom is the small-sample variant.

| Design | Analysis | $\mathrm{cor}\{T(0), T(2)\}$ | $c_{\max}$ | Size of F1A-type selection at 1.96 | Size at $c_{\max}$ with the null bias | Size at fixed $h = 0.2$ |
|------|------|------|------|------|------|------|
| Hybrid | F1 | 0.827 | 2.16 | 0.088 | 0.056 | 0.058 |
| Hybrid | VF | 0.848 | 2.14 | 0.077 | 0.051 | 0.050 |
| Hybrid | S | 0.867 | 2.14 | 0.185 | 0.138 | 0.083 |
| OL+BDC | F1 | 0.783 | 2.19 | 0.095 | 0.057 | 0.053 |
| OL+BDC | VF | 0.811 | 2.17 | 0.081 | 0.050 | 0.050 |
| CO | F1 | 0.976 | 2.04 | 0.061 | 0.051 | 0.050 |
| CO | VF | 0.991 | 2.01 | 0.056 | 0.050 | 0.050 |

Table: Joint null distribution of $T(h)$ over $\mathcal T$, $N = 70$, $t_{1/2} = 0.2$: correlation of the two ends of the grid (the smallest), the 95% critical value of $\max|T|$, and the normal-approximation size of three rules, including the null bias $b_0(h)$ (verified, closed form, `58-tg-maxt-summary.csv`; the multivariate normal integrals are quasi-Monte Carlo, with an error of about 0.005 in $c_{\max}$)

Readings (verified):

- **The five statistics are highly correlated** (0.83 or more in the
  Hybrid design, 0.98 in CO), so the maximum costs little: $c_{\max}$
  is 2.14 to 2.19 rather than the Bonferroni 2.58. In the Hybrid design
  the F1 correlations fall from 0.999 between $h = 0$ and 0.2 to 0.83
  between 0 and 2 weeks.
- **Selection by AIC with a one-test critical value is liberal**: 0.088
  for F1 and 0.077 for VF in the Hybrid design. Without the null bias,
  F1's figure would be 0.080; the selection, not the bias, is most of
  the excess.
- **Calibration repairs the selection but not the bias**: F1's
  calibrated max-$|t|$ test has size 0.056 in the Hybrid design, the
  same as F1 at a fixed $h$, because the null bias of Section 5.2
  remains. VF's calibrated test is at 0.051 (the residual 0.001 is the
  normal approximation of a maximum over five correlated statistics).

**Davies' bound for a continuous range** (derived; formula as
implemented in `segmented::davies.test`, inspected, citing Davies [20]).
For a unit-variance Gaussian process $Z(h)$ the upcrossing argument
gives $\Pr(\sup_h Z(h) > u) \le \Phi(-u) + \mathcal V e^{-u^2/2}/
\sqrt{8\pi}$, with $\mathcal V = \sum_j|Z(h_j) - Z(h_{j-1})|$ the
observed total variation on a fine grid; for $|Z|$ the union of the
two tails gives $2\Phi(-u) + \mathcal V e^{-u^2/2}/\sqrt{2\pi}$. The
expected value of the observed variation is $\sqrt{2/\pi}$ times the
arc length $\sum_j\arccos\mathrm{cor}\{Z(h_j), Z(h_{j-1})\}$, which turns
the bound into $2\Phi(-u) + (\text{arc length})\,e^{-u^2/2}/\pi$
(derived). Over $h \in [0, 2]$ weeks the arc length is 0.65 for F1 and
0.59 for VF in the Hybrid design; the bound puts the critical value at
2.168 and 2.151, against 2.165 and 2.150 from the multivariate normal
on the fine grid, and at 2.22 and 2.17 over $[0, 4]$ weeks
(verified, `58-tg-davies.csv`). The bound is nearly exact here because
the process is smooth and short, and it is a convenient one-line
calibration for a continuous $h$.

**Alternatives to the supremum.** Andrews and Ploberger [21] derive
statistics optimal for weighted average power over the nuisance
parameter, the average and the exponential average of the LR process;
with five highly correlated statistics the supremum and the averages
should differ little (argued), and the supremum is the one whose
rejection region a reader can see.

**What the maximum test reports.** A $p$-value against the null of no
drug effect for any exposure regressor on the grid, and the $h$ at which
$|T|$ is largest. It does not provide an estimate with known coverage;
Section 9.7 gives intervals that account for the maximum.

### 9.3 Randomization tests over path assignments

**Exactness under the sharp null** (derived; [23, 40]). Let the
statistic be any function $\mathcal S(\mathbf Y, \mathbf L)$ of the outcome
vectors and the path labels, including the whole select-then-test
procedure (F1A, the maximum over $\mathcal T$, AIC averaging). Under the
sharp null the outcome vectors are exchangeable across labels (Section
3.1). If $\mathbf L$ was drawn from a randomization distribution
$G$ and $\mathbf L^*$ is drawn from $G$ independently, then
$\mathcal S(\mathbf Y, \mathbf L)$ and
$\mathcal S(\mathbf Y, \mathbf L^*)$ have
the same distribution given $\mathbf Y$, and the $p$-value
$(1 + \#\{m: \mathcal S(\mathbf Y, \mathbf L^*_m) \ge \mathcal S(\mathbf Y,
\mathbf L)\})/(M + 1)$ is super-uniform for any number $M$ of
relabelings. Nothing about the working covariance, the time adjustment,
the half-life or the carryover form enters the argument.

**Which randomization distribution.** The reference set must be the
design's own randomization ("analyze as you randomize"; Rosenberger,
Uschner and Wang [41]).

- *Fixed path counts* (the published allocation, 18, 18, 17, 17):
  uniform over permutations of the labels.
- *The refined process's allocation* (each participant to a path
  uniformly at random, redrawn if a path receives fewer than 2): either
  re-draw labels by the same rule, which is exact unconditionally, or
  permute the observed labels, which is exact conditionally on the
  realized counts and hence also unconditionally (derived; the counts
  are ancillary under the sharp null). The conditional version is
  preferable because the precision of every statistic depends on the
  counts.
- *Stratified or blocked allocation*: permute within strata or
  blocks only.

**A cheap exact version** (derived). Refitting the mixed model for
every relabeling costs about 1 second per relabeling with the grid
(Section 10). It is avoidable. Fit the null model $Y \sim t + \eta_t$
once; its design is the same on every path, so its estimated covariance
$\hat V_0$ does not depend on the labels. A GLS statistic computed with
$\hat V_0$ is then a fixed function of $(\mathbf Y, \mathbf L)$,
and with fixed counts each relabeling changes only the path sums
$\sum_{i \in p}\mathbf Y_i$: $B$ depends on the counts alone. The
studentization by a CR0 sandwich with $\hat V_0$ is equally cheap. The
script `58-common.R` (`perm_setup`, `perm_stats`) computes 999
relabelings of ten statistics in about 2 seconds. The weights are those
of the null model rather than of the fitted alternative, which changes
the estimate little ($-2.699$ against $-2.700$ for F1 at $h = 0$ in the
pilot trial; verified).

**What a randomization test can use** (Proposition 6; derived). Write a
linear statistic with label-invariant weights as

$$
\hat\beta = \underbrace{\frac1N\sum_i\big(\boldsymbol\omega_{p(i)} -
\bar{\boldsymbol\omega}\big)'\mathbf Y_i}_{\hat C}
+ \underbrace{\bar{\boldsymbol\omega}'\,\frac1N\sum_i\mathbf Y_i}_{\hat K},
\qquad \bar{\boldsymbol\omega} = \sum_p\frac{N_p}{N}\boldsymbol\omega_p .
$$

Under permutation with fixed counts $\hat K$ does not change, so the
test is a test of $\hat C$ alone, whose mean is
$\beta_D^{\text{full}}\sum_p\varpi_p(\boldsymbol\omega_p -
\bar{\boldsymbol\omega})'(\boldsymbol\Psi_p - \bar{\boldsymbol\Psi})$: only the
between-path differences of the drug profile carry signal. For VF
$\bar{\boldsymbol\omega} = 0$ and $\hat K = 0$. For F1, $\hat K$
carries the within-path, time-modeled part of the estimate, including
the null bias $b_0$: at $f_m = 0.4$ in the Hybrid design 45% of F1's
limit is in $\hat K$ (verified, `58-cf-randomization-power.csv`).

Two consequences (verified, closed form, normal approximation):

- **The randomization test removes F1's null bias**, since the bias is
  in $\hat K$.
- **Its power is that of the randomized contrast, whatever statistic
  it is applied to**: 0.392, 0.391 and 0.393 for F1, VF and S in the
  Hybrid design at $f_m = 0.4$, against 0.627 for the model-based F1
  test and 0.400 for VF; 0.44 against 0.88 for F1 in OL+BDC; 0.35
  against 0.36 in CO, where the within-path information is already
  balanced by the crossover (Figure 1, panel D). The randomization test
  does not throw information away arbitrarily: the extra power of F1
  is exactly the information whose validity depends on modeling the
  natural course and the expectancy curve.

**The weak null** (Proposition 7; derived). With heterogeneity the
paths have covariances $\Sigma_p = \Sigma + \tau^2\boldsymbol\Psi_p\boldsymbol\Psi_p'$, and the permutation distribution of $\hat C$ has the variance of
the path mixture,
$N^{-1}\sum_p\varpi_p(\boldsymbol\omega_p - \bar{\boldsymbol\omega})'
\bar\Sigma(\boldsymbol\omega_p - \bar{\boldsymbol\omega})$ with
$\bar\Sigma = \sum_r\varpi_r\Sigma_r$, while its sampling variance uses
$\Sigma_p$ on path $p$. The difference is

$$
\frac{\tau^2}{N}\Big[\sum_p\varpi_p\{(\boldsymbol\omega_p -
\bar{\boldsymbol\omega})'\boldsymbol\Psi_p\}^2 - \sum_p\varpi_p\sum_r
\varpi_r\{(\boldsymbol\omega_p - \bar{\boldsymbol\omega})'\boldsymbol\Psi_r\}^2\Big],
$$

positive when the heaviest-weighted paths have the largest exposure,
which makes the unstudentized test liberal. With two equally allocated
paths, $\boldsymbol\omega_2 - \bar{\boldsymbol\omega} =
-(\boldsymbol\omega_1 - \bar{\boldsymbol\omega})$ and the two sums are
equal, so the unstudentized test is asymptotically exact, the
multivariate form of the known result for the two-sample difference in
means with equal group sizes [28] (derived). In the Hybrid design the
variance ratio is 1.007 at $\tau = 2.75$ and 1.027 at $\tau = 5.49$,
which raises the size of the unstudentized test to 0.051 and 0.053; in
OL+BDC and CO it is 1 (verified, `58-tg-weak-null-permutation.csv`).

**Studentization.** A permutation test of a studentized statistic is
asymptotically valid under the weak null when the studentizing variance
is consistent for the statistic's sampling variance and the statistic
is asymptotically linear: Janssen [30] for the Behrens-Fisher problem,
Chung and Romano [28] in general, DiCiccio and Romano [31] for
regression coefficients, and Wu and Ding [29] for Fisher randomization
tests of weak nulls in completely randomized, stratified and clustered
experiments, which they show are exact under the sharp null and
asymptotically conservative under the weak null (quoted in substance
from the full text read). Here the units are participants, the arms are
paths, and the statistic is a linear contrast of arm means of vector
outcomes, so the conditions are those of a multi-arm completely
randomized experiment with a cluster-robust studentization (argued).
Two cautions (derived):

- **Studentize the contrast, not the whole F1 estimate.** For F1,
  $t = (\hat C + \hat K)/\widehat{\mathrm{SE}}$ with $\hat K$ fixed and
  of order $b_0$; under the permutation distribution the standard error
  changes but $\hat K$ does not, so $\hat K/\widehat{\mathrm{SE}}$
  varies with the relabeling and the studentized F1 statistic is not
  pivotal unless $b_0 = 0$. VF, whose $\hat K$ is zero, has no such
  problem; it is the natural statistic for a randomization test.
- **The weak null is not the sharp null with noise.** Under the weak
  null some participants are helped and some harmed. A randomization
  test that rejects says the mean effect is not zero only through the
  asymptotic argument; under the sharp null it says that someone
  responded, which is a weaker conclusion than paper 04 needs (argued;
  cf. Ding, Feller and Miratrix [32] on testing effect variation).

**Selection inside the relabelings.** Because exactness holds for any
function of $(\mathbf Y, \mathbf L)$, the AIC choice, the maximum
over $\mathcal T$ and the averaging can be run inside every relabeling.
The randomization $p$-value of F1A is then exact under the sharp null,
which answers the face-validity request (the data choose the half-life)
without the size cost of Proposition 5.

### 9.4 AIC-weight averaging and the estimand

Model averaging with smoothed AIC weights [36] over $\mathcal T$ has two
problems for the main effect (derived):

- **The averaged quantity must be common to the models.** The
  coefficients of $D_{bc}(h)$ for different $h$ are on different scales
  (Section 4.2): at true $t_{1/2} = 0.2$ the F1 limits are $-2.74$ at
  $h = 0$ and $-3.62$ at $h = 2$ while the model-implied exposure
  contrasts are $-2.74$ and $-2.00$ (Table 4). Averaging coefficients
  averages different estimands. The model-implied exposure contrast
  $\hat\beta_D^{\text{exp}}(h) = \hat\beta_D(h)\bar d(h)$ is common to
  all $h$ and targets $\beta_D^{\text{exp}}$; it is the quantity to
  average, with Buckland's unconditional standard error
  $\sum_h w_h\{\widehat{\mathrm{SE}}_h^2 + (\hat\beta_D^{\text{exp}}(h)
  - \bar\beta)^2\}^{1/2}$.
- **The weights see the tested term entirely.** By Proposition 5,
  $w_h \propto \exp\{\ell_1(h)\} \propto \exp\{\mathrm{LR}(h)/2\}$: the
  weights are a smooth version of F1A's selection, with no null-model
  alternative. Averaging is therefore a smoothed max-type statistic
  referred to a single-$t$ reference. Its size is between that of a
  fixed $h$ and that of F1A, and it is not known in closed form; it is
  checked in Section 11. Hjort and Claeskens [37] show that ignoring
  the selection step gives intervals that are too narrow, and the same
  applies to averaging with data-dependent weights unless the weights
  are accounted for.

Averaging is a reasonable way to report $\beta_D^{\text{exp}}$ across
half-lives; it is not a test with known size.

### 9.5 Parametric bootstrap under the null

For the main effect the null model has no half-life, so the bootstrap
is simpler than in `docs/57`, Section 4.5 (derived): fit $Y \sim t +
\eta_t$ with random intercept and CAR(1) once, simulate $\mathbf
Y^*_i$ from it for the observed participants, visits and paths, and
rerun the whole procedure (fixed $h$, F1A, maximum, averaging) on each
resample. It calibrates the selection exactly under the fitted null
model, and it does not remove a null bias that the null model
reproduces: the bootstrap data have the null model's linear-time and
$\eta_t$ mean, which F1 fits exactly, so F1's bootstrap null limit is
zero while its real null limit is $b_0 = +0.33$ (derived). The
bootstrap of F1 is therefore centered in the wrong place by the same
amount as the model-based test; the bootstrap of VF is not affected.
For a simulation study the warp-speed method of Giacomini, Politis and
White [38], one resample per replicate, estimates the bootstrap test's
rejection rate at about twice the cost of the plain procedure.

### 9.6 Joint estimation of the half-life

**Identification only under the alternative** (derived). In the
partially linear model with mean $\beta_D\,\mathbf x_p(h) +
Z\boldsymbol\xi_Z$, the information about $h$ is proportional to
$\beta_D^2$ [42]: $I_{hh} = \beta_D^2\sum_i \dot{\mathbf x}_{p(i)}'
V^{-1}\dot{\mathbf x}_{p(i)}$ after the other columns are partialled
out, with $\dot x_{pt} = \partial x_{pt}/\partial h = x_{pt}\ln 2\,
t_{sd,pt}/h^2$ off drug. Two singularities follow:

- **At $\beta_D = 0$** the information is zero (Section 9.1).
- **At $h \to 0$** the derivative vanishes faster than any power of
  $h$, because the smallest $t_{sd}$ is one week: every $h$ below about
  0.15 week gives the same regressor. The parameter $r_1 = 2^{-1/h}$
  removes this singularity: $x_{pt} = r_1^{\,t_{sd,pt}}$, a polynomial
  with $\partial x_{pt}/\partial r_1 = 1$ at the first off-drug visit
  of each run and 0 elsewhere when $r_1 = 0$. In $r_1$ the boundary at
  0 is regular, and by Self and Liang [39] the likelihood-ratio
  statistic for $r_1 = 0$ is asymptotically a 50:50 mixture of
  $\chi^2_0$ and $\chi^2_1$, so the one-sided profile interval uses the
  deviance cut 2.71 rather than 3.84 (derived; conditional on
  $\beta_D \ne 0$ and away from $|\beta_D|$ small).

**How much a trial says about $h$.** Table 9 gives the expected ML
deviance of F1, $N\sum_p\varpi_p[\mathbf r_p(h)'V^{-1}\mathbf r_p(h)]$
minus its minimum over $h$, a first-order approximation to the expected
profile deviance, and the best-fitting value $h^*$ (verified,
`58-cf-profile-deviance.csv`).

| True $t_{1/2}$ | $f_m$ | $h^*$ | dev. at $h = 0$ | 0.2 | 0.5 | 1 | 2 | $\infty$ |
|------|------|------|------|------|------|------|------|------|
| 0.2 | 0.4 | 0 | 0 | 0.03 | 0.40 | 1.03 | 1.69 | 4.25 |
| 0.2 | 1 | 0 | 0 | 0.24 | 3.19 | 8.52 | 14.01 | 31.47 |
| 0.5 | 0.4 | 0.5 | 0.10 | 0.07 | 0 | 0.12 | 0.41 | 2.63 |
| 0.5 | 1 | 0.35 | 0.36 | 0.22 | 0.20 | 2.10 | 5.20 | 20.65 |
| 1 | 0.4 | 2.4 | 0.71 | 0.67 | 0.40 | 0.13 | 0.01 | 1.46 |
| 1 | 1 | 1.2 | 2.69 | 2.42 | 0.90 | 0.03 | 0.31 | 11.14 |
| 2 | 1 | 6 | 5.15 | 4.98 | 3.70 | 2.11 | 0.78 | 4.93 |

Table: Expected ML deviance of F1 against the analysis half-life, Hybrid design, $N = 70$, and the best-fitting half-life $h^*$ on a fine grid (verified, first-order approximation)

Readings (verified):

- **At the protocol's alternative ($f_m = 0.4$) a trial of 70 says
  almost nothing about $h$**: no finite $h$ is more than 2 deviance
  units from the best when the truth is short, so every profile
  interval would cover 0 to 2 weeks.
- **At the published ceiling it separates short from long**: a truth
  of 0.2 week excludes $h \ge 1$ (deviance 8.5), as in the interaction
  example of `docs/57`, Section 5.
- **The best-fitting $h$ is a working value that depends on the effect
  size**: for a true one-week half-life, $h^* = 2.4$ at $f_m = 0.4$ and
  1.2 at $f_m = 1$; for a true 2 weeks, $h^* = 6$. Under F1 the null
  bias $b_0(h)$ and the Gompertz onset both pull the fit, and their
  relative weight changes with $f_m$ (derived from Proposition 1). A
  data-chosen $h$ must be reported as a working value, not as the
  drug's carryover half-life.

**Testing $\beta_D$ with $h$ estimated.** The profile likelihood of
$\beta_D$ maximizes over $h$, and at $\beta_D = 0$ the profile is flat
in $h$, so the profile LR statistic for $\beta_D = 0$ is
$\max_h\mathrm{LR}(h)$, the supremum statistic (derived from
Proposition 5). Referring it to $\chi^2_1$, as a profile-likelihood
interval for $\beta_D$ does, is the AIC problem again; the interval is
too short near zero. Away from zero, with $\beta_D$ well determined,
$h$ is identified and the usual profile intervals apply, with $h$
still estimating $h^*$. The Wald test at the joint estimate has the same
defect, and no CR2-type small-sample correction exists for the
nonlinear parameter (`docs/57`, Section 4.6).

### 9.7 Confidence intervals that account for the half-life choice

Three constructions (derived):

1. **Simultaneous intervals over the grid** (in the spirit of PoSI,
   Berk et al. [35]). With $c_{\max}$ of Section 9.2, the intervals
   $\hat\beta_D^{\text{exp}}(h) \pm c_{\max}\,\widehat{\mathrm{SE}}
   \{\hat\beta_D^{\text{exp}}(h)\}$ cover their own targets
   $\operatorname{plim}\hat\beta_D^{\text{exp}}(h)$ simultaneously with
   probability about $1 - \alpha$, so the interval at any data-chosen
   $h$, including the AIC choice, covers its target. Their union
   covers the target of whichever $h$ is reported. They do not remove
   the gap between a target and $\beta_D^{\text{exp}}$ (Table 4).
2. **Inversion of the randomization test** under constant-shape sharp
   hypotheses $H(b, h)$: everyone's drug effect is $b\,D_{bc,pt}(h)$.
   Under $H(b, h)$ the adjusted outcomes $\mathbf Y_i - b\,\mathbf
   x_{p(i)}(h)$ are exchangeable across labels, so the set of $b$ not
   rejected is an exact interval for $b$ given $h$, and its union over
   $h$, converted with $\bar d(h)$, is an interval for the exposure
   contrast that is valid if the true effect has one of the grid's
   shapes. The true Gompertz onset is not of that shape, so the
   coverage under the refined process is approximate (argued).
3. **Maximized $p$-values** (Berger and Boos [33]; Ding, Feller and
   Miratrix [32]): for a hypothesized $b \ne 0$, $h$ is identified, so
   the $p$-value can be maximized over a confidence set for $h$ of
   level $1 - \varsigma$ and $\varsigma$ added. This is valid for
   inverting tests of $b \ne 0$; at $b = 0$ it reduces to the
   randomization test, where $h$ plays no role.

## 10. Worked example: one trial with a drug effect

### 10.1 Set-up

One Hybrid trial was drawn from the refined process of the stress-test
driver (anchored carryover, randomized allocation, separable AR(1) at
$\rho = 0.7$), $N = 70$, path sizes 14, 18, 18 and 20, with a drug
effect $f_m = 0.4$ ($m_{BR} = 4.39$ points), true $t_{1/2} = 0.2$ week
and no heterogeneity. The estimands are $\beta_D^{\text{exp}} = -2.85$
and $\beta_D^{\text{bl}} = -2.99$ points. The seed (5801) was fixed
before the trial was generated, and the trial was not selected for any
property (inspected, `58-randomization-example.R`).

### 10.2 The statistic across the half-life

| $h$ (weeks) | F1 $\hat\beta_D$ | CR2 SE | $T(h)$ | df | $p$ | $\bar d(h)$ | F1 $\hat\beta_D^{\text{exp}}$ | VF $\hat\beta_D$ | VF $T(h)$ | VF $p$ |
|------|------|------|------|------|------|------|------|------|------|------|
| 0 | $-2.604$ | 1.181 | $-2.204$ | 68.9 | 0.031 | 1 | $-2.604$ | $-2.883$ | $-1.635$ | 0.107 |
| 0.2 | $-2.637$ | 1.198 | $-2.201$ | 68.9 | 0.031 | 0.991 | $-2.612$ | $-2.942$ | $-1.631$ | 0.108 |
| 0.5 | $-2.771$ | 1.326 | $-2.090$ | 68.8 | 0.040 | 0.907 | $-2.514$ | $-3.146$ | $-1.504$ | 0.137 |
| 1 | $-2.825$ | 1.538 | $-1.837$ | 68.9 | 0.070 | 0.757 | $-2.138$ | $-2.974$ | $-1.219$ | 0.227 |
| 2 | $-3.084$ | 1.890 | $-1.631$ | 68.9 | 0.107 | 0.554 | $-1.710$ | $-2.894$ | $-0.971$ | 0.335 |

Table: Example: F1 and VF at each grid half-life, REML with CR2 and Satterthwaite degrees of freedom; the model-implied exposure contrast is $\hat\beta_D(h)\bar d(h)$ (verified, `58-ex-grid.csv`; VF degrees of freedom 66.4 to 67.5)

![Example trial ($N = 70$, $f_m = 0.4$, true $t_{1/2} = 0.2$ week, refined process). A: CR2 $t$-statistic of the drug coefficient for F1 and VF over a fine grid of analysis half-lives, with $\pm 1.96$ (grey) and $-c_{\max}$ for each analysis (dashed). B: ML deviance of F1 against $h$, with the cuts 2.71 and 3.84. C: the randomization distribution of $\max_h|T(h)|$ over 999 permutations of the path labels with label-invariant weights, with the observed value (vertical line).](figures/58-fig2-example-trial.png)

Readings (verified):

- **The coefficient grows with $h$ while the exposure contrast
  shrinks**: from $-2.60$ to $-3.08$ against $-2.60$ to $-1.71$, the
  scale effect of Section 4.2. At $h \le 0.5$ the exposure contrast
  is within 0.34 point of the estimand.
- **The profile cannot locate $h$**: the ML estimate is on the boundary
  (0.020 week on a continuous search), the deviance never exceeds 2.1
  at any $h$ (1.40 at 1 week, 1.29 as $h \to \infty$), and the profile
  interval is all of $[0, \infty)$ at either cut (Figure 2B), as Table
  9 predicts for $f_m = 0.4$.
- **The likelihood ratio for the drug term is largest at $h = 0$**
  (4.27), so F1A chooses $h = 0$; in this trial that changes nothing,
  because $T(0)$ and $T(0.2)$ differ by 0.003.

### 10.3 What each procedure reports

| Procedure | $h$ used | $\hat\beta_D^{\text{exp}}$ | Statistic | $p$ |
|------|------|------|------|------|
| F1, fixed $h = 0.2$, CR2 | 0.2 | $-2.61$ | $t = -2.20$ | 0.031 |
| F1A, AIC with the tested term (invalid reference) | 0 | $-2.60$ | $t = -2.20$ | 0.031 |
| F1, AIC-weight average on the $\beta_D^{\text{exp}}$ scale | weighted mean 0.50 | $-2.43$ | $t = -2.01$ | 0.048 |
| F1, max-$\lvert t\rvert$, Bonferroni | 0 | | $2.20$ | 0.154 |
| F1, max-$\lvert t\rvert$, multivariate normal ($c_{\max} = 2.19$) | 0 | | $2.20$ | 0.047 |
| F1, max-$\lvert t\rvert$, Davies' bound on $[0, 4]$ | | | $2.20$ | 0.049 |
| Randomization, F1 at $h = 0.2$, null-model weights (999 permutations; 999 re-randomizations) | 0.2 | | $\lvert t\rvert = 2.13$ | 0.048; 0.050 |
| Randomization, F1 max-$\lvert t\rvert$ with selection inside (same) | | | 2.14 | 0.056; 0.064 |
| Randomization, F1 at $h = 0.2$ and F1A, with refits (49 permutations) | 0.2; 0 | | 2.20 | 0.08; 0.08 |
| Parametric bootstrap from the null model, F1 at 0.2; F1A (49 resamples) | | | 2.20 | 0.08; 0.10 |
| VF, fixed $h = 0.2$, CR2 | 0.2 | $-2.91$ | $t = -1.63$ | 0.108 |
| VF, max-$\lvert t\rvert$, multivariate normal ($c_{\max} = 2.15$); Davies | 0 | | 1.64 | 0.151; 0.186 |
| Randomization, VF at $h = 0.2$ (999; 999) | 0.2 | | 1.53 | 0.121; 0.113 |
| Randomization, VF max-$\lvert t\rvert$ with selection inside (999; 999) | | | 1.53 | 0.182; 0.160 |

Table: Example: what each procedure reports for the same trial (verified, `58-ex-*.csv`; the null-model-weight statistics use CR0 studentization, hence $|t| = 2.13$ rather than 2.20)

| Interval for $\beta_D^{\text{exp}}$ (truth $-2.85$) | Lower | Upper |
|------|------|------|
| F1, fixed $h = 0.2$, CR2 | $-4.98$ | $-0.24$ |
| F1, AIC choice, naive | $-4.96$ | $-0.25$ |
| F1, AIC choice with the max-$\lvert t\rvert$ constant | $-5.19$ | $-0.02$ |
| F1, AIC-weight average (Buckland) | $-4.84$ | $-0.02$ |
| F1, union over $\mathcal T$ with the max-$\lvert t\rvert$ constant | $-5.21$ | $0.58$ |
| VF, fixed $h = 0.2$, CR2 | $-6.48$ | $0.65$ |
| VF, union over $\mathcal T$ with the max-$\lvert t\rvert$ constant | $-6.94$ | $1.95$ |
| VF, randomization inversion at $h = 0.2$ (199 permutations) | $-6.69$ | $0.74$ |
| VF, randomization inversion, union over $\mathcal T$ | $-6.80$ | $1.80$ |

Table: Example: 95% intervals for the exposure contrast that do and do not account for the half-life choice (verified, `58-ex-intervals.csv`; inversion on a 0.25-point grid of $b$)

Readings (verified unless stated):

- **The verdict depends on what the analysis is allowed to assume**,
  not on the half-life. Every F1 procedure that relies on the time and
  expectancy model rejects at or near 0.05 (0.031 at a fixed $h$,
  0.047 to 0.056 once the choice of $h$ is paid for); every procedure
  that relies only on the randomized between-path contrast (VF, its
  calibrated maximum, the VF randomization tests) does not (0.108 to
  0.186). The half-life
  choice moves $p$ by 0.02 at most within either family.
- **The randomization test of the F1 statistic agrees with the
  calibrated F1 tests** (0.048 at a fixed $h$, 0.056 with the maximum
  inside), because its direction comes from the invariant within-path
  part (Section 11.3). The refit version with 49 permutations gives
  0.08, which with 49 relabelings has a resolution of 0.02 and a Monte
  Carlo error of about 0.04 (derived); it was cut from 199 because the
  machine was heavily loaded (Section 10.4), and it is not evidence
  against the cheap version.
- **The bootstrap of F1 is centered in the wrong place by design**
  (Section 9.5); with 49 resamples its $p$-values (0.08, 0.10) are too
  imprecise to compare.
- **Accounting for the choice widens the intervals modestly**: the F1
  interval at the AIC choice goes from $[-4.96, -0.25]$ to
  $[-5.19, -0.02]$ with the max-$|t|$ constant; the union over the
  grid crosses zero. All intervals cover the truth in this trial,
  which is one trial.

### 10.4 Run time

Grids of F1 and VF at five half-lives: 1.7 seconds; the fine profile
(54 half-lives, F1 and VF, REML and ML): 17 seconds; the label-invariant
randomization tests (999 permutations and 999 re-randomizations of ten
statistics): 1.9 seconds; the refit randomization test (49
permutations): 324 seconds and the null-model bootstrap (49): 25
seconds, the same work, so the difference measures the machine load
at the time; the randomization inversion (265 hypotheses, 199
permutations each): 8 seconds; 379 seconds in all on 7 workers
(verified, `58-ex-timing.csv`). A first run with 199 refits and
resamples was stopped after 15 minutes under a load average near 70.

## 11. A modest Monte Carlo check

### 11.1 Set-up

Four Hybrid cells of the refined process at $N = 70$ (anchored
carryover, randomized allocation, separable AR(1) at $\rho = 0.7$):
the sharp null ($f_m = 0$, $\tau = 0$, $t_{1/2} = 0.2$), a weak null
($f_m = 0$, $\tau = 5.49$), and alternatives at $f_m = 0.4$ with true
$t_{1/2} = 0.2$ and 1 week. Per trial: F1 and VF at every $h \in
\mathcal T$ (REML, CR2), the AIC choice with the tested term, AIC
averaging on the $\beta_D^{\text{exp}}$ scale, the maximum over the
grid with Bonferroni and with the influence-based multivariate normal
calibration, and the cheap exact randomization tests of Section 9.3
with 199 relabelings at fixed path counts (inspected,
`58-mc-check.R`). Seeds are per trial, $58{,}000{,}000 + 100{,}000
\times \text{cell} + r$.

**Replicates.** The run was planned at 1,000 null and 500 alternative
trials per cell and was cut to 600 and 300 while running, to respect
the time budget; 1,800 trials were analyzed (600, 600, 300, 300). The
MCSE of a rate near 0.05 is 0.0089 per null cell, and of a power near
0.5 at most 0.029. Differences from the closed form smaller than about
0.02 in size, or 0.06 in power, are not resolved. This is a check of
the closed form, not a size study.

### 11.2 Results

| Quantity | Sharp null | Weak null | Alt., $t_{1/2} = 0.2$ | Alt., $t_{1/2} = 1$ |
|------|------|------|------|------|
| $n_{\text{sim}}$ | 600 | 600 | 300 | 300 |
| F1, $h = 0.2$: mean estimate (closed-form limit) | 0.329 (0.327) | 0.388 (0.327) | $-2.907$ ($-2.783$) | $-1.652$ ($-1.629$) |
| F1, $h = 0.2$: empirical SD (closed-form SE) | 1.245 (1.219) | 1.329 (1.304) | 1.204 (1.219) | 1.283 (1.216) |
| F1, $h = 0.2$: mean CR2 SE; mean CR2 df | 1.216; 68.9 | 1.303; 68.9 | 1.216; 68.9 | 1.216; 68.9 |
| VF, $h = 0.2$: mean estimate (limit) | 0.011 (0) | 0.075 (0) | $-2.907$ ($-2.830$) | $-1.587$ ($-1.680$) |
| VF, $h = 0.2$: empirical SD (closed-form SE) | 1.714 (1.661) | 1.759 (1.714) | 1.701 (1.661) | 1.752 (1.656) |
| F1, fixed $h = 0.2$ | 0.065 (0.058) | 0.060 (0.057) | 0.653 (0.627) | 0.263 (0.269) |
| F1A, AIC with the tested term | 0.093 (0.088) | 0.085 | 0.707 (0.666) | 0.400 (0.407) |
| F1, AIC-weight average | 0.053 | 0.062 | 0.612 | 0.317 |
| F1, max-$\lvert t\rvert$, Bonferroni | 0.020 | 0.027 | 0.428 | 0.187 |
| F1, max-$\lvert t\rvert$, multivariate normal | 0.060 (0.056) | 0.072 | 0.639 (0.590) | 0.340 (0.331) |
| VF, fixed $h = 0.2$ | 0.058 (0.050) | 0.055 (0.050) | 0.400 (0.400) | 0.167 (0.172) |
| VFA, AIC with the tested term | 0.087 (0.077) | 0.083 | 0.420 (0.423) | 0.197 (0.218) |
| VF, AIC-weight average | 0.050 | 0.062 | 0.363 | 0.163 |
| VF, max-$\lvert t\rvert$, multivariate normal | 0.058 (0.051) | 0.068 | 0.383 (0.351) | 0.183 (0.169) |
| Randomization, F1 $\lvert t\rvert$ at $h = 0.2$ | 0.063 | 0.053 | 0.507 | 0.200 |
| Randomization, F1 $\lvert\hat\beta\rvert$ at $h = 0.2$ (unstudentized) | 0.067 | 0.067 | 0.507 | 0.203 |
| Randomization, F1 max-$\lvert t\rvert$ (selection inside) | 0.065 | 0.062 | 0.417 | 0.180 |
| Randomization, VF $\lvert t\rvert$ at $h = 0.2$ | 0.065 | 0.052 | 0.393 (0.391) | 0.163 |
| Randomization, VF $\lvert\hat\beta\rvert$ (unstudentized) | 0.060 | 0.067 (0.053) | 0.420 | 0.173 |
| Randomization, VF max-$\lvert t\rvert$ (selection inside) | 0.055 | 0.060 | 0.370 | 0.167 |

Table: Monte Carlo check, Hybrid, $N = 70$: moments of the drug coefficient and rejection rates at 0.05, with the closed-form value in parentheses where Sections 5 to 9 give one (verified, `58-mc-summary-*.csv`; MCSE 0.0089 for null rates near 0.05, at most 0.029 for power)

| Cell | Analysis | $h = 0$ | 0.2 | 0.5 | 1 | 2 |
|------|------|------|------|------|------|------|
| sharp null | F1A | 0.343 | 0.055 | 0.112 | 0.120 | 0.370 |
| sharp null | VFA | 0.363 | 0.038 | 0.108 | 0.105 | 0.385 |
| alt., $t_{1/2} = 0.2$ | F1A | 0.567 | 0.113 | 0.077 | 0.057 | 0.187 |
| alt., $t_{1/2} = 1$ | F1A | 0.237 | 0.027 | 0.097 | 0.120 | 0.520 |

Table: Share of trials in which AIC (with the tested term) chose each analysis half-life (verified, `58-mc-aic-choice.csv`)

### 11.3 Reading

- **The closed form holds** (verified). F1's mean estimate under the
  sharp null is 0.329 against the limit 0.327, and its empirical SD
  1.245 against 1.219 (2% larger); VF's mean is 0.011 against 0, and
  its SD 3% above the closed form. Power agrees within one to two MCSE
  in every cell checked (0.653 against 0.627, 0.400 against 0.400,
  0.263 against 0.269, 0.167 against 0.172). The mean CR2 standard
  errors match the closed-form SEs to within 0.3%, and the CR2
  Satterthwaite degrees of freedom average 69 for F1 and 66 to 67 for
  VF, so the normal reference is adequate (Section 8).
- **AIC selection with the tested term is liberal, as Proposition 5
  says**: 0.093 and 0.087 under the sharp null (closed form 0.088 and
  0.077), 0.085 and 0.083 under the weak null. Its choices pile up at
  the ends of the grid, 0 and 2 weeks (71% to 75% of null trials), the
  signature of a maximum over correlated statistics whose ends are the
  least correlated. Under the alternative with $t_{1/2} = 1$ week it
  chose 2 weeks in 52% of trials, which is $h^*$ territory (Table 9).
- **F1 at a fixed $h$ is close to, and slightly above, nominal**: 0.065
  (closed form 0.058) under the sharp null and 0.060 under the weak
  null, consistent with the null bias of Section 5.2 within MC error.
  VF at a fixed $h$ is at 0.058 and 0.055.
- **The calibrated maximum and averaging are near nominal; Bonferroni
  is conservative**: max-$|t|$ with the multivariate normal 0.060 (F1)
  and 0.058 (VF) under the sharp null, 0.072 and 0.068 under the weak
  null; averaging 0.053 and 0.050, 0.062 under the weak null;
  Bonferroni 0.020. The calibrated maximum costs 0.014 to 0.017 of
  power relative to the fixed $h$ when the fixed $h$ is right
  ($t_{1/2} = 0.2$) and gains 0.016 to 0.077 when it is wrong
  ($t_{1/2} = 1$).
- **The randomization tests are exact by construction, and the
  rejection rates are consistent with that within MC error**: 0.055
  to 0.067 under the sharp null for the eight randomization
  statistics, 1 to 2 MCSE above 0.05, on the same trials on which the
  model-based VF test rejected 0.058; with 199 relabelings the exact
  size is 0.050 (derived), so these rates indicate a chance excess
  shared by the trials, not a defect, but they cannot show exactness
  more sharply than 0.018.
- **The randomization power matches Proposition 6** for VF (0.393
  against 0.391). For the F1 statistic it is higher than predicted,
  0.507: because $\hat K$ is invariant and carries the within-path
  estimate, $|\hat C + \hat K|$ behaves as a one-sided test of $\hat C$
  in the direction of $\hat K$, at 0.05 one-sided rather than 0.025
  per side; a one-sided test with $\Lambda \approx 1.68$ has power
  0.51 (derived). This is legitimate, since the direction is chosen by
  a label-invariant quantity, and it is the one way the within-path
  information can enter an exact test.
- **The weak null** (single cell, $\tau = 5.49$): the unstudentized VF
  randomization test rejected 0.067 against a predicted 0.053, the
  studentized one 0.052; the difference is within 2 MCSE and does not
  resolve Proposition 7.

### 11.4 Run time

The 1,800 trials took 44,050 worker-seconds as recorded in the chunks,
but that figure is not a property of the procedure: the first launch
ran with a multithreaded BLAS, which oversubscribed the cores, and the
machine was heavily loaded by other processes (load average 50 to 70)
and appears to have been suspended for part of the run. A single trial
with one BLAS thread on an unloaded worker takes 1.7 seconds for the
whole per-trial procedure (verified), so the check would take about
$1{,}800 \times 1.7 / 7 \approx 7$ minutes on 7 workers (derived). The elapsed time
from the first launch to the last chunk was about two hours, well over
the 45-minute budget; the replicate counts were cut as stated above
and no chunk was rerun.

## 12. Comparison and provisional recommendation

### 12.1 Comparison table

"Sharp" and "weak" give validity under the two nulls of Section 4.3.
Size and power are for the Hybrid design at $N = 70$, sharp null and
base alternative ($f_m = 0.4$, $t_{1/2} = 0.2$): closed form (cf,
normal approximation) and the Monte Carlo check (MC, Section 11; MCSE
0.009 for size, 0.028 for power). "Face validity" asks whether the
data choose or check the half-life.

| Method | Sharp null | Weak null | Size, cf / MC | Power, cf / MC | Cost per trial | Estimand | Face validity | Status |
|------|------|------|------|------|------|------|------|------|
| F1, $h$ fixed at 0.2 | approximate: null bias $+0.33$, size grows with $N$ | as sharp; CR2 includes $\tau^2$ | 0.058 / 0.065 | 0.627 / 0.653 | 1 fit | $\beta_D^{\text{exp}}$ via $\bar d(h)$, 3% off | weak | cf, MC |
| VF, $h$ fixed at 0.2 | null limit exactly 0 | as sharp | 0.050 / 0.058 | 0.400 / 0.400 | 1 fit | randomized contrast | weak | cf, MC |
| F1A, AIC with the tested term | invalid (Proposition 5) | invalid | 0.088 / 0.093 | 0.666 / 0.707, inflated | 10 fits | at the chosen $h$ | good | not recommended |
| Max-$\lvert t\rvert$ over $\mathcal T$, F1, multivariate normal | approximate; keeps F1's bias | approximate | 0.056 / 0.060 | 0.590 / 0.639 | 5 fits | none single | moderate | cf, MC |
| Max-$\lvert t\rvert$ over $\mathcal T$, VF, multivariate normal | approximate | approximate | 0.051 / 0.058 | 0.351 / 0.383 | 5 fits | none single | moderate | cf, MC |
| Max-$\lvert t\rvert$, Bonferroni | conservative | conservative | below 0.02 / 0.020 | not computed / 0.428 | 5 fits | none single | moderate | MC |
| Max-$\lvert t\rvert$, Davies' bound, continuous $h$ | approximate, nearly exact here | approximate | bound tight to 0.003 in $c_{\max}$ | as multivariate normal | fine grid of fits | none single | moderate | cf, example |
| Randomization, VF statistic, maximum over $\mathcal T$ inside | exact | asymptotic if studentized | 0.050 by construction / 0.055 | 0.39 / 0.370 | 1 null fit; 2 s for 999 relabelings | randomized contrast | good | cf, MC, example |
| Randomization, F1 statistic | exact | studentized form not pivotal (Section 9.3) | 0.050 by construction / 0.063 | 0.39 two-sided theory / 0.507 | as above | F1's, direction from within-path part | good | MC, example |
| AIC-weight averaging on the $\beta_D^{\text{exp}}$ scale | unknown; weights are $\exp\{\mathrm{LR}/2\}$ | unknown | not derived / 0.053 | not derived / 0.612 | 10 fits | average of $\hat\beta_D^{\text{exp}}(h)$ | good | MC |
| Parametric bootstrap from the null model | calibrates the selection; not F1's bias | approximate | not run at scale | not run | $M$ times 10 fits | as the statistic | good | example ($M = 49$) |
| Joint nonlinear estimation of $h$ | invalid near $\beta_D = 0$ (supremum) | as sharp | not run | not run | 1 `nlme` fit | $h^*$, a working value | good | proposed, not recommended as a test |
| Sample splitting (choose $h$ on one half) | exact for the selection | as fixed $h$ | not run | about a $\sqrt 2$ loss of noncentrality (`docs/57`) | 2 half fits | as fixed $h$ | good | proposed |
| Bayesian marginalization over $h$ | unknown; needs calibration | unknown | not run | not run | MCMC | posterior of $\beta_D$ | good | proposed |
| Intervals: max-$\lvert t\rvert$ constant (PoSI-type) | approximate simultaneous coverage | approximate | | | 5 fits | target at each $h$ | good | example |
| Intervals: randomization inversion, union over $\mathcal T$ | exact if the effect has a grid shape | | | | 265 tests of 199 relabelings, 8 s | exposure contrast | good | example |

Table: Methods for testing the drug main effect with an unknown carryover half-life, Hybrid design, $N = 70$ (cf: closed form, verified; MC: Section 11, verified; validity entries derived or argued as in Sections 5 to 9)

### 12.2 Provisional recommendation for paper 04

Provisional; it rests on closed-form results, 1,800 Monte Carlo
trials in four Hybrid cells and one example trial.

1. **Primary test: F1 with $t_{1/2}^a = 0.2$ week fixed in advance**,
   CR2 with Satterthwaite degrees of freedom, as the protocol proposes.
   A fixed $h$ in the pharmacological range costs almost nothing when
   the truth is in that range (Table 4), and no data-chosen rule can
   gain much power where the data barely identify $h$ (Table 9). Its
   null bias ($+0.33$ points in the Hybrid design) must be reported and
   its size checked at $N = 150$, where the closed form predicts 0.068;
   if the stress test confirms an excess there, VF becomes the primary
   analysis.
2. **Size-exact companion with data-chosen half-life: a randomization
   test of the VF statistic with the maximum over $\mathcal T$ computed
   inside each relabeling**, fixed-count permutation, 999 relabelings,
   label-invariant null-model weights, CR0 studentization. It answers
   the face-validity request (the data choose $h$, and the reader sees
   at which $h$ the evidence is strongest) with exact size under the
   sharp null, and it states plainly how much of the evidence rests on
   the randomization alone.
3. **Report, in every trial**: $\hat\beta_D^{\text{exp}}(h)$ across
   $\mathcal T$ with CR2 intervals; the AIC-weighted average of these
   as a descriptive summary, not a test; the profile deviance of $h$
   in the $r_1$ parameterization with the 2.71 cut; and the statement
   that the best-fitting $h$ is a working value (Section 9.6).
4. **Not to be used as tests**: F1A or VFA with a nominal reference;
   S (linear time without $\eta_t$) in the Hybrid or OL+BDC design;
   the bootstrap of F1 as a correction for its null bias; joint
   estimation with a $\chi^2_1$ reference near zero.

### 12.3 Proposed edits to the protocol (not applied)

The protocol was not edited. Proposed edits, for the user's decision:

- **Section 4.2, estimand.** Add that every $D_{bc}(h)$ analysis
  reports the model-implied exposure contrast $\hat\beta_D(h)\bar d(h)$
  against $\beta_D^{\text{exp}}$, and that bias and coverage are
  computed on that scale, because the coefficient's scale changes with
  $h$ (Section 4.2).
- **Section 5.1, closed form.** Use the ML limit of the random-intercept
  CAR(1) working covariance rather than working independence, include
  the null bias $b_0(h)$ and the mean-misfit term of the sandwich, and
  replace "a nonzero limit is predicted for linear time without
  $\eta_t$" by the Table 3 values, which show a nonzero limit for F1 as
  well.
- **Section 5.2, analysis list. Add:** RV, the randomization test of
  VF with the maximum over $\mathcal T$ inside (999 fixed-count
  permutations, null-model weights, all cells); RF, the same with the
  F1 statistic; MX-F1 and MX-VF, the max-$|t|$ tests over $\mathcal T$
  calibrated by the influence-based multivariate normal; F1A as a
  labeled negative control (its size is predicted at 0.088). **Amend:**
  the "Misspecified half-life" note to use $\mathcal T = \{0, 0.2,
  0.5, 1, 2\}$ in every cell (the cost is five fits).
- **Section 6.1, measures.** Add the $h$ at which $|T|$ is largest and
  the AIC choice, by cell; the profile-interval width in $r_1$.
- **Section 6.4, decision rule.** Add a size screen at $N = 150$ for
  the analyses with a nonzero closed-form null bias (F1, S), since a
  fixed bias inflates size with $N$.
- **Section 8, open decisions. Revise item 12 (L1)** to the
  label-invariant randomization test above, run in every cell rather
  than a subset, since it costs about 2 seconds per trial. **Add:**
  (21) whether paper 04's primary test is model-based (F1, power 0.63
  at the base alternative) or randomization-based (RV, power about
  0.39), or both reported; (22) the half-life policy: fixed at 0.2
  week, or chosen inside a calibrated test (RV, MX); (23) whether the
  weak-null cells (block 2) judge the randomization tests in their
  studentized form only; (24) whether any main-effect result is quoted
  from runs with the compounding carryover recursion (Section 7
  recommends not, at $t_{1/2} \ge 1$ week).

## 13. What the literature says

The protocol's Appendix A reviews 38 simulation and closed-form studies
of main-effect tests in N-of-1, aggregated N-of-1, multi-crossover and
single-case designs; its Section A.3 tables compare their data-generating
processes and analyses. That review is not repeated. This section adds
the literature on the four issues of this note and says how the
development relates to it. The reading extent of each source is in
Section 16.

**Carryover in crossover trials.** Grizzle [7] showed that in the
two-period crossover the direct and residual (carryover) effects are
estimable only with between-participant information, which led to the
two-stage procedure: test carryover first, then analyze both periods or
period 1 only. Freeman [8] showed that the preliminary test is highly
correlated with the period-1 analysis, so the actual significance level
exceeds the nominal level even without carryover and "increasing
sample size only exacerbates the problem" (quoted, abstract); Senn [10,
12] and Fleiss [9] concluded that the procedure should not be used and
that crossover analysis rests on assumptions about carryover, which
washout must make plausible; Kabaila and Vicendese [13] found the same
failure for a four-period design. The two-stage procedure is a
data-chosen carryover model followed by a test that ignores the choice,
the main-effect analog of F1A, and Proposition 5 is its mixed-model
form: a selection over carryover regressors whose criterion is the
tested term's own likelihood ratio. The textbooks of Senn [5] and Jones
and Kenward [6] treat carryover as a nuisance that design should
remove; Kenward and Jones [11] show that baseline-adjusted analyses
need no assumption about the covariance structure, the stance CR2
takes here. The protocol's Appendix A.2.6 and A.2.10 (Cruz et al.;
Sturdevant and Lumley) cover carryover estimation in repeated-measures
crossovers.

**Pharmacodynamic and lag models.** The exponential decay of
$D_{bc}$ is the washout of a first-order effect compartment, the device
of Sheiner et al. [14] (rate constant fitted from concentration and
effect data), whose clinical use Holford and Sheiner [15] describe
(record only). Percha et al. [17] simulate first-order wash-in and
wash-out in single-person trials but do not model carryover in the
analysis; Liao et al. [16] estimate a Bayesian distributed lag of fixed
length chosen from the drug half-life; Gasparrini, Armstrong and
Kenward [18] choose lag bases by an information criterion with no
post-selection inference. Across this work the lag shape is either
fixed in advance or chosen and then treated as known; none calibrates a
test of the main effect for the choice. The working-parameter result of
Section 9.6 adds that, with Gompertz onset and a misspecified non-drug
mean, the best-fitting decay is not the pharmacological one.

**Nuisance parameters present only under the alternative.** Davies
[19, 20] formulated the problem and the upcrossing bound used in
Section 9.2; Andrews and Ploberger [21] derived weighted-average-power
optimal tests (average and exponential statistics); Hansen [22] gave a
simulation method for the $p$-value of the supremum based on the
conditional distribution of the score process, which is the
multiplier version of the influence-based correlation used here. Seber
and Wild [42] treat partially linear nonlinear models, in which the
nonlinear parameter is identified only when its linear coefficient is
not zero (record only). The boundary at $r_1 = 0$ is a case of Self and
Liang [39].

**Randomization tests.** Exactness under the sharp null for any
statistic is the classical randomization-model argument (Edgington and
Onghena [23]; Ernst [40]); Rosenberger, Uschner and Wang [41] stress
that the reference set must follow the actual, possibly restricted,
randomization rule. In single-case research, randomization tests are
designed around the randomization scheme (Onghena and Edgington [24];
Kratochwill and Levin [25]); Michiels et al. [26] wrap a multilevel
model in a randomization test and keep the nominal level where the
Kenward-Roger $F$ test fails, and Bouwmeester and Jongerling [27]
quantify the power cost of gradual onset. For the weak null, Janssen
[30], Chung and Romano [28], DiCiccio and Romano [31] and Wu and Ding
[29] show that permutation tests of unstudentized statistics need not
control the level, while studentized versions are asymptotically valid
(and exact under the sharp null); Ding, Feller and Miratrix [32] test
treatment-effect variation with the mean effect as a nuisance, using
the maximized $p$-value of Berger and Boos [33]. What this note adds is
specific to the Hybrid design: the label-invariant null-model weights
that make a mixed-model randomization test cheap (Section 9.3), the
decomposition showing that such a test uses only the between-path
contrast whatever statistic is permuted (Proposition 6), and the
variance-ratio formula for the weak null with vector outcomes
(Proposition 7).

**Post-selection inference and model averaging.** Leeb and Pötscher
[34] show that inference ignoring selection is invalid and that the
distribution of post-selection estimators cannot be estimated
uniformly; Berk et al. [35] obtain valid intervals after any selection
from a simultaneous max-$|t|$ constant, the construction of Section
9.7, item 1; Buckland, Burnham and Augustin [36] proposed AIC weights
and an unconditional variance, and Hjort and Claeskens [37] show that
naive intervals after selection or averaging are too narrow. For the
main effect the selection is over a single nuisance index, so the PoSI
constant is the $c_{\max}$ of five correlated statistics, near 2.16
rather than the Scheffé bound.

**Small-sample inference and computation.** CR2 with Satterthwaite
degrees of freedom is from Bell and McCaffrey [3] and Pustejovsky and
Tipton [4], who report size "nearly always less than or equal to the
nominal" when the degrees of freedom exceed 4 or 5 (quoted); the
warp-speed method is from Giacomini, Politis and White [38]. The
variance-component results for series of N-of-1 trials (Senn [43];
Araujo, Julious and Senn [44]) give the sharp and weak null
distinction its practical meaning: an analysis that ignores $\tau^2$ is
valid for the first and liberal for the second (protocol, Appendix
A.4).

## 14. Limitations

- **First-order and normal approximations.** The probability limits
  are exact for the stated working covariance; the variances, sizes and
  powers of Sections 5 to 9 are normal approximations at expected path
  shares. The Monte Carlo check (Section 11) covers four Hybrid cells
  at $N = 70$ only.
- **One data-generating process.** The refined process with
  $c_{bm} = 0$, exponential decay, no dropout and no measurement error.
  The null bias of F1 depends on the shapes of the natural course and
  the placebo curve (Section 6), which the refined process fixes at the
  published Gompertz parameters; another shape would give another
  $b_0$.
- **Heterogeneity is the protocol's proposal** ($v_i$ scaling the unit
  profile, independent of everything else), an open decision of the
  protocol (Section 8, item 9). The weak-null results depend on it.
- **The randomization tests use null-model weights.** The refit
  version was run on one trial only (Section 10). Their power under the
  weak null and with dropout is not evaluated; with dropout the cheap
  version needs per-pattern weights, and the refit version is the safe
  choice.
- **Not implemented:** joint nonlinear estimation of $h$ with a
  sandwich, a Bayesian analysis that marginalizes over $h$, sample
  splitting for the main effect, the Andrews-Ploberger average
  statistics, and the coverage of any interval in Section 9.7 beyond
  the single example.
- **The recommendation is provisional** and depends on the paper's
  estimand decision (protocol, Section 8, item 1): a randomization test
  answers a question about the randomized contrast, and paper 04 may
  prefer to report the model-based F1 test with its known small bias.

## 15. How to reproduce

From the repository root, in this order (the figure script reads the
outputs of the others):

```{.bash}
Rscript analysis/scripts/main-effect-half-life/58-plim-and-power.R
Rscript analysis/scripts/main-effect-half-life/58-tgrid-correlation.R
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript analysis/scripts/main-effect-half-life/58-mc-check.R --cores 7 --reps-null 600 --reps-alt 300
Rscript analysis/scripts/main-effect-half-life/58-randomization-example.R --cores 7
Rscript analysis/scripts/main-effect-half-life/58-figures.R
```

All source `58-common.R`, which evaluates the top-level definitions of
`analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R`
(everything before its pre-run checks, skipping the lines that create
output directories), as `57-common.R` does, so the data-generating
process is the driver's own, with anchored carryover and randomized
allocation. The drug effect is scaled by $f_m$ after the driver's
`build()`, which is exact because the carryover rules are linear in the
drug-response mean. Outputs go to
`analysis/data/quick-sim/main-effect-half-life/`; figures to
`docs/figures/58-*.png`.

- `58-common.R`: path moments, design matrices, the GLS limit and
  sandwich, the ML limit of the working covariance, estimands, the
  trial generator with heterogeneity, the main-effect fits with CR2,
  the influence-based correlation, the procedures, and the
  label-invariant randomization tests.
- `58-plim-and-power.R`: Tables 2 to 7 and 9, Figure 1 data
  (`58-cf-*.csv`).
- `58-tgrid-correlation.R`: Section 9.2 tables, Davies' bound, the
  weak-null variance ratio and the max-$|t|$ power (`58-tg-*.csv`).
- `58-mc-check.R`: Section 11 (`58-mc-*.csv`, chunks in
  `58-mc-chunks/`).
- `58-randomization-example.R`: Section 10 (`58-ex-*.csv`).
- `58-figures.R`: Figures 1 and 2.

**Run times** (verified, 4 performance plus 4 efficiency cores, 7
workers). 58-plim-and-power.R 16 seconds; 58-tgrid-correlation.R 15 seconds;
58-mc-check.R 1.7 seconds per trial unloaded, 1,800 trials, but about
two hours elapsed under the conditions of Section 11.4;
58-randomization-example.R 379 seconds with 49 refits and resamples
(Section 10.4); 58-figures.R a few seconds. The time budget of about
45 minutes was exceeded by the Monte Carlo check; its replicate
counts were cut from 1,000 and 500 to 600 and 300, and the example's
refit and bootstrap counts from 199 to 49.


## 16. References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health*. 2020;2:13. doi:10.3389/fdgth.2020.00013.
2. Morris TP, White IR, Crowther MJ. Using simulation studies to
   evaluate statistical methods. *Statistics in Medicine*.
   2019;38(11):2074-2102. doi:10.1002/sim.8086.
3. Bell RM, McCaffrey DF. Bias reduction in standard errors for linear
   regression with multi-stage samples. *Survey Methodology*.
   2002;28(2):169-181. (No DOI; Statistics Canada catalogue 12-001.)
4. Pustejovsky JE, Tipton E. Small-sample methods for cluster-robust
   variance estimation and hypothesis testing in fixed effects models.
   *Journal of Business and Economic Statistics*. 2018;36(4):672-683.
   doi:10.1080/07350015.2016.1247004.
5. Senn S. *Cross-over Trials in Clinical Research*. 2nd ed.
   Chichester: Wiley; 2002. doi:10.1002/0470854596.
6. Jones B, Kenward MG. *Design and Analysis of Cross-Over Trials*.
   3rd ed. Boca Raton: Chapman and Hall/CRC; 2014.
   doi:10.1201/b17537.
7. Grizzle JE. The two-period change-over design and its use in
   clinical trials. *Biometrics*. 1965;21(2):467-480.
   doi:10.2307/2528104.
8. Freeman PR. The performance of the two-stage analysis of
   two-treatment, two-period crossover trials. *Statistics in
   Medicine*. 1989;8(12):1421-1432. doi:10.1002/sim.4780081202.
9. Fleiss JL. A critique of recent research on the two-treatment
   crossover design. *Controlled Clinical Trials*. 1989;10(3):237-243.
   doi:10.1016/0197-2456(89)90065-2.
10. Senn S. The AB/BA crossover: past, present and future? *Statistical
    Methods in Medical Research*. 1994;3(4):303-324.
    doi:10.1177/096228029400300402.
11. Kenward MG, Jones B. The analysis of data from 2 x 2 cross-over
    trials with baseline measurements. *Statistics in Medicine*.
    1987;6(8):911-926. doi:10.1002/sim.4780060806.
12. Senn S. Viewpoint: do not resurrect the two-stage procedure.
    *Pharmaceutical Statistics*. 2022;21(4):808-814.
    doi:10.1002/pst.2224.
13. Kabaila P, Vicendese M. The performance of a two-stage analysis of
    ABAB/BABA crossover trials. *Biometrical Journal*.
    2012;54(3):361-369. doi:10.1002/bimj.201100201.
14. Sheiner LB, Stanski DR, Vozeh S, Miller RD, Ham J. Simultaneous
    modeling of pharmacokinetics and pharmacodynamics: application to
    d-tubocurarine. *Clinical Pharmacology and Therapeutics*.
    1979;25(3):358-371. doi:10.1002/cpt1979253358.
15. Holford NHG, Sheiner LB. Understanding the dose-effect
    relationship: clinical application of
    pharmacokinetic-pharmacodynamic models. *Clinical
    Pharmacokinetics*. 1981;6(6):429-453.
    doi:10.2165/00003088-198106060-00002.
16. Liao Z, Qian M, Kronish IM, Cheung YK. Analysis of N-of-1 trials
    using Bayesian distributed lag model with autocorrelated errors.
    *Statistics in Medicine*. 2023;42(13):2044-2060.
    doi:10.1002/sim.9676.
17. Percha B, Baskerville EB, Johnson M, Dudley JT, Zimmerman N.
    Designing robust N-of-1 studies for precision medicine: simulation
    study and design recommendations. *Journal of Medical Internet
    Research*. 2019;21(4):e12641. doi:10.2196/12641.
18. Gasparrini A, Armstrong B, Kenward MG. Distributed lag non-linear
    models. *Statistics in Medicine*. 2010;29(21):2224-2234.
    doi:10.1002/sim.3940.
19. Davies RB. Hypothesis testing when a nuisance parameter is present
    only under the alternative. *Biometrika*. 1977;64(2):247-254.
    doi:10.1093/biomet/64.2.247.
20. Davies RB. Hypothesis testing when a nuisance parameter is present
    only under the alternative. *Biometrika*. 1987;74(1):33-43.
    doi:10.1093/biomet/74.1.33.
21. Andrews DWK, Ploberger W. Optimal tests when a nuisance parameter
    is present only under the alternative. *Econometrica*.
    1994;62(6):1383-1414. doi:10.2307/2951753.
22. Hansen BE. Inference when a nuisance parameter is not identified
    under the null hypothesis. *Econometrica*. 1996;64(2):413-430.
    doi:10.2307/2171789.
23. Edgington ES, Onghena P. *Randomization Tests*. 4th ed. Boca
    Raton: Chapman and Hall/CRC; 2007. doi:10.1201/9781420011814.
24. Onghena P, Edgington ES. Customization of pain treatments:
    single-case design and analysis. *Clinical Journal of Pain*.
    2005;21(1):56-68. doi:10.1097/00002508-200501000-00007.
25. Kratochwill TR, Levin JR. Enhancing the scientific credibility of
    single-case intervention research: randomization to the rescue.
    *Psychological Methods*. 2010;15(2):124-144.
    doi:10.1037/a0017736.
26. Michiels B, Tanious R, De TK, Onghena P. A randomization test
    wrapper for synthesizing single-case experiments using multilevel
    models: a Monte Carlo simulation study. *Behavior Research
    Methods*. 2020;52(2):654-666. doi:10.3758/s13428-019-01266-6.
27. Bouwmeester S, Jongerling J. Power of a randomization test in a
    single case multiple baseline AB design. *PLoS ONE*.
    2020;15(2):e0228355. doi:10.1371/journal.pone.0228355.
28. Chung E, Romano JP. Exact and asymptotically robust permutation
    tests. *Annals of Statistics*. 2013;41(2):484-507.
    doi:10.1214/13-AOS1090.
29. Wu J, Ding P. Randomization tests for weak null hypotheses in
    randomized experiments. *Journal of the American Statistical
    Association*. 2021;116(536):1898-1913.
    doi:10.1080/01621459.2020.1750415.
30. Janssen A. Studentized permutation tests for non-i.i.d. hypotheses
    and the generalized Behrens-Fisher problem. *Statistics and
    Probability Letters*. 1997;36(1):9-21.
    doi:10.1016/S0167-7152(97)00043-6.
31. DiCiccio CJ, Romano JP. Robust permutation tests for correlation
    and regression coefficients. *Journal of the American Statistical
    Association*. 2017;112(519):1211-1220.
    doi:10.1080/01621459.2016.1202117.
32. Ding P, Feller A, Miratrix L. Randomization inference for
    treatment effect variation. *Journal of the Royal Statistical
    Society, Series B*. 2016;78(3):655-671. doi:10.1111/rssb.12124.
33. Berger RL, Boos DD. P values maximized over a confidence set for
    the nuisance parameter. *Journal of the American Statistical
    Association*. 1994;89(427):1012-1016.
    doi:10.1080/01621459.1994.10476836.
34. Leeb H, Pötscher BM. Model selection and inference: facts and
    fiction. *Econometric Theory*. 2005;21(1):21-59.
    doi:10.1017/S0266466605050036.
35. Berk R, Brown L, Buja A, Zhang K, Zhao L. Valid post-selection
    inference. *Annals of Statistics*. 2013;41(2):802-837.
    doi:10.1214/12-AOS1077.
36. Buckland ST, Burnham KP, Augustin NH. Model selection: an integral
    part of inference. *Biometrics*. 1997;53(2):603-618.
    doi:10.2307/2533961.
37. Hjort NL, Claeskens G. Frequentist model average estimators.
    *Journal of the American Statistical Association*.
    2003;98(464):879-899. doi:10.1198/016214503000000828.
38. Giacomini R, Politis DN, White H. A warp-speed method for
    conducting Monte Carlo experiments involving bootstrap estimators.
    *Econometric Theory*. 2013;29(3):567-589.
    doi:10.1017/S0266466612000655.
39. Self SG, Liang KY. Asymptotic properties of maximum likelihood
    estimators and likelihood ratio tests under nonstandard
    conditions. *Journal of the American Statistical Association*.
    1987;82(398):605-610. doi:10.1080/01621459.1987.10478472.
40. Ernst MD. Permutation methods: a basis for exact inference.
    *Statistical Science*. 2004;19(4):676-685.
    doi:10.1214/088342304000000396.
41. Rosenberger WF, Uschner D, Wang Y. Randomization: the forgotten
    component of the randomized clinical trial. *Statistics in
    Medicine*. 2019;38(1):1-12. doi:10.1002/sim.7901.
42. Seber GAF, Wild CJ. *Nonlinear Regression*. New York: Wiley; 1989.
    doi:10.1002/0471725315.
43. Senn S. Sample size considerations for n-of-1 trials. *Statistical
    Methods in Medical Research*. 2019;28(2):372-383.
    doi:10.1177/0962280217726801.
44. Araujo A, Julious S, Senn S. Understanding variation in sets of
    N-of-1 trials. *PLoS ONE*. 2016;11(12):e0167167.
    doi:10.1371/journal.pone.0167167.
45. Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS*. New
    York: Springer; 2000. doi:10.1007/b98882.

**How the references were checked** (2026-10-09). Bibliographic
details of 2 to 42 were checked against at least one of PubMed, the
Crossref record of the DOI (queried without any contact parameter), the
publisher page or the Zotero item; 1, 43, 44 and 45 are as verified in
the protocol and `docs/57`. Reading extent: full text, at least in
part, for 3, 4, 16, 17, 18, 27, 28, 29, 32, 35, 38, 40 and 41 (16, 17,
18 and 27 through a summarizing reader of the PMC page, so their
quotations may be paraphrased); abstract only for 7 to 14, 19 to 22,
24 to 26, 31, 33, 34, 36, 37 and 39; record only for 5, 6, 15, 23, 30
and 42. The Davies bound of Section 9.2 was
checked against its implementation in `segmented::davies.test`
(inspected), not against Davies' text. Statements attributed to a
source read only in abstract are limited to what the abstract says.

## Appendix A. Notation check

Every symbol used in this note, with its source. "Canonical" means
`analysis/report/NOTATION.md`, Part 1; "protocol" means the working
protocol's notation index; "defined here" means local to this note.

| Symbol | Source |
|------|------|
| $i$, $t$, $N$, $P$, $Y_{it}$, $\mathrm{BL}_i$, $TV$, $PB$, $BR$ | canonical |
| $\eta_t$ | canonical (design expectancy) |
| $D_{it}$, $D_{bc,it}$, $t_{od}$, $t_{sd}$, $t_{pb}$, $t_{1/2}$ | canonical |
| $t_{1/2}^a$ | canonical as used in `docs/56`, `docs/57` and the protocol |
| $\beta_D$, $\rho$, $\alpha$, $\pi$, $\sigma_{BR}$ | canonical |
| $w_t$, $c_1$ | canonical |
| $m_{BR}$, $\tau^2$ | canonical (Gompertz ceiling of $BR$; heterogeneity variance) |
| $\Psi_{it}$, $\boldsymbol\Psi_p$ | canonical (2026-10-09; unit drug-response profile, formerly $q_{it}$, $\mathbf q_p$) |
| $f_m$, $v_i$, $\varpi_p$, $p$ (path) | protocol; the path letter $p$ follows the protocol and `NOTATION.md`'s $q$, $p$ row, and coincides in letter with the $p$-value, which the text distinguishes by context |
| $\beta_D^{\text{full}}$, $\beta_D^{\text{exp}}$, $\beta_D^{\text{bl}}$ | protocol |
| $h$, $T(h)$, $\mathrm{LR}(h)$, $M$ | `docs/57` (local there) |
| $n_{\text{sim}}$, $\Pr$, $\mathbb E$ | canonical |
| $\ell_1(h)$, $\ell_0$ | `docs/57` notation, redefined here ($\ell_0$ is the model without any drug term) |
| $\sigma_u^2$, $\rho_a$ | canonical |
| $\Lambda$ | canonical (noncentrality) |
| $p(i)$, $\mathcal T$, $r_1$, $\boldsymbol\mu_p$, $\boldsymbol\mu^0$, $\Sigma_p$, $\Sigma$, $\bar\Sigma$ | defined here |
| $X_p(h)$, $\mathbf x_p(h)$, $Z$, $\boldsymbol\xi$, $\boldsymbol\xi^*$, $V$, $\hat V_0$, $B$ | defined here |
| $\boldsymbol\omega_p(h)$, $\bar{\boldsymbol\omega}$, $b_0(h)$, $\bar\Psi_\omega(h; t_{1/2})$, $\bar d(h)$, $\mathbf r_p$ | defined here |
| $\hat\beta_D^{\text{exp}}(h)$, $\hat C$, $\hat K$, $c_{\max}$, $\mathcal V$, $U_i(h)$, $h^*$, $I_{hh}$ | defined here |
| $\mathbf L$, $\mathbf L^*$ (label vectors), $G$, $\mathcal S(\cdot)$ (a statistic) | defined here |
| $\varsigma$ | defined here (Berger-Boos allowance; $\gamma$ is reserved) |
| $z_{0.975}$, $t_{0.975}$, $\Phi$, $\chi^2_1$ | standard |

Table: Every symbol used, with its source

(Notation, 2026-10-09: $R$ is now $n_{\text{sim}}$, $u_i(h)$ is now
$U_i(h)$, probability is $\Pr$ and expectation $\mathbb E$, following
`NOTATION.md`; entries that `NOTATION.md` now lists are marked
canonical. The unit drug-response profile, written $\mathbf q_p$ and
$q_{pt}$ in the protocol and here until 2026-10-09, is now
$\boldsymbol\Psi_p$ and $\Psi_{pt}$ in both, because $q$ is the
canonical path index; $\Psi$ is added to `NOTATION.md`. The protocol's
path letter $p$ is kept.)

Reserved symbols checked: $\gamma$, $\delta$, $\kappa$, $\nu$,
$\theta$ and $\lambda$ are not used; $\pi$ is used only for power.
`perl tools/notation-lint.pl` reports the file clean (Section 15).
