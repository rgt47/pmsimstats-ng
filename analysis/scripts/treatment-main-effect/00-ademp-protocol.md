# The Drug Main Effect in Aggregated N-of-1 Designs under the Refined Process: ADEMP Working Protocol for the Rebuild of Paper 04 {.unlisted .unnumbered}
*2026-10-08 15:57 PDT*

**Author.** pmsimstats team

**Manuscript.** `analysis/report/04-treatment-main-effect/report.Rmd`
(to be rebuilt).

**Protocol history.** This file replaces the working protocol of
2026-05-07, which is retained unchanged as
`analysis/scripts/treatment-main-effect/00-ademp-protocol-v1.md`
(byte-identical copy, checked with `cmp`). Version 1 describes the
study behind the current manuscript and its driver
`production-run.R`; it is kept as the record of that study. On
2026-10-09 this file was revised from the main-effect half-life note
`docs/58-main-effect-testing-and-carryover-half-life.md`, Section 12.3
(Revision 1 under "Changes to the working protocol"). Each revised
passage carries a dated note naming the docs/58 section it comes from.
The revision remains a working proposal, made before any result of
this study exists.

**Purpose.** This working protocol sets out a new simulation study of
the treatment (drug) main effect $\beta_D$ in aggregated N-of-1
designs, to replace the study of the current paper 04. It is a
planning document: no code has been written for it and nothing has
been simulated. Section 8 lists the decisions that must be settled
before code is written. The study is exploratory. Its criteria may be
revised as results come in, every revision will be recorded as an
amendment with its reason, and its findings should be confirmed in an
independent study. The structure follows ADEMP: aims, data-generating
mechanisms, estimands, methods and performance measures (Morris, White
and Crowther [1]).

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. Symbols marked *local*
are defined here and are not part of the canonical set.

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $i$, $t$ | participant and visit indices; $t = 0$ is baseline | canonical |
| $w_t$ | calendar week of visit $t$ | local (`docs/40`) |
| $N$ | total participants across all paths or arms | canonical |
| $P$ | number of randomization paths (4 Hybrid, 2 OL+BDC, 2 CO) | canonical |
| $Y_{it}$ | symptom score (CAPS points; code column `Sx`) | canonical |
| $\mathrm{BL}_i$ | baseline level | canonical |
| $BR_{it}$, $PB_{it}$, $TV_{it}$ | drug, placebo-belief and natural-course responses | canonical |
| $\eta_t$ | design expectancy: 1 open label, 0.5 blinded, 0 at baseline (code `De`) | canonical |
| $D_{it}$, $D_{bc,it}$ | binary drug state; exposure-decayed indicator, 0 before any exposure (code `Dbc`) | canonical |
| $t_{od}$, $t_{sd}$, $t_{pb}$ | time on drug, time since discontinuation, time under expectancy | canonical |
| $t_{1/2}$, $t_{1/2}^{a}$ | true carryover half-life; half-life assumed by an analysis (weeks) | canonical; local (`docs/45`) |
| $h$ | a candidate value of $t_{1/2}^{a}$ (weeks); $D_{bc,it}(h)$ is the indicator built with it | local (`docs/57`, `docs/58`) |
| $\mathcal T$ | analysis half-life grid $\{0, 0.2, 0.5, 1, 2\}$ weeks | local (`docs/58`) |
| $r_1$ | $2^{-1/h}$, the share of the drug effect that $D_{bc}(h)$ keeps one week after stopping | local (`docs/58`) |
| $G(t; m, d, r)$ | modified Gompertz curve with ceiling $m$, displacement $d$ and rate $r$ per week, $t$ in weeks (code `maxr`, `disp`, `rate`) | local (`docs/40`) |
| $m_{BR}$, $m_{TV}$, $m_{PB}$ | Gompertz ceilings of the three components (published 10.986, 6.506, 6.506) | local |
| $f_m$ | drug-effect multiplier, $m_{BR}/10.986$ | local |
| $\Psi_{it}$, $\boldsymbol\Psi_p$ | unit drug-response profile: the $BR$ mean of participant $i$'s path divided by $m_{BR}$, anchored carryover included; $\boldsymbol\Psi_p$ its vector over the visits of path $p$; $\bar\Psi_d$ the design's on-minus-off mean (renamed 2026-10-09 from $q_{it}$, $\mathbf q_p$, $\bar q_d$, because $q$ is the canonical path index; capital $\Psi$, distinct from the gating slope $\psi_g$ of paper 03 and from Araujo et al.'s $\psi^2$ quoted in Appendix A.2.1) | canonical |
| $v_i$ | participant-specific deviation of the drug-response ceiling | local |
| $\tau^2$ | variance of $v_i$, the participant-by-treatment variance ($\sigma_v^2$ of `docs/51`) | local |
| $\beta_D$ | drug main effect, the primary estimand | canonical |
| $\beta_D^{\text{full}}$, $\beta_D^{\text{exp}}$, $\beta_D^{\text{bl}}$ | the three candidate estimands of Section 4.1 | local |
| $\hat\beta_D^{\text{exp}}(h)$, $\bar d(h)$ | model-implied exposure contrast $\hat\beta_D(h)\,\bar d(h)$; on-minus-off mean of $D_{bc}(h)$ over post-baseline visits | local (`docs/58`) |
| $b_0(h)$ | null bias: probability limit of the drug coefficient at $f_m = 0$ | local (`docs/58`) |
| $T(h)$, $c_{\max}$ | CR2 $t$-statistic of the drug coefficient at $h$; 95% critical value of $\max_{h \in \mathcal T}\lvert T(h)\rvert$ | local (`docs/57`, `docs/58`) |
| $\beta_T$, $\beta_\eta$ | coefficients of linear calendar time $w_t$ and of the expectancy covariate $\eta_t$ in an analysis model | local |
| $u_i$ | participant random intercept of an analysis model | canonical |
| $\varpi_p$ | expected share of path $p$, $1/P$ ($\pi$ is reserved for power) | local |
| $\sigma_u^2$, $\rho_a$ | random-intercept variance and CAR(1) parameter (per week) of the working covariance | local (`docs/56`) |
| $\rho$, $c_1$ | AR(1) serial correlation per week (0.7); cross-component correlation (0.2) | canonical; local (`docs/40`) |
| $\rho_{cs}$ | compound-symmetry correlation (equal for every pair of visits) | local |
| $\sigma_{BR}$ | SD of $BR$ (8) | canonical |
| $\sigma_{\varepsilon,0}$, $r_{BL}$ | SD of measurement error added to the observed baseline; its reliability $\sigma_{BL}^2/(\sigma_{BL}^2 + \sigma_{\varepsilon,0}^2)$ | local |
| $\sigma_{\varepsilon,\text{OL}}$, $\sigma_{\varepsilon,\text{bl}}$ | SD of measurement error at open-label and blinded visits | local (`docs/46`) |
| $\kappa$ | mean model SE over empirical SD | canonical |
| $\alpha$, $\pi$ | nominal level (0.05); power | canonical |
| $n_{\text{reps}}$ | Monte Carlo replicates per cell | canonical |

Table: Notation index for symbols used in this protocol

*Revised 2026-10-09 from docs/58, Section 12.3, and the notation
decisions of that date:* the rows for $h$, $\mathcal T$, $r_1$, the three estimands,
$\hat\beta_D^{\text{exp}}(h)$, $\bar d(h)$, $b_0(h)$, $T(h)$,
$c_{\max}$, $\beta_T$, $\beta_\eta$, $u_i$, $\varpi_p$, $\sigma_u^2$,
$\rho_a$ and $\rho_{cs}$ were added, and the Gompertz row now names its parameters.
$\eta_t$ is the design expectancy, a covariate, never a coefficient.

Collisions avoided deliberately. `docs/37` uses $\tau^2$ for
the residual variance of the paired contrast given the biomarker; that
quantity is written $\sigma^2_{\Delta \mid B}$ here, and $\tau^2$ is
reserved for the participant-by-treatment variance (treatment-effect
heterogeneity), as in the literature on series of N-of-1 trials
[2, 3]. `docs/57` writes $\tau^2$ for the random-intercept variance and
$\theta$ for $\log t_{1/2}^{a}$; here the random-intercept variance is
$\sigma_u^2$ and $\theta$ is not used (*revised 2026-10-09 from
docs/58, Notation index*). The canonical $\delta$ (standardized bias)
and $\theta$ (generic interaction coefficient) are not used for the
drug effect. $\rho$ is reserved for the AR(1) correlation per week; a
compound-symmetry correlation is written $\rho_{cs}$. $r_1$ is
unrelated to the Gompertz rate $r$ and to the reliability $r_{BL}$.

**Sign convention** (NOTATION.md, rule 5). The components are
non-negative reductions in symptom severity, so an increase in any of
them lowers $Y$. The drug main effect $\beta_D$ is therefore negative
when the drug helps; $m_{BR}$, $f_m$ and $\tau$ are non-negative. Bias
is reported on the $\beta_D$ scale, so a positive bias is a bias
toward no effect.

**Units.** Carryover half-lives are in **weeks** throughout, matching
the refined process and the interaction manuscripts (1 week = 7 days).
*Revised 2026-10-09:* `NOTATION.md`, Part 1, "Units", previously
assigned days to the main-effect manuscript; by the user's decision of
2026-10-09 it now assigns weeks to the rebuilt paper 04, so this
protocol no longer departs from it. Gompertz times ($t_{od}$, $t_{pb}$,
$w_t$) are also in weeks. Effect sizes are given in CAPS points and in
SD units (Section 3.5).

**Evidence labels.** As in `docs/53`: *quoted*, *inspected* (read from
source), *derived* (follows algebraically), *verified* (computed by
running code; here only arithmetic checks in the session scratchpad,
`estimand-check.R`, not a project script), *inferred* and *assumed*.

**Glossary.**

- **The refined process.** The recommended data-generating process
  (DGP) of paper 01, Table `tab:spec`; "Refined" in tables. **The
  published process**: Hendrickson's code at commit 58b32a9
  (`NOTATION.md`, "Named data-generating processes").
- **Hybrid.** The Hendrickson design: open label to week 8, blinded
  discontinuation in weeks 8 to 12, a blinded two-period crossover in
  weeks 12 to 20; visits at weeks 4, 8, 9, 10, 11, 12, 16, 20 [4].
- **Multi-crossover.** The design of the current paper 04: eight
  one-week periods, four on drug and four on placebo in a random order
  per participant. It is not the Hybrid design.
- **Drug-expectancy confounding.** In open-label phases every visit is
  on drug at full expectancy, so the drug state and $\eta_t$ are
  correlated (0.588 over the nine Hybrid rows; `compendium-mathematics.tex`,
  Section "The published Hybrid construct", verified there). An analysis
  that does not model $\eta_t$ credits part of $PB$ to the drug.
- **Heterogeneous null.** A cell with $f_m = 0$ and $\tau > 0$: the mean
  drug effect is zero but individual effects are not.
- **Size-adjusted power, shortfall, minimax regret.** As in `docs/46`
  and `docs/47`: power at the 5% quantile of the analysis's own null
  $p$-values; an analysis's shortfall in a cell is its distance from the
  best surviving analysis there; the selected analysis minimizes its
  largest shortfall.
- **Oracle analysis.** An analysis that uses information the trial
  would not have (a true half-life, the true Gompertz shape, an
  error-free baseline). Reported for reference; never selectable.

## 1. Why the rebuild

The current paper's simulation does not study the design, the process
or the inference that the paper describes. Seven problems were found
(inspected, 2026-10-08):

1. **The design is not the Hybrid design.** The driver generates a
   randomized 4 + 4 multi-crossover of one-week periods
   (`production-run.R`, lines 80 to 85: `seq_base` of four 1s and four
   0s, permuted by `sample()`), while the abstract and Results call the
   N-of-1 arm "hybrid" (`report.Rmd`, line 8; lines 1216 to 1221).
   Methods describes the multi-crossover correctly (lines 515 to 517).
2. **The outcome has no placebo response, natural course, expectancy
   or onset.** The outcome is
   `baseline + Dbc * true_effect + carryover_adj + eps`
   (`production-run.R`, line 103). There is no
   $PB$, no $TV$, no $\eta_t$ and no Gompertz onset, so the
   drug-expectancy confounding of the open-label phase never arises,
   and the analysis's linear time term is fitted to a mean with no
   trend.
3. **No treatment-effect heterogeneity.** Every participant has the
   same effect (`true_effect`, line 103); the paper says that
   between-patient heterogeneity lies "purely" in baseline frequency
   (`report.Rmd`, lines 1362 to 1374), although its own literature
   review rests on $\tau^2$ (lines 153 to 252).
4. **Carryover exclusion uses oracle information.** The analysis drops
   the rows flagged `is_carryover`, a column the DGP writes
   (`production-run.R`, lines 110 and 120 to 127). A trial does not
   know which visits carry residual drug effect.
5. **Carryover is truncated and on a day scale.** Carryover reaches
   only the period immediately after a drug period (lines 94 to 101),
   and half-lives are in days (lines 53 and 95), off the scale of every
   other paper in the compendium.
6. **The inference labeled Kenward-Roger is not.** The test uses the
   `nlme` table's degrees of freedom (lines 172 to 178), which are the
   containment (between-within) degrees of freedom, not Kenward-Roger
   and not Satterthwaite, although the header (line 13) and Methods
   (`report.Rmd`, lines 564 to 568) say otherwise; the reduced model
   fitted for the correction (lines 159 to 170) is never used.
7. **The randomized-trial covariate is measured without error.** The
   RCT endpoint is generated from, and adjusted for, the participant's
   true intercept (`production-run.R`, lines 197 to 204 and 235). An
   ANCOVA on an error-free covariate removes all between-participant
   baseline variance, which favors the RCT in a way no trial can
   reproduce. The RCT arm also has no placebo response or natural
   course.

The compendium's own documents had already recorded parts of this
(`compendium-mathematics.tex`, Section "The treatment main-effect
track", notes of 2026-10-08), and the remediation plan
(`analysis/report/04-treatment-main-effect/00-remediation-plan.md`,
Problems 4 and 5) asked for a correct small-sample test and an
estimand defined under the Gompertz process. Neither can be repaired
inside the current DGP, because the DGP lacks the features that make
the estimand nontrivial. Hence a rebuild rather than a revision.

## 2. Aims (A)

**Primary aim.** Estimate the size and power of tests of the drug main
effect $\beta_D$ in the Hybrid design under the refined process, and
compare them with the OL+BDC design, the two-period crossover (CO) and
a parallel-group randomized controlled trial (RCT), all simulated from
the same component model at the same total $N$.

**Secondary aims.**

1. Estimate the participant-by-treatment variance $\tau^2$, and assess
   how well each design recovers individual drug effects.
2. Recommend an analysis for the drug main effect in each design, by
   the working decision rule of Section 6.4.
3. Quantify, in closed form first, how much of each design's power
   comes from randomized between-path contrasts and how much from
   modeled within-participant contrasts, and how much bias the
   drug-expectancy confounding introduces into each analysis.

**What the study does not aim at.** The biomarker interaction
(papers 01, 02 and the stress test of `docs/46`); optimal design
(`docs/39`); the comparison of response-curve families (paper 07).

**Prior work.** Appendix A reviews the published Monte Carlo and
closed-form studies of the treatment main effect in N-of-1, aggregated
N-of-1, multi-crossover and single-case designs, compares their
data-generating mechanisms and analyses with this protocol, and lists
the proposals that follow from that review (Section 8, items 12 to 20).

## 3. Data-generating mechanisms (D)

### 3.1 The refined process

The base DGP is the refined process of paper 01 (Table `tab:spec`;
Appendix B), with the biomarker uncoupled ($c_{bm} = 0$) in the base
cell. For participant $i$ on path $p$,

$$
Y_{i0} = \mathrm{BL}_i, \qquad
Y_{it} = \mathrm{BL}_i - \bigl(TV_{it} + PB_{it} + BR_{it}\bigr)
  - v_i\, \Psi_{it} + \varepsilon^{\text{me}}_{it}, \quad t \ge 1 .
$$

- **Components and means.** $E(TV_{it}) = G(w_t; m_{TV}, 5, 0.35)$ on
  calendar time; $E(PB_{it}) = \eta_t\, G(t_{pb}; m_{PB}, 5, 0.35)$;
  $E(BR_{it}) = G(t_{od}; m_{BR}, 5, 0.42)$ on drug, with the published
  ceilings $m_{TV} = m_{PB} = 6.506$ and $m_{BR} = 10.986$ unless a cell
  varies them. $\mathrm{BL}_i$ has mean 83.07 and SD 18.48 and is
  independent of the components (`docs/40`, Section 3).
- **Expectancy.** $\eta_t$ scales both the mean and the SD of $PB$
  (SD $10\eta_t$); this is the only route by which blinding enters the
  data (`docs/40`, Section 3).
- **Covariance.** Separable AR(1), $K \otimes A$ with
  $A_{ts} = \rho^{|w_t - w_s|}$, $\rho = 0.7$ per week, and
  cross-component correlation $c_1\rho^{|w_t - w_s|}$, $c_1 = 0.2$; SDs
  10 ($TV$), $10\eta_t$ ($PB$) and 8 ($BR$) at every visit, on and off
  drug (paper 01, "Why separability").
- **Anchored carryover.** Off drug after exposure,
  $E(BR_{it}) = \mu_{BR,\text{last}}\, 2^{-t_{sd}/t_{1/2}}$, decayed once
  from the last on-drug mean, not by the published compounding
  recursion (`docs/53`; paper 01, Appendix B.3). The same decay applies
  to the drug-response profile $\Psi_{it}$.
- **Half-life.** Pharmacological range $t_{1/2} \in \{0.1, 0.2\}$ week
  (17 to 34 hours, about 7 to 15 plasma half-lives of prazosin; paper
  01, "Pharmacological grounding"), base 0.2. A **behavioral carryover
  scenario** at $t_{1/2} = 1$ week is simulated and labeled as such
  everywhere it is reported: it is not a pharmacological claim about
  prazosin.
- **Re-exposure.** As in the published process, $t_{od}$ restarts at zero
  when the drug is resumed (justified by re-titration of prazosin).
  Combined with anchored carryover, this gives a downward step in the
  mean at the first re-exposure visit when the residual exceeds the
  restarted curve (paper 01, Appendix B.3).
- **Validity.** The generator stops on a correlation matrix that is not
  positive definite (paper 01, Appendix B.2). With $c_{bm} = 0$ no cell
  is near a ceiling; the stop is kept for the cell with $c_{bm} > 0$
  (Section 3.6).
- **Allocation.** Each participant is assigned to a path, or to an RCT
  arm, uniformly at random in every trial, redrawn if a path or arm
  receives fewer than 2 participants (the rule of
  `11-strawman-stress-test.R`, `RERUN-anchored.md`, Section 1). $N$ is
  the total (NOTATION.md, rule 4).

### 3.2 New element: a random drug-effect slope

The refined process has no participant-by-treatment variance: the $BR$
component varies between participants with SD 8 equally on and off
drug, so its person-level part cancels in an on-minus-off contrast
(`docs/51`, Feature 8; derived). Treatment-effect heterogeneity is
added as a participant-specific ceiling,

$$
v_i \sim N(0, \tau^2), \quad \text{independent of everything else},
\qquad \Psi_{it} = E(BR_{it})/m_{BR},
$$

that is, $\Psi_{it}$ is the drug-response mean of participant $i$'s path
computed with unit ceiling, including onset and anchored carryover; it
does not depend on $m_{BR}$ and is defined when $m_{BR} = 0$.
Participant $i$'s expected drug response is then
$(m_{BR} + v_i)\, \Psi_{it}$: heterogeneity follows the drug effect, is
zero before exposure and decays with the drug effect, as `docs/51`
recommends ("let the drug-response variance follow the drug effect").
$\tau$ is set in CAPS points, independent of $m_{BR}$, so that
heterogeneous nulls ($m_{BR} = 0$, $\tau > 0$) exist. A normal $v_i$
allows a participant to be harmed when $\tau > m_{BR}$; that is
intended in the heterogeneous null. Whether $v_i$ should instead scale
the binary on-drug state, or be correlated with $\mathrm{BL}_i$, is an
open decision (Section 8, item 9).

### 3.3 Designs

| Design | $P$ | Visit weeks (after baseline) | $\eta_t$ | Paths (drug state by visit) |
|---|---|---|---|---|
| Hybrid | 4 | 4, 8, 9, 10, 11, 12, 16, 20 | 1, 1, then 0.5 | 11110010, 11110001, 11100010, 11100001 |
| OL+BDC | 2 | 4, 8, 12, 16, 17, 18, 19, 20 | 1 to week 16, then 0.5 | 11111100, 11111000 |
| CO | 2 | every 2.5 weeks to week 20 | 0.5 | 11110000, 00001111 |
| RCT (parallel) | 2 arms | 4, 8, 9, 10, 11, 12, 16, 20 | 0.5 | 11111111 or 00000000 |
| Multi-crossover (optional) | random | eight weekly visits | 0.5 | random 4 + 4 per participant |

Table: Designs, randomization paths, schedules and expectancy

The first three rows are the designs of `analysis/scripts/quick-sim/
hendrickson-problems/04-power-simulation.R`, lines 74 to 83
(inspected). All have 8 post-baseline visits over 20 weeks [4].

**The RCT arm.** A parallel-group comparator of the kind that
motivated the Hybrid design [13], double-blind throughout ($\eta_t = 0.5$ from week 0),
so $PB$ and $TV$ follow the same curves as in the N-of-1 designs and are
identical in the two arms. The active arm has $t_{od} = w_t$; the
placebo arm has $BR$ mean 0. Two properties need stating:

- **Inherited $BR$ noise.** Under the refined process the placebo arm
  still carries $BR$ noise with SD 8, as an unexposed crossover path
  does. This is the convention of the published construct, not
  pharmacology (`docs/51`, Feature 8); it adds 64 to the variance of
  every RCT measurement. Whether to keep it is an open decision
  (Section 8, item 10).
- **Baseline measured with error.** The analysis sees
  $Y^{*}_{i0} = \mathrm{BL}_i + \varepsilon_{i0}$,
  $\varepsilon_{i0} \sim N(0, \sigma_{\varepsilon,0}^2)$, at reliability
  $r_{BL}$. Published test-retest reliabilities of the CAPS are high
  but below 1 [14], so an error-free baseline is the optimistic limit.
  The same observed baseline is used in the N-of-1 designs,
  where the baseline row then has its own noise, which also removes the
  degenerate noise-free baseline row noted in `docs/41`, Proposition 5.
  The default and the grid are an open decision (Section 8, item 7).

**The multi-crossover.** Under the refined process with re-titration,
a single on-drug week reaches 3.1% of the ceiling (0.34 points at the
published $m_{BR}$; verified, arithmetic), so the design of the current
paper cannot express the drug effect of a titrated drug except through
runs of consecutive on-drug weeks. It is offered only as an optional
bridge to the current paper (Section 8, item 5).

### 3.4 Factors and grids

| Factor | Levels | Base |
|---|---|---|
| Drug effect $f_m = m_{BR}/10.986$ | 0, 0.2, 0.4, 0.6, 0.8, 1.0 | curve; 0.4 in stress cells |
| Heterogeneity $\tau$ (points) | 0, 2.75, 5.49 ($\tau/10.986$ = 0, 0.25, 0.5) | 0 |
| Total $N$ | 35, 70, 150 | 70 |
| Carryover $t_{1/2}$ (weeks) | 0, 0.1, 0.2 (pharmacological); 1.0 (behavioral, labeled) | 0.2 |
| Decay form | exponential; Weibull shape 0.5 and 2 at $t_{1/2} = 1$ | exponential |
| Natural course $m_{TV}$ | $\times 0$, $\times 1$, $\times 2$ | $\times 1$ (6.506) |
| Placebo response $m_{PB}$ | $\times 0$, $\times 1$, $\times 2$ | $\times 1$ (6.506) |
| Dropout | none; 20% monotone MCAR; 20% monotone MAR on $Y_{i,t-1}$ | none |
| Phase measurement error | none; $\sigma_{\varepsilon,\text{OL}} = 3$, $\sigma_{\varepsilon,\text{bl}} = 6$ | none |
| Baseline reliability $r_{BL}$ | 1; 0.8 ($\sigma_{\varepsilon,0} = 9.24$) | open (Section 8) |

Table: Data-generating factors, levels and base values

The dropout and measurement-error mechanisms are those of the options
patch and the stress test (`analysis/patches/orig-to-refined/
README.md`, options table; `docs/46`, Section 5).

### 3.5 Effect size in outcome and SD units

The drug-effect axis is defined on the DGP, by the ceiling $m_{BR}$,
because each design's estimand takes a different value for the same
drug (Section 4). For an SD scale we propose the SD of a blinded
post-baseline score conditional on baseline,

$$
\sigma_{Y \mid \mathrm{BL}}^2 = \sigma_{TV}^2 + (10\eta_t)^2 + \sigma_{BR}^2
  + 2c_1\bigl(10\sigma_{TV}\eta_t + \sigma_{TV}\sigma_{BR}
  + 10\eta_t\sigma_{BR}\bigr) = 257
$$

at $\eta_t = 0.5$, so $\sigma_{Y \mid \mathrm{BL}} = 16.03$ (derived;
verified, arithmetic; the same value as `docs/41`, Proposition 3,
because marginal variances do not depend on the serial structure). The
unconditional SD is $\sqrt{598.5} = 24.46$.

| $f_m$ | 0 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 |
|---|---|---|---|---|---|---|
| $m_{BR}$ (points) | 0 | 2.20 | 4.39 | 6.59 | 8.79 | 10.99 |
| $m_{BR}/\sigma_{Y \mid \mathrm{BL}}$ | 0 | 0.14 | 0.27 | 0.41 | 0.55 | 0.69 |

Table: Drug-effect grid in multiplier, outcome and SD units

The grid was chosen to keep the N-of-1 designs off the power ceiling.
That rests on rough derived values (an RCT SE near 3.8 points at
$N = 70$ with an error-free baseline, Section 5.1; a Hybrid SE near 1
point, taken from the compound-symmetry results of `docs/40`,
Section 9) and must be checked against the closed form of Section 5.1
before the run.

### 3.6 Draft cell table

Every block varies one factor from the base cell (Hybrid, OL+BDC, CO
and RCT; refined process; $N = 70$; $t_{1/2} = 0.2$; $\tau = 0$; no
dropout or measurement error). Stress blocks have a null twin
($f_m = 0$) and an alternative at $f_m = 0.4$. The RCT is not simulated
in blocks whose factor does not reach it.

| Block | Factor varied | Levels | $f_m$ | Designs | Cells (null / alt) |
|---|---|---|---|---|---|
| 1 | power curves | base | 0 to 1.0 | 4 | 24 (4 / 20) |
| 2 | heterogeneity | $\tau$ = 2.75, 5.49 | 0, 0.4 | 4 | 16 (8 / 8) |
| 3 | sample size | $N$ = 35, 150 | 0, 0.4 | 4 | 16 (8 / 8) |
| 4 | carryover | $t_{1/2}$ = 0; 1.0; Weibull 0.5 and 2 at 1.0 | 0, 0.4 | 3 (no RCT) | 24 (12 / 12) |
| 5 | natural course | $m_{TV} \times 0$, $\times 2$ | 0, 0.4 | 4 | 16 (8 / 8) |
| 6 | placebo response | $m_{PB} \times 0$, $\times 2$ | 0, 0.4 | 4 | 16 (8 / 8) |
| 7 | dropout | MCAR, MAR (20%) | 0, 0.4 | 4 | 16 (8 / 8) |
| 8 | phase measurement error | 3 / 6 | 0, 0.4 | 4 | 8 (4 / 4) |
| 9 | baseline error | $r_{BL} = 0.8$ | 0, 0.4 | 4 | 8 (4 / 4) |
| Total | | | | | 144 (64 / 80) |

Table: Draft cell table: blocks, levels and cell counts

Optional cells, not counted above: $t_{1/2} = 0.1$ (expected to be
indistinguishable from 0.2 at weekly visits); explained heterogeneity
through covariance moderation, $c_{bm} = 0.30$, $\tau = 0$ (the
biomarker then induces a drug-effect SD of $c_{bm}\sigma_{BR} = 2.4$
points on drug; derived from `docs/40`, Section 5); and the
multi-crossover in blocks 1 and 4.

**Misspecified analysis half-life** needs no new data. *Revised
2026-10-09 from docs/58, Section 12.3:* in every cell, the
$D_{bc}$-based analyses are fitted at every $h \in \mathcal T =
\{0, 0.2, 0.5, 1, 2\}$ weeks (Section 5.2, note "Misspecified
half-life"), replacing the earlier plan of $t_{1/2}^{a} \in \{0, 0.5,
1\}$ in the cells of block 1 and the $t_{1/2} = 1$ cells of block 4.

**Replicates and MCSE targets.** Two options (Section 8, item 6):

| Option | Null cells | Alternative cells | MCSE of size at 0.05 | Largest MCSE of power | MCSE of coverage at 0.95 |
|---|---|---|---|---|---|
| R1 | 2,000 | 1,000 | 0.0049 | 0.016 | 0.0049 (null), 0.0069 (alt) |
| R2 (as `docs/46`) | 1,000 | 500 | 0.0069 | 0.022 | 0.0069 (null), 0.0097 (alt) |

Table: Replicate options with Monte Carlo standard error targets

R1 gives 128,000 null and 80,000 alternative trials (208,000), of which
about 165,000 are N-of-1 trials; R2 halves both. The MCSE of bias is
$\mathrm{SD}(\hat\beta_D)/\sqrt{n_{\text{reps}}}$, about 0.03 points
for an empirical SD of 1 at $n_{\text{reps}} = 1{,}000$.

## 4. Estimands (E)

### 4.1 Candidate definitions of $\beta_D$

Under the refined process the drug effect is not one number: it rises
on the Gompertz curve, restarts at re-exposure and decays after
discontinuation, and the open-label visits add full expectancy. Every
analysis with a single drug coefficient targets some weighted average
of a drug response that ranges from 4.3 to 10.2 points across the
Hybrid visits (`docs/40`, Section 9). Three definitions were
considered. Values are at the published ceiling, $t_{1/2} = 0$, equal
path weights (verified, arithmetic, `estimand-check.R`):

1. **Full-effect contrast**, $\beta_D^{\text{full}} = -m_{BR}$: the
   reduction after long exposure, 10.99 points. Design-free, but no
   design observes it: no Hybrid or crossover visit is past 10 weeks of
   continuous exposure, and the re-exposure visits reach 39% of the
   ceiling. Only an analysis that knows the Gompertz shape can target
   it (the oracle analysis OP, Section 5); every other analysis is
   biased toward zero by construction, which makes bias uninformative
   about the analyses.
2. **Contrast at the design's time on drug**,
   $\beta_D^{\text{exp}} = -\bigl[\overline{E(BR)}_{\text{on}} -
   \overline{E(BR)}_{\text{off}}\bigr]$, the on-minus-off difference in
   the mean drug response over the design's own visits, averaged over
   participants, net of $PB$ and $TV$ by definition. It is what the
   drug did at the times the design exposed people. Values: Hybrid
   $-7.22$ (paths 1 and 2, $-7.55$; paths 3 and 4, $-6.89$); OL+BDC
   $-9.35$; CO $-6.70$; RCT at week 20, $-10.97$. Under carryover the
   off-drug mean is positive and the value shrinks; under expectancy
   the definition excludes $PB$, so an analysis that credits
   open-label expectancy to the drug is biased away from zero.
3. **Blinded on-minus-off contrast at matched expectancy**,
   $\beta_D^{\text{bl}}$: the same contrast restricted to blinded
   visits. Hybrid $-7.56$ (paths 1 and 2, $-8.09$; 3 and 4, $-7.04$);
   undefined for OL+BDC except at its two or three blinded visits. Its
   randomized version, the equal-weight average of the between-path
   contrasts at weeks 10, 16 and 20, is $-6.25$ in the Hybrid design,
   the probability limit of the visit-factor analysis (`docs/40`,
   Section 9; `compendium-mathematics.tex`, Equation `eq:bp-fwl`).
   Matching expectancy removes the $PB$ step but not the natural
   course: the blinded within-path contrast without time adjustment has
   a null bias of $+0.37$ points, against $+1.29$ for the contrast over
   all visits (verified, arithmetic).

**Implications under Gompertz onset.** Definitions 2 and 3 depend on
the schedule, so designs differ in their estimand for the same drug:
the RCT's week-20 contrast is 1.5 times the Hybrid's $\beta_D^{\text{exp}}$.
That difference is a real property of the designs (short re-exposures
of a titrated drug), and it will appear in power. It means that power
can be compared across designs at a common $m_{BR}$, but bias and
coverage only against each design's own estimand value.

**Implications under expectancy.** For definition 2, the limits of
`docs/40`, Section 9, at $t_{1/2} = 0$ are linear time $-8.21$, linear
time with $\eta_t$ $-7.28$, spline $-6.22$, visit factor $-6.25$, all
under working independence and the published allocation. Against
$\beta_D^{\text{exp}} = -7.22$ linear time overstates the effect by about
1 point and the expectancy model is close (inferred from those limits
and the value above; the weightings differ slightly, so this must be
recomputed in the closed-form step).

**The test is common to all definitions.** Under $m_{BR} = 0$ every
definition is zero, so size and power refer to the same null
hypothesis in every design. The estimand choice affects bias, coverage
and the interpretation of $\hat\beta_D$, not the rejection rates.

### 4.2 Recommendation (open decision)

We recommend **definition 2**, $\beta_D^{\text{exp}}$, as the primary
estimand for estimation, with definition 3's randomized version
reported as a secondary estimand for the Hybrid design. Reasons: it is
defined for every design including the RCT; it is a pharmacological
quantity, so expectancy confounding shows as bias rather than being
built into the target; and it is attainable by analyses a trial can
run. Its cost is that it is design-specific. This choice is open
(Section 8, item 1). Its value in every cell is computed from the
population means of the cell, in closed form, and stored with the cell.

**The scale of a $D_{bc}$ coefficient.** *(Revised 2026-10-09 from
docs/58, Section 12.3; derivation in docs/58, Section 4.2.)* The
coefficient of $D_{bc,it}(h)$ is an effect per unit of the indicator
built at the analysis half-life $h$, so its scale changes with $h$: a
longer $h$ codes the off-drug visits as partly exposed, and the
coefficient grows in magnitude even when the data do not change. It is
therefore not comparable across $h$, nor with $\beta_D^{\text{exp}}$.
Every analysis with a $D_{bc}(h)$ term also reports the model-implied
exposure contrast
$$
\hat\beta_D^{\text{exp}}(h) = \hat\beta_D(h)\,\bar d(h), \qquad
\bar d(h) = 1 - \frac1P\sum_p \bar x_p^{\,\text{off}}(h),
$$
where $\bar x_p^{\,\text{off}}(h)$ is the mean of $D_{bc}(h)$ over the
post-baseline off-drug visits of path $p$; this is the definition of
$\beta_D^{\text{exp}}$ applied to the fitted drug part of the mean, since
$D_{bc} = 1$ on drug. Bias and coverage against $\beta_D^{\text{exp}}$
(Section 6.1) are computed on this scale, and its interval is the CR2
interval of $\hat\beta_D(h)$ multiplied by $\bar d(h)$. In the Hybrid
design $\bar d(h)$ is 1, 0.991, 0.906, 0.755 and 0.551 at $h$ = 0, 0.2,
0.5, 1 and 2 weeks (verified in docs/58), so at the default
$t_{1/2}^{a} = 0.2$ the two scales differ by 1%. In the closed form of
docs/58 (Hybrid, $f_m = 0.4$, $N = 70$), F1's exposure contrast is
within 0.09 point (3%) of $\beta_D^{\text{exp}}$ when $h = t_{1/2} =
0.2$, including its null bias, and a wrong long $h$ attenuates it by
14% to 30% (docs/58, Table 4).

### 4.3 Secondary estimands

- **$\tau^2$, on the scale of the primary estimand.** Participant $i$'s
  exposure-averaged effect is $-(m_{BR} + v_i)\,\bar\Psi_d$, with
  $\bar\Psi_d$ the design's on-minus-off mean of $\Psi_{it}$; its variance is
  $\tau_d^2 = \tau^2 \bar\Psi_d^2$ (derived). The random-slope analyses of
  Section 5 target $\tau_d^2$, not $\tau^2$; both are reported.
- **Individual effects.** The participant-level value
  $-(m_{BR} + v_i)\,\bar\Psi_{d,p(i)}$, conditional on $v_i$ and excluding
  the transient $BR$ noise.
- **Size.** The rejection rate under $m_{BR} = 0$, with and without
  heterogeneity.

## 5. Methods (M)

### 5.1 Closed form first

Before any simulation, every analysis whose mean model is fixed is
evaluated in closed form, extending `docs/37` and `docs/39` from the
biomarker slope to the drug coefficient.

- **Probability limit.** For an analysis with mean design $X_p$ on path
  $p$ and working covariance $V$,
  $\operatorname{plim}\hat\beta_D = \bigl[(\sum_p \varpi_p X_p'V^{-1}X_p)^{-1}
  \sum_p \varpi_p X_p'V^{-1}\mu_p\bigr]_D$, with $\varpi_p = 1/P$ the
  expected path share (local symbol; $\pi$ is reserved for power) and $\mu_p$ the true mean vector (`compendium-mathematics.tex`,
  "Time adjustment and the probability limit of the drug coefficient",
  extended from working independence to a general $V$; derived). Its
  difference from the estimand value is the bias, and under
  $m_{BR} = 0$ a nonzero limit, the null bias $b_0(h)$, predicts size
  inflation. *(Revised 2026-10-09 from docs/58, Section 12.3; values
  from docs/58, Sections 5.1 and 5.2, Table 3.)* The earlier sentence
  here predicted a nonzero null limit for linear time without $\eta_t$
  only (inferred from `docs/40`, Section 9). The closed form of docs/58
  shows a nonzero limit for F1 as well; Table "Closed-form null bias"
  below gives the values. $b_0(h)$ is exactly 0 for every $h$ and every
  working covariance when the time adjustment contains visit effects
  (VF), because the non-drug mean is then in the span of the non-drug
  columns, which are common to all paths (docs/58, Proposition 3;
  derived); linear time, with or without $\eta_t$, does not satisfy that
  condition.
- **Working covariance.** *(Revised 2026-10-09 from docs/58, Section
  12.3.)* The likelihood analyses are evaluated at the limit of their
  own working covariance, not at working independence: the
  random-intercept CAR(1) member (parameters $\sigma_u^2$, $\rho_a$ and
  the residual variance) closest to the truth in Kullback-Leibler
  distance, found by minimizing the expected ML criterion per
  participant, $\sum_p \varpi_p[\log\lvert V\rvert +
  \operatorname{tr}\{V^{-1}(\Sigma_p + \mathbf m_p\mathbf m_p')\}]$,
  with $\mathbf m_p = \boldsymbol\mu_p - X_p\boldsymbol\xi^*$ the mean
  misfit at the limit $\boldsymbol\xi^*$ of the fixed effects,
  recomputed for each candidate $V$ (REML has the same limit; derived
  in docs/58, Proposition 2). In the Hybrid design under the null this
  gives $\rho_a = 0.646$ per week for F1. The choice matters: F1's
  null bias at $h = 0.2$ is $+0.29$ points under working independence,
  $+0.18$ under the true covariance and $+0.33$ under the
  random-intercept CAR(1) limit, with standard errors 1.63, 1.21 and
  1.22 at $N = 70$ (verified in docs/58).
- **Variance.** The sandwich
  $\mathrm{Var}(\hat\beta_D) = N^{-1}[B^{-1}\{\sum_p \varpi_p X_p'V^{-1}
  (\Sigma_p + \mathbf m_p\mathbf m_p') V^{-1} X_p\}B^{-1}]_{DD}$,
  $B = \sum_p \varpi_p X_p'V^{-1}X_p$,
  with $\Sigma_p$ the true covariance of $\mathbf Y_i$: the refined
  separable covariance plus $\mathrm{BL}$ variance, plus
  $\tau^2 \boldsymbol\Psi_p\boldsymbol\Psi_p'$, plus measurement error.
  *(Revised 2026-10-09 from docs/58, Section 12.3.)* The mean-misfit
  term $\mathbf m_p\mathbf m_p'$ is included because path allocation is
  random: a misspecified mean contributes between-path variance to the
  estimate. Under fixed path counts it drops, and CR2 estimates the form
  that includes it (docs/58, Proposition 2). For the
  summary-measure analyses PW and PR this is exact (the contrast is a
  fixed linear combination, as E9 is in `docs/37`, Section 4); for the
  likelihood analyses it holds at the limit of the working covariance
  above (approximate).
- **Contrast variance by design.** `docs/39`, Section 3, gives the
  contrast variances of the on-minus-off weights by design under CS and
  AR(1) (Hybrid 107.3 to 126.8 under AR(1) at $\rho = 0.8$ without
  measurement error, against 212.8 for CO; to be recomputed at the
  refined $\rho = 0.7$). For the main effect the same quadratic form
  $a'Ca$ applies, with two changes: the drug signal is $a'\mu_{BR}$
  instead of the coupling term, and the weights must also remove the
  non-drug mean, either by orthogonality to the time and expectancy
  terms or by the between-path construction. The heterogeneity adds
  $\tau^2(a'\Psi)^2$.
- **RCT.** For ANCOVA at week 20 with $r_{BL}$,
  $\mathrm{Var}(\hat\beta_D) \approx (2/N)\bigl(2\sigma^2_{20 \mid Y^*_0}
  + \tau^2 \Psi_{20}^2\bigr)$ at equal arms, with
  $\sigma^2_{20 \mid Y^*_0} = 598.5 - 18.48^2\, r_{BL}$ the placebo-arm
  residual variance (derived). At $N = 70$ and $\tau = 0$:
  SE 3.83 at $r_{BL} = 1$, 4.08 at 0.9, 4.31 at 0.8; normal-approximation
  power at the published ceiling 0.82, 0.77, 0.72 (verified,
  arithmetic; ignores the $t$ reference distribution).
- **Analytic power curves.** Power by noncentral $t$ at the analysis's
  reference degrees of freedom, over the $f_m$ grid, for every design.
  The simulation then validates the closed form in the cells of block 1
  and of blocks 2 and 4: agreement within 2 MCSE for rejection rates
  and within 2 MCSE for the empirical SD, as `docs/37`, Section 6 and
  `docs/39`, Section 9 required. Disagreement is reported, not tuned
  away.

| Design | Analysis | $b_0(0)$ | $b_0(0.2)$ | $b_0(1)$ | $b_0(2)$ | SE at $h = 0.2$ | Size at $h = 0.2$ |
|-----------|---------|--------|--------|--------|--------|----------|----------|
| Hybrid | S | $-0.58$ | $-0.61$ | $-1.12$ | $-1.63$ | 1.14 | 0.083 |
| Hybrid | F1 | $+0.32$ | $+0.33$ | $+0.25$ | $+0.13$ | 1.22 | 0.058 |
| Hybrid | VF | 0 | 0 | 0 | 0 | 1.66 | 0.050 |
| OL+BDC | S | $-2.16$ | $-2.23$ | $-3.29$ | $-3.92$ | 1.18 | 0.473 |
| OL+BDC | F1 | $-0.22$ | $-0.22$ | $-0.09$ | $+0.60$ | 1.42 | 0.053 |
| OL+BDC | VF | 0 | 0 | 0 | 0 | 2.35 | 0.050 |
| CO | S | $-0.26$ | $-0.26$ | $-0.33$ | $-0.43$ | 1.48 | 0.053 |
| CO | F1 | 0 | 0 | $-0.06$ | $-0.15$ | 1.59 | 0.050 |
| CO | VF | 0 | 0 | 0 | 0 | 1.59 | 0.050 |

Table: Closed-form null bias $b_0(h)$ (points), standard error and normal-approximation size at $N = 70$, $t_{1/2} = 0.2$, random-intercept CAR(1) working covariance at its limit (quoted from docs/58, Table 3, verified there)

Because a fixed null bias does not shrink while the standard error
falls as $N^{-1/2}$, the size of F1 in the Hybrid design grows with
$N$: 0.054 at $N = 35$, 0.058 at 70, 0.068 at 150 and 0.086 at 300
(derived in docs/58 by scaling the standard error). S is not usable for
the main effect in OL+BDC and is liberal in the Hybrid design. This
motivates the size screen at $N = 150$ in Section 6.4. The null bias of
F1 is the sum of four row-group contributions that nearly cancel
($-0.52$, $+0.61$, $-1.21$, $+1.44$ points; docs/58, Section 6), so its
value is fragile to the shapes of the natural course and the placebo
curve, which the refined process fixes (argued in docs/58).

The closed form is also the gate for the grids: if it shows a design at
power above 0.99 or below 0.10 across most of the grid, the grid is
revised before the run, and the revision is recorded.

### 5.2 Analyses

The analyses adapt the fourteen stress-test analyses
(`compendium-mathematics.tex`, Table `tab:invd-analyses`;
`docs/47`) to the main effect by dropping the biomarker terms, and add
analyses that only the main effect needs. IDs shared with the stress
test denote the main-effect analog of the same model. Notation: $X_{it}$
is $D_{bc,it}$ at $t_{1/2}^{a}$, with the pre-exposure guard (0 at
baseline and before first exposure; paper 01, Appendix B.4), and
$X_{it} = D_{it}$ when $t_{1/2}^{a} = 0$. Default $t_{1/2}^{a} = 0.2$
week. All N-of-1 analyses use the long data with the observed baseline
as a row at $t = 0$ unless stated. RC is `nlme::lme` with a random
intercept and `corCAR1(form = ~ t | ptID)` on weeks, REML; RCH adds
`varIdent(form = ~ 1 | phase)` (baseline, open label, blinded); CR2 is
`clubSandwich::coef_test(vcov = 'CR2', test = 'Satterthwaite')`
clustered on participant [5, 6]; the `nlme` models follow Pinheiro and
Bates [95]. *(Revised 2026-10-09, notation decisions of that date.)* In
the mean-model column, linear time is the calendar week $w_t$, entering
as $\beta_T w_t$; "visit factor" is one indicator per visit; $\eta_t$
is the design expectancy. A model written $w_t + X + \eta_t$ is
$Y_{it} = \beta_0 + \beta_T w_t + \beta_D X_{it} + \beta_\eta\eta_t +
u_i + \varepsilon_{it}$, with $\beta_\eta$ the coefficient of the
expectancy covariate (local). Where a row is evaluated over the
analysis grid, $X_{it}(h)$ denotes the indicator at $h \in \mathcal T$.

| ID | Mean model | Covariance | Inference | df | Purpose |
|---|---|---|---|---|---|
| S | $w_t + X$ | RC | CR2 | Satterthwaite (CR2) | stress-test reference; exposes the expectancy confounding |
| F1 | $w_t + X + \eta_t$ | RC | CR2 | Satterthwaite (CR2) | expectancy-term model; proposed reference (strawman) |
| G1 | as F1 | RCH | CR2 | Satterthwaite (CR2) | phase variances |
| G2 | as F1 | RCH | model-based | containment | phase variances, model-based |
| NS | natural spline in $w_t$ (3 df, `ns(t, 3)`) $+ X + \eta_t$ | RC | CR2 | Satterthwaite (CR2) | flexible time (`docs/42`) |
| VF | visit factor $+ X$ | RC | CR2 | Satterthwaite (CR2) | randomized between-path contrasts only |
| A | $w_t + X$ | unstructured, `mmrm` | Kenward-Roger | Kenward-Roger | stress-test analog; known liberal at small $N$ |
| VA | visit factor $+ X$ | unstructured, `mmrm` | Kenward-Roger | Kenward-Roger | the MMRM convention [7] |
| R | $w_t + D$ | compound symmetry, `mmrm` | model-based | Satterthwaite (`mmrm`) | published-type reference [4] |
| XS | $w_t + D + \eta_t$, first off-drug visit after each discontinuation dropped by schedule | RC | CR2 | Satterthwaite (CR2) | period exclusion without oracle information; versus $D_{bc}$ |
| XO | as XS, dropping visits whose true carryover fraction exceeds 0.05 | RC | CR2 | Satterthwaite (CR2) | oracle; bridge to the v1 analysis |
| PW | per participant, blinded on-minus-off mean; one-sample $t$ | none | exact $t$ | $N - 1$ | E9 analog; time-confounded by construction |
| PR | per participant, $s_i = \sum_t (D_{p(i)t} - \bar D_t)Y_{it}$; path-stratified mean | none | $t$ | $N - P$ | randomization-based summary measure (the FWL estimator) |
| RS | as F1, random slope on $X$ correlated with the intercept | RC + slope | CR2 | Satterthwaite (CR2) | $\tau^2$, individual effects |
| RSL | as RS without $\eta_t$ | RC + slope | CR2 | Satterthwaite (CR2) | $\tau^2$ under linear time |
| OP | $w_t + \eta_t + Q$, $Q_{it} = \Psi_{it}$ the true profile | RC | CR2 | Satterthwaite (CR2) | oracle shape; targets $\beta_D^{\text{full}}$ |
| AN | RCT: week-20 score on arm and $Y^*_{i0}$ | OLS | HC2 | $N - 3$ | RCT primary |
| AN0 | as AN on the error-free $\mathrm{BL}_i$ | OLS | HC2 | $N - 3$ | oracle; reproduces the v1 RCT |
| AM | RCT: all visits, visit-by-arm, $Y^*_{i0}$ covariate, contrast at week 20 | unstructured, `mmrm` | Kenward-Roger | Kenward-Roger | RCT using every visit |
| MX-F1 | F1 at every $h \in \mathcal T$; statistic $\max_h \lvert T(h)\rvert$ | RC | CR2; critical value $c_{\max}$ from the influence-based multivariate normal | multivariate normal over $\mathcal T$ | data-chosen half-life with the selection calibrated; keeps F1's null bias |
| MX-VF | as MX-F1 with VF | RC | as MX-F1 | as MX-F1 | as MX-F1, without a null bias |
| RV | VF at every $h \in \mathcal T$, null-model weights; statistic $\max_h \lvert T(h)\rvert$, CR0-studentized | RC of the null model $w_t + \eta_t$, fitted once | randomization: 999 fixed-count permutations of the path labels, maximum inside each | none (exact under the sharp null) | size-exact companion with a data-chosen half-life; randomized between-path contrast only |
| RF | as RV with the F1 statistic | as RV | as RV | as RV | randomization test of F1; direction set by its label-invariant within-path part |
| F1A | F1 at the AIC-chosen $h \in \mathcal T$ (ML for the choice, REML refit) | RC | CR2, single-test reference | Satterthwaite (CR2) | labeled negative control; predicted size 0.088 (Hybrid, $N = 70$); never selected |

Table: Analyses of the drug main effect: mean model, covariance, inference, degrees of freedom and purpose

Notes on the table:

- **Kenward-Roger** is the `mmrm` implementation. The stress test used
  its linear variant (`Kenward-Roger-Linear`; compendium,
  "Inference"); the same variant is used here for comparability, which
  should be confirmed (Section 8, item 4) [8].
- **CR2 as default.** CR2 with Satterthwaite degrees of freedom is the
  default inference for every random-intercept model, because model-
  based and Kenward-Roger inference were anticonservative with flexible
  time in the Hybrid carryover cell and CR2 was not (`docs/42`, update
  of 2026-10-07; `docs/45`, Section 4.4), and A failed the size screen
  of the stress test (0.087 at $N = 35$; `docs/47`, Section 3).
- **HC2** for the RCT is the singleton-cluster case of CR2; it protects
  against the unequal arm variances that $\tau > 0$ creates.
- **F1 in CO.** With $\eta_t$ constant after baseline, $\eta_t$ is a
  baseline indicator there, so F1 frees the baseline mean.
- **XS against $D_{bc}$.** XS replaces the v1 oracle exclusion by a rule
  a trial can apply. In the Hybrid design it drops weeks 11 (paths 1
  and 2) or 10 (paths 3 and 4) and the first post-crossover off-drug
  visit.
- **PW and PR** are the main-effect counterparts of E9. PW keeps the
  time confounding of a within-path contrast (null bias $+0.37$ points
  over blinded visits; Section 4.1) and is included because its moments
  are exact and it shows the size of the confounding. PR uses only the
  visits at which paths differ in drug state, so its mean is the
  randomized contrast and its variance is exact.
- **Random slopes.** In the Hybrid design each participant has four or
  five on-drug and three or four off-drug visits, so the slope variance
  may be weakly identified and hit the boundary. $\tau_d^2$ is tested
  by the likelihood ratio against RC with the 50:50 mixture of
  $\chi^2_1$ and $\chi^2_2$ [9]; its interval is the `nlme::intervals()`
  approximation (Wald on the log SD scale), whose coverage is itself a
  result. Boundary and nonconvergence rates are reported.
- **Misspecified half-life.** *(Revised 2026-10-09 from docs/58,
  Section 12.3.)* S, F1, VF, RS and XS are fitted at every
  $h \in \mathcal T = \{0, 0.2, 0.5, 1, 2\}$ weeks in every cell, at a
  cost of five fits per analysis; VF is added to the list because MX-VF
  and RV need it. This replaces the earlier plan of refits at
  $t_{1/2}^{a} \in \{0, 0.5, 1\}$ in the cells of block 1 and the
  $t_{1/2} = 1$ cells of block 4 (Section 3.6). $\mathcal T$ is the
  grid of `docs/57` with 0.25 replaced by the default 0.2. Whether RS
  and XS need the whole grid, given their cost, is part of item 22.
- **Analyses added from docs/58** *(revised 2026-10-09 from docs/58,
  Section 12.3; method in docs/58, Sections 9.1 to 9.3).* The half-life
  enters the main-effect model only through the tested column, so it is
  not identified under the null, and any likelihood-based choice of $h$
  looks at the tested term [73, 74, 76]. Under F1, $h$ chosen by AIC
  is the $h$ with the largest likelihood-ratio statistic for the drug
  term, because the model without that term does not depend on $h$;
  referring its $t$ to a single-test critical value makes F1A a max-$|t|$
  test with the wrong critical value (docs/58, Proposition 5; derived).
  - *MX-F1, MX-VF.* The joint null distribution of $T(h)$ over
    $\mathcal T$ is estimated from the participants' influence
    contributions (the multiplier construction of Hansen [76]), and
    $c_{\max}$ is its 95% quantile of $\max_h \lvert T(h)\rvert$
    (`mvtnorm::qmvnorm`; a multivariate $t$ at the smallest CR2 degrees
    of freedom is the small-sample variant). Closed form (Hybrid,
    $N = 70$): $c_{\max} \approx 2.16$ (F1) and 2.14 (VF), against the
    Bonferroni 2.58; size 0.056 and 0.051. MX-F1 calibrates the
    selection but keeps F1's null bias; MX-VF has none. Davies' bound
    [74] gives the same constant for a continuous $h \in [0, 2]$ weeks
    to 0.003 (docs/58, Section 9.2) and is an optional check.
  - *RV, RF.* Exact under the sharp null for any statistic, including
    the selection, when the reference set follows the design's own
    randomization [77, 92, 93]. The null model $w_t + \eta_t$ has the
    same design on every path, so its fitted covariance does not depend
    on the labels; with weights and CR0 studentization taken from it,
    each relabeling changes only the path sums, and 999 relabelings of
    ten statistics take about 2 seconds per trial (docs/58, Section 9.3;
    verified there on one machine). Permutation is of the observed
    labels at the realized path counts, which is exact conditionally
    and hence unconditionally. Such a test uses only the randomized
    between-path contrasts, whatever statistic is permuted (docs/58,
    Proposition 6): its closed-form power at the base alternative is
    about 0.39 for both, against 0.63 for model-based F1, whose extra
    power rests on within-path contrasts protected only by the time and
    expectancy model. In the docs/58 Monte Carlo check RF rejected more
    often than that (0.507), because F1's label-invariant part sets the
    direction; this is legitimate but makes RF behave as a one-sided
    test of the randomized contrast. Under the weak null (block 2) an
    unstudentized permutation test is not exact; the studentized
    version is asymptotically valid [80 to 83], and the studentized F1 statistic is not pivotal unless
    $b_0 = 0$, which is why RV, not RF, is the proposed companion
    (item 23).
  - *F1A.* Included only as a labeled negative control: closed-form
    size 0.088 in the Hybrid design at $N = 70$, Monte Carlo 0.093 (600
    trials, MCSE 0.012; docs/58, Section 11). It is never selected and
    never reported as a test of $\beta_D$. The earlier exclusion rested
    on the interaction stress test (F1A, 0.057 pooled; `docs/47`,
    Section 8), where a selection blind to the tested term exists; for
    the main effect it does not.
  - *Not added.* AIC-weight averaging on the $\hat\beta_D^{\text{exp}}$
    scale [88, 89] is reported descriptively only (Section 6.1); the
    parametric bootstrap from the null model is centered at zero for F1
    while F1's real null limit is $b_0$, so it does not correct the
    null bias (docs/58, Section 9.5), and its rejection rate would need
    the warp-speed method [90] to be affordable; joint nonlinear estimation of $h$ has
    a supremum statistic at $\beta_D = 0$ [91, 94]; sample splitting
    and Bayesian marginalization over $h$ are not evaluated. Background
    on data-chosen carryover models in crossover trials [62 to 69],
    pharmacodynamic and lag models [70 to 72], weighted-average-power
    tests [75], single-case randomization tests [78, 79], randomization
    inference for effect variation [84, 85] and post-selection
    inference [86, 87] is reviewed in docs/58, Section 13.

## 6. Performance measures (P)

### 6.1 Measures

For each analysis and cell, with MCSEs as in Morris, White and
Crowther [1], Table 6:

| Measure | Definition | MCSE |
|---|---|---|
| Size (null cells) | $\hat\pi_0 = n_{\text{reps}}^{-1}\sum I(p_r < 0.05)$ | $\sqrt{\hat\pi_0(1 - \hat\pi_0)/n_{\text{reps}}}$ |
| Power | as above, alternative cells | as above |
| Size-adjusted power | rejection at the 5% quantile of the analysis's null $p$-values in the matching null cell | binomial, conditional on the quantile |
| Bias | $\bar\beta - \beta_D$, against the cell's estimand value | $\mathrm{SD}(\hat\beta)/\sqrt{n_{\text{reps}}}$ |
| Empirical SE | $\mathrm{SD}(\hat\beta)$ | $\mathrm{SD}(\hat\beta)/\sqrt{2(n_{\text{reps}} - 1)}$ |
| Model SE calibration | $\kappa = \overline{\mathrm{SE}}/\mathrm{SD}(\hat\beta)$ | delta method [1] |
| Coverage | share of nominal 95% intervals containing $\beta_D$ | binomial |
| $\tau_d^2$ bias, coverage | as above for the RS and RSL estimates and intervals | as above |
| Individual effects | RMSE and correlation of the BLUP with the true individual value | bootstrap over replicates |
| Failed fits, boundary fits | rates | binomial |
| Half-life of the largest $\lvert T\rvert$ | share of replicates in which each $h \in \mathcal T$ gives the largest $\lvert T(h)\rvert$, for F1 and VF, by cell | binomial, per level |
| AIC choice of $h$ | share of replicates in which AIC (ML, with the tested term) chooses each $h \in \mathcal T$, for F1 and VF, by cell | binomial, per level |
| Profile-interval width in $r_1$ | width of F1's profile-likelihood interval for $r_1 = 2^{-1/h} \in [0, 1)$, deviance cut 2.71 (boundary at $r_1 = 0$ [91]); share of replicates in which it covers all of $[0, 1)$ | $\mathrm{SD}/\sqrt{n_{\text{reps}}}$; binomial |

Table: Performance measures and their Monte Carlo standard errors

*Revised 2026-10-09 from docs/58, Section 12.3:* the last three rows
were added. The half-life measures describe what the data say about
$h$; the best-fitting $h$ is a working value that depends on the effect
size, not an estimate of the drug's carryover half-life (docs/58,
Section 9.6; for a true half-life of 1 week, the best-fitting $h$ was
2.4 weeks at $f_m = 0.4$ and 1.2 at $f_m = 1$, first-order closed
form). At the base alternative a trial of 70 is expected to say almost
nothing about $h$: no finite $h$ is more than 2 deviance units from the
best when the truth is short (docs/58, Table 9). The profile on $r_1$
needs a finer grid than $\mathcal T$ (docs/58 used 54 half-lives; Section
7.3). Bias and coverage of $\beta_D^{\text{exp}}$ for every $D_{bc}$
analysis are computed for $\hat\beta_D^{\text{exp}}(h)$ (Section 4.2).
The AIC-weighted average of $\hat\beta_D^{\text{exp}}(h)$ over
$\mathcal T$, with Buckland's unconditional standard error [88], is
reported as a descriptive summary, not as a test.

Failed fits are excluded from rates and counted, as in the stress test
(`docs/46`, Amendment 2, item 3); an analysis with more than 2% failed
fits in any cell is flagged.

### 6.2 The nominal band

At R1, a single null cell has MCSE 0.0049 and a 2-MCSE band of 0.040
to 0.060; the pooled rate over the 51 N-of-1 null cells of one analysis
across designs has MCSE about 0.0007. At R2 the single-cell band is
0.036 to 0.064, as in `docs/46`.

### 6.3 Paired comparisons

All analyses of a design are fitted to the same trials, and power
differences within a cell are tested by McNemar's test. Designs are not
paired: their schedules differ. Participant-level draws that do not
depend on the schedule ($\mathrm{BL}_i$, $v_i$, $\varepsilon_{i0}$) use
common random numbers across designs; design comparisons are otherwise
treated as independent proportions.

### 6.4 Working decision rule

Within each design, adapted from `docs/46`, Section 7, and stated as a
working rule that may be revised with recorded reasons:

1. **Size screen.** Discard an analysis whose pooled null rejection
   rate over that design's null cells lies outside 2 MCSE of 0.05, or
   that exceeds 0.075 in any null cell. Heterogeneous nulls (block 2)
   count; they test whether ignoring $\tau^2$ costs size.
   *Size screen at $N = 150$ (revised 2026-10-09 from docs/58, Section
   12.3).* An analysis with a nonzero closed-form null bias $b_0$ in a
   design (Section 5.1, Table "Closed-form null bias": F1 and S in the
   Hybrid and OL+BDC designs, S in CO) must in addition lie within the
   single-cell 2-MCSE band of Section 6.2 in that design's $N = 150$
   null cell of block 3, judged on its own and not only through the
   pooled rate or the 0.075 cap, because a fixed bias inflates size as
   $N$ grows (closed form for F1, Hybrid: 0.058 at $N = 70$, 0.068 at
   $N = 150$, above the upper band limits 0.060 at R1 and 0.064 at R2).
   If F1 fails this screen, VF is the candidate replacement reference
   (docs/58, Section 12.2; item 21).
2. **Estimation screen** (new). In the base alternative cells, flag an
   analysis whose coverage of $\beta_D$ is below 0.925 or whose absolute
   bias exceeds 10% of $|\beta_D|$. Flagged analyses may still be
   reported for testing, not for estimation.
3. **Oracle exclusion.** XO, OP and AN0 are never selected. *(Revised
   2026-10-09 from docs/58, Section 12.3.)* Nor is the negative control
   F1A. The randomization tests RV and RF are exact under the sharp
   null by construction; they are screened like the others, with their
   null rejection rates read as a check of the implementation, and
   their heterogeneous-null rates as in item 23.
4. **Selection.** Among survivors, select by minimax regret: the
   smallest largest shortfall in size-adjusted power over the
   alternative cells.
5. **Margin.** The proposed reference F1 is retained unless a survivor
   beats it by more than 0.05 in largest shortfall, with the McNemar
   difference significant at 0.05 in at least half the alternative
   cells. Ties within 0.02 go to the analysis with fewer covariance
   parameters.

**Across designs.** Power curves of each design's selected analysis
are compared at common $m_{BR}$ and $N$; differences are tested as
independent proportions. No design is "selected": the comparison
reports operating characteristics.

## 7. Implementation plan

### 7.1 What to reuse

- **The stress-test driver**
  `analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R`
  (inspected): its `build()` already gives anchored carryover and
  randomized allocation by default (lines 25 to 35, 118 to 133), the
  separable AR(1) construct, Weibull decay, dropout, measurement error,
  per-chunk seeds and checkpointing. It is the basis of the new driver.
- **The patches** `analysis/patches/orig-to-refined/` implement the
  refined DGP and fifteen analyses selectable by `op$analysis`, verified
  against the driver. They are the reference for reproducing the
  refined process in Hendrickson's code, but not the engine: their
  `generateSimulatedResults()` fits the analysis once to the whole
  population, which its README reports as impractical with the refined
  analyses (a single F1 fit on 7,200 participants did not finish in 10
  minutes), and their analyses test `bmc:Dbc` only.
- **Design definitions** from `hendrickson-problems/04-power-simulation.R`
  (lines 74 to 97).

### 7.2 What to add

1. **Closed-form module** (`01-closed-form.R`): estimand values per
   cell; probability limits and sandwich variances for every fixed-
   mean analysis; analytic power curves. Runs first; gates the grid.
2. **Drug-effect scale.** The $BR$, $TV$ and $PB$ ceilings are
   hard-coded in `build()` (10.98604 at line 120; 6.50647 at lines 134
   and 135); make them cell parameters.
3. **Random drug slope.** $v_i\,\Psi_{it}$ added after the draw, with
   $\Psi_{it}$ from the unit-ceiling profile.
4. **Designs.** OL+BDC and CO in the driver (11 uses Hybrid only,
   line 72); the RCT arm; the optional multi-crossover.
5. **Observed baseline with error** for all designs.
6. **Main-effect analyses** of Section 5.2, each returning the estimate,
   SE, df, $p$, interval, and for RS and RSL the variance component, its
   test and the BLUPs.
7. **Unit checks** before the run: the estimand values reproduce
   Section 4.1 at the published ceiling; the RCT closed form matches a
   large draw; $\tau^2$ is recovered from the generated individual
   effects; allocation shares; the closed form against a direct Monte
   Carlo of PW and PR.
8. **Half-life grid, maximum and randomization tests** *(revised
   2026-10-09 from docs/58, Section 12.3).* The fits over $\mathcal T$,
   the influence-based correlation and $c_{\max}$, the AIC choice, and
   the label-invariant randomization tests, adapted from
   `analysis/scripts/main-effect-half-life/58-common.R` (`perm_setup`,
   `perm_stats`), which reproduces the closed form of docs/58 in a
   Monte Carlo check of four Hybrid cells. With dropout the cheap
   randomization test needs per-pattern weights (docs/58, Section 14);
   in block 7 the refit version, or its omission, must be decided.

### 7.3 Compute estimate (inferred)

From `RERUN-anchored.md`, Section 3: the stress test costs about 11
worker-seconds per trial for about twelve analyses, dominated by three
unstructured Kenward-Roger fits; its post hoc runs cost about 0.4
worker-seconds per random-intercept fit with CR2. The list above has
nine random-intercept fits, two unstructured fits, one `mmrm` compound
symmetry fit, two random-slope fits and two summary measures, about 11
worker-seconds per N-of-1 trial if a random-slope fit costs 1 second
(assumed). At R1 that is about 165,000 N-of-1 trials, 500 worker-hours:
about 63 hours on 8 cores, 170 hours at simrun's default of 3 cores. R2
halves this. The $N = 150$ cells add perhaps 10% (assumed). RCT fits
are negligible. A pilot of 20 replicates per cell is run first to time
the fits, as in `docs/46`, Section 9.

*Revised 2026-10-09 from docs/58, Section 12.3 (inferred, not timed).*
The grid $\mathcal T$ in every cell adds four fits to each of S, F1,
VF, RS and XS: 16 random-intercept fits at about 0.4 s and 4
random-slope fits at about 1 s (assumed), about 10 worker-seconds per
N-of-1 trial, which nearly doubles the estimate above. The
randomization tests add about 2 s per trial (docs/58, Section 9.3,
timed on one machine); F1A and the maxima reuse the grid fits apart
from one ML refit per $h$ for the AIC choice. A fine profile in $r_1$
(Section 6.1) would add more: docs/58 timed 17 s per trial for 54
half-lives with F1 and VF under REML and ML (Section 10.4 there), so it
is either restricted to F1 under ML on a coarser grid or run in a
subset of cells. The pilot decides; restricting RS and XS to $h = 0.2$
is the main saving available (item 22).

### 7.4 Execution on the simulation server

Through `simrun`, as described in `RERUN-anchored.md`, Section 4.2:
the job runs the committed tree, so the driver must be committed first,
and `simrun submit` refuses while tracked files have uncommitted
changes; commands are argument vectors without redirection; the job's
`SIM_CORES` sets the worker count; results come back with `simrun
fetch` and `simrun apply`. One job per block keeps jobs short and
restartable. The core count the host can give is unknown (unverified).

### 7.5 Output naming

- Drivers: `analysis/scripts/treatment-main-effect/rebuild/`,
  `01-closed-form.R`, `02-simulate.R`, `03-decision.R`,
  `04-tables-figures.R`. The v1 drivers stay where they are.
- Data: `analysis/data/treatment-main-effect/rebuild/`:
  `closed-form-{estimands,plim,power}.csv`; `cells-{tag}/` with chunk
  files; `replicates-{tag}.rds`; `summary-{tag}.csv`;
  `decision-{screen,estimation,power,minimax}-{tag}.csv`. Tags `full`,
  `pilot20`, and a block suffix (`-b1` to `-b9`) for server jobs.
- Seeds per chunk, $20261009 + 1000 \times \text{cell} + \text{chunk}$,
  so results do not depend on the number of cores (the convention of
  `docs/47`, Section 11).

## 8. Open decisions for the user

1. **Estimand definition.** Adopt $\beta_D^{\text{exp}}$ (exposure-
   averaged on-minus-off drug response at the design's own visits) as
   primary, with the randomized blinded contrast as secondary for the
   Hybrid design, or choose the full-effect or the blinded contrast as
   primary (Section 4).
2. **$\tau^2$ grid.** $\tau \in \{0, 2.75, 5.49\}$ points (0, 0.25 and
   0.5 of the published ceiling), set in absolute units so that
   heterogeneous nulls exist; or a grid relative to $m_{BR}$.
3. **Effect-size grid units.** The multiplier $f_m$ on the drug-response
   ceiling, reported in points and in units of
   $\sigma_{Y \mid \mathrm{BL}} = 16.03$; or another SD (the unconditional
   24.46, or an RCT-based Cohen's $d$).
4. **Analysis list.** The twenty-four analyses of Section 5.2
   (nineteen before the revision of 2026-10-09 from docs/58, Section
   12.3, which added MX-F1, MX-VF, RV, RF and the negative control
   F1A), including the oracle arms XO, OP and AN0, the variant of
   Kenward-Roger, and the default $t_{1/2}^{a} = 0.2$ week; drop or add
   analyses.
5. **Multi-crossover.** Whether to keep the current paper's design as an
   additional design (bridge to v1, expected to have little power under
   re-titration), and if so in which blocks.
6. **Replicate counts.** R1 (2,000 null, 1,000 alternative; about 63
   hours on 8 cores) or R2 (1,000 and 500, as `docs/46`; about 32
   hours).
7. **RCT baseline.** Whether the RCT is analyzed by ANCOVA on a baseline
   measured with error, and how much error: reliability 1 (error-free,
   as v1), 0.9 or 0.8; and whether the same observed baseline is used
   in the N-of-1 designs (proposed: yes).

Further decisions that the protocol has had to assume:

8. **Base half-life** 0.2 week, with 0.1 optional and 1 week as the
   labeled behavioral scenario.
9. **Form of heterogeneity.** $v_i$ scaling the Gompertz profile
   (proposed), or scaling the binary drug state; independent of
   $\mathrm{BL}_i$ (proposed) or correlated with it.
10. **$BR$ noise in unexposed participants.** Keep the inherited SD 8
    off drug and in the RCT placebo arm (the refined process as
    specified), or make the $BR$ SD follow the drug effect.
11. **Decision rule.** The estimation screen of Section 6.4, item 2, is
    new relative to `docs/46`; its thresholds (coverage 0.925, bias
    10%) and the single-cell cap of 0.075 at R1 (5 MCSE) need
    confirming.

**From the literature review (Appendix A).** Proposals only; none is
applied to Sections 1 to 7, except item 12 as revised on 2026-10-09.
Reasons and sources are in Appendix A.5.

12. **Randomization test (L1).** *(Revised 2026-10-09 from docs/58,
    Section 12.3.)* Add the label-invariant randomization tests of
    Section 5.2, RV (VF statistic) and RF (F1 statistic): fixed-count
    permutation of the observed path labels, 999 relabelings, weights
    and CR0 studentization from the null model $w_t + \eta_t$ fitted
    once, and the maximum over $\mathcal T$ computed inside each
    relabeling, run in every cell rather than in a subset, since they
    cost about 2 seconds per trial (docs/58, Section 9.3). This replaces
    the earlier proposal (a re-draw of the path labels applied to the PR
    statistic in every cell and to F1 with refits in a subset of cells).
    RV and RF are entered in Section 5.2 as working proposals; whether
    either is the primary test is item 21.
13. **Bayesian hierarchical arm (L2).** Add a Bayesian random-slope
    analysis in blocks 1 to 3 at $N = 35$ and 70, judged by the same
    frequentist operating characteristics, or record its omission as a
    stated limitation.
14. **Onset exclusion by schedule (L3).** Add an analysis XT that drops
    the first on-drug visit after each start or restart of the drug,
    the onset counterpart of XS.
15. **Literature bridge cell (L4).** Add a reduced cell (constant
    effect, no $PB$ or $TV$, $t_{1/2} = 0$, AR(1) only) in which PW and
    the random-slope model can be checked against the published closed
    forms of Senn [2], Araujo et al. [11] and Yang et al. [16].
16. **Heterogeneous null with carryover (L5).** Cross block 2's null
    cells with the behavioral half-life $t_{1/2} = 1$ week.
17. **Two-point heterogeneity (L6).** Add one cell in which $v_i$ takes
    two values (responders and non-responders) with the same $\tau^2$
    as a normal cell; this extends item 9.
18. **RCT comparator for the design comparison (L7).** Report AM, which
    uses every visit, as the information-matched comparator in the
    across-design comparison, and AN as the conventional comparator.
19. **Additional measures (L8).** Add interval width for $\beta_D$,
    coverage of individual-effect intervals, and a per-participant OLS
    comparator for the BLUPs.
20. **Bounded outcome (L9).** Add one sensitivity cell in which $Y$ is
    rounded and truncated to the CAPS range.

**From the main-effect half-life note** *(revised 2026-10-09 from
docs/58, Section 12.3; reasons in docs/58, Sections 7, 9 and 12.2).*

21. **Model-based or randomization-based primary test.** Whether paper
    04's primary test of $\beta_D$ is model-based (F1, closed-form power
    0.63 at the base alternative, Hybrid, $N = 70$) or
    randomization-based (RV, power about 0.39), or both are reported.
22. **Half-life policy.** Whether the analysis half-life is fixed at
    $t_{1/2}^{a} = 0.2$ week, or chosen from the data inside a
    calibrated test (RV, MX-F1, MX-VF). F1A with a single-test
    reference is not an option (Section 5.2). If the half-life is fixed,
    RS and XS can be fitted at $h = 0.2$ only (Section 7.3).
23. **Randomization tests under the weak null.** Whether the
    heterogeneous-null cells (block 2) judge the randomization tests in
    their studentized form only. An unstudentized permutation test is
    not exact under the weak null; its closed-form size in the Hybrid
    design is 0.053 at $\tau = 5.49$ (docs/58, Proposition 7).
24. **Compounding carryover.** Whether any main-effect result is quoted
    from runs with the compounding carryover recursion of the published
    process (docs/58, Section 7, recommends not, at $t_{1/2} \ge 1$
    week). At $t_{1/2} = 1$ week that recursion inflates
    $\beta_D^{\text{exp}}$ by 11% (Hybrid) and 14% (OL+BDC) and F1's
    power by 0.047 and 0.068, against anchored carryover.

## 9. What carries over from the current paper

**Carried over, with edits** (`report.Rmd`):

- The Introduction's background on N-of-1 designs and sample size
  (lines 76 to 141), including the Senn and Zucker variance-component
  apparatus (lines 142 to 252), which becomes directly relevant once
  $\tau^2$ is in the DGP; references [2, 3, 10, 11] stay.
- The single-case framework, prazosin, regulatory and literature
  sections (lines 290 to 458), with the comparisons of Blackston et al.
  [12] and others, rewritten where they quote the paper's own results.
- The comparison with prior simulation studies in the Discussion
  (lines 1770 to 1922), kept as structure, with every number replaced.

**Replaced.** The abstract; the Objectives (line 460); all of Methods
(lines 472 to 600); all Results (lines 1053 to 1573); the remaining
Discussion and the Conclusions; the cover letter, title page, short
report and bullets. The word "hybrid" is removed from every statement
about the v1 simulation.

**Bibliography fix found while reading.** The key
`araujoComparisonFourMethods2014` in the paper's `references.bib`
(line 22) holds Chen and Chen (2014); the key should be renamed
(inspected).

## 10. References

1. Morris TP, White IR, Crowther MJ. Using simulation studies to
   evaluate statistical methods. *Statistics in Medicine*
   2019;38(11):2074-2102. doi:10.1002/sim.8086.
2. Senn S. Sample size considerations for n-of-1 trials. *Statistical
   Methods in Medical Research* 2019;28(2):372-383.
   doi:10.1177/0962280217726801.
3. Zucker DR, Ruthazer R, Schmid CH. Individual (N-of-1) trials can be
   combined to give population comparative treatment effect estimates:
   methodologic considerations. *Journal of Clinical Epidemiology*
   2010;63(12):1312-1323. doi:10.1016/j.jclinepi.2010.04.020.
4. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020;2:13. doi:10.3389/fdgth.2020.00013.
5. Bell RM, McCaffrey DF. Bias reduction in standard errors for linear
   regression with multi-stage samples. *Survey Methodology*
   2002;28(2):169-181.
6. Pustejovsky JE, Tipton E. Small-sample methods for cluster-robust
   variance estimation and hypothesis testing in fixed effects models.
   *Journal of Business and Economic Statistics* 2018;36(4):672-683.
   doi:10.1080/07350015.2016.1247004.
7. Jones B, Kenward MG. *Design and Analysis of Cross-Over Trials*, 3rd
   edition. Boca Raton: CRC Press, 2014. doi:10.1201/b17537.
8. Kenward MG, Roger JH. Small sample inference for fixed effects from
   restricted maximum likelihood. *Biometrics* 1997;53(3):983-997.
   doi:10.2307/2533558.
9. Stram DO, Lee JW. Variance components testing in the longitudinal
   mixed effects model. *Biometrics* 1994;50(4):1171-1177.
10. Senn S. Mastering variation: variance components and personalised
    medicine. *Statistics in Medicine* 2016;35(7):966-977.
    doi:10.1002/sim.6739.
11. Araujo A, Julious S, Senn S. Understanding variation in sets of
    N-of-1 trials. *PLoS ONE* 2016;11(12):e0167167.
    doi:10.1371/journal.pone.0167167.
12. Blackston JW, Chapple AG, McGree JM, McDonald S, Nikles J.
    Comparison of aggregated N-of-1 trials with parallel and crossover
    randomized controlled trials using simulation studies. *Healthcare*
    2019;7(4):137. doi:10.3390/healthcare7040137.
13. Raskind MA, Peterson K, Williams T, et al. A trial of prazosin for
    combat trauma PTSD with nightmares in active-duty soldiers returned
    from Iraq and Afghanistan. *American Journal of Psychiatry*
    2013;170(9):1003-1010. doi:10.1176/appi.ajp.2013.12081133.
14. Weathers FW, Keane TM, Davidson JR. Clinician-administered PTSD
    scale: a review of the first ten years of research. *Depression and
    Anxiety* 2001;13(3):132-156. doi:10.1002/da.1029.
15. pmsimstats team. `docs/37`, `docs/39`, `docs/40`, `docs/41`,
    `docs/42`, `docs/45`, `docs/46`, `docs/47`, `docs/49`, `docs/51`,
    `docs/53`, `docs/56`, `docs/57`, `docs/58`;
    `analysis/report/compendium-mathematics.tex`; paper 01.
16. Yang J, Steingrimsson JA, Schmid CH. Sample size calculations for
    n-of-1 trials. arXiv preprint 2110.08970, 2021.
    doi:10.48550/arXiv.2110.08970.
17. Hedges LV, Shadish WR, Natesan Batley P. Power analysis for
    single-case designs: computations for (AB)^k designs. *Behavior
    Research Methods* 2023;55(7):3494-3503.
    doi:10.3758/s13428-022-01971-9.
18. Hedges LV, Pustejovsky JE, Shadish WR. A standardized mean
    difference effect size for single case designs. *Research Synthesis
    Methods* 2012;3(3):224-239. doi:10.1002/jrsm.1052.
19. Shi D, Ye T. Behavioral carry-over effect and power consideration in
    crossover trials. *Biometrics* 2024;80(2):ujae023.
    doi:10.1093/biomtc/ujae023.
20. Carrozzo AE, Zimmermann G, Bathke AC, Neunhaeuserer D, Niebauer J,
    Kulnik ST. Two-arm crossover randomized controlled trial versus
    meta-analysis of N-of-1 studies: comparison of statistical
    efficiency in determining an intervention effect. *Biometrical
    Journal* 2025;67(2):e70045. doi:10.1002/bimj.70045.
21. Chen X, Chen P. A comparison of four methods for the analysis of
    N-of-1 trials. *PLoS ONE* 2014;9(2):e87752.
    doi:10.1371/journal.pone.0087752.
22. Gärtner T, Schneider J, Arnrich B, Konigorski S. Comparison of
    Bayesian networks, G-estimation and linear models to estimate causal
    treatment effects in aggregated N-of-1 trials with carry-over
    effects. *BMC Medical Research Methodology* 2023;23:191.
    doi:10.1186/s12874-023-02012-5.
23. Krone T, Boessen R, Bijlsma S, van Stokkum R, Clabbers NDS, Pasman
    WJ. The possibilities of the use of N-of-1 and do-it-yourself trials
    in nutritional research. *PLoS ONE* 2020;15(5):e0232680.
    doi:10.1371/journal.pone.0232680.
24. Tang J, Landes RD. Some t-tests for N-of-1 trials with serial
    correlation. *PLoS ONE* 2020;15(2):e0228077.
    doi:10.1371/journal.pone.0228077.
25. Wang Y, Schork NJ. Power and design issues in crossover-based N-of-1
    clinical trials with fixed data collection periods. *Healthcare
    (Basel)* 2019;7(3):84. doi:10.3390/healthcare7030084.
26. Lin SX, Morrison L, Smith PWF, Hargood C, Weal M, Yardley L.
    Properties of bootstrap tests for N-of-1 studies. *British Journal
    of Mathematical and Statistical Psychology* 2016;69(3):276-290.
    doi:10.1111/bmsp.12071.
27. Percha B, Baskerville EB, Johnson M, Dudley JT, Zimmerman N.
    Designing robust N-of-1 studies for precision medicine: simulation
    study and design recommendations. *Journal of Medical Internet
    Research* 2019;21(4):e12641. doi:10.2196/12641.
28. Schork NJ. Accommodating serial correlation and sequential design
    elements in personalized studies and aggregated personalized
    studies. *Harvard Data Science Review* 2022;(Special Issue 3).
    doi:10.1162/99608f92.f1eef6f4.
29. Liao Z, Qian M, Kronish IM, Cheung YK. Analysis of N-of-1 trials
    using Bayesian distributed lag model with autocorrelated errors.
    *Statistics in Medicine* 2023;42(13):2044-2060. doi:10.1002/sim.9676.
30. Selukar S, Prince DK, May S. A framework for sequential monitoring of
    individual N-of-1 trials and combining results across a series of
    sequentially monitored N-of-1 trials. *Clinical Trials*
    2025;22(3):257-266. doi:10.1177/17407745241304284.
31. Jiang S, Arnold SE, Betensky RA. Design and analysis of N-of-1
    trials that incorporate sequential monitoring. *Statistics in
    Medicine* 2025;44(15-17):e70177. doi:10.1002/sim.70177.
32. Senarathne SGJ, Overstall AM, McGree JM. Bayesian adaptive N-of-1
    trials for estimating population and individual treatment effects.
    *Statistics in Medicine* 2020;39(29):4499-4518. doi:10.1002/sim.8737.
33. Shrestha S, Jain S. A Bayesian-bandit adaptive design for N-of-1
    clinical trials. *Statistics in Medicine* 2021;40(7):1825-1844.
    doi:10.1002/sim.8873.
34. Giai J, Peron J, Roustit M, Cracowski JL, Roy P, Ozenne B, Buyse M,
    Maucort-Boulch D. Individualized Net Benefit estimation and
    meta-analysis using generalized pairwise comparisons in N-of-1
    trials. *Statistics in Medicine* 2023;42(6):878-893.
    doi:10.1002/sim.9648.
35. Cruz NA, Melo OO, Mylona K. Penalized GEE for complex carry-over in
    repeated-measures crossover designs. *Statistics in Medicine*
    2026;45(10-12):e70561. doi:10.1002/sim.70561.
36. Ferron JM, Bell BA, Hess MR, Rendina-Gobioff G, Hibbard ST. Making
    treatment effect inferences from multiple-baseline data: the utility
    of multilevel modeling approaches. *Behavior Research Methods*
    2009;41(2):372-384. doi:10.3758/BRM.41.2.372.
37. Ferron JM, Farmer JL, Owens CM. Estimating individual treatment
    effects from multiple-baseline data: a Monte Carlo study of
    multilevel-modeling approaches. *Behavior Research Methods*
    2010;42(4):930-943. doi:10.3758/BRM.42.4.930.
38. Moeyaert M, Ugille M, Ferron JM, Beretvas SN, Van den Noortgate W.
    Modeling external events in the three-level analysis of
    multiple-baseline across-participants designs: a simulation study.
    *Behavior Research Methods* 2013;45(2):547-559.
    doi:10.3758/s13428-012-0274-1.
39. Petit-Bois M, Baek EK, Van den Noortgate W, Beretvas SN, Ferron JM.
    The consequences of modeling autocorrelation when synthesizing
    single-case studies using a three-level model. *Behavior Research
    Methods* 2016;48(2):803-812. doi:10.3758/s13428-015-0612-1.
40. Declercq L, Jamshidi L, Fernández-Castilla B, Beretvas SN, Moeyaert
    M, Ferron JM, Van den Noortgate W. Analysis of single-case
    experimental count data using the linear mixed effects model: a
    simulation study. *Behavior Research Methods* 2019;51(6):2477-2497.
    doi:10.3758/s13428-018-1091-y.
41. Declercq L, Jamshidi L, Fernández Castilla B, Moeyaert M, Beretvas
    SN, Ferron JM, Van den Noortgate W. Multilevel meta-analysis of
    individual participant data of single-case experimental designs:
    one-stage versus two-stage methods. *Multivariate Behavioral
    Research* 2022;57(2-3):298-317. doi:10.1080/00273171.2020.1822148.
42. Moeyaert M, Yang P, Xu X. The power to explain variability in
    intervention effectiveness in single-case research using
    hierarchical linear modeling. *Perspectives on Behavior Science*
    2022;45(1):13-35. doi:10.1007/s40614-021-00304-z.
43. Chen M, Pustejovsky JE. Multilevel meta-analysis of single-case
    experimental designs using robust variance estimation.
    *Psychological Methods* 2024;29(3):537-560. doi:10.1037/met0000510.
44. Bouwmeester S, Jongerling J. Power of a randomization test in a
    single case multiple baseline AB design. *PLoS ONE*
    2020;15(2):e0228355. doi:10.1371/journal.pone.0228355.
45. Michiels B, Tanious R, De TK, Onghena P. A randomization test wrapper
    for synthesizing single-case experiments using multilevel models: a
    Monte Carlo simulation study. *Behavior Research Methods*
    2020;52(2):654-666. doi:10.3758/s13428-019-01266-6.
46. Levin JR, Ferron JM, Kratochwill TR. Nonparametric statistical tests
    for single-case systematic and randomized ABAB...AB and alternating
    treatment intervention designs: new developments, new directions.
    *Journal of School Psychology* 2012;50(5):599-624.
    doi:10.1016/j.jsp.2012.05.001.
47. Tueller S, Ramirez D, Cance JD, Ye A, Wheeler AC, Fan Z, et al.
    Power analysis for idiographic (within-subject) clinical trials:
    implications for treatments of rare conditions and precision
    medicine. *Behavior Research Methods* 2023;55(8):4175-4199.
    doi:10.3758/s13428-022-02012-1.
48. Natesan Batley P, Hedges LV. Accurate models vs. accurate estimates:
    a simulation study of Bayesian single-case experimental designs.
    *Behavior Research Methods* 2021;53(4):1782-1798.
    doi:10.3758/s13428-020-01522-0.
49. Joo SH, Ferron JM. Application of the within- and between-series
    estimators to non-normal multiple-baseline data: maximum likelihood
    and Bayesian approaches. *Multivariate Behavioral Research*
    2019;54(5):666-689. doi:10.1080/00273171.2018.1564877.
50. Zucker DR, Schmid CH, McIntosh MW, D'Agostino RB, Selker HP, Lau J.
    Combining single patient (N-of-1) trials to estimate population
    treatment effects and to evaluate individual patient responses to
    treatment. *Journal of Clinical Epidemiology* 1997;50(4):401-410.
    doi:10.1016/S0895-4356(96)00429-5.
51. Stunnenberg BC, Woertman W, Raaphorst J, Statland JM, Griggs RC,
    Timmermans J, et al. Combined N-of-1 trials to investigate
    mexiletine in non-dystrophic myotonia using a Bayesian approach;
    study rationale and protocol. *BMC Neurology* 2015;15:43.
    doi:10.1186/s12883-015-0294-4.
52. Stunnenberg BC, Raaphorst J, Groenewoud HM, Statland JM, Griggs RC,
    Woertman W, et al. Effect of mexiletine on muscle stiffness in
    patients with nondystrophic myotonia evaluated using aggregated
    N-of-1 trials. *JAMA* 2018;320(22):2344-2353.
    doi:10.1001/jama.2018.18020.
53. Nikles J, Mitchell GK, Schluter P, Good P, Hardy J, Rowett D, et al.
    Aggregating single patient (n-of-1) trials in populations where
    recruitment and retention was difficult: the case of palliative
    care. *Journal of Clinical Epidemiology* 2011;64(5):471-480.
    doi:10.1016/j.jclinepi.2010.05.009.
54. Schmid C, Yang J. Bayesian models for N-of-1 trials. *Harvard Data
    Science Review* 2022;(Special Issue 3).
    doi:10.1162/99608f92.3f1772ce.
55. Defelippe V, Bartos F, Wagenmakers EJ, Braun KPJ, Jansen FE, Otte
    WM. Continuous outcome estimation in N-of-1 trials for accelerated
    decision-making. *Epilepsia* 2026;67(5):2254-2269.
    doi:10.1002/epi.70119.
56. Diaz FJ. Using population crossover trials to improve the decision
    process regarding treatment individualization in N-of-1 trials.
    *Statistics in Medicine* 2021;40(20):4345-4361. doi:10.1002/sim.9030.
57. Sturdevant SG, Lumley T. Statistical methods for testing carryover
    effects: a mixed effects model approach. *Contemporary Clinical
    Trials Communications* 2021;22:100711.
    doi:10.1016/j.conctc.2021.100711.
58. Cheung K, Mitsumoto H. Evaluating personalized (N-of-1) trials in
    rare diseases: how much experimentation is enough? *Harvard Data
    Science Review* 2022;(Special Issue 3).
    doi:10.1162/99608f92.e11adff0.
59. Moeyaert M, Ugille M, Ferron JM, Beretvas SN, Van den Noortgate W.
    The three-level synthesis of standardized single-subject
    experimental data: a Monte Carlo simulation study. *Multivariate
    Behavioral Research* 2013;48(5):719-748.
    doi:10.1080/00273171.2013.816621.
60. Xie T, Yu Z. N-of-1 design and its applications to personalized
    treatment studies. *Statistics in Biosciences* 2017;9(2):662-675.
    doi:10.1007/s12561-016-9165-9.
61. Swan DM, Pustejovsky JE. A gradual effects model for single-case
    designs. *Multivariate Behavioral Research* 2018;53(4):574-593.
    doi:10.1080/00273171.2018.1466681.
62. Senn S. *Cross-over Trials in Clinical Research*, 2nd edition.
    Chichester: Wiley, 2002. doi:10.1002/0470854596.
63. Grizzle JE. The two-period change-over design and its use in
    clinical trials. *Biometrics* 1965;21(2):467-480.
    doi:10.2307/2528104.
64. Freeman PR. The performance of the two-stage analysis of
    two-treatment, two-period crossover trials. *Statistics in
    Medicine* 1989;8(12):1421-1432. doi:10.1002/sim.4780081202.
65. Fleiss JL. A critique of recent research on the two-treatment
    crossover design. *Controlled Clinical Trials* 1989;10(3):237-243.
    doi:10.1016/0197-2456(89)90065-2.
66. Senn S. The AB/BA crossover: past, present and future? *Statistical
    Methods in Medical Research* 1994;3(4):303-324.
    doi:10.1177/096228029400300402.
67. Kenward MG, Jones B. The analysis of data from 2 x 2 cross-over
    trials with baseline measurements. *Statistics in Medicine*
    1987;6(8):911-926. doi:10.1002/sim.4780060806.
68. Senn S. Viewpoint: do not resurrect the two-stage procedure.
    *Pharmaceutical Statistics* 2022;21(4):808-814.
    doi:10.1002/pst.2224.
69. Kabaila P, Vicendese M. The performance of a two-stage analysis of
    ABAB/BABA crossover trials. *Biometrical Journal*
    2012;54(3):361-369. doi:10.1002/bimj.201100201.
70. Sheiner LB, Stanski DR, Vozeh S, Miller RD, Ham J. Simultaneous
    modeling of pharmacokinetics and pharmacodynamics: application to
    d-tubocurarine. *Clinical Pharmacology and Therapeutics*
    1979;25(3):358-371. doi:10.1002/cpt1979253358.
71. Holford NHG, Sheiner LB. Understanding the dose-effect
    relationship: clinical application of
    pharmacokinetic-pharmacodynamic models. *Clinical
    Pharmacokinetics* 1981;6(6):429-453.
    doi:10.2165/00003088-198106060-00002.
72. Gasparrini A, Armstrong B, Kenward MG. Distributed lag non-linear
    models. *Statistics in Medicine* 2010;29(21):2224-2234.
    doi:10.1002/sim.3940.
73. Davies RB. Hypothesis testing when a nuisance parameter is present
    only under the alternative. *Biometrika* 1977;64(2):247-254.
    doi:10.1093/biomet/64.2.247.
74. Davies RB. Hypothesis testing when a nuisance parameter is present
    only under the alternative. *Biometrika* 1987;74(1):33-43.
    doi:10.1093/biomet/74.1.33.
75. Andrews DWK, Ploberger W. Optimal tests when a nuisance parameter
    is present only under the alternative. *Econometrica*
    1994;62(6):1383-1414. doi:10.2307/2951753.
76. Hansen BE. Inference when a nuisance parameter is not identified
    under the null hypothesis. *Econometrica* 1996;64(2):413-430.
    doi:10.2307/2171789.
77. Edgington ES, Onghena P. *Randomization Tests*, 4th edition. Boca
    Raton: Chapman and Hall/CRC, 2007. doi:10.1201/9781420011814.
78. Onghena P, Edgington ES. Customization of pain treatments:
    single-case design and analysis. *Clinical Journal of Pain*
    2005;21(1):56-68. doi:10.1097/00002508-200501000-00007.
79. Kratochwill TR, Levin JR. Enhancing the scientific credibility of
    single-case intervention research: randomization to the rescue.
    *Psychological Methods* 2010;15(2):124-144.
    doi:10.1037/a0017736.
80. Chung E, Romano JP. Exact and asymptotically robust permutation
    tests. *Annals of Statistics* 2013;41(2):484-507.
    doi:10.1214/13-AOS1090.
81. Wu J, Ding P. Randomization tests for weak null hypotheses in
    randomized experiments. *Journal of the American Statistical
    Association* 2021;116(536):1898-1913.
    doi:10.1080/01621459.2020.1750415.
82. Janssen A. Studentized permutation tests for non-i.i.d. hypotheses
    and the generalized Behrens-Fisher problem. *Statistics and
    Probability Letters* 1997;36(1):9-21.
    doi:10.1016/S0167-7152(97)00043-6.
83. DiCiccio CJ, Romano JP. Robust permutation tests for correlation
    and regression coefficients. *Journal of the American Statistical
    Association* 2017;112(519):1211-1220.
    doi:10.1080/01621459.2016.1202117.
84. Ding P, Feller A, Miratrix L. Randomization inference for
    treatment effect variation. *Journal of the Royal Statistical
    Society, Series B* 2016;78(3):655-671. doi:10.1111/rssb.12124.
85. Berger RL, Boos DD. P values maximized over a confidence set for
    the nuisance parameter. *Journal of the American Statistical
    Association* 1994;89(427):1012-1016.
    doi:10.1080/01621459.1994.10476836.
86. Leeb H, Pötscher BM. Model selection and inference: facts and
    fiction. *Econometric Theory* 2005;21(1):21-59.
    doi:10.1017/S0266466605050036.
87. Berk R, Brown L, Buja A, Zhang K, Zhao L. Valid post-selection
    inference. *Annals of Statistics* 2013;41(2):802-837.
    doi:10.1214/12-AOS1077.
88. Buckland ST, Burnham KP, Augustin NH. Model selection: an integral
    part of inference. *Biometrics* 1997;53(2):603-618.
    doi:10.2307/2533961.
89. Hjort NL, Claeskens G. Frequentist model average estimators.
    *Journal of the American Statistical Association*
    2003;98(464):879-899. doi:10.1198/016214503000000828.
90. Giacomini R, Politis DN, White H. A warp-speed method for
    conducting Monte Carlo experiments involving bootstrap estimators.
    *Econometric Theory* 2013;29(3):567-589.
    doi:10.1017/S0266466612000655.
91. Self SG, Liang KY. Asymptotic properties of maximum likelihood
    estimators and likelihood ratio tests under nonstandard
    conditions. *Journal of the American Statistical Association*
    1987;82(398):605-610. doi:10.1080/01621459.1987.10478472.
92. Ernst MD. Permutation methods: a basis for exact inference.
    *Statistical Science* 2004;19(4):676-685.
    doi:10.1214/088342304000000396.
93. Rosenberger WF, Uschner D, Wang Y. Randomization: the forgotten
    component of the randomized clinical trial. *Statistics in
    Medicine* 2019;38(1):1-12. doi:10.1002/sim.7901.
94. Seber GAF, Wild CJ. *Nonlinear Regression*. New York: Wiley, 1989.
    doi:10.1002/0471725315.
95. Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS*. New
    York: Springer, 2000. doi:10.1007/b98882.

References 2, 3, 8 to 14 were checked against PubMed records on
2026-10-08 (titles, volumes, pages, DOIs; Stram and Lee has no DOI in
its record). References 1, 4 to 6 are as cited in `docs/46` and
`docs/47`; reference 7 is a book, not checked. The CAPS reliability
values of reference 14 were not read from the article; only the
record's statement that reliability is high was checked.

References 16 to 61 were added on 2026-10-08 with Appendix A. Each was
checked against at least one of a PubMed record, Crossref, an arXiv
abstract page or the Zotero record; the source checked and the extent
of reading (full text or abstract only) are listed in Appendix A.1,
Table "Studies screened in full and how each was verified". Reference
16 is a preprint for which no journal version was found; reference 60
was checked against Crossref only and not read.

*Revised 2026-10-09 from docs/58, Section 16:* references 62 to 95
were added, in the order of the docs/58 reference list, which they
reproduce (docs/58 references 5, 7 to 15, 18 to 25, 28 to 42 and 45;
its references 1 to 4, 6, 16, 17, 26, 27, 43 and 44 are this list's 4,
1, 5, 6, 7, 29, 27, 45, 44, 2 and 11). They were not re-checked for
this revision. docs/58 checked 62 to 94 against at least one of
PubMed, the Crossref record of the DOI, the publisher page or the
Zotero item; reference 95 is as verified in `docs/57`. docs/58 records
the extent of reading of each (full text, abstract or record only);
statements here attributed to a source rest on that reading. The DOI
of reference 7 was added from docs/58.

## Changes to the working protocol

Changes made after results are seen will be appended as numbered
amendments, each with its date, scope and reason. Revisions made before
any result of this study exists are recorded here as numbered
revisions.

**Revision 1 (2026-10-09).** *Scope:* the edits proposed in
`docs/58-main-effect-testing-and-carryover-half-life.md`, Section 12.3,
approved by the user on 2026-10-09, and the notation decisions of the
same date. *Reason:* the closed-form and Monte Carlo development of
docs/58 showed that a fixed analysis half-life does not by itself
guarantee size for the main effect (F1 has a null bias of $+0.33$
points in the Hybrid design), that choosing the half-life by AIC with
the tested term is a max-$|t|$ test with a single-test critical value
(size 0.088), and that a randomization test can contain the choice at
about 2 seconds per trial. *Edits:* notation index, collisions and
units (Front matter); misspecified half-life plan (Section 3.6); scale
of a $D_{bc}$ coefficient (Section 4.2); working covariance, null bias
and mean-misfit sandwich term, with Table "Closed-form null bias"
(Section 5.1); $w_t$ and $\beta_T$ in the mean models, analyses MX-F1,
MX-VF, RV, RF and the negative control F1A, the grid $\mathcal T$ in
every cell (Section 5.2); half-life measures (Section 6.1); size screen
at $N = 150$ and exclusion of F1A (Section 6.4); implementation and
compute notes (Sections 7.2 and 7.3); items 4 and 12 revised and items
21 to 24 added (Section 8); references 62 to 95 (Section 10);
cross-references in Appendix A. Every edit is marked in place with
"Revised 2026-10-09". The study remains exploratory, and every revised
item remains a working proposal.

**Revision 2 (2026-10-09).** *Scope:* notation only. *Reason:* the
unit drug-response profile was written $q_{it}$ ($\mathbf q_p$ for a
path's vector, $\bar q_d$ for its on-minus-off mean, $q_{20}$ at week
20), which collides with the canonical path index $q$ of `NOTATION.md`.
*Edits:* the profile is now $\Psi_{it}$, $\boldsymbol\Psi_p$,
$\bar\Psi_d$ and $\Psi_{20}$ throughout (notation index; Sections 3.1,
3.2, 4.3, 5.1, 5.2 and 7.2), with the same rename in `docs/58`, and
$\Psi$ is added to `NOTATION.md`. No definition, value or result
changes.

```{=latex}
\clearpage
```

## Appendix A. Prior simulation and analytic approaches to testing the treatment main effect

This appendix describes how earlier work has evaluated tests and
estimators of an average treatment effect in N-of-1, aggregated N-of-1,
multi-crossover and single-case designs, and compares that work with
Sections 3 to 6. It changes nothing in Sections 1 to 10; its proposals
are collected in A.5 and entered in Section 8 as items 12 to 20.
*(Revised 2026-10-09 from docs/58, Section 12.3:* item 12, the
randomization test, has since been revised and entered in Section 5.2
as RV and RF; the appendix text below is otherwise as written on
2026-10-08.)

**Notation in this appendix.** Each source's model is restated in the
protocol's symbols where that is possible: $Y_{it}$ the outcome,
$\beta_D$ the population mean treatment effect, $v_i \sim N(0, \tau^2)$
the participant-by-treatment deviation, $\rho$ an AR(1) coefficient,
$t_{1/2}$ a carryover half-life and $N$ the number of participants. Two
symbols are local to the appendix: $\sigma_w^2$, the within-participant
residual variance of one occasion, and $k$, the number of AB cycles per
participant; $m$ is the number of observations in a single-person
series. A source's own symbol is given once, in parentheses, when it
differs. Several sources use $\tau$ or $\tau^2$ for something other
than treatment heterogeneity (a treatment effect, or a random-intercept
variance); this is pointed out where it matters, because it is the
most common source of misreading in this literature.

**Evidence labels.** As in the protocol, with one addition. *Quoted*:
the source's words. *Inspected*: read in the source on 2026-10-08.
*Derived*: follows algebraically from the source's model. *Argued*: a
judgment of this team about relevance to the Hybrid setting, not a
finding of the source. Unless a statement is labeled otherwise,
statements about a source's design, model and results are inspected,
from the reading named for that source in Table "Studies screened in
full and how each was verified".

### A.1 Search

**Sources and dates.** All searches were run on 2026-10-08.

- **Zotero, collection `pmsimstats-ng`** (key L29UF4HU): listed in full,
  391 entries in two pages (200 and 191), including duplicate records
  of the same paper. The library has no collection named `nof1`; a
  library-wide keyword search for "N-of-1" returned 84 records, all but
  a handful already in the collection. The library's semantic search
  failed (an embedding-service authorization error) and was not used.
- **PubMed**, eight queries: N-of-1 (title) with simulation or Monte
  Carlo and power or type I error; series, aggregated, combined or
  meta-analysis of N-of-1 with simulation; single-case with multilevel
  or hierarchical and simulation; Ferron, Moeyaert or Van den Noortgate
  with single-case or multiple-baseline and simulation; N-of-1 or
  multiple crossover with sample size or power and autocorrelation,
  serial correlation or carryover; N-of-1 with Bayesian hierarchical
  and simulation or operating characteristics; Zucker, Nikles,
  Stunnenberg or Schmid with N-of-1; multiple, multi-period, replicate
  or repeated crossover with simulation. The queries returned 119
  records, with overlap.
- **Web search** for journals outside PubMed (statistics and
  psychometrics journals, arXiv, Harvard Data Science Review), three
  queries, and the reference lists of every paper read in full.
- **Repository sources**: paper 04's literature sections
  (`report.Rmd`, lines 142 to 458 and 1770 to 1922) and
  `references.bib`; `docs/51` (Section 2 and its reference list);
  `docs/references/nof1-simulation.bib`; `docs/29`.

**Inclusion.** A study was included in A.2 and A.3 if it evaluated, by
Monte Carlo simulation or in closed form, a test or estimator of an
average (or individual) treatment effect in an N-of-1, aggregated
N-of-1, multi-crossover, AB/BA crossover with repeated measures, or
single-case design, and reported at least one of size, power, bias,
coverage or precision. Papers that only applied a method to one data
set, tutorials, protocols with a single sample-size calculation,
reviews and studies of omics "N-of-1" signals were excluded from the
comparison and, where useful, are listed as background in A.2.10.
Hendrickson et al. [4], whose target is the biomarker interaction, is
described because its data-generating process is the ancestor of the
refined process.

**Counts.** About 520 records were screened on title or abstract
(391 Zotero entries, 119 PubMed records, and the web and reference-list
leads, with heavy overlap between sources). Fifty-one were assessed in
full or in abstract: 38 were included in A.2 and A.3, and 13 are kept
as background (A.2.10). By where they are held, 26 of the 38 included
studies are in the `pmsimstats-ng` collection; 11 were found through
PubMed and 1 through the repository's
`docs/references/nof1-simulation.bib`, and these 12 are not in the
collection. Further leads from reference lists (for example on
randomization tests for reversal designs, gradual-effect models and
one-stage IPD meta-analysis) were identified but not assessed; they
are listed in the report that accompanies this revision, not here.

**How each item was verified.** Reading was done by six parallel
extraction passes in this session, each told to verify every record,
read full text where it was available, and say where only an abstract
was read. Full text from PubMed Central loses inline mathematics; where
an equation was reconstructed from surrounding text or a publisher
page, the extraction said so, and such equations are given below in
words or flagged.

| Ref | Study | Record checked | Reading |
|---|---|---|---|
| [2] | Senn 2019 | PubMed, Zotero | full text (Zotero PDF) |
| [11] | Araujo, Julious, Senn 2016 | PubMed | full text (PLOS, PMC) |
| [16] | Yang, Steingrimsson, Schmid 2021 | arXiv | full text (Zotero PDF, v1) |
| [17] | Hedges, Shadish, Natesan Batley 2023 | Crossref, Zotero | full text (Zotero PDF) |
| [18] | Hedges, Pustejovsky, Shadish 2012 | PubMed | full text (ERIC copy) |
| [19] | Shi, Ye 2024 | PubMed | full text (arXiv v2); supplement not read |
| [12] | Blackston et al. 2019 | PubMed | full text (PMC) |
| [20] | Carrozzo et al. 2025 | PubMed | full text (PMC) |
| [21] | Chen, Chen 2014 | PubMed | full text (PMC, PLOS); 4-cycle supplement not read |
| [22] | Gärtner et al. 2023 | PubMed | full text (PMC); supplement not read |
| [23] | Krone et al. 2020 | PubMed | full text (PMC) |
| [24] | Tang, Landes 2020 | PubMed | full text (PMC); figures not read |
| [25] | Wang, Schork 2019 | PubMed | full text (PMC); power figures not read |
| [26] | Lin et al. 2016 | PubMed | full text (PMC); appendix tables not read |
| [27] | Percha et al. 2019 | PubMed | full text (PMC); equations partly reconstructed |
| [28] | Schork 2022 | Crossref, Zotero | full text (Zotero PDF, author manuscript) |
| [29] | Liao et al. 2023 | PubMed | full text (PMC) |
| [30] | Selukar, Prince, May 2025 | PubMed | full text (PMC); supplement not read |
| [31] | Jiang, Arnold, Betensky 2025 | PubMed | full text (PMC) |
| [32] | Senarathne, Overstall, McGree 2020 | Crossref, Zotero | full text (arXiv v3) |
| [33] | Shrestha, Jain 2021 | PubMed, Zotero | **abstract only** |
| [34] | Giai et al. 2023 | PubMed, Crossref, Zotero | **abstract only** |
| [35] | Cruz, Melo, Mylona 2026 | Crossref, arXiv | full text (arXiv v3, not the published version) |
| [36] | Ferron et al. 2009 | PubMed | full text (publisher PDF) |
| [37] | Ferron, Farmer, Owens 2010 | PubMed | full text (publisher PDF) |
| [38] | Moeyaert et al. 2013 (BRM) | PubMed | full text (publisher PDF) |
| [39] | Petit-Bois et al. 2016 | PubMed | full text (publisher PDF) |
| [40] | Declercq et al. 2019 | PubMed | full text (publisher PDF) |
| [41] | Declercq et al. 2022 | PubMed | full text (ERIC copy) |
| [42] | Moeyaert, Yang, Xu 2022 | PubMed, Crossref, Zotero | PMC page through a summarizing reader only |
| [43] | Chen, Pustejovsky 2024 | PubMed | full text (accepted manuscript, PsyArXiv) |
| [44] | Bouwmeester, Jongerling 2020 | PubMed | full text (PMC) |
| [45] | Michiels et al. 2020 | PubMed | full text (author manuscript, KU Leuven) |
| [46] | Levin, Ferron, Kratochwill 2012 | PubMed, Crossref, Zotero | **abstract only** |
| [47] | Tueller et al. 2023 | PubMed | full text (PMC); table values through a summarizing reader |
| [48] | Natesan Batley, Hedges 2021 | PubMed | full text (PMC) |
| [49] | Joo, Ferron 2019 | PubMed | **abstract only** |
| [4] | Hendrickson et al. 2020 | PubMed | full text (local PDF) |

Table: Studies screened in full and how each was verified

The background items of A.2.10 were checked as follows: [3], [50],
[51], [53], [54], [56], [57] and [55] against PubMed and read in full
(for [53], the Zotero full-text cache); [52] against PubMed, abstract
and targeted passages of the PMC page; [58] and [59] against PubMed,
abstract only; [60] against Crossref only, not read (no abstract in
Zotero, publisher page not accessible); [61] against PubMed, not read.

### A.2 The studies

The studies are grouped by what they contribute. Within each, the
extraction terms of the task are given in a fixed order: design; DGP
(effect form, onset, carryover, serial correlation, heterogeneity,
expectancy and trend, missingness, outcome); estimand; analyses;
inference; performance measures and MCSE; scale; comparator;
conclusions. "nr" means not reported in what was read.

#### A.2.1 Closed-form variance-component frameworks

**Senn 2019 [2].** Closed form, no simulation. Design: $N$ participants
each randomized within $k \ge 2$ cycles of a pair of occasions. Model on
within-cycle differences: $d_{ij} = \beta_D + v_i + e_{ij}$,
$v_i \sim N(0, \tau^2)$ (Senn's treatment-by-patient interaction),
$e_{ij} \sim N(0, 2\sigma_w^2)$; all terms independent, no serial
correlation, trend or period effect; carryover is assumed removed by a
washout or by measuring late in each occasion (quoted: "carry-over can
be eliminated either by using a washout or by limiting measurement to
the latter part of any treatment occasion"). Estimands: $\beta_D$,
individual effects (naive and shrunk), and the variance components.
Analyses: the summary-measures random-effects analysis,
$\mathrm{Var}(\hat\beta_D) = (\tau^2 + 2\sigma_w^2/k)/N$ on $N - 1$
degrees of freedom; the fixed-effects analysis for the participants
studied, $2\sigma_w^2/(Nk)$ on $N(k - 1)$ degrees of freedom, with the
remark that software assuming $N - 1$ misstates power; shrinkage
standard errors that treat the weights as known. Worked example: with
$\tau^2 = 1$, $\sigma_w^2 = 4$, $\beta_D = 1$, $k = 3$, the random-effects
analysis needs $N = 31$ and the fixed-effects analysis $N = 22$ for 80%
power. Conclusion (quoted): "very different purposes for n-of-1 trials
may require very different samples sizes."

**Araujo, Julious and Senn 2016 [11].** Mainly theory, with a small
simulation. Model: $Y_{irs} = \lambda_i + \beta_{ir} + Z_{irs}(\beta_D +
v_i) + \varepsilon_{irs}$, $Z = \pm 1/2$, with a patient effect, a
cycle-within-patient effect, $v_i \sim N(0, \tau^2)$ (their $\psi^2$)
and independent error; no carryover, serial correlation, onset or
missingness in the simulation. The central closed-form result
(derived in the paper): the true variance of the mean effect is
$\tau^2/N + 2\sigma_w^2/(Nk)$, while the matched-pairs $t$-test over all
cycles estimates $\{(N-1)/(Nk-1)\}\tau^2/N + 2\sigma_w^2/(Nk)$, which is
too small whenever $\tau^2 > 0$ and $k > 1$; under the strict null
($\tau^2 = 0$) the matched-pairs test is valid. Simulation: $\tau^2$ from
0 to 1, $\sigma_w^2 = 1$, $k = 3$, 1,000 trials per value, comparing
empirical and model variances (figures only; no size, power or MCSE).
They recommend allowing for the random interaction when computing
intervals, and (quoted) "We do not recommend attempting to adjust for
carry-over but instead recommend using washout periods of adequate
length." They also argue that AR(1) errors alone do not allow for a
between-patient main effect and should be applied to within-patient
errors only. The paper is in part a critique of Chen and Chen [21].

**Yang, Steingrimsson and Schmid 2021 [16]** (preprint; no journal
version found). Closed form, no simulation. LMM
$Y_{ijkl} = \mu_{ij} + \delta_{ij}A_{ijkl} + \varepsilon_{ijkl}$ with
fixed or random intercepts crossed with common or random slopes
(random slope variance $\tau^2$, their $\sigma^2_\delta$), and
independent, exchangeable or AR(1) residuals, with several measurements
per period. Standard error by generalized least squares with known
variance components and normal-approximation power; a Shiny tool.
Assumed absent (quoted): "no underlying trend in the outcome or
carryover effects." Findings: random-slope models need more total
measurements than common-slope models, and the residual structure
changes the required sample size. The authors note that a $t$ rather
than normal reference would be needed in small samples.

**Hedges, Shadish and Natesan Batley 2023 [17]; Hedges, Pustejovsky and
Shadish 2012 [18].** Closed-form power ([17]) and the effect-size theory
behind it ([18]) for multi-case $(AB)^k$ reversal designs. Model:
$Y_{it} = \mu_{\text{phase}} + \eta_i + \varepsilon_{it}$ with a random
intercept $\eta_i$ of variance written $\tau^2$ in these papers (an
intercept variance, **not** treatment heterogeneity: the treatment
effect is the same for every case) and stationary AR(1) errors running
across phases. Effect: an immediate, constant mean shift; no onset,
carryover, trend or expectancy. The test statistic is the average
within-person contrast over a pooled between-case SD, referred to a
noncentral $t$ with a Box-approximation degrees of freedom that depends
on $\rho$ and on the intraclass correlation, both treated as known in
[17]. In [17], power is non-monotone in $\rho$, with a minimum near 0.6
to 0.8; with $k = 2$ pairs, 3 observations per phase, $\rho = 0.5$,
intraclass correlation 0.5 and a standardized effect of 0.75, 80%
power needs 7 cases. [18] simulated the $(AB)^1$ case only (8,000
iterations per cell) and found the variance of the estimator
underestimated when the nuisance parameters are estimated; it
evaluates estimation, not test size.

**Shi and Ye 2024 [19].** Closed form with simulation, AB/BA two-period
crossover against a parallel-group design. Potential-outcomes model in
which carryover is a fixed additive shift in period 2 ($\lambda_0$,
$\lambda_1$), with a period (natural-course) change and, in one case,
treatment-by-covariate interactions. Closed-form bias of the basic
crossover estimator, $\tfrac12(\theta_1 + \tilde\theta_2 - \lambda_0 -
\lambda_1)$, and closed-form size and power. Conclusion (quoted):
"when the carry-over effect is positive, the basic estimator
underestimates the treatment effect, which does not inflate the type I
error of one-sided tests but negatively impacts the power"; negative
carryover inflates size. With within-person correlation 0.3, 0.5 and
0.7, the crossover beats a parallel trial only if the mean carryover
is below 0.41, 0.50 and 0.61 times the alternative effect. Inference is
asymptotic $z$ with design-based covariate adjustment; 10,000
replicates per grid point in one case, 2,000 in a resampling study; no
MCSE. "Behavioral" carryover in this paper is a learned change in
behavior that persists, the sense in which the protocol labels its
$t_{1/2} = 1$ week scenario.

#### A.2.2 Aggregated N-of-1 trials against parallel and crossover trials

**Blackston et al. 2019 [12].** Design: aggregated N-of-1 with 3 cycles
of one treatment and one placebo observation (6 per participant);
the crossover comparator uses cycle 1 (2 observations), and the
parallel comparator one observation; all three are cut from the same
simulated data; $N$ from 10 to 200. DGP: $Y_{ick} = \mu_i + \beta_D
Z_{ick} + \varepsilon_{ick}$ with a random intercept $\mu_i$ of SD
$\sigma_\mu \in \{0, 0.1, 0.5\}$ and error SD $\sigma_\varepsilon \in
\{0.5, 1\}$, in four scenarios; $\beta_D = 0.25$ (written $\tau$ in the
paper). The "heterogeneity" of the scenarios is intercept variance;
the treatment effect is common to all participants except in a
selection-bias mixture in which a fraction of participants comes from
a subpopulation with a different effect. Carryover: a constant (0.05,
0.10 or 0.15) added to the next observation after any treatment
observation, regardless of time (quoted: "we do not consider length of
cycle or washout"), present also when $\beta_D = 0$. No serial
correlation, trend, onset or expectancy. Analyses: `lme4::lmer` with a
random intercept, treatment tested by a deviance (likelihood-ratio)
test; ordinary regression for the parallel trial; an extra model with
an indicator for treatment in the previous period. Inference
asymptotic; REML or ML not stated. 5,000 replicates per configuration;
MCSE nr. Results without carryover: the aggregated design reaches 80%
power at $N = 30$ in three scenarios and $N = 100$ in the fourth; the
parallel trial needs $N = 150$ and the crossover $N = 100$; size 0.05
for the aggregated design, 0.05 to 0.07 for the crossover. With
carryover: size of the aggregated design 6% to 8% at 0.05, 14% to 15%
at 0.10, 23% to 27% at 0.15; the crossover 6% to 17%. Under the
selection mixture with half the sample from the other subpopulation,
the aggregated design's size reached 37% to 40% in three scenarios.
Conclusion (quoted): "Carryover effects can result in inflated type I
error probabilities for both crossover and especially N-of-1 designs."
Whether the previous-treatment adjustment removed the inflation is
ambiguous in the text read.

**Carrozzo et al. 2025 [20].** Design: parallel trial, AB/BA
crossover, and an aggregated N-of-1 design of 3 cycles in which each
participant repeats one sequence; 10 to 100 per arm or sequence. DGP:
$Y_{ijk} = \gamma_i + \beta_D Z_{ijk} + \varepsilon_{ijk}$, a random
intercept $\gamma_i$ with variance written $\tau^2 \in \{0, 0.01,
0.25\}$ (intercept variance, cancelling in within-person differences),
and error SD drawn per participant from $U(0.05, 1)$. Carryover: 10% of
$\beta_D$ added to treatment-A observations that follow B; because it
is proportional to $\beta_D$ it vanishes under the null, so it cannot
affect size by construction (derived). No serial correlation, trend,
onset, expectancy or dropout. Analyses: unpaired $t$ (parallel),
period-difference $t$ (crossover), and a two-stage random-effects
meta-analysis of per-participant mean differences with a $z$ interval;
a sequential cumulative meta-analysis with a Bonferroni correction.
5,000 replicates (100 for the sequential analysis); MCSE nr. Results:
the N-of-1 meta-analysis reaches 80% power at about $N = 30$ (40 with
carryover), the crossover at 150, the parallel trial at 300 to 420;
size about 0.05 for all. The paper states that the Kenward-Roger
method was used as "the estimator of the between-study variance";
Kenward-Roger is an inference adjustment, not a variance estimator, so
the method used is unclear.

#### A.2.3 Comparisons of analyses for aggregated trials

**Chen and Chen 2014 [21].** (Two authors; the citation "Chen, Chen and
Xia" in the current paper 04 is wrong.) Design: aggregated N-of-1 with
3 or 4 AB cycles, one observation per period, $N \in \{1, 3, 5, 10, 20,
30\}$. DGP: multivariate normal with unit variances and compound
symmetry ($\rho_{cs}$ = 0, 0.5, 0.8) or AR(1) ($\rho = 0.5$, 0.8) across
periods; effects 2 for B and 2 to 3 for A (differences 0 to 1);
carryover from the previous period only, 0% or 20% of the previous
treatment's absolute effect, so placebo also carries over. No
$\tau^2$, onset, trend, expectancy or dropout. Analyses: M1, paired $t$
over all pairs ignoring the participant; M2, a mixed model of
within-cycle differences with a random participant effect; M3, a mixed
model of the raw responses with period and carryover-indicator terms
and a random intercept; M4, DerSimonian-Laird meta-analysis of
per-participant mean differences. The fitted residual covariance of M2
and M3 was set to the generating one. Degrees-of-freedom method nr
(SAS 9.2; the defaults were presumably used). 3,000 replicates; MCSE nr.
Results (3 cycles): M1 and M3 near nominal size under compound
symmetry; M2 severely conservative at small $N$ (0.004 to 0.006 at
$N = 3$); all methods conservative under AR(1) (M1 0.029, M3 0.015 to
0.021 at $N = 30$); M4 anticonservative at $N = 3$ (0.081 to 0.086) and
biased upward by 26% to 182% of the true difference even without
carryover, which suggests an implementation problem (argued). With 20%
carryover, M1 and M2 were biased by about $-10\%$ and M3 was nearly
unbiased. Conclusion (quoted): "We recommended paired t-test to be
used for normally distributed data of N-of-1 trials because of its
optimal statistical performance. In the presence of carryover effects,
mixed effects model could be used as an alternative." Araujo et al.
[11] show that this recommendation depends on $\tau^2 = 0$.

**Gärtner et al. 2023 [22].** Design: aggregated N-of-1 comparing two
active treatments, daily measurements, two blocks of two 4-week
periods; samples of 5 to 100 drawn from one generated population of
1,000. DGP: a directed acyclic graph of linear structural equations
with time-varying covariates, a Wiener-process drift (natural course),
constant effects ($-2$ and $-4$), and carryover through a recursive
wash-in and wash-out indicator with time constants per treatment (the
Percha form [27]); 10% missing rows with probability rising over time,
plus a 10-day gap per person; an ordinal outcome analyzed as truncated
Gaussian. No random effects or treatment heterogeneity (quoted: "we
only included linear dependencies and fixed effects"), no expectancy,
no explicit AR residual. Analyses: mean difference; OLS with
covariates; a carryover-adjusted parametric model (COAPM) whose decay
regressor's time constants are chosen by grid search; Bayesian
networks; G-estimation via GEE with independence or AR(1) working
correlation. The standard error reported is the SD of estimates over
100 subsamples, which are not independent replicates (argued); no
tests, no coverage, no MCSE. Conclusion (quoted): "When there is
carry-over, COAPM yields unbiased and more efficient estimates while
all other methods show some bias"; the authors concede that COAPM was
close to the generating process.

**Krone et al. 2020 [23].** Estimation only. Aggregated N-of-1 with
$N \in \{20, 30, 40\}$ and 10 to 30 measurements in random order;
$Y_{it} = \beta_{0i} + (\beta_D + v_i)x_{it} + \varepsilon_{it}$ with a
random intercept and random treatment effect (the text and a figure
disagree on whether 0.5 is the variance or the SD of $v_i$); no serial
correlation, carryover or trend. Two-stage meta-analysis against
frequentist and Bayesian random-slope LMMs, judged by the RMSE of
individual effects and the estimated heterogeneity; 1,000 replicates
per scenario; no test of $\beta_D$. Conclusion: the mixed models are
preferable, and the two-stage method overestimates heterogeneity.

#### A.2.4 Single-person series

These studies analyze one participant. They bear on the protocol
through the within-person part of the information and through the
small-sample behavior of serial-correlation estimates.

**Tang and Landes 2020 [24].** Single person, paired (multiple
crossover) and two-sample designs, $m$ from 4 to 12 and 30 to 100;
$Y \sim N(X\beta, \sigma^2 V)$ with $V_{ij} = \rho^{|i-j|}$, $\rho \in
\{-0.33, 0, 0.33, 0.67\}$; level or rate change; no carryover. Four
"serial $t$-tests" using a bias-corrected lag-1 estimate of $\rho$ and
an effective sample size, against the usual $t$-tests; one-sided 5%;
10,000 replicates per cell with a stated Monte Carlo error target
("within 0.01 ... with 95% confidence"). The serial tests improve on
the usual tests, but (quoted) "Type I error is still substantial,
power optimistic, and interval widths biased for small n"; they "do
not account for carryover effects". Paper 04's specific values for
this paper (size about 0.22 against 0.06 at $m = 12$, $\rho = 0.67$)
come from figures that were not read here and are unverified.

**Wang and Schork 2019 [25].** Single person, 400 observations in 2
to 80 alternating periods, fixed order, with or without washout; AR(1)
errors; a constant effect of 0.3 SD; quoted: "we assumed a constant
treatment response, not a graded one ... and no carry-over or
intervention order effects for many of the analyses." Analytic power
from the likelihood-ratio test under AR(1), checked by 1,000
simulations per setting fitting GLS with AR(1); asymptotic $\chi^2_1$.
Power falls with stronger serial correlation, which more frequent
alternation partly mitigates; washouts reduce power slightly when there
is no carryover. The current paper 04 attributes to this study a
finding that washouts reduce effective sample size "in proportion to
the carryover half-life"; that finding is not in the text read.

**Lin et al. 2016 [26].** Single person, 28 daily observations, ABAB or
AB; $Y_t = \alpha + \beta_D D_t + e_t$ with AR(1) $e_t$, $\rho \in \{0,
0.2, 0.5, 0.7\}$, normal or contaminated-normal errors. GLS with AR(1)
by REML, tested by a Wald test with a normal reference, a parametric
bootstrap and a residual bootstrap; 10,000 replicates per cell, with a
stated tolerance interval for size (0.0457 to 0.0543). Wald size 0.061
to 0.113; parametric bootstrap 0.050 to 0.054; REML bias of
$\hat\rho$ about $-0.12$ at $\rho = 0.7$.

**Percha et al. 2019 [27].** Single person, blocks of randomly ordered
periods. The most mechanistic onset-and-carryover DGP found: for each
treatment, the effect state approaches its asymptote with a wash-in
time constant while active and decays with a separate wash-out time
constant when stopped (first-order kinetics); a Wiener-process
baseline drift; mean reversion of the latent state toward drift plus
effects; process noise and measurement error; discrete outcomes by
link transforms. No heterogeneity, no expectancy component. Analysis:
OLS with treatment and block effects, ignoring serial correlation, and
a median-differencing rule; 50 trials per setting; no size, no MCSE.
Treatment order, period length and number of blocks changed both power
and the mean estimate (e.g. a true effect of 10 estimated at 5.8 to
12.0 across four orders). Recommendation (quoted): "our main
recommendation for N-of-1 study designers is to simulate the study".

**Schork 2022 [28].** Single person, 80 to 400 measurements in
alternating blocks with unmeasured washout gaps; AR(1) with $\rho$ from
$-0.99$ to 0.99; no carryover. OLS $t$, `nlme::gls` (whether an AR(1)
term was specified is unclear in the text), Newey-West HAC errors, and
sequential probability ratio tests (standard and bootstrap) for
single and aggregated series; 100 to 10,000 replicates. False
positives rise with positive $\rho$ and are worse with washout gaps;
randomized order mitigates this; HAC errors help.

**Liao et al. 2023 [29].** Single person, 120 daily observations with
several switching schedules and no washout. A distributed-lag outcome
model, $Y_t = \mu + \sum_{l=0}^{7}\beta_l X_{t-l} + \varepsilon_t$ with
AR($p$) errors, and five lag curves (geometric decay, slow absorption
peaking at lag 2, non-monotone curves with sign changes, none). A
Bayesian distributed-lag model with AR errors against three Bayesian
lag models with independent errors; bias and RMSE of immediate,
carryover and total effects; 100 replicates per scenario; no tests, no
MCSE. Under-specifying the AR order hurt more than over-specifying it.

#### A.2.5 Sequential and adaptive designs

**Selukar, Prince and May 2025 [30].** Single-person trials of 3 to 26
blocks of 2 to 6 periods, combined in series of 5 or 10 trials; normal
data with a block random intercept; quoted assumptions "no carryover
effect" and "no outcome drift"; one true effect for all trials of a
series, so no heterogeneity (inferred). An LMM with a block random
intercept, monitored by O'Brien-Fleming or Pocock boundaries, with
naive, bias-adjusted and median-unbiased estimators; inverse-variance
combination across trials. 10,000 trials per setting; MCSE nr. With
13 or fewer blocks, normal boundaries inflated size two to seven
times, while (quoted) "all settings with t-adjusted boundaries have
nominal or smaller than nominal type 1 error rates". Paper 04's
summary ("reaches nominal rates with more blocks or alpha-spending
corrections") should be revised to this.

**Jiang, Arnold and Betensky 2025 [31].** Single-person cycles of two
occasions with washout, monitored by O'Brien-Fleming boundaries on the
cumulative mean difference with a three-way decision (superiority,
equivalence, inferiority) and futility; individual effects from a
three-point mixture $\{-\mu, 0, \mu\}$ with no within-class spread; a
curtailed binomial rule across participants for the proportion of
responders. Known-variance normal theory; 1,000 replicates; MCSE nr.
Restricted by the authors to treatments with "rapid onsets and
cessations of action". The population target is a proportion of
individual classifications, not $\beta_D$.

**Senarathne, Overstall and McGree 2020 [32].** Aggregated N-of-1 with 3
cycles; GLMM with random intercept and random treatment effect
(variance 2.25 or 9, our $\tau^2$), normal, log-normal and Poisson
outcomes; no carryover or serial correlation. Compares Bayesian
information-optimal adaptive allocation, a bandit allocation and fixed
randomization by posterior precision and correct identification of each
participant's better treatment; 20 repetitions per scenario; no
frequentist size or power.

**Shrestha and Jain 2021 [33]** (abstract only). A Bayesian hierarchical
model for individual and population effects with Thompson-sampling
allocation in individual and aggregated N-of-1 trials, evaluated by a
simulation motivated by a completed trial. DGP, replicates and measures
could not be read.

**Giai et al. 2023 [34]** (abstract only). Individual and population
Net Benefit from generalized pairwise comparisons with prioritized
outcomes, in a simulation mimicking a 38-patient series; size was
reported as good and power as low even at 140 days per patient. A
nonparametric probabilistic estimand, not a mean difference.

#### A.2.6 Carryover in crossover designs with repeated measures

**Cruz, Melo and Mylona 2026 [35]** (read as arXiv v3). AB/BA crossover
with 10 to 100 measurements per period and 2 to 50 units per sequence;
a within-period trend and sequence-specific carryover functions
(sine and cosine) centered to mean zero; independent errors; no
heterogeneity. GEE with exchangeable working correlation and sandwich
errors, with and without penalized-spline carryover terms; 500
replicates; MCSE nr. All models estimated the treatment effect without
bias, but models without carryover terms had coverage near 70%, and
the proposed models near 95% as units or measurements increased.
Because the carryover functions have mean zero, they cannot bias the
treatment effect by construction (derived), so the study speaks to
variance and coverage, not to bias.

#### A.2.7 Single-case designs: multilevel Monte Carlo studies

All of these use multiple-baseline AB designs (none simulates ABAB
reversals), an immediate level shift (sometimes with a slope change),
no carryover, no expectancy and no dropout.

**Ferron et al. 2009 [36].** $N \in \{4, 6, 8\}$, series of 10 to 30;
$Y_{ij} = b_{0j} + b_{1j}\,\text{phase}_{ij} + r_{ij}$ with random
intercept and random treatment effect (variance 0.1 or 0.3, our
$\tau^2$) and AR(1) level-1 errors, $\rho$ from 0 to 0.4. REML LMM with
independent or AR(1) errors, crossed with five degrees-of-freedom
methods: containment, residual, between-within, Satterthwaite and
Kenward-Roger. 5,000 replicates per condition; MCSE given once
(0.0031 for coverage 0.95). With AR(1) modeled, Kenward-Roger coverage
averaged 0.952 (range 0.935 to 0.965) and Satterthwaite 0.949;
residual and between-within about 0.917; containment 0.972 with wider
intervals. Without AR(1), coverage fell as $\rho$ rose. The treatment-
effect variance was overestimated (relative bias 0.34, 0.25, 0.21 at
$N = 4$, 6, 8) and its interval covered about 0.83 to 0.85; the
authors call these results "not encouraging".

**Ferron, Farmer and Owens 2010 [37].** Same design and model;
individual effects $b_{1j}$ estimated by per-person OLS and by
empirical Bayes from an AR(1) LMM; 1,000 replicates. Empirical Bayes
estimates were shrunk (mean bias $-0.24$ on a true 1.0) but had lower
MSE. Interval coverage: Kenward-Roger 0.958, Satterthwaite 0.920,
containment 0.862, OLS 0.899. The authors judge Kenward-Roger
"preferable" for individual effects.

**Moeyaert et al. 2013 [38].** Three-level multiple-baseline synthesis
with independent errors and an "external event" that shifts all
participants of a study at the same occasions (constant or fading),
the only time-local shock common to everyone in this group, loosely a
period effect. Ignoring the event biased the effects (coverage of the
trend effect as low as 0.26); modeling it restored coverage of about
0.94 to 0.97. 400 replicates per condition; MCSE nr.

**Petit-Bois et al. 2016 [39].** Three-level synthesis; errors
independent, AR(1) ($\rho = -0.3$, 0.3, 0.6) or MA(1), always analyzed
with AR(1) and Kenward-Roger. Fixed effects nearly unbiased, coverage
0.913 to 0.981 under correct, over- and misspecified error structures;
variance components badly biased. 3,000 replicates; MCSE nr.

**Declercq et al. 2019 [40].** Count outcomes from a Poisson GLMM with
random intercept and random treatment effect; no autocorrelation. A
misspecified normal LMM with Satterthwaite degrees of freedom held
size at 0.04 to 0.05, while the correctly specified GLMM with a Wald
$z$ test was inflated (0.07 to 0.12). 2,000 replicates; MCSE given once.

**Declercq et al. 2022 [41].** One-stage three-level LMM against
two-stage meta-analysis of per-case effects, $z$ and $t$ inference; no
autocorrelation. Estimates nearly identical; coverage 0.951 to 0.964
for the one-stage model with Satterthwaite degrees of freedom against
0.912 to 0.917 for two-stage; size with 5 studies 0.03 to 0.05 against
0.09 to 0.10. 1,000 replicates; MCSE nr.

**Moeyaert, Yang and Xu 2022 [42]** (read through a summarizing reader).
Power to detect participant-level moderators of the treatment effect
with Kenward-Roger degrees of freedom and no autocorrelation; about 20
participants needed for three moderators. Conceptually closest to the
biomarker interaction of papers 01 and 02, not to $\beta_D$.

**Chen and Pustejovsky 2024 [43].** Two-stage meta-analysis of case-
level effect sizes from Poisson AR(1) data ($\rho$ 0, 0.2, 0.4);
model-based errors against robust variance estimation via
`clubSandwich` with a small-sample correction and Satterthwaite
degrees of freedom, the machinery of the protocol's CR2; 4,000
replicates. Robust errors protected coverage only when the case-level
estimator was unbiased (quoted: robust variance estimation "cannot
correct for systematic bias in effect size estimates").

#### A.2.8 Single-case designs: power, randomization tests and Bayesian analysis

**Bouwmeester and Jongerling 2020 [44].** Multiple-baseline AB, 2 to 12
participants, 15 to 60 measurements; AR(1) $\rho$ 0 to 0.5; **gradual
onset** of the effect over 0, 1/4, 1/3 or 1/2 of the intervention
phase; a one-tailed randomization test with 400 sampled permutations.
500 replicates per estimate, each scenario repeated 100 times, with the
SD of power reported. Size 0.046 to 0.05 with 3 or more participants.
Gradual onset cut power sharply: 0.97 with immediate onset, 0.40, 0.26
and 0.15 with onset over 1/4, 1/3 and 1/2 of the phase; $\rho = 0.5$
cut it to 0.80. The authors recommend excluding measurements taken
while the effect is emerging.

**Michiels et al. 2020 [45].** Randomized multiple-baseline, 3 to 6
cases; normal, uniform, AR(1) ($\rho = 0.6$) or bimodal data; a
two-level LMM with random intercept and random treatment effect and
Kenward-Roger $F$, against the same model wrapped in a randomization
test with 1,000 assignments. 1,000 replicates; MCSE nr. The wrapper
held the nominal level throughout; Kenward-Roger was grossly inflated
for bimodal data. The text read is internally inconsistent about the
model's size under the other distributions.

**Levin, Ferron and Kratochwill 2012 [46]** (abstract only).
Randomization and permutation tests for ABAB...AB and alternating-
treatment designs, systematic or randomized phase order, with or
without positive autocorrelation; systematic designs kept size at or
below nominal with appropriately conceived tests.

**Tueller et al. 2023 [47].** AB designs with 20 to 100 participants and
30 to 100 observations; phase jump or gradual linear change; random
intercept and slope, **independent** errors (quoted: "which reflects
no autocorrelation among residuals"); carryover not simulated. Power
only (one-sided), 1,000 replicates per condition with 95% intervals
for power; size not evaluated.

**Natesan Batley and Hedges 2021 [48].** Single-case AB series, phase
length 5 to 15, AR(1) $\rho$ 0 to 0.8; four Bayesian models; 100
datasets per condition. In short series, slope and AR(1) could not be
separated (both credible intervals contained 0 in 66% of datasets);
the authors favor simpler models.

**Joo and Ferron 2019 [49]** (abstract only). Non-normal multiple-
baseline data; within- and between-series estimates by REML and MCMC;
minimal bias, coverage "frequently acceptable, but varied".

#### A.2.9 The ancestor of the refined process

**Hendrickson et al. 2020 [4].** Designs OL, OL+BDC, CO and Hybrid
over 20 weeks with 2.5 weeks of titration at each start. DGP: the
three-component construction with modified Gompertz means and
variances, expectancy 1 open label and 0.5 blinded, exponential decay
of the last on-drug $BR$ after discontinuation, a joint multivariate
normal draw with within-factor autocorrelation, cross-factor
correlations and a biomarker-$BR$ correlation (covariance moderation),
coerced to positive definite; outcome-dependent dropout. Estimand: the
biomarker interaction; the paper states that its power "does not
address whether the treatment is effective as compared with placebo",
and **no main-effect results are reported**. Analysis: `lmer` with a
random intercept and no serial-correlation or carryover term (the text
also calls it an MMRM with unstructured covariance, which is
inconsistent); 1,000 replicates; MCSE nr. Even a 0.1-week half-life
gave "a precipitous decline in power" for the Hybrid and OL+BDC
designs.

#### A.2.10 Screened and kept as background

- **Applications with competing analyses, no simulation.** Zucker et
  al. 1997 [50] (hierarchical Bayesian aggregation and shrinkage,
  fibromyalgia); Zucker, Ruthazer and Schmid 2010 [3] (summary and IPD
  meta-analysis, LMMs with four residual structures chosen by BIC,
  Bayesian models, and design-mimicking subsets on 58 trials; the only
  mention of simulation is a call for it; a precision table for
  patients versus period pairs). The paper 04 text that says this
  study "characterized the within- and between-patient variance
  structure" is accurate; it should not be cited as a simulation.
- **Protocols, reviews and trial reports with a single simulation.**
  Stunnenberg et al. 2015 [51] (a Bayesian expected-posterior-
  probability calculation, 1,000 iterations); Stunnenberg et al. 2018
  [52] (a post hoc simulation showing power 69% and bias $-0.46$ on an
  assumed effect of 1.75 under the trial's response-based individual
  stopping rule); Nikles et al. 2011 [53] (one sentence on a SAS
  simulation for sample size under cycle-wise attrition).
- **Tutorial.** Schmid and Yang 2022 [54]: Bayesian models with trend,
  AR(1) terms, carryover decay regressors and missing-data imputation;
  no simulation of operating characteristics.
- **One simulated trajectory.** Defelippe et al. 2026 [55]: Bayesian
  continuous monitoring of a Poisson outcome with a placebo level
  shift; no replication.
- **Individualization, not the population mean.** Diaz 2021 [56]
  (partial empirical Bayes prediction for new patients, correlated
  random effects, no carryover or period effects); Cheung and Mitsumoto
  2022 [58] (abstract only: analytic power for evaluating a personalized
  trial program).
- **Carryover as target in a parallel trial.** Sturdevant and Lumley
  2021 [57]: linear post-discontinuation carryover in a prevention
  trial, MAR dropout by diagnosis.
- **Not read in full.** Moeyaert et al. 2013 [59] (abstract only: a
  three-level synthesis needing 30 or more studies and 20 or more
  occasions); Xie and Yu 2017 [60] (record only); Swan and Pustejovsky
  2018 [61] (a gradual-effects model for reversal designs, found late;
  record only, high priority to read, because its build-up and
  dissipation of effects resembles the refined process's onset and
  carryover).

### A.3 Comparison tables

Abbreviations: OL, open label; BDC, blinded discontinuation; CO, AB/BA
crossover; MBD, multiple-baseline design; RI, random intercept; RS,
random slope; CS, compound symmetry; KR, Kenward-Roger; Satt.,
Satterthwaite; LRT, likelihood-ratio test; MA, meta-analysis; HTE,
treatment heterogeneity ($\tau^2$); nr, not reported in what was read;
"intercept only" means the source's "heterogeneity" is a random
intercept, not $\tau^2$.

| Study | Design | Carryover in DGP | Serial correlation | $\tau^2$ (HTE) | Expectancy or trend |
|---|---|---|---|---|---|
| Senn [2] | aggregated, $k$ paired cycles | assumed absent (washout) | none | yes, normal | none |
| Araujo [11] | aggregated, 3 cycles | none | none | yes, 0 to 1 | cycle effect |
| Yang [16] | aggregated, periods with repeats | assumed absent | AR(1), CS | yes (RS) | none |
| Hedges [17, 18] | multi-case $(AB)^k$ | none | AR(1), known | no (intercept only) | none |
| Shi, Ye [19] | CO vs parallel | additive shift, signed | none | via covariate interaction | period change |
| Blackston [12] | 3 cycles vs CO vs parallel | constant added after treatment, also under null | none | no (intercept only); selection mixture | none |
| Carrozzo [20] | 3 cycles vs CO vs parallel | 10% of effect, zero under null | none | no (intercept only) | none |
| Chen, Chen [21] | 3 to 4 cycles, 1 obs per period | 0 or 20% of previous effect | CS $\rho_{cs}$ 0 to 0.8, AR(1) 0.5, 0.8 | no | none |
| Gärtner [22] | 2 blocks x 2 periods, daily | wash-in and wash-out recursion | drift only | no | Wiener drift |
| Krone [23] | aggregated, random order | none | none | yes | none |
| Tang, Landes [24] | single person | none | AR(1) $-0.33$ to 0.67 | n/a | rate-change slope |
| Wang, Schork [25] | single person, 400 obs | none (main runs) | AR(1) | n/a | none |
| Lin [26] | single person, 28 obs | none | AR(1) 0 to 0.7 | n/a | none |
| Percha [27] | single person, blocks | first-order wash-in and wash-out | drift plus mean reversion | n/a | Wiener drift |
| Schork [28] | single person, gaps | none | AR(1) $-0.99$ to 0.99 | n/a | slope (one study) |
| Liao [29] | single person, daily | distributed lag, 5 shapes | AR($p$), ARMA | n/a | trend, weekend |
| Selukar [30] | single person, blocks; series | assumed absent | block intercept | no | assumed absent |
| Jiang [31] | single person, cycles | assumed absent (washout) | none | three-point mixture | none |
| Senarathne [32] | aggregated, 3 cycles | none | none | yes | none |
| Shrestha [33] | individual and aggregated | nr (abstract) | nr | nr | nr |
| Giai [34] | series, 38 patients | nr (abstract) | nr | nr | nr |
| Cruz [35] | CO, repeated measures | mean-zero functions | none | no | within-period trend |
| Ferron [36, 37] | MBD, 4 to 8 cases | none | AR(1) 0 to 0.4 | yes, 0.1 to 0.5 | none |
| Moeyaert [38] | MBD, 3 levels | none | none | yes | external event |
| Petit-Bois [39] | MBD, 3 levels | none | AR(1), MA(1) | yes | slope change |
| Declercq [40] | MBD, counts | none | none | yes | none |
| Declercq [41] | AB, 3 levels | none | none | yes | linear, quadratic |
| Moeyaert [42] | MBD | none | none | yes, moderated | none |
| Chen, Pustejovsky [43] | AB series, counts | none | Poisson AR(1) | study and case level | none |
| Bouwmeester [44] | MBD, 2 to 12 | none; gradual onset | AR(1) 0 to 0.5 | no | none |
| Michiels [45] | randomized MBD | none | AR(1) 0.6 (one level) | per-occasion noise | none |
| Levin [46] | ABAB, alternating | nr (abstract) | 0 or positive | nr | nr |
| Tueller [47] | AB, 20 to 100 | none | none | yes (RS) | gradual slope |
| Natesan Batley [48] | single-case AB | none | AR(1) 0 to 0.8 | n/a | slopes |
| Joo, Ferron [49] | MBD, non-normal | nr (abstract) | nr | nr | nr |
| Hendrickson [4] | OL, OL+BDC, CO, Hybrid | exponential decay of $BR$ | within-factor autocorrelation | via $BR$ variance | expectancy 1 / 0.5; trend |
| **This protocol** | Hybrid, OL+BDC, CO, RCT (multi-crossover optional) | anchored exponential, $t_{1/2}$ 0 to 0.2 wk; behavioral 1 wk; Weibull | separable AR(1), $\rho = 0.7$ per week | yes, $\tau$ 0 to 5.49, heterogeneous nulls | $\eta_t$ 1 / 0.5; Gompertz $TV$ |

Table: Prior studies compared with this protocol: designs and
data-generating features

| Study | Analyses compared | Inference, df | Measures | MCSE | Reps per cell | Comparator |
|---|---|---|---|---|---|---|
| Senn [2] | summary-measures $t$; fixed-effects ANOVA; shrinkage | $N - 1$ vs $N(k-1)$ | analytic power, SE | n/a | n/a | none |
| Araujo [11] | paired $t$, mixed models, summary measures, MA, BLUP | df decompositions | empirical vs model variance | no | 1,000 | none |
| Yang [16] | LMM with RI/RS | GLS, known variances, $z$ | analytic power | n/a | n/a | none |
| Hedges [17, 18] | between-case standardized contrast | Box df, known $\rho$ | power ([17]); bias ([18]) | no | 8,000 ([18]) | none |
| Shi, Ye [19] | crossover and parallel $z$, covariate-adjusted | asymptotic | size, power | no | 2,000 to 10,000 | parallel |
| Blackston [12] | `lmer` RI, LRT; OLS | asymptotic LRT | size, power | no | 5,000 | parallel, CO |
| Carrozzo [20] | $t$ tests; two-stage MA | $z$ interval | size, power | no | 5,000 | parallel, CO |
| Chen, Chen [21] | paired $t$; 2 mixed models; DerSimonian-Laird | nr | size, power, bias, MSE | no | 3,000 | single person |
| Gärtner [22] | OLS, decay-regressor OLS, Bayesian networks, G-estimation | empirical SD | bias, SE | no | 100 subsamples | none |
| Krone [23] | two-stage, LMM, Bayesian LMM | none | RMSE of individual effects | no | 1,000 | none |
| Tang, Landes [24] | serial vs usual $t$ | effective df | size, power, interval width | target | 10,000 | none |
| Wang, Schork [25] | GLS AR(1) LRT | asymptotic | power, size | no | 1,000 | none |
| Lin [26] | GLS AR(1) Wald vs bootstrap | normal ref.; bootstrap | size, power | tolerance band | 10,000 | none |
| Percha [27] | OLS with blocks; median rule | OLS | power, estimate | no | 50 | none |
| Schork [28] | OLS, GLS, HAC, SPRT | $t$, HAC | false positives, power | partial | 100 to 10,000 | none |
| Liao [29] | Bayesian lag models | posterior | bias, RMSE | no | 100 | none |
| Selukar [30] | LMM with boundaries; combination | normal vs $t$ boundaries | size, power, bias, MSE | no | 10,000 | fixed design |
| Jiang [31] | group-sequential mean difference | known variance | decision rates, sample size | no | 1,000 | fixed design |
| Senarathne [32] | adaptive vs fixed allocation, Bayesian GLMM | posterior | precision, correct choice | no | 20 | fixed design |
| Shrestha [33] | Bayesian hierarchical bandit | posterior | nr (abstract) | nr | nr | nr |
| Giai [34] | generalized pairwise comparisons | nr | size, power (abstract) | nr | nr | none |
| Cruz [35] | GEE with and without carryover splines | sandwich, no correction | bias, coverage | no | 500 | none |
| Ferron [36] | LMM, independent vs AR(1) | containment, residual, between-within, Satt., KR | coverage, width; $\tau^2$ bias | single value | 5,000 | none |
| Ferron [37] | per-person OLS vs empirical Bayes | containment, Satt., KR | individual bias, MSE, coverage | single value | 1,000 | none |
| Moeyaert [38] | two-stage three-level, with and without event terms | Satt. | bias, coverage | no | 400 | none |
| Petit-Bois [39] | three-level LMM with AR(1) | KR | bias, coverage | no | 3,000 | none |
| Declercq [40] | LMM vs Poisson GLMM | Satt. vs Wald $z$ | size, power, bias | single value | 2,000 | none |
| Declercq [41] | one-stage vs two-stage | Satt. $t$; $z$ | size, coverage, bias | no | 1,000 | none |
| Moeyaert [42] | HLM with moderators | KR | power | no | 1,000 | none |
| Chen, Pustejovsky [43] | multilevel MA, robust variance | CR2-type, Satt. | bias, RMSE, coverage | no | 4,000 | none |
| Bouwmeester [44] | randomization test | permutation | size, power | SD of power | 500 x 100 | none |
| Michiels [45] | LMM KR vs randomization-test wrapper | KR; permutation | size, power | no | 1,000 | none |
| Levin [46] | randomization and permutation tests | permutation | size, power (abstract) | nr | nr | none |
| Tueller [47] | LMM RI+RS, piecewise | nr, one-sided | power only | 95% intervals | 1,000 | none |
| Natesan Batley [48] | four Bayesian models | posterior | RMSE, bias, coverage | no | 100 | none |
| Joo, Ferron [49] | REML vs MCMC, two estimators | nr | bias, coverage (abstract) | nr | nr | none |
| Hendrickson [4] | `lmer` RI (interaction only) | nr | power, bias | no | 1,000 | among four designs |
| **This protocol** | 24 analyses (Section 5.2) incl. RS, summary measures, oracles, randomization and max-$\lvert t\rvert$ tests (revised 2026-10-09 from docs/58, Section 12.3) | CR2 Satt.; KR; containment; exact $t$; HC2 | size, size-adjusted power, bias, empirical SE, $\kappa$, coverage, $\tau_d^2$, BLUP RMSE | every measure | 1,000 to 2,000 null; 500 to 1,000 alt. | parallel RCT (ANCOVA, MMRM) |

Table: Prior studies compared with this protocol: analyses, inference
and performance measures

### A.4 Evaluation

**What the literature does well** (inspected unless marked).

- **The variance algebra is settled.** Senn [2] and Araujo et al. [11]
  give the variance of the mean effect with $\tau^2$, separate the
  strict null (no effect in anyone) from the population null (mean
  zero), and show that an analysis that ignores $\tau^2$ is valid for
  the first and anticonservative for the second. This is the single
  most useful result for the protocol, because its heterogeneous nulls
  (block 2) are exactly the population null.
- **Small-sample inference in LMMs with AR(1) errors has been studied
  carefully,** but in multiple-baseline designs. Ferron et al. [36]
  compared five degrees-of-freedom methods with AR(1) modeled or not;
  Satterthwaite and Kenward-Roger gave coverage near 0.95 only when
  AR(1) was modeled, and containment was conservative. Declercq et al.
  [41] and Selukar et al. [30] show that the choice between normal and
  $t$ references changes size materially at small $N$ or few blocks.
- **Randomization tests** give exact size under the design's own
  randomization and are robust to non-normal data (Michiels et al.
  [45], Bouwmeester and Jongerling [44], Levin et al. [46]).
- **Onset and carryover have mechanistic models** outside the
  aggregated-trial literature: first-order wash-in and wash-out
  (Percha et al. [27], adopted by Gärtner et al. [22]), distributed
  lags with non-monotone shapes (Liao et al. [29]), gradual onset with
  a quantified power cost (Bouwmeester and Jongerling [44]).
- **Individual effects and their intervals** are evaluated (Ferron et
  al. [37], Krone et al. [23], Senn [2]); Bayesian hierarchical and
  adaptive designs exist ([32], [33], [50]); count and ordinal
  outcomes are covered in the single-case literature ([40], [43]).

**Where a DGP or analysis would mislead in a Hybrid, titrated-drug
setting** (argued, from the features in Table "Prior studies compared
with this protocol: designs and data-generating features").

- **Constant, immediate effect.** Every aggregated-trial simulation
  found uses a constant effect from the first on-drug observation.
  Under re-titration a one-week exposure reaches 3.1% of the ceiling
  (Section 3.3), so a constant-effect DGP overstates the information in
  short re-exposures, and its power comparisons with a parallel trial
  favor the N-of-1 design more than a titrated drug would. Bouwmeester
  and Jongerling's onset result [44] (power from 0.97 to 0.40 when the
  effect builds over a quarter of the phase) shows the size of the
  effect in a related design.
- **Carryover that is constant, mean-zero or zero under the null.**
  Blackston et al. [12] add carryover even when there is no effect, so
  their size inflation (up to 27%) measures a carryover with no drug
  behind it; Carrozzo et al. [20] make carryover proportional to the
  effect, so it cannot affect size; Cruz et al. [35] center it to mean
  zero, so it cannot cause bias. In the refined process carryover is
  the decayed drug response, zero under the homogeneous null and
  nonzero for individuals under the heterogeneous null (derived from
  Section 3.2). None of these three choices reproduces that.
- **Carryover excluded with oracle knowledge, or ignored.** Chen and
  Chen [21] adjust with the true carryover indicators; the v1 paper 04
  driver drops the rows the DGP flags (Section 1, problem 4). Both
  remove the problem the design must solve.
- **No expectancy.** No study other than Hendrickson et al. [4] has an
  expectancy component that differs between open-label and blinded
  phases; only Defelippe et al. [55] has a placebo level shift. The
  drug-expectancy confounding of the Hybrid and OL+BDC designs
  (Glossary) is therefore untested in the main-effect literature.
- **No treatment heterogeneity, under a name that suggests it.**
  Blackston et al. [12], Carrozzo et al. [20] and Hedges et al. [17,
  18] vary an intercept variance that cancels in within-person
  contrasts and call it heterogeneity or write it $\tau^2$. Their
  comparisons with parallel trials therefore show the gain from
  removing between-person baseline variance, not the cost of $\tau^2$,
  which in Senn's variance formula [2] becomes the dominant term as $k$
  grows (derived).
- **Comparisons that are not information-matched.** Blackston et al.
  [12] and Carrozzo et al. [20] compare 6 observations per person with
  2 (crossover) and 1 (parallel). The power advantage they report
  combines the design's within-person contrast with three to six times
  as many measurements.
- **Inference.** Asymptotic likelihood-ratio or Wald tests with normal
  references (Blackston et al. [12], Wang and Schork [25], Lin et al.
  [26], Carrozzo et al. [20]), unstated degrees of freedom (Chen and
  Chen [21], Tueller et al. [47]), and sandwich errors without
  small-sample correction at 2 to 10 clusters (Cruz et al. [35]).
  Lin et al. [26] show Wald size up to 0.113 at 28 observations; the
  protocol's CR2 with Satterthwaite degrees of freedom is the kind of
  correction this literature lacks.
- **MCSE.** Of 38 studies, 8 report a Monte Carlo error in some form
  ([24], [26], [28], [36], [37], [40], [44], [47]), most as a single
  value or a band for one measure; none reports it for every measure.
  Several have 20 to 100 replicates per cell ([22], [27], [29], [32],
  [48]), too few to resolve size differences of 0.01 to 0.02.

**What is new in this protocol relative to the literature found**
(argued; "new" means not found in the 38 studies, not a claim about
all published work).

- A main-effect simulation of the Hybrid and OL+BDC designs. The only
  study of these designs, Hendrickson et al. [4], targets the
  interaction and reports no main-effect result.
- An expectancy component that changes between open-label and blinded
  phases, and an analysis that models it (F1), so that drug-expectancy
  confounding appears as bias rather than being built into the target.
- Gompertz onset with restart at re-exposure, combined with anchored
  exponential carryover and a design-specific estimand
  ($\beta_D^{\text{exp}}$, Section 4) whose value is computed per cell.
- Heterogeneous nulls with $\tau$ set in absolute units, and a
  heterogeneity that follows the drug-response profile (zero before
  exposure, decaying after it).
- An RCT comparator generated from the same component model with a
  baseline measured with error, and an all-visits MMRM comparator (AM).
- Carryover handled without oracle information (XS, $D_{bc}$ at a
  misspecified half-life) and the oracle kept only as a labeled bridge.
- Closed-form probability limits and sandwich variances before
  simulation, with the simulation validating them (Section 5.1).
- MCSE targets and MCSEs for every measure, size-adjusted power and a
  working decision rule with an estimation screen.

**What the literature handles better, or the protocol lacks** (argued).

- **Randomization tests.** The protocol has none, although its paths
  are randomized; a permutation test is the one test whose size does
  not depend on the covariance model (Michiels et al. [45]).
  *(Revised 2026-10-09 from docs/58, Section 12.3:* RV and RF are now
  in Section 5.2.)
- **Bayesian hierarchical analysis** of $\tau^2$ and individual effects
  ([50], [32], [33], [54]) is absent; with 35 to 150 participants and
  a weakly identified slope variance, a Bayesian model with a proper
  prior on $\tau$ may estimate $\tau^2$ better than REML at the
  boundary. Ferron et al. [36] and Petit-Bois et al. [39] show how
  poorly REML estimates treatment-effect variances with few cases.
- **Individual-effect inference.** The protocol reports BLUP RMSE and
  correlation (Section 6.1) but not interval coverage for individual
  effects, and has no per-participant OLS comparator; Ferron et al.
  [37] show the degrees-of-freedom choice matters more for these
  intervals than for the mean.
- **Design adaptivity and sequential monitoring** ([30], [31], [32],
  [33]) are out of scope here; the bias from response-based individual
  stopping in a real aggregated trial (Stunnenberg et al. [52]) is a
  reason to keep stopping rules out of the rebuilt study's DGP unless
  they are studied deliberately.
- **Count, ordinal and bounded outcomes** ([40], [43], [22], [27]):
  the protocol's outcome is a normal total score with no range check
  (`docs/51`, Feature 17).
- **Non-monotone or slow-absorption carryover** (Liao et al. [29]) and
  **sign-dependent carryover** (Shi and Ye [19]: positive carryover
  attenuates and is conservative, negative carryover inflates size).
  The refined process has only positive, monotone carryover; a rebound
  after discontinuation is not represented.
- **Common time shocks** (Moeyaert et al. [38]) and **non-representative
  samples** (Blackston et al. [12], selection mixture): neither is in
  the protocol; the first is partly covered by the visit-factor
  analyses, the second is outside the study's aims.

**Corrections to the current paper 04's literature summary** (inspected
against the sources; for Section 9's rewrite):

- Chen and Chen [21] has two authors; "Chen, Chen and Xia" is wrong
  (the bib key problem of Section 9 is separate).
- Blackston et al. [12] and Carrozzo et al. [20] vary intercept
  variance, not treatment heterogeneity; "patient heterogeneity" in
  paper 04 should say so.
- The washout finding attributed to Wang and Schork [25] is not in the
  text read.
- Selukar et al. [30]: nominal size is reached with $t$-adjusted
  boundaries, not "with more blocks or alpha-spending corrections".
- Tang and Landes [24]: the size values quoted in paper 04 come from
  figures not read here and must be checked before reuse.

### A.5 Implications for the protocol

Each item is a proposal for the user's decision. None is applied to
Sections 1 to 7; all are listed in Section 8 as items 12 to 20.
*(Revised 2026-10-09 from docs/58, Section 12.3:* L1 is superseded by
item 12 as revised, which is entered in Section 5.2.)

- **L1. Add a randomization test** (Section 8, item 12). Re-draw the
  path labels under the protocol's allocation rule (uniform, at least 2
  per path) and recompute the statistic. For PR, a closed-form linear
  statistic, 999 permutations per trial cost little. For F1, which
  needs a refit per permutation, run 199 permutations in the cells of
  blocks 1 and 2 only. Reason: exact size under the strict null
  ([45], [44]); a reference against which CR2 and KR size can be read.
  *(Superseded on 2026-10-09 by item 12 as revised from docs/58,
  Section 12.3: null-model weights make the F1 version cheap enough to
  run in every cell without refits.)* Caveat (derived): a permutation test of path labels tests the strict
  null of no effect in anyone, not the population null of the
  heterogeneous cells, so its size there is itself a result.
- **L2. Add a Bayesian random-slope analysis, or record its omission**
  (item 13). Normal priors on fixed effects, a half-normal or
  half-$t$ prior on $\tau$, evaluated by the same frequentist measures in
  blocks 1 to 3 at $N = 35$ and 70. Cost: MCMC per trial is likely
  ten or more times an `nlme` fit (assumed); a pilot would decide.
- **L3. Add an onset-exclusion analysis XT** (item 14): F1 with the
  first on-drug visit after each start or restart dropped by schedule,
  the onset counterpart of XS. Reason: Bouwmeester and Jongerling [44]
  recommend excluding observations while the effect is emerging, and
  under re-titration the first on-drug visits carry little signal.
- **L4. Add a literature bridge cell** (item 15): constant effect, no
  $PB$ or $TV$, $t_{1/2} = 0$, AR(1) only, Hybrid and CO. In it PW and
  RS can be checked against the closed forms of Senn [2], Araujo et
  al. [11] and Yang et al. [16], and the protocol's code against
  published results. Reason: none of the published DGPs can be
  reproduced inside the refined process otherwise.
- **L5. Cross the heterogeneous null with carryover** (item 16): block
  2's null cells at $t_{1/2} = 1$ week. Reason (derived): only there
  does carryover exist under a null of the refined process, so only
  there can the Blackston-type size inflation [12] arise in a form
  consistent with pharmacology.
- **L6. Add a two-point heterogeneity cell** (item 17): $v_i$ taking
  two values with the variance of the $\tau = 5.49$ cell, as a
  responder and non-responder population (cf. the mixtures of [12] and
  [31]). It extends item 9.
- **L7. Designate the RCT comparator for design comparisons** (item
  18): AM, which uses every visit, as the information-matched
  comparator, and AN as the conventional one, so that the design
  comparison does not repeat the unmatched comparisons of [12] and
  [20].
- **L8. Add measures** (item 19): interval width for $\beta_D$ ([36]);
  coverage of individual-effect intervals with Satterthwaite and KR
  degrees of freedom and a per-participant OLS comparator ([37]).
- **L9. Add a bounded-outcome sensitivity cell** (item 20): $Y$ rounded
  and truncated to the CAPS range in the base alternative cell, with
  the fraction of out-of-range scores reported (`docs/51`, Feature 17;
  [40], [43]).

Proposals considered and not made: sequential monitoring and adaptive
allocation ([30] to [33]), because they change the design under study;
non-monotone carryover ([29]) and negative (rebound) carryover ([19]),
because the refined process has no mechanism for them and adding one
is a DGP decision beyond this review; common time shocks ([38]),
because the visit-factor analyses already absorb any shock common to
all paths at a visit.
