# The analysis model: exposure codings, time and baseline

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** E1 to E9 exposure codings; time adjustment; baseline as row or covariate; cLDA.  
**Prerequisites:** ch15

---

## Contents

1. [Notation](#notation)
2. [Summary](#summary)
3. [Construct box](#construct-box)
4. [Setting and argument](#setting-and-argument)
5. [Derivations](#derivations)
6. [Numerical evidence](#numerical-evidence)
7. [Relation to the manuscripts](#relation-to-the-manuscripts)
8. [Limitations and open items](#limitations-and-open-items)
9. [Reproduce](#reproduce)
10. [References](#references)
11. [Change log](#change-log)

---

## Notation

*Local symbols used in this chapter only. For project-wide notation see
`analysis/report/NOTATION.md`.*

<!-- notation entries -->

---

## Summary

<!-- evidence-labeled results go here -->

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N, analysis ID,
t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: CM "The analysis model and its carryover specifications" (6
  subsections); 43 Section 2; 40 Sections 8 to 9; 41 Propositions 1 to 5 and
  Section 7; 42 (pros, cons, update); 45 Section 7; 44 Section 5; 37 Appendix
  A "Visits before first exposure"; 24 "Why the analysis model does not
  match"; `plan-coherent-analysis-comparison.md` Q1 to Q3.

- New to write: one consolidated time-adjustment section (40 Section 9 and 42
  merged).
-->

### The common analysis model

<!-- From analysis/report/compendium-mathematics.tex Section "The
analysis model and its carryover specifications" -->

Across the compendium both architectures are fit identically:

```
Sx ~ bm + t + Dbc + bm:Dbc
```

via `nlme::lme` with `corCAR1(form = ~ t|ptID)` residual correlation
and a participant random intercept. The inferential machinery does not
know which architecture generated the data it is fitting.

The published analysis (commit `58b32a9`) fits the `lmer` model

$$
Y_{it} = \beta_0 + \beta_B B_i + \beta_D D_{it} + \beta_T t
  + \beta_{bm:D} B_i D_{it} + u_i + \varepsilon_{it},
$$

with the binary $D_{it}$ and baseline as a row with $D_{i0} = 0$.
Under carryover a binary indicator estimates an attenuated
on-minus-off contrast, so part of the power loss reflects the
analysis, not the data-generating process.

### Dimensions of an analysis

<!-- From docs/43 Section 2 -->

An analysis is a combination of about a dozen choices. They are: the
estimator family, the exposure coding, the time adjustment, the
expectancy term, the baseline handling, the random effects, the working
covariance, the standard error and degrees of freedom, the path term,
the biomarker coding, the target coefficient, and the design-specific
branches.

**Estimator family.** The main options are: linear mixed model with
`lmer` or `nlme::lme`; marginal MMRM with `mmrm` or `gls`; summary
measures (E9); strict RM-ANOVA with a dichotomized biomarker; GEE.
About 60 distinct approaches have been used, evaluated or proposed.

**Drug-exposure coding.**

| Option | Code | Source |
|---|---|---|
| Binary on-drug indicator `Db` | E1 | published default |
| Exposure-decayed `Dbc`, assumed half-life | E2 | package default |
| Binary plus lagged just-off-drug term $L_{it}$ | E3 | Jones-Kenward |
| `Dbc` with half-life chosen by AIC | E7 | liberal; fails size |

E4 to E6 are retired; E7 is evaluated as liberal in the crossover
(pooled 0.065 for F1A; `docs/43`, Section 2.2).

**Time adjustment.**

| Option | Source | Status |
|---|---|---|
| Linear `t` | published | evaluated; recommended |
| Natural spline in `t`, 3 df | `docs/40`, `docs/42` | excess size in Hybrid under KR |
| Visit as factor | `docs/40`, `docs/42` | as the spline under KR |
| Linear `t` plus expectancy `De` | orig option; F1 | evaluated; selected |

**Expectancy.** The design expectancy `De` (0 at baseline, 1 open
label, 0.5 blinded) is present in the published code but switched off
for the publication. Over the nine rows at the stress-test allocation,
the allocation-weighted correlation of $\eta_t$ with the binary
$D_{it}$ is 0.588 (`docs/56`, Section 3.3; verified). Including it in
the model (F1) is required when the biomarker may be prognostic
through the placebo response.

**Baseline handling.** Options range from the published baseline row
(coded off drug) to ANCOVA (post-baseline visits only, baseline as a
fixed covariate) to cLDA (baseline in the response with a common
mean). Their properties are given in `docs/41` Propositions 1 to 5
and summarized in the section below.

### Building the analysis model, term by term

<!-- From docs/40 Section 8 -->

The simulation's mean structure, for a participant with standardized
biomarker $b_i$, is

$$
\mathbb E(Y_{it} \mid b_i) = \mu_{\mathrm{BL}} - TV(w_t)
  - \eta_t\,PB(t_{pb}) - BR_q(t)
  - [\text{moderation}]_t \, b_i .
$$

The published model replaces each part with something estimable:

- **$\beta_0 + u_i$: baseline level and stable individual
  differences.** The random intercept matches compound symmetry. It
  ignores the change in variance at blinding.
- **$\beta_B B_i$: a prognostic biomarker effect,** zero in the
  simulation. The biomarker is uncentered (mean 124), so $\beta_D$ is
  the drug effect at a biomarker of zero. Centering fixes the
  interpretation and leaves $\beta_{bm:D}$ unchanged.
- **$\beta_D D_{it}$: the drug response $BR_q(t)$.** One number for a
  response that is 4.3 at week 4, 10.2 at week 10 and 4.3 again after
  re-exposure.
- **$\beta_T w_t$: the non-drug trends $TV(w_t) + \eta_t
  PB(t_{pb})$.** A straight line through two sigmoids and a step.
- **$\beta_{bm:D} B_i D_{it}$: the moderation,** the estimand.
- **$\beta_\eta \eta_t$ (`De`): the expectancy step.** Present in the
  code but switched off for the publication because of collinearity
  with $D_{it}$ (correlation 0.59 in the Hybrid long data; verified).

**Time as a factor.** Replacing $\beta_T w_t$ by one free mean per
visit absorbs everything common to all participants at a visit,
including $TV$, $PB$, the expectancy step, and the part of $BR$
every path shares. The drug effect is then identified only by the
between-path comparisons. A spline in time (`ns(t, 3)`) is the
intermediate choice.

**Open-label branch.** When no participant varies in drug state (the
OL design), the code fits `Sx ~ bm + t + bm:t + (1|ptID)` and tests
the biomarker-by-time slope. That slope is fully confounded with
natural course and expectancy.

**Summary measures and RM-ANOVA.** The E9 statistic is the
two-stage summary-measures counterpart: an on-minus-off contrast per
participant, then a regression on the biomarker. It equals the
interaction test of a repeated-measures ANCOVA with drug state as a
two-level within-subject factor. It loses power in the crossover,
where the time trend enters the two sequences' contrasts with opposite
sign (`docs/37`, Section 5.7).

### Does the form of the time adjustment matter?

<!-- From docs/40 Section 9 -->

Three forms of time adjustment were compared on the same simulated
trials: linear `t` (published), `ns(t, df = 3)`, and `factor(t)`.
The setting was the published construct (CS), $N = 70$, strength 0.25,
the published allocation, and three coupling rules at
$t_{1/2} \in \{0, 1, 2\}$ weeks, in Hybrid and CO (500 replicates per
cell, 1000 per null cell; script `07-time-adjustment.R`; verified).

**The interaction does not depend on the time adjustment.** In every
cell the three mean estimates of $\beta_{bm:D}$ agree to within 0.0009.
The biomarker is unrelated to the non-drug components, so their misfit
is irrelevant to the moderation.

**The drug main effect depends on the time adjustment.** Linear time
credits the open-label expectancy to the drug, overstating the effect
by about 2 points. Visit effects target the between-path contrasts
and pay in precision (SD 1.27 against 0.86).

**In Hybrid, flexible adjustments are anticonservative under the
published random-intercept covariance.** Under one-week carryover the
Type I error is 0.075 (spline) and 0.080 (visit factor), against 0.056
for linear time. Their higher power, 0.02 to 0.07, is partly this
inflation. In the crossover all three hold their size and the flexible
forms gain 0.01 to 0.04.

The excess under spline and visit factor is confined to model-based
and Kenward-Roger inference and is absent under CR2 (`docs/45`,
Section 4.4; two independent runs, each 2 to 4 MCSE; verified). Linear
time is never worse in any cell and is the recommended form.

<!-- TODO: one consolidated time-adjustment section merging docs/40
Section 9 and docs/42 into a single narrative with the key table, the
mechanism, and the CR2 verdict -->

### Baseline as a row or a covariate: five propositions

<!-- From docs/41 Sections 3 to 7 -->

The published analysis treats the baseline symptom score as a ninth
response row coded off drug. The alternative is baseline as a
covariate in a model for the post-baseline visits (ANCOVA). A third
option, constrained longitudinal data analysis (cLDA), keeps baseline
in the response but models it with a common mean across groups.

**Proposition 1: the covariate model is the conditional model.**
Conditional on baseline,

$$
\mathbf{Y}_{ip} \mid Y_{i0} \;\sim\; N\!\Bigl(
  X_{ip}\boldsymbol\beta + \boldsymbol\beta_{BL}\,(Y_{i0}
  - \mu_0 - \beta_B B_i),\;\; \Sigma_{pp\cdot 0} \Bigr),
\quad
\boldsymbol\beta_{BL} = \frac{\boldsymbol\sigma_{p0}}{\sigma_{00}}.
$$

ANCOVA with common baseline slope $\beta_{BL}$ is the special case
$\beta_{BL,t} \equiv \beta_{BL}$. Conditioning leaves $\beta_D$
and $\beta_{bm:D}$ unchanged.

In the construct, $\sigma_{t0} = \sigma_{00}$ for every $t$, so
$\boldsymbol\beta_{BL} = \mathbf{1}$ (verified) and the covariate
analysis coincides with change-from-baseline.

**Proposition 2: a covariate adds nothing to a within-participant
contrast when the slopes are equal.** With contrast weights
$\mathbf{a}$ with $\sum_t a_t = 0$, the reduction from conditioning on
baseline is $Q = \sigma_{00}(\mathbf{a}'\boldsymbol\beta_{BL})^2$, which
is zero when slopes are equal. In the construct
$\mathbf{a}'\boldsymbol\sigma_{p0} = 0$ exactly, and the contrast
variance is the same with and without conditioning (verified).

**Proposition 3: between participants, the classic ANCOVA gain.**
The conditional variance $\sigma_{tt} - \sigma_{t0}^2/\sigma_{00}$
is never larger than either the raw or the change-score variance.
In the construct, this is a 57% reduction per blinded visit (verified).
The gain accrues to between-participant terms (the biomarker main
effect and the path means), which are not the target.

**Proposition 4: the baseline row as an off-drug visit.** The
response coding includes baseline as an off-drug observation. With
graded coupling and no carryover ($\bar g_{\text{off}} = 0$), E9 and
E9b coincide. Under the published step rule with any carryover
($\bar g_{\text{on}} = \bar g_{\text{off}} = 1$), the E9 slope is zero;
E9b retains a residual power floor produced entirely by the uncoupled
baseline row. The contrast variance rises from 44.13 to 51.04 (paths
1 and 2) and from 45.03 to 51.77 (paths 3 and 4): increases of 16%
and 15% (verified). In the closed form this costs 0.06 to 0.07 in power.

**Proposition 5: the working covariance and the validity of the
standard errors.** The published random intercept implies
$V = \sigma_u^2 J + \sigma_\varepsilon^2 I$ over the nine rows.
Every row of $V$ carries the same residual variance $\sigma_\varepsilon^2 > 0$,
whereas the baseline row of $\Sigma$ has none. No choice of
$(\sigma_u^2, \sigma_\varepsilon^2)$ matches it, so the model-based standard
errors are not guaranteed correct. Only an unstructured covariance
or a sandwich estimator removes all of the mismatch.

**cLDA (Section 7 of docs/41).** With complete data, an unstructured
covariance and $\beta_B$ freed at baseline, keeping baseline in the
response is equivalent to conditioning on it. The follow-up study
(`docs/45`, Section 6) confirmed this equivalence. Its table in
`docs/41` applied the equivalence to the published row coding, which
shares $\beta_B$: that application is wrong (see Section 6 of `docs/45`).

**Status (2026-10-07).** The stress test selected F1 (a random
intercept with `corCAR1` residuals and CR2 inference, the baseline
row, linear time, and the expectancy terms) over the recommendation of
`docs/41`. The baseline row with a shared $\beta_B$ confers a power
advantage under AR(1) data when the biomarker is purely predictive,
because the baseline row adds a clean off-drug observation. The
analysis closest to the `docs/41` recommendation, EA, held its size
(pooled 0.052) but had the largest maximin shortfall of the survivors
(0.318; `docs/46`, Amendment 2).

### Carryover in the analysis: $D_{bc,it}$

<!-- From docs/45 Section 7 -->

The analysis may model carryover by using the exposure-decayed
indicator $D_{bc,it}$ (E2) in place of the binary $D_{it}$ (E1).

**Size.** $D_{bc,it}$ does not change the size of any analysis
(`RI+CAR1-CR2` is 0.051 to 0.054 throughout; verified).

**Power** (`us-KR`, linear time, raw): at $t_{1/2}^{a} = 0.5$ it costs
0.016 or less without carryover and gains 0.04 to 0.06 at
$t_{1/2} = 1$. In CO it makes little difference at any assumed
half-life.

**Visits before first exposure.** Pre-exposure off-drug visits must
be coded 0, not $e^{-\lambda \cdot 0} = 1$. Manuscript 02's
implementation codes them 1, giving CO's placebo-first path
$D_{bc,it} = 1$ at every visit and costing 0.26 to 0.30 in power
(`docs/45`, Section 8). Fix: `Dbc = ifelse(Db == 1, 1, ifelse(tsd > 0,
carryover_decay(tsd, h), 0))`.

### Analysis implications of each DGP step

<!-- From docs/44 Section 5 -->

| DGP step | Analysis implication |
|---|---|
| Graded coupling, anchored carryover | `Dbc` matches the coupling's decay; binary `Db` remains valid but loses power under carryover |
| AR(1) | Working covariance must allow serial correlation: `corCAR1` or structured with CR2; `us` with KR is slightly liberal and failed the stress test |
| Feasible strengths | Power comparisons at matched, feasible $c_{bm}$ |
| Pharmacodynamic model | `Dbc` from the response half-life |

The exploratory stress test selects F1: random intercept with
`corCAR1` residuals and CR2 inference, the baseline row, linear time,
`Dbc` at an assumed half-life of 0.5 week, and the expectancy term
`De` with its biomarker interaction.

### Why the analysis model does not match the DGP

<!-- From docs/24 "Why the analysis model doesn't match the DGP" -->

Four constraints make a component-matched analysis impractical.

1. **Identifiability.** The BR-PB-TV decomposition is identified by
   phase contrasts. A pure open-label trial cannot distinguish $BR$
   from $PB$, regardless of how the analyst writes the model.
2. **Degrees of freedom.** A component-matched analysis adds on the
   order of a dozen new parameters, half of them non-linear. At trial
   sample sizes of $N \in [30, 150]$, the variance of
   $\hat\beta_{bm:D}$ is dominated by participant and timepoint
   counts; adding parametric structure usually inflates it more than
   it reduces bias.
3. **Convergence.** Non-linear mixed models (`nlme::nlme` with
   self-starting Gompertz forms) are more sensitive to starting values
   and more likely to fail to converge on a single trial dataset than
   the linear-mixed analogue.
4. **Estimand robustness.** For predictive-biomarker validation, the
   regulator wants to know whether the biomarker predicts treatment
   response under trial conditions; a clean BR-versus-PB
   decomposition is not required. The simple `bm:Dbc` interaction is
   the right marginal estimand for the validation question.

The smallest useful step toward component-aware analysis is phase
indicators and their interactions with the drug term, schematically:

```r
Sx ~ bm + t + Dbc + phase + Dbc:phase
     + bm:Dbc + bm:Dbc:phase + (1 | ptID),
     correlation = corCAR1(...)
```

(Note: with three phase levels in the Hybrid design, `Dbc:phase` is
an exact linear combination of the intercept, the `phase` columns and
the `Dbc:phase` columns, so `Dbc:phase` is dropped automatically;
verified, `check-gompertz-rmanova-cycles.R`.)

---

## Derivations

<!-- one row per result: statement | CM label | verification script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the numbers agree? -->

---

## Limitations and open items

<!-- link each to the open-items register at docs/series/registry/open-items.csv -->

---

## Reproduce

<!-- one command per output -->

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | Setting and argument | Initial content move from docs/43, 40, 41, 42, 44, 45, 24, CM | Phase 2 content moves |
