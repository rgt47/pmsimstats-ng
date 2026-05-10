---
title: "Cycle-and-period design sweep: planning note"
author: "pmsimstats team"
date: "2026-05-07"
---

# Cycle-and-period design sweep: planning note

*2026-05-07 06:47 PDT*

This document specifies the parameter grid, simulation plan, and
analysis approach for a sensitivity sweep that quantifies how the
biomarker-treatment interaction test
($H_0: \beta_{\text{bm}:\text{Dbc}} = 0$) responds to changes in
the protocol-level design parameters that govern within-subject
treatment cycling: the number of on/off cycles, the on-drug and
off-drug period lengths within each cycle, the symmetry of the
on/off split, and the optional run-in and run-out phases.

The sweep is intended to inform protocol-design decisions for the
prazosin-PTSD aggregated N-of-1 application that motivates the
pmsimstats project [@hendricksonOptimizing2020]. It is the
**S13** sensitivity block; numbering continues from the existing
S1-S12 blocks documented at
`analysis/scripts/nof1-design-sensitivity/sensitivity/`.

## 1. Scientific motivation

Existing pmsimstats sweeps (S1-S12) examine sensitivity of the
**main drug effect** test under axis-wise perturbations of single
parameters around a reference design (Hendrickson hybrid). The
sweeps holds the protocol structure fixed and varies the data-
generating process. The complementary question, which the
existing infrastructure does not answer, is how the interaction
test responds to changes in the **protocol structure itself**:

1. **Number of cycles** $k$. More within-subject treatment
   switches contribute more within-subject contrasts to the
   estimator of $\beta_{\text{bm}:\text{Dbc}}$. Under the
   simplest balanced design with no carryover, the Fisher
   information on the interaction coefficient is approximately
   linear in $k$. Whether this scaling holds at realistic
   prazosin-PTSD parameters with carryover is the empirical
   question.
2. **Period lengths** $T_{\text{on}}, T_{\text{off}}$. The
   on-drug period must be long enough for the BR Gompertz
   trajectory to develop measurable signal; the off-drug period
   must be long enough for residual drug effect to wash out and
   restore the on-drug-versus-off-drug contrast. The two
   requirements compete at fixed total trial duration.
3. **Within-cycle asymmetry** $T_{\text{on}} \neq T_{\text{off}}$.
   Pharmacokinetic onset and washout typically operate on
   different timescales; the interaction-information harmonic
   mean argument (§4 below) predicts that symmetric periods are
   optimal at fixed total duration absent carryover, but
   carryover should bend the optimum toward longer off-drug
   periods.
4. **Across-cycle asymmetry** ($\Delta_k \neq 1$). Cycle 1 may
   reasonably differ in length from later cycles, particularly
   if it doubles as titration. The empirical question is whether
   such asymmetry materially affects power.
5. **Run-in and run-out phases**. Open-label titration before
   blinding (the OL+BDC pattern) and post-cycle observation
   contribute information differently from within-cycle
   observations. Including these as design parameters allows
   the sweep to characterise their power contribution.

A reviewer of any pmsimstats manuscript can reasonably ask, "if
the trialist reduces the number of cycles to fit within a
12-week window, by how much does the interaction power
decrease?" The S13 sweep is intended to answer such questions
quantitatively.

## 2. Design parameterisation

### 2a. Terminology

We distinguish two concepts that are conflated in much of the
clinical-trial methodology literature:

- **Run-in.** A period before treatment starts. Observations
  during run-in are pre-treatment baseline measurements.
  Drug is *not* given. Run-in observations have $\text{Db} = 0$
  and $\text{Dbc} = 0$ in the analysis model.
- **Open-label.** A period during which active drug is given
  and the subject is aware that the treatment is active.
  Open-label observations have $\text{Db} = 1$ and
  $\text{Dbc} = 1$ in the analysis model.

Run-in and open-label are orthogonal: a trial may have either,
both, or neither. The Hendrickson 2020 designs have open-label
phases (OL1, OL2 at weeks 4, 8 in the Hybrid; OL1-OL4 at
weeks 4-16 in the OL+BDC) but no run-in: every patient receives
active drug from week 0 onward.

The two phases also play different roles in identifying the
interaction coefficient $\beta_{\text{bm}:\text{Dbc}}$:

- Run-in observations contribute *directly* to identifying the
  interaction, because they are off-drug ($\text{Dbc} = 0$)
  observations with biomarker variation, contrasting against the
  later on-drug observations.
- Open-label observations contribute *indirectly* (through
  random-intercept precision, AR(1) parameter estimation, and
  linear-in-$t$ slope absorption) but not directly: their
  $\text{Dbc} = 1$ value is identical to the within-cycle
  on-drug observations, so they do not enlarge the on/off
  contrast.

The S13 sweep treats the two as separate parameters.

### 2b. Sampling density convention

We standardise on **one observation per week** throughout all
phases as the default sampling density for S13. This is dense
enough to track the Gompertz BR trajectory in the
prazosin-PTSD parameter regime, sparse enough to avoid
near-singularity in the AR(1) covariance when $\rho$ is large,
and matches the within-cycle density used in the Hendrickson
Hybrid blinded discontinuation phase. Under this default, the
number of observations per phase is directly proportional to the
phase duration: a 4-week on-drug period has 4 observations, a
2-week off-drug period has 2 observations. A separate Tier 2
block (§3d below) varies the density to characterise the
sampling-cost frontier.

### 2c. Helper

The current `buildtrialdesign()` interface takes a binary
`ondrug` vector at user-specified `timepoints`. This is
expressive (any sequence is encodable) but tedious for
protocol-level sweeps. We introduce a higher-level helper:

```r
build_cycle_design(
  T_runin     = 0,     # weeks of pre-treatment baseline (no drug)
  T_openlabel = 0,     # weeks of on-drug, unblinded titration
  k,                   # number of blinded on/off cycles
  T_on,                # blinded on-drug period (weeks); scalar
                       #   or length-k vector
  T_off,               # blinded off-drug period; scalar or length-k
  T_runout    = 0,     # post-cycle observation weeks
  order       = 'on_first',  # 'on_first' or 'off_first'
  density     = 1,     # observations per week (default 1/wk)
  expectancy  = 0.5    # blinded-phase expectancy
)
```

`build_cycle_design()` constructs the timepoints, on-drug binary
vector, and expectancy vector, then delegates to
`buildtrialdesign()`. The downstream simulation core
(`generate_data`, `lme_analysis`, `prepare_long_data`) consumes
the result unchanged.

A round-trip equivalence check is part of the implementation:
calling `build_cycle_design(T_openlabel = 8, k = 1, T_on = 0,
T_off = 4, T_runout = 8, density = 1)` with appropriate density
adjustment must reproduce the canonical Hendrickson Hybrid
schedule (open-label observations at weeks 4 and 8; blinded
discontinuation observations at weeks 9, 10, 11, 12; brief
crossover at weeks 16 and 20). The Hendrickson schedule uses
phase-specific densities (1/4 wk for open-label, 1/wk for
blinded discontinuation, 1/6 wk for brief crossover); a separate
helper `build_hendrickson_design()` will be provided for exact
Hendrickson reproduction.

## 3. Sweep design

### 3a. Tier 1: principal $k \times T_{\text{on}} \times T_{\text{off}}$ factorial

Symmetric across-cycle ($\Delta_k = 1$), focused on the founding
interaction question. Sampling density held at 1 observation per
week throughout all phases (so number of observations per phase
equals the phase duration in weeks); run-in held at 0 weeks (no
pre-treatment baseline, matching Hendrickson); open-label held
at 8 weeks (matching the Hendrickson Hybrid two-observation
titration).

| Axis | Values | Levels |
|---|---|---|
| Cycles $k$ | 1, 2, 3, 4 | 4 |
| On-drug period $T_{\text{on}}$ (weeks) | 2, 4, 6 | 3 |
| Off-drug period $T_{\text{off}}$ (weeks) | 2, 4, 6 | 3 |
| DGP architecture | A (mean moderation), B (MVN) | 2 |
| Carryover half-life $t_{1/2}$ (weeks) | 0, 0.5, 1.0 | 3 |
| Biomarker moderation $c_{\text{bm}}$ | 0, 0.45 | 2 |

Held fixed: $T_{\text{runin}} = 0$, $T_{\text{openlabel}} = 8$,
$T_{\text{runout}} = 0$, density = 1 / wk, $N = 70$.

Cell count: $4 \times 3 \times 3 \times 2 \times 3 \times 2 = 432$
combinations. **Feasibility filter.** Total subject-time per
participant is $T_{\text{runin}} + T_{\text{openlabel}} +
k \cdot (T_{\text{on}} + T_{\text{off}}) + T_{\text{runout}}$.
With the held-fixed values above, the constraint
"total trial duration $\leq 32$ weeks" eliminates the longest-
duration cells ($k = 4$ at $T_{\text{on}} = T_{\text{off}} = 6$
is $8 + 4 \times 12 = 56$ weeks; $k = 3$ at $T_{\text{on}} =
T_{\text{off}} = 6$ is 44 weeks). The constraint reflects the
upper end of patient-feasibility for N-of-1 trials in chronic
neurobehavioral conditions, where adherence and dropout become
limiting at trial durations beyond eight months
[@hoskins2015ptsd]. The remaining $\sim 264$ cells constitute
Tier 1 of the sweep.

Reference cell, used as the anchor for asymmetric extensions in
Tier 2: $k = 2$, $T_{\text{on}} = T_{\text{off}} = 4$,
Architecture B, $t_{1/2} = 0.5$, $c_{\text{bm}} = 0.45$,
$N = 70$, $T_{\text{openlabel}} = 8$. This is approximately the
Hendrickson-2020 Hybrid configuration, expressed in cycle-period
parameters with the standardised 1-per-week sampling density
extended across the cycles (where the published Hybrid uses
1 / wk during blinded discontinuation but coarser density during
open-label and brief crossover).

Sample size $N$ is fixed at 70 for Tier 1 to keep cell count
manageable; an $N = 35$ sweep is added in Tier 2.

### 3b. Tier 2a: within-cycle asymmetry sweep

At fixed $k = 2$ and matched total period length
$T_{\text{on}} + T_{\text{off}} = 8$ weeks, walk the asymmetry
ratio.

| Axis | Values | Levels |
|---|---|---|
| $(T_{\text{on}}, T_{\text{off}})$ | (2,6), (3,5), (4,4), (5,3), (6,2) | 5 |
| Architecture | A, B | 2 |
| $t_{1/2}$ | 0, 0.5, 1.0 | 3 |
| $c_{\text{bm}}$ | 0, 0.45 | 2 |

Cell count: $5 \times 2 \times 3 \times 2 = 60$. Sample size
fixed at $N = 70$. The expected pattern under the harmonic-mean
prediction (§4) is a maximum at $(4, 4)$ and symmetric decay; the
expected pattern under carryover is a minimum at $(6, 2)$
(insufficient washout) and a non-symmetric maximum biased toward
longer off-drug periods.

### 3c. Tier 2b: across-cycle asymmetry sweep

At fixed $k = 2$, $T_{\text{on}} = T_{\text{off}}$, walk the
relative period length of cycle 1 vs cycle 2 keeping total
duration constant at 16 weeks.

| Axis | Values | Levels |
|---|---|---|
| Cycle-length ratio $T_1 / T_2$ | 0.5, 0.75, 1.0, 1.5, 2.0 | 5 |
| Architecture | A, B | 2 |
| $t_{1/2}$ | 0, 0.5, 1.0 | 3 |
| $c_{\text{bm}}$ | 0.45 | 1 |

Cell count: $5 \times 2 \times 3 \times 1 = 30$.

### 3d. Tier 2c: sampling-density sweep

Vary the standardised sampling density (in observations per week)
at the reference $(k, T_{\text{on}}, T_{\text{off}})$ to
characterise the sampling-cost frontier and the AR(1)-bounded
marginal value of additional observations. The principal sweep
runs at 1/wk; this Tier 2c block extends below and above that
default.

| Axis | Values | Levels |
|---|---|---|
| Density (observations per week) | 0.5, 1, 2, 4 | 4 |
| Architecture | A, B | 2 |
| $t_{1/2}$ | 0, 0.5, 1.0 | 3 |
| $c_{\text{bm}}$ | 0.45 | 1 |

Cell count: $4 \times 2 \times 3 \times 1 = 24$. The 0.5 / wk
level halves the observation count and reflects monthly visits;
the 4 / wk level quadruples it and reflects daily symptom
self-reporting. Held fixed at the Tier 1 reference cell
otherwise.

### 3e. Tier 2d: $N \times k$ trade-off

A two-axis sweep showing whether more cycles per patient can
substitute for more patients at fixed total person-weeks. This
is the design-economics question for trialists choosing between
recruiting an additional participant and adding a cycle to the
existing protocol.

| Axis | Values | Levels |
|---|---|---|
| $N$ | 20, 35, 50, 70 | 4 |
| $k$ | 1, 2, 3 | 3 |
| Architecture | B | 1 |
| $t_{1/2}$ | 0.5 | 1 |
| $c_{\text{bm}}$ | 0.45 | 1 |

Cell count: $4 \times 3 = 12$.

### Sweep totals

| Tier | Cells | Replicates | Cell-replicates |
|---|---|---|---|
| 1: principal $k \times T_{\text{on}} \times T_{\text{off}}$ | 264 | 500 | 132,000 |
| 2a: within-cycle asymmetry | 60 | 500 | 30,000 |
| 2b: across-cycle asymmetry | 30 | 500 | 15,000 |
| 2c: sampling density | 24 | 500 | 12,000 |
| 2d: $N \times k$ trade-off | 12 | 500 | 6,000 |
| Total | 390 | | 195,000 |

At the Tier 1 production-run rate (approximately 13 cell-replicates
per second on the existing 8-core hardware, per the
carryover-sensitivity production timing), the full sweep
completes in approximately 4 hours. Dev-scale (50 replicates per
cell) finishes in approximately 24 minutes.

## 4. Theoretical expectations

Under Architecture A with no carryover, no AR(1) correlation, and
balanced symmetric cycles, the Fisher information on
$\beta_{\text{bm}:\text{Dbc}}$ scales approximately as

$$
I(\beta_{\text{bm}:\text{Dbc}}) \;\propto\; n \cdot k \cdot
\frac{T_{\text{on}} T_{\text{off}}}
     {T_{\text{on}} + T_{\text{off}}} \cdot \text{Var}(B_i)
$$

at leading order, where the harmonic-mean term arises from the
within-subject contrast variance in the matched-decay $D_{it}$
predictor. Three predictions follow:

1. **Linear in $k$.** Doubling $k$ at fixed
   $(T_{\text{on}}, T_{\text{off}})$ approximately doubles
   information, so $\beta_{\text{bm}:\text{Dbc}}$ standard error
   shrinks by $\sqrt{2}$.
2. **Symmetric optimum.** At fixed total period length
   $T_{\text{on}} + T_{\text{off}}$, the harmonic-mean term is
   maximised at $T_{\text{on}} = T_{\text{off}}$.
3. **Carryover bends the optimum.** When residual drug effect
   contaminates the early off-drug period, the effective $D_{it}$
   contrast is reduced. The optimal asymmetry is then
   $T_{\text{off}} > T_{\text{on}}$ to recover usable washout.

The empirical task is to confirm or contradict each prediction
under the Gompertz BR trajectory and AR(1) correlation that
pmsimstats actually implements.

## 5. Performance measures

For each cell, report the Morris-style metrics
[@morris2019simulation, Table 6] now standard in pmsimstats:

1. Power for $H_0: \beta_{\text{bm}:\text{Dbc}} = 0$ at
   $\alpha = 0.05$, with Monte Carlo SE.
2. Coverage of nominal 95 percent Wald CI, with MC SE.
3. Bias on the calibrated effect-size scale
   ($\theta_{\text{true}} = -c_{\text{bm}} \cdot \sigma_{BR}$),
   with MC SE.
4. Empirical SE.
5. Mean squared error.
6. Non-convergence rate per cell.

## 6. Visualisation

- **Tier 1 heatmap matrix.** $T_{\text{on}}$ on $x$,
  $T_{\text{off}}$ on $y$, one panel per $k$, two columns per
  carryover half-life. Color-encodes power.
- **Cycles-per-patient line plot.** Power versus $k$ at fixed
  $T_{\text{on}} = T_{\text{off}} = 4$, faceted by carryover.
- **Asymmetry curve.** Power versus $T_{\text{on}} / T_{\text{off}}$
  ratio at fixed total period length, faceted by carryover.
- **$N$-versus-$k$ trade-off plot.** Power as a function of
  total person-weeks $(N \cdot k \cdot (T_{\text{on}} +
  T_{\text{off}}))$, with $(N, k)$ combinations labelled. The
  question is whether the points fall on a single curve (perfect
  trade-off) or whether $k$ dominates $N$ or vice versa.
- **Lasagna plot of all 375 cells.** Same convention as
  `analysis/figures/02-lasagna-power.pdf`.

## 7. Implementation plan

1. **Helper.** `build_cycle_design()` in
   `analysis/scripts/nof1-design-sensitivity/simulation-core.R`
   (or a shared file `analysis/scripts/_shared/design-helpers.R`
   if multiple sweeps will reuse it). Include unit tests that
   round-trip against the existing OL, CO, Hybrid, and OL+BDC
   presets.
2. **Driver.** New file
   `analysis/scripts/nof1-design-sensitivity/13-cycle-period.R`
   following the pattern of the existing 01-12 sensitivity-block
   drivers. Writes
   `analysis/scripts/nof1-design-sensitivity/output/13-cycle-period.rds`.
3. **Summariser.** Extend
   `analysis/scripts/nof1-design-sensitivity/summarise.R` (or the
   carryover-sensitivity equivalent) to handle the new block.
   Same Morris-compliant output columns as the existing
   summariser.
4. **Figures.** New file
   `analysis/scripts/nof1-design-sensitivity/13-cycle-period-figures.R`
   producing the heatmap matrix, line plot, asymmetry curve,
   trade-off plot, and lasagna plot.
5. **Manuscript prose.** Add §S13 to
   `analysis/report/05-nof1-design-sensitivity/report.Rmd` with
   abstract revision to mention the design-parameter sweep, a new
   Methods subsection, and a Results subsection per the figure
   set.
6. **Reproducibility.** Production run records master seed,
   commit SHA, R version, package versions in the output
   metadata, per the conventions in
   `docs/feedback_morris_simulation_standards.md`.

## 8. Risks and what they would mean

- **Linear-in-$k$ prediction fails.** If doubling $k$ does not
  approximately double information, the Fisher-information
  derivation in §4 is missing a term. Most likely candidate:
  AR(1) correlation reduces the effective independence of
  cycle-level contrasts. Diagnostic: rerun the Tier 1 grid at
  $\rho = 0$ (compound symmetry of within-factor) to isolate the
  AR(1) contribution.
- **Symmetric-optimum prediction fails even without carryover.**
  If the asymmetry sweep at $t_{1/2} = 0$ does not peak at
  $(4,4)$, the Gompertz curvature is contributing a non-
  negligible asymmetric term to the Fisher information. This
  would be a material finding that overturns the textbook
  recommendation for symmetric crossovers.
- **Carryover does not bend the optimum.** If the asymmetry sweep
  at $t_{1/2} > 0$ also peaks near symmetric, the carryover
  effect on the optimum design is smaller than the direct effect
  on power. This would simplify trial-design recommendations:
  symmetric is good across the realistic carryover range.
- **Across-cycle asymmetry has near-zero effect.** If the Tier
  2b sweep is essentially flat, then the choice of cycle-1
  versus cycle-2 length is not load-bearing for the interaction
  test. This would be a simplifying finding and would justify
  treating $T_{\text{on}}$ and $T_{\text{off}}$ as scalars
  throughout the manuscript.
- **$N$-versus-$k$ are not interchangeable.** If the Tier 2d
  trade-off plot shows that $N$ dominates $k$ (or vice versa) at
  fixed total person-weeks, the design-economics conclusion is
  unambiguous. If the curves coincide, the trialist is free to
  choose based on practical considerations (recruitment cost,
  adherence).

## 9. Reporting hooks for the manuscript

The S13 results should support three direct statements in the
nof1-design-sensitivity manuscript:

1. A trade-off curve quantifying how much per-patient duration
   ($k \cdot (T_{\text{on}} + T_{\text{off}})$) the trial must
   accept to reach 80 percent power at the prazosin-PTSD design
   parameters.
2. A recommendation on cycle count: at the prazosin-PTSD
   reference parameters, the marginal value of adding the third
   cycle is X percentage points of power, the fourth cycle is
   Y, and the recommended $k$ is whichever delivers the
   target power most economically.
3. A statement on asymmetry: whether the textbook symmetric
   recommendation is robust to carryover, or whether a specific
   asymmetric form is preferred for the prazosin-PTSD biology.

## 10. Resolved decisions

- Open-label duration as a sweep axis: **resolved
  (2026-05-07).** Held fixed at $T_{\text{openlabel}} = 8$
  weeks, matching the Hendrickson Hybrid two-observation
  titration. Not swept in S13.
- True run-in (pre-treatment baseline) as a sweep axis:
  **resolved (2026-05-07).** Held fixed at $T_{\text{runin}} =
  0$ weeks, matching the Hendrickson convention. Not swept in
  S13. Note that run-in observations contribute *directly* to
  identifying the interaction coefficient (clean off-drug
  observations unaffected by carryover); a future companion
  sweep could quantify that contribution if the question
  becomes relevant for protocol design.
- Analysis-specification axis: **resolved (2026-05-07).** Held
  fixed at A2 (matched-decay) for Tier 1, consistent with the
  carryover-sensitivity manuscript's finding that A2 dominates
  A1 and A3 across the parameter space.
- Parallel-RCT comparator: **resolved (2026-05-07).** Not
  included. S13 is concerned exclusively with the within-subject
  design space (cycle structure, period lengths, run-in/run-out).
  The parallel-RCT versus N-of-1 contrast is the responsibility
  of the broader nof1-design-sensitivity manuscript and is not
  re-evaluated here.
All open questions for S13 are now resolved. The implementation
plan in §7 can proceed without further pre-implementation
decisions.

## References

- [@hendricksonOptimizing2020] Hendrickson RC, Thomas RG,
  Schork NJ, Raskind MA. Optimizing aggregated N-of-1 trial
  designs for predictive biomarker validation: statistical
  methods and practical considerations. *Frontiers in Digital
  Health*. 2020;2:13.
- [@hoskins2015ptsd] Hoskins M, Pearce J, Bethell A, et al.
  Pharmacotherapy for post-traumatic stress disorder: systematic
  review and meta-analysis. *British Journal of Psychiatry*.
  2015;206(2):93-100.
- [@morris2019simulation] Morris TP, White IR, Crowther MJ.
  Using simulation studies to evaluate statistical methods.
  *Statistics in Medicine*. 2019;38(11):2074-2102.
  doi:10.1002/sim.8086.

---

*Rendered on 2026-05-07 at 07:14 PDT.*<br>
*Source: ~/Dropbox/prj/alz/10-pmsimstats-ng/pmsimstats-ng/docs/27-cycle-period-design-sweep-plan.md*
