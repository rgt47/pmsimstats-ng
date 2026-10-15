# Informative dropout and missingness

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Dropout hazard; MCAR, MAR, MNAR; interaction with design;
stress-test dropout cells.
**Prerequisites:** ch02, ch10

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

| Symbol | Meaning |
|---|---|
| $p^{\mathrm{flat}}_{it}$ | time-dependent dropout probability component |
| $p^{\mathrm{bias}}_{it}$ | symptom-trajectory-dependent dropout probability component |
| $\beta_0^{\mathrm{dr}}$ | overall dropout rate |
| $\beta_1^{\mathrm{dr}}$ | fraction of dropout attributable to symptom-biased mechanism |
| $\Delta_{Y,i}(t)$ | symptom change since baseline, $Y_{it} - Y_{i0}$ |
| $e_{b1}$ | exponent on the biased component (fixed at 2) |
| $d_i$ | first missing post-baseline visit ($\infty$ if complete) |

---

## Summary

<!-- evidence-labeled results go here -->

Paper 09 evaluates a censoring mechanism whose hazard depends on the
participant's own symptom trajectory, crossed with the four trial
designs. The marginal per-timepoint censoring probability is the union
of a flat (time-dependent) component and a biased (symptom-change-
dependent) component. The mechanism is MNAR for every dropout pattern
with $\beta_1^{\mathrm{dr}} > 0$ because missingness of $Y_{it}$
depends on the value of $Y_{it}$ itself.

Paper 09 documents a "happy accident" selection effect: participants
with strong open-label response drop out early, enriching the
randomized blocks for participants whose crossover contrasts are most
informative for the interaction test. This effect is design-dependent
and was not quantified per path in the version inspected here.

The stress test uses two simpler mechanisms (MCAR and MAR, 20%
marginal rate each) that are distinct from Paper 09's hazard.

---

## Construct box

*Data source: job 15
(`analysis/data/quick-sim/carryover-closed-form/strawman-stress/`
`summary-full-anchored.csv`). The Dropout column is not in the
standard construct box spec; it is included here because this
chapter's numerical results vary by dropout mechanism.*

| Process | Carryover | Coupling | Allocation | Design | N | Analysis ID | t_half_a (wk) | n_reps | Dropout |
|---|---|---|---|---|---|---|---|---|---|
| refined | anchored | graded | randomized | Hybrid | 70 | F1, S | 0.5 | 500 | none |
| refined | anchored | graded | randomized | Hybrid | 70 | F1, S | 0.5 | 500 | MCAR 20% |
| refined | anchored | graded | randomized | Hybrid | 70 | F1, S | 0.5 | 500 | MAR 20% |

Stored runs (job 15) used $t_{1/2}^a = 0.5$ week; the 0.2 week
policy applies to new runs.

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3.3 for this chapter:
- Sources: CM "Informative dropout" (4 subsections); paper 09 (with
  the matched-N correction); 46 Section 5; 47 dropout cells;
  Appendix I Section 2.
- New to write: most of the chapter (no whitepaper covers dropout).
-->

### The dropout hazard

<!-- From CM sec:dropout, lines 4360-4398 -->

Marginal per-timepoint censoring probability is the union of two
independent components:

$$
\Pr(\mathrm{drop}_{it}) = 1 - \big(1 - p^{\mathrm{flat}}_{it}\big)
  \big(1 - p^{\mathrm{bias}}_{it}\big),
$$

with:

- $p^{\mathrm{flat}}_{it} \propto t$ depending on elapsed time alone
  (calibrated to a marginal dropout fraction of
  $\beta_0^{\mathrm{dr}}(1 - \beta_1^{\mathrm{dr}})$).
- $p^{\mathrm{bias}}_{it} \propto t\,(\Delta_{Y,i}(t)+100)^{e_{b1}}$
  depending on the symptom change $\Delta_{Y,i}(t) = Y_{it} - Y_{i0}$
  since baseline (calibrated to $\beta_0^{\mathrm{dr}}
  \beta_1^{\mathrm{dr}}$), with $e_{b1} = 2$ throughout.

(Sources: Paper 09, `R/censordata.R`; inspected.)

The additive shift of 100 inside the biased component keeps its convex
argument positive across the realized range of symptom scores, so that
large improvement lowers dropout probability while worsening raises it.
Once a participant is censored at a timepoint, every later timepoint
for that participant is also set missing (monotone/absorbing dropout).

### Missingness classification

<!-- From CM sec:dropout, lines 4400-4413 -->

Because the flat component depends only on a fully observed design
variable (elapsed time), it induces ignorable missingness on its own.
Because the biased component's censoring indicator is a function of
$\Delta_{Y,i}(t)$, the quantity realized at, and thereby rendered
missing by, that timepoint, Paper 09 argues the mechanism is MNAR for
every dropout pattern with $\beta_1^{\mathrm{dr}} > 0$: missingness of
$Y_{it}$ depends on the value of $Y_{it}$ itself, the defining
condition for missing-not-at-random data.

No manuscript in the compendium gives a formal probability statement for
MCAR/MAR/MNAR beyond this argument; the classification rests on the
hazard construction alone.

### Interaction with design

<!-- From CM sec:dropout, lines 4415-4436 -->

The estimand is the same interaction test as elsewhere,
$\pi = \Pr(p < 0.05)$ for $H_0: \beta_{bm:D} = 0$, fit post-dropout.
Paper 09 documents a "happy accident" selection effect:

- Participants with unambiguous open-label response are the most
  readily resolved and most prone to drop out early.
- This enriches the later, randomized (crossover/discontinuation)
  blocks for participants whose open-label response was ambiguous,
  the subgroup whose crossover contrasts are most informative for
  the interaction test.
- This selection effect is design-dependent: under the OL design,
  where $D_{bc}$ is identically one throughout, there is no
  within-subject contrast for `bm:Dbc` at all, so dropout there can
  only reduce power.

A quantified, per-path decomposition of the happy-accident effect is
explicitly scoped for a future production-grade rerun and is not
delivered in the version inspected. An earlier mismatched-$N$ run had
produced an apparent robustness finding (the Hybrid design appearing
dropout-immune) that a matched-$N = 70$ rerun showed was a ceiling
effect rather than a property of the design.

### The MAR and MCAR mechanisms of the stress test

<!-- From CM sec:invf-mar, lines 4438-4477 -->

The stress test of `docs/46` uses two monotone mechanisms distinct from
Paper 09's hazard. With $d_i$ the first missing post-baseline visit
($d_i = \infty$ if complete), $t$ indexing post-baseline visits
$1, \dots, 8$, and $U_{it}$ independent uniforms:

$$
\begin{aligned}
\text{MCAR:}\quad & \Pr(d_i < \infty) = 0.2, \qquad d_i \mid d_i < \infty
  \sim \mathrm{Uniform}\{3, \dots, 8\}; \\
\text{MAR:}\quad & \operatorname{logit} p^{\mathrm{MAR}}_{it} = \mu_h
  + \beta_h\, z_{i,t-1}, \quad \beta_h = 1, \quad t = 3, \dots, 8, \\
  & \frac{1}{n_p}\sum_{i \in p}\Big[1 - \prod_{t=3}^{8}
  \big(1 - p^{\mathrm{MAR}}_{it}\big)\Big] = 0.2,
\end{aligned}
$$

where $p^{\mathrm{MAR}}_{it}$ is the hazard at visit $t$ given
presence at $t - 1$, $d_i$ is the first $t$ with
$U_{it} < p^{\mathrm{MAR}}_{it}$, $z_{i,t-1}$ is the score at the
previous visit standardized visit by visit across the $n_p$
participants of randomization path $p$, and $\mu_h$ is solved
numerically for each path of each trial.

Because the hazard depends only on the previously observed score, the
mechanism is MAR and likelihood analysis of the observed data is
ignorable for it, in contrast to the MNAR classification of Paper 09's
mechanism.

Note: the driver fixes $\beta_h = 1$. Because $\mu_h$ is calibrated
within each path of each trial to the realized scores, the 20% rate
holds in expectation conditional on those scores rather than as a
population parameter; $\mu_h$ varies across trials, and a
participant's hazard also depends on the other participants' observed
values: dropout depends only on observed data, but not through a fixed
per-participant model.

(Sources: `docs/46`, Section 5; `docs/51`, Feature 9; inspected,
`11-strawman-stress-test.R`, `dropout_visit()`.)

### New-to-write items

<!-- TODO: power table for the dropout stress-test cells: power by
design by dropout mechanism (none, MCAR, MAR) at the base and peak
cells of the stress test (docs/47, Table 4 or equivalent); with
matched N = 70 throughout -->

<!-- TODO: whether any analysis addresses MNAR adequately or whether
the compendium's analyses are all ignorable (likelihood or
least-squares based); state what the MNAR mechanism implies for bias
in the interaction coefficient under the biased component (qualitative
argument) -->

<!-- TODO: Appendix I Section 2: the dropout mechanism for the
appendix stress test (if different from docs/46) -->

<!-- TODO: reconciliation of the matched-N and unmatched-N runs of
Paper 09: which run is the current authoritative one? What were the
specific N values in each? Were the matched-N results published or
only documented in the compendium? -->

---

## Derivations

| Statement | Source | Script |
|---|---|---|
| MNAR classification for $\beta_1^{\mathrm{dr}} > 0$ | CM `eq:dropout-hazard`; Paper 09 | derived |
| MAR classification for stress-test mechanism | CM `eq:invf-mar`; `docs/46` Section 5 | derived |
| Monotone/absorbing censoring | `R/censordata.R`; inspected | |

---

## Numerical evidence

### Power attenuation under dropout

The stress-test dropout cells (MCAR and MAR at 20% marginal rate, job
15) document the power cost of dropout at $c_{bm} = 0.3$, $N = 70$
Hybrid design, analysis F1. Results below are from
`analysis/data/quick-sim/carryover-closed-form/strawman-stress/summary-full-anchored.csv`,
job 15 [inspected, summary-full-anchored.csv]:

| Dropout mechanism | F1 power | S power | MCSE (approx.) |
|---|---|---|---|
| None (base cell) | 0.424 | 0.499 | 0.022 |
| MCAR (20%) | 0.352 | 0.455 | 0.021 |
| MAR (20%) | 0.319 | 0.418 | 0.021 |

The F1 analysis power drops from 0.424 to 0.352 under MCAR (-17%)
and to 0.319 under MAR (-25%). The MAR attenuation is larger because
the dropout selection mechanism correlates with the outcome: subjects
with higher symptom scores (the ones whose data are most informative
for detecting the interaction) are more likely to drop out, leaving an
analysis sample that under-represents the high-end of the biomarker
distribution.

Under the S (symptom-only) analysis, the pattern is similar but the
absolute attenuation is smaller because S does not rely on the
between-subject biomarker contrast.

---

## Relation to the manuscripts

- **Paper 09:** Paper 09 defines the MNAR dropout hazard; documents
  the happy-accident selection effect; reports the matched-$N = 70$
  rerun correcting an earlier mismatched comparison.
- **`docs/46`:** The stress-test protocol, Section 5 specifies the
  MAR and MCAR mechanisms and the 20% marginal rate.
- **`docs/47`:** Stress-test results for the dropout cells.

---

## Limitations and open items

- The happy-accident effect has not been quantified per path or per
  design (open item, explicitly deferred in Paper 09).
- The matched-$N = 70$ Paper 09 run is confirmed in the compendium;
  whether its results were propagated to the paper is not verified
  here.
- The bias in $\hat\beta_{bm:D}$ under the MNAR mechanism has not
  been derived; the compendium records only the power
  ($\Pr(p < 0.05)$), not the bias.
- A formal MCAR/MAR/MNAR audit of all compendium mechanisms is not
  written; only the two mechanisms above are classified.

---

## Reproduce

```bash
# Stress-test dropout cells
Rscript analysis/scripts/quick-sim/\
  11-strawman-stress-test.R
# Paper 09 main run
bash tools/render.sh \
  analysis/report/09-dropout/report.Rmd
```

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | All | Phase 2 content move from CM sec:dropout, docs/46 Section 5, docs/47 | Whitepaper refactor |
| 2026-10-15 | Numerical evidence | GP-24: replace power-table TODO with dropout stress-test results at c_bm=0.3 N=70 | Phase 3 gap fill |
