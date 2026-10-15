# Biomarker coupling under carryover

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Coupling profile g_t; step, graded and ratio rules;
the pharmacokinetic argument; the coupling gap G.
**Prerequisites:** ch07

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

*Local symbols used in this chapter only. For project-wide notation
see `analysis/report/NOTATION.md`.*

| Symbol | Meaning |
|---|---|
| $g_t$ | coupling: fraction of $c_{bm}$ present at visit $t$ |
| $\bar g_\text{on}$, $\bar g_\text{off}$ | mean coupling over on-drug and off-drug post-baseline visits |
| $\mathcal G = \bar g_\text{on} - \bar g_\text{off}$ | coupling gap |
| $s_f$ | scale factor on the coupling exponent (2024 code) |

---

## Summary

<!-- evidence-labeled results go here -->

The coupling $g_t = \mathrm{Cor}(B_i, BR_{it})/c_{bm}$ measures the
fraction of the on-drug biomarker-response association present at
visit $t$.

The interaction estimand for E9 is:

$$
\beta_{bm:D} = -c_{bm}\,\frac{\sigma_{BR}}{\sigma_{bm}}\,
  (\bar g_\text{on} - \bar g_\text{off}) = -c_{bm}\,
  \frac{\sigma_{BR}}{\sigma_{bm}}\,\mathcal{G}.
$$

Power is a strictly increasing function of the coupling gap $\mathcal
G$ and does not depend on how the coupling is distributed over visits,
only on its average (`docs/37`, Proposition 1 and Section 5.4;
derived).

**Step rule (published):** $g_t = \mathbf{1}[\mu_{BR,t} \ne 0]$. With
any carryover, $\bar g_\text{off} = 1$, so $\mathcal G = 0$ and the
interaction vanishes at every positive half-life (derived; verified,
`docs/55` Section 5).

**Graded coupling (refined process):** $g_t = 1$ on drug,
$\phi_t = 2^{-t_{sd,t}/t_{1/2}}$ off drug after exposure, 0 before
first exposure. The gap falls smoothly from 1 at $t_{1/2} = 0$ to 0
as $t_{1/2} \to \infty$. At the short prazosin pharmacological
half-life (0.1 week), power under graded coupling is within 0.001
of the no-carryover case (verified, `docs/37`).

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3.3 for this chapter:
- Sources: CM "The estimand and the coupling rules", "The
  pharmacokinetic argument for graded coupling", "ratio-graded",
  "Shape of the coupling decay"; 55 Sections 4 to 5; 37 Section
  4.2; 51 Feature 5; Appendix I Section 7.

- New to write: the answer to "is graded decay structurally required"
  (02 page 5, notes.md); the trait-against-state variance question
  as an open modeling choice.
-->

### The four coupling rules

<!-- From CM "The estimand and the coupling rules", lines 1478-1525
and docs/55 Section 5 -->

Four rules set the off-drug interaction (Paper 01, Table 1; `docs/37`
Section 2):

- **Step (published):** $g_t = \mathbf{1}[\mu_{BR,t} \ne 0]$
- **Graded:** $g_t = 1$ on drug, $\phi_t = 2^{-t_{sd,t}/t_{1/2}}$
  off drug after exposure, 0 before first exposure
- **Mean moderation:** $g_t = 0$; $BR_{it}$ is shifted by
  $\beta_{bm}\sigma_{BR} b_i h_t$, $h_t = D_{it}$
- **Decayed mean moderation:** $g_t = 0$; $h_t = 1$ on drug,
  $\phi_t$ off drug

Off-drug visits before first exposure carry coupling and moderation 0
under every rule. This pre-exposure zero is a condition, not a
detail: a rule "coupling equals the decay weight" with
$e^{-\lambda \cdot 0} = 1$ at $t_{sd} = 0$ would give full coupling
at baseline and at placebo-first crossover visits.

### The coupling gap determines the estimand

<!-- From docs/37 Section 4.2 and docs/55 Section 5 -->

For the paired-difference statistic E9, the interaction slope is:

$$
\beta_{bm:D}^{\text{E9}} = -c_{bm}\frac{\sigma_{BR}}{\sigma_{bm}}
  \big(\bar g_\text{on} - \bar g_\text{off}\big).
$$

The interaction slope is the biomarker slope on drug minus the
biomarker slope off drug, proportional to the coupling gap $\mathcal
G = \bar g_\text{on} - \bar g_\text{off}$.

**Step coupling at any positive half-life:** every visit after first
exposure has $g_t = 1$ in the Hybrid and OL+BDC designs, so
$\bar g_\text{on} = \bar g_\text{off} = 1$, the gap is 0, and
$\beta_{bm:D}^{\text{E9}} = 0$. E9b (counting baseline as one more
off-drug visit) retains a remnant of $c_{bm}\sigma_{BR}/(n_\text{off}
+ 1)\sigma_{bm}$, about one-quarter of the full slope on Hybrid paths
1 and 2 (derived; `docs/55` Section 5).

**Graded coupling:** $\bar g_\text{off} = n_\text{off}^{-1}
\sum_\text{off} 2^{-t_{sd,t}/t_{1/2}}$ is strictly increasing in
$t_{1/2}$ and tends to 0 as $t_{1/2} \to 0$. At the short prazosin
half-life of 0.1 week every off-drug visit has $t_{sd} \ge 1$ week,
so $\bar g_\text{off} \le 2^{-10} \approx 0.001$ and the gap is
above 0.999. As $t_{1/2} \to \infty$, $\bar g_\text{off} \to 1$: the
step rule is the infinite-half-life limit of graded coupling.

**Power depends only on the gap.** $\mathrm{Var}(\Delta)$ does not
depend on $t_{1/2}$ (`docs/37`, Proposition 1; derived; verified).
Power is strictly increasing in $|\beta_{bm:D}|$, hence in the gap
(`docs/37`, Section 5.4).

### The pharmacokinetic argument for graded coupling

<!-- From CM "The pharmacokinetic argument for graded coupling",
lines 1617-1692 -->

**The spread argument.** On drug, joint normality gives
$E(BR_{it} \mid b_i = \pm 1) = \mu_{BR,t} \pm c_{bm}\sigma_{BR}$,
a spread $\Delta_\text{on} = 2c_{bm}\sigma_{BR}$ between participants
one standard deviation above and below the biomarker mean. If
carryover is pharmacokinetic, each participant's drug effect decays
by the same factor $\phi_t = 2^{-t_{sd,t}/t_{1/2}}$, so the
off-drug conditional means are $\phi_t(\mu_{BR,t_0} \pm
c_{bm}\sigma_{BR})$ and the spread is $\Delta_\text{off} =
\phi_t \Delta_\text{on}$. Matching this to the conditional-normal
spread $2 c_{bm} g_t \sigma_{BR}$ gives $g_t = \phi_t$ (derived;
`docs/23`, Sections 5.2-5.4).

| Visit, rule | $c_{bm}g_t$ | $b = +1$ | $b = -1$ | Spread/mean |
|---|---|---|---|---|
| On drug (week 10) | 0.50 | 14.19 | 6.19 | 0.785 |
| Off drug, graded ($\lambda = \ln 2/t_{1/2}$) | 0.25 | 7.09 | 3.09 | 0.785 |
| Off drug, zero coupling | 0 | 5.09 | 5.09 | 0 |
| Off drug, full coupling (step rule) | 0.50 | 9.09 | 1.09 | 1.57 |

(Hybrid path 1, $t_{1/2} = 1$ week, $c_{bm} = 0.5$; verified,
`check-historical-results.R` block G.)

Only graded coupling keeps the ratio of the biomarker spread to the
mean response at its on-drug value. Under the step rule at
$t_{1/2} = 0.1$ week, the drug-response mean one week after
discontinuation is $10.19 \cdot 2^{-10} \approx 0.010$, but the step
rule retains $c_{bm} = 0.6$, giving conditional means at $b = \pm 1$
of $+4.81$ and $-4.79$ around a near-zero mean (derived; verified,
block G).

Two qualifications of the spread argument:

1. The premise (a common decay of each participant's conditional mean)
   is decayed mean moderation. The argument fixes the conditional mean
   and says nothing about the conditional variance. The construct
   keeps $\mathrm{Var}(BR_{it} \mid B_i) = \sigma_{BR}^2(1 -
   c_{bm}^2 g_t^2) \approx \sigma_{BR}^2$ off drug, whereas a decay
   of the whole drug-dependent deviation would shrink it by $\phi_t^2$
   (`docs/23`, Section 5.8).
2. The pharmacokinetic argument justifies an exponential coupling
   decay at the same rate as the carryover. Other decay shapes are
   possible (Weibull family); see "Shape of the coupling decay" below.

### The ratio-graded rule (2024 code)

<!-- From CM "The 2024 revision: ratio-graded coupling", lines
1585-1616 -->

In the 2024 upstream code an off-drug visit after exposure receives
the coupling $c_{bm}\,\mu_{BR,t}/\mu_{BR,t-1}$. The history cancels,
giving $g_t = 2^{-s_f\,t_{sd,t}/t_{1/2}}$ at every off-drug visit
of a run -- graded coupling at effective half-life $t_{1/2}/s_f$, with
the default $s_f = 2$ halving the half-life. The first visit is
never assigned ($g_1 = 0$), because the loop starts at the second
visit (Section 3.4 of the source; inspected,
`vendored-hendrickson-generateData.R`).

### Shape of the coupling decay

<!-- From CM "Shape of the coupling decay", lines 2919-2981 -->

If the off-drug coupling followed the Weibull family
($g_t = \exp\{-(\lambda_w t_{sd,t})^k\}$), the shape enters the
ceiling only through the off-drug entries of $g$. The package decays
the correlation exponentially only; Weibull coupling is a
counterfactual.

At $t_{1/2} = 0.5$ week the `covar` ceiling at $\rho = 0.7$ falls
monotonically with $k$ (heavier tail, higher ceiling):

| $k$ | 0.25 | 0.5 | 1 (exp.) | 2 | 4 |
|---|---|---|---|---|---|
| Hybrid, $t_{1/2} = 0.5$ | 0.543 | 0.531 | 0.506 | 0.469 | 0.456 |

(Verified, `docs/34`, Section 5.6.)

### New-to-write items

<!-- TODO: the answer to "is graded decay structurally required":
Paper 02 page 5 and notes.md argue that graded decay is not
structurally required -- a step coupling that detects the pre-
exposure zero correctly (zero before first exposure, one on drug, one
off drug wherever the mean is nonzero) is internally consistent, and
the power collapse is a feature of the step-rule's specific definition
rather than of step functions in general. Graded coupling is preferred
on pharmacokinetic grounds, not on structural grounds. State this
clearly and cite Paper 02 page 5. -->

<!-- TODO: Appendix I Section 7: off-drug coupling and the coupling
gap under AR(1) -- not yet written, leave TODO -->

<!-- TODO: the trait-against-state variance question as an open
modeling choice: the biomarker is treated as a fixed participant
trait (no within-person variation in $B_i$); in reality biomarker
levels fluctuate over time. If the biomarker is a state variable
rather than a trait, its correlation with $BR$ at a given visit
would depend on the joint dynamics of both, not just on the drug
effect. Characterize this as an unresolved modeling choice. -->

---

## Derivations

| Statement | Source | Script |
|---|---|---|
| E9 slope = $-c_{bm}(\sigma_{BR}/\sigma_{bm})\mathcal G$ | `eq:cg-estimand`, `docs/37` Prop. 2 | `check-coupling-closed-form.R` |
| $\mathrm{Var}(\Delta)$ does not depend on $t_{1/2}$ | `docs/37` Prop. 1 | `check-coupling-closed-form.R` |
| Step rule gap = 0 at every positive $t_{1/2}$ | derived from step definition | `docs/55` Section 5 |
| Graded = spread-preserving off-drug coupling | CM `eq:inva-pk` | `check-historical-results.R` block G |

---

## Numerical evidence

From `docs/36` Table 9, the power collapse under step coupling
(configuration A) versus graded coupling (configuration C) in the
Hybrid design at $c_{bm} = 0.25$:

| Config. | $t_{1/2} = 0$ | 0.1 | 0.2 | 0.5 | 1 |
|---|---|---|---|---|---|
| A (step, CS) | 0.59 | 0.08 | 0.08 | 0.09 | 0.10 |
| C (graded, sep. AR(1)) | 0.51 | 0.48 | 0.52 | 0.41 | 0.28 |

(`corCAR1` analysis; verified.)

---

## Relation to the manuscripts

- **Paper 01:** defines the four coupling rules; derives the
  coupling gap and the estimand.
- **Paper 02:** evaluates Weibull shapes for the coupling decay.
- **`docs/37`:** closed-form power as a function of the coupling gap.
- **`docs/55`:** the step-coupling defect in the published code.

---

## Limitations and open items

- The trait-versus-state nature of the biomarker is unresolved:
  within-person biomarker variation is not modeled.
- Whether graded decay is structurally required or merely preferred on
  pharmacokinetic grounds is not explicitly stated in any manuscript
  (open item; see TODO above).
- Appendix I Section 7 (off-drug coupling under AR(1)) is not yet
  written.

---

## Reproduce

```bash
Rscript analysis/scripts/compendium-mathematics/\
  check-coupling-closed-form.R
Rscript analysis/scripts/compendium-mathematics/\
  check-historical-results.R
```

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | All | Phase 2 content move from CM sec:cg-estimand, CM sec:inva-pk, CM sec:inve-ce-shape, docs/37, docs/55, docs/51 | Whitepaper refactor |
