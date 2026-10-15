# Carryover on the mean

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Anchored vs compounding carryover; the pharmacodynamic
half-life; restart on re-exposure.
**Prerequisites:** ch05

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
| $\phi(t_{sd})$ | carryover decay function; $\phi(t_{1/2}) = 0.5$ for three of four forms |
| $\mu_{BR,\text{last}}$ | drug-response mean at the most recent on-drug visit |
| $\lambda = \ln 2 / t_{1/2}$ | exponential decay rate |
| $s_f$ | scale factor on the carryover exponent (2024 upstream code) |
| $k$ | Weibull shape parameter |
| $p$ | power decay exponent |

---

## Summary

<!-- evidence-labeled results go here -->

The carryover half-life in this compendium is a pharmacodynamic
(response) half-life, not a plasma half-life. For prazosin, the plasma
half-life is approximately 2.3 hours; the behavioral or physiological
carryover that matters at weekly visits is a separate, slower process
set by the turnover of the response.

The refined process uses anchored exponential decay:
$\mu_{BR,t} = \mu_{BR,\text{last}} \cdot 2^{-t_{sd,t}/t_{1/2}}$,
decaying from the last on-drug value. The published code compounds the
decay across consecutive off-drug visits, running faster than the
stated half-life from the second off-drug visit onward (derived;
`docs/53`, Sections 4 and 5).

The anchored and compounding forms agree at the first off-drug visit
of each run and diverge only at the second and later visits. The
largest mean shortfall on the Hybrid design at $t_{1/2} = 1$ week is
1.273 points (verified, `docs/53`, Table 5).

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3.3 for this chapter:
- Sources: CM "Carryover on the mean" (3 subsections); 53 Section 5;
  51 Features 3 to 4; 38 Sections 2 to 3; 40 Section 3; 24
  Component 1.

- New to write: the reconciliation of a 2.3-hour plasma half-life
  with a 0.5 to 1 week carryover half-life.
-->

### Two kinds of carryover

<!-- From docs/38 Section 2 and docs/51 Feature 3 -->

**Pharmacokinetic carryover** is drug still in the body. For prazosin
in 24 normotensive subjects, the plasma half-life is approximately
2.3 hours (Hobbs et al. 1978). Five half-lives -- about 12 hours --
remove 97% of the drug. Hendrickson et al.'s half-lives of 0.1 and
0.2 weeks (17 and 34 hours) are on this scale or slightly above it.

**Pharmacodynamic (disease-state) carryover** is improvement that
persists after the drug has gone: consolidated sleep, fewer
nightmares, behavioral learning, recovery of a physiological set
point. Its time course is set by the turnover of the response, not by
drug elimination. Effect-compartment models and indirect-response
models (Dayneka et al. 1993; Sheiner et al. 1979) are the standard
quantitative treatments. In indirect-response models the response
returns to baseline slowly after the drug is gone.

For a symptom outcome measured weekly or less often, pharmacokinetic
carryover is negligible and the carryover that matters is
pharmacodynamic. The half-life parameter of the simulation should
therefore be read as the response half-life $t_{1/2}^R$.

<!-- TODO: reconciliation of the 2.3-hour prazosin plasma half-life
with the 0.5 to 1 week carryover half-life assumed in the
simulation: the plasma compartment equilibrates quickly; the
behavioral consolidation of reduced nightmares and improved sleep
architecture may persist for days; a response half-life of 0.5 week
= 84 hours is 36-fold longer than the plasma half-life; this is
consistent with indirect-response models (e.g. Dayneka et al. 1993)
where the response half-life is governed by the turnover rate of the
target process (here REM-sleep architecture) rather than by drug
exposure; the appropriate value for prazosin in PTSD is empirically
unknown and should be treated as uncertain -->

### Decay families

<!-- From CM "Carryover decay families", lines 1116-1156 -->

Four decay forms are implemented (`carryover_decay()` in
`implementations/nof1power/R/carryover_models.R`). The first three
satisfy $\phi(t_{1/2}) = 0.5$; the power form does not.

- **Exponential:** $\phi_\text{exp}(t_{sd}) = \exp(-\lambda t_{sd})$,
  $\lambda = \ln 2 / t_{1/2}$
- **Linear:** $\phi_\text{lin}(t_{sd}) = \max(0, 1 -
  t_{sd}/(2 t_{1/2}))$
- **Weibull:** $\phi_\text{weib}(t_{sd}) = \exp(-(\lambda_w
  t_{sd})^k)$, $\lambda_w = (\ln 2)^{1/k}/t_{1/2}$
- **Power:** $\phi_\text{pow}(t_{sd}) = \max(0, (1 -
  t_{sd}/(3 t_{1/2}))^p)$

The power form gives $\phi(t_{1/2}) = (2/3)^p$, which equals 0.5
only at $p = \ln 2 / \ln 1.5 \approx 1.71$; at the default $p = 1.5$
it is 0.544 (derived; corrected 2026-10-08).

The Weibull family recovers the exponential at $k = 1$; $k < 1$
produces a heavier tail, $k > 1$ a lighter one. Paper 02 tests
$k \in \{0.25, 0.5, 2.0, 4.0\}$ against the exponential default.
The advantage of the exposure-weighted analysis widens under
light-tailed shapes and narrows to insignificance under the two
heavy-tailed shapes (Paper 02, Abstract).

### Anchored versus compounding carryover

<!-- From CM "Anchored versus compounding carryover", lines
1158-1223 -->

The published process uses a carryover recursion that decays the
previous visit's already-decayed mean by the cumulative time since
discontinuation:

$$
\mu_{BR,t} = \mu_{BR,t-1}\, 2^{-t_{sd,t}/t_{1/2}}.
$$

Because $t_{sd}$ is the cumulative time since discontinuation, the
recursion applies $\prod_{l \le j}\phi(t_{sd,l})$ instead of
$\phi(t_{sd,j})$ at the $j$-th consecutive off-drug visit. It
over-decays from the second off-drug visit onward (Paper 01, Appendix
B.3; derived).

The refined process uses anchored carryover, decaying from the last
on-drug value:

$$
\mu_{BR,t} = \mu_{BR,\text{last}}\, 2^{-t_{sd,t}/t_{1/2}},
$$

with $\mu_{BR,\text{last}}$ the drug-response mean at the most recent
on-drug visit (`R/generateData.R`, inspected; Paper 01, Methods,
refinement step 3).

At the $k$-th consecutive off-drug visit, the compounding form gives
$\mu_\text{last} \cdot 2^{-\sum_{j=1}^k t_{sd,j}/t_{1/2}}$
against the anchored $\mu_\text{last} \cdot 2^{-t_{sd,k}/t_{1/2}}$.
On Hybrid path 1 at $t_{1/2} = 1$ week ($\mu_\text{last} = 10.187$):

| Visit | Week | $t_{sd}$ | Compounding | Anchored | Ratio |
|---|---|---|---|---|---|
| BD3 | 11 | 1 | 5.094 | 5.094 | 1.000 |
| BD4 | 12 | 2 | 1.273 | 2.547 | 0.500 |
| COp | 20 | 4 | 0.268 | 0.268 | 1.000 |

(Verified, `docs/53` Section 6.2; the two forms agree at the first
off-drug visit of each run.)

A paired simulation (Hybrid, $N = 70$, $c_{bm} = 0.25$, 500
replicates) found that replacing the recursion by the anchored form
changes the interaction estimate by at most 0.005 and at most 4 of
500 test decisions at $t_{1/2} \le 2$ weeks (`docs/53`, Section 8.3).

### The 2024 revision: a scale factor on the exponent

<!-- From CM "The 2024 revision of the published code: a scale factor",
lines 1225-1286 -->

The 2024 upstream code (HEAD `06dac83`) adds a scale factor $s_f$
(default 2) to the carryover exponent:

$$
\mu_{BR,t} = \mu_{BR,t-1}\, 2^{-s_f\, t_{sd,t}/t_{1/2}}.
$$

At $s_f = 2$ even the first off-drug visit of a run is understated
by the factor $2^{-t_{sd,1}/t_{1/2}}$. The scale factor acts on the
mean channel; its counterpart in the correlation matrix is the
ratio-graded coupling of ch08.

### Feature 3: restart on re-exposure

<!-- From docs/51 Feature 4 -->

When the drug is restarted, the time on drug $t_{od}$ resets to zero
and the response rebuilds from zero (Gompertz clock restart). The
refined process, the package and the published process all use this
rule. The alternative -- restart from the residual response -- is a
turnover-model feature (ch14, `docs/38` Section 3) and is not
implemented.

For a drug with residual response at re-exposure, the zero-restart
rule understates the on-drug response in the crossover period: on
Hybrid path 1 a participant who was near full response at week 10 is
credited with only $G(4) = 4.281$ of 10.99 at week 16, 39% of the
maximum (verified, `docs/40` Section 4).

---

## Derivations

| Statement | Source | Script |
|---|---|---|
| Compounding at $k$-th off-drug visit: $\mu_\text{last} \cdot 2^{-\sum t_{sd,j}/t_{1/2}}$ | CM `eq:cg-anchor`; Paper 01 Appendix B.3 | `check-historical-results.R` |
| Anchored and compounding agree at first off-drug visit of each run | derived | `docs/53` Section 6.2 |
| Power-form nominal half-life error: $\phi(t_{1/2}) = (2/3)^p \ne 0.5$ at $p = 1.5$ | derived (2026-10-08) | |

---

## Numerical evidence

From `docs/53` Table 5: largest published-to-intended shortfall over
all off-drug visits by design and half-life:

| Design | $t_{1/2}$ | Smallest ratio | Largest shortfall (points) |
|---|---|---|---|
| Hybrid | 0.5 | $2.4 \times 10^{-4}$ | 0.478 |
| Hybrid | 1 | 0.0156 | 1.273 |
| Hybrid | 2 | 0.125 | 2.238 |
| CO | 1 | $3.1 \times 10^{-5}$ | 0.262 |
| CO | 2 | 0.0055 | 1.044 |

(Verified, `docs/53` Table 5.)

---

## Relation to the manuscripts

- **Paper 01 Methods:** states anchored carryover as refinement step 3;
  Appendix B.3 derives the compounding effect.
- **Paper 02:** tests Weibull shapes $k \in \{0.25, 0.5, 2.0, 4.0\}$.
- **Paper 05:** states the carryover recursion; stored results may have
  been produced under the recursion rather than anchored carryover
  (not verified).

---

## Limitations and open items

- The appropriate pharmacodynamic half-life for prazosin in PTSD is
  not known; the simulation should treat it as uncertain.
- Restart from the residual is pharmacologically preferable for
  drugs without re-titration requirements but is not implemented.
- The 2024 upstream code's scale factor is documented but not used
  in the refined process.

---

## Reproduce

```bash
# docs/53 tables (requires pmsimstats checked out at 58b32a9)
Rscript docs/53-carryover-recursion-in-published-code.md  # code block
```

---

## References

- Dayneka NL, Garg V, Jusko WJ (1993). Comparison of four basic
  models of indirect pharmacodynamic responses. J Pharmacokin Biopharm
  21:457-478.
- Hobbs DC, Twomey TM (1978). Prazosin pharmacokinetics in man.
  J Clin Pharmacol 18:402-406.

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | All | Phase 2 content move from docs/38, docs/51, docs/53, docs/40, CM sec:carryover-mean | Whitepaper refactor |
