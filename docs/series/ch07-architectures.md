# Interaction architectures: mean, covariance and dual-channel

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Mean, covariance and dual-channel architectures;
estimands; latent-variable reading.
**Prerequisites:** ch06

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
| $\beta_{bm}$ | mean-moderation parameter: coefficient scaling BR by $b_i$ |
| $c_{bm}$ | covariance-moderation strength: $\mathrm{Cor}(B_i, BR_{it})$ on drug |
| $c_{bm,A}$, $c_{bm,B}$ | mean- and covariance-channel weights (dual-channel) |
| $h_t$ | mean-moderation profile: 1 on drug, decay or 0 off drug |
| $g_t$ | covariance coupling profile (ch08) |
| $b_i$ | standardized biomarker $(B_i - \mu_{bm})/\sigma_{bm}$ |
| $\beta_{bm:D}$ | biomarker-by-treatment interaction (estimand) |
| $\psi_g$ | logistic gating slope (latent-class reading) |

---

## Summary

<!-- evidence-labeled results go here -->

Two architectures encode the biomarker-treatment interaction.

**Mean moderation** applies a post-draw shift $\beta_{bm} b_i
\sigma_{BR} h_t$ to $BR_{it}$; the estimand is $\beta_{bm:D} =
-\beta_{bm}\,\sigma_{BR}/\sigma_{bm}$. The interaction is zero when
the drug is not taken (plain profile $h_t = 0$ off drug), so power
is flat in the carryover half-life.

**Covariance moderation** places the interaction in the MVN
correlation matrix: $\mathrm{Cor}(B_i, BR_{it}) = c_{bm} g_t$;
the estimand is $\beta_{bm:D} = -c_{bm}\,\sigma_{BR}/\sigma_{bm}$.
Conditioning on the biomarker, this is mean moderation with profile
$h_t = g_t$ and a rank-one reduction in the residual covariance
(Paper 01, Discussion; derived). The size of the carryover loss is
therefore set by the off-drug coupling profile $g_t$, not by which
channel carries the interaction.

**The dual-channel architecture** (Paper 11) activates both channels
with independent weights $c_{bm,A}$ (mean) and $c_{bm,B}$
(covariance).

At matched total $N = 70$ and matched effect size the power
difference between mean moderation and covariance moderation under the
refined process is 15.4% for covariance moderation against 6.9% for
mean moderation under carryover (Paper 01 Results; NOTATION.md Part
2). This reflects the graded off-drug coupling of the covariance
channel, which attenuates the interaction signal under carryover.

### Estimands and units by architecture

The interaction estimand is the regression coefficient on the
centered-biomarker-by-$D_{bc}$ term. The formula is the same for the
two single-channel architectures at matched $c_{bm}$.

| Architecture | Interaction parameter | Estimand | Units |
|---|---|---|---|
| Mean moderation | $\beta_{bm}$ | $\beta_{bm:D} = -\beta_{bm}\,\sigma_{BR}/\sigma_{bm}$ | Change in BR per SD of biomarker, per unit $D_{bc}$ |
| Covariance moderation | $c_{bm}$ | $\beta_{bm:D} = -c_{bm}\,\sigma_{BR}/\sigma_{bm}$ | Change in BR per SD of biomarker, per unit $D_{bc}$ |
| Dual-channel | $c_{bm,A}$, $c_{bm,B}$ | Combined; linear in both channels | Change in BR per SD of biomarker, per unit $D_{bc}$ |

When $c_{bm}$ (covariance moderation) is set equal to $\beta_{bm}$
(mean moderation) and the same $\sigma_{BR}$, $\sigma_{bm}$ are used,
the two single-channel architectures produce the same estimand. Their
difference in power arises from the covariance structure of the
data-generating process -- specifically from the off-drug coupling
profile $g_t$ -- not from a difference in the estimand being targeted.
For the dual-channel architecture, the two channel contributions add
on the covariance scale, giving an effective coupling of $c_{bm,A} +
c_{bm,B}$ on drug and a combined estimand that is linear in both
parameters (Paper 11; derived).

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3.3 for this chapter:
- Sources: paper 01; CM "Two architectures", "The dual-channel
  architecture", "A latent-variable reading"; 21 "Biological
  Assumptions"; 23 Sections 1.6, 8.1 to 8.2, 11; 02 Section 4;
  33 Sections 2.3 and 3.2; 36 configurations D and E;
  plan-architecture-c.md; Appendix I Section 5 (calibration).

- New to write: an estimand-and-units table; the current contrast
  (NOTATION Part 2; paper 01 under graded coupling) replacing 10
  to 15% against 40 to 60%; the mean-moderation estimand shift
  (33 Section 2.3).
-->

### Mean moderation

<!-- From CM "Mean moderation", lines 1308-1358 -->

In the mean-moderation architecture (`dgp_architecture =
'mean_moderation'`, default `moderation_scaling = 'constant'`;
`R/generateData.R`, inspected) the package applies the shift

$$
BR_{it} \;\leftarrow\; BR_{it} + \beta_{bm}\, b_i\, \sigma_{BR}\,
  h_t, \qquad h_t = \begin{cases}
  1 & D_{it} = 1 \\
  0 & D_{it} = 0 \text{ (plain)}
\end{cases}
$$

after the multivariate-normal draw. The decayed profile ($h_t =
2^{-t_{sd,t}/t_{1/2}}$ off drug) is Paper 01's comparison arm.
The implied on-drug interaction coefficient on the raw biomarker is

$$
\beta_{bm:D} = -\beta_{bm}\,\frac{\sigma_{BR}}{\sigma_{bm}}
$$

(Paper 01, Estimand). Under the plain profile the moderation stops
when the drug stops, and power is flat in $t_{1/2}$ (0.681 at every
half-life for E9 in Hybrid at $N = 70$, $\beta_{bm} = 0.25$,
compound symmetry; Paper 01, Results). Paper 01 calls this flatness
a coding choice rather than a finding, because the drug effect itself
carries over.

### Covariance moderation

<!-- From CM "Covariance moderation", lines 1360-1409 -->

In the covariance-moderation architecture the biomarker and the
response components are drawn jointly from one multivariate normal
whose biomarker-response correlation depends on exposure state:

$$
\mathrm{Cor}(B_i, BR_{it}) = c_{bm}\, g_t,
$$

with the graded coupling $g_t$ of ch08: 1 on drug,
$e^{-\lambda t_{sd}} = 2^{-t_{sd}/t_{1/2}}$ off drug after exposure,
0 before any exposure. Without carryover the estimand is
$\beta_{bm:D} = -c_{bm}\,\sigma_{BR}/\sigma_{bm}$ (Paper 01,
Estimand).

Conditioning on the biomarker:

$$
E(BR_{it} \mid B_i) = \mu_{BR,t} + c_{bm}\, g_t\, \sigma_{BR}\,
  b_i, \qquad \mathrm{Cov}(BR_{i\cdot} \mid B_i) = \Sigma_{BR} -
  c_{bm}^2 \sigma_{BR}^2\, g g^\top,
$$

so the covariance encoding is mean moderation with exposure profile
$h_t = g_t$ and strength $\beta_{bm} = c_{bm}$, and a drug-response
covariance reduced by a rank-one term (Paper 01, Discussion; derived;
verified to machine precision on Hybrid path 1 under graded coupling
at $t_{1/2} = 1$ week, `check-coupling-closed-form.R`). The
divergence in power between the two architectures reflects the
exposure profile of the interaction, not a difference in kind.

Under the published step rule ($g_t = 1$ wherever the drug-response
mean is nonzero), the coupling gap $1 - \bar g_\text{off}$ is zero
at every positive half-life, and the test collapses to its size
(`docs/36`, Section 5; Paper 01, Results).

### Biological rationale

<!-- From docs/legacy/21 "Biological Assumptions"; docs/36
Section 6.3 -->

Hendrickson et al.'s rationale for covariance moderation: the
biomarker predicts the drug effect only probabilistically. A
high-biomarker participant is likely to have a large response but will
not always; the biomarker enters the model through the correlation
matrix, not the mean.

Mean moderation asserts the opposite: the biomarker fixes each
participant's drug effect deterministically (conditional on the
biomarker, no residual coupling variance remains on drug). This is
the standard precision-medicine model for a predictive biomarker.

From `docs/36` Section 6.3: mean moderation removes the ceiling on
the moderation parameter (no positive-definiteness constraint from the
biomarker row) and addresses both the PD problem and the carryover
sensitivity problem. It changes the model rather than repairing it.

### Configurations D and E

<!-- From docs/36 configurations D and E -->

Configuration D: the published compound-symmetry construct with mean
moderation ($\rho_{cs} = 0.8$, constant cross-factor block, no
covariance coupling, $\beta_{bm}$ shifts the BR mean on drug). No
ceiling. Power is flat across all half-lives in Hybrid at $c_{bm} =
0.30$: 0.78, 0.75, 0.74, 0.77, 0.73 (verified, `docs/36` Table 9).

Configuration E: separable AR(1) with mean moderation and no
covariance coupling ($\rho = 0.7$, separable $c_1 = 0.2$, no
biomarker coupling in the matrix). No ceiling. Power at $c_{bm} =
0.30$: 0.67, 0.65, 0.67, 0.66, 0.65 (verified; note the CO result
differs substantially from D because AR(1) eliminates the
compound-symmetry-driven on/off correlation at 2.5-week spacing).

### The dual-channel architecture

<!-- From CM "The dual-channel architecture", lines 1815-1814 -->

Paper 11 activates the mean and covariance channels simultaneously
with independent weights $c_{bm,A}$ (mean channel) and $c_{bm,B}$
(covariance channel).

**Construction.** The covariance channel is built into the MVN
correlation matrix with $c_{bm,B}$. The mean channel is then applied
as a post-draw additive shift at on-drug timepoints:

$$
\Delta BR_{it} = c_{bm,A}\, b_i\, \sigma_{BR}\,
  \mathbf{1}[D_{it} = 1]
$$

(`dgp_architecture = 'combined'` in `R/generateData.R`, inspected).
The covariance channel carries the graded off-drug coupling; the
mean channel stops when the drug stops.

**The channels add on the covariance scale.** Each channel's on-drug
contribution to $\mathrm{Cov}(B, BR)$ can be computed separately
(covariance channel: $c_{bm,B}\sigma_{bm}\sigma_{BR}$; mean channel:
$c_{bm,A}\sigma_{bm}\sigma_{BR}$), and because the combined-architecture
on-drug response is their sum, the two contributions add, giving
effective coupling $c_{bm,A} + c_{bm,B}$ on drug.

**Boundary reductions.** At $(c_{bm,A} = 0, c_{bm,B} > 0)$ the
architecture reduces to pure covariance moderation; at $(c_{bm,A} >
0, c_{bm,B} = 0)$ it reduces to pure mean moderation. These boundary
cells are needed to confirm that paper 11's DGP reproduces the
single-channel results (open work item; see `CLAUDE.md`).

### A latent-variable reading

<!-- From CM "A latent-variable reading of covariance moderation",
lines 2015-2133 -->

Paper 03 asks what mechanism a covariance-moderation model could
summarize if the true process is a discrete responder/non-responder
dichotomy. With class label $Z_i \mid B_i \sim
\mathrm{Bernoulli}(\pi(B_i))$ and $\pi(B_i) = \mathrm{expit}(\psi_0
+ \psi_g b_i)$, a two-component mixture with class-dependent means
and common within-class covariance produces a marginal
$(B_i, BR_{it})$ distribution whose second moments are reproduced by
a single MVN with $c_{bm}^2 = f_B f_{BR}$ -- the product of the
between-class variance fractions in $B$ and in $BR$ (derived; Paper
01, Discussion).

Under this reading the covariance-moderation architecture is the
second-moment signature of a latent-class mixture. The approximation
breaks down for strongly bimodal responses, a step-function $\pi$, or
class-varying within-class covariance.

### The current power contrast

<!-- From NOTATION.md Part 2; Paper 01 Results -->

At matched total $N = 70$ in the Hybrid design, the relative power
loss under carryover at $t_{1/2} = 0.5$ week is 15.4% for covariance
moderation against 6.9% for mean moderation (Paper 01 Results;
NOTATION.md Part 2). This supersedes the earlier "40-60% versus
8-13%" contrast, which came from a pre-restatement version of the
paper that used unmatched N and the published step-rule coupling
(not graded coupling).

<!-- TODO: the current contrast under graded coupling restated as a
table: architecture, t_1/2, power (no carryover), power (carryover),
relative loss; sourced from Paper 01 production results -->

### The estimand shift under mean moderation

Under mean moderation the package applies the shift

$$
BR_{it} \;\leftarrow\; BR_{it} + \beta_{bm}\, b_i\, \sigma_{BR}\, h_t
$$

after the multivariate-normal draw, where $b_i =
(B_i - \mu_{bm})/\sigma_{bm}$ is the standardized biomarker and $h_t$
is the exposure profile. Under the plain profile ($h_t = 1$ on drug,
$h_t = 0$ off drug), the shift is zero at every off-drug visit. The
regression coefficient on the centered biomarker-by-$D_{bc}$
interaction then recovers $-\beta_{bm}\,\sigma_{BR}/\sigma_{bm}$
(Paper 01, Estimand; derived). This on-drug-only shift means the
moderation signal is present only while the drug is taken and
disappears entirely in the washout period.

Because the shift is zero wherever $h_t = 0$, the off-drug visits
contribute no interaction signal under the plain profile. An on-minus-
off estimand -- which uses information from both on-drug and off-drug
assessments -- therefore loses nothing from off-drug visits: there is
nothing to lose. Power is flat across carryover half-lives, because
the half-life controls only off-drug mean carryover, not the
moderation itself (verified under configuration D, Hybrid, $N = 70$,
$c_{bm} = 0.30$: power 0.78, 0.75, 0.74, 0.77, 0.73 at half-lives
0, 0.1, 0.2, 0.5, 1 week; `docs/36` Table 9). Paper 01 notes this
flatness is a coding choice -- the drug effect carries over, but the
moderation is not coded to carry over with it.

Under the decayed profile ($h_t = 2^{-t_{sd,t}/t_{1/2}}$ off drug),
the shift persists during washout: a high-biomarker participant has a
higher residual response than a low-biomarker participant in the
off-drug period. This makes the moderation profile match the drug
effect profile and gives behavior closer to covariance moderation.
The off-drug contribution to the interaction signal then depends on
$t_{1/2}$, and power falls as $t_{1/2}$ increases and the off-drug
moderation signal approaches zero (`docs/38`, Section 4).

The shift under graded coupling in the covariance architecture is
mechanistically identical: conditioning on $B_i$, the conditional
mean of $BR_{it}$ shifts by $c_{bm}\, g_t\, \sigma_{BR}\, b_i$ with
the graded coupling profile $g_t$. The power difference between the
two architectures at matched $c_{bm}$ is therefore entirely a
consequence of which coupling profile -- plain-stop versus graded
exponential decay -- is applied to the off-drug visits, not of any
difference in what is being estimated.

---

## Derivations

| Statement | Source | Script |
|---|---|---|
| Covariance moderation = mean moderation + rank-one covariance reduction | `eq:archB-conditional` | `check-coupling-closed-form.R` |
| Channels add on the covariance scale: effective coupling $c_{bm,A} + c_{bm,B}$ | CM `eq:dc-covA`, `eq:dc-covB` | derived |
| Latent-class approximation: $c_{bm}^2 = f_B f_{BR}$ | CM `eq:latent-mvn-approx` | derived |
| Estimand under mean moderation: $\beta_{bm:D} = -\beta_{bm}\sigma_{BR}/\sigma_{bm}$ | `eq:archA-kent` | Paper 01 Estimand |

---

## Numerical evidence

From `docs/36` Table 9 (power under configurations D and E; verified):

| Config. | Architecture | $c_{bm}$ | $t_{1/2} = 0$ | 0.1 | 0.2 | 0.5 | 1 |
|---|---|---|---|---|---|---|---|
| D | Mean (CS) | 0.30 | 0.78 | 0.75 | 0.74 | 0.77 | 0.73 |
| E | Mean (sep. AR(1)) | 0.30 | 0.67 | 0.65 | 0.67 | 0.66 | 0.65 |

---

## Relation to the manuscripts

- **Paper 01:** defines both architectures; states the estimands;
  derives the conditional equivalence; reports the power contrast
  under graded coupling.
- **Paper 11:** defines the dual-channel architecture; boundary
  reductions not yet run.
- **Paper 03:** latent-variable reading; pilot results on class-aware
  tests.

---

## Limitations and open items

- The dual-channel boundary cells (pure mean, pure covariance) need
  to be re-run with a driver whose parameters match Paper 01's
  production driver.
- The super-additivity result for the dual channel has not been
  derived analytically.
- The mean-moderation estimand shift under decayed profile has not
  been quantified against the graded-coupling covariance result.

---

## Reproduce

```bash
bash tools/render.sh \
  analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd
```

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | All | Phase 2 content move from Paper 01, CM sec:architectures, CM sec:dualchannel, CM sec:latent, docs/36, NOTATION.md | Whitepaper refactor |
