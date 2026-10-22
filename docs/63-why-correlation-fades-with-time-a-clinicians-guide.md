---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# Compound Symmetry versus Separable AR(1): Why We Changed the Visit-to-Visit Covariance Model {.unlisted .unnumbered}
*2026-10-21 11:30 PDT*

**Author.** pmsimstats team

**Purpose.** This note explains one specific change in how our
simulations model the correlation structure of a patient's repeated
symptom-component scores: replacing a **compound-symmetry** (CS)
correlation matrix with a **separable AR(1)** one, and, cross-factor,
replacing a constant cross-component correlation with one that decays
with the separable AR(1) factor. It is written for a reader with
statistical training (correlation matrices, variance components,
mixed models) but no prior exposure to this project's notation, in
the style of our more code-focused `docs/53`/`docs/55` notes. Unlike
`docs/35-separable-response-covariance.md`, the full technical
derivation this note draws on, this is meant to be read once, stand
alone, and go as far as the Kronecker/separability argument and its
consequence for the feasible biomarker-effect range, without
reproducing every derivation.

```{=latex}
\clearpage
\tableofcontents
\clearpage
```

## Notation

| Symbol | Meaning |
|---|---|
| $w_t$ | calendar week of visit $t$ (Hybrid design: 4, 8, 9, 10, 11, 12, 16, 20) |
| $\rho$ | within-component AR(1) autocorrelation per week (0.7 here) |
| $\rho_{cs}$ | compound-symmetry correlation, constant across all visit pairs (0.8 as published by Hendrickson et al.; 0.7 in the package default) |
| $c_1$ | cross-component correlation at zero lag (same visit), under the separable form (0.2) |
| $c_\times$ | the published process's constant cross-component correlation at every lag, same or different visit (0.1) |
| TV, PB, BR | the three latent response components (natural course, placebo belief, biological drug response) |
| $K_p$ | the $p \times p$ cross-component correlation matrix at a single, shared timepoint ($p=3$ here) |
| $K_t$ | the $n \times n$ within-component correlation matrix across the $n$ visits |
| $\otimes$ | Kronecker product |
| $\Sigma$ | the full $(np) \times (np)$ correlation matrix over every (component, visit) pair |

## 1. Summary

**The question.** How correlated is a patient's score on one
component (say, the drug-response component, BR) at two different
visits, as a function of how many weeks apart those visits are? And
separately: how correlated are *different* components (say TV and
PB) with each other, at the same visit and at different visits?

**Compound symmetry (the published process).** Within a component,
every pair of visits shares one correlation, $\rho_{cs}$, regardless
of lag. Across components, every pair of visits also shares one
constant, $c_\times$, again regardless of lag, except at the same
visit, where it is a different constant $c_1$. Formally, the
within-component block is

$$
K_t^{CS} = (1-\rho_{cs})\,I_n + \rho_{cs}\,J_n,
$$

where $I_n$ is the identity and $J_n$ the all-ones matrix: every
off-diagonal entry equals $\rho_{cs}$, independent of which pair of
visits it connects. This is the classical **equicorrelation** or
**exchangeable** structure: visit order and spacing carry no
information.

**Separable AR(1) (the refined process).** Within a component,
correlation decays geometrically with the gap in weeks:

$$
[K_t]_{t,s} = \rho^{\lvert w_t - w_s \rvert}.
$$

Across components, the correlation at any pair of visits, same or
different, is the cross-component correlation $K_p$ (constant $c_1$
off-diagonal) **multiplied by the same time-decay factor** that
governs the within-component case:

$$
\Sigma_{(c,t),(c',s)} = [K_p]_{c,c'} \times \rho^{\lvert w_t-w_s\rvert}.
$$

This multiplicative structure — one factor for "which pair of
components" and a separate, shared factor for "how far apart in
time" — is what **separable** means here, made precise in Section 3.

**Why it matters.** The compound-symmetry structure imposes a hard
ceiling on how strong the biomarker-response coupling can be before
the joint distribution becomes mathematically infeasible (Section 4).
Under Hendrickson's own published value, that ceiling is 0.256 (on a
correlation scale where 1 is deterministic). Under separable AR(1) at
$\rho=0.7$, it is approximately 0.474-0.480 — nearly double. The
clinical motivation (symptom scores recorded a week apart should
resemble each other more than scores sixteen weeks apart) and the
statistical motivation (a higher, more defensible ceiling; exact
coherence with the fitted analysis model's residual structure) point
the same direction, which is why both were adopted together.

## 2. Compound symmetry and AR(1) as covariance models

Restricting attention to one component for now (say BR, at the $n=8$
Hybrid-design visits), a correlation matrix is a statement about how
*any* two visits for the same patient relate. Two classical structures
dominate the repeated-measures literature:

**Compound symmetry / exchangeable.** $\text{Cor}(Y_t, Y_s) =
\rho_{cs}$ for every $t \ne s$. This is exactly the correlation
structure induced by a **random-intercept-only** model: $Y_{it} = b_i
+ \varepsilon_{it}$, with $b_i \sim N(0,\sigma_b^2)$ a patient-level
random effect and $\varepsilon_{it} \sim N(0,\sigma_e^2)$ i.i.d. Then
$\rho_{cs} = \sigma_b^2 / (\sigma_b^2+\sigma_e^2)$, the intraclass
correlation (ICC). Every pair of visits is equally correlated because
the *only* source of within-patient dependence is a single,
time-invariant patient effect; nothing in the model lets two visits
that happen to be close in time resemble each other any more than two
visits far apart.

**AR(1).** $\text{Cor}(Y_t, Y_s) = \rho^{\lvert w_t - w_s \rvert}$.
This is the correlation induced by a first-order autoregressive
process in continuous time: $Y_{it} = \mu_i + u_{it}$, where $u_{it}$
follows $u_{i,t} = \rho^{\Delta} u_{i,t-1} + \text{noise}$ for a time
gap $\Delta$. Two visits close together inherit most of each other's
deviation from the patient's mean trend; two visits far apart inherit
almost none of it. This is the structure `nlme::corCAR1()` fits (the
continuous-time generalization of discrete-lag AR(1), needed here
because Hybrid-design visits are unevenly spaced).

Neither is more "correct" in the abstract; they encode different
substantive claims about *why* repeated measurements on the same
patient are correlated. Compound symmetry says: there is a stable
patient trait that shifts every visit equally. AR(1) says: deviations
from trend are serially correlated and fade with elapsed time. For a
fluctuating symptom measure recorded over 4 to 20 weeks, the second
claim is the more clinically defensible one — and, not coincidentally,
it is also the residual structure the project's recommended analysis
model (`corCAR1`, random intercept) actually assumes when fitting real
or simulated data. Using CS to *generate* the data while fitting an
AR(1) residual structure to *analyze* it is an avoidable
data-generating-process/analysis-model mismatch; Section 5 returns to
this coherence point.

## 3. What "separable" means, precisely

The full problem is not one component at $n$ visits, but **three**
latent components (TV, PB, BR) each observed at the same $n$ visits,
plus their correlation with a biomarker measured once. Restricting to
the three components, the full correlation structure is an $(np)
\times (np)$ matrix $\Sigma$ indexed by pairs (component, visit). A
**separable** (equivalently, Kronecker-product) covariance structure
is one that factors as

$$
\Sigma = K_p \otimes K_t,
$$

where $K_p$ is $p \times p$ (here $p=3$, the component-to-component
correlation at a single shared timepoint) and $K_t$ is $n \times n$
(the within-component, cross-visit correlation). Written entrywise,
for components $c, c'$ and visits $t, s$:

$$
\Sigma_{(c,t),(c',s)} = [K_p]_{c,c'} \cdot [K_t]_{t,s}.
$$

The substantive claim this encodes is that **the correlation between
any two (component, visit) pairs splits into two independent
factors**: one that depends only on *which pair of components* is
involved, and one that depends only on *which pair of visits* is
involved, with no interaction between the two. Equivalently: the
*ratio* of cross-component to within-component correlation is the
same at every lag. Under the refined process, $K_t$ is the AR(1)
matrix of Section 2 ($[K_t]_{t,s} = \rho^{\lvert w_t-w_s\rvert}$), and
$K_p$ has 1's on the diagonal and $c_1$ off it, so

$$
\Sigma_{(c,t),(c',s)} =
\begin{cases}
\rho^{\lvert w_t-w_s\rvert} & c = c' \\
c_1 \, \rho^{\lvert w_t-w_s\rvert} & c \ne c'.
\end{cases}
$$

Compare this with the published process. There, the within-component
block uses $K_t^{CS}$ (Section 2), and the cross-component block uses
a *constant* $c_\times$ at every lag except zero, where it is $c_1$
instead. That is **not** separable in the sense defined above: the
cross-component correlation does not scale with the *same* decay
function that governs the within-component block (indeed it has no
decay at all once $t \ne s$), so no single $K_p \otimes K_t$
factorization reproduces it. `docs/35-separable-response-covariance.md`
(Section 4) works out the published structure explicitly as a sum of
two terms, $K_p^{\text{person}} \otimes J_n + K_p^{\text{occasion}}
\otimes (I_n - \tfrac{1}{n}J_n)$ in one parameterization — a valid
correlation structure, but a strictly more complicated, non-separable
one, with person-level and occasion-level cross-component coupling
ratios that differ (0.143 and 0.333 at $\rho_{cs}=0.7$; *quoted*,
`docs/35`).

**Why adopt the simpler, separable form.** Four reasons, elaborated in
Section 5: it is valid (positive definite) for *every* visit schedule,
not only ones checked case by case; it needs one cross-component
parameter instead of two; it is the standard model for multivariate
repeated measures, with an established literature on estimation and on
testing separability itself; and, with constant component standard
deviations, it makes the *summed* observed score (the sum of the three
latent components, which is what gets analyzed) exactly a
random-intercept-plus-AR(1) process — the same structure the
recommended analysis model fits.

## 4. A worked numeric comparison

The Hybrid design's eight visits fall at weeks 4, 8, 9, 10, 11, 12, 16
and 20. The table compares, within one component, what each
correlation model predicts for a patient's own score at representative
visit pairs.

| Visit pair | Gap (weeks) | CS: $\rho_{cs}$ | AR(1): $\rho^{\text{gap}}$, $\rho=0.7$ |
|---|---|---|---|
| Week 11 vs. week 12 | 1 | 0.80 | 0.700 |
| Week 4 vs. week 8 | 4 | 0.80 | 0.240 |
| Week 4 vs. week 12 | 8 | 0.80 | 0.058 |
| Week 4 vs. week 20 | 16 | 0.80 | 0.0033 |

Table: Within-component correlation under each structure (*derived*)

Every CS entry is identical by construction; AR(1) shrinks by a factor
of $\rho$ for each additional week of gap. At the largest gap in the
design, CS overstates the correlation by a factor of roughly 240
relative to AR(1).

For the **cross-component** comparison, take BR and PB at the same
visit (lag 0) and at an 8-week gap. Separable AR(1) gives $c_1
\rho^0 = 0.2$ at lag 0 and $c_1\rho^8 = 0.2 \times 0.7^8 = 0.2 \times
0.0576 = 0.0115$ at an 8-week gap. The published process gives $c_1 =
0.2$ at lag 0 (same as separable AR(1), since both use $c_1$ at zero
lag) but a flat $c_\times = 0.1$ at *every* positive lag, including
the 8-week one — 0.1 against separable AR(1)'s 0.0115, an
almost-ninefold overstatement at that gap, and with no further decay
at longer gaps either.

## 5. Why we want separability, in more depth

**Validity for every schedule, in closed form.** A compound-symmetric
or separable correlation matrix's eigenvalues have closed forms that
do not depend on how many visits there are or how they are spaced
(`docs/35`, Section 4.2-4.3): for the separable form, the matrix is
positive definite whenever $K_p$ and $K_t$ individually are, which
holds for any $\rho \in (-1,1)$ and any $c_1$ inside an explicit,
$n$-independent range. The published process's non-separable
structure instead requires checking positive-definiteness
numerically, design by design, and (as `docs/55`/`docs/61` document)
the shipped code's response to a violation was to silently repair the
matrix rather than reject the input.

**An exact, checkable ceiling.** Section 6 of `docs/35` derives the
largest biomarker-response coupling $c_{bm}^{\ast}$ the joint
distribution can support, as a closed-form function of $\rho$ (or
$\rho_{cs}$) and the visit schedule — not a number that must be
rediscovered by trial and error for each new design. Numerically:
0.256 under the published compound-symmetry structure; 0.440 at
$\rho_{cs}=0.8$ under separable AR(1); 0.474 (0.480 in the Hybrid
design specifically) at $\rho=0.7$ (*quoted*,
`docs/series/ch10-refined-process.md`).

**Fewer parameters, standard structure.** One cross-component constant
($c_1$) instead of two ($c_1, c_\times$); Kronecker-product
covariances are an established model for multivariate repeated
measures, with known estimation theory and known tests of the
separability assumption itself (`docs/35`, Section 5, items 4 and 9).

**Coherence with the analysis model.** With the three components'
standard deviations held constant across visits, the *sum* of the
three latent components (what the analysis model actually sees as the
outcome) is, under the separable form, exactly a random intercept plus
an AR(1) residual process — precisely `corCAR1()`. Under the
non-separable published structure, the summed outcome instead carries
an extra white-noise term the analysis model does not represent
(`docs/35`, Section 4.6). Generating data under one correlation model
while analyzing it as though a different one held is an
avoidable mismatch; adopting the separable form removes it.

**Faster computation.** Inverting and Cholesky-factorizing a
Kronecker-product matrix reduces to factorizing the $p \times p$ and
$n \times n$ pieces separately, rather than the full $np \times np$
matrix.

## 6. What this costs

**Old results do not reproduce exactly.** Because the lagged
cross-component correlations change value (not just the within-
component ones), no simulation result computed under the published
structure reproduces exactly under the separable one; this is
intentional, not an error (`docs/35`, Section 5, item 1).

**One decay rate for every component.** The separable form uses a
single $\rho$ for TV, PB and BR alike. It cannot represent, for
instance, a placebo-belief effect that fades on its own, faster or
slower timescale than the drug-response component does. A valid
extension (a sum of several separable terms, one per rate) exists but
gives up the closed-form ceiling (`docs/35`, Section 8).

**No same-visit excess correlation.** Under the separable form, the
ratio of cross-component to within-component correlation is the fixed
constant $c_1$ at *every* lag, including zero; a same-visit "shared
shock" stronger than this ratio implies cannot be represented without
an explicit additional term.

**Separability is itself a testable, not a proven, assumption.**
Whether the three latent components actually satisfy it cannot be
checked directly on real data, since only their sum (the observed
symptom score) is measured; the components themselves are latent.
Adopting separable AR(1) is a stated modeling commitment, consistent
with standard multivariate repeated-measures practice, not an
empirically validated fact about this particular condition
(`docs/35`, Section 5, item 5).

## 7. Where this leaves the simulations

Separable AR(1) (Kronecker-product correlation, one shared time-decay
factor across all three components) replaces the published process's
non-separable compound-symmetry structure, as part of the refined
process paper 01 recommends. This is one of several bundled changes;
`docs/53` covers the carryover-decay fix and `docs/55` the
biomarker-coupling fix, both independent of the covariance-structure
choice discussed here.

## References

- `analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd`, Methods
  and "Room for larger effect sizes" sections, for $\rho_{cs}$ as
  published (0.8) and as implemented (0.7), and the ceiling table.
- `docs/35-separable-response-covariance.md`, the full technical
  derivation: the Kronecker/separability argument (Section 4), the
  advantages-and-disadvantages discussion (Section 5), the ceiling
  calculation (Section 6), and the per-component decay-rate
  alternative (Section 8).
- `docs/series/ch10-refined-process.md`, for the refined process's
  adoption of separable AR(1) and the quoted ceiling figures.
- `docs/53-carryover-recursion-in-published-code.md` and
  `docs/55-step-coupling-in-published-code.md`, companion notes on
  the other two defects in the published code.
- `docs/61-cross-factor-correlation-bug-paper02-driver.md`, on a
  later implementation bug in one driver's cross-component
  correlation parameters, unrelated to the structural choice discussed
  here but concerning the same parameters.
