---
geometry: margin=2.4cm
fontsize: 10pt
---

# A Separable Response Covariance for the Covariance-Moderation DGP {.unlisted .unnumbered}
*2026-09-30 11:12 PDT; expanded 2026-10-01 11:46 PDT*

**Author:** pmsimstats team

**Revision 2026-10-07.**

- **Added:**
  - an intuitive overview ("In brief");
  - the general ceiling with prognostic and baseline couplings
    (Section 4.8);
  - the power cost of the higher ceiling (Section 6);
  - the simulation evidence now available (Section 10).
- **Corrected:**
  - Hendrickson's published $\rho$ is 0.8; the figures quoted for
    `orig` were at the package value 0.7. Both are now given.
  - The coupling pattern is written $g_t$, not $\phi_t$, which was also
    used for the AR(1) step correlation.
  - "The summed response is AR(1)" holds when the component SDs are
    constant over visits.
  - The analysis model that the separable structure matches includes a
    random intercept.
  - The `orig` ceiling depends only on visit count and on-share for
    binary coupling.

**Purpose.** The project is converging on a single covariance-based
data-generating process (DGP) to replace the two now in use, the
package's Architecture B (`covar`) and the published Hendrickson et al.
implementation [1] at commit `58b32a9` (`orig`). An AR(1) within-factor
kernel has been chosen. This note sets out one further structural
choice, that the correlation among the three response factors across
visits be *separable*, and examines it in detail: what it asserts, why
it is coherent where the current construct is not, how it relates to
Hendrickson's construct, what it implies for the largest testable
biomarker correlation, and what it would change across the compendium.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## In brief: what separable AR(1) means

*Added 2026-10-07: an intuitive description, ahead of the formal
treatment.*

**The setting.** Each participant has three hidden responses that move
over the trial:

- the drug response (BR);
- the placebo, or expectancy, response (PB);
- the natural course (TV).

At any visit each one is a little above or below its expected
trajectory. The question is how those deviations relate across visits,
and across the three responses.

**"AR(1)": memory that fades with time.** If a participant's natural
course is running better than expected this week, it will probably
still be better than expected next week, less so in a month, and hardly
at all by the end of the trial. With $\rho = 0.7$ per week, the
correlation between two visits of the same response is:

| Weeks apart | 1 | 2 | 4 | 16 |
|---|---|---|---|---|
| Correlation, $0.7^{\text{weeks}}$ | 0.70 | 0.49 | 0.24 | 0.003 |

Table: AR(1) correlation between visits by weeks apart at $\rho = 0.7$

Each week's deviation is mostly the previous one carried forward and
shrunk, plus something new. Nearby visits look alike; distant visits
are nearly independent.

**"Separable": one clock for everything.** The three responses are
linked to one another at a fixed strength, a correlation of
$c_1 = 0.2$ at the same visit. Separability says that this link fades
with time on exactly the same clock as each response's own memory. The
correlation between two responses at two visits is therefore the
product of two pieces:

$$
\mathrm{Cor}(\text{response } c \text{ at week } t,\ \text{response } c' \text{ at week } s)
= \underbrace{K_{cc'}}_{\text{which responses}} \times
  \underbrace{\rho^{|t - s|}}_{\text{how far apart}} .
$$

For example, placebo at week 8 and natural course at week 12 are
correlated $0.2 \times 0.7^4 \approx 0.05$. The "which responses" part
and the "how far apart" part never interact.

**The generative picture.** Imagine three independent random processes,
each drifting slowly with the same 0.7-per-week memory, and make each
response a fixed blend of the three. Every correlation the model claims
is then produced by something. That is why the correlation matrix can
never be invalid, whatever the visit schedule (Section 4.2).

**How it compares with the alternatives.**

- **Hendrickson et al.'s construct (`orig`, compound symmetry).** Each
  person carries a permanent offset for each response, plus independent
  noise at each visit. Any two visits of a response are correlated at
  RH's $\rho = 0.8$, whether one week or twenty weeks apart. The
  construct is valid, but its memory never fades. Sections 3 and 4.5
  quote its implied quantities at the package value $\rho = 0.7$.
- **The package's current AR(1) form (`covar`, not separable).** It has
  fading memory, but it also asserts an extra same-visit link between
  responses that no underlying process generates. On closely spaced
  visits it asks for correlations no valid covariance can deliver, and
  the matrix fails (Section 4.4).
- **Separable AR(1) (proposed).** Fading memory on one clock, with
  nothing extra, so it is always valid.

**Why choose it.**

- **It is valid for any visit schedule.**
- **The largest testable biomarker effect has an exact formula.** Each
  on/off switch between nearby visits has a cost, so schedules with
  many switches can carry only smaller effects (Section 4.3).
- **It matches the analysis model.** The outcome also contains each
  person's permanent baseline level. The data therefore look like a
  random intercept plus fading AR(1) noise, which is what the
  random-intercept-plus-`corCAR1` analysis assumes. The match is close
  but not exact, because the placebo response's SD scales with
  expectancy (Section 4.6).

**What it gives up.**

- **All three responses forget at the same rate.** A faster-fading
  placebo response, for example, is not allowed (Section 8 gives the
  valid extension).
- **The responses carry almost no long-term memory.** Week 4 and week
  20 are essentially unrelated within a response (correlation 0.003).
  Lasting person-to-person differences come only from the baseline
  level.

## 1. Summary

1. **The proposal.** Replace the response block of the correlation
   matrix with a Kronecker product $M = K \otimes A$: a $3 \times 3$
   factor correlation $K$ times a single AR(1) kernel $A$ in calendar
   time. Lagged cross-factor correlation becomes $c_1\rho^{\text{gap}}$,
   and the parameter `c.cfct` is retired.
2. **All three constructs are sums of separable terms.** Hendrickson's
   `orig` is exactly a person-level term plus an occasion-level term,
   each a valid covariance, and is positive definite for every schedule.
   The current `covar` is an AR(1) term plus an occasion-level term
   whose factor matrix has zero diagonal and is therefore indefinite.
   That term carries cross-factor correlation without carrying any
   variance, which is why `covar` fails on dense schedules. The proposal
   is a single valid term (Section 4).
3. **`covar` has an exact failure condition.** With a common $\rho$, it
   is positive definite if and only if
   $\lambda_{\min}(A)(1 - c_\times) > c_1 - c_\times$. At $\rho = 0.7$
   this fails for half-week spacing at the production values
   $(c_1, c_\times) = (0.2, 0.1)$ and paper 04's $(0.3, 0.2)$, and for
   quarter-week spacing at $(0.1, 0.05)$ (verified, Section 4.4).
4. **The ceiling on $c_{bm}$ becomes exact and interpretable.** Under
   separability it factorizes into a factor term and a schedule term,
   and the schedule term is a sum of per-step costs: on-to-on
   $(1-\phi)/(1+\phi)$, on-to-off $\phi^2/(1-\phi^2)$, off-to-on
   $1/(1-\phi^2)$, with $\phi = \rho^{\text{gap}}$ (verified,
   Section 4.3).
5. **The change reaches most of the compendium.** The same
   non-separable construct is implemented four times: in the package's
   `buildSigma()`, in `implementations/tidyverse`, in
   `implementations/nof1power`, and (in Hendrickson's form) in paper
   04's private generator. The response block is shared by both DGP
   architectures, so mean-moderation results change too. Papers 01-03,
   05-11 and 13-14 would need re-runs (Section 7).
6. **Separability matches the analysis model.** Under separability the
   summed response component is an AR(1) process when the component SDs
   are constant over visits. With the person-level baseline in the
   outcome, the data are then a random intercept plus AR(1), which is
   what the random-intercept-plus-`corCAR1` analysis assumes. Under
   `covar` they also carry a white-noise term the analysis model does
   not include (Section 4.6, derived).

## 2. Review of the original note

The two-page version of 30 September holds up on its central claims,
each of which was rechecked numerically for this revision: the
positive-definiteness guarantee, the factorized ceiling, the factor term
$\sqrt{1 - R^2}$ with $R^2 = 0.067$ at $c_1 = 0.2$, and the ceiling table.
It had five gaps, addressed below.

- **It did not explain why `covar` fails.** It reported the half-week
  failure without a mechanism. Section 4.4 derives the exact condition.
- **It did not contrast Hendrickson's construct.** `orig` turns out to
  be a coherent two-level model, and the comparison clarifies what was
  lost and gained in moving to `covar` (Sections 3 and 4.5).
- **It understated the scope of the change.** It described an edit to
  `buildSigma()`. The construct is duplicated across three other
  engines, affects both architectures, and invalidates the parity
  baselines (Section 7).
- **It did not note a consequence of the AR(1)-only kernel.** Within a
  factor, correlation at a 16-week gap is $0.7^{16} = 0.003$ under
  AR(1), against a constant $\rho$ under `orig` (0.8 as published; 0.7
  at the package value). Person-level persistence
  in the outcome is then carried only by the baseline term
  $\mathrm{BL}_i$ (Section 3).
- **Its main reference was vague.** Section 10 now gives the sources,
  with their verification status.

## 3. Intuition

Each participant's symptom score is $Y_{it} = \mathrm{BL}_i - (TV_{it} +
PB_{it} + BR_{it})$, and the three components are drawn jointly. The
question is how a deviation in one component at one visit relates to
deviations in the components at other visits.

**Hendrickson's construct: traits plus shocks.** `orig` behaves as if
each participant carried a fixed offset for each component that lasts
the whole trial (a trait), plus an independent weekly disturbance (a
shock). Because traits never fade, any two visits of a factor are
correlated at $\rho$, whether one week or twenty weeks apart: 0.8 in the
published construct, 0.7 at the package value. The traits are weakly
correlated across components and the visit-level shocks more strongly:
0.125 and 0.500 at the published values, 0.143 and 0.333 at
$\rho = 0.7$ (verified). This is the classic
random-intercept-plus-noise picture, and it is internally consistent.

**`covar`: memory that fades, with a leftover.** `covar` replaced the
everlasting trait with fading memory: today's deviation carries into
next week at 0.7, into the week after at 0.49, and so on. That is
clinically more plausible. But it kept Hendrickson's two cross-factor
numbers. The same-visit cross-factor correlation of 0.2 is
Hendrickson's shock-level sharing, and in `orig` it is carried by the
shocks' own variance. `covar` has no separate shock variance, so the
same-visit excess (0.2 against 0.1 at other visits) is asserted without
any variance component to produce it. On widely spaced visits this is
harmless. On closely spaced visits, where adjacent measurements of a
factor are nearly identical, the matrix is asked to make two nearly
identical quantities correlate differently with a third, which no valid
covariance can do.

**The separable proposal: one clock.** The proposal lets all variation
flow through a single fading memory shared by the three components. A
shared disturbance today affects all three components, and its
influence fades at the same rate for each component's own history and
for their association with each other. Generatively, the three
components are fixed blends of three independent AR(1) processes that
tick on the same clock. Nothing is asserted that the clock does not
generate, which is why validity cannot fail.

**What the AR(1)-only choice gives up.** Under fading memory alone the
components carry almost no person-level persistence across a 20-week
trial. The outcome still has person-level persistence through
$\mathrm{BL}_i$, which is a participant-level term in $Y_{it}$, so the
analysis model's random intercept remains supported. If persistent
component-level traits are wanted, Section 8 gives the valid extension.

## 4. Derivations

### 4.1 Setup and the feasibility condition

For one participant on one path with $n$ visits, the DGP draws
$X = (B, \mathrm{BL}, TV_{1:n}, PB_{1:n}, BR_{1:n})$ from a multivariate
normal distribution with correlation matrix

$$
R = \begin{pmatrix} 1 & 0 & \tilde r^{\top} \\ 0 & 1 & 0^{\top} \\
\tilde r & 0 & M \end{pmatrix},
\qquad \tilde r = c_{bm}\,(e_3 \otimes u),
$$

where $M$ is the $3n \times 3n$ response block, $e_3 = (0, 0, 1)^\top$
selects $BR$, and $u \in \mathbb{R}^n$ is the coupling pattern ($u_t = 1$
on drug and $u_t = g_t$ off drug, with $g_t$ the residual fraction of
the on-drug effect; $g_t = 0$ before any exposure). Variances scale $R$ to the covariance and do not
affect positive definiteness. By the Schur complement on the biomarker
row, $R$ is positive definite if and only if

$$
M \succ 0 \quad\text{and}\quad 1 - \tilde r^{\top} M^{-1} \tilde r > 0,
\qquad\text{so}\qquad
c_{bm}^{\ast} = \bigl((e_3 \otimes u)^{\top} M^{-1} (e_3 \otimes u)\bigr)^{-1/2}.
$$

The quantity $\tilde r^{\top} M^{-1}\tilde r$ is the $R^2$ of the
biomarker regressed on all response variables, so the condition states
only that the biomarker cannot be more than perfectly predictable.

### 4.2 The separable structure

Let $A_{ts} = \rho^{|w_t - w_s|}$ and
$K = (1 - c_1) I_3 + c_1 J_3$, and set $M = K \otimes A$, so that
$\mathrm{Cor}(X_{c,t}, X_{c',s}) = K_{cc'} A_{ts}$ with factor-major
ordering. Three standard Kronecker identities [15] carry the argument:

$$
(K \otimes A)(k \otimes a) = (\kappa k) \otimes (\lambda a), \qquad
(K \otimes A)^{-1} = K^{-1} \otimes A^{-1}, \qquad
\det(K \otimes A) = (\det K)^n (\det A)^3 .
$$

**Validity.** The eigenvalues of $M$ are the products $\kappa_i\lambda_j$
of those of $K$ and $A$. $K$ has eigenvalues $1 + 2c_1$ (once) and
$1 - c_1$ (twice), so $K \succ 0$ for $-\tfrac12 < c_1 < 1$. $A$ is the
correlation matrix of a stationary Ornstein-Uhlenbeck process sampled at
the visit times; the exponential kernel $e^{-\theta|\Delta|}$ is positive
definite for any distinct times [3, 17]. Hence $M \succ 0$ for every
schedule.

**Generative form.** Writing the response as an $n \times 3$ matrix
$\mathbf{X}$ with one column per factor, separability is the
matrix-normal model $\mathrm{Cov}(\mathrm{vec}\,\mathbf{X}) = K \otimes A$
[8], realized as $\mathbf{X} = L_A \mathbf{Z} L_K^{\top}$ with
$L_A L_A^{\top} = A$, $L_K L_K^{\top} = K$ and $\mathbf{Z}$ standard
normal. Each column of $L_A\mathbf{Z}$ is an independent AR(1) path; the
factors are their fixed mixtures through $L_K$.

### 4.3 The factorized ceiling and its switch-cost form

Substituting $M^{-1} = K^{-1} \otimes A^{-1}$ into Section 4.1:

$$
(e_3 \otimes u)^{\top}(K^{-1} \otimes A^{-1})(e_3 \otimes u)
= [K^{-1}]_{33}\; u^{\top} A^{-1} u,
\qquad
c_{bm}^{\ast} = \bigl([K^{-1}]_{33}\bigr)^{-1/2}\bigl(u^{\top}A^{-1}u\bigr)^{-1/2}.
$$

**Factor term.** For the equicorrelated $K$,
$K^{-1} = (1-c_1)^{-1}\bigl(I_3 - \tfrac{c_1}{1 + 2c_1} J_3\bigr)$, so

$$
[K^{-1}]_{33} = \frac{1 + c_1}{(1 - c_1)(1 + 2c_1)} = \frac{1}{1 - R^2},
\qquad R^2 = \frac{2c_1^2}{1 + c_1},
$$

where $R^2$ is the squared multiple correlation of $BR$ on $TV$ and $PB$
at one visit: 0.0182, 0.0667 and 0.1385 at $c_1 = 0.1, 0.2, 0.3$
(verified).

**Schedule term.** A sampled Ornstein-Uhlenbeck process is Markov:
$X_t = \phi_t X_{t-1} + \sqrt{1-\phi_t^2}\,\varepsilon_t$ with
$\phi_t = \rho^{d_t}$ and $d_t = w_t - w_{t-1}$ [3]. Factorizing the
joint density gives $A^{-1} = L^{\top}L$ with $L$ lower bidiagonal,
$L_{11} = 1$, $L_{tt} = (1-\phi_t^2)^{-1/2}$ and
$L_{t,t-1} = -\phi_t(1-\phi_t^2)^{-1/2}$, so

$$
u^{\top} A^{-1} u = u_1^2 + \sum_{t=2}^{n}
\frac{(u_t - \phi_t u_{t-1})^2}{1 - \phi_t^2}.
$$

Each term is the squared, standardized one-step prediction error of the
coupling pattern under the AR(1) predictor. For a binary pattern the
steps cost:

| Step | Cost | Short gap ($\phi \to 1$) |
|---|---|---|
| on to on | $(1-\phi)/(1+\phi)$ | small |
| on to off | $\phi^2/(1-\phi^2)$ | large |
| off to on | $1/(1-\phi^2)$ | large |
| off to off | $0$ | none |

Table: Cost of each coupling-pattern step under the AR(1) predictor

The ceiling is therefore
$c_{bm}^{\ast\,-2} = [K^{-1}]_{33}\{u_1 + \sum \text{step costs}\}$: a
coupling pattern that changes abruptly between closely spaced visits is
expensive, and one that is constant across closely spaced visits is
cheap. The formula agrees with the full matrix to $9 \times 10^{-16}$
on the paper 01 binding paths and a 16-visit multi-cycle schedule
(verified). It also corrects the approximation in
`docs/34-cbm-feasible-range.md` Section 5.5, which assigned every switch
the off-to-on cost; an on-to-off switch costs exactly one unit less. A
pattern that starts on drug pays $u_1^2 = 1$ up front, so the two
orderings of a single switch cost the same in total.

For a constant pattern ($u = 1$), $u^{\top}A^{-1}u = 1 + \sum_t
(1-\phi_t)/(1+\phi_t) = n_{\text{eff}}$, the effective number of
independent visits, and the ceiling is
$\sqrt{(1 - R^2)/n_{\text{eff}}}$.

### 4.4 Why the current `covar` fails

With a common $\rho$, `covar` sets within-factor entries to $A_{ts}$,
same-visit cross-factor entries to $c_1$, and different-visit
cross-factor entries to $c_\times A_{ts}$. Hence (verified, exact
reconstruction):

$$
M_{\text{covar}} = K_a \otimes A + K_b \otimes I_n, \qquad
K_a = (1 - c_\times) I_3 + c_\times J_3, \qquad
K_b = (c_1 - c_\times)(J_3 - I_3).
$$

$K_b$ has zero diagonal and positive off-diagonal entries, so it is
indefinite, with eigenvalues $2(c_1 - c_\times)$ and $-(c_1 - c_\times)$
(twice). It is the same-visit excess correlation, attached to no
variance. Diagonalizing $A = Q\Lambda Q^{\top}$,

$$
(I_3 \otimes Q)^{\top} M_{\text{covar}} (I_3 \otimes Q)
= \bigoplus_{j=1}^{n} \bigl(\lambda_j K_a + K_b\bigr),
$$

and since $K_a$ and $K_b$ share eigenvectors, the eigenvalues of
$M_{\text{covar}}$ are $\lambda_j(1 + 2c_\times) + 2(c_1 - c_\times)$
and $\lambda_j(1 - c_\times) - (c_1 - c_\times)$. Therefore

$$
M_{\text{covar}} \succ 0 \iff
\lambda_{\min}(A) > \frac{c_1 - c_\times}{1 - c_\times}.
$$

The predicted smallest eigenvalue matches the computed one to
$10^{-15}$ on six schedules and three parameter sets (verified). For
equally spaced visits with $\phi = \rho^d$, the eigenvalues of the AR(1)
Toeplitz matrix lie above the minimum of its spectral density,
$(1 - \phi)/(1 + \phi)$, and approach it as $n$ grows [16]. A sufficient
condition for validity at any $n$ is then
$(1-\phi)/(1+\phi) \geq (c_1 - c_\times)/(1 - c_\times)$, which at
$\rho = 0.7$ requires gaps of at least 0.63 weeks for $(0.2, 0.1)$, 0.71
weeks for $(0.3, 0.2)$, and 0.30 weeks for $(0.1, 0.05)$. Every
production design examined satisfies this; dense or long schedules may
not.

### 4.5 Hendrickson's construct as two valid terms

`orig` sets within-factor entries to $\rho$ at every lag and
cross-factor entries to $c_1$ (same visit) and $c_\times$ (otherwise).
Hence (verified, exact reconstruction):

$$
M_{\text{orig}} = K_o \otimes I_n + K_p \otimes J_n, \qquad
K_p = \rho I_3 + c_\times(J_3 - I_3), \qquad
K_o = (1 - \rho) I_3 + (c_1 - c_\times)(J_3 - I_3).
$$

This is the covariance of $X_{c,t} = P_c + E_{c,t}$, with person-level
effects $P \sim N(0, K_p)$ and occasion-level effects $E_t \sim N(0, K_o)$
independent across visits [2]. The implied cross-factor correlations are
$c_\times/\rho$ at the person level and
$(c_1 - c_\times)/(1 - \rho)$ at the occasion level: 0.143 and 0.333 at
the package value $\rho = 0.7$, and 0.125 and 0.500 at the published
$\rho = 0.8$ (verified). Because $J_n$ has eigenvalues $n$ and $0$,
$M_{\text{orig}} \succ 0$ if and only if $K_o \succ 0$ and
$K_o + nK_p \succ 0$; with $K_p \succeq 0$ this reduces to
$\lambda_{\min}(K_o) = (1 - \rho) - (c_1 - c_\times) > 0$, independent of
$n$ and the schedule. It is 0.20 at $\rho = 0.7$ (verified for $n = 4$
to 32) and 0.10 at the published $\rho = 0.8$ (verified). This recovers the expression given in paper 01
Appendix B.7.

`orig` is a single separable term only when $K_o \propto K_p$, that is
when $c_\times = c_1\rho$, in which case
$M_{\text{orig}} = K \otimes A_{\text{CS}}$ (verified: the difference is
0.040 at $c_\times = 0.1$ and 0 at $c_\times = 0.14$).

### 4.6 The summed response and the analysis model

The analysis model regresses $Y_{it}$, which contains the sum of the
three components, with a random intercept and `corCAR1` residuals. Let
$s_t = (s_{TV,t}, s_{PB,t}, s_{BR,t})^{\top}$ hold the factor standard
deviations at visit $t$ and $S_t = \sum_c s_{c,t} X_{c,t}$ the summed
component. Under the separable structure,

$$
\mathrm{Cov}(S_t, S_s) = (s_t^{\top} K s_s)\, A_{ts},
$$

which is exactly proportional to the AR(1) kernel whenever the factor
standard deviations are constant over visits. In this construct the PB
SD is $10e_t$, so it differs between open-label and blinded visits. The
correlation factor
$s_t^\top K s_s / \sqrt{(s_t^\top K s_t)(s_s^\top K s_s)}$ between an
open-label and a blinded visit is 0.976 (verified), so the departure
from AR(1) is small.

The outcome is $Y_{it} = \mathrm{BL}_i - S_{it}$, and $\mathrm{BL}_i$ is
a person-level constant independent of the responses. Hence
$\mathrm{Cov}(Y_{it}, Y_{is}) = \sigma_{BL}^2 + \mathrm{Cov}(S_t, S_s)$:
a random intercept plus an AR(1) process. That is exactly the residual
structure of the random-intercept-plus-`corCAR1` analysis model, not of
`corCAR1` alone. Under `covar`, by the
decomposition of Section 4.4,

$$
\mathrm{Cov}(S_t, S_s) = (s_t^{\top} K_a s_s)\, A_{ts}
+ (s_t^{\top} K_b s_t)\,\mathbf{1}[t = s],
$$

an AR(1) process plus a white-noise term, since $s^{\top} K_b s =
2(c_1 - c_\times)(s_{TV}s_{PB} + s_{TV}s_{BR} + s_{PB}s_{BR}) > 0$. A
`corCAR1` residual model does not include such a term, so under `covar`
the analysis model's residual structure is misspecified even before the
biomarker enters, whereas under separability it is correctly specified
up to the time-varying placebo-belief scale (which the package ties to
the expectancy weight). Both statements follow algebraically from the
verified decompositions; their effect on test calibration has not been
measured.

### 4.7 One framework for all three

All three constructs are instances of the linear model of
coregionalization, $M = \sum_j K_j \otimes A_j$, in which each term is a
factor matrix times a correlation kernel [13, 14]. A sum of such terms
is guaranteed valid when every $K_j$ is positive semidefinite.

| | `orig` (58b32a9) | `covar` (current) | Proposed |
|---|---|---|---|
| Terms | $K_p \otimes J + K_o \otimes I$ | $K_a \otimes A_{\text{AR}} + K_b \otimes I$ | $K \otimes A_{\text{AR}}$ |
| Factor matrices | both PSD | $K_b$ indefinite | PSD |
| Interpretation | person traits + occasion shocks | fading memory + unbacked same-visit excess | fading memory shared by all factors |
| Within-factor correlation at 16 weeks | $\rho$ (0.8 published, 0.7 package) | 0.003 | 0.003 |
| Valid for every schedule | yes | no (Section 4.4) | yes |
| Ceiling depends on | visit count and on-share only (binary coupling) | full schedule | full schedule, in closed form |
| Cross-factor parameters | $c_1$, $c_\times$ | $c_1$, $c_\times$ | $c_1$ |

Table: Comparison of the `orig`, `covar` and proposed response covariance constructs

### 4.8 The ceiling with prognostic and baseline couplings

The stress test of `docs/46` couples the biomarker to the natural
course, the placebo response or the baseline level as well as to the
drug response. Write the biomarker's correlations with the three
response factors as vectors $v_{TV}$, $v_{PB}$, $v_{BR}$ over the
visits, and its correlation with $\mathrm{BL}$ as $c_{bl}$. Since
$\mathrm{BL}$ is uncorrelated with the responses, the Schur complement
of Section 4.1 gives

$$
R \succ 0 \iff M \succ 0 \quad\text{and}\quad
c_{bl}^2 + \sum_{c, c'} [K^{-1}]_{cc'}\; v_c^\top A^{-1} v_{c'} < 1 ,
$$

using $M^{-1} = K^{-1} \otimes A^{-1}$. The quadratic form agrees with
the full matrix to 10 digits on the Hybrid design (verified,
`analysis/scripts/quick-sim/cbm-ceiling/06-separable-review-checks.R`). With $c_{bm} = c_{bm,TV} = c_{bl} = 0.3$ on Hybrid path A,
the left side is 0.75, and the stress-test cells are feasible.

The off-diagonal elements of $K^{-1}$ are negative,
$-c_1/((1 - c_1)(1 + 2c_1))$. Couplings of the same sign to correlated
factors therefore cost less than the sum of their separate costs. On
Hybrid path A with $c_1 = 0.2$, coupling 0.3 to BR alone costs 0.390,
0.3 to TV alone costs 0.342, and both together cost 0.663: the cross
term is $-0.070$ (verified). A biomarker that is both predictive and
prognostic in the same direction is accordingly easier to accommodate
than either coupling alone would suggest.

## 5. Advantages and disadvantages

**Advantages.**

1. **Validity for every schedule.** No spacing, length or irregularity
   can make the response block indefinite (Section 4.2).
2. **No implementation ambiguity.** There is no per-pair cross-factor
   decay rate, so the loop-order dependence in `buildSigma()` and its
   copies cannot arise.
3. **An exact, interpretable ceiling** that can be checked for every
   simulation cell before data are drawn (Section 4.3).
4. **A recognized structure.** Kronecker-product covariances are the
   standard model for multivariate repeated measures [6, 7]; estimation
   [8] and tests of separability [9, 10] are established.
5. **Fewer parameters.** One cross-factor parameter instead of two.
6. **Faster computation.** Inversion and Cholesky factorization reduce
   to the $3 \times 3$ and $n \times n$ factors.
7. **Coherence with the analysis model.** With constant component SDs,
   the outcome is a random intercept plus AR(1), the structure of the
   random-intercept-plus-`corCAR1` analysis. Under `covar` it carries
   an additional white-noise term (Section 4.6).

**Disadvantages.**

1. **Results change.** Lagged cross-factor correlations double at
   production values ($c_1 \rho^{\text{gap}}$ instead of
   $c_\times \rho^{\text{gap}}$ with $c_\times = c_1/2$). No current
   result reproduces exactly.
2. **One persistence for all factors.** Factor-specific $\rho$ (for
   example faster-fading placebo belief) is excluded. The valid
   extension is a sum of separable terms (Section 8), at the cost of the
   closed-form ceiling.
3. **No same-visit excess.** The ratio of cross-factor to within-factor
   correlation is $c_1$ at every lag, so cross-factor association cannot
   be stronger at the same visit than the within-factor pattern implies.
   If occasion-level shared shocks are substantively important, a valid
   occasion-level term must be added explicitly (Section 8).
4. **No lead-lag structure.** One factor cannot predict another at a
   later visit more strongly than the reverse; neither current construct
   allows this either.
5. **Separability is a testable assumption** [9, 10], and it is
   frequently rejected in applied multivariate and space-time data,
   which has motivated nonseparable alternatives [11, 12]. Whether
   symptom-component data would reject it is unknown, since the
   components are latent. For a simulation DGP this is a modeling
   commitment to state, not a defect.

## 6. Effect on the largest testable $c_{bm}$

Ceilings without carryover at $\rho = 0.7$, $c_1 = 0.2$ (current
structure also $c_\times = 0.1$); closed form and full matrix agree to
$10^{-16}$ (verified):

| Path | Current `covar` | Separable AR(1) |
|---|---|---|
| CO, path A | 0.619 | 0.616 |
| Hybrid, path A | 0.456 | 0.480 |
| Hybrid, path C | 0.465 | 0.491 |
| OL+BDC, path A | 0.452 | 0.474 |
| Multi-cycle, 16 weekly visits | 0.262 | 0.286 |
| Dense, 0.5-week spacing | not positive definite | 0.509 |

Table: Ceilings without carryover under the current and separable AR(1) structures by path

The design conclusions of `docs/34` carry over unchanged in direction,
because they are driven by the schedule term, which separability makes
exact rather than approximate. Their numbers would need recomputing.

**A higher ceiling is not more power.** At the same $c_{bm}$, AR(1)
data carry much less information about the interaction than
compound-symmetric data. The person-level persistence that compound
symmetry provides cancels in a within-participant contrast; fading
memory does not. In the Hybrid design at $c_{bm} = 0.25$, the
closed-form power of the paired-difference statistic is 0.721 under
compound symmetry and 0.335 under the separable structure (`docs/37`,
Section 5.8). Effect sizes under the two structures should therefore
be compared on the scale of the implied interaction slope, not of
$c_{bm}$.

**Graded coupling and the step costs.** Under graded coupling the
off-drug pattern is not binary: $u_t = g_t = 2^{-t_{sd,t}/t_{1/2}}$.
The general term $(u_t - \phi_t u_{t-1})^2/(1 - \phi_t^2)$ of Section 4.3
then replaces the binary on-to-off cost $\phi_t^2/(1 - \phi_t^2)$ with
$(g_t - \phi_t)^2/(1 - \phi_t^2)$, which is smaller when $g_t$ is close
to $\phi_t$. This is why the ceiling under graded coupling rises with
the half-life.

## 7. Impact across the compendium

The response block is shared by both architectures: mean moderation
draws $TV$, $PB$ and $BR$ from the same matrix and adds its shift
afterwards. The change therefore affects every paper whose data come
from one of the engines below, regardless of architecture.

| Engine | Construct | Parameters |
|---|---|---|
| `R/generateData.R` (`buildSigma`) | `covar` | set by drivers |
| `implementations/tidyverse` | `covar`, same loops | set by drivers |
| `implementations/nof1power` | `covar`, same loops | set by drivers |
| Paper 04 private generator (`treatment-main-effect/vig*.R`) | `orig`-like (CS within, constant cross) with the 2024 scale factor | $(0.3, 0.2)$ |

Table: Simulation engines with their response covariance construct and parameters

Per paper, from the drivers (inspected; production designs were not all
re-checked against the validity condition of Section 4.4):

| Paper | Engine | $(c_1, c_\times)$ | Consequence |
|---|---|---|---|
| 01 DGP architectures | package | (0.2, 0.1) | Section 3 re-run already required; Appendix B.3, B.5.2, B.6-B.7 and the ceiling discussion rewritten around one construct |
| 02 Carryover sensitivity | tidyverse | (0.2, 0.1) | Full re-run; decay-shape study should also decay the correlation with the chosen shape |
| 03 Latent class | package | (0.1, 0.05) | Re-run of mixture and comparator DGPs |
| 04 Main effect | private, `orig`-like | (0.3, 0.2) | Unaffected unless migrated; decide whether paper 04 should adopt the unified DGP |
| 05 Design sensitivity | nof1power | (0.1, 0.05) | Re-run; design sweeps over period length should be checked against Section 4.4, since short periods are where `covar` fails |
| 06 Component decomposition | package | (0.1, 0.05) | Re-run; its contaminated-biomarker variant couples $B$ to $PB$ as well, which generalizes the ceiling to $\sum_{cc'}[K^{-1}]_{cc'}\,u_c^{\top}A^{-1}u_{c'}$ |
| 07 Gompertz evaluation | package | (0.2, 0.1) | Re-run |
| 08 Test procedure and design | package | (0.1, 0.05) | Re-run; cycle and period sweeps to be checked against Section 4.4 |
| 09 Informative dropout | package | (0.1, 0.05) | Re-run; its reproduction of Hendrickson's Figure 4A moves further from `orig` |
| 10 Test calibration | tidyverse (via paper 02 core) | (0.2, 0.1) | Re-run; standard-error calibration depends directly on the residual covariance, and part of the `corCAR1` miscalibration it reports may come from the white-noise term of Section 4.6 (not tested) |
| 11 Combined architecture | package | driver-set | Re-run |
| 12 Design efficiency | closed form, single-factor AR(1) | none | Unaffected; its closed forms extend exactly under separability by the factor term |
| 13, 14 Two-component | package, `components` | (0.2, 0.1) | Re-run; with two factors $[K^{-1}]_{22} = 1/(1 - c_1^2)$ |

Table: Engine, cross-factor parameters and consequence of the change for each paper

Two infrastructure items also change. The cross-implementation parity
baselines (`analysis/scripts/parity/`) would need regenerating, and the
Hendrickson comparison arm should be rebased on the `58b32a9` sources,
since the vendored copy now in the repository is the 2024 state.

## 8. Alternatives considered

- **Keep `covar` and fix the loop order.** Retains the failure region of
  Section 4.4 and the unbacked same-visit excess.
- **Separable AR(1) (proposed).** One valid term, one clock.
- **Separable AR(1) plus an occasion-level term,**
  $K_s \otimes A_{\text{AR}} + K_o \otimes I$ with $K_o \succeq 0$. This
  legitimately restores a same-visit excess, at the price of a
  within-factor nugget (within-factor correlation below 1 at lag
  zero-plus) and the loss of the closed-form ceiling.
- **Person-level term plus AR(1),** $K_p \otimes J + K_s \otimes
  A_{\text{AR}}$. Restores Hendrickson's persistent traits within a
  fading-memory model, the decomposition advocated by Diggle [2, 4].
  The AR(1)-only kernel was chosen instead; this remains the natural
  route if component-level persistence is later wanted.
- **Separable compound symmetry,** $K \otimes A_{\text{CS}}$:
  Hendrickson's construct with $c_\times = c_1\rho$. Valid and simple,
  but without decay with lag.

## 9. Recommendation and implementation

Adopt the separable AR(1) structure for the single covariance-moderation
DGP, and use the same response block for mean moderation. Implementation
steps:

1. In `buildSigma()`, replace the within- and cross-factor loops with
   `kronecker(K, A)` over the retained components, keeping factor-major
   labels; retire `c.cfct` (error if supplied, rather than ignore).
2. Add `cbmCeiling()` implementing Section 4.3, and call it before
   simulation in `generateSimulatedResults()` and
   `validateParameterGrid()`, stopping when $c_{bm} \geq c_{bm}^{\ast}$.
   Remove the silent positive-definiteness repair for this construct.
3. Apply the same change to `implementations/tidyverse` and
   `implementations/nof1power`, or retire them in favor of the package.
4. Tests: positive definiteness at quarter-week spacing; closed-form
   ceiling equal to the numerical ceiling; exact equality with the
   current construct when $c_1 = c_\times = 0$, where the two coincide;
   two-component case.
5. Regenerate parity baselines and re-run the affected papers in the
   order of Section 7.

## 10. Evidence status and limitations

- **Verified by computation** (`analysis/scripts/quick-sim/cbm-ceiling/
  05-separable-derivations.R`, run 2026-10-01): the decompositions of
  `covar` and `orig` (exact reconstruction), the validity condition for
  `covar` (eigenvalues to $10^{-15}$), the closed form of
  $[K^{-1}]_{33}$, the switch-cost form of $u^{\top}A^{-1}u$, the
  separability condition for `orig`, and the ceilings of Section 6.
- **Derived, not separately computed:** the equal-spacing gap
  thresholds of Section 4.4, which rest on the Toeplitz eigenvalue
  bound [16] and are consistent with the computed eigenvalues; and the
  covariance of the summed response in Section 4.6, which follows from
  the verified decompositions. Its effect on paper 10's calibration
  results is a hypothesis.
- **Inspected:** the engines and parameter values of Section 7, from
  the driver sources. Paper 04's production engine was not traced in
  full, and the production designs of papers 05 and 08 were not checked
  against the validity condition.
- **Simulated since 2026-10-01.** The separable structure at
  $\rho = 0.7$ (configuration C of `docs/36`) has been used in three
  studies, with no positive-definiteness failure at $c_{bm}$ up to
  0.45:
  - a power arm of `docs/36`;
  - the base construct of the covariance study (`docs/45`, 378,000
    fits);
  - the base construct of the pre-registered stress test (`docs/46`,
    26,500 trials, plus prognostic and baseline couplings).

  Under these data:
  - the random-intercept-plus-`corCAR1` analysis with CR2 standard
    errors was nominal everywhere (pooled null rejection 0.048 to
    0.053);
  - its model-based version with phase-specific variances was nominal
    too (0.047 to 0.052);
  - a compound-symmetry working model was anticonservative in the
    crossover design, with $\kappa = 0.74$ derived in `docs/45`,
    Section 4.3.

  This is consistent with Section 4.6, but it does not test the
  white-noise hypothesis about `covar`.
- **Adopted** as the reference structure in paper 01 (Methods, "Why
  separability"; Appendix A.5).
- **Not implemented in the package.** Steps 1 to 3 of Section 9 are
  pending. The package's sampler still repairs a matrix silently when
  its Cholesky factorization fails (`R/generateData.R`, lines 390 to
  392). The structure exists in the simulation drivers only
  (`04-power-simulation.R` arm C; `08-covariance-study.R`;
  `11-strawman-stress-test.R`).
- **References 4 to 17 are still to be verified** before use in a
  manuscript. Paper 01 now relies on the results that cite them.

## 11. References

Status: [1] from the project bibliography; [2] and [3] verified in
PubMed (PMIDs 3233259 and 2049497); the remainder cited from
bibliographic knowledge and to be verified before use in a manuscript.

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020;2.
2. Diggle PJ. An approach to the analysis of repeated measurements.
   *Biometrics* 1988;44(4):959-971.
3. Jones RH, Boadi-Boateng F. Unequally spaced longitudinal data with
   AR(1) serial correlation. *Biometrics* 1991;47(1):161-175.
4. Diggle PJ, Heagerty P, Liang K-Y, Zeger SL. *Analysis of
   Longitudinal Data*, 2nd ed. Oxford University Press; 2002.
5. Pinheiro JC, Bates DM. *Mixed-Effects Models in S and S-PLUS*.
   Springer; 2000.
6. Galecki AT. General class of covariance structures for two or more
   repeated factors in longitudinal data analysis. *Communications in
   Statistics: Theory and Methods* 1994;23(11):3105-3119.
7. Naik DN, Rao SS. Analysis of multivariate repeated measures data with
   a Kronecker product structured covariance matrix. *Journal of Applied
   Statistics* 2001;28(1):91-105.
8. Dutilleul P. The MLE algorithm for the matrix normal distribution.
   *Journal of Statistical Computation and Simulation*
   1999;64(2):105-123.
9. Lu N, Zimmerman DL. The likelihood ratio test for a separable
   covariance matrix. *Statistics and Probability Letters*
   2005;73(4):449-457.
10. Mitchell MW, Genton MG, Gumpertz ML. A likelihood ratio test for
    separability of covariances. *Journal of Multivariate Analysis*
    2006;97(5):1025-1043.
11. Genton MG. Separable approximations of space-time covariance
    matrices. *Environmetrics* 2007;18(7):681-695.
12. Gneiting T. Nonseparable, stationary covariance functions for
    space-time data. *Journal of the American Statistical Association*
    2002;97(458):590-600.
13. Goulard M, Voltz M. Linear coregionalization model: tools for
    estimation and choice of cross-variogram matrix. *Mathematical
    Geology* 1992;24(3):269-286.
14. Wackernagel H. *Multivariate Geostatistics*, 3rd ed. Springer; 2003.
15. Horn RA, Johnson CR. *Topics in Matrix Analysis*. Cambridge
    University Press; 1991.
16. Grenander U, Szego G. *Toeplitz Forms and Their Applications*.
    University of California Press; 1958.
17. Rasmussen CE, Williams CKI. *Gaussian Processes for Machine
    Learning*. MIT Press; 2006.
