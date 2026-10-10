# A critical evaluation of the Gompertz model in pmsimstats-ng {.unlisted .unnumbered}

*2026-05-07 06:11 PDT*

**Author.** pmsimstats team

**Purpose.** This note evaluates the modified Gompertz function that
pmsimstats-ng uses as the mean trajectory of its three latent response
components. It sets out where the form is defensible, which
assumptions it imposes without defense, how those assumptions bear on
the biomarker-treatment interaction test, and which alternative
families could serve as sensitivity analyses. The argument is from the
literature and from inspection of the package source; no simulation
was run for this note. Several of its recommendations were later
taken up by paper 07 (`analysis/report/07-gompertz-evaluation/`), and
dated notes in Sections 6 to 8 point to that work.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. Symbols marked *local*
are defined here. Code identifiers appear in backticks.

| Symbol | Meaning | Status |
|---|---|---|
| $t$ | time in weeks (the Gompertz argument) | canonical |
| $y(t)$ | value of a response-curve family at time $t$ | local |
| $a$, $g_b$, $g_c$ | asymptote, displacement and rate of the textbook Gompertz form (code `maxr`, `disp`, `rate`) | local |
| $\epsilon_{\mathrm{mach}}$, $\operatorname{expm1}(x)$ | double-precision machine epsilon; $e^x - 1$ evaluated without cancellation | local |
| $E_{\max}$, $EC_{50}$, $h$ | maximum effect, half-maximal time and Hill exponent of the Emax/Hill family | local |
| $BR_{it}$, $PB_{it}$, $TV_{it}$ | biological, placebo-belief and time-variant response components | canonical |
| $Y_{it}$ | observed symptom score (code `Sx`) | canonical |
| $B_i$, $b_i$ | biomarker and standardized biomarker (code `bm`, `bm_z`) | canonical |
| $D_{bc,it}$ | exposure-decayed drug indicator (code `Dbc`) | canonical |
| $\beta_{bm:D}$ | biomarker-by-treatment interaction coefficient (`bm:Dbc`) | canonical |
| $\beta_D$ | treatment main effect | canonical |
| $\beta_0$, $\beta_B$, $\beta_t$ | intercept, biomarker main effect and linear time slope of the analysis model | local |
| $u_i$, $\varepsilon_{it}$ | participant random intercept; residual | canonical |
| $c_{bm}$, $\beta_{bm}$ | covariance-moderation correlation; mean-moderation multiplier | canonical |
| $\sigma_{BR}$ | standard deviation of the $BR$ component | canonical |
| $t_{1/2}$, $t_{od}$ | carryover half-life; time on drug (code `tod`) | canonical |
| $N$ | total sample size across randomization paths | canonical |
| $k$ | Weibull shape parameter | canonical |
| $\lambda_Z$, $\mu_m$ | lag time and maximum specific growth rate of Zwietering's modified Gompertz | local |
| $a_i$, $a_0$, $\tau_a^2$ | participant-level asymptote, its mean and variance (Recommendation 5) | local |

Table: Notation index and glossary of symbols used in this document

- **Modified Gompertz (pmsimstats sense).** The Gompertz curve shifted
  so that $y(0) = 0$ and rescaled so that $y(\infty) = a$ (code
  `modgompertz`).
- **Covariance-moderation architecture.** The interaction lives in the
  on-drug $B_i$-$BR$ correlation $c_{bm}$ (formerly Architecture B).
- **Mean-moderation architecture.** The interaction is an on-drug shift
  $\beta_{bm} b_i \sigma_{BR}$ of the $BR$ mean (formerly
  Architecture A).
- **Inflection.** The time at which the curve's slope is greatest.

## 1. Motivation and scope

The pmsimstats-ng simulation framework uses a modified Gompertz
function as the parametric template for three latent response
components: the biological response (BR), the time-varying
non-drug component (TV), and the placebo-belief component (PB).
The function appears in `R/utilities.R` and in
`implementations/{original,original-extended,tidyverse}/R/`
collections; it is invoked by the simulation core at every cell
of every parameter grid. Because the assumed functional form of
the response trajectory is loaded into the simulated outcome
before any analysis-model question is asked, the choice of family
is not a peripheral implementation detail but a foundational
modeling decision. This document evaluates that decision against
the relevant biological and biostatistical literature, articulates
where the Gompertz form is well-justified for the present setting,
identifies where it imposes assumptions the framework does not
explicitly defend, and outlines the alternative families that
could be used in sensitivity analyses.

## 2. The Gompertz function

The Gompertz function in its original form is a double-exponential
sigmoid

$$
y(t) = a \cdot \exp\bigl(-g_b \cdot \exp(-g_c \cdot t)\bigr),
$$

where $a > 0$ is an asymptote, $g_b > 0$ controls displacement of
the inflection along the time axis, and $g_c > 0$ is a rate
parameter [1]. The function was introduced as a demographic model
of human mortality and was first re-purposed as a growth curve by
Winsor [2], who observed that the sigmoidal shape of biological
size-by-time data could be captured parsimoniously by the same
three-parameter form. Subsequent work established the Gompertz
model as a standard description of tumor growth [3, 4, 5],
bacterial population dynamics [6], and innovation diffusion. A
review by Tjørve and Tjørve [7] catalogs parameterizations and
emphasizes the non-orthogonality of $g_b$ and $g_c$, which
complicates identification when both are estimated from data.

In pharmacodynamic modeling, the dominant family is the Hill or
Emax model, $y = E_{\max} \cdot t^{h} / (EC_{50}^{h} + t^{h})$
[8, 9, 10], which differs from the Gompertz in two respects
relevant here. First, the Hill model is symmetric about its
half-maximal point on the logarithmic time scale (it is a logistic
function of $\ln t$); on the linear time scale it is not symmetric,
and for $h > 1$ its inflection lies at
$EC_{50}\,[(h-1)/(h+1)]^{1/h}$, before $EC_{50}$ (derived). Second,
its half-maximal point $EC_{50}$ has a direct pharmacological
interpretation (half-maximal effect concentration, or time when the
argument is time). The Gompertz inflection occurs at
$t = \ln(g_b)/g_c$ and has no comparable interpretive anchor in the
drug-response literature.

## 3. The modified form used in pmsimstats-ng

`R/utilities.R` implements

```r
y <- maxr * exp(-disp * exp(-rate * t))
vert_offset <- maxr * exp(-disp * exp(-rate * 0))
y <- y - vert_offset
y <- y * maxr / (maxr - vert_offset)
```

This is the Gompertz function shifted vertically so that
$y(0) = 0$ and rescaled so that $y(\infty) \to$ `maxr`. The shift
plus rescale maps the family $\{$ Gompertz with arbitrary $y(0)$
and asymptote `maxr` $\}$ onto the family
$\{$ Origin-passing curves with asymptote `maxr` $\}$. The
adjustment is convenient for symptom-improvement modeling
because clinical outcome scales (CAPS, PCL, HAM-D, BDI) are
typically expressed as change-from-baseline so that 0 is the
correct value at study entry.

Two cautions about this implementation. First, the rescaling
factor `maxr / (maxr - vert_offset)`, which equals
$1/(1 - e^{-g_b})$ because `vert_offset` is $a\,e^{-g_b}$,
grows as $1/g_b$ as `disp` shrinks because `vert_offset` approaches
`maxr`, and the evaluation loses accuracy there (derived). (Corrected
2026-10-08: this sentence previously located the instability in the
rescaling factor itself. The factor is benign in exact arithmetic,
and the curve has the finite limit $a(1 - e^{-g_c t})$ as
$g_b \to 0$. The loss of accuracy is floating-point cancellation in
the numerator $\exp(-g_b e^{-g_c t}) - e^{-g_b}$, which subtracts two
numbers close to 1; the relative error is of order
$\epsilon_{\mathrm{mach}}/g_b$, about $6 \times 10^{-9}$ at
$g_b = 10^{-8}$ and $3 \times 10^{-4}$ at $g_b = 10^{-14}$ (at
$g_c = 0.42$, $t = 5$), and $g_b = 0$ returns `NaN`. It is negligible
for $g_b \gtrsim 10^{-6}$. Writing the curve as
$a\,[\operatorname{expm1}(-g_b e^{-g_c t}) -
\operatorname{expm1}(-g_b)] / [-\operatorname{expm1}(-g_b)]$ is
algebraically identical and accurate to machine precision (verified,
`analysis/scripts/compendium-mathematics/check-gompertz-rmanova-cycles.R`).
See `analysis/report/compendium-mathematics.tex`, Section "Shipped
parameters, landmarks, the displacement, and numerical conditioning"
(`sec:invc-gomp-landmarks`).) The current code guards
against the degenerate `maxr == 0` case but not against a
degenerate `disp` close to zero (inspected, `R/utilities.R`,
2026-10-07). Second, the term 'modified Gompertz' in the
microbiology literature [6] denotes a different
reparameterization, in which a lag time $\lambda_Z$ and a maximum
specific growth rate $\mu_m$ are explicit parameters; the
pmsimstats-ng usage is unrelated to that established
modification, and overloading the name in code comments and
manuscripts invites confusion.

## 4. How Hendrickson et al. (2020) use it

Hendrickson et al. [15] use the modified Gompertz to encode three
components of the simulated symptom trajectory. The biological
response (BR) is a function of time on drug $t_{od}$; the
time-varying non-drug component (TV) is a function of calendar
time since study entry; the placebo-belief component (PB) is a
function of cumulative time under positive expectancy, scaled by
the expectancy weight. Each component receives its own
$(\text{maxr}, \text{disp}, \text{rate})$ triple. (In all current
pmsimstats-ng DGPs $t_{od}$ restarts at 0 on re-exposure; it is
time since the current exposure began, not lifetime cumulative
exposure.) The published calibration, the package dataset
`extracted_rp`, is

| Component | maxr | disp | rate (per week) | sd |
|---|---|---|---|---|
| BR | 10.99 | 5 | 0.42 | 8 |
| TV |  6.51 | 5 | 0.35 | 10 |
| PB |  6.51 | 5 | 0.35 | 10 |

Table: Published (Hendrickson) modified Gompertz parameters for the BR, TV, and PB components, from the package dataset `extracted_rp`

The carryover-sensitivity driver
(`default_resp_param()` in
`analysis/scripts/carryover-sensitivity/simulation-core.R`) departs
from this: it gives TV the BR curve (10.99, 5, 0.42) and uses
standard deviations of 6, 2 and 8 for TV, PB and BR. Both sets were
read from the source on 2026-10-07 (inspected); an earlier version
of this table listed the driver values as the Hendrickson defaults.

The published curves rise from zero at study entry and pass through
inflection near $t = \ln(5)/0.42 \approx 3.8$ weeks for BR and
$\ln(5)/0.35 \approx 4.6$ weeks for TV and PB. At week 12 the BR
curve has reached about 97 percent of its asymptote (derived). The
calibration target is the prazosin-versus-placebo CAPS trajectory
of the Raskind trials [16, 17, 18]. The BR-versus-TV-versus-PB
separation in the Hendrickson model is an architectural choice
rather than a mechanism testable from the published data: the
Raskind trials report the total CAPS difference, not a
decomposition into biological-response and non-drug components.
The package documentation of `extracted_rp` says as much: it calls
the set a 'tabula rasa' in which the time course and maximum of
the expectancy (PB) and time-varying (TV) factors were assumed
equal. This identifiability constraint is acknowledged in the
package documentation but is not, to our knowledge, evaluated
empirically in the published literature.

## 5. Strengths in the symptom-trajectory setting

**Sigmoidal shape.** Drug-induced symptom improvement in chronic
neurobehavioral conditions is widely observed to be sigmoidal
rather than linear: there is a delay before measurable benefit,
followed by a phase of more rapid improvement, followed by a
plateau as the underlying pathophysiology stabilizes [10].
A monotone three-parameter sigmoid is therefore a defensible
parametric template for the mean structure of such trajectories,
and the Gompertz is one of several families that can produce this
shape.

**Parsimony.** Three parameters per component is the minimum
needed to represent magnitude, timing, and rate, which is
attractive for a simulation framework that must enumerate many
parameter grids without combinatorial explosion. Richer families
(four-parameter logistic, Weibull cumulative, splines) introduce
identifiability problems when the simulation is calibrated to a
single observed trajectory.

**Computational tractability.** The Gompertz is closed-form,
infinitely differentiable, and free of the behavior of the Hill
function at $t = 0$, where the slope is unbounded for $h < 1$ and
zero for $h > 1$ (derived). The double-exponential form is
numerically stable across the parameter ranges used in the
framework.

**Established precedent in adjacent fields.** The Gompertz has
been used for decades in tumor-growth modeling [3, 4, 5] and
bacterial-growth modeling [6], giving
the family a level of biological face validity that less commonly
used sigmoids do not enjoy.

## 6. Limitations and concerns

**Asymmetry that may not match the data.** The Gompertz is more
strongly skewed than its biological cousins (logistic, Hill,
Weibull cumulative). Its inflection occurs at $1/e$ of the
asymptote (approximately 37 percent), rather than at 50 percent
as in the logistic; for the origin-shifted pmsimstats form with
$g_b = 5$ the fraction is $(e^{-1} - e^{-5})/(1 - e^{-5}) \approx
0.36$ (derived). Whether prazosin-induced CAPS improvement in
PTSD truly accelerates earlier than it decelerates is an
empirical question; the published Raskind trajectories [19] are
not granular enough to settle the question of skewness, and the
Hendrickson framework does not report sensitivity to the choice
of sigmoid family. (*Note, 2026-10-07:* paper 07 has since run
that sensitivity; see the note in Section 6b.)

**Monotonicity.** The Gompertz cannot accommodate response
followed by relapse, which is documented in PTSD pharmacotherapy
[20, 21] and chronic
depression. Patients who improve early in treatment and then lose
benefit at week 8-10 cannot be simulated under any
$(\text{maxr}, \text{disp}, \text{rate})$ choice; the curve is
monotone non-decreasing by construction. For a simulation
framework whose downstream analyses include the OL+BDC design
(open-label titration followed by blinded discontinuation), the
inability to simulate within-treatment relapse is a meaningful
restriction. The pmsimstats-ng framework currently models
relapse only as off-drug carryover decay, which conflates two
distinct biological mechanisms (treatment escape during continued
dosing, versus pharmacological washout after discontinuation).

**Plateau is fixed.** The asymptote `maxr`, like `disp` and
`rate`, is a population-level parameter; there is no
participant-level random effect on any Gompertz parameter.
Participant heterogeneity enters instead as additive deviations
around the common mean curve: each component is drawn from a
multivariate normal with a constant per-timepoint standard
deviation (`respparam$sd`) and AR(1) serial correlation, and under
the mean-moderation architecture the on-drug $BR$ mean is shifted
by $\beta_{bm}\,b_i\,\sigma_{BR}$ (inspected, `buildSigma()` and
`generateData()` in `R/generateData.R`, 2026-10-07). Heterogeneity
of treatment effect, the phenomenon that motivates
predictive-biomarker work, is therefore grafted on at the
architecture layer (the covariance-moderation correlation $c_{bm}$
or the mean-moderation shift), not as a property of the Gompertz
family. Every participant shares the same curve shape and timing;
individuals differ in level, not in rate or ceiling. Under the
default `moderation_scaling = 'constant'` the biomarker shift has
the same size early in exposure, when the mean $BR$ is near zero,
as at the plateau; the `'trajectory'` option scales it with the
mean $BR$ curve instead. That a single shape with level-only
heterogeneity fits individual prazosin responders is a strong
assumption that is not, to our knowledge, defended empirically.
(*Correction, 2026-10-07:* an earlier version said participants
approach a common ceiling at different rates; the source shows no
rate heterogeneity.)

**Non-orthogonal parameters.** The displacement parameter $g_b$
(`disp`) and the rate parameter $g_c$ (`rate`) are strongly
correlated in the Gompertz likelihood surface [7]. In a
calibration setting this is a
nuisance; in a simulation setting this manifests as
parameter-grid axes that are not really independent. A grid that
varies `disp` while holding `rate` fixed is also implicitly
varying the inflection time and the early-trajectory slope, which
may be undesirable for sensitivity analysis. The pmsimstats-ng
parameter-sensitivity sweeps in
`analysis/scripts/carryover-sensitivity/` and
`analysis/scripts/nof1-design-sensitivity/` do not currently
isolate these geometric features of the Gompertz curve.

**Vertical-offset rescaling can amplify noise.** The
$y - y(0)$ shift plus rescale produces a curve whose value at any
given $t$ is sensitive to small changes in `disp` because both
`vert_offset` and the rescale factor depend on it. The function
is well-behaved in the parameter regime used by Hendrickson et
al. [15] but the framework does not currently flag the regions
where the rescale factor diverges.

**Departure from established pharmacodynamic practice.** The
dominant pharmacodynamic family is the Emax/Hill model [8, 9, 10].
Reviewers familiar with the PK-PD tradition may ask why a Gompertz
was preferred to an Emax model whose half-maximal-time parameter
has a direct interpretation. The pmsimstats-ng default
calibration does not currently report a sensitivity analysis under
an Emax alternative, and the published Hendrickson paper does not
defend the Gompertz choice over a Hill or Emax form on biological
grounds.

**Misspecification relative to the analysis model.** The
analysis model in pmsimstats-ng is a linear mixed-effects model
(`nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc, ...)`) with a continuous
exposure-decayed predictor $D_{bc,it}$ (`Dbc`). Under the Gompertz
DGP, the true mean structure is a sum of three sigmoidal
trajectories; the analysis fits this with a linear-in-time term
plus a binary-or-decayed-exposure term. The mismatch is
intentional, and the simulation studies the consequences for power
and Type I error, but the mismatch does mean that the simulated
power estimates are conditioned on a specific functional form of
the DGP that is itself a modeling choice rather than an
established fact.

## 6b. Implications for the biomarker-treatment interaction test

The founding pmsimstats question is the test of the
biomarker-treatment interaction in the mixed-effects model
`Sx ~ bm + t + Dbc + bm:Dbc` with a random intercept per
participant, that is
$Y_{it} = \beta_0 + \beta_B B_i + \beta_t t + \beta_D D_{bc,it} +
\beta_{bm:D} B_i D_{bc,it} + u_i + \varepsilon_{it}$, with
$H_0: \beta_{bm:D} = 0$ [15]. The Gompertz-induced response
trajectory enters the simulated data through the BR component
mean. The limitations enumerated in Section 6 do not contribute
equally to the validity of this interaction test. Ranking them
by their impact on the bias, power, and coverage of
$\hat\beta_{bm:D}$ (argued, not simulated):

**Most critical for the interaction test.**

1. **Common shape, level-only heterogeneity (no random effect
   on `maxr`).** The covariance-moderation architecture encodes
   the interaction signal in the second-moment $B_i$-$BR$
   correlation; the mean-moderation architecture encodes it in
   a first-moment mean shift. In both, between-participant
   heterogeneity in the biomarker-driven treatment effect is a
   difference in *level* around a common curve: no participant
   reaches a different ceiling through a different `maxr`, and
   none responds faster or slower (Section 6). The
   biomarker-stratification literature typically conceptualizes
   responder-versus-non-responder heterogeneity as different
   asymptotic benefit [23, 24]. A level shift of constant size is
   close to that at the plateau but not early in exposure, so
   the restriction matters most for designs whose on-drug
   observations fall on the rising part of the curve. The
   direction of any resulting bias in simulated power is not
   established; a plateau-heterogeneous DGP would be needed to
   measure it. A reviewer who accepts the latent-subtype framing
   of the covariance-moderation architecture will reasonably ask
   why the response ceiling is treated as a population constant.

2. **Monotonicity (no within-treatment relapse).** The
   Gompertz cannot generate response-then-relapse trajectories.
   In real PTSD pharmacotherapy data [20, 21] a non-trivial
   fraction of participants show early benefit followed by
   decay, often correlated with comorbid depression and chronic
   pain status. If the true symptom trajectory has this shape
   and the simulation does not, the linear-in-time analysis
   model is being asked to fit a smoother surface than it would
   encounter in practice. The interaction test is then
   evaluated under conditions kinder than the application; this
   plausibly biases the simulated power upward (optimistic).
   Because the direction of the point (1) effect is not
   established, the two cannot be assumed to offset.

3. **Misspecification relative to the linear analysis model.**
   The analysis model fits `bm + t + Dbc + bm:Dbc`. The true BR
   mean is sigmoidal, not linear in $t$. The linear-in-$t$ term
   absorbs trajectory curvature into a single slope, and any
   residual curvature is pushed into the random intercept and
   the AR(1) residual structure. The interaction term
   `bm:Dbc` is the residual quantity that survives this
   absorption; it is therefore the term most exposed to the
   choice of sigmoid family. A Hill-shaped DGP (more symmetric
   inflection) and a Gompertz-shaped DGP (asymmetric inflection
   at 1/e of the asymptote) push different amounts of curvature
   into the linear approximation and may yield different
   apparent `bm:Dbc` standard errors at matched calibration
   targets. The framework currently does not report sensitivity
   under an alternative sigmoid family; this is the single most
   actionable Gompertz-related contribution to interaction-test
   validity.

   *Note, 2026-10-07.* Paper 07
   (`analysis/report/07-gompertz-evaluation/`) has since run this
   sensitivity: a 16-cell factorial (four families, both
   architectures, two effect sizes; $N = 70$, OL+BDC,
   $t_{1/2} = 0$, 1,000 replicates per cell). As reported in its
   abstract (inspected, not re-run here), the family is exactly
   inert under covariance moderation, where matched seeds give
   identical rejection decisions; under mean moderation with
   trajectory-scaled moderation the four powers span 0.743 to
   0.782, with the Gompertz lowest. The 144-cell extension with
   carryover remains open.

**Moderate impact on the interaction test.**

4. **Asymmetric inflection.** The Gompertz inflection at
   $t = \ln(g_b)/g_c$ corresponds to $1/e \approx 0.37$ of the
   asymptote, in contrast to the logistic $0.5$. Under matched
   `maxr` and matched inflection time, the Gompertz therefore
   places less of its rise before the inflection and more after
   it than the logistic does. Under the mean-moderation
   architecture with trajectory-scaled moderation, the
   interaction signal follows the $BR$ curve, so the family
   changes the mix of contributions across timepoints and the
   absolute power numbers, though plausibly not the qualitative
   conclusions (argued). Under the covariance-moderation
   architecture the interaction lives in the covariance and is
   mean-independent, so the family does not affect it; paper 07
   confirms this (see the note under point 3). (*Correction,
   2026-10-07:* an earlier version of this point said the
   covariance-moderation signal was weighted toward the
   early-trajectory observations; paper 07 shows it is not.)

5. **Non-orthogonal `disp` and `rate`.** The displacement and
   rate parameters are correlated in the Gompertz likelihood
   surface [7]. For interaction-test
   simulation, this matters when sensitivity sweeps vary one of
   these parameters. A sweep that walks `rate` while holding
   `disp` fixed implicitly co-varies the inflection time and
   the early slope; the resulting effect on the interaction
   power is a mixture of two geometric properties of the curve,
   not a clean single-axis sensitivity. The carryover-sensitivity
   and nof1-design-sensitivity papers do not currently
   decompose this. A cleaner reparameterization (Tjørve's
   alternative form, or fixing the inflection time as the
   primary parameter and deriving `disp` as a function of it)
   would produce sensitivity sweeps with clearer interpretation.

**Low impact on the interaction test.**

6. **Vertical-offset rescaling instability.** A numerical
   concern at extreme parameter values. The default
   $(\text{disp}=5, \text{rate}=0.42)$ regime is far from the
   instability region; the interaction test is not affected at
   the published parameter settings. A guard in the code is
   prudent but the issue does not undermine current results.
   (Corrected 2026-10-08: the instability is cancellation in the
   numerator at very small `disp`, not the rescaling factor; an
   `expm1` form removes it. See Section 3.)

7. **Naming overlap with Zwietering's modified Gompertz.** A
   documentation concern, not a statistical concern. The
   interaction test is unaffected.

8. **Departure from Emax/Hill pharmacodynamic tradition.** A
   reviewer-facing concern about defensibility of the choice. The
   interaction test is unaffected at fixed parameters; it would
   change under an Emax DGP only through the same mechanism as
   point (3) above.

**Summary.** The Gompertz-related issues that most warrant
attention before the next interaction-test manuscript revision
are the fixed plateau (point 1), the monotonicity assumption
(point 2), and the absence of a sensitivity sweep over sigmoid
family (point 3). Points 4 and 5 are second-order. Points 6
through 8 are housekeeping. A direct sensitivity sweep at fixed
$(\text{maxr}, \text{inflection time})$ across Gompertz, logistic
4PL, and Hill sigmoids would address point 3 and is the highest-
value next experiment for the founding interaction-test question.
(*Note, 2026-10-07:* paper 07 ran a version of it with the
logistic, hyperbolic-tangent and piecewise-linear families; the
Hill family was not included.)

## 7. Alternative parametric families to consider

| Family | Strength | Weakness | Reference |
|---|---|---|---|
| Logistic (4PL) | Symmetric inflection at 50% asymptote | Same monotonicity problem as Gompertz | [11] |
| Hill / Emax | Closed-form pharmacological interpretation; standard in PK-PD | Unbounded slope at $t = 0$ when $h < 1$; less natural for time courses than for dose-response | [8] |
| Weibull cumulative | Tunable shape parameter, recovers the monoexponential approach at $k=1$ | Three-parameter and identifiability is similar to Gompertz | [12] |
| Bi-exponential | Two-phase response (rapid plus slow) | Non-monotone rate of approach; estimation can be unstable | [10] |
| Spline / GAM mean structure | Non-parametric, captures relapse and treatment escape | Loses the parsimony argument; harder to calibrate to a single target trajectory | [13] |
| Latent change score / GMM | Allows participant-level trajectory heterogeneity in shape | Heavy machinery; identifiability requires more timepoints than a typical N-of-1 trial provides | [14] |
| Stochastic differential equation | Endogenous treatment-escape and relapse | Loses closed-form simulation; substantially heavier compute | [22] |

Table: Alternative parametric response families with strengths, weaknesses, and references

## 8. Recommendations for pmsimstats-ng

1. **Document the Gompertz choice as a modeling assumption,
   not a mechanism.** The pmsimstats-ng manuscripts (notably the
   companion to this document at
   `analysis/report/01-dgp-mean-moderation-vs-mvn/` and the
   sensitivity papers at `analysis/report/02-carryover-sensitivity/`
   and `analysis/report/05-nof1-design-sensitivity/`) should state
   in the simulation-design section that the Gompertz is chosen
   for parsimony and tractability, that its sigmoidal shape is
   defensible but not the only defensible choice, and that the
   power estimates are conditioned on this functional form.

2. **Add a Gompertz-versus-alternative sensitivity sweep.** A
   direct test of robustness would replace the Gompertz BR with a
   logistic 4PL and a Hill-form sigmoid at matched
   asymptote-and-inflection coordinates and report the resulting
   power changes. The infrastructure for such a sweep already
   exists in the
   `analysis/scripts/{carryover,nof1-design}-sensitivity/`
   pipelines and would require an additional axis in the cell
   tibble plus a small change to `simulation-core.R` to dispatch
   on response family. (*Note, 2026-10-07:* family dispatch now
   exists in the package, `trajectoryShape()` with the `br_family`
   argument of `generateData()`, and paper 07's driver
   `02-faithful-trajectory-sweep.R` (in
   `analysis/scripts/gompertz-evaluation/`) uses it; the Hill form is not among the implemented families.)

3. **Rename `modgompertz` or document the deviation from the
   Zwietering modified form.** Either the function name should be
   changed to reflect the actual modification (vertical-offset and
   rescale to origin) or the in-code roxygen documentation should
   make the deviation from Zwietering explicit. At present, a
   reader who recognizes 'modified Gompertz' from the
   bacterial-growth literature will encounter a different
   parameterization under the same name.

4. **Guard the rescale factor.** Add an explicit check on
   `(maxr - vert_offset)` for near-zero values and either issue a
   warning or fall back to the unmodified Gompertz at low `disp`.

5. **Consider participant-level variation in `maxr`.** A
   random-effects extension in which each simulated participant
   draws an asymptote $a_i \sim N(a_0, \tau_a^2)$ (code `maxr`) would
   accommodate plateau heterogeneity without abandoning the
   Gompertz family. This would also bring the framework into
   better alignment with the latent-class and finite-mixture
   formulations developed in the companion paper at
   `analysis/report/03-latent-class-mixture-application/`.

6. **Acknowledge response-then-relapse as a scope exclusion.**
   Until the framework supports non-monotone response trajectories,
   manuscripts should state explicitly that the simulation does
   not capture mid-treatment loss of benefit, and that the
   reported power estimates may be optimistic for clinical
   populations in which such loss is documented.

## 9. Summary

The Gompertz model is a defensible parametric template for the
simulated symptom trajectories in pmsimstats-ng. It is
parsimonious, numerically tractable, and has a long pedigree in
adjacent biological domains. It is also a strong assumption: it
imposes monotone non-decreasing trajectories with a fixed
population asymptote and a specific asymmetric inflection that may
not match the empirical PTSD or chronic-pain symptom-trajectory
literature. When this note was written (2026-05-07) the
pmsimstats-ng manuscripts used the Gompertz without a sensitivity
comparison against the established Emax / Hill alternative or
against logistic and Weibull-cumulative families. (*Note,
2026-10-07:* paper 07 now supplies a comparison against the
logistic, hyperbolic-tangent and piecewise-linear families at one
design and $t_{1/2} = 0$; Emax / Hill and Weibull remain untested.)
We recommend that future revisions
either defend the Gompertz choice empirically against a candidate
alternative or report a sensitivity sweep that demonstrates the
power conclusions are robust to the choice of sigmoid family. The
Gompertz alone should not bear the weight of a methodological
recommendation about trial design.

## References

Metadata were checked on 2026-10-07 against PubMed (entries 2 to 10,
15 to 21, 23, 24) and Crossref (entries 1, 11 to 13); entries 14 and
22 are books checked against the publisher record.

1. Gompertz B. On the nature of the function expressive of the law
   of human mortality, and on a new mode of determining the value of
   life contingencies. *Philosophical Transactions of the Royal
   Society of London*. 1825;115:513-583. doi:10.1098/rstl.1825.0026.
2. Winsor CP. The Gompertz curve as a growth curve. *Proceedings of
   the National Academy of Sciences*. 1932;18(1):1-8.
3. Laird AK. Dynamics of tumor growth. *British Journal of Cancer*.
   1964;13(3):490-502. doi:10.1038/bjc.1964.55.
4. Norton L, Simon R. Tumor size, sensitivity to therapy, and design
   of treatment schedules. *Cancer Treatment Reports*.
   1977;61(7):1307-1317.
5. Norton L. A Gompertzian model of human breast cancer growth.
   *Cancer Research*. 1988;48(24 Pt 1):7067-7071.
6. Zwietering MH, Jongenburger I, Rombouts FM, van 't Riet K.
   Modeling of the bacterial growth curve. *Applied and
   Environmental Microbiology*. 1990;56(6):1875-1881.
7. Tjørve KMC, Tjørve E. The use of Gompertz models in growth
   analyses, and new Gompertz-model approach: an addition to the
   Unified-Richards family. *PLOS ONE*. 2017;12(6):e0178691.
8. Holford NHG, Sheiner LB. Understanding the dose-effect
   relationship: clinical application of
   pharmacokinetic-pharmacodynamic models. *Clinical
   Pharmacokinetics*. 1981;6(6):429-453.
9. Sheiner LB. Clinical pharmacology and the choice between theory
   and empiricism. *Clinical Pharmacology and Therapeutics*.
   1989;46(6):605-615. doi:10.1038/clpt.1989.195.
10. Holford N. Pharmacodynamic principles and target concentration
    intervention. *Translational and Clinical Pharmacology*.
    2018;26(4):150-154. doi:10.12793/tcp.2018.26.4.150.
11. Richards FJ. A flexible growth function for empirical use.
    *Journal of Experimental Botany*. 1959;10(2):290-301.
    doi:10.1093/jxb/10.2.290.
12. Weibull W. A statistical distribution function of wide
    applicability. *Journal of Applied Mechanics*.
    1951;18(3):293-297. doi:10.1115/1.4010337.
13. Brumback BA, Rice JA. Smoothing spline models for the analysis
    of nested and crossed samples of curves. *Journal of the
    American Statistical Association*. 1998;93(443):961-976.
    doi:10.1080/01621459.1998.10473755.
14. Muthén B. Latent variable mixture modeling. In: Marcoulides GA,
    Schumacker RE, eds. *New Developments and Techniques in
    Structural Equation Modeling*. Lawrence Erlbaum; 2001:1-33.
15. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
    aggregated N-of-1 trial designs for predictive biomarker
    validation: statistical methods and theoretical findings.
    *Frontiers in Digital Health*. 2020;2:13.
    doi:10.3389/fdgth.2020.00013.
16. Raskind MA, Peskind ER, Kanter ED, et al. Reduction of
    nightmares and other PTSD symptoms in combat veterans by
    prazosin: a placebo-controlled study. *American Journal of
    Psychiatry*. 2003;160(2):371-373.
17. Raskind MA, Peskind ER, Hoff DJ, et al. A parallel group
    placebo controlled study of prazosin for trauma nightmares and
    sleep disturbance in combat veterans with post-traumatic stress
    disorder. *Biological Psychiatry*. 2007;61(8):928-934.
18. Raskind MA, Peterson K, Williams T, et al. A trial of prazosin
    for combat trauma PTSD with nightmares in active-duty soldiers
    returned from Iraq and Afghanistan. *American Journal of
    Psychiatry*. 2013;170(9):1003-1010.
19. Raskind MA, Peskind ER, Chow B, et al. Trial of prazosin for
    post-traumatic stress disorder in military veterans. *New
    England Journal of Medicine*. 2018;378(6):507-517.
20. Steckler T, Risbrough V. Pharmacological treatment of PTSD:
    established and new approaches. *Neuropharmacology*.
    2012;62(2):617-627.
21. Hoskins M, Pearce J, Bethell A, et al. Pharmacotherapy for
    post-traumatic stress disorder: systematic review and
    meta-analysis. *British Journal of Psychiatry*.
    2015;206(2):93-100.
22. Øksendal B. *Stochastic Differential Equations: An Introduction
    with Applications*. 6th ed. Springer; 2003.
23. Kent DM, Paulus JK, van Klaveren D, et al. The Predictive
    Approaches to Treatment effect Heterogeneity (PATH) Statement.
    *Annals of Internal Medicine*. 2020;172(1):35-45.
    doi:10.7326/M18-3667.
24. Rothwell PM. Treating individuals 2. Subgroup analysis in
    randomised controlled trials: importance, indications, and
    interpretation. *Lancet*. 2005;365(9454):176-186.
    doi:10.1016/S0140-6736(05)17709-5.
