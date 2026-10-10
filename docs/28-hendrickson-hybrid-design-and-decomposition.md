# The Hendrickson approach, the hybrid design, and when component decomposition pays off {.unlisted .unnumbered}

*2026-06-17 10:40 PDT*

> **SUPERSEDED (2026-09-09). Do not cite or build on this
> document.** Its content is merged into
> `docs/31-component-decomposition-and-reductions.md`, which also
> adds Part IV on the two-component reductions (papers 13 and 14).
> This copy is retained for history only and received a light
> review on 2026-10-07 (notation, obvious errors, US English); its
> numbers were not re-verified. Consult docs/31 instead.

**Author.** pmsimstats team

**Purpose.** Historical record of the June 2026 synthesis of the
original Hendrickson approach, the identifiability the hybrid
design supplies, and the paper-06 conclusion on when component
decomposition repays its cost. Superseded by docs/31.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. Code identifiers
appear in backticks.

| Symbol | Meaning | Status |
|---|---|---|
| $BR$, $PB$, $TV$ | biological, placebo-belief and time-variant response components | canonical |
| $Y$ | observed symptom score (code `Sx`) | canonical |
| $\mathrm{BL}_i$ | participant baseline level | canonical |
| $\eta_t$ | design expectancy at visit $t$: 1 open label, 0.5 blinded, 0 covert (code `e`, `expectancies`) | canonical |
| $D_{bc}$ | exposure-decayed drug indicator (code `Dbc`) | canonical |
| $\lambda$, $t_{1/2}$, $t_{sd}$ | carryover decay rate $\ln 2 / t_{1/2}$; half-life; time since discontinuation | canonical |
| $\beta_{bm}^{\text{lumped}}$, $\beta_{bm}^{BR}$, $\beta_{bm}^{PB}$, $\beta_{bm}^{TV}$ | one-component and component-specific biomarker slopes | canonical |
| $w_{PB}$, $w_{TV}$ | design-dependent weights in the omitted-variable-bias identity | local |
| $m_{BR}$, $m_{PB}$, $m_{TV}$ | Gompertz maxima of the three components (code `max`) | local |
| $r_{BR}$, $d_{BR}$ | Gompertz rate and displacement of $BR$ (code `rate`, `disp`) | local |
| $G_C(t)$ | modified-Gompertz curve of component $C$ divided by its maximum $m_C$; $d$, $r$ without subscript are displacement and rate | local |
| $I_\ell$ | indicator of phase level $\ell$ (open label OL is the reference) | local |
| $c_{bm,PB}$ | biomarker-$PB$ coupling correlation in the contamination studies | local |
| $N$ | total sample size | canonical |
| $\alpha$ | significance level | canonical |

Table: Notation index and glossary of symbols used in this document

- **One-component (lumped) analysis.** The Hendrickson model with a
  single drug indicator, `Sx ~ bm + t + Dbc + bm:Dbc`.
- **Attribution.** Apportioning response among pharmacology
  ($BR$), belief ($PB$) and natural history ($TV$).

## 1. Purpose and scope

This white paper consolidates, in one place, three things that are
otherwise distributed across the package source and the paper-06
manuscript (`analysis/report/06-component-decomposition/`):

1. what the original Hendrickson et al. (2020) approach actually
   does, on both the data-generating and the analysis side;
2. how the hybrid open-label-plus-blinded-discontinuation design
   supplies identifiability, and what it cannot supply;
3. the conclusion the paper-06 analysis draws about when component
   decomposition repays its cost, with explicit, practical advice on
   how and when to field the hybrid design.

The intended reader is an investigator deciding whether to run an
aggregated N-of-1 trial under the hybrid design and whether to layer
a component decomposition on top of the standard lumped analysis. The
treatment is academic and self-contained; the simulation numbers are
quoted from the realized Study A, the Study B pilot, and Study C
reported in the paper-06 manuscript.

## 2. The original Hendrickson approach

The founding pmsimstats question is the detection of a
biomarker-by-treatment interaction in an aggregated N-of-1 trial. The
Hendrickson reference implementation
(`implementations/original/`) is the historical anchor for the whole
program. Its defining feature, and the starting point for paper 06,
is an asymmetry between how it generates data and how it analyzes it.

### 2.1 The data-generating process is fully decomposed

The DGP in `implementations/original/R/generateData.R` builds the
observed symptom score as the sum of three causally distinct
components, each a modified Gompertz trajectory with
participant-specific random effects:

- $BR$ (biological response, the pharmacological component), driven
  by cumulative time on drug, decaying off drug with a carryover
  half-life;
- $PB$ (placebo-belief response), driven by cumulative time on
  belief and scaled by the phase expectancy $\eta_t$;
- $TV$ (time-variant natural-history component), driven by time
  since trial entry.

The components are summed and subtracted from baseline (lines near
99 and 104 of `generateData.R`):

```{.r}
D_t  = rowSums(t.br + t.pb + t.tv)   # total response at timepoint t
Sx_t = BL - D_t                       # symptom = baseline minus total
```

The raw simulated data retain the individual component values
(`t1.br`, `t1.pb`, `t1.tv`, and so on for each timepoint), so the
ground truth is three-component by construction. The prazosin-PTSD
calibration used throughout the program draws its Gompertz
parameters from the Hendrickson et al. (2020) dataset
($m_{BR} \approx 11$, $r_{BR} \approx 0.42$ per week,
$d_{BR} \approx 5$).

### 2.2 The analysis model is a single lumped drug indicator

The analysis in `implementations/original/R/lme_analysis.R` does not
mirror the DGP. It fits one continuous drug indicator, not three
components:

```{.r}
Sx ~ bm + t + Dbc + bm:Dbc
random      = ~ 1 | ptID
correlation = corCAR1(form = ~ t | ptID)
```

where `Dbc` is the continuous drug state $D_{bc}$ (equal to 1 on
drug and to $e^{-\lambda t_{sd}}$ off drug, with
$\lambda = \ln 2 / t_{1/2}$), and
`bm:Dbc` is the biomarker-by-drug interaction. Variant forms in the
same file add an expectancy indicator (`De`) or a carryover term
(`tsd`), but they remain single-drug-indicator models. The model
estimates one lumped drug main effect and one lumped biomarker
slope; it does not recover $BR$, $PB$, and $TV$ separately.

### 2.3 The asymmetry that paper 06 studies

The original therefore applies decomposition only as a
data-generating device. It knows the response is $BR + PB + TV$ and
simulates accordingly, but its inferential model collapses the three
into a single treatment coefficient. This lumped estimator is exactly
the baseline that paper 06 formalizes, and the central question of
that paper is the consequence of the asymmetry: when does the lumped
analysis of three-component data give a biased or uninformative
answer, and when does adding analysis-side decomposition fix it?

## 3. The hybrid design

Component decomposition is meaningless unless the data contain enough
information to pin down each component separately. That information is
supplied by trial design, not by the analysis model. The hybrid
open-label-plus-blinded-discontinuation design is the worked example
throughout paper 06 and the primary identifiability vehicle.

### 3.1 Phase structure

The hybrid design strings together three phases, each contributing a
different mixture of the three components:

| Phase | Drug | Belief | Components present |
|---|---|---|---|
| Open-label on drug | Yes | Full ($\eta_t = 1$) | $BR + PB(\eta_t = 1) + TV$ |
| Blinded discontinuation | Hidden | Hedged ($\eta_t = 0.5$) | $BR$ (on drug) $+ PB(\eta_t = 0.5) + TV$ |
| Blinded crossover | Either | Hedged ($\eta_t = 0.5$) | $BR$ (on drug) $+ PB(\eta_t = 0.5) + TV$ |

Table: Phases of the Hybrid design, with drug state, expectancy and the response components present in each

The 'on drug' qualifier on $BR$ is shorthand for the continuous-decay
drug state $D_{bc}$, not a binary switch. (*Correction, 2026-10-07:*
the third row previously read 'Open-label crossover, knows,
$\eta = 1$'. In the Hybrid design as coded, the expectancy vector is
`c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)`, so the crossover is blinded
with $\eta_t = 0.5$; inspected in
`analysis/scripts/carryover-sensitivity/simulation-core.R` and the
`figure4` drivers. The same row in docs/31 needs the same fix.)

### 3.2 The three identifying contrasts

Three phase-by-allocation contrasts identify the three components:

- **Open-label versus blinded (within drug)** yields
  $PB(\eta_t = 1) - PB(\eta_t = 0.5)$, which is approximately
  $0.5\,PB$, and so
  isolates the belief component, conditional on a known blinded-phase
  expectancy. (Corrected 2026-10-08: the approximation holds only to
  the extent that the component curves change little between week 8
  and week 9. The within-drug change is $m_{BR}[G_{BR}(9) -
  G_{BR}(8)] + m_{PB}[0.5\,G_{PB}(9) - G_{PB}(8)] + m_{TV}[G_{TV}(9) -
  G_{TV}(8)]$; under the published $PB$ parameters ($m_{PB} = 6.51$,
  $d = 5$, $r = 0.35$) the $PB$ term alone is 2.17 points, against
  $0.5\,E[PB]$ of 2.40 at week 8 and 2.62 at week 9, and the $BR$ and
  $TV$ terms do not vanish (derived; recomputed 2026-10-08). See
  `analysis/report/compendium-mathematics.tex`, Section "Time
  arguments of the components, the shared clock, and the phase
  contrasts" (`sec:invb-clocks`).)
- **Blinded on-drug versus blinded placebo** yields $BR$ once drug
  carryover has elapsed, because $PB$ and $TV$ match across allocation
  arms within the blinded phase.
- **Off-drug phases** yield $TV$ plus any residual off-drug $PB$, with
  no active $BR$ contribution beyond carryover.

The lumped LME combines these into one model, but the inference is
driven by the design. The one-component analysis already exploits the
drug-on/off contrast through its `Dbc` term; that is why, as Section 5
shows, it is not as impoverished as its single coefficient suggests.

### 3.3 The expectancy ($\eta$) identification limit

The hybrid design has one principled limitation that deserves
emphasis. With only an open-label and a blinded phase, the data
identify the product $\eta\,m_{PB}$, not the two factors separately. The
working value $\eta$ approximately 0.5 in the blinded phase is a
modeling assumption, not an estimate, and absent a belief measurement
it is not falsifiable from outcome data alone. The identified set for
$m_{PB}$ is the interval obtained by dividing the estimated product by
$\eta$ over its plausible range, collapsing to a point only when $\eta$
is pinned down. Two design features close the gap:

- an open-label placebo phase, which supplies a second equation that
  separates $\eta$ from $m_{PB}$;
- direct measurement of treatment expectancy, as in the open-hidden
  and balanced-placebo designs that manipulate or record belief.

Any trial relying on the decomposition for a PB-versus-BR attribution
should include such a phase or collect an expectancy measure, and
should report sensitivity of the attribution to $\eta$ across its
plausible range.

### 3.4 Identifiability is necessary but not sufficient

Design-contrast identifiability does not guarantee stable estimation
at finite sample sizes. The three Gompertz trajectories share an onset
shape and overlap in calendar time, so their parameters can be
near-collinear even when the identifying contrasts are present. Joint
recovery additionally requires a phase in which each component varies
while the others are held fixed, and an off-drug window long enough,
relative to the carryover half-life, for $TV$ to separate from a still
decaying $BR$. Where these conditions are weak the components are
identified only up to a near-collinear direction; the symptom is the
rank deficiency the Study B pilot encountered when the full
phase-augmented formula was fitted at $N = 35$. (Corrected
2026-10-08: that rank deficiency is not a near-collinearity of the
component trajectories. It is an exact aliasing of the
phase-augmented formula, present at every $N$: because every
open-label visit is on drug, $D_{bc} = 1$ in the reference phase and
$D_{bc} = 1 - \sum_{\ell \ne \mathrm{OL}} I_\ell +
\sum_{\ell \ne \mathrm{OL}} D_{bc} I_\ell$ exactly, so the `Dbc`
column is a linear combination of the intercept, the `phase` columns
and the `Dbc:phase` columns. `nlme::lme` stops with "Singularity in
backsolve" on the full formula (verified), and the paper 06 driver
(`phase_augmented_fit()` in
`analysis/scripts/component-decomposition/00-common.R`) refits
without `Dbc:phase`, and omits the full formula without attempting it
for $N < 70$ (inspected). Verified by the rank of the model matrix,
`analysis/scripts/compendium-mathematics/check-gompertz-rmanova-cycles.R`;
see `analysis/report/compendium-mathematics.tex`, Section "Rank of
the phase-augmented model, and the drug-expectancy correlation by
design" (`sec:invb-collinearity`). The general point of this
subsection about near-collinear Gompertz trajectories stands as an
argument, but this symptom is not evidence for it.) Recoverability is to be
mapped by identifiability diagnostics (convergence rates,
variance-inflation factors on the component-specific parameters, and
sensitivity of the $BR$ estimate to the assumed $\eta$), not asserted.

## 4. The paper-06 analysis

### 4.1 The omitted-variable-bias identity

The analytic core of paper 06 is an exact identity for what the
lumped estimator estimates. By linearity of covariance, the lumped
biomarker-by-treatment slope decomposes as

$$
\beta_{bm}^{\text{lumped}} = \beta_{bm}^{BR} + w_{PB}\,\beta_{bm}^{PB}
+ w_{TV}\,\beta_{bm}^{TV}
$$

so the lumped slope is displaced from the pharmacological slope
$\beta_{bm}^{BR}$ by a term proportional to the biomarker's covariance with
$PB$ and $TV$, weighted by design-dependent weights $w_{PB}$ and $w_{TV}$.
The displacement is in the probability limit, not the variance: a
larger one-component trial converges on the same biased value. The
bias is zero if and only if the biomarker is orthogonal to $PB$ and
$TV$, regardless of the magnitude of those components. (Corrected
2026-10-08: "if and only if" overstates the condition. The bias is
zero exactly when $w_{PB}\beta_{bm}^{PB} + w_{TV}\beta_{bm}^{TV} = 0$;
orthogonality implies this, but a cancellation between the two terms
could also satisfy it, so orthogonality is sufficient, not necessary.
See `analysis/report/compendium-mathematics.tex`, Section
"Partitioning the treatment effect" (`sec:partition`).) The companion
identity for the lumped treatment effect is displaced by the $PB$ and
$TV$ contributions themselves, so a lumped pre-post treatment effect
is biased whenever $PB$ or $TV$ is non-zero.

### 4.2 Three analysis strategies compared

Paper 06 compares three analysis strategies against three-component
data:

1. **One-component (lumped).** The Hendrickson model,
   `Sx ~ bm + t + Dbc + bm:Dbc`.
2. **Phase-augmented LME.** A modest extension,
   `Sx ~ bm + t + Dbc + phase + Dbc:phase + bm:Dbc + bm:Dbc:phase`,
   which tests whether the drug and biomarker-by-drug effects attenuate
   under blinding, without the full Gompertz parametric form.
   (Corrected 2026-10-08: on the Hybrid design this formula is
   rank-deficient at every $N$ (Section 3.4), and every
   phase-augmented result that paper 06 reports, including those in
   Section 4.3 below, is from the fallback formula without
   `Dbc:phase`, `Sx ~ bm + t + Dbc + phase + bm:Dbc + bm:Dbc:phase`.
   In that formula the reported `bm:Dbc` coefficient is an
   open-label versus blinded contrast, in which drug state and
   expectancy change together; see
   `analysis/report/compendium-mathematics.tex`,
   `sec:invb-collinearity`.)
3. **Full component-matched decomposition.** A parametric model that
   mirrors the DGP and recovers component-specific slopes, feasible only
   inside the recoverable region of Section 3.4.

### 4.3 What the simulations show

**Study A (orthogonal-biomarker arm, hybrid design; 1,000 replicates
per alternative cell, 5,000 per null).** With the biomarker coupled to
$BR$ only ($c_{bm,PB} = 0$), the one-component estimate of the
biomarker-by-treatment slope is essentially unbiased across the entire
$(m_{PB}, m_{TV})$ grid: mean absolute bias 0.016 at $N = 35$ and 0.006 at
$N = 150$ against a true value of -2.25, both within one Monte Carlo
standard error of zero, with near-nominal type I error (0.031 to 0.056
for the one-component analysis). The phase-augmented analysis reduces
no bias (0.028 versus 0.016 at $N = 35$) and is markedly less powerful
(mean power 0.35 for the phase-augmented analysis versus 0.59 for
the one-component analysis at $N = 35$; 0.62 versus 0.90 at
$N = 70$; both reach approximately 1.0 by $N = 100$; order of the
pairs corrected 2026-10-07 against paper 06, `report-slim.Rmd`).
Under an
uncontaminated biomarker, decomposition is therefore strictly
dominated.

**Study B (coupling-arm pilot; 100 replicates per cell).** With the
biomarker coupled to $PB$, the one-component bias scales with the
coupling strength $c_{bm,PB}$, not with the magnitude of $PB$, exactly as
the identity predicts (the one-component bias is approximately 0.51 at
the contamination level studied). The phase-augmented analysis offers
no protection: its bias is statistically indistinguishable from the
one-component bias at every contamination level. (Note 2026-10-08:
the phase-augmented results of Studies A and B are from the fallback
formula without `Dbc:phase`; see Section 4.2.)

**Study C (recovery under a belief-decoupling design; 1,000 replicates
per cell, $N$ in {70, 150}, $c_{bm,PB}$ in {0, 0.1, 0.2, 0.3}).** This is
the decisive design comparison. On the standard hybrid design the
drug-state and belief contrasts are collinear (analysis-stage
correlation approximately +0.6; note 2026-10-08: this figure, taken
from paper 06, could not be reproduced. For the 16-visit two-path
design of the recovery study, pooled over visits with equal path
weights, $\mathrm{Cor}(D, \eta) = 0.775$ with binary $D$ and 0.753
with $D_{bc}$ at $t_{1/2} = 1$ (0.785 with a baseline row
$D = \eta = 0$); for the 8-visit four-path Hybrid it is 0.509 and
0.485 (0.645 with the baseline row). The balanced-placebo value is
reproduced: exactly $-2/3$ with binary $D$ and $-0.656$ with
$D_{bc}$. Verified,
`analysis/scripts/compendium-mathematics/check-gompertz-rmanova-cycles.R`;
see `analysis/report/compendium-mathematics.tex`,
`sec:invb-collinearity`), and the belief-covariate,
blinded-stratum, and Gompertz-basis decompositions all fail to recover
$\beta_{bm}^{BR}$, no better than the one-component analysis. A
balanced-placebo (open-hidden) design adds a covert on-drug phase
($\eta_t = 0$) and an open-placebo phase ($\eta_t = 1$), reducing the
drug-state-by-belief correlation to approximately -0.66 and breaking
the collinearity. On that design the belief-covariate decomposition
returns the unbiased pharmacological slope across the feasible
contamination range (bias -0.022, +0.008, +0.021, -0.026 at $N = 150$,
MCSE approximately 0.02; coverage 0.94 to 0.96), whereas the
one-component estimator degrades monotonically (bias +0.16, +0.34,
+0.48; coverage down to 0.86). The grid is capped at $c_{bm,PB} = 0.30$
because beyond approximately 0.45 the implied covariance is
non-positive-definite.

The power comparison between the two estimators on the balanced-placebo
design is more nuanced than a uniform decomposition penalty (rejection
rate for the interaction at $\alpha = 0.05$, computed from the
1,000-replicate per-cell p-values; paired standard error of the gap at
most 0.014):

| $c_{bm,PB}$ | $N$ | One-component | Decomposition | Gap |
|---|---|---|---|---|
| 0.0 | 70 | 0.85 | 0.66 | +0.19 |
| 0.1 | 70 | 0.78 | 0.67 | +0.11 |
| 0.2 | 70 | 0.74 | 0.69 | +0.05 |
| 0.3 | 70 | 0.64 | 0.67 | -0.03 |
| 0.0 | 150 | 0.99 | 0.93 | +0.07 |
| 0.1 | 150 | 0.98 | 0.93 | +0.05 |
| 0.2 | 150 | 0.97 | 0.93 | +0.04 |
| 0.3 | 150 | 0.93 | 0.93 | -0.00 |

Table: Power for the interaction on the balanced-placebo design, one-component versus belief-covariate decomposition (Study C, 1,000 replicates per cell)

At zero contamination the one-component analysis is the more powerful of
the two (by 0.19 at $N = 70$, 0.07 at $N = 150$), reflecting the degrees of
freedom the decomposition spends. That advantage erodes monotonically as
contamination grows and reverses by $c_{bm,PB} = 0.3$ (0.64 versus 0.67 at
$N = 70$), because contamination in this direction drags the biased
one-component estimate toward the null and so costs it power, while the
decomposition holds the estimate at the pharmacological target. The
one-component power advantage is therefore confined to the
low-contamination regime in which decomposition is unnecessary. (The sign
of the gap depends on the direction of the coupling: a coupling that
reinforced the BR sign would inflate the one-component estimate and give
it spuriously higher power, which is power to detect a contaminated
quantity rather than a genuine advantage.)

On the hybrid design the same comparison is uninformative because power
saturates near 1.0 for both estimators (the BR effect of -2.25 is large
relative to its standard error at $N = 70$ to 150), and both are biased
under contamination in any case; there the operative axis is bias and
coverage, not power.

## 5. The conclusion drawn from the 06 analysis

The inferential value of decomposition is conditional, not universal,
and the condition is a joint function of the biomarker's correlation
structure and the trial design, not of either alone. Three points
summarize the result.

1. **The design contrast matters, not the parametric model.** The
   trial must supply the drug-on/off and belief contrasts that identify
   the components. The lumped one-component analysis already exploits
   the drug-on/off contrast through `Dbc`. Adding parametric component
   structure on top is not free: it spends degrees of freedom and, under
   an uncontaminated biomarker, reduces no bias while losing substantial
   power.

2. **For the biomarker-by-treatment slope, the bias is governed by
   covariance, not magnitude.** A biomarker orthogonal to $PB$ and $TV$
   gives an unbiased lumped slope no matter how large $PB$ and $TV$ are.
   A biomarker coupled to $PB$ or $TV$ gives a contaminated slope that no
   amount of data fixes, because the failure is identification, not
   precision.

3. **Recovery under contamination is a design property.** When the
   biomarker is coupled to $PB$, the bias is recoverable only under a
   design that varies belief independently of drug exposure. On the
   standard hybrid design no analysis recovers it; on a balanced-placebo
   design the belief-covariate decomposition does.

### The bottom line specific to the hybrid design

For the biomarker-by-treatment interaction estimand on the hybrid
design, component decomposition never wins:

- if the biomarker is orthogonal to $PB$ and $TV$, the lumped analysis
  is already unbiased and is the more powerful of the two, so
  decomposition only costs power;
- if the biomarker is coupled to $PB$, the lumped analysis is biased
  but decomposition cannot rescue it on the hybrid design, because the
  identifying contrasts are collinear.

Decomposition earns its place on the hybrid design only for a
different estimand: the **attribution** of response to pharmacology
versus belief versus natural history. For attribution the hybrid
design is the primary identification vehicle, subject to the $\eta$
limit of Section 3.3 and the collinearity caveat of Section 3.4.

### Was the decomposition in Hendrickson's DGP unnecessary?

A natural reading of the paper-06 result is that, because the lumped
analysis suffices for the biomarker-by-treatment estimand, the
three-component structure Hendrickson built into the simulator was
itself superfluous. That reading conflates two distinct choices. The
decomposition was unnecessary in the *analysis* model; it was
load-bearing in the *data-generating process*. The two claims are not
interchangeable.

The DGP decomposition was necessary, or at least not removable, for two
reasons.

- It defines the ground truth the whole program tests against. The
  omitted-variable-bias identity of Section 4.1, and every simulation
  that confirms it, can be exhibited only if the simulator generates
  $BR$, $PB$, and $TV$ as separately controllable channels. A DGP that
  produced a single lumped treatment effect would offer no $c_{bm,PB}$
  knob to turn, no way to construct the contamination, and therefore no
  way to show that the lumped analysis is unbiased under orthogonality
  but biased under coupling. The DGP decomposition is what makes the
  bias question askable.
- It encodes genuinely distinct mechanisms with different time-drivers:
  $BR$ by cumulative time on drug with off-drug carryover decay, $PB$ by
  time on belief scaled by the phase expectancy $\eta$, and $TV$ by time
  since entry. These produce different trajectories across trial phases,
  and the design-contrast structure of Section 3.2 (open-label versus
  blinded yielding approximately $0.5\,PB$, blinded on-drug versus
  blinded placebo yielding $BR$) cannot be reproduced from a single
  lumped curve. The components are not cosmetic; they drive the phase
  behavior the whole identification argument relies on.

What paper 06 shows to be unnecessary is only the symmetric move of
mirroring that structure in the analysis. Fitting the components back
out is unnecessary, and under an orthogonal biomarker counterproductive,
for the biomarker-by-treatment estimand, and unrecoverable on the hybrid
design when the biomarker is coupled. Hendrickson's choice to analyze
with a single lumped `Dbc` was, in hindsight, the appropriate one for
that estimand.

The qualifier is the estimand. The analysis-side lumping is validated as
sufficient only for the interaction-slope target under an orthogonal
biomarker. Had the scientific target been attribution, or a $PB$-coupled
biomarker, the lumped analysis would be insufficient, and the DGP's
decomposition would then require both a matching analysis and a
belief-decoupling design. The accurate summary is therefore that the DGP
decomposition was necessary to pose and answer the question, whereas the
analysis-side decomposition was unnecessary for the specific estimand
Hendrickson actually pursued.

## 6. Advice: how and when to use the hybrid design

The hybrid design is the right default for aggregated N-of-1 trials
whose primary estimand is the pharmacological treatment effect or a
biomarker-by-treatment interaction with a mechanistically
pharmacological biomarker. It is not the right design when the
biomarker may couple to placebo responsiveness and the pharmacological
slope must be recovered. The following guidance operationalizes that.

**Use the hybrid design, with the lumped one-component analysis, when:**

- the estimand is the lumped treatment effect or a biomarker-by-drug
  interaction, and
- the candidate biomarker is mechanistically orthogonal to belief and
  to natural history. Pharmacokinetic or receptor-binding biomarkers
  (for example a CYP metabolizer phenotype, or an adrenergic-receptor
  variant affecting prazosin binding) qualify: enzymatic metabolism and
  receptor affinity cannot plausibly couple to belief state or untreated
  disease trajectory. Here the lumped slope is already the target, and a
  parametric component model is counterproductive.

**Use the hybrid design, but add an open-label placebo phase or an
expectancy measure, when:**

- the estimand is attribution (how much of the response is
  pharmacological versus belief versus natural history). The base
  two-phase hybrid identifies only the product $\eta\,m_{PB}$; the
  open-label placebo phase or a recorded expectancy separates them. Report
  attribution sensitivity to $\eta$ across its plausible range, and verify
  the off-drug window is long enough, relative to the carryover
  half-life, for $TV$ to separate from a decaying $BR$.

**Do not rely on the hybrid design; field a balanced-placebo
(open-hidden) design instead, when:**

- the estimand is the pharmacological biomarker slope and the
  candidate biomarker may correlate with placebo responsiveness or with
  natural-history drift. Psychological or symptom-severity biomarkers
  fall here: baseline PCL-5 severity, for example, couples to both
  expectancy and regression-to-the-mean $TV$, and a COMT genotype is a
  documented predictor of placebo-response magnitude. On the hybrid
  design the contamination is unrecoverable; only the covert-drug and
  open-placebo phases of the balanced-placebo design decouple belief from
  drug exposure enough to recover the unbiased slope.

**Design-stage checklist before committing to the hybrid design.**

- Identify the estimand first: lumped effect, attribution, or
  pharmacological biomarker slope. The estimand, not convenience,
  selects the design.
- Classify the biomarker as orthogonal to or potentially coupled with
  $PB$ and $TV$, on mechanistic grounds, before enrollment. If coupling
  cannot be ruled out, the hybrid design is the wrong instrument.
- Verify that every component the estimand requires has at least one
  phase or contrast in which it is uniquely present. A pure open-label
  design cannot separate $BR$ from $PB$; a pure blinded design cannot
  reach full-strength $PB$.
- If attribution is in scope, plan the open-label placebo phase or
  expectancy measurement at design time, and size the off-drug window
  against the carryover half-life.
- Do not plan to buy bias protection with the analysis model. On the
  hybrid design the phase-augmented and full-decomposition analyses do
  not reduce biomarker-slope bias; design is the only lever.

## 7. Limitations and open items

The Study C recovery result is established at 1,000 replicates over
$c_{bm,PB}$ in {0, 0.1, 0.2, 0.3} and $N$ in {70, 150}; recovery beyond
the positive-definite range ($c_{bm,PB}$ above approximately 0.45),
participant-level component recovery, and replication under an
mean-moderation coupling remain open. The Study B
coupling-arm evidence is a 100-replicate pilot; its production run is
set out in the study's working protocol
(`analysis/scripts/component-decomposition/00-ademp-protocol.md`). The $\eta$ working value of 0.5 in the blinded phase
remains a modeling assumption absent a belief measurement. The
framework's participant-level claims are conditional on falling inside
the recoverable region of Section 3.4.

## References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker
   validation: statistical methods and theoretical findings.
   *Frontiers in Digital Health*. 2020;2:13.
   doi:10.3389/fdgth.2020.00013. (Title corrected 2026-10-07
   against PubMed; the prazosin-PTSD reference dataset and source
   of the Gompertz calibration.)
2. pmsimstats team. Paper 06, *Three-component decomposition of
   treatment response in aggregated N-of-1 trials*
   (`analysis/report/06-component-decomposition/report-slim.Rmd`),
   Sections 4 (identifiability), 5 (analysis-side methodology), 6
   (omitted-variable-bias identity), and 7 (Studies A, B, C).
3. pmsimstats team. Component-decomposition pedagogy
   (`docs/24-component-decomposition-pedagogy.md`).
4. pmsimstats team. Component decomposition and reductions
   (`docs/31-component-decomposition-and-reductions.md`), the
   successor to this document.
