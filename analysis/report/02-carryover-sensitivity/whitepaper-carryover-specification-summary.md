---
geometry: margin=0.85in
fontsize: 10pt
---

# White paper: choosing a carryover-mitigation analysis strategy

*2026-07-29 09:10 PDT*

**A two-page summary of 'Choosing a carryover-mitigation analysis
strategy for biomarker-treatment interactions in aggregated N-of-1
trials: a simulation-based tutorial and factorial evaluation'
(pmsimstats team, 2026-07-04).**

## Goal of the paper

In an aggregated N-of-1 trial the drug does not switch off when it is
discontinued; its effect decays over days to weeks, so nominally
off-drug occasions carry residual exposure. Treating them as cleanly
off-drug mis-specifies the design matrix and attenuates the
biomarker-treatment interaction test, which is the estimand these
trials exist to serve.

The analyst therefore has to decide how to represent residual exposure
in the linear predictor. The paper's goal is to settle that decision
empirically rather than by convention. Specifically, it aims to:

1. Compare three analyst-side representations of drug exposure under a
   fixed mixed model, and identify which is most robust to carryover.
2. Test whether the ranking is invariant to the assumed biomarker
   biology, the assumed decay shape, and the design, since a ranking
   that transfers would justify a single recommendation.
3. Quantify how much a wrong assumption actually costs, relative to
   other threats to power.

The motivation is a gap left by the companion paper, which showed that
the carryover penalty depends on how the biomarker effect is encoded
in the data-generating process (1-3% of power under Architecture A,
direct mean moderation; roughly 31% under Architecture B, MVN
differential correlation). That work fixed the decay shape at
exponential and used a single analysis specification, leaving the
joint sensitivity to mis-specification at *both* the data-generating
and the analyst layer uncharacterized.

### The three analysis specifications

All three share one model, `nlme::lme(Sx ~ bm * X + t)` with a random
intercept and `corCAR1` residual correlation, and differ only in how
the exposure regressor $X$ is built:

- **A1 (binary).** $X = D_{it}$, on-drug indicator; carryover is
  ignored. The dominated reference.
- **A2 (exposure-weighted).** $X = \mathrm{Dbc}_{it}$, a continuous
  predictor decaying off-drug under an assumed form and half-life.
  This is the current pmsimstats-ng default.
- **A3 (lagged-augmented).** $X = D_{it}$ plus a lagged just-off-drug
  nuisance term, the classical Jones-Kenward crossover device. No
  decay form is assumed.

### Simulation setup

An ADEMP-structured factorial crossing three DGP decay forms (linear,
exponential, Weibull at $k \in \{0.7, 1.0, 1.5\}$) with the three
specifications, under both DGP architectures, three designs (CO,
Hybrid, OL+BDC), two sample sizes ($N = 35, 70$) and three interaction
strengths: 540 cells at 500 replicates each, plus five single-axis
sensitivity blocks and a cluster-robust recalibration block. Common
random numbers within each cell make the A2-versus-A1 contrast a
paired comparison. Monte Carlo standard error on power is at most
0.022 and is reported with every estimate.

## Results

**A2 does not dominate.** This is the paper's central and deliberately
negative finding. The exposure-weighted specification attains the
highest power in only 68 of the 360 non-null cells (19%).

**Where it does win, it wins clearly.** Under Architecture B in the
two designs with a blinded-discontinuation phase (Hybrid, OL+BDC), A2
leads A1 and A3 by about 10 percentage points at the reference cell
(0.860 against 0.766 and 0.770 at $t_{1/2} = 1.0$), roughly five
MCSEs. A higher-precision 24-cell rerun at 600 trials confirms this
with a paired McNemar test: 13.7 points at the strongest OL+BDC cell,
$p = 6.6 \times 10^{-17}$.

**Where it loses, it loses badly.** Under the classical crossover
design A2 is markedly inferior (0.488 against 0.830 for A1). Under
Architecture A it is marginally dominated throughout (Hybrid,
$t_{1/2} = 1.0$: 0.960 against 0.972).

**Mechanism.** A2 recovers signal that carryover has moved into the
attenuated off-drug correlation. Under Architecture A the treatment
effect stays proportional at every occasion, so there is no attenuated
correlation to recover and the simpler binary indicator is at least as
good. Under the crossover design the within-subject correlation
already carries most of that signal, so A2 has nothing to gain and its
down-weighting of off-drug points only costs.

**A3 is essentially A1.** Its lagged term enters only as a main-effect
nuisance covariate, so the tested interaction coefficient is
identical to A1's; the two differ only through absorbed variance and
track closely throughout. The Jones-Kenward expectation that a
shape-free lagged term should beat a parametric predictor under
mis-specification is not supported in the range tested.

**Where A2 leads, the lead is robust.** Wrong assumed decay shape
never pushes A2 below A1 or A3 (worst case: a linear truth). A two- to
four-fold error in the assumed half-life moves A2 only from about 0.83
to 0.77. Stronger autocorrelation lifts all three specifications
equally, leaving the 5-6 point lead $\rho$-independent; the joint
$\rho \times t_{1/2}$ check finds the effects additive, not
interactive.

**Dropout dominates the specification choice.** At 30% MCAR dropout
all three methods collapse to 0.12-0.25 and become indistinguishable;
the A2-minus-A1 gap shrinks from roughly 0.10 to 0.03. Severity-related
(MAR) dropout is slightly worse than MCAR because the
highest-signal patients leave first.

**The test is mildly conservative, and fixable.** Mean null rejection
is 0.039 across the 540 null cells, an artefact of the `corCAR1`
working covariance that over-estimates the standard error by 6-10%. A
CR2 cluster-robust standard error restores $\kappa$ to about one and
size to nominal (0.045-0.049), recovers two to five points of power,
and leaves the ranking intact, at the cost of far fewer degrees of
freedom (about 480 down to 23). This was validated only in the
Architecture B / Hybrid / $N = 70$ cells.

## What follows for practice

The widely applied default 'always use the exposure-weighted
predictor' is only conditionally correct. Use A2 when the design has a
discontinuation phase *and* the biomarker is believed to act through a
correlation channel (Architecture B), and report it with a CR2
cluster-robust standard error. Under a crossover design, or when the
biomarker mechanically scales the mean effect (Architecture A), the
simple binary indicator is at least as good and considerably simpler
to defend. A roughly guessed pharmacokinetic half-life is adequate
wherever A2 is used.

The blunt practical message is that none of this competes with
retention: anticipated dropout governs power far more than the
carryover specification does, and should dominate design effort
accordingly.

## Scope and limitations

The numbers come from a single parameter set calibrated to a
prazosin/PTSD trajectory, so absolute power values are illustrative
and only the ranking is portable, and even the ranking is conditional
on architecture and design. The sensitivity blocks vary one factor at
a time and can miss factor interactions; the single joint check (S5)
showed none. The CR2 recalibration is validated in one stratum only.

One caveat internal to the headline table deserves emphasis: bias,
MSE, and coverage there are scored against a single calibrated target
($\theta_{\mathrm{true}} = -3.6$), which is A2's estimand, while A1
and A3 estimate a different coefficient. Those columns are therefore
not comparable across specifications; only power and Type I error
support like-for-like comparison, and the paper says so explicitly.

---
*Rendered on 2026-07-29 at 09:10 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/02-carryover-sensitivity/whitepaper-carryover-specification-summary.md*
