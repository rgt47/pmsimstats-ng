---
geometry: margin=0.75in
fontsize: 10pt
header-includes:
  - \linespread{0.96}
  - \setlength{\parskip}{0.45em}
---

# White paper: DGP architecture choice and carryover power loss

*2026-07-29 09:02 PDT*

**Summary of 'Three data-generating architectures for
biomarker-treatment interactions in aggregated N-of-1 trials, and why
the choice changes power under carryover: a tutorial' (pmsimstats
team, 2026-06-17).**

## Goal of the paper

Aggregated N-of-1 trials pool many single-patient on/off series to
validate predictive biomarkers. Power for the biomarker-by-treatment
interaction has no closed form in these designs, so it is obtained by
Monte Carlo simulation. Every such simulation must encode how the
biomarker moderates the treatment effect in its data-generating
process (DGP), and there is more than one way to do so. That encoding
is normally treated as an implementation detail.

The paper's goal is to show that it is not a detail but a substantive
scientific assumption. Specifically, it aims to:

1. Formalize three DGP architectures as co-equal modeling choices.
2. Quantify, by simulation with the analysis model held fixed, how
   much each architecture's estimated power erodes as carryover
   half-life increases, across three within-subject designs.
3. Supply guidance on which architecture matches which biomarker
   biology, and on what simulation studies must report.

The work is motivated by a concrete discrepancy: the original
Hendrickson et al. (2020) code used one encoding, a re-implementation
in the same software project quietly used another, and the two gave
different accounts of how much carryover costs.

### The three architectures

- **Architecture A (direct mean moderation).** The biomarker scales
  the treatment effect in the conditional mean,
  $Y_{it} = \mu_0 + \beta_1 D_{it}(1 + \beta_{bm} B_i) + \epsilon_{it}$.
  The interaction is deterministic and lives in the first moment.
  This is the near-universal convention in the trial-simulation
  literature.
- **Architecture B (differential correlation, MVN).** Biomarker and
  biological response are drawn jointly from a multivariate normal
  whose BM-BR correlation is treatment-state dependent: $c_{bm}$ on
  drug, $c_{bm}e^{-\lambda t_{sd}}$ off drug. The interaction is
  probabilistic and lives in the second moment. This is the
  Hendrickson N-of-1 approach.
- **Architecture C (combined).** Both channels act at once through
  independent parameters $c_{bm,A}$ (mean) and $c_{bm,B}$
  (covariance); A and B are its single-channel limits, enforced
  algebraically rather than by convention.

### Simulation setup

Identical configurations were run under each architecture, differing
only in the BM-BR linkage. Designs: crossover (CO), Hybrid, and
open-label with blinded discontinuation (OL+BDC). Primary run
$N = 35$ per randomization path, $c_{bm} = \beta_{bm} = 0.45$,
DGP carryover half-life $t_{1/2} \in \{0, 0.5, 1.0\}$ weeks with a
matched analysis half-life, 1,000 replicates per cell (Monte Carlo
SE near 0.016), 36,000 fits, 100% convergence. The analysis model
was fixed throughout: `nlme::lme` with `corCAR1` residual
correlation and the interaction term `bm:Dbc`, where `Dbc` equals 1
on drug and $e^{-\lambda t_{sd}}$ off drug.

## Results

**Architecture A preserves power.** Moving from no carryover to a
one-week half-life costs only 1.9% relative power in CO, 1.4% in
Hybrid, and 2.8% in OL+BDC (0.751 to 0.730). The biomarker
proportionality is maintained at every occasion, so the interaction
signal contracts in amplitude without losing identifiability.

**Architecture B loses substantially, and design-dependently.** The
same carryover costs 3.0% in CO (not distinguishable from zero),
5.1% in Hybrid, and 30.6% in OL+BDC, where power falls from 0.777
to 0.539. The OL+BDC loss is roughly eleven times the matched
Architecture A loss.

**The gap is far beyond Monte Carlo error.** A difference-of-
differences $z$-test on the within-architecture losses gives
$z = 7.64$ ($p < 10^{-13}$) at OL+BDC and $z = 3.96$ ($p < 10^{-4}$)
at Hybrid; the CO contrast is null ($z = 0.29$), because the
crossover design's within-subject AR(1) structure already dominates
the correlation channel that Architecture B uses to carry the signal.

**Mechanism.** Both architectures suffer mean blurring, since
carryover compresses the `Dbc` contrast. Architecture B suffers a
second, unique attack: the BM-BR correlation that generates its
second-moment signal decays exponentially off drug. Under
Architecture A the drug-mediated and drug-independent components of
BR variance move together, so their ratio, and hence identifiability,
survives.

**Architecture C is intermediate.** Across the $3 \times 3$ grid of
$c_{bm,A} \times c_{bm,B} \in \{0, 0.22, 0.45\}^2$ at
$t_{1/2} = 1.0$, power increases strictly with each parameter, and
the boundary cells reproduce the pure architectures exactly. The
paper reads the interior cells as super-additive: Hybrid
$(0.22, 0.22)$ reaches 0.299 against an additive-on-the-probability-
scale benchmark of 0.157. At full signal the Hybrid design retains
0.870, exceeding either pure architecture alone.

**Auxiliary checks.** Type I error across the 18 null cells sat in
[0.010, 0.074] (mean 0.043) with no architecture-specific inflation.
An $N = 70$ run on OL+BDC preserved the B-versus-A ordering (9.3
versus 1.7 percentage points lost), compressed by ceiling
saturation.

## What follows for practice

The architecture encodes *why* the biomarker predicts response, so
choosing one is a biological claim. Use Architecture A when the
biomarker mechanically sets the size of the drug effect (renal
clearance, receptor density, metabolizer genotype). Use Architecture
B when the biomarker is only a proxy for a hidden responder state
(resting blood pressure in PTSD, CRP or IL-6, polygenic risk
scores). The motivating prazosin/PTSD application is an Architecture
B story, which makes its carryover sensitivity a clinically relevant
warning rather than an artifact.

Because the literature uses Architecture A almost exclusively,
published power figures may be systematically optimistic for
biomarkers that genuinely operate through a correlation channel. The
distinction bites hardest in treatment-switching designs (crossover,
N-of-1, basket, adaptive enrichment) and is immaterial in
parallel-group designs with no off-drug phase. Recommendations:
state which architecture a simulation adopts and why; when the
mechanism is uncertain, sweep Architecture C across a mixing weight
and report the range; treat Architecture B as the conservative
default; and reduce exposure to the question by designing adequate
washout.

## Does Architecture C earn its place?

Architecture C is worth keeping, but for narrower reasons than the
paper claims, and its current execution understates the case. The
defensible arguments are two. First, a verification argument: the
boundary cells $(c_{bm,A} = 0.45, c_{bm,B} = 0)$ and
$(0, 0.45)$ recover the pure architectures algebraically rather than
by convention, so any C-based simulation carries a built-in
correctness gate on the DGP implementation. Second, a design
argument: when the mechanism is unknown, C with a mixing weight
$\alpha$ is the only construct spanning the A-to-B continuum in one
parameterization, which is what a sensitivity analysis requires.
The power argument, by contrast, is weak as presented. Power is a
convex function of effect size below 0.5, so additivity on the
probability scale is the wrong benchmark: a normal approximation in
which the two channels merely add on the effect-size scale already
predicts 0.254 for the Hybrid $(0.22, 0.22)$ cell against 0.157 for
the probability-scale benchmark and 0.299 observed, leaving a
residual of roughly three Monte Carlo standard errors rather than
the reported gap of nine. Three further caveats constrain what C can
deliver. The $c_{bm,A}$ and $c_{bm,B}$ parameters are not separately
identifiable from a fitted `bm:Dbc` coefficient, so C is a
DGP-side construct that can be simulated but never estimated from
trial data, which makes the recommendation to obtain independent
evidence on each channel difficult to act on. The C grid uses
$N = 35$ as a total rather than per path, so its cells are not
readable against the headline tables. And the carryover comparison
within the grid is reported only for the crossover design, which is
precisely the design in which Section 3 found no
architecture-by-carryover interaction; the informative panel would
be OL+BDC. Net assessment: retain C as a scoped sensitivity device
with the effect-size-matched parameterization
($c_{bm,A} = (1-\alpha)c$, $c_{bm,B} = \alpha c$) that the paper
defines but declines to use, rather than presenting it as a
co-equal third architecture.

## Scope and limitations

The paper is a tutorial and simulation study, not an empirical
validation against trial data. Section 6's power-recovery strategies
(on-drug-only analysis, exposure weighting, within-subject
contrasts, two-stage random slopes) are qualitative hypotheses with
no simulation support here; a companion paper addresses them.
Architecture B is acknowledged as a correlation-level approximation
to what would more faithfully be a finite-mixture responder model.
The effect-size-scale recomputation in the preceding section is a
normal-approximation check performed for this summary, not a
re-simulation.

\vspace{0.5em}
*Rendered on 2026-07-29 at 09:12 PDT. Source:
~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/01-dgp-mean-moderation-vs-mvn/whitepaper-dgp-architectures-summary.md*
