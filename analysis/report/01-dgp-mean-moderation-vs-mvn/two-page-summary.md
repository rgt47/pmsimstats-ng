---
geometry: margin=0.75in
fontsize: 10pt
header-includes:
  - \linespread{0.96}
  - \setlength{\parskip}{0.45em}
---

# Two-page summary

*2026-08-12 11:02 PDT*

**Mean moderation and covariance moderation as data-generating
architectures for biomarker-treatment interactions in aggregated
N-of-1 trials, and their divergent power under carryover**
(pmsimstats team).

## Goal

Simulation is the standard tool for power and sample-size
calculation in aggregated N-of-1 trials, designs in which each
participant is repeatedly crossed on and off treatment and many
single-patient series are pooled to validate a predictive biomarker.
Any such simulation must specify how the biomarker moderates the
treatment effect in its data-generating process (DGP), and that
specification is a substantive scientific choice rather than an
implementation detail. Each of the two established parameterizations
carries its own advantages and costs. The paper's goal is to
formalize both, and to show that the choice between them materially
changes the estimated power loss under carryover.

Under \emph{mean moderation}, the biomarker scales the treatment
effect in the conditional mean: a deterministic, proportional
relationship appropriate when the biomarker directly determines
drug pharmacokinetics or pharmacodynamics (for example, kidney
function setting the clearance of a renally cleared drug). Under
\emph{covariance moderation}, the interaction is instead
represented as a correlation between biomarker and biological
response that depends on drug-exposure condition (on drug versus
off drug): a probabilistic relationship appropriate when the
biomarker is only an imperfect proxy for a latent responder
subtype. Covariance moderation is often the more biologically
defensible choice in that latter setting, and it is the
specification Hendrickson et al. (2020) selected for the paper's
motivating application, a prazosin trial for PTSD in which resting
blood pressure is thought to index noradrenergic tone rather than
directly drive drug levels. That defensibility, the paper argues,
does not come free: because the covariance-moderation signal lives
in a correlation that is sustained only while the drug is active,
any off-drug erosion of drug effect directly erodes the statistical
signal itself.

## Method

Both architectures are embedded in an otherwise identical Monte
Carlo simulation, holding the analysis model fixed at a linear
mixed model with continuous-time AR(1) residual correlation. Power
to detect the biomarker-by-treatment interaction is estimated as
the carryover half-life increases from zero to one week, across
three within-subject designs matched to the prazosin and PTSD
reference application: a classical crossover, a Hybrid N-of-1
design, and an open-label design followed by blinded discontinuation
(OL+BDC). Sample size is $N = 70$ throughout, matched across
designs so the comparison is a design comparison, not a
sample-size one. The study follows the ADEMP reporting framework,
with 1,000 replicates per cell and every performance measure
reported with its Monte Carlo standard error.

## Results

Under mean moderation, power is nearly preserved as carryover
increases, because carryover does not disrupt the biomarker's
proportional relationship to the treatment effect: in the OL+BDC
design, power falls only from 0.751 to 0.730 (a 2.8% relative
loss), and losses across all three designs range from 1.9% to
6.9%. Under covariance moderation, the same carryover erodes the
differential correlation that carries the signal, and power falls
much further, from 0.777 to 0.539 in OL+BDC (a 30.6% relative
loss), roughly eleven times the matched mean-moderation loss and
far beyond Monte Carlo error ($z = 7.64$, $p < 10^{-13}$).

The gap between the two architectures is strongly design-dependent.
It is largest in designs with a long blinded-discontinuation phase
(OL+BDC, then Hybrid at 15.4% loss) and negligible in the classical
crossover (3.0%, not statistically distinguishable from zero),
because the within-subject AR(1) correlation in a crossover design
already carries most of the signal that covariance moderation
depends on, leaving carryover little left to erode. A robustness
check at $N = 140$ on the OL+BDC design confirms that this ordering
survives a larger sample, though absolute magnitudes are compressed
by ceiling saturation. Type I error across all null cells sat
uniformly below the nominal 5% level for both architectures alike,
a known conservatism of the working covariance structure that does
not disturb the architecture comparison.

## Conclusion

The DGP architecture is a substantive scientific assumption, not a
parameterization convenience, and it determines how severely
carryover appears to reduce power. The choice should be guided by
the hypothesized biological mechanism: mean moderation for
biomarkers that directly determine drug pharmacokinetics or
pharmacodynamics, covariance moderation for biomarkers that
statistically predict treatment response through latent subtype
membership. The existing clinical trial simulation literature uses
mean moderation almost exclusively, so published power estimates
for biomarkers that actually operate through a correlation-based
mechanism may be systematically optimistic. Biomarker-validation
simulation studies should state which architecture they adopt and
report power under both.
