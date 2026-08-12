---
geometry: margin=0.75in
fontsize: 10pt
header-includes:
  - \linespread{0.96}
  - \setlength{\parskip}{0.45em}
---

# Two-page summary

*2026-08-12 11:02 PDT*

**Unadjusted, lag-adjusted, and exposure-weighted carryover
specifications under differential-correlation biomarker
interactions: a factorial simulation study** (pmsimstats team).

## Goal

A companion paper establishes that under covariance moderation, the
multivariate-normal differential-correlation architecture in which
the biomarker-treatment interaction lives in a treatment-state
correlation, carryover (residual drug effect that decays after
discontinuation) costs up to roughly 31% of power in relative
terms. That earlier work, and the original Hendrickson et al.
(2020) publication it builds on, fixed the pharmacokinetic decay
shape at exponential and evaluated only a single way of
representing residual exposure in the analysis model. Given a
trial analyzed under covariance moderation, the analyst must still
decide how to code that residual exposure in the linear predictor,
and it is not established which representation is most robust to
the ordinarily unknown decay shape. This paper takes up that
question directly, restricting attention to covariance moderation
and the exponential and Weibull decay forms.

Three analysis-model specifications for the drug-exposure regressor
are compared. \emph{Unadjusted} is the binary on-drug indicator,
ignoring carryover entirely; it is what an analyst obtains by
default and serves as the benchmark any remedy must beat.
\emph{Lag-adjusted} augments that indicator with a lagged
just-off-drug term that estimates the magnitude of carryover
non-parametrically, following the classical crossover-trial device
of Jones and Kenward. \emph{Exposure-weighted} is a continuous
indicator that commits to a parametric decay form and half-life,
the current default in the pmsimstats-ng package. Unadjusted and
Lag-adjusted estimate the same coefficient and differ by exactly
one nuisance term, so their contrast isolates that term's
contribution alone.

## Method

All three specifications are fitted to the same simulated datasets
under common random numbers, so every pairwise contrast is paired
and estimated with reduced Monte Carlo variance. The principal
comparison is an ADEMP-structured factorial crossing two
data-generating decay forms (exponential; Weibull at three shape
values) with the three specifications, three trial designs
(crossover, Hybrid, OL+BDC), two sample sizes ($N = 35, 70$), and
three biomarker-moderation strengths (216 cells, 500 replicates
each). Four further single-axis sensitivity blocks perturb
within-factor autocorrelation, the analyst's assumed half-life
against the true one, dropout, and biomarker effect size, and a
cluster-robust (CR2) recalibration checks whether the model-based
test's known mild conservatism disturbs the ranking among
specifications.

## Results

Two findings dominate. First, the lagged-treatment term buys
nothing: across all 216 cells, Lag-adjusted and Unadjusted are
statistically indistinguishable, differing by a mean of 0.003 in
power against a Monte Carlo standard error of 0.018, and agreeing
to three decimal places on bias, mean squared error, and coverage.
Because the two specifications share an estimand, are fitted to the
same datasets, and differ by exactly one term, this is a clean
null rather than an underpowered comparison. The mechanism follows
from where the term sits: it enters only as a nuisance main effect,
so it can absorb mean-level carryover but cannot repair the
attenuation of the biomarker interaction itself.

Second, Exposure-weighted attains the highest power only in the
two designs with a blinded-discontinuation phase. At the reference
cell (Hybrid, $N = 70$), its advantage over Unadjusted is about ten
percentage points (0.860 versus 0.766), a difference far beyond
Monte Carlo error. Under the classical crossover design, however,
it is markedly inferior (0.488 versus 0.830), because the
within-subject correlation there already carries most of the
signal the decaying predictor was designed to recover, so its
down-weighting of off-drug observations works against it. Where
Exposure-weighted leads, the lead is robust to incorrect
decay-shape and half-life assumptions and to autocorrelation
strength. The model-based interaction test is mildly conservative
for all three specifications (an artifact of the working
covariance, common to both architectures), and a CR2 cluster-robust
standard error restores nominal Type I error without disturbing
the ranking among specifications. Across every condition examined,
anticipated dropout reduces power far more than the choice among
specifications: at 30% dropout, all three specifications collapse
to between 0.12 and 0.25, essentially indistinguishable.

## Conclusion

Under covariance moderation, we recommend the exposure-weighted
predictor, reported with a cluster-robust standard error, for the
Hybrid and OL+BDC designs; under a crossover design it should not
be used, and the unadjusted binary indicator is materially better.
The classical lagged-treatment remedy cannot be recommended in any
design examined here, since it never improves on ignoring carryover
altogether. An approximate half-life is adequate where exposure
weighting is used, with the error asymmetric in the safer
direction (under-assuming). Preventing dropout is a more effective
use of design effort than any choice among these three
specifications.
