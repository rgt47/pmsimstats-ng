---
geometry: margin=0.85in
fontsize: 10pt
---

# White paper: should biomarker moderation be encoded as a latent class?

*2026-07-29 09:14 PDT*

**A two-page summary of 'Latent-class and mixture-model formulations
for biomarker-treatment interaction in N-of-1 trials' (pmsimstats
team, 2026-06-13).**

## Goal of the paper

Aggregated N-of-1 trials exist to detect a biomarker-by-treatment
interaction. Published practice encodes that interaction one way only:
as a smooth continuous term, $\beta_{bm:D} \cdot B \cdot D$, in a
single-component regression. Every participant is assumed to respond
to some degree, and the biomarker is assumed to scale that response
linearly.

There is a second, equally defensible encoding. Participants may
belong to unobserved responder and non-responder classes with
qualitatively distinct response curves, with the biomarker predicting
class membership without determining it. This is a finite mixture, and
the psychometric, econometric, and marketing-science literatures have
developed mature machinery for it (latent class and profile analysis,
growth mixture models, factor mixture models, regression mixtures and
mixtures of experts) that the N-of-1 biostatistics literature has not
engaged.

The paper's claim is that the choice between the two encodings is a
substantive scientific question, not a parameterization convenience,
and it pursues four aims:

1. Formalize the latent-class data-generating process for aggregated
   N-of-1 trials and derive when a single multivariate normal can
   stand in for it.
2. Characterize estimation and the identifiability conditions at
   trial-relevant sample sizes, $N \in [30, 150]$.
3. Benchmark, by simulation, whether a class-aware analysis beats the
   standard linear mixed model.
4. Re-analyze the Hendrickson et al. (2020) prazosin-PTSD data under
   the latent-class formulation.

Aims 1 and 2 are delivered in full. Aim 3 is delivered only as a
240-replicate pilot of the last of five pre-registered studies. Aim 4
has not been carried out.

### Why the encoding matters under carryover

The two encodings diverge precisely where N-of-1 trials are weakest.
Carryover enters as a decaying drug indicator $D_{it} =
\exp(-\lambda t_{sd})$. Under continuous heterogeneity the moderation
enters as a product with $D_{it}$, so carryover attenuates the signal
directly, which is the mechanism behind the 40-60% covariance-only
power loss reported in the companion work. Under a latent-class DGP
the class label is fixed across drug states, so some moderation
information survives washout and is in principle recoverable from the
full trajectory. The gap between the two power curves upper-bounds
what a class-aware analysis could recover. The paper is explicit that
this bound is asymptotic and can shrink or reverse at finite $N$.

The theory section also establishes when the MVN shortcut is safe. A
two-component mixture is reproduced to second order by a single MVN
with $c_{bm}^2$ equal to the product of the between-class variance
fractions in $B$ and in $BR$, so the achievable marginal correlation
is bounded by the geometric mean of the two class separations. The
approximation fails in three named regimes: bimodal response
distributions, a near-step gating function (where a thresholded
mean-moderation model is the better approximation), and class-varying
covariance.

## Results

**The pilot.** Study 5, 240 replicates per cell, $N = 70$, eight
timepoints, a heterogeneous-random-slopes DGP with logistic gating in
the standardized biomarker at gating slopes $\{0, 1.0, 1.5\}$, run on
four cores in 558 s. Three tests per replicate: the linear-mixed
`bm:Dbc` Wald, an unconditional Wald on the `lcmm::hlme` gating
coefficient, and a BIC-conditional gating procedure that tests the
gating coefficient only if BIC prefers $ng = 2$. Monte Carlo standard
error on power is about 0.030; `lcmm` convergence was 92.9% at the
null and 98.3-98.8% elsewhere.

**The unconditional class-aware test is badly anti-conservative.**
Type I error 0.197 (95% CI [0.145, 0.250]) against nominal 0.05, more
than five MCSE of separation. The diagnosis is a finite-sample
identifiability artifact, not bias: the gating coefficient is centered
at zero, but its empirical standard deviation across replicates (63.9)
is roughly fifty times the median within-replicate Wald standard error
(0.86). Forced to fit two classes to data with none, the optimizer
latches onto random class boundaries and the Wald standard error does
not absorb the labeling uncertainty. This is exactly the overextraction
hazard predicted by Bauer and Curran (2003) and Lubke and Muthen
(2005). The naive likelihood-ratio test is anti-conservative for a
different reason, rejecting at 0.146, because the true reference
distribution at the $ng = 1$ boundary is chi-bar-squared.

**The BIC-conditional test controls size but is inert.** Type I error
0.004 under the null, but power of only 0.025 at gating slope 1.0 and
0.004 at 1.5. At this $N$ and DGP the BIC step selects a single class
in essentially every replicate, so the procedure approximates a
never-reject rule. Power is non-monotone in the gating slope, which
the paper reads as the single-class LME absorbing the
heterogeneous-slopes signal well enough that BIC never prefers two
classes.

**The standard analysis dominates.** The linear-mixed `bm:Dbc` Wald
test held nominal Type I (0 rejections of 240 under the null) and
reached 0.542 power at gating slope 1.5, against 0.360 for the
unconditional Wald and 0.004 for the BIC-conditional procedure. The
unconditional test buys nothing for its size inflation: its power is
below the LME baseline at both non-null cells.

**A prior claim is corrected.** The paper's own earlier proposition,
that class-aware analyses are unusable as primary tests at small $N$,
is narrowed. It holds for an unconditional Wald test on a single
mixture coefficient, which requires a model-comparison wrapper (BIC or
a bootstrapped LRT) to control size. The BIC-conditional procedure is
size-valid; its problem is sensitivity, not validity.

## What follows for practice

At $N \in [30, 150]$, a mixture model should not be the primary
analysis for the biomarker-treatment interaction. Its legitimate roles
are as a sensitivity check and as a DGP for stress-testing the
pre-specified linear mixed model. If the simulation infrastructure is
extended, the heterogeneous-random-slopes form is the right target,
because it admits latent-class generative structure while leaving the
analysis model a standard `nlme::lme` fit and preserving the audit
chain.

The paper also warns that a sensitivity analysis run only under the
covariance-only DGP will overstate carryover-induced power loss, since
that specification puts the entire signal in the one channel carryover
erodes by construction. Any signal flowing through the mean structure
survives. Report power at three points on the spectrum: pure mean
moderation, pure covariance moderation, and combined.

## Scope and limitations

This is predominantly a theory and taxonomy paper with a pilot
attached, and the paper says so. The five-study production programme
(approximately 2,300 pre-registered ADEMP cells; estimated at about a
week of background compute on 8-16 cores) has not been run, and the
prazosin-PTSD re-analysis has not been attempted.

The empirical claims therefore rest on 240 replicates at a single
sample size ($N = 70$), a single design, and one generative
mechanism. The conclusion is a non-competitiveness result in that
cell, not a general one: the regime where a class-aware test should
win, large class separation with genuinely bimodal response, was never
tested. Indeed the pilot's own DGP is one where the single-class model
performs well, which is the paper's stated explanation for BIC's
inertia and is a reason to expect the finding to be regime-specific.
A canonical 1,000-replicate run at $N \in \{35, 70, 100\}$ over a wider
gating-slope grid is required before the negative result can be
generalized.

Two further cautions. The class-aware analysis costs roughly 32 times
the linear-mixed fit per replicate, which is what confines the pilot
to 240 replicates and rules out the bootstrapped LRT that would
correctly calibrate the likelihood-ratio test. And the reported Type I
error of 0.000 for the LME baseline is a resolution artifact of 240
replicates rather than evidence of exact size; it bounds the rate
below roughly 0.015, no tighter.

---
*Rendered on 2026-07-29 at 09:14 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/03-latent-class-mixture-application/whitepaper-latent-class-mixture-summary.md*
