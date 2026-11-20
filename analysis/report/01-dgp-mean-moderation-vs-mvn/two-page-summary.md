---
geometry: margin=0.75in
fontsize: 10pt
header-includes:
  - \linespread{0.96}
  - \setlength{\parskip}{0.45em}
---

# Two-page summary

*2026-11-19 20:12 PST*

**Toward an optimal data-generating process for simulating
predictive-biomarker interactions in aggregated N-of-1 trials**
(pmsimstats team).

## Goal

A predictive biomarker is validated by an interaction test. In an
aggregated N-of-1 design every participant is observed on and off
drug, so the test can be made within participants. Its power is
estimated by simulation. A simulation can show only what its
data-generating process (DGP) allows. Paper 01 asks what that DGP
should be.

## Approach

The DGP is a multivariate normal model. The baseline biomarker and
three latent response components (drug response, placebo response,
natural history) at every visit form one vector. The design assembles
its correlation matrix. The predictive biomarker enters as a
correlation $c_{bm}$ with the drug response at on-drug visits
(Hendrickson et al. 2020). We take the elements of the DGP one at a
time and compare candidate choices for each. We judge them against
three requirements.

1. **A valid joint distribution.** The correlation matrix must be
   positive definite.
2. **An interaction that survives carryover.** The
   biomarker-by-treatment interaction must remain a contrast between
   on-drug and off-drug visits.
3. **Flexibility for stress conditions.** Dropout, measurement error
   and a different decay form are available as options, off by
   default.

Pharmacological plausibility breaks ties. Every parameter is labeled
as a default value, a choice made here, or uncalibrated. We derive the
closed-form power of a paired-difference interaction statistic (E9).
We also derive a Schur-complement bound on the largest feasible
interaction, the feasibility bound. Study 1 checks both by Monte Carlo
simulation. Study 2 stress-tests ten analyses under the optional
features. The biomarker is assumed predictive only. That is, it is
correlated with the drug response and with nothing else (the
prognostic case is paper 15). The Hybrid design is the primary case.
All simulations use the default parameter values of Hendrickson et
al.

## Results

**Infeasible matrices.** Above the feasibility bound the correlation
matrix is not positive definite. A DGP can repair such a matrix or
refuse it. A repaired matrix simulates a different DGP from the one
stated. The recommended DGP therefore refuses it (requirement 1).

**The biomarker-response correlation under carryover.** This is the
element that matters most. Under the constant rule the correlation
stays at its on-drug value after first exposure. The interaction is
then absent among the post-baseline visits. Power falls nearly to the
test size at every positive half-life. Under the decaying rule the
correlation decays at the carryover half-life. Power then stays near
its no-carryover value at pharmacological half-lives of 0.1 to 0.2
week. It declines smoothly at longer ones. Power depends on the rule
only through the on/off correlation difference (requirement 2).

**Carryover on the mean.** We decay the drug-response mean from the
last on-drug value rather than from the previous visit. The half-life
then means what it says. The effect of this choice on power is at
most 0.006.

**Serial correlation and cross-component structure.** The correlation
structure bounds the largest interaction the DGP can hold. The bound
is 0.256 under compound symmetry and 0.474 under separable AR(1) at a
serial correlation of $\rho = 0.7$. The separable form is positive
definite for any visit schedule and gives an exact bound. The value
0.7 is a convention. A PTSD symptom series suggests 0.9. There the
bound falls to 0.34, and the Hybrid design needs 154 participants at
$c_{bm} = 0.25$ instead of 278.

**Optional features.** In study 2, measurement error lowered power by
0.05 to 0.09 and dropout of 20% by 0.01 to 0.06. A Weibull decay moved
it by up to 0.10. The planning sample sizes assume none of these
features.

**Planning.** At the pharmacological half-life the Hybrid design
needs about 280 participants at $c_{bm} = 0.25$ and 84 at 0.45 for
80% power of E9 in closed form. At a 1-week half-life it needs 480 and
146. These are exact values for the paired-difference statistic.
Paper 02 gives sample sizes by simulation for the mixed-model analysis
a trial would use. Paper 02 calls that analysis A4, and this paper
calls it G10. Paper 02 uses the same parameter values. In its runs A4
had at least 80% power at the closed-form sample size wherever that
size was at most 400, so the closed form is conservative for A4.

**Differences from the DGP of Hendrickson et al.** A separate section
lists every difference and reruns their code. As published, power
fell from 0.668 without carryover to 0.126 at a 0.1-week half-life. A
4,000-replicate rerun gave 0.090. With only the correlation lines
replaced by the decaying rule, power was 0.638. In their grid, 20 of
54 matrices at $c_{bm} \in \{0.3, 0.6\}$ were repaired without a
warning.

## The recommended DGP

In the recommended DGP the biomarker-response correlation decays
after discontinuation.
Carryover on the mean decays from the last on-drug value. The serial
correlation is separable AR(1) at $\rho = 0.7$, with a
cross-component correlation of 0.2. A check refuses infeasible effect
sizes. Time on drug restarts at re-exposure. Measurement error, the
decay form and dropout are options. Every parameter is labeled as a
default value, a choice or uncalibrated.

## Conclusion

First, the
biomarker-response correlation must decay after discontinuation. A
constant correlation removes the interaction under any carryover.
Second, the serial correlation structure sets the largest interaction
the DGP can hold. The required sample size depends strongly on $\rho$.
Third, a matrix outside the feasibility bound should be refused rather
than repaired. We propose the recommended DGP as a reference for
simulation studies of these designs when the biomarker is predictive
only. Calibrating the serial correlation, the drug-response SD and the
symptom half-life is the next step toward an optimal DGP.
