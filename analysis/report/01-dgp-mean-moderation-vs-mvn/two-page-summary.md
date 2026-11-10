---
geometry: margin=0.75in
fontsize: 10pt
header-includes:
  - \linespread{0.96}
  - \setlength{\parskip}{0.45em}
---

# Two-page summary

*2026-11-09 17:58 PST*

**Toward an optimal data-generating process for simulating
biomarker-treatment interactions in aggregated N-of-1 trials**
(pmsimstats team).

## Goal

A predictive biomarker is validated by an interaction test. With two
equal subgroups and an interaction as large as the main effect, that
test needs about four times the sample of a main-effect test
(Brookes et al. 2004). Aggregated N-of-1 designs have been proposed
to make it affordable. Their power is estimated by simulation, and a
simulation can only show what its data-generating process (DGP)
allows. Paper 01 asks what that DGP should be.

## Approach

We start from the one DGP we found that encodes a predictive
biomarker, the multivariate normal model of Hendrickson et al.
(2020). We call it the starting DGP. In it the biomarker acts as a
correlation with the drug response. We take the elements of the DGP
one at a time, compare the candidate choices for each, and judge them
by five criteria:

1. existence (the joint distribution must be valid);
2. identifiability under carryover;
3. fidelity to the pharmacology;
4. analyzability;
5. calibration status.

The choices are collected as the recommended DGP, a step toward an
optimal one. We derive the closed-form power of a paired-difference
interaction statistic (E9) and a Schur-complement bound on the
largest feasible interaction. We check both by Monte Carlo
simulation. Three simulation studies verify the closed forms, rerun
the starting code, and stress-test ten analyses under a prognostic
biomarker. The Hybrid design is the primary case, and all
simulations use the starting response and baseline parameters.

## Results

**Feasibility.** On the grid of the earlier study, the starting code
repaired 20 correlation matrices (25%) that were not positive
definite, without a warning. The recommended DGP refuses such a
matrix instead. On the extended grid it finds 75% of cells infeasible
at an interaction correlation of $c_{bm} = 0.6$. This choice rests on
the existence criterion.

**The biomarker-response correlation under carryover.** This is the
element that matters most. Under the constant rule, the correlation
stays at its on-drug value after first exposure. The interaction is
then absent among the post-baseline visits, and power falls nearly to
the test size at every positive half-life. Under the decaying rule,
the correlation decays at the carryover half-life. Power then stays
near its no-carryover value at pharmacological half-lives of 0.1 to
0.2 week and declines smoothly at longer ones. Power depends on the
rule only through the on/off correlation difference. The starting
code shows the same pattern. At 0.1 week its power was 0.09 with the
constant correlation and 0.64 with the decaying one, against 0.67
without carryover. The choice of the decaying rule rests on
identifiability under carryover.

**Carryover on the mean.** The starting code compounds the decay of
the drug-response mean across off-drug visits. The recommended DGP
anchors it at the last on-drug mean. The choice rests on fidelity to
the pharmacology, as a judgment; its effect on power is at most 0.006
either way.

**Serial correlation.** The correlation structure bounds the largest
interaction the DGP can hold: 0.256 under compound symmetry and 0.474
under separable AR(1) at a serial correlation of $\rho = 0.7$. Across
$\rho$ from 0.4 to 0.95 the bound is largest at 0.7 and falls steeply
above it. The value 0.7 is a convention. A PTSD symptom series
suggests 0.9, where the bound falls to 0.34 and the Hybrid design
needs 154 participants at $c_{bm} = 0.25$ instead of 278.

**Prognostic biomarker.** If a prognostic biomarker acts through a
pathway absent at baseline, analyses that share its effect with
baseline are anticonservative. A prognostic correlation also shares
the biomarker row with $c_{bm}$, so each lowers the feasible room for
the other.

**Planning.** We give the total sample size for 80% power of E9 and
the simulated power of the recommended analysis G10 (a
random-intercept CAR(1) mixed model with an analysis half-life chosen by AIC
on the model without the interaction, and CR2 standard errors). G10
rejected the null at 0.046 to 0.053. At the E9 sample size, G10 power
was 0.844 to 0.936 in the 20 uncapped cells, so the closed-form rule
is conservative for G10. In the Hybrid design the G10 sample size for
80% power is 0.71 to 0.83 times the E9 value. At $\rho = 0.9$, 8 to
12% of G10 fits failed in each crossover cell at $N = 70$. The fitted
serial correlation reached its limit of 1 and could not be told apart
from the random intercept.

## The recommended DGP

The recommended DGP has a biomarker-response correlation that decays
after discontinuation, carryover on the mean anchored at the last
on-drug value, a separable AR(1) structure at $\rho = 0.7$ with a
cross-component correlation of 0.2, and a check that refuses
infeasible effect sizes. Path allocation is unchanged from the
starting DGP. Every parameter is labeled as a starting value, a
choice or uncalibrated.

## Conclusion

Three points are to be emphasized. First, the biomarker-response
correlation must decay after discontinuation; a constant correlation
removes the interaction under any carryover. Second, the serial
correlation structure sets the largest interaction the DGP can hold,
and the required sample size depends strongly on $\rho$ above 0.7.
Third, a matrix outside the feasibility bound should be refused
rather than repaired. We propose the recommended DGP as a reference
for simulation studies of these designs. Calibrating the serial
correlation, the drug-response SD and the symptom half-life is the
next step toward an optimal DGP.
