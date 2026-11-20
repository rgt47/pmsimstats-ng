---
geometry: margin=0.75in
fontsize: 10pt
header-includes:
  - \linespread{0.96}
  - \setlength{\parskip}{0.45em}
---

# Two-page summary

*2026-11-19 19:59 PST*

**Testing a predictive biomarker under drug carryover in aggregated
N-of-1 trials: a staged simulation study** (pmsimstats team).

## Goal

A predictive biomarker is validated by a biomarker-by-treatment
interaction test. In an aggregated N-of-1 trial that test is made
within patients, from repeated on-drug and off-drug contrasts. Drug
carryover weakens the contrast. Observations taken after
discontinuation still carry part of the drug effect, and a binary
on/off regressor codes them as if they carried none. The analyst must
therefore decide three things. These are how to represent residual
exposure, how to size the test, and whether the half-life must be
known. Paper 02 asks which analysis answers all three.

## Approach

We simulate under the structure of the data-generating process of the
companion paper (paper 01). There the interaction is a correlation
$c_{bm}$ between the biomarker and the drug response that is present
on drug and decays off drug at the carryover half-life. The parameter
values are the starting values of Hendrickson et al., which paper 01
also uses (parameter set 1). The whole study was repeated under a
second set with a weaker placebo pathway (parameter set 2). The
biomarker is predictive only. We compare five analyses arranged as a
ladder. Each rung changes one thing about the analysis, except the
step from rung 0 to rung 1, which changes three.

- **A0**, the starting analysis of Hendrickson et al. (2020): random
  intercept, binary drug indicator, no residual correlation.
- **A1**: the binary indicator with `corCAR1` residual correlation
  and the working-covariance standard error.
- **A2**: A1 with a CR2 cluster-robust standard error and
  Satterthwaite degrees of freedom.
- **A3**: an exposure-weighted predictor at the true half-life, with
  CR2. No analyst knows the true half-life, so A3 is an oracle
  benchmark.
- **A4**: the exposure-weighted predictor at the half-life chosen by
  AIC over four candidates. The selection uses the model without the
  interaction. The interaction model is then refitted by REML and
  tested with CR2. Papers 01 and 15 call this analysis G10.

The ladder is compared in six scenario families. Family 1 is the base
grid of 30 cells. It covers three designs (crossover, OL+BDC and
Hybrid) and $N$ of 35 and 70. It sets $c_{bm}$ to 0, 0.30 and 0.45 at
a 1-week half-life, and the half-life to 0 and 0.5 week at
$c_{bm} = 0.45$. Families 2 to 5 add a misspecified half-life and
decay shape, serial correlation of 0.5 and 0.9, patient dropout, and
a patient slope variance. Family S is a Hybrid run with independent
seeds for every cell. Within a cell every analysis is fitted to the
same trials, so the comparisons are paired. A separate factorial of
exponential and Weibull decay compares the binary indicator with the
exposure-weighted predictor (rung 3). The study is exploratory. Its
rungs and scenarios were not fixed in advance.

## Results

**Rungs 0 and 1: residual correlation.** Pooled over the three
designs in the base grid, A0 rejected the null at 0.107. It is liberal
in every design, from 0.074 to 0.076 under Hybrid to 0.128 to 0.150
under the crossover design. Isolation analyses attribute this to its
missing residual correlation. A1 rejected at 0.038. That is, the
working-covariance test is conservative. A0 departs from the analysis
as published in two respects: it adds an expectancy term and it has
no baseline row. The analysis exactly as published rejected at 0.097,
and at 0.045 once `corCAR1` was added. Its baseline row, an off-drug
outcome at time zero, carried 13 to 17 points of power. No rung of
the ladder uses that row, so every rung's absolute power is probably
understated.

**Rung 2: cluster-robust inference.** CR2 moved the pooled size to
0.054. It brings every analysis to within about half a point of
nominal. It changed power by a mean of $+3.4$ points under Hybrid,
$+4.9$ under OL+BDC and $-3.4$ under the crossover design.

**Rung 3: exposure weighting.** At a 1-week half-life the
exposure-weighted predictor at the true half-life recovers 20% to 22%
of the power that carryover removes under Hybrid, and about 12%
under OL+BDC. The rest has to be met by the design, for example with
more patients. Where it does not gain, it is within Monte Carlo error
of the binary indicator.

**Rung 4: choosing the half-life.** Paired against the oracle, AIC
selection changed size by $-0.005$ to $+0.005$ and power by $-0.011$
to $+0.002$ across the six families. The selected half-life is spread
over the candidates. The small cost therefore reflects how little the
predictor depends on the half-life, not recovery of the true value.
Leaving the interaction out of the selection matters. Selecting with
it in the model raised size in every family. In the run with
independent seeds, the size of A4 minus that of the oracle was
$+0.0008$ (90% interval $-0.0020$ to $+0.0035$). That is equivalent
within $\pm 0.01$. The same equivalence held in 2 of the other five
families; in the other three the interval reached just past the
margin.

**A4 against the binary indicator.** In the Hybrid design A4 exceeded
A2 by only 0.4 points on average over 40 alternative cells, and by up
to 4.1 points at long half-lives. It was below A2 in 13 cells, by at
most 1.5 points. Under dropout of 10% to 30% of patients, completely at
random or dependent on baseline severity, A4 kept its size. In all but
one cell it also kept its power advantage.

**Sample size.** Halving the effect roughly quadruples the required
$N$. Under parameter set 1, a Hybrid trial with a 1-week half-life
analyzed with A4 needs about 154, 273 and 613 patients at $c_{bm}$ of
0.4, 0.3 and 0.2 for 80% power. The binary indicator needs about 15%
to 18% more. A4 held its size at every $N$ checked, up to 1,012.
Under parameter set 2 the requirement is about half as large.

## Recommendation

For a Hybrid trial with a biomarker known to be predictive only, use
A4 to test the interaction. It is correctly sized and needs no
half-life. Three qualifications apply. First, A4 is for the test, not
for an estimate. Its estimate is biased when the selected half-life is
wrong, and its interval covers below the nominal rate. Second, under
the crossover design use A2 instead. Third, whether a biomarker is
also prognostic is rarely known before a trial. A biomarker that
predicts the natural-history trend biases A4 as fitted here. For a
biomarker of unknown role, the default is the analysis of the
companion paper on prognostic biomarkers (paper 15). This paper
supports the parts of it that A4 shares: the exposure-weighted
predictor, the choice of half-life without the interaction, and the
CR2 test.

## Limits

The results are conditional on the correlation encoding, on two
parameter sets that gave the same ranking, on one design and on a
predictive-only biomarker. Size was evaluated at $N$ of 35 and 70 for every analysis,
and up to 1,012 for A2 and A4. The analysis model is close to the
conditional mean of the data-generating process. Part of what the
ladder shows is therefore that a correctly specified regressor
outperforms a misspecified one. Even the oracle recovers only part of
the power that carryover removes. The rest has to be met by the
design, with more patients and an effort to retain them.
