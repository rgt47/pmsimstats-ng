# Referee report: paper 02, carryover representation and inference for biomarker-treatment interaction tests (Hybrid design, covariance moderation)
*2026-10-27 11:45 PDT*

Manuscript: `analysis/report/02-carryover-sensitivity/report.Rmd`
(revised 2026-10-27). Reviewed against the standard of a Statistics in
Medicine or Trials methodology article. Evidence labels used below:
VERIFIED (re-derived or executed on the stored data), INSPECTED (read
in the source and judged inconsistent), SUSPECTED (plausible, not
confirmed).

## 1. Summary of the manuscript

The manuscript reports a staged simulation study of how to analyze the
biomarker-by-treatment interaction in an aggregated N-of-1 trial when
drug carryover is present, under the covariance-moderation
data-generating process (the 'refined process' of companion paper 01).
Stage 1 compares a binary on-drug indicator (G1), that indicator plus a
lagged-treatment nuisance term (G2) and an exposure-weighted predictor
at the true half-life (G3) in a 216-cell factorial. Stage 2 replaces
the oracle half-life by an AIC-selected one (G4). Stage 3 attaches a
CR2 cluster-robust standard error with Satterthwaite degrees of
freedom (G6 to G9). Stage 4 compares G9 with the published analysis
(E0) and an expectancy-term family (F1, F1A0B) on scenario cells, and
a gap-fill study adds scenario grids and patient-level dropout.
Stage 5 stresses G9 with four prognostic-biomarker pathways, traces a
power loss at confound strength 0.6 to a silently repaired
non-positive-definite covariance matrix, derives a closed-form
feasibility ceiling, and runs a feasible confound ladder and large-$N$
null cells. The authors recommend G9 (exposure-weighted predictor,
half-life chosen by AIC, CR2 inference) for the Hybrid design.

## 2. Overall assessment and recommendation

Recommendation: MAJOR REVISION.

The question is relevant, the staged logic is easy to follow, the
honest reporting of two code defects and of the repair artifact is
commendable, and the feasibility-ceiling derivation is a sound,
reusable contribution (the Schur-complement condition is checked by
eigenvalues on either side of each ceiling in the stored data). The
core qualitative findings are plausible and largely reproduce from the
stored summaries: the lagged term is inert, exposure weighting gains
power under Hybrid and OL+BDC, and CR2 repairs the conservative
model-based test. However, the manuscript's central size claim for the
recommended method is contradicted by its own per-replicate data; the
description of the selection step does not match the code; the
evidence for 'robustness' rests on 500-replicate cells that share
random-number streams; several counts and ranges in the abstract and
text are wrong; the prognostic stress test contains a vacuous pathway
and omits the pathway most likely to bias the estimator; and the
comparisons conflate size and power. None of these is fatal, and most
can be repaired with reanalysis of existing per-replicate outputs plus
a modest number of additional runs, but together they mean that the
headline statements ('nominal size including after half-life
selection', 'at or above G1 and F1 at every alternative but one',
'liberal at every null cell') cannot stand as written.

## 3. Significance and novelty

What is new, stated precisely: (i) an evaluation, under one DGP, of
carryover representations for the interaction test (not the main
effect) in the Hybrid aggregated N-of-1 design; (ii) the combination
of an exponentially decayed exposure regressor, AIC half-life
selection on a four-point grid, and CR2/Satterthwaite inference; (iii)
a demonstration that the published analysis (no residual correlation)
is liberal under this DGP; (iv) the closed-form PD feasibility ceiling
for prognostic couplings. Item (iv) is the most transferable.

Competing and neighboring literature the paper should engage with:

- Gartner, Schneider, Arnrich and Konigorski (BMC Med Res Methodol
  2023;23:191) compare analyses of aggregated N-of-1 trials with
  carryover and introduce a carry-over adjusted parametric model with
  exponential wash-in and wash-out terms. This is the closest prior art
  to an exponentially decayed exposure regressor and is not cited. Its
  estimand is the treatment effect, not the biomarker interaction, and
  that distinction is where the present paper's contribution lies, but
  the paper must say so.
- Percha et al. (J Med Internet Res 2019;21:e12641) simulate N-of-1
  designs with carryover and wash-in and give design recommendations.
- Blackston et al. (Healthcare 2019;7:137) show inflated Type I error
  for aggregated N-of-1 analyses that ignore moderate-to-strong
  carryover, directly relevant to the E0 and G1 findings.
- Chen and Chen (PLoS ONE 2014;9:e87752) compare four N-of-1 analyses
  (currently cited under the wrong author names, comment M9).
- Small-sample cluster-robust inference: Pustejovsky and Tipton (J Bus
  Econ Stat 2018;36:672-683) is the reference for CR2 with
  Satterthwaite df as implemented in `clubSandwich`; it is not cited.
  Kenward and Roger (Biometrics 1997) is the standard mixed-model
  small-sample alternative and a script for it exists in the repository
  but is not mentioned.
- Post-selection inference: Leeb and Potscher (Econometric Theory 2005)
  and Kuchibhotla, Kolassa and Kuffner (Annu Rev Stat Appl 2022) show
  why a post-selection test that is correctly sized under one
  configuration need not be uniformly valid.

Fit. For Trials (methodology article) the topic fits, subject to the
format points in comment m14. For a top-tier statistics journal the
statistical novelty is modest: the estimators are standard components
assembled and compared on one DGP, and the generality of the
conclusions is limited by comments M3, M5 and M7. Contemporariness:
the reference list contains no work after 2022 and omits the
2019 to 2023 simulation literature on aggregated N-of-1 analysis cited
above.

## 4. Major comments

Correctness problems come first.

**M1. The half-life selection is described as using the null model; the
code uses the model with the interaction (INSPECTED, formula read in
source).** The Abstract (Conclusions), Section 2.4.4 ('F1A0 selects it
on the null model and refits with the interaction, as G4 does') and
Section 4.2 ('ML-AIC ... on the null model') state that G4 and G9
select the half-life on the null model. In `simulation-core.R`,
`select_half_life_aic()` fits `Sx ~ bm + t + Dbc + bm:Dbc` at every
candidate and takes the minimum AIC, so the interaction is in the
selection model. In the project's own F-family vocabulary this is the
F1A construct, not F1A0 (`FM0_F` omits `bm:Dbc`). The distinction is
not cosmetic: the project CLAUDE.md records that selecting by AIC with
the interaction present and testing as if the half-life were fixed is
liberal, and the selected half-life then depends on the very
coefficient being tested. Remedy: correct the text, or change the code
to select on the model without `bm:Dbc` and rerun; ideally report both
variants (the null-model variant is the one with a clean argument for
approximate independence of selection and test under the null).

**M2. The claim that selection 'costs nothing detectable in size' is
contradicted by the per-replicate data (VERIFIED, executed on
`g9-gap-fill/job-31/all.rds`).** Pooling the four 2,000-replicate null
cells (Table 4 in the manuscript), G9 rejects in 0.0555 of 7,987 fits
(MCSE 0.0024, about 2.3 MCSE above 0.05), whereas G8, the same
predictor at the oracle half-life, rejects in 0.0493. In the paired
data G9 rejects where G8 does not in 77 replicates and the reverse in
30 (difference +0.0059, paired SE 0.0013). Conditional on the selected
half-life, the null rejection rate is 0.051 when 0.25 week is chosen
and 0.058 to 0.059 when 0.5, 1 or 2 weeks is chosen. Across all 12
non-ladder null cells of the gap-fill summary G9 exceeds G8 in 10 (mean
+0.006). The effect is small in practical terms (about 0.055 against
0.05), but 'nothing detectable' (Section 3.3, Conclusion 3: CR2
'absorbs the post-selection risk') is false as stated. The same
selection effect appears as optimistic standard errors under the
alternative (M3). Remedy: report the pooled and paired comparison,
state the size as approximately 0.055 under selection, and either
accept it explicitly or add a selection-aware correction (a parametric
bootstrap as in F1A0B, or a null-model selection step per M1).

**M3. Under alternatives the G9 standard error is optimistic because of
selection, and this is misattributed to dropout (VERIFIED from
`g9-gap-fill/summary.rds`).** Mean calibration ratio $\kappa$ at
alternative cells: G9 0.934 ($N = 35$, 15 cells) and 0.963 ($N = 70$,
22 cells), against G8 1.02 and 1.01, G6 1.00 and 1.02, and F1 0.99 and
1.01 (minimum for G9 0.875). The CR2 estimator without selection is
calibrated; the shortfall is specific to the data-chosen half-life.
Section 3.4 (Dropout) states that 'the part of that shortfall that
dropout causes is not separated from the part that is there already'
and the Abstract says that G9's standard error 'becomes optimistic'
under heavy missingness. The G8 comparison separates them: the
optimism is a selection effect that is present without dropout (for
instance $\kappa = 0.87$ at $N = 35$ with $c_{bm} = 0.45$ and no
dropout in the 1,000-replicate Stage 5 reference cell). Consequences:
(a) nominal CIs for G9 will under-cover; the manuscript never reports
coverage of the CR2 interval for any CR2 specification (the Stage 1
coverage column is model-based and has an estimand mismatch for
G1/G2); (b) the post-selection claim in Conclusion 3 is supported only
at the null. Remedy: report coverage and bias of G8 and G9 at
alternatives on the calibrated scale, and either qualify the claim or
adopt a selection-aware interval.

**M4. Cell-level evidence is much weaker than the 'k of n cells inside
the band' language suggests (VERIFIED arithmetic and INSPECTED seeds).**
(a) In `32-run-g9-gap-fill.R` the seed is `set.seed(26102700L + r)`
for replicate $r$ in every cell, and the other drivers use the same
pattern (`24102600L + r`, `25102600L + r`, `23102600L + r`). Every cell
at a given $N$ and design therefore reuses the same random-number
stream. This is useful for paired contrasts across cells, but it means
that '14 of 14 null cells inside the band', '16 of 16', and the seven
ladder cells are not independent confirmations; Monte Carlo noise is
largely common to them. (The prog-TV-grow rows are bit-identical to the
Hybrid null, see M7.) (b) A 500-replicate band of $[0.031, 0.069]$ has
little power to detect inflation. If the true size is 0.065 a cell lands
inside the band with probability 0.64, at 0.070 with probability 0.47,
and at 0.075 with probability 0.31 (normal approximation). 'Inside the
band' therefore cannot support 'holds nominal size'. (c) The same
Hybrid null and OL+BDC rows come from 500 replicates, and the Stage 1
duplication of the exponential cells as Weibull $k = 1$ shows how noisy
single cells are: the 54 'duplicate' cells differ from the exponential
cells by up to 0.048 in power (VERIFIED; the text says they duplicate
the exponential 'numerically'). Remedy: report size as an estimate with
a confidence interval, pooled over cells where the design allows; run
at least 5,000 to 10,000 replicates for the null cells that carry the
recommendation (Hybrid, $N \in \{35, 70\}$, two half-lives, with and
without a prognostic confound); use independent seeds across cells or
state the common-stream design explicitly; use per-replicate paired
standard errors (the `all.rds` files are in the repository) for every
ordering claim that rests on differences of one or two points.

**M5. Several counts and ranges in the Abstract and text do not match
the stored data (VERIFIED by execution).**

- 'Eight null cells' (Abstract, Section 3.3, caption of Table
  'stage3-null'): the 30-cell grid contains six null cells (3 designs
  by 2 sample sizes at $t_{1/2} = 1$); the table itself shows six rows.
  The Stage 3 null evidence also contains no CR2 null cell at
  $t_{1/2} \in \{0, 0.5\}$ for CO or OL+BDC.
- 'E0 is liberal at every null cell (0.06 to 0.15)' (Abstract) and
  'never at or below nominal' (Section 3.4): E0 rates of 0.060 (Hybrid
  $N = 70$, Stage 3), 0.060, 0.062 and 0.066 (Stage 4) lie inside the
  500-replicate band, and the text itself says 'at its edge under
  Hybrid'. E0 is clearly liberal under CO (0.128, 0.150) and OL+BDC
  (0.106, 0.144, 0.100) and marginally so under Hybrid (0.062 to 0.064
  at 2,000 replicates). 'Every' overstates.
- 'at or above G1 and F1 at every alternative but one' (Abstract,
  Section 3.4, Section 4.1): over the 37 non-dropout alternative cells
  of the gap-fill summary G9 is below G1 at three: hl-t0.2 $N = 35$
  ($-0.010$, the one cited), Weibull $k = 2$, $c_{bm} = 0.25$,
  $N = 35$ ($-0.008$) and Weibull $k = 2$, $c_{bm} = 0.45$, $N = 70$
  (a tie to four decimals). The two $c_{bm} = 0.25$ shape cells are
  not shown in the printed tables.
- 'G9 lost 17 to 26 points' at strength 0.6 (Section 3.5): against the
  stated no-confound references (0.759, 0.969) the losses are 16 and
  12 points for prog-PB and 26 and 20 for prog-TV-const, i.e. 12 to 26.
- 'F1A0B's bootstrap calibration costs a few points of power relative to
  F1 at every alternative' (Section 3.4): at the half-life
  misspecification cell F1A0B (0.228) exceeds F1 (0.192), and at
  hybrid_alt the difference is 0.002; in the 30-cell grid F1A0B exceeds
  F1 in several Hybrid cells (for example 0.513 against 0.442).
- 'power at strength 0.3 is at or above the no-confound reference on
  every pathway': prog-PB and prog-TV-const at $N = 70$ are 0.968
  against 0.969 (trivial, but the statement is literal).
- Background and Introduction: 'at a half-life of one week the relative
  loss is roughly 31% under OL+BDC but only 3 to 5% under crossover and
  Hybrid'. In this manuscript's own corrected Stage 1 grid, the G1
  relative loss at $t_{1/2} = 1$ is 47% to 66% (OL+BDC), 19% to 56%
  (Hybrid) and 7% to 16% (CO); the quoted figures resemble the
  $t_{1/2} = 0.5$ column (OL+BDC 16% to 34%; Hybrid $N = 70$,
  $c_{bm} = 0.45$ 4.8%). Either the companion paper used a different
  metric or the half-life label is wrong; as written the Background
  understates the problem by a factor of two relative to Stage 1.
- 'Across all 216 cells' (Stage 1): 54 of the 216 are the Weibull
  $k = 1$ duplicates; there are 162 distinct configurations. 'G3 is
  never worse than G1': G3 is below G1 in 17 of 162 configurations
  (minimum $-0.016$ at Hybrid), which is within noise but not 'never'.
- Section 4.3 and Conclusion 5: 'holds nominal size ... at every
  feasible strength examined' is qualified by the $N = 280$ null cells,
  where G8, G9 and F1 reject at 0.062 to 0.067 (band $[0.036, 0.064]$;
  F1 is above it at both strengths).

Remedy: regenerate every in-text number from the data objects (most
already are, through inline R) and replace the hard-coded ones; add an
automated check that fails the build when a quoted count differs.

**M6. Comparisons mix size and power, and two comparators are not
like-for-like (INSPECTED and VERIFIED).**

(a) G9 (CR2, size about 0.05) is compared with G1 and G3 (model-based,
size 0.03 to 0.04) throughout Stages 3 to 5. A conservative test has
lower power by construction, so a part of 'exposure weighting recovers
power' is size restoration. Using the stored gap-fill alternatives,
G9 minus G1 averages +3.9 points, whereas G9 minus G6 (G1 with CR2, the
equal-size comparator) averages +1.8 points, and G9 trails G6 in four of
37 cells. In the 30-cell grid CR2 raises power by 2.1 to 2.4 points on
average under Hybrid and OL+BDC, but lowers it under CO (mean $-1.4$,
minimum $-5.6$ points for G4 to G9), so 'CR2 recovers two to five
points of power' (Conclusion 3) holds for two of three designs.
Remedy: make G6 and G8 the primary comparators for G9 in every power
table, and add size-adjusted power (critical values calibrated to the
cell's null) for E0, which is liberal and for which the manuscript
rightly declines to compare power.

(b) E0 differs from G1 in more than the residual correlation. In
`fit_spec_e0` it adds the expectancy main effect `De` whenever
expectancy varies, uses `lme4`/`lmerTest` Satterthwaite rather than
`nlme` with `clubSandwich`, and drops patients never on drug. The
statement that 'dropping the residual correlation structure ... is
itself a source of Type I error inflation' (Section 3.4) is an
untested attribution. Remedy: add the two single-factor contrasts (E0
plus `corCAR1`; G1 without `corCAR1`), or soften the sentence. The
manuscript's own Section 2.4.4 describes E0 without the `De` term.

(c) F1 fixes the analysis half-life at 0.2 week while the true values
are 0.5, 1 and 2, so the 'G9 beats F1' contrasts are largely a
half-life misspecification contrast by construction; F1A0B uses a
different candidate grid ({0, 0.2, 0.5, 1, 2}) and, per the protocol,
a one-resample-per-replicate pooled ('warp-speed') bootstrap whose
critical value is estimated from the 500 resamples of the same cell.
None of this is disclosed in Section 2.4.4, which says only 'a
parametric bootstrap that repeats the selection on every resample'.
Remedy: describe the calibration exactly, and include an F-family
member with the same grid and a selected half-life without the
bootstrap (F1A0) in the headline tables.

**M7. The prognostic stress test contains a vacuous pathway and omits
the most plausible one (INSPECTED code; SUSPECTED consequence).**

(a) prog-TV-grow is implemented in `implementations/tidyverse/R/
functions.R` as `b_bm_tv * bm_z_tv * tv_sd` added at every timepoint,
that is, a constant shift proportional to the biomarker. Despite the
name it does not grow with time; it is absorbed exactly by the `bm`
main effect, which is why its rows are bit-identical to the null row
(for example 0.0461 at both strengths, and identical to
`hybrid_null` across all methods in Stage 4). It therefore tests
nothing, and the claim of 'four prognostic pathways' should read three.

(b) The mechanism most likely to confound the interaction in an
actual trial, a biomarker that predicts the rate of natural-history
change (a $bm \times t$ or $bm \times$ trend component), is not
simulated. The analysis model has no `bm:t` term, and in the Hybrid
design $D_{bc}$ is strongly negatively related to time. A back-of-
envelope projection (my own calculation, not a simulation: regress a
normalized Gompertz natural-history shape on $D_{bc}$ over the 32
path-occasions at $t_{1/2} = 0.5$) gives a slope of about $-0.19$ per
unit of shift, so a biomarker-proportional trend shift of
$0.3 \times 6$ points (the stated $b_{tv}$ and TV SD) would bias
`bm:Dbc` by roughly $-0.33$ (about 9% of the $-3.6$ target and a
non-trivial fraction of the standard error at $N = 70$). This is
suspected, not demonstrated; it should be simulated before 'no
confound effect inside the feasible region' is asserted.

(c) The 'feasible region' is a property of the correlation-matrix
parameterization: couplings enter as correlations of one latent
biomarker with the components, so PD limits them (0.306 for $c_{pb}$
at $c_{bm} = 0.45$). Clinically strong prognostic effects expressed as
mean shifts have no such bound. 'Nothing moves inside the feasible
region' is thus a statement about this construction, not about
confounding generally.

(d) The strength-0.6 null rows on prog-PB and prog-TV-const (and the
abstract's 'holds nominal size ... at both confound strengths') are
generated under a repaired matrix, that is, a different DGP with
attenuated coupling. They should not enter a robustness claim.

(e) The Discussion reports a population contrast of $-0.14$ at
$c_{pb} = 0.3$, but the stored G9 mean estimates at $N = 280$ are
$-0.013$ ($c_{pb} = 0.3$) and $-0.039$ (0.45), and at $N = 70$ about
$-0.13$ (0.45). These are on different scales (the estimand is a
regression coefficient on $D_{bc}$, not the raw on/off contrast) and
may be reconcilable, but the manuscript does not reconcile them, and
the explanation offered for the $N = 280$ excess ('a property of the
small-sample correction') is conjecture (SUSPECTED): rejection rates
rise with $N$ (0.045 at 70, 0.055 to 0.059 at 140, 0.062 to 0.063 at
280) for all CR2 tests, and a small estimand bias of 0.2 to 0.3 standard
errors would produce about 0.06.

**M8. Generality and circularity (INSPECTED).** The analysis model is
essentially the conditional mean of the Gaussian DGP. Under covariance
moderation the regression of the response on the biomarker at an
occasion is linear in $c_{bm}\phi(t_{sd})$, so $D_{bc}$ at the true
$\phi$ is the correctly specified regressor; Stage 1 is therefore
partly a demonstration that the correct model beats an incorrect one.
In the Stage 1 Weibull cells the analysis uses the true Weibull shape
(`prepare_long_data` defaults `analysis_form` to the DGP form, and
`01-run-factorial.R` passes no override), which contradicts Section 2.3
('the form every analysis specification assumes' is exponential) and the
Stage 1 sentence that only the half-life is set to truth. Shape
misspecification enters only through the gap-fill $k \in \{0.5, 2\}$
cells at $c_{bm} \in \{0.25, 0.45\}$. Everything uses one parameter set
(prazosin/PTSD trajectories, $\rho = 0.7$ except the $\rho$ grid,
$c_1 = 0.2$, a single response SD vector), one design for the
recommendation, $N \le 70$ for power, and factor-at-a-time variation
(acknowledged). The words 'robust' and 'holds' in the Abstract and
Conclusions are stronger than this. In addition, the covariance-
moderation premise is an unvalidated modeling assumption of the same
project; the recommendation is conditional on it, and the Discussion
should say so in the first paragraph, not only in the last sentence.
Remedy: vary the response SDs and trajectory parameters on a small
space-filling design; add one misspecified-mean-structure DGP (for
example a latent-class or mean-moderation-like pattern is excluded by
scope, but a different decay form together with a different $\rho$
and dropout is not); report AIC selection frequencies.

**M9. The selection step is noisier than the narrative implies, and the
selected half-life is never reported (VERIFIED from per-replicate
data).** At a true half-life of 0.5 week G9 selects 0.25 in 61% to 63%
and 0.5 in 28% to 32% of replicates; at 1 week it selects 0.25, 0.5, 1,
and 2 in about 12% to 20%, 25%, 28% to 39% and 22% to 27%. The power of
G4/G9 matches the oracle despite this, which is informative (the
exposure weight is insensitive to the half-life over this range), but
'AIC selection recovers the half-life' would be wrong and the reader
should see the distribution. If the weight is insensitive, a fixed
compromise value (a CR2 test at one pre-set half-life, as F1 does at
0.2) removes the post-selection issue altogether; G9 should be
compared with G8 at a fixed, deliberately wrong half-life (for example
0.25 and 1 against truth 0.5 and 2) to show what selection buys.

**M10. Heavy-missingness $\kappa$ is not interpretable (VERIFIED).** For
per-occasion MCAR 0.3 at $N = 35$, $\kappa = 0.49$ is produced by a
heavy-tailed estimate distribution (median $|\hat\beta|$ 3.4, maximum
130.7; SD 12.4 against a MAD-based SD of 5.1). The Type I error in that
cell is 0.022, below nominal, which is inconsistent with the
statement that 'G9's CR2 standard error understates its variability at
the smaller sample'; a robust ratio (median SE over MAD-based SD) is
0.86. The patient-level alternative $\kappa$ of 0.875 does have a
robust counterpart of 0.85 and is genuine (see M3). Remedy: define
$\kappa$ with the root mean squared SE in the numerator and a robust
scale in the denominator, or report coverage; drop the claim about
understated variability or support it with coverage.

**M11. Dropout evidence is narrower than the conclusion.** The MAR
mechanism depends on baseline severity (BL), which is not among the
columns passed to the analysis (`prepare_long_data` selects `ptID`,
`path`, `bm` and the outcome columns), so the mechanism is MAR given
the analyzed data only insofar as early outcomes proxy BL (INFERRED).
Dropout is monotone, begins after the first occasion, and is independent
of treatment and of the outcome trajectory; informative dropout
(adverse events on drug, non-response) is the clinically relevant
case for an N-of-1 design and is flagged only as not simulated. The
Abstract and Conclusion 7 ('keeps its size and its power advantage')
should say 'under MCAR and baseline-severity-dependent dropout'. With
$N = 35$ the dropout cells give G9 size 0.050 to 0.062 and 0.026 to
0.054 at $N = 70$ with an unexplained range (0.026 at MCAR 0.2,
$N = 70$) from 500 replicates; this again illustrates M4.

**M12. Random-effects check does not close the question.** The
augmentation was run on four no-confound and eight strength-0.6 cells,
eight of which are under the repaired matrix, in a DGP that contains no
patient-specific slopes. 'Adds nothing' and 'closes the random-effects
question for G9 on this design' follow from the DGP, not from evidence
about trials. At least a DGP with a patient-specific slope variance
should be included, and the 2.5% non-convergence of the three-component
fit should be shown to be uninformative (it is in the first-stage
count only). The text also states paired per-replicate differences of
at most 0.008, which I could not check from the summary (the summary
supports 'within one point').

## 5. Minor comments

**m1.** Section 2.4 introduces $\hat{t}_{1/2}$ for both the assumed and
the selected half-life; NOTATION.md (Part on estimands) reserves
$t_{1/2}^a$ for the analysis half-life. Use $t_{1/2}^a$ for the
assumed value and $\hat{t}_{1/2}^a$ for the selected one.

**m2.** $\kappa$ is the mean SE divided by the empirical SD; by
Jensen's inequality the mean of SEs is below the root mean squared SE,
so $\kappa$ is biased downward by a few percent at small $N$. Use the
root mean squared SE (or report the variance ratio).

**m3.** 'CR2 slightly overcorrects ($\kappa$ 0.91 to 0.93)' at $N = 35$:
$\kappa < 1$ means the SE is too small, that is, the CR2 estimator
under-corrects relative to the empirical variability, while Type I error
sits near 0.06. Reword, and note that this is G9 only (G8 is 0.98).

**m4.** Stage 1: the 'about four MCSEs' comparison (8 points against
0.0177) uses the unpaired MCSE of a single power estimate; the
contrast is paired, so report the paired SE from `all.rds`.

**m5.** The Stage 1 text says 'mean bias, MSE and coverage agree to
three decimal places'. Cell-level differences between G1 and G2 reach
0.018 (bias), 0.081 (MSE) and 0.026 (coverage); only the means agree.
Write 'on average'.

**m6.** Section 2.6 describes the factorial as 540 cells, then 216 and
162; the stored grid also contains a 'linear' decay form and an `E4`
specification that the text does not mention. State which are excluded
and why.

**m7.** The stored `g-vs-f-comparison` and `g9-gap-fill` summaries still
contain `meanmod_*` and `B3-meanmod` cells (and a stale `mean_moderation`
slice in the 02 grid). They are correctly unused, but a reader of the
deposited data will find them; label them as out of scope in a data
dictionary.

**m8.** Stage 5 text says the ladder reaches 'within 0.001 to 0.05 of
the ceiling on each pathway'; for prog-BL at $c_{bm} = 0.3$ the tested
0.8 is 0.015 below 0.815, for prog-PB at 0.45 it is 0.001, and for
prog-BL at $c_{bm} = 0.45$ the tested 0.45 is 0.046 below 0.496. This is
accurate but not 'each pathway' at the strongest tested level in every
$c_{bm}$ setting, and prog-PB and prog-TV-const are not run at
$c_{bm} = 0.45$ below the 0.306 ceiling beyond 0.3. State the grid.

**m9.** The ceiling for $c_{bm}$ alone (0.518 at $\rho = 0.7$, $t_{1/2} =
0.5$ week) differs from the 0.474 (Hybrid 0.480) recorded in the project
CLAUDE.md for the refined process. The difference may come from the
carryover setting, but it should be reconciled in the text, because
readers will compare the two.

**m10.** The Stage 1 and 2 tables carry no CR2 null results at
$t_{1/2} \in \{0, 0.5\}$ outside Hybrid; say so, or add them.

**m11.** The Methods say candidate fits are compared by ML and the
winner refitted by REML. Report how often the REML refit fails or
differs in convergence from the ML fit (the non-convergence column
combines both).

**m12.** Stage 4 and 5 compare cells across jobs and hosts (the
no-confound reference is job 29; the ladder is job 40 on a second host,
which the project CLAUDE.md notes is not draw-for-draw identical). These
comparisons are unpaired and carry the full Monte Carlo error; the
text treats 'power does not fall' as an exact check.

**m13.** The Weibull decay formula in Section 2.3 is correct, but the
text says all forms are 'matched to a common half-life' while the
Stage 1 'linear' form in the data is not described. Remove it from the
data or document it.

**m14.** Trials formatting (INSPECTED): the Abstract is about 840 words
including headings (Trials limits abstracts to 350 words); there is no
Keywords line; Declarations lack 'Authors' contributions',
'Acknowledgements' and 'Authors' information'; 'Availability of data
and materials' points to a project repository without a persistent
identifier (deposit a tagged release with a DOI); the PDF carries line
numbers and a table of contents that the journal will remove. An
ADEMP summary table (aims, data-generating mechanisms, estimands,
methods, performance measures) would shorten Sections 2.1 to 2.7 and
satisfy the Morris et al. reporting expectation more directly.

**m15.** The statement that the study is exploratory is correct and
should be repeated in the Abstract, because the headline sentence
('we recommend G9') reads as confirmatory.

**m16.** Reference-key and attribution problems are listed in
Section 6. In addition the in-text names in Section 4.5 ('Senn,
Julious and Araujo's finding') do not match the cited entry.

## 6. Missing and questionable references

| Location | Issue | Suggested source |
|---|---|---|
| Sec. 4.5 and `senn2014four` | In-text attribution to 'Senn, Julious and Araujo' but the bib entry is Chen and Chen (PLoS ONE); the entry's note says 'Chen, Chen and Xia'. Two authors per the journal record. The key and text must be corrected, and the Senn work cited separately if intended | Chen X, Chen P. A comparison of four methods for the analysis of N-of-1 trials. PLoS ONE 2014;9(2):e87752. doi:10.1371/journal.pone.0087752. For the Senn claim: Araujo A, Julious S, Senn S. Understanding variation in sets of N-of-1 trials. PLoS ONE 2016;11(12):e0167167 |
| Sec. 2.4.5, Abstract (CR2, Satterthwaite) | Only Bell and McCaffrey cited; the implemented estimator (`clubSandwich`) is described by Pustejovsky and Tipton | Pustejovsky JE, Tipton E. Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. J Bus Econ Stat 2018;36(4):672-683 |
| Sec. 2.4.5 (small-sample inference in mixed models) | No mention of the standard alternative (Kenward-Roger), for which a diagnostic script exists | Kenward MG, Roger JH. Small sample inference for fixed effects from restricted maximum likelihood. Biometrics 1997;53:983-997 |
| Sec. 2.4.4 and Sec. 3.3 (post-selection) | Berk et al. cited, but no evidence about uniform validity of post-selection tests | Leeb H, Potscher BM. Model selection and inference: facts and fiction. Econometric Theory 2005;21:21-59; Kuchibhotla AK, Kolassa JE, Kuffner TA. Post-selection inference. Annu Rev Stat Appl 2022;9:505-527 |
| Introduction, novelty | Prior work on exponential wash-in and wash-out carryover models in aggregated N-of-1 trials is not cited | Gartner T, Schneider J, Arnrich B, Konigorski S. Comparison of Bayesian networks, G-estimation and linear models to estimate causal treatment effects in aggregated N-of-1 trials with carry-over effects. BMC Med Res Methodol 2023;23:191 |
| Introduction, simulation of N-of-1 designs with carryover | Not cited | Percha B, Baskerville EB, Johnson M, Dudley JT, Zimmerman N. Designing robust N-of-1 studies for precision medicine: simulation study and design recommendations. J Med Internet Res 2019;21(4):e12641 |
| Sec. 3.1 and 3.4 (Type I error when carryover is ignored) | Not cited; consistent with the E0 and G1 findings | Blackston JW, Chapple AG, McGree JM, McDonald S, Nikles J. Comparison of aggregated N-of-1 trials with parallel and crossover randomized controlled trials using simulation studies. Healthcare 2019;7(4):137 |
| Introduction (carryover in crossover designs; the lagged-term device) | The case against testing or adjusting for carryover with a single lag term is classical | Senn S. Cross-over Trials in Clinical Research. 2nd ed. Wiley; 2002; Freeman PR. The performance of the two-stage analysis of two-treatment, two-period crossover trials. Stat Med 1989;8:1421-1432 |
| Introduction (pooling N-of-1 trials) | Foundational aggregation methodology absent | Zucker DR, Ruthazer R, Schmid CH. Individual (N-of-1) trials can be combined to give population comparative treatment effect estimates: methodologic considerations. J Clin Epidemiol 2010;63:1312-1323 |
| Sec. 2.7 (simulation reporting) | Morris et al. cited; fuller reporting templates and cautions on common random numbers are not | Siepe BS, Bartos F, Morris TP, Boulesteix AL, Heck DW, Pawel S. Simulation studies for methodological research in psychology: a standardized template for planning, preregistration, and reporting. Psychol Methods 2024; Pawel S, Kook L, Reeve K. Pitfalls and potentials in simulation studies. Biom J 2024 |
| Introduction, estimand | Power for treatment-by-covariate interactions and the sample-size cost relative to the main effect | Brookes ST, Whitely E, Egger M, Smith GD, Mulheran PA, Peters TJ. Subgroup analyses in randomized trials: risks of subgroup-specific analyses; power and sample size for the interaction test. J Clin Epidemiol 2004;57:229-236 |
| Reporting, Trials | CENT is the reporting guideline for N-of-1 trials; the Declarations do not state its relevance (real-trial reporting of the recommended analysis) | Vohra S, Shamseer L, Sampson M, et al. CONSORT extension for reporting N-of-1 trials (CENT) 2015. BMJ 2015;350:h1738 |
| Sec. 2.4.1, mixed-model AIC | AIC from ML fits of LMMs with differing fixed effects: no discussion of the conditional-AIC literature | Vaida F, Blanchard S. Conditional Akaike information for mixed-effects models. Biometrika 2005;92:351-370 |
| Companion references | Both companion papers are 'unpublished'; Sections 1 and 2.4.5 rest on figures (31%, 3 to 5%, calibration findings) that a reader cannot check | Provide preprint DOIs or include the figures in a supplement |

Not cited in the manuscript and unverified by me as to relevance: the
Bell-McCaffrey result for models with correlation structures. I did not
verify that `clubSandwich` CR2 for `lme` fits with `corCAR1` and a
random intercept reproduces the stated theoretical properties (its
adjustment matrices use the fitted working covariance); the stored
size results at 2,000 replicates are the only evidence offered.

## 7. Suggestions for strengthening the paper

- Reanalyze existing `all.rds` files: pooled and paired size of G9 and
  G8, selection frequencies, coverage on the calibrated scale, robust
  $\kappa$, paired standard errors for every ordering claim.
- Add a null-model selection variant of G9 and a fixed-half-life CR2
  variant (two or three wrong values) so that the value of selection is
  quantified against its size cost.
- Simulate a genuine $bm \times$ trend prognostic pathway as a mean
  shift (no PD limit), and a patient-level slope variance, and delete or
  rename prog-TV-grow.
- Use independent seeds per cell, run the null cells that carry the
  recommendation at 5,000 to 10,000 replicates, and tabulate size with
  confidence intervals.
- Replace G1/G3 as power comparators by G6/G8, report size-adjusted
  power for E0, and report separate single-factor contrasts for the E0
  differences (`corCAR1`, `De`, software).
- Present a short decision table: which analysis for which design, with
  the evidence tier (verified in cells, extrapolated) for each cell.
- Tighten the abstract to the Trials limit, replace 'robust' by the
  specific conditions examined, and add a limitations sentence on single
  parameter set and covariance-moderation premise.
- Deposit the repository snapshot with a DOI.

## 8. Scope of this review

Verified (executed read-only R on stored data): Stage 1 grid
(`02-grid-summary.rds`): G1/G2 differences, G3 versus G1 by design,
null rates, relative-loss figures, exponential versus Weibull $k = 1$
duplicates; `tier1-g-vs-f-e0` (cell counts, null table, CR2 power
changes, G9 leadership); `g-vs-f-comparison` (E0 range, F1 versus F1A0B,
G9 versus F family); `g9-gap-fill` summary (size diagnostic values and
$\kappa$, scenario grids, counts of cells inside the band, G9 versus
G1/G6/G8/F1, dropout ranges, $\kappa$ by method and $N$); per-replicate
`all.rds` for jobs 31, 32, 33, 40 (pooled and paired size, selection
frequencies, heavy-tail $\kappa$); `g9-prognostic-stress` and
`g9-random-effects` summaries (16 null cells, losses, nonconvergence,
reference $\kappa$); `pd-check.rds` and `coupling-ceilings.rds`
(contrasts, repaired paths, eigenvalue signs at the ceilings). Monte
Carlo band arithmetic ([0.031, 0.069], [0.036, 0.064], [0.040, 0.060])
and the probability of landing inside the band under inflated size were
recomputed.

Inspected (source read, not executed): `simulation-core.R` (design
presets, `exposure_weight`, `select_half_life_aic`, `cr2_extract`,
`fit_spec_s7`, `fit_spec_e0`, F-family code, `apply_dropout`,
`simulate_cell*`), `32-run-g9-gap-fill.R`, the seed lines of drivers
27, 28, 29 and 06, the prognostic-pathway code in
`implementations/tidyverse/R/functions.R`, the gap-fill protocol,
`references.bib`, NOTATION.md (specification labels, $\kappa$,
$t_{1/2}^a$), the project CLAUDE.md. `33-assemble-gap-fill.R` was not
read line by line.

Not checked: whether the figures (`02xs-*.pdf`) match the data; the
rendered PDF; `clubSandwich` behavior for `lme` with `corCAR1`
(theory); paired per-replicate differences quoted for the random-effects
fits; the Stage 1 factorial per-replicate data; the companion papers;
that 'the published analysis' reproduces commit 58b32a9 (taken from the
code comment); the full bootstrap implementation in
`06-cr2-recalibration.R` and the supplement. Full re-execution of any
simulation was not attempted; no repository file other than this report
was modified.

Literature searches run (web search, no personal identifiers sent to any
service): the authorship and content of Chen and Chen 2014; Pustejovsky
and Tipton 2018 (CR2, Satterthwaite, `clubSandwich`); Percha et al.
2019; Blackston et al. 2019; Gartner et al. 2023 (page fetched and
summarized); and a search for 2023 to 2025 work on biomarker-by-
treatment interaction analyses with carryover in aggregated N-of-1
trials, which returned nothing newer than Hendrickson 2020 on that exact
topic. References listed from memory without a search in this session
(Kenward-Roger 1997, Leeb-Potscher 2005, Kuchibhotla et al. 2022,
Zucker et al. 2010, Araujo et al. 2016, Freeman 1989, Senn 2002,
Vaida-Blanchard 2005, Brookes et al. 2004, Vohra et al. 2015, Siepe et
al. 2024, Pawel et al. 2024) should have their volume and page details
checked before they are added to the bibliography. Confidence: high for
the items labeled VERIFIED; moderate for the INSPECTED items (a reading
of code, not a rerun); low for the SUSPECTED consequences in M7.
