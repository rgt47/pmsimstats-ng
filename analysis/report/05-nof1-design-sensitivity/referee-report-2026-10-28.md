# Referee report: paper 05, hybrid N-of-1 design sensitivity
*2026-10-28 07:36 PDT*

**Manuscript.** `analysis/report/05-nof1-design-sensitivity/report.Rmd`,
"Power sensitivity of a hybrid open-label-plus-blinded-discontinuation
N-of-1 design under carryover, serial correlation, and
variance-component variation" (pmsimstats team; cover letter addressed
to *Statistics in Medicine*).

**Evidence labels.** VERIFIED: re-derived or executed by the referee.
INSPECTED: read in source or driver code and judged. SUSPECTED: a
plausible reading not confirmed by execution. UNVERIFIED: not checked.

## 1 Summary of the manuscript

The manuscript reports twelve one-axis sensitivity sweeps (S1 to S12)
of Monte Carlo power for the treatment main effect, estimated as the
`Dbc` coefficient of an `nlme::lme` fit with a random intercept and a
`corCAR1` residual, under a multivariate-normal data-generating process
with modified-Gompertz component means. The primary design is described
as the Hendrickson et al. (2020) Hybrid (open-label titration, blinded
discontinuation, brief crossover). Comparators are an alternating
A/B/A/B aggregated N-of-1 design, a two-period crossover and a
parallel-group RCT, all at a stated baseline of $N = 35$, eight weekly
visits, a drug-response ceiling of 2 units, a 3-day carryover
half-life, 500 replicates per cell (2,000 for S12). In the authors'
words, the contributions are: (i) a characterization of the hybrid's
power profile along effect size, carryover half-life, decay form,
between- and within-patient variance, cycles, participants, AR(1)
correlation, period length, biomarker strength by DGP architecture, a
misspecification factorial and a null calibration; (ii) Morris-style
MCSE for every performance measure; (iii) identification of the regimes
in which the hybrid is preferred to the RCT and where that ordering
reverses; (iv) the finding that the two-period crossover analyzed with
`corCAR1` has inflated Type I error.

## 2 Overall assessment and recommendation

**Recommendation: reject in the present form.** A resubmission could be
considered only after the simulation is rerun under a corrected driver
and the manuscript is rewritten from the regenerated artifacts. Major
revision would understate what is required, because the defects reach
the design being evaluated, the analysis model being described, the
comparator against which the design is judged, and the correspondence
between the prose and the shipped data.

Justification, in order of consequence:

1. The numbers in Results for S2 to S8 and S10 to S12 do not reproduce
   from the artifacts the manuscript says they were read from, and the
   text contains no code that reads any artifact (VERIFIED, Section 4,
   comment 1).
2. The analysis never models carryover. Every sweep fits a binary drug
   indicator (`carryover_t1half = 0` in the analysis options), so the
   statements in S3 ("the analysis model always assumes exponential
   decay") and S11 ("the decaying `Dbc` in the analysis absorbs
   carryover correctly") describe a model that was not fitted
   (VERIFIED, comment 2).
3. Sweep S9 never varies the carryover half-life: the driver passes only
   the first row of its parameter grid, so all fifteen cells were run at
   a 0.5-day half-life and relabeled (VERIFIED, comment 3).
4. The design called "the Hendrickson 2020 hybrid" is a two-path,
   weekly, eight-visit schedule with two open-label visits and no
   placebo-first crossover sequence. It is not the published four-path
   Hybrid at weeks 4, 8, 9, 10, 11, 12, 16, 20, and the manuscript's
   claims of reproducing Hendrickson "by construction" and of a
   "four-week open-label titration phase" are false (INSPECTED, comment
   4).
5. The parallel-RCT comparator is biased under the null in every null
   cell (mean estimate about $-0.12$ when the true effect is 0, 14 to 25
   standard errors from zero), and its estimate drifts with the AR(1)
   parameter. The manuscript's explanations of S8 and S11 rest on this
   artifact without recognizing it (VERIFIED bias; SUSPECTED mechanism;
   comment 5).

Several further findings would each on their own require major
revision: the four designs estimate different quantities (comment 6);
non-positive-definite covariance matrices were silently repaired in
S8 and S10 and the current package would refuse them (comment 7); the
within-subject designs show a consistent mild Type I inflation that
the text describes as nominal (comment 8); the "between-patient
heterogeneity" sweep does not vary treatment heterogeneity, and the DGP
contains none (comment 9); a "prototype confirmation" section describes
a different study whose results contradict S7 (comment 10).

The study question is worthwhile and the sweep structure is sensible.
A corrected version with a faithful Hybrid, a sound RCT analysis, a
defined estimand, carryover-aware and carryover-naive analyses crossed
with the DGP, and prose generated from the data would be a useful
contribution.

## 3 Significance and novelty

**What is claimed.** A systematic multi-axis sensitivity study of a
hybrid aggregated N-of-1 design for the main effect, with MCSE, where
prior work varied one or two axes.

**Competing literature.** The axes are each covered by prior simulation
or analytic work, and the manuscript cites most of it, but the paper's
positioning against that work is inaccurate in several places
(Section 6). The relevant sources, checked on PubMed on 2026-10-28:

- Blackston JW, Chapple AG, McGree JM, McDonald S, Nikles J. Comparison
  of aggregated N-of-1 trials with parallel and crossover randomized
  controlled trials using simulation studies. *Healthcare* 2019;
  7(4):137. doi:10.3390/healthcare7040137. Three cycles, four variance
  structures, constant additive carryover, selection bias; parallel RCT
  uses one observation per patient.
- Wang Y, Schork NJ. Power and design issues in crossover-based N-of-1
  clinical trials with fixed data collection periods. *Healthcare*
  2019; 7(3):84. doi:10.3390/healthcare7030084. Single-subject
  analytic power under AR(1), period stacking, heteroscedasticity,
  washout; no parallel-group arm.
- Percha B, Baskerville EB, Johnson M, Dudley JT, Zimmerman N.
  Designing robust N-of-1 studies for precision medicine: simulation
  study and design recommendations. *J Med Internet Res* 2019;
  21(4):e12641. doi:10.2196/12641. Block count, ordering, treatment
  duration and sampling frequency with carryover and wash-in. Not
  cited; directly overlaps S6 and S9.
- Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
  aggregated N-of-1 trial designs for predictive biomarker validation.
  *Front Digit Health* 2020; 2:13. doi:10.3389/fdgth.2020.00013.
- Senn S. Sample size considerations for n-of-1 trials. *Stat Methods
  Med Res* 2019; 28(2):372-383. doi:10.1177/0962280217726801. Variance
  components: within-cycle within-patient variation and
  patient-by-treatment interaction.
- Chen X, Chen P. A comparison of four methods for the analysis of
  N-of-1 trials. *PLoS ONE* 2014; 9(2):e87752.
  doi:10.1371/journal.pone.0087752. CS versus AR(1), 20% carryover,
  paired t versus mixed models.
- Araujo A, Julious S, Senn S. Understanding variation in sets of N-of-1
  trials. *PLoS ONE* 2016; 11(12):e0167167.
  doi:10.1371/journal.pone.0167167. The critique of Chen and Chen that
  the models omit the patient-by-treatment interaction.
- Gärtner T, Schneider J, Arnrich B, Konigorski S. Comparison of
  Bayesian networks, G-estimation and linear models to estimate causal
  treatment effects in aggregated N-of-1 trials with carry-over
  effects. *BMC Med Res Methodol* 2023; 23:191.
  doi:10.1186/s12874-023-02012-5. Analysis-method comparison under
  carryover; not a decay-form misspecification study.
- Carrozzo AE, Zimmermann G, Bathke AC, Neunhaeuserer D, Niebauer J,
  Kulnik ST. Two-arm crossover randomized controlled trial versus
  meta-analysis of N-of-1 studies. *Biom J* 2025; 67(2):e70045.
  doi:10.1002/bimj.70045. The comparator is a two-arm *crossover* RCT.
- Lin SX, Morrison L, Smith PWF, Hargood C, Weal M, Yardley L.
  Properties of bootstrap tests for N-of-1 studies. *Br J Math Stat
  Psychol* 2016; 69(3):276-290. doi:10.1111/bmsp.12071. AB versus ABAB
  under autocorrelation. Not cited.
- Tang J, Landes RD. Some t-tests for N-of-1 trials with serial
  correlation. *PLoS ONE* 2020; 15(2):e0228077. In the bibliography,
  not cited in the text.
- Shi D, Ye T. Behavioral carry-over effect and power consideration in
  crossover trials. *Biometrics* 2024; 80(2):ujae023.
  doi:10.1093/biomtc/ujae023. Crossover-versus-parallel power trade-off
  under carryover. In the bibliography, not cited.
- Jones B, Kenward MG. *Design and Analysis of Cross-Over Trials*, 3rd
  ed. Chapman and Hall/CRC, 2014; Senn S. *Cross-over Trials in Clinical
  Research*, 2nd ed. Wiley, 2002. The standard sources for period
  effects, carryover and the two-period design. Neither is cited.
- The project's own `docs/39-nof1-design-optimization.md` and
  `docs/40-hendrickson-hybrid-design-primer.md` contain closed-form and
  Monte Carlo design comparisons (switch count, spacing, open-label
  cost, block type) that overlap S2, S6 and S9 for the interaction
  estimand and reach conclusions the present manuscript does not engage
  with (for example, that under AR(1) correlation more switches help
  and that open-label visits cost power).

**What would be new if the study were sound.** A main-effect power
profile for an open-label-plus-blinded-discontinuation design under a
component-based DGP with expectancy, with carryover and AR(1) crossed,
would be new. Percha et al. (2019) and Wang and Schork (2019) cover
cycles, period length and serial correlation for alternating designs;
Blackston (2019) and Carrozzo (2025) cover the N-of-1-versus-RCT
comparison with simpler DGPs. The distinctive element is the hybrid
structure, and that is precisely the element the implementation does
not deliver (comment 4). As submitted the paper does not clear the bar
of *Statistics in Medicine* or *Trials*. To clear it, the study would
need the faithful Hybrid, a defined estimand common to the designs, a
correct RCT analysis, and the carryover-aware analysis the text already
claims.

**Contemporariness.** The bibliography reaches 2025 (Carrozzo, Jiang,
Selukar). Percha 2019, Lin 2016 and Araujo 2016 are missing; Senn 2024
(paired cycles) and Shi and Ye 2024 are in the `.bib` but unused.

## 4 Major comments

Numbering is by consequence, correctness first. Line numbers refer to
`report.Rmd`.

### 4.1 Results prose does not reproduce from the shipped artifacts

*Label.* VERIFIED.

*Location.* Results, lines 352-358 ("the numerical values quoted below
are read directly from `analysis/data/derived_data/sensitivity/`"),
Abstract lines 60-62 ("compute MCSE at render time"), TODO box lines
329-350, and every sweep subsection.

*Problem.* No chunk in `report.Rmd` reads an `.rds` file; all twelve
chunks are `knitr::include_graphics()`. Every number and every MCSE in
the prose is typed by hand. I recomputed each cell from the twelve
artifacts (script in my scratchpad; summaries below). S1 and S9 match
the prose. S2 to S8 and S10 to S12 do not:

| Sweep | Cell | Prose | Artifact |
|---|---|---|---|
| S2 | Hybrid, 7 / 10 / 14 d | 0.976 / 0.900 / 0.736 | 0.970 / 0.870 / 0.692 |
| S2 | Nof1-alt, 0 / 3 / 14 d | 0.896 / 0.778 / 0.304 | 0.882 / 0.734 / 0.280 |
| S2 | RCT range | 0.92 to 0.96 | 0.900 to 0.924 |
| S3 | Nof1-alt exp / linear / W0.8 | 0.816 / 0.904 / 0.754 | 0.718 / 0.898 / 0.718 |
| S4 | RCT, tau 1.0 / 2.0 | 0.978 / 0.874 | 0.962 / 0.866 |
| S5 | RCT, sigma 1.0 / 1.5 | 0.942 / 0.814 | 0.888 / 0.766 |
| S6 | RCT, k = 2 / 8 | 0.292 / 0.940 | 0.224 / 0.902 |
| S8 | RCT, rho 0 / 0.9 | 0.984 / 0.746 | 0.970 / 0.736 |
| S10 | Nof1-alt range | 0.682 to 0.916 | 0.616 to 0.884 |
| S11 | RCT null, AR(1) high | 0.112 to 0.126 | 0.092 to 0.116 |
| S12 | Hybrid / Nof1-alt / XO | [0.056, 0.062] / [0.057, 0.064] / [0.100, 0.130] | [0.054, 0.057] / [0.058, 0.066] / [0.108, 0.122] |

The TODO box states that the artifacts are stale at $N = 40$ and the
code is at $N = 35$. The reverse is true: the stored `N` column is 35
in S2 to S8 and S10 to S12 (VERIFIED), and 40 only in S1 and S9, whose
drivers hard-code `N = 40`. The prose is the stale object. The figures
(rendered 2026-05-08 from the artifacts) therefore disagree with the
text printed beside them.

*Why it matters.* A simulation paper whose reported numbers cannot be
regenerated from its own data is not reproducible, and the discrepancy
is not uniform in sign or size (S3 exponential: 0.816 reported versus
0.718 stored, 5 MCSE apart).

*Remedy.* Read the artifacts in chunks, compute every performance
measure and MCSE in code, and use inline R for every number in the
text. Replace `N = 40` in `01-effect-size.R` and `09-period.R` with the
common baseline. Remove the TODO box after the rerun and state the
replicate count per sweep in a table.

### 4.2 The analysis never models carryover

*Label.* VERIFIED by inspection of the driver and package source.

*Location.* Methods lines 253-263; S3 lines 466-470 and 486-489; S11
lines 737-743; Conclusions line 1024.

*Problem.* `default_analysisparams()` in `00-common.R` sets
`carryover_t1half = 0`, and `lme_analysis()`
(`implementations/nof1power/R/analysis.R`, lines 118-126) then sets
`Dbc = 0` at every off-drug visit. Every sweep therefore fits the
binary indicator (specification E1 in the project's vocabulary), not a
decaying $D_{bc,it}$. The S3 claim that "the analysis model always
assumes exponential decay in constructing `Dbc`" and the S11 claim
that "the decaying `Dbc` in the analysis absorbs carryover correctly"
are both false. The S11 driver's own header says "never-carryover
`Dbc`". The fitted formula is also not the one printed: with
`target_coef = "main"` the model is `Sx ~ bm + t + Dbc` (plus `De`
when expectancy varies), without `bm:Dbc` (analysis.R lines 171-181).

*Why it matters.* S2, S3 and S11 are presented as tests of a
carryover-aware analysis against carryover in the DGP. They are tests
of a carryover-naive analysis. The interpretation of S3 (decay-form
"misspecification" relative to an exponential assumption) has no
referent, because no decay form is assumed on the analysis side. The
Discussion's contrast with Blackston's Type I inflation "when not
modeled" is then a contrast with the same unmodeled situation.

*Remedy.* Cross the analysis half-life with the DGP half-life (at
least $t_{1/2}^a \in \{0, \text{true}\}$) in S2, S3 and S11, and report
both. Print the formula actually fitted. If the authors intend E1
throughout, say so and rewrite S3 and S11 accordingly.

### 4.3 S9 never varies the carryover half-life

*Label.* VERIFIED.

*Location.* S9, lines 655-676; driver `09-period.R` line 44
(`modelparams = model_params[1, , drop = FALSE]`).

*Problem.* The driver builds a fifteen-row grid but passes only its
first row, so all fifteen design cells ran at $t_{1/2} = 0.5$ days. The
stored `carryover_t1half` column takes the single value 0.0714 weeks
(VERIFIED). The `halflife_days` label is then attached from the design
index, so cells labeled 3 and 7 days are relabeled 0.5-day cells. The
"finding" that "at fixed period length, the carryover half-life has
only a modest effect" is the absence of any variation plus Monte Carlo
noise at 300 replicates. The sweep also uses $N = 40$ and 300
replicates, not the stated 35 and 500.

*Remedy.* Pass the full grid (fifteen model-parameter rows, or a
design-by-parameter pairing), rerun, and rewrite S9. Add an assertion
in every driver that the number of distinct stored parameter values
equals the grid size.

### 4.4 The design evaluated is not the Hendrickson Hybrid

*Label.* INSPECTED.

*Location.* Abstract lines 54-57; Introduction lines 107-110; Pipeline
lines 177-179; S2 lines 442-447 ("four-week open-label titration
phase"); Discussion lines 855-858 ("direct tidyverse re-implementation
... reproduces Hendrickson by construction"); Limitations lines 917-923.

*Problem.* `design_hybrid()` in `00-common.R` builds eight weekly
visits at weeks 1 to 8 with expectancy 1 at weeks 1 and 2 and 0.5
thereafter, and two paths: on at weeks 1-4, off 5-6, on 7-8; and on at
weeks 1-6, off 7-8. The published Hybrid has visits at weeks 4, 8, 9,
10, 11, 12, 16, 20, an eight-week open-label phase, four paths from two
randomizations, and a crossover in which one sequence is placebo first
(`docs/40`, Section 2, verified against commit 58b32a9). The
implemented design has a two-visit open-label phase, no path that is
off drug and then on drug in a blinded crossover of equal periods, and
re-exposure in only one path. The statement that the open-label phase
is four weeks is wrong either way (two visits at weeks 1-2 in the
code; eight weeks in the publication). The claim that the framework
"reproduces Hendrickson by construction" cannot hold for a different
schedule, and the companion-paper language about the "Hendrickson 2020
minimum" of $N = 35$ is unsupported (Hendrickson's figures show
$N = 35$ panels; no minimum is defined).

The comparators carry a second confound. `design_rct()` and
`design_nof1_alternating()` set expectancy 1 at every visit, while the
hybrid is blinded (0.5) from week 3. In this DGP expectancy scales both
the mean and the standard deviation of the placebo-response component
(`generate_data.R`, lines 232-234), so the RCT and the alternating
design run with twice the placebo-response standard deviation of the
hybrid's blinded visits. Part of the hybrid's apparent advantage is
built into the comparators' noise.

*Why it matters.* The paper's unit of contribution is the hybrid
design. If the schedule is not the Hybrid, the title, the abstract and
the comparison with Hendrickson are unsupported, and the design's
power is a property of a schedule the reader is not told about.

*Remedy.* Either implement the published Hybrid (four paths, the
published visit weeks, the published expectancy pattern; a
`build_hendrickson_design()` helper is already planned in
`planning/27-...md`) or rename the design and describe it visit by
visit in a table with path, week, drug state and expectancy. Put all
comparators on the same blinding assumption, or state the expectancy
difference and quantify its effect.

### 4.5 The parallel-RCT comparator is biased under the null

*Label.* VERIFIED (the bias and the SE miscalibration); SUSPECTED (the
mechanism).

*Location.* S1 line 414; S8 lines 628-645; S11 lines 729-743; S12;
Limitations lines 931-946; Synthesis lines 813-822.

*Problem.* In every null cell the RCT's mean estimate is far from zero:
S1 ($\delta = 0$) $-0.115$ (empirical SE 0.220, 1,000 replicates, 17
standard errors); S11 null cells $-0.108$ to $-0.136$; S12 $-0.070$ to
$-0.140$ at 2,000 replicates (25 standard errors). The hybrid's null
means are $-0.014$ to $-0.030$. The RCT's model-based SE is also
miscalibrated: 0.29 against an empirical 0.22 at $\rho \le 0.3$
(conservative, Type I 0.014 to 0.034), and 0.24 against 0.26 to 0.28 at
$\rho = 0.7$ (anticonservative, 0.09 to 0.12). Under the alternative
the RCT's mean estimate falls from $-1.00$ at $\rho = 0$ to $-0.59$ at
$\rho = 0.9$ (S8), which the manuscript reports without comment.

A plausible mechanism (SUSPECTED, not tested): `lme_analysis()`
prepends a baseline row at $t = 0$ with `Dbc = 0` and `Sx` equal to
the noise-free baseline, so in the RCT `Dbc` equals $\mathbf 1[t > 0]$
in the drug arm only. The mean trajectory is Gompertz in both arms but
the model has one shared linear `t`. The drug arm has a free parameter
for the baseline-to-follow-up step and the placebo arm does not, so the
misfit of the shared line is loaded onto `Dbc`. That produces a
non-zero expectation under the null and makes the estimate depend on
the GLS weights, hence on $\rho$. The manuscript's explanation of the
S11 inflation ("the variance estimator still treats each participant as
providing approximately independent observations") does not account
for the bias or for the conservatism at $\rho = 0$.

*Why it matters.* Every hybrid-versus-RCT conclusion (sample-size
ratio, robustness to variance components, AR(1) asymmetry, Type I
asymmetry) is a comparison against a comparator whose estimator is
biased and whose standard error is wrong in both directions. The
Limitations' promise that "a better-specified parallel-RCT analysis
... would likely restore nominal Type I error" is untested and may be
false if the bias is structural.

*Remedy.* Analyze the RCT with a standard arm contrast: an MMRM with
visit as a factor and a visit-by-arm interaction, or ANCOVA on the
follow-up mean adjusted for baseline; drop the baseline row from the
`Dbc` coding or code `Dbc` as the arm indicator at every row. Verify
zero mean under the null and calibrated SE before any comparison. Report
bias and SE calibration for every design in every sweep (the ADEMP
protocol promises both).

### 4.6 The designs estimate different quantities

*Label.* VERIFIED. The "effect size" axis is a Gompertz ceiling, not
the estimand.

*Location.* S1 lines 392-405; S2; S6 lines 573-585; Methods lines
219-221; Discussion throughout.

*Problem.* Time on drug restarts from zero at each re-exposure
(`trial_design.R`, lines 50-57), and the drug response follows
$G(t_{od}; m = 2, d = 1, r = 0.3)$. The realized on-drug mean is
therefore 0.34 after one week and 1.73 after eight (VERIFIED by
evaluating `mod_gompertz`). The stored mean estimates at the baseline
cell are: alternating N-of-1 $-0.31$ to $-0.37$; RCT $-0.89$; hybrid
$-0.92$ to $-1.06$; crossover $-1.26$. The four designs are being
compared on power to detect four different true effects, none of which
is the stated "$-2$ nightmares per week" (the Gompertz ceiling $m$,
never attained in eight weeks). The S1 "$\delta$" grid is a grid of
ceilings. The S6 RCT mean estimate rises from $-0.43$ at $k = 2$ to
$-1.32$ at $k = 16$ because longer trials accumulate more Gompertz
response: the "cycles" axis is confounded with trial duration and with
the effect itself.

No true value is defined anywhere, so bias, MSE and coverage (promised
in the ADEMP protocol and in Methods line 322-323) are not reported and
cannot be. The `ETbeta` column stored in every artifact (the pooled
large-sample fit) could serve as a target but is unused.

*Why it matters.* Power comparisons across designs with different
estimands are not design comparisons. Senn's and Blackston's
comparisons hold the effect fixed.

*Remedy.* Define the estimand (for example, the mean on-minus-off
difference over the design's visits, or the Gompertz value at a fixed
time on drug) and either match designs on it or report it per design
beside the power. Report bias, empirical SE, model SE, coverage and
non-convergence per cell. Consider a DGP option in which response does
not restart from zero (`docs/38`).

### 4.7 Non-positive-definite matrices were silently repaired

*Label.* VERIFIED by execution. No flag was stored, and the current
package would refuse the same cells.

*Location.* Methods lines 250-253; S3 lines 470-475; S8; S10.

*Problem.* The stored runs predate the 2026-10-26 change of the default
to `stop_on_nonpd = TRUE`. I rebuilt the S10 and S8 covariance matrices
with the current package: at $c_{bm} = 0.6$ every path of every design
except the RCT placebo arm is non-PD; at $c_{bm} = 0.3$, $\rho = 0.9$
the hybrid's paths are non-PD (minimum eigenvalue of the repaired
matrix 0.002, the `corpcor` tolerance). Those S10 and S8 cells were
therefore run on repaired matrices whose realized $c_{bm}$ is not the
nominal value, which is the attenuation documented in the project's
`docs/65`. The manuscript says a "per-cell flag" is returned; no such
flag exists in the stored results (`make_result_row` in `simulation.R`
does not carry `was_pd_corrected`). Under the current default the S10
and S8 drivers would stop with an error, so the rerun demanded by the
TODO box cannot complete as written.

A second defect: `carryover_decay(form = "power")` evaluates
`pmax(0, (1 - tsd / (3 t_{1/2}))^{1.5})`, which is `NaN` when the base
is negative. For the hybrid at $t_{1/2} = 3/7$ week and $t_{sd} = 2$
weeks the base is negative, the correlation matrix contains `NaN`, and
`build_sigma_matrix()` errors in `eigen()` (VERIFIED: my check script
halted there). The S3 text attributes this to non-PD matrices; it is a
coding error, which also affects the component-mean carryover.

*Remedy.* Compute the feasibility bound (`docs/34`, `docs/35`) for
every design-cell before running; drop or flag infeasible cells; store
the repair flag; rerun S8 and S10 under refusal. Fix the power form
(`pmax(0, 1 - tsd/max_time)^power`). State in Methods which rule was in
force for every reported cell.

### 4.8 Type I error claims for the within-subject designs

*Label.* VERIFIED (arithmetic on the stored rates).

*Location.* S1 lines 412-419; S11 lines 718-720; S12 lines 757-765;
Synthesis lines 823-828.

*Problem.* At 2,000 replicates, two MCSE at nominal is 0.0097, so the
largest rate "within 2 MCSE of nominal" is 0.060. The text's alternating
range [0.057, 0.064] and the stored 0.066 both exceed it. All eight
S12 within-subject cells lie above 0.05 (stored 0.054 to 0.066; pooled
0.058 over 16,000 replicates, MCSE 0.0018, 4.5 MCSE above nominal).
S1 gives 0.066 at 1,000 replicates and S11 gives 0.066 in one hybrid
cell. The pattern is a consistent anticonservatism of about one
percentage point, with model SE below empirical SE (S1 hybrid: 0.159
against 0.168), which the manuscript describes as "near nominal" and
"controlled." The S12 claim about the RCT (three conservative, one
nominal) is correct for the stored data but is a symptom of comment
4.5, not of good calibration.

*Remedy.* Report rates with MCSE and state the test used for
"nominal" (a binomial test or a 2-MCSE band computed at the nominal
rate). Report the ratio of mean model SE to empirical SD ($\kappa$ in
the project's notation) per cell. Discuss the small-sample
anticonservatism of `lme` Wald tests at $N = 35$ and whether a
Satterthwaite or Kenward-Roger correction or CR2 standard errors (used
elsewhere in this compendium) are needed.

### 4.9 The variance-component sweeps are mislabeled

*Label.* INSPECTED.

*Location.* S4 lines 499-522; S5 lines 532-552; S6 lines 562-585;
Introduction lines 130-134; Synthesis lines 795-802.

*Problem.* S4 varies the standard deviation of the baseline level
`BL`, a random intercept that enters every visit identically. Senn's
$\tau^2$ (and NOTATION.md's) is the between-patient variance of the
treatment effect. The DGP has no patient-by-treatment heterogeneity:
every participant shares the Gompertz ceiling and individual
differences in response enter only through the biomarker correlation
and the noise. Senn's decomposition, in which the N-of-1 design is
needed to estimate $\tau^2$ at all, is not tested here. The RCT's loss
of power with baseline SD is the familiar unadjusted-analysis loss of a
parallel design, not an N-of-1 result. S5 varies only the
drug-response component SD, not the residual or total within-patient
SD. S6 varies the number of weekly visits (and so the trial length and
the realized effect; comment 4.6), not cycles; the text alternates
between "cycles" and "timepoints."

*Remedy.* Add a participant-level random slope on the drug response
(`v_i \Psi_{it}` with $\mathrm{Var}(v_i) = \tau^2$, as in the paper 04
protocol) and sweep it; rename S4 to "baseline SD" if it is retained.
Define the within-patient variance swept in S5 in terms of the outcome.
Re-specify S6 as cycles at fixed trial length or at fixed per-cycle
response.

### 4.10 The "prototype confirmation" section is a different study

*Label.* VERIFIED. Its results contradict S7.

*Location.* Discussion lines 975-1014; Reproducibility lines 1067-1095.

*Problem.* The section reports a $c_{br} \times N$ grid at
`analysis/data/quick-sim/05-design-replicates.rds`. That file comes
from `analysis/scripts/quick-sim/05-design-prototype.R`, which uses
different response parameters (ceilings 10.99, 6.51, $-2$; SDs 5, 2,
5; baseline SD 10), a one-week half-life with analysis half-life 1.0,
and `c_br` as the AR(1) correlation of the drug-response component.
Its stored summary gives power 0.17 to 0.34 at $N = 25$ to 50 and
positive mean estimates (1.0 to 1.3), the opposite sign to every sweep.
The section is not a confirmation of anything in S1 to S12, its symbol
$c_{br}$ is undefined in the paper, and its results (no cell reaches
80% at $N = 50$) flatly contradict S7 (hybrid 0.95 at $N = 16$). The
Reproducibility section's seeding description (`set.seed(20260507)`,
`mclapply`, "9-core") is this prototype's, not the sweeps' (which use
`set.seed(42 + id)` and `furrr` with `seed = TRUE`).

*Remedy.* Delete the section or move it to a clearly labeled appendix
with its own parameters, and rewrite Reproducibility for the sweeps
actually run.

### 4.11 Protocol, drivers and manuscript disagree

*Label.* INSPECTED. The disagreement covers grids, replicate counts and
provenance.

*Location.* ADEMP protocol `00-ademp-protocol.md`; Methods lines
265-325; README; Reproducibility.

*Problem.* The protocol's table gives grids that match none of the
drivers (S2 five half-lives versus nine; S6 $k \in \{1,2,3,4,6\}$ versus
$\{2,...,16\}$; S7, S8, S10 all differ; S11 is described as a $4 \times
4$ half-life factorial but is a $2 \times 2$ carryover-by-AR(1)
factorial; S12 "single cell" versus four). Replicates are 1,000 in S1
and 300 in S9, not 500. The protocol cites `01-base-seed.txt` and
`02-deviations-log.md`; neither exists in the sweep directory. Methods
cites `docs/sensitivity-sweep-plan-2026-04-17.md`, which does not exist
in the repository. The artifact modification times (18:15 to 19:07 on
2026-05-07) span about 55 minutes, not the stated 32. The cover letter
calls the study ADEMP-compliant; the deliverables listed in the
protocol (bias, coverage, per-sweep ADEMP tables) are absent.

*Remedy.* Bring the protocol into agreement with the drivers (with a
deviations log), or state that the protocol was superseded and by what.
Record replicate counts per sweep in a table.

### 4.12 Misattributed and overstated use of the prior literature

*Label.* INSPECTED against the retrieved sources.

*Location.* Methods lines 174-176 and 204-206; S8 lines 641-645; S2
lines 448-452; S3 lines 485-489; Discussion 840-896.

*Problem.* (a) The mean-moderation architecture is attributed to
Hendrickson et al. (2020); that paper contains only the covariance
encoding, and the mean-shift encoding is this project's addition
(paper 01). (b) `nof1power` is cited to Hendrickson; it is a vendored
project package. (c) Wang and Schork (2019) studied one subject with
400 observations and no parallel group; they found that serial
correlation *reduces* the within-subject design's power and that
stacking mitigates it. The S8 attribution ("the N-of-1 design benefits
from strong serial correlation ... while the parallel-group design
remains relatively sensitive") is not in that paper. The S9 attribution
("more frequent alternation benefits correlation-modeled analyses but
is problematic for absolute-effect inference") is at best a loose
paraphrase of their remark that the practicality of frequent
alternation under carryover "is an open question." The Table 1 values
quoted at lines 876-882 could not be checked (tables are not in the
PMC text) and are UNVERIFIED. (d) Blackston's parallel RCT used one
observation per patient, so its immunity to carryover is by
construction, not an "observation." (e) Gärtner et al. (2023) compare
analysis methods (including a carryover-adjusted parametric model)
under carryover; they do not study decay-form misspecification, so the
S3 precedent is misplaced. (f) The Hendrickson power collapse at
$t_{1/2} = 0.1$ weeks is, by the project's own `docs/40` Section 5 and
`docs/36`, an artifact of the published step rule under which the
interaction estimand is exactly zero; the manuscript presents it as a
genuine finding ("Both findings are correct"). (g) Carrozzo's
comparator is a two-arm *crossover* RCT, not a parallel RCT.

*Remedy.* Correct each attribution; cite Percha (2019) for S6 and S9;
cite Jones and Kenward or Senn (2002) for the two-period crossover
analysis recommendation in S12.

## 5 Minor comments

1. Lines 2-4, title: the design is "hybrid" in lower case here and
   "Hybrid" (capitalized, per NOTATION.md Part 3) nowhere.
2. Line 45: "moderate treatment main effects" is not defined; the
   baseline ceiling is 2 outcome units with component SDs 0.5 to 1.0.
3. Lines 52-53 and README: "twelve sweeps against the pmsimstats-ng
   tidyverse pipeline" but the engine is `nof1power`; say which.
4. Line 64 and NOTATION.md rule 4: $N = 35$ is a total split 18/17
   across two paths (`simulation.R` lines 142-147); state the per-path
   and per-arm counts, since an RCT with 17 per arm is unusually small
   for the claims made.
5. Lines 67-68, Abstract: "robust to AR(1) serial correlation that
   inflates RCT Type I error" presupposes comment 4.5.
6. Lines 170-171: "MCSE below 1%" is not met in S1 to S11 at 500
   replicates for power near 0.5 (MCSE 0.022).
7. Lines 212-213: the three components are labeled "time-variant,
   pharm-biomarker, bio-response." The `pb` component is the
   placebo-belief (expectancy) response, scaled by expectancy
   (`generate_data.R` lines 152 and 233); "pharm-biomarker" is wrong.
   Use NOTATION.md's $TV$, $PB$, $BR$ names and STYLE.md's
   natural-history, placebo-response and drug-response terms.
8. Lines 215-217: the Gompertz curve is written with code names
   (`maxr`, `disp`, `rate`); NOTATION.md requires $m$, $d$, $r$ in
   equations. Lines 392, 688-689, 754 use `max[br]`, `dgp\_architecture`,
   `halflife\_days` inside mathematics (rule 3).
9. Lines 222-223: the carryover recursion compounds (the previous
   decayed mean is decayed again with cumulative $t_{sd}$), so the
   realized half-life is shorter than the nominal one (`docs/53`).
   State this, since the day-scale axes are interpreted as
   pharmacokinetic.
10. Line 233 onward, Units: the exemption in `notation-lint.pl` was
    written for the retired paper 04 v1 day scale; confirm with the
    compendium whether paper 05 is to stay on days.
11. Lines 258-260: with `useDE = TRUE`, `De` is in the hybrid model but
    not the RCT model (expectancy constant); say so, since the hybrid's
    open-label visits then carry an expectancy adjustment the RCT lacks.
12. Line 307: "$c_{bm} = 0$" for the S12 null; the driver sets the
    drug-response ceiling to 0 and keeps $c_{bm} = 0.3$.
13. Lines 360-367 and the drivers: the hybrid is labeled `Nof1` in S1,
    S4, S5, S7, S8, S11 and `Hybrid` in S2, S10, S12. Use one label.
14. Line 398: "$|\delta| \ge 2$" with $\delta$ reserved in NOTATION.md
    for standardized bias; $\tau$ (line 500) and $\sigma$ collide with
    the heterogeneity variance and the residual SD; $k$ is used for
    timepoints and is also the Weibull shape in the same paper.
15. Line 413: "95% Wald CI [0.050, 0.082]" covers 0.05 only at its
    lower endpoint; "at the margin" understates it.
16. Lines 449-452: "largely immaterial to carryover" is Blackston by
    construction (comment 4.12 d).
17. Lines 470-475: "interacts poorly with the 'power' form" is the `NaN`
    bug of comment 4.7.
18. Line 520: "$\tau^2$ is part of $\mathrm{Var}(\hat\Delta)$" uses
    $\Delta$, reserved for the E9 contrast.
19. Line 576: "dead heat"; line 631 "warrants a closer look"; line 698
    "reassuring rather than surprising"; lines 438, 735, 811 "immune";
    line 416 "It appears that." STYLE.md forbids metaphor and
    characterization without the number.
20. Lines 639-640: $1/(1 + (k-1)\rho)$ is the compound-symmetry design
    effect, not the AR(1) one.
21. Lines 680 and 686-695: the S10 figure caption says
    "biomarker-interaction power" but the sweep stores only the main
    effect (`target_coef = "main"`), and the text says so. Fix the
    caption and the `plot-all.R` title.
22. Lines 701 and 845: "coupling"; line 910: "ceiling." STYLE.md
    requires the standard terms in prose and a Glossary section, which
    the manuscript lacks.
23. Lines 767-775: the two-period crossover has three rows per subject
    (baseline included), so "two per subject" is inaccurate; the
    recommendation to drop `corCAR1` is sound but untested here. A
    one-line diagnostic (refit S12 XO without `corCAR1`) would settle
    it.
24. Lines 842-853: the Hendrickson numbers 0.25/0.94/0.96/0.95 and
    0.95 to 0.18 are UNVERIFIED (figure values); give the figure and
    panel.
25. Line 877: "$1p \times 2i$" is Wang and Schork's notation for one
    period per intervention; define it.
26. Lines 925-929: the proposed `op$use_corCAR1` option is a one-line
    change; either do it or drop the limitation.
27. Lines 950-973: "three N-of-1 design variants" but only two are
    N-of-1 (hybrid, alternating); the crossover is the third. The
    upper-bound/lower-bound claims for unevaluated designs are
    assertions, and `docs/39` Section 4 shows more switches can lower
    power under CS with carryover, contradicting the lower-bound claim
    for multi-cycle designs.
28. Lines 1058-1063 and title page: the repository URL is a
    placeholder; "openly available" is not yet true.
29. README build line renders with `rmarkdown::render()`, which the
    compendium forbids (`tools/render.sh`).
30. `claudecode.tex` is loaded in the header (line 16) and is an inert
    leftover per CLAUDE.md; `lineno` with `mathlines` is a submission
    convenience that should not ship in the compendium copy.
31. Sentences over 30 words occur in the Abstract (lines 64-69) and
    throughout; semicolons join clauses at lines 143, 380, 1009.
32. US spelling is used in the Rmd; the title page and cover letter use
    "organisation," "favoured," "programme."

## 6 Missing and questionable references

| Location | Issue | Suggested source |
|---|---|---|
| Lines 174-176, 204-206 | `nof1power` and the mean-moderation architecture attributed to Hendrickson 2020 | Cite paper 01 of this compendium for mean moderation; cite the package by name and version with no external author |
| S6, S9 (lines 288-298) | Precedents for cycles and period length omit the closest prior simulation | Percha B, et al. Designing robust N-of-1 studies for precision medicine. *J Med Internet Res* 2019;21(4):e12641. doi:10.2196/12641 |
| S8 (lines 641-645) | Claim about parallel versus within-subject sensitivity to AR(1) attributed to Wang and Schork, who had no parallel arm | Remove the attribution; for ABAB versus AB under autocorrelation cite Lin SX, et al. *Br J Math Stat Psychol* 2016;69:276-290. doi:10.1111/bmsp.12071 |
| S3 (lines 485-489) | Gärtner 2023 cited for decay-form misspecification; it compares analysis methods under carryover | Keep Gärtner for carryover-adjusted analysis; for decay families cite the project's `docs/62` or drop the precedent |
| S2 (lines 449-452) | Blackston RCT immunity presented as a finding | State that the RCT had one observation per patient; cite Blackston accurately |
| S12 (lines 770-775) | Recommendation for two-period crossover analysis has no citation | Jones B, Kenward MG. *Design and Analysis of Cross-Over Trials*, 3rd ed. CRC Press, 2014; Senn S. *Cross-over Trials in Clinical Research*, 2nd ed. Wiley, 2002 |
| Lines 128-131, 283-285, 512, 543, 570 | "Senn 2018" in prose; the entry is 2019 (SMMR 28(2):372-383, epub 2017) | Use 2019 consistently |
| Lines 606-613 | Carrozzo's comparator described as "two-arm RCT" | It is a two-arm crossover RCT; revise the sentence |
| Lines 842-853 | Hendrickson power values quoted without figure reference; carryover collapse presented as genuine | Add figure and panel; cite `docs/36`/`docs/40` for the step-rule artifact |
| Line 932 and Limitations | Small-sample Wald anticonservatism in `lme` uncited | Kenward MG, Roger JH. *Biometrics* 1997;53:983-997; or the CR2 reference used in paper 10 |
| `references.bib` entries `araujoComparisonFourMethods2014`, `kravitzSinglepatientNof1Trials2014`, `chenMethodologicalReviewRandomised2024`, and nine others | Keys name one author, fields name another, with "retained for compatibility" notes; `araujo...` note says "Chen, Chen and Xia" but PubMed lists two authors (Chen X, Chen P) | Rename keys to match authors; add Araujo A, Julious S, Senn S. *PLoS ONE* 2016;11(12):e0167167 as its own entry |
| `references.bib` | 30 of 47 entries are uncited (Raskind, Khazaie, Miller, Daza, Hussey, Matthews, PSU lesson, Clinical Leader news item, Tang and Landes, Shi and Ye, Senn 2024, Hawksworth 2024, Chatters 2024, others) | Prune to cited entries; two web pages and a news report are not citable sources for a statistics journal |
| `references.bib` `giaiIndividualizedNet2023` | No volume or pages | Complete the record |
| `kravitzIntroductionNof1Trials2014`, `schmidStatisticalDesignAnalytic2014` | Notes say authorship "not checked against the PDF" | Check and remove the notes |

## 7 Suggestions for strengthening the paper

1. Rebuild the study around a defined estimand and a faithful Hybrid,
   then rerun under refusal of non-PD matrices with the feasibility
   bound computed per cell.
2. Cross the analysis half-life with the DGP half-life in every
   carryover sweep, so that each sweep separates the cost of unmodeled
   carryover from the cost of carryover itself.
3. Replace the RCT analysis with a standard MMRM arm contrast and show
   its null calibration before any comparison.
4. Report bias, empirical and model SE, coverage and $\kappa$ for every
   cell; this is what distinguishes an ADEMP study from a power table.
5. Lower the baseline effect so that the hybrid is off the ceiling, as
   the Limitations already recommend; most of S4, S5, S8 and S10 are
   uninformative at power 1.000.
6. Add a treatment-heterogeneity component so that Senn's decomposition
   is actually exercised, and sweep it.
7. Connect the main-effect findings to the project's interaction-focused
   design work (`docs/39`): the switch-count and spacing results there
   give closed-form predictions that the main-effect sweeps could test.
8. Generate the entire Results section from the artifacts; keep the
   drivers' grids, replicate counts and seeds in one machine-readable
   table that the manuscript reads.
9. Add the Glossary section STYLE.md requires and run the STYLE.md grep
   and `notation-lint.pl` before rendering (the linter passes now; the
   grep does not).

## 8 Scope of this review

**Verified (executed).** Per-cell power, MCSE, mean estimate,
empirical SE, mean model SE and singular-fit fraction for all twelve
stored artifacts (`sensitivity-01` to `-12.rds`), compared cell by cell
with the prose; the stored `N`, `carryover_t1half` and `use_DE`
columns per sweep; positive-definiteness and minimum eigenvalue of the
driver's covariance matrices for the S10 ($c_{bm} \in \{0.4, 0.6\}$,
$t_{1/2} \in \{0, 3\}$ days) and S8 ($\rho \in \{0.7, 0.9\}$) cells
under the current package; the `NaN` failure of the power decay form on
the hybrid; the realized Gompertz means at the baseline parameters; the
2-MCSE arithmetic for S12; `perl tools/notation-lint.pl` (clean) and
the STYLE.md grep (three hits); the existence or absence of files cited
by the manuscript and protocol; artifact modification times.

**Inspected (read).** `report.Rmd` in full, including YAML and chunks;
`README.md`, `title-page.md`, `cover-letter.md`, `references.bib`,
`planning/27-cycle-period-design-sweep-plan.md`; `00-common.R`,
`00-ademp-protocol.md`, `00-smoke-test.R`, drivers `01` to `12`,
`run-all.R`, `plot-all.R`, `05-design-prototype.R` and its summary;
`implementations/nof1power/R/analysis.R`, `generate_data.R`,
`simulation.R`, `trial_design.R`, `carryover_models.R`,
`data_processing.R`, `parallel_utilities.R`; `NOTATION.md`,
`STYLE.md`, `docs/39`, `docs/40`, the project `CLAUDE.md`.

**Not checked.** The twelve figure PDFs were not opened. No simulation
was rerun. The RCT bias mechanism (comment 4.5) was not tested by
refitting. The Wang and Schork Table 1 values and the Hendrickson
figure values quoted in the Discussion were not verified (tables and
figures are not in the PMC text). `bullets.Rmd` and the rendered
`report.pdf`/`report.tex` were not read. S13 (`13-cycle-period.R`) is
outside the manuscript and was not reviewed.

**Literature searches run (PubMed, 2026-10-28, no email or identifier
sent).** Citation lookups for Blackston 2019, Wang and Schork 2019,
Senn 2019, Zucker 2010, Morris 2019, Hendrickson 2020, Carrozzo 2025,
Gärtner 2023, Jiang 2025, Selukar 2025, Chen and Chen 2014, Shi and Ye
2024, Senn 2024, Hawksworth 2024, Araujo 2016, Percha 2019 (all
resolved; Jiang and Shi required a title search). Full text retrieved
for Wang and Schork 2019, Blackston 2019 and Hendrickson 2020. Keyword
searches: "(N-of-1 OR single-patient trial) AND simulation AND (power
OR sample size) AND (carryover OR autocorrelation OR serial
correlation)" from 2015 (9 hits, all screened); "(aggregated N-of-1 OR
series of N-of-1) AND (parallel group OR randomized controlled trial)
AND simulation AND (power OR efficiency) AND comparison" from 2015 (2
hits); "N-of-1 trial design optimization number of cycles period
length washout simulation power" from 2015 (0 hits). Web search was not
used; the Zotero MCP server was unavailable.
