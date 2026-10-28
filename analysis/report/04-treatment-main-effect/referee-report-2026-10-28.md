# Referee Report: Paper 04, N-of-1 versus RCT Power
*2026-10-28 07:35 PDT*

**Manuscript.** `analysis/report/04-treatment-main-effect/report.Rmd`
(2,081 lines, last substantive commit 2026-10-09). Driver:
`analysis/scripts/treatment-main-effect/production-run.R`. Stored
output: `analysis/data/derived_data/04-mainfx-production.rds` (run
2026-05-19, R 4.5.3). Target venue assumed: *Statistics in Medicine*
or *Trials*.

**Evidence labels.** VERIFIED: re-derived by hand or executed in R on
the stored output or the driver. INSPECTED: read in source and judged.
SUSPECTED: plausible from context, not confirmed. Every finding below
carries one of these labels.

## 1 Summary of the manuscript

The paper reports a Monte Carlo comparison of the power to detect a
treatment main effect $\beta_D$ in two designs at a matched participant
count of $N = 35$: an "aggregated N-of-1 hybrid" in which each
participant completes eight one-week periods (four on prazosin, four on
placebo, in a random balanced order), analyzed by `nlme::lme` with a
random intercept and `corCAR1` residual correlation; and a
parallel-group RCT with a single end-of-study measurement analyzed by
ANCOVA on baseline. The data-generating process (DGP) is linear:
outcome = baseline + $\beta_D D$ + carryover + AR(1) noise, with
between-participant SD 1.5, within-participant SD 2.3, lag-one
correlation 0.3, and an additive exponential carryover with half-life
3 days reaching only the placebo period that immediately follows a drug
period. Carryover-contaminated placebo periods are removed before
fitting, using a flag written by the DGP. The effect-size grid is
$\beta_D \in \{0, -0.5, -1, -2\}$ nightmares per week; a carryover
sweep varies the half-life over $\{0.5, 1, 3, 5, 7\}$ days; 2,000
replicates per cell.

Stated contributions, in the manuscript's words: an "ADEMP-compliant"
simulation showing that "the aggregated N-of-1 hybrid dominated the
parallel-group RCT across the full effect-size range" (power 1.000
versus 0.675 at $\beta_D = -2$; 0.830 versus 0.230 at $-1$), that bias
is negligible, that coverage is 0.940 to 0.949 (N-of-1) and 0.939
(RCT), and that "the period-exclusion strategy for carryover
management proved robust" because N-of-1 power, bias and SE were
constant across the five half-life cells. The Introduction adds a
tutorial-style survey of closed-form power for series of N-of-1 trials
(Senn, Zucker) and a review of prior simulation comparisons.

## 2 Overall assessment and recommendation

**Recommendation: reject in its present form.** Resubmission after the
rebuild the authors' own working protocol describes
(`00-ademp-protocol.md`, version of 2026-10-08 and 2026-10-09) would be
a different paper and should be reviewed as one.

Justification, in order of weight.

1. The study the paper describes is not the study the code runs. The
   N-of-1 arm is labeled a Hendrickson hybrid in the title-level claims
   but is an eight-period multi-crossover; the inference is labeled
   Satterthwaite but uses `nlme` containment degrees of freedom with no
   correction applied; the carryover sweep is reported at $\beta_D =
   -2$ but was run at $-1$; the half-life "robustness" finding is an
   algebraic identity of the oracle exclusion rule, not a result.
   Every one of these is VERIFIED below.
2. The headline comparison is reproducible in closed form from the
   stated variance components in under a second (VERIFIED, Section 4
   and Major comment 6). The simulation therefore adds no information,
   and the comparison it makes is unfair to the RCT in two ways the
   authors only partly acknowledge (one measurement versus about six,
   and an error-free baseline covariate that no trial can observe).
3. The Results and Discussion quote numbers from sensitivity sweeps
   (S6, S7, S8, S11, S12) that the production run never executed and
   whose only source is the workspace the authors themselves have
   marked "SUPERSEDED, NOT AUTHORITATIVE. Do not cite any number from
   this document" (VERIFIED).
4. Several numerical statements contradict the stored output
   (VERIFIED): the RCT Type I interval does not contain 5%; the stated
   null-cell empirical SE (0.28) is 0.2455; the stated intraclass
   correlation (0.69) is 0.30; the "paired MCSE" that is said to be
   smaller than the independent-sum formula is computed with the
   independent-sum formula.
5. The authors' own current protocol (`00-ademp-protocol.md`,
   Section 1) lists seven defects of this simulation and concludes
   that a rebuild, not a revision, is required. I agree with that
   assessment and with its reasons. The manuscript cites that protocol
   as the source of its design.
6. The literature review contains material errors: the null NEJM trial
   of prazosin (Raskind et al. 2018) is cited as evidence of efficacy;
   the numbers quoted from Tang and Landes (2020) are taken from the
   wrong rows and the wrong table; one single-author paper is cited as
   two-author; two companion documents cited by path do not exist in
   the repository.

The paper has real strengths that the rebuild should keep: the ADEMP
scaffolding, MCSEs on every performance measure, a single uniform
replicate count, the L'Ecuyer stream discipline, the Morris-style
tables, and a candid Limitations paragraph on the observation-count
confound. These are matters of form; the substance must be redone.

## 3 Significance and novelty

**What the paper claims is new.** "Empirical power estimates across a
clinically calibrated parameter range" for prazosin in PTSD, extending
Blackston et al. (2019), with carryover handled by period exclusion and
serial correlation by `corCAR1`.

**What is actually new.** Nothing that survives scrutiny. The DGP has a
constant treatment effect, no treatment-by-participant variance, no
time trend, no placebo response, and carryover that is removed by an
oracle flag before analysis. Once the contaminated rows are gone the
analysis model is a correctly specified Gaussian linear mixed model
with known structure, and its standard error is a deterministic
function of the design: averaging the GLS information over the 70
balanced sequences gives SE 0.2455 with all eight periods (the null
cell) and 0.336 after exclusion (the alternative cells), against 0.2455
and 0.342 observed (VERIFIED, scratch script, under one second). The
RCT standard error is $2.3\sqrt{1/17 + 1/18} = 0.778$ (observed 0.783).
Noncentral-$t$ power from these SEs is 0.315, 0.840, 1.000 (N-of-1) and
0.096, 0.239, 0.703 (RCT) at $\beta_D = -0.5, -1, -2$, against the
simulated 0.312, 0.830, 1.000 and 0.105, 0.230, 0.675. The paper's
central table is a noisy rendering of this arithmetic. The result that
a within-participant contrast with six observations per person beats a
between-participant contrast with one observation per person, at
matched $N$, is a textbook consequence of the variance decomposition
the Introduction itself presents (Senn 2019, equation for
$\mathrm{Var}(\hat\Delta)$).

**Competing and overlapping literature.** The same comparison, with
more content, exists in several places:

- Blackston JW, Chapple AG, McGree JM, McDonald S, Nikles J.
  Comparison of aggregated N-of-1 trials with parallel and crossover
  randomized controlled trials using simulation studies. *Healthcare*
  2019;7(4):137. doi:10.3390/healthcare7040137.
  Four variance structures, carryover of three sizes, selection bias,
  5,000 replicates, and (unlike the present paper) a Type I error
  study under unmodeled carryover and a carryover-adjusted analysis.
  Note that Blackston's RCT also uses one observation per participant
  drawn from the first cycle, so the present paper inherits rather
  than resolves the observation-count confound.
- Carrozzo AE, Zimmermann G, Bathke AC, Neunhaeuserer D, Niebauer J,
  Kulnik ST. Two-arm crossover randomized controlled trial versus
  meta-analysis of N-of-1 studies: comparison of statistical efficiency
  in determining an intervention effect. *Biometrical Journal*
  2025;67(2):e70045. doi:10.1002/bimj.70045.
  Two-arm RCT, crossover and N-of-1 meta-analysis; heterogeneity and
  carryover varied; sequential aggregation.
- Chen X, Chen P. A comparison of four methods for the analysis of
  N-of-1 trials. *PLoS ONE* 2014;9(2):e87752.
  doi:10.1371/journal.pone.0087752.
  Compound-symmetry and AR(1) covariance, 0% and 20% carryover, 3,000
  replicates; includes the mixed model with a carryover regressor that
  the present paper defers to "sweep S11".
- Percha B, Baskerville EB, Johnson M, Dudley JT, Zimmerman N.
  Designing robust N-of-1 studies for precision medicine: simulation
  study and design recommendations. *J Med Internet Res*
  2019;21(4):e12641. doi:10.2196/12641.
  A stochastic time-series simulator with carryover, wash-in and
  autocorrelation, and design recommendations on block count, order,
  period length and sampling frequency. Not cited.
- Wang Y, Schork NJ. Power and design issues in crossover-based N-of-1
  clinical trials with fixed data collection periods. *Healthcare*
  2019;7(3):84. doi:10.3390/healthcare7030084.
  Analytical power under AR(1), heteroscedasticity, washout and
  carryover. Cited, but it already supplies the closed forms this
  paper replaces with simulation.
- Senn S. Sample size considerations for n-of-1 trials. *Stat Methods
  Med Res* 2019;28(2):372-383. doi:10.1177/0962280217726801; and Yang
  J, Steingrimsson JA, Schmid CH. Sample size calculations for
  n-of-1 trials. arXiv:2110.08970 (2021). Closed-form planning for
  series of N-of-1 trials, including the individual-effect precision
  this paper's DGP cannot address because $\tau^2 = 0$. The latter is
  not cited.
- Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
  aggregated N-of-1 trial designs for predictive biomarker validation.
  *Front Digit Health* 2020;2:13. doi:10.3389/fdgth.2020.00013.
  The design the paper claims to simulate, with the three-component
  Gompertz DGP the v1 protocol promised and the driver omits.

**Contemporariness.** The reference list reaches 2025 (Carrozzo; Jiang,
Arnold and Betensky, *Stat Med* 2025, doi:10.1002/sim.70177; Selukar,
Prince and May, *Clin Trials* 2025, doi:10.1177/17407745241304284),
so the survey is current. The simulation is not: it is a pilot-grade
linear model whose own protocol has superseded it.

**What would clear a top-tier bar.** The rebuild in the current
protocol, executed: the Hendrickson hybrid, OL+BDC, CO and a
multi-visit RCT simulated from one component DGP with placebo
response, natural course, Gompertz onset, treatment-effect
heterogeneity ($\tau^2 > 0$) and measured-with-error baseline; carryover
handled by rules a trial can apply (schedule-based exclusion, or
$D_{bc}$ at a pre-specified analysis half-life), never by a DGP flag;
closed-form probability limits and sandwich variances computed first
and the simulation used to validate them; estimands defined on the
exposure scale so that bias and coverage are meaningful; and a
matched-burden RCT arm so the design-family effect is separated from
the observation-count effect. That study would be publishable in
*Statistics in Medicine*. The present one would not be publishable in
*Trials*.

## 4 Major comments

Each comment gives the location, the problem, why it matters, and a
concrete remedy. Correctness items come first.

### M1. The N-of-1 arm is not the hybrid design the paper names (VERIFIED)

*Location.* Abstract ("aggregated N-of-1 hybrid", three times); Results
"Headline findings", "Power curves", "Optimal design conditions";
Conclusions; 21 occurrences of "hybrid" in the Rmd. Driver lines 80 to
85.

*Problem.* `generate_nof1_patient()` builds `seq_base <- c(rep(1L, 4),
rep(0L, 4))` and permutes it with `sample()`. That is a randomized
4 + 4 multi-crossover of one-week periods. The Hendrickson hybrid
(open-label titration, blinded discontinuation, brief crossover; visits
at weeks 4, 8, 9, 10, 11, 12, 16, 20) has no open-label phase here, no
expectancy covariate, and no titration. Methods (Section "Simulation
design") describes the multi-crossover correctly, so the paper
contradicts itself. The v1 protocol (`00-ademp-protocol-v1.md`, "Trial
design (N-of-1 arm)") specified the hybrid; the driver never
implemented it.

*Why it matters.* The hybrid's defining property, that most of its
power comes from within-path contrasts protected only by the time and
expectancy model, is exactly what the compendium's main-effect work
(`docs/58`) is about. A reader who takes the paper at its word will
attribute the multi-crossover's power to the hybrid.

*Remedy.* Either simulate the hybrid (the protocol's Table "Designs"
gives the paths and schedule, and `hendrickson-problems/04-power-
simulation.R` has them in code) or remove "hybrid" from every claim
and call the design what it is. The former is the rebuild; the latter
is the minimum for honesty.

### M2. Oracle carryover exclusion; vacuous S1; `Dbc` misdescribed (VERIFIED)

*Location.* Methods "Statistical analysis models"; Results "Carryover
half-life sensitivity (S1)", "Bias and coverage", "Optimal design
conditions"; Discussion "Comparison of carryover handling strategies";
Conclusions, second point. Driver lines 94 to 110, 126 to 127.

*Problem, part (a).* The analysis drops rows where `is_carryover` is
`TRUE`, a column the DGP sets from the true carryover term. A trial
cannot observe this. Because exclusion removes every contaminated row
whatever the half-life, the half-life cannot enter the estimating
equations, and the five S1 cells are identical to the fourth decimal
(bias 0.01059, SE 0.3419, coverage 0.9395 in all five; stored output).
The paper states this mechanism itself ("the true half-life no longer
enters the estimating equations") and then presents the invariance as
a finding about "robustness to carryover half-life misspecification".
It is not a finding; it is the definition of the rule. At $t_{1/2} =
0.5$ day the residual carryover is $2^{-14}$ of the effect (VERIFIED,
multiplier 0.00006) and the period is still discarded.

*Problem, part (b).* Under $\beta_D = 0$ the carryover term is zero, so
`carryover_adj != 0` is `FALSE` everywhere and nothing is excluded.
The null cell therefore analyzes all eight periods per participant
while every alternative cell analyzes about six (mean excluded periods
2.000 over the 70 sequences; VERIFIED), which is why the null-cell
empirical SE is 0.2455 and the alternative-cell SE is 0.3419. The Type
I error and coverage reported for "the" N-of-1 analysis are the size of
one procedure and the coverage of another.

*Problem, part (c).* Methods calls `Dbc` "the continuous
exposure-decayed drug indicator" and the Discussion says the model
"already carries a continuous drug-exposure regressor" that "absorbs
exponential-decay carryover". In the driver `Dbc` is the binary
sequence vector (line 84) and is never decayed. The compendium's
$D_{bc}$ with pre-exposure guard and decay at an analysis half-life
(`docs/56`, NOTATION.md) is not in this code.

*Why it matters.* Carryover handling is the one methodological question
the paper poses that the closed form does not answer, and the paper
answers it with a rule that uses the truth. The Discussion's claim that
"period exclusion ... prevents bias" is true only because the DGP tells
the analyst which rows to drop. The Conclusions' second point is
unsupported.

*Remedy.* Implement exclusion by schedule (the protocol's XS analysis:
drop the first off-drug period after each discontinuation, whatever the
data say), apply it identically in null and alternative cells, and
compare it with a $D_{bc}$ analysis at a fixed $t_{1/2}^a$ and with an
explicit carryover regressor (Chen and Chen's M3). Report the oracle
rule, if at all, as a labeled oracle (XO). Define `Dbc` correctly or
rename the column `D`.

### M3. The inference is not Satterthwaite; no correction occurred (VERIFIED)

*Location.* Methods "Statistical analysis models" ("REML Wald
$t$-statistic with Satterthwaite denominator degrees of freedom from
`nlme`, which corrects the small-sample anti-conservatism we observed
in the preliminary runs (Type I = 0.061 ...)"); table footnotes;
Conclusions third point. Driver lines 158 to 183.

*Problem.* `nlme::lme` reports containment (between-within) degrees of
freedom; it has no Satterthwaite option. The driver reads
`st["Dbc", "DF"]`, computes a $p$-value from it, and stores it as
`p_kr`. That is the same reference distribution `nlme` already used, so
`nof1_power_wald` equals `nof1_power_kr` in all nine cells (stored
output, VERIFIED), and the "correction" is the identity. The driver's
header still says "KR correction applied ... via pbkrtest", loads
`pbkrtest`, fits a reduced model (lines 159 to 170) and never uses it.
The 0.061 the Methods attribute to the pre-correction N-of-1 test is,
in the stored output, the RCT's Type I rate; the N-of-1 null rate is
0.051 with or without the "correction".

*Why it matters.* With about 170 containment df the practical effect is
small here, which is why the Type I rate is near nominal. But the
claim is false as written, the remediation plan's Problem 4 asked for
a real small-sample correction, and the compendium's standard for
every other paper is CR2 with Satterthwaite df via `clubSandwich`.

*Remedy.* State that the test uses `nlme` containment df, or implement
CR2 with Satterthwaite df and report both. Delete the sentence about
correcting 0.061. Remove `pbkrtest` and the unused reduced fit.

### M4. The carryover sweep ran at $\beta_D = -1$, not $-2$ (VERIFIED)

*Location.* Methods ("sweep S1 varies $t_{1/2}$ ... at fixed $\beta_D
= -2.0$"); Results "Carryover Decay Model Sensitivity Analysis" ("at
the primary cell that model yields N-of-1 power of `r ...`"); Figure
`power-heatmap`; driver line 354 (`run_cell(-1.0, hl, ...)`).

*Problem.* The stored `s1_carryover` table has `true_effect = -1` in
every row. The inline "primary cell" number in the decay-sensitivity
paragraph is therefore 83.0%, not the primary cell's 100%. The heatmap
(`p2b`) is built from `s1_carry` alone: one column ($-1$) and five
rows that all read "+60", so the caption's "blue regions ... purple
regions ... white regions" describe a figure that cannot exist. The
object `power_wide_es` that would supply the other columns is built
and never used. The v1 protocol, moreover, named the effect-size sweep
S1 and the carryover sweep S2; the driver and the paper call the
carryover sweep S1.

*Why it matters.* A reader cannot tell which effect size any S1 number
refers to, and the one figure that promised a two-way view is empty.

*Remedy.* Correct the Methods to $-1$ (or rerun at $-2$, where the
N-of-1 arm is saturated and the sweep would show nothing), fix the
inline reference, rebuild the heatmap from both tables or drop it, and
reconcile the sweep labels with the protocol in force.

### M5. One RNG stream for all cells; false paired-MCSE claim (VERIFIED)

*Location.* Methods "Simulation size and random-number discipline";
Results "Headline findings" ("the paired MCSE is smaller than the
square-root-sum of the two marginal MCSEs, reflecting the variance
reduction gained from sharing the patient pool"); Rmd chunk
`load-simulation-results`, lines 650 to 651. Driver lines 35 to 36,
278 to 283.

*Problem, part (a).* `parallel::mclapply(..., mc.set.seed = TRUE)` does
not advance the master L'Ecuyer stream between calls. I ran two
consecutive `mclapply` calls after one `set.seed(2024L)` and obtained
identical draws (scratch script). The driver seeds once and calls
`mclapply` nine times, so every cell used the same 2,000 random-number
streams. The stored output shows the symptom: `rct_bias` is 0.01461 in
all nine cells including the null, and all alternative-cell N-of-1
statistics are identical to four decimals. Common random numbers are a
legitimate variance-reduction device, but here they are undisclosed,
unintended, and they make the cross-cell comparisons (and the "exceeded
five combined MCSEs" argument in the Conclusions) dependent rather
than independent.

*Problem, part (b).* The Rmd computes `paired_power_diff_mcse <-
sqrt(n_of_1_power_mcse^2 + rct_power_mcse^2)`, which is the
independence formula, and the prose asserts that this number "is
smaller than the square-root-sum". It is the square-root-sum. The
within-replicate covariance cannot be computed because per-replicate
rejections are not saved (see M9).

*Why it matters.* The reproducibility section is the one part of the
paper that reads as careful, and it is wrong in both directions: the
stream discipline described was not achieved, and the paired-design
benefit described was not computed.

*Remedy.* Advance the stream per cell (`.Random.seed <-
parallel::nextRNGStream(.Random.seed)` before each `mclapply`, or seed
each cell from a recorded table of seeds as the protocol's Section 7.5
proposes). Save per-replicate results. Compute the paired MCSE from the
replicate-level agreement of the two rejection indicators, and delete
the sentence until it is true.

### M6. Unfair RCT comparator; the result is a closed form (VERIFIED)

*Location.* Methods "Simulation design" (RCT arm); Discussion
"Limitations", fourth paragraph; Discussion "Statistical efficiency".
Driver lines 197 to 227.

*Problem.* (i) The ANCOVA covariate is the participant's true random
intercept `baselines[i]`, measured without error. With this covariate
the between-participant variance is removed entirely from the RCT, a
limit no trial attains; the protocol (Section 3.3, "Baseline measured
with error") already flags this. (ii) The RCT has one post-baseline
measurement against about six analyzed N-of-1 periods. The paper's
Limitation 4 concedes (ii) candidly but then reports the unmatched
comparison as the headline. (iii) Neither arm has a placebo response or
natural course, so the one feature that distinguishes a blinded RCT
from an open-label-containing design is absent.

Given (i) to (iii), the entire Table 1 and Figure "power curves" are a
deterministic function of $\sigma_u = 1.5$, $\sigma_w = 2.3$, $\rho =
0.3$, $N = 35$ and the exclusion rule. I reproduced every power entry
within about one MCSE from the GLS information averaged over the 70
sequences (Section 3 above; scratch script, under one second).

*Why it matters.* A simulation whose answer is known before it runs
adds nothing to the closed form, and the closed form shows that the
"N-of-1 advantage" here is the ratio $0.778 / 0.336 = 2.3$ of two
standard errors, which is entirely an observation-count and
error-free-baseline effect. The paper's clinical framing ("N-of-1
trials achieve substantially higher statistical power ... when
carryover is managed by period exclusion") attributes to a design and
a carryover rule what is due to six measurements versus one.

*Remedy.* Add the matched-burden RCT (eight post-baseline visits, MMRM
or ANCOVA on the final visit) and the measured-with-error baseline
($r_{BL} = 0.8$ as the protocol proposes). Present the closed form
first and use the simulation to validate it, as the protocol's Section
5.1 requires. State in the Abstract that the comparison is at matched
participant count, not matched observation count.

### M7. Statements contradicted by the stored output (VERIFIED)

*Location and problem, item by item.*

1. Results "Type I Error Control": "Both confidence intervals include
   the nominal 5% level." The RCT rate is 0.0615 with MCSE 0.0054; the
   95% interval is 0.051 to 0.072. It excludes 0.05. The Conclusions'
   third point says the opposite of the Results ("marginally elevated
   Type I rate of 0.061") and attributes it to "known small-sample
   effects in `lm`". The ANCOVA $t$-test is exact under the DGP's
   Gaussian errors; there is no small-sample effect. The excess is 2.1
   MCSE and is most plausibly chance (SUSPECTED), compounded by the
   shared stream in M5. The paper should say so, not invent a
   mechanism.
2. Results "Power curves": "The N-of-1 empirical SE is calibrated to
   approximately 0.28 ... at the null cell." Stored: 0.2455. The same
   0.28 appears in "Preliminary validation run" as the achieved
   calibration; it was the target, not the result. The driver's own
   comment derives 0.276 from $\sqrt{2/(35 \times 4)}$, a formula that
   ignores the AR(1) structure and the time term.
3. Results "Between-patient variance": "$\sigma_\tau = 1.5$ ...
   calibrated to yield an intraclass correlation near 0.69." With
   $\sigma_w = 2.3$ the ICC is $2.25 / (2.25 + 5.29) = 0.30$. The value
   0.69 is what $\sigma_w = 1.0$ gives, the within-SD of the superseded
   short report. Methods adds "calibrated to yield $I^2 \approx 50\%$
   in the parallel analysis"; $I^2$ is a meta-analysis statistic that
   this pipeline never computes (the Rmd sets `i_squared_mean <-
   NA_real_`).
4. Coverage. Stored N-of-1 coverage at the alternative cells is 0.9395
   (MCSE 0.0053, 1.97 MCSE below 0.95) and RCT coverage 0.9385 (2.14
   MCSE below). The Results say coverage "tracks the nominal level";
   the Figure `fig-morris-coverage` caption attributes the shortfall
   to "small-sample Wald test conservatism" (conservatism would give
   over-coverage); the "Primary Power Comparison" paragraph attributes
   the N-of-1 shortfall to "systematic carryover bias ... (see Section
   bias)", and that section says bias is negligible (0.011). Three
   incompatible explanations of a 1 percentage-point shortfall that is
   probably noise.
5. Results "Bias and coverage": "relative error within $\pm 2\%$ for
   both arms". Stored: N-of-1 $-2.1\%$, RCT $-2.5\%$.
6. Results "Headline findings" and "Primary Power Comparison" quote
   the same two numbers twice in consecutive paragraphs with different
   wording and rounding; "Comprehensive Performance Summary" then
   quotes them a third time with a "Mean estimate" row computed from
   500 synthetic normal draws (`rnorm(500, ...)`, chunk
   `prepare-plot-data`), not from replicates.

*Remedy.* Regenerate every inline number from the replicate-level
output after M5 is fixed, delete the mechanistic explanations that the
data do not support, and remove the synthetic "distribution of
estimates" figure and its derived row (see also m12).

### M8. Results from sweeps never run; missing companion files (VERIFIED)

*Location.* Discussion "Comparison with Prior Simulation Studies":
"our hybrid reaches 80% power at $N \le 16$ while the parallel RCT
requires $N \approx 40$" (S7); "our S12 null calibration 0.10-0.13";
"the hybrid maintains Type I in [0.056, 0.062]" (S12); "the
40-participant aggregated analysis" (S8); "requires $k \ge 6$ to match
RCT power (S6)". Introduction "Prior Simulation Comparisons" and the
Discussion cite `docs/literature-review-nof1-vs-parallel-2026-04-17.md`
and `docs/06-blackston.md` by path.

*Problem.* The production output contains the M1 grid and the
carryover sweep only (two data frames, 1,068 bytes). Sweeps S6, S7, S8,
S11 and S12 exist only in `sim_workspace.RData`, the input of
`report_short.Rmd`, whose banner says "Do not cite any number from this
document" and whose generating script does not exist. The quoted
numbers describe a 40-participant design and a 1,500-replicate run,
neither of which is this paper's. Neither cited companion document is
in the repository (`find` over the whole tree returns nothing; the
legacy index has no entry). The subsection also promises "five
completed head-to-head comparisons" and delivers three.

*Why it matters.* These are the paper's only quantitative links to the
prior literature, and they rest on retracted numbers and missing
files.

*Remedy.* Delete every S6 to S12 number, or run those sweeps in the
production driver and cite the result. Restore the two companion
documents or remove the references.

### M9. Reproducibility and data-availability statements are false (VERIFIED)

*Location.* "Reproducibility" ("Each simulation driver sets
`set.seed(20260507)` at the top of the replicate loop and uses
`parallel::mc.set.seed` ... Per-replicate seeds derive from the master
seed by `+ rep_idx`"); "Data availability statement" ("Per-replicate
Monte Carlo outputs ... are versioned under `analysis/data/`");
Results "Between-patient variance" ("Per-replicate estimates of
$\hat\sigma_u^2$ are available in the `.rds` output"); Methods ("We
capture per-replicate `.Random.seed` states").

*Problem.* The driver seeds with `set.seed(2024L)` once, as Methods
says; the Reproducibility section describes a different scheme. The
driver captures `seed_state` per replicate (line 232) and then
discards it with everything else replicate-level: `run_cell()` returns
one summary row and `reps` is dropped. The `.rds` has no $\hat\sigma_u^2$,
no seeds, no per-replicate estimates. The "8-core `mclapply`" is
`detectCores() - 1`. The GitHub URL is a placeholder.

*Remedy.* Save the replicate-level list (it is small), write the
Reproducibility section from the code, and remove the three false
availability statements.

### M10. The literature review misstates its sources (VERIFIED unless marked)

1. Introduction "Clinical Application": "Prazosin ... has demonstrated
   efficacy for reducing trauma-related nightmares in multiple RCTs
   [@raskind2013; @raskind2018]." Raskind et al. 2018 (*N Engl J Med*
   378:507-517, doi:10.1056/NEJMoa1507598)
   randomized 304 veterans and found no difference from placebo on any
   of its three primary outcomes (PubMed abstract, VERIFIED). Citing
   it as evidence of efficacy is a factual error, and the null result
   is directly relevant to the paper's premise that the effect is
   $-2$ nightmares per week.
2. "Senn and Julious [-@sennAnalysisContinuousData2024]": the paper is
   single-authored (Senn S, *Trials* 2024;25:128,
   doi:10.1186/s13063-024-07964-7; VERIFIED).
3. "Chen, Chen, and Xia (2014)": the paper has two authors, Chen X and
   Chen P (doi:10.1371/journal.pone.0087752; VERIFIED). The bib entry's
   note says "verified authorship is Chen,
   Chen and Xia", which is itself wrong.
4. Tang and Landes (2020), Discussion subsection: "at $m = 4$, $\rho =
   0$ the detectable $\mu_{A-B}$ is $1.65\sigma$; at $\rho = 0.33$ it
   is $2.81\sigma$; at $\rho = 0.67$ it is $13.73\sigma$ ... At $m =
   12$, the same points are $0.52\sigma$, $1.25\sigma$, and
   $5.43\sigma$." Their Table 2 is indexed by $\rho \in \{0, 0.2, 0.4,
   0.6, 0.8\}$, not 0.33 and 0.67. At $m = 4$: 1.65 ($\rho = 0$), 4.08
   ($\rho = 0.4$), 13.73 ($\rho = 0.6$); the quoted 2.81 is the $m =
   5$, $\rho = 0.4$ entry. At $m = 12$: 0.77 ($\rho = 0$), 1.25 ($\rho
   = 0.4$), 1.86 ($\rho = 0.6$); the quoted 0.52 is a Table 1 margin
   of error and 5.43 is the $\rho = 0.8$ entry (PMC full text,
   VERIFIED). The "48-fold inflation" is therefore 8.3-fold at the
   rows that exist. The stated "one-sided 0.10 level" is their
   "one-sided 0.05". The Type I rates "roughly 0.22 ... roughly 0.06"
   are read from a figure I could not check (SUSPECTED).
5. Blackston et al.: the quoted 92/32/50% at $n = 30$, the Type I rates
   15% and 25% under carryover, and the crossover 0.07 are correct
   (PMC full text, VERIFIED). "Adjusting for carryover restored
   nominal Type I" overstates their sentence that results "were nearly
   identical when controlling for carryover effects" (SUSPECTED
   overstatement). The 150 and 100 sample sizes could not be checked
   in the extracted text.
6. "Zucker et al. [-@zucker2010] ... Their precision formula,
   $\mathrm{Var}(\hat\Delta) = (\tau^2 + \sigma_w^2/k)/n$": Zucker,
   Ruthazer and Schmid 2010 is a comparison of models on one published
   series; the abstract does not present a precision formula
   (SUSPECTED misattribution). The formula is Senn's with $2\sigma^2$
   renamed $\sigma_w^2$, and the paper presents the two as separate
   results.
7. Wang and Schork: "washout periods reduce effective sample size in
   proportion to the carryover half-life relative to the period
   length" is not in their abstract and reads like this paper's own
   model (SUSPECTED).
8. The 3-day primary half-life is justified by "prazosin's
   neuroadaptive effects on sleep architecture [@Khazaie2016Prazosin]".
   Khazaie et al. 2016 is a narrative piece in the *Iranian Journal of
   Psychiatry and Behavioral Sciences*; I found no half-life estimate
   in it (SUSPECTED unsupported). The compendium's own pharmacological
   range is 0.1 to 0.2 week (17 to 34 hours), and NOTATION.md requires
   longer values to be labeled behavioral scenarios. Three days is
   outside that range and is not labeled.
9. The calibration "8.0 nightmares per week, SD 1.5, within SD 2.3,
   effect $-2$" is attributed to "Raskind 2013 / Hendrickson 2020"
   (driver comment). Raskind 2013 used the CAPS recurrent distressing
   dreams item, not weekly counts; no source for these four numbers is
   given anywhere (SUSPECTED unsupported).

*Remedy.* Correct items 1 to 4 and 6 to 9 in the text; add a
calibration table with sources; label the 3-day half-life a behavioral
scenario or move the primary cell to 0.2 week.

### M11. The paper's own protocol disowns the study (VERIFIED)

*Location.* Methods "Estimands and study protocol" cites
`analysis/scripts/treatment-main-effect/00-ademp-protocol.md` as the
protocol that "records the parameter grid for the primary comparison
(M1) and the twelve sensitivity sweeps (S1-S12)". Results "Power
curves": "A thorough comparison against the full DGP of the working
protocol is deferred to the complete production run, and the absolute
power magnitudes reported here should be read in that light." Six
further deferrals to "the full production program" (S3, S11).

*Problem.* The file at that path is the 2026-10-08 rebuild protocol.
Its Section 1, "Why the rebuild", lists seven defects of this exact
driver (design, outcome model, no heterogeneity, oracle exclusion,
truncated day-scale carryover, mislabeled inference, error-free RCT
covariate) and concludes "a rebuild rather than a revision". My
findings M1 to M3 and M6 agree with it independently. The v1 protocol
it replaced promised a three-component Gompertz DGP, the hybrid, a
16-week five-visit RCT and twelve sweeps; the driver delivers a linear
DGP, a multi-crossover, a one-visit RCT and one sweep under a changed
label. The v1 protocol's promised `02-deviations-log.md` does not exist
for this study.

*Why it matters.* A manuscript cannot cite as its protocol a document
whose first section says the manuscript's simulation is not fit for
purpose. And a paper that calls its central result provisional pending
"the complete production run" is a pilot report.

*Remedy.* Run the rebuild. If a bridging paper is wanted first, say in
the title and Abstract that it is a linear-DGP pilot, cite the v1
protocol as its protocol, and record every deviation.

## 5 Minor comments

Anchored to section or line of `report.Rmd`.

- m1. Line 74: an empty `# Abstract` heading follows the YAML abstract
  and renders an empty section.
- m2. Line 8 (Abstract) and throughout: "ADEMP-compliant" is a claim
  about form; drop the adjective and let the structure speak.
- m3. Abstract: three sentences exceed 30 words (31, 33, 36), against
  STYLE.md Section 5.
- m4. Lines 163 to 173 (Senn's model): $\varepsilon_{ij} \sim N(0,
  2\sigma^2)$ and "the factor of 2 arises because $d_{ij}$ is the
  difference of two independent within-period observations" is fine,
  but $\sigma^2$ is then reused at line 149 for "within-patient,
  between-period variance" and $\sigma_w$ at line 519 for the
  residual SD. Three meanings of $\sigma$ in one paper.
- m5. Lines 167 and 508: $\tau^2$ is defined as the
  treatment-by-patient interaction variance (Senn's usage, and
  NOTATION.md's), then $\sigma_\tau$ is used for the between-patient
  baseline SD, which is not treatment heterogeneity. Rename the
  baseline SD $\sigma_u$ (the model already uses $u_i$) and state
  explicitly that $\tau^2 = 0$ in this DGP.
- m6. Line 519: $\rho = 0.3$ is called a "decay rate"; it is the
  lag-one correlation per week, and NOTATION.md reserves $\rho$ for
  exactly that. Say so.
- m7. Lines 541 to 547, Units: NOTATION.md (Units, revised 2026-10-09)
  now assigns weeks to paper 04; the day scale is tolerated only for
  the v1 study. Either convert or state that this is the v1 study.
- m8. Line 218 to 220: "standard sample size software assumes $n - 1$
  degrees of freedom and therefore slightly underestimates power for
  the fixed-effects case" is asserted without a source or a number.
- m9. Line 487: the "primary inferential output" is power at $p <
  0.05$, but the Estimands paragraph calls the coefficient the
  "primary estimand" and power a "secondary estimand" in the v1
  protocol; align the two.
- m10. Lines 604 to 671: the loading chunk is `cache = TRUE`; a cached
  chunk that reads an `.rds` will silently keep stale numbers after a
  rerun of the driver. Drop the cache or key it on the file's mtime.
- m11. Lines 53 to 71: the "fallback stubs" that let the document knit
  with `NA` placeholders are the pattern the remediation plan (Problem
  7) asked to replace with `stopifnot()` guards; they remain.
- m12. Lines 752 to 760 and Figure `effect-distributions`: the
  "distribution of treatment effect estimates" is 500 draws from a
  normal with the cell mean and SD, seeded with `set.seed(1L)`. The
  caption presents it as simulation output. Remove it, or save the
  replicates and plot them.
- m13. Lines 762 to 786 and Figure `sample-trial`: the "representative
  individual N-of-1 trial" is a separate synthetic patient with noise
  SD 1 (not 2.3) and the Results report its estimate as "close to the
  true effect size". It is illustrative at best; label it so, use the
  DGP's SD, and drop the inference.
- m14. Lines 788 to 794 and Figure `p_power_sensitivity`: a one-bar
  bar chart ("Exponential" only) with a dashed 80% line. Remove until
  S3 exists.
- m15. Line 1057 to 1061: a chunk comment says the headline objects
  are "produced by the simulation script sourced above"; nothing is
  sourced.
- m16. Lines 1164 to 1207 and 1379 to 1441: two Morris-style tables
  report the same primary cell; the second adds MSE computed as
  $\mathrm{bias}^2 + \mathrm{SE}^2$ and a "Mean estimate" from the
  synthetic draws of m12. Keep one.
- m17. Line 1315 and 1486: "analysable" is British; use "analyzable".
  `report_short.Rmd` has "Generalisation".
- m18. Line 134: "But how large is the advantage, and how does it
  depend on the carryover regime?" is a rhetorical question
  (STYLE.md, Section 1).
- m19. Fifteen sentences begin "We note", "We see that", "We should
  note", "We think" or "We may"; STYLE.md asks for the fact, not the
  announcement of it.
- m20. No Glossary section (STYLE.md, Section 4). "Hybrid", "OL",
  "BDC", "$D_{bc}$", "mean-moderation architecture" and "Architecture
  A" appear without definition; "Architecture A" is retired notation
  (project CLAUDE.md).
- m21. Line 2018: "Author Contributions: [To be filled in]"; line
  2079: "[GitHub repository URL]"; lines 2075 to 2082 duplicate the
  competing-interest, data-availability and funding statements already
  given as sections.
- m22. The `in_header: claudecode.tex` line is an inert leftover
  (project CLAUDE.md, "Paragraph structure").
- m23. Line 1788: Blackston's effect is written $\tau = 0.25$;
  $\tau$ is this paper's heterogeneity symbol. Write $\delta = 0.25$
  or "their treatment effect of 0.25".
- m24. Line 1802: "our S12 null calibration 0.10-0.13" refers to a
  two-period crossover arm this paper never simulates.
- m25. `report_short.Rmd` versus `report.Rmd`: the README names
  `report.Rmd` as canonical and the short report carries a
  superseded banner, so the submission is clear. But the short report
  still reports RCT power 0.189 against this paper's 0.675 for "the
  same" comparison (the short report used $N = 100$ for the RCT and
  10 participants for the N-of-1 arm, so the two are not the same
  comparison at all). Delete `report_short.Rmd`, or move it to an
  archive directory with its workspace, so no reader meets both.
- m26. Bibliography hygiene. Twelve cite-keys name an author who is
  not an author of the entry (`araujoComparisonFourMethods2014` is
  Chen and Chen; `kratochwillSinglecaseExperimentalDesigns2013` is
  Smith 2012; `kravitzSinglepatientNof1Trials2014` is Duan et al.
  2013; `ledfordSinglecaseDesignAnalysis2018` is Lobo et al. 2017;
  `michielsMultiSCEDToolMetaanalyzing2020` is Declercq et al. 2019;
  `patsopoulosUnderstandingControlledTrials2011` is Sibbald and
  Roberts 1998; `yapBehavioralCarryoverEffect2024` is Shi and Ye;
  `pequignatIdiographicPower2023` is Tueller et al.;
  `chenMethodologicalReviewRandomised2024` is Hawksworth et al.;
  `mitchellDevelopmentSetKey2024` is Chatters et al.;
  `mcdonaldMakingSwitchCase2019` is Margolis and Giuliano;
  `nam2021expanding` is Davidson et al. 2018). Keys are invisible in
  the rendered PDF, but they have already produced one in-text error
  (M10.3) and will produce more. `schork2015` has author "News, NJ
  Schork - Nature and 2015, undefined". `fda2021clinical` is a trade
  news page cited, at line 332, as "recognition of N-of-1 trials by
  regulatory agencies"; cite the FDA guidance itself.
  `hussey2020sced` and `nof1sced2024analysis` have no venue or URL.
  `Khazaie2016Prazosin`, `psu2024crossover` and
  `matthewsSinglesubjectsDesign2023` are not sources a *Statistics in
  Medicine* reader expects for methodological claims. `schmidStatist
  icalDesignAnalytic2014` carries the note "Authors as reported by a
  web search; not checked against the PDF". `senn2002`, `araujo2016`,
  `senn2016mastering`, `sturdevant2021carryover`, `punjaNof1TrialsCan
  2016` and `niklesAggregatingSinglePatient2011` sit in the pool
  uncited although each is more relevant than several cited entries.

## 6 Missing and questionable references

| Claim location | Issue | Suggested source (full citation) |
|---|---|---|
| Intro "Clinical Application", efficacy "in multiple RCTs [@raskind2018]" | Cited trial was null (VERIFIED) | Raskind MA, Peskind ER, Chow B, et al. Trial of prazosin for post-traumatic stress disorder in military veterans. *N Engl J Med* 2018;378:507-517. doi:10.1056/NEJMoa1507598. Cite as the null trial and discuss. |
| Intro "N-of-1 trials ... each patient serves as their own control" | Foundational source absent | Guyatt G, Sackett D, Taylor DW, Chong J, Roberts R, Pugsley S. Determining optimal therapy: randomized trials in individual patients. *N Engl J Med* 1986;314:889-892 (not verified in this session). |
| Intro "aggregation ... requires meta-analytic approaches [@zucker2010]" | Original aggregation paper absent | Zucker DR, Schmid CH, McIntosh MW, D'Agostino RB, Selker HP, Lau J. Combining single patient (N-of-1) trials to estimate population treatment effects and to evaluate individual patient responses to treatment. *J Clin Epidemiol* 1997;50:401-410. doi:10.1016/S0895-4356(96)00429-5 (VERIFIED on PubMed). |
| Intro "Zucker et al. ... their precision formula" | Misattribution (SUSPECTED); formula is Senn's | Senn S. Sample size considerations for n-of-1 trials. *Stat Methods Med Res* 2019;28:372-383. doi:10.1177/0962280217726801. |
| Intro "Sample Size and Power Considerations" | Current closed-form planning work absent | Yang J, Steingrimsson JA, Schmid CH. Sample size calculations for n-of-1 trials. arXiv:2110.08970, 2021. |
| Intro "Sample Size and Power Considerations" | Variance-components view of personalized medicine absent (in bib, uncited) | Araujo A, Julious S, Senn S. Understanding variation in sets of N-of-1 trials. *PLoS ONE* 2016;11:e0167167. Senn S. Mastering variation: variance components and personalised medicine. *Stat Med* 2016;35:966-977. |
| Intro "Prior Simulation Comparisons" | Directly overlapping simulator not cited | Percha B, Baskerville EB, Johnson M, Dudley JT, Zimmerman N. Designing robust N-of-1 studies for precision medicine. *J Med Internet Res* 2019;21:e12641. doi:10.2196/12641 (VERIFIED). |
| Methods, crossover carryover handling | Standard crossover references absent | Jones B, Kenward MG. *Design and Analysis of Cross-Over Trials*, 3rd ed. Chapman and Hall/CRC, 2014. Senn S, Lambrou D. Robust and realistic approaches to carry-over. *Stat Med* 1998;17:2849-2864. Senn S. *Cross-over Trials in Clinical Research*, 2nd ed. Wiley, 2002 (in bib, uncited). |
| Methods "ADEMP" | Fine; also cite the companion on MCSE for paired designs if the paired difference is retained | Morris TP, White IR, Crowther MJ. *Stat Med* 2019;38:2074-2102 (cited). |
| Methods, 3-day half-life "[@Khazaie2016Prazosin]" | Source does not give a half-life (SUSPECTED) | A pharmacokinetic source for prazosin (plasma half-life 2 to 3 h), and a stated rationale for any behavioral carryover scenario. |
| Methods, calibration to "8/week, SD 1.5, SD 2.3" | No source | Calibration table with the published trial statistics used (Raskind 2003, 2013) and the conversion from CAPS item scores. |
| Discussion, Tang and Landes numbers | Wrong rows and table (VERIFIED) | Tang J, Landes RD. Some t-tests for N-of-1 trials with serial correlation. *PLoS ONE* 2020;15:e0228077, Table 2, rows $\rho \in \{0, 0.4, 0.6\}$. |
| Discussion, "Senn and Julious" | Single author (VERIFIED) | Senn S. The analysis of continuous data from n-of-1 trials using paired cycles: a simple tutorial. *Trials* 2024;25:128. doi:10.1186/s13063-024-07964-7. |
| Discussion, "Chen, Chen, and Xia" | Two authors (VERIFIED) | Chen X, Chen P. *PLoS ONE* 2014;9:e87752. doi:10.1371/journal.pone.0087752. |
| Intro "regulatory agencies [@fda2021clinical]" | Trade-press page | The FDA draft guidance itself, or omit. |
| Discussion "Implications for precision medicine" | A real aggregated N-of-1 program that replaced an RCT would anchor the clinical claim | Stunnenberg BC, Raaphorst J, Groenewoud HM, et al. Effect of mexiletine on muscle stiffness in patients with nondystrophic myotonia evaluated using aggregated N-of-1 trials. *JAMA* 2018;320:2344-2353 (not verified in this session). |
| Discussion "Recent design refinements" | A 2026 review with sample-size material | Konigorski S, Ried-Larsen M, Schmid CH. Digital N-of-1 trials and their application in experimental physiology. *Exp Physiol* 2026;111:3525-3545. doi:10.1113/EP092753 (VERIFIED on PubMed). |

## 7 Suggestions for strengthening the paper

1. Execute the rebuild protocol. Its closed-form-first plan (Section
   5.1 there) is the right design for this question, and the
   compendium already has the machinery (`11-strawman-stress-test.R`,
   `58-common.R`, the CR2 analyses).
2. Put the closed form in the paper. For the linear DGP it is a
   half-page: GLS information averaged over sequences, noncentral $t$.
   Readers of *Statistics in Medicine* will want it, and it exposes the
   observation-count effect cleanly: the N-of-1 to RCT SE ratio is
   $0.336/0.778$, and a reader can see what changes when the RCT gets
   eight visits.
3. Make carryover the actual subject. Compare schedule-based
   exclusion, $D_{bc}$ at a fixed analysis half-life, and an explicit
   carryover regressor, under true half-lives that include the
   pharmacological range and a labeled behavioral scenario, with size
   checked under the null for each rule. That comparison is the one
   this paper keeps deferring to "S11", and it is the only one with
   novelty.
4. Add treatment-effect heterogeneity ($\tau^2 > 0$). Without it the
   paper cannot address the Senn and Zucker trade-off its Introduction
   spends three pages on, and the N-of-1 design's distinctive output,
   individual effects, is undefined.
5. Save replicate-level output, advance the RNG stream per cell, and
   compute the paired MCSE from the replicate-level agreement.
6. Replace the synthetic figures (estimate density, sample patient,
   one-bar sensitivity chart, degenerate heatmap) with figures from
   the replicates: power by effect size with MCSE bars; size by
   analysis rule; coverage by rule.
7. Reframe the clinical anchor. Prazosin's evidence base is now mixed
   (2013 positive, 2018 null); a paper that uses it as the motivating
   example should say so and should use CAPS units or state a defended
   conversion.
8. Cut the Introduction's tutorial on closed-form power to the
   equations the paper uses, and move the twelve-paragraph literature
   summary to a short table with one row per study (design, DGP,
   analysis, replicates, main finding).
9. Apply STYLE.md: declarative sentences, no "We note", a Glossary,
   and the compendium's symbols (NOTATION.md lint is clean; the
   collisions in m4 to m6 are semantic, not typographic).

## 8 Scope of this review

**Read in full.** `report.Rmd` (all 2,081 lines including every code
chunk); `report_short.Rmd`; `README.md`; `00-remediation-plan.md`;
`PLAN_power_methods_section.md`; `references.bib` (all 70 entries by
key, title, author and year; full text of 25 entries);
`production-run.R` (all 383 lines); `00-ademp-protocol.md` (lines 1 to
1100 of 2,637: notation, Sections 1 to 8); `00-ademp-protocol-v1.md`;
`docs/58` Sections 1, 2, 12.2, 12.3 and its heading list; STYLE.md;
NOTATION.md (units and main-effect rows); the project CLAUDE.md; the
2026-08-20 pub_review whitepaper's issues 2.6 and 2.9.

**Executed (VERIFIED).** Read the stored `.rds` and printed every cell
of `primary`, `s1_carryover` and `meta`. Ran a scratch script that (a)
showed two consecutive `mclapply` calls after one `set.seed` return
identical draws; (b) computed the closed-form N-of-1 SE with and
without exclusion over the 70 balanced sequences, the mean number of
excluded periods, the ICC, the RCT SE, and noncentral-$t$ power for
both arms at each effect size; (c) computed the carryover multipliers,
the RCT Type I interval and the coverage $z$-scores. Ran the abstract
sentence-length count, the project notation linter (clean), and greps
for em dashes (none), British spellings, hedging phrases, and
"hybrid". Confirmed by `find` that the two cited companion documents
and a `02-deviations-log.md` for this study do not exist.

**Inspected, not executed.** The driver's DGP, exclusion rule,
inference code and RCT generator; the Rmd's inline computations; the
bib entries. Did not run the production simulation (20 to 40 minutes).
Did not read `simulation.R` (the retired DerSimonian-Laird driver) or
the `nof1power` package source, since `production-run.R` imports
neither; the README's statement that the paper depends on `nof1power`
is therefore inaccurate for this driver (INSPECTED).

**Not checked.** Whether `sim_workspace.RData` reproduces the numbers
in `report_short.Rmd` (its generating script is absent, as the banner
says). The content of `docs/58` Sections 3 to 11 beyond what the
protocol quotes. The Blackston 150/100 sample sizes (lost in text
extraction). The Tang and Landes Type I rates read from their figures.
Whether Khazaie 2016 contains any half-life statement (full text not
retrieved). Guyatt 1986 and Stunnenberg 2018 are cited from memory
and marked as such.

**Literature searches run.** PubMed (via the PubMed MCP tool, no email
or identifier sent): metadata for PMIDs 38365817, 28882093, 29414272,
24503561, 40071868, 30932871, 39744919, 40690192, 9179098, 20863658,
32017772, 31698799, 31269712, 34713026, 42318774, 38563531; searches
for Carrozzo 2025, Selukar 2025, Jiang 2025, Raskind 2018, Percha
2019, "sample size n-of-1 trials" with Schmid, Zucker 1997, Guyatt
1986 (no useful hit), Shi and Ye 2024, and "N-of-1 trials AND (power
OR sample size) AND simulation AND (parallel group OR crossover)" from
2022 (7 hits, all inspected; two were biosimilar complete-N-of-1
designs, not relevant). PMC full text retrieved for Blackston et al.
2019 (PMC6955665) and Tang and Landes 2020 (PMC6999905). Web search
for the Yang, Steingrimsson and Schmid preprint (arXiv:2110.08970
confirmed; no journal version found). Zotero MCP was unavailable
(connection refused), so the project's `nof1` collection was not
consulted.
