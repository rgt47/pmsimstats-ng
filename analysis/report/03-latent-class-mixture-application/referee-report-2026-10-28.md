# Referee report: Latent-class and mixture-model formulations for biomarker-treatment interaction in N-of-1 trials
*2026-10-28 07:32 PDT*

Manuscript: `analysis/report/03-latent-class-mixture-application/report.Rmd`
(with `supplement.Rmd`, `references.bib`, `simulation-plan.md`,
`README.md`). Target venue assumed: *Statistics in Medicine* (methods
article). Reviewer role: statistical referee. Evidence labels used
throughout: VERIFIED (re-derived or executed), INSPECTED (read and
judged), SUSPECTED (plausible but not confirmed).

## 1 Summary of the manuscript

The manuscript argues that biomarker moderation of treatment response
in aggregated N-of-1 trials can be encoded either as a continuous
interaction (the convention in the trial literature and in
Hendrickson et al. 2020) or as a discrete latent-class structure (the
convention in psychometrics and marketing science), and that the
choice is a scientific question rather than a parameterization
convenience. In its own words the paper has four aims: (1) formalize
the latent-class and finite-mixture formulations with identifiability
conditions specific to the N-of-1 sample-size regime; (2) characterize
estimation strategies and diagnostics for distinguishing genuine from
spurious class structure at small $N$; (3) report a factorial
simulation benchmarking mixture-DGP and continuous-DGP power under
linear-mixed and class-aware analyses with carryover; (4) re-analyze
the prazosin-PTSD data (deferred).

The delivered content is: a literature narrative (main text plus
supplement Section 1); a two-component mixture DGP with logistic
gating in the standardized biomarker; a second-moment formula for the
marginal correlation induced by a two-class mixture; three prose
"failure regimes" of the MVN approximation; a four-row taxonomy table
and a five-row "moderation spectrum" table; a software list; and a
240-replicate-per-cell pilot at $N = 70$ with gating slope
$\psi_g \in \{0, 1.0, 1.5\}$ comparing a linear-mixed test, an
unconditional Wald test on the `lcmm::hlme` gating coefficient, and a
BIC-conditional version of that test. The headline results are an
unconditional-Wald Type I error of 0.197, a BIC-conditional Type I
error of 0.004 with power 0.025 and 0.004, and linear-mixed power of
0.325 and 0.542 with Type I error 0 of 240.

## 2 Overall assessment and recommendation

**Recommendation: reject in the present form; a resubmission after
the production study is run would be considered as a new submission.**

Justification. The manuscript is explicit that the simulation program
on which its contribution rests has not been executed, and the pilot
that stands in for it does not test what the paper says it tests.
Specifically: (i) on the all-on-drug open-label design used for every
pilot cell, the package's `lme_analysis()` does not fit the
`bm:Dbc` model at all; it fits `Sx ~ bm + t + bm:t`, so the
"linear-mixed `bm:Dbc` Wald test" reported throughout is a
biomarker-by-linear-time test, and it is compared against an `hlme`
fit with a different fixed-effects structure (Major comment 1,
VERIFIED by execution); (ii) that test's Type I error of 0 of 240 is
not "nominal", it is a severely conservative test (model SE exceeds
the empirical SD by a factor of 1.67; 1 rejection in 2000 null
replicates in an unreported production run), so its "power" is not
comparable to any other arm (Major comment 2, VERIFIED); (iii) the
inflation of the unconditional gating Wald test is a textbook
consequence of a parameter that is unidentified under the null
(Davies 1977, 1987; Hansen 1996), not of the Bauer-Curran
non-normality mechanism the paper invokes, and it would not go away
at any sample size (Major comment 4, INSPECTED); (iv) the DGP the
paper calls "heterogeneous random slopes" has no random slope; it
zeros the drug-response component for non-responders (Major comment
5, VERIFIED); (v) a 4000-replicate production run of the Type I cells
exists on disk and contradicts the manuscript's statement that the
production driver has not been run (Major comment 3, VERIFIED).

Aims 1 and 2 are not delivered: no identifiability condition is
derived and no diagnostic is developed beyond listing posterior
entropy and BIC. Aim 3 is not delivered. The remaining material is a
well-read position piece. A position piece with an unexecuted
simulation plan does not meet the bar for a methods article in
*Statistics in Medicine* or *Trials*. The literature narrative is
competent and parts of it would survive into a resubmission, but the
closest competing work (growth mixture modeling of treatment effects
in randomized trials, Muthen et al. 2002; the variance-components
treatment of individual response in crossover and N-of-1 designs,
Senn 2004, 2016) is not cited, and it weakens the novelty claim that
"the biostatistical literature on N-of-1 trials has not engaged"
mixture heterogeneity.

## 3 Significance and novelty

**What is new, precisely.** (a) The explicit statement that the
covariance-only (correlation) encoding of Hendrickson et al. (2020)
is one point on a moment-by-moment spectrum of moderation encodings,
with a mixture DGP at the discrete end. (b) The closed-form product
formula for the squared marginal correlation induced by a two-class
mixture with within-class independence (Section "MVN approximation";
VERIFIED, see Section 8). (c) Empirical evidence, at one sample size
and one design, that an unconditional Wald test on an `hlme` gating
coefficient is anti-conservative and that BIC does not select two
classes under a deterministic two-class DGP at $N = 70$.

**What is not new.** Item (a) is a framing contribution. Item (b) is
the standard between/within decomposition of covariance under a
latent grouping; it is an exercise, not a result. Item (c) is a known
consequence of non-regular testing in mixtures and has been
documented in the class-enumeration literature at sample sizes
including the N-of-1 range.

**Competing and foundational literature the paper does not engage.**

- Muthen B, Brown CH, Masyn K, Jo B, Khoo ST, Yang CC, Wang CP,
  Kellam SG, Carlin JB, Liao J. General growth mixture modeling for
  randomized preventive interventions. *Biostatistics* 2002;
  3(4): 459-475. doi:10.1093/biostatistics/3.4.459. This is the
  direct precedent: a finite-mixture random-effects model used to
  estimate treatment effects within latent trajectory classes in a
  randomized longitudinal trial. The manuscript's "latent-class
  encoding of treatment-effect heterogeneity" is this model with the
  biomarker as a class-membership covariate.
- Senn S. Individual response to treatment: is it a valid assumption?
  *BMJ* 2004; 329(7472): 966-968. doi:10.1136/bmj.329.7472.966; and
  Senn S. Mastering variation: variance components and personalised
  medicine. *Statistics in Medicine* 2016; 35(7): 966-977.
  doi:10.1002/sim.6739 (PMID 26415869, located in this review's
  PubMed search). These papers are the N-of-1 and crossover
  literature's treatment of exactly the responder/non-responder
  dichotomy the manuscript says that literature lacks: they show that
  the patient-by-treatment interaction variance is identifiable only
  with repeated within-patient crossover, and that apparent
  responder classes are usually artifacts of within-patient
  variation. Araujo A, Julious S, Senn S. Understanding variation in
  sets of N-of-1 trials. *PLoS ONE* 2016; 11(12): e0167167, extends
  this with simulation.
- Lubke GH, Neale MC. Distinguishing between latent classes and
  continuous factors: resolution by maximum likelihood? *Multivariate
  Behavioral Research* 2006; 41(4): 499-532.
  doi:10.1207/s15327906mbr4104_4. A simulation study of precisely the
  continuous-versus-discrete question the manuscript poses, with
  class separation and within-class sample size as the governing
  factors. Bauer DJ, Curran PJ. The integration of continuous and
  discrete latent variable models: potential problems and promising
  opportunities. *Psychological Methods* 2004; 9(1): 3-29, is the
  companion position piece.
- Davies RB. Hypothesis testing when a nuisance parameter is present
  only under the alternatives. *Biometrika* 1977; 64(2): 247-254 and
  1987; 74(1): 33-43; Hansen BE. Inference when a nuisance parameter
  is not identified under the null hypothesis. *Econometrica* 1996;
  64(2): 413-430. The correct explanation of the gating-Wald
  inflation (Major comment 4).
- McLachlan GJ. On bootstrapping the likelihood ratio test statistic
  for the number of components in a normal mixture. *Applied
  Statistics* 1987; 36(3): 318-324; Lo Y, Mendell NR, Rubin DB.
  Testing the number of components in a normal mixture. *Biometrika*
  2001; 88(3): 767-778; Chen H, Chen J, Kalbfleisch JD. A modified
  likelihood ratio test for homogeneity in finite mixture models.
  *JRSS B* 2001; 63(1): 19-29. The calibrated tests of one versus two
  classes; the manuscript cites only Nylund et al. (2007) for the
  bootstrap LRT.
- Dziak JJ, Lanza ST, Tan X. Effect size, statistical power, and
  sample size requirements for the bootstrap likelihood ratio test in
  latent class analysis. *Structural Equation Modeling* 2014; 21(4):
  534-552. Directly addresses power to detect a second class as a
  function of $N$ and class separation, which is the manuscript's
  open question.
- Hipp JR, Bauer DJ. Local solutions in the estimation of growth
  mixture models. *Psychological Methods* 2006; 11(1): 36-53.
  Relevant to the `gridsearch(rep = 5, maxiter = 30)` choice and to
  the degenerate gating estimates observed (Major comment 4).
- Lubke GH, Muthen BO. Performance of factor mixture models as a
  function of model size, covariate effects, and class-specific
  parameters. *Structural Equation Modeling* 2007; 14(1): 26-47. The
  Monte Carlo study the manuscript attributes to Lubke and Muthen
  (2005).
- Vermunt JK. Latent class modeling with covariates: two improved
  three-step approaches. *Political Analysis* 2010; 18(4): 450-469;
  Asparouhov T, Muthen B. Auxiliary variables in mixture modeling:
  three-step approaches using Mplus. *Structural Equation Modeling*
  2014; 21(3): 329-341. The standard treatment of covariates (here,
  the biomarker) in class membership, which bears on whether the
  one-step `classmb = ~bm` specification is the right class-aware
  test.
- Proust-Lima C, Sene M, Taylor JMG, Jacqmin-Gadda H. Joint latent
  class models for longitudinal and time-to-event data: a review.
  *Statistical Methods in Medical Research* 2014; 23(1): 74-90.
- Percha B, Baskerville EB, Johnson M, Dudley JT, Zimmerman N.
  Designing robust N-of-1 studies for precision medicine: simulation
  study and design recommendations. *Journal of Medical Internet
  Research* 2019; 21(4): e12641. doi:10.2196/12641. Post-2015 N-of-1
  simulation work with carryover and wash-in, which the manuscript's
  review of the N-of-1 methods literature omits.
- Nylund-Gibson K, Choi AY. Ten frequently asked questions about
  latent class analysis. *Translational Issues in Psychological
  Science* 2018; 4(4): 440-461; van de Schoot R, Sijbrandij M,
  Winter SD, Depaoli S, Vermunt JK. The GRoLTS-checklist: guidelines
  for reporting on latent trajectory studies. *Structural Equation
  Modeling* 2017; 24(3): 451-467. Post-2015 reporting standards for
  mixture studies that a resubmission should follow.

**Does it clear a top-tier bar?** Not in this form. It would if the
production program (or a disciplined subset: Study 5 at three sample
sizes with calibrated tests, plus Study 1 or 3 at one design with an
off-drug contrast) were run, the analysis arms were matched on
fixed-effects structure, the tests were calibrated, and the effect
sizes were placed on a common scale.

**Contemporariness.** The bibliography stops, with two exceptions
(Kent 2020, Morris 2019), at 2017. The three-step, BLRT-power,
GRoLTS, and N-of-1 simulation literatures since 2014 are absent.

## 4 Major comments

Ordered with correctness first.

**M1. The pilot's linear-mixed arm does not fit the model the paper
says it fits, and the two arms fit different models (VERIFIED).**
Location: Abstract (Methods, Results), Section "Pilot results",
Section "Reproducibility"; drivers `quick-sim-driver-parallel.R`,
`03-study5-pilot.R`, `04-pilot-power-random-slopes.R`. Every pilot
cell uses `buildtrialdesign(... ondrug = list(pathA = rep(1, 8)))`,
an open-label design with all eight post-baseline visits on drug. In
`R/lme_analysis.R` lines 130-147, when no participant has
within-person variation in drug state (`varInDb == FALSE`), the
fixed-effects formula becomes `Sx ~ bm + t + bm:t` and the reported
coefficient is `bm:t`. I ran one null replicate with the pilot's
exact design, parameters and seed scheme and obtained
`FORMULA FIT BY lme_analysis: Sx ~ bm + t + bm:t` with the `tTable`
row `bm:t` returned as `beta`. The `Dbc` column is 0 at baseline and
1 at all eight visits. Meanwhile `lcmm_analysis()` hard-codes
`fixed = Sx ~ bm + t + Dbc + bm:Dbc` with `mixture = ~Dbc`, so in the
class-aware arm `Dbc` is a post-baseline indicator and `bm:Dbc` is a
biomarker-by-post-baseline term. The paper's "lme `bm:Dbc` Wald"
versus "lcmm gating Wald" comparison is therefore between a
biomarker-by-linear-time test in a model without `Dbc` and a gating
test in a model with a post-baseline step. Why it matters: the paper's
estimand is the biomarker-by-treatment interaction $\beta_{bm:D}$,
and NOTATION.md Part 3 states that the OL design "supplies no
within-subject interaction information". Nothing in the pilot speaks
to that estimand, and the carryover discussion that frames the pilot
(Sections "Implications for power under carryover", "Discussion") is
inapplicable to a design with no discontinuation (the
`carryover_t1half = 1.0` setting in `op` has no effect because
`tsd = 0` everywhere). Remedy: rerun on a design with an off-drug
contrast (OL+BDC or Hybrid, as the protocol specifies), state the
fitted formula of each arm in the manuscript, and make the
fixed-effects structure identical across arms so that the only
difference is the mixture.

**M2. The linear-mixed test is not "nominal"; it is severely
conservative, which invalidates the power comparison (VERIFIED).**
Location: Abstract ("holds nominal Type I"), Section "Pilot results"
("holds nominal Type I (0 of 240 under the null)"). Under a
calibrated 0.05 test, 0 rejections in 240 has probability
$0.95^{240} = 4.5 \times 10^{-6}$. From the stored replicates
(`03-latent-replicates.rds`, null cell): empirical SD of
$\hat\beta$ is 0.049, mean model SE is 0.082, ratio 1.67; the
smallest of 240 null p-values is 0.090; the median is 0.69. The
100-replicate pre-flight shows the same (0 of 100; ratio 1.65). The
unreported production run (M3) shows 0.0025 at $N = 35$ and 0.0005 at
$N = 70$ over 2000 replicates each. Why it matters: a test whose size
is about 0.001 cannot have its power compared with tests of size 0.20
or 0.004; the "lme dominates" conclusion is a comparison of three
mis-calibrated tests. The conservativeness is itself a finding about
the `bm:t` specification with `corCAR1` on a Gompertz trajectory and
deserves diagnosis, but it must not be described as nominal. Remedy:
report $\kappa$ (model SE over empirical SD, as paper 10 does), state
the realized size with MCSE, and either calibrate (parametric
bootstrap or size-adjusted critical values) or compare power at
matched realized size.

**M3. A production Type I error run exists on disk and contradicts
the manuscript (VERIFIED).** Location: Abstract ("which we have not
yet executed"), Section "Reproducibility" ("`10-study5-production.R`
(committed, not yet run)"), Discussion. The file
`analysis/data/study5/study5-progress.rds` (2026-05-09) holds 4000
replicates: cell 1 ($N = 35$) and cell 2 ($N = 70$), both null,
2000 each. Computed from it: `lcmm` convergence 0.942 and 0.947;
unconditional gating Wald rejection 0.166 (MCSE 0.009) and 0.240
(MCSE 0.010); BIC selects two classes in 2.0 percent and 0.9 percent;
BIC-conditional rejection 0.006 and 0.006; within-class `bm:Dbc`
rejection 0.000 and 0.000; linear-mixed rejection 0.0025 and 0.0005.
Why it matters: these are the paper's primary Type I estimands at
eight times the pilot's precision and at a second sample size, and
they are more informative than the 223-replicate figure quoted. The
statement that the run has not happened is false as written (it may
be incomplete; the file name says "progress"). Remedy: report these
cells, say what the run's status is, and reconcile the deviations log
(`02-deviations-log.md` D1 to D3) with the manuscript.

**M4. The explanation of the unconditional-Wald inflation is wrong,
and the remedy proposed follows from the wrong explanation
(INSPECTED; the diagnostic numbers are VERIFIED).** Location:
Abstract ("consistent with the Bauer-Curran and Lubke-Muthen
overextraction predictions"), Section "Pilot results", Supplement
Section 6 ("finite-sample identifiability artifact predicted by the
Bauer-Curran ... literature"). Bauer and Curran (2003) concern
spurious classes extracted from a single non-normal population. The
null DGP here is a single multivariate normal, so the pilot does not
test, and cannot be "consistent with", that mechanism. What the pilot
shows is the non-regular testing problem: under one class, the gating
slope $\psi_g$ is a parameter that exists only under the alternative,
the information matrix is singular, and the Wald statistic has no
$\chi^2_1$ limit at any $N$ (Davies 1977, 1987; Hansen 1996). The
stored estimates bear this out: across the 223 converged null fits
the median $|\hat\psi_g|$ is 1.4, 18 percent exceed 10 in absolute
value, the SD is 70, and the z-statistic SD is 1.44; these are
boundary solutions with a near-empty class, not an SE that "does not
absorb labeling uncertainty". The same pattern persists at gating
slopes 1.0 and 1.5 (z SD 1.67 and 1.77), so it is not confined to the
null. Why it matters: the paper's proposed fix (gate the Wald test on
BIC) treats a symptom; the principled tests of one versus two classes
are the bootstrap LRT (McLachlan 1987; Nylund 2007), the
Lo-Mendell-Rubin adjusted LRT (2001), or the modified LRT of Chen,
Chen and Kalbfleisch (2001), and the test of whether the biomarker
predicts class membership given two classes is a separate, regular
problem. Also, `lrt_df = 4` in the wrapper is the raw parameter-count
difference, and the naive $\chi^2_4$ comparison reported in the
supplement (0.146) is known to be wrong for the reason the supplement
itself gives. Remedy: reframe the inflation as a nuisance-parameter
non-identification problem with the citations above; replace the
"BIC-conditional" test with a BLRT-gated or LMR-gated procedure and
report its size; retain BIC selection frequency as a separate
estimand.

**M5. The DGP is mislabeled: there is no random slope in the
"heterogeneous-random-slopes DGP" (VERIFIED).** Location: Abstract,
Table `tab:spectrum` row 4, Section "Pilot results", Supplement
Section 6. `02-heterogeneous-slopes-dgp.R` multiplies the `BR`
component at on-drug visits by 1.0 (class 1) or 0.0 (class 0). The
non-responder class therefore has its drug-response component, mean
and random variation alike, deleted; the responder class is the
unmodified MVN draw. This is a deterministic two-class scaling, not
the Verbeke-Lesaffre mixture over random slopes that the protocol
(`00-ademp-protocol.md`, Study 3) specifies and that
`07-true-heterogeneous-slopes-dgp.R` implements with
$\beta_R \sim \mathcal N(\mu_R, \sigma_R^2)$. The header of file 07
states this discrepancy explicitly. Two consequences follow. First,
the class-conditional residual variance differs between classes (the
`BR` variance, $\sigma_{BR}^2 = 25$ at each visit, is absent in class
0), so the pilot DGP is the paper's own "class-varying covariance"
failure regime, while the `hlme` analysis assumes a common residual
variance; the paper does not notice this. Second, the gating slopes
1.0 and 1.5 are never placed on the common effect-size scale the
protocol mandates; the stored `03-calibration-table.rds` is cited in
the Reproducibility section but is not used by the pilot and does not
contain this DGP (its `marginal_mixture` rows come from
`05-marginal-mixture-dgp.R`). Why it matters: the reader cannot tell
what effect size the 0.542 power corresponds to, and the "Study 5"
label is wrong (the protocol's Study 5 has no gating-slope axis; the
gating-slope cells are Study 3's). Remedy: name the DGP accurately
(deterministic class-specific drug-response scaling), state the
class-conditional variances, report the marginal on-drug
$\mathrm{Cor}(B, BR)$ or Cohen's $f^2$ for each cell from the
calibration machinery, and relabel the pilot by protocol study.

**M6. The Reproducibility section misstates the provenance of the
pilot and nothing in the manuscript is computed from the stored data
(VERIFIED).** Location: Section "Reproducibility". (a) The
240-replicate run was produced by `quick-sim-driver-parallel.R`
(it alone matches the run log: `mclapply mc.cores=4`, 558.6 s,
1440 rows), not by `03-study5-pilot.R`, which is the 100-replicate
sequential null-only pre-flight writing to `output/03-study5-pilot.rds`.
(b) The seed scheme is `SEED_BASE + rep_idx + 1000 * cell_index` with
`mc.set.seed = FALSE`, not "`set.seed(20260507)` at the top of the
replicate loop, with per-replicate seeds derived by `+ rep_idx` under
`parallel::mc.set.seed`". (c) Four drivers (`quick-sim-driver.R`,
`-10min.R`, `-30min.R`, `-parallel.R`) write to the same
`03-latent-replicates.rds`, `03-latent-summary.txt` and
`03-latent-runlog.txt` paths; the stored file's identity depends on
which ran last. (d) `report.Rmd` has one chunk (`setup`) and reads
no data; every number in the Results is typed by hand. (e) The
manuscript contains no table and no figure of the pilot; the driver's
figure (`analysis/figures/quick-sim/03-latent-power.pdf`) is not
included. (f) `fig.path = '../../figures/04-'` belongs to paper 04.
Remedy: a results chunk that reads the summary file and renders a
table with MCSE and the converged-replicate denominators; one driver
per stored artifact; correct the seed description.

**M7. Misattributed and superseded quantitative claims about
carryover power loss (INSPECTED; the supersession is VERIFIED against
the compendium's NOTATION.md Part 2).** Location: Section
"Implications for power under carryover" ("40--60% across designs in
@hendrickson2020"); Introduction lines 133-137 ("Recent methodological
work has documented substantial power losses ... with the loss varying
sharply according to how the moderation effect is represented in the
data-generating process [@hendrickson2020]"); Section "Why this
matters", third point (cites @hendrickson2020 for differential power
profiles and then "the companion carryover-sensitivity work" without a
citation). Hendrickson et al. (2020) implement only the correlation
encoding and do not compare encodings; the comparison and the
"40-60 percent" figure are this compendium's own earlier results, and
NOTATION.md records that the matched-$N$ rerun replaced them with 15.4
percent (correlation encoding) against 6.9 percent (mean-shift
encoding) in the Hybrid design. Remedy: cite the compendium paper or
NOTATION.md Part 2 for the current figures, attribute Hendrickson et
al. only for what they report, and remove "40--60%".

**M8. Over-interpretation of differences inside Monte Carlo error,
and the BIC-conditional power is reported as a test result when it is
a selection-rate result (VERIFIED from the stored data).** Location:
Section "Pilot results", paragraph on non-monotonicity; Abstract.
BIC-conditional rejection is 0.0253 (6 of 237) at $\psi_g = 1.0$ and
0.0042 (1 of 236) at 1.5; the difference is 0.021 with standard error
0.011 ($z = 1.9$). The paragraph that explains "the non-monotonicity"
by the single-class model "capturing the signal" is an explanation
of noise. The quantity that carries the information is the BIC
selection rate itself: two classes are selected in 0.4, 5.1 and 1.7
percent of fits at the three gating slopes, and $\Delta_{BIC} > 6$
never occurs. The conditional test's power is bounded by that rate and
the paper should say so directly. Likewise "sits below the lme
baseline" at $\psi_g = 1.0$ compares 0.304 with 0.325 (difference
within one MCSE). Remedy: report selection frequencies as the
estimand, give differences with MCSE, and delete causal narratives
for differences under two MCSE.

**M9. A fourth test form was run and is omitted from the manuscript
(VERIFIED).** Location: Abstract ("We record three test forms"),
Section "Pilot results". The stored replicates contain the
within-class `bm:Dbc` Wald test from `hlme` (`within_p`). From the
converged fits it rejects 0.000, 0.430 and 0.521 at the three gating
slopes; at $\psi_g = 1.5$ this is within one MCSE of the linear-mixed
arm (0.542) with Type I error 0 of 223, and in the production run
0 of about 1890 at each $N$. The supplement lists it in the pre-flight
table and the wrapper documents it as the "secondary signal", but the
main text drops it. Why it matters: it is the class-aware test with
the best observed operating characteristics in the pilot, and its
omission changes the conclusion that "neither class-aware test is
competitive". Its Type I error of exactly zero also needs the same
calibration diagnosis as M2. Remedy: report all recorded test forms
or justify the exclusion in the protocol deviations log.

**M10. Aims 1 and 2 are stated as contributions but are not delivered
(INSPECTED).** Location: Abstract, Section "Aims and structure",
Sections "Identifiability at trial-relevant sample sizes" and
"Finite-mixture taxonomy". No identifiability condition is derived;
the identifiability section is a paragraph asserting that posterior
entropy "requires $N$ several hundred to several thousand" with three
citations and no analysis. Entropy is a classification-certainty
index, not an identifiability diagnostic, and the sentence "the EM
algorithm for finite mixtures of mixed-effects models requires either
informative priors on class-membership probabilities or strong
separation" is unsupported (`hlme` uses a Marquardt algorithm, not
EM, and no priors enter). The "five-component spectrum" the abstract
says the paper "develops" is a five-row table with one sentence per
row. Remedy: either derive something (for example, the Fisher
information for $\psi_g$ as a function of class separation and $N$
under the paper's own DGP, or the noncentrality of the BLRT via Dziak
et al. 2014) or rewrite the abstract to claim a framing plus a pilot.

**M11. Design and sample-size claims are internally inconsistent with
the pilot (INSPECTED).** Location: Abstract ("trial-relevant sample
sizes ($N = 70$, eight measurement timepoints)"), Section "Why this
matters" ("$N \in [30, 150]$"). The pilot uses nine rows per
participant (baseline plus eight), one randomization path, no
blinding, no discontinuation and constant expectancy, so the `De`
term is dropped by `lme_analysis()`. This is not the Hybrid design of
Hendrickson et al. nor any design in the protocol's design axis
except OL, which the compendium's own glossary says carries no
interaction information. Remedy: state the design and its limitation
in the abstract, or rerun on a design with an off-drug contrast (M1).

## 5 Minor comments

1. Notation, Section "Notation" and "Faithful latent-class DGP":
   $\pi$ is used for the mixing proportion and $\pi(B_i)$ for the
   class probability; NOTATION.md reserves $\pi$ for power. Use
   $\pi_z$ or $\omega$ for mixing and state the departure.
2. Section "Notation": the definition of $D_{bc,it}$ omits the
   pre-exposure value (0 at off-drug visits before first exposure),
   which NOTATION.md requires (2026-10-09 revision). Irrelevant in the
   pilot's design but wrong as a definition.
3. Table `tab:spectrum`, last row: "Architecture A" is a retired
   label (NOTATION.md Part 3). Write "the mean-shift encoding".
4. `tools/notation-lint.pl` reports `bare-E-expect` five times
   (report lines 590, 593, 600; supplement 386, 424): write
   $\mathbb E[\cdot]$.
5. STYLE.md requires a Glossary section before the Abbreviations;
   the manuscript has neither. Terms needing a glossary entry:
   gating slope, `gs`, `ng`, BIC-conditional procedure, OL,
   `bm:Dbc`, `Dbc`, `De`, entropy (as computed in the wrapper).
6. STYLE.md voice: the Abstract has sentences over 30 words (the
   "We commit to ..." sentence is about 45 words; the "We treat the
   choice ..." sentence is about 30; the Results sentence beginning
   "In a 240-replicate-per-cell pilot" is about 35). Split them.
7. STYLE.md terminology: "coupling" appears in the supplement (line
   595); "covariance-only", "covariance moderation" and "mean
   moderation" are on the Avoid list; use "the correlation encoding"
   and "the mean-shift encoding" with one definition each.
8. US English: "Behaviour" (table header), "catalogue" (Abstract),
   "systematised", "finalise", "Parallelisation", "Generalised"
   (supplement), and "labelling", "favour", "parameterisation" in
   code comments that are quoted. Spearman's title is a verbatim
   quotation and may keep its form.
9. Section "MVN approximation": the formula is correct (re-derived:
   $\mathrm{Cov}(B, BR) = \pi(1-\pi)\Delta_B\Delta_{BR}$ under
   within-class independence), but the sentence "bounded by the
   geometric mean of the two intra-class correlations" should read
   "equals the geometric mean of the two between-class variance
   fractions"; these are $\eta^2$ quantities, not intraclass
   correlations, and the relation is an equality. Also state that the
   formula conditions on $Z$ and so requires $B \perp BR \mid Z$,
   which is a restriction beyond the Section 3 DGP (where $B$ enters
   only through the gating).
10. Section "Failure regimes", first regime: "Power under the true
    mixture DGP will then exceed the MVN-LMM prediction" is asserted
    without evidence and the pilot shows the opposite at $N = 70$.
    Hedge or support.
11. Section "Related work" and Supplement 1.4: Zhang and Davidian
    (2001) use a seminonparametric (Gallant-Nychka) random-effects
    density, not a "semi-parametric mixture specification"
    (SUSPECTED mischaracterization; I did not reread the paper).
12. Supplement 1.1: Bartholomew's (1987) contribution is the fourfold
    classification by continuous/categorical manifest and latent
    variables; "graded versus unfolding" is an item-response
    distinction and is not the continuous-versus-discrete question
    (SUSPECTED misattribution).
13. Supplement 1.3: "Lubke and Muthen (2005) conducted the first
    systematic Monte Carlo evaluation of factor mixture models". The
    2005 paper is primarily expository; the Monte Carlo evaluation is
    Lubke and Muthen (2007) (SUSPECTED).
14. Section "Pilot results": state that rejection rates use the
    converged-replicate denominator (223, 237, 236 of 240) and
    report non-convergence as an outcome; under the null, 7 percent
    non-convergence is itself informative, and treating
    non-convergence as non-rejection gives 0.183.
15. Section "Progress to date" quotes per-replicate `lcmm` timing of
    1.57 s from the 100-replicate run; the 240-replicate run averaged
    2.9 s per replicate (four workers sharing BLAS threads). Say which
    run the budget extrapolation uses.
16. Section "Progress to date": "0.250 against 0.05" is from 96
    converged pre-flight fits (MCSE 0.044); say so.
17. The within-class `bm:Dbc` test and the BIC-conditional procedure
    have Type I error exactly zero in every run; a test with zero
    size is not "controlled", it is uninformative about calibration.
18. The supplement's smoke-test numbers (46 percent responders,
    $\mathrm{Cor}(B, BR) = 0.63$, $p = 0.002$, entropy 0.86, BIC 170
    lower) come from a single draw in `99-smoke-test.R` and are not
    stored; they should not be quoted as findings (not re-run here).
19. `references.bib`: key `guyatt1988` is a 1990 paper (Ann Intern
    Med 112:293-299); `kent2020path` and `clark2013` use "and
    others"; most entries lack DOIs; `mplus` year "1998-2017";
    "Raftery 1995" is cited in the supplement and wrapper with no
    bib entry; `vermuntmagidson2002` is unused.
20. README.md is stale: it says the driver directory "is not yet
    created" and lists four scope tasks that no longer match the
    paper. `in_header: ../sim-preamble.tex` plus the inert
    `claudecode.tex` leftover should be checked for the line-number
    package interaction with `longtable`.
21. The manuscript's `date` field is `Sys.Date()`; a submission
    should carry a fixed date.
22. The pilot wrapper fits `hlme(ng = 2)` with `maxiter = 30` and
    `gridsearch(rep = 5)`; report the convergence criterion and the
    fraction of fits at the iteration cap, since Hipp and Bauer
    (2006) show local solutions are common at this budget.

## 6 Missing and questionable references

| Claim location | Issue | Suggested source (full citation) |
|---|---|---|
| Abstract, Background: "the biostatistical literature on N-of-1 trials has not engaged" mixture heterogeneity | Overstated; the responder/variance-components literature exists | Senn S. Individual response to treatment: is it a valid assumption? BMJ 2004;329(7472):966-968. Senn S. Mastering variation: variance components and personalised medicine. Stat Med 2016;35(7):966-977. Araujo A, Julious S, Senn S. Understanding variation in sets of N-of-1 trials. PLoS ONE 2016;11(12):e0167167 |
| Introduction, "Related work": GMM lineage | Omits the GMM-for-randomized-trials precedent | Muthen B et al. General growth mixture modeling for randomized preventive interventions. Biostatistics 2002;3(4):459-475 |
| Introduction, "Two formal answers": continuous versus discrete | Omits the direct simulation treatments | Lubke GH, Neale MC. Distinguishing between latent classes and continuous factors: resolution by maximum likelihood? Multivariate Behav Res 2006;41(4):499-532. Bauer DJ, Curran PJ. The integration of continuous and discrete latent variable models. Psychol Methods 2004;9(1):3-29 |
| Abstract and Pilot results: inflation "consistent with Bauer-Curran" | Wrong mechanism; needs non-regular testing references | Davies RB. Biometrika 1977;64(2):247-254; Davies RB. Biometrika 1987;74(1):33-43; Hansen BE. Econometrica 1996;64(2):413-430 |
| Supplement 6: naive LRT, chi-bar-squared | Cites only McLachlan and Peel (2000) ch. 6 and Nylund (2007) | McLachlan GJ. Appl Stat 1987;36(3):318-324. Lo Y, Mendell NR, Rubin DB. Biometrika 2001;88(3):767-778. Chen H, Chen J, Kalbfleisch JD. JRSS B 2001;63(1):19-29 |
| Discussion, open question 1 (usability at small $N$) | Power of class enumeration at small $N$ is studied | Dziak JJ, Lanza ST, Tan X. Struct Equ Modeling 2014;21(4):534-552 |
| Supplement 6: gridsearch, degenerate gating estimates | Local solutions literature | Hipp JR, Bauer DJ. Psychol Methods 2006;11(1):36-53 |
| Supplement 1.3: "first systematic Monte Carlo evaluation" of FMM | Wrong paper | Lubke GH, Muthen BO. Struct Equ Modeling 2007;14(1):26-47 |
| Wrapper and Supplement 6: `classmb = ~bm` one-step | Covariate-in-membership estimation has a known bias literature | Vermunt JK. Polit Anal 2010;18(4):450-469. Asparouhov T, Muthen B. Struct Equ Modeling 2014;21(3):329-341 |
| Section "Implications for power under carryover": "40--60% ... in @hendrickson2020" | Misattributed and superseded | Cite the compendium's paper 01 or NOTATION.md Part 2 (15.4 vs 6.9 percent at matched $N = 70$) |
| Introduction: N-of-1 methods literature since 2015 | Missing simulation and reporting work | Percha B et al. J Med Internet Res 2019;21(4):e12641. Shamseer L et al. CONSORT extension for N-of-1 trials (CENT) 2015. BMJ 2015;350:h1793 |
| Supplement 6: "Raftery 1995" BIC thresholds | No bibliography entry | Raftery AE. Bayesian model selection in social research. Sociol Methodol 1995;25:111-163 |
| Supplement 1.4: Zhang and Davidian as "mixture" | Mischaracterized (SUSPECTED) | Keep the citation; describe it as a seminonparametric random-effects density |
| Supplement 1.1: Bartholomew "graded/unfolding" | Misattributed (SUSPECTED) | Bartholomew DJ. Latent Variable Models and Factor Analysis. Griffin; 1987, cited for the fourfold classification |
| Reporting of the pilot | No mixture-reporting standard followed | van de Schoot R et al. Struct Equ Modeling 2017;24(3):451-467 (GRoLTS) |

## 7 Suggestions for strengthening the paper

1. Run the pilot on a design with an off-drug contrast (OL+BDC or
   Hybrid) at $N \in \{35, 70, 140\}$ so that $\beta_{bm:D}$ is the
   estimand in both arms, and make the fixed-effects structure
   identical across arms.
2. Replace the three ad hoc test forms with a principled set: (i) the
   linear-mixed `bm:Dbc` Wald with CR2 and Satterthwaite df, as the
   compendium's F1 analysis uses; (ii) a BLRT or LMR test of one
   versus two classes; (iii) conditional on two classes, a Wald or
   LRT on $\psi_g$; (iv) the within-class `bm:Dbc` test. Report size
   and $\kappa$ for every test before power.
3. Use the protocol's actual random-slopes DGP (file 07) and the
   deterministic scaling DGP (file 02) as two distinct cells, and
   calibrate both to a common marginal correlation or $f^2$ with the
   existing calibration sub-study.
4. Deliver on Aim 1 with one derivation: the expected information for
   $\psi_g$ under the two-class DGP as a function of class separation,
   $N$ and visits, which yields a sample-size curve the pilot can
   check.
5. Report the production Type I error cells already on disk, and
   present all results in a table generated from the stored summary
   with MCSE and denominators.
6. Trim the literature narrative to what the argument uses; move the
   software list to a supplement; add the competing literature of
   Section 3 and position the contribution relative to Muthen et al.
   (2002) and Senn (2016).
7. Decide whether the paper is a methods paper (then the simulation
   program is the paper) or a tutorial/position paper (then it should
   not claim simulation results in the abstract), and submit to a
   venue matching that choice.

## 8 Scope of this review

**Verified by execution.** (a) Recomputation of every rejection rate,
MCSE, convergence rate and denominator in Section "Pilot results"
from `analysis/data/quick-sim/03-latent-replicates.rds`; all quoted
figures (0.197, CI [0.145, 0.250], 0.004, 0.025, 0.004, 0.325, 0.542,
0 of 240, 92.9 percent, 98.3 to 98.8 percent, 558 s, MCSE about
0.030) agree with the stored data. (b) Recomputation of the
pre-flight numbers in Supplement Section 6 from
`output/03-study5-pilot.rds`: 0.250, 0.000, 0.000, 0.010, 0.146,
0.000, 96 percent convergence, 0.05 s and 1.6 s timings, median
gating SE 0.83, median gating coefficient 0.13 (the text says
$-0.10$), sign split 43/53 (the text says 46/54), gating SD 82 (the
text says 63.9), z SD 1.60 (the text says 1.47); the last four differ
from the text, plausibly because the stored object (dated
2026-05-09) is a rerun of the 2026-05-07 pilot the text describes,
which the manuscript should state. (c) One null replicate of the
pilot's design run through `lme_analysis()` to confirm the fitted
formula and regressors (M1). (d) Recomputation of Type I error and
selection rates from `analysis/data/study5/study5-progress.rds` (M3).
(e) Re-derivation of the mixture-induced correlation formula. (f)
`tools/notation-lint.pl` on both Rmd files. (g) Grep for STYLE.md
terms and British spellings.

**Inspected.** `R/lme_analysis.R`, the five driver scripts and the
wrapper in `analysis/scripts/latent-class-mixture-application/`, the
ADEMP protocol and deviations log, `simulation-plan.md`, `README.md`,
`NOTATION.md`, `STYLE.md`, `CLAUDE.md`, and the full text of both
Rmd files and the bibliography. Bibliographic details of the 48
entries were checked against my own knowledge of the works; I did not
resolve DOIs.

**Not checked.** `R/generateData.R` and `buildSigma()` internals; the
`lcmm` package's BIC and convergence definitions; the smoke-test
numbers in Supplement Section 6 (not rerun); the production driver
beyond its header; `bullets.Rmd` and the whitepaper summary; whether
the rendered PDFs match the sources; the `05-`, `06-`, `09-` and
`11-` scripts.

**Literature searches run.** PubMed (via the PubMed MCP tool, no
identifiers sent): "growth mixture model randomized trial treatment
effect latent class heterogeneity" (12 hits; the Muthen 2002 paper
was then located by web search); "N-of-1 trials responder
heterogeneity latent class OR mixture OR individual response" (query
too broad to be useful); "latent class mixed model small sample class
enumeration bootstrap likelihood ratio test simulation" (0 hits);
"Senn mastering variation variance components personalised medicine"
(1 hit, PMID 26415869). Web searches: Davies (1977, 1987) and Hansen
(1996) on nuisance parameters present only under the alternative;
Lubke and Neale (2006); Muthen et al. (2002); Percha et al. (2019).
The remaining references in Sections 3 and 6 are cited from
knowledge of the literature and should be checked against the
originals before use.
