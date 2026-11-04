# Referee report: "Simulating predictive-biomarker interactions in aggregated N-of-1 trials: an extended data-generating process building on Hendrickson et al."
*2026-11-01 15:55 PST*

Manuscript: `analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd`
(rendered 2026-11-01, after the G9-to-G10 switch and the prose pass).
Target journal orientation: *Statistics in Medicine* or *Trials*
(simulation methodology for clinical trial design).

## 1. Summary of the manuscript

The manuscript studies the data-generating process (DGP) of
Hendrickson et al. (2020) for aggregated N-of-1 trials, in which a
predictive biomarker is encoded as a correlation $c_{bm}$ between the
baseline biomarker and a latent drug-response component inside one
multivariate normal vector assembled from each participant's exposure
history. Using a paired-difference interaction statistic (E9, E9b)
whose distribution is available in closed form, it shows that the
original rule, which holds the correlation constant after first
exposure, removes the interaction at any positive carryover half-life,
so that power equals the test size; and, via a Schur complement, that
the original compound-symmetry structure admits no joint distribution
for $c_{bm} > 0.256$, so the original effect size of 0.6 relied on a
silent nearest-positive-definite repair. It proposes an "extended DGP":
a correlation that decays at the carryover half-life, a stop instead of
a repair, carryover decaying from the last on-drug mean, separable
(Kronecker) AR(1) correlation at $\rho = 0.7$, and optional prognostic
pathways, measurement error, Weibull decay and dropout. It verifies the
closed form by Monte Carlo, re-runs the original code with only the
correlation lines changed, runs an exploratory stress test of fourteen
analyses under a prognostic biomarker, and gives closed-form planning
sample sizes together with simulated power of a recommended mixed-model
analysis (G10, from the companion paper 02).

## 2. Overall assessment and recommendation

**Recommendation: major revision.**

The core of the paper is sound and useful. The closed-form treatment
of the paired-difference statistic is correct (I re-derived the
moments, the slope, the noncentrality and the monotonicity argument;
Section 8), and it turns a puzzling published finding (power collapsing
under trivial carryover) into a transparent consequence of one coding
rule. The Schur-complement feasibility bound and its Kronecker
factorization are correct and give a practical existence check that
simulation papers in this area lack. I reproduced every feasibility
bound and on/off correlation difference quoted in the text.

Three problems prevent acceptance as it stands. First, and most
serious, the simulated G10 results in the planning section and in the
G10 panel of Figure 4B were generated under response parameters that
differ from the DGP the paper specifies and from the parameters used
for the closed-form E9 numbers they are compared with (Major comment
1). The paper's headline planning conclusion, that the closed-form
sample sizes are conservative for the recommended analysis, is largely
an artifact of that mismatch. Second, the stress-test table reports the
first run under the original carryover rule, not the extended DGP the
paper recommends, and its caption's claim that the rerun changed no
verdict is contradicted by the stored rerun (Major comment 2). Third,
the motivating biomarker association did not replicate in the
project's own reanalysis of the largest prazosin trial, and the paper
does not say so (Major comment 8). Fourth, the novelty and literature
positioning rest on a single PubMed query, one Background sentence is
contradicted by a reference already in the bibliography, and several
citations do not support the claims attached to them (Sections 3 and
6). For a top-tier statistics journal
(JASA, Biometrika) the methodological contribution is narrow: the
theory is elementary multivariate-normal algebra applied to one DGP. It
fits *Statistics in Medicine* or a trial-methodology journal well,
provided the correctness issues are fixed.

## 3. Significance and novelty

Closest prior work. The paper's own Discussion shows that,
conditional on the biomarker, the correlation encoding is a
random-coefficient (cross-level interaction) model with a reduced
residual variance; the model class is standard and what is new is the
parameterization, a correlation with fixed marginal variance indexed by
exposure state. Ding, Feller and Miratrix (2019) define the same
explained-variance quantity as $c_{bm}^2$. Grenet et al. (2020)
simulate crossover trials with a predictive biomarker, but with a
mean-term interaction. Copula and NORTA methods generate covariates and
outcomes jointly from target correlations, and the NORTA literature
treats infeasible target matrices (Ghosh and Henderson 2003). No
exposure-indexed covariate-response correlation was found in any N-of-1,
crossover or single-case simulation, so the narrow claim is plausible;
the search behind it is not adequate to assert it (M7).

What is new:

- A closed-form account of why the original DGP's power collapses
  under carryover, isolating the correlation rule as the cause
  (Propositions 1 to 3). This is a new result about a published
  simulator, not new statistical theory.
- The observation that the original effect sizes exceed the
  compound-symmetry feasibility bound and were silently repaired, with
  a count of repaired matrices. This is a useful and verifiable
  correction of the record.
- An exact feasibility bound for structured correlation matrices with
  a biomarker row, including the Kronecker factorization
  $c^{*}_{bm} = \sqrt{1 - R^2}\,(g^\top A^{-1} g)^{-1/2}$ and the joint
  ellipsoidal region with prognostic correlations. The algebra is
  standard (Schur complement, Kronecker eigenstructure); the application
  to simulation design is the contribution.
- A specified, labeled reference DGP with a provenance table
  (Table `tab:spec`).

What is not new, and the paper says so: multivariate-normal simulation
of repeated measures with carryover, the pharmacometric decomposition,
separable covariances, and nearest-correlation-matrix repair.

Contemporariness. The paper engages ADEMP (Morris et al. 2019) and
CR2 inference, but reads as positioned inside one project rather than
the current N-of-1 methods literature. A current referee would expect
engagement with Bayesian and multilevel N-of-1 analysis and its
simulation literature, with the heterogeneity-of-treatment-effect
literature on whether replicated crossovers identify individual
response, and with the copula and latent-variable approaches to joint
simulation of covariates and outcomes, which are the natural
alternatives to placing the biomarker in the response covariance. See
Section 6.

## 4. Major comments

**M1. The G10 simulations use a different DGP from the one the paper
specifies and from the closed-form numbers they are compared with
(planning section; Table `tab:planning` last column; Table
`tab:planningg10`; Figure `fig:fig4bg10`; Limitations item 1;
Conclusions).** *Verified.*

The closed-form E9 numbers come from
`hendrickson-problems/04-power-simulation.R`, whose `build()` uses the
original parameters the Methods describe: natural-history (TV) and
placebo (PB) components on the PB Gompertz curve, SDs 10 and
$10\eta_t$, BR SD 8, biomarker SD 15.4, baseline SD 18.5. The G10
column and the G10 Figure 4B panel come from
`19-g9-planning-and-fig4b.R`, which draws data with
`simulation-core.R::default_resp_param()` and
`default_baseline_param()`: TV on the BR curve (max 10.99, rate 0.42)
with SD 6, PB SD $2\eta_t$, biomarker SD 1, baseline mean 70 and SD 10
(`build_sigma_matrix()` reads these SDs directly). TV and PB together
carry far less noise.

I ran a Monte Carlo of E9 itself (1,000 replicates; post-baseline on
minus off mean, OLS on the biomarker, paths pooled) with the generator
of driver 19 at Hybrid, $N = 70$, $c_{bm} = 0.25$, $t_{1/2} = 0.2$
week. E9 power is 0.549 (MCSE 0.016). The paper compares G10's 0.649
with the closed-form 0.285. On a common DGP the G10 advantage over E9
at this cell is about 0.10, not 0.36. Consequences:

- "At $N = 70$ G10 has higher power than E9 in every cell" and the
  quoted contrasts (0.649 vs 0.285; 0.403 vs 0.202) compare two DGPs.
- "The closed-form E9 rule is therefore conservative for G10"
  (power 0.984 to 1.000 at the E9 sample size) is confounded. The
  conservatism measured is mainly that of lower noise variance.
- The Figure 4B G10 paragraph says the two Figure 4B versions differ
  "in its correlation structure ... so the two figures are not matched
  cell by cell". They also differ in the TV curve, the TV and PB SDs,
  and the biomarker and baseline scales.
- The same parameter set is used by the companion paper 02, so the two
  papers' "extended DGP" are not the same process in their response
  parameters, although paper 02 cites paper 01 for it.

Remedy: rerun driver 19 with paper 01's stated parameters (the
`build()` parameters), or state and justify a single parameter set
used throughout both papers and recompute the closed form under it.
Then rewrite the planning results, the G10 Figure 4B text, the
Limitations bullet and the Conclusions sentence. Add a check that the
driver's generator reproduces the closed-form moments of E9 (the
existing study 1 check does this for `build()`, not for
`simulation-core.R`).

**M2. The stress-test table reports the original carryover rule, and
the rerun does change a verdict (Stage 6, Table `tab:prognostic` and
its caption; Methods, "Which form of the DGP each study used";
Appendix B.5).** *Verified from stored output.*

The table is the first run (`strawman-stress/summary-full.csv`), which
uses the original compounding carryover rule that Stage 3 replaces.
Its "Original analysis" and "CAR(1), CR2" rows match that file exactly.
The rerun under the extended DGP's own rule
(`summary-full-anchored.csv`) gives, for the original analysis,
0.077 under the growing natural-history pathway (first run 0.053) and
0.119, 0.167, 0.069 in the other three columns (first run 0.107,
0.159, 0.075); for the CAR(1) CR2 model, 0.114 under the placebo
pathway (first run 0.100). The original analysis under the growing
pathway therefore moves from inside the band to above both the band
and the table's own bold threshold. The caption's "changed no
analysis's verdict" is not correct at the cell level. The paper also
recommends a DGP and then reports its key validity table under a rule
it has rejected.

Remedy: report the anchored rerun as the primary table (it is the
extended DGP), give the first run in an appendix if wanted, and replace
the caption claim with the per-cell changes. Note that the rerun also
randomized allocation, so it changes two things at once (Methods); say
which change is responsible, or rerun with one change at a time.

**M3. Monotonicity and "returns to the no-carryover value" are stated
more strongly in the Abstract than proved (Abstract, Results;
Appendix A.3).** *Inspected.*

Proposition 3 proves strict monotonicity for the path-stratified
statistic. The pooled statistic is checked on a 120-point grid, and
the mixed model only at the simulated half-lives (the text says so).
The Abstract says power "declines monotonically at longer ones" without
qualification. Either prove the pooled case (the between-path term
depends on $t_{1/2}$ only through the carryover means, so a bound on its
derivative may suffice) or qualify the Abstract.

**M4. The Abstract's existence claim is unconditional but the result is
conditional (Abstract, Results sentence 4; Stage 1).** *Inspected.*

"Under the original compound-symmetry structure no joint distribution
exists for correlations above 0.256" holds without carryover, or under
the decaying rule. Under the original rule at positive half-lives the
original DGP passes many such cells (Table `tab:nonpd`: 12 of 40
repaired at $c_{bm} = 0.3$ on the extended grid, not 40), as Stage 1
itself explains. State the condition in the Abstract.

**M5. The planning table is the practical deliverable but gives sample
sizes for a statistic nobody would analyze with (planning section;
Limitations).** *Inspected.*

The paper gives closed-form $N$ for 80% power of E9 and then reports
that the recommended analysis has power near 1 at those $N$ and that
"we did not search" for the G10 sample size. After M1 is fixed the gap
will be smaller, but a trialist still needs $N$ for the analysis they
will run. Remedy: search for the G10 $N_{80}$ in at least the Hybrid
cells (a bisection on $N$ with a few thousand replicates per step is
cheap relative to the runs already made), and present E9 as the
analytical approximation it is, with its measured bias.

**M6. An unexplained discrepancy between two implementations of the
original DGP (Stage 2, paragraph on the original analysis).**
*Inspected.*

The reconstruction and the separately seeded run of the original code
differ at 0.1 week by about 2.4 combined Monte Carlo SE, and the text
says the source was not traced. Since the paper's argument rests on
reproducing the original code exactly (Table `tab:nonpd` says the
repair counts reproduce "exactly"), a nominally significant difference
between two supposed implementations of the same process should be
explained or the comparison rerun with more replicates. It does not
affect the conclusion (both are near the test size), but it is the kind
of loose end a referee will ask about.

**M7. The novelty claim rests on one database query (Background, "The
priority claim above rests on a PubMed search").** *Inspected.*

A single PubMed query returning one record does not establish that
the encoding "has no precedent in the N-of-1 literature". Joint
multivariate-normal simulation of covariates with latent outcomes, and
moderation expressed through a covariance, appear in other literatures
(Section 6). Remedy: report the search strings and databases, extend
to Web of Science or Scopus and to the simulation-methods literature,
and soften the claim to "in the N-of-1 literature we searched".

**M8. The motivating biomarker association did not replicate in the
project's own reanalysis of the largest prazosin trial (Background
l.115-119; planning section l.1326-1345; "What a DGP must satisfy",
criterion 5).** *Inspected (`docs/69` of the repository; the
reanalysis itself was not rerun).*

The paper motivates the biomarker with the 67-participant secondary
analysis of Raskind et al. (2016), 14 CAPS points per 10 mm Hg, and
calibrates the planning discussion against it. The project's
reanalysis of CSP 563 (304 veterans; `docs/69`, 2026-10-30) estimates
the standing-SBP by arm interaction at +0.25 points per 10 mm Hg
(95% CI -2.72 to 3.21), which excludes -14. The same document gives a
biomarker SD of 14.2 mm Hg (manuscript: 15.4), a baseline SD of 16.0
(manuscript: 18.5), a placebo-arm improvement about 1.8 times the
DGP's, an SD of change from baseline of 17 to 24 that grows over time
(DGP: 16.0, constant), and a total-score correlation decaying at about
0.975 per week with about 10 percent measurement noise. None of this is
in the manuscript, which states that no parameter has been estimated
from prazosin trial data.

Remedy: report the failure to replicate in the Background and the
planning section; replace "no parameter has been estimated" with what
CSP 563 does and does not identify (the total-score variance and
correlation are estimable; the TV/PB split, the expectancy weight and
the carryover half-life are not); and either run the CSP 563-calibrated
sensitivity arm that `docs/69` Section 9 lays out, or state it as
future work.

## 5. Minor comments

1. **Table `tab:estimands`.** The stress-test row says these analyses
   use $D_{bc,it}$ at the true half-life, but the stress test includes
   analyses with an AIC-chosen half-life (Table `tab:prognostic` rows
   3, 5 and 9). Split the row.
2. **Methods, the original DGP.** The text gives TV "the same curve as
   PB and SD 10". State whether the companion paper uses these values
   (see M1); at present a reader of both papers cannot tell.
3. **Stage 2, step-function threshold.** "Above about 0.007 week ...
   below that the off-drug mean underflows". This is a property of
   double-precision arithmetic in the original code, not of the DGP;
   say so in one clause, or move it to Appendix B.
4. **Stage 4.** "We take $\rho = 0.7$ as a convention ... the original
   0.8 lowered by 0.1, which makes the original effect size of 0.3
   feasible under the original structure". At $\rho_{cs} = 0.7$ the CS
   bound is 0.346, so 0.3 is feasible; but the adopted structure is
   separable AR(1), where 0.3 is feasible at 0.8 already (bound 0.440).
   The stated reason for 0.7 therefore does not apply to the structure
   adopted. Give the actual reason or call it a convention only.
5. **Stage 5, property 4.** "This follows from the algebra; its effect
   on test calibration has not been measured." Give the derivation in
   Appendix A.5 explicitly (the summed response is AR(1) only up to the
   expectancy scaling of PB, so it is not exactly AR(1) in the Hybrid
   design).
6. **Stage 6, size more generally.** "Both values of 0.057 lie above
   the pooled band": the band is 0.046 to 0.054 for 11,000 replicates
   at nominal 0.05; state the interval used.
7. **Planning section.** "about a quarter of the sample [@brookes2004]"
   and "about four times" in the Background: Brookes et al. give a
   factor of about four for an interaction of the same size as the
   main effect under specific conditions (two equal subgroups). State
   the condition.
8. **Glossary.** "Original rule (`step`)" and "Decaying correlation
   (`graded`)" use code labels that the figures print ("Step",
   "Graded"); the figure captions explain this, but the figures
   themselves would be clearer relabeled.
9. **Declarations.** Competing interests and the consultation of the
   original authors are placeholders; they must be completed before
   submission, especially since the paper corrects the original study.
10. **Reproducibility.** Several tables are typed into the Rmd
    (Tables `tab:closedform`, `tab:patch`, `tab:rhostage`,
    `tab:ceilings`, `tab:prognostic`, `tab:options`, and the E9 columns
    of `tab:planning`). I checked the E9 planning values against
    `planning-n80.csv` (all 36 match) and the first-run stress-test
    rows against `summary-full.csv` (match). Reading every table from
    its stored output, as the G10 tables already do, would remove the
    risk of drift and would have exposed M2.

## 6. Missing and questionable references

All 35 cited keys have bib entries, and the 29 DOIs on cited entries
resolve with matching metadata. The problems are with what some
citations are used to support, and with literature that is absent.

**Citations that do not support the attached claim, or only partly**

| Location | Issue | Remedy |
|---|---|---|
| Background, l.236 (li2023) | Li et al. model HTE as a mean-model interaction; their correlations are ICCs of the modifier and of the outcome given the modifier. There is no modifier-outcome correlation, so this is not a precursor of the correlation encoding (full text read). | Reword to "designs for HTE analysis in cluster trials" or drop from the novelty paragraph. |
| Background, l.119 (raskind2003, 2013, 2018) | Efficacy trials, cited for "blood pressure is therefore a candidate predictive biomarker". They do not analyze BP moderation, and Raskind 2018 (NEJM) was null overall, which the paper does not say. | Cite raskind2016bp alone for the BP finding; state the 2018 null result. |
| Methods, l.363 (higham2002) | `corpcor::make.positive.definite` implements Higham's 1988 nearest symmetric positive-definite algorithm, which does not keep a unit diagonal; it is not the nearest-correlation-matrix method of Higham 2002. | Cite Higham NJ. Linear Algebra Appl. 1988;103:103-118, doi:10.1016/0024-3795(88)90223-6, and correct the description. |
| Background, l.86 (sargent2005) | No "four times" factor found in it (abstract only). | Keep brookes2004 for the factor, with its condition (Minor 7). |
| Planning, l.1394 (brookes2004) | Carrying the factor of 4 to a continuous-biomarker slope "of the same standardized size" extrapolates beyond Brookes's binary-subgroup setting. | State the extrapolation. |
| l.156, l.177 (zucker1997, senn2019samplesize) | Cited for what simulation studies usually do; neither is a simulation study (abstracts only). | Cite simulation papers (chen2014, wang2019, percha2019) for that claim. |
| A.4, l.2048 (liu2009baseline) | Supports the analogy (baseline as covariate or response), not the specific constancy constraint stated. | Say "analogous to". |
| Discussion, l.1552 (pmsimstats-paper03) | Unpublished, cited for mixture models needing larger samples. | Cite Bauer and Curran 2003, Lubke and Muthén 2005, Nylund et al. 2007 (already in the bib, uncited). |

**Missing references** (full citations in the audit file; metadata
checked on Crossref unless noted)

| Claim location | Missing source |
|---|---|
| Background l.131: "None of them represents a baseline biomarker that modifies the drug response" | Grenet G et al. Basic Clin Pharmacol Toxicol. 2020;126(1):59-64, doi:10.1111/bcpt.13293 (in the bib, uncited). A crossover simulation with a predictive biomarker (mean-term encoding). The sentence is false as written. |
| Carryover throughout; Stage 3 | Jones B, Kenward MG. Design and Analysis of Cross-Over Trials, 3rd ed. CRC; 2014 (in bib, uncited). Senn S. Stat Med. 1992;11(6):715-726, doi:10.1002/sim.4780110603. Senn S, Lambrou D. Stat Med. 1998;17(24):2849-2864. |
| Individual response, l.114, l.156, Discussion | Senn S. BMJ. 2004;329:966-968, doi:10.1136/bmj.329.7472.966. Senn S, Rolfe K, Julious SA. Stat Methods Med Res. 2011;20(6):657-666 (in bib, uncited). |
| Definition of a moderator / predictive biomarker | Kraemer HC et al. Arch Gen Psychiatry. 2002;59(10):877-883, doi:10.1001/archpsyc.59.10.877. |
| Discussion: $c_{bm}^2$ as the share of drug-response variance explained | Ding P, Feller A, Miratrix L. J Am Stat Assoc. 2019;114(525):304-317, doi:10.1080/01621459.2017.1407322. Their systematic-variation $R^2$ is the same quantity. |
| Stage 6, prognostic biomarker | Haller B, Ulm K. Trials. 2018;19:128 (in bib, uncited). Kent DM et al. (PATH). Ann Intern Med. 2020;172(1):35-45 (in bib, uncited). |
| N-of-1 methodology after 2015 | Schmid CH, Yang J. Harvard Data Sci Rev. 2022, doi:10.1162/99608f92.3f1772ce. Kravitz RL, Duan N. Harvard Data Sci Rev. 2022, doi:10.1162/99608f92.901255e7. Daza EJ. Methods Inf Med. 2018;57(1):e10-e21. Mirza RD et al. J R Soc Med. 2017;110(8):330-340. Stunnenberg BC et al. JAMA. 2018;320(22):2344-2353. Tang J, Landes RD. PLoS ONE. 2020;15(2):e0228077 and Tueller S et al. Behav Res Methods. 2023;55(8):4175-4199 (both in bib, uncited). |
| Single-case multilevel models with moderators (novelty claim) | Van den Noortgate W, Onghena P. Behav Res Methods Instrum Comput. 2003;35(1):1-10. Moeyaert M et al. J Sch Psychol. 2014;52(2):191-211. |
| Simulation reporting; post hoc additions to study 3 | Siepe BS et al. Psychol Methods. 2024, doi:10.1037/met0000695. Pawel S, Kook L, Reeve K. Biom J. 2024;66(1):2200091, doi:10.1002/bimj.202200091. |
| Stage 1 repair versus stop | Higham NJ, Strabić N, Šego V. SIAM Rev. 2016;58(2):245-263, doi:10.1137/140996112, a repair that keeps a fixed block and shrinks the rest, which is the real alternative to stopping; the paper should say why stopping is preferred. Ghosh S, Henderson SG. ACM TOMACS. 2003;13(3):276-294 (feasibility of target correlations in NORTA). Joe H. J Multivar Anal. 2006;97(10):2177-2189 (positive definiteness by construction). |
| Separable covariance, l.234 and A.5 | Naik DN, Rao SS. J Appl Stat. 2001;28(1):91-105. |
| Symptom half-life longer than plasma half-life, l.283-286 | Holford NHG, Sheiner LB. Clin Pharmacokinet. 1981;6(6):429-453 (effect-compartment models). |

**Bibliography hygiene.** DOIs missing on raskind2003
(10.1176/appi.ajp.160.2.371) and raskind2013
(10.1176/appi.ajp.2013.12081133). Key `tueller2022` carries year
2023. `sturdevant2021carryover` gives the first author as "Sarah G.";
Crossref has "S. Gwynn". Journal names vary (PLoS ONE / PLoS One;
Healthcare / Healthcare (Basel)). gartner2023 lacks an issue number.
38 of 73 entries are uncited (dropped by pandoc, so harmless), but
`pmsimstats-paper11` should be removed from this paper's bib.

## 7. Suggestions for strengthening the paper

- **One parameter set, one generator.** Make the closed form, study 3
  and the G10 runs draw from one function with one parameter table,
  and add a unit test that the generator reproduces the closed-form
  mean and variance of E9 at a few cells. This would have caught M1.
- **Mixed-model power without the paired-difference proxy.** The
  closed form is exact for E9 only. A small simulation of G10 power
  over the half-life grid on the extended DGP would let the Stage 2
  monotonicity claim be stated for the analysis a trialist would use.
- **The sensitivity model.** The Discussion proposes a
  random-coefficient (individual sensitivity) DGP that removes the
  feasibility bound and the titration asymmetry. Implementing it for
  one cell, and comparing power with the correlation encoding at
  matched implied slope, would turn a suggestion into a result and
  answer the obvious referee question of whether the bound is a
  property of the encoding or of the trial.
- **Calibration.** The paper is candid that no parameter is
  estimated. A short appendix using any available weekly PTSD symptom
  series (the project's own CSP 563 analysis, `docs/69`, appears to be
  one) to estimate $\rho$ and $\sigma_{BR}$ would address the most
  sensitive inputs.

## 8. Scope of this review

**Verified (re-derived or executed):**

- Appendix A.1 variance formula under compound symmetry, term by term;
  Proposition 2 slope and its E9 and E9b forms; A.2 noncentrality and
  its $r/\sqrt{1 - r^2}$ form; A.3 derivative and Proposition 3; the
  small- and large-half-life limits.
- The Schur-complement bound (A.4) and the Kronecker factor
  $\sqrt{1 - R^2}$ with $R^2 = 2c_1^2/(1 + c_1)$; the conditional
  mean and covariance in the Discussion; $c_{bm}^2 = f_B f_{BR}$ for
  a two-class mixture with within-class independence.
- Feasibility bounds by bisection on the smallest eigenvalue with the
  closed-form `build()`: CS 0.256, 0.256, 0.262; two-constant 0.291
  ($\rho = 0.8$) and 0.452 ($\rho = 0.7$); separable 0.440 and 0.474
  (Hybrid 0.480). The Kronecker closed form gives 0.480 for Hybrid at
  $\rho = 0.7$, matching the bisection.
- On/off correlation differences 0.9997, 0.991, 0.906, 0.754, 0.551
  (Hybrid, allocation-weighted) and the $1/\mathcal G^2$ factors 1.8
  and 3.3.
- Planning arithmetic (slope per SD, 0.13 points per mm Hg, 21.5 points
  per SD, $\sigma_{BR} \ge 45$).
- The E9 planning values against `planning-n80.csv` (36 of 36).
- M1: parameters of both generators read from code; E9 Monte Carlo
  under driver 19's generator (0.549, 1,000 replicates).
- M2: first-run and anchored-rerun stress-test summaries read from
  `strawman-stress/`.

**Inspected, not verified:** Appendices B.1 to B.5 (code patches and
package changes), Table `tab:patch` and Figure `fig:fig4bmc` numbers
(not regenerated), Table `tab:options`, the repair counts of Table
`tab:nonpd`, and the stress-test analyses' implementation.

**Not checked:** the Monte Carlo study 1 agreement count (10 of 278
cell-variants); the CSP 563 material; the companion paper's results
beyond the parameter question in M1.

**Literature searches** (run by a delegated audit, results checked
for consistency, not re-run): 22 web searches covering crossover and
N-of-1 data simulation, carryover models, interaction power,
heterogeneity of treatment effect in replicated crossover, nearest
correlation matrices and positive-definiteness of structured
matrices, Kronecker covariances, pharmacometric disease-progression
models, covariance- or correlation-based moderation in SEM and
psychometrics, copula and NORTA generation, post-2015 N-of-1 methods,
simulation reporting, CR2 inference and prazosin PTSD trials; about 25
Crossref bibliographic queries and DOI checks of 32 DOIs; full text
read for Li et al. 2023 and Chen and Chen 2014. Verdicts on Sargent
2005, Wang 2019, Senn 2019, Zucker 1997, Liu 2009 and Raskind 2016 rest
on abstracts or summaries. No systematic database search was run, so
precedents for the novelty claim cannot be ruled out.
