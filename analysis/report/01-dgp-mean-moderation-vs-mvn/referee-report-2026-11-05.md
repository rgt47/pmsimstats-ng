# Referee report (re-review): "Simulating predictive-biomarker interactions in aggregated N-of-1 trials: an extended data-generating process building on Hendrickson et al."
*2026-11-05 11:57 PST*

Manuscript: `analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd`
(3,159 lines; rendered `report.tex` and `report.pdf`, 75 pages, of
2026-11-05 11:20). Target journal: *Statistics in Medicine*, with
JASA/Biometrika standards noted where relevant. Context documents read:
`referee-report-2026-11-01.md`, `revisions/referee-2026-11-01-response-
audit.md`, `revisions/full-read-2026-11-05.md`. Line numbers refer to
`report.Rmd` unless stated.

Evidence labels: **VERIFIED** (re-derived by hand or recomputed in R by
this referee), **INSPECTED** (read in source, code or data and judged;
not recomputed), **SUSPECTED** (plausible, unconfirmed).

## 1. Summary of the manuscript

The manuscript revisits the data-generating process (DGP) of
Hendrickson et al. (2020) for aggregated N-of-1 trials, in which a
baseline biomarker and three latent response components (natural
history, placebo response, drug response) at every visit are drawn
from one multivariate normal distribution, and the predictive
biomarker is a correlation $c_{bm}$ between the biomarker and the drug
response that depends on each visit's exposure state. Using a
paired-difference statistic (E9, E9b) with closed-form moments and
power, it shows that the original rule, which keeps the correlation
constant after first exposure, removes the post-baseline interaction
at any positive carryover half-life, and, through a Schur complement,
that the original compound-symmetry structure admits no joint
distribution above $c_{bm} = 0.256$ without carryover, so that part of
the original grid relied on a silent positive-definite repair. It
proposes an "extended DGP": a correlation that decays at the carryover
half-life, a stop instead of a repair, carryover decaying from the last
on-drug mean, a separable (Kronecker) AR(1) structure at $\rho = 0.7$,
and optional prognostic pathways, measurement error, Weibull decay and
dropout. It verifies the closed form by Monte Carlo, reruns the
original code with only the correlation lines changed, reports an
exploratory stress test of ten analyses under a prognostic biomarker,
and gives closed-form planning sample sizes for E9 together with
simulated power and simulated $N_{80}$ for a mixed-model analysis
(G10) recommended by a companion paper. In its own words, the paper
"extends the original DGP and judges each extension against five
criteria" (existence, identifiability under carryover, fidelity to the
pharmacology, analyzability, calibration status) and "gives the
sample sizes a trialist needs to plan an aggregated N-of-1 trial
under it."

## 2. Overall assessment and recommendation

**Recommendation: major revision.**

The mathematical core is correct, and I re-derived or recomputed it
independently of the project's code (Section 8): the E9 moments, the
slope and noncentrality, the monotonicity argument, the Schur and
Kronecker feasibility bounds, the two-constant eigenvalue condition,
the summed-response factor $q_{ts} = 0.976$, and the E9 planning
sample sizes. The diagnosis of the original simulator (constant
correlation after first exposure; silent repair) is clear and useful
as a correction of the record. The revision has also repaired the most
serious defect of the previous round: the G10 simulations now use the
same response parameters as the closed form, and the rendered G10
numbers agree with the stored output.

What prevents acceptance now is different from last time. First, the
planning deliverable rests on $\rho = 0.7$, which the paper calls a
convention, while the paper's own CSP 563 reanalysis points to a much
slower decay ($\varphi \approx 0.975$ per week on the total score;
the project's document estimates that matching it within the DGP needs
$\rho \approx 0.87$ to 0.91). I recomputed that at $\rho = 0.9$ the
Hybrid feasibility bound falls from 0.480 to 0.340, so two of the three
planning effect sizes become infeasible, and the E9 sample size at
$c_{bm} = 0.25$ falls from 278 to 154 (M1). Second, the sample sizes are
given for G10, an analysis whose close relatives fail the paper's own
prognostic-biomarker stress test (null rejection up to 0.13), without
any caveat in the planning section (M2). Third, the primary validity
table is still not computed under the DGP as the paper defines it
(fixed allocation); the rerun is in progress (M3). Fourth, the
motivating argument for the correlation encoding overstates its
advantages over a random-coefficient mean model, which the paper's own
Discussion shows to be its conditional equivalent and which would
remove the feasibility bound (M4). The literature positioning and
reliance on non-citable internal documents remain to be fixed (M5,
M6). For JASA or Biometrika the methodological contribution would be
too narrow (standard multivariate-normal algebra applied to one
simulator); for *Statistics in Medicine* it is adequate in kind if the
planning results are made robust to $\rho$ and to the analysis.

**Changes since the 2026-11-01 report.** I checked each earlier major
comment against the current text and data, not against the audit's
verdicts.

- M1 (G10 under a different DGP): resolved. The published-set run is
  gated by `PARAM_SET`, the G10 values in the text match
  `planning-g10.csv` and `n80-g10.csv`, and my independent E9 code
  reproduces the closed-form column (VERIFIED).
- M2 (stress test under the rejected carryover rule): partly
  resolved. The anchored run is primary and the caption claim is gone,
  but the run also randomized allocation, which the extended DGP does
  not; the fixed-allocation rerun is pending (carried forward as M3).
- M3 (monotonicity overstated) and M4 (unconditional existence claim):
  resolved by qualification in the Abstract (INSPECTED).
- M5 (sample size for the analysis actually used): resolved for the
  Hybrid design by a probit search (VERIFIED against
  `n80-g10.csv`); not done for the other designs.
- M6 (unexplained 0.126 vs 0.080): resolved by a 4,000-replicate rerun
  (0.0898, MCSE 0.0045; VERIFIED from the stored chunks).
- M7 (novelty search): partly resolved; Scopus and Web of Science not
  searched, NORTA and copula literature not engaged.
- M8 (non-replication of the biomarker association): reported, but the
  implication of the same reanalysis for $\rho$ is not drawn (new M1).
- Declarations remain placeholders, and ten of the previously
  suggested references are still uncited (M5, M6).

## 3. Significance and novelty

What is new, stated narrowly:

1. A closed-form explanation of a published simulator's behavior: the
   power collapse under trivial carryover is entirely a property of the
   correlation rule, through the on/off correlation difference
   $\mathcal G = 1 - \bar g_{\text{off}}$ (Propositions 1 to 3). This is
   a result about one simulator, not new statistical theory.
2. A documented correction of the record: 20 of the 54 nonzero-effect
   matrices of the original Figure 4 grid were silently repaired, and
   the repair at $c_{bm} = 0.6$ changed the simulated DGP.
3. An exact feasibility bound for a structured correlation matrix with
   one covariate row, its Kronecker factorization
   $c^{*}_{bm} = \sqrt{1 - R^2}\,(g^\top A^{-1} g)^{-1/2}$, the per-step
   (Markov innovation) form of $g^\top A^{-1} g$, and the joint
   ellipsoidal region with prognostic correlations. The algebra is
   textbook (Schur complement; Kronecker eigenstructure; AR(1)
   precision); the contribution is its use as an existence check in
   simulation design.
4. A specified, provenance-labeled reference DGP (Table `tab:spec`).

Competing and closely related literature. The paper's Discussion
shows that, conditional on the biomarker, the encoding is a
random-coefficient model with a mean shift $c_{bm} g_t \sigma_{BR} b_i$
and a rank-one reduction of residual covariance. That model class is
standard for individual response in replicated crossover and N-of-1
trials (Senn 2016; Araujo et al. 2016; Senn, Rolfe and Julious 2011)
and for moderators in single-case multilevel models (Van den Noortgate
and Onghena 2003; Moeyaert et al. 2014, 2022), all now cited except
Senn et al. 2011. The quantity $c_{bm}^2$ is the systematic share of
treatment-effect variation of Ding, Feller and Miratrix (2019, JASA).
Grenet et al. (2020) simulate crossover trials with a predictive
biomarker by a mean term. The feasibility problem is the same one
treated for target correlations in NORTA generation (Ghosh and
Henderson 2003), and the real alternative to "stop" is a
structure-preserving repair that keeps a fixed block (Higham, Strabić
and Šego 2016), not the unstructured nearest-matrix repair the paper
argues against. My own searches (Section 8) found no simulation of
N-of-1, crossover or single-case trials that encodes a moderator as an
exposure-indexed covariate-response correlation, so the narrow
priority claim is plausible; the searches were not systematic and do
not establish it.

Bar for the target journal. Items 1 and 2 are valuable to users of
`pmsimstats` and to readers of Hendrickson et al., but are of limited
general interest on their own. Item 3 is the most transferable part.
Item 4 would be a significant contribution only if the reference DGP
were calibrated or its planning outputs shown robust to the
uncalibrated inputs, which is not yet the case (M1, M2). What would
clear the bar: (a) planning results presented as functions of $\rho$
over a data-supported range; (b) a direct comparison of the
correlation encoding with the random-coefficient ("sensitivity")
model at matched implied slope, showing when the two give different
design conclusions; (c) a substantially shorter paper with the
project-specific material moved to supplementary files.

## 4. Major comments

**M1. The planning results depend strongly on $\rho$, the chosen value
of 0.7 is contradicted by the paper's own data source, and the
dependence is not shown (Stage 4, l.1264 to 1332, Table `tab:rhostage`;
planning section l.1843 to 1860 and Tables `tab:planning`,
`tab:planningg10`, `tab:n80g10`; Discussion l.2266 to 2275).**
*VERIFIED (computation) and INSPECTED (`docs/69`).*

Stage 4 adopts $\rho = 0.7$ "as a convention only", adding that 0.7
"has the largest separable feasibility bound (0.474)" of the four
values examined (l.1279 to 1280). Choosing an uncalibrated parameter
partly because it maximizes the range of effect sizes that can be
simulated is not a defensible basis for a reference DGP. The planning
section then reports that CSP 563's total-score correlation decays by
about 0.975 per week with 10 percent noise, and the Discussion that
compound symmetry fits it worst. The repository document the paper
cites (`docs/69`, Section 4.2) goes further: matching the observed 4-
and 20-week correlations within the DGP's structure "would require
$\rho$ of about 0.87 to 0.91 per week", and the on-minus-off contrast
noise under the CSP 563 structure is about one quarter of that under
the extended DGP for the Hybrid paths (0.115 to 0.128 against 0.451 to
0.531). The manuscript does not report either consequence. Stage 4
also says "We have no PTSD symptom series with which to decide whether
[compound symmetry] overstates long-lag correlation" (l.1265 to 1267),
which the planning section and the Discussion contradict.

I recomputed, with my own implementation of the stated DGP
(separable AR(1), $c_1 = 0.2$, decaying correlation, original response
parameters), the Hybrid feasibility bound and the E9 sample size at
$c_{bm} = 0.25$, $t_{1/2} = 0.2$ week:

| $\rho$ | Bound, Hybrid | Bound, OL+BDC | E9 power, $N = 70$ | E9 $N_{80}$ |
|---|---|---|---|---|
| 0.70 | 0.480 | 0.474 | 0.285 | 278 |
| 0.80 | 0.440 | 0.468 | 0.329 | 234 |
| 0.85 | 0.400 | 0.445 | 0.376 | 200 |
| 0.90 | 0.340 | 0.395 | 0.465 | 154 |

The 0.70 row reproduces the paper exactly. At a $\rho$ consistent with
CSP 563, the planning rows at $c_{bm} = 0.35$ and 0.45 are infeasible
in the Hybrid design, and the required $N$ at 0.25 is about 45 percent
lower. Why it matters: the planning table is the paper's practical
deliverable, and Table `tab:rhostage` (0.5 to 0.8) stops short of the
only data-supported region. The caveat that results "should be read
with the sensitivity range of that table" is insufficient when the
range excludes the relevant values.

Remedy: (a) extend Table `tab:rhostage` and the planning table to
$\rho \in \{0.7, 0.8, 0.85, 0.9\}$, or present $N_{80}$ as a function of
$\rho$ in a figure, with infeasible cells marked; (b) state in Stage 4
and the Abstract that the only available PTSD series supports a slower
decay than 0.7, and remove the "no PTSD symptom series" sentence; (c)
drop "largest feasibility bound" as a reason for 0.7, or state plainly
that 0.7 was chosen to keep the original effect sizes simulable; (d)
note that component-level $\rho$ is not identified by a total-score
series, so the CSP 563 value bounds rather than fixes it.

**M2. Sample sizes are given for an analysis that the paper's own
stress test indicates is not valid under the prognostic scenarios the
paper argues are important (planning section l.1991 to 2043; Table
`tab:prognostic` l.1577 to 1587; "Power and robustness" l.1766 to 1780;
Recommendations item 5, l.2295 to 2298).** *INSPECTED (table values
VERIFIED against `summary-full*-anchored.csv` via the rendered text);
the inference about G10 itself is SUSPECTED, because G10 was not run in
the prognostic cells.*

The planning section presents G10 as "the analysis a trial would use"
and gives its power and $N_{80}$. G10 is a random-intercept CAR(1)
model with the half-life chosen by AIC and CR2 inference, with no
expectancy terms and with baseline as an outcome row. In Table
`tab:prognostic`, its two nearest relatives reject the null at 0.099
and 0.114 (fixed half-life, row S) and 0.127 and 0.130 (half-life by
AIC with the interaction in the model, row DS, which corresponds to
G9; the paper reports G9 and G10 within 0.040 of each other in every
planning cell) under the constant natural-history and placebo
pathways. Only the expectancy-term family holds its size, at a cost of
about 13 percent in size-adjusted power (`docs/47`). The planning
section does not mention this; a reader sizing a trial from Table
`tab:n80g10` would not learn that the analysis is anticonservative if
the biomarker is prognostic through a pathway absent at baseline. CSP
563 shows no prognostic association with baseline severity, but says
nothing about association with the placebo response or natural history.

Remedy: run G10 in the four prognostic null cells (the driver
exists); add a sentence to the planning section and the Limitations
stating the condition under which the G10 sample sizes apply; and
either give planning sizes for the expectancy-term analysis (F1 or G1)
alongside G10, or give the inflation factor in $N$ that robustness
costs.

**M3. The primary validity table is still not computed under the
extended DGP as defined (Methods l.635 to 651; Table `tab:spec` row
"Path allocation"; Appendix B.5 l.3099 to 3110; Appendix C).**
*INSPECTED.*

The second run of study 3 changed both the carryover rule and the
allocation (to randomized), and the extended DGP is defined with fixed
allocation. The paper discloses this, and the per-cell changes between
runs are small (largest 0.031; one cell crosses the bold threshold, the
original analysis under the growing pathway, 0.053 to 0.077), so I do
not expect the qualitative conclusions to change. The rerun with
anchored carryover and fixed allocation is reported to be running. It
should replace Tables `tab:prognostic` and `tab:options` before
acceptance, and the two earlier runs should be used to attribute each
change to the carryover rule or the allocation. Two related points:
the stress-test script now defaults to an analysis half-life of 0.2
week (`TH_A <- 0.2`, `11-strawman-stress-test.R` l.244), while the
stored and reported runs used 0.5 week; the Availability section
(l.2558 to 2559) should give the exact invocation that reproduces the
reported tables. And the "growing natural-history" column is a mean
shift of size 0.3 × 10 points per SD at the end of the TV curve
(`b_tv`, l.204 to 208 of the driver), not a correlation, although the
table caption calls every column a "prognostic correlation 0.3".

**M4. The case for the correlation encoding is overstated, and the
paper's own results show it to be a reparameterization of a model that
has none of its limitations (Background, "Why the joint distribution
is built from the design", l.236 to 289; l.277 to 279; Discussion
l.2156 to 2256).** *INSPECTED; the conditional moments VERIFIED.*

The four reasons given for building the joint distribution from the
design do not hold as stated against a properly specified mean-model
alternative:

- Reason 2 says that with a mean term "the residual variance of the
  drug response is the same at every biomarker value" and that the
  simple form "has no parameter" for a biomarker-response correlation.
  Under the correlation encoding the conditional covariance
  $\Sigma_{BR} - c_{bm}^2\sigma_{BR}^2 g g^\top$ also does not depend on
  the biomarker (l.2163 to 2165), so homoscedasticity is shared. A mean
  model with a random participant-by-treatment slope, the form the
  paper itself cites (Senn 2016; Araujo et al. 2016), implies a
  biomarker-response correlation
  $\beta\sigma_B/\sqrt{\beta^2\sigma_B^2 + \sigma_u^2}$; the parameter
  exists.
- Reason 1 (carryover of moderation "would be whatever the analyst
  specifies") applies equally to the decaying correlation, which is
  also specified by the analyst. A mean term $\beta B_i D_{bc,it}$
  represents exactly the same decay.
- Reason 3 (prognostic pathways cannot be generated) is contradicted
  by the stress-test driver, which generates the growing
  natural-history pathway as a mean shift outside the correlation
  matrix (Appendix B.5, l.3084 to 3085).
- Reason 4 treats the absence of a feasibility constraint in the mean
  form as a defect. The bound of the correlation encoding arises
  because the marginal drug-response variance is held fixed and,
  as the Discussion concedes (l.2223 to 2226), the drug-response
  component keeps full variance and serial correlation at off-drug and
  pre-exposure visits. The statement that the bound "is a property of
  the design" (l.277 to 279) is therefore inaccurate: it is a property
  of the design together with these two modeling choices.

Why it matters: the reader is asked to adopt a reference DGP whose
main structural limits (an effect-size ceiling below the motivating
association, an asymmetry between titration and washout, off-drug
noise in a "drug response") are absent from the random-coefficient
model the Discussion proposes (l.2232 to 2256) but leaves unimplemented.
Remedy: rewrite the four reasons as a statement of what the
parameterization buys (an effect size on the scale of explained
variance with fixed marginal variance) and what it costs; correct the
"property of the design" sentence; and implement the sensitivity model
for at least the Hybrid design at matched implied slope, reporting
whether power, feasibility and the design ranking differ. If that is
out of scope, say explicitly that the extended DGP is recommended for
continuity with the original study, not because the encoding is
preferable.

**M5. Literature positioning remains incomplete (Background l.159 to
232; Stage 1 l.707 to 738; Stage 3; Stage 6; Discussion).**
*INSPECTED; existence of the suggested sources checked by web search
(Section 8).*

Ten references recommended in the previous round are still absent,
several attached to specific claims: carryover in crossover trials
(Jones and Kenward 2014, already in the bib; Senn 1992; Senn and
Lambrou 1998); individual response (Senn 2004; Senn, Rolfe and Julious
2011, in the bib); the moderator definition (Kraemer et al. 2002);
$c_{bm}^2$ as an explained share of treatment-effect variation (Ding,
Feller and Miratrix 2019); prognostic biomarkers in interaction tests
(Haller and Ulm 2018; Kent et al. 2020, both in the bib); post-2015
N-of-1 methodology (Schmid and Yang 2022; Kravitz and Duan 2022;
Tueller et al. 2023, in the bib); simulation reporting and post hoc
changes (Siepe et al. 2024; Pawel, Kook and Reeve 2024); and, most
substantively, the repair alternative. Stage 1 argues against an
opaque nearest-matrix repair, which is fair, but not against a
fixed-block shrinking repair (Higham, Strabić and Šego 2016) that keeps
the response block and shrinks only the biomarker row, nor does it
connect the feasibility bound to the NORTA feasibility literature
(Ghosh and Henderson 2003). The novelty search (l.213 to 232) still
omits Scopus and Web of Science and says it does not cover joint
covariate-outcome simulation; one sentence engaging that literature is
needed. Remedy: add the references at the locations listed in Section
6; add a paragraph in Stage 1 explaining why stopping is preferred to a
block-preserving shrink (for instance, that shrinking the biomarker row
lowers the effective $c_{bm}$, which is what the stop reports
honestly); and run and report the Scopus and Web of Science searches
or state their absence as a limitation.

**M6. Load-bearing claims rest on non-citable internal documents, an
unpublished companion paper and non-public data, and the Declarations
are incomplete (l.152, 667, 1773, 1832, 2271, 2755, 2885, 2996;
Declarations l.2563 to 2573).** *INSPECTED.*

The non-replication of the motivating association and every CSP 563
quantity rest on `docs/69`, a repository document describing a
reanalysis of participant-level VA Cooperative Studies Program data
that are not public. The size-adjusted power comparison (0.431 against
0.496; 0.510 against 0.490) rests on `docs/47`; the derivations behind
Appendix A.4's correction claims on a project compendium; the full
derivation of A.5 on `docs/35`; the definition and validation of G10 on
an unpublished companion paper. A journal will require these as
supplementary material, and the CSP 563 reanalysis will need a short
methods description (model, analysis set of 279, covariates) and
confirmation that publishing it is permitted under the data use
agreement. The competing-interests statement is a placeholder; one
author of the original study (Thomas RG) appears to be associated with
this work, and since the paper corrects the original study, the
overlap and whether the original authors were consulted must be
stated. Remedy: move the cited internal documents into a supplement,
or remove the claims that depend on them; complete the Declarations.

## 5. Minor comments

1. **A.2 and A.3, pooled power is an approximation (l.2631 to 2652).**
   The per-path result is exact, but the pooled statistic treats
   between-path mean and slope differences as extra residual variance
   and ignores heteroscedasticity across paths (the code does this
   explicitly, `01-e9-closed-form.R` l.142 to 162). Call the pooled
   value an approximation. My direct Monte Carlo of pooled E9 at
   Hybrid, $\rho = 0.7$, $t_{1/2} = 1$ week (4,000 replicates) gave
   0.177 (MCSE 0.006) against the closed form 0.184, so the
   approximation is adequate (VERIFIED).
2. **Stage 4 contradicts the planning section (l.1265 to 1267 versus
   l.1852 to 1853, l.2271).** See M1 (INSPECTED).
3. **Crossover G10 power at 0 and 0.2 week (Table `tab:planningg10`).**
   At $c_{bm} = 0.25$, $N = 70$, G10 power is 0.187 without carryover and
   0.243 at 0.2 week, 3.1 combined MCSE apart, between two DGPs that are
   practically identical in the crossover (off-drug visits are 2.5 or
   more weeks after stopping). Figure 4B's independent run shows the
   same direction (0.20 against 0.27). The replicate mean estimates are
   $-1.93$ and $-2.09$ against a true $-2.0$ (both within 1.6 SE), so
   chance is the likely explanation, not a defect (VERIFIED from
   `published-planning/trials/`). Common random numbers across
   half-lives would remove such noise from comparisons the text draws.
4. **Multiplicity in single-cell band classifications (Stage 6,
   l.1589 to 1603, l.1703 to 1743).** With ten analyses and fifteen null
   cells, about seven cell-analysis pairs fall outside a two-MCSE band
   by chance. Single-cell exceedances (the free-baseline model at 0.073;
   the AIC expectancy model at 0.078 under Weibull 0.5) should be read
   against that rate, or a multiplicity-adjusted band given (INSPECTED).
5. **Original analysis in study 3 (Table `tab:prognostic`, row 1).**
   The driver fits it with `mmrm` and a compound-symmetry covariance
   with Satterthwaite inference, not the `lmer` random-intercept model
   of the original code. These coincide only when the within-person
   correlation is non-negative; state the implementation (INSPECTED).
6. **Provenance.** The rendered PDF cites commit 091953725 (l.2531),
   but the rendered source was committed afterwards in 717baca00, so the
   cited commit does not contain this version. The G10 runs' `meta.rds`
   has an empty `git_sha` (git exit status 128) (VERIFIED). Record the
   commit from the host at submission time.
7. **Abstract (l.64 to 99).** "Test ten analyses" understates the
   fourteen fitted (l.589); say ten are reported. The Abstract also
   gives the separable bound of 0.474 without the $\rho$ caveat of M1.
8. **$N_{80}$ for G10 (Table `tab:n80g10`).** The probit estimates
   carry Monte Carlo error that the table does not show (disclosed in
   Limitations). Add a parametric-bootstrap or delta-method interval;
   with 500 replicates at six $N$ values this is cheap.
9. **Table `tab:nonpd` is still typed (l.740 to 756).** The counts match
   `rho-stage/nonpd-by-rho.csv` (VERIFIED), but the table should be
   read from that file like the others.
10. **"Property of the design" (l.277 to 279).** See M4; also qualify
    "The bound is low because the correlation profile switches between
    visits that the structure ties together" (l.703 to 705) by noting
    the role of off-drug BR variance.
11. **Estimands (l.469 to 479).** The claim that analyses with
    $D_{bc}$ at the true half-life "recover the on-drug moderation"
    holds because the decaying correlation and $D_{bc}$ share one decay
    function; say so, since under the Weibull option or a separate
    correlation half-life it fails.
12. **Abbreviations (l.2504 to 2518).** SD, SE, CI, ML, REML, AR(1),
    CSP and VA are used but not listed; PD is listed but spelled out
    in the text.
13. **Length and project-specific material.** At 75 pages the paper is
    long for *Statistics in Medicine*. Appendix B code, driver names,
    repository document numbers, code labels in the Glossary and the
    stress-test history (l.657 to 670, Appendix C) belong in
    supplementary material.
14. **Exploratory status (l.657 to 670).** The paragraph discloses post
    hoc additions well; cite Pawel, Kook and Reeve (2024) there, and
    Siepe et al. (2024) beside Morris et al. (2019) at l.564.
15. **Companion-paper title (`references.bib`, `pmsimstats-paper02`).**
    It still uses "covariance moderation", vocabulary the paper has
    retired; align with the current title or with the paper's terms.
16. **Stage 6 feasibility paragraph (l.1639 to 1643).** "A correlation
    with baseline severity has a looser bound, since baseline severity
    is uncorrelated with the response components" is correct, but the
    following A.4 sentence that a biomarker correlated with the drug
    response "already implies a small correlation with the other
    components" (l.2711 to 2712) is not: the DGP sets those correlations
    to zero. What is negative is the off-diagonal of $A$, which reflects
    partial, not marginal, association. Reword (INSPECTED; the bounds
    themselves VERIFIED against `prognostic-joint-ceilings.csv`).

## 6. Missing and questionable references

All 41 cited keys have bibliography entries. The DOIs I checked by
web search (Ding et al. 2019; Higham et al. 2016; Siepe et al. 2024;
Joe 2006) resolve to the stated works. Raskind et al. (2016): the
14-point-per-10-mm-Hg standing SBP figure is confirmed by a secondary
summary of the article, not the full text.

| Claim location | Issue | Suggested source |
|---|---|---|
| Stage 1, l.707 to 738 (stop vs repair) | Only the unstructured repair is argued against | Higham NJ, Strabić N, Šego V. Restoring definiteness via shrinking, with an application to correlation matrices with a fixed block. SIAM Rev. 2016;58(2):245-263. doi:10.1137/140996112 |
| Methods l.543 to 560; Stage 1 (feasibility) | Feasibility of target correlations is studied in NORTA generation | Ghosh S, Henderson SG. Behavior of the NORTA method for correlated random vector generation as the dimension increases. ACM Trans Model Comput Simul. 2003;13(3):276-294 |
| Stage 1 (PD by construction) | Alternative construction that guarantees PD | Joe H. Generating random correlation matrices based on partial correlations. J Multivar Anal. 2006;97(10):2177-2189 |
| Discussion l.2187 to 2189 ($c_{bm}^2$ as explained share) | Uncited; same quantity as systematic variation | Ding P, Feller A, Miratrix L. Decomposing treatment effect variation. J Am Stat Assoc. 2019;114(525):304-317. doi:10.1080/01621459.2017.1407322 |
| Background l.244 to 249; Stage 3 (carryover) | No crossover-carryover methodology cited | Jones B, Kenward MG. Design and Analysis of Cross-Over Trials. 3rd ed. CRC; 2014. Senn S. Is the 'simple carry-over' model useful? Stat Med. 1992;11(6):715-726. Senn S, Lambrou D. Robust and realistic approaches to carry-over. Stat Med. 1998;17(24):2849-2864 |
| Background l.197 to 201; Discussion l.2242 to 2245 (individual response) | Key sources absent | Senn S. Individual response to treatment: is it a valid assumption? BMJ. 2004;329:966-968. Senn S, Rolfe K, Julious SA. Stat Methods Med Res. 2011;20(6):657-666 |
| Background l.206 to 208 (predictive biomarker as moderator) | Definition uncited | Kraemer HC, Wilson GT, Fairburn CG, Agras WS. Mediators and moderators of treatment effects in randomized clinical trials. Arch Gen Psychiatry. 2002;59(10):877-883 |
| Stage 6 l.1460 to 1470 (prognostic biomarker inflates interaction tests) | Uncited though in bib | Haller B, Ulm K. Trials. 2018;19:128. Kent DM et al. The PATH statement. Ann Intern Med. 2020;172(1):35-45 |
| Background l.159 to 176 (N-of-1 methods) | No post-2015 methodological reviews | Schmid CH, Yang J. Bayesian models for N-of-1 trials. Harv Data Sci Rev. 2022 (Special Issue 3). Kravitz RL, Duan N. Harv Data Sci Rev. 2022. Tueller S et al. Behav Res Methods. 2023;55(8):4175-4199 (in bib as `tueller2022`) |
| Methods l.564; l.657 to 670 | Simulation reporting; post hoc changes | Siepe BS, Bartoš F, Morris TP, Boulesteix AL, Heck DW, Pawel S. Psychol Methods. 2024. doi:10.1037/met0000695. Pawel S, Kook L, Reeve K. Biom J. 2024;66(1):2200091 |
| Background l.295 to 297; A.5 (separable covariance) | One source only | Naik DN, Rao SS. Analysis of multivariate repeated measures data with a Kronecker product structured covariance matrix. J Appl Stat. 2001;28(1):91-105 |
| Background l.348 to 354 (symptom vs plasma half-life) | Assumed range has no pharmacodynamic citation | Holford NHG, Sheiner LB. Understanding the dose-effect relationship. Clin Pharmacokinet. 1981;6(6):429-453 |
| l.152, l.1832 (CSP 563 non-replication) | Internal document on non-public data | Supplementary methods note; cite Raskind et al. 2018 for the trial |
| l.1773 (`docs/47`), l.2755 (compendium), l.2885 (`docs/35`) | Non-citable | Supplementary material |

Bibliography hygiene: key `tueller2022` carries year 2023 (cosmetic
until cited); `bell2002` has no DOI (none exists for that journal
issue, acceptable); about 30 entries are uncited and dropped at render,
which is harmless but should be pruned before submission.

## 7. Suggestions for strengthening the paper

- **Make $\rho$ an explicit design input.** A single figure of $N_{80}$
  against $\rho$ (0.6 to 0.95) for the three designs, with the
  feasibility bound overlaid, would replace much of Stage 4 and answer
  M1 directly. Use the CSP 563 structure (exponential plus nugget) as a
  further sensitivity arm, as the project's protocol already proposes.
- **Implement the sensitivity model for one comparison.** The
  random-coefficient model of l.2232 to 2251 removes the bound and the
  titration asymmetry. Comparing it with the correlation encoding at
  matched implied slope in the Hybrid design would turn the Discussion
  into a result and would let the paper recommend an encoding on
  evidence rather than continuity.
- **Plan for the robust analysis too.** Report $N_{80}$ for the
  expectancy-term analysis beside G10 so a trialist can see the price
  of robustness to a prognostic pathway.
- **One generator.** Three generators remain (the closed-form
  `build()`, the stress-test driver, `simulation-core.R`). Converting
  `42-check-param-sets.R` into a `tinytest` unit test that checks E9
  moments for all three would protect the planning numbers against
  future drift.
- **Shorten.** Move code, repository paths and project history to a
  supplement and keep the main text to the closed form, the
  feasibility bound, the extended DGP and a $\rho$-robust planning
  table.

## 8. Scope of this review

**Verified (re-derived by hand or recomputed in R, with code written
for this review from the paper's stated DGP, not sourced from the
project):**

- Gompertz landmarks (inflection 3.8, half-maximum 4.7, 90 percent 9.2
  weeks); re-exposure mean 4.3 and residual 0.16 points.
- A.1 slope, residual variance and noncentrality, including the
  $r/\sqrt{1 - r^2}$ form; A.3 derivative and the small- and
  large-half-life limits; E9b slope $1/(n_{\text{off}} + 1)$ under the
  original rule, consistent with the stored mean estimates.
- A.4 Schur bound; A.5 $[K^{-1}]_{33} = 1/(1 - R^2)$, the Markov
  innovation form of $g^\top A^{-1} g$ and its per-step table; the
  two-constant eigenvalues and threshold 0.111; the compound-symmetry
  decomposition and smallest eigenvalue 0.10; $s^\top K s$ values 368,
  257, 300 and $q_{ts} = 0.976$; A.6 $c_{bm}^2 = f_B f_{BR}$.
- Feasibility bounds: compound symmetry 0.256 (Hybrid, crossover),
  0.262 (OL+BDC), 0.346 at 0.7; separable 0.440 (0.8) and 0.474
  (OL+BDC at 0.7), Hybrid 0.480; compound symmetry under the decaying
  correlation 0.307 (1 week) and 0.369 (2 weeks); eigenvalue sign change
  at 0.999 and 1.001 times the bound; smallest eigenvalues 0.122, 0.193,
  0.28, 0.43; 20 weekly visits 0.1118; daily 0.016 and 0.026.
- On/off correlation differences 0.9997, 0.991, 0.906, 0.754, 0.551 and
  factors 1.8 and 3.3.
- E9 closed-form power 0.290, 0.285, 0.184 (Hybrid, $\rho = 0.7$);
  0.721, 0.335, 0.628; crossover 0.202, OL+BDC 0.214; E9 $N_{80}$ 272,
  278, 480, 84, 146; a direct Monte Carlo of pooled E9 (Minor 1); the
  $\rho$ sensitivity table of M1.
- Planning arithmetic (2.0, 2.8, 3.6, 3.8 points per SD; 21.5 points
  per SD; $\sigma_{BR} \ge 45$; CSP 563 interval $-0.57$ to 0.48;
  1.4 to 2.5 points per 10 mm Hg).
- Stored data against the text: `graded-patch-reps500.csv` (0.668,
  0.638, 0.656; 0.126), the 4,000-replicate chunks (0.0898),
  `nonpd-by-rho.csv` against Table `tab:nonpd`,
  `prognostic-joint-ceilings.csv` against Table `tab:jointceiling`,
  `planning-g10.csv` and `n80-g10.csv`; replicate-level G10 output for
  four planning cells (Minor 3); `meta.rds` `git_sha`.

**Inspected, not recomputed:** the stress-test driver
(`11-strawman-stress-test.R`, cells, DGP, analyses), the stress-test
tables as rendered against their summary files (via the Rmd code that
reads them), Appendix B code patches, the repair magnitudes at
$c_{bm} = 0.3$ and 0.6, the G10 implementation in
`simulation-core.R` (half-life selection only), `docs/69` Sections 1
to 4.

**Not checked:** the study 1 Monte Carlo agreement count (10 of 278),
the original-code reruns themselves, Figures 3 to 6 against their
data, the CSP 563 reanalysis (data not available to me), the companion
paper's results, and whether the pending fixed-allocation stress-test
rerun changes any cell.

**Literature searches run (web search, no email or API identifiers
sent):** simulation of aggregated N-of-1 trials with a predictive
biomarker and carryover (returned only Hendrickson et al. 2020);
maximum feasible covariate-outcome correlation and Schur-complement
feasibility in simulation; NORTA feasibility (Ghosh and Henderson);
Bayesian hierarchical models for series of N-of-1 trials with
covariate heterogeneity (Schmid and Yang 2022 and theses; no
covariate-moderator simulation found); crossover simulation with
biomarker-by-treatment interaction and correlated random effects
(found gene-by-treatment crossover simulation and carryover-modeling
papers, none with a correlation-encoded moderator); Higham-Strabić-Šego
fixed-block shrinking; Ding-Feller-Miratrix; Siepe et al. simulation
template; Raskind et al. 2016 effect size; sample size for moderators
in series of N-of-1 trials (found N-of-1 sample-size work and
cluster-trial moderator power formulas, none for N-of-1 moderators);
Joe 2006. None of these searches found a prior exposure-indexed
correlation encoding of a moderator; they were not systematic database
searches, so precedents cannot be excluded. Scopus and Web of Science
were not searched by this referee either.
