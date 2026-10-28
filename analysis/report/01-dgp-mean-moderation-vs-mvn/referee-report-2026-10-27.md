# Referee report: paper 01, refined data-generating process for aggregated N-of-1 trials
*2026-10-27 11:34 PDT*

Manuscript: 'Simulating predictive-biomarker interactions in aggregated
N-of-1 trials: a refined data-generating process following Hendrickson
et al.' (`report.Rmd`, with its Appendices A and B and the companion
`appendix-b.Rmd`). Target journal as instructed: Trials (Methodology
article). Reviewed as it stands after the 2026-10-27 removal of all
mean-moderation material; the absence of that comparison is not
criticized. The studies are treated as exploratory throughout.

Evidence labels used below: VERIFIED (I re-derived or executed it),
INSPECTED (I read the source and judged it), SUSPECTED (plausible, not
confirmed).

## 1. Summary of the manuscript

The paper revisits the simulation framework of Hendrickson et al.
(2020; 'RH'), which encodes a biomarker-by-treatment interaction as a
correlation c_bm between the biomarker and the drug-response component
of a multivariate normal symptom model. The authors identify two
problems in RH's code at commit 58b32a9: a coupling rule that sets the
biomarker correlation to its full on-drug value at every off-drug visit
with a nonzero residual drug mean, and a silent repair of correlation
matrices that are not positive definite. Using a paired-difference
statistic (E9, E9b) they derive closed-form moments and noncentral-t
power (Propositions 1 to 3), show that the step rule removes the
interaction at every positive carryover half-life, and show that
'graded' coupling (decaying with the residual drug effect) restores a
smooth decline. A Schur-complement argument gives an exact ceiling for
c_bm (0.256 under RH's compound symmetry, 0.440 to 0.480 under a
separable AR(1) structure). Patched RH code is rerun, and a 41-cell
exploratory stress test of 14 analysis models, including prognostic
biomarker scenarios, ends in a recommendation of a model with
expectancy terms. The paper proposes a 'refined process' as the
reference DGP for such simulations.

## 2. Overall assessment and recommendation

Recommendation: major revision. For Trials specifically: not suitable
in its present form; suitable after the reframing and repairs listed in
Major Comments 1 and 2 (the first concerning fit and audience, the
second concerning consistency of the deliverable).

Justification. The central technical finding is sound and, in my
assessment, correct: with step coupling the contrast weights on the
coupling profile sum to zero for E9 only if the off-drug coupling is
below one, so the interaction vanishes for any t_half above numerical
underflow. I re-derived Proposition 2 and the ceiling formula and
reproduced every closed-form number in Tables 'closedform' and
'ceilings' and the c_bm = 0.2 values quoted in the text. The
PD-ceiling and separability algebra (Appendix A.5) is also correct,
with one stated number wrong (Minor 1). The paper is candid about its
exploratory status in the stress test. The weaknesses are of a different
kind. (i) The contribution is a correction to one prior simulation
framework; its wider value depends on guidance that trialists can use,
and that guidance is thin (Major 1). (ii) The refined process as
specified in Table 'spec' is not what the package implements, and the
reported stress-test results derive from earlier variants of the
process (Major 2). (iii) At the recommended pharmacological half-lives
the refined process is numerically indistinguishable from no carryover,
which weakens the framing of 'restoring' a decline (Major 3). (iv) The
stress test is post hoc and partly circular (Major 4). (v) The
manuscript contains internal revision notes, unfinished declarations and
a stale companion appendix, so it is not submission-ready (Major 5 and
Minor comments).

## 3. Significance and novelty

What is new, stated precisely:

- A closed-form identification of why RH's power collapses under
  carryover (coupling gap), with exact moments and power for E9/E9b.
  The algebra is elementary but it is correct and it explains a
  published, surprising result. Clears the bar for a methods note; it
  is not a deep result.
- An exact feasibility ceiling for c_bm (Schur complement) and a
  factored form under separability, with a per-step cost table. The
  Schur-complement condition itself is standard linear algebra; the
  application to a named simulation framework and the Kronecker
  factorization (c* = sqrt(1 - R^2) (g' A^-1 g)^(-1/2), VERIFIED
  numerically to 1e-16 against the direct computation and against the
  zero smallest eigenvalue at the ceiling) are useful.
- Code patches (Appendix B) and a reference specification.
- An exploratory comparison of 14 analyses under prognostic scenarios.

Competing or adjacent literature that the paper does not engage with:

- Percha et al. (JMIR 2019) built a stochastic simulator for N-of-1
  studies with exponential wash-in and wash-out, baseline drift and
  design recommendations; it is the closest prior work on carryover
  modeling in N-of-1 simulation and is not cited.
- Gärtner et al. (BMC Medical Research Methodology 2023) simulate
  aggregated N-of-1 data with carryover and covariate-dependent effects
  from DAG-specified mechanisms and compare Bayesian networks,
  G-estimation and linear models. Directly relevant to DGP choice and to
  effect-modifier treatment.
- Liao and colleagues (arXiv 2112.13991; hierarchical extension in a
  dissertation) on Bayesian distributed-lag models with AR errors for
  carryover in N-of-1 series.
- Classic and current work on power for treatment-by-covariate
  interactions (Brookes et al. 2004; Sargent et al. 2005; Simon 2012;
  Ballman 2015) and the standard result that interaction tests need
  roughly four times the sample of a main effect. The paper never
  relates its power curves to this benchmark.
- Higham (2002) on nearest correlation matrices, relevant to RH's silent
  repair and to the choice of stopping as the alternative.

Assessment. For Trials, which publishes methodology of direct use to
trial design, the contribution is a plausible fit as a 'simulation
methods' article if reframed for trialists. For a top-tier statistics
journal (Statistics in Medicine, Biometrics) it would be below the bar:
the mathematics is elementary, the object of study is one unpublished-in-
practice simulator, and the model is Gaussian with a single biomarker.
Contemporariness: the focus on a 2020 framework is narrow, and the
post-2019 N-of-1 simulation literature (above) is absent.

## 4. Major comments

Correctness and consistency problems first, then framing.

**1. Fit for Trials and audience (framing; INSPECTED).** Trials
Methodology articles are expected to advance methods that trialists use.
Here the stated audience is 'simulation studies that build on RH', and
the recommendations (Discussion, items 1 to 6) are addressed to
simulation developers. A trialist designing an aggregated N-of-1 study
cannot take from the paper a sample size, a defensible c_bm, or a rule
for choosing a design. Specific gaps: (a) c_bm is a correlation whose
clinical meaning is not translated; the paper concedes that the slope
(-0.130 outcome points per mmHg at c_bm = 0.25, VERIFIED) cannot be
compared with the 14 CAPS points per 10 mmHg of Raskind et al.; (b) no
worked example shows the number of participants needed for 80% power in
Hybrid, crossover and OL+BDC designs under the refined process; the
closed-form power at N = 70 is the only operating point and never
exceeds 0.72 (E9, compound symmetry) or 0.34 (separable AR(1), c_bm =
0.25), so the reader learns that the designs are underpowered at
c_bm = 0.25 without being told what is adequate; (c) no software
route is described for a trialist (the package function is noted in
Appendix B.5, but the measurement-error option is not in the package).
Remedy: add a section 'Using the refined process to plan a study' with
a sample-size table by design, c_bm (with implied slope) and half-life,
and relate the results to the standard four-fold interaction penalty;
state in the Background what a trialist should do differently as a
result.

**2. The reference specification is not the object that was analyzed
or implemented (INSPECTED, partly VERIFIED).**

- Appendix B.5 states that `buildSigma()` builds the cross-component
  correlation from `c.cf1t` at the same visit and `c.cfct * rho^tg` at
  other visits, not `c_1 rho^|w_t - w_s|` at every lag as Table 'spec'
  and Section 'Why separability' state. The Kronecker factorization,
  the ceiling formula and the validity-for-every-schedule claim hold
  only for the latter. The authors leave this 'not resolved'. I verified
  that the separable form gives the quoted ceilings; I did not run the
  package, so I cannot say which numbers the package would give. Either
  the package or the specification is wrong, and the paper's central
  deliverable depends on which.
- The stress-test tables come from a first run with the published
  recursion and fixed allocation. The rerun used anchored carryover and
  randomized path allocation, and randomization has since been removed
  from the process definition (Results, last section; Table 'spec').
  No reported stress-test number therefore comes from the process as
  now defined. The text says verdicts 'did not change within Monte
  Carlo error' but then reports that two analyses changed screen status
  and the leading model changed by 0.026 (0.214 against 0.240).
- The 'Estimand' section states that closed-form results and RH's
  analysis use the binary indicator while the stress test uses D_bc at
  the true half-life; the F1 quick check uses a fixed assumed
  half-life of 0.2. Three analyses, three estimands, one conclusion.
- Measurement error (a stress-test default: SD 3 and 6) and dropout are
  in Table 'spec' as options but exist only in the drivers.

Why it matters: a methods paper proposing a reference DGP must make the
specification, the code and the reported numbers one object. Remedy:
settle the cross-component question; rerun the stress test (28,000
trials is inexpensive at the stated 2-core host) under the final
definition, with fixed allocation, and replace the first-run tables;
tabulate which estimand each analysis targets.

**3. At the pharmacological half-lives the refined process is
effectively no carryover, so the 'restoration' is the no-carryover
result (VERIFIED computations; judgment INSPECTED).** Visits are at
least one week after discontinuation, so mean off-drug coupling is at
most 2^(-10) at t_half = 0.1 and the coupling gap is 0.991 at 0.2 week
(my recomputation: pooled slope ratios 0.9906, 0.906, 0.754 at 0.2,
0.5 and 1 week, matching the text). The abstract's 'power raised from
0.13 to 0.64' is thus a return to 0.67, the no-carryover value
(Table 'patch'). That is a valid and useful correction of RH, but it
means (a) the default scenario of the refined process (0.1 to 0.2 week)
cannot show any carryover effect, so the process is not a test bed for
carryover until half-lives of 0.5 week or more, for which the paper
offers no pharmacological support and calls 'behavioral'; (b) the
claim that 'power declines smoothly' is a statement about the
non-default scenarios. Further, the 0.1 to 0.2 week range is argued
from a plasma half-life of 2.3 hours (Hobbs et al. 1978, plausible) by
assuming 7 to 15 plasma half-lives; the text says concentration 'falls
by orders of magnitude within a day', which if anything argues for
shorter symptom half-lives. The range is an assumption and should be
labeled so. A second conceptual point: graded coupling makes the
biomarker association fade as the drug response fades off drug, but on
drug the moderation is the full c_bm sigma_BR b at every visit, even
during titration when the Gompertz mean is near zero (mean at 1 week
on drug is a small fraction of the maximum). The pharmacological
justification given ('it predicts the size of the drug effect, so it
cannot predict a drug effect that is no longer present') applies
equally to the early on-drug weeks. The Discussion's proposed
random-sensitivity model acknowledges this implicitly. Remedy: state
the asymmetry as a limitation, show power under a coupling proportional
to the current mean, and reword the abstract to say the correction
returns the no-carryover result at pharmacological half-lives.

**4. Propositions and scope of the closed form (VERIFIED algebra;
scope INSPECTED).**

- Proposition 3 gives strict monotonicity for the path-stratified
  statistic. The pooled statistic's residual variance also contains the
  between-path variation of mean contrasts, which depends on t_half
  through the carryover means, so monotonicity of the pooled power is
  only checked on a grid (as the text says). Strictness also fails for
  any path whose off-drug visits all precede first exposure (the
  placebo-first crossover path), where g_off is fixed at 0. The wording
  'strictly decreasing' in Appendix A.3 and the abstract's 'monotone'
  should be qualified.
- Proposition 1 holds for unrepaired matrices only, as stated; under the
  refined process (stop on non-PD) this is satisfied, but the
  proposition is stated for compound symmetry while Table 'spec' uses
  separable AR(1) with expectancy-scaled PB variance; the variance
  formula in A.1 is given only for compound symmetry. A one-line
  extension to the separable case is needed since the recommended
  process is separable.
- Closed-form power is exact for E9, not for the mixed model. The
  agreement for the Hybrid design is good (E9b 0.654 against lmer
  0.632 at t_half = 0, 500 reps, MCSE about 0.02) but poor in the
  crossover design (step rule, 0.279 against 0.306 to 0.354,
  `e9-vs-lmer.csv`). The Limitations mention this but the abstract
  does not.
- The statement that the reconstruction and RH-code results differ by
  'about 2.4 combined MCSE, within the spread of the reconstruction's
  positive-half-life cells' (0.08 to 0.11) is not accurate: 0.126 lies
  outside that range. A 2.4 SE discrepancy is nominally significant and
  should be discussed (analysis-model difference, E9b versus lmer, is
  the likely reason; SUSPECTED).

**5. Attribution and numerical claims about RH (INSPECTED).**

- The abstract states that 'the published effect sizes of 0.3 and 0.6
  relied on silently repaired matrices'. The Background states that at
  0.3 the repair changed no correlation by more than 0.007, so results
  are 'essentially those of the stated process'. The abstract should say
  that only 0.6 is materially affected.
- 'RH reported that power ... fell from 0.74 to 0.12 at a half-life of
  0.1 week' (Background). The file `published-power.csv` behind these
  numbers is a 1,000-replicate rerun by the present authors of RH's
  code, not a figure read off RH's publication. State this, and say
  whether RH's printed values agree.
- Table 'nonpd' reports the refined process failing PD at 26% of the
  published grid (21 of 81) against 25% for the published process, and
  75% at c_bm = 0.6. The text frames the published pass at 0.6 under
  carryover as 'cells that have lost their signal'; this is correct
  (VERIFIED by the constant-coupling argument), but the table also
  shows that the refined process is more restrictive in the grid RH
  actually used. The caption specifies rho = 0.8, whereas Table 'spec'
  recommends rho = 0.7; give the 0.7 counts too.

**6. The stress test is post hoc, partly circular and its headline
selection is not robust (INSPECTED).**

- The study is declared exploratory (to its credit), but then yields a
  recommendation ('selects the expectancy-term model') from a rule
  that was revised after seeing results, with the analysis set growing
  from 11 to 14 and the 'reference' analysis replaced. The selected
  model and its two phase-variance variants tie within 0.005 on the
  chosen metric (0.216, 0.221, 0.221), so the 'selection' is a
  tie-break on parameter count.
- The prognostic scenarios are constructed so that the biomarker is
  associated with TV or PB but not with baseline (BL is independent of
  all other variables in the DGP). That is the exact mechanism by which
  analyses sharing the biomarker main effect across baseline and
  follow-up fail (Appendix A.4, which labels its account 'heuristic, not
  a derivation'). The only scenario with a plausible baseline
  association (baseline severity) produced little harm. The abstract
  statement that a prognostic biomarker 'makes analyses that share the
  biomarker effect with baseline anticonservative' should therefore be
  conditioned on a prognostic effect that is absent at baseline; the
  plausibility of that scenario and of the 0.3 magnitude is not argued.
- Robustness costs power: the selected model gives 0.431 against 0.496
  for the standard model at c_bm = 0.30 (13% less pooled; 29% in the
  first run). In an exploratory comparison with MCSE around 0.02 per
  cell this difference, and the opposite sign at compound symmetry
  (0.510 against 0.490), cannot support a general recommendation.
- Preliminary F1 check (Results). Half-lives of 0, 0.1 and 0.2 are
  exactly where theory predicts a flat power curve, so the check is
  uninformative on whether power 'survives' carryover; the reported
  rise (0.262 to 0.294) is about 1.1 combined SE at 500 replicates and
  is not evidence of anything. Either extend to 0.5 to 2 weeks or remove
  the section.
- Scenario counts in the ADEMP table (41 cells, 14 analyses, 28,000
  trials plus post hoc additions) are a lot of researcher freedom for a
  paper whose message is about a DGP. Consider moving the stress test
  to a companion paper (the project already has several) and keeping
  this one to the DGP correction, which is cleaner and more citable.

**7. Not submission-ready as a manuscript (INSPECTED).**

- Six inline editorial notes of the form 'Notation, 2026-10-09: ...' in
  the Methods and Appendix A, 'Updated 2026-10-09 / 2026-10-15' notes in
  Appendix B.3, and an unresolved caveat in B.5. A stray '(lambda is
  the carryover decay rate)' sits after the non-separable PD condition
  in A.5 without a lambda in the display.
- `appendix-b.Rmd` is a different, stale document under the same
  heading: it refers to 'Architecture B', 'covar' and 'orig', Sections
  2.2.2 and 2.2.4 that no longer exist, a re-run 'flagged' in a section
  that does not exist, a September 2026 correction, a visible
  '[TO RE-VERIFY]' block, and line numbers. Its content (correlation
  structure and ceiling discussion) is not the Appendix B ('Code
  changes') of `report.Rmd`. It also claims that the compound-symmetry
  ceiling is tighter 'than AR(1) does', which contradicts the main
  text's finding that AR(1) with decaying cross-correlation lowers the
  ceiling in the open-label design. It should be rewritten or removed.
- Declarations: Competing interests, Funding, Authors' contributions
  and Acknowledgements are placeholders. RH's author list
  (Hendrickson, Thomas, Schork, Raskind) overlaps with a project that
  describes itself as 'pmsimstats team'; if any author of this paper is
  an author of RH, the manuscript's third-person framing ('to our
  knowledge, the first ...') and claims of defects in RH's code should
  be accompanied by disclosure and by a statement that RH's authors
  agree with the diagnosis or were consulted.
- Availability of data and materials points to file paths in a
  repository without URL, DOI or archived release; several load-bearing
  arguments are deferred to `docs/46`, `docs/53`, `docs/56` and
  `compendium-mathematics.tex`, and a key claim (responder-class
  equivalence, c_bm^2 = f_B f_BR) rests solely on an unpublished
  manuscript (paper 03). Trials expects data and code deposited in a
  persistent public archive with an identifier.
- The Trials format: the abstract has the four headings (310 words,
  within the 350-word limit as I recall it; I could not retrieve
  Trials' current guideline page, so the limit and the exact Methodology
  article requirements are UNVERIFIED) but lacks a 'Trial registration:
  not applicable' line, which Trials normally expects; keywords are
  present. The abstract's Results section is largely about the code
  correction and does not give the stress-test result or any figure
  usable by a trialist.

## 5. Minor comments

1. Appendix A.5, 'RH's structure': the positive-definiteness margin
   (1 - rho_cs) - (c_1 - c_x) is 0.10 at RH's values, not 0.20 as
   written. VERIFIED: the smallest eigenvalue of RH's 8-visit response
   block is 0.1000.
2. Same paragraph, 'which fails on dense schedules': for the Hybrid
   schedule at rho = 0.8 the smallest eigenvalue of A is 0.1218, above
   the threshold (c_1 - c_x)/(1 - c_x) = 0.111, so the non-separable
   form is valid there (VERIFIED); give a schedule where it fails (CO at
   2.5-week spacing has 0.281, also valid). The claim is untested for
   the paper's designs.
3. Table 'ceilings': the caption says the minimum is over three designs;
   the 0.474 at rho = 0.7 is attained in OL+BDC (VERIFIED: Hybrid 0.480,
   OL+BDC 0.474, CO 0.616). Say so, and state that CO is never binding.
4. Methods, 'The published data-generating process': the modified
   Gompertz G is introduced with offset term e^{-d}; state that the
   code subtracts the same offset and rescales so that G(0) = 0 and
   G(infinity) = m.
5. 'Number of replicates': 'chance alone produces [flags] in about
   4.6% of cells' is correct for 2 MCSE; the observed 10 of 278
   (3.6%) is consistent with this. Say that the 278 cell-variants are
   not independent.
6. The size screens report 0.057 for Kenward-Roger and 0.057 for the
   AIC expectancy model, both above the 0.046 to 0.054 pooled band, but
   only the second is described as a failure; state the rule uniformly.
7. Table 'prognostic' lacks an MCSE note; with 1,000 replicates the
   single-cell band is 0.036 to 0.064, so 0.063 and 0.062 for 'free
   baseline slope' are at the boundary and are called 'nominal'.
8. Results (Room for larger effect sizes): 'A biomarker correlation of
   0.45 is feasible in every design under the separable form at rho =
   0.7' holds without carryover (VERIFIED, minimum 0.474); state that
   the check was at t_half = 0.
9. Several items of notation: Methods defines the coupling gap as
   1 - g_off, NOTATION.md as g_on - g_off (equal here since g_on = 1);
   NOTATION.md describes c_x multiplied by rho^lag in the separable
   form, whereas the separable form in this paper has no c_x. Align.
10. Figure captions: 'Hybrid design, N = 70, compound symmetry' on
    Figures 1 and 3 should also state E9 or E9b and the number of Monte
    Carlo replicates (5,000). Figure fig4bmc has 250 replicates (MCSE up
    to 0.032), too few to read differences of 0.05 to 0.10 between
    half-lives; the quoted numbers (e.g. 0.08 against Table 'patch'
    0.126 at 0.1 week for the same design and seeds) illustrate this.
11. Table 'refinements', step 3: anchored carryover is a correction
    only if the compounding in the published recursion is unintended;
    the paper gives no evidence of RH's intent. Describe it as a
    different specification that the authors judge more coherent, and
    note its effect is at most 0.006 in power.
12. 'Re-exposure' is retained as in RH (time on drug restarts). The
    paper itself notes that, combined with anchored carryover, this
    produces a downward step at re-exposure; show the mean trajectory
    for one path so readers can judge plausibility.
13. Exploratory status paragraph is clear and should be kept; however
    the text elsewhere uses 'recommended analysis' and 'selected'.
    Use 'preferred in this comparison'.
14. Abbreviations list omits E9, E9b, ADEMP, MCSE, OL, BDC (OL+BDC is
    there), CAPS.
15. References: the Statistics in Medicine CSL is used for a BMC
    journal; reformat to Trials' numbered style at submission.
16. Language: 'Notation' notes and 'to our knowledge, the first'
    claims (twice in the manuscript) need a citation search statement
    (see Section 8) or removal.

## 6. Missing and questionable references

| Location | Issue | Suggested source |
|---|---|---|
| Background, 'Closed-form power ... only simplified' | Omits interaction-power literature | Brookes ST et al. J Clin Epidemiol 2004;57:229-236; Ballman KV. J Clin Oncol 2015;33:3968-3971 |
| Background, DGP choice | No prior N-of-1 simulators discussed | Percha B et al. J Med Internet Res 2019;21:e12641 (corrected 2020, 9:e16179) |
| Background, carryover | Recent aggregated N-of-1 simulation with carryover and effect modification | Gärtner et al. BMC Med Res Methodol 2023 (PMC10440905; details to confirm) |
| Pharmacological grounding | Carryover modeling alternatives | Liao et al. arXiv:2112.13991 (Bayesian distributed lag, AR errors) |
| Priority claim ('first in N-of-1 literature') | Unsupported search statement | Add search protocol; also Zucker 2010, Duan 2013 already in bib but not discussed here |
| Predictive vs prognostic | Definitions and design context | Sargent DJ et al. J Clin Oncol 2005;23:2020-2027; Simon R. Stat Med 2012;31 (predictive medicine trials) |
| Silent repair | Alternative repair methods | Higham NJ. IMA J Numer Anal 2002;22:329-343 |
| ADEMP | Cited (Morris 2019) but the template is richer | Siepe BS et al. Psychol Methods 2024 (simulation-study template; check details) |
| Discussion, c_bm^2 = f_B f_BR | Rests on unpublished paper 03 | Publish or give the derivation in an appendix |
| Pharmacological grounding | Plasma t1/2 does not give symptom t1/2 | Add PK/PD source on alpha-1 receptor occupancy or state as assumption |
| Baseline as covariate or outcome | Cited (Liu 2009) | Also Vickers and Altman BMJ 2001 (ANCOVA) and Kenward-Roger-based longitudinal analyses |
| Satterthwaite, CR2 | Correct as cited | Kuznetsova 2017; Bell and McCaffrey 2002; Pustejovsky and Tipton 2018 (INSPECTED, plausible) |

Reference details marked 'check' or given without article numbers were
recalled from search results and should be verified before use.
Citations present in the bibliography (Hobbs 1978, Raskind 2016, Wang
2019, Morris 2019, Liu 2009) plausibly support the claims attached to
them (INSPECTED against titles only; I did not read the full texts).

## 7. Suggestions for strengthening the paper

- Split the paper: (a) the DGP correction (Sections on the published
  process, closed form, patches, ceiling) as the Trials methodology
  article; (b) the stress test and analysis recommendation as a separate
  paper with a fresh, fixed-protocol confirmation run.
- Add a trialist-facing section: implied slope in outcome units, power
  tables by design, a decision flow for choosing c_bm below the ceiling,
  and how to elicit t_half and rho.
- State in the abstract that the correction returns the no-carryover
  result at pharmacological half-lives, and give the stress-test result
  with its exploratory label.
- Provide a package vignette that reproduces Tables 'closedform',
  'patch' and 'ceilings' in a few lines, deposit a tagged release with a
  DOI (Zenodo), and cite it.
- Add a figure of the correlation matrix before and after the repair
  at c_bm = 0.6, and one of mean trajectories by path, to make the
  pathologies visible.
- Test the closed form beyond E9 with a simulation under the refined
  process using the mixed model of RH across all three designs.
- Consider reporting power as a function of the implied slope rather
  than of c_bm, so that the compound-symmetry and separable results are
  commensurate (the paper's own advice in the Ceiling section).

## 8. Scope of this review

Verified (executed or re-derived):

- Propositions 1 and 2 and the noncentral-t power formula (A.1, A.2),
  by reading the closed-form script and reproducing every entry of
  Table 'closedform' (step and graded, E9 and E9b, six half-lives) with
  the authors' functions run from a scratch copy.
- Coupling-gap values 0.991, 0.906, 0.754; c_bm = 0.2 E9b values; CS
  ceiling 0.2555 (Hybrid), 0.256 (CO), 0.262 (OL+BDC); graded ceilings
  0.307 (1 week), 0.369 (2 weeks); AR(1) decaying 0.291; separable 0.440
  (rho 0.8), 0.480 and 0.474 (rho 0.7).
- The separable ceiling factorization, using my own independent
  implementation, including a smallest eigenvalue of zero at the ceiling.
- The cost table of A.5 and the non-separable PD condition (derived by
  hand; the eigenvalue formulas are correct).
- RH's compound-symmetry margin (0.10, not 0.20); non-separable
  validity for the Hybrid and CO schedules at rho = 0.8.
- Table 'patch' entries against `graded-patch-reps500.csv`; the
  N-of-1 0.74/0.12/0.10 values against `published-power.csv`.

Inspected only (read, not executed): all prose, tables of stress-test
rejection rates, the reported 28,000-trial results and the
post hoc pooling (stored outputs not recomputed), the `lmer`-versus-E9
comparison file, Appendix B code listings (not run against commit
58b32a9), the package `buildSigma()` (not read or run), NOTATION.md
for consistency, `appendix-b.Rmd`.

Not checked: the 20-of-54 and 21-of-81 repair counts (the count files
exist and the totals are internally consistent); the 0.007 and 0.10
maximum correlation changes; Monte Carlo replicates; the F1 quick
check; figure files themselves; the earlier referee report of
2026-06-13 (not read, as instructed it is out of date); the full texts
of the cited papers.

Searches run (web search only; no email or identifiers passed to any
service): Trials methodology-article submission guidelines (the
journal's own page redirected and could not be read, so format
requirements are partly UNVERIFIED); aggregated N-of-1 power with
biomarker interaction and carryover; Bayesian hierarchical N-of-1
models with carryover and covariate interaction; predictive versus
prognostic biomarker interaction power; Percha et al. N-of-1 simulator.
I did not search for a prior statement of the step-coupling defect in
other work, nor for post-2020 critiques of RH specifically; the
'first in the literature' claims are therefore neither confirmed nor
refuted here. PubMed tools were not used.
