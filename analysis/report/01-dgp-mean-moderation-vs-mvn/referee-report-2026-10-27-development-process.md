# Second referee report: paper 01 as a development process toward a DGP
*2026-10-27 11:48 PDT*

Manuscript: 'Simulating predictive-biomarker interactions in aggregated
N-of-1 trials: a refined data-generating process following Hendrickson
et al.' (`report.Rmd` with Appendices A and B; companion
`appendix-b.Rmd`). Target journal: Trials (Methodology article). This
is an independent second review. The first report of today
(`referee-report-2026-10-27.md`) was read before the review began and
is engaged with in Section 4; judgments here were formed and checked
separately. The frame for this review is the author's stated intent:
that the paper present a logical, step-by-step development from the
published process of Hendrickson et al. (2020; 'RH', commit 58b32a9)
to a recommended (the author's word: optimal) data-generating process
for simulating aggregated N-of-1 trials of the RH type, in particular
the Hybrid design. Mean moderation was removed from this paper on
2026-10-27 by design and is not recommended back in. The studies are
treated as exploratory throughout.

Evidence labels: VERIFIED (re-derived or executed by me), INSPECTED
(read and judged), SUSPECTED (plausible, not confirmed). Anything not
labeled is a judgment about exposition.

## 1. Summary of the manuscript

The paper takes RH's simulation framework for aggregated N-of-1 trials
of prazosin in PTSD, in which the biomarker-by-treatment interaction is
a correlation c_bm between a baseline biomarker and the drug-response
component of a multivariate normal symptom model, and identifies two
defects in the published code: a step coupling rule that keeps the full
biomarker correlation at every off-drug visit once any drug mean
remains, and a silent repair of correlation matrices that are not
positive definite. A paired-difference statistic (E9, E9b) gives
closed-form moments and noncentral-t power (Propositions 1 to 3), from
which the step rule's removal of the interaction at every positive
half-life follows, and a Schur-complement bound gives the largest
feasible c_bm (0.256 under RH's compound symmetry; 0.440 to 0.480
under a separable AR(1) structure). Graded coupling is shown to restore
a smooth power decline, in closed form and in a patched rerun of RH's
code. A reference specification (Table 'spec') then adds AR(1) within
components at rho = 0.7, a separable cross-component block with c_1 =
0.2, anchored carryover, a stop on non-PD matrices, retention of RH's
fixed path allocation and re-exposure rule, and options for prognostic
couplings, measurement error, decay form and dropout. An exploratory
stress test of 14 analyses under that process (41 cells, 28,000 trials
plus post hoc additions) finds that analyses sharing the biomarker main
effect with baseline are anticonservative under a prognostic biomarker
and prefers an expectancy-term model.

## 2. Overall assessment and recommendation

Recommendation: major revision. The technical core is correct (I
re-derived the closed-form power, the coupling-gap values, every entry
of Table 'ceilings', and the patch results against the stored files;
Section 9), and the two defects in RH's code are real and
consequential. Read against the author's intent, however, the
manuscript does not yet read as a development process. It presents the
decisions (Table 'refinements', Table 'spec', 'Why separability') in
the Methods before any evidence for them appears in the Results; it
never states what makes one DGP better than another, so 'refined' has
no criterion behind it; two of the six elements are justified by
evidence, one by a derivation, and the rest by assertion or
preference; and the rho = 0.7 default, which is never argued for,
accounts for most of the headline rise in the feasible effect size
that the paper attributes to separability (VERIFIED, Section 5,
comment 2). The recommended process also exists in three slightly
different forms (specification, stress-test driver, package). For
Trials the paper is a plausible Methodology article only after it is
reorganized around explicit criteria and stages, with a short
planning section that tells a trialist what to simulate and what
remains uncalibrated. The first report's recommendation (major
revision; not suitable in present form) is one I share, for partly
different reasons.

## 3. The development process: ideal chain against the manuscript

The ideal chain for each element is: defect or limitation in the
published process; diagnosis (why it happens); candidate remedy;
evidence that the remedy works and does not cost more than it gains;
decision, with the reason stated. The table below maps the manuscript
onto that chain for every element of Table 'spec'.

| Element | Defect stated | Diagnosis | Remedy | Evidence | Decision reason | Verdict |
|---|---|---|---|---|---|---|
| Graded coupling | Yes: power collapse under any carryover (Background 2) | Yes, exact: Prop. 2, gap = 0 under step rule (Results 1) | Yes: g = 2^(-tsd/t1/2) off drug | Strong: closed form, 5,000-rep MC, patched RH code (Table 'patch', Fig. fig4bmc) | Pharmacological reading given (Methods, 'Coupling') | Complete. The one well-developed chain. Alternatives to the specific profile (coupling proportional to the current mean; a separate coupling half-life) are not considered, and at the default t1/2 the remedy is numerically indistinguishable from zero off-drug coupling (gap 0.9997 at 0.1 week, VERIFIED) |
| Stop on non-PD | Yes: 20 of 54 matrices repaired; correlations moved up to 0.10 at 0.6 (Background 1) | Yes: Schur ceiling (Methods, A.4) | Yes: stop (B.2) | Partial: repair counts and ceiling; no evidence of what the repair did to RH's reported power at 0.6 | Implicit (a repaired matrix is a different process) | Diagnosed and remedied; consequence for the user left open. Under the stop, RH's 0.6 cells become infeasible (75% of the extended grid, Table 'nonpd') and the paper offers no implemented route to effect sizes above the ceiling (the sensitivity model is 'not yet implemented') |
| Anchored carryover | Asserted: recursion 'compounds' the decay | Yes, mechanical (B.3) | Yes | Effect at most 0.006 in power (B.3) | None: RH's intent is not established | A specification choice presented as a correction. Harmless, but it dilutes the correction/choice distinction the paper itself draws |
| AR(1) within components | Asserted: 'compound symmetry overstates long-lag correlation' (Table 'refinements', step 4) | None: no data, no citation, no diagnostic | rho^|w_t - w_s| | None for the defect; cost shown (power 0.721 to 0.335 at c_bm = 0.25, rho = 0.8) | None given | Asserted. The claim about long-lag correlation is a statement about PTSD symptom data that the paper does not support |
| rho = 0.7 default | Not stated as a defect of RH's 0.8 | None | Lower rho | None; rho = 0.7 appears only in Table 'spec', the ADEMP table and the ceiling numbers | None | Missing step. Lowering rho from 0.8 to 0.7 alone raises the compound-symmetry ceiling from 0.256 to 0.346 and the separable ceiling from 0.440 to 0.474, and lowers E9 power at c_bm = 0.25 from 0.335 to 0.290 (VERIFIED). The paper never separates this change from separability |
| Separable cross-component block, c_1 = 0.2 | Partial: non-separable form can be non-PD on dense schedules; ceiling too low | Yes: A.5 (eigenvalues, factored ceiling, matched residual structure) | c_1 rho^|lag| | Derivation only; the PD failure 'on dense schedules' is not exhibited for any design in the paper (the first report verified it does not fail for Hybrid or CO at rho = 0.8) | Four properties (Methods, 'Why separability') | Justified as a modeling choice by derivation, not by evidence of fit. At rho = 0.7 the move from the package's two-constant form to the separable form raises the minimum ceiling only from 0.452 to 0.474 (VERIFIED). c_1 = 0.2 is RH's same-visit value carried over without comment |
| Fixed path allocation | Not a defect; retained | n/a | n/a | None; the stress-test rerun used random allocation and the definition was changed afterward (Results, last section; driver default `--allocation random`, INSPECTED) | Stated: 'as in RH' | Retained, but the only reported run under the refined process did not use it |
| Re-exposure (clock restarts) | Not a defect; retained | Argued: re-titration of prazosin | Option to rebuild from residual for non-titrated drugs | None; the downward step at re-exposure it produces with anchored carryover is acknowledged (B.3) but not shown | Pharmacological argument | Retained by argument. The mean trajectory under this rule is never plotted |
| Prognostic couplings (option) | Limitation: 'not representable' in RH | n/a | Correlations of the biomarker with TV or PB | Yes: prognostic null cells show which analyses fail (Results, 'Consequences') | Stated: a simulation that omits it cannot reveal the failures | Justified by the paper's own evidence as an option. The 0.3 magnitude and the absence of a baseline association in three of four scenarios are not argued |
| Measurement error, decay form, dropout (options) | Limitation only | n/a | Yes, in drivers | None reported: the stress test varied them, but the paper reports no result about their effect | None | Asserted. They appear in Table 'spec' with stress-test defaults and nowhere else; measurement error is not in the package (B.5) |

Taken together: one element has a complete chain (graded coupling),
one is diagnosed and remedied but its consequence is left open (the
stop), one is justified by derivation (separability), and the
remaining elements are retained, asserted, or absent from the argument
(AR(1), rho = 0.7, c_1, anchored carryover, options). The manuscript's
own distinction between 'corrections' (steps 1 to 3) and 'modeling
choices' (steps 4 to 6) is the right one, but it is drawn in a table
in the Methods and not carried through: anchored carryover sits among
the corrections on no evidence of RH's intent, and rho = 0.7, the
choice with the largest numerical effect on the headline ceiling, is
not in the table at all.

A second structural problem is order. The Methods present the
Refinements table, the four properties of separability and the full
reference specification before the Results show the first defect. A
reader following a development process meets the conclusions first
and the evidence later, and has to reconstruct which result supports
which decision. The estimand section, placed before the closed form,
tells the reader that three analyses target three different
quantities (binary indicator in the closed form and RH's analysis;
D_bc at the true half-life in the stress test; F1 at a fixed 0.2
week), which is true but is not then used to reconcile the three sets
of results.

## 4. Relation to the first report

Where I agree (and have checked independently where marked):

- Major 2 (the specification, the code and the reported numbers are
  not one object). Agreed and partly settled by inspection: the
  stress-test driver (`11-strawman-stress-test.R`, line 184) builds
  the cross-component block as `C1 * A`, which is the separable form
  of Table 'spec', so the reported stress-test numbers come from the
  specified structure; the package `buildSigma()` (`R/generateData.R`,
  lines 389 to 398) uses `c.cf1t` at lag zero and `c.cfct * rho^tg`
  otherwise, which is not. The divergence is between specification and
  package, not between specification and reported results. That
  narrows the repair: fix or document the package. The allocation
  discrepancy (rerun under random allocation; specification fixed)
  stands as the first report describes (INSPECTED, driver default).
- Major 3 (at pharmacological half-lives the refined process is
  effectively no carryover). Agreed; my recomputed gaps are 0.9997,
  0.9906, 0.906, 0.754 and 0.551 at 0.1, 0.2, 0.5, 1 and 2 weeks
  (VERIFIED). The on-drug asymmetry (full coupling during titration
  when the mean is near zero) is a real limitation and belongs in the
  development chain as an alternative that was not tried.
- Major 5 (attribution: `published-power.csv` is a 1,000-replicate
  rerun, not RH's printed figure; the abstract overstates the repair's
  effect at 0.3). Agreed (INSPECTED: the file carries `reps = 1000`).
- Major 7 and the Minor comments on submission readiness, Minor 1 (the
  margin is 0.10, not 0.20; VERIFIED, smallest eigenvalue 0.1000),
  Minor 3 (0.474 is attained in OL+BDC; VERIFIED), Minor 8 (0.45 is
  feasible at t1/2 = 0; VERIFIED), Minor 11 (anchored carryover is a
  choice, not a correction), Minor 12 (plot the re-exposure trajectory).
- The stale companion `appendix-b.Rmd` (INSPECTED: it still refers to
  'Architecture B', `covar` and `orig`, and reproduces a different
  Appendix B).

Where I disagree or weigh differently:

- Major 1 and Suggestion 1 (split the paper; move the stress test
  out). Under the development-process framing the stress test is not
  a detachable study: it is the only evidence in the paper for the
  'prognostic coupling' option and for the analyzability criterion
  (which DGP features change which analysis is valid). I would keep a
  reduced stress test as the evidence for that one element and drop
  the selection of an analysis model, which is the part that is post
  hoc, tie-broken on parameter count and more naturally a companion
  paper. The first report's objection to the selection stands; its
  remedy of removing the whole study would leave the prognostic option
  unsupported.
- Major 4 on Proposition 3: agreed on the qualification, but I would
  not call the pooled-statistic monotonicity a gap in the argument.
  The paper says it was checked on a grid, and the grid check is in
  the stored curves (INSPECTED). The CO placebo-first path, where the
  off-drug coupling is identically zero at every half-life (VERIFIED),
  is the case that needs an explicit exclusion in the statement.
- Major 6, third bullet (robustness costs power, so the recommendation
  is unsupported). The 13% cost is a reason to state the trade-off,
  not a reason to withhold a preference: a DGP paper should say which
  analyses are valid under which features and leave the choice of
  power against robustness to the trialist.

What this report adds:

- The rho = 0.7 confound (Section 5, comment 2), which the first
  report does not raise and which bears directly on the paper's
  central feasibility claim.
- The absence of stated criteria for 'better' (comment 1).
- The order problem (decisions before evidence) and a concrete
  reorganization (Section 8).
- Trialist-facing numbers the paper could supply from its own closed
  form (comment 4): the implied slope is 2.0 symptom points per SD of
  biomarker at c_bm = 0.25 (3.6 at 0.45; VERIFIED), and the E9
  sample size for 80% power under the refined process at t1/2 = 0.2
  week is about 280 at c_bm = 0.25, 140 at 0.35 and 84 at 0.45
  (VERIFIED for E9; the mixed model will differ somewhat).
- A PubMed check of the priority claim (Section 7).

## 5. Major comments

**1. 'Optimal' and 'refined' have no stated criterion (Background,
'Why a follow-up is needed'; Methods, 'Refinements'; INSPECTED).**
Problem: the aims paragraph names three desiderata (consistent with
the pharmacology of short-half-life titrated drugs; intuitive
parameters; flexible features) and the rest of the paper uses at least
three more without naming them: existence of the joint distribution
(the ceiling), identifiability of the interaction under carryover (the
coupling gap), and analyzability (which features change which
analysis is valid). Calibration to trial data is deferred in every
section where it arises. Why it matters: without criteria the reader
cannot judge whether a step is an improvement; the AR(1) step, for
example, is an improvement in feasibility and a loss in power at
matched c_bm, and the paper says both without saying which criterion
governs. Remedy: state the criteria once, in the Background, in a
numbered list (fidelity to the pharmacology; existence of the
distribution for clinically plausible effects; identifiability of the
interaction; matched analyzability; calibration status), and end each
stage of the Results with a sentence stating which criterion the
decision rests on. Replace 'refined' by 'recommended' or keep
'refined' but define it as 'meeting criteria 1 to 4; criterion 5
open'.

**2. The rho = 0.7 default is never justified and is confounded with
the separability argument (Table 'spec'; Results, 'Room for larger
effect sizes'; Table 'ceilings'; VERIFIED).** Problem: the text
attributes the rise in the feasible c_bm to the separable form ('The
ceiling rises to 0.440 at RH's rho = 0.8 and to 0.474 at rho = 0.7').
Table 'ceilings' itself shows, in its compound-symmetry column, that
RH's own structure at rho = 0.7 has a ceiling of 0.346, so RH's effect
size 0.3 becomes feasible by lowering rho alone, with no change of
structure. My recomputation across rho: compound symmetry 0.256, 0.346,
0.410, 0.458 at rho = 0.8, 0.7, 0.6, 0.5; separable 0.440, 0.474,
0.462, 0.448 (minimum over designs). The separable form's advantage is
large at rho = 0.8 and nearly gone by rho = 0.6, so the paper's
preferred combination (separable, rho = 0.7) gets most of its room
from the rho change, which no section argues for. Lowering rho also
costs power: E9 at c_bm = 0.25, Hybrid, N = 70 is 0.721 (CS, 0.8),
0.335 (separable, 0.8), 0.290 (separable, 0.7); the paper quotes the
first two and not the value at its own default. Why it matters: a
development process must show each change separately with its reason;
here a parameter change is folded into a structural change and
credited to the structure. Remedy: add a stage 'Serial correlation'
that states where 0.7 comes from (an estimate from PTSD symptom series,
a sensitivity range, or a convention), shows the ceiling and the
power at matched implied slope for rho in {0.5, 0.6, 0.7, 0.8} under
both structures, and only then introduces separability as a separate
step with its own, smaller, gain and its validity and matched-model
properties as the main reasons.

**3. Decisions are presented before evidence, so the paper reads as a
specification with a justification appended, not as a development
(Methods, 'Refinements', 'Why separability', 'Reference specification';
INSPECTED).** Problem: by the end of the Methods the reader has the
full Table 'spec', the four properties of separability and the
pharmacological grounding, but has not yet seen that the step rule
removes the interaction, that the repair changes correlations, or what
the ceiling is under each structure. The Results then supply the
evidence out of the order in which the decisions were made. Why it
matters: this is the author's stated intent and the present order
defeats it; a referee cannot tell which result licenses which row of
Table 'spec'. Remedy: see Section 8. In brief, keep the Methods for
tools (the closed-form statistic, the ceiling, the simulation designs)
and move every decision into a staged Results section where each
stage ends with the decision and its reason; Table 'spec' becomes the
last table of the Results, with a column citing the stage that
justifies each row.

**4. A trialist cannot use the paper to plan a study (Discussion,
'Recommendations'; INSPECTED; numbers VERIFIED).** Problem: the
recommendations are addressed to simulation developers. Nothing tells a
reader planning an RH-type Hybrid trial what c_bm means in outcome
units, what N is needed, which parameters are RH's values, which are
the authors' choices, and which are uncalibrated. Why it matters: Trials
publishes methodology for trialists; and the paper's own closed form
can supply the numbers in a few lines. Remedy: add a short section
'Using the refined process to plan an aggregated N-of-1 trial' with (a)
the implied slope -c_bm sigma_BR per SD of biomarker (2.0 points at
0.25, 3.6 at 0.45, with sigma_BR = 8), (b) a table of N for 80% power
by design and c_bm under the refined process at t1/2 = 0.2 week (my
E9 values for Hybrid: about 280, 140 and 84 at c_bm = 0.25, 0.35 and
0.45; the authors should give mixed-model values), (c) a list of
parameters by status (RH's; chosen here; uncalibrated), and (d) the
one-line rule for choosing c_bm below the ceiling of every design in
the study. State that the interaction test needs roughly four times
the sample of a main-effect test of the same size and relate the
table to that benchmark.

**5. The recommended process exists in three forms (Table 'spec';
Appendix B.5; `11-strawman-stress-test.R`; `R/generateData.R`;
INSPECTED).** Problem: the specification says c_1 rho^|lag| at every
lag; the stress-test driver implements that; the package implements two
constants (`c.cf1t`, `c.cfct`). The specification says fixed
allocation; the stress-test rerun used random allocation. The
specification lists measurement error; the package lacks it. Why it
matters: a reader who installs the package to reproduce or extend the
paper gets a different process from the one specified and run. Remedy:
make the package match Table 'spec' (or add a `separable = TRUE`
argument that does) and state in B.5 which form each reported table
used; rerun the stress-test cells that feed the paper's tables under
fixed allocation, or change Table 'spec' back to random allocation and
say why; move measurement error into the package or drop it from Table
'spec'.

**6. The stop on non-PD matrices removes RH's larger effect size without
offering a replacement (Results, Table 'nonpd'; Discussion;
INSPECTED).** Problem: under the stop, c_bm = 0.6 is infeasible in
every design in the Hybrid family (75% of the extended grid at rho =
0.8; every Hybrid cell at rho = 0.7 by Table 'ceilings'). The paper's
answer, the sensitivity model, is 'not yet implemented'. A development
process that ends with 'the published effect sizes cannot be
simulated' needs to say what a user should do now. Why it matters: the
clinical association the project cites (14 CAPS points per 10 mmHg) is
large, and a process that caps the effect at a correlation of 0.47
without translating that cap into outcome units leaves the user unable
to tell whether the cap binds in practice. Remedy: translate the
ceiling into an implied slope per SD of biomarker at the recommended
parameters, compare it with the observed association on a common
scale as far as possible (or state exactly why it is not possible), and
either implement the sensitivity model or state the covariance
encoding's range of use as a limitation in the abstract.

**7. The stress test's analysis selection does not belong in a DGP
paper, but its prognostic result does (Results, 'Consequences for the
analysis'; Methods, 'Exploratory status'; INSPECTED).** Problem: the
selection rule was revised after results, the analysis set grew from
11 to 14, the reference analysis was replaced, and the final choice is
a tie broken on parameter count. The paper is candid about all of
this. Why it matters: the selection is the part of the stress test
that attracts the first report's Major 6 and that a reader will
distrust; the prognostic null result (analyses sharing the biomarker
effect with baseline reject at 0.10 to 0.51 under a prognostic
biomarker) is the part that supports a DGP element and is robust in
direction across the first run and the rerun. Remedy: keep the
prognostic table and the size comparison as the evidence for the
prognostic option, state them as 'which analyses are valid under this
feature', and move the selection, the post hoc 1,500-trial pooling and
the shortfall metric to the companion paper. The F1 quick check
(Results, 'A preliminary check') covers half-lives where the closed
form predicts no change and should be removed or extended to 0.5 to 2
weeks.

**8. Several claims in the chain are asserted with no evidence or
citation (Table 'refinements' step 4; Methods, 'Pharmacological
grounding'; INSPECTED).** 'Compound symmetry overstates long-lag
correlation' is a claim about PTSD symptom series with no source. The
0.1 to 0.2 week symptom half-life is derived from a plasma half-life
by an unexplained factor of 7 to 15 (the first report makes the same
point). 'Which fails on dense schedules' (non-separable form) is not
exhibited for any schedule in the paper. Remedy: cite or show data for
the first; label the second an assumption and give the sensitivity
range; for the third, either show a schedule where the condition
fails or delete the phrase.

## 6. Minor comments

1. Abstract, Results: 'raised power in the published code from 0.13 to
   0.64' is correct against `graded-patch-reps500.csv` (0.126 to 0.638;
   VERIFIED), but the no-carryover value is 0.668, so say 'restored to
   within Monte Carlo error of the no-carryover value'.
2. Abstract, Results: 'the published effect sizes of 0.3 and 0.6 relied
   on silently repaired matrices' overstates the case at 0.3 (maximum
   change 0.007, Background); say 'at 0.6, and to a negligible degree
   at 0.3'.
3. Methods, 'The published data-generating process': the Gompertz
   offset is subtracted and rescaled in the code so that G(0) = 0 and
   G(infinity) = m (INSPECTED in `04-power-simulation.R`); the display
   should match.
4. Methods, 'Closed-form analysis': E9b 'reproduces [RH's mixed model]
   mean estimates to within 0.003 and its power closely' holds for the
   Hybrid design (0.654 against 0.632 at t1/2 = 0); in the crossover
   design the closed form is 0.645 against 0.708 and 0.279 against
   0.306 to 0.354 under the step rule (`e9-vs-lmer.csv`, INSPECTED).
   State the design.
5. Table 'ceilings' caption: name the design attaining each minimum
   (OL+BDC at rho = 0.7 separable, 0.474; Hybrid at rho = 0.8 CS,
   0.256; VERIFIED) and state that the crossover design is never
   binding.
6. Table 'nonpd' is at rho = 0.8; Table 'spec' recommends 0.7. Give the
   counts at 0.7 or explain why 0.8 is shown.
7. Appendix A.5, 'RH's structure': the margin at RH's values is 0.10,
   not 0.20 (VERIFIED; same as the first report's Minor 1).
8. Appendix A.3, Proposition 3: exclude paths whose off-drug visits all
   precede first exposure (CO placebo-first; VERIFIED that g_off is
   identically 0 there), where the slope is constant in t1/2.
9. Appendix A.5, non-separable form: delete the stray '(lambda is the
   carryover decay rate)'.
10. Results, 'Graded coupling restores a smooth decline': the sample
    size scaling '1/G^2' is a small-r approximation; say so and give
    the exact statement (noncentrality proportional to r/sqrt(1-r^2)).
11. Table 'refinements': move rho = 0.7 into the table as its own step
    (see Major 2), and move anchored carryover out of the 'corrections'
    group (see Section 3).
12. The six 'Notation, 2026-10-09' notes and the 'Updated 2026-10-09 /
    2026-10-15' notes in B.3 are revision history, not manuscript text.
13. Declarations are placeholders; 'Availability of data and materials'
    cites repository paths with no URL, DOI or release; 'Competing
    interests' should state whether any author is an author of RH.
14. `appendix-b.Rmd` is a different, stale document and should be
    removed or rewritten to match Appendix B of `report.Rmd`.
15. Abbreviations: add E9, E9b, ADEMP, MCSE, OL, BDC, CAPS, PD.
16. Keywords and abstract headings are in Trials' form; a 'Trial
    registration: not applicable' line is customary for BMC journals
    (UNVERIFIED for Trials specifically; the journal's guideline page
    redirects to a login and could not be read).
17. The CSL is Statistics in Medicine; Trials uses numbered Vancouver
    style.

## 7. Missing and questionable references

| Location | Issue | Source (status) |
|---|---|---|
| Background, 'Closed-form power calculations exist only for simplified versions' | Prior simulators of N-of-1 designs with explicit wash-in and wash-out, not discussed | Percha B et al. J Med Internet Res 2019;21(4):e12641, doi:10.2196/12641 (PubMed metadata retrieved; relevance INSPECTED from abstract) |
| Background and Discussion, DGP with carryover and covariate-dependent effects | Closest recent aggregated N-of-1 simulation with carryover and treatment dependence on covariates | Gärtner T et al. BMC Med Res Methodol 2023;23:191, doi:10.1186/s12874-023-02012-5 (PubMed metadata retrieved) |
| Table 'refinements', step 4 (CS against AR(1)) | A prior N-of-1 simulation that generated under both CS and AR(1) with 20% carryover and compared analyses | Chen X, Chen P. PLoS One 2014;9:e87752, doi:10.1371/journal.pone.0087752 (PubMed metadata retrieved) |
| Methods, serial correlation | Serial correlation and washout in personalized studies, by an RH co-author | Schork NJ. Harv Data Sci Rev 2022, doi:10.1162/99608f92.f1eef6f4 (PubMed metadata retrieved) |
| Background, carryover in aggregated N-of-1 simulation | Already in the bibliography (`blackston2019`) but not cited in this paper's text at the point where carryover's effect on Type I error is relevant | Blackston JW et al. Healthcare 2019;7:137, doi:10.3390/healthcare7040137 |
| Background, 'to our knowledge, the first' (twice) | A PubMed query for N-of-1 AND simulation AND biomarker/effect modifier AND power returned only RH itself (one record, 2026-10-27). This weakly supports the claim; add a search statement | PubMed search, see Section 9 |
| Comment 4, interaction sample size | The four-fold penalty for interaction tests and predictive-biomarker design | Brookes ST et al. J Clin Epidemiol 2004;57:229-236; Sargent DJ et al. J Clin Oncol 2005;23:2020-2027 (recalled; not retrieved today) |
| Appendix B.2, repair | Nearest-correlation-matrix repair, the method RH's `make.positive.definite` approximates | Higham NJ. IMA J Numer Anal 2002;22:329-343 (recalled; not retrieved today) |
| Discussion, c_bm^2 = f_B f_BR | Rests solely on an unpublished companion (`pmsimstats-paper03`) | Give the two-line derivation in Appendix A |
| Methods, 'Pharmacological grounding' | Hobbs 1978 supports the plasma half-life only; the symptom half-life range is an assumption | Label as assumption or add a PK/PD source |
| Methods, simulation reporting | Morris et al. 2019 is cited correctly (INSPECTED, bibliography entry matches the Stat Med article) | No change |
| Analysis citations (Kuznetsova 2017; Bell and McCaffrey 2002; Pustejovsky and Tipton 2018; Kenward and Roger 1997; Liu et al. 2009) | Plausible for the claims attached (INSPECTED against bibliography entries, full texts not read) | No change |

## 8. A proposed restructuring

The outline below keeps the Trials Methodology skeleton (Background,
Methods, Results, Discussion, Conclusions) and makes the Results a
sequence of stages, each closing with a decision. Existing material is
reused wherever named.

**Background.**

- Aggregated N-of-1 trials and the prazosin example (existing).
- RH's framework and its encoding (existing, shortened).
- What a DGP for these trials must satisfy: the criteria list of
  Major 1, numbered, each with one sentence on how it is measured
  (ceiling; coupling gap; agreement of the residual structure with the
  analysis model; calibration status).
- Aims: to take the published process through those criteria, repair
  what fails them, choose what is open, and state what remains
  uncalibrated.

**Methods.** Tools only, no decisions.

- The published process (existing, with the Gompertz display
  corrected).
- Estimands (existing), with a short table of which analysis targets
  which quantity.
- The paired-difference statistic and Propositions 1 to 3 (existing).
- The ceiling (existing).
- Simulation designs in ADEMP form (existing Table 'ademp'), with the
  stress test reduced to the prognostic and size cells that the paper
  uses.
- Software and reproducibility (existing), plus the package version
  and a statement of which form of the cross-component block each
  study uses.

**Results.** Stages, each with: what fails or is open; why; the
candidate; the evidence; the decision and the criterion it rests on.

- Stage 1, Validity of the joint distribution: repair counts
  (Background 1 material, Table 'nonpd'), the ceiling under RH's
  structure, the decision to stop (B.2), and the consequence (RH's 0.6
  infeasible). Criterion: existence.
- Stage 2, The coupling rule: Proposition 2, the step rule's zero gap,
  Table 'closedform', Figures powerfig and gapfig, the patched RH code
  (Table 'patch', Figure fig4bmc), graded coupling as the decision, and
  the limitation that at pharmacological half-lives the gap is within
  0.01 of 1. Criterion: identifiability under carryover; fidelity.
- Stage 3, Carryover on the mean: anchored against compounding, effect
  at most 0.006, decision stated as a specification choice; one plotted
  mean trajectory per path showing the re-exposure step. Criterion:
  fidelity (stated as a judgment).
- Stage 4, Serial correlation: rho as a parameter; ceiling and power
  at matched implied slope for rho in {0.5, 0.6, 0.7, 0.8} under CS;
  the source or sensitivity range for 0.7. Criterion: existence against
  power; calibration status open.
- Stage 5, Cross-component structure: the non-separable form's
  validity condition and the separable form's four properties (A.5),
  the ceilings at rho = 0.7 for both (0.452 against 0.474), the
  matched-residual-structure argument; decision. Criterion:
  analyzability; existence.
- Stage 6, Features the published process cannot represent: prognostic
  couplings with the prognostic null table as evidence that they change
  which analyses are valid; measurement error, decay form and dropout
  listed as options with one sentence each on what the stress test
  showed (or removed if nothing is reported). Criterion:
  analyzability.
- Stage 7, The recommended process: Table 'spec' with a 'Justified in'
  column pointing to the stage; a second column marking each row
  'repair', 'choice' or 'option'.
- Using the process to plan a trial: the section of Major 4.

**Discussion.** The encoding's conditional structure and strengths
(existing), the ceiling as a limitation in outcome units, the
sensitivity model as the route above the ceiling, what remains
uncalibrated, and limitations (existing, plus the on-drug coupling
asymmetry).

**Conclusions.** Existing, reworded so that 'power declines smoothly'
is attributed to half-lives of 0.5 week and above and 'negligibly at
pharmacological half-lives' is given its numerical meaning.

**Appendices.** A (derivations, with the Proposition 3 exclusion and
the 0.10 margin corrected); B (code changes, with B.5 stating the
package's current form of the cross-component block and the
measurement-error status). Remove `appendix-b.Rmd`.

## 9. Scope of this review

Verified (executed from an independent implementation written for this
review, `check.R` in the review scratch directory, using the design
definitions and published parameters read from
`analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R`):

- Ceilings under graded coupling at t1/2 = 0 for CS, the two-constant
  AR(1) form and the separable form at rho = 0.8, 0.7, 0.6 and 0.5,
  for Hybrid, OL+BDC and CO; every entry of Table 'ceilings' agrees
  (0.256, 0.291, 0.440; 0.346, 0.452, 0.474; Hybrid 0.480).
- The smallest eigenvalue of RH's response block at rho_cs = 0.8
  (0.1000).
- Coupling gaps in the Hybrid design, allocation-weighted (0.9997,
  0.9906, 0.9060, 0.7542, 0.5509 at 0.1, 0.2, 0.5, 1, 2 weeks).
- Closed-form E9 and E9b power at c_bm = 0.25, N = 70 (0.721, 0.654
  at t1/2 = 0; E9 0.473 at 1 week; separable 0.335 at rho = 0.8 and
  0.290 at rho = 0.7), by my own noncentral-t implementation.
- Proposition 2's slope at t1/2 = 0 (-0.1302 per unit biomarker; 2.0
  points per SD at c_bm = 0.25).
- The CO placebo-first path's off-drug coupling (identically 0).
- Approximate E9 sample sizes for 80% power under the refined process
  at t1/2 = 0.2 week (280, 140, 84 at c_bm = 0.25, 0.35, 0.45).
- Table 'patch' against `graded-patch-reps500.csv`, and the 0.74, 0.12,
  0.10 values against `published-power.csv` (which records 1,000
  replicates).

Inspected (read, not executed): all prose and tables; the
closed-form script `01-e9-closed-form.R`; the stored curves, per-design
ceilings, factorization check and `e9-vs-lmer.csv`; the stress-test
driver's covariance construction and allocation default; the package
`buildSigma()` cross-component code; NOTATION.md's named-process
entries; the bibliography entries for the sources named in Section 7;
the head of `appendix-b.Rmd`.

Not checked: the stress-test rejection rates and shortfall metrics
(stored outputs not recomputed); the 20-of-54 and 21-of-81 repair
counts; the 0.007 and 0.10 maximum correlation changes; the figure
files; the F1 quick check; the docs/46 protocol; full texts of cited
papers; RH's printed power values; Trials' current author guidelines
(the journal's page redirected to a login; the format points in Minor
16 and 17 are UNVERIFIED for Trials and rest on the BMC house style).

Searches run (no email address or identifier passed to any service):
web search for Percha et al. 2019; for Gärtner et al. 2023; for
Trials methodology-article guidelines (not retrievable); for RH 2020
and its citations; for N-of-1 carryover simulation 2022 to 2024.
PubMed (via the PubMed tool): (i) N-of-1 OR single patient OR
aggregated N-of-1 AND simulation AND biomarker OR effect modifier OR
treatment-by-covariate AND power OR sample size: one record, RH
itself (doi:10.3389/fdgth.2020.00013); (ii) N-of-1 AND simulation AND
carryover OR carry-over OR washout AND design OR power: ten records,
of which Percha 2019, Gärtner 2023, Chen and Chen 2014, Blackston 2019,
Wang and Schork 2019 and Schork 2022 are relevant and are listed in
Section 7 with their DOIs. I did not search for post-2020 critiques of
RH specifically.
