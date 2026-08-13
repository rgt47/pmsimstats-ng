---
geometry: margin=0.85in
fontsize: 10pt
header-includes:
  - \linespread{0.97}
  - \setlength{\parskip}{0.4em}
---

# White paper: the pmsimstats-ng manuscript compendium

*2026-07-29 09:32 PDT, updated 2026-08-12*

**A five-page review of the eleven manuscripts under
`analysis/report/`, with an assessment of the contribution each makes
and a prioritized set of next steps.**

## Scope and method of this review

The compendium comprises eleven manuscripts, each in its own numbered
directory, all targeting *Statistics in Medicine* style and sharing a
common preamble, bibliography apparatus, and simulation substrate.
Several directories carry `report-slim`, `report-extra-slim`,
`report-mini`, or `report_short` variants; these are length-reduced
renderings of the same study, not separate investigations, and are
assessed here as one paper each.

This review is source-level. Every `report.Rmd` was read for its
abstract, section structure, and results, together with the project
READMEs, the six committed referee reports, the paper-04 remediation
plan, and the three existing per-paper white papers (01, 02, 03). No
simulation was re-executed and no numerical result was independently
recomputed for this document, so all quantitative statements below are
reported as the manuscripts report them. Assessments of contribution
and of weakness are the reviewer's judgment; claims about draft state
(placeholder counts, missing results, deferred sweeps) are verified by
inspection of the sources.

**Update, 2026-08-12.** This revision reflects work committed to the
compendium between 2026-07-29 and 2026-08-12, established by
`git diff` against commit `fdfa40e` (the commit that produced the
original text below) and by reading the current working-tree sources.
The comparison is source-level in the same sense as the original: no
simulation was re-executed for this revision either, so a report that
a number is unchanged means the manuscript still states that number,
not that it was independently re-derived. Every paragraph below that
changed is marked inline; paragraphs with no marker were checked
against the current source and found unchanged. Papers 04, 05, and 08
received only copy-edits (British-to-US spelling, an added units
footnote reconciling this series' weekly $t_{1/2}$ convention with 04
and 05's day-denominated one) and are otherwise as originally
assessed.

## What the programme is

The compendium is a single coherent research programme rather than a
collection of related papers. Its object is one estimand: the
biomarker-by-treatment interaction in aggregated N-of-1 trials, the
quantity such trials exist to estimate when they are used to validate a
predictive biomarker. Its motivating application throughout is the
prazosin/PTSD nightmare-frequency setting of Hendrickson et al. (2020),
and its method throughout is Monte Carlo simulation reported under the
ADEMP framework of Morris, White, and Crowther (2019).

The eleven papers partition into five layers.

- **Generative substrate.** 06 (three-component BR/PB/TV
  decomposition) and 07 (Gompertz trajectory template) justify the
  shape of the simulated data.
- **Interaction encoding.** 01 (Architecture A vs B), 11 (Architecture
  C, combined), and 03 (latent-class and mixture formulations) ask how
  the biomarker's moderating effect should enter the data-generating
  process at all.
- **Analyst-side inference.** 02 (carryover specification A1/A2/A3),
  08 (RM-ANOVA vs mixed model vs GEE), and 10 (working-covariance
  miscalibration and cluster-robust repair) ask what the analyst should
  fit.
- **Design.** 05 (twelve-sweep design sensitivity), 09 (informative
  dropout by design family), and the design axis of 08.
- **Clinical framing.** 04 (N-of-1 versus parallel-group RCT for the
  treatment main effect).

The programme's intellectual center is paper 01. Its finding, that the
carryover penalty is an order of magnitude larger when the interaction
is encoded in the covariance (Architecture B) than in the mean
(Architecture A), propagates as a conditioning variable into 02, 03,
07, and 11, and it is the reason those papers report
architecture-stratified rather than pooled conclusions. That is a
genuine and unusual structural virtue: the programme has identified a
modeling choice that the surrounding literature treats as an
implementation detail and has made it an explicit axis everywhere
downstream.

## Paper-by-paper assessment

**01. Two data-generating architectures.** *Useful; the keystone
paper, and materially strengthened since 2026-07-29.* [Updated.]
Establishes that Architecture B loses 30.6% relative power at a
one-week carryover half-life in the OL+BDC design against 2.8% for
Architecture A, with a difference-of-differences $z = 7.64$, while the
crossover design shows no architecture effect; these headline numbers
are unchanged. What changed is the paper around them. It was retitled
("Mean moderation and covariance moderation as data-generating
architectures..."), rebuilt on `bookdown::pdf_document2` so
cross-references resolve, and rewritten end to end: the three-part
`bullets` / `rgt` / `orig` scaffold is gone entirely, replaced by
continuous author prose with zero remaining placeholder blocks, the
first of the eleven papers to reach that state (with 02 and 03,
below). The $N$ convention was also corrected to a uniform $N = 70$
total throughout, resolving the per-path-versus-total ambiguity the
prior draft carried. Most consequentially, the Architecture C material
was fully extracted rather than folded back as this review had
recommended (see "Recommended consolidation," updated below); paper 01
now contains no Architecture C content at all, so the duplication this
review flagged no longer exists, but the super-additivity
benchmark-scale problem travels with the material to paper 11 and is
still unresolved there. Section 6's power-recovery strategies remain
qualitative hypotheses with no simulation behind them. The directory's
`README.md` is now stale in the other direction: it still asserts the
`rgt` blocks carry `rgt to complete.` placeholders, which is no longer
true of `report.Rmd`.

**02. Carryover-mitigation analysis strategy.** *Useful; the strongest
negative result in the compendium, and the scope has narrowed to its
sharpest claim.* [Updated.] The five drafts this review counted
(`report`, `-slim`, `-extra-slim`, `-extra-slim-rgt`, `-devresults`)
were consolidated on 2026-07-30 into a single master `report.Rmd`,
with the others moved to `archive/`; the scaffold is fully retired
here too, with zero placeholder blocks remaining. The paper's scope
narrowed with it: rather than the full 540-cell, three-specification,
three-architecture factorial, the master now restricts to Architecture
B, exponential and Weibull decay, and a head-to-head between the
exposure-weighted specification (E2, formerly A2) and the
lagged-treatment specification (E3, formerly A3) against an unadjusted
baseline (E1, formerly A1), filtering to 216 of the original 540
cells. The "19% of non-null cells" framing this review quoted is gone
along with the broader grid; in its place the paper states directly
that the exposure-weighted advantage is markedly inferior under the
classical crossover design (0.488 against 0.830, unchanged from the
number this review cited) and should not be extended there. Two of the
three weaknesses flagged here are resolved: the CR2 cluster-robust
recalibration is now validated across all five S6 stress cells rather
than one, holding or widening the exposure-weighted advantage at four
of them and narrowing only under high autocorrelation and heavy
dropout; and the previously noted inconsistency with paper 01's
carryover-loss figures was reconciled before this review was written
and remains reconciled, since 01's headline numbers are unchanged. The bias/MSE/
coverage cross-specification comparability caveat persists as a
structural feature of the design (E1/E3 target a different coefficient
than E2), not a defect to be fixed.

**03. Latent-class and mixture formulations.** *Partially useful; a
theory paper with a pilot attached.* Aims 1 and 2 (formalizing the
latent-class DGP, deriving when a single MVN substitutes for it,
characterizing identifiability) are delivered and are the paper's real
contribution, particularly the result that the achievable marginal
correlation is bounded by the geometric mean of the two class
separations. Aim 3 is a 240-replicate pilot of the last of five
pre-registered studies; aim 4, the prazosin re-analysis, has not been
attempted. The pilot's findings are sharp (the unconditional
class-aware Wald test rejects at 0.197 against nominal 0.05; the
BIC-conditional test controls size but is essentially inert at 0.025
power), but they rest on one sample size, one design, and one
generative mechanism, and the regime where a class-aware test should
win, large separation with genuinely bimodal response, was never
tested. [Updated.] The scaffold is now fully written out here as
well, zero placeholder blocks remain, and the paper took the
reframe-rather-than-rerun path this review offered as an alternative:
the abstract and Results now state plainly, in the author's own
prose, that the finding holds "in a 240-replicate-per-cell pilot" and
"within the 240-replicate pilot," rather than presenting the pilot as
if it carried the weight of the pre-registered production run. The
canonical 1,000-replicate run at $N \in \{35, 70, 100\}$ over a wider
gating-slope grid, including a genuinely bimodal separation regime,
has not been executed, so the substantive limitation this review
raised is unchanged; what changed is that the manuscript now discloses
it correctly rather than needing to.

**04. N-of-1 versus parallel-group RCT.** *Useful clinically, but
carries an unaddressed comparison-fairness problem.* Reports that the
aggregated N-of-1 hybrid dominates the parallel RCT across the entire
effect-size range, with roughly a threefold power ratio at the
moderate-effect cell, negligible bias, and near-nominal coverage in
both arms. The ADEMP apparatus is the most complete in the compendium
and the 2026-05-19 remediation plan's critical items (dual
implementations, RNG discipline) appear to have been carried out; the
DerSimonian-Laird estimator now survives only as literature discussion.
The problem is the comparison itself. Both arms are matched at $N = 35$
participants, but the N-of-1 arm contributes eight periods per
participant against the RCT's single post-baseline measurement, so the
result substantially restates that more measurements per subject buy
more information. The Limitations section covers compliance, single
outcome, and decay form, but does not mention observation count or
participant burden at all. A referee at a general medical venue will
raise this immediately. Either add a measurement-matched or
burden-matched sensitivity arm, or state the asymmetry explicitly and
defend it as the decision-relevant comparison for a fixed recruitable
population.

**05. Twelve-sweep design sensitivity.** *Useful as a reference
appendix; weak as a standalone paper.* The twelve sweeps (S1-S12) are
well chosen, each anchored to a literature precedent, and the null
calibration sweep at 2,000 replicates is good practice. The synthesis
is defensible: the Hendrickson hybrid is a competent default; ordering
flips to the parallel RCT only when carryover exceeds the
within-subject period length. But the paper is organized as twelve
one-factor-at-a-time sweeps with a single factorial (S11), so it
characterizes the axes rather than their interactions, and its target
estimand is the main effect, which makes it a companion to 04 rather
than to the interaction papers that constitute the programme's core.
Its most likely fate is as supplementary material to 04.

**06. Three-component decomposition.** *Useful; the most intellectually
substantial paper after 01.* The omitted-variable-bias identity is the
contribution: the one-component biomarker slope is displaced from the
pharmacological slope by a term proportional to the biomarker's
covariance with the placebo-belief and natural-history components, not
by the magnitude of those components. Study A confirms that an
uncontaminated biomarker suffers essentially no bias and that the cheap
phase-augmented remedy is strictly dominated; the contaminated-biomarker
pilot exhibits the failure mode at +0.51 bias; Study B shows recovery
is possible only under a balanced-placebo design that decouples
exposure from belief. That last result, that the cure is a design
property rather than an analysis property, is the paper's real payload
and is what the fourth referee report was asking for. Weaknesses: the
paper is 26,000 words of source, by a wide margin the longest in the
compendium, and needs the slim variant promoted to canonical; the
contaminated-biomarker result rests on a 100-replicate pilot while the
uncontaminated case gets 1,000; and the practical recommendation
(balanced-placebo designs) is rarely feasible in the clinical settings
the programme targets, which the Discussion should concede more
squarely.

**07. Gompertz trajectory evaluation.** *Useful, and materially
improved by revision.* The first-round referee identified a fatal
construction flaw, that the four trajectory families were not four
data-generating processes but a post-hoc imposition on Gompertz-
generated data. The current draft reports an architecture-dependent
answer instead: under Architecture B the trajectory family is exactly
inert because the interaction is mean-independent, so matched seeds
return identical rejection decisions; under a faithfully multiplicative
Architecture A the four families separate by 0.039 in power with
Gompertz the mildly conservative member. This is a better paper than
the one reviewed, and the exact-inertness result is a clean structural
finding worth stating on its own. It remains a narrow paper: one design
(OL+BDC), one sample size, zero carryover, 16 cells. The regimes the
README promises to identify (non-saturating growth, biphasic patterns,
breakpoint behavior) are discussed but not simulated, and the
zero-carryover choice is awkward in a programme whose central finding
is about carryover. [Checked, unchanged.] The 0.039 power-separation
figure, the "exactly inert" result, and the paper's scope (one design,
one sample size, zero carryover, 16 cells) are unchanged in the
current source; the `rgt` layer still carries 51 placeholder blocks.

**08. Test procedure and design.** *Useful; the cleanest practical
message in the compendium.* The closed-form RM-ANOVA derivation is
correct per the referee's independent re-derivation, and the revised
simulation now separates the two effects that the first draft
confounded: dichotomizing the biomarker costs 0.08 to 0.21 in power,
while the further step from a continuous-biomarker contrast to the full
mixed model adds at most 0.04. Naive GEE inflates Type I error to
1.5-1.9 times nominal at all four sample sizes and the Mancl-DeRouen
correction repairs it at roughly ten points of power. The
dichotomization result is the finding most likely to change applied
practice, and it deserves to be the paper's headline rather than a
component of it. The unresolved problem is scope: the cycle-by-period
design grid is specified in detail and then deferred, which leaves the
title's second half unearned. Either run the grid or retitle the paper
around the test-procedure axis.

**09. Informative dropout by design.** *Promising, and closer to a
paper than it was.* [Updated.] The gap is real and correctly
identified: Hendrickson's Figure 4A implies a design-by-dropout
coupling that the N-of-1 methods literature has not pursued, and no
other manuscript here centers dropout. The abstract's Results
paragraph no longer reads "[Forthcoming]" (one internal methods note
still uses that word, but the reader-facing Results and Discussion do
not), and a `cover-letter.md` has been added. The ceiling-saturation
problem this review raised is partly addressed: the design grid was
re-run at a matched $N = 70$ and now reports all four designs rather
than two, hybrid 0.968, OL+BDC 0.950, traditional crossover 0.854, and
open-label 0.088, so the ordering has real resolution between crossover
and the top two designs even though hybrid and OL+BDC still sit close
to the ceiling and are explicitly flagged in the text as
statistically inseparable ($z = 1.44$). The manuscript itself now
argues that the hybrid design's earlier near-immunity to dropout was a
ceiling artifact of running at half the matched sample size, which is
the same diagnosis this review made independently. [Verified
2026-08-12.] Two follow-up items are now checked directly. First, the
verbatim-duplicate paragraph this review originally flagged is gone:
Results (report.Rmd:947-961) and Discussion's "Prototype Monte Carlo
confirmation" subsection (report.Rmd:1304-1328) cover the same
prototype run but are no longer the same text; the Discussion version
is longer and adds the per-path $N$ breakdown and the contrast with
the earlier, unmatched 35-per-path run, so this is no longer a defect.
Second, Architecture B and a wider $c_{bm}$ range were not added: the
prototype remains single-architecture (mean-moderation only, at one
value, $c_{bm} = 0.45$) and covariance-moderation coverage for the OL
row is discussed only as a scoped future extension
(report.Rmd:1423, 1439, 1561-1575), not run. The `rgt` layer still
carries 50 placeholder blocks (still Lorem ipsum throughout the
`.rgt` divs, e.g. report.Rmd:944, 1301).

**10. Interaction-test calibration.** *Useful, and the most portable
paper in the compendium.* It is the only manuscript whose contribution
is not specific to N-of-1 trials: the finding is that a misspecified
parametric working covariance biases the model-based standard error of
a mixed-model interaction asymptotically, in either direction. The
diagnostic staging is convincing, ruling out degrees-of-freedom rules,
positive-definiteness correction, and point-estimate bias, then
localizing the 6-10% standard-error inflation to the corCAR1 residual
layer, demonstrating that the bias does not attenuate from $N = 35$ to
$N = 280$, and repairing calibration with a CR2 cluster-robust standard
error and independently from the marginal side with bias-corrected GEE.
The reusable staged protocol is a real deliverable. Two gaps: the
worked example is the same prazosin cell as everything else, so the
generality claim rests on argument rather than on a second application;
and the paper has no README and no cover letter, indicating it has not
been through the same submission-preparation pass as its neighbors.
[Checked, unchanged.] Only cosmetic edits since 2026-07-29: US
spelling corrections and relabeling the analysis specifications from
A1/A2/A3 to E1/E2/E3 to stay synchronized with paper 02's new naming.
Both gaps persist, the `rgt` layer still carries 37 placeholder
blocks, and the directory still has no `README.md` or cover letter.

**11. Combined DGP architecture (Architecture C).** *Still not
assessable, and the extraction this review questioned is now
complete rather than reversed.* [Updated.] The author did not adopt
this review's recommendation to fold 11 back into 01; instead the
split was finished. Paper 01's `report.Rmd` no longer contains any
Architecture C material, so the duplication this review originally
flagged is gone, but that also means paper 11 now carries the entire
burden of the material this review found weakest, with none of the
consolidation risk addressed. The manuscript's own abstract still
states "[To be finalized after the common-driver boundary re-run]",
and the super-additivity claim still rests on the language "will
confirm agreement to within Monte Carlo standard error" rather than a
completed comparison, so the benchmark-scale problem (probability
scale versus the effect-size scale this review argued for) has not
been resolved. Twelve `rgt` placeholder blocks remain. The
crossover-versus-OL+BDC panel question and the non-identifiability of
$c_{bm,A}$ and $c_{bm,B}$ from a fitted `bm:Dbc` coefficient were not
independently re-checked in this pass.

## What the compendium establishes collectively

Read together, the eleven papers deliver a coherent and defensible set
of claims about aggregated N-of-1 biomarker-validation trials.

1. How the interaction is encoded in the data-generating process is a
   scientific assumption with large consequences under carryover, not
   an implementation detail (01, with 03 and 11 mapping the
   neighborhood).
2. Analyst-side defaults are conditional, not universal. The
   exposure-weighted carryover predictor helps only under one
   architecture and one class of design (02); the mixed model beats
   GEE and RM-ANOVA, but most of that gap is dichotomization rather
   than the test (08); the model-based standard error is not
   trustworthy when the working covariance is a convenient parametric
   guess (10).
3. Design dominates analysis. Retention beats specification choice
   (02); the balanced-placebo contrast, not any estimator, is what
   recovers an unbiased pharmacological slope under a contaminated
   biomarker (06); design family governs dropout robustness (09).
4. Two modeling choices that might have been feared are largely
   inert: the trajectory family (07, exactly inert under Architecture
   B) and mixture versus continuous encoding at trial-relevant $N$
   (03, where the standard mixed model dominates a class-aware test).

Negative and inertness results of this kind are the compendium's
comparative advantage, and it is worth being deliberate about that in
the cover letters. The programme has repeatedly tested its own defaults
and reported when they failed.

## Cross-cutting weaknesses

**Author prose was absent everywhere; it no longer is, in three
papers.** [Updated.] At the original baseline all eleven manuscripts
used the three-part `bullets` / `rgt` / `orig` scaffold with no author
text in any `rgt` block, 645 placeholder blocks in total. As of
2026-08-12, papers 01, 02, and 03 have been rewritten as continuous
author prose with the scaffold removed entirely, zero placeholder
blocks in each. The remaining eight still carry placeholders: 04
(110), 06 (101), 08 (59), 05 (54), 07 (51), 09 (50), 10 (37), and 11
(12), 474 blocks in total. This is progress on the compendium's single
binding submission constraint, but eight of the eleven papers,
including paper 10, which this review's next-steps list prioritized
third, are unchanged on this dimension.

**One reference parameter set.** Nearly every simulation is calibrated
to the same prazosin/PTSD trajectory. Absolute power values are
therefore illustrative throughout, and only orderings port. Paper 10 is
the case where this matters most, because its claim is explicitly
general.

**Ceiling saturation.** Papers 09 and, in places, 04 and 05 report
power at or near 1.00, where designs and methods cannot be
discriminated. Cells should be chosen to sit near 0.5 for the contrast
of interest.

**Prototype-versus-production gap.** Papers 03, 06 (contaminated
pilot), 08 (deferred design grid), and 09 present pilots or partial
grids as if they carried the weight of the pre-registered production
runs. Each states this, but the abstracts do not always reflect it.

**Duplicate and stale scaffolding.** [Updated.] The READMEs that
asserted driver directories "do not yet exist" have been corrected;
none remain as of this pass. Paper 09 gained a `cover-letter.md`.
Papers 10 and 11 still lack READMEs and cover letters. Paper 01's
README is now stale in a new way (see paper 01, above): it still
describes the `rgt` layer as unfilled placeholders after the paper was
rewritten past that state. Paper 09's Results/Discussion duplication is
confirmed resolved: the two sections cover the same prototype run in
different prose, not a verbatim repeat.

## Recommended consolidation

Eleven manuscripts on a single estimand, calibrated throughout to one
prazosin/PTSD parameter set, invites the charge of over-partitioning,
and a referee who encounters two of them will notice. Six papers carry
a defensible independent claim: 01 (the architecture distinction), 02
(the conditional failure of the exposure-weighted default), 06 (the
omitted-variable identity and the design-side cure), 08 (the
dichotomization cost), 10 (working-covariance miscalibration, the only
contribution not specific to N-of-1 trials), and 04 (the clinical
comparison, conditional on repairing its comparison basis). Three do
not, and the objection in each case is to the packaging rather than the
content, all of which is worth keeping.

- **Fold 11 back into 01** as a scoped sensitivity section rather than
  completing the extraction. Architecture C is the common
  generalization of A and B, not a rival to them, and its two channel
  parameters are not separately identifiable from a fitted `bm:Dbc`
  coefficient, so it is a simulation-side construct with no estimation
  counterpart. Its two legitimate uses, a boundary-cell correctness
  gate on the DGP implementation and a mixing-weight sweep when the
  mechanism is unknown, are arguments about how to use 01's framework
  and read naturally as part of 01. **Status as of 2026-08-12: not
  adopted.** The author completed the extraction instead, and paper 01
  is now clean of Architecture C material while paper 11 stands alone,
  still unfinished, still carrying the unresolved benchmark-scale
  problem. This recommendation is renewed rather than withdrawn: paper
  11 is 3,800-odd words against a compendium median near 13,000, has no
  README or cover letter, and its abstract still reads as
  "to be finalized." If the author's judgment is that Architecture C
  earns independent-manuscript status on scientific grounds, that case
  should be made explicitly in 11's Introduction, since nothing in the
  current source argues for it.
- **Publish 05 as supplementary material to 04.** Its twelve sweeps
  target the treatment main effect, not the interaction, which makes it
  a companion to 04 rather than to the programme's core, and its
  one-factor-at-a-time structure characterizes axes rather than their
  interactions.
- **Merge 07 into 06 as a robustness section on the generative
  substrate.** Both concern how the simulated trajectories are shaped;
  neither is large on its own; and 07's cleanest result, that
  trajectory family is exactly inert under a covariance-channel DGP,
  is a statement about the substrate that 06 formalizes. If 06 is
  slimmed for submission as recommended above, 07's sixteen cells fit
  comfortably in the technical companion.

If adopted, the result would be seven or eight submittable papers of
more even weight: 01 (absorbing 11), 02, 04 (with 05 appended), 06
(absorbing 07), 08, 10, and 03 and 09 once their outstanding work
completes. As of 2026-08-12 only the 03/09 half of that plan has moved
forward on its own terms (03 by reframing, 09 by re-running toward
resolvable cells); the 01/11 fold was not adopted, and 05-into-04 and
07-into-06 have not been acted on either.

## Recommended next steps, in priority order

1. **Write the `rgt` prose for paper 01, then 02, then 10.** *Status:
   01 and 02 done (2026-08-12); 03 done as well, ahead of its place in
   this list, since it was rewritten in the same pass. Papers 06, 04,
   08, 05, 07, 09, and 10 remain, in descending order of remaining
   placeholder count (101, 110, 59, 54, 51, 50, 37). 10 is still next
   in priority even though it is not the largest remaining count,
   because it is the compendium's only non-N-of-1-specific
   contribution and this review's next-most-portable paper after 01.*
   Run `~/bin/strip-claudecode` on each as its `rgt` layer closes; 01
   and 02 have not yet had it run, since neither has a `README.md`
   confirming submission-preparation completion in the sense the other
   finished papers do.
2. **Settle the Architecture C question and re-verify 01's boundary
   cells.** *Status: not settled, and settled in the opposite direction
   from this review's recommendation.* The extraction was completed
   rather than reversed (see "Recommended consolidation," above). What
   remains outstanding is unchanged: recompute the super-additivity
   comparison on the effect-size scale, replace the crossover carryover
   panel with OL+BDC, and confirm that $(c_{bm,A} = 0.45, c_{bm,B} = 0)$
   and $(0, 0.45)$ reproduce 01's published power figures exactly
   before accepting the grid. The 01/02 carryover-loss consistency has
   already been re-checked and holds, since 01's published numbers did
   not move.
3. **Fix paper 04's comparison basis.** Add either a
   measurement-matched or a burden-matched RCT arm, or a stated,
   defended asymmetry, and put observation count in the Limitations.
   This is a one-sweep change and it removes the paper's most likely
   rejection reason.
4. **Re-scope paper 09 around resolvable cells.** *Status: partially
   done, now verified.* The design-by-dropout grid was re-run at a
   matched $N = 70$ and now includes all four designs rather than two,
   giving traditional crossover (0.854) and open-label (0.088) genuine
   separation from the hybrid/OL+BDC pair, and the "[Forthcoming]"
   abstract placeholder is gone. The duplicated-paragraph defect this
   review originally flagged is confirmed removed: Results and
   Discussion cover the same prototype run in different text, not a
   verbatim repeat. Still open, and now confirmed rather than merely
   flagged: Architecture B has not been added (the prototype remains
   single-architecture, mean-moderation only) and $c_{bm}$ is still
   evaluated at one value, 0.45, with no range.
5. **Decide paper 08's scope.** Unchanged; no commits to this paper's
   `report.Rmd` since 2026-07-29 beyond a spelling pass.
6. **Promote slim variants to canonical for 06 and 02.** *Status: done
   for 02, open for 06.* Paper 02's consolidation folded its slim
   variant into the single master `report.Rmd` along with the scope
   narrowing described above, which should bring it within journal
   budget; this was not independently checked by page count in this
   pass. Paper 06 is unchanged, still roughly 26,000 words of source
   against `report-slim.Rmd`'s roughly 18,700, and the slim variant has
   not been promoted.
7. **Add a second application to paper 10.** Unchanged; no commits to
   this paper's content since 2026-07-29 beyond a spelling and
   specification-label pass.
8. **Run the canonical replication for paper 03, or reframe it.**
   *Status: reframed, not re-run.* The abstract and Results now state
   directly, in the author's own prose, that the finding holds "in a
   240-replicate-per-cell pilot," rather than presenting it as the
   pre-registered production result. The 1,000-replicate run at
   $N \in \{35, 70, 100\}$ over a wider gating-slope grid, including a
   genuinely bimodal separation regime, has not been executed and would
   still be required before the non-competitiveness result generalizes
   beyond the pilot.
9. **Housekeeping.** *Status: partially done.* The stale
   "do not yet exist" README language has been refreshed compendium-
   wide, and paper 09 gained a cover letter. Still open: READMEs and
   cover letters for 10 and 11 (11 also still needs the compendium
   table entry in `analysis/report/README.md` extended to cover it,
   not independently re-checked this pass), and paper 01's own README
   now needs refreshing in the opposite direction, since it undersells
   a paper that no longer carries `rgt` placeholders.

A reasonable submission sequence, given that the consolidation above
was not adopted for paper 11, is 01 first as the keystone (it no
longer depends on 11 to be complete), then 02 and 10, both close to
complete and 02 already carrying a clean practical message; then 08
once its scope decision is made; then 06 in slim form with 07 as a
robustness section, if that consolidation is adopted; then 04 with 05
as supplementary material, once the comparison basis is repaired; with
03, 09, and 11 last, as their outstanding work, reframing accepted for
03, resolution and duplication cleanup for 09, and the benchmark-scale
and completion problems for 11, is finished or explicitly closed.
