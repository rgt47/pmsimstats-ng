---
geometry: margin=0.85in
fontsize: 10pt
header-includes:
  - \linespread{0.97}
  - \setlength{\parskip}{0.4em}
---

# White paper: the pmsimstats-ng manuscript compendium

*2026-07-29 09:32 PDT*

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

**01. Two (three) data-generating architectures.** *Useful; the
keystone paper.* Establishes that Architecture B loses 30.6% relative
power at a one-week carryover half-life in the OL+BDC design against
2.8% for Architecture A, with a difference-of-differences $z = 7.64$,
while the crossover design shows no architecture effect. The mechanism
is argued cleanly: Architecture B's signal lives in a correlation that
decays off drug, so carryover attacks it twice. Type I error is checked
across 18 null cells and an $N = 70$ robustness run preserves the
ordering. This paper should be submitted first; the rest of the
compendium cites it. Two weaknesses persist. The Architecture C
material now duplicated in paper 11 uses an additive-on-the-probability
-scale benchmark to claim super-additivity, which is the wrong
benchmark for a convex power function and inflates the apparent effect
from roughly three Monte Carlo standard errors to nine. Section 6's
power-recovery strategies are qualitative hypotheses with no simulation
behind them.

**02. Carryover-mitigation analysis strategy.** *Useful; the strongest
negative result in the compendium.* A 540-cell factorial establishes
that the exposure-weighted specification A2, the framework's own
default, attains the highest power in only 19% of non-null cells. It
wins by about ten points under Architecture B in designs with a blinded
discontinuation phase and loses badly under the classical crossover
(0.488 against 0.830). A3, the Jones-Kenward lagged-term device, is
shown to be numerically almost identical to A1 because its lagged term
enters only as a main-effect nuisance. The paper also reports that 30%
dropout collapses all three specifications to indistinguishability,
which is the most actionable single sentence in the compendium.
Remaining weaknesses: bias, MSE, and coverage in the headline table are
scored against A2's estimand and are therefore not comparable across
specifications (the paper says so, but the table invites the
comparison anyway); the CR2 recalibration is validated in one stratum
only; and the previously noted inconsistency with paper 01's
carryover-loss figures should be re-verified after 01's final numbers
settle.

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
tested. As drafted the paper claims less than it appears to; it should
either be reframed explicitly as theory-plus-pilot or the canonical
1,000-replicate run should be executed before submission.

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
is about carryover.

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

**09. Informative dropout by design.** *Promising but not yet a paper.*
The gap is real and correctly identified: Hendrickson's Figure 4A
implies a design-by-dropout coupling that the N-of-1 methods literature
has not pursued, and no other manuscript here centers dropout. But the
abstract's Results paragraph still reads "[Forthcoming]", the empirical
content is a 500-replicate 16-cell prototype at one sample size, one
value of $c_{bm}$, and Architecture A only, and the same prototype
paragraph appears verbatim in both Results and Discussion. More
seriously, the reported power ordering is measured at the ceiling:
hybrid at approximately 1.00 and OL+BDC at approximately 0.95 cannot
discriminate designs, and OL sits at nominal alpha by construction
because it supplies no off-drug contrast. The design ordering claim
therefore has almost no resolution in the cells that produced it. The
"happy accident" randomization-path decomposition, which is the paper's
most original idea, needs a cell where power is near 0.5 to be
visible at all.

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

**11. Combined DGP architecture (Architecture C).** *Not yet
assessable; currently an extraction in progress.* The manuscript is
3,762 words of source against a compendium median of about 13,000, has
no rendered PDF, and its abstract states "[To be finalized after the
common-driver boundary re-run]". Its content is the Architecture C
material currently also present in paper 01, being split out per
`01-.../architecture-c-extraction-plan.md`. Whether the split is
warranted is doubtful: C is a sensitivity device rather than a co-equal
third architecture, and the case for it is better made inside 01 than
as an independent manuscript (see the consolidation section below).
Two substantive issues must be fixed either way. The super-additivity
claim needs the effect-size-scale
benchmark rather than the probability-scale one, under which the
observed excess shrinks to roughly three Monte Carlo standard errors.
And the grid's carryover comparison is currently reported only for the
crossover design, the one design in which paper 01 found no
architecture-by-carryover interaction; the informative panel is OL+BDC.
The paper should also state plainly that $c_{bm,A}$ and $c_{bm,B}$ are
not separately identifiable from a fitted `bm:Dbc` coefficient, which
makes C a simulation-side construct only.

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

**Author prose is absent everywhere.** All eleven manuscripts use the
three-part `bullets` / `rgt` / `orig` scaffold, and no `rgt` block in
any paper contains author text: papers 01, 02, and 11 hold a "to
complete" placeholder and the remaining eight hold Lorem ipsum, 645
blocks in total. This is the single binding constraint on submission and it is
not a formatting task; the `rgt` layer is where the author's voice and
argument are supposed to live.

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

**Duplicate and stale scaffolding.** Several READMEs still assert that
driver directories "do not yet exist" when they do; paper 09 repeats a
paragraph verbatim across two sections; papers 10 and 11 lack READMEs
and cover letters. These are cheap to fix and they signal draft state
to a referee.

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
  and read naturally as part of 01.
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

The result is seven or eight submittable papers of more even weight:
01 (absorbing 11), 02, 04 (with 05 appended), 06 (absorbing 07), 08,
10, and 03 and 09 once their outstanding simulation work completes.

## Recommended next steps, in priority order

1. **Write the `rgt` prose for paper 01, then 02, then 10.** Nothing
   can be submitted until this is done, and 01 gates the citation chain
   for the rest. Papers 02 and 10 are the next most nearly complete.
   Run `~/bin/strip-claudecode` on each as its `rgt` layer closes.
2. **Settle the Architecture C question and re-verify 01's boundary
   cells.** Adopt the consolidation above (fold 11 into 01) or defend
   the extraction explicitly. Either way, recompute the
   super-additivity comparison on the effect-size scale, replace the
   crossover carryover panel with OL+BDC, and confirm that
   $(c_{bm,A} = 0.45, c_{bm,B} = 0)$ and $(0, 0.45)$ reproduce 01's
   published power figures exactly before accepting the grid. Re-check
   the 01/02 carryover-loss consistency afterward.
3. **Fix paper 04's comparison basis.** Add either a
   measurement-matched or a burden-matched RCT arm, or a stated,
   defended asymmetry, and put observation count in the Limitations.
   This is a one-sweep change and it removes the paper's most likely
   rejection reason.
4. **Re-scope paper 09 around resolvable cells.** Re-run the
   design-by-dropout grid at parameter values placing power near 0.5,
   add Architecture B, extend beyond a single $c_{bm}$, and write the
   Results paragraph of the abstract. Drop the duplicated prototype
   paragraph.
5. **Decide paper 08's scope.** Either run the cycle-by-period grid or
   retitle around the test-procedure axis and promote the
   dichotomization result to the headline.
6. **Promote slim variants to canonical for 06 and 02.** Both exceed
   plausible journal budgets in their full form. The full versions
   become technical companions in `docs/`.
7. **Add a second application to paper 10.** A non-N-of-1 longitudinal
   example, even a small one, converts its generality claim from
   argument to demonstration and makes the paper submittable to a
   broader methods venue.
8. **Run the canonical replication for paper 03, or reframe it.** A
   1,000-replicate run at $N \in \{35, 70, 100\}$ over a wider
   gating-slope grid, including a genuinely bimodal separation regime,
   is required before the non-competitiveness result generalizes. If
   that compute is not available, reframe as theory plus pilot and say
   so in the abstract.
9. **Housekeeping.** Refresh the stale READMEs, add READMEs and cover
   letters for 10 and 11, update the compendium table in
   `analysis/report/README.md` to cover papers 10 and 11, and record
   the intended submission order.

A reasonable submission sequence, assuming the consolidation above, is
01 (absorbing 11) and 10 first, being the keystone and the most
portable; 02 and 08 next, both close to complete and both carrying
clean practical messages; then 06 in slim form with 07 as a robustness
section; then 04 with 05 as supplementary material, once the
comparison basis is repaired; with 03 and 09 last, as their
outstanding simulation work completes.

---
*Rendered on 2026-07-29 at 09:32 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/whitepaper-compendium-summary.md*
