# White paper: referee-style review of "Mean moderation and covariance
moderation as data-generating architectures for biomarker-treatment
interactions in aggregated N-of-1 trials, and their divergent power
under carryover"

*Review date: 2026-08-26 08:22 PDT*

Scope note: this review is deliberately limited to the manuscript in
this directory (`analysis/report/01-dgp-mean-moderation-vs-mvn/`),
by the author's explicit choice, and does not re-review the other ten
report series in the pmsimstats-ng compendium. Portfolio-wide reviews
already exist at `docs/pub_review_whitepaper_2026-08-16.md` and
`docs/pub_review_whitepaper_2026-08-20.md` at the repository root. An
older, differently structured referee report for an earlier draft of
this manuscript exists at `../referee-report-2026-06-13.md`; it is
noted here as prior context but its point-by-point items are not
individually re-litigated, since the master manuscript has since been
substantially rewritten (the tutorial-style consolidation of
2026-07-31 described in `../README.md`).

## 1. Summary of the work under review

The manuscript (`report.Rmd`, rendering to a 17-page PDF) compares two
architectures for encoding a biomarker-treatment interaction in the
data-generating process (DGP) of a Monte Carlo power simulation for
aggregated N-of-1 trials: mean moderation, in which the biomarker
scales the treatment effect directly in the conditional mean, and
covariance moderation, in which the interaction is represented as a
drug-exposure-dependent correlation in a joint multivariate normal
(MVN) distribution, following Hendrickson et al. (2020). Holding the
analysis model fixed at a linear mixed model with continuous-time
AR(1) residual correlation, the paper runs a 36-cell factorial (2
architectures x 3 designs x 3 carryover half-lives x 2 biomarker
effect sizes, 1,000 replicates per cell, N = 70 total per design) and
finds that covariance moderation loses far more power to carryover
than mean moderation, most sharply in the open-label-plus-blinded-
discontinuation (OL+BDC) design (30.6% relative loss versus 2.8%,
$z = 7.64$). A secondary N = 140 run on OL+BDC alone is offered as a
robustness check that the ordering survives power saturation. The
paper closes with a discussion of which biological mechanisms favor
each architecture, reporting recommendations for future simulation
studies, and untested mitigation strategies.

During this review, at the author's request, Section 4.1 was also
extended in place with material drawn from the companion latent-class
manuscript (`analysis/report/03-latent-class-mixture-application/`),
deepening the discussion of which drug and biomarker properties, and
which trial-design features, favor each architecture, and of the
boundary regimes (near-deterministic biomarkers, bimodal responder
subpopulations, class-varying covariance) in which neither pure
architecture is adequate. That edit is described in Section 4 below
and has already been applied to `report.Rmd`; it is not a pending
recommendation.

## 2. Major issues

### 2.1 The N = 140 robustness check (Section 3.3) has no traceable generating script or output data in the repository

**Location.** `report.Rmd`, Section 3.3, "Robustness check at
N = 140" (the table reporting covariance-moderation power 0.978,
0.945, 0.885 and mean-moderation power 0.976, 0.973, 0.959 across
$t_{1/2} \in \{0, 0.5, 1.0\}$, described as an 800-replicate-per-cell,
12-cell, 9,600-fit run restricted to the OL+BDC design).

**Problem.** Unlike the Section 3.1 production table, which
reproduces exactly from `analysis/data/quick-sim/01-dgp-summary.txt`
(verified by direct comparison; see Section 2.2 below for the
verification), no script in `analysis/scripts/`, and no data file in
`analysis/data/`, references an N = 140, OL+BDC-only, 800-replicate
run. A targeted search covered: every `*summary*.txt` file under
`analysis/data/` (none contains an N = 140 row); every script under
`analysis/scripts/` for the literal string "N = 140" (only an
unrelated file in `treatment-main-effect/plan1.md` matches); and the
git history for a deleted script matching that description (none
found). One N = 140 dataset does exist, in
`analysis/.n-migration-backup/quick-sim/01-dgp-replicates-unmatched.rds`,
but inspection (`table(x140$design)`) shows it covers the Hybrid
design only, with 12,000 replicate-level rows, and is a leftover from
the pre-N-matching parameterization era referenced in the
Reproducibility section of this same manuscript — it is not the
OL+BDC N = 140 dataset the manuscript describes.

**Why it matters to a referee.** This is precisely the failure mode
the pub_review criteria ask reviewers to distinguish: "verified
results (code exists and runs)" versus "asserted results (numbers in
prose with no reproducible source)." A specific numbered results
subsection, with a $z$-statistic and a stated fit count, currently has
no reproducible source anywhere in the repository. This is the same
class of defect that a companion review in this compendium (paper 07,
"29-mixnormalmri" and "30-missforest" in the portfolio-wide whitepaper)
flagged as fabricated-appearing content when a results section
described numbers that did not trace to a live run. This review found
no evidence the numbers were fabricated (they are internally
consistent and plausible given the Section 3.1 pattern, and the
manuscript is explicit that this is a smaller, less powered check), but
that is an inference from plausibility, not a verification, and the
distinction matters for a referee.

**Remediation required.** Either (a) locate the original driver script
and output (check other machines, prior branches, or the `share/`
PDF history for a version whose commit predates any deletion), or (b)
re-run the described N = 140, OL+BDC, 800-replicate-per-cell factorial
from the existing `01-dgp-prototype.R` driver with an N = 140
parameterization and confirm the reported numbers reproduce within
Monte Carlo error, or (c) revise Section 3.3 to state plainly that this
check has not yet been re-executed under the current N-matching
convention and either remove the specific numbers or label them
explicitly as a legacy result predating the N-matching correction
documented in the Reproducibility section. Option (c) is the fastest
path to correctness; option (b) is preferable for a final submission
because it restores a genuine, currently-absent robustness check.

### 2.2 Section 3.1 and 3.4 numeric results are verified reproducible from committed data (positive finding, stated for calibration)

**Location.** `report.Rmd`, Sections 3.1 and 3.4.

**Finding.** Every power figure in the two Section 3.1 tables (CO
0.739/0.739/0.725 and 0.745/0.754/0.723; Hybrid 0.842/0.844/0.784 and
0.840/0.795/0.711; OL+BDC 0.751/0.748/0.730 and 0.777/0.694/0.539)
reproduces exactly, to three decimal places, from
`analysis/data/quick-sim/01-dgp-summary.txt` (36 data rows, matching
the stated 36-cell design). The Type I error range [0.010, 0.074] and
mean 0.042 in Section 3.4, and the $\hat\beta_{bm:D}$ range -0.230 to
-0.281, both reproduce from the same file by direct recomputation
(verified: independently aggregated the 18 null cells and the 18
alternative cells and recomputed the summary statistics). This is
stated here, rather than only in Section 3, because it calibrates the
severity of finding 2.1: the paper's primary results (Section 3.1,
which the discussion explicitly identifies as "the canonical reference
for the narrative magnitude") are solidly verified, and the
reproducibility gap is confined to the secondary robustness check.

### 2.3 No systematic verification of the analytical claims in Section 2.2.3

**Location.** `report.Rmd`, Section 2.2.3, "Analytical comparison."

**Problem.** The algebraic claims in this subsection (that the
mean-moderation biomarker-effect ratio is invariant to proportional
scaling of $D_{it}$, and that the covariance-moderation conditional
expectation converges to $\mu_{BR,\text{off}}$ as $t_{sd} \to \infty$)
are correct by direct inspection of the stated equations and are not
in question. What is not established anywhere in the paper is a
quantitative link between this closed-form mechanism and the specific
magnitudes reported in Section 3 (2.8% versus 30.6% relative loss).
The paper is explicit that "no such closed form is available" for the
full mixed-model power (Introduction), which is a defensible position,
but the consequence is that Section 2.2.3's analytical argument
establishes only the qualitative direction of the effect (covariance
moderation loses more), not its size, and the paper should say so
plainly rather than let the juxtaposition of Section 2.2.3 and the
Section 3 tables imply a tighter connection than exists.

**Remediation required.** Add one sentence at the end of Section 2.2.3
stating explicitly that the closed-form comparison establishes
direction, not magnitude, and that the magnitude is established solely
by the Monte Carlo results of Section 3. This is a low-cost
clarification, not a new derivation.

## 3. Minor issues

### 3.1 `README.md` contains two stale claims

**Location.** `../README.md` (this directory).

The README's "Master source" bullet list states, as of 2026-08-16,
that `report.Rmd` "carries zero remaining `rgt` placeholder markers"
and has "reached zero-placeholder state" — this is correct; the master
carries no `bullets`/`rgt`/`orig` scaffolding tags at all (verified by
grep). But the same file's "Known defects in the master" section, a
few lines below, still lists "The `rgt` blocks throughout are
placeholders reading `rgt to complete.`" as an open defect, directly
contradicting the zero-placeholder claim above it. The same "Known
defects" section also lists a raw-citekey rendering defect for
`[@pmsimstats-paper08]`; this citekey is registered in
`references.bib` (`@unpublished{pmsimstats-paper08, ...}`) and appears
in Section 1 as plain markdown syntax, not inside a raw-LaTeX
scaffolding environment, so the defect as described no longer applies
to the current master (it may be a holdover from when the comprehensive
archived draft, which does retain scaffolding tags, was the reference
point). Both stale claims should be removed from the "Known defects"
section.

### 3.2 README page-count claim does not match the current render

**Location.** `../README.md`, "Master source" section: "The current
render is 21 pages."

The manuscript as rendered during this review (including the Section
4.1 extension applied below) is 17 pages (`pdfinfo report.pdf`). This
may simply predate a prior trim of the manuscript; either way the
figure should be refreshed the next time the README is touched.

### 3.3 Type I error median in Section 3.4 rounds inconsistently

**Location.** `report.Rmd`, Section 3.4: "mean 0.042 and median
0.035."

Independent recomputation from the 18 null cells in
`01-dgp-summary.txt` gives a median of 0.0345 (the two central sorted
values are 0.034 and 0.035). Reporting 0.035 rather than 0.034 or
0.0345 is a trivial rounding choice with no bearing on the argument,
noted only for completeness.

### 3.4 Mitigation strategies in Section 4.3 are explicitly untested, which is appropriate but easy to overlook on a skim

Not a defect: Section 4.3 already states plainly, "None of these
strategies was tested by simulation here." This is good epistemic
practice and is flagged here only to record, for the "what remains to
be done" checklist below, that this is real, not yet started work
rather than a completed contribution.

## 4. What remains to be done

**(a) Required for correctness.**

1. Resolve finding 2.1: either regenerate or clearly re-label the
   Section 3.3 N = 140 robustness check. This is the single blocking
   item; everything else in this manuscript is either verified or is
   an appropriately labeled limitation.
2. Add the one-sentence magnitude-versus-direction clarification to
   Section 2.2.3 (finding 2.3).

**(b) Required for acceptance.**

3. Clean up `README.md`'s self-contradictory scaffolding claim and its
   stale citekey-rendering defect note (finding 3.1); refresh the page
   count (finding 3.2).
4. Section 4.1 has been extended during this review (see Section 5
   below) to connect the mean-versus-covariance-moderation choice to
   the boundary regimes and design-dependence mechanisms developed in
   the companion latent-class paper. A referee reading this extension
   alongside the original text should confirm no redundancy was
   introduced with Section 4.4's existing "Relation to the wider
   literature" discussion of trial-design dependence; a light pass
   found none, but the author should re-read both together once more
   before submission, since both sections now discuss trial-design
   dependence from adjacent angles.

**(c) Desirable polish.**

5. Section 4.3's mitigation strategies remain untested by design (a
   companion paper is deferred to take this up); no action needed for
   this manuscript, only confirmation the companion paper is tracked
   elsewhere (it is: `analysis/report/` series numbering suggests this
   is future, unassigned work, not yet a numbered report).
6. Consider adding a one-line data-availability pointer in the
   Reproducibility section directly to
   `analysis/data/quick-sim/01-dgp-summary.txt` and
   `analysis/scripts/quick-sim/01-dgp-prototype.R`, since this review
   had to locate them by search; an explicit path would let a referee
   verify Section 3.1 in under a minute rather than requiring a
   repository-wide search.
7. `analysis/scripts/quick-sim/run-01-dgp.sh` hard-codes
   `cd /Users/zenn/Dropbox/prj/alz/10-pmsimstats-ng/pmsimstats-ng`, a
   path outside the current repository location
   (`.../prj/res/36-pmsimstats-ng/pmsimstats-ng`). This does not affect
   the manuscript's correctness, since `01-dgp-prototype.R` itself uses
   relative paths and can be run directly with
   `Rscript analysis/scripts/quick-sim/01-dgp-prototype.R` from the
   repository root, but the wrapper script would fail if invoked as
   written and should be corrected or removed the next time this
   simulation is re-run.

## 5. Recommended framing

**(a) Plausible framings.** (i) A methodological caution paper: "DGP
architecture choice materially changes power estimates under
carryover; report both, or justify your choice." (ii) A software/tools
paper documenting the pmsimstats-ng dual-architecture capability. (iii)
A biology-first argument paper: "covariance moderation is the correct
model for the prazosin-PTSD biomarker, and here is what it costs."

**(b) Recommendation.** Framing (i), which is close to what the
manuscript already does, remains the right choice, and this review
found no reason to change that recommendation from what the
manuscript's own framing already assumes. Framing (ii) undersells the
scientific content (this is not primarily a software paper; the
software is instrumental). Framing (iii) is too narrow: the paper's
strongest, most exportable claim is the general one, that
mean-moderation power estimates in the literature are optimistic
whenever the true mechanism has a covariance-moderation component and
the design has off-drug phases, not merely a claim about one trial.
The prazosin-PTSD case is best kept as the motivating running example
it already is, not elevated to the paper's thesis.

**(c) Implications of framing (i).** Title, abstract, and introduction
are already well aligned with this framing and need no structural
change. The Section 4.1 extension applied during this review
strengthens framing (i) specifically, because it broadens the paper's
claim from a two-architecture dichotomy to a spectrum with identified
boundary regimes, which is a more defensible and more general
methodological caution than a strict either/or choice. Comparators
remain appropriate: Simon (2010), Freidlin and Korn (2014), and Kent
et al.'s PATH framework are the right anchors for a
methodology-in-simulation paper; no additional comparator class is
needed. Target journal: a methods-in-simulation-friendly outlet
(Statistics in Medicine, given the CSL already in use, or
Pharmaceutical Statistics) fits framing (i) better than a
disease-specific or software-focused venue.

**(d) Material to emphasize or de-emphasize.** Emphasize Section 3.2's
mechanistic explanation (mean blurring versus correlation erosion) and
the newly extended Section 4.1 boundary-regime discussion; these are
the paper's most generalizable content. De-emphasize, or move to a
supplement if length becomes a constraint, the prazosin-PTSD-specific
argumentation in the second half of Section 4.1 (the "three further
considerations" paragraph); it is well-argued but is application-
specific support for a claim (covariance moderation fits this one
trial) that is secondary to the paper's general methodological
contribution.

## 6. Assessment

**Verdict: minor revision.** The manuscript's central empirical claim
(Section 3.1, the architecture-by-carryover interaction, $z = 7.64$)
is solidly verified against committed data, the analytical framing is
correct, the literature positioning is accurate and appropriately
narrow (mean moderation dominates the cited literature; covariance
moderation is a genuine minority position with a clear biological
rationale), and the paper's own limitations section is honest about
what remains untested. The manuscript is publication-track modulo one
concrete, bounded fix: the Section 3.3 robustness check needs either a
regenerated data source or an explicit downgrade to "not yet
re-verified under the current N-matching convention." Everything else
identified in this review is either already resolved (as verified
during this pass), cosmetic (README staleness), or explicitly and
appropriately flagged by the authors themselves as future work.

## 7. Revision history

- 2026-08-26: First pub_review-branded whitepaper for this manuscript
  directory (prior review record for an earlier draft exists at
  `../referee-report-2026-06-13.md` under a different naming
  convention and predates the 2026-07-31 tutorial-style consolidation;
  its items are not individually carried forward here). Verified
  Section 3.1/3.4 numeric results against committed simulation output
  (`analysis/data/quick-sim/01-dgp-summary.txt`); found Section 3.3's
  N = 140 robustness check has no traceable generating script or data
  in the repository (Major issue 2.1); found two stale claims in
  `README.md` (Minor issues 3.1-3.2). At the author's request, applied
  one remediation during this review pass rather than only
  recommending it: extended Section 4.1 with material on boundary
  regimes and trial-design mechanism drawn from the companion
  latent-class manuscript (`03-latent-class-mixture-application`),
  re-rendered the manuscript successfully (17 pages, no errors), and
  recorded the change under "What remains to be done," item (b)4, for
  a final author read-through.
