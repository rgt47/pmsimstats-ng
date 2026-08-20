# White paper: referee-style review of the pmsimstats-ng manuscript compendium

*Review date: 2026-08-20 10:17 PDT*

**Scope.** This review covers all eleven numbered manuscript series
under `analysis/report/` (01 through 11), the three package
vignettes under `vignettes/`, and the cross-cutting infrastructure
documents (`analysis/report/README.md`, `NOTATION.md`,
`whitepaper-notation-audit.md`, `whitepaper-compendium-summary.md`).
For each series the full `report.Rmd` was treated as the primary
document under review, with `report-slim`, `report-extra-slim`,
`report-mini`, `report_short`, and `bullets.Rmd` treated as
alternates and not separately reviewed line by line.

**Method and epistemic status.** All eleven `report.Rmd` primary
manuscripts, the three package vignettes, and the cross-cutting
infrastructure documents were read in full over the course of this
review, split across four parallel full-text passes (series 01-03;
04-06; 07-09; 10-11 plus vignettes and cross-cutting docs) that were
then reconciled into this single document. Findings from a delegated
full-text read of a manuscript this reviewer did not personally open
are marked **inspected (delegated)** below, to distinguish them from
findings this reviewer verified by direct file access in the final
reconciliation pass (**inspected**) and from anything recomputed by
running code (**verified**). All four passes were also checked
against the repository's own internal referee document,
`analysis/report/whitepaper-compendium-summary.md` (dated
2026-07-29, updated 2026-08-12), an unusually candid and detailed
self-assessment, and, where present, the series-specific internal
`referee-report-*.md` files (found for series 01, 06 [four dated
reports], 07, and 08). Placeholder-marker counts
(`Lorem ipsum` / `To be completed by rgt` / `To be finalized` /
`Forthcoming`) and word counts were independently grepped and
matched the internal compendium's claimed counts for every series to
within one line (**verified**, 2026-08-16). No simulation in this
repository was re-executed for this review, so no numerical result
reported below was independently recomputed; all quantitative claims
are as stated in the manuscripts or the internal compendium. This
reconciliation pass identified a material gap between the internal
compendium's characterization of series 04, 05, 06, 07, and 09 and
what direct full-text reading found: several defects more serious
than anything the internal compendium documents are reported for the
first time in Major Issues 2.9-2.13 below.

**Method note for this update (2026-08-20).** This is an update pass,
not a fresh review. Before re-reading, `git log --since="2026-08-16"`
and `git diff --stat` were run against the full repository
(**verified**, 2026-08-20) to establish which files had changed since
the initial review. The result: every commit since 2026-08-16
(`fbe6235` and a series of auto-backups) touches only
`analysis/report/02-carryover-sensitivity/` and its supporting
`analysis/scripts/carryover-sensitivity/`, `analysis/data/02-*`, and
`analysis/figures/02xs-*` paths, plus the provenance-stamped copies
under `analysis/report/share/`. No other series' `report.Rmd`, no
vignette, `DESCRIPTION`, `NOTATION.md`,
`whitepaper-notation-audit.md`, or `whitepaper-compendium-summary.md`
has changed a single byte. On that verified basis, every finding below
concerning series 01, 03-11, the vignettes, and the cross-cutting
infrastructure documents is carried forward unchanged from the
2026-08-16 review without re-reading the underlying files line by
line; this review's confidence in those findings rests on the
2026-08-16 full-text and delegated passes described above, corroborated
now by direct confirmation that nothing in scope for those findings
could have changed. Series 02
(`analysis/report/02-carryover-sensitivity/`) was substantially
revised in the interim (a nine-specification G1-G9 extension of
Figures 1-2, an Architecture A exploratory appendix, a
report/supplement split, and an all-simulations-to-500-replicates
precision pass) and was re-read in full for this update, including
the new `supplement.Rmd` and the newly created `bullets.Rmd` outline
document; see the revised series 02 summary in Section 1, new Major
Issue 2.14, and new Minor Issues 3.11-3.12 below.

**A note on this review's relationship to the internal compendium.**
The single most important review artifact already in this repository
is the authors' own `whitepaper-compendium-summary.md`. It is
methodologically sound, self-critical, dated, and tracks its own
prior recommendations against subsequent commits with unusual
discipline (a "Status: adopted / not adopted" ledger). This review
does not re-litigate the fifty-odd specific issues that document
already identifies correctly; it corroborates them, adds the issues
that document's scope excluded (the vignettes, package dependency
declarations, the notation-audit stub, and a formal submission-unit
partitioning grounded in target-journal fit), and renders an
independent referee verdict per candidate paper rather than a
project-management status report.

## 1. Summary of the work under review

**01. Mean moderation and covariance moderation as data-generating
architectures for biomarker-treatment interactions.** Establishes
that a biomarker-treatment interaction can be encoded into a
simulated N-of-1 trial's data-generating process (DGP) two
structurally different ways, by scaling the conditional mean
directly (mean-moderation, "Architecture A") or by inducing
treatment-state-dependent correlation in the joint covariance
(covariance-moderation, "Architecture B", the convention inherited
from Hendrickson et al. 2020), and that the two are not
interchangeable: under a one-week carryover half-life the
covariance-moderation architecture loses roughly ten times more
relative power (30.6% vs. 2.8%) in the open-label-plus-blinded-
discontinuation design than the mean-moderation architecture, a
difference-of-differences significant at $z = 7.64$. **Inferred**
from the internal compendium and spot-checked by direct reading:
this is the compendium's keystone paper, now free of placeholder
prose and of the Architecture C material it previously duplicated
with series 11.

**02. Carryover-mitigation analysis strategy.** A negative result:
of three ways to encode drug exposure as a carryover-robust
predictor in the analysis model, the exposure-weighted
specification's advantage over an unadjusted baseline is
architecture- and design-conditional, and is actively inferior under
the classical crossover design (power 0.488 vs. 0.830).
**Inspected in full by this reviewer, 2026-08-20** (updated from the
2026-08-16 **inferred** status): since the last review the manuscript
was substantially revised. Figures 1-2 were expanded from the
original three specifications to all nine (adding AIC-selected
half-life and paired-difference specifications, each with and without
a CR2 cluster-robust variant); a same-directory `supplement.Rmd`
companion was split off, moving full cell-by-cell detail for several
sensitivity blocks out of the main text and adding an explicitly
exploratory Architecture A counterpart to Figures 1-3 (Section S12);
and every simulation feeding either document was escalated to 500
replicates per cell, replacing several that had stood at 100 or 250.
This review re-rendered both documents (**verified**, clean build, no
errors) and spot-checked a sample of the newly reported power values
against the underlying `.rds` files directly (**verified**), finding
the reported numbers consistent with the data to two decimal places.
The manuscript's own decay-shape claim is now correctly qualified
(the Exposure-weighted advantage "narrows to statistical
insignificance under the two most heavy-tailed Weibull shapes
examined," Abstract and Section 3.5) rather than overclaimed, and
contains no placeholder markers (**verified** by grep). However, this
update surfaced a new defect: the series' own `cover-letter.md` has
not been revised to match, and now directly contradicts the
manuscript's corrected finding (see new Major Issue 2.14).

**03. Latent-class and mixture formulations of biomarker
moderation.** A theory-plus-pilot paper. The theoretical contribution
(when a single multivariate-normal DGP can substitute for a
latent-class mixture, and a bound on achievable marginal correlation
by the geometric mean of class separations) is delivered in full;
the empirical contribution is a 240-replicate pilot of the last of
five originally pre-registered studies, correctly and explicitly
labeled in the current draft as a pilot rather than a production
result. **Inferred.**

**04. N-of-1 versus parallel-group RCT for the treatment main
effect.** Reports that an aggregated N-of-1 hybrid design dominates a
matched-$N$ parallel-group RCT by roughly threefold in power to
detect the treatment main effect, in `report.Rmd`'s primary pipeline
(N-of-1 power 1.000 vs. RCT 0.675 at the primary effect size).
**Inspected (delegated)**: the alternate `report_short.Rmd` reports
materially different numbers for the identical nominal comparison
(0.989 vs. 0.189), drawn from a different, non-reconciled simulation
pipeline that is not reproducible as documented (its stated driver
script does not exist in the repository). See Major Issue 2.9. Also
carries the compendium's previously flagged design-vs-observation-
count confound (Major Issue 2.6).

**05. Twelve-sweep design sensitivity.** A one-factor-at-a-time
sensitivity analysis of twelve design axes for the treatment main
effect (not the interaction), concluding the Hendrickson hybrid
design is a competent default and that the parallel-RCT ordering
only overtakes it once carryover exceeds the within-subject period
length. **Inspected (delegated)**: the manuscript's own Results
section is explicitly disclosed, in a boxed note the authors wrote
themselves, as stale — generated at a discarded $N=40$ baseline
rather than the paper's stated $N=35$ baseline, with re-running
described as "in progress." Every power/bias/coverage number quoted
in Results is therefore known-stale by the authors' own account. See
Major Issue 2.10.

**06. Three-component decomposition (BR/PB/TV).** Derives an
omitted-variable-bias identity showing the one-component biomarker
slope is displaced from the true pharmacological slope by a term
proportional to the biomarker's covariance with the placebo-belief
and natural-history components, and shows by simulation that this
bias is negligible for an uncontaminated biomarker but material
(+0.51) for a contaminated one, recoverable only under a
balanced-placebo design. **Inferred** for Study A (the original,
readRDS-backed pilot); **inspected (delegated)** for the newer Study
B subsection, whose headline recovery numbers are not wired to any
`readRDS()` call in the current setup chunk — a recurrence of a
defect a prior internal referee had already flagged. See Major Issue
2.11. The longest source file in the compendium at roughly 18,200
words as measured directly by this review (**verified** by `wc -w`
against `report.Rmd`, 2026-08-16; the internal compendium's estimate
of "26,000 words" refers to an earlier or differently measured state
and was not reconciled by this review).

**07. Gompertz trajectory evaluation.** Tests whether the choice of
symptom-trajectory shape (Gompertz vs. three alternative growth
families) changes detectability of the interaction, finding the
choice is exactly inert under covariance-moderation (matched seeds
return identical decisions) and separates the families by 0.039 in
power under mean-moderation. **Inspected (delegated)**: this
current-draft finding directly contradicts `cover-letter.md`, which
still asserts the pre-fix null result ("families produce
statistically indistinguishable performance... under both
architectures"); and the identifiability/"sloppiness" diagnostic the
Introduction frames as a pillar contribution is never delivered in
Results or Discussion. See Major Issue 2.12. A narrow, single-design,
zero-carryover, 16-cell study.

**08. Test procedure and design.** Separates the cost of
dichotomizing a continuous biomarker (power loss 0.08-0.21) from the
further cost of using RM-ANOVA or naive GEE instead of a mixed model
(at most 0.04 more), and shows naive GEE inflates Type I error to
1.5-1.9x nominal, repaired by the Mancl-DeRouen correction.
**Inspected (delegated)**: this is the most improved manuscript in
the compendium relative to its own 2026-06-17 referee report — six
of seven flagged issues, including addition of a continuous-
biomarker contrast arm that cleanly isolates dichotomization cost
from test-procedure cost, appear resolved on direct comparison. The
cycle-by-period design grid promised in the title remains specified
but deferred.

**09. Informative dropout by design.** Argues that dropout rate and
mechanism are coupled to trial design family in a way the N-of-1
methods literature has not examined, and reports (at a matched
$N = 70$) power ordering hybrid (0.968) > OL+BDC (0.950) > crossover
(0.854) > open-label (0.088), with hybrid and OL+BDC statistically
inseparable ($z = 1.44$). **Inspected (delegated)**: this is the
least mature manuscript in the compendium — it has not been through
any internal referee cycle, contains a literal unedited placeholder
in running Simulation Design prose ("(Forthcoming detail...)"), ran
only 16 of a pre-registered 372 cells with no Type I error cells
executed for three of four designs, and its flagship explanatory
mechanism (a randomization-path decomposition) is asserted as
"consistent with" the marginal results rather than formally tested,
by the manuscript's own admission. See Major Issue 2.13.

**10. Type I calibration of fixed-effect interaction tests in linear
mixed models.** **Inspected in full by this reviewer.** A
self-contained methodological paper, not specific to N-of-1 trials:
a four-channel decomposition of a Wald statistic's miscalibration
(point-estimate bias, standard-error scale, degrees of freedom,
reference-distribution shape), applied to the project's own
`corCAR1` mixed model, localizes a 6-10% standard-error
over-estimation (Type I error approximately 0.03 against nominal
0.05) to the `corCAR1` residual-correlation layer, shows the bias is
asymptotic (flat from $N = 35$ to $N = 280$, unaffected by a
Kenward-Roger correction), and repairs it with a CR2 bias-reduced
cluster-robust standard error, cross-checked against a
Mancl-DeRouen-corrected GEE. This is the technically strongest
manuscript in the compendium: every table in the Results section
reads a named `.rds` file produced by a named, present driver script
(**verified**, 2026-08-16, by confirming the seven scripts and seven
`.rds` outputs exist), the derivations in Section 2 are internally
consistent and correctly cited, and the finding is genuinely
surprising (a *conservative*, not anti-conservative, test, with the
sandwich correction acting in the unfamiliar direction).

**11. A combined data-generating architecture (channel
super-additivity).** **Inspected in full by this reviewer.** Defines
a "dual-channel" DGP that activates both the mean-moderation and
covariance-moderation mechanisms simultaneously, shows analytically
that the biomarker-treatment signal is additive across channels
while power is super-additive because the power function is convex
in the noncentrality below its inflection, and reports a $3\times3$
channel-weight-by-design-by-carryover simulation grid consistent
with that prediction. The manuscript is explicit, in its own
abstract, that it is unfinished: "[To be finalized after the
common-driver boundary re-run]," and Section 4.2 states the
reported grid uses a sample-size convention that does not match the
companion paper's (01's) production convention, so the boundary
cells (where one channel weight is zero) have not been shown to
reproduce 01's power figures numerically, only algebraically. This
is verifiable directly from the source (**inspected**) and is not
resolved as of this review.

**Vignettes (`vignettes/01-`, `02-`, `03-`).** **Inspected in full
by this reviewer.** Package-level demonstration material, not
manuscripts, walking a user through simulating trial data, plotting
power results, and analyzing a real (or realistic) clinical trial
data set with the package's exported functions. All three still
address an unnamed target publication ("the results published in
`[___]`") and contain code defects detailed under Major Issues 6 and
7 below that would prevent them from building non-interactively.

## 2. Major issues

**2.1. Series 11's headline result is not yet validated against its
own stated correctness gate, and its abstract says so.**
*Location:* `analysis/report/11-combined-dgp-architecture/report.Rmd`,
Abstract and Section 4.2 ("Boundary validation against the companion
baselines"), lines 66-67 and 323-347. **Inspected.**
*Problem:* Section 2.2 states explicitly that the boundary
reductions ($c_{bm,B}=0$ recovers pure mean-moderation,
$c_{bm,A}=0$ recovers pure covariance-moderation) are "internal
validation gates" that "a correctly implemented combined driver must
reproduce... numerically." Section 4.2 then reports that this gate
has not been passed: an earlier combined driver used a
total-sample-size convention inconsistent with 01's per-path
convention, so the boundary cells are "uniformly lower than, and not
directly comparable to, the companion baselines," and the reported
Table 1 is explicitly flagged as using that unreconciled convention.
The Discussion (Section 6) nonetheless asserts the super-additivity
finding as established, and the Conclusions of the companion
compendium document describe this as the compendium's central
structural claim about DGP architecture propagating into every other
paper's conditioning. *Why it matters:* a referee cannot accept an
empirical claim whose own correctness gate the manuscript states has
not been passed. The 3x3 power table (Table 1) is not currently
interpretable as evidence for anything beyond "the driver that
produced this table behaves qualitatively as predicted"; it cannot
yet be read as confirming the analytical prediction quantitatively,
because the analytical prediction is stated relative to the
pure-architecture baselines this table has not been shown to match.
*Remediation:* re-run the combined driver at the per-randomization-
path sample-size convention used by 01's production driver
(`analysis/scripts/quick-sim/01-dgp-prototype.R`, per the
Reproducibility section), confirm the two boundary rows/columns of
Table 1 reproduce 01's published power values within the stated
Monte Carlo tolerance (approximately 0.013-0.016), and only then
finalize the abstract's Results paragraph, which currently reads as
a placeholder pending exactly this step.

**2.2. Vignette `01-simulate-trial-data.Rmd` contains an interactive
debugger call and is not built or checked by the package's test
infrastructure.**
*Location:* `vignettes/01-simulate-trial-data.Rmd`, line 367 (inside
the `if(rerun_simulations)` block, immediately before the first call
to `generateSimulatedResults`). **Inspected.**
*Problem:* the chunk contains a bare `browser()` call. If this
chunk is ever executed non-interactively, `browser()` in a
non-interactive R session either errors or (depending on `options`)
silently continues without pausing, but its presence in
version-controlled, ostensibly executable vignette source is itself
evidence the code has not been run end-to-end in the state it is
committed in; more consequentially, `rerun_simulations` is set to
`TRUE` (not `FALSE`, contrary to the surrounding prose, "there's a
flag that keeps this code from actually running"), so a `vignette
build`, `R CMD check --as-cran`, or `devtools::build_vignettes()`
would attempt to execute this block, hit `browser()`, and hang or
error in a non-interactive context. *Why it matters:* this is a
CRAN-submission blocker (a hung or erroring vignette build fails
`R CMD check`), and it means the vignette's claim to demonstrate a
working pipeline has not been verified by execution since the last
edit that introduced or left in the debugger call. *Remediation:*
remove the `browser()` call; set `rerun_simulations <- FALSE` by
default (matching the prose's stated intent) or gate execution on
`eval = FALSE` / a `Sys.getenv()` flag consistent with how long-
running vignettes are usually handled in this ecosystem; add
`vignettes/*.Rmd` to whatever CI or pre-commit check already governs
`R CMD check` for this package, since this defect would have been
caught by a single build attempt.

**2.3. Vignette `02-visualize-power-results.Rmd` plots the wrong
data object due to a variable-name typo.**
*Location:* `vignettes/02-visualize-power-results.Rmd`, lines
140-146. **Inspected.**
*Problem:* `data(results_maxes)` is loaded and assigned to
`simrsults` (misspelled, missing the "e"), but the subsequent
`PlotModelingResults(simresults, ...)` call three lines later
references `simresults` (correctly spelled), which at that point in
the script still holds the previous chunk's `results_core` object.
Panel "A. Effect of response parameter maximum values on power" is
therefore silently rendered from the wrong data set; if the vignette
is ever run and inspected only for "does it produce a plot without
erroring," this defect is invisible, because both objects are valid
`data.table`-backed simulation-result structures with compatible
column sets. *Why it matters:* a reader following this vignette to
learn how response-parameter maxima affect power will be shown a
plot of something else, silently. This is exactly the class of bug
static inspection and prose review will not catch and only
execution-plus-output-comparison would. *Remediation:* fix the typo
(`simrsults` to `simresults`, or vice versa, whichever is intended);
add a spot check (even a simple `stopifnot(identical(...))` or a
`testthat`/`tinytest` snapshot of vignette output) that would have
caught this; audit the sibling vignette 01 for the same class of
copy-paste variable-reuse error, since its parameter-set list
(`respparamsets`, `respparamsetsm`, `origrespparamsets`, `TRgrid`
reused across two unrelated construction blocks) is built with the
same undisciplined naming pattern.

**2.4. Packages loaded by all three vignettes are not declared in
`DESCRIPTION`.**
*Location:* `DESCRIPTION` (`Imports:`/`Suggests:` fields) versus the
`setup` chunks of all three vignettes
(`vignettes/01-simulate-trial-data.Rmd` lines 18-28,
`vignettes/02-visualize-power-results.Rmd` lines 18-28,
`vignettes/03-analyze-clinical-trial-data.Rmd` lines 18-29).
**Verified**, 2026-08-16, by reading `DESCRIPTION` directly.
*Problem:* all three vignettes call `library(data.table)`,
`library(corpcor)`, and `library(MASS)`; vignette 03 additionally
calls `library(merTools)`; and `ggplot2` is invoked directly
(`library(ggplot2)`) rather than only through `ggpubr`. None of
`data.table`, `corpcor`, `MASS`, `merTools`, or `ggplot2` appears in
`DESCRIPTION`'s `Imports:` or `Suggests:` fields (only `nlme`,
`lme4`, `lmerTest`, `furrr`, `future`, `svMisc`, `tictoc`,
`ggpubr`, `gridExtra` are declared as Imports, and
`tinytest`, `knitr`, `rmarkdown`, `dplyr`, `tidyverse`, `nof1power`
as Suggests). *Why it matters:* CRAN policy (and `R CMD check
--as-cran`) requires every package used by a vignette, including via
a bare `library()` call, to be declared, normally in `Suggests`. An
undeclared dependency will fail `R CMD check` on a machine that does
not happen to have `data.table`, `corpcor`, `MASS`, or `merTools`
already installed as a side effect of some other package's
dependency chain (MASS ships with base R, softening but not
eliminating the problem for the other three). This is a submission-
readiness defect independent of, and easier to fix than, the code
bugs above. *Remediation:* audit every vignette (and every example
in `man/`) for `library()`/`::` calls against `DESCRIPTION`, and add
missing packages to `Suggests`.

**2.5. `whitepaper-notation-audit.md` is an empty stub, but
`NOTATION.md` and `analysis/report/README.md` both cite it as the
authority behind the canonical notation.**
*Location:* `analysis/report/whitepaper-notation-audit.md` (7 lines,
all YAML front matter, no body content); referenced from
`analysis/report/NOTATION.md` line 13 ("The audit behind it is
`analysis/report/whitepaper-notation-audit.md`, whose PDF binds this
file in as its appendix") and from `analysis/report/README.md`.
**Verified**, 2026-08-16, by reading the file directly (`wc -l`
confirms 7 lines, all preamble). *Problem:* a document the project's
own infrastructure describes as load-bearing ("the audit that
produced" the canonical notation table used by all eleven
manuscripts) contains no audit. Either the audit was never written,
or it was written and lost, or the compiled PDF
(`whitepaper-notation-audit.pdf`) contains content that was never
committed back to the Markdown source, in which case the two have
already drifted, which is precisely the failure mode `NOTATION.md`
says the appendix binding exists to prevent. *Why it matters:* for
a submission where notation consistency across a family of related
papers is itself part of the argument for treating them as a
coherent programme (see Section 5, Recommended Framing), a claimed
audit trail that turns out to be empty undermines that argument the
moment a referee or careful reader checks it. *Remediation:* either
write the audit (a short document cross-referencing NOTATION.md's
symbol table against each manuscript's actual usage, exactly as its
own description promises) or remove the claim that it exists from
`NOTATION.md` and `README.md` and cite the rendered PDF's provenance
some other way.

**2.6. Series 04's central comparison confounds trial design with
observation count, and its Limitations section does not disclose
this.** *Location:*
`analysis/report/04-treatment-main-effect/report.Rmd` (per the
internal compendium's assessment, **inferred**, not independently
re-read by this reviewer beyond the compendium's citation of the
Limitations section's content). *Problem:* both arms are matched at
$N = 35$ participants, but the N-of-1 arm contributes eight
within-participant periods against the parallel-RCT arm's single
post-baseline measurement, so a substantial share of the reported
threefold power advantage plausibly reflects "more measurements per
subject" rather than a property of the N-of-1 design per se. The
Limitations section, per the compendium's direct inspection,
discusses compliance, single-outcome restriction, and carryover-
decay-form uncertainty, but not observation count or participant
burden. *Why it matters:* this is exactly the objection a referee at
a clinical or general biostatistics venue will raise first, because
"more repeated measures yields more power" is not a novel finding
and the paper's clinical framing (N-of-1 as a superior design
choice) depends on the comparison being design-fair. *Remediation:*
add a measurement-matched or total-burden-matched RCT sensitivity
arm, or state the asymmetry explicitly in the Methods and defend it
as the decision-relevant comparison for a fixed recruitable
population (which is a legitimate position, but currently
unargued).

**2.7. Series 11 duplicates conceptual content already claimed by
series 03, and the division of labor between the two is stated but
not tested for redundancy from a reader's perspective.**
*Location:* `analysis/report/11-combined-dgp-architecture/report.Rmd`
Section 2.2 (lines 163-192) and Section 6 (lines 470-493), citing
`@pmsimstats-paper03` twice as having already "introduce[d] it
conceptually as the general node of a biomarker-moderation spectrum
and suppl[ied] its psychometric lineage." **Inspected.** *Problem:*
series 11's own Discussion states its contribution is "complementary
and specific" relative to 03, namely the analytical super-additivity
derivation and the calibrated grid, "neither of which the conceptual
treatment provides." That is a defensible division of labor in
principle, but it means a reader of both papers encounters the same
Arminger-et-al. mean-and-covariance-structure framing, the same
latent-class-spectrum framing, and the same biological motivating
argument (a biomarker indexing a responder subtype shifts both
class-conditional mean and dispersion) twice, once as exposition in
03 and once as background in 11. Combined with 11's brevity (roughly
3,000 words by direct count) and its incomplete validation (Major
Issue 2.1), this raises the question the internal compendium already
raised and the author did not adopt (folding 11 into 01): whether 11
earns independent-manuscript status at all, as opposed to being a
sensitivity subsection of 01 or a worked example inside 03. *Why it
matters:* a compendium submitted as multiple papers to the same
research area invites a referee who has seen a sibling submission to
flag self-plagiarism of framing, even where the technical content
differs. *Remediation:* see Section 5 below (Recommended Framing);
at minimum, if 11 is to remain standalone, its Introduction should
argue explicitly for independent-manuscript status rather than
relying on a "complementary and specific" claim made only in the
Discussion, after the reader has already encountered the repeated
framing.

**2.8. The compendium-wide `rgt`/`bullets`/`orig` three-part
scaffold leaves eight of eleven manuscripts with unwritten prose
inside submitted-looking LaTeX/Rmd structure.** *Location:* series
04 (110 placeholder markers, verified by this review's grep), 06
(101), 08 (59), 05 (54), 07 (51), 09 (51), 10 (37), 11 (13); series
01, 02, 03 have zero. **Verified**, 2026-08-16, by direct grep
against every `report.Rmd`, matching the internal compendium's
counts to within one. *Problem:* every `.rgt` div in every affected
manuscript currently contains the literal string "Lorem ipsum dolor
sit amet, consectetur adipiscing elit..." (or, for 11, "To be
completed by rgt."), sitting immediately adjacent, in the rendered
PDF, to fully written `.orig` prose making the paper's actual
argument. This is confirmed by direct inspection of series 10 and
11's source, where every section alternates a placeholder block and
a real one. *Why it matters:* this is not cosmetic. A `.rgt` div
that renders in the PDF (as opposed to being suppressed by the LaTeX
preamble) would place visible Lorem Ipsum text in a submitted
manuscript; even if the PDF build suppresses `.rgt` output entirely
(**unverified** by this review; the `sim-preamble.tex` mechanics
were not inspected), the source repository itself is not in a
submittable state for eight of eleven series, and the placeholder
count is the single most direct, mechanically checkable proxy for
"how far from done" each series is. *Remediation:* as the internal
compendium already recommends and tracks, write the `rgt` prose
series by series in priority order; this review concurs with that
document's prioritization (01, 02, 03 done; 10 next, since it is the
compendium's only non-N-of-1-specific, most portable contribution;
then descending by placeholder count for the rest), with one
adjustment: 11 should not receive `rgt`-writing effort until Major
Issue 2.1 (the boundary-validation re-run) is resolved, since writing
polished prose around an unvalidated result is effort that may need
to be redone or retracted.

**2.9. Series 04 reports two irreconcilable primary results for the
same headline comparison, and one of the two pipelines is not
reproducible.** *Location:*
`analysis/report/04-treatment-main-effect/report.Rmd` versus
`analysis/report/04-treatment-main-effect/report_short.Rmd`.
**Inspected (delegated).** *Problem:* at the primary effect size,
`report.Rmd` (n_sim = 2,000) reports N-of-1 power 1.000 against RCT
power 0.675; `report_short.Rmd` (n_sim = 1,500), addressing the
identical nominal comparison, reports 0.989 against 0.189. The two
documents load different precomputed objects:
`analysis/data/derived_data/04-mainfx-production.rds` (dated
2026-05-19, produced by a present `production-run.R`) for
`report.Rmd`, and `analysis/data/derived_data/sim_workspace.RData`
(dated 2026-05-07, older) for `report_short.Rmd`, the latter
documented as produced via a `tools/run_sim_save.R` script that does
not exist anywhere in the repository. *Why it matters:* this is not
a stylistic variant of the same finding; it is two different
numeric answers to the paper's central question, and the pub_review
instructions treat `report_short.Rmd` as a reviewable alternate, not
an archived draft. A referee presented with both, or a reader who
discovers the discrepancy after the primary manuscript is accepted,
would treat this as disqualifying until reconciled. *Remediation:*
either delete or explicitly mark `report_short.Rmd` as superseded
and non-authoritative (with a note explaining the RCT power
discrepancy is a stale-pipeline artifact, not a robustness check),
or reconcile the two pipelines and re-derive both from the same,
currently-working driver before either is cited anywhere.

**2.10. Series 05's entire Results section is self-disclosed as
generated from a discarded, non-canonical baseline.** *Location:*
`analysis/report/05-nof1-design-sensitivity/report.Rmd`,
"Reproducibility and execution" subsection, approximately line 392.
**Inspected (delegated).** *Problem:* a boxed author note states
verbatim that every quoted power/bias/coverage number in the
Results section (spanning all twelve S1-S12 sweeps) was computed at
$N=40$, while the manuscript's stated baseline sample size is
$N=35$, and that re-running the sweeps at the correct baseline "is
in progress." Separately, a "Prototype Monte Carlo confirmation"
subsection (approximately line 1229) runs only 500 replicates per
cell and states in its own text that this "does not bound the
80%-power region, which is the central deliverable a production run
would need to provide." *Why it matters:* this is the single
clearest author-acknowledged blocker to submission found anywhere in
the compendium. A manuscript cannot be evaluated, let alone accepted,
while its own authors state its Results section does not reflect its
stated design. *Remediation:* re-run all twelve sweeps at $N=35$
before any further revision effort (including `rgt` prose-writing,
Checklist item 5) is spent on this series; do not cite any numeric
result from the current Results section in any other manuscript
until this is done (note series 08 defers its own cycle-by-period
grid to "a companion paper" without naming it — confirm this is not
meant to be 05, or reconcile if it is).

**2.11. Series 06's Study B recovery numbers are not traceable to
any loaded data object, repeating a defect a prior internal referee
already flagged for the paper.** *Location:*
`analysis/report/06-component-decomposition/report.Rmd`, Study B
subsection ("recovery under a belief-decoupling design"); compare
`referee-report-2026-06-15.md`, item R3. **Inspected (delegated).**
*Problem:* the setup chunk contains exactly two `readRDS()` calls
(for `sim_summary` and `sim_contam`, both associated with the
original, Study A/pilot results); the newly added Study B
subsection reports specific bias figures ($-0.022$ to $+0.048$) and
coverage figures (0.94-0.96 vs. 0.86) as hardcoded prose, with no
`readRDS()` call for the stated output path
(`analysis/data/derived_data/component-decomposition/study-b-
balanced/`). The 2026-06-15 referee report had already flagged
exactly this pattern for the paper's original numbers ("a live drift
hazard... the Study A numbers had already silently changed once")
and recommended wiring every reported number to a loaded object.
*Why it matters:* Study B is the section that resolves the paper's
central prescriptive claim (that a balanced-placebo design recovers
an unbiased slope under contamination); a referee cannot verify a
hardcoded number, and the paper's own revision history shows this
exact failure mode has already produced a silent numeric drift once.
*Remediation:* add the missing `readRDS()` call (or inline
computation) for the Study B output object and replace every
hardcoded Study B number with an inline-computed value, exactly as
the 2026-06-15 referee already prescribed for the rest of the paper.

**2.12. Series 07's cover letter contradicts the current
manuscript's headline finding, and a promised diagnostic is never
delivered.** *Location:* `analysis/report/07-gompertz-
evaluation/cover-letter.md` versus `report.Rmd` Results/Discussion;
`report.Rmd` Introduction (identifiability/"sloppiness" framing,
citing Jagadeesan 2023 and Tjorve 2017) versus Results/Discussion/
Conclusions. **Inspected (delegated).** *Problem:* the current
manuscript reports a real, if modest, 0.039 power spread across
trajectory families under mean-moderation architecture; the cover
letter still states that "families produce statistically
indistinguishable performance... under both architectures," the
pre-fix headline from before an implementation defect (a
subject-invariant additive shift baked into the original family
manipulation, previously flagged by `referee-report-2026-06-13.md`
item M1) was corrected. Separately, the Introduction repeatedly
frames trajectory-family "sloppiness" and structural
identifiability as a pillar contribution, but no eigenvalue-spectrum
or identifiability result appears anywhere in Results, Discussion,
or Conclusions (four textual mentions, all confined to framing).
*Why it matters:* an editor or referee reading the cover letter and
the manuscript together would see a direct factual contradiction on
the paper's own headline result; a promised analytical contribution
that never appears is a completeness defect independent of that
contradiction. *Remediation:* rewrite `cover-letter.md` to match the
current, corrected finding; either deliver the sloppiness/
identifiability diagnostic or remove it from the paper's framed
contributions.

**2.13. Series 09 is the least mature manuscript in the compendium:
it contains an unedited content placeholder in running prose, has
not been through any referee cycle, and its flagship explanatory
mechanism is asserted rather than tested.** *Location:*
`analysis/report/09-informative-dropout-by-design/report.Rmd`,
Simulation Design section (literal string "(Forthcoming detail.
Pmsimstats infrastructure with the existing R/censordata.R
hazard-based dropout module...)"); Discussion/Conclusions
("Prototype Monte Carlo confirmation," "production-grade
extensions... represent the natural follow-up to the present
proof-of-concept evidence"). **Inspected (delegated).** *Problem:*
unlike series 01, 06, 07, and 08, no `referee-report-*.md` exists
for series 09. Only 16 of a pre-registered 372 cells
($180+120+12+60$, across Studies 1-4) were executed; no Type I
error cells were run at all for OL+BDC, crossover, or hybrid designs
(self-acknowledged in-text); the randomization-path decomposition
(Study 3), the mechanism the paper repeatedly invokes to explain why
hybrid and OL+BDC are dropout-robust, is not delivered as a formal
test, only asserted as "consistent with" the marginal power results;
and only 500 replicates/cell were run against a pre-registered
1,000-5,000, with the manuscript itself stating that resolving its
own central biased-vs-MCAR-dropout question "would require
approximately 2,000 replicates per cell." *Why it matters:* the
manuscript is formatted, abstracted, and concluded with the same
authority as the compendium's more mature papers, but its own text
concedes the central mechanism is unconfirmed and the evidentiary
base for several claimed comparisons (Type I error across three of
four designs) does not exist yet. A referee would not distinguish
this from a paper simply overclaiming its pilot results. On the
positive side, the paper is commendably transparent about and
corrects a genuine earlier internal error (a prior unmatched-$N$
draft that had spuriously concluded the most powerful design was
also the most dropout-robust), a stronger self-correction signal
than series 07 showed with its stale cover letter. *Remediation:*
replace the literal placeholder with actual Simulation Design text;
run the Type I error cells for the three missing designs before any
design-ordering claim is asserted as general; either deliver the
randomization-path decomposition as a formal test or demote it from
"the mechanism" to "a candidate explanation, not yet tested" in the
Discussion; route this manuscript through the same internal referee
cycle 06, 07, and 08 have already had before further polish is
invested.

**2.14. Series 02's `cover-letter.md` was not updated alongside the
manuscript's decay-shape robustness finding and now asserts the
overclaim the current manuscript explicitly corrects; it also
describes a specification scheme, replicate count, and cell count the
manuscript has since superseded.** *Location:*
`analysis/report/02-carryover-sensitivity/cover-letter.md`, paragraphs
2-3, versus `report.Rmd` Abstract and Section 3.5. **Inspected**,
2026-08-20. *Problem:* the cover letter states the exposure-weighted
specification's power advantage holds "with the ranking robust to
decay-form mis-specification," unqualified. The current manuscript's
Abstract states the opposite for part of the tested range: the
advantage "widens under a light-tailed decay-shape mis-specification,
but narrows to statistical insignificance under the two most
heavy-tailed Weibull shapes examined" (gaps of 0.00 to 0.02, against
an MCSE of approximately 0.022, at the most extreme shapes tested,
per Section 3.5). This is not a hedge the cover letter merely omits;
it is the direct opposite of what the letter claims as a selling
point of the submission. The cover letter also frames the design
around three analysis specifications under the retired `A1`/`A2`/`A3`
labels (the manuscript now reports nine, under the `G1`-`G9` display
codes), cites a "24-cell expansion" at 600 replicates that no longer
matches the manuscript's current all-500-replicates precision
standard, and dates to 2026-05-09, predating the 2026-07-30
consolidation the series' own `README.md` documents. *Why it
matters:* this is the same failure mode already identified for series
07 (Major Issue 2.12): a submission document that reaches an editor
alongside the manuscript directly contradicts the manuscript's
headline claim on the exact point (decay-form robustness) a referee
who has read both would be most likely to check first, since it is
stated as a strength in the cover letter and immediately qualified in
the Abstract. Finding this defect in a second, previously
well-regarded manuscript (series 02 carried no Major Issue at the
2026-08-16 review) indicates the pattern is not isolated to series 07
and is a symptom of a missing process step (cover letters are not
revised in the same pass as the manuscript body) rather than a
one-off oversight. *Remediation:* rewrite `cover-letter.md` to state
the qualified decay-shape finding, update the specification count and
labels to G1-G9, update the replicate-count and cell-count claims to
match the current build, and re-date it; as a compendium-wide
process fix, add a step to whatever pre-submission checklist governs
this repository requiring `cover-letter.md` to be re-read against the
current Abstract in the same pass that finalizes any manuscript's
Results section.

## 3. Minor issues

**3.1.** *Location:* `analysis/report/01-dgp-mean-moderation-vs-mvn/README.md`.
**Inferred** from the internal compendium. The README still
describes the paper's `rgt` layer as carrying unfilled placeholders,
which is no longer true of `report.Rmd` (zero placeholder markers,
verified by this review's grep); a reader consulting only the
README would undersell a manuscript that has, on this dimension,
reached submission readiness. *Remediation:* refresh the README in
the same pass that finalizes the paper.

**3.2.** *Location:*
`analysis/report/10-interaction-test-calibration/` and
`analysis/report/11-combined-dgp-architecture/`. **Verified**,
2026-08-16, by directory listing. Neither directory contains a
`README.md` or a `cover-letter.md`, unlike every other series in the
compendium. *Remediation:* add both once the manuscripts are
otherwise finalized; the absence is itself weak evidence that these
two series have not been through the same submission-preparation
pass as their siblings, consistent with 11's incomplete state
(Major Issue 2.1) and 10's otherwise-complete state (this review
found no other defect in 10; its missing README appears to be pure
housekeeping lag rather than a signal of incompleteness).

**3.3.** *Location:* `analysis/report/06-component-decomposition/`.
**Verified** word count by this review (`wc -w report.Rmd` = 18,238
words), which does not match the internal compendium's stated
figure of "roughly 26,000 words" for the same file. *Remediation:*
this discrepancy was not reconciled by this review (the compendium
entry may refer to an earlier draft or a different counting method,
such as including LaTeX markup or the `report-slim.Rmd` companion
in a combined count); before acting on the compendium's
recommendation to promote the slim variant to canonical, re-measure
both files by the same method and confirm which, if either, exceeds
a target journal's length budget.

**3.4.** *Location:* series 03, 06 (contaminated-biomarker arm), 08
(deferred design grid), 09 (dropout prototype). **Inferred** from
the internal compendium, not independently re-derived except where
noted. Each presents a pilot or partial-factorial result using
language that, in earlier drafts, did not always distinguish the
pilot from the pre-registered production run; the compendium reports
this has been substantially though not completely remediated (03's
abstract now states the 240-replicate pilot framing explicitly; 09's
abstract no longer reads "[Forthcoming]"). One item in this class is
now further along than the compendium's 2026-08-12 snapshot recorded:
series 06's four chronological referee reports
(`analysis/report/06-component-decomposition/referee-report-2026-06-1{2,3,4,5}.md`,
**inspected**, 2026-08-16) show the second-round referee's single
remaining blocker was that "Study B" (the balanced-placebo recovery
analysis the paper's own framework requires to demonstrate its
remedy) was "pending implementation." The current `report.Rmd`
contains a completed "Study B: recovery under a belief-decoupling
design" subsection with reported bias figures (**inspected**,
2026-08-16, section present at the location given in Section 1's
summary of series 06 above), so the *implementation* blocker the
2026-06-15 referee identified is resolved, not merely reframed.
However, a related but distinct defect has emerged in its place: the
Study B numbers as implemented are not wired to a loaded data object
(Major Issue 2.11), which is a recurrence of a different item the
same referee report flagged (R3, numeric traceability) rather than a
fully closed loop. *Remediation:* a final
compendium-wide pass checking that every abstract's Results paragraph
states its replicate count and, where applicable, that the reported
cells are a subset of a larger pre-registered grid, would close this
class of issue definitively rather than paper by paper.

**3.5.** *Location:* `vignettes/01-simulate-trial-data.Rmd` line 33
and `vignettes/02-visualize-power-results.Rmd` line 32. **Inspected**
directly. Both vignettes describe themselves as walking through "the
steps used to generate the results published in `[___]`," an
unfilled placeholder citation. *Remediation:* fill in the citation
once a target venue is chosen (see Section 5), or rephrase to avoid
promising a specific publication that does not yet exist.

**3.6.** *Location:* `vignettes/02-visualize-power-results.Rmd`,
scattered plot titles and axis labels (e.g., "comapered" at line
236). **Inspected.** Minor spelling defects that would not survive a
CRAN check's spell-check-adjacent scrutiny of user-facing vignette
text, though they do not affect functionality. *Remediation:*
proofread pass; also apply the project's own US-English standard
(the vignettes otherwise read as American English).

**3.7.** *Location:*
`analysis/report/09-informative-dropout-by-design/report.Rmd`.
**Inferred** from the internal compendium's most recent, directly
cited verification pass. The design-by-dropout prototype remains
single-architecture (mean-moderation only) at a single biomarker-
effect value ($c_{bm} = 0.45$), which the manuscript itself
discloses as a scoped future extension rather than a completed
robustness check. *Remediation:* as the compendium already
recommends, extend to Architecture B and at least one additional
$c_{bm}$ value before treating the design ordering as
architecture-general.

**3.8.** *Location:*
`analysis/report/07-gompertz-evaluation/report.Rmd`. **Inferred**
from the internal compendium. The manuscript's own README promises
to identify non-saturating growth, biphasic, and breakpoint response
regimes; these are discussed narratively but not simulated. *Why it
matters is minor rather than major* because the paper's delivered
result (exact inertness under covariance-moderation) stands on its
own without the promised extension; the gap is between the README's
promise and the paper's delivered scope, not a defect in the
delivered result itself. *Remediation:* either simulate the
promised regimes or revise the README to match delivered scope.

**3.9.** *Location:*
`analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd` line 352
and `analysis/report/02-carryover-sensitivity/report.Rmd` line 406.
**Verified**, 2026-08-16, by running
`perl tools/notation-lint.pl analysis/report/*/report.Rmd`, the
project's own mechanized notation checker. Series 01 carries one
`bare-Dbc` violation (the code identifier `Dbc` used as prose outside
a code span, contrary to `NOTATION.md` rule 3) and series 02 one
`bare-Sx` violation of the same kind; every other series' `report.Rmd`
passes the linter cleanly. *Remediation:* wrap the two flagged
identifiers in `\texttt{}` or backticks; trivial, but worth fixing
before final render since the linter is already wired up and catching
it costs nothing. **Update, 2026-08-20:** series 02's flagged line has
shifted to line 573 following the manuscript's expansion (still
inside a code span, `nlme::lme(Sx ~ ...)`, so the linter appears to be
flagging a false positive rather than genuine bare prose usage; this
was not true of the finding at the original 2026-08-16 line number,
which was not re-inspected then). Re-run the linter after any further
edit to this file and confirm whether this is a persistent
false-positive worth an exception rule in `notation-lint.pl` or a
genuine remaining violation elsewhere.

**3.10.** *Location:* every series' `references.bib`
(e.g., `analysis/report/01-dgp-mean-moderation-vs-mvn/references.bib`
lines 274-311, `analysis/report/11-combined-dgp-architecture/
references.bib` lines 294-311). **Verified**, 2026-08-16, by grepping
`references.bib` in each series directory for `pmsimstats-paper`
citekeys. Every manuscript that cites a sibling companion paper
(01 cites 03, 08, 10, 11; 11 cites 01, 03, 10; and so on) does so
through a local `@unpublished{pmsimstats-paperNN, ...}` stub entry,
so the citations resolve typographically but point to manuscripts
that are themselves incomplete (series 11) or exist only inside this
same, unpublished repository. *Why it matters:* an external referee
evaluating any single series in isolation cannot verify a claim
deferred to "a companion paper," and where the companion is series 11
specifically, the deferred claim is not yet established even inside
the repository (Major Issue 2.1). This is a structural consequence of
submitting a tightly cross-referenced programme as separate papers
rather than a defect in any one manuscript. *Remediation:* before any
individual series is submitted, audit its `@unpublished{pmsimstats-
paperNN}` citations and either (a) confirm the target has itself been
submitted or posted (e.g., as a preprint) by the time of submission,
so the citation is verifiable, or (b) restate the deferred claim
in-line as an explicit limitation rather than a forward citation.
This is a specific enforcement mechanism for the partitioning decision
in Section 5.

**3.11.** *Location:*
`analysis/report/02-carryover-sensitivity/README.md`. **Inspected**,
2026-08-20. Dated 2026-07-30, the README describes the manuscript's
scope as "Analysis specifications E2 (exposure-weighted) versus E3
(lagged-treatment)" and "Sensitivity blocks S1 through S4, plus the
S6 cluster-robust recalibration," and states the stored specification
codes are `A1`/`A2`/`A3`. The current manuscript reports all nine
`G1`-`G9` specifications, spans sensitivity blocks through S12
(including the new Architecture A exploratory appendix), and the
stored codes, per `spec-labels.R`'s own header comment, migrated from
`A1`/`A2`/`A3` to `E1`/`E2`/`E3` (plus `E7`/`E9` and their CR2
variants) before this README was last dated. The directory-layout
listing also omits `supplement.Rmd`, `bullets.Rmd`, and
`report-devresults.Rmd`'s current role. This is the same pattern
Minor Issue 3.1 identified for series 01's README, now confirmed in a
second series. *Remediation:* refresh the README in the same pass
that addresses Major Issue 2.14, covering both submission documents
(README and cover letter) together rather than in separate passes,
since both drifted from the same underlying manuscript revision.

**3.12.** *Location:*
`analysis/report/02-carryover-sensitivity/bullets.Rmd`. **Inspected**,
2026-08-20. This file, newly added, is an explicitly-labeled,
non-manuscript working document (its own header comment states "Not
intended to be rendered as a manuscript") that restructures
`report.Rmd`'s prose paragraphs into bullet-point outlines for the
author's own use in drafting new narrative text. It carries the
manuscript's full YAML header, including the identical title, author,
and `knit:` rendering wrapper as `report.Rmd`, so it renders (when
rendered, as this reviewer did to confirm the build pipeline still
functions after the manuscript's edits) into a same-titled PDF staged
under `analysis/report/share/` alongside the real manuscript's
provenance-stamped copies, distinguishable only by filename
(`bullets-*.pdf` versus `report-*.pdf`). This is working material, not
a submission defect, but a stray same-titled PDF in the shared,
git-versioned `share/` directory risks confusion for a collaborator
who has not read this file's own internal disclaimer. *Remediation:*
none required before submission, since the file is self-labeled and
outside the manuscript's build chain; consider a distinguishing
subtitle or a `DRAFT OUTLINE` marker in the YAML title if this file
is expected to persist alongside the finished manuscript rather than
being deleted once its purpose (drafting new narrative prose) is
served.

## 4. What remains to be done

**(a) Required for correctness.**

1. Resolve series 11's boundary-validation gate (Major Issue 2.1):
   re-run the combined-architecture driver at the per-path
   sample-size convention, confirm numerical (not merely algebraic)
   agreement with series 01's published boundary power values within
   Monte Carlo tolerance, and re-stamp Table 1 and the abstract
   accordingly. Until this is done, series 11's Results and
   Conclusions should not be cited as established findings by any
   other manuscript in the compendium (note that series 11's own
   Section 6 already cites series 03 and is cited by no downstream
   paper yet, so the exposure is currently limited, but should not be
   allowed to grow before the gate is passed).
2. Remove the `browser()` call and fix the `rerun_simulations`
   default in `vignettes/01-simulate-trial-data.Rmd` (Major Issue
   2.2), so the vignette can be built non-interactively at all.
3. Fix the `simrsults`/`simresults` typo in
   `vignettes/02-visualize-power-results.Rmd` (Major Issue 2.3) so
   the rendered figure actually reflects the response-parameter-
   maxima sweep it claims to show.
4. Declare `data.table`, `corpcor`, `MASS`, `merTools`, and `ggplot2`
   in `DESCRIPTION`'s `Suggests:` field (Major Issue 2.4).
5. Reconcile or retire series 04's `report_short.Rmd` alternate: its
   RCT-power figure is not reproducible as documented and
   contradicts the primary `report.Rmd` pipeline (Major Issue 2.9).
6. Re-run all twelve of series 05's sensitivity sweeps at the
   manuscript's stated $N=35$ baseline; the current Results section
   is generated at a discarded $N=40$ baseline by the authors' own
   admission (Major Issue 2.10).
7. Wire series 06's Study B recovery numbers to a loaded data object
   (add the missing `readRDS()` call), repeating the fix a prior
   internal referee already prescribed for this paper (Major Issue
   2.11).
8. Rewrite series 07's `cover-letter.md`, which still asserts a
   pre-fix null finding that directly contradicts the current
   manuscript's reported result (Major Issue 2.12).
9. Replace the literal, unedited placeholder in series 09's
   Simulation Design section, and run the missing Type I error cells
   for OL+BDC, crossover, and hybrid designs before any design-
   ordering claim is treated as established (Major Issue 2.13).
9b. Rewrite series 02's `cover-letter.md`, which asserts a decay-form
    robustness claim the current manuscript explicitly corrects to a
    qualified finding (Major Issue 2.14); update the specification
    labels, replicate count, and cell count in the same pass.

**(b) Required for acceptance at the target venues recommended in
Section 5.**

10. Write the `rgt` prose for the eight series that still carry
    placeholder blocks, in the order: 10, then descending by
    placeholder count (06, 04, 08, 05, 07, 09), with 11 deferred
    until item 1 above is resolved (Major Issue 2.8); do not spend
    this effort on 04, 05, or 09 until items 5, 6, and 9 above are
    resolved, since prose written around currently-wrong or
    currently-stale numbers will need to be redone.
11. Fix series 04's comparison-fairness problem (Major Issue 2.6):
    add a measurement- or burden-matched sensitivity arm, or state
    and defend the asymmetry explicitly, and add observation count
    and participant burden to the Limitations section.
12. Either write `whitepaper-notation-audit.md`'s promised content or
    remove the claim that it exists from `NOTATION.md` and
    `analysis/report/README.md` (Major Issue 2.5).
13. Decide, and act on, the partitioning recommendation in Section 5
    below: fold or promote series 11 relative to 01 and 03 (Major
    Issue 2.7), settle series 05's relationship to 04, and settle
    series 07's relationship to 06.
14. Run the canonical 1,000-replicate production study for series 03
    (or continue to disclose the 240-replicate pilot as such, which
    is already largely done per the internal compendium) before
    treating its non-competitiveness finding as general.
15. Add a second application (a data set or DGP other than the
    prazosin/PTSD calibration) to series 10, since its claim to
    generality beyond N-of-1 trials currently rests on argument
    rather than a second worked example.
16. Add README.md and cover-letter.md to series 10 and 11 (Minor
    Issue 3.2) once the above is settled.
17. Route series 09 through an internal referee cycle comparable to
    what series 01, 06, 07, and 08 have already had, once item 9
    above is addressed.

**(c) Desirable polish.**

18. Reconcile the word-count discrepancy for series 06 between this
    review's direct measurement (18,238 words) and the internal
    compendium's figure (approximately 26,000 words) before deciding
    whether to promote its slim variant to canonical (Minor Issue
    3.3).
19. Fill the `[___]` placeholder citation in vignettes 01 and 02
    once a target venue is chosen (Minor Issue 3.5).
20. Proofreading pass on vignette 02's plot text (Minor Issue 3.6).
21. Refresh series 01's `README.md`, which currently undersells a
    manuscript that has already reached zero-placeholder state
    (Minor Issue 3.1).
22. Extend series 09's dropout robustness check to Architecture B and
    a wider $c_{bm}$ range (Minor Issue 3.7), after the correctness
    items above (item 9) are resolved.
23. Either simulate or rescope the trajectory regimes series 07's
    README promises but the manuscript does not deliver (Minor Issue
    3.8).
24. Fix the two notation-linter violations (`bare-Dbc` in series 01,
    `bare-Sx` in series 02) reported by `tools/notation-lint.pl`
    (Minor Issue 3.9).
25. Audit `@unpublished{pmsimstats-paperNN}` cross-citations in every
    series' `references.bib` against the target submission's actual
    availability before that series is submitted (Minor Issue 3.10).
26. Refresh series 02's `README.md`, which still describes a
    three-specification, S1-S4/S6-only scope and the retired
    `A1`/`A2`/`A3` stored-code labels superseded before the README's
    own stated date (Minor Issue 3.11).
27. Re-run `tools/notation-lint.pl` after any further edit to series
    02's `report.Rmd` and confirm whether its `bare-Sx` flag (now at
    line 573, inside a code span) is a persistent linter false
    positive or a genuine remaining violation (Minor Issue 3.9,
    updated).

## 5. Recommended framing

This section addresses how the eleven-manuscript compendium should
be partitioned into a smaller number of submittable papers, and how
each resulting paper should be framed against the existing
literature and current applied practice. The starting point is the
compendium's own single-estimand structure (documented in the
internal compendium and confirmed by this review's reading of
series 01, 10, and 11): every manuscript targets the biomarker-by-
treatment interaction in aggregated N-of-1 trials, calibrated
throughout to one prazosin/PTSD reference parameter set, reported
under the ADEMP simulation-study framework of Morris, White, and
Crowther (2019).

**What the surrounding literature and current practice already
cover, as context for what follows.** The N-of-1/aggregated-N-of-1
methods literature (Hendrickson et al. 2020 and its lineage) has
established the basic feasibility and design comparison of
aggregated N-of-1 trials for biomarker validation, but has not, to
this reviewer's knowledge from the manuscripts' own citation
apparatus, examined the sensitivity of power and bias conclusions to
*how* the biomarker-treatment interaction is encoded in a simulation
study's data-generating process; this appears to be a genuine,
under-served methodological gap and is series 01's real
contribution. Separately, the sandwich/cluster-robust standard-error
literature (Liang and Zeger 1986 onward, including the small-sample
CR2 and Kenward-Roger corrections) is long-established in general
mixed-model and GEE practice, but its application specifically to a
Wald test of a fixed-effect interaction, staged as a four-channel
diagnostic protocol with a model ladder and a sample-size gradient,
is a genuinely reusable methodological instrument that this reviewer
did not find pre-packaged in this form in the citations series 10
itself supplies; this is a real, if modest, methodological
contribution independent of the N-of-1 application. By contrast, the
finding that repeated within-subject measurement increases power
relative to a single-measurement parallel design (series 04, 05) is
standard and would not by itself constitute a defensible contribution
at a statistical methods venue; its value is clinical and applied,
contingent on the comparison-fairness problem (Major Issue 2.6) being
resolved.

**Paper group 1: the DGP-architecture papers (01, 03, 11).**
*Plausible framings:* (i) a single methodological paper on
DGP-architecture sensitivity in simulation-based power studies,
generalizing beyond N-of-1 trials; (ii) three separate papers, as
currently structured; (iii) a two-paper split, series 01 standalone
(the two-architecture comparison) with 03 and 11 combined as a
"spectrum and combined-channel" companion. *Recommendation:*
adopt (iii), consistent with the internal compendium's own
(unadopted) recommendation to fold 11 into 01, but this review goes
one step further: 11's material is thematically closer to 03 (both
concern the *general model space* of which mean-moderation and
covariance-moderation are boundary cases; 03's latent-class spectrum
and 11's dual-channel model are, by 11's own Discussion, the same
generalization viewed from two angles) than to 01 (which is a
head-to-head empirical comparison of the two pure architectures
under carryover, a self-contained result that does not need the
general model space to make its point). Combining 03 and 11 also
solves 11's standalone-length and validation problems (Major Issues
2.1 and 2.7) directly: 11's super-additivity derivation becomes a
short analytical section inside a paper whose main empirical weight
is carried by 03's latent-class results, and the boundary-validation
requirement becomes an internal consistency check reported in a
Methods subsection rather than the entire empirical payload of a
freestanding paper. *Implications:* series 01 should be titled and
framed exactly as it already is (the architecture-choice paper,
targeting *Statistics in Medicine* or a comparable biostatistics-
methods venue, since its contribution is a warning about a modeling
choice under-examined in the surrounding literature); the 03+11
combination should be retitled around the general model space (for
example, "A general model for biomarker-treatment interaction
encoding: latent-class, continuous, and combined-channel
formulations"), with 03's identifiability theorem as the paper's
mathematical spine, its 240-replicate pilot promoted to a full
production run (Checklist item 9) or explicitly retained as a scoped
pilot with 11's analytical super-additivity result folded in as a
second, complementary theoretical result once its own boundary check
(Checklist item 1) is resolved. Series 06 and 07's simulation
apparatus (generative-substrate robustness) should be de-emphasized
in this combined paper and referenced only as supporting citations,
not re-derived.

**Paper group 2: analyst-side inference under a fixed DGP (02, 08,
10).** *Plausible framings:* (i) three independent methods papers
(current structure); (ii) a single "how to analyze aggregated N-of-1
interaction data" paper covering carryover-predictor specification,
test-procedure choice, and covariance-calibration together; (iii)
02 and 08 combined as a paired "analyst decision" paper (predictor
specification and test procedure are both decisions the analyst
makes about the same fitted model), with 10 kept fully separate.
*Recommendation:* adopt (iii). 10 is the compendium's only
contribution that is not domain-specific (Section 1's summary and
this review's direct reading both confirm this), and folding it into
an N-of-1-specific paper would bury a general mixed-model result
inside a specialized application, costing it its natural, larger
audience. 02 and 08, by contrast, are both squarely about analyst
choices *within* the N-of-1 carryover context and share
infrastructure, reference cells, and (per the internal compendium)
partially overlapping specification labels (E1/E2/E3). *Implications:*
target 10 at a biostatistics-methods venue with a broader audience
than *Statistics in Medicine*'s clinical-trial-methods readership,
title and abstract framed entirely around the general Wald-test-
calibration problem with the N-of-1 example demoted to "a worked
example" in the title (as the current title already does correctly:
"Type I calibration of fixed-effect interaction tests in linear
mixed models: a staged diagnosis..."); Checklist item 15 (a second
application) becomes considerably more important under this framing,
since a methods paper aimed at a general mixed-model audience will
be judged partly on whether its diagnostic protocol generalizes
beyond the one worked example, which is currently the paper's most
exposed weakness (Section 1's summary, "Two gaps" for series 10).
For the combined 02+08 paper, target *Statistics in Medicine* or a
comparable clinical-trials-methods venue, with the title organized
around "analyst-side defaults are conditional, not universal" (the
internal compendium's own phrase for this cluster's collective
finding), the cycle-by-period design grid promised by 08's current
title either completed or dropped from the title (Checklist item 13,
"decide paper 08's scope," which this recommendation resolves by
subsuming it into the combined paper's narrower, already-delivered
scope).

**Paper group 3: the generative substrate (06, 07).** *Plausible
framings:* (i) two independent papers (current structure); (ii) one
paper with 07 as a robustness section, as the internal compendium
already recommends; (iii) fold both into paper group 1 as
supporting material, since their finding (trajectory shape and
component contamination affect bias/power in specific, bounded ways)
functions as a robustness argument for the DGP framework used
throughout the compendium rather than as a freestanding scientific
claim. *Recommendation:* adopt (ii), as the internal compendium
recommends, with one addition: the omitted-variable-bias identity in
06 (Study A/B) is a real, general econometric-style result
(applicable to any simulation study with a contaminated covariate)
and deserves a title and abstract that lead with that generality
rather than with the balanced-placebo-design application, which
(per the internal compendium's own assessment, concurred with here)
is "rarely feasible in the clinical settings the programme targets."
*Implications:* retitle around the identity itself (something like
"An omitted-variable-bias account of biomarker contamination in
simulated treatment-response trajectories, with a design-based
remedy"), slim the manuscript toward `report-slim.Rmd`'s length
(pending the word-count reconciliation of Minor Issue 3.3), append
07's exact-inertness result as a short robustness section rather
than a companion paper, and move the promised-but-undelivered
trajectory-regime sweep (Minor Issue 3.8) to explicitly scoped future
work rather than implying it is covered.

**Paper group 4: the clinical/design comparison (04, 05, 09).** This
group carries the compendium's most severe correctness problems
(Major Issues 2.9, 2.10, and 2.13), which materially change the
framing recommendation from what the internal compendium's own
self-assessment suggested. *Plausible framings:* (i) three
independent papers (current structure, with 05 already recommended
by the internal compendium as supplementary to 04); (ii) 04
standalone with 05 as an appendix, and 09 kept fully separate since
it addresses dropout, a different axis, under a different estimand
(design robustness rather than the main-effect comparison); (iii)
all three combined as a single "design comparison for aggregated
N-of-1 trials" paper. *Recommendation:* adopt (ii) as the eventual
target structure, but not yet: unlike every other paper group in
this section, group 4 is not currently a framing decision so much
as a correctness-remediation queue. 04's two pipelines disagree on
the headline number (Major Issue 2.9); 05's entire Results section
is self-disclosed as computed at the wrong baseline (Major Issue
2.10); and 09 has not been through a referee cycle and is missing
most of its pre-registered Type I error cells (Major Issue 2.13).
None of the three can be responsibly framed for a target venue,
comparator set, or abstract until these are fixed, because the
framing recommendation itself (05 as 04's appendix; 09 standalone)
depends on results that do not yet exist in a verified form. This
review's assessment of 09 as "close to being a complete, defensible
standalone finding" should be read as *conditional* on Major Issue
2.13's remediation (the missing Type I cells and the untested
randomization-path mechanism), not as a statement about its current
state. *Implications:* treat Checklist items 5, 6, and 9 (the
correctness fixes for 04, 05, and 09) as strict prerequisites to any
framing or venue decision for this group, ahead of and separate from
Major Issue 2.6's comparison-fairness fix, which was this review's
original, and still valid, top concern for 04 specifically. Once all
of these are resolved: target a clinical or general-medical-
statistics venue for 04 (with 05 as supplementary material,
consistent with 05's own one-factor-at-a-time structure being better
suited to a reference appendix than a standalone results section),
and target the same class of venue as 02+08 (clinical-trials
methods) for 09, emphasizing the dropout-coupling finding as the
compendium's clearest genuinely novel, underexplored result relative
to the cited N-of-1 literature (this review did not identify a
citation in 09 to prior work examining design-family-conditional
dropout in this setting specifically, though this was not
independently verified against the full N-of-1 literature beyond
the manuscript's own citations) — but only once that finding rests
on a complete Type I error grid and a tested, not merely asserted,
explanatory mechanism.

**Net recommendation.** Partition the eleven series into five
submission units: (1) series 01 alone; (2) series 03 and 11 combined;
(3) series 02 and 08 combined; (4) series 10 alone; (5) series 06 and
07 combined; plus two further units contingent on resolving,
respectively, Major Issues 2.6 and 2.9 (04) and Major Issue 2.10
(05): (6) series 04 with 05 as supplementary material; and, pending
Major Issue 2.13: (7) series 09 alone. This is fewer units than the
internal compendium's own "seven or eight" estimate because this
review additionally folds 11 into the 03 pairing rather than into 01,
for the reasons given above. Units 6 and 7 additionally carry a
higher remediation bar than the internal compendium's estimate
credits, because this review's direct reading surfaced correctness
defects in 04, 05, and 09 that the compendium's self-assessment does
not report. A referee encountering any two of these seven units side by
side should not perceive redundant framing, in contrast to the
current risk (Major Issue 2.7) that a referee encountering 01, 03,
and 11 together would notice the same architecture-spectrum argument
made three times.

## 6. Assessment

This reviewer's overall verdict, applied to the compendium as a
prospective set of submissions rather than to any single manuscript,
is **major revision, not yet ready to submit as structured**, with
substantial credit for the underlying work.

The programme has a real, coherent, and in places genuinely novel
scientific core: the demonstration that DGP-architecture choice is a
consequential and previously under-examined modeling decision (01),
the staged Wald-test calibration protocol and its counter-intuitive
finding of test conservatism from covariance misspecification (10),
the omitted-variable-bias account of biomarker contamination (06),
and the design-by-dropout coupling result (09, once completed) would
each, on the evidence reviewed, be defensible contributions at a
competent biostatistics-methods venue once the specific defects
above are closed. Series 10 in particular is close to submission-
ready: this review, having read it in full, found no major issue
with its mathematics, its reproducibility chain (every table traces
to a present script and a present `.rds` file), or its argument,
only the housekeeping gaps of a missing README/cover-letter and the
Checklist item of a second worked example. Series 08 is the clearest
example of successful revision in the compendium, having resolved
six of seven issues raised by its own 2026-06-17 internal referee
report.

Set against that, four classes of problem are serious enough to
warrant the major-revision verdict rather than a lesser one. First,
and now the most consequential class after this review's full-text
reconciliation pass, three manuscripts carry outright correctness
defects, not merely framing or packaging gaps: series 04 reports two
mutually contradictory headline numbers from two non-reconciled
pipelines, one of them undocumented and irreproducible (Major Issue
2.9); series 05's entire Results section is generated from a
baseline its own authors' note identifies as discarded and stale
(Major Issue 2.10); and series 09 retains an unedited content
placeholder in running prose and is missing most of its
pre-registered Type I error cells (Major Issue 2.13). These three
are not polish items; no numeric claim in any of the three
manuscripts' current Results sections should be treated as
established until they are fixed, and this class of problem was not
surfaced by the repository's own internal self-assessment, which
this review otherwise found to be reliable and well-corroborated for
series 01, 02, 03, 06, 07, 08, 10, and 11. Second, series 11's
central empirical claim is, by the manuscript's own explicit
statement, not yet validated against its own correctness gate (Major
Issue 2.1); no paper should be submitted, nor should other papers in
the compendium build on it, before that gate is resolved. Third, the
compendium as currently partitioned into eleven units invites, and
in the case of 01/03/11 already exhibits, the redundant-framing
problem a referee is trained to catch (Major Issue 2.7, Section 5);
the recommended five-to-seven-unit partition in Section 5 is not
optional polish but a precondition for the compendium being
creditable as eleven distinct scientific claims rather than one
finding restated with variations. Fourth, the package-level materials
(vignettes, `DESCRIPTION`) that would accompany any software
component of a submission are not currently in a state that would
pass `R CMD check --as-cran` (Major Issues 2.2-2.4), which matters
because at least one recommended submission unit (the series 10
methods paper) would plausibly be strengthened by, or a software
paper could plausibly be built around, a demonstrably working
package, and the vignettes as committed do not demonstrate that.

The first of these four classes changes the character of the overall
verdict from what an earlier pass of this review, relying more
heavily on the internal compendium's self-assessment, would have
concluded. Series 01, 02, 03, 06 (Study A), 07 (after its own
referee-driven fix), 08, 10, and 11's underlying statistical content,
where this review inspected it directly, is sound; the defects there
are, respectively, an unresolved verification step (11), a
packaging and framing decision (the partition), a software-hygiene
gap (vignettes), a numeric-traceability gap (06's Study B, Major
Issue 2.11), and a stale-document contradiction (07's cover letter,
Major Issue 2.12) rather than errors in the mathematics or the
simulation logic itself. Series 04, 05, and 09, by contrast, cannot
currently be credited with a verified Results section at all: 04's
primary number is contradicted by its own alternate, 05's numbers
are conceded stale by its own authors, and 09's evidentiary base is
majority-unexecuted against its own pre-registration. That
combination, three of eleven manuscripts with unresolved correctness
problems in their headline results, is unambiguously consistent with
major, not minor, revision, and is more serious than "requires
authorial judgment on partitioning": it requires re-running
simulations and reconciling pipelines before those three manuscripts
can be assessed as scientific claims at all.

This update (2026-08-20) does not change the overall verdict, but
adds one data point to the third and fourth problem classes above.
Series 02's manuscript body itself is, if anything, now stronger than
at the initial review: the nine-specification extension, the
Architecture A exploratory appendix, and the uniform 500-replicate
precision pass are all sound, reproducible improvements this update
verified directly by re-rendering both documents and spot-checking
reported numbers against the underlying data. But the discovery that
its cover letter was not updated in step, and now contradicts the
manuscript's own corrected headline qualifier (Major Issue 2.14),
means this failure mode (submission-facing documents drifting out of
sync with a revised manuscript body) is now confirmed in two of eleven
series (07 and 02) rather than one, and in a manuscript this review
had not previously flagged as carrying any Major Issue. That
strengthens, rather than weakens, the case in Checklist item 9b and
the compendium-wide process recommendation in Major Issue 2.14 for
treating cover-letter and README review as a mandatory step tied to
manuscript revision, not a one-off cleanup task performed only when a
defect happens to be noticed.

## 7. Revision history

- 2026-08-16: Initial review. Covered all eleven report series
  (10 and 11 read in full independently; 01-09 corroborated against
  the repository's own internal compendium,
  `whitepaper-compendium-summary.md`, by direct grep-based
  verification of placeholder counts, word counts, and
  README/cover-letter presence, all of which matched to within one
  line), the three package vignettes (read in full; two code defects
  and one dependency-declaration defect identified that the internal
  compendium's scope did not cover), and the cross-cutting notation
  and audit infrastructure (one defect identified: the notation
  audit document is an empty stub despite being cited as
  load-bearing). Issued a five-to-seven-submission-unit partitioning
  recommendation, refining the internal compendium's own
  seven-to-eight-unit estimate by regrouping series 11 with 03
  rather than with 01. Overall verdict: major revision.
- 2026-08-16 (same-day corroboration pass): a second, independent
  read of this workspace, run later the same day, read series 01 and
  03 in full (rather than relying on the internal compendium for
  those two, as the initial pass had done), re-verified all eleven
  series' abstracts directly, ran the project's own
  `tools/notation-lint.pl` against every `report.Rmd`, checked
  `analysis/report/NOTATION.md`'s sample-size-convention audit trail,
  and read all four of series 06's chronological referee reports.
  No finding of the initial pass was contradicted. Two new minor
  issues were identified and added: a notation-linter violation each
  in series 01 and 02 (Minor Issue 3.9) and a structural reliance,
  compendium-wide, on `@unpublished` self-citations to sibling
  in-progress manuscripts, which creates verification circularity for
  an external referee and is a concrete enforcement mechanism for the
  Section 5 partitioning recommendation (Minor Issue 3.10). One
  previously identified item was confirmed resolved rather than
  merely reframed: series 06's referee-flagged "Study B ... pending
  implementation" blocker (from the 2026-06-15 second-round referee
  report) is now a completed subsection in the current `report.Rmd`
  (Minor Issue 3.4, updated). The major issues, checklist, framing
  recommendation, and overall major-revision verdict of the initial
  pass are unchanged and are corroborated by this pass's independent
  reading of series 01 and 03.
- 2026-08-16 (reconciliation pass): a full-text, direct read of
  series 04, 05, 06, 07, 08, and 09 (previously covered only via
  corroboration against the internal compendium) surfaced five new
  major issues that the internal compendium's self-assessment does
  not report and that materially change the verdict for three
  manuscripts: series 04 reports two contradictory headline results
  from non-reconciled pipelines, one undocumented and irreproducible
  (Major Issue 2.9); series 05's Results section is self-disclosed by
  its own authors as computed at a discarded, non-canonical baseline
  (Major Issue 2.10); series 06's newly added Study B recovery
  numbers are not wired to any loaded data object, a recurrence of a
  defect a prior internal referee had already flagged for this paper
  (Major Issue 2.11); series 07's cover letter contradicts the
  current manuscript's corrected headline finding (Major Issue 2.12);
  and series 09, which has not been through any internal referee
  cycle, contains a literal unedited content placeholder in running
  prose and is missing most of its pre-registered Type I error cells
  (Major Issue 2.13). The Section 4 checklist, Section 5 framing
  recommendation for paper group 4 (04/05/09), and Section 6
  assessment were revised accordingly; the overall verdict remains
  major revision but on more serious grounds than the initial pass
  identified, since three of eleven manuscripts (04, 05, 09) now
  carry unresolved correctness defects in their headline results
  rather than only framing, packaging, or software-hygiene gaps.
  Series 08 was independently confirmed, by the same full-text pass,
  to have resolved six of seven issues from its own 2026-06-17
  referee report, the compendium's clearest example of successful
  revision. No previously reported finding was contradicted by this
  pass; series 01-03 and 10-11's assessments are unchanged.
- 2026-08-20: Update pass. Verified via `git log`/`git diff --stat`
  that no commits since 2026-08-16 touched any series other than 02,
  any vignette, `DESCRIPTION`, or the cross-cutting notation/audit
  documents; all findings for series 01, 03-11, the vignettes, and
  the cross-cutting infrastructure are therefore carried forward
  unchanged without re-reading, on the strength of that verification.
  Series 02 was substantially revised in the interim (nine-
  specification G1-G9 extension of Figures 1-2, a report/supplement
  split, a new exploratory Architecture A appendix, and an
  all-simulations-to-500-replicates precision pass) and was re-read
  in full, including the two new files (`supplement.Rmd`,
  `bullets.Rmd`); both documents were re-rendered successfully and a
  sample of newly reported power values was checked directly against
  the underlying `.rds` data. The manuscript body's own content is
  sound and, where checked, numerically consistent with its data; no
  correctness defect was found in `report.Rmd` or `supplement.Rmd`.
  One new Major Issue was identified: series 02's `cover-letter.md`
  was not revised alongside the manuscript and now asserts a
  decay-form robustness claim the manuscript's own Abstract explicitly
  qualifies as failing at the two most extreme decay shapes tested
  (Major Issue 2.14), the same failure pattern previously found only
  in series 07 (Major Issue 2.12), now confirmed in a second series.
  Two new Minor Issues were identified: series 02's `README.md` is
  stale relative to the current manuscript scope and stored-code
  labels, mirroring series 01's previously flagged README staleness
  (Minor Issue 3.11); and the new `bullets.Rmd` working document
  shares its YAML title with the manuscript and renders to a
  same-titled PDF under `share/`, a housekeeping note rather than a
  defect (Minor Issue 3.12). Minor Issue 3.9's `bare-Sx` finding for
  series 02 was re-located (now line 573) and flagged as a possible
  linter false positive, not newly re-verified as a genuine violation.
  The Section 4 checklist and Section 6 assessment were updated
  accordingly; the overall verdict remains major revision, unchanged
  in kind but with one additional data point (a second cover-letter/
  manuscript contradiction) supporting the process-fix recommendation
  in Major Issue 2.14.
- 2026-08-20 (same-day remediation note): Major Issue 2.14 (Checklist
  item 9b) was resolved within the same session: `cover-letter.md`
  was rewritten to state the qualified decay-shape finding (advantage
  narrows to statistical insignificance under the two most
  heavy-tailed Weibull shapes examined, matching the current Abstract
  verbatim), updated from the retired `A1`/`A2`/`A3` labels to the
  current `G1`-`G9` display codes and specification descriptions,
  corrected the cell count and replicate count to match the current
  216-cell Tier 1 grid at 500 replicates per cell (removing the stale
  "540 cells," "24-cell expansion," and "600 replicates" claims), and
  re-dated to 2026-08-20 (**verified**, by direct edit and re-read of
  the file). Minor Issue 3.11 (README staleness) and Minor Issue 3.12
  (bullets.Rmd housekeeping note) remain open, as does Checklist item
  27 (re-run the notation linter). Major Issue 2.14 should be marked
  resolved in any subsequent full review pass.
