---
geometry: left=0.85in,top=0.8in,bottom=0.8in,right=0.85in
fontsize: 10pt
header-includes:
  - \linespread{1.15}
  - \setlength{\parskip}{0.3em}
---

# Submission readiness of the fifteen manuscripts

*2026-11-09 19:00 PST*

**Author.** pmsimstats team

**Purpose.** A paper-by-paper statement of how far each manuscript under
`analysis/report/` is from journal submission as of 2026-11-09, what
blocks it, and what it would cost to clear. This document supersedes
`docs/66-manuscript-readiness-assessment.md` (2026-10-29, fourteen
papers). Paper 15 is new since then.

## Method and evidence

The method is that of `docs/66`. For each paper the source is the newest
internal referee report in its directory. Each major finding was checked
against the commits made after that report and against the current
source. A finding is marked addressed only where the fix was seen in the
source or a commit diff. Five assessment agents worked in parallel, three
related papers each, and their grades were then checked against one
another and against the scale.

The evidence base differs by paper:

- **Current** (no rewrite since the latest referee report): 02 (report
  2026-11-08), 03, 12, 13 and 14.
- **Rebuilt or rerun since the report, not re-refereed:** 01 (reframed
  in 1619ceb9d after the 11-05 report), 04, 05, 06, 07, 08, 09, 10
  and 11 (late-October reports).
- **Never refereed:** 15 (first draft, d07c71374, 2026-11-07).

No full referee pass was run for this assessment (`--referee` was not
given). No simulation was rerun and no manuscript was rendered. Headline
numbers were recomputed from the stored result files where the agents
could do so with Rscript. Each claim below is verified (executed and
observed), inspected (read in source or diff) or inferred, and the
labels are given where they matter. Effort estimates are rough, and
grades are the assessor's judgment.

**Grades.** **A** submittable after copy-edit; **B** minor revision;
**C** major revision, no new simulation; **D** needs new or rerun
simulation; **E** incomplete draft (results section empty or pending).
A grade in parentheses applies if an optional or contested analysis is
required.

## Summary table

| Paper | Short title | Target | Grade | Previous | Effort to submit |
|---|---|---|---|---|---|
| 01 | Toward an optimal DGP | *Trials*, Methodology | B | B | 1-2 days + re-review |
| 02 | Predictive biomarker under carryover | *Trials* | B (C) | C (D) | 2-4 days + re-review |
| 03 | Latent-class / mixture | SiM (format: *Trials*) | E | E | 1.5-2 weeks |
| 04 | Main effect, linear-DGP pilot | undecided | B (D) | C (D) | 1-2 days + re-review |
| 05 | Hybrid power, main effect | SiM | D | D | 3-4 days |
| 06 | Three-component decomposition | SiM | C (D) | C (D) | 4-7 days |
| 07 | Response-curve family | SiM, Methods Note | B (D) | B | 1 day + re-review |
| 08 | Test-procedure choice | SiM | B | B | 1 day |
| 09 | Informative dropout | *Trials* (SiM: D) | B (D) | B | 1-2 days + re-review |
| 10 | Interaction-test calibration | none stated | C (D) | D | 1-2 days |
| 11 | Dual-channel DGP | SiM (assumed) | C | C | 2-3 days, or fold |
| 12 | Closed-form design efficiency | SiM, short comm. | D/E | D/E | 1-2 weeks, or retire |
| 13 | Two-component, drop PB | none named | D/E | D/E | 3-4 days with 14 |
| 14 | Two-component, drop TV | none named | D | D | (with 13) |
| 15 | Prognostic biomarker | none named | D | new | 1-1.5 weeks |

Three grades moved up: 02 (the bm:t check and the repaired-matrix rows
left with the paper 15 split), 04 (every finding addressed after the
pilot rebuild) and 10 (the four pending blocks are cleared). None moved
down. Papers 07 and 09 gained a parenthetical D for contested reruns.
Three papers (03, 12, 13) have new problems that do not change their
grade.

**Common to every paper.** The author roster is a placeholder
("Temperance Persons") or the team name. Declarations are incomplete or
missing. No paper has an archived, DOI-bearing code release. Papers 01
and 02 additionally need your decisions on authors, the competing-
interests statement (a shared authorship with Hendrickson et al. 2020)
and the build settings. None of these takes long, but they block every
submission. The journal word and abstract limits quoted below were not
checked against the live author pages.

## Paper-by-paper assessment

**01. Toward an optimal DGP. Grade B.** The 11-05 referee asked for
major revision on six points, and five are addressed in substance.
Planning is shown at serial correlations of 0.5 and 0.9 (ce827a4ef). G10
is run under the prognostic pathways (1619ceb9d). The stress test is
rerun under the DGP as defined (6ad5d4948). The literature is filled in
(d80ae7ab0). Every spot-checked number reproduces (verified), including
the bound of 0.474 at $\rho = 0.7$ and 0.340 at 0.9 and the Hybrid E9
sample sizes of 314, 278 and 154. What remains is mostly administrative:

- 6ad5d4948 removed the Funding, Authors' contributions and
  Acknowledgements headings, which *Trials* requires (inferred).
- Competing interests reads "none" and omits the shared authorship.
- The data-use-agreement confirmation is a placeholder.
- `title-page.md` (retired title) and `README.md` (SiM style) are
  stale; `appendix-b.Rmd` is still an orphan.
- The bib title for paper 02 is stale.
- The main text is about 14,000 words; the abstract is about 345 to
  365 words against a limit thought to be 350.

The reframing in 1619ceb9d added results no referee has read, so a
short re-review should precede submission. One discrepancy with paper
15 must be reconciled first (cross-cutting finding 2). The science is
settled. The risk is in presentation: "optimal" in the title sits
beside text saying the DGP is not yet optimal.

**02. Predictive biomarker under carryover. Grade B (C if the
published-set reruns must be reported).** The 11-08 referee asked for
major revision, and 8d491ccd4 answered nearly every point from stored
output. It reports the share of the carryover loss G8 recovers (15-44%
in Hybrid, 19-25% in OL+BDC; verified), adds a parameter table and the
fixed-half-life comparator, and adds G10 coverage (0.934-0.955) and a
G10-versus-G6 equivalence test. The random-slope fits (0d16536db,
d56339fb3) close m17 without moving any number. The bm:t question and
the repaired-matrix cells moved to paper 15, which resolves two
`docs/66` items. Still open:

- The Recommendation still uses the "within one paired SE in every
  family" criterion the referee rejected, and the family MCSEs ignore
  shared random-number streams.
- The hard-coded "(Section 2.8)", the build settings and the unused
  bib entries (m7, m20, m21).
- `cover-letter.md` and `title-page.md` address SiM under an old title
  and cite a retracted finding; `README.md` tells the reader to call
  `rmarkdown::render` directly.
- About 15,100 words of main text and a 370-word abstract.

The paper runs on the simcore2026 parameter set (SDs 6/2/8), while
papers 01 and 15 use the published set (10/10/8). This is disclosed.
Paper 15, however, says a paper-02-era conclusion does not hold under
the published values, which invites a referee to ask why 02 does not
report its stored published-set reruns. Paper 02 hands the case of a
biomarker of unknown role to paper 15, so its submission is tied to
that grade-D companion.

**03. Latent-class and mixture formulations. Grade E.** No commit has
touched the folder since 245e4369f (10-28), but the situation has
changed. Both blocks of the corrected study have now run under the
published set (null on 11-03, alternative on 11-04; verified). The
outputs sit in an untracked folder, and the Rmd reads a `summary.rds`
that does not exist, so the Results still print the "has not
completed" fallback and the `RESULTS-PENDING` marker remains (verified).
The new data do not support the conclusion the draft anticipated. Under
the null on the Hybrid design, the single-resample BLRT rejects 0.088,
0.136 and 0.200 at N = 35, 70 and 140, and the within-class Wald test
rejects about 0.25 (verified). The linear mixed model with CR2 holds
its size (0.046-0.052) and has the best power. Other problems:

- M7 has regressed: the Background cites paper 01 for the 15.4% and
  6.9% figures, which paper 01 no longer contains.
- Paper 01 Appendix A.6 now holds the paper's only theoretical result.
- The Methods do not state the parameter set, and the record cells use
  older sets.

A size-adjusted comparison can probably be made from the stored
statistics without new simulation (inferred). The paper should be
reframed as an empirical study of class-aware tests at N-of-1 sample
sizes. It moves to C once the Results are written from the stored runs.

**04. Treatment main effect, linear-DGP pilot. Grade B (D if the
Hybrid rebuild is required).** The 10-28 referee recommended rejection.
Commits 270a729d8 and 33fb58b1b rebuilt it as a pilot that presents
its closed form first, and all eleven major findings are now addressed
in the source. The 15 abstract guards pass against the stored run
(verified): the exclusion-by-schedule rule has power 0.84 and 0.83, and
null rates lie between 0.043 and 0.054. What remains is contribution and
venue. The simulation confirms a closed form, and the Hybrid, treatment
heterogeneity and placebo response are not simulated. The cover letter
(May) still describes the retracted study and addresses SiM; the title
page suggests *Trials* or *Clinical Trials*; the README says the tables
are pending. At about 3,500 words it reads as a short methods note. It
no longer shares an engine with paper 05, though the two share an
estimand.

**05. Hybrid power for the main effect. Grade D.** Nothing has changed
since `docs/66`. The `nof1power` engine still lacks a baseline-to-
follow-up shift term (verified), so the alternating design and the
crossover over-reject (0.094 and 0.148) and their power is not
comparable. The paper discloses this, but its abstract still reports
their power. The engine also still compounds carryover, and its
cross-component correlation is 0.1/0.05. The comparison of the Hybrid
with the parallel group (0.778 against 0.798) is sound as it stands.
Clearing it needs 1-2 days of engine work, about 3 hours on the host
and a day of rewriting. A further problem: the engine version 0.1.0
did not change with the censoring sign fix (24f187bc6), and the
artifacts record an empty `git_rev`, so a run cannot be tied to code.

**06. Three-component decomposition. Grade C (D if $c_1 = 0.2$ or a
Study 3 rerun is required).** d0070b189 and eeda569ef dealt with the
substance of the 10-28 report: seven of ten major findings are
addressed, and the identity predicts the sign of the displacement on
every design (0.780, 0.612 and 1.203 of the predicted size; verified).
The paper has since been overtaken by the parameter remediation. The
published-set reruns finished on 11-04 and 11-05 but are untracked and
not read by the Rmd (verified). The Rmd types the fixture values (true
slope $-2.25$, $\sigma_{PB} = 2$), and one render guard would fail on
the new data, so integration is more than repointing paths. The new
numbers change the size statements: balanced-design null size is 0.073
to 0.080 under every analysis (verified). The June Study 2 and Study 3
runs predate anchored carryover, which the text does not state. Length
is the largest cost: about 17,400 words of main text against roughly
8,000. `report-slim.Rmd` still carries the withdrawn headline.

**07. Response-curve family. Grade B (D if $c_1 = 0.2$ or a G10/CR2
arm is required).** The common-random-numbers rebuild (0c0f1e260)
addressed all nine major findings, and its stored numbers match the
referee's independent run exactly (verified). One new error was found:
the pooled null pools four near-identical fits per seed, so its MCSE
(0.0020) should be about 0.004 and its $z$ of $-6.7$ about $-3.3$
(verified). That changes one abstract number and one sentence. Other
issues:

- Companion work is cited by repository path in prose, with no bib
  entries.
- It describes paper 01 as the home of the architecture comparison.
- Its conjecture for the size deficit (working-covariance mismatch)
  contradicts paper 10's diagnosis.
- The test is the historical lme with model-based SEs, on OL+BDC, with
  cross-component correlations of 0.2/0.1.

Each optional rerun is about 15 minutes locally. It has not been
re-refereed since the rewrite.

**08. Test-procedure choice. Grade B.** 0bd472f4c addressed all seven
major findings, and 091953725 moved every number to the published-set
run. The setup chunk refuses to render if the parameter set or any
narrated direction changes, and that guard passes today (verified).
Remaining issues:

- Its mechanism sentence ("a working-covariance calibration defect")
  cites paper 10 for the opposite of what paper 10 now shows.
- It names G9 as the series' recommendation; papers 01 and 02 now use
  G10.
- The README (CR2 arm "pending", wrong path B week), title page and
  cover letter (five procedures, not six) trail the manuscript, and
  the author line differs between files.
- `report-slim.Rmd` and `report-mini.Rmd` are unmaintained variants.
- Novelty is narrow, and it overlaps paper 10.

At about 7,000 words it fits the SiM guidance.

**09. Informative dropout. Grade B for *Trials* (D if a referee
requires $c_1 = 0.2$, a G9/G10 arm or F1); D for SiM.** The rebuild
and the published-set rerun (8918e5c6e through 091953725) addressed
every correctness finding of the 10-28 report: the mechanism, labels,
the OL row, bias against a stated true value of $-3.6$, calibration,
pairing, stale numbers and the seed statement. Every number is read
from stored output, and the cover-letter figures match the files
(verified). The rerun reversed two of the rebuild's own conclusions.
For example, CO's larger loss is no longer explained by its lower
starting power. No referee has read either version. Other issues:

- One parameter point, cross-component correlations of 0.1/0.05, no
  shared-parameter dropout.
- An abstract of about 405 words, over either journal's limit.
- The title page says 500 null replicates; the run has 1,000.
- The README is stale, and `09-dropout-prototype-10min.R` has an
  uncommitted change.

**10. Interaction-test calibration. Grade C (D if a G10 arm is
required, which consistency with 01 and 02 will likely demand).** The
rewrite (6d801ed8f) and zeus job 55 (b611d4998) addressed all ten major
findings, and the four `RESULTS-PENDING` blocks counted in `docs/66` are
now inert fallbacks (verified by running the chunks). The paper is
behind again, in two ways. First, job 55 ran under simcore2026. The
published-set rerun finished on athena on 11-04 but is untracked and
not read. In it CR2 sizes fall from about 0.050 to 0.044 (verified),
so "restores the size" needs re-reading. Second, the paper calls G9
the recommendation and reports G9's selection increment (+0.010),
which is the reason papers 01 and 02 moved to G10, but it has no G10
arm. Other issues:

- "Each run records the git commit" is false; the field is empty.
- Section 5.2 quotes paper 08's old values.
- No cover letter or declarations; 20 uncited bib entries; a tracked
  `tmp-pdfcrop-50901.tex`.

Paper 10 explains the CR2 story better than any other paper in the
series and is the natural home for it.

**11. Dual-channel DGP. Grade C.** The 10-29 rerun (6ea37805e) cleared
every correctness finding: it uses the refined settings, a draw-for-
draw validation gate passes over 1,800 comparisons, the super-additivity
claim is withdrawn, and the result is reframed as one Gaussian law
(verified). The G9 arm still has 200 replicates, and its null rejection
reaches 0.085 (CO) and 0.080 (Hybrid) at $t_{1/2} = 1$ without comment.
It uses "original" and "extended DGP" and calls the parameters a
calibration, which paper 01 now says they are not. There is no front
matter. The open question is contribution. Paper 01's Discussion now
states the analytical core of paper 11: the correlation encoding is a
mean shift with a rank-one smaller residual covariance. Paper 11 adds
the simulation along the equal-total line. It also supplies a current
carryover-loss comparison: Hybrid relative loss 0.206 for the
correlation encoding against 0.085 for the mean shift (verified),
which could replace the 15.4% and 6.9% figure that has no live source.

**12. Closed-form design efficiency. Grade D/E.** Nothing has changed
since the 10-29 reject, and all ten major findings are open. The
centered $R^2$ is $-0.49$ (reproduced), and the ranking inverts for CO.
The evidence base has also moved under it. Paper 11's rerun replaced
`01-combined-summary.rds` with a file at N = 70 only, without the
`spec` column or the E1 arm, so the validation script can no longer run
(verified schema). Paper 01's closed-form paired-difference power now
meets the paper's aim. The one question paper 12 raises that no other
paper answers is the precision cost of the expectancy terms under an
open-label lead-in. That fits as a paragraph in paper 01's planning
section, computed with the GLS form.

**13. Two-component DGP, drop PB. Grade D, bordering E.** Nothing has
moved since the referee pass. The only commit, 588156b5d, adds a
`select` binding for CI; it did not change the data path. The Rmd still
reads `13-two-component-grid.rds`, which does not exist, and the
purled chunks fail at the first table (verified), so the document
cannot render. The text still reports paper 14's grid, and four of its
findings are false against its own data (verified). Its framing still
rests on a paper 01 comparison that paper 01 dropped. The driver still
uses `c.cfct = 0.1`, repair-enabled settings and model-based SEs. The
mean-shift estimate is biased about 13-20% at $t_{1/2} = 1$. The commit
message "fix CI renders of papers ... 13" overstates what it did.

**14. Two-component DGP, drop TV. Grade D.** Every number reproduces
from its file (verified). The argument has been overtaken in three
ways. The biomarker-TV parameters added on 2026-10-22 remove the safety
asymmetry that is half its case for dropping TV. The mean-shift arm's
small carryover cost is partly a misspecification bias of 13-22%
(verified). The paper 01 result it claims to reproduce is gone. Text
copied from paper 13 still argues the opposite choice, and the null
rejection rate (0.061 CO, 0.021 OL+BDC) is not diagnosed.

**15. Prognostic biomarker. Grade D (never refereed).** Paper 15 was
cut out of paper 02 on 11-07 and built on the published-set reruns. Its
strongest result is clean: under a placebo-response correlation the
G-family estimate is biased, and at N = 280 G9 rejects at 0.092 and
G10 at 0.089, while F1A0 rejects at 0.058 (verified). The expectancy
terms remove the bias. All 29 render-time claims hold (verified). Three
problems remain:

- The biomarker-by-time remedy, half the title's promise, is a
  `RESULTS PENDING` box, and G10 and F1A0 were not fitted under the
  trend pathway.
- Conclusion 2 overstates. At 5,000 replicates its own block N puts G9
  above the band in all four cells (0.059-0.074) and G10 in two of four
  (verified). The render-time checks test only the 500-replicate runs.
- F1A0B is offered as a recommendation, yet in the stress job it is
  conservative (0.030-0.034) and gives up 6 to 13 points of power to
  G10 at N = 70 (verified). It uses a single bootstrap resample.

A referee would also likely ask whether the expectancy remedy only
undoes a mechanism the simulator builds in, and why the paper is
separate from paper 01's Stage 6 (inferred). The literature is thin
(11 references). Every published-set run records an empty `git_sha`.
It moves to C once the biomarker-by-time block and the trend-pathway
refit are stored.

## Cross-cutting findings

1. **The recommended analysis differs across papers (checked).** Paper
   01 plans sample sizes for G10. Paper 02 recommends G10. Paper 15
   recommends F1A0 or F1A0B and does not recommend F1. Papers 08 and
   10 name G9. CLAUDE.md names F1, with F1A0B provisional. One
   recommendation, stated identically, is needed before 01, 02 and 15
   go out as a set.
2. **G10 null rates under a prognostic correlation disagree between 01
   and 15 (checked; cause unknown).** Paper 01's Stage 6 reports 0.108
   (constant natural history) and 0.110 (placebo) at 1,000 replicates.
   Paper 15 gives 0.038 to 0.062 for what appear to be the same
   pathways and settings. The papers use different drivers
   (`quick-sim/carryover-closed-form/11-strawman-stress-test.R` and
   `carryover-sensitivity/simulation-core.R`), and "prog-TV-grow"
   means a growing shift in one and a constant shift in the other.
   Paper 01's inflation factor for the robust model rests on its
   number. This must be resolved first.
3. **The mechanism of the LME conservatism is told two ways
   (checked).** Papers 07 and 08 attribute it to the working
   covariance. Paper 10 shows it is mainly mean-model misfit absorbed
   by REML. Paper 08 cites 10 for the opposite of what 10 says.
4. **Parameter sets and integration (checked).** Published-set reruns
   exist for 03, 06 and 10 but are untracked and not read by their
   manuscripts. Paper 02 stays on simcore2026 (disclosed). Papers 04,
   05, 07, 13 and 14 are outside the remediation.
5. **The cross-component correlation is still inconsistent
   (checked).** 05, 06, 08 and 09 use 0.1/0.05; 07, 13 and 14 use
   0.2/0.1. Only 01, 02, 10, 11 and 15 use the refined $c_1 = 0.2$.
6. **Run provenance (checked).** Container runs record an empty git
   commit (papers 10 and 15). Paper 05's artifacts record an empty
   `git_rev`, and the `nof1power` version was not bumped with the
   censoring fix.
7. **Paper 01 now contains the analytical core of 03, 11 and 12
   (inspected).** Appendix A.6 has 03's formula. The Discussion has
   11's conditional-structure result. The closed-form power covers
   12's aim.
8. **Mean-shift bias under a decaying exposure indicator (checked in
   09, 13 and 14).** It is handled in 09 by a binary-indicator arm and
   unhandled in 13 and 14.
9. **Stale companion references and side files (checked).** Bib titles
   for papers 01, 02 and 15 are out of date in every paper that cites
   them. Papers 03, 07, 11, 12, 13 and 14 still describe paper 01 as
   comparing the two encodings. Cover letters and title pages trail
   the manuscripts in 01, 02, 04, 06, 08 and 09. Paper 07 cites
   companions by repository path.

## Recommended consolidation

These are proposals, inferred from overlap in design, engine or
argument. None has been agreed.

- **11 into paper 01**, as an Additional file under a short subsection
  on how the biomarker enters. This changes `docs/66`, which proposed
  01 or 03.
- **12 retired**, with its one open question moved to paper 01's
  planning section.
- **13 and 14 merged** and rerun once with full, drop-PB and drop-TV
  arms. Given how little would remain, a supplement to paper 06 or 11,
  or a compendium section, may be better value than a standalone
  paper.
- **08 into 10**, making 10 the single calibration reference. This
  looks stronger now than on 10-29.
- **07 as a robustness appendix of paper 01**, rather than a note with
  08.
- **04 as a closed-form section of 05**, since they share an estimand.

Together these would reduce fifteen manuscripts to about nine.

## Suggested order of work

1. **Resolve finding 2 and settle one recommended analysis** (finding
   1). Both gate papers 01, 02 and 15.
2. **Papers 01, 02 and 15 as a set to *Trials*.**
   - 01: declarations and headings, side files, `appendix-b.Rmd`, the
     bib, a short re-review.
   - 02: decide whether to report the published-set reruns; replace
     the "one paired SE" criterion; side files.
   - 15: run the biomarker-by-time block and the trend-pathway refit;
     correct Conclusion 2; justify or drop F1A0B; fill in the
     literature.
3. **Commit and wire in the untracked published-set runs:** 10 (about
   a day, plus a G10 null arm of about 1.5 CPU-hours), 06 (with the
   cut to about 8,000 words), 03 (Results from the stored runs).
4. **The near-ready B papers:** 08 and 09 (to *Trials*), 04 and 07.
   Fix the mechanism sentences in 07 and 08 once 10 is settled.
5. **Reruns:** the `nof1power` engine fix and rerun for 05; the merged
   13/14 grid.
6. **Scope decisions:** 11 into 01, 12 retired, 07's placement.
7. **Paper 03** last, reframed around what the corrected study shows.

## Limitations of this assessment

The status rests on source and commit inspection and on recomputing
headline numbers from stored files, not on re-execution of any
simulation. No manuscript was rendered, so LaTeX and cross-reference
failures were not checked. Paper 15 has never been refereed, and nine
other papers were rebuilt or rerun after their latest referee report;
their grades may move on a fresh read. Several agents delegated part
of their work to helper agents and relayed the results; those checks
were not repeated. Journal word limits, abstract limits and required
declaration headings were not checked against the live author pages.
The cause of the paper 01 and paper 15 discrepancy (finding 2) was not
investigated. Word counts are rough.
