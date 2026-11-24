---
geometry: left=0.85in,top=0.8in,bottom=0.8in,right=0.85in
fontsize: 10pt
header-includes:
  - \linespread{1.15}
  - \setlength{\parskip}{0.3em}
---

# Submission readiness of the fourteen manuscripts

*2026-10-29 16:55 PDT*

**Author.** pmsimstats team

**Purpose.** A paper-by-paper statement of how far each manuscript under
`analysis/report/` is from journal submission as of 2026-10-29, what
blocks it, and what it would cost to clear. This document supersedes
`analysis/report/whitepaper-compendium-summary.md` (2026-08-12, eleven
papers) as the portfolio status record. For papers 01 and 02, the
step-level plan in `docs/60` remains the working detail; its status
column predates the 2026-10-27 to 2026-10-29 changes summarized here.

## Method and evidence

For papers 01 to 12, the source is the most recent internal referee
report in each manuscript directory (2026-10-27 to 2026-10-29). Each
major finding was checked against the commits made after that report.
A finding is marked addressed only where the change was seen in the
source or commit diff, and open only where no such change exists.
Papers 13 and 14 had never been refereed. Each received a full referee
pass today (`referee-report-2026-10-29.md` in each directory). That
pass included recomputing inline numbers from the stored result files.
No simulation was rerun and no manuscript was rendered for this
document. Effort estimates are rough, and grades are the reviewer's
judgment.

**Grades.** **A** submittable after copy-edit; **B** minor revision;
**C** major revision, no new simulation; **D** needs new or rerun
simulation; **E** incomplete draft (results section empty).

## Summary table

| Paper | Short title | Target | Grade | Effort to submit |
|---|---|---|---|---|
| 01 | Extended DGP | *Trials*, Methodology | B | 1-3 days |
| 02 | Carryover representation, G9 | *Trials* | C (D) | 1-2 weeks |
| 03 | Latent-class / mixture | SiM (assumed) | E | 2-3 weeks |
| 04 | Main effect, pilot | undecided | C (D) | 2-4 days |
| 05 | Hybrid power, main effect | SiM | D | 3-5 days |
| 06 | Three-component decomposition | SiM | C (D) | 3-5 days |
| 07 | Response-curve family | SiM, Methods Note | B | 1 day + re-review |
| 08 | Test-procedure choice | SiM | B | 1 day |
| 09 | Informative dropout | *Trials* (SiM: D) | B | 1 day |
| 10 | Interaction-test calibration | SiM | D | 2 days |
| 11 | Dual-channel DGP | SiM (assumed) | C | 1-2 days |
| 12 | Closed-form design efficiency | SiM, short comm. | D/E | 1 week |
| 13 | Two-component, drop PB | SiM (CSL only) | D | 3-4 days |
| 14 | Two-component, drop TV | SiM (CSL only) | D | 3-4 days |

A grade in parentheses applies if an optional or contested analysis
is required.

**Common to every paper.** The author roster is a placeholder
("Temperance Persons"). Affiliations, declarations and suggested
reviewers are incomplete. No paper has an archived, DOI-bearing code
release. None of these takes long, but they block every submission.

## Paper-by-paper assessment

**01. Extended DGP building on Hendrickson et al. Grade B.**
The 2026-10-27 referee asked for major revision. Revisions on
10-27 and 10-28 (3efaadf89, 23a144aa4) addressed trialist guidance by
adding a planning section and G9 sample sizes. They also corrected the
no-carryover wording and the Proposition scope. Partly open: the Stage
6 stress-test tables still come from the first run under recursive
carryover. The rerun used random path allocation, which is not part of
the refined process, so a cheap rerun under the final definition is
advisable. Administrative blockers:

- The declarations read "To be completed".
- There is no code DOI.
- `title-page.md` still carries the retired "Two Architectures" title.
- `appendix-b.Rmd` is an orphan file.

The G9 Figure 4B panel (job 41) predates the 2026-10-29
`censordata()` sign fix. The stored panel and the text agree, but the
committed driver would no longer reproduce it. The *Trials* author
instructions have not been checked against the live page. The science
is settled.

**02. Carryover representation and inference. Grade C (D if the
bm:t check is required).** The referee (10-27) asked for major
revision. Since then:

- **Addressed:** the size cost of half-life selection (about 0.055),
  coverage and calibration, the equal-size comparators (G6, G8) and
  isolation of E0 (job 42). The CR2 calibration in OL+BDC was
  corrected on 10-29 (8e810ca9d).
- **Partly open:** the prognostic-trend pathway inflates size to about
  0.07 for every analysis. The proposed remedy, a bm:t term, has not
  been evaluated.
- **Still tabled:** the strength-0.6 cells made under a repaired
  matrix (`docs/65`), with caveats.
- **Disclosed only:** generality.

Other blockers:

- The manuscript is about 1.5 times the length of paper 01.
- The cover letter and title page address SiM under an old title.
- The companion bib entry for paper 01 is stale.
- It cites the unpublished paper 10.

Paper 02 depends on paper 01 as its DGP companion.

**03. Latent-class and mixture formulations. Grade E.**
The referee (10-28) recommended rejection, with a new submission
considered after the production study. Commit 245e4369f addressed the
textual findings: DGP naming, inline numbers, the non-regular-testing
framing and attribution. The study itself is not done:

- The Results section is a `RESULTS-PENDING` marker (verified).
- Job 38's null-block output was never applied to `analysis/data`. It
  sits at a different path from the one the report reads.
- The alternative block has not been run.

Five major findings depend on that study. The paper names no journal.
Paper 01 cites paper 03, which should be softened while 03 is
incomplete.

**04. Treatment main effect, multi-crossover versus parallel group.
Grade C as a pilot note (D as originally scoped).** The referee
(10-28) recommended rejection. Commits 270a729d8 and 33fb58b1b (the
production run, job 45) addressed the design labels, the df and CR2
statement, the RNG streams, the deleted S6-S12 numbers and the
Raskind 2018 description. What remains is structural:

- The simulation mainly confirms a closed form (Senn).
- The Hybrid rebuild is declared not run.

A decision is needed on whether a short closed-form-confirming note is
worth submitting, and where. The cover letter (SiM) and title page
(*Trials* or *Clinical Trials*) disagree. The paper overlaps paper 05
in engine and estimand.

**05. Power of the Hybrid design for the main effect. Grade D.**
The referee (10-28) recommended rejection. The 10-29 full rerun (zeus
job 51, e238025cf) addressed the following:

- the S9 grid
- the design schedule
- RCT null bias (ANCOVA, size 0.054)
- non-PD repair
- disclosure of Type I error (alternating 0.094, crossover 0.148)

The headline comparison is still invalid. The `nof1power` engine
lacks a baseline-to-follow-up shift term, so alternating and crossover
power is wrong until the engine is fixed and those sweeps are rerun.
The engine also compounds carryover rather than anchoring it, and
Discussion numbers are hard-coded. Estimated effort is 1-2 days of
engine work, about 3 hours on the host, and 1-2 days of rewriting.

**06. Three-component decomposition. Grade C (D if Studies 1 and 4
are rerun).** The referee (10-28) asked for major revision. The repair
artifact and identity prediction are corrected (job 49, d0070b189,
eeda569ef). The displacement is 0.78 and 0.61 of the prediction on the
Hybrid designs and 1.20 on the balanced design, with no refused cells.
Hybrid recovery is partly addressed, and the confounding of N with
design remains open. Studies 1, 3 and 4 are June runs that predate the
anchored-carryover change, so the paper mixes two DGPs without saying
so (suspected to be immaterial to the bias results). The main cost is
length: about 25,000 words against roughly 8,000 for SiM. The README
still lists job 49 as not run.

**07. Response-curve family. Grade B.** The referee (10-28) recommended
rejection pending a common-random-numbers rerun. The rebuild
(0c0f1e260) retracted the old headline:

- Gompertz and logistic give identical decisions in all 1,000 seeds.
- tanh is 0.010 lower.
- Type I error is 0.034.

The numbers are current. The run uses $t_{1/2} = 0$, so the October
DGP changes do not apply. What weakens it is scope:

- a single near-null result
- the historical lme test rather than F1 with CR2
- the two-path OL+BDC design rather than the Hybrid

It has not been re-refereed since the rewrite. It is a candidate for
merging with 08.

**08. Test-procedure choice for the interaction. Grade B.**
The referee (10-28) asked for major revision. Commit 0bd472f4c fixed
the Kenward-Roger direction and Mancl-DeRouen with a $t$ reference
(null rates 0.047-0.059). It also added the missing literature and
narrated a CR2 arm (size 0.043-0.059, $\kappa$ 0.96-0.98). Remaining
issues:

- The model-based arm is G1 plus an expectancy main effect, not F1.
- With no carryover, the G9 recommendation goes untested.
- The cost of dichotomization is a known result, so novelty is narrow.
- It overlaps substantially with paper 10.
- The README still says the CR2 arm is pending.

**09. Informative dropout and design choice. Grade B for *Trials*
(D for SiM).** The referee (10-28) recommended rejection. A full
rebuild (8918e5c6e, then b13e1ec84 and 8e810ca9d) addressed the
following:

- the hazard formalization
- the worsening, improvement and return conditions
- removal of the OL design
- the true value, $\kappa$, coverage and CR2
- common random numbers
- inline numbers and the seed
- attribution

The `censordata()` sign fix does not affect its results, because its
driver does not call that function. It is still one parameter point.
Shared-parameter dropout is disclosed as not modeled. That is adequate
for *Trials*, as the referee suggested, but thin for SiM. The rebuild
has not been re-refereed.

**10. Interaction-test calibration. Grade D.** The referee (10-28)
recommended rejection. Commit 6d801ed8f rewrote the argument:

- the mechanism is mean-model misfit with the exact population
  $\kappa$
- the null-bias claim is withdrawn
- paired MCSEs, GEE efficiency, G9 and F1 are added
- the literature is added

The drivers now use the current DGP (`paper10-common.R`), but the
rerun has not been done. Every simulated table is from the
pre-correction run and is labeled as such. There are four
`RESULTS-PENDING` blocks (verified). Clearing it means running
`run-paper10-refined.sh` (about 10 CPU-hours) and re-reading the
conditional prose against the new numbers. There is no cover letter.
Paper 10 is the natural single methods reference for the CR2 story
that 02, 08 and 09 each repeat.

**11. Dual-channel DGP. Grade C.** The referee (10-28) recommended
rejection, with a reframed resubmission considered. Commit 6ea37805e
reran it under the current DGP, with a pre-exposure guard and a
draw-for-draw validation gate over 1,800 comparisons. It withdrew the
super-additivity claim and reframed the result as one Gaussian law. No
new simulation is strictly required. The G9 arm has only 200
replicates and should have 1,000. The main run's metadata records
uncommitted driver changes (disclosed). The open question is
contribution: once reframed, the result is "two parameterizations of
one slope". That may sit better as a section of paper 01 or 03 than
as a standalone paper. There is no front matter.

**12. Closed-form design efficiency. Grade D, bordering E.**
The referee report (10-29) recommended rejection until the closed form
is rebuilt on the GLS variance. Nothing has changed since (last edit
10-08), and all ten major findings are open:

- The reported $R^2$ is uncentered; the centered value is $-0.49$.
- The ranking inverts for CO.
- The attenuation factor models carryover the DGP lacks.
- The validation script is broken, because its July data were
  archived.
- The cover letter still claims "95%".

Paper 11's rerun cannot simply be swapped in as validation data.
Paper 01's exact E9 power and the compendium's GLS separation largely
supersede it. Folding it into one of those is likely cheaper than a
one-week rebuild.

**13. Two-component DGP, drop PB. Grade D (bordering E).**
The first referee pass found the following:

- **Wrong data file (verified).** The Rmd reads
  `13-two-component-grid.rds`, which no longer exists; the data now
  sits in `13-drop-pb-grid.rds`. Every hard-coded result matches
  paper 14's drop-TV grid, not paper 13's own.
- **Four findings false (verified).** Against its own data, four
  stated findings do not hold. The power floor at $c_{bm} = 0.45$ is
  0.790. CO does lose power to carryover. Mean moderation does lose
  power (about 0.08, $z \approx 2.7$). The cost ratio is about one
  half, not one sixth.
- **Stale settings.** The driver uses `c.cfct = 0.1` (the `docs/61`
  defect) and repair-enabled settings, though no cell was actually
  repaired. It also uses model-based SEs, with $\kappa$ between 0.91
  and 1.25.
- **Biased estimate.** The mean-moderation estimate is biased 13-20%
  at $t_{1/2} = 1$. The cause is that the mean shift is on drug only,
  while $D_{bc}$ assumes carryover.

Clearing it needs a rerun of the 108-cell grid (2-4 compute hours),
then rewriting Sections 3 and 5. The framing also needs re-anchoring:
paper 01 no longer makes the architecture comparison.

**14. Two-component DGP, drop TV. Grade D.** Every inline number
reproduces from `14-drop-tv-grid.rds`. The argument, however, has been
overtaken:

- **Premise no longer holds.** The paper argues that dropping TV is
  unconditionally safe because no biomarker-TV parameter exists. Three
  have existed since 2026-10-22, and `docs/65` shows the biomarker-TV
  and biomarker-PB feasibility bounds are identical.
- **Biased estimate (verified).** The mean-moderation estimate is
  biased 13.6-22.2% at $t_{1/2} = 1$ ($z$ from $-5.9$ to $-16.6$). So
  part of that architecture's small carryover cost is bias, which
  contradicts Appendix A.4.
- **Stale run.** It predates the current settings (`c.cfct = 0.1`,
  model-based SEs). The null rejection rate is 0.021 under OL+BDC and
  0.067 under CO, and neither is diagnosed.
- **Copied text.** Text taken from paper 13 contradicts this paper.

Papers 13 and 14 share a driver, grid, seed scheme and bibliography.
They should be merged before either is rerun.

## Cross-cutting findings

1. **Inconsistent cross-component correlation.** Papers 05, 06 and 08
   use $c_{cf1t}/c_{cfct} = 0.1/0.05$, and 07 uses $0.2/0.1$. Papers
   13 and 14 carry the `docs/61` value 0.1. None uses the refined
   process's single $c_1 = 0.2$. The values are disclosed, but the
   series is not internally consistent.
2. **Mean-moderation bias under the continuous drug indicator.**
   The bias found in 13 and 14 (mean shift on drug only, $D_{bc}$
   with carryover) likely applies wherever a mean-moderation arm is
   analyzed with $t_{1/2} > 0$. This is inferred, not checked for
   other papers.
3. **The CR2 story is told four times** (02, 08, 09, 10).
4. **Companion citations** point to unpublished or incomplete papers.
   Paper 01 cites 02, 03, 06, 08 and 10. Paper 02 cites 10.
5. **Stale side files.** README, title-page and cover-letter files
   trail the manuscripts in 01, 02, 04, 06 and 08.

## Recommended consolidation

The status reviews suggest five merges. Each is inferred from overlap
in design, engine or argument and has not been agreed:

- **13 + 14** into one two-component paper (rerun once).
- **08 into 10**, making 10 the single calibration reference.
- **07 + 08** as one robustness note, if 08 is not folded into 10.
- **12** folded into paper 01 or the compendium.
- **11** as a section of 01 or 03.

Together these would reduce fourteen manuscripts to about nine.

## Suggested order of work

1. **Papers 01 and 02 to *Trials* as companions.** First:
   - Clear 01's administrative items.
   - Decide on 02's bm:t check.
   - Remove 02's repaired-matrix rows.
   - Soften 01's citation of 03.
2. **Paper 09 to *Trials*** after front matter and a re-read (about
   1 day).
3. **Run the queued reruns on the simulation host:** paper 10 (about
   10 CPU-hours), then the merged 13/14 grid. Fix the 05 engine before
   its rerun.
4. **Decide scope for 04, 07, 08, 11 and 12** (merge, short note or
   shelve) before investing further rewriting.
5. **Paper 03** last. Its production study is the largest outstanding
   piece of work.

**Limitations of this assessment.** The 01-12 status rests on source
and commit inspection, not re-execution. Several items carry a
"suspected" or "inferred" label in the underlying agent reports, for
example the effect of the mixed DGP in paper 06. The *Trials* author
instructions remain unverified. No paper has been re-refereed after
its post-review rewrite. The grades for 07, 08, 09 and 11 in
particular may move on a fresh read.
