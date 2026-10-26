# Gap analysis and remediation plan: papers 01 and 02 for Trials

*2026-10-20 16:55 PDT*

**Author:** pmsimstats team

**Update, 2026-10-26.** Two further `simrun` jobs (28 and 30) and a
direct covariance check, all bearing on step 11's prognostic-
biomarker work and on the `G9`-versus-F1A0B decision.

- **Job 28 (`28-run-g9-prognostic-stress.R`, 32 cells x 500 reps)**
  crossed all four prognostic pathways with two confound strengths
  (0.3, 0.6), two $c_{bm}$ levels (0, 0.45) and $N \in \{35, 70\}$
  on Hybrid, fitting `G1`-`G9`, `E0` and the `F`-family on the same
  trials. Verified: `G9`'s Type I error is 0.042-0.050 at every null
  cell, including strength 0.6, and at strength 0.3 no method loses
  power relative to the no-confound alternative. At strength 0.6 on
  the prog-PB and prog-TV-const pathways every method, `G` and `F`
  alike, lost 15-25 points of power (`G9`: 0.597/0.853 on prog-PB
  and 0.502/0.772 on prog-TV-const at $N=35/70$, against 0.759/0.969
  with no confound). `G9` is at least as powerful as F1A0B in all 16
  alternative cells of this job.
- **Job 30 (`29-run-g9-random-effects.R`, 12 cells x 1000 reps)**
  was launched on the hypothesis that this loss was added patient-
  level variance that a richer random-effects structure might
  absorb. It was not: `G9` with a random slope on `t`, on `De`, or
  both matches plain `G9` to within one point of power and of SE in
  all 12 cells (paired per-replicate gains at most 0.008), at the
  cost of up to 2.5% non-convergence.
- **The strength-0.6 loss is a DGP artifact, not a confound
  effect.** The new check script
  `30-check-prognostic-confound-pd.R` (output
  `analysis/data/carryover-sensitivity/g9-prognostic-stress/pd-check.rds`)
  shows, with no simulation, that every strength-0.6 cell on prog-PB
  and prog-TV-const, and prog-BL at $c_{bm}=0.45$, has a non-
  positive-definite $\Sigma$ on all four Hybrid paths (minimum
  eigenvalue $-0.11$ to $-0.85$). Both drivers run with
  `stop_on_nonpd = FALSE` (the default retained in step 3), so
  `build_sigma_matrix()` repaired the matrix silently. The repair
  shrinks the population on-drug minus off-drug slope contrast of
  $S_x$ on the biomarker, the quantity `bm:Dbc` estimates, from
  $-3.53$ to $-2.38$ (prog-PB), from $-3.26$ to $-1.81$
  (prog-TV-const) and from $-3.26$ to $-3.00$ (prog-BL); the mean
  estimates jobs 28 and 30 report at those cells match the repaired
  values. Strength-0.3 cells are PD on every path, which is why they
  show no loss. prog-TV-grow is a mean shift and never affects PD.

Consequences. Job 28's strength-0.6 power numbers must not be cited
as the cost of a strong confound, and no analysis-side remedy for
them should be pursued; the premise of job 30 was wrong. `G9`'s
nominal size at those cells stands as a result but was obtained
under a repaired matrix. These cells are infeasible under the
refined process, which is exactly what `stop_on_nonpd = TRUE`
exists to refuse, and a valid "stronger confound" test needs a
PD-feasible strength ladder run with `stop_on_nonpd = TRUE` (not
yet done). One side observation, verified by the same script:
prog-PB induces a small nonzero population contrast even at
$c_{bm}=0$ ($-0.14$ at strength 0.3, $-0.27$ at 0.6, against $-3.3$
at the alternative), because the placebo-response SD scales with
expectancy and so differs on and off drug in Hybrid; `G9` showed no
detectable size inflation from it at these $N$, but it is the
mechanism the `F`-family's `De` term exists for and may matter at
larger $N$ or stronger confounds. For the `G9`-versus-F1A0B
decision: `G9` continues to hold nominal size under every
prognostic pathway tested, now at doubled strength, and was never
less powerful than F1A0B on the same trials in job 28. Neither
manuscript yet states any of this.

**Update, 2026-10-25.** Two compute jobs launched to close gaps this
note flagged as open on 2026-10-24: both are complete, fetched, and
applied.

- **Step 4 (the headline rerun) is now done.** `simrun` job 26 reran
  the unchanged `01-run-factorial.R` (540 cells, 500 reps) under the
  corrected `c.cf1t`/`c.cfct = 0.2` constants, superseding the
  2026-07-31 file (`01-factorial.rds`, now 1,080,000 rows, verified
  structurally sound after fetch/apply). The previously reported
  `G9`-best headline ranking has not yet been recomputed from this
  file; that recheck, and regenerating the headline table/figures
  from it, is the next step before this can be marked fully closed.
- **The E0 Type I error finding (2026-10-24, below) is now
  generalized across design and N, and is larger than first
  reported.** `simrun` job 27 ran the new driver
  `27-run-tier1-g-vs-f-e0.R` on the existing 30-cell Tier 1 figure
  grid (design x {CO, Hybrid, OLBDC}, N x {35, 70}, both Panel A and
  B cells), fitting `E0` alongside all of `G1`-`G9` and the `F`-family
  on the same simulated trials. At the eight null cells ($c_{bm}=0$,
  Panel A), `E0`'s rejection rate is:

  | Design | $N=35$ | $N=70$ |
  |---|---|---|
  | CO | 0.150 | 0.128 |
  | Hybrid | 0.070 | 0.060 |
  | OLBDC | 0.144 | 0.106 |

  This is a wider range than the 0.06-0.10 reported from S12's
  15-cell slice on 2026-10-24: CO and OLBDC run at roughly 2-3x
  nominal $\alpha=0.05$, while Hybrid is close to nominal. By
  contrast, every `G`-family and `F`-family specification stays in
  the 0.03-0.07 range across the same eight null cells (full table:
  `analysis/data/carryover-sensitivity/tier1-g-vs-f-e0/summary.rds`).
  The design-dependence (worst in CO and OLBDC, mild in Hybrid) is
  itself new information not in the 2026-10-24 note, and should be
  stated alongside the overall inflation once step 13 is drafted --
  it suggests the effect interacts with the design's expectancy
  structure, not just with the absence of `corCAR1`. Not yet decided:
  whether job 27's results become a new report block, an extension of
  S12, or are folded directly into step 13's synthesis discussion.

**Update, 2026-10-24.** This note's body (Parts 1-4) is unchanged
from 2026-10-21 below and is now stale on two Part 2 items; read
this update first, then the body for everything else, which still
holds. Changes since 2026-10-21, both in paper 02's Part 2 plan
table (updated in place below) and summarized here:

- **Step 11, partially done.** `docs/64` and paper 02's new Block S12
  (`report.Rmd`, added 2026-10-23/24) add at least one prognostic-
  biomarker null cell -- in fact all four of paper 01 Study 3's
  pathways (prog-BL, prog-TV-const, prog-TV-grow, prog-PB) -- to a
  grid that also fits every `G1`-`G9` specification. All nine hold
  nominal Type I error under every pathway without needing a `De`
  term, a result paper 02 did not have before. The other half of
  step 11, an expectancy-term variant of the `G1`-`G9` family
  itself, was not built; S12 compares the existing `G`-family
  (no `De`) against the separate `F`-family (which has `De`) rather
  than adding `De` to `G1`-`G9`.
- **Step 12, substantially done.** `F1A0B` (and `F1`/`F1A`/`F1A0`)
  are now tested alongside all nine `G`-family specifications and
  `E0` (Hendrickson's actual original analysis, newly added; see
  below), on the same simulated trials each replicate, using the
  corrected `c.cf1t`/`c.cfct = 0.2` constants from step 1 -- i.e. on
  paper 02's own corrected driver, closing the main concern in this
  step (F1A0B is no longer an untested, borrowed conclusion). Not yet
  done: the Methods text (S12) names `docs/62` by its internal
  document number rather than citing `docs/57`/`ch23` as the
  method's source, which is adequate for an internal working draft
  but is not a resolvable reference for a journal reviewer and needs
  converting to a proper citation (or removing the internal
  cross-reference) before submission. The same applies to every other
  `docs/NN` in-text reference S12 added.
- **Step 13, still open.** No explicit paper-01/paper-02/F1 synthesis
  statement has been written in either paper.
- **Step 4 (the headline rerun), still open and not to be confused
  with the above.** The data files behind paper 02's headline
  `G1`-`G9` table and figures (`analysis/scripts/carryover-
  sensitivity/output/01-factorial.rds` etc.) are dated 2026-07-31,
  before the step-1 `c.cf1t`/`c.cfct` fix (2026-10-21). The headline
  ranking (`G9` best-performing) is still reported from the pre-fix
  DGP. S12's new comparison is a separate grid (15 cells, not the
  full factorial) and does not supersede this gap.
- **New finding, not anticipated by this note's original gap list.**
  S12 found that `E0`, Hendrickson's actual original analysis
  (fetched directly from GitHub commit `58b32a9`, which omits
  `corCAR1` entirely and differs from this project's own vendored
  reference copy), runs Type I error at $0.06$-$0.10$ against
  nominal $0.05$ across every null cell checked, while every
  `G`-family and `F`-family specification in this project stays at
  $0.03$-$0.06$. This is evidence that dropping the residual-
  correlation structure this project added on top of the published
  method is itself consequential, a point neither paper currently
  states and that belongs in paper 02's Discussion once step 13 is
  drafted.

**Scope.** This note assesses how far paper 01
(`analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd`) and
paper 02 (`analysis/report/02-carryover-sensitivity/report.Rmd`) are
from being ready to submit to *Trials* (BioMed Central), and lays out
a sequenced plan to close the gap. It assumes each paper's stated
goal is the actual submission goal:

- **Paper 01's goal:** investigate defects in Hendrickson's published
  data-generating process (DGP) and recommend an improved DGP (the
  refined process) as the reference for aggregated N-of-1 simulation
  studies of this kind.
- **Paper 02's goal:** investigate candidate analysis specifications
  for the carryover problem and recommend one (currently
  Exposure-weighted, `G3`/`G8`, or possibly `G4`/`G9`) for use with
  that DGP.

Evidence labels follow the project convention: *verified* (observed
by running code or inspecting output this session), *inspected*
(read from source code or the manuscript text), *inferred*
(consistent with surrounding material but not directly checked), and
*unverified* (not checked at all). Items below are inspected/verified
during this session's review unless marked otherwise.

A caveat on the journal-specific checklist in Part 3: I was not able
to retrieve *Trials*' current Instructions for Authors page during
this session (it redirected to an authenticated Springer login), so
the journal-specific items are from general knowledge of BMC *Trials*
practice and should be confirmed against the live page before
submission.

---

## Part 1: Paper 01 (DGP paper)

### Where it stands

Paper 01 largely achieves its stated goal. The refined process is
fully specified (Table `tab:spec`), matches the project's canonical
definition in `docs/series/ch10-refined-process.md` exactly, and the
four defect-and-repair pairs (step vs. graded coupling, compound
symmetry vs. separable AR(1), compounding vs. anchored carryover,
silent repair vs. stop-on-non-PD) are each documented with a
closed-form or code-level explanation. The Recommendations and
Conclusions sections state the refined process as the paper's
proposed reference DGP. This is a coherent, complete piece of work
on its own terms.

### Gaps against a Trials submission

1. **Missing headline result the rest of the project cites as
   coming from this paper.** [verified] The project's own
   `CLAUDE.md` and `NOTATION.md` Part 2 cite a matched-N (N=70,
   Hybrid) relative power loss under carryover of 15.4% (covariance
   moderation) versus 6.9% (mean moderation) as "paper 01 results."
   That comparison does not currently exist anywhere in
   `report.Rmd`. It survives only in an archived pre-restatement
   draft (`archive/report-2026-10-07-pre-restatement.Rmd`) and a
   review document, after the manuscript's scope was narrowed on
   2026-10-07 to "specification of the compendium's reference DGP."
   Three other project documents point to a result that is not in
   the paper. This is the single highest-priority item: a referee
   who checks the cited figure against the manuscript will not find
   it.

2. **The package does not implement its own top recommendation.**
   [inspected, confirmed this session] `R/generateData.R` still
   silently repairs non-positive-definite covariance matrices via
   `corpcor::make.positive.definite()` rather than stopping, which is
   recommendation #1 in the paper's own Recommendations list
   ("Always apply the corrections: stop on non-positive-definite
   matrices..."). The paper's Appendix B code patches presumably
   implement the stop behavior for the paper's own analysis scripts,
   but the shipped package default does not. A *Trials* reviewer
   who reads the code-availability statement and checks the
   repository will find the package contradicting its own
   recommendation. Needs either a package-level fix or an explicit,
   visible statement in the paper that the stop behavior is
   implemented at the analysis-script level only, pending a package
   change.

3. **Companion-paper citations need a fabrication check.**
   [unverified] `report.Rmd` cites `@pmsimstats_dgp_architectures2026`
   and other in-house companion references. Given that this
   portfolio's referee-review campaign found fabricated-appearing
   content (invented seeds, unimplemented methods reported as
   implemented) in at least three other workspaces this year, every
   self-citation to an in-house companion paper or dataset should be
   checked for existence and content match before submission, not
   assumed correct because it is internal.

4. **Correction (2026-10-21): this item was wrong.** The original
   version of this gap analysis claimed the paper never names ADEMP
   explicitly. On a fuller read, `report.Rmd`'s "Simulation study
   design" section (around line 437) states plainly: "We designed
   and report the simulations following the ADEMP structure of
   Morris et al. [@morris2019simulation]: aims, data-generating
   mechanisms, estimands, methods and performance measures," with
   Table `tab:ademp` immediately after. This gap does not exist; no
   action is needed. Left here, struck through in substance rather
   than deleted, so the correction is visible rather than silently
   erasing the earlier error.

5. **Parameter calibration is disclosed as a limitation but not
   mitigated.** [inspected] The Limitations section already states
   plainly that Gompertz, correlation, and coupling parameters are
   not estimated from real trial data. This is honest and adequate
   as a limitation statement, but a *Trials* reviewer may ask for at
   least a sensitivity range or a cited external plausibility
   argument (e.g., prazosin pharmacokinetic half-life literature) 
   rather than a bare disclosure. Low priority relative to items 1-3.

### Plan

| Step | Action | Depends on | Status |
|---|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---|---|
| 1 | Decide: reinstate the matched-N architecture-loss section in paper 01 (rerun if the archived numbers are stale), or formally redirect `CLAUDE.md`/`NOTATION.md`/`ch07-architectures.md` to cite wherever that analysis now belongs. This is a scope decision, not a writing task, and blocks everything else that treats "paper 01 results" as settled. | none | **Done (2026-10-21).** User decided: redirect. The 15.4%/6.9% figure is out of scope for paper 01, whose goal is the covariance-moderation DGP, not the architecture contrast. |
| 2 | If reinstating: pull the matched-N comparison and z-statistics from `archive/report-2026-10-07-pre-restatement.Rmd`, verify the numbers against a fresh run under the current refined-process code (package and parameters may have drifted since 2026-10-07), and reintegrate as a results section. | 1 | **Not applicable.** Redirect was chosen over reinstatement. |
| 3 | Resolve the non-PD stop-vs-repair mismatch: either patch `R/generateData.R` to add a `stop_on_nonpd` argument used by default in the paper's own driver scripts, or add an explicit statement in Methods naming this as a package-vs-analysis-script distinction. | none | **Done (2026-10-21).** User decided: patch the package. `buildSigma()`/`generateData()`/`generateSimulatedResults()` now accept `stopOnNonPD` (default `FALSE`, backward compatible); when `TRUE` it raises an error instead of repairing, at both the explicit PD check and the previously-silent Cholesky-fallback repair. Tests added (`inst/tinytest/test_sigma.R`); full suite passes (79/79). Appendix B.5 and `CLAUDE.md`'s "Known open item" updated accordingly. |
| 4 | Verify every in-house companion citation (`@pmsimstats_dgp_architectures2026`, etc.) resolves to an existing, content-matched document. | none | **Done for paper 01 (2026-10-21).** Every `\@`-citation in `report.Rmd` (excluding bookdown's `\@ref()` cross-references, a false positive in the original check) resolves to a real `references.bib` entry; no fabricated citations found in paper 01 itself. `@pmsimstats_dgp_architectures2026`/`@pmsimstats_calibration2026` are cited only in paper 02, not paper 01; that remains paper 02's open item (Part 2, step 8) — both resolve to `@unpublished` placeholder entries for companion manuscripts that do not exist as actual report directories, which paper 02's remediation should address. |
| 5 | Add an explicit ADEMP-conformance statement in Methods, mapping the paper's existing table structure onto ADEMP's five elements by name. | none | **Not needed — already present.** See the correction to gap item 4 above; this was an error in the original gap analysis. |
| 6 | Full read-through against *Trials*' current Instructions for Authors (word limit, structured abstract, declarations, data/code availability statement, reporting-guideline checklist) once confirmed live. | 1-5 | **Blocked.** `trialsjournal.biomedcentral.com`'s submission-guidelines page redirects to an authenticated Springer login; could not be fetched this session. Needs a manual check against the live page. |

Paper 01's remaining gaps (steps 1, 3, 4, 5) are resolved as of
2026-10-21. Step 6 (the journal-specific read-through) is the only
item still open, and it is blocked on access to *Trials*' live
Instructions for Authors rather than on any further work in this
repository.

---

## Part 2: Paper 02 (analysis-method paper)

### Where it stands

Paper 02 is a substantial, well-scoped comparison: nine analysis
specifications (`G1`-`G9`), explicit scoping rationale for what was
and was not tested, a clear negative result (the lagged-treatment
term does nothing), a higher-precision paired-McNemar confirmation of
the headline ranking, and an explicit (if incomplete) check of the
post-selection inference risk in its data-driven half-life
specification. This is further from submission-ready than paper 01,
for reasons specific to this session's review plus the DGP
inheritance problem raised by the user's prompt ("02 should use the
canonical DGP").

### Gaps against a Trials submission

1. **The simulation driver does not fully implement the refined
   process paper 01 is recommending as the project's reference DGP.**
   [verified this session, see prior turns] Checked against the five
   refined-process properties:
   - Graded coupling: present, correctly implemented.
   - Anchored carryover: present, correctly implemented.
   - Separable AR(1) at $\rho=0.7$: present, but the **cross-factor
     correlation constant is numerically wrong**. The canonical form
     is a single $c_1 = 0.2$ scaled by $\rho^{|w_t-w_s|}$ at every
     lag; the driver (`simulation-core.R:444,839`) instead uses two
     different constants, 0.2 at same-timepoint and
     $0.1 \times \rho^{\text{lag}}$ at every other lag. This is a
     real numeric deviation from the documented spec, not an
     undocumented-but-correct case.
   - Stop on non-PD matrices: absent. The driver calls with
     `make_positive_definite = TRUE`, i.e., silent repair, matching
     the published process, not the refined one.
   - Randomized path allocation: absent. Patients are split across
     paths deterministically (`floor(N/n_paths)` plus remainder to
     the first paths), with no `sample()` call anywhere in the
     driver.

   This means paper 02's simulations are not run under the DGP that
   its own companion paper is about to recommend as the project's
   reference process. This is the direct answer to "02 should use
   the canonical DGP": it does not, on two of five properties, and
   is numerically off on a third.

2. **Stale architecture terminology throughout.** [inspected] The
   project retired the "Architecture A/B" labels for
   covariance-moderation/mean-moderation on 2026-09-30-ish per
   `NOTATION.md` Part 3, but paper 02 still uses "Architecture B
   (covariance moderation) only" in its DGP-architecture table row
   and in multiple figure captions ("Architecture A", "Architecture
   B"). A referee reading both papers side by side, as a *Trials*
   reviewer handling companion submissions likely would, will notice
   the terminology mismatch immediately.

3. **No explicit DGP-specification section.** [inspected] Unlike
   paper 01's Table `tab:spec`, paper 02's "Data-generating
   mechanism" section (lines 505-575) never names "refined process,"
   "anchored," "$c_1$," "graded coupling," or "randomized
   allocation." Even the two properties that are correctly
   implemented (graded coupling, anchored carryover) are described
   without using the project's own established vocabulary for them,
   so a reader cannot tell from the text alone whether this is the
   same DGP as paper 01's recommendation.

4. **The data-driven half-life specification's Type I error claim is
   under-supported and inconsistent with the project's own closer
   look at the same question.** [verified this session] Block S11
   checks `G4` (AIC-selected half-life, uncalibrated) and `G9`
   (`G4`+CR2) for Type I error inflation at one cell (Hybrid, N=70,
   true $t_{1/2}=1.0$ week, $n_{sim}=2000$) and finds none (0.038 and
   0.052 respectively), explicitly noting this "does not generalize
   beyond the single cell tested." The project's own deeper
   investigation of this exact question
   (`docs/57-data-chosen-carryover-half-life.md`,
   `docs/series/ch23-data-chosen-half-life.md`) found pooled
   inflation of 0.067-0.074 across 11 null cells using a richer
   analysis model (the full recommended `F1` model with expectancy
   terms) and a wider AIC grid, and explicitly states that even at
   its own $N=70$, $t_{1/2}=1$ cell, the naive-AIC excess was "not
   borne out cleanly." Paper 02's single-cell, simpler-model check
   disagrees with the closest comparable result in the project's own
   prior work, and the disagreement is not acknowledged. Either `G4`
   and `G9` need a broader Type I error check matching docs/57's
   pooled design before the power ranking that currently favors them
   (`G9` is reported as the single best-performing specification) can
   be trusted, or paper 02 should adopt the bootstrap-calibrated
   `F1A0B`/`F1AB` selection from `ch23` instead of the naive
   AIC selection, which was developed specifically to resolve this
   problem.

5. **Companion-paper citations need the same fabrication check as
   paper 01.** [unverified] `@pmsimstats_dgp_architectures2026` and
   `@pmsimstats_calibration2026` are both cited; verify both exist
   and match what is claimed of them.

6. **Correction (2026-10-21): this item was wrong, same error as
   paper 01's original item 4.** `report.Rmd` names ADEMP explicitly
   and repeatedly: "We embed these three in an ADEMP-structured"
   (line 268), "the ADEMP reporting framework of Morris, White, and
   Crowther" (line 410), "ADEMP order" (line 417), "ADEMP-compliant
   cell-level summaries" (line 2163), and "Morris ADEMP protocol"
   (line 2252). This gap does not exist; no action is needed. Left
   here, struck through in substance rather than deleted, so the
   correction is visible.

7. **Paper 02's "best method" claim is narrower than it reads, and
   F1A0B's provenance must not be attributed to paper 02.** [verified
   this session, see conversation record 2026-10-21] Three related
   problems, surfaced by discussing paper 01's Study 3 and the
   project's F1 analysis alongside paper 02:
   - Paper 02's `G1`-`G9` family never includes an expectancy term
     (`De`/`bmc:De`). Paper 01's own Study 3 (the 14-analysis
     comparison, Part 1 above) found that term necessary for
     robustness to a prognostic biomarker: every analysis that shared
     the biomarker's baseline association with its off-drug
     association failed at least one prognostic null cell; only the
     expectancy-term family held nominal size throughout. Paper 02
     has never tested a prognostic-biomarker cell at all, so its
     recommendation of Exposure-weighted (`G3`/`G8`) as the best
     method is unverified under exactly the condition paper 01 showed
     matters most.
   - The project's own named "recommended analysis," F1
     (`Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De`, CR2 cluster-robust
     SEs, Satterthwaite df; `CLAUDE.md`), is in effect the synthesis
     of paper 02's winning carryover representation (`Dbc` is
     identical to `G3`'s Exposure-weighted predictor) with paper 01's
     necessary expectancy terms. **Neither paper currently states
     this synthesis or cites the other for it.** A referee reading
     both companion papers, and then `CLAUDE.md`'s recommended
     analysis, would have no way to see that F1 is where their two
     separate conclusions actually meet.
   - `F1A0B` (the bootstrap-calibrated, data-driven half-life
     selection developed in `docs/57`/`ch23` specifically to resolve
     the Type I error risk in gap 4 above) has **never been tested in
     paper 02 or cited by it**. It must not be described, in either
     paper or in any summary of this work, as a conclusion paper 02
     reached after investigating `G1`-`G9` — it wasn't, and doing so
     would repeat the same kind of misattributed-result problem
     identified for the 15.4%/6.9% figure in Part 1. If F1A0B is
     eventually recommended, that claim needs its own test on paper
     02's (or a successor paper's) own grid before either paper states
     it.

### Plan

| Step | Action | Depends on |
|---|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---|
| Step | Action | Depends on | Status |
|---|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---|---|
| 1 | Fix the cross-factor correlation bug in `simulation-core.R`/`implementations/tidyverse/R/functions.R` so a single $c_1=0.2$ scales by $\rho^{\text{lag}}$ at every lag, matching the documented refined-process spec. | none | **Done (2026-10-21).** `simulation-core.R`'s two model-param constructors (`simulate_cell`, `simulate_cell_s7`) and `26-co-dbc-fix-paired.R` changed `c.cfct` from `0.1` to `0.2`, matching `c.cf1t`, so both constants now carry the documented $c_1=0.2$ at every lag. Full derivation, worked numeric example, and verified matrix-level effect in `docs/61-cross-factor-correlation-bug-paper02-driver.md`. Left unfixed, out of scope for this pass: the older one-off diagnostic scripts (`diagnose-null-type1*.R`, `check-null-architecture-equivalence.R`, `06-cr2-recalibration.R`, `diagnose-alt-cell-cr2.R`) and the earlier `quick-sim-prototype*.R` exploratory scripts still use the old constants (`0.1`/`0.05`); none of these feed the manuscript's reported tables or figures directly, but should be swept in a later pass if any of their output is ever cited. |
| 2 | Add randomized path allocation (replace the deterministic `floor`+remainder split with a `sample()`-based assignment) in the same driver. | none | **Done (2026-10-21).** `generate_data_multi_path()` in `simulation-core.R` now assigns the `N %% n_paths` remainder patients to a `sample.int()`-chosen subset of paths each call, instead of always the first `remainder` paths. |
| 3 | Decide, in coordination with paper 01's step 3, whether to fix non-PD handling (stop vs. repair) for this paper's simulations, or carry it forward as a shared, explicitly documented limitation in both papers. | Paper 01 step 3 | **Done (2026-10-21).** Paper 01 chose "patch the package"; the same patch was now mirrored into `implementations/tidyverse/R/functions.R` (`build_sigma_matrix()`/`generate_data()` accept `stop_on_nonpd`, default `FALSE`) and threaded through `simulation-core.R`'s `generate_data_multi_path()`/`simulate_cell()`/`simulate_cell_s7()` (same parameter name, same default). One behavioral note: `simulate_cell()`/`simulate_cell_s7()` wrap `generate_data_multi_path()` in a per-replicate `tryCatch`, so a `stop_on_nonpd = TRUE` error is caught and recorded as a non-converged replicate rather than halting the run — the appropriate adaptation for a Monte Carlo loop, documented inline, but different from the package's harder `stop()` semantics. Smoke-tested directly (confirmed: a known non-PD cell at `c_bm=0.95` returns `converged = FALSE` under `stop_on_nonpd = TRUE`); package `tinytest` suite still 79/79; parity tests (`test-parity-extended-tidyverse.R`) show the same 144-cell regression with or without these changes, confirmed pre-existing and unrelated (see note below). |
| 4 | Re-run the affected cells under the corrected settings from steps 1-3. Scope this to the cells that matter for the headline ranking (`G1`-`G9` at the Hybrid reference grid) rather than the full 216-cell grid if time-boxing is needed; flag any cell not rerun as a known limitation rather than silently leaving stale numbers in place. | 1, 2, 3 | **Done (2026-10-25).** The full factorial (`01-run-factorial.R`, unchanged script, 540 cells, 500 reps) was rerun on the simulation host under the corrected `c.cf1t`/`c.cfct = 0.2` constants (`simrun` job 26, 72a09e3, 5h39m; 1,080,000 rows, verified structurally after fetch/apply). This supersedes the 2026-07-31 file; the headline `G1`-`G9` table and figures should be regenerated from it and re-checked against the previously reported `G9`-best ranking before submission -- not yet done as of this update. |
| 5 | Add a DGP-specification section/table to paper 02 modeled on paper 01's Table `tab:spec`, naming the refined process explicitly and cross-referencing paper 01 as the paper that establishes and justifies it. This also resolves gap 3. | 1, 2, 3 (so the table describes what the rerun code actually does) | Open; code dependencies (1-3) are now satisfied, so this can be drafted before step 4's rerun completes, since it describes the DGP specification, not rerun results. |
| 6 | Replace "Architecture A/B" with "covariance moderation"/"mean moderation" throughout (table cells, all figure captions, prose). | none, can run in parallel | **Done (2026-10-21).** All 25 occurrences in `report.Rmd` and all 24 in `supplement.Rmd` (including the bare "under B"/"under A" shorthand found alongside them) replaced; both documents re-rendered cleanly. |
| 7 | Expand Block S11 to a pooled multi-cell check matching docs/57's design (multiple true half-lives, larger pooled n), or replace `G4`/`G9` with the bootstrap-calibrated `F1A0B`/`F1AB` method from ch23. The second option avoids re-deriving a size guarantee from scratch but requires implementing the bootstrap procedure in `simulation-core.R` and re-running `G4`/`G9`'s cells. | 4 (rerun baseline should be settled first) | Superseded by steps 11-12, which add F1A0B as a tested candidate directly; keep this step only if F1A0B is not pursued. |
| 8 | Verify companion citations (shared with paper 01 step 4). | none | **Done (2026-10-21).** Every `\@`-citation in paper 02's `report.Rmd` resolves to a real `references.bib` entry (no missing-key errors). But `@pmsimstats_dgp_architectures2026` and `@pmsimstats_calibration2026` are `@unpublished` placeholder entries for companion manuscripts that do not exist as actual report directories in this repo — a dangling-citation problem, not a missing-key one. See gap item 7: this is the same provenance issue as the F1A0B misattribution risk. These two citations need either a real companion manuscript to point to, or removal/rewording until one exists. |
| 9 | Name ADEMP conformance explicitly (shared pattern with paper 01 step 5). | none | **Not needed — already present.** See the correction to gap item 6 above; this was an error in the original gap analysis. |
| 10 | Full read-through against *Trials* Instructions for Authors once confirmed live. | 1-9 | Blocked, same as paper 01 step 6 — journal page not reachable this session. |
| 11 | Add an expectancy-term (`De`/`bmc:De`) variant to the `G1`-`G9` family, and add at least one prognostic-biomarker null cell to the grid (reusing paper 01 Study 3's cell design), before paper 02's "best method" claim is extended to cover a prognostic biomarker. | 4 | **Partially done (2026-10-23/24).** All four of paper 01 Study 3's prognostic pathways (prog-BL, prog-TV-const, prog-TV-grow, prog-PB) were added and run against all nine existing `G1`-`G9` specifications plus `E0`/the `F`-family (paper 02's new Block S12). All nine `G`-family specs hold nominal Type I error under every pathway *without* a `De` term -- a result paper 01 Study 3 did not establish for this specific family. The `De`/`bmc:De` variant of `G1`-`G9` itself was not built; S12 compares the unmodified `G`-family against the separate `De`-bearing `F`-family rather than adding `De` to `G1`-`G9` directly. **Extended (2026-10-26):** job 28 doubled the confound strength to 0.6 and added co-occurring-interaction cells; `G9` stays nominal throughout, but the strength-0.6 cells on prog-PB and prog-TV-const are non-PD and were silently repaired, so their power results are a DGP artifact (see the 2026-10-26 update and `30-check-prognostic-confound-pd.R`). A PD-feasible strength ladder under `stop_on_nonpd = TRUE` is still needed. |
| 12 | Add `F1A0B` (and/or `F1AB`) as a tested candidate alongside `G1`-`G9`, run on paper 02's own corrected grid (post steps 1-3). Implement it directly from `docs/57-data-chosen-carryover-half-life.md` Section 4.5 (the parametric null-model bootstrap that repeats the selection step in every resample) and `docs/series/ch23-data-chosen-half-life.md`'s construct box and candidate grid ($\mathcal{T}_a = \{0, 0.25, 0.5, 1, 2\}$ weeks); cite both explicitly in paper 02's Methods as the source of the method, not just of the idea. If F1A0B is recommended, this makes paper 02 itself the source of that conclusion, on its own tested grid, rather than an uncited internal document the recommendation is currently borrowed from without attribution. | 1, 2, 3, 11 | **Substantially done (2026-10-23/24), extended (2026-10-25).** `F1A0B` (plus `F1`/`F1A`/`F1A0` and `E0`, Hendrickson's actual original analysis) is tested alongside all nine `G`-family specifications on the same simulated trials each replicate, using the corrected `c.cf1t`/`c.cfct = 0.2` constants from step 1. S12's own 15-cell comparison (driver `06-g-vs-f-comparison.R`, `simrun` job 21) is now joined by a second, larger comparison on the existing Tier 1 figure grid (driver `27-run-tier1-g-vs-f-e0.R`, `simrun` job 27, 30 cells x 500 reps x design {CO, Hybrid, OLBDC} x N {35, 70}, 4h03m, fetched and verified: `analysis/data/carryover-sensitivity/tier1-g-vs-f-e0/summary.rds`). This generalizes the E0 finding below across design and N rather than the single 15-cell slice. Not done: S12's Methods text cites `docs/62` by internal document number rather than `docs/57`/`ch23` as a proper reference; every `docs/NN` in-text pointer S12 added needs converting to a citation (or removing) before submission; and job 27's results are not yet narrated in either manuscript (see new finding below). |
| 13 | Once 11-12 land, state the paper-01/paper-02/F1 synthesis explicitly in both papers: paper 02's Discussion should name `Dbc` as the carryover representation its own comparison selected, now combined with paper 01's expectancy terms into F1; paper 01 should re-verify, at least as a confirmatory check, that its headline refined-process power-vs-half-life results are consistent under F1/F1A0B, not only under `E9`/RH's `lmer`. This is a small addition to paper 01 (already otherwise resolved, Part 1) in addition to a paper 02 change. | 12 | Open. 11-12 are now far enough along (prognostic cells run, F1A0B tested) that this can be drafted; also fold in S12's new `E0`/`corCAR1` Type I error finding (see 2026-10-24 update above), which belongs in this synthesis discussion. |

Paper 02's gap is substantially larger than paper 01's. Steps 1, 2,
and 6 are mechanical and low-risk. Step 7 is the long pole: it either
requires a materially larger S11 rerun or a new bootstrap-calibration
implementation, and the choice determines whether `G4`/`G9` can keep
their current "best-performing specification" claim at all.

---

## Part 3: Cross-cutting dependencies

- **Paper 02 cannot credibly claim to use "the canonical DGP" until
  paper 01's own DGP recommendation is stable.** If paper 01's step 1
  (restoring or redirecting the matched-N architecture result)
  changes any parameter defaults or the stated ceiling values, paper
  02's rerun (its step 4) should happen after, not before, that
  settles, to avoid rerunning twice.
- **The non-PD stop-vs-repair decision (paper 01 step 3, paper 02
  step 3) is one decision, not two.** Whatever is decided should
  apply identically to both papers' code paths, since they currently
  share the same underlying `functions.R`/`generateData.R`
  implementations.
- **Terminology consistency.** Once paper 02's "Architecture A/B"
  labels are retired (step 6), confirm no other document in
  `docs/series/` or `NOTATION.md` still points to paper 02 using the
  old labels.
- **Citation verification (paper 01 step 4, paper 02 step 8)** should
  be done once, covering both papers' reference lists together, since
  they likely share several companion citations.
- **New finding, out of scope for this session: a pre-existing
  cross-implementation parity regression, unrelated to today's
  changes.** [verified 2026-10-21] Running
  `analysis/scripts/parity/test-parity-extended-tidyverse.R` (which
  checks numerical agreement between `original-extended` and
  `tidyverse` under matching seeds) reports 144 cells with a sigma-
  means mismatch between the two implementations, all under
  `mean_moderation` architecture. Confirmed via `git stash` that this
  regression is present identically with or without this session's
  code changes (same 144-cell count, same pattern), so it predates
  this work and is not something introduced here. It is a genuine
  open correctness question for the `tidyverse`/`original-extended`
  parity guarantee the project otherwise relies on, worth a dedicated
  investigation, but was not pursued further in this session since it
  sits outside both papers' plans above. Flagging it here so it is
  not lost.

## Part 4: Trials-journal checklist (unverified against the live page)

The following are general BMC *Trials* expectations from prior
experience with the journal, not confirmed against its current
Instructions for Authors during this session. Confirm each before
submission:

- Structured abstract (BMC journals typically cap abstracts around
  350 words; *Trials* may have its own structure requirement for
  methodology articles).
- A data and code availability statement naming the public
  repository and, if applicable, a release tag or commit hash
  corresponding to the submitted results, not a moving `main` branch.
- A declarations section (funding, competing interests, ethics/
  consent, as applicable; likely "not applicable" for a simulation
  study, but the journal will expect that stated explicitly rather
  than omitted).
- Explicit identification of the reporting guideline followed for
  simulation studies (ADEMP), if *Trials* lists this as an expected
  guideline for methodology submissions in its current instructions.
- Trial registration: not applicable for a simulation/methods paper,
  but confirm the journal's submission system allows marking this
  "not applicable" rather than requiring a registration number.
- Exploratory-study framing: per this project's own established
  convention, neither paper should describe its simulation grid as
  pre-registered or pre-specified; confirm final text keeps this
  framing (both papers appeared to do so as of this review, but this
  was not re-checked line by line during this session).

## Summary

**Paper 01, updated 2026-10-21: all five in-repository action items
(steps 1, 3, 4, 5, and the correction to the ADEMP finding) are
resolved.** The matched-N architecture-loss result was formally
scoped out of paper 01 (it belongs to the mean-vs-covariance
architecture question, not paper 01's actual goal of specifying the
best covariance-moderation DGP), with `CLAUDE.md` and
`ch07-architectures.md` corrected to stop citing "paper 01" for it.
The package now implements its own stop-on-non-PD recommendation via
an opt-in `stopOnNonPD` argument, tested and passing the full suite.
Paper 01's own citations were checked and all resolve to real
entries. The ADEMP-conformance gap in the original version of this
analysis was itself an error; the paper already names ADEMP
explicitly. Only step 6 (the live *Trials* Instructions-for-Authors
read-through) remains, blocked on access to that page rather than on
further work here.

Paper 02 was further out and unchanged as of the 2026-10-21 revision:
it had a real, numeric deviation from the DGP its companion paper
recommends, an under-powered and apparently-contradicted Type I
error claim supporting its currently-best-performing specification,
and pervasive stale terminology (see Part 2). Paper 02's own
companion citations (`@pmsimstats_dgp_architectures2026`,
`@pmsimstats_calibration2026`) were confirmed in that session to be
placeholder `@unpublished` bib entries for companion manuscripts that
do not exist as actual report directories — the fabrication-risk
check flagged as needed for both papers was therefore done for
paper 01 (clean) and partly done for paper 02 (two citations found
dangling; full code-level fabrication check still pending).

**Correction/update (2026-10-24): "unchanged" above is now stale.**
Steps 11 and 12 have moved substantially since 2026-10-21 (see the
2026-10-24 update at the top of this document) -- the DGP deviation
(step 1) was fixed that same day, the Type I error concern behind
step 7 now has a direct, shared-grid `G`-vs-`F`-vs-`E0` comparison
(paper 02's new Block S12) rather than remaining an unresolved
single-cell disagreement with `docs/57`, and prognostic-biomarker
robustness is now checked for every `G`-family specification. The
two dangling companion citations and the stale-terminology sweep
status are unchanged from 2026-10-21. The one item that has *not*
moved and is now the clear long pole is step 4: the headline
`G1`-`G9` table and figures still run on the pre-fix DGP.
