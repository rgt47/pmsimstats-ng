---
geometry: margin=2.4cm
fontsize: 10pt
---

# A Plan for Refactoring the pmsimstats-ng Whitepapers into One Series {.unlisted .unnumbered}
*2026-10-11 11:30 PDT; Section 10 decisions recorded 2026-10-11 16:50 PDT*

**Author.** pmsimstats team

**Purpose.** This document reviews the whole of `docs/` (58 numbered
whitepapers, 5 planning files, 6 other files) for consistency of
notation, methods and conclusions, and for redundancy, and sets out a
plan for rebuilding the set as one series that can be read from start
to finish. The series is meant to carry more detail than the
manuscripts under `analysis/report/`: every derivation, every table
and every number, each traceable to a script and an output file. The
plan states the design principles of the series, its Parts and
chapters, where every current file goes, the conflicts to settle
before content is moved, the gaps to fill, and the order of work. It
incorporates the state of the simulation work on 2026-10-11 and three
user decisions of the same day (Section 2.11) and the ten decisions
recorded in Section 10.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation and evidence

Notation follows `analysis/report/NOTATION.md` (revision 2026-10-09).
This plan introduces no new mathematical symbol. It uses four kinds of
local label, none of which is a symbol: chapter identifiers `ch01` to
`ch28`, front-matter identifiers `fm0` to `fm2` and appendix letters
`A` to `J` (Section 3); conflict identifiers CR-01 onward (Section 5);
gap identifiers GP-01 onward (Section 8); and decision identifiers
DQ-1 onward (Section 10). Old documents are cited as `docs/NN`, with
section numbers as they stand on 2026-10-11.

| Label | Meaning in this plan |
|------------|----------------------------------------------------|
| observed | read in the source document, a file listing or a command output on 2026-10-11 by the author of this plan |
| spot-checked | a survey claim re-read in the source at the cited section or line on 2026-10-11 |
| survey | reported by one of the four survey files listed below and not re-checked here |
| inferred | a judgment drawn from the pattern of what was read, not checked |
| decided | a user decision of 2026-10-11 (Sections 2.11 and 10) |

Table: Evidence labels used in this plan

The evidence base is four survey files written on 2026-10-11 that
between them read every file in `docs/` in full: survey a (docs/00 to
docs/20, including the `.tex` set, 09a and 19), survey b (docs/21 to
docs/33 and the planning files), survey c (docs/34 to docs/45) and
survey d (docs/46 to docs/58, the guides, the `pub_review` files and
the slide deck). Read directly for this plan: `NOTATION.md` in full;
the `\part`, `\section` and `\subsection` outline of
`analysis/report/compendium-mathematics.tex` (written CM below; 8,723
lines, 4 parts, 23 sections); the YAML titles of the 14 manuscripts;
the headings of the paper 04 working protocol
`analysis/scripts/treatment-main-effect/00-ademp-protocol.md` (written
P04 below); the project `CLAUDE.md`; `docs/57` Sections 1 and 8.1 to
8.2; `RERUN-anchored.md` Sections 1, 2 and 7; and the summary files of
the anchored reruns. Section 11 lists what was not checked.

## 1. Findings summary

### 1.1 State of the corpus

The top level of `docs/` holds 116 files besides this plan: 58
numbered documents, 5 planning files, 6 other files, and their PDFs
(observed, `find`). Four
subfolders (`archive/`, `figures/`, `references/`, `reviews/`) and a
`.claude/` settings folder sit beside them.

| Status (as judged by the surveys) | Count | Documents |
|----------------------|-----|---------------------------------------------|
| Current within its stated scope | 13 | 22 (parity workflow), 26, 35, 36, 37, 39, 47, 51, 53, 55, 56, 57, 58 |
| Partly superseded | 16 | 04, 06, 08, 19, 24, 25, 31, 32, 34, 40, 41, 42, 45, 48, 49, 54 |
| Superseded or historical record | 22 | 00, 01, 02, 03, 07, 07a, 09, 09a, 13, 14, 15, 16, 17, 18, 20, 21, 23, 28, 33, 46, 50, 52 |
| Planning | 4 | 22 (journal targets), 27, 38, 44 |
| Reference | 3 | 29, 30, 43 |
| Total numbered | 58 | numbers 05, 10, 11, 12 are absent; 22 is used twice; 07a and 09a are lettered |

Table: Numbered documents by status (survey; classification merged
across the four surveys)

| Format | Count | Notes |
|------------|-----|-------------------------------------------------------|
| `.tex` only | 16 | 00, 01, 03, 04, 06, 07, 07a, 08, 09, 13, 14, 15, 16, 17, 18, 20; none reviewed on 2026-10-07 or swept by the 2026-10-09 notation pass; only 07a has a PDF, and it predates the repository (December 2025) |
| `.md` and `.tex` | 2 | 09a and 19: the `.md` was reviewed and corrected, the `.tex` was not, and the two now disagree |
| `.md` | 39 | 21 to 58; each has a PDF at least as new as its source (survey, mtimes) |
| PDF only | 1 | 02, an annotated review copy of a deleted `.tex` |

Table: Numbered documents by format (observed, `ls`; review status
from the surveys)

Four facts about the simulation work frame the plan (observed on
disk, 2026-10-11, and from the coordinator's brief):

- The anchored stress-test main run (`simrun` job 15, server) has
  finished: `summary-full-anchored.csv` holds 451 rows, the 41 cells
  of `docs/46` times its 11 original analyses (S, C1, C2, A, ES, EA,
  F1, F2, R, DS, DA). By `RERUN-anchored.md` Section 1 the driver's
  default is anchored carryover with randomized path allocation, so
  this is the stress test under the refined process (inspected). No
  document reports it.
- The local anchored `--extra` (S, G1, G2) and `--extra2` (S, F1A)
  arms are not paired with job 15 (different BLAS and LAPACK draws; 0
  of 28,000 trials match). They are being rerun on mercury as jobs
  21 and 22 to pair with job 15 (decided).
- The local anchored main pass of the covariance study has finished
  (`summary-alt250-null500-anchored.csv`, 320 rows; fixed allocation,
  per `RERUN-anchored.md` Section 1). Its follow-up and $D_{bc,it}$
  passes were cancelled. No document reports the anchored main pass.
- The data-chosen half-life extension (job 14; 58 cells, 40,000
  trials, refined process) is reported in `docs/57` Section 8, with
  F1A0B as the provisional data-driven primary (spot-checked, `docs/57`
  Section 1, lines 209 to 223).

### 1.2 Classes of problem

**Notation.** The 16 `.tex` documents and the `.tex` halves of 09a and
19 predate every notation decision. They use the retired labels
"Architecture A" and "Architecture B", "N-of-1" for the Hybrid design,
"ER" and "TR" for $PB$ and $TV$, $S_{i,t}$ or `Sx` for $Y_{it}$,
$\rho$ for the compound-symmetry correlation (canonical $\rho_{cs}$),
$c_{bm}$ for the mean-moderation multiplier (canonical $\beta_{bm}$,
rule 2), and $D_{bc}$ for the binary indicator (rule 1) (survey a).
The documents 34 to 58, swept on 2026-10-09, are clean in substance
but carry regressions: reserved letters reused ($q$ in 35 and 36, $f$
in 38, 39, 41, 42, 45, $\Delta$ in 35, $M$, $m$, $g$, $h$, $\omega$
in 57, $B$, $G$, $K$, $Z$ in 58), plain $E(\cdot)$ in 37, $N(\cdot)$ in
41, one letter with several meanings across documents ($v$, $\omega$,
$C$, $\mu$), and one object with several symbols ($\tilde r$, $r$,
$\tilde g$ for the padded coupling vector; $n_{\text{on}}$ against
$n_2$ for row counts; $\sigma_v^2$ against $\tau^2$) (surveys c and
d). Section 7 lists them.

**Methods described differently.** The same quantity is quoted at three
parameter points without saying which: the Hybrid ceiling is 0.456
(`docs/34`, configuration F, $\rho = 0.7$), 0.291 (`docs/36`,
configuration F, $\rho = 0.8$), 0.440 (`docs/36`, configuration C,
$\rho = 0.8$) and 0.480 (`docs/35`, configuration C, $\rho = 0.7$),
each correct for its construct (survey c). The stress-test and
covariance-study numbers in 45 to 56 rest on the recursive carryover
and fixed allocation of the stored runs, while 57 Section 8 and job 15
use the refined process (survey d; observed). At least three
incompatible lists of differences between the package and the
published process exist (`docs/23` Section 2, `docs/32` Section 2,
`docs/54` Section 1), and NOTATION Part 3 gives a fourth (survey b).
Size bands differ (0.043 to 0.057 over 8 cells in 45; 0.0458 to 0.0542
in 43; 0.046 to 0.054 over 12 cells in 57), and evidence-label
vocabularies differ in at least six ways across 46 to 58 (surveys c
and d).

**Conclusion conflicts.** The most consequential, all detailed in
Section 5:

1. The recommended interaction analysis is stated five ways: the plain
   `corCAR1` model (24, 26, 31, `CLAUDE.md`), `us-KR` with
   $D_{bc,it}$ (45), F1 at a fixed half-life (49 Section 7, 50 Section
   2, 56), and F1A0B as data-driven primary with F1 as fallback (57
   Section 1, revised 2026-10-10).
2. The mechanism of the published carryover collapse: standard-error
   inflation with an unbiased estimate (03, 04, 09, 02) against
   attenuation because the coupling gap is zero (36 Section 5.3, 37
   Section 5.1, 55 Section 5).
3. The size of the architecture contrast: 10 to 15% against 40 to 60%
   relative loss (02, 19, `CLAUDE.md`) against 6.9% against 15.4% at
   matched total $N = 70$ (NOTATION Part 2).
4. Coupling ceilings: "0.25 to 0.49+", "about 0.6", "about 0.45"
   (03, 07a, 19.tex, `CLAUDE.md`) against 0.256 under compound
   symmetry and 0.474 (0.480 in Hybrid) under the separable form at
   $\rho = 0.7$ (34, 35, 52 Section 2).
5. The cost of F1 relative to S: 29% at the reference cell of the
   stored run (47 Section 1, line 123, spot-checked: S 0.514, F1 0.364),
   never re-estimated under the refined process, where F1 alone has
   0.438 (57 Section 8.3, line 1544, spot-checked). Job 15 now holds
   both.

**Redundancy.** Thirty-three clusters of overlap were found across the four
surveys (Section 6). The worst are: `docs/24` against `docs/31` Parts
I to III (near verbatim, each holding corrections the other lacks);
`docs/40` Section 9 against `docs/42` Sections 2 and 4 (identical
tables and figures); `docs/46` Amendment 2 against `docs/47`; the
"superseded by F1" paragraph repeated in six documents (35, 41, 42,
43, 44, 45); the description of the published process in nine places
(24, 25, 28, 31, 36, 37, 40, 44, 53); the ceiling formula and tables in
eight (32, 34, 35, 36, 37, 44, 54, CM); package implementation status
maintained in at least five (34, 35, 36, 44, 52), each with its own
inspection date; and `docs/23`, which copies 21 earlier documents.

**Convention lapses.** Every `.tex` document has a hand-written
"Rendered on" footer with the obsolete path `~/prj/alz/10-pmsimstats-ng`
and no notation index, Purpose or evidence labels (survey a). Further
lapses: "Working draft" in the header of `docs/36`, the most-cited
document (cited by 20 others); a word the user's style rules forbid in
`plan-coherent-analysis-comparison.md` line 20 and in the `pub_review`
files; British spellings in `reorg-plan-2026-04-22.md` and in
`CLAUDE.md`; "patient(s)" in 37, 38 and 39 against "participant";
"hybrid" in lower case in 28 and 31; the file name of `docs/51`
("ten") against its seventeen features; two numbers 22; swapped figure
files in `docs/40`; one orphan figure (`36-fig4-contrast-reps100.png`,
cited nowhere) (surveys a to d).

**Stale PDFs and `.tex`.** No current PDF exists for any of the 16
`.tex`-only documents; 07a's PDF is from 2025; dated duplicate PDFs of
29 and 30 from June remain; `docs/18` cannot render (its figures under
`analysis/figure5/` are gone); `docs/02` is a PDF of a deleted source
(survey a; survey b; observed).

## 2. Design principles for the refactored series

### 2.1 Audience and purpose

The series is written for a methodologist who has to understand,
reproduce or extend the work: a statistician joining the project, a
referee of a manuscript, or one of the original authors. It assumes
graduate statistics and no prior knowledge of the project. It is the
long form of the program. Where a manuscript states a result and
cites its derivation, the series gives the derivation steps, the
construct, every cell of every table, the Monte Carlo error and the
command that regenerates the number. It is not a manuscript: it
records the history of decisions in one appendix rather than
throughout, and it can be longer than any journal allows.

Three reading paths are offered in the reading guide (`fm0`): the full
path (Parts I to VI); the practitioner's path (`ch01`, `ch02`, `ch10`,
`ch22`, `ch23`, `ch25`); and the original authors' path (`ch03`,
`ch04`, `ch26` to `ch28`).

### 2.2 Relation to the manuscripts

Each chapter names the manuscripts that use it and the manuscript
sections that depend on it. The rule is one direction of detail: a
manuscript may summarize a chapter, and a chapter never depends on a
manuscript for a number it does not itself tabulate. Every number in a
chapter carries a source line (script, output file, commit) in its
table caption or in the chapter's "Reproduce" section. Where a
manuscript reports a result the series reproduces, the chapter states
whether the numbers agree and, if not, why (for example the published
power at $t_{1/2} = 0.1$ week, quoted as 0.126, 0.078, 0.100 and 0.076
in 37, 53, 55 and paper 01; reconciled in 55 Section 7.3, survey d).

### 2.3 Relation to the mathematics compendium

Three placements were considered.

1. **The compendium as the series' mathematical volume (decided
   2026-10-11, DQ-1).**
   CM stays one LaTeX document and becomes Volume II of the series.
   Chapters state each result, give the derivation in enough steps to
   follow the argument, and cite the CM subsection label for the full
   proof and its verification script. Derivations that exist only in a
   whitepaper (57 Sections 3 to 4, 58 Propositions 1 to 7, 55 Section
   5) are moved into CM once, and the chapter keeps the argued form.
2. **A separate reference.** CM stays as it is, cited ad hoc. This is
   the status quo; it leaves the reader without a path through CM and
   lets the two drift (CM's main-effect section still describes
   version 1 of paper 04, and CM cites none of 55 to 58; survey d).
3. **Folded into the chapters.** Each CM section moves into the
   chapter it serves. This gives one document per topic but destroys
   CM's LaTeX structure, its equation labels and its pairing with the
   `analysis/scripts/compendium-mathematics/check-*.R` verification
   scripts, and it is the largest and riskiest single move of the
   whole refactor.

Option 1 was recommended, and adopted, because the main defect the
surveys found is
redundancy, and two full copies of each derivation would recreate it;
because CM already has a reading order (its Section "Scope and reading
order") and verification scripts that a Markdown chapter cannot hold;
and because 46 of the 118 files that cite `docs/NN` outside `docs/`
are under `analysis/report/` and already cite CM labels (observed,
grep). The chapter template (Section 2.10) gives each chapter a
"Derivations" section with a row per result: statement, CM label,
verification script. The cost of option 1 is that CM must be brought up
to date (GP-28) and its Part order aligned with the series Parts
(Section 3 shows that the existing order already matches closely).

### 2.4 Relation to NOTATION.md

`analysis/report/NOTATION.md` stays the single source for symbols,
code identifiers and named processes. The series adds nothing to it
without a NOTATION revision. Each chapter's notation block lists only
its local symbols, with a status column (canonical or local), and links
every canonical symbol to NOTATION by its table name. A local symbol
may not use a letter on NOTATION's reserved list. The front-matter
chapter `fm2` holds the project vocabulary that NOTATION Part 3 holds
for the manuscripts (named processes, configurations, analysis IDs,
evidence labels), so that the lexicon (Appendix E) can stay a
cross-field glossary.

### 2.5 One evidence-label vocabulary

The series uses seven labels, defined once in `fm2` and used in every
chapter:

| Label | Use |
|------------|-----------------------------------------------------|
| verified | computed or run, and the output observed; cite the script |
| derived | proved analytically; cite the CM label or the chapter section |
| inspected | read in code or in a source document; cite file and line or commit |
| quoted | taken from another document; cite document and section |
| argued | a reasoned position, not a proof |
| inferred | consistent with the surrounding pattern, not checked |
| unverified | not checked at all |

Table: The evidence-label vocabulary of the series (merging the sets
of `docs/51`, `docs/53`, `docs/54`, `docs/55`, `docs/56`, `docs/57`
and `docs/58`)

### 2.6 One Monte Carlo methodology

`fm2` states, and every chapter applies, the following conventions,
collected from 45 Section 2, 46 Section 3, 57 Section 8.1 and CM
"Shared performance measures":

- **Reporting frame.** ADEMP for every study (aims, data-generating
  mechanisms, estimands, methods, performance measures), with the
  construct box of Section 2.10.
- **Size band.** For a pooled null rate over $n$ null trials the band
  is $0.05 \pm 2\sqrt{0.05 \times 0.95 / n}$, stated in the caption
  with $n$; for 12,000 trials, 0.046 to 0.054. Per-cell cap 0.075;
  a prognostic cell above 0.075 raises a flag. These are the screens of
  `docs/46` as amended (quoted, 57 Section 8.1).
- **Distances from nominal.** "$k$ MCSE above 0.05" uses the MCSE at
  0.05, $\sqrt{0.05 \times 0.95 / n_{\text{reps}}}$, as `docs/42`
  Section 6 does; the convention of `docs/40` Section 9 (MCSE at the
  observed rate) is retired (resolves CR-26).
- **MCSE with every number.** Every rejection rate carries its binomial
  MCSE $\sqrt{\hat\pi(1-\hat\pi)/n_{\text{reps}}}$; every pooled rate
  its pooled MCSE; every difference between independent runs
  $\sqrt{\text{MCSE}_1^2 + \text{MCSE}_2^2}$ (`RERUN-anchored.md`
  Section 7); paired comparisons within a run use McNemar.
- **Warp-speed bootstrap.** Rates from one resample per trial
  (Giacomini, Politis and White, 2013) are labeled warp-speed; their
  MCSE is computed by resampling trials and is larger than the
  binomial one (0.0097 to 0.0112 against 0.0069 per cell at
  $n_{\text{sim}} = 1{,}000$; quoted, 57 Section 8.1). The binomial
  MCSE is never quoted alone for a warp-speed rate.
- **Size-adjusted power** follows CM "Size-adjusted power"; maximin
  shortfall states its cell set in the caption (non-prognostic
  alternative cells, or all 26 alternative cells), so the two readings
  of `docs/43` and `docs/47` Section 7 stay apart (CR-29).
- **Replicate bands.** Null cells $n_{\text{reps}} \ge 1{,}000$ and
  alternative cells $\ge 500$ for any number that enters a decision;
  smaller runs are labeled "illustration" or "pilot" in the caption.
- **Sample size.** $N$ is always the total (NOTATION rule 4), and every
  study table states the allocation (fixed or randomized).

### 2.7 One analysis-ID scheme

The stored interaction-track labels (S, C1, C2, A, ES, EA, F1, F2, R,
DS, DA, G1, G2, F1A, F1A0, F1A0B, F1AB, MA, MaxT and their suffixes) are
keys in the stored stress-test files, and NOTATION Part 2 forbids a
display-time mapping layer. The scheme (decided 2026-10-11, DQ-5)
therefore keeps them and changes what is new or local. A document
names the track at first use ("interaction analysis F1", "main-effect
analysis `M-F1`").

| Class | Form | Examples | Change |
|----------------------|----------|--------------------|----------------------|
| Interaction track (stress test, covariance study selection, data-chosen half-life) | bare label | S, F1, F1A0B, G2 | none; the default track of Parts III and IV |
| Main-effect track (P04, `docs/58`) | `M-` prefix | `M-S`, `M-F1`, `M-VF`, `M-RV`, `M-F1A` | rename in P04, `docs/58` content and `ch24` before paper 04's stored data exist (DQ-5, decided) |
| Covariance study analyses (`docs/45`, `docs/43` Sections 3.7 and 3.8) | `CV-` prefix | `CV-V1` to `CV-V7`, `CV-F1` to `CV-F7` | removes the collision of catalog F1 with stress-test F1 and configuration F |
| Paper-local catalog IDs | paper prefix, as NOTATION's "05-S2" rule | `02-G1` to `02-G9`, `06-C1` to `06-C4`, `08-T1`, `10-M0` | removes collisions with G1, G2, C1, C2 and configuration C |
| Published analyses (`docs/43` Section 3.1) | `PUB-` prefix | `PUB-O1` to `PUB-O5` | none in substance |
| Summary measures and time adjustments (`docs/43` S1 to S3, D1 to D7) | `SM-`, `TA-` prefix | `SM-1`, `TA-3` | removes collisions with sensitivity sweeps S1, S2 and with $D_{it}$ |
| Exposure codings | E1 to E9, E9b, E`n`cr2 | unchanged | NOTATION extended E series |
| Constructs | "configuration A" to "configuration F", always spelled out, with $\rho$ | "configuration C at $\rho = 0.7$" | never a bare letter |

Table: The analysis-ID scheme proposed for the series

The registry of all IDs, with formula, track, source and stored-data
key, is a machine-readable table (`docs/series/registry/analyses.csv`)
from which Appendix A is generated.

### 2.8 Canonical names

- **The published process**: Hendrickson's `generateData()` and
  `buildtrialdesign()` at commit `58b32a9` (NOTATION Part 3). Never
  "the original code", "orig" or `42ac030` for the published process;
  `42ac030` is the initial commit, whose `generateData.R` equals
  `58b32a9`'s, but analysis defaults may differ (survey a, memory note).
- **The 2024 revision**: commits `8609f12`, `f6ee86d`, `f70f86d`,
  `325314f` (scale factor $s_f$, ratio coupling, timepoint-1 guard,
  `Dbc` in the analysis). **The vendored copy**: the 2024 state copied
  as `06dac83`, which `docs/54` calls `orig`.
- **The refined process**: as NOTATION Part 3; with the user decision
  of Section 2.11 its validated fixed-half-life analysis uses
  $t_{1/2}^a = 0.2$ week.
- **Configurations A to F**: as NOTATION Part 3, always with $\rho$.
  `covar` is replaced in prose by "configuration F at the package's
  parameters"; `orig` by "the published construct" or "the vendored
  copy", whichever is meant.
- **The published construct**: configuration A only (compound symmetry
  with $\rho_{cs} = 0.8$, constant cross block). AR(1) variants of the
  `04-power-simulation.R` construct are named by configuration
  (resolves CR-30).
- **Architectures**: mean-moderation, covariance-moderation,
  dual-channel. **Designs**: OL, CO, OL+BDC, Hybrid (capitalized).
- **The estimand**: the biomarker-treatment interaction $\beta_{bm:D}$
  (docs/30 D3); "biomarker-by-drug" in the titles of 43 to 57 becomes
  "biomarker-treatment" (DQ-7, decided).
- **Participant**, never patient; **visit**, never occasion;
  **off-drug period**, never washout for the span after discontinuation.

### 2.9 Superseded material: the dated-note policy

The current documents carry dozens of inline "(Superseded 2026-10-NN
...)" paragraphs; one paragraph about F1 recurs in six documents
(survey c, cluster 14), and `docs/57` holds three successive
recommendations. The series adopts a different rule.

1. **Chapters state the current answer only.** When a chapter is
   written, superseded statements are rewritten, not annotated. A
   chapter contains no dated correction paragraph at release.
2. **History lives in two appendices.** Appendix B (decision history)
   records each decision once: date, what was recommended, what
   replaced it, the evidence and the old document. Appendix C (errata)
   records each correction to an old document or a manuscript once,
   with its status.
3. **After release, changes are logged, not annotated.** Each chapter
   ends with a change log (date, section, change, reason). A dated
   inline note is allowed only where a result a reader may already have
   used changes sign or magnitude materially, and it points to the
   change log.
4. **Old documents are frozen.** Archived documents are not edited
   further, except for a banner at the top naming the chapter that
   replaces them.

### 2.10 Chapter template

Every chapter follows one template:

1. Title; timestamp line; **Author.** pmsimstats team; **Purpose**
   (question answered, audience); **Prerequisites** (chapters).
2. Contents block, as in this document.
3. Notation: local symbols only, with status; a pointer to NOTATION.
4. **Summary**: the results, each with its evidence label.
5. **Construct box**: a fixed table stating, for every numeric result
   in the chapter, the process (published, refined, configuration X at
   $\rho$), carryover (anchored or compounding), coupling rule,
   allocation (fixed or randomized), design, $N$ (total), analysis ID,
   $t_{1/2}^a$ and $n_{\text{reps}}$.
6. Setting and argument (the body).
7. **Derivations**: one row per result with its statement, CM label
   and verification script.
8. **Numerical evidence**: tables with MCSE; each caption gives the
   script and output file.
9. **Relation to the manuscripts**: which manuscript sections depend on
   the chapter, and whether their numbers agree.
10. **Limitations and open items**, each linked to the open-items
    register of Appendix D.
11. **Reproduce**: commands, one per output.
12. References (from the series bibliography); change log.

### 2.11 User decisions of 2026-10-11

Three decisions were made while this plan was written and are treated
as settled:

1. **Fixed analysis half-life.** The validated fixed-half-life fallback
   uses $t_{1/2}^a = 0.2$ week in both tracks, interaction and main
   effect. The value 0.5 week is described as the historical base
   setting of the stress test and the covariance study. Evidence for
   F1 at 0.2 week comes from the job 14 extension (F1h02: pooled null
   rate 0.050 over 12 non-prognostic cells, worst cell 0.058, worst
   prognostic cell 0.064; quoted, 57 Section 8.2). The 41-cell stress
   test (stored run and job 15) and the covariance study evaluated F1
   and $D_{bc,it}$ at 0.5 week only (survey; observed for job 15's
   file layout), so chapters 20 to 22 must say which half-life each
   number used. Code and driver defaults set at 0.5 week (survey d,
   `docs/56` Section 8) are aligned in Phase 2 (Section 9).
2. **The older `.tex` whitepapers** (00 to 20, with 02, 09a and 19)
   are consolidated into one document that holds only results not
   covered by the newer whitepapers, with the suspect references
   checked or dropped; the originals then go to the archive. This
   document is Appendix I, the legacy record (Section 3.8).
3. **The `--extra` and `--extra2` arms** of the stress test are being
   rerun on the server mercury (jobs 21 and 22) to pair with job 15.
   (Done 2026-10-15: jobs 21 and 22 are paired with job 15 and reported
   in `docs/47`, with the post hoc job 37; F1 retained.)

## 3. Target structure

The series has front matter, six Parts with 28 chapters, ten
appendices and the mathematics compendium as Volume II. Lengths are
estimates in letter-size PDF pages at the current rendering (for
scale, `docs/37` is 1,673 lines and renders to 43 pages; observed).
The total is about 670 pages, against roughly 50,800 lines of `.md`
and `.tex` in `docs/` today (observed, `wc`). Much of the reduction in
duplicated text is offset by tables that are now split across
documents being given in full once. In the tables below, "Needs" lists
the prerequisite chapters and "pp." the estimated length in pages.

### 3.1 Front matter

| ID | Title | Needs | pp. |
|----|--------------------------------------------------|--------------|----|
| `fm0` | Reading guide and index | none | 6 |
| `fm1` | Executive summary | none | 4 |
| `fm2` | Notation, vocabulary and conventions | none | 8 |

Table: Front matter of the series

**`fm0`. Reading guide and index.**

- Scope: three reading paths; status of every chapter; old-to-new map of every
  `docs/NN`.
- Sources: `docs/00`, `docs/23` Section 0, CM "Scope and reading order",
  Section 4 of this plan.
- New to write: all; generated from the registry
  `docs/series/registry/chapters.csv`.

**`fm1`. Executive summary.**

- Scope: the program's results in four pages: the published process and its
  defects, the refined process, the closed form, the recommended analysis, the
  main effect.
- Sources: `docs/50` (structure), 36 Section 1, 37 summary, 56 Section 1, 57
  Section 1, NOTATION Part 2 (paper 01 at matched $N$).
- New to write: all, after Phase 1.

**`fm2`. Notation, vocabulary and conventions.**

- Scope: pointer to NOTATION; named processes, configurations, analysis IDs;
  evidence labels; Monte Carlo conventions; construct box.
- Sources: NOTATION Parts 1 to 3, `docs/29` (project terms), `docs/30`,
  Sections 2.5 to 2.8 of this plan, CM "Shared performance measures".
- New to write: the project glossary and the conventions.


### 3.2 Part I. The problem and the published process

| ID | Title | Needs | pp. |
|----|--------------------------------------------------|--------------|----|
| `ch01` | The question: a predictive biomarker in aggregated N-of-1 trials | `fm2` | 12 |
| `ch02` | Trial designs, randomization paths and exposure clocks | `ch01` | 14 |
| `ch03` | The published process at commit `58b32a9` | `ch02` | 22 |
| `ch04` | Three defects of the published process | `ch03` | 18 |

Table: Part I of the series

**`ch01`. The question: a predictive biomarker in aggregated N-of-1 trials.**

- Scope: the estimand $\beta_{bm:D}$; predictive against prognostic; the sign
  convention; why carryover matters; a one-page history of the program.
- Sources: 24 "Big picture"; 40 Sections 1 to 2; 31 Part IV (Sections 14 to
  16); 02 Sections 1 and 7; 21 "Biological Assumptions"; papers 01 and 06
  introductions.
- New to write: the program history; a framing that answers the co-author's
  comments on 02 pages 2 and 3.

**`ch02`. Trial designs, randomization paths and exposure clocks.**

- Scope: OL, CO, OL+BDC, Hybrid, balanced placebo, cycle designs; schedules
  $w_t$; paths $P$; $\eta_t$; $t_{od}$, $t_{sd}$, $t_{pb}$ and their restart
  rules; $N$ as total and the 2026-08-08 matched-$N$ correction.
- Sources: 13 Section 2 (after a check against `58b32a9`); 33 Section 2.5
  (design table); 40 Sections 2 to 3; 24 Component 2; 37 Section 2; 27
  Sections 2a and 2c (cycle designs, Hybrid round trip corrected); NOTATION
  Part 2; CM "Trial designs and treatment exposure", "Time arguments of the
  components".
- New to write: one table per design; a Hybrid schedule figure.

**`ch03`. The published process at commit `58b32a9`.**

- Scope: `generateData()` and `buildtrialdesign()` as published; parameters;
  compound symmetry $\rho_{cs} = 0.8$; constant cross block; step coupling;
  compounding recursion; silent repair; the published `lmer` analysis with
  binary $D_{it}$; the censoring grid; Figures 4 and 5; the commit attribution
  table.
- Sources: 40 Sections 3 to 6 and 8; 36 Section 3; 53 Sections 3 to 4; 55
  Sections 3 to 4; 44 Section 2; 25 Section 4; 23 Section 3.1; CM "The
  published compound-symmetry construct" and the three "2024 revision"
  subsections; Appendix I Sections 1 to 3.
- New to write: the commit attribution table (`42ac030`, `58b32a9`, `8609f12`,
  `f6ee86d`, `f70f86d`, `325314f`, `06dac83`); what
  `implementations/original/` is and is not.

**`ch04`. Three defects of the published process.**

- Scope: matrices above the compound-symmetry ceiling, silently repaired; step
  coupling and the collapse under negligible carryover (attenuation, zero
  coupling gap); the compounding recursion (mean shortfall; negligible effect
  on power at 0.1 to 0.2 week).
- Sources: 36 Sections 1, 4, 5 (spine); 55 Sections 1, 5 to 8 (summary only);
  53 Sections 1, 8 (summary only); 37 Section 5.1; 34 Section 4.6; 50 Section
  1; 07 counts via Appendix I.
- New to write: one table reconciling the published power at $t_{1/2} = 0.1$
  week (0.126, 0.078, 0.100, 0.076: cell, code, seeds), from 55 Section 7.3;
  one reference cell quoted throughout the series.


### 3.3 Part II. The data-generating process

| ID | Title | Needs | pp. |
|----|--------------------------------------------------|--------------|----|
| `ch05` | Response components and the modified Gompertz trajectory | `ch02` | 16 |
| `ch06` | Carryover on the mean | `ch05` | 10 |
| `ch07` | Interaction architectures: mean, covariance and dual-channel | `ch06` | 20 |
| `ch08` | Biomarker coupling under carryover | `ch07` | 12 |
| `ch09` | The response covariance and the coupling ceiling | `ch08` | 30 |
| `ch10` | The refined process: specification, the adjustments evaluated, and the roadmap | `ch04`, `ch09` | 18 |
| `ch11` | Features of an N-of-1 data-generating process | `ch10` | 26 |
| `ch12` | Component decomposition, identification and reduction | `ch05`, `ch07` | 22 |
| `ch13` | Informative dropout and missingness | `ch02`, `ch10` | 10 |
| `ch14` | Beyond the refined process: a pharmacodynamic sensitivity model | `ch09`, `ch10` | 10 |

Table: Part II of the series

**`ch05`. Response components and the modified Gompertz trajectory.**

- Scope: $BR$, $PB$, $TV$, $\mathrm{BL}_i$; their clocks; $G(t_{od}; m, d,
  r)$, landmarks, conditioning; alternative families; paper 07's family
  result.
- Sources: 24 Components 1 to 3; 25 Sections 2 to 3, 5 to 8 (literature,
  alternatives, ranked impact); 40 Sections 3 to 4; CM "The response
  trajectory" (4 subsections); paper 07; 18 via Appendix I.
- New to write: a family-overlay figure; 25 Section 6 and Recommendation 4
  rewritten to the corrected account of Section 3; the link from 25
  Recommendation 5 to P04's $v_i\Psi_{it}$.

**`ch06`. Carryover on the mean.**

- Scope: decay families; the half-life as a pharmacodynamic response
  half-life; anchored against compounding carryover; restart on re-exposure;
  the 2024 scale factor.
- Sources: CM "Carryover on the mean" (3 subsections); 53 Section 5; 51
  Features 3 to 4; 38 Sections 2 to 3; 40 Section 3; 24 Component 1.
- New to write: the reconciliation of a 2.3-hour plasma half-life with a 0.5
  to 1 week carryover half-life.

**`ch07`. Interaction architectures: mean, covariance and dual-channel.**

- Scope: definitions; estimand and units per architecture; the latent-variable
  reading; biological rationale; paper 01's contrast; paper 11's dual channel.
- Sources: paper 01; CM "Two architectures", "The dual-channel architecture",
  "A latent-variable reading"; 21 "Biological Assumptions"; 23 Sections 1.6,
  8.1 to 8.2, 11; 02 Section 4; 33 Sections 2.3 and 3.2; 36 configurations D
  and E; `plan-architecture-c.md`; Appendix I Section 5 (calibration).
- New to write: an estimand-and-units table; the current contrast (NOTATION
  Part 2; paper 01 under graded coupling) replacing 10 to 15% against 40 to
  60%; the mean-moderation estimand shift (33 Section 2.3).

**`ch08`. Biomarker coupling under carryover.**

- Scope: coupling profile $g_t$; step, graded and ratio rules; the
  pharmacokinetic argument and its two qualifications; the coupling gap
  $\mathcal G$; shape of the decay.
- Sources: CM "The estimand and the coupling rules", "The pharmacokinetic
  argument for graded coupling", "ratio-graded", "Shape of the coupling
  decay"; 55 Sections 4 to 5; 37 Section 4.2; 51 Feature 5; Appendix I Section
  7.
- New to write: the answer to "is graded decay structurally required" (02 page
  5, `notes.md`); the trait-against-state variance question as an open
  modeling choice.

**`ch09`. The response covariance and the coupling ceiling.**

- Scope: joint correlation matrix; the exact ceiling $c_{bm}^{\ast}$; the
  compound-symmetry construct, configuration F and the separable form;
  Kronecker decompositions; factorized ceiling and per-step costs; schedule
  dependence; silent repair; components and the ceiling; prognostic couplings
  in the ceiling; the factored generator; the outcome covariance.
- Sources: 35 "In brief" and Sections 4.1 to 4.8 (spine); 34 Sections 3 to 5;
  36 Sections 4, 6.4, 6.5, 6.7, Appendix A; 37 Section 5.8; 32 Sections 3 to 5
  (numbers labeled configuration F); 54 Section 5 (crossed decomposition, as
  corrected) if CM lacks it; 31 Section 20; 44 Steps 5 to 6; CM "Validity of
  covariance binding" (15 subsections).
- New to write: one full ceiling table (construct by $\rho \in \{0.7, 0.8\}$
  by design by $t_{1/2}$); a box "a higher ceiling is not more power" (35
  Section 6).

**`ch10`. The refined process: specification, the adjustments evaluated, and
  the roadmap.**

- Scope: the authoritative definition; the one table of differences from the
  published process; adjustments 1 to 3 and separability over configurations A
  to F; power under the adjustments; the step table of 44; a faithful
  published reference arm; operating characteristics of the refined process as
  a whole.
- Sources: 36 Sections 6.1 to 6.6, 6.8; 44 Sections 1, 3, 5, 6; 32 Section 2
  (repair, choice, consequence); NOTATION "Named data-generating processes";
  paper 01 `tab:spec`; 50 Section 1; 51 Section 20; 23 Section 2 (history
  only).
- New to write: the difference table; the refined-process operating
  characteristics from job 14 and job 15; the unreported common-grid contrast
  behind `36-fig4-contrast-reps100.png`.

**`ch11`. Features of an N-of-1 data-generating process.**

- Scope: the seventeen-feature framework and the five-DGP comparison.
- Sources: 51 Sections 1 to 19, 22 (spine); CM "Further features" (8
  subsections).
- New to write: the evidence base stated per feature (the stress-test driver,
  now under the refined process).

**`ch12`. Component decomposition, identification and reduction.**

- Scope: identification of $BR$, $PB$, $TV$ by the Hybrid design; three
  contrasts; the $\eta$ limit; rank of the phase-augmented model; paper 06's
  studies; the two-component reductions.
- Sources: 24 identification and "Why the analysis does not match"; 31 Parts
  II, III (paper 06 studies, totals), IV, V; 28 Sections 3 to 4 (2026-10-08
  corrections); CM "Partitioning the treatment effect" (4 subsections); papers
  06, 13, 14.
- New to write: 31 Sections 18.1, 18.3 and 22.2 brought up to date (the
  prognostic couplings of the stress-test driver exist; F1 carries `De`).

**`ch13`. Informative dropout and missingness.**

- Scope: the dropout hazard; MCAR, MAR, MNAR; interaction with design; the
  stress test's MAR and MCAR mechanisms; the published censoring grid.
- Sources: CM "Informative dropout" (4 subsections); paper 09 (with the
  matched-$N$ correction); 46 Section 5; 47 dropout cells; Appendix I Section
  2.
- New to write: most of the chapter (no whitepaper covers dropout).

**`ch14`. Beyond the refined process: a pharmacodynamic sensitivity model.**

- Scope: turnover model; sensitivity binding; consequences for the test; a
  proposal, not implemented.
- Sources: 38 (spine, renamed $f(t) \to h(t)$, $F \to \Delta_h$); 44 Steps 7
  to 8; CM "A pharmacodynamic sensitivity model" (4 subsections); 36 Section
  6.7.
- New to write: $\sigma_v^2$ rewritten as $\tau^2$ (DQ-6, decided).


### 3.4 Part III. Theory of the interaction test

| ID | Title | Needs | pp. |
|----|--------------------------------------------------|--------------|----|
| `ch15` | Closed-form theory: E9, the coupling gap and power | `ch08`, `ch09` | 30 |
| `ch16` | Design efficiency and the choice of schedule | `ch15` | 20 |

Table: Part III of the series

**`ch15`. Closed-form theory: E9, the coupling gap and power.**

- Scope: E9 and E9b; Propositions 1 to 3; power through the coupling gap; step
  coupling gives zero gap; graded coupling declines monotonically; anchored
  against compounding means; covariance binding as mean moderation; other
  designs; AR(1) structures; validation; the baseline row.
- Sources: 37 Sections 3 to 6 (spine); 41 Proposition 4 and Section 6; 39
  Section 2; CM Part "Closed-form theory" (13 subsections); 17 via Appendix I
  (origin of E9).
- New to write: the reconciled anchored-change figures (CR-31).

**`ch16`. Design efficiency and the choice of schedule.**

- Scope: the criterion; contrast variance by design; switches, spacing,
  open-label visits, BDC and CO blocks; the exhaustive search; the design
  procedure; cycle count and period length.
- Sources: 39 (spine); 34 Section 5 (feasibility levers); 27 Section 4
  (information formula, corrected optimum (2, 6)); CM "Design efficiency and
  time adjustment" (16 subsections); papers 05 and 12.
- New to write: the criterion renamed (39's $r_{\Delta B}^2$ is not a squared
  correlation); the construct stated in the opening.


### 3.5 Part IV. Analysis of the interaction

| ID | Title | Needs | pp. |
|----|--------------------------------------------------|--------------|----|
| `ch17` | The analysis model: exposure codings, time and baseline | `ch15` | 22 |
| `ch18` | Inference and calibration | `ch17` | 16 |
| `ch19` | The prognostic biomarker and the shared biomarker effect | `ch17`, `ch18` | 16 |
| `ch20` | The covariance study: size of candidate analyses | `ch18`, `ch19` | 16 |
| `ch21` | The stress test: protocol, results and the refined-process rerun | `ch19`, `ch20` | 26 |
| `ch22` | The recommended analysis | `ch21` | 22 |
| `ch23` | Choosing the analysis half-life from the data | `ch22` | 30 |

Table: Part IV of the series

**`ch17`. The analysis model: exposure codings, time and baseline.**

- Scope: the common model; E1 to E9 and the extended E series; the open-label
  branch; the pre-exposure value of $D_{bc,it}$; time adjustment (linear,
  spline, visit factor, phase); baseline as a row or a covariate (Propositions
  1 to 5); cLDA; the published analysis term by term; why the analysis does
  not mirror the DGP.
- Sources: CM "The analysis model and its carryover specifications" (6
  subsections); 43 Section 2; 40 Sections 8 to 9; 41 Propositions 1 to 5 and
  Section 7; 42 (pros, cons, update); 45 Section 7; 44 Section 5; 37 Appendix
  A "Visits before first exposure"; 24 "Why the analysis model does not
  match"; `plan-coherent-analysis-comparison.md` Q1 to Q3.
- New to write: one consolidated time-adjustment section (40 Section 9 and 42
  merged).

**`ch18`. Inference and calibration.**

- Scope: what a random intercept absorbs under AR(1); a decomposition of
  realized Type I error; $\kappa$ by design and architecture; Kenward-Roger
  and CR2 repair different failures; Satterthwaite degrees of freedom; the
  cost of dichotomizing.
- Sources: CM "Inference" (6 subsections); 45 Section 4; 33 Section 2.4; 56
  Section 3.5; 36 Section 6.8; Appendix I Section 4 (the `lmer` evidence);
  Appendix H (RM-ANOVA).
- New to write: a $\kappa$ section (33 Section 2.4 has the only analysis); the
  0.14 to 0.16 against 0.14 to 0.23 ranges labeled by construct (CR-28).

**`ch19`. The prognostic biomarker and the shared biomarker effect.**

- Scope: the shared slope; the leak $(1-f)s_{\text{prog}}$; freeing the
  baseline slope; the expectancy terms; the pathways table; leak factors; the
  power cost.
- Sources: 48 Sections 1 to 3 (lead-in); 49 Sections 2 to 4; 56 Section 3.3
  (rigorous core); 45 Section 6; 41 Section 3; 51 Feature 6; CM "Baseline, a
  prognostic biomarker, and the shared biomarker effect" (6 subsections),
  "Prognostic pathways".
- New to write: F1 stated as approximately unbiased (CR-21); leak factors
  under the refined process if job 15 supports them.

**`ch20`. The covariance study: size of candidate analyses.**

- Scope: 63 analyses, 16 cells, 378,000 fits; size table; Kenward-Roger
  against CR2; size-adjusted power; the shared $\beta_B$; $D_{bc,it}$ as
  insurance; the manuscript 02 coding defect.
- Sources: 45 Sections 2 to 8 (spine); 43 Sections 3.7 to 3.8; the anchored
  main pass (`summary-alt250-null500-anchored.csv`, fixed allocation).
- New to write: the anchored-pass results; a statement of what the cancelled
  follow-up and $D_{bc,it}$ passes would have supplied.

**`ch21`. The stress test: protocol, results and the refined-process rerun.**

- Scope: criteria; the 14 analyses; the 41 cells; results of the stored run;
  results under the refined process (job 15, 11 analyses); G1, G2 and F1A
  paired with job 15 (jobs 21, 22); the symmetry screen; failed fits; the
  protocol of record as an annex.
- Sources: 47 (spine); 46 (annex, with Amendment 2 item 4 cut to a pointer);
  48 Section 4; 49 Section 5; 50 Section 2; CM "The stress-test construct and
  its fourteen analyses" (5 subsections), "Size, power screens and the
  selection criteria".
- New to write: the job 15 results section; the job 21 and 22 results; the
  comparison of the two runs as independent proportions (`RERUN-anchored.md`
  Section 7).

**`ch22`. The recommended analysis.**

- Scope: F1 term by term; the fixed half-life of 0.2 week; worked examples;
  the decision framework; alternatives; the recipe in R; limitations; a
  summary of F1A0B with a pointer to `ch23`.
- Sources: 56 Sections 1 to 5, 7 to 8 (spine); 49 Section 6; 57 Sections 1 and
  7.2 (summary); 57 Section 8.2 (F1h02 evidence).
- New to write: the status, half-life and limitations sections rewritten;
  Section 5 of 56 rebased on job 15 and job 14.

**`ch23`. Choosing the analysis half-life from the data.**

- Scope: why AIC with the interaction inflates size; what the data say; eight
  methods and permutation; worked examples; the decision table; the extension
  results; the cost.
- Sources: 57 (spine; Sections 6 and 7.3 moved to an annex); 47 Section 8; 56
  Section 6; CM "AIC selection of the analysis half-life", "The assumed
  half-life of $D_{bc}$"; Appendix I Section 9 (`characterize_carryover`).
- New to write: F1A0, the bootstrap, max-$\lvert t\rvert$ and the
  orthogonality result moved into CM; the full-bootstrap check (GP-15) when
  run.


### 3.6 Part V. The treatment main effect and application

| ID | Title | Needs | pp. |
|----|--------------------------------------------------|--------------|----|
| `ch24` | Testing the drug main effect when the half-life is unknown | `ch15`, `ch23` | 30 |
| `ch25` | Analyzing a real trial | `ch22` to `ch24` | 12 |

Table: Part V of the series

**`ch24`. Testing the drug main effect when the half-life is unknown.**

- Scope: estimands $\beta_D^{\text{full}}$, $\beta_D^{\text{exp}}$,
  $\beta_D^{\text{bl}}$; probability limits and null bias; row-group weights;
  anchored against compounding; closed-form power; the half-life not
  identified under the null; max-$\lvert t\rvert$; randomization tests;
  averaging; joint estimation.
- Sources: 58 (spine); P04 (pointer for operations); CM "The treatment
  main-effect track" (to be rewritten); paper 04.
- New to write: the `M-` labels; the fixed 0.2 week policy stated for both
  tracks with 57's reconciliation of size (57 Section 1, Section 7.2 item 7).

**`ch25`. Analyzing a real trial.**

- Scope: `analyze_trial_extended()`, `characterize_carryover()`; predicted
  effect and threshold on the fitted (centered) scale; which visits carry the
  carryover signal; reporting (CENT 2015, PATH); two-stage random slopes for
  real data.
- Sources: 14 via Appendix I Section 9; vignette 03; 56 Section 8; 57 Section
  5; 58 Section 10; 23 Section 10; 02 page 13 comment.
- New to write: most of the chapter.


### 3.7 Part VI. Notes for the original authors

| ID | Title | Needs | pp. |
|----|--------------------------------------------------|--------------|----|
| `ch26` | Note: the carryover recursion in the published code | none (standalone) | 20 |
| `ch27` | Note: step coupling in the published code | none (standalone) | 24 |
| `ch28` | Note: patching `58b32a9` | `ch26`, `ch27` | 8 |

Table: Part VI of the series

**`ch26`. Note: the carryover recursion in the published code.**

- Scope: 53 as a standalone note.
- Sources: 53 (whole).
- New to write: deduplicate Section 8.3 and Section 9.5 against `ch28`.

**`ch27`. Note: step coupling in the published code.**

- Scope: 55 as a standalone note, plus the upstream never-treated `Dbc` coding
  defect.
- Sources: 55 (whole); 54 Section 3.6.2.
- New to write: the 54 Section 3.6.2 item.

**`ch28`. Note: patching `58b32a9`.**

- Scope: the eight-line graded patch (power at 0.1 week from 0.126 to 0.638),
  the anchored fix, and their joint effect, with one power table.
- Sources: 37 Appendix A; 53 Section 9; 55 Sections 9 to 10 (Table 8, runs A
  to C); 36 patch material; the patch README under `analysis/scripts/`.
- New to write: the single power table.


### 3.8 Appendices and Volume II

| ID | Title | Sources | New to write | pp. |
|---|--------------|------------------------|----------------|---|
| A | Catalog of analyses and the ID registry | 43 (regenerated from `analyses.csv`, extended to 56 to 58) | the registry; the global IDs of Section 2.7 | 20 |
| B | Decision history | the dated notes of 35, 41 to 45, 49, 50, 56, 57; 52 Section 2; 30 execution record; 45 Section 9; 57 superseded plans | all | 8 |
| C | Errata to earlier documents and manuscripts | 36 Section 8; 40 Section 11; 45 Sections 8 and 10; CM main-effect section | the status of each item | 8 |
| D | Implementation status and open-items register | 35 Sections 7, 9; 44 Section 4; 52 Section 4; 23 Section 14; 34 Sections 4.7, 6.2; 36 Section 9 | one maintained register with one inspection date | 8 |
| E | Lexicon | 29 (with the $\tau^2$ entries fixed) | a project-terms table pointing to `fm2` | 16 |
| F | Style guide | 30 (extended to the whitepapers); Sections 2 and 7 of this plan | the whitepaper scope | 6 |
| G | Bibliography | 23 Section 17 (promoted); `docs/references/*.bib` | verification of the suspect entries (Appendix I Section 11) | 14 |
| H | Classical test equivalents: repeated-measures ANOVA | 26 (whole; reserved letters renamed) | none beyond renaming | 14 |
| I | The legacy record: results of docs/00 to docs/20 not carried elsewhere | see the outline below | the consolidation (decided) | 40 |
| J | Software and parity workflow (pointer only) | 22 (parity), moved to `analysis/scripts/parity/README.md` | a one-page pointer | 1 |
| Vol. II | Mathematics compendium | CM | updates of GP-28 | (existing) |

Table: Appendices of the series

**Outline of Appendix I, the legacy record (decided).** One document
that holds only what docs/00 to docs/20 (with 02, 09a and 19) contain
and no newer whitepaper covers. Each section opens with a provenance
line (source document, section, date, construct) and a status label;
none restates a superseded mechanism except to name it as superseded.

| Section | Content | Sources | Notes |
|----|----------------------------------------|----------------------|------------------------------|
| 0 | Purpose, scope, and a map from each original to where its superseded content now lives | survey a recommendations | lists what was deliberately dropped (07a entire; 00, 01 entire; 03 and 04 mechanisms) |
| 1 | Commit history of the published code and the 2024 revision: per-commit agreement with published Figure 4 ($r$), PD failure rates (40 of 162 at publication, 102 of 162 after `8609f12`), condition numbers (median 229 against $3.1 \times 10^{18}$), eigenvalue range | 04 Appendix A.7; 07 Section 3; 23 Section 3.1 | commit attribution corrected to `58b32a9` against `42ac030` and the 2024 commits |
| 2 | The published simulation grid: the four `buildtrialdesign` calls, the censoring model and its five patterns, 72 cells by 1,000 replicates by 5 censoring conditions | 13 Sections 2, 4, 6 | each entry checked against `58b32a9`; "last observation carried forward" corrected to monotone dropout; dropout form checked against `censordata` (17 Section 7.1 item 4 disagrees) |
| 3 | Reproductions of published Figures 4 and 5 | 04 Section 2 ($r = 0.994$ at 200 replicates); 18 Sections 1 to 5 | 18's figures are missing; its comparison is reported as a historical result under the April construct, or regenerated (GP-26) |
| 4 | The analysis-model switch: `lmer` Type I error under AR(1) data (200 replicates; OL 0.13/0.09, CO 0.17/0.13 at two sizes) and `corCAR1` for unequal spacing | 06 Sections 1, 3.3, 6 | 06's residual-correlation formula replaced by CM `eq:inva-ri`; "always understates" qualified |
| 5 | Mean-moderation formulations that failed (three attempts) and the calibration check (3.36 and 3.48 per SD against 3.6); positive definiteness under covariance moderation against mean moderation; the published grid under mean moderation | 19.md (reviewed); 19.tex Sections 6 to 7 | 19.tex ceilings ("about 0.45") replaced by `ch09` values; the analysis variable of each table stated |
| 6 | The constant-exposure study and the two-period crossover bias | 09a.md Sections 2 to 7 (corrected 2026-10-08) | the Shi and Ye (2024) expression with corrected signs; the withdrawn conclusion 4 stated as withdrawn |
| 7 | Open modeling questions: trait against state variance off drug; clinical plausibility of a negative $BR$ under the step rule | 08 Section 8; 09 Section 5 | carried to `ch08` as open items |
| 8 | The co-author's ten questions on the April architecture paper, and where the series answers each | 02 annotations (pages 2, 3, 5, 7, 13); `notes.md` | a table: question, answer, chapter |
| 9 | The March 2026 collaborator goals and their status: Goal 2b (residual analysis, never implemented); the threshold and predicted-effect formulas (defective, fixed in code 2026-10-08); `characterize_carryover()` grid (0.1 or 0.25 to 8 weeks, to be checked) | 14; 03 Section 6 | links to `ch23`, `ch25` |
| 10 | Abandoned designs and software history: the simplification plan (not adopted) as the origin of E9 and of the RM-ANOVA comparator; the 2025 tidyverse alignment and its unused decay functions | 17 Sections 3.3, 8; 15; 16 item 1 | one paragraph each |
| 11 | Bibliography corrections: the misattributions of 20 and their corrections; checks of Chen, Li and Yang (2020, *PLoS One* e0216949) and Eng and Konings (2018, *Trials* 19:72); also Pequignat et al. (2023) and De Bock et al. (2021) | 20; 17 Section 7; 23 Section 17; 02 references | each entry verified against Crossref or PubMed or dropped (GP-29) |

Table: Outline of Appendix I, the legacy record (decided 2026-10-11)

## 4. Disposition of every current file

"Archive" means `docs/legacy/` with the file name unchanged and a
one-line banner naming the replacing chapter (Section 9.3). "Remove from
series" means the file moves out of `docs/` (destination named) and is
not part of the read-through. Where the surveys disagreed, the row says
so and gives the resolution.

| File(s) | Destination | Reason |
|----------------------------------|------------------------------------|------------------------------|
|  00-documentation-index.tex | archive | superseded; replaced by the generated index `fm0`  |
|  01-codebase-overview.tex | archive | repository description belongs in `CLAUDE.md` and package docs; nothing unique  |
|  02-dgp-mean-moderation-vs-mvn-rchcoments.pdf | `docs/reviews/`, renamed review-comments-2026-04-29.pdf; Appendix I Section 8 | record of the co-author's comments; the questions enter Appendix I  |
|  03-audit-and-revision-report.tex | archive; Appendix I Sections 1, 9 | claims superseded; corrected commit attribution in `ch03`  |
|  04-revised-power-analysis.tex | archive; Appendix I Sections 1, 3 | commit history and Figure 4 reproduction are unique; mechanism and ceilings superseded  |
|  06-ar1-residual-correlation.tex | archive; Appendix I Section 4 | `lmer` Type I evidence unique; derivation superseded by CM  |
|  07-positive-definiteness-failures.tex | archive; Appendix I Section 1 | counts and condition numbers kept; root cause refuted by `docs/36` Section 8  |
|  07a-pd-constraints-theory.tex, .pdf | archive | superseded and partly wrong; nothing carried (surveys a and c agree)  |
|  08-biomarker-correlation-decay.tex | archive; `ch08` credit; Appendix I Section 7 | origin of graded coupling; trait-against-state question unique  |
|  09-carryover-correlation-artifact.tex | archive; Appendix I Section 7 | superseded by 55 and 36; clinical-plausibility point kept  |
|  09a-carryover-bias-variance.md, .pdf, .tex | archive all three; Appendix I Section 6 | survey a offered "keep the `.md` as an appendix"; the user decision places its unique content in Appendix I; the `.tex` contradicts the corrected `.md`  |
|  13-figure4-walkthrough.tex | archive; Appendix I Section 2 | grid and censoring specification unique; DGP sections mislabel the 2024 code  |
|  14-trial-analysis-implementation-guide.tex | archive; Appendix I Section 9; `ch25` | goals record; formulas superseded by the 2026-10-08 code fix  |
|  15-tidyverse-alignment-plan.tex | archive; Appendix I Section 10 | planning record  |
|  16-tidyverse-code-review.tex | archive; Appendix I Section 10 | planning record  |
|  17-simplification-plan.tex | archive; Appendix I Sections 10, 11 | origin credit for E9 and RM-ANOVA; a suspect reference  |
|  18-response-parameter-sensitivity.tex | archive; Appendix I Section 3 | figures missing; result kept as historical  |
|  19-mean-moderation-implementation-notes.md, .pdf, .tex | archive all three; Appendix I Section 5; `ch07` pointer | survey a recommended a merged appendix to the architectures chapter; the user decision places it in Appendix I  |
|  20-annotated-bibliography.tex | archive; Appendix G; Appendix I Section 11 | misattributions; corrected list in 23 Section 17  |
|  21-architecture-comparison-summary.md, .pdf | archive; `ch07` (rationale) | superseded; 50-replicate tables unreproducible  |
|  22-journal-target-recommendation.md, .pdf | remove from series: the paper 01 folder, `planning/` | planning, not methods  |
|  22-tidyverse-development-parity-workflow.md, .pdf | remove from series: a README in `analysis/scripts/parity/`; Appendix J pointer | current developer guide, not methods; frees the duplicate number  |
|  23-documentation-compendium.md, .pdf | archive; salvage 3.1 to `ch03`, 8.1 to 8.2 to `ch07`, 14 to Appendix D, 17 to Appendix G | every section has a fuller current home  |
|  24-component-decomposition-pedagogy.md, .pdf | `ch05` (components) and `ch12` (identification) | current pedagogy; analysis advice rewritten to point to `ch22`  |
|  25-gompertz-model-evaluation.md, .pdf | `ch05` | literature, alternatives and ranked impact kept; numerics from CM  |
|  26-rm-anova-interaction-test.md, .pdf | Appendix H | current derivation, reference character  |
|  27-cycle-period-design-sweep-plan.md, .pdf | remove from series: the paper 05 folder, `planning/`; Section 4 to `ch16` | planning record; corrected information result kept  |
|  28-hendrickson-hybrid-design-and-decomposition.md, .pdf | archive after porting its 2026-10-08 corrections to `ch12` | superseded by 31 (its own banner)  |
|  29-nof1-precision-medicine-lexicon.md, .pdf | Appendix E | reference; fix $\tau^2$  |
|  29-nof1-precision-medicine-lexicon_2026-06-18_1719.pdf | delete | stale dated duplicate  |
|  30-terminology-consistency-plan.md, .pdf | Appendix F (extended); execution record to Appendix B | prose authority beside NOTATION  |
|  30-terminology-consistency-plan_2026-06-18_0834.pdf | delete | stale dated duplicate  |
|  31-component-decomposition-and-reductions.md, .pdf | split: Parts II to V to `ch12`; Parts I and VI removed (pointer to `ch05`); then archive | executes `plan-doc31-revision` Phase 0  |
|  32-why-covar-departs-from-orig.md, .pdf | `ch10` (framing), `ch09` (numbers labeled configuration F); then archive | survey b recommended merging into 54; 54 is itself largely superseded (survey d), so the framing goes to `ch10` instead  |
|  33-paper01-dgp-architecture-review.md, .pdf | remove from series: the paper 01 folder, `revisions/`; findings 2.3, 2.4 to `ch07`, `ch18`; Section 2.5 to `ch02` | manuscript review, not a methods chapter  |
|  34-cbm-feasible-range.md, .pdf | `ch09` (formula, schedule dependence), `ch16` (Section 5 levers); Sections 6.1, 7 archive; 6.2 to Appendix D | partly superseded at $\rho_{cs} = 0.7$  |
|  35-separable-response-covariance.md, .pdf | `ch09` spine; Sections 7, 9 to Appendix D | current  |
|  36-hendrickson-pd-and-carryover.md, .pdf | split: Sections 1 to 5 to `ch04`; 6.1 to 6.6, 6.8 to `ch10`; 6.7 and Appendix A.7 to `ch09`; 8 to Appendix C; 9 to Appendix D | current; drop "Working draft"  |
|  37-carryover-power-closed-form.md, .pdf | `ch15` spine; Appendix A to `ch28` | current  |
|  38-pharmacodynamic-carryover-model.md, .pdf | `ch14` spine | proposal; symbols renamed to 44's  |
|  39-nof1-design-optimization.md, .pdf | `ch16` spine | current within its construct  |
|  40-hendrickson-hybrid-design-primer.md, .pdf | split: Sections 2 to 3 to `ch02`; 3 to 6 and 8 to `ch03`, `ch05`; 9 to `ch17`; 11 to Appendix C | primer; duplicated by 42  |
|  41-baseline-as-covariate.md, .pdf | `ch17` (Propositions 1 to 5); Proposition 4 to `ch15`; recommendation to Appendix B | partly superseded  |
|  42-flexible-time-adjustment.md, .pdf | `ch17` (time adjustment, merged with 40 Section 9) | duplicate tables removed  |
|  43-analysis-approaches-catalog.md, .pdf | Appendix A (regenerated) | reference, incomplete since 56  |
|  44-dgp-roadmap.md, .pdf | `ch10` (step table); Steps 7 to 8 to `ch14`; status to Appendix D | planning, current as status board  |
|  45-correctly-sized-hybrid-analysis.md, .pdf | `ch20` spine; Section 6 to `ch19`; 7 to `ch17`; 8 to Appendix C; 9 to Appendix B | covariance-study record  |
|  46-strawman-stress-test-protocol.md, .pdf | `ch21` annex (protocol of record) | historical protocol; Amendment 2 item 4 cut to a pointer  |
|  47-stress-test-results.md, .pdf | `ch21` spine | current record of the stored run  |
|  48-shared-biomarker-effect-and-stress-test-issues.md, .pdf | Sections 1 to 3 to `ch19`; 4 to `ch21`; 5 dropped | partly superseded; Section 4 table corrected (CR-20)  |
|  49-predictive-prognostic-biomarker-expectancy-term.md, .pdf | Sections 2 to 4 to `ch19`; 6 to `ch22`; 5 and 7 to Appendix B | recommendation superseded  |
|  50-dgp-and-analysis-status-brief.md, .pdf | archive; regenerated as `fm1` | superseded brief  |
|  51-ten-features-of-an-nof1-dgp.md, .pdf | `ch11` (renamed "seventeen features"); Section 20 to `ch10`; 21 to paper 01 | current framework  |
|  52-whitepaper-review-2026-10-07.md, .pdf | archive; Section 2 to Appendix B; Section 4 to Appendix D | review log  |
|  53-carryover-recursion-in-published-code.md, .pdf | `ch26` (kept standalone, DQ-3) | current note for the original authors  |
|  54-covar-vs-orig-dgp.md, .pdf | archive; Section 3.6.2 to `ch27`; Section 5 to `ch09` if CM lacks it | largely superseded; assessment of the ratio rule overtaken by 55  |
|  55-step-coupling-in-published-code.md, .pdf | `ch27` (kept standalone, DQ-3) | current note for the original authors  |
|  56-recommended-analysis-expectancy-cr2.md, .pdf | `ch22` spine; Section 6 to `ch23` | current, internally stale in Section 9  |
|  57-data-chosen-carryover-half-life.md, .pdf | `ch23` spine; Sections 6 and 7.3 to an annex | current  |
|  58-main-effect-testing-and-carryover-half-life.md, .pdf | `ch24` spine | current  |
|  59-whitepaper-refactoring-plan.md, .pdf (this plan) | archive at release; its Section 4 becomes `fm0`'s old-to-new map | planning  |
|  plan-architecture-c.md | remove from series: the paper 11 folder, `planning/` | superseded as to placement  |
|  plan-coherent-analysis-comparison.md | remove from series: the paper 02 folder, `revisions/`; Q1 to Q3 to `ch17` | planning; the forbidden word on line 20 is left in the moved copy only  |
|  plan-doc31-revision.md, .pdf | archive | executed by this refactor  |
|  notes.md | `docs/reviews/` with a header naming its source and date; Appendix I Section 8 | reviewer comments  |
|  reorg-plan-2026-04-22.md | archive | executed in a different form  |
|  ZZCOLLAB_USER_GUIDE.md | remove from series: `docs/tooling/` | framework guide, no project content  |
|  pub_review_whitepaper_2026-08-16.md, pub_review_whitepaper_2026-08-20.md, pub_review_remediation_2026-08-20.md | `docs/reviews/` | manuscript reviews of the August campaign  |
|  Architecting_Precision_Medicine.pdf | moved to `docs/archive/`, labeled unreviewed (DQ-4, decided) | 12.2 MB AI-generated deck, no text layer, cited nowhere, unreviewed  |
|  pmsimstats-preamble.tex | archive, beside the archived `.tex` set | needed to compile the archived `.tex`  |
|  figures/ (44 files) | stays; new figures named `chNN-figK-*.png`; `40-fig2` and `40-fig3` renamed to print order | archived documents keep their paths  |
|  archive/ (4 files) | stays, for non-whitepaper historical material; receives the slide deck (DQ-4) | already historical; the old whitepapers go to `docs/legacy/` instead |
|  references/ (13 files) | stays; `nof1-simulation.bib` becomes the series bibliography source | `pmsimstats-nof1-whitepaper.md` there was not surveyed (Section 11)  |
|  reviews/ (1 file) | stays; receives 02's PDF, `notes.md`, the `pub_review` files | review records  |
|  .claude/ | unchanged | settings, not documentation  |

Table: Disposition of every file in `docs/` (survey recommendations;
disagreements resolved in the rows for 09a, 19 and 32)

## 5. Conflict register

Each conflict is listed once. "Best answer" names the source of the
current answer; "Action" names where the series states it.

(Updated 2026-10-11, Phase 1: CR-31, CR-43, CR-44 and CR-51 are
resolved in Section 5.2 by inspection and by comparing stored outputs,
and GP-29 in Section 8. The evidence, with commands and outputs, is in
the Phase 1 log `docs/series/notes/phase1/phase1-cr-log.md`, with its
check scripts; it is to be carried into Appendix I. Other `SP/` paths
below name scratch outputs that were not kept in the repository. Corrections for archived
documents (13, 17, 20, 23 and 03 in `docs/legacy/`) are recorded there,
not made in the archived files.)

### 5.1 Resolved by later work: needs propagation only

| ID | Conflict | Documents and sections | Best answer (source) | Action |
|-------|----------------|---------------------------|------------------------------|--------------|
| CR-01 | Mechanism of the published collapse: standard-error inflation with an unbiased estimate, or attenuation | 03 Section 4.2; 04 Section 3.3; 09 Sections 4.4, 6.3; 02 Section 6 against 36 Section 5.3, 37 Section 5.1, 55 Section 5 | zero coupling gap; the estimate attenuates by about three quarters (36 Section 5.3; 55 Section 5) | `ch04`; old documents archived |
| CR-02 | What the step rule does off drug | 03 Section 2.2 ("correlation drops to zero instantly", lines 57 to 58, spot-checked) against 04 Appendix, 08 Section 1, 09 Section 3.1, 55 Section 4 | full $c_{bm}$ wherever the off-drug mean is nonzero (55 Section 4, code quoted) | `ch03`, `ch04` |
| CR-03 | Cause and pattern of the published PD failures | 03 Section 3.2; 04 Section 4.2; 07 Sections 3.4, 4.1 against 36 Sections 1, 4, 8 | $c_{bm}$ above the compound-symmetry ceiling 0.256 to 0.262; diagonal scaling cannot change definiteness (36 Section 8) | `ch04`; Appendix I Section 1 keeps the counts only |
| CR-04 | Coupling ceilings | 01 Section 6; 03 Sections 2.1, 2.7; 04 Table A.5; 07a Section 8; 19.tex Section 7.1; `CLAUDE.md` against 34, 35, 36, 52 Section 2 | by construct: 0.256 (compound symmetry, $\rho_{cs} = 0.8$); 0.440 (separable, $\rho = 0.8$); 0.474, Hybrid 0.480 (separable, $\rho = 0.7$) (52 Section 2) | `ch09` full table; `CLAUDE.md` (Section 9.5) |
| CR-05 | Size of the architecture contrast | 02 Sections 3.1, 5.1; 19.md Background; 21; 23 Sections 1.4 to 1.6, 18; `CLAUDE.md` against NOTATION Part 2 and paper 01 | Hybrid 6.9% (mean) against 15.4% (covariance) at matched total $N = 70$; under graded coupling losses nearly equal (NOTATION Part 2; paper 01 as restated) | `ch07`; `CLAUDE.md` |
| CR-06 | Commit attribution of the scale factor, the timepoint-1 guard and `Dbc` | 03 Sections 2.4, 2.5; 06 Section 1.1; 13 throughout against 04 Appendix, 36 Section 8, the memory note | 2024 commits `8609f12` and `325314f`; the published analysis used `lmer` with binary $D_{it}$ | `ch03` attribution table |
| CR-07 | Whether ignoring serial correlation always understates standard errors; OLS inflating Type I error | 06 Section 3.1; 17 Method 5 against CM `sec:inva-ri` | direction depends on the contrast; OLS is conservative for the within-participant interaction (CM) | `ch18` |
| CR-08 | Residual correlation under random intercept plus AR(1) | 06 Section 3.2 against CM `eq:inva-ri` | CM (rescaled by $\sigma_\varepsilon^2$) | `ch18`; Appendix I Section 4 |
| CR-09 | Threshold and predicted-effect formulas | 14 Section 4.2 against the code fix of 2026-10-08 (52 Section 4) | computed on the fitted scale (code after 2026-10-08) | `ch25`; Appendix C |
| CR-10 | Recommended interaction analysis | 24 "Four reasons"; 26 Section 10; 28 and 31 Sections 12 to 13; 40 Section 10; 41 Section 8; 45 Section 9; 49 Section 7; 50 Section 2; 56; `CLAUDE.md` against 57 Sections 1, 7.2 | F1A0B primary data-driven candidate, F1AB alternative, F1 at a fixed half-life the fallback (57 Section 1, lines 209 to 223, spot-checked); the fallback at 0.2 week (decided) | `ch22`, `ch23`; Appendix B |
| CR-11 | Hendrickson et al. (2020) subtitle and author order; Haller and Ulm against "Eng and Konings" for one DOI; the other misattributions of 20 | 02 references; 20 Section 1; 17 Section 7 against 09a.md, 23 Section 17 | "Statistical methods and theoretical findings", Hendrickson first (PubMed-checked in 09a.md and 23 Section 17) | Appendix G; Appendix I Section 11 |
| CR-12 | Hybrid crossover blinded or open-label | 31 Sections 3, 5, 7.1 against 24, 28 Section 3.1, NOTATION | blinded, $\eta_t = 0.5$ (NOTATION) | `ch02`, `ch12` |
| CR-13 | Rank of the phase-augmented model | 31 Sections 7.4, 8 against 24, 28 Sections 3.4, 4.2, CM "Rank of the phase-augmented model" | exact aliasing at every $N$; five identified parameters (24, 28, CM) | `ch12` |
| CR-14 | $\mathrm{Cor}(D, \eta)$ on the Hybrid design | 31 Section 11.3 (about +0.6) against 28 Section 4.3 note (0.509 Hybrid, 0.775 recovery design), 56 Section 3 (0.588) | the reproduced values of 28 and 56, with the design and row set stated | `ch12`, `ch19` |
| CR-15 | Paper 06 Study A sample sizes | 28 Section 4.3 (150, 100) against 31 Section 11.3 (300, 200) | totals (31; NOTATION Part 2 migration) | `ch12` |
| CR-16 | Gompertz parameters per participant or population constants | 28 Section 2.1; 31 Section 6.1 against 24, 25 Section 6, 31 Section 9 | population constants with MVN deviations (24, 25) | `ch05` |
| CR-17 | Whether a biomarker-$TV$ coupling exists | 31 Section 18.1 ("no biomarker-TV counterpart anywhere", line 1837, spot-checked) against NOTATION ($c_{bm,TV}$), 46 to 49 | exists in the stress-test driver (`c_tv`), not in the package | `ch12`, `ch19` |
| CR-18 | Restoring the recursion in the package | 32 Section 6 against paper 01 | anchored carryover (paper 01; 32's dated note) | `ch10` |
| CR-19 | Mean moderation "immune to carryover" | 21; 23 Sections 1.4 to 1.6, 18 against 33 Section 2.3 and paper 01 | the estimand grows with carryover (13% Hybrid, 22% OL+BDC; 33 Section 2.3); losses not of the same kind | `ch07` |
| CR-20 | ES against F1 on the non-prognostic cells | 48 Section 4 table, row "Symmetry screen dropped" (line 384, spot-checked: "ES on non-prognostic cells") against 47 Section 7 and 46 Amendment 1 | F1 ahead (0.242 against 0.278 maximin shortfall; 47 Section 7) | `ch21` |
| CR-21 | F1 exactness under prognostic pathways | 49 Section 1 item 3 against 49 Section 3 note, 47 Section 4, 56 Section 3.3 | approximately unbiased; leak factors $-0.051$ to 0.171 (56 Section 3.3) | `ch19` |
| CR-22 | Status statements: anchored reruns "pending", data-chosen half-life "not yet evaluated" | 38 Section 7; 43 header; 44 Section 1, Step 3, Section 4; 45 header; 56 Section 9 (lines 1004, 1022, spot-checked); 57 Section 1 (line 144, spot-checked) | job 15 finished (stress test, refined process); covariance anchored main pass finished; follow-up and $D_{bc,it}$ passes cancelled; job 14 reported in 57 Section 8; jobs 21, 22 running (Section 1.1) | `ch20`, `ch21`, `ch23` |
| CR-23 | Citation of the 14-analysis decision rerun | 35 Section 10 ("Amendment 1") against 46 lines 401 to 447, 45, 52 | Amendment 2 | `ch09` |
| CR-24 | Symbols of the pharmacodynamic model | 38 ($f(t)$, $F$, $\omega_i$, $\zeta$) against 44 index ($\theta_i$, $\kappa$, $z$; $h(t)$, $\Delta_h$) | NOTATION: $h(t)$; 38's current $\omega_i$, $\zeta$, $b_i$ | `ch14` |
| CR-25 | Flexible-time remedy | 40 Section 10 recommendation 2 (`varIdent` plus `corCAR1`, untested) against 42 update, 45 Section 4.4 | CR2 removes the size excess (42 update) | `ch17` |
| CR-26 | MCSE multiples of the flexible-time excess | 40 Section 9 ("about 3 and 3.5") against 42 Section 6 (3.6 and 4.3, at the MCSE at 0.05) | 42 (Section 2.6 convention) | `ch17` |
| CR-27 | The 2024 ratio coupling rule | 54 Sections 3.3.1 to 3.3.2, 6 ("Keep ... the interaction gate replacement", lines 876 to 880, spot-checked) against 55 Sections 4.5, 9.4, 10 | 55: graded at half the stated rate, visit 1 never coupled, E9 power 0.495 against 0.721 at $t_{1/2} = 0$ (55 lines 444, 464, 685, spot-checked) | `ch08`, `ch27`; 54 archived |
| CR-28 | The published analysis on AR(1) data in CO | 43 Summary item 4 (0.14 to 0.16) against 43 O1 row and 36 Section 6.8 (0.14 to 0.23) | both correct for their constructs ($\rho = 0.7$ covariance study; $\rho = 0.8$ in 36) | `ch18` with construct boxes |
| CR-29 | Two maximin readings | 43 (non-prognostic cells: G2 0.232, F1 0.242) against 47 Section 7 and 52 (all 26 alternative cells: G1 0.308, G2 0.312, F1 0.358) | both correct; label the cell set | `ch21` |
| CR-30 | "Published construct" for AR(1) variants | 48 Section 3 note; 49 Section 3 note against NOTATION, 56 Section 3.3 | configuration A only (NOTATION) | `ch19` |
| CR-56 | Published power at $t_{1/2} = 0.1$ week | 37 Appendix A (0.126); 51; 53 (0.078); 55 run C (0.100); paper 01, 50 (0.076, 0.08) | reconciled in 55 Section 7.3 | `ch04` one table |
| CR-32 | Lists of differences from the published process | 23 Section 2 (seven); 32 Section 2 (four); 54 Section 1 (five); NOTATION; paper 01 Appendix B | NOTATION's definitions plus 53, 55 and 36; resolved in substance, never tabulated once | `ch10` one table (GP-02) |
| CR-33 | What "original Hendrickson" means | 28 Section 2; 31 Section 6 (`implementations/original/`); 54 (`orig` = `06dac83`) against 32, NOTATION, the memory note | the published process is `58b32a9`; `implementations/original/` is a corrected data.table reference; 54's `orig` is the vendored 2024 copy | `ch03` |
| CR-34 | E2 dominance | 27 Section 10 against paper 02 grid | paper 02 (27's note) | `ch16` |
| CR-35 | Meaning of $\tau^2$ in the lexicon | 29 (between-patient variance; also "tau" for effect size) against NOTATION decision 1 | NOTATION | Appendix E |
| CR-36 | The main-effect model in CM | CM "The treatment main-effect track" (linear DGP, days, $\rho = 0.3$) against P04 and 58 | P04 and 58 (refined process, weeks, $\rho = 0.7$) | CM update (GP-28); `ch24` |
| CR-37 | Fixed analysis half-life policy | 46, 47, 49, 50, 56, 57 (0.5 week); 56 Section 8 recipe (0.1 to 0.2 week); 57 Section 7.2 item 3 (0.2 no gain); 58 Section 12 and P04 (0.2 week, line 222 of 58, spot-checked); `CLAUDE.md` (0.5 to 1.0) | 0.2 week in both tracks (decided); 0.5 week named the historical base setting | `fm2`, `ch22`, `ch24`; `CLAUDE.md` |
| CR-38 | Hybrid round trip in the cycle-design plan | 27 Section 2c against 24 Component 2, 33 Section 2.5 | 33 Section 2.5 table | `ch02` |
| CR-39 | $TV$ identification in the Hybrid design | 24 Q&A and summary table ("off-drug timepoints with high efficiency") against 24 identification section | weakly identified; $TV$ never observed without $PB$ (24 identification section) | `ch12` |
| CR-40 | The Gompertz rescale-factor account | 25 Section 6 and Recommendation 4 against 25 Section 3 (corrected 2026-10-08) | numerator cancellation, `expm1` remedy (25 Section 3) | `ch05` |
| CR-41 | Whether `lmer` was correctly sized for the published results | 06 Section 7 ("correct under CS") against 18 Section 5.2 ("anti-conservative") | correctly specified for compound-symmetry data (06; inferred consistent with 36 Section 6.8, which finds excess only on AR(1) data); 36 open item 4 (OL+BDC under CS) remains open (CR-46) | `ch18` |

Table: Conflicts resolved by later work (propagation only)

### 5.2 Genuinely open: needs analysis, inspection or a decision

| ID | Conflict | Documents and sections | Current best answer | Action |
|-------|----------------|---------------------------|------------------------------|--------------|
| CR-31 | Anchored-change figures inside 37 | header note (closed form at most 0.017 in CO; Monte Carlo 0.021; lines 16 to 17, spot-checked) against Appendix A final bullet (closed form up to 0.041 in CO; Monte Carlo 0.018; line 1662, spot-checked); 39 header repeats 0.017 | Resolved 2026-10-11: both closed-form figures are correct for different grids. Closed-form change in CO is at most 0.017 at the tabulated half-lives of 0 to 2 weeks (compound symmetry, mean moderation, E9, $t_{1/2} = 2$: 0.548 to 0.531) and 0.041 on the plotted curves to 4 weeks (same cell, $t_{1/2} = 4$: 0.505 to 0.464); Hybrid at most 0.0016, OL+BDC 0.0020. Monte Carlo change on the same seeds (carryover only): CO 0.018, Hybrid 0.003, OL+BDC 0.003. The header's 0.021 was a stored Monte Carlo change of a reseeded Hybrid cell (0.0216), not a carryover effect. 39's 0.017 is the CO change in `design-power.csv` (grid 0, 0.5, 1, 2), correct. Source: pre-anchored backups against the anchored outputs, and the recursive sandbox rerun on the same seeds (verified; `docs/series/notes/phase1/phase1-cr-log.md`). 37 header and Appendix A and 39 header corrected with dated notes | quote in `ch15` with the grid stated (GP-34 done) |
| CR-42 | Design power ranking | 13 Section 10.1, 14 Section 6.3, 18 (Hybrid > CO > OL+BDC > OL) against 04 Section 5.3, 19.tex (Hybrid > OL+BDC > CO > OL at $c_{bm} = 0.45$) against the published abstract via 23 Section 17.1 and 17 Section 7.1 (CO above Hybrid) | depends on construct, $c_{bm}$ and $N$ allocation (inferred) | compute the ranking under the refined process at matched total $N$ (GP-25); `ch02`, `ch16` |
| CR-43 | Whether `characterize_carryover()` is the liberal F1A procedure, and its grid | 03 Section 6.2, 14 Section 3, 23 Section 9.1 ("validated") against 57 (F1A liberal); 01 (0.1 to 8 weeks) against 03 Section 6.2 (0.25 to 8 weeks) | Resolved 2026-10-11: yes, it is an F1A-type procedure, and a weaker one. `characterize_carryover()` (`R/carryover_analysis.R` lines 27 to 132) fits `lme_analysis()` at each candidate half-life with the interaction in the model, picks the lowest AIC (lines 116 to 118) and returns the model-based $t$-test p-values of the drug effect and the interaction at that half-life as if it were fixed (lines 76, 90; printed at lines 318 to 325). Unlike F1A of 57 it compares REML likelihoods (`nlme::lme` default, `R/lme_analysis.R` lines 167 to 177), which are not comparable across half-lives because $D_{bc,it}$ changes the fixed-effect design; it uses model-based, not CR2, standard errors; and a failed `corCAR1` fit silently falls back to no residual correlation (lines 174 to 177), so the AIC can compare different covariance models. Default grid: 0.1, 0.25, 0.5, 1, 2, 4, 8 weeks (lines 28 to 29), so 01 is right and 03 Section 6.2 (0.25 to 8) is wrong; 0 is not in the grid. Its p-values are liberal (57, 58); "validated" in 23 Section 9.1 is not supported (inspected; code change proposed in `docs/series/notes/phase1/phase1-cr-log.md`, not implemented) | `ch23`, `ch25`; Appendix D (code change) |
| CR-44 | The published dropout form | 13 Section 4 against 17 Section 7.1 item 4 ($\beta_0 + \beta_1(\text{change})^2$) | Resolved 2026-10-11: neither is exact. At `58b32a9` (`R/censordata.R` lines 53 to 86, identical on GitHub and in the local copy) each participant-visit is flagged by two independent uniform draws: a flat flag with probability $\min\{1,\ f_1\,\Delta t_k/\overline{\Delta t}\}$ and a biased flag with probability $\min\{1,\ f_2\,p_{ik}/\bar p\}$, where $p_{ik} = \Delta t_k\,(D_{ik} + 100)^{e_{b1}}$, $\Delta t_k$ is the visit interval (`t_wk`), $D_{ik}$ the summed component improvement $TV + PB + BR$ at visit $k$ (`D_` columns, `generateData.R` line 174), $f_1 = \beta_0(1-\beta_1)$, $f_2 = \beta_0\beta_1$, and the bars are means over all participant-visits of the trial; a visit is dropped if either flag is set, and dropout is carried forward (monotone). 13's $(\text{change} + 100)^{e_{b1}}$ term is right but its direction is wrong: because $D_{ik}$ is improvement, participants who improve more drop out more, contrary to the roxygen text and to 13; 17's "responders more likely to drop out" has the direction right but the form wrong (no additive $\beta_0$, no square of the change alone). $\beta_0$ is the per-visit flag rate, not the total dropout (inspected; direction and per-visit rate verified on synthetic data, `docs/series/notes/phase1/cr44-censor-check.R`) | inspect `censordata` at `58b32a9` (GP-36); Appendix I Section 2, `ch13` |
| CR-45 | $\sigma_v^2$ against $\tau^2$ | 38 and 51 ($v_i$ with $\sigma_v^2$) against NOTATION, P04, 58 ($\mathrm{Var}(v_i) = \tau^2$) | one concept (DQ-6, decided): $\tau^2$ everywhere | NOTATION revision; `ch11`, `ch14` |
| CR-46 | Unexplained OL+BDC Type I error under compound symmetry | 36 Section 9 item 4; 37 Section 6 | open | analysis (GP-32) |
| CR-47 | Silent repair still in the package sampler | 23 Section 18 ("zero silent PD corrections") against 34 Section 4.7, 52 Section 4 (`R/generateData.R:391-393`), 51 Section 20 | the Cholesky fallback still repairs silently (inspected by 34 and 52; not re-checked) | code fix; Appendix D |
| CR-48 | Parity report mode | 22 (parity) Section 9 against the committed quick-mode `parity_report.rds` (52 Section 4) | mismatch stands (52; not re-checked) | regenerate in default mode; Appendix D |
| CR-49 | Manuscript 02 pre-exposure coding of $D_{bc,it}$ and its CO rerun | 45 Sections 8, 12 question 3; 44; 52 | defect confirmed (45 Section 8: CO power 0.484 to 0.744 after correction, compound symmetry); rerun of paper 02 unconfirmed | rerun (GP-31); Appendix C |
| CR-50 | Driver cross-component values in the April comparison; the 19 against 21 Hybrid row | 21 (0.1, 0.05) against 23 Section 6.1 (0.2); 19 (0.74/0.68/0.72) against 21 (0.74/0.72/0.68) | not resolvable (driver archived, output missing) | drop the tables; mention in Appendix I Section 5 only |
| CR-51 | Failed-fit count and DS shares | 48 Purpose (27 failed fits in 291,500) against 47 Section 9 (sum 28); 50 Section 2 (DS: 82%, 86%) against 56 Section 6.1 and 57 Section 3.5 (F1A: 87%, 77%) | Resolved 2026-10-11: no statement is wrong; they use different denominators and analyses. Failed fits: 27 of 291,500 is the 39-cell run of 48 (26,500 trials by 11 analyses); the final 41-cell main run has 28 of 308,000 (the two `prog-BL` cells add one F1 failure); 47 Section 9's per-analysis counts are exact over 14 analyses and 28,000 trials each (55 of 392,000; at most 0.6% per cell, G1 and G2 in cell 1). Selection shares on the 1,000 null trials of each cell: DS (S with AIC) chose 0 in 81.6% at true $t_{1/2} = 0$ and 1 or 2 weeks in 85.7% at true 1 week (50's 82%, 86%); F1A chose 87.0% and 76.7% (56's and 57's 87%, 77%). Source: `replicates-full.rds`, `-extra.rds`, `-extra2.rds` (verified) | state analysis and denominator in `ch21`, `ch23` |

Table: Conflicts that need analysis, inspection or a decision

### 5.3 Needs the pending reruns, or the reporting of finished ones

| ID | Conflict or open number | Documents and sections | State on 2026-10-11 | Action |
|-------|----------------|---------------------------|------------------------------|--------------|
| CR-52 | Cost of F1 relative to S | 47 Section 1, 48 Section 3, 49 Section 1.5, 56 Section 1 (29%, stored run) against 57 Section 8.3 (F1 0.438 under the refined process) | job 15 holds S and F1 under the refined process (observed, 451 rows); not yet tabulated | report in `ch21`, `ch22` (writing only) |
| CR-53 | G1, G2 and F1A under the refined process | 46 Amendment 2, 47 Section 7, 49 Section 5 | local anchored arms unpaired with job 15; jobs 21 and 22 running (decided) | report when jobs 21 and 22 finish; `ch21`. Done 2026-10-15: reported in `docs/47` (rerun section); F1A fails size, G1 and G2 tie with F1 |
| CR-54 | The covariance study under anchored carryover | 45 (recursive, fixed allocation); 35 Section 10 | anchored main pass finished (320 rows, observed); follow-up and $D_{bc,it}$ passes cancelled | report the main pass in `ch20`; state that the follow-up and $D_{bc,it}$ passes stay cancelled and that 45 Sections 6, 7 rest on the recursive passes (DQ-8, decided) |
| CR-55 | The 41-cell selection under the refined process | 46 Amendments 1, 2; 47 Section 7 (F1 by the tie rule) | needs job 15 plus jobs 21, 22 through `12-strawman-stress-decision.R --tag -anchored` | rerun the decision script when 21 and 22 finish; `ch21`. Done 2026-10-15: on jobs 15, 21, 22 G2 led F1 by 0.026; with job 37 pooled (`15-extra4-run.sh`) F1, G1, G2 tie within 0.02 and the tie-break selects F1 (`docs/47`) |

Table: Conflicts that wait on reruns or on reporting finished runs

## 6. Redundancy map

| Cluster | Documents now | Single home | What the others become |
|----------------|----------------------------------------|------------|------------------------|
| R-01 Description of the published process | 24, 25 Section 4, 28 Section 2, 31 Sections 2, 6, 36 Section 3, 37 Section 2, 40 Sections 2 to 6, 41 Section 3, 42 Section 2, 44 Section 2, 53 Section 3, 13 | `ch03` | pointers; glossaries of 53, 55 to 58 point to `fm2` |
| R-02 The April audit (what changed relative to the published code) | 03 Section 2, 04 Section 3 and Appendix A, 13 scope box, 14 Section 4, 23 Section 2, CM 2024-revision subsections | `ch03` (attribution) and Appendix I Section 1 (per-commit numbers) | archived |
| R-03 Differences between the package and the published process | 23 Section 2, 32 Section 2, 33 Section 3.4, 54 Section 1, 53, 55, 36, NOTATION | `ch10` difference table | archived or pointer |
| R-04 Positive definiteness and the ceiling | 01 Section 6, 03 Section 3, 04 Section 4 and Table A.5, 07, 07a, 19.tex Section 7, 23 Section 3, 31 Section 20, 32 Sections 3.3, 4.3, 8, 34, 35, 36 Sections 4, 6.5, Appendix A.7.4, 37 Sections 5.5, 5.8, 44 Step 6, 46 Section 3, 50 Section 1, 51 Feature 5c, 52 Section 2, 54 Sections 4 to 5, 55 Section 8.2, CM | `ch09`; derivations in CM | removed; one-line pointers |
| R-05 Kronecker decompositions and generative forms | 35 Sections 4.2, 4.4 to 4.7, 36 Sections 6.4, 6.7, Appendix A.7.5 to A.7.6, 37 Section 5.8, 38 Section 5, 44 Step 5 | `ch09` | removed |
| R-06 Step coupling and the collapse | 04 Section 3, 08 Sections 1, 5, 09, 09a.md Section 5.5 note, 13 scope box, 23 Section 5.7, 34 Section 4.6, 36 Section 5, 37 Sections 5.1 to 5.4 and Appendix A, 40 Section 5, 44 Step 2, 50 Section 1, 51 Feature 5, 55, CM "Step coupling" | `ch04` (diagnosis), `ch27` (note) | others archived |
| R-07 Graded coupling from a pharmacokinetic argument | 04 Section 3.5, 08, 14 Section 3.1, 23 Section 5, 55 Section 3, CM `sec:inva-pk` | `ch08` | 08 credited, archived |
| R-08 The compounding recursion | 08 Section 2.2, 09 Section 2.2, 13 Section 3.4, 51 Feature 3, 53, 54 Section 3.5, CM "Anchored versus compounding" | `ch06` (mean), `ch26` (note) | removed |
| R-09 The two power tables of 53 and 55 | 53 Section 8.3 and Table 8, 55 Table 8 and Section 10 Table 11 | `ch28` (one table) | 53 and 55 point to it |
| R-10 Mean against covariance moderation | 02, 19 (.md and .tex), 01 Section 1, 03 Section 5.4, 21, 22 (journal), 23 Sections 1, 8, `plan-architecture-c.md`, 33, 36 configurations D and E, paper 01, CM Part II | `ch07`; Appendix I Section 5 (history) | archived |
| R-11 Gompertz curve and parameters | 23 Section 2.6, 24, 25, 28 Section 2.1, 31 Section 2, 40 Section 4, CM | `ch05` | removed; three notations reduced to $m$, $d$, $r$ |
| R-12 Component decomposition pedagogy | 24 (whole) = 31 Parts I, III Section 8, VI; 40 Sections 3 to 4, 11; CM; paper 06 | `ch05` and `ch12` | 31 Parts I and VI deleted |
| R-13 Hybrid design, schedule and expectancy | 24 Component 2, 27 Section 2c, 28 Section 3, 31 Sections 5, 7, 33 Section 2.5, 37 Section 2, 40 Sections 2 to 3, glossaries of 53, 55 to 58, CM | `ch02` | pointers |
| R-14 Paper 06 results | 28 Section 4, 31 Section 11, paper 06 | `ch12` | 28 archived |
| R-15 Pharmacodynamic carryover and the response half-life | 38, 44 Steps 7 to 8, 40 Section 4, 39, 38 Section 2 = 40 Section 3 = 45 Section 7 (prazosin) | `ch06` (half-life), `ch14` (model) | removed |
| R-16 E9 and E9b closed form | 17 Section 3.3, 37 Sections 3 to 6, 39 Sections 2 to 3, 40 Section 8, 41 Proposition 4 and Section 6, 43 Section 3.5, CM | `ch15` | removed |
| R-17 Design efficiency and cycles | 23 Section 13.1 item 4, 27 Section 4, 34 Section 5, 39, CM | `ch16` | removed |
| R-18 Analysis model and its evolution | 23 Section 4, 24, 25 Section 6b, 26 Section 9, 27, 28, 31, 06, 43, 45 to 49, 56, 57 | `ch17` (model), Appendix B (evolution) | removed |
| R-19 Time adjustment | 40 Sections 9 to 10, 42 (whole), 45 Section 4.4, 43 time-adjustment rows | `ch17` | 42's tables and figures not duplicated |
| R-20 Baseline coding | 41, 43 Section 2.5, 45 Sections 6, 10, CM | `ch17` | removed |
| R-21 Alternative analyses and test methods | 02 Section 6, 17 Section 3.3, 23 Sections 10, 13, 26, 43, `plan-coherent-analysis-comparison.md`, CM | Appendix A (catalog), Appendix H (RM-ANOVA), `ch17` | removed |
| R-22 Shared slope and prognostic leak | 41, 45, 48 Sections 1 to 3, 49 Sections 2 to 3, 51 Feature 6, 56 Section 3.3, CM | `ch19` | 48 Sections 1 to 2 kept as the lead-in only |
| R-23 Expectancy-drug confounding | 40 Section 8, 49 Section 2, 51 Feature 1, 56 Section 3.3, 58 Section 6, CM | `ch19` | removed |
| R-24 Stress-test design and results | 46, 47, 48 Section 5, 49 Sections 1, 4, 5, 50 Section 2, 56 Sections 5, 7, CM | `ch21` | 46 as annex; others removed |
| R-25 Symmetry screen | 46 Amendment 1, 47 Section 5, 48 Section 4, CM "Sign symmetry of power" | `ch21` (one section) | removed |
| R-26 AIC choice of the half-life | 47 Section 8, 49 Section 5, 50 Section 2, 56 Section 6, 57, 58 Section 9, CM | `ch23` (interaction), `ch24` (main effect) | removed |
| R-27 Recommendation and recipe | 49 Section 7, 50 Section 2, 56 Sections 7 to 8, 57 Section 7.2 | `ch22` | removed |
| R-28 "Superseded by F1" notes | 35 Section 10, 41 Sections 1, 8, 42 update, 43 Sections 1, 5, 44 Sections 1, 5, 45 header and Sections 1, 9 | Appendix B (one entry) | removed |
| R-29 Package implementation status | 34 Section 4.7, 35 Sections 9, 10, 36 Sections 6.7, 9, 44 Sections 3, 4, 50, 51 Section 20, 52 Section 4 | Appendix D | removed; chapters link to it |
| R-30 Repository, parity and tooling history | 00, 01, 03 Section 5, 15, 16, 22 (parity), 23 Sections 7, 14, 15, `CLAUDE.md` | `CLAUDE.md` and `analysis/scripts/parity/README.md`; Appendix I Section 10 | archived |
| R-31 Bibliography | 20, 23 Section 17, reference lists of 02, 09a, 17 and of 34 to 58, `docs/references/` | Appendix G and one `.bib` | per-chapter lists generated from the `.bib` |
| R-32 Terminology authorities | NOTATION, 29, 30, `CLAUDE.md` lexicon section, every notation index | NOTATION (symbols), `fm2` (project terms), Appendix E (lexicon), Appendix F (style) | `CLAUDE.md` section points to them |
| R-33 Main-effect half-life mathematics | 58, P04 Revision 1 tables, CM main-effect section | `ch24`; P04 keeps operations only | P04 points to `ch24` |

Table: Redundancy clusters and their single home

## 7. Notation and convention fixes

### 7.1 Documents before 34: wholesale

These documents are archived rather than edited (Section 2.9), so the
fixes apply when their content is carried into a chapter or into
Appendix I.

| Fix | Where it occurs (survey) | Canonical form |
|----------------------|----------------------------------------|------------------------------|
| Architecture letters | all `.tex`; 21, 23, 32, `plan-architecture-c.md` | mean-moderation, covariance-moderation, dual-channel |
| "N-of-1" for the Hybrid design | 13, 14, 18 | Hybrid |
| ER, TR, "pharm_biomarker" for components | 07, 07a, 09, 13, 15, 18 | $PB$, $TV$ |
| $S_{i,t}$, `Sx` as a symbol | `.tex` set | $Y_{it}$ |
| $\rho$ for compound symmetry | 07a, 18, 23 Sections 2.1 to 3.3, 31 Section 20, 32 ($\rho_c$) | $\rho_{cs}$ |
| $\rho$ or $\rho_t$ for the biomarker-response correlation | 04 Section 3.1, 08, 09 | $c_{bm}\,g_t$ |
| $c_{bm}$ as the mean-moderation multiplier | 01, 19.tex | $\beta_{bm}$ |
| $B_i^{\ast}$ for the standardized biomarker | 01, 19.tex | $b_i$ |
| $D_{bc}$ for a binary indicator | 06 Section 4.1, 09 | $D_{it}$ |
| $c_{tv}$, $c_{pb}$, $c_{br}$ for autocorrelations; $c_{cf1t}$, $c_{cfct}$ | 07, 07a, 24, 31 | $\rho_{cs}$ or $\rho$; $c_1$, $c_\times$ |
| $e_t$ for expectancy | 07 | $\eta_t$ |
| $\alpha$ for the decay factor or a power-law exponent; $\lambda_{cor}$, $\lambda_{drug}$ | 08, 09 | $\phi_t$; $\lambda$ |
| $f_G$, $a$, $g_b$, $g_c$ for the Gompertz curve | 08, 25 | $G(t_{od}; m, d, r)$ |
| $\beta_t$, $\beta_1$, $\beta_2$ for time and biomarker coefficients | 23 Section 4.1, 24, 25, 31 | $\beta_T$, $\beta_B$ |
| $K$ for blocks, visits or components; $n$ for occasions | 17, 23 Section 13, 24, 31 Section 20 | $n_c$; $n_{\text{on}}$, $n_{\text{off}}$; $K$ only for the cross-component matrix |
| $g$, $\phi_j$, $\tau_t$ in RM-ANOVA | 26 | another stratum letter; $\omega_j$ or similar for phase effects; another letter for time effects |
| $v_i$ for a random slope | 24 | another letter ($v_i$ is the drug-effect deviation) |
| $\tilde r$, $v$, $r$ for the coupling vector | 31, 32, 34, 36 | $c_{bm}\tilde g$, $\tilde g$ |
| $g_t$ as the off-drug factor only | 32 | $g_t$ is the whole profile (1 on drug) |
| $N(\cdot)$ for the normal law | 09a.md, 25, 26 | $\mathcal N$ |
| local labels for now-canonical symbols ($\phi_t$, $g_t$, $\mu_{bm}$, $w_t$, $c_1$, $M$, $m$, $d$, $r$) | 23, 24, 32, 33 | marked canonical |
| "biomarker-by-treatment", "washout" for the off-drug span, lower-case "hybrid" | 21, 23, 24, 28, 31 | biomarker-treatment interaction; off-drug period; Hybrid |
| hand-written "Rendered on" footers with `~/prj/alz/...` | every `.tex`; `reorg-plan-2026-04-22.md` | removed; `tools/render.sh` owns provenance |
| British spellings | `reorg-plan-2026-04-22.md` | US English |

Table: Notation and convention fixes for content carried from
documents before docs/34

### 7.2 Regressions in 34 to 58

| Document and location (survey) | Problem | Fix |
|-------------------------|------------------------------|------------------------------|
| 34 Sections 3.1 to 3.2 | $\tilde r$; unpadded $g$ in the ceiling formula; $C$ as cross block | $\tilde g$; $(\tilde g^\top M^{-1}\tilde g)^{-1/2}$; name the block explicitly |
| 35 | $q_K$, $q_A$ (eigenvectors); $\Delta$ as a time lag; $\mu$ as an eigenvalue; $v_{TV}$, $v_{PB}$, $v_{BR}$; $s_t$ as factor SDs; $c_{bm,BL}$ marked local | letters chosen against NOTATION's reserved list; $c_{bm,BL}$ marked canonical |
| 36 Section 6.7, Appendix A.7.5 | $q$ for the quadratic form; $r$ for the coupling vector; $K$, $A_c$ marked local | another letter; $c_{bm}\tilde g$; canonical |
| 37 Sections 4, 4.4 (lines 492, 624, 635) | plain $E(\cdot)$; $\mu$, $\mu_p$, $\bar\mu$ local; $v_t$ factor SDs; $e_i$ residual; "patient(s)" | $\mathbb E$; another letter; participant |
| 38 | $f(t)$, $\bar f$, $F$; $G$, $H$ in the Section 6 proof; $\sigma_v^2$ | $h(t)$, $\Delta_h$ (as 44); other letters; $\tau^2$ (DQ-6) |
| 39 | $f(d)$; $\omega$ (measurement-error share); $C$ (response covariance); $v$; $d$ (spacing); $w_t$ marked local; $r_{\Delta B}^2$ for a non-correlation; "Values changed by at most 0"; "patient(s)" | other letters; rename the criterion; fix the phrase; participant |
| 41 line 205 | $N(\cdot)$; $f(t)$; $d_{it}$ | $\mathcal N$; another letter; another letter |
| 42, 45 | $f(t)$ for the time adjustment | another letter (for example $s(w_t)$) |
| 44 index | stale description of 38's symbols | NOTATION's |
| 45 index | $B_i$ "centered" | $B^c_i$ |
| 48 Section 3 note; 49 Section 3 note | "published construct" for AR(1) variants | configuration name with $\rho$ |
| 49, 56 | $n_{\text{on}}$, $n_{\text{off}}$, $n_v$ (49) against $n_0$ to $n_3$ (56) | NOTATION's $n_{\text{on}}$, $n_{\text{off}}$, $n_v$ |
| 51 | $\sigma_v^2$; $E(t)$, $E_{it}$ beside E1 to E9; $L$, $U$ beside $L_{it}$; file name "ten" | $\tau^2$ (DQ-6); other letters; rename |
| 53 | $m$ (run position), $L$ (run length) | other letters |
| 54 | $w_t$ marked local; "occasion"; $D_{bc}$ without subscripts; $n$ for occasions | canonical; visit |
| 56 Section 1 | AR(1) exponent on the visit index | $\lvert w_t - w_s\rvert$ |
| 57 | $M$, $m$ (bootstrap); $g(Y)$; $\gamma$; $h$ (candidate half-life); $\omega$ ($\log t_{1/2}^a$); $\hat\xi$; $e_i$; grid $\{0, 0.25, 0.5, 1, 2\}$ | letters chosen against NOTATION's reserved list (for example $R$ resamples, $\mathcal S(Y)$ for the selection rule); one grid across 57 and 58 (DQ-10) |
| 58 | $B$ (GLS bread), $G$ (randomization distribution) beside the Gompertz $G$, $\hat K$, $Z$, $r_p$, $f_m$, $M$, $h$; grid $\{0, 0.2, 0.5, 1, 2\}$; analysis labels shared with the interaction track | other letters; `M-` labels (Section 2.7) |

Table: Notation regressions in docs/34 to docs/58

### 7.3 Cross-cutting rules for the series

1. NOTATION is the only symbol authority; a chapter adds only local
   symbols, never on a reserved letter.
2. Every result states its construct in the chapter's construct box
   (Section 2.10); a ceiling, size or power number without $\rho$, the
   construct, carryover, coupling, allocation and $N$ is not quoted.
3. $N$ is the total; $t_{1/2}$ in weeks; the fixed analysis half-life
   is 0.2 week unless a historical result is quoted (decided).
4. Evidence labels from Section 2.5 only.
5. Analysis IDs from the registry only (Section 2.7).
6. Each table has a caption line, each figure a caption with its
   generating script; figure files named in print order.
7. No dated correction paragraphs at release (Section 2.9).
8. US English; no em dashes; prose vocabulary per Appendix F; this
   project's studies are exploratory and its protocols are working
   protocols, never described with the vocabulary of registered trials.
9. The linter `tools/notation-lint.pl` runs over every chapter; two
   checks are added for the series: plain $E($ and $N($ for the
   normal law, and reserved letters reused in a chapter's local table
   (GP-37).

## 8. Gaps

(Updated 2026-10-11: GP-29, GP-34, GP-35 and GP-36 are done; see the
rows and Section 5.2.)

Merged from the four surveys' gap lists and deduplicated. "Writing" means
the content exists in outputs or derivations and needs only to be
written; "computation" means new runs or code; "blocked" means it waits
on a running job or a decision. Costs are rough working time with
assistance; server time is stated separately.

| ID | Gap | Status | Rough cost |
|-------|-----------------------------------------------------|----------------------|--------------|
| GP-01 | Front matter: reading guide, executive summary, project glossary | writing | 1 day |
| GP-02 | One table of differences between the published and the refined process | writing | 0.5 day |
| GP-03 | Commit attribution table (`42ac030` to `06dac83`), checked by `gh api` | writing, light inspection | 0.5 day |
| GP-04 | One design table per design, including balanced placebo and cycle designs; Hybrid schedule figure | writing | 1 day |
| GP-05 | The full ceiling grid (construct by $\rho$ by design by $t_{1/2}$) | computation, closed form (existing `cbm-ceiling` scripts) | 2 hours |
| GP-06 | Estimands and units per architecture; the mean-moderation estimand shift | writing (33 Sections 2.3, 3.2) | 0.5 day |
| GP-07 | $\kappa$ by design and architecture | writing from the 33 Appendix 1 reanalysis | 0.5 day |
| GP-08 | Plasma against pharmacodynamic half-life | writing | 2 hours |
| GP-09 | Answers to the co-author's ten questions | writing | 2 hours |
| GP-10 | Trait against state variance off drug, as an open choice | writing | 2 hours |
| GP-11 | Operating characteristics of the refined process as a whole | writing from jobs 14 and 15 | 1 day |
| GP-12 | Stress test under the refined process, all 14 analyses | done 2026-10-15 (`docs/47`, rerun section; jobs 15, 21, 22 and the post hoc job 37): F1 retained, G1 and G2 tie as post hoc alternatives | |
| GP-13 | Covariance study under anchored carryover | writing (main pass); the follow-up and $D_{bc,it}$ passes stay cancelled (DQ-8, decided) | 0.5 day |
| GP-14 | F1 cost relative to S under the refined process | writing (job 15) | 2 hours |
| GP-15 | Full bootstrap ($M = 999$) check of F1A0B in a few cells | computation: about 17 worker-minutes per trial (57 Section 8.4), so 3 cells of 200 trials is about 170 worker-hours | server |
| GP-16 | Sample-size planning for F1 at 0.2 week and F1A0B (power by $N$, $\rho$, $c_{bm}$) | computation | 2 days plus server |
| GP-17 | The interaction analysis in CO and OL+BDC (expectancy terms constant or confounded differently) | computation | 2 days plus server |
| GP-18 | Prognostic couplings at strengths other than 0.3 and combined pathways | computation | 1 day plus server |
| GP-19 | Simulation of Features 11 to 17 of 51 | computation (large) | defer |
| GP-20 | Calibration of effect sizes and half-life to prazosin trial data | blocked on data access | unknown |
| GP-21 | Main-effect stress test under the rebuilt P04 | computation (paper 04's plan) | per P04 Section 7.3 |
| GP-22 | Coverage of profile-likelihood and PoSI-type intervals | computation | 1 day plus server |
| GP-23 | Monte Carlo methodology section | writing | 0.5 day |
| GP-24 | Dropout chapter under the refined process | writing (CM, paper 09, stress-test MAR and MCAR cells) | 1 day |
| GP-25 | Design ranking under the refined process at matched total $N$ | computation (closed form plus a small Monte Carlo) | 1 day |
| GP-26 | Figures: Gompertz family, Hybrid schedule, coupling profiles, ceiling against $\rho$; optional Figure 5 regeneration | writing and plotting from existing scripts | 1 day |
| GP-27 | Implementation status appendix (separable `buildSigma`, factored generator, Cholesky fallback, `cbmCeiling()`, anchored carryover, prognostic couplings in the package) | writing; the code work is separate | 0.5 day |
| GP-28 | CM updates: main-effect track (version 1 to P04 and 58); 57 procedures; 55 cited; anchored stress-test DGP; 58 Propositions 1 to 7 | writing (LaTeX) | 2 days |
| GP-29 | Verification of four suspect references (Chen, Li and Yang 2020; Eng and Konings 2018; Pequignat et al. 2023; De Bock et al. 2021) | done 2026-10-11 (PubMed and Crossref): all four misattributed, none exists as cited. "Chen, Li and Yang": the title is Blackston, Chapple, McGree, McDonald and Nikles (2019), *Healthcare* 7(4):137, doi:10.3390/healthcare7040137; the cited DOI 10.1371/journal.pone.0216949 is Ao et al. (2019), an Ebola-glycoprotein HIV vaccine paper, *PLoS One* 14(5). "Eng and Konings": the title and DOI are Haller and Ulm (2018), *Trials* 19:128 (not 19:72), doi:10.1186/s13063-018-2491-0. "Pequignat et al.": the title and DOI are Tueller, Ramirez, Cance et al. (2022), *Behav Res Methods* 55(8):4175 to 4199, doi:10.3758/s13428-022-02012-1 (no Pequignat author in PubMed for this topic). "De Bock et al.": the title and DOI are Natesan Batley and Hedges (2021), *Behav Res Methods* 53(4):1782 to 1798, doi:10.3758/s13428-020-01522-0. The annotations of 20 for these entries describe the fabricated attributions and are not to be carried | Appendix G corrections |
| GP-30 | The unreported common-grid contrast (`36-fig4-contrast-reps100.png`) | writing, if `13-summarize-contrast.R` outputs suffice (inferred) | 2 hours |
| GP-31 | Paper 02 CO rerun with the corrected pre-exposure coding | computation | 0.5 day plus server |
| GP-32 | The OL+BDC Type I error under compound symmetry | analysis | 1 day |
| GP-33 | Paper 11 boundary cells, analytic super-additivity, MCSE columns (`CLAUDE.md`) | computation | 1 day plus server |
| GP-34 | Reconcile 37's anchored-change figures | done 2026-10-11 (CR-31) | |
| GP-35 | `characterize_carryover()` against F1A; its grid | done 2026-10-11 (CR-43); code change proposed, not implemented | |
| GP-36 | The published dropout form | done 2026-10-11 (CR-44) | |
| GP-37 | Linter checks for $E($, $N($ and reserved letters | code | 2 hours |
| GP-38 | A single open-items register | writing | 0.5 day |
| GP-39 | F1 at 0.2 week in the 41-cell stress test with the full screens of 46 (only the job 14 extension evaluated F1h02) | computation, optional: the job 14 evidence covers size in 58 cells including the prognostic ones (quoted, 57 Section 8.2); maximin under the 46 rule at 0.2 week is not reported | 0.5 day plus server |

Table: Gaps the series must fill

## 9. Migration plan

### 9.1 Phases

| Phase | Tasks | Depends on | Verification | Effort |
|----------|------------------------------------|------------|------------------------|--------|
| 0. Freeze and archive | create `docs/legacy/` and `docs/series/`; tag the repository (`whitepapers-pre-series`); move the archived files with `git mv`, unchanged except for a banner; move the files removed from the series; move the slide deck to `docs/archive/` with an unreviewed label; delete the two dated duplicate PDFs; write `docs/legacy/INDEX.md` (old number to new home) | DQ-1, DQ-2, DQ-4, DQ-9 (decided) | every file in Section 4 accounted for (a script lists `docs/` before and after); `git log --follow` works for a sample | 0.5 day |
| 1. Resolve conflicts and report finished runs | settle CR-31, CR-43, CR-44, CR-51 by inspection; tabulate job 15 and the covariance anchored pass; report jobs 21 and 22 when they finish; rerun the decision script on the refined process (done 2026-10-15: `docs/47`, rerun section; F1 retained after the additional trials of job 37); GP-29 references | Phase 0; jobs 21, 22 | each CR row has a stated answer and source; decision files `*-anchored.csv` produced | 3 days plus job time |
| 2. Skeleton and content moves | write the chapter template and registries (`chapters.csv`, `analyses.csv`); create every chapter with its outline; move content from the spines; write Appendix I; align driver defaults with the 0.2 week policy (inspect `t_half_a` defaults first) | Phase 1 for Parts IV and V; Phases 0 and 1 for Appendix I | every source section of Section 3 mapped to a chapter section (a checklist generated from the tables of Section 3) | 8 days -- **Done 2026-10-15**: 41 chapter files (fm0-fm2, ch01-ch28, appA-appJ) plus `chapters.csv` and `analyses.csv` committed; 13 spine chapters have verbatim spine text; all others have source content or TODO stubs; `TH_A` driver default aligned to 0.2 week. Appendix I: 11-section legacy consolidation (1,107 lines). Remaining TODOs (new writing, unverified entries) pass to Phase 3. |
| 3. New writing | GP-01 to GP-14, GP-23, GP-24, GP-26 to GP-28, GP-30, GP-38 | Phase 2 | each chapter passes the template checklist | 10 days -- **Substantially done 2026-10-15**: GP-01 to GP-14, GP-23, GP-24, GP-27, GP-28 (partial), GP-30, GP-38 done (commits fd766c0a to dbf3558e). Deferred: GP-09 (co-author PDF-blocked), GP-26 (figures, R plotting), GP-28 full (docs/57 procedures, docs/58 propositions 1-7), GP-37 (linter code). |
| 4. Notation and convention pass | Section 7; linter with the new checks (GP-37); evidence labels; construct boxes; analysis IDs; prose per Appendix F | Phase 3 | `perl tools/notation-lint.pl docs/series/*.md` clean; a grep for retired labels returns nothing outside quotations | 3 days -- **In progress 2026-10-15**: linter clean across all 43 series files (commits 3a10bfef); construct boxes filled in ch10, ch13, ch18, ch20, ch21, ch22 (confirmed or in-flight for ch23/ch24); evidence labels pass started on ch21. Remaining: evidence labels in ch20, ch22 and other numerical chapters; construct boxes for all other chapters; GP-37 (linter extensions). |
| 5. Render, link check, retire | render every chapter through `tools/render.sh`; check every `docs/NN`, `series/chNN` and CM label reference resolves; update manuscripts, CM, NOTATION and `CLAUDE.md` (Sections 9.4 to 9.5); archive this plan | Phase 4 | link checker over `docs/`, `analysis/report/`, `analysis/scripts/`; PDFs current | 2 days |

Table: Migration phases (computational gaps GP-15 to GP-22, GP-25,
GP-31 to GP-33 and GP-39 run alongside, on the server)

### 9.2 File layout and naming

Decided (DQ-2, 2026-10-11): a new folder `docs/series/` with identifiers that
cannot be confused with the old numbers.

```
docs/
  series/
    fm0-reading-guide.md        fm1-executive-summary.md
    fm2-notation-conventions.md
    ch01-the-question.md ... ch28-patching-58b32a9.md
    appA-analysis-catalog.md ... appJ-software-pointer.md
    registry/chapters.csv  registry/analyses.csv
    registry/open-items.csv
  legacy/
    INDEX.md                    (old number -> new home)
    00-documentation-index.tex ... 59-whitepaper-refactoring-plan.md
  figures/                      (shared; new files chNN-figK-*.png)
  archive/                      (non-whitepaper history; the slide deck)
  references/  reviews/  tooling/
```

The alternative, keeping the old numbers and adding new documents as
60 and up, was rejected: numbers would not follow reading order, and
the series would be interleaved with the documents it replaces.
Renumbering old documents in place was rejected because 118 files
outside `docs/` cite `docs/NN` (observed, grep: 71 scripts, 46 files
under `analysis/report/`, and `CLAUDE.md`).

### 9.3 How old numbers keep resolving

Old numbers are permanent identifiers. `docs/NN` keeps meaning "legacy
document NN"; the file lives at `docs/legacy/NN-*.md` (or `.tex`) and
`docs/legacy/INDEX.md` maps each number, and each section where a
document was split, to its new home. The archived file carries one
banner line ("Archived 2026-MM-DD; current treatment: `series/chNN`
Section K; see `docs/legacy/INDEX.md`"). No redirect stubs are left in
`docs/`; a stub per old file would put 58 files back into the folder
the refactor is clearing. `fm0` prints the same map for readers.

### 9.4 Updating citations in manuscripts, CM and scripts

- **Scripts (71 files).** Their `docs/NN` comments are provenance and
  keep pointing to the legacy document, which is correct by
  Section 9.3. No edit is needed; a one-line note in the script
  header convention says how to resolve them.
- **Manuscripts (46 files under `analysis/report/`, including
  `bullets*.Rmd`, `report-slim.Rmd`, the `.tex` outputs and the
  `revisions/` files).** Citations to current content are changed to
  `series/chNN` with section; citations that refer to history (for
  example a referee report's quotation) keep `docs/NN`. The `.tex`
  files are regenerated by rendering, not edited. The edit list is
  produced by a grep and applied with the Edit tool (the repository
  path may resolve under cloud storage; in-place stream editors are
  not used, per the user's rule).
- **CM.** Its citations of docs 46 to 54 move to chapter references;
  it gains citations of `ch23`, `ch24`, `ch27` (GP-28).
- **NOTATION.** The "Used in" columns are rewritten from `docs/NN` to
  `series/chNN` in one revision with a dated entry; "Prose
  terminology" points to Appendix F.
- **P04.** Its 76 citations of `docs/58` (survey d) become `ch24`; its
  analysis labels take the `M-` prefix (DQ-5).

### 9.5 Keeping the project `CLAUDE.md` in sync

`CLAUDE.md` is project configuration, and this plan does not edit it.
The thirteen stale or missing statements this plan found are listed
below. Items 1 to 12 were corrected in `CLAUDE.md` on 2026-10-11 by the
user (decided; not re-read here). Item 13 is open: it points to
`docs/series/`, which does not yet exist, so `CLAUDE.md` carries for
now a "Documentation map" paragraph pointing to this plan instead; the
paragraph is replaced by the pointer of item 13 in Phase 5.

1. Corrected 2026-10-11. "Key architectural feature": "40-60% for
   Architecture B vs. 8-13% for Architecture A" was stale (CR-05).
2. Corrected 2026-10-11. "Architecture A", "Architecture B",
   "Architecture C" throughout: retired labels (NOTATION Part 3).
3. Corrected 2026-10-11. "expanded feasible range for `c.bm` from max
   0.25 under CS to 0.49+ under AR(1)": stale (CR-04), as `docs/34`
   Section 6.2 had asked.
4. Corrected 2026-10-11. "Analysis model ... Fixed formula:
   `Sx ~ bm + t + Dbc + bm:Dbc`" with a random intercept: no longer the
   recommended analysis (CR-10).
5. Corrected 2026-10-11. "`carryover_t1half` ... (default 0.5-1.0 for
   the PTSD/prazosin context)": conflicted with the 0.2 week analysis
   policy and did not distinguish the data-generating half-life from
   $t_{1/2}^a$ (CR-37).
6. Corrected 2026-10-11. The repository layout listed
   `analysis/reports/` with three unnumbered slugs, and "Render a
   manuscript" named three slugs and called `rmarkdown::render()`
   directly, which the user's rule for `~/prj/res` forbids.
7. Corrected 2026-10-11. "`docs/figures/` PNG assets referenced by
   `docs/04-*.tex`": stale.
8. Corrected 2026-10-11. "See
   `docs/19-mean-moderation-implementation-notes.tex`": pointed to the
   unreviewed `.tex`; after Phase 0 the reference is Appendix I
   Section 5.
9. Corrected 2026-10-11. "Related external documents": the
   reorganization plan overstated as producing the current layout;
   `docs/23` described as current; `docs/plan-architecture-c.md`.
10. Corrected 2026-10-11. "Hendrickson et al. reported 2.5-5.8%":
    from `docs/03`, unverified against the publication.
11. Corrected 2026-10-11. British spellings.
12. Corrected 2026-10-11. The "Last updated" date.
13. Open. A short section pointing to
    `docs/series/fm0-reading-guide.md` and to `docs/legacy/INDEX.md`,
    added in Phase 5 in place of the interim "Documentation map"
    paragraph.

## 10. Decisions

All ten choices this plan put to the user were decided on 2026-10-11,
in addition to the three decisions of Section 2.11. Each is recorded
here with what it settles; the rest of the plan is written to them.

1. **DQ-1. Placement of the mathematics compendium.** Decided
   2026-10-11: CM is the series' mathematical reference, Volume II,
   the single home of derivations, cited by label from each chapter's
   "Derivations" section (Section 2.3).
2. **DQ-2. Naming and layout.** Decided 2026-10-11: `docs/series/`
   with `fm`, `ch` and `app` identifiers, and `docs/legacy/` holding
   the old whitepapers under their old names, with an old-to-new index
   (Sections 9.2 to 9.3).
3. **DQ-3. Docs 53 and 55.** Decided 2026-10-11: both stay standalone
   notes for the original authors (`ch26`, `ch27`), and a series
   chapter, `ch28`, summarizes both, with the patches and the single
   power table; their duplicated tables are removed.
4. **DQ-4. The slide deck.** Decided 2026-10-11:
   `Architecting_Precision_Medicine.pdf` is moved out of the top level
   of `docs/` into `docs/archive/`, not deleted from the repository,
   and is labeled there as unreviewed (an AI-generated deck of
   12.2 MB with no text layer, cited nowhere).
5. **DQ-5. Analysis IDs.** Decided 2026-10-11: the `M-` prefix marks
   main-effect analyses only (`M-F1`, `M-VF`, `M-RV`, ...), in P04,
   `docs/58` and `ch24`, before paper 04's stored data exist; the
   interaction labels stay as stored; documents name the track at
   first use ("interaction analysis F1"). `M-` with a hyphen does not
   collide with the paper-local `10-M0` to `10-M3` labels of Section
   2.7.
6. **DQ-6. $\sigma_v^2$ and $\tau^2$.** Decided 2026-10-11: they are
   the same concept, the between-participant variance of the drug
   effect; $\tau^2$ is written everywhere (in `ch11`, `ch14` and every
   carried passage of 38 and 51), and the identification is recorded
   in NOTATION.
7. **DQ-7. Name of the estimand.** Decided 2026-10-11:
   "biomarker-treatment interaction" everywhere, replacing
   "biomarker-by-drug" in the titles and text carried from 43 to 57.
8. **DQ-8. The covariance-study passes.** Decided 2026-10-11: `ch20`
   reports the anchored main pass; the follow-up and $D_{bc,it}$
   passes stay cancelled, and `ch20` states this, together with the
   fact that 45 Sections 6 and 7 rest on the recursive passes and that
   the recommendation (F1 at 0.2 week, F1A0B) does not depend on them.
9. **DQ-9. Location of the archive.** Decided 2026-10-11:
   `docs/legacy/` inside the repository, so that `docs/NN` keeps
   resolving and history stays in git.
10. **DQ-10. Half-life grid.** Decided 2026-10-11: new runs use one
    grid for both tracks, $\{0, 0.2, 0.5, 1, 2\}$; the stored
    interaction runs keep $\{0, 0.25, 0.5, 1, 2\}$, stated wherever
    their results are quoted (in the construct box).

## 11. Evidence basis and limits of this plan

- **What rests on the surveys.** Most per-document findings (status,
  section contents, notation problems, citation counts) are taken from
  the four survey files and were not re-checked; they are labeled
  survey where cited.
- **What was spot-checked** against the source on 2026-10-11: `docs/03`
  Section 2.2 (lines 57 to 58); `docs/31` Section 18.1 (line 1837);
  `docs/37` header (lines 16 to 17) and Appendix A (line 1662);
  `docs/47` Section 1 (line 123); `docs/48` Section 4 (line 384);
  `docs/54` Section 6 (lines 876 to 880); `docs/55` (lines 444, 464,
  685); `docs/56` Section 9 (lines 1004, 1022); `docs/57` Sections 1
  (lines 136 to 234), 8.1 and 8.2, and line 1544; `docs/58` line 222;
  `RERUN-anchored.md` Sections 1, 2 and 7; the analyses and row counts
  of `summary-full-anchored.csv`, `summary-full-extra-anchored.csv`,
  `summary-full-extra2-anchored.csv` and
  `summary-alt250-null500-anchored.csv`.
- **Not checked.** The contents of `docs/references/` (including
  `pmsimstats-nof1-whitepaper.md`), `docs/archive/` and the 10 pages of
  the slide deck after the first; whether CM already holds 54 Section
  5 in full; the code defaults for $t_{1/2}^a$; whether job 15 used the
  same cells and seeds as the stored run beyond the file layout; the
  page estimates, which are judgments.
- **Counts.** File counts are from `ls` of `docs/` (115 top-level
  files); the citing-file counts are from a grep for `docs/NN` over
  `analysis/report/`, `analysis/scripts/`, `R/` and `CLAUDE.md`.
- **Nothing was edited** besides this file, apart from pruning the
  earlier staged copies of this plan in `analysis/report/share/` and
  their `MANIFEST.md` rows; no simulation was run.

## References

Giacomini R, Politis DN, White H (2013). A warp-speed method for
conducting Monte Carlo experiments involving bootstrap estimators.
*Econometric Theory* 29(3):567-589.

Morris TP, White IR, Crowther MJ (2019). Using simulation studies to
evaluate statistical methods. *Statistics in Medicine*
38(11):2074-2102.
