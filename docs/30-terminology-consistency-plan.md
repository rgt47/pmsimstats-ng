# Manuscript Terminology Consistency Plan {.unlisted .unnumbered}
*2026-06-17 19:37 PDT*

**Author.** pmsimstats team

**Purpose.** A plan to align the terminology in all `report.Rmd` and
`report-slim.Rmd` manuscripts with the project lexicon
(`docs/29-nof1-precision-medicine-lexicon.md`) and with each other,
together with the record of its execution on 2026-07-30. The plan
governs prose surface forms and acronym expansion; mathematical
symbols and code identifiers are governed by
`analysis/report/NOTATION.md`.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

This plan is about prose terminology, so it uses almost no
mathematics. The symbols and labels that do occur are listed here
with their status against `analysis/report/NOTATION.md`.

| Symbol or label | Meaning | Status |
|---|---|---|
| $D_{it}$ | binary drug state, 1 on drug and 0 off drug | canonical |
| $D_{bc}$ | continuous exposure-decayed drug indicator (code `Dbc`) | canonical |
| $c_{bm}$ | covariance-moderation correlation (code `c.bm`) | canonical |
| $\beta_{bm}^{\text{lumped}}$ | one-component biomarker slope of paper 06 | canonical |
| E1, E2, E3 | carryover analysis specifications: binary indicator; exposure-weighted $D_{bc}$; binary plus lagged nuisance term | canonical labels |
| `A1`-`A3`, `M1`-`M3` | earlier labels of E1-E3, retained only where this record quotes them | historical |
| D1-D7 | terminology decisions of Section 1 | local |

Table: Symbols and labels used in this plan, with status against NOTATION.md

**Glossary.**

- *Surface form*: the exact spelling or wording of a term as printed
  (for example `Type I error` against `type-I error`).
- *Tier 1, 2, 3*: project-wide mechanical fixes, per-paper consistency
  fixes, and first-use glosses respectively (Sections 2 to 4).
- *First-use expansion*: spelling out an acronym at its first
  reader-facing occurrence in each paper.
- *Off-drug period*: the post-discontinuation span in which carryover
  decays; distinguished from a *washout*, an inserted treatment-free
  clearance gap (decision D5).
- *MCSE*: Monte Carlo standard error.
- *HTE*: heterogeneity of treatment effects.

## Execution status (2026-07-30)

**Executed** with D1 locked to US English. Tier 1 was applied to all 21
`report*.Rmd` files (the plan was written against 13; six variants have
been added since, and paper 11 did not exist). Tier 2, Tier 3, and the
Section 5 full-versus-slim pass were applied where the targets still
existed; most had already converged through intervening revisions.
Section 6 items remain open and are for the author.

What changed, and what was already clean, is recorded at the end of
this file. All 21 manuscripts were re-rendered and built.

**Later changes (note added 2026-10-07).** Two statements below are
superseded. First, the analysis specifications were relabeled again on
2026-07-31, from `M1`/`M2`/`M3` to E1/E2/E3, in the stored data,
drivers and manuscripts (`analysis/report/NOTATION.md`, Part 3); this
plan's references to the specification now use E2, quoting the older
label only where the record requires it. Second, the three-part
`bullets` / `rgt` / `orig` paragraph scaffolding was dropped from all
manuscripts on 2026-10-07 (project `CLAUDE.md`, "Paragraph
structure"), so the remarks below about Lorem ipsum `rgt` placeholders
and a second terminology pass over them no longer apply: the
manuscripts are plain prose.

## Scope and method

Thirteen manuscript files across ten papers were reviewed against the
lexicon (prose, headings, captions, and abstracts only; R code chunks and
in-model variable names such as `c.bm`, `Dbc`, `carryover_t1half` were
excluded). Papers 01, 02, and 06 have both a full and a slim variant; the
variants were reviewed together for cross-file consistency. The `rgt`
paragraph blocks are still Lorem ipsum placeholders and were not
reviewable; a second pass will be required once the author writes that
prose.

**Finding in one line.** No statistical term is *wrong*; the issues are
(i) inconsistent surface forms for the same concept across and within
papers, (ii) acronyms used before first-use expansion, and (iii) a few
colloquial coinages and one genuinely narrow lexicon term ("washout")
used loosely. Three *substantive* (non-terminological) issues surfaced
incidentally and are listed separately in Section 6; they require author
judgment and must not be folded into a mechanical edit.

**Important execution constraint.** These files live under
`~/Library/CloudStorage/` (Dropbox). Per the project guardrails, in-place
stream editors (`sed -i`, `perl -i`, etc.) must NOT be used on
cloud-mounted paths. All edits must go through the Edit tool (which is
path-aware) or the copy-to-`/tmp`, edit, `mv`-back pattern. Global
substitutions should be applied per file with `replace_all`, never with a
recursive `sed -i`.

## 1. Decisions to lock before editing

These choices propagate everywhere, so confirm them first. Recommended
defaults are evidence-based (current corpus counts in parentheses).

| # | Decision | Recommended default | Evidence |
|---|---|---|---|
| D1 | Spelling convention | **US / `-ize`, `-or`, `modeling`** (locked 2026-07-30; overrides the frequency-based recommendation below) | The corpus leans British (modelling 73 vs modeling 43; characterise 26; behaviour 18 vs 3; generalised 16 vs 8), but frequency records drift, not a decision. The project-wide standard is US English, so the direction of travel is British -> US. |
| D2 | Casing of "Type I error" | `Type I error` (cap T, cap I, no hyphen) | Type I error 121 vs type I 61 vs type-I 45 |
| D3 | Biomarker-interaction surface form | `biomarker-treatment interaction` | 91 vs biomarker-by-treatment 45; matches lexicon §10 headword |
| D4 | Monte Carlo SE form | spell `Monte Carlo standard error (MCSE)` at first use, then `MCSE` | MCSE 281, MC SE 59, spelled 57, "Monte Carlo SE" 25 |
| D5 | "washout" usage | reserve `washout` for an inserted treatment-free clearance gap; use `off-drug period / occasions` for the general post-discontinuation span | off-drug 279 vs washout 92; the papers' own thesis is that carryover contaminates these spans, so "washout" is self-contradictory there (lexicon §2) |
| D6 | Precision vs personalized medicine | `precision medicine` (primary), "personalized medicine" only as a glossed synonym | precision 11 vs personalized (corpus form `personalised`) 4; lexicon §10 headword |
| D7 | Participant unit noun | `participant` in prose; use "subject"/"cluster" only where theory requires (e.g. sandwich estimator) | papers 05/08 use participant; 10 mixes cluster/subject/participant |

Table: Terminology decisions with recommended defaults and supporting corpus evidence

## 2. Tier 1: project-wide standardizations (safe, high value)

Apply across all 13 files. Each is a low-risk, presentational fix that
does not alter meaning. Apply with the Edit tool (`replace_all` per file),
checking each hit in context (especially D5, which is contextual, not a
blind replace).

1. **Type I / Type II error casing** (D2). Normalize `type I error`,
   `type-I error`, `type-I`, `Type-I` -> `Type I error`; likewise
   `Type II`. Affects 02, 04, 05, 07, 08, 09, 10. (~106 non-canonical
   instances.)
2. **Biomarker-treatment interaction** (D3). `biomarker-by-treatment
   interaction` -> `biomarker-treatment interaction`; also normalize
   `treatment-by-biomarker`, `biomarker:treatment`, "predictive
   interaction" where they denote the same estimand. Affects 02, 03, 06,
   08, 09, 10. (~45 instances.)
3. **MCSE** (D4). `MC SE` and `Monte Carlo SE` -> `MCSE`; ensure one
   spelled-out first use per file. Affects 01, 02, 07 (and any ADEMP
   paper). (~84 instances.)
4. **HTE first-use + label**. Spell `heterogeneity of treatment effects
   (HTE)` at first use in every paper that uses the concept, then `HTE`.
   Replace loose paraphrases ("treatment-response heterogeneity",
   "individual treatment heterogeneity", "between-person response
   variability") with HTE or with "patient-by-treatment interaction"
   where the variance-component meaning is intended. Critical in 06
   (concept central but acronym never introduced) and 01 (used before
   expansion). Also 03, 04.
5. **"washout" -> "off-drug period/occasions"** (D5). Contextual: change
   only where the text means the post-discontinuation span in which
   carryover decays; keep "washout" where a clearance gap is genuinely
   inserted. Affects 01, 02, 04. (~92 "washout" hits to triage.)
6. **Acronym first-use expansion.** Ensure each is spelled out at first
   reader-facing use, then abbreviated: `DGP` (data-generating
   process/mechanism), `MVN` (multivariate normal), `LME` (linear
   mixed-effects model), `MCAR`/`MAR`/`MNAR`, `OL`/`BDC`/`OL+BDC`, `CO`
   (crossover), `PK/PD`, `BLUP`, `SEM`, `IL-6`, `RM-ANOVA`, `OFAT`. The
   three-part `bullets`/`orig` structure aggravates this because a
   `bullets` block often leads with a bare acronym that the later `orig`
   block expands; fix by expanding in whichever block appears first.
7. **Spelling normalization** (D1, direction reversed on 2026-07-30).
   `generalised` -> `generalized`, `modelling` -> `modeling`,
   `behaviour` -> `behavior`, `parameterise` -> `parameterize`,
   `characterise` -> `characterize`, `normalise` -> `normalize`,
   `whilst` -> `while`, `artefact` -> `artifact`, `grey` -> `gray`,
   etc. Preserve British spelling inside verbatim quotations and in
   published citation or reference titles. Do not touch the US-English
   false friends (`exercise`, `revise`, `raise`, `advise`, `promise`,
   `compromise`, `analyses`, `characteristic`, `programmer`).
8. **Code tokens leaking into prose.** Replace reader-facing code
   identifiers with their prose terms: `OLBDC` -> `OL+BDC` (02),
   `modGompertz` -> "modified Gompertz function" in narrative (07; keep
   `modGompertz` only when referring to the R function). Keep all genuine
   model-symbol references (`$D_{bc}$`, `$c.bm$`, `bm:Dbc`) as code.
   (Note added 2026-10-07: `NOTATION.md` rule 3 now separates the two
   registers, so mathematics uses $D_{bc}$ and $c_{bm}$ while code
   contexts use `Dbc` and `c.bm`; the mixed form `$c.bm$` should not
   be used.)

## 3. Tier 2: per-paper consistency fixes

| Paper | Fix | Notes |
|---|---|---|
| 01 | Standardize the drug-indicator naming ("treatment indicator" / "drug indicator" / "drug exposure variable" -> one prose term for the DGP indicator and one for the analysis variable). | Within-paper inconsistency for `$D_{it}$` / `Dbc`. |
| 02 | Settle on ONE label for analysis spec E2 (`A2` when this plan was written; "matched-decay" / "exposure-weighted" / "continuous-$D_{bc}$"); reserve "matched-decay" for the matched DGP cell. Fix "censored" where it means dropped/missing (lexicon §8 vs §11). Standardize "Hybrid N-of-1 design" and "crossover (CO)" (drop "classical/traditional" alternation). | E2 label multiplicity (~25) is the biggest clarity win. |
| 03 | Disambiguate unobserved-group vocabulary: use `latent class` for unobserved groups; reserve `subgroup`/`stratum` for observed splits (lexicon §10). Gloss `subadditivity` and `frailty` at first use; change "diagnostic biomarker" -> `predictive biomarker` where differential benefit is meant. | "frailty" is an unglossed third synonym for random effect. |
| 06 | Replace colloquial `lumped` -> `one-component` / `single-indicator` (both already used). Define `response expectancy` once, then use "expectancy" rather than rotating belief/suggestion/suggestibility/optimism. Gloss `balanced-placebo (open-hidden) design` at first use. Standardize "between-participant". | See also the Study B/C cross-file bug in Section 6. |
| 07 | Capitalize `Gompertz` consistently in prose (lowercase only as a factor level/code). Remove `modGompertz` from narrative (Tier 1.8). Add a one-line gloss tying Hill / sigmoid-Emax / indirect-response to the lexicon "Emax / dose-response model" cluster. | Mostly cosmetic but pervasive. |
| 08 | Disambiguate `main effect`: use `average treatment effect (ATE)` for the population mean effect; keep "main effect" only for the ANOVA-factor sense. Introduce `anti-conservative (inflated Type I error)` once, then use one term. | "main effect" overload is a genuine ambiguity. |
| 09 | Consolidate `informative dropout` (canonical) over "biased dropout"/"attrition" (keep "biased" only for cell labels). Use `ignorable` (not "non-informative") for the MAR component to avoid collision with "noninformative prior". Replace colloquial `estimand drift` -> "bias in the estimand" and `design surface area` -> a defined phrase. Normalize `MNAR` spelling. | See MCAR_25 flag in Section 6. |
| 10 | Choose one orthography for `mis-specified`/`misspecified`. Introduce `cluster-robust (sandwich) standard error` once, then one term. Keep `working correlation structure` for the GEE object only; use `working covariance` for the LME structure. Fix participant noun (D7). | Most careful paper; fixes are orthographic. |

Table: Per-paper terminology consistency fixes with notes

## 4. Tier 3: glossing and definition additions

Add a one-line first-use gloss (no statistical change) for load-bearing
terms absent from the lexicon: `subadditivity` (03, 01), `frailty` (03),
`balanced-placebo / open-hidden design` (06), `vulnerability window` and
`design surface area` (09), `lasagna plot` (02). Consider adding
`subadditivity`, `balanced-placebo design`, and `gating function` as new
entries to `docs/29-...-lexicon.md` so future papers inherit a fixed
definition.

## 5. Cross-file consistency (full vs slim)

For papers 01, 02, 06, run a final diff pass so the slim variant uses the
same canonical forms and acronym expansions as the full report
(e.g. IL-6 is defined in 01 full but not slim; "RM-ANOVA" vs
"repeated-measures analogue" differs across 02 full/slim).

## 6. Out of scope: substantive issues to flag to the author

These are NOT terminology and must not be auto-edited; each needs author
judgment because a "fix" could change meaning or a stated result.

1. **Paper 02, `report.Rmd` Conclusion #1.** States E2 (then labeled
   `A2`) is "the recommended default ... across the full factorial",
   which conflicts with the paper's own data (E2 leads in ~19% of
   cells) and with the slim variant's correctly hedged conclusion. (Consistent with the prior
   project note that the E2-dominance headline is not universal.)
2. **Paper 06, Study B/C naming divergence.** The full report labels the
   belief-decoupling study "Study B" (sec 7.5) and the coupling pilot the
   "contaminated-biomarker pilot" (7.4); the slim report calls them
   "Study C" and "Study B pilot". Any reader cross-referencing the two
   files is misdirected. Decide one scheme and apply to both.
3. **Paper 09, `MCAR_25` cell label.** Its flat dropout component scales
   with elapsed time (a design variable), which is MAR, not MCAR, under
   the lexicon definitions. Renaming would change a stated result label,
   so confirm the mechanism before any change.

*Status spot-check (2026-10-07, inspected, not exhaustive).* The
quoted "across the full factorial" wording no longer occurs in paper
02's `report.Rmd`; paper 06's `report-slim.Rmd` still uses "Study B
pilot" and "Study C", so item 2 is plausibly still open; item 3 was
not re-checked.

## 7. Execution workflow

1. Lock decisions D1-D7 (Section 1).
2. Apply Tier 1 globally, one file at a time, via the Edit tool with
   `replace_all`; triage D5 ("washout") hit-by-hit. Do not use `sed -i`
   on these cloud paths.
3. Apply Tier 2 per paper; Tier 3 glosses; Section 5 full-vs-slim diff.
4. Re-render each changed manuscript with the appropriate R environment
   (note: paper 04 requires the Framework R / kableExtra environment per
   the project render notes) and confirm it still builds and that the
   PDF render-timestamp footer is updated (since superseded: the render
   wrapper `tools/render.sh` now stamps provenance on every page).
5. Address Section 6 items separately with the author.

**Effort estimate.** Tier 1 is mechanical and dominated by D2/D3/D4
(~235 substitutions) plus contextual D5 triage (~92 hits). Tier 2 is a
few dozen targeted edits per paper. The re-render and full-vs-slim diff
are the main time cost. None of the work changes a number, estimand, or
model.

## 8. Execution record (2026-07-30)

**Tier 1 applied** across 21 files by a chunk-aware transform (prose
rules skip fenced R chunks so that ggplot's `colour=` and dplyr's
`summarise()` are untouched; display-string rules such as `MC SE` are
applied inside chunks as well, because those strings are table
headers):

- D2 Type I/II casing, D3 `biomarker-by-treatment` ->
  `biomarker-treatment`, D4 `MC SE`/`MC SEs` -> `MCSE`/`MCSEs`,
  Tier 1.8 `OLBDC` -> `OL+BDC`.
- D1 in the US direction: `modelling`, `generalis*`, `behaviour`,
  `characteris*` (guarding `characteristic`), `normalis*`,
  `summaris*` (prose only), `artefact`, `randomis*`, `standardis*`,
  `dichotomis*`, `ageing`, `sceptic*`, `favour*`, `enrolment`,
  `signalling`, `programme`, `labelled`, `modelled`, `analyse`,
  `centre*`, `colour` (prose only), `whilst`, `amongst`.
- 16 of 21 files changed; residual British-looking tokens are all R
  identifiers (`dplyr::summarise`) and are correct as they stand.

**Tier 1 items already satisfied:** HTE is expanded at first use in
01 and 06; the acronym sweep (DGP, MVN, LME, MCAR/MAR/MNAR, OL/BDC,
CO, PK/PD, BLUP, SEM, RM-ANOVA, OFAT, GEE, ICC) found only five true
misses, fixed here: MCAR and MAR in 02, BLUP in 04, DGP and MVN in 11,
plus RM-ANOVA in 01-slim and MCAR/MAR in 02-slim for Section 5
parity. ADEMP is treated as a named, always-cited framework and is
not expanded.

**D5 (`washout`) triage.** 35 hits reviewed individually. All but
three are genuine inserted clearance gaps and were kept. The three
changed are in paper 02's toy example, where the span is the
carryover-contaminated post-discontinuation window: `three-week
washout` -> `three-week off-drug period`, `the second washout week` ->
`the second off-drug week`, `throughout the washout` -> `throughout
the off-drug period`.

**Tier 2, applied:** 02, one residual `matched-decay (M2)` reference to
the specification changed to `exposure-weighted (M2)` (the
matched-decay *cell* retains the term, as the plan specifies).

**Tier 2, already converged:** 06 (`lumped` survives only as the math
superscript $\beta_{bm}^{\text{lumped}}$, not as prose); 07
(`modGompertz` absent from narrative; lowercase `gompertz` occurs only
in citation keys, file paths, and equation labels); 08 (`main effect`
absent; `anti-conservative` already single-form); 09 (`attrition`,
`estimand drift`, `design surface area`, `non-informative` all absent;
`biased dropout` survives only as cell-label description, which the
plan permits); 10 (`misspecified` already uniform, 13/13).

**Tier 3 glosses added:** `subadditivity` (03), `frailty` (03),
`balanced-placebo design` (06). `lasagna plot` has no occurrences.

**Not done:** the plan's suggestion to add `subadditivity`,
`balanced-placebo design`, and `gating function` to
`docs/29-...-lexicon.md`; Section 6 substantive items; and the second
terminology pass over the `rgt` blocks, which remain placeholders.

**Interaction with the notation pass.** This execution followed the
symbol-level pass recorded in `analysis/report/NOTATION.md`. The
analysis specifications were then `M1`/`M2`/`M3` in all reader-facing
text, which is why this record says M2 where the 2026-06-17 draft said
`A2`. They have since been relabeled E1/E2/E3 (2026-07-31; see the
note under the execution status), which is the label used in Sections
3 and 6 above.
