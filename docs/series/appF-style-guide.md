# Style guide

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team
**Purpose:** Prose, notation, and figure conventions for the series
and for manuscripts.
**Prerequisites:** None

---

## Contents

1. [Notation](#notation)
2. [Summary](#summary)
3. [Construct box](#construct-box)
4. [Setting and argument](#setting-and-argument)
5. [Derivations](#derivations)
6. [Numerical evidence](#numerical-evidence)
7. [Relation to the manuscripts](#relation-to-the-manuscripts)
8. [Limitations and open items](#limitations-and-open-items)
9. [Reproduce](#reproduce)
10. [References](#references)
11. [Change log](#change-log)

---

## Notation

*Local symbols used in this chapter only. For project-wide notation
see `analysis/report/NOTATION.md`.*

| Symbol or label | Meaning | Status |
|---|---|---|
| $D_{it}$ | binary drug state, 1 on drug and 0 off drug | canonical |
| $D_{bc}$ | continuous exposure-decayed drug indicator (`Dbc`) | canonical |
| $c_{bm}$ | covariance-moderation correlation (`c.bm`) | canonical |
| $\beta_{bm}^{\text{lumped}}$ | one-component biomarker slope, paper 06 | canonical |
| E1, E2, E3 | exposure codings: binary indicator; `Dbc`; binary plus lagged nuisance | canonical labels |
| D1-D7 | terminology decisions, Section 1 of `docs/30` | local |

Table: Symbols and labels used in this appendix

---

## Summary

This appendix combines two sources:

1. **Terminology consistency plan** (`docs/30`): decisions D1-D7,
   Tier 1-3 standardizations, cross-file consistency rules, and the
   execution record (2026-07-30). Governs prose surface forms for the
   manuscripts.
2. **Whitepaper series notation rules** (`docs/59` Sections 7.1-7.3):
   notation fixes when carrying content from pre-34 documents, known
   regressions in docs/34-58, and nine cross-cutting rules for all
   series chapters.

The two sources govern different things: `docs/30` governs manuscript
prose; `docs/59` Section 7 governs notation and presentation of
derived content in series chapters. Both apply when writing chapters
that quote or summarize manuscript results.

---

## Construct box

*No Monte Carlo simulation results in this chapter.*

---

## Setting and argument

### Part I: Terminology consistency plan

<!-- From docs/30-terminology-consistency-plan.md (2026-06-17;
execution record added 2026-07-30) -->

#### Scope and method

Thirteen manuscript files across ten papers were reviewed against the
lexicon (prose, headings, captions, and abstracts only; R code chunks
and in-model variable names such as `c.bm`, `Dbc`, `carryover_t1half`
were excluded). Papers 01, 02, and 06 have both a full and a slim
variant; the variants were reviewed together for cross-file
consistency.

**Finding in one line.** No statistical term is wrong; the issues are
(i) inconsistent surface forms for the same concept across and within
papers, (ii) acronyms used before first-use expansion, and (iii) a
few colloquial coinages and one genuinely narrow lexicon term
("washout") used loosely. Three substantive (non-terminological)
issues surfaced incidentally and are listed separately in Section 6 of
`docs/30`; they require author judgment and must not be folded into a
mechanical edit.

**Execution status (2026-07-30).** Tier 1 was applied to all 21
`report*.Rmd` files. The analysis specifications were relabeled from
`M1`/`M2`/`M3` to E1/E2/E3 on 2026-07-31. The three-part
`bullets`/`rgt`/`orig` paragraph scaffolding was dropped from all
manuscripts on 2026-10-07; manuscripts are now plain prose and D5
triage of `rgt` blocks no longer applies.

#### 1. Decisions locked before editing (D1-D7)

| # | Decision | Locked default | Evidence |
|---|---|---|---|
| D1 | Spelling convention | **US English** (`-ize`, `-or`, `modeling`); British -> US direction; locked 2026-07-30 | Project-wide standard; overrides corpus-frequency recommendation |
| D2 | Casing of "Type I error" | `Type I error` (cap T, cap I, no hyphen) | `Type I error` 121 vs `type I` 61 vs `type-I` 45 |
| D3 | Biomarker-interaction surface form | `biomarker-treatment interaction` | 91 vs `biomarker-by-treatment` 45; matches lexicon Section 10 headword |
| D4 | Monte Carlo SE form | Spell `Monte Carlo standard error (MCSE)` at first use, then `MCSE` | MCSE 281; MC SE 59; spelled 57 |
| D5 | "washout" usage | Reserve `washout` for an inserted treatment-free clearance gap; use `off-drug period / occasions` for the general post-discontinuation span | off-drug 279 vs washout 92; "washout" self-contradictory in carryover context (Appendix E Section 2) |
| D6 | Precision vs personalized medicine | `precision medicine` (primary); "personalized medicine" only as a glossed synonym | precision 11 vs personalized (British `personalised`) 4; Appendix E Section 10 headword |
| D7 | Participant unit noun | `participant` in prose; use "subject"/"cluster" only where theory requires (e.g., sandwich estimator) | papers 05/08 use participant; 10 mixes cluster/subject/participant |

Table: Terminology decisions with locked defaults and supporting
corpus evidence (source: `docs/30` Section 1)

#### 2. Tier 1: project-wide standardizations (safe, high value)

Apply across all manuscript files. Each is a low-risk, presentational
fix that does not alter meaning. Apply with the Edit tool
(`replace_all` per file), checking each hit in context (especially D5,
which is contextual, not a blind replace).

1. **Type I / Type II error casing** (D2). Normalize `type I error`,
   `type-I error`, `type-I`, `Type-I` to `Type I error`; likewise
   `Type II`. Affects 02, 04, 05, 07, 08, 09, 10.
2. **Biomarker-treatment interaction** (D3). `biomarker-by-treatment
   interaction` to `biomarker-treatment interaction`; normalize
   `treatment-by-biomarker`, `biomarker:treatment`, "predictive
   interaction" where they denote the same estimand. Affects 02, 03,
   06, 08, 09, 10.
3. **MCSE** (D4). `MC SE` and `Monte Carlo SE` to `MCSE`; ensure one
   spelled-out first use per file. Affects 01, 02, 07.
4. **HTE first-use and label.** Spell `heterogeneity of treatment
   effects (HTE)` at first use in every paper that uses the concept,
   then `HTE`. Replace loose paraphrases with HTE or with
   "patient-by-treatment interaction" where the variance-component
   meaning is intended.
5. **"washout" to "off-drug period/occasions"** (D5). Contextual:
   change only where the text means the post-discontinuation span in
   which carryover decays; keep "washout" where a clearance gap is
   genuinely inserted. Affects 01, 02, 04.
6. **Acronym first-use expansion.** Ensure each is spelled out at
   first reader-facing use: `DGP`, `MVN`, `LME`, `MCAR`/`MAR`/`MNAR`,
   `OL`/`BDC`/`OL+BDC`, `CO` (crossover), `PK/PD`, `BLUP`, `SEM`,
   `IL-6`, `RM-ANOVA`, `OFAT`.
7. **Spelling normalization** (D1). `generalised` to `generalized`,
   `modelling` to `modeling`, `behaviour` to `behavior`,
   `parameterise` to `parameterize`, `characterise` to `characterize`,
   `normalise` to `normalize`, `whilst` to `while`, `artefact` to
   `artifact`, `grey` to `gray`, etc. Preserve British spelling inside
   verbatim quotations and in published citation or reference titles.
   Do not touch US-English false friends (`exercise`, `revise`,
   `raise`, `advise`, `promise`, `compromise`, `analyses`,
   `characteristic`, `programmer`).
8. **Code tokens in prose.** Replace reader-facing code identifiers:
   `OLBDC` to `OL+BDC`; `modGompertz` to "modified Gompertz
   function" in narrative (keep `modGompertz` only when referring to
   the R function). Keep genuine model-symbol references (`$D_{bc}$`,
   `$c_{bm}$`, `bm:Dbc`) as code. Do not use the mixed form `$c.bm$`
   (NOTATION.md rule 3).

#### 3. Tier 2: per-paper consistency fixes

| Paper | Fix | Notes |
|---|---|---|
| 01 | Standardize the drug-indicator naming ("treatment indicator" / "drug indicator" / "drug exposure variable" to one prose term for the DGP indicator and one for the analysis variable) | Within-paper inconsistency for `$D_{it}$` / `Dbc` |
| 02 | Settle on ONE label for analysis spec E2; reserve "matched-decay" for the matched DGP cell only. Fix "censored" where it means dropped/missing. Standardize "Hybrid N-of-1 design" and "crossover (CO)" | E2 label multiplicity the biggest clarity win |
| 03 | Disambiguate unobserved-group vocabulary: use `latent class` for unobserved groups; reserve `subgroup`/`stratum` for observed splits. Gloss `subadditivity` and `frailty` at first use; change "diagnostic biomarker" to `predictive biomarker` where differential benefit is meant | "frailty" is an unglossed synonym for random effect |
| 06 | Replace colloquial `lumped` to `one-component` / `single-indicator`. Define `response expectancy` once, then use "expectancy". Gloss `balanced-placebo (open-hidden) design`. Standardize "between-participant" | Study B/C naming divergence (item 2 of Section 6 in `docs/30`) still plausibly open as of 2026-10-07 |
| 07 | Capitalize `Gompertz` consistently in prose. Remove `modGompertz` from narrative. Add one-line gloss tying Hill / sigmoid-Emax / indirect-response to the lexicon | Mostly cosmetic but pervasive |
| 08 | Disambiguate `main effect`: use `average treatment effect (ATE)` for the population mean effect; keep "main effect" only for the ANOVA-factor sense. Introduce `anti-conservative (inflated Type I error)` once | "main effect" overload is a genuine ambiguity |
| 09 | Consolidate `informative dropout` over "biased dropout"/"attrition". Use `ignorable` (not "non-informative") for the MAR component. Replace colloquial `estimand drift` to "bias in the estimand". Normalize `MNAR` spelling | See `MCAR_25` flag in `docs/30` Section 6 |
| 10 | Choose one orthography for `mis-specified`/`misspecified`. Introduce `cluster-robust (sandwich) standard error` once. Keep `working correlation structure` for the GEE object only; use `working covariance` for the LME structure. Fix participant noun (D7) | Most careful paper; fixes are orthographic |

Table: Per-paper terminology consistency fixes (source: `docs/30`
Section 3)

#### 4. Tier 3: glossing and definition additions

Add a one-line first-use gloss for load-bearing terms absent from
the lexicon: `subadditivity` (03, 01), `frailty` (03),
`balanced-placebo / open-hidden design` (06), `vulnerability window`
and `design surface area` (09), `lasagna plot` (02). Consider adding
`subadditivity`, `balanced-placebo design`, and `gating function` as
new entries to Appendix E so future papers inherit a fixed definition.

#### 5. Cross-file consistency (full vs slim)

For papers 01, 02, 06, run a final diff pass so the slim variant uses
the same canonical forms and acronym expansions as the full report
(e.g., `IL-6` is defined in 01 full but not slim; "RM-ANOVA" vs
"repeated-measures analogue" differs across 02 full/slim).

#### 6. Out of scope: substantive issues to flag to the author

These are NOT terminology and must not be auto-edited; each needs
author judgment because a "fix" could change meaning or a stated
result.

1. **Paper 02, `report.Rmd` Conclusion no. 1.** States E2 is "the
   recommended default ... across the full factorial", which
   conflicts with the paper's own data (E2 leads in approximately 19%
   of cells) and with the slim variant's correctly hedged conclusion.
   Status (2026-10-07, inspected): the quoted "across the full
   factorial" wording no longer occurs in paper 02's `report.Rmd`.
2. **Paper 06, Study B/C naming divergence.** The full report labels
   the belief-decoupling study "Study B" and the coupling pilot the
   "contaminated-biomarker pilot"; the slim report calls them "Study C"
   and "Study B pilot". Status (2026-10-07, inspected): paper 06's
   `report-slim.Rmd` still uses "Study B pilot" and "Study C";
   plausibly still open.
3. **Paper 09, `MCAR_25` cell label.** Its flat dropout component
   scales with elapsed time (a design variable), which is MAR, not
   MCAR, under the lexicon definitions. Status: not re-checked
   2026-10-07.

#### 7. Execution workflow

1. Lock decisions D1-D7 (Section 1).
2. Apply Tier 1 globally, one file at a time, via the Edit tool with
   `replace_all`; triage D5 ("washout") hit-by-hit. Do not use
   `sed -i` on cloud paths.
3. Apply Tier 2 per paper; Tier 3 glosses; Section 5 full-vs-slim
   diff.
4. Re-render each changed manuscript with `tools/render.sh` (not
   `rmarkdown::render()` directly).
5. Address Section 6 items separately with the author.

#### 8. Execution record (2026-07-30)

**Tier 1 applied** across 21 files:

- D2 Type I/II casing, D3 `biomarker-by-treatment` to
  `biomarker-treatment`, D4 `MC SE`/`MC SEs` to `MCSE`/`MCSEs`,
  Tier 1.8 `OLBDC` to `OL+BDC`.
- D1 in the US direction: `modelling`, `generalis*`, `behaviour`,
  `characteris*` (guarding `characteristic`), `normalis*`,
  `summaris*` (prose only), `artefact`, `randomis*`, `standardis*`,
  `dichotomis*`, `ageing`, `sceptic*`, `favour*`, `enrolment`,
  `signalling`, `programme`, `labelled`, `modelled`, `analyse`,
  `centre*`, `colour` (prose only), `whilst`, `amongst`.
- 16 of 21 files changed; residual British-looking tokens are all R
  identifiers (`dplyr::summarise`) and are correct as they stand.

**D5 (`washout`) triage.** 35 hits reviewed individually. Three
changed: in paper 02's toy example, `three-week washout` to
`three-week off-drug period`, `the second washout week` to `the
second off-drug week`, `throughout the washout` to `throughout the
off-drug period`.

**Not done (still open):** adding `subadditivity`, `balanced-placebo
design`, and `gating function` to Appendix E; Section 6 substantive
items; the second terminology pass over `rgt` blocks (now moot:
the scaffolding was dropped 2026-10-07).

---

### Part II: Whitepaper series notation rules

<!-- From docs/59 Sections 7.1, 7.2, 7.3 -->

These rules apply to all series chapters (front matter, parts, and
appendices), not just manuscripts.

#### 7.1 Documents before docs/34: notation fixes on carry-in

These documents are archived rather than edited, so the fixes apply
when their content is carried into a chapter or into Appendix I.

| Fix | Where it occurs | Canonical form |
|---|---|---|
| Architecture letters | all `.tex`; 21, 23, 32, `plan-architecture-c.md` | mean-moderation, covariance-moderation, dual-channel |
| "N-of-1" for the Hybrid design | 13, 14, 18 | Hybrid |
| ER, TR, "pharm_biomarker" for components | 07, 07a, 09, 13, 15, 18 | $PB$, $TV$ |
| $S_{i,t}$, `Sx` as a symbol | `.tex` set | $Y_{it}$ |
| $\rho$ for compound symmetry | 07a, 18, 23 Sections 2.1-3.3, 31 Section 20, 32 | $\rho_{cs}$ |
| $\rho$ or $\rho_t$ for the biomarker-response correlation | 04 Section 3.1, 08, 09 | $c_{bm}\,g_t$ |
| $c_{bm}$ as the mean-moderation multiplier | 01, 19.tex | $\beta_{bm}$ |
| $B_i^{\ast}$ for the standardized biomarker | 01, 19.tex | $b_i$ |
| $D_{bc}$ for a binary indicator | 06 Section 4.1, 09 | $D_{it}$ |
| $e_t$ for expectancy | 07 | $\eta_t$ |
| $\alpha$ for the decay factor or a power-law exponent; $\lambda_{cor}$, $\lambda_{drug}$ | 08, 09 | $\phi_t$; $\lambda$ |
| $\beta_t$, $\beta_1$, $\beta_2$ for time and biomarker coefficients | 23 Section 4.1, 24, 25, 31 | $\beta_T$, $\beta_B$ |
| "biomarker-by-treatment", "washout" for the off-drug span, lower-case "hybrid" | 21, 23, 24, 28, 31 | biomarker-treatment interaction; off-drug period; Hybrid |
| hand-written "Rendered on" footers | every `.tex` | removed; `tools/render.sh` owns provenance |
| British spellings | `reorg-plan-2026-04-22.md` | US English |

Table: Notation and convention fixes for content carried from
documents before docs/34 (source: `docs/59` Section 7.1)

#### 7.2 Known notation regressions in docs/34-58

| Document and location | Problem | Fix |
|---|---|---|
| 34 Sections 3.1-3.2 | $\tilde r$; unpadded $g$ in the ceiling formula; $C$ as cross block | $\tilde g$; $(\tilde g^\top M^{-1}\tilde g)^{-1/2}$; name block explicitly |
| 35 | $q_K$, $q_A$ (eigenvectors); $\Delta$ as a time lag; $\mu$ as an eigenvalue; $v_{TV}$, $v_{PB}$, $v_{BR}$; $s_t$ as factor SDs | Letters chosen against NOTATION's reserved list; fix each to a non-reserved letter |
| 37 Sections 4, 4.4 | plain `$E(\cdot)$`; $\mu$, $\mu_p$, $\bar\mu$ local; $v_t$ factor SDs; $e_i$ residual; "patient(s)" | $\mathbb{E}(\cdot)$; another letter; participant |
| 38 | $f(t)$, $\bar f$, $F$; $G$, $H$ in Section 6 proof; $\sigma_v^2$ | $h(t)$, $\Delta_h$ (as in docs/44); other letters; $\tau^2$ (DQ-6) |
| 41 line 205 | `$N(\cdot)$`; $f(t)$; $d_{it}$ | $\mathcal{N}(\cdot)$; another letter; another letter |
| 42, 45 | $f(t)$ for the time adjustment | another letter (e.g., $s(w_t)$) |
| 45 index | $B_i$ "centered" | $B^c_i$ |
| 51 | $\sigma_v^2$; `$E(t)$`, `$E_{it}$` beside E1-E9; $L$, $U$ beside $L_{it}$ | $\tau^2$ (DQ-6); other letters; rename |
| 56 Section 1 | AR(1) exponent on the visit index | $\lvert w_t - w_s\rvert$ |
| 57 | $M$, $m$ (bootstrap); $g(Y)$; $\gamma$; $h$ (candidate half-life); $\omega$; $\hat\xi$; $e_i$; grid $\{0, 0.25, 0.5, 1, 2\}$ | letters not on reserved list; one half-life grid across docs/57 and docs/58 (DQ-10) |
| 58 | $B$ (GLS bread), $G$ (randomization distribution) beside Gompertz $G$; analysis labels shared with interaction track | other letters; `M-` labels (Appendix A, Section 2.7) |

Table: Known notation regressions in docs/34-58 (source: `docs/59`
Section 7.2)

#### 7.3 Cross-cutting rules for all series chapters

1. NOTATION is the only symbol authority; a chapter adds only local
   symbols, never on a reserved letter.
2. Every result states its construct in the chapter's construct box;
   a ceiling, size, or power number without $\rho$, the construct,
   carryover, coupling, allocation, and $N$ is not quoted.
3. $N$ is the total; $t_{1/2}$ in weeks; the fixed analysis half-life
   is 0.2 week unless a historical result is quoted.
4. Evidence labels from `docs/59` Section 2.5 only (verified,
   inspected, inferred, unverified).
5. Analysis IDs from the registry only (Appendix A).
6. Each table has a caption line; each figure has a caption with its
   generating script; figure files named in print order.
7. No dated correction paragraphs at release.
8. US English; no em dashes; prose vocabulary per this appendix; this
   project's studies are exploratory and its protocols are working
   protocols, never described with the vocabulary of registered trials.
9. The linter `tools/notation-lint.pl` runs over every chapter; two
   additional checks for the series: plain `$E($` and `$N($` for the
   normal law; and reserved letters reused in a chapter's local table
   (GP-37).

---

## Derivations

<!-- one row per result: statement | CM label | verification
script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

This appendix is the single reference for prose and notation
conventions. Manuscripts cite D1-D7 by number for locked decisions.
Series chapters cite rules 1-9 by number.

---

## Limitations and open items

- Section 6 substantive issues (paper 02 conclusion, paper 06 Study
  B/C naming, paper 09 `MCAR_25` cell) were not resolved in the Tier
  1-3 pass.
- Tier 3 additions to Appendix E (subadditivity, balanced-placebo
  design, gating function) are not yet written.
- docs/59 Section 7.2 regression entries for docs/36 Sections 6.7
  and A.7.5, docs/39, docs/44 index, docs/48 Section 3 note, docs/49
  Section 3 note, docs/49 and 56 `n_on`/`n_off`/`n_v`, docs/53, and
  docs/54 have been omitted from the compact table above; consult
  `docs/59` Section 7.2 for the full list.

---

## Reproduce

The source documents are:

```bash
cat docs/30-terminology-consistency-plan.md
```

```bash
grep -n "^### 7\." docs/59-whitepaper-refactoring-plan.md
```

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | Setting and argument | Content move from docs/30 (all sections) and docs/59 Sections 7.1-7.3 | Phase 2 content move |
