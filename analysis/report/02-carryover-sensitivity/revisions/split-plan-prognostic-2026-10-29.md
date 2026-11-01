# Split plan: moving the prognostic material and the F-family out of paper 02

*2026-10-29 19:16 PDT*

**Author.** pmsimstats team

**Status.** Plan only. No manuscript, driver or data file has been
changed. Line numbers refer to `report.Rmd` at commit 8e810ca9d
(2026-10-29). They will drift as soon as the first cut is made, so
work from the section labels.

## Notation

Symbols follow `analysis/report/NOTATION.md`. Code names are in
`monospace`. Two symbols have no NOTATION.md entry yet. They are
marked *proposed* and should be added when the split is adopted.
Paper 02's current text uses retired forms ($c_{pb}$, $c_{tv}$,
$c_{bl}$, $D_e$). The cut should rename them to the forms below.
Do not declare them as departures.

| Symbol | Meaning | Code |
|---|---|---|
| $B_i$, $B^c_i$ | biomarker of participant $i$; centered at the trial's baseline sample mean | `bm`, `bmc` |
| $w_t$ | calendar week of visit $t$ | `t` |
| $c_{bm}$ | biomarker-response correlation on drug: the interaction under covariance moderation | `c_bm`, `c.bm` |
| $c_{bm}^{\ast}$ | ceiling: the largest $c_{bm}$ for which every path's correlation matrix is positive definite | |
| $c_{bm,PB}$, $c_{bm,TV}$, $c_{bm,BL}$ | prognostic correlations: of the biomarker with the placebo-response component, the natural-history component, and the baseline level | `c_pb`, `c_tv`, `c_bl`; `c.bm.pb` |
| $b_{TV}$ (*proposed*) | constant natural-history mean shift per SD of the biomarker; absorbed by the biomarker main effect | `b_tv`, `b.bm.tv` |
| $b_{TV}^{\mathrm{tr}}$ (*proposed*) | trend-pathway coefficient: shift of $b_{TV}^{\mathrm{tr}}\,\sigma_{TV}\,z_i\,g(w_t)$, with $z_i$ the standardized biomarker and $g$ the TV mean course scaled to peak 1 | `b_tv_trend`, `b.bm.tv.trend` |
| $c_1$ | cross-component correlation at the same visit (refined process 0.2) | `c.cf1t`, `c.cfct` |
| $t_{1/2}$, $t_{1/2}^a$ | carryover half-life of the DGP; half-life assumed by the analysis | `t1half`; `h`, `hl_analysis` |
| $D_{bc,it}$ | exposure-weighted drug indicator built at $t_{1/2}^a$ | `Dbc` |
| $\eta_t$ | design expectancy: 1 open label, 0.5 blinded, 0 at baseline | `e`; regressor `De` |
| $\beta_\eta$, $\beta_{bm:\eta}$ | coefficients of the expectancy terms $\eta_t$ and $B^c_i\eta_t$ | `De`, `bmc:De` |
| $\kappa$ | mean model-based SE divided by the empirical SD of the estimates | `kappa` |
| $\hat\pi$, MCSE | estimated rejection rate; its Monte Carlo SE, $\sqrt{\hat\pi(1-\hat\pi)/n_{\text{reps}}}$ | `reject`, `mcse` |

Paper 02's `bm` is already centered by the DGP (mean 0), so `bm` and
`bmc` coincide in its drivers. The remedy term is written `bmc:t`
($B^c_i w_t$) throughout this plan, and its quadratic variant adds
$B^c_i w_t^2$.

## Glossary

- **Covariance moderation.** The data-generating convention in which
  the interaction is a correlation between the biomarker and the drug
  response, present on drug and decaying off drug. Paper 02 calls it
  the "correlation encoding"; code labels it `mvn`.
- **Refined process.** The recommended DGP: graded coupling,
  separable AR(1) with $\rho = 0.7$ and $c_1 = 0.2$, anchored
  carryover, and a stop on non-positive-definite matrices. Paper 02
  calls it the "extended DGP". **Published process:** Hendrickson's
  code at commit 58b32a9.
- **Predictive biomarker.** A biomarker that changes the drug
  response ($c_{bm} \ne 0$). **Prognostic biomarker:** one associated
  with the outcome through a channel other than the drug response.
  A biomarker can be both.
- **Prognostic pathway.** One channel of prognostic association:
  - **prog-BL:** $c_{bm,BL}$, baseline level.
  - **prog-TV-const:** $c_{bm,TV}$, natural-history component.
  - **prog-PB:** $c_{bm,PB}$, placebo-response component.
  - **constant shift:** $b_{TV}$, stored as `prog-TV-grow` despite
    being constant.
  - **trend pathway:** $b_{TV}^{\mathrm{tr}}$, a natural-history mean
    shift that grows over the trial.

  The first three are correlations and are bounded by feasibility.
  The last two are mean shifts and are not.
- **Feasibility bound.** The largest prognostic correlation, given
  $c_{bm}$, for which the joint correlation matrix is positive
  definite on every path. It has a closed form via the Schur
  complement (`docs/65`).
- **Silent repair; repair artifact.** Before 2026-10-26 the drivers
  replaced a non-positive-definite matrix with a nearby
  positive-definite one without warning. In the strength-0.6
  prognostic cells this attenuated the interaction and showed up as a
  power loss in every analysis (`docs/65`). Repair is now refused by
  default.
- **G-family.** The analysis ladder of paper 02, stored with `E`
  codes. The members this plan names:
  - **G1** (`E1`): binary drug indicator, `corCAR1`.
  - **G6** (`E1cr2`): G1 with CR2.
  - **G8** (`E2cr2`): $D_{bc}$ at the true half-life, with CR2.
  - **G9** (`E7cr2`): $D_{bc}$ at a half-life chosen by AIC on the
    interaction model, with CR2. Slightly liberal (about +0.005);
    the side analysis after the split.
  - **G10** (`E7null_cr2`): half-life chosen by AIC on the model
    without the interaction, then the interaction model refitted
    and tested with CR2. **The recommended analysis of paper 02**
    (Section 2.8). It was called "G9-null" in earlier drafts and the
    referee-variant runs (decision 5).
  - **Oracle benchmark.** An analysis given the true half-life (G3,
    G8). It is not available to a trialist. It serves only as the
    reference that prices the selection step (decision 4).
  - **G8 at a fixed half-life** (`E3fix_cr2`): a prespecified value
    of 0.25 or 1 week. **Dropped from paper 02** (decision 4).
- **E0.** Hendrickson's original analysis: `lmer`, random intercept,
  binary indicator, and the expectancy main effect where expectancy
  varies. It stays in paper 02.
- **F-family.** The G3 formula plus $\eta_t$ and $B^c_i\eta_t$, with
  CR2. Members:
  - **F1:** $t_{1/2}^a$ fixed at the 0.2-week policy value.
  - **F1A:** chosen by AIC with the interaction in the model.
  - **F1A0:** chosen by AIC without it.
  - **F1A0B:** F1A0 with a parametric-bootstrap calibration that
    repeats the selection.

  The family moves to paper 15.
- **Joint recommendation.** The single analysis advice that papers
  02 and 15 give together (Section 4).
- **CR2.** The Bell-McCaffrey bias-reduced cluster-robust standard
  error, used with Satterthwaite degrees of freedom.
- **Scenario family.** A set of simulation cells sharing one
  departure from paper 02's original grid (Section 2.3).
- **Block (N, A, T, S).** The four blocks of
  `36-run-referee-variants.R`:
  - **N:** null cells.
  - **A:** alternative cells.
  - **T:** trend pathway.
  - **S:** patient slope variance.
- **Patient slope variance.** A patient-specific linear time trend
  unrelated to the biomarker (`sd_slope`). It is a heterogeneity
  stress, not a prognostic one, and stays in paper 02.
- **Two-MCSE band.** The interval of 0.05 plus or minus two MCSE:
  $[0.036, 0.064]$ at 1,000 replicates.
- **Hybrid, OL+BDC, CO.** The designs: open-label titration then
  blinded discontinuation then a brief crossover; open-label
  titration then blinded discontinuation; two-period crossover.

## 1. Decision

Paper 02 is narrowed to one question: how to represent carryover and
choose the analysis half-life for the biomarker-treatment interaction
when the biomarker is predictive only. The data come from the refined
process under covariance moderation, meaning the biomarker-response
correlation is the interaction mechanism.

Everything about a biomarker that is also prognostic moves to a new
companion paper, provisionally paper 15. So does the F-family (F1,
F1A, F1A0, F1A0B), whose expectancy terms exist to protect against a
prognostic biomarker.

**Recommended analysis (decided 2026-10-30).** Paper 02 recommends
**G10**:

- $D_{bc}$, with the analysis half-life chosen by AIC over
  $\{0.25, 0.5, 1, 2\}$ weeks on the model *without* the interaction
- the interaction model refitted by REML at that half-life
- the test of `bm:Dbc` by CR2 with Satterthwaite df

G9, which selects on the model *with* the interaction, becomes the
comparator that shows what that choice costs. The reason is size.
In the paired block N evidence, G9 is liberal against the oracle in
every cell, while G10 cannot be distinguished from it (Section
2.8). G10 gives up 1 to 3 points of power for that. Among
analyses that do not need the true half-life, it is the most powerful
one that holds nominal size. That is the sense in which it is "the
best model" in paper 02. It is not the most powerful analysis, and
the paper must say so. Section 2.8 lists what has to be run before
the paper can carry this claim.

A note on wording: the narrowed paper 02 does not use the published
process at 58b32a9. Step coupling in that process removes the signal
at any positive half-life (`docs/55`). Describe the scope as "the
refined process under covariance moderation, with a predictive-only
biomarker".

**Rationale.**

- The split removes two of paper 02's submission blockers outright:
  - the repaired strength-0.6 rows
  - the untested `bmc:t` remedy

  Both become properly scoped results in paper 15, where the 0.6
  cells are the repair-artifact exhibit and the `bmc:t` check becomes a
  tested analysis.
- It cuts about 15 to 20% of the manuscript.
- It gives the prognostic material a coherent paper of its own. That
  paper is largely drafted already across `docs/47`, `docs/49`,
  `docs/56`, `docs/64`, `docs/65` and `docs/series/ch19`.

## 1a. Decisions needed before work starts

Five points emerged while mapping the split. All five are resolved;
the detail is in the section cited.

1. **Which F-family member carries the joint recommendation
   (Section 4). Resolved 2026-10-30.** With G10 as paper 02's
   analysis, the expectancy-term counterpart is F1A0. That is
   G10's selection rule plus $\eta_t$ and $B^c_i\eta_t$, in its
   bootstrap-calibrated form F1A0B. This agrees with CLAUDE.md's
   provisional recommendation.

   F1, which fixes the half-life at 0.2 week, is not the
   recommendation. In paper 02's cells it loses up to 10 points of
   power to a misspecified half-life.
2. **Large-N size evidence (Section 2.6). Resolved 2026-10-30:
   option (a).** Family 6, the $N \in \{140, 280\}$ cells, is entirely
   prog-PB. After the split, paper 02 has no size evidence above
   $N = 70$. The paper states the $N \le 70$ limit, and no large-N run
   is made.
3. **Overlap between paper 01 Study 3 and paper 15 (Section 5.2).
   Resolved 2026-10-31: cite, do not reopen.** Paper 01 already
   contains a prognostic stress test and a joint feasibility table
   (`tab:prognostic`). Paper 15 cites paper 01 for the DGP-level
   facts and presents only the analysis-level consequences. Paper 01,
   which is near submission at grade B, is not reopened for the
   prognostic material. **Revised 2026-10-31: paper 01 is reopened
   narrowly for G9 to G10.** Paper 01 simulates G9 in its planning
   tables (`tab:planning`, `tab:planningg9`) and its Figure 4B panel
   (`fig:fig4bg9`) and calls G9 the analysis the companion paper
   recommends. Those cells are rerun with G10 on the same trials
   (Section 5.1) and the text is changed to match; nothing else in
   paper 01 changes.
4. **Role of the models that fix the analysis half-life. Resolved
   2026-10-30.** Models that assume a known half-life lack face
   validity as candidate analyses. The three groups are handled as
   follows:
   - **(i) Oracle analyses at the true half-life (G3, G8). Kept**,
     labeled "oracle benchmark, not a candidate analysis" throughout.
     G8 is the reference that prices the selection step. The size
     evidence of Section 2.8 is stated against it.
   - **(ii) Analyses at a prespecified value (G8 at 0.25 or 1 week,
     `E3fix_cr2`). Dropped** from paper 02, together with every result
     and sentence about them. F1 leaves with the F-family. The cut-list
     consequences are in Section 2.2.
   - **(iii) Binary-indicator analyses (E0, G1, G6). Kept.** These fix
     the half-life at 0 implicitly. They stay as the published and
     default reference analyses.

   **Cost of dropping (ii).** Paper 02 loses its evidence that a fixed
   1-week value stays within 2 points of the oracle at true half-lives
   of 0.5 to 2 weeks. A referee may still ask why the analyst should
   not simply fix the half-life. The answer the paper keeps is the
   binary indicator, which fixes it at 0. In family 1 it reaches
   0.563 power against the oracle's 0.640 (G1 against G8). Prespecified nonzero values
   are not evaluated, and the Limitations should say so in one
   sentence.
5. **Name for G9-null. Resolved 2026-10-30: renamed G10.**
   - The stored code stays `E7null_cr2`. The display mapping in the
     Rmd (`m_lab`, `rv_lab`) and the ladder script map it to G10, as
     they already map `E` codes to `G` codes.
   - G9 is described as "G10 with the half-life chosen on the
     interaction model".
   - Record G10 in NOTATION.md and in paper 02's glossary.
   - The name "G9-null" survives only in stored data labels and in the
     referee-variant text being replaced.

## 2. Paper 02 after the split

### 2.1 Scope statement (Introduction and Methods)

Add one paragraph to the Introduction. Suggested text:

> The biomarker in this study is predictive only: it changes the
> drug response and nothing else. A biomarker that is also
> prognostic raises a different question, namely whether the
> interaction test confuses an association with the natural history
> or the placebo response with moderation of the drug effect. A
> constant prognostic shift is absorbed exactly by the biomarker main
> effect of every analysis here and cannot move the interaction
> estimate. Prognostic associations that vary over time, and the
> design-expectancy terms that protect against them, are the subject
> of a companion paper [paper 15].

Keep the absorption check for the constant shift ($b_{TV}$,
"prog-TV-grow" in the stored data). It takes one sentence and
answers the first question a referee will ask.

### 2.2 Section-by-section cut list

Status values:

- **Delete:** remove the section.
- **Trim:** remove the prognostic or F1 content only.
- **Rewrite:** recast the text without it.
- **Keep:** leave as it is.
- **Move:** relocate the material to paper 15.

| Section (label) | Lines | Action |
|---|---|---|
| Setup chunk: `s28`, `s30`, `pdc`, `ceil` loads | 122-129 | Delete. These are the prognostic stress, random-effects, PD-check and ceiling data. |
| Setup chunk: `m_lab`, `focus`, `focus_cols`, `rv_lab`, `rv_method` | 140-164 | Trim. Remove F1, F1A, F1A0 and F1A0B. Remove `E3fix_cr2` ('G8-fixed'). Map `E7null_cr2` to 'G10' (decision 5). |
| Setup chunk: `rvT` | 168 | Delete. Keep `rvS`; see the slope row below. |
| Abstract | 230-273 | Rewrite: (a) Methods: drop "or add expectancy terms (F1, F1A0B)". (b) Results: drop the F1 sentence, the trend-bias sentence and the strength-0.6 artifact sentence. (c) "Eight scenario families" becomes six; see 2.3. (d) Conclusions: recommend G10, not G9, and add "for a biomarker that is predictive only". (e) Methods: rung 4 is G10; G9 is the side analysis that selects on the interaction model. (f) Results: report the G10 size against the oracle, and G9's half-point excess as the cost of selecting on the tested model. (g) Methods and Conclusions: drop "fix it" from the side analyses and the sentence "A fixed 1-week analysis half-life with CR2 holds nominal size for at most two points of power" (decision 4). |
| Introduction, ladder overview | 385-400 | Trim. Remove the F-family sentence and the prognostic, trend and large-N items from the list of stresses. Add the scope paragraph from 2.1. |
| Aims (`sec-aims`) | 457-480 | Trim. Remove "prognostic pathways, large samples, a biomarker that predicts the natural-history trend" from the second aim. Keep "patient-specific slope variance". |
| Table `tab-families` | 482-500 | Rewrite. Delete rows 5 (prognostic) and 6 (large N; every cell there is prog-PB). Row 7 becomes slope variance only. Renumber. |
| `sec-ladder-def`, table `tab-ladder-def` | 643-672 | Rewrite. Remove the pointers to the F-family. Redefine rung 4 as G10 (half-life chosen by AIC without the interaction). Label G3 and G8 "oracle benchmark, not a candidate analysis" (decision 4). Rename G9-null to G10 (decision 5). |
| `sec-g4` (rung 4) | 769-802 | Rewrite. Lead with null-model selection and why it keeps the test honest. Then present selection on the interaction model (G4, G9) as the variant whose size cost is measured. The existing F1A/F1A0 explanation (793-794) carries over. |
| `sec-e0` | 673-688 | Keep. E0 contains the expectancy regressor `De` because Hendrickson's analysis does. That is part of the original analysis, not an F-family feature. Add one clause noting `De` is constant, and so dropped, under CO. |
| `sec-comparators` | 803-839 | Rewrite. Delete the F-family paragraph (807-831). G10 moves out of the side analyses into rung 4. Delete the "G8 at a fixed half-life" paragraph (decision 4). The only side analysis left is G9 (selection on the interaction model). Retitle "Side analysis: selecting on the interaction model". |
| `sec-pathways` | 841-866 | Delete, except one sentence on constant-shift absorption, which goes in the scope paragraph. |
| `sec-choice` | 868-886 | Keep. The earlier `bm:L` and `bm:tsd` models are a carryover question. |
| `sec-grid` | 924-960 | Trim. Remove the prognostic-cell clauses (924, 929, 950-955). |
| `sec-corrections` | 1008-1058 | Trim. The PD-refusal correction stays as a statement about the driver. Its prognostic-scenario evidence moves to 15 (1041-1044). |
| Ladder-at-a-glance (`sec-glance`) and Fig. 1-3 captions | 1072-1274 | Rewrite after the ladder script is rerun (2.4). Remove F1 and F1A0B from the captions (1127, 1159) and from the prose (1206-1207, 1271). |
| `sec-rung1`, E0 isolation | 1282-1580 | Keep. The expectancy mentions here describe E0's `De` term, which stays. |
| Rung 4 null variants (`stage3-variants-null`) | 2077-2125 | Rewrite. Remove the F1 column and the 'G8 at 0.25 wk' and 'G8 at 1 wk' columns. The table becomes G8 (oracle), G9 and G10 in the eight Hybrid null cells. |
| `sec-alternatives` | 2172-2341 | Rewrite. Delete the fixed-half-life alternative (decision 4). Delete the F-family evidence (2176-2183, 2279-2339, the 'G9>F1' paired step). Move the G9-versus-F1A0B comparison to 15. What remains (null-model selection) has moved into rung 4, so fold any surviving text into `sec-g4` and delete the section. |
| `sec-scen-grids` | 2353-2457 | Trim. Remove the F1 columns (2381-2394, 2436-2448). |
| `sec-dropout` | 2458-2585 | Trim. Remove the F1 columns (2472-2483, 2499-2500, 2514, 2524-2531, 2565-2566). Line 2533 cites the no-dropout reference of the prognostic stress job. Re-anchor it on a non-prognostic reference or drop the comparison. |
| `sec-prognostic` (all subsections, including `sec-feasibility`, `sec-feasible-grid` and "Larger samples") | 2586-2939 | Move to 15. |
| `sec-trend`: trend pathway | 2940-3012 | Move to 15 (table `stage5-trend`, Fig. 4). |
| `sec-trend`: patient slope variance | 3014-3045 | Keep as its own short subsection, "Patient slope variance", with G8 and G9 only. It is a heterogeneity stress, not a prognostic one: the slope is unrelated to the biomarker. Delete the fig4 trend contrast (3028-3040). |
| `sec-random` | 3046-3090 | Move to 15. It was fitted on the prognostic stress cells, 8 of 12 of them the repaired strength-0.6 cells, and was motivated by the artifact. |
| `sec-calibration` | 3091-3178 | Trim. Remove the F1 $\kappa$ values (3142-3143). |
| `sec-which`, ladder score narrative | 3179-3262 | Rewrite from the regenerated table B, with G10 as the reference analysis ("percent of G10 power"). "G9 has the largest cell size ... under the trend pathway" no longer applies. Remove the F1 power-loss sentence (3250). "For a Hybrid trial G9 is the analysis to use" becomes G10, stated together with its power cost against G9. |
| Discussion, `sec-bought` | 3265-3364 | Trim. Remove the F-family and prognostic summaries (3331-3359). |
| `sec-recommend` | 3365-3414 | Rewrite. Recommend G10. Delete the fixed-1-week-half-life sentences (decision 4). Give the trade in one sentence: half a point of size bought back for 1 to 3 points of power. State the scope condition (predictive-only biomarker). Insert the joint-recommendation sentence (Section 4). Remove the trend caveat (3404-3406); it moves to 15. |
| "What the expectancy terms are for" | 3415-3442 | Move to 15. In 02, replace it with the one-sentence expectancy note in 2.5. |
| Limitations | 3443-3497 | Trim. Remove the repaired-cell, prognostic-grid and large-N F-family sentences (3469-3488). Add the large-N gap (2.6). Add one sentence: analyses at a prespecified nonzero half-life were not evaluated (decision 4). |
| Conclusions | 3534-3595 | Rewrite items 4, 5, 6, 7 and 8. Item 4 becomes "choose the half-life on the model without the interaction". Item 7 recommends G10. Item 5 keeps only the interaction-model-selection finding (G9's size excess); the fixed-half-life finding is deleted (decision 4). Item 6 (feasibility) moves to 15. Item 8 loses the prog-PB and trend sentences and gains the scope condition. |
| Glossary | 3596-3715 | Trim. Remove F1, F1A, F1A0, F1A0B, the prognostic labels, the feasible-grid note and the half-life policy value (3598, 3625, 3642). Replace the "G9-null; G8 at a fixed half-life" entry with a G10 entry and an "oracle benchmark" entry. |
| Reproducibility | 3751-3878 | Trim. Remove the jobs 28 and 30, scripts 30 and 31, and block T entries (3786-3792, 3829). |

### 2.3 Scenario families after the split

| New | Old | Content |
|---|---|---|
| 1 | 1 | Original 30-cell grid |
| 2 | 2 | Misspecified half-life and decay shape |
| 3 | 3 | Serial correlation |
| 4 | 4 | Patient-level dropout |
| 5 | 7 (block S only) | Patient slope variance |
| S | S (block N: $c_{bm,PB} = 0$ cells only; block A) | Independent-seed Hybrid run |

Block N of the referee variants holds four null cells at
$c_{bm,PB} = 0.3$. They leave family S. The pooled G9 null figure in the
abstract ("7,987 null fits ... 0.055") must be recomputed without
them. Check whether `m2-pooled-size.csv` (script 34) already pools
only $c_{bm,PB} = 0$; the ladder script sets
`pooled_null = c_pb == 0` for block N, but script 34 has not been
checked.

### 2.4 Ladder script (`37-ladder-tables-figures.R`)

The script reads stored data only, so regenerating it is cheap. The
changes are:

1. Drop family 5 (job 28 and the gap-fill sub `C1-ladder`), family 6
   (sub `C2-largeN`), block T, and the block N cells with
   $c_{bm,PB} > 0$.
2. Drop F1 and F1A0B from `analysis_cols`, the paired steps (the
   'G9>F1' row) and `side_x`. Drop them from the Fig. 2 and Fig. 3
   analysis lists. Drop 'G8-0.25' and 'G8-1' the same way
   (decision 4), and relabel 'G9-null' as 'G10' (decision 5).
3. Make G10 rung 4 and the reference analysis:
   - Change `mean_power_pct_of_g9` to percent of G10.
   - Set the paired steps to G8 to G10 (rung 3 to 4) and G10
     to G9 (side).
   - Move G9 to `side_x`.

   This needs the backfill of Section 2.8 first. Until then G10
   exists in family S only.
4. Delete Fig. 4 (trend). Paper 15 will need its own version.
5. Rerun the script and check that tables A, B and C and Fig. 1-3
   still satisfy the Rmd's `stopifnot` guards.
6. Re-read every inline `lv()`, `inc()` and `bv()` call against the
   new tables.

Rather than deleting code, add a `--scope predictive` switch that
filters these out. Paper 15 can then run the same script with
`--scope all` for its G-versus-F comparisons.

### 2.5 Expectancy note for 02

A referee will ask why the G-family omits a design covariate that
varies in Hybrid and OL+BDC. Answer in one sentence in
`sec-recommend`:

> With a predictive-only biomarker the G-family holds nominal size
> without the design-expectancy terms in every family examined; the
> terms matter when the biomarker may also be prognostic, which the
> companion paper [15] examines.

Do not claim more than that. The claim that F1 matched G9's power at
every alternative belongs to paper 15, which carries the evidence.

### 2.6 Evidence gap the split opens

Family 6 (large N) consists entirely of prog-PB cells. After the
split, paper 02 has no size evidence above $N = 70$ in the
predictive-only case.

**Decision (2026-10-30): state the limit; no large-N run.** Paper 02
reports size at $N \in \{35, 70\}$ only. It makes no claim about
larger trials. The text changes are:

- **Limitations.** Add: "Size was evaluated at $N = 35$ and
  $N = 70$. Behavior at larger samples was not examined for a
  predictive-only biomarker." This replaces the present partial
  statement ("measured only up to $N = 70$ with complete data") and
  the $N = 280$ sentences that move to paper 15.
- **Conclusions and `sec-recommend`.** Qualify the recommendation as
  holding "at $N \le 70$". Most aggregated N-of-1 trials fall in this
  range.
- **Do not argue that size problems shrink at larger $N$.** The
  stored data do not support it. Mean model-based size moves from
  0.033 to 0.028 (G1) and from 0.032 to 0.027 (G3) between $N = 35$
  and $N = 70$, consistent with conservatism that comes from model
  misfit rather than small samples. CR2's small-sample refinements
  would fade at larger $N$, but the need for a robust standard error
  would not. One sentence in the Discussion can say this, without a
  claim about where the boundary lies.

The run considered and declined was a Hybrid null at
$N \in \{140, 280\}$ with G8, G9 and G10, and model-based and CR0
arms of G10. Estimated at 1 to 3 hours on zeus. It remains available
if a referee asks.

### 2.7 Effort for 02

- 1 to 2 days of rewriting.
- About 1 hour for the ladder script and checks.
- The G10 backfill (2.8):
  - half a day for the driver
  - 7.7 local core-hours, timed 2026-10-31 (about 3 h on zeus)
- No large-N run (decision 2).

### 2.8 G10 as the recommended analysis: evidence and backfill

**Evidence in hand.** Block N, job 34: Hybrid, 5,000 replicates per
null cell, with G8, G9 and G10 fitted to the same trials. The
table shows the four cells with $c_{bm,PB} = 0$, the ones paper 02
keeps.

| Cell | G8 (oracle) | G9 | G10 | G9 $-$ G8 (paired SE) | G10 $-$ G8 (paired SE) |
|---|---|---|---|---|---|
| $t_{1/2} = 0.5$, $N = 35$ | 0.048 | 0.052 | 0.047 | +0.004 (0.0016) | $-$0.001 (0.0015) |
| $t_{1/2} = 0.5$, $N = 70$ | 0.051 | 0.054 | 0.051 | +0.002 (0.0014) | $-$0.000 (0.0014) |
| $t_{1/2} = 1$, $N = 35$ | 0.051 | 0.060 | 0.051 | +0.009 (0.0019) | $-$0.000 (0.0017) |
| $t_{1/2} = 1$, $N = 70$ | 0.046 | 0.051 | 0.045 | +0.004 (0.0017) | $-$0.001 (0.0016) |

The two-MCSE band at 5,000 replicates is $[0.044, 0.056]$. The
results are:

- **G9 is liberal against the oracle in every cell.** It is in all 8
  block N cells, including the prognostic ones, and 7 of those 8 are
  above 2 paired SE.
- **G10 is within 1 SE of the oracle everywhere.**
- **The standard errors agree with this.** $\kappa$ is 0.89 to 0.97
  for G9 and 0.94 to 0.99 for G10.
- **Power.** G10's cost is 1 to 3 points, from block A: 8
  alternative cells at 1,000 replicates, Hybrid only. In ladder
  table B it reaches 96.6% of G9's power.

**What is missing.** G10 was fitted only in family S (blocks N
and A, Hybrid, $N \le 70$). A recommended analysis must be shown in
every family the paper reports:

- family 1, the 30-cell grid over CO, Hybrid and OL+BDC (job 27)
- family 2, misspecified half-life and decay shape (job 21 cells and
  gap-fill B1, B2)
- family 3, serial correlation (gap-fill B4)
- family 4, patient-level dropout (gap-fill A1, plus the A2 size
  diagnostic)
- the slope-variance cells (block S)

**Backfill run.** Driver `41-run-g10-backfill.R` (written
2026-10-31; number 40 was taken), with the generator
`simulate_cell_g10()` added to `simulation-core.R`:

- It loops over the cell grids of `27-run-tier1-g-vs-f-e0.R` (job
  27), `06-g-vs-f-comparison.R` (the four Stage 4 cells, job 21),
  `32-run-g9-gap-fill.R` (A1 patient-level from job 33, A2 from job
  31, B1, B2 and B4 from job 32) and block S of
  `36-run-referee-variants.R` (job 37). B3-meanmod and job 31's
  per-occasion A1 cells are excluded; neither is in a paper 02
  family.
- For each (cell, replicate) it fits `E7null_cr2` and `E7cr2`
  through `fit_spec_s7()`.
- **Pairing.** Each replicate is regenerated with its stored job's
  seed (`base + r`, or `label_seed()` for block S). G9 is refitted
  only as a check: its estimate must match the stored value to
  $10^{-6}$. The smoke run (2 replicates of all 96 cells, local)
  reproduced every stored trial, with refitted-G9 differences up to
  $4.5 \times 10^{-8}$. That is optimizer and BLAS rounding on
  identical data; a different trial would differ by about 0.1. The
  tolerance was widened from the $10^{-8}$ first proposed, which
  flagged 5 of 192 replicates on rounding alone.
- **Platform.** Because the local machine reproduces the stored
  zeus trials, the run is made locally and needs no commit.
- **If the check fails** for a source, the stored G8 and G9 cannot be
  paired. Rerun G8, G9 and G10 together for it and pair within the
  new run.
- **Size:** 55,000 replicates (96 cells).
- **Cost (timed 2026-10-31):** 0.31 to 0.69 s per replicate on one
  local worker; 7.7 local core-hours in all. Projected at 3.0 h on
  zeus's 8 cores, calibrated against job 37.

**Result (run 2026-10-31, local, 181 min).** All 55,000 replicates
completed. Output: `analysis/data/carryover-sensitivity/g10-backfill/`.

- **Guard.** 54,989 of 55,000 refitted G9 estimates match the stored
  value to $10^{-6}$ (99.9th percentile of the difference
  $1.2 \times 10^{-7}$). Of the 11 that do not, 8 differ by
  $1.2$ to $6.7 \times 10^{-6}$ with the same half-life and p-value,
  1 failed to converge locally, and 2 selected a different half-life:
  one because a candidate fit is singular on this platform, one
  because a candidate's ML fit reached a different optimum (an AIC
  gap of 0.985 locally; inferred, not checked on zeus). None changes
  a rejection decision. The 11 are excluded from the paired
  contrasts.
- **Size.** Pooled over 25,915 paired null replicates: G8 0.0451,
  G9 0.0509, G10 0.0456. G9 $-$ G8 $= +0.0057$ (paired SE 0.0007);
  G10 $-$ G8 $= +0.0004$ (0.0007). G10 is within about 1.5 paired SE
  of the oracle in every family, while G9 is above G8 in every
  family. One null cell is outside its two-MCSE band for G10 (family
  4, MCAR 20%, $N = 70$: 0.026), on the conservative side, and G9 is
  outside it too.
- **Power.** G10 reaches 98.5% of G9's power over the 57 alternative
  cells (97.9% to 99.7% by family). The worst single-cell loss is
  2.8 points (family 1, Hybrid).
- **CO.** The design flagged below behaves like the others: null
  0.054 against G8 0.051 (paired difference $+0.003$, SE 0.003);
  power 98.5% of G9.

Neither overturning outcome occurred. The recommendation stands.

| Family | Null reps | G8 | G9 | G10 | G10 $-$ G8 (SE) | G10 power, % of G9 |
|---|---|---|---|---|---|---|
| 1 CO | 999 | 0.051 | 0.057 | 0.054 | +0.003 (0.003) | 98.5 |
| 1 Hybrid | 998 | 0.049 | 0.059 | 0.052 | +0.003 (0.004) | 98.1 |
| 1 OL+BDC | 998 | 0.050 | 0.053 | 0.051 | +0.001 (0.003) | 99.0 |
| 2 | 5,975 | 0.043 | 0.050 | 0.043 | +0.000 (0.001) | 97.9 |
| 3 | 1,994 | 0.044 | 0.047 | 0.043 | $-$0.002 (0.002) | 98.6 |
| 4 | 5,973 | 0.040 | 0.045 | 0.041 | +0.001 (0.001) | 98.6 |
| A2 size diagnostic | 7,978 | 0.049 | 0.055 | 0.049 | +0.000 (0.001) | n/a |
| Slope (block S) | 1,000 | 0.042 | 0.043 | 0.042 | +0.000 (0.003) | 99.7 |

**What would have overturned the recommendation.** Two outcomes would:

- a family in which G10 leaves the size band while G9 does not
- a family in which G10 loses much more than 3 points of power

The CO design deserves attention. There the candidate half-lives
differ least in their fixed-effect fit, so null-model AIC carries
little information about the half-life.

Expected grade after the cut: B, with the front matter (cover letter,
title page, declarations) still outstanding as in `docs/66`.

## 3. Paper 15 outline: prognostic biomarkers and the interaction test

Working title: *When the biomarker is also prognostic: bias in the
biomarker-treatment interaction test of aggregated N-of-1 trials and
two remedies.* Directory: `analysis/report/15-prognostic-biomarker/`.
Target: the same venue as 01 and 02, as a third companion, or a
methods journal on its own.

| Section | Content | Source (existing) |
|---|---|---|
| 1. Introduction | Predictive versus prognostic biomarkers; why the interaction test can confuse the two in N-of-1 designs with varying expectancy | `docs/49`, `ch19`, lexicon `docs/29` |
| 2. Taxonomy of pathways | Constant shift (absorbed); correlation pathways prog-BL, prog-TV-const and prog-PB; mean-shift trend pathway (unbounded) | 02 `sec-pathways`; `tidyverse/R/functions.R` trend block |
| 3. Feasibility | Closed-form joint bound on the prognostic correlation and $c_{bm}$ (Schur complement); feasible grid | `docs/65`; scripts 30 and 31; 02 `sec-feasibility`. **Paper 01 already has a joint-ceiling table** (`tab:prognostic`, `18-prognostic-joint-ceiling.R`). Cite 01 for it rather than duplicating; see 5.2. |
| 4. The repair artifact | Strength-0.6 cells: power loss from silent PD repair, not confounding | `docs/65`; jobs 28 and 30; 02 `sec-prognostic` "strength-0.6" subsection; `sec-random` |
| 5. Analyses | G-family (G8, G9) against the F-family (F1, F1A, F1A0, F1A0B); expectancy terms $\eta_t$, $B^c_i\eta_t$ (regressors `De`, `bmc:De`); why they absorb the expectancy-scaled on/off contrast under prog-PB | `docs/56`, `docs/57`, `docs/62`, `docs/64`; 02 `sec-comparators` and "What the expectancy terms are for" |
| 6. Correlation pathways: results | Size and power at strength 0.3 and on the feasible grid; large N under prog-PB | Job 28; gap-fill `C1-ladder` and `C2-largeN`; 02 tables `stage5-*` |
| 7. Trend pathway: results and the `bmc:t` remedy | Bias of every analysis; linear and quadratic `bmc:t`; cost when there is no trend | Block T (job 36); **new run** (3.1) |
| 8. G-family against F-family without a prognostic effect | The F-family's cost when unneeded; F1's fixed-half-life penalty versus F1A0B | Job 27 (`tier1-g-vs-f-e0`), job 21 (`g-vs-f-comparison`); 02 `sec-alternatives` F-material |
| 9. Recommendation | Joint recommendation (Section 4) | |
| Appendix | Constant-shift absorption proof; repair mechanics | `docs/65` |

### 3.1 New simulation for paper 15

**`bmc:t` block.**

- Add specifications `E2cr2_bmt`, `E7cr2_bmt` and `F1_bmt` to
  `simulation-core.R`: the existing formula plus `bmc:t`. For G9, add
  the term to both the AIC-selection formula and the refit formula.
  Optionally add quadratic variants.
- Add a new block to `36-run-referee-variants.R`. Cross the trend
  pathway $b_{TV}^{\mathrm{tr}} \in \{0, 0.3\}$ with $c_{bm} \in \{0, 0.45\}$
  at Hybrid, $N = 70$, $t_{1/2} = 0.5$, 1,000 replicates, plus a
  5,000-replicate null trend cell for size.
- Estimate about 30 to 45 minutes on zeus. Block T took 8 minutes for
  2 cells and 3 specifications.

A local pilot on 2026-10-29 (scratch script, 300 replicates, G8 only)
gave the following:

- The null bias of +0.26 fell to +0.10 with linear `bmc:t` and to +0.00
  with quadratic.
- Power fell from 0.96 to 0.94 and 0.92.

The DGP's trend follows the TV Gompertz course, not a straight line,
so linear `bmc:t` is expected to be only a partial remedy.

### 3.2 Effort for paper 15

About 1 to 2 weeks:

- Drafting is mostly assembly from the sources above.
- One new simulation block (about 1 hour of compute plus half a day
  of code).
- The ladder script runs with `--scope all` for the G-versus-F
  figures.

## 4. Joint recommendation (both papers, CLAUDE.md, NOTATION.md)

The two papers must give one answer between them. Present the
F-family as the G-family plus $\eta_t$ and $B^c_i\eta_t$ (regressors `De` and `bmc:De`):

- **Paper 02** settles the carryover layer with G10. That means
  an exposure-weighted predictor, with the half-life chosen by AIC on
  the model without the interaction, and CR2 with Satterthwaite df.
- **Paper 15** adds the expectancy layer when the biomarker may be
  prognostic, and a `bmc:t` layer if the biomarker may predict the
  natural-history trend.

**Resolved 2026-10-30: the F-family member is F1A0, and its
calibrated form is F1A0B.** F1A0 is G10's selection rule with
$\eta_t$ and $B^c_i\eta_t$ added. Choosing it makes the two layers
nest exactly: G10 and F1A0 differ only in the expectancy terms.
It also matches CLAUDE.md's provisional F1A0B.

Two things remain for paper 15 to show:

- whether F1A0B's bootstrap calibration is needed once the selection
  excludes the interaction (`docs/64` found no Type I inflation)
- what F1A0 costs against G10 when the biomarker is not
  prognostic

F1 at the 0.2-week policy value is not the recommendation. It loses
up to 10 points of power where that value is wrong (family 1).

Recommended wording:

> For a biomarker believed to be predictive only, use G10: the
> exposure-weighted predictor with the half-life chosen by AIC on the
> model without the interaction, tested with CR2 (paper 02). Where
> the biomarker may also be prognostic, add the design-expectancy
> terms $\eta_t$ and $B^c_i\eta_t$ (F1A0, paper 15). Where it may
> predict the rate of natural-history change, also add a
> biomarker-by-time term.

Once decided, update:

- the CLAUDE.md "Recommended analysis for the interaction (F1)"
  paragraph
- NOTATION.md (analysis IDs; the new paper number)
- `docs/66` (portfolio table)
- `docs/60` Part 2 (paper 02 plan)

## 5. Cross-paper effects

### 5.1 Paper 01

Paper 01 cites paper 02. Check that none of its citations points at
material that moves. Grep for "companion" in 01. As of 8e810ca9d,
01's own prognostic content is its Study 3 stress test, which is
independent of 02.

Paper 01 also runs the recommended analysis itself. Driver
`quick-sim/carryover-closed-form/19-g9-planning-and-fig4b.R` fits G9
(`E7cr2`) for the planning block P (63 cells, 72,000 replicates) and
the Figure 4B block F (18 estimable cells, 250 replicates, five
censoring settings each). These were stored by simrun on 2026-10-28.
The paper 02 backfill did not cover them. The driver now takes
`--specs E7null_cr2`, which fits G10 on the same trials (seeds are
`label_seed(cell_id-rep-r)`), and `--guard-reps 50`, which refits G9
on the first 50 replicates of each cell and compares the refit with
the stored run. A two-replicate check matched the stored G9 on all
126 refits. The full rerun writes to
`analysis/data/quick-sim/carryover-closed-form/g10-planning/`
(block F in `block-F/`). Paper 01's `planningg9data` and
`fig4bg9data` chunks, `summary.Rmd`, `21-fig4b-g9-heatmap.R`, the
text at the G9 passages and the glossary entry then move to G10.

### 5.2 Overlap between paper 01 Study 3 and paper 15

Paper 01 already contains:

- a stress test under a prognostic biomarker (Study 3, `tab:prognostic`)
- a joint feasibility table (`18-prognostic-joint-ceiling.R`)
- the statement that the stop-on-non-PD default came from "a
  prognostic scenario" with correlation 0.6 (lines 612-618)

Paper 15 must cite 01 for the DGP-level facts: the bound exists, and
repair is refused. It should present only the analysis-level
consequences. Paper 01 is at grade B and is not reopened for this
(decision 3, resolved 2026-10-31).

### 5.3 Portfolio count

The split adds a fifteenth manuscript. `docs/66` recommends
consolidating toward about nine. The split is justified because the
boundary is scientific, not editorial. It should be recorded
alongside the merges already proposed there: 13 with 14, 08 into 10,
12 into 01, and 11 into 01 or 03.

### 5.4 Deferred: bootstrap-calibrated selection (G9-WB), a paper 10 pilot

**Status: deferred (2026-10-30).** This is not part of the paper 02
work. Paper 02 proceeds with G10 whatever the outcome.

**Question.** Can one calibration step fix both size errors that
rungs 2 and 4 fix separately?

- the conservatism of the model-based standard error (CR2's job)
- the liberal excess from choosing the half-life on the interaction
  model (null-model selection's job)

If it can, G9's selection rule, which is 1 to 3 points more powerful
than G10, could be used at nominal size.

**Why a parametric bootstrap will not do.** Resampling from the
fitted working model (random intercept plus `corCAR1`) reproduces the
working model's error structure. That model is misspecified relative
to the DGP, which is the source of the conservatism, and under it the
model-based SE is correct by construction. Such a bootstrap would
therefore fix only the selection excess, which G10 already fixes
more cheaply.

**Candidate: G9-WB.** A null-imposed wild cluster bootstrap of the
whole procedure (Cameron, Gelbach and Miller 2008). Each resample:

1. Fits the model without the interaction.
2. Multiplies each patient's residual vector by one random weight.
   Use Webb six-point weights, given $N = 35$ clusters.
3. Rebuilds the outcome.
4. Repeats G9's AIC half-life selection on the interaction model.
5. Refits and records the model-based $t$ for `bm:Dbc`.

The observed $t$ is compared with the resampled distribution.
Patient-level resampling keeps the true within-patient covariance,
so no sandwich estimator is needed.

**Pilot design.**

- **Cells:** family S. That is block N (8 null cells, 5,000
  replicates) and block A (8 alternative cells, 1,000 replicates):
  Hybrid, $N \in \{35, 70\}$.
- **Arms:** G9-WB alongside G8, G9 and G10, on the same trials.
- **Bootstrap:** the warp-speed method (Giacomini, Politis and White
  2013). Each Monte Carlo replicate draws one bootstrap resample, and
  the critical value comes from the pooled resamples of the cell, as
  F1A0B does. This validates the procedure's size and power on
  average, not single-trial p-values.
- **Driver:** a new specification in `simulation-core.R` and a new
  block of `36-run-referee-variants.R`.
- **Cost:** half a day of code. About 1 to 3 hours on zeus, at
  roughly twice the cost of G9 per replicate (estimate, untimed).

**Adoption rule (exploratory, stated before the run).** Adopt G9-WB
only if both of these hold:

- its null rejection lies in the two-MCSE band in all 8 block N
  cells
- its mean power exceeds G10's by at least 2 points over block A

Otherwise record the negative result in paper 10 and stop.

**Reasons for low expectations.**

- F1A0B, the nearest existing bootstrap-calibrated analysis, reached
  92.9% of G9's power in paper 02's ladder, with a largest cell size
  of 0.070.
- `docs/64` found no power gain for F1A0B over F1.
- The most G9-WB could win is the 1 to 3 points that separate G9
  from G10.

**Home.** Paper 10 (interaction-test calibration). Its current
argument is mean-model misfit with CR2 as the remedy, and G9-WB is a
competing remedy. A trialist would run the full bootstrap (about 999
resamples times 5 fits, a few minutes per trial), so a positive
result would also need a single-trial calibration check that the
warp-speed pilot cannot give.

## 6. Order of work

All five decisions in Section 1a are settled (the last on
2026-10-31), so work can start at step 1.

1. Write and time `41-run-g10-backfill.R` (2.8) and run it. Check
   the G9 pairing guard before using any result. **Done 2026-10-31;
   guard passed (11 of 55,000 numerical exceptions, excluded).**
2. Confirm that G10 holds the size band in every family. If a
   family overturns the recommendation (2.8), stop and revisit
   before rewriting. **Done 2026-10-31; no family overturns it.**
2a. Rerun paper 01's planning and Figure 4B cells with G10 (5.1),
   check the G9 guard, then move paper 01's G9 tables, figure, text
   and glossary entry to G10 and re-render 01 and its summary.
   **Done 2026-11-01.** Outputs are in `g10-planning/` (block P) and
   `g10-planning/block-F/`. The G9 refit on the first 50 replicates
   of every cell matched the stored G9 in 3,143 of 3,150 block P fits
   and 4,492 of 4,500 block F fits, with no rejection flips; the rest
   are missing in the stored run, plus one block F fit at 1.02e-6.
   The block F guard has to read zeus job 56. The local
   `g9-planning/block-F/all.rds` predates the `censordata()` sign fix
   (24f187bc6), and against it only the uncensored fits match. Driver
   19 now reads job 56. G10 power is at most 0.022 below G9 in any
   cell; at the E9 sample size it is 0.984 to 1.000 in the uncapped
   cells, the same as G9. No paper 01 conclusion changed. One
   sentence of the paper 01 Figure 4B text was corrected regardless
   of G10: OL+BDC power at 0.2 week (0.46) is not within Monte Carlo
   error of its no-carryover value (0.58), as the G9 text had said.
3. Add `--scope` and the G10 reference to the ladder script, then
   regenerate the predictive-only tables and figures. **Done
   2026-10-31.** `--scope full` reproduces the stored CSV files byte
   for byte. `--scope predictive` writes to `ladder-predictive/` and
   to Figures 1 to 3 with a `-predictive` suffix. It uses 46,989
   guard-matched G10 backfill replicates; the other 8,000 are the
   A2-size diagnostic cells, which are not a ladder family. G10 is
   in the 0.04 to 0.06 band in 6 of 6 families, against 5 of 6 for
   G8. Pooled by family, G10 power is 96.6 (family S) to 99.7
   (family 7) percent of G9's. G10 has no
   coverage entry (the stored coverage covers G8 and G9 only). The
   Rmd still reads `ladder/`. 22 inline lookups name dropped steps
   or analyses (`G8>G9`, `G9>F1`, `G9>G9-null`, `G8>G8-1`,
   `G8>G8-0.25`, `F1`), and `tab-ladder-inc` lists families 5 and
   6. All of these move in step 4.
4. Cut and rewrite paper 02 per 2.2. Re-verify every inline number.
   Run `perl tools/notation-lint.pl`. Render through
   `tools/render.sh`.
5. Create `analysis/report/15-prognostic-biomarker/` and move the
   excised Rmd sections into a first draft, keeping their chunks
   working against the same data. Cite paper 01 for the feasibility
   bound and the refusal of repair; do not edit paper 01 beyond
   step 2a.
6. Implement and run the `bmc:t` block (3.1). Include F1A0 alongside
   G10.
7. Update CLAUDE.md (the recommended-analysis paragraph: G10,
   with F1A0/F1A0B as the prognostic-safe form), NOTATION.md,
   `docs/60` and `docs/66`.

## 7. Unverified assumptions in this plan

- The line ranges come from a grep of section headers and F1 and
  prognostic mentions, not a line-by-line read. Some cross-references
  outside the listed ranges will also need cutting.
- I have not checked whether `m2-pooled-size.csv` excludes the
  $c_{bm,PB} = 0.3$ null cells (2.3).
- The conservatism-from-misfit reading in 2.6 rests on two sample
  sizes ($N = 35$, $70$), mostly at 500 replicates per cell. It is
  consistent with the data but not a demonstrated trend.
- The cost of the `bmc:t` block is scaled from block T's runtime. Its
  findings rest on a 300-replicate G8-only local pilot.
- G10's standing as the best-sized analysis rests on the 8 block N
  null cells and 8 block A alternative cells (Hybrid, $N \le 70$).
  Its behavior in CO, OL+BDC, dropout, misspecified-shape and
  serial-correlation cells is unknown until the backfill runs.
- The backfill's pairing assumes that rerunning the stored seeds
  regenerates identical trials. The smoke run verified this on 2
  replicates of every cell (2.8); the full run checks every
  replicate.
- The deferred G9-WB pilot (5.4) has never been tried in this
  project. Two things about it are unverified:
  - that a wild cluster bootstrap with Webb weights behaves well at
    $N = 35$ with eight visits per patient
  - its cost
- The deferred pilot is named "G9-WB", not G10, because G10 is now
  the recommended analysis (decision 5).
