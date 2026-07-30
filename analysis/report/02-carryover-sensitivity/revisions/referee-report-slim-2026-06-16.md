# Referee Report: Paper 02 (slim / undergraduate version)
*2026-06-16 18:30 PDT*

**Manuscript.** "Which way of handling carryover gives the most
powerful test for a biomarker-by-treatment interaction? A student's
guide to a simulation study"
(`analysis/report/02-carryover-sensitivity/report-slim.Rmd`).

**Nature of the document.** This is an explicitly pedagogical, condensed
version of the full research manuscript (`report.Rmd`), retitled "a
student's guide" and pitched at an undergraduate statistics major. It is
not a journal submission and should not be judged as one; the review
below holds the *mathematics and empirical claims* to full rigor, but
assesses *exposition* against the document's own stated audience rather
than against a top-tier journal's. Where a finding also applies to the
full manuscript, I say so.

---

## 1. Summary

The document compares three ways of coding drug exposure in the analysis
model for an aggregated N-of-1 trial -- a binary on/off indicator (A1), a
continuous exposure-weighted predictor that decays with an assumed
half-life (A2), and the binary indicator plus a lagged "just-stopped"
flag (A3) -- and asks which yields the most power to detect a
biomarker-by-treatment interaction under carryover. Using a factorial
simulation (three decay shapes, two data-generating architectures, three
designs, two sample sizes, three biomarker strengths; 500 replicates per
cell), it concludes that A2's advantage is real but *conditional*: A2
leads only under Architecture B in the two designs with a withdrawal
phase, is strongly dominated under the crossover design, and is slightly
dominated under Architecture A; and that preventing dropout matters more
than the choice of method. The condensation faithfully preserves the
corrected, conditional conclusion of the full study.

## 2. Overall assessment

**Recommendation: minor-to-major revision (as a teaching document).** The
science is correctly represented and the central numbers are right (I
re-verified them by recompute). Three issues keep it from being a clean
teaching artifact: (i) a residual overstatement in the results narration
that contradicts the corrected conditional conclusion (Major M1); (ii) a
Type-I-error argument that assumes an independence the design does not
have (Major M2); and (iii) a register split -- the abstract, introduction,
and discussion are genuinely undergraduate, but the methods (§3) and
results-table narration (§4) are lifted verbatim from the research version
and are dense with undefined terms, which defeats the document's stated
purpose (Major M3-M4). None is a deep error; all are mechanical to fix.

## 3. Significance and novelty

Not assessed in the usual sense: this is a derivative exposition, not a
new contribution. The underlying study's contribution -- a
design-and-architecture-conditional comparison of three concrete
carryover specifications, with the finding that the popular
"always use the exposure-weighted predictor" advice is only conditionally
correct -- is sound and worth teaching. The condensation conveys it
accurately.

## 4. Major comments (correctness first)

**M1. The results narration still overstates the A2 advantage and
contradicts the paper's own corrected conclusion.** *Verified by
recompute.* In §"Power by analysis specification" the text reads "the A2
advantage is not an artefact of the reference cell: it emerges at
$t_{1/2}=0.5$ weeks and grows through $t_{1/2}=1.0$ weeks in each design
stratum," and earlier "A2 retains a modest but consistent advantage."
This paragraph was carried over unedited from the research version; it
contradicts the corrected abstract and discussion, which (correctly) state
that under the crossover (CO) design A2 is *dominated* -- A2 = 0.488 versus
A1 = 0.830 at the matched CO cell ($t_{1/2}=1.0$, Architecture B),
confirmed from `02-grid-summary.rds`. "In each design stratum" is simply
false for CO. **This finding applies to the full manuscript as well** (the
same §4 narration was not corrected there). *Remedy:* scope the claim to
"Architecture B, Hybrid and OL+BDC designs" and add the one-clause CO
exception, matching the discussion.

**M2. The Type-I multiplicity argument assumes independence the design
violates.** *Inspected.* §"Type I error control" argues the cell maximum
of 0.084 is acceptable because it falls within "the expected maximum among
540 independent binomial draws of size 500." The 540 null cells are not
independent: within a cell the three specifications share common random
numbers (stated in §3.6), and at $t_{1/2}=0$ the five decay-form cells are
exact duplicates (the decay multiplier is zero for every form). The
conclusion -- no specification inflates Type I -- is probably correct, but
the stated justification is invalid. **Also applies to the full
manuscript.** *Remedy:* either drop the "independent draws" sentence and
report the per-cell binomial check over the unique cells, or explicitly
note the dependence understates the expected maximum.

**M3. Register split defeats the stated audience.** *Inspected.* The
abstract, introduction, and discussion are well pitched for an
undergraduate and define their terms. But §3 (Estimands, Methods,
$D_{bc}$ construction) and the §4 table narration are the research prose
verbatim, and introduce, without definition, terms the stated reader will
not know: "estimand," "ADEMP," "corCAR1 working covariance," "MVN
differential correlation," "Wald confidence interval," "Satterthwaite
degrees of freedom," "MCSE," "second moment," "lasagna plot," "Tier 1 /
Tier 2." The companion slim (paper 01) solved exactly this with a short
"Key terms" box; this document has none. *Remedy:* add the same glossary
box after the introduction, and lightly translate the §3-4 prose into the
register of the abstract. This is the single highest-value revision.

**M4. The headline table is research-grade and unexplained for the
audience.** *Inspected.* Table 1 presents seven Morris performance
measures (power, bias, empirical SE, coverage, MSE, their MCSEs,
non-convergence). The undergraduate text never defines bias, coverage,
empirical SE, or MSE, and -- per the (correct) caveat now in the prose --
the bias/MSE/coverage columns are *not comparable across specifications*
because A1/A3 estimate a different coefficient than A2. Presenting four
non-comparable, undefined columns to a student invites exactly the
misreading the caveat warns against. *Remedy:* for the teaching version,
reduce the headline table to the two comparable, defined measures (power
and Type I error), or add a plain-language key and visually de-emphasise
the non-comparable columns.

## 5. Minor comments

1. **Dangling cross-reference.** §"Decay-form mis-specification" ends
   "This is taken up further in §5.2." The discussion was condensed and
   its subsections renamed, so §5.2 no longer denotes the intended
   content. *Suspected* (the bookdown build will not error on a prose
   "§5.2"). Re-point or delete.
2. **Architecture B estimand description vs scale.** §Estimands calls the
   Architecture B estimand "the on-drug biomarker-response correlation
   $c_{bm}$, dimensionless," then the table uses
   $\theta_{\text{true}} = -c_{bm}\sigma_{BR} = -3.6$ in BR units. The
   change of scale and the negative sign are never explained; for this
   audience, one sentence is needed.
3. **Repo-internal detail not needed for the reader.** The "Implementation
   note" (paths, `build_sigma_matrix()`, `original-extended` supporting
   only exponential) and parts of §3.6 reference code internals
   irrelevant to a student; cut or move to a footnote.
4. **"Lasagna plot"** is shown (Fig.) and named but never explained; for
   the audience, one sentence on how to read it would help, or drop the
   figure (it is the least self-explanatory of the nine).
5. **Orphan bibliography entry.** `sturdevant2021carryover` remains in
   `references.bib` but is no longer cited (the condensed discussion
   dropped the Sturdevant paragraph). Harmless but should be removed for
   hygiene.
6. **Unpublished-companion dependency.** The intro's "31%" and "1-3%"
   figures cite `pmsimstats_dgp_architectures2026`, and the conservatism
   discussion cites `pmsimstats_calibration2026`; both are unpublished. A
   student cannot follow these up. Acceptable for an internal teaching
   note, but worth a one-line "companion papers, in preparation."
7. **Math spot-checks pass.** I verified by hand that the toy-$D_{bc}$
   table values equal $(1/2)^{t_{sd}/\hat t_{1/2}}$ for all nine entries,
   and that all three decay forms (linear, exponential, Weibull) return
   exactly $0.5$ at $t_{sd}=t_{1/2}$, i.e. the half-life calibration is
   correct and Weibull reduces to exponential at $k=1$. No errors found.

## 6. Missing and questionable references

| Location | Issue | Action |
|---|---|---|
| `references.bib` | `sturdevant2021carryover` present but uncited after condensation | Remove the entry |
| Introduction ("31%", "1-3%") | Cites unpublished companion `pmsimstats_dgp_architectures2026` | Mark "in preparation"; no other source needed |
| §Type I / §S6 | Conservatism (6-10%, CR2) cites unpublished `pmsimstats_calibration2026` | Mark "in preparation" |

No incorrect attributions were found in the citations that remain; the
Sturdevant-Lumley mischaracterisation flagged in the full-manuscript
review does not appear here because that paragraph was cut.

## 7. Suggestions for strengthening

1. Add the "Key terms" glossary box (mirror the companion slim) -- the
   biggest single improvement for the audience (M3).
2. Add the corrective "where A2 loses" exhibit (the by-design /
   by-architecture advantage figure already generated under
   `readability-review/figures/fig2-advantage-by-design.png`). It makes
   the conditional finding visual and directly forecloses the M1
   overstatement.
3. Translate the §3 methods and §4 table narration into the abstract's
   register, and reduce the headline table to comparable measures (M4).
4. Consider dropping one or two of the nine figures (the lasagna plot in
   particular) for a 20-page teaching document; the figure budget is the
   main remaining space cost.

## 8. Scope of this review

**Verified by recompute** (from `analysis/data/02-grid-summary.rds`, this
session): reference-cell power A2/A1/A3 = 0.860/0.766/0.770 at
$t_{1/2}=1.0$; the CO/Architecture B counter-example A2 = 0.488 vs
A1 = 0.830; Architecture A Hybrid power 0.976 -> 0.960; A2 strictly highest
in 68/360 non-null cells (19%); mean null Type I 0.039.
**Verified by hand:** the nine toy-$D_{bc}$ table entries and the
half-life calibration of all three decay forms.
**Inspected:** the full slim source, the rendered structure, the citation
keys against `references.bib`, and the register/voice of each section.
**Not checked:** I did not re-run any simulation, did not re-derive the
`docs/19` calibration ($\sigma_{BR}=8$), and did not run fresh literature
searches for this slim (the full-manuscript referee report
`referee-report-2026-06-16.md` records the literature searches; the slim
cites a strict subset of the same works).

---
*Rendered on 2026-06-16 at 18:30 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/02-carryover-sensitivity/referee-report-slim-2026-06-16.md*
