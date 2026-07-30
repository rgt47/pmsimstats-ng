# Referee Report

*2026-06-16 15:55 PDT*

**Manuscript.** Robustness of carryover-mitigation analysis strategies
for biomarker-treatment interaction: a factorial simulation study
(`analysis/report/02-carryover-sensitivity/report.Rmd`).

**Target venue (assumed).** Statistics in Medicine (per cover letter
and CSL).

---

## 1. Summary of the manuscript

The manuscript reports an ADEMP-structured factorial Monte Carlo study
of how the analyst's choice of carryover specification affects power to
detect a biomarker-by-treatment interaction in aggregated N-of-1 and
crossover trials. Three analysis-side specifications are compared: a
binary on-drug indicator that ignores carryover (A1), an
exposure-weighted continuous predictor committing to an assumed decay
form and half-life (A2, the current `pmsimstats-ng` default), and a
binary indicator augmented by a lagged-treatment nuisance term (A3,
attributed to the Jones-Kenward crossover tradition). The principal
factorial crosses three data-generating decay forms (linear,
exponential, Weibull) against the three specifications under two
data-generating architectures, three designs, two sample sizes, three
biomarker-effect levels, and three carryover half-lives (a stated 540
cells at 500 replicates), supplemented by six single-axis sensitivity
blocks (S1-S6) and a high-precision 24-cell rerun. The headline claim
is that A2 dominates A1 and A3 by roughly ten percentage points of
power at a prazosin-calibrated reference cell, that this ranking is
invariant across the grid and to half-life and decay-form
mis-specification, and that a CR2 cluster-robust standard error restores
nominal size and recovers a few points of power without disturbing the
ranking.

## 2. Overall assessment and recommendation

**Recommendation: Major revision.**

The study is competently engineered, the simulation pipeline is
reproducible, the ADEMP scaffolding is appropriate, and the central
power finding (A2 attains higher power than A1 and A3 under carryover)
is real: I reproduced the headline numbers exactly from the archived
cell summaries (A2 = 0.860, A1 = 0.766, A3 = 0.770 at the reference
cell, t<sub>1/2</sub> = 1.0; mean null rejection 0.039, cell maximum
0.084 across 540 null cells). The work does not, in its present form,
clear the bar for a top-tier methodological journal for four reasons,
each developed below: (i) the bias, MSE, and coverage performance
measures are scored against a single estimand that is correct for only
one of the three specifications, rendering those cross-specification
comparisons uninterpretable; (ii) the A3 comparator, as implemented,
tests the same interaction coefficient as A1 and differs only by a
main-effect nuisance regressor, so the "A1 approx A3" result is close to
a structural near-certainty rather than an empirical finding, and the
"contradicts Jones-Kenward" framing is not supported by what was
actually computed; (iii) a load-bearing citation (Sturdevant and
Lumley 2021) is misattributed, and the novelty claim is not demarcated
against the existing aggregated-N-of-1 carryover-simulation literature;
and (iv) the submitted source is an unfinished internal draft in which
every author-rewrite ("rgt") block is Lorem ipsum placeholder. None of
these is fatal to the underlying science, but together they require
substantive revision and re-analysis before the manuscript is
reviewable as a finished contribution.

## 3. Significance and novelty

The narrow contribution -- a factorial comparison of three concrete
analysis-side carryover specifications, under both DGP architectures and
under decay-form and half-life mis-specification, with MCSE-reported
performance measures -- is a reasonable and useful methodological
exercise for the precision-medicine N-of-1 niche. It is plausibly novel
in the specific sense that I did not locate a prior study contrasting
exactly these three analysis specifications under decay-form
mis-specification in the biomarker-interaction setting.

However, the manuscript overstates the gap in the literature. The claim
(Introduction; cover letter) that "no systematic evaluation of this
question exists for the within-subject biomarker-interaction setting" is
not demarcated against directly adjacent work:

- Hendrickson et al. (2020), the acknowledged base case, already
  studies carryover-induced power loss across these designs.
- The aggregated-N-of-1 simulation literature has examined
  carryover-induced power loss across designs with multiple variance
  structures (e.g. the Stat-Med-adjacent "Comparison of Aggregated
  N-of-1 Trials with Parallel and Crossover Randomized Controlled
  Trials Using Simulation Studies," PMC6955665, 2019). This is not
  cited.
- Sturdevant and Lumley (2021) is itself a carryover-testing
  methods paper that is far more relevant to the substantive question
  than the reporting-convention use to which it is put here (see M3).

The contribution is, on balance, an applied/methodological sensitivity
study of solid but incremental significance. As written it fits a
specialist clinical-trials methods venue (Statistics in Medicine is
defensible; Contemporary Clinical Trials, Clinical Trials, or BMC
Medical Research Methodology are equally appropriate). To rise to the
top tier would require either a transportable theoretical result (e.g.
an analytic characterization of why the exposure-weighted predictor
dominates, rather than a single-parameterization simulation) or a
markedly broader and better-positioned empirical claim.

The reporting is contemporary in its adoption of ADEMP and MCSE
discipline, which is a genuine strength and is executed correctly.

## 4. Major comments

Correctness problems first.

**M1. Bias, MSE, and coverage are scored against an estimand that is
correct for only one of the three specifications (verified).**
Section 3.5 and Table 1 (`tab:headline`) report bias, MSE, and coverage
for all three specifications against a single calibrated target
theta_true = -c_bm * sigma_BR = -3.6. I confirmed from
`analysis/data/02-grid-summary.rds` that `true_value` is identical
(-3.6) across A1, A2, A3 in every one of the 324 non-null cells. But the
three specifications do not estimate the same parameter: A1 and A3 test
the coefficient on `bm:Db` (interaction with the binary on-drug
indicator), whereas A2 tests the coefficient on `bm:Dbc` (interaction
with the continuous decayed-exposure predictor). These are different
population quantities. The consequence is visible in the archived means
at the reference cell, t<sub>1/2</sub> = 1.0: A1 and A3 return mean
estimates near -2.51 and are scored as having bias +1.09 and coverage
0.82, while A2 returns -3.68 with bias -0.08 and coverage 0.97. The A1
and A3 "bias" is not estimator bias; it is the difference between the
`bm:Db` coefficient and a target defined for A2's `bm:Dbc`
parameterization. As reported, the bias/MSE/coverage columns therefore
do not measure what their column headers claim for A1 and A3, and the
cross-specification comparison on those three measures is invalid.

This matters because the abstract and Discussion present the headline
table as part of the evidence for A2's superiority, and a reader will
naturally read "coverage 0.83 for A1 vs 0.97 for A2" as A1 being a
poorly calibrated estimator, when in fact A1 is being measured against
the wrong target. Power and Type I error are unaffected -- each tests
that specification's own coefficient against zero, and those comparisons
are valid -- which is why the core power finding survives.

*Remedy.* Either (a) define a matched estimand for each specification's
coefficient and report bias/coverage/MSE against the matched target
(this requires deriving the population value of the `bm:Db` coefficient
under each DGP, which is non-trivial under decay because the binary
indicator is a coarsened exposure), or (b) restrict the cross-spec
comparison to power and Type I error and either drop the
bias/MSE/coverage columns or present them per-specification with an
explicit statement that the targets differ and the columns are not
comparable across rows. Option (b) is the honest minimum.

**M2. The A3 comparator, as implemented, tests the same interaction as
A1 and the "A1 approx A3" finding is near-structural (inspected).**
In `simulation-core.R` (`fit_spec`), A3 is
`Sx ~ bm + t + Db + bm:Db + L`: the lagged indicator `L` enters only as
a main-effect nuisance, and the interaction extracted is `bm:Db`, the
identical term to A1. A3 can therefore differ from A1 only through `L`
absorbing residual variance in the off-drug timepoints; it does not
model a biomarker-moderated carryover at all. This near-fully explains
the "A1 approx A3" result reported throughout Tier 1 and the sensitivity
blocks -- it is close to a mathematical near-certainty given the model,
not a discovered empirical regularity. It also undercuts the framing
(Sections 3.3.3, 4.3, 5.3, Conclusion 3) that the data "contradict the
Jones-Kenward expectation" that a form-free carryover indicator should
outperform a parametric model under mis-specification: the implemented
A3 is not the parametric-vs-nonparametric contrast that argument is
about, because the carryover term it adds is not in the interaction
being tested.

*Remedy.* State precisely what A3 estimates and why it is structurally
close to A1. If the intended contrast is a model-free carryover handling
of the interaction, implement it as such (e.g. a `bm:L` term, a
period-stratified within-subject contrast, or exclusion of contaminated
observations) and re-run. At minimum, temper every claim that the
results "contradict" or "do not support" Jones-Kenward; the present
design cannot adjudicate that question.

**M3. The Sturdevant-Lumley (2021) citation is misattributed
(inspected against the source).** Section 4.6 states that "Sturdevant
and Lumley developed the two-tier reporting convention we have adopted:
a principal factorial covering the main contrast, supplemented by
marginal sensitivity blocks anchored at a reference cell." That paper
(Contemporary Clinical Trials Communications, 2021) develops a
mixed-effects-model *test for carryover effects after treatment
cessation* in the threshold-crossing / open-label TROPHY hypertension
setting; it is not a methodology for structuring or reporting simulation
studies. The two-tier "principal factorial plus marginal sensitivity
blocks" structure is ordinary ADEMP practice (Morris et al. 2019, Sec.
7) and is not attributable to Sturdevant-Lumley. This is a citation used
for a claim the source does not make.

*Remedy.* Remove the attribution. Separately, engage Sturdevant-Lumley
on its actual content -- it is a genuinely relevant prior treatment of
carryover testing and its omission as a substantive comparator is a gap.

**M4. The "540-cell factorial" overstates unique design coverage; ~27%
of cells are exact duplicates (verified in code).** `carryover_decay()`
returns `rep(0, ...)` for every decay form when `halflife == 0`
(`implementations/tidyverse/R/functions.R:54`). At t<sub>1/2</sub> = 0
the five decay-shape levels (linear; exponential; Weibull at k in
{0.7, 1.0, 1.5}) therefore generate identical data-generating processes
and identical analyst predictors. Of the 180 cells at t<sub>1/2</sub> =
0, only 36 are unique; the remaining 144 (27% of the 540) are exact
replications differing only by an inert decay-form label. This does not
bias any result, but the repeatedly emphasized "540 cells" figure
inflates the apparent breadth of the design and the redundant cells
consume roughly a quarter of the compute.

*Remedy.* Report the unique-cell count (396) alongside the nominal grid
size, or drop the decay-shape factor at t<sub>1/2</sub> = 0, and note
the collapse explicitly so readers do not over-read the factorial's
resolution.

**M5. Central interpretive claims depend on an unpublished companion
(inspected).** The interpretation that the mean null rejection of 0.039
is `corCAR1` working-covariance conservatism rather than ordinary
sampling behavior, the quantitative claims that the standard error is
over-estimated by six to ten percent, and the entire S6 CR2 recovery
narrative are sourced to `@pmsimstats_calibration2026`, an `@unpublished`
companion (paper 10 in the repository). A referee cannot verify the
load-bearing calibration result, yet the manuscript uses it to recast
all reported power figures as "mildly conservative lower bounds" and to
justify the headline recommendation ("A2 with a CR2 standard error").

*Remedy.* Make the key calibration result self-contained in this
manuscript (a short subsection deriving or demonstrating the
conservatism on the present grid), or defer the S6 recalibration and the
"lower bound" framing until the companion is available to reviewers.
Co-submission of the companion would also suffice if the venue permits
linked review.

**M6. The multiplicity argument for the Type-I maximum assumes
independence that the design violates (inspected).** Section 4.5 argues
the cell maximum of 0.084 is acceptable because "the expected maximum
among 540 independent binomial draws of size 500 at p = 0.05 exceeds
0.07." The 540 null cells are not independent: the three specifications
within a cell share common random numbers (stated in Section 3.6), and
the t<sub>1/2</sub> = 0 decay-form duplicates (M4) are exact
replications. The expected-maximum heuristic is therefore optimistic.
The conclusion (no specification inflates Type I) is probably still
correct, but the stated justification should either acknowledge the
dependence or be replaced with a direct check (e.g. the per-cell
binomial test with a multiplicity correction over the unique cells).

**M7. The submitted source is an unfinished draft (verified).**
`claudecode.tex` renders three versions of every paragraph -- a bullet
summary, a blue-italic author-rewrite ("rgt") block, and a gray "orig"
block -- and every single `rgt` block in `report.Rmd` contains "Lorem
ipsum dolor sit amet ..." placeholder text. The rendered `report.pdf`
(299 KB) accordingly interleaves placeholder Latin throughout. The
manuscript is an internal drafting artifact, not a submittable
document. I reviewed the `.orig` prose as the manuscript of record. The
`strip-claudecode` tooling exists in the repository to collapse the
scaffold to a single clean prose stream; this must be run, and the body
proof-read as continuous prose, before the paper is reviewable as
finished. The current word count (10,772) is also long for the venue and
will fall once the duplication is removed.

## 5. Minor comments

1. **Estimand sign not motivated.** Section 3.2 / Table 1 use
   theta_true = -c_bm * sigma_BR = -3.6. The negative sign is never
   explained in the text; state the convention (the on-drug correlation
   maps to a negative symptom-scale coefficient under the calibration in
   `docs/19`).

2. **Architecture-B estimand described as "dimensionless correlation"
   but scored in BR units.** Section 3.2 calls the Architecture-B
   estimand the "on-drug biomarker-response correlation c_bm,
   dimensionless," yet bias is reported on the calibrated BR-unit scale.
   Reconcile the two descriptions explicitly.

3. **High-precision rerun is on a different design than the headline
   cell.** Section 4.7 confirms the ranking on the OL+BDC design, whereas
   the Tier-1 reference cell is Hybrid. This is fine but should be stated
   so readers do not infer the McNemar p = 6.6e-17 attaches to the Hybrid
   reference cell.

4. **MCSE formula and McNemar usage are correct.** The
   sqrt(p(1-p)/n_sim) power MCSE (Morris Table 6) and the paired McNemar
   on common datasets are appropriate and correctly applied. No change
   needed; noted for completeness.

5. **Coverage above nominal at t<sub>1/2</sub> = 0 (0.990).** Consistent
   with the documented conservatism; worth one sentence connecting it to
   the S6 discussion.

6. **"500 replicates" vs the rerun's 600 and the dev run's 50.** The
   replicate counts are internally consistent across abstract, Section
   3.5, and Section 4.7; no action required, but ensure the abstract's
   "500 replicates per cell" is not read as applying to the rerun.

7. **Bibliography hygiene.** Entries are clean and complete. Jones and
   Kenward (2014) is a book and correctly typed. `@unpublished`
   companion entry lacks an institution/URL; add a stable identifier so
   reviewers can locate it (see M5). Verify the Sturdevant-Lumley DOI
   resolves once its role is corrected (M3).

8. **Author/affiliation placeholders.** Title page and cover letter
   carry an unresolved roster ("Temperance Persons", gmail address,
   "[organisation]" GitHub placeholder). Acceptable at review stage but
   must be resolved before submission.

9. **Section cross-references.** Section 4.6 cites companion paper 05 as
   "confirmed in Section 4.4," and the architecture-dependence subsection
   cross-refers to paper 01. Verify these resolve to the intended
   subsections after the scaffold is stripped (numbering shifts when the
   bullets/rgt environments are removed).

## 6. Missing and questionable references

| Claim location | Issue | Suggested source |
|---|---|---|
| Sec. 4.6: "two-tier reporting convention" attributed to Sturdevant-Lumley | Misattribution; the source is a carryover-*testing* methods paper, not a simulation-reporting convention | Drop attribution; for the convention cite Morris, White & Crowther (2019), *Stat Med* 38(11):2074-2102, Sec. 7. Engage Sturdevant & Lumley (2021), *Contemp Clin Trials Commun* 22:100711 on its actual content |
| Introduction / cover letter: "no systematic evaluation ... exists" | Novelty claim not demarcated against adjacent N-of-1 carryover simulations | Selukar et al. / "Comparison of Aggregated N-of-1 Trials with Parallel and Crossover RCTs Using Simulation Studies," PMC6955665 (2019); cite and distinguish |
| Sec. 3.3.3 / 4.3: Jones-Kenward "expectation" that form-free indicators beat parametric models under mis-specification | Claim is used to frame the A3 result, but A3 as implemented does not realize that contrast (M2); the specific Jones-Kenward proposition should be quoted | Quote the precise passage in Jones & Kenward (2014), 3rd ed., or remove the framing |
| Sec. 4.4 / Discussion: carryover erodes biomarker-response correlation "directly" under Architecture B | Mechanistic claim asserted without citation to the companion's derivation | Cite the specific section of paper 01 / `docs/19` |
| Sec. 4.5 / S6: corCAR1 conservatism, SE over-estimated 6-10% | Sourced only to an unpublished companion | Make self-contained or provide a citable companion (M5) |

## 7. Suggestions for strengthening the paper

1. **Derive, do not only simulate, the A2 advantage.** A first-order
   argument for why a correctly weighted exposure predictor recovers more
   interaction information than a binary indicator under exponential
   decay (an attenuation/efficiency calculation) would lift the
   contribution from "one calibration's simulation" toward a
   transportable result and would directly address the
   single-parameterization limitation the authors concede.

2. **Replace or supplement A3 with a genuine model-free carryover
   handling of the interaction** (M2), so the Jones-Kenward comparison is
   actually tested. As it stands the third arm adds little.

3. **Report a matched-estimand bias/coverage analysis** (M1), or move to
   a power-and-size-only comparison and reallocate the table space to a
   per-design power table with MCSEs, which is what the argument needs.

4. **Add at least one non-prazosin trajectory parameterization** to
   substantiate the repeated claim that the *ranking* (not the absolute
   power) is portable. This is currently asserted as expectation, not
   shown.

5. **Probe one genuine joint perturbation** beyond S5 -- high rho with
   high dropout is the obvious candidate the authors themselves flag --
   since the OFAT structure cannot detect the interactions most likely to
   matter operationally.

6. **Reproducibility.** Pin a single R version (the manuscript states R
   4.5.3 in Section 3.6 and R 4.4.x in the Reproducibility section -- a
   contradiction to fix), and report the unique-cell count (M4) in the
   grid description.

## 8. Scope of this review

**Verified (executed or recomputed from archived outputs).** Headline
reference-cell power, bias, coverage, MSE, and non-convergence
(`02-grid-summary.rds`): A2/A1/A3 = 0.860/0.766/0.770 at t<sub>1/2</sub>
= 1.0, reproduced exactly. Identity of `true_value` (-3.6) across all
three specifications in all 324 non-null cells (basis of M1). Mean null
rejection 0.039 and cell maximum 0.084 across 540 null cells. Grid size
arithmetic (5x3x2x3x3x2x3 = 540) and the t<sub>1/2</sub> = 0 decay-form
collapse (basis of M4). The `carryover_decay()` branches and the
`halflife == 0` early return (`functions.R:50-66`). The A1/A2/A3 model
formulae and the A3 `L`-as-main-effect construction
(`simulation-core.R`, `fit_spec`). Existence and non-emptiness of the
companion calibration manuscript (paper 10).

**Inspected (read and judged, not independently recomputed).** The
ADEMP pre-registration and `simulation-grid-plan.md` (S1-S6 match the
manuscript). The CR2 / GEE robust fitters in `simulation-core.R`. The
S6 table-generating chunk (`tab-s6`). The rendering scaffold
(`claudecode.tex`) and the presence of Lorem ipsum in every `rgt` block.

**Not checked.** I did not re-run the 540-cell production factorial, the
sensitivity blocks, or the S6 CR2 computation from scratch; I relied on
the archived `.rds`/`.txt` summaries. I did not verify the `docs/19`
calibration derivation (sigma_BR = 8, theta_true = -3.6) beyond
confirming the code uses br sd = 8. I did not verify the companion
calibration paper's internal numerical claims (the 6-10% SE
over-estimation). I did not inspect the figure PDFs pixel-by-pixel; I
read their captions and the underlying summary values.

**Literature searches run.** (1) "Sturdevant Lumley statistical methods
testing carryover effects mixed effects model" -- established the M3
misattribution. (2) "N-of-1 trial aggregated biomarker treatment
interaction simulation carryover analysis specification" -- surfaced
Hendrickson (2020) and the 2019 aggregated-N-of-1 carryover-simulation
comparison (PMC6955665) used in Section 3 and M4. I did not run a
dedicated search for exposure-weighted / dose-on-board predictor theory
or for recent (2022-2026) N-of-1 analysis-model developments; the
authors should close that gap before resubmission.

---
*Rendered on 2026-06-16 at 15:55 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/02-carryover-sensitivity/referee-report-2026-06-16.md*
