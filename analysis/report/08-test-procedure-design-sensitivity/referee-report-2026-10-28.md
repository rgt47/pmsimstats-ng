# Referee report: Test-procedure choice for the biomarker-treatment interaction in aggregated N-of-1 trials
*2026-10-28 17:40 PDT*

Manuscript: `analysis/report/08-test-procedure-design-sensitivity/report.Rmd`
(revision of 2026-10-08). Reviewed as for *Statistics in Medicine*.
Prior report: `referee-report-2026-06-17.md` in the same directory.
`report-slim.Rmd`, `report-mini.Rmd` and `bullets.Rmd` exist in the
directory and were not reviewed.

Evidence labels used below: VERIFIED (I executed code or recomputed
numbers), INSPECTED (I read the source and judged it), SUSPECTED (a
discrepancy I could not fully confirm), UNVERIFIED (not checked).

## 1. Summary

The manuscript compares five test procedures for the
biomarker-by-treatment interaction in an aggregated N-of-1 design: a
strict split-plot RM-ANOVA on a median-dichotomized biomarker with
phase pre-averaging, a continuous-biomarker within-subject contrast
regressed on the biomarker, the package's `nlme::lme` analysis with
`corCAR1`, GEE with the Liang-Zeger sandwich, and GEE with the
Mancl-DeRouen (2001) bias-corrected sandwich. Section 2 derives the
closed-form RM-ANOVA interaction $F$ and its reduction to a two-sample
$t$-test on within-subject differences. The simulation runs a 40-cell
factorial (5 procedures, $\beta_{bm} \in \{0, 0.45\}$,
$N \in \{25, 35, 70, 100\}$) on the package DGP under mean moderation,
with no carryover, 2000 null and 1000 alternative replicates per cell,
reported with bias, empirical and model SE, and coverage. The headline
results are that dichotomization accounts for 0.08 to 0.21 of absolute
power and the test procedure for at most 0.04; that the uncorrected
GEE sandwich over-rejects at 1.5 to 1.9 times nominal at every $N$;
that the Mancl-DeRouen correction brings the rate near nominal at a
7 to 10 point power cost; and that the `lme` analysis is conservative
(size 0.02 to 0.03). The cycle-by-period design grid is deferred.

## 2. Overall assessment and recommendation

**Recommendation: Major revision.**

The revision responds well to the earlier report's main items. The
continuous-biomarker contrast arm now isolates dichotomization from
the test procedure (prior M2), the misattributed 13 to 17 percent
figure is corrected and explicitly labeled as the `lmer` result (prior
M3, VERIFIED against `docs/legacy/06`), the sample-size grid is
restored (prior M4), the secondary estimands are computed and
tabulated (prior M6), the seeding statement now matches the code in
outline (prior M7), the inline Mancl-DeRouen SE is benchmarked against
`geesmv` (prior suggestion 4), and the title no longer promises the
design grid (prior M5, partly). The stored summary reproduces every
number in Tables 2 to 4 and in the Results prose (VERIFIED), the
covariance matrices used by the driver are positive definite at the
driver's parameters so the stored run does not depend on the
now-disabled non-PD repair (VERIFIED), and convergence was 100 percent
across 60,000 fits (VERIFIED from the replicate file).

Four problems remain that I regard as blocking. First, the paper's
diagnosis of, and remedy for, the `lme` conservatism is wrong in
direction: a Kenward-Roger or Satterthwaite correction can only make a
test more conservative, and the over-coverage the paper documents is a
variance-calibration defect that the compendium's own companion papers
have already diagnosed and remedied with a cluster-robust standard
error. Second, the Mancl-DeRouen arm uses a normal reference
distribution, and the paper's claim that the correction "restores
nominal Type I control" is contradicted by its own null cell at
$N = 25$ (6.35 percent, $p = 0.008$ against 0.05). Third, the
linear-mixed analysis evaluated here is the compendium's historical
analysis (`bm:Dbc` with model-based SE, binary exposure), not the
analysis paper 02 now recommends (G9: AIC-selected half-life with CR2
Satterthwaite inference) or the F1 model, so the paper's "preferred
default" recommendation is out of step with the series it belongs to.
Fourth, the Introduction and Section 2 still argue the thesis of the
pre-revision draft (sphericity as the mechanism, the joint
design-by-test study as the contribution, five DGP violations that the
executed study does not contain), so the paper contradicts its own
Results in several places. None of these requires new large
simulations; the CR2 arm is the only new computation I would ask for.

## 3. Significance and novelty

The defensible contribution is narrow: a clean decomposition, in the
aggregated N-of-1 setting with a two-path open-label plus blinded
discontinuation design, of the power loss from dichotomizing a
continuous moderator versus the power loss from the test procedure,
together with a documented small-sample GEE comparison and a
calibration table. That is useful, and the decomposition is the right
way to answer the question the earlier draft asked badly.

Two things limit the significance. The cost of a median split is one
of the oldest quantitative results in applied statistics: Cohen (1983)
gave the efficiency loss of dichotomizing a normal predictor as about
36 percent, MacCallum et al. (2002) and Royston, Altman and Sauerbrei
(2006) documented the power loss in regression at length, and none of
this literature is cited. The paper's 0.538 to 0.358 drop at $N = 25$
is that known result observed once more. The paper should present it
as a confirmation in this design, not as its first finding. The GEE
result is likewise known in outline (Mancl and DeRouen 2001; Pan and
Wall 2002; Li and Redden 2015), and the paper's own numbers reproduce
the standard caveat that a bias-corrected sandwich with a normal
reference still over-rejects at small cluster counts. What is new is
the N-of-1 interaction estimand and the side-by-side calibration of
the four continuous-biomarker arms. The novelty would be materially
greater if the paper included the analysis the compendium actually
recommends (CR2) so that a reader could see, in one table, how the
historical `lme`, the CR2-robust `lme`, and the corrected GEE compare.

## 4. Major comments

**M1. The Kenward-Roger/Satterthwaite remedy for the `lme`
conservatism is wrong in direction (Sections 5, 6 and 7.1).**
Location: Results paragraph "Fourth, the comparison is not
contaminated...", recommendation table row "$N \le 35$, concern about
LME conservatism", Discussion item four. Problem: the paper attributes
the sub-nominal size (0.021 to 0.034, VERIFIED, $z$ between $-3.4$ and
$-6.1$ against 0.05 in every null cell) to "`corCAR1`
degrees-of-freedom accounting" and recommends a Kenward-Roger or
Satterthwaite correction to "address the conservatism and recover the
small power increment it forgoes". Both corrections inflate the
estimated variance (Kenward-Roger) and reduce the denominator degrees
of freedom; each change can only lower the rejection rate. They are
remedies for anti-conservative tests, not conservative ones. The
paper's own Table 4 shows the mechanism: the model-based SE exceeds the
empirical SE by 11 to 16 percent (VERIFIED), which is a misspecified
working covariance (random intercept plus CAR(1) fitted to a
three-component DGP with a separate baseline level), not a
degrees-of-freedom effect. `nlme` assigns within-participant terms
such as `bm:Dbc` the inner-level degrees of freedom, which here exceed
150 even at $N = 25$, so no df choice can produce a 2 percent size. Why
it matters: the recommendation table tells a reader to apply a
correction that moves the test further from nominal. The compendium's
companion papers already diagnose this exact defect (paper 10,
`10-interaction-test-calibration`, the $\kappa$ ratio of NOTATION.md)
and remedy it with a CR2 cluster-robust standard error with
Satterthwaite degrees of freedom (paper 02, rung G6, which moved size
from the G1 value to near nominal). Remedy: delete the Kenward-Roger
recommendation; state the conservatism as a working-covariance
calibration defect with $\kappa = 1.11$ to $1.16$; add a sixth arm,
`lme` with CR2 Satterthwaite inference (`clubSandwich`, as in paper
02's driver), and report its size, power and coverage alongside the
others. INSPECTED (lme_analysis.R, nlme df convention) and VERIFIED
(summary file).

**M2. The Mancl-DeRouen arm uses a normal reference and the paper
overstates its size control (Abstract, Sections 5, 6, 8).**
Location: `mancl_derouen_se()` returns `p = 2 * pnorm(-abs(z))`
(driver lines 317 to 318, INSPECTED); Abstract "the Mancl-DeRouen
correction removes this"; Conclusions "restores nominal Type I
control". Problem: the corrected null rates are 0.0635, 0.053, 0.0625
and 0.059. Against 0.05 with 2000 replicates, the $N = 25$ cell is at
$z = 2.77$ ($p = 0.008$) and the $N = 70$ cell at $z = 2.56$
($p = 0.014$) (VERIFIED by exact binomial test). "Within or close to
MCSE of nominal" is not accurate for those two cells; 0.0635 is 1.27
times nominal and 2.7 MCSE away. The standard pairing of the
Mancl-DeRouen estimator with a $t$ reference on $K - p$ degrees of
freedom (Li and Redden 2015 conclude that the Wald $z$ test "should be
avoided ... even when bias-corrected sandwich estimators are used") is
exactly the step this arm omits, and the residual over-rejection is
the expected consequence at $K = 25$ clusters and $p = 5$ parameters.
Why it matters: the GEE recommendation row and the Abstract rest on
the claim of restored size. Remedy: use $t_{K - p}$ (or report both),
rerun the GEE cells (cheap: the fits are unchanged, only the reference
distribution changes, so the stored `beta` and `se` suffice), and
revise the claim to "near nominal" with the numbers beside it. Also
rename "naive sandwich": in the GEE literature "naive" denotes the
model-based variance, and what the paper calls naive is the standard
Liang-Zeger robust sandwich. "Uncorrected sandwich" is unambiguous.
VERIFIED and INSPECTED.

**M3. The linear-mixed analysis evaluated is the compendium's
superseded analysis, and the paper does not say so (Sections 4, 6).**
Location: "the standard pmsimstats linear-mixed analysis with
`corCAR1`" (Abstract, Section 4), recommendation table "Default".
Problem: `fit_lme()` calls `lme_analysis()` with `carryover_t1half =
0` and `useDE = TRUE`, which fits `Sx ~ bm + t + Dbc + bm:Dbc + De`
with a random intercept, `corCAR1`, model-based SE and the `nlme` Wald
$t$ (INSPECTED). Because the DGP has no carryover and the analysis
half-life is zero, `Dbc` is the binary `D` at every visit, so this is
specification E1 (paper 02's G1) at the true half-life. The
compendium's current recommendation is G9 (exposure-weighted
predictor, half-life by AIC, CR2 Satterthwaite inference; paper 02
Conclusions), with F1 (`Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De`,
CR2) as the validated fixed-half-life fallback (CLAUDE.md, "Recommended
analysis for the interaction"). Nothing in this study exercises the
half-life question, which is defensible given $t_{1/2} = 0$, but the
inference method differs, and the paper's historical analysis is the
one paper 02 shows to be miscalibrated. Why it matters: a reader of
the series will find paper 08 recommending as "preferred default" an
analysis that paper 02 and paper 10 replace. Remedy: name the fitted
model and its inference explicitly (formula, random structure,
model-based SE, `nlme` df), label it G1 with a cross-reference to
paper 02's Table of codes, add the CR2 arm (M1), and scope the
recommendation to what was tested. State also that the DGP here is
the package default configuration under mean moderation (AR(1) with
$\rho = 0.7$, cross-component correlations 0.1 and 0.05) and is
neither "the original DGP" nor "the extended DGP" of NOTATION.md
($c_1 = 0.2$), so a reader can place it. INSPECTED.

**M4. The design is mislabeled and the deviations table is incomplete
(Sections 1, 3, 4, 6).** Location: "the Hendrickson hybrid
open-label-plus-blinded-discontinuation design", "Hendrickson hybrid
OL+BDC" (Sections 1.1, 3, 6), Table 1 (deviations). Problem: the
driver builds a two-path OL+BDC design with visits at weeks 4, 8, 12,
16, 17, 18, 19, 20 and discontinuation at visit 7 or 6 (INSPECTED).
NOTATION.md distinguishes OL+BDC ($P = 2$) from the Hendrickson Hybrid
($P = 4$, weeks 4, 8, 9, 10, 11, 12, 16, 20, with a crossover block),
and the working protocol specifies "8-week OL + 4-week BDC + 4-week +
4-week crossover blocks". The executed study has no crossover block
and a different visit schedule. Table 1 does not record this change,
nor the substrate change (`implementations/simple` in the protocol and
README; the package in the paper), nor the change in the RM-ANOVA
arm's form (the protocol's "Greenhouse-Geisser-corrected $F$" on the
multi-timepoint series; the executed arm pre-averages to a $2 \times
2$), nor the seed mechanism (protocol: "per-replicate seed via
cell-descriptor hash"; executed: `mc.set.seed`). Section 6 says the
cycle count and period length are "held at the reference values
defined in Section 3", but Section 3 now defines nothing. Why it
matters: the design determines the number of off-drug visits per
participant (3 or 4 here) and therefore the variance of the phase
contrast that drives every power number; and the ADEMP discipline the
paper invokes requires that protocol departures be listed. Remedy: give
the executed design in full (paths, visit weeks, on/off pattern,
expectancies) in Section 4, call it OL+BDC throughout, and add the
missing rows to Table 1 with a one-line rationale each. INSPECTED.

**M5. The Introduction and Section 2 argue the pre-revision thesis and
contradict the Results.** Location: Section 1.1 ("dominates the strict
classical RM-ANOVA test in essentially every parameter cell where
within-participant residuals depart from compound symmetry"); Section
1.1 opening ("an analyst's two principal choices ... should be treated
jointly"); Section 1.5 contributions two and four (a design grid in
Section 3 and a recommendation table "indexed by the carryover
half-life and the target effect size"); Section 2.4 (five assumptions
"violated to varying degrees in the pmsimstats DGP", including
carryover at $t_{1/2} = 1$ week with "attenuation by 30 to 35 percent
in the production data" and missing data); Section 2.4 item 1 ("SE
roughly 1.4 to 1.7 times ... a 30 to 50 percent power penalty").
Problem: the executed study has no carryover, no dropout, equal
spacing within the blinded block only, and a pre-averaged RM-ANOVA
whose Type I error is nominal, so the compound-symmetry mechanism of
Section 1.1 and three of the five violations of Section 2.4 are not
present in the data analyzed. The 30 to 35 percent and 1.4 to 1.7
figures cite no stored artifact and are inconsistent with the paper's
own decomposition (the dichotomization cost is 0.08 to 0.21 absolute;
the ratio of contrast to RM-ANOVA power at $N = 25$ is 1.5, not a
30 to 50 percent "penalty" in any stated sense). Contributions two and
four are not delivered. Why it matters: a reader who stops at Section
2 takes away the opposite of what Section 5 shows. Remedy: rewrite
Section 1.1 and 1.5 around the delivered decomposition; cut Section 3
to one paragraph of future work or remove it; in Section 2.4 separate
"assumptions of the closed form" from "which of them the present study
exercises" (only the categorical biomarker); delete the unverifiable
production-data figures or cite the artifact that produced them;
remove the parenthetical "(corrected 2026-10-08 ...)" change-log notes
from the manuscript text, which belong in a revision log, not in a
submission. INSPECTED; the 30 to 35 and 1.4 to 1.7 figures are
UNVERIFIED and I could find no source for them.

**M6. The dichotomization and small-sample GEE literatures are
missing, and one citation is wrong.** Location: Sections 1.3, 2.4, 5,
7. Problem: no citation to Cohen (1983), MacCallum et al. (2002),
Royston, Altman and Sauerbrei (2006) or Altman and Royston (2006) for
the headline finding; `@pan2001` (Pan 2001, "Sample size and power
calculations with correlated binary data") is cited as a bias-corrected
sandwich reference, but the relevant paper is Pan and Wall (2002,
*Statistics in Medicine* 21:1429-1441), and the canonical
corrections alongside Mancl-DeRouen (Kauermann and Carroll 2001; Fay
and Graubard 2001) and their comparison (Li and Redden 2015) are
absent. Why it matters: the novelty claim in Section 1.3 ("nor does it
isolate the dichotomization cost") is not sustainable against the
dichotomization literature, and the GEE discussion cannot explain its
own residual over-rejection without the $t$-reference literature.
Remedy: see the table in Section 6 below. INSPECTED (bib) and VERIFIED
(PubMed metadata for Pan and Wall 2002, Li and Redden 2015, Royston
et al. 2006).

**M7. Reproducibility and implementation statements need
correction.** Location: Reproducibility section; Section 4; Abstract.
Problems, each INSPECTED or VERIFIED:

- "each replicate draws from an independent, reproducible
  random-number stream" is not what `mclapply(mc.set.seed = TRUE,
  mc.preschedule = TRUE)` does. Each forked worker receives one
  L'Ecuyer stream and processes a prescheduled block of jobs
  sequentially from it; reproduction therefore requires the same
  worker count (7, stated in Results but not in the Reproducibility
  section) and the same prescheduling. Say so.
- The driver passes `makePositiveDefinite = TRUE` and does not set
  `stopOnNonPD`. Under the package as of 2026-10-26 the default
  refuses non-PD matrices. I rebuilt both paths' matrices at the
  driver's parameters under the refusal default: minimum eigenvalue
  0.168 for both paths and both $\beta_{bm}$ values, and `c.bm` does
  not enter $\Sigma$ under mean moderation (VERIFIED). The stored run
  therefore did not rely on repair. The manuscript should state this
  and the driver should drop the obsolete argument.
- The GEE working correlation is `geepack` `corstr = 'ar1'` on
  observation order within participant, while the visits are at weeks
  0, 4, 8, 12, 16, 17, 18, 19, 20. This is a discrete-lag AR(1) on
  unequally spaced occasions, not a continuous-time AR(1). Section 4
  does not state the working correlation at all. State it, and note
  the mismatch with the `corCAR1` arm.
- The Mancl-DeRouen benchmark is reported as "within approximately
  2.5 percent". The mean ratio is 1.0245 (SD 0.021), but the
  per-replicate ratio ranges from 0.999 to 1.092 (VERIFIED from
  `08-test-mancl-benchmark.txt`). Report the range.
- The contrast and GEE arms include the baseline visit as an off-drug
  observation and the GEE and `lme` arms include `t` and (for `lme`)
  `De`; the contrast arm includes neither. The "contrast to
  linear-mixed step" therefore also changes the mean model, not only
  the covariance treatment. Say so, or drop `De` from the `lme` arm
  for this comparison.

## 5. Minor comments

1. Abstract: "The closed-form strict RM-ANOVA test agrees with the
   linear-mixed Wald test under sphericity and balanced design"
   overstates Section 2.3, which gives three conditions (dichotomized
   biomarker, compound symmetry, binary `Dbc`) and asymptotic
   equivalence only. Align the two.
2. Section 2.3.1: "Under the alternative $H_1: \hat L = \Delta_0$"
   should read $\mathbb E[\hat L] = \Delta_0$; $\hat L$ is the
   estimate. The non-centrality and variance algebra are correct
   (VERIFIED by re-derivation: $SS_{AB} = nT\hat L^2/4$,
   $\mathrm{Var}(\hat L) = 4\sigma_e^2/(nT)$, $t^2 = F$).
3. Section 2.2.2: the exactness claim "for any within-participant
   covariance" needs the $\Delta_{il}$ to share a variance across
   participants. In the executed design path A has 6 on-drug and 3
   off-drug visits and path B has 5 and 4, so the contrast variance
   differs by path and the pooled $t$ is approximate. The observed
   nominal size shows the approximation is harmless here; say so.
4. Section 2.2: with $N$ odd the median split gives unequal strata
   (13 and 12 at $N = 25$); the $2(n-1)$ df statement assumes balance.
5. Table 2: the contrast arm's null rate at $N = 70$ is 0.0645
   ($z = 2.98$, $p = 0.004$ by binomial test, VERIFIED). With 20 null
   cells one such excursion is not alarming, but the text calls the
   arm "nominal" without comment. Note it.
6. Results: "its model-based SE is 9 to 12 percent smaller than the
   empirical SE" for uncorrected GEE mixes null and alternative cells;
   at $\beta_{bm} = 0.45$ the deficits are 7 and 9 percent at $N = 25$
   and 35 (VERIFIED). Give the cell.
7. STYLE.md compliance (grep run, VERIFIED): "machinery" (three
   times), "the price of admission", "bought at the price", "the
   honest decomposition", "the honest statement", "Unfortunately",
   "We shall" (twice), and "Architecture-A"/"Architecture-B" (retired
   labels; NOTATION.md Part 3). Six of 19 Abstract sentences exceed 30
   words (maximum 48). Semicolons join clauses throughout. There is no
   Glossary section, which STYLE.md requires. `perl
   tools/notation-lint.pl` is clean (VERIFIED).
8. Cross-references: `docs/27-cycle-period-design-sweep-plan.md` does
   not exist anywhere in the repository (VERIFIED by `find`);
   `docs/06-ar1-residual-correlation.tex` now lives at
   `docs/legacy/06-...`; `docs/26` exists. Fix the paths or drop them.
9. Section 2.3 cites the `lmer` figures from `docs/06` as 13 to 17
   percent "for the crossover (CO) and open-label (OL) designs". The
   source gives CO 0.17 and 0.13, OL 0.13 and 0.09 at $N = 35$ and 70
   (VERIFIED); "9 to 17 percent" is the accurate range.
10. Front matter is stale: `cover-letter.md` describes four
    procedures, $N \in \{25, 35\}$, 800 replicates, a $2.4\times$
    inflation and $c_{bm}$; `title-page.md` gives "800 reps x 16
    cells" in the data statement and British spellings
    ("generalised"); `README.md` still names the `implementations/
    simple` substrate and $t_{1/2} = 1$ week. All three contradict the
    manuscript.
11. `claudecode.tex` is still included via `in_header`; CLAUDE.md
    records it as an inert leftover. Remove the include before
    submission.
12. Bibliography hygiene: `morris2019simulation` and
    `mancl2001covariance` lack DOIs (10.1002/sim.8086 and
    10.1111/j.0006-341X.2001.00126.x); `hedges2022power` still carries
    a 2022 key against `year = {2023}` (prior minor 7, not addressed;
    PubMed gives online 2022-10-12, print 55(7) 2023);
    `tong2024stepwedge` is print volume 22(1), 2025, with the online
    date 2024-10-08 (VERIFIED), so the year and note should say so;
    `pmsimstats-paper10` is an unpublished self-citation and should
    give the paper's working title and directory; `guyatt1986`,
    `lillie2011`, `jones2014crossover`, `levin2012nonparametric`,
    `grajeda2016` and `schlain1993car1` are in the bib but never cited
    (VERIFIED by key grep).
13. The prazosin/PTSD biomarker is described as resting SBP in Section
    2.4 but the DGP draws `bm` with mean 0 and SD 1 and the paper never
    gives units; either state that the biomarker is standardized or
    give the clinical scale.
14. "$T = 8$" is used for the number of post-baseline visits, while
    Section 2 uses $T$ for within-phase timepoints. Use distinct
    symbols.

## 6. Missing and questionable references

| Claim location | Issue | Suggested source |
|---|---|---|
| Sec. 1.3, 5, 8: cost of dichotomizing the biomarker presented as a finding | Classical result, uncited | Cohen J (1983) The cost of dichotomization. *Appl Psychol Meas* 7:249-253; MacCallum RC et al. (2002) *Psychol Methods* 7:19-40; Royston P, Altman DG, Sauerbrei W (2006) *Stat Med* 25:127-141, doi:10.1002/sim.2331 (PubMed, VERIFIED); Altman DG, Royston P (2006) *BMJ* 332:1080 |
| Sec. 1.3: `@pan2001` cited for bias-corrected small-cluster GEE | Wrong paper (Pan 2001 is on sample size for correlated binary data) | Pan W, Wall MM (2002) Small-sample adjustments in using the sandwich variance estimator in GEE. *Stat Med* 21:1429-1441, doi:10.1002/sim.1142 (PubMed, VERIFIED) |
| Sec. 4, 5: Mancl-DeRouen with normal reference; residual over-rejection unexplained | Companion corrections and the $t$-reference evidence absent | Kauermann G, Carroll RJ (2001) *JASA* 96:1387-1396; Fay MP, Graubard BI (2001) *Biometrics* 57:1198-1206; Li P, Redden DT (2015) *Stat Med* 34:281-296, doi:10.1002/sim.6344 (PubMed, VERIFIED) |
| Sec. 5, 6: Kenward-Roger/Satterthwaite recommended to fix conservatism | Direction wrong (M1); the relevant remedy is CR2 | Pustejovsky JE, Tipton E (2018) Small-sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. *J Bus Econ Stat* 36:672-683; the compendium's paper 02 (G6, G9) and paper 10 |
| Sec. 1.3, 1.4: analysis of sets of N-of-1 trials | The mixed-model versus randomization-based distinction for aggregated N-of-1 analysis is not engaged | Araujo A, Julious S, Senn S (2016) Understanding variation in sets of N-of-1 trials. *PLoS ONE* 11:e0167167, doi:10.1371/journal.pone.0167167 (PubMed, VERIFIED) |
| Sec. 2.3: `docs/06` figures | Path moved to `docs/legacy/`; range misquoted (minor 9) | Correct path and range |
| Sec. 1.5, 3, 8: `docs/27` | File does not exist | Remove or restore |

Citations spot-checked and found to support their claims: Tong et al.
(the 6.6 percent figure is exact; PubMed, VERIFIED); Natesan Batley et
al. 2023 (verified in the prior round); Hedges, Shadish and Natesan
Batley (title and venue correct; PubMed, VERIFIED); Haverkamp and
Beauducel 2017 and Blanca et al. 2023 (the gloss is now balanced and
the 0.75 threshold is used consistently; prior minors 2 and 3
addressed, INSPECTED); Ma et al. 2012 (now cited; prior table row
addressed); Blackston et al. 2019 (now cited); Mancl and DeRouen 2001
(the inline $(I - H_{ii})^{-1} r_i$ algebra matches the estimator,
INSPECTED; the `geesmv` agreement is VERIFIED).

## 7. Suggestions for strengthening

1. Add the CR2 Satterthwaite `lme` arm (M1, M3). It is the single
   computation that would connect this paper to papers 02 and 10 and
   give the recommendation table a defensible default. The stored
   replicate file does not contain the fitted objects, so this needs a
   rerun of the `lme` cells only (about 12,000 fits, under ten
   minutes at the stated 7 cores).
2. Rerun the Mancl-DeRouen reference distribution as $t_{K-p}$ from
   the stored `beta` and `se` (M2); no new fits are needed.
3. Report $\kappa$ (model SE over empirical SE) for every arm and
   every cell in Table 4 rather than two sample sizes, since it is the
   quantity that explains all four arms' size behavior.
4. Reframe the dichotomization result as a confirmation of Cohen's
   loss in the N-of-1 interaction setting and compare the observed
   power ratio with the theoretical efficiency factor $2/\pi \approx
   0.64$ for a median split of a normal moderator. Agreement or
   disagreement with that factor is a more informative statement than
   "0.18 of power".
5. Add one carryover cell ($t_{1/2} = 0.5$ week, the series' standard
   grid) so that the "no carryover model" limitation of Section 2.4 is
   exercised rather than asserted, or drop that limitation from the
   list of what this study shows.
6. Either deliver a minimal design slice (prior suggestion 2) or cut
   Section 3 and every reference to a "joint design-by-test" study to
   a single future-work sentence.

## 8. Scope of this review

- **Verified (executed or recomputed):** every number in Tables 2, 3
  and 4 and in the Results prose against `08-test-summary.txt`; the
  replicate counts, 100 percent convergence and 60,000-fit total
  against `08-test-replicates.rds`; exact binomial tests of every
  null cell against 0.05; the Mancl-DeRouen benchmark mean, SD and
  range; positive definiteness of both paths' covariance matrices at
  the driver's parameters under the current `stopOnNonPD = TRUE`
  default, and that `c.bm` does not enter $\Sigma$ under mean
  moderation; `generateData()` runs under the current default with
  the driver's arguments; the Section 2 algebra (sums of squares,
  $t^2 = F$, variance, non-centrality) by re-derivation; the STYLE.md
  grep, the Abstract sentence-length check, the notation linter, the
  citation-key usage grep, and the existence of `docs/06`, `docs/26`
  and `docs/27`; the `docs/06` Type I figures; the timestamps of the
  stored artifacts (consistent with the stated 587 s run).
- **Inspected (read and judged):** the driver
  `01-study1-test-procedure.R` in full; `R/lme_analysis.R` and the PD
  branch of `R/generateData.R`; the working protocol; paper 02's
  Abstract and analysis-code table; CLAUDE.md, STYLE.md, NOTATION.md;
  the README, cover letter, title page and the 2026-06-17 report.
- **Not checked:** I did not rerun any simulation cell; I did not fit
  the CR2 arm I ask for; I did not verify the `nlme` degrees of
  freedom numerically for these data (the claim in M1 rests on the
  `nlme` inner-outer rule); I did not inspect the rendered PDF; I did
  not verify the Cohen (1983) and MacCallum (2002) bibliographic
  details against a database (they are from memory and should be
  checked); I did not locate the source of the Section 2.4 production-
  data figures.
- **Searches run:** PubMed (no user identifiers passed): Tong et al.
  stepped-wedge review; Pan and Wall 2002; Li and Redden 2015;
  Royston, Altman and Sauerbrei 2006; Hedges, Shadish and Natesan
  Batley; Araujo, Julious and Senn 2016; Kenward-Roger/Satterthwaite
  small-sample comparisons since 2015; aggregated N-of-1 mixed-model
  simulations since 2016 (the last two returned little of direct
  relevance). Web: Kauermann and Carroll 2001 and Fay and Graubard
  2001 bibliographic details; degrees-of-freedom practice for the
  Mancl-DeRouen estimator.
