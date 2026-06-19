# Referee Report: Test-procedure and trial-design choices for the biomarker-treatment interaction in aggregated N-of-1 trials
*2026-06-17 17:09 PDT*

Manuscript: `analysis/report/08-test-procedure-design-sensitivity/report.Rmd`
Target orientation: *Statistics in Medicine* (methodological, applied
biostatistics). Reviewed as for a top-tier methods venue.

---

## 1. Summary of the manuscript

The paper addresses inference about the biomarker-by-treatment
interaction in aggregated N-of-1 trials, framed around two analyst
choices: the test procedure (strict repeated-measures ANOVA, the
linear mixed-effects model with continuous-time AR(1) residual
correlation, and generalised estimating equations with naive and with
Mancl-DeRouen small-sample sandwich variance) and the trial-design
parameters (cycle count, on/off period lengths, symmetry, run-in/out).
Its stated contributions are fourfold: (i) a closed-form strict
RM-ANOVA F-statistic for the biomarker-stratified treatment
interaction in a balanced split-plot design, and its relation to the
linear-mixed Wald test; (ii) specification of a cycle-by-period design
grid calibrated to a prazosin/PTSD reference application; (iii) a
Monte Carlo study, reported under the Morris et al. (2019) ADEMP
framework, crossing the test procedures against the
biomarker-moderation coefficient and sample size; and (iv) a
design-and-test recommendation table. The empirical work delivered is
a 16-cell factorial (4 procedures x c_bm in {0, 0.45} x N in {25, 35},
800 replicates per cell). The cycle-by-period sweep is specified but
deferred to a companion paper.

## 2. Overall assessment and recommendation

**Recommendation: Major revision** (and, absent the deferred design
sweep, the contribution presently fits a specialist methods journal
rather than a top-tier general venue).

The closed-form derivation in Section 2 is correct and cleanly
presented; I re-derived it and it holds (see Section 4, comment M0).
The simulation is competently engineered and reproducible. However,
the paper has five problems that are individually serious and jointly
disqualifying at the current draft: (a) the headline test-procedure
comparison is confounded with an analyst preprocessing choice
(forced median dichotomisation of the biomarker), so the reported
"linear-mixed dominates RM-ANOVA" gap does not isolate the test
procedure; (b) the theoretical motivation built around sphericity
violation is disconnected from the mechanism actually operating in the
simulation, because the simulated RM-ANOVA pre-averages to one value
per phase and therefore never exercises the sphericity problem; (c)
the executed design deviates from the paper's own pre-registration on
nearly every factor level, undocumented; (d) the title and abstract
promise a joint design-by-test study that the paper does not deliver;
and (e) secondary estimands that the Methods section pre-specifies
(bias, CI coverage) are never computed or reported. None of these is
fatal to the underlying programme, but all must be resolved before the
claims are defensible.

## 3. Significance and novelty

The genuine and defensible novelty is narrow but real: the
specialisation of the test-procedure comparison to the
biomarker-by-treatment *interaction* estimand in the aggregated N-of-1
/ Hendrickson hybrid setting, with continuous-time AR(1) residual
correlation and an inline Mancl-DeRouen small-sample sandwich. The
closed-form strict RM-ANOVA interaction F and its reduction to a
two-sample t-test on within-subject treatment differences is a tidy
specialisation, though it is textbook split-plot algebra
(Maxwell-Delaney-Kelley) rather than new theory.

The broader claim, asserted repeatedly (Section 1.3, Section 1.4),
that "no test-procedure comparison addresses the cross-over or N-of-1
setting directly" and that the field lacks a procedure comparison, is
overstated. The generic comparison of RM-ANOVA versus GEE versus mixed
models for longitudinal data is well established, with quantitative
findings strikingly close to this paper's own. In particular:

- Ma, Mazumdar & Memtsoudis (2012), *Beyond repeated measures ANOVA*,
  *Reg Anesth Pain Med* 37:99-105, compares RM-ANOVA, GEE, and mixed
  models on power and Type I error for longitudinal endpoints and
  reports RM-ANOVA at roughly 30% and 50% lower power than GEE and
  mixed models respectively. This is the same ordering and nearly the
  same magnitude the present paper reports as a finding. It is not
  cited.
- The aggregated-N-of-1 simulation literature (e.g. the
  parallel/crossover/aggregated-N-of-1 simulation comparison in
  *PLoS ONE* / PMC6955665, 2020) bears directly on the design half of
  the framing and is not engaged.

The contribution would clear a specialist-methods bar once the
confound (M2) is fixed and the comparison genuinely isolates the test
procedure. It does not clear a top-tier general-statistics bar as
written, principally because the empirical content reduces to a
16-cell confirmatory exercise whose central comparison is not clean,
and because the joint design-by-test question, which is the paper's
stated reason to exist, is deferred in its entirety.

## 4. Major comments

**M0 (verification, positive).** I re-derived the Section 2 closed
form. The 2x2 reduction SS_AB = nT L-hat^2 / 4 (eq. 2.7-2.9), the
identity t_AB^2 = F_AB with Var(L-hat) = 4 sigma_e^2 / (nT)
(eqs. 2.11-2.14), and the non-centrality lambda = nT delta^2 /
(4 sigma_e^2) (eq. 2.15) are internally consistent and correct under
sphericity and balance. The numerical claims in Section 5 are also
internally consistent: the N=25 power gap 0.533 - 0.340 = 0.193 at
~8 SE of the difference, the N=35 gap 0.284 at ~12 SE, and the GEE
N-level gap 0.043 at ~2.9 SE all reproduce on recomputation. No action
needed; stated so the authors know these were checked.

**M1 (correctness; motivation disconnected from evidence).** The
sphericity argument is the spine of Sections 1.2 and 2.3: chronic-
condition autocorrelation violates sphericity, strict RM-ANOVA inflates
Type I error, the mixed model with `corCAR1` avoids the
degrees-of-freedom penalty. But the simulated RM-ANOVA procedure does
not exercise sphericity at all. Inspection of the driver
(`08-test-procedure-prototype.R`, lines 157-191) shows it median-splits
the biomarker and then **pre-averages to a single value per subject per
phase** (`agg <- long2[, .(Sx = mean(Sx, ...)), by = .(ptID, bm_high,
Db_f)]`) before `aov(Sx ~ bm_high * Db_f + Error(ptID/Db_f))`. With one
value per phase, the within-subject covariance is 2x2 and sphericity is
automatic; there is no epsilon < 1 to correct. Consistently, the
simulated RM-ANOVA holds Type I error at nominal (5.25%, 5.00%;
Section 5 table), exactly the opposite of the 13-17% inflation the
motivation predicts. The power deficit RM-ANOVA suffers in the
simulation is therefore driven by information loss (dichotomisation
plus pre-averaging), not by sphericity. The paper cannot motivate the
result with a mechanism the experiment switches off. *Remedy:* either
(a) reconcile the two by analysing the un-averaged repeated measures
with strict RM-ANOVA so sphericity actually bites and is corrected via
Greenhouse-Geisser, and report that arm; or (b) drop the sphericity
framing as the explanation for the simulated gap and rewrite Sections
1.2/2.3 so the stated mechanism (dichotomisation + averaging) matches
what the code does.

**M2 (correctness; confounded comparison).** The four procedures do not
share an estimand, and the paper concedes this (Section 5,
"related but not identical estimands"). RM-ANOVA is fed a
median-dichotomised biomarker; LME and GEE are fed the continuous
biomarker. By the paper's own Section 2.4 accounting, dichotomisation
alone inflates the interaction SE by 1.4-1.7x and costs 30-50% of
power. That is the same order as the entire reported LME-over-RM-ANOVA
gap. The headline claim that "the linear-mixed analysis dominates
strict RM-ANOVA" (Abstract; Section 6) therefore largely measures the
cost of an analyst preprocessing decision, not a property of the
F-test machinery. Calling this "the more honest comparison"
(Section 5) does not rescue it: an honest *procedure* comparison must
hold the data representation fixed. *Remedy:* add an arm that isolates
the test procedure from the preprocessing, e.g. a continuous-biomarker
within-subject contrast t-test (the exact analogue of eq. 2.11 without
dichotomisation), or feed all four procedures the dichotomised
biomarker. Only then can the power gap be attributed to the test.

**M3 (correctness; misattributed Type I figure).** Section 2.3
(orig block, around line 1092) states "the strict RM-ANOVA Type I
error sits at 13 to 17 percent for the CO and OL designs absent
correction," citing the audit `docs/06-ar1-residual-correlation.tex`.
Per that audit, the 13-17% figures are the Type I inflation of a
**misspecified `lmer` model with no residual correlation structure**
(0.17 at N=35, 0.13 at N=70 for CO; 0.13/0.09 for OL), which is then
fixed by switching to `nlme::lme` with `corCAR1`. Those numbers are
for a naive mixed model, not for strict RM-ANOVA. As written the
sentence attributes one procedure's failure to a different procedure.
This is inspected via the source document rather than re-executed, so I
flag it as a strong discrepancy to verify, not a certainty; but if it
holds it is a factual error in the central motivation. *Remedy:* verify
against `docs/06`, correct the attribution, and report the actual
strict-RM-ANOVA Type I behaviour (which, per M1, is nominal in the
pre-averaged form).

**M4 (methodological integrity; undocumented pre-registration
deviations).** The paper foregrounds ADEMP discipline and a
pre-registered plan (`00-ademp-pre-registration.md`). The executed
Study 1 departs from that plan on essentially every level: procedures
3 -> 4 (a naive-sandwich GEE arm added); c_bm {0, 0.3, 0.6} ->
{0, 0.45}; N {35, 70, 100} -> {25, 35}; replicates 1000 (5000 for Type
I) -> 800. Notably, N = 25, the sample size at which the marquee GEE
inflation (2.4x) and the largest framing claims land, was not in the
pre-registered grid, while N = 70 and 100 were dropped. A
pre-registration whose every cell is changed, with the headline
condition added post hoc and the deviations log referenced but not
shown, undercuts the ADEMP claim rather than supporting it. *Remedy:*
include an explicit deviations table (pre-registered vs executed, with
rationale for each change), and either justify the N restriction to
{25, 35} or restore the pre-registered N levels.

**M5 (scope; title and abstract overpromise).** The title names
"the cycle-by-period protocol grid" and the abstract frames a *joint*
characterisation of "both axes." The paper delivers neither: the design
grid is specified (Section 3) but swept nowhere, and the recommendation
(Section 6) is narrowed to test procedure at a single fixed design. The
body is candid about this (Section 3, Section 6), but a title and
abstract that promise the joint study while the contribution is a
single-design procedure comparison will mislead readers and reviewers.
*Remedy:* retitle to reflect the test-procedure scope (drop or
subordinate "cycle-by-period protocol grid"); rewrite the abstract so
the joint study is described as motivation/future work, not as
delivered method and result.

**M6 (ADEMP completeness; unreported estimands).** Section 4
pre-specifies secondary estimands: "bias of beta_bm:D ..., nominal-95%
CI coverage, and convergence rate." Code inspection confirms only the
rejection rate, its MCSE, and convergence are computed and saved
(`08-test-procedure-prototype.R`, lines 431-436, 479-482); bias and CI
coverage are never produced (the saved `beta` is unused, and RM-ANOVA
returns `NA` for beta). A paper invoking the Morris et al. framework as
its discipline should not define performance measures it does not
report. *Remedy:* compute and tabulate bias and coverage for the LME
and GEE arms, or remove them from the estimand list and explain the
omission.

**M7 (reproducibility statement inaccurate).** The Reproducibility
section states "Per-replicate seeds derive from the master seed by
+ rep_idx" and that the driver "uses `parallel::mc.set.seed` inside
`mclapply`." Code inspection found a single `set.seed(20260507)` and no
visible per-replicate `+ rep_idx` or `mc.set.seed` stream management. A
reproducibility statement that does not match the code is worse than
none. *Remedy:* either implement the described per-replicate seeding
(needed for exact reproducibility under `mclapply`, which otherwise
does not guarantee a fixed stream per worker) or correct the statement
to describe what the code actually does.

## 5. Minor comments

1. **Three vs four procedures (internal inconsistency).** Section 1
   contribution 3 (orig block) says "three test procedures (strict
   RM-ANOVA, ... linear-mixed ..., GEE with bias-corrected sandwich
   variance per @wang2016geesmv)," while the Abstract, Section 4, and
   Section 5 describe four (naive + Mancl-DeRouen, cited
   @mancl2001covariance). Reconcile the count and the citation: the
   executed work uses Mancl-DeRouen (2001), implemented inline;
   `wang2016geesmv` is the package (`geesmv`) that wraps it, not the
   estimator's source. State both correctly.

2. **Blanca et al. (2023) epsilon threshold.** Section 1.2 and Section
   2.3 attribute "Greenhouse-Geisser when epsilon < 0.6" to Blanca et
   al. (2023). The conventional Girden threshold, and the one the paper
   itself uses at line 1039 ("Huynh-Feldt ... recommended when epsilon
   > 0.75"), is 0.75. The 0.6 figure appears to be a misquote; verify
   against the source and reconcile the two thresholds the paper cites.

3. **Haverkamp & Beauducel (2017) over-read.** The paper cites this as
   showing RM-ANOVA "inflates Type I error rapidly" as occasions
   increase, implying mixed models are uniformly safer. The source is
   more balanced: corrected RM-ANOVA and MLM perform comparably, and
   MLM is in some conditions more liberal. Soften the gloss to match.

4. **"27 cells" vs 16 cells.** Section 4 (orig block, ~line 1418)
   refers to "Study 1: test-procedure comparison at fixed design,
   27 cells," but the executed Study 1 is 16 cells. Reconcile in text
   (this is the same pre-registration deviation as M4 and should be
   cross-referenced).

5. **Manuscript not in final narrative form.** Every section carries a
   `.rgt` block of Lorem ipsum placeholder text alongside duplicated
   `bullets` and `orig` prose. The reviewable content is coherent, but
   the paper cannot be submitted in this tri-partite scaffold; the
   final single-narrative version should be the object of the next
   review round.

6. **Kenward-Roger recommended but uncited.** Sections 5 and 6
   recommend a Kenward-Roger or Satterthwaite correction for the LME
   conservatism but cite neither Kenward & Roger (1997) nor Satterthwaite.
   Add the references for any method the recommendation table names.

7. **Bibliography hygiene.** `hedges2022power` carries a 2022 cite key
   but `year = {2023}` in the entry (Behavior Research Methods 55(7));
   align the key with the year, or note the 2022 online/2023 print
   distinction. Otherwise the bib is clean and DOIs are present.

8. **MCSE resolution stated correctly.** The text's handling of the
   N=35 Mancl-DeRouen rate (0.045, "covering the nominal") respects the
   ~0.008 MCSE; no overclaiming there. Good.

9. **Geometry.** `right=5cm` margin is unusually wide; presumably for
   annotation during review. Cosmetic only.

## 6. Missing and questionable references

| Claim location | Issue | Suggested source |
|---|---|---|
| Sec. 1.3-1.4: "no test-procedure comparison addresses ... the N-of-1 setting"; the RM-ANOVA/GEE/MEM power ordering reported as new | Generic RM-ANOVA vs GEE vs mixed-model power/Type-I comparison is established prior art with near-identical magnitudes; not engaged | Ma Y, Mazumdar M, Memtsoudis SG (2012), *Beyond repeated measures ANOVA: advanced statistical methods for the analysis of longitudinal data in anesthesia research.* Reg Anesth Pain Med 37(1):99-105. doi:10.1097/AAP.0b013e31823ebc74 |
| Sec. 1.1, 1.4: aggregated-N-of-1 design framing and novelty | Aggregated-N-of-1 simulation comparisons exist and bear on the design axis | Aggregated N-of-1 vs parallel/crossover simulation study, *PLoS ONE* 2020 (PMC6955665); verify full citation and cite |
| Sec. 2.3 (~line 1092): "strict RM-ANOVA Type I error ... 13 to 17 percent" | Per `docs/06`, those figures are for a misspecified `lmer` (no corstr), not strict RM-ANOVA | Correct attribution; cite the internal audit accurately or provide an in-paper RM-ANOVA Type I figure |
| Sec. 1.2 / 2.3: GG threshold "epsilon < 0.6 per Blanca et al. (2023)" | Threshold likely 0.75, not 0.6 | Re-verify Blanca et al. (2023) doi:10.3389/fpsyg.2023.1192453; reconcile with the 0.75 used at line 1039 |
| Sec. 5-6: Kenward-Roger / Satterthwaite recommendation | Method recommended without citation | Kenward MG, Roger JH (1997), Biometrics 53:983-997; Satterthwaite FE (1946), Biometrics Bull 2:110-114 |

Citations I spot-checked and found accurately supporting their claims:
Natesan Batley et al. (2023) "83.8% ignore autocorrelation" (verified
against the source, exact); Hendrickson et al. (2020); Morris et al.
(2019); Mancl & DeRouen (2001) (the inline `(I - H_ii)^{-1} r_i`
algebra matches the estimator). Bibliography entries are complete and
correctly venued.

## 7. Suggestions for strengthening the paper

1. **Make the comparison clean (highest value).** Add the
   continuous-biomarker within-subject contrast arm (M2). If, with the
   data representation held fixed, LME still dominates, the paper's
   central claim becomes defensible and considerably stronger; if the
   gap collapses, that is itself the important finding.

2. **Deliver at least a slice of the joint study.** The paper's reason
   to exist is the design-by-test interaction. Even a reduced
   cycle-by-period sweep (e.g. k in {2, 4} x one period contrast) at
   the two strongest procedures would convert the title from a promise
   to a result and lift the contribution toward a general venue.

3. **Restore the pre-registered N range.** Reporting N in {25, 35, 70,
   100} would show whether the GEE inflation and the LME-RM-ANOVA gap
   behave monotonically in N, and would remove the appearance of having
   selected the two smallest N to maximise the effect.

4. **Externally validate the inline Mancl-DeRouen sandwich.** Benchmark
   against `geesmv::GEE.var.md` on a few cells (the paper already flags
   this as a follow-up); doing it removes a stated caveat cheaply.

5. **Report the full ADEMP table.** Bias, empirical SE, model SE, and
   coverage per cell (M6), with MCSEs, is the expected deliverable for
   a Morris-framework paper and would let readers see the
   conservatism/efficiency trade-off directly rather than inferring it.

6. **Raise replicate count for the Type I cells.** The pre-registration
   itself called for 5000 replicates at the null; 800 gives MCSE ~0.008
   at p=0.05, marginal for the 0.045-vs-0.05 distinctions the paper
   draws. The null cells are cheap; run them deeper.

## 8. Scope of this review

- **Verified (re-derived or recomputed):** the Section 2 closed-form
  derivation (SS_AB reduction, t^2 = F identity, variance, NCP); the
  internal arithmetic of the Section 5 power gaps and SE-of-difference
  multipliers; the cell-count reconciliation (4 cells x 800 x 4
  procedures = 12,800 fits).
- **Inspected (read source/code, judged):** the driver
  `08-test-procedure-prototype.R` (RM-ANOVA dichotomise + pre-average;
  LME `corCAR1` formula; inline Mancl-DeRouen `(I - H_ii)^{-1} r_i`;
  parameters c_bm {0,0.45}, N {25,35}, 800 reps, seed 20260507;
  rejection rate as proportion p<0.05; secondary estimands not
  computed; AR(1) GEE working correlation). The pre-registration
  (`00-ademp-pre-registration.md`) and the audits `docs/26`,
  `docs/27`, `docs/06` were read via a sub-agent; the M3 and M4
  discrepancies rest on that reading and are flagged for the authors to
  confirm against the primary files.
- **Not checked:** I did not execute the simulation or reproduce the
  reported rates from `08-test-replicates.rds`; the per-replicate seed
  mechanism (M7) was not run; the absolute calibration of the inline
  Mancl-DeRouen magnitude was not benchmarked; the rendered PDF layout
  was not inspected.
- **Literature searches run:** Blanca et al. (2023) sphericity
  threshold; Natesan Batley et al. (2023) autocorrelation statistic
  (confirmed exact); Haverkamp & Beauducel (2017) actual conclusions;
  prior RM-ANOVA/GEE/mixed-model procedure comparisons for longitudinal
  data (located Ma et al. 2012 and an aggregated-N-of-1 simulation
  study bearing on novelty). I did not run a PubMed systematic sweep of
  the N-of-1 design-optimisation literature beyond these targeted
  queries.

---
*Rendered on 2026-06-17 at 17:09 PDT.*<br>
*Source: ~/prj/res/36-pmsimstats-ng/pmsimstats-ng/analysis/report/08-test-procedure-design-sensitivity/referee-report-2026-06-17.md*
