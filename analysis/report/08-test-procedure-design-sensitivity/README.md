# 08-test-procedure-design-sensitivity

*2026-10-28 18:30 PDT*

## Working title

Test-procedure choice for the biomarker-treatment interaction in
aggregated N-of-1 trials: dichotomization, classical RM-ANOVA,
mixed-effects, and generalized estimating equations.

## Origin

The manuscript grew from two project notes:

- `docs/legacy/26-rm-anova-interaction-test.md` derives the
  closed-form strict repeated-measures ANOVA test of the
  biomarker-treatment interaction in a balanced split-plot design
  and relates it to the linear-mixed Wald test.
- A planned cycle-by-period design grid (Studies 2 and 3 of the
  working protocol). The grid was not executed and the paper no
  longer promises it; it is recorded as future work.

## Scope (as executed)

1. Derive the strict RM-ANOVA $F$-statistic for the
   biomarker-stratified treatment interaction and show that it is a
   two-sample $t$-test on within-participant phase differences,
   exact for any within-participant covariance once each phase is
   averaged.
2. Separate the power cost of dichotomizing the biomarker from the
   cost of the test procedure, using a continuous-biomarker
   within-participant contrast that differs from strict RM-ANOVA
   only in keeping the biomarker continuous.
3. Tabulate size, power, coverage and the calibration ratio
   $\kappa$ for the continuous-biomarker procedures.

## Study as executed

- **Procedures (stored run, 2026-10-08).** Strict RM-ANOVA on a
  median-split biomarker (pre-averaged $2 \times 2$); the
  continuous-biomarker contrast (E9b); LME with `corCAR1` and
  model-based inference (`lme_analysis()`, G1 formula plus `De`);
  GEE with the Liang-Zeger sandwich (`geepack`, `corstr = 'ar1'`);
  GEE with the inline Mancl-DeRouen sandwich. The Mancl-DeRouen
  arm's $t_{K-5}$ reference was recomputed from the stored
  estimates (`02-mancl-t-reference.R`).
- **Pending arm.** LME plus CR2 (clubSandwich, Satterthwaite df),
  implemented in the driver as `lme_cr2` and smoke-tested; needs a
  rerun (about 11 minutes on 7 local cores).
- **Grid.** $\beta_{bm} \in \{0, 0.45\}$, $N \in \{25, 35, 70,
  100\}$ (total, split across two paths); 2000 replicates per null
  cell and 1000 per power cell.
- **Design.** Two-path OL+BDC: open-label visits at weeks 4, 8, 12,
  16; blinded weekly visits at weeks 17 to 20; discontinuation after
  week 18 (path A) or week 17 (path B). Not the four-path Hybrid
  design.
- **Substrate.** The pmsimstats package (`generateData()` under
  `dgp_architecture = 'mean_moderation'`, `lme_analysis()`), with
  AR(1) $\rho = 0.7$ within components, cross-component correlations
  0.1 and 0.05, no carryover, no dropout. This is neither the
  original DGP nor the extended DGP ($c_1 = 0.2$) of NOTATION.md.
  The `implementations/simple` substrate named in the working
  protocol was not used.

## Relationship to other papers in the series

- `01-dgp-mean-moderation-vs-mvn`: shares the mean-shift encoding of
  the interaction.
- `02-carryover-sensitivity`: defines the G1 to G9 analysis codes
  and recommends G9 (CR2 inference) for a trial with carryover. The
  present paper's LME arm is G1's formula with model-based
  inference; the pending CR2 arm is G6's inference.
- `10-interaction-test-calibration`: diagnoses the model-based LME's
  conservatism as working-covariance misspecification and remedies
  it with CR2; the present paper's $\kappa$ table shows the same
  defect.
- `05-nof1-design-sensitivity`: design-parameter sensitivity for the
  main effect rather than the interaction.

## Files

- `report.Rmd`: the manuscript. Render with
  `bash tools/render.sh analysis/report/08-test-procedure-design-sensitivity/report.Rmd`.
- `references.bib`, `statistics-in-medicine.csl`.
- `cover-letter.md`, `title-page.md`.
- `referee-report-2026-06-17.md`, `referee-report-2026-10-28.md`.
- `report-slim.Rmd`, `report-mini.Rmd`, `bullets.Rmd`: earlier
  condensed variants, not maintained with the main manuscript.
- Driver: `analysis/scripts/test-procedure-design-sensitivity/`
  (`00-ademp-protocol.md`, `01-study1-test-procedure.R`,
  `02-mancl-t-reference.R`).
- Data: `analysis/data/quick-sim/08-test-*.txt|rds`.

## Author

pmsimstats team
