# Analyzing a real trial

*2026-10-15 09:59 PDT*

**Author:** pmsimstats team  
**Purpose:** analyze_trial_extended(); characterize_carryover(); CENT
2015 reporting.  
**Prerequisites:** ch22, ch23, ch24

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

*Local symbols used in this chapter only. For project-wide notation see
`analysis/report/NOTATION.md`.*

<!-- notation entries -->

---

## Summary

<!-- evidence-labeled results go here -->

---

## Construct box

<!-- table: process, carryover, coupling, allocation, design, N,
analysis ID, t_half_a, n_reps -->

---

## Setting and argument

<!-- SOURCES FROM docs/59 Section 3 for this chapter:
- Sources: 14 via Appendix I Section 9; vignette 03; 56 Section 8;
  57 Section 5; 58 Section 10; 23 Section 10.

- New to write: most of the chapter.
-->

### Package entry points for a real trial

<!-- From docs/legacy/14-trial-analysis-implementation-guide.tex
lines 1-120 -->

Two functions export the analysis pipeline to applied users:

**`analyze_trial_extended(td, dat, op, threshold = 5)`** runs the
full F1 analysis on a single real or simulated dataset and returns a
summary list including the LME model object, the coefficient table,
the biomarker threshold on the fitted biomarker scale, and the
predicted response trajectories. The `threshold` argument sets the
minimum fitted response improvement used to classify a participant as
a likely responder.

**`characterize_carryover(td, dat, half_lives = c(0.25, 0.5, 1, 2,
4))`** fits the analysis model at each candidate half-life, extracts
the AIC, and returns the AIC-preferred value as `cc$best_t_half`.
This is a descriptive procedure; `docs/57` and `docs/58` document why
testing at the AIC-preferred half-life without calibration is liberal.

Both are demonstrated in `vignettes/03-analyze-clinical-trial-data.Rmd`.

### Worked analysis of the included CTdata

<!-- From vignettes/03-analyze-clinical-trial-data.Rmd lines 1-137 -->

The vignette uses the package's included `CTdata` dataset (synthetic
N-of-1 trial with $N = 12$ participants, Hybrid design). The analysis
proceeds in three steps.

**Step 1: define the trial design.**

```r
td <- buildtrialdesign(
  n_reps          = 1,
  design          = 'Hybrid',
  carryover_t1half = 0.5
)
```

The `carryover_t1half` argument sets the assumed half-life for the
analysis model, i.e., the $t_{1/2}^a$ at which `Dbc` is built. Use
`characterize_carryover()` to guide this choice from the data.

**Step 2: set analysis options.**

```r
op <- list(
  dgp_architecture = 'mvn',
  use_dbc          = TRUE,
  use_expectancy   = TRUE
)
```

Setting `use_expectancy = TRUE` adds `De` and `bm:De` to the model
(the F1 specification). Setting `use_dbc = TRUE` uses the
exposure-decayed indicator $D_{bc}$ rather than the binary $D$.

**Step 3: run the analysis.**

```r
result <- lme_analysis(td, CTdata, op, full_model_out = TRUE)
print(result$summary['bm:Dbc', ])
```

The `full_model_out = TRUE` flag retains the `nlme` model object for
downstream diagnostics.

**Predicted trajectories.**

```r
preds <- predict_trajectories(result, CTdata, td)
plot_predicted_response(preds)
```

The trajectory plot shows per-participant fitted symptom curves
overlaid on observed values, with shaded carryover windows.

### Practical recipe for F1 with R code

<!-- From docs/56 Section 8 -->

For a researcher applying F1 to their own trial dataset:

```r
library(nlme)
library(clubSandwich)

# Fit F1 (random intercept + corCAR1 + CR2 inference)
fit <- lme(
  Sx ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De,
  random = ~ 1 | ptID,
  correlation = corCAR1(form = ~ t | ptID),
  data = dat,
  method = "REML"
)

# CR2 + Satterthwaite inference for bmc:Dbc
coef_test(
  fit,
  vcov    = "CR2",
  cluster = dat$ptID,
  coefs   = "bmc:Dbc",
  test    = "Satterthwaite"
)
```

Notes on the recipe:

- `bmc` is the biomarker centered at its sample mean.
- `De` is the design expectancy $\eta_t$: 1 open label, 0.5 blinded,
  0 placebo.
- `Dbc` uses an assumed half-life; derive it via
  `characterize_carryover()` or fix at 0.2 weeks for the validated
  fallback.
- `cluster = dat$ptID` identifies the 50 to 100 participants as
  independent clusters for CR2.

### Null trial worked example: profile likelihood at each half-life

<!-- From docs/57 Section 5 -->

A null trial with $c_{bm} = 0$ ($N = 70$, Hybrid, 500 replicates)
was analyzed at each half-life in the grid $\{0, 0.2, 0.5, 1, 2\}$
weeks. The profile likelihood table (`docs/57` Table 5) shows:

- At $t_{1/2}^a = 0$: null rejection 0.050.
- At $t_{1/2}^a = 0.5$ (correct assumed value): null rejection 0.053.
- At the AIC-selected $t_{1/2}^a$ (F1A, no bootstrap calibration):
  null rejection 0.065 (inspected; pooled estimate, all models liberal
  by at least 0.013).

The AIC selection is anticonservative without bootstrap calibration.
The recommended fallback is to fix $t_{1/2}^a = 0.2$ weeks unless
there is strong evidence for a longer half-life in the specific
compound.

### Drug effect trial worked example

<!-- From docs/58 Section 10 -->

A trial with a true drug main effect ($c_{bm} = 0$, drug response
mean = 5 points, $N = 70$) was analyzed at each grid half-life. The
F1 analysis power at $t_{1/2}^a = 0.2$ weeks was 0.71; at
$t_{1/2}^a = 1$ week it fell to 0.58 (inspected from `docs/58`
Table 10). The visit-factor (VF) model achieved 0.73 and 0.67
respectively, at the cost of size inflation of 0.015 under
`ri-CAR1-KR` inference.

The data-chosen half-life procedure (F1A0B: choose by AIC on a model
without the interaction, then calibrate by a parametric bootstrap) held
size at 0.052 and reached 0.69 power across the grid, close to the
best fixed half-life without requiring the researcher to know the
pharmacology in advance.

### Alternative analysis strategies when carryover is severe

<!-- From docs/legacy/23-documentation-compendium.md Section 10 -->

When carryover is expected to be strong (e.g., behavioral or
psychometric mechanisms), several complementary strategies are
available beyond adjusting $t_{1/2}^a$:

1. **Restrict to on-drug observations only** and compare participants
   at matched on-drug timepoints. This eliminates the off-drug
   information but avoids contamination. It is equivalent to an ANCOVA
   with drug state removed and substantially reduces power.
2. **Weighted analysis** downweighting observations in the estimated
   carryover window, with weights $w_t = D_{bc,t}$ or a threshold
   rule. Protects against carryover contamination in the control arm
   at the cost of effective sample size.
3. **Within-subject contrast** collapsing each participant's data to
   a single on-minus-off value, then regressing on the biomarker.
   This is the two-stage E9 estimate and is discussed in Appendix H.
4. **Two-stage random slopes:** in the first stage, fit a
   mixed model and extract the random slope for the biomarker-by-drug
   interaction; in the second stage, test the mean of the empirical
   Bayes slopes. More complex but accounts for individual variation
   in treatment response.
5. **Design-level solutions:** longer washout periods, separate open-
   label and blinded phases, participant-level crossover randomization.
   These are planning decisions and do not help post-hoc.

<!-- TODO: write the full chapter. The sources above provide anchors
for the practical recipe, the two worked examples (null trial and
drug effect trial), and the alternative strategies section. The
following sections remain to be written:
  - CENT 2015 reporting checklist mapped to the pmsimstats functions
  - How to interpret the interaction estimate for a clinical audience
  - Threshold analysis and responder classification via
    analyze_trial_extended()
  - characterize_carryover() in detail: AIC table, profile plot,
    bootstrapped calibration for F1A0B
  - Sensitivity analysis: vary t_half_a across the grid and report
    all six estimates
  - Example report table formatted for a manuscript Methods section
-->

---

## Derivations

<!-- one row per result: statement | CM label | verification script -->

---

## Numerical evidence

<!-- tables with MCSE; each caption cites script and output file -->

---

## Relation to the manuscripts

<!-- which manuscript sections depend on this chapter; do the numbers
agree? -->

---

## Limitations and open items

<!-- link each to the open-items register at
docs/series/registry/open-items.csv -->

---

## Reproduce

<!-- one command per output -->

---

## References

<!-- from docs/references/nof1-simulation.bib -->

---

## Change log

| Date | Section | Change | Reason |
|------|---------|--------|--------|
| 2026-10-15 | Setting and argument | Initial content move from docs/legacy/14, vignette 03, docs/56, 57, 58, 23 | Phase 2 content moves |
