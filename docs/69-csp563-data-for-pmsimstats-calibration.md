---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# What the VA Prazosin Trial CSP 563 Can Tell the pmsimstats-ng Simulations: Empirical Checks and Candidate Calibrations {.unlisted .unnumbered}
*2026-10-30 11:47 PDT*

**Author.** pmsimstats team

**Purpose.** The pmsimstats-ng simulations rest on a data-generating
process (DGP) whose parameters were set from the published code of
Hendrickson et al. (2020), which in turn calibrated them to an earlier
positive trial of prazosin in active-duty soldiers. This paper asks
which of those parameters can be checked or estimated against
participant-level data from a second, larger trial: the VA Cooperative
Studies Program trial CSP 563 (Raskind et al. 2018, 304 veterans,
negative overall). It reports what was estimated, what could not be,
and what each result implies for the simulations and the manuscripts.

**Source.** All analyses were run in the private compendium
`~/prj/alz/14-murray-nejm/murraynejm`, branch `sbp-caps-replication`,
directory `analysis/report/04-sbp-caps-replication/`. Two analysis
plans were committed before any model in them was fitted (`plan.md`,
`plan-components.md`); a measurement-error follow-up and the parameter
summaries were added afterwards and are labeled as such. Only
aggregate results appear here. The participant-level data stay in that
compendium. The studies are exploratory.

**Evidence labels.** *Computed* means produced by a script in this
session from the CSP 563 data or from the pmsimstats-ng formulas.
*Quoted* means taken from a publication. *Inferred* means argued from
the computed results and not itself tested.

```{=latex}
\clearpage
\tableofcontents
\clearpage
```

## Notation

Symbols follow `analysis/report/NOTATION.md`.

| Symbol | Meaning |
|------------|----------------------------------------------------------------|
| $Y_{it}$ | CAPS-IV total symptom severity of participant $i$ at visit $t$ |
| $BL_i$ | baseline severity; the DGP writes $Y_{it} = BL_i - (TV_{it} + PB_{it} + BR_{it})$ |
| $TV$, $PB$, $BR$ | natural-history, placebo-response and drug-response components |
| $\eta_t$ | design expectancy: 1 open label, 0.5 blinded |
| $B_i$ | biomarker: standing systolic blood pressure (SBP), mm Hg |
| $\sigma_{BR}$, $\sigma_{bm}$ | DGP standard deviations of the drug response (8) and the biomarker |
| $c_{bm}$ | correlation between the biomarker and the drug response (the interaction strength) |
| $c_{bm,BL}$ | correlation between the biomarker and baseline severity |
| $\rho$, $c_1$ | DGP AR(1) correlation per week within a component (0.7) and between components (0.2) |
| $\varphi$, $\nu$ | fitted serial correlation per week and measurement-error share of the CAPS total |
| $\tau$ | SD of individual drug effects (treatment heterogeneity) |
| $a^\top A a$ | time part of the variance of a participant's on-minus-off contrast |

Table: Notation

## Glossary

- **CSP 563.** VA Cooperative Studies Program trial of prazosin for
  PTSD in veterans; 304 randomized, 10-week primary endpoint, 26-week
  follow-up; no benefit on its primary outcomes (Raskind et al. 2018).
- **CAPS-IV total.** Clinician-Administered PTSD Scale for DSM-IV,
  total symptom severity; the outcome scale of the DGP. A secondary
  outcome of CSP 563.
- **Correlation encoding.** The DGP's encoding of the
  biomarker-by-treatment interaction as a correlation between the
  biomarker and the drug-response component.
- **CR2.** Bias-reduced cluster-robust standard error with
  Satterthwaite degrees of freedom.
- **Extended DGP.** The pmsimstats-ng recommended DGP, with separable
  AR(1) correlation ($\rho = 0.7$, $c_1 = 0.2$).
- **Nugget.** Visit-specific measurement error: variance not shared
  with any other visit.
- **Predictive, prognostic.** A biomarker that modifies the drug
  effect, and one associated with the outcome regardless of drug.

## 1. Summary

1. **The blood-pressure biomarker effect of Raskind et al. (2016) does
   not replicate.** In 279 participants the arm-by-SBP interaction on
   CAPS total change was +0.25 points per 10 mm Hg (95% CI −2.72 to
   3.21; computed). The prazosin-arm slope excludes the published −14.
2. **CSP 563 cannot test the effect sizes the simulations use.** On the
   simulation scale the interval is $c_{bm}$ −0.57 to 0.48. The
   project's grid, $c_{bm}$ 0.25 to 0.45, implies only 1.4 to 2.5
   points per 10 mm Hg, inside the interval. The published −14 would
   need $c_{bm} \approx 2.5$, which is impossible for a correlation
   (computed).
3. **The biomarker is not prognostic.** SBP and baseline CAPS correlate
   at $r = 0.02$ (−0.09 to 0.13; computed). The DGP's
   $c_{bm,BL} = 0$ is supported.
4. **CAPS correlation decays slowly, with about 10 percent measurement
   error.** Compound symmetry fits worst. The total-score correlation
   implied by the extended DGP falls too fast at short lags and then
   flattens (computed).
5. **Under that covariance the Hybrid design's advantage over the
   crossover is larger** than under the extended DGP: contrast noise
   0.37 to 0.41 of the crossover's, against 0.59 to 0.69 (computed).
6. **Placebo-arm improvement is about 1.8 times the DGP's, with the
   same time course.** There is no mean drug response, and the SD of
   individual drug effects is at most about 8 points (computed).
7. **Missingness was unrelated to prior worsening**, unlike the DGP's
   biased dropout (computed).
8. **Natural history and placebo response cannot be separated** in a
   fully blinded parallel trial. Neither can the carryover half-life
   or the per-component correlations.

## 2. The data and what the design identifies

CSP 563 randomized 304 veterans with chronic PTSD to prazosin or
placebo, blinded throughout (Raskind et al. 2018, quoted). CAPS total
was recorded at screening and weeks 6, 10, 18 and 26 (303, 273, 263,
249 and 245 records). Standing and supine blood pressure were recorded
at screening, on randomization day and at most later visits.

The design fixes what can be learned (inferred from the design):

- Everyone was blinded at every visit, so $\eta_t = 0.5$ throughout and
  time under expectancy equals calendar time. $TV$ and $PB$ share one
  clock and one weight and are **not separately identifiable**.
- No one stopped drug during the trial, so there is **no carryover
  information**.
- It is a parallel trial: the drug response is identified only as the
  difference between arms, and its variability only through the
  difference in outcome variance between arms.
- The data export has no expectancy or treatment-guess measure and no
  study-drug dose records.

## 3. The biomarker

### 3.1 Predictive effect

Mixed model of CAPS total change at weeks 6 to 26 with arm, SBP and
week and their interactions, adjusted for baseline CAPS, the two
randomization stratifiers and site; unstructured covariance; CR2
inference. SBP is the mean of the screening and randomization-day
readings (computed).

| Quantity, CAPS points per 10 mm Hg | Estimate | 95% CI |
|---|---|---|
| Interaction (prazosin minus placebo slope) | +0.25 | −2.72 to 3.21 |
| Prazosin-arm slope | −0.89 | −3.14 to 1.35 |
| Placebo-arm slope | −1.14 | −3.23 to 0.95 |
| Interaction, week 10 | +0.09 | −4.85 to 5.03 |

Table: SBP-by-arm interaction on CAPS total, 279 participants

The result did not depend on the analysis choices (computed). The
interaction was −0.51 (−3.06 to 2.04) with screening SBP only, +1.22
(−1.41 to 3.86) with randomization-day SBP only, −0.25 (−3.80 to 3.29)
with supine SBP, and +0.35 (−3.00 to 3.70) in a week-10 ANCOVA.

### 3.2 Translation to the simulation scale

Under the correlation encoding the arm-by-biomarker slope is
$-c_{bm}\,\sigma_{BR}/\sigma_{bm}$ per biomarker unit (`docs/39`). With
$\sigma_{BR} = 8$ and the CSP 563 biomarker SD of 14.2 mm Hg
(computed):

| Quantity | Value |
|---|---|
| CSP 563 interval for $c_{bm}$ | −0.57 to 0.48 |
| Slope implied by $c_{bm} = 0.25$ | −1.41 per 10 mm Hg |
| Slope implied by $c_{bm} = 0.45$ | −2.54 per 10 mm Hg |
| $c_{bm}$ implied by Raskind 2016 (−14 per 10 mm Hg) | about 2.5 |

Table: The interaction on the $c_{bm}$ scale

The project's effect sizes are an order of magnitude below the
published one and are not excluded by CSP 563. The published effect
cannot be represented by the DGP at all, because a correlation cannot
exceed 1 (inferred from the formula). With a smaller drug-response SD
the interval in $c_{bm}$ widens: at $\sigma_{BR} = 5$ it is −0.91 to
0.77.

**Implication.** Manuscripts that cite the 14-points-per-10-mm-Hg
association as motivation should state that it was not replicated and
that the simulated effects are far smaller. The simulated grid stays
defensible as a range of hypothetical biomarker effects. No trial has
yet estimated $c_{bm}$ for SBP with useful precision.

### 3.3 Prognostic association

| Comparison | n | Pearson $r$ (95% CI) |
|---|---|---|
| Standing SBP vs screening CAPS | 303 | 0.02 (−0.09 to 0.13) |
| Supine SBP vs screening CAPS | 303 | 0.01 (−0.10 to 0.12) |
| Placebo arm: standing SBP vs week-10 CAPS change | 132 | −0.07 (−0.24 to 0.10) |

Table: Biomarker and baseline severity (computed)

Correlations larger than about 0.13 in absolute value are excluded.
The DGP's $c_{bm,BL} = 0$ is consistent with these data, and the
`prog-BL` scenarios of paper 02 and `docs/64` remain sensitivity
analyses rather than expected conditions.

## 4. Covariance of the CAPS total

### 4.1 Estimates

Placebo arm, 152 participants, CAPS total at weeks 0, 6, 10, 18 and 26
(computed):

| | 0 | 6 | 10 | 18 | 26 |
|---|---|---|---|---|---|
| 0 | 1 | 0.63 | 0.60 | 0.55 | 0.46 |
| 6 | | 1 | 0.82 | 0.77 | 0.64 |
| 10 | | | 1 | 0.84 | 0.73 |
| 18 | | | | 1 | 0.81 |
| 26 | | | | | 1 |

Table: Empirical correlation of CAPS total, placebo arm

| Model (REML, week-factor mean) | ΔAIC | Parameters |
|---|---|---|
| Unstructured | 0 | as the table above |
| Exponential decay + nugget | 39 | $\nu = 0.10$, $\varphi = 0.975$ per week |
| Random intercept + exponential | 50 | intercept share 0.40, $\varphi = 0.924$ |
| Exponential (AR(1)) | 52 | $\varphi = 0.959$ |
| Compound symmetry | 110 | intraclass correlation 0.68 |

Table: Covariance models for CAPS total, placebo arm (computed; the
nugget model was added after the plan)

The nugget model implies correlations of 0.88, 0.82, 0.74 and 0.55 at
lags of 1, 4, 8 and 20 weeks, against observed 0.82 to 0.84 at 4
weeks, 0.77 at 8 and 0.64 at 20. Screening correlates less with every
later visit than follow-up visits do with each other. Selection on
screening CAPS for eligibility and regression to the mean are
plausible explanations; this was not tested.

### 4.2 Comparison with the extended DGP

The extended DGP's total-score correlation, with $BL$ SD 18.48,
component SDs 10, $10\eta_t$ and 8, $\rho = 0.7$ and $c_1 = 0.2$, is
0.87 at 1 week, 0.67 at 4 weeks and about 0.57 from 8 weeks on
(computed). It decays faster than CAPS at short lags and then
flattens, where CAPS keeps decaying. Matching the observed 4-week and
20-week correlations within the DGP's structure would require $\rho$
of about 0.87 to 0.91 per week (computed).

The variance also differs. The DGP's SD of change from baseline is
16.0 at every visit. In the CSP 563 placebo arm it is 17.0, 19.4, 22.1
and 24.4 at weeks 6, 10, 18 and 26 (computed).

### 4.3 Consequence for design comparisons

The time part of the on-minus-off contrast variance, $a^\top A a$ with
$A = (1 - \nu)\varphi^{\lvert w_t - w_s\rvert} + \nu I$ (computed):

| Structure | Crossover | Hybrid paths | Alternating | Hybrid / crossover |
|---|---|---|---|---|
| Extended DGP, $\rho = 0.7$, no noise | 0.770 | 0.451 to 0.531 | 0.235 | 0.59 to 0.69 |
| `docs/39`, $\rho = 0.8$, 20% noise | 0.763 | 0.389 to 0.466 | 0.232 | 0.51 to 0.61 |
| CSP 563, $\varphi = 0.975$, 10% noise | 0.311 | 0.115 to 0.128 | 0.074 | 0.37 to 0.41 |

Table: Contrast noise by design and covariance structure

Under the CSP 563 structure every within-patient contrast is less
noisy, and interleaving on and off visits helps more. The design
rankings of `docs/39` and `docs/68` Section 9 under AR(1) are, if
anything, understated (inferred). Two limits apply. CSP 563 assessed
CAPS every 4 to 8 weeks, so correlations at the Hybrid design's
one-week lags are extrapolated from the fitted form. And the
calculation covers the noise of the contrast only, not carryover or
the proportional-moderation question of
`analysis/scripts/quick-sim/carryover-closed-form/22-proportional-moderation-designs.R`.

## 5. Response components

### 5.1 Mean trajectories

Adjusted mean CAPS change, CR2 95% CIs, against the DGP's means under
blinding (computed):

| Week | Placebo ($TV + PB$) | DGP | Prazosin minus placebo ($BR$) | DGP |
|---|---|---|---|---|
| 6 | −9.7 (−13.1 to −6.3) | −5.3 | −1.1 (−5.0 to 2.8) | −7.3 |
| 10 | −13.0 (−16.8 to −9.2) | −8.4 | +1.1 (−3.3 to 5.5) | −10.2 |
| 18 | −17.7 (−21.9 to −13.6) | −9.7 | +4.4 (−0.6 to 9.3) | −11.0 |
| 26 | −17.5 (−22.1 to −12.9) | −9.8 | +1.5 (−4.1 to 7.2) | −11.0 |

Table: Mean CAPS change (negative is improvement)

A modified Gompertz curve fitted to the placebo means, with the
displacement fixed at 5, gives a maximum of 17.2 points (SE 0.9) and a
rate of 0.34 per week (SE 0.03). The DGP's $TV + 0.5\,PB$ has maximum
9.8 and rate 0.35 (computed). The shape matches and the size is about
1.8 times larger. How the larger improvement divides between $TV$ and
$PB$ is not identifiable here.

### 5.2 Drug-effect heterogeneity

Residual SD of CAPS change by arm (week-by-week ANCOVA, 2,000
bootstrap resamples within arm; computed):

| Week | SD placebo | SD prazosin | Ratio (95% CI) | $\hat\tau$ upper 95% limit |
|---|---|---|---|---|
| 6 | 16.1 | 16.1 | 1.00 (0.81 to 1.24) | 10.2 |
| 10 | 19.0 | 17.5 | 0.92 (0.75 to 1.12) | 8.3 |
| 18 | 20.7 | 18.2 | 0.88 (0.71 to 1.08) | 7.0 |
| 26 | 23.3 | 21.7 | 0.93 (0.76 to 1.14) | 11.0 |

Table: Outcome SD by arm and the implied bound on $\tau$

A longitudinal test of arm-specific variances gave $\chi^2 = 3.87$, 4 df,
p = 0.42. The bound assumes individual drug effects are uncorrelated
with the other components.

### 5.3 Blood-pressure fall as a marker of drug action

Prazosin lowered supine SBP by 5.7 mm Hg relative to placebo at week
10 (−8.8 to −2.6). Within the placebo arm a larger fall in SBP went
with greater CAPS improvement (3.1 points per 10 mm Hg, 0.4 to 5.8).
The relation was not stronger in the prazosin arm (interaction −2.0,
−5.7 to 1.7; computed). A fall in blood pressure tracks a shared state,
not individual drug response. SBP change is post-randomization; the
analysis is descriptive.

## 6. Dropout

The share of participants seen at one visit but not the next was 10
to 11 percent from screening to week 6 and 3 to 8 percent in later
intervals, similar between arms. In a logistic model over the
post-baseline intervals, CAPS worsening so far did not predict
missing the next visit: odds ratio 0.96 per 10 points of worsening
(0.82 to 1.13; computed). The DGP's `censordata()` makes dropout
depend on worsening (corrected 2026-10-29). CSP 563 shows no such
dependence. The count includes intermittent missed visits as well as
permanent dropout.

## 7. Parameter table

| DGP parameter | Current value | CSP 563 | Status |
|---|---|---|---|
| $BL$ mean, SD | 83.1, 18.5 | 81.1, 16.0 | estimated |
| Biomarker mean, SD | 124.3, 15.4 | 129.8, 14.2 | estimated |
| $c_{bm,BL}$ | 0 | 0.02 (−0.09 to 0.13) | estimated; supports current value |
| $c_{bm}$ | grid 0.25 to 0.45 | −0.57 to 0.48 | not informative at the grid's scale |
| $BR$ mean | Gompertz maximum 11.0 | about 0 | estimated for this trial |
| $TV + PB$ under blinding | maximum 9.8, rate 0.35 | maximum 17.2, rate 0.34 | estimated in combination |
| $TV$ versus $PB$ split | equal (assumption) | none | not identifiable |
| Expectancy weight 0.5 | assumption | none | not identifiable |
| Component SDs (10, $10\eta$, 8) | fixed | SD of change 17 to 24, growing | total only |
| $\rho$, $c_1$ | 0.7, 0.2 | total score: $\varphi = 0.975$, $\nu = 0.10$ | total only |
| $\tau$ (paper 04) | grid | at most about 8 at week 10 | bounded |
| Dropout dependence on worsening | yes | none detected | estimated |
| Carryover half-life | 0.1 to 0.2 week | none | not identifiable |
| Proportional moderation | not used | none | not identifiable |

Table: DGP parameters against CSP 563

## 8. Implications for the manuscripts

- **Papers 01 and 02 (DGP and carryover).** The cited 14 points per
  10 mm Hg did not replicate; say so, and note that the simulated
  $c_{bm}$ grid is an order of magnitude smaller and is not excluded.
  $c_{bm,BL} = 0$ now has empirical support.
- **Paper 04 (main effect).** The drug-response mean of 11 points
  describes the earlier positive trial; CSP 563 found none. Check the
  heterogeneity grid against the bound of about 8 points.
- **Paper 05 and `docs/68` (design).** The CSP 563 covariance
  favors interleaving more than the extended DGP does. A sensitivity
  arm at $\varphi \approx 0.97$ per week with 10 percent noise would
  test whether the design rankings hold.
- **Paper 09 (dropout).** CSP 563 gives no support for
  worsening-dependent dropout. A missing-completely-at-random arm is
  the empirically closer condition for veterans in a blinded trial.
- **Paper 10 (calibration).** Paper 10 traced $\kappa$ above 1 mainly
  to mean-model misfit. The growing variance of change seen here is a
  second, untested candidate (inferred).

## 9. A candidate CSP 563-calibrated sensitivity arm

If the project adds one sensitivity DGP, these values follow from the
estimates (inferred; none has been simulated):

1. $BL$ SD 16; biomarker SD 14.
2. $TV$ and $PB$ maxima each scaled by 1.8, keeping the equal split
   (an assumption this trial cannot test) and the rate 0.35.
3. Component correlation decaying at about 0.97 per week, plus a
   measurement-error component of about 10 percent of total variance.
   The per-component split remains an assumption.
4. Component variances growing over time so that the SD of change
   rises from about 17 to 24 over 26 weeks.
5. Dropout independent of worsening, about 3 to 11 percent per
   interval, highest in the first six weeks.
6. $c_{bm}$ grid unchanged; $c_{bm,BL} = 0$.
7. $BR$ mean at 0 for a CSP 563-like scenario and at 11 for the
   earlier trial, reported side by side.

## 10. Limitations

- One trial, negative overall, in chronic veterans; the DGP's
  calibration source enrolled active-duty soldiers. The two may differ
  in placebo response and drug sensitivity.
- CAPS total was a secondary outcome of CSP 563.
- Visits 4 to 8 weeks apart limit what can be said about one-week
  lags.
- The heterogeneity bound and the component means assume additivity
  and independence of the drug response from the other components.
- The plans were committed before fitting but not registered; the
  nugget models, the parameter summaries and the dropout analysis were
  added afterwards.

## References

- Hendrickson RC, Thomas RG, Schork NJ, Raskind MA (2020). Optimizing
  aggregated N-of-1 trial designs for predictive biomarker validation:
  statistical methods and theoretical findings. *Frontiers in Digital
  Health* 2:13. doi:10.3389/fdgth.2020.00013.
- Raskind MA, Millard SP, Petrie EC, et al. (2016). Higher pretreatment
  blood pressure is associated with greater posttraumatic stress
  disorder symptom reduction in soldiers treated with prazosin.
  *Biological Psychiatry* 80(10):736-742.
  doi:10.1016/j.biopsych.2016.03.2108.
- Raskind MA, Peskind ER, Chow B, et al. (2018). Trial of prazosin for
  post-traumatic stress disorder in military veterans. *New England
  Journal of Medicine* 378(6):507-517. doi:10.1056/NEJMoa1507598.
- Project documents: `docs/39` (design optimization), `docs/40`
  (Hybrid design primer), `docs/64` and paper 02 (prognostic
  scenarios), `docs/68` (design tutorial), papers 04, 09 and 10.
- Analysis sources: `~/prj/alz/14-murray-nejm/murraynejm/analysis/report/04-sbp-caps-replication/`
  (`plan.md`, `results.md`, `plan-components.md`,
  `results-components.md`, scripts `01` to `04`).
