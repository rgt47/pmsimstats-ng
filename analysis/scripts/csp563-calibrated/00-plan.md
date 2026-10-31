# G9 Under a CSP 563-Calibrated DGP: Working Protocol
*2026-10-30 17:32 PDT*

pmsimstats team. Exploratory working protocol, written before any
code or run. It follows `docs/69`, Section 9, and the closed-form
result of
`analysis/scripts/quick-sim/carryover-closed-form/23-csp563-covariance-designs.R`
(commit f1a13be93).

## Notation

| Symbol | Meaning |
|---|---|
| $\rho$ | within-component AR(1) correlation per week (`c.tv`, `c.pb`, `c.br`) |
| $c_1$ | cross-component correlation, $c_1 \rho^{\lvert w_t - w_s \rvert}$ |
| $\sigma_e$ | SD of the added measurement error (new) |
| $g(w)$ | multiplier on the $TV$ and $PB$ SDs at week $w$ (new) |
| $\varphi, \nu$ | CSP 563 exponential decay per week and nugget share |
| $c_{bm}$ | biomarker-$BR$ correlation |
| $t_{1/2}$ | data-generating carryover half-life |
| $\kappa$ | SD of the estimate over the mean reported SE |

## 1. Question

Does G9 keep its size, and does the Hybrid's power advantage over the
one-cycle crossover hold, when the DGP is calibrated to CSP 563?
Script 23 answered the second part for E9 in closed form: the
Hybrid/CO power ratio at $t_{1/2} = 0$ rises from 1.29 under
$\rho = 0.7$ to 1.61 under the CSP covariance. That calculation cannot
address three things this run can:

1. G9's size when its residual model (random intercept, `corCAR1`, no
   nugget, constant variance) is misspecified for the CSP structure.
2. The effect of the means on G9, which models the mean through `t`
   and `Dbc` rather than absorbing it. Paper 10 traced $\kappa$ above
   1 mainly to mean-model misfit.
3. Dropout.

## 2. A calibration point found while planning

`docs/69` compared CSP 563 with the parameters of the closed-form
scripts (`hendrickson-problems/04-power-simulation.R`: $BL$ SD 18.5;
$TV$, $PB$, $BR$ SDs 10, $10\eta$, 8; $TV$ maximum 6.5). The G9 driver
of job 40 (`19-g9-planning-and-fig4b.R`) uses `default_resp_param()`
and `default_baseline_param()` of `simulation-core.R` instead: $BL$
SD 10; $TV$, $PB$, $BR$ SDs 6, $2\eta$, 8; $TV$ maximum 11.0 at rate
0.42; $PB$ maximum 6.5 at rate 0.35. The factor of 1.8 on the placebo
improvement and the comparison of SDs in `docs/69` therefore do not
transfer to this driver. Stage 1 recalibrates against the driver's
own parameters. Which parameter set is Hendrickson's is not settled
here; the two files disagree on the $TV$ curve and should be checked
against 58b32a9 separately.

## 3. DGP arms

All arms keep the refined process otherwise: graded coupling,
anchored carryover, $c_1 = 0.2$, stop on non-PD, fixed allocation to
paths, $c_{bm,BL} = 0$.

| Arm | Covariance | Variance growth | Means | Dropout |
|---|---|---|---|---|
| R | $\rho = 0.7$ (job 40) | none | driver default | none |
| C1 | calibrated $\rho$, $\sigma_e$ | none | driver default | none |
| C2 | as C1 | $g(w)$ | driver default | none |
| C3 | as C1 | $g(w)$ | $TV$, $PB$ scaled; $BR$ maximum 11.0 | none |
| C4 | as C1 | $g(w)$ | as C3; $BR$ maximum 0 | none |
| C3d | as C3 | | | MCAR, 20% of patients |
| C4d | as C4 | | | MCAR, 20% of patients |

Arm R is not rerun: job 40's stored block P supplies it, and one cell
is rerun as a cross-check (Stage 3). The steps C1, C2 and C3 add one
feature at a time, so each difference is attributable.

Grid per arm: designs Hybrid, CO, OL+BDC; $t_{1/2} \in \{0, 0.2,
1\}$ week; $c_{bm} \in \{0, 0.25, 0.35\}$; N = 70. Replicates: 2,000
at $c_{bm} = 0$ (MCSE of size 0.005), 1,000 otherwise. That is 27
cells per arm and 162 in all. The $c_{bm}$ values may change after
the pre-flight check (Stage 2): a larger $\rho$ changes the
feasibility region, and 0.45 is left out for that reason.

The alternating designs of script 23 are not included. Their
advantage there rests on correlations at 2.5-week lags, which CSP 563
did not observe. They can be added as presets in a later run.

## 4. Analyses

- **Primary:** G9, `fit_spec_s7(dat_long, 'E7cr2')` unchanged (as job
  40).
- **Secondary:** G9 with `corExp(form = ~ t | ptID, nugget = TRUE)`
  in place of `corCAR1`, both in the AIC selection and the REML refit.
  It asks whether modeling the nugget changes size or power. It is
  dropped if `--time` projects more than 12 hours of host time with
  it included.
- **Reference:** E9 closed form under each arm's covariance (Stage 1),
  for comparison with script 23.

## 5. Stages

### Stage 0: generator changes

In `implementations/tidyverse/R/functions.R`, two gated additions that
change nothing when absent (the same convention as `sd_slope` and
`b.bm.tv.trend`):

- `model_param$sd_error`: after the MVN draw and outcome construction,
  add independent $\mathcal{N}(0, \sigma_e^2)$ error to every post-baseline
  outcome column. The $BL$ column is not perturbed: in the DGP it is
  the patient's latent level, and the analysis does not use it as an
  observation.
- `model_param$sd_growth`: a per-timepoint multiplier vector applied
  to the $TV$ and $PB$ SDs in `build_sigma_matrix()`. $BR$ is not
  scaled, so the interaction's size in outcome units is unchanged.
  The CSP 563 growth was observed in the placebo arm, which has no
  $BR$ component.

`simulation-core.R`'s `generate_data_multi_path()` passes both
through. `R/generateData.R` is not changed; the parity scripts are
unaffected because both fields are absent in every existing caller.

Tests in `inst/tinytest/` (or the tidyverse collection's tests):

- with both fields absent, `generate_data()` output is identical to
  the current output under a fixed seed;
- with `sd_error` set, the variance of the outcome minus its
  error-free value is within Monte Carlo error of $\sigma_e^2$;
- with `sd_growth` set, the $TV$ and $PB$ diagonal of `sigma` equals
  the unscaled diagonal times $g^2$, and the $BR$ block is unchanged.

### Stage 1: calibration (closed form, no simulation)

Script `01-calibrate.R`. Choose $\rho$, $\sigma_e$ and $g(w)$ for the
driver's parameters so that, in a blinded placebo-only pseudo-design
with visits at weeks 6, 10, 18 and 26:

- the correlations of the post-baseline outcome match the CSP 563
  exponential-plus-nugget fit (0.82, 0.74 and 0.55 at lags of 4, 8
  and 20 weeks) within 0.05;
- under C2 and later, the SD of change from baseline matches CSP 563
  (17.0, 19.4, 22.1, 24.4) within 1 point at each week.

The calibration fits the outcome's total-score structure. The driver
adds a random intercept ($BL$, SD 10) that CSP 563's fitted model did
not have, so the correlation floor differs; if no $\rho$, $\sigma_e$
pair meets the first criterion with $BL$ SD 10, the $BL$ SD becomes a
third calibration parameter and the change is reported.

The $TV$ and $PB$ scaling for C3 is the factor that makes the
driver's blinded placebo-arm mean improvement, $TV + 0.5\,PB$, reach
the CSP 563 Gompertz maximum of 17.2 points. The equal treatment of
$TV$ and $PB$ is kept, as `docs/69` notes, as an assumption.

Outputs: the calibrated values, the implied correlation and SD tables
against CSP 563, and the E9 closed-form power for every cell of every
arm. The calibrated values are fixed in a committed CSV before Stage
3.

### Stage 2: pre-flight

Every (arm, design, $c_{bm}$, $t_{1/2}$) covariance matrix is checked
for positive definiteness, as in script 19. Infeasible $c_{bm}$
values are replaced by the largest feasible value on a 0.05 grid,
recorded before any replicate is run.

### Stage 3: smoke test, timing and cross-check

Driver `02-run-g9.R`, modeled on script 19 (per-replicate seeds from
`label_seed()`, per-replicate checkpoints, `--smoke`, `--time`,
`--preflight`, `--cross-check`). The cross-check reruns job 40's
Hybrid, $c_{bm} = 0.25$, $t_{1/2} = 0$, N = 70 cell under arm R with
200 replicates and compares it with the stored value; agreement within
2 MCSE is required before the full run. Results from the container
and a local session are not draw-for-draw identical, so the
cross-check runs on the host that runs the job.

### Stage 4: run

`simrun` on zeus (8 cores per job, worker count from
`parallelly::availableCores()`), or on `sim-mercury-ts` if zeus is
busy. BLAS threads set to 1.

### Stage 5: summaries and write-up

Per cell: size or power with MCSE, mean estimate, $\kappa$, the
distribution of the selected half-life, non-convergence. Comparisons:

- each arm against R, and adjacent arms against each other (C1-R,
  C2-C1, C3-C2, C4-C3, C3d-C3);
- the Hybrid/CO power ratio by arm, against script 23's E9 ratio;
- G9 against its nugget variant.

A results note in this directory, then updates to `docs/69` Section
9 and, if the design conclusions change, `docs/68` Section 9.2.

## 6. Expectations stated before the run

These are stated so the results can be read against them; they are
not tests with a decision rule.

1. The Hybrid/CO G9 power ratio is larger under C1 than under R, as
   for E9 in script 23.
2. G9 size stays within 0.05 plus or minus 2 MCSE (0.040 to 0.060)
   under C1 and C2. CR2 protects the test against a misspecified
   correlation model; this is the expectation most at risk.
3. C3 changes G9's size or $\kappa$ relative to C2, through mean-model
   misfit (paper 10). The direction is not predicted.
4. C4 and C3 give similar power. Under the correlation encoding, the
   interaction acts through the variability of $BR$, not its mean;
   but $BR$ mean 0 removes the carryover signal that guides the AIC
   choice of half-life, so G9 may select differently.
5. Dropout of 20 percent lowers power without changing size.

## 7. Limitations known in advance

- The CSP 563 structure is estimated for the total score; applying a
  single $\rho$ to every component and adding the nugget as
  measurement error are assumptions.
- CSP 563 visits were 4 to 8 weeks apart; the Hybrid's one-week lags
  are extrapolated.
- `apply_dropout(per = 'patient')` never removes the first occasion,
  whereas CSP 563's highest per-interval rate was before week 6. The
  20 percent figure approximates CSP 563's cumulative missingness by
  week 26 and includes intermittent misses, which this mechanism
  treats as permanent.
- One trial, negative overall, in chronic veterans.

## 8. Cost

Script 19's `--time` estimate will be repeated for this grid. A rough
expectation before timing: 162 cells, about 216,000 replicates at
about 1.5 core-seconds each for G9, or 90 core-hours; about 6 hours of
host time on 8 cores at the host's per-core speed, and roughly half
again with the nugget variant. These figures are guesses until Stage 3
measures them.
