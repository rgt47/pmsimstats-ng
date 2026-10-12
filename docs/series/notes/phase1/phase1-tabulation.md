# Phase 1 tabulation: anchored reruns of the stress test and covariance study
*2026-10-11 22:41 PDT*

Author: pmsimstats team

Repository: `/Users/zenn/prj/res/36-pmsimstats-ng/pmsimstats-ng` (REPO).
Working files: `phase1/` in this scratchpad (scripts, logs, CSVs).
Nothing was committed.

## 0. Evidence status and procedure

- **Verified (run and observed).** `12-strawman-stress-decision.R` was
  run from REPO with single-threaded BLAS (`OPENBLAS_NUM_THREADS`,
  `OMP_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `MKL_NUM_THREADS` all 1)
  as `--tag -anchored` and `--strawman DS --no-symmetry --tag
  -anchored`, without `--with-extras` (`phase1/run12.sh`; logs
  `phase1/dec-s-anchored.log`, `phase1/dec-ds-anchored.log`). Outputs in
  `strawman-stress/`: `decision-{screen,sensitivity,power,maximin}-anchored.csv`
  and `decision-{...}-ds-nosym-anchored.csv` (8 new, untracked files).
- **Stored-run baseline recomputed.** The stored untagged
  `decision-*-ds-nosym.csv` files (2026-10-07 10:20) predate the
  addition of the baseline-prognostic cells 40 and 41 to
  `replicates-full.rds` (11:50): the screen file has no `prog-BL`
  column. For a like-for-like comparison the 11-analysis decision was
  recomputed from the current `replicates-full.rds` (41 cells) into
  `phase1/stored-recomputed/` (`--dir`, so no stored file was touched;
  `phase1/run12-stored.sh`). Check: the recomputed `-all14` files are
  identical to the four stored `decision-*-ds-nosym-all14.csv` files
  (verified, `all.equal` TRUE for each). The stored untagged
  `-ds-nosym` and plain `decision-*.csv` files were left stale as they
  are.
- **Comparison script.** `phase1/compare.R` (output
  `phase1/compare.out`; per-entry tables `phase1/power-old-new.csv`,
  `phase1/summary-old-new.csv`).
- **Covariance study.** `09-covariance-study-tables.R --tag -anchored
  --reps 250,500` failed as written, because it requires the follow-up
  and Dbc passes, which were cancelled for the anchored run (DQ-8). A
  minimal change to 09 (uncommitted) skips any pass with no output for
  the tag and prints which; defaults and the untagged behavior are
  unchanged (verified: rerunning 09 on the stored recursive passes in a
  scratch `--dir` reproduces `summary-alt250-null500-combined.csv` and
  `size-adjusted-power-alt250-null500.csv` byte for byte, and the
  stored files' checksums are unchanged). Outputs:
  `covariance-study/summary-alt250-null500-combined-anchored.csv` and
  `size-adjusted-power-alt250-null500-anchored.csv` (main pass only).
  Comparison: `phase1/compare-cov.R`, `phase1/compare-cov.out`.
- **Convention.** The decision script counts a failed fit as a
  non-rejection (stored decision files likewise). `summary-full*.csv`
  excludes failed fits. The two differ by at most 0.002 in any rate.

## 1. Stress test: what changed

| Feature | Stored run (2026-10-07) | Rerun (job 15, `e077526`) |
|---|---|---|
| Carryover on the BR mean | published recursion (elapsed time counted twice from the 2nd off-drug visit) | anchored, $\mu_{\text{last}} 2^{-t_{sd}/t_{1/2}}$ (Weibull weight in the Weibull cells) |
| Path allocation | fixed published allocation | uniform random per participant, redrawn if a path has fewer than 2 |
| Cells, trials, seeds | 41, 28,000, per-chunk seeds | same |
| Analyses | 11 main + G1, G2 (extra) + F1A (extra2) | 11 main only; G1, G2, F1A pending (jobs 21, 22) |
| Pairing with the other run | | none: means and random streams differ; independent estimates |

Failed fits: 26 of 308,000 (stored 28). By analysis, rerun (stored): S
3 (6), DS 3 (2), C1 4 (6), C2 4 (6), F1 4 (3), F2 3 (2), ES 5 (3), EA,
A, DA, R 0 (0).

## 2. Size (11 non-prognostic null cells, 11,000 trials; band 0.0458 to 0.0542; cell cap 0.075)

| Analysis | Pooled old | Pooled new | $z$ | Largest old | Largest new | Nominal old / new |
|---|---|---|---|---|---|---|
| S | 0.0483 | 0.0516 | 1.14 | 0.059 MCAR | 0.061 N35 | yes / yes |
| DS | 0.0528 | 0.0569 | 1.39 | 0.068 MCAR | 0.067 Weibull 0.5 | yes / **no** (pooled) |
| C1 | 0.0531 | 0.0529 | -0.06 | 0.069 N35 | 0.065 meas. error | yes / yes |
| C2 | 0.0519 | 0.0492 | -0.93 | 0.066 N35 | 0.055 meas. error | yes / yes |
| F1 | 0.0511 | 0.0514 | 0.09 | 0.062 AR(1) curve | 0.064 thalf0 | yes / yes |
| F2 | 0.0486 | 0.0519 | 1.11 | 0.059 MCAR | 0.060 N35 | yes / yes |
| ES | 0.0478 | 0.0490 | 0.40 | 0.059 AR(1) curve | 0.060 Weibull 0.5 | yes / yes |
| EA | 0.0517 | 0.0575 | 1.98 | 0.070 N35 | 0.076 N35 | yes / **no** (pooled and cell) |
| A | 0.0573 | 0.0571 | -0.06 | 0.087 N35 | 0.084 N35 | no / no |
| DA | 0.0653 | 0.0656 | 0.12 | 0.103 N35 | 0.093 N35 | no / no |
| R | 0.0687 | 0.0686 | -0.03 | 0.081 N35 | 0.084 N35 | no / no |

$z$ uses $\sqrt{2 \times 0.05 \times 0.95 / 11000}$.

Prognostic null cells (1,000 trials each), old / new:

| Analysis | TV const | TV grow | Placebo | Baseline | Flag old / new |
|---|---|---|---|---|---|
| S | 0.099 / 0.099 | 0.049 / 0.055 | 0.100 / 0.114 | 0.050 / 0.056 | yes / yes |
| DS | 0.115 / 0.127 | 0.058 / 0.068 | 0.112 / 0.130 | 0.055 / 0.057 | yes / yes |
| C1 | 0.382 / 0.367 | 0.173 / 0.196 | 0.205 / 0.191 | 0.054 / 0.056 | yes / yes |
| C2 | 0.393 / 0.382 | 0.172 / 0.190 | 0.193 / 0.178 | 0.050 / 0.051 | yes / yes |
| F1 | 0.047 / 0.050 | 0.046 / 0.050 | 0.040 / 0.048 | 0.040 / 0.052 | no / no |
| F2 | 0.090 / 0.094 | 0.048 / 0.052 | 0.098 / 0.112 | 0.048 / 0.054 | yes / yes |
| ES | 0.049 / 0.045 | 0.063 / 0.073 | 0.049 / 0.063 | 0.050 / 0.057 | no / no |
| EA | 0.049 / 0.049 | 0.057 / 0.057 | 0.050 / 0.063 | 0.066 / 0.061 | no / no |
| A | 0.416 / 0.411 | 0.195 / 0.210 | 0.222 / 0.193 | 0.051 / 0.063 | yes / yes |
| DA | 0.509 / 0.515 | 0.245 / 0.276 | 0.256 / 0.229 | 0.064 / 0.065 | yes / yes |
| R | 0.107 / 0.119 | 0.053 / 0.077 | 0.159 / 0.167 | 0.075 / 0.069 | yes / yes |

Cell-by-cell null rates (summary files): 7 of 165 (15 cells by 11
analyses) have $|z| > 2$ (about 7.5 expected): EA thalf0 0.040 to
0.069 (2.86), DA MAR 0.075 to 0.049 (-2.41), S thalf0 0.036 to 0.058
(2.32), EA Weibull 2 0.042 to 0.064 (2.20), R TV grow 0.053 to 0.077
(2.18), EA meas. error 0.043 to 0.065 (2.18), F2 thalf0 0.037 to
0.057 (2.12). Full $z$ grid in `phase1/compare.out`, section 1b.

## 3. Sensitivity and survivors

All curves monotone in both constructs, both runs. Symmetry (used
only in the original reading) flips between runs: AR(1) curve S, DS,
F2 asymmetric in the stored run, symmetric in the rerun; CS curve S,
DS, C2, F2 symmetric in the stored run, asymmetric in the rerun; C1,
ES, EA, A the reverse on the CS curve.

| Reading | Stored (11 analyses, 41 cells) | Rerun |
|---|---|---|
| Amendment 1 (DS strawman, no symmetry): survivors | C1, C2, S, F2, DS, F1, ES, EA | C2, C1, F2, S, F1, ES |
| of which unflagged | F1, ES, EA | F1, ES |
| strawman | DS survives (flagged) | DS fails size: strawman rule not reached |
| Original (S strawman, symmetry on): survivors | C2, F1 | C1, F1, ES |
| strawman | S fails symmetry | S fails symmetry (CS curve) |

Detectable coupling (power 0.80 on the AR(1) curve), old / new: S
0.446 / 0.433, DS NA / 0.433, C1 0.403 / 0.387, C2 0.403 / 0.394, F2
NA / 0.439, A 0.430 / 0.406, DA 0.434 / 0.412; F1, ES, EA, R not
reached in either. Steepness (0.15 to 0.30), old / new: F1 1.80 /
2.09, ES 1.71 / 1.63, EA 1.48 / 1.64, S 2.32 / 2.35.

## 4. Maximin shortfall (survivors; best per cell over all survivors)

| Analysis | Non-prog. old | Non-prog. new | worst new | All old | All new | worst new |
|---|---|---|---|---|---|---|
| C2 | 0.040 | 0.022 | curve-C07 0.45 | 0.040 | 0.058 | prog-TV-const |
| C1 | 0.034 | 0.032 | meas-error | 0.064 | 0.032 | meas-error |
| F2 | 0.160 | 0.164 | weibull2 | 0.244 | 0.208 | prog-TV-grow |
| S | 0.144 | 0.172 | weibull2 | 0.220 | 0.172 | weibull2 |
| DS | 0.182 | (fails size) | | 0.214 | | |
| F1 | 0.242 | 0.240 | weibull0.5 | 0.358 | 0.348 | prog-PB |
| ES | 0.278 | 0.266 | weibull0.5 | 0.544 | 0.498 | prog-TV-grow |
| EA | 0.318 | (fails size) | | 0.528 | | |

Stored `-all14` additionally had G2 0.232, G1 0.252 (non-prognostic)
and G1 0.308, G2 0.312 (all cells); the 11-analysis stored rows are
identical to the `-all14` rows (G1, G2 were never the best in a cell).
Original reading (rerun): C1 0.022, F1 0.234, ES 0.260.

Among the 11 analyses F1 leads ES by 0.026 (outside the 0.02 tie
margin). The stored-run selection of F1 over G2 (0.010 apart, decided
by the last tie-break) cannot be rechecked until G1 and G2 are rerun.

## 5. Size-adjusted power (unadjusted in parentheses)

| Analysis | AR(1) 0.30 old | AR(1) 0.30 new | $z$ | CS 0.25 old | CS 0.25 new | $z$ |
|---|---|---|---|---|---|---|
| S | 0.514 (0.508) | 0.504 (0.506) | -0.32 | 0.554 (0.554) | 0.490 (0.498) | -2.03 |
| DS | 0.486 (0.516) | 0.502 (0.522) | 0.51 | 0.480 (0.476) | 0.372 (0.404) | -3.47 |
| C1 | 0.578 (0.602) | 0.612 (0.594) | 1.10 | 0.660 (0.658) | 0.624 (0.632) | -1.19 |
| C2 | 0.606 (0.606) | 0.606 (0.602) | 0.00 | 0.656 (0.630) | 0.616 (0.604) | -1.32 |
| F1 | 0.364 (0.410) | 0.448 (0.424) | 2.71 | 0.522 (0.544) | 0.510 (0.534) | -0.38 |
| F2 | 0.516 (0.496) | 0.482 (0.490) | -1.08 | 0.552 (0.552) | 0.488 (0.494) | -2.03 |
| ES | 0.364 (0.412) | 0.376 (0.392) | 0.39 | 0.662 (0.634) | 0.602 (0.572) | -1.97 |
| EA | 0.306 (0.392) | 0.390 (0.392) | 2.80 | 0.602 (0.596) | 0.564 (0.548) | -1.22 |
| A | 0.534 (0.566) | 0.556 (0.562) | 0.70 | 0.618 (0.618) | 0.566 (0.556) | -1.68 |
| DA | 0.498 (0.592) | 0.540 (0.568) | 1.33 | 0.572 (0.596) | 0.524 (0.542) | -1.53 |
| R | 0.398 (0.476) | 0.448 (0.482) | 1.60 | 0.596 (0.594) | 0.488 (0.546) | -3.45 |

$z$ uses the binomial MCSE of each size-adjusted rate (it ignores the
noise of the estimated critical value, so it slightly overstates
significance).

**F1 against S.**

| Quantity | Stored | Rerun |
|---|---|---|
| AR(1), $c_{bm} = 0.30$: S, F1 | 0.514, 0.364 | 0.504, 0.448 |
| F1 deficit vs S at that cell | 29.2% | 11.1% |
| CS, $c_{bm} = 0.25$: S, F1 | 0.554, 0.522 | 0.490, 0.510 |
| F1 deficit at CS cell | 5.8% | -4.1% (F1 ahead) |
| Mean F1/S ratio, 20 non-prognostic alt cells with $c_{bm} > 0$ | 0.831 | 0.877 |
| McNemar F1 vs S, 22 non-prog. cells: F1 lower / higher | 13 / 0 | 12 / 0 |
| McNemar F1 vs DS: F1 lower / higher | 15 / 0 | 12 / 3 |
| McNemar ES vs S: ES lower / higher | 15 / 3 | 14 / 5 |
| McNemar C1 vs S: C1 higher / lower | 20 / 0 | 20 / 0 |
| McNemar F1 vs ES: F1 lower / higher | 5 / 3 | 5 / 4 |

F1/S ratio by cell, old / new: N150 0.90 / 0.90; N35 0.86 / 0.81;
C07 0.075 0.76 / 1.06; 0.15 0.57 / 0.88; 0.225 0.75 / 0.91; 0.30 0.71 /
0.89; 0.375 0.81 / 0.94; 0.45 0.83 / 0.93; CS 0.05 0.97 / 1.00; 0.10
0.89 / 0.82; 0.15 0.88 / 0.91; 0.20 0.98 / 0.92; 0.25 0.94 / 1.04; MAR
0.91 / 0.75; MCAR 0.89 / 0.80; meas. error 0.80 / 0.71; thalf0 0.84 /
0.80; thalf1 0.70 / 0.73; Weibull 0.5 0.76 / 0.74; Weibull 2 0.86 /
1.01. The AR(1) curve moved toward F1 throughout; the stress cells
moved slightly against it.

**Stress cells, size-adjusted power at 0.30, old / new.**

| Cell | S | C1 | DS | F1 | ES | EA |
|---|---|---|---|---|---|---|
| N150 | 0.828 / 0.840 | 0.918 / 0.906 | 0.822 / 0.832 | 0.748 / 0.758 | 0.704 / 0.746 | 0.698 / 0.728 |
| N35 | 0.222 / 0.258 | 0.272 / 0.360 | 0.232 / 0.264 | 0.192 / 0.208 | 0.206 / 0.188 | 0.196 / 0.198 |
| thalf0 | 0.594 / 0.498 | 0.666 / 0.556 | 0.588 / 0.508 | 0.496 / 0.398 | 0.504 / 0.428 | 0.490 / 0.380 |
| thalf1 | 0.388 / 0.374 | 0.502 / 0.454 | 0.412 / 0.368 | 0.272 / 0.274 | 0.288 / 0.276 | 0.294 / 0.250 |
| Weibull 0.5 | 0.404 / 0.392 | 0.514 / 0.522 | 0.376 / 0.358 | 0.308 / 0.288 | 0.270 / 0.262 | 0.230 / 0.260 |
| Weibull 2 | 0.588 / 0.494 | 0.682 / 0.666 | 0.588 / 0.518 | 0.504 / 0.498 | 0.458 / 0.466 | 0.464 / 0.412 |
| meas. error | 0.400 / 0.396 | 0.444 / 0.462 | 0.402 / 0.398 | 0.318 / 0.282 | 0.336 / 0.286 | 0.324 / 0.274 |
| MCAR | 0.426 / 0.446 | 0.542 / 0.556 | 0.426 / 0.438 | 0.378 / 0.356 | 0.360 / 0.352 | 0.346 / 0.316 |
| MAR | 0.418 / 0.458 | 0.518 / 0.508 | 0.424 / 0.464 | 0.382 / 0.344 | 0.356 / 0.356 | 0.322 / 0.308 |
| prog TV const | 0.606 / 0.598 | 0.774 / 0.752 | 0.590 / 0.610 | 0.482 / 0.456 | 0.370 / 0.402 | 0.378 / 0.396 |
| prog TV grow | 0.516 / 0.544 | 0.672 / 0.684 | 0.522 / 0.542 | 0.378 / 0.430 | 0.192 / 0.206 | 0.208 / 0.258 |
| prog PB | 0.652 / 0.586 | 0.736 / 0.746 | 0.610 / 0.576 | 0.446 / 0.398 | 0.498 / 0.448 | 0.438 / 0.386 |
| prog BL | 0.444 / 0.490 | 0.552 / 0.582 | 0.438 / 0.502 | 0.396 / 0.416 | 0.336 / 0.406 | 0.290 / 0.358 |

**Run-to-run variation.** Unadjusted power: 15 of 286 entries with
$|z| > 2$; mean change by analysis between -0.001 and +0.008. Mean
change by cell (over analyses) from -0.048 (CS 0.25) to +0.054 (prog
BL); SD about 0.028, against about 0.032 expected for the difference
of two independent 500-trial rates near 0.5. Size-adjusted power: 52
of 286 entries with $|z| > 2$, but entries within a cell share trials,
so they move together (thalf0 lower for every analysis; AR(1) 0.45
higher for every analysis). The thalf0 shift (S 0.594 to 0.498) is in
a cell whose means anchoring does not change, so it reflects the
allocation change or chance; not separated (`--allocation fixed` would
separate it).

**AIC selection shares (all cells).**

| Analysis | Run | 0 | 0.25 | 0.5 | 1 | 2 |
|---|---|---|---|---|---|---|
| DS | stored | 0.170 | 0.144 | 0.299 | 0.208 | 0.179 |
| DS | rerun | 0.216 | 0.148 | 0.250 | 0.173 | 0.213 |
| DA | stored | 0.188 | 0.144 | 0.294 | 0.181 | 0.193 |
| DA | rerun | 0.234 | 0.143 | 0.249 | 0.156 | 0.219 |

## 6. Covariance study, main pass (paired: fixed allocation, same seeds)

120,000 fits, 0 failed, both runs. Same seven analyses by three time
forms. $t_{1/2} = 0$ cells: every rejection rate and mean estimate
identical (verified). $t_{1/2} = 1$ cells: null rates move at most
0.004, raw power at most 0.012, size-adjusted power at most 0.028, mean
estimate at most 0.00045. No size class changes in any time form.

Size over the eight null cells, linear time, recursive / anchored:

| Analysis | Mean | Min | Max CS | Max C07 Hybrid | Max C07 CO | Min $\kappa$ | Class |
|---|---|---|---|---|---|---|---|
| `cs-Satt` | 0.079 / 0.079 | 0.042 / 0.042 | 0.058 / 0.058 | 0.066 / 0.066 | 0.158 / 0.158 | 0.73 / 0.73 | liberal |
| `csh-Satt` | 0.077 / 0.077 | 0.040 / 0.040 | 0.060 / 0.060 | 0.054 / 0.052 | 0.164 / 0.162 | 0.74 / 0.74 | liberal |
| `sp_exp-Satt` | 0.042 / 0.042 | 0.028 / 0.028 | 0.062 / 0.062 | 0.048 / 0.048 | 0.034 / 0.034 | 0.95 / 0.95 | conservative |
| `us-Satt` | 0.090 / 0.090 | 0.072 / 0.072 | 0.092 / 0.092 | 0.102 / 0.102 | 0.104 / 0.104 | 0.84 / 0.84 | liberal |
| `us-KR` | 0.058 / 0.058 | 0.052 / 0.052 | 0.062 / 0.062 | 0.066 / 0.068 | 0.056 / 0.056 | 0.94 / 0.94 | slightly liberal |
| `BL+cs-Satt` | 0.077 / 0.077 | 0.038 / 0.038 | 0.052 / 0.052 | 0.070 / 0.072 | 0.160 / 0.160 | 0.73 / 0.73 | liberal |
| `BL+us-Satt` | 0.087 / 0.087 | 0.072 / 0.074 | 0.094 / 0.094 | 0.092 / 0.092 | 0.094 / 0.094 | 0.86 / 0.86 | liberal |

Spline and visit forms: means within 0.001 (us-KR visit 0.058 to
0.057); all classes unchanged (`phase1/compare-cov.out`).

Size-adjusted power, linear time, recursive / anchored (Hybrid):

| Analysis | CS 0 | CS 1 | C07 0 | C07 1 |
|---|---|---|---|---|
| `us-KR` | 0.688 / 0.688 | 0.424 / 0.408 | 0.940 / 0.940 | 0.756 / 0.764 |
| `us-Satt` | 0.680 / 0.680 | 0.424 / 0.408 | 0.940 / 0.940 | 0.756 / 0.764 |
| `cs-Satt` | 0.664 / 0.664 | 0.456 / 0.452 | 0.844 / 0.844 | 0.676 / 0.676 |
| `csh-Satt` | 0.684 / 0.684 | 0.456 / 0.456 | 0.888 / 0.888 | 0.704 / 0.704 |
| `sp_exp-Satt` | 0.496 / 0.496 | 0.312 / 0.312 | 0.920 / 0.920 | 0.628 / 0.628 |
| `BL+us-Satt` | 0.704 / 0.704 | 0.408 / 0.412 | 0.784 / 0.784 | 0.456 / 0.464 |
| `BL+cs-Satt` | 0.752 / 0.752 | 0.448 / 0.452 | 0.744 / 0.744 | 0.420 / 0.428 |

CO: `us-KR` 0.720 / 0.720, 0.748 / 0.752, 0.804 / 0.804, 0.768 /
0.764; other analyses change by at most 0.004.

Rows of docs/45 quoted from the main pass:

| docs/45 location | Recursive | Anchored |
|---|---|---|
| 4.1 `us-KR` max C07 Hybrid | 0.066 | 0.068 |
| 4.4 `cs-Satt` CS Hybrid $t_{1/2}=1$ linear / spline / visit | 0.054 / 0.070 / 0.072 | unchanged |
| 4.4 `us-KR` linear / visit | 0.056 / 0.072 | 0.056 / 0.070 |
| 5 `us-KR` Hybrid CS 1 | 0.424 | 0.408 |
| 5 `us-KR` Hybrid C07 1 | 0.756 | 0.764 |
| 5 `cs-Satt` Hybrid CS 1 | 0.456 | 0.452 |
| 5 `us-KR` CO CS 1 / C07 1 | 0.748 / 0.768 | 0.752 / 0.764 |
| 6 `us-KR` visit raw, C07 Hybrid 1 | 0.796 | 0.804 (all other cells unchanged) |
| 7 $D_{it}$ row (us-KR linear raw), CS Hybrid 1 | 0.436 | 0.440 (all other cells unchanged) |
| 7 estimate CS Hybrid 1, $D_{it}$ | $-0.1000$ | $-0.0999$ |

Rows from the follow-up and Dbc passes (KR/CR2 of cs, csh, sp_exp;
RI+CAR1-CR2; us-CR2; BL+ KR and CR2; us-KR-sepbm; all Dbc rows) remain
recursive (DQ-8).

## 7. Surprises and caveats

1. **The 29% F1 deficit was a single-cell extreme.** At the AR(1)
   reference cell F1 rose 0.364 to 0.448 ($z = 2.7$) while S barely
   moved; the average F1/S ratio moved only 0.83 to 0.88. Inferred:
   mostly Monte Carlo error in one 500-trial cell, plus possibly a
   real effect of the DGP change on the AR(1) curve (the whole curve
   moved toward F1). Not separable without the `--allocation fixed`
   run.
2. **DS and EA now fail size** (pooled 0.0569 and 0.0575; EA 0.076 at
   N35), so the Amendment 1 strawman is out before its rule is reached.
   Both are boundary crossings within about 1.4 to 2 $z$ of the stored
   values.
3. **ES is close to the prognostic cap** in the TV-growing cell (0.073
   against 0.075).
4. **The symmetry screen changes its verdict between runs** (S fails it
   on the CS curve in the rerun; survivors of the original reading
   change from C2, F1 to C1, F1, ES), consistent with docs/48's
   argument that it tests noise.
5. **The covariance study is almost insensitive to anchoring** (paired:
   estimates move by at most 0.0005). This suggests (inferred) that in
   the stress test too, most old-versus-new differences come from the
   allocation change and from independent Monte Carlo error rather
   than from anchoring.
6. **Stale stored files.** `decision-*-ds-nosym.csv` and
   `decision-*.csv` (untagged) in `strawman-stress/` are from 39 cells
   (no `prog-BL`); docs/47 quotes the `-all14` files, which are
   current. Not changed.
7. **09 failed without the cancelled passes;** fixed by a small
   uncommitted change (skip missing passes).
