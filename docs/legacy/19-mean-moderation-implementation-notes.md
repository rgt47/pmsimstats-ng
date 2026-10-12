# Mean Moderation Implementation Notes {.unlisted .unnumbered}
*Written about 2026-04-08; first committed to this repository 2026-05-01 16:23 PDT*

*Archived 2026-10-11; current treatment: Appendix I (legacy record) Section 5, with a pointer from `series/ch07`; see `docs/legacy/INDEX.md`.*

**Author.** pmsimstats team

**Purpose.** This note records how the mean-moderation architecture
(Architecture A) was implemented in `pmsimstats-ng`: three incorrect
formulations, the calibrated formulation that was adopted, and a
small pilot comparison of its power with the covariance-moderation
architecture (Architecture B) under carryover. It is a quick-reference
companion to the `.tex` version. *Review note (2026-10-07):* code
claims were re-checked against the current `R/generateData.R`. The
pilot power figures (50 replicates) have been superseded by paper 01
and by the closed-form results of `docs/37`. Dated notes mark these
points below.

> **Note**: This is a quick-reference summary. For the full
> mathematical derivation, the positive-definiteness ceiling
> discussion, the `validateParameterGrid()` known limitation,
> and the empirical evaluation of the Hendrickson (2020)
> publication parameter grid under Architecture A, see the
> comprehensive document `19-mean-moderation-implementation-notes.tex`
> (rendered as `19-mean-moderation-implementation-notes.pdf`).

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

The notation follows `analysis/report/NOTATION.md`. *Canonical*
symbols have the meaning given there. *Local* symbols are specific to
this note. Code identifiers appear only in code spans. The package
passes the moderation strength as `c.bm` under both architectures. In
the mathematics it is $\beta_{bm}$ under mean moderation and $c_{bm}$
under covariance moderation. The pilot table reports both at the
shared value $c_{bm} = \beta_{bm} = 0.45$ (NOTATION rule 2).

| Symbol | Meaning | Status |
|---|---|---|
| $BR$ | biologic (drug) response component | canonical |
| $B$ | raw biomarker (code `bm`) | canonical |
| $b$ | standardized biomarker $(B - \mu_{bm})/\sigma_{bm}$ (code `bm_z`) | canonical |
| $\mu_{bm}$, $\sigma_{bm}$ | biomarker mean and SD (124.33 and 15.36 in `extracted_bp`) | local ($\mu_{bm}$); canonical ($\sigma_{bm}$) |
| $\mu_{BR}$, $\sigma_{BR}$ | mean and SD of $BR$ ($\sigma_{BR} = 8$ in `extracted_rp`) | local ($\mu_{BR}$); canonical ($\sigma_{BR}$) |
| $\beta_{bm}$ | mean-moderation multiplier | canonical |
| $c_{bm}$ | covariance-moderation correlation | canonical |
| $t_{1/2}$ | carryover half-life, weeks | canonical |
| $t_{sd}$, $t_{od}$ | time since discontinuation, time on drug (code `tsd`, `tod`) | canonical |
| $N$ | total number of participants (convention of the pilot not re-verified) | canonical |

Table: Notation index: symbols used in this note

**Glossary.** *Architecture A* is the former label for the
mean-moderation architecture (`dgp_architecture = 'mean_moderation'`),
and *Architecture B* the former label for the covariance-moderation
architecture (`'mvn'`). *N-of-1* in the table is the Hendrickson
Hybrid design. *Calibration* means choosing the mean-moderation shift
so that a one-SD biomarker difference moves $BR$ by the same amount
as in the covariance-moderation conditional expectation.

## Summary

This document records the implementation of Architecture A
(direct mean moderation) for the biomarker-treatment
interaction in `pmsimstats-ng`, including two failed
approaches and the final correct formulation derived from
the original 2025 codebase.

## Background

The original R package code (`R/generateData.R`) uses
Architecture B: the BM-BR interaction is encoded as a
differential correlation in the covariance matrix of the
joint MVN draw. Under carryover, this correlation erodes as
the drug effect decays, producing the 40-60% power drops
documented in `05-revision-summary.pdf`. (That file is no longer in
`docs/`; flagged 2026-10-07. The 40-60% figure is a 50-replicate
pilot result. See the note after the power table below.)

The 2025 tidyverse codebase (pre-alignment version in
`~/prj/alz/01-pmsimstats/pmsimstats2025`, initial commit
957c8bd) used Architecture A: additive mean moderation
applied post-MVN draw. That version showed full power
preservation under carryover.

## Implementation Attempts

### Attempt 1: Multiplicative scaling (incorrect)

```r
dat[, (br_col) := get(br_col) * (1 + beta_bm * bm_centered)]
```

**Problem:** When the base BR value is near zero (off-drug
timepoints with minimal carryover), multiplying by
`(1 + small)` produces negligible moderation signal. The
analysis model sees strong biomarker signal on-drug but
almost none off-drug, mimicking Architecture B's correlation
erosion. Power dropped about 60% under carryover, no
improvement over MVN.

### Attempt 2: Additive with drug exposure decay (incorrect)

```r
drug_exposure <- (1/2)^(tsd / t1half)  # on-drug=1, decays off-drug
dat[, (br_col) := get(br_col) + beta_bm * bm_centered *
  drug_exposure * br_sd]
```

**Problem:** Drug exposure decays rapidly during off-drug
periods (e.g., 0.18 at first off-drug timepoint, 0.001 by
the fourth, at 2.5-week spacing with $t_{1/2} = 1$; computed). The
additive signal is proportional to this decay, so the off-drug
moderation is still negligible relative to noise
($\sigma_{BR} = 8$; corrected 2026-10-07 from `br_sd = 5`, against
`extracted_rp`). The `bm:Dbc` interaction
term sees minimal contrast. Power remained low under
carryover.

### Attempt 3: Additive with time-on-drug (incorrect scaling)

```r
dat[, (br_col) := get(br_col) + bm * tod * beta_bm]
```

**Source:** Original 2025 codebase, initial commit 957c8bd,
`analysis/scripts/pm_functions.R` lines ~650-670:

```r
treatment_status <- trial_design$tod[timepoint_idx]
participant_data[[br_col]] <- participant_data[[br_col]] +
  biomarker * treatment_status * interaction_strength
```

**Problem:** With uncentered biomarker (mean about 124) and
$t_{od}$ scaling with weeks (up to 10-20), the moderation shift at
interaction strength 0.45 was about $124 \times 10 \times 0.45 = 558$
units at $t_{od} = 10$, and larger at longer exposure (the earlier
"~650" corresponds to $t_{od} \approx 11.6$; made consistent with
`docs/23` on 2026-10-07), against BR noise SD of 8. Power was 1.00
everywhere. This was an artifact of a grossly overpowered signal, not
genuine resilience to carryover. The original
2025 codebase had this same scaling issue.

(Corrected 2026-10-08: the figure of about 558 units is the total
shift at the mean biomarker, and almost all of it is biomarker-free
drift, not moderation. Writing $B_i = \mu_{bm} + \sigma_{bm} b_i$, the
shift is $\beta_{bm}\mu_{bm}t_{od} + \beta_{bm}\sigma_{bm}t_{od}\,b_i$.
At $\beta_{bm} = 0.45$, $\mu_{bm} = 124.33$, $\sigma_{bm} = 15.36$ and
$t_{od} = 10$ the drift is 559.5, which the drug and time terms absorb
and which carries no interaction information; the moderation that the
interaction term sees is 69.1 per SD of the biomarker, about 8.6
$\sigma_{BR}$ (and 19.2 times the calibrated 3.6), growing with
$t_{od}$. The comparison of 558 with the noise SD therefore overstates
the signal-to-noise ratio by roughly a factor of 8; the conclusion, a
grossly overpowered signal, stands. Verified by direct computation. See
`analysis/report/compendium-mathematics.tex`, subsection "Historical
mean-moderation formulations and the calibration check", Equation
`eq:inva-attempt3`.)

### Attempt 4: Centered, on-drug binary, calibrated (correct)

```r
bm_z <- (dat$bm - bm_mean) / bm_sd
for (tp in 1:nP) {
  if (d[tp]$onDrug) {
    dat[, (br_col) := get(br_col) + c.bm * bm_z * br_sd]
  }
}
```

In mathematical form, the shift applied at every on-drug visit is
$\beta_{bm}\, b\, \sigma_{BR}$, with $\beta_{bm}$ passed in code as
`c.bm`.

**Why this works:** At $\beta_{bm} = c_{bm}$ the formulation matches
Architecture B's conditional expectation from the MVN model:

$$E[BR \mid B, \text{on drug}] = \mu_{BR} + c_{bm} \cdot
\frac{\sigma_{BR}}{\sigma_{bm}} \cdot (B - \mu_{bm})
= \mu_{BR} + c_{bm}\, \sigma_{BR}\, b$$

The standardized biomarker $b$ (code `bm_z`) centers the signal
(mean moderation = 0 in the population). The $\sigma_{BR}$ (code
`br_sd`) scaling produces comparable effect sizes between
architectures: a 1-SD biomarker shift produces $\beta_{bm}\sigma_{BR}$
units of BR shift, 3.6 at $\beta_{bm} = 0.45$ and $\sigma_{BR} = 8$.
The binary on-drug indicator ensures moderation is present whenever
drug is active and absent otherwise.

**Calibration ($c_{bm} = \beta_{bm} = 0.45$, extracted parameters):**

- Architecture A: 3.36 BR units per SD of biomarker
- Architecture B: 3.48 BR units per SD of biomarker
- Signal-to-noise ratio: about 0.65-0.9

Both empirical shifts are within sampling error of the theoretical 3.6
(the `.tex` version gives the 500-participant crossover draw behind
them). The signal-to-noise definition behind "0.65-0.9" is not stated,
and the range was not reproduced. The ratio of the per-SD shift to
$\sigma_{BR}$ is 0.45 (flagged 2026-10-07).

**Power (N-of-1, $N = 70$, $c_{bm} = \beta_{bm} = 0.45$, 50 replicates):**

| $t_{1/2}$ | Architecture A | Architecture B |
|:---:|:---:|:---:|
| 0.0 | 0.74 | 0.82 |
| 0.5 | 0.68 | 0.64 |
| 1.0 | 0.72 | 0.50 |

Table: N-of-1 power by $t_{1/2}$ for Architecture A and Architecture B

Architecture A power is stable across carryover levels
(variation is Monte Carlo noise at 50 replicates, MCSE about 0.065).
Architecture B shows the characteristic 40% relative decline
(0.82 to 0.50).

*Note (2026-10-07).* The 0.5 and 1.0 Architecture A cells disagree
with `docs/21`, which reports 0.72 and 0.68 where this note (and its
`.tex` version) report 0.68 and 0.72. The source run was not located,
so the discrepancy is unresolved. Either ordering lies within Monte
Carlo error of flat power. Whether $N = 70$ here is a total or a
per-path count was not re-verified (NOTATION rule 4). Paper 01's later
re-run at a matched total of $N = 70$ gives Hybrid carryover losses of
6.9 percent (mean moderation) and 15.4 percent (covariance
moderation), as recorded in `analysis/report/NOTATION.md` (Part 2).
`docs/37` derives the result that, under on-drug-only mean
moderation, power does not depend on the half-life. It also shows
that under graded covariance coupling power declines monotonically
with the half-life.

## Architectural Difference

The key distinction between Architectures A and B is where
the interaction signal lives:

- **Architecture B (MVN):** Signal is in the second moment
  (covariance). The analysis model infers the interaction
  from conditional correlations between bm and BR. Carryover
  erodes the off-drug correlation, compressing the contrast.

- **Architecture A (mean moderation):** Signal is in the
  first moment (mean). The moderation is a fixed additive
  shift to BR during on-drug timepoints. Off-drug timepoints
  carry no moderation signal regardless of carryover state.
  The analysis model detects the interaction from the
  differential mean between high-biomarker and
  low-biomarker participants when on drug.

## Implementation Details

The `dgp_architecture` parameter is threaded through three
functions:

- `buildSigma()`: when
  `dgp_architecture = "mean_moderation"`, the BM-BR
  correlation block is omitted from the covariance matrix.
- `generateData()`: after the MVN draw, the additive
  moderation `c.bm * bm_z * br_sd` is applied to the BR
  column at every on-drug timepoint (`tod > 0`), where `bm_z` is the
  standardized biomarker and `br_sd` is the BR standard
  deviation from the response parameters.
- `generateSimulatedResults()`: passes the parameter to
  both the sigma cache and the per-paramset data generation.

Default is `"mvn"` to preserve backward compatibility with
the Hendrickson et al. (2020) publication. All 67 existing
tests passed under the default at the time of the merge.

*Note (2026-10-07, inspected `R/generateData.R`).* The package has
since gained a third value, `dgp_architecture = "combined"` (paper
11). It applies the mean shift with weight `c.bm_a` and the covariance
coupling with weight `c.bm_b`. A `moderation_scaling` argument was
also added. Its default `"constant"` reproduces the shift above, and
`"trajectory"` scales it by the on-drug BR mean curve, normalized to
mean magnitude $\sigma_{BR}$. `validateParameterGrid()` is still not
threaded. It calls `buildSigma()` with the default `"mvn"`, so under
mean moderation it can report spurious PD failures (the known
limitation described in the `.tex` version).

## Files Modified

- `R/generateData.R`: `buildSigma()` and `generateData()`
- `R/generateSimulatedResults.R`: parameter threading

## Status

Merged to `master` on 2026-04-08 in commit `5ed7b41` (the
`feature/mean-moderation-dgp` branch). That commit is not present in
the current repository's history, which begins on 2026-05-01, so this
statement is unverified here. The same merge also
introduced the `implementations/` directory with three
parallel collections (`original`, `original-extended`,
`tidyverse`); the active `R/` package directory mirrors
`implementations/original-extended/`.
