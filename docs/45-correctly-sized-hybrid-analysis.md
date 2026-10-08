# Choosing a Correctly Sized Analysis for the Biomarker-by-Drug Interaction in a Hybrid N-of-1 Trial {.unlisted .unnumbered}
*2026-10-06 17:36 PDT*

**Author.** pmsimstats team

**Purpose.** This paper chooses a pre-specifiable primary analysis, and a
set of sensitivity analyses, for the biomarker-by-drug interaction in a
Hybrid N-of-1 trial of prazosin with blood pressure as the candidate
predictive biomarker. The criterion is that the test must hold its
nominal size before power is considered.

The evidence is the covariance study (`08-covariance-study.R`), run in
three passes on the same simulated trials: 63 analyses, 16 cells, and
378,000 model fits, with one failed fit. It also uses a separate check
of the drug-exposure coding used in manuscript 02
(`10-dbc-preexposure-check.R`). The paper consolidates the discussion
in `docs/40` to `docs/43` and corrects several statements made there
before the results were in.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. Symbols marked *local*
are defined here or in `docs/41`.

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $Y_{it}$ | symptom score of participant $i$ at visit $t$ (code `Sx`); $t = 0$ is baseline | canonical |
| $\mathrm{BL}_i$ | baseline symptom score | canonical |
| $B_i$ | biomarker (code `bm`); centered at its sample mean in every analysis | canonical |
| $D_{it}$ | binary drug state (code `Db`); 0 at baseline | canonical |
| $D_{bc,it}$ | exposure-decayed drug indicator (code `Dbc`): 1 on drug, $2^{-t_{sd}/t_{1/2}^{a}}$ off drug after exposure, 0 before any exposure | canonical; the pre-exposure value is made explicit here |
| $t_{sd}$ | time since discontinuation, in weeks; 0 before first exposure | canonical |
| $t_{1/2}$ | carryover half-life of the data-generating process, in weeks | canonical |
| $t_{1/2}^{a}$ | half-life assumed by the analysis when it builds $D_{bc,it}$ | local |
| $\lambda$ | decay rate, $\ln 2 / t_{1/2}$ | canonical |
| $c_{bm}$ | covariance-moderation strength: on-drug correlation of the biomarker with the drug response | canonical |
| $\rho$ | AR(1) serial correlation per week | canonical |
| $\beta_B$ | biomarker main effect | local (`docs/41`) |
| $\beta_D$, $\beta_{bm:D}$ | drug main effect; biomarker-by-drug interaction, the estimand | canonical |
| $\alpha$ | nominal level, 0.05 | canonical |
| $\pi$ | power, $\Pr(p < \alpha)$ | canonical |
| $\kappa$ | mean model-based standard error over the empirical standard deviation of the estimate; below 1 means anticonservative | canonical |
| $\nu$ | denominator degrees of freedom | canonical |
| $\Sigma$ | covariance matrix of a participant's outcomes | canonical |
| $N$ | total sample size, 70 | canonical |

Table: Notation index and glossary of symbols, meanings and status

**Sign convention.** The response components reduce symptom severity,
so $\beta_D$ and $\beta_{bm:D}$ are negative and $c_{bm}$ is positive.

**Analysis labels.** Each analysis is named by its baseline coding,
working covariance and inference, for example `BL+us-KR`.

- **Baseline coding.** No prefix: baseline is a ninth row of the
  response, coded off drug at $t = 0$ (the published coding; `docs/41`
  calls it the response coding). `BL+`: baseline is a covariate and the
  eight post-baseline visits are the response.
- **Working covariance.**
  - `cs`: compound symmetry; with a row-coded baseline and linear time
    this is the published random-intercept model.
  - `csh`: heterogeneous compound symmetry.
  - `us`: unstructured; 45 parameters over 9 rows, 36 over 8.
  - `sp_exp`: exponential spatial correlation in continuous time.
  - `RI+CAR1`: a random intercept plus `corCAR1` residuals, fitted with
    `nlme`. This is manuscript 02's analysis, G6 with $D_{it}$ and G8
    with $D_{bc,it}$.
- **Inference.**
  - `Satt`: asymptotic model-based covariance of the estimates, with
    Satterthwaite $\nu$.
  - `KR`: Kenward-Roger [2], with the `Kenward-Roger-Linear` covariance
    for `us`.
  - `CR2`: the bias-reduced cluster-robust (sandwich) covariance, with
    Satterthwaite $\nu$ (Bell and McCaffrey [3]).
- **Suffixes.** `-Dbc0.5`, `-Dbc1` and `-Dbc2` replace $D_{it}$ with
  $D_{bc,it}$ at $t_{1/2}^{a} = 0.5$, 1 and 2 weeks. `-sepbm` frees the
  biomarker main effect at baseline (Section 6).
- **Time adjustment.** Linear $t$ unless stated; `spline` is a natural
  spline with 3 degrees of freedom, and `visit` is a visit factor.

**Glossary.**

- **Size.** The rejection rate of the interaction test when
  $\beta_{bm:D} = 0$.
- **Correctly sized.** Defined in Section 2.
- **Size-adjusted power.** Power when each analysis rejects at the
  5% quantile of its own null p-values in the matching null cell, so
  that analyses are compared at a common realized size.
- **Construct.** The data-generating process. *CS* is the published
  compound-symmetry construct with graded coupling at $c_{bm} = 0.25$.
  *C07* is configuration C (separable AR(1), $\rho = 0.7$) with graded
  coupling at $c_{bm} = 0.45$ (`docs/36`, `docs/44`).
- **Shared biomarker effect.** In the row coding, a single $\beta_B$
  describes the biomarker's association with the outcome at baseline
  and at every post-baseline off-drug visit.
- **Constancy assumption.** The assumption that this association is in
  fact the same at baseline and at post-baseline off-drug visits
  (Section 6).

## 1. Summary

1. **The working covariance and the small-sample correction determine
   the size.** The time adjustment, the baseline coding and the exposure
   coding move the mean Type I error by at most about 0.01 (Section 4).
2. **Kenward-Roger and CR2 repair different failures, so neither
   substitutes for the other.**
   - Kenward-Roger corrects for the cost of estimating a covariance that
     is correctly specified. With `us` it brings the mean size from 0.090
     to 0.058. With a misspecified `cs` or `csh` it does nothing: 0.15 to
     0.16 in CO under AR(1) data.
   - CR2 corrects for misspecification of a parsimonious working
     covariance. With `cs`, `csh`, `sp_exp` or `RI+CAR1` it gives a mean
     size of 0.044 to 0.054 in every construct. With `us` it fails: 0.086,
     no better than no correction.
3. **Correctly sized in every one of the eight null cells (Section 4):**
   - nominal: `cs-CR2`, `csh-CR2`, `sp_exp-CR2`, `RI+CAR1-CR2` and
     `BL+cs-CR2`, mean 0.044 to 0.054;
   - slightly liberal: `us-KR` and `BL+us-KR`, mean 0.058 to 0.059, no
     cell above 0.072.

   The published analysis (`cs-Satt`) is correctly sized under CS data
   and in Hybrid. Under AR(1) data in CO it rejects 0.15 to 0.16 of the
   time.
4. **At a common size, power differs little under CS data.** In Hybrid
   without carryover, every correctly sized analysis except `sp_exp` and
   `RI+CAR1` has size-adjusted power of 0.66 to 0.72. Under AR(1) data
   the row-coded `us`, `RI+CAR1` and `sp_exp` analyses reach 0.94. The
   correctly sized baseline-covariate analyses reach 0.76 to 0.78
   (Section 5).
5. **The row coding's advantage under AR(1) comes from the shared
   biomarker effect, not from the covariance.** Freeing $\beta_B$ at
   baseline (`us-KR-sepbm`) gives the same power as baseline as a
   covariate: 0.780 against 0.792 in Hybrid. The advantage is valid only
   under the constancy assumption, which the null cells cannot test
   because the construct satisfies it by design (Section 6).
6. **$D_{bc,it}$ at a short assumed half-life is a near-free insurance
   against carryover.** With `us-KR` in Hybrid at $t_{1/2}^{a} = 0.5$:
   - without carryover it costs at most 0.016 in power;
   - at $t_{1/2} = 1$ week it gains 0.04 to 0.06;
   - its size is unchanged.

   A longer assumed half-life gains more under long carryover but costs
   up to 0.14 without it. With $t_{1/2}^{a} = t_{1/2}$ the estimate
   recovers the no-carryover interaction (Section 7).
7. **Manuscript 02's $D_{bc,it}$ codes off-drug visits before the first
   exposure as 1.** In CO this makes the indicator constant for every
   participant randomized to placebo first. On the same data, the
   corrected coding restores power from 0.484 to 0.744 under CS and from
   0.256 to 0.556 under AR(1) (verified, Section 8). This is the likely
   cause of manuscript 02's finding that the exposure-weighted analysis
   collapses under CO.
8. **Recommendation (Section 9).**
   - **Primary analysis:** an MMRM with unstructured covariance and
     Kenward-Roger inference, baseline as a response row, linear time,
     and $D_{bc,it}$ at a pre-specified $t_{1/2}^{a} = 0.5$ week.
   - **Pre-specified sensitivity analyses:** the same model with
     $\beta_B$ freed at baseline (equivalently, baseline as a
     covariate), and `RI+CAR1-CR2`.
   - **If the protocol cannot defend the constancy assumption,** the
     freed model should be primary.

## 2. The question and the criterion

A confirmatory test of a predictive biomarker is useful only if its
Type I error is what the protocol states. In this setting, an analysis
that rejects 8% of the time under the null cannot be rescued by its
power, so the analyses are first screened on size, and power is
compared only among those that pass.

**Monte Carlo precision.** Each null cell has 500 replicates, so the
Monte Carlo standard error (MCSE) of a rejection rate of 0.05 is 0.0097.
A single cell between 0.031 and 0.069 is consistent with nominal size.
With 8 null cells, the pooled rate rests on 4,000 replicates (MCSE
0.0034), with a 95% band of 0.043 to 0.057.

**Criterion.** An analysis is called:

- **nominal** if its mean size over the eight null cells is within
  0.043 to 0.057 and no cell exceeds 0.075;
- **slightly liberal** if the mean is 0.057 to 0.065 and no cell exceeds
  0.075;
- **conservative** if the mean is below 0.043;
- **liberal** otherwise.

An analysis is *correctly sized* if it is nominal or slightly liberal.
The upper cell limit of 0.075 is 2.6 MCSE above nominal. With 63
analyses and 8 cells, a few individual cells above 0.069 are expected by
chance, so the cell limit is set above the single-cell band.

## 3. The study

**Constructs.** Both constructs use the published means and standard
deviations of `04-power-simulation.R` (which reproduces the `58b32a9`
matrices), $N = 70$, and the published allocation (18/18/17/17 in
Hybrid, 35/35 in CO). Both use graded coupling: full coupling on drug,
$c_{bm} 2^{-t_{sd}/t_{1/2}}$ off drug after exposure, zero before
exposure (`docs/37`, Appendix A).

- **CS:** compound symmetry within each response component (0.8 at
  every lag), cross-component correlations 0.2 and 0.1, $c_{bm} = 0.25$
  (below the ceiling of 0.256).
- **C07:** separable AR(1) at $\rho = 0.7$ per week, $c_{bm} = 0.45$
  (below the ceiling of 0.480 in Hybrid and 0.616 in CO).

**Designs.**

- **Hybrid:** visits at weeks 4, 8, 9, 10, 11, 12, 16 and 20 after
  baseline, with four paths: discontinuation after week 9 or week 10,
  crossed with drug-first or placebo-first crossover.
- **CO:** eight visits 2.5 weeks apart, with drug-first and
  placebo-first paths.

**Cells.** 2 constructs × 2 designs × $t_{1/2} \in \{0, 1\}$ week ×
{alternative, null}, 16 cells. The null twin of each cell sets
$c_{bm} = 0$. Alternative cells have 250 replicates (power MCSE at most
0.032); null cells have 500.

**Analyses.** The mean model is $\beta_0 + \beta_B B_i + \beta_D X_{it}
+ f(t) + \beta_{bm:D} B_i X_{it}$, with $X_{it}$ either $D_{it}$ or
$D_{bc,it}$ and $f(t)$ linear, spline or visit factor. The baseline
coding, working covariance and inference vary as in the glossary. All
analyses are fitted with `mmrm`, except `RI+CAR1` (`nlme` plus
`clubSandwich`).

**Three passes, paired.**

| Pass | Analyses | Fits | Failed fits |
|----------|------------------------------------------------------|---------|-------|
| Main | `cs`, `csh`, `us`, `sp_exp` with Satterthwaite; `us-KR`; `BL+cs`, `BL+us` with Satterthwaite; each time form | 120,000 | 0 |
| Follow-up | KR and CR2 versions of the above; `us-KR-sepbm` | 138,000 | 0 |
| Exposure | $D_{bc,it}$ at three $t_{1/2}^{a}$ under `us-KR`, `us-CR2` and `RI+CAR1-CR2`; $D_{it}$ under `RI+CAR1-CR2` | 120,000 | 1 |

Table: Covariance study passes with analyses, model fits and failed fits

Every pass uses the same per-chunk seeds, so the analyses are fitted to
identical simulated trials. Differences between analyses are therefore
paired comparisons, and are more precise than the MCSE of each rate
suggests.

## 4. Size

### 4.1 The size table

Mean, minimum and maximum size over the eight null cells (linear time),
with the maximum in each group of cells: the four CS cells, the two C07
Hybrid cells and the two C07 CO cells. The last numeric column is the
smallest $\kappa$ over the eight cells.

| Analysis | Mean | Min | Max, CS | Max, C07 Hybrid | Max, C07 CO | Min $\kappa$ | Class |
|-------------|------|------|------|------|------|------|----------|
| `cs-Satt` | 0.079 | 0.042 | 0.058 | 0.066 | 0.158 | 0.73 | liberal |
| `cs-KR` | 0.079 | 0.042 | 0.058 | 0.066 | 0.158 | 0.73 | liberal |
| `csh-KR` | 0.073 | 0.040 | 0.058 | 0.048 | 0.154 | 0.75 | liberal |
| `cs-CR2` | 0.051 | 0.042 | 0.064 | 0.044 | 0.060 | 0.95 | nominal |
| `csh-CR2` | 0.054 | 0.048 | 0.064 | 0.054 | 0.056 | 0.93 | nominal |
| `sp_exp-Satt` | 0.042 | 0.028 | 0.062 | 0.048 | 0.034 | 0.95 | conservative |
| `sp_exp-CR2` | 0.048 | 0.040 | 0.064 | 0.050 | 0.044 | 0.92 | nominal |
| `RI+CAR1-CR2` | 0.051 | 0.040 | 0.064 | 0.056 | 0.050 | 0.94 | nominal |
| `us-Satt` | 0.090 | 0.072 | 0.092 | 0.102 | 0.104 | 0.84 | liberal |
| `us-KR` | 0.058 | 0.052 | 0.062 | 0.066 | 0.056 | 0.94 | slightly liberal |
| `us-CR2` | 0.086 | 0.070 | 0.092 | 0.100 | 0.088 | 0.83 | liberal |
| `BL+cs-KR` | 0.077 | 0.038 | 0.052 | 0.070 | 0.158 | 0.73 | liberal |
| `BL+cs-CR2` | 0.048 | 0.036 | 0.052 | 0.050 | 0.054 | 0.97 | nominal |
| `BL+us-Satt` | 0.087 | 0.072 | 0.094 | 0.092 | 0.094 | 0.86 | liberal |
| `BL+us-KR` | 0.059 | 0.050 | 0.068 | 0.060 | 0.066 | 0.96 | slightly liberal |
| `BL+us-CR2` | 0.080 | 0.068 | 0.086 | 0.090 | 0.086 | 0.86 | liberal |

Table: Null size and minimum $\kappa$ by analysis over the eight null cells

`cs-Satt` is the published analysis. All figures are verified
(`09-covariance-study-tables.R`). The visit-factor and spline versions
fall in the same classes, with mean sizes within 0.005 of the linear
versions. `us-KR` is slightly liberal by about 0.008, a little over 2
MCSE of the pooled rate. The excess appears in every `us-KR` variant,
but they are fitted to the same trials and are not independent
confirmations.

### 4.2 Why Kenward-Roger and CR2 are not interchangeable

Kenward-Roger [2] makes two adjustments. It inflates the model-based
covariance of the fixed effects for the uncertainty in the estimated
covariance parameters, and it matches the denominator degrees of
freedom to the moments of the resulting Wald statistic. Both
adjustments assume the working covariance is correct.

- **`us` is correct by construction.** Its only failure is the cost of
  estimating 45 parameters from 70 participants. That is what
  Kenward-Roger repairs: $\kappa$ rises from 0.84-0.90 to 0.94-1.01, and
  the mean size falls from 0.090 to 0.058.
- **`cs` and `csh` under AR(1) data are misspecified.** Their failure is
  bias in the model-based variance, which Kenward-Roger does not
  address: $\kappa$ stays at 0.73-0.75 in CO.

CR2 [3] replaces the model-based variance with a sandwich whose
residuals are adjusted by the working model's hat matrix. It is
consistent whatever the working covariance.

- **With a parsimonious working model** (2 to 3 parameters), the
  weights are nearly fixed and the bias reduction works. $\kappa$ is
  0.92-1.04, and every such analysis is nominal.
- **With `us`, the weights are themselves estimated from 45
  parameters.** They adapt to the noise in each replicate, and the
  sandwich does not account for this. $\kappa$ is 0.83-0.89, the same as
  the uncorrected model-based variance (an inferred mechanism; the
  $\kappa$ values are verified).

The practical rule is to pair `us` with Kenward-Roger and a
parsimonious working covariance with CR2.

### 4.3 Why the published model fails in CO under AR(1) data

The CO contrast compares two blocks of four visits, 2.5 weeks apart
within a block, whose centers are ten weeks apart. Under AR(1) at
$\rho = 0.7$ per week, the average correlation within a block is 0.272
and between blocks 0.069.

A compound-symmetry fit estimates one common correlation, the average
over all 28 pairs, 0.156. For the difference of block means, in units
of the per-visit variance:

- the true variance is $2(1 + 3 \times 0.272)/4 - 2 \times 0.069 =
  0.770$;
- the compound-symmetry variance is $2(1 - 0.156)/4 = 0.422$.

The implied $\kappa$ is $\sqrt{0.422/0.770} = 0.740$, against the
observed 0.73 to 0.75 (derived; approximate, because it ignores the
baseline row, the time term and the cross-component structure). In
Hybrid most off-drug visits are one week from an on-drug visit, so the
averaging does less harm, and $\kappa$ is 0.92 to 0.97.

The failure therefore depends on the design. It would appear in any
design that contrasts long, well-separated blocks under serial
correlation.

### 4.4 The time adjustment

`docs/40`, `docs/42` and `docs/43` reported an unexplained excess size,
about 0.07, for flexible time adjustments in Hybrid under carryover.

| Analysis | Linear | Spline | Visit |
|---|---|---|---|
| `cs-Satt` | 0.054 | 0.070 | 0.072 |
| `us-KR` | 0.056 | 0.068 | 0.072 |
| `BL+us-KR` | 0.056 | 0.060 | 0.064 |
| `cs-CR2` | 0.058 | | 0.058 |
| `csh-CR2` | 0.060 | | 0.060 |

Table: Null rejection rate in the CS Hybrid cell by time adjustment

*Null rejection rate in the CS Hybrid cell at $t_{1/2} = 1$, by time
adjustment (500 replicates, MCSE 0.0097).*

The excess is confined to this one cell. It is 2.0 to 2.3 MCSE above
nominal and absent under CR2. Averaged over the eight cells, the visit
factor and the spline have the same mean size as linear time in every
analysis (0.058 for `us-KR`). The evidence for a systematic effect of
the time form on size is therefore weak. Linear time is never worse in
any cell, and it is the recommended form.

## 5. Power among the correctly sized analyses

Size-adjusted power, linear time, in Hybrid (verified):

| Analysis | CS, $t_{1/2}=0$ | CS, $t_{1/2}=1$ | C07, $t_{1/2}=0$ | C07, $t_{1/2}=1$ |
|---|---|---|---|---|
| `us-KR` | 0.688 | 0.424 | 0.940 | 0.756 |
| `us-KR-Dbc0.5` | 0.700 | 0.472 | 0.908 | 0.812 |
| `RI+CAR1-CR2` | 0.616 | 0.456 | 0.940 | 0.692 |
| `cs-CR2` | 0.664 | 0.428 | 0.828 | 0.676 |
| `csh-CR2` | 0.692 | 0.392 | 0.840 | 0.712 |
| `sp_exp-CR2` | 0.476 | 0.332 | 0.936 | 0.664 |
| `BL+us-KR` | 0.704 | 0.404 | 0.784 | 0.464 |
| `BL+cs-CR2` | 0.724 | 0.476 | 0.760 | 0.468 |
| `cs-Satt` (published) | 0.664 | 0.456 | 0.844 | 0.676 |

Table: Size-adjusted power of correctly sized analyses in the Hybrid design

And in CO:

| Analysis | CS, $t_{1/2}=0$ | CS, $t_{1/2}=1$ | C07, $t_{1/2}=0$ | C07, $t_{1/2}=1$ |
|---|---|---|---|---|
| `us-KR` | 0.720 | 0.748 | 0.804 | 0.768 |
| `us-KR-Dbc0.5` | 0.716 | 0.748 | 0.796 | 0.776 |
| `RI+CAR1-CR2` | 0.684 | 0.704 | 0.688 | 0.732 |
| `cs-CR2` | 0.712 | 0.724 | 0.592 | 0.640 |
| `csh-CR2` | 0.688 | 0.688 | 0.608 | 0.688 |
| `sp_exp-CR2` | 0.336 | 0.312 | 0.660 | 0.668 |
| `BL+us-KR` | 0.708 | 0.756 | 0.492 | 0.520 |
| `BL+cs-CR2` | 0.784 | 0.792 | 0.504 | 0.524 |

Table: Size-adjusted power of correctly sized analyses in the CO design

The size-adjusted figures rest on a critical value estimated from 500
null replicates. That adds about 0.02 to their uncertainty beyond the
power MCSE of up to 0.032, so differences under about 0.05 should not
be read as real.

- **Under CS data the choice hardly matters,** with two exceptions.
  `sp_exp` loses a third of the power in Hybrid and half in CO, because
  it imposes a decay the data do not have. `RI+CAR1` loses about 0.07
  in Hybrid.
- **Under AR(1) data the row-coded `us` analysis is the most powerful,
  or tied for it, in every cell.** `RI+CAR1` and `sp_exp` match it without carryover in
  Hybrid but fall behind with carryover and in CO.
- **Kenward-Roger costs no power.** At a common size, `us-KR` and
  `us-Satt` agree to within 0.008 in every cell (verified).
  Kenward-Roger changes the standard error and degrees of freedom, not
  the estimate, so the ordering of the replicates barely changes. The
  higher raw power of `us-Satt` (for example 0.956 against 0.944) is its
  excess size.
- **The baseline-covariate analyses lose 0.15 to 0.31 under AR(1)
  data.** Section 6 explains why.

## 6. The baseline row and the shared biomarker effect

In the row coding, the baseline row and the post-baseline off-drug
visits share one biomarker main effect $\beta_B$. In the covariate
coding, baseline is conditioned on and the biomarker's association with
it is not used. The follow-up pass added a single analysis to separate
the two explanations for the row coding's power: `us-KR-sepbm`, the row
coding with an extra term $B_i \cdot \mathbf{1}\{t = 0\}$ that frees
the biomarker effect at baseline.

| Analysis | CS Hybrid 0 | CS Hybrid 1 | C07 Hybrid 0 | C07 Hybrid 1 | C07 CO 0 | C07 CO 1 |
|----------------------|-----------|-----------|-----------|-----------|-----------|-----------|
| `us-KR` | 0.704 | 0.468 | 0.944 | 0.796 | 0.800 | 0.788 |
| `us-KR-sepbm` | 0.720 | 0.460 | 0.780 | 0.504 | 0.592 | 0.528 |
| `BL+us-KR` | 0.704 | 0.456 | 0.792 | 0.484 | 0.600 | 0.572 |

Table: Raw power of `us-KR` with shared versus freed baseline biomarker effect

*Raw power with the visit-factor time adjustment, 250 replicates per
cell (verified). Column labels give construct, design and $t_{1/2}$.
`us-KR` shares $\beta_B$ with the baseline row; `us-KR-sepbm` frees it.*

Freeing $\beta_B$ at baseline removes the whole advantage. The freed
row model and the covariate model then agree to within Monte Carlo
error in every cell. Their size is the same (0.058 and 0.059).

**Correction to `docs/41`.** Section 7 of `docs/41` derived that, with
complete data, an unstructured covariance and $\beta_B$ free at
baseline, keeping baseline in the response is equivalent to conditioning
on it. That derivation is confirmed. Its table, however, applied the
equivalence to the published row coding, which shares $\beta_B$. That
application is wrong: the shared form is a different, more constrained
model, and under AR(1) data it is substantially more powerful.

**What the constraint buys.** In Hybrid, the post-baseline off-drug
visits are few and late (two or three blinded-discontinuation visits and
one crossover visit). The interaction is the difference between the
biomarker's on-drug and off-drug associations. With the constraint,
the baseline row adds a further, well-measured off-drug observation of
the biomarker association. Under CS the participant effect absorbs most
of what this adds, which is consistent with Proposition 4 of `docs/41`.
Under AR(1) it does not. The mechanism is verified empirically by the
`sepbm` arm; an analytic account has not been derived.

**What the constraint assumes.** The constancy assumption holds in both
constructs by design: the biomarker is coupled only to the drug
response, never to the baseline or to the natural-course and
expectancy components. The null cells therefore cannot reveal what
happens when it fails.

In a prazosin trial it fails if blood pressure is prognostic for the
symptom course off drug: for example, if participants with higher
pressure improve more, or less, on placebo or with time. The shared
$\beta_B$ would then be a compromise between two different
associations. The interaction estimate would absorb the difference,
and its null distribution would no longer be centered at zero
(inferred, not simulated).

## 7. Carryover in the analysis: $D_{bc,it}$

`docs/43` asked whether the analysis should model carryover. The
exposure pass answers this for `us-KR`, `us-CR2` and `RI+CAR1-CR2`.

**Size.** $D_{bc,it}$ does not change the size of any analysis. `us-KR`
is 0.058 with $D_{it}$, and 0.057, 0.055 and 0.059 at
$t_{1/2}^{a} = 0.5$, 1 and 2. `RI+CAR1-CR2` is 0.051 to 0.054
throughout. `us-CR2` stays liberal (0.085) (verified).

**Power** (`us-KR`, linear time, raw; verified):

| Exposure coding | CS Hybrid 0 | CS Hybrid 1 | C07 Hybrid 0 | C07 Hybrid 1 | CO, all cells |
|----------------------------|-----------|-----------|-----------|-----------|-------------|
| $D_{it}$ | 0.700 | 0.436 | 0.944 | 0.812 | 0.752-0.808 |
| $D_{bc,it}$, $t_{1/2}^{a} = 0.5$ | 0.708 | 0.492 | 0.928 | 0.852 | 0.748-0.808 |
| $D_{bc,it}$, $t_{1/2}^{a} = 1$ | 0.632 | 0.504 | 0.892 | 0.884 | 0.744-0.808 |
| $D_{bc,it}$, $t_{1/2}^{a} = 2$ | 0.564 | 0.468 | 0.840 | 0.876 | 0.736-0.796 |

Table: Raw `us-KR` power by exposure coding and assumed half-life $t_{1/2}^{a}$

- **In Hybrid it is a trade-off.** An assumed half-life longer than the
  true one spreads the contrast onto off-drug visits that carry no
  effect. One shorter than the true one leaves residual effect in the
  off-drug reference.
- **At $t_{1/2}^{a} = 0.5$ the trade-off is favorable in every cell.**
  It costs 0.016 or less without carryover and gains 0.04 to 0.06 at
  $t_{1/2} = 1$. The size-adjusted figures agree (Section 5).
- **In CO it makes little difference** (0.736 to 0.808 in every cell
  and at every assumed half-life). The first off-drug visit is 2.5 weeks
  after discontinuation, where $D_{bc,it}$ is 0.03 at
  $t_{1/2}^{a} = 0.5$ but 0.42 at $t_{1/2}^{a} = 2$. Why even the longer
  assumptions cost nothing in CO was not analyzed.

**The estimand.** With $D_{it}$, carryover attenuates the estimate, for
example from $-0.135$ to $-0.100$ in CS Hybrid at $t_{1/2} = 1$.
$D_{bc,it}$ with $t_{1/2}^{a} = t_{1/2}$ recovers the full on-drug
interaction ($-0.129$). A longer assumed half-life overshoots it
($-0.151$ at $t_{1/2}^{a} = 2$; $-0.168$ without carryover). Under
$D_{bc,it}$, $\beta_{bm:D}$ is the interaction per unit of modeled
exposure, so its scale depends on $t_{1/2}^{a}$. The assumed half-life
must be pre-specified, and the estimate is interpretable as the on-drug
interaction only if the assumption is close to the truth.

**For prazosin.** The pharmacokinetic half-life of prazosin is a few
hours, so a pharmacological carryover half-life of 0.5 week is
generous. Carryover much longer than that would have to be behavioral
or psychological (inferred). An assumed half-life of 0.5 week therefore
covers the plausible pharmacological range, and costs almost nothing if
there is no carryover.

**`RI+CAR1-CR2` with $D_{bc,it}$** shows the same pattern at lower
power, with a steeper cost at long assumed half-lives (0.412 at
$t_{1/2}^{a} = 2$ in CS Hybrid without carryover, against 0.616 with
$D_{it}$).

## 8. The pre-exposure coding of $D_{bc,it}$ in manuscript 02

Manuscript 02 builds the exposure-weighted indicator as 1 on drug and
`carryover_decay(tsd, t_half)` off drug
(`analysis/scripts/carryover-sensitivity/simulation-core.R`, lines 179
to 183; the same pattern at lines 662 and 688). Two facts combine
(inspected):

- the package sets $t_{sd} = 0$ before the first exposure
  (`R/buildtrialdesign.R`, line 105);
- `carryover_decay(0, h)` returns $e^{0} = 1$ for any $h > 0$
  (`implementations/tidyverse/R/functions.R`, line 56, which the
  manuscript 02 driver uses; the `nof1power` copy behaves the same).

Off-drug visits that precede any exposure are therefore coded as fully
exposed. In Hybrid every path starts on drug, so the defect does not
arise; manuscript 02 also has no baseline row. In CO the placebo-first
path is off drug for its first four visits, and these are coded 1. The
indicator is then 1 at every visit for half the participants, and they
contribute nothing to the within-participant contrast.

**Check.** `10-dbc-preexposure-check.R` fits the same trials with both
codings, using manuscript 02's model (random intercept, `corCAR1`, no
baseline row, linear time) at $t_{1/2} = t_{1/2}^{a} = 1$, with 250
replicates per cell (verified):

| Construct, design | $D_{it}$ | $D_{bc,it}$ corrected | $D_{bc,it}$ as in 02 | Constant indicator |
|--------------------|----------|--------------|--------------|------------|
| CS, Hybrid | 0.428 | 0.496 | 0.496 | 0 |
| CS, CO | 0.744 | 0.744 | 0.484 | 35 |
| C07, Hybrid | 0.476 | 0.576 | 0.576 | 0 |
| C07, CO | 0.548 | 0.556 | 0.256 | 35 |

Table: Power under corrected and manuscript 02 pre-exposure coding of $D_{bc,it}$

*Power with CR2 standard errors. The last column counts participants
whose indicator is constant over all visits under manuscript 02's
coding.*

- **In Hybrid** the two codings give identical fits, as they must.
- **In CO** the defect costs 0.26 to 0.30 in power, and it raises the
  empirical standard deviation of the estimate by a factor of 1.44 to
  1.47 (from 0.049 to 0.070 under CS). That is close to the $\sqrt{2}$
  expected from losing half the participants' contrasts.
- **The corrected coding performs like $D_{it}$ in CO,** as in the
  covariance study.

Manuscript 02 reports exposure-weighted power of 0.488 against 0.830 for
the binary coding in CO (its Section 3.4). It attributes this to "the
continuous decayed predictor's mismatch with CO's design". The size of
the loss matches the defect, and the corrected coding shows no mismatch.
The reported finding is very likely an artifact of the coding (inferred:
manuscript 02's own pipeline has not been rerun with the correction).

The same applies to its finding that AIC selection of the half-life
picks the largest candidate in CO. That selection was made among
models all affected by the defect.

**Fix.** Code pre-exposure off-drug visits as 0, for example
`Dbc = ifelse(Db == 1, 1, ifelse(tsd > 0, carryover_decay(tsd, h), 0))`.
Then rerun the CO cells of manuscript 02, its blocks S7 to S10 in
particular.

## 9. Recommendation for a Hybrid trial of prazosin with blood pressure

### 9.1 Primary analysis

An MMRM for all nine rows (baseline as a row at $t = 0$, coded off
drug):

$$
Y_{it} = \beta_0 + \beta_B B_i + \beta_T t + \beta_D D_{bc,it}
  + \beta_{bm:D} B_i D_{bc,it} + \varepsilon_{it},
\qquad \operatorname{Var}(\boldsymbol\varepsilon_i) = \Sigma
  \text{ unstructured},
$$

with $B_i$ centered at its sample mean, $D_{bc,it}$ built at
$t_{1/2}^{a} = 0.5$ week with pre-exposure visits coded 0, REML
estimation, and Kenward-Roger inference. The test is the two-sided
Wald test of $\beta_{bm:D} = 0$ at $\alpha = 0.05$.

```{.r}
library(mmrm)
dat$Dbc <- ifelse(dat$Db == 1, 1,
                  ifelse(dat$tsd > 0, 0.5^(dat$tsd / 0.5), 0))
dat$bmc <- dat$bm - mean(dat$bm[dat$t == 0])
fit <- mmrm(Sx ~ bmc + t + Dbc + bmc:Dbc + us(visit | ptID),
            data = dat, method = 'Kenward-Roger',
            vcov = 'Kenward-Roger-Linear')
summary(fit)$coefficients['bmc:Dbc', ]
```

**Why this analysis.**

- **Size:** slightly liberal, 0.057 on average, never above 0.064 in
  any cell with linear time.
- **Robustness:** its validity does not depend on the serial
  correlation structure.
- **Power:** it is the most powerful correctly sized analysis under
  AR(1) data and as powerful as any under CS data.
- **Carryover:** at $t_{1/2}^{a} = 0.5$ it insures against
  pharmacological carryover at almost no cost.

**What the protocol must state.**

- the assumed half-life and its justification;
- the constancy assumption (Section 6), with the reason to believe it;
- linear time as the time adjustment;
- that the unstructured covariance is estimated over all nine rows;
- the handling of missing visits: MMRM under missing at random, which
  this study did not test (Section 11).

### 9.2 Pre-specified sensitivity analyses

1. **$\beta_B$ freed at baseline.** Add `bmc:bl0`, where `bl0` is 1 at
   baseline. Equivalently, analyze the eight post-baseline visits with
   baseline as a covariate (`BL+us-KR`). Its size is the same; it gives
   up the constancy assumption at a power cost of up to 0.15 under
   AR(1) data. A material disagreement between the primary and this
   analysis is evidence against the constancy assumption.
2. **`RI+CAR1-CR2` with the same $D_{bc,it}$.** A nominally sized
   analysis on a different basis of inference (a robust variance about
   a parsimonious working model). Agreement shows the conclusion does
   not rest on the unstructured covariance.
3. **$D_{it}$ in place of $D_{bc,it}$.** The published exposure coding,
   for comparability with Hendrickson et al. [1].

### 9.3 When to change the primary analysis

- **If the constancy assumption cannot be defended:** for example, if
  pilot data suggest blood pressure predicts placebo response or the
  natural course. Make sensitivity analysis 1 primary.
- **If the sample is much smaller than 70,** or visits are often
  missed, so that 45 covariance parameters are poorly determined: use
  `csh-CR2` or `RI+CAR1-CR2` as primary. Both are nominal in every
  cell, at a size-adjusted power cost of up to about 0.20 against the
  recommended primary under AR(1) data (least in Hybrid without
  carryover, most in CO). This study did not test smaller samples.
- **If the design is CO rather than Hybrid,** the same recommendation
  holds, and $D_{bc,it}$ makes no difference.

### 9.4 What not to use

- **The published random-intercept model with model-based inference
  (`cs-Satt`, equivalently `lmer`).** It is anticonservative under
  serial correlation in CO, 0.15 to 0.16.
- **`us` without Kenward-Roger, and `us` with CR2.** Both are liberal,
  0.080 to 0.090.
- **`cs` or `csh` with Kenward-Roger.** Same failure as the published
  model.
- **`sp_exp` with model-based inference.** It is conservative and loses
  a third of the power under CS data.
- **Manuscript 02's $D_{bc,it}$ coding** in any design with off-drug
  visits before the first exposure.
- **A long assumed half-life chosen for safety.** It costs up to 0.14 in
  power when carryover is short.

## 10. Corrections to earlier documents

| Document | Statement | Correction |
|----------------|------------------------------|----------------------------------------|
| `docs/41`, Section 7 table | cLDA and ANCOVA coincide for the unstructured row model | True only with $\beta_B$ free at baseline; the published row coding shares it and differs (Section 6) |
| `docs/41`, Section 8 | Recommends baseline as a covariate for the trial analysis | Superseded by Section 9: the covariate form is the sensitivity analysis, or the primary if constancy cannot be defended |
| `docs/42` and `docs/43`, Section 4 | The flexible-time size excess persists under every covariance and its cause is unexplained | Confined to one cell, 2 MCSE, absent under CR2; no systematic effect (Section 4.4) |
| `docs/43`, Sections 1 and 5 | `us-KR` is sized at 0.046 to 0.072 | Slightly liberal, mean 0.058; the CR2 analyses with parsimonious covariances are nominal (Section 4) |
| `docs/43`, Section 6 | The follow-up questions are open | Answered here (Sections 4 to 7) |
| Manuscript 02, Section 3.4 and blocks S7 to S10 | The exposure-weighted analysis collapses under CO | Likely an artifact of the pre-exposure coding (Section 8) |

Table: Corrections to statements in earlier documents

## 11. Limitations and evidence status

- **Two constructs, both satisfying the constancy assumption.** No null
  cell has a prognostic biomarker, so the row coding's size when the
  assumption fails is unknown. This is the most important open
  question for the recommendation.
- **One sample size, no missing data.** $N = 70$ and complete data
  throughout. Kenward-Roger with `us` is known to degrade in small
  samples, and missing visits change both the covariance estimation
  and the cLDA-ANCOVA comparison (`docs/41`, Section 7).
- **Two half-lives.** $t_{1/2} \in \{0, 1\}$ week; the published values
  of 0.1 and 0.2 lie between them. The claim that $D_{bc,it}$ at 0.5
  costs almost nothing at short carryover is interpolated.
- **Untested combinations.** $D_{bc,it}$ was crossed only with row-coded
  `us-KR`, `us-CR2` and `RI+CAR1-CR2`. Its size under the freed
  ($\beta_B$ at baseline) and covariate codings is inferred from its
  invariance in the row coding. The expectancy term, random slopes and
  GEE were not part of this study (`docs/43`).
- **Monte Carlo precision.** Power MCSE is up to 0.032 per cell. Size
  MCSE is 0.0097 per cell and 0.0034 pooled. Classes near a boundary
  (`us-KR` at 0.058) could change with more replicates.
- **Mechanisms.** The account of the CR2 failure with `us` is inferred.
  The CS-under-AR(1) account (Section 4.3) is an approximate derivation
  that agrees with the observed $\kappa$. The row-coding mechanism is
  verified empirically but not derived.
- **Manuscript 02.** The defect is verified by inspection and its effect
  by simulation in this construct; manuscript 02's own pipeline has not
  been rerun.

## 12. Open questions

1. **Size of the row coding under a prognostic biomarker.** Couple the
   biomarker to the natural-course or expectancy component in the null
   cells. This decides between the primary and sensitivity analysis 1.
2. **The recommended primary and sensitivity analyses in the remaining
   conditions.** Run them at $t_{1/2} \in \{0.1, 0.2, 0.5\}$, at
   $N \in \{30, 50\}$, and with the published dropout mechanisms.
3. **Manuscript 02 with the corrected coding.** Rerun its CO cells.
4. **An analytic account of the shared-$\beta_B$ gain under AR(1),**
   along the lines of `docs/41`, Proposition 4.

## 13. Reproducibility

Run from the repository root:

```bash
Rscript analysis/scripts/quick-sim/carryover-closed-form/08-covariance-study.R \
  --reps-alt 250 --reps-null 500 --cores 8
Rscript analysis/scripts/quick-sim/carryover-closed-form/08-covariance-study.R \
  --reps-alt 250 --reps-null 500 --cores 8 --arms followup
Rscript analysis/scripts/quick-sim/carryover-closed-form/08-covariance-study.R \
  --reps-alt 250 --reps-null 500 --cores 8 --arms dbc
Rscript analysis/scripts/quick-sim/carryover-closed-form/09-covariance-study-tables.R
Rscript analysis/scripts/quick-sim/carryover-closed-form/10-dbc-preexposure-check.R \
  --reps 250 --cores 4
```

- **Study passes (the first three commands).** Each writes
  `summary-alt250-null500{,-followup,-dbc}.csv` and the matching
  replicate files to
  `analysis/data/quick-sim/carryover-closed-form/covariance-study/`.
  Chunks are checkpointed under `cells-*`, so an interrupted pass
  resumes.
- **Tables (`09`).** Combines the passes into
  `summary-alt250-null500-combined.csv` and computes
  `size-adjusted-power-alt250-null500.csv`.
- **Pre-exposure check (`10`).** Writes to
  `analysis/data/quick-sim/carryover-closed-form/dbc-preexposure/`.
- **Run time.** On 8 cores the three passes took about 3, 25 and 23
  hours of wall time. The last two include periods when the machine
  slept. The check took about 15 minutes.

## 14. References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020; 2:13. doi:10.3389/fdgth.2020.00013.
2. Kenward MG, Roger JH. Small sample inference for fixed effects from
   restricted maximum likelihood. *Biometrics* 1997; 53(3):983-997.
   doi:10.2307/2533558.
3. Bell RM, McCaffrey DF. Bias reduction in standard errors for linear
   regression with multi-stage samples. *Survey Methodology* 2002;
   28(2):169-181. Pustejovsky JE, Tipton E. Small-sample methods for
   cluster-robust variance estimation and hypothesis testing in fixed
   effects models. *Journal of Business and Economic Statistics* 2018;
   36(4):672-683. doi:10.1080/07350015.2016.1247004.
4. Sabanés Bové D, et al. *mmrm: Mixed Models for Repeated Measures*. R
   package, CRAN.
5. Liu GF, Lu K, Mogg R, Mallick M, Mehrotra DV. Should baseline be a
   covariate or dependent variable in analyses of change from baseline
   in clinical trials? *Statistics in Medicine* 2009; 28(20):2509-2530.
   doi:10.1002/sim.3639.
6. pmsimstats team. `docs/37-carryover-power-closed-form.md`,
   `docs/40-hendrickson-hybrid-design-primer.md`,
   `docs/41-baseline-as-covariate.md`,
   `docs/42-flexible-time-adjustment.md`,
   `docs/43-analysis-approaches-catalog.md`, `docs/44-dgp-roadmap.md`.
