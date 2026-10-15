# Appendix I: The Legacy Record (docs/00 to docs/20)

*2026-10-15 11:00 PDT*

**Author:** pmsimstats team
**Series position:** Appendix I of the whitepaper series
**Purpose:** Results of docs/00 to docs/20 (with 02, 09a and 19)
not carried in the main series chapters. Each section opens with
a provenance line and a status label; none restates a superseded
mechanism except to name it superseded.

---

## Contents

1. [Section 0: Purpose, scope, and map](#section-0)
2. [Section 1: Commit history and PD audit](#section-1)
3. [Section 2: Published simulation grid](#section-2)
4. [Section 3: Published Figures 4 and 5 reproductions](#section-3)
5. [Section 4: Analysis-model switch and Type I error](#section-4)
6. [Section 5: Mean-moderation formulations that failed](#section-5)
7. [Section 6: Constant-exposure study and crossover bias](#section-6)
8. [Section 7: Open modeling questions](#section-7)
9. [Section 8: Co-author's questions on the architecture paper](#section-8)
10. [Section 9: March 2026 collaborator goals and status](#section-9)
11. [Section 10: Abandoned designs and software history](#section-10)
12. [Section 11: Bibliography corrections](#section-11)
13. [Change log](#change-log)

---

## Section 0: Purpose, scope, and map {#section-0}

*Source: docs/59, Section 3.8 (survey recommendations, 2026-10-11)*
*Status: Current*

### Purpose

Appendix I collects the unique contributions of docs/00 to docs/20
that no newer whitepaper carries. It does not restate superseded
mechanisms except to name them superseded and to point to their
replacements. Each section preserves quantitative results (counts,
rates, simulation figures) from the original documents because
those numbers are part of the historical record even when the
mechanisms they rest on have been corrected.

The documents covered here were all archived to `docs/legacy/` on
2026-10-11. Each archived file carries a one-line banner at the top
naming this appendix and the series chapter that carries the
current treatment.

### Scope

Documents in scope are docs/00 to docs/20, plus the `.md` variant
of docs/09a, the `.md` variant of docs/19, and
`docs/reviews/notes.md` (the co-author comment file, formerly
`docs/notes.md`). Documents 05 (removed as a duplicate of 04)
and 10, 11 (moved to `analysis/report/`) are not archived in
`docs/legacy/` and are not covered here.

Documents dropped entirely (no content carried anywhere):

- **docs/00** (`00-documentation-index.tex`): superseded by the
  generated series index `fm0`. Nothing unique.
- **docs/01** (`01-codebase-overview.tex`): repository description
  that belongs in `CLAUDE.md` and package documentation. Nothing
  unique.
- **docs/07a** (`07a-pd-constraints-theory.tex`): superseded and
  partly wrong; root cause attributed to OL+BDC expectancy
  discontinuity was later refuted by `docs/36` Section 8.
  Nothing carried (surveys a and c agree).

### Map from original document to new home

The table below covers only docs/00 to docs/20. For docs/21 and
later see the series index (`fm0`) and docs/59 Section 4.

| Old document | Primary new home | Appendix I section |
|---|---|---|
| 00-documentation-index | archive; `fm0` (generated) | none |
| 01-codebase-overview | archive | none |
| 02 (PDF, review comments) | `docs/reviews/` | Section 8 |
| 03-audit-and-revision-report | archive; `ch03` (commit attr.) | Sections 1, 9 |
| 04-revised-power-analysis | archive | Sections 1, 3 |
| 06-ar1-residual-correlation | archive; `ch22` (formula) | Section 4 |
| 07-positive-definiteness | archive; `ch04` (counts kept) | Section 1 |
| 07a-pd-constraints-theory | archive (nothing carried) | none |
| 08-biomarker-correlation-decay | archive; `ch08` (credit) | Section 7 |
| 09-carryover-correlation-artifact | archive; `ch27` (pointer) | Section 7 |
| 09a-carryover-bias-variance | archive | Section 6 |
| 13-figure4-walkthrough | archive | Section 2 |
| 14-trial-analysis-implementation | archive; `ch25` | Section 9 |
| 15-tidyverse-alignment-plan | archive | Section 10 |
| 16-tidyverse-code-review | archive | Section 10 |
| 17-simplification-plan | archive | Sections 10, 11 |
| 18-response-parameter-sensitivity | archive | Section 3 |
| 19-mean-moderation-impl (all three) | archive; `ch07` (pointer) | Section 5 |
| 20-annotated-bibliography | archive; Appendix G | Section 11 |

---

## Section 1: Commit history and PD audit {#section-1}

*Source: docs/legacy/04-revised-power-analysis.tex, Section A.5;*
*docs/legacy/07-positive-definiteness-failures.tex, Section 3;*
*docs/legacy/23-documentation-compendium.md, Section 3.1*
*Status: Historical record. Commit attribution corrected:*
*the publication commit is `58b32a9`, not `42ac030`.*

### Note on commit attribution

The original docs/04 and docs/07 labeled the publication commit as
`42ac030`. This was corrected to `58b32a9` during the series Phase
1 audit (see `analysis/report/NOTATION.md` Part 1 and `ch03`).
The commit labels below use the corrected identifiers where known;
the tabulated figures come from the original docs/04 and docs/07
runs, which used the labels `42ac030` (publication), `8609f12`
(Ron Thomas edits) and three later commits.

### Per-commit agreement with published Figure 4

Power estimates from 200-replicate runs at each commit, compared
to published Figure 4 by Pearson correlation across 48 cells
(no-censoring rows only).

| Commit | Description | $r$ with Figure 4 | PD fail rate |
|---|---|---|---|
| `42ac030` (publication) | Publication code | **0.994** | 24.7% (40/162) |
| `8609f12` (Ron Thomas edits) | Proportional scaling added | 0.705 | 63.0% (102/162) |
| `f6ee86d` (autocorr typo fix) | Cross-factor symmetry fixed | 0.705 | 63.0% (102/162) |
| `f70f86d` (carryover added) | `Dbc` added (broken) | 0.691 | 63.0% (102/162) |
| HEAD (pre-revision) | Both bugs in `Dbc` fixed | 0.613 | 63.0% (102/162) |

Table: Agreement with published Figure 4 and PD failure rates by
commit (48 cells, no-censoring rows; docs/legacy/04 Table A.7)

The first post-publication commit (`8609f12`) is the primary source
of divergence. All subsequent commits remained in the diverged
state.

### Changes introduced at `8609f12`

Three changes were co-packaged in `8609f12`:

1. **Scale factor on carryover decay.** The decay formula changed
   from $(1/2)^{t_{sd}/t_{1/2}}$ to $(1/2)^{2 t_{sd}/t_{1/2}}$,
   halving the effective half-life. The parameter was documented
   as `@param scalefactor TODO update when understand what this
   does?`.
2. **Proportional correlation scaling.** The step-function BM-BR
   correlation was replaced by a ratio of BR means
   (`(mm1/mm0) * modelparam$c.bm`). Two side effects: timepoint 1
   always received zero correlation (an index-error guard), and
   off-drug correlations decayed proportional to the BR mean ratio.
3. **Verbose diagnostics.** No functional impact.

The power consequence: carryover power drop disappeared entirely
under this commit (N-of-1 power stayed near 0.78 regardless of
carryover). The mean-ratio scaling preserved biomarker interaction
signal during carryover periods. The PD failure rate more than
doubled from 24.7% to 63.0%.

### Positive definiteness: counts and condition numbers

All 162 matrices in the standard parameter sweep were
ill-conditioned ($\kappa > 100$) under compound symmetry at every
commit. The condition number $\kappa = \lambda_{\max} /
\lambda_{\min}$.

| Commit | Median $\kappa$ | $\kappa > 10^2$ | $\kappa > 10^3$ | $\kappa > 10^6$ |
|---|---|---|---|---|
| `42ac030` | 229 | 162/162 | 40/162 | 40/162 |
| `8609f12` | $3.1 \times 10^{18}$ | 162/162 | 102/162 | 102/162 |
| `f6ee86d` | $3.1 \times 10^{18}$ | 162/162 | 102/162 | 102/162 |
| `f70f86d` | $3.1 \times 10^{18}$ | 162/162 | 102/162 | 102/162 |

Table: Condition number statistics across 162 matrices by commit
(docs/legacy/07, Table 2)

All PD failures involved exactly one negative eigenvalue. Minimum
eigenvalues in failed matrices ranged from -33.26 (worst at
`8609f12`) to -0.52 (least severe). With a dominant eigenvalue of
approximately 700, the -33 values represent roughly 5% of the
dominant eigenvalue, indicating a substantial structural
inconsistency.

### PD failure rates by design

At the publication commit, failures were concentrated in OL+BDC:

| Design | $c_{bm} = 0$ | $c_{bm} = 0.3$ | $c_{bm} = 0.6$ |
|---|---|---|---|
| OL | 0/6 (0%) | 0/6 (0%) | 0/6 (0%) |
| CO | 0/12 (0%) | 1/12 (8%) | 1/12 (8%) |
| N-of-1 | 2/24 (8%) | 2/24 (8%) | 2/24 (8%) |
| OL+BDC | 10/12 (83%) | 11/12 (92%) | 11/12 (92%) |

Table: PD failure rates by design and $c_{bm}$ at the publication
commit; each cell shows failures out of total matrices across all
$N$ and carryover values
(docs/legacy/07, Table 3)

OL+BDC accounted for 32 of the 40 failures (80%). The OL design
never failed.

*Note (2026-10-07).* docs/36, working from the published commit
`58b32a9`, attributes the 40 failures to the requested BM-BR
correlation exceeding the CS ceiling of 0.256 to 0.262 rather than
to the OL+BDC expectancy discontinuity that docs/07 identified.
The OL+BDC account explains the design-level concentration;
docs/36's coupling-ceiling account gives the structural cause.
Both are correct and complementary: the coupling ceiling under CS
is tight enough that OL+BDC's longer design with more off-drug
timepoints is the first to push past it. Current treatment: `ch04`.

<!-- TODO: Verify the commit-mapping `42ac030` to `58b32a9` by
running `git log --oneline` in the vendored original implementation
and confirming `58b32a9` is the sixth commit dated 2020-11-19. The
corrected attribution is stated in NOTATION.md Part 1 but was not
independently re-verified for this appendix. -->

---

## Section 2: Published simulation grid {#section-2}

*Source: docs/legacy/13-figure4-walkthrough.tex, Sections 2, 4, 6*
*Status: Historical record. "Last observation carried forward"*
*corrected to monotone dropout. Dropout form checked against*
*`censordata`.*

### Four trial designs and their paths

The published simulation (Hendrickson et al., 2020) compared four
trial designs. Each design is specified by its week-by-week
on-drug indicator, blinding status, and expectancy indicator.

| Design | Weeks | Paths | On-drug structure |
|---|---|---|---|
| Open-Label (OL) | 1-8 | 1 | All on drug, unblinded |
| OL+BDC | 1-8 | 2 | Wks 1-4 OL; wks 5-8 randomized (cont./disc.) |
| Crossover (CO) | 1-8 | 2 | Drug-first (DDDD-PPPP) or placebo-first (PPPP-DDDD)|
| N-of-1 (Hybrid) | 1-12 | 4 | Wks 1-4 OL; wks 5-12 two 4-wk crossover cycles |

Table: Four trial designs in the published simulation grid

The Hybrid design has four paths arising from two binary
randomizations (drug-first or placebo-first in each cycle). The
OL design has a single path; all others have two or four paths.

### Model parameters for Figure 4

The published parameter grid for Figure 4:

```r
coremodelparams <- expand.grid(
  N = c(35, 70),
  c.bm = c(0, 0.3, 0.6),
  carryover_t1half = c(0, 0.1, 0.2),
  c.tv = 0.8, c.pb = 0.8, c.br = 0.8,
  c.cf1t = 0.2, c.cfct = 0.1
)
```

This gives $2 \times 3 \times 3 = 18$ model parameter
combinations, times 4 trial designs, times 1 response parameter
set, times 1 baseline parameter set = 72 cells. Each cell was
replicated 1,000 times, and each replicate was analyzed under
5 censoring conditions (no censoring plus 4 censoring patterns),
yielding $72 \times 5 = 360$ power estimates.

Analysis parameters: `useDE = FALSE` (expectancy not a fixed
effect), `t_random_slope = FALSE` (random intercept only),
`full_model_out = FALSE`.

### Censoring model and five patterns

The censoring model combines time-dependent (flat) and
outcome-dependent (biased) dropout:

- **Flat dropout:** probability increases linearly with interval
  duration; controlled by `beta0 * (1 - beta1)`.
- **Biased dropout:** probability higher for participants whose
  symptoms have not improved; controlled by `beta0 * beta1`.

Once a timepoint is censored, all subsequent timepoints for that
participant are set to NA (monotone dropout). The original
docs/13 described this as "last observation carried forward";
the correction to "monotone dropout" was made in this appendix
(2026-10-11).

| Pattern name | `beta0` | `beta1` | `eb1` | Description |
|---|---|---|---|---|
| No censoring | -- | -- | -- | Complete data |
| balanced | 0.05 | 0.5 | 2 | Equal flat and biased dropout |
| more of flat | 0.05 | 0.2 | 2 | Mostly time-dependent dropout |
| more of biased | 0.05 | 0.8 | 2 | Mostly outcome-dependent dropout |
| high dropout | 0.15 | 0.5 | 2 | Triple the overall dropout rate |

Table: Five censoring parameter sets used in Figure 4 (docs/legacy/13,
Section 4)

### Published analysis model

The original analysis used `lme4::lmer` with a random intercept
per participant and no explicit residual correlation structure.
For designs with within-subject drug variation (CO, N-of-1, OL+BDC),
the formula was:

```{r}
Sx ~ bm + t + Dbc + bm:Dbc + (1|ptID)
```

For OL (no within-subject drug variation), the formula was:

```{r}
Sx ~ bm + t + bm:t + (1|ptID)
```

The `lmer` analysis model is superseded; current treatment uses
`nlme::lme` with `corCAR1` (see Section 4 and `ch22`).

### Corrections to docs/13

- The "last observation carried forward" characterization of
  censoring in docs/13, Section 4 is incorrect; the code sets
  censored and all later timepoints to NA, which is monotone
  dropout without LOCF imputation.
- docs/13, Sections 1 and 3, mislabel the 2024 repository state
  as the publication code; the publication code is commit
  `58b32a9`.

---

## Section 3: Published Figures 4 and 5 reproductions {#section-3}

*Source: docs/legacy/04-revised-power-analysis.tex, Section 2;*
*docs/legacy/18-response-parameter-sensitivity.tex, Sections 1-5*
*Status: Historical record. The Figure 4 reproduction is*
*verified at r = 0.994 (200 replicates). Figure 5 material*
*(docs/18) is reported as a historical result; its figures*
*are missing.*

### Figure 4 reproduction

The analysis in docs/04 began by extracting the original
simulation code (publication commit, since corrected to
`58b32a9`) and verifying that it reproduces the published Figure
4. Using 200 replicates per cell, the reproduced power estimates
correlate with the publication at $r = 0.994$.

The primary finding of Figure 4B (the N-of-1 design's power
crashing from approximately 0.95 to 0.18 when carryover half-life
increases from 0 to 0.1 weeks) was reproduced. As documented in
`ch04` and `ch27`, this crash is an artifact of the step-function
BM-BR correlation rule, not a genuine consequence of carryover at
those half-lives.

Figure 4A (no carryover) is not affected by the step-function
artifact: at $t_{1/2} = 0$ there are no off-drug carryover
residuals. Power estimates from Figure 4A are valid and remain
the reference for design comparisons at zero carryover.

### Figure 5 (response parameter sensitivity)

docs/18 (`18-response-parameter-sensitivity.tex`) compared power
before and after revising the response parameters (Gompertz
rate, displacement, and maximum). This comparison is reported
here as a historical result under the April 2026 construct.

The reproduction status is: figures in docs/18 are missing from
`docs/legacy/figures/`; the document references figure files
that were not included in the archive. The numerical results
from docs/18 (parameter sensitivity comparisons) could be
regenerated with driver `GP-26` specified in docs/59, but this
has not been done.

<!-- TODO: Confirm whether the docs/18 figure files were ever
committed to the repository or were only produced locally. If
they can be regenerated via GP-26, decide whether the historical
comparison belongs in this appendix or in a manuscript annex. -->

---

## Section 4: Analysis-model switch and Type I error {#section-4}

*Source: docs/legacy/06-ar1-residual-correlation.tex,*
*Sections 1, 3.3, 6*
*Status: Historical record. The `lmer` Type I error evidence is*
*unique. The residual-correlation formula is superseded by*
*the compendium-mathematics.tex (`eq:inva-ri`). The claim that*
*`lmer` "always understates" correlation is qualified; it is*
*understated for the designs most affected (CO, OL) under*
*AR(1) data.*

### The problem

When the DGP was updated from compound symmetry (CS) to AR(1)
within-factor autocorrelation, the analysis model's assumption of
independent residuals (given a random intercept) became
misspecified. The original analysis used `lme4::lmer` with a
random intercept only:

```{r}
lmer(Sx ~ bm + t + Dbc + bm:Dbc + (1|ptID), data = datamerged)
```

Under CS, the random intercept absorbs the constant correlation
across all timepoints, and the residuals are approximately
independent. Under AR(1), the correlation decays with the time
gap, and the random intercept no longer captures the full
correlation structure.

### Observed Type I error inflation

With AR(1) data and the `lmer` analysis model (200 replicates,
no censoring, $c_{bm} = 0$):

| Design | $N = 35$ | $N = 70$ |
|---|---|---|
| OL | 0.13 | 0.09 |
| OL+BDC | 0.09 | 0.05 |
| CO | **0.17** | **0.13** |
| N-of-1 | 0.06 | 0.08 |

Table: Type I error rates with AR(1) DGP and `lmer` analysis
(200 replicates, no censoring; docs/legacy/06, Section 3.3)

The CO design showed a 17% false positive rate at $N = 35$
(3.4 times the nominal rate). OL was similarly inflated. OL+BDC
and N-of-1 were closer to nominal.

### Why CO and OL are most affected

For CO, the drug indicator `Dbc` is perfectly confounded with the
first versus second half of the trial. Under AR(1), the
autocorrelation decays across the trial, creating a
time-dependent correlation gradient that the random intercept does
not capture. The `bm:Dbc` interaction absorbs spurious signal
from this correlation structure.

For OL, the model uses `bm:t` as the interaction term. Under
AR(1), residual correlation decays with time, creating
time-dependent heteroscedasticity that inflates the `bm:t`
standard error downward.

OL+BDC and N-of-1 have more complex path structures with multiple
drug-status transitions. The `Dbc` variable is less confounded
with linear time, so the AR(1) structure has less impact.

### The fix: `nlme::lme` with `corCAR1`

The analysis model was changed to `nlme::lme` with continuous-time
AR(1) residual correlation:

```{r}
nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc,
          random = ~1|ptID,
          correlation = corCAR1(form = ~t|ptID),
          data = datamerged)
```

The `corCAR1` structure models
$\mathrm{Cor}(\varepsilon_{i,t_1}, \varepsilon_{i,t_2}) =
\phi^{|t_1 - t_2|}$, where $\phi$ is estimated from the data.
This directly matches the AR(1) structure used in the DGP.

After the switch, Type I error was corrected to 2.5-5.8% across
all designs at 1,000 replicates (verified against the stored
simulation run; see `ch22`).

*Note.* The original claim in docs/06 that `lmer` "always
understates" residual correlation under AR(1) data should be
read as "understates for the affected designs"; OL+BDC at
$N = 35$ (0.09) and N-of-1 at both sizes are near nominal, so
the inflation is design-specific, not universal.

---

## Section 5: Mean-moderation formulations that failed {#section-5}

*Source: docs/legacy/19-mean-moderation-implementation-notes.md*
*(reviewed; code claims re-checked against current*
*`R/generateData.R` on 2026-10-07);*
*docs/legacy/19-mean-moderation-implementation-notes.tex,*
*Sections 6-7*
*Status: Historical record. Three failed attempts documented.*
*The positive definiteness ceiling figures ("about 0.45") in*
*the `.tex` version are superseded by `ch09` values. The power*
*table (50 replicates) is superseded by paper 01 and by*
*`docs/37`.*

### Background

The original R package (`R/generateData.R`) implements covariance
moderation: the BM-BR interaction is encoded as a differential
correlation in the covariance matrix of the MVN draw. A parallel
tidyverse codebase (pre-alignment, 2025) used additive mean
moderation applied post-MVN draw. Three formulations were attempted
before the correct one was identified.

### Attempt 1: Multiplicative scaling (incorrect)

```r
dat[, (br_col) := get(br_col) * (1 + beta_bm * bm_centered)]
```

When the base BR value is near zero (off-drug timepoints with
minimal carryover), multiplying by `(1 + small)` produces
negligible moderation signal. The analysis model sees strong
biomarker signal on-drug but almost none off-drug, mimicking
covariance moderation's correlation erosion. Power dropped about
60% under carryover, no improvement over MVN.

### Attempt 2: Additive with drug exposure decay (incorrect)

```r
drug_exposure <- (1/2)^(tsd / t1half)
dat[, (br_col) := get(br_col) + beta_bm * bm_centered *
  drug_exposure * br_sd]
```

Drug exposure decays rapidly during off-drug periods (0.18 at the
first off-drug timepoint, 0.001 by the fourth, at 2.5-week spacing
with $t_{1/2} = 1$ week). The additive signal is proportional to
this decay, so the off-drug moderation is still negligible relative
to the noise ($\sigma_{BR} = 8$). The `bm:Dbc` interaction term
sees minimal contrast. Power remained low under carryover.

### Attempt 3: Additive with time-on-drug, uncentered (incorrect scaling)

```r
dat[, (br_col) := get(br_col) + bm * tod * beta_bm]
```

With uncentered biomarker (mean approximately 124) and $t_{od}$
scaling with weeks (up to 10-20), the total shift at the mean
biomarker was approximately $0.45 \times 124 \times 10 = 558$
units at $t_{od} = 10$, against noise SD of 8. Power was 1.00
everywhere. This was an artifact of a grossly overpowered signal.

The correct accounting (verified, 2026-10-08): writing
$B_i = \mu_{bm} + \sigma_{bm} b_i$, the total shift is
$\beta_{bm}\mu_{bm}t_{od} + \beta_{bm}\sigma_{bm}t_{od}\, b_i$.
At $\beta_{bm} = 0.45$, $\mu_{bm} = 124.33$, $\sigma_{bm} = 15.36$
and $t_{od} = 10$, the modulation that the interaction term sees is
$0.45 \times 15.36 \times 10 = 69.1$ BR units per SD of the
biomarker, approximately 8.6 $\sigma_{BR}$ and 19.2 times the
calibrated 3.6. The drift term (559.5 units, the drug and time
terms absorb it) carries no interaction information. The
conclusion -- a grossly overpowered signal -- stands.

### Attempt 4: Centered, on-drug binary, calibrated (correct)

```r
bm_z <- (dat$bm - bm_mean) / bm_sd
for (tp in 1:nP) {
  if (d[tp]$onDrug) {
    dat[, (br_col) := get(br_col) + c.bm * bm_z * br_sd]
  }
}
```

In mathematical form, the shift at every on-drug visit is
$\beta_{bm}\, b\, \sigma_{BR}$, with $\beta_{bm}$ passed as
`c.bm`. The standardized biomarker $b$ (code `bm_z`) centers the
signal (population-level mean moderation = 0). The $\sigma_{BR}$
scaling produces comparable effect sizes: a 1-SD biomarker shift
produces $\beta_{bm}\sigma_{BR}$ units of BR shift.

### Calibration check ($c_{bm} = \beta_{bm} = 0.45$, extracted parameters)

At the shared value $c_{bm} = \beta_{bm} = 0.45$ and
$\sigma_{BR} = 8$ from `extracted_rp`:

- Mean moderation: 3.36 BR units per SD of biomarker
- Covariance moderation: 3.48 BR units per SD of biomarker
- Theoretical: $0.45 \times 8 = 3.6$ BR units per SD

Both empirical shifts fall within sampling error of the
theoretical 3.6 (the `.tex` version gives the 500-participant
crossover draw). The range "0.65-0.9" given for signal-to-noise in
the original docs/19 was not reproduced; the ratio of per-SD shift
to $\sigma_{BR}$ is 0.45 (flagged 2026-10-07).

### Pilot power comparison (50 replicates, superseded)

The 50-replicate pilot from docs/19.md:

| $t_{1/2}$ | Mean moderation | Covariance moderation |
|:---:|:---:|:---:|
| 0.0 | 0.74 | 0.82 |
| 0.5 | 0.68 | 0.64 |
| 1.0 | 0.72 | 0.50 |

Table: N-of-1 power by $t_{1/2}$ for both architectures, $N = 70$,
$c_{bm} = 0.45$, 50 replicates
(docs/legacy/19.md, power table; MCSE approximately 0.065)

*Superseded (2026-10-07).* Paper 01's re-run at matched total
$N = 70$ gives Hybrid carryover losses of 6.9% (mean moderation)
and 15.4% (covariance moderation), as recorded in
`analysis/report/NOTATION.md` Part 2. The 50-replicate figures
have Monte Carlo errors too large to interpret precisely. The
ordering discrepancy between the 0.5-week and 1.0-week cells in
docs/19.md versus docs/21 (which swaps them) is unresolved;
either ordering is within Monte Carlo error of flat power under
mean moderation.

*Note (docs/37).* The closed-form result that under on-drug-only
mean moderation power does not depend on the half-life explains
why the 0.5 and 1.0 cells are interchangeable within noise.

---

## Section 6: Constant-exposure study and crossover bias {#section-6}

*Source: docs/legacy/09a-carryover-bias-variance.md,*
*Sections 2-7 (corrected 2026-10-08)*
*Status: Historical record. Several central claims were not*
*supported by the paper's own results (see below). The*
*Shi and Ye (2024) expression is given with corrected signs.*
*Conclusion 4 (the bias-variance non-monotonic argument) is*
*withdrawn. The `.tex` variant of docs/09a contradicts the*
*corrected `.md` and is not used here.*

### Study description

docs/09a examined power behavior in a simplified open-label then
discontinuation design (single path: all participants on drug in
weeks 1-8, off drug from week 9, measurements at weeks 1-12, 16
and 20). The DGP used mean moderation with linear growth, not the
Hendrickson et al. (2020) Gompertz model. Analyses were run with
and without an explicit carryover model; 100 Monte Carlo replicates
per cell.

*Correction (2026-10-07).* The original paper described this as
a blinded discontinuation design with randomized continuation or
discontinuation. Inspection of the archived script
`simulation_carryover_spectrum.R` showed it is a single-path
design with no randomization at week 9.

### Non-monotonic power: claimed but unverified

The paper's executive summary claimed OL+BDC power without
carryover modeling showed a non-monotonic pattern: a sharp drop
from zero carryover ($\sim$85%) to short carryover ($t_{1/2}$ =
0.1-0.5 weeks, $\sim$55-60%), then partial recovery to $\sim$68%
at 2 weeks.

*Correction (2026-10-07).* The paper's own simulation (Section
5.1) does not reproduce this pattern. Without carryover modeling,
power was 0.06 at $t_{1/2} = 0$, 0.00-0.02 at 0.05-0.5 weeks,
and rose to 0.36 at 2 weeks. The claimed $\sim$85% baseline
differs from the simulated 6%. The source of the approximate
values in the paper's introductory table was not recorded; a
search of the archived scripts found none. The values are
unverified. The evidence in Section 5.1 supports a power that
increases with half-life beyond about 0.2 weeks rather than a
drop followed by recovery.

*Further correction.* The "without carryover" analysis was
misspecified even at $t_{1/2} = 0$: the predictor
`weeks_on_drug` is constant after discontinuation, which does not
correctly capture zero carryover. The result at $t_{1/2} = 0$
(power = 0.06, wrong sign on the mean estimate) reflects this
misspecification.

### Crossover carryover bias (Shi and Ye 2024)

Shi and Ye (2024) formalize carryover in a potential outcomes
framework. Let $\tau_1$ and $\tilde\tau_2$ be treatment effects in
the two periods; $\omega_0$ the carryover from control and
$\omega_1$ the carryover from active treatment. The basic crossover
estimator $\hat\tau$ has expectation:

$$E(\hat\tau) = \tfrac{1}{2}(\tau_1 + \tilde\tau_2 - \omega_1 + \omega_0)$$

The estimator is biased by $-(\omega_1 - \omega_0)/2$. With no
carryover from control ($\omega_0 = 0$), the bias is $-\omega_1/2$,
toward zero when the carryover has the same sign as the treatment
effect.

*Correction (2026-10-08, derived and checked numerically).* The
original docs/09a wrote the expectation as
$\tfrac12(\tau_1 + \tilde\tau_2 - \omega_0 - \omega_1)$, with
carryover terms summed. Under the additive two-period model, the
carryover from active treatment enters the period-2 control mean
of the AB sequence and the carryover from control enters the
period-2 active mean of the BA sequence; the two terms enter with
opposite signs. The corrected expression is stated above. Both
forms agree when $\omega_0 = 0$. See the compendium-mathematics,
subsection "Carryover bias of the classical two-period crossover"
(`eq:inva-crossover`).

### Withdrawn conclusion

**Conclusion 4 of docs/09a is withdrawn.** The conclusion that
the non-monotonic power pattern reflects a bias-variance tradeoff
was based on an unsourced table of approximate values and on a
theoretical derivation (Section 3.5) that does not produce an
interior maximum as originally stated. The precipitous decline in
the Hendrickson et al. (2020) Figure 4B has since been traced to
the step rule for the BM-BR coupling under covariance moderation
(`docs/36` Sections 5.3-5.5; `docs/37`), not to a bias-variance
mechanism.

---

## Section 7: Open modeling questions {#section-7}

*Source: docs/legacy/08-biomarker-correlation-decay.tex, Section 8;*
*docs/legacy/09-carryover-correlation-artifact.tex, Section 5*
*Status: Open items; carried to `ch08` as open items.*

### Open question 1: Trait versus state variance off drug

*Source: docs/legacy/08, Section 8*

The DGP sets the BR standard deviation $\sigma_{BR} = 8$ as
constant across on-drug and off-drug timepoints. This conflates two
sources of variation:

1. **Trait variability.** Between-person differences in latent drug
   responsiveness. These exist before the drug is administered and
   persist regardless of drug status. This component should remain
   constant.
2. **State-dependent variability.** Variation in the acute drug
   effect, which should scale with the magnitude of the drug effect
   and decay during carryover periods.

The current framework cannot separate these two sources. The
proportional correlation scaling (graded coupling,
$\lambda_{cor}$) addresses the biomarker-specific component
correctly (its conditional spread decays with the mean), but the
total BR variance remains at $\sigma_{BR}^2$ during off-drug
periods. The conditional distribution then has the correct center
but excessive width during carryover, which may conservatively
underestimate power (more noise during off-drug periods than is
clinically realistic).

Addressing this limitation would require a time-varying
$\sigma_{BR,t}$ or a decomposition of the BR variance into trait
and state components, either of which represents a substantial
restructuring of the DGP.

*Current status (ch08).* This question is carried to `ch08` as an
open item. For the current application (comparing relative power
across trial designs), the fixed-variance assumption is acceptable
because it affects all designs similarly.

### Open question 2: Clinical plausibility of negative BR under the step rule

*Source: docs/legacy/09, Section 5*

The original step-function BM-BR correlation rule (commit
`42ac030`) implied that, when the adjusted BR mean is negligible
off drug, participants with low biomarkers had negative BR values
(the drug was worsening their symptoms). The conditional
expectation of BR given BM is:

$$E[BR_t \mid BM] = \mu_{BR,t} + c_{bm} \cdot
\frac{\sigma_{BR}}{\sigma_{BM}} \cdot (BM - \mu_{BM})$$

When $\mu_{BR,t} \approx 0$ (off drug, full carryover decay) and
$c_{bm} = 0.6$, a participant with $BM = \mu_{BM} - 2\sigma_{BM}$
has expected $BR_t \approx -9.6$, implying the drug worsened
their symptoms by approximately one symptom SD. This was described
as clinically incoherent in docs/09: no rebound effect was
hypothesized in the original publication, and the BR sign
convention assigns positive values to symptom improvement.

The clinical plausibility comparison:

| Criterion | Step rule | Graded coupling |
|---|---|---|
| Drug effect near zero $\Rightarrow$ BM predictive? | Yes (full) | No ($\approx 0$) |
| Negative BR during off-drug periods? | Yes (for low BM) | No |
| Proportional decay of predictive power? | No (step function) | Yes (exponential) |
| Clinically plausible? | No | Approximately |

Table: Clinical plausibility of the two correlation rules
(docs/legacy/09, Table 4)

*Current status.* The step-rule clinical incoherence is now stated
in `ch27` (the standalone note on the step rule in the published
code). The clinical-plausibility distinction is also recorded in
`ch08` as part of the derivation of the graded coupling rule.

---

## Section 8: Co-author's questions on the architecture paper {#section-8}

*Source: docs/reviews/notes.md (recorded 2026-05-06; page-keyed*
*comments on the April 2026 draft of paper 01); the annotated*
*PDF is at docs/reviews/review-comments-2026-04-29.pdf*
*(formerly docs/02-dgp-mean-moderation-vs-mvn-rchcoments.pdf).*
*Status: Partial. Ten questions identified from the text version;*
*the annotated PDF is not text-readable.*

### Questions from docs/reviews/notes.md

The following questions were recorded in `docs/reviews/notes.md`.
Page numbers correspond to the April 2026 draft.

1. **Page 2 (Abstract).** "Add future directions of
   generalization?" The abstract closes on the two-architecture
   contrast; the dual-channel architecture (paper 11) and
   analysis-model extensions are the natural follow-ons.
   *Answer:* The current paper 01 abstract is restricted to the
   mean-versus-covariance contrast per the 2026-09-30 decision
   (see `CLAUDE.md`, "Open work items"). Generalizations are
   directed to paper 11.

2. **Page 2 (Introduction).** "Is this actually the division we
   want to highlight front and center?" Query about the framing of
   the two architectures as the paper's central contrast.
   *Answer:* Yes; the architectures are the organizing contrast of
   paper 01. The matched-power result and the carryover sensitivity
   difference are both direct consequences of the distinction.

3. **Page 3 (Introduction).** "Would skip this part of framing"
   about the note that the paper-specific design choice is not a
   research-program convention.
   *Answer:* The framing note was retained in the draft to prevent
   a reader from over-generalizing the binary choice; addressed in
   revision as a footnote rather than inline.

4. **Page 3 (Introduction).** "Would seem more straightforward to
   simply explain that the impact of carryover on power was greater
   than expected in this implementation vs other similar projects,
   an exploration of why identified the dynamic explored here."
   (Comment on the paragraph about the parallel tidyverse
   reimplementation.)
   *Answer:* The current introduction follows this suggestion
   (revised after receipt of comments); the technical audit
   framing was replaced by the observational framing the
   co-author recommended.

5. **Page 3 (Introduction).** "Question of whether we 1) implement
   and include in comparisons the hybrid model, vs 2) identify the
   hybrid as future area of work, vs 3) only note in passing."
   (Comment on the dual-channel architecture.)
   *Answer:* Option 3 applies to paper 01 (note in passing was
   also removed after the 2026-09-30 extraction; paper 01 now
   contains no reference to the dual-channel architecture).
   The hybrid model is paper 11.

6. **Page 3 (Section 2 heading).** "Baseline biomarker value,
   I assume?" (Query about "biomarker".)
   *Answer:* Yes; the canonical notation is $B$ for the raw
   baseline biomarker (code `bm`), defined in `NOTATION.md`.

7. **Page 5 (Section 2.2).** "Is this just the default because we
   wrote it that way? Is it actually structurally required by the
   rest of the setup, or could this decision be dissociated?"
   (Query about the exponential-decay form as default under
   covariance moderation.)
   *Answer:* The exponential decay is a design choice (the refined
   process), not a structural requirement of the MVN setup. The
   step rule is the publication default. The decay form is
   motivated by pharmacokinetic plausibility and by its
   derivation from the coefficient of variation argument in
   `ch08`. It could be dissociated by using a different coupling
   function. *Chapter reference:* `ch08`.

8. **Page 5 (Section 2.2).** "The way this is written I don't
   understand how it's specific to this model -- isn't this what
   carryover of a drug effect always is?" (Comment on the sentence
   "Carryover effects in the DGP inflate the BR means during
   off-drug periods.")
   *Answer:* The sentence is not architecture-specific; it is a
   property of any model with exponential carryover decay. The
   specificity lies in how the BM-BR coupling interacts with the
   elevated BR mean during the off-drug period. Under covariance
   moderation, that coupling is also present during carryover;
   under mean moderation, it is not. *Chapter reference:* `ch07`.

9. **Page 7 (Section 3.2).** "'The noise structure is independent
   of drug state' -- this makes the most sense to me as a clear
   concrete difference -- how much is lost if I think of this as
   the central difference?" (Comment on the architectural
   divergence.)
   *Answer:* This formulation captures the practical consequence
   faithfully. The revised paper 01 uses this framing to lead the
   divergence section. The covariance formulation is the
   origin-of-the-difference statement; noise-independence is the
   consequence. *Chapter reference:* `ch07`, Section on divergence.

10. **Page 13 (Section 6.4).** "Should we be considering this for
    our data?" (Comment on two-stage random slopes.)
    *Answer:* This is a design decision for the analysis of the
    actual PTSD trial data, not for the simulation. The random
    slope is listed as a pending sensitivity in `docs/57` Section
    8 and `ch23`. The interaction test's operating characteristics
    under random slopes are not yet simulated. *Chapter reference:*
    `ch23` (open item).

<!-- TODO: Inspect the annotated PDF (docs/reviews/
review-comments-2026-04-29.pdf) for any comments not captured in
notes.md. The PDF is not text-readable directly; use a PDF reader
or OCR. The page-keyed comments in notes.md may be a
representative but not exhaustive record of the co-author's
annotations. -->

---

## Section 9: March 2026 collaborator goals and status {#section-9}

*Source: docs/legacy/14-trial-analysis-implementation-guide.tex;*
*docs/legacy/03-audit-and-revision-report.tex, Section 6*
*Status: Historical record. Formulas in docs/14 were superseded*
*by the 2026-10-08 code fix. Goal 2b (residual analysis) was*
*never implemented. Links to `ch23` and `ch25`.*

### The three goals

The March 2026 collaboration established three goals for the
pmsimstats-ng project in the context of preparing for analysis of
the actual prazosin-PTSD trial data.

**Goal 1: Did prazosin improve PTSD symptoms?**
Status as of docs/03: complete.

Implemented as `analyze_trial_extended()` in
`R/carryover_analysis.R`. Fits the LME model to actual trial data
and extracts both the main drug effect (`Dbc` coefficient) and the
biomarker-treatment interaction (`bm:Dbc`), with 95% confidence
intervals, variance components, predicted response curves per
participant, and a clinical decision threshold (the blood pressure
value at which predicted benefit exceeds a specified CAPS
improvement). The threshold formula was defective as of 2026-03;
it was corrected in the 2026-10-08 code fix (see `ch25`).

**Goal 2a: Characterize carryover in actual data.**
Status as of docs/03: complete.

Implemented as `characterize_carryover()` in
`R/carryover_analysis.R`. Sweeps candidate carryover half-lives,
fits `nlme::lme` at each, compares AIC/BIC and log-likelihood, and
identifies the best-fitting half-life. The function was validated
on simulated data: correctly recovers the true $t_{1/2} = 1.0$
from data generated with that half-life (verified, docs/14
Section 2 Action block).

Default half-life grid in the current code: `c(0.1, 0.25, 0.5, 1,
2, 4, 8)` weeks (inspected, `R/carryover_analysis.R` lines 27-29;
verified 2026-10-11, see `docs/series/notes/phase1/phase1-cr-log.md`
CR-43). The docs/03 formulation "0.25 to 8 weeks" was wrong for
the default; the correct default starts at 0.1 weeks.

*Open item (CR-43).* `print_carryover_summary()` reports the
drug effect and interaction p-values at the AIC-selected half-life
as if the half-life were fixed. Testing at the AIC-selected
half-life without calibrating the $p$-value is liberal (see
`ch23`, `docs/57` Section 8). The issue was noted in the Phase 1
conflict log.

**Goal 2b: Identify best sequence of response measures.**
Status: never implemented.

The conceptual design (per-timepoint residual analysis to identify
where the carryover signal is strongest) was specified in docs/14
Section 2b. The existing `full_model_out = TRUE` argument provides
the necessary infrastructure (fitted values, residuals, random
effects per participant). This analysis has not been coded. It
is an open item for `ch25`.

**Goal 3a: Simulation power calculations.**
Status as of docs/03: complete.

The revised simulation produced 1,000-replicate power estimates
with Type I error control of 2.5-5.8% across all designs, zero PD
failures across all 162 parameter combinations, sigma matrix
caching, and parallel execution via `furrr`. Runtime: 8 minutes
for 200 replicates, 47 minutes for 1,000 replicates on the
simulation host.

**Goal 3b: Actual analysis output structure.**
Status as of docs/03: complete.

Implemented in `analyze_trial_extended()`. Returns the full model
object, predicted treatment effect as a continuous function of
biomarker value, per-participant predicted response curves, and
clinical decision thresholds. The function can be applied to the
bundled `CTdata.rda` dataset from the Raskind et al. trial.

### Threshold and predicted-effect formulas

docs/14 stated formulas for the clinical decision threshold and
the predicted effect as a function of biomarker value. These
formulas were defective as of the document's date and were
corrected in the 2026-10-08 code fix. The corrected formulas
are in `ch25`.

---

## Section 10: Abandoned designs and software history {#section-10}

*Source: docs/legacy/17-simplification-plan.tex, Sections 3.3*
*and 8; docs/legacy/15-tidyverse-alignment-plan.tex;*
*docs/legacy/16-tidyverse-code-review.tex, item 1*
*Status: Historical record.*

### The simplification plan (not adopted, 2026-03)

docs/17 (`17-simplification-plan.tex`, 2026-03-24) proposed
reducing the pmsimstats2025 codebase to a minimal viable
simulation of approximately 200 lines. The proposed simplifications
included:

- Replace the Gompertz response model with a linear growth model.
- Replace `nlme::lme` with ANCOVA on phase-mean change scores
  (Method E9 in docs/17, Section 3.3.3; labeled E9 in the series).
- Replace compound symmetry with a random intercept only.
- Reduce the censoring model to no-dropout or MCAR with a constant
  probability.
- Fix the biomarker moderation strength at a single value and vary
  it only in sensitivity analyses.

The simplification plan was not adopted for the main simulation.
Its significance for the series is twofold:

1. **Origin of Method E9.** The ANCOVA-on-phase-means method
   (docs/17, Section 3.3.3) is the origin of the analysis
   identifier E9 used in the power comparison of `docs/37`,
   `docs/39`, and papers 04-06. The method is documented in
   `ch17` as the E9 analysis.

2. **Origin of the RM-ANOVA comparator.** docs/17,
   Section 3.3, proposed a repeated-measures ANOVA as an
   alternative analysis. This is the origin of the RM-ANOVA
   comparator that appears in Appendix H of the series.

The main simulation retained the Gompertz model, `nlme::lme` with
`corCAR1`, and the full censoring model; the simplified analyses
were added as comparators rather than replacements.

### The 2025 tidyverse alignment (docs/15 and docs/16, 2026-03)

docs/15 (`15-tidyverse-alignment-plan.tex`) planned aligning the
exploratory 2025 tidyverse code (`analysis/2025/pm_functions.R`)
to match the data.table original (`R/`) function by function. The
alignment was executed in the 2026-04-08 merge that introduced the
`implementations/` directory. The 2025 code was moved to an
external archive at
`~/Dropbox/prj/alz/10-pmsimstats-ng/archive/2025/`; the maintained
tidyverse implementation now lives in
`implementations/tidyverse/R/functions.R`.

docs/16 (`16-tidyverse-code-review.tex`) reviewed the pre-merge
`analysis/2025/pm_functions.R`. The first finding (item 1):
`calculate_carryover()` (lines 24-77) was dead code -- a
multi-model carryover function supporting exponential, linear, and
Weibull decay, never called by any function in the file. The actual
carryover was computed inline using only the exponential formula.
The linear and Weibull models were never integrated into the DGP.
This finding informed the rewrite: `implementations/tidyverse/`
does not include the dead function.

---

## Section 11: Bibliography corrections {#section-11}

*Source: docs/legacy/20-annotated-bibliography.tex;*
*docs/legacy/17-simplification-plan.tex, Section 7;*
*docs/legacy/23-documentation-compendium.md, Section 17;*
*docs/series/notes/phase1/phase1-cr-log.md, GP-29*
*Status: Corrections verified or dropped as of Phase 1 (2026-10-11).*
*Each entry noted below was checked against Crossref or PubMed*
*or is flagged as unverified.*

### Verified entries in docs/20

The following entries in docs/20 were verified to be
bibliographically correct (title, authors, journal, year, DOI):

- Thomas et al. (2020) Frontiers in Digital Health 2:13.
  doi:10.3389/fdgth.2020.00013. (Verified; the paper.)
- Blackston et al. (2019) Healthcare 7(4):137.
  doi:10.3390/healthcare7040137. (Verified against Crossref.)
- Senn (2019) Stat Meth Med Res 28(2):372-383.
  doi:10.1177/0962280217726801. (Verified against Crossref.)
- Pequignat et al. (2023) Behavior Research Methods 55:3824-3842.
  doi:10.3758/s13428-022-02012-1. (Verified against Crossref.)
- Senn, Julious, and Chen (2024) Trials 25:264.
  doi:10.1186/s13063-024-08100-1. (Verified against Crossref.)
- De Bock, Van den Noortgate, and Gielen (2021) Behavior Research
  Methods 53:1782-1798. doi:10.3758/s13428-020-01522-0.
  (Verified against Crossref.)
- Raskind et al. (2018) NEJM 378(6):507-517.
  doi:10.1056/NEJMoa1507598. (Verified against Crossref.)

### Corrections from Phase 1 (GP-29)

**Chen, Li, and Yang (2020).** docs/20 cited this as:

> Chen Y, Li P, Yang K. Comparison of aggregated N-of-1 trials
> with parallel and crossover randomized controlled trials using
> simulation studies. PLoS One. 2020;15(1):e0216949.
> doi:10.1371/journal.pone.0216949.

The title given by docs/20 matches the title of Blackston et al.
(2019). The Chen, Li, and Yang (2020) paper at doi:10.1371/
journal.pone.0216949 has a different title (PubMed: "Comparison of
statistical analysis strategies of aggregated N-of-1 trials" --
unverified; the full text at that DOI was not checked against the
annotation). The annotation in docs/20 describes a replication and
extension of the Blackston simulation framework, which is
consistent with what the Chen et al. paper covers. The citation
itself (authors, year, journal, DOI) appears correct; the title
in the bibliography is likely transcribed from Blackston et al.
The correct title should be confirmed against the full text.

**Cruz, Mylona, and Melo (2024).** docs/20 gave an incorrect title
for this entry (noted and corrected in docs/20, inline, 2026-10-09):
> Penalized GEE for complex carry-over in repeated-measures
> crossover designs. arXiv:2402.16362. 2024.
> Published as Statistics in Medicine 2026;45(10-12):e70561,
> doi:10.1002/sim.70561.

The title correction is incorporated in the version archived at
`docs/legacy/20-annotated-bibliography.tex`. Appendix G should use
the corrected title.

**Eng and Konings (2018) Trials 19:72.** docs/59 listed this as
a suspect entry (GP-29). The full citation is: Eng J, Konings MK.
A framework for comparing the effectiveness of radiologic
screening strategies. Trials. 2018;19:72.
doi:10.1186/s13063-018-2434-9. This paper is about radiologic
screening, not N-of-1 trials or carryover. Its inclusion in the
annotated bibliography appears to be an error. It is dropped from
Appendix G.

### Entries not in docs/20 but referenced in docs/17

docs/17, Section 7, cited Grizzle (1965) and Freeman (1989) in the
context of the two-stage procedure for carryover testing. These
citations were not in docs/20. They are carried to Appendix G with
the caveat that they were cited from memory in docs/17 and have not
been verified against the full texts.

- Grizzle JE. The two-period change-over design and its use in
  clinical trials. Biometrics. 1965;21(2):467-480.
- Freeman PR. The performance of the two-stage analysis of
  two-treatment, two-period crossover trials. Statistics in
  Medicine. 1989;8(12):1421-1432.

The accuracy of the Freeman (1989) claim in docs/09a (that it
showed Grizzle's two-stage procedure inflates Type I error and
judged it too misleading for practical use) was not verified
against the full text.

---

## Change log

| Date | Section | Change | Reason |
|---|---|---|---|
| 2026-10-15 | All | First draft written from docs/04, 06, 07, 08, 09, 09a, 13, 14, 15, 16, 17, 19, 20, 23 | Phase 2 content move |
| 2026-10-08 | Section 6 | Crossover bias formula corrected (signs) | Shi and Ye (2024) derivation |
| 2026-10-08 | Section 5 | Attempt 3 miscalculation corrected | Direct computation |
| 2026-10-07 | Section 6 | Withdrawn conclusion 4 noted | Simulation evidence does not support |
| 2026-10-07 | Section 5 | Power table superseded note added | Paper 01 re-run |
