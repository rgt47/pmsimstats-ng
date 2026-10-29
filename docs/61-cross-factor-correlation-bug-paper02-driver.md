---
geometry: margin=1in
header-includes:
  - \AtBeginDocument{\fvset{fontsize=\footnotesize}}
---

# The Cross-Factor Correlation Mismatch in Paper 02's Simulation Driver: `c.cfct` Set to Half of `c.cf1t` Where the Refined Process Specifies One Constant {.unlisted .unnumbered}
*2026-10-21 10:23 PDT*

**Author.** pmsimstats team

**Purpose.** This note describes one specific defect found on
2026-10-21 in `analysis/scripts/carryover-sensitivity/simulation-core.R`
and one other driver script
(`analysis/scripts/carryover-sensitivity/26-co-dbc-fix-paired.R`), the
data-generating drivers behind paper 02 ("A factorial comparison of
nine carryover representations for aggregated N-of-1 trials under
covariance moderation"). The refined process specifies the
correlation between different response components, at any pair of
timepoints including the same visit, as a single separable form,
$c_1\,\rho^{\lvert w_t - w_s\rvert}$, with $c_1 = 0.2$ (paper 01,
Table `tab:spec`; `docs/series/ch10-refined-process.md`). The
driver's two model parameters that jointly implement this formula,
one for same-timepoint pairs and one for different-timepoint pairs,
were set to two different values, so every cross-factor correlation
beyond the same visit was exactly half the documented value, for the
entire production history of paper 02 until this session's fix. This
note states the intended model, quotes the code as it stood, derives
the exact size of the discrepancy, works a numeric example with
verified values, states what is and is not affected, gives the fix
verbatim, and records what remains open. It is written for a reader
with statistical training (correlation matrices, variance components)
but no prior exposure to this project's code, in the style of
`docs/53` and `docs/55`, which document analogous defects in
Hendrickson et al.'s published code; this defect is instead in this
project's own code, not in that external reference. Section 3 uses the
Kronecker-product (separability) framing developed at length in
`docs/63-why-correlation-fades-with-time-a-clinicians-guide.md`; a
reader unfamiliar with that framing may want to read Section 3 of
that note first, though this note restates what is needed to follow
the bug itself.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows the project's canonical table
(`analysis/report/NOTATION.md`) where applicable. Symbols marked
*local* are defined here for this note only.

| Symbol | Meaning | Status |
|------------|--------------------------------------------------|----------|
| $\rho$ | within-factor AR(1) autocorrelation per week (code `c.tv`, `c.pb`, `c.br`, here always 0.7 in the affected cells) | canonical |
| $w_t$, $w_s$ | calendar week of visit $t$, $s$ | canonical |
| $c_1$ | the refined process's single cross-factor correlation constant (paper 01 Table `tab:spec`) | canonical |
| `c.cf1t` | code parameter: cross-factor correlation at the same timepoint (lag 0) | local (code) |
| `c.cfct` | code parameter: cross-factor correlation scale at different timepoints, multiplying $\rho^{\text{lag}}$ | local (code) |
| TV, PB, BR | the three response components: natural course, placebo belief, biological drug response | canonical |
| $c_{bm}$ | on-drug biomarker-BR correlation (code `c.bm`); unaffected by this defect | canonical |
| "the affected driver" | `simulation-core.R`'s `simulate_cell()` and `simulate_cell_s7()`, and `26-co-dbc-fix-paired.R`'s `run_cell()` | local |
| "the shared function" | `build_correlation_matrix()` in `implementations/tidyverse/R/functions.R`, mirrored by `buildSigma()` in the package's `R/generateData.R` | local |

Table: Notation index for this note

**Glossary.**

- **Cross-factor correlation.** The correlation between two
  *different* response components (e.g. TV and PB) at a given pair
  of visits, as distinct from the within-factor AR(1) correlation
  between the same component at two visits.
- **The refined process.** The pmsimstats-ng data-generating process
  recommended by paper 01 as the reference DGP for simulation studies
  of this kind (`NOTATION.md`; paper 01 Table `tab:spec`).
- **The affected driver.** The caller-side code that builds the
  `model_param` list passed into the shared function, not the shared
  function itself (Section 4 explains why this distinction matters
  here).

**Evidence labels.** Each claim carries one of four labels: *quoted*
(text reproduced from a source), *inspected* (read from source code),
*derived* (follows algebraically) or *verified* (computed by running
code, with the exact commands given in Section 9).

## 1. Summary for the busy reader

**What the code did.** `simulation-core.R`'s two model-parameter
constructors (`simulate_cell()`, line 467 as now fixed, and
`simulate_cell_s7()`, line 878 as now fixed) and
`26-co-dbc-fix-paired.R`'s `run_cell()` (line 130) set
`c.cf1t = 0.2` but `c.cfct = 0.1` (inspected, prior to this session's
fix). The shared function that consumes these,
`build_correlation_matrix()`, uses `c.cf1t` directly for same-
timepoint cross-factor pairs and `c.cfct * rho^lag` for every other
pair (inspected, `implementations/tidyverse/R/functions.R` lines
236-258).

**What it should have done.** The refined process specifies a single
constant $c_1 = 0.2$ scaling $\rho^{\lvert w_t-w_s\rvert}$ at *every*
lag, including lag 0 where $\rho^0 = 1$ makes the formula collapse to
exactly $c_1$ (quoted, paper 01 Table `tab:spec`: "Cross-component
correlation \textbar{} separable, $c_1\rho^{\lvert w_t-w_s\rvert}$
[...] \textbar{} $c_1 = 0.2$"; `docs/series/ch10-refined-process.md`
line 115 states the same form for the refined process against the
published process's constant $c_\times = 0.1$ at every lag). The
code's two-parameter structure is general enough to represent this
single-constant form as a special case — it simply requires
`c.cf1t` and `c.cfct` to be set to the same value. They were not.

**How large.** Because the intended and the actual value both scale
by the same $\rho^{\text{lag}}$ factor, the ratio of actual to
intended is exactly `c.cfct` / $c_1$ = 0.1 / 0.2 = **0.5** at every
lag greater than zero, independent of $\rho$ and of the specific lag
(derived). At the same timepoint (lag 0), the correlation was
correct throughout, since that case uses `c.cf1t` directly, which was
set to 0.2 as documented. Concretely, at a 4-week gap with
$\rho = 0.7$: the intended TV-PB correlation is
$0.2 \times 0.7^4 = 0.04802$; the code produced $0.1 \times 0.7^4 =
0.02401$ (verified, Section 6).

**What this is not.** This is not a defect in the shared function
`build_correlation_matrix()` (or the package's `buildSigma()`, which
implements the identical two-parameter structure). Both functions
correctly implement a general separable cross-factor form; the
defect is entirely in what the caller passed for `c.cfct`. No change
to the shared function was needed or made.

**What it does and does not touch.** The within-factor AR(1)
correlations (`c.tv`, `c.pb`, `c.br`) and the biomarker-BR
differential coupling `c.bm` are generated by separate code blocks
and are unaffected (inspected). What is affected is the background
correlation among the three nuisance/signal components (TV, PB, BR)
with each other across different visits — every cross-timepoint,
cross-factor entry of the covariance matrix, in every design with
more than one visit, which is every design this manuscript uses
(Hybrid, OL+BDC, CO all have multiple visits). At the one reference
cell checked (Hybrid, $c_{bm} = 0.45$, $t_{1/2} = 1.0$), the
resulting covariance matrix remained positive definite under both the
old and the corrected parameterization, though the minimum eigenvalue
was smaller under the old, weaker-correlation parameterization
(0.1025) than under the corrected one (0.1834) at this cell — a
verified but not further explained fact (verified, Section 6).

**The one-line-per-site fix.** Change `c.cfct` from `0.1` to `0.2` at
each of the three sites (Section 8), so both constants carry $c_1$
and the formula collapses to the documented single-constant separable
form at every lag.

**What remains open.** Whether this defect changed any of paper 02's
reported power, Type I error, or ranking conclusions has not been
assessed. Every stored simulation output and every table or figure in
`report.Rmd` and `supplement.Rmd` generated before this session's fix
was computed under the halved cross-factor correlation. A rerun is a
separate, larger task, tracked as Part 2 step 4 of
`docs/60-trials-submission-gap-analysis.md`, and has not been
performed as of this writing.

## 2. Scope

This note covers only the cross-factor correlation parameterization
defect. Two other defects were found and fixed in the same driver in
the same session (deterministic rather than randomized path
allocation; silent rather than optional-stop non-positive-definite
matrix repair), and are documented separately in
`docs/60-trials-submission-gap-analysis.md` and in `CLAUDE.md`'s
"Resolved 2026-10-21" entries, not here. This note also does not
re-derive the refined process itself, which is paper 01's subject;
it takes paper 01's Table `tab:spec` and `docs/series/ch10-refined-
process.md` as the specification against which the code is checked.

The defect affects every call into `generate_data_multi_path()` from
`simulate_cell()` or `simulate_cell_s7()`, which is to say every
cell of paper 02's Tier 1 grid (`G1`-`G9` at the factorial grid of
trial design, $N$, $c_{bm}$, $t_{1/2}$), every Tier 2 sensitivity
block (S1-S11), the decay-shape sensitivity runs, and the
higher-precision 24-cell paired rerun, from paper 02's first
production run through 2026-10-21. It does not affect paper 01's own
reported results, which do not use this driver (paper 01's own
reproducibility pipeline is separate; see paper 01 Appendix B and
its Reproducibility section).

## 3. The intended model

Paper 01's Table `tab:spec` (quoted):

> Cross-component correlation \textbar{} separable,
> $c_1\rho^{\lvert w_t-w_s\rvert}$; the Kronecker form applies to the
> correlation, since the PB SD scales with the expectancy \textbar{}
> $c_1 = 0.2$

`docs/series/ch10-refined-process.md` states the same contrast
against the published process (quoted, line 115):

> Cross-factor correlation ($t \ne s$) \textbar{} constant
> $c_\times = 0.1$ \textbar{} separable: $c_1 \rho^{\lvert t-s
> \rvert}$

and gives the separable form's defining property explicitly (line
215, quoted): "Separable AR(1), $\rho = 0.7$ per week,
$c_1 = 0.2$ \textbar{} choice: raises coupling ceiling from 0.256 to
0.480 (Hybrid)". The form is separable in the sense that the
cross-factor correlation at any pair of visits, same or different,
is one constant $c_1$ times the same decay $\rho^{\lvert
w_t-w_s\rvert}$ that governs the within-factor AR(1) structure; at
$w_t = w_s$ this is $c_1 \times 1 = c_1$ exactly, so a single
parameter $c_1$ fully determines the cross-factor correlation at
every timepoint pair.

**In matrix terms** (see
`docs/63-why-correlation-fades-with-time-a-clinicians-guide.md`,
Section 3, for the full derivation): the correlation structure over
every (component, visit) pair is meant to be the Kronecker product
$\Sigma = K_p \otimes K_t$, with $K_t$ the AR(1) matrix
($[K_t]_{t,s}=\rho^{\lvert w_t-w_s\rvert}$) and $K_p$ the
$3\times 3$ cross-component correlation matrix (TV, PB, BR), with 1's
on the diagonal and the single constant $c_1$ off it. Entrywise,
$\Sigma_{(c,t),(c',s)} = [K_p]_{c,c'}\cdot[K_t]_{t,s}$, which for
$c\ne c'$ is exactly $c_1\,\rho^{\lvert w_t-w_s\rvert}$. The code's
two parameters, `c.cf1t` and `c.cfct`, are not two independent
structural choices; they are two different code paths computing the
*same* entry of $K_p$ (the off-diagonal $c_1$) at two different
ranges of lag ($\ell=0$ and $\ell>0$ respectively). Separability
requires both code paths to return the same number. Section 4 shows
they did not.

## 4. The code as it stood

`build_correlation_matrix()`, `implementations/tidyverse/R/functions.R`,
lines 236-258 (inspected; unchanged by this fix, quoted here to show
the two-parameter structure the caller must supply consistently):

```r
for (factor2_idx in setdiff(seq_along(factors), factor_idx)) {
  other_factor <- factors[factor2_idx]

  n1_same <- paste(trial_design$timepoint_name,
                   current_factor, sep = ".")
  n2_same <- paste(trial_design$timepoint_name,
                   other_factor, sep = ".")
  idx1_same <- cbind(match(n1_same, labels), match(n2_same, labels))
  idx2_same <- cbind(idx1_same[, 2], idx1_same[, 1])
  correlations[idx1_same] <- model_param$c.cf1t
  correlations[idx2_same] <- model_param$c.cf1t

  if (num_timepoints > 1) {
    n1_lag <- paste(trial_design$timepoint_name[upper_pairs$p1],
                    current_factor, sep = ".")
    n2_lag <- paste(trial_design$timepoint_name[upper_pairs$p2],
                    other_factor, sep = ".")
    lag_vals <- model_param$c.cfct * ar_vals
    idx1_lag <- cbind(match(n1_lag, labels), match(n2_lag, labels))
    idx2_lag <- cbind(idx1_lag[, 2], idx1_lag[, 1])
    correlations[idx1_lag] <- lag_vals
    correlations[idx2_lag] <- lag_vals
  }
}
```

Here `ar_vals <- autocorrelation^time_gaps` with `autocorrelation`
the within-factor $\rho$ computed earlier in the same loop iteration
(line 222), so `lag_vals` is exactly `c.cfct * rho^lag` as stated in
Section 1. The package's `buildSigma()` (`R/generateData.R`, lines
322-344) implements the identical structure: `correlations[n1,n2] <-
modelparam$c.cf1t` for same-timepoint pairs (lines 333-334) and
`cfval <- modelparam$c.cfct * rho^tg` for different-timepoint pairs
(line 342). Neither function was changed by this fix.

What was changed is the caller. `simulate_cell()`
(`simulation-core.R`, before this session's fix, inspected from the
version control history):

```r
model_param <- list(
  N = cell$N,
  carryover_t1half = cell$t1half,
  carryover_form = cell$carryover_form,
  weibull_shape = cell$weibull_shape,
  dgp_architecture = cell$dgp_arch,
  c.tv = rho, c.pb = rho, c.br = rho,
  c.cf1t = 0.2, c.cfct = 0.1
)
```

`simulate_cell_s7()` and `26-co-dbc-fix-paired.R`'s `run_cell()` set
the identical pair, `c.cf1t = 0.2, c.cfct = 0.1` (inspected).

## 5. Why the discrepancy is exactly a factor of 2

Write the intended cross-factor correlation at lag $\ell =
\lvert w_t - w_s \rvert$ as $c_1 \rho^\ell$, and the code's actual
value as `c.cfct` $\times\, \rho^\ell$. The ratio is

$$
\frac{\text{actual}}{\text{intended}} =
\frac{\texttt{c.cfct} \times \rho^\ell}{c_1 \times \rho^\ell} =
\frac{\texttt{c.cfct}}{c_1} = \frac{0.1}{0.2} = 0.5,
$$

independent of $\rho$ and of $\ell$, for every $\ell > 0$ (derived).
At $\ell = 0$ the same-timepoint branch of the code uses `c.cf1t`
directly, which was set to $0.2 = c_1$ correctly, so same-visit
cross-factor correlations were never affected. The defect is
therefore not a diffuse numerical error but a clean, uniform halving
of every cross-timepoint, cross-factor correlation entry, for as long
as the driver has existed.

## 6. Worked example: Hybrid design, path 1

The Hybrid design's eight visits fall at weeks 4, 8, 9, 10, 11, 12,
16 and 20 (named `OL1, OL2, BD1, BD2, BD3, BD4, COd, COp` in the
`tidyverse` implementation's trial design; verified).

**Same timepoint (lag 0).** TV-PB correlation at `OL1`: 0.2 under
both the old and the corrected parameterization (verified) — unaffected,
as expected.

**Different timepoints (lag 4 weeks, `OL1` to `OL2`, weeks 4 to 8).**
At $\rho = 0.7$: intended $0.2 \times 0.7^4 = 0.04802$; code (before
the fix) produced $0.1 \times 0.7^4 = 0.02401$; ratio 0.5 exactly
(verified, matching Section 5's derivation).

**General pattern across lags**, $\rho = 0.7$ (derived from the
verified ratio and the known decay):

| Lag (weeks) | $\rho^{\text{lag}}$ | Intended ($c_1 = 0.2$) | Code (before fix, `c.cfct = 0.1`) |
|---|---|---|---|
| 0 | 1.0000 | 0.2000 | 0.2000 (unaffected; uses `c.cf1t`) |
| 1 | 0.7000 | 0.1400 | 0.0700 |
| 4 | 0.2401 | 0.0480 | 0.0240 |
| 8 | 0.0576 | 0.0115 | 0.0058 |
| 12 | 0.0138 | 0.0028 | 0.0014 |

Table: Intended versus actual cross-factor correlation by lag,
$\rho = 0.7$; the lag-4 row is directly verified, the others derived
from the same exact-factor-of-2 relationship

**Matrix-level effect.** At the reference cell (Hybrid, $N = 70$,
$c_{bm} = 0.45$, $t_{1/2} = 1.0$ week), the Frobenius norm of the
difference between the old and corrected correlation matrices is
0.625 (verified). Both the old and the corrected covariance matrices
are positive definite at this cell; the minimum eigenvalue is 0.1025
under the old (weaker cross-correlation) parameterization and 0.1834
under the corrected one (verified). This is reported as an observed
fact at one cell, not a general claim about how cross-factor
correlation strength relates to positive-definiteness across the
grid; no further mechanism is derived here.

## 7. Consequences

**What changes.** Every cross-timepoint, cross-factor correlation
entry (TV-PB, TV-BR, PB-BR, at every pair of distinct visits) in
every design the manuscript uses. In an eight-visit design this is
28 of the $\binom{8}{2} = 28$ possible cross-timepoint pairs per
factor pair, all of them halved from their documented value; only the
8 same-timepoint pairs per factor pair were correct. With three
factor pairs (TV-PB, TV-BR, PB-BR), the overwhelming majority of the
model's cross-factor correlation structure was affected.

**What does not change.** The within-factor AR(1) correlations
(`c.tv`, `c.pb`, `c.br`, all $\rho = 0.7$ in the affected cells), the
biomarker-BR differential coupling (`c.bm`, the core
interaction-carrying parameter under covariance moderation), graded
coupling's off-drug decay, carryover on the drug-response mean, and
the analysis model itself are all generated or specified by code
paths independent of `c.cf1t`/`c.cfct` (inspected). The defect is
confined to the background correlation structure among the
nuisance/signal components with each other, not to the biomarker-
treatment interaction the manuscript's power comparisons target
directly.

**Unknown downstream impact.** Whether a halved cross-factor
correlation among TV, PB and BR measurably changes any of paper 02's
reported power, Type I error, or specification rankings is not
established by this note. The covariance misspecification this
introduces relative to the documented refined process could in
principle affect the residual structure the `corCAR1` analysis model
assumes, and therefore standard errors and power, in either
direction; no simulation comparing old and corrected parameterizations
head to head on the actual test statistics has been run as of this
writing. This is the natural next step and is tracked separately
(`docs/60-trials-submission-gap-analysis.md`, Part 2, step 4).

**Scope of affected output.** Every `.rds` output file and every
table or figure in `analysis/report/02-carryover-sensitivity/report.Rmd`
and `supplement.Rmd` produced by a script that calls `simulate_cell()`
or `simulate_cell_s7()` was generated under the halved cross-factor
correlation, from paper 02's first production run through 2026-10-21.

## 8. The fix

Three sites, each changing `c.cfct` from `0.1` to `0.2`:

**`simulation-core.R`, `simulate_cell()`** (model_param list,
previously ending `c.cf1t = 0.2, c.cfct = 0.1`):

```r
c.tv = rho, c.pb = rho, c.br = rho,
c.cf1t = 0.2, c.cfct = 0.2
```

**`simulation-core.R`, `simulate_cell_s7()`** (same change, same
pattern).

**`26-co-dbc-fix-paired.R`, `run_cell()`** (same change):

```r
c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.2
```

No change was made to `build_correlation_matrix()` or `buildSigma()`;
both already correctly implement the general two-constant structure
and require only that the caller supply matching constants to recover
the documented single-$c_1$ form.

**Verification that the fix restores the documented ratio of 1** (not
0.5) at lag $> 0$: see Section 9's reproduction commands, which
confirm the lag-4, $\rho = 0.7$ cell now equals $0.2 \times 0.7^4 =
0.04802$ exactly, and that the lag-0 value is unchanged at 0.2 under
both the old and new parameterizations.

Several older one-off diagnostic scripts in the same directory
(`diagnose-null-type1*.R`, `check-null-architecture-equivalence.R`,
`06-cr2-recalibration.R`, `diagnose-alt-cell-cr2.R`) and the earlier
`quick-sim-prototype*.R` exploratory scripts still construct
`model_param` directly with the old, or in some cases a different
(`c.cf1t = 0.1, c.cfct = 0.05`), pair of constants. None of these were
found to feed the manuscript's reported tables or figures directly,
and they were left unchanged in this pass; this was not exhaustively
verified for every script in the directory, so a later sweep should
confirm none of their output is cited before relying on that claim.
Correction (2026-10-29): the `diagnose-null-type1*.R`,
`diagnose-alt-cell-cr2.R` and `check-null-architecture-equivalence.R`
scripts do feed paper 10
(`analysis/report/10-interaction-test-calibration/`), whose tables
were produced with `c.cfct = 0.1`, and they now take `c.cfct = 0.2`
from `paper10-common.R`.

## 9. How to verify independently

Run in a checkout of this repository, from the repository root:

```r
library(corpcor)
source("implementations/tidyverse/R/functions.R")
source("analysis/scripts/carryover-sensitivity/simulation-core.R")

design_set <- build_design_set("Hybrid")
td <- design_set[[1]]
rp <- default_resp_param(); bp <- default_baseline_param()

mp_new <- list(N = 70, carryover_t1half = 1.0,
               carryover_form = "exponential", weibull_shape = 1,
               dgp_architecture = "mvn", c.bm = 0.45,
               c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
               c.cf1t = 0.2, c.cfct = 0.2)       # corrected
mp_old <- mp_new; mp_old$c.cfct <- 0.1           # as the driver stood before this fix

s_new <- build_sigma_matrix(mp_new, rp, bp, td, dgp_architecture = "mvn")
s_old <- build_sigma_matrix(mp_old, rp, bp, td, dgp_architecture = "mvn")

tdat <- prepare_trial_data(td)
nm <- td$timepoint_name

# lag 0: unaffected
i1 <- paste0(nm[1], ".tv"); i2 <- paste0(nm[1], ".pb")
s_old$correlations[i1, i2]; s_new$correlations[i1, i2]   # both 0.2

# lag 4 weeks (OL1 to OL2): halved before the fix
j1 <- paste0(nm[1], ".tv"); j2 <- paste0(nm[2], ".pb")
s_old$correlations[j1, j2]   # 0.02401
s_new$correlations[j1, j2]   # 0.04802

norm(s_old$correlations - s_new$correlations, type = "F")  # 0.625
is.positive.definite(s_old$sigma); is.positive.definite(s_new$sigma)  # both TRUE
min(eigen(s_old$sigma, only.values = TRUE)$values)   # 0.1025
min(eigen(s_new$sigma, only.values = TRUE)$values)   # 0.1834
```

All values above reproduce exactly (verified, this session).

## 10. Limitations

- **No downstream power or Type I error impact has been assessed.**
  This note documents the DGP-level defect and its exact size, not
  its effect on any of paper 02's specific reported results. See
  Section 7.
- **Positive-definiteness and eigenvalue behavior were checked at one
  reference cell only** (Hybrid, $c_{bm} = 0.45$, $t_{1/2} = 1.0$),
  not swept across the full grid of designs, sample sizes, effect
  sizes and half-lives paper 02 uses.
- **Older diagnostic and prototype scripts were not exhaustively
  audited**; Section 8 names the ones found still using stale
  constants and states only that they were not found to feed the
  manuscript directly, not that this was verified for every script.
- **This note does not address the driver's other two defects**
  (deterministic path allocation; silent non-PD repair), which were
  found and fixed in the same session but are out of scope here; see
  `docs/60-trials-submission-gap-analysis.md` and `CLAUDE.md`.

## 11. Reproducibility

Section 9 gives the exact commands used to produce every verified
number in this note. No stored data files are required; the
comparison is deterministic given the model parameters shown (no
random seed is involved, since `build_sigma_matrix()` only
constructs the covariance matrix, not simulated draws from it).

## 12. References

- Paper 01, Table `tab:spec`
  (`analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd`).
- `docs/series/ch10-refined-process.md`, the refined process's
  canonical specification.
- `docs/60-trials-submission-gap-analysis.md`, Part 2, which first
  identified this defect (step 1) as part of a broader session
  auditing paper 02 against the refined process.
- `CLAUDE.md` (pmsimstats-ng project file), "Resolved 2026-10-21"
  entries, for this defect's sibling fixes (path allocation, non-PD
  handling) in the same driver.
- `docs/53-carryover-recursion-in-published-code.md` and
  `docs/55-step-coupling-in-published-code.md`, whose structure this
  note follows, though both describe defects in Hendrickson et al.'s
  published code rather than in this project's own.
