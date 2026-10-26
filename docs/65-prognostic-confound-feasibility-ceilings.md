# Prognostic Confounds and the Coupling Ceiling: A Silent Repair Artifact and Its Analytic Boundary

*2026-10-26 16:56 PDT*

**Author.** pmsimstats team

**Purpose.** Two paper 02 simulation jobs (jobs 28 and 30) reported
that every analysis specification, `G`-family and `F`-family alike,
lost 15 to 25 points of power when a strong prognostic confound was
added to the data-generating process alongside the true
biomarker-treatment interaction. This note shows that the loss was
not a property of any analysis model but of the data-generating
process itself: the confound strength used was infeasible, the joint
covariance matrix was not positive definite, and the drivers' silent
repair step removed part of the interaction signal before any model
saw the data. It then derives, in closed form, the exact boundary of
the feasible region for a prognostic coupling alongside the
predictive coupling $c_{bm}$, so that future stress cells can be
chosen inside it. Both halves are verified by code in the
repository. The studies are exploratory.

## 1. Summary

- At confound strength 0.6 on the prog-PB and prog-TV-const pathways
  (and prog-BL at $c_{bm} = 0.45$), the joint correlation matrix of
  the refined-process DGP is non-positive-definite on every Hybrid
  path, with smallest eigenvalue between $-0.11$ and $-0.85$.
- The paper 02 drivers ran with silent repair enabled, so
  `corpcor::make.positive.definite` replaced the matrix by a nearby
  positive-definite one. The repair shrinks the population
  on-drug/off-drug slope contrast of the outcome on the biomarker,
  which is the quantity the interaction term estimates, from
  $-3.5$ to $-2.4$ (prog-PB) and from $-3.3$ to $-1.8$
  (prog-TV-const). The fitted estimates from the simulations match
  these repaired values, and standard errors do not grow.
- No analysis model can recover a signal the DGP did not generate.
  Job 30's three random-effects augmentations of `G9` changed power
  by at most one point in every cell, which is the expected result
  once the mechanism is understood.
- The compendium's Schur-complement ceiling for $c_{bm}$ extends
  directly to a biomarker row with several couplings. The feasible
  region is the interior of an ellipsoid $c^{\top} A c < 1$, with
  $A$ computed from the inverse of the correlation matrix with the
  biomarker removed. At $c_{bm} = 0.45$ on Hybrid, the largest
  feasible prog-PB or prog-TV-const coupling is $0.306$ and the
  largest prog-BL coupling is $0.496$; at $c_{bm} = 0.30$ they are
  $0.451$ and $0.815$. A coupling of 0.6 on PB or TV is infeasible
  at every $c_{bm}$.
- Consequences: job 28's strength-0.6 power numbers must not be
  cited as a confound effect; every valid cell (strength 0.3) shows
  no power cost at all; the project's non-PD default is now stop,
  not repair, in every collection; and stress grids should be
  chosen from the ceiling table in Section 5 and checked before
  submission.

## 2. Background

### 2.1 The prognostic pathways

Paper 01's Study 3 and `docs/57` Section 7.3 define four ways a
biomarker can be prognostic without being predictive, each a
separate channel in the DGP:

| Pathway | Parameter | Mechanism |
|---|---|---|
| prog-BL | `c.bm.bl` ($c_{bl}$) | correlation between the biomarker and baseline severity BL |
| prog-TV-const | `c.bm.tv` ($c_{tv}$) | correlation between the biomarker and the natural-history component TV at every timepoint |
| prog-TV-grow | `b.bm.tv` ($b_{tv}$) | an additive shift of the TV mean proportional to the standardized biomarker |
| prog-PB | `c.bm.pb` ($c_{pb}$) | correlation between the biomarker and the placebo-response component PB at every timepoint |

Three of the four are entries in the biomarker row of the joint
correlation matrix. prog-TV-grow is a mean shift and never touches
the covariance; it is absorbed exactly by the biomarker main effect
and plays no part in what follows.

### 2.2 The predictive coupling and its ceiling

Under covariance moderation (the refined process), the
biomarker-treatment interaction is encoded as a correlation $c_{bm}$
between the biomarker and the drug-response component BR, equal to
$c_{bm}$ on drug and decaying off drug. `docs/34` and the mathematics
compendium show that $c_{bm}$ is bounded above by a ceiling
$c_{bm}^{\ast}$ beyond which the joint matrix is not positive
definite, and give the ceiling in closed form by a Schur complement
on the biomarker row. Every production grid in papers 01 and 02 sits
below it.

### 2.3 Silent repair

Hendrickson's published code repairs a non-positive-definite
covariance matrix with `corpcor::make.positive.definite(tol = 1e-3)`
and continues. The refined process (paper 01, Table `tab:spec` and
Appendix B.2) recommends stopping instead. A stop option
(`stopOnNonPD` in the package, `stop_on_nonpd` in the tidyverse
collection) was added on 2026-10-21 with a backward-compatible
default of `FALSE`, so every caller that did not opt in kept
repairing silently. The paper 02 drivers did not opt in.

## 3. The issue: what jobs 28 and 30 actually measured

### 3.1 Job 28, the prognostic stress test

Driver `28-run-g9-prognostic-stress.R` crossed the four pathways with
confound strength $\{0.3, 0.6\}$, $c_{bm} \in \{0, 0.45\}$ and
$N \in \{35, 70\}$ on the Hybrid design, 500 replicates per cell,
fitting `G1`-`G9`, `E0` and `F1`/`F1A`/`F1A0`/`F1A0B` to the same
simulated trial each replicate. Its two questions were whether a
doubled confound inflates `G9`'s Type I error and whether a confound
co-occurring with a true interaction costs power.

Type I error was nominal everywhere: `G9` rejected between 0.042 and
0.050 at all sixteen null cells. Power at $c_{bm} = 0.45$ (`G9`,
$N = 35$ / $N = 70$):

| Pathway | No confound | Strength 0.3 | Strength 0.6 |
|---|---|---|---|
| prog-BL | 0.759 / 0.969 | 0.782 / 0.976 | 0.738 / 0.950 |
| prog-PB | 0.759 / 0.969 | 0.798 / 0.968 | 0.597 / 0.853 |
| prog-TV-const | 0.759 / 0.969 | 0.782 / 0.968 | 0.502 / 0.772 |
| prog-TV-grow | 0.759 / 0.969 | 0.778 / 0.970 | 0.778 / 0.970 |

(The no-confound baseline is job 30's, 1000 replicates.) At strength
0.3 every pathway is at or above baseline within Monte Carlo error.
At 0.6 the PB and TV-const pathways drop by 17 to 26 points, and so
does every other specification in the job: F1A0B, for instance, falls
from 0.741 to 0.562 on prog-PB at $N = 35$. The uniformity across
fourteen specifications of three different families was the first
hint that the loss was upstream of the analysis.

### 3.2 Job 30, a hypothesis that was wrong

The initial reading of job 28 was that the confound added
patient-level variance that the models' random-intercept structure
could not absorb. Driver `29-run-g9-random-effects.R` therefore fit
`G9` with three richer random-effects structures (random slope on
time, on expectancy, and both) on the same twelve cells at 1000
replicates. None recovered anything: all three variants matched
plain `G9` to within one point of power and of standard error in
every cell, with paired per-replicate differences of at most 0.008,
at the cost of up to 2.5% non-convergence.

The per-cell estimates in the same output already contradicted the
variance hypothesis. The mean `G9` estimate of the interaction at
$N = 35$ was $-3.61$ with no confound, $-2.50$ under prog-PB 0.6 and
$-1.97$ under prog-TV-const 0.6, while the empirical standard error
fell from 1.38 to 1.16 to 1.02. An added-variance mechanism would
leave the estimate centered and inflate its spread; what was
observed was attenuation with a tighter spread.

### 3.3 The direct check

Script `30-check-prognostic-confound-pd.R` builds the covariance
matrix for every (pathway, strength, $c_{bm}$) combination the jobs
used, on each of the four Hybrid paths, with no simulation, and
reports two things: whether the raw matrix is positive definite
(`build_sigma_matrix()`'s own `was_pd_corrected` flag and the
smallest eigenvalue of $DRD$ with $D$ the standard deviations), and
the population slope of the outcome on the biomarker at each
timepoint, computed from the matrix directly. Since the outcome is
$S_t = \mathrm{BL} - TV_t - PB_t - BR_t$,

$$
\beta_t = \frac{\mathrm{Cov}(bm, \mathrm{BL}) - \mathrm{Cov}(bm, TV_t)
 - \mathrm{Cov}(bm, PB_t) - \mathrm{Cov}(bm, BR_t)}{\mathrm{Var}(bm)},
$$

and the on-drug minus off-drug average of $\beta_t$ is the population
value of the interaction contrast. Comparing it before and after
repair isolates what the repair did to the signal.

| Pathway | Strength | $c_{bm}$ | Paths repaired | Min. eigenvalue | Contrast, raw | Contrast, repaired |
|---|---|---|---|---|---|---|
| none | 0 | 0.45 | 0 | 0.180 | $-3.26$ | $-3.26$ |
| prog-BL | 0.3 | 0.45 | 0 | 0.145 | $-3.26$ | $-3.26$ |
| prog-BL | 0.6 | 0.45 | 4 | $-0.110$ | $-3.26$ | $-3.00$ |
| prog-PB | 0.3 | 0.45 | 0 | 0.010 | $-3.40$ | $-3.40$ |
| prog-PB | 0.6 | 0 | 4 | $-0.240$ | $-0.27$ | $-0.24$ |
| prog-PB | 0.6 | 0.45 | 4 | $-0.642$ | $-3.53$ | $-2.38$ |
| prog-TV-const | 0.3 | 0.45 | 0 | 0.011 | $-3.26$ | $-3.26$ |
| prog-TV-const | 0.6 | 0 | 4 | $-0.354$ | 0.00 | 0.00 |
| prog-TV-const | 0.6 | 0.45 | 4 | $-0.853$ | $-3.26$ | $-1.81$ |
| prog-TV-grow | 0.6 | 0.45 | 0 | 0.180 | $-3.26$ | $-3.26$ |

Three things follow. First, every strength-0.3 cell is positive
definite, which is why no power was lost there. Second, the two
strength-0.3 cells on PB and TV-const have a smallest eigenvalue of
0.010, almost exactly on the boundary; Section 5 shows why. Third,
the repaired contrasts, $-2.38$ and $-1.81$, are what the fitted
models reported ($-2.50$ / $-2.32$ and $-1.97$ / $-1.84$ at
$N = 35$ / $70$). The power loss is the repair removing 33% to 45% of
the interaction from the generated data. There was nothing for any
analysis model to recover.

### 3.4 Two side observations

**prog-PB is not a pure confound of the interaction.** Even at
$c_{bm} = 0$, the raw contrast under prog-PB is nonzero ($-0.14$ at
strength 0.3, $-0.27$ at 0.6, against $-3.3$ at the alternative).
The placebo-response amplitude scales with expectancy, which in the
Hybrid design differs between on-drug and off-drug visits, so a
biomarker correlated with PB acquires a small on/off slope
difference by construction. `G9` showed no detectable Type I error
inflation from it at $N \le 70$, but this is exactly the mechanism
the `F`-family's expectancy term exists to absorb, and it may matter
at larger $N$ or stronger couplings.

**Zero-variance components are a separate case.** Flipping the
repair default (Section 6) exposed a second, unrelated use of the
repair. Because PB's amplitude scales with expectancy, a timepoint
with expectancy 0 (a run-in phase) has exactly zero PB variance and
makes the covariance matrix singular by construction. That is a
degenerate component, not an infeasible parameter set. The stop path
now checks and factors only the positive-variance sub-block and draws
the degenerate components at their means; no correlation is modified.

## 4. The analytic boundary

### 4.1 Setup

For one participant on one path with $n$ visits, the DGP draws
$(B, \mathrm{BL}, TV_{1:n}, PB_{1:n}, BR_{1:n})$ from a multivariate
normal distribution. Write $R$ for its correlation matrix and $R_0$
for $R$ with the biomarker row and column removed. $R_0$ contains
the within-component AR(1) blocks, the cross-component correlations
and the uncorrelated BL, and depends on no biomarker coupling. The
biomarker row is a sum of fixed patterns scaled by the couplings:

$$
r = c_{bm}\, g + c_{pb}\, p + c_{tv}\, q + c_{bl}\, b,
$$

where $g$ is the predictive pattern (1 on the BR positions while on
drug, $2^{-t_{sd}/t_{1/2}}$ off drug after first exposure, 0 before
exposure, 0 elsewhere), $p$ is 1 on every PB position, $q$ is 1 on
every TV position and $b$ is 1 on BL. The variances rescale $R$ to
the covariance and do not affect positive definiteness.

### 4.2 The condition

By the Schur complement on the biomarker row, $R \succ 0$ if and
only if $R_0 \succ 0$ and

$$
r^{\top} R_0^{-1} r < 1 .
$$

The left side is the biomarker's squared multiple correlation on all
other variables, so the condition says only that the biomarker cannot
be more than perfectly predictable from the rest. Substituting the
pattern expansion,

$$
r^{\top} R_0^{-1} r = c^{\top} A\, c, \qquad
A_{jk} = x_j^{\top} R_0^{-1} x_k, \quad x \in \{g, p, q, b\},
$$

with $c = (c_{bm}, c_{pb}, c_{tv}, c_{bl})$. The feasible region is
the interior of the ellipsoid $c^{\top} A c < 1$ in the positive
orthant. The compendium's single-coupling ceiling is the special case
$c_{bm}^{\ast} = A_{gg}^{-1/2}$.

### 4.3 The boundary for one coupling given the others

Fix all couplings but one, say $c_{pb}$. Then $c^{\top} A c = 1$ is a
quadratic in $c_{pb}$,

$$
A_{pp}\, c_{pb}^2 + 2\Bigl(\sum_{k \ne p} A_{pk} c_k\Bigr) c_{pb}
 + \Bigl(\sum_{j,k \ne p} c_j A_{jk} c_k - 1\Bigr) = 0,
$$

and the ceiling is its positive root. If the discriminant is
negative the fixed couplings are already infeasible on their own. A
design's ceiling is the minimum over its paths, since all paths share
one parameter set.

Because BL is uncorrelated with every response component, $b$ is an
eigenvector of $R_0$ with $A_{bb} = 1$ and $A_{bk} = 0$ for
$k \ne b$, so $c_{bl}$ enters only as $c_{bl}^2$: the prog-BL ceiling
given $c_{bm}$ is $\sqrt{1 - c_{bm}^2 A_{gg}}$. The PB and TV
patterns are structurally identical (a constant on every position of
one component with the same AR(1) block and the same cross-component
correlations), so $A_{pp} = A_{qq}$ and $A_{gp} = A_{gq}$, and
prog-PB and prog-TV-const have identical ceilings. Expectancy
scaling of PB lives in the standard deviations and does not enter.

### 4.4 Implementation and verification

Script `31-prognostic-coupling-ceiling.R` obtains $R_0$ and the four
patterns from the paper 02 driver's own `build_sigma_matrix()`
(building the correlation matrix at zero couplings, then at each
unit coupling in turn), forms $A$ per path, and solves the quadratic.
Each closed-form ceiling is checked against the smallest eigenvalue
of the full $R$: positive at $0.999$ times the ceiling, negative at
$1.001$ times it. All fourteen cases in the table below pass.

## 5. Ceilings for the paper 02 configuration

Hybrid design, minimum over its four paths, $t_{1/2} = 0.5$ week,
$\rho = 0.7$, $c_1 = 0.2$ at every lag (the corrected constants of
`docs/61`).

| Coupling | Alone | Given $c_{bm} = 0.30$ | Given $c_{bm} = 0.45$ |
|---|---|---|---|
| $c_{bm}$ | 0.518 | | |
| $c_{pb}$ | 0.513 | 0.451 | 0.306 |
| $c_{tv}$ | 0.513 | 0.451 | 0.306 |
| $c_{bl}$ | 1.000 | 0.815 | 0.496 |

In the other direction, given $c_{pb}$ or $c_{tv} = 0.3$ the ceiling
on $c_{bm}$ is 0.454; given 0.6, no $c_{bm}$ is feasible.

The coefficient matrix for path 1 (the others differ in the third
decimal) is

$$
A = \begin{pmatrix}
 3.725 & -0.400 & -0.400 & 0 \\
-0.400 &  3.797 & -0.633 & 0 \\
-0.400 & -0.633 &  3.797 & 0 \\
 0 & 0 & 0 & 1
\end{pmatrix}
$$

in the order $(g, p, q, b)$. The negative off-diagonals mean the
couplings are mildly complementary rather than purely competing: the
response components are positively correlated with one another, so
a biomarker tied to BR already implies a little correlation with PB
and TV, and an explicit $c_{pb}$ of the same sign uses up less of the
budget than it would if the components were independent.

Reading job 28 against this table: at $c_{bm} = 0.45$, strength 0.3
on PB or TV-const sat at 0.306, inside the region by 0.006, which is
the eigenvalue of 0.010 in Section 3.3; strength 0.6 was outside it
at every $c_{bm}$; and prog-BL at 0.6 was outside it at
$c_{bm} = 0.45$ (ceiling 0.496) but inside at $c_{bm} = 0$. The
$c_{bm}$-alone ceiling of 0.518 is higher than the 0.480 Hybrid
figure in `docs/34` because that figure used the pre-`docs/61`
cross-lag constant of 0.1 and $t_{1/2} = 0$.

## 6. Consequences and recommendations

1. **Do not cite job 28's strength-0.6 power numbers as a confound
   effect.** The valid finding from jobs 28 and 30 is narrower and
   cleaner: `G9` holds nominal size under every prognostic pathway at
   both strengths, and at every feasible strength tested it loses no
   power relative to the no-confound alternative. The job 30
   random-effects variants had no artifact-free motivation and are
   not a direction worth pursuing.
2. **Non-positive-definite matrices are now refused, never repaired,
   by default** (commit `4aa4d0ef4`): in the package, the tidyverse
   and original-extended collections, the nof1power engine, and
   paper 02's driver core, where a refused cell appears as a
   non-converged replicate. Silent repair survives only as an
   explicit `FALSE` opt-in for the scripts that study the repair
   itself and for reproducing runs made before 2026-10-26.
   `implementations/original` is untouched, since silent repair is
   part of the published process it reproduces.
3. **The stop is a backstop, not a design tool.** Parameters must be
   chosen inside the feasible region before a run, from the table in
   Section 5 or by running script 31 for a new configuration. A
   pre-flight eigenvalue check over the whole grid, as the package's
   `generateSimulatedResults()` already performs, should be added to
   the paper 02 drivers so an infeasible cell is rejected before the
   simulation host spends time on it.
4. **A valid "stronger confound" test is still open.** At
   $c_{bm} = 0.45$ there is almost no room above 0.3. The question
   can be asked at $c_{bm} = 0.30$, where confounds up to 0.45 are
   feasible on PB and TV-const and up to 0.81 on BL, or at
   $c_{bm} = 0.45$ with prog-BL up to 0.49. Whether a ladder of this
   kind is worth running depends on whether the manuscript needs to
   claim robustness beyond strength 0.3; it currently does not.
5. **The prog-PB expectancy mechanism (Section 3.4) deserves a
   targeted size check** at larger $N$ before `G9`'s robustness to a
   prognostic placebo responder is stated without qualification.

## 7. Reproducibility

- `analysis/scripts/carryover-sensitivity/30-check-prognostic-confound-pd.R`
  writes `pd-check.rds` in
  `analysis/data/carryover-sensitivity/g9-prognostic-stress/`
  (Section 3.3).
- `analysis/scripts/carryover-sensitivity/31-prognostic-coupling-ceiling.R`
  writes `.../g9-prognostic-stress/coupling-ceilings.rds` (Section 5).
- Job 28: `28-run-g9-prognostic-stress.R`, results in
  `.../g9-prognostic-stress/` (commit `a69050073`). Job 30:
  `29-run-g9-random-effects.R`, results in
  `.../g9-random-effects/` (commit `3d1ee857a`). Both were generated
  under the old repair default; their strength-0.6 rows are the
  artifact described here and are retained for the record.
- The stop-by-default change and the zero-variance handling: commit
  `4aa4d0ef4`; tests in `inst/tinytest/test_sigma.R`.
- Both scripts run from the repository root in seconds on the host R
  session and need no simulation.

## 8. Caveats, unsolicited

- The ceilings in Section 5 are for one configuration (Hybrid,
  $t_{1/2} = 0.5$, $\rho = 0.7$, $c_1 = 0.2$). They move with every
  one of those settings; the compendium shows the $c_{bm}$ ceiling
  alone varies from 0.25 to 0.52 across the project's constructs.
  Script 31 must be rerun for any other configuration.
- The attenuation argument in Section 3.3 is a population
  calculation. It explains the direction and size of the observed
  loss and matches the fitted means, but the simulated estimates are
  also subject to the usual finite-sample behavior of `G9`'s
  AIC-selected half-life and CR2 standard errors, which this note
  does not separate out.
- `G9`'s nominal Type I error at the strength-0.6 null cells was
  measured under the repaired matrix, so it is a result about a
  modified DGP, not about the intended one. The strength-0.3 null
  results are unaffected.
- The claim that no analysis model can recover the lost signal is
  exact for the population contrast but was checked empirically only
  for `G9` and its three random-effects variants, and for the
  specifications in job 28; it is not a theorem about every
  conceivable estimator of a different estimand.
- Parity tests between the implementation collections were not
  rerun after the default change.

## References

- `analysis/report/compendium-mathematics.tex`, section "Validity of
  covariance binding: the coupling ceiling", equation `ce-ceiling`.
- `docs/34-cbm-feasible-range.md` and
  `docs/35-separable-response-covariance.md`, the single-coupling
  ceiling and its dependence on the correlation construct.
- `docs/57-data-chosen-carryover-half-life.md`, Section 7.3, the four
  prognostic pathways.
- `docs/60-trials-submission-gap-analysis.md`, 2026-10-26 update,
  the chronology of jobs 28 and 30 and the default change.
- `docs/61-cross-factor-correlation-bug-paper02-driver.md`, the
  $c_1 = 0.2$ correction that the ceilings in Section 5 assume.
- `docs/64-f1a0b-first-pass-evaluation.md`, the strength-0.3
  prognostic cells for the `F`-family, which Section 3.3 confirms are
  feasible.
- `analysis/scripts/quick-sim/cbm-ceiling/01-cbm-ceiling.R`, the
  single-coupling ceiling script this note's script 31 generalizes.
