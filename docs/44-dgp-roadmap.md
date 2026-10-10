# From the Published Process to a Defensible Data-Generating Process: A Step-by-Step Roadmap {.unlisted .unnumbered}
*2026-10-04 17:51 PDT*

**Author.** pmsimstats team

**Purpose.** This paper sets out, in a recommended order, the changes
that take the data-generating process (DGP) of Hendrickson et al. [1]
(`orig`, commit `58b32a9`) to the refined process, the DGP this project
now recommends (paper 01, Table `tab:spec`; named in
`analysis/report/NOTATION.md`), and beyond it to the pharmacodynamic model proposed in `docs/38`. For
each step it gives:

- what the step changes and which defect it removes;
- the evidence and its status;
- what it costs, in feasible effect size, comparability and analysis;
- whether the package (`R/generateData.R`) already implements it.

The rationale for individual changes is developed elsewhere: `docs/32`,
`docs/34`, `docs/35`, `docs/36`, `docs/37`, `docs/38` and
`docs/54`. This paper puts them in sequence and
says where the project stands.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

Notation follows `analysis/report/NOTATION.md`. Symbols marked *local*
are defined here. (Notation, 2026-10-09: the compound-symmetry
correlation is $\rho_{cs}$, with $\rho$ kept for the AR(1) correlation
per week; $\mu_{\text{stop}}$ is now the canonical $\mu_{\text{last}}$.)

| Symbol | Meaning | Status |
|----------------|----------------------------------------|----------|
| $BR$, $PB$, $TV$ | drug, expectancy and natural-course responses | canonical |
| $B_i$, $b_i$ | biomarker; standardized biomarker | canonical |
| $c_{bm}$ | covariance-moderation strength (correlation of $B$ with $BR$ where coupled) | canonical |
| $\beta_{bm}$ | mean-moderation multiplier | canonical |
| $c_{bm,TV}$, $c_{bm,PB}$ | prognostic biomarker-$TV$ and biomarker-$PB$ correlations (`docs/46`) | canonical |
| $\beta_{bm}^{TV}$ | mean moderation of the natural course by $b_i$ | canonical (component slope) |
| $N$ | total sample size across paths | canonical |
| $\rho_{cs}$ | common within-component correlation under CS, the same at every lag (0.8 in the published process) | canonical |
| $\rho$ | AR(1) within-component correlation per week, $\rho^{\lvert w_t - w_s\rvert}$ | canonical |
| $c_1$, $c_\times$ | cross-component correlation at the same visit and at different visits (code `c.cf1t`, `c.cfct`) | canonical |
| $t_{1/2}$, $\lambda$ | carryover half-life; $\lambda = \ln 2 / t_{1/2}$ | canonical |
| $t_{od}$, $t_{sd}$ | time on drug; time since discontinuation | canonical |
| $g_t$ | coupling at visit $t$: the fraction of $c_{bm}$ applied | canonical |
| $c^\ast_{bm}$ | ceiling: the largest $c_{bm}$ with every path's matrix positive definite | canonical |
| $w_t$, $w_s$ | calendar weeks of visits $t$ and $s$ | canonical |
| $K$, $A$ | $3 \times 3$ component correlation; common time correlation of the separable block $K \otimes A$ | canonical |
| $\mu_{\text{last}}$ | $BR$ mean at the last on-drug visit, the anchor of anchored carryover (written $\mu_{\text{stop}}$ before 2026-10-09) | canonical |
| $h(t)$ | population drug-response time course on $[0, 1]$ (`docs/38`'s $f(t)$) | local |
| $\Delta_h$, $\bar h_{\text{on}}$, $\bar h_{\text{off}}$ | contrast-weighted exposure $\bar h_{\text{on}} - \bar h_{\text{off}}$ (`docs/38`'s $F$); on- and off-drug means of $h(t)$ | local |
| $\omega_i$, $\bar\omega$ | participant drug sensitivity (`docs/38`'s $\theta_i$); its population mean | local |
| $\zeta$ | sensitivity moderation by the biomarker (`docs/38`'s $\kappa$) | local |
| $v_i$, $\sigma_v$ | unexplained sensitivity; its SD | local |
| $\xi_{it}$ | serial noise in $BR_{it}$ in the `docs/38` model | local |

Table: Notation index and glossary of symbols used in this document

$\omega_i$ and $\zeta$ are renamed from `docs/38`, whose $\theta$,
$\kappa$ and $z$ collide with canonical symbols (a generic interaction
coefficient, the standard-error ratio, and the standardized biomarker
$b_i$). $v_i$ and $\xi_{it}$ were written $u_i$ and $\varepsilon_{it}$
in this document before 2026-10-07; they were renamed because NOTATION
reserves $u_i$ for the participant random intercept and
$\varepsilon_{it}$ for the observation-level residual. $h(t)$ and
$\Delta_h$ are renamed from `docs/38`'s $f(t)$ and $F$: $f(t)$ is the
time adjustment of the analysis model in `docs/41`, `docs/42` and
`docs/45`, and $F$ would be read as configuration F.

- **Coupling rule.** How $g_t$ is set. *Step* (orig): $g_t = 1$ wherever
  the mean of $BR$ is nonzero. *Graded*: 1 on drug,
  $e^{-\lambda t_{sd}}$ off drug after exposure, 0 before any exposure.
- **CS.** Compound symmetry: correlation $\rho_{cs}$ at every lag.
- **AR(1).** Correlation $\rho^{|w_t - w_s|}$ between visits in weeks
  $w_t$, $w_s$.
- **Separable.** The response block is $K \otimes A$: a $3 \times 3$
  component correlation $K$ times a common time correlation $A$, so
  cross-component correlations decay with time like within-component
  ones (`docs/35`).
- **Configuration.** The lettered constructs of `docs/36`:
  - A, published;
  - B, AR(1) with the decaying (`covar`) cross block and step coupling;
  - C, separable AR(1) with graded coupling;
  - D, the published compound-symmetry construct (constant cross
    block) with mean moderation;
  - E, separable AR(1) (AR(1) within component, cross block
    $c_1\rho^{|w_t - w_s|}$) with mean moderation and no covariance
    coupling;
  - F, AR(1) with the decaying (`covar`) cross block and graded coupling.

  (Clarified 2026-10-09: `docs/36` formerly labeled E "Adj. 1 + 2 + 3"
  and `NOTATION.md` "all three". Graded coupling, adjustment 2, has no
  effect in E, whose biomarker coupling vector is zero; E is
  adjustment 1 with the separable cross block plus adjustment 3, as
  arm E of `hendrickson-problems/04-power-simulation.R` builds it.)
- **The refined process.** The recommended DGP (paper 01, Table
  `tab:spec`): configuration C at $\rho = 0.7$ and $c_1 = 0.2$, with
  anchored carryover, a stop on matrices that are not positive definite
  and randomized path allocation. Configuration C alone is only its
  covariance and coupling construct.
- **Problem 1 and problem 2.** Correlation matrices that are not
  positive definite and are silently repaired; the collapse of power
  under negligible carryover (`docs/36`).

## 1. Summary

1. **Eight steps take the published process to a generative model.**
   They are:
   1. fail loudly instead of repairing;
   2. graded coupling;
   3. anchored carryover;
   4. AR(1) correlation;
   5. separable cross-component structure with a factored generator;
   6. choosing $\rho$ and the strength inside the ceiling;
   7. restart from the residual response;
   8. the pharmacodynamic sensitivity model of `docs/38`.

   Steps 1 to 6 give the recommended near-term construct, the refined
   process: configuration C with anchored carryover and a stop on
   matrices that are not positive definite, at $\rho = 0.7$ (paper 01
   adds randomized path allocation and options). Steps 7 and 8 are a proposal awaiting implementation and
   validation (Section 3).
2. **The order matters.** Graded coupling alone, under compound
   symmetry, converts problem 2 into problem 1: the ceiling stays at
   0.256 while the published effect sizes are 0.3 and 0.6. AR(1) alone
   leaves the step rule's collapse in place and, on the published grid,
   leaves more matrices to repair (48% against 37%). Each step is
   placed after the steps it depends on (`docs/36`, Section 6.6).
3. **The package already implements steps 2, 3 and 4.** It does not
   implement step 1 fully: the Cholesky step still repairs silently.
   It also lacks step 5 (it uses the decaying cross-factor form,
   configuration F) and steps 7 and 8 (Section 4).
4. **The near-term target is the refined process** (configuration C
   with anchored carryover). That is separable AR(1), graded coupling
   and anchored carryover, generated by the factored generator, which
   cannot simulate an invalid model. At $\rho = 0.7$ it admits
   $c_{bm} = 0.45$ in every design and half-life with a margin of at
   least 0.02. At $\rho = 0.8$, the published value of $\rho_{cs}$, the
   Hybrid ceiling of configuration C is 0.440 (`docs/36`; ceilings recomputed for this
   paper).
5. **Each step needs an analysis that keeps pace.** AR(1) data need a
   serial-correlation or unstructured working covariance with a
   small-sample correction. Graded coupling and anchored carryover make
   an exposure-weighted drug indicator the natural analysis; the
   covariance-study runs now in progress evaluate both (`docs/43`).
   *Update (2026-10-07):* those runs are complete (`docs/45`), and a
   stress test on configuration C at $\rho = 0.7$ (the covariance
   construct of the refined process, with the published recursive
   carryover; see the status note) has selected the analysis
   (`docs/46`, `docs/49`; Section 5).

**Status note (2026-10-07).** The step statuses below were rechecked on
2026-10-07. The package's data-generating code is unchanged since this
paper was written (`R/generateData.R` last changed on 2026-09-08,
inspected with the git log), so the package column of Section 4 still
holds. Configuration C at $\rho = 0.7$, as used by the covariance
study and by the stress test (`11-strawman-stress-test.R`), still keeps
the published recursive carryover and the fixed published allocation
(inspected), so neither study ran the refined process in full. (Update
2026-10-09: both drivers now default to anchored carryover, and the
stress test to randomized allocation; the stored results are still the
recursive runs, and the anchored reruns are pending; see Step 3.) The
changes since 2026-10-04 are on
the analysis side (Section 5) and in the stress-test driver's new
prognostic mechanisms (Section 4).

## 2. The starting point: orig

The published process (`docs/36`, Section 3; `docs/40`, Sections 3
to 6):

| Feature | orig |
|----------------------------|----------------------------------------|
| Joint distribution | one 26-variable multivariate normal per path |
| Within-component correlation | CS, $\rho_{cs} = 0.8$ at every lag |
| Cross-component correlation | $c_1 = 0.2$ same visit, $c_\times = 0.1$ otherwise |
| Biomarker coupling | step rule: $c_{bm}$ wherever the $BR$ mean is nonzero |
| Carryover on the $BR$ mean | recursive, cumulative $t_{sd}$, over-decays |
| Re-exposure | time on drug restarts at zero |
| Off-drug $BR$ variance | full (SD 8) at every visit |
| Invalid matrices | silently repaired (`make.positive.definite`) |

Table: Features of the published process (orig)

Its two known defects share one cause, the step rule. Without carryover
the coupling has an on/off contrast, and the matrix is invalid for any
$c_{bm}$ above 0.256. With any carryover the coupling is constant, the
matrix is valid, and the interaction disappears (`docs/36`, Section 5.4;
`docs/40`, Section 5).

## 3. The steps

### Step 1: fail loudly instead of repairing

**Change.** Remove every silent repair. A cell whose matrix is not
positive definite stops, or is reported as infeasible, instead of being
simulated from a nearby matrix.

**Defect removed.** The repair changes the model that is simulated. At
$c_{bm} = 0.6$ the on-drug correlation actually simulated is 0.54 to
0.57, off-drug visits acquire 0.03 to 0.04, and other correlations move
by up to 0.10 (`docs/36`, Section 4.3). The published no-carryover
results at 0.6 describe a model other than the one stated.

**Cost.** Some published cells become infeasible. That is a statement
about the construct, not a loss.

**Package status: partial.** `validateParameterGrid()` checks
feasibility. Repairs are still silent in two places:

- with `makePositiveDefinite = TRUE` the repair warns only when
  `verbose = TRUE`;
- whatever that setting, the Cholesky factorization falls back to
  `make.positive.definite()` on failure, so a matrix that is not
  positive definite is always simulated.

The project's quick-sim drivers (`04-power-simulation.R`) skip
infeasible cells instead.

### Step 2: graded coupling

**Change.** Couple the biomarker to $BR$ in proportion to the drug
effect that remains: $g_t = 1$ on drug, $e^{-\lambda t_{sd}}$ off drug
after exposure, and 0 before any exposure.

**Defect removed.** Problem 2. Power now declines smoothly and
monotonically with the half-life, and at the published half-lives (0.1
and 0.2 weeks) the loss is under 0.01 (`docs/37`, Proposition 3,
closed form validated by simulation).

**Cost.** Under compound symmetry the ceiling barely moves (0.256
without carryover, 0.307 at one week in Hybrid), so above 0.256 graded
coupling alone turns problem 2 into problem 1 in every cell. This step
must be followed by steps 4 and 5 if larger effects are to be studied.

**Package status: implemented** (`buildSigma`, coupling
$c_{bm}e^{-\lambda t_{sd}}$ off drug; pre-exposure visits keep 0).

### Step 3: anchored carryover on the drug-response mean

**Change.** Off drug, decay the $BR$ mean from its value at
discontinuation, $\mu_{\text{last}} \, 2^{-t_{sd}/t_{1/2}}$, instead of
recursively decaying the previous visit's already-decayed mean with a
cumulative $t_{sd}$.

**Defect removed.** The recursion counts elapsed time twice, so the
mean decays faster than the stated half-life. At $t_{1/2} = 1$ week, two
weeks after stopping, the published mean is 1.27 where anchored decay
gives 2.55 (`docs/40`, Section 3; `docs/32`, Section 6).

**Cost.** Results at a nominal half-life differ from the published ones,
since the stated half-life now means what it says.

**Package status: implemented** (`buildSigma`, the `last_on` anchor).
On 2026-10-07 the quick-sim reproductions of orig kept the recursion
deliberately, and so did the configuration C drivers
(`04-power-simulation.R`, `08-covariance-study.R`,
`11-strawman-stress-test.R`; inspected 2026-10-07). *Update
(2026-10-09):* since 2026-10-08 those three drivers default to anchored
carryover, and `--carryover recursive` reproduces the earlier runs
(`docs/53`; `analysis/scripts/quick-sim/carryover-closed-form/RERUN-anchored.md`).
The stored covariance-study and stress-test results are still those of
the recursive runs; their anchored reruns are pending on the
simulation server.

### Step 4: AR(1) within-component correlation

**Change.** Within each component, correlation $\rho^{|w_t - w_s|}$, with
$\rho$ per week, in place of a constant $\rho_{cs}$ at every lag.

**Defect removed.** Compound symmetry says that symptoms one week apart
and sixteen weeks apart are equally correlated, which is implausible
for symptom trajectories. It is also the source of the low ceiling:
under CS the smallest eigenvalue of the response block is
$(1 - \rho_{cs}) - (c_1 - c_\times) = 0.1$ for every design (`docs/36`,
Section 4.2).

**Cost.**

- **Applied without the separable cross-component form, AR(1) helps
  little and can hurt.** It raises the ceiling where the coupling has an
  on/off contrast and lowers it where the coupling is constant (OL
  falls from 0.893 to 0.576). On the published grid it leaves more
  matrices to repair (48% against 37%).
  (Corrected 2026-10-08: this bullet first said "applied with the
  published cross-component form". The values 0.576 (OL) and 0.291
  (Hybrid, $t_{1/2} = 0$), and the 48% repair count, are those of
  `docs/36` configuration B, which pairs AR(1) with the decaying
  cross-component block $c_\times\rho^{|w_t - w_s|}$ (`docs/36`,
  Appendix A.6), not the published constant block $c_\times$. With
  AR(1) and the published constant block the ceilings are 0.535 (OL)
  and 0.363 (Hybrid, $t_{1/2} = 0$), still below the compound-symmetry
  OL value 0.893, so the direction of the statement holds. "AR(1)
  alone" in this paper and in `docs/36` therefore means AR(1) with the
  decaying cross block; in `docs/54` it means
  AR(1) with the constant block. See
  `analysis/report/compendium-mathematics.tex`, Section "Ceilings by
  serial correlation and by cross-component form".)
- **The analysis must change too.** On AR(1) data the published
  random-intercept analysis rejects 14% to 23% of null crossover trials
  (`docs/36`, Section 6.8; covariance study, 0.14 to 0.16).

**Package status: implemented** (`buildSigma`: `rho^tg` within
components, and the decaying cross-component form $c_\times\rho^{|w_t -
w_s|}$, configuration F).

### Step 5: a separable cross-component structure and the factored generator

**Change.** Make the response block separable, $K \otimes A$, with
cross-component correlations $c_1 \rho^{|w_t - w_s|}$ that decay with
time exactly as the within-component ones do. Generate data from the
factorization, without forming the full matrix:

1. three AR(1) recursions;
2. a $3 \times 3$ mixing step;
3. a conditional draw of the biomarker (`docs/36`, Section 6.7).

**Defect removed.**

- **The response block is valid for every schedule.** A Kronecker
  product of valid matrices is valid. The decaying cross-factor form
  of step 4 has no such guarantee: its smallest eigenvalue falls to
  0.010.
- **The ceiling rises.** The Hybrid ceiling without carryover goes from
  0.291 (configuration F) to 0.440 at $\rho = 0.8$ (`docs/35`; `docs/36`,
  Section 6.4).
- **An infeasible effect size stops the simulation instead of being
  repaired.** The factored generator's conditional variance is positive
  exactly below the ceiling.

**Cost.** Separability changes feasibility, not power: with graded
coupling, configurations C and F give the same power within $\pm 0.04$
at $\rho = 0.8$ once each test is judged against its own null
(`docs/36`, Section 6.8).
It is a modeling assumption: all components share one time correlation.

**Package status: not implemented.** The package uses configuration F.
Configuration C exists in the quick-sim drivers (`04-power-simulation.R`,
arm C) and was used by the covariance study at $\rho = 0.7$.

### Step 6: choose the correlation and the effect size inside the ceiling

**Change.** A decision, not code: set $\rho$ and the grid of $c_{bm}$
so that every cell is feasible with a margin.

**Evidence.** Exact graded-coupling ceilings, minimum over designs,
recomputed for this paper (`ceilings-graded-by-rho.csv`):

| $\rho$ ($\rho_{cs}$ for CS) | CS | AR(1), decaying cross (F) | Separable AR(1) (C) |
|--------|--------|--------|--------|
| 0.8 | 0.256 | 0.291 | 0.440 |
| 0.7 | 0.346 | 0.452 | 0.474 |
| 0.6 | 0.410 | 0.457 | 0.462 |
| 0.5 | 0.458 | 0.448 | 0.448 |

Table: Exact graded-coupling ceilings by $\rho$ (by $\rho_{cs}$ for CS) for CS, AR(1) and separable AR(1) constructs

These are the binding cells, at $t_{1/2} = 0$. At $\rho = 0.8$ and
$0.7$ carryover raises or keeps every ceiling. At lower $\rho$ a few
long-half-life cells fall slightly: at $\rho = 0.5$ the CO ceilings at
$t_{1/2} = 2$ are 0.003 to 0.004 below their no-carryover values, under
both AR(1) forms (verified in `ceilings-graded-by-rho.csv`, 2026-10-07;
the earlier sentence said carryover raises every ceiling). The minimum
over designs is still at $t_{1/2} = 0$ in every row of the table.

**Recommendation.** Configuration C at $\rho = 0.7$, the covariance
construct of the refined process, admits $c_{bm} = 0.45$,
the compendium's reference value (NOTATION.md, rule 2), everywhere with
a margin of at least 0.02. Alternatively keep $\rho = 0.8$, the
published value, and cap $c_{bm}$ at 0.4. Avoid running within about 0.02 of
a ceiling: there the conditional variance of $BR$ given the biomarker
approaches zero and results become fragile.

**Cost.** $\rho$ is itself an unestimated assumption about PTSD symptom
autocorrelation. Changing it changes power at a given strength, so
results across $\rho$ are not directly comparable.

### Step 7: restart from the residual response on re-exposure

**Change.** When the drug is restarted, let the response rebuild from
its residual level, not from zero. The simplest version keeps the
Gompertz onset curve and starts it at the level corresponding to the
residual response. The principled version is the turnover model of
step 8 (`docs/38`, Section 3).

**Defect removed.** Time on drug resets at every re-exposure, so the
Hybrid crossover's four-week on-drug period recovers only 39% of the
ceiling (4.3 of 11 points) even when the previous exposure ended weeks
earlier with a near-full response (`docs/40`, Section 4). For a drug
whose benefit persists, that understates the crossover's on-drug
response.

**Cost.** It interacts with carryover: a long response half-life raises
both the residual off-drug level and the restart level. The
pharmacology is uncertain, so the step should be reported as a
sensitivity, not imposed.

**Package status: not implemented** (`buildtrialdesign` resets $t_{od}$).

### Step 8: the pharmacodynamic sensitivity model (`docs/38`)

**Change.** Replace the matrix specification of $BR$ with a generative
model:

$$
BR_{it} = \omega_i\, h(t) + \xi_{it},
\qquad
\omega_i = \bar\omega\,(1 + \zeta b_i) + v_i ,
$$

with $h(t)$ the turnover time course (anchored decay, restart from the
residual), $\omega_i$ the participant's drug sensitivity moderated by
the standardized biomarker, $v_i$ unexplained sensitivity, and
$\xi_{it}$ serial noise. The mean, the variances and every
correlation follow (`docs/38`, Section 5).

**Defects removed.**

- **Moderation is proportional to the response:** small early in
  titration, decaying with the response off drug, never a step.
- **The off-drug drug-response variance shrinks with the response.**
  Every current construct gives $BR$ its full SD of 8 at off-drug
  visits with zero mean. `docs/38` judges this the least defensible
  feature of the construct, and it adds noise to every on-off contrast.
- **Validity by construction.** The covariance is a rank-one sensitivity
  term plus a valid noise covariance. There is no ceiling on effect
  size, because correlations are outputs, not inputs.

**Consequences for the interaction test** (`docs/38`, Section 6; derived,
not simulated). Power is strictly increasing in the contrast-weighted
exposure $\Delta_h = \bar h_{\text{on}} - \bar h_{\text{off}}$. It is capped by
the share of sensitivity variation the biomarker explains. And the
contrast variance now depends on carryover.

**Cost.**

- **A new architecture.** All results through `generateData()` would
  need comparison before any are replaced.
- **New parameters:** $\bar\omega$, $\zeta$ (or the peak correlation),
  $\sigma_v$, the serial noise, and the response half-life. The last is
  unknown for prazosin in PTSD.
- **A modeling choice:** multiplicative moderation assumes the biomarker
  acts on sensitivity. A biomarker that acted on onset or offset rates
  would need a different model.

**Package status: not implemented; proposal.** Required checks before
use (`docs/38`, Section 9):

- size holds when $\zeta = 0$;
- results reduce to the no-carryover case as $t_{1/2} \to 0$;
- simulated moments match the derivations.

## 4. Status and the recommended target

| Step | orig | Package (`covar`) | Configuration C (quick-sim) | `docs/38` model |
|--------------------|----------|----------------|----------------|----------|
| 1 Fail loudly | no | partial (Cholesky fallback repairs) | yes (cells skipped) | yes (no matrix) |
| 2 Graded coupling | no (step) | yes | yes | yes (proportional) |
| 3 Anchored carryover | no | yes | yes, driver default since 2026-10-08; stored results recursive | yes |
| 4 AR(1) | no (CS) | yes | yes | yes (serial noise) |
| 5 Separable, factored | no | no (F) | yes | not needed |
| 6 $\rho$, strength in range | 0.3 and 0.6 exceed 0.256 | grid-dependent | $\rho = 0.7$, 0.45 | no ceiling |
| 7 Restart from residual | no | no | no | yes |
| 8 Sensitivity model | no | no | no | yes |
| Hybrid ceiling, $t_{1/2} = 0$ | 0.256 | 0.291 at $\rho = 0.8$ | 0.440 at 0.8; 0.480 at 0.7 | none |

Table: Status of each roadmap step for orig, covar, Configuration C and the docs/38 model

**Near-term target: the refined process (configuration C with
anchored carryover, at a chosen $\rho$; $\rho = 0.7$ in paper 01).**
This needs four changes:

1. **Add the separable form to `buildSigma`**, as a new value of the
   cross-component option alongside the current decaying form.
2. **Add the factored generator** as the sampler for the separable form.
   It stops on an infeasible effect size.
3. **Remove the Cholesky fallback,** or make it an error.
4. **Use anchored carryover in configuration C.** As of 2026-10-07 the
   quick-sim version kept the published recursion for comparability
   with orig; the package version already anchors.

*Status (2026-10-07):* none of the four changes has been made; the
table above was rechecked against the package and the configuration C
drivers on this date and is unchanged (inspected).

*Update (2026-10-09):* change 4 has since been made in the quick-sim
drivers: `04-power-simulation.R`, `08-covariance-study.R` and
`11-strawman-stress-test.R` default to anchored carryover (the stress
test also to randomized allocation) from 2026-10-08, with the
recursion kept behind `--carryover recursive` (`RERUN-anchored.md`).
The stored covariance-study and stress-test results were produced
with the recursion, and their anchored reruns are pending. Changes 1
to 3 concern the package and remain open.

**Longer-term target: the `docs/38` model,** implemented as a new
`dgp_architecture` value after the checks of step 8. It would be run
alongside the refined process until the two have been compared.

**Branches that remain available.** Mean moderation (`'mean_moderation'`,
constant or trajectory-scaled), the dual-channel architecture
(`'combined'`, manuscript 11), the contaminated-biomarker option
(`c.bm.pb`, manuscript 06), reduced component sets (manuscripts 13 and
14) and alternative response-curve families (manuscript 07) are
orthogonal to steps 1 to 6 and can be combined with the refined process.
One extension is missing: there is no parameter for a biomarker related
to natural course, so a biomarker that predicts $TV$ cannot be simulated
(`docs/31`, Section 18.1). *Update (2026-10-07):* this is still true of
the package, but the stress-test driver (`11-strawman-stress-test.R`)
now simulates prognostic biomarkers on configuration C at
$\rho = 0.7$: a constant
biomarker-$TV$ correlation $c_{bm,TV}$, a mean moderation of $TV$ that
grows along its curve ($\beta_{bm}^{TV}$), a biomarker-$PB$ correlation
$c_{bm,PB}$, and a biomarker-baseline correlation (`docs/46`, Section 5;
inspected). These are candidates for the package.

## 5. Analysis that keeps pace with the DGP

| DGP step | Analysis implication | Evidence |
|-------------------------|----------------------------------------|----------------|
| 2, 3 (graded coupling, anchored carryover) | an exposure-weighted drug indicator `Dbc` matches the coupling's decay; binary `Db` remains valid but loses power under carryover | `docs/37`; manuscript 02 (Hybrid, OL+BDC); `docs/45`, Section 7 (the `Dbc` run, completed) |
| 4 (AR(1)) | the working covariance must allow serial correlation: `corCAR1` or a structured covariance with CR2; the unstructured covariance with Kenward-Roger is slightly liberal and failed the stress test | covariance study (`docs/45`); stress test (`docs/46`); manuscripts 02 and 10 |
| 5, 6 (feasible strengths) | larger effects become testable; power comparisons should be made at feasible, matched strengths | `docs/36`, Section 6.8 |
| 7, 8 (pharmacodynamic model) | `Dbc` built from the response half-life; a nonlinear turnover model could estimate that half-life | `docs/38`, Section 8 |

Table: Analysis implications and supporting evidence for each DGP step

*Update (2026-10-07).* The table's evidence column was revised on this
date. The exploratory stress test on configuration C at $\rho = 0.7$,
with the published recursive carryover (`docs/46`;
41 cells, 28,000 trials, 14 analyses) selects F1: a random intercept
with `corCAR1` residuals and CR2 inference, the baseline row, linear
time, `Dbc` at an assumed half-life of 0.5 week, and the expectancy
term `De` with its biomarker interaction (`docs/49`). The unstructured
covariance with Kenward-Roger, which `docs/45` recommended, rejected
0.057 of null trials pooled and 0.087 at $N = 35$, and failed all three
original prognostic cells. The prognostic cells are the reason
the analysis needs the expectancy term: they are DGP extensions the
package does not yet have (Section 4).

**A coding warning for `Dbc`.** Off-drug visits before any exposure must
be coded 0, not $e^{-\lambda \cdot 0} = 1$. Manuscript 02's implementation
codes them 1, which gives CO's placebo-first path `Dbc` = 1 at every
visit. That plausibly explains its finding that exposure weighting is
markedly inferior in CO (inspected in
`analysis/scripts/carryover-sensitivity/simulation-core.R`; to be
corrected and rerun). *Update (2026-10-07):* `docs/45` (Section 8)
confirmed the effect on the same data (CO power 0.484 to 0.744 under
CS). The driver has since been corrected: `exposure_weight()` now codes
pre-exposure off-drug visits 0 (inspected). A rerun of manuscript 02's
CO cells was not confirmed in this review.

## 6. Comparability with the published results

Every step departs from the published process, so no corrected result
can be set directly beside a published figure. Keep a faithful
reproduction of orig as a reference arm (the published reference arm,
configuration A without repair; `docs/36` called it the strawman, a
word `docs/46` now uses for a provisional analysis), report every
corrected result alongside
it, and attribute each difference to a named step. Steps 1, 2 and 3
correct errors and need no defense beyond the evidence cited. Steps 4 to
8 are modeling choices and should be presented as such.

## 7. Evidence status and limitations

- **Verified.** The ceilings (exact computations, recomputed for this
  paper; Step 6 table rechecked against `ceilings-graded-by-rho.csv` on
  2026-10-07); the effects of the repair, the step rule, AR(1) and
  separability (`docs/36`); the closed-form results of `docs/37`; the
  covariance-study results.
- **Inspected.** Package status, read in `R/generateData.R` at the
  current working tree (rechecked 2026-10-07: the Cholesky fallback to
  `make.positive.definite()` is still present, and there is no
  separable cross-component option); the manuscript 02 `Dbc` coding.
- **Derived, not simulated.** Step 8 and its consequences.
- **Not evaluated.** Step 7 (time on drug restarts at zero on
  re-exposure in every current DGP); the refined process itself, that
  is, the combination of anchored carryover with configuration C.
- **Limitations.**
  - The roadmap concerns the covariance-moderation architecture. Mean
    moderation needs only steps 3, 4, 7 and 8.
  - The choice of $\rho$ and of the response half-life are open
    empirical questions for prazosin in PTSD.

## 8. References

1. Hendrickson RC, Thomas RG, Schork NJ, Raskind MA. Optimizing
   aggregated N-of-1 trial designs for predictive biomarker validation:
   statistical methods and theoretical findings. *Frontiers in Digital
   Health* 2020; 2:13. doi:10.3389/fdgth.2020.00013.

Project documents: `docs/31`, `docs/32`, `docs/34`, `docs/35`, `docs/36`,
`docs/37`, `docs/38`, `docs/40`, `docs/43`,
`docs/54`.
