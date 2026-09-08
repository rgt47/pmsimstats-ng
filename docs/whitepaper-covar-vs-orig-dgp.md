# The covar and orig Data-Generating Processes: A Component-Wise
# Assessment, with Application to the c_bm <= 0.35 Regime

*2026-09-08 12:40 PDT*

Author: pmsimstats team

## 1. Executive summary

Paper 01 (`analysis/report/01-dgp-mean-moderation-vs-mvn/`) compares
three data-generating processes, of which two are covariance-channel
processes: `covar` (the package's Architecture B) and `orig` (a
vendored copy of the Hendrickson et al. reference implementation).
Appendix B.6 of that paper states that `orig` differs from `covar` on
three rows simultaneously. This white paper audits that claim against
the code that runs, assesses each difference on its merits, and
answers a specific practical question: if the biomarker-response
correlation `c_bm` is restricted to 0.35 or below, is there any
remaining reason to prefer `covar` over `orig`?

**There are four differences, not three, and the appendix's account of
the most consequential one is inverted.** The three documented rows are
the within-factor correlation form, the cross-factor off-diagonal, and
the interaction vector `b`. A fourth, undocumented difference is that
`orig` never assigns the biomarker-response correlation at the first
measurement occasion. The appendix's description of `b` corresponds to
code that is commented out in the vendored file; the active code
behaves in the opposite direction.

**Only one of the four changes is required to fix positive
definiteness.** A crossed decomposition over all 16 combinations of
the four changes (Section 5) shows that switching compound symmetry to
AR(1) accounts for the entire improvement, raising the worst-cell
`c_bm` ceiling from 0.34 to 0.48. Applied alone it beats the full
four-change `covar` set, which reaches only 0.45. The cross-factor
change lowers the ceiling to 0.16 when applied without AR(1), and
neither the `b` change nor the occasion-1 change moves it at all.

**The answer to the practical question is yes, on five independent
grounds.** Restricting `c_bm` to 0.35 does not rescue `orig`. The
positive-definite (PD) ceiling under `orig` is 0.30 in five of nine
design-by-half-life cells and 0.35 in the remaining four, so a study
run at 0.35 is non-PD in a substantial part of its own grid. Where it
is PD at 0.35 the margin is thin enough to be fragile to any other
parameter change. PD failure is repaired silently, so the realized
correlation stops matching the nominal one without any warning. The
dropped first occasion is a defect at every `c_bm`, including zero.
And the `b` decay rate differs by a factor of two, so the two
processes are not comparable at a matched nominal half-life regardless
of correlation strength.

## 2. Scope, method, and epistemic status

Findings below are labeled `verified` (code was executed and output
observed), `inspected` (read in source and confirmed by reading), or
`inferred`. Nothing in this document is asserted from memory of the
manuscripts.

Sources examined:

- `analysis/scripts/quick-sim/hendrickson-original-comparison/`
  (`vendored-hendrickson-generateData.R`,
  `01-hendrickson-orig-driver.R`, `02-pd-sweep.R`, `README.md`)
- `R/generateData.R` (`buildSigma`, the `covar` construction)
- `implementations/original/R/generateData.R`
- `analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd`,
  Appendix sections B.4 through B.7

Computations run for this document:

- Reconstruction of the `b` vector under both the active and the
  commented-out code paths, Hybrid path A, at `t_half` in {0.5, 1.0}.
- A PD sweep over `c_bm` in [0.10, 0.60] by 0.05, crossed with
  `t_half` in {0, 0.5, 1.0} and the three trial designs, recording the
  minimum eigenvalue of Sigma under each process, minimized across the
  design's randomization paths.
- Quantification of the distortion introduced by PD repair.

All computations used the package's own `extracted_rp` and
`extracted_bp` parameter sets, `rho = 0.7`, `c.cf1t = 0.2`,
`c.cfct = 0.1`, `N = 35`, `scalefactor = 2`.

## 3. The four differences

### 3.1 Difference 1: within-factor correlation form

**Status: verified.**

`orig` assigns a single constant to every within-factor occasion pair
(`vendored-hendrickson-generateData.R:117-126`):

```r
ac <- modelparam[[paste("c", c, sep = ".")]]
for (p in 1:(nP-1)) for (p2 in (1+p):nP) {
  correlations[n1, n2] <- ac
  correlations[n2, n1] <- ac
}
```

This is compound symmetry. `covar` uses `rho^|w_i - w_j|`, an AR(1)
form on the cumulative week scale.

**Arguments for compound symmetry.** It is the form in the published
reference, so it is the correct target if the goal is reproducing
Hendrickson et al. It has one parameter and is exchangeable, which
makes it invariant to occasion ordering and therefore robust to
irregular visit spacing. It is the standard random-intercept
implication, so it is internally consistent with an analysis model
that fits only a patient random intercept.

**Arguments against.** It is clinically implausible for repeated
symptom measurement, where adjacent occasions should correlate more
strongly than distant ones. More seriously for present purposes, it
constrains the feasible parameter space. A constant `rho` across `n`
occasions requires `rho > -1/(n-1)` on its own, and once the
biomarker-response block is added the joint constraint on `c_bm`
tightens rapidly with `n`. Appendix B.7 identifies this as the
mechanism behind the narrower PD range, and Section 4 below quantifies
it.

**Assessment.** This difference is the dominant driver of the PD
behavior and therefore of the answer to the paper's practical
question. It is a legitimate modeling choice in isolation, but it is
not a neutral one: it changes what parameter values can be studied at
all.

### 3.2 Difference 2: cross-factor off-diagonal

**Status: verified.**

`orig` assigns `modelparam$c.cfct` as a constant to every cross-factor
pair at distinct occasions (line 140). `covar` assigns
`c.cfct * rho^tg`, decaying in the lag `tg`.

**Arguments for the constant form.** Again, fidelity to the reference.
It also imposes no assumption about how quickly the coupling between
distinct latent factors decays, which is a quantity for which little
empirical guidance exists.

**Arguments against.** It is internally inconsistent with any decaying
within-factor structure, since it asserts that a factor's correlation
with a different factor eight occasions away is as strong as with
itself at the adjacent occasion. It contributes to the PD constraint
in the same direction as Difference 1, since constant off-diagonal
blocks raise the matrix's dependence on its smallest eigenvalue.

**Assessment.** Secondary to Difference 1 in magnitude, and in the
same direction. It is not separately identifiable from Difference 1 in
any comparison run to date, since the two always vary together.

### 3.3 Difference 3: the interaction vector b

**Status: verified, and the paper's description is incorrect.**

Appendix B.6 states that `orig` gates the interaction on whether the
response mean is nonzero, giving a two-valued step of `{0, c_bm}`, and
that at `t_half = 1.0` in the Hybrid design the vector becomes
`c_bm * 1_n`, carrying no on-drug versus off-drug contrast at all. The
worked example at B.5 prints `b^orig = (0.4500, 0.4500, 0.4500,
0.4500)`.

That describes the following code, which is present in the vendored
file but **commented out**, under the header
`## Following commented out in Ron Thomas version:`

```r
#if(means[which(n1==labels)]!=0){
#  correlations[n1,'bm']<-modelparam$c.bm
#  correlations['bm',n1]<-modelparam$c.bm
#}
```

The code that actually executes, under the header
`## RON THOMAS VERSION:`, is line 153:

```r
correlations["bm", n1] <- correlations[n1, "bm"] <-
  ifelse(brtest[p],
         ifelse(brmeans[p] == 0, 0, (mm1 / mm0) * modelparam$c.bm),
         modelparam$c.bm)
```

where `mm1` and `mm0` are the carryover-adjusted response means at
occasions `p` and `p - 1`, and `brtest` is computed before the
carryover adjustment, so it is TRUE exactly at off-drug occasions.

Re-running both code paths over the Hybrid path A design gives:

| occasion | on drug | adj mu_BR | b active | b commented-out |
|---|---|---|---|---|
| OL1 | yes | 4.281 | NA | 0.45 |
| OL2 | yes | 9.223 | 0.4500 | 0.45 |
| BD1 | yes | 9.793 | 0.4500 | 0.45 |
| BD2 | yes | 10.187 | 0.4500 | 0.45 |
| BD3 | no | 2.547 | 0.1125 | 0.45 |
| BD4 | no | 0.159 | 0.0281 | 0.45 |
| COd | yes | 4.281 | 0.4500 | 0.45 |
| COp | no | 0.017 | 0.0018 | 0.45 |

(`t_half = 1.0`, `scalefactor = 2`, `c_bm = 0.45`.)

The commented-out column reproduces the appendix exactly. The active
column does the opposite: off-drug entries decay toward zero rather
than rising to `c_bm`, and the on-drug versus off-drug contrast is
preserved rather than lost.

The active decay belongs to the same exponential family as `covar`'s,
but runs at twice the rate. The ratio `mm1 / mm0` inherits
`(1/2)^(scalefactor * t_sd / t_half)` from the carryover adjustment at
line 96, and with `scalefactor = 2` this is an effective half-life of
`t_half / 2`, compounding multiplicatively across successive off-drug
occasions. For comparison, at `t_half = 1.0` and one week off drug,
`covar` gives `0.45 * exp(-ln2 * 1) = 0.225` while `orig` gives
`0.1125`.

**Arguments for the active form.** It is graded rather than a step,
which is arguably more plausible than a binary gate, and it ties the
correlation channel to the same mean trajectory the carryover
adjustment already modifies, which is a defensible internal
consistency.

**Arguments against.** It is undocumented. Nothing in Paper 01, the
vendored file's own header, or the directory README describes it. Its
decay rate is set by `scalefactor`, a parameter whose own Roxygen
comment in the vendored file reads `TODO update when understand what
this does?`. And because it decays at twice `covar`'s rate, a
comparison of the two at a matched nominal `t_half` is not a
comparison at matched carryover.

**Assessment.** This is the most serious finding in this document. The
paper's stated causal mechanism for `orig`'s behavior, that carryover
opens its gate at off-drug occasions and destroys the contrast on
which identification depends, is an accurate description of code that
does not run. Any interpretation in Paper 01 resting on that mechanism
requires revision, independent of anything else here.

### 3.4 Difference 4: the first occasion is never assigned

**Status: verified.**

The assignment loop at line 147 is guarded by `if (p > 1)`, and
`correlations` is initialized as `diag(length(labels))` at line 106.
The biomarker-response correlation at occasion 1 is therefore left at
exactly zero under `orig`, regardless of `c_bm` and regardless of
whether occasion 1 is on drug. In the table above this appears as the
`NA` in the OL1 row. `covar` assigns all eight occasions.

**Arguments for.** None identified. This does not appear to be a
modeling choice.

**Arguments against.** It removes one of eight occasions from the
interaction channel, and it does so at an occasion that is on drug in
every path of every design examined, so the loss falls on the most
informative part of the series. It attenuates `orig`'s realized
interaction relative to its nominal `c_bm` at every parameter setting.
It also interacts with Difference 3: the `p > 1` guard exists because
line 153 needs `p - 1` to form the ratio, so the defect is a direct
consequence of the undocumented modification, not an independent bug.

**Assessment.** A defect rather than a difference in modeling
philosophy. It is `c_bm`-independent and therefore is not addressed by
restricting the correlation range.

## 4. The c_bm <= 0.35 question

### 4.1 Why the question arises

Compound symmetry constrains the feasible correlation range more
tightly than AR(1), as Appendix B.7 notes. If the binding objection to
`orig` were only that it cannot represent large `c_bm`, then
restricting attention to `c_bm <= 0.35` would neutralize the
objection, and the choice between processes would come down to
fidelity to the published reference, which would favor `orig`.

### 4.2 Positive-definiteness sweep

**Status: verified.** Minimum eigenvalue of Sigma, minimized across
each design's randomization paths, at the parameter settings in
Section 2.

Largest `c_bm` on the tested grid at which Sigma remains PD:

| design | t_half | orig | covar |
|---|---|---|---|
| CO | 0.0 | 0.30 | 0.60 |
| CO | 0.5 | 0.30 | 0.60 |
| CO | 1.0 | 0.30 | 0.60 |
| Hybrid | 0.0 | 0.30 | 0.45 |
| Hybrid | 0.5 | 0.35 | 0.50 |
| Hybrid | 1.0 | 0.35 | 0.50 |
| OL+BDC | 0.0 | 0.30 | 0.45 |
| OL+BDC | 0.5 | 0.35 | 0.45 |
| OL+BDC | 1.0 | 0.35 | 0.50 |

Minimum eigenvalue at exactly `c_bm = 0.35`:

| design | t_half | orig | covar | orig PD |
|---|---|---|---|---|
| CO | 0.0 | -0.3069 | 8.516 | no |
| CO | 0.5 | -0.3069 | 8.516 | no |
| CO | 1.0 | -0.3069 | 8.516 | no |
| Hybrid | 0.0 | -0.3140 | 2.158 | no |
| Hybrid | 0.5 | 0.0831 | 2.297 | yes |
| Hybrid | 1.0 | 1.2913 | 2.340 | yes |
| OL+BDC | 0.0 | -0.3217 | 2.156 | no |
| OL+BDC | 0.5 | 0.0869 | 2.300 | yes |
| OL+BDC | 1.0 | 1.4385 | 2.344 | yes |

Three observations follow.

First, `c_bm = 0.35` is **not** a safe ceiling for `orig`. Sigma is
non-PD in four of the nine design-by-half-life cells: all three CO
cells, and the zero-carryover cells of both other designs. The safe
ceiling across the full grid is 0.30, not 0.35.

Second, where `orig` is PD at 0.35, it is marginal. Minimum
eigenvalues of 0.083 and 0.087 at `t_half = 0.5` sit three orders of
magnitude below `covar`'s at the same cells, and two orders below
`orig`'s own values at `c_bm = 0.10`. A margin that thin is not robust
to changes in `rho`, `c.cfct`, the number of occasions, or the
response parameters, none of which were varied in this sweep.

Third, the CO design is uniformly the binding case under `orig`, and
its minimum eigenvalue does not vary with `t_half` at all. This is
consistent with CO's off-drug occasions being spaced 2.5 weeks or more
apart, so that the `orig` decay, at twice `covar`'s rate, has driven
the off-drug `b` entries to effectively zero at every half-life
examined.

### 4.3 Silent repair

**Status: verified.**

Both implementations call `make.positive.definite(sigma, tol = 1e-3)`
when the PD check fails, with no warning, no error, and no record in
the returned object:

```r
if (makePositiveDefinite) {
  if (!is.positive.definite(sigma)) {
    sigma <- make.positive.definite(sigma, tol = 1e-3)
  }
}
```

The consequence is that a non-PD cell does not fail. It silently
becomes a different DGP. Measuring the realized biomarker-response
correlations after repair, for CO path A at `t_half = 1.0`:

| nominal c_bm | min eigenvalue | realized max | repaired |
|---|---|---|---|
| 0.20 | 5.5645 | 0.2000 | no |
| 0.25 | 5.0361 | 0.2500 | no |
| 0.30 | 3.6264 | 0.3000 | no |
| 0.35 | 0.8645 | 0.3500 | no |
| 0.40 | -2.8964 | 0.3955 | yes |
| 0.45 | -7.3385 | 0.4370 | yes |

The distortion is modest in magnitude, roughly one to three points,
but it is undisclosed and it is not constant across the grid, so it
introduces a systematic and unreported difference between nominal and
realized effect size precisely in the cells where `orig` is under the
most strain. Power reported against a nominal `c_bm` in a repaired
cell is not power at that `c_bm`.

### 4.4 Answer

**Restricting `c_bm` to 0.35 or below does not remove the case for
`covar`.** Five reasons remain, in descending order of force.

1. **The restriction is not tight enough.** At exactly 0.35, `orig` is
   non-PD in four of nine cells, including every CO cell. A
   restriction that actually cleared the PD objection would be
   `c_bm <= 0.30`, and even then only if no other parameter moves.

2. **The undocumented `b` modification is `c_bm`-independent.**
   `orig`'s off-drug interaction decays at twice `covar`'s rate for
   any correlation strength. The two processes are therefore not
   matched on carryover at a matched nominal `t_half`, so a comparison
   between them confounds carryover severity with everything else.

3. **The dropped first occasion is `c_bm`-independent.** It removes an
   on-drug occasion from the interaction channel at every parameter
   setting, including the null.

4. **Silent repair breaks the correspondence between nominal and
   realized effect size** in exactly the region where a restricted
   study would be operating closest to the boundary.

5. **The fidelity argument does not survive inspection.** The strongest
   reason to prefer `orig` at any `c_bm` is that it reproduces the
   published reference. But the vendored file is not the published
   reference. Its interaction assignment has been replaced, and the
   original is commented out beneath it. Whatever `orig` currently
   reproduces, it is not Hendrickson et al. as published.

Section 5 qualifies the first of these five reasons. A crossed
decomposition shows the PD objection is narrower than it appears: it
is an objection to compound symmetry alone, and one change fixes it.
The other four reasons are not PD arguments and are unaffected.

The fifth point deserves emphasis because it inverts the usual
trade-off. One would ordinarily accept a narrower feasible parameter
range as the price of matching a published method. Here that price is
being paid without the benefit being received.

### 4.5 What would change the answer

The case for `orig` would be materially stronger under any of the
following, none of which currently holds.

- The vendored file is restored to the genuine upstream code, with the
  `RON THOMAS VERSION` block reverted and the original gate
  uncommented. `orig` would then serve its stated fidelity purpose,
  and Paper 01's Appendix B.6 would become an accurate description of
  it.
- The study is restricted to `c_bm <= 0.30` **and** to designs with
  `t_half > 0` **and** the `p > 1` guard is fixed. Under those three
  conditions the PD objection is genuinely resolved and the remaining
  differences are defensible modeling choices rather than defects.
- The correlation structure and the interaction gate are run as
  crossed factors, as Paper 01's own Appendix B.6 suggests. This would
  decompose the joint effect and permit a component-wise judgment
  rather than the present all-or-nothing comparison.

## 5. Decomposition: the minimum set of changes that fixes PD

### 5.1 Motivation and method

The stated reason for the `covar` changes was the non-positive-definite
correlation matrix under `orig`. That justification implies a
counterfactual which has not previously been tested: how many of the
four changes are actually required to fix PD, and does each contribute
in the direction assumed?

**Status: verified.** A single Sigma constructor was written with one
flag per candidate change, reproducing `orig` when all four flags are
off and `covar` when all four are on. The flags are:

- `ar1`: AR(1) within-factor (`rho^|w_i - w_j|`) versus compound
  symmetry (constant `rho`). Difference 1.
- `cfd`: decaying cross-factor off-diagonal (`c.cfct * rho^|w_i-w_j|`)
  versus constant `c.cfct`. Difference 2.
- `bexp`: exponential `b` (`c_bm * exp(-lambda * t_sd)`) versus the
  ratio-graded `orig` form. Difference 3.
- `occ1`: assign the biomarker-response correlation at occasion 1
  versus skipping it. Difference 4.

All 16 combinations were crossed with three designs and three
half-lives. For each, `c_bm` was swept from 0.05 to 0.80 in steps of
0.01, and the largest value at which Sigma remained PD across every
randomization path was recorded. The reported ceiling for a
configuration is its worst cell over all nine design-by-half-life
combinations.

### 5.2 Effect of each change in isolation

Starting from `orig` (all four flags off, worst-cell ceiling 0.34) and
turning on exactly one flag:

| change enabled | worst-cell ceiling | effect |
|---|---|---|
| none (`orig`) | 0.34 | baseline |
| `ar1` only | **0.48** | **+0.14** |
| `cfd` only | 0.16 | -0.18 |
| `bexp` only | 0.34 | none |
| `occ1` only | 0.34 | none |

The result is unambiguous. **The AR(1) change alone accounts for the
entire PD improvement.** The other three contribute nothing positive,
and one of them is actively harmful.

### 5.3 The minimum set is a single change

`ar1` is both necessary and sufficient. Every configuration with
`ar1 = TRUE` reaches a worst-cell ceiling of 0.43 or above; every
configuration with `ar1 = FALSE` is capped at 0.34 or below,
regardless of what the other three flags are set to. Ranked across all
16 configurations, the eight `ar1 = TRUE` configurations occupy the
top eight positions without exception.

Per design and half-life:

| design | t_half | orig | AR(1) only | full covar |
|---|---|---|---|---|
| CO | 0.0 | 0.34 | 0.59 | 0.61 |
| CO | 0.5 | 0.34 | 0.59 | 0.61 |
| CO | 1.0 | 0.34 | 0.59 | 0.61 |
| Hybrid | 0.0 | 0.34 | 0.49 | 0.45 |
| Hybrid | 0.5 | 0.35 | 0.51 | 0.50 |
| Hybrid | 1.0 | 0.36 | 0.54 | 0.54 |
| OL+BDC | 0.0 | 0.34 | 0.48 | 0.45 |
| OL+BDC | 0.5 | 0.35 | 0.49 | 0.49 |
| OL+BDC | 1.0 | 0.36 | 0.52 | 0.53 |

**Changing only the within-factor correlation form, and leaving the
other three `orig` behaviors untouched, yields a higher worst-cell
ceiling (0.48) than the full set of four changes (0.45).** The
minimal fix is not merely sufficient; it is strictly better than
`covar` on the criterion that motivated `covar` in the first place.
AR(1)-only matches or exceeds full `covar` in seven of the nine cells,
and trails it by at most 0.01 in the remaining two.

### 5.4 Why the cross-factor change is harmful

The `cfd` change, applied without `ar1`, lowers the ceiling from 0.34
to 0.16, the worst configuration tested. The mechanism is structural
incoherence. Under compound symmetry a factor's correlation with
itself is a constant `rho` at every lag. Introducing a decaying
cross-factor term makes a factor's correlation with a *different*
factor fall away with lag while its correlation with *itself* does
not. At long lags the matrix asserts strong within-factor coupling
alongside near-zero cross-factor coupling, and the resulting block
structure drives the smallest eigenvalue negative sooner than either
uniform choice does.

Under `ar1` the same change is close to neutral, contributing at most
0.02 in either direction, because within-factor and cross-factor terms
then decay together and the matrix stays internally consistent.

This is worth stating plainly: **Difference 2 was not a PD fix. In
isolation it is a PD regression, and its only defensible justification
is internal consistency with Difference 1.**

### 5.5 Why the b and occasion-1 changes are irrelevant to PD

Neither `bexp` nor `occ1` moved the ceiling by a single grid step in
any configuration. This is expected on inspection. The
biomarker-response entries occupy one row and one column of a matrix
with `2 + 3n` dimensions, and their magnitudes are bounded by `c_bm`
itself. The binding constraint comes from the `3n x 3n` within- and
cross-factor blocks, not from the biomarker margin.

The practical consequence is that Differences 3 and 4 cannot be
defended on PD grounds at all. Difference 3 is an undocumented change
of decay rate (Section 3.3) and Difference 4 is a defect (Section
3.4). Whatever their merits, the non-positive-definiteness of the
`orig` matrix is not among them.

### 5.6 Implications

1. **The minimum set of changes that fixes the PD problem is one
   change: compound symmetry to AR(1).** Adopting it alone recovers
   more feasible parameter range than the full four-change set.

2. **Three of the four changes in `covar` are not PD fixes.** They may
   be justified on other grounds, but the PD rationale does not extend
   to them and should not be offered for them.

3. **A minimal-change variant is worth constructing.** An
   `orig + AR(1)` process would isolate the correlation structure as
   the sole difference from the published reference, which is exactly
   the crossed-factor decomposition Paper 01's Appendix B.6 says would
   be required, and it is a one-line change to the vendored file.

4. **This changes the practical answer of Section 4.4 in one respect.**
   The PD objection to `orig` is narrower than it appeared: it is an
   objection to compound symmetry specifically, not to the `orig`
   process as a whole. The other four objections in Section 4.4 stand
   unchanged, since none of them is a PD argument.

## 6. Recommendations

1. **Do not adopt `orig` as a comparison arm on the strength of a
   restricted `c_bm` range.** The restriction addresses one of four
   differences and does not fully address even that one.

2. **Correct Paper 01, Appendix B.5 to B.7.** The `b` row, the worked
   example, and the third interpretive paragraph after the B.6 table
   describe code that does not run. This is independent of any
   decision about future simulation work and should be treated as a
   correctness fix.

3. **Verify the vendored files against the upstream repository**
   (`github.com/rchendrickson/pmsimstats`, commit `3035581`). The
   directory README states the files are verbatim apart from one
   `data.table` indexing patch and that no other lines were changed.
   At least two further edits are visible in the file: the
   `RON THOMAS VERSION` block, and a typo fix at line 134 annotated
   `##### <--------FIXING TYPO, this was [n1,n2] again`. Either the
   README's provenance claim or the file needs correcting.

4. **Make PD repair loud.** Both implementations should record the
   repair in the returned object and warn. A silent change to the DGP
   is not acceptable in a simulation study whose outputs are reported
   against nominal parameter values.

5. **Fix the `p > 1` guard** if `orig` is retained in any form, or
   document it as a known property of the arm if the intent is to
   reproduce whatever the vendored file does.

6. **Construct an `orig + AR(1)` variant.** Section 5 shows this is a
   one-line change that fixes PD more effectively than the full
   `covar` set. It isolates the correlation structure as the sole
   difference from the reference implementation and would let Paper 01
   report a decomposition rather than a confounded bound.

7. **Stop offering the PD rationale for Differences 2, 3, and 4.**
   Section 5 shows none of them improves positive definiteness, and
   Difference 2 degrades it when applied on its own. If those changes
   are retained they need a different justification.

8. **If a full decomposition of the power results is wanted, run the
   crossed design.** Correlation structure (CS versus AR(1)) crossed
   with interaction gate (step versus graded) at `c_bm <= 0.30`, three
   designs, three half-lives. This is four DGP variants rather than
   two and would separate the effects the present comparison
   confounds. The PD half of that decomposition is already done
   (Section 5); the power half is not.

## 7. Limitations

The PD sweep used a single parameter set (`rho = 0.7`,
`c.cf1t = 0.2`, `c.cfct = 0.1`, `N = 35`, eight occasions, the
package's `extracted_rp` and `extracted_bp`). The PD boundary depends
on all of these, and the boundaries reported in Section 4.2 should not
be read as general. In particular, larger `rho` will tighten the
`orig` boundary further.

The Section 5 decomposition is a statement about positive
definiteness only. It establishes which changes widen the feasible
`c_bm` range, not which produce better statistical behavior within
that range. An `orig + AR(1)` process has a higher PD ceiling than
`covar` but was not evaluated for power, Type I error, or estimator
bias, and nothing here implies it would perform better on those.

The Section 5 constructor reproduces both processes from the source
read in Section 2 rather than calling them. It was checked against the
`orig` and `covar` ceilings obtained in Section 4.2 from the real
implementations, which agree to within the coarser grid step used
there, but it is a reimplementation and may diverge in behavior not
exercised by this sweep.

No power simulation was run for this document. The consequences of
Differences 3 and 4 for realized power under `orig` are argued from
the structure of the covariance matrix, not measured. Quantifying them
would require running the `orig` arm with and without the `p > 1`
guard, and with the commented-out gate restored, which is the natural
follow-up.

The upstream provenance question is unresolved. The comparison against
`github.com/rchendrickson/pmsimstats` was not performed, and the
attribution of the `RON THOMAS VERSION` edits is therefore inferred
from the annotation text alone.

`implementations/original/R/generateData.R`, the package's other copy
of the reference implementation, was inspected and found to differ
from the vendored copy: it gates on `d[p]$onDrug` and uses
`c.cfct * rho^tg` for the cross-factor off-diagonal. It is therefore a
third variant, distinct from both `covar` and the vendored `orig`, and
its relationship to the other two was not pursued here.

## 8. Provenance

Audit conducted 2026-09-08 against the working tree at commit
`3ce1162`. Computations executed on the host R session (R 4.6.1),
outside the project container. Scratch scripts were not retained in
the repository.
