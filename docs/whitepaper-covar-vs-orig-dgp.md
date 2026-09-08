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

**There are five differences, not three, and the appendix describes
the wrong version of the most consequential one.** The three
documented rows are the within-factor correlation form, the
cross-factor off-diagonal, and the interaction vector `b`. Two further
differences are undocumented: `orig` never assigns the
biomarker-response correlation at the first measurement occasion, and
it decays the response mean at twice the nominal rate through a
`scalefactor` parameter that `covar` lacks. The appendix's description
of `b` corresponds to code that is commented out in the vendored file;
the active code behaves in the opposite direction.

**The active `b` code is the better of the two, contrary to first
appearance.** The commented-out published version raises the off-drug
biomarker-response correlation to full strength as carryover grows,
inverting the effect residual exposure should have. The 2024 revision
repairs this. Its defects are the `scalefactor` default and the
dropped occasion, not the change of gate itself (Section 3.3.1).

**Only one of the four covariance-matrix changes is required to fix
positive definiteness.** A crossed decomposition over all 16
combinations of those four (Section 5) shows that switching compound
symmetry to AR(1) accounts for the entire improvement, raising the
worst-cell `c_bm` ceiling from 0.34 to 0.48. Applied alone it beats the full
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
And the `scalefactor` default halves the effective carryover
half-life, so the two processes are not comparable at a matched
nominal half-life regardless of correlation strength.

**The vendored analysis model was also replaced, and carries a
separate correctness bug.** The `orig` arm does not run the published
analysis either: the binary `Db` regressor became the continuous
`Dbc` in 2024. In that code every never-treated occasion, baseline
included, is coded as fully on drug whenever `carryover_t1half > 0`,
so the crossover design's never-discontinuing path contributes no
on-drug versus off-drug contrast at all. The package's own
`R/lme_analysis.R` already guards this case, so no result computed
through the package pipeline is affected (Section 3.6).

**The remedy is a corrected `orig`, not a choice between the two.**
Three of the five 2024 DGP changes are sound and should be kept; the
`scalefactor` default and the dropped occasion should be fixed; and
the single AR(1) change of Section 5 resolves positive definiteness
outright (Section 4.5).

## 2. Scope, method, and epistemic status

Findings below are labeled `verified` (code was executed and output
observed), `inspected` (read in source and confirmed by reading), or
`inferred`. Nothing in this document is asserted from memory of the
manuscripts.

Sources examined:

- `analysis/scripts/quick-sim/hendrickson-original-comparison/`
  (`vendored-hendrickson-generateData.R`,
  `01-hendrickson-orig-driver.R`, `02-pd-sweep.R`, `README.md`)
- `R/generateData.R` (`buildSigma`, the `covar` construction) and
  `R/lme_analysis.R` (the package's own analysis model)
- `implementations/original/R/generateData.R`
- `analysis/report/01-dgp-mean-moderation-vs-mvn/report.Rmd`,
  Appendix sections B.4 through B.7
- The upstream repository `github.com/rchendrickson/pmsimstats`,
  cloned in full (21 commits, 2020-02-22 to 2024-11-27)

Computations run for this document:

- Reconstruction of the `b` vector under both the active and the
  commented-out code paths, Hybrid path A, at `t_half` in {0.5, 1.0}.
- A PD sweep over `c_bm` in [0.10, 0.60] by 0.05, crossed with
  `t_half` in {0, 0.5, 1.0} and the three trial designs, recording the
  minimum eigenvalue of Sigma under each process, minimized across the
  design's randomization paths.
- Quantification of the distortion introduced by PD repair.
- A crossed decomposition over all 16 combinations of the four
  covariance-matrix changes (Section 5).
- Comparison of the carryover mean adjustment at `scalefactor` 1 and
  2 against a correct exponential decay from the last on-drug value
  (Section 3.5).
- Diffs of the vendored files against the upstream initial commit
  (`42ac030`), the last pre-acceptance commit (`3035581`), and HEAD
  (`06dac83`).
- Evaluation of the `Dbc` exposure regressor across all occasions of
  the CO and Hybrid designs, under both the vendored and the package
  analysis code, at `t_half` in {0, 1.0} (Section 3.6).

**Provenance correction.** The directory README states the vendored
files come from `3035581` (dated there as 2026-06-28; the commit is
actually 2020-06-27). They do not. Diffing establishes the base as
HEAD, `06dac83`, 2024-11-27, from which the vendored copy differs by
exactly the one documented `data.table` line. The README's claim that
no other lines were changed is therefore true of the file that was
vendored and false of the file it names. The `RON THOMAS`
annotations are upstream, introduced by RC Hendrickson in commit
`8609f12` (2024-05-06, 'Ron Thomas' edits, currated into working
code'), not local edits.

All computations used the package's own `extracted_rp` and
`extracted_bp` parameter sets, `rho = 0.7`, `c.cf1t = 0.2`,
`c.cfct = 0.1`, `N = 35`, `scalefactor = 2`.

## 3. The five DGP differences, and the analysis model

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

#### 3.3.1 Which of the two is correct

**Status: verified.** The natural reading of the paragraphs above is
that the commented-out code is the correct target and the active code
an unwarranted deviation from it. On assessment that reading is wrong,
and it should be stated plainly before the arguments are weighed.

The published gate tests `means[...] != 0`, where `means` is the
**carryover-adjusted** response mean. The adjustment runs before the
correlation loop, so at any off-drug occasion carrying residual
effect, the adjusted mean is nonzero, the gate opens, and the occasion
receives the **full** `c_bm`. The consequence is that increasing
carryover progressively converts off-drug occasions into
full-strength on-drug occasions in the covariance channel, until at
`t_half = 1.0` in the Hybrid design the vector is `c_bm * 1_n` and no
on-drug versus off-drug contrast survives at all.

That behavior is backwards. Carryover should attenuate the
biomarker-response coupling at off-drug occasions, not raise it to
on-drug strength. The published DGP therefore has the property that
more carryover means a *stronger*, not weaker, off-drug interaction
signal, which is not a defensible representation of residual drug
exposure.

The active code gates on `brtest`, computed before the carryover
adjustment (Difference 5, Section 3.5), so it identifies genuinely
off-drug occasions, and then grades the correlation downward rather
than holding it at full strength. **The direction is correct and the
published version was wrong.** Whatever its implementation defects,
the 2024 revision repairs a real error rather than introducing one.

#### 3.3.2 Assessment

**Arguments for the active form.** It fixes the inversion described
above, which is the substantive point. It is graded rather than a
step, which is the more plausible representation. And it ties the
correlation channel to the same mean trajectory the carryover
adjustment already modifies, which is internally consistent.

**Arguments against.** Its decay rate is set by `scalefactor`, whose
default of 2 is itself an error (Section 3.5). It drops occasion 1
(Section 3.4). And it is undocumented: nothing in Paper 01, the
vendored file's header, or the directory README describes it, so a
reader of the manuscript cannot know which behavior the `orig` arm
has.

**Assessment.** The serious finding here is documentary rather than
methodological. Paper 01's Appendix B.5 and B.6 describe the 2020
code, which the `orig` arm does not run, and the appendix's causal
story for `orig`, that carryover opens its gate and destroys the
identifying contrast, is an accurate account of a defect that the
running code has already fixed. Any interpretation in Paper 01 resting
on that mechanism requires revision. But the fix should not be to
restore the 2020 behavior; it should be to describe the 2024 behavior
accurately and correct its two remaining implementation flaws.

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

### 3.5 Difference 5: the scalefactor on the carryover mean

**Status: verified.**

Paper 01's Appendix B.6 records the carryover adjustment to the
response mean as an 'identical line in all three' processes. It is
not. `orig` carries a `scalefactor` multiplier that `covar` does not:

```r
## orig
brmeans[p] <- brmeans[p] +
  brmeans[p-1]*(1/2)^(scalefactor * d$tsd[p]/carryover_t1half)

## covar
brmeans[p] <- brmeans[p] +
  brmeans[p-1]*(1/2)^(d$tsd[p]/carryover_t1half)
```

`scalefactor` defaults to 2, so `orig` decays the response mean at
twice the nominal rate. This is a difference in the **mean** channel,
distinct from the covariance-channel differences above, and it affects
any power comparison between the two processes even though it does not
enter Sigma and so does not appear in the Section 5 decomposition.

The default is an error, not merely a choice. Comparing the adjustment
against a correct exponential decay from the last on-drug value, at
`t_half = 1.0` in the Hybrid design:

| occasion | t_sd | sf = 1 | sf = 2 | correct |
|---|---|---|---|---|
| BD3 | 1 | 5.094 | 2.547 | 5.094 |
| BD4 | 2 | 1.273 | 0.159 | 2.547 |
| COp | 4 | 0.268 | 0.017 | 0.637 |

As a ratio to the correct value: at `sf = 1` the first off-drug
occasion is exact (1.000x) and later occasions over-decay (0.500x,
0.420x). At `sf = 2` even the first occasion is wrong (0.500x), and
the fourth is understated by a factor of 38 (0.026x). A dose one week
old at a one-week half-life should retain half its effect; the shipped
default gives it a quarter.

The parameter is also undocumented in the strict sense: its Roxygen
entry reads `@param scalefactor TODO update when understand what this
does?`, so it was shipped with a behavior-changing default while
recorded as not understood.

The residual over-decay at `sf = 1` is a separate, pre-existing bug
inherited from the published code. The recursion multiplies the
previous **already-adjusted** value by `(1/2)^(cumulative t_sd)`,
double-counting elapsed time. The correct form uses either the
interval rather than cumulative `t_sd`, or cumulative `t_sd` applied
to the last on-drug mean.

**Arguments for.** None identified for the default of 2. The
parameter itself is a harmless generalization if defaulted to 1.

**Arguments against.** It breaks a case the published code got right,
it compounds an existing over-decay rather than fixing it, and it
silently changes the DGP for every caller that does not override it.

**Assessment.** This is the change that should be removed. Setting
`scalefactor = 1` restores correct behavior at the first off-drug
occasion at no cost; fixing the recursion would correct the rest.

### 3.6 The analysis model: lme_analysis

**Status: verified.**

The five differences above concern `generateData.R`. The comparison
arm vendors a second file, `vendored-hendrickson-lme_analysis.R`,
which supplies the analysis model. It is byte-identical to upstream
HEAD, so the README's statement that it is unpatched is exactly true.
Five commits changed it between `3035581` and HEAD, altering 110
lines. Three findings follow.

#### 3.6.1 The analysis model itself was replaced

The published analysis and the running analysis test different
coefficients:

| | 2020 (`3035581`) | 2024 (HEAD, what runs) |
|---|---|---|
| exposure regressor | `Db`, binary `tod > 0` | `Dbc`, continuous decay |
| formula term | `bm*Db` | `bm*Dbc` |
| coefficient extracted | `bm:DbTRUE` | `bm:Dbc` |

The change came in `f70f86d` (2024-05-06), the same day as the DGP
changes. The `orig` arm therefore does not run the published analysis
either. Its exposure regressor is structurally the same object this
program calls Exposure-weighted (`G3`).

This **corroborates** the retraction in Paper 02, Section 3.3. That
section argues the reference analysis is effectively Unadjusted
because the drivers leave the analysis-side half-life at its default
of zero, at which the decayed predictor collapses onto the binary
indicator. Verified directly: at `carryover_t1half = 0` the off-drug
branch evaluates `(1/2)^(t_sd/0) = 0` and the on-drug branch 1, which
is exactly `Db`. The reasoning in that section holds.

#### 3.6.2 Never-treated occasions are coded as fully on drug

`Dbc` is assigned to every `Db == FALSE` row as
`(1/2)^(scalefactor * t_sd / t_half)`. That expression is meaningful
only after discontinuation. For occasions before a patient's first
dose, `buildtrialdesign` zeroes `t_sd` (multiplying by `everondrug`),
so the expression returns `(1/2)^0 = 1`, the **maximum** exposure
value.

CO path B, the never-discontinuing arm, at `t_half = 1.0`:

| occasion | on drug | t_sd | Dbc |
|---|---|---|---|
| BL | FALSE | 0 | **1.0** |
| COa1 - COa4 | FALSE | 0 | **1.0** |
| COb1 - COb4 | TRUE | 0 | 1.0 |

All nine rows are coded `Dbc = 1`. Four pre-treatment occasions and
baseline are indistinguishable from fully on drug, the path
contributes no on-drug versus off-drug contrast at all, and untreated
baseline data enter the model as treated. Under Hybrid the damage is
one row per patient rather than five, but it is present in every
design whenever `carryover_t1half > 0`.

At `carryover_t1half = 0` the same rows evaluate `(1/2)^(0/0)` and
return `NaN`, so baseline is silently dropped by `na.omit` instead of
miscoded. Different symptom, same root cause. The current driver
passes no `carryover_t1half` and so hits the `NaN` branch, which is
why this has not surfaced as a visible failure.

This is a genuine correctness bug, upstream, and distinct from
anything in `generateData.R`.

#### 3.6.3 The package's own lme_analysis is not affected

**Status: verified.** `R/lme_analysis.R` in this package already
guards the case, and its comment names the failure mode precisely:

```r
# Guard against 0/0 = NaN when carryover_t1half == 0 or
# tsd == 0 for never-on-drug subjects; both collapse to Dbc=0.
if (op$carryover_t1half == 0) {
  data.m2[Db==FALSE, Dbc:=0]
} else {
  data.m2[Db==FALSE & tsd<=0, Dbc:=0]
  data.m2[Db==FALSE & tsd>0,
          Dbc:=((1/2)^(op$carryover_scalefactor*tsd/op$carryover_t1half))]
}
```

Comparing the two on CO path B:

| `t_half` | package | vendored |
|---|---|---|
| 0 | 0 0 0 0 0 1 1 1 1 | NaN NaN NaN NaN NaN 1 1 1 1 |
| 1.0 | 0 0 0 0 0 1 1 1 1 | 1 1 1 1 1 1 1 1 1 |

Finding 3.6.2 is therefore confined to the vendored comparison arm.
No result in Paper 01 or Paper 02 that uses the package's own
analysis pipeline is affected by it.

#### 3.6.4 Minor

`error(...)` is called in the guard against combining `simplecarryover`
with a half-life-based carryover. No such function exists in base R
(verified), so the guard halts with `could not find function "error"`
rather than its intended message. It was evidently never exercised.

Assessed and found unimportant: the `union(timeptnames, "BL")` change,
which prevents `BL` being selected twice with no behavioral effect
otherwise; and the refactor of the model formula from nested `if`
blocks to string concatenation with `eval(parse(...))`, which is
self-described in the source as 'poor programming form!' but is
equivalent. Commit `325314f` changes only whitespace in this file, and
`baf4127`, titled 'first commit. edit to lme_analysis', does not touch
the file at all.

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

2. **The `scalefactor` default is `c_bm`-independent and wrong.**
   `orig` decays the response mean at twice the nominal rate at every
   correlation strength (Section 3.5), understating residual exposure
   by a factor of two at the first off-drug occasion and by up to 38
   at the fourth. The two processes are therefore not matched on
   carryover at a matched nominal `t_half`, so any comparison between
   them confounds carryover severity with everything else.

3. **The dropped first occasion is `c_bm`-independent.** It removes an
   on-drug occasion from the interaction channel at every parameter
   setting, including the null.

4. **Silent repair breaks the correspondence between nominal and
   realized effect size** in exactly the region where a restricted
   study would be operating closest to the boundary.

5. **The fidelity argument does not survive inspection, but not in
   the way it first appears.** The strongest reason to prefer `orig`
   at any `c_bm` is that it reproduces the published reference. The
   vendored file does not: it is upstream HEAD (`06dac83`,
   2024-11-27), four years after the paper, not the pre-acceptance
   commit the directory README names. So `orig` as run is not
   Hendrickson et al. as published.

   The remedy is not to restore the published version. Section 3.3.1
   shows the 2020 interaction gate is defective: it raises the
   off-drug biomarker-response correlation to full strength as
   carryover grows, inverting the effect residual exposure should
   have. Pinning to the pre-acceptance commit would buy documentary
   fidelity at the cost of reinstating a scientifically wrong DGP.

Section 5 qualifies the first of these five reasons. A crossed
decomposition shows the PD objection is narrower than it appears: it
is an objection to compound symmetry alone, and one change fixes it.
The other four reasons are not PD arguments and are unaffected.

The fifth point deserves the most emphasis, because it dissolves the
trade-off rather than resolving it. One would ordinarily accept a
narrower feasible parameter range as the price of matching a published
method. Here neither branch is attractive: the current arm does not
match the published method, and the published method is not worth
matching. The way out is a corrected process rather than a choice
between the two existing ones, which is the subject of Section 4.5.

### 4.5 What would change the answer

Rather than choosing between the two existing processes, the better
course is a corrected `orig`. Assessing the five 2024 changes
individually (Sections 3.3 to 3.5) shows they are not of a piece:
three are sound and two are defective.

**Keep.** The `brtest` / `rawbrmeans` capture (Difference 5's
infrastructure), which correctly separates pre-carryover state from
post-carryover mean. The interaction gate replacement, which repairs
the inversion in the published code (Section 3.3.1). The `verbose`
diagnostic blocks, which are inert with respect to returned values and
are the mechanism by which the `scalefactor` problem would have been
caught.

**Fix.** Set `scalefactor = 1`, or remove the parameter (Section 3.5).
Remove the `p > 1` guard so occasion 1 is assigned (Section 3.4),
which requires handling `p = 1` explicitly since the ratio has no
predecessor there.

**Fix separately.** The carryover recursion double-counts elapsed time
by applying `(1/2)^(cumulative t_sd)` to an already-adjusted value.
This is inherited from the published code and is present in `covar`
as well, so correcting it changes both processes.

A process so corrected would be better than the published original,
better than current HEAD, and would isolate the correlation structure
as the substantive remaining difference from `covar`. Combined with
the single AR(1) change of Section 5, it would also resolve the PD
objection outright.

Two further conditions would independently strengthen any comparison:

- Restricting the study to `c_bm <= 0.30` and to designs with
  `t_half > 0`, under which the PD objection does not arise even
  without the AR(1) change.
- Running the correlation structure and the interaction gate as
  crossed factors, as Paper 01's own Appendix B.6 suggests. Section 5
  does this for positive definiteness; the power half remains
  undone.

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

3. **Correct the directory README's provenance paragraph.** The
   vendored base is upstream HEAD `06dac83` (2024-11-27), not
   `3035581`, and the date given for that commit (2026-06-28) is
   wrong in both year and day (it is 2020-06-27). The 'no other lines
   were changed' claim is accurate once the correct base is named.
   The `RON THOMAS` annotations are upstream, from `8609f12`, and
   should not be described as local modifications.

3b. **Do not revert to the published version to resolve this.**
   Section 3.3.1 shows the 2020 interaction gate is defective.
   Pinning to `3035581` would buy documentary fidelity at the cost of
   a scientifically wrong DGP. Correct the description and the two
   remaining implementation flaws instead (Section 4.5).

3c. **Report the `Dbc` never-treated bug upstream** (Section 3.6.2).
   It is a genuine correctness defect in `rchendrickson/pmsimstats`
   HEAD, it is independent of everything else in this document, and
   the fix already exists in this package's `R/lme_analysis.R` and
   can be offered directly. Until it is fixed, the `orig` arm must
   not be run with `carryover_t1half > 0`, since its crossover
   results would be meaningless.

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

The upstream provenance question is now settled by direct comparison
(Section 2), but one part of it is not. Commit `8609f12` is authored
and committed by RC Hendrickson under the message 'Ron Thomas' edits,
currated into working code'. The history cannot show which parts of
that commit originated with Ron Thomas and which arose in the
curation, so responsibility for the `scalefactor` default and the
`p > 1` guard cannot be assigned from the record. This matters only
for deciding with whom to raise them.

The upstream code has not been executable on a current stack for some
time: the pre-2021 `data.table` indexing idiom is present in all 21
commits including HEAD, so `generateData()` fails outright without the
vendored patch. The `orig` arm therefore reproduces upstream source
but not upstream behavior, since upstream does not run unmodified.

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
