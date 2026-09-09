# 13-two-component-dgp

*2026-09-08 16:30 PDT*

## Working title

Data generation machinery for N-of-1 clinical trial simulations
under a two-component response decomposition.

## Origin

Paper 13 is a reduced-decomposition counterpart to paper 01
(`analysis/report/01-dgp-mean-moderation-vs-mvn/`). It asks paper
01's question, comparing mean moderation (Architecture A) against
covariance moderation (Architecture B), under a data-generating
process that carries only the pharmacological (BR) and time-variant
natural-history (TV) components. The placebo-belief component (PB)
is removed.

The motivation is twofold. The primary one is simplification of the
correlation matrix, from `(2+3n)` to `(2+2n)`. The second, which
determines *which* component to drop, is conformity with the
literature: a trend term is conventional in crossover and N-of-1
models, whereas a modelled placebo-belief trajectory is an
innovation of this program rather than an inherited convention.
The placebo literature (Hrobjartsson and Gotzsche) further regards
the belief and natural-history channels as separable only under
three-arm designs that most trials do not field, so a model that
keeps natural history and drops modelled belief makes the weaker
assumption. Section 1.1 of the manuscript develops this.

**Asymmetry worth knowing.** The two candidate reductions are not
equally safe. No biomarker-TV parameter exists, so dropping TV would
have been safe unconditionally. Dropping PB is safe only while the
optional `c.bm.pb` contamination parameter is zero, which is the
default and is what this program's simulations use. A study setting
it nonzero cannot use the reduced DGP.

## The two reductions differ on three axes

| axis | TV-drop (paper 14) | PB-drop (paper 13) |
|---|---|---|
| Safety | unconditional: no `c.bm.tv` parameter exists | conditional: holds only while `c.bm.pb` is 0 |
| Design structure | nothing becomes inert | expectancy weight becomes inert |
| Convention | deletes the trend term the literature expects | matches crossover / N-of-1 practice |

The first two favor paper 14, the third favors paper 13. Neither
dominates. Section 4.5 of both manuscripts sets them against each
other and declines to pick a winner: TV-drop for work internal to
this program, where contamination is live and the blinding contrast
is used; PB-drop for work aimed at the wider N-of-1 community. What
matters is that neither be adopted without noticing a choice was
made.

Note that dropping TV is safe because of a *gap in the software*, not
a property of natural history. Under paper 06's identity a biomarker
correlated with drift would bias the slope exactly as a PB-correlated
one does; there is simply no parameter to express it.

## Differences from paper 01

1. **Two components rather than three.** The covariance matrix falls
   from `(2 + 3n)` to `(2 + 2n)`: 26x26 to 18x18 at eight occasions.
   PB is the component removed; BR and TV are retained.
2. **Two arms rather than three.** The `orig` vendored-Hendrickson
   arm is dropped. Paper 13 compares `mean` against `covar` only,
   both through the package's own `generateData()`.
3. **Corrected carryover recursion.** The BR mean carryover is
   anchored to the mean at discontinuation rather than recursing on
   the already-adjusted previous value. See Section 2.4 of the
   manuscript and `02-deviations.md` in the driver directory.
4. **Numbers read from data.** Paper 01 hardcodes its results in
   prose. Paper 13's manuscript reads every quoted figure from the
   driver output at render time.

## Scope

1. Restate the two architectures under the reduced decomposition.
2. Cross architecture x design x N x `c_bm` x `t_half` (108 cells)
   and compare the carryover-response profile of each architecture.
3. Verify Type I error at the null under the reduced DGP.
4. Set out when the reduction is and is not appropriate.

## What this paper does not claim

That PB is dispensable. It remains necessary when the estimand is
the attribution of response to pharmacology versus belief (paper
06), when a biomarker may correlate with placebo responsiveness, and
for any study of blinding or expectancy.

**Important consequence.** The design's expectancy weight `e` enters
the DGP only through the PB component (its mean and its standard
deviation). With PB removed, `e` has no effect at all, so the
open-label versus blinded contrast is invisible to data generation
and the three designs differ only through their on-drug patterns.
This does not affect the carryover comparison, which acts through
`tod` and `tsd`, but it bounds what else the reduced DGP can be used
for.

## Driver

`analysis/scripts/two-component-dgp/01-run-grid.R`

```bash
# smoke test
Rscript analysis/scripts/two-component-dgp/01-run-grid.R --reps 3
# production (~75 min)
Rscript analysis/scripts/two-component-dgp/01-run-grid.R --reps 500
```

Writes `analysis/data/13-two-component-grid.rds`.

## Render

```bash
bash tools/render.sh analysis/report/13-two-component-dgp/report.Rmd
```

No simulation code runs during knitting.

## Package support

The reduction is implemented by the `components` argument added to
`generateData()` and `buildSigma()` in `R/generateData.R`. It
defaults to `c('tv','pb','br')`, so existing calls are unaffected.
Passing `c('tv','br')` produces the reduced matrix. The argument
requires that `br` be retained, since BR carries the interaction
channel in both architectures.

## Relationship to other papers

- `01-dgp-mean-moderation-vs-mvn`: the three-component original;
  paper 13 shares its question, notation, designs, and analysis
  model.
- `06-component-decomposition`: the methodological basis for the
  three-component decomposition and the source of the
  omitted-variable identity paper 13 relies on.
- `02-carryover-sensitivity`: shares the carryover parameterization
  and the corrected recursion.

## Author

pmsimstats team
