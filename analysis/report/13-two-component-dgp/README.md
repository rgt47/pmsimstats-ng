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
process that carries only the pharmacological (BR) and
placebo-belief (PB) components. The time-variant natural-history
component (TV) is removed.

The motivation is that TV carries no biomarker coupling in any
architecture. The package provides `c.bm.pb` for a biomarker-PB
correlation but no equivalent for TV, so the biomarker-TV entry of
the covariance matrix is structurally zero in every simulation this
program has run. Paper 06's omitted-variable identity then implies
TV cannot displace the biomarker-treatment slope. Paper 13 tests
whether the architecture comparison survives its removal, and finds
that it does.

## Differences from paper 01

1. **Two components rather than three.** The covariance matrix falls
   from `(2 + 3n)` to `(2 + 2n)`: 26x26 to 18x18 at eight occasions.
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

That TV is dispensable. It remains necessary when the estimand is
the attribution of response to its causal sources (paper 06), when a
biomarker may correlate with natural-history drift (not currently
simulable), and when absolute power is reported as a design input
rather than as a comparison.

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
Passing `c('pb','br')` produces the reduced matrix. The argument
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
