# 14-two-component-drop-tv

*2026-09-08 16:30 PDT*

## Working title

Data generation machinery for N-of-1 clinical trial simulations
under a two-component response decomposition.

## Origin

Paper 14 is the companion to paper 13. The two papers evaluate the
two candidate two-component reductions of the three-component DGP:

| paper | components | dropped |
|---|---|---|
| 13 | BR + TV | placebo-belief (PB) |
| 14 | BR + PB | natural history (TV) |

**This paper drops TV.** That reduction is *unconditionally* safe:
the implementation provides no biomarker-TV parameter, so the
deleted block cannot carry any part of the biomarker coupling under
any setting. It also preserves every feature by which the designs
differ, since the expectancy weight continues to act through PB and
the calendar-time argument continues to drive the AR(1) week gaps.

The cost is conventional rather than structural: a trend term is
standard in crossover and N-of-1 models, so this reduction deletes
the component that literature would most expect a model to keep.
Paper 13 makes the opposite trade. Section 4.5 of the manuscript
sets the two arguments against each other rather than declaring a
winner.

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
   TV is the component removed; BR and PB are retained.
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
biomarker may correlate with natural-history drift, and when
absolute power is reported as a design input rather than as a
comparison.

Note also that the safety of this reduction rests on a *gap in the
software* (no `c.bm.tv` parameter) rather than on a property of
natural history. Under paper 06's identity, a biomarker correlated
with drift would bias the slope exactly as a PB-correlated one does.

## Driver

`analysis/scripts/two-component-dgp/01-run-grid.R`

```bash
# smoke test
Rscript analysis/scripts/two-component-dgp/01-run-grid.R --drop tv --reps 3
# production (~75 min)
Rscript analysis/scripts/two-component-dgp/01-run-grid.R --drop tv --reps 500
```

Writes `analysis/data/14-drop-tv-grid.rds`.

## Render

```bash
bash tools/render.sh analysis/report/14-two-component-drop-tv/report.Rmd
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
