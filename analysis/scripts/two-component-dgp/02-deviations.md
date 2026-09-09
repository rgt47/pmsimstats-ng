# Paper 13 driver: deviations and implementation notes

*2026-09-08 16:30 PDT*

## Deviation 1: corrected carryover recursion

The BR mean carryover in `R/generateData.R` and
`implementations/tidyverse/R/functions.R` previously computed

```r
brmeans[p] <- brmeans[p] +
  brmeans[p-1] * (1/2)^(tsd[p] / carryover_t1half)
```

`brmeans[p-1]` is the already-adjusted value from the preceding
occasion, while `tsd[p]` is cumulative time since discontinuation.
Elapsed time is therefore counted twice, and every off-drug occasion
after the first in an uninterrupted run is over-decayed. At a
one-week half-life the second off-drug occasion received a quarter
of the pre-discontinuation effect where it should receive a half.

Both files now anchor the decay to the mean at discontinuation:

```r
last_on <- 0
for (p in seq_len(nP)) {
  if (onDrug[p]) {
    last_on <- brmeans[p]
  } else if (tsd[p] > 0) {
    brmeans[p] <- brmeans[p] + last_on * (1/2)^(tsd[p] / t_half)
  }
}
```

**Effect on published results.** Checked at three reference cells
(Hybrid, OL+BDC, CO; N = 70, `c_bm` = 0.45, `t_half` = 1.0,
exponential, 400 replicates, common seed). Largest power movement
0.005 against an MCSE of about 0.025; point estimates agreed to two
or three decimals. Paper 02's published numbers stand. Paper 06 was
not checked and is more exposed, since its estimands include
component means directly.

**Not yet corrected.** `implementations/original-extended/` and
`implementations/nof1power/` carry the same recursion.
`implementations/original/` carries it and should keep it, since its
role is back-compatibility testing against the historical reference.
Correcting `original-extended` will break the parity tests in
`analysis/scripts/parity/` until the baselines are regenerated.

## Deviation 2: the components argument

`generateData()` and `buildSigma()` gained a `components` argument,
defaulting to `c('tv','pb','br')`. Existing calls are unaffected.

Validation: the argument must be a subset of `c('tv','pb','br')` and
must contain `'br'`, since BR carries the interaction channel under
both architectures. The supplied set is reordered to the canonical
`tv, pb, br` sequence so that label and standard-deviation ordering
is stable regardless of how the caller writes it.

Downstream assembly of the outcome columns already iterated over the
component list generically, so no further change was required.

The package test suite (70 assertions across six files) passes
unchanged with the default argument.

## Deviation 3: no matched three-component run

Section 3.3 of the manuscript compares absolute power against paper
01's published figures rather than against a three-component grid run
under identical conditions. Paper 01's figures were produced under
the earlier carryover recursion and a partly different grid, so the
comparison is qualitative only. A precise decomposition of the power
difference into its variance and recursion parts would require a
matched run and is not available.

## Grid

108 cells: 2 architectures x 3 designs x 2 sample sizes x 3
biomarker levels x 3 carryover half-lives, 500 replicates each.
Master seed 20260908, per-cell seed derived as `SEED + i`.

`N` is the total number of participants, allocated as evenly as
possible across randomization paths.

## Known limitation of the seeding scheme

The two architectures are run in separate replicate loops from
different cell seeds, so their replicates are not paired. This is
unavoidable, since the DGP itself differs by architecture and the
same random draw does not correspond to the same trial. Comparisons
between architectures at a common nominal `c_bm` are therefore
unpaired, and in any case not interpretable on an absolute scale
(the parameter is not calibrated to equal power across channels).
The interpretable comparison is the carryover profile within each
architecture.
