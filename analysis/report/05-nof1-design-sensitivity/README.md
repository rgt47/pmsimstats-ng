# nof1-design-sensitivity
*2026-10-28 17:50 PDT*

**Author.** pmsimstats team

## Working title

Power of the Hybrid N-of-1 design under carryover, serial correlation
and variance components.

## Status

The manuscript (`report.Rmd`) was rewritten on 2026-10-28 after the
referee report of the same date (`referee-report-2026-10-28.md`). It
reads every Results number from the stored artifacts through
`results-helpers.R`.

2026-10-29: the full rerun (host job 51, commit 6485001, `nof1power`
0.1.0, 2 h 45 min) replaced the 2026-05 artifacts. Every sweep
artifact has schema-2 metadata, no `*-UNCHECKED.rds` exists, no fit
failed, and `pd_repaired` is `FALSE` throughout. The `git_rev` field
is empty because the host container had no git metadata. S13
(`13-cycle-period.R`) is not in `run-all.R`; its 2026-05 artifact is
stale and the manuscript does not use it. The rerun showed that the
alternating design and the crossover have a shifted mean under the
null (Type I 0.09 and 0.15). `18-null-shift-check.R` (run locally)
traces it to the missing term for the change from baseline to
follow-up in the analysis model. The fix is an engine change and has
not been made. See the manuscript's Open items.

## What changed from the 2026-05 version

- The Hybrid is the published schedule (weeks 4, 8, 9, 10, 11, 12,
  16, 20; four paths). The old two-path weekly schedule is kept as
  `Compact2` for S1 only.
- All designs share the visit calendar and the blinding assumption
  (expectancy 0.5, except the open-label visits of the Hybrid).
  `RCT-OL` keeps the old open-label assumption for comparison.
- A common estimand, `delta_on` (mean drug response over on-drug
  follow-up visits). The engine sets the drug-response asymptote per
  design with `target_delta`. Baseline `delta_on = -0.5`.
- The parallel group is analyzed by baseline-adjusted ANCOVA. The old
  mixed-model analysis is `parallel_analysis = 'legacy'`.
- The analysis half-life is crossed with the data-generating
  half-life in S2, S3, S11 and S12.
- S4 varies treatment heterogeneity `tau` (new in the DGP) and,
  separately, the baseline SD. S5 scales all three component SDs. S6
  fixes trial length and varies periods. S9 passes its whole grid.
- Cells with a non-positive-definite covariance matrix are dropped
  and listed in the artifact. The repair is never used.
  `pd_repaired` is stored.
- The power form of `carryover_decay()` no longer returns `NaN`.
- Every driver asserts that the stored parameter values equal its
  grid. A failed assertion saves `*-UNCHECKED.rds`.

## Directory contents

| File | Role |
|------|------|
| `report.Rmd` | Manuscript source |
| `results-helpers.R` | Reads artifacts, builds Results text and tables |
| `references.bib` | Cited entries only |
| `statistics-in-medicine.csl` | Citation style |
| `referee-report-2026-10-28.md` | Referee report |
| `title-page.md`, `cover-letter.md` | Submission material (drafts) |
| `bullets.Rmd`, `claudecode.tex` | Inert scaffolding of the retired three-part paragraph layout |

Drivers: `analysis/scripts/nof1-design-sensitivity/`. Engine:
`implementations/nof1power/` (version 0.1.0). Protocol:
`analysis/scripts/nof1-design-sensitivity/00-ademp-protocol.md`.

## Rerun commands

Run from the repository root. The engine must be installed from the
working tree, because the parallel workers load the installed
package.

```bash
Rscript -e "devtools::install_local(
  'implementations/nof1power', force = TRUE)"
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript -e "tinytest::run_test_dir('implementations/nof1power/inst/tinytest')"
Rscript analysis/scripts/nof1-design-sensitivity/00-smoke-test.R
```

Then either all sweeps in sequence:

```bash
Rscript analysis/scripts/nof1-design-sensitivity/run-all.R
```

or each driver (the order in `run-all.R` is S14, S1, S12, S12b, S2,
S6, S7, S4, S5, S8, S9, S3, S10, S11):

```bash
Rscript analysis/scripts/nof1-design-sensitivity/14-rct-null-check.R
Rscript analysis/scripts/nof1-design-sensitivity/17-feasibility.R
Rscript analysis/scripts/nof1-design-sensitivity/01-effect-size.R
Rscript analysis/scripts/nof1-design-sensitivity/02-carryover-halflife.R
Rscript analysis/scripts/nof1-design-sensitivity/03-decay-form.R
Rscript analysis/scripts/nof1-design-sensitivity/04-tau2.R
Rscript analysis/scripts/nof1-design-sensitivity/05-sigma2.R
Rscript analysis/scripts/nof1-design-sensitivity/06-cycles.R
Rscript analysis/scripts/nof1-design-sensitivity/07-patients.R
Rscript analysis/scripts/nof1-design-sensitivity/08-ar1.R
Rscript analysis/scripts/nof1-design-sensitivity/09-period.R
Rscript analysis/scripts/nof1-design-sensitivity/10-biomarker.R
Rscript analysis/scripts/nof1-design-sensitivity/11-misspec.R
Rscript analysis/scripts/nof1-design-sensitivity/12-morris-null.R
Rscript analysis/scripts/nof1-design-sensitivity/15-xo-corcar1.R
```

Finally:

```bash
Rscript analysis/scripts/nof1-design-sensitivity/16-summarize-all.R
Rscript analysis/scripts/nof1-design-sensitivity/plot-all.R
bash tools/render.sh analysis/report/05-nof1-design-sensitivity/report.Rmd
perl tools/notation-lint.pl analysis/report/05-nof1-design-sensitivity/report.Rmd
```

Long runs go to the simulation host. See
`analysis/scripts/quick-sim/carryover-closed-form/RERUN-anchored.md`
for the `simrun` procedure. Results from the container and from a
local session are not draw-for-draw identical, so run all sweeps on
one platform.

## Expected time

Timings were measured on the development laptop, single core,
200 replicates per cell and analysis, with the pooled large-sample fit
included: Hybrid 18 s, Alt 14 s, RCT 5 s, XO 7 s, Alt with 16 visits
32 s. The cost is linear in the replicates and in the number of
analyses. The host speed has not been measured. The host column
assumes 2.5 effective workers (2 physical cores, 4 threads) at 0.7
times the speed of the laptop core. Treat it as an estimate within a
factor of two.

| Sweep | Replicates | Cells | Laptop core-hours | Host hours (est.) |
|---|---|---|---|---|
| S1 | 1,000 | 42 | 0.63 | 0.36 |
| S2 | 500 | 36 | 0.56 | 0.32 |
| S3 | 500 | 15 | 0.40 | 0.23 |
| S4 | 500 | 30 | 0.26 | 0.15 |
| S5 | 500 | 15 | 0.13 | 0.08 |
| S6 | 500 | 5 | 0.10 | 0.06 |
| S7 | 500 | 18 | 0.25 | 0.14 |
| S8 | 500 | 15 (2 or 3 dropped) | 0.12 | 0.07 |
| S9 | 500 | 15 | 0.15 | 0.09 |
| S10 | 500 | 48 (6 dropped) | 0.37 | 0.21 |
| S11 | 500 | 24 | 0.64 | 0.36 |
| S12 | 2,000 | 16 | 1.00 | 0.57 |
| S12b | 2,000 | 2 | 0.12 | 0.07 |
| Total | | | 4.7 | 2.7 |

S14 (1,000 replicates, 24 cells) took 7 minutes on 7 workers locally
and has been run. Script 17 (feasibility table, no simulation) takes
seconds and has been run. Rendering takes a few seconds.

## Checks after the rerun

1. No artifact is named `*-UNCHECKED.rds`.
2. `pd_repaired` is `FALSE` in every artifact, and the dropped cells
   in S8 and S10 are listed in the Results.
3. In S1, the RCT and RCT-OL Type I error is within two MCSE of 0.05
   and $\kappa$ is within 5% of 1, as in the S14 check.
4. In S12, the pooled Type I error of every design is reported with
   its distance from 0.05 in MCSE.
5. The rendered PDF has no RESULTS-PENDING note and no `??`.

## Open decisions

- Units: the manuscript keeps days. NOTATION.md puts papers 01 to 04
  on weeks.
- The baseline effect `delta_on = -0.5` is a choice made so that the
  Hybrid is off power 1.000. Change `BASE_DELTA` in `00-designs.R` to
  move it.
- The author list in `title-page.md` and `cover-letter.md`.
- Small-sample corrections (Kenward-Roger, CR2) are named in the
  Limitations and not implemented.
