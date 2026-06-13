# How to run the Study A simulation

*2026-06-13*

## What it does

Runs the Study A grid from paper 06 (component decomposition):
80 alternative cells (m_PB x m_TV x N) and 4 null cells, each
at N_REPS replicates, fitting two analysis models
(one_component, phase_augmented) per replicate.

## Prerequisites

Host R must have these packages:

```r
install.packages(c('pkgload', 'data.table', 'nlme'))
# parallel ships with base R
```

`pkgload` is a lightweight substitute for `devtools::load_all()`.
If `devtools` is already installed, both work.

## Configure and run

1. Open `run-study-a-prod.R` and set the one parameter at the top:

   ```r
   N_REPS <- 150L   # change to any integer
   ```

2. From the repository root:

   ```
   cd /path/to/pmsimstats-ng
   Rscript analysis/scripts/component-decomposition/run-study-a-prod.R
   ```

   Or to run in the background and capture a log:

   ```
   nohup Rscript analysis/scripts/component-decomposition/run-study-a-prod.R \
     > run-study-a-N150.log 2>&1 &
   ```

## Where results land

All output goes under:

```
analysis/data/derived_data/component-decomposition/study-a-N<nreps>/
```

Each cell writes one file as it completes:

```
cell-001.rds  ...  cell-084.rds   # raw replicates, one per cell
study-a-all-reps.rds              # collated raw data (written at end)
study-a-summary.rds               # Morris ADEMP summary table
study-a-summary.txt               # same, tab-delimited
```

N_REPS is encoded in the directory name so runs at different rep
counts do not overwrite each other.

## Resuming an interrupted run

The script is resume-safe. If interrupted, re-run the same command.
Cells with an existing `cell-NNN.rds` are skipped; only incomplete
cells are re-run. The collation step at the end rebuilds
`study-a-all-reps.rds` and the summary from whatever checkpoints
are present.

## Running on a different machine

1. Clone or sync the repo to the new machine.
2. Install prerequisites (see above).
3. Set N_REPS, then run as above from the repo root.
4. Copy the output directory back:

   ```
   rsync -av \
     othermachine:/path/to/pmsimstats-ng/analysis/data/derived_data/component-decomposition/study-a-N1000/ \
     analysis/data/derived_data/component-decomposition/study-a-N1000/
   ```

5. If you want to split the cell grid across machines, the simplest
   approach is to run the full script on each machine and take the
   results from whichever finished first. Because each cell uses a
   deterministic seed (`STUDY_SEED + 1000 * cell_id + rep_idx`),
   results are reproducible and identical across machines given the
   same R version and package versions.

## Timing reference (8-core MacBook Pro M-series)

| N_REPS | Estimated wall time |
|--------|-------------------|
| 50     | ~22 min           |
| 150    | ~67 min           |
| 500    | ~3.7 hr           |
| 1000   | ~8.8 hr           |

Estimates scale linearly; 7 cores are used by default
(`detectCores() - 1`). N=100 and N=150 cells take 5-10x longer
per replicate than N=35 cells because they use the full 16-visit
design.

## Reading the summary

```r
library(data.table)
s <- readRDS(
  'analysis/data/derived_data/component-decomposition/
  study-a-N150/study-a-summary.rds')
# columns: analysis, regime, cell_id, m_PB, m_TV, N, c_bm,
#          n_reps, n_converged, conv_rate,
#          power, mcse_power,
#          mean_beta, bias, mcse_bias,
#          coverage
s[analysis == 'one_component' & regime == 'alt'][order(N, m_PB)]
```
