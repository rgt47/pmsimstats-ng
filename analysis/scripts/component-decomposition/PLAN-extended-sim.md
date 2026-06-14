# Plan: Study A extended simulation (pre-registered spec)

*2026-06-13 17:43 PDT*

## Objective

Run Study A at the pre-registered rep counts:

- **80 alternative cells** (c_bm = 0.45): 1,000 reps each
- **4 null cells** (c_bm = 0): 5,000 reps each

The 250-rep pilot (completed 2026-06-13) confirmed infrastructure
and established power estimates, but has insufficient precision to
demonstrate the load-bearing bias claim (signal/noise ~1.6 at
m_PB=10, N=70) or meet the pre-registered Type I error target
(MCSE 0.003 requires 5,000 null reps).

## Estimated wall time

Based on the 250-rep timing run (0.317 sec/rep on 7 cores):

| Segment          | Cells | Reps  | Est. time    |
|------------------|-------|-------|--------------|
| Alt N=35 (x20)   | 20    | 1,000 | ~0.3 hr      |
| Alt N=70 (x20)   | 20    | 1,000 | ~0.5 hr      |
| Alt N=100 (x20)  | 20    | 1,000 | ~2.6 hr      |
| Alt N=150 (x20)  | 20    | 1,000 | ~3.8 hr      |
| Null (x4)        | 4     | 5,000 | ~1.8 hr      |
| **Total**        | **84**|       | **~9.0 hr**  |

Scales inversely with core count. On a machine with 14 cores
(e.g. M3 Pro), expected time ~5.5 hr.

## Prerequisites on the other laptop

R must be installed (4.4+). Install three packages:

```r
install.packages(c('pkgload', 'data.table', 'nlme'),
                 repos = 'https://cloud.r-project.org')
# parallel is a base package -- no install needed
```

Verify with:

```r
library(pkgload); library(data.table); library(nlme)
```

## Sync the repo

The repo syncs automatically via Dropbox. Confirm the script is
present on the other machine at:

```
~/Library/CloudStorage/Dropbox/prj/res/36-pmsimstats-ng/
  pmsimstats-ng/analysis/scripts/component-decomposition/
  run-study-a-prod.R
```

If Dropbox is not available, clone via git:

```bash
git clone <repo-url> pmsimstats-ng
cd pmsimstats-ng
git checkout main
```

## Run the simulation

From the repo root, in a terminal that will stay open (or via
`nohup` to survive logout):

```bash
cd ~/Library/CloudStorage/Dropbox/prj/res/36-pmsimstats-ng/pmsimstats-ng

nohup Rscript \
  analysis/scripts/component-decomposition/run-study-a-prod.R \
  > /tmp/study-a-prod.log 2>&1 &

echo "PID: $!"
```

Monitor progress:

```bash
tail -f /tmp/study-a-prod.log
```

The script prints one line per cell on start and completion:

```
[   0.0s] cell 01 (m_PB= 0 m_TV=-1 N= 35) start
[  42.3s] cell 01 (m_PB= 0 m_TV=-1 N= 35) done (1000/1000 ok)
```

## Where results land

All output goes to:

```
analysis/data/derived_data/component-decomposition/
  study-a-alt1000-null5000/
    cell-001.rds  ...  cell-084.rds   # per-cell raw replicates
    study-a-all-reps.rds              # collated (written at end)
    study-a-summary.rds               # Morris ADEMP summary
    study-a-summary.txt               # tab-delimited summary
```

If the run is via Dropbox, results sync back automatically as
each cell completes (one checkpoint file per cell). If via git,
copy the output directory back manually:

```bash
rsync -av \
  /path/on/other-machine/analysis/data/derived_data/component-decomposition/study-a-alt1000-null5000/ \
  analysis/data/derived_data/component-decomposition/study-a-alt1000-null5000/
```

## Resume after interruption

The script is resume-safe. If it is interrupted (power, sleep,
crash), re-run the same command. Cells with an existing checkpoint
file are skipped; only incomplete cells re-run. The collation step
at the end rebuilds the summary from all present checkpoints.

## After the run completes

Update the data path in `report.Rmd`:

```r
# In the load-sim-data chunk, change:
sim_summary <- readRDS(here::here(
  'analysis/data/derived_data/component-decomposition',
  'study-a-N250/study-a-summary.rds'))

# To:
sim_summary <- readRDS(here::here(
  'analysis/data/derived_data/component-decomposition',
  'study-a-alt1000-null5000/study-a-summary.rds'))
```

Then re-render:

```bash
cd analysis/report/06-component-decomposition
Rscript -e "rmarkdown::render('report.Rmd')"
```
