## analysis/scripts/carryover-sensitivity/33-assemble-gap-fill.R
##
## Assemble the gap-fill study (32-run-g9-gap-fill.R, protocol
## 32-gap-fill-protocol.md) from its simrun jobs into one summary.
## The three blocks ran as separate jobs on two hosts, and block A1
## ran twice (per-occasion and patient-level dropout conventions), so
## each job's summary.rds lands in its own fetched directory under
## analysis/remote_runs/. This script copies the small outputs of
## each job (summary, preflight, meta, all) into
## analysis/data/carryover-sensitivity/g9-gap-fill/job-<id>/ and
## writes the combined summary.rds at the top level with `job` and
## `host` columns. Checkpoints stay in remote_runs (gitignored).
##
## Usage: Rscript analysis/scripts/carryover-sensitivity/33-assemble-gap-fill.R
suppressMessages({library(dplyr); library(purrr)})

jobs <- tribble(
  ~job, ~host,      ~dir,
  31L,  'zeus',     'analysis/remote_runs/job-31-ca2eef2',
  32L,  'zeus',     'analysis/remote_runs/job-32-ca2eef2',
  40L,  'mercury',  'analysis/remote_runs/job-40-64f1e70-simhost',
  33L,  'zeus',     'analysis/remote_runs/job-33-fc7ee08'
)
rel <- 'results/analysis/data/carryover-sensitivity/g9-gap-fill'
out_dir <- 'analysis/data/carryover-sensitivity/g9-gap-fill'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

present <- jobs |> mutate(src = file.path(dir, rel)) |>
  filter(file.exists(file.path(src, 'summary.rds')))
missing <- setdiff(jobs$job, present$job)
if (length(missing) > 0)
  message('Not yet fetched: job ', paste(missing, collapse = ', '))

combined <- pmap_dfr(present, function(job, host, dir, src) {
  dest <- file.path(out_dir, sprintf('job-%02d', job))
  dir.create(dest, showWarnings = FALSE)
  for (f in c('summary.rds', 'preflight.rds', 'meta.rds', 'all.rds'))
    file.copy(file.path(src, f), file.path(dest, f), overwrite = TRUE)
  readRDS(file.path(src, 'summary.rds')) |>
    mutate(job = job, host = host, .before = 1)
})

## The two A1 runs are told apart by dropout_per (job 31 predates the
## column and is per-occasion).
if (!'dropout_per' %in% names(combined)) combined$dropout_per <- NA_character_
combined <- combined |>
  mutate(dropout_per = if_else(is.na(dropout_per) & dropout_rate > 0,
                               'occasion', dropout_per),
         dropout_per = if_else(dropout_rate == 0, 'none', dropout_per))
stopifnot(!anyDuplicated(combined[c('cell_id', 'method')]))

saveRDS(combined, file.path(out_dir, 'summary.rds'))
cat(sprintf('%d rows from jobs %s -> %s/summary.rds\n',
            nrow(combined), paste(present$job, collapse = ', '), out_dir))
print(count(combined, job, host, block, sub))
