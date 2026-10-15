## Build decision-script inputs for job 37 (11 --extra4) alone and pooled
## with jobs 15 and 21, in the format 12-strawman-stress-decision.R reads.
## Job 37 rep ids are offset by 10,000 so pooled trials stay distinct.
## The inputs are large and are written to a work directory outside the
## repo (argument 1, or the EXTRA4_WORK environment variable); see
## 15-extra4-run.sh. Post hoc (docs/47, 'Fresh trials for the F1, G1, G2
## choice').
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
work <- if (length(args) >= 1) args[1] else Sys.getenv('EXTRA4_WORK')
stopifnot(nzchar(work))
d <- 'analysis/data/quick-sim/carryover-closed-form/strawman-stress'
x4 <- readRDS(file.path(d, 'replicates-full-extra4-anchored.rds'))
m <- readRDS(file.path(d, 'replicates-full-anchored.rds'))
e1 <- readRDS(file.path(d, 'replicates-full-extra-anchored.rds'))
stopifnot(setequal(unique(x4$analysis), c('S', 'C1', 'C2', 'ES', 'F1', 'F2', 'G1', 'G2')),
          uniqueN(x4$cell) == 41, all(x4[, .N, by = .(cell, analysis)]$N == 1500))
cat('job 37 failed fits:', sum(is.na(x4$p)), 'of', nrow(x4), '\n')
x4[, rep := rep + 10000L]
keep8 <- c('S', 'C1', 'C2', 'ES', 'F1', 'F2')
write_set <- function(dir, main, extra) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(main, file.path(dir, 'replicates-full-x.rds'))
  saveRDS(extra, file.path(dir, 'replicates-full-extra-x.rds'))
  saveRDS(extra[0], file.path(dir, 'replicates-full-extra2-x.rds'))
}
write_set(file.path(work, 'new'), x4[analysis %in% keep8], x4[analysis %in% c('G1', 'G2')])
pooled_main <- rbindlist(list(m[analysis %in% keep8], x4[analysis %in% keep8]), fill = TRUE)
pooled_ex <- rbindlist(list(e1[analysis %in% c('G1', 'G2')], x4[analysis %in% c('G1', 'G2')]),
                       fill = TRUE)
write_set(file.path(work, 'pooled'), pooled_main, pooled_ex)
fwrite(rbindlist(list(pooled_main, pooled_ex), fill = TRUE)[, .(cell, label, cbm, hyp, rep, analysis, p)],
       file.path(work, 'pooled-long.csv'))
fwrite(x4[, .(cell, label, cbm, hyp, rep, analysis, p)], file.path(work, 'new-long.csv'))
