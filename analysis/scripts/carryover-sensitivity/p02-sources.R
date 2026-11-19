## analysis/scripts/carryover-sensitivity/p02-sources.R
##
## Data folders of paper 02 by parameter set. Paper 02 and its derived
## scripts read every stored run through these helpers, so one switch
## selects the set:
##
##   PAPER02_PARAMS=simcore2026  (default) the stored runs under
##                               analysis/data/carryover-sensitivity/
##   PAPER02_PARAMS=published    the runs under Hendrickson et al.'s
##                               values, in .../carryover-sensitivity/
##                               published/ (docs/74, wave 1)
##
## The helpers fail closed. A published-set run that is missing, or
## whose PARAM_SET marker does not read 'published', stops with an
## error; nothing falls back to the other set. Derived outputs (ladder
## tables, reanalysis files, figures) take a '-published' suffix under
## the published set, so the two sets never overwrite each other.
##
## Under the published set A4 (E7null_cr2) was fitted in the runs
## themselves, so the G10 backfill (script 41) does not apply.

p02_set <- function() {
  s <- Sys.getenv('PAPER02_PARAMS', 'simcore2026')
  if (!s %in% c('simcore2026', 'published'))
    stop("PAPER02_PARAMS must be 'simcore2026' or 'published', not '", s, "'")
  s
}
p02_published <- function() p02_set() == 'published'
p02_cs <- function() here::here('analysis/data/carryover-sensitivity')

## A published-set run folder, checked against its marker.
p02_check_marker <- function(d) {
  m <- file.path(d, 'PARAM_SET')
  if (!dir.exists(d)) stop('no published-set run at ', d)
  if (!file.exists(m) || !identical(readLines(m, n = 1L), 'published'))
    stop('folder lacks a published PARAM_SET marker: ', d)
  d
}

## A stored simulation run, by its folder name under the set-1 layout
## (for example 'tier1-g-vs-f-e0', 'referee-variants/block-N').
p02_run <- function(name) {
  if (p02_published()) p02_check_marker(file.path(p02_cs(), 'published', name))
  else file.path(p02_cs(), name)
}

## A derived output folder or file name under analysis/data/
## carryover-sensitivity/ (for example 'ladder-predictive').
p02_derived <- function(name) {
  file.path(p02_cs(), if (p02_published()) paste0(name, '-published') else name)
}

## A figure base name (no extension).
p02_fig <- function(base) if (p02_published()) paste0(base, '-published') else base

## Whether A4 comes from the backfill (set 1) or from the runs (published).
p02_backfill <- function() !p02_published()

## The supplementary (gap-fill) study. Under set 1 it is jobs 31, 32, 33
## and 40, each with all.rds and a combined summary.rds; under the
## published set it is the parts B, C, dropout and size. `job` holds the
## job number or part name, so cell keys stay unique within a set.
p02_gapfill_parts <- function() {
  root <- file.path(p02_cs(), if (p02_published()) 'published/g9-gap-fill' else 'g9-gap-fill')
  if (p02_published()) {
    parts <- c('B', 'C', 'dropout', 'size')
    for (p in parts) p02_check_marker(file.path(root, p))
    setNames(file.path(root, parts), parts)
  } else {
    jobs <- c('31', '32', '33', '40')
    setNames(file.path(root, paste0('job-', jobs)), jobs)
  }
}
p02_gapfill_summary <- function() {
  if (!p02_published())
    return(readRDS(file.path(p02_cs(), 'g9-gap-fill', 'summary.rds')))
  parts <- p02_gapfill_parts()
  dplyr::bind_rows(lapply(names(parts), function(p) {
    s <- readRDS(file.path(parts[[p]], 'summary.rds'))
    cbind(data.frame(job = p, host = NA_character_), s)
  }))
}
p02_gapfill_all <- function() {
  parts <- p02_gapfill_parts()
  dplyr::bind_rows(lapply(names(parts), function(p) {
    a <- readRDS(file.path(parts[[p]], 'all.rds'))
    a$job <- p
    a
  }))
}
