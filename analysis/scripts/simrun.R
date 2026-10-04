## simrun.R -- progress reporting for jobs offloaded with `simrun`.
##
## Copy into a project (simrun init does this) and source it from a
## simulation script. Both functions do nothing unless the script runs
## under simrun, so the same script runs unchanged on a workstation.
##
##   simrun_progress_init(n_reps, 'factorial grid')
##   res <- parallel::mclapply(seq_len(n_reps), function(i) {
##     out <- run_one(i)
##     simrun_tick()
##     out
##   }, mc.cores = as.integer(Sys.getenv('SIM_CORES', '1')))
##
## Each tick appends one line to the progress file. Appends of a single
## short line are atomic, so forked workers (mclapply, future multicore)
## can tick concurrently without coordination.

simrun_progress_init <- function(total, label = '') {
  f <- Sys.getenv('SIMRUN_PROGRESS_FILE')
  if (!nzchar(f)) return(invisible(FALSE))
  cat(sprintf('#total %d %s\n', as.integer(total), label),
      file = f, append = TRUE)
  invisible(TRUE)
}

simrun_tick <- function(n = 1L) {
  f <- Sys.getenv('SIMRUN_PROGRESS_FILE')
  if (!nzchar(f)) return(invisible(FALSE))
  cat(strrep('.\n', as.integer(n)), file = f, append = TRUE)
  invisible(TRUE)
}
