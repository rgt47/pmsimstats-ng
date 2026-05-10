#' run-all.R
#'
#' Orchestrator for the paper 06 (component-decomposition)
#' production simulation programme. Sources each study script in its
#' own local() environment so state does not leak between studies.
#'
#' Order: smoke test, Study A (largest), Study B (full nonlinear,
#' currently stubbed), Study C (subadditivity).
#'
#' Run from the repository root, e.g.:
#'   Rscript analysis/scripts/component-decomposition/run-all.R
#' or, in a Docker container via `make r`:
#'   Rscript analysis/scripts/component-decomposition/run-all.R

script_dir <- file.path('analysis', 'scripts',
                        'component-decomposition')

order <- c(
  '06-smoke-test.R',
  '02-study-A.R',
  '03-study-B.R',
  '04-study-C.R'
)

total_start <- Sys.time()

for (script in order) {
  path <- file.path(script_dir, script)
  if (!file.exists(path)) {
    warning('Missing script: ', path)
    next
  }
  cat(sprintf('\n=== %s ===\nStart: %s\n',
              script, format(Sys.time(), '%Y-%m-%d %H:%M:%S')))
  t0 <- Sys.time()
  local({ source(path, local = TRUE) })
  elapsed <- difftime(Sys.time(), t0, units = 'mins')
  cat(sprintf('Done: %.2f min\n', as.numeric(elapsed)))
}

total_elapsed <- difftime(Sys.time(), total_start, units = 'mins')
cat(sprintf('\n== Paper 06 production programme: %.2f min total ==\n',
            as.numeric(total_elapsed)))
