#' run-all.R
#'
#' Orchestrator for the paper 09 (informative-dropout-by-design)
#' production simulation programme. Sources each study script in its
#' own local() environment so state does not leak between studies.
#'
#' Order: smoke test, Study 1 (largest), Study 2 (two-mode bias),
#' Study 3 (path-conditional decomposition), Study 4 (analysis-side
#' misspecification, MI stub).
#'
#' Run from the repository root:
#'   Rscript analysis/scripts/informative-dropout-by-design/run-all.R

script_dir <- file.path('analysis', 'scripts',
                        'informative-dropout-by-design')

order <- c(
  '99-smoke-test.R',
  '01-study1.R',
  '02-study2.R',
  '03-study3.R',
  '04-study4.R'
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
cat(sprintf('\n== Paper 09 production programme: %.2f min total ==\n',
            as.numeric(total_elapsed)))
