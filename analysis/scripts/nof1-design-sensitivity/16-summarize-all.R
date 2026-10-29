#' 16-summarize-all.R
#'
#' Write one table of performance measures for every cell of every
#' schema-2 artifact: power or Type I error, bias against the
#' large-sample value and against the design-independent estimand,
#' empirical SE, model SE, kappa, coverage, failed fits, and their
#' Monte Carlo standard errors. Output:
#' analysis/data/derived_data/sensitivity/cell-summaries.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(here)
})
source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-summarize.R'
))

data_dir <- Sys.getenv(
  'SENS_OUT_DIR',
  here::here('analysis', 'data', 'derived_data', 'sensitivity')
)
files <- list.files(data_dir, pattern = '^sensitivity-\\d+-.*\\.rds$',
                    full.names = TRUE)
files <- files[!grepl('UNCHECKED', files)]

keys <- c('axis', 'delta_label', 'config', 'design_name',
          'dgp_architecture', 'trialdesign', 'respparamset',
          'blparamset', 'modelparamset', 'analysis_row',
          'analysis_t1half', 'delta_on', 'N', 'c.bm', 'c.br', 'tau',
          'br_sd', 'bl_sd', 'carryover_t1half', 'carryover_form',
          'weibull_shape', 'n_visits', 'n_switch', 'min_eig',
          'pd_repaired')

out <- lapply(files, function(f) {
  obj <- readRDS(f)
  if (!identical(obj$run_meta$schema, 2L)) return(NULL)
  r <- obj$results
  by <- intersect(keys, names(r))
  smry <- summarize_cells(r, by)
  smry$sweep <- obj$run_meta$sweep_id
  smry$n_reps_planned <- obj$run_meta$n_reps
  smry
})
out <- bind_rows(out)
dest <- file.path(data_dir, 'cell-summaries.csv')
utils::write.csv(out, dest, row.names = FALSE)
message('Wrote ', nrow(out), ' cells to ', dest)
