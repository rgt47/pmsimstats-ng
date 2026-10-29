#' 17. Feasibility of the S8 and S10 covariance matrices.
#'
#' No simulation. For every design and every cell of the S8 (AR(1)
#' correlation) and S10 (biomarker correlation by half-life, correlation
#' encoding) grids, builds the covariance matrix of every path with the
#' engine, without repair, and records the smallest eigenvalue of the
#' correlation matrix of the components with positive variance. A
#' negative value means the cell is infeasible and the sweep drops it.
#' The mean-shift encoding adds no correlation and is feasible whenever
#' the baseline is. Output: feasibility-s8-s10.rds in the artifact
#' directory.

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-common.R'
))

designs <- list(
  Hybrid = design_hybrid(),
  Alt = design_nof1_alternating(),
  RCT = design_rct()
)
rp <- baseline_resp_param()
bp <- baseline_bl_param()

min_eig <- function(td, mp) {
  vals <- vapply(td$trialpaths, function(path) {
    sg <- build_sigma_matrix(mp, rp, bp, path, stop_on_nonpd = FALSE)
    nof1power:::.min_eigenvalue(sg)
  }, numeric(1))
  min(vals)
}

s8 <- expand.grid(design = names(designs), rho = c(0, 0.3, 0.5, 0.7, 0.9),
                  stringsAsFactors = FALSE)
s8$min_eig <- mapply(function(d, rho) {
  min_eig(designs[[d]], baseline_model_param(
    N = BASE_N, c.tv = rho, c.pb = rho, c.br = rho))
}, s8$design, s8$rho)
s8$sweep <- 'S8'

s10 <- expand.grid(design = names(designs), c_bm = c(0, 0.2, 0.4, 0.6),
                   halflife_days = c(0, 3), stringsAsFactors = FALSE)
s10$min_eig <- mapply(function(d, cb, hd) {
  min_eig(designs[[d]], baseline_model_param(
    N = BASE_N, c.bm = cb, carryover_t1half = hd / 7))
}, s10$design, s10$c_bm, s10$halflife_days)
s10$sweep <- 'S10'

out <- dplyr::bind_rows(s8, s10)
out$feasible <- out$min_eig > 0
res <- list(
  table = out,
  run_meta = list(
    schema = 2L, script = '17-feasibility.R',
    engine_version = as.character(utils::packageVersion('nof1power')),
    date = format(Sys.time(), '%Y-%m-%d %H:%M %Z'))
)
saveRDS(res, file.path(SENS_DATA_DIR, 'feasibility-s8-s10.rds'))
print(as.data.frame(out[!out$feasible | out$min_eig < 0.05, ]),
      digits = 3)
