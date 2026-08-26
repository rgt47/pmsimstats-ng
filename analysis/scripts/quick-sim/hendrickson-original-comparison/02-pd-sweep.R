suppressPackageStartupMessages({
  library(devtools); load_all('/Users/zenn/Library/CloudStorage/Dropbox/prj/res/36-pmsimstats-ng/pmsimstats-ng', quiet = TRUE)
  library(data.table); library(corpcor); library(MASS)
})

pd_log <- new.env(parent = emptyenv())
pd_log$rows <- list()

original_is_pd <- corpcor::is.positive.definite
instrumented_is_pd <- function(m, tol = 1e-8, method = c('eigen','chol')) {
  result <- original_is_pd(m, tol = tol)
  pd_log$rows[[length(pd_log$rows) + 1]] <- data.frame(pd = result)
  result
}
environment(instrumented_is_pd) <- asNamespace('corpcor')
assignInNamespace('is.positive.definite', instrumented_is_pd, ns = 'corpcor')
assign('is.positive.definite', instrumented_is_pd, envir = .GlobalEnv)

orig_env <- new.env(parent = globalenv())
sys.source(file.path(Sys.getenv('SCRATCH'), 'orig_generateData.R'), envir = orig_env)
environment(orig_env$generateData) <- orig_env

data(extracted_bp); data(extracted_rp)
bp <- extracted_bp; rp <- extracted_rp
td_OLBDC <- buildtrialdesign(
  name_shortform = 'OL+BDC', timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptnames = c('OL1','OL2','OL3','OL4','BD1','BD2','BD3','BD4'),
  expectancies = c(1,1,1,1,.5,.5,.5,.5),
  ondrug = list(pathA = c(1,1,1,1,1,1,0,0), pathB = c(1,1,1,1,1,0,0,0))
)
path1 <- td_OLBDC$trialpaths[[1]]

cbm_grid <- seq(0.1, 0.8, by = 0.05)
t1half <- 1.0
out <- data.frame()
for (cbm in cbm_grid) {
  pd_log$rows <- list()
  mp_orig <- data.table(N = 35L, c.bm = cbm, carryover_t1half = t1half,
                         c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.1)
  invisible(orig_env$generateData(mp_orig, rp, bp, data.table::copy(path1), FALSE, TRUE))
  pd_hendrickson <- pd_log$rows[[1]]$pd

  pd_log$rows <- list()
  mp_cur <- list(N = 35L, c.bm = cbm, carryover_t1half = t1half, c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.1)
  invisible(buildSigma(modelparam = mp_cur, respparam = rp, blparam = bp,
                        trialdesign = path1, makePositiveDefinite = TRUE,
                        lambda_cor = log(2)/t1half, dgp_architecture = 'mvn'))
  pd_archB <- pd_log$rows[[1]]$pd

  out <- rbind(out, data.frame(c.bm = cbm, pd_hendrickson_cs = pd_hendrickson, pd_archB_ar1 = pd_archB))
}
print(out, row.names = FALSE)
