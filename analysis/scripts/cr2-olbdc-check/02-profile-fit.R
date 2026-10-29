## Time the pieces of one fit_analysis() call (lme vs CR2).
suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table)
  loadNamespace('clubSandwich')
})
Sys.setenv(CHK_REPS = '0')
env09 <- new.env(parent = globalenv())
exprs <- parse('analysis/scripts/informative-dropout-by-design/20-paired-dropout.R',
               keep.source = FALSE)
stop_at <- which(vapply(exprs, function(e) {
  is.call(e) && identical(e[[1]], as.name('<-')) &&
    identical(e[[2]], as.name('units'))
}, logical(1)))[1]
for (e in exprs[seq_len(stop_at - 1L)]) eval(e, envir = env09)

set.seed(1)
td <- env09$designs$OLBDC
paths <- td$trialpaths
mp <- env09$make_mp(0)
templates <- lapply(paths, env09$path_long_template)
yl <- bl <- bml <- fdl <- list()
for (g in 1:2) {
  mp$N <- 35
  di <- generateData(mp, env09$rp, env09$bp, paths[[g]], empirical = FALSE,
                     makePositiveDefinite = FALSE,
                     cached_sigma = env09$sigma_cache$OLBDC[[g]],
                     dgp_architecture = 'mean_moderation', stopOnNonPD = TRUE)
  tn <- paths[[g]]$timeptnames
  yl[[g]] <- as.matrix(di[, ..tn]); bl[[g]] <- di$BL; bml[[g]] <- di$bm
  fdl[[g]] <- rep(NA_integer_, 35)
}
long <- env09$build_long(yl, bl, bml, fdl, templates)
long[, Dbc := as.numeric(Db)]
print(system.time(fit <- nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc + De,
  random = ~1 | ptID, correlation = nlme::corCAR1(form = ~t | ptID),
  data = long, control = nlme::lmeControl(opt = 'optim', maxIter = 200,
                                          msMaxIter = 200))))
print(system.time(ct <- clubSandwich::coef_test(fit, vcov = 'CR2',
                                                test = 'Satterthwaite')))
print(system.time(v <- clubSandwich::vcovCR(fit, type = 'CR2')))
print(system.time(ct2 <- clubSandwich::coef_test(fit, vcov = v,
                                                 test = 'Satterthwaite',
                                                 coefs = 'bm:Dbc')))
print(ct); print(ct2)
print(packageVersion('clubSandwich')); print(packageVersion('nlme'))
print(sessionInfo()$BLAS)
