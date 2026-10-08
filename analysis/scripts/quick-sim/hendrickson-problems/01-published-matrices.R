## Problem 1 evidence. Run Hendrickson's own generateData() at 58b32a9
## with the published parameters, capture the covariance matrix it
## passes to mvrnorm, and test it with the corpcor check the published
## code uses; then apply the published repair and measure the change.
##
## Requires 00-fetch-58b32a9.sh. Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/01-published-matrices.R
## Whitepaper: docs/36-hendrickson-pd-and-carryover.md
suppressMessages({library(data.table); library(corpcor)})
src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
out_dir <- 'analysis/data/quick-sim/hendrickson-problems'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

env <- new.env()
for (f in c('utilities.R', 'buildtrialdesign.R', 'generateData.R'))
  sys.source(file.path(src, 'R', f), envir = env)
## Intercept the draw: record Sigma and return the means. A non-PD
## Sigma cannot be sampled, so no real draw is attempted.
captured <- NULL
env$mvrnorm <- function(n, mu, Sigma, empirical = FALSE) {
  captured <<- Sigma
  matrix(rep(mu, each = n), nrow = n, dimnames = list(NULL, colnames(Sigma)))
}

e <- new.env(); load(file.path(src, 'data/results_core.rda'), envir = e)
ps <- e$results_core$parameterselections
rp <- as.data.table(ps$respparamsets[[1]]$param)
bp <- as.data.table(ps$blparamsets[[1]]$param)
tds <- ps$trialdesigns
dnames <- sapply(tds, function(td) td$metadata$name_shortform)

res <- list()
for (k in seq_along(tds)) {
  paths <- tds[[k]]$trialpaths
  for (g in seq_along(paths)) {
    for (i in seq_len(nrow(ps$modelparams))) {
      ## The published code indexes modelparam as a data.frame row.
      mp <- ps$modelparams[i, ]
      mp$N <- 5
      invisible(env$generateData(mp, rp, bp, paths[[g]], FALSE, FALSE))
      S <- captured
      R <- cov2cor(S)
      pd <- is.positive.definite(S)
      Rrep <- if (pd) R else cov2cor(make.positive.definite(S, tol = 1e-3))
      br <- grep('\\.br$', colnames(S))
      on <- paths[[g]]$tod > 0
      res[[length(res) + 1]] <- data.table(
        design = dnames[k], path = g, N = ps$modelparams$N[i],
        c.bm = mp$c.bm, t_half = mp$carryover_t1half, pd = pd,
        min_eig_cor = min(eigen(R, TRUE, TRUE)$values),
        r_on_spec = mean(R['bm', br][on]),
        r_on_used = mean(Rrep['bm', br][on]),
        r_off_spec = if (any(!on)) mean(R['bm', br][!on]) else NA_real_,
        r_off_used = if (any(!on)) mean(Rrep['bm', br][!on]) else NA_real_,
        max_cor_change = max(abs(Rrep - R)))
    }
  }
}
out <- rbindlist(res)
fwrite(out, file.path(out_dir, 'published-matrices.csv'))

cat('Failures of the published PD check, by design (design x path x 18 cells)\n')
print(out[, .(checks = .N, failures = sum(!pd)), by = design])
cat(sprintf('total: %d of %d\n\n', sum(!out$pd), nrow(out)))
cat('Failing cells (N omitted; it does not enter the covariance)\n')
print(unique(out[pd == FALSE, !'N'])[, lapply(.SD, function(x)
  if (is.numeric(x)) round(x, 3) else x)])
