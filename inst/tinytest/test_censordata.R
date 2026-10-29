## censordata(): the biased term must lower dropout with improvement
## (Hendrickson et al. 2020, Methods), and beta0 is a per-timepoint
## hazard. Before 2026-10-29 the term raised dropout with improvement.

make_censor_data <- function(N, seed = 20261029) {
  td <- buildtrialdesign(
    name_longform = 'open label',
    name_shortform = 'OL',
    timepoints = cumulative(rep(2.5, 8)),
    timeptnames = paste0('OL', 1:8),
    expectancies = rep(1, 8),
    ondrug = list(pathA = rep(1, 8))
  )
  tp <- td$trialpaths[[1]]
  labs <- tp$timeptnames
  set.seed(seed)
  dat <- data.table(ptID = seq_len(N), BL = 70)
  imp <- rnorm(N, 15, 12)
  for (i in seq_along(labs)) {
    dat[, (paste0('D_', labs[i])) := imp * i / length(labs)]
    dat[, (labs[i]) := BL - get(paste0('D_', labs[i]))]
  }
  list(dat = dat, tp = tp, labs = labs, imp = imp)
}

{
  s <- make_censor_data(N = 40000)
  set.seed(1)
  out <- censordata(copy(s$dat), s$tp,
                    list(beta0 = 0.05, beta1 = 0.999, eb1 = 4))
  dropped <- rowSums(is.na(as.matrix(out[, s$labs, with = FALSE]))) > 0
  q <- cut(s$imp, quantile(s$imp, c(0, 0.25, 0.75, 1)),
           include.lowest = TRUE, labels = c('low', 'mid', 'high'))
  rates <- tapply(dropped, q, mean)
  expect_true(rates[['low']] > rates[['high']] + 0.03,
              info = paste('biased dropout is higher among the least',
                           'improved than among the most improved'))
  expect_true(cor(dropped, s$imp) < 0,
              info = 'dropout is negatively correlated with improvement')
}

{
  s <- make_censor_data(N = 40000, seed = 7)
  set.seed(2)
  out <- censordata(copy(s$dat), s$tp,
                    list(beta0 = 0.05, beta1 = 0.5, eb1 = 2))
  miss <- is.na(as.matrix(out[, s$labs, with = FALSE]))
  first <- apply(miss, 1, function(r) if (any(r)) which(r)[1] else NA)
  new_drop <- tabulate(first[!is.na(first)], nbins = length(s$labs))
  at_risk <- nrow(miss) - c(0, cumsum(new_drop)[-length(s$labs)])
  ## The flat and biased draws can coincide, so the union is slightly
  ## below beta0 (about 0.049 here): approximately, not exactly, beta0.
  expect_equal(mean(new_drop / at_risk), 0.05, tolerance = 0.05,
               info = 'beta0 is approximately the per-timepoint hazard')
  expect_true(all(apply(miss, 1, function(r)
                    !any(r) || all(r[which(r)[1]:length(r)]))),
              info = 'dropout is monotone: once censored, always censored')
}
