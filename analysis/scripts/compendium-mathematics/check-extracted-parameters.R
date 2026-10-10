## Copied 2026-10-08 from the review scratch script inv_c_check2.R;
## run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-extracted-parameters.R
# Supplementary checks for inventory part c: exact Gompertz parameters and the
# working-independence variants of docs/33 Section 2.6.
load('data/extracted_rp.rda')
obj <- get(ls()[1])
print(obj)

gomp <- function(t, m, d, r) {
  m * (exp(-d * exp(-r * t)) - exp(-d)) / (1 - exp(-d))
}
w <- c(4, 8, 9, 10, 11, 12, 16, 20)
cat(sprintf('TV/PB (6.51) at weeks: %s\n',
  paste(sprintf('%.3f', gomp(w, 6.51, 5, 0.35)), collapse = ' ')))
cat(sprintf('BR (10.99) at 4, 8, 9, 10: %s\n',
  paste(sprintf('%.3f', gomp(c(4, 8, 9, 10), 10.99, 5, 0.42)), collapse = ' ')))

tsd_of <- function(w, on) {
  out <- rep(NA_real_, length(w)); last_on <- NA
  for (k in seq_along(w)) {
    if (on[k] == 1) {
      out[k] <- 0; last_on <- w[k]
    } else if (!is.na(last_on)) out[k] <- w[k] - last_on
  }
  out
}
dbc <- function(w, on, th) {
  tsd <- tsd_of(w, on)
  g <- numeric(length(w)); g[on == 1] <- 1
  k <- on == 0 & !is.na(tsd)
  if (th > 0) g[k] <- 2^(-tsd[k] / th)
  g
}
designs <- list(
  CO = list(seq(2.5, 20, 2.5), list(c(1, 1, 1, 1, 0, 0, 0, 0),
    c(0, 0, 0, 0, 1, 1, 1, 1)), c(35, 35)),
  Hybrid = list(w, list(c(1, 1, 1, 1, 0, 0, 1, 0), c(1, 1, 1, 1, 0, 0, 0, 1),
    c(1, 1, 1, 0, 0, 0, 1, 0), c(1, 1, 1, 0, 0, 0, 0, 1)), c(18, 18, 17, 17)),
  'OL+BDC' = list(c(4, 8, 12, 16, 17, 18, 19, 20), list(c(1, 1, 1, 1, 1, 1, 0, 0),
    c(1, 1, 1, 1, 1, 0, 0, 0)), c(35, 35)))
variant <- function(dn, th, pooled, weighted, base) {
  d <- designs[[dn]]
  xs <- lapply(d[[2]], function(on) {
    x <- dbc(d[[1]], on, th); if (base) c(0, x) else x
  })
  wt <- if (weighted) d[[3]] else rep(1, length(xs))
  if (pooled) {
    allx <- unlist(xs); allw <- rep(wt, sapply(xs, length))
    m <- sum(allw * allx) / sum(allw)
    sum(allw * (allx - m)^2)
  } else {
    sum(wt * sapply(xs, function(x) sum((x - mean(x))^2)))
  }
}
for (pooled in c(FALSE, TRUE)) for (weighted in c(FALSE, TRUE)) for (base in c(FALSE, TRUE)) {
  out <- sapply(names(designs), function(dn)
    1 - variant(dn, 1, pooled, weighted, base) / variant(dn, 0, pooled, weighted, base))
  cat(sprintf('pooled %s weighted %s baseline %s: %s\n', pooled, weighted, base,
    paste(sprintf('%s %.1f%%', names(out), 100 * out), collapse = '; ')))
}
