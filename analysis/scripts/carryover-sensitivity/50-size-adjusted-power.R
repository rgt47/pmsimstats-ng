## analysis/scripts/carryover-sensitivity/50-size-adjusted-power.R
##
## Size-adjusted power of A0 to A4 on the 30-cell grid (job 27; A4 from
## the backfill refits of the same trials), so that the liberal A0 can
## be compared with the correctly sized analyses (referee report
## 2026-11-14, M3). For each analysis, design and N the critical value
## is the empirical 5 percent quantile of the p-values in the null cell
## (c_bm = 0, t1/2 = 1 week); size-adjusted power is the share of
## alternative replicates (t1/2 = 1 week, c_bm in {0.30, 0.45}) with a
## p-value at or below it. Null cells exist only at t1/2 = 1 week, so
## only those alternatives are adjusted. A bootstrap over the null and
## alternative replicates (2,000 resamples) gives the standard error,
## which includes the error of the estimated critical value. Nothing is
## re-simulated.
##
## Run from the repository root:
##   Rscript analysis/scripts/carryover-sensitivity/50-size-adjusted-power.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

root <- here::here()
cs <- file.path(root, 'analysis/data/carryover-sensitivity')
out_dir <- file.path(cs, 'referee-reanalysis')

specs <- c('E0', 'E1', 'E1cr2', 'E2cr2', 'E7null_cr2')
a27 <- readRDS(file.path(cs, 'tier1-g-vs-f-e0/all.rds')) |>
  filter(method %in% specs, t1half == 1) |>
  select(cell_id, rep, method, p_value, converged, design, N, c_bm)
key27 <- a27 |> distinct(cell_id, design, N, c_bm)
g27 <- readRDS(file.path(cs, 'g10-backfill/all.rds')) |>
  filter(source == 't1', method == 'E7null_cr2') |>
  select(cell_id, rep, method, p_value, converged) |>
  inner_join(key27, by = 'cell_id')
d <- bind_rows(a27, g27) |>
  filter(converged, !is.na(p_value))
stopifnot(setequal(unique(d$method), specs))

adj_one <- function(p0, p1, B = 2000) {
  q <- unname(quantile(p0, 0.05, type = 7))
  est <- mean(p1 <= q)
  boot <- replicate(B, {
    qb <- unname(quantile(sample(p0, replace = TRUE), 0.05, type = 7))
    mean(sample(p1, replace = TRUE) <= qb)
  })
  c(crit_p = q, size = mean(p0 < 0.05), power = mean(p1 < 0.05),
    power_adj = est, power_adj_se = sd(boot))
}

set.seed(20261114)
res <- d |>
  distinct(design, N, method) |>
  arrange(design, N, method)
out <- list()
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  p0 <- d$p_value[d$design == r$design & d$N == r$N &
                    d$method == r$method & d$c_bm == 0]
  stopifnot(length(p0) >= 450)
  for (cb in c(0.30, 0.45)) {
    p1 <- d$p_value[d$design == r$design & d$N == r$N &
                      d$method == r$method & abs(d$c_bm - cb) < 1e-9]
    stopifnot(length(p1) >= 450)
    out[[length(out) + 1]] <- tibble(design = r$design, N = r$N,
                                     method = r$method, c_bm = cb,
                                     n_null = length(p0), n_alt = length(p1),
                                     !!!as.list(adj_one(p0, p1)))
  }
}
tab <- bind_rows(out) |>
  mutate(analysis = c(E0 = 'A0', E1 = 'A1', E1cr2 = 'A2', E2cr2 = 'A3',
                      E7null_cr2 = 'A4')[method])
print(as.data.frame(tab |> filter(design == 'Hybrid') |>
  select(N, c_bm, analysis, size, power, power_adj, power_adj_se)),
  digits = 3)
write.csv(tab, file.path(out_dir, 'size-adjusted-power-f1.csv'),
          row.names = FALSE)
cat('wrote', file.path(out_dir, 'size-adjusted-power-f1.csv'), '\n')

## Paired differences: A4 against A0 and A2 on the same trials. The
## null and alternative replicates are resampled jointly across the
## analyses, so the bootstrap SE reflects the pairing.
wide <- d |>
  select(design, N, c_bm, rep, method, p_value) |>
  pivot_wider(names_from = method, values_from = p_value) |>
  filter(!is.na(E0), !is.na(E1cr2), !is.na(E7null_cr2))
pair_one <- function(w0, w1, a, b, B = 2000) {
  stat <- function(i0, i1) {
    qa <- quantile(w0[[a]][i0], 0.05, names = FALSE)
    qb <- quantile(w0[[b]][i0], 0.05, names = FALSE)
    mean(w1[[a]][i1] <= qa) - mean(w1[[b]][i1] <= qb)
  }
  est <- stat(seq_len(nrow(w0)), seq_len(nrow(w1)))
  boot <- replicate(B, stat(sample(nrow(w0), replace = TRUE),
                            sample(nrow(w1), replace = TRUE)))
  c(diff = est, se = sd(boot))
}
pairs <- list()
for (dz in unique(wide$design)) for (n in unique(wide$N)) {
  w0 <- wide |> filter(design == dz, N == n, c_bm == 0)
  for (cb in c(0.30, 0.45)) {
    w1 <- wide |> filter(design == dz, N == n, abs(c_bm - cb) < 1e-9)
    for (b in c('E0', 'E1cr2')) {
      pairs[[length(pairs) + 1]] <- tibble(
        design = dz, N = n, c_bm = cb, vs = b,
        !!!as.list(pair_one(w0, w1, 'E7null_cr2', b)))
    }
  }
}
pairs <- bind_rows(pairs) |>
  mutate(vs = c(E0 = 'A0', E1cr2 = 'A2')[vs])
print(as.data.frame(pairs), digits = 3)
write.csv(pairs, file.path(out_dir, 'size-adjusted-power-f1-paired.csv'),
          row.names = FALSE)
cat('wrote', file.path(out_dir, 'size-adjusted-power-f1-paired.csv'), '\n')

## Pooled difference per design: the mean of the four (N, c_bm) cell
## differences (referee report 2026-11-15, M2). Within a job, the null
## and alternative cells at one design and N share replicate r's
## random-number stream, so the bootstrap resamples replicate indices
## jointly across the three cells at each N; the two N are resampled
## independently.
pool_one <- function(wd, a, b, B = 2000) {
  ns <- sort(unique(wd$N))
  by_n <- lapply(ns, function(n) {
    w <- wd[wd$N == n, ]
    list(reps = sort(unique(w$rep)),
         cells = split(w, sprintf('%.2f', w$c_bm)))
  })
  stat <- function(draw) {
    mean(unlist(lapply(seq_along(by_n), function(k) {
      cl <- by_n[[k]]$cells
      pick <- function(x) x[match(draw[[k]], x$rep, nomatch = 0), ]
      w0 <- pick(cl[['0.00']])
      vapply(c('0.30', '0.45'), function(cb) {
        w1 <- pick(cl[[cb]])
        mean(w1[[a]] <= quantile(w0[[a]], 0.05, names = FALSE)) -
          mean(w1[[b]] <= quantile(w0[[b]], 0.05, names = FALSE))
      }, numeric(1))
    })))
  }
  est <- stat(lapply(by_n, `[[`, 'reps'))
  boot <- replicate(B, stat(lapply(by_n, function(x)
    sample(x$reps, replace = TRUE))))
  c(diff = est, se = sd(boot))
}
pooled <- list()
for (dz in unique(wide$design)) for (b in c('E0', 'E1cr2')) {
  pooled[[length(pooled) + 1]] <- tibble(
    design = dz, vs = b,
    !!!as.list(pool_one(wide |> filter(design == dz), 'E7null_cr2', b)))
}
pooled <- bind_rows(pooled) |>
  mutate(vs = c(E0 = 'A0', E1cr2 = 'A2')[vs])
print(as.data.frame(pooled), digits = 3)
write.csv(pooled, file.path(out_dir, 'size-adjusted-power-f1-pooled.csv'),
          row.names = FALSE)
cat('wrote', file.path(out_dir, 'size-adjusted-power-f1-pooled.csv'), '\n')
