## analysis/scripts/carryover-sensitivity/49-coverage-family1-s.R
##
## Coverage of the nominal 95 percent CR2 interval for A3 (stored E2cr2)
## and A4 (E7null_cr2) in the two families that 43-g10-coverage.R does
## not cover (referee report 2026-11-14, M2):
##   family 1, the alternative cells of the 30-cell grid (job 27; A3 from
##     tier1-g-vs-f-e0/all.rds, A4 from the backfill refits, source t1);
##   family S, the alternative cells of the independent-seed Hybrid run
##     (referee-variants/block-A/all.rds).
## Nothing is re-simulated. The method is that of 34-referee-reanalysis.R
## and 43-g10-coverage.R: the Satterthwaite df is recovered from the
## stored (t, p) pair where |t| >= 0.5, and the remaining replicates are
## classified as covering, not covering or ambiguous, which gives lower
## and upper bounds. The target is theta = -c_bm * sigma_BR, sigma_BR = 8.
##
## Run from the repository root:
##   Rscript analysis/scripts/carryover-sensitivity/49-coverage-family1-s.R

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

root <- here::here()
cs <- file.path(root, 'analysis/data/carryover-sensitivity')
## Parameter set from PAPER02_PARAMS (p02-sources.R). Under the published
## set A4 was fitted in job 27's counterpart, so it is read from there.
source(file.path(root, 'analysis/scripts/carryover-sensitivity',
                 'p02-sources.R'))
out_dir <- p02_derived('referee-reanalysis')

solve_df <- function(t, p) {
  if (is.na(t) || is.na(p) || abs(t) < 0.5 || p <= 0 || p >= 1) {
    return(NA_real_)
  }
  f <- function(ldf) 2 * pt(-abs(t), exp(ldf)) - p
  lo <- f(log(1)); hi <- f(log(1e6))
  if (sign(lo) == sign(hi)) return(NA_real_)
  exp(uniroot(f, c(log(1), log(1e6)), tol = 1e-8)$root)
}

coverage_block <- function(d, theta) {
  d <- d |> filter(!is.na(estimate), !is.na(std_error))
  z <- abs(d$estimate - theta) / d$std_error
  df <- map2_dbl(d$t, d$p_value, solve_df)
  crit <- ifelse(is.na(df), NA_real_, qt(0.975, df))
  cov <- ifelse(!is.na(crit), z < crit,
                ifelse(z < 1.96, TRUE,
                       ifelse(z > qt(0.975, 5), FALSE, NA)))
  tibble(n = nrow(d), n_df_solved = sum(!is.na(df)),
         n_ambiguous = sum(is.na(cov)),
         coverage_lower = mean(cov %in% TRUE),
         coverage_upper = mean(cov %in% TRUE | is.na(cov)),
         emp_sd = sd(d$estimate), mean_se = mean(d$std_error))
}

sigma_br <- 8

## Family 1: job 27 alternatives. A3 from the stored per-replicate file,
## A4 from the backfill refits of the same trials.
t1_all <- readRDS(file.path(p02_run('tier1-g-vs-f-e0'), 'all.rds'))
a27 <- t1_all |>
  filter(method == 'E2cr2', c_bm > 0, converged) |>
  select(cell_id, rep, method, estimate, std_error, p_value, t,
         design, N, c_bm, t1half)
key27 <- a27 |> distinct(cell_id, design, N, c_bm, t1half)
g27_src <- if (p02_backfill()) {
  readRDS(file.path(cs, 'g10-backfill/all.rds')) |> filter(source == 't1')
} else t1_all
g27 <- g27_src |>
  filter(method == 'E7null_cr2', c_bm > 0, converged) |>
  select(cell_id, rep, method, estimate, std_error, p_value, t) |>
  inner_join(key27, by = 'cell_id')
stopifnot(nrow(g27) > 0, setequal(unique(g27$cell_id), unique(a27$cell_id)))
f1 <- bind_rows(a27, g27) |> mutate(family = '1')

## Family S: independent-seed alternatives (block A).
bA <- readRDS(file.path(p02_run('referee-variants/block-A'), 'all.rds')) |>
  filter(method %in% c('E2cr2', 'E7null_cr2'), converged) |>
  mutate(design = 'Hybrid',
         N = as.integer(sub('.*-N(\\d+)-.*', '\\1', cell_id)),
         c_bm = as.numeric(sub('.*-cbm', '', cell_id)),
         t1half = as.numeric(sub('^alt-t([0-9.]+)-.*', '\\1', cell_id)),
         family = 'S') |>
  select(cell_id, rep, method, estimate, std_error, p_value, t,
         design, N, c_bm, t1half, family)
stopifnot(!anyNA(bA$N), !anyNA(bA$c_bm), !anyNA(bA$t1half))

cov_tab <- bind_rows(f1, bA) |>
  group_by(family, cell_id, design, N, c_bm, t1half, method) |>
  group_modify(~ coverage_block(.x, theta = -unique(.y$c_bm) * sigma_br)) |>
  ungroup() |>
  mutate(kappa = mean_se / emp_sd,
         method = c(E2cr2 = 'G8', E7null_cr2 = 'G10')[method])

cov_sum <- cov_tab |>
  group_by(family, method) |>
  summarise(cells = n(), mean_cov_lower = mean(coverage_lower),
            mean_cov_upper = mean(coverage_upper),
            min_cov_lower = min(coverage_lower),
            min_kappa = min(kappa), .groups = 'drop')
print(as.data.frame(cov_sum))
print(as.data.frame(cov_tab |>
  filter(family == '1', design == 'Hybrid', N == 70, c_bm == 0.45,
         t1half == 1) |>
  select(method, coverage_lower, coverage_upper, kappa)))
write.csv(bind_rows(cov_tab |> mutate(level = 'cell'),
                    cov_sum |> mutate(level = 'summary')) |>
            select(level, everything()),
          file.path(out_dir, 'm6-coverage-family1-s.csv'), row.names = FALSE)
cat('wrote', file.path(out_dir, 'm6-coverage-family1-s.csv'), '\n')
