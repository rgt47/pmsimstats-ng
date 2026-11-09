## analysis/scripts/carryover-sensitivity/43-g10-coverage.R
##
## Coverage of the nominal 95 percent CR2 interval for G10, and for G9
## as a check, over the predictive-only alternative cells of the
## supplementary study (referee report 2026-11-08, M5). Reads the
## per-replicate G9/G10 refits of 41-run-g10-backfill.R; nothing is
## re-simulated. The method is that of 34-referee-reanalysis.R: the
## Satterthwaite df is recovered from the stored (t, p) pair where
## |t| >= 0.5, and the remaining replicates are classified as
## covering, not covering or ambiguous, which gives lower and upper
## bounds. The G9 rows must reproduce m3-coverage.csv.
##
## Run from the repository root:
##   Rscript analysis/scripts/carryover-sensitivity/43-g10-coverage.R

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

root <- here::here()
cs <- file.path(root, 'analysis/data/carryover-sensitivity')
out_dir <- file.path(cs, 'referee-reanalysis')

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
         df_median = median(df, na.rm = TRUE),
         n_ambiguous = sum(is.na(cov)),
         coverage_lower = mean(cov %in% TRUE),
         coverage_upper = mean(cov %in% TRUE | is.na(cov)),
         coverage_normal_approx = mean(z < 1.96))
}

## sigma_BR = 8 in the simcore2026 parameter set, as in
## 34-referee-reanalysis.R (grid$sigma_BR).
sigma_br <- 8
bf <- readRDS(file.path(cs, 'g10-backfill/all.rds')) |>
  filter(source == 'gf', c_bm > 0, converged)
cov_tab <- bf |>
  group_by(cell_id, sub, N, c_bm, method) |>
  group_modify(~ coverage_block(.x, theta = -unique(.y$c_bm) * sigma_br)) |>
  ungroup() |>
  mutate(method = c(E7cr2 = 'G9', E7null_cr2 = 'G10')[method])

## Check: the G9 bounds reproduce the stored reanalysis.
m3 <- read.csv(file.path(out_dir, 'm3-coverage.csv')) |>
  filter(level == 'cell', method == 'G9') |>
  select(cell_id, lo_ref = coverage_lower, hi_ref = coverage_upper)
chk <- cov_tab |> filter(method == 'G9') |> inner_join(m3, by = 'cell_id')
stopifnot(nrow(chk) == sum(cov_tab$method == 'G9'),
          max(abs(chk$coverage_lower - chk$lo_ref)) < 1e-12,
          max(abs(chk$coverage_upper - chk$hi_ref)) < 1e-12)

cov_sum <- cov_tab |>
  group_by(method, N) |>
  summarise(cells = n(), mean_cov_lower = mean(coverage_lower),
            mean_cov_upper = mean(coverage_upper),
            min_cov_lower = min(coverage_lower), .groups = 'drop')
print(as.data.frame(cov_sum))
write.csv(bind_rows(cov_tab |> mutate(level = 'cell'),
                    cov_sum |> mutate(level = 'summary')) |>
            select(level, everything()),
          file.path(out_dir, 'm5-coverage-g10.csv'), row.names = FALSE)
cat('wrote', file.path(out_dir, 'm5-coverage-g10.csv'), '\n')
