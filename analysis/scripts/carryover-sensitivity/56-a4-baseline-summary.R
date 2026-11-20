## analysis/scripts/carryover-sensitivity/56-a4-baseline-summary.R
##
## Tables for docs/77 (A4 with and without the baseline row), from the
## output of 55-run-a4-baseline.R and, for the published analysis with
## corCAR1 (E0pnc) on the same trials, 54-run-a0-published.R. Prints
## Markdown tables and the numbers quoted in the text; writes nothing.
##
## Usage (from the repository root):
##   Rscript analysis/scripts/carryover-sensitivity/56-a4-baseline-summary.R
suppressMessages({library(dplyr); library(tidyr); library(knitr)})
D <- 'analysis/data/carryover-sensitivity/published'
res   <- readRDS(file.path(D, 'a4-baseline/all.rds'))
guard <- readRDS(file.path(D, 'a4-baseline/guard.rds'))
meta  <- readRDS(file.path(D, 'a4-baseline/meta.rds'))
stopifnot(meta$param_set == 'published')
a0    <- readRDS(file.path(D, 'a0-published/all.rds'))
a0g   <- readRDS(file.path(D, 'a0-published/guard.rds'))
SIGMA_BR <- 8

cat(sprintf('Run %s, %d replicates per cell, guard tol %g\n',
            format(meta$date), meta$reps, meta$guard_tol))
cat(sprintf('Guard: %d of %d (cell, rep) pass; A4 half-life mismatches %d\n',
            sum(guard$match), nrow(guard), sum(!guard$h_match)))

ok <- guard |> filter(match) |> select(cell_id, rep)
## E0pnc on the same trials, only where script 54's own guard held.
a0_ok <- a0g |> filter(match) |> select(cell_id, rep)
e0 <- a0 |> filter(method == 'E0pnc') |>
  semi_join(a0_ok, by = c('cell_id', 'rep')) |>
  select(cell_id, rep, method, estimate, std_error, p_value, converged,
         panel, design, N, c_bm, t1half)
d <- res |> semi_join(ok, by = c('cell_id', 'rep')) |>
  select(-best_t_half) |>
  bind_rows(e0 |> semi_join(ok, by = c('cell_id', 'rep'))) |>
  mutate(rej = as.numeric(p_value < 0.05))
## Paired comparisons keep only replicates in which both arms converged.
both_conv <- function(dd, a, b) dd |> filter(method %in% c(a, b)) |>
  group_by(cell_id, rep) |> filter(n() == 2, all(converged)) |> ungroup()
nonconv <- d |> group_by(method) |>
  summarise(nonconv = sum(!converged), n = n(), .groups = 'drop')
cat('\nNon-converged fits by method:\n'); print(as.data.frame(nonconv))

rate_tab <- function(dd, methods) dd |> filter(method %in% methods,
                                              converged) |>
  group_by(design, N, method) |>
  summarise(r = mean(rej), n = n(), .groups = 'drop') |>
  mutate(cell = sprintf('%.3f', r)) |>
  select(design, N, method, cell) |>
  pivot_wider(names_from = method, values_from = cell)

## Pooled paired difference b - a with its SE, pooling by replicate
## (cells share random-number streams through the seed).
paired <- function(dd, a, b) {
  x <- both_conv(dd, a, b) |>
    select(cell_id, rep, method, rej) |>
    pivot_wider(names_from = method, values_from = rej) |>
    mutate(dif = .data[[b]] - .data[[a]])
  ## Replicate r of every cell shares seed r, so the per-replicate
  ## means are independent and their SE carries the cross-cell
  ## correlation.
  per_rep <- x |> group_by(rep) |> summarise(m = mean(dif),
                                             .groups = 'drop')
  c(diff = mean(per_rep$m), se = sd(per_rep$m) / sqrt(nrow(per_rep)),
    n = nrow(x))
}
pooled <- function(dd, m) {
  x <- dd |> filter(method == m, converged)
  r <- mean(x$rej)
  c(rate = r, n = nrow(x))
}

M_ALL <- c('A2', 'A2bl', 'A3', 'A3bl', 'A4', 'A4bl', 'E0pnc')
null <- d |> filter(c_bm == 0)
cat('\n## Size: null cells (t1/2 = 1 week)\n\n')
print(kable(rate_tab(null, M_ALL), format = 'pipe', align = 'llrrrrrrr'))
cat('\nPooled null rejection:\n')
for (m in M_ALL) {
  p <- pooled(null, m)
  cat(sprintf('  %-6s %.4f (n %d, MCSE %.4f)\n', m, p['rate'], p['n'],
              sqrt(p['rate'] * (1 - p['rate']) / p['n'])))
}
for (pr in list(c('A4', 'A4bl'), c('A2', 'A2bl'), c('A3', 'A3bl'),
                c('A4bl', 'E0pnc'))) {
  p <- paired(null, pr[1], pr[2])
  cat(sprintf('  size %s - %s: %+.4f (paired SE %.4f, n %d)\n',
              pr[2], pr[1], p['diff'], p['se'], p['n']))
}

alt <- d |> filter(c_bm > 0)
ref <- alt |> filter(c_bm == 0.45, t1half == 1)
cat('\n## Power at the six reference alternatives (c_bm 0.45, t1/2 1)\n\n')
print(kable(rate_tab(ref, M_ALL), format = 'pipe', align = 'llrrrrrrr'))
cat('\nPooled power, reference alternatives:\n')
for (m in M_ALL) {
  p <- pooled(ref, m)
  cat(sprintf('  %-6s %.4f (n %d)\n', m, p['rate'], p['n']))
}
for (pr in list(c('A4', 'A4bl'), c('A2', 'A2bl'), c('A3', 'A3bl'),
                c('A2bl', 'A4bl'), c('A4bl', 'E0pnc'))) {
  p <- paired(ref, pr[1], pr[2])
  cat(sprintf('  power %s - %s: %+.4f (paired SE %.4f, n %d)\n',
              pr[2], pr[1], p['diff'], p['se'], p['n']))
}

cat('\n## Power, all 24 alternative cells: A4, A4bl and the gain\n\n')
pw <- alt |> filter(method %in% c('A4', 'A4bl'), converged) |>
  group_by(design, N, c_bm, t1half, method) |>
  summarise(r = mean(rej), .groups = 'drop') |>
  pivot_wider(names_from = method, values_from = r) |>
  mutate(gain = A4bl - A4) |> arrange(design, c_bm, t1half, N)
print(kable(pw |> mutate(across(c(A4, A4bl, gain), ~ sprintf('%.3f', .x))),
            format = 'pipe'))
cat(sprintf('\nGain range %.3f to %.3f; mean %.3f; negative in %d of %d cells\n',
            min(pw$gain), max(pw$gain), mean(pw$gain), sum(pw$gain < 0),
            nrow(pw)))
gd <- pw |> group_by(design) |> summarise(mean_gain = mean(gain),
                                          .groups = 'drop')
print(as.data.frame(gd))
gt <- pw |> group_by(t1half) |> summarise(mean_gain = mean(gain),
                                          .groups = 'drop')
print(as.data.frame(gt))
cat('\nPaired pooled power difference A4bl - A4, all alternatives:\n')
p <- paired(alt, 'A4', 'A4bl')
cat(sprintf('  %+.4f (paired SE %.4f, n %d)\n', p['diff'], p['se'], p['n']))

cat('\n## Estimation at c_bm = 0.45 (theta = -c_bm * 8 = -3.6)\n\n')
est <- alt |> filter(c_bm == 0.45, method %in% c('A3', 'A3bl', 'A4', 'A4bl'),
                     converged) |>
  group_by(method, t1half) |>
  summarise(mean_est = mean(estimate), emp_sd = sd(estimate),
            mean_se = mean(std_error), .groups = 'drop') |>
  mutate(bias = mean_est + 0.45 * SIGMA_BR, kappa = mean_se / emp_sd)
print(kable(est, format = 'pipe', digits = 3))
cat('\nBy design at the reference alternatives (t1/2 = 1):\n')
est_d <- ref |> filter(method %in% c('A3', 'A3bl', 'A4', 'A4bl'),
                       converged) |>
  group_by(design, method) |>
  summarise(mean_est = mean(estimate), emp_sd = sd(estimate),
            mean_se = mean(std_error), .groups = 'drop') |>
  mutate(kappa = mean_se / emp_sd)
print(kable(est_d, format = 'pipe', digits = 3))
cat('\nkappa under the null, by method:\n')
kn <- null |> filter(method %in% c('A2', 'A2bl', 'A3', 'A3bl', 'A4', 'A4bl'),
                     converged) |>
  group_by(method) |>
  summarise(emp_sd = sd(estimate), mean_se = mean(std_error),
            .groups = 'drop') |> mutate(kappa = mean_se / emp_sd)
print(kable(kn, format = 'pipe', digits = 3))

cat('\n## Selected half-life (A4 against A4bl), by true half-life\n\n')
hs <- res |> semi_join(ok, by = c('cell_id', 'rep')) |>
  filter(method %in% c('A4', 'A4bl'), converged) |>
  count(t1half, method, best_t_half) |>
  group_by(t1half, method) |> mutate(p = n / sum(n)) |> ungroup() |>
  select(-n) |>
  pivot_wider(names_from = best_t_half, values_from = p, values_fill = 0)
print(kable(hs, format = 'pipe', digits = 2))
same_h <- res |> semi_join(ok, by = c('cell_id', 'rep')) |>
  filter(method %in% c('A4', 'A4bl')) |>
  select(cell_id, rep, method, best_t_half) |>
  pivot_wider(names_from = method, values_from = best_t_half)
cat(sprintf('\nSame half-life selected in %.1f%% of replicates\n',
            100 * mean(same_h$A4 == same_h$A4bl, na.rm = TRUE)))
