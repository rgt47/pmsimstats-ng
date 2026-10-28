## analysis/scripts/carryover-sensitivity/34-referee-reanalysis.R
##
## Reanalysis of the stored per-replicate and summary data of paper 02
## (analysis/report/02-carryover-sensitivity/report.Rmd) to answer the
## 2026-10-27 referee comments. Every number is computed from stored
## data; nothing is re-simulated. Outputs are written as CSV files to
## analysis/data/carryover-sensitivity/referee-reanalysis/ and printed
## to stdout.
##
## Stored method codes: E1 = G1, E3 = G2, E2 = G3, E7 = G4, E9 = G5,
## E1cr2 = G6, E3cr2 = G7, E2cr2 = G8, E7cr2 = G9; E0 is the published
## analysis; F1, F1A, F1A0, F1A0B the expectancy-term family (F1A0B is
## summary-only: it is calibrated from F1A0's bootstrap draws and has
## no per-replicate row).
##
## Conventions shared with the drivers (32-run-g9-gap-fill.R):
## reject = p_value < 0.05 among replicates with a non-missing
## p-value; kappa = mean(std_error) / sd(estimate). The robust kappa
## here is median(std_error) / (1.4826 * MAD(estimate)); note that
## stats::mad() already includes the 1.4826 constant.
##
## Run from the repository root:
##   Rscript analysis/scripts/carryover-sensitivity/34-referee-reanalysis.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})

root <- here::here()
cs <- file.path(root, 'analysis/data/carryover-sensitivity')
out_dir <- file.path(cs, 'referee-reanalysis')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(root,
  'analysis/scripts/carryover-sensitivity/spec-labels.R'))
m_lab <- c(E0 = 'E0', E1 = 'G1', E3 = 'G2', E2 = 'G3', E7 = 'G4',
           E9 = 'G5', E1cr2 = 'G6', E3cr2 = 'G7', E2cr2 = 'G8',
           E7cr2 = 'G9', F1 = 'F1', F1A = 'F1A', F1A0 = 'F1A0',
           F1A0B = 'F1A0B')
g <- function(code) {
  ifelse(code %in% names(m_lab), unname(m_lab[code]), code)
}

emit <- function(df, name) {
  df <- as.data.frame(df)
  path <- file.path(out_dir, name)
  write.csv(df, path, row.names = FALSE)
  cat('\n==========', name, '(', nrow(df), 'rows ) ==========\n')
  print(df, row.names = FALSE, digits = 4)
  invisible(df)
}
options(width = 220)

## ---------------------------------------------------------------------
## Load stored data
## ---------------------------------------------------------------------

sgf <- readRDS(file.path(cs, 'g9-gap-fill/summary.rds'))
gf_jobs <- c(31, 32, 33, 40)
gf_all <- map_dfr(gf_jobs, function(j) {
  readRDS(file.path(cs, sprintf('g9-gap-fill/job-%d/all.rds', j))) |>
    mutate(job = j)
})
## Cell-level covariates for the per-replicate gap-fill rows come
## from the summary; cell_id is unique within the gap-fill study.
gf_cells <- sgf |>
  distinct(cell_id, job, block, sub, N, c_bm, t1half, rho, dgp_arch,
           carryover_form, weibull_shape, dropout_rate, dropout_mech,
           dropout_per, c_pb, c_tv, c_bl)
stopifnot(!anyDuplicated(gf_cells$cell_id))
gf_all <- gf_all |>
  select(-job, -block, -sub) |>
  left_join(gf_cells, by = 'cell_id')
stopifnot(!any(is.na(gf_all$N)))

s21 <- readRDS(file.path(cs, 'g-vs-f-comparison/summary.rds'))
a21 <- readRDS(file.path(cs, 'g-vs-f-comparison/all.rds'))
s27 <- readRDS(file.path(cs, 'tier1-g-vs-f-e0/summary.rds'))
a27 <- readRDS(file.path(cs, 'tier1-g-vs-f-e0/all.rds'))
s28 <- readRDS(file.path(cs, 'g9-prognostic-stress/summary.rds'))
a28 <- readRDS(file.path(cs, 'g9-prognostic-stress/all.rds'))
s30 <- readRDS(file.path(cs, 'g9-random-effects/summary.rds'))
a30 <- readRDS(file.path(cs, 'g9-random-effects/all.rds'))
grid <- readRDS(file.path(root, 'analysis/data/02-grid-summary.rds'))$summary
raw1 <- readRDS(file.path(root,
  'analysis/scripts/carryover-sensitivity/output/01-factorial.rds'))$results

## Stage 4 cell covariates, transcribed from 06-g-vs-f-comparison.R
## (base cell: Hybrid, N = 70, t1half = 0.5; t1hifhi cells t1half = 2;
## olbdc cells the OL+BDC design).
s4_cells <- tibble(cell = unique(s21$cell)) |>
  mutate(design = ifelse(grepl('^olbdc', cell), 'OLBDC', 'Hybrid'),
         N = 70,
         t1half = ifelse(grepl('^t1hifhi', cell), 2, 0.5),
         null = grepl('_null$', cell))
s4_null8 <- c('hybrid_null', 'olbdc_null', 'progbl_null',
              'progtvconst_null', 'progtvgrow_null', 'progpb_null',
              't1hifhi_null', 'shapemis_null')

## ---------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------

rej <- function(p) as.integer(p < 0.05)
mcse <- function(p, n) sqrt(p * (1 - p) / n)

## Paired comparison of rejection indicators for two methods fitted to
## the same replicate; rows with a missing p-value on either side are
## dropped from the pair set.
paired_rej <- function(d, m_a, m_b, keys = 'rep') {
  a <- d |> filter(method == m_a, !is.na(p_value)) |>
    transmute(across(all_of(keys)), ra = rej(p_value))
  b <- d |> filter(method == m_b, !is.na(p_value)) |>
    transmute(across(all_of(keys)), rb = rej(p_value))
  j <- inner_join(a, b, by = keys)
  dd <- j$ra - j$rb
  tibble(method_a = g(m_a), method_b = g(m_b), n_pairs = nrow(j),
         rate_a = mean(j$ra), rate_b = mean(j$rb),
         a_only = sum(j$ra == 1 & j$rb == 0),
         b_only = sum(j$ra == 0 & j$rb == 1),
         diff = mean(dd), paired_se = sd(dd) / sqrt(nrow(j)),
         naive_se = sqrt(mean(j$ra) * (1 - mean(j$ra)) / nrow(j) +
                         mean(j$rb) * (1 - mean(j$rb)) / nrow(j)))
}

robust_kappa <- function(se, est) {
  ok <- !is.na(se) & !is.na(est)
  median(se[ok]) / mad(est[ok], constant = 1.4826)
}

## Satterthwaite df is not stored; it is recovered from the stored
## (t, p) pair by solving p = 2 * pt(-|t|, df). The solution is
## ill-determined when |t| is small, so it is attempted only for
## |t| >= 0.5; for the rest the CI is classified as covering (|z| <
## 1.96, covers for any df), not covering (|z| > qt(0.975, 5)), or
## ambiguous.
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

## ---------------------------------------------------------------------
## Item 1: pooled size at the four 2,000-replicate Hybrid null cells
## ---------------------------------------------------------------------

a2 <- gf_all |> filter(sub == 'A2-size')
stopifnot(length(unique(a2$cell_id)) == 4)
pooled <- a2 |>
  filter(method %in% c('E7cr2', 'E2cr2', 'E1cr2'), !is.na(p_value)) |>
  group_by(method) |>
  summarise(n = n(), reject = mean(rej(p_value)), .groups = 'drop') |>
  mutate(mcse = mcse(reject, n), item = 'pooled', detail = 'all four cells',
         method = g(method)) |>
  select(item, detail, method, n, reject, mcse)
pr <- paired_rej(a2, 'E7cr2', 'E2cr2', keys = c('cell_id', 'rep'))
paired_row <- tibble(item = 'paired G9 vs G8', detail = sprintf(
  'n_pairs=%d; G9 only=%d; G8 only=%d; naive_se=%.4f',
  pr$n_pairs, pr$a_only, pr$b_only, pr$naive_se),
  method = 'G9-G8', n = pr$n_pairs, reject = pr$diff, mcse = pr$paired_se)
cond <- a2 |> filter(method == 'E7cr2', !is.na(p_value)) |>
  group_by(h) |>
  summarise(n = n(), reject = mean(rej(p_value)), .groups = 'drop') |>
  mutate(mcse = mcse(reject, n), item = 'G9 conditional on selected h',
         detail = sprintf('h=%g', h), method = 'G9') |>
  select(item, detail, method, n, reject, mcse)
per_cell <- a2 |>
  filter(method %in% c('E7cr2', 'E2cr2'), !is.na(p_value)) |>
  group_by(cell_id, method) |>
  summarise(n = n(), reject = mean(rej(p_value)), .groups = 'drop') |>
  mutate(mcse = mcse(reject, n), item = 'per cell', detail = cell_id,
         method = g(method)) |>
  select(item, detail, method, n, reject, mcse)
nl <- sgf |>
  filter(c_bm == 0, sub != 'C1-ladder', method %in% c('E7cr2', 'E2cr2')) |>
  select(cell_id, sub, method, reject) |>
  pivot_wider(names_from = method, values_from = reject) |>
  mutate(d = E7cr2 - E2cr2)
nl_rows <- nl |>
  transmute(item = 'non-ladder null cells: G9 - G8', detail = cell_id,
            method = sub, n = NA_integer_, reject = d, mcse = NA_real_)
nl_sum <- tibble(
  item = 'non-ladder null cells: summary',
  detail = c('cells', 'G9 > G8', 'G9 < G8', 'tie', 'mean diff',
             'cells excl. C2-largeN', 'G9 > G8 excl. C2-largeN'),
  method = 'G9-G8', n = NA_integer_,
  reject = c(nrow(nl), sum(nl$d > 0), sum(nl$d < 0), sum(nl$d == 0),
             mean(nl$d), sum(nl$sub != 'C2-largeN'),
             sum(nl$d > 0 & nl$sub != 'C2-largeN')),
  mcse = NA_real_)
emit(bind_rows(pooled, paired_row, cond, per_cell, nl_sum, nl_rows),
     'm2-pooled-size.csv')

## ---------------------------------------------------------------------
## Item 2: AIC-selected half-life frequencies
## ---------------------------------------------------------------------

sel_src <- bind_rows(
  gf_all |> filter(method %in% c('E7', 'E7cr2')) |>
    transmute(source = sprintf('gap-fill job %d', job), scenario = sub,
              cell = cell_id, true_t1half = t1half, design = 'Hybrid',
              N, method, h),
  a21 |> filter(method %in% c('E7', 'E7cr2')) |>
    left_join(s4_cells, by = 'cell') |>
    transmute(source = 'g-vs-f (Stage 4)', scenario = cell, cell,
              true_t1half = t1half, design, N, method, h),
  a27 |> filter(method %in% c('E7', 'E7cr2')) |>
    transmute(source = 'tier1 (Stage 3)', scenario = sprintf(
      'cbm%.2f', c_bm), cell = cell_id, true_t1half = t1half,
      design, N, method, h),
  a28 |> filter(method %in% c('E7', 'E7cr2')) |>
    transmute(source = 'prognostic-stress (Stage 5)',
              scenario = sprintf('%s s=%.1f cbm%.2f', pathway, strength, c_bm),
              cell = cell_id, true_t1half = 0.5, design = 'Hybrid', N,
              method, h),
  a30 |> filter(method == 'G9') |>
    transmute(source = 'random-effects (Stage 5)',
              scenario = sprintf('%s s=%.1f cbm%.2f', pathway, strength, c_bm),
              cell = cell_id, true_t1half = 0.5, design = 'Hybrid', N,
              method = 'E7cr2', h)) |>
  filter(!is.na(h)) |>
  mutate(method = g(method))
sel_by <- function(d, ...) {
  d |> group_by(...) |>
    summarise(n = n(), p_0.25 = mean(h == 0.25), p_0.5 = mean(h == 0.5),
              p_1 = mean(h == 1), p_2 = mean(h == 2), .groups = 'drop')
}
sel_tab <- bind_rows(
  sel_by(sel_src, level = 'by true half-life x design x N (all sources)',
         true_t1half, design, N, method),
  sel_by(sel_src, level = 'by source x scenario', source, scenario,
         true_t1half, design, N, method)) |>
  arrange(level, true_t1half, design, N, method)
emit(sel_tab, 'm9-selection-frequencies.csv')

## ---------------------------------------------------------------------
## Item 3: kappa at alternative cells; coverage of CR2 intervals
## ---------------------------------------------------------------------

gf_alt_sets <- list(
  'non-dropout' = sgf |> filter(c_bm > 0, sub != 'A1-dropout'),
  'non-dropout excl. meanmod' =
    sgf |> filter(c_bm > 0, sub != 'A1-dropout', sub != 'B3-meanmod'),
  'non-dropout excl. meanmod + patient-level dropout' =
    sgf |> filter(c_bm > 0, sub != 'B3-meanmod',
                  !(sub == 'A1-dropout' & dropout_per == 'occasion')),
  'all alternatives incl. both dropout conventions' =
    sgf |> filter(c_bm > 0))
k_sum <- map_dfr(names(gf_alt_sets), function(nm) {
  gf_alt_sets[[nm]] |>
    group_by(method, N) |>
    summarise(cells = n(), mean_kappa = mean(kappa), min_kappa = min(kappa),
              max_kappa = max(kappa), .groups = 'drop') |>
    mutate(cell_set = nm, method = g(method))
})
rk_cell <- gf_all |>
  filter(c_bm > 0) |>
  group_by(cell_id, sub, N, dropout_per, method) |>
  summarise(n = sum(!is.na(estimate)),
            kappa_recomputed = mean(std_error, na.rm = TRUE) /
              sd(estimate, na.rm = TRUE),
            robust_kappa = robust_kappa(std_error, estimate),
            .groups = 'drop') |>
  mutate(method = g(method))
rk_sum <- map_dfr(names(gf_alt_sets), function(nm) {
  ids <- unique(gf_alt_sets[[nm]]$cell_id)
  rk_cell |> filter(cell_id %in% ids) |>
    group_by(method, N) |>
    summarise(mean_robust_kappa = mean(robust_kappa),
              min_robust_kappa = min(robust_kappa), .groups = 'drop') |>
    mutate(cell_set = nm)
})
k_tab <- k_sum |> left_join(rk_sum, by = c('cell_set', 'method', 'N')) |>
  select(cell_set, method, N, cells, mean_kappa, min_kappa, max_kappa,
         mean_robust_kappa, min_robust_kappa) |>
  arrange(cell_set, method, N)
emit(k_tab, 'm3-kappa-alternatives.csv')
emit(rk_cell |> arrange(sub, cell_id, method), 'm3-kappa-per-cell.csv')
g9_min <- sgf |> filter(c_bm > 0, method == 'E7cr2') |>
  arrange(kappa) |> select(cell_id, sub, N, kappa) |> head(5)
cat('\nFive smallest G9 kappa among all gap-fill alternative cells:\n')
print(as.data.frame(g9_min))

## Coverage: theta_true = -c_bm * sigma_BR with sigma_BR = 8, the
## calibration used by 02-summarise-grid.R (grid$true_value); the
## estimand of G8/G9 is bm:Dbc, whose calibrated target is this value
## under the exponential DGP at the true half-life. Mean-moderation
## (B3) cells are reported but flagged.
sigma_br <- unique(grid$sigma_BR)
stopifnot(length(sigma_br) == 1)
cov_tab <- gf_all |>
  filter(c_bm > 0, method %in% c('E2cr2', 'E7cr2')) |>
  group_by(cell_id, sub, N, c_bm, dgp_arch, dropout_per, method) |>
  group_modify(~ coverage_block(.x, theta = -unique(.y$c_bm) * sigma_br)) |>
  ungroup() |>
  mutate(method = g(method), theta_true = -c_bm * sigma_br,
         note = ifelse(dgp_arch == 'mean_moderation',
                       'mean-moderation DGP: estimand mapping assumed', ''))
cov_sum <- cov_tab |> filter(dgp_arch != 'mean_moderation') |>
  group_by(method, N) |>
  summarise(cells = n(), mean_cov_lower = mean(coverage_lower),
            mean_cov_upper = mean(coverage_upper),
            mean_cov_normal = mean(coverage_normal_approx),
            min_cov_lower = min(coverage_lower), .groups = 'drop')
cat('\nCoverage summary (mvn cells only):\n'); print(as.data.frame(cov_sum))
emit(bind_rows(cov_tab |> mutate(level = 'cell'),
               cov_sum |> mutate(level = 'summary (mvn cells)')) |>
       select(level, everything()), 'm3-coverage.csv')

## ---------------------------------------------------------------------
## Item 4: heavy missingness, per-occasion MCAR 0.3, N = 35
## ---------------------------------------------------------------------

hm <- gf_all |>
  filter(sub == 'A1-dropout', dropout_per == 'occasion',
         dropout_mech == 'MCAR', dropout_rate == 0.3, N == 35) |>
  group_by(cell_id, c_bm, method) |>
  summarise(n = sum(!is.na(estimate)), nonconv = mean(!converged),
            sd_est = sd(estimate, na.rm = TRUE),
            mad_sd_est = mad(estimate, constant = 1.4826, na.rm = TRUE),
            median_abs_est = median(abs(estimate), na.rm = TRUE),
            max_abs_est = max(abs(estimate), na.rm = TRUE),
            mean_se = mean(std_error, na.rm = TRUE),
            median_se = median(std_error, na.rm = TRUE),
            kappa = mean_se / sd_est,
            robust_kappa = robust_kappa(std_error, estimate),
            reject = mean(rej(p_value), na.rm = TRUE),
            .groups = 'drop') |>
  mutate(method = g(method), cell_type = ifelse(c_bm == 0, 'null', 'alt')) |>
  arrange(c_bm, match(method, m_lab))
emit(hm, 'm10-heavy-missingness.csv')

## ---------------------------------------------------------------------
## Item 5: counts and claims
## ---------------------------------------------------------------------

row5 <- function(claim, detail, quantity, value) {
  tibble(claim = claim, detail = detail, quantity = quantity, value = value)
}
m5 <- list()

## (a)
null_cells_s3 <- s27 |> filter(c_bm == 0) |> distinct(cell_id)
m5$a <- bind_rows(
  row5('a', 'Stage 3 30-cell grid', 'distinct cells',
       n_distinct(s27$cell_id)),
  row5('a', 'Stage 3 30-cell grid', 'null cells (c_bm = 0)',
       nrow(null_cells_s3)))

## (b)
e0_s3 <- s27 |> filter(c_bm == 0, method == 'E0') |>
  transmute(claim = 'b', detail = sprintf('Stage 3 %s', cell_id),
            quantity = 'E0 null rejection', value = reject)
e0_s4 <- s21 |> filter(method == 'E0', grepl('_null$', cell)) |>
  transmute(claim = 'b', detail = sprintf('Stage 4 %s%s', cell,
              ifelse(cell %in% s4_null8, '', ' (not tabulated in paper)')),
            quantity = 'E0 null rejection', value = reject)
e0_all <- bind_rows(e0_s3, e0_s4)
m5$b <- bind_rows(e0_all,
  row5('b', 'Stage 3 + Stage 4 (8 tabulated null cells)', 'E0 min',
       min(c(e0_s3$value, e0_s4$value[!grepl('not tabulated', e0_s4$detail)]))),
  row5('b', 'Stage 3 + Stage 4 (8 tabulated null cells)', 'E0 max',
       max(c(e0_s3$value, e0_s4$value[!grepl('not tabulated', e0_s4$detail)]))),
  row5('b', 'Stage 3 only', 'E0 min', min(e0_s3$value)),
  row5('b', 'Stage 3 only', 'E0 max', max(e0_s3$value)),
  row5('b', 'Stage 4 8 null cells only', 'E0 min',
       min(e0_s4$value[!grepl('not tabulated', e0_s4$detail)])),
  row5('b', 'Stage 4 8 null cells only', 'E0 max',
       max(e0_s4$value[!grepl('not tabulated', e0_s4$detail)])))

## (c)
alt_nd <- gf_alt_sets[['non-dropout']] |>
  filter(method %in% c('E1', 'E1cr2', 'E2cr2', 'E7cr2', 'F1')) |>
  select(cell_id, sub, N, method, reject) |>
  pivot_wider(names_from = method, values_from = reject) |>
  mutate(d_G1 = E7cr2 - E1, d_G6 = E7cr2 - E1cr2, d_G8 = E7cr2 - E2cr2,
         d_F1 = E7cr2 - F1)
c_rows <- map_dfr(c('d_G1', 'd_G6', 'd_G8', 'd_F1'), function(v) {
  comp <- sub('d_', '', v)
  below <- alt_nd |> filter(.data[[v]] < 0)
  ties <- alt_nd |> filter(.data[[v]] == 0)
  bind_rows(
    row5('c', sprintf('non-dropout alternatives (%d cells)', nrow(alt_nd)),
         sprintf('cells where G9 < %s', comp), nrow(below)),
    row5('c', sprintf('non-dropout alternatives (%d cells)', nrow(alt_nd)),
         sprintf('cells where G9 == %s', comp), nrow(ties)),
    if (nrow(below) > 0) below |>
      transmute(claim = 'c', detail = cell_id,
                quantity = sprintf('G9 - %s where G9 below', comp),
                value = .data[[v]]) else NULL,
    if (nrow(ties) > 0) ties |>
      transmute(claim = 'c', detail = cell_id,
                quantity = sprintf('G9 - %s tie', comp), value = 0) else NULL)
})
m5$c <- c_rows

## (d), (f)
ref30 <- s30 |> filter(pathway == 'none', c_bm == 0.45, method == 'G9') |>
  select(N, ref = reject)
pw28 <- s28 |> filter(c_bm == 0.45, method == 'E7cr2') |>
  left_join(ref30, by = 'N') |>
  mutate(loss = reject - ref)
m5$d <- pw28 |> filter(strength == 0.6, pathway %in% c('progpb', 'progtvconst')) |>
  transmute(claim = 'd', detail = sprintf('%s s=0.6 N=%d (ref %.3f)', pathway, N, ref),
            quantity = 'G9 power minus no-confound reference', value = loss)
m5$f <- pw28 |> filter(strength == 0.3) |>
  transmute(claim = 'f', detail = sprintf('%s s=0.3 N=%d (ref %.3f)', pathway, N, ref),
            quantity = 'G9 power minus no-confound reference', value = loss)

## (e)
f_diff <- function(s, label, keys) {
  s |> filter(method %in% c('F1', 'F1A0B')) |>
    select(all_of(keys), method, reject) |>
    pivot_wider(names_from = method, values_from = reject) |>
    filter(!is.na(F1), !is.na(F1A0B), !is.nan(F1)) |>
    mutate(d = F1A0B - F1, src = label)
}
e_gf <- f_diff(sgf |> filter(c_bm > 0), 'gap-fill', c('cell_id')) |>
  rename(cell = cell_id)
e_21 <- f_diff(s21 |> filter(grepl('_alt', cell)), 'Stage 4', c('cell'))
e_27 <- f_diff(s27 |> filter(c_bm > 0), 'Stage 3', c('cell_id')) |>
  rename(cell = cell_id)
e_28 <- f_diff(s28 |> filter(c_bm > 0), 'prognostic-stress', c('cell_id')) |>
  rename(cell = cell_id)
e_all <- bind_rows(e_gf, e_21, e_27, e_28)
m5$e <- bind_rows(
  e_all |> transmute(claim = 'e', detail = sprintf('%s %s', src, cell),
                     quantity = 'F1A0B - F1 power', value = d),
  row5('e', 'all alternative cells with both', 'cells', nrow(e_all)),
  row5('e', 'all alternative cells with both', 'mean F1A0B - F1', mean(e_all$d)),
  row5('e', 'all alternative cells with both', 'min', min(e_all$d)),
  row5('e', 'all alternative cells with both', 'max', max(e_all$d)),
  row5('e', 'all alternative cells with both', 'cells F1A0B > F1',
       sum(e_all$d > 0)))

## (g)
g1 <- grid |> filter(carryover_form == 'exponential', dgp_arch == 'mvn',
                     spec == 'E1', c_bm > 0) |>
  select(design, N, c_bm, t1half, power) |>
  pivot_wider(names_from = t1half, values_from = power, names_prefix = 'p') |>
  mutate(rel_loss_1 = (p0 - p1) / p0, rel_loss_0.5 = (p0 - p0.5) / p0)
m5$g <- bind_rows(
  g1 |> transmute(claim = 'g', detail = sprintf('%s N=%d cbm=%.2f', design, N, c_bm),
                  quantity = 'G1 relative power loss, t1/2 = 1 vs 0',
                  value = rel_loss_1),
  g1 |> transmute(claim = 'g', detail = sprintf('%s N=%d cbm=%.2f', design, N, c_bm),
                  quantity = 'G1 relative power loss, t1/2 = 0.5 vs 0',
                  value = rel_loss_0.5),
  g1 |> group_by(design) |>
    summarise(min1 = min(rel_loss_1), max1 = max(rel_loss_1),
              min05 = min(rel_loss_0.5), max05 = max(rel_loss_0.5),
              .groups = 'drop') |>
    pivot_longer(-design, names_to = 'q', values_to = 'value') |>
    transmute(claim = 'g', detail = design,
              quantity = c(min1 = 'min rel loss t1/2=1',
                           max1 = 'max rel loss t1/2=1',
                           min05 = 'min rel loss t1/2=0.5',
                           max05 = 'max rel loss t1/2=0.5')[q], value))

## (h)
cfg_keys <- c('carryover_form', 'weibull_shape', 'dgp_arch', 't1half',
              'design', 'N', 'c_bm')
grid_d <- grid |> filter(!(carryover_form == 'weibull' & weibull_shape == 1))
n_cfg <- grid_d |> distinct(across(all_of(cfg_keys))) |> nrow()
w31 <- grid_d |> filter(spec %in% c('E1', 'E2')) |>
  select(all_of(cfg_keys), spec, power) |>
  pivot_wider(names_from = spec, values_from = power) |>
  mutate(d = E2 - E1)
twin <- grid |>
  filter(carryover_form %in% c('exponential', 'weibull'),
         weibull_shape == 1) |>
  select(all_of(cfg_keys), spec, power) |>
  pivot_wider(names_from = carryover_form, values_from = power) |>
  mutate(ad = abs(weibull - exponential))
m5$h <- bind_rows(
  row5('h', 'Stage 1 grid', 'cells before removing Weibull k=1 duplicates',
       grid |> distinct(across(all_of(cfg_keys))) |> nrow()),
  row5('h', 'Stage 1 grid', 'distinct configurations after removal', n_cfg),
  row5('h', 'Stage 1 grid (distinct configs)', 'configs with G3 power < G1',
       sum(w31$d < 0)),
  row5('h', 'Stage 1 grid (distinct configs)', 'configs with G3 power == G1',
       sum(w31$d == 0)),
  row5('h', 'Stage 1 grid (distinct configs)', 'min (G3 - G1)', min(w31$d)),
  row5('h', 'Stage 1 grid (distinct configs)', 'min (G3 - G1) among alt cells',
       min(w31$d[w31$c_bm > 0])),
  row5('h', 'Stage 1 grid (distinct configs)', 'configs with G3 < G1 among alt cells',
       sum(w31$d[w31$c_bm > 0] < 0)),
  row5('h', 'Stage 1 grid (distinct configs)', 'mean (G3 - G1) among alt cells',
       mean(w31$d[w31$c_bm > 0])),
  w31 |> filter(d < -0.05) |> arrange(d) |>
    transmute(claim = 'h',
              detail = sprintf('%s k=%g %s t1/2=%g %s N=%d cbm=%.2f',
                               carryover_form, weibull_shape, dgp_arch,
                               t1half, design, N, c_bm),
              quantity = 'G3 - G1 where below -0.05', value = d),
  w31 |> filter(c_bm > 0) |>
    group_by(dgp_arch, carryover_form) |>
    summarise(n_below = sum(d < 0), n_cells = n(), .groups = 'drop') |>
    transmute(claim = 'h', detail = sprintf('%s %s (alt cells)', dgp_arch,
                                            carryover_form),
              quantity = sprintf('configs with G3 < G1 of %d', n_cells),
              value = n_below),
  row5('h', 'Weibull k=1 vs exponential twin, all specs',
       'max |power difference|', max(twin$ad)),
  row5('h', 'Weibull k=1 vs exponential twin, all specs',
       'twin pairs with nonzero difference', sum(twin$ad > 0)),
  row5('h', 'Weibull k=1 vs exponential twin, all specs',
       'twin pairs', nrow(twin)))

## (i)
m5$i <- sgf |> filter(sub == 'C2-largeN', N == 280,
                      method %in% c('E2cr2', 'E7cr2', 'F1', 'E1', 'F1A0B')) |>
  transmute(claim = 'i', detail = cell_id,
            quantity = sprintf('%s null rejection (n=%d)', g(method), n),
            value = reject)
emit(bind_rows(m5), 'm5-counts.csv')

## ---------------------------------------------------------------------
## Item 6: size-adjusted comparisons
## ---------------------------------------------------------------------

m6_gf <- map_dfr(names(gf_alt_sets), function(nm) {
  w <- gf_alt_sets[[nm]] |>
    filter(method %in% c('E1', 'E1cr2', 'E2cr2', 'E7cr2')) |>
    select(cell_id, method, reject) |>
    pivot_wider(names_from = method, values_from = reject)
  tibble(part = 'gap-fill alternatives', cell_set = nm, design = 'Hybrid',
         cells = nrow(w),
         mean_G9_minus_G1 = mean(w$E7cr2 - w$E1),
         mean_G9_minus_G6 = mean(w$E7cr2 - w$E1cr2),
         mean_G9_minus_G8 = mean(w$E7cr2 - w$E2cr2),
         cells_G9_below_G6 = sum(w$E7cr2 < w$E1cr2),
         cells_G9_below_G1 = sum(w$E7cr2 < w$E1),
         cells_G9_below_G8 = sum(w$E7cr2 < w$E2cr2))
})
cr2_s3 <- function(d, label) {
  d |> filter(method %in% c('E1', 'E1cr2', 'E2', 'E2cr2', 'E7', 'E7cr2')) |>
    select(cell_id, design, method, reject) |>
    pivot_wider(names_from = method, values_from = reject) |>
    mutate(d16 = E1cr2 - E1, d38 = E2cr2 - E2, d49 = E7cr2 - E7) |>
    group_by(design) |>
    summarise(cells = n(),
              mean_G6_minus_G1 = mean(d16), min_G6_minus_G1 = min(d16),
              mean_G8_minus_G3 = mean(d38), min_G8_minus_G3 = min(d38),
              mean_G9_minus_G4 = mean(d49), min_G9_minus_G4 = min(d49),
              mean_all_three = mean(c(d16, d38, d49)),
              min_all_three = min(c(d16, d38, d49)), .groups = 'drop') |>
    mutate(part = label, cell_set = 'Stage 3 grid')
}
m6_s3 <- bind_rows(cr2_s3(s27 |> filter(c_bm > 0), 'Stage 3 alternatives (24 cells)'),
                   cr2_s3(s27, 'Stage 3 all 30 cells'))
emit(bind_rows(m6_gf, m6_s3), 'm6-size-adjusted.csv')

## ---------------------------------------------------------------------
## Item 7: paired standard errors
## ---------------------------------------------------------------------

p7 <- list()
p7$a <- paired_rej(gf_all |> filter(cell_id == 'hl-t0.2-cbm0.25-N35'),
                   'E7cr2', 'E1') |>
  mutate(source = 'gap-fill job 32', cell = 'hl-t0.2-cbm0.25-N35')
s1cell <- raw1 |> filter(carryover_form == 'exponential', dgp_arch == 'mvn',
                         design == 'Hybrid', N == 70, c_bm == 0.45,
                         t1half == 1) |>
  rename(method = spec)
p7$b <- paired_rej(s1cell, 'E2', 'E1') |>
  mutate(source = 'Stage 1 raw (01-factorial.rds)',
         cell = 'exponential mvn Hybrid N70 cbm0.45 t1/2=1')
p7$b2 <- paired_rej(a27 |> filter(design == 'Hybrid', N == 70, c_bm == 0.45,
                                  t1half == 1), 'E2', 'E1') |>
  mutate(source = 'tier1 (Stage 3) same cell',
         cell = 'panelA Hybrid N70 cbm0.45 t1/2=1')
for (h0 in c(0, 0.5)) {
  p7[[sprintf('b_h%g', h0)]] <- paired_rej(
    raw1 |> filter(carryover_form == 'exponential', dgp_arch == 'mvn',
                   design == 'Hybrid', N == 70, c_bm == 0.45,
                   t1half == h0) |> rename(method = spec), 'E2', 'E1') |>
    mutate(source = 'Stage 1 raw (01-factorial.rds)',
           cell = sprintf('exponential mvn Hybrid N70 cbm0.45 t1/2=%g', h0))
}
alt_ids <- unique(gf_alt_sets[['non-dropout']]$cell_id)
p7$c <- map_dfr(alt_ids, function(id) {
  d <- gf_all |> filter(cell_id == id)
  map_dfr(c('E1', 'E1cr2', 'E2cr2', 'F1'), function(m) {
    if (!any(!is.na(d$p_value[d$method == m]))) return(NULL)
    paired_rej(d, 'E7cr2', m) |>
      mutate(source = sprintf('gap-fill job %d', unique(d$job)), cell = id)
  })
})
m7 <- bind_rows(p7) |>
  select(source, cell, method_a, method_b, n_pairs, rate_a, rate_b, diff,
         paired_se, naive_se, a_only, b_only)
emit(m7, 'm4-paired-se.csv')

## ---------------------------------------------------------------------
## Item 8: Stage 1, G1 vs G2
## ---------------------------------------------------------------------

g12 <- function(d, label) {
  w <- d |> filter(spec %in% c('E1', 'E3')) |>
    select(all_of(cfg_keys), spec, power, bias, mse, coverage) |>
    pivot_wider(names_from = spec, values_from = c(power, bias, mse, coverage))
  tibble(cell_set = label, cells = nrow(w),
         max_abs_diff_bias = max(abs(w$bias_E3 - w$bias_E1)),
         max_abs_diff_mse = max(abs(w$mse_E3 - w$mse_E1)),
         max_abs_diff_coverage = max(abs(w$coverage_E3 - w$coverage_E1)),
         max_abs_diff_power = max(abs(w$power_E3 - w$power_E1)),
         mean_diff_bias = mean(w$bias_E3 - w$bias_E1),
         mean_diff_mse = mean(w$mse_E3 - w$mse_E1),
         mean_diff_coverage = mean(w$coverage_E3 - w$coverage_E1),
         mean_diff_power = mean(w$power_E3 - w$power_E1),
         mean_abs_diff_bias = mean(abs(w$bias_E3 - w$bias_E1)),
         mean_abs_diff_mse = mean(abs(w$mse_E3 - w$mse_E1)),
         mean_abs_diff_coverage = mean(abs(w$coverage_E3 - w$coverage_E1)),
         mean_bias_G1 = mean(w$bias_E1), mean_bias_G2 = mean(w$bias_E3),
         mean_mse_G1 = mean(w$mse_E1), mean_mse_G2 = mean(w$mse_E3),
         mean_cov_G1 = mean(w$coverage_E1), mean_cov_G2 = mean(w$coverage_E3))
}
m8 <- bind_rows(
  g12(grid |> filter(dgp_arch == 'mvn',
                     carryover_form %in% c('exponential', 'weibull')),
      'mvn, exponential + weibull (paper slice)'),
  g12(grid, 'all 540 cells'),
  g12(grid_d, 'all distinct configurations (432)'))
emit(m8, 'm5-stage1-g1-g2.csv')

## ---------------------------------------------------------------------
## Item 9: ML / REML convergence
## ---------------------------------------------------------------------

m9 <- tibble(
  item = 'ML vs REML convergence',
  status = 'not stored',
  note = paste(
    'Per-replicate files record a single logical `converged` per',
    '(cell, rep, method) and a `reason` string in 01-factorial.rds;',
    'the ML comparison fits used by G4/G9 to choose the half-life',
    '(simulation-core.R, AIC by ML then refit by REML) are not',
    'recorded separately. Overall non-convergence rates are in the',
    '`nonconv` column of the summaries.'))
emit(m9, 'm11-reml.csv')

## ---------------------------------------------------------------------
## Item 10: random-effects augmentation, paired per replicate
## ---------------------------------------------------------------------

m10 <- map_dfr(unique(a30$cell_id), function(id) {
  d <- a30 |> filter(cell_id == id)
  map_dfr(c('G9rs', 'G9re', 'G9rse'), function(v) {
    paired_rej(d, v, 'G9') |>
      mutate(cell = id,
             nonconv_variant = mean(!d$converged[d$method == v]),
             nonconv_G9 = mean(!d$converged[d$method == 'G9']))
  })
}) |>
  select(cell, variant = method_a, n_pairs, rate_variant = rate_a,
         rate_G9 = rate_b, diff, paired_se, variant_only = a_only,
         G9_only = b_only, nonconv_variant, nonconv_G9)
m10_sum <- m10 |> group_by(variant) |>
  summarise(cell = 'ALL CELLS: max |diff|', n_pairs = sum(n_pairs),
            diff = max(abs(diff)), paired_se = max(paired_se),
            nonconv_variant = max(nonconv_variant), .groups = 'drop')
emit(bind_rows(m10, m10_sum), 'm12-random-effects.csv')

cat('\nDone. Outputs in', out_dir, '\n')
