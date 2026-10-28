## analysis/scripts/carryover-sensitivity/35-check-new-variants.R
##
## Local checks for the additions made for paper 02's referee report
## of 2026-10-27 (M1, M7, M9, M12):
##
##   1. Backward compatibility, bit for bit: with b_tv_trend = 0 and
##      sd_slope = 0 the current generator and driver core reproduce,
##      for the same seed, the output of the versions committed at
##      HEAD (fetched with `git show`, sourced into a separate
##      environment). Covers generate_data(), simulate_cell_combined()
##      and the G9 selection path (select_half_life_aic defaults).
##   2. The new specifications (E7null_cr2, E3fix_cr2 at 0.25 and 1
##      week, F1 via fit_f1_cr2) return estimates on 4 replicates of
##      one Hybrid cell at N = 35, c_bm = 0.45, t1/2 = 0.5.
##   3. b_tv_trend > 0 shifts the TV component in proportion to bm_z
##      and to the normalized TV time course, so the shift grows with
##      time; the constant b_tv shift does not.
##   4. sd_slope > 0 adds a patient slope whose contribution to the
##      outcome variance grows as t^2; sd_slope = 0 is unchanged.
##
## Usage (from the repository root):
##   Rscript analysis/scripts/carryover-sensitivity/35-check-new-variants.R
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

## ---------------- HEAD versions in a separate environment ----------------
old_dir <- file.path(tempdir(), 'head-versions')
dir.create(old_dir, showWarnings = FALSE)
for (f in c('implementations/tidyverse/R/functions.R',
            'analysis/scripts/carryover-sensitivity/simulation-core.R')) {
  writeLines(system2('git', c('show', shQuote(paste0('HEAD:', f))),
                     stdout = TRUE),
             file.path(old_dir, basename(f)))
}
old <- new.env()
sys.source(file.path(old_dir, 'functions.R'), envir = old)
sys.source(file.path(old_dir, 'simulation-core.R'), envir = old)
old_generate_data <- get('generate_data', envir = old)
old_simulate_cell_combined <- get('simulate_cell_combined', envir = old)
old_fit_spec_s7 <- get('fit_spec_s7', envir = old)

cell <- list(design = 'Hybrid', N = 35, t1half = 0.5,
             carryover_form = 'exponential', weibull_shape = 1,
             dgp_arch = 'mvn', c_bm = 0.45)

cat('\n== 1. Backward compatibility against HEAD ==\n')
design_set <- build_design_set('Hybrid')
mp <- list(N = 9, carryover_t1half = 0.5, carryover_form = 'exponential',
           weibull_shape = 1, dgp_architecture = 'mvn', c.bm = 0.45,
           c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.2,
           c.bm.bl = 0, c.bm.tv = 0, b.bm.tv = 0, c.bm.pb = 0,
           b.bm.tv.trend = 0, sd_slope = 0)
set.seed(3101)
d_new <- generate_data(mp, default_resp_param(), default_baseline_param(),
                       design_set[[1]], empirical = FALSE,
                       make_positive_definite = TRUE)
set.seed(3101)
d_old <- old_generate_data(mp, default_resp_param(), default_baseline_param(),
                           design_set[[1]], empirical = FALSE,
                           make_positive_definite = TRUE)
cat('generate_data identical to HEAD (trend = 0, slope = 0):',
    identical(d_new, d_old), '\n')

set.seed(3102)
r_new <- simulate_cell_combined(cell, n_reps = 1, n_boot = 0)
set.seed(3102)
r_old <- old_simulate_cell_combined(cell, n_reps = 1, n_boot = 0)
cat('simulate_cell_combined identical to HEAD (one replicate, 13 methods):',
    identical(r_new, r_old), '\n')
print(as.data.frame(r_new |> filter(method %in% c('E2cr2', 'E7cr2', 'F1')) |>
  select(method, h, estimate, std_error, p_value)), row.names = FALSE)

## select_half_life_aic default path on a shared trial
set.seed(3103)
dat <- generate_data_multi_path(c(mp, list()) |> modifyList(list(N = 35)),
                                default_resp_param(), default_baseline_param(),
                                design_set)
dat_long <- prepare_long_data(dat, design_set, carryover_t1half = 0.5)
cat('fit_spec_s7(E7cr2) identical to HEAD on a shared trial:',
    identical(fit_spec_s7(dat_long, 'E7cr2'),
              old_fit_spec_s7(dat_long, 'E7cr2')), '\n')

cat('\n== 2. New specifications on 4 replicates (Hybrid, N = 35, c_bm = 0.45, t1/2 = 0.5) ==\n')
set.seed(3104)
t0 <- Sys.time()
ref <- simulate_cell_referee(cell, n_reps = 4)
cat(sprintf('4 replicates in %.1f s\n',
            as.numeric(Sys.time() - t0, units = 'secs')))
print(as.data.frame(ref |>
  mutate(across(c(estimate, std_error, p_value, t), ~ round(.x, 4)))),
  row.names = FALSE)
stopifnot(all(ref$converged), !anyNA(ref$estimate))
cat('null-model selection agrees with interaction-model selection in',
    sum((ref |> filter(method == 'E7cr2'))$h ==
        (ref |> filter(method == 'E7null_cr2'))$h), 'of 4 replicates\n')

cat('\n== 3. b_tv_trend shifts TV progressively; b_tv does not ==\n')
bp <- default_baseline_param(); rp <- default_resp_param()
gen <- function(extra) {
  set.seed(3105)
  generate_data(modifyList(mp, c(list(N = 400), extra)), rp, bp,
                design_set[[1]], empirical = FALSE,
                make_positive_definite = TRUE)
}
d0 <- gen(list())
d_tr <- gen(list(b.bm.tv.trend = 0.3))
d_ct <- gen(list(b.bm.tv = 0.3))
tv_cols <- paste(design_set[[1]]$timepoint_name, 'tv', sep = '.')
week <- cumsum(design_set[[1]]$t_wk)
sig <- build_sigma_matrix(mp, rp, bp, design_set[[1]])
g_t <- sig$means[match(tv_cols, sig$labels)]
g_t <- g_t / max(abs(g_t))
bm_z <- (d0$bm - bp$m[bp$cat == 'bm']) / bp$sd[bp$cat == 'bm']
slope_per_sd <- function(d) {
  map_dbl(tv_cols, ~ coef(lm(I(d[[.x]] - d0[[.x]]) ~ bm_z))[2])
}
cat('bm identical across the three draws:',
    identical(d0$bm, d_tr$bm) && identical(d0$bm, d_ct$bm), '\n')
print(tibble(timepoint = design_set[[1]]$timepoint_name, week = week,
             g_t = round(g_t, 4),
             expected_trend = round(0.3 * 6 * g_t, 4),
             observed_trend = round(slope_per_sd(d_tr), 4),
             observed_const = round(slope_per_sd(d_ct), 4)) |>
        as.data.frame(), row.names = FALSE)
cat('observed trend shift equals 0.3 * sigma_TV * g(t) to 1e-10:',
    isTRUE(all.equal(slope_per_sd(d_tr), 0.3 * 6 * g_t, tolerance = 1e-10)),
    '\n')
cat('trend shift is non-decreasing in time:',
    all(diff(slope_per_sd(d_tr)) >= 0), '\n')
cat('BR columns untouched by the trend:',
    identical(d0$OL1.br, d_tr$OL1.br), '\n')

cat('\n== 4. sd_slope adds a patient-specific linear trend ==\n')
d_sl <- gen(list(sd_slope = 0.3))
d_sl0 <- gen(list(sd_slope = 0))
cat('sd_slope = 0 identical to the default:', identical(d_sl0, d0), '\n')
delta <- sapply(tv_cols, function(k) d_sl[[k]] - d0[[k]])
u_hat <- delta[, 8] / week[8]
print(tibble(timepoint = design_set[[1]]$timepoint_name, week = week,
             sd_delta = round(apply(delta, 2, sd), 3),
             expected = round(0.3 * week, 3)) |> as.data.frame(),
      row.names = FALSE)
cat('delta is exactly u_i * t at every timepoint:',
    isTRUE(all.equal(delta, outer(u_hat, week), tolerance = 1e-10,
                     check.attributes = FALSE)), '\n')
cat(sprintf('empirical SD of the recovered slopes: %.3f (target 0.3)\n',
            sd(u_hat)))
cat('sd_slope changes the outcome column Sx:',
    !identical(d0$OL1, d_sl$OL1), '\n')

cat('\n== 5. simulate_cell_referee with trend and slope cells runs ==\n')
set.seed(3106)
r_tr <- simulate_cell_referee(modifyList(cell, list(c_bm = 0, N = 70,
                                                    b_tv_trend = 0.3)),
                              n_reps = 1,
                              specs = c('E2cr2', 'E7cr2', 'F1'))
set.seed(3106)
r_sl <- simulate_cell_referee(modifyList(cell, list(N = 70, sd_slope = 0.3)),
                              n_reps = 1,
                              specs = c('E2cr2', 'E7cr2', 'F1'))
print(as.data.frame(bind_rows(mutate(r_tr, cell = 'trend-0.3'),
                              mutate(r_sl, cell = 'slope-0.3')) |>
  mutate(across(c(estimate, std_error, p_value, t), ~ round(.x, 4)))),
  row.names = FALSE)
cat('\nAll checks ran.\n')
