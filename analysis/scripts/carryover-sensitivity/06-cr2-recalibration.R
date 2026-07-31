## analysis/scripts/carryover-sensitivity/06-cr2-recalibration.R
##
## Sensitivity block S6 for the carryover-sensitivity manuscript:
## cluster-robust recalibration of the three analysis specifications.
##
## The companion calibration study (analysis/report/
## 10-interaction-test-calibration/) shows that the model-based Wald
## test of the biomarker-by-treatment interaction in the lme + corCAR1
## model is conservative: the model-based standard error is
## over-estimated by 6 to 10 percent and the realised type-I error sits
## near 0.03 rather than 0.05. The conservatism is asymptotic and is
## removed by a CR2 bias-reduced cluster-robust standard error. Because
## the carryover-sensitivity study evaluates that same model-based test,
## its absolute power figures are mildly pessimistic and its reported
## null rejection rate (mean 0.039) is the same conservatism surfacing
## across the grid.
##
## This block re-fits all three specifications under both the
## model-based test and the CR2 cluster-robust test at the reference
## cell and four stress cells (small N, high autocorrelation, half-life
## mis-specification, heavy dropout). It answers two questions the
## conservatism raises:
##   (1) does the recalibration preserve the A2 > A1, A3 ranking, or
##       does it reorder the specifications?
##   (2) what is the correctly-sized power of the recommended analysis
##       (A2), and what does the CR2 degrees-of-freedom penalty cost in
##       the sparse-information cells?
##
## Per cell and specification it records, across n_reps replicates, the
## interaction estimate, the model-based and CR2 standard errors and
## p-values, the CR2 Satterthwaite degrees of freedom, and derives the
## standard-error ratio kappa = mean(SE)/sd(estimate) under each
## inference, the power under each inference, and the A2 - A1 gap.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/06-cr2-recalibration.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(nlme)
  library(clubSandwich)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

args <- commandArgs(trailingOnly = TRUE)
reps_idx <- which(args == '--reps')
N_REPS <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else 500L
SEED <- 20260611L

## --- Stress cells (all Hybrid, so A1 and A3 are well defined; the
## design axis is covered by the main grid). c_bm = 0.45 throughout,
## the reference non-null effect, so the metric is power. -----------

base_cell <- list(
  design = 'Hybrid', N = 70, c_bm = 0.45, t1half = 1.0,
  carryover_form = 'exponential', weibull_shape = 1, rho = 0.7,
  dgp_arch = 'mvn')

cells <- list(
  modifyList(base_cell, list(label = 'Reference (N=70)')),
  modifyList(base_cell, list(label = 'Small N (N=35)', N = 35)),
  modifyList(base_cell, list(label = 'High autocorr (rho=0.9)',
                             rho = 0.9)),
  modifyList(base_cell, list(label = 'Half-life mis-spec (assume 0.25)',
                             analysis_t1half = 0.25)),
  modifyList(base_cell, list(label = '30% MCAR dropout',
                             dropout_rate = 0.30, dropout_mech = 'MCAR')))

## --- Null cells (c_bm = 0) mirroring each stress cell, so the
## rejection rate is the type-I error. A correctly-sized test is the
## precondition for any "best analysis" recommendation: comparing
## power is only meaningful once each specification holds nominal
## alpha. The model-based test is conservative (type-I near 0.03);
## CR2 should restore it to ~0.05. -----------------------------------
null_cells <- lapply(cells, function(cl) modifyList(cl, list(c_bm = 0)))

resp_param     <- default_resp_param()
baseline_param <- default_baseline_param()

## --- one long dataset for a cell (mirrors simulate_cell's data path,
## including the analyst-assumed half-life for the A2 predictor and the
## optional dropout block) --------------------------------------------

gen_cell_long <- function(cell) {
  design_set <- build_design_set(cell$design)
  rho <- default_val(cell[['rho']], 0.7)
  analysis_t1half <- default_val(cell[['analysis_t1half']], cell$t1half)
  analysis_form   <- default_val(cell[['analysis_form']],
                                 cell$carryover_form)
  analysis_shape  <- default_val(cell[['analysis_shape']],
                                 cell$weibull_shape)
  dropout_rate <- default_val(cell[['dropout_rate']], 0)
  dropout_mech <- default_val(cell[['dropout_mech']], 'MCAR')

  model_param <- list(
    N = cell$N, c.bm = cell$c_bm,
    carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    dgp_architecture = cell$dgp_arch,
    c.tv = rho, c.pb = rho, c.br = rho,
    c.cf1t = 0.2, c.cfct = 0.1)

  dat <- tryCatch(
    generate_data_multi_path(model_param, resp_param,
                             baseline_param, design_set),
    error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  if (dropout_rate > 0)
    dat <- apply_dropout(dat, design_set, rate = dropout_rate,
                         mechanism = dropout_mech)
  dl <- prepare_long_data(
    dat, design_set,
    carryover_t1half = cell$t1half,
    carryover_form = cell$carryover_form,
    weibull_shape = cell$weibull_shape,
    analysis_t1half = analysis_t1half,
    analysis_form = analysis_form,
    analysis_shape = analysis_shape)
  dl
}

spec_formula <- function(spec) switch(spec,
  E1 = Sx ~ bm + t + Db  + bm:Db,
  E2 = Sx ~ bm + t + Dbc + bm:Dbc,
  E3 = Sx ~ bm + t + Db  + bm:Db + L)

target_names <- function(spec)
  if (spec == 'E2') c('bm:Dbc', 'Dbc:bm') else c('bm:Db', 'Db:bm')

na_row <- function(spec) tibble(
  spec = spec, estimate = NA_real_, mod_se = NA_real_,
  mod_p = NA_real_, cr2_se = NA_real_, cr2_p = NA_real_,
  cr2_df = NA_real_, converged = FALSE)

## model-based AND CR2 inference for one specification on one dataset.
fit_spec_both <- function(dl, spec) {
  if (spec %in% c('E1', 'E3') && !any(dl$Db == 0)) return(na_row(spec))
  if (spec == 'E3' && !any(dl$L == 1)) return(na_row(spec))
  fit <- tryCatch(
    nlme::lme(spec_formula(spec), random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dl,
              control = nlme::lmeControl(opt = 'optim',
                                         maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
  if (is.null(fit)) return(na_row(spec))
  cc  <- summary(fit)$tTable
  tgt <- intersect(target_names(spec), rownames(cc))
  if (length(tgt) == 0) return(na_row(spec))
  tgt <- tgt[1]

  ct <- tryCatch(
    clubSandwich::coef_test(fit, vcov = 'CR2', cluster = dl$ptID,
                            test = 'Satterthwaite'),
    error = function(e) NULL)
  grab <- function(row, patterns) {
    for (p in patterns) {
      hit <- grep(p, names(row), value = TRUE, ignore.case = TRUE)
      if (length(hit)) return(as.numeric(row[[hit[1]]]))
    }
    NA_real_
  }
  cr2_se <- cr2_p <- cr2_df <- NA_real_
  if (!is.null(ct)) {
    cf  <- if (!is.null(ct$Coef)) as.character(ct$Coef) else rownames(ct)
    idx <- which(cf == tgt)
    if (length(idx)) {
      row    <- ct[idx[1], , drop = FALSE]
      cr2_se <- grab(row, '^SE$')
      cr2_p  <- grab(row, c('^p_Satt', '^p_val', '^p\\.', '^p$', '^p'))
      cr2_df <- grab(row, c('df_Satt', '^df'))
    }
  }
  tibble(
    spec = spec,
    estimate = cc[tgt, 'Value'],
    mod_se   = cc[tgt, 'Std.Error'],
    mod_p    = cc[tgt, 'p-value'],
    cr2_se   = cr2_se, cr2_p = cr2_p, cr2_df = cr2_df,
    converged = TRUE)
}

one_rep <- function(cell, i) {
  dl <- gen_cell_long(cell)
  if (is.null(dl)) return(NULL)
  map_dfr(c('E1', 'E2', 'E3'), function(sp) fit_spec_both(dl, sp)) |>
    mutate(rep = i, .before = 1)
}

run_cell <- function(cell) {
  future_map_dfr(seq_len(N_REPS), function(i) one_rep(cell, i),
                 .options = furrr_options(seed = SEED)) |>
    mutate(cell = cell$label, .before = 1)
}

cat(sprintf('S6 CR2 recalibration: %d cells, %d reps each\n',
            length(cells), N_REPS))
plan(multicore, workers = max(1, parallel::detectCores() - 1))
t0 <- Sys.time()
res      <- map_dfr(cells, run_cell)
res_null <- map_dfr(null_cells, run_cell)
elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
cat(sprintf('Completed in %.0f s.\n\n', elapsed))

summ <- res |>
  filter(converged) |>
  group_by(cell, spec) |>
  summarise(
    n            = n(),
    emp_se       = sd(estimate, na.rm = TRUE),
    mod_se_ratio = mean(mod_se, na.rm = TRUE) / sd(estimate, na.rm = TRUE),
    cr2_se_ratio = mean(cr2_se, na.rm = TRUE) / sd(estimate, na.rm = TRUE),
    power_mod    = mean(mod_p < 0.05, na.rm = TRUE),
    power_cr2    = mean(cr2_p < 0.05, na.rm = TRUE),
    cr2_df       = mean(cr2_df, na.rm = TRUE),
    .groups      = 'drop')

cat('=== S6: model-based vs CR2 across stress cells ===\n')
cat('mod_se_ratio : model SE / empirical SD (>1 = conservative)\n')
cat('cr2_se_ratio : CR2 SE / empirical SD (~1 = calibrated)\n')
cat('power_mod/power_cr2 : rejection rate at c_bm = 0.45\n')
print(as.data.frame(summ), digits = 3)

## A2 - A1 power gap under each inference, per cell
gaps <- summ |>
  select(cell, spec, power_mod, power_cr2) |>
  pivot_wider(names_from = spec,
              values_from = c(power_mod, power_cr2)) |>
  transmute(cell,
            gap_mod = power_mod_A2 - power_mod_A1,
            gap_cr2 = power_cr2_A2 - power_cr2_A1)
cat('\n=== A2 - A1 power gap (ranking check) ===\n')
print(as.data.frame(gaps), digits = 3)

## --- Type-I error at the null (c_bm = 0) under each inference ------
type1 <- res_null |>
  filter(converged) |>
  group_by(cell, spec) |>
  summarise(
    n         = n(),
    type1_mod = mean(mod_p < 0.05, na.rm = TRUE),
    type1_cr2 = mean(cr2_p < 0.05, na.rm = TRUE),
    .groups   = 'drop')
cat('\n=== S6: type-I error at c_bm = 0 (model vs CR2) ===\n')
cat('Nominal alpha = 0.05; model-based is conservative, CR2 ~ 0.05\n')
print(as.data.frame(type1), digits = 3)

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
saveRDS(list(results = res, results_null = res_null,
             summary = summ, gaps = gaps, type1 = type1,
             cells = cells, null_cells = null_cells,
             n_reps = N_REPS, seed = SEED,
             elapsed_secs = elapsed),
        file.path(out_dir, 'diag-s6-cr2.rds'))
cat(sprintf('\nWrote %s\n', file.path(out_dir, 'diag-s6-cr2.rds')))
