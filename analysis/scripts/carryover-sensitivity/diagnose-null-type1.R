## analysis/scripts/carryover-sensitivity/diagnose-null-type1.R
##
## Paper 10 (analysis/report/10-interaction-test-calibration/),
## Section "The conservatism is standard-error over-estimation".
## Stages 1 to 3 of type1-examination-plan.md: the first-order
## decomposition of the null rejection rate of the biomarker-by-drug
## interaction test under the production model M3 (lme, random
## intercept plus corCAR1).
##
## Stage 1  Per-replicate output: estimate, model SE, lme degrees of
##          freedom, corCAR1 rho_a, random-intercept and residual SDs.
## Stage 2  Bias, kappa, SD of the statistic, CV of the SE, size.
## Stage 3  Normal-reference p-values (degrees-of-freedom channel).
## Oracle   The same M3 fit to Sx minus its true path-by-visit mean
##          (correct mean model), on the same datasets. Isolates the
##          share of the excess kappa due to mean-model misfit.
##
## Corrections of 2026-10-29 (referee report 2026-10-28):
##   - refined DGP (paper10-common.R: c.cfct = 0.2, carryover decaying
##     from the last on-drug mean, stop_on_nonpd = TRUE);
##   - no fallback refit without corCAR1: a failed M3 fit is recorded
##     as not converged (the earlier version counted 32 to 40 random-
##     intercept-only fits as M3).
##
## Primary null cell: c_bm = 0, exponential decay, Hybrid design,
## N = 70, t1half = 1.0 week, rho = 0.7. The seed is shared with the
## stage-4, CR2, GEE and KR drivers, so replicate i is the same
## dataset in every one of them.
##
## Usage (thread variables set to 1, see run-paper10-refined.sh):
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-null-type1.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(nlme)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/paper10-common.R'))

args <- commandArgs(trailingOnly = TRUE)
reps_idx <- which(args == '--reps')
N_REPS <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else 5000L
SEED <- 20260520L

cell <- list(
  c_bm = 0, carryover_form = 'exponential', weibull_shape = 1,
  t1half = 1.0, design = 'Hybrid', N = 70, rho = 0.7,
  dgp_arch = 'mvn')

design_set  <- build_design_set(cell$design)
model_param <- p10_model_param(cell$N, cell$c_bm, cell$t1half, cell$rho)
true_means  <- p10_true_means(model_param, design_set)

## Every path builds under stop_on_nonpd = TRUE; record the minimum
## eigenvalue of each.
pd_min_eig <- map_dbl(design_set, function(td) {
  s <- build_sigma_matrix(model_param, default_resp_param(),
                          default_baseline_param(), td,
                          dgp_architecture = cell$dgp_arch,
                          stop_on_nonpd = TRUE)$sigma
  pos <- diag(s) > 0
  min(eigen(s[pos, pos], only.values = TRUE)$values)
})

na_spec_row <- function(spec, mean_model) {
  tibble(spec = spec, mean_model = mean_model, estimate = NA_real_,
         std_error = NA_real_, df = NA_real_, t_value = NA_real_,
         p_value = NA_real_, rho_a = NA_real_, sd_ranef = NA_real_,
         sd_resid = NA_real_, converged = FALSE)
}

fit_spec_diag <- function(dat_long, spec, mean_model = 'fitted') {
  dl <- dat_long
  if (mean_model == 'oracle') dl$Sx <- dl$Sx - dl$mu_true
  fit <- p10_lme(p10_spec_formula(spec), dl,
                 correlation = nlme::corCAR1(form = ~t | ptID))
  if (is.null(fit)) return(na_spec_row(spec, mean_model))
  cc <- summary(fit)$tTable
  target <- intersect(p10_target(spec), rownames(cc))
  if (length(target) == 0) return(na_spec_row(spec, mean_model))
  rho_a <- as.numeric(coef(fit$modelStruct$corStruct,
                           unconstrained = FALSE))[1]
  sd_ranef <- suppressWarnings(
    as.numeric(nlme::VarCorr(fit)['(Intercept)', 'StdDev']))
  tibble(spec = spec, mean_model = mean_model,
         estimate  = cc[target[1], 'Value'],
         std_error = cc[target[1], 'Std.Error'],
         df        = cc[target[1], 'DF'],
         t_value   = cc[target[1], 't-value'],
         p_value   = cc[target[1], 'p-value'],
         rho_a = rho_a, sd_ranef = sd_ranef, sd_resid = fit$sigma,
         converged = TRUE)
}

grid <- expand_grid(spec = c('E1', 'E2', 'E3'),
                    mean_model = c('fitted', 'oracle'))

one_rep <- function(i) {
  dl <- p10_gen_long(model_param, design_set, cell$t1half)
  if (is.null(dl)) {
    return(pmap_dfr(grid, na_spec_row) |>
             mutate(rep = i, gen_ok = FALSE, .before = 1))
  }
  dl$path_i <- as.integer(sub('_.*', '', dl$ptID))
  dl <- dl |> left_join(true_means,
                        by = c(path_i = 'path', 'timepoint_name'))
  pmap_dfr(grid, function(spec, mean_model)
    fit_spec_diag(dl, spec, mean_model)) |>
    mutate(rep = i, gen_ok = TRUE, .before = 1)
}

cat(sprintf('Null cell: %s, N=%d, c_bm=0, t1half=%.1f, c.cfct=%.1f, %d reps\n',
            cell$design, cell$N, cell$t1half, P10_CFCT, N_REPS))

plan(multicore, workers = max(1, parallelly::availableCores() - 1))
t0 <- Sys.time()
res <- future_map_dfr(seq_len(N_REPS), one_rep,
                      .options = furrr_options(seed = SEED))
elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
cat(sprintf('Completed in %.0f s.\n\n', elapsed))

res <- res |> mutate(p_value_z = 2 * pnorm(-abs(t_value)))

decomp <- res |>
  filter(converged) |>
  group_by(spec, mean_model) |>
  summarise(
    n_conv    = n(),
    mean_est  = mean(estimate),
    mcse_mean = sd(estimate) / sqrt(n()),
    emp_se    = sd(estimate),
    mod_se    = mean(std_error),
    se_ratio  = mean(std_error) / sd(estimate),
    cv_se     = sd(std_error) / mean(std_error),
    sd_tstat  = sd(t_value),
    type1_01  = mean(p_value < 0.01),
    type1_05  = mean(p_value < 0.05),
    type1_10  = mean(p_value < 0.10),
    type1_05_z = mean(p_value_z < 0.05),
    mean_df   = mean(df),
    ri_boundary = mean(sd_ranef < 1),
    .groups = 'drop') |>
  mutate(mcse_05 = p10_mcse(type1_05, n_conv))

conv <- res |>
  group_by(spec, mean_model) |>
  summarise(n = n(), n_conv = sum(converged), gen_fail = sum(!gen_ok),
            .groups = 'drop')

print(as.data.frame(conv))
print(as.data.frame(decomp), digits = 4)

out_dir <- p10_out_dir(repo_root)
saveRDS(
  list(results = res, decomp = decomp, conv = conv, cell = cell,
       n_reps = N_REPS, seed = SEED, pd_min_eig = pd_min_eig,
       true_means = true_means, dgp = p10_dgp_meta(repo_root),
       elapsed_secs = elapsed),
  file.path(out_dir, 'diag-null-type1.rds'))
cat(sprintf('Wrote %s\n', file.path(out_dir, 'diag-null-type1.rds')))
