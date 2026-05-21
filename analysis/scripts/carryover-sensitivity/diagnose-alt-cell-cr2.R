## analysis/scripts/carryover-sensitivity/diagnose-alt-cell-cr2.R
##
## Alternative-cell robustness check for the type-I examination
## (manuscript 10). Stages 1 to 6 and diagnose-null-type1-cr2.R
## established, at the null (c_bm = 0), that the production model
## M3 over-estimates the interaction standard error and that a CR2
## cluster-robust standard error restores calibration.
##
## This script verifies the one link not yet checked: that at an
## ALTERNATIVE cell (c_bm > 0, where power is measured) the CR2
## recalibration preserves the headline ranking A2 > A1 ~ A3. The
## cell is the manuscript-02 reference cell: Hybrid design,
## Architecture B, N = 70, c_bm = 0.45, exponential DGP,
## t1half = 1.0.
##
## For each replicate the production M3 model is fitted for A1, A2,
## A3, and the interaction is tested both with the model-based Wald
## standard error and with the CR2 cluster-robust standard error
## (clubSandwich, clustered on participant). Power is the rejection
## rate under each.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-alt-cell-cr2.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(ggplot2)
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
  as.integer(args[reps_idx + 1]) else 3000L
SEED <- 20260521L

## manuscript-02 reference cell, at the alternative (c_bm = 0.45)
cell <- list(
  c_bm = 0.45, carryover_form = 'exponential', weibull_shape = 1,
  t1half = 1.0, design = 'Hybrid', N = 70, rho = 0.7,
  dgp_arch = 'mvn')

resp_param     <- default_resp_param()
baseline_param <- default_baseline_param()
design_set     <- build_design_set(cell$design)

model_param <- list(
  N = cell$N, c.bm = cell$c_bm,
  carryover_t1half = cell$t1half,
  carryover_form   = cell$carryover_form,
  weibull_shape    = cell$weibull_shape,
  dgp_architecture = cell$dgp_arch,
  c.tv = cell$rho, c.pb = cell$rho, c.br = cell$rho,
  c.cf1t = 0.2, c.cfct = 0.1)

gen_long <- function() {
  dat <- tryCatch(
    generate_data_multi_path(model_param, resp_param,
                             baseline_param, design_set),
    error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  prepare_long_data(
    dat, design_set,
    carryover_t1half = cell$t1half,
    carryover_form   = cell$carryover_form,
    weibull_shape    = cell$weibull_shape)
}

spec_formula <- function(spec) switch(spec,
  A1 = Sx ~ bm + t + Db  + bm:Db,
  A2 = Sx ~ bm + t + Dbc + bm:Dbc,
  A3 = Sx ~ bm + t + Db  + bm:Db + L)

target_names <- function(spec)
  if (spec == 'A2') c('bm:Dbc', 'Dbc:bm') else c('bm:Db', 'Db:bm')

na_row <- function(spec) tibble(
  spec = spec, estimate = NA_real_, mod_se = NA_real_,
  mod_p = NA_real_, cr2_se = NA_real_, cr2_p = NA_real_,
  cr2_df = NA_real_, converged = FALSE)

cr2_extract <- function(fit, tgt, cluster_vec) {
  ct <- tryCatch(
    clubSandwich::coef_test(fit, vcov = 'CR2', cluster = cluster_vec,
                            test = 'Satterthwaite'),
    error = function(e) NULL)
  if (is.null(ct)) return(list(se = NA_real_, p = NA_real_,
                               df = NA_real_))
  cf <- if (!is.null(ct$Coef)) as.character(ct$Coef) else rownames(ct)
  idx <- which(cf == tgt)
  if (length(idx) == 0) return(list(se = NA_real_, p = NA_real_,
                                    df = NA_real_))
  row <- ct[idx[1], , drop = FALSE]
  grab <- function(patterns) {
    for (p in patterns) {
      hit <- grep(p, names(row), value = TRUE, ignore.case = TRUE)
      if (length(hit)) return(as.numeric(row[[hit[1]]]))
    }
    NA_real_
  }
  list(se = grab('^SE$'),
       p  = grab(c('^p_Satt', '^p_val', '^p\\.', '^p$', '^p')),
       df = grab(c('df_Satt', '^df')))
}

fit_M3_cr2 <- function(dl, spec) {
  fit <- tryCatch(
    nlme::lme(spec_formula(spec), random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dl,
              control = nlme::lmeControl(
                opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
  if (is.null(fit)) return(na_row(spec))
  cc <- summary(fit)$tTable
  tgt <- intersect(target_names(spec), rownames(cc))
  if (length(tgt) == 0) return(na_row(spec))
  tgt <- tgt[1]
  cr2 <- cr2_extract(fit, tgt, dl$ptID)
  tibble(
    spec = spec, estimate = cc[tgt, 'Value'],
    mod_se = cc[tgt, 'Std.Error'], mod_p = cc[tgt, 'p-value'],
    cr2_se = cr2$se, cr2_p = cr2$p, cr2_df = cr2$df,
    converged = TRUE)
}

one_rep <- function(i) {
  dl <- gen_long()
  if (is.null(dl)) return(NULL)
  map_dfr(c('A1', 'A2', 'A3'),
          function(sp) fit_M3_cr2(dl, sp)) |>
    mutate(rep = i, .before = 1)
}

cat(sprintf('Alternative-cell CR2 check: %s, Arch B, N=%d, c_bm=%.2f,\n',
            cell$design, cell$N, cell$c_bm))
cat(sprintf('  exponential DGP, t1half=%.1f, %d replicates\n',
            cell$t1half, N_REPS))

plan(multicore, workers = max(1, parallel::detectCores() - 1))
t0 <- Sys.time()
res <- future_map_dfr(seq_len(N_REPS), one_rep,
                      .options = furrr_options(seed = SEED))
elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
cat(sprintf('Completed in %.0f s.\n\n', elapsed))

summ <- res |>
  filter(converged, !is.na(cr2_p)) |>
  group_by(spec) |>
  summarise(
    n            = n(),
    emp_se       = sd(estimate),
    mod_se_ratio = mean(mod_se) / sd(estimate),
    cr2_se_ratio = mean(cr2_se) / sd(estimate),
    power_mod    = mean(mod_p < 0.05),
    power_cr2    = mean(cr2_p < 0.05),
    mc_se_mod    = sqrt(mean(mod_p < 0.05) *
                        (1 - mean(mod_p < 0.05)) / n()),
    mc_se_cr2    = sqrt(mean(cr2_p < 0.05) *
                        (1 - mean(cr2_p < 0.05)) / n()),
    .groups = 'drop')

cat('=== Power at the alternative cell: model-based vs CR2 ===\n')
cat('mod_se_ratio : production SE / empirical SE (expect ~1.05-1.10)\n')
cat('cr2_se_ratio : CR2 robust SE / empirical SE (expect ~1.00)\n')
cat('power_mod    : rejection rate, production Wald test\n')
cat('power_cr2    : rejection rate, CR2 cluster-robust test\n\n')
print(as.data.frame(summ), digits = 4)

pull1 <- function(s, col) summ[[col]][summ$spec == s]
gap <- function(col, hi, lo) pull1(hi, col) - pull1(lo, col)

cat('\n=== Ranking check ===\n')
cat(sprintf('Model-based:  A2=%.3f  A1=%.3f  A3=%.3f\n',
            pull1('A2','power_mod'), pull1('A1','power_mod'),
            pull1('A3','power_mod')))
cat(sprintf('CR2 robust :  A2=%.3f  A1=%.3f  A3=%.3f\n',
            pull1('A2','power_cr2'), pull1('A1','power_cr2'),
            pull1('A3','power_cr2')))
cat(sprintf('A2 - A1 advantage: model-based %+.3f, CR2 %+.3f\n',
            gap('power_mod','A2','A1'), gap('power_cr2','A2','A1')))
cat(sprintf('A2 - A3 advantage: model-based %+.3f, CR2 %+.3f\n',
            gap('power_mod','A2','A3'), gap('power_cr2','A2','A3')))

ranking_holds <- pull1('A2','power_cr2') >= pull1('A1','power_cr2') &&
                 pull1('A2','power_cr2') >= pull1('A3','power_cr2')
cat(sprintf('\nVERDICT: under CR2 the A2 >= A1, A2 >= A3 ranking %s\n',
    if (ranking_holds) 'HOLDS.' else 'does NOT hold - inspect above.'))

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')

theme_paper <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'top')

p_fig <- summ |>
  select(spec, power_mod, power_cr2) |>
  pivot_longer(c(power_mod, power_cr2), names_to = 'method',
               values_to = 'power') |>
  mutate(method = recode(method,
    power_mod = 'Model-based Wald', power_cr2 = 'CR2 cluster-robust')) |>
  ggplot(aes(spec, power, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf('%.3f', power)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = c(
    'Model-based Wald' = '#e31a1c',
    'CR2 cluster-robust' = '#1f78b4')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = 'Analysis specification', y = 'Power', fill = NULL,
       title = 'Power at the reference alternative cell',
       subtitle = sprintf(
         'Arch B, Hybrid, N=70, c_bm=0.45, exp DGP, %d reps',
         N_REPS)) +
  theme_paper

ggsave(file.path(out_dir, 'diag-alt-cell-cr2.pdf'),
       p_fig, width = 6, height = 3.8)

saveRDS(list(results = res, summary = summ, cell = cell,
             n_reps = N_REPS, seed = SEED, elapsed_secs = elapsed),
        file.path(out_dir, 'diag-alt-cell-cr2.rds'))

cat(sprintf('\nWrote %s\n',
    file.path(out_dir, 'diag-alt-cell-cr2.pdf')))
cat(sprintf('Wrote %s\n',
    file.path(out_dir, 'diag-alt-cell-cr2.rds')))
