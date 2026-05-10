## Study 5 production-results analysis: Phase B kill-switch
## evaluation and Morris-format reporting tables.
##
## Consumes the per-replicate output from 10-study5-production.R
## and produces:
##   - Table 5.1: Type I error by N and test (Gaussian null cell)
##   - Table 5.2: Class-count selection, single-component DGPs
##   - Table 5.3: Class-count selection, true two-class DGP
##   - Table 5.4: Convergence and timing summary
##   - Phase B kill-switch evaluation against pre-registered
##     thresholds.
##
## Output: analysis/data/study5/study5-tables.rds (list of tables)
##         analysis/tables/study5-table-{N}.csv (one CSV per table)
##         analysis/data/study5/study5-killswitch.txt
##
## Run from the package root:
##   Rscript analysis/scripts/latent-class-mixture-application/11-study5-analysis.R

suppressPackageStartupMessages({
  library(data.table)
})

DATA_DIR <- file.path('analysis', 'data', 'study5')
TBL_DIR  <- file.path('analysis', 'tables')
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

results <- readRDS(file.path(DATA_DIR, 'study5-results.rds'))

cat(sprintf('Study 5 results loaded: %d rows, %d cells\n',
            nrow(results),
            length(unique(results$cell_id))))

binom_mcse <- function(rate, n) {
  ifelse(is.na(rate) | is.na(n) | n == 0, NA_real_,
         sqrt(rate * (1 - rate) / n))
}

## ---- Table 5.1: Type I error under Gaussian null ----
##
## Pre-registered prediction: nominal Type I at alpha = 0.05 for
## the BIC-conditional gating Wald and the lme bm:Dbc Wald. The
## unconditional gating Wald is on the boundary of the parameter
## space and is documented as miscalibrated; we report it for
## completeness.
table_5_1 <- results[
  dgp_family == 'gaussian_null',
  {
    n_lme  <- sum(!is.na(lme_p))
    n_lcmm <- sum(!is.na(lcmm_gating_p))
    rej_lme <- mean(lme_p < 0.05, na.rm = TRUE)
    rej_lcmm_uncond <- mean(lcmm_gating_p < 0.05, na.rm = TRUE)
    rej_lcmm_within <- mean(lcmm_within_p < 0.05, na.rm = TRUE)
    rej_lcmm_cond <- mean(
      (lcmm_delta_bic > 0) & (lcmm_gating_p < 0.05),
      na.rm = TRUE
    )
    list(
      n_reps = n_lcmm,
      conv_lme = mean(!is.na(lme_p)),
      conv_lcmm = mean(lcmm_conv == 1, na.rm = TRUE),
      rej_lme = rej_lme,
      mcse_lme = binom_mcse(rej_lme, n_lme),
      rej_lcmm_uncond = rej_lcmm_uncond,
      mcse_lcmm_uncond = binom_mcse(rej_lcmm_uncond, n_lcmm),
      rej_lcmm_within = rej_lcmm_within,
      mcse_lcmm_within = binom_mcse(rej_lcmm_within, n_lcmm),
      rej_lcmm_cond = rej_lcmm_cond,
      mcse_lcmm_cond = binom_mcse(rej_lcmm_cond, n_lcmm)
    )
  },
  by = .(N)
]
setorder(table_5_1, N)

## ---- Table 5.2: Class-count selection, single-component DGPs ----
##
## Three DGPs: Gaussian, t5, skew-normal. Truth in all three is
## K = 1. Pre-registered prediction: BIC selects K = 1 in
## >= 80% of replicates at N >= 70.
table_5_2 <- results[
  dgp_family %in% c('gaussian_null', 't5_null', 'skewnorm4_null'),
  {
    n_lcmm <- sum(!is.na(lcmm_delta_bic))
    p_k1 <- mean(lcmm_delta_bic <= 0, na.rm = TRUE)
    p_k2 <- mean(lcmm_delta_bic > 0, na.rm = TRUE)
    p_k2_strong <- mean(lcmm_delta_bic > 6, na.rm = TRUE)
    list(
      n_reps = n_lcmm,
      conv_lcmm = mean(lcmm_conv == 1, na.rm = TRUE),
      p_select_k1 = p_k1,
      p_select_k2 = p_k2,
      p_select_k2_strong = p_k2_strong,
      mcse_k1 = binom_mcse(p_k1, n_lcmm)
    )
  },
  by = .(dgp_family, N)
]
setorder(table_5_2, dgp_family, N)

## ---- Table 5.3: Class-count selection, true two-class ----
##
## Class separation Delta in {0, 0.5, 1.0, 1.5, 2.0}. Truth is
## K = 2 (well, K = 1 at Delta = 0 since there is no class
## structure left). Pre-registered prediction: P(BIC selects K=2)
## increases with Delta.
table_5_3 <- results[
  dgp_family == 'two_class',
  {
    n_lcmm <- sum(!is.na(lcmm_delta_bic))
    p_k2 <- mean(lcmm_delta_bic > 0, na.rm = TRUE)
    p_k2_strong <- mean(lcmm_delta_bic > 6, na.rm = TRUE)
    list(
      n_reps = n_lcmm,
      conv_lcmm = mean(lcmm_conv == 1, na.rm = TRUE),
      p_select_k2 = p_k2,
      p_select_k2_strong = p_k2_strong,
      mcse_k2 = binom_mcse(p_k2, n_lcmm)
    )
  },
  by = .(delta, N)
]
setorder(table_5_3, delta, N)

## ---- Table 5.4: Convergence and timing ----
table_5_4 <- results[, {
  n_total <- .N
  list(
    n_reps = n_total,
    conv_lme = mean(!is.na(lme_p)),
    conv_lcmm = mean(lcmm_conv == 1, na.rm = TRUE),
    mean_lme_t = mean(lme_t, na.rm = TRUE),
    mean_lcmm_t = mean(lcmm_t, na.rm = TRUE),
    median_lcmm_t = median(lcmm_t, na.rm = TRUE),
    p95_lcmm_t = quantile(lcmm_t, 0.95, na.rm = TRUE)
  )
}, by = .(dgp_family, delta, N)]
setorder(table_5_4, dgp_family, delta, N)

## ---- Save tables ----
saveRDS(list(table_5_1 = table_5_1, table_5_2 = table_5_2,
              table_5_3 = table_5_3, table_5_4 = table_5_4),
        file.path(DATA_DIR, 'study5-tables.rds'))
fwrite(table_5_1, file.path(TBL_DIR, 'study5-table-5_1.csv'))
fwrite(table_5_2, file.path(TBL_DIR, 'study5-table-5_2.csv'))
fwrite(table_5_3, file.path(TBL_DIR, 'study5-table-5_3.csv'))
fwrite(table_5_4, file.path(TBL_DIR, 'study5-table-5_4.csv'))

cat('\n=== Table 5.1: Type I error under Gaussian null ===\n')
print(table_5_1)
cat('\n=== Table 5.2: Class-count selection, single-component ===\n')
print(table_5_2)
cat('\n=== Table 5.3: Class-count selection, true two-class ===\n')
print(table_5_3)
cat('\n=== Table 5.4: Convergence and timing ===\n')
print(table_5_4)

## ---- Phase B kill-switch evaluation ----
##
## Per simgo-paper03-phasing.md:
##   Criterion 1. Type I error within +/- 0.01 of nominal at
##                N >= 35.
##   Criterion 2. Convergence rate >= 90% in cells 1 and 2.
##   Criterion 3. Class-count selection correctly identifies
##                K = 1 at single-component DGPs in >= 80% of
##                replicates at N >= 70.
##
## The Type I criterion is evaluated on the BIC-conditional
## procedure (procedure of record per the pre-registration); the
## unconditional gating Wald is reported for completeness but is
## known to be miscalibrated under boundary nulls.
killswitch_eval <- function(table_5_1, table_5_2) {
  c1_pass <- all(
    abs(table_5_1$rej_lcmm_cond - 0.05) <= 0.01 &
    abs(table_5_1$rej_lme - 0.05) <= 0.01,
    na.rm = TRUE
  )
  c2_cells <- table_5_4[
    dgp_family %in% c('gaussian_null', 't5_null', 'skewnorm4_null') &
    N >= 35
  ]
  c2_pass <- all(c2_cells$conv_lcmm >= 0.90, na.rm = TRUE)
  c3_cells <- table_5_2[N >= 70]
  c3_pass <- all(c3_cells$p_select_k1 >= 0.80, na.rm = TRUE)
  list(
    type_I_pass = c1_pass,
    convergence_pass = c2_pass,
    k_selection_pass = c3_pass,
    overall_pass = c1_pass && c2_pass && c3_pass
  )
}

ks <- killswitch_eval(table_5_1, table_5_2)
ks_path <- file.path(DATA_DIR, 'study5-killswitch.txt')
writeLines(c(
  'Phase B kill-switch evaluation',
  paste0('Date: ', format(Sys.time(), '%Y-%m-%d %H:%M %Z')),
  '',
  sprintf('Criterion 1 (Type I within +/- 0.01 nominal): %s',
          ifelse(ks$type_I_pass, 'PASS', 'FAIL')),
  sprintf('Criterion 2 (Convergence >= 90%% at N >= 35): %s',
          ifelse(ks$convergence_pass, 'PASS', 'FAIL')),
  sprintf('Criterion 3 (K=1 selection >= 80%% at N >= 70): %s',
          ifelse(ks$k_selection_pass, 'PASS', 'FAIL')),
  '',
  sprintf('Overall Phase B decision: %s',
          ifelse(ks$overall_pass,
                 'PASS - proceed to Phase C',
                 'FAIL - resubmit paper 03 as path (a)'))
), ks_path)
cat(sprintf('\nKill-switch evaluation -> %s\n', ks_path))
cat(readLines(ks_path), sep = '\n')
