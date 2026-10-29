## analysis/scripts/carryover-sensitivity/diagnose-null-type1-gee.R
##
## Paper 10 (analysis/report/10-interaction-test-calibration/),
## Section "The marginal-model alternative". Fits the GEE counterpart
## of the analysis model (geepack::geeglm, Gaussian, identity link) at
## the primary null cell and tests the interaction with
##   (a) the Liang-Zeger sandwich SE reported by geeglm, and
##   (b) the Mancl-DeRouen bias-corrected sandwich, computed inline
##       from the same geeglm fit (estimate and SE from one fit),
## each referred both to a normal reference and to t on K - p df
## (K clusters, p fixed effects), so the SE estimator and the
## reference can be compared one at a time.
##
## Working correlation. geepack's corstr = 'ar1' indexes the lag by
## visit order within cluster, alpha^|j - k|, and ignores the spacing
## of the visits (weeks 4, 8, 9, 10, 11, 12, 16, 20) even when waves
## are supplied (checked 2026-10-29: the geeglm SE equals an inline
## sandwich with order lags, not with calendar lags). The inline
## Mancl-DeRouen sandwich uses the same order-lag matrix so that both
## SEs refer to one working correlation. A calendar-time working
## correlation is not available in geepack and is not evaluated.
##
## Corrections of 2026-10-29 (referee report 2026-10-28): refined DGP
## (paper10-common.R); estimate and SE from one fit (the earlier MD SE
## came from a separate geesmv refit); both references reported for
## both SEs.
##
## The seed matches diagnose-null-type1.R, so the datasets are
## identical and the GEE estimate can be compared with the lme M3
## estimate replicate by replicate (relative efficiency).
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-null-type1-gee.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(geepack)
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

na_row <- function(spec) tibble(
  spec = spec, estimate = NA_real_, alpha = NA_real_,
  lz_se = NA_real_, md_se = NA_real_, df = NA_real_,
  converged = FALSE)

## Mancl-DeRouen (2001) sandwich for a Gaussian identity-link GEE with
## working correlation alpha^|j - k| in visit order (rows sorted by
## time within cluster). The scale cancels.
md_vcov <- function(X, y, beta, id, alpha) {
  res <- as.numeric(y - X %*% beta)
  ids <- unique(id)
  blocks <- lapply(ids, function(k) {
    rows <- which(id == k)
    j <- seq_along(rows)
    R <- alpha^abs(outer(j, j, '-'))
    list(X = X[rows, , drop = FALSE], r = res[rows], Ri = solve(R))
  })
  bread <- Reduce(`+`, lapply(blocks, function(b)
    crossprod(b$X, b$Ri) %*% b$X))
  bread_inv <- solve(bread)
  meat <- Reduce(`+`, lapply(blocks, function(b) {
    H <- b$X %*% bread_inv %*% crossprod(b$X, b$Ri)
    u <- crossprod(b$X, b$Ri) %*% solve(diag(nrow(H)) - H, b$r)
    tcrossprod(u)
  }))
  bread_inv %*% meat %*% bread_inv
}

fit_gee <- function(dl, spec) {
  form <- p10_spec_formula(spec)
  fit <- tryCatch(
    geepack::geeglm(form, family = gaussian, data = dl, id = id,
                    corstr = 'ar1'),
    error = function(e) NULL)
  if (is.null(fit)) return(na_row(spec))
  sm <- summary(fit)$coefficients
  tgt <- intersect(p10_target(spec), rownames(sm))
  if (length(tgt) == 0) return(na_row(spec))
  tgt <- tgt[1]
  alpha <- as.numeric(fit$geese$alpha)
  X <- model.matrix(form, data = dl)
  y <- dl$Sx
  vc <- tryCatch(md_vcov(X, y, coef(fit), dl$id, alpha),
                 error = function(e) NULL)
  md_se <- if (is.null(vc)) NA_real_ else
    sqrt(vc[match(tgt, names(coef(fit))), match(tgt, names(coef(fit)))])
  tibble(spec = spec, estimate = sm[tgt, 'Estimate'], alpha = alpha,
         lz_se = sm[tgt, 'Std.err'], md_se = md_se,
         df = length(unique(dl$id)) - length(coef(fit)),
         converged = TRUE)
}

one_rep <- function(i) {
  dl <- p10_gen_long(model_param, design_set, cell$t1half)
  if (is.null(dl)) return(map_dfr(c('E1', 'E2', 'E3'), na_row) |>
                            mutate(rep = i, .before = 1))
  dl <- dl[order(dl$ptID, dl$t), ]
  dl$id <- as.integer(factor(dl$ptID))
  dl <- as.data.frame(dl)
  map_dfr(c('E1', 'E2', 'E3'), function(sp) fit_gee(dl, sp)) |>
    mutate(rep = i, .before = 1)
}

cat(sprintf('GEE run: %s, N=%d, c_bm=0, c.cfct=%.1f, %d reps\n',
            cell$design, cell$N, P10_CFCT, N_REPS))

plan(multicore, workers = max(1, parallelly::availableCores() - 1))
t0 <- Sys.time()
res <- future_map_dfr(seq_len(N_REPS), one_rep,
                      .options = furrr_options(seed = SEED))
elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
cat(sprintf('Completed in %.0f s.\n\n', elapsed))

res <- res |>
  mutate(lz_p_z = 2 * pnorm(-abs(estimate / lz_se)),
         lz_p_t = 2 * pt(-abs(estimate / lz_se), df),
         md_p_z = 2 * pnorm(-abs(estimate / md_se)),
         md_p_t = 2 * pt(-abs(estimate / md_se), df))

summ <- res |>
  filter(converged, !is.na(md_se)) |>
  group_by(spec) |>
  summarise(
    n           = n(),
    emp_se      = sd(estimate),
    lz_se_ratio = mean(lz_se) / sd(estimate),
    md_se_ratio = mean(md_se) / sd(estimate),
    rej_lz_z    = mean(lz_p_z < 0.05),
    rej_lz_t    = mean(lz_p_t < 0.05),
    rej_md_z    = mean(md_p_z < 0.05),
    rej_md_t    = mean(md_p_t < 0.05),
    alpha_med   = median(alpha),
    .groups = 'drop')

fails <- res |> group_by(spec) |>
  summarise(n = n(), fail = sum(!converged | is.na(md_se)),
            .groups = 'drop')

print(as.data.frame(summ), digits = 4)
print(as.data.frame(fails))

out_dir <- p10_out_dir(repo_root)
saveRDS(list(results = res, summary = summ, failures = fails,
             cell = cell, n_reps = N_REPS, seed = SEED,
             dgp = p10_dgp_meta(repo_root), elapsed_secs = elapsed),
        file.path(out_dir, 'diag-null-type1-gee.rds'))
cat(sprintf('Wrote %s\n', file.path(out_dir, 'diag-null-type1-gee.rds')))
