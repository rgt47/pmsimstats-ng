## analysis/scripts/dgp-combined/03-analytic-projection.R
##
## Closed-form quantities for paper 11 (referee report 2026-10-28, M4,
## M5 and suggestion 3), computed from the covariance matrices the
## extended DGP builds. No simulation.
##
## For each design, carryover half-life and channel-weight pair:
##   1. The slope profile of the outcome on the biomarker, visit by
##      visit, including the baseline visit and the off-drug visits:
##        s_t = -(c_bm_a * D_t + c_bm_b * g_t) * sigma_BR / sigma_bm
##      in outcome points per biomarker unit. The mean channel acts at
##      on-drug visits only (D_t); the covariance channel follows the
##      correlation profile g_t (1 on drug, the carryover decay off
##      drug after exposure, 0 before first exposure).
##   2. The conditional covariance V of the outcome vector given the
##      biomarker. Only the covariance channel changes it.
##   3. The population value of the E2 interaction coefficient: the
##      generalized-least-squares projection of s onto the biomarker
##      regressors (bm, bm:Dbc) with Dbc at the true half-life, pooled
##      over paths with weights n_p V_p^{-1}.
##   4. The model-free standard error of that GLS estimator,
##      sqrt([(sigma_bm^2 sum_p n_p X_p' V_p^{-1} X_p)^{-1}]_22), and
##      the noncentrality Lambda = theta / SE.
##   5. The on-drug conditional variance of BR and of the outcome.
##
## Output: analysis/data/dgp-combined/03-analytic.rds (a list with
## `cells` and `visits`).

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))
e <- new.env()
load(file.path(repo_root, 'data/extracted_rp.rda'), envir = e)
load(file.path(repo_root, 'data/extracted_bp.rda'), envir = e)
rp <- e$extracted_rp
bp <- e$extracted_bp

N_TOTAL <- 70L
sd_bm <- bp$sd[bp$cat == 'bm']
sd_br <- rp$sd[rp$cat == 'br']

allocate <- function(n, k) rep(n %/% k, k) + c(rep(1L, n %% k),
                                               rep(0L, k - n %% k))

path_moments <- function(td, t1half, c_bm_b) {
  mp <- list(N = 10, carryover_t1half = t1half,
             carryover_form = 'exponential', weibull_shape = 1,
             dgp_architecture = 'combined',
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
             c.cf1t = 0.2, c.cfct = 0.2, c.bm_a = 0, c.bm_b = c_bm_b)
  sg <- build_sigma_matrix(mp, rp, bp, td, dgp_architecture = 'combined',
                           stop_on_nonpd = TRUE)
  S <- sg$sigma
  lab <- sg$labels
  tp <- td$timepoint_name
  n_y <- length(tp) + 1L
  C <- matrix(0, n_y, length(lab))
  C[, match('BL', lab)] <- 1
  for (k in seq_along(tp)) {
    C[k + 1L, match(paste0(tp[k], c('.tv', '.pb', '.br')), lab)] <- -1
  }
  ib <- match('bm', lab)
  syy <- C %*% S %*% t(C)
  syb <- drop(C %*% S[, ib])
  V <- syy - tcrossprod(syb) / S[ib, ib]
  D <- c(0, as.numeric(td$tod > 0))
  tsd <- c(0, td$tsd)
  dbc <- exposure_weight(D, tsd, t1half)
  ibr <- match(paste0(tp, '.br'), lab)
  list(
    slope_cov = syb / S[ib, ib],
    D = D, dbc = dbc, V = V,
    min_eig = min(eigen(S[sg$standard_deviations > 0,
                          sg$standard_deviations > 0],
                        only.values = TRUE)$values),
    var_br_marg = diag(S)[ibr],
    var_br_cond = diag(S)[ibr] - S[ibr, ib]^2 / S[ib, ib],
    var_y_cond = diag(V)[-1]
  )
}

grid <- tidyr::expand_grid(
  design = c('CO', 'Hybrid', 'OLBDC'),
  t1half = c(0, 0.2, 0.5, 1),
  c_bm_a = c(0, 0.225, 0.45),
  c_bm_b = c(0, 0.225, 0.45)
)

cells <- pmap_dfr(grid, function(design, t1half, c_bm_a, c_bm_b) {
  ds <- build_design_set(design)
  n_p <- allocate(N_TOTAL, length(ds))
  pm <- map(ds, path_moments, t1half = t1half, c_bm_b = c_bm_b)
  info <- matrix(0, 2, 2)
  score <- c(0, 0)
  for (p in seq_along(pm)) {
    X <- cbind(1, pm[[p]]$dbc)
    Vi <- solve(pm[[p]]$V)
    s <- pm[[p]]$slope_cov - c_bm_a * sd_br / sd_bm * pm[[p]]$D
    info <- info + n_p[p] * t(X) %*% Vi %*% X
    score <- score + n_p[p] * drop(t(X) %*% Vi %*% s)
  }
  th <- solve(info, score)
  se <- sqrt(solve(sd_bm^2 * info)[2, 2])
  tibble(design, t1half, c_bm_a, c_bm_b,
         theta = th[2], theta_bm = th[1], se_gls = se,
         Lambda = abs(th[2]) / se,
         nominal = -(c_bm_a + c_bm_b) * sd_br / sd_bm,
         min_eig = min(map_dbl(pm, 'min_eig')))
})

## Per-unit coefficients: theta = c_bm_a * kappa_A + c_bm_b * kappa_B.
kappas <- cells |>
  group_by(design, t1half) |>
  summarise(
    kappa_A = theta[c_bm_a == 0.45 & c_bm_b == 0] / 0.45,
    kappa_B = theta[c_bm_a == 0 & c_bm_b == 0.45] / 0.45,
    max_dev_linear = max(abs(theta - (c_bm_a * kappa_A +
                                      c_bm_b * kappa_B))),
    .groups = 'drop'
  ) |>
  mutate(nominal_per_unit = -sd_br / sd_bm,
         overshoot_A = kappa_A / nominal_per_unit - 1,
         overshoot_B = kappa_B / nominal_per_unit - 1)

## Visit-level moments for the Hybrid design, path 1, at t1/2 = 1
visits <- map_dfr(c(0, 0.45), function(b) {
  td <- build_design_set('Hybrid')[[1]]
  pm <- path_moments(td, 1, b)
  tibble(c_bm_b = b, visit = td$timepoint_name, on_drug = td$tod > 0,
         var_br_marg = pm$var_br_marg, var_br_cond = pm$var_br_cond,
         var_y_cond = pm$var_y_cond,
         slope_cov = pm$slope_cov[-1])
})

out <- list(cells = cells, kappas = kappas, visits = visits,
            sd_bm = sd_bm, sd_br = sd_br,
            created = format(Sys.time(), '%Y-%m-%d %H:%M %Z'))
saveRDS(out, file.path(repo_root,
  'analysis/data/dgp-combined/03-analytic.rds'))

print(as.data.frame(kappas), digits = 4)
print(as.data.frame(cells |> filter(t1half %in% c(0, 1),
                                    c_bm_a + c_bm_b > 0.44,
                                    c_bm_a + c_bm_b < 0.46)),
      digits = 4)
print(as.data.frame(visits), digits = 4)
cat('Smallest eigenvalue over all cells:', min(cells$min_eig), '\n')
