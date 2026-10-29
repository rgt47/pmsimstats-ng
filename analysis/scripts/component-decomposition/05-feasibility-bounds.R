## 05-feasibility-bounds.R  (paper 06, 2026-10-28)
##
## Feasibility of the contaminated-biomarker cells, by the Schur-
## complement argument of docs/65 and the compendium (eq. ce-ceiling),
## applied to the three trial designs paper 06 uses:
##
##   simple    eight-visit single-path OL+BDC (Study 1 at N <= 70,
##             Study 2 at N = 70)
##   full      sixteen-visit two-path OL+BDC+crossover (Study 1 at
##             N >= 200, Study 2 at N = 300, Study 3)
##   balanced  OL / covert on-drug / concealed placebo (Study 4)
##
## With the biomarker row removed, the correlation matrix R_0 does not
## depend on any biomarker correlation. The biomarker row is
##     r = c_bm * g + c_pb * p,
## with g the BR pattern (1 on drug, decaying off drug) and p the PB
## pattern (1 at every PB position). R is positive definite iff R_0 is
## and r' R_0^{-1} r < 1, a quadratic in c_pb at fixed c_bm. The design
## bound is the minimum over its paths.
##
## For every cell of the stored contamination run (study-contam-
## alt100-null100) and hybrid recovery run (study-b-recovery-pilot)
## the script also reports what the June 2026 repair
## (corpcor::make.positive.definite, tol = 1e-3, the package's old
## repair path) did to the matrix: the repaired on-drug biomarker-BR
## correlation, the repaired population value of the bm:Dbc contrast,
## and the largest change in any correlation. Those cells drew from
## the repaired matrix, so bias in them must be read against the
## repaired population value, not against -2.25. The repair is
## applied here directly with corpcor on a copy of the raw matrix; the
## package's own repair branch is never enabled.
##
## The identity prediction (manuscript Section 6.1, corollary) is
## reported two ways: the one-line approximation
##     -c_pb * sigma_PB * (eta_on - eta_off)
## with eta_on and eta_off the design's mean expectancy over on-drug
## and off-drug visits (binary drug state, pooled over paths), and the
## exact population on-minus-off contrast of the biomarker slope,
## computed from the raw covariance matrix as in docs/65 Section 3.3.
##
## Usage (from the repository root, host R):
##   Rscript analysis/scripts/component-decomposition/05-feasibility-bounds.R
## Output:
##   analysis/data/derived_data/component-decomposition/feasibility-bounds.rds
## Runs in seconds; no simulation.

suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table)
})

OUT_DIR <- file.path('analysis', 'data', 'derived_data',
                     'component-decomposition')
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------------ ##
## Designs (identical to the drivers)
## ------------------------------------------------------------------ ##

design_simple <- function() buildtrialdesign(
  name_longform = 'OL + BDC (simple)', name_shortform = 'OLBDC',
  timepoints = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20),
  timeptnames = c(paste0('OL', 1:5), paste0('BDC', 1:3)),
  expectancies = c(rep(1, 5), rep(0.5, 3)),
  ondrug = list(pathA = c(rep(1, 5), 1, 0, 0)))

design_full <- function() buildtrialdesign(
  name_longform = 'OL + BDC + crossover (full)', name_shortform = 'OLBDCxo',
  timepoints = seq(1, 20, by = 1.25),
  timeptnames = c(paste0('OL', 1:8), paste0('BDC', 1:4),
                  paste0('XO1_', 1:2), paste0('XO2_', 1:2)),
  expectancies = c(rep(1, 8), rep(0.5, 4), rep(0.5, 4)),
  ondrug = list(AB = c(rep(1, 8), rep(0, 4), rep(1, 2), rep(0, 2)),
                BA = c(rep(1, 8), rep(0, 4), rep(0, 2), rep(1, 2))))

design_balanced <- function() buildtrialdesign(
  name_longform = 'balanced-placebo (OL/covert/concealed placebo)',
  name_shortform = 'BAL',
  timepoints = seq(1.5, 22.5, by = 1.5),
  timeptnames = c(paste0('OL', 1:3), paste0('COV', 1:6), paste0('OLP', 1:6)),
  expectancies = c(rep(1, 3), rep(0, 6), rep(1, 6)),
  ondrug = list(p = c(rep(1, 3), rep(1, 6), rep(0, 6))))

designs <- list(simple = design_simple(), full = design_full(),
                balanced = design_balanced())

make_resp_params <- function(m_PB = 6, m_TV = 0, m_BR = 10.98604)
  data.table(cat = c('tv', 'pb', 'br'), max = c(m_TV, m_PB, m_BR),
             disp = c(5, 5, 5), rate = c(0.42, 0.35, 0.42), sd = c(5, 2, 5))
make_bl_params <- function()
  data.table(cat = c('bm', 'BL'), m = c(0, 70), sd = c(1, 10))
make_model_params <- function(c_bm = 0.45, c_bm_pb = 0)
  data.table(N = 50, c.bm = c_bm, c.bm.pb = c_bm_pb, carryover_t1half = 1,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.1, c.cfct = 0.05)

SIGMA_PB <- 2
SIGMA_BR <- 5

## Raw covariance (never repaired): stopOnNonPD = FALSE and
## makePositiveDefinite = FALSE return the matrix as built. Without
## options(pmsimstats.allow_repair) the package repairs nothing.
raw_sigma <- function(path, c_bm, c_bm_pb, rp = make_resp_params()) {
  s <- buildSigma(make_model_params(c_bm, c_bm_pb), rp, make_bl_params(),
                  path, makePositiveDefinite = FALSE, stopOnNonPD = FALSE)
  s$sigma
}

## Positive-variance sub-block (the balanced design has exactly zero
## PB variance at its covert visits, which is degenerate, not non-PD).
pos_block <- function(sig) {
  pos <- sqrt(diag(sig)) > 0
  sig[pos, pos, drop = FALSE]
}

min_eig <- function(sig)
  min(eigen(pos_block(sig), symmetric = TRUE, only.values = TRUE)$values)

## ------------------------------------------------------------------ ##
## Population bm:Dbc contrast from a covariance matrix (docs/65, 3.3)
## ------------------------------------------------------------------ ##
## Sx_t = BL - TV_t - PB_t - BR_t, so the slope of Sx_t on bm is
##   (Cov(bm,BL) - Cov(bm,TV_t) - Cov(bm,PB_t) - Cov(bm,BR_t)) / Var(bm)
## and the on-minus-off average of that slope is the population value
## the one-component bm:Dbc coefficient targets (binary drug state).

slope_contrast <- function(sig, path) {
  lbl <- rownames(sig)
  v_bm <- sig['bm', 'bm']
  tp <- path$timeptnames
  on <- path$tod > 0
  comp_cov <- function(cc) vapply(tp, function(nm) {
    key <- paste(nm, cc, sep = '.')
    if (key %in% lbl) sig['bm', key] else 0
  }, numeric(1))
  cv_br <- comp_cov('br'); cv_pb <- comp_cov('pb'); cv_tv <- comp_cov('tv')
  v_br_on <- vapply(tp[on], function(nm)
    sig[paste(nm, 'br', sep = '.'), paste(nm, 'br', sep = '.')], numeric(1))
  ## Slope of Sx_t on bm, component by component (sign: reductions).
  ## The one-component analysis carries the drug state as D_bc, which
  ## absorbs the off-drug decay of the BR correlation, so its target
  ## under the DGP is the on-drug BR slope plus the on-minus-off
  ## displacement from any component the biomarker is correlated
  ## with (Section 6.1 of the manuscript; docs/65 Section 3.3).
  disp <- function(cv) -(mean(cv[on]) - mean(cv[!on])) / v_bm
  list(beta_br_on = -mean(cv_br[on]) / v_bm,
       pb_disp = disp(cv_pb), tv_disp = disp(cv_tv),
       cor_bm_br_on = mean(cv_br[on] / sqrt(v_bm * v_br_on)))
}

## ------------------------------------------------------------------ ##
## Schur-complement bound on c_pb given c_bm, per design
## ------------------------------------------------------------------ ##

path_pieces <- function(path) {
  R_full <- cov2cor(pos_block(raw_sigma(path, 0, 0)))
  i_bm <- which(rownames(R_full) == 'bm')
  R0 <- R_full[-i_bm, -i_bm]
  stopifnot(min(eigen(R0, symmetric = TRUE, only.values = TRUE)$values) > 0)
  pattern <- function(c_bm, c_pb)
    cov2cor(pos_block(raw_sigma(path, c_bm, c_pb)))[i_bm, -i_bm]
  X <- cbind(g = pattern(1, 0), p = pattern(0, 1))
  list(R0 = R0, X = X, A = t(X) %*% solve(R0, X))
}

bound_given <- function(A, y, fixed) {
  other <- setdiff(colnames(A), y)
  a <- A[y, y]
  bq <- 2 * A[y, other] * fixed
  cq <- A[other, other] * fixed^2 - 1
  disc <- bq^2 - 4 * a * cq
  if (disc < 0) return(NA_real_)
  (-bq + sqrt(disc)) / (2 * a)
}

bounds <- rbindlist(lapply(names(designs), function(dn) {
  td <- designs[[dn]]
  pieces <- lapply(td$trialpaths, path_pieces)
  per_path <- function(f) min(vapply(pieces, f, numeric(1)))
  data.table(
    design = dn,
    c_bm_alone = per_path(function(p) 1 / sqrt(p$A['g', 'g'])),
    c_pb_alone = per_path(function(p) 1 / sqrt(p$A['p', 'p'])),
    c_pb_given_cbm_0.30 = per_path(function(p) bound_given(p$A, 'p', 0.30)),
    c_pb_given_cbm_0.45 = per_path(function(p) bound_given(p$A, 'p', 0.45)),
    c_bm_given_cpb_0.30 = per_path(function(p) bound_given(p$A, 'g', 0.30)),
    c_bm_given_cpb_0.45 = per_path(function(p) bound_given(p$A, 'g', 0.45)))
}))

## Eigenvalue verification of each bound at c_bm = 0.45: positive just
## below, negative just above.
verify <- rbindlist(lapply(names(designs), function(dn) {
  cstar <- bounds[design == dn, c_pb_given_cbm_0.45]
  me <- function(cpb) min(vapply(designs[[dn]]$trialpaths, function(p)
    min_eig(raw_sigma(p, 0.45, cpb)), numeric(1)))
  data.table(design = dn, bound = cstar,
             min_eig_below = me(0.999 * cstar),
             min_eig_above = me(1.001 * cstar))
}))
stopifnot(all(verify$min_eig_below > 0 & verify$min_eig_above < 0))

## ------------------------------------------------------------------ ##
## Expectancy loading and the one-line identity prediction
## ------------------------------------------------------------------ ##

loading <- rbindlist(lapply(names(designs), function(dn) {
  td <- designs[[dn]]
  eta_on <- unlist(lapply(td$trialpaths, function(p) p$e[p$tod > 0]))
  eta_off <- unlist(lapply(td$trialpaths, function(p) p$e[!(p$tod > 0)]))
  data.table(design = dn, eta_on = mean(eta_on), eta_off = mean(eta_off),
             w_pb = mean(eta_on) - mean(eta_off))
}))

## ------------------------------------------------------------------ ##
## Per-cell table: raw eigenvalue, repaired matrix, population values
## ------------------------------------------------------------------ ##

grid <- CJ(design = names(designs),
           c_bm_pb = c(0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.45),
           sorted = FALSE)

cells <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  dn <- grid$design[i]; cpb <- grid$c_bm_pb[i]
  td <- designs[[dn]]
  rbindlist(lapply(seq_along(td$trialpaths), function(g) {
    path <- td$trialpaths[[g]]
    sig <- raw_sigma(path, 0.45, cpb)
    sub <- pos_block(sig)
    ev <- min(eigen(sub, symmetric = TRUE, only.values = TRUE)$values)
    raw <- slope_contrast(sig, path)
    ## What the old repair path did to a non-PD matrix.
    rep_sub <- if (ev > 0) sub else
      corpcor::make.positive.definite(sub, tol = 1e-3)
    rep_sig <- sig
    rep_sig[rownames(sub), colnames(sub)] <- rep_sub
    rep <- slope_contrast(rep_sig, path)
    data.table(
      design = dn, path = names(td$trialpaths)[g], c_bm = 0.45,
      c_bm_pb = cpb, min_eig = ev, pd = ev > 0,
      cor_bm_br_on_raw = raw$cor_bm_br_on,
      cor_bm_br_on_repaired = rep$cor_bm_br_on,
      beta_br_raw = raw$beta_br_on,
      beta_br_repaired = rep$beta_br_on,
      pb_disp_raw = raw$pb_disp,
      pb_disp_repaired = rep$pb_disp,
      tv_disp_repaired = rep$tv_disp,
      target_raw = raw$beta_br_on + raw$pb_disp + raw$tv_disp,
      target_repaired = rep$beta_br_on + rep$pb_disp + rep$tv_disp,
      max_abs_cor_change = if (ev > 0) 0 else
        max(abs(cov2cor(rep_sub) - cov2cor(sub))),
      w_pb = loading[design == dn, w_pb],
      identity_approx = -cpb * SIGMA_PB * loading[design == dn, w_pb])
  }))
}))

## Design-level view: minimum eigenvalue over paths, path-mean of the
## population quantities (the two paths of the full design differ in
## the third decimal).
cells_design <- cells[, .(
  min_eig = min(min_eig), pd = all(pd),
  cor_bm_br_on_raw = mean(cor_bm_br_on_raw),
  cor_bm_br_on_repaired = mean(cor_bm_br_on_repaired),
  beta_br_raw = mean(beta_br_raw),
  beta_br_repaired = mean(beta_br_repaired),
  pb_disp_raw = mean(pb_disp_raw),
  pb_disp_repaired = mean(pb_disp_repaired),
  tv_disp_repaired = mean(tv_disp_repaired),
  target_raw = mean(target_raw),
  target_repaired = mean(target_repaired),
  max_abs_cor_change = max(max_abs_cor_change),
  w_pb = w_pb[1], identity_approx = identity_approx[1]),
  by = .(design, c_bm, c_bm_pb)]

## ------------------------------------------------------------------ ##
## Study 1 grid: feasibility does not depend on m_PB or m_TV, which
## enter the means only; checked at the grid corners and at the
## reference cell on every path.
## ------------------------------------------------------------------ ##

study1 <- rbindlist(lapply(c('simple', 'full'), function(dn) {
  td <- designs[[dn]]
  rbindlist(lapply(list(c(0, -1), c(0, 0), c(10, 2), c(6, 0)), function(mm) {
    rp <- make_resp_params(m_PB = mm[1], m_TV = mm[2])
    data.table(design = dn, m_PB = mm[1], m_TV = mm[2],
               min_eig = min(vapply(td$trialpaths, function(p)
                 min_eig(raw_sigma(p, 0.45, 0, rp)), numeric(1))))
  }))
}))
stopifnot(all(study1$min_eig > 0))
stopifnot(all(study1[, diff(range(min_eig)) < 1e-12, by = design]$V1))

## ------------------------------------------------------------------ ##
## Report and save
## ------------------------------------------------------------------ ##

options(width = 160)
cat('Feasibility bounds (c_bm = 0.45, sigma_BR = 5, sigma_PB = 2, rho = 0.7,\n',
    'c_1 = 0.1, c_x = 0.05, t1/2 = 1 week), minimum over paths:\n')
print(bounds[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat('\nBound verification (eigenvalue just below / above):\n')
print(verify[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
cat('\nExpectancy loading w_PB = eta_on - eta_off:\n')
print(loading[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
cat('\nPer-cell table (design level):\n')
print(cells_design[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat('\nStudy 1 grid corners (min eigenvalue, independent of m_PB, m_TV):\n')
print(study1)

out <- list(bounds = bounds, verify = verify, loading = loading,
            cells = cells, cells_design = cells_design, study1 = study1,
            meta = list(date = Sys.time(), c_bm = 0.45, sigma_pb = SIGMA_PB,
                        sigma_br = SIGMA_BR, rho = 0.7, c1 = 0.1,
                        cx = 0.05, t1half = 1,
                        repair = 'corpcor::make.positive.definite(tol = 1e-3) on the positive-variance block'))
saveRDS(out, file.path(OUT_DIR, 'feasibility-bounds.rds'))
cat(sprintf('\nWritten: %s\n', file.path(OUT_DIR, 'feasibility-bounds.rds')))
