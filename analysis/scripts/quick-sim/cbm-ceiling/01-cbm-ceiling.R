## Exact c_bm ceilings under crossed correlation structures.
##
## c* = (v' M^{-1} v)^{-1/2}, from the Schur complement on the biomarker
## row (Appendix B.7). Crossed factors: within-factor block (CS, AR1),
## cross-factor off-diagonal (constant, decaying), off-drug coupling
## (graded, step), carryover half-life, design path.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/cbm-ceiling/01-cbm-ceiling.R
## Outputs: analysis/data/quick-sim/cbm-ceiling/{design,path,rho}.csv
## Whitepaper: docs/34-cbm-feasible-range.md

suppressMessages(pkgload::load_all('.', quiet = TRUE))
library(data.table)
out_dir <- 'analysis/data/quick-sim/cbm-ceiling'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
data(extracted_bp)
data(extracted_rp)

designs <- list(
  CO = buildtrialdesign(
    name_longform = 'CO', name_shortform = 'CO',
    timepoints = cumulative(rep(2.5, 8)),
    timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
    expectancies = rep(.5, 8),
    ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
                  pathB = c(0, 0, 0, 0, 1, 1, 1, 1))
  ),
  Hybrid = buildtrialdesign(
    name_longform = 'Hybrid', name_shortform = 'Hybrid',
    timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
    timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4',
                    'COd', 'COp'),
    expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
    ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
                  pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
                  pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
                  pathD = c(1, 1, 1, 0, 0, 0, 0, 1))
  ),
  `OL+BDC` = buildtrialdesign(
    name_longform = 'OLBDC', name_shortform = 'OLBDC',
    timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
    timeptnames = c('OL1', 'OL2', 'OL3', 'OL4', 'BD1', 'BD2',
                    'BD3', 'BD4'),
    expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
    ondrug = list(pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
                  pathB = c(1, 1, 1, 1, 1, 0, 0, 0))
  )
)

#' Coupling pattern v on the BR positions, with c_bm factored out.
coupling_pattern <- function(tp, t_half, coupling) {
  on <- tp$tod > 0
  phi <- if (t_half > 0) (1 / 2)^(tp$tsd / t_half) else rep(0, nrow(tp))
  phi[on | tp$tsd <= 0] <- 0
  off_val <- if (coupling == 'graded') phi else as.numeric(phi > 0)
  ifelse(on, 1, off_val)
}

#' Correlation matrix R over (bm, BL, tv_1..n, pb_1..n, br_1..n).
build_R <- function(tp, c_bm, t_half, rho, c1, cx,
                    within, cross, coupling) {
  n <- nrow(tp)
  w <- cumulative(tp$t_wk)
  gap <- abs(outer(w, w, '-'))
  A <- if (within == 'ar1') rho^gap else {
    m <- matrix(rho, n, n)
    diag(m) <- 1
    m
  }
  C <- if (cross == 'decay') cx * rho^gap else matrix(cx, n, n)
  diag(C) <- c1
  M <- rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
  v <- c(rep(0, 2 * n), coupling_pattern(tp, t_half, coupling))
  R <- diag(2 + 3 * n)
  R[3:(2 + 3 * n), 3:(2 + 3 * n)] <- M
  R[1, 3:(2 + 3 * n)] <- R[3:(2 + 3 * n), 1] <- c_bm * v
  list(R = R, M = M, v = v)
}

ceiling_exact <- function(M, v) {
  min_eig_M <- min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)
  if (min_eig_M <= 0) return(c(c_star = NA_real_, min_eig_M = min_eig_M))
  c(c_star = 1 / sqrt(drop(t(v) %*% solve(M, v))), min_eig_M = min_eig_M)
}

ceiling_bisect <- function(tp, t_half, rho, c1, cx, within, cross,
                           coupling) {
  min_eig <- function(cb) {
    R <- build_R(tp, cb, t_half, rho, c1, cx, within, cross,
                 coupling)$R
    min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  }
  if (min_eig(1) > 0) return(Inf)
  uniroot(min_eig, c(0, 1), tol = 1e-10)$root
}

base <- list(rho = 0.7, c1 = 0.2, cx = 0.1)

## ---------------------------------------------------------------
## 1. Validate the builder against buildSigma() (covar configuration)
## ---------------------------------------------------------------
cat('== Validation against buildSigma() ==\n')
max_err <- 0
for (dn in names(designs)) {
  for (tp in designs[[dn]]$trialpaths) {
    for (th in c(0, 0.5, 1)) {
      mp <- data.table(N = 10, c.bm = 0.45, carryover_t1half = th,
        c.tv = base$rho, c.pb = base$rho, c.br = base$rho,
        c.cf1t = base$c1, c.cfct = base$cx)
      pkg <- cov2cor(buildSigma(mp, extracted_rp, extracted_bp, tp,
        makePositiveDefinite = FALSE, dgp_architecture = 'mvn')$sigma)
      own <- build_R(tp, 0.45, th, base$rho, base$c1, base$cx,
                     'ar1', 'decay', 'graded')$R
      max_err <- max(max_err, max(abs(pkg - own)))
    }
  }
}
cat(sprintf('max |R_builder - R_package| over all paths and t1/2: %.2e\n',
            max_err))

## ---------------------------------------------------------------
## 2. Crossed ceiling table at production parameters
## ---------------------------------------------------------------
grid <- CJ(within = c('ar1', 'cs'), cross = c('decay', 'const'),
           coupling = c('graded', 'step'), t_half = c(0, 0.5, 1),
           design = names(designs))
path_rows <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i]
  paths <- designs[[g$design]]$trialpaths
  for (k in seq_along(paths)) {
    b <- build_R(paths[[k]], 1, g$t_half, base$rho, base$c1, base$cx,
                 g$within, g$cross, g$coupling)
    ce <- ceiling_exact(b$M, b$v)
    path_rows[[length(path_rows) + 1]] <- data.table(
      g, path = k, c_star = ce[['c_star']], min_eig_M = ce[['min_eig_M']])
  }
}
paths_dt <- rbindlist(path_rows)
design_dt <- paths_dt[, .(c_star = min(c_star), min_eig_M = min(min_eig_M),
                          binding_path = path[which.min(c_star)]),
                      by = .(within, cross, coupling, t_half, design)]

## ---------------------------------------------------------------
## 3. Bisection cross-check of the exact formula
## ---------------------------------------------------------------
chk <- paths_dt[!is.na(c_star)]
chk[, c_bisect := mapply(function(d, k, th, wi, cr, co)
  ceiling_bisect(designs[[d]]$trialpaths[[k]], th, base$rho, base$c1,
                 base$cx, wi, cr, co),
  design, path, t_half, within, cross, coupling)]
cat(sprintf('\n== Bisection check ==\nmax |c*_formula - c*_bisection|: %.2e over %d path configurations\n',
            max(abs(chk$c_star - chk$c_bisect)), nrow(chk)))

cat('\n== Design-level ceiling (min over paths), rho = 0.7, c1 = 0.2, cx = 0.1 ==\n')
wide <- dcast(design_dt, within + cross + coupling + design ~ t_half,
              value.var = 'c_star')
setnames(wide, c('0', '0.5', '1'), c('t0', 't0.5', 't1'))
setorder(wide, within, cross, coupling, design)
print(wide, digits = 3)

cat('\n== Smallest eigenvalue of M (response block alone) ==\n')
print(unique(design_dt[, .(min_eig_M = round(min(min_eig_M), 4)),
                       by = .(within, cross, design)]))

cat('\n== Binding path (covar configuration) ==\n')
print(design_dt[within == 'ar1' & cross == 'decay' & coupling == 'graded',
                .(design, t_half, c_star = round(c_star, 4), binding_path)])

## ---------------------------------------------------------------
## 4. Sensitivity to rho (design-level ceilings)
## ---------------------------------------------------------------
sens <- list()
for (rho in c(0.3, 0.5, 0.7, 0.9)) {
  for (cfg in list(c('ar1', 'decay', 'graded'), c('cs', 'const', 'step'),
                   c('cs', 'const', 'graded'))) {
    for (th in c(0, 1)) {
      for (dn in names(designs)) {
        cs <- sapply(designs[[dn]]$trialpaths, function(tp) {
          b <- build_R(tp, 1, th, rho, base$c1, base$cx,
                       cfg[1], cfg[2], cfg[3])
          ceiling_exact(b$M, b$v)[['c_star']]
        })
        sens[[length(sens) + 1]] <- data.table(rho = rho,
          config = paste(cfg, collapse = '/'), t_half = th,
          design = dn, c_star = min(cs))
      }
    }
  }
}
sens_dt <- dcast(rbindlist(sens), config + t_half + design ~ rho,
                 value.var = 'c_star')
cat('\n== Ceiling by rho (columns) ==\n')
print(sens_dt, digits = 3)

## ---------------------------------------------------------------
## 5. Silent PD repair just below and above the covar ceiling
## ---------------------------------------------------------------
cat('\n== makePositiveDefinite repair around c* (covar, t1/2 = 0) ==\n')
for (dn in names(designs)) {
  cstar <- design_dt[within == 'ar1' & cross == 'decay' &
                     coupling == 'graded' & t_half == 0 &
                     design == dn, c_star]
  k <- design_dt[within == 'ar1' & cross == 'decay' &
                 coupling == 'graded' & t_half == 0 &
                 design == dn, binding_path]
  tp <- designs[[dn]]$trialpaths[[k]]
  for (cb in c(0.99, 1.01) * cstar) {
    mp <- data.table(N = 10, c.bm = cb, carryover_t1half = 0,
      c.tv = base$rho, c.pb = base$rho, c.br = base$rho,
      c.cf1t = base$c1, c.cfct = base$cx)
    raw <- buildSigma(mp, extracted_rp, extracted_bp, tp,
      makePositiveDefinite = FALSE, dgp_architecture = 'mvn')
    adj <- buildSigma(mp, extracted_rp, extracted_bp, tp,
      makePositiveDefinite = TRUE, dgp_architecture = 'mvn')
    br <- grep('\\.br$', raw$labels)
    r_raw <- cov2cor(raw$sigma)[1, br]
    r_adj <- cov2cor(adj$sigma)[1, br]
    cat(sprintf('%-7s path %d  c_bm = %.4f  max|r_adj - r_spec| = %.4f  on-drug r_adj = %s\n',
      dn, k, cb, max(abs(r_adj - r_raw)),
      paste(sprintf('%.3f', r_adj[tp$tod > 0]), collapse = ' ')))
  }
}

fwrite(design_dt, file.path(out_dir, 'design.csv'))
fwrite(paths_dt, file.path(out_dir, 'path.csv'))
fwrite(rbindlist(sens), file.path(out_dir, 'rho.csv'))
