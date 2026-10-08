## How often would each configuration's correlation matrix be silently
## repaired? Counts, for the strawman (Hendrickson's 58b32a9 code) and for
## configurations A-F of 04-power-simulation.R, the matrices that fail the
## test the published code uses before repairing,
## corpcor::is.positive.definite() on the covariance matrix.
##
## One matrix per design path x c.bm x t1/2: N does not enter the matrix,
## so the published grid's 162 builds (N = 35 and 70) are 81 distinct
## matrices. Two grids are counted:
##   published  the Figure 4 grid: OL, OL+BDC, CO, N-of-1 (Hybrid);
##              c.bm {0, .3, .6}; t1/2 {0, .1, .2}
##   extended   the 04-power-simulation.R grid: OL+BDC, CO, Hybrid;
##              c.bm {0, .25, .3, .6}; t1/2 {0, .1, .2, .5, 1}
## The strawman column comes from 01-published-matrices.csv (58b32a9
## generateData() itself); configuration A re-derives the same matrices
## with the 04 generator, so the two must agree.
##
## Configurations D and E carry no biomarker coupling in the matrix, so
## they cannot fail on account of c.bm; they are counted anyway.
##
## Ceilings. A failure has one of two sources: the response block M
## (TV, PB, BR visits) is itself not positive definite, or M is fine and
## the biomarker coupling exceeds what it can support. For each
## configuration with coupling (A, B, F, C), design and t1/2 the script
## reports M's smallest eigenvalue and the c.bm ceiling
##   c* = min over paths (u' M^-1 u)^(-1/2),
## u the coupling pattern at c.bm = 1 (BL is uncorrelated with the rest,
## so it drops out). c.bm^2 u' M^-1 u is the squared multiple correlation
## of the biomarker on the response trajectory, so it must stay below 1.
## The ceiling is checked against the counts: a path fails exactly when
## c.bm exceeds its own ceiling (corpcor's tolerance can move a cell
## lying within ~1e-3 of the ceiling; any such cell is listed).
##
## Run from the repository root (after 01-published-matrices.R):
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/11-pd-repair-count.R
suppressMessages({library(data.table); library(corpcor)})
here <- 'analysis/scripts/quick-sim/hendrickson-problems'
out_dir <- 'analysis/data/quick-sim/hendrickson-problems'

## Take the DGP definitions from 04 verbatim, without running its
## simulation: evaluate only the top-level assignments named here.
keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build')
dgp <- new.env()
## Lines such as 'RHO <- 0.8; C1 <- 0.2' parse as separate expressions.
for (e in parse(file.path(here, '04-power-simulation.R'))) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
missing_defs <- setdiff(keep, ls(dgp))
if (length(missing_defs)) stop('not found in 04: ', paste(missing_defs, collapse = ', '))
designs <- dgp$designs
designs$OL <- list(t_wk = rep(2.5, 8), e = rep(1, 8), paths = list(rep(1, 8)))

count_arm <- function(arm, dz, cbms, halves) {
  rows <- list()
  for (dn in dz) for (cbm in cbms) for (th in halves) {
    d <- designs[[dn]]
    for (g in seq_along(d$paths)) {
      pi <- dgp$path_info(d$t_wk, d$e, d$paths[[g]])
      b <- dgp$build(pi, arm, cbm, th)
      rows[[length(rows) + 1]] <- data.table(
        arm = arm, design = dn, path = g, c.bm = cbm, t_half = th,
        repaired = !is.positive.definite(b$Sigma),
        min_eig = min(eigen(b$R, TRUE, TRUE)$values))
    }
  }
  rbindlist(rows)
}

arms <- c('A', 'B', 'F', 'C', 'D', 'E')
arm_lab <- c(A = 'A published construct', B = 'B AR(1)',
             F = 'F AR(1) + graded, covar', C = 'C AR(1) + graded, separable',
             D = 'D mean moderation', E = 'E all three')

pub_d <- c('OL', 'OL+BDC', 'CO', 'Hybrid')
grid_pub <- rbindlist(lapply(arms, count_arm, pub_d, c(0, .3, .6), c(0, .1, .2)))
grid_ext <- rbindlist(lapply(arms, count_arm, c('OL+BDC', 'CO', 'Hybrid'),
                             c(0, .25, .3, .6), c(0, .1, .2, .5, 1)))

## Strawman reference: the matrices 58b32a9 generateData() actually built.
sm <- fread(file.path(out_dir, 'published-matrices.csv'))
sm <- unique(sm[, .(design = fifelse(design == 'N-of-1', 'Hybrid', design),
                    path, c.bm, t_half, repaired = !pd)])
stopifnot(nrow(sm) == 81)
chk <- merge(sm, grid_pub[arm == 'A'], by = c('design', 'path', 'c.bm', 't_half'),
             suffixes = c('_58b32a9', '_A'))
stopifnot(nrow(chk) == 81)
n_disagree <- sum(chk$repaired_58b32a9 != chk$repaired_A)
cat(sprintf('Check: configuration A vs 58b32a9 generateData(), 81 matrices: %d disagree\n',
            n_disagree))

tab <- function(d, by) {
  d[, .(matrices = .N, repaired = sum(repaired),
        pct = round(100 * mean(repaired), 1)), by = by]
}
fmt <- function(d) {
  d[, cell := sprintf('%d/%d (%.0f%%)', repaired, matrices, pct)]
  d
}

strawman_rows <- sm[, .(arm = 'S', design, c.bm, repaired)]
pub_all <- rbind(strawman_rows, grid_pub[, .(arm, design, c.bm, repaired)])
lab <- c(S = 'Strawman (58b32a9 code)', arm_lab)
pub_all[, config := factor(lab[arm], levels = lab[c('S', arms)])]

cat('\nPublished grid: matrices that would be silently repaired',
    '(c.bm > 0 only; c.bm = 0 never fails)\n')
by_design <- fmt(tab(pub_all[c.bm > 0], c('config', 'design')))
print(dcast(by_design, config ~ design, value.var = 'cell'))
overall <- fmt(tab(pub_all, 'config'))
overall_cbm <- fmt(tab(pub_all[c.bm > 0], 'config'))
cat('\nPublished grid, all designs\n')
print(merge(overall[, .(config, `all c.bm` = cell)],
            overall_cbm[, .(config, `c.bm > 0` = cell)], by = 'config'))

cat('\nExtended grid (04 power simulation), all designs\n')
grid_ext[, config := factor(arm_lab[arm], levels = arm_lab[arms])]
ext <- fmt(tab(grid_ext, c('config', 'c.bm')))
print(dcast(ext, config ~ c.bm, value.var = 'cell'))

## ---------------- ceilings ----------------
ceil_rows <- list()
for (arm in c('A', 'B', 'F', 'C')) for (dn in pub_d) for (th in c(0, .1, .2, .5, 1)) {
  d <- designs[[dn]]
  for (g in seq_along(d$paths)) {
    b <- dgp$build(dgp$path_info(d$t_wk, d$e, d$paths[[g]]), arm, 1, th)
    M <- b$R[-(1:2), -(1:2)]
    u <- b$R[1, -(1:2)]
    me <- min(eigen(M, TRUE, TRUE)$values)
    ## A response block that is not positive definite fails at any c.bm.
    cs <- if (me > 0) 1 / sqrt(drop(crossprod(u, solve(M, u)))) else 0
    ceil_rows[[length(ceil_rows) + 1]] <- data.table(
      arm = arm, design = dn, path = g, t_half = th, min_eig_M = me, ceiling = cs)
  }
}
ceil_path <- rbindlist(ceil_rows)
ceil <- ceil_path[, .(min_eig_M = min(min_eig_M), ceiling = min(ceiling)),
                  by = .(arm, design, t_half)]

## Consistency with the counts: a path fails exactly when c.bm > ceiling.
counted <- rbind(grid_pub[arm %in% c('A', 'B', 'F', 'C') & c.bm > 0],
                 grid_ext[arm %in% c('A', 'B', 'F', 'C') & c.bm > 0], fill = TRUE)
counted <- unique(counted[, .(arm, design, path, c.bm, t_half, repaired)])
cc <- merge(counted, ceil_path, by = c('arm', 'design', 'path', 't_half'))
stopifnot(nrow(cc) == nrow(counted))
cc[, predicted := c.bm > ceiling]
off <- cc[predicted != repaired]
cat(sprintf('\nCheck: ceiling predicts the repair count for %d of %d matrices\n',
            nrow(cc) - nrow(off), nrow(cc)))
if (nrow(off)) print(off[, .(arm, design, path, c.bm, t_half, ceiling, repaired)])

cat('\nResponse block M: smallest eigenvalue over paths and t1/2 (no c.bm involved)\n')
print(dcast(ceil[, .(m = round(min(min_eig_M), 3)), by = .(config = arm_lab[arm], design)],
            config ~ design, value.var = 'm'))
cat('\nc.bm ceiling: largest c.bm with every path positive definite (columns: t1/2)\n')
print(dcast(ceil[, .(config = factor(arm_lab[arm], levels = arm_lab[arms]), design,
                     t_half, c = round(ceiling, 3))],
            config + design ~ t_half, value.var = 'c'))
fwrite(ceil[, .(config = arm_lab[arm], design, t_half, min_eig_M, ceiling)],
       file.path(out_dir, 'pd-ceilings.csv'))

fwrite(rbind(pub_all[, .(grid = 'published', config, design, c.bm, repaired)],
             grid_ext[, .(grid = 'extended', config, design, c.bm, repaired)]),
       file.path(out_dir, 'pd-repair-matrices.csv'))
fwrite(rbind(cbind(grid = 'published', by_design),
             cbind(grid = 'published-overall', design = 'all', overall),
             cbind(grid = 'extended', design = 'all', ext), fill = TRUE),
       file.path(out_dir, 'pd-repair-count.csv'))
cat('\nWrote', file.path(out_dir, 'pd-repair-count.csv'), 'and',
    file.path(out_dir, 'pd-ceilings.csv'), '\n')
