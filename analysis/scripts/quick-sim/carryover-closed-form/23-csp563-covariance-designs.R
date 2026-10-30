## Closed-form design comparison under the CAPS covariance estimated
## from the VA trial CSP 563 (docs/69, Sections 4 and 9, items 3 and 4).
## Reuses the machinery of script 22 and adds three response structures:
##   ar0.7      separable AR(1), rho = 0.7 per week, no measurement error
##              (the extended DGP)
##   csp        exponential decay 0.975 per week plus 10% measurement
##              error, applied to every response component (the CSP 563
##              placebo-arm estimate for the total score; the component
##              split is an assumption)
##   csp_grow   csp, with every component SD scaled by g(w) so the SD of
##              the response sum follows the CSP 563 placebo-arm SD of
##              change (17.0, 19.4, 22.1, 24.4 at weeks 6, 10, 18, 26;
##              linear fit), against the DGP's constant 16.0
## The other items of docs/69 Section 9 do not enter this closed form:
## the stratified E9 absorbs the TV, PB and BR means in the sequence
## intercepts, and its power does not depend on the biomarker's scale.
## Dropout is not modeled.
##
## sigma_v 7.5 is near the CSP 563 upper bound on the SD of individual
## drug effects (about 8 at week 10; with c_bm sd_BR = 2 the total is
## sqrt(7.5^2 + 2^2) = 7.8).
##
## Run from anywhere:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/23-csp563-covariance-designs.R \
##     [--mc-reps R] [--cores K]
suppressMessages({library(data.table); library(MASS); library(parallel)})
root <- here::here()
src22 <- readLines(file.path(root, 'analysis/scripts/quick-sim/carryover-closed-form/22-proportional-moderation-designs.R'))
eval(parse(text = src22[1:(grep('^grid <- CJ', src22)[1] - 1)]))
mc_from <- grep('^mc_cell <- function', src22)
mc_to <- mc_from - 1 + which(src22[mc_from:length(src22)] == '}')[1]
eval(parse(text = src22[mc_from:mc_to]))
mc_reps <- arg('--mc-reps', 4000L)
out_dir <- 'analysis/data/quick-sim/carryover-closed-form/csp563-covariance'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

STRUCTS <- c('CS', 'nug0.2', 'ar0.7', 'csp', 'csp_grow')
SIGMA_V <- c(0, 2, 7.5)
PHI_CSP <- 0.975
NU_CSP <- 0.10
csp_sd <- data.frame(w = c(6, 10, 18, 26), sd = c(17.0, 19.4, 22.1, 24.4))
sd_fit <- coef(lm(sd ~ w, data = csp_sd))
dgp_sd_sum <- sqrt(sum(outer(c(10, 5, sd_BR), c(10, 5, sd_BR)) * K3))
stopifnot(abs(dgp_sd_sum - 16.0) < 0.05)
growth <- function(w) (sd_fit[[1]] + sd_fit[[2]] * w) / dgp_sd_sum

resp_cor <- function(pi, structure) {
  dw <- abs(outer(pi$w, pi$w, '-'))
  switch(structure,
    CS = dgp$build(pi, 'A', 0, 0)$R[-(1:2), -(1:2)],
    ar0.7 = kronecker(K3, 0.7^dw),
    csp = , csp_grow = kronecker(K3, (1 - NU_CSP) * PHI_CSP^dw + NU_CSP * diag(n)),
    {
      nu <- as.numeric(sub('nug', '', structure))
      kronecker(K3, (1 - nu) * dgp$RHO^dw + nu * diag(n))
    })
}
resp_cov <- function(pi, structure) {
  g <- if (structure == 'csp_grow') growth(pi$w) else rep(1, n)
  sds <- c(rep(10, n), 10 * pi$e, rep(sd_BR, n)) * rep(g, 3)
  L <- cbind(diag(n), diag(n), diag(n))
  L %*% (outer(sds, sds) * resp_cor(pi, structure)) %*% t(L)
}

grid <- CJ(design = names(DES), rule = RULES_P, t_half = HALVES,
           structure = STRUCTS, sigma_v = SIGMA_V, sorted = FALSE)
res <- rbindlist(mclapply(seq_len(nrow(grid)), function(k)
  eval_cell(grid$design[k], grid$rule[k], grid$t_half[k], grid$structure[k],
            grid$sigma_v[k]), mc.cores = cores))
stopifnot(nrow(res) == nrow(grid), all(res$power > 0 & res$power < 1))
fwrite(res, file.path(out_dir, 'closed-form-power.csv'))

## Validation 1: script 22's structures are unchanged (docs/39: CS,
## current rule, no carryover, Hybrid 0.663, CO 0.765).
old <- fread(file.path(root, 'analysis/data/quick-sim/carryover-closed-form/proportional-moderation/closed-form-power.csv'))
old <- old[structure %in% c('CS', 'nug0.2')]
cmp <- merge(res[structure %in% c('CS', 'nug0.2') & sigma_v %in% c(0, 2)],
             old, by = c('design', 'rule', 't_half', 'structure', 'sigma_v'))
stopifnot(nrow(cmp) == nrow(old), max(abs(cmp$power.x - cmp$power.y)) < 1e-10)
cat(sprintf('Validation 1: %d cells match script 22\n', nrow(cmp)))

## Validation 2: visit-level Monte Carlo (script 22's mc_cell, which
## draws from resp_cov and so uses the structures defined here).
mc_grid <- CJ(design = c('CO', 'Hybrid'), rule = c('current', 'prop_zero'),
              t_half = c(0.2, 1), structure = c('csp', 'csp_grow'),
              sigma_v = c(0, 7.5), sorted = FALSE)
mc_grid[, seed := 20261031L + .I]
mc <- rbindlist(mclapply(seq_len(nrow(mc_grid)), function(k) {
  g <- mc_grid[k]
  data.table(g, mc_power = mc_cell(g$design, g$rule, g$t_half, g$structure,
                                   g$sigma_v, mc_reps, g$seed))
}, mc.cores = cores))
mc <- merge(mc, res, by = c('design', 'rule', 't_half', 'structure', 'sigma_v'))
mc[, mcse := sqrt(mc_power * (1 - mc_power) / mc_reps)]
mc[, z := (mc_power - power) / mcse]
fwrite(mc, file.path(out_dir, 'monte-carlo-check.csv'))
cat(sprintf('Validation 2: %d cells, %d reps each; |z| max %.2f, mean z %+.2f, %d beyond 2\n',
            nrow(mc), mc_reps, max(abs(mc$z)), mean(mc$z), sum(abs(mc$z) > 2)))

## Summaries.
cat(sprintf('\nGrowth factor g(w): %s at weeks 2.5, 10, 20\n',
            paste(round(growth(c(2.5, 10, 20)), 3), collapse = ', ')))
noise <- res[rule == 'current' & t_half == 0 & sigma_v == 0,
             .(design, structure, aCa = round(aCa_mean, 1))]
cat('\nMean contrast variance a\'Ca (current rule, no carryover)\n')
print(dcast(noise, design ~ structure, value.var = 'aCa'))
for (sv in SIGMA_V) {
  x <- res[rule == 'current' & sigma_v == sv,
           .(design, structure, t_half, p = round(power, 3))]
  cat(sprintf('\nPower, current rule, sigma_v %g (columns: t1/2 in weeks)\n', sv))
  print(dcast(x, structure + design ~ t_half, value.var = 'p'))
}
hc <- dcast(res[design %in% c('Hybrid', 'CO')],
            rule + structure + sigma_v + t_half ~ design, value.var = 'power')
hc[, ratio := round(Hybrid / CO, 2)]
fwrite(hc, file.path(out_dir, 'hybrid-vs-co.csv'))
cat('\nHybrid / CO power ratio, sigma_v 0 (columns: t1/2)\n')
print(dcast(hc[sigma_v == 0], rule + structure ~ t_half, value.var = 'ratio'))
rank_tab <- res[, .(best = design[which.max(power)],
                    order = paste(design[order(-power)], collapse = ' > ')),
                by = .(rule, structure, sigma_v, t_half)]
fwrite(rank_tab, file.path(out_dir, 'design-ranking.csv'))
cat('\nBest design by condition, sigma_v 0\n')
print(dcast(rank_tab[sigma_v == 0], structure + t_half ~ rule, value.var = 'best'))
flips <- rank_tab[, .(n_orders = uniqueN(order)), by = .(rule, sigma_v, t_half)]
cat('\nConditions where the full ordering differs across structures:',
    sum(flips$n_orders > 1), 'of', nrow(flips), '\n')
