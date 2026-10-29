## 21-summarise-paired.R
##
## Summaries of the paired paper 09 run written by 20-paired-dropout.R.
## Run from the repository root:
##   Rscript analysis/scripts/informative-dropout-by-design/21-summarise-paired.R
##
## Outputs (analysis/data/informative-dropout-by-design/):
##   paired-cells.csv     per design x regime x condition x analysis x
##                        test: rejection rate, bias against the true
##                        value, empirical SD, mean SE, kappa, 95%
##                        coverage, each with its MCSE
##   paired-loss.csv      alt regime, each dropout condition against
##                        'none' on the same datasets: power loss,
##                        coefficient shift and SE ratio with paired
##                        MCSE (and the unpaired MCSE for comparison)
##   paired-contrast.csv  differences in power loss between designs
##                        (independent datasets, SEs add in quadrature)
##   paired-dropout.csv   dropout diagnostics per design x regime x
##                        condition

suppressPackageStartupMessages(library(data.table))

dir <- Sys.getenv('P09_OUT', 'analysis/data/informative-dropout-by-design')
tag <- Sys.getenv('P09_TAG', 'paired')
reps <- readRDS(file.path(dir, paste0(tag, '-replicates.rds')))
drop <- readRDS(file.path(dir, paste0(tag, '-dropout.rds')))
stopifnot(!('error' %in% names(reps)) || all(is.na(reps$error)))

long <- rbind(
  reps[, .(design, regime, rep, condition, analysis, test = 'mb',
           est, se = se_mb, df = df_mb, p = p_mb, true_beta, car1,
           fit_ok)],
  reps[, .(design, regime, rep, condition, analysis, test = 'cr2',
           est, se = se_cr2, df = df_cr2, p = p_cr2, true_beta, car1,
           fit_ok)]
)
long[, rej := as.numeric(p < 0.05)]
long[, cover := as.numeric(abs(est - true_beta) <=
                             qt(0.975, df) * se)]

mcse_p <- function(x) sqrt(mean(x) * (1 - mean(x)) / length(x))

cells <- long[, {
  ok <- !is.na(p) & !is.na(se)
  n <- sum(ok)
  e <- est[ok]
  list(n_units = .N, n_ok = n, car1_rate = mean(car1[ok]),
       rej = mean(rej[ok]), mcse_rej = mcse_p(rej[ok]),
       mean_est = mean(e), true_beta = true_beta[1],
       bias = mean(e) - true_beta[1], mcse_bias = sd(e) / sqrt(n),
       emp_sd = sd(e), mcse_emp_sd = sd(e) / sqrt(2 * (n - 1)),
       mean_se = mean(se[ok]),
       kappa = mean(se[ok]) / sd(e),
       coverage = mean(cover[ok]), mcse_coverage = mcse_p(cover[ok]))
}, by = .(design, regime, condition, analysis, test)]

base <- long[condition == 'none',
             .(design, regime, rep, analysis, test,
               rej0 = rej, est0 = est, se0 = se)]
pair <- merge(long[condition != 'none'], base,
              by = c('design', 'regime', 'rep', 'analysis', 'test'))
pair <- pair[!is.na(p) & !is.na(rej0)]
loss <- pair[regime == 'alt', {
  n <- .N
  d <- rej0 - rej
  list(n = n, power_none = mean(rej0), power = mean(rej),
       loss = mean(d), mcse_loss = sd(d) / sqrt(n),
       mcse_loss_unpaired = sqrt(var(rej0) / n + var(rej) / n),
       shift = mean(est - est0), mcse_shift = sd(est - est0) / sqrt(n),
       se_ratio = mean(se / se0), mcse_se_ratio = sd(se / se0) / sqrt(n),
       sd_ratio = sd(est) / sd(est0))
}, by = .(design, condition, analysis, test)]

size_shift <- pair[regime == 'null', .(
  size_change = mean(rej - rej0),
  mcse_size_change = sd(rej - rej0) / sqrt(.N)
), by = .(design, condition, analysis, test)]

des <- unique(loss$design)
pairs <- combn(sort(des), 2, simplify = FALSE)
contrast <- rbindlist(lapply(pairs, function(pp) {
  a <- loss[design == pp[1]]
  b <- loss[design == pp[2]]
  m <- merge(a, b, by = c('condition', 'analysis', 'test'),
             suffixes = c('_a', '_b'))
  m[, .(design_a = pp[1], design_b = pp[2], condition, analysis, test,
        loss_a, loss_b, diff = loss_a - loss_b,
        se_diff = sqrt(mcse_loss_a^2 + mcse_loss_b^2),
        z = (loss_a - loss_b) / sqrt(mcse_loss_a^2 + mcse_loss_b^2),
        se_ratio_a, se_ratio_b)]
}))

dsum <- drop[, {
  m <- function(x) mean(x, na.rm = TRUE)
  s <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  list(n = .N,
       ever_drop = m(ever_drop), obs_missing = m(obs_missing),
       frac_drop_offdrug = m(frac_drop_offdrug),
       cor_drop_improve = m(cor_drop_improve),
       mcse_cor_drop_improve = s(cor_drop_improve),
       cor_drop_bm = m(cor_drop_bm), mcse_cor_drop_bm = s(cor_drop_bm),
       drop_p1 = m(drop_p1), drop_p2 = m(drop_p2),
       drop_p3 = m(drop_p3), drop_p4 = m(drop_p4))
}, by = .(design, regime, condition)]

fwrite(cells, file.path(dir, paste0(tag, '-cells.csv')))
fwrite(merge(loss, size_shift, by = c('design', 'condition', 'analysis',
                                      'test'), all.x = TRUE),
       file.path(dir, paste0(tag, '-loss.csv')))
fwrite(contrast, file.path(dir, paste0(tag, '-contrast.csv')))
fwrite(dsum, file.path(dir, paste0(tag, '-dropout.csv')))
cat('Wrote summaries to', dir, '\n')
