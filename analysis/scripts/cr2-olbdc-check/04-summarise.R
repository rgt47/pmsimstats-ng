## Summarise every production rds in analysis/data/cr2-olbdc-check/:
## empirical SD, mean model-based and CR2 SE, CR2 df quantiles,
## rejection rates at 0.05 with MCSE and exact two-sided binomial
## p-values against 0.05, and paired differences between variants.
suppressPackageStartupMessages(library(data.table))
dir <- 'analysis/data/cr2-olbdc-check'
fs <- list.files(dir, pattern = '^(p09|p02|iso_[a-z]+)-.*-N[0-9]+(-from[0-9]+)?\\.rds$',
                 full.names = TRUE)
d <- rbindlist(lapply(fs, readRDS), fill = TRUE)
if ('error' %in% names(d)) {
  print(d[!is.na(error), .N, by = .(setting, design, N, error)])
  d <- d[is.na(error)]
}
binp <- function(x, n) binom.test(x, n, 0.05)$p.value
s <- d[fit_ok == TRUE, {
  k_cr2 <- sum(p_cr2 < 0.05, na.rm = TRUE); n_cr2 <- sum(!is.na(p_cr2))
  k_mb <- sum(p_mb < 0.05); n_mb <- .N
  q <- quantile(df_cr2, c(0.05, 0.5, 0.95), na.rm = TRUE)
  .(n = .N, car1_frac = mean(car1),
    emp_sd = sd(est), mean_se_mb = mean(se_mb),
    mean_se_cr2 = mean(se_cr2, na.rm = TRUE),
    ratio_mb = mean(se_mb) / sd(est),
    ratio_cr2 = mean(se_cr2, na.rm = TRUE) / sd(est),
    df_q05 = q[[1]], df_q50 = q[[2]], df_q95 = q[[3]],
    rej_mb = k_mb / n_mb, mcse_mb = sqrt(k_mb / n_mb *
                                         (1 - k_mb / n_mb) / n_mb),
    binp_mb = binp(k_mb, n_mb),
    rej_cr2 = k_cr2 / n_cr2, mcse_cr2 = sqrt(k_cr2 / n_cr2 *
                                             (1 - k_cr2 / n_cr2) / n_cr2),
    binp_cr2 = binp(k_cr2, n_cr2))
}, by = .(setting, design, N, variant)]
setorder(s, design, setting, N, variant)
fwrite(s, file.path(dir, 'summary.csv'))
print(s[, .(setting, design, N, variant, n, car1_frac,
            emp_sd = round(emp_sd, 3), se_mb = round(mean_se_mb, 3),
            se_cr2 = round(mean_se_cr2, 3),
            r_mb = round(ratio_mb, 3), r_cr2 = round(ratio_cr2, 3))])
print(s[, .(setting, design, N, variant,
            df = sprintf('%.1f/%.1f/%.1f', df_q05, df_q50, df_q95),
            rej_mb = sprintf('%.4f (%.4f) p=%.3g', rej_mb, mcse_mb, binp_mb),
            rej_cr2 = sprintf('%.4f (%.4f) p=%.3g', rej_cr2, mcse_cr2,
                              binp_cr2))])

w <- dcast(d[fit_ok == TRUE & !is.na(p_cr2)],
           setting + design + N + rep ~ variant,
           value.var = 'p_cr2')
pairs <- list(c('ta0', 'g6form'), c('ta0', 'ta0_noDe'),
              c('ta0', 'ta0_noBL'), c('g6', 'p09form'))
pd <- rbindlist(lapply(pairs, function(p) {
  if (!all(p %in% names(w))) return(NULL)
  x <- w[!is.na(get(p[1])) & !is.na(get(p[2]))]
  x[, .(a = p[1], b = p[2], n = .N,
        diff = mean((get(p[1]) < 0.05) - (get(p[2]) < 0.05)),
        se = sd((get(p[1]) < 0.05) - (get(p[2]) < 0.05)) / sqrt(.N)),
    by = .(setting, design, N)]
}))
if (nrow(pd)) {
  print(pd)
  fwrite(pd, file.path(dir, 'paired-differences.csv'))
}
