## Which design characteristics drive the power of the paired-difference
## interaction test? Uses the closed form of 01-e9-closed-form.R (its
## definitions are evaluated, nothing is simulated).
##
## Without carryover the interaction slope is the same in every design,
## so power depends on the design only through the noise of the on-off
## contrast, Var(Delta) = a' C a, and (with carryover) through the
## off-drug coupling that survives, gbar_off. The script
##   1. decomposes Var(Delta) under CS into its terms by design path;
##   2. reports the AR(1) time part a' A a and the expectancy confound
##      sum(a e) by design path;
##   3. evaluates hypothetical designs that change one characteristic at
##      a time, under CS and separable AR(1), graded coupling, for
##      t1/2 in {0, 0.5, 1, 2}, reporting the c_bm ceiling, the per-patient
##      signal-to-noise ratio and power at c_bm = 0.25 (E9 as specified,
##      and with the path term removed, as path intercepts would).
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/02-design-characteristics.R
suppressMessages(library(data.table))
src <- readLines('analysis/scripts/quick-sim/carryover-closed-form/01-e9-closed-form.R')
stop_at <- grep('^## ---------------- per-path table', src)[1] - 1
eval(parse(text = src[1:stop_at]))
out_dir <- 'analysis/data/quick-sim/carryover-closed-form'

## ---------------- hypothetical designs ----------------
## Every design has 8 post-baseline visits over 20 weeks, so only the
## arrangement changes. Expectancy 1 marks open-label visits.
hyb_t <- c(4, 8, 9, 10, 11, 12, 16, 20)
eq_t <- seq(2.5, 20, by = 2.5)
add_design <- function(name, times, e, paths) {
  dgp$designs[[name]] <- list(t_wk = diff(c(0, times)), e = e, paths = paths)
}
add_design('Hybrid, blinded', hyb_t, rep(0.5, 8), dgp$designs$Hybrid$paths)
add_design('CO, open-label first period', eq_t, c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
           dgp$designs$CO$paths)
add_design('ABAB (2 crossovers)', eq_t, rep(0.5, 8),
           list(c(1, 1, 0, 0, 1, 1, 0, 0), c(0, 0, 1, 1, 0, 0, 1, 1)))
add_design('ABAB, active first only', eq_t, rep(0.5, 8),
           list(c(1, 1, 0, 0, 1, 1, 0, 0)))
add_design('Alternating every visit', eq_t, rep(0.5, 8),
           list(c(1, 0, 1, 0, 1, 0, 1, 0), c(0, 1, 0, 1, 0, 1, 0, 1)))
add_design('Hybrid timing, ABAB blocks', hyb_t, rep(0.5, 8),
           list(c(1, 1, 0, 0, 1, 1, 0, 0), c(0, 0, 1, 1, 0, 0, 1, 1)))
designs_eval <- c('Hybrid', 'Hybrid, blinded', 'OL+BDC', 'CO',
                  'CO, open-label first period', 'ABAB (2 crossovers)',
                  'ABAB, active first only', 'Alternating every visit',
                  'Hybrid timing, ABAB blocks')

## ---------------- 1-2. structure of the contrast noise ----------------
## CS terms of Var(Delta) (docs/37, Section 4.1), E9.
cs_terms <- function(pi) {
  a <- weights(pi, 'E9'); e <- pi$e
  rho <- dgp$RHO; c1 <- dgp$C1; cx <- dgp$CX
  s_tv <- 10; s_pb <- 10; s_br <- sd_BR
  A2 <- sum(a^2); Ae <- sum(a * e); Aee <- sum(a^2 * e^2); A2e <- sum(a^2 * e)
  data.table(
    balance = (1 - rho) * (s_tv^2 + s_br^2) * A2,
    pb_occasion = (1 - rho) * s_pb^2 * Aee,
    expectancy_confound = rho * s_pb^2 * Ae^2,
    cross_factor = 2 * (c1 - cx) * (s_tv * s_br * A2 + (s_tv + s_br) * s_pb * A2e),
    sum_a2 = A2, sum_ae = Ae,
    a_A_a = drop(t(a) %*% (dgp$RHO^abs(outer(pi$w, pi$w, '-'))) %*% a))
}
struct <- rbindlist(lapply(designs_eval, function(dn) {
  use_design(dn)
  rbindlist(lapply(seq_len(nP), function(p) {
    pi <- paths[[p]]
    cbind(design = dn, path = p, n_on = sum(pi$on), n_off = sum(!pi$on),
          switches = sum(diff(c(0, pi$on)) != 0), cs_terms(pi),
          var_cs = path_moments(pi, 'step', 0, 'E9', construct = 'CS')$var_D,
          var_sep = path_moments(pi, 'step', 0, 'E9', construct = 'sep')$var_D)
  }))
}))
fwrite(struct, file.path(out_dir, 'design-contrast-structure.csv'))
cat('Contrast noise by design path (E9; CS terms sum to var_cs)\n')
print(struct[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

## ---------------- 3. power of the designs ----------------
snr <- function(m) {
  ## Per-patient squared contrast-biomarker correlation, patient-weighted,
  ## without the path term: the design's information per patient.
  w <- Ns / N
  sum(w * (m$gamma * sd_B)^2 / m$var_D)
}
res <- rbindlist(lapply(designs_eval, function(dn) {
  rbindlist(lapply(c('CS', 'sep'), function(cn) {
    rbindlist(lapply(c(0, 0.5, 1, 2), function(th) {
      ceil <- ceiling_of('graded', th, cn, dn)
      use_design(dn)
      m <- rbindlist(lapply(paths, path_moments, rule = 'graded', th = th,
                            variant = 'E9', construct = cn))
      ms <- copy(m)[, mu := 0]
      data.table(design = dn, construct = cn, t_half = th, ceiling = ceil,
                 snr = snr(m), path_term = sum(Ns / N * (m$mu - sum(Ns / N * m$mu))^2),
                 power = if (CBM < ceil) pooled(m)$power else NA_real_,
                 power_strat = if (CBM < ceil) pooled(ms)$power else NA_real_)
    }))
  }))
}))
use_design('Hybrid')
fwrite(res, file.path(out_dir, 'design-power.csv'))
cat('\nGraded coupling, c_bm = 0.25, E9: ceiling, per-patient SNR (rho^2), power',
    '(NA = above the ceiling), power with the path term removed\n')
print(dcast(res[, .(design, construct, t_half,
                    v = sprintf('c*%.3f snr%.4f p%s/%s', ceiling, snr,
                                ifelse(is.na(power), ' NA ', sprintf('%.3f', power)),
                                ifelse(is.na(power_strat), 'NA', sprintf('%.3f', power_strat))))],
            design + construct ~ t_half, value.var = 'v'))
