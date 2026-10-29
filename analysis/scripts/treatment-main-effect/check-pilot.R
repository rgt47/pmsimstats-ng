# check-pilot.R
#
# Fast checks for the paper 04 pilot (about 30 seconds). Run from the
# project root:
#   Rscript analysis/scripts/treatment-main-effect/check-pilot.R
#
# Sections
#   1  Structural checks of the analysis rules (no leakage of the
#      truth into schedule-based rules; pre-exposure guard).
#   2  Closed form against the stored 2026-05-19 run (v1 driver).
#   3  Closed form against a short Monte Carlo of the v2 driver.
#   4  RNG: per-replicate streams are reproducible and distinct.
#
# Sections 1 and 4 are structural or exact checks. Sections 2 and 3 are
# numerical agreement checks with stated tolerances. None of them
# replaces the production run.

suppressPackageStartupMessages({
  library(nlme)
  library(clubSandwich)
})
source(file.path('analysis', 'scripts', 'treatment-main-effect',
                 'pilot-core.R'))

n_fail <- 0L
check <- function(ok, label) {
  cat(sprintf('%s  %s\n', if (isTRUE(ok)) 'PASS' else 'FAIL', label))
  if (!isTRUE(ok)) n_fail <<- n_fail + 1L
}

# ---- 1. Structure --------------------------------------------------
seqs <- all_sequences()
check(length(seqs) == 70L, 'there are 70 balanced 4-of-8 sequences')

first_off_counts <- vapply(seqs, function(d) {
  sum(exposure_structure(d)$first_off)
}, numeric(1))
check(abs(mean(first_off_counts) - 2) < 1e-12,
      'XS drops 2.000 periods on average over the 70 sequences')

keep_xs_h02 <- lapply(seqs, function(d) analysis_design(d, 0.2)$XS$keep)
keep_xs_h2 <- lapply(seqs, function(d) analysis_design(d, 2)$XS$keep)
check(identical(keep_xs_h02, keep_xs_h2),
      'XS retained rows do not depend on the true half-life')

keep_xo_h02 <- lapply(seqs, function(d) analysis_design(d, 0.2)$XO$keep)
keep_xo_h2 <- lapply(seqs, function(d) analysis_design(d, 2)$XO$keep)
check(!identical(keep_xo_h02, keep_xo_h2),
      'oracle XO retained rows depend on the true half-life (as labeled)')

check(!('beta_d' %in% names(formals(analysis_design))),
      'the oracle rule has no effect-size argument, so null and alternative cells keep the same rows')

pre <- c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 0L)
ad <- analysis_design(pre, 0.2)
check(all(ad$DBC$x[1:2] == 0),
      'Dbc is 0 before first exposure (pre-exposure guard)')
check(all(ad$DBC$x >= 0 & ad$DBC$x <= 1), 'Dbc lies in [0, 1]')
check(abs(ad$DBC$x[5] - carry_fraction(0, 0.2)) < 1e-12 &&
        abs(ad$DBC$x[6] - carry_fraction(1, 0.2)) < 1e-12,
      'Dbc off-drug values equal the carryover fraction at k = 0, 1')
check(identical(ad$DBCL$l1, c(0, 0, 0, 0, 1, 0, 0, 1)),
      'lag regressor flags the first off-drug period after each drug run')
check(abs(carry_fraction(0, 0.2) - 0.2796) < 1e-3,
      'period-mean carryover fraction is 0.280 at t_half = 0.2 week')
check(all(carry_fraction(0:3, 0) == 0), 't_half = 0 gives no carryover')

# ---- 2. Closed form against the stored v1 run -------------------------
old_file <- file.path('analysis', 'data', 'derived_data',
                      '04-mainfx-production.rds')
if (file.exists(old_file)) {
  old <- readRDS(old_file)$primary
  cf_null <- cf_nof1('ALL', 0, 0.2)
  cf_alt <- cf_nof1('XS', -1, 0.2)
  check(abs(cf_null$se - old$nof1_emp_se[old$true_effect == 0]) < 0.005,
        sprintf('null-cell N-of-1 SE: closed form %.4f, v1 run %.4f',
                cf_null$se, old$nof1_emp_se[old$true_effect == 0]))
  check(abs(cf_alt$se - old$nof1_emp_se[old$true_effect == -1]) < 0.01,
        sprintf('alternative-cell N-of-1 SE: closed form %.4f, v1 %.4f',
                cf_alt$se, old$nof1_emp_se[old$true_effect == -1]))
  rc <- cf_rct(-1, 1, 1)
  check(abs(rc$se - old$rct_emp_se[1]) < 0.03,
        sprintf('RCT SE: closed form %.4f, v1 empirical %.4f',
                rc$se, old$rct_emp_se[1]))
  for (b in c(-0.5, -1, -2)) {
    p_cf <- cf_nof1('XS', b, 0.2)$power
    p_old <- old$nof1_power_kr[old$true_effect == b]
    check(abs(p_cf - p_old) < 0.03,
          sprintf('N-of-1 power at %.1f: closed form %.3f, v1 %.3f',
                  b, p_cf, p_old))
    r_cf <- cf_rct(b, 1, 1)$power
    r_old <- old$rct_power[old$true_effect == b]
    check(abs(r_cf - r_old) < 0.04,
          sprintf('RCT power at %.1f: closed form %.3f, v1 %.3f',
                  b, r_cf, r_old))
  }
} else {
  cat('SKIP  stored v1 run not found\n')
}

# ---- 3. Closed form against a short Monte Carlo ------------------------
n_mc <- 60L
streams <- make_streams(2L * n_mc, 77L)
reps <- do.call(rbind, lapply(seq_len(n_mc), function(i) {
  run_replicate(i, streams[[i]], -1, 0.5)
}))
smry <- summarize_cell(reps, -1)
for (a in c('XS', 'DBC', 'DBCL', 'ALL')) {
  cf <- cf_nof1(a, -1, 0.5)
  s <- smry[smry$arm == a, ]
  check(abs(s$mean_se_c / cf$se - 1) < 0.08,
        sprintf('%s model SE %.3f vs closed form %.3f', a, s$mean_se_c,
                cf$se))
  check(abs(s$bias - cf$bias) < 4 * s$bias_mcse + 0.02,
        sprintf('%s bias %.3f vs closed form %.3f', a, s$bias, cf$bias))
}
n_rct <- 400L
sr <- make_streams(n_rct, 78L)
rct_se <- vapply(seq_len(n_rct), function(i) {
  assign('.Random.seed', sr[[i]], envir = globalenv())
  u <- rnorm(PILOT$n_patients, 0, PILOT$sigma_u)
  tw <- seq_len(8)
  cw <- chol(PILOT$sigma_w^2 * outer(tw, tw, function(a, b) {
    PILOT$rho^abs(a - b)
  }))
  c(fit_rct(generate_rct(u, -1, cw, 1, 1))$beta,
    fit_rct(generate_rct(u, -1, cw, PILOT$r_bl, 1))$beta,
    fit_rct(generate_rct(u, -1, cw, PILOT$r_bl, 8))$beta)
}, numeric(3))
emp <- apply(rct_se, 1L, sd)
cf_r <- c(cf_rct(-1, 1, 1)$se, cf_rct(-1, PILOT$r_bl, 1)$se,
          cf_rct(-1, PILOT$r_bl, 8)$se)
check(all(abs(emp / cf_r - 1) < 0.12),
      sprintf('RCT empirical SEs %s vs closed form %s',
              paste(round(emp, 3), collapse = '/'),
              paste(round(cf_r, 3), collapse = '/')))

# ---- 4. RNG ------------------------------------------------------------
a1 <- run_replicate(1, streams[[1]], -1, 0.5)
a2 <- run_replicate(1, streams[[1]], -1, 0.5)
check(isTRUE(all.equal(a1$beta, a2$beta)),
      'a replicate reproduces exactly from its stored stream')
check(!isTRUE(all.equal(reps$beta[reps$rep == 1 & reps$arm == 'XS'],
                        reps$beta[reps$rep == 2 & reps$arm == 'XS'])),
      'adjacent replicates use different streams')
cell2 <- run_replicate(1, streams[[n_mc + 1L]], -1, 0.5)
check(!isTRUE(all.equal(a1$beta, cell2$beta)),
      'the same replicate index in two cells draws different data')

cat(sprintf('\n%d check(s) failed\n', n_fail))
if (n_fail > 0L) quit(status = 1L)
