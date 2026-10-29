library(tinytest)

#---------------------------------------------------------------------
# Engine changes of version 0.1.0 (paper 05 referee report, 2026-10-28)
#---------------------------------------------------------------------

# Power decay form: a non-integer exponent must not give NaN
phi <- nof1power::carryover_decay(
  c(0, 0.5, 1, 2, 5), halflife = 0.5, form = 'power', power = 1.5
)
expect_false(anyNA(phi), info = 'power form returns no NaN past max_time')
expect_equal(phi[1], 1, info = 'multiplier is 1 at tsd = 0')
expect_equal(phi[5], 0, info = 'multiplier is 0 beyond max_time')

# Parallel-group analysis: ANCOVA arm contrast
td_rct <- nof1power::build_trial_design(
  name_longform = 'RCT', name_shortform = 'RCT',
  timepoints = c(4, 8, 12), expectancies = rep(0.5, 3),
  ondrug = list(rep(1, 3), rep(0, 3))
)
resp <- data.frame(
  cat = c('tv', 'pb', 'br'), max = c(0.5, 1, 0), disp = 1,
  rate = 0.3, sd = c(0.5, 0.5, 1), stringsAsFactors = FALSE
)
bl <- data.frame(cat = c('bm', 'BL'), m = c(120, 8), sd = c(15, 1.5),
                 stringsAsFactors = FALSE)
mp <- data.frame(
  N = 40, c.bm = 0.3, carryover_t1half = 0, c.tv = 0.3, c.pb = 0.3,
  c.br = 0.3, c.cf1t = 0.1, c.cfct = 0.05, stringsAsFactors = FALSE
)
simp <- list(Nreps = 300, progressiveSave = FALSE, basesavename = 't',
             nRep2save = .Machine$integer.max, saveunit2start = 1,
             savedir = tempdir())
ap <- data.frame(
  useDE = TRUE, t_random_slope = FALSE, full_model_out = FALSE,
  carryover_t1half = 0, simplecarryover = FALSE,
  carryover_scalefactor = 1, target_coef = 'main',
  stringsAsFactors = FALSE
)
set.seed(11)
res <- nof1power::generate_simulated_results(
  list(td_rct), list(list(name = 'null', param = resp, target_delta = 0)),
  list(list(name = 'bl', param = bl)), NA, mp, simp, ap,
  n_cores = 1
)
r <- res$results
mc_se <- stats::sd(r$beta) / sqrt(nrow(r))
expect_true(abs(mean(r$beta)) < 3 * mc_se,
            info = 'ANCOVA arm contrast has mean zero under the null')
expect_true(
  abs(mean(r$betaSE) / stats::sd(r$beta) - 1) < 0.1,
  info = 'ANCOVA model SE within 10% of the empirical SD (kappa)'
)
expect_equal(unique(r$delta_on), 0, info = 'null estimand is zero')
expect_false(any(r$pd_repaired), info = 'repair flag is stored and FALSE')
expect_true(all(grepl('OLS', r$fit_formula)),
            info = 'fitted formula is recorded')

# Infeasible cells are skipped and reported, never repaired
mp_bad <- rbind(mp, transform(mp, c.cf1t = 0.95, c.cfct = 0.9))
res2 <- nof1power::generate_simulated_results(
  list(td_rct), list(list(name = 'null', param = resp)),
  list(list(name = 'bl', param = bl)), NA, mp_bad,
  modifyList(simp, list(Nreps = 3)), ap,
  n_cores = 1, skip_infeasible = TRUE
)
expect_equal(nrow(res2$infeasible), 1L,
             info = 'one infeasible cell reported')
expect_equal(unique(res2$results$modelparamset), 1L,
             info = 'only the feasible cell ran')
expect_error(
  nof1power::generate_simulated_results(
    list(td_rct), list(list(name = 'null', param = resp)),
    list(list(name = 'bl', param = bl)), NA, mp_bad,
    modifyList(simp, list(Nreps = 3)), ap, n_cores = 1
  ),
  info = 'without skip_infeasible a non-PD cell stops the run'
)

# Matched estimand: target_delta sets the ceiling per design
td_alt <- nof1power::build_trial_design(
  name_longform = 'Alt', name_shortform = 'Alt',
  timepoints = c(4, 8, 12), expectancies = rep(0.5, 3),
  ondrug = list(c(1, 0, 1), c(0, 1, 0))
)
res3 <- nof1power::generate_simulated_results(
  list(td_rct, td_alt),
  list(list(name = 'm', param = transform(resp, max = c(0.5, 1, 2)),
            target_delta = -0.5)),
  list(list(name = 'bl', param = bl)), NA, mp,
  modifyList(simp, list(Nreps = 2)), ap, n_cores = 1
)
expect_equal(unique(round(res3$results$delta_on, 8)), -0.5,
             info = 'every design has the target delta_on')
expect_true(length(unique(round(res3$results$br_max, 6))) == 2,
            info = 'ceiling differs between designs')

# Treatment heterogeneity: tau adds variance to the drug arm only
set.seed(5)
d0 <- nof1power::generate_data(
  transform(mp, N = 4000), transform(resp, max = c(0.5, 1, 2)), bl,
  nof1power:::ensure_timepoint_name(td_rct$trialpaths[[1]]),
  empirical = FALSE, make_positive_definite = TRUE
)
d1 <- nof1power::generate_data(
  transform(mp, N = 4000, tau = 1), transform(resp, max = c(0.5, 1, 2)),
  bl, nof1power:::ensure_timepoint_name(td_rct$trialpaths[[1]]),
  empirical = FALSE, make_positive_definite = TRUE
)
expect_true(stats::var(d1$V3.br) > stats::var(d0$V3.br) + 0.5,
            info = 'tau > 0 raises the variance of the drug response')
expect_true(abs(mean(d1$V3.br) - mean(d0$V3.br)) < 0.1,
            info = 'tau leaves the mean drug response unchanged')

# use_corCAR1 = FALSE fits independent residuals and says so
td_alt3 <- nof1power::build_trial_design(
  name_longform = 'Alt', name_shortform = 'Alt',
  timepoints = c(4, 8, 12, 16), expectancies = rep(0.5, 4),
  ondrug = list(c(1, 0, 1, 0), c(0, 1, 0, 1))
)
ap2 <- transform(ap, use_corCAR1 = FALSE)
set.seed(3)
res4 <- nof1power::generate_simulated_results(
  list(td_alt3), list(list(name = 'm', param = resp)),
  list(list(name = 'bl', param = bl)), NA, mp,
  modifyList(simp, list(Nreps = 3)), ap2, n_cores = 1
)
expect_true(all(grepl('no serial correlation', res4$results$fit_formula)),
            info = 'fitted formula records the absence of corCAR1')
