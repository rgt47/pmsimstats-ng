library(tinytest)

#---------------------------------------------------------------------
# Fixtures (pmsimstats short-form conventions)
#---------------------------------------------------------------------

trial_path <- data.frame(
  timepoint_name = paste0("V", 1:4),
  t_wk = c(2, 2, 2, 2),
  e = c(0.5, 0.5, 0, 0),
  tod = c(2, 4, 4, 4),
  tsd = c(0, 0, 2, 4),
  tpb = c(2, 4, 4, 4),
  stringsAsFactors = FALSE
)

model_param <- data.frame(
  N = 50, c.bm = 0.5, carryover_t1half = 4,
  c.tv = 0.6, c.pb = 0.6, c.br = 0.6,
  c.cf1t = 0.3, c.cfct = 0.1,
  stringsAsFactors = FALSE
)

resp_param <- data.frame(
  cat = c("tv", "pb", "br"),
  max = c(10, 5, 8), disp = c(1, 1, 1),
  rate = c(0.5, 0.5, 0.5), sd = c(2, 1, 3),
  stringsAsFactors = FALSE
)

baseline_param <- data.frame(
  cat = c("bm", "BL"),
  m = c(120, 50), sd = c(15, 10),
  stringsAsFactors = FALSE
)

#---------------------------------------------------------------------
# 1. match.arg rejects invalid architecture names
#---------------------------------------------------------------------

expect_error(
  nof1power::generate_data(
    model_param, resp_param, baseline_param, trial_path,
    FALSE, TRUE, dgp_architecture = "invalid"
  )
)

#---------------------------------------------------------------------
# 2. build_sigma_matrix returns correct structure
#---------------------------------------------------------------------

sig <- nof1power::build_sigma_matrix(
  model_param, resp_param, baseline_param, trial_path
)

expect_true(is.list(sig))
expect_true(all(c("sigma", "labels", "means", "chol_sigma",
                  "was_pd_corrected") %in% names(sig)))
expect_true(is.matrix(sig$sigma))
expect_equal(nrow(sig$sigma), ncol(sig$sigma))
expect_equal(length(sig$labels), nrow(sig$sigma))

#---------------------------------------------------------------------
# 3. MVN architecture: non-zero BM-BR correlations
#---------------------------------------------------------------------

sig_mvn <- nof1power::build_sigma_matrix(
  model_param, resp_param, baseline_param, trial_path,
  dgp_architecture = "mvn"
)

br_cols <- grep("\\.br$", sig_mvn$labels, value = TRUE)
bm_br_cors <- sig_mvn$correlations["bm", br_cols]
expect_true(any(bm_br_cors != 0))

#---------------------------------------------------------------------
# 4. Mean moderation: zero BM-BR correlations
#---------------------------------------------------------------------

sig_mm <- nof1power::build_sigma_matrix(
  model_param, resp_param, baseline_param, trial_path,
  dgp_architecture = "mean_moderation"
)

bm_br_cors_mm <- sig_mm$correlations["bm", br_cols]
expect_true(all(bm_br_cors_mm == 0))

#---------------------------------------------------------------------
# 5. Mean moderation: high-biomarker participants show higher BR
#    at on-drug timepoints
#---------------------------------------------------------------------

set.seed(42)
dat_mm <- nof1power::generate_data(
  model_param, resp_param, baseline_param, trial_path,
  FALSE, TRUE, seed = 42, dgp_architecture = "mean_moderation"
)

cor_bm_br <- cor(dat_mm$bm, dat_mm[["V1.br"]])
expect_true(cor_bm_br > 0.1)

#---------------------------------------------------------------------
# 6. lambda_cor auto-computed from carryover_t1half
#---------------------------------------------------------------------

sig_auto <- nof1power::build_sigma_matrix(
  model_param, resp_param, baseline_param, trial_path,
  lambda_cor = NA, dgp_architecture = "mvn"
)

expected_lambda <- log(2) / model_param$carryover_t1half
sig_manual <- nof1power::build_sigma_matrix(
  model_param, resp_param, baseline_param, trial_path,
  lambda_cor = expected_lambda, dgp_architecture = "mvn"
)

expect_equal(sig_auto$sigma, sig_manual$sigma)

#---------------------------------------------------------------------
# 7. Both architectures produce same column structure
#---------------------------------------------------------------------

set.seed(42)
dat_mvn_check <- nof1power::generate_data(
  model_param, resp_param, baseline_param, trial_path,
  FALSE, TRUE, seed = 42, dgp_architecture = "mvn"
)

expect_equal(names(dat_mvn_check), names(dat_mm))
expect_equal(nrow(dat_mvn_check), nrow(dat_mm))

#---------------------------------------------------------------------
# 8. Cached sigma produces identical results
#---------------------------------------------------------------------

cached <- nof1power::build_sigma_matrix(
  model_param, resp_param, baseline_param, trial_path,
  dgp_architecture = "mvn"
)

set.seed(99)
dat_fresh <- nof1power::generate_data(
  model_param, resp_param, baseline_param, trial_path,
  FALSE, TRUE, seed = 99, dgp_architecture = "mvn"
)

set.seed(99)
dat_cached <- nof1power::generate_data(
  model_param, resp_param, baseline_param, trial_path,
  FALSE, TRUE, seed = 99, dgp_architecture = "mvn",
  cached_sigma = cached
)

# Strip PD attributes before comparison (cached path sets
# was_pd_corrected = FALSE by design)
attributes(dat_fresh)[c("non_positive_definite_count",
                        "non_positive_definite_rate")] <- NULL
attributes(dat_cached)[c("non_positive_definite_count",
                         "non_positive_definite_rate")] <- NULL
expect_equal(dat_fresh, dat_cached)

#---------------------------------------------------------------------
# 9. Default parameters are backward-compatible
#---------------------------------------------------------------------

set.seed(123)
dat_default <- nof1power::generate_data(
  model_param, resp_param, baseline_param, trial_path,
  FALSE, TRUE, seed = 123
)

expect_true(is.data.frame(dat_default))
expect_true("ptID" %in% names(dat_default))
expect_true("bm" %in% names(dat_default))
expect_true("BL" %in% names(dat_default))
expect_equal(nrow(dat_default), model_param$N)

#---------------------------------------------------------------------
# 10. PD-stats attributes are recorded on fresh (non-cached) draws
#---------------------------------------------------------------------

expect_true(!is.null(attr(dat_default, "sigma_count")))
expect_true(!is.null(attr(dat_default,
                          "non_positive_definite_count")))
expect_true(!is.null(attr(dat_default,
                          "non_positive_definite_rate")))

#---------------------------------------------------------------------
# 11. carryover_decay has correct boundary behavior
#---------------------------------------------------------------------

# Exponential: phi(0) = 1, phi(halflife) = 0.5
expect_equal(nof1power::carryover_decay(0, halflife = 1), 1)
expect_equal(nof1power::carryover_decay(1, halflife = 1), 0.5)

# Linear: phi(0) = 1, phi(2 * halflife) = 0
expect_equal(
  nof1power::carryover_decay(0, halflife = 1, form = "linear"), 1
)
expect_equal(
  nof1power::carryover_decay(2, halflife = 1, form = "linear"), 0
)

# Weibull k = 1 recovers exponential
expect_equal(
  nof1power::carryover_decay(1, halflife = 1, form = "weibull",
                             shape = 1),
  nof1power::carryover_decay(1, halflife = 1, form = "exponential")
)
