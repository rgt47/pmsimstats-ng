# Toy analysis ladder for docs/67 (undergraduate tutorial).
#
# A deliberately simple data-generating process, not the extended DGP
# of papers 01 and 02: a random intercept, a linear time trend, AR(1)
# residuals in calendar weeks, and a drug effect whose size depends on
# the biomarker and decays after discontinuation. It exists so that
# students can watch each rung of the ladder change Type I error and
# power. Its numbers are not the paper 02 numbers.
#
# Usage: Rscript analysis/scripts/tutorial-g9/01-toy-ladder.R [n_reps]

suppressPackageStartupMessages({
  library(nlme)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args)) as.integer(args[1]) else 1000
n_cores <- max(1, parallelly::availableCores() - 1)

weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)

# Drug state and time since discontinuation of the four Hybrid paths
# (docs/40, Table "Drug state, t_od and t_sd by visit").
paths <- list(
  list(on = c(1, 1, 1, 1, 0, 0, 1, 0), tsd = c(0, 0, 0, 0, 1, 2, 0, 4)),
  list(on = c(1, 1, 1, 1, 0, 0, 0, 1), tsd = c(0, 0, 0, 0, 1, 2, 6, 0)),
  list(on = c(1, 1, 1, 0, 0, 0, 1, 0), tsd = c(0, 0, 0, 1, 2, 3, 0, 4)),
  list(on = c(1, 1, 1, 0, 0, 0, 0, 1), tsd = c(0, 0, 0, 1, 2, 3, 7, 0))
)

# Every path starts on drug, so no visit precedes first exposure and
# the pre-exposure rule (Dbc = 0) is never needed here.
exposure_weight <- function(on, tsd, half_life) {
  ifelse(on == 1, 1, 0.5^(tsd / half_life))
}

simulate_trial <- function(n, beta_int, true_half_life = 1,
                           beta_drug = -5, beta_time = -0.2,
                           sd_patient = 4, sd_resid = 3, rho = 0.7) {
  path_id <- rep_len(1:4, n)
  lag <- abs(outer(weeks, weeks, '-'))
  chol_r <- chol(sd_resid^2 * rho^lag)
  rows <- map(seq_len(n), function(i) {
    p <- paths[[path_id[i]]]
    b <- rnorm(1)
    expo <- exposure_weight(p$on, p$tsd, true_half_life)
    e <- drop(rnorm(length(weeks)) %*% chol_r)
    y <- 50 + rnorm(1, 0, sd_patient) + beta_time * weeks +
      (beta_drug + beta_int * b) * expo + e
    data.frame(ptID = i, t = weeks, on = p$on, tsd = p$tsd,
               bm = b, Sx = y)
  })
  dat <- do.call(rbind, rows)
  dat$bm <- dat$bm - mean(dat$bm[!duplicated(dat$ptID)])
  dat
}

p_model <- function(fit, term) summary(fit)$tTable[term, 'p-value']

p_cr2 <- function(fit, term, cluster) {
  ct <- clubSandwich::coef_test(fit, vcov = 'CR2', cluster = cluster,
                                test = 'Satterthwaite')
  ct$p_Satt[ct$Coef == term]
}

fit_lme <- function(form, dat, corr, method = 'REML') {
  nlme::lme(form, random = ~ 1 | ptID, correlation = corr, data = dat,
            method = method,
            control = nlme::lmeControl(opt = 'optim'))
}

analyze <- function(dat, true_half_life = 1,
                    grid = c(0.25, 0.5, 1, 2)) {
  car1 <- nlme::corCAR1(form = ~ t | ptID)
  out <- c()

  # Rung A: ordinary least squares, every row treated as independent.
  ols <- lm(Sx ~ bm + t + on + bm:on, data = dat)
  out['ols'] <- summary(ols)$coefficients['bm:on', 'Pr(>|t|)']

  # Rung B: random intercept, independent residuals.
  ri <- fit_lme(Sx ~ bm + t + on + bm:on, dat, NULL)
  out['ri'] <- p_model(ri, 'bm:on')

  # Rung C (G1): random intercept plus continuous-time AR(1).
  g1 <- fit_lme(Sx ~ bm + t + on + bm:on, dat, car1)
  out['g1'] <- p_model(g1, 'bm:on')

  # Rung D (G6): G1 with CR2 inference.
  out['g6'] <- p_cr2(g1, 'bm:on', dat$ptID)

  # Rung E (G8): exposure-weighted predictor at the true half-life.
  dat$Dbc <- exposure_weight(dat$on, dat$tsd, true_half_life)
  g3 <- fit_lme(Sx ~ bm + t + Dbc + bm:Dbc, dat, car1)
  out['g3'] <- p_model(g3, 'bm:Dbc')
  out['g8'] <- p_cr2(g3, 'bm:Dbc', dat$ptID)

  # Rung F (G9): half-life chosen by ML AIC, REML refit, CR2.
  aics <- map_dbl(grid, function(h) {
    dat$Dbc <- exposure_weight(dat$on, dat$tsd, h)
    AIC(fit_lme(Sx ~ bm + t + Dbc + bm:Dbc, dat, car1, 'ML'))
  })
  h_hat <- grid[which.min(aics)]
  dat$Dbc <- exposure_weight(dat$on, dat$tsd, h_hat)
  g4 <- fit_lme(Sx ~ bm + t + Dbc + bm:Dbc, dat, car1)
  out['g9'] <- p_cr2(g4, 'bm:Dbc', dat$ptID)
  out['h_hat'] <- h_hat
  out
}

one_rep <- function(seed, n, beta_int) {
  set.seed(seed)
  res <- tryCatch(analyze(simulate_trial(n, beta_int)),
                  error = function(e) NULL)
  if (is.null(res)) return(NULL)
  c(seed = seed, beta_int = beta_int, res)
}

run_cell <- function(n, beta_int, seeds) {
  rows <- parallel::mclapply(seeds, one_rep, n = n, beta_int = beta_int,
                             mc.cores = n_cores)
  do.call(rbind, compact(rows))
}

seeds <- 20261030 + seq_len(n_reps)
cells <- list(null = run_cell(40, 0, seeds),
              alt = run_cell(40, -0.8, seeds))

methods <- c('ols', 'ri', 'g1', 'g6', 'g3', 'g8', 'g9')
summary_tab <- do.call(rbind, map(names(cells), function(k) {
  m <- cells[[k]]
  rej <- colMeans(m[, methods] < 0.05)
  data.frame(cell = k, method = methods, n_ok = nrow(m),
             reject = round(rej, 3),
             mcse = round(sqrt(rej * (1 - rej) / nrow(m)), 3))
}))
rownames(summary_tab) <- NULL
print(summary_tab)
cat('\nSelected half-life under the alternative:\n')
print(table(cells$alt[, 'h_hat']))

out_dir <- here::here('analysis', 'data', 'tutorial-g9')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_tab, file.path(out_dir, 'toy-ladder-summary.csv'),
          row.names = FALSE)
saveRDS(cells, file.path(out_dir, 'toy-ladder-reps.rds'))
