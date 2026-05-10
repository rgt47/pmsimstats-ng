## Paper 04 prototype: aggregated N-of-1 vs parallel-group RCT
## Compares power for the average treatment effect at matched N.
## Morris-compatible per-replicate output.

suppressPackageStartupMessages({
  library(devtools)
  library(data.table)
  library(ggplot2)
  library(nlme)
})

devtools::load_all('.', quiet = TRUE)

set.seed(20260507)

out_dir_data <- 'analysis/data/quick-sim'
out_dir_fig  <- 'analysis/figures/quick-sim'
dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_fig,  recursive = TRUE, showWarnings = FALSE)

## --------------------------------------------------------------- ##
## 1. N-of-1 hybrid trial design                                   ##
## --------------------------------------------------------------- ##

td_hybrid <- buildtrialdesign(
  name_longform  = 'hybrid',
  name_shortform = 'Hyb',
  timepoints     = c(2.5, 5, 7.5, 10, 12, 14, 16, 18, 20),
  timeptnames    = c(paste0('OL', 1:5),
                     paste0('BDC', 1:2),
                     paste0('CO',  1:2)),
  expectancies   = c(rep(1,   5), rep(0.5, 2), rep(0.5, 2)),
  ondrug         = list(pathA = c(rep(1, 5), rep(1, 2), 1, 0))
)

## --------------------------------------------------------------- ##
## 2. Modified lme_analysis that returns the Dbc main-effect row   ##
##    rather than the bm:Dbc interaction row.                      ##
## --------------------------------------------------------------- ##

lme_main_effect <- function(trialdesign_set, dat) {
  ## reuse pmsimstats-ng setup but extract Dbc rather than bm:Dbc
  res <- lme_analysis(trialdesign_set, dat,
                      op = list(useDE          = TRUE,
                                t_random_slope = FALSE,
                                full_model_out = TRUE,
                                carryover_t1half = 0,
                                simplecarryover  = FALSE,
                                carryover_scalefactor = 1))
  if (is.null(res$fit)) {
    return(data.table(beta_main = NA_real_, betaSE_main = NA_real_,
                      p_main = NA_real_, converged = FALSE))
  }
  tt <- summary(res$fit)$tTable
  if (!('Dbc' %in% rownames(tt))) {
    return(data.table(beta_main = NA_real_, betaSE_main = NA_real_,
                      p_main = NA_real_, converged = FALSE))
  }
  data.table(
    beta_main   = tt['Dbc', 'Value'],
    betaSE_main = tt['Dbc', 'Std.Error'],
    p_main      = tt['Dbc', 'p-value'],
    converged   = TRUE
  )
}

## --------------------------------------------------------------- ##
## 3. Single-rep functions                                         ##
## --------------------------------------------------------------- ##

run_nof1_rep <- function(true_effect, N = 35) {
  ## true_effect maps to br$max: br component grows toward max while
  ## on drug. Sx = BL - sum(components), so a positive br$max yields
  ## a negative Dbc coefficient. We pass abs(true_effect) and the
  ## sign comes out via the analysis.
  br_max <- abs(true_effect)
  mp <- data.table(
    N = N, c.bm = 0,
    carryover_t1half = 0,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.1, c.cfct = 0.05
  )
  rp <- data.table(
    cat  = c('tv', 'pb', 'br'),
    max  = c(2.0,  1.5,  br_max),
    disp = c(5,    5,    5),
    rate = c(0.42, 0.35, 0.42),
    sd   = c(2,    1,    2)
  )
  bp <- data.table(
    cat = c('bm', 'BL'),
    m   = c(0,    15),
    sd  = c(1,    5)
  )
  dat <- generateData(mp, rp, bp, td_hybrid$trialpaths[[1]],
                      empirical = FALSE,
                      makePositiveDefinite = TRUE,
                      dgp_architecture = 'mvn')
  dat[, path := 1]
  out <- tryCatch(lme_main_effect(td_hybrid$trialpaths, dat),
                  error = function(e) data.table(
                    beta_main = NA_real_, betaSE_main = NA_real_,
                    p_main = NA_real_, converged = FALSE))
  out
}

run_rct_rep <- function(true_effect, N_total = 70,
                        residual_sd = 5, baseline_mean = 15,
                        baseline_sd = 5) {
  ## Simple parallel-group RCT, 1:1 randomisation. Endpoint = wk-12
  ## measurement. We model Sx_endpoint ~ baseline + treat*group + eps.
  ## At matched per-participant N=35 in N-of-1, we double for RCT to
  ## keep total measurements roughly comparable would inflate; the
  ## brief specifies matched per-participant, so N_total = 70 RCT pts
  ## is NOT what we want -- we want N=35 per design (per-pt match).
  ## The user's request says 'N=35 per design (matched per-participant)'.
  ## Interpretation: each design uses 35 participants total. Override
  ## default and split 17/18.
  n_treat <- floor(N_total / 2)
  n_ctrl  <- N_total - n_treat
  bl <- rnorm(N_total, baseline_mean, baseline_sd)
  group <- c(rep(1, n_treat), rep(0, n_ctrl))
  ## true_effect is the mean change in nightmares/week for treated.
  ## Negative effect = improvement.
  y <- bl + true_effect * group + rnorm(N_total, 0, residual_sd)
  fit <- tryCatch(lm(y ~ bl + group),
                  error = function(e) NULL)
  if (is.null(fit)) {
    return(data.table(beta_main = NA_real_, betaSE_main = NA_real_,
                      p_main = NA_real_, converged = FALSE))
  }
  s <- summary(fit)$coefficients
  if (!('group' %in% rownames(s))) {
    return(data.table(beta_main = NA_real_, betaSE_main = NA_real_,
                      p_main = NA_real_, converged = FALSE))
  }
  data.table(
    beta_main   = s['group', 'Estimate'],
    betaSE_main = s['group', 'Std. Error'],
    p_main      = s['group', 'Pr(>|t|)'],
    converged   = TRUE
  )
}

## --------------------------------------------------------------- ##
## 4. Calibration: confirm one rep of each works                   ##
## --------------------------------------------------------------- ##

t0 <- Sys.time()
cat('Calibration runs:\n')
calib_n <- run_nof1_rep(-2.0, N = 35)
cat('  N-of-1 (effect=-2): beta=', round(calib_n$beta_main, 3),
    ' p=', round(calib_n$p_main, 4), '\n', sep = '')
calib_r <- run_rct_rep(-2.0, N_total = 35)
cat('  RCT    (effect=-2): beta=', round(calib_r$beta_main, 3),
    ' p=', round(calib_r$p_main, 4), '\n', sep = '')
cat('Calibration time:',
    round(as.numeric(Sys.time() - t0, units = 'secs'), 2), 's\n')

## --------------------------------------------------------------- ##
## 5. Time a small batch (10 reps each cell) to set rep budget     ##
## --------------------------------------------------------------- ##

cells <- list(
  list(design = 'Nof1', true_effect = 0,    fn = run_nof1_rep,
       args = list(N = 35)),
  list(design = 'Nof1', true_effect = -2.0, fn = run_nof1_rep,
       args = list(N = 35)),
  list(design = 'RCT',  true_effect = 0,    fn = run_rct_rep,
       args = list(N_total = 35)),
  list(design = 'RCT',  true_effect = -2.0, fn = run_rct_rep,
       args = list(N_total = 35))
)

t1 <- Sys.time()
cat('Timing 10 reps per cell to project budget...\n')
for (cell in cells) {
  for (i in 1:10) {
    do.call(cell$fn, c(list(true_effect = cell$true_effect),
                       cell$args))
  }
}
elapsed_pilot <- as.numeric(Sys.time() - t1, units = 'secs')
cat('Pilot (10 reps x 4 cells) took ', round(elapsed_pilot, 2),
    's\n', sep = '')

## Decide reps. Target: complete in ~3.5 min after calibration.
## Estimated cost per rep across 4 cells = elapsed_pilot / 10.
## Budget: 210s. n_reps = floor(210 / (elapsed_pilot / 10)).
budget_remaining_s <- 210
reps_target <- floor(budget_remaining_s / (elapsed_pilot / 10))
n_reps <- max(50, min(200, reps_target))
cat('Reps per cell:', n_reps, '\n')

## --------------------------------------------------------------- ##
## 6. Main loop                                                    ##
## --------------------------------------------------------------- ##

t2 <- Sys.time()
all_results <- vector('list', length(cells) * n_reps)
idx <- 1L
for (cell in cells) {
  for (i in seq_len(n_reps)) {
    r <- tryCatch(
      do.call(cell$fn, c(list(true_effect = cell$true_effect),
                         cell$args)),
      error = function(e) data.table(
        beta_main = NA_real_, betaSE_main = NA_real_,
        p_main = NA_real_, converged = FALSE))
    r[, design      := cell$design]
    r[, true_effect := cell$true_effect]
    r[, rep_idx     := i]
    all_results[[idx]] <- r
    idx <- idx + 1L
  }
  cat('  finished', cell$design, 'effect =', cell$true_effect,
      '\n')
}
elapsed_main <- as.numeric(Sys.time() - t2, units = 'secs')
cat('Main loop:', round(elapsed_main, 2), 's\n')

reps <- rbindlist(all_results, fill = TRUE)
setcolorder(reps, c('design', 'true_effect', 'rep_idx',
                    'beta_main', 'betaSE_main', 'p_main',
                    'converged'))
saveRDS(reps, file.path(out_dir_data, '04-mainfx-replicates.rds'))
cat('Saved replicates:', nrow(reps), 'rows\n')

## --------------------------------------------------------------- ##
## 7. Morris-style summary                                          ##
## --------------------------------------------------------------- ##

alpha <- 0.05
summary_dt <- reps[, {
  conv <- converged & !is.na(p_main)
  n_used <- sum(conv)
  rej   <- sum(p_main[conv] < alpha)
  rate  <- if (n_used > 0) rej / n_used else NA_real_
  ## binomial MCSE for rate
  mcse_rate <- if (n_used > 0)
    sqrt(rate * (1 - rate) / n_used) else NA_real_
  mb <- mean(beta_main[conv])
  sb <- sd(beta_main[conv])
  mcse_mb <- if (n_used > 1) sb / sqrt(n_used) else NA_real_
  list(
    n_reps           = .N,
    n_converged      = n_used,
    conv_rate        = n_used / .N,
    power_or_typeI   = rate,
    mcse_rate        = mcse_rate,
    mean_beta        = mb,
    sd_beta          = sb,
    mcse_mean_beta   = mcse_mb
  )
}, by = .(design, true_effect)]

setorder(summary_dt, design, true_effect)
fwrite(summary_dt,
       file.path(out_dir_data, '04-mainfx-summary.txt'),
       sep = '\t')
cat('\nSummary table:\n')
print(summary_dt)

## --------------------------------------------------------------- ##
## 8. Figure                                                       ##
## --------------------------------------------------------------- ##

plot_dt <- copy(summary_dt)
plot_dt[, scenario := ifelse(true_effect == 0,
                             'Type I error (effect = 0)',
                             'Power (effect = -2)')]
plot_dt[, ci_lo := pmax(0, power_or_typeI -
                          1.96 * mcse_rate)]
plot_dt[, ci_hi := pmin(1, power_or_typeI +
                          1.96 * mcse_rate)]

p <- ggplot(plot_dt,
            aes(x = design, y = power_or_typeI,
                fill = design)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.2) +
  geom_hline(data = data.frame(
    scenario = 'Type I error (effect = 0)', y = 0.05),
    aes(yintercept = y), linetype = 'dashed') +
  facet_wrap(~ scenario) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = c('Nof1' = 'steelblue',
                               'RCT'  = 'darkorange'),
                    guide = 'none') +
  labs(x = 'Design',
       y = 'Rejection rate',
       title = 'Paper 04 prototype: average treatment effect',
       subtitle = paste0('N=35 per design, ',
                         summary_dt$n_reps[1],
                         ' reps per cell')) +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir_fig, '04-mainfx-power.pdf'),
       p, width = 7, height = 4)

cat('\nWrote:\n  ', file.path(out_dir_data,
    '04-mainfx-replicates.rds'),
    '\n  ', file.path(out_dir_data, '04-mainfx-summary.txt'),
    '\n  ', file.path(out_dir_fig, '04-mainfx-power.pdf'), '\n',
    sep = '')

cat('Total elapsed:',
    round(as.numeric(Sys.time() - t0, units = 'secs'), 2),
    's\n')
