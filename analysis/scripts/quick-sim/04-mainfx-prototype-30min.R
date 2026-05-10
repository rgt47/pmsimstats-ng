## Paper 04 prototype (30-min budget): aggregated N-of-1 vs
## parallel-group RCT power curves for the average treatment effect.
##
## Cells: design x true_effect = 2 x 4 = 8 cells.
##   design      in {Nof1_hybrid, parallel_RCT}
##   true_effect in {0, -0.5, -1.0, -2.0}
## Target: 1500 reps per cell (12000 total). Wall-time cap 1700 s.
## Per-replicate output written to:
##   analysis/data/quick-sim/04-mainfx-replicates.rds

suppressPackageStartupMessages({
  library(devtools)
  library(data.table)
  library(ggplot2)
  library(nlme)
  library(parallel)
})

devtools::load_all('.', quiet = TRUE)

set.seed(20260507)

out_dir_data <- 'analysis/data/quick-sim'
out_dir_fig  <- 'analysis/figures/quick-sim'
dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_fig,  recursive = TRUE, showWarnings = FALSE)

WALLTIME_S <- 1700      # 28.3 min hard ceiling
TARGET_REPS <- 1500
N_PER_DESIGN <- 35

n_cores <- max(1L, parallel::detectCores() - 2L)
cat('Cores available:', parallel::detectCores(),
    ' workers used:', n_cores, '\n')

## --------------------------------------------------------------- ##
## 1. Trial design                                                 ##
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
## 2. Helpers                                                      ##
## --------------------------------------------------------------- ##

lme_main_effect <- function(trialdesign_set, dat) {
  res <- lme_analysis(trialdesign_set, dat,
                      op = list(useDE          = TRUE,
                                t_random_slope = FALSE,
                                full_model_out = TRUE,
                                carryover_t1half = 0,
                                simplecarryover  = FALSE,
                                carryover_scalefactor = 1))
  if (is.null(res$fit)) {
    return(data.table(beta_main = NA_real_,
                      betaSE_main = NA_real_,
                      p_main = NA_real_, converged = FALSE))
  }
  tt <- summary(res$fit)$tTable
  if (!('Dbc' %in% rownames(tt))) {
    return(data.table(beta_main = NA_real_,
                      betaSE_main = NA_real_,
                      p_main = NA_real_, converged = FALSE))
  }
  data.table(
    beta_main   = tt['Dbc', 'Value'],
    betaSE_main = tt['Dbc', 'Std.Error'],
    p_main      = tt['Dbc', 'p-value'],
    converged   = TRUE
  )
}

run_nof1_rep <- function(true_effect, N = N_PER_DESIGN) {
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
                    beta_main = NA_real_,
                    betaSE_main = NA_real_,
                    p_main = NA_real_, converged = FALSE))
  out
}

run_rct_rep <- function(true_effect, N_total = N_PER_DESIGN,
                        residual_sd = 5, baseline_mean = 15,
                        baseline_sd = 5) {
  n_treat <- floor(N_total / 2)
  n_ctrl  <- N_total - n_treat
  bl <- rnorm(N_total, baseline_mean, baseline_sd)
  group <- c(rep(1, n_treat), rep(0, n_ctrl))
  y <- bl + true_effect * group + rnorm(N_total, 0, residual_sd)
  fit <- tryCatch(lm(y ~ bl + group),
                  error = function(e) NULL)
  if (is.null(fit)) {
    return(data.table(beta_main = NA_real_,
                      betaSE_main = NA_real_,
                      p_main = NA_real_, converged = FALSE))
  }
  s <- summary(fit)$coefficients
  if (!('group' %in% rownames(s))) {
    return(data.table(beta_main = NA_real_,
                      betaSE_main = NA_real_,
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
## 3. Calibration                                                  ##
## --------------------------------------------------------------- ##

t0 <- Sys.time()
cat('Calibration:\n')
calib_n <- run_nof1_rep(-2.0)
cat('  Nof1_hybrid (eff=-2): beta=',
    round(calib_n$beta_main, 3),
    ' p=', round(calib_n$p_main, 4), '\n', sep = '')
calib_r <- run_rct_rep(-2.0)
cat('  parallel_RCT (eff=-2): beta=',
    round(calib_r$beta_main, 3),
    ' p=', round(calib_r$p_main, 4), '\n', sep = '')

## --------------------------------------------------------------- ##
## 4. Build the work plan                                          ##
## --------------------------------------------------------------- ##

cells <- list(
  list(design = 'Nof1_hybrid',  true_effect =  0.0),
  list(design = 'Nof1_hybrid',  true_effect = -0.5),
  list(design = 'Nof1_hybrid',  true_effect = -1.0),
  list(design = 'Nof1_hybrid',  true_effect = -2.0),
  list(design = 'parallel_RCT', true_effect =  0.0),
  list(design = 'parallel_RCT', true_effect = -0.5),
  list(design = 'parallel_RCT', true_effect = -1.0),
  list(design = 'parallel_RCT', true_effect = -2.0)
)

run_one <- function(seed_offset, design, true_effect) {
  set.seed(20260507L + seed_offset)
  if (design == 'Nof1_hybrid') {
    out <- run_nof1_rep(true_effect)
  } else {
    out <- run_rct_rep(true_effect)
  }
  out[, design      := design]
  out[, true_effect := true_effect]
  out[, rep_idx     := seed_offset]
  out
}

## Time a small parallel batch (40 reps, 5/cell) to project cost.
pilot_args <- list()
for (cell in cells) {
  for (i in 1:5) {
    pilot_args[[length(pilot_args) + 1L]] <-
      list(seed_offset = i, design = cell$design,
           true_effect = cell$true_effect)
  }
}
t1 <- Sys.time()
pilot_out <- parallel::mclapply(
  pilot_args,
  function(a) run_one(a$seed_offset, a$design, a$true_effect),
  mc.cores = n_cores)
elapsed_pilot <- as.numeric(Sys.time() - t1, units = 'secs')
n_pilot <- length(pilot_args)
sec_per_rep_parallel <- elapsed_pilot / n_pilot
cat('Pilot (', n_pilot, ' reps in parallel) took ',
    round(elapsed_pilot, 2), ' s; ',
    round(sec_per_rep_parallel * 1000, 1),
    ' ms per rep wall-time.\n', sep = '')

## --------------------------------------------------------------- ##
## 5. Run in chunks of 200, checking the wall clock between chunks ##
## --------------------------------------------------------------- ##

CHUNK <- 200L              # also defines progress cadence

results_acc <- vector('list', 0L)
total_done <- 0L
aborted <- FALSE

t_main_start <- Sys.time()
elapsed_total <- function() {
  as.numeric(Sys.time() - t0, units = 'secs')
}

next_seed <- 1L

reps_per_cell_completed <- setNames(
  rep(0L, length(cells)),
  vapply(cells, function(c) paste(c$design, c$true_effect),
         character(1)))

while (any(reps_per_cell_completed < TARGET_REPS) && !aborted) {

  ## Build the next chunk: round-robin across cells whose
  ## counter is < TARGET_REPS, to keep cells balanced if we abort.
  chunk_args <- list()
  for (k in seq_along(cells)) {
    cell <- cells[[k]]
    key <- paste(cell$design, cell$true_effect)
    remaining <- TARGET_REPS - reps_per_cell_completed[[key]]
    if (remaining <= 0L) next
    take <- min(CHUNK %/% length(cells), remaining)
    if (take <= 0L) take <- 1L
    for (i in seq_len(take)) {
      chunk_args[[length(chunk_args) + 1L]] <- list(
        seed_offset = next_seed,
        design = cell$design,
        true_effect = cell$true_effect,
        cell_key = key
      )
      next_seed <- next_seed + 1L
    }
  }
  if (length(chunk_args) == 0L) break

  ## Pre-flight time check.
  if (elapsed_total() > WALLTIME_S) {
    aborted <- TRUE; break
  }

  chunk_out <- parallel::mclapply(
    chunk_args,
    function(a) run_one(a$seed_offset, a$design, a$true_effect),
    mc.cores = n_cores
  )

  ## Bind successful results, record per-cell counts.
  good <- !vapply(chunk_out,
                  function(x) inherits(x, 'try-error') ||
                              is.null(x),
                  logical(1))
  chunk_dt <- rbindlist(chunk_out[good], fill = TRUE)
  results_acc[[length(results_acc) + 1L]] <- chunk_dt

  for (a in chunk_args) {
    reps_per_cell_completed[[a$cell_key]] <-
      reps_per_cell_completed[[a$cell_key]] + 1L
  }

  total_done <- total_done + nrow(chunk_dt)

  cat('  progress: ', total_done, ' fits;  elapsed ',
      round(elapsed_total(), 1), ' s;  per-cell counts: ',
      paste(reps_per_cell_completed, collapse = '/'),
      '\n', sep = '')

  if (elapsed_total() > WALLTIME_S) {
    aborted <- TRUE
  }
}

reps <- rbindlist(results_acc, fill = TRUE)
elapsed_main <- as.numeric(Sys.time() - t_main_start, units = 'secs')
cat('Main loop:', round(elapsed_main, 1),
    's; total reps recorded:', nrow(reps), '\n')
if (aborted) {
  cat('NOTE: wall-time budget reached; writing partial output.\n')
}

setcolorder(reps, c('design', 'true_effect', 'rep_idx',
                    'beta_main', 'betaSE_main', 'p_main',
                    'converged'))
saveRDS(reps,
        file.path(out_dir_data, '04-mainfx-replicates.rds'))
cat('Wrote replicates: ', nrow(reps), ' rows to ',
    file.path(out_dir_data, '04-mainfx-replicates.rds'),
    '\n', sep = '')

## --------------------------------------------------------------- ##
## 6. Morris-style summary                                         ##
## --------------------------------------------------------------- ##

alpha <- 0.05
summary_dt <- reps[, {
  conv <- converged & !is.na(p_main)
  n_used <- sum(conv)
  rej   <- sum(p_main[conv] < alpha)
  rate  <- if (n_used > 0) rej / n_used else NA_real_
  mcse_rate <- if (n_used > 0)
    sqrt(rate * (1 - rate) / n_used) else NA_real_
  mb <- mean(beta_main[conv])
  sb <- sd(beta_main[conv])
  mcse_mb <- if (n_used > 1) sb / sqrt(n_used) else NA_real_
  ## Coverage of true_effect by 95% Wald CI from beta_main
  ## +/- 1.96 * betaSE_main.
  bs <- betaSE_main[conv]
  bm <- beta_main[conv]
  ci_lo <- bm - 1.96 * bs
  ci_hi <- bm + 1.96 * bs
  cov_in <- (true_effect[1] >= ci_lo) & (true_effect[1] <= ci_hi)
  coverage <- mean(cov_in, na.rm = TRUE)
  mcse_cov <- if (n_used > 0)
    sqrt(coverage * (1 - coverage) / n_used) else NA_real_
  list(
    n_reps         = .N,
    n_converged    = n_used,
    conv_rate      = n_used / .N,
    rejection_rate = rate,
    mcse_rate      = mcse_rate,
    mean_beta      = mb,
    sd_beta        = sb,
    mcse_mean_beta = mcse_mb,
    coverage       = coverage,
    mcse_coverage  = mcse_cov
  )
}, by = .(design, true_effect)]

setorder(summary_dt, design, true_effect)
fwrite(summary_dt,
       file.path(out_dir_data, '04-mainfx-summary.txt'),
       sep = '\t')

cat('\nSummary table:\n')
print(summary_dt)

## --------------------------------------------------------------- ##
## 7. Power-curve figure                                           ##
## --------------------------------------------------------------- ##

plot_dt <- copy(summary_dt)
plot_dt[, ci_lo := pmax(0, rejection_rate -
                          1.96 * mcse_rate)]
plot_dt[, ci_hi := pmin(1, rejection_rate +
                          1.96 * mcse_rate)]

p <- ggplot(plot_dt,
            aes(x = true_effect, y = rejection_rate,
                colour = design, group = design)) +
  geom_hline(yintercept = 0.05, linetype = 'dashed',
             colour = 'grey40') +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.05, alpha = 0.7) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  scale_x_reverse(breaks = c(0, -0.5, -1.0, -2.0)) +
  scale_colour_manual(values = c('Nof1_hybrid'  = 'steelblue',
                                 'parallel_RCT' = 'darkorange')) +
  labs(x = 'True average treatment effect (negative = improvement)',
       y = 'Rejection rate at alpha = 0.05',
       colour = 'Design',
       title = paste0('Paper 04 prototype: ',
                      'power curves for the main treatment effect'),
       subtitle = paste0('N = ', N_PER_DESIGN,
                         ' per design; up to ', TARGET_REPS,
                         ' reps per cell;',
                         ' dashed line = nominal alpha')) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom')

ggsave(file.path(out_dir_fig, '04-mainfx-power.pdf'),
       p, width = 7.5, height = 5)

cat('\nWrote:\n  ',
    file.path(out_dir_data, '04-mainfx-replicates.rds'),
    '\n  ',
    file.path(out_dir_data, '04-mainfx-summary.txt'),
    '\n  ',
    file.path(out_dir_fig, '04-mainfx-power.pdf'),
    '\n', sep = '')

cat('\nTotal elapsed:',
    round(as.numeric(Sys.time() - t0, units = 'secs'), 1),
    's; aborted_due_to_walltime =', aborted, '\n')
