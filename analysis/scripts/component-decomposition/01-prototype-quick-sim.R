## 06-component-decomposition: prototype quick simulation (10-min budget)
## Purpose: Morris-compatible prototype for paper 06 (Phase 1, Study A).
## Design: 3 pb_strength x 3 N x 2 analyses = 18 cells at OL+BDC,
##   c.bm = 0.45, t1/2 = 1.
## Target: 500 reps per cell, n_fits = 500 * 18 = 9000.
## Parallelism: HARD REQUIREMENT mclapply with mc.cores = 8.
## Time budget: 540 s with graceful abort, progress every 100 reps.
## See analysis/scripts/component-decomposition/00-ademp-pre-registration.md
## and 02-deviations-log.md.

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(nlme)
  library(ggplot2)
  library(parallel)
})

set.seed(20260507)

out_data_dir <- 'analysis/data/quick-sim'
out_fig_dir  <- 'analysis/figures/quick-sim'
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir,  showWarnings = FALSE, recursive = TRUE)

TIME_BUDGET_SEC <- 540
PROGRESS_EVERY  <- 100L
TARGET_REPS     <- 500L
N_CORES         <- 8L
SCRIPT_T0       <- Sys.time()

cat('Cores requested:', N_CORES,
    '  detectCores =', parallel::detectCores(), '\n')

elapsed_sec <- function() as.numeric(difftime(Sys.time(), SCRIPT_T0,
                                              units = 'secs'))

## ---- Trial design: Hybrid OL+BDC -----------------------------------------
## Single-path 8 visit OL+BDC; on-drug indicator 1,1,1,1,1, 1,0,0.
td <- buildtrialdesign(
  name_longform  = 'OL+BDC',
  name_shortform = 'OLBDC',
  timepoints   = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20),
  timeptnames  = c(paste0('OL', 1:5), paste0('BDC', 1:3)),
  expectancies = c(rep(1, 5), rep(0.5, 3)),
  ondrug = list(pathA = c(rep(1, 5), 1, 0, 0))
)

## ---- Canonical parameter setup -------------------------------------------
## carryover_t1half = 1 (per task spec).
make_params <- function(N = 70, c.bm = 0.45, pb_max = 6.50647) {
  mp <- data.table(
    N = N, c.bm = c.bm,
    carryover_t1half = 1,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.1, c.cfct = 0.05
  )
  rp <- data.table(
    cat  = c('tv', 'pb', 'br'),
    max  = c(10.98604, pb_max, 10.98604),
    disp = c(5, 5, 5),
    rate = c(0.42, 0.35, 0.42),
    sd   = c(5, 2, 5)
  )
  bp <- data.table(
    cat = c('bm', 'BL'),
    m   = c(0, 70),
    sd  = c(1, 10)
  )
  list(mp = mp, rp = rp, bp = bp)
}

## ---- True coefficient (Architecture B) -----------------------------------
## With c.bm = 0.45, sigma_br = 5, sigma_bm = 1, the population conditional
## moderation of BR by bm has slope 2.25; on Sx = BL - (BR + PB + TV) the
## bm:Dbc coefficient is the negative.
TRUE_BETA <- -0.45 * 5 / 1

## ---- Phase-augmented analysis --------------------------------------------
## At N >= 70 attempt the FULL pre-registered formula including Dbc:phase.
## If lme fails because of rank deficiency, fall back to the reduced
## formula. At N = 35 go directly to the reduced formula. Record the
## dropped term (or 'none').
phase_augmented_fit <- function(td, dat, N) {
  op_prep <- list(useDE = FALSE, t_random_slope = FALSE,
                  full_model_out = TRUE, simplecarryover = FALSE,
                  carryover_t1half = 1, carryover_scalefactor = 1)
  prep <- tryCatch(lme_analysis(td$trialpaths, dat, op_prep),
                   error = function(e) NULL)
  if (is.null(prep) || is.null(prep$datamerged)) {
    return(data.table(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
                      converged = FALSE,
                      formula_dropped = 'prep_failed'))
  }
  d <- copy(prep$datamerged)
  d <- d[!is.na(Sx) & !is.na(bm) & !is.na(t) & !is.na(Dbc) & !is.na(De)]
  d <- d[t > 0]
  d[, phase := factor(ifelse(De >= 0.99, 'OL', 'BDC'),
                      levels = c('OL', 'BDC'))]
  if (length(unique(d$phase)) < 2) {
    return(data.table(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
                      converged = FALSE,
                      formula_dropped = 'single_phase'))
  }

  ctl <- nlme::lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)
  formula_dropped <- 'none'
  fit <- NULL

  if (N >= 70) {
    fit <- tryCatch(
      nlme::lme(Sx ~ bm + t + Dbc + phase + Dbc:phase
                     + bm:Dbc + bm:Dbc:phase,
                random = ~1 | ptID,
                correlation = nlme::corCAR1(form = ~t | ptID),
                data = d, control = ctl),
      error = function(e) NULL
    )
    if (is.null(fit)) formula_dropped <- 'Dbc:phase'
  } else {
    formula_dropped <- 'Dbc:phase'
  }

  if (is.null(fit)) {
    fit <- tryCatch(
      nlme::lme(Sx ~ bm + t + Dbc + phase + bm:Dbc + bm:Dbc:phase,
                random = ~1 | ptID,
                correlation = nlme::corCAR1(form = ~t | ptID),
                data = d, control = ctl),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) {
    return(data.table(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
                      converged = FALSE,
                      formula_dropped = paste0(formula_dropped, '+all')))
  }
  ct <- summary(fit)$tTable
  target <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(ct))
  if (length(target) == 0) {
    return(data.table(beta = NA_real_, betaSE = NA_real_, p = NA_real_,
                      converged = FALSE,
                      formula_dropped = paste0(formula_dropped, '+bmDbc')))
  }
  data.table(beta            = unname(ct[target[1], 'Value']),
             betaSE          = unname(ct[target[1], 'Std.Error']),
             p               = unname(ct[target[1], 'p-value']),
             converged       = TRUE,
             formula_dropped = formula_dropped)
}

## ---- One-component analysis wrapper --------------------------------------
one_component_fit <- function(td, dat) {
  op <- list(useDE = FALSE, t_random_slope = FALSE,
             full_model_out = FALSE, simplecarryover = FALSE,
             carryover_t1half = 1, carryover_scalefactor = 1)
  res <- tryCatch(lme_analysis(td$trialpaths, dat, op),
                  error = function(e) NULL)
  if (is.null(res)) {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = 'fit_failed'))
  }
  data.table(beta            = res$beta,
             betaSE          = res$betaSE,
             p               = res$p,
             converged       = !is.na(res$beta),
             formula_dropped = 'none')
}

## ---- One replicate -------------------------------------------------------
one_rep <- function(rep_idx, pb_strength_label, pb_max, N_val) {
  set.seed(20260507L + rep_idx + 1000L * N_val
           + as.integer(round(100 * as.numeric(pb_strength_label))))
  params <- make_params(N = N_val, c.bm = 0.45, pb_max = pb_max)
  dat <- tryCatch(
    generateData(params$mp, params$rp, params$bp, td$trialpaths[[1]],
                 empirical = FALSE, makePositiveDefinite = TRUE),
    error = function(e) NULL
  )
  if (is.null(dat)) {
    na_row <- data.table(beta = NA_real_, betaSE = NA_real_,
                         p = NA_real_, converged = FALSE,
                         formula_dropped = 'gen_failed')
    return(rbind(
      cbind(pb_strength = pb_strength_label, N = N_val,
            analysis = 'one_component', rep_idx = rep_idx, na_row),
      cbind(pb_strength = pb_strength_label, N = N_val,
            analysis = 'phase_augmented', rep_idx = rep_idx, na_row)
    ))
  }
  dat[, path := 1]
  oc <- one_component_fit(td, dat)
  pa <- phase_augmented_fit(td, dat, N = N_val)
  rbind(
    cbind(pb_strength = pb_strength_label, N = N_val,
          analysis = 'one_component', rep_idx = rep_idx, oc),
    cbind(pb_strength = pb_strength_label, N = N_val,
          analysis = 'phase_augmented', rep_idx = rep_idx, pa)
  )
}

## ---- Cell grid -----------------------------------------------------------
pb_levels <- c(0, 6.50647, 13.01294)
pb_labels <- c('0', '6.5', '13')
N_levels  <- c(35L, 70L, 100L)

cells <- CJ(pb_idx = seq_along(pb_levels),
            N      = N_levels,
            sorted = FALSE)
cells[, pb_label := pb_labels[pb_idx]]
cells[, pb_max   := pb_levels[pb_idx]]
n_cells <- nrow(cells)

cat(sprintf('Cells (pb_strength x N): %d (analyses x cells = %d).\n',
            n_cells, n_cells * 2L))
cat(sprintf('Target reps/cell: %d; total fits = %d.\n',
            TARGET_REPS, TARGET_REPS * n_cells * 2L))
cat(sprintf('Wall-time budget: %.0f s; mc.cores = %d.\n\n',
            TIME_BUDGET_SEC, N_CORES))

## ---- Pilot timing (tiny parallel batch) ----------------------------------
pilot_args <- list()
for (cell_i in seq_len(nrow(cells))) {
  pilot_args[[length(pilot_args) + 1L]] <- list(
    rep_idx = 0L,
    pb_strength_label = cells$pb_label[cell_i],
    pb_max = cells$pb_max[cell_i],
    N_val = cells$N[cell_i]
  )
}
t_pilot <- Sys.time()
pilot_out <- parallel::mclapply(
  pilot_args,
  function(a) one_rep(a$rep_idx, a$pb_strength_label,
                      a$pb_max, a$N_val),
  mc.cores = N_CORES
)
elapsed_pilot <- as.numeric(difftime(Sys.time(), t_pilot, units = 'secs'))
cat(sprintf('[%6.1fs] Pilot: %d cells in %.2f s wall (%.0f ms/cell).\n',
            elapsed_sec(), length(pilot_args), elapsed_pilot,
            1000 * elapsed_pilot / length(pilot_args)))

## ---- Chunked round-robin main loop with graceful abort -------------------
## Each chunk dispatches CHUNK args across all cells (round-robin).
CHUNK_PER_CELL <- 25L

results_acc <- vector('list', 0L)
reps_done <- setNames(rep(0L, n_cells),
                      paste0('cell_', seq_len(n_cells)))
total_reps_done <- 0L
last_progress <- 0L
total_target <- TARGET_REPS * n_cells
abort_flag <- FALSE

while (any(reps_done < TARGET_REPS) && !abort_flag) {

  if (elapsed_sec() > TIME_BUDGET_SEC) {
    cat(sprintf('[%6.1fs] Wall-time hit before chunk; aborting.\n',
                elapsed_sec()))
    abort_flag <- TRUE
    break
  }

  chunk_args <- list()
  for (cell_i in seq_len(nrow(cells))) {
    key <- paste0('cell_', cell_i)
    remaining <- TARGET_REPS - reps_done[[key]]
    if (remaining <= 0L) next
    take <- min(CHUNK_PER_CELL, remaining)
    next_start <- reps_done[[key]] + 1L
    for (i in seq_len(take)) {
      chunk_args[[length(chunk_args) + 1L]] <- list(
        rep_idx = next_start + i - 1L,
        pb_strength_label = cells$pb_label[cell_i],
        pb_max = cells$pb_max[cell_i],
        N_val = cells$N[cell_i],
        cell_key = key
      )
    }
  }
  if (length(chunk_args) == 0L) break

  ## Estimate remaining time. If we cannot complete the chunk, scale down.
  if (length(results_acc) > 0L && elapsed_pilot > 0) {
    remaining_budget <- TIME_BUDGET_SEC - elapsed_sec() - 5
    ## crude per-rep estimate from elapsed and reps so far
    if (total_reps_done > 0L) {
      sec_per_rep <- elapsed_sec() / total_reps_done
      max_reps_remaining <- floor(remaining_budget / sec_per_rep)
      if (max_reps_remaining < length(chunk_args)) {
        if (max_reps_remaining <= 0L) {
          cat(sprintf('[%6.1fs] No budget remaining; aborting.\n',
                      elapsed_sec()))
          abort_flag <- TRUE
          break
        }
        chunk_args <- chunk_args[seq_len(max_reps_remaining)]
        cat(sprintf('[%6.1fs] Trimming chunk to %d reps (budget).\n',
                    elapsed_sec(), max_reps_remaining))
      }
    }
  }

  chunk_out <- parallel::mclapply(
    chunk_args,
    function(a) one_rep(a$rep_idx, a$pb_strength_label,
                        a$pb_max, a$N_val),
    mc.cores = N_CORES
  )

  good <- !vapply(chunk_out,
                  function(x) inherits(x, 'try-error') || is.null(x),
                  logical(1))
  chunk_dt <- rbindlist(chunk_out[good], fill = TRUE)
  results_acc[[length(results_acc) + 1L]] <- chunk_dt

  for (a in chunk_args) {
    reps_done[[a$cell_key]] <- reps_done[[a$cell_key]] + 1L
  }
  total_reps_done <- sum(reps_done)

  if (total_reps_done - last_progress >= PROGRESS_EVERY) {
    last_progress <- total_reps_done
    cat(sprintf('[%6.1fs] Progress: %d/%d reps; per-cell min/max = %d/%d\n',
                elapsed_sec(), total_reps_done, total_target,
                min(reps_done), max(reps_done)))
  }

  if (elapsed_sec() > TIME_BUDGET_SEC) {
    cat(sprintf('[%6.1fs] Wall-time hit after chunk; aborting.\n',
                elapsed_sec()))
    abort_flag <- TRUE
  }
}

cat(sprintf('[%6.1fs] Loop ended; total replicates produced = %d.\n',
            elapsed_sec(), total_reps_done))
cat('Per-cell rep counts:\n')
print(data.table(cell = seq_len(n_cells),
                 pb = cells$pb_label, N = cells$N,
                 reps_done = unname(reps_done)))

## ---- Persist per-replicate data ------------------------------------------
all_rows <- rbindlist(results_acc, fill = TRUE)
setnames(all_rows, c('beta', 'betaSE', 'p'),
         c('beta_bmDbc', 'betaSE_bmDbc', 'p_bmDbc'))
setcolorder(all_rows, c('pb_strength', 'N', 'analysis', 'rep_idx',
                        'beta_bmDbc', 'betaSE_bmDbc', 'p_bmDbc',
                        'converged', 'formula_dropped'))
saveRDS(all_rows, file.path(out_data_dir, '06-decomp-replicates.rds'))
cat(sprintf('Wrote %s (n = %d rows).\n',
            file.path(out_data_dir, '06-decomp-replicates.rds'),
            nrow(all_rows)))

## ---- Morris summary table ------------------------------------------------
## Performance measures with MCSE per Morris, White & Crowther (2019):
##   bias  = mean_beta - TRUE_BETA;  MCSE_bias = sd_beta / sqrt(n)
##   power = mean(p < 0.05);         MCSE_power = sqrt(p(1-p)/n)
summarise_cell <- function(d) {
  conv <- d[converged == TRUE]
  n <- nrow(conv)
  if (n == 0) {
    return(data.table(
      n_reps_done    = nrow(d),
      n_converged    = 0L,
      conv_rate      = 0,
      power          = NA_real_,
      mcse_power     = NA_real_,
      mean_beta      = NA_real_,
      sd_beta        = NA_real_,
      mcse_mean_beta = NA_real_,
      bias           = NA_real_,
      mcse_bias      = NA_real_,
      n_fmla_full    = 0L,
      n_fmla_dropped = 0L,
      most_common_drop = NA_character_
    ))
  }
  power <- mean(conv$p_bmDbc < 0.05, na.rm = TRUE)
  mcse_power <- sqrt(power * (1 - power) / n)
  mean_beta <- mean(conv$beta_bmDbc, na.rm = TRUE)
  sd_beta   <- sd(conv$beta_bmDbc,  na.rm = TRUE)
  mcse_mean_beta <- sd_beta / sqrt(n)
  bias <- mean_beta - TRUE_BETA
  mcse_bias <- mcse_mean_beta
  drops <- d$formula_dropped
  drops_non_none <- drops[drops != 'none' & !is.na(drops)]
  drop_tab <- if (length(drops_non_none) == 0) NULL else
    sort(table(drops_non_none), decreasing = TRUE)
  data.table(
    n_reps_done    = nrow(d),
    n_converged    = n,
    conv_rate      = mean(d$converged, na.rm = TRUE),
    power          = power,
    mcse_power     = mcse_power,
    mean_beta      = mean_beta,
    sd_beta        = sd_beta,
    mcse_mean_beta = mcse_mean_beta,
    bias           = bias,
    mcse_bias      = mcse_bias,
    n_fmla_full    = sum(drops == 'none', na.rm = TRUE),
    n_fmla_dropped = length(drops_non_none),
    most_common_drop = if (is.null(drop_tab)) 'none' else
      names(drop_tab)[1]
  )
}

summary_tbl <- all_rows[,
  summarise_cell(.SD),
  by = .(pb_strength, N, analysis)]

setorder(summary_tbl, analysis, N, pb_strength)

write.table(summary_tbl,
            file = file.path(out_data_dir, '06-decomp-summary.txt'),
            sep = '\t', quote = FALSE, row.names = FALSE)
cat(sprintf('Wrote %s.\n',
            file.path(out_data_dir, '06-decomp-summary.txt')))
cat(sprintf('TRUE_BETA (Architecture B; -c.bm * br_sd / bm_sd): %.4f\n',
            TRUE_BETA))
print(summary_tbl)

## ---- Figure: bias vs PB strength faceted by analysis x N -----------------
plot_dat <- copy(summary_tbl)
plot_dat[, pb_strength_num := as.numeric(pb_strength)]
plot_dat[, analysis := factor(analysis,
                              levels = c('one_component', 'phase_augmented'),
                              labels = c('One-component',
                                         'Phase-augmented'))]
plot_dat[, N_lab := factor(paste0('N = ', N),
                           levels = paste0('N = ', sort(unique(N))))]

p_bias <- ggplot(plot_dat,
                 aes(x = pb_strength_num, y = bias)) +
  geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey40') +
  geom_line(linewidth = 0.7, colour = '#2c7fb8') +
  geom_point(size = 2.4, colour = '#2c7fb8') +
  geom_errorbar(aes(ymin = bias - 1.96 * mcse_bias,
                    ymax = bias + 1.96 * mcse_bias),
                width = 0.4, linewidth = 0.5, colour = '#2c7fb8') +
  facet_grid(N_lab ~ analysis) +
  scale_x_continuous(breaks = as.numeric(pb_labels)) +
  labs(x = 'Population PB component strength (m_PB)',
       y = expression('Bias of '*hat(beta)[bm:Dbc]*
                      '  (mean - true; truth = -2.25)'),
       title = 'Paper 06 prototype: bias of bm:Dbc by analysis and N',
       subtitle = sprintf(
         'OL+BDC, c.bm = 0.45, t1/2 = 1; target n_reps = %d/cell; truth = %.2f',
         TARGET_REPS, TRUE_BETA)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = 'grey92'))

ggsave(file.path(out_fig_dir, '06-decomp-bias.pdf'),
       p_bias, width = 7.2, height = 6.0)
cat(sprintf('Wrote %s.\n',
            file.path(out_fig_dir, '06-decomp-bias.pdf')))

cat(sprintf('Total elapsed: %.1f s. abort_flag = %s. Done.\n',
            elapsed_sec(), abort_flag))
