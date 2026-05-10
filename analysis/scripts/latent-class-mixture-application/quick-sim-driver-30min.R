## Parallel quick-sim driver for paper 03-latent-class-mixture-
## application, 30-minute wall-time budget.
##
## Cells: gating_slope in {0, 1.0, 1.5} crossed with N in {35, 70}
## (6 cells), Hybrid OL+BDC design, 600 reps/cell. Two analyses
## per replicate (lme + lcmm) with three test channels recorded
## (lme bm:Dbc, lcmm gating Wald unconditional, lcmm gating
## BIC-conditional).
##
## Parallelism: mclapply with mc.cores = 8 primary, with fall-back
## escalation to mc.cores = 4, then sequential at n_reps = 200/cell
## if serialisation issues hit lcmm.
##
## Time budget: aborts gracefully and writes whatever is finished
## once elapsed > 1700 s. Progress printed every 100 reps.
##
## Run from the package root:
##   Rscript analysis/scripts/latent-class-mixture-application/\
##           quick-sim-driver-30min.R

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
  library(parallel)
})

source('analysis/scripts/latent-class-mixture-application/01-lcmm-wrapper.R')
source('analysis/scripts/latent-class-mixture-application/02-heterogeneous-slopes-dgp.R')

DATA_DIR <- file.path('analysis', 'data', 'quick-sim')
FIG_DIR  <- file.path('analysis', 'figures', 'quick-sim')
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

SEED_BASE <- 20260507L
GATING_SLOPES <- c(0, 1.0, 1.5)
N_VALUES <- c(35L, 70L)
N_REPS_PER_CELL <- 600L
MC_CORES_PRIMARY <- 8L
MC_CORES_FALLBACK <- 4L
SEQ_REPS <- 200L
WALL_BUDGET_SEC <- 1700        # graceful abort threshold
PROGRESS_EVERY  <- 100L

td <- buildtrialdesign(
  name_longform  = 'open label',
  name_shortform = 'OL',
  timepoints     = cumulative(rep(2.5, 8)),
  timeptnames    = paste0('OL', 1:8),
  expectancies   = rep(1, 8),
  ondrug         = list(pathA = rep(1, 8))
)

rp <- data.table(
  cat  = c('tv', 'pb', 'br'),
  max  = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd   = c(5, 2, 5)
)
bp <- data.table(
  cat = c('bm', 'BL'),
  m   = c(0, 70),
  sd  = c(1, 10)
)

op <- list(
  useDE = TRUE, t_random_slope = FALSE, full_model_out = FALSE,
  carryover_t1half = 1.0, simplecarryover = FALSE,
  carryover_scalefactor = 1
)

make_modelparam <- function(gating_slope, N_size) {
  if (gating_slope == 0) {
    data.table(
      N = N_size, c.bm = 0,
      carryover_t1half = 1.0,
      c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
      c.cf1t = 0.1, c.cfct = 0.05
    )
  } else {
    data.table(
      N = N_size, c.bm = 0,
      carryover_t1half = 1.0,
      c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
      c.cf1t = 0.1, c.cfct = 0.05,
      beta_R_factor = 1.0, beta_NR_factor = 0.0,
      gating_intercept = 0.0, gating_slope = gating_slope
    )
  }
}

simulate_dat <- function(gating_slope, N_size, rep_idx) {
  ## Seed strategy: combine slope cell index and N cell index so the
  ## same rep_idx differs across cells but is reproducible.
  slope_idx <- which(GATING_SLOPES == gating_slope)
  n_idx <- which(N_VALUES == N_size)
  set.seed(SEED_BASE + rep_idx +
           1000L * slope_idx +
           250000L * n_idx)
  mp <- make_modelparam(gating_slope, N_size)
  if (gating_slope == 0) {
    dat <- generateData(
      modelparam = mp, respparam = rp, blparam = bp,
      trialdesign = td$trialpaths[[1]],
      empirical = FALSE, makePositiveDefinite = TRUE,
      dgp_architecture = 'mvn'
    )
  } else {
    dat <- generate_random_slopes_data(
      modelparam = mp, respparam = rp, blparam = bp,
      trialdesign = td$trialpaths[[1]],
      empirical = FALSE, makePositiveDefinite = TRUE
    )
  }
  dat[, path := 1]
  dat
}

run_one <- function(gating_slope, N_size, rep_idx) {
  dat <- simulate_dat(gating_slope, N_size, rep_idx)

  t_lme0 <- Sys.time()
  out_lme <- tryCatch(
    lme_analysis(td$trialpaths, dat, op),
    error = function(e) data.table(
      beta = NA_real_, betaSE = NA_real_, p = NA_real_,
      issingular = NA, warning = paste0('err: ', conditionMessage(e))
    )
  )
  t_lme <- as.numeric(difftime(Sys.time(), t_lme0, units = 'secs'))

  t_lcmm0 <- Sys.time()
  out_lcmm <- tryCatch(
    lcmm_analysis(td$trialpaths, dat, op),
    error = function(e) data.table(
      beta = NA_real_, betaSE = NA_real_, p = NA_real_,
      beta_within = NA_real_, betaSE_within = NA_real_,
      p_within = NA_real_, delta_bic = NA_real_,
      lrt_stat = NA_real_, lrt_df = NA_integer_,
      bic_ng1 = NA_real_, bic_ng2 = NA_real_,
      issingular = NA, warning = paste0('err: ', conditionMessage(e)),
      entropy = NA_real_, bic = NA_real_, conv = 0L, conv_ng1 = 0L
    )
  )
  t_lcmm <- as.numeric(difftime(Sys.time(), t_lcmm0, units = 'secs'))

  list(
    lme = data.table(
      gating_slope = gating_slope, N = N_size,
      analysis = 'lme', rep_idx = rep_idx,
      beta = out_lme$beta, betaSE = out_lme$betaSE, p = out_lme$p,
      gating_p = NA_real_, within_p = NA_real_,
      delta_bic = NA_real_,
      conv = as.integer(!is.na(out_lme$p)), elapsed_s = t_lme
    ),
    lcmm = data.table(
      gating_slope = gating_slope, N = N_size,
      analysis = 'lcmm', rep_idx = rep_idx,
      beta = out_lcmm$beta, betaSE = out_lcmm$betaSE,
      p = out_lcmm$p,
      gating_p = out_lcmm$p, within_p = out_lcmm$p_within,
      delta_bic = out_lcmm$delta_bic,
      conv = as.integer(out_lcmm$conv == 1),
      elapsed_s = t_lcmm
    )
  )
}

worker <- function(job) {
  res <- tryCatch(
    run_one(job$gating_slope, job$N, job$rep_idx),
    error = function(e) {
      list(
        lme = data.table(
          gating_slope = job$gating_slope, N = job$N,
          analysis = 'lme', rep_idx = job$rep_idx,
          beta = NA_real_, betaSE = NA_real_, p = NA_real_,
          gating_p = NA_real_, within_p = NA_real_,
          delta_bic = NA_real_, conv = 0L, elapsed_s = NA_real_
        ),
        lcmm = data.table(
          gating_slope = job$gating_slope, N = job$N,
          analysis = 'lcmm', rep_idx = job$rep_idx,
          beta = NA_real_, betaSE = NA_real_, p = NA_real_,
          gating_p = NA_real_, within_p = NA_real_,
          delta_bic = NA_real_, conv = 0L, elapsed_s = NA_real_
        )
      )
    }
  )
  res
}

make_jobs <- function(gating_slopes, N_values, n_reps) {
  jobs <- list()
  for (g in gating_slopes) {
    for (n_size in N_values) {
      for (r in seq_len(n_reps)) {
        jobs[[length(jobs) + 1L]] <-
          list(gating_slope = g, N = n_size, rep_idx = r)
      }
    }
  }
  jobs
}

is_structural_failure <- function(x) {
  inherits(x, 'try-error') || is.null(x) ||
    !is.list(x) || is.null(x$lme) || is.null(x$lcmm)
}

## Time-bounded chunked executor: feed jobs to mclapply in chunks
## so that elapsed time can be checked between chunks. Each chunk
## takes roughly chunk_size / mc_cores * mean_per_rep seconds.
run_parallel_chunked <- function(jobs, mc_cores, t_global,
                                 wall_budget, chunk_size = 80L,
                                 progress_every = PROGRESS_EVERY) {
  results <- vector('list', length(jobs))
  total <- length(jobs)
  done <- 0L
  next_progress <- progress_every
  i <- 1L
  while (i <= total) {
    elapsed <- as.numeric(difftime(Sys.time(), t_global,
                                   units = 'secs'))
    if (elapsed > wall_budget) {
      cat(sprintf(
        '  TIME ABORT: elapsed %.1f s > %.1f s budget; %d/%d jobs complete\n',
        elapsed, wall_budget, done, total))
      break
    }
    j_end <- min(i + chunk_size - 1L, total)
    out_chunk <- mclapply(jobs[i:j_end], worker,
                          mc.cores = mc_cores,
                          mc.preschedule = FALSE,
                          mc.set.seed = FALSE)
    fail_idx <- which(vapply(out_chunk, is_structural_failure,
                             logical(1)))
    if (length(fail_idx) > 0) {
      cat(sprintf(
        '  STRUCTURAL FAIL: %d/%d in chunk starting at %d\n',
        length(fail_idx), length(out_chunk), i))
      return(list(results = NULL, done = done))
    }
    for (k in seq_along(out_chunk)) {
      results[[i + k - 1L]] <- out_chunk[[k]]
    }
    done <- done + length(out_chunk)
    if (done >= next_progress) {
      pct <- 100 * done / total
      el <- as.numeric(difftime(Sys.time(), t_global,
                                units = 'secs'))
      eta <- if (done > 0) el * (total - done) / done else NA_real_
      cat(sprintf(
        '  progress: %d/%d (%.1f%%) elapsed=%.1f s eta=%.1f s\n',
        done, total, pct, el, eta))
      next_progress <- next_progress + progress_every
    }
    i <- j_end + 1L
  }
  list(results = results[seq_len(done)], done = done)
}

run_sequential_capped <- function(jobs, t_global, wall_budget,
                                  progress_every = PROGRESS_EVERY) {
  results <- vector('list', length(jobs))
  done <- 0L
  next_progress <- progress_every
  for (i in seq_along(jobs)) {
    elapsed <- as.numeric(difftime(Sys.time(), t_global,
                                   units = 'secs'))
    if (elapsed > wall_budget) {
      cat(sprintf(
        '  TIME ABORT (sequential): elapsed %.1f s; %d/%d jobs done\n',
        elapsed, done, length(jobs)))
      break
    }
    results[[i]] <- worker(jobs[[i]])
    done <- done + 1L
    if (done >= next_progress) {
      el <- as.numeric(difftime(Sys.time(), t_global,
                                units = 'secs'))
      cat(sprintf('  progress: %d/%d elapsed=%.1f s\n',
                  done, length(jobs), el))
      next_progress <- next_progress + progress_every
    }
  }
  list(results = results[seq_len(done)], done = done)
}

t_global <- Sys.time()
cat(sprintf('30-min driver started at %s\n', format(t_global)))
cat(sprintf('Cells: gating_slope = %s; N = %s\n',
            paste(GATING_SLOPES, collapse = ', '),
            paste(N_VALUES, collapse = ', ')))
cat(sprintf('Target reps per cell: %d (total jobs %d)\n',
            N_REPS_PER_CELL,
            length(GATING_SLOPES) * length(N_VALUES) *
              N_REPS_PER_CELL))
cat(sprintf('Wall budget: %.1f s\n', WALL_BUDGET_SEC))

## Probe: a single rep at the most-expensive cell (gating=1.5, N=70)
## to surface any load-order or hlme errors before paying parallel
## startup cost.
cat('Probe: 1 rep at gating=1.5, N=70 (sequential)...\n')
t_probe0 <- Sys.time()
probe_res <- worker(list(gating_slope = 1.5, N = 70L, rep_idx = 1L))
t_probe <- as.numeric(difftime(Sys.time(), t_probe0,
                               units = 'secs'))
cat(sprintf('  probe elapsed: %.2f s; lcmm conv = %d\n',
            t_probe, probe_res$lcmm$conv))

## Build the full job list. Order: interleave by rep_idx within
## (gating, N) so an early time-abort still yields some reps in
## every cell. We keep the cell-grouped order to make progress
## reporting per-cell-friendly; the time-abort behaviour is
## handled by the chunk loop.
jobs <- make_jobs(GATING_SLOPES, N_VALUES, N_REPS_PER_CELL)
cat(sprintf('Total jobs: %d\n', length(jobs)))

mode_used <- NA_character_
results_list <- NULL
n_done <- 0L

attempt <- function(mc_cores, label) {
  cat(sprintf('Attempting %s (mc.cores = %d, %d jobs)...\n',
              label, mc_cores, length(jobs)))
  t0 <- Sys.time()
  res <- run_parallel_chunked(jobs, mc_cores, t_global,
                              WALL_BUDGET_SEC)
  t1 <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
  cat(sprintf('  %s elapsed: %.1f s; %d jobs complete\n',
              label, t1, res$done))
  res
}

primary <- attempt(MC_CORES_PRIMARY,
                   sprintf('mc.cores = %d', MC_CORES_PRIMARY))
if (!is.null(primary$results) && primary$done > 0L) {
  results_list <- primary$results
  n_done <- primary$done
  mode_used <- sprintf('mclapply mc.cores=%d', MC_CORES_PRIMARY)
} else {
  fb <- attempt(MC_CORES_FALLBACK,
                sprintf('mc.cores = %d', MC_CORES_FALLBACK))
  if (!is.null(fb$results) && fb$done > 0L) {
    results_list <- fb$results
    n_done <- fb$done
    mode_used <- sprintf('mclapply mc.cores=%d',
                         MC_CORES_FALLBACK)
  } else {
    cat(sprintf('Falling back to sequential at %d reps/cell...\n',
                SEQ_REPS))
    jobs_seq <- make_jobs(GATING_SLOPES, N_VALUES, SEQ_REPS)
    seq_res <- run_sequential_capped(jobs_seq, t_global,
                                     WALL_BUDGET_SEC)
    results_list <- seq_res$results
    n_done <- seq_res$done
    mode_used <- 'sequential (fallback)'
  }
}

t_total <- as.numeric(difftime(Sys.time(), t_global,
                               units = 'secs'))
cat(sprintf('\nTotal wall time: %.1f s (%.2f min)\n',
            t_total, t_total / 60))
cat(sprintf('Mode used: %s\n', mode_used))
cat(sprintf('Jobs complete: %d / %d\n', n_done, length(jobs)))

if (n_done == 0L) {
  stop('No replicates completed; nothing to summarise.')
}

## Stitch per-replicate data.
per_rep <- rbindlist(c(
  lapply(results_list, `[[`, 'lme'),
  lapply(results_list, `[[`, 'lcmm')
), use.names = TRUE, fill = TRUE)
setorder(per_rep, gating_slope, N, analysis, rep_idx)

saveRDS(per_rep, file.path(DATA_DIR, '03-latent-replicates.rds'))
cat(sprintf('Per-replicate output -> %s (rows=%d)\n',
            file.path(DATA_DIR, '03-latent-replicates.rds'),
            nrow(per_rep)))

## Morris summary table -----------------------------------------
binom_mcse <- function(rate, n) {
  if (is.na(rate) || is.na(n) || n == 0) return(NA_real_)
  sqrt(rate * (1 - rate) / n)
}

summarise_cell <- function(df_lme, df_lcmm, gating_slope, N_size) {
  dgp_label <- if (gating_slope == 0) 'null' else 'random-slopes'

  n_lme <- sum(!is.na(df_lme$p))
  conv_lme <- mean(df_lme$conv == 1, na.rm = TRUE)
  rej_lme <- mean(df_lme$p < 0.05, na.rm = TRUE)

  n_lcmm <- sum(!is.na(df_lcmm$gating_p))
  conv_lcmm <- mean(df_lcmm$conv == 1, na.rm = TRUE)
  rej_uncond <- mean(df_lcmm$gating_p < 0.05, na.rm = TRUE)

  rej_cond <- mean(
    (df_lcmm$delta_bic > 0) & (df_lcmm$gating_p < 0.05),
    na.rm = TRUE
  )

  rbind(
    data.table(
      dgp = dgp_label, gating_slope = gating_slope, N = N_size,
      analysis_test = 'lme_bmDbc',
      n_reps = n_lme, rejection_rate = rej_lme,
      mcse_rate = binom_mcse(rej_lme, n_lme),
      conv_rate = conv_lme
    ),
    data.table(
      dgp = dgp_label, gating_slope = gating_slope, N = N_size,
      analysis_test = 'lcmm_gating_unconditional',
      n_reps = n_lcmm, rejection_rate = rej_uncond,
      mcse_rate = binom_mcse(rej_uncond, n_lcmm),
      conv_rate = conv_lcmm
    ),
    data.table(
      dgp = dgp_label, gating_slope = gating_slope, N = N_size,
      analysis_test = 'lcmm_gating_BIC_conditional',
      n_reps = n_lcmm, rejection_rate = rej_cond,
      mcse_rate = binom_mcse(rej_cond, n_lcmm),
      conv_rate = conv_lcmm
    )
  )
}

cell_grid <- CJ(gating_slope = GATING_SLOPES, N = N_VALUES,
                sorted = FALSE)
summary_tab <- rbindlist(lapply(seq_len(nrow(cell_grid)),
                                function(i) {
  g <- cell_grid$gating_slope[i]
  n_size <- cell_grid$N[i]
  df_lme  <- per_rep[gating_slope == g & N == n_size &
                       analysis == 'lme']
  df_lcmm <- per_rep[gating_slope == g & N == n_size &
                       analysis == 'lcmm']
  summarise_cell(df_lme, df_lcmm, g, n_size)
}))

summary_path <- file.path(DATA_DIR, '03-latent-summary.txt')
fwrite(summary_tab, summary_path, sep = '\t')
cat(sprintf('Morris summary -> %s\n', summary_path))
print(summary_tab)

## Run log -------------------------------------------------------
log_path <- file.path(DATA_DIR, '03-latent-runlog.txt')
reps_per_cell_tab <- per_rep[analysis == 'lcmm',
                             .(reps = .N), by = .(gating_slope, N)]
conv_per_cell_tab <- per_rep[analysis == 'lcmm',
                             .(conv_rate = mean(conv == 1,
                                                na.rm = TRUE)),
                             by = .(gating_slope, N)]
log_lines <- c(
  sprintf('run_started: %s', format(t_global)),
  sprintf('mode_used: %s', mode_used),
  sprintf('total_wall_sec: %.1f', t_total),
  sprintf('jobs_target: %d', length(jobs)),
  sprintf('jobs_done: %d', n_done),
  sprintf('time_was_binding: %s',
          if (t_total >= WALL_BUDGET_SEC * 0.97) 'TRUE'
          else 'FALSE'),
  sprintf('n_replicate_rows: %d', nrow(per_rep)),
  '',
  'reps per (gating_slope, N):',
  capture.output(print(reps_per_cell_tab)),
  '',
  'lcmm convergence rate per (gating_slope, N):',
  capture.output(print(conv_per_cell_tab))
)
writeLines(log_lines, log_path)
cat(sprintf('Run log -> %s\n', log_path))

## Figure --------------------------------------------------------
plot_dat <- copy(summary_tab)
plot_dat[, ci_lo := pmax(0, rejection_rate - 1.96 * mcse_rate)]
plot_dat[, ci_hi := pmin(1, rejection_rate + 1.96 * mcse_rate)]
plot_dat[, analysis_test := factor(
  analysis_test,
  levels = c('lme_bmDbc',
             'lcmm_gating_unconditional',
             'lcmm_gating_BIC_conditional'),
  labels = c('lme bm:Dbc',
             'lcmm gating Wald (unconditional)',
             'lcmm gating, BIC-conditional')
)]
plot_dat[, N_label := factor(sprintf('N = %d', N),
                             levels = sprintf('N = %d',
                                              sort(N_VALUES)))]

n_per_cell_actual <- max(per_rep[analysis == 'lcmm',
                                 .N, by = .(gating_slope, N)]$N)

p_fig <- ggplot(plot_dat,
                aes(x = gating_slope, y = rejection_rate)) +
  geom_hline(yintercept = 0.05, linetype = 'dashed',
             colour = 'grey50') +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.08) +
  geom_point(size = 2.3) +
  geom_line() +
  facet_grid(N_label ~ analysis_test) +
  scale_x_continuous(breaks = GATING_SLOPES) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = 'Gating slope (heterogeneous-slopes DGP intensity)',
    y = 'Rejection rate (alpha = 0.05)',
    title = sprintf(
      paste0('30-min quick-sim: class-aware vs lme moderation tests',
             ' (n_reps up to %d/cell, %s)'),
      n_per_cell_actual, mode_used
    ),
    caption = sprintf(
      paste0('Reference line at 0.05; error bars are 95%% binomial',
             ' CIs. Total wall: %.1f min.'),
      t_total / 60
    )
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        plot.caption = element_text(size = 8))

fig_path <- file.path(FIG_DIR, '03-latent-power.pdf')
ggsave(fig_path, p_fig, width = 9, height = 5)
cat(sprintf('Figure -> %s\n', fig_path))

cat('\nDone.\n')
