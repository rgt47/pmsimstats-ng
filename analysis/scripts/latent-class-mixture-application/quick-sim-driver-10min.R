## Parallel quick-sim driver for paper 03-latent-class-mixture-
## application. 10-minute hard wall budget (540 s working budget),
## six cells = 3 gating slopes x 2 N values, three test forms per
## cell (lme bm:Dbc, lcmm gating Wald unconditional, lcmm BIC-
## conditional). Uses parallel::mclapply with mc.cores = 8 as a
## hard requirement.
##
## Run from the package root:
##   Rscript analysis/scripts/latent-class-mixture-application/\
##           quick-sim-driver-10min.R

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

SEED_BASE         <- 20260508L
GATING_SLOPES     <- c(0, 1.0, 1.5)
N_LEVELS          <- c(35L, 70L)
N_REPS_PER_CELL   <- 150L
N_REPS_FALLBACK   <- 100L
MC_CORES_PRIMARY  <- 8L
MC_CORES_FALLBACK <- 4L
WALL_BUDGET_SEC   <- 540

## --------------------------------------------------------------
## Trial design: OL+BDC (open-label run-in then blinded
## discontinuation), N split across two paths. The user prompt
## describes this as 'Hybrid OL+BDC'; the codebase OLBDC preset
## is the matching definition.
## --------------------------------------------------------------

td <- buildtrialdesign(
  name_longform  = 'OL+BDC',
  name_shortform = 'OLBDC',
  timepoints     = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptnames    = c('OL1', 'OL2', 'OL3', 'OL4',
                     'BD1', 'BD2', 'BD3', 'BD4'),
  expectancies   = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
  ondrug         = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
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

## --------------------------------------------------------------
## Per-cell modelparam constructor. Splits N across the two paths
## for the OLBDC design.
## --------------------------------------------------------------

split_n <- function(N) {
  c(ceiling(N / 2), floor(N / 2))
}

make_modelparam <- function(gating_slope, N_path) {
  if (gating_slope == 0) {
    data.table(
      N = N_path, c.bm = 0,
      carryover_t1half = 1.0,
      c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
      c.cf1t = 0.1, c.cfct = 0.05
    )
  } else {
    data.table(
      N = N_path, c.bm = 0,
      carryover_t1half = 1.0,
      c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
      c.cf1t = 0.1, c.cfct = 0.05,
      beta_R_factor = 1.0, beta_NR_factor = 0.0,
      gating_intercept = 0.0, gating_slope = gating_slope
    )
  }
}

simulate_dat <- function(gating_slope, N, rep_idx) {
  cell_idx <- which(GATING_SLOPES == gating_slope)
  N_idx    <- which(N_LEVELS == N)
  set.seed(SEED_BASE + rep_idx +
           1000L * cell_idx + 100000L * N_idx)
  N_per_path <- split_n(N)

  mk_path_dat <- function(path_i) {
    mp <- make_modelparam(gating_slope, N_per_path[path_i])
    if (gating_slope == 0) {
      d <- generateData(
        modelparam = mp, respparam = rp, blparam = bp,
        trialdesign = td$trialpaths[[path_i]],
        empirical = FALSE, makePositiveDefinite = TRUE,
        dgp_architecture = 'mvn'
      )
    } else {
      d <- generate_random_slopes_data(
        modelparam = mp, respparam = rp, blparam = bp,
        trialdesign = td$trialpaths[[path_i]],
        empirical = FALSE, makePositiveDefinite = TRUE
      )
    }
    d[, path := path_i]
    d
  }

  d1 <- mk_path_dat(1)
  d2 <- mk_path_dat(2)
  d2[, ptID := ptID + max(d1$ptID)]
  rbind(d1, d2, fill = TRUE)
}

run_one <- function(gating_slope, N, rep_idx) {
  dat <- simulate_dat(gating_slope, N, rep_idx)

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
      gating_slope = gating_slope, N = N, analysis = 'lme',
      rep_idx = rep_idx,
      beta = out_lme$beta, betaSE = out_lme$betaSE, p = out_lme$p,
      gating_p = NA_real_, within_p = NA_real_,
      delta_bic = NA_real_,
      conv = as.integer(!is.na(out_lme$p)), elapsed_s = t_lme
    ),
    lcmm = data.table(
      gating_slope = gating_slope, N = N, analysis = 'lcmm',
      rep_idx = rep_idx,
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
  tryCatch(
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
}

## Build the job list as gating_slope x N x rep_idx tuples.
make_jobs <- function(gating_slopes, N_levels, n_reps) {
  jobs <- list()
  for (g in gating_slopes) {
    for (N in N_levels) {
      for (r in seq_len(n_reps)) {
        jobs[[length(jobs) + 1L]] <- list(
          gating_slope = g, N = N, rep_idx = r
        )
      }
    }
  }
  jobs
}

t_global <- Sys.time()
cat(sprintf('Latent-class quick-sim 10-min driver started at %s\n',
            format(t_global)))
cat(sprintf('Cells: gating_slope = %s ; N = %s\n',
            paste(GATING_SLOPES, collapse = ', '),
            paste(N_LEVELS, collapse = ', ')))
cat(sprintf('Wall budget: %d s\n', WALL_BUDGET_SEC))

## Probe: one rep at gating=0, N=35 (cheapest cell), sequential.
cat('Probe: 1 rep at gating=0, N=35 (sequential)...\n')
t_probe0 <- Sys.time()
probe_res <- worker(list(gating_slope = 0, N = 35L, rep_idx = 1L))
t_probe <- as.numeric(difftime(Sys.time(), t_probe0, units = 'secs'))
cat(sprintf('  probe elapsed: %.2f s; lcmm conv = %d\n',
            t_probe, probe_res$lcmm$conv))

## Estimate runtime. The probe cell is the cheapest (gating=0,
## N=35); N=70 cells with gating>0 cost more. Apply a 1.6x safety
## multiplier so we under-promise on rep counts and exit with
## time to spare for plotting/IO.
COST_SAFETY_MULT <- 1.6
target_reps <- N_REPS_PER_CELL
cell_count  <- length(GATING_SLOPES) * length(N_LEVELS)
predicted_serial_s <- COST_SAFETY_MULT * t_probe *
                      target_reps * cell_count
predicted_parallel_s <- predicted_serial_s / MC_CORES_PRIMARY
cat(sprintf(
  'Predicted parallel runtime at %d reps/cell, mc.cores=%d: %.0f s (%.1fx safety)\n',
  target_reps, MC_CORES_PRIMARY, predicted_parallel_s,
  COST_SAFETY_MULT
))
budget_remaining <- WALL_BUDGET_SEC - t_probe -
                    as.numeric(difftime(Sys.time(), t_global, 'secs'))
if (predicted_parallel_s > 0.85 * budget_remaining) {
  cat(sprintf(
    '  prediction exceeds 85%% of remaining budget; falling back to %d reps/cell\n',
    N_REPS_FALLBACK
  ))
  target_reps <- N_REPS_FALLBACK
}

mode_used <- NA_character_
results <- NULL

attempt_parallel <- function(jobs, mc_cores, label) {
  cat(sprintf('Attempting %s (mc.cores = %d, %d jobs)...\n',
              label, mc_cores, length(jobs)))
  t0 <- Sys.time()
  out <- tryCatch(
    mclapply(jobs, worker,
             mc.cores = mc_cores,
             mc.preschedule = FALSE,
             mc.set.seed = FALSE),
    error = function(e) {
      cat(sprintf('  ERROR: %s\n', conditionMessage(e)))
      NULL
    }
  )
  t1 <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
  cat(sprintf('  %s elapsed: %.1f s (%.2f min)\n',
              label, t1, t1 / 60))
  if (is.null(out)) return(NULL)
  err_idx <- which(vapply(out, function(x) {
    inherits(x, 'try-error') || is.null(x) ||
      !is.list(x) || is.null(x$lme) || is.null(x$lcmm)
  }, logical(1)))
  if (length(err_idx) > 0) {
    cat(sprintf('  WARN: %d/%d jobs had structural failures\n',
                length(err_idx), length(out)))
    cat(sprintf('  example: %s\n',
                paste(deparse(out[[err_idx[1]]])[1:3],
                      collapse = ' ')))
    return(NULL)
  }
  list(out = out, elapsed = t1)
}

## Build job list and run with mc.cores=8 (primary) -> 4 (fallback)
## -> sequential at reduced rep count if both fail.
jobs <- make_jobs(GATING_SLOPES, N_LEVELS, target_reps)
cat(sprintf('Total jobs: %d\n', length(jobs)))

primary <- attempt_parallel(jobs, MC_CORES_PRIMARY,
                             sprintf('mc.cores = %d', MC_CORES_PRIMARY))
if (!is.null(primary)) {
  results <- primary$out
  mode_used <- sprintf('mclapply mc.cores=%d', MC_CORES_PRIMARY)
} else {
  cat('Primary attempt failed; trying fallback core count...\n')
  fallback <- attempt_parallel(jobs, MC_CORES_FALLBACK,
                                sprintf('mc.cores = %d',
                                        MC_CORES_FALLBACK))
  if (!is.null(fallback)) {
    results <- fallback$out
    mode_used <- sprintf('mclapply mc.cores=%d', MC_CORES_FALLBACK)
  } else {
    cat('Falling back to sequential at reduced reps...\n')
    SEQ_REPS <- 30L
    jobs_seq <- make_jobs(GATING_SLOPES, N_LEVELS, SEQ_REPS)
    t_seq0 <- Sys.time()
    results <- lapply(jobs_seq, worker)
    t_seq <- as.numeric(difftime(Sys.time(), t_seq0, units = 'secs'))
    cat(sprintf('  sequential elapsed: %.1f s\n', t_seq))
    mode_used <- 'sequential (fallback)'
  }
}

t_total <- as.numeric(difftime(Sys.time(), t_global, units = 'secs'))
cat(sprintf('\nTotal wall time: %.1f s (%.2f min)\n',
            t_total, t_total / 60))
cat(sprintf('Mode used: %s\n', mode_used))

## --------------------------------------------------------------
## Stitch per-replicate data
## --------------------------------------------------------------

per_rep <- rbindlist(c(
  lapply(results, `[[`, 'lme'),
  lapply(results, `[[`, 'lcmm')
), use.names = TRUE, fill = TRUE)
setorder(per_rep, gating_slope, N, analysis, rep_idx)

saveRDS(per_rep, file.path(DATA_DIR, '03-latent-replicates.rds'))
cat(sprintf('Per-replicate output -> %s (rows=%d)\n',
            file.path(DATA_DIR, '03-latent-replicates.rds'),
            nrow(per_rep)))

## --------------------------------------------------------------
## Morris summary (per gating x N x test form)
## --------------------------------------------------------------

binom_mcse <- function(rate, n) {
  if (is.na(rate) || is.na(n) || n == 0) return(NA_real_)
  sqrt(rate * (1 - rate) / n)
}

summarise_cell <- function(df_lme, df_lcmm, gating_slope, N) {
  dgp_label <- if (gating_slope == 0) 'null' else 'random-slopes'

  n_lme <- sum(!is.na(df_lme$p))
  conv_lme <- mean(df_lme$conv == 1, na.rm = TRUE)
  rej_lme  <- mean(df_lme$p < 0.05, na.rm = TRUE)

  n_lcmm <- sum(!is.na(df_lcmm$gating_p))
  conv_lcmm <- mean(df_lcmm$conv == 1, na.rm = TRUE)
  rej_uncond <- mean(df_lcmm$gating_p < 0.05, na.rm = TRUE)

  rej_cond <- mean(
    (df_lcmm$delta_bic > 0) & (df_lcmm$gating_p < 0.05),
    na.rm = TRUE
  )

  rbind(
    data.table(
      dgp = dgp_label, gating_slope = gating_slope, N = N,
      analysis_test = 'lme_bmDbc',
      n_reps = n_lme, rejection_rate = rej_lme,
      mcse_rate = binom_mcse(rej_lme, n_lme),
      conv_rate = conv_lme
    ),
    data.table(
      dgp = dgp_label, gating_slope = gating_slope, N = N,
      analysis_test = 'lcmm_gating_unconditional',
      n_reps = n_lcmm, rejection_rate = rej_uncond,
      mcse_rate = binom_mcse(rej_uncond, n_lcmm),
      conv_rate = conv_lcmm
    ),
    data.table(
      dgp = dgp_label, gating_slope = gating_slope, N = N,
      analysis_test = 'lcmm_gating_BIC_conditional',
      n_reps = n_lcmm, rejection_rate = rej_cond,
      mcse_rate = binom_mcse(rej_cond, n_lcmm),
      conv_rate = conv_lcmm
    )
  )
}

cells <- CJ(g = GATING_SLOPES, N = N_LEVELS, sorted = FALSE)
summary_tab <- rbindlist(lapply(seq_len(nrow(cells)), function(i) {
  g <- cells$g[i]; N_i <- cells$N[i]
  df_lme  <- per_rep[gating_slope == g & N == N_i & analysis == 'lme']
  df_lcmm <- per_rep[gating_slope == g & N == N_i & analysis == 'lcmm']
  summarise_cell(df_lme, df_lcmm, g, N_i)
}))

summary_path <- file.path(DATA_DIR, '03-latent-summary.txt')
fwrite(summary_tab, summary_path, sep = '\t')
cat(sprintf('Morris summary -> %s\n', summary_path))
print(summary_tab)

## --------------------------------------------------------------
## Run log
## --------------------------------------------------------------

reps_per_cell_obs <- per_rep[analysis == 'lcmm',
                              .(n = .N), by = .(gating_slope, N)]

log_path <- file.path(DATA_DIR, '03-latent-runlog.txt')
writeLines(c(
  sprintf('run_started: %s', format(t_global)),
  sprintf('mode_used: %s', mode_used),
  sprintf('target_reps_per_cell: %d', target_reps),
  sprintf('reps_per_cell_observed: %s',
          paste(sprintf('g=%g N=%d:%d',
                        reps_per_cell_obs$gating_slope,
                        reps_per_cell_obs$N,
                        reps_per_cell_obs$n),
                collapse = '; ')),
  sprintf('total_wall_sec: %.1f', t_total),
  sprintf('n_replicate_rows: %d', nrow(per_rep)),
  sprintf('probe_sec: %.2f', t_probe),
  sprintf('lcmm_conv_overall: %.3f',
          mean(per_rep[analysis == 'lcmm']$conv, na.rm = TRUE)),
  sprintf('lme_conv_overall: %.3f',
          mean(per_rep[analysis == 'lme']$conv, na.rm = TRUE))
), log_path)
cat(sprintf('Run log -> %s\n', log_path))

## --------------------------------------------------------------
## Figure: rejection-rate panels by analysis_test x N
## --------------------------------------------------------------

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
plot_dat[, N_label := factor(paste0('N = ', N),
                              levels = paste0('N = ', N_LEVELS))]

n_per_cell_actual <- max(per_rep[analysis == 'lcmm',
                                 .N, by = .(gating_slope, N)]$N)

p_fig <- ggplot(plot_dat,
                aes(x = gating_slope, y = rejection_rate)) +
  geom_hline(yintercept = 0.05, linetype = 'dashed',
             colour = 'grey50') +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.08) +
  geom_point(size = 2.5) +
  geom_line() +
  facet_grid(N_label ~ analysis_test) +
  scale_x_continuous(breaks = GATING_SLOPES) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = 'Gating slope (heterogeneous-slopes DGP intensity)',
    y = 'Rejection rate (alpha = 0.05)',
    title = sprintf(
      'Quick-sim: class-aware vs lme moderation tests (n_reps up to %d/cell, %s)',
      n_per_cell_actual, mode_used
    ),
    caption = sprintf(
      'Reference line at 0.05; error bars are 95%% binomial CIs. Total wall: %.1f min. OL+BDC design.',
      t_total / 60
    )
  ) +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        plot.caption = element_text(size = 8),
        strip.text = element_text(size = 9))

fig_path <- file.path(FIG_DIR, '03-latent-power.pdf')
ggsave(fig_path, p_fig, width = 9, height = 5.5)
cat(sprintf('Figure -> %s\n', fig_path))

cat('\nDone.\n')
