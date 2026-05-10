## quick-sim-prototype-parallel.R
##
## Parallel quick-sim driver for paper 02-carryover-sensitivity,
## 10-minute wall-time budget (HARD).
##
## Cells: trial = OL+BDC; analysis_spec in {A1_binary, A2_matched,
## A3_lagged}; c.bm in {0, 0.45}; N in {35, 70}; t1half_dgp in
## {0.5, 1.0}. 24 cells, target 600 reps/cell, 14,400 fits.
##
## Parallelism: parallel::mclapply with mc.cores = 8 primary, with
## fall-back escalation to mc.cores = 4, then sequential at 300
## reps/cell on the 6 (analysis_spec x c.bm) cells if forking
## fails. Pattern matches sister papers 04, 05, 08.
##
## Time budget: aborts gracefully and writes whatever is finished
## once elapsed > 540 s (90% of 600 s budget). Progress printed
## every 500 fits.
##
## Run from the package root:
##   Rscript analysis/scripts/carryover-sensitivity/\
##           quick-sim-prototype-parallel.R

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
  library(parallel)
})

DATA_DIR <- file.path('analysis', 'data', 'quick-sim')
FIG_DIR  <- file.path('analysis', 'figures', 'quick-sim')
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

SEED_BASE        <- 20260508L
ANALYSIS_SPECS   <- c('A1_binary', 'A2_matched', 'A3_lagged')
C_BM_LEVELS      <- c(0, 0.45)
N_LEVELS         <- c(35L, 70L)
T1HALF_LEVELS    <- c(0.5, 1.0)
N_REPS_PER_CELL  <- 600L
MC_CORES_PRIMARY <- 8L
MC_CORES_FALLBACK <- 4L
SEQ_REPS         <- 300L
WALL_BUDGET_SEC  <- 540          # graceful abort threshold
PROGRESS_EVERY   <- 500L         # in fits, not reps

t_global <- Sys.time()
cat(sprintf('Parallel driver started at %s\n', format(t_global)))
cat(sprintf('Cells: %d analysis_spec x %d c.bm x %d N x %d t1half_dgp = %d\n',
            length(ANALYSIS_SPECS), length(C_BM_LEVELS),
            length(N_LEVELS), length(T1HALF_LEVELS),
            length(ANALYSIS_SPECS) * length(C_BM_LEVELS) *
              length(N_LEVELS) * length(T1HALF_LEVELS)))
cat(sprintf('Target reps per cell: %d\n', N_REPS_PER_CELL))
cat(sprintf('Wall budget: %.0f s\n', WALL_BUDGET_SEC))

## --------------------------------------------------------------
## Trial design and fixed parameter blocks
## --------------------------------------------------------------

td <- buildtrialdesign(
  name_longform  = 'OL+BDC',
  name_shortform = 'OLBDC',
  timepoints     = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptnames    = c('OL1', 'OL2', 'OL3', 'OL4',
                     'BD1', 'BD2', 'BD3', 'BD4'),
  expectancies   = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)

resp_param <- data.table(
  cat  = c('tv', 'pb', 'br'),
  max  = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd   = c(5, 2, 5)
)

baseline_param <- data.table(
  cat = c('bm', 'BL'),
  m   = c(0, 70),
  sd  = c(1, 10)
)

## A2's analysis-side carryover half-life tracks the DGP value
## within each cell so 'matched' remains correctly matched as
## t1half_dgp varies.
make_op_specs <- function(t1half_dgp) {
  list(
    A1_binary = list(useDE = TRUE, t_random_slope = FALSE,
                     full_model_out = FALSE,
                     carryover_t1half = 0,
                     simplecarryover = FALSE,
                     carryover_scalefactor = 1),
    A2_matched = list(useDE = TRUE, t_random_slope = FALSE,
                      full_model_out = FALSE,
                      carryover_t1half = t1half_dgp,
                      simplecarryover = FALSE,
                      carryover_scalefactor = 1),
    A3_lagged = list(useDE = TRUE, t_random_slope = FALSE,
                     full_model_out = FALSE,
                     carryover_t1half = 0,
                     simplecarryover = TRUE,
                     carryover_scalefactor = 1)
  )
}

## --------------------------------------------------------------
## Worker: simulate one dataset, fit all three analysis specs
## --------------------------------------------------------------

simulate_one_dataset <- function(c_bm, N_total, t1half_dgp, seed) {
  set.seed(seed)
  N_per_path <- c(ceiling(N_total / 2), floor(N_total / 2))
  mp1 <- data.table(N = N_per_path[1], c.bm = c_bm,
                    carryover_t1half = t1half_dgp,
                    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
                    c.cf1t = 0.1, c.cfct = 0.05)
  mp2 <- data.table(N = N_per_path[2], c.bm = c_bm,
                    carryover_t1half = t1half_dgp,
                    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
                    c.cf1t = 0.1, c.cfct = 0.05)
  d1 <- generateData(mp1, resp_param, baseline_param,
                     td$trialpaths[[1]],
                     empirical = FALSE,
                     makePositiveDefinite = TRUE)
  d1[, path := 1]
  d2 <- generateData(mp2, resp_param, baseline_param,
                     td$trialpaths[[2]],
                     empirical = FALSE,
                     makePositiveDefinite = TRUE)
  d2[, path := 2]
  d2[, ptID := ptID + max(d1$ptID)]
  rbind(d1, d2)
}

run_one_rep <- function(c_bm, N_size, t1half_dgp, rep_idx) {
  cb_idx <- which(C_BM_LEVELS == c_bm)
  n_idx  <- which(N_LEVELS == N_size)
  th_idx <- which(T1HALF_LEVELS == t1half_dgp)
  ## Same dataset seed across the three analysis specs so within-rep
  ## comparisons are paired (enables McNemar across A1/A2/A3).
  seed_i <- SEED_BASE + rep_idx +
            10000L * cb_idx + 250000L * n_idx +
            5000000L * th_idx
  dat <- tryCatch(
    simulate_one_dataset(c_bm, N_size, t1half_dgp, seed_i),
    error = function(e) NULL
  )
  op_specs <- make_op_specs(t1half_dgp)
  out_list <- vector('list', length(op_specs))
  names(out_list) <- names(op_specs)
  for (spec_name in names(op_specs)) {
    if (is.null(dat)) {
      out_list[[spec_name]] <- data.table(
        analysis_spec = spec_name, c.bm = c_bm, N = N_size,
        t1half_dgp = t1half_dgp, rep_idx = rep_idx,
        beta = NA_real_, betaSE = NA_real_, p = NA_real_,
        converged = FALSE
      )
      next
    }
    fit <- tryCatch(
      lme_analysis(td$trialpaths, dat, op_specs[[spec_name]]),
      error = function(e) NULL
    )
    if (is.null(fit) || is.na(fit$beta)) {
      out_list[[spec_name]] <- data.table(
        analysis_spec = spec_name, c.bm = c_bm, N = N_size,
        t1half_dgp = t1half_dgp, rep_idx = rep_idx,
        beta = NA_real_, betaSE = NA_real_, p = NA_real_,
        converged = FALSE
      )
    } else {
      out_list[[spec_name]] <- data.table(
        analysis_spec = spec_name, c.bm = c_bm, N = N_size,
        t1half_dgp = t1half_dgp, rep_idx = rep_idx,
        beta = fit$beta, betaSE = fit$betaSE, p = fit$p,
        converged = TRUE
      )
    }
  }
  rbindlist(out_list)
}

worker <- function(job) {
  tryCatch(
    run_one_rep(job$c_bm, job$N, job$t1half_dgp, job$rep_idx),
    error = function(e) {
      ## Synthetic 3-row failure record so structural shape is stable
      data.table(
        analysis_spec = ANALYSIS_SPECS,
        c.bm = job$c_bm, N = job$N,
        t1half_dgp = job$t1half_dgp, rep_idx = job$rep_idx,
        beta = NA_real_, betaSE = NA_real_, p = NA_real_,
        converged = FALSE
      )
    }
  )
}

## --------------------------------------------------------------
## Job list: interleave by rep_idx so an early time-abort still
## yields some reps in every cell.
## --------------------------------------------------------------

make_jobs <- function(n_reps) {
  jobs <- list()
  for (r in seq_len(n_reps)) {
    for (cb in C_BM_LEVELS) {
      for (n_size in N_LEVELS) {
        for (th in T1HALF_LEVELS) {
          jobs[[length(jobs) + 1L]] <-
            list(c_bm = cb, N = n_size, t1half_dgp = th,
                 rep_idx = r)
        }
      }
    }
  }
  jobs
}

is_structural_failure <- function(x) {
  inherits(x, 'try-error') || is.null(x) ||
    !is.data.table(x) || nrow(x) != 3L
}

## Time-bounded chunked executor; checks elapsed time between
## chunks. With ~14400 fits over ~24 DGP cells, chunk_size = 60
## means ~2.5 seconds of work per core per chunk at primary 8.
run_parallel_chunked <- function(jobs, mc_cores, t_anchor,
                                 wall_budget,
                                 chunk_size = 60L,
                                 progress_every = PROGRESS_EVERY) {
  results <- vector('list', length(jobs))
  total <- length(jobs)
  done <- 0L
  fits_done <- 0L
  next_progress <- progress_every
  i <- 1L
  while (i <= total) {
    elapsed <- as.numeric(difftime(Sys.time(), t_anchor,
                                   units = 'secs'))
    if (elapsed > wall_budget) {
      cat(sprintf(
        '  TIME ABORT: elapsed %.1f s > %.0f s budget; %d/%d jobs done\n',
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
    if (length(fail_idx) > 0L) {
      cat(sprintf(
        '  STRUCTURAL FAIL: %d/%d in chunk starting at %d\n',
        length(fail_idx), length(out_chunk), i))
      return(list(results = NULL, done = done))
    }
    for (k in seq_along(out_chunk)) {
      results[[i + k - 1L]] <- out_chunk[[k]]
    }
    done <- done + length(out_chunk)
    fits_done <- fits_done + 3L * length(out_chunk)
    if (fits_done >= next_progress) {
      pct <- 100 * done / total
      el <- as.numeric(difftime(Sys.time(), t_anchor,
                                units = 'secs'))
      eta <- if (done > 0) el * (total - done) / done
             else NA_real_
      cat(sprintf(
        '  progress: jobs=%d/%d (%.1f%%) fits=%d elapsed=%.1f s eta=%.1f s\n',
        done, total, pct, fits_done, el, eta))
      next_progress <- next_progress + progress_every
    }
    i <- j_end + 1L
  }
  list(results = results[seq_len(done)], done = done)
}

run_sequential_capped <- function(jobs, t_anchor, wall_budget,
                                  progress_every = PROGRESS_EVERY) {
  results <- vector('list', length(jobs))
  done <- 0L
  fits_done <- 0L
  next_progress <- progress_every
  for (i in seq_along(jobs)) {
    elapsed <- as.numeric(difftime(Sys.time(), t_anchor,
                                   units = 'secs'))
    if (elapsed > wall_budget) {
      cat(sprintf(
        '  TIME ABORT (sequential): elapsed %.1f s; %d/%d jobs done\n',
        elapsed, done, length(jobs)))
      break
    }
    results[[i]] <- worker(jobs[[i]])
    done <- done + 1L
    fits_done <- fits_done + 3L
    if (fits_done >= next_progress) {
      el <- as.numeric(difftime(Sys.time(), t_anchor,
                                units = 'secs'))
      cat(sprintf(
        '  progress: jobs=%d/%d fits=%d elapsed=%.1f s\n',
        done, length(jobs), fits_done, el))
      next_progress <- next_progress + progress_every
    }
  }
  list(results = results[seq_len(done)], done = done)
}

## --------------------------------------------------------------
## Probe: surface load-order or DGP errors before paying parallel
## startup cost.
## --------------------------------------------------------------

cat('Probe: 1 rep at c.bm=0.45, N=70, t1half=1.0 (sequential)...\n')
t_probe0 <- Sys.time()
probe_res <- worker(list(c_bm = 0.45, N = 70L,
                         t1half_dgp = 1.0, rep_idx = 1L))
t_probe <- as.numeric(difftime(Sys.time(), t_probe0,
                               units = 'secs'))
cat(sprintf('  probe elapsed: %.2f s; rows = %d; conv = %d/3\n',
            t_probe, nrow(probe_res), sum(probe_res$converged)))

## --------------------------------------------------------------
## Build job list and run with primary -> fallback -> sequential
## --------------------------------------------------------------

jobs <- make_jobs(N_REPS_PER_CELL)
cat(sprintf('Total jobs (datasets): %d -> %d fits\n',
            length(jobs), 3L * length(jobs)))

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
  mode_used <- sprintf('mclapply mc.cores=%d',
                       MC_CORES_PRIMARY)
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
    jobs_seq <- make_jobs(SEQ_REPS)
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
cat(sprintf('Datasets complete: %d / %d\n', n_done,
            length(jobs)))

if (n_done == 0L) {
  stop('No replicates completed; nothing to summarise.')
}

## Stitch per-replicate data.
replicates <- rbindlist(results_list, use.names = TRUE,
                        fill = TRUE)
setorder(replicates, c.bm, N, t1half_dgp, analysis_spec, rep_idx)

saveRDS(replicates,
        file.path(DATA_DIR, '02-carryover-replicates.rds'))
cat(sprintf('Per-replicate output -> %s (rows=%d)\n',
            file.path(DATA_DIR, '02-carryover-replicates.rds'),
            nrow(replicates)))

## --------------------------------------------------------------
## Morris ADEMP summary
## --------------------------------------------------------------

binom_mcse <- function(rate, n) {
  if (is.na(rate) || is.na(n) || n == 0) return(NA_real_)
  sqrt(rate * (1 - rate) / n)
}

summary_dt <- replicates[, {
  conv <- !is.na(beta)
  n_total <- .N
  n_conv  <- sum(conv)
  pwr     <- if (n_conv > 0) mean(p[conv] < 0.05, na.rm = TRUE)
             else NA_real_
  mcse_pwr <- binom_mcse(pwr, n_conv)
  m_b     <- if (n_conv > 0) mean(beta[conv], na.rm = TRUE)
             else NA_real_
  s_b     <- if (n_conv > 1) sd(beta[conv], na.rm = TRUE)
             else NA_real_
  mcse_b  <- if (n_conv > 0) s_b / sqrt(n_conv) else NA_real_
  list(n_reps = n_total,
       power = pwr, mcse_power = mcse_pwr,
       mean_beta = m_b, sd_beta = s_b,
       mcse_mean_beta = mcse_b,
       conv_rate = if (n_total > 0) n_conv / n_total
                   else NA_real_)
}, by = .(analysis_spec, c.bm, N, t1half_dgp)]

setorder(summary_dt, c.bm, N, t1half_dgp, analysis_spec)

## Write a Morris-style summary with header and run metadata.
summary_path <- file.path(DATA_DIR, '02-carryover-summary.txt')
header_lines <- c(
  sprintf('# Morris ADEMP summary: 02-carryover-sensitivity'),
  sprintf('# run_started: %s', format(t_global)),
  sprintf('# total_wall_sec: %.1f', t_total),
  sprintf('# mode_used: %s', mode_used),
  sprintf('# datasets_target: %d', length(jobs)),
  sprintf('# datasets_done: %d', n_done),
  sprintf('# fits_done: %d', 3L * n_done),
  sprintf('# time_was_binding: %s',
          if (t_total >= WALL_BUDGET_SEC * 0.97) 'TRUE'
          else 'FALSE'),
  '#'
)
writeLines(header_lines, summary_path)
fwrite(summary_dt, summary_path, sep = '\t', append = TRUE,
       col.names = TRUE)
cat(sprintf('Morris summary -> %s\n', summary_path))
print(summary_dt)

## --------------------------------------------------------------
## Figure: power across analysis_spec faceted by N x t1half_dgp,
## c.bm = 0.45 only.
## --------------------------------------------------------------

plot_dt <- copy(summary_dt[c.bm == 0.45])
plot_dt[, ci_lo := pmax(0, power - 1.96 * mcse_power)]
plot_dt[, ci_hi := pmin(1, power + 1.96 * mcse_power)]
plot_dt[, N_lbl := factor(sprintf('N = %d', N),
                          levels = sprintf('N = %d',
                                           sort(N_LEVELS)))]
plot_dt[, th_lbl := factor(sprintf('DGP t1/2 = %.1f wk',
                                   t1half_dgp),
                           levels = sprintf('DGP t1/2 = %.1f wk',
                                            sort(T1HALF_LEVELS)))]
plot_dt[, analysis_spec := factor(analysis_spec,
                                  levels = ANALYSIS_SPECS)]

n_per_cell_actual <- max(replicates[c.bm == 0.45,
                                    .N,
                                    by = .(analysis_spec, N,
                                           t1half_dgp)]$N)

p_fig <- ggplot(plot_dt,
                aes(x = analysis_spec, y = power,
                    ymin = ci_lo, ymax = ci_hi)) +
  geom_col(fill = 'grey75', colour = 'grey30') +
  geom_errorbar(width = 0.2) +
  geom_text(aes(label = sprintf('%.3f', power)),
            vjust = -0.6, size = 3.0) +
  facet_grid(N_lbl ~ th_lbl) +
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  labs(x = 'Analysis-side carryover specification',
       y = 'Power',
       title = paste0(
         'Power by analysis-side carryover specification ',
         '(c.bm = 0.45)'),
       subtitle = sprintf(
         paste0('OL+BDC; faceted by N x DGP carryover ',
                't1/2; reps per cell up to %d (%s)'),
         n_per_cell_actual, mode_used),
       caption = sprintf(
         paste0('Error bars are 95%% binomial CIs. ',
                'Total wall: %.1f min.'),
         t_total / 60)) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11),
        plot.caption = element_text(size = 8))

fig_path <- file.path(FIG_DIR, '02-carryover-power.pdf')
ggsave(fig_path, p_fig, width = 8, height = 6)
cat(sprintf('Figure -> %s\n', fig_path))

cat('\nDone.\n')

invisible(NULL)
