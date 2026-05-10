# 07-gompertz-evaluation.R
#
# Prototype simulation for paper 07-gompertz-evaluation.
#
# Compares four BR-trajectory families (gompertz, logistic,
# hyperbolic_tangent, piecewise_linear_breakpoint) crossed with
# two DGP architectures (Architecture B = MVN moderation;
# Architecture A = mean moderation).  All families share the same
# asymptote (maxr) and are matched (calibrated) to the Gompertz
# value at the anchor week 8.
#
# Trajectory families:
#   gompertz:   y = maxr * exp(-disp*exp(-rate*t))  (rescaled, baseline)
#   logistic:   y = maxr / (1 + exp(-log_rate*(t - t0)))
#                 implemented previously, kept here for parity
#   hyperbolic_tangent:
#               y = (maxr/2) * (1 + tanh(rate*(t - t0)))
#                 NEWLY IMPLEMENTED in this driver
#   piecewise_linear_breakpoint:
#               y = pmin(slope*t, maxr)  with slope = maxr / t_bp
#                 NEWLY IMPLEMENTED in this driver
#
# All three non-Gompertz families share maxr with Gompertz and are
# calibrated to satisfy y(t_anchor=8) = gompertz_anchor; this matches
# their effective on-drug BR mean at the canonical Hendrickson anchor.
# The shape of the trajectory differs between families, which is the
# scientific quantity of interest.
#
# Cells: family {4} x architecture {2} = 8.  Target 600 reps/cell
# (4800 fits) within a 540-s wall-clock budget.  Forked parallelism
# via parallel::mclapply with mc.cores = 8 (HARD requirement, fixed,
# not auto-detected).  Graceful abort: residual cells run with
# whatever budget remains.
#
# Author: pmsimstats team
# Last updated: 2026-05-08

suppressMessages({
  devtools::load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
  library(parallel)
})

set.seed(20260507)

out_data_dir <- 'analysis/data/quick-sim'
out_fig_dir  <- 'analysis/figures/quick-sim'
dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig_dir,  recursive = TRUE, showWarnings = FALSE)

# -- Trial design (OL+BDC hybrid) ------------------------------------

td <- buildtrialdesign(
  name_longform  = 'open label+blinded discontinuation',
  name_shortform = 'OL+BDC',
  timepoints     = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptnames    = c('OL1', 'OL2', 'OL3', 'OL4',
                     'BD1', 'BD2', 'BD3', 'BD4'),
  expectancies   = c(1, 1, 1, 1, .5, .5, .5, .5),
  ondrug         = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)
n_paths <- length(td$trialpaths)

# -- Parameter sets (Hendrickson canonical) --------------------------

data(extracted_bp); data(extracted_rp)
bp <- extracted_bp
rp <- extracted_rp

mp <- data.table(
  N = 35, c.bm = 0.45,
  carryover_t1half = 0,
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
  c.cf1t = 0.2, c.cfct = 0.1
)
op <- list(useDE = FALSE, t_random_slope = FALSE,
           full_model_out = FALSE, carryover_t1half = 0,
           simplecarryover = FALSE, carryover_scalefactor = 1)

# -- Trajectory family definitions and calibration -------------------
#
# Anchor at t = 8 weeks, matching gomp_br(8).  All families share
# the Gompertz asymptote (maxr).  Each non-Gompertz family is a
# 3-parameter form whose free parameter is solved analytically to
# satisfy the anchor constraint.

br_pars <- as.list(rp[cat == 'br'])
maxr_v  <- br_pars$max

gomp_br <- function(t) {
  modgompertz(t, br_pars$max, br_pars$disp, br_pars$rate)
}
t_anchor <- 8
g_anchor <- gomp_br(t_anchor)

# Logistic: anchor solves t0 given log_rate
log_rate <- 0.5
log_t0   <- t_anchor + log(maxr_v / g_anchor - 1) / log_rate
logistic_br <- function(t) {
  maxr_v / (1 + exp(-log_rate * (t - log_t0)))
}

# Hyperbolic tangent: y = (maxr/2)*(1 + tanh(rate*(t-t0)))
# Anchor: tanh(rate*(t_anchor-t0)) = 2*g_anchor/maxr - 1
htan_rate <- 0.5
htan_t0   <- t_anchor - atanh(2 * g_anchor / maxr_v - 1) / htan_rate
htan_br <- function(t) {
  (maxr_v / 2) * (1 + tanh(htan_rate * (t - htan_t0)))
}

# Piecewise-linear breakpoint: linear ramp from 0 at t=0 up to maxr
# at t = t_bp, then flat at maxr.  Anchor sets slope, hence t_bp.
pl_slope <- g_anchor / t_anchor
pl_tbp   <- maxr_v / pl_slope
pl_br <- function(t) {
  pmin(pl_slope * t, maxr_v)
}

cat('-- Trajectory family calibration --\n')
cat(sprintf('  Anchor week %.1f, gompertz BR(anchor) = %.4f, maxr = %.4f\n',
            t_anchor, g_anchor, maxr_v))
cat(sprintf('  logistic:           log_rate = %.3f, t0 = %.3f\n',
            log_rate, log_t0))
cat(sprintf('  hyperbolic_tangent: rate = %.3f, t0 = %.3f\n',
            htan_rate, htan_t0))
cat(sprintf('  piecewise_linear:   slope = %.4f, t_bp = %.3f wks\n',
            pl_slope, pl_tbp))

stopifnot(abs(logistic_br(t_anchor) - g_anchor) < 1e-8)
stopifnot(abs(htan_br(t_anchor) - g_anchor) < 1e-8)
stopifnot(abs(pl_br(t_anchor) - g_anchor) < 1e-8)

family_fns <- list(
  gompertz                    = gomp_br,
  logistic                    = logistic_br,
  hyperbolic_tangent          = htan_br,
  piecewise_linear_breakpoint = pl_br
)

# Per-path BR mean shifts for each family.  Each shift is a vector
# of length nP.  Family 'gompertz' has shift identically zero.
br_shift_by_family_path <- list()
for (fam in names(family_fns)) {
  br_shift_by_family_path[[fam]] <- vector('list', n_paths)
  for (g in seq_len(n_paths)) {
    tp_path <- td$trialpaths[[g]]
    tod_g   <- tp_path$tod
    shift_g <- family_fns[[fam]](tod_g) - gomp_br(tod_g)
    names(shift_g) <- tp_path$timeptnames
    br_shift_by_family_path[[fam]][[g]] <- shift_g
  }
}

# -- Sigma cache: one per architecture x path ------------------------
# Architecture changes the BM-BR cross-correlation block, so caches
# are not reusable across architectures.

architectures <- c('mvn', 'mean_moderation')
cs_by_arch_path <- list()
for (arch in architectures) {
  cs_by_arch_path[[arch]] <- vector('list', n_paths)
  for (g in seq_len(n_paths)) {
    cs_by_arch_path[[arch]][[g]] <- buildSigma(
      mp, rp, bp, td$trialpaths[[g]],
      makePositiveDefinite = TRUE,
      dgp_architecture = arch
    )
  }
}

# -- Single-replicate driver -----------------------------------------

run_one <- function(family, architecture, seed) {
  set.seed(seed)
  shifts <- br_shift_by_family_path[[family]]
  cs_set <- cs_by_arch_path[[architecture]]

  dats <- vector('list', n_paths)
  for (g in seq_len(n_paths)) {
    tp_path <- td$trialpaths[[g]]
    dat <- generateData(mp, rp, bp, tp_path,
                        empirical = FALSE,
                        makePositiveDefinite = TRUE,
                        cached_sigma = cs_set[[g]],
                        dgp_architecture = architecture)
    if (family != 'gompertz') {
      shift_g <- shifts[[g]]
      for (tp in tp_path$timeptnames) {
        if (shift_g[tp] != 0) {
          br_col <- paste0(tp, '.br')
          tv_col <- paste0(tp, '.tv')
          pb_col <- paste0(tp, '.pb')
          d_col  <- paste0('D_', tp)
          dat[, (br_col) := get(br_col) + shift_g[tp]]
          dat[, (d_col)  := get(tv_col) + get(pb_col) + get(br_col)]
          dat[, (tp)     := BL - get(d_col)]
        }
      }
    }
    dat[, path := g]
    dats[[g]] <- dat
  }
  dat <- rbindlist(dats, fill = TRUE)
  res <- lme_analysis(td$trialpaths, dat, op)

  data.table(
    family       = family,
    architecture = architecture,
    rep_idx      = seed,
    beta_bmDbc   = res$beta,
    betaSE_bmDbc = res$betaSE,
    p_bmDbc      = res$p,
    issingular   = res$issingular,
    converged    = !is.na(res$beta)
  )
}

# -- Run plan: cells x reps with graceful abort ----------------------

reps_target  <- 600
families     <- names(family_fns)
n_cells      <- length(families) * length(architectures)

# Distinct seed offsets per (family, architecture) cell so that
# replicates are independent across cells but reproducible.
seed_base <- 700000L
seed_offsets <- list()
k <- 0L
for (fam in families) for (arch in architectures) {
  k <- k + 1L
  seed_offsets[[paste(fam, arch, sep = '|')]] <-
    seed_base + (k - 1L) * 100000L
}

# Time budget and parallelism.
# HARD REQUIREMENT: mc.cores fixed at 8 (sister papers 04, 05, 08 used
# 8-9 core mclapply and finished cleanly; the previous detectCores()-1
# strategy caused org-usage exhaustion on a 30-min run).
t_start    <- Sys.time()
budget_s   <- 540
deadline   <- t_start + budget_s
n_workers  <- 8L

cat(sprintf('\n-- Run plan --\n'))
cat(sprintf('  Cells: %d (%d families x %d architectures)\n',
            n_cells, length(families), length(architectures)))
cat(sprintf('  Target reps/cell: %d  -> total fits target: %d\n',
            reps_target, reps_target * n_cells))
cat(sprintf('  Wall-clock budget: %d s\n', budget_s))
cat(sprintf('  Forked workers: %d\n', n_workers))

run_log <- character()
log_msg <- function(msg) {
  ts <- format(Sys.time(), '%H:%M:%S')
  line <- sprintf('[%s] %s', ts, msg)
  run_log[[length(run_log) + 1L]] <<- line
  cat(line, '\n')
}

log_msg(sprintf('Start at %s, deadline %s',
                format(t_start, '%H:%M:%S'),
                format(deadline, '%H:%M:%S')))

# Time per fit estimate (rough): use first cell's first 8 fits
# sequentially to size remaining work.
log_msg('Pilot timing (8 fits, sequential, gompertz/mvn) ...')
t_pilot_start <- Sys.time()
pilot <- rbindlist(lapply(seq_len(8), function(i)
  run_one('gompertz', 'mvn', seed_base + i)))
t_pilot <- as.numeric(Sys.time() - t_pilot_start, units = 'secs') / 8
log_msg(sprintf('Per-fit (sequential): %.3f s', t_pilot))

# Round-robin schedule: chunk reps in batches of `chunk_size` per
# cell, dispatch via mclapply, then check deadline.  This keeps a
# graceful abort tractable: if the deadline is hit mid-cell, we
# drop the unfinished tail of that batch.

chunk_size  <- 50L
all_results <- list()

cells <- expand.grid(family = families,
                     architecture = architectures,
                     stringsAsFactors = FALSE)
n_cells <- nrow(cells)

# Track completed reps per cell
done_per_cell <- integer(n_cells)
names(done_per_cell) <- paste(cells$family, cells$architecture,
                              sep = '|')

cell_results <- vector('list', n_cells)
names(cell_results) <- names(done_per_cell)

n_fits_done   <- 0L
last_progress <- 0L
abort         <- FALSE

# Loop over batches (round-robin across cells) until either
# reps_target is reached for all cells or deadline is hit.
batch_idx <- 0L
repeat {
  batch_idx <- batch_idx + 1L
  any_remaining <- FALSE
  for (ci in seq_len(n_cells)) {
    fam   <- cells$family[ci]
    arch  <- cells$architecture[ci]
    cellk <- paste(fam, arch, sep = '|')
    rem   <- reps_target - done_per_cell[ci]
    if (rem <= 0) next
    any_remaining <- TRUE
    take <- min(chunk_size, rem)
    seed_lo <- seed_offsets[[cellk]] + done_per_cell[ci] + 1L
    seeds   <- seed_lo:(seed_lo + take - 1L)

    t_now <- Sys.time()
    if (t_now >= deadline) {
      abort <- TRUE; break
    }
    # Estimated time for this batch given pilot timing and
    # parallelism.  If insufficient time even for one chunk,
    # shrink the chunk to fit.
    est_batch_s <- t_pilot * take / n_workers
    if ((t_now + est_batch_s) > deadline) {
      max_take <- max(1L,
        floor(as.numeric(deadline - t_now, units = 'secs') *
              n_workers / t_pilot))
      take <- min(take, max_take)
      seeds <- seed_lo:(seed_lo + take - 1L)
      if (take < 1L) { abort <- TRUE; break }
    }

    out <- mclapply(seeds, function(s) run_one(fam, arch, s),
                    mc.cores = n_workers, mc.preschedule = TRUE)
    out <- rbindlist(out)
    cell_results[[cellk]] <- rbind(cell_results[[cellk]], out)
    done_per_cell[ci] <- done_per_cell[ci] + take
    n_fits_done <- n_fits_done + take

    if ((n_fits_done - last_progress) >= 200L) {
      el <- as.numeric(Sys.time() - t_start, units = 'secs')
      log_msg(sprintf(
        'progress: %d fits done (cell %s @ %d, elapsed %.0f s)',
        n_fits_done, cellk, done_per_cell[ci], el))
      last_progress <- n_fits_done
    }

    if (Sys.time() >= deadline) { abort <- TRUE; break }
  }
  if (abort) break
  if (!any_remaining) break
}

t_end <- Sys.time()
elapsed <- as.numeric(t_end - t_start, units = 'secs')
log_msg(sprintf('Run complete. Elapsed %.1f s. Aborted=%s',
                elapsed, abort))
log_msg(sprintf('Reps achieved per cell:'))
for (cellk in names(done_per_cell)) {
  log_msg(sprintf('  %-45s : %d', cellk, done_per_cell[[cellk]]))
}

replicates <- rbindlist(cell_results)
saveRDS(replicates, file.path(out_data_dir, '07-gompertz-replicates.rds'))

# -- Morris ADEMP summary --------------------------------------------

alpha <- 0.05
summarise_cell <- function(df) {
  conv <- df[converged == TRUE]
  n_conv <- nrow(conv)
  n_total <- nrow(df)
  if (n_conv == 0L) {
    return(data.table(
      n_reps = n_total, power = NA_real_, mcse_power = NA_real_,
      mean_beta = NA_real_, sd_beta = NA_real_,
      mcse_mean_beta = NA_real_, conv_rate = 0
    ))
  }
  power <- mean(conv$p_bmDbc < alpha)
  mcse_power <- sqrt(power * (1 - power) / n_conv)
  mean_beta <- mean(conv$beta_bmDbc)
  sd_beta <- sd(conv$beta_bmDbc)
  mcse_mean_beta <- sd_beta / sqrt(n_conv)
  data.table(
    n_reps         = n_total,
    power          = round(power, 4),
    mcse_power     = round(mcse_power, 4),
    mean_beta      = round(mean_beta, 4),
    sd_beta        = round(sd_beta, 4),
    mcse_mean_beta = round(mcse_mean_beta, 4),
    conv_rate      = round(n_conv / n_total, 4)
  )
}

summary_dt <- replicates[, summarise_cell(.SD),
                         by = .(family, architecture)]
setcolorder(summary_dt, c('family', 'architecture', 'n_reps',
                          'power', 'mcse_power', 'mean_beta',
                          'sd_beta', 'mcse_mean_beta', 'conv_rate'))

# Write summary text and run log together
sink(file.path(out_data_dir, '07-gompertz-summary.txt'))
cat('# 07-gompertz-evaluation Morris summary\n')
cat(sprintf('# Generated %s\n',
            format(Sys.time(), '%Y-%m-%d %H:%M %Z')))
cat(sprintf('# Wall-clock elapsed: %.1f s (budget %d s, abort=%s)\n',
            elapsed, budget_s, abort))
cat(sprintf('# Total fits: %d  |  workers: %d\n',
            n_fits_done, n_workers))
cat('# Calibration anchor: week 8\n')
cat(sprintf(
  '#   logistic:           log_rate=%.3f t0=%.3f\n', log_rate, log_t0))
cat(sprintf(
  '#   hyperbolic_tangent: rate=%.3f t0=%.3f\n', htan_rate, htan_t0))
cat(sprintf(
  '#   piecewise_linear:   slope=%.4f t_bp=%.3f\n', pl_slope, pl_tbp))
cat('#\n')
cat('# Run log:\n')
for (line in run_log) cat('#   ', line, '\n', sep = '')
cat('#\n')
cat('# Cell summary (alpha = 0.05):\n')
print(summary_dt)
sink()

cat('\n--- Morris summary ---\n')
print(summary_dt)

# -- Power figure ----------------------------------------------------

summary_dt[, ci_lo := pmax(0, power - 1.96 * mcse_power)]
summary_dt[, ci_hi := pmin(1, power + 1.96 * mcse_power)]

# Order families for plotting clarity
fam_levels <- c('gompertz', 'logistic',
                'hyperbolic_tangent', 'piecewise_linear_breakpoint')
summary_dt[, family := factor(family, levels = fam_levels)]
summary_dt[, architecture := factor(architecture,
                                    levels = c('mvn', 'mean_moderation'),
                                    labels = c('B: MVN moderation',
                                               'A: mean moderation'))]

p <- ggplot(summary_dt,
            aes(x = family, y = power, fill = family)) +
  geom_col(width = 0.6, colour = 'black') +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_text(aes(label = sprintf('%.3f', power)),
            vjust = -0.6, size = 3) +
  facet_wrap(~ architecture, ncol = 2) +
  scale_fill_brewer(palette = 'Set2') +
  scale_y_continuous(limits = c(0, 1.08),
                     breaks = seq(0, 1, 0.2)) +
  labs(
    x = NULL,
    y = sprintf('Power at alpha = %.2f', alpha),
    title = 'Power for bm:Dbc by trajectory family x DGP architecture',
    subtitle = sprintf(
      'OL+BDC hybrid, N=%d, c.bm=%.2f, t1/2=0, anchored at week %.0f',
      mp$N, mp$c.bm, t_anchor)
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = 'none',
        plot.title = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(out_fig_dir, '07-gompertz-power.pdf'),
       p, width = 9, height = 4.6)

cat('\nWrote:\n  ',
    file.path(out_data_dir, '07-gompertz-replicates.rds'), '\n  ',
    file.path(out_data_dir, '07-gompertz-summary.txt'),    '\n  ',
    file.path(out_fig_dir,  '07-gompertz-power.pdf'),       '\n',
    sep = '')
