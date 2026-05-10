# 01-architecture-a-trajectory-sweep.R
#
# Tier-F item 27, paper 07 (gompertz-evaluation): Architecture A
# trajectory comparators. Replaces the earlier 200-rep, two-family
# Architecture-B-only prototype with a 16-cell sweep that crosses
# DGP architecture x trajectory family x biomarker effect size, with
# 1000 replicates per cell (MCSE ~ 0.016 at p=0.5).
#
# Cells (16 total):
#   architecture in {mvn, mean_moderation}                           (2)
#   family       in {gompertz, logistic, hyperbolic_tangent,
#                    piecewise_linear_breakpoint}                    (4)
#   c.bm         in {0 (null), 0.45 (alternative)}                   (2)
#
# Trajectory families. The package exports four 3-parameter monotone-
# saturating trajectories via trajectoryShape():
#
#   gompertz                    p2=disp,         p3=rate
#   logistic                    p2=log_rate,     p3=t0
#   hyperbolic_tangent          p2=rate,         p3=t0
#   piecewise_linear_breakpoint p2=t_breakpoint, p3=post_slope
#
# Cross-family calibration. All four families share the asymptote
# maxr from extracted_rp[cat=='br'] and are calibrated to produce a
# common value at the anchor week t_anchor=5. The non-Gompertz
# families fix one shape parameter at a default and solve the other
# analytically:
#
#   logistic:           log_rate = 0.5, t0 from anchor closed form
#   hyperbolic_tangent: rate     = 0.5, t0 from anchor closed form
#   piecewise_linear:   post_slope = 0, t_breakpoint = anchor *
#                                                     maxr / anchor_y
#
# Family is varied via post-hoc additive shifts on the BR columns of
# the simulated data. The shift is family(tod) - gompertz(tod) at each
# timepoint, applied after the MVN draw (and after, if applicable,
# the Architecture-A mean-moderation shift). This decouples family
# from architecture and from the biomarker correlation block, both of
# which depend on the package's gompertz-derived BR mean.
#
# Parallelism. parallel::mclapply with mc.cores = 8 (HARD requirement,
# fixed; not auto-detected). Time budget 1400 s with graceful abort.
#
# Author: pmsimstats team
# Last updated: 2026-05-09

suppressMessages({
  devtools::load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
  library(parallel)
})

set.seed(20260509)

t_master_start <- Sys.time()
master_budget_s <- 1400

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

op <- list(useDE = FALSE, t_random_slope = FALSE,
           full_model_out = FALSE, carryover_t1half = 0,
           simplecarryover = FALSE, carryover_scalefactor = 1)

# -- Trajectory family calibration -----------------------------------
# Anchor at t=5 weeks. modGompertz BR(5) defines the common target
# value the other three families are calibrated to match. All four
# families share maxr.

br_pars <- as.list(rp[cat == 'br'])
maxr_v  <- br_pars$max
g_disp  <- br_pars$disp
g_rate  <- br_pars$rate

t_anchor <- 5
g_anchor <- modgompertz(t_anchor, maxr_v, g_disp, g_rate)

# Logistic: shape_logistic uses vertical-offset adjustment. Solve t0
# given log_rate so that shape_logistic(t_anchor) = g_anchor.
# Closed form below uses Newton search since the vertical-offset
# rescale couples both ends of the curve.
# The root-finding interval is restricted to t0 in [-5, 20] so that
# the vertical-offset rescale stays finite (t0 -> -infty drives the
# offset to maxr and the rescale denominator to zero).
solve_logistic_t0 <- function(target, anchor, log_rate, maxr) {
  f <- function(t0) {
    shape_logistic(anchor, maxr, log_rate, t0) - target
  }
  uniroot(f, interval = c(-5, 20), tol = 1e-10)$root
}
log_rate <- 0.5
log_t0   <- solve_logistic_t0(g_anchor, t_anchor, log_rate, maxr_v)

# Hyperbolic tangent: same vertical-offset convention. Solve t0 given
# rate.
solve_htan_t0 <- function(target, anchor, rate, maxr) {
  f <- function(t0) {
    shape_hyperbolic_tangent(anchor, maxr, rate, t0) - target
  }
  uniroot(f, interval = c(-5, 20), tol = 1e-10)$root
}
htan_rate <- 0.5
htan_t0   <- solve_htan_t0(g_anchor, t_anchor, htan_rate, maxr_v)

# Piecewise-linear breakpoint: y(t) = maxr*min(1,t/t_bp) (post_slope=0)
# implies y(t_anchor)=g_anchor when t_bp = t_anchor * maxr / g_anchor.
pl_post_slope <- 0
pl_tbp        <- t_anchor * maxr_v / g_anchor

family_param <- list(
  gompertz = list(p2 = g_disp,    p3 = g_rate),
  logistic = list(p2 = log_rate,  p3 = log_t0),
  hyperbolic_tangent = list(p2 = htan_rate, p3 = htan_t0),
  piecewise_linear_breakpoint = list(p2 = pl_tbp, p3 = pl_post_slope)
)

# Verify the calibration. Anchor must match across families to <1e-4.
# (uniroot's analytic-equivalent tolerance after the offset rescale is
# of order tol; we keep the assertion at 1e-4 to allow for the
# composition of root finder and rescale numerical error.)
for (fam in names(family_param)) {
  pp <- family_param[[fam]]
  v <- trajectoryShape(fam, t_anchor, maxr_v, pp$p2, pp$p3)
  stopifnot(abs(v - g_anchor) < 1e-4)
}

cat('-- Trajectory family calibration --\n')
cat(sprintf('  Anchor week %.1f, gompertz BR(anchor) = %.4f, maxr = %.4f\n',
            t_anchor, g_anchor, maxr_v))
cat(sprintf('  logistic:           log_rate = %.3f, t0 = %.3f\n',
            log_rate, log_t0))
cat(sprintf('  hyperbolic_tangent: rate = %.3f, t0 = %.3f\n',
            htan_rate, htan_t0))
cat(sprintf('  piecewise_linear:   t_bp = %.3f, post_slope = %.3f\n',
            pl_tbp, pl_post_slope))

# Per-path BR mean shift relative to gompertz, for each family.
# Family 'gompertz' has shift identically zero by construction.
families <- names(family_param)
br_shift_by_family_path <- list()
for (fam in families) {
  pp <- family_param[[fam]]
  br_shift_by_family_path[[fam]] <- vector('list', n_paths)
  for (g in seq_len(n_paths)) {
    tp_path <- td$trialpaths[[g]]
    tod_g   <- tp_path$tod
    fam_y   <- trajectoryShape(fam, tod_g, maxr_v, pp$p2, pp$p3)
    gomp_y  <- modgompertz(tod_g, maxr_v, g_disp, g_rate)
    shift_g <- fam_y - gomp_y
    names(shift_g) <- tp_path$timeptnames
    br_shift_by_family_path[[fam]][[g]] <- shift_g
  }
}

# -- Sigma cache: one per (architecture, c.bm, path) -----------------
# c.bm enters the BM-BR cross-correlation under Architecture B and
# the post-hoc mean shift under Architecture A. Caching at the
# (arch, c.bm, path) granularity avoids per-replicate sigma rebuilds.

architectures <- c('mvn', 'mean_moderation')
cbm_levels    <- c(0, 0.45)
N_subjects    <- 35

cs_cache <- list()
for (arch in architectures) {
  for (cv in cbm_levels) {
    key <- paste(arch, cv, sep = '|')
    cs_cache[[key]] <- vector('list', n_paths)
    mp_local <- data.table(
      N = N_subjects, c.bm = cv,
      carryover_t1half = 0,
      c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
      c.cf1t = 0.2, c.cfct = 0.1
    )
    for (g in seq_len(n_paths)) {
      cs_cache[[key]][[g]] <- buildSigma(
        mp_local, rp, bp, td$trialpaths[[g]],
        makePositiveDefinite = TRUE,
        dgp_architecture = arch
      )
    }
  }
}

# -- Single-replicate driver -----------------------------------------

run_one <- function(family, architecture, cbm, seed) {
  set.seed(seed)
  shifts <- br_shift_by_family_path[[family]]
  cs_set <- cs_cache[[paste(architecture, cbm, sep = '|')]]

  mp_local <- data.table(
    N = N_subjects, c.bm = cbm,
    carryover_t1half = 0,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.2, c.cfct = 0.1
  )

  dats <- vector('list', n_paths)
  for (g in seq_len(n_paths)) {
    tp_path <- td$trialpaths[[g]]
    dat <- generateData(mp_local, rp, bp, tp_path,
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
    c.bm         = cbm,
    rep_idx      = seed,
    beta_bmDbc   = res$beta,
    betaSE_bmDbc = res$betaSE,
    p_bmDbc      = res$p,
    issingular   = ifelse(is.na(res$issingular), NA, res$issingular),
    converged    = !is.na(res$beta)
  )
}

# -- Run plan: 16 cells x 1000 reps ---------------------------------

reps_target <- 1000L
cells <- expand.grid(
  family = families,
  architecture = architectures,
  c.bm = cbm_levels,
  stringsAsFactors = FALSE
)
n_cells <- nrow(cells)

# Distinct seed offsets per cell so that replicates are independent
# across cells but reproducible.
seed_base <- 800000L
cells$cellk <- with(cells, paste(family, architecture, c.bm, sep = '|'))
seed_offsets <- setNames(
  seed_base + (seq_len(n_cells) - 1L) * 200000L,
  cells$cellk
)

t_start    <- Sys.time()
budget_s   <- 1400 - as.numeric(t_start - t_master_start, units = 'secs')
deadline   <- t_master_start + 1400
n_workers  <- 8L

cat(sprintf('\n-- Run plan --\n'))
cat(sprintf('  Cells: %d (4 families x 2 arch x 2 c.bm)\n', n_cells))
cat(sprintf('  Target reps/cell: %d  -> total fits target: %d\n',
            reps_target, reps_target * n_cells))
cat(sprintf('  Wall-clock budget: %.0f s (master deadline)\n',
            budget_s))
cat(sprintf('  Forked workers: %d\n', n_workers))

run_log <- character()
log_msg <- function(msg) {
  ts <- format(Sys.time(), '%H:%M:%S')
  line <- sprintf('[%s] %s', ts, msg)
  run_log[[length(run_log) + 1L]] <<- line
  cat(line, '\n')
}

log_msg(sprintf('Start at %s, master deadline %s',
                format(t_start, '%H:%M:%S'),
                format(deadline, '%H:%M:%S')))

# Pilot timing on a few sequential fits so we can estimate total time.
log_msg('Pilot timing (8 fits, sequential, gompertz/mvn/c.bm=0.45) ...')
t_pilot_start <- Sys.time()
pilot <- rbindlist(lapply(seq_len(8), function(i)
  run_one('gompertz', 'mvn', 0.45, seed_base + i)))
t_pilot <- as.numeric(Sys.time() - t_pilot_start, units = 'secs') / 8
log_msg(sprintf(paste0('Per-fit (sequential): %.3f s; ',
                       '8-core projection for %d fits: %.0f s'),
                t_pilot, reps_target * n_cells,
                t_pilot * reps_target * n_cells / n_workers))

# Round-robin chunked dispatch with per-batch deadline check.
chunk_size  <- 100L
done_per_cell <- integer(n_cells)
names(done_per_cell) <- cells$cellk
cell_results <- vector('list', n_cells)
names(cell_results) <- cells$cellk

n_fits_done   <- 0L
last_progress <- 0L
abort         <- FALSE
batch_idx     <- 0L

repeat {
  batch_idx <- batch_idx + 1L
  any_remaining <- FALSE
  for (ci in seq_len(n_cells)) {
    fam   <- cells$family[ci]
    arch  <- cells$architecture[ci]
    cbm   <- cells$c.bm[ci]
    cellk <- cells$cellk[ci]
    rem   <- reps_target - done_per_cell[ci]
    if (rem <= 0) next
    any_remaining <- TRUE
    take <- min(chunk_size, rem)
    seed_lo <- seed_offsets[[cellk]] + done_per_cell[ci] + 1L
    seeds   <- seed_lo:(seed_lo + take - 1L)

    t_now <- Sys.time()
    if (t_now >= deadline) { abort <- TRUE; break }
    est_batch_s <- t_pilot * take / n_workers
    if ((t_now + est_batch_s) > deadline) {
      max_take <- max(1L,
        floor(as.numeric(deadline - t_now, units = 'secs') *
              n_workers / t_pilot))
      take <- min(take, max_take)
      seeds <- seed_lo:(seed_lo + take - 1L)
      if (take < 1L) { abort <- TRUE; break }
    }

    out <- mclapply(seeds,
                    function(s) run_one(fam, arch, cbm, s),
                    mc.cores = n_workers, mc.preschedule = TRUE)
    out <- rbindlist(out)
    cell_results[[cellk]] <- rbind(cell_results[[cellk]], out)
    done_per_cell[ci] <- done_per_cell[ci] + take
    n_fits_done <- n_fits_done + take

    if ((n_fits_done - last_progress) >= 1000L) {
      el <- as.numeric(Sys.time() - t_master_start, units = 'secs')
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
elapsed <- as.numeric(t_end - t_master_start, units = 'secs')
log_msg(sprintf('Run complete. Elapsed %.1f s. Aborted=%s',
                elapsed, abort))
log_msg(sprintf('Reps achieved per cell:'))
for (cellk in names(done_per_cell)) {
  log_msg(sprintf('  %-55s : %d', cellk, done_per_cell[[cellk]]))
}

replicates <- rbindlist(cell_results)
saveRDS(replicates,
        file.path(out_data_dir, '07-gompertz-replicates.rds'))

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
                         by = .(family, architecture, c.bm)]
setcolorder(summary_dt,
            c('family', 'architecture', 'c.bm', 'n_reps',
              'power', 'mcse_power', 'mean_beta',
              'sd_beta', 'mcse_mean_beta', 'conv_rate'))

sink(file.path(out_data_dir, '07-gompertz-summary.txt'))
cat('# 07-gompertz-evaluation Morris ADEMP summary\n')
cat(sprintf('# Generated %s\n',
            format(Sys.time(), '%Y-%m-%d %H:%M %Z')))
cat(sprintf(paste0('# Wall-clock elapsed: %.1f s ',
                   '(master budget %d s, abort=%s)\n'),
            elapsed, master_budget_s, abort))
cat(sprintf('# Total fits: %d  |  workers: %d\n',
            n_fits_done, n_workers))
cat(sprintf('# Trial design: OL+BDC, N=%d, t1/2=0\n', N_subjects))
cat(sprintf(paste0('# Calibration anchor: week %.1f, ',
                   'gompertz BR(anchor)=%.4f\n'),
            t_anchor, g_anchor))
cat(sprintf('#   logistic:           log_rate=%.3f t0=%.3f\n',
            log_rate, log_t0))
cat(sprintf('#   hyperbolic_tangent: rate=%.3f t0=%.3f\n',
            htan_rate, htan_t0))
cat(sprintf('#   piecewise_linear:   t_bp=%.3f post_slope=%.3f\n',
            pl_tbp, pl_post_slope))
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

fam_levels <- c('gompertz', 'logistic',
                'hyperbolic_tangent', 'piecewise_linear_breakpoint')
summary_dt[, family := factor(family, levels = fam_levels)]
summary_dt[, architecture := factor(architecture,
                                    levels = c('mvn', 'mean_moderation'),
                                    labels = c('B: MVN moderation',
                                               'A: mean moderation'))]
summary_dt[, c.bm_lab := factor(c.bm,
                                levels = c(0, 0.45),
                                labels = c('c.bm = 0 (null)',
                                           'c.bm = 0.45 (alt)'))]

p <- ggplot(summary_dt,
            aes(x = family, y = power, fill = family)) +
  geom_col(width = 0.6, colour = 'black') +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
  geom_text(aes(label = sprintf('%.3f', power)),
            vjust = -0.6, size = 2.7) +
  facet_grid(c.bm_lab ~ architecture) +
  scale_fill_brewer(palette = 'Set2') +
  scale_y_continuous(limits = c(0, 1.08),
                     breaks = seq(0, 1, 0.2)) +
  labs(
    x = NULL,
    y = sprintf('Rejection rate at alpha = %.2f', alpha),
    title = paste('Power and Type I error for bm:Dbc by',
                  'trajectory family,'),
    subtitle = sprintf(
      'OL+BDC hybrid, N=%d, t1/2=0, anchored at week %.0f',
      N_subjects, t_anchor)
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = 'none',
        plot.title = element_text(face = 'bold'),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(out_fig_dir, '07-gompertz-power.pdf'),
       p, width = 11, height = 7.0)

cat('\nWrote:\n  ',
    file.path(out_data_dir, '07-gompertz-replicates.rds'), '\n  ',
    file.path(out_data_dir, '07-gompertz-summary.txt'),    '\n  ',
    file.path(out_fig_dir,  '07-gompertz-power.pdf'),       '\n',
    sep = '')
