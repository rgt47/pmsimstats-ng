# 09-dropout-prototype.R
#
# Small-scale prototype simulation for paper
# 09-informative-dropout-by-design.
#
# Goal: contrast two trial designs (OL, OL+BDC) under two
# dropout patterns (none, symptom-worsening-driven hazard) at
# N=35 with c.bm=0.45 under Architecture A (mean moderation).
# 4 cells, 200 reps target, ~4 minute wall-clock budget.
#
# Outputs:
#   analysis/data/quick-sim/09-dropout-replicates.rds
#   analysis/data/quick-sim/09-dropout-summary.txt
#   analysis/figures/quick-sim/09-dropout-power.pdf

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
})

t_start <- Sys.time()
budget_secs <- 4 * 60     # hard abort at 4 minutes
probe_budget <- 25        # if first 25 reps in cell 1 exceed 25s, halve

set.seed(20260507)

# ---------------------------------------------------------------
# Trial designs
#
# OL: all 8 timepoints on drug. Because Db has no variation,
# lme_analysis falls back to Sx ~ bm + t + bm:t and returns the
# bm:t coefficient. We extract whichever interaction the
# fitter returned (bm:Dbc or bm:t) and mark the term used.
#
# OL+BDC: prompt's literal spec. 5 OL on-drug timepoints,
# then BDC1 on-drug, BDC2 and BDC3 off-drug. Db varies, so
# lme_analysis returns bm:Dbc.
# ---------------------------------------------------------------

td_OL <- buildtrialdesign(
  name_longform = 'open label',
  name_shortform = 'OL',
  timepoints = cumulative(rep(2.5, 8)),
  timeptnames = paste0('OL', 1:8),
  expectancies = rep(1, 8),
  ondrug = list(pathA = rep(1, 8))
)

td_OLBDC <- buildtrialdesign(
  name_longform = 'open label + blinded discontinuation',
  name_shortform = 'OLBDC',
  timepoints = cumulative(c(rep(2.5, 5), rep(1, 3))),
  timeptnames = c(paste0('OL', 1:5), paste0('BDC', 1:3)),
  expectancies = c(rep(1, 5), rep(0.5, 3)),
  ondrug = list(pathA = c(rep(1, 5), 1, 0, 0))
)

designs <- list(OL = td_OL, OLBDC = td_OLBDC)

# ---------------------------------------------------------------
# DGP parameters (matched to extracted defaults / typical
# pmsimstats prototype scale)
# ---------------------------------------------------------------

bp <- data.table(
  cat = c('bm', 'BL'),
  m = c(0, 70),
  sd = c(1, 10)
)

rp <- data.table(
  cat = c('tv', 'pb', 'br'),
  max = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd = c(5, 2, 5)
)

# ---------------------------------------------------------------
# Inline dropout module
#
# At each visit t > BL, the per-participant hazard is
#   h(t) = beta0 + beta1 * delta_sx(t)^2
# where delta_sx(t) = (Sx(t) - Sx(BL)) standardised so that
# squared values are positive and on a comparable scale across
# replicates. Once a participant drops, all subsequent
# measurements are NA. The intercept and slope are then rescaled
# globally so the achieved overall dropout rate by the last
# visit is approximately `target_rate`.
#
# This is the simple operationalisation requested by the prompt;
# it differs from R/censordata.R, which uses a multiplicative
# (beta0, beta1) parameterisation tied to the per-timepoint
# delta-from-baseline raised to a shape exponent. We honour the
# prompt's quadratic hazard formulation directly.
# ---------------------------------------------------------------

apply_hazard_dropout <- function(dat, trialdesign,
                                 beta0 = 0.05, beta1 = 0.5,
                                 target_rate = 0.25,
                                 scale_to_target = TRUE) {
  d <- data.table(trialdesign)
  tnames <- d$timeptnames
  n_pt <- nrow(dat)

  # delta_sx(t) = Sx(t) - BL; positive => worsening (since lower
  # post-baseline Sx means improvement when BL > Sx)
  # In pmsimstats, BL is the pre-treatment outcome, and Sx[t] is
  # generated as BL minus accumulated response components, so
  # smaller Sx means larger improvement. Worsening relative to
  # baseline corresponds to Sx[t] - BL > 0. For the hazard, we
  # take the squared standardised deviation.
  bl <- dat$BL
  sx_mat <- as.matrix(dat[, mget(tnames)])
  delta_mat <- sx_mat - bl
  # standardise per visit by SD across participants to keep the
  # squared term on a stable scale
  v_sd <- apply(delta_mat, 2, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s < 1e-6) 1 else s
  })
  delta_z <- sweep(delta_mat, 2, v_sd, '/')

  # raw hazard at each visit (n_pt x n_visit)
  haz_raw <- beta0 + beta1 * (delta_z^2)
  haz_raw[haz_raw < 0] <- 0
  haz_raw[haz_raw > 1] <- 1

  # If requested, find a global multiplier alpha so that the
  # expected proportion ever-dropped by last visit = target_rate.
  # We use a one-pass survival approximation:
  #   P(survive to end) ~= prod_t (1 - alpha * haz_raw)
  # Solve for alpha by univariate search.
  alpha <- 1
  if (scale_to_target) {
    obj <- function(a) {
      surv <- exp(rowSums(log(pmax(1 - a * haz_raw, 1e-9))))
      mean(1 - surv) - target_rate
    }
    a_lo <- 1e-6
    a_hi <- 1
    if (obj(a_lo) > 0) {
      alpha <- a_lo
    } else if (obj(a_hi) < 0) {
      alpha <- a_hi
    } else {
      alpha <- uniroot(obj, c(a_lo, a_hi), tol = 1e-4)$root
    }
  }
  haz <- pmin(pmax(alpha * haz_raw, 0), 1)

  # Bernoulli draws -> first-failure index per participant
  draws <- matrix(runif(n_pt * length(tnames)),
                  nrow = n_pt, ncol = length(tnames))
  drop_event <- draws < haz
  first_drop <- apply(drop_event, 1, function(x) {
    w <- which(x)
    if (length(w) == 0) NA_integer_ else w[1]
  })

  # Apply: NA out from first_drop onwards
  for (i in seq_len(n_pt)) {
    fi <- first_drop[i]
    if (!is.na(fi)) {
      for (j in fi:length(tnames)) {
        set(dat, i = i, j = tnames[j], value = NA_real_)
      }
    }
  }

  dropout_fraction <- mean(!is.na(first_drop))
  list(dat = dat, dropout_fraction = dropout_fraction,
       alpha = alpha)
}

# ---------------------------------------------------------------
# One replicate
# ---------------------------------------------------------------

run_one_rep <- function(design_name, dropout_pattern) {
  td <- designs[[design_name]]
  paths <- td$trialpaths
  c_bm <- 0.45
  t1half <- 0.5

  mp <- data.table(
    N = 35L, c.bm = c_bm,
    carryover_t1half = t1half,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.1, c.cfct = 0.05
  )

  dat_list <- vector('list', length(paths))
  drop_fracs <- numeric(length(paths))
  for (g in seq_along(paths)) {
    di <- generateData(
      modelparam = mp, respparam = rp, blparam = bp,
      trialdesign = paths[[g]],
      empirical = FALSE, makePositiveDefinite = TRUE,
      dgp_architecture = 'mean_moderation'
    )
    di[, path := g]

    if (dropout_pattern == 'biased') {
      out <- apply_hazard_dropout(
        di, paths[[g]],
        beta0 = 0.05, beta1 = 0.5,
        target_rate = 0.25,
        scale_to_target = TRUE
      )
      di <- out$dat
      drop_fracs[g] <- out$dropout_fraction
    } else {
      drop_fracs[g] <- 0
    }
    dat_list[[g]] <- di
  }
  dat <- rbindlist(dat_list, fill = TRUE)

  op <- list(
    useDE = TRUE,
    t_random_slope = FALSE,
    full_model_out = FALSE,
    carryover_t1half = t1half,
    simplecarryover = FALSE,
    carryover_scalefactor = 1
  )

  res <- tryCatch(
    lme_analysis(td$trialpaths, dat, op),
    error = function(e) {
      data.table(beta = NA_real_, betaSE = NA_real_,
                 p = NA_real_, issingular = NA,
                 warning = conditionMessage(e))
    }
  )
  res[, dropout_fraction := mean(drop_fracs)]
  res
}

# ---------------------------------------------------------------
# Cell grid
# ---------------------------------------------------------------

cells <- expand.grid(
  design = c('OL', 'OLBDC'),
  dropout_pattern = c('none', 'biased'),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
cat('Cells:', nrow(cells), '\n')
print(cells)

n_reps <- 200L

# ---------------------------------------------------------------
# Probe: time first 25 reps of cell 1 to see if we need to halve
# ---------------------------------------------------------------

probe_n <- 25L
probe_cell <- cells[1, ]
cat(sprintf('\nProbe: cell 1 (%s, dropout=%s), %d reps...\n',
            probe_cell$design, probe_cell$dropout_pattern,
            probe_n))
t_probe <- Sys.time()
probe_results <- vector('list', probe_n)
for (k in seq_len(probe_n)) {
  probe_results[[k]] <- run_one_rep(probe_cell$design,
                                    probe_cell$dropout_pattern)
}
probe_secs <- as.numeric(difftime(Sys.time(), t_probe,
                                  units = 'secs'))
cat(sprintf('Probe time: %.1f s for %d reps\n',
            probe_secs, probe_n))

# Extrapolate to total cost across 4 cells x n_reps. If the
# extrapolation exceeds 200s (matches the prompt's >100s/100reps
# rule of thumb), halve.
extrapolated <- probe_secs * (n_reps / probe_n) * nrow(cells)
cat(sprintf('Extrapolated total: %.0f s\n', extrapolated))
if (extrapolated > 200) {
  n_reps <- 100L
  cat(sprintf('Halving to n_reps=%d\n', n_reps))
} else {
  cat(sprintf('Within budget; keeping n_reps=%d\n', n_reps))
}

# ---------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------

all_rows <- list()
aborted <- FALSE

for (i in seq_len(nrow(cells))) {
  ci <- cells[i, ]
  cat(sprintf('\nCell %d/%d: design=%s dropout=%s\n',
              i, nrow(cells), ci$design, ci$dropout_pattern))
  cell_rows <- vector('list', n_reps)

  start_k <- 1L
  if (i == 1L) {
    n_reuse <- min(probe_n, n_reps)
    cell_rows[seq_len(n_reuse)] <- probe_results[seq_len(n_reuse)]
    start_k <- n_reuse + 1L
  }

  if (start_k <= n_reps) {
    for (k in seq.int(start_k, n_reps)) {
      cell_rows[[k]] <- run_one_rep(ci$design,
                                    ci$dropout_pattern)
      if (k %% 50 == 0) {
        elapsed <- as.numeric(difftime(Sys.time(), t_start,
                                       units = 'secs'))
        cat(sprintf('  rep %d, total elapsed %.1fs\n',
                    k, elapsed))
        if (elapsed > budget_secs) {
          cat('  Budget exceeded; aborting after this cell.\n')
          aborted <- TRUE
          cell_rows <- cell_rows[seq_len(k)]
          break
        }
      }
    }
  }

  cell_dt <- rbindlist(cell_rows, fill = TRUE)
  cell_dt[, `:=`(design = ci$design,
                 dropout_pattern = ci$dropout_pattern,
                 rep_idx = seq_len(.N))]
  all_rows[[i]] <- cell_dt

  if (aborted) break
}

reps_dt <- rbindlist(all_rows, fill = TRUE)
setnames(reps_dt,
         old = c('beta', 'betaSE', 'p'),
         new = c('beta_bmDbc', 'betaSE_bmDbc', 'p_bmDbc'))
setcolorder(reps_dt,
            c('design', 'dropout_pattern', 'rep_idx',
              'beta_bmDbc', 'betaSE_bmDbc', 'p_bmDbc',
              'dropout_fraction'))
cat(sprintf('\nTotal reps collected: %d across %d cells\n',
            nrow(reps_dt),
            uniqueN(reps_dt[, .(design, dropout_pattern)])))

t_total <- as.numeric(difftime(Sys.time(), t_start,
                               units = 'secs'))
cat(sprintf('Total wall-clock: %.1f s\n', t_total))

# ---------------------------------------------------------------
# Persist replicates
# ---------------------------------------------------------------

dir.create('analysis/data/quick-sim', showWarnings = FALSE,
           recursive = TRUE)
dir.create('analysis/figures/quick-sim', showWarnings = FALSE,
           recursive = TRUE)

saveRDS(reps_dt, 'analysis/data/quick-sim/09-dropout-replicates.rds')

# ---------------------------------------------------------------
# Morris-style summary table with MCSE columns
# ---------------------------------------------------------------

summary_dt <- reps_dt[, .(
  n_reps = .N,
  n_converged = sum(!is.na(p_bmDbc)),
  conv_rate = mean(!is.na(p_bmDbc)),
  power = mean(p_bmDbc < 0.05, na.rm = TRUE),
  mean_dropout_fraction = mean(dropout_fraction, na.rm = TRUE),
  mean_beta = mean(beta_bmDbc, na.rm = TRUE),
  sd_beta = sd(beta_bmDbc, na.rm = TRUE)
), by = .(design, dropout_pattern)]

summary_dt[, mcse_power := sqrt(power * (1 - power) /
                                pmax(n_converged, 1))]
summary_dt[, mcse_mean_beta := sd_beta /
                               sqrt(pmax(n_converged, 1))]
setcolorder(summary_dt,
            c('design', 'dropout_pattern', 'n_reps',
              'power', 'mcse_power',
              'mean_dropout_fraction',
              'mean_beta', 'sd_beta', 'mcse_mean_beta',
              'conv_rate'))
setorder(summary_dt, design, dropout_pattern)

write.table(summary_dt,
            file = 'analysis/data/quick-sim/09-dropout-summary.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)

cat('\n--- Summary ---\n')
print(summary_dt)

# ---------------------------------------------------------------
# Power figure
# ---------------------------------------------------------------

plot_dt <- copy(summary_dt)
plot_dt[, ci_lo := pmax(0, power - 1.96 * mcse_power)]
plot_dt[, ci_hi := pmin(1, power + 1.96 * mcse_power)]
plot_dt[, design_lab := factor(design, levels = c('OL', 'OLBDC'),
                               labels = c('OL', 'OL+BDC'))]
plot_dt[, drop_lab := factor(
  dropout_pattern,
  levels = c('none', 'biased'),
  labels = c('No dropout', 'Symptom-worsening hazard'))]

p_power <- ggplot(plot_dt,
                  aes(x = design_lab, y = power,
                      fill = design_lab)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.2) +
  geom_hline(yintercept = 0.05, linetype = 'dashed',
             alpha = 0.4) +
  facet_wrap(~ drop_lab) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = 'Design', y = 'Power (P(p < 0.05))',
       fill = 'Design',
       title = paste0('Power for biomarker x treatment ',
                      'interaction by design and dropout'),
       subtitle = sprintf(
         paste0('Architecture A, N=35, c.bm=0.45, ',
                '%d reps per cell'), n_reps),
       caption = paste0('Error bars: 1.96 x binomial MCSE. ',
                        'Dashed line: nominal alpha = 0.05.')) +
  theme_minimal(base_size = 11) +
  theme(legend.position = 'none')

ggsave('analysis/figures/quick-sim/09-dropout-power.pdf',
       p_power, width = 8, height = 4.5)

cat('\nWrote:\n')
cat('  analysis/data/quick-sim/09-dropout-replicates.rds\n')
cat('  analysis/data/quick-sim/09-dropout-summary.txt\n')
cat('  analysis/figures/quick-sim/09-dropout-power.pdf\n')

if (aborted) {
  cat('\nNOTE: budget aborted; cells may have fewer reps.\n')
}
