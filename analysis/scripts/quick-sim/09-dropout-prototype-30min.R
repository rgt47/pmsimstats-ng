## 09-dropout-prototype-30min.R
##
## 30-minute prototype simulation for paper
## 09-informative-dropout-by-design.
##
## Grid: design {OL, OLBDC, CO, Hybrid}
##       x dropout_pattern {none, MCAR_25, biased_25, biased_40}
##       = 16 cells.
## Target: 800 reps per cell (12,800 fits) within 1700 s wall clock,
## with graceful abort. Architecture A (mean moderation), N = 35,
## c.bm = 0.45, carryover_t1half = 0.5.
##
## Outputs (overwrite):
##   analysis/data/quick-sim/09-dropout-replicates.rds
##   analysis/data/quick-sim/09-dropout-summary.txt
##   analysis/figures/quick-sim/09-dropout-power.pdf

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
})

t_start <- Sys.time()
budget_secs <- 1700              # graceful abort target

set.seed(20260507)

## -----------------------------------------------------------------
## Trial designs (canonical Hendrickson presets reproduced from
## analysis/scripts/carryover-sensitivity/simulation-core.R)
## -----------------------------------------------------------------

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
  timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptnames = c('OL1', 'OL2', 'OL3', 'OL4',
                  'BD1', 'BD2', 'BD3', 'BD4'),
  expectancies = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)

td_CO <- buildtrialdesign(
  name_longform = 'traditional crossover',
  name_shortform = 'CO',
  timepoints = cumulative(rep(2.5, 8)),
  timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
  expectancies = rep(0.5, 8),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
)

td_Hybrid <- buildtrialdesign(
  name_longform = 'hybrid N-of-1',
  name_shortform = 'Hybrid',
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptnames = c('OL1', 'OL2', 'BD1', 'BD2',
                  'BD3', 'BD4', 'COd', 'COp'),
  expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
  )
)

designs <- list(OL = td_OL, OLBDC = td_OLBDC,
                CO = td_CO, Hybrid = td_Hybrid)

## -----------------------------------------------------------------
## DGP parameters (matched to prior 5-min driver)
## -----------------------------------------------------------------

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

## -----------------------------------------------------------------
## Dropout module
##
## Three flavours:
##   none      -- no dropout
##   MCAR_<r>  -- per-visit hazard independent of outcome trajectory,
##                rescaled so mean ever-dropout by last visit ~= r
##   biased_<r>-- per-visit hazard h(t) = beta0 + beta1 * delta_sx(t)^2,
##                with delta_sx(t) standardised across participants;
##                rescaled so mean ever-dropout by last visit ~= r
##
## In all cases, once a participant drops, all subsequent
## measurements are NA (last-observation-carried-forward into NA).
## -----------------------------------------------------------------

apply_dropout <- function(dat, trialdesign,
                          mode = c('none', 'MCAR', 'biased'),
                          target_rate = 0.25,
                          beta0 = 0.05, beta1 = 0.5) {
  mode <- match.arg(mode)
  if (mode == 'none' || target_rate <= 0) {
    return(list(dat = dat, dropout_fraction = 0, alpha = 0))
  }

  d <- data.table(trialdesign)
  tnames <- d$timeptnames
  n_pt <- nrow(dat)
  n_v <- length(tnames)

  if (mode == 'MCAR') {
    haz_raw <- matrix(1, nrow = n_pt, ncol = n_v)
  } else {
    bl <- dat$BL
    sx_mat <- as.matrix(dat[, mget(tnames)])
    delta_mat <- sx_mat - bl
    v_sd <- apply(delta_mat, 2, function(x) {
      s <- sd(x, na.rm = TRUE)
      if (!is.finite(s) || s < 1e-6) 1 else s
    })
    delta_z <- sweep(delta_mat, 2, v_sd, '/')
    haz_raw <- beta0 + beta1 * (delta_z^2)
    haz_raw[!is.finite(haz_raw)] <- beta0
    haz_raw[haz_raw < 0] <- 0
    haz_raw[haz_raw > 1] <- 1
  }

  ## rescale so that expected mean(1 - prod(1 - alpha*h)) = target_rate
  obj <- function(a) {
    surv <- exp(rowSums(log(pmax(1 - a * haz_raw, 1e-9))))
    mean(1 - surv) - target_rate
  }
  a_lo <- 1e-6
  a_hi <- 1
  alpha <- if (obj(a_lo) > 0) {
    a_lo
  } else if (obj(a_hi) < 0) {
    a_hi
  } else {
    uniroot(obj, c(a_lo, a_hi), tol = 1e-4)$root
  }
  haz <- pmin(pmax(alpha * haz_raw, 0), 1)

  draws <- matrix(runif(n_pt * n_v), nrow = n_pt, ncol = n_v)
  drop_event <- draws < haz
  first_drop <- apply(drop_event, 1, function(x) {
    w <- which(x)
    if (length(w) == 0) NA_integer_ else w[1]
  })

  for (i in seq_len(n_pt)) {
    fi <- first_drop[i]
    if (!is.na(fi)) {
      for (j in fi:n_v) {
        set(dat, i = i, j = tnames[j], value = NA_real_)
      }
    }
  }

  list(dat = dat,
       dropout_fraction = mean(!is.na(first_drop)),
       alpha = alpha)
}

## -----------------------------------------------------------------
## One replicate
## -----------------------------------------------------------------

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

  ## decode dropout_pattern into (mode, target_rate)
  spec <- switch(
    dropout_pattern,
    none      = list(mode = 'none',   rate = 0),
    MCAR_25   = list(mode = 'MCAR',   rate = 0.25),
    biased_25 = list(mode = 'biased', rate = 0.25),
    biased_40 = list(mode = 'biased', rate = 0.40),
    stop('Unknown dropout pattern: ', dropout_pattern)
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

    out <- apply_dropout(di, paths[[g]],
                         mode = spec$mode,
                         target_rate = spec$rate)
    dat_list[[g]] <- out$dat
    drop_fracs[g] <- out$dropout_fraction
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

## -----------------------------------------------------------------
## Cell grid
## -----------------------------------------------------------------

cells <- expand.grid(
  design = c('OL', 'OLBDC', 'CO', 'Hybrid'),
  dropout_pattern = c('none', 'MCAR_25', 'biased_25', 'biased_40'),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
n_cells <- nrow(cells)
cat('Cells:', n_cells, '\n')
print(cells)

n_reps_target <- 800L

## -----------------------------------------------------------------
## Main loop with graceful abort
##
## We interleave reps across cells (round-robin) so that, if the
## budget is exhausted, every cell has a comparable number of reps
## rather than the early cells being saturated and the late cells
## empty.
## -----------------------------------------------------------------

cell_rows <- vector('list', n_cells)
for (i in seq_len(n_cells)) {
  cell_rows[[i]] <- vector('list', n_reps_target)
}
cell_counts <- integer(n_cells)

aborted <- FALSE
total_fits <- 0L
last_progress <- 0L

for (k in seq_len(n_reps_target)) {
  for (i in seq_len(n_cells)) {
    ci <- cells[i, ]
    cell_rows[[i]][[k]] <- run_one_rep(ci$design,
                                       ci$dropout_pattern)
    cell_counts[i] <- k
    total_fits <- total_fits + 1L

    if (total_fits - last_progress >= 200L) {
      elapsed <- as.numeric(difftime(Sys.time(), t_start,
                                     units = 'secs'))
      cat(sprintf('  fits %d, rep round %d/%d, elapsed %.1fs\n',
                  total_fits, k, n_reps_target, elapsed))
      last_progress <- total_fits
      if (elapsed > budget_secs) {
        cat('  Budget exceeded; stopping after this fit.\n')
        aborted <- TRUE
        break
      }
    }
  }
  if (aborted) break

  ## also check budget at end of each round
  elapsed <- as.numeric(difftime(Sys.time(), t_start,
                                 units = 'secs'))
  if (elapsed > budget_secs) {
    cat(sprintf('  Round %d complete; budget reached at %.1fs.\n',
                k, elapsed))
    aborted <- TRUE
    break
  }
}

t_total <- as.numeric(difftime(Sys.time(), t_start,
                               units = 'secs'))
cat(sprintf('\nTotal fits: %d   Wall clock: %.1f s\n',
            total_fits, t_total))
cat('Per-cell rep counts:\n')
print(setNames(cell_counts, paste(cells$design, cells$dropout_pattern,
                                  sep = ':')))

## -----------------------------------------------------------------
## Assemble replicates table
## -----------------------------------------------------------------

all_rows <- vector('list', n_cells)
for (i in seq_len(n_cells)) {
  ci <- cells[i, ]
  rows_i <- cell_rows[[i]]
  rows_i <- rows_i[!vapply(rows_i, is.null, logical(1))]
  if (length(rows_i) == 0L) next
  cell_dt <- rbindlist(rows_i, fill = TRUE)
  cell_dt[, `:=`(design = ci$design,
                 dropout_pattern = ci$dropout_pattern,
                 rep_idx = seq_len(.N))]
  all_rows[[i]] <- cell_dt
}

reps_dt <- rbindlist(all_rows, fill = TRUE)
setnames(reps_dt,
         old = c('beta', 'betaSE', 'p'),
         new = c('beta_bmDbc', 'betaSE_bmDbc', 'p_bmDbc'))
reps_dt[, converged := !is.na(p_bmDbc)]
setcolorder(reps_dt,
            c('design', 'dropout_pattern', 'rep_idx',
              'beta_bmDbc', 'betaSE_bmDbc', 'p_bmDbc',
              'converged', 'dropout_fraction'))

cat(sprintf('\nReplicates table: %d rows across %d cells\n',
            nrow(reps_dt),
            uniqueN(reps_dt[, .(design, dropout_pattern)])))

## -----------------------------------------------------------------
## Persist replicates
## -----------------------------------------------------------------

dir.create('analysis/data/quick-sim', showWarnings = FALSE,
           recursive = TRUE)
dir.create('analysis/figures/quick-sim', showWarnings = FALSE,
           recursive = TRUE)

saveRDS(reps_dt,
        'analysis/data/quick-sim/09-dropout-replicates.rds')

## -----------------------------------------------------------------
## Morris-style summary table with MCSE columns
## -----------------------------------------------------------------

summary_dt <- reps_dt[, .(
  n_reps = .N,
  n_converged = sum(converged),
  conv_rate = mean(converged),
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
            c('design', 'dropout_pattern',
              'n_reps', 'n_converged', 'conv_rate',
              'power', 'mcse_power',
              'mean_dropout_fraction',
              'mean_beta', 'sd_beta', 'mcse_mean_beta'))
setorder(summary_dt, design, dropout_pattern)

write.table(summary_dt,
            file = 'analysis/data/quick-sim/09-dropout-summary.txt',
            sep = '\t', quote = FALSE, row.names = FALSE)

cat('\n--- Morris-style summary ---\n')
print(summary_dt)

## -----------------------------------------------------------------
## Power figure: power across dropout_pattern, faceted by design
## -----------------------------------------------------------------

plot_dt <- copy(summary_dt)
plot_dt[, ci_lo := pmax(0, power - 1.96 * mcse_power)]
plot_dt[, ci_hi := pmin(1, power + 1.96 * mcse_power)]
plot_dt[, design_lab := factor(
  design,
  levels = c('OL', 'OLBDC', 'CO', 'Hybrid'),
  labels = c('OL', 'OL+BDC', 'Crossover', 'Hybrid'))]
plot_dt[, drop_lab := factor(
  dropout_pattern,
  levels = c('none', 'MCAR_25', 'biased_25', 'biased_40'),
  labels = c('No dropout', 'MCAR 25%',
             'Biased 25%', 'Biased 40%'))]

p_power <- ggplot(plot_dt,
                  aes(x = drop_lab, y = power, fill = drop_lab)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25) +
  geom_hline(yintercept = 0.05, linetype = 'dashed', alpha = 0.4) +
  facet_wrap(~ design_lab, nrow = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = 'Dropout pattern', y = 'Power (P(p < 0.05))',
       fill = 'Dropout pattern',
       title = paste0('Power for biomarker x treatment ',
                      'interaction by design and dropout'),
       subtitle = sprintf(
         paste0('Architecture A, N=35, c.bm=0.45, ',
                'target %d reps per cell'), n_reps_target),
       caption = paste0('Error bars: 1.96 x binomial MCSE. ',
                        'Dashed line: nominal alpha = 0.05.')) +
  theme_minimal(base_size = 11) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave('analysis/figures/quick-sim/09-dropout-power.pdf',
       p_power, width = 11, height = 4.5)

cat('\nWrote:\n')
cat('  analysis/data/quick-sim/09-dropout-replicates.rds\n')
cat('  analysis/data/quick-sim/09-dropout-summary.txt\n')
cat('  analysis/figures/quick-sim/09-dropout-power.pdf\n')

if (aborted) {
  cat(sprintf('\nNOTE: budget aborted at %.1fs; cells may have ',
              t_total))
  cat('fewer than target reps.\n')
} else {
  cat('\nCompleted full target rep count within budget.\n')
}
