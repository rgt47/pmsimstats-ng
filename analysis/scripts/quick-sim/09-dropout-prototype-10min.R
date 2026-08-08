## 09-dropout-prototype-10min.R
##
## 10-minute prototype simulation for paper
## 09-informative-dropout-by-design.
##
## Grid: design          {OL, OLBDC, traditional_crossover, hybrid}
##       x dropout_pattern {none, MCAR_25, biased_25, biased_40}
##       = 16 cells.
## Target: 500 reps per cell (8,000 fits) within a 540 s wall-clock
## budget, parallel::mclapply with mc.cores = 8 (HARD REQUIREMENT,
## matched to sister papers 04, 05, 08). Architecture A (mean
## moderation), N = 35, c.bm = 0.45, carryover_t1half = 0.5.
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
  library(parallel)
})

t_start <- Sys.time()
# Graceful abort target. Default 540 s (the 10-minute budget this
# driver was written to); override with PMSIM_BUDGET_SECS for a full
# re-run on a slower host. An aborted run writes to *-partial.rds and
# leaves the canonical outputs untouched.
budget_secs <- as.integer(Sys.getenv('PMSIM_BUDGET_SECS', '540'))

set.seed(20260508)

n_cores <- 8L
cat('Cores requested:', n_cores,
    '  detected:', parallel::detectCores(), '\n')

## -----------------------------------------------------------------
## Trial designs
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

designs <- list(
  OL                    = td_OL,
  OLBDC                 = td_OLBDC,
  traditional_crossover = td_CO,
  hybrid                = td_Hybrid
)

## -----------------------------------------------------------------
## Sample size. N is the TOTAL across randomization paths
## (NOTATION.md rule 4) and is allocated across them. This driver
## holds the per-path count at 35, so the total varies with the
## design's path count: 35 (OL), 70 (OLBDC, crossover), 140
## (hybrid). Earlier versions passed 35 to every path and recorded
## that as N. The draws are unchanged; only the label is corrected.
##
## The designs are consequently NOT matched on N, so the cross-design
## power ordering this driver produces is confounded with sample
## size. The manuscript states this in the discussion.
## -----------------------------------------------------------------

n_per_path_target <- 35L
n_total_for_design <- vapply(
  designs, function(td) n_per_path_target * length(td$trialpaths),
  integer(1))

allocate_across_paths <- function(n_total, n_paths) {
  base <- n_total %/% n_paths
  rep(base, n_paths) +
    c(rep(1L, n_total %% n_paths),
      rep(0L, n_paths - n_total %% n_paths))
}

## -----------------------------------------------------------------
## DGP parameters (matched to prior 5-min and 30-min drivers)
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
##   MCAR      -- per-visit hazard independent of outcome trajectory
##   biased    -- per-visit hazard h(t) = beta0 + beta1*delta_sx(t)^2
##                with delta_sx(t) standardised across participants
##
## In all dropout cases an alpha multiplier is solved so that the
## expected ever-dropout fraction by the last visit matches
## target_rate. Once a participant drops, all subsequent measurements
## are NA.
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
## One replicate (seeded for reproducibility under mclapply)
## -----------------------------------------------------------------

run_one_rep <- function(design_name, dropout_pattern, seed) {
  set.seed(seed)

  td <- designs[[design_name]]
  paths <- td$trialpaths
  c_bm <- 0.45
  t1half <- 0.5

  n_per_path <- allocate_across_paths(
    n_total_for_design[[design_name]], length(paths))

  mp <- data.table(
    N = NA_integer_, c.bm = c_bm,
    carryover_t1half = t1half,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.1, c.cfct = 0.05
  )

  spec <- switch(
    dropout_pattern,
    none      = list(mode = 'none',   rate = 0),
    MCAR_25   = list(mode = 'MCAR',   rate = 0.25),
    biased_25 = list(mode = 'biased', rate = 0.25),
    biased_40 = list(mode = 'biased', rate = 0.40),
    stop('Unknown dropout pattern: ', dropout_pattern)
  )

  out_na <- function(reason) {
    data.table(
      beta = NA_real_, betaSE = NA_real_, p = NA_real_,
      issingular = NA, warning = reason,
      dropout_fraction = NA_real_
    )
  }

  tryCatch(
    {
      dat_list <- vector('list', length(paths))
      drop_fracs <- numeric(length(paths))
      for (g in seq_along(paths)) {
        mp$N <- n_per_path[[g]]
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

      res <- lme_analysis(td$trialpaths, dat, op)
      res[, dropout_fraction := mean(drop_fracs)]
      res
    },
    error = function(e) out_na(conditionMessage(e))
  )
}

## -----------------------------------------------------------------
## Cell grid
## -----------------------------------------------------------------

cells <- expand.grid(
  design = c('OL', 'OLBDC',
             'traditional_crossover', 'hybrid'),
  dropout_pattern = c('none', 'MCAR_25',
                      'biased_25', 'biased_40'),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
n_cells <- nrow(cells)
n_reps_target <- 500L
chunk_size <- n_cores * 10L     # 80 fits per mclapply call

cat('Cells:', n_cells, '  target reps/cell:', n_reps_target, '\n')
cat('Total target fits:', n_cells * n_reps_target, '\n')
cat('Chunk size (per mclapply call):', chunk_size, '\n\n')

## -----------------------------------------------------------------
## Main loop: per-cell, mclapply in chunks, abort on time budget
## -----------------------------------------------------------------

cell_rows <- vector('list', n_cells)
cell_counts <- integer(n_cells)
total_fits <- 0L
last_progress <- 0L
aborted <- FALSE
abort_cell <- NA_integer_

for (i in seq_len(n_cells)) {
  ci <- cells[i, ]
  remaining <- budget_secs -
    as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
  if (remaining <= 0) {
    aborted <- TRUE
    abort_cell <- i
    cat(sprintf('\nBudget exhausted before cell %d.\n', i))
    break
  }

  cat(sprintf('\n[%4.0fs left] Cell %d/%d: design=%s dropout=%s\n',
              remaining, i, n_cells,
              ci$design, ci$dropout_pattern))

  cell_seeds <- 100000L * i + seq_len(n_reps_target)
  cell_results <- vector('list', n_reps_target)
  k_done <- 0L

  while (k_done < n_reps_target) {
    elapsed_now <- as.numeric(difftime(Sys.time(), t_start,
                                       units = 'secs'))
    if (elapsed_now > budget_secs) {
      cat(sprintf('  Budget exceeded mid-cell at %d/%d reps\n',
                  k_done, n_reps_target))
      aborted <- TRUE
      abort_cell <- i
      break
    }

    chunk_end <- min(k_done + chunk_size, n_reps_target)
    idx <- (k_done + 1L):chunk_end
    seeds_chunk <- cell_seeds[idx]

    chunk_res <- mclapply(
      seeds_chunk,
      function(s) run_one_rep(ci$design, ci$dropout_pattern, s),
      mc.cores = n_cores,
      mc.preschedule = TRUE
    )

    cell_results[idx] <- chunk_res
    k_done <- chunk_end
    total_fits <- total_fits + length(idx)

    if (total_fits - last_progress >= 200L) {
      el <- as.numeric(difftime(Sys.time(), t_start,
                                units = 'secs'))
      cat(sprintf('  ...%d fits done (cell %d, %d/%d), %.1fs\n',
                  total_fits, i, k_done, n_reps_target, el))
      last_progress <- total_fits
    }
  }

  cell_results <- cell_results[
    !vapply(cell_results, is.null, logical(1))]
  cell_counts[i] <- length(cell_results)

  if (length(cell_results) > 0L) {
    cell_dt <- rbindlist(cell_results, fill = TRUE)
    cell_dt[, `:=`(design = ci$design,
                   dropout_pattern = ci$dropout_pattern,
                   rep_idx = seq_len(.N))]
    cell_rows[[i]] <- cell_dt
  }

  if (aborted) break
}

t_total <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
cat(sprintf('\nTotal fits: %d   Wall clock: %.1f s (budget %d)\n',
            total_fits, t_total, budget_secs))
cat('Per-cell rep counts:\n')
print(setNames(cell_counts,
               paste(cells$design, cells$dropout_pattern,
                     sep = ':')))

## -----------------------------------------------------------------
## Assemble replicates table
## -----------------------------------------------------------------

reps_dt <- rbindlist(cell_rows, fill = TRUE)
setnames(reps_dt,
         old = c('beta', 'betaSE', 'p'),
         new = c('beta_bmDbc', 'betaSE_bmDbc', 'p_bmDbc'))
reps_dt[, converged := !is.na(p_bmDbc)]
reps_dt[, N := n_total_for_design[as.character(design)]]
setcolorder(reps_dt,
            c('design', 'N', 'dropout_pattern', 'rep_idx',
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

reps_out <- reps_dt[, .(design, N, dropout_pattern, rep_idx,
                        beta_bmDbc, betaSE_bmDbc, p_bmDbc,
                        converged)]

# An aborted run covers only a prefix of the cell grid; writing it to
# the canonical path would replace a complete run with a truncated
# one.
reps_path <- if (aborted) {
  'analysis/data/quick-sim/09-dropout-replicates-partial.rds'
} else {
  'analysis/data/quick-sim/09-dropout-replicates.rds'
}
saveRDS(reps_out, reps_path)
if (aborted) {
  cat(sprintf(paste0(
    '\nRUN INCOMPLETE: wrote %s and left the canonical replicates\n',
    'untouched. Re-run with a larger PMSIM_BUDGET_SECS.\n'), reps_path))
}

## -----------------------------------------------------------------
## Morris-style summary table with MCSE columns
## -----------------------------------------------------------------

alpha <- 0.05
summary_dt <- reps_dt[, {
  conv <- converged & !is.na(p_bmDbc)
  n_used <- sum(conv)
  rej <- if (n_used > 0) sum(p_bmDbc[conv] < alpha) else 0L
  rate <- if (n_used > 0) rej / n_used else NA_real_
  mcse_rate <- if (n_used > 0)
    sqrt(rate * (1 - rate) / n_used) else NA_real_
  mb <- if (n_used > 0) mean(beta_bmDbc[conv]) else NA_real_
  sb <- if (n_used > 1) sd(beta_bmDbc[conv]) else NA_real_
  mcse_mb <- if (n_used > 1) sb / sqrt(n_used) else NA_real_
  mdrop <- mean(dropout_fraction, na.rm = TRUE)
  list(
    n_reps         = .N,
    n_converged    = n_used,
    conv_rate      = n_used / .N,
    power          = rate,
    mcse_power     = mcse_rate,
    mean_dropout_fraction = mdrop,
    mean_beta      = mb,
    sd_beta        = sb,
    mcse_mean_beta = mcse_mb
  )
}, by = .(design, N, dropout_pattern)]

setorder(summary_dt, design, dropout_pattern)

write.table(summary_dt,
            file = if (aborted)
              'analysis/data/quick-sim/09-dropout-summary-partial.txt'
            else 'analysis/data/quick-sim/09-dropout-summary.txt',
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
  levels = c('OL', 'OLBDC',
             'traditional_crossover', 'hybrid'),
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
                'target %d reps per cell, 8-core mclapply'),
         n_reps_target),
       caption = paste0('Error bars: 1.96 x binomial MCSE. ',
                        'Dashed line: nominal alpha = 0.05.')) +
  theme_minimal(base_size = 11) +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 30, hjust = 1))

if (!aborted) {
  ggsave('analysis/figures/quick-sim/09-dropout-power.pdf',
         p_power, width = 11, height = 4.5)
}

cat('\nWrote:\n')
cat(' ', reps_path, '\n')
if (aborted) {
  cat('  analysis/data/quick-sim/09-dropout-summary-partial.txt\n')
  cat('  (figure not written: run incomplete)\n')
} else {
  cat('  analysis/data/quick-sim/09-dropout-summary.txt\n')
  cat('  analysis/figures/quick-sim/09-dropout-power.pdf\n')
}

if (aborted) {
  cat(sprintf('\nNOTE: budget aborted at %.1fs (cell %d); ',
              t_total, abort_cell))
  cat('some cells have fewer than target reps.\n')
} else {
  cat('\nCompleted full target rep count within budget.\n')
}
