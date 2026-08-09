# 01-dgp-prototype.R
#
# Tier-F item 25: paper 01 (dgp-mean-moderation-vs-mvn) at a total
# of N = 70, to escape the saturation regime that compressed the
# covariance-versus-mean-moderation power-loss-under-carryover gap
# in the earlier prototype at twice that size. N is the TOTAL across
# randomization paths and every design runs at the same total, so
# the design comparison is matched on N.
#
# Cells (24 total): architecture {mvn, mean_moderation} x
# design {CO, Hybrid, OL+BDC} x t1half {0, 0.5, 1.0} x
# c.bm {0, 0.45}.  Reps/cell = 1000 (MCSE 0.016 at p = 0.5).
# Total fits = 24,000.
#
# Wall budget: PMSIM_BUDGET_SECS, default 5400 s. Abort 50 s before
# it. The original 900 s prototype budget truncates the full 36-cell
# grid at about cell 16 on this host, so it is no longer the default.
# An aborted run writes to *-partial.rds and leaves the canonical
# outputs untouched, so a truncated run cannot overwrite a complete
# one.
#
# Parallelism policy (HARD REQUIREMENT):
#   - Primary: mclapply with mc.cores = 8.
#   - On per-cell catastrophic failure (zero converged fits in
#     a parallel chunk), automatic step-down to mc.cores = 4.
#   - On second catastrophic failure, step-down to sequential at
#     200 reps per cell.
#
# Outputs (overwritten):
#   analysis/data/quick-sim/01-dgp-replicates.rds
#   analysis/data/quick-sim/01-dgp-summary.txt
#   analysis/figures/quick-sim/01-dgp-power.pdf

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
  library(parallel)
})

t_start <- Sys.time()
budget_secs <- as.integer(Sys.getenv('PMSIM_BUDGET_SECS', '5400'))
abort_secs <- budget_secs - 50L
target_reps <- 1000L

set.seed(20260509)

# Parallelism state. Mutable across cells: a catastrophic chunk
# failure steps the whole run down for all subsequent cells.
n_cores <- 8L
fallback_seq <- FALSE
fallback_seq_reps <- 200L

cat(sprintf(
  'Parallelism: mclapply with mc.cores = %d (host has %d cores).\n',
  n_cores, parallel::detectCores()))

# ---------------------------------------------------------------
# Trial designs: CO, Hybrid (N-of-1), OL+BDC.  Templates match
# analysis/scripts/figure4/01-generate-power-sweep.R so the
# present run is comparable with the figure-4 reproduction.
# ---------------------------------------------------------------

td_CO <- buildtrialdesign(
  name_longform = 'traditional crossover',
  name_shortform = 'CO',
  timepoints = cumulative(rep(2.5, 8)),
  timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
  expectancies = rep(.5, 8),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
)

td_Hybrid <- buildtrialdesign(
  name_longform = 'hybrid N-of-1 design',
  name_shortform = 'Hybrid',
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptnames = c('OL1', 'OL2', 'BD1', 'BD2',
                  'BD3', 'BD4', 'COd', 'COp'),
  expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
  )
)

td_OLBDC <- buildtrialdesign(
  name_longform = 'open label+blinded discontinuation',
  name_shortform = 'OL+BDC',
  timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptnames = c('OL1', 'OL2', 'OL3', 'OL4',
                  'BD1', 'BD2', 'BD3', 'BD4'),
  expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)

td_lookup <- list(CO = td_CO, Hybrid = td_Hybrid, `OL+BDC` = td_OLBDC)

# ---------------------------------------------------------------
# Sample size. N is the TOTAL across randomization paths (see
# analysis/report/NOTATION.md rule 4) and is allocated across them,
# remainder spread one per path.
#
# All three designs run at the SAME total, so the cross-design
# comparison is matched on N. Earlier versions fixed the per-path
# count at 35 instead, which gave the four-path Hybrid design twice
# the sample of the two-path designs and confounded any design
# difference with sample size. At a total of 70 the two-path designs
# still put 35 on each path, so their cells are unchanged; Hybrid
# now splits 70 as 18/18/17/17 rather than running 35 per path.
# ---------------------------------------------------------------

N_TOTAL <- 70L
n_total_for_design <- c(CO = N_TOTAL, Hybrid = N_TOTAL,
                        `OL+BDC` = N_TOTAL)

allocate_across_paths <- function(n_total, n_paths) {
  base <- n_total %/% n_paths
  rep(base, n_paths) +
    c(rep(1L, n_total %% n_paths),
      rep(0L, n_paths - n_total %% n_paths))
}

# ---------------------------------------------------------------
# Parameter sets: package's extracted reference parameters
# (prazosin/PTSD calibration).
# ---------------------------------------------------------------

data(extracted_bp)
data(extracted_rp)
bp <- extracted_bp
rp <- extracted_rp

# ---------------------------------------------------------------
# Cell grid (24 cells)
# ---------------------------------------------------------------

cells <- expand.grid(
  architecture = c('mvn', 'mean_moderation'),
  design       = c('CO', 'Hybrid', 'OL+BDC'),
  t1half       = c(0, 0.5, 1.0),
  c.bm         = c(0, 0.45),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
n_cells <- nrow(cells)
cat('Cells:', n_cells, '\n')
print(cells)

# ---------------------------------------------------------------
# Single-rep worker. Returns a one-row data.table with the
# canonical lme_analysis output (beta, betaSE, p, ...).
# ---------------------------------------------------------------

run_one_rep <- function(architecture, design, c_bm, t1half, seed) {
  set.seed(seed)

  out_na <- function(reason) {
    data.table(
      beta = NA_real_, betaSE = NA_real_, p = NA_real_,
      issingular = NA, warning = reason
    )
  }

  tryCatch({
    td <- td_lookup[[design]]
    paths <- td$trialpaths
    n_per_path <- allocate_across_paths(n_total_for_design[[design]],
                                        length(paths))
    mp <- data.table(
      N = NA_integer_, c.bm = c_bm,
      carryover_t1half = t1half,
      c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
      c.cf1t = 0.2, c.cfct = 0.1
    )
    dat_list <- vector('list', length(paths))
    for (g in seq_along(paths)) {
      mp$N <- n_per_path[[g]]
      di <- generateData(
        modelparam = mp, respparam = rp, blparam = bp,
        trialdesign = paths[[g]],
        empirical = FALSE, makePositiveDefinite = TRUE,
        dgp_architecture = architecture
      )
      di[, path := g]
      dat_list[[g]] <- di
    }
    dat <- rbindlist(dat_list, fill = TRUE)

    op <- list(
      useDE = FALSE,
      t_random_slope = FALSE,
      full_model_out = FALSE,
      carryover_t1half = t1half,
      simplecarryover = FALSE,
      carryover_scalefactor = 1
    )

    res <- lme_analysis(td$trialpaths, dat, op)
    if (!('beta' %in% names(res))) return(out_na('no beta col'))
    res
  }, error = function(e) out_na(conditionMessage(e)))
}

# ---------------------------------------------------------------
# Main loop. Outer loop over cells; inner mclapply chunks for
# mid-cell budget checks. Step down on catastrophic chunk failure.
# ---------------------------------------------------------------

all_rows <- list()
total_fits <- 0L
last_progress <- 0L
aborted <- FALSE
abort_cell <- NA_integer_

for (i in seq_len(n_cells)) {
  ci <- cells[i, ]

  elapsed_now <- as.numeric(difftime(Sys.time(), t_start,
                                     units = 'secs'))
  remaining <- abort_secs - elapsed_now
  if (remaining <= 0) {
    aborted <- TRUE
    abort_cell <- i
    cat(sprintf('\nAbort threshold reached before cell %d.\n', i))
    break
  }

  cat(sprintf(
    paste0('\n[%4.0fs left] Cell %d/%d: arch=%s design=%s ',
           't1half=%.2f c.bm=%.2f\n'),
    remaining, i, n_cells,
    ci$architecture, ci$design, ci$t1half, ci$c.bm))

  n_reps <- if (fallback_seq) fallback_seq_reps else target_reps
  # mclapply chunks of n_cores * 25 so we can probe budget mid-cell
  chunk_size <- if (fallback_seq) {
    min(n_reps, 50L)
  } else {
    min(n_reps, n_cores * 25L)
  }

  cell_rows <- vector('list', n_reps)
  k_done <- 0L
  cell_seed_base <- 100000L * i

  while (k_done < n_reps) {
    elapsed_now <- as.numeric(difftime(Sys.time(), t_start,
                                       units = 'secs'))
    if (elapsed_now > abort_secs) {
      cat(sprintf('  abort threshold mid-cell at %d/%d reps\n',
                  k_done, n_reps))
      aborted <- TRUE
      abort_cell <- i
      break
    }

    chunk_end <- min(k_done + chunk_size, n_reps)
    idx <- (k_done + 1L):chunk_end
    seeds_chunk <- cell_seed_base + idx

    if (fallback_seq) {
      chunk_res <- lapply(seeds_chunk, function(s)
        run_one_rep(ci$architecture, ci$design,
                    ci$c.bm, ci$t1half, s))
    } else {
      chunk_res <- mclapply(
        seeds_chunk,
        function(s) run_one_rep(ci$architecture, ci$design,
                                ci$c.bm, ci$t1half, s),
        mc.cores = n_cores,
        mc.preschedule = TRUE
      )

      # Catastrophic failure check: zero converged fits in the
      # entire chunk. Step down once, retry the chunk.
      n_conv <- sum(vapply(chunk_res, function(r) {
        is.data.table(r) && !is.na(r$p[[1]])
      }, logical(1)))
      if (n_conv == 0L && length(chunk_res) >= 8L) {
        if (n_cores > 4L) {
          cat('  Catastrophic chunk failure; stepping down ',
              'mc.cores from ', n_cores, ' to 4.\n', sep = '')
          n_cores <- 4L
          chunk_res <- mclapply(
            seeds_chunk,
            function(s) run_one_rep(ci$architecture, ci$design,
                                    ci$c.bm, ci$t1half, s),
            mc.cores = n_cores,
            mc.preschedule = TRUE
          )
          n_conv2 <- sum(vapply(chunk_res, function(r) {
            is.data.table(r) && !is.na(r$p[[1]])
          }, logical(1)))
          if (n_conv2 == 0L) {
            cat('  Still failing at mc.cores = 4; falling ',
                'back to sequential at ', fallback_seq_reps,
                ' reps/cell.\n', sep = '')
            fallback_seq <- TRUE
            n_reps <- fallback_seq_reps
            chunk_size <- min(n_reps, 50L)
            cell_rows <- vector('list', n_reps)
            k_done <- 0L
            next
          }
        } else {
          cat('  mc.cores = 4 also failed; falling back ',
              'sequential at ', fallback_seq_reps, ' reps/cell.\n',
              sep = '')
          fallback_seq <- TRUE
          n_reps <- fallback_seq_reps
          chunk_size <- min(n_reps, 50L)
          cell_rows <- vector('list', n_reps)
          k_done <- 0L
          next
        }
      }
    }

    # Defensive: any non-data.table (e.g., mc.preschedule error
    # objects) gets coerced to NA row.
    chunk_res <- lapply(chunk_res, function(r) {
      if (is.data.table(r)) r else
        data.table(beta = NA_real_, betaSE = NA_real_,
                   p = NA_real_,
                   issingular = NA, warning = 'mc error')
    })

    cell_rows[idx] <- chunk_res
    k_done <- chunk_end
    total_fits <- total_fits + length(idx)

    if (total_fits - last_progress >= 500L) {
      el <- as.numeric(difftime(Sys.time(), t_start,
                                units = 'secs'))
      cat(sprintf('  ...%d fits done (cell %d, %d/%d), %.1fs\n',
                  total_fits, i, k_done, n_reps, el))
      last_progress <- total_fits
    }
  }

  cell_rows <- cell_rows[seq_len(k_done)]
  cell_rows <- cell_rows[!vapply(cell_rows, is.null, logical(1))]
  if (length(cell_rows) == 0L) next

  cell_dt <- rbindlist(cell_rows, fill = TRUE)
  cell_dt[, `:=`(architecture = ci$architecture,
                 design       = ci$design,
                 N            = n_total_for_design[[ci$design]],
                 c.bm         = ci$c.bm,
                 t1half       = ci$t1half,
                 rep_idx      = seq_len(.N))]
  all_rows[[i]] <- cell_dt

  if (aborted) break
}

reps_dt <- rbindlist(all_rows, fill = TRUE)
reps_dt[, converged := !is.na(p)]
setcolorder(reps_dt, c('architecture', 'design', 'N', 't1half', 'c.bm',
                       'rep_idx', 'beta', 'betaSE', 'p',
                       'converged'))

t_total <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
cat(sprintf('\nTotal wall-clock: %.1f s (budget %d s)\n',
            t_total, budget_secs))
cat(sprintf('Total reps stored: %d (aborted=%s)\n',
            nrow(reps_dt), aborted))
cat('Reps per cell:\n')
print(reps_dt[, .N, by = .(architecture, design, N, t1half, c.bm)])

# ---------------------------------------------------------------
# Persist replicates
# ---------------------------------------------------------------

dir.create('analysis/data/quick-sim', showWarnings = FALSE,
           recursive = TRUE)
dir.create('analysis/figures/quick-sim', showWarnings = FALSE,
           recursive = TRUE)

reps_out <- reps_dt[, .(architecture, design, N, t1half, c.bm, rep_idx,
                        beta, betaSE, p, converged)]

# A run that hit the wall budget covers only a prefix of the cell
# grid. Writing it to the canonical path would silently replace a
# complete run with a truncated one, so it goes to a partial file
# instead and the canonical outputs are left alone.
reps_path <- if (aborted) {
  'analysis/data/quick-sim/01-dgp-replicates-partial.rds'
} else {
  'analysis/data/quick-sim/01-dgp-replicates.rds'
}
saveRDS(reps_out, reps_path)
if (aborted) {
  cat(sprintf(paste0(
    '\nRUN INCOMPLETE: %d of %d cells covered. Wrote %s and left the\n',
    'canonical replicates untouched. Re-run with a larger\n',
    'PMSIM_BUDGET_SECS to refresh them.\n'),
    abort_cell - 1L, n_cells, reps_path))
}

# ---------------------------------------------------------------
# Morris-style summary table
# ---------------------------------------------------------------

alpha <- 0.05
summary_dt <- reps_dt[, {
  conv <- converged & !is.na(p)
  n_used <- sum(conv)
  rej   <- if (n_used > 0L) sum(p[conv] < alpha) else 0L
  rate  <- if (n_used > 0L) rej / n_used else NA_real_
  mcse_rate <- if (n_used > 0L)
    sqrt(rate * (1 - rate) / n_used) else NA_real_
  mb <- if (n_used > 0L) mean(beta[conv]) else NA_real_
  sb <- if (n_used > 1L) sd(beta[conv]) else NA_real_
  mcse_mb <- if (n_used > 1L) sb / sqrt(n_used) else NA_real_
  list(
    n_reps         = .N,
    n_converged    = n_used,
    conv_rate      = n_used / .N,
    power          = rate,
    mcse_power     = mcse_rate,
    mean_beta      = mb,
    sd_beta        = sb,
    mcse_mean_beta = mcse_mb
  )
}, by = .(architecture, design, N, t1half, c.bm)]

setorder(summary_dt, architecture, design, c.bm, t1half)
fwrite(summary_dt,
       if (aborted) 'analysis/data/quick-sim/01-dgp-summary-partial.txt'
       else 'analysis/data/quick-sim/01-dgp-summary.txt',
       sep = '\t')

cat('\n--- Summary ---\n')
print(summary_dt)

# ---------------------------------------------------------------
# Power figure: power vs t1half, faceted by architecture x design,
# two effect sizes overlaid.
# ---------------------------------------------------------------

plot_dt <- copy(summary_dt)
plot_dt[, ci_lo := pmax(0, power - 1.96 * mcse_power)]
plot_dt[, ci_hi := pmin(1, power + 1.96 * mcse_power)]
plot_dt[, arch_lab := ifelse(architecture == 'mvn',
                              'B: MVN (correlation)',
                              'A: mean moderation')]
plot_dt[, cbm_lab := factor(paste0('c.bm = ', c.bm),
                            levels = c('c.bm = 0',
                                       'c.bm = 0.45'))]
plot_dt[, design_lab := factor(design,
                               levels = c('CO', 'Hybrid', 'OL+BDC'))]

p_power <- ggplot(plot_dt,
                  aes(x = t1half, y = power,
                      colour = cbm_lab, group = cbm_lab)) +
  geom_hline(yintercept = 0.05, linetype = 'dashed', alpha = 0.4) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.04) +
  facet_grid(arch_lab ~ design_lab) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0)) +
  labs(x = 'Carryover half-life t1/2 (weeks)',
       y = 'Power (P(p < 0.05))',
       colour = 'Effect size',
       title = 'Power for biomarker x treatment interaction',
       subtitle = sprintf(
         paste0('N=35, target %d reps/cell ',
                '(mclapply mc.cores=%d)'),
         target_reps, n_cores),
       caption = paste0('Error bars: 1.96 x MCSE (binomial). ',
                        'Dashed = nominal alpha = 0.05.')) +
  theme_bw(base_size = 11) +
  theme(legend.position = 'bottom')

if (!aborted) {
  ggsave('analysis/figures/quick-sim/01-dgp-power.pdf',
         p_power, width = 9, height = 6)
}

cat('\nWrote:\n')
cat(' ', reps_path, '\n')
if (aborted) {
  cat('  analysis/data/quick-sim/01-dgp-summary-partial.txt\n')
  cat('  (figure not written: run incomplete)\n')
} else {
  cat('  analysis/data/quick-sim/01-dgp-summary.txt\n')
  cat('  analysis/figures/quick-sim/01-dgp-power.pdf\n')
}

# ---------------------------------------------------------------
# Run-mode reporting (honest)
# ---------------------------------------------------------------

cat('\n--- Run mode ---\n')
cat(sprintf('Parallelism mode used: %s\n',
            if (fallback_seq) 'sequential (fallback)'
            else sprintf('mclapply mc.cores=%d', n_cores)))
cat(sprintf('Total fits attempted: %d\n', total_fits))
cat(sprintf('Total fits stored:    %d\n', nrow(reps_dt)))
cat(sprintf('Mean conv rate across cells: %.3f\n',
            mean(summary_dt$conv_rate)))
cat(sprintf('Wall-clock: %.1f s of %d s budget (abort %d s)\n',
            t_total, budget_secs, abort_secs))
if (aborted) {
  cat(sprintf('NOTE: aborted at cell %d on time budget.\n',
              abort_cell))
}
