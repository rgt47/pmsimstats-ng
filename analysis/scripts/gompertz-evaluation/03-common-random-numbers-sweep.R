# 03-common-random-numbers-sweep.R
#
# Paper 07 production driver: the trajectory-family comparison under
# common random numbers.
#
# The earlier driver (02-faithful-trajectory-sweep.R) gave every one of
# its 16 cells its own seed stream, so a between-family contrast was a
# difference of two independent binomial proportions (SE about 0.019
# at 1,000 replicates). The 2026-10-28 referee report showed that the
# 0.039 spread it reported was Monte Carlo noise. This driver runs the
# same design, calibration and parameters but uses the same seed for
# the same replicate index in every family and both architectures, so
# the families are compared on identical random inputs and the
# contrasts are paired.
#
# Cells.
#   Alternative (c.bm = 0.45): 4 families x 2 architectures, 1,000
#     matched seeds each (8,000 fits).
#   Null (c.bm = 0): 4 families, 2,000 matched seeds each (8,000 fits).
#     At c.bm = 0 the mean-moderation and covariance-moderation DGPs
#     coincide (no mean shift is applied and no biomarker-response
#     correlation is built), so the null is run once per family and
#     labeled architecture 'both'. A 20-seed check that the two
#     architectures return identical results at the null is saved with
#     the summary.
#
# Design: two-path OL+BDC, visits at weeks 4, 8, 12, 16, 17, 18, 19,
# 20; discontinuation at week 18 (path A) or 19 (path B); N = 70 total
# (35 per path); t_{1/2} = 0.
#
# Outputs (analysis/data/quick-sim/):
#   07-gompertz-crn-replicates.rds   per-replicate decisions, p-values
#                                    and estimates
#   07-gompertz-crn-paired.rds       paired contrasts against Gompertz
#   07-gompertz-crn-summary.txt      cell summary and run metadata
#   07-gompertz-crn-meta.rds         run metadata (elapsed, workers,
#                                    calibration, reference value)
#
# Usage (from the repository root, about 8 minutes on 8 cores):
#   Rscript analysis/scripts/gompertz-evaluation/03-common-random-numbers-sweep.R
#   Rscript ... --smoke        (20 seeds per cell; no files written)
#
# Author: pmsimstats team
# Last updated: 2026-10-28

suppressMessages({
  devtools::load_all('.', quiet = TRUE)
  library(data.table)
  library(parallel)
})

args  <- commandArgs(trailingOnly = TRUE)
smoke <- '--smoke' %in% args

n_workers <- 8L
n_alt  <- if (smoke) 20L else 1000L
n_null <- if (smoke) 20L else 2000L
n_chk  <- 20L

t_master_start <- Sys.time()
out_data_dir <- 'analysis/data/quick-sim'
dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)

# -- Trial design (two-path OL+BDC) ---------------------------------

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

data(extracted_bp); data(extracted_rp)
bp <- extracted_bp
rp <- extracted_rp

op <- list(useDE = FALSE, t_random_slope = FALSE,
           full_model_out = FALSE, carryover_t1half = 0,
           simplecarryover = FALSE, carryover_scalefactor = 1)

# -- Trajectory-family calibration (anchor at week 5, matched m) ----

br_pars <- as.list(rp[cat == 'br'])
maxr_v  <- br_pars$max
g_disp  <- br_pars$disp
g_rate  <- br_pars$rate
sd_br   <- br_pars$sd
sd_bm   <- bp[cat == 'bm']$sd

t_anchor <- 5
g_anchor <- modgompertz(t_anchor, maxr_v, g_disp, g_rate)

solve_t0 <- function(shape_fun, target, anchor, rate, maxr) {
  uniroot(function(t0) shape_fun(anchor, maxr, rate, t0) - target,
          interval = c(-5, 20), tol = 1e-10)$root
}
log_rate  <- 0.5
log_t0    <- solve_t0(shape_logistic, g_anchor, t_anchor, log_rate, maxr_v)
htan_rate <- 0.5
htan_t0   <- solve_t0(shape_hyperbolic_tangent, g_anchor, t_anchor,
                      htan_rate, maxr_v)
pl_post_slope <- 0
pl_tbp        <- t_anchor * maxr_v / g_anchor

family_param <- list(
  gompertz = list(p2 = g_disp,    p3 = g_rate),
  logistic = list(p2 = log_rate,  p3 = log_t0),
  hyperbolic_tangent = list(p2 = htan_rate, p3 = htan_t0),
  piecewise_linear_breakpoint = list(p2 = pl_tbp, p3 = pl_post_slope)
)
families <- names(family_param)

for (fam in families) {
  pp <- family_param[[fam]]
  v <- trajectoryShape(fam, t_anchor, maxr_v, pp$p2, pp$p3)
  stopifnot(abs(v - g_anchor) < 1e-4)
}

# Closed-form reference for the interaction coefficient on the raw
# biomarker scale with D_bc binary at t_{1/2} = 0 (referee M5, 2026-10-28):
# the conditional mean of BR on drug given B_i has slope
# c_bm * sigma_BR / sigma_bm, and the sign convention (NOTATION rule 5)
# makes the fitted coefficient negative.
c_bm_alt  <- 0.45
beta_ref  <- -c_bm_alt * sd_br / sd_bm

cat(sprintf('-- Calibration: anchor week %.1f, BR(anchor)=%.4f, m=%.4f\n',
            t_anchor, g_anchor, maxr_v))
cat(sprintf('-- Reference interaction: -c_bm sigma_BR/sigma_bm = %.4f\n',
            beta_ref))

architectures <- c('mvn', 'mean_moderation')

N_total <- 70L
allocate_across_paths <- function(n_total, n_paths) {
  base <- n_total %/% n_paths
  rep(base, n_paths) +
    c(rep(1L, n_total %% n_paths),
      rep(0L, n_paths - n_total %% n_paths))
}
n_per_path <- allocate_across_paths(N_total, n_paths)

mp_for <- function(cv, n) data.table(
  N = n, c.bm = cv, carryover_t1half = 0,
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.1)

# -- Sigma cache keyed by (family, architecture, c.bm, path) --------

build_cache <- function(fam) {
  pp <- family_param[[fam]]
  cache <- list()
  for (arch in architectures) {
    for (cv in c(0, c_bm_alt)) {
      key <- paste(arch, cv, sep = '|')
      cache[[key]] <- lapply(seq_len(n_paths), function(g)
        buildSigma(mp_for(cv, n_per_path[[g]]), rp, bp,
                   td$trialpaths[[g]],
                   makePositiveDefinite = TRUE, dgp_architecture = arch,
                   br_family = fam, br_p2 = pp$p2, br_p3 = pp$p3,
                   stopOnNonPD = TRUE))
    }
  }
  cache
}
fam_caches <- setNames(lapply(families, build_cache), families)

# -- Single replicate ------------------------------------------------
# moderation_scaling = 'trajectory' under mean moderation (the
# multiplicative form of paper 01); the argument is inert under 'mvn'.

run_one <- function(family, architecture, cbm, seed) {
  set.seed(seed)
  ms <- if (architecture == 'mean_moderation') 'trajectory' else 'constant'
  cs_set <- fam_caches[[family]][[paste(architecture, cbm, sep = '|')]]
  dats <- vector('list', n_paths)
  for (g in seq_len(n_paths)) {
    dat <- generateData(mp_for(cbm, n_per_path[[g]]), rp, bp,
                        td$trialpaths[[g]],
                        empirical = FALSE, makePositiveDefinite = TRUE,
                        cached_sigma = cs_set[[g]],
                        dgp_architecture = architecture,
                        moderation_scaling = ms)
    dat[, path := g]
    dats[[g]] <- dat
  }
  dat <- rbindlist(dats, fill = TRUE)
  res <- lme_analysis(td$trialpaths, dat, op)
  data.table(family = family, architecture = architecture, N = N_total,
             c.bm = cbm, seed = seed,
             beta_bmDbc = res$beta, se_bmDbc = res$betaSE,
             p_bmDbc = res$p, reject = isTRUE(res$p < 0.05),
             converged = !is.na(res$beta))
}

# One seed, every cell that shares it. The seed is reset inside
# run_one, so every family and architecture draws the same Z matrix.
run_seed_alt <- function(seed) rbindlist(lapply(architectures, function(a)
  rbindlist(lapply(families, function(f) run_one(f, a, c_bm_alt, seed)))))

run_seed_null <- function(seed) rbindlist(lapply(families, function(f) {
  r <- run_one(f, 'mvn', 0, seed)
  r[, architecture := 'both']
  r
}))

# -- Alternative cells: 1,000 common seeds -------------------------

seeds_alt <- 810001L + seq_len(n_alt) - 1L
cat(sprintf('\n-- Alternative: %d seeds x %d families x %d architectures\n',
            n_alt, length(families), length(architectures)))
alt <- rbindlist(mclapply(seeds_alt, run_seed_alt, mc.cores = n_workers))
cat(sprintf('   done (elapsed %.0f s)\n',
            as.numeric(Sys.time() - t_master_start, units = 'secs')))

# -- Null cells: 2,000 common seeds, run once per family -------------

seeds_null <- 820001L + seq_len(n_null) - 1L
cat(sprintf('\n-- Null: %d seeds x %d families\n', n_null, length(families)))
nul <- rbindlist(mclapply(seeds_null, run_seed_null, mc.cores = n_workers))
cat(sprintf('   done (elapsed %.0f s)\n',
            as.numeric(Sys.time() - t_master_start, units = 'secs')))

# -- Null architecture-equivalence check ------------------------------
# At c.bm = 0 the two architectures must return identical fits from
# the same seed.

seeds_chk <- 830001L + seq_len(n_chk) - 1L
chk <- rbindlist(mclapply(seeds_chk, function(s) rbindlist(lapply(families,
  function(f) rbindlist(lapply(architectures,
    function(a) run_one(f, a, 0, s))))), mc.cores = n_workers))
chk_w <- dcast(chk, family + seed ~ architecture, value.var = 'p_bmDbc')
chk_wb <- dcast(chk, family + seed ~ architecture, value.var = 'beta_bmDbc')
null_identical <- isTRUE(all.equal(chk_w$mvn, chk_w$mean_moderation,
                                   tolerance = 0)) &&
  isTRUE(all.equal(chk_wb$mvn, chk_wb$mean_moderation, tolerance = 0))
cat(sprintf('\n-- Null equivalence check (%d seeds x 4 families): %s\n',
            n_chk, if (null_identical) 'identical' else 'DIFFER'))

replicates <- rbind(alt, nul)

# -- Cell summary ------------------------------------------------------

summarise_cell <- function(df, cbm) {
  conv <- df[converged == TRUE]
  n <- nrow(conv)
  power <- mean(conv$reject)
  bias <- if (cbm > 0) mean(conv$beta_bmDbc) - beta_ref else
    mean(conv$beta_bmDbc)
  data.table(n_reps = nrow(df), power = power,
             mcse_power = sqrt(power * (1 - power) / n),
             mean_beta = mean(conv$beta_bmDbc),
             bias = bias,
             mcse_bias = sd(conv$beta_bmDbc) / sqrt(n),
             sd_beta = sd(conv$beta_bmDbc),
             mean_se = mean(conv$se_bmDbc),
             conv_rate = n / nrow(df))
}
summary_dt <- replicates[, summarise_cell(.SD, .BY$c.bm),
                         by = .(family, architecture, c.bm)]
setorder(summary_dt, c.bm, architecture, family)

# -- Paired contrasts against Gompertz -------------------------------
# Decisions: McNemar-type SE, sqrt((b + c))/n, where b and c are the
# discordant counts. Estimates: paired-difference SE.

paired_vs_gompertz <- function(df) {
  w  <- dcast(df, seed ~ family, value.var = 'reject')
  wp <- dcast(df, seed ~ family, value.var = 'p_bmDbc')
  wb <- dcast(df, seed ~ family, value.var = 'beta_bmDbc')
  n  <- nrow(w)
  rbindlist(lapply(setdiff(families, 'gompertz'), function(f) {
    b <- sum(w[[f]] & !w$gompertz)
    cc <- sum(!w[[f]] & w$gompertz)
    d <- mean(w[[f]]) - mean(w$gompertz)
    se <- sqrt((b + cc)) / n
    db <- wb[[f]] - wb$gompertz
    data.table(family = f, n = n,
               power_diff = d, discordant_up = b, discordant_down = cc,
               paired_se = se, z = if (se > 0) d / se else NA_real_,
               max_abs_p_diff = max(abs(wp[[f]] - wp$gompertz)),
               beta_diff = mean(db), beta_diff_se = sd(db) / sqrt(n))
  }))
}
paired_dt <- replicates[, paired_vs_gompertz(.SD),
                        by = .(architecture, c.bm)]

# Seeds on which any family disagrees with any other, per architecture.
discord <- replicates[, {
  w <- dcast(.SD, seed ~ family, value.var = 'reject')
  R <- as.matrix(w[, -1])
  .(n_seeds = nrow(R),
    n_discordant = sum(apply(R, 1, function(x) length(unique(x)) > 1)))
}, by = .(architecture, c.bm)]

# Paired difference between architectures within family at the
# alternative (same seed, same family).
arch_pair <- alt[, {
  w <- dcast(.SD, seed ~ architecture, value.var = 'reject')
  b <- sum(w$mvn & !w$mean_moderation)
  cc <- sum(!w$mvn & w$mean_moderation)
  d <- mean(w$mvn) - mean(w$mean_moderation)
  .(power_diff = d, discordant = b + cc, paired_se = sqrt(b + cc) / nrow(w))
}, by = family]

elapsed <- as.numeric(Sys.time() - t_master_start, units = 'secs')

meta <- list(
  generated = format(Sys.time(), '%Y-%m-%d %H:%M %Z'),
  elapsed_s = elapsed, workers = n_workers,
  total_fits = nrow(replicates) + nrow(chk),
  n_alt = n_alt, n_null = n_null, n_chk = n_chk,
  N_total = N_total, n_per_path = n_per_path,
  t_anchor = t_anchor, g_anchor = g_anchor, maxr = maxr_v,
  family_param = family_param, beta_ref = beta_ref,
  sd_br = sd_br, sd_bm = sd_bm, c_bm_alt = c_bm_alt,
  null_identical = null_identical,
  seeds = list(alt = range(seeds_alt), null = range(seeds_null),
               chk = range(seeds_chk)),
  r_version = R.version.string,
  platform = R.version$platform,
  blas = extSoftVersion()[['BLAS']],
  smoke = smoke)

if (!smoke) {
  saveRDS(replicates, file.path(out_data_dir, '07-gompertz-crn-replicates.rds'))
  saveRDS(list(paired = paired_dt, discordant = discord,
               arch_pair = arch_pair, summary = summary_dt),
          file.path(out_data_dir, '07-gompertz-crn-paired.rds'))
  saveRDS(meta, file.path(out_data_dir, '07-gompertz-crn-meta.rds'))
  sink(file.path(out_data_dir, '07-gompertz-crn-summary.txt'))
  cat('# 07-gompertz-evaluation common-random-numbers summary\n')
  cat(sprintf('# Generated %s\n', meta$generated))
  cat(sprintf('# Wall-clock elapsed: %.1f s | total fits: %d | workers: %d\n',
              elapsed, meta$total_fits, n_workers))
  cat(sprintf('# %s | %s | BLAS %s\n', meta$r_version, meta$platform,
              meta$blas))
  cat(sprintf(paste0('# Trial design: two-path OL+BDC, N=%d total across %d ',
                     'paths (%s per path), t1/2=0; anchor week %.1f\n'),
              N_total, n_paths, paste(n_per_path, collapse = '/'), t_anchor))
  cat(sprintf('# Alternative: c.bm=%.2f, %d common seeds per family and architecture\n',
              c_bm_alt, n_alt))
  cat(sprintf('# Null: c.bm=0, %d common seeds per family; architectures coincide\n',
              n_null))
  cat(sprintf('# Null equivalence check (%d seeds): %s\n', n_chk,
              if (null_identical) 'identical across architectures' else 'DIFFER'))
  cat(sprintf('# Reference interaction at c.bm=%.2f: %.4f points per biomarker unit\n',
              c_bm_alt, beta_ref))
  cat('#\n# Cell summary (alpha = 0.05):\n')
  print(summary_dt, digits = 4)
  cat('\n# Paired contrasts against Gompertz (same seeds):\n')
  print(paired_dt, digits = 4)
  cat('\n# Seeds with any between-family disagreement:\n')
  print(discord)
  cat('\n# Paired architecture contrast within family (mvn minus mean):\n')
  print(arch_pair, digits = 4)
  sink()
}

cat('\n--- Cell summary ---\n')
print(summary_dt, digits = 4)
cat('\n--- Paired contrasts against Gompertz ---\n')
print(paired_dt, digits = 4)
cat('\n--- Discordant seeds ---\n')
print(discord)
cat('\n--- Architecture contrast within family ---\n')
print(arch_pair, digits = 4)
cat(sprintf('\nElapsed %.1f s; smoke = %s\n', elapsed, smoke))
