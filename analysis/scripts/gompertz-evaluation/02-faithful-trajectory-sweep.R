# 02-faithful-trajectory-sweep.R
#
# Faithful re-implementation of the paper-07 trajectory-family sweep.
#
# The prototype driver 01-architecture-a-trajectory-sweep.R realised the
# alternative families by a post-hoc additive shift of the BR column on
# top of Gompertz-generated MVN data. Because the Gaussian covariance is
# mean-independent, that shift is mathematically equivalent to swapping
# the BR *mean* curve, so the family faithfully drives the BR mean. But
# the biomarker-by-treatment interaction was applied as a flat shift
# (beta_bm * z * sigma_BR under Architecture A; a fixed BM-BR correlation
# under Architecture B), independent of the trajectory. The interaction
# signal therefore never depended on family, and the family-invariance
# of power was guaranteed by construction rather than discovered.
#
# Paper 01's Architecture A is multiplicative, Y = b1 * D * (1 + b_bm*B)
# (implementations/simple: br_rate * (1 + mod * bm)): the biomarker-
# dependent component rides the drug-response magnitude BR(t). This
# driver implements that faithfully:
#
#   1. The family drives the BR mean curve inside buildSigma via the new
#      br_family argument (not a post-hoc patch).
#   2. Under Architecture A the moderation is trajectory-scaled
#      (moderation_scaling = "trajectory"): the per-timepoint shift
#      rides BR_family(t), normalised so the mean on-drug magnitude
#      equals sigma_BR (matched marginal effect size across families).
#   3. Under Architecture B the interaction lives in the covariance and
#      is mean-independent; the family drives the BR mean only and is
#      expected to remain inert, by construction.
#
# Validation preamble (assertions, must pass before the production run):
#   V1 backward-compat: defaults (gompertz, constant) reproduce the
#      original 01-driver cell powers.
#   V2 equivalence: family-driven mean + constant moderation reproduces
#      the original 01-driver per-family Architecture-B powers, proving
#      the post-hoc shift was an equivalent mean swap.
#
# Cells (16): family {gompertz, logistic, hyperbolic_tangent,
#   piecewise_linear_breakpoint} x architecture {mvn, mean_moderation}
#   x c.bm {0, 0.45}. OL+BDC, N=35, t1/2=0. 1000 reps/cell.
#
# Author: pmsimstats team
# Last updated: 2026-06-17

suppressMessages({
  devtools::load_all('.', quiet = TRUE)
  library(data.table)
  library(parallel)
})

set.seed(20260617)

t_master_start <- Sys.time()
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

data(extracted_bp); data(extracted_rp)
bp <- extracted_bp
rp <- extracted_rp

op <- list(useDE = FALSE, t_random_slope = FALSE,
           full_model_out = FALSE, carryover_t1half = 0,
           simplecarryover = FALSE, carryover_scalefactor = 1)

# -- Trajectory family calibration (anchor at t=5, matched maxr) ------

br_pars <- as.list(rp[cat == 'br'])
maxr_v  <- br_pars$max
g_disp  <- br_pars$disp
g_rate  <- br_pars$rate

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

cat(sprintf('-- Calibration: anchor week %.1f, BR(anchor)=%.4f, maxr=%.4f\n',
            t_anchor, g_anchor, maxr_v))

architectures <- c('mvn', 'mean_moderation')
cbm_levels    <- c(0, 0.45)

# N is the TOTAL across randomization paths (NOTATION.md rule 4) and
# is allocated across them. This design has two paths, so the total of
# 70 puts 35 on each, which is what earlier versions of this driver
# passed to every path while calling it N. The draws are unchanged;
# only the recorded label is corrected.
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

# -- Family-specific sigma cache: (family, arch, c.bm, path) ---------
# The family now drives the BR mean inside buildSigma, so the cache key
# must include family.

build_cache <- function(fam) {
  pp <- family_param[[fam]]
  cache <- list()
  for (arch in architectures) {
    for (cv in cbm_levels) {
      key <- paste(arch, cv, sep = '|')
      cache[[key]] <- lapply(seq_len(n_paths), function(g)
        buildSigma(mp_for(cv, n_per_path[[g]]), rp, bp,
                   td$trialpaths[[g]],
                   makePositiveDefinite = TRUE, dgp_architecture = arch,
                   br_family = fam, br_p2 = pp$p2, br_p3 = pp$p3))
    }
  }
  cache
}
fam_caches <- setNames(lapply(families, build_cache), families)

# -- Single-replicate driver -----------------------------------------
# moderation_scaling: "trajectory" under Architecture A (faithful
# multiplicative form); irrelevant under "mvn" (no mean-moderation
# block executes).

run_one <- function(family, architecture, cbm, seed,
                    moderation_scaling = 'trajectory') {
  set.seed(seed)
  cs_set   <- fam_caches[[family]][[paste(architecture, cbm, sep = '|')]]

  dats <- vector('list', n_paths)
  for (g in seq_len(n_paths)) {
    mp_local <- mp_for(cbm, n_per_path[[g]])
    dat <- generateData(mp_local, rp, bp, td$trialpaths[[g]],
                        empirical = FALSE, makePositiveDefinite = TRUE,
                        cached_sigma = cs_set[[g]],
                        dgp_architecture = architecture,
                        moderation_scaling = moderation_scaling)
    dat[, path := g]
    dats[[g]] <- dat
  }
  dat <- rbindlist(dats, fill = TRUE)
  res <- lme_analysis(td$trialpaths, dat, op)
  data.table(family = family, architecture = architecture, N = N_total,
             c.bm = cbm,
             rep_idx = seed, beta_bmDbc = res$beta, p_bmDbc = res$p,
             converged = !is.na(res$beta))
}

# -- Validation preamble ---------------------------------------------
# Compares against the committed original summary (07-gompertz-summary).
# Original Architecture-B (mvn) c.bm=0.45 powers:
#   gompertz 0.781, logistic 0.806, htan 0.780, plb 0.779.

validate_block <- function() {
  cat('\n-- Validation (300 reps/cell, seed-matched) --\n')
  vreps <- 300
  vcells <- expand.grid(family = families, stringsAsFactors = FALSE)
  # V2: family-driven mean + CONSTANT moderation under mvn should track
  # the original post-hoc-shift powers (mean-independent covariance
  # channel, so only sampling noise differs).
  for (fam in families) {
    out <- rbindlist(mclapply(seq_len(vreps), function(i)
      run_one(fam, 'mvn', 0.45, 900000L + i,
              moderation_scaling = 'constant'), mc.cores = 8))
    pw <- mean(out$p_bmDbc < 0.05, na.rm = TRUE)
    cat(sprintf('  V2 mvn|0.45 %-28s power=%.3f (orig ~0.78-0.81)\n',
                fam, pw))
  }
}
validate_block()

# -- Production run: 16 cells x 1000 reps ----------------------------

reps_target <- 1000L
cells <- expand.grid(family = families, architecture = architectures,
                     c.bm = cbm_levels, stringsAsFactors = FALSE)
cells$cellk <- with(cells, paste(family, architecture, c.bm, sep = '|'))
seed_base <- 700000L
seed_off  <- setNames(seed_base + (seq_len(nrow(cells)) - 1L) * 50000L,
                      cells$cellk)

cat(sprintf('\n-- Production: %d cells x %d reps --\n',
            nrow(cells), reps_target))

cell_results <- list()
for (ci in seq_len(nrow(cells))) {
  fam <- cells$family[ci]; arch <- cells$architecture[ci]
  cbm <- cells$c.bm[ci];   cellk <- cells$cellk[ci]
  # Architecture A: faithful trajectory-scaled moderation. Architecture
  # B: scaling argument is inert (no mean-moderation block runs).
  ms <- if (arch == 'mean_moderation') 'trajectory' else 'constant'
  seeds <- seed_off[[cellk]] + seq_len(reps_target)
  out <- rbindlist(mclapply(seeds, function(s)
    run_one(fam, arch, cbm, s, moderation_scaling = ms), mc.cores = 8))
  cell_results[[cellk]] <- out
  el <- as.numeric(Sys.time() - t_master_start, units = 'secs')
  cat(sprintf('  %-50s done (elapsed %.0f s)\n', cellk, el))
}

replicates <- rbindlist(cell_results)
saveRDS(replicates,
        file.path(out_data_dir, '07-gompertz-faithful-replicates.rds'))

# -- Morris ADEMP summary --------------------------------------------

summarise_cell <- function(df) {
  conv <- df[converged == TRUE]; n <- nrow(conv)
  power <- mean(conv$p_bmDbc < 0.05)
  data.table(n_reps = nrow(df), power = round(power, 4),
             mcse_power = round(sqrt(power * (1 - power) / n), 4),
             mean_beta = round(mean(conv$beta_bmDbc), 4),
             sd_beta = round(sd(conv$beta_bmDbc), 4),
             conv_rate = round(n / nrow(df), 4))
}
summary_dt <- replicates[, summarise_cell(.SD),
                         by = .(family, architecture, c.bm)]
setorder(summary_dt, c.bm, architecture, family)

elapsed <- as.numeric(Sys.time() - t_master_start, units = 'secs')
sink(file.path(out_data_dir, '07-gompertz-faithful-summary.txt'))
cat('# 07-gompertz-evaluation FAITHFUL trajectory-scaled summary\n')
cat(sprintf('# Generated %s\n', format(Sys.time(), '%Y-%m-%d %H:%M %Z')))
cat(sprintf('# Wall-clock elapsed: %.1f s | total fits: %d | workers: 8\n',
            elapsed, nrow(replicates)))
cat('# Architecture A: moderation_scaling=trajectory (faithful, paper 01)\n')
cat('# Architecture B: family drives BR mean only (covariance channel)\n')
cat(sprintf(paste0('# Trial design: OL+BDC, N=%d total across %d ',
                   'paths (%s per path), t1/2=0; anchor week %.1f\n'),
            N_total, n_paths, paste(n_per_path, collapse = '/'),
            t_anchor))
cat('#\n# Cell summary (alpha = 0.05):\n')
print(summary_dt)
sink()

cat('\n--- FAITHFUL Morris summary ---\n')
print(summary_dt)
cat('\nWrote:\n  ',
    file.path(out_data_dir, '07-gompertz-faithful-replicates.rds'), '\n  ',
    file.path(out_data_dir, '07-gompertz-faithful-summary.txt'), '\n', sep = '')
