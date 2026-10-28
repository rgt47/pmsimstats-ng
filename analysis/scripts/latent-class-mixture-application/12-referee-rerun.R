## analysis/scripts/latent-class-mixture-application/12-referee-rerun.R
##
## Rerun of the paper 03 pilot requested by the referee report of
## 2026-10-28 (analysis/report/03-latent-class-mixture-application/
## referee-report-2026-10-28.md; Major comments M1, M2, M4, M5, M6, M9
## and Suggestions 1 to 5). Modeled on
## analysis/scripts/carryover-sensitivity/36-run-referee-variants.R.
## Replaces the open-label pilot (quick-sim-driver-parallel.R,
## 03-study5-pilot.R, 10-study5-production.R), which fitted
## `Sx ~ bm + t + bm:t` in the linear-mixed arm because the all-on-drug
## design carries no within-person variation in drug state (M1).
##
## Design (M1, M11). The Hybrid design of papers 01 and 02
## (design_preset('Hybrid') in carryover-sensitivity/simulation-core.R:
## open-label titration at weeks 4 and 8, blinded discontinuation at
## weeks 9 to 12, a brief crossover at weeks 16 and 20; four paths), built
## with the package's buildtrialdesign() and N allocated across the four
## paths in the package's fixed proportions (floor(N / 4) each, the
## remainder to the first paths, as generateSimulatedResults() does; no
## randomization of path allocation). With drug state and expectancy
## both varying within person, lme_analysis() fits
##   Sx ~ bm + t + Dbc + bm:Dbc + De
## (verified by 13-check-referee-rerun.R), and Dbc is 1 on drug,
## (1/2)^(tsd / 0.2) off drug after first exposure, and 0 before exposure.
##
## Data-generating process. The package generateData() with the extended
## settings of the compendium: separable AR(1) at rho = 0.7
## (c.tv = c.pb = c.br = 0.7), cross-factor correlations
## c.cf1t = c.cfct = 0.2, the default decaying biomarker-response
## correlation (lambda = ln 2 / t_half), anchored carryover on the mean
## (package default), stopOnNonPD = TRUE (package default: a non-PD
## matrix is refused, never repaired), DGP carryover half-life 0.5 week.
## Response and baseline parameters are the compendium defaults
## (default_resp_param(), default_baseline_param() in simulation-core.R:
## BR sd 8, TV sd 6, PB sd 2; biomarker N(0, 1); baseline N(70, 100)),
## not the pilot's BR sd 5. Every analysis builds Dbc at the project's
## policy analysis half-life of 0.2 week.
##
## DGP cells (M5), each at N in {35, 70, 140}:
##   null   the single-Gaussian null, c.bm = 0;
##   cs     the class-scaling DGP of 02-heterogeneous-slopes-dgp.R (a
##          deterministic class-specific scaling of the drug-response
##          component: responders keep BR, non-responders have it
##          deleted at on-drug visits; no random slope), gating slopes
##          {1.0, 1.5};
##   rs     the true random-slopes DGP of
##          07-true-heterogeneous-slopes-dgp.R at its defaults
##          (mu_R = 1, sigma_R = 0.3, mu_NR = 0, sigma_NR = 0.3), gating
##          slopes {1.0, 1.5};
##   mvn    the correlation-encoding MVN DGP (generateData() with
##          c.bm > 0) at the c.bm matched to the marginal on-drug
##          Cor(B, BR) of the class-scaling cell at gating slope 1.5,
##          computed inside this driver from a 20,000-participant draw.
## For every DGP configuration the same large draw also records the
## marginal on-drug Cor(B, BR) at the week-8 visit (OL2, on drug on
## every path), Cohen's f^2, and the class-conditional on-drug
## variances and means of BR (calibration.rds), so the paper can put
## the cells on one effect-size scale.
##
## Analyses, all with the same fixed-effects formula (M1), on the same
## simulated trial each replicate:
##   lme_model    lme_analysis() as the package fits it: random
##                intercept, corCAR1(~t | ptID); bm:Dbc Wald with the
##                model-based SE and nlme's t reference;
##   lme_cr2      the same fit with CR2 cluster-robust SE and
##                Satterthwaite df (cr2_extract(); the compendium's
##                recommended inference);
##   hlme_gating  lcmm::hlme, ng = 2, fixed effects identical to the lme
##                fit, mixture = ~Dbc, classmb = ~bm, random intercept,
##                gridsearch(rep = 5, maxiter = 30, minit = ng-1 fit):
##                the gating coefficient 'bm class1' and its Wald z
##                (the unconditional gating test);
##   hlme_within  the within-class bm:Dbc coefficient of the same ng = 2
##                fit and its Wald z (the fourth test form, M9);
##   lrt          2 (loglik_2 - loglik_1) for ng = 2 against ng = 1.
## Per replicate the ng = 2 fit also records BIC for ng = 1 and 2 with
## the BIC-selection indicator, posterior entropy, the convergence code,
## whether the final fit hit the iteration cap (HLME_MAXITER), the
## iteration count and the elapsed time.
##
## Calibration of the LRT (M4). The naive chi-square reference for the
## one-versus-two-class LRT is invalid: under one class the gating slope
## and the class-2 parameters exist only under the alternative, the
## information matrix is singular, and the statistic has no chi-square
## limit at any N (Davies 1977, 1987; Hansen 1996; McLachlan 1987). The
## column lrt_p_naive is kept only to document the size of that error.
## The calibration is a single-resample parametric bootstrap, the
## pattern of boot_stat_f1a0() / sim_null_f() in simulation-core.R: for
## each replicate ONE dataset is simulated from the fitted ng = 1 model
## (same N, same design, same Dbc, same biomarker values; random
## intercept plus independent residuals, as hlme specifies), ng = 1 and
## ng = 2 are fitted to it by the same procedure, and its LRT statistic
## is stored (boot_lrt_stat). Pooling the resamples within a cell gives
## the bootstrap LRT (BLRT) reference: per replicate
##   p* = (1 + #{boot >= observed}) / (n_boot + 1),
## reject when p* < 0.05 (ws_size() in simulation-core.R). The
## BLRT-gated gating test rejects only when the BLRT rejects at 0.05 and
## the gating Wald then rejects at 0.05. The pooled p* is appended to
## all.rds as blrt_p_pooled once every replicate of a cell is in.
##
## Seeds (M6). Every (cell, replicate) job seeds from
## label_seed(sprintf('%s-rep-%d', cell_id, r)) inside the worker with
## mc.set.seed = FALSE; no stream is shared between cells, and the
## bootstrap resample draws from the replicate's own stream after the
## observed fits. The calibration draws seed from
## label_seed(sprintf('calib-%s', dgp_id)).
##
## Replicates: REPS_NULL (default 2,000) per null cell, REPS_ALT (default
## 1,000) per alternative cell; 6,000 null + 15,000 alternative
## replicates in all. See --time for the measured cost.
##
## Outputs (OUT_DIR, default
## analysis/data/latent-class-mixture-application/referee-rerun/):
##   trials/           one .rds per (cell, rep); resumable checkpoint,
##                     gitignored like the carryover-sensitivity caches
##   all.rds           one row per replicate and analysis
##   summary.rds       one row per (cell, test)
##   calibration.rds   large-draw correlations and class-conditional
##                     variances per DGP configuration
##   preflight.rds     smallest eigenvalue of the DGP covariance per
##                     cell (minimum over the four paths)
##   meta.rds          replicates, seed scheme, package versions, elapsed
##
## Usage (from the repo root):
##   Rscript analysis/scripts/latent-class-mixture-application/12-referee-rerun.R \
##     [--block null|alt|all] [--reps-null R] [--reps-alt R] [--cores K] \
##     [--out-dir PATH] [--smoke] [--preflight] [--time]
##   --cores defaults to the SIM_CORES environment variable (set by
##   simrun) and then to 4; --smoke runs two replicates of every cell;
##   --preflight checks every cell's covariance matrix and exits; --time
##   runs TIME_REPS replicates of one alternative cell at each N on one
##   worker and prints the per-replicate cost with a host projection.

Sys.setenv(VECLIB_MAXIMUM_THREADS = '1', OPENBLAS_NUM_THREADS = '1',
           OMP_NUM_THREADS = '1')

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(parallel)
  library(lcmm)
  library(clubSandwich)
})

## simulation-core.R supplies design_preset(), cr2_extract(),
## label_seed(), ws_size(), default_resp_param() and
## default_baseline_param(); its data-generation helpers (which need
## implementations/tidyverse) are not used here.
source('analysis/scripts/carryover-sensitivity/simulation-core.R')
source('analysis/scripts/latent-class-mixture-application/02-heterogeneous-slopes-dgp.R')
source('analysis/scripts/latent-class-mixture-application/07-true-heterogeneous-slopes-dgp.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
BLOCK     <- argv('--block', 'all')
REPS_NULL <- as.integer(argv('--reps-null', '2000'))
REPS_ALT  <- as.integer(argv('--reps-alt', '1000'))
CORES     <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '4')))
PREFLIGHT <- '--preflight' %in% args
SMOKE     <- '--smoke' %in% args
TIME      <- '--time' %in% args
OUT_DIR <- argv('--out-dir',
  if (SMOKE) 'analysis/data/latent-class-mixture-application/referee-rerun-smoke'
  else 'analysis/data/latent-class-mixture-application/referee-rerun')
ck_dir <- file.path(OUT_DIR, 'trials')
dir.create(ck_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(BLOCK %in% c('null', 'alt', 'all'))

## ---------------- fixed settings ----------------
N_VALUES        <- c(35, 70, 140)
GATING_SLOPES   <- c(1.0, 1.5)
RHO             <- 0.7
C_CF            <- 0.2
T_HALF_DGP      <- 0.5
T_HALF_ANALYSIS <- 0.2
N_CALIB         <- 20000L
REF_VISIT       <- 'OL2'
HLME_MAXITER    <- 100L
GRID_REP        <- 5L
GRID_MAXITER    <- 30L
TIME_REPS       <- 3L
HOST_SPEEDUP    <- 2
HOST_CORES      <- 8
ALPHA           <- 0.05

RP <- as.data.table(default_resp_param())
BP <- as.data.table(default_baseline_param())
OP <- list(useDE = TRUE, t_random_slope = FALSE, full_model_out = TRUE,
           carryover_t1half = T_HALF_ANALYSIS, simplecarryover = FALSE,
           carryover_scalefactor = 1)

spec <- design_preset('Hybrid')
TD <- buildtrialdesign(
  name_longform = 'Hybrid', name_shortform = 'Hybrid',
  timepoints = spec$timepoints, timeptnames = spec$timeptnames,
  expectancies = spec$expectancies, ondrug = spec$ondrug
)
DESIGN_SET <- TD$trialpaths
N_PATHS <- length(DESIGN_SET)

base_modelparam <- function(N, c_bm) {
  data.table(N = N, c.bm = c_bm, carryover_t1half = T_HALF_DGP,
             c.tv = RHO, c.pb = RHO, c.br = RHO,
             c.cf1t = C_CF, c.cfct = C_CF)
}

## The package's fixed allocation of N across paths
## (generateSimulatedResults(): floor, remainder to the first paths).
path_allocation <- function(N) {
  n_each <- N %/% N_PATHS
  n_each + c(rep(1, N %% N_PATHS), rep(0, N_PATHS - N %% N_PATHS))
}

## ---------------- DGP configurations and cells ----------------
## A DGP configuration is N-free; a cell is a configuration at one N.
dgp_config <- function(dgp, gating_slope = NA_real_, c_bm = 0) {
  id <- switch(dgp,
    null = 'null',
    cs   = sprintf('cs-gs%.1f', gating_slope),
    rs   = sprintf('rs-gs%.1f', gating_slope),
    mvn  = sprintf('mvn-cbm%.3f', c_bm))
  list(dgp_id = id, dgp = dgp, gating_slope = gating_slope, c_bm = c_bm)
}

DGP_LABELS <- c(null = 'single-gaussian-null', cs = 'class-scaling',
                rs = 'true-random-slopes', mvn = 'correlation-encoding-mvn')

## One simulated trial: N participants across the four Hybrid paths.
generate_trial <- function(cfg, N) {
  n_path <- path_allocation(N)
  rbindlist(lapply(seq_len(N_PATHS), function(g) {
    if (n_path[g] == 0) return(NULL)
    mp <- base_modelparam(n_path[g], cfg$c_bm)
    d <- switch(cfg$dgp,
      null = ,
      mvn = generateData(mp, RP, BP, DESIGN_SET[[g]], empirical = FALSE,
                         makePositiveDefinite = TRUE,
                         dgp_architecture = 'mvn'),
      cs = {
        mp[, `:=`(beta_R_factor = 1.0, beta_NR_factor = 0.0,
                  gating_intercept = 0.0, gating_slope = cfg$gating_slope)]
        generate_random_slopes_data(mp, RP, BP, DESIGN_SET[[g]])
      },
      rs = {
        mp[, `:=`(gating_intercept = 0.0, gating_slope = cfg$gating_slope)]
        generate_true_random_slopes_data(mp, RP, BP, DESIGN_SET[[g]])
      })
    d[, path := g]
    d
  }), fill = TRUE)
}

## ---------------- calibration: large draw per DGP configuration ----------------
## Marginal on-drug Cor(B, BR) at the week-8 visit (OL2.br, after any
## class scaling), Cohen's f^2 = r^2 / (1 - r^2), and the
## class-conditional on-drug variance and mean of BR. Under the MVN DGP
## Cor(bm, OL2.br) equals c.bm by construction, so the matched c.bm for
## the mvn cell is the class-scaling value at gating slope 1.5.
calibrate_config <- function(cfg) {
  set.seed(label_seed(sprintf('calib-%s', cfg$dgp_id)))
  d <- generate_trial(cfg, N_CALIB)
  br <- d[[paste0(REF_VISIT, '.br')]]
  r <- cor(d$bm, br)
  has_class <- 'latent_class' %in% names(d)
  z <- if (has_class) d$latent_class else rep(NA_integer_, nrow(d))
  stat_class <- function(k, f) if (has_class) f(br[z == k]) else NA_real_
  ## Under the class-scaling DGP the non-responder class has BR deleted
  ## at on-drug visits (zero variance), so its within-class correlation
  ## is undefined.
  cor_class <- function(k) {
    if (!has_class || var(br[z == k]) == 0) return(NA_real_)
    cor(d$bm[z == k], br[z == k])
  }
  data.table(
    dgp_id = cfg$dgp_id, dgp = DGP_LABELS[[cfg$dgp]],
    gating_slope = cfg$gating_slope, c_bm = cfg$c_bm,
    n_draw = nrow(d), ref_visit = REF_VISIT,
    cor_bm_br = r, f2 = r^2 / (1 - r^2),
    var_br_all = var(br), mean_br_all = mean(br),
    prop_class1 = if (has_class) mean(z) else NA_real_,
    var_br_class0 = stat_class(0, var), var_br_class1 = stat_class(1, var),
    mean_br_class0 = stat_class(0, mean), mean_br_class1 = stat_class(1, mean),
    cor_bm_br_class0 = cor_class(0), cor_bm_br_class1 = cor_class(1)
  )
}

configs <- list(dgp_config('null'))
for (gs in GATING_SLOPES) configs <- c(configs, list(dgp_config('cs', gs)))
for (gs in GATING_SLOPES) configs <- c(configs, list(dgp_config('rs', gs)))
calib <- rbindlist(lapply(configs, calibrate_config))
c_bm_matched <- calib[dgp_id == sprintf('cs-gs%.1f', max(GATING_SLOPES)), cor_bm_br]
cfg_mvn <- dgp_config('mvn', c_bm = c_bm_matched)
configs <- c(configs, list(cfg_mvn))
calib <- rbind(calib, calibrate_config(cfg_mvn))
attr(calib, 'c_bm_matched') <- c_bm_matched
attr(calib, 'matched_to') <- sprintf('cs-gs%.1f', max(GATING_SLOPES))
saveRDS(calib, file.path(OUT_DIR, 'calibration.rds'))
cat(sprintf('Calibration (N = %d per configuration), week-8 Cor(B, BR); matched c.bm = %.4f\n',
            N_CALIB, c_bm_matched))
print(as.data.frame(calib[, .(dgp_id, gating_slope, c_bm = round(c_bm, 4),
                              cor_bm_br = round(cor_bm_br, 4), f2 = round(f2, 4),
                              prop_class1 = round(prop_class1, 3),
                              var_br_class0 = round(var_br_class0, 2),
                              var_br_class1 = round(var_br_class1, 2))]),
      row.names = FALSE)

grid <- rbindlist(lapply(configs, function(cfg) {
  data.table(block = if (cfg$dgp == 'null') 'null' else 'alt',
             cell_id = sprintf('%s-N%d', cfg$dgp_id, N_VALUES),
             dgp_id = cfg$dgp_id, dgp = DGP_LABELS[[cfg$dgp]],
             gating_slope = cfg$gating_slope, c_bm = cfg$c_bm,
             N = N_VALUES,
             reps = if (cfg$dgp == 'null') REPS_NULL else REPS_ALT,
             cfg = rep(list(cfg), length(N_VALUES)))
}))
if (BLOCK != 'all') grid <- grid[block == BLOCK]
stopifnot(nrow(grid) > 0, !anyDuplicated(grid$cell_id))
if (SMOKE) grid[, reps := 2L]

## ---------------- pre-flight: the DGP covariance must be PD on every path ----------------
## The class scaling and the random slopes act on the MVN draw after it
## is taken and do not enter the covariance, so the check is of the
## backbone matrix (c.bm = 0 for those cells, the matched c.bm for mvn).
preflight_cell <- function(cfg, N) {
  mp <- base_modelparam(N, cfg$c_bm)
  rbindlist(lapply(seq_len(N_PATHS), function(g) {
    b <- buildSigma(mp, RP, BP, DESIGN_SET[[g]], makePositiveDefinite = FALSE,
                    stopOnNonPD = FALSE)
    sds <- sqrt(diag(b$sigma))
    pos <- sds > 0
    ev <- eigen(b$sigma[pos, pos], symmetric = TRUE, only.values = TRUE)$values
    data.table(path = g, min_eigenvalue = min(ev), feasible = min(ev) > 0)
  }))
}
pre <- grid[, {
  pf <- preflight_cell(cfg[[1]], N)
  .(min_eigenvalue = min(pf$min_eigenvalue), feasible = all(pf$feasible))
}, by = .(block, cell_id, c_bm)]
saveRDS(pre, file.path(OUT_DIR, 'preflight.rds'))
cat(sprintf('\nPre-flight: %d cells, %d feasible, min eigenvalue %.4f\n',
            nrow(pre), sum(pre$feasible), min(pre$min_eigenvalue)))
if (!all(pre$feasible)) {
  print(as.data.frame(pre[feasible == FALSE]), row.names = FALSE)
  stop('Infeasible cells present; adjust the grid before running.')
}
if (PREFLIGHT) {
  print(as.data.frame(pre[, .(block, cell_id, c_bm = round(c_bm, 4),
                              min_eigenvalue = round(min_eigenvalue, 4),
                              feasible)]), row.names = FALSE)
  quit(save = 'no')
}

## ---------------- per-replicate analyses ----------------
na_row <- function(analysis) {
  data.table(analysis = analysis, estimate = NA_real_, std_error = NA_real_,
             statistic = NA_real_, p_value = NA_real_, converged = FALSE)
}

## hlme's variance-covariance of the estimates: VarCov() when it works,
## otherwise the vech unpacking used by 01-lcmm-wrapper.R.
hlme_vcov <- function(fit) {
  est <- fit$best
  V <- tryCatch(VarCov(fit), error = function(e) NULL)
  if (is.matrix(V) && nrow(V) == length(est)) {
    dimnames(V) <- list(names(est), names(est))
    return(V)
  }
  n_par <- length(est)
  V <- matrix(0, n_par, n_par)
  V[upper.tri(V, diag = TRUE)] <- as.numeric(fit$V)
  V <- V + t(V) - diag(diag(V))
  dimnames(V) <- list(names(est), names(est))
  V
}

fit_hlme_pair <- function(form, df) {
  m1 <- tryCatch(
    hlme(fixed = form, random = ~1, subject = 'ptID', ng = 1, data = df,
         verbose = FALSE, maxiter = HLME_MAXITER),
    error = function(e) NULL)
  if (is.null(m1)) return(list(m1 = NULL, m2 = NULL))
  m2 <- tryCatch(
    gridsearch(
      m = hlme(fixed = form, mixture = ~Dbc, random = ~1, subject = 'ptID',
               ng = 2, classmb = ~bm, data = df, verbose = FALSE,
               maxiter = HLME_MAXITER),
      rep = GRID_REP, maxiter = GRID_MAXITER, minit = m1),
    error = function(e) NULL)
  list(m1 = m1, m2 = m2)
}

hlme_ok <- function(m) !is.null(m) && isTRUE(m$conv == 1)
hlme_cap <- function(m) {
  if (is.null(m)) return(NA)
  isTRUE(m$conv == 2) || isTRUE(m$niter >= HLME_MAXITER)
}
lrt_of <- function(m1, m2) {
  if (is.null(m1) || is.null(m2)) return(NA_real_)
  2 * (as.numeric(m2$loglik) - as.numeric(m1$loglik))
}

## One parametric resample from the fitted ng = 1 hlme model: the same
## design matrix (participants, visits, Dbc, De, biomarker values), a
## random intercept N(0, varcov 1) per participant and independent
## N(0, stderr^2) residuals, which is exactly the ng = 1 model hlme fits.
sim_null_hlme <- function(m1, form, df) {
  est <- estimates(m1, cholesky = FALSE)
  X <- model.matrix(form[-2], df)
  colnames(X)[colnames(X) == '(Intercept)'] <- 'intercept'
  stopifnot(all(colnames(X) %in% names(est)))
  mu <- drop(X %*% est[colnames(X)])
  tau <- sqrt(max(as.numeric(est['varcov 1']), 0))
  sig <- as.numeric(est['stderr'])
  ids <- unique(df$ptID)
  b <- rnorm(length(ids), 0, tau)[match(df$ptID, ids)]
  ds <- df
  ds$Sx <- mu + b + rnorm(nrow(df), 0, sig)
  ds
}

run_one_rep <- function(cfg, N, cell_id, r) {
  seed <- label_seed(sprintf('%s-rep-%d', cell_id, r))
  set.seed(seed)
  dat <- generate_trial(cfg, N)

  t0 <- Sys.time()
  prep <- tryCatch(lme_analysis(DESIGN_SET, dat, OP), error = function(e) NULL)
  t_lme <- as.numeric(Sys.time() - t0, units = 'secs')
  fit <- if (is.null(prep)) NULL else prep$fit
  form <- if (is.null(prep)) NULL else prep$form
  dl <- if (is.null(prep)) NULL else prep$datamerged
  lme_formula <- if (is.null(form)) NA_character_ else paste(deparse(form), collapse = '')
  lme_corcar1 <- if (is.null(fit)) NA else
    inherits(fit$modelStruct$corStruct, 'corCAR1')

  rows <- list()
  if (is.null(fit)) {
    rows$lme_model <- na_row('lme_model')
    rows$lme_cr2 <- na_row('lme_cr2')
  } else {
    tt <- summary(fit)$tTable
    tgt <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(tt))
    rows$lme_model <- if (length(tgt) == 0) na_row('lme_model') else
      data.table(analysis = 'lme_model', estimate = tt[tgt[1], 'Value'],
                 std_error = tt[tgt[1], 'Std.Error'],
                 statistic = tt[tgt[1], 't-value'],
                 p_value = tt[tgt[1], 'p-value'], converged = TRUE)
    cr <- cr2_extract(fit, dl, c('bm:Dbc', 'Dbc:bm'))
    rows$lme_cr2 <- if (is.null(cr) || !is.finite(cr$p_value)) na_row('lme_cr2') else
      data.table(analysis = 'lme_cr2', estimate = cr$estimate,
                 std_error = cr$std_error, statistic = cr$estimate / cr$std_error,
                 p_value = cr$p_value, converged = TRUE)
  }

  ## hlme on the identical fixed-effects formula and the identical rows.
  t0 <- Sys.time()
  pair <- if (is.null(form)) list(m1 = NULL, m2 = NULL) else {
    df <- as.data.frame(dl)
    df <- df[complete.cases(df[, all.vars(form)]), ]
    fit_hlme_pair(form, df)
  }
  t_hlme <- as.numeric(Sys.time() - t0, units = 'secs')
  m1 <- pair$m1; m2 <- pair$m2
  m2_ok <- hlme_ok(m2)
  wald_row <- function(name, target) {
    if (!m2_ok) return(na_row(name))
    est <- m2$best
    tg <- intersect(target, names(est))
    if (length(tg) == 0) return(na_row(name))
    V <- hlme_vcov(m2)
    se <- suppressWarnings(sqrt(V[tg[1], tg[1]]))
    if (!is.finite(se) || se <= 0) return(na_row(name))
    z <- as.numeric(est[tg[1]]) / se
    data.table(analysis = name, estimate = as.numeric(est[tg[1]]), std_error = se,
               statistic = z, p_value = 2 * pnorm(-abs(z)), converged = TRUE)
  }
  rows$hlme_gating <- wald_row('hlme_gating', c('bm class1', 'bm class2'))
  rows$hlme_within <- wald_row('hlme_within', c('bm:Dbc', 'Dbc:bm'))

  lrt_stat <- if (hlme_ok(m1) && m2_ok) lrt_of(m1, m2) else NA_real_
  lrt_df <- if (hlme_ok(m1) && m2_ok) length(m2$best) - length(m1$best) else NA_integer_
  rows$lrt <- data.table(analysis = 'lrt', estimate = lrt_stat, std_error = NA_real_,
                         statistic = lrt_stat, p_value = NA_real_,
                         converged = hlme_ok(m1) && m2_ok)

  entropy <- NA_real_
  if (m2_ok) {
    P <- as.matrix(m2$pprob[, grep('^prob', names(m2$pprob)), drop = FALSE])
    entropy <- 1 - mean(-rowSums(ifelse(P > 0, P * log(P), 0))) / log(2)
  }

  ## One parametric-bootstrap null resample from the ng = 1 fit.
  t0 <- Sys.time()
  boot <- list(m1 = NULL, m2 = NULL)
  if (hlme_ok(m1)) {
    ds <- tryCatch(sim_null_hlme(m1, form, df), error = function(e) NULL)
    if (!is.null(ds)) boot <- fit_hlme_pair(form, ds)
  }
  t_boot <- as.numeric(Sys.time() - t0, units = 'secs')
  boot_ok <- hlme_ok(boot$m1) && hlme_ok(boot$m2)

  out <- rbindlist(rows)
  out[, `:=`(
    cell_id = cell_id, rep = r, seed = seed, N = N,
    bic_ng1 = if (is.null(m1)) NA_real_ else as.numeric(m1$BIC),
    bic_ng2 = if (is.null(m2)) NA_real_ else as.numeric(m2$BIC),
    conv_ng1 = if (is.null(m1)) 0L else as.integer(m1$conv),
    conv_ng2 = if (is.null(m2)) 0L else as.integer(m2$conv),
    niter_ng1 = if (is.null(m1)) NA_integer_ else as.integer(m1$niter),
    niter_ng2 = if (is.null(m2)) NA_integer_ else as.integer(m2$niter),
    hit_cap_ng1 = hlme_cap(m1), hit_cap_ng2 = hlme_cap(m2),
    entropy = entropy, lrt_stat = lrt_stat, lrt_df = lrt_df,
    lrt_p_naive = if (is.na(lrt_stat)) NA_real_ else
      pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE),
    boot_lrt_stat = if (boot_ok) lrt_of(boot$m1, boot$m2) else NA_real_,
    boot_conv_ng1 = if (is.null(boot$m1)) 0L else as.integer(boot$m1$conv),
    boot_conv_ng2 = if (is.null(boot$m2)) 0L else as.integer(boot$m2$conv),
    boot_hit_cap_ng2 = hlme_cap(boot$m2),
    lme_formula = lme_formula, lme_corcar1 = lme_corcar1,
    t_lme = t_lme, t_hlme = t_hlme, t_boot = t_boot
  )]
  out[, `:=`(bic_select2 = bic_ng2 < bic_ng1,
             delta_bic = bic_ng1 - bic_ng2)]
  setcolorder(out, c('cell_id', 'rep', 'seed', 'N', 'analysis'))
  out
}

## ---------------- --time: per-replicate cost at each N ----------------
if (TIME) {
  cfg_t <- grid[dgp_id == sprintf('cs-gs%.1f', max(GATING_SLOPES))]
  if (nrow(cfg_t) == 0) cfg_t <- grid[block == 'alt']
  if (nrow(cfg_t) == 0) cfg_t <- grid
  timing <- rbindlist(lapply(seq_len(nrow(cfg_t)), function(i) {
    t0 <- Sys.time()
    res <- rbindlist(lapply(seq_len(TIME_REPS), function(r)
      run_one_rep(cfg_t$cfg[[i]], cfg_t$N[i], cfg_t$cell_id[i], r)))
    secs <- as.numeric(Sys.time() - t0, units = 'secs')
    data.table(cell_id = cfg_t$cell_id[i], N = cfg_t$N[i],
               secs_per_rep = secs / TIME_REPS,
               t_lme = mean(res$t_lme), t_hlme = mean(res$t_hlme),
               t_boot = mean(res$t_boot),
               n_rows = nrow(res) / TIME_REPS)
  }))
  cat(sprintf('\nTiming on one local worker (%d replicates per N):\n', TIME_REPS))
  print(as.data.frame(timing[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)]),
        row.names = FALSE)
  full <- grid[, .(reps = sum(reps)), by = N][timing, on = 'N', nomatch = 0]
  full[, worker_secs := reps * secs_per_rep]
  local_h <- sum(full$worker_secs) / 3600
  cat(sprintf('Grid as configured (block %s): %d replicates, %.1f local worker-hours\n',
              BLOCK, sum(grid$reps), local_h))
  cat(sprintf('Host projection (%d cores, %.1fx faster per core): %.1f h wall\n',
              HOST_CORES, HOST_SPEEDUP, local_h / HOST_SPEEDUP / HOST_CORES))
  quit(save = 'no')
}

## ---------------- run ----------------
cat(sprintf('\nReferee rerun, block %s: %d cells, %d (cell, rep) jobs, %d workers\n',
            BLOCK, nrow(grid), sum(grid$reps), CORES))
jobs <- grid[, .(rep = seq_len(reps)), by = cell_id]

one_job <- function(i) {
  cell_id <- jobs$cell_id[i]; r <- jobs$rep[i]
  f <- file.path(ck_dir, sprintf('%s-rep-%04d.rds', cell_id, r))
  if (file.exists(f)) return(invisible(NULL))
  k <- match(cell_id, grid$cell_id)
  res <- tryCatch(run_one_rep(grid$cfg[[k]], grid$N[k], cell_id, r),
                  error = function(e) NULL)
  if (!is.null(res)) saveRDS(res, f)
  invisible(NULL)
}

t0 <- Sys.time()
invisible(mclapply(seq_len(nrow(jobs)), one_job, mc.cores = CORES,
                   mc.preschedule = FALSE, mc.set.seed = FALSE))
elapsed_min <- as.numeric(Sys.time() - t0, units = 'mins')

files <- list.files(ck_dir, full.names = TRUE)
files <- files[sub('-rep-\\d+\\.rds$', '', basename(files)) %in% grid$cell_id]
res <- rbindlist(lapply(files, readRDS), fill = TRUE)
cell_meta <- grid[, .(block, cell_id, dgp_id, dgp, gating_slope, c_bm)]
res <- cell_meta[res, on = 'cell_id']
n_done <- uniqueN(res[, .(cell_id, rep)])
cat(sprintf('%d of %d (cell, rep) jobs present in %.1f min (this run)\n',
            n_done, nrow(jobs), elapsed_min))

## Pooled single-resample BLRT within each cell: p* for each observed
## LRT against every finite bootstrap LRT of the same cell.
pooled_p <- function(obs, boot) {
  s <- boot[is.finite(boot)]
  if (length(s) == 0) return(rep(NA_real_, length(obs)))
  vapply(obs, function(a) {
    if (!is.finite(a)) return(NA_real_)
    (1 + sum(s >= a)) / (length(s) + 1)
  }, numeric(1))
}
res[, blrt_p_pooled := pooled_p(lrt_stat, boot_lrt_stat[analysis == 'lrt']),
    by = cell_id]
res[, n_boot_pooled := sum(is.finite(boot_lrt_stat[analysis == 'lrt'])), by = cell_id]
setcolorder(res, c('block', 'cell_id', 'dgp_id', 'dgp', 'gating_slope', 'c_bm',
                   'N', 'rep', 'seed', 'analysis'))
saveRDS(res, file.path(OUT_DIR, 'all.rds'))

## ---------------- summary ----------------
## Per (cell, test): rejection rate over converged replicates with its
## binomial MCSE, non-convergence over all replicates, mean estimate,
## mean SE, empirical SD, kappa = mean SE / empirical SD, the fraction
## of ng = 2 fits at the iteration cap, the BIC two-class selection
## rate, the pooled-BLRT rejection rate and the mean entropy (the last
## four over the replicates whose ng = 2 fit converged).
wald_tests <- c('lme_model', 'lme_cr2', 'hlme_gating', 'hlme_within')
rep_level <- res[analysis == 'lrt']
rep_level[, blrt_reject := is.finite(blrt_p_pooled) & blrt_p_pooled < ALPHA]
gating <- res[analysis == 'hlme_gating', .(cell_id, rep, gating_p = p_value,
                                            gating_conv = converged)]
rep_level <- gating[rep_level, on = c('cell_id', 'rep')]
rep_level[, gated_reject := blrt_reject & gating_conv & gating_p < ALPHA]

mcse <- function(p, n) sqrt(p * (1 - p) / n)
summ_wald <- res[analysis %in% wald_tests, {
  ok <- converged & is.finite(p_value)
  n <- sum(ok)
  rej <- mean(p_value[ok] < ALPHA)
  .(n = n, n_total = .N, reject = rej, mcse = mcse(rej, n),
    nonconv = mean(!ok),
    mean_est = mean(estimate[ok]), mean_se = mean(std_error[ok]),
    emp_sd = sd(estimate[ok]), kappa = mean(std_error[ok]) / sd(estimate[ok]),
    hit_cap_frac = if (analysis[1] %in% c('hlme_gating', 'hlme_within'))
      mean(hit_cap_ng2[conv_ng2 != 0], na.rm = TRUE) else NA_real_)
}, by = .(cell_id, test = analysis)]

summ_lrt <- rep_level[, {
  ok <- converged & is.finite(lrt_stat)
  rbind(
    data.table(test = 'blrt', n = sum(ok & is.finite(blrt_p_pooled)), n_total = .N,
               reject = mean(blrt_reject[ok]),
               mcse = mcse(mean(blrt_reject[ok]), sum(ok)),
               nonconv = mean(!ok), mean_est = mean(lrt_stat[ok]),
               mean_se = NA_real_, emp_sd = sd(lrt_stat[ok]), kappa = NA_real_,
               hit_cap_frac = mean(hit_cap_ng2[ok])),
    data.table(test = 'blrt_gated_gating', n = sum(ok), n_total = .N,
               reject = mean(gated_reject[ok]),
               mcse = mcse(mean(gated_reject[ok]), sum(ok)),
               nonconv = mean(!ok), mean_est = NA_real_, mean_se = NA_real_,
               emp_sd = NA_real_, kappa = NA_real_,
               hit_cap_frac = mean(hit_cap_ng2[ok])),
    data.table(test = 'lrt_naive_chisq', n = sum(ok), n_total = .N,
               reject = mean(lrt_p_naive[ok] < ALPHA),
               mcse = mcse(mean(lrt_p_naive[ok] < ALPHA), sum(ok)),
               nonconv = mean(!ok), mean_est = mean(lrt_stat[ok]),
               mean_se = NA_real_, emp_sd = sd(lrt_stat[ok]), kappa = NA_real_,
               hit_cap_frac = mean(hit_cap_ng2[ok]))
  )
}, by = cell_id]

cell_diag <- rep_level[, {
  ok <- conv_ng2 == 1
  .(bic_select2_rate = mean(bic_select2[ok], na.rm = TRUE),
    blrt_reject_rate = mean(blrt_reject[ok]),
    entropy_mean = mean(entropy[ok], na.rm = TRUE),
    n_boot_pooled = n_boot_pooled[1],
    boot_nonconv = mean(!is.finite(boot_lrt_stat)),
    mean_t_lme = mean(t_lme), mean_t_hlme = mean(t_hlme), mean_t_boot = mean(t_boot))
}, by = cell_id]

summary_tab <- rbind(summ_wald, summ_lrt)
summary_tab <- cell_diag[summary_tab, on = 'cell_id']
summary_tab <- cell_meta[summary_tab, on = 'cell_id']
summary_tab <- grid[, .(cell_id, N)][summary_tab, on = 'cell_id']
test_order <- c(wald_tests, 'blrt', 'blrt_gated_gating', 'lrt_naive_chisq')
summary_tab[, test := factor(test, levels = test_order)]
setorder(summary_tab, block, dgp_id, N, test)
summary_tab[, test := as.character(test)]
setcolorder(summary_tab, c('block', 'cell_id', 'dgp_id', 'dgp', 'gating_slope',
                           'c_bm', 'N', 'test', 'n', 'n_total', 'reject', 'mcse',
                           'nonconv', 'mean_est', 'mean_se', 'emp_sd', 'kappa',
                           'hit_cap_frac', 'bic_select2_rate', 'blrt_reject_rate',
                           'entropy_mean'))
saveRDS(summary_tab, file.path(OUT_DIR, 'summary.rds'))

cat('\nRejection rate per (cell, test); kappa = mean SE / empirical SD:\n')
print(as.data.frame(summary_tab[, .(cell_id, test, n, reject = round(reject, 3),
                                    mcse = round(mcse, 3), nonconv = round(nonconv, 3),
                                    kappa = round(kappa, 2),
                                    bic2 = round(bic_select2_rate, 3),
                                    blrt = round(blrt_reject_rate, 3),
                                    entropy = round(entropy_mean, 3))]),
      row.names = FALSE, max = 5000)

meta <- list(
  script = '12-referee-rerun.R', block = BLOCK, date = Sys.time(),
  reps_null = REPS_NULL, reps_alt = REPS_ALT, smoke = SMOKE, cores = CORES,
  seed_scheme = 'label_seed(sprintf("%s-rep-%d", cell_id, r)) per (cell, rep); mc.set.seed = FALSE',
  design = 'Hybrid (design_preset) via buildtrialdesign; fixed path allocation',
  t_half_dgp = T_HALF_DGP, t_half_analysis = T_HALF_ANALYSIS, rho = RHO, c_cf = C_CF,
  c_bm_matched = c_bm_matched, n_calib = N_CALIB,
  hlme = list(maxiter = HLME_MAXITER, grid_rep = GRID_REP, grid_maxiter = GRID_MAXITER),
  elapsed_min_this_run = elapsed_min, n_jobs_present = n_done, n_jobs = nrow(jobs),
  r_version = R.version.string,
  packages = sapply(c('lcmm', 'nlme', 'clubSandwich', 'data.table', 'pmsimstats'),
                    function(p) tryCatch(as.character(packageVersion(p)),
                                         error = function(e) NA_character_)),
  git_sha = tryCatch(system('git rev-parse --short HEAD', intern = TRUE),
                     error = function(e) NA_character_)
)
saveRDS(meta, file.path(OUT_DIR, 'meta.rds'))
cat(sprintf('\nWritten: %s\n', OUT_DIR))
