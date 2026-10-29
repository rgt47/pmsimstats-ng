## 01-run-cr2-null.R
##
## Null calibration of the CR2 cluster-robust test (Satterthwaite df)
## for the biomarker-by-drug interaction, under paper 09's driver and
## paper 02's simulation core. Run through 01-run.sh from the
## repository root (it pins BLAS and OpenMP threads to 1).
##
## Data generation and model fitting reuse the drivers' own code:
##   p09: the function and object definitions of
##        analysis/scripts/informative-dropout-by-design/
##        20-paired-dropout.R (every top-level expression before the
##        `units <-` line, evaluated in their own environment), i.e.
##        generateData() with the cached Sigma, path_long_template(),
##        build_long() and fit_analysis().
##   p02: analysis/scripts/carryover-sensitivity/simulation-core.R
##        (with implementations/tidyverse/R/functions.R), i.e.
##        generate_data_multi_path() and prepare_long_data() at the
##        tier-1 null cell (mvn, c_bm = 0, t1half = 1, exponential).
## Every analysis variant is fitted with p09's fit_analysis(), which
## for (t_a = 0, use_de = FALSE) on p02's long data is paper 02's G6
## (E1cr2) model: Sx ~ bm + t + Db + bm:Db, random intercept,
## corCAR1(~t | ptID), CR2 + Satterthwaite (smoke test checks the
## agreement with fit_spec_s7(., 'E1cr2')). fit_analysis() falls back
## to a fit without corCAR1 on error; `car1` records which was used.
##
## Analysis variants (all on the same dataset within a setting):
##   p09 data (baseline row included by build_long()):
##     ta0        paper 09 'ta0': binary drug, + De if De varies, BL row
##     ta0.5      paper 09 'ta0.5': Dbc at t_a = 0.5, + De, BL row
##     ta0_noDe   binary, no De, BL row
##     ta0_noBL   binary, + De, BL row dropped
##     g6form     binary, no De, BL row dropped (paper 02's G6 form)
##   p02 data (no baseline row in prepare_long_data()):
##     g6         paper 02's G6 (binary, no De, no BL)
##     g6_De      binary, + De, no BL
##     g6_BL      binary, no De, BL row added (Sx = BL, t = 0, De = 0)
##     p09form    binary, + De, BL row added (paper 09's ta0 form)
## Control designs (CO, Hybrid) fit only ta0 and ta0.5.
##
## Environment: CHK_SETTING (p09 | p02), CHK_DESIGN (OLBDC | CO |
## Hybrid), CHK_N, CHK_REPS, CHK_CORES, CHK_OUT, CHK_SMOKE (1 = also
## check G6 agreement).

suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table)
  library(parallel)
  loadNamespace('clubSandwich')
})

setting <- Sys.getenv('CHK_SETTING', 'p09')
design <- Sys.getenv('CHK_DESIGN', 'OLBDC')
n_total <- as.integer(Sys.getenv('CHK_N', '70'))
n_reps <- as.integer(Sys.getenv('CHK_REPS', '20'))
n_cores <- as.integer(Sys.getenv('CHK_CORES', '8'))
out_dir <- Sys.getenv('CHK_OUT', 'analysis/data/cr2-olbdc-check')
smoke <- Sys.getenv('CHK_SMOKE', '0') == '1'

env09 <- new.env(parent = globalenv())
local({
  exprs <- parse('analysis/scripts/informative-dropout-by-design/20-paired-dropout.R',
                 keep.source = FALSE)
  stop_at <- which(vapply(exprs, function(e) {
    is.call(e) && identical(e[[1]], as.name('<-')) &&
      identical(e[[2]], as.name('units'))
  }, logical(1)))[1]
  for (e in exprs[seq_len(stop_at - 1L)]) eval(e, envir = env09)
})

env02 <- new.env(parent = globalenv())
sys.source('implementations/tidyverse/R/functions.R', envir = env02)
sys.source('analysis/scripts/carryover-sensitivity/simulation-core.R',
           envir = env02)

fit_analysis <- env09$fit_analysis

## Setting 'iso' (targeted isolation): paper 09's generator and
## analyses with some of paper 02's DGP settings swapped in, chosen by
## CHK_ISO: 'all' (component SDs tv 6 / br 8, cross-factor
## correlations 0.2 / 0.2, carryover half-life 1 week), 'sd', 'cf' or
## 'hl' (one group only). The Sigma cache is rebuilt to match.
iso <- Sys.getenv('CHK_ISO', 'all')
if (setting == 'iso') {
  if (iso %in% c('all', 'sd')) {
    env09$rp[cat == 'tv', sd := 6]
    env09$rp[cat == 'br', sd := 8]
  }
  mp_cf <- if (iso %in% c('all', 'cf')) c(0.2, 0.2) else c(0.1, 0.05)
  mp_hl <- if (iso %in% c('all', 'hl')) 1.0 else 0.5
  env09$make_mp <- function(beta_bm) {
    data.table(N = NA_integer_, c.bm = beta_bm,
               carryover_t1half = mp_hl,
               c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
               c.cf1t = mp_cf[1], c.cfct = mp_cf[2])
  }
  env09$sigma_cache <- lapply(env09$designs, function(td) {
    lapply(td$trialpaths, function(path) {
      buildSigma(env09$make_mp(0.45), env09$rp, env09$bp, path,
                 makePositiveDefinite = FALSE,
                 dgp_architecture = 'mean_moderation',
                 stopOnNonPD = TRUE)
    })
  })
  setting_tag <- paste0('iso_', iso)
} else {
  setting_tag <- setting
}

## 'iso' reuses the p09 seeds, so its datasets share the normal draws
## of the p09 cell (paired comparison).
setting_idx <- if (setting == 'p02') 2L else 1L
design_idx <- match(design, c('OLBDC', 'CO', 'Hybrid'))
seed_of <- function(rep) {
  70000000L + 10000000L * setting_idx + 1000000L * design_idx +
    100000L * as.integer(n_total == 70L) + rep
}

gen_p09 <- function(rep) {
  set.seed(seed_of(rep))
  td <- env09$designs[[design]]
  paths <- td$trialpaths
  n_pp <- env09$allocate_across_paths(n_total, length(paths))
  mp <- env09$make_mp(0)
  templates <- lapply(paths, env09$path_long_template)
  y_list <- bl_list <- bm_list <- fd_list <- vector('list', length(paths))
  for (g in seq_along(paths)) {
    mp$N <- n_pp[g]
    di <- generateData(mp, env09$rp, env09$bp, paths[[g]],
                       empirical = FALSE, makePositiveDefinite = FALSE,
                       cached_sigma = env09$sigma_cache[[design]][[g]],
                       dgp_architecture = 'mean_moderation',
                       stopOnNonPD = TRUE)
    tn <- paths[[g]]$timeptnames
    y_list[[g]] <- as.matrix(di[, ..tn])
    bl_list[[g]] <- di$BL
    bm_list[[g]] <- di$bm
    fd_list[[g]] <- rep(NA_integer_, n_pp[g])
  }
  long <- env09$build_long(y_list, bl_list, bm_list, fd_list, templates)
  use_de <- length(unique(paths[[length(paths)]]$e)) > 1
  list(long = long, use_de = use_de)
}

p02_cell <- list(design = design, N = n_total, c_bm = 0, t1half = 1.0,
                 carryover_form = 'exponential', weibull_shape = 1.0,
                 dgp_arch = 'mvn')

## Same model_param as simulate_cell_s7()/simulate_cell_combined()
## for this cell (rho 0.7, c.cf1t = c.cfct = 0.2, no confounds).
p02_model_param <- list(
  N = n_total, carryover_t1half = 1.0, carryover_form = 'exponential',
  weibull_shape = 1.0, dgp_architecture = 'mvn', c.bm = 0,
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.2,
  c.bm.bl = 0, c.bm.tv = 0, b.bm.tv = 0, b.bm.tv.trend = 0,
  sd_slope = 0, c.bm.pb = 0)

gen_p02 <- function(rep) {
  set.seed(seed_of(rep))
  ds <- env02$build_design_set(design)
  dat <- env02$generate_data_multi_path(
    p02_model_param, env02$default_resp_param(),
    env02$default_baseline_param(), ds, stop_on_nonpd = TRUE)
  dl <- env02$prepare_long_data(dat, ds, carryover_t1half = 1.0,
                                carryover_form = 'exponential',
                                weibull_shape = 1.0)
  dl <- as.data.table(dl)
  bl <- as.data.table(dat)[, .(ptID = paste(path, ptID, sep = '_'),
                               bm, Sx = BL)]
  bl[, `:=`(t = 0, Db = 0, Dbc = 0, L = 0, tsd = 0, De = 0)]
  list(dl = dl, bl = bl, dat = dat)
}

as_fit_long <- function(d) {
  d <- copy(d)
  d[, Db := as.logical(Db)]
  d[, .(ptID, bm, Sx, t, De, Db, tsd)]
}

row_of <- function(variant, r) {
  data.table(variant = variant, as.data.table(r))
}

## CHK_VARIANTS (comma-separated) restricts the variants fitted; the
## machine was shared with other jobs, so the full runs fit only the
## variants needed for the comparison.
variants <- strsplit(Sys.getenv('CHK_VARIANTS', ''), ',')[[1]]
want <- function(v) !length(variants) || v %in% variants

run_rep <- function(rep) {
  out <- list()
  if (setting %in% c('p09', 'iso')) {
    g <- gen_p09(rep)
    long <- g$long
    nobl <- long[t > 0]
    ude <- g$use_de
    if (want('ta0')) out$ta0 <- row_of('ta0', fit_analysis(long, 0, ude))
    if (want('ta0.5')) {
      out$ta05 <- row_of('ta0.5', fit_analysis(long, 0.5, ude))
    }
    if (want('ta0_noDe')) {
      out$nd <- row_of('ta0_noDe', fit_analysis(long, 0, FALSE))
    }
    if (want('ta0_noBL')) {
      out$nb <- row_of('ta0_noBL', fit_analysis(nobl, 0, ude))
    }
    if (want('g6form')) {
      out$g6 <- row_of('g6form', fit_analysis(nobl, 0, FALSE))
    }
    out <- rbindlist(out, fill = TRUE)
  } else {
    g <- gen_p02(rep)
    a <- as_fit_long(g$dl)
    b <- as_fit_long(rbind(g$dl, g$bl, use.names = TRUE, fill = TRUE))
    setorder(b, ptID, t)
    if (want('g6')) out$g6 <- row_of('g6', fit_analysis(a, 0, FALSE))
    if (want('g6_De')) out$gd <- row_of('g6_De', fit_analysis(a, 0, TRUE))
    if (want('g6_BL')) out$gb <- row_of('g6_BL', fit_analysis(b, 0, FALSE))
    if (want('p09form')) {
      out$pf <- row_of('p09form', fit_analysis(b, 0, TRUE))
    }
    out <- rbindlist(out, fill = TRUE)
    if (smoke) {
      ref <- env02$fit_spec_s7(g$dl, 'E1cr2')
      out[variant == 'g6', `:=`(ref_est = ref$estimate,
                                ref_se = ref$std_error,
                                ref_p = ref$p_value)]
    }
  }
  out[, `:=`(setting = setting_tag, design = design, N = n_total,
             rep = rep)]
  out
}

cat(sprintf('%s %s N=%d reps=%d cores=%d\n', setting, design, n_total,
            n_reps, n_cores))
t0 <- Sys.time()
## CHK_REP0 offsets the replicate indices (and so the seeds) for
## top-up runs; their output file carries a '-from<REP0+1>' suffix.
rep0 <- as.integer(Sys.getenv('CHK_REP0', '0'))
res <- mclapply(rep0 + seq_len(n_reps), function(r) {
  tryCatch(run_rep(r), error = function(e) {
    data.table(setting = setting, design = design, N = n_total, rep = r,
               error = conditionMessage(e))
  })
}, mc.cores = n_cores, mc.preschedule = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))
res <- rbindlist(res, fill = TRUE)
cat(sprintf('Elapsed %.1f s (%.2f s per dataset per core)\n', elapsed,
            elapsed * n_cores / n_reps))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fn <- file.path(out_dir, sprintf('%s%s-%s-N%d%s.rds',
                                 if (smoke) 'smoke-' else '', setting_tag,
                                 design, n_total,
                                 if (rep0 > 0) sprintf('-from%d', rep0 + 1L)
                                 else ''))
saveRDS(res, fn)
attr_meta <- list(elapsed = elapsed, n_reps = n_reps, n_cores = n_cores,
                  finished = Sys.time(), session = sessionInfo())
saveRDS(attr_meta, sub('\\.rds$', '-meta.rds', fn))

if ('error' %in% names(res)) print(res[!is.na(error), .N, by = error])
print(res[, .(n = sum(fit_ok, na.rm = TRUE),
              emp_sd = sd(est, na.rm = TRUE),
              se_mb = mean(se_mb, na.rm = TRUE),
              se_cr2 = mean(se_cr2, na.rm = TRUE),
              rej_mb = mean(p_mb < 0.05, na.rm = TRUE),
              rej_cr2 = mean(p_cr2 < 0.05, na.rm = TRUE),
              df_cr2 = median(df_cr2, na.rm = TRUE)),
          by = variant])
if (smoke && setting == 'p02') {
  print(res[variant == 'g6', .(max_abs_dest = max(abs(est - ref_est)),
                               max_abs_dse = max(abs(se_cr2 - ref_se)),
                               max_abs_dp = max(abs(p_cr2 - ref_p)))])
}
