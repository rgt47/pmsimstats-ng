## analysis/scripts/dgp-combined/01-run-combined-factorial.R
##
## Factorial simulation driver for paper 11, the dual-channel
## architecture (data value `combined`). A 3 x 3 panel of channel
## weights c_bm_a (mean channel) x c_bm_b (covariance channel) in
## {0, 0.225, 0.45}, crossed with three designs (CO, Hybrid, OL+BDC)
## and the carryover half-lives in T1HALF, at a total of N = 70
## participants allocated across randomization paths.
##
## Rebuilt 2026-10-29 (paper 11 referee report 2026-10-28, M1, M2, M6,
## M8). Differences from the 2026-07-31 run, which is archived under
## analysis/data/dgp-combined/archive-2026-07-31/:
##
##   1. Extended DGP throughout: cross-component correlation c_1 = 0.2
##      at every lag (c.cfct = 0.2, was 0.1), carryover that decays
##      from the last on-drug mean (was the original compounding rule),
##      and a stop on a covariance matrix that is not positive definite
##      (stop_on_nonpd = TRUE; repair was enabled before but was never
##      triggered). All three come from the current
##      implementations/tidyverse/R/functions.R.
##   2. Pre-exposure guard on Dbc: 0 at off-drug visits before first
##      exposure (exposure_weight() in simulation-core.R). Without it a
##      placebo-first CO participant was coded exposed at every visit.
##   3. The baseline visit (t = 0, off drug, Dbc = 0, Sx = BL) enters
##      the analysis data, as in the package lme_analysis(), which
##      produced the single-channel results of paper 01.
##   4. Common random numbers: within a replicate (design, half-life,
##      replicate index) the same seed is set before every weight
##      cell's data generation. The multivariate-normal deviates and
##      the biomarker draws are therefore shared by all nine weight
##      cells, and every contrast between cells is paired.
##   5. Validation gate (implemented, not just promised): for the
##      first GATE_REPS replicates of every design and half-life, the
##      boundary cells are regenerated under the single-channel
##      architectures (`mean_moderation` with c.bm = c_bm_a, `mvn` with
##      c.bm = c_bm_b) from the same seed. The gate passes only if the
##      generated data and the fitted E2 estimates are identical; the
##      script stops otherwise.
##   6. Three analyses of every simulated trial:
##        E2     Sx ~ bm + t + Dbc + bm:Dbc, nlme::lme, random
##               intercept, corCAR1(~ t | ptID), REML, Dbc at the true
##               half-life, model-based standard error and t test.
##        E2cr2  the same fit, CR2 cluster-robust standard error and
##               Satterthwaite t test (clubSandwich); G8 in paper 02.
##        G9     half-life chosen by ML AIC over {0.25, 0.5, 1, 2}
##               weeks on the interaction model, REML refit, CR2 and
##               Satterthwaite (E7cr2 in simulation-core.R).
##      No fallback fit: a failed corCAR1 fit is recorded as not
##      converged rather than refitted without the correlation.
##
## Usage (thread variables are set by run-combined.sh):
##   bash analysis/scripts/dgp-combined/run-combined.sh [options]
##   options: --smoke             2 replicates, CO and Hybrid, t1/2 = 1
##            --reps N            replicates per design and half-life
##            --gate-reps N       replicates checked by the gate (50)
##            --workers N         parallel workers (6)
##            --methods e2,g9     analyses: e2 = E2 and E2cr2; g9 = G9
##            --cells all|g9      weight cells: all nine, or the null
##                                cell and the equal-total line
##                                (0,0), (0.45,0), (0.225,0.225),
##                                (0,0.45)
##            --t1half a,b,...    half-lives in weeks (0,0.2,0.5,1)
##            --out-tag TAG       suffix for the output file names
##   Because the seed depends only on (design, half-life, replicate),
##   a G9 pass run separately with --methods g9 --cells g9 analyzes
##   exactly the trials of the main pass's first N replicates.
##
## Runs made for the manuscript (2026-10-29):
##   main: --methods e2 --cells all --reps 1000 --t1half 0,0.2,0.5,1
##   G9:   --methods g9 --cells g9 --reps 200 --t1half 0,0.5,1
##         --gate-reps 0 --out-tag -g9
##
## Outputs (analysis/data/dgp-combined/):
##   01-combined-replicates<TAG>.rds  one row per replicate x cell x method
##   01-combined-summary<TAG>.rds     cell-level power, MCSE, estimates
##   01-combined-gate<TAG>.rds        validation-gate comparisons
##   01-combined-meta<TAG>.rds        commit, versions, settings, timing

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(parallel)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

local({
  e <- new.env()
  load(file.path(repo_root, 'data/extracted_rp.rda'), envir = e)
  load(file.path(repo_root, 'data/extracted_bp.rda'), envir = e)
  assign('pkg_rp', e$extracted_rp, envir = .GlobalEnv)
  assign('pkg_bp', e$extracted_bp, envir = .GlobalEnv)
})

## -----------------------------------------------------------------
## Arguments
## -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) && i < length(args)) args[i + 1] else default
}
smoke_mode <- '--smoke' %in% args
methods_in <- strsplit(arg_value('--methods', 'e2,g9'), ',')[[1]]
stopifnot(all(methods_in %in% c('e2', 'g9')))
run_e2     <- 'e2' %in% methods_in
run_g9     <- 'g9' %in% methods_in
cell_set   <- match.arg(arg_value('--cells', 'all'), c('all', 'g9'))
n_reps     <- as.integer(arg_value('--reps', if (smoke_mode) 2L else 1000L))
gate_reps  <- if (run_e2) as.integer(arg_value('--gate-reps',
                if (smoke_mode) 2L else 50L)) else 0L
n_workers  <- as.integer(arg_value('--workers', 6L))
t1half_set <- as.numeric(strsplit(
  arg_value('--t1half', if (smoke_mode) '1' else '0,0.2,0.5,1'), ',')[[1]])
out_tag    <- arg_value('--out-tag', if (smoke_mode) '-smoke' else '')
designs    <- if (smoke_mode) c('CO', 'Hybrid') else c('CO', 'Hybrid', 'OLBDC')

N_TOTAL   <- 70L
W_LEVELS  <- c(0, 0.225, 0.45)
weights   <- tidyr::expand_grid(c_bm_a = W_LEVELS, c_bm_b = W_LEVELS)
if (cell_set == 'g9') {
  weights <- dplyr::filter(weights,
    (c_bm_a == 0 & c_bm_b == 0) | abs(c_bm_a + c_bm_b - 0.45) < 1e-9)
}
stopifnot(nrow(weights) %in% c(4L, 9L))

units <- tidyr::expand_grid(design = designs, t1half = t1half_set,
                            rep = seq_len(n_reps))
cat(sprintf(paste0('Units (design x half-life x replicate): %d; ',
                   'weight cells per unit: %d; G9: %s; workers: %d\n'),
            nrow(units), nrow(weights), run_g9, n_workers))

## -----------------------------------------------------------------
## Data generation and analysis
## -----------------------------------------------------------------

unit_seed <- function(design, t1half, rep) {
  label_seed(sprintf('paper11|%s|%.3f|%d', design, t1half, rep))
}

model_param_for <- function(arch, t1half, a, b) {
  mp <- list(
    N = N_TOTAL, carryover_t1half = t1half,
    carryover_form = 'exponential', weibull_shape = 1,
    dgp_architecture = arch,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.2, c.cfct = 0.2
  )
  if (arch == 'combined') {
    mp$c.bm_a <- a
    mp$c.bm_b <- b
  } else if (arch == 'mean_moderation') {
    mp$c.bm <- a
  } else {
    mp$c.bm <- b
  }
  mp
}

generate_unit <- function(seed, arch, t1half, a, b, design_set) {
  set.seed(seed)
  generate_data_multi_path(model_param_for(arch, t1half, a, b),
                           pkg_rp, pkg_bp, design_set,
                           stop_on_nonpd = TRUE)
}

## Baseline visit as an analysis row: t = 0, off drug, never exposed,
## design expectancy 0, outcome the baseline level BL.
analysis_data <- function(dat, design_set, t1half) {
  dl <- prepare_long_data(dat, design_set, carryover_t1half = t1half)
  bl <- dat |>
    transmute(ptID = paste(path, ptID, sep = '_'), path, bm,
              timepoint_name = 'BL', Sx = BL, t = 0, Db = 0, Dbc = 0,
              L = 0, tsd = 0, De = 0)
  bind_rows(bl, dl) |> arrange(ptID, t)
}

na_rows <- function(methods) {
  tibble(method = methods, estimate = NA_real_, std_error = NA_real_,
         p_value = NA_real_, h = NA_real_, converged = FALSE)
}

fit_e2_pair <- function(dl, t1half) {
  fit <- tryCatch(
    nlme::lme(Sx ~ bm + t + Dbc + bm:Dbc, random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dl,
              control = nlme::lmeControl(
                opt = 'optim', maxIter = 200, msMaxIter = 200)),
    error = function(e) NULL)
  if (is.null(fit)) return(na_rows(c('E2', 'E2cr2')))
  cc <- summary(fit)$tTable
  e2 <- tibble(method = 'E2', estimate = cc['bm:Dbc', 'Value'],
               std_error = cc['bm:Dbc', 'Std.Error'],
               p_value = cc['bm:Dbc', 'p-value'], h = t1half,
               converged = TRUE)
  r <- cr2_extract(fit, dl, c('bm:Dbc', 'Dbc:bm'))
  cr2 <- if (is.null(r)) na_rows('E2cr2') else
    tibble(method = 'E2cr2', estimate = r$estimate,
           std_error = r$std_error, p_value = r$p_value, h = t1half,
           converged = TRUE)
  bind_rows(e2, cr2)
}

fit_g9 <- function(dl) {
  r <- fit_spec_s7(dl, 'E7cr2')
  tibble(method = 'G9', estimate = r$estimate, std_error = r$std_error,
         p_value = r$p_value, h = r$best_t_half, converged = r$converged)
}

methods_all <- c(if (run_e2) c('E2', 'E2cr2'), if (run_g9) 'G9')

run_unit <- function(u) {
  design_set <- build_design_set(u$design)
  seed <- unit_seed(u$design, u$t1half, u$rep)
  cells <- map_dfr(seq_len(nrow(weights)), function(k) {
    a <- weights$c_bm_a[k]
    b <- weights$c_bm_b[k]
    dat <- tryCatch(generate_unit(seed, 'combined', u$t1half, a, b,
                                  design_set),
                    error = function(e) NULL)
    res <- if (is.null(dat)) na_rows(methods_all) else {
      dl <- analysis_data(dat, design_set, u$t1half)
      bind_rows(if (run_e2) fit_e2_pair(dl, u$t1half),
                if (run_g9) fit_g9(dl))
    }
    mutate(res, c_bm_a = a, c_bm_b = b, .before = 1)
  })
  cells <- mutate(cells, design = u$design, t1half = u$t1half,
                  rep = u$rep, .before = 1)

  gate <- NULL
  if (u$rep <= gate_reps) {
    checks <- tibble(arch = c('mean_moderation', 'mean_moderation',
                              'mvn', 'mvn'),
                     c_bm_a = c(0.225, 0.45, 0, 0),
                     c_bm_b = c(0, 0, 0.225, 0.45))
    gate <- map_dfr(seq_len(nrow(checks)), function(k) {
      ck <- checks[k, ]
      d1 <- generate_unit(seed, 'combined', u$t1half, ck$c_bm_a,
                          ck$c_bm_b, design_set)
      d2 <- generate_unit(seed, ck$arch, u$t1half, ck$c_bm_a,
                          ck$c_bm_b, design_set)
      num <- names(d1)[vapply(d1, is.numeric, logical(1))]
      data_diff <- max(abs(as.matrix(d1[num]) - as.matrix(d2[num])))
      e1 <- fit_e2_pair(analysis_data(d1, design_set, u$t1half),
                        u$t1half)
      e2 <- fit_e2_pair(analysis_data(d2, design_set, u$t1half),
                        u$t1half)
      tibble(design = u$design, t1half = u$t1half, rep = u$rep,
             single_arch = ck$arch, c_bm_a = ck$c_bm_a,
             c_bm_b = ck$c_bm_b, max_abs_data_diff = data_diff,
             max_abs_est_diff = max(abs(e1$estimate - e2$estimate)),
             max_abs_p_diff = max(abs(e1$p_value - e2$p_value)))
    })
  }
  list(cells = cells, gate = gate)
}

## -----------------------------------------------------------------
## Run, in chunks, with checkpoints so an interrupted run resumes
## -----------------------------------------------------------------

out_dir <- file.path(repo_root, 'analysis/data/dgp-combined')
ckpt_dir <- file.path(out_dir, paste0('checkpoints', out_tag))
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

chunk_size <- 10L
chunk_id <- ceiling(seq_len(nrow(units)) / chunk_size)
chunks <- split(seq_len(nrow(units)), chunk_id)
settings_key <- sprintf('reps%d_gate%d_%s_%s_hl%s', n_reps, gate_reps,
                        paste(methods_in, collapse = '-'), cell_set,
                        paste(t1half_set, collapse = '-'))

t_start <- Sys.time()
done <- mclapply(seq_along(chunks), function(ci) {
  f <- file.path(ckpt_dir, sprintf('%s_chunk%05d.rds', settings_key, ci))
  if (file.exists(f)) return(f)
  out <- map(chunks[[ci]], function(i) run_unit(units[i, ]))
  saveRDS(list(cells = bind_rows(map(out, 'cells')),
               gate = bind_rows(map(out, 'gate'))), f)
  f
}, mc.cores = n_workers, mc.preschedule = FALSE)

failed <- vapply(done, function(x) inherits(x, 'try-error'), logical(1))
if (any(failed)) {
  print(done[failed][[1]])
  stop(sprintf('%d chunk(s) failed; rerun to resume from checkpoints.',
               sum(failed)))
}
elapsed <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
cat(sprintf('Simulation complete. Elapsed: %.1f s.\n', elapsed))

parts <- map(unlist(done), readRDS)
results <- bind_rows(map(parts, 'cells'))
gate <- bind_rows(map(parts, 'gate'))
stopifnot(nrow(results) == nrow(units) * nrow(weights) *
            length(methods_all))

## -----------------------------------------------------------------
## Validation gate
## -----------------------------------------------------------------

gate_pass <- if (gate_reps == 0) NA else nrow(gate) > 0 &&
  all(gate$max_abs_data_diff == 0) &&
  all(gate$max_abs_est_diff < 1e-10, na.rm = TRUE)
if (gate_reps > 0) {
  cat(sprintf(paste0('Validation gate: %d comparisons; max data diff ',
                     '%.3g; max estimate diff %.3g; %s\n'),
              nrow(gate), max(gate$max_abs_data_diff),
              max(gate$max_abs_est_diff, na.rm = TRUE),
              if (gate_pass) 'PASSED' else 'FAILED'))
} else {
  cat('Validation gate not run (--gate-reps 0 or no E2 analysis).\n')
}
if (isFALSE(gate_pass)) {
  saveRDS(gate, file.path(out_dir, paste0('01-combined-gate', out_tag,
                                          '.rds')))
  stop('Validation gate failed: the boundary cells do not reproduce ',
       'the single-channel architectures draw for draw.')
}

## -----------------------------------------------------------------
## Save
## -----------------------------------------------------------------

summary_dt <- results |>
  group_by(design, t1half, c_bm_a, c_bm_b, method) |>
  summarise(
    n_reps     = n(),
    n_conv     = sum(converged),
    power      = mean(p_value[converged] < 0.05),
    mcse_power = sqrt(power * (1 - power) / max(n_conv, 1L)),
    mean_est   = mean(estimate[converged]),
    sd_est     = sd(estimate[converged]),
    mean_se    = mean(std_error[converged]),
    kappa      = mean_se / sd_est,
    .groups = 'drop'
  ) |>
  mutate(N = N_TOTAL, .after = t1half)

git_head <- tryCatch(system2('git', c('-C', repo_root, 'rev-parse',
                                      'HEAD'), stdout = TRUE),
                     error = function(e) NA_character_)
git_dirty <- tryCatch(system2('git', c('-C', repo_root, 'status',
                                       '--porcelain', '--',
                                       'implementations/tidyverse/R',
                                       'analysis/scripts/carryover-sensitivity',
                                       'analysis/scripts/dgp-combined'),
                              stdout = TRUE),
                      error = function(e) NA_character_)
meta <- list(
  created = format(Sys.time(), '%Y-%m-%d %H:%M %Z'),
  git_head = git_head, git_dirty_relevant = git_dirty,
  r_version = R.version.string, platform = R.version$platform,
  blas = La_library(),
  pkg_versions = c(nlme = as.character(packageVersion('nlme')),
                   clubSandwich = as.character(packageVersion('clubSandwich'))),
  n_reps = n_reps, gate_reps = gate_reps, workers = n_workers,
  t1half = t1half_set, designs = designs, N = N_TOTAL,
  weights = W_LEVELS, methods = methods_all,
  dgp = list(rho = 0.7, c.cf1t = 0.2, c.cfct = 0.2,
             carryover = 'decays from the last on-drug mean',
             stop_on_nonpd = TRUE, resp_param = pkg_rp,
             baseline_param = pkg_bp),
  elapsed_secs = elapsed, gate_pass = gate_pass
)

saveRDS(results, file.path(out_dir, paste0('01-combined-replicates',
                                           out_tag, '.rds')))
saveRDS(summary_dt, file.path(out_dir, paste0('01-combined-summary',
                                              out_tag, '.rds')))
saveRDS(gate, file.path(out_dir, paste0('01-combined-gate', out_tag,
                                        '.rds')))
saveRDS(meta, file.path(out_dir, paste0('01-combined-meta', out_tag,
                                        '.rds')))
cat('Outputs written to', out_dir, '\n')

summary_dt |>
  dplyr::filter(method == methods_all[1]) |>
  dplyr::select(design, t1half, c_bm_a, c_bm_b, power, mcse_power,
                mean_est) |>
  print(n = 30)

invisible(NULL)
