## analysis/scripts/carryover-sensitivity/51-refit-failed-fits.R
##
## Diagnoses the failed fits of the sample-size runs
## (46-run-sample-size.R; referee report 2026-11-15, minor 12). Every
## lme fit in simulation-core.R uses lmeControl(opt = 'optim'). On a
## failing replicate (published OL+BDC, N = 740, replicate 1) the REML
## fit with optim stops with 'system is computationally singular',
## while the default nlminb optimizer fits the same data with ordinary
## estimates (phi 0.75). This script regenerates every trial with a
## failed A2 or A4 fit from its stored seed, and
##
##   1. reruns the stored fits (fit_spec_s7) as the pairing guard: a
##      fit that failed must fail again, and a fit that succeeded must
##      reproduce its stored estimate to GUARD_TOL;
##   2. refits the failed analysis with the same model and an nlminb
##      fallback whenever an optim fit fails;
##   3. writes the refits and, per affected cell, the rejection rate
##      as stored (failures dropped), with failures counted as
##      non-rejections, and with the refits.
##
## Nothing in the stored runs is changed.
##
## Usage (from the repository root; one call per parameter set):
##   Rscript analysis/scripts/carryover-sensitivity/51-refit-failed-fits.R \
##     --params simcore2026|published [--cores K]
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
source('implementations/tidyverse/R/functions.R')
source('analysis/scripts/carryover-sensitivity/simulation-core.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
CORES  <- as.integer(argv('--cores', '6'))
PARAMS <- Sys.getenv('PMSIM_PARAMS')
stopifnot(PARAMS %in% c('simcore2026', 'published'))
DATA_DIR <- 'analysis/data/carryover-sensitivity'
OUT_DIR  <- file.path(DATA_DIR, 'referee-reanalysis')
GUARD_TOL <- 1e-6
SPECS <- c('E1cr2', 'E7null_cr2')

## ---------------- trial regeneration (as 46-run-sample-size.R) ----------------
parse_cell <- function(cell_id) {
  m <- regmatches(cell_id, regexec(
    '^ss-([A-Za-z]+)-t([0-9.]+)-cbm([0-9.]+)-N([0-9]+)$', cell_id))[[1]]
  list(design = m[2], t1half = as.numeric(m[3]), c_bm = as.numeric(m[4]),
       N = as.numeric(m[5]), carryover_form = 'exponential',
       weibull_shape = 1, dgp_arch = 'mvn')
}

regen <- function(cell_id, r) {
  cell <- parse_cell(cell_id)
  set.seed(label_seed(sprintf('%s-rep-%d', cell_id, r)))
  design_set <- build_design_set(cell$design)
  mp <- list(N = cell$N, carryover_t1half = cell$t1half,
             carryover_form = cell$carryover_form,
             weibull_shape = cell$weibull_shape,
             dgp_architecture = cell$dgp_arch, c.bm = cell$c_bm,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2,
             c.cfct = 0.2, c.bm.bl = 0, c.bm.tv = 0, b.bm.tv = 0,
             b.bm.tv.trend = 0, sd_slope = 0, c.bm.pb = 0)
  dat <- generate_data_multi_path(mp, default_resp_param(),
                                  default_baseline_param(), design_set,
                                  stop_on_nonpd = TRUE)
  prepare_long_data(dat, design_set, carryover_t1half = cell$t1half,
                    carryover_form = cell$carryover_form,
                    weibull_shape = cell$weibull_shape,
                    analysis_t1half = cell$t1half,
                    analysis_form = cell$carryover_form,
                    analysis_shape = cell$weibull_shape)
}

## ---------------- fits with an nlminb fallback ----------------
fit_fb <- function(form, dl, method = 'REML') {
  one <- function(opt) tryCatch(
    nlme::lme(form, random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = dl, method = method,
              control = nlme::lmeControl(opt = opt, maxIter = 200,
                                         msMaxIter = 200)),
    error = function(e) NULL)
  f <- one('optim')
  if (!is.null(f)) return(list(fit = f, opt = 'optim'))
  f <- one('nlminb')
  if (!is.null(f)) return(list(fit = f, opt = 'nlminb'))
  NULL
}

refit_fb <- function(dl, spec) {
  na <- tibble(estimate = NA_real_, std_error = NA_real_,
               p_value = NA_real_, h = NA_real_, opt = NA_character_)
  if (spec == 'E1cr2') {
    f <- fit_fb(Sx ~ bm + t + Db + bm:Db, dl)
    if (is.null(f)) return(na)
    r <- cr2_extract(f$fit, dl, c('bm:Db', 'Db:bm'))
    if (is.null(r)) return(na)
    return(tibble(estimate = r$estimate, std_error = r$std_error,
                  p_value = r$p_value, h = NA_real_, opt = f$opt))
  }
  ## A4: ML selection over the grid on the model without the
  ## interaction, REML refit with it, CR2 (as select_half_life_aic).
  cands <- map(c(0.25, 0.5, 1, 2), function(hl) {
    d2 <- dl
    d2$Dbc <- exposure_weight(d2$Db, d2$tsd, hl)
    f <- fit_fb(E7_NULL_SELECTION, d2, 'ML')
    list(hl = hl, dl = d2, aic = if (is.null(f)) NA_real_ else AIC(f$fit),
         opt = if (is.null(f)) NA_character_ else f$opt)
  })
  aics <- map_dbl(cands, 'aic')
  if (all(is.na(aics))) return(na)
  best <- cands[[which.min(aics)]]
  f <- fit_fb(Sx ~ bm + t + Dbc + bm:Dbc, best$dl)
  if (is.null(f)) return(na)
  r <- cr2_extract(f$fit, best$dl, c('bm:Dbc', 'Dbc:bm'))
  if (is.null(r)) return(na)
  opts <- c(map_chr(cands, 'opt'), f$opt)
  tibble(estimate = r$estimate, std_error = r$std_error,
         p_value = r$p_value, h = best$hl,
         opt = if (any(opts == 'nlminb', na.rm = TRUE)) 'nlminb' else 'optim')
}

## ---------------- failed fits of this parameter set ----------------
dirs <- list.files(DATA_DIR, pattern = sprintf('^sample-size-%s-', PARAMS),
                   full.names = TRUE)
stored <- bind_rows(lapply(dirs, function(d)
  readRDS(file.path(d, 'all.rds')) |> mutate(run = basename(d))))
failed <- stored |>
  filter(!converged | is.na(p_value)) |>
  distinct(run, cell_id, rep)
cat(sprintf('%s: %d trials with a failed fit\n', PARAMS, nrow(failed)))

one_trial <- function(i) {
  f <- failed[i, ]
  dl <- regen(f$cell_id, f$rep)
  st <- stored |> filter(run == f$run, cell_id == f$cell_id, rep == f$rep)
  map_dfr(SPECS, function(s) {
    s_row <- st |> filter(method == s)
    stored_ok <- nrow(s_row) == 1 && s_row$converged && !is.na(s_row$p_value)
    again <- fit_spec_s7(dl, s)
    again_ok <- isTRUE(again$converged) && !is.na(again$p_value)
    guard <- if (stored_ok) {
      again_ok && abs(again$estimate - s_row$estimate) < GUARD_TOL
    } else {
      !again_ok
    }
    rf <- if (stored_ok) {
      tibble(estimate = s_row$estimate, std_error = s_row$std_error,
             p_value = s_row$p_value, h = NA_real_, opt = 'stored')
    } else {
      refit_fb(dl, s)
    }
    tibble(run = f$run, cell_id = f$cell_id, rep = f$rep, method = s,
           stored_ok = stored_ok, guard = guard) |>
      bind_cols(rf)
  })
}

t0 <- Sys.time()
refits <- bind_rows(mclapply(seq_len(nrow(failed)), one_trial,
                             mc.cores = CORES, mc.preschedule = FALSE))
cat(sprintf('%.1f min\n', as.numeric(Sys.time() - t0, units = 'mins')))
cat(sprintf('guard held in %d of %d fits\n', sum(refits$guard),
            nrow(refits)))
cat('refits of failed fits, by optimizer that succeeded:\n')
print(table(refits$opt[!refits$stored_ok], useNA = 'ifany'))
write.csv(refits, file.path(OUT_DIR,
                            sprintf('sample-size-refits-%s.csv', PARAMS)),
          row.names = FALSE)

## ---------------- rejection rates in the affected cells ----------------
cells <- stored |>
  semi_join(failed, by = c('run', 'cell_id')) |>
  left_join(refits |> filter(!stored_ok) |>
              select(run, cell_id, rep, method, p_refit = p_value),
            by = c('run', 'cell_id', 'rep', 'method')) |>
  mutate(p_fill = if_else(is.na(p_value), p_refit, p_value)) |>
  group_by(run, cell_id, method) |>
  summarise(n_reps = n(), n_failed = sum(!converged | is.na(p_value)),
            n_still_failed = sum(is.na(p_fill)),
            reject_stored = mean(p_value < 0.05, na.rm = TRUE),
            reject_fail_as_zero = mean(!is.na(p_value) & p_value < 0.05),
            reject_refit = mean(p_fill < 0.05, na.rm = TRUE),
            .groups = 'drop') |>
  mutate(params = PARAMS)
write.csv(cells, file.path(OUT_DIR,
                           sprintf('sample-size-refit-cells-%s.csv', PARAMS)),
          row.names = FALSE)
cat(sprintf('largest change in a cell rejection rate from the refits: %.4f\n',
            max(abs(cells$reject_refit - cells$reject_stored))))
