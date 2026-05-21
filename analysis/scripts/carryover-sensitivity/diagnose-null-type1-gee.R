## analysis/scripts/carryover-sensitivity/diagnose-null-type1-gee.R
##
## GEE verification for manuscript 10. The lme-plus-CR2 remedy adds
## a cluster-robust standard error to a conditional (mixed) model
## after the fact. A generalized estimating equations analysis is
## the marginal-model counterpart, with the sandwich variance built
## in. This script fits the GEE marginal model at the primary null
## cell and checks type-I calibration of the interaction test under
##   (a) the naive GEE sandwich standard error (geepack default),
##   (b) the Mancl-DeRouen bias-corrected sandwich (geesmv).
##
## The expectation: the naive sandwich is mildly anti-conservative
## at a moderate number of clusters, and the Mancl-DeRouen
## correction restores type-I to approximately 0.05.
##
## Primary null cell: c_bm = 0, exponential DGP, Hybrid design,
## N = 70, t1half = 1.0, rho = 0.7. The working correlation is
## AR(1); by the sandwich principle the inference is valid whatever
## the working correlation, which is the point of the comparison.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-null-type1-gee.R [--reps N]

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(ggplot2)
  library(geepack)
  library(geesmv)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

args <- commandArgs(trailingOnly = TRUE)
reps_idx <- which(args == '--reps')
N_REPS <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else 5000L
SEED <- 20260520L

cell <- list(
  c_bm = 0, carryover_form = 'exponential', weibull_shape = 1,
  t1half = 1.0, design = 'Hybrid', N = 70, rho = 0.7,
  dgp_arch = 'mvn')

resp_param     <- default_resp_param()
baseline_param <- default_baseline_param()
design_set     <- build_design_set(cell$design)

model_param <- list(
  N = cell$N, c.bm = cell$c_bm,
  carryover_t1half = cell$t1half,
  carryover_form   = cell$carryover_form,
  weibull_shape    = cell$weibull_shape,
  dgp_architecture = cell$dgp_arch,
  c.tv = cell$rho, c.pb = cell$rho, c.br = cell$rho,
  c.cf1t = 0.2, c.cfct = 0.1)

gen_long <- function() {
  dat <- tryCatch(
    generate_data_multi_path(model_param, resp_param,
                             baseline_param, design_set),
    error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  dl <- prepare_long_data(
    dat, design_set,
    carryover_t1half = cell$t1half,
    carryover_form   = cell$carryover_form,
    weibull_shape    = cell$weibull_shape)
  ## geeglm and geesmv require cluster rows contiguous and, for an
  ## AR(1) working correlation, time-ordered within cluster; geesmv
  ## additionally needs a numeric cluster id and a plain data frame.
  dl <- dl[order(dl$ptID, dl$t), ]
  dl$id <- as.integer(factor(dl$ptID))
  as.data.frame(dl)
}

spec_formula <- function(spec) switch(spec,
  A1 = Sx ~ bm + t + Db  + bm:Db,
  A2 = Sx ~ bm + t + Dbc + bm:Dbc,
  A3 = Sx ~ bm + t + Db  + bm:Db + L)

target_names <- function(spec)
  if (spec == 'A2') c('bm:Dbc', 'Dbc:bm') else c('bm:Db', 'Db:bm')

na_row <- function(spec) tibble(
  spec = spec, estimate = NA_real_, naive_se = NA_real_,
  naive_p = NA_real_, md_se = NA_real_, md_p = NA_real_,
  converged = FALSE)

fit_gee <- function(dl, spec) {
  form <- spec_formula(spec)

  fit <- tryCatch(
    geepack::geeglm(form, family = gaussian, data = dl,
                    id = id, corstr = 'ar1'),
    error = function(e) NULL)
  if (is.null(fit)) return(na_row(spec))

  sm <- summary(fit)$coefficients
  tgt <- intersect(target_names(spec), rownames(sm))
  if (length(tgt) == 0) return(na_row(spec))
  tgt <- tgt[1]
  est      <- sm[tgt, 'Estimate']
  naive_se <- sm[tgt, 'Std.err']
  naive_p  <- sm[tgt, 'Pr(>|W|)']

  ## Mancl-DeRouen bias-corrected sandwich (geesmv).
  md <- tryCatch(
    {
      utils::capture.output(
        m <- geesmv::GEE.var.md(form, id = 'id', family = gaussian,
                                data = dl, corstr = 'AR-M'))
      m
    },
    error = function(e) NULL)
  if (is.null(md)) {
    md_se <- NA_real_; md_p <- NA_real_
  } else {
    idx   <- match(tgt, names(coef(fit)))
    md_se <- sqrt(md$cov.beta[idx])
    n_id  <- length(unique(dl$ptID))
    df    <- n_id - length(coef(fit))
    md_p  <- 2 * pt(-abs(est / md_se), df)
  }

  tibble(spec = spec, estimate = est, naive_se = naive_se,
         naive_p = naive_p, md_se = md_se, md_p = md_p,
         converged = TRUE)
}

one_rep <- function(i) {
  dl <- gen_long()
  if (is.null(dl)) return(NULL)
  map_dfr(c('A1', 'A2', 'A3'), function(sp) fit_gee(dl, sp)) |>
    mutate(rep = i, .before = 1)
}

cat(sprintf('GEE verification: %s, N=%d, c_bm=0, exp DGP, %d reps\n',
            cell$design, cell$N, N_REPS))

plan(multicore, workers = max(1, parallel::detectCores() - 1))
t0 <- Sys.time()
res <- future_map_dfr(seq_len(N_REPS), one_rep,
                      .options = furrr_options(seed = SEED))
elapsed <- as.numeric(Sys.time() - t0, units = 'secs')
cat(sprintf('Completed in %.0f s.\n\n', elapsed))

summ <- res |>
  filter(converged, !is.na(md_p)) |>
  group_by(spec) |>
  summarise(
    n              = n(),
    emp_se         = sd(estimate),
    naive_se       = mean(naive_se),
    md_se          = mean(md_se),
    naive_se_ratio = mean(naive_se) / sd(estimate),
    md_se_ratio    = mean(md_se) / sd(estimate),
    type1_naive    = mean(naive_p < 0.05),
    type1_md       = mean(md_p < 0.05),
    .groups = 'drop')

mcse <- sqrt(0.05 * 0.95 / unique(summ$n)[1])

cat('=== GEE marginal model: naive vs Mancl-DeRouen sandwich ===\n')
cat('naive_se_ratio : geepack sandwich SE / empirical SE\n')
cat('md_se_ratio    : Mancl-DeRouen SE / empirical SE (expect ~1.00)\n')
cat('type1_naive    : type-I, naive GEE sandwich (z reference)\n')
cat('type1_md       : type-I, Mancl-DeRouen sandwich (t reference)\n')
cat(sprintf('MC SE on a type-I estimate: ~%.4f\n\n', mcse))
print(as.data.frame(summ), digits = 4)

verdict <- nrow(summ) == 3L && all(abs(summ$type1_md - 0.05) < 3 * mcse)
cat(sprintf('\nVERDICT: Mancl-DeRouen GEE calibration %s\n',
    if (verdict) 'CONFIRMED (type-I ~0.05).'
    else 'NOT confirmed - inspect table above.'))

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')

p_ecdf <- res |>
  filter(converged, !is.na(md_p)) |>
  select(spec, naive_p, md_p) |>
  pivot_longer(c(naive_p, md_p), names_to = 'source',
               values_to = 'p') |>
  mutate(source = recode(source,
    naive_p = 'Naive GEE sandwich',
    md_p = 'Mancl-DeRouen sandwich')) |>
  ggplot(aes(p, colour = source)) +
  stat_ecdf(geom = 'step', linewidth = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed',
              colour = 'grey50') +
  facet_wrap(~ spec) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c(
    'Naive GEE sandwich' = '#e31a1c',
    'Mancl-DeRouen sandwich' = '#1f78b4')) +
  labs(x = 'Null p-value', y = 'Empirical CDF', colour = NULL,
       title = 'GEE null p-values: naive vs bias-corrected sandwich',
       subtitle = sprintf(
         'Marginal model, Hybrid, N=70, c_bm=0, exp DGP, %d reps',
         N_REPS)) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = 'top')

ggsave(file.path(out_dir, 'diag-null-type1-gee-ecdf.pdf'),
       p_ecdf, width = 7.5, height = 3.4)

saveRDS(list(results = res, summary = summ, cell = cell,
             n_reps = N_REPS, seed = SEED, elapsed_secs = elapsed),
        file.path(out_dir, 'diag-null-type1-gee.rds'))

cat(sprintf('\nWrote %s\n',
    file.path(out_dir, 'diag-null-type1-gee-ecdf.pdf')))
cat(sprintf('Wrote %s\n',
    file.path(out_dir, 'diag-null-type1-gee.rds')))
