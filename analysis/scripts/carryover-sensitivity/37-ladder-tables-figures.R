## analysis/scripts/carryover-sensitivity/37-ladder-tables-figures.R
##
## Analysis-ladder tables and figures for paper 02, built only from
## stored per-replicate and summary data (no simulation). Outputs are
## described in 37-ladder-README.md.
##
## Run from the repository root:
##   Rscript analysis/scripts/carryover-sensitivity/37-ladder-tables-figures.R
##
## Environment overrides (used for testing the isolation branch
## without touching the repository):
##   LADDER_ISOLATION_DIR  directory holding the e0-isolation all.rds
##   LADDER_OUT_DATA       output directory for the CSV files
##   LADDER_OUT_FIG        output directory for the figures

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})

data_dir <- here::here('analysis/data/carryover-sensitivity')
iso_dir <- Sys.getenv('LADDER_ISOLATION_DIR',
  file.path(data_dir, 'e0-isolation'))
out_data <- Sys.getenv('LADDER_OUT_DATA', file.path(data_dir, 'ladder'))
out_fig <- Sys.getenv('LADDER_OUT_FIG', here::here('analysis/figures'))
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

alpha_level <- 0.05
band <- c(0.04, 0.06)
iso_present <- file.exists(file.path(iso_dir, 'all.rds'))
message('Isolation data present: ', iso_present)

## Analyses ---------------------------------------------------------

analysis_meta <- tibble::tribble(
  ~analysis, ~label, ~rung_type, ~inference,
  'E0', 'E0', 'original', 'model-based',
  'E0n', 'E0n', 'original', 'model-based',
  'E0nc', 'E0nc', 'original', 'model-based',
  'E1nc', 'E1nc', 'original', 'model-based',
  'G1', 'G1', 'binary', 'model-based',
  'G6', 'G6', 'binary', 'CR2',
  'G3', 'G3', 'exposure-weighted', 'model-based',
  'G8', 'G8', 'exposure-weighted', 'CR2',
  'G4', 'G4', 'selected', 'model-based',
  'G9', 'G9', 'selected', 'CR2',
  'G9-null', 'G9-null', 'alternative', 'CR2',
  'G8-0.25', 'G8 at 0.25 wk', 'alternative', 'CR2',
  'G8-1', 'G8 at 1 wk', 'alternative', 'CR2',
  'F1', 'F1', 'alternative', 'CR2',
  'F1A0B', 'F1A0B', 'alternative', 'bootstrap-calibrated CR2'
)

analysis_cols <- c(
  E0 = '#4D4D4D', E0n = '#737373', E0nc = '#8C8C8C', E1nc = '#A6A6A6',
  G1 = '#56B4E9', G6 = '#0072B2',
  G3 = '#7FCDB0', G8 = '#009E73',
  G4 = '#E69F00', G9 = '#D55E00',
  `G9-null` = '#CC79A7', `G8-0.25` = '#6A5ACD', `G8-1` = '#A89FD8',
  F1 = '#7B2D8E', F1A0B = '#B5679B')

analysis_shapes <- c(
  E0 = 15, E0n = 15, E0nc = 15, E1nc = 15,
  G1 = 1, G6 = 16, G3 = 1, G8 = 16, G4 = 1, G9 = 16,
  `G9-null` = 17, `G8-0.25` = 17, `G8-1` = 17, F1 = 17, F1A0B = 17)

recode_analysis <- function(method, h = NA_real_) {
  case_when(
    method == 'E0' ~ 'E0',
    method == 'E0n' ~ 'E0n',
    method == 'E0nc' ~ 'E0nc',
    method == 'E1nc' ~ 'E1nc',
    method == 'E1' ~ 'G1',
    method == 'E1cr2' ~ 'G6',
    method == 'E2' ~ 'G3',
    method == 'E2cr2' ~ 'G8',
    method == 'E7' ~ 'G4',
    method == 'E7cr2' ~ 'G9',
    method == 'E7null_cr2' ~ 'G9-null',
    method == 'E3fix_cr2' & abs(h - 0.25) < 1e-8 ~ 'G8-0.25',
    method == 'E3fix_cr2' & abs(h - 1) < 1e-8 ~ 'G8-1',
    method == 'F1' ~ 'F1',
    .default = NA_character_)
}

## Families ---------------------------------------------------------

fam_levels <- c('1', '2', '3', '4', '5', '6', '7', 'S')
fam_labels <- c(
  '1' = '1 Original grid', '2' = '2 Misspecified decay',
  '3' = '3 Autocorrelation', '4' = '4 Dropout',
  '5' = '5 Prognostic', '6' = '6 Large N',
  '7' = '7 Trend and slope', 'S' = 'S Referee Hybrid run')

## Data loading -----------------------------------------------------

as_reps <- function(d, cell_key) {
  d |>
    mutate(
      cell_key = cell_key,
      analysis = recode_analysis(method, if ('h' %in% names(d)) h else NA_real_),
      valid = (if ('converged' %in% names(d)) converged else TRUE) &
        !is.na(p_value),
      reject = if_else(valid, as.numeric(p_value < alpha_level), NA_real_)) |>
    filter(!is.na(analysis)) |>
    transmute(cell_key, rep, analysis, reject,
              std_error = if_else(valid, std_error, NA_real_),
              estimate = if_else(valid, estimate, NA_real_))
}

make_cells <- function(d) {
  d |>
    mutate(
      label = sprintf('%s | N=%d | cbm=%.2f | t=%g', scenario, as.integer(N),
                      c_bm, t_half),
      pooled = if_else(role == 'null', pooled_null, ref)) |>
    select(cell_key, family, scenario, role, ref, pooled, N, c_bm, t_half,
           label)
}

## Original grid (family 1) ------------------------------------------

t1_all <- readRDS(file.path(data_dir, 'tier1-g-vs-f-e0/all.rds'))
t1_sum <- readRDS(file.path(data_dir, 'tier1-g-vs-f-e0/summary.rds'))
design_lab <- c(CO = 'CO', Hybrid = 'Hybrid', OLBDC = 'OL+BDC')

t1_cells_raw <- t1_all |>
  distinct(cell_id, design, N, c_bm, t1half) |>
  mutate(
    scenario = design_lab[design],
    role = if_else(c_bm == 0, 'null', 'alt'),
    ref = c_bm == 0.45 & t1half == 1,
    pooled_null = TRUE,
    t_half = t1half)

t1_cells <- t1_cells_raw |>
  mutate(cell_key = paste0('t1|', cell_id), family = '1') |>
  make_cells()
t1_reps <- as_reps(t1_all, paste0('t1|', t1_all$cell_id))

iso_reps <- NULL
iso_cells <- NULL
if (iso_present) {
  iso_all <- readRDS(file.path(iso_dir, 'all.rds'))
  stopifnot(all(c('cell_id', 'rep', 'method', 'p_value') %in% names(iso_all)))
  iso_reps <- as_reps(iso_all, paste0('iso|', iso_all$cell_id))
  iso_cells <- t1_cells_raw |>
    filter(cell_id %in% unique(iso_all$cell_id)) |>
    mutate(cell_key = paste0('iso|', cell_id), family = '1iso') |>
    make_cells()
}

## Stage 4 scenario cells (family 2, stored g-vs-f-comparison) ---------

gvf_all <- readRDS(file.path(data_dir, 'g-vs-f-comparison/all.rds'))
gvf_sum <- readRDS(file.path(data_dir, 'g-vs-f-comparison/summary.rds'))
gvf_cells_raw <- tibble::tribble(
  ~cell, ~scenario, ~c_bm, ~t_half,
  'shapemis_null', 'Weibull k=2 (Stage 4)', 0, 0.5,
  'shapemis_alt', 'Weibull k=2 (Stage 4)', 0.25, 0.5,
  't1hifhi_null', 't1/2=2 wk (Stage 4)', 0, 2,
  't1hifhi_alt', 't1/2=2 wk (Stage 4)', 0.25, 2) |>
  mutate(cell_key = paste0('gvf|', cell), family = '2', N = 70,
         role = if_else(c_bm == 0, 'null', 'alt'),
         ref = role == 'alt', pooled_null = TRUE)
gvf_cells <- make_cells(gvf_cells_raw)
gvf_reps <- gvf_all |>
  filter(cell %in% gvf_cells_raw$cell) |>
  (\(d) as_reps(d, paste0('gvf|', d$cell)))()

## Supplementary study (families 2 to 6) --------------------------------

gf_sum <- readRDS(file.path(data_dir, 'g9-gap-fill/summary.rds'))
gf_cells_raw <- gf_sum |>
  distinct(job, sub, cell_id, N, c_bm, t1half, rho, weibull_shape,
           dropout_rate, dropout_mech, dropout_per, c_pb, c_tv, c_bl) |>
  mutate(path = case_when(c_pb > 0 ~ 'PB', c_tv > 0 ~ 'TV-const',
                          c_bl > 0 ~ 'BL', .default = 'none'),
         strength = pmax(c_pb, c_tv, c_bl))

gf_def <- gf_cells_raw |>
  mutate(
    family = case_when(
      sub == 'A1-dropout' & dropout_per == 'patient' ~ '4',
      sub %in% c('B1-shape', 'B2-halflife') ~ '2',
      sub == 'B4-rho' ~ '3',
      sub == 'C1-ladder' ~ '5',
      sub == 'C2-largeN' ~ '6',
      .default = NA_character_),
    scenario = case_when(
      family == '4' ~ sprintf('%s %d%%', dropout_mech,
                              as.integer(round(100 * dropout_rate))),
      sub == 'B1-shape' ~ sprintf('Weibull k=%g', weibull_shape),
      sub == 'B2-halflife' ~ sprintf('t1/2=%g wk', t1half),
      family == '3' ~ sprintf('rho=%g', rho),
      family == '5' ~ sprintf('%s %.2f (ladder)', path, strength),
      family == '6' ~ sprintf('%s %.2f (large N)', path, strength),
      .default = NA_character_),
    ref_cbm = case_when(family == '4' ~ 0.45, family == '2' ~ 0.25,
                        family == '3' ~ 0.30, .default = NA_real_),
    role = if_else(c_bm == 0, 'null', 'alt'),
    ref = role == 'alt' & !is.na(ref_cbm) & abs(c_bm - ref_cbm) < 1e-8,
    pooled_null = TRUE,
    t_half = t1half,
    cell_key = paste('gf', job, cell_id, sep = '|')) |>
  filter(!is.na(family))
gf_cells <- make_cells(gf_def)

gf_files <- list.files(file.path(data_dir, 'g9-gap-fill'),
                       pattern = '^all\\.rds$', recursive = TRUE,
                       full.names = TRUE)
gf_reps <- map(gf_files, \(f) {
  job <- as.integer(sub('^job-', '', basename(dirname(f))))
  d <- readRDS(f) |> mutate(job = job)
  as_reps(d, paste('gf', d$job, d$cell_id, sep = '|'))
}) |>
  list_rbind() |>
  filter(cell_key %in% gf_cells$cell_key)

## Prognostic stress test, strength 0.3 (family 5) -----------------------

j28_all <- readRDS(file.path(data_dir, 'g9-prognostic-stress/all.rds'))
j28_sum <- readRDS(file.path(data_dir, 'g9-prognostic-stress/summary.rds'))
path_lab <- c(progpb = 'PB', progtvconst = 'TV-const', progbl = 'BL')
j28_cells_raw <- j28_sum |>
  distinct(cell_id, pathway, strength, c_bm, N) |>
  filter(pathway %in% names(path_lab), abs(strength - 0.3) < 1e-8) |>
  mutate(cell_key = paste0('j28|', cell_id), family = '5',
         scenario = sprintf('%s %.2f', path_lab[pathway], strength),
         role = if_else(c_bm == 0, 'null', 'alt'),
         ref = role == 'alt' & abs(c_bm - 0.45) < 1e-8,
         pooled_null = TRUE, t_half = 0.5)
j28_cells <- make_cells(j28_cells_raw)
j28_reps <- j28_all |>
  filter(cell_id %in% j28_cells_raw$cell_id) |>
  (\(d) as_reps(d, paste0('j28|', d$cell_id)))()

## Referee-variant blocks (family 7 and the side-analysis family S) --------

rv_dir <- file.path(data_dir, 'referee-variants')
rv_load <- function(block) {
  s <- readRDS(file.path(rv_dir, paste0('block-', block), 'summary.rds'))
  a <- readRDS(file.path(rv_dir, paste0('block-', block), 'all.rds'))
  cells <- s |>
    distinct(cell_id, N, c_bm, t1half, c_pb, b_tv_trend, sd_slope) |>
    mutate(
      cell_key = paste0('rv', block, '|', cell_id),
      family = if_else(block %in% c('T', 'S'), '7', 'S'),
      scenario = switch(block,
        N = sprintf('Hybrid null, cpb=%.1f', c_pb),
        A = rep('Hybrid alt', n()),
        T = rep('Trend 0.3', n()),
        S = rep('Patient slope 0.3', n())),
      role = if_else(c_bm == 0, 'null', 'alt'),
      ref = role == 'alt' & abs(c_bm - 0.45) < 1e-8,
      pooled_null = if (block == 'N') c_pb == 0 else TRUE,
      t_half = t1half)
  list(cells = make_cells(cells),
       reps = as_reps(a, paste0('rv', block, '|', a$cell_id)))
}
rv <- map(c('N', 'A', 'T', 'S'), rv_load)
rv_cells <- list_rbind(map(rv, 'cells'))
rv_reps <- list_rbind(map(rv, 'reps'))

## Assemble ---------------------------------------------------------

cells <- bind_rows(t1_cells, gvf_cells, gf_cells, j28_cells, rv_cells,
                   iso_cells)
reps <- bind_rows(t1_reps, gvf_reps, gf_reps, j28_reps, rv_reps, iso_reps)

stopifnot(
  'duplicate (cell, rep, analysis) rows' =
    !anyDuplicated(reps[, c('cell_key', 'rep', 'analysis')]),
  'replicate rows without a registered cell' =
    all(reps$cell_key %in% cells$cell_key))

cell_stats_reps <- reps |>
  filter(!is.na(reject)) |>
  group_by(cell_key, analysis) |>
  summarise(n = n(), k = sum(reject), kappa = mean(std_error) / sd(estimate),
            .groups = 'drop')

f1a0b_sum <- bind_rows(
  t1_sum |> transmute(cell_key = paste0('t1|', cell_id), method, reject, n),
  gvf_sum |> transmute(cell_key = paste0('gvf|', cell), method, reject, n),
  gf_sum |> transmute(cell_key = paste('gf', job, cell_id, sep = '|'),
                      method, reject, n),
  j28_sum |> transmute(cell_key = paste0('j28|', cell_id), method, reject, n)
) |>
  filter(method == 'F1A0B', n > 0, !is.na(reject),
         cell_key %in% cells$cell_key) |>
  transmute(cell_key, analysis = 'F1A0B', n = as.integer(n),
            k = round(reject * n), kappa = NA_real_)

cell_stats <- bind_rows(cell_stats_reps, f1a0b_sum) |>
  mutate(value = k / n, mcse = sqrt(value * (1 - value) / n)) |>
  inner_join(cells, by = 'cell_key')

pool_stats <- function(cs, group_vars = c('family', 'role', 'analysis')) {
  cs |>
    filter(pooled) |>
    group_by(across(all_of(group_vars))) |>
    summarise(n = sum(n), k = sum(k), n_cells = n(), .groups = 'drop') |>
    mutate(value = k / n, mcse = sqrt(value * (1 - value) / n))
}

pooled <- pool_stats(cell_stats)

kappa_fam <- cell_stats |>
  filter(!is.na(kappa)) |>
  group_by(family, analysis) |>
  summarise(kappa = mean(kappa), n_cells = n(), .groups = 'drop')
kappa_fam_role <- cell_stats |>
  filter(!is.na(kappa)) |>
  group_by(family, role, analysis) |>
  summarise(kappa = mean(kappa), n_cells = n(), .groups = 'drop')

## Isolation family: report the new analyses under family 1 ---------------

pooled_iso <- pooled |>
  filter(family == '1iso', analysis %in% c('E0n', 'E0nc', 'E1nc')) |>
  mutate(family = '1')
pooled <- pooled |>
  filter(family != '1iso') |>
  bind_rows(pooled_iso)
kappa_fam <- kappa_fam |>
  mutate(family = if_else(family == '1iso', '1', family)) |>
  distinct(family, analysis, .keep_all = TRUE)
kappa_fam_role <- kappa_fam_role |>
  mutate(family = if_else(family == '1iso', '1', family)) |>
  distinct(family, role, analysis, .keep_all = TRUE)

## Paired differences -----------------------------------------------------

wide <- reps |>
  pivot_wider(id_cols = c(cell_key, rep), names_from = analysis,
              values_from = reject)

paired_cells <- function(from, to) {
  if (!all(c(from, to) %in% names(wide))) {
    return(NULL)
  }
  wide |>
    filter(!is.na(.data[[from]]), !is.na(.data[[to]])) |>
    mutate(d = .data[[to]] - .data[[from]]) |>
    group_by(cell_key) |>
    summarise(n = n(), diff = mean(d), se = sd(d) / sqrt(n),
              .groups = 'drop') |>
    filter(n >= 2) |>
    inner_join(cells, by = 'cell_key')
}

paired_family <- function(from, to, family_override = NULL) {
  pc <- paired_cells(from, to)
  if (is.null(pc) || nrow(pc) == 0) {
    return(NULL)
  }
  pc |>
    filter(pooled) |>
    group_by(family, role) |>
    summarise(
      diff = sum(n * diff) / sum(n),
      se = sum(n * se) / sum(n),
      n = sum(n), n_cells = n(), .groups = 'drop') |>
    mutate(from = from, to = to)
}

steps <- tibble::tribble(
  ~step_id, ~from, ~to, ~kind, ~what,
  'E0>G1', 'E0', 'G1', 'core',
  'Software and df method, residual correlation, expectancy term and patient dropping change together (three respects at once)',
  'G1>G6', 'G1', 'G6', 'core',
  'CR2 cluster-robust standard errors with Satterthwaite df replace the model-based standard error',
  'G6>G8', 'G6', 'G8', 'core',
  'Binary drug indicator replaced by the exposure-weighted predictor at the true half-life',
  'G8>G9', 'G8', 'G9', 'core',
  'True half-life replaced by the half-life chosen by AIC over the candidate grid',
  'G9>G9-null', 'G9', 'G9-null', 'side',
  'AIC selection on the model without the interaction term, then refit with it',
  'G9>F1', 'G9', 'F1', 'side',
  'Selected half-life replaced by the fixed 0.2 week policy value, and the design-expectancy terms added',
  'G1>G3', 'G1', 'G3', 'optional',
  'Model-based: binary indicator replaced by the exposure-weighted predictor',
  'G3>G8', 'G3', 'G8', 'optional',
  'Exposure-weighted predictor: model-based replaced by CR2 inference',
  'G3>G4', 'G3', 'G4', 'optional',
  'Model-based: true half-life replaced by the AIC-selected half-life',
  'G4>G9', 'G4', 'G9', 'optional',
  'AIC-selected half-life: model-based replaced by CR2 inference',
  'G8>G8-0.25', 'G8', 'G8-0.25', 'side',
  'True half-life replaced by a fixed analysis half-life of 0.25 week',
  'G8>G8-1', 'G8', 'G8-1', 'side',
  'True half-life replaced by a fixed analysis half-life of 1 week')

iso_steps <- tibble::tribble(
  ~step_id, ~from, ~to, ~kind, ~what,
  'E0>E0n', 'E0', 'E0n', 'isolation',
  'Software and degrees-of-freedom method only (isolation run)',
  'E0n>E0nc', 'E0n', 'E0nc', 'isolation',
  'Residual correlation (corCAR1) only (isolation run)',
  'E0n>E1nc', 'E0n', 'E1nc', 'isolation',
  'Expectancy term and patient dropping removed (isolation run)',
  'E1nc>G1', 'E1nc', 'G1', 'isolation',
  'Residual correlation added to the binary analysis without the original expectancy term (isolation run)')

step_results <- steps |>
  mutate(res = map2(from, to, paired_family)) |>
  unnest(res, names_sep = '_') |>
  transmute(step_id, kind, what, family = res_family, role = res_role,
            from = res_from, to = res_to, diff = res_diff, se = res_se,
            n = res_n, n_cells = res_n_cells)

if (iso_present) {
  iso_res <- iso_steps |>
    mutate(res = map2(from, to, \(a, b) {
      pc <- paired_cells(a, b)
      if (is.null(pc)) {
        return(NULL)
      }
      pc |>
        filter(family == '1iso', pooled) |>
        group_by(family, role) |>
        summarise(diff = sum(n * diff) / sum(n), se = sum(n * se) / sum(n),
                  n = sum(n), n_cells = n(), .groups = 'drop')
    })) |>
    unnest(res, names_sep = '_') |>
    transmute(step_id, kind, what, family = '1', role = res_role,
              from, to, diff = res_diff, se = res_se, n = res_n,
              n_cells = res_n_cells)
  step_results <- bind_rows(step_results, iso_res)
}

step_results <- step_results |>
  mutate(quantity = if_else(role == 'null', 'size', 'power'))

## Coverage (stored in the referee reanalysis) ------------------------------

cov_raw <- read.csv(file.path(data_dir, 'referee-reanalysis/m3-coverage.csv'),
                    stringsAsFactors = FALSE) |>
  filter(level == 'cell')
coverage_fam <- cov_raw |>
  mutate(family = case_when(
    sub == 'A1-dropout' & dropout_per == 'patient' ~ '4',
    sub %in% c('B1-shape', 'B2-halflife') ~ '2',
    sub == 'B4-rho' ~ '3',
    sub == 'C1-ladder' ~ '5',
    .default = NA_character_),
    analysis = method) |>
  filter(!is.na(family), analysis %in% c('G8', 'G9')) |>
  group_by(family, analysis) |>
  summarise(lower = mean(coverage_lower), upper = mean(coverage_upper),
            n_cells = n(), .groups = 'drop') |>
  mutate(value = (lower + upper) / 2)

## Tidy long file -------------------------------------------------------

long_cell <- cell_stats |>
  filter(family != '1iso' | analysis %in% c('E0n', 'E0nc', 'E1nc')) |>
  mutate(family = if_else(family == '1iso', '1', family)) |>
  transmute(
    family, cell_group = label, scenario, level = 'cell', analysis,
    quantity = if_else(role == 'null', 'size', 'power'),
    value, mcse, n, from = NA_character_, reference = ref,
    in_family_pool = pooled)

long_pooled <- pooled |>
  transmute(
    family, cell_group = 'ALL', scenario = NA_character_, level = 'family',
    analysis, quantity = if_else(role == 'null', 'size', 'power'),
    value, mcse, n, from = NA_character_, reference = role == 'alt',
    in_family_pool = TRUE)

long_kappa <- kappa_fam_role |>
  transmute(
    family, cell_group = 'ALL', scenario = NA_character_, level = 'family',
    analysis, quantity = paste0('kappa_', role), value = kappa,
    mcse = NA_real_, n = n_cells, from = NA_character_, reference = NA,
    in_family_pool = TRUE) |>
  bind_rows(
    kappa_fam |>
      transmute(family, cell_group = 'ALL', scenario = NA_character_,
                level = 'family', analysis, quantity = 'kappa',
                value = kappa, mcse = NA_real_, n = n_cells,
                from = NA_character_, reference = NA, in_family_pool = TRUE))

long_diff <- step_results |>
  transmute(
    family, cell_group = 'ALL', scenario = NA_character_, level = 'family',
    analysis = to, quantity = paste0(quantity, '_diff'), value = diff,
    mcse = se, n, from, reference = role == 'alt', in_family_pool = TRUE)

long_cov <- coverage_fam |>
  pivot_longer(c(value, lower, upper), names_to = 'part', values_to = 'v') |>
  transmute(
    family, cell_group = 'ALL', scenario = NA_character_, level = 'family',
    analysis,
    quantity = if_else(part == 'value', 'coverage', paste0('coverage_', part)),
    value = v, mcse = NA_real_, n = n_cells, from = NA_character_,
    reference = TRUE, in_family_pool = TRUE)

ladder_long <- bind_rows(long_cell, long_pooled, long_kappa, long_diff,
                         long_cov) |>
  arrange(match(family, fam_levels), level, quantity, analysis, cell_group)

readr::write_csv(ladder_long, file.path(out_data, 'ladder-long.csv'))

## Table A ----------------------------------------------------------------

fmt_pm <- function(d, se) sprintf('%+.3f (%.3f)', d, se)

table_a_numeric <- step_results |>
  mutate(family = factor(family, levels = fam_levels)) |>
  select(step_id, kind, from, to, family, quantity, diff, se, n) |>
  pivot_wider(names_from = quantity, values_from = c(diff, se, n),
              names_glue = '{quantity}_{.value}')

steps_all <- bind_rows(steps, if (iso_present) iso_steps) |>
  mutate(order = row_number())
if (iso_present) {
  steps_all <- steps_all |>
    mutate(order = case_when(
      kind == 'isolation' ~ 1.5 + 0.1 * (row_number() - nrow(steps)),
      .default = order)) |>
    arrange(order)
} else {
  steps_all <- steps_all |>
    mutate(what = if_else(
      step_id == 'E0>G1',
      paste0(what, '; isolation data not available'), what))
}

table_a_cells <- step_results |>
  mutate(txt = fmt_pm(diff, se)) |>
  select(step_id, family, quantity, txt) |>
  pivot_wider(names_from = c(family, quantity), values_from = txt,
              names_glue = 'family{family}_{quantity}')

fam_cols <- as.vector(t(outer(paste0('family', fam_levels), c('size', 'power'),
                              paste, sep = '_')))
for (nm in setdiff(fam_cols, names(table_a_cells))) {
  table_a_cells[[nm]] <- NA_character_
}

table_a <- steps_all |>
  select(step_id, from, to, kind, what) |>
  left_join(table_a_cells, by = 'step_id') |>
  select(step_id, from, to, kind, what, all_of(fam_cols))
table_a[is.na(table_a)] <- '--'
readr::write_csv(table_a, file.path(out_data, 'ladder-table-A.csv'))
readr::write_csv(
  table_a_numeric |> arrange(match(step_id, steps_all$step_id)),
  file.path(out_data, 'ladder-table-A-numeric.csv'))

## Table B ----------------------------------------------------------------

pooled_size <- pooled |> filter(role == 'null')
pooled_power <- pooled |> filter(role == 'alt')

cell_size <- cell_stats |>
  filter(role == 'null', family != '1iso')
cell_alt <- cell_stats |>
  filter(role == 'alt', family != '1iso')

g9_cell_power <- cell_alt |>
  filter(analysis == 'G9') |>
  select(cell_key, g9 = value)

pool_g9_on <- function(fams_cells, role_) {
  cell_stats |>
    filter(analysis == 'G9', role == role_, pooled,
           cell_key %in% fams_cells) |>
    group_by(family) |>
    summarise(k = sum(k), n = sum(n), .groups = 'drop') |>
    mutate(g9 = k / n)
}

table_b_row <- function(a) {
  sz <- pooled_size |> filter(analysis == a)
  pw <- pooled_power |> filter(analysis == a)
  in_band <- sum(sz$value >= band[1] & sz$value <= band[2])

  cs_null <- cell_size |> filter(analysis == a)
  worst_hi <- cs_null |> slice_max(value, n = 1, with_ties = FALSE)
  worst_lo <- cs_null |> slice_min(value, n = 1, with_ties = FALSE)

  pw_on <- cell_stats |>
    filter(analysis == a, role == 'alt', pooled, family != '1iso')
  rel <- if (nrow(pw_on) > 0 && a != 'G9') {
    x <- pw_on |>
      group_by(family) |>
      summarise(kx = sum(k), nx = sum(n), keys = list(cell_key),
                .groups = 'drop')
    g9_on <- map_dfr(seq_len(nrow(x)), \(i) {
      pool_g9_on(x$keys[[i]], 'alt') |> filter(family == x$family[i])
    })
    if (nrow(g9_on) > 0) {
      tibble(fam = x$family, px = x$kx / x$nx) |>
        inner_join(g9_on |> select(fam = family, g9), by = 'fam') |>
        filter(g9 > 0) |>
        summarise(rel = 100 * mean(px / g9), loss = max(100 * (g9 - px)),
                  loss_family = fam[which.max(g9 - px)], m = n())
    } else {
      tibble(rel = NA_real_, loss = NA_real_, loss_family = NA_character_,
             m = 0L)
    }
  } else if (a == 'G9') {
    tibble(rel = 100, loss = 0, loss_family = NA_character_,
           m = nrow(pw_on |> distinct(family)))
  } else {
    tibble(rel = NA_real_, loss = NA_real_, loss_family = NA_character_,
           m = 0L)
  }

  cell_loss <- cell_alt |>
    filter(analysis == a) |>
    inner_join(g9_cell_power, by = 'cell_key') |>
    mutate(loss = 100 * (g9 - value)) |>
    slice_max(loss, n = 1, with_ties = FALSE)

  cov <- coverage_fam |> filter(analysis == a)

  tibble(
    analysis = a,
    families_size_available = nrow(sz),
    families_size_in_band = in_band,
    size_in_band = sprintf('%d of %d', in_band, nrow(sz)),
    max_size = if (nrow(worst_hi)) worst_hi$value else NA_real_,
    max_size_cell = if (nrow(worst_hi)) worst_hi$label else NA_character_,
    min_size = if (nrow(worst_lo)) worst_lo$value else NA_real_,
    min_size_cell = if (nrow(worst_lo)) worst_lo$label else NA_character_,
    families_power_available = rel$m,
    mean_power_pct_of_g9 = rel$rel,
    worst_family_power_loss_pp = rel$loss,
    worst_family = rel$loss_family,
    worst_cell_power_loss_pp = if (nrow(cell_loss)) cell_loss$loss else NA_real_,
    worst_cell = if (nrow(cell_loss)) cell_loss$label else NA_character_,
    coverage_mean = if (nrow(cov)) mean(cov$value) else NA_real_,
    coverage_families = nrow(cov))
}

table_b_analyses <- c('E0', 'G1', 'G6', 'G3', 'G8', 'G4', 'G9', 'G9-null',
                      'G8-1', 'G8-0.25', 'F1', 'F1A0B')
table_b <- map(table_b_analyses, table_b_row) |>
  list_rbind() |>
  left_join(analysis_meta |> select(analysis, label), by = 'analysis') |>
  relocate(label, .after = analysis)
readr::write_csv(table_b, file.path(out_data, 'ladder-table-B.csv'))

## Table C ----------------------------------------------------------------

table_c <- tibble::tribble(
  ~analysis, ~rung, ~formula, ~inference, ~assumes,
  'E0', '0',
  'Sx ~ bm + t + Db + bm:Db (+ De where expectancy varies); random intercept',
  'lme4 with lmerTest Satterthwaite df, no residual correlation',
  'Binary drug indicator; independent residuals given the random intercept; patients never on drug dropped',
  'G1', '1',
  'Sx ~ bm + t + Db + bm:Db; random intercept, corCAR1 residuals',
  'nlme model-based standard error',
  'Binary drug indicator; AR(1) residual correlation at the stated working structure',
  'G6', '2',
  'Sx ~ bm + t + Db + bm:Db; random intercept, corCAR1 residuals',
  'CR2 cluster-robust standard error, Satterthwaite df',
  'Binary drug indicator; working covariance may be wrong',
  'G3', 'optional (between 2 and 3)',
  'Sx ~ bm + t + Dbc + bm:Dbc; random intercept, corCAR1 residuals',
  'nlme model-based standard error',
  'Exponential carryover at the true half-life; working covariance correct',
  'G8', '3',
  'Sx ~ bm + t + Dbc + bm:Dbc; random intercept, corCAR1 residuals',
  'CR2 cluster-robust standard error, Satterthwaite df',
  'Exponential carryover at the true half-life',
  'G4', 'optional (between 3 and 4)',
  'G3 formula with the half-life chosen by AIC over {0.25, 0.5, 1, 2} weeks',
  'nlme model-based standard error, selection ignored',
  'Exponential carryover; true half-life near the candidate grid',
  'G9', '4',
  'G4 formula; half-life chosen by AIC over {0.25, 0.5, 1, 2} weeks',
  'CR2 cluster-robust standard error, Satterthwaite df, selection ignored',
  'Exponential carryover; true half-life near the candidate grid',
  'G9-null', 'side',
  'G9 formula; half-life chosen by AIC on the model without bm:Dbc, then refit with it',
  'CR2 cluster-robust standard error, Satterthwaite df',
  'As G9; selection does not use the interaction term',
  'G8-0.25', 'side',
  'G8 formula at a fixed analysis half-life of 0.25 week',
  'CR2 cluster-robust standard error, Satterthwaite df',
  'Exponential carryover; the fixed half-life need not equal the truth',
  'G8-1', 'side',
  'G8 formula at a fixed analysis half-life of 1 week',
  'CR2 cluster-robust standard error, Satterthwaite df',
  'Exponential carryover; the fixed half-life need not equal the truth',
  'F1', 'side',
  'Sx ~ bm + t + Dbc + bm:Dbc + De + bm:De; random intercept, corCAR1 residuals; half-life fixed at 0.2 week',
  'CR2 cluster-robust standard error, Satterthwaite df',
  'Exponential carryover; design expectancy De varies within patient (not identifiable under CO)',
  'F1A0B', 'side',
  'F1 formula; half-life chosen by AIC on the model without the interaction over {0, 0.2, 0.5, 1, 2} weeks',
  'Single-resample parametric bootstrap repeating the selection under the null',
  'Exponential carryover; design expectancy varies within patient (not identifiable under CO)')

if (iso_present) {
  table_c <- bind_rows(
    table_c,
    tibble::tribble(
      ~analysis, ~rung, ~formula, ~inference, ~assumes,
      'E0n', 'isolation',
      'E0 formula', 'Software and degrees-of-freedom method of the new code',
      'Isolation run: differs from E0 in software and df only',
      'E0nc', 'isolation',
      'E0 formula with corCAR1 residuals', 'As E0n',
      'Isolation run: adds residual correlation to E0n',
      'E1nc', 'isolation',
      'Sx ~ bm + t + Db + bm:Db; random intercept, no residual correlation',
      'As E0n',
      'Isolation run: removes the expectancy term and patient dropping from E0n')) |>
    arrange(match(analysis, c('E0', 'E0n', 'E0nc', 'E1nc', 'G1', 'G6', 'G3',
                              'G8', 'G4', 'G9', 'G9-null', 'G8-0.25',
                              'G8-1', 'F1', 'F1A0B')))
}
readr::write_csv(table_c, file.path(out_data, 'ladder-table-C.csv'))

## Figures ----------------------------------------------------------------

core_x <- c('E0', if (iso_present) c('E0n', 'E0nc', 'E1nc'),
            'G1', 'G6', 'G8', 'G9')
side_x <- c('G9-null', 'G8-0.25', 'G8-1', 'F1', 'F1A0B')
x_pos <- c(set_names(seq_along(core_x), core_x),
           set_names(length(core_x) + 1.8 + seq_along(side_x) - 1, side_x))
x_labs <- analysis_meta |>
  filter(analysis %in% names(x_pos)) |>
  mutate(label = c('G8-0.25' = 'G8 0.25 wk', 'G8-1' = 'G8 1 wk')[analysis] |>
           coalesce(label))

fig_theme <- theme_classic(base_size = 7.5) +
  theme(
    legend.position = 'none',
    strip.background = element_blank(),
    strip.text = element_text(size = 7, face = 'bold', hjust = 0),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 6.5),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    panel.spacing = unit(0.5, 'lines'))

x_scale <- scale_x_continuous(
  breaks = unname(x_pos[x_labs$analysis]),
  labels = x_labs$label,
  expand = expansion(add = 0.5))

prep <- function(d) {
  d |>
    filter(analysis %in% names(x_pos)) |>
    mutate(x = unname(x_pos[analysis]),
           family = factor(family, levels = fam_levels,
                           labels = unname(fam_labels[fam_levels])),
           analysis = factor(analysis, levels = names(analysis_cols)))
}

point_layers <- function(lo_mult = 2) {
  list(
    geom_errorbar(aes(ymin = pmax(value - lo_mult * mcse, -Inf),
                      ymax = value + lo_mult * mcse, colour = analysis),
                  width = 0, linewidth = 0.4),
    geom_point(aes(colour = analysis, shape = analysis), size = 1.7,
               stroke = 0.6),
    scale_colour_manual(values = analysis_cols),
    scale_shape_manual(values = analysis_shapes))
}

gap_line <- function(x) geom_vline(xintercept = x, linewidth = 0.25,
                                   colour = 'grey70', linetype = 'dotted')
gap_x <- length(core_x) + 0.9

fig1_data <- pooled |>
  filter(role == 'null') |>
  prep()

fig1 <- ggplot(fig1_data, aes(x, value)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = band[1], ymax = band[2],
           fill = 'grey85', alpha = 0.6) +
  geom_hline(yintercept = alpha_level, linewidth = 0.3, colour = 'grey40') +
  gap_line(gap_x) +
  point_layers() +
  x_scale +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 0.2, 0.05),
                     expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~family, ncol = 4) +
  labs(x = NULL, y = 'Null rejection rate') +
  fig_theme

for (ext in c('pdf', 'png')) {
  ggsave(file.path(out_fig, paste0('02-ladder-fig1-size.', ext)), fig1,
         width = 6.5, height = 4.6, dpi = 300)
}

fig2_power <- pooled |>
  filter(role == 'alt') |>
  mutate(panel = 'Power at reference alternative') |>
  prep()

step_prev <- tibble(
  analysis = c(if (iso_present) c('E0n', 'E0nc', 'E1nc'), 'G1', 'G6', 'G8',
               'G9', 'G9-null', 'G8-0.25', 'G8-1', 'F1'),
  from = c(if (iso_present) c('E0', 'E0n', 'E0n'),
           if (iso_present) 'E1nc' else 'E0', 'G1', 'G6', 'G8', 'G9', 'G8',
           'G8', 'G9'))

fig2_diff_src <- step_results |>
  filter(kind %in% c('core', 'side', 'isolation')) |>
  inner_join(step_prev, by = c('to' = 'analysis', 'from' = 'from')) |>
  filter(role == 'alt') |>
  transmute(family, analysis = to, value = diff, mcse = se,
            panel = 'Change from previous rung') |>
  prep()

empty_note <- tibble(
  family = factor(fam_labels['6'], levels = unname(fam_labels)),
  panel = factor('Power at reference alternative',
                 levels = c('Power at reference alternative',
                            'Change from previous rung')),
  x = mean(range(x_pos)), value = 0.5, text = 'no alternative cell')

fig2_data <- bind_rows(fig2_power, fig2_diff_src) |>
  mutate(panel = factor(panel, levels = c('Power at reference alternative',
                                          'Change from previous rung')))

fig2_block <- function(fams, show_strip_y = TRUE) {
  fl <- unname(fam_labels[fams])
  d <- fig2_data |> filter(family %in% fl) |>
    mutate(family = factor(family, levels = fl))
  en <- empty_note |> mutate(family = factor(as.character(family),
                                             levels = fl)) |>
    filter(!is.na(family))
  ggplot(d, aes(x, value)) +
    geom_hline(data = tibble(panel = factor('Change from previous rung',
                                            levels = levels(d$panel)),
                             yint = 0),
               aes(yintercept = yint), linewidth = 0.3, colour = 'grey40') +
    gap_line(gap_x) +
    point_layers() +
    geom_text(data = en, aes(x = x, y = value, label = text), size = 2.2,
              colour = 'grey30', inherit.aes = FALSE) +
    x_scale +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
    facet_grid(panel ~ family, scales = 'free_y', switch = 'y',
               labeller = labeller(panel = c(
                 'Power at reference alternative' = 'Power',
                 'Change from previous rung' = 'Change'))) +
    labs(x = NULL, y = NULL) +
    fig_theme +
    theme(strip.placement = 'outside',
          strip.text.y.left = element_text(angle = 90, face = 'plain'))
}

fig2 <- fig2_block(c('1', '2', '3', '4')) / fig2_block(c('5', '6', '7', 'S'))

for (ext in c('pdf', 'png')) {
  ggsave(file.path(out_fig, paste0('02-ladder-fig2-power.', ext)), fig2,
         width = 6.5, height = 7.4, dpi = 300)
}

## Figure 3: calibration and coverage --------------------------------------

fig3_analyses <- c('E0', 'G1', 'G6', 'G3', 'G8', 'G4', 'G9', 'F1')
fam_x <- c(fam_levels)

fig3a_data <- kappa_fam |>
  filter(analysis %in% fig3_analyses) |>
  mutate(family = factor(family, levels = fam_x),
         analysis = factor(analysis, levels = names(analysis_cols)))

dodge <- position_dodge(width = 0.8)

fig3a <- ggplot(fig3a_data, aes(family, kappa, colour = analysis,
                                shape = analysis)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = 'grey40') +
  geom_point(position = dodge, size = 1.5, stroke = 0.6) +
  scale_colour_manual(values = analysis_cols, breaks = fig3_analyses,
                      drop = FALSE) +
  scale_shape_manual(values = analysis_shapes, breaks = fig3_analyses,
                     drop = FALSE) +
  scale_x_discrete(labels = fam_x, drop = FALSE) +
  labs(x = 'Scenario family', y = 'Calibration ratio (model SE / empirical SD)',
       colour = NULL, shape = NULL) +
  fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7),
        legend.position = 'bottom', legend.key.size = unit(0.7, 'lines'),
        legend.text = element_text(size = 6.5))

fig3b_data <- coverage_fam |>
  mutate(family = factor(family, levels = fam_x),
         analysis = factor(analysis, levels = names(analysis_cols)))

fig3b <- ggplot(fig3b_data, aes(family, value, colour = analysis,
                                shape = analysis)) +
  geom_hline(yintercept = 0.95, linewidth = 0.3, colour = 'grey40') +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = dodge,
                width = 0, linewidth = 0.5) +
  geom_point(position = dodge, size = 1.5, stroke = 0.6) +
  scale_colour_manual(values = analysis_cols, breaks = fig3_analyses,
                      drop = FALSE) +
  scale_shape_manual(values = analysis_shapes, breaks = fig3_analyses,
                     drop = FALSE) +
  scale_x_discrete(labels = fam_x, drop = FALSE) +
  scale_y_continuous(limits = c(0.88, 0.97)) +
  labs(x = 'Scenario family', y = 'Coverage of the nominal 95% interval',
       colour = NULL, shape = NULL) +
  fig_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 7))

fig3 <- fig3a + (fig3b + theme(legend.position = 'none')) +
  plot_layout(widths = c(2.2, 1))

for (ext in c('pdf', 'png')) {
  ggsave(file.path(out_fig, paste0('02-ladder-fig3-calibration.', ext)), fig3,
         width = 6.5, height = 3.4, dpi = 300)
}

## Figure 4: size under the trend pathway and the patient-slope cell -------

fig4_cells <- cell_stats |>
  filter(family == '7', role == 'null') |>
  mutate(group = if_else(grepl('^Trend', scenario),
                         'Trend pathway (0.3 per week)',
                         'Patient slope SD 0.3'),
         analysis = factor(analysis, levels = names(analysis_cols)))

fig4_data <- fig4_cells |>
  mutate(group = factor(group, levels = c('Trend pathway (0.3 per week)',
                                          'Patient slope SD 0.3')),
         label_y = value + 2 * mcse)

fig4 <- ggplot(fig4_data, aes(analysis, value, fill = analysis,
                              alpha = group)) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = band[1], ymax = band[2],
           fill = 'grey85', alpha = 0.6) +
  geom_hline(yintercept = alpha_level, linewidth = 0.3, colour = 'grey40') +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = value - 2 * mcse, ymax = value + 2 * mcse,
                    group = group),
                position = position_dodge(width = 0.8), width = 0,
                linewidth = 0.4, alpha = 1) +
  geom_text(aes(y = label_y, label = sprintf('%.3f', value), group = group),
            position = position_dodge(width = 0.8), vjust = -0.5, size = 2.3,
            alpha = 1) +
  scale_fill_manual(values = analysis_cols, guide = 'none') +
  scale_alpha_manual(values = c(1, 0.45), name = NULL) +
  guides(alpha = guide_legend(override.aes = list(fill = 'grey30'))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = 'Null rejection rate') +
  theme_classic(base_size = 7.5) +
  theme(legend.position = 'bottom', axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3))

for (ext in c('pdf', 'png')) {
  ggsave(file.path(out_fig, paste0('02-ladder-fig4-trend.', ext)), fig4,
         width = 6.5, height = 3.2, dpi = 300)
}

## Console summary ----------------------------------------------------------

message('\nPooled size by family and rung')
pooled |>
  filter(role == 'null', analysis %in% c('G1', 'G6', 'G8', 'G9')) |>
  arrange(match(family, fam_levels), match(analysis, names(analysis_cols))) |>
  mutate(value = round(value, 4), mcse = round(mcse, 4)) |>
  select(family, analysis, n, n_cells, value, mcse) |>
  as.data.frame() |>
  print()

message('\nPooled power by family and rung')
pooled |>
  filter(role == 'alt', analysis %in% c('G1', 'G6', 'G8', 'G9')) |>
  arrange(match(family, fam_levels), match(analysis, names(analysis_cols))) |>
  mutate(value = round(value, 4), mcse = round(mcse, 4)) |>
  select(family, analysis, n, n_cells, value, mcse) |>
  as.data.frame() |>
  print()

message('\nFiles written to ', out_data, ' and ', out_fig)
