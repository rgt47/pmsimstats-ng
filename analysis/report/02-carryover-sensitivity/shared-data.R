## analysis/report/02-carryover-sensitivity/shared-data.R
##
## Data loading and helpers shared by report.Rmd, summary.Rmd and
## supplement.Rmd of paper 02. Reads stored data only. Each document
## sources it after loading tibble, dplyr, tidyr and knitr.


## Display names. Stored values stay E1/E2/E3/E7/E9 and their cr2
## variants everywhere; only the display layer maps them to G1-G9.
source(here::here(
  'analysis/scripts/carryover-sensitivity/spec-labels.R'))

## Data loaded once for the tables below.
grid <- readRDS(here::here('analysis/data/02-grid-summary.rds'))$summary
s27  <- readRDS(here::here(
  'analysis/data/carryover-sensitivity/tier1-g-vs-f-e0/summary.rds'))
s21  <- readRDS(here::here(
  'analysis/data/carryover-sensitivity/g-vs-f-comparison/summary.rds'))
sgf  <- readRDS(here::here(
  'analysis/data/carryover-sensitivity/g9-gap-fill/summary.rds'))

## G10 fitted afterward to the stored trials of families 1 to 5
## (41-run-g10-backfill.R). The backfill refit of G9 reproduces the
## stored G9 summaries, so its G10 rows sit beside the stored rows.
bf_dir <- here::here('analysis/data/carryover-sensitivity/g10-backfill')
bf_sum <- readRDS(file.path(bf_dir, 'summary.rds'))
bf_all <- readRDS(file.path(bf_dir, 'all.rds'))
bf_guard <- readRDS(file.path(bf_dir, 'guard.rds'))
bf_rej <- function(src, m, ...) {
  d <- bf_sum |> filter(source == src, method == m, ...)
  stopifnot(nrow(d) == 1)
  d$reject
}
bf_mean_est <- bf_all |>
  filter(method == 'E7null_cr2') |>
  group_by(source, cell_id) |>
  summarise(mean_est = mean(estimate, na.rm = TRUE), .groups = 'drop')
bf_gvf <- bf_sum |>
  filter(source == 'gvf', method == 'E7null_cr2') |>
  transmute(cell = cell_id, method, reject)
## One rejection rate of the scenario job (job 21), stored or backfilled.
g21 <- function(cell_id, m) s21$reject[s21$cell == cell_id & s21$method == m]
g21b <- function(cell_id, m) bf_gvf$reject[bf_gvf$cell == cell_id &
                                            bf_gvf$method == m]
sgf_g10 <- bf_sum |>
  filter(source == 'gf', method == 'E7null_cr2') |>
  select(cell_id, method, reject, n, nonconv, mean_se, emp_se, kappa) |>
  left_join(bf_mean_est |> filter(source == 'gf') |> select(-source),
            by = 'cell_id') |>
  left_join(sgf |> filter(method == 'E7cr2') |>
              select(-method, -reject, -n, -nonconv, -mean_est, -mean_se,
                     -emp_se, -kappa),
            by = 'cell_id')
stopifnot(!anyNA(sgf_g10$sub))
sgf <- bind_rows(sgf, sgf_g10)
## G10 (A4) on the job-27 trials (backfill source t1), added beside the
## stored analyses of the 30-cell grid. Replicates whose G9 refit failed
## the pairing guard are dropped, as in 37-ladder-tables-figures.R.
bf_ok <- bf_guard |> filter(match) |> select(source, cell_id, rep)
cell_meta27 <- s27 |> distinct(cell_id, panel, design, N, c_bm, t1half)
a27_g10 <- bf_all |>
  filter(source == 't1', method == 'E7null_cr2') |>
  semi_join(bf_ok, by = c('source', 'cell_id', 'rep')) |>
  select(cell_id, rep, method, h, estimate, std_error, p_value,
         converged, t) |>
  left_join(cell_meta27, by = 'cell_id') |>
  mutate(family = 'G')
s27_g10 <- a27_g10 |>
  filter(converged, !is.na(p_value)) |>
  group_by(cell_id, panel, design, N, c_bm, t1half, family, method) |>
  summarise(reject = mean(p_value < 0.05), n = n(), .groups = 'drop')
stopifnot(nrow(s27_g10) == 30, !anyNA(s27_g10$t1half))
s27 <- bind_rows(s27, s27_g10)
## One rejection rate from the gap-fill summary; the dots are filter
## conditions and must identify exactly one cell.
gf <- function(sub_, m, ...) {
  d <- sgf |> filter(sub == sub_, method == m, ...)
  stopifnot(nrow(d) == 1)
  d$reject
}

## One column of the 2,000-replicate size diagnostic (A2-size cells)
## for one method, and a formatted range (referee report 2026-11-08,
## m14: inline instead of hard-coded).
a2v <- function(m, col, n = c(35, 70)) {
  x <- sgf[[col]][sgf$sub == 'A2-size' & sgf$method == m & sgf$N %in% n]
  stopifnot(length(x) > 0)
  x
}
rngf <- function(x, d = 3) sprintf(paste0('%.', d, 'f to %.', d, 'f'),
                                   min(x), max(x))
## G1 and G3 (working-covariance standard error) in the rho = 0.9
## cells of the serial-correlation grid, at c_bm = cb.
rho9 <- function(col, cb) {
  x <- sgf[[col]][sgf$sub == 'B4-rho' & sgf$rho == 0.9 &
                    sgf$c_bm == cb & sgf$method %in% c('E1', 'E2')]
  stopifnot(length(x) == 4)
  x
}

## Method display labels shared by the tables.
m_lab <- c(E0 = 'E0', E1 = 'G1', E3 = 'G2', E2 = 'G3', E7 = 'G4',
           E9 = 'G5', E1cr2 = 'G6', E3cr2 = 'G7', E2cr2 = 'G8',
           E7cr2 = 'G9', E7null_cr2 = 'G10')
## Tables show the paper 02 short list (A0-A4) only.
focus <- c('E0', 'E1', 'E1cr2', 'E2cr2', 'E7null_cr2')
focus_cols <- a_order
design_lab <- c(CO = 'CO', Hybrid = 'Hybrid', OLBDC = 'OL+BDC')
fmt3 <- function(x) ifelse(is.na(x), '--', sprintf('%.3f', x))
## Reanalysis of the per-replicate output for the 2026-10-27 referee
## report (34-referee-reanalysis.R); the tables that use it are
## skipped if the files are absent.
rr_dir <- here::here('analysis/data/carryover-sensitivity/referee-reanalysis')
has_rr <- file.exists(file.path(rr_dir, 'm2-pooled-size.csv'))
g9_pooled <- if (has_rr) read.csv(file.path(rr_dir, 'm2-pooled-size.csv')) else NULL
## G9 null rejection by the half-life it selected (four A2 cells).
g9h <- if (has_rr) g9_pooled[g9_pooled$item == 'G9 conditional on selected h', ] else NULL
## The four referee-variant runs (36-run-referee-variants.R, simrun
## jobs 34 to 37 on zeus, 2026-10-27), one directory per block.
rv_dir <- here::here('analysis/data/carryover-sensitivity/referee-variants')
rv_read <- function(block, what) {
  f <- file.path(rv_dir, paste0('block-', block), paste0(what, '.rds'))
  if (file.exists(f)) readRDS(f) else NULL
}
has_rv <- !is.null(rv_read('N', 'summary'))
rv_lab <- c(E2cr2 = 'G8', E7null_cr2 = 'G10', E7cr2 = 'G9')
## E3fix_cr2 is the G8 formula at a fixed analysis half-life (0.25 or
## 1 week) instead of the true one (referee report 2026-11-08, M3).
rv_name <- function(m, h) ifelse(m == 'E3fix_cr2',
                                 sprintf('Fixed %g', h), rv_lab[m])
rv_ord <- c('G8', 'G10', 'G9', 'Fixed 0.25', 'Fixed 1')
## Paper 02 tables show the short list only: A3 (G8) and A4 (G10).
rv_tab_m <- c('E2cr2', 'E7null_cr2')
## Block N also holds null cells with a prognostic placebo pathway
## (c_pb = 0.3), which belong to the companion paper.
rvN <- rv_read('N', 'summary')
if (!is.null(rvN)) rvN <- rvN |> filter(c_pb == 0)
rvA <- rv_read('A', 'summary'); rvS <- rv_read('S', 'summary')
rvNp <- rv_read('N', 'paired')
if (!is.null(rvNp)) rvNp <- rvNp |> filter(cell_id %in% rvN$cell_id)
rvAp <- rv_read('A', 'paired')
rv_cell <- function(d, m, ...) {
  x <- d |> filter(method == m, ...)
  stopifnot(nrow(x) == 1)
  x
}

## Isolation block (38-run-e0-isolation.R, checked by
## 39-check-e0-isolation.R; simrun job 42). Paired step differences
## are computed from the per-replicate file.
iso_dir <- here::here('analysis/data/carryover-sensitivity/e0-isolation')
iso <- readRDS(file.path(iso_dir, 'summary.rds'))
iso_all <- readRDS(file.path(iso_dir, 'all.rds'))
stopifnot(nrow(iso) == 150, nrow(iso_all) == 150000)

## Analysis-ladder tables and figures (37-ladder-tables-figures.R).
## The script reads stored data only. This document reads its CSV
## files and recomputes nothing.
lad_dir <- here::here('analysis/data/carryover-sensitivity/ladder-predictive')
lad_long <- read.csv(file.path(lad_dir, 'ladder-long.csv'),
                     stringsAsFactors = FALSE)
lad_A  <- read.csv(file.path(lad_dir, 'ladder-table-A.csv'),
                   stringsAsFactors = FALSE, check.names = FALSE)
lad_An <- read.csv(file.path(lad_dir, 'ladder-table-A-numeric.csv'),
                   stringsAsFactors = FALSE)
lad_B  <- read.csv(file.path(lad_dir, 'ladder-table-B.csv'),
                   stringsAsFactors = FALSE)
lad_C  <- read.csv(file.path(lad_dir, 'ladder-table-C.csv'),
                   stringsAsFactors = FALSE)
lad_fam <- lad_long |> filter(level == 'family')
## One family-level value (size, power, kappa, coverage) of an
## analysis. The filter must identify exactly one row.
lv <- function(a, q, f = '1', w = 'value') {
  x <- lad_fam |> filter(analysis == a, quantity == q, family == f)
  stopifnot(nrow(x) == 1)
  x[[w]]
}
## One paired increment between adjacent analyses, from the numeric
## copy of Table A. `w` is size_diff, size_se, power_diff or power_se.
inc <- function(step, f = '1', w = 'power_diff') {
  x <- lad_An |> filter(step_id == step, family == f)
  stopifnot(nrow(x) == 1)
  x[[w]]
}
## The range over families of one paired increment, formatted.
inc_rng <- function(step, w = 'power_diff') {
  x <- (lad_An |> filter(step_id == step))[[w]]
  stopifnot(length(x) > 0, !anyNA(x))
  sprintf('%+.3f to %+.3f', min(x), max(x))
}
## Cell-level power of an analysis in family 1, pooled over the cells
## of one design (weighted by replicates).
cell_pw <- function(d, a) {
  x <- lad_long |> filter(level == 'cell', family == '1',
                          quantity == 'power', analysis == a,
                          sub(' .*', '', scenario) == d)
  stopifnot(nrow(x) > 0)
  weighted.mean(x$value, x$n)
}
fam_ids <- c('1', '2', '3', '4', '5', 'S')
stopifnot(setequal(setdiff(unique(lad_fam$family), '1iso'), fam_ids))
fam_size <- function(a) {
  x <- lad_fam |> filter(analysis == a, quantity == 'size')
  setNames(x$value, x$family)
}
f3 <- function(x) sprintf('%.3f', x)
fs3 <- function(x) sprintf('%+.3f', x)
lag_mean_g2 <- grid |>
  filter(dgp_arch == 'mvn', carryover_form %in% c('exponential', 'weibull'),
         !(carryover_form == 'weibull' & weibull_shape == 1),
         spec %in% c('E1', 'E3')) |>
  select(design, N, c_bm, t1half, carryover_form, weibull_shape, spec,
         power) |>
  pivot_wider(names_from = spec, values_from = power) |>
  summarise(m = mean(E3 - E1)) |> pull(m)
## Hybrid-only headline (referee report 2026-11-08, M7): G10 against
## G6, the simplest correctly sized analysis that needs no half-life,
## over the Hybrid alternative cells of families 1 to 4 (families 5
## and S have no G6). Same trials, so the cell differences are paired.
hy_cells <- lad_long |>
  filter(level == 'cell', quantity == 'power', family %in% c('1', '2', '3', '4'),
         analysis %in% c('G6', 'G10'),
         family != '1' | sub(' .*', '', scenario) == 'Hybrid') |>
  select(family, cell_group, analysis, value) |>
  pivot_wider(names_from = analysis, values_from = value) |>
  mutate(d = G10 - G6)
stopifnot(!anyNA(hy_cells$d))
hy_d_mean <- mean(hy_cells$d)
hy_below <- sum(hy_cells$d < 0)
hy_min <- min(hy_cells$d)
hy_max <- max(hy_cells$d)
## Size of G10 against the oracle in family S (independent seeds), with
## the 90 percent interval of the equivalence test (M6).
eqS <- lad_An |> filter(step_id == 'G8>G10', family == 'S')
eqS_lo <- eqS$size_diff - qnorm(0.95) * eqS$size_se
eqS_hi <- eqS$size_diff + qnorm(0.95) * eqS$size_se
stopifnot(nrow(eqS) == 1, eqS_lo > -0.01, eqS_hi < 0.01)
## OL+BDC CR2 check at 5,000 null datasets per cell, this paper's
## parameter set (analysis/data/cr2-olbdc-check/; the -published/
## directory beside it holds the same check under the starting values).
olb <- read.csv(here::here('analysis/data/cr2-olbdc-check/summary.csv')) |>
  filter(setting == 'p02', design == 'OLBDC', variant == 'g6')
stopifnot(nrow(olb) == 2, setequal(olb$N, c(35, 70)))
## kappa with the root-mean-square SE against the mean-SE kappa used
## throughout, CR2 analyses, by family and role (37-ladder-tables-
## figures.R, ladder-kappa-rms.csv; referee report 2026-11-08, m13).
krms <- read.csv(file.path(lad_dir, 'ladder-kappa-rms.csv')) |>
  filter(analysis %in% c('G6', 'G8', 'G9', 'G10'))
krms_d <- krms$kappa_rms - krms$kappa
## The text states that G9 stays the lowest of the four in every family
## and role, and below one at the alternatives, under either summary.
krms_g9 <- krms |>
  group_by(family, role) |>
  summarise(g9_low = all(kappa[analysis == 'G9'] <= kappa) &&
              all(kappa_rms[analysis == 'G9'] <= kappa_rms),
            .groups = 'drop')
krms_g10_alt_rms <- krms$kappa_rms[krms$analysis == 'G10' & krms$role == 'alt']
stopifnot(all(krms_d > 0), max(krms_d) < 0.04, all(krms_g9$g9_low),
          all(krms$kappa_rms[krms$analysis == 'G9' & krms$role == 'alt'] < 1))


## ---- scen-setup (formerly a report.Rmd chunk) ----
null_lab <- c(hybrid_null = 'Hybrid', olbdc_null = 'OL+BDC',
              t1hifhi_null = '$t_{1/2} = 2$ wk',
              shapemis_null = 'Weibull $k = 2$ data')
alt_lab <- c(hybrid_alt = 'Hybrid, $c_{bm} = 0.45$',
             hybrid_alt_lo = 'Hybrid, $c_{bm} = 0.25$',
             olbdc_alt = 'OL+BDC, $c_{bm} = 0.45$',
             t1hifhi_alt = '$t_{1/2} = 2$ wk, $c_{bm} = 0.25$',
             shapemis_alt = 'Weibull $k = 2$ data, $c_{bm} = 0.25$')
e0_null <- s21$reject[s21$method == 'E0' & s21$cell %in% names(null_lab)]
n6 <- s27 |> filter(c_bm == 0)
e0_6 <- n6$reject[n6$method == 'E0']
g1_6 <- n6$reject[n6$method == 'E1']

## ---- iso-numbers (formerly a report.Rmd chunk) ----
iso_lab <- c(E0 = 'A0', E0n = 'E0n', E0nc = 'E0nc', E1nc = 'E1nc',
             E1 = 'A1')
iso_rep <- unique(iso$n[iso$method != 'E0nc'])
stopifnot(length(iso_rep) == 1,
          all(iso$nonconv[iso$method != 'E0nc'] == 0),
          setequal(as.character(unique(iso$method)), names(iso_lab)))
iso_nc_max <- max(iso_rep - iso$n)
iso_rw <- iso_all |>
  mutate(rej = ifelse(converged & !is.na(p_value), p_value < 0.05, NA)) |>
  select(cell_id, design, N, c_bm, t1half, rep, method, rej) |>
  pivot_wider(names_from = method, values_from = rej)
iso_pw <- iso_all |>
  select(cell_id, rep, method, p_value) |>
  pivot_wider(names_from = method, values_from = p_value)
## Software and df method: largest p-value difference between E0 and
## E0n over all replicates.
iso_dp_sw <- max(abs(iso_pw$E0 - iso_pw$E0n), na.rm = TRUE)
## The exclusion rule of E0 drops no patient in this grid.
iso_no_drop <- all(iso_all$n_pat == iso_all$N)
stopifnot(iso_no_drop)
pd <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  d <- as.numeric(y[ok]) - as.numeric(x[ok])
  c(mean(d), sd(d) / sqrt(sum(ok)))
}
iso_st <- iso_rw |>
  group_by(design, N, c_bm, t1half) |>
  summarise(car_d = pd(E0n, E0nc)[1], car_se = pd(E0n, E0nc)[2],
            de_d = pd(E0n, E1nc)[1], de_se = pd(E0n, E1nc)[2],
            car1_d = pd(E1nc, E1)[1], car1_se = pd(E1nc, E1)[2],
            .groups = 'drop')
iso_de_null <- iso_st |> filter(c_bm == 0, design != 'CO')
iso_de_alt <- iso_st |> filter(c_bm > 0, design != 'CO')
iso_de_co <- iso_st |> filter(design == 'CO')
stopifnot(all(iso_de_co$de_d == 0))
iso_alt <- iso_st |> filter(c_bm > 0)
iso_gain <- iso_alt |> filter(car_d > 0)
iso_t1 <- iso_alt |> filter(t1half == 1)
## Pooled null rejection rate by design and analysis (over N).
iso_null <- iso_rw |>
  filter(c_bm == 0) |>
  pivot_longer(all_of(names(iso_lab)), names_to = 'method',
               values_to = 'rej') |>
  filter(!is.na(rej)) |>
  group_by(design, method) |>
  summarise(p = mean(rej), n = n(), .groups = 'drop') |>
  mutate(mcse = sqrt(p * (1 - p) / n))
inull <- function(d, m, w = 'p') {
  x <- iso_null |> filter(design == d, method == m)
  stopifnot(nrow(x) == 1)
  x[[w]]
}
iso_k <- function(m, ...) iso |> filter(method == m, c_bm == 0, ...) |>
  pull(kappa)
iso_cell <- function(m, d, nn, cb, th, w) {
  x <- iso |> filter(method == m, design == d, N == nn,
                     abs(c_bm - cb) < 1e-9, abs(t1half - th) < 1e-9)
  stopifnot(nrow(x) == 1)
  x[[w]]
}
iso_ratio <- iso |>
  filter(method %in% c('E0n', 'E0nc')) |>
  select(design, N, c_bm, t1half, method, mean_est, mean_se, emp_sd) |>
  pivot_wider(names_from = method,
              values_from = c(mean_est, mean_se, emp_sd)) |>
  mutate(se_ratio = mean_se_E0nc / mean_se_E0n,
         sd_ratio = emp_sd_E0nc / emp_sd_E0n,
         shrink = 1 - mean_est_E0nc / mean_est_E0n)
iso_shr <- function(d, th) {
  x <- iso_ratio |> filter(c_bm > 0, t1half == th, design %in% d)
  100 * range(x$shrink)
}
pct <- function(x) sprintf('%.0f', x)

## ---- stage4-grid (formerly a report.Rmd chunk) ----
sB <- sgf |>
  filter(block == 'B', sub %in% c('B1-shape', 'B2-halflife', 'B4-rho'),
         method %in% focus) |>
  mutate(level = case_when(
           sub == 'B1-shape' ~ sprintf('Weibull $k = %g$ data', weibull_shape),
           sub == 'B2-halflife' ~ sprintf('$t_{1/2} = %g$ wk', t1half),
           sub == 'B4-rho' ~ sprintf('$\\rho = %g$', rho)),
         Method = to_a(method)) |>
  arrange(match(sub, c('B1-shape', 'B2-halflife', 'B4-rho')),
          weibull_shape, t1half, rho, N)
alt_level <- c('B1-shape' = 0.45, 'B2-halflife' = 0.25, 'B4-rho' = 0.30)
grid_tab <- function(d) d |>
  select(Scenario = level, N, Method, reject) |>
  pivot_wider(names_from = Method, values_from = reject) |>
  select(Scenario, N, all_of(focus_cols)) |>
  mutate(across(-c(Scenario, N), fmt3))
## G10 against G1, G6, G8 and G9 over every alternative cell of the
## three scenario grids.
alt_gfB <- sgf |>
  filter(block == 'B', sub %in% c('B1-shape', 'B2-halflife', 'B4-rho'),
         c_bm > 0,
         method %in% c('E1', 'E1cr2', 'E2cr2', 'E7cr2', 'E7null_cr2')) |>
  select(cell_id, sub, N, t1half, weibull_shape, rho, method, reject) |>
  pivot_wider(names_from = method, values_from = reject) |>
  mutate(d_g1 = E7null_cr2 - E1, d_g6 = E7null_cr2 - E1cr2,
         d_g8 = E7null_cr2 - E2cr2, d_g9 = E7cr2 - E7null_cr2)
stopifnot(!anyNA(alt_gfB$E7null_cr2))

## ---- stage3-null computations (null cells of job 27) ----
## Null rejection by method at the six null cells, for the text.
st3n <- function(m) s27$reject[s27$c_bm == 0 & s27$method %in% m]
st3_g6max <- s27 |> filter(c_bm == 0, method == 'E1cr2') |>
  slice_max(reject, n = 1)
st3_g9g8 <- s27 |> filter(c_bm == 0, method %in% c('E2cr2', 'E7cr2')) |>
  select(design, N, method, reject) |>
  pivot_wider(names_from = method, values_from = reject) |>
  with(max(abs(E7cr2 - E2cr2)))

## ---- stage3-selection computations (A4 selected half-life) ----
## Half-life selected by A4 (stored E7null_cr2) per replicate, from the
## backfill refits of the job-27 grid and the supplementary study,
## Hybrid design, guarded replicates only.
sel_meta <- bind_rows(
  sgf |> distinct(cell_id, t1half),
  cell_meta27 |> filter(design == 'Hybrid') |> distinct(cell_id, t1half))
stopifnot(!anyDuplicated(sel_meta$cell_id))
m9 <- bf_all |>
  filter(method == 'E7null_cr2', source %in% c('gf', 't1'),
         unname(design) == 'Hybrid', converged, !is.na(h)) |>
  semi_join(bf_ok, by = c('source', 'cell_id', 'rep')) |>
  inner_join(sel_meta, by = 'cell_id') |>
  group_by(true_t1half = t1half, N = unname(N)) |>
  summarise(n = n(), p_0.25 = mean(h == 0.25), p_0.5 = mean(h == 0.5),
            p_1 = mean(h == 1), p_2 = mean(h == 2), .groups = 'drop') |>
  arrange(true_t1half, N)
stopifnot(all(abs(m9$p_0.25 + m9$p_0.5 + m9$p_1 + m9$p_2 - 1) < 1e-9))
m9_05 <- m9 |> filter(true_t1half == 0.5)

## ---- sample-size study (formerly the ss-setup chunk) ----
## 46-run-sample-size.R (six simrun jobs, 2026-11-11) and
## 47-sample-size-curves.R; 48-reconcile-sample-size-paper01.R compares
## the published Hybrid run with the companion paper's runs.
ss_dir <- here::here('analysis/data/carryover-sensitivity')
ss_n80 <- read.csv(file.path(ss_dir, 'sample-size-n80.csv'))
ss_fit <- read.csv(file.path(ss_dir, 'sample-size-fit.csv'))
ss_sep <- read.csv(file.path(ss_dir, 'sample-size-sepcheck.csv'))
ss_rec <- read.csv(file.path(ss_dir, 'sample-size-reconcile.csv'))
ss_reps <- min(vapply(Sys.glob(file.path(ss_dir, 'sample-size-*-*',
                                         'meta.rds')),
                      function(f) readRDS(f)$reps, numeric(1)))
ssN <- function(p, d, h, cb, a = 'A4', w = 'N80') {
  x <- ss_n80[ss_n80$params == p & ss_n80$design == d &
                ss_n80$t1half == h & abs(ss_n80$c_bm - cb) < 1e-9 &
                ss_n80$analysis == a, ]
  stopifnot(nrow(x) == 1)
  x[[w]]
}
## A2 needs a constant multiple of A4's N within a (set, design, t1/2),
## because the curve model scales with c_bm sqrt(N).
ss_ratio <- ss_n80 |>
  filter(abs(c_bm - 0.3) < 1e-9) |>
  select(params, design, t1half, analysis, N80) |>
  pivot_wider(names_from = analysis, values_from = N80) |>
  mutate(r = A2 / A4)
ssr <- function(d, h) {
  x <- ss_ratio |> filter(design == d, t1half == h)
  lo <- sprintf('%.0f', 100 * (min(x$r) - 1))
  hi <- sprintf('%.0f', 100 * (max(x$r) - 1))
  if (lo == hi) lo else paste(lo, 'to', hi)
}
## The text says A2 needs more N than A4 at a 1-week half-life in
## Hybrid and OL+BDC, and that the two are within a few percent
## elsewhere.
stopifnot(all(ss_ratio$r[ss_ratio$t1half == 1 & ss_ratio$design != 'CO'] > 1.1),
          all(abs(ss_ratio$r[!(ss_ratio$t1half == 1 & ss_ratio$design != 'CO')] - 1) < 0.06))
## The companion paper's A4 runs by serial correlation (Hybrid,
## starting values; 19-g9-planning-and-fig4b.R block N at rho 0.5, 0.7
## and 0.9). c_bm = 0.25 is the value common to all three.
rho_root <- here::here('analysis/data/quick-sim/carryover-closed-form',
                       'published-planning')
rho_n80 <- bind_rows(lapply(c(0.5, 0.7, 0.9), function(r) {
  f <- file.path(rho_root, if (r == 0.7) 'block-N' else
                   file.path(sprintf('rho%.1f', r), 'block-N'), 'n80-g10.csv')
  read.csv(f) |> mutate(rho = r)
})) |>
  filter(abs(c_bm - 0.25) < 1e-9, !extrapolated)
rho_r <- rho_n80 |>
  select(rho, t_half, N80_G) |>
  pivot_wider(names_from = rho, values_from = N80_G, names_prefix = 'r') |>
  mutate(r09 = r0.9 / r0.7, r05 = r0.5 / r0.7)
stopifnot(nrow(rho_r) == 3, !anyNA(rho_r))
ss_p_ratio <- ss_n80 |>
  filter(analysis == 'A4') |>
  select(params, design, t1half, c_bm, N80) |>
  pivot_wider(names_from = params, values_from = N80) |>
  mutate(r = published / simcore2026)
