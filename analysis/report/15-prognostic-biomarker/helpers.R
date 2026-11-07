## analysis/report/15-prognostic-biomarker/helpers.R
##
## Data loading and formatting helpers for report.Rmd. This file reads
## stored simulation summaries only. It runs no simulation.
##
## Two kinds of stored run are read (docs/74, decision D5):
##   published-set runs  analysis/data/carryover-sensitivity/published/
##                       Hendrickson et al.'s parameter values,
##                       stop on non-positive-definite matrices.
##   historical exhibit  analysis/data/carryover-sensitivity/<job>/
##                       the 'simcore2026' parameter set, run with
##                       silent repair before 2026-10-26.
## Every file is optional. A missing file sets has_* to FALSE, and the
## tables and inline values that need it print a pending note.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

if (!exists('repo_root')) {
  repo_root <- local({
    d <- normalizePath(getwd())
    while (!file.exists(file.path(d, 'DESCRIPTION')) &&
             d != dirname(d)) d <- dirname(d)
    d
  })
}
cs_dir  <- file.path(repo_root, 'analysis', 'data',
                     'carryover-sensitivity')
pub_dir <- file.path(cs_dir, 'published')

rd <- function(path) if (file.exists(path)) readRDS(path) else NULL
rd_csv <- function(path) {
  if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE) else NULL
}

src <- list(
  stress_p = file.path(pub_dir, 'g9-prognostic-stress'),
  gapC_p   = file.path(pub_dir, 'g9-gap-fill', 'C'),
  blkT_p   = file.path(pub_dir, 'referee-variants', 'block-T'),
  blkN_p   = file.path(pub_dir, 'referee-variants', 'block-N'),
  blkA_p   = file.path(pub_dir, 'referee-variants', 'block-A'),
  tier1_p  = file.path(pub_dir, 'tier1-g-vs-f-e0'),
  scen_p   = file.path(pub_dir, 'g-vs-f-comparison'),
  stress_h = file.path(cs_dir, 'g9-prognostic-stress'),
  re_h     = file.path(cs_dir, 'g9-random-effects'),
  rr_h     = file.path(cs_dir, 'referee-reanalysis')
)

## Published-set runs.
st_p     <- rd(file.path(src$stress_p, 'summary.rds'))
st_p_m   <- rd(file.path(src$stress_p, 'meta.rds'))
gfC      <- rd(file.path(src$gapC_p, 'summary.rds'))
gfC_all  <- rd(file.path(src$gapC_p, 'all.rds'))
gfC_m    <- rd(file.path(src$gapC_p, 'meta.rds'))
rvT      <- rd(file.path(src$blkT_p, 'summary.rds'))
rvTp     <- rd(file.path(src$blkT_p, 'paired.rds'))
rvT_m    <- rd(file.path(src$blkT_p, 'meta.rds'))
rvN      <- rd(file.path(src$blkN_p, 'summary.rds'))
rvNp     <- rd(file.path(src$blkN_p, 'paired.rds'))
rvN_m    <- rd(file.path(src$blkN_p, 'meta.rds'))
rvA      <- rd(file.path(src$blkA_p, 'summary.rds'))
rvA_m    <- rd(file.path(src$blkA_p, 'meta.rds'))
s27      <- rd(file.path(src$tier1_p, 'summary.rds'))
s27_all  <- rd(file.path(src$tier1_p, 'all.rds'))
s27_m    <- rd(file.path(src$tier1_p, 'meta.rds'))
s21      <- rd(file.path(src$scen_p, 'summary.rds'))
s21_m    <- rd(file.path(src$scen_p, 'meta.rds'))

## Historical exhibit (simcore2026, silent repair).
st_h     <- rd(file.path(src$stress_h, 'summary.rds'))
st_h_m   <- rd(file.path(src$stress_h, 'meta.rds'))
pdc      <- rd(file.path(src$stress_h, 'pd-check.rds'))
ceil     <- rd(file.path(src$stress_h, 'coupling-ceilings.rds'))
re_h     <- rd(file.path(src$re_h, 'summary.rds'))
re_h_m   <- rd(file.path(src$re_h, 'meta.rds'))
m12      <- rd_csv(file.path(src$rr_h, 'm12-random-effects.csv'))

has_st_p <- !is.null(st_p)
has_gfC  <- !is.null(gfC)
has_gfCa <- !is.null(gfC_all)
has_rvT  <- !is.null(rvT) && !is.null(rvTp)
has_rvN  <- !is.null(rvN) && !is.null(rvNp)
has_rvA  <- !is.null(rvA)
has_s27  <- !is.null(s27)
has_s27a <- !is.null(s27_all)
has_s21  <- !is.null(s21)
has_st_h <- !is.null(st_h)
has_pdc  <- !is.null(pdc)
has_ceil <- !is.null(ceil)
has_re_h <- !is.null(re_h)
has_m12  <- !is.null(m12)

## The parameter set a stored run records. Runs made before the
## parameter switch (2026-11-03) record none.
param_set <- function(meta) {
  ps <- meta$params$param_set
  if (is.null(ps)) 'not recorded' else ps
}

## A published-set run must say so in its metadata; a run that does
## not is reported as a provenance failure in the consistency checks.
pub_ok <- function(meta) identical(param_set(meta), 'published')

## The drug-response SD of the published set, for the true value of
## the interaction, -c_bm * sigma_BR per SD of the biomarker.
sigma_br <- local({
  m <- if (!is.null(gfC_m)) gfC_m else st_p_m
  rp <- m$params$resp_param
  if (is.null(rp)) NA_real_ else rp$sd[rp$cat == 'br']
})
sigma_pb <- local({
  m <- if (!is.null(gfC_m)) gfC_m else st_p_m
  rp <- m$params$resp_param
  if (is.null(rp)) NA_real_ else rp$sd[rp$cat == 'pb']
})

## The simcore2026 parameter set is not recorded in the historical
## runs' metadata. It is read from its definition in the simulation
## core, by evaluating that one function and nothing else.
simcore_rp <- local({
  f <- file.path(repo_root, 'analysis', 'scripts',
                 'carryover-sensitivity', 'simulation-core.R')
  out <- NULL
  if (file.exists(f)) {
    env <- new.env()
    for (e in parse(f)) {
      if (is.call(e) && identical(e[[1]], as.name('<-')) &&
            identical(e[[2]], as.name('simcore2026_resp_param'))) {
        eval(e, env)
      }
    }
    if (exists('simcore2026_resp_param', envir = env)) {
      out <- env$simcore2026_resp_param()
    }
  }
  out
})
sigma_pb_simcore <- if (is.null(simcore_rp)) NA_real_ else
  simcore_rp$sd[simcore_rp$cat == 'pb']

## Display labels. Stored codes stay E/F; the display layer maps them.
m_lab <- c(E0 = 'E0', E1 = 'G1', E2 = 'G3', E7 = 'G4', E1cr2 = 'G6',
           E2cr2 = 'G8', E7cr2 = 'G9', E7null_cr2 = 'G10', F1 = 'F1',
           F1A = 'F1A', F1A0 = 'F1A0', F1A0B = 'F1A0B')
show_m <- c('E2cr2', 'E7cr2', 'E7null_cr2', 'F1', 'F1A0', 'F1A0B')
show_c <- unname(m_lab[show_m])
path_lab <- c(progbl = '$c_{bm,BL}$', progtvconst = '$c_{bm,TV}$',
              progpb = '$c_{bm,PB}$',
              progtvgrow = '$b_{TV}$ (constant)')
path_ord <- c('progpb', 'progtvconst', 'progbl', 'progtvgrow')

## Formatting.
f2  <- function(x) ifelse(is.na(x), '--', sprintf('%.2f', x))
f3  <- function(x) ifelse(is.na(x), '--', sprintf('%.3f', x))
f4  <- function(x) ifelse(is.na(x), '--', sprintf('%.4f', x))
fs2 <- function(x) ifelse(is.na(x), '--', sprintf('%+.2f', x))
fs3 <- function(x) ifelse(is.na(x), '--', sprintf('%+.3f', x))
fs4 <- function(x) ifelse(is.na(x), '--', sprintf('%+.4f', x))
fn  <- function(x) format(x, big.mark = ',', scientific = FALSE,
                          trim = TRUE)
rng <- function(x, f = f3) {
  lo <- f(min(x))
  hi <- f(max(x))
  if (lo == hi) lo else paste0(lo, ' to ', hi)
}

## Two-MCSE band around a nominal 0.05 at n replicates.
band <- function(n) 0.05 + c(-2, 2) * sqrt(0.05 * 0.95 / n)
band_txt <- function(n) {
  b <- band(n)
  sprintf('$[%.3f, %.3f]$', b[1], b[2])
}
in_band <- function(x, n) x >= band(n)[1] & x <= band(n)[2]

## Inline values. Each is set in a chunk when its data exist; V()
## returns a pending marker otherwise.
V_env <- new.env()
setv <- function(k, x) assign(k, x, envir = V_env)
V <- function(k) {
  if (exists(k, envir = V_env, inherits = FALSE)) {
    get(k, envir = V_env)
  } else {
    '[pending]'
  }
}

## Qualitative statements in the prose are checked against the data.
## A failed check does not stop the render. It is listed in the
## appendix so that the text can be revised.
claim_env <- new.env()
claim_env$n <- 0L
claim_env$fail <- character(0)
claim <- function(ok, text) {
  claim_env$n <- claim_env$n + 1L
  if (!isTRUE(all(ok))) claim_env$fail <- c(claim_env$fail, text)
  invisible(ok)
}

## One row of a summary; the filter must identify exactly one row.
row1 <- function(d, ...) {
  x <- dplyr::filter(d, ...)
  stopifnot(nrow(x) == 1)
  x
}

## Paired difference in rejection between two analyses fitted to the
## same replicates of one cell (m2 minus m1), from a per-replicate
## file. Rejection is p < 0.05 on the stored p-value. F1A0B is not
## supported, because its critical value is a pooled bootstrap value.
paired_diff <- function(all, cell, m1, m2) {
  a <- all |>
    dplyr::filter(cell_id == cell, method %in% c(m1, m2),
                  !is.na(p_value)) |>
    dplyr::transmute(rep, method, r = as.numeric(p_value < 0.05)) |>
    tidyr::pivot_wider(names_from = method, values_from = r)
  a <- a[stats::complete.cases(a[, c(m1, m2)]), ]
  d <- a[[m2]] - a[[m1]]
  c(diff = mean(d), se = stats::sd(d) / sqrt(length(d)),
    n = length(d))
}

## A table or a pending note.
pending_note <- function(what) {
  knitr::asis_output(paste0('\\textit{Pending: ', what,
                            ' is not in the stored data.}\n'))
}
