## The strawman: a fresh run of Hendrickson et al. (2020) exactly as
## published, using the 58b32a9 ('orig') code itself. It should show both
## problems the adjustments are meant to fix:
##   1. correlation matrices that are not positive definite, silently
##      repaired by corpcor::make.positive.definite(tol = 1e-3);
##   2. power that collapses under small carryover (t1/2 = 0.1, 0.2).
##
## What is run: their generateSimulatedResults() -> generateData() (with
## the repair switched on, as in their driver) -> censordata() ->
## lme_analysis(), on their Figure 4 grid (4 designs x N {35, 70} x
## c.bm {0, .3, .6} x t1/2 {0, .1, .2} x no censoring + 4 censoring
## patterns), with the parameter sets saved in their results_core.rda and
## the analysis options of their vignette. Only three things are added
## around their code: (a) each of the 72 parameter sets runs as its own
## call, in parallel, with a fixed seed; (b) tic()/toc() progress printing
## is silenced; (c) make.positive.definite() is wrapped so every repair
## is logged with how far it moved the correlations.
##
## Requires 00-fetch-58b32a9.sh. censordata.R and plottingfunctions.R are
## read from the cache, or from the directory in H58_EXTRA. Run from the
## repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/10-strawman-58b32a9.R \
##     [--reps R] [--cores K]
suppressMessages({
  library(data.table); library(lme4); library(lmerTest); library(corpcor)
  library(MASS); library(ggplot2); library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
reps <- arg('--reps', 100L)
cores <- arg('--cores', 8L)

src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
extra <- Sys.getenv('H58_EXTRA', file.path(src, 'R'))
find_src <- function(f) {
  p <- c(file.path(src, 'R', f), file.path(extra, f))
  p <- p[file.exists(p)]
  if (!length(p)) stop('missing 58b32a9 source ', f)
  p[1]
}
out_dir <- 'analysis/data/quick-sim/hendrickson-problems/strawman'
fig_dir <- 'analysis/figures/quick-sim/hendrickson-problems'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

env <- new.env()
for (f in c('utilities.R', 'buildtrialdesign.R', 'generateData.R',
            'censordata.R', 'lme_analysis.R', 'generateSimulatedResults.R'))
  sys.source(find_src(f), envir = env)
env$tic <- function(...) invisible(NULL)
env$toc <- function(...) invisible(NULL)
env$.repairs <- list()
env$make.positive.definite <- function(m, tol) {
  r <- corpcor::make.positive.definite(m, tol = tol)
  R0 <- cov2cor(m); R1 <- cov2cor(r)
  br <- grep('\\.br$', colnames(m))
  env$.repairs[[length(env$.repairs) + 1]] <- data.table(
    min_eig = min(eigen(R0, TRUE, TRUE)$values),
    max_cor_change = max(abs(R1 - R0)),
    bm_br_spec = max(R0['bm', br]), bm_br_used = max(R1['bm', br]))
  r
}

pub <- new.env(); load(file.path(src, 'data', 'results_core.rda'), envir = pub)
ps <- pub$results_core$parameterselections
## Their vignette's analysis options; the saved copy lacks full_model_out,
## which lme_analysis() requires.
ap <- expand.grid(useDE = FALSE, t_random_slope = FALSE, full_model_out = FALSE)

## The designs saved in results_core come from an earlier code state
## (column 'timeptname'; the 58b32a9 lme_analysis() needs 'timeptnames'),
## so rebuild them with the 58b32a9 buildtrialdesign() and the calls in
## its vignette, then check the schedules match the saved ones.
btd <- env$buildtrialdesign
tds <- list(
  btd(name_longform = 'open label', name_shortform = 'OL',
      timepoints = env$cumulative(rep(2.5, 8)), timeptnames = paste0('OL', 1:8),
      expectancies = rep(1, 8), ondrug = list(pathA = rep(1, 8))),
  btd(name_longform = 'open label+blinded discontinuation', name_shortform = 'OL+BDC',
      timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
      timeptnames = c('OL1', 'OL2', 'OL3', 'OL4', 'BD1', 'BD2', 'BD3', 'BD4'),
      expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
      ondrug = list(pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
                    pathB = c(1, 1, 1, 1, 1, 0, 0, 0))),
  btd(name_longform = 'traditional crossover', name_shortform = 'CO',
      timepoints = env$cumulative(rep(2.5, 8)),
      timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
      expectancies = rep(.5, 8),
      ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
                    pathB = c(0, 0, 0, 0, 1, 1, 1, 1))),
  btd(name_longform = 'primary N-of-1 design', name_shortform = 'N-of-1',
      timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
      timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
      expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
      ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0), pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
                    pathC = c(1, 1, 1, 0, 0, 0, 1, 0), pathD = c(1, 1, 1, 0, 0, 0, 0, 1))))
num_cols <- c('t_wk', 'e', 'tod', 'tsd', 'tpb')
for (k in seq_along(tds)) {
  stopifnot(identical(tds[[k]]$metadata$name_shortform,
                      ps$trialdesigns[[k]]$metadata$name_shortform),
            length(tds[[k]]$trialpaths) == length(ps$trialdesigns[[k]]$trialpaths))
  for (g in seq_along(tds[[k]]$trialpaths))
    stopifnot(isTRUE(all.equal(
      as.data.frame(tds[[k]]$trialpaths[[g]])[num_cols],
      as.data.frame(ps$trialdesigns[[k]]$trialpaths[[g]])[num_cols],
      check.attributes = FALSE)))
}
ps$trialdesigns <- tds
grid <- CJ(td = seq_along(ps$trialdesigns), mp = seq_len(nrow(ps$modelparams)))
grid[, id := .I]

run_set <- function(k) {
  g <- grid[k]
  env$.repairs <- list()
  set.seed(20200 + g$id)
  t0 <- proc.time()[['elapsed']]
  res <- env$generateSimulatedResults(
    trialdesigns = ps$trialdesigns[g$td], respparamsets = ps$respparamsets[1],
    blparamsets = ps$blparamsets[1], censorparams = ps$censorparams,
    modelparams = ps$modelparams[g$mp, , drop = FALSE],
    ## Their loop reads saveunit2start even when not saving.
    simparam = list(Nreps = reps, progressiveSave = FALSE, saveunit2start = 1),
    analysisparams = ap, rawdataout = FALSE)
  r <- res$results
  r[, `:=`(trialdesign = g$td, modelparamset = g$mp,
           irep = (g$id - 1L) * reps + (irep - 1L) %% reps + 1L)]
  r[, warning := NULL]
  rp <- rbindlist(env$.repairs)
  list(results = r,
       repairs = data.table(id = g$id, trialdesign = g$td, modelparamset = g$mp,
                            n_repairs = nrow(rp),
                            min_eig = if (nrow(rp)) min(rp$min_eig) else NA_real_,
                            max_cor_change = if (nrow(rp)) max(rp$max_cor_change) else 0,
                            bm_br_spec = if (nrow(rp)) max(rp$bm_br_spec) else NA_real_,
                            bm_br_used = if (nrow(rp)) max(rp$bm_br_used) else NA_real_,
                            seconds = proc.time()[['elapsed']] - t0))
}

cat(sprintf('strawman: %d parameter sets x %d reps, %d cores\n',
            nrow(grid), reps, cores))
t_all <- proc.time()[['elapsed']]
out <- mclapply(grid$id, function(k) tryCatch(run_set(k), error = function(e)
  list(error = conditionMessage(e), id = k)), mc.cores = cores,
  mc.preschedule = FALSE)
err <- Filter(function(o) !is.null(o$error), out)
if (length(err)) {
  for (e in err) cat('set', e$id, 'failed:', e$error, '\n')
  stop(length(err), ' parameter sets failed')
}
results <- rbindlist(lapply(out, `[[`, 'results'))
repairs <- rbindlist(lapply(out, `[[`, 'repairs'))
cat(sprintf('done in %.1f min\n', (proc.time()[['elapsed']] - t_all) / 60))

ps_out <- ps
ps_out$analysisparams <- ap
ps_out$simparam <- list(Nreps = reps, seed_base = 20200)
strawman <- list(results = results, parameterselections = ps_out)
saveRDS(strawman, file.path(out_dir, sprintf('strawman-reps%d.rds', reps)))

## ---------------- problem 1: silent repairs ----------------
dn <- sapply(ps$trialdesigns, function(td) td$metadata$name_shortform)
mp <- as.data.table(ps$modelparams)[, mp := .I]
repairs <- merge(repairs, mp[, .(mp, N, c.bm, carryover_t1half)],
                 by.x = 'modelparamset', by.y = 'mp')
repairs[, design := dn[trialdesign]]
fwrite(repairs[order(trialdesign, modelparamset)],
       file.path(out_dir, sprintf('strawman-repairs-reps%d.csv', reps)))
cat('\nParameter sets whose matrices were silently repaired:\n')
print(repairs[, .(sets = .N, repaired = sum(n_repairs > 0),
                  bm_br_spec = if (any(n_repairs > 0)) max(bm_br_spec, na.rm = TRUE) else NA_real_,
                  bm_br_used = if (any(n_repairs > 0)) min(bm_br_used, na.rm = TRUE) else NA_real_),
              by = .(design, c.bm)][order(design, c.bm)])

## ---------------- problem 2: power, against the published ----------------
pw <- function(d) d[, .(power = mean(p < .05), n = .N),
                    by = .(trialdesign, N, c.bm, carryover_t1half, censorparamset)]
cmp <- merge(pw(results), pw(as.data.table(pub$results_core$results)),
             by = c('trialdesign', 'N', 'c.bm', 'carryover_t1half', 'censorparamset'),
             suffixes = c('_new', '_pub'))
cmp[, design := dn[trialdesign]]
cmp[, mcse := sqrt(power_pub * (1 - power_pub) / n_new + power_pub * (1 - power_pub) / n_pub)]
cmp[, z := (power_new - power_pub) / pmax(mcse, 1e-9)]
fwrite(cmp, file.path(out_dir, sprintf('strawman-power-vs-published-reps%d.csv', reps)))
cat(sprintf('\nPower vs published (%d cells): max |diff| %.3f; |z| > 2 in %d cells (%.0f%%; ~5%% expected)\n',
            nrow(cmp), max(abs(cmp$power_new - cmp$power_pub)), sum(abs(cmp$z) > 2),
            100 * mean(abs(cmp$z) > 2)))
cat('\nNo censoring, c.bm = 0.6, N = 70: power by carryover half-life\n')
print(dcast(cmp[censorparamset == 0 & c.bm == 0.6 & N == 70],
            design ~ carryover_t1half, value.var = 'power_new'))

## ---------------- Figure 4 layout, from their plotting code ----------------
penv <- new.env(); sys.source(find_src('plottingfunctions.R'), envir = penv)
prov <- sprintf('Strawman rerun of 58b32a9: %d replicates per cell (published: 1000)', reps)
pA <- penv$PlotModelingResults(strawman, 'power',
  c('trialdesign', 'N', 'c.bm', 'censorparamset'),
  data.table(blparamset = 1, respparamset = 1, carryover_t1half = 0),
  c('c.cfct', 'modelparamset')) +
  ggtitle('A. Effect of trial design, censoring and biomarker effect on power',
          subtitle = prov) + theme(plot.title = element_text(hjust = 0))
pB <- penv$PlotModelingResults(strawman, 'power',
  c('trialdesign', 'N', 'carryover_t1half', 'censorparamset'),
  data.table(blparamset = 1, respparamset = 1, c.bm = .6),
  c('c.cfct', 'modelparamset')) +
  ggtitle('B. Effect of nonzero carryover term on power', subtitle = prov) +
  theme(plot.title = element_text(hjust = 0))
ggsave(file.path(fig_dir, sprintf('strawman-fig4a-reps%d.png', reps)), pA,
       width = 9, height = 4.5, dpi = 200)
ggsave(file.path(fig_dir, sprintf('strawman-fig4b-reps%d.png', reps)), pB,
       width = 9, height = 4.5, dpi = 200)
cat('Wrote', file.path(fig_dir, sprintf('strawman-fig4[ab]-reps%d.png', reps)), '\n')
