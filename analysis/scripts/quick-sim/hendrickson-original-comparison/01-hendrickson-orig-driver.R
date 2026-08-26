suppressPackageStartupMessages({
  library(devtools); load_all('/Users/zenn/Library/CloudStorage/Dropbox/prj/res/36-pmsimstats-ng/pmsimstats-ng', quiet = TRUE)
  library(data.table); library(corpcor); library(MASS); library(lme4); library(lmerTest); library(parallel)
})

orig_env <- new.env(parent = globalenv())
sys.source(file.path(Sys.getenv('SCRATCH'), 'orig_generateData.R'), envir = orig_env)
sys.source(file.path(Sys.getenv('SCRATCH'), 'orig_lme_analysis.R'), envir = orig_env)
environment(orig_env$generateData) <- orig_env
environment(orig_env$lme_analysis) <- orig_env

data(extracted_bp); data(extracted_rp)
bp <- extracted_bp; rp <- extracted_rp

td_CO <- buildtrialdesign(name_shortform='CO', timepoints=cumulative(rep(2.5,8)),
  timeptnames=c(paste0('COa',1:4),paste0('COb',1:4)), expectancies=rep(.5,8),
  ondrug=list(pathA=c(1,1,1,1,0,0,0,0), pathB=c(0,0,0,0,1,1,1,1)))
td_Hybrid <- buildtrialdesign(name_shortform='Hybrid', timepoints=c(4,8,9,10,11,12,16,20),
  timeptnames=c('OL1','OL2','BD1','BD2','BD3','BD4','COd','COp'), expectancies=c(1,1,.5,.5,.5,.5,.5,.5),
  ondrug=list(pathA=c(1,1,1,1,0,0,1,0), pathB=c(1,1,1,1,0,0,0,1), pathC=c(1,1,1,0,0,0,1,0), pathD=c(1,1,1,0,0,0,0,1)))
td_OLBDC <- buildtrialdesign(name_shortform='OL+BDC', timepoints=c(4,8,12,16,17,18,19,20),
  timeptnames=c('OL1','OL2','OL3','OL4','BD1','BD2','BD3','BD4'), expectancies=c(1,1,1,1,.5,.5,.5,.5),
  ondrug=list(pathA=c(1,1,1,1,1,1,0,0), pathB=c(1,1,1,1,1,0,0,0)))

td_lookup <- list(CO=td_CO, Hybrid=td_Hybrid, `OL+BDC`=td_OLBDC)

allocate_across_paths <- function(n_total, n_paths) {
  base <- n_total %/% n_paths
  rep(base, n_paths) + c(rep(1L, n_total %% n_paths), rep(0L, n_paths - n_total %% n_paths))
}

N_TOTAL <- 70L
cells <- expand.grid(design = c('CO','Hybrid','OL+BDC'), t1half = c(0, 0.5, 1.0),
                      c.bm = c(0, 0.45), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
n_cells <- nrow(cells)
n_reps <- as.integer(Sys.getenv('NREPS', '300'))
cat('Cells:', n_cells, ' reps/cell:', n_reps, '\n')

run_one_rep <- function(design, c_bm, t1half, seed) {
  set.seed(seed)
  out_na <- function(reason) data.table(beta=NA_real_, betaSE=NA_real_, p=NA_real_, issingular=NA, warning=reason)
  tryCatch({
    td <- td_lookup[[design]]
    paths <- td$trialpaths
    n_per_path <- allocate_across_paths(N_TOTAL, length(paths))
    mp <- data.table(N = NA_integer_, c.bm = c_bm, carryover_t1half = t1half,
                      c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = 0.1)
    dat_list <- vector('list', length(paths))
    for (g in seq_along(paths)) {
      mp$N <- n_per_path[[g]]
      di <- orig_env$generateData(mp, rp, bp, paths[[g]], FALSE, TRUE)  # default scalefactor=2
      di[, path := g]
      dat_list[[g]] <- di
    }
    dat <- rbindlist(dat_list, fill = TRUE)
    # Faithful to Hendrickson's Produce_Publication_Results vignette:
    # analysisparams <- expand.grid(useDE=FALSE, t_random_slope=FALSE, full_model_out=FALSE)
    # -- no carryover_t1half in op, so lme_analysis's own fallback sets it to 0
    # (fixed/mis-specified Dbc treatment, not tracking the true DGP carryover).
    op <- list(useDE = FALSE, t_random_slope = FALSE, full_model_out = FALSE)
    res <- orig_env$lme_analysis(td$trialpaths, dat, op)
    if (!('beta' %in% names(res))) return(out_na('no beta col'))
    res
  }, error = function(e) out_na(conditionMessage(e)))
}

all_rows <- list()
for (i in seq_len(n_cells)) {
  ci <- cells[i, ]
  cat(sprintf('Cell %d/%d: design=%s t1half=%.2f c.bm=%.2f\n', i, n_cells, ci$design, ci$t1half, ci$c.bm))
  seeds <- 100000L*i + seq_len(n_reps)
  chunk_res <- mclapply(seeds, function(s) run_one_rep(ci$design, ci$c.bm, ci$t1half, s),
                         mc.cores = max(1, detectCores() - 1), mc.preschedule = TRUE)
  cr <- rbindlist(chunk_res, fill = TRUE)
  cr[, `:=`(design = ci$design, t1half = ci$t1half, c.bm = ci$c.bm, rep_idx = seq_len(.N))]
  all_rows[[i]] <- cr
}
replicates <- rbindlist(all_rows, fill = TRUE)
saveRDS(replicates, file.path(Sys.getenv('SCRATCH'), 'hendrickson_orig_replicates.rds'))

summary_dt <- replicates[, .(
  n_reps = .N, n_converged = sum(!is.na(p)),
  conv_rate = mean(!is.na(p)),
  power = mean(p < 0.05, na.rm = TRUE),
  mcse_power = sqrt(mean(p < 0.05, na.rm=TRUE)*(1-mean(p < 0.05, na.rm=TRUE))/sum(!is.na(p))),
  mean_beta = mean(beta, na.rm = TRUE), sd_beta = sd(beta, na.rm = TRUE)
), by = .(design, t1half, c.bm)]
setorder(summary_dt, design, c.bm, t1half)
print(summary_dt)
fwrite(summary_dt, file.path(Sys.getenv('SCRATCH'), 'hendrickson_orig_summary.txt'), sep = '\t')
cat('DONE\n')
