## Reproduce stored paper 09 null replicates (OLBDC, condition 'none',
## analysis ta0) draw for draw: same seed formula and the same RNG
## sequence as run_unit() in 20-paired-dropout.R (generateData, then
## the runif() dropout draws, per path). Fits only ta0 on the
## no-dropout data, then compares with the stored rows.
suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table)
  library(parallel)
  loadNamespace('clubSandwich')
})
env09 <- new.env(parent = globalenv())
exprs <- parse('analysis/scripts/informative-dropout-by-design/20-paired-dropout.R',
               keep.source = FALSE)
stop_at <- which(vapply(exprs, function(e) {
  is.call(e) && identical(e[[1]], as.name('<-')) &&
    identical(e[[2]], as.name('units'))
}, logical(1)))[1]
for (e in exprs[seq_len(stop_at - 1L)]) eval(e, envir = env09)

n_reps <- as.integer(Sys.getenv('CHK_REPS', '40'))
n_cores <- as.integer(Sys.getenv('CHK_CORES', '2'))

one <- function(r) local({
  rep_i <- r
  design <- 'OLBDC'
  d_idx <- match(design, names(designs))
  set.seed(P09_SEED + 100000L * d_idx + rep_i + 50000L)
  paths <- designs[[design]]$trialpaths
  n_pp <- allocate_across_paths(N_TOTAL, length(paths))
  mp <- make_mp(0)
  templates <- lapply(paths, path_long_template)
  y_list <- bl_list <- bm_list <- fd_list <- vector('list', length(paths))
  for (g in seq_along(paths)) {
    mp$N <- n_pp[g]
    di <- generateData(mp, rp, bp, paths[[g]], empirical = FALSE,
                       makePositiveDefinite = FALSE,
                       cached_sigma = sigma_cache[[design]][[g]],
                       dgp_architecture = 'mean_moderation',
                       stopOnNonPD = TRUE)
    tn <- paths[[g]]$timeptnames
    y_list[[g]] <- as.matrix(di[, ..tn])
    bl_list[[g]] <- di$BL
    bm_list[[g]] <- di$bm
    invisible(runif(length(y_list[[g]])))
    fd_list[[g]] <- rep(NA_integer_, n_pp[g])
  }
  long <- build_long(y_list, bl_list, bm_list, fd_list, templates)
  data.table(rep = rep_i, as.data.table(fit_analysis(long, 0, TRUE)))
}, envir = list2env(list(r = r), parent = env09))

new <- rbindlist(mclapply(seq_len(n_reps), one, mc.cores = n_cores))
old <- as.data.table(readRDS(
  'analysis/data/informative-dropout-by-design/paired-replicates.rds'))
old <- old[design == 'OLBDC' & regime == 'null' & condition == 'none' &
           analysis == 'ta0' & rep <= n_reps]
m <- merge(new, old[, .(rep, est_old = est, se_cr2_old = se_cr2,
                        p_cr2_old = p_cr2)], by = 'rep')
print(m[, .(n = .N, max_d_est = max(abs(est - est_old)),
            max_d_se = max(abs(se_cr2 - se_cr2_old)),
            max_d_p = max(abs(p_cr2 - p_cr2_old)))])
print(head(m[, .(rep, est, est_old, p_cr2, p_cr2_old)], 8))

all_old <- as.data.table(readRDS(
  'analysis/data/informative-dropout-by-design/paired-replicates.rds'))[
  design == 'OLBDC' & regime == 'null' & condition == 'none' &
    analysis == 'ta0']
z <- all_old$est / all_old$se_cr2
cat(sprintf('stored: n %d, rej %.3f, sd(t_cr2) %.3f, kurtosis(est) %.2f\n',
            nrow(all_old), mean(all_old$p_cr2 < 0.05), sd(z),
            mean(scale(all_old$est)^4)))
cat(sprintf('stored rejections in reps 1-250: %d, 251-500: %d\n',
            sum(all_old[rep <= 250]$p_cr2 < 0.05),
            sum(all_old[rep > 250]$p_cr2 < 0.05)))
