## analysis/scripts/carryover-sensitivity/39-check-e0-isolation.R
##
## Checks for 38-run-e0-isolation.R:
##   1. the formula and model call each method fits, shown for one
##      replicate of the Hybrid N = 70 cell;
##   2. E0 and E1 from the chain reproduce job 27's stored rejection
##      rates (tier1-g-vs-f-e0/summary.rds, 500 replicates, seeds
##      24102600 + r) in the same cells within Monte Carlo error. The
##      two runs use different seeds, so the comparison is between
##      independent estimates: a z statistic with the two MCSEs
##      combined, flagged when |z| > 3;
##   3. the chain's paired internal consistency (every method has one
##      row per replicate).
## Usage: Rscript analysis/scripts/carryover-sensitivity/39-check-e0-isolation.R \
##          [--reps 200] [--cores K]
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(parallel)
})
E0_ISOLATION_LIB_ONLY <- TRUE
source('analysis/scripts/carryover-sensitivity/38-run-e0-isolation.R')

args <- commandArgs(trailingOnly = TRUE)
argv <- function(k, d) if (k %in% args) args[which(args == k) + 1] else d
CHECK_REPS  <- as.integer(argv('--reps', '200'))
CHECK_CORES <- as.integer(argv('--cores', Sys.getenv('SIM_CORES', '4')))

cat('--- 1. Fitted models, one Hybrid N = 70 replicate ---\n')
hyb_id <- 'panelA-designHybrid-N70-cbm0.45-t1.00'
hyb <- as.list(grid[match(hyb_id, grid$cell_id), ])
set.seed(label_seed(sprintf('%s-rep-%d', hyb_id, 1)))
dl <- draw_trial(hyb)
s <- e0_setup(dl)
cat('E0   : lmerTest::lmer(', deparse(update(s$form, . ~ . + (1 | ptID))),
    ')  [patients: ', n_distinct(s$dd$ptID), ' of ', n_distinct(dl$ptID),
    ']\n', sep = '')
cat('E0n  : nlme::lme(', deparse(s$form),
    ', random = ~ 1 | ptID)  [same data as E0]\n', sep = '')
cat('E0nc : nlme::lme(', deparse(s$form),
    ', random = ~ 1 | ptID, correlation = corCAR1(~ t | ptID))',
    '  [same data as E0]\n', sep = '')
cat('E1nc : nlme::lme(', deparse(E1_FORM),
    ', random = ~ 1 | ptID)  [all ', n_distinct(dl$ptID), ' patients]\n',
    sep = '')
cat('E1   : nlme::lme(', deparse(E1_FORM),
    ', random = ~ 1 | ptID, correlation = corCAR1(~ t | ptID))  ',
    '(fit_spec_s7)\n', sep = '')
print(as.data.frame(fit_chain(dl)), digits = 4)

cat('\n--- 2. Reproduction of job 27 (E0, E1), ', CHECK_REPS,
    ' replicates ---\n', sep = '')
stored <- readRDS('analysis/data/carryover-sensitivity/tier1-g-vs-f-e0/summary.rds') |>
  filter(method %in% c('E0', 'E1')) |>
  select(cell_id, method, reject_27 = reject, n_27 = n)
check_ids <- c('panelA-designHybrid-N70-cbm0.45-t1.00',
               'panelA-designHybrid-N35-cbm0.00-t1.00',
               'panelA-designCO-N70-cbm0.45-t1.00')
jobs <- expand_grid(cell_id = check_ids, rep = seq_len(CHECK_REPS))
res <- bind_rows(mclapply(seq_len(nrow(jobs)), function(i) {
  cid <- jobs$cell_id[i]
  one_trial(as.list(grid[match(cid, grid$cell_id), ]), cid, jobs$rep[i])
}, mc.cores = CHECK_CORES, mc.preschedule = FALSE))

cmp <- res |>
  filter(method %in% c('E0', 'E1')) |>
  group_by(cell_id, method) |>
  summarise(n = sum(!is.na(p_value)), reject = mean(p_value < 0.05, na.rm = TRUE),
            .groups = 'drop') |>
  left_join(stored, by = c('cell_id', 'method')) |>
  mutate(mcse = sqrt(reject * (1 - reject) / n),
         mcse_27 = sqrt(reject_27 * (1 - reject_27) / n_27),
         z = (reject - reject_27) / sqrt(mcse^2 + mcse_27^2),
         ok = abs(z) <= 3)
print(as.data.frame(cmp |> mutate(across(where(is.double), ~ round(.x, 3)))),
      row.names = FALSE)
stopifnot(all(!is.na(cmp$reject_27)))
cat(if (all(cmp$ok)) 'E0 and E1 reproduce job 27 within 3 combined MCSE.\n'
    else 'WARNING: at least one cell differs from job 27 by more than 3 MCSE.\n')

cat('\n--- 3. Structure of the chain ---\n')
per_rep <- res |> count(cell_id, rep)
stopifnot(all(per_rep$n == length(METHODS)))
cat('Every (cell, rep) has exactly', length(METHODS), 'method rows.\n')
cat('Rejection rates in the checked cells, all methods:\n')
print(as.data.frame(res |> group_by(cell_id, method) |>
  summarise(reject = round(mean(p_value < 0.05, na.rm = TRUE), 3),
            nonconv = round(mean(!converged), 3), .groups = 'drop') |>
  pivot_wider(names_from = method, values_from = c(reject, nonconv))),
  row.names = FALSE)
