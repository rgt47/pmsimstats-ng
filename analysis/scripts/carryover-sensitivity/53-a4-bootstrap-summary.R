## analysis/scripts/carryover-sensitivity/53-a4-bootstrap-summary.R
##
## Combines the per-replicate output of 52-a4-bootstrap-interval.R and
## scores five nominal 95 percent intervals for A4 at the reference
## cell against theta = -c_bm * sigma_BR = -3.6:
##
##   cr2_normal   estimate +/- 1.96 CR2 SE (the stored interval uses
##                Satterthwaite df, which are not stored; Table 7 of
##                the paper gives its coverage bounds)
##   percentile   2.5 and 97.5 percent quantiles of the bootstrap
##                estimates
##   basic        2 * estimate minus the same quantiles
##   boot_se      estimate +/- 1.96 SD of the bootstrap estimates
##   boot_t       estimate - q * CR2 SE, with q the 97.5 and 2.5
##                percent quantiles of (bootstrap estimate - estimate)
##                / bootstrap CR2 SE
##
## Only replicates whose A4 fit passed the pairing guard are used.
## Writes <data>/a4-bootstrap/all.rds and
## <data>/referee-reanalysis/a4-bootstrap-coverage.csv.
##
## Usage (from the repository root):
##   Rscript analysis/scripts/carryover-sensitivity/53-a4-bootstrap-summary.R
## With PMSIM_PARAMS=published it reads a4-bootstrap-published/ and
## writes referee-reanalysis-published/a4-bootstrap-coverage.csv. The
## true value is the same under both sets (sigma_BR = 8).
suppressMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
})
DATA_DIR <- 'analysis/data/carryover-sensitivity'
PUB <- Sys.getenv('PMSIM_PARAMS', 'simcore2026') == 'published'
BOOT_DIR <- file.path(DATA_DIR, if (PUB) 'a4-bootstrap-published' else 'a4-bootstrap')
RR_DIR <- file.path(DATA_DIR, if (PUB) 'referee-reanalysis-published' else 'referee-reanalysis')
THETA <- -0.45 * 8

files <- list.files(file.path(BOOT_DIR, 'trials'), pattern = '^rep-.*\\.rds$',
                    full.names = TRUE)
reps <- lapply(files, readRDS)
obs <- map_dfr(reps, ~ tibble(rep = .x$rep, guard = .x$guard, est = .x$est,
                              se = .x$se, p = .x$p, h = .x$h))
boot <- map_dfr(reps, function(x)
  if (is.null(x$boot)) NULL else mutate(x$boot, rep = x$rep, .before = 1))
saveRDS(list(obs = obs, boot = boot), file.path(BOOT_DIR, 'all.rds'))
cat(sprintf('%d replicates, guard held in %d; %d resamples, %d failed\n',
            nrow(obs), sum(obs$guard), nrow(boot), sum(is.na(boot$est))))

ok <- obs |> filter(guard)
iv <- boot |>
  filter(!is.na(est), !is.na(se), se > 0) |>
  inner_join(ok |> select(rep, est0 = est, se0 = se), by = 'rep') |>
  group_by(rep) |>
  summarise(
    n_boot = n(),
    est0 = first(est0), se0 = first(se0),
    q_lo = quantile(est, 0.025, names = FALSE),
    q_hi = quantile(est, 0.975, names = FALSE),
    sd_b = sd(est),
    t_lo = quantile((est - est0) / se, 0.025, names = FALSE),
    t_hi = quantile((est - est0) / se, 0.975, names = FALSE),
    .groups = 'drop')

lims <- bind_rows(
  iv |> transmute(rep, interval = 'cr2_normal',
                  lo = est0 - 1.96 * se0, hi = est0 + 1.96 * se0),
  iv |> transmute(rep, interval = 'percentile', lo = q_lo, hi = q_hi),
  iv |> transmute(rep, interval = 'basic',
                  lo = 2 * est0 - q_hi, hi = 2 * est0 - q_lo),
  iv |> transmute(rep, interval = 'boot_se',
                  lo = est0 - 1.96 * sd_b, hi = est0 + 1.96 * sd_b),
  iv |> transmute(rep, interval = 'boot_t',
                  lo = est0 - t_hi * se0, hi = est0 - t_lo * se0))

cov <- lims |>
  group_by(interval) |>
  summarise(n = n(), coverage = mean(lo <= THETA & THETA <= hi),
            mcse = sqrt(coverage * (1 - coverage) / n),
            miss_low = mean(hi < THETA), miss_high = mean(lo > THETA),
            mean_width = mean(hi - lo), .groups = 'drop') |>
  mutate(min_boot = min(iv$n_boot), median_boot = median(iv$n_boot),
         sd_est = sd(ok$est), mean_se = mean(ok$se),
         mean_sd_boot = mean(iv$sd_b))
print(as.data.frame(cov), digits = 3)
write.csv(cov, file.path(RR_DIR, 'a4-bootstrap-coverage.csv'),
          row.names = FALSE)
