## Summarize the stored null replicates of paper 09 and paper 02
## (OL+BDC and controls) before the new runs.
suppressPackageStartupMessages(library(data.table))
d <- as.data.table(readRDS(
  'analysis/data/informative-dropout-by-design/paired-replicates.rds'))
print(d[regime == 'null' & condition == 'none',
        .(n = .N, emp_sd = sd(est, na.rm = TRUE),
          se_mb = mean(se_mb, na.rm = TRUE),
          se_cr2 = mean(se_cr2, na.rm = TRUE),
          rej_mb = mean(p_mb < 0.05, na.rm = TRUE),
          rej_cr2 = mean(p_cr2 < 0.05, na.rm = TRUE),
          df_cr2_med = median(df_cr2, na.rm = TRUE),
          car1 = mean(car1)),
        by = .(design, analysis)])

t1 <- as.data.table(readRDS(
  'analysis/data/carryover-sensitivity/tier1-g-vs-f-e0/all.rds'))
print(t1[c_bm == 0 & method %in% c('E1', 'E1cr2', 'E2', 'E2cr2'),
         .(n = sum(!is.na(p_value)), emp_sd = sd(estimate, na.rm = TRUE),
           se = mean(std_error, na.rm = TRUE),
           rej = mean(p_value < 0.05, na.rm = TRUE)),
         by = .(design, N, method)][order(design, N, method)])
