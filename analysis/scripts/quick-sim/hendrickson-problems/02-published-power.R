## Problem 2 evidence. Power from Hendrickson's own saved replicate-level
## results at 58b32a9 (results_core: 4 designs x 18 model cells x 5
## censoring conditions x 1,000 replicates), without dropout.
##
## Requires 00-fetch-58b32a9.sh. Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/02-published-power.R
## Whitepaper: docs/36-hendrickson-pd-and-carryover.md
suppressMessages(library(data.table))
src <- 'analysis/data/quick-sim/hendrickson-58b32a9'
out_dir <- 'analysis/data/quick-sim/hendrickson-problems'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

e <- new.env(); load(file.path(src, 'data/results_core.rda'), envir = e)
r <- as.data.table(e$results_core$results)
dnames <- sapply(e$results_core$parameterselections$trialdesigns,
                 function(td) td$metadata$name_shortform)
r[, design := dnames[trialdesign]]
pw <- r[censorparamset == 0,
        .(power = mean(p < 0.05, na.rm = TRUE), reps = .N,
          mean_beta = mean(beta, na.rm = TRUE),
          mean_se = mean(betaSE, na.rm = TRUE)),
        by = .(design, N, c.bm, carryover_t1half)]
fwrite(pw, file.path(out_dir, 'published-power.csv'))
cat(sprintf('replicates per cell: %s\n\n', paste(unique(pw$reps), collapse = ' ')))
cat('Published power, no dropout (columns: carryover half-life, weeks)\n')
print(dcast(pw, design + N + c.bm ~ carryover_t1half, value.var = 'power'),
      digits = 2)
cat('\nMean interaction estimate (columns: half-life)\n')
print(dcast(pw[c.bm > 0], design + N + c.bm ~ carryover_t1half,
            value.var = 'mean_beta'), digits = 3)
cat('\nResidual fraction one week after stopping, by half-life:\n')
print(setNames(0.5^(1 / c(0.1, 0.2, 0.5, 1)), c(0.1, 0.2, 0.5, 1)), digits = 3)
