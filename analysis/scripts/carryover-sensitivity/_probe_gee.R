suppressPackageStartupMessages({
  library(dplyr); library(geepack); library(geesmv)
})
repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

set.seed(1)
design_set <- build_design_set('Hybrid')
mp <- list(N = 70, c.bm = 0, carryover_t1half = 1.0,
           carryover_form = 'exponential', weibull_shape = 1,
           dgp_architecture = 'mvn',
           c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
           c.cf1t = 0.2, c.cfct = 0.1)
dat <- generate_data_multi_path(mp, default_resp_param(),
                                default_baseline_param(), design_set)
dl <- prepare_long_data(dat, design_set, carryover_t1half = 1.0,
                        carryover_form = 'exponential',
                        weibull_shape = 1)
dl <- dl[order(dl$ptID, dl$t), ]
dl$id <- as.integer(factor(dl$ptID))
cat('n clusters:', length(unique(dl$id)), '\n')

fit <- geepack::geeglm(Sx ~ bm + t + Db + bm:Db, family = gaussian,
                       data = dl, id = id, corstr = 'ar1')
cat('geeglm coef:', paste(names(coef(fit)), collapse = ', '), '\n')

md <- tryCatch(
  geesmv::GEE.var.md(Sx ~ bm + t + Db + bm:Db, id = 'id',
    family = gaussian, data = as.data.frame(dl), corstr = 'AR-1'),
  error = function(e) conditionMessage(e))
cat('\n--- GEE.var.md result ---\n')
str(md)
if (is.list(md)) {
  cat('cov.beta:', paste(sprintf('%.4f', md$cov.beta),
                         collapse = ', '), '\n')
  cat('naive SEs:',
      paste(sprintf('%.4f', summary(fit)$coefficients[,'Std.err']),
            collapse = ', '), '\n')
}
