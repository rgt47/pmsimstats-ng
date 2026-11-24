## analysis/scripts/quick-sim/biomarker-scaling/01-invariance-check.R
##
## docs/78, Section 2. Fits the interaction model to 20 synthetic
## data sets (N = 70, eight visits, biomarker N(124.3, 15.4^2)) with
## the biomarker raw and standardized, under lmer (the published
## analysis) and nlme corCAR1 with CR2 SEs (the paper 02 core), and
## prints the largest differences in t, df and p. Synthetic data, not
## the papers' DGP.
##
## Usage (from the repository root):
##   Rscript analysis/scripts/quick-sim/biomarker-scaling/\
##     01-invariance-check.R
suppressMessages({library(lmerTest); library(nlme); library(clubSandwich)})
set.seed(1)
res <- list()
for (r in 1:20) {
  N <- 70; nv <- 8
  d <- expand.grid(v = 1:nv, ptID = 1:N)
  bm <- rnorm(N, 124.3, 15.4)
  d$bm <- bm[d$ptID]
  d$t <- d$v * 2.5
  d$Db <- as.numeric((d$v %in% c(3:5)) == (d$ptID %% 2 == 0))
  u <- rnorm(N, 0, 8)
  d$Sx <- 80 - 0.5 * d$t - 3 * d$Db - 0.13 * (d$bm - 124.3) * d$Db +
    u[d$ptID] + rnorm(nrow(d), 0, 6)
  d$bz <- (d$bm - 124.3) / 15.4
  ## A0 as in 58b32a9: lmer, random intercept, Satterthwaite
  f1 <- lmer(Sx ~ bm + Db + t + bm:Db + (1 | ptID), data = d)
  f2 <- lmer(Sx ~ bz + Db + t + bz:Db + (1 | ptID), data = d)
  c1 <- coef(summary(f1))['bm:Db', ]; c2 <- coef(summary(f2))['bz:Db', ]
  ## nlme corCAR1 + CR2, as in the paper 02 core
  g1 <- lme(Sx ~ bm + t + Db + bm:Db, random = ~1 | ptID,
            correlation = corCAR1(form = ~ t | ptID), data = d)
  g2 <- lme(Sx ~ bz + t + Db + bz:Db, random = ~1 | ptID,
            correlation = corCAR1(form = ~ t | ptID), data = d)
  k1 <- coef_test(g1, vcov = 'CR2', coefs = 'bm:Db')
  k2 <- coef_test(g2, vcov = 'CR2', coefs = 'bz:Db')
  res[[r]] <- data.frame(
    lmer_ratio = c2['Estimate'] / c1['Estimate'],
    lmer_dt = c1['t value'] - c2['t value'],
    lmer_dp = c1['Pr(>|t|)'] - c2['Pr(>|t|)'],
    lme_dt = k1$tstat - k2$tstat,
    cr2_ddf = k1$df_Satt - k2$df_Satt,
    cr2_dp = k1$p_Satt - k2$p_Satt)
}
x <- do.call(rbind, res)
cat('estimate ratio (expect 15.4):', range(x$lmer_ratio), '\n')
print(sapply(x[-1], function(v) max(abs(v))))
