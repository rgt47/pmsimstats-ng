## analysis/scripts/quick-sim/biomarker-scaling/02-free-baseline-check.R
##
## docs/78, Section 3. Usage (from the repository root):
##   Rscript analysis/scripts/quick-sim/biomarker-scaling/\
##     02-free-baseline-check.R
## Synthetic check of the free-baseline model (ES in paper 01's stress
## test): Sx ~ b + t + X + b:X + b:bl0, with no bl0 main effect, where
## b is the biomarker raw, centered at the sample baseline mean, or
## standardized by the sample. Also counts lmer scale warnings.
suppressMessages({library(nlme); library(clubSandwich); library(lme4)})
one <- function(seed, beta, off) {
  set.seed(seed)
  N <- 70; nv <- 8
  d <- expand.grid(v = 0:nv, ptID = 1:N)
  bm <- rnorm(N, 124.3, 15.4)
  d$bm <- bm[d$ptID]
  d$t <- d$v * 2.5
  d$X <- ifelse(d$v == 0, 0,
    as.numeric((d$v %in% 3:5) == (d$ptID %% 2 == 0)))
  d$bl0 <- as.numeric(d$v == 0)
  u <- rnorm(N, 0, 8)
  d$Sx <- 80 - 0.5 * d$t - 3 * d$X + off * d$bl0 +
    beta * (d$bm - 124.3) * d$X + u[d$ptID] + rnorm(nrow(d), 0, 6)
  mb <- mean(bm); sb <- sd(bm)
  d$braw <- d$bm; d$bc <- d$bm - mb; d$bz <- (d$bm - mb) / sb
  fit <- function(v) {
    f <- as.formula(sprintf('Sx ~ %s + t + X + %s:X + %s:bl0', v, v, v))
    g <- lme(f, random = ~1 | ptID, data = d,
             correlation = corCAR1(form = ~ t | ptID))
    k <- coef_test(g, vcov = 'CR2', coefs = sprintf('%s:X', v))
    s <- if (v == 'bz') sb else 1
    c(est = k$beta / s, p = k$p_Satt)
  }
  w <- 0L
  withCallingHandlers(
    lmer(Sx ~ braw * X + t + (1 | ptID), data = d),
    warning = function(e) { w <<- w + 1L; invokeRestart('muffleWarning') },
    message = function(e) { w <<- w + 1L; invokeRestart('muffleMessage') })
  c(raw = fit('braw'), cen = fit('bc'), z = fit('bz'), lmer_warn = w)
}
R <- 300
for (cfg in list(c(0, 0), c(0, 5), c(-0.13, 5))) {
  x <- t(sapply(1:R, one, beta = cfg[1], off = cfg[2]))
  cat(sprintf('\nbeta %.2f, baseline offset %g, %d reps\n',
              cfg[1], cfg[2], R))
  cat(sprintf('  reject raw %.3f  centered %.3f  standardized %.3f\n',
      mean(x[, 'raw.p'] < .05), mean(x[, 'cen.p'] < .05),
      mean(x[, 'z.p'] < .05)))
  cat(sprintf('  mean est (per raw unit) raw %.4f  centered %.4f  z %.4f\n',
      mean(x[, 'raw.est']), mean(x[, 'cen.est']), mean(x[, 'z.est'])))
  cat(sprintf('  max |p diff| centered vs z %.2e; raw vs centered %.3f\n',
      max(abs(x[, 'cen.p'] - x[, 'z.p'])),
      max(abs(x[, 'raw.p'] - x[, 'cen.p']))))
  cat(sprintf('  lmer raw-scale warnings in %d of %d fits\n',
      sum(x[, 'lmer_warn'] > 0), R))
}
