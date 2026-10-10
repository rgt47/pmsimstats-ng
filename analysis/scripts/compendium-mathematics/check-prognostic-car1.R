## Leak of S and F1 under the RI + CAR1 working covariance to which REML
## converges (approximated by one large fit), separable AR(1) arm, rho 0.7.
suppressMessages({library(MASS); library(nlme)})
dgp <- new.env()
keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build')
for (e in parse('analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R')) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
dgp$RHO <- 0.7
h <- dgp$designs$Hybrid
Ns <- c(18, 18, 17, 17)
set.seed(20261008)
arm <- 'C'
long <- do.call(rbind, lapply(1:4, function(k) {
  pi <- dgp$path_info(h$t_wk, h$e, h$paths[[k]])
  b <- dgp$build(pi, arm, 0, 0)
  n <- b$n; m <- 1500
  X <- mvrnorm(m, b$mu, b$Sigma)
  Y <- cbind(X[, 2], X[, 2] - (X[, 3:(2 + n)] + X[, (3 + n):(2 + 2 * n)] +
                                 X[, (3 + 2 * n):(2 + 3 * n)]))
  data.frame(ptID = rep(paste(k, 1:m), each = n + 1),
             t = rep(c(0, pi$w), m), Y = as.vector(t(Y)),
             bm = rep(X[, 1], each = n + 1), path = k,
             Dbc = rep(c(0, ifelse(pi$on, 1, ifelse(pi$tsd > 0,
                                       0.5^(pi$tsd / 0.5), 0))), m),
             De = rep(c(0, pi$e), m))
}))
long$bmc <- long$bm - mean(long$bm)
fit <- lme(Y ~ bmc + t + Dbc + bmc:Dbc + De + bmc:De, random = ~ 1 | ptID,
           correlation = corCAR1(form = ~ t | ptID), data = long,
           method = 'ML')
su2 <- as.numeric(VarCorr(fit)[1, 1]); se2 <- fit$sigma^2
phi <- coef(fit$modelStruct$corStruct, unconstrained = FALSE)
cat(sprintf('RI var %.1f, resid var %.1f, CAR1 phi %.3f\n', su2, se2, phi))

tt <- c(0, cumsum(h$t_wk))
Vw <- su2 + se2 * phi^abs(outer(tt, tt, '-'))
dbc <- function(pi) c(0, ifelse(pi$on, 1, ifelse(pi$tsd > 0, 0.5^(pi$tsd / 0.5), 0)))
leak <- function(Zfun, sfun, V = Vw) {
  A <- 0; r <- 0
  for (k in 1:4) {
    pi <- dgp$path_info(h$t_wk, h$e, h$paths[[k]])
    Z <- Zfun(pi); s <- sfun(pi); Vi <- solve(V)
    A <- A + Ns[k] * t(Z) %*% Vi %*% Z; r <- r + Ns[k] * t(Z) %*% Vi %*% s
  }
  drop(solve(A, r))[2]
}
one <- function(pi) rep(1, 9)
ZS <- function(pi) cbind(one(pi), dbc(pi))
ZF <- function(pi) cbind(one(pi), dbc(pi), c(0, pi$e))
ZE <- function(pi) cbind(one(pi), dbc(pi), c(1, rep(0, 8)))
sTV <- function(pi) c(0, rep(1, 8))
sPB <- function(pi) c(0, pi$e)
out <- data.frame(model = c('S', 'F1', 'ES'),
  TV_const = c(leak(ZS, sTV), leak(ZF, sTV), leak(ZE, sTV)),
  PB_eta = c(leak(ZS, sPB), leak(ZF, sPB), leak(ZE, sPB)))
print(out, digits = 3)
