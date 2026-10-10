## Checks for the bp fragment, section 1: shared biomarker effect,
## prognostic leak, and the expectancy-term cancellation.
dgp <- new.env()
keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build')
for (e in parse('analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R')) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
h <- dgp$designs$Hybrid
Ns <- c(18, 18, 17, 17)

## Implied covariance of (Y_0, ..., Y_8) for a path under an arm.
ycov <- function(k, arm) {
  pi <- dgp$path_info(h$t_wk, h$e, h$paths[[k]])
  b <- dgp$build(pi, arm, 0, 0)
  n <- b$n
  L <- matrix(0, n + 1, 2 + 3 * n); L[, 2] <- 1
  for (t in 1:n) L[t + 1, c(2 + t, 2 + n + t, 2 + 2 * n + t)] <- -1
  list(V = L %*% b$Sigma %*% t(L), pi = pi)
}

## Probability limit of the biomarker block of a GLS fit. Because the
## centered biomarker is independent of path and of every non-biomarker
## column, the normal equations separate and the biomarker coefficients
## are the GLS regression of the true slope vector s on Z, pooled over
## paths with allocation weights.
plim_bm <- function(arm, Zfun, sfun, working) {
  A <- 0; r <- 0
  for (k in 1:4) {
    yc <- ycov(k, arm)
    V <- switch(working,
      indep = diag(nrow(yc$V)),
      true  = yc$V,
      ri    = {
        off <- mean(yc$V[upper.tri(yc$V)]); d <- mean(diag(yc$V))
        m <- matrix(off, nrow(yc$V), nrow(yc$V)); diag(m) <- d; m
      })
    Z <- Zfun(yc$pi); s <- sfun(yc$pi)
    Vi <- solve(V)
    A <- A + Ns[k] * t(Z) %*% Vi %*% Z
    r <- r + Ns[k] * t(Z) %*% Vi %*% s
  }
  drop(solve(A, r))
}

dbc <- function(pi, ha = 0.5) {
  x <- ifelse(pi$on, 1, ifelse(pi$tsd > 0, 0.5^(pi$tsd / ha), 0))
  c(0, x)
}
Dbin <- function(pi) c(0, as.numeric(pi$on))
eta <- function(pi) c(0, pi$e)
bl0 <- function(pi) c(1, rep(0, length(pi$on)))
one <- function(pi) rep(1, length(pi$on) + 1)

Zs <- list(
  S_bin  = function(pi) cbind(one(pi), Dbin(pi)),
  F1_bin = function(pi) cbind(one(pi), Dbin(pi), eta(pi)),
  ES_bin = function(pi) cbind(one(pi), Dbin(pi), bl0(pi)),
  S_dbc  = function(pi) cbind(one(pi), dbc(pi)),
  F1_dbc = function(pi) cbind(one(pi), dbc(pi), eta(pi)),
  ES_dbc = function(pi) cbind(one(pi), dbc(pi), bl0(pi)))

## Prognostic pathways, per unit of prognostic slope.
paths_s <- list(
  TV_const = function(pi) c(0, rep(1, length(pi$on))),
  PB_eta   = function(pi) c(0, pi$e))

res <- expand.grid(model = names(Zs), pathway = names(paths_s),
                   arm = c('A', 'C'), working = c('indep', 'ri', 'true'),
                   stringsAsFactors = FALSE)
res$interaction <- mapply(function(m, p, a, w)
  plim_bm(a, Zs[[m]], paths_s[[p]], w)[2],
  res$model, res$pathway, res$arm, res$working)
res$interaction <- round(res$interaction, 4)
print(reshape(res, idvar = c('pathway', 'arm', 'working'),
              timevar = 'model', direction = 'wide'), row.names = FALSE)

## Closed-form checks under independence, binary D (four visit groups).
on_bl <- sapply(1:4, function(k) sum(h$paths[[k]][3:8] == 1))
off_bl <- sapply(1:4, function(k) sum(h$paths[[k]][3:8] == 0))
w0 <- sum(Ns); w1 <- 2 * sum(Ns); w2 <- sum(Ns * on_bl); w3 <- sum(Ns * off_bl)
cat(sprintf('\nweights: baseline %d, OL %d, blinded on %d, blinded off %d\n',
            w0, w1, w2, w3))
cat(sprintf('S leak (1-f) = w0/(w0+w3) = %.4f\n', w0 / (w0 + w3)))
K <- 1 / w0 + 1 / w1 + 1 / w2 + 1 / w3
cat(sprintf('F1 leak = (1/w3 - 1/w2)/K = %.4f\n', (1 / w3 - 1 / w2) / K))
cat('Per-path 1/(n_off+1):', 1 / (sapply(h$paths, function(p) sum(p == 0)) + 1), '\n')
