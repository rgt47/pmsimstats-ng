## Copied 2026-10-08 from the review scratch script inv_f_check.R;
## run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-baseline-design-search.R
## Independent checks for the inventory part f (docs/39, 41, 42, 45,
## 48, 49, 51). Run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-baseline-design-search.R
## Everything below is recomputed from first principles (own construct
## builder, own power functions); stored CSVs are read only to compare.
suppressMessages({library(data.table); library(MASS); library(nlme)})
dd <- 'analysis/data/quick-sim/carryover-closed-form'
hdr <- function(s) cat('\n==================', s, '==================\n')

## ---------------- construct (published values) ----------------
gomp <- function(t, m, d, r) {
  y <- m * exp(-d * exp(-r * t)); off <- m * exp(-d)
  (y - off) * m / (m - off)
}
RHO <- 0.8; C1 <- 0.2; CX <- 0.1
sTV <- 10; sPB <- 10; sBR <- 8; sbm <- 15.36159; sBL <- 18.48267
BETA <- 0.25; ALPHA <- 0.05

## Path information from absolute visit times (weeks), expectancy and
## on/off pattern: tod and tsd accumulate within a run; tsd is zero
## before first exposure.
pinfo <- function(times, e, od) {
  n <- length(od); gaps <- diff(c(0, times))
  tod <- tsd <- numeric(n); exposed <- FALSE
  for (i in seq_len(n)) {
    if (od[i] == 1) {
      tod[i] <- gaps[i] + if (i > 1 && od[i - 1] == 1) tod[i - 1] else 0
      exposed <- TRUE
    } else if (exposed) {
      tsd[i] <- gaps[i] + if (i > 1 && od[i - 1] == 0) tsd[i - 1] else 0
    }
  }
  list(w = times, e = e, on = od == 1, tod = tod, tsd = tsd, n = n)
}
hprof <- function(p, th) {
  off <- if (th > 0) 0.5^(p$tsd / th) else rep(0, p$n)
  off[p$tsd <= 0] <- 0
  ifelse(p$on, 1, off)
}
## Response correlation, component-major (TV, PB, BR).
rcor <- function(p, structure, rho = RHO) {
  n <- p$n
  if (structure == 'CS') {
    A <- matrix(rho, n, n); diag(A) <- 1
    X <- matrix(CX, n, n); diag(X) <- C1
    return(rbind(cbind(A, X, X), cbind(X, A, X), cbind(X, X, A)))
  }
  nu <- as.numeric(sub('nug', '', structure))
  K <- (1 - C1) * diag(3) + C1
  kronecker(K, (1 - nu) * rho^abs(outer(p$w, p$w, '-')) + nu * diag(n))
}
rcov <- function(p, structure, rho = RHO) {
  s <- c(rep(sTV, p$n), sPB * p$e, rep(sBR, p$n))
  outer(s, s) * rcor(p, structure, rho)
}
CSum <- function(p, structure, rho = RHO) {
  L <- cbind(diag(p$n), diag(p$n), diag(p$n))
  L %*% rcov(p, structure, rho) %*% t(L)
}
wts <- function(p, use = rep(TRUE, p$n), extra_off = 0) {
  on <- p$on & use; off <- !p$on & use
  ifelse(on, 1 / sum(on), ifelse(off, -1 / (sum(off) + extra_off), 0))
}
alloc <- function(P, N = 70L) N %/% P + c(rep(1L, N %% P), rep(0L, P - N %% P))

## Decayed mean moderation moments of one path (E9 on the used visits).
mm_moments <- function(p, th, structure, use = rep(TRUE, p$n)) {
  a <- wts(p, use); h <- hprof(p, th)
  aCa <- drop(t(a) %*% CSum(p, structure) %*% a)
  sh <- BETA * sBR * sum(a * h)
  list(gamma = -BETA * sBR / sbm * sum(a * h), varD = aCa + sh^2, tau2 = aCa,
       hoff = mean(h[!p$on & use]))
}
## Stratified (sequence-specific intercept) power: pooled slope over P
## strata, S ~ sbm^2 chi^2_{N-P}, df N - P - 1; slope differences enter
## the residual to first order. Gauss-Laguerre-free: quantile quadrature
## on 4000 points.
strat_power <- function(m, Ns) {
  P <- length(Ns); Nt <- sum(Ns); w <- (Ns - 1) / sum(Ns - 1)
  g <- sapply(m, `[[`, 'gamma'); t2 <- sapply(m, `[[`, 'tau2')
  vD <- sapply(m, `[[`, 'varD')
  gb <- sum(w * g); res <- sum(w * (t2 + sbm^2 * (g - gb)^2))
  df <- Nt - P - 1; q <- qchisq((seq_len(4000) - 0.5) / 4000, Nt - P)
  ncp <- gb * sbm * sqrt(q) / sqrt(res); tc <- qt(1 - ALPHA / 2, df)
  list(power = mean(pt(-tc, df, ncp = ncp) + pt(tc, df, ncp = ncp, lower.tail = FALSE)),
       info = sum(Ns / Nt * (g * sbm)^2 / vD))
}
## Common-intercept pooled power (E9 as specified), with path means mu.
pooled_power <- function(g, t2, mu, Ns) {
  N <- sum(Ns); w <- Ns / N; gb <- sum(w * g); mb <- sum(w * mu)
  res <- sum(w * (t2 + (mu - mb)^2 + sbm^2 * (g - gb)^2))
  q <- qchisq((seq_len(4000) - 0.5) / 4000, N - 1); tc <- qt(1 - ALPHA / 2, N - 2)
  ncp <- gb * sbm * sqrt(q) / sqrt(res)
  list(power = mean(pt(-tc, N - 2, ncp = ncp) + pt(tc, N - 2, ncp = ncp, lower.tail = FALSE)),
       mean = gb, sd = sqrt(res / ((N - 3) * sbm^2)), res = res)
}

hyb_t <- c(4, 8, 9, 10, 11, 12, 16, 20)
hyb_e <- c(1, 1, rep(.5, 6))
hyb_p <- list(c(1,1,1,1,0,0,1,0), c(1,1,1,1,0,0,0,1), c(1,1,1,0,0,0,1,0), c(1,1,1,0,0,0,0,1))
olb_t <- c(4, 8, 12, 16, 17, 18, 19, 20); olb_e <- c(1, 1, 1, 1, rep(.5, 4))
olb_p <- list(c(1,1,1,1,1,1,0,0), c(1,1,1,1,1,0,0,0))
co_t <- 2.5 * (1:8)
pat <- function(s) as.numeric(strsplit(s, '')[[1]])
comp <- function(s) 1 - pat(s)

## ================= docs/39, Section 3: components of Var(Delta) =========
hdr('docs/39 S3: CS components of Var(Delta) and AR(1) time part')
cs_parts <- function(p) {
  a <- wts(p); e <- p$e
  m1 <- sum(a); m2 <- sum(a^2); me <- sum(a * e); mee <- sum(a^2 * e^2); m2e <- sum(a^2 * e)
  c(balance = (sTV^2 + sBR^2) * (RHO * m1^2 + (1 - RHO) * m2),
    pb_occ = sPB^2 * (1 - RHO) * mee,
    exp_conf = sPB^2 * RHO * me^2,
    cross = 2 * sTV * sBR * (CX * m1^2 + (C1 - CX) * m2) +
      2 * (sTV + sBR) * sPB * (CX * m1 * me + (C1 - CX) * m2e))
}
designs <- list(
  'Hybrid p1' = pinfo(hyb_t, hyb_e, hyb_p[[1]]), 'Hybrid p2' = pinfo(hyb_t, hyb_e, hyb_p[[2]]),
  'Hybrid p3' = pinfo(hyb_t, hyb_e, hyb_p[[3]]), 'Hybrid p4' = pinfo(hyb_t, hyb_e, hyb_p[[4]]),
  'Hybrid blinded p1' = pinfo(hyb_t, rep(.5, 8), hyb_p[[1]]),
  'Hybrid blinded p3' = pinfo(hyb_t, rep(.5, 8), hyb_p[[3]]),
  'OL+BDC p1' = pinfo(olb_t, olb_e, olb_p[[1]]), 'OL+BDC p2' = pinfo(olb_t, olb_e, olb_p[[2]]),
  'CO p1' = pinfo(co_t, rep(.5, 8), pat('11110000')),
  'CO OL-first p1' = pinfo(co_t, c(1,1,1,1,.5,.5,.5,.5), pat('11110000')),
  'ABAB 2.5wk' = pinfo(co_t, rep(.5, 8), pat('11001100')),
  'Alternating 2.5wk' = pinfo(co_t, rep(.5, 8), pat('10101010')),
  'ABAB at Hybrid times' = pinfo(hyb_t, rep(.5, 8), pat('11001100')))
tab3 <- rbindlist(lapply(names(designs), function(nm) {
  p <- designs[[nm]]; a <- wts(p); pr <- cs_parts(p)
  A <- RHO^abs(outer(p$w, p$w, '-'))
  data.table(design = nm, on = sum(p$on), off = sum(!p$on),
             balance = pr[1], pb_occ = pr[2], exp_conf = pr[3], cross = pr[4],
             sum_parts = sum(pr), qf_CS = drop(t(a) %*% CSum(p, 'CS') %*% a),
             aAa = drop(t(a) %*% A %*% a),
             var_AR1 = drop(t(a) %*% CSum(p, 'nug0') %*% a))
}))
print(tab3[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat('vKv blinded (TV,PB=5,BR):', sum(outer(c(10, 5, 8), c(10, 5, 8)) *
                                    ((1 - C1) * diag(3) + C1)), '\n')
ref <- fread(file.path(dd, 'design-contrast-structure.csv'))
cat('stored var_cs range', range(ref$var_cs), ' var_sep range', range(ref$var_sep), '\n')

## Power: Hybrid E9 published vs blinded; CO vs CO with an OL first
## period (stratified E9, decayed mean moderation, CS, t1/2 = 0).
hdr('docs/39 S3: power consequences (CS, t1/2 = 0)')
pw <- function(paths, Ns, th = 0, structure = 'CS', use = NULL)
  strat_power(lapply(paths, function(p) mm_moments(p, th, structure,
                                                  if (is.null(use)) rep(TRUE, p$n) else use)), Ns)$power
cat(sprintf('Hybrid stratified %.3f, blinded %.3f\n',
            pw(lapply(hyb_p, function(o) pinfo(hyb_t, hyb_e, o)), alloc(4)),
            pw(lapply(hyb_p, function(o) pinfo(hyb_t, rep(.5, 8), o)), alloc(4))))
cat(sprintf('CO stratified %.3f, OL-first %.3f\n',
            pw(list(pinfo(co_t, rep(.5, 8), pat('11110000')),
                    pinfo(co_t, rep(.5, 8), comp('11110000'))), c(35L, 35L)),
            pw(list(pinfo(co_t, c(1,1,1,1,.5,.5,.5,.5), pat('11110000')),
                    pinfo(co_t, c(1,1,1,1,.5,.5,.5,.5), comp('11110000'))), c(35L, 35L))))

## ================= docs/39 S4, S6: exhaustive search ==================
hdr('docs/39 S4/S6: exhaustive search (own implementation)')
pats <- as.matrix(expand.grid(rep(list(0:1), 8)))
pats <- pats[rowSums(pats) >= 2 & rowSums(pats) <= 6 & pats[, 1] == 1, ]
cat('designs:', nrow(pats), '\n')
STR <- c('CS', 'nug0.2', 'nug0.5'); TH <- c(0.5, 1, 2); DG <- c(0.5, 1, 1.5, 2, 2.5, 3, 4)
eval_seq <- function(pt, d, th, s, k = 0, use_ol = TRUE) {
  nb <- length(pt); n <- nb + k; times <- d * seq_len(n)
  e <- c(rep(1, k), rep(.5, nb)); use <- c(rep(use_ol, k), rep(TRUE, nb))
  ps <- list(pinfo(times, e, c(rep(1, k), pt)), pinfo(times, e, c(rep(1, k), 1 - pt)))
  strat_power(lapply(ps, mm_moments, th = th, structure = s, use = use), c(35L, 35L))
}
res <- rbindlist(lapply(seq_len(nrow(pats)), function(i) {
  pt <- pats[i, ]
  rbindlist(lapply(DG, function(d) rbindlist(lapply(STR, function(s) rbindlist(lapply(TH,
    function(th) { r <- eval_seq(pt, d, th, s)
      data.table(pattern = paste(pt, collapse = ''), sw = sum(diff(pt) != 0), n_on = sum(pt),
                 d = d, s = s, th = th, info = r$info, power = r$power) }))))))
}))
best <- res[, .SD[which.max(info)], by = .(s, th)]
print(best[, .(s, th, pattern, sw, d, power = round(power, 3))])
stb <- fread(file.path(dd, 'design-search/design-search-best.csv'))
cmpb <- merge(best, stb, by.x = c('s', 'th'), by.y = c('structure', 't_half'))
cat('best: same pattern', all(cmpb$pattern.x == sprintf('%08s', cmpb$pattern.y) |
                               cmpb$pattern.x == as.character(cmpb$pattern.y)),
    ' max |power diff|', signif(max(abs(cmpb$power.x - cmpb$power.y)), 3), '\n')
bysw <- res[, .(p = max(power)), by = .(s, th, sw)]
print(dcast(bysw[, .(s, th, sw, p = round(p, 3))], s + th ~ sw, value.var = 'p'))
## CS landscape bands at t1/2 = 0.5, by number of on-drug visits.
print(res[s == 'CS' & th == 0.5, .(p = max(power)), by = .(pattern, n_on)][
  , .(min = round(min(p), 3), max = round(max(p), 3)), by = n_on][order(n_on)])
rel <- merge(res, best[, .(s, th, bp = power)], by = c('s', 'th'))
rob <- rel[, .(min_ratio = min(power / bp), mean_ratio = mean(power / bp)), by = .(pattern, d)][
  order(-min_ratio)]
print(head(rob[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], 6))
cat('alternating at 4 weeks:', unlist(rob[pattern == '10101010' & d == 4, .(min_ratio, mean_ratio)]), '\n')

## ================= docs/39 S5: spacing ==================
hdr('docs/39 S5: spacing (common intercept, path means removed, as in 03)')
arr <- list(alternating = '10101010', ABAB = '11001100', CO = '11110000')
dgr <- seq(0.25, 6, by = 0.25)
sp <- rbindlist(lapply(names(arr), function(an) rbindlist(lapply(c('CS', 'nug0', 'nug0.2', 'nug0.5'),
  function(s) rbindlist(lapply(c(0.25, 0.5, 1, 2), function(th) rbindlist(lapply(dgr, function(d) {
    times <- d * (1:8)
    ps <- list(pinfo(times, rep(.5, 8), pat(arr[[an]])), pinfo(times, rep(.5, 8), comp(arr[[an]])))
    m <- lapply(ps, mm_moments, th = th, structure = s)
    g <- sapply(m, `[[`, 'gamma'); t2 <- sapply(m, `[[`, 'tau2'); vD <- sapply(m, `[[`, 'varD')
    data.table(arr = an, s = s, th = th, d = d, info = mean((g * sbm)^2 / vD),
               power = pooled_power(g, t2, c(0, 0), c(35L, 35L))$power)
  }))))))))
opt <- sp[, .SD[which.max(info)], by = .(arr, s, th)]
d95 <- sp[, .(d95 = min(d[power >= 0.95 * max(power)])), by = .(arr, s, th)]
print(dcast(merge(opt, d95)[, .(arr, s, th, v = sprintf('d*=%.2f d95=%.2f', d, d95))],
            arr + s ~ th, value.var = 'v'))

## ================= docs/39 S7: open-label lead-in ==================
hdr('docs/39 S7: open-label lead-in (best power over patterns and spacings)')
olr <- list()
for (k in 0:4) {
  nb <- 8 - k
  bp <- as.matrix(expand.grid(rep(list(0:1), nb)))
  bp <- bp[rowSums(bp) >= 1 & rowSums(bp) <= nb - 1 & bp[, 1] == 1, , drop = FALSE]
  for (u in c(TRUE, FALSE)) {
    if (k == 0 && !u) next
    olr[[length(olr) + 1]] <- rbindlist(lapply(seq_len(nrow(bp)), function(i)
      rbindlist(lapply(DG, function(d) rbindlist(lapply(STR, function(s) rbindlist(lapply(TH,
        function(th) { r <- eval_seq(bp[i, ], d, th, s, k, u)
          data.table(k = k, inc = u, s = s, th = th, info = r$info, power = r$power) }))))))))
  }
}
ol <- rbindlist(olr)
olb <- ol[, .SD[which.max(info)], by = .(k, inc, s, th)]
print(dcast(olb[, .(s, th, col = sprintf('k%d%s', k, ifelse(inc, 'in', 'out')),
                    p = round(power, 3))], s + th ~ col, value.var = 'p'))

## ================= docs/39 S8: blocks ==================
hdr('docs/39 S8: BD and CO blocks')
eval_des <- function(times, e, paths, th, s, ol_in) {
  ps <- lapply(paths, function(o) pinfo(times, e, o))
  use <- c(rep(ol_in, 2), rep(TRUE, length(times) - 2))
  m <- lapply(ps, mm_moments, th = th, structure = s, use = use)
  r <- strat_power(m, alloc(length(paths)))
  data.table(power = r$power, var_mean = mean(sapply(m, `[[`, 'tau2')),
             h_off = mean(sapply(m, `[[`, 'hoff')))
}
ol <- c(4, 8); blk <- c(9, 10, 11, 12)
blocks <- list(
  BD109 = list(c(1,1,1,1,0,0), c(1,1,1,0,0,0)), BD119 = list(c(1,1,1,1,1,0), c(1,1,1,0,0,0)),
  CO = list(c(1,1,1,1,0,0), c(1,1,0,0,1,1)), ALT = list(c(1,1,1,0,1,0), c(1,1,0,1,0,1)))
mk <- function(bdp, cop) { o <- list(); for (b in bdp) for (cc in cop) o[[length(o) + 1]] <- c(1, 1, b, cc); o }
bd <- list(c(1,1,0,0), c(1,0,0,0)); co2 <- list(c(1,0), c(0,1)); co4 <- list(c(1,1,0,0), c(0,0,1,1))
bd2 <- list(c(1,0), c(0,0))
hyb <- list(
  Published = list(c(ol, blk, 16, 20), mk(bd, co2), c(1, 1, rep(.5, 6))),
  EqualVisits = list(c(ol, 10, 12, 14, 16, 18, 20), mk(bd2, co4), c(1, 1, rep(.5, 6))),
  COweekly = list(c(ol, blk, 13, 14, 15, 16), mk(bd, co4), c(1, 1, rep(.5, 8))),
  CO2weekly = list(c(ol, blk, 14, 16, 18, 20), mk(bd, co4), c(1, 1, rep(.5, 8))))
br <- list()
for (s in c('CS', 'nug0.2')) for (th in c(0, .5, 1)) for (oi in c(TRUE, FALSE)) {
  for (nm in names(blocks)) br[[length(br) + 1]] <- cbind(set = 'block', name = nm, s = s, th = th,
    inc = oi, eval_des(c(ol, blk), c(1, 1, rep(.5, 4)), blocks[[nm]], th, s, oi))
  for (nm in names(hyb)) br[[length(br) + 1]] <- cbind(set = 'hybrid', name = nm, s = s, th = th,
    inc = oi, eval_des(hyb[[nm]][[1]], hyb[[nm]][[3]], hyb[[nm]][[2]], th, s, oi))
}
br <- rbindlist(br)
print(dcast(br[inc == TRUE, .(name, col = paste(s, th), p = round(power, 3))], name ~ col,
            value.var = 'p'))
print(br[inc == TRUE & s == 'CS' & th == 1, .(name, h_off = round(h_off, 3))])
print(br[inc == TRUE & s == 'CS' & th == 0, .(name, var_mean = round(var_mean, 2))])
print(br[s == 'nug0.2' & th == 0 & name %in% c('COweekly', 'EqualVisits'),
         .(name, inc, p = round(power, 3))])
print(br[s == 'CS' & th == 0.5 & name == 'ALT', .(inc, p = round(power, 3))])

## ================= docs/39 S10 item 6: the optimal contrast ==========
hdr('docs/39 S10.6: optimal contrast among sum-zero contrasts')
opt_info <- function(C, h) {
  Ci <- solve(C); one <- rep(1, length(h))
  drop(t(h) %*% Ci %*% h - (t(one) %*% Ci %*% h)^2 / (t(one) %*% Ci %*% one))
}
e9_info <- function(C, h, a) drop((sum(a * h))^2 / (t(a) %*% C %*% a))
set.seed(1)
for (case in list(list('Hybrid p1, CS, t=0', pinfo(hyb_t, hyb_e, hyb_p[[1]]), 'CS', 0),
                  list('Hybrid p1, nug0.2, t=0', pinfo(hyb_t, hyb_e, hyb_p[[1]]), 'nug0.2', 0),
                  list('Alternating 2.5wk, nug0.2, t=1', pinfo(co_t, rep(.5, 8), pat('10101010')), 'nug0.2', 1),
                  list('CO 2.5wk, nug0.2, t=1', pinfo(co_t, rep(.5, 8), pat('11110000')), 'nug0.2', 1),
                  list('CO 2.5wk, CS, t=1', pinfo(co_t, rep(.5, 8), pat('11110000')), 'CS', 1))) {
  p <- case[[2]]; C <- CSum(p, case[[3]]); h <- hprof(p, case[[4]])
  oi <- opt_info(C, h); ei <- e9_info(C, h, wts(p))
  ## GLS on the biomarker block (columns 1 and h) gives the same information.
  Z <- cbind(1, h); gls <- 1 / solve(t(Z) %*% solve(C) %*% Z)[2, 2]
  ## Random sum-zero contrasts never exceed the optimum.
  rnd <- max(replicate(20000, { a <- rnorm(p$n); a <- a - mean(a); e9_info(C, h, a) }))
  cat(sprintf('%-32s opt %.5f  GLS %.5f  E9 %.5f  E9/opt %.3f  max random %.5f\n',
              case[[1]], oi, gls, ei, ei / oi, rnd))
}

## ================= docs/41: construct covariance ==================
hdr('docs/41: construct covariance, propositions 3 and 4')
p1 <- pinfo(hyb_t, hyb_e, hyb_p[[1]])
CS1 <- CSum(p1, 'CS')
SY <- matrix(sBL^2, 9, 9); SY[2:9, 2:9] <- SY[2:9, 2:9] + CS1
cat(sprintf('Var(BL) %.1f; OL var %.1f cov %.1f; blinded var %.1f cov %.1f; OL-blinded %.1f\n',
            SY[1, 1], SY[2, 2], SY[2, 3], SY[4, 4], SY[4, 5], SY[2, 4]))
cat(sprintf('cov with baseline: %s\n', paste(round(SY[1, 2:9], 1), collapse = ' ')))
cat(sprintf('cor BL-OL %.3f BL-blinded %.3f; OL-OL %.3f blinded-blinded %.3f\n',
            SY[1, 2] / sqrt(SY[1, 1] * SY[2, 2]), SY[1, 4] / sqrt(SY[1, 1] * SY[4, 4]),
            SY[2, 3] / SY[2, 2], SY[4, 5] / SY[4, 4]))
cond <- SY[2:9, 2:9] - outer(SY[2:9, 1], SY[1, 2:9]) / SY[1, 1]
cat(sprintf('conditional: OL var %.1f cov %.1f; blinded var %.1f cov %.1f\n',
            cond[1, 1], cond[1, 2], cond[3, 3], cond[3, 4]))
cat('phi:', round(SY[2:9, 1] / SY[1, 1], 6), '\n')
cat(sprintf('Prop 3 blinded: unadj %.1f change %.1f conditional %.1f\n', SY[4, 4],
            SY[4, 4] + SY[1, 1] - 2 * SY[1, 4], SY[4, 4] - SY[1, 4]^2 / SY[1, 1]))
for (j in c(1, 3)) {
  p <- pinfo(hyb_t, hyb_e, hyb_p[[j]]); C <- CSum(p, 'CS')
  a9 <- wts(p); a9b <- wts(p, extra_off = 1); noff <- sum(!p$on)
  b <- ifelse(p$on, 0, 1 / noff)
  v9 <- drop(t(a9) %*% C %*% a9); v9b <- drop(t(a9b) %*% C %*% a9b)
  vSoff <- drop(t(b) %*% C %*% b); cv <- drop(t(a9) %*% C %*% b)
  cat(sprintf(paste('path %d: Var E9 %.2f E9b %.2f; identity %.2f;',
                    'Var(Sbar_off) %.2f Cov(D9,Sbar_off) %.2f; added %.2f offset %.2f (%.0f%%);',
                    'approx 200/(noff+1)^2 = %.2f, offset share %.0f%%\n'),
              j, v9, v9b, v9 + 2 * cv / (noff + 1) + vSoff / (noff + 1)^2,
              vSoff, cv, vSoff / (noff + 1)^2, -2 * cv / (noff + 1),
              100 * (-2 * cv / (noff + 1)) / (vSoff / (noff + 1)^2),
              200 / (noff + 1)^2, 100 * (200 / (noff + 1)^2 - (v9b - v9)) / (200 / (noff + 1)^2)))
}

## E9 and E9b pooled power and SD (published construct means; common
## intercept; covariance binding under step rule, plain mean moderation).
hdr('docs/41 S6: E9 and E9b power and SD (pooled, common intercept)')
mu_path <- function(p, th) {
  br <- gomp(p$tod, 10.98604, 5, 0.42)
  for (i in 2:p$n) if (!p$on[i] && p$tsd[i] > 0) br[i] <- br[i] + br[i - 1] * 0.5^(p$tsd[i] / th)
  gomp(p$w, 6.50647, 5, 0.35) + gomp(p$w, 6.50647, 5, 0.35) * p$e + br
}
e9_pool <- function(rule, th, variant) {
  ps <- lapply(hyb_p, function(o) pinfo(hyb_t, hyb_e, o)); Ns <- alloc(4)
  m <- lapply(ps, function(p) {
    a <- wts(p, extra_off = as.numeric(variant == 'E9b')); C <- CSum(p, 'CS')
    aCa <- drop(t(a) %*% C %*% a)
    mus <- if (th > 0) mu_path(p, th) else { x <- mu_path(p, 1e-9); x }
    if (rule == 'step') {
      g <- as.numeric(p$on | (p$tsd > 0 & th > 0)); gam <- -BETA * sBR / sbm * sum(a * g)
      c(gam, aCa - gam^2 * sbm^2, -sum(a * mus))
    } else {
      h <- as.numeric(p$on); gam <- -BETA * sBR / sbm * sum(a * h)
      c(gam, aCa, -sum(a * mus))
    }
  })
  m <- do.call(rbind, m)
  pooled_power(m[, 1], m[, 2], m[, 3], Ns)
}
for (cs in list(c('step', 0), c('step', 0.1), c('mean', 0), c('mean', 1)))
  for (v in c('E9', 'E9b')) {
    r <- e9_pool(cs[1], as.numeric(cs[2]), v)
    cat(sprintf('%-5s t=%-4s %-4s power %.3f mean %.4f sd %.4f\n', cs[1], cs[2], v,
                r$power, r$mean, r$sd))
  }

## ================= docs/42 ==================
hdr('docs/42: population means and Db-De correlation')
tv <- gomp(hyb_t, 6.50647, 5, 0.35)
print(round(rbind(TV = tv, PB = tv * hyb_e, sum = tv * (1 + hyb_e),
                  BR_path1 = gomp(p1$tod, 10.98604, 5, 0.42)), 2))
cat('drop at week 9:', round(sum(tv[2] * 2) - tv[3] * 1.5, 2), '\n')
cat('G(10) =', round(gomp(10, 10.98604, 5, 0.42), 2), '\n')
Ns <- alloc(4)
rows <- rbindlist(lapply(1:4, function(j) data.table(w = Ns[j], Db = c(0, hyb_p[[j]]),
                                                    De = c(0, hyb_e), t = c(0, hyb_t))))
wc <- function(x, y, w) { mx <- sum(w * x) / sum(w); my <- sum(w * y) / sum(w)
  sum(w * (x - mx) * (y - my)) / sqrt(sum(w * (x - mx)^2) * sum(w * (y - my)^2)) }
r9 <- wc(rows$Db, rows$De, rows$w)
r8 <- rows[t > 0, wc(Db, De, w)]
cat(sprintf('cor(Db, De): 9 rows %.3f (VIF %.2f); 8 rows %.3f (VIF %.2f)\n',
            r9, 1 / (1 - r9^2), r8, 1 / (1 - r8^2)))
X <- rows[, cbind(Db, t)]; fitv <- lm.wfit(cbind(1, X), rows$De, rows$w)
r2 <- 1 - sum(rows$w * fitv$residuals^2) / sum(rows$w * (rows$De - weighted.mean(rows$De, rows$w))^2)
cat(sprintf('VIF of De given Db and t (9 rows): %.2f\n', 1 / (1 - r2)))

## ================= docs/45 S4.3 ==================
hdr('docs/45 S4.3: kappa of a CS fit to AR(1) data in CO')
A <- 0.7^abs(outer(co_t, co_t, '-'))
b1 <- 1:4; b2 <- 5:8
rw <- mean(c(A[b1, b1][upper.tri(A[b1, b1])], A[b2, b2][upper.tri(A[b2, b2])]))
rb <- mean(A[b1, b2]); rall <- mean(A[upper.tri(A)])
vt <- 2 * (1 + 3 * rw) / 4 - 2 * rb; a <- c(rep(.25, 4), rep(-.25, 4))
cat(sprintf('r_within %.3f r_between %.3f r_all %.3f; true %.3f (qf %.3f); CS %.3f; kappa %.3f\n',
            rw, rb, rall, vt, drop(t(a) %*% A %*% a), 2 * (1 - rall) / 4,
            sqrt(2 * (1 - rall) / 4 / vt)))
cat(sprintf('Dbc at 2.5 wk: t_a=0.5 %.3f, t_a=2 %.3f\n', 0.5^(2.5 / 0.5), 0.5^(2.5 / 2)))
cat('MCSE 500 reps', round(sqrt(.05 * .95 / 500), 4), ' 4000 reps', round(sqrt(.05 * .95 / 4000), 4),
    ' band', round(.05 + c(-1, 1) * 1.96 * sqrt(.05 * .95 / 4000), 3), '\n')

## ================= docs/45 S7: plim of beta_bm:D under Dbc ==========
hdr('docs/45 S7: plim of bm:Dbc under an assumed half-life (Hybrid, CS, graded)')
## Rows: baseline (Y0 = BL) then 8 visits. Under covariance binding the
## true slope of Y on B at visit t is -c_bm sBR g_t / sbm; baseline 0.
## Working V = conditional covariance of Y given B (what an unstructured
## fit converges to).
plim_bmD <- function(th_true, th_a, binary = FALSE, V_type = 'true') {
  num <- 0; den <- 0
  ps <- lapply(hyb_p, function(o) pinfo(hyb_t, hyb_e, o)); Ns <- alloc(4)
  M <- matrix(0, 2, 2); v <- c(0, 0)
  for (j in 1:4) {
    p <- ps[[j]]; g <- hprof(p, th_true)
    s <- c(0, -BETA * sBR * g / sbm)
    x <- c(0, if (binary) as.numeric(p$on) else hprof(p, th_a))
    SY <- matrix(sBL^2, 9, 9); SY[2:9, 2:9] <- SY[2:9, 2:9] + CSum(p, 'CS')
    u <- c(0, BETA * sBR * g)
    V <- if (V_type == 'true') SY - outer(u, u) else diag(9)
    Z <- cbind(1, x); Vi <- solve(V)
    M <- M + Ns[j] * t(Z) %*% Vi %*% Z; v <- v + Ns[j] * t(Z) %*% Vi %*% s
  }
  drop(solve(M, v))[2]
}
cat(sprintf('estimand %.4f\n', -BETA * sBR / sbm))
for (cfg in list(c(0, NA, 1), c(1, NA, 1), c(0, 0.5, 0), c(1, 0.5, 0), c(1, 1, 0), c(1, 2, 0),
                 c(0, 2, 0), c(0, 1, 0)))
  cat(sprintf('t_true %.1f  %s  plim %.4f (independence %.4f)\n', cfg[1],
              if (cfg[3] == 1) 'binary D' else sprintf('Dbc t_a=%.1f', cfg[2]),
              plim_bmD(cfg[1], cfg[2], cfg[3] == 1, 'true'),
              plim_bmD(cfg[1], cfg[2], cfg[3] == 1, 'ind')))

## ================= docs/48 S4: reflection symmetry ==================
hdr('docs/48 S4: reflection symmetry')
## (i) Joint covariance of (B, BR block) under -c equals D Sigma(+c) D.
p <- p1; g <- hprof(p, 1); n <- p$n
build_joint <- function(cbm) {
  R <- diag(1 + 3 * n); R[-1, -1] <- rcor(p, 'CS')
  R[1, 1 + 2 * n + (1:n)] <- R[1 + 2 * n + (1:n), 1] <- cbm * g
  s <- c(sbm, rep(sTV, n), sPB * p$e, rep(sBR, n)); outer(s, s) * R
}
Dm <- diag(c(-1, rep(1, 3 * n)))
cat('max |Sigma(-c) - D Sigma(+c) D| =', max(abs(build_joint(-0.2) - Dm %*% build_joint(0.2) %*% Dm)), '\n')
## (ii) Fits with the reflected biomarker give identical p-values.
set.seed(48)
ps <- lapply(hyb_p, function(o) pinfo(hyb_t, hyb_e, o)); Ns <- alloc(4)
dat <- rbindlist(lapply(1:4, function(j) {
  pj <- ps[[j]]
  gj <- hprof(pj, 1); S2 <- diag(1 + 3 * n); S2[-1, -1] <- rcor(pj, 'CS')
  S2[1, 1 + 2 * n + (1:n)] <- S2[1 + 2 * n + (1:n), 1] <- 0.2 * gj
  sj <- c(sbm, rep(sTV, n), sPB * pj$e, rep(sBR, n)); S2 <- outer(sj, sj) * S2
  X <- mvrnorm(Ns[j], rep(0, 1 + 3 * n), S2)
  bl <- rnorm(Ns[j], 0, sBL)
  Y <- bl - (X[, 1 + 1:n] + X[, 1 + n + 1:n] + X[, 1 + 2 * n + 1:n])
  data.table(id = rep(sprintf('%d-%02d', j, seq_len(Ns[j])), each = 9),
             t = rep(c(0, hyb_t), Ns[j]), B = rep(X[, 1] + 124.3, each = 9),
             Db = rep(c(0, hyb_p[[j]]), Ns[j]), Y = as.vector(t(cbind(bl, Y))))
}))
dat[, bmc := B - mean(B)]; dat[, bmr := (2 * mean(B) - B) - mean(2 * mean(B) - B)]
f1 <- lme(Y ~ bmc * Db + t, random = ~ 1 | id, correlation = corCAR1(form = ~ t | id), data = dat)
f2 <- lme(Y ~ bmr * Db + t, random = ~ 1 | id, correlation = corCAR1(form = ~ t | id), data = dat)
tt1 <- summary(f1)$tTable; tt2 <- summary(f2)$tTable
cat(sprintf('lme: est %.5f vs %.5f; p %.8f vs %.8f; logLik diff %.2e\n',
            tt1['bmc:Db', 1], tt2['bmr:Db', 1], tt1['bmc:Db', 5], tt2['bmr:Db', 5],
            logLik(f1) - logLik(f2)))
if (requireNamespace('clubSandwich', quietly = TRUE)) {
  v1 <- clubSandwich::coef_test(f1, vcov = 'CR2', test = 'Satterthwaite')
  v2 <- clubSandwich::coef_test(f2, vcov = 'CR2', test = 'Satterthwaite')
  cat(sprintf('CR2: p %.8f vs %.8f\n', v1$p_Satt[rownames(tt1) == 'bmc:Db'],
              v2$p_Satt[rownames(tt2) == 'bmr:Db']))
}
cat('\nDone.\n')
