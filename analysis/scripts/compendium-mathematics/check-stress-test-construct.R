## Copied 2026-10-08 from the review scratch script inv_d_check.R;
## run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-stress-test-construct.R
## Inventory part d: checks of the equations and numbers carried over
## from docs/43, 44, 46, 47, 50, 52, docs/54 and
## plan-architecture-c.md. Run from the repository root:
##   Rscript analysis/scripts/compendium-mathematics/check-stress-test-construct.R
## Reads analysis/data/quick-sim/carryover-closed-form/ceilings-graded-by-rho.csv
## for comparison only; every value is recomputed here.

## ------------------------------------------------------------------
## Designs (58b32a9 schedules, as in check-ceilings.R) and path info
## ------------------------------------------------------------------
designs <- list(
  CO = list(w = seq(2.5, 20, 2.5),
            paths = list(c(1,1,1,1,0,0,0,0), c(0,0,0,0,1,1,1,1))),
  Hybrid = list(w = c(4, 8, 9, 10, 11, 12, 16, 20),
                paths = list(c(1,1,1,1,0,0,1,0), c(1,1,1,1,0,0,0,1),
                             c(1,1,1,0,0,0,1,0), c(1,1,1,0,0,0,0,1))),
  `OL+BDC` = list(w = c(4, 8, 12, 16, 17, 18, 19, 20),
                  paths = list(c(1,1,1,1,1,1,0,0), c(1,1,1,1,1,0,0,0))))

path_tsd <- function(w, od) {
  dw <- diff(c(0, w))
  n <- length(od)
  tod <- dw * od
  tsd <- dw * (od != 1)
  for (i in 2:n) {
    if (tod[i] > 0) tod[i] <- tod[i] + tod[i - 1] else tsd[i] <- tsd[i] + tsd[i - 1]
  }
  tsd <- (cumsum(od) > 0) * tsd
  list(tod = tod, tsd = tsd)
}

## Coupling patterns. graded: 1 on, 2^(-tsd/th) off after exposure, 0
## before. ratio (2024 vendored): 1 on, 2^(-sf*tsd/th) off after
## exposure, 0 before; at th = 0 both are 0 off drug.
coupling <- function(od, tsd, th, rule = 'graded', sf = 2, occ1 = TRUE) {
  off <- if (th <= 0) rep(0, length(od)) else
    if (rule == 'graded') 2^(-tsd / th) else 2^(-sf * tsd / th)
  g <- ifelse(od == 1, 1, ifelse(tsd > 0, off, 0))
  if (!occ1) g[1] <- 0
  g
}

## Response block M (component-major TV, PB, BR).
build_M <- function(w, rho, c1 = 0.2, cx = 0.1, within = 'ar1', cross = 'sep') {
  n <- length(w)
  lag <- abs(outer(w, w, '-'))
  A <- if (within == 'ar1') rho^lag else { a <- matrix(rho, n, n); diag(a) <- 1; a }
  C <- switch(cross,
    sep = c1 * rho^lag,
    decay = { m <- cx * rho^lag; diag(m) <- c1; m },
    const = { m <- matrix(cx, n, n); diag(m) <- c1; m })
  rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
}
ceil <- function(M, g) {
  u <- c(rep(0, 2 * length(g)), g)
  1 / sqrt(drop(crossprod(u, solve(M, u))))
}
min_eig <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

design_ceiling <- function(dn, rho, th, within, cross, rule = 'graded',
                           sf = 2, occ1 = TRUE, c1 = 0.2, cx = 0.1) {
  d <- designs[[dn]]
  M <- build_M(d$w, rho, c1, cx, within, cross)
  vals <- sapply(d$paths, function(od) {
    ps <- path_tsd(d$w, od)
    ceil(M, coupling(od, ps$tsd, th, rule, sf, occ1))
  })
  c(ceiling = min(vals), min_eig_M = min_eig(M))
}

## ------------------------------------------------------------------
## 1. docs/44 Step 6: graded ceilings by rho for CS, F (decaying cross),
##    C (separable); minimum over designs; compared with the stored CSV.
## ------------------------------------------------------------------
cat('== 1. Graded-coupling ceilings by rho (docs/44, Step 6) ==\n')
cons <- list(`CS` = c('cs', 'const'), `F` = c('ar1', 'decay'), `C` = c('ar1', 'sep'))
ths <- c(0, 0.5, 1, 2)
res <- expand.grid(rho = c(0.8, 0.7, 0.6, 0.5), con = names(cons),
                   design = names(designs), th = ths, stringsAsFactors = FALSE)
res$ceiling <- NA_real_
res$minEigM <- NA_real_
for (i in seq_len(nrow(res))) {
  cc <- cons[[res$con[i]]]
  v <- design_ceiling(res$design[i], res$rho[i], res$th[i], cc[1], cc[2])
  res$ceiling[i] <- v['ceiling']
  res$minEigM[i] <- v['min_eig_M']
}
for (r in c(0.8, 0.7, 0.6, 0.5)) {
  row <- sapply(names(cons), function(k) min(res$ceiling[res$rho == r & res$con == k]))
  argm <- sapply(names(cons), function(k) {
    s <- res[res$rho == r & res$con == k, ]
    paste(s$design[which.min(s$ceiling)], s$th[which.min(s$ceiling)])
  })
  cat(sprintf('rho %.1f  CS %.4f (%s)  F %.4f (%s)  C %.4f (%s)\n', r,
              row['CS'], argm['CS'], row['F'], argm['F'], row['C'], argm['C']))
}
cat('\nPer design, t_half = 0:\n')
for (r in c(0.8, 0.7, 0.6, 0.5)) for (k in names(cons)) {
  s <- res[res$rho == r & res$con == k & res$th == 0, ]
  cat(sprintf('rho %.1f %-2s  CO %.4f  Hybrid %.4f  OL+BDC %.4f\n', r, k,
              s$ceiling[s$design == 'CO'], s$ceiling[s$design == 'Hybrid'],
              s$ceiling[s$design == 'OL+BDC']))
}
cat('\nrho 0.5, CO, carryover vs none (F and C):\n')
for (k in c('F', 'C')) {
  s <- res[res$rho == 0.5 & res$con == k & res$design == 'CO', ]
  cat(sprintf('%s: %s\n', k, paste(sprintf('th=%g %.4f', s$th, s$ceiling), collapse = '  ')))
}
cat('\nCells where carryover lowers the ceiling below t_half = 0:\n')
for (i in seq_len(nrow(res))) if (res$th[i] > 0) {
  b <- res$ceiling[res$rho == res$rho[i] & res$con == res$con[i] &
                   res$design == res$design[i] & res$th == 0]
  if (res$ceiling[i] < b - 1e-9)
    cat(sprintf('  rho %.1f %s %s th %g: %.4f vs %.4f (diff %.4f)\n', res$rho[i],
                res$con[i], res$design[i], res$th[i], res$ceiling[i], b, b - res$ceiling[i]))
}
csv <- 'analysis/data/quick-sim/carryover-closed-form/ceilings-graded-by-rho.csv'
if (file.exists(csv)) {
  st <- read.csv(csv)
  key <- c(`CS (published)` = 'CS', `AR(1) covar (F)` = 'F', `Separable AR(1) (C)` = 'C')
  st$con <- key[st$construct]
  mm <- merge(st, res, by.x = c('rho', 'con', 'design', 't_half'),
              by.y = c('rho', 'con', 'design', 'th'))
  cat(sprintf('\nStored CSV rows matched: %d of %d; max |recomputed - stored| %.2e\n',
              nrow(mm), nrow(st), max(abs(mm$ceiling.x - mm$ceiling.y))))
}
cat(sprintf('Stress-test cells (Hybrid, t_half 0.5): CS rho 0.8 %.4f; C rho 0.7 %.4f\n',
            res$ceiling[res$rho == 0.8 & res$con == 'CS' & res$design == 'Hybrid' & res$th == 0.5],
            res$ceiling[res$rho == 0.7 & res$con == 'C' & res$design == 'Hybrid' & res$th == 0.5]))

## OL design, constant pattern (docs/44 Step 4: 0.893 -> 0.576)
cat('\n== OL constant pattern: CS vs AR(1) decaying cross (docs/44 Step 4) ==\n')
for (wn in c('Hybrid', 'OL+BDC')) for (r in c(0.8, 0.7)) {
  w <- designs[[wn]]$w
  g <- rep(1, 8)
  cat(sprintf('weeks of %-6s rho %.1f: CS %.4f  F %.4f  C %.4f; minEig F %.4f\n', wn, r,
              ceil(build_M(w, r, within = 'cs', cross = 'const'), g),
              ceil(build_M(w, r, within = 'ar1', cross = 'decay'), g),
              ceil(build_M(w, r, within = 'ar1', cross = 'sep'), g),
              min_eig(build_M(w, r, within = 'ar1', cross = 'decay'))))
}
cat('Smallest eigenvalue of the F block over designs and rho:\n')
for (dn in names(designs)) for (r in c(0.8, 0.7))
  cat(sprintf('  %-6s rho %.1f  %.4f\n', dn, r,
              min_eig(build_M(designs[[dn]]$w, r, within = 'ar1', cross = 'decay'))))

## ------------------------------------------------------------------
## 2. Covar whitepaper Section 5: 16-configuration decomposition at
##    rho 0.7, c1 0.2, cx 0.1, t_half in {0, 0.5, 1}; occasion 1.
## ------------------------------------------------------------------
cat('\n== 2. Covar whitepaper decomposition (rho 0.7) ==\n')
cfg <- expand.grid(ar1 = c(FALSE, TRUE), cfd = c(FALSE, TRUE),
                   bexp = c(FALSE, TRUE), occ1 = c(FALSE, TRUE))
cell_ceil <- function(ar1, cfd, bexp, occ1, dn, th) {
  within <- if (ar1) 'ar1' else 'cs'
  cross <- if (cfd) 'decay' else 'const'
  rule <- if (bexp) 'graded' else 'ratio'
  design_ceiling(dn, 0.7, th, within, cross, rule, sf = 2, occ1 = occ1)['ceiling']
}
cfg$worst <- NA_real_
for (i in seq_len(nrow(cfg))) {
  v <- c()
  for (dn in names(designs)) for (th in c(0, 0.5, 1))
    v <- c(v, cell_ceil(cfg$ar1[i], cfg$cfd[i], cfg$bexp[i], cfg$occ1[i], dn, th))
  cfg$worst[i] <- min(v)
}
cfg <- cfg[order(-cfg$worst), ]
print(transform(cfg, worst = round(worst, 4), grid001 = floor(worst * 100 + 1e-9) / 100),
      row.names = FALSE)
cat('\nPer design and t_half: orig(all off) | AR1 only | AR1+occ1 | full covar\n')
for (dn in names(designs)) for (th in c(0, 0.5, 1)) {
  cat(sprintf('%-6s th %.1f  %.4f  %.4f  %.4f  %.4f\n', dn, th,
              cell_ceil(FALSE, FALSE, FALSE, FALSE, dn, th),
              cell_ceil(TRUE, FALSE, FALSE, FALSE, dn, th),
              cell_ceil(TRUE, FALSE, FALSE, TRUE, dn, th),
              cell_ceil(TRUE, TRUE, TRUE, TRUE, dn, th)))
}
cat('\nSmallest eigenvalue of the response block M per configuration (rho 0.7):\n')
for (ar in c(FALSE, TRUE)) for (cf in c(FALSE, TRUE)) {
  me <- sapply(names(designs), function(dn)
    min_eig(build_M(designs[[dn]]$w, 0.7, within = if (ar) 'ar1' else 'cs',
                    cross = if (cf) 'decay' else 'const')))
  cat(sprintf('  ar1 %-5s cfd %-5s  %s\n', ar, cf,
              paste(sprintf('%s %.4f', names(me), me), collapse = '  ')))
}
cat('\nCO under orig (CS, const cross, ratio coupling, occ1 off), per path:\n')
for (th in c(0, 0.5, 1)) {
  d <- designs$CO
  M <- build_M(d$w, 0.7, within = 'cs', cross = 'const')
  v <- sapply(d$paths, function(od) {
    ps2 <- path_tsd(d$w, od)
    ceil(M, coupling(od, ps2$tsd, th, 'ratio', 2, FALSE))
  })
  cat(sprintf('  th %.1f  drug-first %.4f  placebo-first %.4f\n', th, v[1], v[2]))
}
cat('\nOne flag on from orig (worst cell):\n')
for (f in c('none', 'ar1', 'cfd', 'bexp', 'occ1')) {
  fl <- c(ar1 = FALSE, cfd = FALSE, bexp = FALSE, occ1 = FALSE)
  if (f != 'none') fl[f] <- TRUE
  w <- cfg$worst[cfg$ar1 == fl['ar1'] & cfg$cfd == fl['cfd'] &
                 cfg$bexp == fl['bexp'] & cfg$occ1 == fl['occ1']]
  cat(sprintf('  %-5s %.4f\n', f, w))
}

## ------------------------------------------------------------------
## 3. Covar whitepaper Sections 3.3 and 3.5: ratio-graded coupling and
##    the scalefactor recursion on the BR mean, Hybrid path A, t_half 1.
## ------------------------------------------------------------------
cat('\n== 3. Ratio-graded coupling and scalefactor recursion ==\n')
gomp <- function(t, maxr, disp, rate) {
  y <- maxr * exp(-disp * exp(-rate * t))
  off <- maxr * exp(-disp)
  (y - off) * maxr / (maxr - off)
}
od <- designs$Hybrid$paths[[1]]
ps <- path_tsd(designs$Hybrid$w, od)
raw <- gomp(ps$tod, 10.98604, 5, 0.42)
recur <- function(sf, th) {
  m <- raw
  for (p in 2:8) if (od[p] == 0 && ps$tsd[p] > 0) m[p] <- m[p] + m[p - 1] * 2^(-sf * ps$tsd[p] / th)
  m
}
anch <- raw
last <- NA
for (p in 1:8) if (od[p] == 1) last <- raw[p] else if (ps$tsd[p] > 0) anch[p] <- last * 2^(-ps$tsd[p] / 1)
m1 <- recur(1, 1); m2 <- recur(2, 1)
cpl <- ifelse(od == 1, 0.45, ifelse(m2 == 0, 0, c(NA, m2[-1] / m2[-8]) * 0.45))
cpl[1] <- NA
tab <- data.frame(visit = c('OL1','OL2','BD1','BD2','BD3','BD4','COd','COp'),
                  on = od, tsd = ps$tsd, sf1 = round(m1, 3), sf2 = round(m2, 3),
                  anchored = round(anch, 3), ratio_sf2 = round(m2 / anch, 4),
                  ratio_sf1 = round(m1 / anch, 4), coupling = round(cpl, 4),
                  closed = round(ifelse(od == 1, 0.45, 0.45 * 2^(-2 * ps$tsd)), 4))
print(tab, row.names = FALSE)
cat(sprintf('graded at 1 week, t_half 1: %.4f\n', 0.45 * exp(-log(2))))
cat('General ratio formula 2^-(sf*sum tsd - tsd_k)/th at BD4 (sf 1):',
    2^(-(1 * (1 + 2) - 2)), ' COp (sf 2):', 2^(-(2 * 4 - 4)), '\n')

## ------------------------------------------------------------------
## 4. Covar whitepaper: CS exact ceilings at rho 0.7 under the 2024
##    coupling with occasion 1 uncoupled (0.346, 0.351).
## ------------------------------------------------------------------
cat('\n== 4. 2024 orig under CS (occ1 off, ratio coupling, const cross) ==\n')
for (dn in names(designs)) for (th in c(0, 0.5, 1))
  cat(sprintf('%-6s th %.1f  %.4f\n', dn, th, cell_ceil(FALSE, FALSE, FALSE, FALSE, dn, th)))

## ------------------------------------------------------------------
## 5. Phase-specific variances implied by configuration C (stress test),
##    Hybrid path A, rho 0.7; RI + varIdent(phase) + CAR1 comparison.
## ------------------------------------------------------------------
cat('\n== 5. Implied outcome covariance by phase, configuration C ==\n')
sd_bl <- 18.48267
w <- designs$Hybrid$w
e <- c(1, 1, .5, .5, .5, .5, .5, .5)
K <- matrix(0.2, 3, 3); diag(K) <- 1
A <- 0.7^abs(outer(w, w, '-'))
S <- cbind(rep(10, 8), 10 * e, rep(8, 8))
CovY <- matrix(0, 9, 9)
CovY[1, ] <- CovY[, 1] <- sd_bl^2
for (t in 1:8) for (s in 1:8)
  CovY[t + 1, s + 1] <- sd_bl^2 + drop(S[t, ] %*% K %*% S[s, ]) * A[t, s]
cat(sprintf('sigma_BL^2 %.1f; Var(Y) baseline %.1f, OL %.1f, blinded %.1f\n', sd_bl^2,
            CovY[1, 1], CovY[2, 2], CovY[4, 4]))
q_ol <- drop(S[1, ] %*% K %*% S[1, ]); q_bl <- drop(S[3, ] %*% K %*% S[3, ])
q_x <- drop(S[1, ] %*% K %*% S[3, ])
cat(sprintf('s^T K s: OL %.1f blinded %.1f cross %.1f; factor %.4f\n', q_ol, q_bl, q_x,
            q_x / sqrt(q_ol * q_bl)))
cat(sprintf('with measurement error (3, 6): OL %.1f blinded %.1f\n',
            CovY[2, 2] + 9, CovY[4, 4] + 36))
## The RI + varIdent + CAR1 structure with sigma_u^2 = sigma_BL^2,
## sigma_e,0 -> 0, sigma_e,OL^2 = q_ol, sigma_e,bl^2 = q_bl, rho_a = 0.7:
sig <- c(0, sqrt(q_ol), sqrt(q_ol), rep(sqrt(q_bl), 6))
tt <- c(0, w)
Vm <- sd_bl^2 + outer(sig, sig) * 0.7^abs(outer(tt, tt, '-'))
cat(sprintf('max |CovY - V_model| over all pairs %.3f; over pairs not OL x blinded %.2e\n',
            max(abs(CovY - Vm)),
            max(abs((CovY - Vm)[c(1:3), c(1:3)]), abs((CovY - Vm)[4:9, 4:9]),
                abs((CovY - Vm)[1, ]))))

## Leak factors plim(beta_bm:D)/s_prog by the separated GLS normal
## equations (MAIN eq:bp-plim), configuration C, rho 0.7, Hybrid,
## N = 70 (18/18/17/17), X = D_bc at t_half^a = 0.5.
cat('\nLeak factors, configuration C (rho 0.7), working V:\n')
alloc0 <- c(18, 18, 17, 17)
Vind <- diag(9)
Vlist <- list(true = CovY, RCH_limit = Vm, independence = Vind)
leak <- function(V, zfun, s) {
  Aacc <- 0; bacc <- 0
  for (k in 1:4) {
    od <- designs$Hybrid$paths[[k]]
    tsd <- c(0, path_tsd(w, od)$tsd)
    Xk <- c(0, ifelse(od == 1, 1, ifelse(tsd[-1] > 0, 0.5^(tsd[-1] / 0.5), 0)))
    Z <- zfun(Xk)
    Vi <- solve(V)
    Aacc <- Aacc + alloc0[k] * t(Z) %*% Vi %*% Z
    bacc <- bacc + alloc0[k] * t(Z) %*% Vi %*% s
  }
  drop(solve(Aacc, bacc))[2]
}
eta9 <- c(0, e); iota <- c(1, rep(0, 8))
zS <- function(X) cbind(1, X)
zE <- function(X) cbind(1, X, iota)
zF <- function(X) cbind(1, X, eta9)
s_tv <- c(0, rep(1, 8)); s_pb <- eta9
for (nm in names(Vlist)) cat(sprintf(
  '  %-12s TV const: S %.3f ES %.3f F1 %.3f | placebo: S %.3f ES %.3f F1 %.3f\n', nm,
  leak(Vlist[[nm]], zS, s_tv), leak(Vlist[[nm]], zE, s_tv), leak(Vlist[[nm]], zF, s_tv),
  leak(Vlist[[nm]], zS, s_pb), leak(Vlist[[nm]], zE, s_pb), leak(Vlist[[nm]], zF, s_pb)))

## ------------------------------------------------------------------
## 6. Correlation of the design expectancy with the binary drug state,
##    Hybrid, 9 rows per participant, allocation 18/18/17/17 (docs/40
##    via docs/43: 0.59).
## ------------------------------------------------------------------
cat('\n== 6. Cor(eta_t, D_it), Hybrid ==\n')
alloc <- c(18, 18, 17, 17)
etas <- c(); ds <- c()
for (k in 1:4) {
  etas <- c(etas, rep(c(0, e), alloc[k]))
  ds <- c(ds, rep(c(0, designs$Hybrid$paths[[k]]), alloc[k]))
}
cat(sprintf('with baseline row %.4f; post-baseline only %.4f\n', cor(etas, ds),
            cor(etas[rep(c(FALSE, rep(TRUE, 8)), 70)], ds[rep(c(FALSE, rep(TRUE, 8)), 70)])))

## ------------------------------------------------------------------
## 7. Dual channel: exact conditional moments given the biomarker.
## ------------------------------------------------------------------
cat('\n== 7. Dual-channel conditional moments ==\n')
sb <- 8; cA <- 0.2; cB <- 0.25
g <- coupling(od, ps$tsd, 1, 'graded')
SigBR <- sb^2 * 0.7^abs(outer(w, w, '-'))
## joint of (b, BR draw): Cov(b, BR_t) = cB g_t sb; BR_C = BR + cA b sb D_t
covbB <- cB * g * sb + cA * sb * od
VarBRC <- SigBR + sb^2 * (cA^2 * outer(od, od) + cA * cB * (outer(od, g) + outer(g, od)))
condC <- VarBRC - outer(covbB, covbB)
cat(sprintf('slope check max |Cov(b,BR_C) - sb(cA D + cB g)| %.2e\n',
            max(abs(covbB - sb * (cA * od + cB * g)))))
cat(sprintf('conditional cov check max |cond - (SigBR - cB^2 sb^2 gg^T)| %.2e\n',
            max(abs(condC - (SigBR - cB^2 * sb^2 * outer(g, g))))))
set.seed(11)
nn <- 4e5
bz <- rnorm(nn)
Lc <- chol(SigBR - cB^2 * sb^2 * outer(g, g))
brd <- outer(bz, cB * sb * g) + matrix(rnorm(nn * 8), nn) %*% Lc
brc <- brd + outer(bz, cA * sb * od)
cat(sprintf('simulated on-drug slope visit 1 %.3f (target %.3f); off visit 5 %.3f (target %.3f)\n',
            coef(lm(brc[, 1] ~ bz))[2], sb * (cA + cB),
            coef(lm(brc[, 5] ~ bz))[2], sb * cB * g[5]))
cat(sprintf('marginal on-drug variance %.2f, formula %.2f\n', var(brc[, 1]),
            sb^2 * (1 + cA^2 + 2 * cA * cB)))

## ------------------------------------------------------------------
## 8. Sign symmetry: reflected biomarker gives identical p-values.
## ------------------------------------------------------------------
cat('\n== 8. Reflection invariance of S and F1 p-values ==\n')
suppressMessages({library(nlme); library(clubSandwich)})
set.seed(3)
np <- 40
id <- rep(seq_len(np), each = 9)
tt9 <- rep(c(0, w), np)
path <- rep(rep(1:4, length.out = np), each = 9)
Db <- unlist(lapply(rep(1:4, length.out = np), function(k) c(0, designs$Hybrid$paths[[k]])))
tsdv <- unlist(lapply(rep(1:4, length.out = np), function(k)
  c(0, path_tsd(designs$Hybrid$w, designs$Hybrid$paths[[k]])$tsd)))
X <- ifelse(Db == 1, 1, ifelse(tsdv > 0, 0.5^(tsdv / 0.5), 0))
De <- rep(c(0, e), np)
bm <- rep(rnorm(np, 124, 15), each = 9)
Y <- 80 + rep(rnorm(np, 0, 15), each = 9) - 5 * X - 0.1 * (bm - 124) * X + rnorm(9 * np, 0, 8)
dat <- data.frame(ptID = factor(id), t = tt9, X = X, De = De, Sx = Y,
                  bmc = bm - mean(bm[tt9 == 0]))
dref <- dat; dref$bmc <- -dat$bmc
fitp <- function(d, f) {
  m <- lme(f, random = ~ 1 | ptID, correlation = corCAR1(form = ~ t | ptID), data = d)
  ct <- coef_test(m, vcov = 'CR2', test = 'Satterthwaite')
  k <- which(ct$Coef == 'bmc:X')
  c(est = ct$beta[k], p_cr2 = ct$p_Satt[k], p_mb = summary(m)$tTable['bmc:X', 'p-value'])
}
for (f in list(Sx ~ bmc + t + X + bmc:X, Sx ~ bmc + t + X + bmc:X + De + bmc:De)) {
  a <- fitp(dat, f); b <- fitp(dref, f)
  cat(sprintf('est %.5f vs %.5f; CR2 p %.6f vs %.6f; model p %.6f vs %.6f\n',
              a['est'], b['est'], a['p_cr2'], b['p_cr2'], a['p_mb'], b['p_mb']))
}

## ------------------------------------------------------------------
## 9. Monte Carlo bands and MCSEs of docs/46 and docs/47.
## ------------------------------------------------------------------
cat('\n== 9. MCSE bands ==\n')
band <- function(R) 0.05 + c(-2, 2) * sqrt(0.05 * 0.95 / R)
cat(sprintf('cell R=1000: MCSE %.4f\n', sqrt(0.05 * 0.95 / 1000)))
cat(sprintf('pooled 11,000: band %.4f to %.4f\n', band(11000)[1], band(11000)[2]))
cat(sprintf('docs/45 8 x 500 = 4,000: band %.4f to %.4f\n', band(4000)[1], band(4000)[2]))
cat(sprintf('power, R=500, worst case: %.4f; difference of two cells: %.4f\n',
            sqrt(0.25 / 500), sqrt(2 * 0.25 / 500)))
