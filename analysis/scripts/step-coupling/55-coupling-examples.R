## docs/55: step coupling of the biomarker in the published pmsimstats code.
##
## Regenerates every table, figure and verified number in
## docs/55-step-coupling-in-published-code.md except the lmer power
## simulation, which is 55-power-check.R (its stored output is read here
## for Figure 2 when present). Run from the repository root:
##   Rscript analysis/scripts/step-coupling/55-coupling-examples.R [out_dir]
##
## Sources: the 58b32a9 files under analysis/data/quick-sim/
## hendrickson-58b32a9/R (byte-identical to github.com/rchendrickson/
## pmsimstats at 58b32a9), and the vendored 2024 upstream generateData()
## (upstream master 06dac83 up to one semantically identical line).
##
## Method. Each version of generateData() is run with mvrnorm() replaced
## by a function that records the mean vector and covariance matrix it is
## given, so every coupling, mean and matrix below is the code's own.

suppressMessages({
  library(data.table)
  library(MASS)
  library(corpcor)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else file.path(tempdir(), 'docs55')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- 'docs/figures'
src <- 'analysis/data/quick-sim/hendrickson-58b32a9/R'
head_file <- paste0('analysis/scripts/quick-sim/',
  'hendrickson-original-comparison/vendored-hendrickson-generateData.R')
sim_dir <- 'analysis/data/quick-sim/step-coupling'
stored14 <- paste0('analysis/data/quick-sim/hendrickson-problems/',
  'graded-patch/graded-patch-reps500.csv')
stopifnot(file.exists(file.path(src, 'generateData.R')), file.exists(head_file))
cat('md5 generateData.R (58b32a9):',
  tools::md5sum(file.path(src, 'generateData.R')), '\n')

## ---- versions of generateData() ---------------------------------------
gd <- readLines(file.path(src, 'generateData.R'))

## Graded coupling: lines 134-137 replaced (report 01, Appendix B.1;
## identical to 14-graded-patch-58b32a9.R).
graded_patch <- c(
  '      if(d[p]$tod>0){',
  '        g<-1',
  '      }else if(d[p]$tsd>0 && modelparam$carryover_t1half>0){',
  '        g<-(1/2)^(d[p]$tsd/modelparam$carryover_t1half)',
  '      }else{',
  '        g<-0',
  '      }',
  '      correlations[n1,\'bm\']<-modelparam$c.bm*g',
  '      correlations[\'bm\',n1]<-modelparam$c.bm*g')
apply_graded <- function(lines) {
  i <- grep('if(means[which(n1==labels)]!=0){', lines, fixed = TRUE)
  stopifnot(length(i) == 1, grepl('^\\s*}\\s*$', lines[i + 3]))
  c(lines[seq_len(i - 1)], graded_patch, lines[(i + 4):length(lines)])
}
## Anchored carryover of the mean (docs/53, Section 9.2).
apply_anchored <- function(lines) {
  i <- grep('brmeans[p]<-brmeans[p]+brmeans[p-1]*', lines, fixed = TRUE)
  j <- grep('brmeans<-modgompertz(d$tod', lines, fixed = TRUE)
  stopifnot(length(i) == 1, length(j) == 1, i - j == 5)
  new <- c(
    '      brmeans<-modgompertz(d$tod,rp$max,rp$disp,rp$rate)',
    '      if(nP>1){',
    '        last_on<-0',
    '        for(p in 1:nP){',
    '          if(d[p]$onDrug){',
    '            last_on<-brmeans[p]',
    '          }else if(d[p]$tsd>0){',
    paste0('            brmeans[p]<-brmeans[p]+last_on*',
      '(1/2)^(d$tsd[p]/modelparam$carryover_t1half)'),
    '          }',
    '        }',
    '      }')
  c(lines[seq_len(j - 1)], new, lines[(j + 10):length(lines)])
}
stopifnot(grepl('if(means[which(n1==labels)]!=0){', gd[134], fixed = TRUE))
gd_graded <- apply_graded(gd)

## The two fixes touch disjoint lines; check that the order of application
## does not matter.
order_ok <- identical(apply_anchored(apply_graded(gd)),
  apply_graded(apply_anchored(gd)))
cat('graded+anchored identical in either order:', order_ok, '\n')

## The 2024 upstream file, and the same file with anchored means (to test
## whether its ratio coupling survives the docs/53 fix).
up <- readLines(head_file)
iu <- grep('brmeans[p]<-brmeans[p]+brmeans[p-1]*(1/2)^(scalefactor',
  up, fixed = TRUE)
stopifnot(length(iu) == 1, grepl('if(nP>1){', up[iu - 4], fixed = TRUE),
  grepl('^\\s*}\\s*$', up[iu + 4]))
up_anch <- c(up[seq_len(iu - 5)],
  '      if(nP>1){',
  '        last_on<-0',
  '        for(p in 1:nP){',
  '          if(d[p]$onDrug){',
  '            last_on<-brmeans[p]',
  '          }else if(d[p]$tsd>0){',
  paste0('            brmeans[p]<-brmeans[p]+last_on*',
    '(1/2)^(scalefactor*d$tsd[p]/modelparam$carryover_t1half)'),
  '          }',
  '        }',
  '      }',
  up[(iu + 5):length(up)])

write_src <- function(lines, name) {
  f <- file.path(out_dir, name)
  writeLines(lines, f)
  f
}
files <- c(
  published = file.path(src, 'generateData.R'),
  graded = write_src(gd_graded, 'generateData-graded.R'),
  upstream = head_file,
  upstream_anchored = write_src(up_anch, 'generateData-2024-anchored.R'))
writeLines(c(gd[131:138]), file.path(out_dir, 'published-lines-131-138.txt'))

make_env <- function(gen_file) {
  env <- new.env()
  for (f in c('utilities.R', 'buildtrialdesign.R')) {
    sys.source(file.path(src, f), envir = env)
  }
  sys.source(gen_file, envir = env)
  env$.cap <- NULL
  env$mvrnorm <- function(n, mu, Sigma, empirical = FALSE, ...) {
    env$.cap <- list(mu = mu, Sigma = Sigma)
    matrix(0, n, length(mu))
  }
  env
}
envs <- lapply(files, make_env)

load('data/extracted_rp.rda')
load('data/extracted_bp.rda')
sd_BR <- extracted_rp[cat == 'br']$sd
sd_bm <- extracted_bp[cat == 'bm']$sd
mu_bm <- extracted_bp[cat == 'bm']$m
cat('sd_BR', sd_BR, 'sd_bm', sd_bm, 'mu_bm', mu_bm, '\n')

bt <- envs$published$buildtrialdesign
designs <- list(
  Hybrid = bt(name_longform = 'primary N-of-1 design',
    name_shortform = 'N-of-1', timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
    timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
    expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
    ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
      pathB = c(1, 1, 1, 1, 0, 0, 0, 1), pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
      pathD = c(1, 1, 1, 0, 0, 0, 0, 1)))$trialpaths,
  'OL+BDC' = bt(name_longform = 'open label+blinded discontinuation',
    name_shortform = 'OL+BDC', timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
    timeptnames = c('OL1', 'OL2', 'OL3', 'OL4', 'BD1', 'BD2', 'BD3', 'BD4'),
    expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
    ondrug = list(pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
      pathB = c(1, 1, 1, 1, 1, 0, 0, 0)))$trialpaths,
  CO = bt(name_longform = 'traditional crossover', name_shortform = 'CO',
    timepoints = envs$published$cumulative(rep(2.5, 8)),
    timeptnames = c(paste('COa', 1:4, sep = ''), paste('COb', 1:4, sep = '')),
    expectancies = rep(.5, 8),
    ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
      pathB = c(0, 0, 0, 0, 1, 1, 1, 1)))$trialpaths)

## Published correlation parameters (results_core modelparams).
mparams <- function(cbm, th, rho = 0.8) {
  data.frame(N = 2, c.bm = cbm, carryover_t1half = th, c.tv = rho,
    c.pb = rho, c.br = rho, c.cf1t = 0.2, c.cfct = 0.1)
}

## Mean vector and covariance matrix of one path, as the code builds them
## (no positive-definiteness repair).
capture <- function(version, td, cbm, th, rho = 0.8, ...) {
  env <- envs[[version]]
  env$.cap <- NULL
  invisible(suppressWarnings(env$generateData(mparams(cbm, th, rho),
    extracted_rp, extracted_bp, td, empirical = FALSE,
    makePositiveDefinite = FALSE, ...)))
  env$.cap
}
idx <- function(n) list(bm = 1, BL = 2, tv = 3:(2 + n), pb = (3 + n):(2 + 2 * n),
  br = (3 + 2 * n):(2 + 3 * n))

## Coupling g_t: Cor(bm, BR_t) / c_bm.
coupling_of <- function(cap, cbm) {
  n <- (length(cap$mu) - 2) / 3
  ii <- idx(n)
  R <- cov2cor(cap$Sigma)
  R[1, ii$br] / cbm
}
br_mean_of <- function(cap) {
  n <- (length(cap$mu) - 2) / 3
  cap$mu[idx(n)$br]
}

## ---- closed form of E9 / E9b from the captured matrix ------------------
## Delta = sum_t a_t Y_t, Y_t = BL - (TV_t + PB_t + BR_t); E9b adds the
## baseline (Y_0 = BL) as one more off-drug row. Exact per-path moments,
## pooled as in docs/37, Section 4.4.
path_moments <- function(cap, td, variant) {
  n <- nrow(td)
  ii <- idx(n)
  on <- td$tod > 0
  n_on <- sum(on)
  n_off <- sum(!on) + (variant == 'E9b')
  a <- ifelse(on, 1 / n_on, -1 / n_off)
  a0 <- if (variant == 'E9b') -1 / n_off else 0
  w <- numeric(length(cap$mu))
  w[ii$BL] <- sum(a) + a0
  w[ii$tv] <- -a
  w[ii$pb] <- -a
  w[ii$br] <- -a
  S <- cap$Sigma
  var_D <- drop(t(w) %*% S %*% w)
  cov_DB <- drop(t(w) %*% S[, 1])
  gamma <- cov_DB / S[1, 1]
  data.table(gamma = gamma, var_D = var_D, tau2 = var_D - gamma^2 * S[1, 1],
    mu = sum(w * cap$mu), n_on = n_on, n_off = sum(!on))
}
ALPHA <- 0.05
N <- 70L
pooled <- function(m) {
  nP <- nrow(m)
  Ns <- N %/% nP + c(rep(1L, N %% nP), rep(0L, nP - N %% nP))
  w <- Ns / N
  chi_q <- qchisq((seq_len(400) - 0.5) / 400, N - 1)
  gbar <- sum(w * m$gamma)
  mubar <- sum(w * m$mu)
  resid <- sum(w * (m$tau2 + (m$mu - mubar)^2 + sd_bm^2 * (m$gamma - gbar)^2))
  tc <- qt(1 - ALPHA / 2, N - 2)
  ncp <- gbar * sd_bm * sqrt(chi_q) / sqrt(resid)
  power <- mean(pt(-tc, N - 2, ncp = ncp) +
    pt(tc, N - 2, ncp = ncp, lower.tail = FALSE))
  data.table(mean_slope = gbar, power = power)
}
closed <- function(version, design, cbm, th, variant, ...) {
  tds <- designs[[design]]
  m <- rbindlist(lapply(tds, function(td)
    path_moments(capture(version, td, cbm, th, ...), td, variant)))
  pooled(m)
}

## Ceiling: the largest c_bm for which the correlation matrix is positive
## definite, min over paths of (u' M^-1 u)^(-1/2) (docs/36, Section 4).
ceiling_of <- function(version, design, th, rho = 0.8, ...) {
  min(vapply(designs[[design]], function(td) {
    cap <- capture(version, td, 0.25, th, rho, ...)
    R <- cov2cor(cap$Sigma)
    u <- R[-1, 1] / 0.25
    M <- R[-1, -1]
    if (!any(u != 0)) return(Inf)
    1 / sqrt(drop(crossprod(u, solve(M, u))))
  }, numeric(1)))
}

CBM <- 0.25
k_slope <- CBM * sd_BR / sd_bm
cat(sprintf('c_bm sigma_BR / sigma_bm = %.6f\n', k_slope))

## ---- 1. The published lines and the mean trigger ------------------------
cat('\nPublished lines 131-138:\n')
writeLines(sprintf('%4d  %s', 131:138, gd[131:138]))

## ---- 2. Worked example 1: Hybrid path 1 --------------------------------
td1 <- designs$Hybrid[[1]]
ex1 <- rbindlist(lapply(c(0.1, 1), function(th) {
  cp <- capture('published', td1, CBM, th)
  cg <- capture('graded', td1, CBM, th)
  cu <- capture('upstream', td1, CBM, th)
  stopifnot(identical(cp$mu, cg$mu))
  data.table(t_half = th, visit = td1$timeptnames,
    week = envs$published$cumulative(td1$t_wk), D = as.integer(td1$tod > 0),
    t_od = td1$tod, t_sd = td1$tsd, mu_pub = br_mean_of(cp),
    nonzero = br_mean_of(cp) != 0, g_step = coupling_of(cp, CBM),
    g_graded = coupling_of(cg, CBM), g_2024 = coupling_of(cu, CBM),
    mu_2024 = br_mean_of(cu))
}))
ex1[, slope_step := k_slope * g_step]
ex1[, slope_graded := k_slope * g_graded]
ex1[, g_formula := ifelse(D == 1, 1, ifelse(t_sd > 0, 2^(-t_sd / t_half), 0))]
stopifnot(max(abs(ex1$g_graded - ex1$g_formula)) < 1e-12)
fwrite(ex1, file.path(out_dir, 'example1.csv'))
cat('\nExample 1 (Hybrid path 1)\n')
print(ex1[, .(t_half, visit, week, D, t_sd, mu_pub = signif(mu_pub, 4),
  nonzero, g_step, g_graded = signif(g_graded, 4),
  slope_step = round(slope_step, 4), slope_graded = signif(slope_graded, 4),
  g_2024 = signif(g_2024, 4))])

## E9 and E9b slopes of the on-minus-off contrast, path 1, from the
## formula and directly from the captured matrix.
e9_formula <- function(g, on, variant) {
  n_off <- sum(!on)
  gbar_off <- mean(g[!on])
  if (variant == 'E9') -k_slope * (1 - gbar_off)
  else -k_slope * (1 - n_off / (n_off + 1) * gbar_off)
}
ex1_slopes <- rbindlist(lapply(c(0, 0.1, 1), function(th) {
  rbindlist(lapply(c('published', 'graded', 'upstream'), function(v) {
    cap <- capture(v, td1, CBM, th)
    g <- coupling_of(cap, CBM)
    on <- td1$tod > 0
    rbindlist(lapply(c('E9', 'E9b'), function(vr) {
      data.table(t_half = th, version = v, variant = vr,
        gbar_on = mean(g[on]), gbar_off = mean(g[!on]),
        gap = mean(g[on]) - mean(g[!on]),
        slope_formula = if (v == 'upstream') NA_real_ else e9_formula(g, on, vr),
        slope_direct = path_moments(cap, td1, vr)$gamma)
    }))
  }))
}))
fwrite(ex1_slopes, file.path(out_dir, 'example1-slopes.csv'))
cat('\nExample 1 contrast slopes (path 1)\n')
print(ex1_slopes, digits = 5)
stopifnot(with(ex1_slopes[version != 'upstream'],
  max(abs(slope_formula - slope_direct))) < 1e-12)

## ---- 3. Every path: does the graded patch reproduce the intended rule? -
halves <- c(0.1, 0.2, 0.5, 1, 2)
chk <- rbindlist(lapply(names(designs), function(dn) {
  rbindlist(lapply(seq_along(designs[[dn]]), function(k) {
    td <- designs[[dn]][[k]]
    rbindlist(lapply(halves, function(th) {
      cp <- capture('published', td, CBM, th)
      cg <- capture('graded', td, CBM, th)
      on <- td$tod > 0
      intended <- ifelse(on, 1, ifelse(td$tsd > 0, 2^(-td$tsd / th), 0))
      gs <- coupling_of(cp, CBM)
      gg <- coupling_of(cg, CBM)
      ## response block (everything but the bm row and column)
      same_block <- max(abs(cp$Sigma[-1, -1] - cg$Sigma[-1, -1]))
      data.table(design = dn, path = k, t_half = th,
        graded_vs_intended = max(abs(gg - intended)),
        step_is_mean_indicator = max(abs(as.numeric(br_mean_of(cp) != 0) -
          gs)) < 1e-12,
        step_gbar_off = if (any(!on)) mean(gs[!on]) else NA_real_,
        step_gbar_off_postexp = if (any(!on & td$tsd > 0))
          mean(gs[!on & td$tsd > 0]) else NA_real_,
        same_means = identical(cp$mu, cg$mu), max_block_diff = same_block)
    }))
  }))
}))
fwrite(chk, file.path(out_dir, 'check-graded-patch.csv'))
cat('\nGraded patch against the intended rule, 8 paths x 5 half-lives\n')
print(chk[, .(cases = .N, max_graded_vs_intended = max(graded_vs_intended),
  step_rule_is_mean_indicator = all(step_is_mean_indicator),
  step_off_coupling_after_exposure = paste(unique(na.omit(step_gbar_off_postexp)),
    collapse = ','),
  means_identical = all(same_means), max_response_block_diff = max(max_block_diff))])

## Package buildSigma(): same rule.
pkg_ok <- tryCatch({
  suppressMessages(pkgload::load_all('.', quiet = TRUE, export_all = TRUE))
  res <- rbindlist(lapply(names(designs), function(dn) {
    rbindlist(lapply(seq_along(designs[[dn]]), function(k) {
      td <- designs[[dn]][[k]]
      rbindlist(lapply(halves, function(th) {
        b <- buildSigma(mparams(CBM, th), extracted_rp, extracted_bp, td,
          makePositiveDefinite = FALSE)
        R <- cov2cor(b$sigma)
        g <- R[1, grep('\\.br$', b$labels)] / CBM
        on <- td$tod > 0
        intended <- ifelse(on, 1, ifelse(td$tsd > 0, 2^(-td$tsd / th), 0))
        data.table(design = dn, path = k, t_half = th,
          diff = max(abs(g - intended)))
      }))
    }))
  }))
  fwrite(res, file.path(out_dir, 'check-package-buildSigma.csv'))
  max(res$diff)
}, error = function(e) {
  cat('package check skipped:', conditionMessage(e), '\n')
  NA_real_
})
cat('package buildSigma() vs intended coupling, max abs diff:', pkg_ok, '\n')

## ---- 4. Closed-form power against the half-life -------------------------
th_tab <- c(0, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.0125, 0.02, 0.05, 0.1,
  0.2, 0.5, 1, 2, 4)
pw <- rbindlist(lapply(th_tab, function(th) {
  rbindlist(lapply(c('published', 'graded', 'upstream'), function(v) {
    rbindlist(lapply(c('E9', 'E9b'), function(vr)
      cbind(data.table(t_half = th, version = v, variant = vr),
        closed(v, 'Hybrid', CBM, th, vr))))
  }))
}))
fwrite(pw, file.path(out_dir, 'closed-form-power.csv'))
cat('\nClosed-form power, Hybrid, N = 70, c_bm = 0.25\n')
print(dcast(pw, t_half ~ version + variant, value.var = 'power'), digits = 3)
print(dcast(pw, t_half ~ version + variant, value.var = 'mean_slope'),
  digits = 3)

## Coupling gap by half-life (patient-weighted over paths).
gap_tab <- rbindlist(lapply(th_tab, function(th) {
  rbindlist(lapply(c('published', 'graded', 'upstream'), function(v) {
    tds <- designs$Hybrid
    Ns <- N %/% 4 + c(rep(1L, N %% 4), rep(0L, 4 - N %% 4))
    gp <- vapply(tds, function(td) {
      g <- coupling_of(capture(v, td, CBM, th), CBM)
      on <- td$tod > 0
      mean(g[on]) - mean(g[!on])
    }, numeric(1))
    data.table(t_half = th, version = v, gap = sum(Ns / N * gp))
  }))
}))
fwrite(gap_tab, file.path(out_dir, 'coupling-gap.csv'))
cat('\nCoupling gap (on minus off mean coupling), Hybrid\n')
print(dcast(gap_tab, t_half ~ version, value.var = 'gap'), digits = 4)

## ---- 5. Boundary: tiny half-lives and floating-point underflow ---------
## Restated lines 83-92 (docs/53 checked the restatement against the code).
pub_mean <- function(td, th) {
  rp <- extracted_rp[cat == 'br']
  m <- envs$published$modgompertz(td$tod, rp$max, rp$disp, rp$rate)
  for (p in 2:nrow(td)) {
    if (td$tod[p] == 0 && td$tsd[p] > 0) m[p] <- m[p] + m[p - 1] *
      (1/2)^(td$tsd[p] / th)
  }
  m
}
under <- rbindlist(lapply(names(designs), function(dn) {
  rbindlist(lapply(seq_along(designs[[dn]]), function(k) {
    td <- designs[[dn]][[k]]
    offs <- which(td$tod == 0 & td$tsd > 0)
    rbindlist(lapply(offs, function(p) {
      lo <- 1e-6; hi <- 1
      for (it in 1:200) {
        mid <- sqrt(lo * hi)
        if (pub_mean(td, mid)[p] == 0) lo <- mid else hi <- mid
      }
      data.table(design = dn, path = k, visit = td$timeptnames[p],
        t_sd = td$tsd[p], threshold = hi)
    }))
  }))
}))
## check the bisection against the code itself, just above and below
under[, code_zero_below := mapply(function(dn, k, v, h) {
  td <- designs[[dn]][[k]]
  br_mean_of(capture('published', td, CBM, h * 0.999))[td$timeptnames == v] == 0
}, design, path, visit, threshold)]
under[, code_nonzero_above := mapply(function(dn, k, v, h) {
  td <- designs[[dn]][[k]]
  br_mean_of(capture('published', td, CBM, h * 1.001))[td$timeptnames == v] != 0
}, design, path, visit, threshold)]
fwrite(under, file.path(out_dir, 'underflow-thresholds.csv'))
cat('\nHalf-life below which the published off-drug mean underflows to 0\n')
print(under, digits = 4)

bnd <- rbindlist(lapply(c(0.01, 0.005, 0.002, 0.001, 0.0005), function(th) {
  rbindlist(lapply(seq_along(designs$Hybrid), function(k) {
    td <- designs$Hybrid[[k]]
    g <- coupling_of(capture('published', td, CBM, th), CBM)
    on <- td$tod > 0
    data.table(t_half = th, path = k,
      step_off = paste(format(g[!on]), collapse = ' '),
      graded_off_max = max(coupling_of(capture('graded', td, CBM, th),
        CBM)[!on]))
  }))
}))
fwrite(bnd, file.path(out_dir, 'boundary-coupling.csv'))
cat('\nStep coupling at the off-drug visits, tiny half-lives\n')
print(bnd)
cat(sprintf('2^(-1/0.001) = %g; 2^(-7/0.001) = %g; 2^(-1/0.0005) = %g\n',
  2^(-1 / 0.001), 2^(-7 / 0.001), 2^(-1 / 0.0005)))

## ---- 6. Ceilings ---------------------------------------------------------
ceil <- rbindlist(lapply(names(designs), function(dn) {
  rbindlist(lapply(c(0, 0.1, 0.2, 0.5, 1, 2), function(th) {
    data.table(design = dn, t_half = th,
      step = ceiling_of('published', dn, th),
      graded = ceiling_of('graded', dn, th),
      upstream = ceiling_of('upstream', dn, th))
  }))
}))
fwrite(ceil, file.path(out_dir, 'ceilings.csv'))
cat('\nCeilings, rho = 0.8 (min over paths)\n')
print(ceil, digits = 4)
ceil7 <- rbindlist(lapply(c('Hybrid', 'OL+BDC', 'CO'), function(dn) {
  rbindlist(lapply(c(0, 0.5, 1), function(th)
    data.table(design = dn, t_half = th,
      upstream_rho07 = ceiling_of('upstream', dn, th, rho = 0.7))))
}))
cat('\n2024 construct ceilings at rho = 0.7 (cross-check of docs/54)\n')
print(ceil7, digits = 4)
fwrite(ceil7, file.path(out_dir, 'ceilings-2024-rho07.csv'))

## ---- 7. The 2024 revision: ratio coupling ------------------------------
up_tab <- rbindlist(lapply(list(c('Hybrid', 1), c('CO', 1)), function(x) {
  td <- designs[[x[1]]][[as.integer(x[2])]]
  rbindlist(lapply(c(2, 1), function(sf) {
    cu <- capture('upstream', td, CBM, 1, scalefactor = sf)
    ca <- capture('upstream_anchored', td, CBM, 1, scalefactor = sf)
    on <- td$tod > 0
    data.table(design = x[1], path = as.integer(x[2]), sf = sf,
      visit = td$timeptnames, D = as.integer(on), t_sd = td$tsd,
      g_2024 = coupling_of(cu, CBM), g_2024_anchored = coupling_of(ca, CBM),
      g_graded_stated = ifelse(on, 1, ifelse(td$tsd > 0, 2^(-td$tsd), 0)),
      g_graded_half = ifelse(on, 1, ifelse(td$tsd > 0, 2^(-sf * td$tsd), 0)))
  }))
}))
fwrite(up_tab, file.path(out_dir, 'upstream-2024-coupling.csv'))
cat('\n2024 ratio coupling at t_half = 1, compounding vs anchored means\n')
print(up_tab[, .(design, sf, visit, D, t_sd, g_2024 = signif(g_2024, 4),
  g_2024_anchored = signif(g_2024_anchored, 4),
  g_graded_stated = signif(g_graded_stated, 4))])
## identity check: 2024 coupling = graded at half-life t_half / sf, except
## visit 1 (never assigned), over all 8 paths and 5 half-lives.
up_chk <- rbindlist(lapply(names(designs), function(dn) {
  rbindlist(lapply(seq_along(designs[[dn]]), function(k) {
    td <- designs[[dn]][[k]]
    rbindlist(lapply(halves, function(th) {
      g <- coupling_of(capture('upstream', td, CBM, th), CBM)
      on <- td$tod > 0
      target <- ifelse(on, 1, ifelse(td$tsd > 0, 2^(-2 * td$tsd / th), 0))
      target[1] <- 0
      data.table(design = dn, path = k, t_half = th,
        diff = max(abs(g - target)))
    }))
  }))
}))
cat('2024 coupling vs graded at t_half/2 with g_1 = 0, 40 cases, max diff:',
  max(up_chk$diff), '\n')

## Example of the 58b32a9 code's own repair status at the article's 0.6.
pd_06 <- rbindlist(lapply(c(0, 0.1), function(th) {
  rbindlist(lapply(seq_along(designs$Hybrid), function(k) {
    cap <- capture('published', designs$Hybrid[[k]], 0.6, th)
    data.table(t_half = th, path = k,
      positive_definite = is.positive.definite(cap$Sigma))
  }))
}))
cat('\nHybrid, c_bm = 0.6: is the published matrix positive definite?\n')
print(dcast(pd_06, path ~ t_half, value.var = 'positive_definite'))

## ---- 8. Main drug effect: marginal distribution unchanged ---------------
marg <- rbindlist(lapply(names(designs), function(dn) {
  rbindlist(lapply(seq_along(designs[[dn]]), function(k) {
    td <- designs[[dn]][[k]]
    rbindlist(lapply(halves, function(th) {
      cp <- capture('published', td, CBM, th)
      cg <- capture('graded', td, CBM, th)
      data.table(design = dn, path = k, t_half = th,
        mean_diff = max(abs(cp$mu - cg$mu)),
        block_diff = max(abs(cp$Sigma[-1, -1] - cg$Sigma[-1, -1])),
        e9_mean_pub = path_moments(cp, td, 'E9')$mu,
        e9_mean_gr = path_moments(cg, td, 'E9')$mu)
    }))
  }))
}))
cat('\nMarginal (biomarker-free) distribution: max |difference| over 40 cases\n')
print(marg[, .(mean = max(mean_diff), response_block = max(block_diff),
  e9_contrast_mean = max(abs(e9_mean_pub - e9_mean_gr)))])

## ---- 9. Figures ----------------------------------------------------------
pal <- c('Graded (intended)' = '#2a78d6', 'Step (58b32a9)' = '#eb6834',
  '2024 upstream' = '#6b6b6b')
lty <- c('Graded (intended)' = 'solid', 'Step (58b32a9)' = 'dashed',
  '2024 upstream' = 'dotted')
shp <- c('Graded (intended)' = 16, 'Step (58b32a9)' = 17, '2024 upstream' = 15)
on_spans <- function(td) {
  wk <- c(0, envs$published$cumulative(td$t_wk))
  od <- td$tod > 0
  data.table(xmin = wk[-length(wk)][od], xmax = wk[-1][od])
}
f1 <- rbind(
  ex1[, .(t_half, week, g = g_step, rule = 'Step (58b32a9)')],
  ex1[, .(t_half, week, g = g_graded, rule = 'Graded (intended)')])
f1 <- rbind(f1, data.table(t_half = c(0.1, 1, 0.1, 1), week = 0, g = 0,
  rule = rep(c('Step (58b32a9)', 'Graded (intended)'), each = 2)))
f1[, panel := factor(sprintf('t[1/2] == %s ~ "week%s"', t_half,
  ifelse(t_half == 1, '', 's')), levels = c('t[1/2] == 0.1 ~ "weeks"',
  't[1/2] == 1 ~ "week"'))]
setorder(f1, rule, t_half, week)
p1 <- ggplot(f1, aes(week, g, colour = rule, linetype = rule, shape = rule)) +
  geom_rect(data = on_spans(td1), inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = 'grey92') +
  geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
  facet_wrap(~panel, labeller = label_parsed) +
  scale_colour_manual(values = pal[1:2]) +
  scale_linetype_manual(values = lty[1:2]) +
  scale_shape_manual(values = shp[1:2]) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  labs(x = 'Week (shaded: on drug)',
    y = expression('Coupling ' * g[t] * ' = Cor(bm, ' * BR[t] * ') / ' * c[bm]),
    colour = NULL, linetype = NULL, shape = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = 'bottom', panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, '55-fig1-coupling-by-visit.png'), p1, width = 7.5,
  height = 3.8, dpi = 200, bg = 'white')

## Figure 2: closed-form power (E9b, the proxy for the published analysis,
## and E9) on a log half-life axis, t_half = 0 drawn at the left edge,
## with lmer power from Hendrickson's code.
th_fine <- sort(unique(c(10^seq(log10(2e-4), log10(4), length.out = 140),
  0.0005, 0.001, 0.002, 0.01, 0.0125, 0.1, 0.2, 0.5, 1, 2)))
x0 <- 1e-4
cf <- rbindlist(lapply(c(0, th_fine), function(th) {
  rbindlist(lapply(c('published', 'graded', 'upstream'), function(v) {
    rbindlist(lapply(c('E9b', 'E9'), function(vr)
      cbind(data.table(t_half = th, version = v, variant = vr),
        closed(v, 'Hybrid', CBM, th, vr))))
  }))
}))
fwrite(cf, file.path(out_dir, 'closed-form-power-fine.csv'))
lab_v <- c(published = 'Step (58b32a9)', graded = 'Graded (intended)',
  upstream = '2024 upstream')
cf[, rule := factor(lab_v[version], levels = names(pal))]
cf[, x := ifelse(t_half == 0, x0, t_half)]
cf[, variant := factor(variant, levels = c('E9b', 'E9'),
  labels = c('E9b (proxy for the published analysis)', 'E9'))]
sims <- list()
if (file.exists(stored14)) {
  s14 <- fread(stored14)
  sims$s14 <- s14[, .(t_half, version, power, source = 'script 14')]
}
f55 <- file.path(sim_dir, 'power-check-reps500.csv')
if (file.exists(f55)) {
  s55 <- fread(f55)
  sims$s55 <- s55[, .(t_half, version, power, source = '55-power-check')]
}
sims <- rbindlist(sims)
if (nrow(sims)) {
  sims[, rule := factor(lab_v[version], levels = names(pal))]
  sims[, x := ifelse(t_half == 0, x0, t_half)]
  sims[, variant := factor('E9b (proxy for the published analysis)',
    levels = levels(cf$variant))]
}
cf_line <- cf[t_half > 0]
cf_zero <- cf[t_half == 0]
p2 <- ggplot(cf_line, aes(x, power, colour = rule, linetype = rule)) +
  geom_vline(xintercept = 1.6e-4, colour = 'grey70', linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(data = cf_zero, aes(shape = rule), size = 2.4) +
  facet_wrap(~variant) +
  scale_x_log10(breaks = c(x0, 1e-3, 1e-2, 0.1, 1),
    labels = c('0', '0.001', '0.01', '0.1', '1')) +
  scale_colour_manual(values = pal) + scale_linetype_manual(values = lty) +
  scale_shape_manual(values = shp, guide = 'none') +
  labs(x = 'Carryover half-life (weeks, log scale; 0 at left)',
    y = 'Power', colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = 'bottom', panel.grid.minor = element_blank())
if (nrow(sims)) {
  p2 <- p2 + geom_point(data = sims, aes(x, power, colour = rule),
    inherit.aes = FALSE, shape = 1, size = 2.6, stroke = 0.9)
}
ggsave(file.path(fig_dir, '55-fig2-power-vs-halflife.png'), p2, width = 8,
  height = 4, dpi = 200, bg = 'white')
cat('\nwrote figures and', out_dir, '\n')
