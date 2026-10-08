## Blinded-discontinuation (BD) blocks against crossover (CO) blocks,
## everything else equal; and the Hybrid design with its crossover block
## run as 4 weekly observations instead of 2 observations over 8 weeks.
## Closed form of docs/37 with the design-search machinery of
## 04-design-search.R (its definitions evaluated; nothing is simulated).
## Whitepaper: docs/39-nof1-design-optimization.md, Section 8.
##
## Every design is also simulated directly (stratified E9, both
## analyses of the open-label visits from the same simulated trials).
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/05-block-comparison.R \
##     [--mc-reps R] [--cores K]
suppressMessages({library(data.table); library(MASS); library(parallel)})
src04 <- readLines('analysis/scripts/quick-sim/carryover-closed-form/04-design-search.R')
stop_at <- grep('^## ---------------- design space', src04)[1] - 1
eval(parse(text = src04[1:stop_at]))
## The evaluated 04 header reads --mc-reps and --cores into mc_reps and
## cores.
out_dir <- 'analysis/data/quick-sim/carryover-closed-form/design-search'

## A design: visit times, expectancy, paths (on/off patterns), and the
## visits used in the on/off contrast. N = 70 split across paths as in
## the published code.
design <- function(times, e, paths, use = rep(TRUE, length(times))) {
  list(t_wk = diff(c(0, times)), e = e, paths = paths, use = use)
}
## Response covariance of one path (factor-major TV, PB, BR blocks).
resp_cov <- function(pi, structure) {
  nv <- length(pi$w)
  sds <- c(rep(10, nv), 10 * pi$e, rep(sd_BR, nv))
  R <- if (structure == 'CS') {
    rho <- dgp$RHO; c1 <- dgp$C1; cx <- dgp$CX
    Aw <- matrix(rho, nv, nv); diag(Aw) <- 1
    Cr <- matrix(cx, nv, nv); diag(Cr) <- c1
    rbind(cbind(Aw, Cr, Cr), cbind(Cr, Aw, Cr), cbind(Cr, Cr, Aw))
  } else {
    nu <- as.numeric(sub('nug', '', structure))
    kronecker(K3, (1 - nu) * dgp$RHO^abs(outer(pi$w, pi$w, '-')) + nu * diag(nv))
  }
  outer(sds, sds) * R
}
alloc <- function(P) N %/% P + c(rep(1L, N %% P), rep(0L, P - N %% P))
eval_general <- function(des, th, structure) {
  nv <- length(des$t_wk)
  n <<- nv
  P <- length(des$paths)
  Ns_d <- alloc(P)
  ps <- lapply(des$paths, function(od) dgp$path_info(des$t_wk, des$e, od))
  m <- rbindlist(lapply(ps, function(pi) {
    a <- weights_use(pi, des$use)
    h <- mod_profile(pi, 'mean_decayed', th)
    L <- cbind(diag(nv), diag(nv), diag(nv))
    aCa <- drop(t(a) %*% L %*% resp_cov(pi, structure) %*% t(L) %*% a)
    shift <- CBM * sd_BR * sum(a * h)
    data.table(gamma = -CBM * sd_BR / sd_B * sum(a * h), var_D = aCa + shift^2,
               tau2 = aCa, gbar_off = mean(h[!pi$on & des$use]),
               n_on = sum(pi$on & des$use), n_off = sum(!pi$on & des$use))
  }))
  r <- strat_power(m, Ns_d)
  data.table(power = r$power, info = r$info, var_mean = mean(m$var_D),
             h_off = mean(m$gbar_off))
}

ol <- c(4, 8)
blk <- c(9, 10, 11, 12)
e_ol <- c(1, 1)
e_blk <- rep(0.5, 4)
## Like for like: open-label lead-in, then one 4-week block of 4 weekly
## visits run as blinded discontinuation, crossover, or alternating.
blocks <- list(
  'BD block (stop at week 10 or 9)' = list(c(1, 1, 1, 1, 0, 0), c(1, 1, 1, 0, 0, 0)),
  'BD block (stop at week 11 or 9)' = list(c(1, 1, 1, 1, 1, 0), c(1, 1, 1, 0, 0, 0)),
  'CO block (on-on-off-off / off-off-on-on)' = list(c(1, 1, 1, 1, 0, 0), c(1, 1, 0, 0, 1, 1)),
  'Alternating block (on-off-on-off / off-on-off-on)' = list(c(1, 1, 1, 0, 1, 0),
                                                              c(1, 1, 0, 1, 0, 1)))
## The Hybrid design as published, and with its crossover block run as 4
## weekly observations (weeks 13-16) instead of 2 over 8 weeks (16, 20).
bd <- list(c(1, 1, 0, 0), c(1, 0, 0, 0))
co2 <- list(c(1, 0), c(0, 1))
co4 <- list(c(1, 1, 0, 0), c(0, 0, 1, 1))
mk_paths <- function(bdp, cop) {
  out <- list()
  for (b in bdp) for (cc in cop) out[[length(out) + 1]] <- c(1, 1, b, cc)
  out
}
## Equal visits and duration: BD block shortened to 2 visits (weeks 10,
## 12; discontinuation after week 10 or week 8, so one BD path is on at
## week 10 and off at 12, the other off at both) to make room for a
## 4-visit CO block every 2 weeks (14-20).
bd2 <- list(c(1, 0), c(0, 0))
hybrids <- list(
  'Hybrid as published (8 visits, 20 weeks)' =
    list(times = c(ol, blk, 16, 20), paths = mk_paths(bd, co2), e = c(e_ol, rep(0.5, 6))),
  'Hybrid, BD 2 visits + CO 4 visits every 2 weeks (8 visits, 20 weeks)' =
    list(times = c(ol, 10, 12, 14, 16, 18, 20), paths = mk_paths(bd2, co4),
         e = c(e_ol, rep(0.5, 6))),
  'Hybrid, CO block 4 weekly visits (10 visits, 16 weeks)' =
    list(times = c(ol, blk, 13, 14, 15, 16), paths = mk_paths(bd, co4), e = c(e_ol, rep(0.5, 8))),
  'Hybrid, CO block 4 visits every 2 weeks (10 visits, 20 weeks)' =
    list(times = c(ol, blk, 14, 16, 18, 20), paths = mk_paths(bd, co4), e = c(e_ol, rep(0.5, 8))))

## All designs as (set, name, times, e, paths).
all_des <- c(
  lapply(names(blocks), function(nm) list(set = 'block', name = nm, times = c(ol, blk),
                                          e = c(e_ol, e_blk), paths = blocks[[nm]])),
  lapply(names(hybrids), function(nm) c(list(set = 'hybrid', name = nm), hybrids[[nm]])))
rows <- list()
for (structure in c('CS', 'nug0.2')) for (th in c(0, 0.5, 1)) for (ol_in in c(TRUE, FALSE)) {
  for (dd in all_des) {
    des <- design(dd$times, dd$e, dd$paths,
                  c(rep(ol_in, 2), rep(TRUE, length(dd$times) - 2)))
    rows[[length(rows) + 1]] <- cbind(set = dd$set, name = dd$name, structure = structure,
                                      t_half = th, ol_in_contrast = ol_in,
                                      eval_general(des, th, structure))
  }
}
res <- rbindlist(rows)

## ---------------- Monte Carlo ----------------
## One simulated set of trials per design, structure and t1/2; the
## stratified E9 is computed with the open-label visits in and out of the
## contrast from the same data.
mc_general <- function(dd, th, structure, reps, seed) {
  set.seed(seed)
  nv <- length(dd$times)
  n <<- nv
  t_wk <- diff(c(0, dd$times))
  ps <- lapply(dd$paths, function(od) dgp$path_info(t_wk, dd$e, od))
  Ns_d <- alloc(length(ps))
  uses <- list(in_ = rep(TRUE, nv), out = c(FALSE, FALSE, rep(TRUE, nv - 2)))
  Bc <- NULL; Dc <- list(in_ = NULL, out = NULL)
  for (j in seq_along(ps)) {
    pi <- ps[[j]]
    X <- mvrnorm(Ns_d[j] * reps, rep(0, 3 * nv), resp_cov(pi, structure))
    z <- rnorm(Ns_d[j] * reps)
    br <- (2 * nv + 1):(3 * nv)
    X[, br] <- X[, br] + CBM * sd_BR * outer(z, mod_profile(pi, 'mean_decayed', th))
    S <- X[, 1:nv] + X[, (nv + 1):(2 * nv)] + X[, br]
    mt <- function(v) matrix(v, nrow = reps, byrow = TRUE)
    Bm <- mt(z * sd_B)
    Bc <- cbind(Bc, Bm - rowMeans(Bm))
    for (u in names(uses)) {
      Dm <- mt(-drop(S %*% weights_use(pi, uses[[u]])))
      Dc[[u]] <- cbind(Dc[[u]], Dm - rowMeans(Dm))
    }
  }
  df <- N - length(ps) - 1
  sxx <- rowSums(Bc^2)
  vapply(names(uses), function(u) {
    slope <- rowSums(Bc * Dc[[u]]) / sxx
    rv <- (rowSums(Dc[[u]]^2) - slope^2 * sxx) / df
    mean(2 * pt(-abs(slope / sqrt(rv / sxx)), df) < ALPHA)
  }, numeric(1))
}
mc_jobs <- CJ(k = seq_along(all_des), structure = c('CS', 'nug0.2'), t_half = c(0, 0.5, 1))
mc_file <- file.path(out_dir, sprintf('block-comparison-mc-reps%d.csv', mc_reps))
if (file.exists(mc_file)) {
  mc <- fread(mc_file)
} else {
  mc <- rbindlist(mclapply(seq_len(nrow(mc_jobs)), function(r) {
    j <- mc_jobs[r]; dd <- all_des[[j$k]]
    cell <- sprintf('%s_%s_t%.1f', dd$name, j$structure, j$t_half)
    seed <- sum(utf8ToInt(cell) * seq_along(utf8ToInt(cell))) %% 1e6 + 23
    p <- mc_general(dd, j$t_half, j$structure, mc_reps, seed)
    data.table(name = dd$name, structure = j$structure, t_half = j$t_half,
               ol_in_contrast = c(TRUE, FALSE), mc_power = unname(p))
  }, mc.cores = cores, mc.preschedule = FALSE))
  fwrite(mc, mc_file)
}
res <- merge(res, mc, by = c('name', 'structure', 't_half', 'ol_in_contrast'))
res[, z := (mc_power - power) / sqrt(power * (1 - power) / mc_reps)]
cat(sprintf('\nMonte Carlo check (%d reps): %d design-conditions, max |z| %.2f, mean z %.2f, beyond 2: %d\n',
            mc_reps, nrow(res), max(abs(res$z)), mean(res$z), sum(abs(res$z) > 2)))
fwrite(res, file.path(out_dir, 'block-comparison.csv'))
for (s in c('block', 'hybrid')) {
  cat(sprintf('\n%s comparison: power (closed form), c_bm = 0.25, N = 70, stratified E9\n', s))
  print(dcast(res[set == s, .(name, structure,
                              col = sprintf('t%.1f-%s', t_half, ifelse(ol_in_contrast, 'olin', 'olout')),
                              p = sprintf('%.3f (%.3f)', power, mc_power))],
              name + structure ~ col, value.var = 'p'))
}
cat('\nContrast variance and surviving off-drug moderation (OL in contrast, t1/2 = 1)\n')
print(res[ol_in_contrast == TRUE & t_half == 1,
          .(name, structure, var_mean = round(var_mean, 1), h_off = round(h_off, 3))])
