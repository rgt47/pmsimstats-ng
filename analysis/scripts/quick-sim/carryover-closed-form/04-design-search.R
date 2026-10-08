## Exhaustive search for the most powerful N-of-1 design for the
## biomarker-treatment interaction, using the closed form of docs/37.
## Whitepaper: docs/39-nof1-design-optimization.md.
##
## Design space: two counterbalanced sequences (an on/off pattern over 8
## blinded visits and its complement, 35 patients each), 2 to 6 on-drug
## visits per sequence, equal spacing d. A pattern and its complement are
## the same design, so 119 designs. Spacing d in {0.5, ..., 4} weeks.
## Rule: decayed mean moderation (anchored decay, no c_bm ceiling).
## Response correlation structures:
##   CS      the published compound-symmetry construct (configuration A)
##   nug0.2  separable AR(1), rho = 0.8 per week, with 20% of each
##           factor's variance white measurement error
##   nug0.5  the same with 50%
## Analysis: E9 with sequence-specific intercepts (stratified), so that
## sequence and period differences do not enter the residual.
## Criterion: per-patient information, the patient-weighted squared
## correlation between the contrast and the biomarker; power at N = 70,
## c_bm = 0.25. Response half-life t1/2 in {0.5, 1, 2} weeks.
## Checks: direct Monte Carlo of the stratified E9 for the best design in
## each condition and for reference designs.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/04-design-search.R \
##     [--mc-reps R] [--cores K]
suppressMessages({library(data.table); library(MASS); library(ggplot2); library(parallel)})
src <- readLines('analysis/scripts/quick-sim/carryover-closed-form/01-e9-closed-form.R')
stop_at <- grep('^## ---------------- per-path table', src)[1] - 1
eval(parse(text = src[1:stop_at]))
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, dflt) if (k %in% args) as.integer(args[which(args == k) + 1]) else dflt
mc_reps <- arg('--mc-reps', 5000L)
cores <- arg('--cores', 6L)
out_dir <- 'analysis/data/quick-sim/carryover-closed-form/design-search'
fig_dir <- 'docs/figures'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

STRUCTS <- c('CS', 'nug0.2', 'nug0.5')
HALVES <- c(0.5, 1, 2)
D_GRID <- c(0.5, 1, 1.5, 2, 2.5, 3, 4)
K3 <- (1 - dgp$C1) * diag(3) + dgp$C1

## ---------------- moments of one sequence ----------------
## Response correlation of a sequence, factor-major (TV, PB, BR blocks).
resp_cor <- function(pi, structure) {
  if (structure == 'CS') return(dgp$build(pi, 'A', 0, 0)$R[-(1:2), -(1:2)])
  nu <- as.numeric(sub('nug', '', structure))
  kronecker(K3, (1 - nu) * dgp$RHO^abs(outer(pi$w, pi$w, '-')) + nu * diag(n))
}
## Contrast weights over the visits used in the contrast (use = TRUE):
## +1/n_on on drug, -1/n_off off drug, 0 for unused visits (an
## open-label lead-in that the analysis leaves out).
weights_use <- function(pi, use = rep(TRUE, n)) {
  on <- pi$on & use; off <- !pi$on & use
  ifelse(on, 1 / sum(on), ifelse(off, -1 / sum(off), 0))
}
seq_moments <- function(pi, th, structure, use = rep(TRUE, n)) {
  a <- weights_use(pi, use)
  h <- mod_profile(pi, 'mean_decayed', th)
  sds <- c(rep(10, n), 10 * pi$e, rep(sd_BR, n))
  C <- outer(sds, sds) * resp_cor(pi, structure)
  L <- cbind(diag(n), diag(n), diag(n))
  aCa <- drop(t(a) %*% L %*% C %*% t(L) %*% a)
  shift <- CBM * sd_BR * sum(a * h)
  data.table(gamma = -CBM * sd_BR / sd_B * sum(a * h), var_D = aCa + shift^2,
             tau2 = aCa)
}

## Stratified pooled slope over P sequences: within-sequence sums of
## squares S_p ~ sd_B^2 chi^2_{n_p - 1}, so the pooled slope has
## N - P - 1 residual df and S ~ sd_B^2 chi^2_{N - P}. Slope differences
## between sequences act as extra residual variance (first order).
strat_power <- function(m, Ns) {
  P <- length(Ns); Nt <- sum(Ns)
  w <- (Ns - 1) / sum(Ns - 1)
  gbar <- sum(w * m$gamma)
  resid <- sum(w * (m$tau2 + sd_B^2 * (m$gamma - gbar)^2))
  df <- Nt - P - 1
  q <- qchisq((seq_len(400) - 0.5) / 400, Nt - P)
  ncp <- gbar * sd_B * sqrt(q) / sqrt(resid)
  tc <- qt(1 - ALPHA / 2, df)
  list(power = mean(pt(-tc, df, ncp = ncp) + pt(tc, df, ncp = ncp, lower.tail = FALSE)),
       info = sum(Ns / Nt * (m$gamma * sd_B)^2 / m$var_D), gbar = gbar)
}

## A design: k open-label on-drug visits for everyone (expectancy 1),
## then a blinded pattern over the remaining visits for sequence 1 and
## its complement for sequence 2 (expectancy 0.5).
make_paths <- function(pattern, d, k = 0L) {
  t_wk <- rep(d, n)
  e <- c(rep(1, k), rep(0.5, n - k))
  p1 <- c(rep(1, k), pattern); p2 <- c(rep(1, k), 1 - pattern)
  list(dgp$path_info(t_wk, e, p1), dgp$path_info(t_wk, e, p2))
}
eval_design <- function(pattern, d, th, structure, k = 0L, use_ol = TRUE) {
  ps <- make_paths(pattern, d, k)
  use <- c(rep(use_ol, k), rep(TRUE, n - k))
  m <- rbindlist(lapply(ps, seq_moments, th = th, structure = structure, use = use))
  strat_power(m, c(35L, 35L))
}

## ---------------- design space ----------------
pats <- as.matrix(expand.grid(rep(list(0:1), n)))
pats <- pats[rowSums(pats) >= 2 & rowSums(pats) <= 6, , drop = FALSE]
## A pattern and its complement are the same design: keep the one whose
## first visit is on drug.
pats <- pats[pats[, 1] == 1, , drop = FALSE]
pat_str <- apply(pats, 1, paste, collapse = '')
cat(sprintf('%d designs x %d spacings x %d structures x %d half-lives\n',
            nrow(pats), length(D_GRID), length(STRUCTS), length(HALVES)))
switches <- function(p) sum(diff(p) != 0)

jobs <- CJ(i = seq_len(nrow(pats)), d = D_GRID, structure = STRUCTS, t_half = HALVES)
res <- rbindlist(mclapply(split(jobs, by = 'i'), function(jb) {
  p <- pats[jb$i[1], ]
  rbindlist(lapply(seq_len(nrow(jb)), function(k) {
    r <- eval_design(p, jb$d[k], jb$t_half[k], jb$structure[k])
    data.table(pattern = pat_str[jb$i[k]], n_on = sum(p), switches = switches(p),
               d = jb$d[k], structure = jb$structure[k], t_half = jb$t_half[k],
               info = r$info, power = r$power)
  }))
}, mc.cores = cores))
fwrite(res, file.path(out_dir, 'design-search-all.csv'))

## ---------------- reference designs ----------------
refs <- data.table(
  name = c('CO', 'ABAB', 'Alternating'),
  pattern = c('11110000', '11001100', '10101010'))
ref_res <- res[pattern %in% refs$pattern][refs, on = 'pattern']

## ---------------- results ----------------
best <- res[, .SD[which.max(info)], by = .(structure, t_half)]
fwrite(best, file.path(out_dir, 'design-search-best.csv'))
cat('\nBest design (pattern of the first sequence; the second is its complement)\n')
print(best[, .(structure, t_half, pattern, n_on, switches, d, info = round(info, 4),
               power = round(power, 3))])
cat('\nReference designs at their best spacing\n')
print(ref_res[, .SD[which.max(info)], by = .(name, structure, t_half)][
  , .(name, structure, t_half, d, power = round(power, 3))][order(structure, t_half, name)])
## Best design at each spacing, to show where the optimum lies.
best_d <- res[, .SD[which.max(info)], by = .(structure, t_half, d)]
fwrite(best_d, file.path(out_dir, 'design-search-best-by-spacing.csv'))
cat('\nBest power at each spacing\n')
print(dcast(best_d[, .(structure, t_half, d, p = sprintf('%.3f %s', power, pattern))],
            structure + t_half ~ d, value.var = 'p'))
## Robustness: the design-spacing pair maximizing the minimum, over
## structures and half-lives, of power relative to the condition's best.
rel <- merge(res, best[, .(structure, t_half, best_power = power)],
             by = c('structure', 't_half'))
rel[, ratio := power / best_power]
robust <- rel[, .(min_ratio = min(ratio), mean_ratio = mean(ratio)), by = .(pattern, d)][
  order(-min_ratio)]
fwrite(robust, file.path(out_dir, 'design-search-robust.csv'))
cat('\nMost robust design-spacing pairs (minimum relative power over the 9 conditions)\n')
print(head(robust, 8))
## Characteristics of good designs: the best power by number of switches.
by_sw <- res[, .(best_power = max(power)), by = .(structure, t_half, switches)]
fwrite(by_sw, file.path(out_dir, 'design-search-by-switches.csv'))
cat('\nBest power by number of switches (first sequence)\n')
print(dcast(by_sw[, .(structure, t_half, switches, p = round(best_power, 3))],
            structure + t_half ~ switches, value.var = 'p'))

## ---------------- open-label lead-in ----------------
## k = 0-4 open-label on-drug visits (expectancy 1) for everyone, then
## counterbalanced blinded sequences over the remaining 8 - k visits
## (at least one on and one off visit each). Two analyses: the open-label
## visits included in the on/off contrast (as E9 and the published
## analysis do) or left out of it.
ol_rows <- list()
for (k in 0:4) {
  nb <- n - k
  bp <- as.matrix(expand.grid(rep(list(0:1), nb)))
  bp <- bp[rowSums(bp) >= 1 & rowSums(bp) <= nb - 1 & bp[, 1] == 1, , drop = FALSE]
  for (use_ol in c(TRUE, FALSE)) {
    if (k == 0 && !use_ol) next
    ol_jobs <- CJ(i = seq_len(nrow(bp)), d = D_GRID, structure = STRUCTS, t_half = HALVES)
    ol_rows[[length(ol_rows) + 1]] <- rbindlist(mclapply(split(ol_jobs, by = 'i'), function(jb) {
      p <- bp[jb$i[1], ]
      rbindlist(lapply(seq_len(nrow(jb)), function(r) {
        v <- eval_design(p, jb$d[r], jb$t_half[r], jb$structure[r], k, use_ol)
        data.table(lead_in = k, ol_in_contrast = use_ol, pattern = paste(p, collapse = ''),
                   d = jb$d[r], structure = jb$structure[r], t_half = jb$t_half[r],
                   info = v$info, power = v$power)
      }))
    }, mc.cores = cores))
  }
}
ol <- rbindlist(ol_rows)
fwrite(ol, file.path(out_dir, 'design-search-open-label.csv'))
ol_best <- ol[, .SD[which.max(info)], by = .(lead_in, ol_in_contrast, structure, t_half)]
fwrite(ol_best, file.path(out_dir, 'design-search-open-label-best.csv'))
cat('\nOpen-label lead-in: best power by number of open-label visits',
    '(in = open-label visits in the contrast; out = left out)\n')
print(dcast(ol_best[, .(structure, t_half, col = sprintf('k%d-%s', lead_in,
                                                          ifelse(ol_in_contrast, 'in', 'out')),
                        p = round(power, 3))],
            structure + t_half ~ col, value.var = 'p'))

## ---------------- Monte Carlo checks ----------------
mc_one <- function(pattern, d, th, structure, reps, seed, k = 0L, use_ol = TRUE) {
  set.seed(seed)
  p <- as.integer(strsplit(pattern, '')[[1]])
  ps <- make_paths(p, d, k)
  use <- c(rep(use_ol, k), rep(TRUE, n - k))
  Bs <- Ds <- NULL
  for (j in seq_along(ps)) {
    pi <- ps[[j]]
    sds <- c(rep(10, n), 10 * pi$e, rep(sd_BR, n))
    C <- outer(sds, sds) * resp_cor(pi, structure)
    X <- mvrnorm(35L * reps, rep(0, 3 * n), C)
    z <- rnorm(35L * reps)
    h <- mod_profile(pi, 'mean_decayed', th)
    br <- (2 * n + 1):(3 * n)
    X[, br] <- X[, br] + CBM * sd_BR * outer(z, h)
    S <- X[, 1:n] + X[, (n + 1):(2 * n)] + X[, br]
    D <- -drop(S %*% weights_use(pi, use))
    mt <- function(v) matrix(v, nrow = reps, byrow = TRUE)
    ## Centered within the sequence: the stratified analysis.
    Bm <- mt(z * sd_B); Dm <- mt(D)
    Bs <- cbind(Bs, Bm - rowMeans(Bm)); Ds <- cbind(Ds, Dm - rowMeans(Dm))
  }
  sxx <- rowSums(Bs^2); slope <- rowSums(Bs * Ds) / sxx
  df <- 70 - 2 - 1
  rv <- (rowSums(Ds^2) - slope^2 * sxx) / df
  pval <- 2 * pt(-abs(slope / sqrt(rv / sxx)), df)
  mean(pval < ALPHA)
}
chk <- unique(rbind(best[, .(pattern, d, structure, t_half)],
                    ref_res[d == 2.5, .(pattern, d, structure, t_half)]))
chk <- merge(chk, res[, .(pattern, d, structure, t_half, power)],
             by = c('pattern', 'd', 'structure', 't_half'))
chk[, `:=`(lead_in = 0L, ol_in_contrast = TRUE)]
## Plus the best two-visit lead-in designs, both analyses, at t1/2 = 1.
chk <- rbindlist(list(chk, ol_best[lead_in == 2 & t_half == 1,
                                   .(pattern, d, structure, t_half, power, lead_in,
                                     ol_in_contrast)]), use.names = TRUE)
## 'cell' rather than 'key': data.table() reads an argument named key as
## a sort key.
chk[, cell := sprintf('k%d%s_%s_d%.1f_%s_t%.1f', lead_in, ifelse(ol_in_contrast, 'in', 'out'),
                      pattern, d, structure, t_half)]
mc_file <- file.path(out_dir, sprintf('design-search-mc-reps%d.csv', mc_reps))
done <- if (file.exists(mc_file)) fread(mc_file) else
  data.table(cell = character(), mc_power = numeric())
todo <- chk[!cell %in% done$cell]
if (nrow(todo)) {
  new <- mclapply(seq_len(nrow(todo)), function(k) {
    j <- todo[k]
    seed <- sum(utf8ToInt(j$cell) * seq_along(utf8ToInt(j$cell))) %% 1e6 + 11
    tryCatch(data.table(cell = j$cell,
                        mc_power = mc_one(j$pattern, j$d, j$t_half, j$structure, mc_reps,
                                          seed, j$lead_in, j$ol_in_contrast)),
             error = function(e) paste(j$cell, conditionMessage(e)))
  }, mc.cores = cores, mc.preschedule = FALSE)
  errs <- Filter(is.character, new)
  if (length(errs)) stop('Monte Carlo failed: ', paste(head(unlist(errs), 3), collapse = '; '))
  new <- rbindlist(new)
  done <- rbind(done, new)
  fwrite(done, mc_file)
}
chk <- merge(chk, done, by = 'cell')
chk[, z := (mc_power - power) / sqrt(power * (1 - power) / mc_reps)]
fwrite(chk, file.path(out_dir, 'design-search-mc-check.csv'))
cat(sprintf('\nMonte Carlo check (%d reps): %d designs, max |z| %.2f, mean z %.2f, beyond 2: %d\n',
            mc_reps, nrow(chk), max(abs(chk$z)), mean(chk$z), sum(abs(chk$z) > 2)))
print(chk[, .(lead_in, ol_in_contrast, pattern, d, structure, t_half,
              power = round(power, 3), mc_power = round(mc_power, 3),
              z = round(z, 2))][order(structure, t_half)])

## ---------------- figures ----------------
ink <- '#0b0b0b'; ink2 <- '#52514e'; surface <- '#fcfcfb'
theme_fig <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = surface, colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = '#e6e5e1', linewidth = 0.3),
          axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink2),
          strip.text = element_text(colour = ink, face = 'bold', hjust = 0),
          plot.title = element_text(colour = ink, face = 'bold'),
          plot.subtitle = element_text(colour = ink2),
          plot.caption = element_text(colour = ink2, hjust = 0),
          plot.title.position = 'plot', plot.caption.position = 'plot',
          legend.position = 'top', legend.justification = 'left',
          legend.title = element_blank(), legend.text = element_text(colour = ink))
}
s_lab <- c(CS = 'Compound symmetry', nug0.2 = 'AR(1) + 20% measurement error',
           nug0.5 = 'AR(1) + 50% measurement error')
## Figure 2: every design at its best spacing, against switches.
at_best <- res[, .SD[which.max(info)], by = .(pattern, structure, t_half)]
at_best[, structure := factor(s_lab[structure], levels = s_lab)]
at_best[, th := factor(sprintf('t1/2 = %s wk', format(t_half)),
                       levels = sprintf('t1/2 = %s wk', format(HALVES)))]
bb <- copy(best)[, structure := factor(s_lab[structure], levels = s_lab)]
bb[, th := factor(sprintf('t1/2 = %s wk', format(t_half)), levels = levels(at_best$th))]
rr <- ref_res[, .SD[which.max(info)], by = .(name, structure, t_half)]
rr[, structure := factor(s_lab[structure], levels = s_lab)]
rr[, th := factor(sprintf('t1/2 = %s wk', format(t_half)), levels = levels(at_best$th))]
g2 <- ggplot(at_best, aes(switches, power)) +
  geom_jitter(width = 0.15, height = 0, colour = '#9ec5f0', size = 1.4) +
  geom_point(data = rr, aes(shape = name), colour = '#17447e', size = 2.6, stroke = 1) +
  geom_point(data = bb, colour = '#eb6834', size = 3.2) +
  scale_shape_manual(values = c(Alternating = 1, ABAB = 2, CO = 0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = 1:7) +
  facet_grid(structure ~ th) +
  labs(x = 'Switches between on and off drug (first sequence)', y = 'Power at the best spacing',
       title = 'All 119 counterbalanced designs at their best spacing',
       subtitle = 'Decayed mean moderation, N = 70, c_bm = 0.25, stratified E9, closed form',
       caption = paste('Light points: designs. Open dark symbols: reference designs.',
                       'Orange: the best design in each condition.')) +
  theme_fig()
ggsave(file.path(fig_dir, '39-fig2-design-landscape.png'), g2, width = 10, height = 8.6,
       dpi = 200, bg = surface)
cat('\nWrote', file.path(fig_dir, '39-fig2-design-landscape.png'), '\n')
