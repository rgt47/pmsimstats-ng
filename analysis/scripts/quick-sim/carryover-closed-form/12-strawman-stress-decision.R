## Applies the decision rule of docs/46 (Section 7) to the output of
## 11-strawman-stress-test.R. Interpretations of points the
## pre-registration leaves open are marked 'interpretation' and are
## reported in docs/47.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/12-strawman-stress-decision.R \
##     [--strawman S|DS] [--no-symmetry] [--with-extras]
## The defaults reproduce the rule as pre-registered. Amendment 1 to
## docs/46 (2026-10-07) uses --strawman DS --no-symmetry: the strawman
## selects its analysis half-life by AIC, and the symmetry screen, which
## tested a property that holds by construction, is dropped (docs/48,
## Section 4). Outputs carry a suffix naming the amendment.
## --with-extras adds the three analyses fitted to the same trials after
## the pre-registration (G1, G2 and F1A from the --extra and --extra2
## runs of 11-strawman-stress-test.R) to every screen and to the
## selection.
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
straw <- if ('--strawman' %in% args) args[which(args == '--strawman') + 1] else 'S'
use_sym <- !('--no-symmetry' %in% args)
extras <- '--with-extras' %in% args
suffix <- paste0(if (straw == 'S' && use_sym) '' else
  paste0('-', tolower(straw), if (use_sym) '' else '-nosym'), if (extras) '-all14' else '')
d <- 'analysis/data/quick-sim/carryover-closed-form/strawman-stress'
x <- readRDS(file.path(d, 'replicates-full.rds'))
an <- c('S', 'A', 'C1', 'C2', 'DS', 'DA', 'ES', 'EA', 'F1', 'F2', 'R')
if (extras) {
  ex <- rbindlist(lapply(c('replicates-full-extra.rds', 'replicates-full-extra2.rds'),
                         function(f) readRDS(file.path(d, f))[analysis != 'S']), fill = TRUE)
  x <- rbindlist(list(x, ex), fill = TRUE)
  an <- c(an, 'G1', 'G2', 'F1A')
}
## Failed fits count as non-rejections (docs/46, Section 6).
x[, rej := !is.na(p) & p < 0.05]
x[, p1 := fifelse(is.na(p), 1, p)]
## The prognostic-baseline cells (40, 41) were added after the
## pre-registration; they are judged with the other prognostic cells.
prog <- c('prog-TV-const', 'prog-TV-grow', 'prog-PB', 'prog-BL')

## ---- 1. nominal screen: non-prognostic null cells ----
nul <- x[hyp == 'null', .(reject = mean(rej), n = .N, fails = sum(is.na(p))),
         by = .(label, analysis)]
main <- nul[!label %in% prog]
pooled <- main[, .(pooled = sum(reject * n) / sum(n), N = sum(n), max_cell = max(reject),
                   worst = label[which.max(reject)]), by = analysis]
pooled[, band := 2 * sqrt(0.05 * 0.95 / N)]
pooled[, nominal := abs(pooled - 0.05) <= band & max_cell <= 0.075]

## ---- 2. prognostic screen ----
pg <- dcast(nul[label %in% prog], analysis ~ label, value.var = 'reject')
pg[, prog_flag := Reduce(`|`, lapply(intersect(prog, names(pg)), function(v) get(v) > 0.075))]
screen <- merge(pooled, pg, by = 'analysis')
screen <- screen[match(an, analysis)]
cat('\nScreens (pooled over', main[, uniqueN(label)], 'non-prognostic null cells)\n')
print(screen, digits = 3)

## ---- size-adjusted power ----
## The matching null cell: same label for stress cells; the zero cell of
## the same curve for curve cells.
x[, null_label := label]
crit <- x[hyp == 'null', .(crit = quantile(p1, 0.05)), by = .(null_label = label, analysis)]
alt <- merge(x[hyp == 'alt'], crit, by = c('null_label', 'analysis'))
alt[, rej_adj := p1 < crit]
pw <- alt[, .(power = mean(rej), power_adj = mean(rej_adj), n = .N), by = .(label, cbm, analysis)]
pw[, mcse := sqrt(power_adj * (1 - power_adj) / n)]

## ---- 3. sensitivity screen: monotone and symmetric power curves ----
## interpretation: a 'decrease larger than 2 MCSE' uses the MCSE of the
## difference between two independent cells; the zero point of each
## curve is the size-adjusted rejection at the null quantile, 0.05.
curve <- rbind(pw[label %in% c('curve-C07', 'curve-CS')],
               unique(x[hyp == 'null' & label %in% c('curve-C07', 'curve-CS'),
                        .(label, cbm, analysis, power = mean(rej), power_adj = 0.05,
                          n = .N), by = .(label, analysis)][, .(label, cbm, analysis, power,
                                                               power_adj, n)])[
                 , mcse := sqrt(power_adj * (1 - power_adj) / n)], fill = TRUE)
sens <- curve[, {
  pos <- .SD[cbm >= 0][order(cbm)]
  dd <- diff(pos$power_adj); sd_d <- sqrt(head(pos$mcse, -1)^2 + tail(pos$mcse, -1)^2)
  mono <- all(dd >= -2 * sd_d)
  neg <- .SD[cbm < 0]; mirror <- pos[abs(cbm - abs(neg$cbm)) < 1e-9]
  sym <- abs(neg$power_adj - mirror$power_adj) <= 2 * sqrt(neg$mcse^2 + mirror$mcse^2)
  s15 <- pos[abs(cbm - 0.15) < 1e-9, power_adj]; s30 <- pos[abs(cbm - 0.30) < 1e-9, power_adj]
  steep <- if (length(s15) && length(s30)) (s30 - s15) / 0.15 else NA_real_
  k <- which(pos$power_adj >= 0.8)[1]
  detect <- if (is.na(k)) NA_real_ else if (k == 1) pos$cbm[1] else
    pos$cbm[k - 1] + (0.8 - pos$power_adj[k - 1]) * (pos$cbm[k] - pos$cbm[k - 1]) /
      (pos$power_adj[k] - pos$power_adj[k - 1])
  .(monotone = mono, symmetric = sym, power_neg = neg$power_adj, power_pos = mirror$power_adj,
    steepness = steep, detectable = detect)
}, by = .(label, analysis)]
cat('\nSensitivity (size-adjusted power curves)\n')
print(dcast(sens, analysis ~ label, value.var = c('monotone', 'symmetric', 'detectable',
                                                  'steepness')), digits = 3)
sens_flag <- if (use_sym) sens[, .(sens_ok = all(monotone & symmetric)), by = analysis] else
  sens[, .(sens_ok = all(monotone)), by = analysis]
cat(sprintf('\nStrawman: %s; symmetry screen: %s\n', straw, if (use_sym) 'on' else 'off'))

## ---- 4. maximin shortfall among survivors ----
## interpretation: the maximin set is every alternative cell except the
## three prognostic cells (judged by the prognostic screen); a version
## including them is also reported.
surv <- intersect(screen[nominal == TRUE, analysis], sens_flag[sens_ok == TRUE, analysis])
maximin <- function(cells_in) {
  p <- pw[analysis %in% surv & label %in% cells_in]
  p[, best := max(power_adj), by = .(label, cbm)]
  p[, shortfall := best - power_adj]
  p[, .(maximin = max(shortfall), worst_cell = paste(label, cbm)[which.max(shortfall)],
        mean_power_adj = mean(power_adj)), by = analysis][order(maximin)]
}
alt_main <- setdiff(unique(pw$label), prog)
mm_main <- maximin(alt_main)
mm_all <- maximin(unique(pw$label))
cat('\nSurvivors:', paste(surv, collapse = ', '), '\n')
cat('\nMaximin shortfall, non-prognostic alternative cells\n'); print(mm_main, digits = 3)
cat('\nMaximin shortfall, all alternative cells\n'); print(mm_all, digits = 3)

## ---- 5. strawman rule: McNemar on paired size-adjusted rejections ----
mcn <- function(a, b) {
  w <- dcast(alt[analysis %in% c(a, b) & label %in% alt_main], label + cbm + rep ~ analysis,
             value.var = 'rej_adj')
  w[, .(n10 = sum(get(a) & !get(b)), n01 = sum(!get(a) & get(b))), by = .(label, cbm)][
    , p_mcnemar := mapply(function(u, v) if (u + v == 0) 1 else
        binom.test(u, u + v)$p.value, n10, n01)][]
}
if (straw %in% surv) {
  sS <- mm_main[analysis == straw, maximin]
  for (a in setdiff(surv, straw)) {
    m <- mcn(a, straw)
    beats <- sS - mm_main[analysis == a, maximin]
    cat(sprintf('\n%s vs %s: maximin advantage %.3f; McNemar p < 0.05 in %d of %d cells (%s ahead in %d)\n',
                a, straw, beats, sum(m$p_mcnemar < 0.05), nrow(m), a,
                sum(m$p_mcnemar < 0.05 & m$n10 > m$n01)))
  }
} else cat(sprintf('\nStrawman %s did not survive the screens.\n', straw))
cat('\nSize-adjusted power, non-prognostic alternative cells\n')
print(dcast(pw[label %in% alt_main], label + cbm ~ analysis, value.var = 'power_adj'),
      digits = 3, nrows = 100)
cat('\nSize-adjusted power, prognostic alternative cells\n')
print(dcast(pw[label %in% prog], label + cbm ~ analysis, value.var = 'power_adj'), digits = 3)
cat('\nAIC selections (DS, DA), share by half-life, all cells\n')
print(dcast(x[analysis %in% c('DS', 'DA'), .N, by = .(analysis, h)][, share := N / sum(N),
            by = analysis], analysis ~ h, value.var = 'share'), digits = 2)
fwrite(screen, file.path(d, sprintf('decision-screen%s.csv', suffix)))
fwrite(sens, file.path(d, sprintf('decision-sensitivity%s.csv', suffix)))
fwrite(pw, file.path(d, sprintf('decision-power%s.csv', suffix)))
fwrite(rbind(mm_main[, set := 'non-prognostic'], mm_all[, set := 'all']),
       file.path(d, sprintf('decision-maximin%s.csv', suffix)))
