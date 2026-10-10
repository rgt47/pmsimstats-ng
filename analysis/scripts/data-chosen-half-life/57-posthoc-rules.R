## docs/57, Section 6: post hoc rules computed from the stored grid fits
## of 57-null-illustration.R (no refitting). Defined after the
## illustration's results were seen, and labeled post hoc in docs/57.
##   F1G  guarded F1A0: keep t1/2a = 0.5 unless the null-model profile
##        rejects it (2 * (max ll0 - ll0(0.5)) > qchisq(0.95, 1)); then
##        use the F1A0 choice.
##   MA0  model averaging with weights from the null-model likelihood
##        (ll0) instead of ll1.
## Also: paired discordance of F1A and F1A0 against F1, and the share
## of trials whose null-model profile rejects 0.5 week.
## Run from the repository root after 57-null-illustration.R.
suppressMessages(library(data.table))
out <- 'analysis/data/quick-sim/data-chosen-half-life'
x <- readRDS(file.path(out, '57-mc-replicates.rds'))
G <- x$grid; M <- x$methods
G[, dev05 := 2 * (max(ll0, na.rm = TRUE) - ll0[h == 0.5]), by = .(cell, rep)]
g_rule <- G[, {
  k <- if (isTRUE(dev05[1] > qchisq(0.95, 1))) which.max(ll0) else which(h == 0.5)
  .(method = 'F1G', h = h[k], p = p[k], switched = isTRUE(dev05[1] > qchisq(0.95, 1)))
}, by = .(cell, rep)]
ma0 <- G[, {
  ok <- is.finite(ll0) & is.finite(est); ll0 <- ll0[ok]; est <- est[ok]; se <- se[ok]; df <- df[ok]; h <- h[ok]
  w <- exp(ll0 - max(ll0)); w <- w / sum(w)
  eb <- sum(w * est); sb <- sum(w * sqrt(se^2 + (est - eb)^2)); dfb <- sum(w * df)
  .(method = 'MA0', h = sum(w * h), p = 2 * pt(-abs(eb / sb), dfb), switched = NA)
}, by = .(cell, rep)]
R <- rbind(g_rule, ma0)
summ <- R[, .(reject = mean(p < 0.05), R = .N, switched = mean(switched)), by = .(cell, method)]
summ[, mcse := sqrt(reject * (1 - reject) / R)]
print(summ)
fwrite(summ, file.path(out, '57-mc-posthoc-rules.csv'))
pooled <- rbind(R[cell %in% c(2, 20), .(reject = mean(p < .05), .N), by = method],
                M[cell %in% c(2, 20), .(reject = mean(p < .05, na.rm = TRUE), .N), by = method])
pooled[, mcse := sqrt(reject * (1 - reject) / N)]
print(pooled)
fwrite(pooled, file.path(out, '57-mc-pooled-null.csv'))
P <- dcast(rbind(M[, .(cell, rep, method, p)], R[, .(cell, rep, method, p)]),
           cell + rep ~ method, value.var = 'p')
disc <- rbindlist(lapply(c('F1A', 'F1A0', 'F1G', 'MA', 'MA0'), function(m)
  P[cell %in% c(2, 20, 25), .(method = m, only_m = sum(get(m) < .05 & F1 >= .05),
                               only_F1 = sum(F1 < .05 & get(m) >= .05)), by = cell]))
print(disc)
fwrite(disc, file.path(out, '57-mc-discordance.csv'))
