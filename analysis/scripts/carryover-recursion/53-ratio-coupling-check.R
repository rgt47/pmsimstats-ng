# Upstream 8609f12 ratio coupling c.bm * mu_t / mu_{t-1} under compounding
# vs anchored means, against graded coupling 2^(-t_sd / t_half).
# Crossover path 1 off-drug run: t_sd = 2.5, 5, 7.5, 10 weeks.
tsd <- c(2.5, 5, 7.5, 10)
mu_last <- 4.283
for (sf in c(1, 2)) {
  for (th in c(1, 2)) {
    comp <- numeric(4); prev <- mu_last
    for (k in 1:4) { comp[k] <- prev * 0.5^(sf * tsd[k] / th); prev <- comp[k] }
    anch <- mu_last * 0.5^(sf * tsd / th)
    ratio_comp <- comp / c(mu_last, comp[-4])
    ratio_anch <- anch / c(mu_last, anch[-4])
    graded_eff <- 0.5^(tsd / (th / sf))
    cat(sprintf('s_f=%d t_half=%g\n', sf, th))
    cat('  ratio rule, compounding means:', signif(ratio_comp, 4), '\n')
    cat('  graded at t_half/s_f         :', signif(graded_eff, 4), '\n')
    cat('  ratio rule, anchored means   :', signif(ratio_anch, 4), '\n')
    cat('  intended graded (t_half)     :', signif(0.5^(tsd / th), 4), '\n')
  }
}
