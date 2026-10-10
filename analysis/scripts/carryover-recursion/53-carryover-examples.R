## docs/53: the carryover recursion in the published pmsimstats code.
##
## Regenerates every number, table and figure in
## docs/53-carryover-recursion-in-published-code.md. Run from the
## repository root:
##   Rscript <path>/53-carryover-examples.R [out_dir]
##
## Sources: the 58b32a9 files under
## analysis/data/quick-sim/hendrickson-58b32a9/R (byte-identical to
## github.com/rchendrickson/pmsimstats at 58b32a9) and the vendored 2024
## upstream generateData() under analysis/scripts/quick-sim/
## hendrickson-original-comparison/.

suppressMessages({
  library(data.table)
  library(MASS)
  library(corpcor)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else file.path(tempdir(), 'docs53')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- 'docs/figures'
src <- 'analysis/data/quick-sim/hendrickson-58b32a9/R'
head_file <- paste0('analysis/scripts/quick-sim/',
  'hendrickson-original-comparison/vendored-hendrickson-generateData.R')
stopifnot(file.exists(file.path(src, 'generateData.R')), file.exists(head_file))

cat('md5 generateData.R (58b32a9):',
  tools::md5sum(file.path(src, 'generateData.R')), '\n')

pub <- new.env()
for (f in c('utilities.R', 'buildtrialdesign.R', 'generateData.R')) {
  sys.source(file.path(src, f), envir = pub)
}

## Anchored patch applied to the published text: lines 83-92 replaced.
gd <- readLines(file.path(src, 'generateData.R'))
stopifnot(grepl('brmeans<-modgompertz', gd[83]),
  grepl('brmeans[p-1]', gd[88], fixed = TRUE), grepl('^      }$', gd[92]))
fix_lines <- c(
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
gd_fix <- c(gd[1:82], fix_lines, gd[93:length(gd)])
fix_file <- file.path(out_dir, 'generateData-anchored.R')
writeLines(gd_fix, fix_file)
fixed <- new.env()
for (f in c('utilities.R', 'buildtrialdesign.R')) {
  sys.source(file.path(src, f), envir = fixed)
}
sys.source(fix_file, envir = fixed)

## One-line alternative: decay by the increment in t_sd since the
## previous visit. Exact for exponential decay only (memoryless).
gd_inc <- gd
gd_inc[88] <- paste0('              brmeans[p]<-brmeans[p]+brmeans[p-1]*',
  '(1/2)^((d$tsd[p]-d$tsd[p-1])/modelparam$carryover_t1half)')
inc_file <- file.path(out_dir, 'generateData-increment.R')
writeLines(gd_inc, inc_file)
incr <- new.env()
for (f in c('utilities.R', 'buildtrialdesign.R')) {
  sys.source(file.path(src, f), envir = incr)
}
sys.source(inc_file, envir = incr)

upstream <- new.env()
for (f in c('utilities.R', 'buildtrialdesign.R')) {
  sys.source(file.path(src, f), envir = upstream)
}
sys.source(head_file, envir = upstream)

## The repository's data/extracted_rp.rda and extracted_bp.rda are
## identical to the 58b32a9 data/ files (checked with identical()).
load('data/extracted_rp.rda')
load('data/extracted_bp.rda')
rp_br <- extracted_rp[cat == 'br']
cat('BR Gompertz parameters: max', rp_br$max, 'disp', rp_br$disp,
  'rate', rp_br$rate, '\n')

bt <- pub$buildtrialdesign
designs <- list(
  Hybrid = bt(name_longform = 'primary N-of-1 design',
    name_shortform = 'N-of-1', timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
    timeptnames = c('OL1', 'OL2', 'BD1', 'BD2', 'BD3', 'BD4', 'COd', 'COp'),
    expectancies = c(1, 1, .5, .5, .5, .5, .5, .5),
    ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
      pathB = c(1, 1, 1, 1, 0, 0, 0, 1), pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
      pathD = c(1, 1, 1, 0, 0, 0, 0, 1))),
  OLBDC = bt(name_longform = 'open label+blinded discontinuation',
    name_shortform = 'OL+BDC', timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
    timeptnames = c('OL1', 'OL2', 'OL3', 'OL4', 'BD1', 'BD2', 'BD3', 'BD4'),
    expectancies = c(1, 1, 1, 1, .5, .5, .5, .5),
    ondrug = list(pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
      pathB = c(1, 1, 1, 1, 1, 0, 0, 0))),
  CO = bt(name_longform = 'traditional crossover', name_shortform = 'CO',
    timepoints = pub$cumulative(rep(2.5, 8)),
    timeptnames = c(paste('COa', 1:4, sep = ''), paste('COb', 1:4, sep = '')),
    expectancies = rep(.5, 8),
    ondrug = list(pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
      pathB = c(0, 0, 0, 0, 1, 1, 1, 1))))

## Lines 83-92 of the published generateData(), restated as a function;
## sf is the 2024 scale factor (1 reproduces 58b32a9).
published_mean <- function(td, t_half, sf = 1) {
  d <- data.table(td)
  d[, onDrug := (tod > 0)]
  nP <- nrow(d)
  brmeans <- pub$modgompertz(d$tod, rp_br$max, rp_br$disp, rp_br$rate)
  for (p in 2:nP) {
    if (!d[p]$onDrug && d[p]$tsd > 0) {
      brmeans[p] <- brmeans[p] + brmeans[p - 1] *
        (1/2)^(sf * d$tsd[p] / t_half)
    }
  }
  brmeans
}

## The intended model: mu_last * 2^(-t_sd / t_half), written directly.
intended_mean <- function(td, t_half) {
  g <- pub$modgompertz(td$tod, rp_br$max, rp_br$disp, rp_br$rate)
  out <- g
  last <- 0
  for (p in seq_along(g)) {
    if (td$tod[p] > 0) {
      last <- g[p]
    } else if (td$tsd[p] > 0) {
      out[p] <- last * 2^(-td$tsd[p] / t_half)
    }
  }
  out
}

## Empirical draw: with empirical = TRUE, mvrnorm returns a sample whose
## column means equal mu exactly, so the code's own means can be read off.
code_mean <- function(env, td, t_half, ...) {
  mp <- data.frame(N = 200, c.bm = 0, carryover_t1half = t_half,
    c.tv = .8, c.pb = .8, c.br = .8, c.cf1t = .2, c.cfct = .1)
  dat <- env$generateData(mp, extracted_rp, extracted_bp, td,
    empirical = TRUE, makePositiveDefinite = TRUE, ...)
  as.numeric(colMeans(dat[, paste(td$timeptnames, 'br', sep = '.'),
    with = FALSE]))
}

halves <- c(0.1, 0.2, 0.5, 1, 2)
checks <- list()
for (dn in names(designs)) {
  for (k in seq_along(designs[[dn]]$trialpaths)) {
    td <- designs[[dn]]$trialpaths[[k]]
    for (h in halves) {
      checks[[length(checks) + 1]] <- data.table(design = dn, path = k,
        t_half = h,
        replica_vs_code = max(abs(published_mean(td, h) -
          code_mean(pub, td, h))),
        fixed_vs_intended = max(abs(intended_mean(td, h) -
          code_mean(fixed, td, h))),
        oneline_vs_intended = max(abs(intended_mean(td, h) -
          code_mean(incr, td, h))),
        upstream_vs_sf2 = max(abs(published_mean(td, h, sf = 2) -
          suppressWarnings(code_mean(upstream, td, h)))),
        same_nonzero = identical(published_mean(td, h) != 0,
          intended_mean(td, h) != 0))
    }
  }
}
checks <- rbindlist(checks)
cat('\nVerification of the replica, the fix and the 2024 code\n')
print(checks[, .(max_replica = max(replica_vs_code),
  max_fixed = max(fixed_vs_intended),
  max_oneline = max(oneline_vs_intended), max_upstream = max(upstream_vs_sf2),
  all_same_nonzero = all(same_nonzero))])
fwrite(checks, file.path(out_dir, 'verification.csv'))

f3 <- function(x) formatC(x, format = 'f', digits = 3)
fr <- function(x) ifelse(is.na(x), '', formatC(x, format = 'f', digits = 3))

path_table <- function(td, design_name, path_label, hs) {
  base <- data.table(visit = td$timeptnames, week = pub$cumulative(td$t_wk),
    D = as.integer(td$tod > 0), t_od = td$tod, t_sd = td$tsd,
    gompertz = pub$modgompertz(td$tod, rp_br$max, rp_br$disp, rp_br$rate))
  for (h in hs) {
    pm <- published_mean(td, h)
    am <- intended_mean(td, h)
    base[, paste0('pub_', h) := pm]
    base[, paste0('anch_', h) := am]
    base[, paste0('ratio_', h) := ifelse(D == 0 & t_sd > 0, pm / am, NA)]
  }
  base[, design := design_name]
  base[, path := path_label]
  base
}

ex1 <- path_table(designs$Hybrid$trialpaths[[1]], 'Hybrid', 'path 1',
  c(1, 0.5))
ex2 <- path_table(designs$CO$trialpaths[[1]], 'CO', 'path 1 (drug first)',
  c(1, 2))
fwrite(ex1, file.path(out_dir, 'example1.csv'))
fwrite(ex2, file.path(out_dir, 'example2.csv'))

md_table <- function(dt, hs) {
  hdr <- c('Visit', 'Week', '$D$', '$t_{od}$', '$t_{sd}$', '$G(t_{od})$')
  for (h in hs) {
    hdr <- c(hdr, sprintf('Pub., %s', h), sprintf('Anch., %s', h),
      sprintf('Ratio, %s', h))
  }
  rows <- apply(dt, 1, function(r) NULL)
  lines <- paste('|', paste(hdr, collapse = ' | '), '|')
  lines <- c(lines, paste0('|', paste(rep('---', length(hdr)),
    collapse = '|'), '|'))
  for (i in seq_len(nrow(dt))) {
    v <- c(dt$visit[i], format(dt$week[i]), dt$D[i], format(dt$t_od[i]),
      format(dt$t_sd[i]), f3(dt$gompertz[i]))
    for (h in hs) {
      v <- c(v, f3(dt[[paste0('pub_', h)]][i]),
        f3(dt[[paste0('anch_', h)]][i]), fr(dt[[paste0('ratio_', h)]][i]))
    }
    lines <- c(lines, paste('|', paste(v, collapse = ' | '), '|'))
  }
  lines
}

sink(file.path(out_dir, 'tables.md'))
cat('## Example 1\n')
writeLines(md_table(ex1, c(1, 0.5)))
cat('\n## Example 2\n')
writeLines(md_table(ex2, c(1, 2)))
sink()

## Step-by-step arithmetic for Example 1 at t_half = 1 and 0.5.
cat('\nExample 1 arithmetic\n')
td1 <- designs$Hybrid$trialpaths[[1]]
g1 <- pub$modgompertz(td1$tod, rp_br$max, rp_br$disp, rp_br$rate)
print(data.table(visit = td1$timeptnames, gompertz = sprintf('%.6f', g1)))
for (h in c(1, 0.5)) {
  cat(sprintf('t_half = %s: BD3 %.4f * 2^-(%g/%g) = %.4f\n', h, g1[4],
    td1$tsd[5], h, g1[4] * 2^(-td1$tsd[5] / h)))
  cat(sprintf('  BD4 published %.4f * 2^-(%g/%g) = %.4f; anchored %.4f',
    g1[4] * 2^(-td1$tsd[5] / h), td1$tsd[6], h,
    g1[4] * 2^(-(td1$tsd[5] + td1$tsd[6]) / h), g1[4] * 2^(-td1$tsd[6] / h)))
  cat(sprintf('; COp %.4f * 2^-(%g/%g) = %.4f\n', g1[7], td1$tsd[8], h,
    g1[7] * 2^(-td1$tsd[8] / h)))
}

## Summary over every path of the three designs with an off-drug run.
summ <- list()
for (dn in names(designs)) {
  for (k in seq_along(designs[[dn]]$trialpaths)) {
    td <- designs[[dn]]$trialpaths[[k]]
    off <- td$tod == 0 & td$tsd > 0
    if (!any(off)) next
    for (h in halves) {
      pm <- published_mean(td, h)[off]
      am <- intended_mean(td, h)[off]
      summ[[length(summ) + 1]] <- data.table(design = dn, path = k,
        t_half = h, n_off = sum(off), min_ratio = min(pm / am),
        max_abs_diff = max(am - pm), max_anch = max(am))
    }
  }
}
summ <- rbindlist(summ)
fwrite(summ, file.path(out_dir, 'summary.csv'))
cat('\nSummary by design and half-life (worst path)\n')
summ_d <- summ[, .(min_ratio = min(min_ratio),
  max_abs_diff = max(max_abs_diff)), by = .(design, t_half)]
print(summ_d, digits = 4)
fwrite(summ_d, file.path(out_dir, 'summary-by-design.csv'))

## Upstream 2024 code (scale factor 2), Hybrid path 1.
cat('\nUpstream (scalefactor = 2) on Hybrid path 1\n')
for (h in c(1, 0.5)) {
  print(data.table(visit = td1$timeptnames, t_half = h,
    pub_sf1 = round(published_mean(td1, h), 4),
    pub_sf2 = round(published_mean(td1, h, sf = 2), 4),
    anchored = round(intended_mean(td1, h), 4)))
}

## Re-exposure: Hybrid paths 3 and 4 have a three- or four-visit run; path
## 3 is re-exposed at week 16 with t_od restarting at zero.
cat('\nRe-exposure step on Hybrid path 3 (week 12 off, week 16 on)\n')
td3 <- designs$Hybrid$trialpaths[[3]]
for (h in c(0.5, 1, 2, 4)) {
  cat(sprintf('t_half %g: week 12 published %.3f anchored %.3f; week 16 %.3f\n',
    h, published_mean(td3, h)[6], intended_mean(td3, h)[6],
    intended_mean(td3, h)[7]))
}

## Figures.
long_fig <- function(td, hs, design_label) {
  rbindlist(lapply(hs, function(h) {
    wk <- c(0, pub$cumulative(td$t_wk))
    rbind(
      data.table(week = wk, mean = c(0, published_mean(td, h)),
        process = 'Published recursion (58b32a9)',
        panel = sprintf('t[1/2] == %s ~ "week%s"', h,
          ifelse(h == 1, '', 's'))),
      data.table(week = wk, mean = c(0, intended_mean(td, h)),
        process = 'Intended (anchored) decay',
        panel = sprintf('t[1/2] == %s ~ "week%s"', h,
          ifelse(h == 1, '', 's'))))
  }))
}

on_spans <- function(td) {
  wk <- c(0, pub$cumulative(td$t_wk))
  od <- td$tod > 0
  data.table(xmin = wk[-length(wk)][od], xmax = wk[-1][od])
}

pal <- c('Intended (anchored) decay' = '#2a78d6',
  'Published recursion (58b32a9)' = '#eb6834')
lty <- c('Intended (anchored) decay' = 'solid',
  'Published recursion (58b32a9)' = 'dashed')
shp <- c('Intended (anchored) decay' = 16,
  'Published recursion (58b32a9)' = 17)

plot_paths <- function(dt, spans, xlab, file) {
  dt[, panel := factor(panel, levels = unique(panel))]
  p <- ggplot(dt, aes(week, mean, colour = process, linetype = process,
      shape = process)) +
    geom_rect(data = spans, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = 'grey92') +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.2) +
    facet_wrap(~panel, labeller = label_parsed) +
    scale_colour_manual(values = pal) +
    scale_linetype_manual(values = lty) +
    scale_shape_manual(values = shp) +
    labs(x = xlab, y = 'Expected drug response, BR mean (points)',
      colour = NULL, linetype = NULL, shape = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = 'bottom', panel.grid.minor = element_blank())
  ggsave(file, p, width = 7.5, height = 3.8, dpi = 200, bg = 'white')
}

plot_paths(long_fig(td1, c(1, 0.5)), on_spans(td1),
  'Week (shaded: on drug)', file.path(fig_dir, '53-fig1-hybrid-path1.png'))
td_co <- designs$CO$trialpaths[[1]]
plot_paths(long_fig(td_co, c(1, 2)), on_spans(td_co),
  'Week (shaded: on drug)', file.path(fig_dir, '53-fig2-crossover-path1.png'))
cat('\nwrote figures and', out_dir, '\n')
