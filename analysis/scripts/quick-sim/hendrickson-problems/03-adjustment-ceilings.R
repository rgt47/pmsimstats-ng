## Ceiling of the biomarker-BR correlation at Hendrickson's published
## parameters (rho = 0.8, c1 = 0.2, cx = 0.1, t1/2 in {0, 0.1, 0.2})
## under each combination of adjustment 1 (CS -> AR(1)) and
## adjustment 2 (step -> graded coupling). Adjustment 3 (mean
## moderation) writes no biomarker correlation into the matrix, so it
## has no ceiling; only the response block must be PD.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/03-adjustment-ceilings.R
## Whitepaper: docs/36-hendrickson-pd-and-carryover.md
library(data.table)
out_dir <- 'analysis/data/quick-sim/hendrickson-problems'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

designs <- list(
  `OL+BDC` = list(w = c(4, 8, 12, 16, 17, 18, 19, 20),
                  paths = list(c(1,1,1,1,1,1,0,0), c(1,1,1,1,1,0,0,0))),
  CO = list(w = cumsum(rep(2.5, 8)),
            paths = list(c(1,1,1,1,0,0,0,0), c(0,0,0,0,1,1,1,1))),
  Hybrid = list(w = c(4, 8, 9, 10, 11, 12, 16, 20),
                paths = list(c(1,1,1,1,0,0,1,0), c(1,1,1,1,0,0,0,1),
                             c(1,1,1,0,0,0,1,0), c(1,1,1,0,0,0,0,1))))

build_M <- function(w, rho, c1, cx, within, cross) {
  n <- length(w); g <- abs(outer(w, w, '-'))
  A <- if (within == 'AR1') rho^g else { m <- matrix(rho, n, n); diag(m) <- 1; m }
  C <- switch(cross,
    const = { m <- matrix(cx, n, n); diag(m) <- c1; m },
    decay = { m <- cx * rho^g; diag(m) <- c1; m },
    separable = c1 * A)
  rbind(cbind(A, C, C), cbind(C, A, C), cbind(C, C, A))
}
## tsd: cumulative time since discontinuation, reset on re-exposure,
## zero before first exposure.
tsd_of <- function(w, on) {
  gaps <- c(w[1], diff(w)); tsd <- numeric(length(on)); ever <- FALSE
  for (t in seq_along(on)) {
    if (on[t] == 1) { tsd[t] <- 0; ever <- TRUE }
    else tsd[t] <- if (ever) (if (t > 1) tsd[t - 1] else 0) + gaps[t] else 0
  }
  tsd
}
coupling <- function(w, on, t_half, rule) {
  tsd <- tsd_of(w, on)
  phi <- if (t_half > 0) 0.5^(tsd / t_half) else rep(0, length(on))
  phi[on == 1 | tsd == 0] <- 0
  off <- if (rule == 'graded') phi else as.numeric(phi > 0)
  ifelse(on == 1, 1, off)
}
ceil <- function(M, u) {
  if (min(eigen(M, TRUE, TRUE)$values) <= 0) return(NA_real_)
  v <- c(rep(0, 2 * length(u)), u)
  1 / sqrt(drop(t(v) %*% solve(M, v)))
}

cfg <- data.table(
  label = c('Published: CS, const cross, step',
            'Adj 1 only: AR(1), decay cross, step',
            'Adj 2 only: CS, const cross, graded',
            'Adj 1+2 (covar): AR(1), decay cross, graded',
            'Adj 1+2, separable: AR(1), c1*A, graded'),
  within = c('CS', 'AR1', 'CS', 'AR1', 'AR1'),
  cross = c('const', 'decay', 'const', 'decay', 'separable'),
  rule = c('step', 'step', 'graded', 'graded', 'graded'))

out <- rbindlist(lapply(seq_len(nrow(cfg)), function(i) {
  rbindlist(lapply(names(designs), function(d) {
    rbindlist(lapply(c(0, 0.1, 0.2), function(th) {
      M <- build_M(designs[[d]]$w, 0.8, 0.2, 0.1, cfg$within[i], cfg$cross[i])
      cs <- sapply(designs[[d]]$paths, function(on)
        ceil(M, coupling(designs[[d]]$w, on, th, cfg$rule[i])))
      data.table(config = cfg$label[i], design = d, t_half = th,
                 ceiling = min(cs))
    }))
  }))
}))
w <- dcast(out, config + design ~ t_half, value.var = 'ceiling')
setnames(w, c('0', '0.1', '0.2'), c('t0', 't0.1', 't0.2'))
w[, config := factor(config, levels = cfg$label)]
setorder(w, config, design)
cat('Ceiling at rho = 0.8, c1 = 0.2, cx = 0.1 (published c.bm: 0.3, 0.6)\n')
print(w, digits = 3)

fwrite(out, file.path(out_dir, 'adjustment-ceilings.csv'))
