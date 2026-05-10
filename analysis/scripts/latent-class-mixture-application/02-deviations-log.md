# Deviations log: latent-class-mixture-application

*2026-05-09 19:40 PDT*

This document records every deviation from the ADEMP
pre-registration in `00-ademp-pre-registration.md`, with date,
reason, and (where applicable) reviewer signoff. Per the
pre-registration §"Deviations from this pre-registration",
deviations are recorded before the affected production runs
are committed.

## D1. Phase B reduced replicate counts (kill-switch budget)

**Date.** 2026-05-09.

**Pre-registered specification.** "$n_{\text{reps}} = 5000$ per
cell for type I error and class-count recovery estimands"
(common infrastructure §"Replicates").

**Deviation.** Phase B production driver
`10-study5-production.R` uses 2000 replicates per Type I error
cell (Cell 1) and 500 replicates per class-count recovery cell
(Cells 2 and 3), rather than the pre-registered 5000 per cell.

**Reason.** Phase B is the kill-switch decision per
`~/prj/res/00-portfolio/simgo-paper03-phasing.md`: a feasibility
check that establishes whether `lcmm::hlme` is usable at trial-
relevant $N$ at all. Reducing replicates from 5000 to 2000 (Type
I) and 500 (class-count) keeps wall-clock under two hours at
8 cores while preserving sufficient MCSE precision to evaluate
the kill-switch criteria:

- Type I MCSE at $\alpha = 0.05$ with 2000 reps: $\sqrt{0.05
  \times 0.95 / 2000} = 0.0049$, well under the
  $\pm 0.01$ kill-switch tolerance.
- Class-count selection MCSE at $p = 0.80$ with 500 reps:
  $\sqrt{0.80 \times 0.20 / 500} = 0.0179$, sufficient to
  resolve a ten-percentage-point gap to the 0.80 threshold.

**Plan.** If Phase B passes the kill-switch criteria, the full
5000-rep Type I and class-count cells will be re-run during
Phase C alongside the Study 1 production runs, recovering the
pre-registered replicate counts.

## D2. Phase B class-count grid restricted to $K \in \{1, 2\}$

**Date.** 2026-05-09.

**Pre-registered specification.** "$K \in \{1, 2, 3\}$ classes
for cells 2 and 3" (Study 5 §M).

**Deviation.** The Phase B production driver fits `lcmm::hlme`
at $K = 1$ and $K = 2$ only, not $K = 3$, via the
`lcmm_analysis()` wrapper (which is committed at $K \in \{1,
2\}$).

**Reason.** $K = 2$ versus $K = 1$ is the load-bearing
identifiability decision: if `hlme` cannot reliably distinguish
single-component from two-class data, then class-aware analyses
of any kind are unworkable at trial-relevant $N$, and $K = 3$
provides no additional Phase B signal. Adding $K = 3$ would
roughly double the per-replicate fit time and is more
informatively run during Phase C alongside Study 1's mixture-
DGP cells, where $K = 3$ overfit-rate is the substantive
question.

**Plan.** The $K = 3$ extension is folded into Phase C of
`simgo-paper03-phasing.md`. The lcmm wrapper will be extended
to $K \in \{1, 2, 3\}$ at the start of Phase C, before the
Study 1 production runs commit.

## D3. Phase B Cell 2 implementation: post-hoc residual transform

**Date.** 2026-05-09.

**Pre-registered specification.** "Heavy-tailed and mildly
skewed residual distributions on $BR$: $t_5$ distributed
standardised residuals, and (separately) skew-normal with
shape parameter 4" (Study 5 §D, cell 2).

**Deviation.** The Phase B production driver implements the
non-Gaussian residual cells via a post-hoc transform of the
MVN draw rather than via a fully non-Gaussian DGP. Specifically,
the BR component at each timepoint is standardised, mapped to
its $t_5$ or skew-normal quantile via probability-integral
transform (with the skew-normal approximated by a delta-tilted
absolute-normal mixture), then re-standardised to preserve the
marginal SD.

**Reason.** A fully non-Gaussian DGP would require either a
multivariate $t$ or multivariate skew-normal draw with the
existing AR(1) structure preserved, neither of which is in the
package source. The post-hoc transform preserves the marginal
distribution targets that the pre-registration specifies and is
sufficient for the Phase B diagnostic question (does mild
non-Gaussianity break `hlme` identifiability or class-count
selection?). The transform alters the full multivariate
distribution at each timepoint independently, so the
within-subject correlation structure is no longer preserved
exactly, but the diagnostic interpretation of class-count
selection at the null is unaffected because BIC and Wald
calibration depend on the marginal distribution.

**Plan.** If Phase B exposes meaningful sensitivity to
non-Gaussianity (e.g. Type I error inflation under $t_5$), a
fully non-Gaussian multivariate DGP will be implemented for
Phase C. If Phase B shows little sensitivity, the post-hoc
transform is retained for production.

---

*Source: ~/prj/alz/10-pmsimstats-ng/pmsimstats-ng/analysis/scripts/latent-class-mixture-application/02-deviations-log.md*
