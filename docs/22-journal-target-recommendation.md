# Journal Target Recommendation: '02' DGP Architecture Comparison {.unlisted .unnumbered}

*2026-04-08*

**Author.** pmsimstats team

**Purpose.** This note recommends target journals for the DGP
architecture comparison manuscript (then numbered '02', now paper 01,
`analysis/report/01-dgp-mean-moderation-vs-mvn/`). It was generated on
2026-04-08 from the Zotero nof1 collection (50 items) and a reading of
the manuscript draft of that date
(`docs/02-dgp-mean-moderation-vs-mvn.tex`, which is no longer in the
repository; only an annotated PDF of that draft remains in `docs/`).
It was
reviewed on 2026-10-07 for notation, reference accuracy and project
conventions.

```{=latex}
\clearpage
\tableofcontents
\listoftables
\clearpage
```

## Notation index and glossary

This note uses no mathematical symbols. Architecture names follow
`analysis/report/NOTATION.md`, Part 3.

| Term | Meaning | Status |
|---|---|---|
| mean-moderation architecture | the biomarker scales the treatment effect in the conditional mean (formerly Architecture A; data value `mean_moderation`) | canonical |
| covariance-moderation architecture | the biomarker-response correlation differs by treatment state within a multivariate normal (MVN) draw (formerly Architecture B; data value `mvn`) | canonical |
| DGP | data-generating process | canonical |
| Path B | the strengthened revision plan for the manuscript as of 2026-04-08 (1,000 or more replicates, evaluated alternative analysis strategies, power curves) | local |
| collection precedent | a paper in the Zotero nof1 collection published in the named journal | local |

Table: Terms used in this note, with meaning and status

**Glossary.** *N-of-1 trial*: a multiple-crossover trial in a single
participant; *aggregated N-of-1 trial*: many such series pooled in one
mixed model. *Carryover*: persistence of a drug effect after
discontinuation.

## Manuscript summary

*Note (2026-10-07).* This section describes the manuscript as it stood
on 2026-04-08. The power-loss figures below are superseded and should
not be quoted. When the comparison was re-run at a matched total of
$N = 70$ for all designs (2026-08-08), the Hybrid-design carryover loss
was 6.9 percent under mean moderation and 15.4 percent under covariance
moderation (`analysis/report/NOTATION.md`, Part 2). Paper 01 has since
been re-scoped as a refined data-generating process following
Hendrickson et al. (2020): the published step coupling makes the
coupling gap zero, so the interaction vanishes at any positive
half-life, and graded coupling ($2^{-t_{sd}/t_{1/2}}$ off drug, 0 before
exposure) removes the collapse (docs/37); under graded coupling,
covariance-based binding and mean moderation that fades with the drug
effect lose power almost identically (paper 01, abstract). The journal
reasoning below should be revisited against the current manuscript
before it is acted on.

'Two Architectures for Simulating Biomarker-Treatment Interactions:
Implications for Statistical Power Under Carryover' formalizes the
distinction between direct mean moderation (the mean-moderation
architecture) and differential correlation within an MVN draw (the
covariance-moderation architecture) as data-generating processes for
biomarker-treatment interaction simulations. The central finding, as
stated in the 2026-04-08 draft, was that the choice of DGP
architecture determines whether carryover effects reduce statistical
power modestly (~10-15%, mean moderation) or substantially (~40-60%,
covariance moderation). The manuscript provides mathematical
formalization, biological rationale for each architecture, simulation
evidence from N-of-1 trial designs, and guidance for DGP selection.

## Methodological positioning

The manuscript sits at the intersection of three domains:

- Statistical simulation methodology (DGP specification and its
  consequences)
- N-of-1 and crossover trial design (carryover, power analysis)
- Biomarker-treatment interaction modeling (predictive biomarker
  validation)

## Most similar papers in the Zotero nof1 collection

| Paper | Similarity | Journal |
|---|---|---|
| Shi & Ye (2024): behavioral carryover and power in crossover trials | Highest | *Biometrics* |
| Blackston et al. (2019): simulation comparison, aggregated N-of-1 vs. RCT | High | *Healthcare* |
| Wang & Schork (2019): power and design in crossover N-of-1 | High | *Healthcare* |
| Senn (2016): variance components and personalized medicine | High | *Statistics in Medicine* |
| Schork (2022): serial correlation in personalized studies | High | *Harvard Data Science Review* |
| Haller & Ulm (2018): biomarker-treatment interaction estimation | High | *Trials* |
| Senn (2019): sample size for N-of-1 (see note) | Moderate | *Statistical Methods in Medical Research* |
| Tang & Landes (2020): t-tests for N-of-1 with serial correlation | Moderate | *PLoS ONE* |
| Hawksworth et al. (2024): methodological review of randomized N-of-1 trials | Moderate | *Trials* |

Table: Most similar papers in the Zotero nof1 collection with similarity and journal

*Attribution corrections (2026-10-07).* Six citations in the original
table and in the ranked list below were misattributed or could not be
found when checked against PubMed, Crossref and the Zotero library.
Three of them match Zotero attachment file names rather than the item
metadata, which is the likely source of the error. "Yap, Wang & Harhay (2024)"
in *Biometrics* could not be found; the collection's 2024 *Biometrics*
paper on carryover and power in crossover trials is Shi and Ye (2024),
and the entry is taken to be that paper. "Chen et al. (2024)" is
Hawksworth et al. (2024). "Senn & Julious (2024)" is Senn (2024), a
single-author paper. "Mitchell et al. (2024)" is Chatters et al.
(2024). "Arends et al. (2019)" in *Statistical Methods in Medical
Research* is Senn (2019): the project bibliography entry behind it
carries the title of Senn's paper with an invented author list
(verified against PubMed 2026-10-07). "Liu et al. (2021)" in
*Contemporary Clinical Trials Communications* does not exist: its
bibliography entry carried the title of Lim and In (2021), "Considerations
for crossover design in clinical study", published in the *Korean
Journal of Anesthesiology*, with invented authors and a DOI that
resolves to an unrelated paper in *Contemporary Clinical Trials
Communications* (verified against PubMed and Crossref 2026-10-07). It
is therefore not a precedent for that journal.

## Ranked journal targets

### Rank 1: *Statistics in Medicine*

- **Collection precedent:** Senn (2016)
- **Rationale:** Premier journal for statistical methodology in
  clinical trial design. The core contribution (in the 2026-04-08
  draft, that DGP architecture choice produces 40-60% divergent power
  estimates under carryover; superseded, see the note under the
  manuscript summary) is a methodological finding with implications
  for how the field designs and reports simulation studies. The
  journal's readership (biostatisticians designing trials) is the
  primary audience who would act on this finding. The Path B
  strengthened version (1,000+ replicates, evaluated alternative
  analysis strategies, power curves) meets the journal's
  expectations for rigor.
- **Risk:** Competitive acceptance rate. Must demonstrate the
  finding generalizes beyond the specific pmsimstats
  parameterization.

### Rank 2: *Biometrics*

- **Collection precedent:** Shi & Ye (2024)
- **Rationale:** Shi & Ye (2024) is the closest comparator paper
  in the collection: carryover's impact on power in crossover
  designs, published in *Biometrics* in 2024. The '02' manuscript
  extends this by showing that the mechanism by which the
  biomarker-treatment interaction is generated determines the
  magnitude of power loss. Natural companion to Shi & Ye.
- **Risk:** Expects stronger theoretical results (formal proofs,
  asymptotic analysis). May view the scope as too narrow for a
  generalist biometrics readership. Would need broader framing in
  terms of covariance-vs-mean signal detection.

### Rank 3: *Trials*

- **Collection precedent:** Haller & Ulm (2018), Hawksworth et al.
  (2024), Senn (2024), Chatters et al. (2024)
- **Rationale:** Highest concentration of N-of-1 methodology papers
  in the collection (4 of 50). Haller & Ulm (2018) is essentially
  the same methodological genre as '02' (simulation study of
  biomarker-treatment interaction estimation). Open access,
  practical orientation, readership includes both methodologists
  and trialists.
- **Risk:** Lower impact factor. Theoretical depth may exceed the
  journal's typical level, though this is minor.

### Rank 4: *Statistical Methods in Medical Research*

- **Collection precedent:** Senn (2019) (listed originally as
  "Arends et al. (2019)"; see the attribution note)
- **Rationale:** Publishes simulation-based methodological work for
  clinical trial design. Senn (2019) is the same genre
  (power/sample size for N-of-1). Occupies middle ground between
  theoretical focus of *Biometrics* and practical focus of *Trials*.
  Accepts longer methodological papers with extensive simulation.
- **Risk:** Slower review times. Less visibility in the N-of-1
  community.

### Rank 5: *Contemporary Clinical Trials Communications*

- **Collection precedent:** Sturdevant & Lumley (2021). (The
  original "Liu et al. (2021)" precedent was a misattributed entry;
  see the attribution note.)
- **Rationale:** Relevant work on carryover testing (Sturdevant &
  Lumley). Good fallback. Practical framing and clinical
  application (prazosin-PTSD) fit well.
- **Risk:** Lower visibility and impact than journals above.

### Rank 6: *Harvard Data Science Review*

- **Collection precedent:** Schork (2022)
- **Rationale:** Schork (2022) on related statistical machinery was
  published here. HDSR publishes methodological perspectives.
- **Risk:** More of a review/perspective venue. Heavy simulation
  results in the Path B version may not fit format expectations.

### Rank 7: *Healthcare* (MDPI)

- **Collection precedent:** Blackston et al. (2019), Wang & Schork
  (2019)
- **Rationale:** Two of the most similar simulation-comparison
  papers were published here. Fast review, open access.
- **Risk:** Lower standing of MDPI journals among methodologists. If
  strengthened to Path B standards, the manuscript would suit a
  higher-impact venue.

## Primary recommendation

The rankings in this note are editorial judgments (argued), based on
collection precedent and fit; no acceptance rates, impact factors or
review times were looked up.

**Target: *Statistics in Medicine*.** The manuscript's central claim
is about simulation methodology itself. *Statistics in Medicine*
reaches the audience most likely to change their practice.

**Backup: *Trials*.** If *Statistics in Medicine* considers the
scope too specialized, *Trials* has the most active N-of-1
methodology pipeline and strong fit with collection precedents.

## References

Metadata verified on 2026-10-07 against PubMed or Crossref, except
where marked.

- Blackston, J. W., Chapple, A. G., McGree, J. M., McDonald, S. and
  Nikles, J. (2019). Comparison of aggregated N-of-1 trials with
  parallel and crossover randomized controlled trials using
  simulation studies. *Healthcare* 7(4):137.
  doi:10.3390/healthcare7040137.
- Chatters, R., Hawksworth, O., Julious, S. and Cook, A. (2024). The
  development of a set of key points to aid clinicians and
  researchers in designing and conducting n-of-1 trials. *Trials*
  (published 2024-07-11; PMID 38992786; a correction appeared
  2024-09-09).
- Haller, B. and Ulm, K. (2018). A simulation study on estimating
  biomarker-treatment interaction effects in randomized trials with
  prognostic variables. *Trials* 19:128.
  doi:10.1186/s13063-018-2491-0.
- Hawksworth, O., Chatters, R., Julious, S., Cook, A., Biggs, K.,
  Solaiman, K., Quah, M. C. H. and Cheong, S. C. (2024). A
  methodological review of randomised n-of-1 trials. *Trials* 25:263.
  doi:10.1186/s13063-024-08100-1.
- Hendrickson, R. C., Thomas, R. G., Schork, N. J. and Raskind, M. A.
  (2020). Optimizing aggregated N-of-1 trial designs for predictive
  biomarker validation: statistical methods and theoretical findings.
  *Frontiers in Digital Health* 2:13. doi:10.3389/fdgth.2020.00013.
- Lim, C.-Y. and In, J. (2021). Considerations for crossover design
  in clinical study. *Korean Journal of Anesthesiology* 74(4):293-299.
  doi:10.4097/kja.21165. (The paper behind the misattributed
  "Liu et al. (2021)" entry.)
- Schork, N. J. (2022). Accommodating serial correlation and
  sequential design elements in personalized studies and aggregated
  personalized studies. *Harvard Data Science Review*, Special
  Issue 3. doi:10.1162/99608f92.f1eef6f4.
- Senn, S. (2016). Mastering variation: variance components and
  personalised medicine. *Statistics in Medicine* 35:966-977.
  doi:10.1002/sim.6739.
- Senn, S. (2019). Sample size considerations for n-of-1 trials.
  *Statistical Methods in Medical Research* (online 2017).
  doi:10.1177/0962280217726801. The paper behind the misattributed
  "Arends et al. (2019)" entry (verified by title).
- Senn, S. (2024). The analysis of continuous data from n-of-1 trials
  using paired cycles: a simple tutorial. *Trials* (published
  2024-02-16; PMID 38365817).
- Shi, D. and Ye, T. (2024). Behavioral carry-over effect and power
  consideration in crossover trials. *Biometrics* 80(2):ujae023.
  doi:10.1093/biomtc/ujae023.
- Sturdevant, S. G. and Lumley, T. (2021). Statistical methods for
  testing carryover effects: a mixed effects model approach.
  *Contemporary Clinical Trials Communications* 22:100711.
  doi:10.1016/j.conctc.2021.100711.
- Tang, J. and Landes, R. D. (2020). Some t-tests for N-of-1 trials
  with serial correlation. *PLoS ONE* 15(2):e0228077.
  doi:10.1371/journal.pone.0228077.
- Wang, Y. and Schork, N. J. (2019). Power and design issues in
  crossover-based N-of-1 clinical trials with fixed data collection
  periods. *Healthcare* 7(3):84. doi:10.3390/healthcare7030084.
