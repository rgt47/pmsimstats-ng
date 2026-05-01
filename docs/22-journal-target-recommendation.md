# Journal Target Recommendation: '02' DGP Architecture Comparison

## Manuscript summary

'Two Architectures for Simulating Biomarker-Treatment Interactions:
Implications for Statistical Power Under Carryover' formalizes the
distinction between direct mean moderation (Architecture A) and
differential correlation via MVN (Architecture B) as
data-generating processes for biomarker-treatment interaction
simulations. The central finding is that the choice of DGP
architecture determines whether carryover effects reduce
statistical power modestly (~10-15%, Architecture A) or
substantially (~40-60%, Architecture B). The manuscript provides
mathematical formalization, biological rationale for each
architecture, simulation evidence from N-of-1 trial designs, and
guidance for DGP selection.

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
| Yap, Wang & Harhay (2024) -- carryover and power in crossover | Highest | *Biometrics* |
| Blackston et al. (2019) -- simulation comparison, aggregated N-of-1 vs. RCT | High | *Healthcare* |
| Wang & Schork (2019) -- power and design in crossover N-of-1 | High | *Healthcare* |
| Senn (2016) -- variance components and personalized medicine | High | *Statistics in Medicine* |
| Schork (2022) -- serial correlation in personalized studies | High | *Harvard Data Science Review* |
| Haller & Ulm (2018) -- biomarker-treatment interaction estimation | High | *Trials* |
| Arends et al. (2019) -- sample size for N-of-1 | Moderate | *Statistical Methods in Medical Research* |
| Tang & Landes (2020) -- t-tests for N-of-1 with serial correlation | Moderate | *PLoS ONE* |
| Chen et al. (2024) -- methodological review of N-of-1 | Moderate | *Trials* |

## Ranked journal targets

### Rank 1: *Statistics in Medicine*

- **Collection precedent:** Senn (2016)
- **Rationale:** Premier journal for statistical methodology in
  clinical trial design. The core contribution -- that DGP
  architecture choice produces 40-60% divergent power estimates
  under carryover -- is a methodological finding with implications
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

- **Collection precedent:** Yap, Wang & Harhay (2024)
- **Rationale:** Yap et al. (2024) is the closest comparator paper
  in the collection -- carryover's impact on power in crossover
  designs, published in *Biometrics* in 2024. The '02' manuscript
  extends this by showing that the mechanism by which the
  biomarker-treatment interaction is generated determines the
  magnitude of power loss. Natural companion to Yap et al.
- **Risk:** Expects stronger theoretical results (formal proofs,
  asymptotic analysis). May view the scope as too narrow for a
  generalist biometrics readership. Would need broader framing in
  terms of covariance-vs-mean signal detection.

### Rank 3: *Trials*

- **Collection precedent:** Haller & Ulm (2018), Chen et al.
  (2024), Senn & Julious (2024), Mitchell et al. (2024)
- **Rationale:** Highest concentration of N-of-1 methodology papers
  in the collection (4 of 50). Haller & Ulm (2018) is essentially
  the same methodological genre as '02' (simulation study of
  biomarker-treatment interaction estimation). Open access,
  practical orientation, readership includes both methodologists
  and trialists.
- **Risk:** Lower impact factor. Theoretical depth may exceed the
  journal's typical level, though this is minor.

### Rank 4: *Statistical Methods in Medical Research*

- **Collection precedent:** Arends et al. (2019)
- **Rationale:** Publishes simulation-based methodological work for
  clinical trial design. Arends et al. (2019) is the same genre
  (power/sample size for N-of-1). Occupies middle ground between
  theoretical focus of *Biometrics* and practical focus of *Trials*.
  Accepts longer methodological papers with extensive simulation.
- **Risk:** Slower review times. Less visibility in the N-of-1
  community.

### Rank 5: *Contemporary Clinical Trials Communications*

- **Collection precedent:** Sturdevant & Lumley (2021), Liu et al.
  (2021)
- **Rationale:** Two relevant papers on carryover testing and
  crossover design. Good fallback. Practical framing and clinical
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
- **Risk:** MDPI stigma. If strengthened to Path B standards, the
  manuscript deserves a higher-impact venue.

## Primary recommendation

**Target: *Statistics in Medicine*.** The manuscript's central claim
is about simulation methodology itself. *Statistics in Medicine*
reaches the audience most likely to change their practice.

**Backup: *Trials*.** If *Statistics in Medicine* considers the
scope too specialized, *Trials* has the most active N-of-1
methodology pipeline and strong fit with collection precedents.

---

*Generated 2026-04-08. Based on Zotero nof1 collection (50 items)
and full reading of docs/02-dgp-mean-moderation-vs-mvn.tex.*
