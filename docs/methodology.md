# Methodology Notes

A reference for the statistical genetics concepts underlying this
pipeline — useful both for your own review before interviews and as
documentation for anyone evaluating the project.

## What is a Polygenic Risk Score?

A PRS aggregates the small, additive effects of many genetic variants
(typically thousands to millions of SNPs) into a single score
estimating an individual's genetic predisposition to a trait or
disease. For each SNP, a genome-wide association study (GWAS) estimates
an effect size (how much that variant shifts disease risk or a trait
value). The PRS is the weighted sum of an individual's genotype at each
SNP, weighted by these effect sizes.

```
PRS_i = Σ (genotype_ij × effect_weight_j)   for all SNPs j
```

## Why PRS calculation is harder than it sounds

1. **SNP matching**: the GWAS discovery SNPs must be matched to SNPs
   present in the target genotype data, accounting for different
   genome builds, strand orientation, and naming conventions.
2. **Linkage disequilibrium (LD)**: nearby SNPs are correlated, so
   simply summing effects across all GWAS-significant SNPs would
   double-count signal from correlated variants. Tools like PRSice-2
   and LDpred2 apply clumping/thresholding or Bayesian shrinkage to
   account for this.
3. **Population/ancestry mismatch**: PRS are typically derived from
   GWAS in one ancestry group (historically, overwhelmingly European-
   ancestry cohorts) and can perform substantially worse — both in
   discrimination and calibration — when applied to individuals from
   different ancestries. This is one of the most important current
   limitations and active research areas in the PRS field.
4. **Effect size shrinkage**: raw GWAS effect sizes are often
   overestimated (winner's curse), especially for variants just
   crossing genome-wide significance. Bayesian methods like LDpred2
   address this directly.

This project's `calculate_prs.py` implements the core weighted-sum
calculation (step in the list above) but does NOT implement LD
clumping or ancestry adjustment — these are flagged as known
simplifications, both in code comments and in `qc_checks.py`, rather
than silently ignored.

## Interpreting PRS results responsibly

A PRS that is statistically significantly associated with disease
status is not automatically clinically useful. Three numbers matter
together:

- **Odds ratio**: how much does risk increase per unit increase in PRS?
  Can look impressive (e.g. 3x per SD) while still corresponding to
  modest absolute risk change if the baseline disease prevalence is low.
- **AUC**: how well does the PRS distinguish cases from controls at
  the individual level? Most published PRS for common complex diseases
  have AUCs in the 0.55–0.70 range — modest discrimination, not
  diagnostic-grade.
- **Variance explained (R²)**: how much of the total phenotypic
  variance does the PRS account for? Typically a few percent to ~20%
  for the best-powered PRS currently available, since common complex
  diseases also involve substantial environmental and unmeasured
  genetic contributions.

The honest, defensible claim for most current PRS is: *useful for
population-level risk stratification and research, not sufficient
alone for individual diagnostic decisions.*

## The liability-threshold model used in the simulation

`simulate_genotypes.py` uses a standard liability-threshold model:
an underlying continuous genetic + environmental "liability" is
simulated, and individuals above a chosen quantile threshold are
labeled as cases. This is a standard simplification used in
quantitative and statistical genetics to simulate complex,
binary disease traits with realistic polygenic architecture, and
matches the conceptual model underlying tools like GCTA-GWAS-simulation
and PLINK's simulation utilities.

## Honest limitations of this project

- Uses simulated, not real, genotype/GWAS data — chosen deliberately to
  avoid any data privacy or consent ambiguity, but means results don't
  reflect the messiness of real genotyping data (batch effects, real
  LD structure, population stratification).
- Does not implement LD clumping/thresholding or Bayesian PRS methods
  (LDpred2-style) — the weighted sum here assumes pre-selected,
  independent causal variants, which is a simplification.
- Ancestry/population matching is flagged but not implemented, since
  it requires reference population panels (e.g. 1000 Genomes) not
  included in this simplified version.
- The PGS Catalog API integration (`pgs_catalog_client.py`) has not
  been verified against a live network call in this development
  environment — see the caveat in that file's docstring.
