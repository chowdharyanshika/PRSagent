"""
simulate_genotypes.py

Generates simulated SNP genotype data, effect sizes, and a phenotype --
entirely synthetic, so the rest of the pipeline can be developed and
tested without touching any real (and potentially sensitive) patient
genotype data.

This mirrors a simplified version of how PRS simulation tools like
PLINK's --simulate or GCTA work, scaled down for a portfolio project.

Model:
  - n_snps additive variants, each with a small effect size drawn from
    a normal distribution (consistent with the polygenic, small-effect
    architecture typical of common complex traits)
  - genotypes coded 0/1/2 (number of risk alleles), drawn from a
    binomial distribution based on a random minor allele frequency (MAF)
    per SNP -- this mimics real population-level allele frequency variation
  - phenotype = weighted sum of genotypes (the "true" genetic liability)
    + environmental noise, then thresholded into a binary case/control
    label -- analogous to a liability-threshold model commonly used in
    quantitative/complex trait genetics
"""

import argparse
import numpy as np
import pandas as pd


def simulate_data(
    n_individuals: int = 2000,
    n_snps: int = 500,
    n_causal_snps: int = 50,
    heritability: float = 0.3,
    seed: int = 42,
):
    """
    Returns:
        genotypes_df: (n_individuals x n_snps) DataFrame, values in {0,1,2}
        effect_sizes_df: (n_snps,) DataFrame with true effect sizes
                          (0 for non-causal SNPs)
        phenotype_df: (n_individuals,) DataFrame with continuous liability
                       and binarized case/control status
    """
    rng = np.random.default_rng(seed)

    # --- Minor allele frequencies, one per SNP, realistic range ---
    mafs = rng.uniform(0.05, 0.5, size=n_snps)

    # --- Genotypes: each SNP ~ Binomial(2, maf) per individual ---
    genotypes = np.column_stack([
        rng.binomial(2, maf, size=n_individuals) for maf in mafs
    ])
    snp_ids = [f"SNP_{i+1}" for i in range(n_snps)]
    genotypes_df = pd.DataFrame(genotypes, columns=snp_ids)
    genotypes_df.insert(0, "individual_id", [f"IND_{i+1}" for i in range(n_individuals)])

    # --- Effect sizes: only n_causal_snps are non-zero (sparse, polygenic) ---
    causal_idx = rng.choice(n_snps, size=n_causal_snps, replace=False)
    effect_sizes = np.zeros(n_snps)
    effect_sizes[causal_idx] = rng.normal(0, 0.2, size=n_causal_snps)
    effect_sizes_df = pd.DataFrame({
        "snp_id": snp_ids,
        "effect_weight": effect_sizes,
        "maf": mafs,
        "is_causal": np.isin(np.arange(n_snps), causal_idx),
    })

    # --- Genetic liability = genotypes @ effect sizes ---
    genetic_liability = genotypes @ effect_sizes

    # --- Scale environmental noise so heritability matches the target ---
    var_genetic = np.var(genetic_liability)
    var_env = var_genetic * (1 - heritability) / heritability if var_genetic > 0 else 1.0
    env_noise = rng.normal(0, np.sqrt(var_env), size=n_individuals)

    total_liability = genetic_liability + env_noise

    # --- Binarize into case/control using a top-quantile threshold ---
    # (mimics a disease with ~20% prevalence in this simulated cohort)
    threshold = np.quantile(total_liability, 0.80)
    case_control = (total_liability >= threshold).astype(int)

    phenotype_df = pd.DataFrame({
        "individual_id": genotypes_df["individual_id"],
        "genetic_liability": genetic_liability,
        "total_liability": total_liability,
        "case_control": case_control,
    })

    return genotypes_df, effect_sizes_df, phenotype_df


def main():
    parser = argparse.ArgumentParser(
        description="Simulate genotype, effect size, and phenotype data for PRS development"
    )
    parser.add_argument("--n_individuals", type=int, default=2000)
    parser.add_argument("--n_snps", type=int, default=500)
    parser.add_argument("--n_causal_snps", type=int, default=50)
    parser.add_argument("--heritability", type=float, default=0.3,
                         help="Proportion of phenotypic variance explained by genetics")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out_dir", type=str, default="data")
    args = parser.parse_args()

    genotypes_df, effect_sizes_df, phenotype_df = simulate_data(
        n_individuals=args.n_individuals,
        n_snps=args.n_snps,
        n_causal_snps=args.n_causal_snps,
        heritability=args.heritability,
        seed=args.seed,
    )

    genotypes_df.to_csv(f"{args.out_dir}/genotypes.csv", index=False)
    effect_sizes_df.to_csv(f"{args.out_dir}/effect_sizes.csv", index=False)
    phenotype_df.to_csv(f"{args.out_dir}/phenotype.csv", index=False)

    print(f"Simulated {args.n_individuals} individuals x {args.n_snps} SNPs")
    print(f"  Causal SNPs: {args.n_causal_snps}")
    print(f"  Target heritability: {args.heritability}")
    print(f"  Case prevalence: {phenotype_df['case_control'].mean():.1%}")
    print(f"\nSaved to {args.out_dir}/genotypes.csv, effect_sizes.csv, phenotype.csv")


if __name__ == "__main__":
    main()
