"""
calculate_prs.py

Calculates a Polygenic Risk Score for each individual as the weighted
sum of their genotypes, weighted by per-SNP effect sizes.

    PRS_i = sum over all SNPs j of (genotype_ij * effect_weight_j)

This is the standard "clumping + thresholding" or simple weighted-sum
approach used by tools like PRSice-2 and PLINK --score, simplified here
for clarity. In a full implementation you would also need to:
  - Match SNPs between the target genotype data and the effect size
    (scoring) file by ID, harmonizing strand/allele coding
  - Apply LD clumping/thresholding if using raw GWAS summary stats
    rather than a pre-clumped scoring file
  - Standardize the resulting PRS (z-score) for interpretability

This script assumes genotypes and effect sizes are already aligned on
SNP ID -- the QC step (qc_checks.py) is what would normally flag a
mismatch before scoring proceeds.
"""

import argparse
import pandas as pd
import numpy as np


def calculate_prs(genotypes_df: pd.DataFrame, effect_sizes_df: pd.DataFrame) -> pd.DataFrame:
    """
    genotypes_df: must contain 'individual_id' plus one column per SNP
    effect_sizes_df: must contain 'snp_id' and 'effect_weight'

    Returns a DataFrame with individual_id, raw PRS, and standardized PRS.
    """
    snp_cols = [c for c in genotypes_df.columns if c != "individual_id"]

    # Align effect sizes to the SNP columns present in the genotype data
    weights = effect_sizes_df.set_index("snp_id")["effect_weight"]
    weights_aligned = weights.reindex(snp_cols)

    missing = weights_aligned.isna().sum()
    if missing > 0:
        print(f"Warning: {missing} SNPs in genotype data have no matching "
              f"effect weight and will be treated as zero-effect.")
        weights_aligned = weights_aligned.fillna(0)

    genotype_matrix = genotypes_df[snp_cols].to_numpy()
    raw_prs = genotype_matrix @ weights_aligned.to_numpy()

    prs_df = pd.DataFrame({
        "individual_id": genotypes_df["individual_id"],
        "prs_raw": raw_prs,
    })
    prs_df["prs_zscore"] = (prs_df["prs_raw"] - prs_df["prs_raw"].mean()) / prs_df["prs_raw"].std()

    return prs_df


def main():
    parser = argparse.ArgumentParser(description="Calculate PRS from genotype + effect size data")
    parser.add_argument("--genotypes", type=str, default="data/genotypes.csv")
    parser.add_argument("--effect_sizes", type=str, default="data/effect_sizes.csv")
    parser.add_argument("--out", type=str, default="results/prs_scores.csv")
    args = parser.parse_args()

    genotypes_df = pd.read_csv(args.genotypes)
    effect_sizes_df = pd.read_csv(args.effect_sizes)

    prs_df = calculate_prs(genotypes_df, effect_sizes_df)
    prs_df.to_csv(args.out, index=False)

    print(f"Calculated PRS for {len(prs_df)} individuals")
    print(f"  Raw PRS range: [{prs_df['prs_raw'].min():.3f}, {prs_df['prs_raw'].max():.3f}]")
    print(f"  Mean (raw): {prs_df['prs_raw'].mean():.3f}, SD: {prs_df['prs_raw'].std():.3f}")
    print(f"\nSaved to {args.out}")


if __name__ == "__main__":
    main()
