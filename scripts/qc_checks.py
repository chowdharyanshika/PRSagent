"""
qc_checks.py

Runs basic quality control checks that should ALWAYS be performed before
trusting a calculated PRS. This is the step most tutorials skip, and
exactly the kind of judgment call an autonomous agent should be making
rather than blindly reporting a number.

Checks implemented:
  1. SNP overlap: what fraction of SNPs in the scoring/effect file were
     actually found in the target genotype data? Low overlap means the
     PRS is built on incomplete information and should not be trusted.
  2. Sample size: is the cohort large enough for the PRS-phenotype
     association to be statistically meaningful?
  3. Missingness: are there individuals with excessive missing genotype
     data who should be excluded?
  4. Variance sanity check: does the PRS have non-zero variance? (a
     constant PRS usually indicates a SNP-matching bug upstream)

In a real-world pipeline you would also check ancestry/population match
between the GWAS discovery cohort and the target sample (a major source
of PRS bias and a key limitation to disclose in any interpretation) --
flagged here as a placeholder since it requires ancestry reference data
not included in this simplified simulation.
"""

import argparse
import pandas as pd
import numpy as np


def run_qc(
    genotypes_df: pd.DataFrame,
    effect_sizes_df: pd.DataFrame,
    prs_df: pd.DataFrame,
    min_snp_overlap: float = 0.8,
    min_sample_size: int = 500,
) -> dict:
    snp_cols = [c for c in genotypes_df.columns if c != "individual_id"]
    effect_snp_ids = set(effect_sizes_df["snp_id"])
    genotype_snp_ids = set(snp_cols)

    overlap = effect_snp_ids & genotype_snp_ids
    snp_overlap_fraction = len(overlap) / len(effect_snp_ids) if effect_snp_ids else 0.0

    n_individuals = len(genotypes_df)

    missingness_per_individual = genotypes_df[snp_cols].isna().mean(axis=1)
    high_missingness_count = (missingness_per_individual > 0.05).sum()

    prs_variance = prs_df["prs_raw"].var()

    checks = {
        "snp_overlap_fraction": round(snp_overlap_fraction, 4),
        "snp_overlap_pass": snp_overlap_fraction >= min_snp_overlap,
        "sample_size": n_individuals,
        "sample_size_pass": n_individuals >= min_sample_size,
        "individuals_with_high_missingness": int(high_missingness_count),
        "missingness_pass": high_missingness_count == 0,
        "prs_variance": round(float(prs_variance), 6),
        "prs_variance_pass": prs_variance > 1e-8,
        "ancestry_match_checked": False,  # placeholder -- see module docstring
    }

    checks["overall_pass"] = all([
        checks["snp_overlap_pass"],
        checks["sample_size_pass"],
        checks["missingness_pass"],
        checks["prs_variance_pass"],
    ])

    return checks


def print_qc_report(checks: dict):
    print("=" * 50)
    print("PRS QUALITY CONTROL REPORT")
    print("=" * 50)

    print(f"\n1. SNP overlap: {checks['snp_overlap_fraction']:.1%} "
          f"{'PASS' if checks['snp_overlap_pass'] else 'FAIL (below threshold)'}")

    print(f"2. Sample size: {checks['sample_size']} individuals "
          f"{'PASS' if checks['sample_size_pass'] else 'FAIL (too small for reliable association testing)'}")

    print(f"3. High-missingness individuals: {checks['individuals_with_high_missingness']} "
          f"{'PASS' if checks['missingness_pass'] else 'FAIL (consider excluding these individuals)'}")

    print(f"4. PRS variance check: {checks['prs_variance']:.6f} "
          f"{'PASS' if checks['prs_variance_pass'] else 'FAIL (PRS has no variance -- check SNP matching)'}")

    print(f"5. Ancestry/population match: NOT CHECKED "
          f"(requires reference population data -- see module docstring)")

    print("\n" + "=" * 50)
    verdict = "PASS - results can be reported with standard caveats" \
        if checks["overall_pass"] else \
        "FAIL - do not report results without addressing the failed check(s) above"
    print(f"OVERALL: {verdict}")
    print("=" * 50)


def main():
    parser = argparse.ArgumentParser(description="Run QC checks on a calculated PRS")
    parser.add_argument("--genotypes", type=str, default="data/genotypes.csv")
    parser.add_argument("--effect_sizes", type=str, default="data/effect_sizes.csv")
    parser.add_argument("--prs", type=str, default="results/prs_scores.csv")
    args = parser.parse_args()

    genotypes_df = pd.read_csv(args.genotypes)
    effect_sizes_df = pd.read_csv(args.effect_sizes)
    prs_df = pd.read_csv(args.prs)

    checks = run_qc(genotypes_df, effect_sizes_df, prs_df)
    print_qc_report(checks)


if __name__ == "__main__":
    main()
