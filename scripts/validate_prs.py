"""
validate_prs.py

Tests whether the calculated PRS actually predicts the phenotype, and
reports results the way they should be reported in practice: with an
effect size AND a measure of how much variance is explained (since a
statistically significant PRS can still have very limited individual-
level predictive power -- a frequently misunderstood point about PRS
that's worth demonstrating explicit understanding of).

Metrics reported:
  - Odds ratio per SD increase in PRS (from logistic regression)
  - AUC (discrimination ability)
  - Nagelkerke's pseudo-R^2 (variance explained, the standard PRS metric)
  - Comparison across PRS quantiles (e.g. top 10% vs bottom 10% risk)
"""

import argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


def nagelkerke_r2(model: LogisticRegression, X: np.ndarray, y: np.ndarray) -> float:
    """
    Nagelkerke's pseudo-R^2: the standard way variance explained is
    reported for PRS in the literature, since ordinary R^2 doesn't
    apply directly to logistic regression.
    """
    n = len(y)
    log_likelihood_full = np.sum(
        y * np.log(model.predict_proba(X)[:, 1] + 1e-10) +
        (1 - y) * np.log(1 - model.predict_proba(X)[:, 1] + 1e-10)
    )
    p_null = y.mean()
    log_likelihood_null = np.sum(
        y * np.log(p_null) + (1 - y) * np.log(1 - p_null)
    )
    cox_snell_r2 = 1 - np.exp((2 / n) * (log_likelihood_null - log_likelihood_full))
    max_r2 = 1 - np.exp((2 / n) * log_likelihood_null)
    return cox_snell_r2 / max_r2


def validate_prs(prs_df: pd.DataFrame, phenotype_df: pd.DataFrame) -> dict:
    merged = prs_df.merge(phenotype_df, on="individual_id")

    X = merged[["prs_zscore"]].to_numpy()
    y = merged["case_control"].to_numpy()

    model = LogisticRegression()
    model.fit(X, y)

    odds_ratio_per_sd = np.exp(model.coef_[0][0])
    predicted_probs = model.predict_proba(X)[:, 1]
    auc = roc_auc_score(y, predicted_probs)
    r2 = nagelkerke_r2(model, X, y)

    # Risk by PRS decile -- intuitive way to communicate results
    merged["prs_decile"] = pd.qcut(merged["prs_zscore"], 10, labels=False) + 1
    risk_by_decile = merged.groupby("prs_decile")["case_control"].mean()
    top_vs_bottom_decile_ratio = (
        risk_by_decile.iloc[-1] / risk_by_decile.iloc[0]
        if risk_by_decile.iloc[0] > 0 else float("inf")
    )

    return {
        "n_individuals": len(merged),
        "odds_ratio_per_sd": round(float(odds_ratio_per_sd), 3),
        "auc": round(float(auc), 3),
        "nagelkerke_r2": round(float(r2), 4),
        "risk_top_decile": round(float(risk_by_decile.iloc[-1]), 3),
        "risk_bottom_decile": round(float(risk_by_decile.iloc[0]), 3),
        "top_vs_bottom_decile_ratio": round(float(top_vs_bottom_decile_ratio), 2),
        "risk_by_decile": risk_by_decile.to_dict(),
    }


def print_interpretation(results: dict):
    print("=" * 55)
    print("PRS VALIDATION & INTERPRETATION")
    print("=" * 55)
    print(f"\nSample size: {results['n_individuals']} individuals")
    print(f"\nOdds ratio per 1 SD increase in PRS: {results['odds_ratio_per_sd']}x")
    print(f"AUC (discrimination): {results['auc']}")
    print(f"Variance explained (Nagelkerke R2): {results['nagelkerke_r2']:.1%}")
    print(f"\nRisk comparison:")
    print(f"  Bottom 10% of PRS distribution: {results['risk_bottom_decile']:.1%} disease prevalence")
    print(f"  Top 10% of PRS distribution:    {results['risk_top_decile']:.1%} disease prevalence")
    print(f"  Relative risk (top vs bottom decile): {results['top_vs_bottom_decile_ratio']}x")

    print("\n" + "-" * 55)
    print("HONEST INTERPRETATION:")
    r2_pct = results["nagelkerke_r2"] * 100
    print(f"  This PRS explains approximately {r2_pct:.1f}% of the variance")
    print(f"  in disease status in this simulated cohort. An AUC of "
          f"{results['auc']} indicates")
    if results["auc"] < 0.6:
        print("  POOR discrimination -- this PRS has limited individual-level")
        print("  predictive value, though it may still be useful for")
        print("  population-level risk stratification research.")
    elif results["auc"] < 0.75:
        print("  MODEST discrimination -- typical of PRS for many common")
        print("  complex diseases. Useful for risk stratification at the")
        print("  population level, but should NOT be used alone for")
        print("  individual clinical decisions.")
    else:
        print("  STRONG discrimination for a PRS, though still likely")
        print("  insufficient as a standalone diagnostic tool.")
    print("-" * 55)


def main():
    parser = argparse.ArgumentParser(description="Validate PRS against phenotype data")
    parser.add_argument("--prs", type=str, default="results/prs_scores.csv")
    parser.add_argument("--phenotype", type=str, default="data/phenotype.csv")
    parser.add_argument("--out", type=str, default="results/validation_report.csv")
    args = parser.parse_args()

    prs_df = pd.read_csv(args.prs)
    phenotype_df = pd.read_csv(args.phenotype)

    results = validate_prs(prs_df, phenotype_df)
    print_interpretation(results)

    # Save the decile breakdown for plotting later
    decile_df = pd.DataFrame(
        list(results["risk_by_decile"].items()),
        columns=["prs_decile", "disease_prevalence"]
    )
    decile_df.to_csv(args.out, index=False)
    print(f"\nDecile breakdown saved to {args.out}")


if __name__ == "__main__":
    main()
