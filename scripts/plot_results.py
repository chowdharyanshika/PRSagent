"""
plot_results.py

Generates the standard plots used to communicate PRS results:
  1. PRS distribution by case/control status (density plot)
  2. Disease prevalence by PRS decile (bar chart)
  3. ROC curve

These mirror the figures typically shown in PRS publications and are
far more interpretable to a non-technical audience than a table of
statistics alone.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, roc_auc_score

sns.set_theme(style="whitegrid")


def plot_prs_distribution(merged: pd.DataFrame, ax):
    sns.kdeplot(data=merged, x="prs_zscore", hue="case_control",
                fill=True, alpha=0.4, ax=ax, palette=["#4C72B0", "#C44E52"])
    ax.set_xlabel("PRS (z-score)")
    ax.set_ylabel("Density")
    ax.set_title("PRS Distribution by Case/Control Status")
    ax.legend(title="Status", labels=["Case", "Control"])


def plot_risk_by_decile(merged: pd.DataFrame, ax):
    merged["prs_decile"] = pd.qcut(merged["prs_zscore"], 10, labels=False) + 1
    risk_by_decile = merged.groupby("prs_decile")["case_control"].mean()
    ax.bar(risk_by_decile.index, risk_by_decile.values, color="#55A868")
    ax.set_xlabel("PRS Decile (1 = lowest risk, 10 = highest risk)")
    ax.set_ylabel("Disease Prevalence")
    ax.set_title("Disease Prevalence by PRS Decile")
    ax.set_xticks(range(1, 11))


def plot_roc(merged: pd.DataFrame, ax):
    from sklearn.linear_model import LogisticRegression
    X = merged[["prs_zscore"]].to_numpy()
    y = merged["case_control"].to_numpy()
    model = LogisticRegression().fit(X, y)
    probs = model.predict_proba(X)[:, 1]
    fpr, tpr, _ = roc_curve(y, probs)
    auc = roc_auc_score(y, probs)
    ax.plot(fpr, tpr, color="#4C72B0", label=f"PRS (AUC = {auc:.3f})")
    ax.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Chance")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curve")
    ax.legend()


def main():
    parser = argparse.ArgumentParser(description="Generate PRS result visualizations")
    parser.add_argument("--prs", type=str, default="results/prs_scores.csv")
    parser.add_argument("--phenotype", type=str, default="data/phenotype.csv")
    parser.add_argument("--out", type=str, default="results/prs_plots.png")
    args = parser.parse_args()

    prs_df = pd.read_csv(args.prs)
    phenotype_df = pd.read_csv(args.phenotype)
    merged = prs_df.merge(phenotype_df, on="individual_id")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    plot_prs_distribution(merged, axes[0])
    plot_risk_by_decile(merged, axes[1])
    plot_roc(merged, axes[2])

    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    print(f"Saved plots to {args.out}")


if __name__ == "__main__":
    main()
