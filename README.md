# PRS Agent: Autonomous Polygenic Risk Score Analysis

An agentic AI system that calculates, validates, and interprets Polygenic
Risk Scores (PRS) using public genomic data sources — combining classical
statistical genetics with modern LLM-based tool orchestration.

## Why this project

Polygenic Risk Scores aggregate the effect of many common genetic variants
to estimate an individual's genetic predisposition to a disease or trait.
Calculating a *trustworthy* PRS is not a one-step process — it requires:

1. Finding the right GWAS summary statistics or existing validated PRS weights
2. Checking SNP overlap between the GWAS source and the target genotype data
3. Accounting for linkage disequilibrium (LD) between variants
4. Checking for population/ancestry mismatch (a major source of PRS bias)
5. Validating the score against a phenotype, if available
6. Interpreting the result in a way that is statistically honest
   (PRS explain *some* variance, not a diagnosis)

This project builds an **agent** that makes these decisions autonomously,
rather than a fixed script — given a trait name and a genotype file, it
decides which data source to use, runs the appropriate quality checks,
and only then produces an interpretation.

## Project phases

- **Phase 1 (this stage):** Core, deterministic pipeline. Search PGS
  Catalog for existing PRS models, calculate a score on (simulated)
  genotype data, run basic QC, and report results. No agent yet — this
  proves each tool works correctly in isolation.
- **Phase 2:** Wrap each step as a LangGraph tool and let an LLM agent
  decide the workflow (e.g. compare multiple candidate PRS models, choose
  the one with best SNP overlap, decide whether to fall back to raw GWAS
  summary statistics).
- **Phase 3:** Add a logistic regression layer testing whether the PRS
  predicts a simulated phenotype, plus ancestry-mismatch warnings.
- **Phase 4:** Streamlit front-end + polished documentation.

## Data sources (all public)

| Source | Purpose | Link |
|---|---|---|
| PGS Catalog | Pre-computed, published PRS weight files | https://www.pgscatalog.org |
| GWAS Catalog | Raw summary statistics if no PGS exists | https://www.ebi.ac.uk/gwas |
| 1000 Genomes Project | Reference genotype data for testing | https://www.internationalgenome.org |

No real patient data is used. For development and demonstration, genotype
data is **simulated** (`scripts/simulate_genotypes.py`) to avoid any
consent or privacy ambiguity — this mirrors best practice for methods
development in statistical genetics.

## Repo structure

```
prs-agent/
├── data/                   # downloaded/simulated data (gitignored)
├── scripts/
│   ├── pgs_catalog_client.py   # search & download PRS models from PGS Catalog
│   ├── simulate_genotypes.py   # generate realistic toy genotype + phenotype data
│   ├── calculate_prs.py        # core PRS scoring logic
│   └── qc_checks.py            # SNP overlap, sample size, basic validation
├── notebooks/
│   └── 01_pipeline_walkthrough.ipynb
├── results/                 # output scores, QC reports
├── docs/
│   └── methodology.md       # explains the statistical genetics behind PRS
├── requirements.txt
└── README.md
```

## Quick start

```bash
pip install -r requirements.txt --break-system-packages

# 1. Generate toy data (simulated SNPs + effect sizes + phenotype)
python scripts/simulate_genotypes.py

# 2. Search PGS Catalog for a real trait's existing PRS models
python scripts/pgs_catalog_client.py --trait "type 2 diabetes"

# 3. Calculate a PRS on the simulated data using simulated weights
python scripts/calculate_prs.py

# 4. Run QC checks on the result
python scripts/qc_checks.py
```

