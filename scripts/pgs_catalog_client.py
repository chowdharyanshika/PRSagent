"""
pgs_catalog_client.py

Searches the PGS Catalog (https://www.pgscatalog.org) for existing,
published Polygenic Score models matching a given trait.

The PGS Catalog is a curated, free database of PRS models from
peer-reviewed publications. Checking here FIRST (before building a PRS
from raw GWAS summary statistics) is standard practice in the field --
it avoids reinventing a score that has already been validated, and it's
exactly the kind of decision an agent should make autonomously in
Phase 2.

API docs: https://www.pgscatalog.org/rest/

NOTE on this script's reliability:
The exact JSON response shape below (`results`, `id`, `variants_number`,
etc.) is based on documented PGS Catalog REST conventions, but it has
NOT been verified against a live network call in this environment (no
outbound internet access here). Before relying on this in a demo or
interview, run it with real network access and:
  1. Print the raw `trait_results` / `scores_results` dicts first
  2. Confirm the key names match what's parsed below
  3. Adjust field names if the live API differs
Treat the parsing logic as a first draft to verify, not ground truth --
good practice for any external API integration.
"""

import argparse
import json
import sys
import requests

BASE_URL = "https://www.pgscatalog.org/rest"


def search_trait(trait_term: str, debug: bool = False) -> dict:
    """
    Search the PGS Catalog for an EFO trait matching the search term,
    then return all PGS (PRS) models associated with that trait.

    Set debug=True the first time you run this with real network access
    -- it dumps the raw JSON so you can confirm/fix the field names
    parsed below against the live schema.
    """
    # Step 1: find the matching trait entry (EFO ontology term)
    trait_search_url = f"{BASE_URL}/trait/search"
    resp = requests.get(trait_search_url, params={"term": trait_term}, timeout=30)
    resp.raise_for_status()
    trait_results = resp.json()

    if debug:
        print("--- RAW trait search response ---")
        print(json.dumps(trait_results, indent=2)[:2000])
        print("--- end raw response ---\n")

    if not trait_results.get("results"):
        return {"trait_term": trait_term, "matched_trait": None, "scores": []}

    # Take the first (best) matching trait
    matched_trait = trait_results["results"][0]
    efo_id = matched_trait["id"]

    # Step 2: fetch PGS scores associated with that trait
    scores_url = f"{BASE_URL}/score/search"
    resp = requests.get(scores_url, params={"trait_id": efo_id}, timeout=30)
    resp.raise_for_status()
    scores_results = resp.json()

    if debug:
        print("--- RAW score search response ---")
        print(json.dumps(scores_results, indent=2)[:2000])
        print("--- end raw response ---\n")

    scores = []
    for s in scores_results.get("results", []):
        scores.append({
            "pgs_id": s.get("id"),
            "name": s.get("name"),
            "num_variants": s.get("variants_number"),
            "publication": s.get("publication", {}).get("title"),
            "ftp_scoring_file": s.get("ftp_scoring_file"),
        })

    return {
        "trait_term": trait_term,
        "matched_trait": {
            "efo_id": efo_id,
            "label": matched_trait.get("label"),
        },
        "n_scores_found": len(scores),
        "scores": scores,
    }


def rank_by_snp_coverage(scores: list) -> list:
    """
    Sort candidate PRS models by number of variants (more variants is
    generally -- though not always -- a proxy for a more comprehensive,
    well-powered score). This is a simple heuristic an agent can use
    when choosing between multiple available models.
    """
    return sorted(
        [s for s in scores if s.get("num_variants")],
        key=lambda s: s["num_variants"],
        reverse=True,
    )


def main():
    parser = argparse.ArgumentParser(description="Search PGS Catalog for a trait")
    parser.add_argument(
        "--trait", type=str, required=True,
        help="Trait name to search for, e.g. 'type 2 diabetes'"
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="Print raw API responses (use on first real run to verify schema)"
    )
    args = parser.parse_args()

    print(f"Searching PGS Catalog for trait: '{args.trait}'...\n")

    try:
        result = search_trait(args.trait, debug=args.debug)
    except requests.exceptions.RequestException as e:
        print(f"Error contacting PGS Catalog API: {e}", file=sys.stderr)
        sys.exit(1)

    if result["matched_trait"] is None:
        print(f"No matching trait found for '{args.trait}'. Try a different term.")
        sys.exit(0)

    print(f"Matched trait: {result['matched_trait']['label']} "
          f"(EFO: {result['matched_trait']['efo_id']})")
    print(f"Found {result['n_scores_found']} existing PRS model(s):\n")

    ranked = rank_by_snp_coverage(result["scores"])
    for s in ranked[:5]:
        print(f"  {s['pgs_id']}: {s['name']}")
        print(f"    Variants: {s['num_variants']}")
        print(f"    Publication: {s['publication']}")
        print()

    # Save full result for downstream use
    out_path = "data/pgs_search_result.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Full results saved to {out_path}")


if __name__ == "__main__":
    main()
