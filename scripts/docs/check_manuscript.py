"""Check local citation and cross-reference identities in the technical report.

This checks the embedded bibliography rather than fetching literature or
substituting for TeX compilation and visual review.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manuscript", type=Path,
        default=Path(__file__).resolve().parents[2] / "paper/production-report/main.tex",
    )
    args = parser.parse_args()
    source = args.manuscript.read_text()
    citation_groups = re.findall(
        r"\\cite(?:p|t|author|year)?\*?(?:\[[^\]]*\])*\{([^}]+)\}", source
    )
    cited = {key.strip() for group in citation_groups for key in group.split(",")}
    bibliography = re.findall(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", source)
    labels = re.findall(r"\\label\{([^}]+)\}", source)
    references = set(re.findall(r"\\(?:ref|eqref|pageref|autoref)\{([^}]+)\}", source))
    duplicates = lambda items: sorted(k for k, n in Counter(items).items() if n > 1)
    result = {
        "citation_groups": len(citation_groups),
        "distinct_cited_sources": len(cited),
        "bibliography_entries": len(bibliography),
        "missing_citation_keys": sorted(cited - set(bibliography)),
        "uncited_bibliography_keys": sorted(set(bibliography) - cited),
        "duplicate_bibliography_keys": duplicates(bibliography),
        "missing_reference_labels": sorted(references - set(labels)),
        "duplicate_labels": duplicates(labels),
    }
    print(json.dumps(result, indent=2))
    errors = [result[k] for k in (
        "missing_citation_keys", "duplicate_bibliography_keys",
        "missing_reference_labels", "duplicate_labels",
    )]
    raise SystemExit(any(errors))


if __name__ == "__main__":
    main()
