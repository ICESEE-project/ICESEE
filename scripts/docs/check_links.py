"""Check repository-relative Markdown destinations without importing ICESEE.

This validates the handbook and the principal documentation entry points.
It checks local files and source line anchors; remote URLs are not fetched.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root.resolve()
    paths = [root / "README.md", root / "docs/model-adapter.md",
             root / "docs/execution-mode-3-design.md"]
    for folder in ["docs/handbook", "docs/architecture", "paper/production-report"]:
        paths.extend(sorted((root / folder).glob("*.md")))
    failures = []
    checked = 0
    for page in paths:
        if not page.is_file():
            failures.append({"page": str(page.relative_to(root)), "reason": "missing page"})
            continue
        for match in re.finditer(r"\]\(([^\s)]+)\)", page.read_text()):
            raw = match.group(1).strip("<>")
            link = urlsplit(raw)
            if link.scheme or link.netloc or not link.path:
                continue
            destination = (page.parent / unquote(link.path)).resolve()
            checked += 1
            reason = None
            if not destination.exists():
                reason = "missing destination"
            elif re.fullmatch(r"L\d+", link.fragment) and destination.is_file():
                line = int(link.fragment[1:])
                if not 1 <= line <= len(destination.read_text().splitlines()):
                    reason = "source line anchor out of range"
            if reason:
                failures.append({"page": str(page.relative_to(root)), "link": raw, "reason": reason})
    print(json.dumps({"pages": len(paths), "local_links": checked, "failures": failures}, indent=2))
    raise SystemExit(bool(failures))


if __name__ == "__main__":
    main()
