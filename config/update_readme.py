#!/usr/bin/env python3
"""
update_flags_readme.py

Refreshes the ICESEE flags section in config/README.md using extract_flags/generate_flags_markdown.
- Detects repo root via git (fallback: walk up until 'src' exists).
- Replaces existing flags block (case-insensitive header variants) or the anchored block.
- Writes atomically; creates config/ if missing.
- Stages the file if running inside a git repo.
- Idempotent: no rewrite when content is unchanged.

Usage:
  python update_flags_readme.py [path/to/script_to_scan.py]

Defaults to: <repo_root>/config/_utility_imports.py
@author: Brian Kyanjo
@date: 2025-09-15
"""

from __future__ import annotations

import os
import re
import sys
import math
import tempfile
import subprocess
from pathlib import Path
from typing import Optional

from extract_flags import extract_flags, generate_flags_markdown


# ---- Constants ----------------------------------------------------------------

SECTION_CANONICAL_TITLE = "## All Main Flags Used in ICESEE"

# Header-based block match (fallback when anchors aren't present).
# Matches any H2 exactly for the flags title regardless of casing on "Used/used",
# and replaces until next H2 or EOF.
FLAGS_BLOCK_REGEX = re.compile(
    r"(?ims)^\s*##\s+All\s+Main\s+Flags\s+(?:Used|used)\s+in\s+ICESEE\s*\n.*?(?=^\s*##\s|\Z)"
)

# Preferred: anchored block (stable and explicit).
ANCHOR_BLOCK_REGEX = re.compile(
    r"(?s)<!--\s*BEGIN:\s*ICESEE-FLAGS\s*-->.*?<!--\s*END:\s*ICESEE-FLAGS\s*-->"
)

# Destination path (relative to repo root)
README_REL_PATH = Path("config/README.md")

# Default script to scan (relative to repo root)
DEFAULT_SCRIPT_REL_PATH = Path("config/_utility_imports.py")


# ---- Subprocess helpers -------------------------------------------------------

def run(cmd: list[str], cwd: Optional[Path] = None, check: bool = False) -> subprocess.CompletedProcess:
    """Run a command; capture stdout/stderr. Set check=True to raise on non-zero."""
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        check=check,
    )


def get_repo_root_via_git(start: Path) -> Optional[Path]:
    """Return git toplevel directory or None if not a git repo."""
    try:
        cp = run(["git", "rev-parse", "--show-toplevel"], cwd=start)
        if cp.returncode == 0 and cp.stdout.strip():
            return Path(cp.stdout.strip())
    except Exception:
        pass
    return None


# ---- Filesystem helpers -------------------------------------------------------

def get_project_root() -> Path:
    """
    Prefer the Git toplevel; fall back to walking up until a 'src' directory is found.
    As a last resort, return the script directory.
    """
    here = Path(__file__).resolve().parent
    root = get_repo_root_via_git(here)
    if root:
        return root

    cur = here
    while cur != cur.root:
        if (cur / "src").exists():
            return cur
        cur = cur.parent
    return here


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return "# ICESEE FLAGS README\n\n"


def write_atomic(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8", dir=str(path.parent)) as tmp:
        tmp.write(content)
        tmp_path = Path(tmp.name)
    os.replace(tmp_path, path)  # atomic on POSIX


# ---- Content normalization ----------------------------------------------------

def strip_existing_flags_header(block: str) -> str:
    """Remove any 'All Main Flags Used/used in ICESEE' header at the start of a block."""
    return re.sub(
        r"(?i)^\s*##\s+All\s+Main\s+Flags\s+(?:Used|used)\s+in\s+ICESEE\s*\n+",
        "",
        block.strip()
    ).strip()


def canonicalize_flags_doc(flags_doc: str) -> str:
    """
    Ensure generated section:
    - has canonical title,
    - contains no duplicate header inside,
    - ends with a single newline,
    - collapses excessive blank lines.
    """
    body = strip_existing_flags_header(flags_doc)
    section = f"{SECTION_CANONICAL_TITLE}\n\n{body.strip()}\n"
    section = re.sub(r"\n{3,}", "\n\n", section)
    return section


def wrap_with_anchors(section: str) -> str:
    """Wrap a section with BEGIN/END anchors (preferred replacement target)."""
    return f"<!-- BEGIN: ICESEE-FLAGS -->\n{section}<!-- END: ICESEE-FLAGS -->"


def update_section(existing: str, new_section: str) -> str:
    """
    Replace the flags section in `existing` with `new_section`:
    1) If anchors exist: replace the anchored block.
    2) Else if header-based block exists: replace all matches.
    3) Else append once at the end (with an extra blank line if needed).
    Normalize redundant blank lines and ensure trailing newline.
    """
    anchored = wrap_with_anchors(new_section)

    if ANCHOR_BLOCK_REGEX.search(existing):
        updated = ANCHOR_BLOCK_REGEX.sub(anchored, existing)
    elif FLAGS_BLOCK_REGEX.search(existing):
        updated = FLAGS_BLOCK_REGEX.sub(anchored, existing)  # replace ALL header matches
    else:
        base = existing.rstrip()
        updated = (base + "\n\n" if base else "") + anchored

    updated = re.sub(r"\n{3,}", "\n\n", updated).rstrip() + "\n"
    return updated


# ---- Git staging --------------------------------------------------------------

def stage_if_git_repo(repo_root: Path, path: Path) -> None:
    """Stage the file if inside a git repo; ignore failures."""
    if get_repo_root_via_git(repo_root):
        try:
            run(["git", "add", str(path)], cwd=repo_root)
        except Exception:
            pass


# ---- Main update logic --------------------------------------------------------

def update_readme(script_path: Path) -> bool:
    """
    Extract flags from `script_path`, generate canonical section, and update README.md
    iff changes occur. Returns True if a write happened, else False.
    """
    repo_root = get_project_root()
    readme_path = repo_root / README_REL_PATH

    flags = extract_flags(str(script_path))
    raw_flags_doc = generate_flags_markdown(flags)
    flags_doc = canonicalize_flags_doc(raw_flags_doc)

    existing = read_text(readme_path)
    updated = update_section(existing, flags_doc)

    if existing == updated:
        print("[flags-readme] No changes detected; README unchanged.")
        return False

    write_atomic(readme_path, updated)
    stage_if_git_repo(repo_root, readme_path)
    rel = readme_path.relative_to(repo_root) if readme_path.is_absolute() else readme_path
    print(f"[flags-readme] Updated {rel}")
    return True


def main(argv: list[str]) -> int:
    repo_root = get_project_root()
    default_script = (repo_root / DEFAULT_SCRIPT_REL_PATH).resolve()
    script_path = Path(argv[1]).resolve() if len(argv) > 1 else default_script

    if not script_path.exists():
        print(f"[flags-readme] ERROR: script not found: {script_path}", file=sys.stderr)
        return 2

    print(f"[flags-readme] Updating README with flags from {script_path}")
    try:
        updated = update_readme(script_path)
    except Exception as e:
        print(f"[flags-readme] ERROR: {e}", file=sys.stderr)
        return 1

    # success (even if no changes were needed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))