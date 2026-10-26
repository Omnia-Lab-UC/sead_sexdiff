#!/usr/bin/env python3
"""Move superseded output directories under ``<results>/superseded/<change>/``.

A MOTCO pin move that changes the measurement re-issues every MOTCO-dependent
output (spec ``motco-release-migration``). Before the reruns, the superseded
directories are moved, never deleted: each keeps its path relative to the
results root, every file is hashed before and after the move (they must be
equal), and ``<results>/superseded/<change>/MOVED.json`` records old path, new
path and sha256 per file. Nothing under ``superseded/`` is an input of any
step.

Usage (from the repository root; standard library only, so it runs with any
python3, e.g. on the cluster login node)::

    python3 scripts/move_superseded.py --change mtg-rank-floor-primary \\
        motco_analysis motco_rank_sensitivity ... microglia/motco_analysis ...
    python3 scripts/move_superseded.py --change mtg-astrocyte-latent-adjustment-sensitivity \\
        motco_analysis_0.8.0=motco_analysis

Each argument is ``old[=new]``, relative to ``--results`` (default
``results``); ``new`` defaults to ``old``. A directory that does not exist is
skipped and listed (a machine may lack some outputs); an existing target stops
the script before anything moves. ``--dry-run`` hashes and prints the plan.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_tree(root: Path) -> dict[str, str]:
    """Relative path → sha256 of every regular file under ``root``; symlinks are an error."""
    out = {}
    for dirpath, dirnames, filenames in os.walk(root):
        for name in sorted(dirnames + filenames):
            p = Path(dirpath) / name
            if p.is_symlink():
                raise SystemExit(f"{p} is a symlink; move it by hand, nothing was moved")
        for name in sorted(filenames):
            p = Path(dirpath) / name
            out[str(p.relative_to(root))] = sha256_file(p)
    return dict(sorted(out.items()))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--change", required=True, help="Name of the change that supersedes the outputs")
    parser.add_argument("--results", default="results", help="Results root (default: results)")
    parser.add_argument("--dry-run", action="store_true", help="Hash and print the plan; move nothing")
    parser.add_argument("dirs", nargs="+", help="old[=new], relative to --results")
    args = parser.parse_args(argv)

    results = Path(args.results)
    target_root = results / "superseded" / args.change
    plan, skipped = [], []
    for spec in args.dirs:
        old_rel, _, new_rel = spec.partition("=")
        new_rel = new_rel or old_rel
        old, new = results / old_rel, target_root / new_rel
        if Path(old_rel).is_absolute() or ".." in Path(old_rel).parts or "superseded" in Path(old_rel).parts:
            raise SystemExit(f"{spec}: give a path relative to {results} outside superseded/")
        if not old.exists():
            skipped.append(str(old))
            continue
        if not old.is_dir() or old.is_symlink():
            raise SystemExit(f"{old} is not a directory; nothing was moved")
        if new.exists():
            raise SystemExit(f"{new} exists; nothing was moved")
        plan.append((old, new))
    if not plan:
        print(f"[move-superseded] nothing to move (absent: {skipped})")
        return 0
    if os.stat(results).st_dev != os.stat(plan[0][0]).st_dev:
        raise SystemExit("the directories are not on the results root's filesystem; nothing was moved")

    before = {}
    for old, new in plan:
        before[old] = hash_tree(old)
        print(f"[move-superseded] {old} -> {new} ({len(before[old])} files)")
    for path in skipped:
        print(f"[move-superseded] absent, skipped: {path}")
    if args.dry_run:
        return 0

    manifest_path = target_root / "MOVED.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {
        "change": args.change,
        "note": (
            "superseded outputs moved here, never deleted; each keeps its path relative to the results root; "
            "never an input of any step"
        ),
        "moves": [],
    }
    for old, new in plan:
        new.parent.mkdir(parents=True, exist_ok=True)
        os.rename(old, new)
        after = hash_tree(new)
        if after != before[old]:
            changed = sorted(set(after.items()) ^ set(before[old].items()))
            raise SystemExit(f"{new}: hashes differ after the move ({changed[:5]} ...); investigate before going on")
        manifest["moves"].append({
            "old": str(old), "new": str(new), "moved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "host": os.uname().nodename, "n_files": len(after),
            "files": [{"path": rel, "sha256": h} for rel, h in after.items()],
        })
        print(f"[move-superseded] moved {old} -> {new}; {len(after)} files hash as before")
    manifest["absent"] = sorted(set(manifest.get("absent", [])) | set(skipped))
    target_root.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[move-superseded] wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
