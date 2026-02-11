"""
Download all experimentally resolved Cas13-family PDB structures from RCSB.

Strategy:
  1. Start with a curated baseline of well-known Cas13 PDB IDs.
  2. Hit the RCSB Search API (full-text "Cas13") to discover any new deposits.
  3. Download every PDB file to  data/structure_pipeline/references/
  4. Write  cas13_reference_catalog.json  with PDB ID, title, subtype, resolution.

Usage:
    python visualization/download_cas13_references.py          # full run
    python visualization/download_cas13_references.py --dry-run # list only
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

# ---------------------------------------------------------------------------
# Curated baseline – guaranteed Cas13 structures (as of early 2026)
# Format: PDB_ID -> (subtype_label, short_description)
# ---------------------------------------------------------------------------
CURATED: Dict[str, tuple] = {
    # --- Cas13a  (Type VI-A) ---
    "5W1H": ("Cas13a", "LshCas13a–crRNA binary"),
    "5WTK": ("Cas13a", "LbuCas13a ternary"),
    "5XWP": ("Cas13a", "LshCas13a ternary"),
    "6LTU": ("Cas13a", "LwaCas13a–crRNA"),
    "6NAR": ("Cas13a", "HheCas13a ternary"),
    "7DMQ": ("Cas13a", "Cas13a crRNA-target"),
    "7OS0": ("Cas13a", "LbuCas13a activated"),
    "7OS1": ("Cas13a", "LbuCas13a pre-activation"),
    # --- Cas13b  (Type VI-B) ---
    "6DTD": ("Cas13b", "PbuCas13b ternary"),
    "6AAX": ("Cas13b", "BzCas13b apo"),
    "6AAY": ("Cas13b", "BzCas13b binary"),
    # --- Cas13d  (Type VI-D / CasRx) ---
    "6IV9": ("Cas13d", "RfxCas13d–crRNA binary"),
    "6E9E": ("Cas13d", "RfxCas13d ternary"),
    "6E9F": ("Cas13d", "RfxCas13d apo"),
    "8BPO": ("Cas13d", "EsCas13d"),
    "8I3Q": ("Cas13d", "Cas13d variant"),
    # --- Cas13bt / Cas13X / Cas13Y (compact Type VI) ---
    "7WAJ": ("Cas13bt", "Cas13bt1 ternary"),
    "7WJP": ("Cas13X", "Cas13X.1 crRNA-target"),
}

RCSB_DOWNLOAD = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_SUMMARY  = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_SEARCH   = "https://search.rcsb.org/rcsbsearch/v2/query"

# Search terms that cover all Cas13 naming conventions
SEARCH_TERMS = ["Cas13", "CRISPR type VI", "C2c2"]


def _fetch_json(url: str, data: bytes | None = None, timeout: int = 30) -> Optional[dict]:
    """GET or POST JSON, return parsed dict or None on failure."""
    req = urllib.request.Request(url)
    req.add_header("Content-Type", "application/json")
    try:
        with urllib.request.urlopen(req, data=data, timeout=timeout) as resp:
            return json.loads(resp.read().decode())
    except (urllib.error.URLError, json.JSONDecodeError, TimeoutError):
        return None


def search_rcsb_cas13() -> List[str]:
    """Query RCSB full-text search for Cas13-related PDB entries."""
    all_ids: set = set()
    for term in SEARCH_TERMS:
        body = json.dumps({
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {"value": term},
            },
            "return_type": "entry",
            "request_options": {
                "results_content_type": ["experimental"],
                "paginate": {"start": 0, "rows": 1000},
            },
        }).encode()
        result = _fetch_json(RCSB_SEARCH, data=body)
        if result and "result_set" in result:
            for entry in result["result_set"]:
                pdb_id = entry.get("identifier", "").upper()
                if pdb_id:
                    all_ids.add(pdb_id)
        time.sleep(0.3)  # be polite to RCSB
    return sorted(all_ids)


def get_pdb_summary(pdb_id: str) -> dict:
    """Fetch title + resolution from RCSB REST API."""
    url = RCSB_SUMMARY.format(pdb_id=pdb_id.upper())
    info = _fetch_json(url)
    if not info:
        return {"title": "", "resolution": None}
    title = ""
    struct = info.get("struct", {})
    if isinstance(struct, dict):
        title = struct.get("title", "")
    resolution = None
    rcsb = info.get("rcsb_entry_info", {})
    if isinstance(rcsb, dict):
        resolution = rcsb.get("resolution_combined", [None])[0] if rcsb.get("resolution_combined") else None
    return {"title": title, "resolution": resolution}


def _is_cas13_title(title: str) -> bool:
    """Heuristic: does the PDB title look like it belongs to the Cas13 family?"""
    t = title.lower()
    cas13_terms = ["cas13", "c2c2", "type vi", "type-vi", "casrx",
                   "hepn nuclease vi", "vi-a", "vi-b", "vi-c", "vi-d"]
    return any(term in t for term in cas13_terms)


def classify_subtype(pdb_id: str, title: str) -> str:
    """Guess Cas13 subtype from PDB title."""
    if pdb_id in CURATED:
        return CURATED[pdb_id][0]
    t = title.lower()
    if "cas13a" in t or "c2c2" in t or "vi-a" in t:
        return "Cas13a"
    if "cas13b" in t or "vi-b" in t:
        return "Cas13b"
    if "cas13c" in t or "vi-c" in t:
        return "Cas13c"
    if "cas13d" in t or "casrx" in t or "vi-d" in t:
        return "Cas13d"
    if "cas13x" in t:
        return "Cas13X"
    if "cas13y" in t:
        return "Cas13Y"
    if "cas13bt" in t:
        return "Cas13bt"
    return "Cas13"


def download_pdb(pdb_id: str, out_dir: Path) -> Path:
    """Download PDB file from RCSB; return local path."""
    out_path = out_dir / f"{pdb_id.upper()}.pdb"
    if out_path.exists():
        return out_path
    url = RCSB_DOWNLOAD.format(pdb_id=pdb_id.upper())
    urllib.request.urlretrieve(url, str(out_path))
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Download all Cas13-family PDB references from RCSB.")
    parser.add_argument(
        "--out-dir",
        default="data/structure_pipeline/references",
        help="Directory to save PDB files (default: data/structure_pipeline/references)",
    )
    parser.add_argument(
        "--catalog",
        default="data/structure_pipeline/references/cas13_reference_catalog.json",
        help="Output JSON catalog path",
    )
    parser.add_argument("--dry-run", action="store_true", help="List PDB IDs without downloading")
    parser.add_argument("--skip-search", action="store_true", help="Only use curated list, skip RCSB API search")
    args = parser.parse_args()

    proj = Path(__file__).resolve().parent.parent
    out_dir = (proj / args.out_dir).resolve()
    catalog_path = (proj / args.catalog).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Collect PDB IDs ----
    pdb_ids: set = set(CURATED.keys())
    print(f"[*] Curated baseline: {len(pdb_ids)} Cas13 PDBs")

    if not args.skip_search:
        print("[*] Searching RCSB for additional Cas13 structures...")
        discovered = search_rcsb_cas13()
        print(f"    Raw search results: {len(discovered)} entries")
        # Filter to only Cas13-related by checking title
        verified = set()
        for pdb_id in discovered:
            if pdb_id in pdb_ids:
                verified.add(pdb_id)
                continue
            summary = get_pdb_summary(pdb_id)
            if _is_cas13_title(summary.get("title", "")):
                verified.add(pdb_id)
            time.sleep(0.15)  # rate limit
        new_ids = verified - pdb_ids
        if new_ids:
            print(f"    Discovered {len(new_ids)} additional Cas13 PDBs: {', '.join(sorted(new_ids))}")
            pdb_ids.update(new_ids)
        else:
            print("    No additional Cas13 PDBs found beyond curated list.")

    pdb_ids_sorted = sorted(pdb_ids)
    print(f"\n[*] Total Cas13 PDB IDs to process: {len(pdb_ids_sorted)}")

    # ---- Build catalog & download ----
    catalog: List[dict] = []
    for pdb_id in pdb_ids_sorted:
        summary = get_pdb_summary(pdb_id) if pdb_id not in CURATED else {}
        title = summary.get("title", "")
        resolution = summary.get("resolution")

        # Use curated info if available
        if pdb_id in CURATED:
            subtype, desc = CURATED[pdb_id]
            if not title:
                title = desc
        else:
            subtype = classify_subtype(pdb_id, title)

        entry = {
            "pdb_id": pdb_id,
            "subtype": subtype,
            "title": title,
            "resolution_A": resolution,
        }

        if args.dry_run:
            print(f"  {pdb_id}  [{subtype:8s}]  {title[:70]}")
        else:
            try:
                path = download_pdb(pdb_id, out_dir)
                entry["file"] = str(path.name)
                size_kb = path.stat().st_size / 1024
                print(f"  [+] {pdb_id}  [{subtype:8s}]  {size_kb:6.0f} KB  {title[:55]}")
            except Exception as e:
                print(f"  [!] {pdb_id}  FAILED: {e}")
                entry["file"] = None
            time.sleep(0.2)

        catalog.append(entry)

    # ---- Write catalog ----
    catalog_path.parent.mkdir(parents=True, exist_ok=True)
    with open(catalog_path, "w", encoding="utf-8") as f:
        json.dump(catalog, f, indent=2)
    print(f"\n[DONE] {len(catalog)} entries written to {catalog_path}")
    if not args.dry_run:
        pdb_count = sum(1 for e in catalog if e.get("file"))
        print(f"       {pdb_count} PDB files in {out_dir}")

    # ---- Summary by subtype ----
    from collections import Counter
    counts = Counter(e["subtype"] for e in catalog)
    print("\n  Subtype breakdown:")
    for sub, n in sorted(counts.items()):
        print(f"    {sub:10s}: {n}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
