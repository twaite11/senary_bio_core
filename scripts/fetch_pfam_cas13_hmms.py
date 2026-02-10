#!/usr/bin/env python3
"""
Download Pfam Cas13/HEPN HMMs into data/hmm/ for USE_HMMER=1 mining.
Requires: HMMER (hmmfetch) in PATH, or run the manual steps in data/hmm/README.md.
"""
import os
import sys
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path

# Pfam accessions for Cas13/Type VI mining (exact denotations)
CAS13_PF_IDS = ["PF05168"]  # HEPN

PFAM_HMM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
ROOT = Path(__file__).resolve().parent.parent
HMM_DIR = ROOT / "data" / "hmm"


def main():
    HMM_DIR.mkdir(parents=True, exist_ok=True)
    if not shutil.which("hmmfetch"):
        print(
            "[!] hmmfetch not found. Install HMMER (http://hmmer.org) or run manually:\n"
            "    See data/hmm/README.md for wget + hmmfetch steps.",
            file=sys.stderr,
        )
        sys.exit(1)

    print("[*] Downloading Pfam-A.hmm.gz...")
    with tempfile.NamedTemporaryFile(suffix=".hmm.gz", delete=False) as tmp_gz:
        tmp_gz_path = tmp_gz.name
    try:
        import urllib.request
        urllib.request.urlretrieve(PFAM_HMM_URL, tmp_gz_path)
    except Exception as e:
        print(f"[!] Download failed: {e}", file=sys.stderr)
        sys.exit(1)

    print("[*] Decompressing...")
    tmp_hmm = tmp_gz_path.replace(".gz", "")
    try:
        with gzip.open(tmp_gz_path, "rb") as f_in:
            with open(tmp_hmm, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
    except Exception as e:
        print(f"[!] Decompress failed: {e}", file=sys.stderr)
        os.unlink(tmp_gz_path)
        sys.exit(1)
    os.unlink(tmp_gz_path)

    for pf in CAS13_PF_IDS:
        out_path = HMM_DIR / f"{pf}.hmm"
        print(f"[*] Fetching {pf} -> {out_path}...")
        try:
            with open(out_path, "wb") as f:
                subprocess.run(
                    ["hmmfetch", tmp_hmm, pf],
                    stdout=f,
                    check=True,
                )
        except subprocess.CalledProcessError as e:
            print(f"[!] hmmfetch {pf} failed: {e}", file=sys.stderr)
        else:
            print(f"    [+] {out_path}")

    try:
        os.unlink(tmp_hmm)
    except OSError:
        pass
    print(f"[+] Done. HMMs in {HMM_DIR}")


if __name__ == "__main__":
    main()
