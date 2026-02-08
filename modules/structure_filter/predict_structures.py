"""
Run OmegaFold (or ESMFold) on a FASTA; output one PDB per sequence.
Used by the pipeline for in-silico structure filtering.
"""
import os
import subprocess
import argparse
from pathlib import Path
from typing import Optional


def run_omegafold(
    input_fasta: str,
    output_dir: str,
    num_cycle: int = 3,
    device: Optional[str] = None,
    omegafold_repo: Optional[str] = None,
) -> None:
    """Run OmegaFold on input FASTA. Writes PDBs to output_dir."""
    input_path = Path(input_fasta).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    repo = omegafold_repo or os.environ.get("OMEGAFOLD_REPO")
    if repo:
        main_py = Path(repo).resolve() / "main.py"
        if not main_py.exists():
            raise FileNotFoundError(f"OmegaFold main.py not found: {main_py}")
        cmd = ["python", str(main_py), str(input_path), str(out_path)]
    else:
        cmd = ["omegafold", str(input_path), str(out_path)]

    cmd.extend(["--num_cycle", str(num_cycle)])
    if device:
        cmd.extend(["--device", device])
    subprocess.run(cmd, check=True, cwd=repo if repo else None)


def main():
    parser = argparse.ArgumentParser(description="Predict structures with OmegaFold for pipeline filter.")
    parser.add_argument("--input", default="data/design/drift_variants.fasta", help="Input FASTA")
    parser.add_argument("--output-dir", default="data/structure_pipeline/structures/omegafold", help="Output PDB dir")
    parser.add_argument("--num-cycle", type=int, default=3)
    parser.add_argument("--device", default=None, help="cuda or cpu")
    parser.add_argument("--omegafold-repo", default=None, help="Path to OmegaFold repo (or OMEGAFOLD_REPO)")
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"[!] Input not found: {args.input}")
        return 1
    run_omegafold(
        args.input,
        args.output_dir,
        num_cycle=args.num_cycle,
        device=args.device,
        omegafold_repo=args.omegafold_repo,
    )
    print(f"[+] Structures written to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
