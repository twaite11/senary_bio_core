"""
Run OmegaFold (or ESMFold) on a FASTA; output one PDB per sequence.
Used by the pipeline for in-silico structure filtering.
Processes sequences in small batches to avoid GPU VRAM OOM.
"""
import os
import subprocess
import argparse
import tempfile
from pathlib import Path
from typing import Optional

from Bio import SeqIO


def _find_omegafold_venv_python(repo_path: Path) -> Optional[str]:
    """
    Find OmegaFold venv Python executable.
    Checks common venv locations relative to the OmegaFold repo.
    Returns path to Python executable if found, None otherwise.
    """
    # Common venv names to check
    venv_names = [
        "venv_omegafold1",
        "venv_omegafold",
        "omegafold_venv",
        "venv",
        "env",
    ]
    
    # Check parent directory (e.g., /workspace/venv_omegafold1 when repo is /workspace/OmegaFold)
    parent_dir = repo_path.parent
    for venv_name in venv_names:
        venv_path = parent_dir / venv_name
        if venv_path.exists() and venv_path.is_dir():
            # Check for Python executable (Linux/Mac)
            python_exe = venv_path / "bin" / "python"
            if python_exe.exists():
                return str(python_exe)
            # Check for Python executable (Windows)
            python_exe = venv_path / "Scripts" / "python.exe"
            if python_exe.exists():
                return str(python_exe)
    
    # Also check if OMEGAFOLD_VENV is explicitly set
    omegafold_venv = os.environ.get("OMEGAFOLD_VENV")
    if omegafold_venv:
        venv_path = Path(omegafold_venv).resolve()
        python_exe = venv_path / "bin" / "python"
        if python_exe.exists():
            return str(python_exe)
        python_exe = venv_path / "Scripts" / "python.exe"
        if python_exe.exists():
            return str(python_exe)
    
    return None


def _run_omegafold_single_fasta(
    input_fasta: str,
    output_dir: str,
    num_cycle: int,
    device: Optional[str],
    omegafold_repo: Optional[str],
    python_exe: str,
    main_py: Path,
    repo: Optional[str],
) -> None:
    """Run OmegaFold once on a single FASTA file."""
    cmd = [python_exe, str(main_py), input_fasta, output_dir]
    cmd.extend(["--num_cycle", str(num_cycle)])
    if device:
        cmd.extend(["--device", device])
    subprocess.run(cmd, check=True, cwd=repo if repo else None)


def run_omegafold(
    input_fasta: str,
    output_dir: str,
    num_cycle: int = 3,
    device: Optional[str] = None,
    omegafold_repo: Optional[str] = None,
    batch_size: int = 5,
) -> None:
    """Run OmegaFold on input FASTA. Writes PDBs to output_dir. Processes in batches to limit VRAM."""
    input_path = Path(input_fasta).resolve()
    out_path = Path(output_dir).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    repo = omegafold_repo or os.environ.get("OMEGAFOLD_REPO")
    if repo:
        repo_path = Path(repo).resolve()
        # Check if repo path is valid (not a placeholder)
        if not repo_path.exists() or "path/to" in str(repo_path).lower() or repo_path == Path("/path/to/OmegaFold"):
            raise FileNotFoundError(
                f"Invalid OmegaFold repository path: {repo}\n"
                f"Please set OMEGAFOLD_REPO environment variable to the actual OmegaFold repository path, "
                f"or use --omegafold-repo argument, or use --skip-structure to skip structure prediction."
            )
        main_py = repo_path / "main.py"
        if not main_py.exists():
            raise FileNotFoundError(
                f"OmegaFold main.py not found at: {main_py}\n"
                f"Repository path: {repo_path}\n"
                f"Please verify that OmegaFold is installed at this location."
            )
        venv_python = _find_omegafold_venv_python(repo_path)
        python_exe = venv_python if venv_python else "python"
        if venv_python:
            print(f"[*] Using OmegaFold venv Python: {venv_python}")
        else:
            print(f"[*] Using system Python (OmegaFold venv not found, using: {python_exe})")
    else:
        python_exe = "python"
        main_py = None
        repo_path = None

    records = list(SeqIO.parse(input_path, "fasta"))
    total = len(records)
    if total == 0:
        print("[!] No sequences in input FASTA.")
        return

    if batch_size < 1:
        batch_size = total
    n_batches = (total + batch_size - 1) // batch_size
    print(f"[*] OmegaFold: {total} sequences in {n_batches} batch(es) of up to {batch_size} (to reduce VRAM).")

    for i in range(0, total, batch_size):
        chunk = records[i : i + batch_size]
        batch_num = i // batch_size + 1
        fd, temp_path = tempfile.mkstemp(suffix=".fasta", prefix="omegafold_batch_")
        try:
            os.close(fd)
            SeqIO.write(chunk, temp_path, "fasta")
            print(f"[*] OmegaFold batch {batch_num}/{n_batches} ({len(chunk)} sequences)...")
            if repo and main_py is not None:
                _run_omegafold_single_fasta(
                    temp_path,
                    str(out_path),
                    num_cycle,
                    device,
                    omegafold_repo,
                    python_exe,
                    main_py,
                    repo,
                )
            else:
                subprocess.run(
                    ["omegafold", temp_path, str(out_path), "--num_cycle", str(num_cycle)]
                    + (["--device", device] if device else []),
                    check=True,
                )
        finally:
            try:
                os.unlink(temp_path)
            except OSError:
                pass

    print(f"[+] Structures written to {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Predict structures with OmegaFold for pipeline filter.")
    parser.add_argument("--input", default="data/design/drift_variants.fasta", help="Input FASTA")
    parser.add_argument("--output-dir", default="data/structure_pipeline/structures/omegafold", help="Output PDB dir")
    parser.add_argument("--num-cycle", type=int, default=3)
    parser.add_argument("--device", default=None, help="cuda or cpu")
    parser.add_argument("--omegafold-repo", default=None, help="Path to OmegaFold repo (or OMEGAFOLD_REPO)")
    parser.add_argument("--batch-size", type=int, default=5, help="Sequences per batch to limit VRAM (default 5)")
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
        batch_size=args.batch_size,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
