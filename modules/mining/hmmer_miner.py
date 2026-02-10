"""
HMMER-based screening for Cas13/Type VI ORFs using Pfam HMMs (e.g. PF05168 HEPN).
Used only when USE_HMMER=1; otherwise mining uses ESM-2 only.
"""
import os
import tempfile
from pathlib import Path

import pyhmmer
from pyhmmer import easel, hmmer, plan7


def load_hmms(hmm_dir):
    """
    Load all .hmm files from hmm_dir. Returns a list of HMM objects
    (pyhmmer.plan7.HMM) for use with scan_sequences.
    """
    hmm_dir = Path(hmm_dir)
    if not hmm_dir.is_dir():
        return []
    hmms = []
    for path in sorted(hmm_dir.glob("*.hmm")):
        try:
            with plan7.HMMFile(path) as hmm_file:
                for hmm in hmm_file:
                    hmms.append(hmm)
                    break  # one HMM per file expected
        except Exception as e:
            # Log but continue with other files
            if os.getenv("HMMER_VERBOSE"):
                print(f"[hmmer_miner] Skip {path}: {e}")
    return hmms


def scan_sequences(sequences, hmm_dir=None, hmms=None, e_value_cutoff=1e-5):
    """
    Run hmmscan: sequences vs HMMs. Only sequences with at least one hit
    below e_value_cutoff are considered passing.

    sequences: list of str (protein sequences)
    hmm_dir: path to dir with .hmm files (used if hmms is None)
    hmms: optional pre-loaded list of HMMs from load_hmms()
    e_value_cutoff: max E-value to count as hit (default 1e-5)

    Returns: list of (best_name, best_evalue) or None per sequence.
    Same length as sequences; None if no hit below cutoff.
    """
    if hmms is None:
        hmms = load_hmms(hmm_dir or os.getenv("HMM_DIR", "data/hmm"))
    if not hmms:
        return [None] * len(sequences)

    alphabet = easel.Alphabet.protein()
    # Build digital sequences in order (hmmscan yields in same order)
    digitized = []
    for i, seq_str in enumerate(sequences):
        s = str(seq_str).strip()
        if not s:
            digitized.append(None)
            continue
        try:
            seq = easel.TextSequence(name=f"seq_{i}".encode(), sequence=s)
            seq_d = seq.digitize(alphabet)
            digitized.append(seq_d)
        except Exception:
            digitized.append(None)

    valid_indices = [i for i, d in enumerate(digitized) if d is not None]
    if not valid_indices:
        return [None] * len(sequences)

    valid_digitized = [digitized[i] for i in valid_indices]
    results = [None] * len(sequences)
    try:
        for k, top_hits in enumerate(hmmer.hmmscan(valid_digitized, hmms, cpus=1, E=e_value_cutoff)):
            if not top_hits:
                continue
            best = top_hits[0]
            idx = valid_indices[k]
            name = best.accession.decode() if best.accession else best.name.decode()
            evalue = getattr(best, "evalue", None) or getattr(best, "e_value", None)
            if evalue is None:
                evalue = 0.0
            if evalue <= e_value_cutoff:
                results[idx] = (name, evalue)
    except Exception as e:
        if os.getenv("HMMER_VERBOSE"):
            print(f"[hmmer_miner] hmmscan error: {e}")
    return results


def scan_single(sequence, hmm_dir=None, hmms=None, e_value_cutoff=1e-5):
    """
    Convenience: scan one sequence. Returns (best_name, best_evalue) or None.
    """
    out = scan_sequences([sequence], hmm_dir=hmm_dir, hmms=hmms, e_value_cutoff=e_value_cutoff)
    return out[0] if out else None
