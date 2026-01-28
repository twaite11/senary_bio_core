import re
from Bio import SeqIO

class HEPNFilter:
    def __init__(self):
        # The classic R-X4-H motif essential for RNA cleavage
        # "R" = Arginine, "." = Any, "{4}" = 4 times, "H" = Histidine
        self.hepn_motif = re.compile(r"R.{4}H")

    def filter_candidates(self, fasta_path):
        """
        Scans a FASTA file and returns only enzymes with TWO HEPN domains.
        """
        valid_enzymes = []
        print(f"[*] Scanning {fasta_path} for molecular scissors...")

        for record in SeqIO.parse(fasta_path, "fasta"):
            seq_str = str(record.seq).upper()
            
            # Find all HEPN domains
            matches = self.hepn_motif.findall(seq_str)
            
            # Cas13d needs TWO HEPN domains to cut
            if len(matches) >= 2:
                valid_enzymes.append(record)
                print(f"   [+] Valid Candidate: {record.id} (Found {len(matches)} HEPN motifs)")
        
        print(f"[*] Retained {len(valid_enzymes)} functional candidates.")
        return valid_enzymes