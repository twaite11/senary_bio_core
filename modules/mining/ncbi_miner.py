import os
from Bio import Entrez, SeqIO
from datetime import datetime

# REQUIRED: Set your email for NCBI tracking
Entrez.email = "founder@collateralbio.com" 

class EnzymeMiner:
    def __init__(self, output_dir="data/raw_sequences"):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
    
    def search_and_fetch(self, query="Cas13d", max_results=50):
        """
        1. Search NCBI Protein DB
        2. Fetch FASTA sequences
        3. Save with IP Timestamp
        """
        print(f"[*] Searching NCBI for: '{query}'...")
        
        # Search
        handle = Entrez.esearch(db="protein", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        if not id_list:
            print("[!] No candidates found.")
            return []

        # Fetch
        print(f"[*] Fetching {len(id_list)} sequences...")
        handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text")
        sequences = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        
        # Save (The IP Event)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"{self.output_dir}/search_{timestamp}.fasta"
        SeqIO.write(sequences, filename, "fasta")
        
        print(f"[SUCCESS] Saved {len(sequences)} candidates to {filename}")
        return sequences