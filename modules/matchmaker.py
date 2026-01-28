import pandas as pd
from Bio import SeqIO
import random
import os
import sys

class MasterMatchmaker:
    def __init__(self, enzyme_fasta, target_file):
        self.enzyme_fasta = enzyme_fasta
        self.target_file = target_file
        self.results = []

    def _load_enzymes(self):
        enzymes = []
        # 1. Check if file exists
        if not os.path.exists(self.enzyme_fasta):
            print(f"[!] Warning: Enzyme file '{self.enzyme_fasta}' not found.")
            print("    -> Switching to 'Mock Enzyme' mode for simulation.")
            return [{'id': 'Cas13d_Mock_Variant_01', 'pfs_rule': 'NOT_G'}]
            
        # 2. Try to parse with Biopython
        try:
            for record in SeqIO.parse(self.enzyme_fasta, "fasta"):
                enzymes.append({'id': record.id, 'pfs_rule': 'NOT_G'})
        except Exception as e:
            print(f"[!] Error parsing enzyme FASTA: {e}")
            return []

        if not enzymes:
            print("[!] Enzyme file is empty. Using Mock Enzyme.")
            return [{'id': 'Cas13d_Mock_Variant_01', 'pfs_rule': 'NOT_G'}]
            
        return enzymes

    def _load_targets(self):
        print(f"[*] Loading Targets from {self.target_file}...")
        
        if not os.path.exists(self.target_file):
            print(f"[!] Critical Error: Target file '{self.target_file}' not found.")
            print("    -> Did you rename the downloaded CSV to 'known_fusions.csv'?")
            return pd.DataFrame()
        
        try:
            # Read CSV. The snippet shows comma-separated values.
            df = pd.read_csv(self.target_file)
            
            # --- THE FIX FOR YOUR SPECIFIC FILE ---
            # Your file header is: fusionsss, count, Pub supported, KB supported
            if 'fusionsss' in df.columns:
                df.rename(columns={'fusionsss': 'Fusion_Name'}, inplace=True)
            
            # Basic Validation
            if 'Fusion_Name' not in df.columns:
                print(f"[!] Columns found: {list(df.columns)}")
                print("[!] Critical Error: Could not find 'Fusion_Name' or 'fusionsss' column.")
                return pd.DataFrame()
                
            print(f"   [+] Successfully loaded {len(df)} targets.")
            return df
            
        except Exception as e:
            print(f"[!] Error reading file: {e}")
            return pd.DataFrame()

    def _get_target_sequence(self, fusion_name):
        # SIMULATION: In Phase 2, we connect to NCBI/RefSeq.
        # For the Pitch Deck, we simulate the RNA context (60nt).
        bases = ['A', 'C', 'T', 'G']
        return "".join(random.choice(bases) for _ in range(60))

    def run_matching(self):
        enzymes = self._load_enzymes()
        targets_df = self._load_targets()

        if targets_df.empty or not enzymes:
            print("[!] Aborting matchmaker due to missing inputs.")
            return

        print(f"[*] Screening {len(enzymes)} Enzymes against Targets...")
        
        # We screen the top 50 targets for speed
        subset = targets_df.head(50)

        for index, row in subset.iterrows():
            target_name = row['Fusion_Name']
            patient_count = row.get('count', 'N/A') # Matches your file's 'count' column
            
            # 1. Get RNA Sequence (Simulation)
            rna_seq = self._get_target_sequence(target_name)
            
            # 2. Test Enzymes
            for enzyme in enzymes:
                score, cut_sites = self._calculate_cut_efficiency(enzyme, rna_seq)
                
                if score > 0:
                    self.results.append({
                        'Target': target_name,
                        'Patient_Count': patient_count,
                        'Enzyme_ID': enzyme['id'],
                        'Valid_Cut_Sites': score,
                        'Best_Spacer': cut_sites[0]
                    })

    def _calculate_cut_efficiency(self, enzyme, rna_seq):
        """
        Calculates how many 22nt spacers exist in the target sequence
        that respect the Cas13d PFS rule (No G at 3' end).
        """
        valid_spacers = []
        spacer_len = 22
        
        for i in range(len(rna_seq) - spacer_len):
            spacer = rna_seq[i : i+spacer_len]
            pfs_index = i + spacer_len
            
            # Check boundary and PFS rule
            if pfs_index < len(rna_seq):
                if rna_seq[pfs_index] != 'G':
                    valid_spacers.append(spacer)
                    
        return len(valid_spacers), valid_spacers

    def save_leads(self, filename="lead_candidates.csv"):
        if not self.results:
            print("[!] No matches found.")
            return
        
        df = pd.DataFrame(self.results)
        
        # Sort by Cut Sites (Druggability) then Patient Count (Market Size)
        df = df.sort_values(by=['Valid_Cut_Sites', 'Patient_Count'], ascending=[False, False])
        
        # Save
        df.to_csv(filename, index=False)
        print(f"\n[SUCCESS] Generated '{filename}' with {len(df)} candidates.")
        print("--- TOP 5 MATCHES ---")
        print(df[['Target', 'Patient_Count', 'Valid_Cut_Sites']].head())

if __name__ == "__main__":
    # --- CONFIGURATION ---
    
    # 1. Enzyme File (Update with your real timestamped file if you have it)
    ENZYME_FILE = "data/raw_sequences/search_20260128_103022.fasta" 
    
    # 2. Target File (Defaulting to the Validation Set you just uploaded)
    # CHANGE THIS to "data/novel_fusions.csv" when you want to run discovery.
    TARGET_FILE = "data/known_fusions.csv"
    
    print(f"--- Collateral Bio: Matchmaker v2.0 ---")
    matcher = MasterMatchmaker(ENZYME_FILE, TARGET_FILE)
    matcher.run_matching()
    matcher.save_leads()