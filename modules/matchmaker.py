import pandas as pd
from Bio import SeqIO
import random
import os
import sys

class MasterMatchmaker:
    def __init__(self, enzyme_fasta, target_file, disease_matrix_file=None):
        self.enzyme_fasta = enzyme_fasta
        self.target_file = target_file
        self.disease_matrix_file = disease_matrix_file
        self.results = []
        self.disease_map = {}

    def _load_enzymes(self):
        enzymes = []
        if not os.path.exists(self.enzyme_fasta):
            print(f"[!] Warning: Enzyme file '{self.enzyme_fasta}' not found.")
            print("    -> Using 'Mock Enzymes' for simulation.")
            # Generating 3 mock enzymes to show the loop functionality
            return [
                {'id': 'Cas13d_Yellowstone_V1', 'pfs_rule': 'NOT_G'},
                {'id': 'Cas13d_DeepSea_V4', 'pfs_rule': 'NOT_G'},
                {'id': 'Cas13d_SaltLake_V9', 'pfs_rule': 'NOT_G'}
            ]
            
        try:
            for record in SeqIO.parse(self.enzyme_fasta, "fasta"):
                enzymes.append({'id': record.id, 'pfs_rule': 'NOT_G'})
        except Exception as e:
            print(f"[!] Error parsing enzyme FASTA: {e}")
            return []

        if not enzymes:
            return [{'id': 'Cas13d_Mock_V1', 'pfs_rule': 'NOT_G'}]
            
        return enzymes

    def _load_disease_map(self):
        """
        Loads the Matrix CSV (Rows=Cancer, Cols=Fusion) to map Fusions -> Diseases.
        """
        if not self.disease_matrix_file or not os.path.exists(self.disease_matrix_file):
            print("[!] Warning: Disease matrix file not found. Diseases will be 'Unknown'.")
            return

        print(f"[*] Loading Disease Map from {self.disease_matrix_file}...")
        try:
            # Load matrix
            df = pd.read_csv(self.disease_matrix_file)
            
            # Assuming format: First column is Cancer Type, subsequent cols are Fusions
            # We want to iterate columns (Fusions) and find which rows (Cancers) have a count > 0
            
            # Set the first column (Cancer) as index
            df.set_index(df.columns[0], inplace=True)
            
            # Iterate through fusions (columns)
            for fusion in df.columns:
                # Find cancers where this fusion has a non-zero count
                associated_cancers = df.index[df[fusion] > 0].tolist()
                if associated_cancers:
                    self.disease_map[fusion] = ", ".join(associated_cancers)
                    
            print(f"   [+] Mapped diseases for {len(self.disease_map)} fusions.")
            
        except Exception as e:
            print(f"[!] Error loading disease map: {e}")

    def _load_targets(self):
        print(f"[*] Loading Targets from {self.target_file}...")
        
        if not os.path.exists(self.target_file):
            print(f"[!] Critical: Target file '{self.target_file}' not found.")
            return pd.DataFrame()
        
        try:
            df = pd.read_csv(self.target_file)
            
            if 'fusionsss' in df.columns:
                df.rename(columns={'fusionsss': 'Fusion_Name'}, inplace=True)
            
            # Deduplicate! Your file has multiple rows for the same fusion.
            # We keep the first occurrence which usually has the count.
            initial_len = len(df)
            df = df.drop_duplicates(subset=['Fusion_Name'])
            print(f"   [+] Loaded {len(df)} unique targets (deduplicated from {initial_len}).")
            
            return df
            
        except Exception as e:
            print(f"[!] Error reading file: {e}")
            return pd.DataFrame()

    def _get_target_sequence(self, fusion_name):
        # SIMULATION: Generating random RNA Junction (60nt)
        # Real pipeline: Call NCBI Entrez API here.
        bases = ['A', 'C', 'T', 'G']
        return "".join(random.choice(bases) for _ in range(60))

    def run_matching(self):
        # 1. Load Data
        enzymes = self._load_enzymes()
        targets_df = self._load_targets()
        self._load_disease_map()

        if targets_df.empty or not enzymes:
            print("[!] Aborting: Missing Data.")
            return

        print(f"[*] Screening {len(enzymes)} Enzymes against top 50 Targets...")
        
        # Select top 50 unique targets
        subset = targets_df.head(50)

        for index, row in subset.iterrows():
            target_name = row['Fusion_Name']
            patient_count = row.get('count', 'N/A')
            
            # Lookup Disease
            disease = self.disease_map.get(target_name, "Unknown")
            
            # Get Sequence (Simulation of Junction)
            rna_seq = self._get_target_sequence(target_name)
            
            # Test EVERY Enzyme
            for enzyme in enzymes:
                score, cut_sites = self._calculate_cut_efficiency(enzyme, rna_seq)
                
                # We save it if it has at least one valid cut site
                if score > 0:
                    self.results.append({
                        'Target_Fusion': target_name,
                        'Associated_Disease': disease,
                        'Patient_Count': patient_count,
                        'Enzyme_Variant': enzyme['id'],
                        'Valid_Cut_Sites': score,
                        'Junction_Sequence_Sim': rna_seq[:15] + "..." # Preview
                    })

    def _calculate_cut_efficiency(self, enzyme, rna_seq):
        valid_spacers = []
        spacer_len = 22
        for i in range(len(rna_seq) - spacer_len):
            spacer = rna_seq[i : i+spacer_len]
            pfs_index = i + spacer_len
            if pfs_index < len(rna_seq):
                if rna_seq[pfs_index] != 'G':
                    valid_spacers.append(spacer)
        return len(valid_spacers), valid_spacers

    def save_leads(self, filename="lead_candidates.csv"):
        if not self.results:
            print("[!] No matches found.")
            return
        
        df = pd.DataFrame(self.results)
        
        # Sort by Disease, then Score
        df = df.sort_values(by=['Associated_Disease', 'Valid_Cut_Sites'], ascending=[True, False])
        
        df.to_csv(filename, index=False)
        print(f"\n[SUCCESS] Generated '{filename}' with {len(df)} candidates.")
        print("--- TOP MATCHES PREVIEW ---")
        print(df[['Target_Fusion', 'Associated_Disease', 'Enzyme_Variant', 'Valid_Cut_Sites']].head(10))

if __name__ == "__main__":
    # --- CONFIGURATION ---
    
    # 1. Enzyme File
    ENZYME_FILE = "data/raw_sequences/search_20260128_103022.fasta"
    
    # 2. Target File (The List)
    TARGET_FILE = "data/known_fusions.csv"
    
    # 3. Disease Map File (The Matrix)
    # RENAME 'Recurrent_table.xlsx - KB_and_Pub_Recur_per_cancer.csv' to this:
    DISEASE_FILE = "data/disease_matrix_known.csv"
    
    print(f"--- Collateral Bio: Matchmaker v2.1 ---")
    matcher = MasterMatchmaker(ENZYME_FILE, TARGET_FILE, DISEASE_FILE)
    matcher.run_matching()
    matcher.save_leads()