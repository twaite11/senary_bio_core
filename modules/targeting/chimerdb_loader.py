import pandas as pd

class ChimerDBLoader:
    def __init__(self, db_path="data/ChimerDB4_Recurrence.csv"):
        # You download this from the ChimerDB website
        self.db_path = db_path
        self.df = None

    def load_db(self):
        print(f"[*] Loading ChimerDB from {self.db_path}...")
        try:
            # ChimerDB often uses tab or comma. Adjust 'sep' accordingly.
            self.df = pd.read_csv(self.db_path, sep='\t')
            print(f"[SUCCESS] Loaded {len(self.df)} fusion records.")
        except FileNotFoundError:
            print("[!] Database file not found. Please download ChimerDB4.")

    def find_tumor_exclusive_fusions(self, tumor_type):
        """
        Filters for fusions that are High in Tumor, Low/Zero in other tissues.
        """
        if self.df is None:
            self.load_db()
            if self.df is None: return

        print(f"[*] Hunting for exclusive fusions in: {tumor_type}...")
        
        # 1. Filter by Disease (Cancer Type)
        # Note: Column names vary by DB version. Assuming standard headers.
        tumor_hits = self.df[self.df['Cancer_Type'].str.contains(tumor_type, case=False, na=False)]
        
        # 2. Sort by Frequency (How common is it?)
        # 'Frequency' usually indicates how many samples had it.
        top_hits = tumor_hits.sort_values(by='Frequency', ascending=False)
        
        return top_hits[['Fusion_Pair', 'Cancer_Type', 'Frequency', 'Kinase_Fusion']]

# Usage
if __name__ == "__main__":
    loader = ChimerDBLoader()
    # Example: Find fusions in Brain Cancer (GBM)
    loader.find_tumor_exclusive_fusions("GBM")