import pandas as pd

class CCLEFilter:
    def __init__(self, metadata_path="data/ccle_metadata.csv"):
        # You download this csv from DepMap portal (free)
        self.metadata = pd.read_csv(metadata_path)

    def find_cell_lines(self, disease=None, mutation=None):
        """
        Filters cell lines by disease type and/or mutation status.
        """
        print(f"[*] Filtering CCLE for Disease='{disease}' / Mutation='{mutation}'...")
        
        df = self.metadata
        
        # 1. Filter by Disease (e.g., "Lung Adenocarcinoma")
        if disease:
            # Flexible case-insensitive search
            df = df[df['primary_disease'].str.contains(disease, case=False, na=False)]
        
        # 2. Filter by Mutation (requires mutation column in metadata or separate merge)
        # Note: DepMap metadata usually has 'Subtype' or 'Lineage'. 
        # For specific mutations (KRAS G12C), you merge with the Mutation csv.
        
        print(f"[*] Found {len(df)} matching cell lines.")
        return df[['DepMap_ID', 'cell_line_name', 'primary_disease', 'sex', 'age']]

# Usage Example
if __name__ == "__main__":
    # Download 'sample_info.csv' from https://depmap.org/portal/download/
    ccle = CCLEFilter("data/sample_info.csv")
    
    # Find all Lung Cancer lines to screen for our targets
    targets = ccle.find_cell_lines(disease="Lung")
    print(targets.head())
    
    # Save the ID list to feed into the AWS downloader
    targets['DepMap_ID'].to_csv("data/lung_cancer_ids.txt", index=False)