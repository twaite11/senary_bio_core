import sys
import os

# Ensure Python can find the 'modules' folder
sys.path.append(os.getcwd())

from modules.targeting.archs4_loader import ARCHS4Loader

# Path to the data file
h5_path = "data/expression_data/human_matrix.h5"

if not os.path.exists(h5_path):
    print(f"Error: Could not find file at {h5_path}")
    print("Did you download human_matrix.h5 from ARCHS4 and put it in data/expression_data/?")
else:
    print(f"[*] Loading database from {h5_path}...")
    loader = ARCHS4Loader(h5_path)
    
    # Check expression of a known cancer gene (KRAS)
    print("[*] Querying KRAS expression...")
    result = loader.get_gene_expression("KRAS")
    
    if result is not None:
        print("\n[SUCCESS] Data loaded! Here are the first 5 tissues:")
        print(result.head())