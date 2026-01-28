import pandas as pd
import h5py
import numpy as np

class ARCHS4Loader:
    def __init__(self, h5_path="data/expression_data/human_matrix.h5"):
        self.h5_path = h5_path

    def get_gene_expression(self, gene_symbol):
        """
        Extracts expression data for a single gene, auto-detecting matrix orientation.
        """
        try:
            with h5py.File(self.h5_path, 'r') as f:
                # 1. Locate Gene Symbols
                # Try 'meta/genes/gene_symbol' (v2) first, then 'meta/genes/genes' (v1)
                if 'meta/genes/gene_symbol' in f:
                    gene_dset = f['meta']['genes']['gene_symbol']
                else:
                    gene_dset = f['meta']['genes']['genes']
                
                # Decode all genes to strings for searching
                all_genes = [g.decode('utf-8') for g in gene_dset[:]]
                
                if gene_symbol not in all_genes:
                    print(f"[!] Gene {gene_symbol} not found in database.")
                    return None
                
                # Get the index of the gene
                gene_idx = all_genes.index(gene_symbol)
                
                # 2. Get Sample Metadata (Tissues)
                # Try 'meta/samples/source_name_ch1'
                if 'meta/samples/source_name_ch1' in f:
                    sample_dset = f['meta']['samples']['source_name_ch1']
                else:
                    print("[!] Could not find sample source metadata.")
                    return None
                    
                tissues = [t.decode('utf-8') for t in sample_dset[:]]
                
                # 3. Smart Extraction (Check Shape)
                expression_dset = f['data']['expression']
                shape = expression_dset.shape
                
                # Case A: Matrix is (Genes, Samples) -> Standard for Bio
                if shape[0] == len(all_genes):
                    # print(f"[*] Detected (Genes, Samples) orientation: {shape}")
                    expression_values = expression_dset[gene_idx, :]
                    
                # Case B: Matrix is (Samples, Genes) -> Compressed format
                elif shape[1] == len(all_genes):
                    # print(f"[*] Detected (Samples, Genes) orientation: {shape}")
                    expression_values = expression_dset[:, gene_idx]
                    
                else:
                    print(f"[!] Shape mismatch! Matrix: {shape}, Genes: {len(all_genes)}, Samples: {len(tissues)}")
                    return None

                # 4. Return Data
                return pd.DataFrame({
                    'Tissue': tissues,
                    'Expression': expression_values
                })
                
        except FileNotFoundError:
            print(f"[!] Database not found at {self.h5_path}.")
            return None
        except Exception as e:
            print(f"[!] Error reading H5 file: {e}")
            return None