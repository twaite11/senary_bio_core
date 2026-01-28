import vcf
import pandas as pd
import os

class MutationMiner:
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
    
    def find_specific_mutations(self, target_gene, min_depth=10):
        """
        Parses a VCF to find mutations in a specific gene.
        Useful for validating if 'KRAS' has the specific 'G12C' mutation in your samples.
        """
        print(f"[*] Mining VCF for mutations in {target_gene}...")
        mutations = []
        
        try:
            vcf_reader = vcf.Reader(filename=self.vcf_path)
            
            for record in vcf_reader:
                # VCFs often use CHROM:POS, so we filter by gene name in INFO field if available
                # Or you must know the genomic coordinates of your gene.
                # Here we assume the VCF is annotated (like SnpEff output)
                
                info = record.INFO
                if 'GENE' in info and info['GENE'][0] == target_gene:
                    for sample in record.samples:
                        if sample['GT'] != '0/0': # If sample has the mutation
                            mutations.append({
                                'Sample': sample.sample,
                                'Position': record.POS,
                                'Ref': record.REF,
                                'Alt': str(record.ALT[0]),
                                'Depth': sample['DP']
                            })
                            
            return pd.DataFrame(mutations)
            
        except FileNotFoundError:
            print("[!] VCF file not found.")
            return None