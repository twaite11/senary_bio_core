import pandas as pd
import os

class ExcelSplitter:
    def __init__(self, excel_path, output_dir="data"):
        self.excel_path = excel_path
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def run_split(self):
        print(f"[*] Opening Excel file: {self.excel_path}...")
        
        try:
            # Read all sheets at once
            xls = pd.ExcelFile(self.excel_path)
            sheet_names = xls.sheet_names
            print(f"[*] Found sheets: {sheet_names}")
            
            # Define specific renaming rules for your project
            # Sheet Name -> Clean Filename
            rename_map = {
                'HRS_Recurrent_known': 'known_fusions.csv',  # Validation Set
                'HRS_Recurrent_novel': 'novel_fusions.csv',  # Discovery Set
            }

            for sheet in sheet_names:
                print(f"   -> Processing sheet: '{sheet}'...")
                df = pd.read_excel(xls, sheet_name=sheet)
                
                # Determine output filename
                if sheet in rename_map:
                    filename = rename_map[sheet]
                    print(f"      [!] Renaming to critical file: {filename}")
                else:
                    # Clean up other sheet names to be safe filenames
                    safe_name = sheet.replace(" ", "_").replace("&", "and")
                    filename = f"{safe_name}.csv"
                
                output_path = os.path.join(self.output_dir, filename)
                df.to_csv(output_path, index=False)
                print(f"      [+] Saved to {output_path}")

            print("\n[SUCCESS] Split complete. Your data folder is ready.")

        except FileNotFoundError:
            print(f"[!] Error: Could not find file '{self.excel_path}'")
            print("    Make sure the Excel file is inside your 'data' folder or provide the full path.")
        except Exception as e:
            print(f"[!] Error processing file: {e}")

if __name__ == "__main__":
    # Update this path to where your .xlsx file actually lives
    # Assuming you dropped it in the data folder
    EXCEL_FILE = "data/Recurrent_table.xlsx"
    
    splitter = ExcelSplitter(EXCEL_FILE)
    splitter.run_split()