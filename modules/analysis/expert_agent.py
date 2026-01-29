import pandas as pd
import json
import os
import time
import google.generativeai as genai
import sys
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Optional Import for Safety Checks
sys.path.append(os.getcwd())
try:
    from modules.targeting.archs4_loader import ARCHS4Loader
except ImportError:
    ARCHS4Loader = None

class ExpertAgent:
    def __init__(self, candidates_file, prompt_file):
        self.candidates_file = candidates_file
        self.prompt_file = prompt_file
        self.safety_loader = None
        self.model = None
        
        # Load API key from .env file
        api_key = os.getenv("GEMINI_API_KEY")
        if not api_key:
            print("Warning: GEMINI_API_KEY not found in environment variables. AI analysis will be disabled.")
        else:
            try:
                genai.configure(api_key=api_key)
                # Use gemini-flash-latest for free tier (faster and free)
                # This is the latest flash model available in the free tier
                self.model = genai.GenerativeModel('gemini-flash-latest')
                print("[OK] Gemini API configured successfully")
            except Exception as e:
                print(f"Error configuring Gemini API: {e}")
                self.model = None
        
        if ARCHS4Loader:
            h5_path = "data/expression_data/human_matrix.h5"
            if os.path.exists(h5_path):
                self.safety_loader = ARCHS4Loader(h5_path)

    def load_candidates(self, max_candidates=10):
        """Load candidates from CSV file, limiting to max_candidates for faster processing"""
        if not os.path.exists(self.candidates_file):
            return []
        return pd.read_csv(self.candidates_file).head(max_candidates).to_dict('records')

    def check_safety(self, target_name):
        if not self.safety_loader: return 50, ["Safety Check Skipped"]
        
        parts = target_name.replace("--", "-").split("-")
        genes = [p for p in parts if len(p) > 2]
        risks = []
        
        for gene in genes:
            df = self.safety_loader.get_gene_expression(gene)
            if df is not None:
                summary = df.groupby('Tissue')['Expression'].mean()
                for organ in ['brain', 'heart', 'liver', 'lung', 'kidney']:
                    matches = summary[summary.index.str.contains(organ, case=False, na=False)]
                    if not matches.empty and matches.mean() > 100:
                        risks.append(f"{gene} in {organ}")
        
        return (100 - len(risks)*25, risks) if risks else (100, ["Clean Profile"])

    def load_expert_prompt(self):
        """Load the expert persona prompt from file"""
        if not os.path.exists(self.prompt_file):
            return ""
        try:
            with open(self.prompt_file, 'r', encoding='utf-8') as f:
                return f.read()
        except Exception as e:
            print(f"Warning: Could not load prompt file: {e}")
            return ""
    
    def analyze_ai(self, candidate, risks):
        """Analyze a candidate using Gemini AI"""
        if not self.model:
            return {
                "Verdict": "N/A",
                "Rationale": "AI model not configured",
                "Strategy": "N/A"
            }
        
        # Load expert persona prompt
        expert_prompt = self.load_expert_prompt()
        
        # Build the analysis prompt
        prompt = f"""{expert_prompt}

TARGET DATA:
- Target Fusion: {candidate.get('Target_Fusion', 'Unknown')}
- Associated Disease: {candidate.get('Associated_Disease', 'Unknown')}
- Patient Count: {candidate.get('Patient_Count', 0)}
- Safety Risks: {', '.join(risks) if risks else 'None identified'}

Please analyze this target and return your evaluation as JSON with the following structure:
{{
    "Verdict": "GO" | "NO-GO" | "HOLD",
    "Scientific_Rationale": "Brief explanation",
    "Screening_Strategy": "Recommended approach",
    "Risk_Factors": "Key concerns or advantages"
}}
"""
        
        try:
            # Use modern Gemini API with structured output
            generation_config = {
                "response_mime_type": "application/json",
                "temperature": 0.7,
            }
            
            response = self.model.generate_content(
                prompt,
                generation_config=generation_config
            )
            
            # Parse JSON response
            if response.text:
                return json.loads(response.text)
            else:
                print(f"Warning: Empty response from Gemini API")
                return {}
                
        except json.JSONDecodeError as e:
            print(f"Error parsing JSON response: {e}")
            print(f"Response text: {response.text if 'response' in locals() else 'N/A'}")
            return {}
        except Exception as e:
            print(f"Error calling Gemini API: {e}")
            return {}

    def run(self, max_candidates=10):
        """Run analysis on candidates"""
        print("--- Collateral Bio: AI Analysis ---")
        candidates = self.load_candidates(max_candidates=max_candidates)
        
        if not candidates:
            print("[!] No candidates found to analyze.")
            return
        
        print(f"[*] Analyzing {len(candidates)} candidates...")
        results = []
        
        for i, cand in enumerate(candidates, 1):
            target = cand.get('Target_Fusion', 'Unknown')
            print(f"[{i}/{len(candidates)}] Analyzing {target}...")
            
            safety_score, risks = self.check_safety(target)
            ai_res = self.analyze_ai(cand, risks)
            
            results.append({
                'target': target,
                'disease': cand.get('Associated_Disease', 'Unknown'),
                'count': cand.get('Patient_Count', 0),
                'risks': risks,
                'ai_verdict': ai_res.get('Verdict', 'N/A'),
                'ai_rationale': ai_res.get('Scientific_Rationale', ai_res.get('Rationale', 'Analyzing...')),
                'ai_strategy': ai_res.get('Screening_Strategy', 'N/A'),
                'ai_risks': ai_res.get('Risk_Factors', 'N/A'),
                'score': safety_score
            })
            
            # Rate limiting: wait between API calls to avoid hitting free tier limits
            if i < len(candidates):
                time.sleep(2)  # Increased delay for free tier (5-15 RPM limit)
            
        # Here we would call generate_dashboard(results) - see previous turn for HTML logic
        print(f"\n[SUCCESS] Analysis Complete. Processed {len(results)} candidates.")
        return results

if __name__ == "__main__":
    ExpertAgent("lead_candidates.csv", "prompts/expert_persona.txt").run()