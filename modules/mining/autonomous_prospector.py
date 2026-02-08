import time
import json
import os
import re
import sqlite3
import csv
import requests
import logging
import random
from datetime import datetime
import sys
from Bio import SeqIO, Entrez

# Ensure imports work for VPS environment
sys.path.append(os.getcwd())
try:
    from modules.mining.sra_scout import SRAScout
    from modules.mining.deep_miner_utils import DeepEngine, NeighborhoodWatch
    from modules.mining.full_orf_checks import full_orf_passes, get_full_orf_config
except ImportError:
    print("[!] Critical: Modules missing. Ensure 'sra_scout.py' and 'deep_miner_utils.py' exist.")
    sys.exit(1)

# --- CONFIGURATION ---
DB_PATH = "data/prospector.db"
LOG_FILE = "prospector.log"
Entrez.email = "founder@senarybio.com"

# Broad, diverse SRA search queries to maximize coverage and minimize overlap.
# Covers extreme environments, marine, terrestrial, gut, plant, industrial, symbionts.
BROAD_SEARCH_QUERIES = [
    "acid mine drainage metagenome",
    "soda lake metagenome",
    "serpentinizing spring metagenome",
    "geothermal hot spring metagenome",
    "Antarctic soil metagenome",
    "Arctic permafrost metagenome",
    "deep sea sediment metagenome",
    "hydrothermal vent metagenome",
    "cold seep metagenome",
    "mangrove sediment metagenome",
    "coral reef microbiome metagenome",
    "ocean gyre metagenome",
    "permafrost metagenome",
    "desert soil metagenome",
    "volcanic soil metagenome",
    "cave microbiome metagenome",
    "compost metagenome",
    "rumen metagenome",
    "termite gut metagenome",
    "insect gut metagenome",
    "fish gut microbiome metagenome",
    "rhizosphere metagenome",
    "phyllosphere metagenome",
    "plant endophyte metagenome",
    "wastewater treatment metagenome",
    "anaerobic digestor metagenome",
    "oil reservoir metagenome",
    "mine tailings metagenome",
    "hypersaline lake metagenome",
    "salt marsh sediment metagenome",
    "freshwater sediment metagenome",
    "sponge symbiont metagenome",
    "coral symbiont metagenome",
    "lichen microbiome metagenome",
    "biofilm metagenome",
    "activated sludge metagenome",
    "salt flat halophile metagenome",
    "brine pool metagenome",
    "sulfuric spring metagenome",
    "alkaline lake metagenome",
    "glacier ice metagenome",
    "deep subsurface metagenome",
    "peat bog metagenome",
    "rice paddy soil metagenome",
    "marine sponge metagenome",
    "hydrothermal plume metagenome",
    "wood decay metagenome",
    "soda ash lake metagenome",
    "saline spring metagenome",
    "iron mine metagenome",
    "salt cave metagenome",
    "thermal spring sediment metagenome",
]

# Very broad, generic queries used when stuck seeing same SRA files 2–3 times in a row.
# Searched without wgs[Prop] for a different slice of NCBI.
SUPER_BROAD_QUERIES = [
    "metagenome",
    "environmental metagenome",
    "shotgun metagenome",
    "whole genome shotgun metagenome",
    "uncultured microbiome",
    "microbial metagenome",
    "environmental DNA sequencing",
    "metagenomic assembly",
    "WGS metagenome",
    "marine metagenome",
    "soil metagenome",
    "gut metagenome",
]

# Setup Logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.FileHandler(LOG_FILE), logging.StreamHandler()]
)

class LLMClient:
    """Abstracts the AI backend. Optimized for Local Llama-3-Bio."""
    def __init__(self):
        self.provider = os.getenv("LLM_PROVIDER", "local").lower() 
        self.local_url = os.getenv("LLM_LOCAL_URL", "http://localhost:11434/api/chat")
        self.model_name = os.getenv("LLM_MODEL", "llama3-bio") 
        logging.info(f"Initialized LLM Backend: {self.provider} | Model: {self.model_name}")

    def generate(self, system_prompt, user_prompt):
        try:
            payload = {
                "model": self.model_name,
                "messages": [
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt}
                ],
                "temperature": 0.8, # Slightly higher temp for more creative search terms
                "stream": False
            }
            response = requests.post(self.local_url, json=payload, timeout=120)
            
            if response.status_code == 200:
                return response.json()['message']['content'].strip()
            else:
                logging.error(f"LLM Error {response.status_code}: {response.text}")
                return ""
        except Exception as e:
            logging.error(f"Local LLM Connection Failed: {e}")
            return ""

class AutonomousProspector:
    def __init__(self):
        self.llm = LLMClient()
        self.scout = SRAScout()
        
        logging.info("Initializing Deep Learning Engines (ESM-2 & CRISPR-Context)...")
        self.deep_engine = DeepEngine()       
        self.context_engine = NeighborhoodWatch() 
        
        self._setup_db()
        self.failure_streak = 0
        self._last_raw_ids = set()
        self._repeat_streak = 0

    def _setup_db(self):
        os.makedirs("data", exist_ok=True)
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute('''CREATE TABLE IF NOT EXISTS history 
                     (id INTEGER PRIMARY KEY, query TEXT, hits INTEGER, strategy TEXT, timestamp DATETIME)''')
        c.execute('''CREATE TABLE IF NOT EXISTS visited_ids 
                     (nucleotide_id TEXT PRIMARY KEY, timestamp DATETIME)''')
        c.execute('''CREATE TABLE IF NOT EXISTS query_cycle 
                     (id INTEGER PRIMARY KEY CHECK (id=1), idx INTEGER)''')
        c.execute('INSERT OR IGNORE INTO query_cycle (id, idx) VALUES (1, 0)')
        conn.commit()
        conn.close()

    def get_next_broad_query(self):
        """Cycle through BROAD_SEARCH_QUERIES to maximize diversity and minimize overlap."""
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute("SELECT idx FROM query_cycle WHERE id=1")
        row = c.fetchone()
        idx = (row[0] if row else 0) % len(BROAD_SEARCH_QUERIES)
        query = BROAD_SEARCH_QUERIES[idx]
        c.execute("UPDATE query_cycle SET idx = ? WHERE id=1", ((idx + 1) % len(BROAD_SEARCH_QUERIES),))
        conn.commit()
        conn.close()
        strategy = f"Broad[{idx + 1}/{len(BROAD_SEARCH_QUERIES)}]: {query[:40]}..."
        return query, strategy

    def get_random_super_broad_query(self):
        """Pick a random very-broad query to break out of repeated same-SRA loops."""
        query = random.choice(SUPER_BROAD_QUERIES)
        strategy = f"VeryBroad(random): {query[:40]}..."
        return query, strategy

    def get_visited_ids(self):
        """Return set of nucleotide IDs already mined to avoid re-processing."""
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute("SELECT nucleotide_id FROM visited_ids")
        ids = {row[0] for row in c.fetchall()}
        conn.close()
        return ids

    def mark_ids_visited(self, id_list):
        """Persist nucleotide IDs as visited after processing."""
        if not id_list:
            return
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        now = datetime.now()
        c.executemany(
            "INSERT OR IGNORE INTO visited_ids (nucleotide_id, timestamp) VALUES (?, ?)",
            [(str(uid), now) for uid in id_list]
        )
        conn.commit()
        conn.close()

    def log_history(self, query, hits, strategy):
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute("INSERT INTO history (query, hits, strategy, timestamp) VALUES (?, ?, ?, ?)",
                  (query, hits, strategy, datetime.now()))
        conn.commit()
        conn.close()

    def get_recent_history(self, limit=5):
        conn = sqlite3.connect(DB_PATH)
        c = conn.cursor()
        c.execute("SELECT query, hits, strategy FROM history ORDER BY id DESC LIMIT ?", (limit,))
        data = c.fetchall()
        conn.close()
        return "\n".join([f"- Strategy: {row[2]} | Query: '{row[0]}' | Hits: {row[1]}" for row in data])

    def fetch_metadata(self, id_list):
        if not id_list: return {}
        try:
            handle = Entrez.esummary(db="nucleotide", id=",".join(str(x) for x in id_list))
            records = Entrez.read(handle)
            handle.close()
            return {str(r['Id']): r.get('Title', r.get('Description', 'Unknown')) for r in records}
        except Exception as e:
            logging.error(f"Metadata fetch failed: {e}")
            return {}

    def semantic_filter(self, metadata_dict):
        if not metadata_dict: return []

        deep_mine_max = int(os.getenv("DEEP_MINE_MAX", "15"))
        meta_keys = list(metadata_dict.keys())
        descriptions = "\n".join([f"ID {k}: {v}" for k, v in list(metadata_dict.items())[:20]])

        system = "You are an expert computational biologist specializing in metagenomics."
        prompt = f"""
        We are mining for **Type VI CRISPR-Cas13d** systems.
        Review these dataset descriptions. Select the Top {min(deep_mine_max, 20)} that are most likely to contain **uncultured CRISPR loci**, **viral defense islands**, or **diverse microbial communities**.

        Look for high-value terms: 'microbial mat', 'biofilm', 'hot spring', 'hypersaline', 'uncultured', 'phage defense', 'rumen', 'permafrost', 'deep sea'.
        Avoid: 'human', 'clinical', 'mitochondrion', 'chloroplast'.

        DATASETS:
        {descriptions}

        OUTPUT FORMAT:
        Return ONLY a JSON list of IDs. Do not write anything else.
        Example: ["12345", "67890"]
        """

        response = self.llm.generate(system, prompt)
        selected_ids = []

        try:
            if response and '[' in response and ']' in response:
                clean_json = response[response.find('['):response.rfind(']')+1]
                parsed = json.loads(clean_json)
                selected_ids = [str(i) for i in parsed if str(i) in metadata_dict]
        except (json.JSONDecodeError, TypeError):
            pass

        if not selected_ids and response:
            raw_ids = re.findall(r'\d+', response)
            selected_ids = [uid for uid in raw_ids if uid in metadata_dict]
            if selected_ids:
                logging.info("Semantic filter: extracted IDs via regex fallback.")

        if len(selected_ids) < 3:
            fallback_count = min(deep_mine_max, len(meta_keys))
            selected_ids = meta_keys[:fallback_count]
            logging.warning(f"Semantic filter returned {len(selected_ids) if selected_ids else 0} valid IDs. Using first {fallback_count} from metadata.")

        return selected_ids[:deep_mine_max]

    def formulate_strategy(self, use_broad_list=False):
        """
        The 'Ralph' Loop: Adapts strategy based on success/failure.
        When use_broad_list=True (e.g. all IDs visited or high failure streak),
        cycles through BROAD_SEARCH_QUERIES to maximize diversity and minimize overlap.
        """
        if use_broad_list:
            return self.get_next_broad_query()

        history = self.get_recent_history()

        system = "You are the Autonomous Director of Discovery. Plan the next NCBI mining operation."
        prompt = f"""
        GOAL: Discover novel Cas13d enzymes in NCBI Whole Genome Shotgun (WGS) data.

        STATUS:
        - Failure Streak: {self.failure_streak}
        - Recent History:
        {history}

        TASK:
        Generate the next specific NCBI Query.
        1. If Failure Streak > 2: PIVOT to a COMPLETELY DIFFERENT environment (not similar to recent).
        2. If succeeding: REFINE the current query.
        3. DIVERSITY IS CRITICAL. Pick environments NOT in recent history.
        4. Avoid repeating: hypersaline, salt flat, sponge symbiont if already tried recently.

        DIVERSE ENVIRONMENT IDEAS (pick one not recently used):
        - Acid mine, soda lake, serpentinizing spring, geothermal
        - Antarctic/Arctic, permafrost, glacier, deep subsurface
        - Cold seep, mangrove, coral reef, ocean gyre
        - Termite gut, fish gut, rhizosphere, phyllosphere
        - Wastewater, oil reservoir, mine tailings, activated sludge
        - Lichen, peat bog, rice paddy, wood decay, hydrothermal plume

        OUTPUT FORMAT:
        JSON with keys: "strategy_name", "query".
        Return ONLY environmental keywords for "query" - no filters like wgs[Prop].
        Example: {{ "strategy_name": "Pivot to Acid Mine", "query": "acid mine drainage metagenome" }}
        """

        try:
            response = self.llm.generate(system, prompt)
            clean_json = response[response.find('{'):response.rfind('}')+1]
            data = json.loads(clean_json)
            return data.get('query'), data.get('strategy_name')
        except Exception:
            return self.get_next_broad_query()

    def deep_mine(self, id_list):
        require_crispr = os.getenv("REQUIRE_CRISPR", "1").lower() in ("1", "true", "yes")
        require_full_structure = os.getenv("REQUIRE_FULL_STRUCTURE", "0").lower() in ("1", "true", "yes")
        min_repeat_count = int(os.getenv("MIN_REPEAT_COUNT", "1"))
        esm_threshold = float(os.getenv("ESM_THRESHOLD", "0.75"))
        use_diversity_band = os.getenv("ESM_SIMILARITY_CEILING", "").strip() != ""
        full_orf_cfg = get_full_orf_config()
        hits = []
        for ncbi_id in id_list:
            logging.info(f"   -> Scanning ID {ncbi_id}...")
            contigs_scanned = 0
            contigs_with_crispr = 0
            orfs_scored = 0
            best_score = 0.0
            try:
                handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
                records = list(SeqIO.parse(handle, "fasta"))
                handle.close()
                for record in records:
                    dna_seq = record.seq
                    contigs_scanned += 1
                    has_crispr = self.context_engine.has_crispr_array(dna_seq)
                    if has_crispr:
                        contigs_with_crispr += 1
                    if require_crispr and not has_crispr:
                        continue
                    if has_crispr:
                        logging.info(f"      [!] CRISPR Array detected in {ncbi_id}. Analyzing proteins...")

                    repeat_domains = self.context_engine.get_repeat_domains(dna_seq)
                    contig_len = len(dna_seq)

                    frames = [dna_seq[i:].translate(to_stop=False) for i in range(3)]
                    frames += [dna_seq.reverse_complement()[i:].translate(to_stop=False) for i in range(3)]
                    for frame_index, frame in enumerate(frames):
                        orfs = str(frame).split("*")
                        for orf_index, orf in enumerate(orfs):
                            if 600 <= len(orf) <= 1400:
                                if not full_orf_passes(
                                    orf,
                                    contig_len,
                                    frame_index,
                                    orf_index,
                                    orfs,
                                    require_m=full_orf_cfg["require_m"],
                                    min_tail=full_orf_cfg["min_tail"],
                                    boundary_margin=full_orf_cfg["boundary_margin"],
                                ):
                                    continue
                                score = self.deep_engine.score_candidate(orf)
                                orfs_scored += 1
                                if score > best_score:
                                    best_score = score
                                if score > esm_threshold and self.deep_engine.passes_diversity_band(score):
                                    if (require_full_structure or require_crispr) and len(repeat_domains) < min_repeat_count:
                                        continue
                                    logging.info(f"      [***] DISCOVERY: Novel Enzyme (Score: {score:.4f})")
                                    hits.append((f"{ncbi_id}_ORF", orf, score, ncbi_id, repeat_domains))

                if contigs_scanned > 0:
                    logging.info(f"      ID {ncbi_id}: {contigs_scanned} contigs, {contigs_with_crispr} with CRISPR, {orfs_scored} ORFs scored, best={best_score:.3f}")
                if orfs_scored > 0 and best_score < esm_threshold:
                    logging.info(f"      (Best score {best_score:.3f} below threshold {esm_threshold})")
            except Exception as e:
                logging.error(f"Error mining {ncbi_id}: {e}")

        return hits

    def run_loop(self):
        require_crispr = os.getenv("REQUIRE_CRISPR", "1").lower() in ("1", "true", "yes")
        esm_threshold = float(os.getenv("ESM_THRESHOLD", "0.75"))
        deep_mine_max = int(os.getenv("DEEP_MINE_MAX", "15"))
        logging.info(f"--- Senary Bio: Deep Prospector (Model: {self.llm.model_name}) ---")
        diversity_ceiling = os.getenv("ESM_SIMILARITY_CEILING", "")
        require_full_structure = os.getenv("REQUIRE_FULL_STRUCTURE", "0").lower() in ("1", "true", "yes")
        min_repeat_count = os.getenv("MIN_REPEAT_COUNT", "1")
        esm_ref_fasta = os.getenv("ESM_REFERENCE_FASTA", "").strip()
        logging.info(f"Config: DEEP_MINE_MAX={deep_mine_max}, ESM_THRESHOLD={esm_threshold}, REQUIRE_CRISPR={require_crispr}, DIVERSITY_CEILING={diversity_ceiling or 'off'}")
        if esm_ref_fasta:
            logging.info(f"ESM_REFERENCE_FASTA set: mining for closest similarity to RfxCas13d/PspCas13a (leave ESM_SIMILARITY_CEILING unset to keep closest matches)")
        logging.info(f"Full enzyme: REQUIRE_START_M={os.getenv('REQUIRE_START_M','1')}, MIN_CTERM_TAIL={os.getenv('MIN_CTERM_TAIL','15')}, REQUIRE_FULL_STRUCTURE={require_full_structure}, MIN_REPEAT_COUNT={min_repeat_count}")
        logging.info(f"ORF size: 600-1400 aa | Broad queries: {len(BROAD_SEARCH_QUERIES)} | SRA_MAX_RECORDS={os.getenv('SRA_MAX_RECORDS', '100')}")

        # Threshold: same/similar SRA set this many times in a row -> force very broad random search
        repeat_streak_threshold = int(os.getenv("REPEAT_STREAK_THRESHOLD", "2"))
        stale_overlap_ratio = float(os.getenv("STALE_OVERLAP_RATIO", "0.6"))
        stale_new_ids_max = int(os.getenv("STALE_NEW_IDS_MAX", "5"))
        pagination_max_pages = int(os.getenv("SRA_PAGINATION_MAX_PAGES", "10"))

        while True:
            visited = self.get_visited_ids()
            broad_threshold = int(os.getenv("BROAD_LIST_THRESHOLD", "500"))
            # Stuck seeing same SRA 2–3 times in a row -> one shot of very broad random search (no wgs[Prop])
            use_very_broad = self._repeat_streak >= repeat_streak_threshold
            if use_very_broad:
                self._repeat_streak = 0
                logging.info("Repeat streak reached: forcing very broad random search to break loop.")
            # Use broad diverse list when: failure streak, or already visited many IDs (reduces overlap)
            use_broad = self.failure_streak >= 1 or len(visited) >= broad_threshold
            if use_very_broad:
                query, strategy = self.get_random_super_broad_query()
                logging.info(f"PLAN: {strategy} | Query: {query}")
                max_records = min(200, int(os.getenv("SRA_MAX_RECORDS", "100")) * 2)
            else:
                query, strategy = self.formulate_strategy(use_broad_list=use_broad)
                logging.info(f"PLAN: {strategy} | Query: {query}")
                max_records = int(os.getenv("SRA_MAX_RECORDS", "100"))

            # Search with pagination: if first page yields only visited IDs, try deeper pages
            raw_ids = []
            ids = []
            for page in range(pagination_max_pages):
                retstart = page * max_records
                if use_very_broad:
                    raw_ids = self.scout.search_very_broad(query, max_records=max_records, retstart=retstart)
                else:
                    raw_ids = self.scout.search_wgs(query, max_records=max_records, retstart=retstart)

                ids = [uid for uid in raw_ids if uid not in visited]
                if visited:
                    logging.info(f"Skipping {len(visited)} previously visited IDs. {len(ids)} new candidates (page {page + 1}).")

                if ids:
                    break
                if len(raw_ids) < max_records:
                    logging.info("Reached end of results (fewer than max_records). No more pages.")
                    break
                logging.info(f"Page {page + 1} all visited. Trying next page...")
                time.sleep(2)

            # Detect "same SRA over and over": high overlap with last iteration or very few new IDs
            if raw_ids:
                curr_set = set(raw_ids)
                overlap = len(curr_set & self._last_raw_ids) / len(curr_set) if curr_set else 0
                if overlap >= stale_overlap_ratio or len(ids) <= stale_new_ids_max:
                    self._repeat_streak += 1
                    logging.info(f"Stale iteration (overlap={overlap:.2f}, new={len(ids)}). Repeat streak={self._repeat_streak}.")
                else:
                    self._repeat_streak = 0
                self._last_raw_ids = curr_set

            if not ids:
                logging.warning("All IDs already visited or none found. Using broad list next iteration.")
                self.failure_streak += 1
                time.sleep(5)
                continue

            logging.info(f"Analyzing {len(ids)} candidates with {self.llm.model_name}...")
            metadata = self.fetch_metadata(ids)
            best_ids = self.semantic_filter(metadata)
            
            logging.info(f"FILTER: Selected {len(best_ids)} high-probability targets.")

            if not best_ids:
                logging.warning("No targets to mine (semantic filter returned empty). Skipping iteration.")

            hits = 0
            if best_ids:
                new_discoveries = self.deep_mine(best_ids)
                if new_discoveries:
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    fasta_path = f"data/raw_sequences/deep_hits_{timestamp}.fasta"
                    meta_path = f"data/raw_sequences/deep_hits_{timestamp}_metadata.csv"
                    with open(fasta_path, "w") as f:
                        for name, seq, score, sra_accession, repeat_domains in new_discoveries:
                            seq_id = f"{name}_Score_{score:.3f}"
                            f.write(f">{seq_id}\n{seq}\n")
                    with open(meta_path, "w", newline="", encoding="utf-8") as m:
                        writer = csv.writer(m)
                        writer.writerow(["sequence_id", "sra_accession", "repeat_domains", "score"])
                        for name, seq, score, sra_accession, repeat_domains in new_discoveries:
                            seq_id = f"{name}_Score_{score:.3f}"
                            repeat_str = "|".join(repeat_domains) if repeat_domains else ""
                            writer.writerow([seq_id, sra_accession, repeat_str, f"{score:.3f}"])
                    logging.info(f"SUCCESS: Saved {len(new_discoveries)} Deep Learning Hits to {fasta_path} and {meta_path}.")
                    hits = len(new_discoveries)
                    self.failure_streak = 0
                else:
                    self.failure_streak += 1
            
            self.log_history(query, hits, strategy)

            # Mark all IDs from this search as visited to avoid re-processing
            self.mark_ids_visited(raw_ids)

            logging.info("Sleeping for 30s...")
            time.sleep(30)

if __name__ == "__main__":
    prospector = AutonomousProspector()
    prospector.run_loop()