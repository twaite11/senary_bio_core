# Walkthrough: Screen → Autonomous Prospector

From installing **GNU Screen** to running the autonomous prospector in a detached session. Use this on a Linux VPS (e.g. Ubuntu) or WSL.

---

## 1. Install Screen

**Ubuntu / Debian:**
```bash
sudo apt-get update
sudo apt-get install -y screen
```

**macOS (Homebrew):**
```bash
brew install screen
```

**Check it’s installed:**
```bash
screen --version
```

---

## 2. One-time project setup (if not done yet)

From your project root (e.g. `~/collateral_bio_core` or `/opt/collateral_bio_core`):

```bash
cd /path/to/collateral_bio_core

# Python venv
python3.10 -m venv venv
source venv/bin/activate   # Linux/macOS

# Dependencies (includes torch, transformers, biopython for prospector)
pip install -r requirements.txt

# Config (optional but recommended)
cp config/pipeline.env.example .env
# Edit .env if you want: ESM_THRESHOLD, PROSPECTOR_WORKERS, EXHAUSTIVE_MODE, etc.
```

**Optional for LLM-driven strategy:** Install and run Ollama (or set `LLM_PROVIDER` / `LLM_LOCAL_URL` in `.env`). The prospector can run without an LLM in **exhaustive mode** (`EXHAUSTIVE_MODE=1`).

---

## 3. Start a Screen session and run the prospector

**Create a named session (e.g. `prospector`):**
```bash
screen -S prospector
```

You’re now inside that session. Your shell prompt is unchanged; the window is managed by Screen.

**Activate venv and start the prospector:**
```bash
cd /path/to/collateral_bio_core
source venv/bin/activate

# Optional env vars (or set in .env)
export ESM_SIMILARITY_CEILING=0.82
# export EXHAUSTIVE_MODE=1
# export PROSPECTOR_WORKERS=4

python modules/mining/autonomous_prospector.py
```

The prospector will log to the terminal and to `prospector.log`. Outputs go to `data/raw_sequences/deep_hits_*.fasta` and `*_metadata.csv`.

**Detach from the session (leave it running in the background):**  
Press:

```
Ctrl+A, then D
```

(Hold **Ctrl**, press **A**, release, then press **D**.)

You’ll see something like `[detached from 12345.prospector]`. The prospector keeps running inside that session.

---

## 4. Reattach to the same session later

**List sessions:**
```bash
screen -ls
```

Example:
```
12345.prospector    (Detached)
```

**Reattach by name:**
```bash
screen -r prospector
```

Or by PID if the name is unique:
```bash
screen -r 12345
```

You’ll see the same terminal and logs. To detach again: **Ctrl+A**, then **D**.

---

## 5. Useful Screen commands (quick reference)

| Action | Keys / Command |
|--------|-----------------|
| Detach (leave running) | **Ctrl+A**, then **D** |
| List sessions | `screen -ls` |
| Reattach by name | `screen -r prospector` |
| Kill a session (from inside) | **Ctrl+A**, then **K** (confirm with **y**) |
| Kill a session from outside | `screen -S prospector -X quit` |

---

## 6. Optional: run in exhaustive mode with workers

To mine all bacteria/archaea WGS contigs with multiple workers (good for a 32‑core box):

```bash
screen -S prospector
cd /path/to/collateral_bio_core
source venv/bin/activate

export EXHAUSTIVE_MODE=1
export FILTER_23_HEPN=1
export PROSPECTOR_WORKERS=4
export ESM_BATCH_SIZE=32

python modules/mining/autonomous_prospector.py
```

Then detach with **Ctrl+A**, **D**.

---

## 7. Check that the prospector is running and writing output

From another terminal (no need to attach to Screen):

```bash
# Latest log lines
tail -f /path/to/collateral_bio_core/prospector.log

# New FASTA/metadata files
ls -la /path/to/collateral_bio_core/data/raw_sequences/deep_hits_*
```

---

## End-to-end checklist

1. Install Screen: `sudo apt-get install -y screen` (or `brew install screen` on macOS).
2. Clone repo, create venv, `pip install -r requirements.txt`, optionally `cp config/pipeline.env.example .env`.
3. Start session: `screen -S prospector`.
4. In the session: `cd` to repo, `source venv/bin/activate`, set env vars if needed, run `python modules/mining/autonomous_prospector.py`.
5. Detach: **Ctrl+A**, then **D**.
6. Later: `screen -r prospector` to reattach; **Ctrl+A**, **D** to detach again.
