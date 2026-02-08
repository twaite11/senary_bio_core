"""
Generate Structure Dashboard HTML: 3D viewer + domain coloring + homology table.
Reads: input FASTA, InterPro domains.json, homology_scores.json, OmegaFold PDB paths.
Optional: synthesis_metadata.csv or family_manifest for SRA + CRISPR repeats.
Predicts NUC/REC lobes for structure display.
Serve the HTML via HTTP (e.g. python -m http.server) to load PDBs.
"""
import argparse
import csv
import json
import re
from pathlib import Path
from typing import List, Tuple

from Bio import SeqIO

HEPN_REGEX = re.compile(r"R.{4,6}H")

# Pfam: HEPN = nuclease lobe; HEL/WYL = recognition lobe
NUC_SIGNATURES = {"PF05168"}  # HEPN
REC_SIGNATURES = {"PF01228", "PF18456"}  # HEL, WYL


def count_hepn(seq_str: str) -> int:
    return len(HEPN_REGEX.findall(str(seq_str)))


def get_hepn_positions(seq_str: str) -> List[Tuple[int, int]]:
    """Return (start, end) 0-based for each HEPN motif."""
    return [(m.start(), m.end()) for m in HEPN_REGEX.finditer(str(seq_str))]


def predict_nuc_rec_spans(
    sequence: str,
    domains: list,
    hepn_padding: int = 50,
) -> dict:
    """
    Predict NUC (nuclease) and REC (recognition) lobe residue spans.
    NUC = HEPN domains (PF05168) or expanded HEPN motifs if no InterPro.
    REC = HEL/WYL domains or remaining regions.
    Returns {"nuc_spans": [(start,end),...], "rec_spans": [(start,end),...], "nuc_summary": "100-250, 400-550"}.
    """
    nuc_spans = []
    rec_spans = []
    seq_len = len(sequence)

    if domains:
        for d in domains:
            sig = (d.get("signature") or "").split(".")[0]
            start = int(d.get("start", 0))
            end = int(d.get("end", 0))
            if sig in NUC_SIGNATURES:
                nuc_spans.append((start, end))
            elif sig in REC_SIGNATURES:
                rec_spans.append((start, end))

    if not nuc_spans:
        # Fallback: expand HEPN motifs ±hepn_padding aa
        for start, end in get_hepn_positions(sequence):
            s = max(0, start - hepn_padding)
            e = min(seq_len, end + hepn_padding)
            nuc_spans.append((s, e))

    # Merge overlapping spans and sort
    def merge_spans(spans: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        if not spans:
            return []
        spans = sorted(spans, key=lambda x: x[0])
        merged = [spans[0]]
        for s, e in spans[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        return merged

    nuc_spans = merge_spans(nuc_spans)
    rec_spans = merge_spans(rec_spans)

    # REC = regions not in NUC (if we have InterPro REC domains, use them; else infer)
    if not rec_spans and nuc_spans:
        nuc_set = set()
        for s, e in nuc_spans:
            nuc_set.update(range(s, e))
        pos = 0
        while pos < seq_len:
            if pos not in nuc_set:
                end = pos
                while end < seq_len and end not in nuc_set:
                    end += 1
                rec_spans.append((pos, end))
                pos = end
            else:
                pos += 1
        rec_spans = merge_spans(rec_spans)

    nuc_summary = ", ".join(f"{s+1}-{e}" for s, e in nuc_spans) if nuc_spans else "-"
    rec_summary = ", ".join(f"{s+1}-{e}" for s, e in rec_spans) if rec_spans else "-"

    return {
        "nuc_spans": nuc_spans,
        "rec_spans": rec_spans,
        "nuc_summary": nuc_summary,
        "rec_summary": rec_summary,
    }


def parse_family_id(seq_id: str) -> str:
    parts = str(seq_id).split("_")
    return parts[0] if parts else seq_id


def load_synthesis_metadata(proj: Path, metadata_path: str = None) -> dict:
    """Load new_id -> {sra_accession, repeat_domains} from synthesis_metadata or family_manifest."""
    lookup = {}
    paths = []
    if metadata_path:
        p = proj / metadata_path if not Path(metadata_path).is_absolute() else Path(metadata_path)
        if p.exists():
            paths.append(p)
    if not paths:
        mined = proj / "data" / "mined_sequences"
        if mined.exists():
            paths.extend(sorted(mined.glob("synthesis_metadata.csv")))
            paths.extend(sorted(mined.glob("family_manifest_*.csv"))[-1:])  # latest manifest
    for p in paths:
        try:
            with open(p, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    nid = row.get("new_id", "").strip()
                    if nid:
                        lookup[nid] = {
                            "sra_accession": row.get("sra_accession", "").strip(),
                            "repeat_domains": row.get("repeat_domains", "").strip(),
                        }
        except Exception as e:
            print(f"[!] Warning: could not read metadata {p}: {e}")
    return lookup


def find_pdb_paths(structures_dir: Path, proj: Path) -> dict:
    """Scan OmegaFold (flat *.pdb) or ColabFold (subdirs) output for {seq_id: relative_path}."""
    result = {}
    # OmegaFold: flat *.pdb in output dir
    for pdb in structures_dir.glob("*.pdb"):
        rel = (structures_dir / pdb.name).relative_to(proj)
        result[pdb.stem] = str(rel).replace("\\", "/")
    # ColabFold fallback: subdirs with *_rank_001*.pdb
    if not result:
        for sub in structures_dir.iterdir():
            if not sub.is_dir():
                continue
            for pdb in sub.glob("*_rank_001*.pdb"):
                name = pdb.stem.replace("_unrelaxed", "").replace("_relaxed", "").split("_rank")[0]
                rel = pdb.relative_to(proj)
                result[name] = str(rel).replace("\\", "/")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Generate structure dashboard with 3D viewer and homology table."
    )
    parser.add_argument(
        "--input",
        default="data/structure_pipeline/passed_structures.fasta",
        help="Input FASTA (pipeline: passed_structures.fasta; or input_2-3_hepn.fasta)",
    )
    parser.add_argument(
        "--domains",
        default="data/structure_pipeline/interpro/domains.json",
        help="InterPro domains JSON",
    )
    parser.add_argument(
        "--homology",
        default="data/structure_pipeline/homology_scores.json",
        help="Homology scores JSON (from pipeline or run_tmscore)",
    )
    parser.add_argument(
        "--filter-results",
        default="data/structure_pipeline/structure_filter_results.json",
        help="Fallback: structure filter results (has cas13a/cas13b/cas13d)",
    )
    parser.add_argument(
        "--structures-dir",
        default="data/structure_pipeline/structures/omegafold",
        help="OmegaFold output directory",
    )
    parser.add_argument(
        "--metadata",
        default=None,
        help="synthesis_metadata.csv or family_manifest CSV (new_id, sra_accession, repeat_domains); searches mined_sequences/ if not set",
    )
    parser.add_argument(
        "--output",
        default="visualization/structure_dashboard.html",
        help="Output HTML path",
    )
    args = parser.parse_args()

    proj = Path(__file__).resolve().parent.parent
    input_path = proj / args.input
    domains_path = proj / args.domains
    homology_path = proj / args.homology
    filter_results_path = proj / args.filter_results
    struct_dir = proj / args.structures_dir
    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = proj / out_path

    # Load sequences
    seqs = {}
    if input_path.exists():
        for rec in SeqIO.parse(str(input_path), "fasta"):
            seq_str = str(rec.seq)
            seqs[rec.id] = {
                "sequence": seq_str,
                "length": len(seq_str),
                "hepn_count": count_hepn(seq_str),
                "family": parse_family_id(rec.id),
            }
    else:
        seqs = {}

    # Load synthesis metadata (SRA, CRISPR repeats)
    synth_meta = load_synthesis_metadata(proj, args.metadata)

    # Load domains
    domains = {}
    if domains_path.exists():
        with open(domains_path) as f:
            domains = json.load(f)

    # Load homology: pipeline writes homology_scores.json; fallback to structure_filter_results (has cas13a/cas13b/cas13d)
    homology = {}
    if homology_path.exists():
        with open(homology_path) as f:
            homology = json.load(f)
    if not homology and filter_results_path.exists():
        with open(filter_results_path) as f:
            raw = json.load(f)
        for sid, r in raw.items():
            h = {k: r[k] for k in ("cas13a", "cas13b", "cas13d") if k in r and r[k] is not None}
            if h:
                homology[sid] = h

    # PDB paths (relative for fetch from project root)
    pdb_paths = find_pdb_paths(struct_dir, proj) if struct_dir.exists() else {}

    # Build rows for table
    all_ids = sorted(set(seqs.keys()) | set(domains.keys()) | set(homology.keys()) | set(pdb_paths.keys()))
    rows = []
    for seq_id in all_ids:
        s = seqs.get(seq_id, {})
        d = domains.get(seq_id, [])
        h = homology.get(seq_id, {})
        cas13a = h.get("cas13a")
        cas13b = h.get("cas13b")
        cas13d = h.get("cas13d")
        best = None
        best_name = None
        if cas13a is not None and (best is None or (cas13a or 0) > best):
            best = cas13a
            best_name = "Cas13a"
        if cas13b is not None and (best is None or (cas13b or 0) > best):
            best = cas13b
            best_name = "Cas13b"
        if cas13d is not None and (best is None or (cas13d or 0) > best):
            best = cas13d
            best_name = "RfxCas13d"

        # NUC/REC lobe prediction
        seq_str = s.get("sequence", "")
        nuc_rec = predict_nuc_rec_spans(seq_str, d) if seq_str else {}

        meta = synth_meta.get(seq_id, {})
        rows.append({
            "id": seq_id,
            "length": s.get("length", 0),
            "hepn_count": s.get("hepn_count", 0),
            "family": s.get("family", parse_family_id(seq_id)),
            "domains": d,
            "cas13a": round(cas13a * 100, 1) if cas13a is not None else None,
            "cas13b": round(cas13b * 100, 1) if cas13b is not None else None,
            "cas13d": round(cas13d * 100, 1) if cas13d is not None else None,
            "best_match": best_name,
            "best_pct": round((best or 0) * 100, 1) if best is not None else None,
            "pdb_path": pdb_paths.get(seq_id),
            "sra_accession": meta.get("sra_accession", ""),
            "repeat_domains": meta.get("repeat_domains", ""),
            "nuc_summary": nuc_rec.get("nuc_summary", "-"),
            "rec_summary": nuc_rec.get("rec_summary", "-"),
            "nuc_spans": nuc_rec.get("nuc_spans", []),
            "rec_spans": nuc_rec.get("rec_spans", []),
        })

    # Sort by best homology desc
    rows.sort(key=lambda r: (r["best_pct"] or 0, r["id"]), reverse=True)

    # Domain color map for 3Dmol; NUC/REC for lobe coloring
    domain_colors = {
        "PF05168": "#22d3ee",   # HEPN (nuclease)
        "PF01228": "#f97316",   # HEL
        "PF18456": "#94a3b8",   # WYL
    }
    lobe_colors = {"nuc": "#22d3ee", "rec": "#f97316"}  # NUC=cyan, REC=orange

    data_json = json.dumps({
        "rows": rows,
        "domain_colors": domain_colors,
        "lobe_colors": lobe_colors,
    })

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Structure Dashboard - 2-3 HEPN Cas13</title>
  <link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&display=swap" rel="stylesheet">
  <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
  <style>
    :root {{ --bg:#0f172a; --card:#1e3a5f; --accent:#38bdf8; --hepn:#22d3ee; --text:#f8fafc; --muted:#94a3b8; }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; padding: 1.5rem; font-family: 'DM Sans', system-ui; background: var(--bg); color: var(--text); min-height: 100vh; }}
    h1 {{ font-size: 1.5rem; margin: 0 0 0.5rem; }}
    .sub {{ color: var(--muted); font-size: 0.9rem; margin-bottom: 1rem; }}
    .layout {{ display: grid; grid-template-columns: 1fr 400px; gap: 1.5rem; }}
    @media (max-width: 900px) {{ .layout {{ grid-template-columns: 1fr; }} }}
    .viewer {{ width: 100%; height: 480px; background: #1e293b; border-radius: 8px; }}
    .table-wrap {{ overflow-x: auto; max-height: 70vh; overflow-y: auto; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.85rem; }}
    th, td {{ padding: 0.4rem 0.6rem; text-align: left; border-bottom: 1px solid #334155; }}
    th {{ color: var(--muted); font-weight: 500; position: sticky; top: 0; background: var(--bg); }}
    tr:hover {{ background: #1e3a5f; }}
    tr.selected {{ background: #1e3a5f; outline: 1px solid var(--accent); }}
    .pct {{ font-variant-numeric: tabular-nums; }}
    .no-pdb {{ color: var(--muted); font-style: italic; }}
    .load-hint {{ color: var(--muted); font-size: 0.8rem; margin-top: 0.5rem; }}
  </style>
</head>
<body>
  <h1>Structure Dashboard – 2-3 HEPN Cas13</h1>
  <p class="sub">Folding homology vs Cas13a/13b · Domain-colored 3D structure · Serve via HTTP to load PDBs</p>

  <div class="layout">
    <div>
      <div id="molviewer" class="viewer"></div>
      <p class="load-hint">Select a row to load structure. Serve from project root: python -m http.server 8000</p>
      <div id="domain-legend" style="margin-top:0.5rem;font-size:0.85rem;color:var(--muted)"></div>
      <div id="detail-panel" style="margin-top:0.75rem;padding:0.5rem;background:#1e293b;border-radius:6px;font-size:0.8rem;color:var(--muted)"></div>
    </div>
    <div class="table-wrap">
      <table>
        <thead>
          <tr>
            <th>ID</th>
            <th>Family</th>
            <th>Len</th>
            <th>HEPN</th>
            <th>NUC</th>
            <th>REC</th>
            <th>Cas13a %</th>
            <th>Cas13b %</th>
            <th>RfxCas13d %</th>
            <th>Best</th>
            <th>SRA</th>
          </tr>
        </thead>
        <tbody id="tbody"></tbody>
      </table>
    </div>
  </div>

  <script>
    const DATA = {data_json};

    function renderTable() {{
      const tbody = document.getElementById("tbody");
      tbody.innerHTML = "";
      DATA.rows.forEach((r, i) => {{
        const tr = document.createElement("tr");
        tr.dataset.id = r.id;
        tr.dataset.index = i;
        const sra = (r.sra_accession || "").slice(0, 12);
        tr.innerHTML = `<td>${{r.id}}</td><td>${{r.family}}</td><td>${{r.length}}</td><td>${{r.hepn_count}}</td>
          <td title="${{r.nuc_summary || ''}}">${{(r.nuc_summary || '-').slice(0,20)}}${{(r.nuc_summary || '').length > 20 ? '…' : ''}}</td>
          <td title="${{r.rec_summary || ''}}">${{(r.rec_summary || '-').slice(0,20)}}${{(r.rec_summary || '').length > 20 ? '…' : ''}}</td>
          <td class="pct">${{r.cas13a != null ? r.cas13a + '%' : '-'}}</td>
          <td class="pct">${{r.cas13b != null ? r.cas13b + '%' : '-'}}</td>
          <td class="pct">${{r.cas13d != null ? r.cas13d + '%' : '-'}}</td>
          <td>${{r.best_match || '-'}}</td>
          <td title="${{r.sra_accession || ''}}">${{sra || '-'}}</td>`;
        tr.addEventListener("click", () => selectRow(i));
        tbody.appendChild(tr);
      }});
    }}

    function selectRow(index) {{
      document.querySelectorAll("#tbody tr").forEach(tr => tr.classList.remove("selected"));
      const tr = document.querySelector(`#tbody tr[data-index="${{index}}"]`);
      if (tr) tr.classList.add("selected");
      const row = DATA.rows[index];
      loadStructure(row);
      const dp = document.getElementById("detail-panel");
      const parts = [];
      if (row.sra_accession) parts.push(`<b>SRA:</b> ${{row.sra_accession}}`);
      if (row.repeat_domains) parts.push(`<b>CRISPR repeats:</b> ${{row.repeat_domains}}`);
      if (parts.length) dp.innerHTML = parts.join(" &nbsp;|&nbsp; ");
      else dp.innerHTML = "";
    }}

    function loadStructure(row) {{
      const viewer = document.getElementById("molviewer");
      const div = document.createElement("div");
      div.className = "viewer";
      div.style.width = "100%";
      div.style.height = "480px";
      viewer.parentNode.replaceChild(div, viewer);
      div.id = "molviewer";

      const view = $3Dmol.createViewer(div, {{ backgroundColor: "#1e293b" }});

      if (!row.pdb_path) {{
        view.addLabel("No PDB available", {{ position: {{ x: 0, y: 0, z: 0 }}, backgroundColor: "transparent" }});
        view.zoomTo();
        view.render();
        return;
      }}

      // Fetch PDB - path relative to project root when served via http.server
      const base = window.location.pathname.includes("/visualization/") ? "../" : "./";
      const url = base + (row.pdb_path || "");
      fetch(url).then(res => res.text()).then(pdb => {{
        view.addModel(pdb, "pdb");
        view.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});

        let legendParts = [];
        if ((row.domains || []).length > 0) {{
          // Domain coloring from InterPro
          (row.domains || []).forEach(d => {{
            const sig = (d.signature || "").split(".")[0];
            const color = DATA.domain_colors[sig] || "#64748b";
            view.setStyle({{ resi: d.start + "-" + d.end }}, {{ cartoon: {{ color: color }} }});
          }});
          const legMap = {{ "PF05168": "HEPN (NUC)", "PF01228": "HEL (REC)", "PF18456": "WYL (REC)" }};
          const used = [...new Set((row.domains || []).map(d => (d.signature || "").split(".")[0]))];
          legendParts = used.map(s => `<span style="color:${{DATA.domain_colors[s] || '#64748b'}}">&#9679;</span> ${{legMap[s] || s}}`);
        }} else if ((row.nuc_spans || []).length > 0 || (row.rec_spans || []).length > 0) {{
          // NUC/REC lobe coloring (1-based resi)
          (row.nuc_spans || []).forEach(([s,e]) => view.setStyle({{ resi: (s+1) + "-" + e }}, {{ cartoon: {{ color: DATA.lobe_colors.nuc }} }}));
          (row.rec_spans || []).forEach(([s,e]) => view.setStyle({{ resi: (s+1) + "-" + e }}, {{ cartoon: {{ color: DATA.lobe_colors.rec }} }}));
          legendParts = [`<span style="color:${{DATA.lobe_colors.nuc}}">&#9679;</span> NUC (nuclease)`, `<span style="color:${{DATA.lobe_colors.rec}}">&#9679;</span> REC (recognition)`];
        }}

        view.zoomTo();
        view.render();

        const legend = document.getElementById("domain-legend");
        legend.innerHTML = legendParts.join(" ") || "";
      }}).catch(() => {{
        view.addLabel("Failed to load PDB", {{ position: {{ x: 0, y: 0, z: 0 }}, backgroundColor: "transparent" }});
        view.zoomTo();
        view.render();
      }});
    }}

    renderTable();
    if (DATA.rows.length) selectRow(0);
  </script>
</body>
</html>
'''

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[SUCCESS] Wrote {out_path}")
    print(f"   [*] Serve via: python -m http.server 8000  (from project root)")
    print(f"   [*] Open: http://localhost:8000/visualization/structure_dashboard.html")


if __name__ == "__main__":
    main()
