"""Day 4: report generation layer.

Combines the Day-2 UMAP figures, the per-cluster interpretations (interpret.py),
and the VDJdb epitope matches (epitope.py) into a clean, self-contained report.

Owns the shared "build the story" logic (used by both the console tool and the
HTML generator) so the heavy compute lives in one place:
  - rank marker genes per cluster
  - annotate predicted epitopes
  - pick the most clonally-expanded clusters
  - build evidence + interpretation for each

Every interpretation carries its source ("Claude API" / "deterministic
fallback") so the report can tag it unambiguously.
"""

from __future__ import annotations

import base64
import html
from dataclasses import dataclass
from pathlib import Path

import scanpy as sc

from . import epitope, interpret

N_MARKERS = 8
N_CLUSTERS_TO_REPORT = 3
MIN_TCR_CELLS = 20  # skip clusters with too few TCR+ cells to say anything

# The four Day-2 UMAPs, in display order, with captions.
FIGURES = [
    ("umap_clusters.png", "Leiden clusters on the transcriptome (UMAP)."),
    ("umap_celltype.png", "Marker-based cell-type / T-cell subset labels."),
    ("umap_clonesize.png", "TCR clone size (log1p of cells per clonotype)."),
    ("umap_expansion.png", "Clonal-expansion category per cell."),
]


@dataclass
class ClusterReport:
    evidence: interpret.ClusterEvidence
    interpretation: interpret.Interpretation


@dataclass
class ReportData:
    clusters: list[ClusterReport]
    epitope_overview: list[str]   # dataset-wide predicted specificities
    n_tcr_matched: int            # TCR+ cells with any predicted match


# --------------------------------------------------------------------------- #
# Build the story
# --------------------------------------------------------------------------- #
def _top_markers(adata, group):
    rg = adata.uns["rank_genes_groups"]
    names = rg["names"][group][:N_MARKERS]
    scores = rg["scores"][group][:N_MARKERS]
    return list(zip(names.tolist(), (float(s) for s in scores)))


def _build_evidence(adata, adata_tcr, group, epitope_cols):
    sub = adata.obs[adata.obs["leiden"] == group]
    barcodes = sub.index
    n_cells = len(sub)
    n_with_tcr = int((sub["clone_size"] > 0).sum())
    n_expanded = int((sub["clone_size"] > 1).sum())
    pct_with_tcr = 100 * n_with_tcr / n_cells if n_cells else 0.0
    pct_expanded = 100 * n_expanded / n_with_tcr if n_with_tcr else 0.0
    top_clone_sizes = sorted(
        sub.loc[sub["clone_size"] > 1, "clone_size"].astype(int), reverse=True
    )[:5]
    mode = sub["cell_type"].mode()
    cell_type = str(mode.iat[0]) if not mode.empty else "Unknown"
    hits = epitope.hits_for_cells(adata_tcr, barcodes, epitope_cols)

    return interpret.ClusterEvidence(
        cluster_id=group,
        cell_type=cell_type,
        n_cells=n_cells,
        top_markers=_top_markers(adata, group),
        n_with_tcr=n_with_tcr,
        pct_with_tcr=pct_with_tcr,
        n_expanded_cells=n_expanded,
        pct_expanded=pct_expanded,
        top_clone_sizes=top_clone_sizes,
        epitope_hits=hits,
    )


def _select_notable(adata):
    """Return the leiden groups with the most clonally-expanded cells."""
    obs = adata.obs
    stats = []
    for group in obs["leiden"].cat.categories:
        sub = obs[obs["leiden"] == group]
        n_tcr = int((sub["clone_size"] > 0).sum())
        if n_tcr < MIN_TCR_CELLS:
            continue
        n_exp = int((sub["clone_size"] > 1).sum())
        stats.append((group, n_exp))
    return [g for g, _ in sorted(stats, key=lambda t: -t[1])][:N_CLUSTERS_TO_REPORT]


def build_report_data(adata, adata_tcr) -> ReportData:
    """Compute everything the report needs from the processed objects."""
    print("[report] ranking marker genes per cluster...")
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)

    adata_tcr, epitope_cols = epitope.annotate_epitopes(adata_tcr)

    notable = _select_notable(adata)
    print(f"[report] most-expanded clusters (by # expanded cells): {notable}")

    clusters = []
    for group in notable:
        ev = _build_evidence(adata, adata_tcr, group, epitope_cols)
        interp = interpret.interpret_cluster(ev)
        print(f"[report]   cluster {group} ({ev.cell_type}) -> {interp.source}")
        clusters.append(ClusterReport(ev, interp))

    # Dataset-wide predicted specificities (all TCR+ cells, not just expanded).
    overview = epitope.hits_for_cells(
        adata_tcr, adata_tcr.obs.index, epitope_cols, expanded_only=False
    )
    n_matched = 0
    if epitope_cols:
        n_matched = int(adata_tcr.obs[epitope_cols].notna().any(axis=1).sum())

    return ReportData(clusters=clusters, epitope_overview=overview, n_tcr_matched=n_matched)


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #
def _img_data_uri(path: Path) -> str:
    b64 = base64.b64encode(path.read_bytes()).decode()
    return f"data:image/png;base64,{b64}"


def _source_badge(source: str) -> str:
    kind = "claude" if source == "Claude API" else "fallback"
    return f'<span class="badge badge-{kind}">{html.escape(source)}</span>'


def _interp_html(text: str) -> str:
    paras = [p.strip() for p in text.split("\n") if p.strip()]
    return "".join(f"<p>{html.escape(p)}</p>" for p in paras) or "<p></p>"


def _evidence_rows(ev: interpret.ClusterEvidence) -> str:
    markers = ", ".join(f"{html.escape(g)} ({s:+.2f})" for g, s in ev.top_markers)
    clones = ", ".join(str(c) for c in ev.top_clone_sizes) or "none > 1"
    rows = [
        ("Provisional label", html.escape(ev.cell_type)),
        ("Cells in cluster", f"{ev.n_cells:,}"),
        ("Top marker genes (rank_genes_groups score)", markers),
        ("TCR coverage", f"{ev.n_with_tcr:,} / {ev.n_cells:,} ({ev.pct_with_tcr:.0f}%)"),
        (
            "Clonal expansion",
            f"{ev.n_expanded_cells:,} cells in clones &gt; 1 "
            f"({ev.pct_expanded:.0f}% of TCR+ cells)",
        ),
        ("Largest clone sizes here", clones),
    ]
    return "".join(
        f"<tr><th>{label}</th><td>{value}</td></tr>" for label, value in rows
    )


def _epitope_block(ev: interpret.ClusterEvidence) -> str:
    if ev.epitope_hits:
        items = "".join(f"<li>{html.escape(h)}</li>" for h in ev.epitope_hits)
        return (
            '<div class="epitope epitope-hit"><strong>Predicted epitope matches '
            "(VDJdb, exact CDR3 aa identity)</strong> &mdash; predicted specificity, "
            f"not confirmed reactivity:<ul>{items}</ul></div>"
        )
    return (
        '<div class="epitope epitope-none"><strong>No VDJdb match</strong> for this '
        "cluster&rsquo;s expanded clones. Absence of a match is a normal, valid "
        "outcome &mdash; not a failure.</div>"
    )


def _cluster_card(cr: ClusterReport) -> str:
    ev = cr.evidence
    return f"""
    <section class="card">
      <div class="card-head">
        <h3>Cluster {html.escape(ev.cluster_id)}: {html.escape(ev.cell_type)}</h3>
        {_source_badge(cr.interpretation.source)}
      </div>
      <table class="evidence">{_evidence_rows(ev)}</table>
      <div class="interp">
        <div class="interp-label">Interpretation
          <span class="src-note">generated by: {html.escape(cr.interpretation.source)}</span>
        </div>
        {_interp_html(cr.interpretation.text)}
      </div>
      {_epitope_block(ev)}
    </section>"""


def render_html(data: ReportData, figures_dir: Path, meta: dict) -> str:
    figs = []
    for name, caption in FIGURES:
        path = figures_dir / name
        if not path.exists():
            continue
        figs.append(
            f'<figure><img src="{_img_data_uri(path)}" alt="{html.escape(caption)}">'
            f"<figcaption>{html.escape(caption)}</figcaption></figure>"
        )
    figures_html = "\n".join(figs)

    if data.epitope_overview:
        items = "".join(f"<li>{html.escape(h)}</li>" for h in data.epitope_overview)
        overview_html = (
            f"<p>{data.n_tcr_matched} TCR+ cell(s) carry a CDR3 that exactly matches "
            "a VDJdb entry. These are <em>predicted</em> specificities (database "
            "hits), not confirmed reactivity:</p>"
            f"<ul class='overview'>{items}</ul>"
        )
    else:
        overview_html = (
            "<p>No exact VDJdb CDR3 matches were found dataset-wide. This is a "
            "normal, valid outcome for a healthy-donor repertoire.</p>"
        )

    cards = "\n".join(_cluster_card(cr) for cr in data.clusters)
    meta_rows = "".join(
        f"<div class='stat'><span class='stat-num'>{v}</span>"
        f"<span class='stat-lab'>{html.escape(k)}</span></div>"
        for k, v in meta.items()
    )

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Clonal Compass &mdash; Report</title>
<style>{_CSS}</style></head><body>
<div class="wrap">
  <header>
    <h1>Clonal Compass</h1>
    <p class="sub">Single-cell immune-repertoire co-pilot &mdash; paired scRNA-seq + TCR-seq.
    This is an <strong>exploratory research</strong> analysis, not a diagnostic one.</p>
    <div class="stats">{meta_rows}</div>
  </header>

  <h2>1. Overview figures</h2>
  <div class="figgrid">{figures_html}</div>

  <h2>2. Interpretation of the most clonally-expanded clusters</h2>
  <p class="note">Each interpretation is tagged with the path that produced it.
  Claims are hedged and cite the specific metric behind them.</p>
  {cards}

  <h2>3. Predicted epitope specificities (VDJdb)</h2>
  <div class="card">{overview_html}
    <p class="fine">Matches are exact CDR3 amino-acid hits against VDJdb and indicate
    <em>possible/predicted</em> antigen specificity only; they require functional
    validation. A &ldquo;no match&rdquo; result is expected and valid.</p>
  </div>

  <footer>Generated by the Clonal Compass pipeline. Research/exploratory use only &mdash;
  no patient-level or clinical claims.</footer>
</div></body></html>"""


def render_markdown(data: ReportData) -> str:
    out = ["# Clonal Compass — cluster interpretations\n"]
    for cr in data.clusters:
        ev = cr.evidence
        out.append(f"## Cluster {ev.cluster_id}: {ev.cell_type} ({ev.n_cells} cells)")
        out.append(f"_Interpretation source: {cr.interpretation.source}_\n")
        out.append("```")
        out.append(ev.as_prompt().rstrip())
        out.append("```\n")
        out.append(cr.interpretation.text + "\n")
    return "\n".join(out)


_CSS = """
:root{--fg:#1a1f2b;--muted:#5b6472;--line:#e3e7ee;--bg:#f7f8fb;--card:#fff;
--accent:#2b6cb0;--claude:#1f7a4d;--claude-bg:#e6f4ec;--fallback:#9a6a00;--fallback-bg:#fdf2dc}
*{box-sizing:border-box}
body{margin:0;font:16px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
color:var(--fg);background:var(--bg)}
.wrap{max-width:940px;margin:0 auto;padding:40px 24px 64px}
header h1{margin:0 0 4px;font-size:34px;letter-spacing:-.02em}
.sub{color:var(--muted);margin:0 0 20px;max-width:70ch}
.stats{display:flex;gap:14px;flex-wrap:wrap;margin-bottom:8px}
.stat{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 16px;min-width:120px}
.stat-num{display:block;font-size:22px;font-weight:700}
.stat-lab{display:block;font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}
h2{margin:40px 0 12px;font-size:22px;border-bottom:2px solid var(--line);padding-bottom:6px}
.note,.fine{color:var(--muted);font-size:14px}
.figgrid{display:grid;grid-template-columns:1fr 1fr;gap:18px}
.figgrid figure{margin:0;background:var(--card);border:1px solid var(--line);border-radius:12px;padding:12px}
.figgrid img{width:100%;height:auto;border-radius:6px;display:block}
.figgrid figcaption{font-size:13px;color:var(--muted);margin-top:8px}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:20px 22px;margin:16px 0}
.card-head{display:flex;align-items:center;justify-content:space-between;gap:12px;flex-wrap:wrap}
.card-head h3{margin:0;font-size:19px}
.badge{font-size:12px;font-weight:700;padding:4px 10px;border-radius:999px;white-space:nowrap}
.badge-claude{color:var(--claude);background:var(--claude-bg)}
.badge-fallback{color:var(--fallback);background:var(--fallback-bg)}
table.evidence{width:100%;border-collapse:collapse;margin:14px 0;font-size:14px}
table.evidence th{text-align:left;color:var(--muted);font-weight:600;width:230px;
vertical-align:top;padding:6px 12px 6px 0;border-bottom:1px solid var(--line)}
table.evidence td{padding:6px 0;border-bottom:1px solid var(--line)}
.interp{margin:14px 0}
.interp-label{font-weight:700;font-size:13px;text-transform:uppercase;letter-spacing:.04em;color:var(--accent)}
.src-note{font-weight:600;color:var(--muted);text-transform:none;letter-spacing:0;margin-left:8px}
.interp p{margin:8px 0}
.epitope{border-radius:10px;padding:12px 14px;font-size:14px;margin-top:12px}
.epitope-hit{background:#eef4fb;border:1px solid #cfe0f4}
.epitope-none{background:#f2f4f7;border:1px solid var(--line);color:var(--muted)}
.epitope ul,.overview{margin:8px 0 0;padding-left:20px}
ul.overview li{margin:2px 0}
footer{margin-top:40px;padding-top:16px;border-top:1px solid var(--line);color:var(--muted);font-size:13px}
@media(max-width:680px){.figgrid{grid-template-columns:1fr}table.evidence th{width:auto}}
"""
