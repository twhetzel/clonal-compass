"""Day 4: generate the judge-facing report (self-contained HTML + markdown).

Combines the Day-2 UMAPs, per-cluster interpretations (tagged Claude API vs
deterministic fallback), and VDJdb epitope results into one clean HTML file
with the figures embedded (no external assets to carry around).

Run:  .venv/bin/python scripts/generate_report.py
(Set ANTHROPIC_API_KEY first to tag interpretations as Claude API.)
"""

from pathlib import Path

from clonal_compass._warnings import silence_demo_warnings

silence_demo_warnings()

import scanpy as sc  # noqa: E402  (after warning config)

from clonal_compass import report  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"
REPORTS = ROOT / "reports"
REPORTS.mkdir(parents=True, exist_ok=True)


def main() -> None:
    adata = sc.read_h5ad(PROC / "gex_annotated.h5ad")
    adata_tcr = sc.read_h5ad(PROC / "tcr_clonotypes.h5ad")

    data = report.build_report_data(adata, adata_tcr)

    meta = {
        "Cells": f"{adata.n_obs:,}",
        "Genes": f"{adata.n_vars:,}",
        "Clusters": adata.obs["leiden"].nunique(),
        "Cell types": adata.obs["cell_type"].nunique(),
        "TCR+ cells": f"{int((adata.obs['clone_size'] > 0).sum()):,}",
    }

    html_path = REPORTS / "clonal_compass_report.html"
    html_path.write_text(report.render_html(data, FIGURES, meta))

    md_path = REPORTS / "cluster_interpretations.md"
    md_path.write_text(report.render_markdown(data))

    sources = {cr.interpretation.source for cr in data.clusters}
    print(f"\n[report] interpretation source(s): {', '.join(sorted(sources))}")
    print(f"[report] wrote {html_path}")
    print(f"[report] wrote {md_path}")


if __name__ == "__main__":
    main()
