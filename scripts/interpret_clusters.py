"""Day 3 console tool: print per-cluster evidence + interpretation.

Inspection front-end over the shared report layer. For the judge-facing
artifact (HTML), use scripts/generate_report.py instead.

Run:  .venv/bin/python scripts/interpret_clusters.py
(Set ANTHROPIC_API_KEY to use Claude; otherwise a deterministic fallback runs.)
"""

from pathlib import Path

from clonal_compass._warnings import silence_demo_warnings

silence_demo_warnings()

import scanpy as sc  # noqa: E402

from clonal_compass import report  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"


def main() -> None:
    adata = sc.read_h5ad(PROC / "gex_annotated.h5ad")
    adata_tcr = sc.read_h5ad(PROC / "tcr_clonotypes.h5ad")

    data = report.build_report_data(adata, adata_tcr)

    for cr in data.clusters:
        ev = cr.evidence
        header = f"Cluster {ev.cluster_id}: {ev.cell_type} ({ev.n_cells} cells)"
        print(f"\n{'=' * 70}\n{header}  [{cr.interpretation.source}]\n{'=' * 70}")
        print(ev.as_prompt())
        print(cr.interpretation.text)

    if data.epitope_overview:
        print(f"\n{'=' * 70}\nDataset-wide predicted epitope specificities (VDJdb)")
        print("=" * 70)
        for hit in data.epitope_overview:
            print(f"  - {hit}")


if __name__ == "__main__":
    main()
