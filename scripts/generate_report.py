"""Day 4: generate the judge-facing report (self-contained HTML + markdown).

Combines the Day-2 UMAPs, per-cluster interpretations (tagged Claude API vs
deterministic fallback), and VDJdb epitope results into one clean HTML file
with the figures embedded (no external assets to carry around).

Two datasets share this exact reporting layer (select with --dataset); each
reads its own suffixed processed objects and writes its own suffixed artifacts
(PBMC keeps the original unsuffixed names, e.g. reports/cluster_evidence.json;
cancer writes reports/cluster_evidence_cancer.json).

Run:  .venv/bin/python scripts/generate_report.py [--dataset pbmc|cancer]
(Set ANTHROPIC_API_KEY first to tag interpretations as Claude API.)
"""

import argparse
from pathlib import Path

from clonal_compass._warnings import silence_demo_warnings

silence_demo_warnings()

import scanpy as sc  # noqa: E402  (after warning config)

from clonal_compass import io, report  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"
REPORTS = ROOT / "reports"
REPORTS.mkdir(parents=True, exist_ok=True)


def main(dataset: str = "pbmc") -> None:
    if dataset not in io.DATASETS:
        raise ValueError(
            f"unknown dataset {dataset!r}; choose one of {sorted(io.DATASETS)}"
        )
    spec = io.DATASETS[dataset]
    sfx = spec.suffix

    adata = sc.read_h5ad(PROC / f"gex_annotated{sfx}.h5ad")
    adata_tcr = sc.read_h5ad(PROC / f"tcr_clonotypes{sfx}.h5ad")

    data = report.build_report_data(adata, adata_tcr)

    meta = {
        "Cells": f"{adata.n_obs:,}",
        "Genes": f"{adata.n_vars:,}",
        "Clusters": adata.obs["leiden"].nunique(),
        "Cell types": adata.obs["cell_type"].nunique(),
        "TCR+ cells": f"{int((adata.obs['clone_size'] > 0).sum()):,}",
    }

    html_path = REPORTS / f"clonal_compass_report{sfx}.html"
    html_path.write_text(report.render_html(data, FIGURES, meta, suffix=sfx))

    md_path = REPORTS / f"cluster_interpretations{sfx}.md"
    md_path.write_text(report.render_markdown(data))

    # Compact per-cluster evidence for the chat interface (never reads .h5ad).
    report.write_evidence_json(
        adata, adata_tcr, data,
        REPORTS / f"cluster_evidence{sfx}.json",
        dataset_name=spec.display_name,
    )

    sources = {cr.interpretation.source for cr in data.clusters}
    print(f"\n[report] dataset={spec.key!r} ({spec.display_name})")
    print(f"[report] interpretation source(s): {', '.join(sorted(sources))}")
    print(f"[report] wrote {html_path}")
    print(f"[report] wrote {md_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset", default="pbmc", choices=sorted(io.DATASETS),
        help="which dataset's processed objects to report on (default: pbmc)",
    )
    args = parser.parse_args()
    main(args.dataset)
