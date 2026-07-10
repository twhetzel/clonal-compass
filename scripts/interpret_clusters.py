"""Day 3 orchestrator: interpret notable clusters + cross-reference VDJdb.

Loads the processed objects from Day 2, computes per-cluster marker genes,
annotates TCRs with predicted VDJdb epitopes, then writes a hedged
plain-language interpretation for the most clonally-expanded clusters.

Run:  .venv/bin/python scripts/interpret_clusters.py
(Set ANTHROPIC_API_KEY to use Claude; otherwise a deterministic fallback runs.)
"""

from pathlib import Path

import scanpy as sc

from clonal_compass import epitope, interpret

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
REPORTS = ROOT / "reports"
REPORTS.mkdir(parents=True, exist_ok=True)

N_MARKERS = 8
N_CLUSTERS_TO_REPORT = 3
MIN_TCR_CELLS = 20  # skip clusters with too few TCR+ cells to say anything


def _top_markers(adata, group):
    """Top marker genes + scores for one leiden group from rank_genes_groups."""
    rg = adata.uns["rank_genes_groups"]
    names = rg["names"][group][:N_MARKERS]
    scores = rg["scores"][group][:N_MARKERS]
    return list(zip(names.tolist(), (float(s) for s in scores)))


def build_evidence(adata, adata_tcr, group, epitope_cols):
    sub = adata.obs[adata.obs["leiden"] == group]
    barcodes = sub.index
    n_cells = len(sub)
    n_with_tcr = int((sub["clone_size"] > 0).sum())
    n_expanded = int((sub["clone_size"] > 1).sum())
    pct_with_tcr = 100 * n_with_tcr / n_cells if n_cells else 0.0
    pct_expanded = 100 * n_expanded / n_with_tcr if n_with_tcr else 0.0
    top_clone_sizes = (
        sorted(sub.loc[sub["clone_size"] > 1, "clone_size"].astype(int), reverse=True)[:5]
    )
    cell_type_mode = sub["cell_type"].mode()
    cell_type = str(cell_type_mode.iat[0]) if not cell_type_mode.empty else "Unknown"
    hits = epitope.hits_for_cells(adata_tcr, barcodes, epitope_cols)

    return interpret.ClusterEvidence(
        cluster_id=group,
        cell_type=str(cell_type),
        n_cells=n_cells,
        top_markers=_top_markers(adata, group),
        n_with_tcr=n_with_tcr,
        pct_with_tcr=pct_with_tcr,
        n_expanded_cells=n_expanded,
        pct_expanded=pct_expanded,
        top_clone_sizes=top_clone_sizes,
        epitope_hits=hits,
    )


def main() -> None:
    adata = sc.read_h5ad(PROC / "gex_annotated.h5ad")
    adata_tcr = sc.read_h5ad(PROC / "tcr_clonotypes.h5ad")

    # Per-cluster marker genes (log-normalized data lives in adata.raw).
    print("[interpret] ranking marker genes per cluster...")
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)

    # Predicted epitope annotation for all TCRs.
    adata_tcr, epitope_cols = epitope.annotate_epitopes(adata_tcr)

    # Rank clusters by number of clonally-expanded cells; report the top few.
    # (Number, not percent — a small cluster with 3% expansion isn't notable.)
    obs = adata.obs
    stats = []
    for group in obs["leiden"].cat.categories:
        sub = obs[obs["leiden"] == group]
        n_tcr = int((sub["clone_size"] > 0).sum())
        if n_tcr < MIN_TCR_CELLS:
            continue
        n_exp = int((sub["clone_size"] > 1).sum())
        stats.append((group, n_exp, n_tcr))
    notable = [g for g, _, _ in sorted(stats, key=lambda t: -t[1])][:N_CLUSTERS_TO_REPORT]
    print(f"[interpret] most-expanded clusters (by # expanded cells): {notable}")

    report_lines = ["# Clonal Compass — cluster interpretations\n"]
    for group in notable:
        ev = build_evidence(adata, adata_tcr, group, epitope_cols)
        text = interpret.interpret_cluster(ev)
        header = f"## Cluster {group}: {ev.cell_type} ({ev.n_cells} cells)"
        print(f"\n{'=' * 70}\n{header}\n{'=' * 70}")
        print(ev.as_prompt())
        print(text)
        report_lines += [header, "", "```", ev.as_prompt().rstrip(), "```", "", text, ""]

    out = REPORTS / "cluster_interpretations.md"
    out.write_text("\n".join(report_lines))
    print(f"\n[interpret] wrote {out}")


if __name__ == "__main__":
    main()
