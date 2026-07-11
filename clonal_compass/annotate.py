"""Stage 3: marker-based cell-type / T-cell subset annotation.

Transparent, score-based labelling: each signature is scored per cell with
`sc.tl.score_genes`, averaged per Leiden cluster, and each cluster is assigned
the top-scoring label. Clusters that score as "T cell" are refined to a T-cell
subset. This is intentionally simple and auditable rather than a black-box
classifier (see interpretation guardrails in CLAUDE.md).
"""

import scanpy as sc

from .markers import LINEAGE_MARKERS, TCELL_SUBSET_MARKERS


def _present(adata, genes):
    return [g for g in genes if g in adata.var_names]


def _score_signatures(adata, signatures, prefix, seed=0):
    """Add one score column per signature: f'{prefix}:{name}'."""
    for name, genes in signatures.items():
        present = _present(adata, genes)
        if not present:
            print(f"[annotate] WARNING: no marker genes found for {name}")
            continue
        sc.tl.score_genes(
            adata, present, score_name=f"{prefix}:{name}",
            use_raw=False, random_state=seed,
        )


def annotate_clusters(
    adata,
    cluster_key: str = "leiden",
    lineage_markers: dict | None = None,
    tcell_subset_markers: dict | None = None,
):
    """Assign a lineage and a refined cell_type label to each cluster.

    Writes `.obs['lineage']` (broad) and `.obs['cell_type']` (T subsets where
    applicable). Returns the same AnnData.

    The signature dicts are parameters so a caller can swap in a dataset-specific
    marker set (see `markers.MARKER_SETS`); they default to the PBMC set. The
    scoring/labelling logic itself is unchanged and marker-set-agnostic.
    """
    if lineage_markers is None:
        lineage_markers = LINEAGE_MARKERS
    if tcell_subset_markers is None:
        tcell_subset_markers = TCELL_SUBSET_MARKERS

    _score_signatures(adata, lineage_markers, "lineage")
    _score_signatures(adata, tcell_subset_markers, "tsub")

    obs = adata.obs
    lineage_cols = [c for c in obs.columns if c.startswith("lineage:")]
    tsub_cols = [c for c in obs.columns if c.startswith("tsub:")]

    means = obs.groupby(cluster_key, observed=True)[lineage_cols + tsub_cols].mean()
    lineage_call = (
        means[lineage_cols].idxmax(axis=1).str.replace("lineage:", "", regex=False)
    )
    tsub_call = (
        means[tsub_cols].idxmax(axis=1).str.replace("tsub:", "", regex=False)
    )

    final = {
        cl: (tsub_call[cl] if lineage_call[cl] == "T cell" else lineage_call[cl])
        for cl in means.index
    }

    adata.obs["lineage"] = obs[cluster_key].map(lineage_call).astype("category")
    adata.obs["cell_type"] = obs[cluster_key].map(final).astype("category")

    summary = (
        adata.obs.groupby("cell_type", observed=True).size().sort_values(ascending=False)
    )
    print("[annotate] cells per cell_type:")
    for label, n in summary.items():
        print(f"    {label:<24} {n:>6,}")
    return adata
