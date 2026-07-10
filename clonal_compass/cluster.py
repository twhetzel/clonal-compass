"""Stage 2: normalization, dimensionality reduction, and clustering.

Standard Scanpy workflow. Raw counts are preserved in `.layers['counts']`
and the log-normalized matrix in `.raw` so downstream marker scoring and
plotting use interpretable values (not the scaled matrix used for PCA).
"""

import scanpy as sc


def normalize_and_cluster(
    adata,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    resolution: float = 1.0,
    seed: int = 0,
):
    """Normalize -> HVG -> scale -> PCA -> neighbors -> UMAP -> Leiden.

    Returns the same AnnData with `.obsm['X_pca']`, `.obsm['X_umap']`, and a
    `leiden` cluster label in `.obs`. Re-runnable on a QC'd AnnData.
    """
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # log-normalized snapshot for scoring/plotting

    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)

    # Scale + PCA on HVGs only, in a copy, so adata.X stays log-normalized.
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=n_pcs, random_state=seed)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    # flavor="igraph" avoids the optional leidenalg dependency.
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    n_clusters = adata.obs["leiden"].nunique()
    print(f"[cluster] {n_clusters} Leiden clusters at resolution={resolution}")
    return adata
