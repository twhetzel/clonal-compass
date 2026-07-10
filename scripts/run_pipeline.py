"""Day 2 orchestrator: run the full core pipeline and save sanity-check plots.

Stages (each a re-runnable function in the `clonal_compass` package):
    0. load        GEX + TCR
    1. qc          filter low-quality cells / rare genes
    2. cluster     normalize -> PCA -> UMAP -> Leiden      [saves UMAP by cluster]
    3. annotate    marker-based cell-type labels           [saves UMAP by cell type]
    4. clonal      clonotypes + expansion, merged on cells [saves UMAP by clone size]

Intermediate objects are written to data/processed/ so any later stage can be
re-run without redoing earlier ones. Figures land in figures/.

Run:  .venv/bin/python scripts/run_pipeline.py
"""

from pathlib import Path

import scanpy as sc

from clonal_compass import annotate, clonal, cluster, io, plots, qc

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
FIG = ROOT / "figures"
PROC.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 1


def main() -> None:
    # --- Stage 0: load ---
    adata = io.load_gex()
    adata_tcr = io.load_tcr()

    # --- Stage 1: QC ---
    adata = qc.run_qc(adata)

    # --- Stage 2: normalize + cluster ---
    adata = cluster.normalize_and_cluster(adata)
    adata.write(PROC / "gex_clustered.h5ad")
    plots.save_umap(
        adata, "leiden", FIG / "umap_clusters.png",
        title="Leiden clusters (transcriptome)",
    )

    # --- Stage 3: annotate ---
    adata = annotate.annotate_clusters(adata)
    plots.save_umap(
        adata, "cell_type", FIG / "umap_celltype.png",
        title="Cell type (marker-based)",
    )

    # --- Stage 4: clonal expansion ---
    adata_tcr = clonal.compute_clonal_expansion(adata_tcr)
    adata = clonal.merge_clone_size(adata, adata_tcr)
    plots.save_umap(
        adata, "log_clone_size", FIG / "umap_clonesize.png",
        title="TCR clone size (log1p of cells per clonotype)",
        color_map="viridis",
    )
    plots.save_umap(
        adata, "clonal_expansion", FIG / "umap_expansion.png",
        title="Clonal expansion category",
    )

    # --- Persist final objects ---
    adata.write(PROC / "gex_annotated.h5ad")
    adata_tcr.write(PROC / "tcr_clonotypes.h5ad")
    print("\n[done] processed objects in data/processed/, figures in figures/")


if __name__ == "__main__":
    main()
