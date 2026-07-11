"""Day 2 orchestrator: run the full core pipeline and save sanity-check plots.

Stages (each a re-runnable function in the `clonal_compass` package):
    0. load        GEX + TCR
    1. qc          filter low-quality cells / rare genes
    2. cluster     normalize -> PCA -> UMAP -> Leiden      [saves UMAP by cluster]
    3. annotate    marker-based cell-type labels           [saves UMAP by cell type]
    4. clonal      clonotypes + expansion, merged on cells [saves UMAP by clone size]

Intermediate objects are written to data/processed/ so any later stage can be
re-run without redoing earlier ones. Figures land in figures/.

Two datasets share this exact pipeline (select with --dataset):
    pbmc   (default)  10x PBMC 5' demo (healthy donor)
    cancer            Wu et al. 2020 tumor-infiltrating T cells (wu2020_3k)
Each dataset's artifacts get a filename suffix so both coexist on disk
(PBMC keeps the original unsuffixed names).

Run:  .venv/bin/python scripts/run_pipeline.py [--dataset pbmc|cancer]
"""

import argparse
from pathlib import Path

from clonal_compass._warnings import silence_demo_warnings

silence_demo_warnings()

import scanpy as sc  # noqa: E402

from clonal_compass import annotate, clonal, cluster, io, markers, plots, qc  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
FIG = ROOT / "figures"
PROC.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 1


def main(dataset: str = "pbmc") -> None:
    # --- Stage 0: load (dataset chosen from the io registry) ---
    adata, adata_tcr, spec = io.load_dataset(dataset)
    sfx = spec.suffix
    print(f"[pipeline] dataset={spec.key!r} ({spec.display_name}); suffix={sfx!r}")

    # --- Stage 1: QC ---
    adata = qc.run_qc(adata)

    # --- Stage 2: normalize + cluster ---
    adata = cluster.normalize_and_cluster(adata)
    adata.write(PROC / f"gex_clustered{sfx}.h5ad")
    plots.save_umap(
        adata, "leiden", FIG / f"umap_clusters{sfx}.png",
        title="Leiden clusters (transcriptome)",
    )

    # --- Stage 3: annotate (marker set chosen per dataset) ---
    mset = markers.MARKER_SETS[spec.marker_set]
    print(f"[annotate] using {spec.marker_set!r} marker set")
    adata = annotate.annotate_clusters(
        adata,
        lineage_markers=mset.lineage,
        tcell_subset_markers=mset.tcell_subset,
    )
    plots.save_umap(
        adata, "cell_type", FIG / f"umap_celltype{sfx}.png",
        title="Cell type (marker-based)",
    )

    # --- Stage 4: clonal expansion ---
    adata_tcr = clonal.compute_clonal_expansion(adata_tcr)
    adata = clonal.merge_clone_size(adata, adata_tcr)
    plots.save_umap(
        adata, "log_clone_size", FIG / f"umap_clonesize{sfx}.png",
        title="TCR clone size (log1p of cells per clonotype)",
        color_map="viridis",
    )
    plots.save_umap(
        adata, "clonal_expansion", FIG / f"umap_expansion{sfx}.png",
        title="Clonal expansion category",
    )

    # --- Persist final objects ---
    adata.write(PROC / f"gex_annotated{sfx}.h5ad")
    adata_tcr.write(PROC / f"tcr_clonotypes{sfx}.h5ad")
    print("\n[done] processed objects in data/processed/, figures in figures/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset", default="pbmc", choices=sorted(io.DATASETS),
        help="which dataset to run through the pipeline (default: pbmc)",
    )
    args = parser.parse_args()
    main(args.dataset)
