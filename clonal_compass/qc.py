"""Stage 1: standard QC filtering.

Computes per-cell QC metrics and applies standard thresholds. Thresholds are
parameters so they can be tuned per dataset without editing code.
"""

import scanpy as sc


def run_qc(
    adata,
    min_genes: int = 200,
    min_cells: int = 3,
    max_pct_mt: float = 15.0,
    max_genes: int = 6000,
):
    """Filter low-quality cells and rarely-detected genes.

    - min_genes:  drop cells expressing fewer than this many genes
    - min_cells:  drop genes detected in fewer than this many cells
    - max_pct_mt: drop cells above this % mitochondrial counts (dying cells)
    - max_genes:  drop cells above this many genes (likely doublets)

    Returns a new, filtered AnnData (input is not modified in place).
    """
    adata = adata.copy()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    n_start = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy()
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes].copy()

    print(
        f"[qc] {n_start:,} -> {adata.n_obs:,} cells kept "
        f"({adata.n_vars:,} genes); "
        f"filters: min_genes={min_genes}, max_pct_mt={max_pct_mt}, "
        f"max_genes={max_genes}"
    )
    return adata
