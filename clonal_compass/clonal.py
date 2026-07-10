"""Stage 4: clonal expansion metrics from the VDJ data, merged onto clusters.

Runs scirpy's clonotype-definition workflow on the TCR AnnData, then merges
per-cell clone size and an expansion category onto the GEX AnnData by barcode.
Cells with no TCR get clone_size 0 / NaN clone_id.
"""

import numpy as np
import scirpy as ir


def compute_clonal_expansion(adata_tcr, clonotype_key: str = "cc_aa_identity"):
    """Index chains -> QC -> distance -> clonotype clusters -> expansion.

    Groups cells into clonotypes by **CDR3 amino-acid identity** (via
    `ir.tl.define_clonotype_clusters`, default `sequence='aa'`). This is more
    permissive than exact nucleotide matching: synonymous-codon variants of the
    same CDR3 protein collapse into one clonotype. The resulting cluster column
    is `cc_aa_identity`, aliased to `clone_id` for a uniform downstream merge.

    Adds `clone_id`, `clonal_expansion`, and `clone_size` to
    `adata_tcr.obs`. Re-runnable on a freshly loaded TCR AnnData.
    """
    ir.pp.index_chains(adata_tcr)
    ir.tl.chain_qc(adata_tcr)
    # Amino-acid CDR3 identity: ir_dist must be computed on the aa sequence.
    ir.pp.ir_dist(adata_tcr, metric="identity", sequence="aa")
    ir.tl.define_clonotype_clusters(
        adata_tcr, sequence="aa", metric="identity",
        receptor_arms="all", dual_ir="primary_only",
    )
    adata_tcr.obs["clone_id"] = adata_tcr.obs[clonotype_key]
    ir.tl.clonal_expansion(adata_tcr, target_col="clone_id")

    counts = adata_tcr.obs["clone_id"].value_counts()
    adata_tcr.obs["clone_size"] = (
        adata_tcr.obs["clone_id"].map(counts).astype("float")
    )

    n_clones = adata_tcr.obs["clone_id"].nunique()
    n_expanded = int((adata_tcr.obs["clone_size"] > 1).sum())
    print(
        f"[clonal] {n_clones:,} distinct clonotypes; "
        f"{n_expanded:,} cells in expanded clones (size > 1)"
    )
    return adata_tcr


def merge_clone_size(adata_gex, adata_tcr):
    """Join clone_id / clone_size / clonal_expansion onto GEX obs by barcode."""
    cols = ["clone_id", "clone_size", "clonal_expansion"]
    # Drop any prior merge so re-running the stage is idempotent.
    adata_gex.obs = adata_gex.obs.drop(
        columns=[c for c in cols + ["log_clone_size"] if c in adata_gex.obs.columns]
    )
    adata_gex.obs = adata_gex.obs.join(adata_tcr.obs[cols], how="left")
    adata_gex.obs["clone_size"] = adata_gex.obs["clone_size"].fillna(0.0)
    # log1p makes the mostly-size-1 distribution legible on a UMAP color scale.
    adata_gex.obs["log_clone_size"] = np.log1p(adata_gex.obs["clone_size"])

    n_paired = int((adata_gex.obs["clone_size"] > 0).sum())
    print(
        f"[clonal] merged onto GEX: {n_paired:,}/{adata_gex.n_obs:,} cells "
        f"have a TCR ({100 * n_paired / adata_gex.n_obs:.1f}%)"
    )
    return adata_gex
