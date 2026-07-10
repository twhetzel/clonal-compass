"""Data loading (stage 0).

Thin, re-runnable loaders for the 10x 5' GEX + VDJ demo dataset. Returns
plain AnnData objects; pairing GEX<->TCR happens later by barcode.
"""

from pathlib import Path

import scanpy as sc
import scirpy as ir

DATA_RAW = Path(__file__).resolve().parent.parent / "data" / "raw"
GEX_H5 = DATA_RAW / "sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5"
VDJ_CSV = DATA_RAW / "sc5p_v2_hs_PBMC_10k_t_filtered_contig_annotations.csv"


def load_gex(path: Path = GEX_H5):
    """Load the filtered gene-expression matrix into AnnData."""
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    return adata


def load_tcr(path: Path = VDJ_CSV):
    """Load the 10x T-cell VDJ contig annotations into an AIRR AnnData."""
    return ir.io.read_10x_vdj(path)
