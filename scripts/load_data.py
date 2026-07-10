"""Day 1 sanity check: load the 10x 5' GEX + VDJ demo dataset into AnnData.

Dataset: 10k Human PBMCs, 5' v2 (healthy donor), from 10x Genomics.
  - GEX:  filtered feature-barcode matrix (.h5)
  - VDJ:  filtered T-cell contig annotations (.csv)

This does NO analysis. It only confirms the environment is working and the
data loads cleanly, then prints basic shape info as a sanity check.
"""

from pathlib import Path

import mudata as md
import scanpy as sc
import scirpy as ir

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"
GEX_H5 = DATA_DIR / "sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5"
VDJ_CSV = DATA_DIR / "sc5p_v2_hs_PBMC_10k_t_filtered_contig_annotations.csv"


def main() -> None:
    print("scanpy", sc.__version__, "| scirpy", ir.__version__)

    # --- GEX ---
    adata = sc.read_10x_h5(GEX_H5)
    adata.var_names_make_unique()
    print(f"\nGEX AnnData loaded: n_obs (cells) = {adata.n_obs:,}, "
          f"n_vars (genes) = {adata.n_vars:,}")

    # --- VDJ (TCR) ---
    # scirpy reads the 10x contig CSV into an AIRR-formatted AnnData/MuData.
    adata_tcr = ir.io.read_10x_vdj(VDJ_CSV)
    print(f"TCR AnnData loaded: n_obs (cells with a contig) = "
          f"{adata_tcr.n_obs:,}")

    # --- Combine GEX + VDJ into a single MuData object ---
    # scirpy >= 0.13 pairs the two modalities via MuData (merge_with_ir was
    # removed). Barcodes are matched across the "gex" and "airr" modalities.
    mdata = md.MuData({"gex": adata, "airr": adata_tcr})
    n_paired = len(adata.obs_names.intersection(adata_tcr.obs_names))
    print(f"\nMuData built: {mdata['gex'].n_obs:,} GEX cells + "
          f"{mdata['airr'].n_obs:,} TCR cells")
    print(f"Cells with both GEX and a TCR: {n_paired:,} "
          f"({100 * n_paired / adata.n_obs:.1f}% of GEX cells)")

    print("\nOK: dataset loaded cleanly into AnnData / MuData objects.")


if __name__ == "__main__":
    main()
