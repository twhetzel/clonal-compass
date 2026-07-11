"""Data loading (stage 0).

Thin, re-runnable loaders that return a paired ``(gex, tcr)`` of plain AnnData
objects; pairing GEX<->TCR happens later by barcode. Two datasets are wired up
behind a small registry (see ``DATASETS`` / ``load_dataset``):

- ``pbmc``   — the 10x 5' GEX + VDJ demo dataset (healthy donor). Default.
- ``cancer`` — scirpy's built-in "3k T cells from cancer" (Wu et al. 2020,
  tumor-infiltrating T cells across multiple cancer types), fetched via
  ``ir.datasets.wu2020_3k()`` and cached by scirpy.

Both feed the *identical* downstream pipeline (QC -> cluster -> annotate ->
clonal); the registry only differs in how the raw objects are obtained and in
a filename ``suffix`` so the two datasets' artifacts coexist on disk.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

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


# --------------------------------------------------------------------------- #
# Dataset registry
# --------------------------------------------------------------------------- #
def _load_pbmc():
    """10x PBMC demo: separate GEX (.h5) and TCR (.csv) files on disk."""
    return load_gex(), load_tcr()


def _load_cancer():
    """Wu et al. 2020 tumor-infiltrating T cells (scirpy's ``wu2020_3k``).

    Returns the two modalities of scirpy's bundled MuData as standalone
    AnnData objects, shaped exactly like ``_load_pbmc`` returns them:
    ``gex`` carries raw integer counts; ``airr`` carries the per-cell TCR
    chains in ``.obsm['airr']`` ready for ``ir.pp.index_chains``. All 3,000
    cells are TCR-containing and already barcode-paired across modalities.
    """
    mdata = ir.datasets.wu2020_3k()
    gex = mdata.mod["gex"].copy()
    gex.var_names_make_unique()
    tcr = mdata.mod["airr"].copy()
    return gex, tcr


@dataclass(frozen=True)
class DatasetSpec:
    """One selectable dataset: how to load it and how to name its artifacts."""

    key: str
    display_name: str
    suffix: str  # appended to every output filename ("" keeps PBMC paths stable)
    loader: Callable  # -> (adata_gex, adata_tcr)
    marker_set: str  # key into markers.MARKER_SETS for annotation signatures


DATASETS: dict[str, DatasetSpec] = {
    "pbmc": DatasetSpec(
        key="pbmc",
        display_name="10k Human PBMC 5' (GEX + VDJ)",
        suffix="",
        loader=_load_pbmc,
        marker_set="pbmc",
    ),
    "cancer": DatasetSpec(
        key="cancer",
        display_name="3k tumor-infiltrating T cells (Wu et al. 2020)",
        suffix="_cancer",
        loader=_load_cancer,
        marker_set="cancer",
    ),
}


def load_dataset(key: str = "pbmc"):
    """Resolve ``key`` in the registry and load it.

    Returns ``(adata_gex, adata_tcr, spec)`` where ``spec`` is the matching
    :class:`DatasetSpec` (carries the display name + filename suffix the
    orchestration scripts use to keep both datasets' outputs side by side).
    """
    if key not in DATASETS:
        raise ValueError(
            f"unknown dataset {key!r}; choose one of {sorted(DATASETS)}"
        )
    spec = DATASETS[key]
    gex, tcr = spec.loader()
    return gex, tcr, spec
