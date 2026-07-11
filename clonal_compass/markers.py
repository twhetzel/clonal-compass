"""Marker-gene signatures for score-based, per-cluster labelling in `annotate.py`.

These are deliberately small, well-established markers so the annotation is
transparent and auditable (see interpretation guardrails in CLAUDE.md).

Two signature sets are provided behind a small registry (`MARKER_SETS`),
mirroring the dataset registry in `io.py`:

- ``pbmc``   — broad PBMC lineages + resting/circulating T-cell subsets. Default.
- ``cancer`` — tuned for **tumor-infiltrating T cells** (all cells are T cells):
  a single-entry ``T cell`` lineage layer (see the note on CANCER_LINEAGE_MARKERS
  for why a multi-lineage argmax fails on a pure-T population), plus a TIL-state
  subset panel covering exhaustion, cytotoxicity, Treg, Tfh, tissue-resident
  memory, and proliferation — states the PBMC panel doesn't resolve.

The annotation *logic* is unchanged and dataset-agnostic (`annotate.py`); only
which set it scores is selected per dataset via `DatasetSpec.marker_set`.
"""

from dataclasses import dataclass

# --------------------------------------------------------------------------- #
# PBMC (default) — healthy circulating repertoire
# --------------------------------------------------------------------------- #
# Broad PBMC lineages — used first to separate T cells from everything else.
LINEAGE_MARKERS = {
    "T cell": ["CD3D", "CD3E", "CD3G", "TRAC"],
    "NK": ["GNLY", "NKG7", "KLRD1", "NCAM1", "KLRF1"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "Monocyte": ["CD14", "LYZ", "FCGR3A", "S100A8", "S100A9"],
    "DC": ["FCER1A", "CST3", "CLEC10A", "CLEC4C"],
    "Platelet": ["PPBP", "PF4"],
}

# T-cell subsets — used to refine clusters that scored as "T cell".
TCELL_SUBSET_MARKERS = {
    "CD4 T naive": ["CD3D", "CD4", "CCR7", "SELL", "TCF7", "LEF1"],
    "CD4 T memory": ["CD3D", "CD4", "IL7R", "S100A4", "AQP3"],
    "CD8 T naive": ["CD3D", "CD8A", "CD8B", "CCR7", "SELL", "TCF7"],
    "CD8 T effector/memory": ["CD3D", "CD8A", "CD8B", "GZMK", "GZMB", "CCL5", "NKG7"],
    "Treg": ["CD3D", "CD4", "FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "MAIT": ["CD3D", "SLC4A10", "KLRB1", "TRAV1-2"],
    "gamma-delta T": ["CD3D", "TRDC", "TRGC1", "TRGC2"],
}


# --------------------------------------------------------------------------- #
# Cancer — tumor-infiltrating T cells (Wu et al. 2020)
# --------------------------------------------------------------------------- #
# Lineage layer: a single "T cell" entry, on purpose.
#
# This dataset is curated to TCR-containing cells — every cluster is T cells.
# `sc.tl.score_genes` measures a signature's enrichment *against a background of
# similarly-expressed genes*, so in a pure-T population the pan-T score has no
# non-T contrast to rise above and lands near zero (often negative). A
# multi-lineage argmax then gets decided by noise and mislabels real T cells as
# B/NK/myeloid. Asserting the one lineage the data actually is (T cell) is both
# honest and routes every cluster into the TIL-state panel below — which is
# where the real discrimination happens. (Genuine contaminant/doublet clusters,
# if any, are surfaced by the per-cluster rank_genes_groups in the report, not
# by this layer.) The shared annotate logic is unchanged: idxmax over a
# single-column score frame simply always returns "T cell".
CANCER_LINEAGE_MARKERS = {
    "T cell": ["CD3D", "CD3E", "CD3G", "CD2"],
}

# TIL differentiation states — the value-add over the PBMC subset panel. Each
# retains a T-lineage anchor (CD3D + CD4/CD8) so the score reflects a T-cell
# state, then adds the state-defining programme (exhaustion checkpoints,
# cytotoxic effectors, Treg, Tfh, tissue-resident memory, cell cycle).
CANCER_TCELL_SUBSET_MARKERS = {
    "CD4 T naive/central-memory": ["CD3D", "CD4", "CCR7", "SELL", "TCF7", "LEF1", "IL7R"],
    "CD4 Treg": ["CD3D", "CD4", "FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TNFRSF4"],
    "CD4 Tfh-like (CXCL13+)": ["CD3D", "CD4", "CXCL13", "PDCD1", "CXCR5", "ICOS", "BCL6"],
    "CD8 cytotoxic effector": ["CD3D", "CD8A", "CD8B", "GZMB", "GZMH", "PRF1", "GNLY", "NKG7", "FGFBP2"],
    "CD8 exhausted/terminal": ["CD3D", "CD8A", "PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "ENTPD1", "LAYN"],
    "CD8 tissue-resident memory": ["CD3D", "CD8A", "ITGAE", "ZNF683", "CXCR6", "ITGA1", "PRDM1"],
    # STMN1 deliberately excluded: broadly expressed in activated/effector T
    # cells, so it inflates the cycling score for non-dividing effectors. Keep
    # the cell-cycle-specific markers (G2/M: MKI67/TOP2A; S-phase: TYMS/PCLAF).
    "Proliferating (cycling)": ["MKI67", "TOP2A", "TYMS", "PCLAF"],
    "MAIT/innate-like": ["CD3D", "SLC4A10", "KLRB1", "TRAV1-2", "RORA"],
}


# --------------------------------------------------------------------------- #
# Registry — selected per dataset via DatasetSpec.marker_set
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MarkerSet:
    """A named pair of signature dicts for the two-stage annotation."""

    lineage: dict
    tcell_subset: dict


MARKER_SETS: dict[str, MarkerSet] = {
    "pbmc": MarkerSet(lineage=LINEAGE_MARKERS, tcell_subset=TCELL_SUBSET_MARKERS),
    "cancer": MarkerSet(
        lineage=CANCER_LINEAGE_MARKERS, tcell_subset=CANCER_TCELL_SUBSET_MARKERS
    ),
}
