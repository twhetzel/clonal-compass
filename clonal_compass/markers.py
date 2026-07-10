"""Standard marker-gene signatures for PBMC lineages and T-cell subsets.

These are deliberately small, well-established markers so the annotation is
transparent and auditable (see interpretation guardrails in CLAUDE.md). They
are used for score-based, per-cluster labelling in `annotate.py`.
"""

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
