"""Plotting helpers.

Every saved figure gets a one-line caption/title stating the metric shown
(a project convention — see CLAUDE.md).
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless; save to file, never open a window
import matplotlib.pyplot as plt  # noqa: E402
import scanpy as sc  # noqa: E402


def save_umap(adata, color, path, title=None, **kwargs):
    """Save a UMAP colored by `color` to `path` (PNG)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ax = sc.pl.umap(
        adata, color=color, title=title or color, show=False, **kwargs
    )
    fig = ax.figure if hasattr(ax, "figure") else plt.gcf()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {path}")
    return path
