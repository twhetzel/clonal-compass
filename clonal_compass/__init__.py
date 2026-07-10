"""Clonal Compass core analysis pipeline.

Each stage is a separate, re-runnable function operating on AnnData/scirpy
objects (see the individual modules). The intended order is:

    io -> qc -> cluster -> annotate -> clonal

`scripts/run_pipeline.py` chains them and writes figures + processed data.
"""

from . import annotate, clonal, cluster, io, markers, plots, qc

__all__ = ["io", "qc", "cluster", "annotate", "clonal", "markers", "plots"]
