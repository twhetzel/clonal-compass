"""Silence known-benign warnings so the demo console is clean.

These are all deferred-from-Day-2 noise, not real problems:
  - mudata `pull_on_update` FutureWarning (scirpy wraps AnnData in MuData
    internally) — adopt the new behaviour explicitly.
  - anndata awkward-array ExperimentalFeatureWarning (VDJ data lives in an
    awkward array in .obsm["airr"]).
  - scanpy DataFrame-fragmentation PerformanceWarning from rank_genes_groups.

Call `silence_demo_warnings()` once at the top of each entry-point script.
"""

from __future__ import annotations

import warnings


def silence_demo_warnings() -> None:
    # Adopt the future mudata behaviour (and silence its FutureWarning).
    try:
        import mudata

        mudata.set_options(pull_on_update=False)
    except Exception:  # noqa: BLE001 - best effort, never fail the demo
        pass
    warnings.filterwarnings("ignore", message=".*pull_on_update.*")

    # awkward-array experimental support in anndata.
    try:
        from anndata._warnings import ExperimentalFeatureWarning

        warnings.filterwarnings("ignore", category=ExperimentalFeatureWarning)
    except Exception:  # noqa: BLE001
        warnings.filterwarnings("ignore", message=".*Awkward.*")

    # scanpy rank_genes_groups fragments its results DataFrame (harmless).
    try:
        from pandas.errors import PerformanceWarning

        warnings.filterwarnings("ignore", category=PerformanceWarning)
    except Exception:  # noqa: BLE001
        pass

    # scirpy ir_query uses a pandas idiom that emits a deprecation FutureWarning.
    warnings.filterwarnings(
        "ignore", message=".*Series.__getitem__ treating keys as positions.*"
    )
