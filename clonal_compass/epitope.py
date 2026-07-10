"""Day 3: epitope cross-referencing against VDJdb.

For TCRs (focusing on expanded clones), look up exact CDR3 amino-acid matches
in the VDJdb reference database via scirpy. Matches are *predicted* specificity,
not confirmed reactivity; the absence of a match is a normal, valid outcome.
"""

from __future__ import annotations

import tempfile
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd
import scanpy as sc
import scirpy as ir
from scirpy.io import AirrCell, from_airr_cells

# VDJdb reference columns we surface. These carry the predicted antigen.
REF_COLS = ["antigen.species", "antigen.epitope"]

# Cache the processed reference under the (gitignored) data dir.
_VDJDB_CACHE = Path(__file__).resolve().parent.parent / "data" / "reference" / "vdjdb.h5ad"
_LATEST_URL = "https://raw.githubusercontent.com/antigenomics/vdjdb-db/master/latest-version.txt"
_META_FIELDS = [
    "species", "mhc.a", "mhc.b", "mhc.class", "antigen.epitope",
    "antigen.gene", "antigen.species", "reference.id", "vdjdb.score",
]


def load_vdjdb():
    """Load VDJdb as a scirpy AnnData, building + caching if needed.

    scirpy 0.17's built-in ``ir.datasets.vdjdb()`` reads ``vdjdb_full.txt`` from
    the archive root, but current VDJdb releases nest it under a dated
    subdirectory, so the built-in loader raises FileNotFoundError. This
    reimplements the same processing with a glob for the real file path.
    """
    if _VDJDB_CACHE.exists():
        vdjdb = sc.read_h5ad(_VDJDB_CACHE)
        vdjdb.uns["DB"] = {"name": "VDJDB"}  # scirpy keys query results on this
        return vdjdb

    with urllib.request.urlopen(_LATEST_URL) as fh:
        release_url = fh.read().decode().split()[0]

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        urllib.request.urlretrieve(release_url, tmp / "vdjdb.zip")
        with zipfile.ZipFile(tmp / "vdjdb.zip") as zf:
            zf.extractall(tmp)
        matches = list(tmp.rglob("vdjdb_full.txt"))
        if not matches:
            raise FileNotFoundError("vdjdb_full.txt not found in VDJdb release")
        df = pd.read_csv(matches[0], sep="\t", low_memory=False)

    cells = []
    # to_dict("records") avoids building a pandas Series per row — much faster
    # than iterrows() over the ~140k VDJdb entries.
    for idx, row in enumerate(df.to_dict("records")):
        cell = AirrCell(cell_id=str(idx))
        if not pd.isnull(row["cdr3.alpha"]):
            chain = AirrCell.empty_chain_dict()
            chain.update({"locus": "TRA", "junction_aa": row["cdr3.alpha"],
                          "v_call": row["v.alpha"], "j_call": row["j.alpha"],
                          "consensus_count": 0, "productive": True})
            cell.add_chain(chain)
        if not pd.isnull(row["cdr3.beta"]):
            chain = AirrCell.empty_chain_dict()
            chain.update({"locus": "TRB", "junction_aa": row["cdr3.beta"],
                          "v_call": row["v.beta"], "d_call": row["d.beta"],
                          "j_call": row["j.beta"], "consensus_count": 0,
                          "productive": True})
            cell.add_chain(chain)
        for f in _META_FIELDS:
            cell[f] = row.get(f)
        cells.append(cell)

    vdjdb = from_airr_cells(cells)
    ir.pp.index_chains(vdjdb)
    vdjdb.uns["DB"] = {"name": "VDJDB"}  # scirpy keys query results on this
    _VDJDB_CACHE.parent.mkdir(parents=True, exist_ok=True)
    vdjdb.write_h5ad(_VDJDB_CACHE)
    print(f"[epitope] built + cached VDJdb reference ({vdjdb.n_obs:,} entries)")
    return vdjdb


def annotate_epitopes(adata_tcr):
    """Annotate each TCR cell with predicted VDJdb epitope matches.

    Requires that clonotypes have already been defined (index_chains + chain_qc
    run in clonal.compute_clonal_expansion). Adds VDJdb reference columns to
    `adata_tcr.obs` where an exact CDR3 aa match exists. Returns the annotated
    AnnData plus the list of column names that were added.
    """
    print("[epitope] loading VDJdb reference...")
    vdjdb = load_vdjdb()

    # Distance + query against the reference on CDR3 amino-acid identity.
    ir.pp.ir_dist(adata_tcr, vdjdb, metric="identity", sequence="aa")
    ir.tl.ir_query(
        adata_tcr, vdjdb, metric="identity", sequence="aa",
        receptor_arms="all", dual_ir="any",
    )
    before = set(adata_tcr.obs.columns)
    # unique-only: annotate a cell only when its predicted epitope is
    # unambiguous. Cells matching multiple conflicting entries stay NaN
    # (reported as "no match") rather than the noisy literal "ambiguous".
    ir.tl.ir_query_annotate(
        adata_tcr, vdjdb, metric="identity", sequence="aa",
        include_ref_cols=REF_COLS, strategy="unique-only",
    )
    added = [c for c in adata_tcr.obs.columns if c not in before]

    n_hit = int(adata_tcr.obs[added].notna().any(axis=1).sum()) if added else 0
    print(
        f"[epitope] {n_hit:,}/{adata_tcr.n_obs:,} TCR+ cells have a predicted "
        f"VDJdb match (columns: {added})"
    )
    return adata_tcr, added


def hits_for_cells(adata_tcr, barcodes, epitope_cols, expanded_only=True):
    """Summarize predicted epitope hits for a set of cells.

    Returns a list of 'species: epitope (n cells)' strings. When
    expanded_only is True, restricts to cells in clones of size > 1 (the clones
    the report should flag). An empty list is a valid 'no match' result.
    """
    if not epitope_cols:
        return []
    obs = adata_tcr.obs.loc[adata_tcr.obs.index.intersection(barcodes)]
    if expanded_only and "clone_size" in obs.columns:
        obs = obs[obs["clone_size"] > 1]
    if obs.empty:
        return []

    species_col = next((c for c in epitope_cols if "species" in c), None)
    epitope_col = next((c for c in epitope_cols if "epitope" in c), None)
    if epitope_col is None:
        return []

    # Drop NaN and scirpy's literal "ambiguous" (CDR3 matched multiple
    # conflicting epitopes) — neither is a single predicted specificity.
    hits = obs.dropna(subset=[epitope_col])
    hits = hits[hits[epitope_col].astype(str).str.lower() != "ambiguous"]
    if hits.empty:
        return []

    counts = {}
    for _, row in hits.iterrows():
        species = row[species_col] if species_col else "unknown"
        label = f"{species}: {row[epitope_col]}"
        counts[label] = counts.get(label, 0) + 1
    return [
        f"{label} ({n} cell{'s' if n > 1 else ''})"
        for label, n in sorted(counts.items(), key=lambda kv: -kv[1])
    ]
