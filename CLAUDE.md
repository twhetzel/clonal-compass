# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Clonal Compass

A single-cell immune repertoire co-pilot, built for Built with Claude: Life Sciences (Gladstone / Cerebral Valley hackathon).

## Environment & commands

Day 1 is done: a dedicated venv exists with a pinned, known-good stack.

- **Python:** 3.11.9 (via pyenv). The system `python3` is 3.8 and too old —
  always use the venv.
- **Setup / reinstall:**
  `~/.pyenv/versions/3.11.9/bin/python3.11 -m venv .venv && .venv/bin/python -m pip install -r requirements.txt`
- **Run anything in the env:** `.venv/bin/python <script>` (or `source .venv/bin/activate`).
- **Pinned deps** live in `requirements.txt`: numpy 1.26.4, pandas 2.2.3,
  anndata 0.10.9, scanpy 1.10.3, scirpy 0.17.2. Don't let these auto-upgrade.
- No test/lint tooling yet — add commands here when introduced.

### Apple Silicon gotcha (parasail)
`scirpy` depends on `parasail`, which has **no arm64 macOS wheel** and must
build from source. That build needs GNU autotools — install once with
`brew install automake libtool` (`autoconf` was already present). Without
them the `pip install` fails on `parasail` with "autoreconf -fi failed".

## Data

10x Genomics **10k Human PBMCs, 5' v2 (GEX + VDJ)** demo dataset. Files live
in `data/raw/` (gitignored). Re-download with:

```
curl -o data/raw/sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5 \
  https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.h5
curl -o data/raw/sc5p_v2_hs_PBMC_10k_t_filtered_contig_annotations.csv \
  https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_t_filtered_contig_annotations.csv
```

Note the VDJ contig CSV lives under the `.../sc5p_v2_hs_PBMC_10k/` sample
directory (not a `_t` subdirectory — that path 403s).

`scripts/load_data.py` is the Day 1 sanity check: loads GEX via
`sc.read_10x_h5`, TCR via `ir.io.read_10x_vdj`, and pairs them into a
`MuData`. Expected output: 10,548 cells × 36,601 genes, 5,328 TCR contigs,
5,308 paired (~50%). Note: scirpy's old `merge_with_ir` was removed in
v0.13 — use `MuData({"gex": ..., "airr": ...})`.

## Project brief

Input: paired scRNA-seq + TCR-seq data (10x Genomics 5' GEX + VDJ format).

Pipeline: run QC, clustering, and T-cell subset annotation on the GEX data;
call clonotypes and compute clonal expansion metrics from the VDJ data; merge
the two; use Claude to interpret expanded/notable clusters, cross-reference
against TCR-epitope databases (VDJdb, IEDB), and write a plain-language,
appropriately hedged report a working immunologist could actually use.

Output: a reproducible notebook/script + a short written report (and,
if time allows, a lightweight dashboard).

## Timeframe: 4 days. Hard scope discipline required.

- **Day 1** — Environment + data loading ONLY. No analysis code.
  This absorbs the single biggest risk (scirpy/Scanpy dependency conflicts),
  so it gets the full day even though the "work" itself is light.

  Starter prompt:
  > "Set up a Python environment for this project using conda or venv.
  > Install scanpy, scirpy, anndata, pandas, and numpy, pinning versions to
  > ones known to work well together. Then download and load 10x Genomics'
  > public 5' GEX + VDJ demo dataset (help me find the right one), and
  > confirm it loads cleanly into an AnnData object. Show me basic dataset
  > info (n_obs, n_vars) as a sanity check."

  **Checkpoint (end of day 1)**: if this isn't clean and working, fall back
  to simpler clone-counting directly from `filtered_contig_annotations.csv`
  instead of fighting scirpy further.

- **Day 2** — Core pipeline: QC → clustering → T-cell subset annotation
  (marker-gene based) → clonal expansion metrics merged onto clusters.

  Starter prompt:
  > "Now build the core analysis pipeline as a set of separate, re-runnable
  > functions (not one monolithic notebook): (1) standard QC filtering,
  > (2) normalization/clustering with Scanpy, (3) T-cell subset annotation
  > using standard marker genes, (4) clonal expansion metrics from the VDJ
  > data merged onto the clusters. Show me a UMAP colored by cluster and by
  > clone size as we go so I can sanity-check each stage."

  Split into two sessions if needed (QC+clustering, then
  annotation+expansion) rather than debugging one huge diff at once.

- **Day 3** — Claude interpretation layer: this is the differentiator —
  protect this time even if day 2 ran long.

  Starter prompt:
  > "Write a function that takes a cluster's marker genes and clonal
  > expansion stats and produces a plain-language interpretation. Follow
  > these guardrails strictly: distinguish 'consistent with X' from 'is X'
  > or 'proves X'; cite the specific metric behind every claim; never state
  > a conclusion more strongly than the stats support. Then add epitope
  > cross-referencing against VDJdb for expanded clones — report matches as
  > 'possible/predicted,' not confirmed, and explicitly note when no match
  > is found as a normal, valid outcome."

- **Day 4** — Report layer + demo prep. No new features.

  Starter prompt:
  > "Generate a short markdown or HTML report combining the plots and
  > Claude's interpretations for 2-3 interesting clusters/clones from our
  > dataset. Keep it clean and readable — this is the main thing I'll show
  > judges. Don't add new analysis, just present what we already have well."

  Reserve the last few hours purely for demo rehearsal and confirming the
  whole thing runs end-to-end from a clean state — a broken live demo is a
  much bigger loss than a missing feature.

**Explicit fallback rule**: if you're not through the Day 1 checkpoint by
end of day 1, or Day 2's pipeline isn't producing sane clusters by midday
on day 3, that's the trigger to cut scope (drop epitope matching, keep just
cluster interpretation) rather than push the whole timeline back.

Cut entirely at this timeframe (not stretch goals, just off the table):
HLA typing, neoantigen matching, fancy interactive dashboards.

## Environment

- Python env dedicated to this project (venv or conda), not shared with
  anything else.
- Core deps: `scanpy`, `scirpy`, `anndata`, `pandas`, `numpy`. Pin versions
  explicitly once a working combination is found — don't let anything
  auto-upgrade mid-week.
- Data: start with 10x Genomics' public 5' GEX + VDJ demo dataset
  (support.10xgenomics.com/single-cell-vdj/datasets). Swap in a richer
  dataset (e.g. a public tumor-infiltrating T-cell dataset) only after the
  pipeline is solid, if time allows.

## Interpretation guardrails (critical — read before writing prompts)

Claude's summaries of biological results must stay scientifically cautious:

- Distinguish "consistent with X" from "proves X" or "is X." Never state a
  cluster/clone conclusion more strongly than the underlying stats support.
- Every claim about a cluster or clone should cite the specific metric
  behind it (expansion fold-change, marker gene expression level, etc.),
  not just an impression.
- Epitope matches from VDJdb/IEDB are a database hit, not proof of
  functional reactivity. Always report them as "predicted/possible match,"
  include a rough confidence framing, and note when no match is found
  (this is common and should be presented as a normal, valid result).
- No patient-level or clinical claims. This is a research/exploratory
  tool, not a diagnostic one — say so explicitly in report output.
- When summarizing, prefer more, shorter hedged claims over fewer,
  confident-sounding ones.

## Conventions

- Keep pipeline steps as separate, re-runnable functions/scripts, not one
  monolithic notebook — makes it easier to swap datasets on day 3–4.
- Every plot should have a one-line caption stating what metric it shows.
- Commit early and often; use `/clear` between major pipeline stages if a
  session gets long, and paste a short "here's where we are" summary after.
