# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Clonal Compass

A single-cell immune repertoire co-pilot, built for Built with Claude: Life Sciences (Gladstone / Cerebral Valley hackathon).

## Current status

Days 1–4 of the plan (below) are **done and verified end-to-end**:

- **Day 1** — venv + pinned stack; the 10x demo data loads cleanly.
- **Day 2** — core pipeline (QC → cluster → annotate → clonal expansion) with UMAPs.
- **Day 3** — Claude interpretation layer + VDJdb epitope cross-referencing.
  Both interpretation paths are verified working: the **Claude API path runs
  live** (needs `ANTHROPIC_API_KEY`), and a deterministic, guardrail-compliant
  **fallback** runs when no key is set. Each interpretation records which path
  produced it (`Interpretation.source` → `"Claude API"` / `"deterministic fallback"`).
- **Day 4** — self-contained HTML + Markdown report, plus a compact per-cluster
  **evidence bundle** (`reports/cluster_evidence.json`, ~20 KB).

**Second dataset added (additive, verified):** the same pipeline now runs on
scirpy's built-in **"3k T cells from cancer"** (Wu et al. 2020, tumor-infiltrating
T cells across lung/renal/colon/endometrial). Select with `--dataset {pbmc,cancer}`
on both pipeline scripts (default `pbmc`); each dataset's artifacts get a filename
suffix (`""` for PBMC, `_cancer` for cancer) so both coexist. No stage logic
changed. Expansion is decisively different — cancer: **22.7%** of TCR+ cells in
expanded clones (largest clone 20); PBMC baseline: **2.8%** (largest clone 9) —
i.e. real tumor-driven clonal expansion, not healthy-donor noise. See
"Dataset selection" below.

**Active stretch goal:** a Streamlit chat interface (see "Chat interface" below).
**Stretch goal — DONE:** the Streamlit chat interface is built and verified
(see "Chat interface" below). It grounds every answer in the evidence bundle
via Claude's tool-use loop, tags each answer with its source path, and now
carries conversation history so follow-up references ("that cluster", "it")
resolve across turns.

## Architecture

Data flows one direction. Each analysis stage is a re-runnable function in the
`clonal_compass/` package; scripts in `scripts/` orchestrate them:

```
raw 10x files ──run_pipeline.py──▶ data/processed/*.h5ad + figures/*.png
data/processed ──generate_report.py──▶ reports/clonal_compass_report.html
                                       reports/cluster_interpretations.md
                                       reports/cluster_evidence.json
```

> A rendered **Mermaid architecture diagram** (both datasets → shared registry →
> pipeline → report → chat) lives in [`docs/architecture.md`](docs/architecture.md).
> Keep it in sync when the data flow changes.

> **`docs/` is also the GitHub Pages site** (`twhetzel.github.io/clonal-compass/`).
> Unlike the gitignored working `/reports/`, `docs/index.html` and the **committed
> snapshots** under `docs/reports/` (`clonal_compass_report.html`,
> `clonal_compass_report_cancer.html`) are the *published* example reports the
> README links to. Regenerate with `generate_report.py` and copy the fresh HTML
> into `docs/reports/` when you want the public links to reflect new output.

**Package (`clonal_compass/`):**
- `io` — load GEX (`sc.read_10x_h5`) + TCR (`ir.io.read_10x_vdj`). Also a small
  **dataset registry** (`DatasetSpec`, `DATASETS`, `load_dataset(key)`) that
  returns a uniform `(gex, tcr, spec)` for either `pbmc` (files on disk) or
  `cancer` (scirpy's `wu2020_3k()` MuData, split into GEX + AIRR AnnData). `spec`
  carries the display name + filename `suffix` the scripts use to keep both
  datasets' outputs side by side.
- `qc` — standard QC filtering (thresholds are function params, not hard-coded)
- `cluster` — normalize → HVG → PCA → UMAP → Leiden (`flavor="igraph"`; leidenalg is **not** installed)
- `markers` + `annotate` — marker-gene signatures → per-cluster lineage / T-subset labels via `sc.tl.score_genes`. Signatures live in a **registry** (`markers.MARKER_SETS`: `pbmc` broad-lineage/subset panel, `cancer` TIL-state panel) selected per dataset via `DatasetSpec.marker_set`. `annotate_clusters()` takes the two signature dicts as params (default = PBMC); its scoring/labelling logic is unchanged and marker-set-agnostic.
- `clonal` — clonotypes by **CDR3 amino-acid identity** (`ir.tl.define_clonotype_clusters`, aliased to `clone_id`) → clone size + expansion, merged onto GEX by barcode
- `plots` — UMAP PNGs (captioned)
- `interpret` — `ClusterEvidence`, `Interpretation`, `interpret_cluster()` (Claude API + fallback). Guardrails live in `SYSTEM_PROMPT`. Model `claude-opus-4-8`, adaptive thinking, effort high.
- `epitope` — VDJdb load/annotate + `hits_for_cells`
- `report` — `build_report_data()` (rank markers + annotate epitopes + pick notable clusters + interpret), `render_html` / `render_markdown`, `write_evidence_json`
- `_warnings.silence_demo_warnings()` — call first in every entry script for a clean demo console

**Scripts (`scripts/`):** `load_data.py` (Day-1 sanity), `run_pipeline.py`
(Day-2 pipeline), `interpret_clusters.py` (Day-3 console inspection),
`generate_report.py` (Day-4 artifacts + evidence bundle).

**Run the whole thing:**
```
.venv/bin/python scripts/run_pipeline.py       # → data/processed + figures  (PBMC, default)
.venv/bin/python scripts/generate_report.py    # → reports/*                  (PBMC, default)
.venv/bin/python scripts/run_pipeline.py    --dataset cancer   # Wu 2020 tumor T cells
.venv/bin/python scripts/generate_report.py --dataset cancer   # → reports/*_cancer.*
```
Set `ANTHROPIC_API_KEY` first for live Claude interpretations (else fallback).
Outputs are all gitignored: `data/processed/`, `figures/`, `reports/`, and
`data/reference/vdjdb.h5ad` (cached VDJdb reference). (The published example
reports under `docs/reports/` are the exception — those are committed snapshots;
see the GitHub Pages note above.)

### Dataset selection
- Both `run_pipeline.py` and `generate_report.py` take `--dataset {pbmc,cancer}`
  (default `pbmc`). The choice flows through `io.load_dataset()`; **no stage
  function changed** — QC → cluster → annotate → clonal → evidence bundle are
  identical for both.
- Artifacts are suffixed so datasets coexist: PBMC keeps the original unsuffixed
  names (`gex_annotated.h5ad`, `cluster_evidence.json`, `umap_*.png`, …); cancer
  writes the `_cancer` variants (`cluster_evidence_cancer.json`, etc.).
- **Per-dataset marker sets.** The cancer path uses a **TIL-state signature
  panel** (`markers.CANCER_*`), not the PBMC lineage panel. Two design points
  worth knowing:
  - Its lineage layer is a **single `T cell` entry** on purpose. Every cell in
    this dataset is a T cell, and `score_genes` measures enrichment *against a
    similarly-expressed background*, so in a pure-T population the pan-T score
    has no contrast and lands near zero — a multi-lineage argmax then gets
    decided by noise and mislabels real T cells as NK/B/myeloid (the old
    behavior). Asserting the one true lineage routes every cluster into the
    subset panel, where the real discrimination happens.
  - The subset panel resolves tumor-relevant states the PBMC panel can't:
    cytotoxic effector, exhausted/terminal, Treg, Tfh (CXCL13+), tissue-resident
    memory, proliferating. Validated against canonical markers (Treg→FOXP3,
    Tfh→CXCL13, TRM→ZNF683/ITGAE, cytotoxic→GZMB/GNLY, cycling→MKI67/TOP2A).
    `STMN1` is deliberately **excluded** from the cycling signature — it's broadly
    expressed in activated effectors and inflates the score for non-dividing cells.
  - Nice consequence: the most clonally-expanded clusters now correctly label as
    **CD8 cytotoxic effector** (they were mislabeled "NK" under the PBMC panel).
  - Weakly-differentiated clusters can still land on a narrow-margin call (an
    inherent property of score-based argmax, same as the PBMC path); the hedged
    report + per-cluster `rank_genes_groups` are the backstop.
- `wu2020_3k` has **0 `MT-` genes**, so the `max_pct_mt` QC filter is a harmless
  no-op on it (2,992/3,000 cells kept on default thresholds).

### Gotchas & decisions already made
- **Clonotypes are amino-acid based.** We switched from exact-nt
  `define_clonotypes` to `define_clonotype_clusters` (CDR3 aa identity) — more
  defensible for expansion claims (won't over-count silent mutations).
- **scirpy's `ir.datasets.vdjdb()` is broken** on current VDJdb releases (reads
  `vdjdb_full.txt` from the archive root, but it's now nested in a dated
  subdir). `epitope.load_vdjdb()` reimplements the loader with a recursive glob
  and caches to `data/reference/vdjdb.h5ad`.
- **Epitope matches use `strategy="unique-only"`**; scirpy's literal
  `"ambiguous"` (CDR3 matched multiple epitopes) is treated as "no match".
- **Notable clusters are chosen by number of expanded cells**, not percent.
- The VDJdb layer genuinely finds public epitopes (GILGFVFTL flu, NLVPMVATV
  CMV, LLWNGPMAV YFV) — but in this healthy donor they're all in singleton
  clones, so the expanded-clone report correctly shows "no match".

## Chat interface (DONE)

A lightweight **Streamlit** app for asking natural-language questions about the
clusters/clones. Built and verified end-to-end.

**Run it:**
```
.venv/bin/streamlit run app.py
```
Set `ANTHROPIC_API_KEY` first for live Claude answers (else the deterministic
fallback runs). `streamlit==1.59.1` is pinned in `requirements.txt`.

**Pieces:**
- `clonal_compass/chat.py` — the grounded Q&A layer. `load_evidence_bundle()`
  reads the JSON; `ask_question(question, evidence_bundle, history=None)` runs
  Claude's tool-use loop and returns an `Interpretation(text, source)` just like
  `interpret.interpret_cluster` (source → `"Claude API"` / `"deterministic
  fallback"`).
- `scripts/ask.py` — one-shot CLI to ask a question from the terminal.
- `app.py` — the Streamlit app, now a **single fused surface: report visuals +
  chat**. The **sidebar** holds the dataset `selectbox`, the dataset **stat
  chips** (cells / clusters / TCR+ / VDJdb-matched), and a **tabbed "Repertoire
  visuals" panel**: the four Day-2 UMAPs (loaded as base64 data-URIs straight
  from `figures/`) plus an **Evidence tab** that re-renders the HTML report's
  source-tagged cluster cards + dataset-wide epitope overview natively. The
  **main pane** is a clean chat under a **sticky hero header** (title + "not a
  diagnostic tool" disclaimer + which dataset is active), with example-question
  hints and a per-answer **source badge** (green = `Claude API`, amber =
  `deterministic fallback`). Nothing is recomputed — UMAP PNGs and every
  stat/interpretation/epitope hit come from disk; bundles + figure data-URIs are
  `@st.cache_data`-cached by path. Styling is **ported from `report.py` under a
  `.cc` scope** so the live app and the static report read as one product, and
  Streamlit's dev chrome (menu/footer/header/toolbar) is hidden via CSS for a
  clean demo.
- **Branding.** The app's logo lives in `assets/`: `clonal_compass_logo.png`
  (512×512) renders inline next to the "Clonal Compass" title in the sidebar
  identity block (base64 data-URI via the same `_img_data_uri` helper as the
  UMAPs), and `assets/favicon.png` (a 64×64 downscale of the logo) is the
  browser-tab `page_icon`. Both fall back to the 🧭 emoji if the files are
  missing. `assets/` is committed (not gitignored).
- **Example-question copy buttons.** The hero's "Try:" hints each carry an inline
  copy icon; a small `components.html` script (using `execCommand('copy')` on the
  parent document, so no clipboard-write permission is needed) copies the plain
  question text on click and flashes a ✓.
- **`.streamlit/config.toml`** — a `light` theme keyed to the report's design
  tokens (`--accent #2b6cb0`, `--bg #f7f8fb`, `--card #fff`, `--fg #1a1f2b`) so
  the app matches the report out of the box, plus `toolbarMode = "minimal"` to
  drop Streamlit's dev toolbar in recorded demos.

**Multi-dataset selection.** `app.py` discovers which `reports/cluster_evidence{suffix}.json`
bundles exist (one per `io.DATASETS` entry) and offers them in a sidebar
`selectbox` labelled by `display_name`. The active dataset's `name` + stats show
in the sidebar, and the sticky hero header names the active dataset, so it's
always clear which dataset is loaded. Because each chat is grounded in ONE
dataset, **switching datasets clears the conversation** (`st.session_state.messages`)
— otherwise a follow-up like "that cluster" would resolve against the wrong
dataset's clusters. Discovery is uncached (cheap `stat()`) so a newly-generated
bundle appears on the next rerun.

**Design constraints (all honored):**
- **Reads only the evidence bundle(s)** (`reports/cluster_evidence*.json`) — never
  opens the `.h5ad` objects (those are hundreds of MB; each JSON is ~20 KB and
  loads instantly). Regenerate with `generate_report.py [--dataset …]`.
- Uses the Anthropic **tool-calling** loop: `chat.TOOLS` exposes
  `get_dataset_overview`, `list_clusters`, `get_cluster_evidence`,
  `get_epitope_matches` so Claude grounds every answer in the data rather than
  free-associating. Bounded by `_MAX_TOOL_ROUNDS = 6`.
- **Epitope tool has two scopes — don't confuse them.** `get_epitope_matches`
  called **without** `cluster_id` returns the **dataset-wide** matches (all TCR+
  cells: the total count + every predicted specificity); called **with** a
  `cluster_id` it returns only that cluster's **expanded-clone** hits (a strict
  subset). A general "are there any epitope matches?" must use the no-cluster
  form — the most-expanded clusters often have zero expanded-clone hits even when
  the dataset has many (e.g. cancer: clusters 4/5 empty, but 42 matched cells
  dataset-wide). The per-cluster empty note and a system-prompt line both warn
  against generalizing a single cluster's empty list to the whole dataset. (This
  was a real bug: the tool used to be per-cluster-only, so the chat wrongly
  answered "no epitope matches".)
- **Reuses the exact same guardrails** — `CHAT_SYSTEM_PROMPT` is
  `interpret.SYSTEM_PROMPT` + tool-use operating instructions only; no second,
  looser prompt. Same model/effort as `interpret.py`.
- Additive only: does not modify the working Day-1–4 pipeline.

**Conversation history.** `ask_question` takes an optional `history` — the prior
turns as `[{"role": "user"|"assistant", "content": str}, ...]`, **not** including
the current question — and seeds the tool-use loop with them so follow-ups like
"what about that cluster's marker genes?" resolve "that"/"it" to what was
discussed. Only the plain **text** of each prior turn is carried forward, not the
earlier tool_use/tool_result blocks (enough to resolve references, and avoids
re-pairing tool blocks across turns). `app.py` passes
`st.session_state.messages` (captured *before* appending the new question) as
`history`. The deterministic fallback ignores history — without the model it
can't resolve references, and it says so.

## Environment & commands

Day 1 is done: a dedicated venv exists with a pinned, known-good stack.

- **Python:** 3.11.9 (via pyenv). The system `python3` is 3.8 and too old —
  always use the venv.
- **Setup / reinstall:**
  `~/.pyenv/versions/3.11.9/bin/python3.11 -m venv .venv && .venv/bin/python -m pip install -r requirements.txt && .venv/bin/python -m pip install -e .`
- **Run anything in the env:** `.venv/bin/python <script>` (or `source .venv/bin/activate`).
  The `pip install -e .` step makes `clonal_compass` importable — **no `PYTHONPATH=.` prefix needed**.
- **Pinned deps** live in `requirements.txt`: numpy 1.26.4, pandas 2.2.3,
  anndata 0.10.9, scanpy 1.10.3, scirpy 0.17.2, anthropic 0.116.0. `pyproject.toml`
  reads its dependencies from `requirements.txt`, so it stays the single pin source.
  Don't let these auto-upgrade.
- No test/lint tooling yet — add commands here when introduced.

### Apple Silicon gotcha (parasail)
`scirpy` depends on `parasail`, which has **no arm64 macOS wheel** and must
build from source. That build needs GNU autotools — install once with
`brew install automake libtool` (`autoconf` was already present). Without
them the `pip install` fails on `parasail` with "autoreconf -fi failed".

## Data

Two datasets are wired up (select with `--dataset`; see "Dataset selection").

**1. `pbmc` (default)** — 10x Genomics **10k Human PBMCs, 5' v2 (GEX + VDJ)**
demo dataset (healthy donor). Files live in `data/raw/` (gitignored). Re-download
with:

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

**2. `cancer`** — scirpy's built-in **"3k T cells from cancer"** (Wu et al.
2020): 3,000 tumor-infiltrating T cells downsampled across lung/renal/colon/
endometrial patients (Tumor/NAT/Blood). No download step — fetched and cached
by scirpy via `ir.datasets.wu2020_3k()` (a MuData with `gex` raw-count +
`airr` modalities; `io._load_cancer()` splits them). All 3,000 cells are
TCR-containing and already barcode-paired.

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

> **Status: Days 1–4 are all complete** (see "Current status" up top). The
> per-day plan and starter prompts below are retained as reference for what was
> asked and why. New work is the chat-interface stretch goal, not these days.

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
