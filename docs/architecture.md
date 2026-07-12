# Clonal Compass — architecture

Data flows one direction: both datasets converge on a single loader/registry,
run through the *same* dataset-agnostic pipeline (a filename `suffix` keeps their
outputs side by side), and the Streamlit chat reads only the compact ~20 KB
evidence bundle — never the multi-hundred-MB `.h5ad` objects.

```mermaid
flowchart TD
    subgraph sources["Data sources"]
        PBMC["PBMC · 10x 5' GEX+VDJ<br/>data/raw/*.h5 + *.csv"]
        CANCER["Cancer · Wu 2020 TILs<br/>scirpy wu2020_3k (cached)"]
    end

    REG["io.load_dataset(key)<br/>DATASETS registry → DatasetSpec<br/>display_name · suffix · loader · marker_set"]
    PBMC --> REG
    CANCER --> REG

    MARKERS["markers.MARKER_SETS<br/>pbmc | cancer signatures"]

    subgraph pipe["run_pipeline.py --dataset pbmc/cancer"]
        QC[qc.run_qc] --> CLU[cluster.normalize_and_cluster]
        CLU --> ANN["annotate.annotate_clusters<br/>(marker set per spec)"]
        ANN --> CLO["clonal.compute_clonal_expansion<br/>+ merge_clone_size"]
    end
    REG -->|"gex, tcr, spec"| QC
    MARKERS -.->|"selected by spec.marker_set"| ANN

    PROC["data/processed/*{suffix}.h5ad<br/>figures/umap_*{suffix}.png"]
    CLO --> PROC

    subgraph rpt["generate_report.py --dataset …"]
        BUILD["report.build_report_data<br/>rank markers · epitopes · notable · interpret"]
    end
    PROC --> BUILD
    VDJ["epitope.py · VDJdb (cached)"] -.->|"epitope hits"| BUILD
    INT["interpret.py<br/>Claude API / deterministic fallback"] -.->|"hedged text"| BUILD

    ART["reports/*{suffix}<br/>report.html · interpretations.md · cluster_evidence.json"]
    BUILD --> ART

    subgraph chatui["Streamlit chat · app.py"]
        SEL["sidebar dataset selector<br/>discovers cluster_evidence*.json"]
        ASK["chat.ask_question<br/>tool-use loop over evidence bundle"]
        SEL --> ASK
    end
    ART -->|"cluster_evidence*.json only (~20 KB)"| SEL
    INT -.->|"same guardrails · Claude / fallback"| ASK

    classDef ds fill:#e6f0ff,stroke:#3b6fb0;
    classDef out fill:#eef7ee,stroke:#3a8a3a;
    class PBMC,CANCER ds;
    class PROC,ART out;
```

## How to read it

- **Two datasets, one path.** `PBMC` and `Cancer` both resolve through
  `io.load_dataset` into a single `DatasetSpec`; the pipeline stages
  (`qc → cluster → annotate → clonal`) are identical for both. The spec's
  `suffix` (`""` for PBMC, `_cancer` for cancer) is what keeps their processed
  objects, figures, and reports coexisting on disk.
- **Marker set is injected, not branched.** `annotate.annotate_clusters` scores
  whichever `markers.MARKER_SETS[spec.marker_set]` the dataset selects — the
  labelling logic itself is unchanged and dataset-agnostic.
- **The chat is bundle-only.** `app.py` discovers the per-dataset
  `cluster_evidence*.json` bundles, lets the user pick one, and answers grounded
  in that ~20 KB JSON via `chat.ask_question`'s tool-use loop. It never opens the
  `.h5ad` objects.
- **Solid arrows** are data written/read; **dotted arrows** are cross-cutting
  helpers (marker signatures, VDJdb epitopes, the Claude/​fallback interpretation
  layer with its shared guardrails).

> Rendered by any Mermaid-aware viewer (GitHub, VS Code, mermaid.live). To export
> an image: `mmdc -i docs/architecture.md -o docs/architecture.svg` (needs
> `@mermaid-js/mermaid-cli`).
