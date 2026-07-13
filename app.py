"""Chat interface, step 3 — now a single fused surface: report visuals + chat.

Discovers which per-dataset evidence bundles exist under ``reports/``
(``cluster_evidence{suffix}.json``, one per entry in ``io.DATASETS``), lets the
user pick one, and lays out one page that brings the two demo surfaces together:

  * left  — the Day-2 UMAPs (loaded straight from ``figures/``) plus the same
    source-tagged evidence cards the HTML report shows, re-rendered natively in
    the report's palette;
  * right — the grounded chat, answering each question via ``chat.ask_question``
    against the selected bundle (Claude tool-use loop, deterministic fallback
    when no API key is set).

Nothing is recomputed: the UMAP PNGs already exist on disk and every stat/
interpretation/epitope hit is read from the ~20 KB evidence bundle. The styling
is ported from ``clonal_compass/report.py`` so the live app and the static
report read as one product. Switching datasets clears the conversation, since
each chat is grounded in a single dataset.

Run:
  .venv/bin/streamlit run app.py

Regenerate a bundle (and its figures) with:
  scripts/generate_report.py [--dataset pbmc|cancer]
"""

from __future__ import annotations

import base64
import html
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from clonal_compass._warnings import silence_demo_warnings

silence_demo_warnings()

from clonal_compass import chat, io  # noqa: E402

ROOT = Path(__file__).resolve().parent
REPORTS_DIR = ROOT / "reports"
FIGURES_DIR = ROOT / "figures"

# The Day-2 UMAPs, in display order: (base filename, short tab label, caption).
# Mirrors report.FIGURES; kept local so the lightweight chat app never imports
# the scanpy-heavy report module.
FIGURES = [
    ("umap_clusters.png", "Clusters", "Leiden clusters on the transcriptome."),
    ("umap_celltype.png", "Cell types", "Marker-based cell-type / T-cell subset labels."),
    ("umap_clonesize.png", "Clone size", "TCR clone size (log1p of cells per clonotype)."),
    ("umap_expansion.png", "Expansion", "Clonal-expansion category per cell."),
]


# --------------------------------------------------------------------------- #
# Bundle discovery (unchanged behaviour)
# --------------------------------------------------------------------------- #
def _bundle_path(spec: io.DatasetSpec) -> Path:
    """Path to a dataset's evidence bundle (mirrors generate_report.py's suffix)."""
    return REPORTS_DIR / f"cluster_evidence{spec.suffix}.json"


def _available_datasets() -> dict[str, Path]:
    """Registry keys -> bundle path, restricted to bundles that exist on disk."""
    return {
        key: _bundle_path(spec)
        for key, spec in io.DATASETS.items()
        if _bundle_path(spec).exists()
    }


@st.cache_data(show_spinner=False)
def _load_bundle(path: str) -> dict:
    """Load one compact evidence bundle (cached per path across reruns)."""
    return chat.load_evidence_bundle(path)


@st.cache_data(show_spinner=False)
def _img_data_uri(path: str) -> str:
    """Base64 data-URI for a figure PNG (cached; mirrors report._img_data_uri)."""
    return "data:image/png;base64," + base64.b64encode(Path(path).read_bytes()).decode()


def _figures_for(suffix: str) -> list[tuple[str, str, str]]:
    """Resolve (data_uri, tab_label, caption) for each figure that exists on disk."""
    out: list[tuple[str, str, str]] = []
    for name, label, caption in FIGURES:
        p = Path(name)
        fp = FIGURES_DIR / f"{p.stem}{suffix}{p.suffix}"
        if fp.exists():
            out.append((_img_data_uri(str(fp)), label, caption))
    return out


# --------------------------------------------------------------------------- #
# HTML fragments — ported from report.py's palette, scoped under `.cc`
# --------------------------------------------------------------------------- #
def _badge_html(source: str) -> str:
    kind = "claude" if source == "Claude API" else "fallback"
    return f'<span class="cc-badge cc-badge-{kind}">{html.escape(source)}</span>'


def _figure_card(data_uri: str, caption: str) -> str:
    return (
        f'<div class="cc"><figure class="cc-fig">'
        f'<img src="{data_uri}" alt="{html.escape(caption)}">'
        f"<figcaption>{html.escape(caption)}</figcaption></figure></div>"
    )


def _stat_chips(ds: dict) -> str:
    def fmt(v):
        return f"{v:,}" if isinstance(v, int) else "?"

    chips = [
        ("Cells", fmt(ds.get("n_cells"))),
        ("Clusters", fmt(ds.get("n_clusters"))),
        ("TCR+ cells", fmt(ds.get("n_tcr_cells"))),
        ("VDJdb-matched", fmt(ds.get("n_tcr_epitope_matched"))),
    ]
    inner = "".join(
        f'<div class="cc-stat"><span class="cc-stat-num">{v}</span>'
        f'<span class="cc-stat-lab">{html.escape(k)}</span></div>'
        for k, v in chips
    )
    return f'<div class="cc"><div class="cc-stats">{inner}</div></div>'


def _evidence_card(cluster: dict) -> str:
    interp = cluster.get("interpretation") or {}
    source = interp.get("source", "deterministic fallback")
    markers = ", ".join(
        f"{html.escape(g)} ({s:+.2f})" for g, s in cluster.get("top_markers", [])[:6]
    )
    clones = ", ".join(str(c) for c in cluster.get("top_clone_sizes", [])) or "none > 1"
    rows = [
        ("Cells", f"{cluster.get('n_cells', 0):,}"),
        (
            "TCR coverage",
            f"{cluster.get('n_with_tcr', 0):,} / {cluster.get('n_cells', 0):,} "
            f"({cluster.get('pct_with_tcr', 0):.0f}%)",
        ),
        (
            "Clonal expansion",
            f"{cluster.get('n_expanded_cells', 0):,} cells in clones &gt; 1 "
            f"({cluster.get('pct_expanded', 0):.0f}% of TCR+)",
        ),
        ("Top markers", markers or "&mdash;"),
        ("Largest clones", clones),
    ]
    row_html = "".join(f"<tr><th>{lab}</th><td>{val}</td></tr>" for lab, val in rows)

    interp_text = interp.get("text", "")
    paras = "".join(
        f"<p>{html.escape(p.strip())}</p>" for p in interp_text.split("\n") if p.strip()
    )

    hits = cluster.get("epitope_hits") or []
    if hits:
        items = "".join(f"<li>{html.escape(h)}</li>" for h in hits)
        epi = (
            '<div class="cc-epi cc-epi-hit"><strong>Predicted epitope matches '
            "(VDJdb, exact CDR3)</strong> &mdash; predicted specificity, not "
            f"confirmed reactivity:<ul>{items}</ul></div>"
        )
    else:
        epi = (
            '<div class="cc-epi cc-epi-none"><strong>No VDJdb match</strong> for '
            "this cluster's expanded clones. Absence of a match is a normal, valid "
            "outcome &mdash; not a failure.</div>"
        )

    return f"""<div class="cc"><section class="cc-card">
      <div class="cc-card-head">
        <h3>Cluster {html.escape(str(cluster.get('cluster_id', '?')))}:
            {html.escape(str(cluster.get('cell_type', 'Unknown')))}</h3>
        {_badge_html(source)}
      </div>
      <table class="cc-evidence">{row_html}</table>
      <div class="cc-interp">
        <div class="cc-interp-label">Interpretation
          <span class="cc-src">generated by: {html.escape(source)}</span>
        </div>{paras or '<p>&mdash;</p>'}
      </div>
      {epi}
    </section></div>"""


def _epitope_overview(bundle: dict) -> str:
    ds = bundle.get("dataset", {})
    overview = bundle.get("epitope_overview") or []
    n_matched = ds.get("n_tcr_epitope_matched", 0)
    if overview:
        items = "".join(f"<li>{html.escape(h)}</li>" for h in overview)
        body = (
            f"<p>{n_matched} TCR+ cell(s) carry a CDR3 that exactly matches a VDJdb "
            "entry. These are <em>predicted</em> specificities (database hits), not "
            f"confirmed reactivity:</p><ul class='cc-overview'>{items}</ul>"
        )
    else:
        body = (
            "<p>No exact VDJdb CDR3 matches were found dataset-wide &mdash; a normal, "
            "valid outcome for a healthy-donor repertoire.</p>"
        )
    return f'<div class="cc"><section class="cc-card">{body}</section></div>'


# Injected once. Ports the report tokens/classes under a `.cc` scope, restyles
# the native chat bubbles, and hides Streamlit's menu/footer/header chrome.
_STYLE = """
<style>
:root{--fg:#1a1f2b;--muted:#5b6472;--line:#e3e7ee;--bg:#f7f8fb;--card:#fff;
--accent:#2b6cb0;--claude:#1f7a4d;--claude-bg:#e6f4ec;--fallback:#9a6a00;--fallback-bg:#fdf2dc}

/* --- hide Streamlit dev chrome for a clean demo --- */
#MainMenu,footer,header[data-testid="stHeader"]{visibility:hidden;height:0}
[data-testid="stToolbar"]{display:none}
.block-container{padding-top:1.4rem;padding-bottom:6rem;max-width:900px}

/* --- scoped report styles (only affect our injected HTML) --- */
.cc,.cc *{box-sizing:border-box}
.cc h3{margin:0;font-size:18px;color:var(--fg)}
.cc-stats{display:flex;gap:12px;flex-wrap:wrap;margin:2px 0 6px}
.cc-stat{background:var(--card);border:1px solid var(--line);border-radius:10px;
padding:10px 16px;min-width:110px}
.cc-stat-num{display:block;font-size:21px;font-weight:700;color:var(--fg)}
.cc-stat-lab{display:block;font-size:11px;color:var(--muted);text-transform:uppercase;
letter-spacing:.04em}
.cc-fig{margin:0;background:var(--card);border:1px solid var(--line);border-radius:12px;
padding:12px}
.cc-fig img{max-width:100%;max-height:340px;width:auto;height:auto;border-radius:6px;display:block;margin:0 auto}
.cc-fig figcaption{font-size:13px;color:var(--muted);margin-top:8px}
.cc-card{background:var(--card);border:1px solid var(--line);border-radius:12px;
padding:18px 20px;margin:14px 0}
.cc-card-head{display:flex;align-items:center;justify-content:space-between;gap:12px;
flex-wrap:wrap}
.cc-badge{font-size:12px;font-weight:700;padding:4px 10px;border-radius:999px;white-space:nowrap}
.cc-badge-claude{color:var(--claude);background:var(--claude-bg)}
.cc-badge-fallback{color:var(--fallback);background:var(--fallback-bg)}
.cc-evidence{width:100%;border-collapse:collapse;margin:12px 0;font-size:14px}
.cc-evidence th{text-align:left;color:var(--muted);font-weight:600;width:180px;
vertical-align:top;padding:6px 12px 6px 0;border-bottom:1px solid var(--line)}
.cc-evidence td{padding:6px 0;border-bottom:1px solid var(--line);color:var(--fg)}
.cc-interp{margin:12px 0}
.cc-interp-label{font-weight:700;font-size:13px;text-transform:uppercase;
letter-spacing:.04em;color:var(--accent)}
.cc-src{font-weight:600;color:var(--muted);text-transform:none;letter-spacing:0;margin-left:8px}
.cc-interp p{margin:8px 0;color:var(--fg)}
.cc-epi{border-radius:10px;padding:12px 14px;font-size:14px;margin-top:10px}
.cc-epi-hit{background:#eef4fb;border:1px solid #cfe0f4}
.cc-epi-none{background:#f2f4f7;border:1px solid var(--line);color:var(--muted)}
.cc-epi ul,.cc-overview{margin:8px 0 0;padding-left:20px}

/* --- product-y header band --- */
.cc-hero h1{margin:0 0 2px;font-size:30px;letter-spacing:-.02em;color:var(--fg)}
.cc-hero p{margin:0 0 4px;color:var(--muted);font-size:15px;max-width:80ch}
.cc-panel-label{font-size:12px;font-weight:700;text-transform:uppercase;
letter-spacing:.05em;color:var(--muted);margin:2px 0 8px}

/* --- header hero + panel labels --- */
.cc-hero h1{margin:0 0 2px;font-size:28px}

/* --- widen the sidebar so the UMAP reads well; the sidebar stays fixed and
       scrolls on its own natively, so the visuals sit still while the chat
       scrolls — no sticky / calc / overflow hacks needed --- */
section[data-testid="stSidebar"]{width:440px!important;min-width:440px!important}
[data-testid="stSidebar"] .cc-fig{padding:8px}
[data-testid="stSidebar"] .cc-fig img{max-height:300px}
[data-testid="stSidebar"] .cc-stats{gap:8px}
[data-testid="stSidebar"] .cc-stat{min-width:0;flex:1 1 44%;padding:8px 12px}
[data-testid="stSidebar"] .cc-stat-num{font-size:18px}
[data-testid="stSidebar"] .cc-card{padding:14px 16px}

/* --- prompt hints under the "Ask the repertoire" label --- */
.cc-hint{font-size:13px;color:var(--muted);line-height:1.4;margin:-2px 0 12px}
.cc-hint em{font-style:italic}
.cc-ex{white-space:nowrap}
.cc-copy{display:inline-flex;align-items:center;cursor:pointer;color:var(--muted);
margin-left:5px;vertical-align:middle;opacity:.65}
.cc-copy:hover{color:var(--accent);opacity:1}

/* --- sticky main header (hero + label + hint) stays put while the chat
       scrolls beneath it; solid bg + top z-index so messages mask cleanly,
       and a divider line marks where the chat scrolls under --- */
[data-testid="stElementContainer"]:has(.cc-mainhead){position:sticky;top:0;z-index:50;
background:var(--bg);padding-top:.6rem;padding-bottom:26px;border-bottom:1px solid var(--line)}
.cc-mainhead .cc-hero h1{font-size:24px}
.cc-mainhead .cc-hero p{font-size:13.5px}
.cc-mainhead .cc-hint{margin:0}

/* --- restyle native chat bubbles onto the report palette --- */
[data-testid="stChatMessage"]{background:var(--card);border:1px solid var(--line);
border-radius:12px;padding:12px 16px}
</style>
"""


# --------------------------------------------------------------------------- #
# Page
# --------------------------------------------------------------------------- #
st.set_page_config(page_title="Clonal Compass", page_icon="🧭", layout="wide")
st.markdown(_STYLE, unsafe_allow_html=True)

available = _available_datasets()
if not available:
    st.error(
        "No evidence bundles found under `reports/`. Generate at least one with "
        "`scripts/generate_report.py [--dataset pbmc|cancer]`."
    )
    st.stop()

# --- Dataset selector (sidebar) -------------------------------------------- #
st.sidebar.header("Dataset")
dataset_key = st.sidebar.selectbox(
    "Active dataset",
    options=list(available),
    format_func=lambda k: io.DATASETS[k].display_name,
    key="dataset_key",
)
st.sidebar.caption(
    "Set `ANTHROPIC_API_KEY` before launching for live Claude answers; "
    "otherwise a deterministic fallback runs."
)

# Each chat is grounded in ONE dataset, so switching invalidates prior turns.
if st.session_state.get("active_dataset") != dataset_key:
    st.session_state.active_dataset = dataset_key
    st.session_state.messages = []

bundle = _load_bundle(str(available[dataset_key]))
ds = bundle.get("dataset", {})
ds_name = ds.get("name", io.DATASETS[dataset_key].display_name)

# --- Sidebar: dataset stats + visuals -------------------------------------- #
# The sidebar is Streamlit's one natively-fixed, self-scrolling region, so the
# UMAP / evidence sit still while the chat scrolls — no sticky/calc/overflow.
figs = _figures_for(io.DATASETS[dataset_key].suffix)
interp_clusters = [c for c in bundle.get("clusters", []) if c.get("interpretation")]
with st.sidebar:
    st.markdown(_stat_chips(ds), unsafe_allow_html=True)
    st.markdown('<div class="cc-panel-label">Repertoire visuals</div>', unsafe_allow_html=True)
    tab_labels = [label for _, label, _ in figs] + ["Evidence"]
    tabs = st.tabs(tab_labels)
    for tab, (data_uri, _label, caption) in zip(tabs, figs):
        with tab:
            st.markdown(_figure_card(data_uri, caption), unsafe_allow_html=True)
    with tabs[-1]:  # Evidence
        if interp_clusters:
            st.caption("Most clonally-expanded clusters — each tagged with the path "
                       "that produced it.")
            for cluster in interp_clusters:
                st.markdown(_evidence_card(cluster), unsafe_allow_html=True)
        st.markdown(_epitope_overview(bundle), unsafe_allow_html=True)

# --- Main area: hero + a clean, ordinary chat ------------------------------ #
st.session_state.setdefault("messages", [])

_COPY_SVG = (
    '<svg viewBox="0 0 24 24" width="13" height="13" fill="none" stroke="currentColor" '
    'stroke-width="2" stroke-linecap="round" stroke-linejoin="round">'
    '<rect x="9" y="9" width="13" height="13" rx="2"/>'
    '<path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg>'
)


def _example(q: str) -> str:
    """One example question with an inline copy icon (copies the plain text)."""
    return (
        f'<span class="cc-ex"><em>“{html.escape(q)}”</em>'
        f'<span class="cc-copy" data-copy="{html.escape(q, quote=True)}" '
        f'title="Copy question">{_COPY_SVG}</span></span>'
    )


_EXAMPLES = [
    "Which cluster is most clonally expanded, and why?",
    "Do any clones match a known epitope?",
]

st.markdown(
    '<div class="cc cc-mainhead">'
    '<div class="cc-hero"><h1>🧭 Clonal Compass</h1>'
    '<p>Single-cell immune-repertoire co-pilot &mdash; paired scRNA-seq + TCR-seq. '
    'An <strong>exploratory research</strong> tool, not a diagnostic one. '
    f'Answering about: <strong>{html.escape(ds_name)}</strong>.</p></div>'
    '<div class="cc-panel-label">Ask the repertoire</div>'
    '<div class="cc-hint">Try: '
    + " &nbsp;·&nbsp; ".join(_example(q) for q in _EXAMPLES)
    + "</div></div>",
    unsafe_allow_html=True,
)

# Clipboard handler: copies each example's text on click. Runs in the page (via
# the parent document) using execCommand so it works inside the click gesture
# without needing clipboard-write permissions.
components.html(
    """<script>
const d = window.parent.document;
d.querySelectorAll('.cc-copy').forEach(function(el){
  el.onclick = function(){
    const ta = d.createElement('textarea');
    ta.value = el.getAttribute('data-copy');
    ta.style.position='fixed'; ta.style.opacity='0';
    d.body.appendChild(ta); ta.focus(); ta.select();
    try { d.execCommand('copy'); } catch(e){}
    d.body.removeChild(ta);
    const orig = el.innerHTML;
    el.innerHTML = '✓'; el.style.color = '#1f7a4d';
    setTimeout(function(){ el.innerHTML = orig; el.style.color=''; }, 1000);
  };
});
</script>""",
    height=0,
)

question = st.chat_input("Ask about a cluster, clone, or epitope…")

# Existing turns, then the new turn at the bottom (spinner inside the bubble).
for msg in st.session_state.messages:
    with st.chat_message(msg["role"]):
        st.markdown(msg["content"])
        if msg.get("source"):
            st.markdown(_badge_html(msg["source"]), unsafe_allow_html=True)
if question:
    history = list(st.session_state.messages)  # resolves "that cluster"/"it"
    with st.chat_message("user"):
        st.markdown(question)
    with st.chat_message("assistant"):
        with st.spinner("Consulting the evidence bundle…"):
            answer = chat.ask_question(question, bundle, history=history)
        st.markdown(answer.text)
        st.markdown(_badge_html(answer.source), unsafe_allow_html=True)
    st.session_state.messages.append({"role": "user", "content": question})
    st.session_state.messages.append(
        {"role": "assistant", "content": answer.text, "source": answer.source}
    )
