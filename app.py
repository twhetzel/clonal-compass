"""Chat interface, step 3: a Streamlit chat UI over the evidence bundle.

Discovers which per-dataset evidence bundles exist under `reports/`
(`cluster_evidence{suffix}.json`, one per entry in `io.DATASETS`), lets the user
pick one in the sidebar, and answers each question via `chat.ask_question` —
which grounds every answer in the selected bundle through Claude's tool-use loop
(or a deterministic fallback when no API key is set). Switching datasets clears
the conversation, since each chat is grounded in a single dataset.

Run:
  .venv/bin/streamlit run app.py

No deployment needed — this is a local demo tool. It reads only the ~20 KB
evidence bundles, never the .h5ad objects. Regenerate them with
`scripts/generate_report.py [--dataset pbmc|cancer]`.
"""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from clonal_compass._warnings import silence_demo_warnings

silence_demo_warnings()

from clonal_compass import chat, io  # noqa: E402

REPORTS_DIR = Path(__file__).resolve().parent / "reports"

# Source badge styling: green = grounded live Claude answer, amber = fallback.
_BADGE = {
    "Claude API": ("#e6f4ea", "#137333", "Claude API"),
    "deterministic fallback": ("#fef7e0", "#b06000", "deterministic fallback"),
}


def _bundle_path(spec: io.DatasetSpec) -> Path:
    """Path to a dataset's evidence bundle (mirrors generate_report.py's suffix)."""
    return REPORTS_DIR / f"cluster_evidence{spec.suffix}.json"


def _available_datasets() -> dict[str, Path]:
    """Registry keys → bundle path, restricted to bundles that exist on disk.

    Not cached: it's a cheap stat() and we want newly-generated bundles to appear
    on the next rerun without clearing any cache.
    """
    return {
        key: _bundle_path(spec)
        for key, spec in io.DATASETS.items()
        if _bundle_path(spec).exists()
    }


@st.cache_data(show_spinner=False)
def _load_bundle(path: str) -> dict:
    """Load one compact evidence bundle (cached per path across reruns)."""
    return chat.load_evidence_bundle(path)


def _badge_html(source: str) -> str:
    bg, fg, label = _BADGE.get(source, ("#eee", "#333", source))
    return (
        f'<span style="background:{bg};color:{fg};padding:2px 8px;'
        f'border-radius:10px;font-size:0.75rem;font-weight:600;'
        f'white-space:nowrap;">{label}</span>'
    )


st.set_page_config(page_title="Clonal Compass — chat", page_icon="🧭")
st.title("🧭 Clonal Compass")
st.caption(
    "Ask natural-language questions about the clusters and clones. "
    "Answers are grounded in the selected dataset's evidence bundle — a research "
    "tool, not a diagnostic one."
)

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

# Each chat is grounded in ONE dataset, so switching invalidates prior turns
# ("that cluster" would point at the other dataset's clusters). Reset on change.
if st.session_state.get("active_dataset") != dataset_key:
    st.session_state.active_dataset = dataset_key
    st.session_state.messages = []

bundle = _load_bundle(str(available[dataset_key]))

# Orientation: name the loaded dataset and its headline stats.
ds = bundle.get("dataset", {})


def _fmt(key: str) -> str:
    v = ds.get(key)
    return f"{v:,}" if isinstance(v, int) else "?"


st.sidebar.caption(ds.get("name", io.DATASETS[dataset_key].display_name))
st.sidebar.write(f"**Cells:** {_fmt('n_cells')}")
st.sidebar.write(f"**Clusters:** {_fmt('n_clusters')}")
st.sidebar.write(f"**TCR+ cells:** {_fmt('n_tcr_cells')}")
st.sidebar.caption(
    "Set `ANTHROPIC_API_KEY` before launching for live Claude answers; "
    "otherwise a deterministic fallback runs."
)

# Make the active dataset unmistakable in the main pane too.
st.info(f"Answering about: **{ds.get('name', io.DATASETS[dataset_key].display_name)}**", icon=":material/database:")

# Conversation history: list of {"role", "content", "source"?}.
st.session_state.setdefault("messages", [])

# Replay the scrollable conversation history.
for msg in st.session_state.messages:
    with st.chat_message(msg["role"]):
        st.markdown(msg["content"])
        if msg.get("source"):
            st.markdown(_badge_html(msg["source"]), unsafe_allow_html=True)

# Chat input at the bottom.
if question := st.chat_input("Ask about a cluster, clone, or epitope…"):
    # Prior turns (everything before this new question) become the context that
    # lets Claude resolve follow-ups like "that cluster" or "it".
    history = list(st.session_state.messages)

    st.session_state.messages.append({"role": "user", "content": question})
    with st.chat_message("user"):
        st.markdown(question)

    with st.chat_message("assistant"):
        with st.spinner("Consulting the evidence bundle…"):
            answer = chat.ask_question(question, bundle, history=history)
        st.markdown(answer.text)
        st.markdown(_badge_html(answer.source), unsafe_allow_html=True)

    st.session_state.messages.append(
        {"role": "assistant", "content": answer.text, "source": answer.source}
    )
