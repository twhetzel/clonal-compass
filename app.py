"""Chat interface, step 3: a Streamlit chat UI over the evidence bundle.

Loads `reports/cluster_evidence.json` once on startup and answers each question
via `chat.ask_question`, which grounds every answer in the bundle through
Claude's tool-use loop (or a deterministic fallback when no API key is set).

Run:
  .venv/bin/streamlit run app.py

No deployment needed — this is a local demo tool. It reads only the ~20 KB
evidence bundle, never the .h5ad objects. Regenerate the bundle with
`scripts/generate_report.py`.
"""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from clonal_compass._warnings import silence_demo_warnings

silence_demo_warnings()

from clonal_compass import chat  # noqa: E402

BUNDLE_PATH = Path(__file__).resolve().parent / "reports" / "cluster_evidence.json"

# Source badge styling: green = grounded live Claude answer, amber = fallback.
_BADGE = {
    "Claude API": ("#e6f4ea", "#137333", "Claude API"),
    "deterministic fallback": ("#fef7e0", "#b06000", "deterministic fallback"),
}


@st.cache_data(show_spinner=False)
def _load_bundle(path: str) -> dict:
    """Load the compact evidence bundle once (cached across reruns)."""
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
    "Answers are grounded in `reports/cluster_evidence.json` — a research "
    "tool, not a diagnostic one."
)

if not BUNDLE_PATH.exists():
    st.error(
        f"Evidence bundle not found at `{BUNDLE_PATH}`. "
        "Generate it first with `scripts/generate_report.py`."
    )
    st.stop()

bundle = _load_bundle(str(BUNDLE_PATH))

# Small orientation line so the demo shows what's loaded.
ds = bundle.get("dataset", {})
st.sidebar.header("Dataset")
st.sidebar.write(f"**Cells:** {ds.get('n_cells', '?')}")
st.sidebar.write(f"**Clusters:** {ds.get('n_clusters', '?')}")
st.sidebar.write(f"**TCR+ cells:** {ds.get('n_tcr_cells', '?')}")
st.sidebar.caption(
    "Set `ANTHROPIC_API_KEY` before launching for live Claude answers; "
    "otherwise a deterministic fallback runs."
)

# Conversation history: list of {"role", "content", "source"?}.
if "messages" not in st.session_state:
    st.session_state.messages = []

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
