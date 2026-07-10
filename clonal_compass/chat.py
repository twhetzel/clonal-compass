"""Chat interface, step 2: grounded Q&A over the evidence bundle.

`ask_question(question, evidence_bundle)` answers a natural-language question
about the clusters/clones using Claude with **tool use**: instead of dumping the
whole ~20 KB bundle into every prompt, Claude is handed a few small tools
(`get_dataset_overview`, `list_clusters`, `get_cluster_evidence`,
`get_epitope_matches`) and looks up only what a given question needs. Every
answer is therefore grounded in the bundle rather than free-associated.

Guardrails are NOT re-invented here — we import `interpret.SYSTEM_PROMPT`
verbatim and append only operational tool-use instructions on top of it (a chat
turn can span several clusters and vary in length, unlike the single-cluster
interpreter). The scientific hedging rules are unchanged.

Like `interpret.interpret_cluster`, every answer is returned as an
`Interpretation(text, source)` so the caller can tag which path produced it:
`"Claude API"` when a key is set, `"deterministic fallback"` otherwise.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

from .interpret import MODEL, SYSTEM_PROMPT, Interpretation

# Reuse the interpreter's guardrails exactly; add only tool-use operating
# instructions. This does NOT loosen any scientific rule above it.
CHAT_SYSTEM_PROMPT = SYSTEM_PROMPT + textwrap.dedent(
    """

    ---
    You are now answering an immunologist's free-form question about ONE
    processed dataset. You cannot see the data directly — you must call the
    provided tools to look up cluster evidence, and ground every quantitative
    claim in what the tools return. Do not invent cluster ids, markers, counts,
    or epitopes that a tool did not give you.

    - Call `list_clusters` or `get_dataset_overview` first when you need to find
      the relevant cluster(s); then `get_cluster_evidence` /
      `get_epitope_matches` for the specifics.
    - If the tools do not contain what is needed to answer, say so plainly
      rather than guessing.
    - All the hedging rules above still apply: cite the metric behind each
      claim, say "consistent with"/"predicted" rather than "is"/"confirmed",
      and treat a missing epitope match as a normal, valid result.
    - Answer the question actually asked, as briefly as it allows. Plain prose,
      no headings.
    """
)

_MAX_TOOL_ROUNDS = 6  # safety bound on the tool-use loop


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def load_evidence_bundle(path: str | Path = "reports/cluster_evidence.json") -> dict:
    """Load the compact evidence bundle written by `report.write_evidence_json`."""
    return json.loads(Path(path).read_text())


# --------------------------------------------------------------------------- #
# Tools over the evidence bundle
# --------------------------------------------------------------------------- #
# Claude never touches the .h5ad objects; it only queries this in-memory dict.
# Each handler takes (bundle, **input) and returns a small JSON-able object.

TOOLS = [
    {
        "name": "get_dataset_overview",
        "description": (
            "Top-level summary of the dataset: cell/gene/cluster counts, number "
            "of TCR+ cells, and the dataset-wide predicted VDJdb epitope "
            "specificities. Call this to orient before drilling into clusters."
        ),
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "list_clusters",
        "description": (
            "List every cluster with its provisional cell-type label and "
            "expansion stats (n_cells, n_with_tcr, n_expanded_cells, "
            "pct_expanded). Use to find which cluster a question is about. "
            "Set expanded_only=true to restrict to clusters with expanded clones."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "expanded_only": {
                    "type": "boolean",
                    "description": "Only return clusters with >0 expanded cells.",
                }
            },
        },
    },
    {
        "name": "get_cluster_evidence",
        "description": (
            "Full evidence for one cluster by id: provisional cell type, cell "
            "counts, TCR coverage, clonal-expansion stats, top marker genes with "
            "rank_genes_groups scores, epitope hits, and the stored hedged "
            "interpretation if one was generated for this cluster."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "cluster_id": {
                    "type": "string",
                    "description": "The cluster id, e.g. '6'.",
                }
            },
            "required": ["cluster_id"],
        },
    },
    {
        "name": "get_epitope_matches",
        "description": (
            "Predicted VDJdb epitope matches (exact CDR3 aa identity) for one "
            "cluster's expanded clones. Returns an empty list when there is no "
            "match, which is a normal, valid result — not an error."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "cluster_id": {
                    "type": "string",
                    "description": "The cluster id, e.g. '6'.",
                }
            },
            "required": ["cluster_id"],
        },
    },
]


def _find_cluster(bundle: dict, cluster_id: str) -> dict | None:
    cid = str(cluster_id)
    for c in bundle.get("clusters", []):
        if str(c.get("cluster_id")) == cid:
            return c
    return None


def _tool_get_dataset_overview(bundle: dict) -> dict:
    return {
        "dataset": bundle.get("dataset", {}),
        "epitope_overview": bundle.get("epitope_overview", []),
        "n_clusters_in_bundle": len(bundle.get("clusters", [])),
    }


def _tool_list_clusters(bundle: dict, expanded_only: bool = False) -> list[dict]:
    out = []
    for c in bundle.get("clusters", []):
        if expanded_only and not c.get("n_expanded_cells"):
            continue
        out.append(
            {
                "cluster_id": c.get("cluster_id"),
                "cell_type": c.get("cell_type"),
                "n_cells": c.get("n_cells"),
                "n_with_tcr": c.get("n_with_tcr"),
                "n_expanded_cells": c.get("n_expanded_cells"),
                "pct_expanded": c.get("pct_expanded"),
            }
        )
    return out


def _tool_get_cluster_evidence(bundle: dict, cluster_id: str) -> dict:
    c = _find_cluster(bundle, cluster_id)
    if c is None:
        return {
            "error": f"No cluster with id {cluster_id!r}.",
            "available_ids": [x.get("cluster_id") for x in bundle.get("clusters", [])],
        }
    return c


def _tool_get_epitope_matches(bundle: dict, cluster_id: str) -> dict:
    c = _find_cluster(bundle, cluster_id)
    if c is None:
        return {
            "error": f"No cluster with id {cluster_id!r}.",
            "available_ids": [x.get("cluster_id") for x in bundle.get("clusters", [])],
        }
    hits = c.get("epitope_hits", [])
    return {
        "cluster_id": c.get("cluster_id"),
        "epitope_hits": hits,
        "note": (
            "Predicted specificities (VDJdb CDR3 aa hits), not confirmed "
            "reactivity." if hits else
            "No VDJdb match for this cluster's expanded clones — a normal, "
            "valid outcome, not a failure."
        ),
    }


_DISPATCH = {
    "get_dataset_overview": _tool_get_dataset_overview,
    "list_clusters": _tool_list_clusters,
    "get_cluster_evidence": _tool_get_cluster_evidence,
    "get_epitope_matches": _tool_get_epitope_matches,
}


def _run_tool(bundle: dict, name: str, tool_input: dict) -> str:
    handler = _DISPATCH.get(name)
    if handler is None:
        return json.dumps({"error": f"Unknown tool {name!r}."})
    try:
        result = handler(bundle, **(tool_input or {}))
    except Exception as exc:  # noqa: BLE001 - surface, don't crash the loop
        result = {"error": f"{name} failed: {exc}"}
    return json.dumps(result)


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #
def ask_question(question: str, evidence_bundle: dict, use_llm: bool = True) -> Interpretation:
    """Answer a question about the dataset, grounded in the evidence bundle.

    Runs Claude's tool-use loop so the model looks up only the clusters it needs
    (via `TOOLS`) instead of receiving the whole bundle. Returns an
    `Interpretation` whose `.source` records the path taken. Falls back to a
    deterministic, guardrail-compliant response when no API key is available.
    """
    import os

    if use_llm and os.environ.get("ANTHROPIC_API_KEY"):
        try:
            return Interpretation(_ask_with_claude(question, evidence_bundle), "Claude API")
        except Exception as exc:  # noqa: BLE001 - never crash the chat UI
            print(f"[chat] Claude call failed ({exc}); using fallback.")
    return Interpretation(_ask_fallback(question, evidence_bundle), "deterministic fallback")


def _ask_with_claude(question: str, bundle: dict) -> str:
    import anthropic

    client = anthropic.Anthropic()
    messages = [{"role": "user", "content": question}]

    for _ in range(_MAX_TOOL_ROUNDS):
        response = client.messages.create(
            model=MODEL,
            max_tokens=2048,
            thinking={"type": "adaptive"},
            output_config={"effort": "high"},
            system=CHAT_SYSTEM_PROMPT,
            tools=TOOLS,
            messages=messages,
        )
        # Preserve the full assistant turn (thinking + tool_use blocks) verbatim.
        messages.append({"role": "assistant", "content": response.content})

        if response.stop_reason != "tool_use":
            return "".join(b.text for b in response.content if b.type == "text").strip()

        tool_results = []
        for block in response.content:
            if block.type == "tool_use":
                out = _run_tool(bundle, block.name, block.input)
                tool_results.append(
                    {
                        "type": "tool_result",
                        "tool_use_id": block.id,
                        "content": out,
                    }
                )
        messages.append({"role": "user", "content": tool_results})

    # Exhausted the tool-round budget without a final answer.
    return (
        "I could not settle on an answer within the tool-lookup budget for this "
        "question. Try asking about a specific cluster id or a narrower metric."
    )


def _ask_fallback(question: str, bundle: dict) -> str:
    """Deterministic, guardrail-compliant response when no API key is set.

    Without the Claude tool-use loop we cannot parse an arbitrary question, so
    we stay honest: state the limitation and surface the bundle's contents so
    the answer is still useful for orientation.
    """
    ds = bundle.get("dataset", {})
    clusters = bundle.get("clusters", [])
    expanded = [c for c in clusters if c.get("n_expanded_cells")]
    top = sorted(expanded, key=lambda c: -(c.get("n_expanded_cells") or 0))[:3]
    top_desc = (
        "; ".join(
            f"cluster {c.get('cluster_id')} ({c.get('cell_type')}, "
            f"{c.get('n_expanded_cells')} expanded cells)"
            for c in top
        )
        or "none"
    )
    overview = bundle.get("epitope_overview", [])
    epi_desc = "; ".join(overview) if overview else "no dataset-wide VDJdb matches"
    return (
        "Interactive Q&A needs the Claude API (set ANTHROPIC_API_KEY); without "
        "it I can only report the bundle directly rather than answer the "
        f"question '{question.strip()}'. This dataset has "
        f"{ds.get('n_cells', '?')} cells across {ds.get('n_clusters', '?')} "
        f"clusters, with {ds.get('n_tcr_cells', '?')} TCR+ cells. The most "
        f"clonally-expanded clusters are: {top_desc}. Predicted (not confirmed) "
        f"dataset-wide epitope specificities: {epi_desc}. These are exploratory "
        "research observations, not diagnostic or clinical claims."
    )
