"""Day 3: Claude interpretation layer.

Turns a cluster's marker genes + clonal expansion stats into a plain-language,
appropriately-hedged interpretation an immunologist could actually use.

The scientific guardrails (see CLAUDE.md) are enforced in the system prompt AND
mirrored in a deterministic fallback so the pipeline still runs without an API
key. Guardrails:
  - distinguish "consistent with X" from "is X" / "proves X"
  - cite the specific metric behind every claim
  - never state a conclusion more strongly than the stats support
  - prefer more, shorter hedged claims over fewer confident ones
  - research/exploratory tool, not diagnostic; no patient/clinical claims
"""

from __future__ import annotations

import os
import textwrap
from dataclasses import dataclass, field

MODEL = "claude-opus-4-8"

SYSTEM_PROMPT = textwrap.dedent(
    """\
    You are a cautious computational-immunology assistant. You are given
    quantitative evidence about a single cell cluster from a paired
    scRNA-seq + TCR-seq experiment and must write a short, plain-language
    interpretation a working immunologist could use.

    Follow these rules without exception:
    - Distinguish "consistent with X" or "suggestive of X" from "is X" or
      "proves X". Never state a conclusion more strongly than the numbers
      support.
    - Cite the specific metric behind every claim (a marker gene and its
      score, an expansion fraction, a clone size, a percentage). No claim
      without a number attached.
    - Prefer several short, individually-hedged statements over one
      confident-sounding paragraph.
    - Epitope/VDJdb matches are database hits, not proof of function: call
      them "predicted" or "possible", never "confirmed". If no match is
      given, say so plainly — absence of a match is a normal, valid result,
      not a failure.
    - This is a research/exploratory analysis, not a diagnostic one. Make no
      patient-level or clinical claims. State this framing if you draw any
      functional conclusion.
    - If the evidence is ambiguous, say it is ambiguous. Do not invent
      markers, numbers, or cell identities beyond what you are given.

    Write 4-8 sentences. No preamble, no headings — just the interpretation.
    """
)


@dataclass
class Interpretation:
    """An interpretation plus which path produced it (for unambiguous tagging)."""

    text: str
    source: str  # "Claude API" or "deterministic fallback"


@dataclass
class ClusterEvidence:
    """Quantitative evidence for one cluster, fed to the interpreter."""

    cluster_id: str
    cell_type: str
    n_cells: int
    top_markers: list[tuple[str, float]]  # (gene, score) high -> low
    n_with_tcr: int
    pct_with_tcr: float
    n_expanded_cells: int          # cells whose clone has size > 1
    pct_expanded: float            # of TCR+ cells in this cluster
    top_clone_sizes: list[int]     # largest clone sizes in this cluster
    epitope_hits: list[str] = field(default_factory=list)  # "species: epitope (n cells)"

    def as_prompt(self) -> str:
        markers = ", ".join(f"{g} ({s:+.2f})" for g, s in self.top_markers)
        clones = ", ".join(str(s) for s in self.top_clone_sizes) or "none > 1"
        epitopes = (
            "\n".join(f"    - {h}" for h in self.epitope_hits)
            if self.epitope_hits
            else "    - none (no VDJdb match for this cluster's expanded clones)"
        )
        return textwrap.dedent(
            f"""\
            Cluster: {self.cluster_id}
            Marker-based label (provisional): {self.cell_type}
            Cells in cluster: {self.n_cells}
            Top marker-gene scores (sc.tl.rank_genes_groups, higher = more specific):
                {markers}
            TCR coverage: {self.n_with_tcr}/{self.n_cells} cells have a TCR ({self.pct_with_tcr:.0f}%)
            Clonal expansion: {self.n_expanded_cells} cells in clones of size > 1
                ({self.pct_expanded:.0f}% of TCR+ cells in this cluster)
            Largest clone sizes here: {clones}
            Predicted epitope matches (VDJdb, exact CDR3 aa identity):
            {epitopes}
            """
        )


def interpret_cluster(evidence: ClusterEvidence, use_llm: bool = True) -> Interpretation:
    """Return a hedged plain-language interpretation plus its source path.

    Uses the Claude API when an API key is available; otherwise returns a
    deterministic, guardrail-compliant fallback so the pipeline still runs.
    The `.source` field records which path produced the text.
    """
    if use_llm and os.environ.get("ANTHROPIC_API_KEY"):
        try:
            return Interpretation(_interpret_with_claude(evidence), "Claude API")
        except Exception as exc:  # noqa: BLE001 - never break the pipeline
            print(f"[interpret] Claude call failed ({exc}); using fallback.")
    return Interpretation(_interpret_fallback(evidence), "deterministic fallback")


def _interpret_with_claude(evidence: ClusterEvidence) -> str:
    import anthropic

    client = anthropic.Anthropic()
    response = client.messages.create(
        model=MODEL,
        max_tokens=1024,
        thinking={"type": "adaptive"},
        output_config={"effort": "high"},
        system=SYSTEM_PROMPT,
        messages=[{"role": "user", "content": evidence.as_prompt()}],
    )
    return "".join(b.text for b in response.content if b.type == "text").strip()


def _interpret_fallback(evidence: ClusterEvidence) -> str:
    """Deterministic, guardrail-compliant summary (no LLM).

    Every sentence cites a metric and hedges; used when no API key is set.
    """
    top = ", ".join(f"{g} ({s:+.2f})" for g, s in evidence.top_markers[:5])
    lines = [
        f"Cluster {evidence.cluster_id} ({evidence.n_cells} cells) shows highest "
        f"marker-gene scores for {top}, which is consistent with the "
        f"provisional label '{evidence.cell_type}' rather than proving it.",
    ]
    if evidence.pct_expanded >= 10:
        lines.append(
            f"{evidence.pct_expanded:.0f}% of the {evidence.n_with_tcr} TCR+ cells "
            f"here sit in expanded clones (size > 1; largest "
            f"{max(evidence.top_clone_sizes, default=1)}), a level of clonal "
            f"expansion consistent with an antigen-experienced population."
        )
    else:
        lines.append(
            f"Only {evidence.pct_expanded:.0f}% of the {evidence.n_with_tcr} TCR+ "
            f"cells are in expanded clones, so there is little evidence of clonal "
            f"expansion in this cluster."
        )
    if evidence.epitope_hits:
        lines.append(
            "VDJdb reports predicted (not confirmed) matches for some expanded "
            f"clones: {'; '.join(evidence.epitope_hits)}. These are database "
            "hits by CDR3 amino-acid identity and require functional validation."
        )
    else:
        lines.append(
            "No VDJdb epitope match was found for this cluster's expanded clones, "
            "which is a normal and valid outcome, not a failure."
        )
    lines.append(
        "This is an exploratory research analysis, not a diagnostic one; no "
        "patient-level or clinical conclusions should be drawn."
    )
    return " ".join(lines)
