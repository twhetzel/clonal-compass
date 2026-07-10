"""Chat interface step 2: ask a natural-language question about the dataset.

Reads only reports/cluster_evidence.json (never the .h5ad objects) and answers
via Claude's tool-use loop, grounding every claim in the evidence bundle.

Run:
  .venv/bin/python scripts/ask.py "Which cluster is the most clonally expanded?"
(Set ANTHROPIC_API_KEY to use Claude; otherwise a deterministic fallback runs.)
"""

import sys
from pathlib import Path

from clonal_compass._warnings import silence_demo_warnings

silence_demo_warnings()

from clonal_compass import chat  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BUNDLE = ROOT / "reports" / "cluster_evidence.json"

DEFAULT_Q = "Which cluster is the most clonally expanded, and what is it likely to be?"


def main() -> None:
    question = " ".join(sys.argv[1:]).strip() or DEFAULT_Q
    bundle = chat.load_evidence_bundle(BUNDLE)

    print(f"Q: {question}\n")
    answer = chat.ask_question(question, bundle)
    print(answer.text)
    print(f"\n[source: {answer.source}]")


if __name__ == "__main__":
    main()
