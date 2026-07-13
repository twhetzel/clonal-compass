# Clonal Compass

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)

Clonal Compass is a single-cell immune-repertoire co-pilot for 
paired scRNA-seq + TCR-seq data. It turns standard repertoire-analysis outputs 
into a reviewable report with a grounded chat interface that helps immunologists 
ask: which T-cell clones expanded, what transcriptional states are they 
associated with, and what evidence supports that interpretation?

This project was created for the [Build with Claude: Life Sciences](https://cerebralvalley.ai/e/built-with-claude-life-sciences) hackathon.

## Requirements

- **Python 3.11+** (developed and pinned against 3.11.9).

> ⚠️ On many systems `python3` / `python` points to an older interpreter
> (e.g. macOS may ship 3.8, which is too old for this project's scanpy/scirpy
> stack). Always create the venv with an explicit `python3.11` so you get the
> right version. Check with `python3.11 --version` before starting.

## Setup

```bash
# 1. Create a dedicated virtual environment with Python 3.11
python3.11 -m venv .venv

# 2. Activate it
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 3. Install the pinned dependencies
pip install -r requirements.txt

# 4. Install clonal_compass as an editable package (so `import clonal_compass`
#    works from anywhere — no PYTHONPATH needed). Dependencies are read from
#    requirements.txt, so this step is safe to run right after it.
pip install -e .
```

If `python3.11` isn't on your PATH, point at it explicitly, e.g. via pyenv:
`~/.pyenv/versions/3.11.9/bin/python3.11 -m venv .venv`.

**Apple Silicon note:** `scirpy` depends on `parasail`, which has no arm64
macOS wheel and builds from source. If the install fails on `parasail`, run
`brew install automake libtool` once, then retry.

## Sanity check

With the venv active, load the demo dataset and confirm it reads into AnnData:

```bash
python scripts/load_data.py
```


## How Claude was used in this project
This project used Claude across three different surfaces, each for a distinct purpose:

<b>Claude Chat</b> — used for project scoping, architecture planning, and risk assessment before writing any code: comparing hackathon project ideas, breaking the build into a day-by-day plan with explicit checkpoints and fallback triggers, identifying and de-risking the biggest technical unknowns (environment/dependency setup, scientific overclaiming) ahead of time, and reviewing/interpreting results at each major milestone (UMAP sanity checks, clonotype definition tradeoffs, live-API guardrail verification) to catch issues early rather than after the fact.

<b>Claude Code</b> — used to build the entire pipeline and application: environment setup and dependency debugging (including diagnosing a Python version mismatch and fixing two real upstream bugs — a broken scirpy VDJdb loader path, and an epitope-matching tool that generalized a per-cluster "no match" into an incorrect dataset-wide "no match"), the Scanpy/scirpy analysis pipeline (QC, clustering, T-cell subset annotation, clonal expansion metrics), the report generation layer, and the Streamlit chat interface — including diagnosing a dataset-specific marker annotation failure and building a corrected, per-dataset marker registry as the fix.

<b>Claude API</b> (via the Anthropic Python SDK, claude-opus-4-8 with adaptive thinking) — the core scientific differentiator of the tool itself, called directly from the pipeline code (not just used to build it) to generate plain-language, guardrail-compliant interpretations of T-cell clusters and clonal expansion patterns, and to power a tool-using conversational agent that looks up only the specific evidence needed to answer a given question (rather than receiving the full dataset every time), grounded against VDJdb epitope data. A deterministic, guardrail-compliant fallback path ensures the pipeline runs end-to-end even without an API key.

## Why this is designed to be cautious
Clonal Compass is built around a specific concern: LLMs summarizing biological results can sound confident even when the underlying evidence is weak or ambiguous. Several things are built in to guard against that:

- Every interpretation is grounded in a specific cited metric — no claim is made without a number behind it (expansion percentage, marker gene score, TCR coverage).
- Hedged language is enforced, not just encouraged — the system prompt explicitly requires "consistent with" rather than "is" or "proves," and this was stress-tested against a deliberately leading question ("does this data suggest infection?") designed to provoke overclaiming; the tool held its hedged framing under that pressure.
- "No match" is treated as a normal, valid result, not a failure — epitope database misses are common (VDJdb's coverage is limited) and are explicitly framed that way, rather than as a gap in the analysis.
- No clinical or patient-level claims, ever — every interpretation explicitly frames itself as exploratory research, and the tool refuses direct requests for treatment/clinical recommendations.
- A deterministic, non-LLM fallback path exists — if the API is unavailable, the tool still produces guardrail-compliant output rather than failing or guessing.

## Known limitations
- Epitope matching is exact-CDR3 only — misses functionally similar TCRs with slightly different sequences, and VDJdb's coverage is sparse and skewed toward well-studied viral epitopes. Absence of a match often reflects database gaps, not biology.
- Cell-type annotation is marker-gene-based and unvalidated against ground truth — reasonably reliable for the two dataset types tested here, but not guaranteed to generalize.
- No HLA typing — epitope matches aren't filtered for whether they're even biologically plausible given a sample's genetic background.
- Research tool only — not validated for, and not intended for, any clinical or diagnostic use.
- Tested on two public datasets only (10x PBMC, Wu et al. 2020 tumor-infiltrating T cells); behavior on other data is unverified.
- Live interpretation requires an Anthropic API key; without one, the tool uses a deterministic fallback that's guardrail-compliant but less nuanced.
- No automated test suite — verification has been manual (data sanity checks, targeted question testing), not unit/integration tests.
- LLM outputs should be spot-checked, not trusted blindly — guardrails were tested extensively during development, but Claude's interpretations aren't immune to occasional inconsistency, especially on data outside our testing.

# Demo Materials
- [Demo slides](docs/slides/clonal-compass-demo.pdf)
- [Demo video]()

## Sample Reports

Reviewers can inspect example outputs without running the pipeline:

- [PBMC baseline report](https://twhetzel.github.io/clonal-compass/reports/clonal_compass_report.html)
- [Cancer TIL report](https://twhetzel.github.io/clonal-compass/reports/clonal_compass_report_cancer.html)

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE.md) for details.
