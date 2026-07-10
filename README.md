# clonal-compass

Build with Claude: Life Sciences hackathon project — a single-cell immune
repertoire co-pilot for paired scRNA-seq + TCR-seq (10x 5' GEX + VDJ) data.

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

See `CLAUDE.md` for the full project brief, environment details, and the
4-day plan.
