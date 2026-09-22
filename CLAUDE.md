# CLAUDE.md

## Context

This is a secondary research project using the MICrONS dataset (mouse visual cortex: an EM connectome co-registered with calcium imaging). Access is through the `caveclient` Python package, datastack `minnie65_public`.

- **First result (December):** how the like-to-like wiring rule changes with reconstruction completeness (proofreading level).
- **Later (Feb–Mar):** a connectome-constrained model of activity, compared against shuffled-wiring and data-driven controls.

The dataset changes over time, so never assume table names or version numbers. Query the server for them, and check against https://tutorial.microns-explorer.org when unsure.

## Rules

- Pin one materialization version for all analyses; record it in `config.py`.
- Always report the shuffled-wiring control next to any real-wiring result.
- Never commit data or tokens.
- Keep the core model to co-registered neurons with good connectivity and good functional data.

## Layout

- `notebooks/` — analysis scripts/notebooks, numbered `01_`, `02_`, ... (`# %%` cell format)
- `data/` — downloaded tables (gitignored)
- `figures/` — output figures
- `scripts/` — helpers
- Environment: `.venv/` (Jupyter kernel `microns`); pinned package versions in `requirements.txt`
