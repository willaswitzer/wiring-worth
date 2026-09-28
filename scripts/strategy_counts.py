"""Session 4, task 1: proofreading strategy columns for the usable set, by axon level.

Run from the repo root: `.venv/bin/python scripts/strategy_counts.py`
proofreading_status_and_strategy has two strategy columns (strategy_axon, strategy_dendrite)
and two boolean status columns (status_axon, status_dendrite). axon_level in the usable set
is strategy_axon, so its own crosstab is the diagonal.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import VERSION  # noqa: E402

DATA = ROOT / "data"
AXON_ORDER = ["axon_fully_extended", "axon_partially_extended", "axon_interareal", "none"]

proof = pd.read_parquet(DATA / f"proofreading_status_and_strategy_v{VERSION}.parquet")
usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
m = usable[["root_id", "axon_level"]].merge(proof, left_on="root_id", right_on="pt_root_id", validate="1:1")
assert len(m) == len(usable) and (m["axon_level"] == m["strategy_axon"]).all()

cols = [c for c in proof.columns if c.startswith(("strategy_", "status_"))]
for c in cols:
    t = pd.crosstab(m[c].astype(str), m["axon_level"], margins=True, margins_name="All")
    t = t.reindex(columns=[*AXON_ORDER, "All"], fill_value=0)
    print(f"\n### {c}\n\n{t.to_markdown()}")
