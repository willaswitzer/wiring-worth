"""v1926 add-on: write the complete-axon "set" that like_to_like_trunc.R needs, for a release without truncation runs.

Run from the repo root after scripts/prepare_like_to_like.py (with MICRONS_VERSION set):
    MICRONS_VERSION=1926 .venv/bin/python scripts/make_full_set.py
Writes data/truncation<VSFX>/sets_v<VERSION>/radial_1.00_0.parquet (rows aligned with l2l_pairs: Ld_um, Nsyn as
measured) and index.csv listing it, so the shuffles (cell-type and distance-matched), per-cell-type and
per-cell-type-pair fits can run on complete axons exactly as for the pinned release.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import VERSION, VSFX  # noqa: E402

DATA = ROOT / "data"


def main():
    pairs = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet")
    set_dir = DATA / f"truncation{VSFX}" / f"sets_v{VERSION}"
    set_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"row": np.arange(len(pairs), dtype=np.int32), "Ld_um": pairs["Ld_um"].to_numpy(np.float32),
                  "Nsyn": pairs["Nsyn"].to_numpy(np.int16)}).to_parquet(set_dir / "radial_1.00_0.parquet")
    pd.DataFrame([("radial_1.00_0", "radial", 1.0, 0)], columns=["set", "cond", "level", "rep"]).to_csv(
        set_dir / "index.csv", index=False)
    print(f"v{VERSION}: {len(pairs):,} pairs, {pairs['pre_root_id'].nunique()} complete axons, "
          f"{int(pairs['Nsyn'].sum()):,} synapses -> {set_dir.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
