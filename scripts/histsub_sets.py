"""Session 12: count-matched random loss -- does losing each axon's synapses at random reproduce the history?

Run from the repo root after scripts/edit_history.py build:
    .venv/bin/python scripts/histsub_sets.py [--reps 20] [--seed 12]
then fit with
    TRUNC_INDEX=index_histsub.csv TRUNC_TAG=histsub .r-env/bin/Rscript scripts/like_to_like_trunc.R fit poisson 4

For each history state (Nov 2025, Feb 2026) and repeat, every complete axon keeps a random subset of its l2l
synapses of the same size as it had attached at that state (Ld left complete). Within an axon the loss is random;
across axons it matches the real history exactly. If the real history's weaker orientation rule comes from which
axons were incomplete (not which synapses within them), these sets reproduce it.
Writes data/truncation/sets_v1822/histsub_<state>_<rep>.parquet and index_histsub.csv (condition histsub,
level = state, rep).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import VERSION  # noqa: E402

DATA = ROOT / "data"
SET_DIR = DATA / "truncation" / f"sets_v{VERSION}"
STATES = [202511, 202602]
MAX_DIST_UM = 3.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--seed", type=int, default=12)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    usable["root_id"] = usable["root_id"].astype("int64")
    pairs = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        pairs[c] = pairs[c].astype("int64")
    pres = pairs["pre_root_id"].unique()
    syn = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        syn[c] = syn[c].astype("int64")
    syn = syn[(syn["pre_root_id"] != syn["post_root_id"]) & (syn["pre_axon_dist_um"] <= MAX_DIST_UM)
              & syn["pre_root_id"].isin(pres) & syn["post_root_id"].isin(usable["root_id"])].reset_index(drop=True)
    pair_idx = pd.Series(np.arange(len(pairs)), index=pd.MultiIndex.from_frame(pairs[["pre_root_id", "post_root_id"]]))
    s_row = pair_idx.reindex(pd.MultiIndex.from_frame(syn[["pre_root_id", "post_root_id"]])).to_numpy()
    assert not np.isnan(s_row).any()
    s_row = s_row.astype(int)
    full = np.bincount(s_row, minlength=len(pairs))
    assert np.array_equal(full, pairs["Nsyn"].to_numpy()), "l2l synapses do not reproduce l2l_pairs Nsyn"

    groups = syn.groupby("pre_root_id").indices
    index = [("radial_1.00_0", "radial", 1.0, 0)]
    summary = []
    for v in STATES:
        att = pd.read_parquet(DATA / "history" / f"attached_{v}.parquet").set_index("synapse_id")["attached"]
        a = att.reindex(syn["synapse_id"]).fillna(False).to_numpy(bool)
        hist_ns = np.bincount(s_row[a], minlength=len(pairs))
        ref = pd.read_parquet(SET_DIR / f"hist_{v}_0.parquet")
        assert np.array_equal(hist_ns, ref["Nsyn"].to_numpy()), f"attached synapses do not reproduce hist_{v}_0 Nsyn"
        k = {r: int(a[ix].sum()) for r, ix in groups.items()}
        for rep in range(args.reps):
            keep = np.zeros(len(syn), bool)
            for r, ix in groups.items():
                keep[rng.choice(ix, k[r], replace=False)] = True
            ns = np.bincount(s_row[keep], minlength=len(pairs))
            name = f"histsub_{v}_{rep}"
            pd.DataFrame({"row": np.arange(len(pairs), dtype=np.int32), "Ld_um": pairs["Ld_um"].to_numpy(np.float32),
                          "Nsyn": ns.astype(np.int16)}).to_parquet(SET_DIR / f"{name}.parquet")
            index.append((name, "histsub", float(v), rep))
        summary.append((v, int(a.sum()), len(a)))
        print(f"{v}: {a.sum()} of {len(a)} synapses kept per state; {args.reps} reps written", flush=True)
    pd.DataFrame(index, columns=["set", "cond", "level", "rep"]).to_csv(SET_DIR / "index_histsub.csv", index=False)


if __name__ == "__main__":
    main()
