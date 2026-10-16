"""Session 12: roll back the postsynaptic side too.

Session 8 rolled back only the presynaptic axon: a synapse counted in a past state if its presynaptic supervoxel
was attached to the presynaptic soma then. Here it must also have been attached to the postsynaptic neuron's soma
(its postsynaptic supervoxel had the same root as that neuron's nucleus supervoxel at that time).

Stage 1 -- fetch (needs CAVE; ~5 min):
    .venv/bin/python scripts/postsyn_history.py fetch
  postsynaptic supervoxel of every l2l synapse (synapse table @ v1822), nucleus supervoxel of every usable
  neuron (nucleus_detection_v0), and their roots at Apr 2023, Nov 2025, Feb 2026 and v1822.
  Writes data/history/attached_post_<state>.parquet.
Stage 2 -- build (local):
    .venv/bin/python scripts/postsyn_history.py build
  Sets histpp_<state>_0: Nsyn counts synapses attached on both sides; Ld is the presynaptic-rollback Ld of
  hist_<state>_0 (dendrite cable is not rolled back). Writes index_histpp.csv and
  notebooks/12_postsyn_history.md (synapses attached on each side, per state).
Fit with TRUNC_INDEX=index_histpp.csv TRUNC_TAG=histpp (see scripts/run_session12.sh).
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import DATASTACK, VERSION  # noqa: E402

DATA = ROOT / "data"
HIST = DATA / "history"
SET_DIR = DATA / "truncation" / f"sets_v{VERSION}"
STATES = [202304, 202511, 202602, VERSION]
CHUNK = 500
CHUNK_SV = 20000
MAX_DIST_UM = 3.0


def l2l_synapses():
    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    usable["root_id"] = usable["root_id"].astype("int64")
    pairs = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        pairs[c] = pairs[c].astype("int64")
    syn = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        syn[c] = syn[c].astype("int64")
    syn = syn[(syn["pre_root_id"] != syn["post_root_id"]) & (syn["pre_axon_dist_um"] <= MAX_DIST_UM)
              & syn["pre_root_id"].isin(pairs["pre_root_id"].unique())
              & syn["post_root_id"].isin(usable["root_id"])].reset_index(drop=True)
    return usable, pairs, syn


def fetch(args):
    from caveclient import CAVEclient
    client = CAVEclient(DATASTACK, version=VERSION)
    syn_table = client.info.get_datastack_info()["synapse_table"]
    nuc_table = client.info.get_datastack_info().get("soma_table") or "nucleus_detection_v0"
    usable, _, syn = l2l_synapses()
    t0 = time.time()
    sv_path = HIST / "post_sv.parquet"
    if not sv_path.exists():
        ids = syn["synapse_id"].astype("int64").tolist()
        parts = []
        for i in range(0, len(ids), CHUNK):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                df = client.materialize.query_table(syn_table, filter_in_dict={"id": ids[i:i + CHUNK]},
                                                    select_columns=["id", "post_pt_root_id", "post_pt_supervoxel_id"])
            assert not any("Limited query" in str(x.message) for x in w), f"chunk {i} truncated"
            parts.append(df)
        sv = pd.concat(parts, ignore_index=True).rename(columns={"id": "synapse_id", "post_pt_root_id": "post_root_id",
                                                                 "post_pt_supervoxel_id": "sv"})
        assert len(sv) == len(ids), f"{len(sv)} of {len(ids)} synapses returned"
        sv.to_parquet(sv_path)
    sv = pd.read_parquet(sv_path)
    nuc_path = HIST / "post_nucleus_sv.parquet"
    if not nuc_path.exists():
        nuc = client.materialize.query_table(nuc_table, filter_in_dict={"id": usable["nucleus_id"].astype("int64").tolist()},
                                             select_columns=["id", "pt_supervoxel_id", "pt_root_id"])
        nuc = nuc.rename(columns={"id": "nucleus_id", "pt_supervoxel_id": "nuc_sv", "pt_root_id": "root_now"})
        nuc.to_parquet(nuc_path)
    nuc = pd.read_parquet(nuc_path)
    u = usable[["root_id", "nucleus_id"]].astype("int64").merge(nuc, on="nucleus_id", validate="1:1")
    bad = (u["root_now"].astype("int64") != u["root_id"]).sum()
    if bad:
        print(f"WARNING: {bad} postsynaptic nuclei have a different root at v{VERSION} than the usable set")
    nuc_sv = u.set_index("root_id")["nuc_sv"]
    meta = pd.read_csv(HIST / "versions.csv")
    meta["time_stamp"] = pd.to_datetime(meta["time_stamp"], format="ISO8601", utc=True)
    post_of = sv["post_root_id"].astype("int64").to_numpy()
    sv_ids = sv["sv"].to_numpy(np.uint64)
    for v in STATES:
        out = HIST / f"attached_post_{v}.parquet"
        if out.exists():
            continue
        ts = meta.loc[meta["version"] == v, "time_stamp"].iloc[0].to_pydatetime()
        r_syn = np.concatenate([client.chunkedgraph.get_roots(sv_ids[i:i + CHUNK_SV], timestamp=ts)
                                for i in range(0, len(sv_ids), CHUNK_SV)])
        r_nuc = pd.Series(client.chunkedgraph.get_roots(nuc_sv.to_numpy(np.uint64), timestamp=ts), index=nuc_sv.index)
        attached = r_syn == r_nuc.reindex(post_of).to_numpy(np.uint64)
        pd.DataFrame({"synapse_id": sv["synapse_id"].to_numpy(), "post_root_id": post_of,
                      "attached": attached}).to_parquet(out)
        print(f"  {v}: {attached.mean():.3f} of l2l synapses attached to their postsynaptic soma; "
              f"{time.time() - t0:.0f}s", flush=True)


def build(args):
    usable, pairs, syn = l2l_synapses()
    pair_idx = pd.Series(np.arange(len(pairs)), index=pd.MultiIndex.from_frame(pairs[["pre_root_id", "post_root_id"]]))
    s_row = pair_idx.reindex(pd.MultiIndex.from_frame(syn[["pre_root_id", "post_root_id"]])).to_numpy().astype(int)
    n = pd.read_parquet(DATA / f"l2l_neurons_v{VERSION}.parquet")
    n["root_id"] = n["root_id"].astype("int64")
    area = n.set_index("root_id")["area"]
    v1v1 = ((area.reindex(syn["pre_root_id"]).to_numpy() == "V1") & (area.reindex(syn["post_root_id"]).to_numpy() == "V1"))
    index, rows = [("radial_1.00_0", "radial", 1.0, 0)], []
    for v in STATES:
        pre_att = pd.read_parquet(HIST / f"attached_{v}.parquet").set_index("synapse_id")["attached"]
        post_att = pd.read_parquet(HIST / f"attached_post_{v}.parquet").set_index("synapse_id")["attached"]
        a_pre = pre_att.reindex(syn["synapse_id"]).fillna(False).to_numpy(bool)
        a_post = post_att.reindex(syn["synapse_id"]).fillna(False).to_numpy(bool)
        both = a_pre & a_post
        hist = pd.read_parquet(SET_DIR / f"hist_{v}_0.parquet")
        assert np.array_equal(np.bincount(s_row[a_pre], minlength=len(pairs)), hist["Nsyn"].to_numpy())
        name = f"histpp_{v}_0"
        pd.DataFrame({"row": hist["row"].to_numpy(), "Ld_um": hist["Ld_um"].to_numpy(),
                      "Nsyn": np.bincount(s_row[both], minlength=len(pairs)).astype(np.int16)}).to_parquet(SET_DIR / f"{name}.parquet")
        index.append((name, "histpp", float(v), 0))
        rows.append(dict(state=v, synapses=len(syn), presynaptic_attached=int(a_pre.sum()),
                         postsynaptic_attached=int(a_post.sum()), both_attached=int(both.sum()),
                         V1_V1_pre=int(a_pre[v1v1].sum()), V1_V1_both=int(both[v1v1].sum()),
                         V1_V1_total=int(v1v1.sum())))
        print(rows[-1], flush=True)
    pd.DataFrame(index, columns=["set", "cond", "level", "rep"]).to_csv(SET_DIR / "index_histpp.csv", index=False)
    (ROOT / "notebooks" / "12_postsyn_history.md").write_text(
        f"# Session 12: synapses attached on each side, per history state (auto-generated, v{VERSION})\n\n"
        "A synapse counts in `histpp` sets only if attached to both the presynaptic and the postsynaptic soma.\n\n"
        + pd.DataFrame(rows).to_markdown(index=False) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["fetch", "build"])
    args = ap.parse_args()
    {"fetch": fetch, "build": build}[args.stage](args)


if __name__ == "__main__":
    main()
