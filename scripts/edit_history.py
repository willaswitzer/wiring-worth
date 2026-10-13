"""Session 8: the real proofreading history of the fully extended axons (edit-log natural experiment).

Suggested by C. Schneider-Mizell (Allen Institute, 2026-10-13): partial axon extension in MICrONS
followed axons to a fixed radial distance around the proofread column, so the bias of real
incompleteness is best measured on the same cells' own past states, not simulated.

Stage 1 -- fetch (needs CAVE access; your Mac):
    .venv/bin/python scripts/edit_history.py fetch [--limit N]
  - presynaptic supervoxel of every output synapse of the 193 fully extended axons (synapses_pni_2 @ v1822)
  - nucleus supervoxel of each of those neurons (nucleus_detection_v0 @ v1822)
  - for the timestamp of every materialization version of the datastack (expired ones included; from
    get_versions_metadata), the root id of each synapse supervoxel and each nucleus supervoxel
    (chunkedgraph.get_roots). An output synapse is "attached at t" if its supervoxel's root at t equals
    the nucleus supervoxel's root at t.
  - for every version that can still be queried, the axon proofreading label of each neuron at that
    version (proofreading table joined through the nucleus id).
  Resumable: one parquet per timestamp in data/history/.

Stage 2 -- build (local, no CAVE):
    .venv/bin/python scripts/edit_history.py build
  For each timestamp, every axon is cut down to its state at t: an axon piece (1 um, as in
  compute_proximity.py) is kept if its nearest output synapse on that axon was attached at t;
  Ld is recomputed from the kept pieces and Nsyn counts the kept l2l synapses. Sets are written in the
  truncation format as condition `hist` (one set per version, level = version number), with
  index_hist.csv; fit them with
    TRUNC_INDEX=index_hist.csv TRUNC_TAG=hist .r-env/bin/Rscript scripts/like_to_like_trunc.R fit poisson 4
    TRUNC_INDEX=index_hist.csv TRUNC_TAG=hist .r-env/bin/Rscript scripts/like_to_like_trunc.R boot poisson 200 4
  Writes data/history/history_axons_v1822.parquet (per axon x version: synapses / cable attached, label)
  and notebooks/08_history_summary.md.
Caveats: synapse attribution on the postsynaptic side is today's (only the presynaptic axon is rolled
back); a piece's state is inferred from its nearest output synapse, so branches without synapses
follow their neighbours.
"""
import argparse
import itertools
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from config import DATASTACK, VERSION  # noqa: E402

DATA = ROOT / "data"
HIST = DATA / "history"
SET_DIR = DATA / "truncation" / f"sets_v{VERSION}"
CHUNK_ROOTS = 10
CHUNK_SV = 20000
PROOF_TABLES = ["proofreading_status_and_strategy", "proofreading_status_public_release"]


def ref_pres():
    pairs = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet", columns=["pre_root_id"])
    pres = pairs["pre_root_id"].astype("int64").unique()
    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    usable["root_id"] = usable["root_id"].astype("int64")
    u = usable.set_index("root_id").loc[pres]
    return pd.DataFrame({"root_id": pres, "nucleus_id": u["nucleus_id"].astype("int64").to_numpy(),
                         "cell_type": u["cell_type"].to_numpy()})


# ------------------------------------------------------------------ stage 1
def fetch(args):
    global HIST
    from caveclient import CAVEclient
    if args.limit or args.versions:
        HIST = DATA / "history_test"  # a test run never mixes with the real outputs
    HIST.mkdir(parents=True, exist_ok=True)
    client = CAVEclient(DATASTACK, version=VERSION)
    syn_table = client.info.get_datastack_info()["synapse_table"]
    nuc_table = client.info.get_datastack_info().get("soma_table") or "nucleus_detection_v0"
    pres = ref_pres()
    if args.limit:
        pres = pres.head(args.limit)
    t0 = time.time()

    # supervoxels of output synapses and nuclei, at the pinned version
    sv_path = HIST / "output_sv.parquet"
    if not sv_path.exists():
        parts = []
        roots = pres["root_id"].tolist()
        for i in range(0, len(roots), CHUNK_ROOTS):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                df = client.materialize.query_table(
                    syn_table, filter_in_dict={"pre_pt_root_id": roots[i:i + CHUNK_ROOTS]},
                    select_columns=["id", "pre_pt_root_id", "pre_pt_supervoxel_id"])
            assert not any("Limited query" in str(x.message) for x in w), f"chunk {i} truncated"
            parts.append(df)
            print(f"  output synapses: {i + CHUNK_ROOTS}/{len(roots)} roots, {time.time() - t0:.0f}s", flush=True)
        sv = pd.concat(parts, ignore_index=True).rename(columns={"id": "synapse_id", "pre_pt_root_id": "pre_root_id",
                                                                 "pre_pt_supervoxel_id": "sv"})
        sv.to_parquet(sv_path)
    sv = pd.read_parquet(sv_path)
    nuc_path = HIST / "nucleus_sv.parquet"
    if not nuc_path.exists():
        nuc = client.materialize.query_table(nuc_table, filter_in_dict={"id": pres["nucleus_id"].tolist()},
                                             select_columns=["id", "pt_supervoxel_id", "pt_root_id"])
        nuc = nuc.rename(columns={"id": "nucleus_id", "pt_supervoxel_id": "nuc_sv", "pt_root_id": "root_now"})
        nuc.to_parquet(nuc_path)
    nuc = pd.read_parquet(nuc_path)
    pres = pres.merge(nuc, on="nucleus_id", validate="1:1")
    bad = (pres["root_now"].astype("int64") != pres["root_id"]).sum()
    if bad:
        print(f"WARNING: {bad} nuclei have a different root at v{VERSION} than the usable set; check before trusting them")
    print(f"{len(sv):,} output synapses of {len(pres)} axons; {time.time() - t0:.0f}s", flush=True)

    # timestamps of every version, expired ones included
    meta = pd.DataFrame(client.materialize.get_versions_metadata(expired=True))
    meta = meta[["version", "time_stamp"]].drop_duplicates("version").sort_values("version")
    meta["time_stamp"] = pd.to_datetime(meta["time_stamp"], utc=True)
    meta = meta[meta["version"] <= VERSION]
    if args.versions:
        meta = meta[meta["version"].isin([int(v) for v in args.versions.split(",")])]
    meta.to_csv(HIST / "versions.csv", index=False)
    print(f"{len(meta)} versions: {meta['version'].tolist()}", flush=True)

    # roots at each timestamp
    sv_ids = sv["sv"].to_numpy(np.uint64)
    pre_of = sv["pre_root_id"].to_numpy(np.int64)
    nuc_sv = pres.set_index("root_id")["nuc_sv"]
    for _, m in meta.iterrows():
        out = HIST / f"attached_v{int(m.version):04d}.parquet"
        if out.exists():
            continue
        ts = m.time_stamp.to_pydatetime()
        r_syn = np.concatenate([client.chunkedgraph.get_roots(sv_ids[i:i + CHUNK_SV], timestamp=ts)
                                for i in range(0, len(sv_ids), CHUNK_SV)])
        r_nuc = pd.Series(client.chunkedgraph.get_roots(nuc_sv.to_numpy(np.uint64), timestamp=ts), index=nuc_sv.index)
        attached = r_syn == r_nuc.reindex(pre_of).to_numpy(np.uint64)
        pd.DataFrame({"synapse_id": sv["synapse_id"].to_numpy(), "pre_root_id": pre_of,
                      "attached": attached}).to_parquet(out)
        frac = pd.Series(attached).groupby(pre_of).mean()
        print(f"  v{int(m.version)} ({ts:%Y-%m-%d}): median {frac.median():.2f} of each axon's outputs attached; "
              f"{time.time() - t0:.0f}s", flush=True)

    # proofreading labels at each still-queryable version
    lab_path = HIST / "labels_by_version.parquet"
    labs = []
    for v in meta["version"]:
        try:
            n_v = client.materialize.query_table(nuc_table, filter_in_dict={"id": pres["nucleus_id"].tolist()},
                                                 select_columns=["id", "pt_root_id"], materialization_version=int(v))
        except Exception as e:  # noqa: BLE001  (expired or missing)
            print(f"  labels v{v}: nucleus table unavailable ({repr(e)[:80]})")
            continue
        for tbl in PROOF_TABLES:
            try:
                p = client.materialize.query_table(tbl, filter_in_dict={"pt_root_id": n_v["pt_root_id"].tolist()},
                                                   materialization_version=int(v))
            except Exception:  # noqa: BLE001
                continue
            col = next((c for c in ("strategy_axon", "status_axon") if c in p.columns), None)
            if col is None:
                continue
            m_ = n_v.merge(p[["pt_root_id", col]], on="pt_root_id", how="left")
            labs.append(pd.DataFrame({"version": v, "nucleus_id": m_["id"], "axon_label": m_[col].astype(str),
                                      "table": tbl}))
            print(f"  labels v{v}: {tbl}.{col}, {m_[col].notna().sum()} of {len(m_)} labelled")
            break
    if labs:
        pd.concat(labs, ignore_index=True).to_parquet(lab_path)
    print(f"fetch done; {time.time() - t0:.0f}s")


# ------------------------------------------------------------------ stage 2
def build(args):
    from scipy.spatial import cKDTree
    from compute_proximity import AXON, DEND, PROX, load_skeleton, resample
    t0 = time.time()
    pres = ref_pres()
    meta = pd.read_csv(HIST / "versions.csv")
    files = {int(v): HIST / f"attached_v{int(v):04d}.parquet" for v in meta["version"]}
    files = {v: f for v, f in files.items() if f.exists()}
    att = {v: pd.read_parquet(f).set_index("synapse_id")["attached"] for v, f in files.items()}
    versions = sorted(att)

    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    usable["root_id"] = usable["root_id"].astype("int64")
    root_idx = pd.Series(np.arange(len(usable)), index=usable["root_id"])
    pairs = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        pairs[c] = pairs[c].astype("int64")
    io = pd.read_parquet(DATA / f"synapses_io_usable_v{VERSION}.parquet")
    io = io[(io["direction"] == "out") & io["root_id"].astype("int64").isin(pres["root_id"])]
    io["root_id"] = io["root_id"].astype("int64")
    syn = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        syn[c] = syn[c].astype("int64")
    syn = syn[(syn["pre_root_id"] != syn["post_root_id"]) & (syn["pre_axon_dist_um"] <= 3.0)
              & syn["pre_root_id"].isin(pres["root_id"]) & syn["post_root_id"].isin(usable["root_id"])]

    d_pts, d_len, d_owner = [], [], []
    for i, r in enumerate(usable["root_id"]):
        sk = load_skeleton(r)
        if sk is None:
            continue
        pts, ln, _ = resample(sk, DEND)
        d_pts.append(pts), d_len.append(ln), d_owner.append(np.full(len(ln), i))
    d_pts, d_len, d_owner = np.vstack(d_pts), np.concatenate(d_len), np.concatenate(d_owner)
    dtree = cKDTree(d_pts)
    n_post = len(usable)
    pre_ids = pres["root_id"].to_numpy()
    ld = {v: np.zeros((len(pre_ids), n_post), np.float32) for v in versions}
    ns = {v: np.zeros((len(pre_ids), n_post), np.int16) for v in versions}
    rows = []
    for j, r in enumerate(pre_ids):
        sk = load_skeleton(r)
        a_pts, a_len, _ = resample(sk, AXON)
        atree = cKDTree(a_pts)
        nbrs = atree.query_ball_tree(dtree, PROX)
        lens = np.fromiter(map(len, nbrs), int, len(nbrs))
        a_idx = np.repeat(np.arange(len(nbrs)), lens)
        d_idx = np.fromiter(itertools.chain.from_iterable(nbrs), int, lens.sum())
        other = d_owner[d_idx] != root_idx[r]
        a_idx, d_idx = a_idx[other], d_idx[other]
        uniq_d, inv = np.unique(d_idx, return_inverse=True)
        owner_u, len_u = d_owner[uniq_d], d_len[uniq_d]
        outs = io[io["root_id"] == r]
        o_xyz = outs[["ctr_x_nm", "ctr_y_nm", "ctr_z_nm"]].to_numpy(float)
        piece_out = cKDTree(o_xyz).query(a_pts)[1] if len(outs) else None   # nearest output synapse per piece
        s = syn[syn["pre_root_id"] == r]
        s_post = root_idx.reindex(s["post_root_id"]).to_numpy().astype(int)
        for v in versions:
            a_out = att[v].reindex(outs["synapse_id"]).fillna(False).to_numpy(bool)
            keep_a = a_out[piece_out] if piece_out is not None else np.zeros(len(a_len), bool)
            hit = np.zeros(len(uniq_d), bool)
            hit[inv[keep_a[a_idx]]] = True
            ld[v][j] = np.bincount(owner_u[hit], weights=len_u[hit], minlength=n_post) / 1e3
            keep_s = att[v].reindex(s["synapse_id"]).fillna(False).to_numpy(bool)
            ns[v][j] = np.bincount(s_post[keep_s], minlength=n_post)
            rows.append((r, v, a_out.mean() if len(a_out) else np.nan, a_len[keep_a].sum() / a_len.sum(),
                         keep_s.sum(), len(s)))
        if (j + 1) % 25 == 0:
            print(f"  axon {j + 1}/{len(pre_ids)}: {time.time() - t0:.0f}s", flush=True)

    pj = pd.Series(np.arange(len(pre_ids)), index=pre_ids)[pairs["pre_root_id"]].to_numpy()
    pk = root_idx[pairs["post_root_id"]].to_numpy()
    idx = [("radial_1.00_0", "radial", 1.0, 0)]
    for v in versions:
        name = f"hist_{v:04d}_0"
        pd.DataFrame({"row": np.arange(len(pairs), dtype=np.int32), "Ld_um": ld[v][pj, pk],
                      "Nsyn": ns[v][pj, pk]}).to_parquet(SET_DIR / f"{name}.parquet")
        idx.append((name, "hist", float(v), 0))
    pd.DataFrame(idx, columns=["set", "cond", "level", "rep"]).to_csv(SET_DIR / "index_hist.csv", index=False)

    ax = pd.DataFrame(rows, columns=["pre_root_id", "version", "outputs_attached", "cable_attached",
                                     "l2l_syn_attached", "l2l_syn_full"])
    ax = ax.merge(pres.rename(columns={"root_id": "pre_root_id"}), on="pre_root_id")
    lab_path = HIST / "labels_by_version.parquet"
    if lab_path.exists():
        ax = ax.merge(pd.read_parquet(lab_path)[["version", "nucleus_id", "axon_label"]], on=["version", "nucleus_id"],
                      how="left")
    ax.to_parquet(HIST / f"history_axons_v{VERSION}.parquet")
    meta = meta.set_index("version")
    summ = ax.groupby("version").agg(axons_complete=("cable_attached", lambda x: (x > 0.95).sum()),
                                     median_cable=("cable_attached", "median"),
                                     synapses_kept=("l2l_syn_attached", "sum"), synapses_full=("l2l_syn_full", "sum"))
    summ["date"] = pd.to_datetime(meta.loc[summ.index, "time_stamp"]).dt.strftime("%Y-%m-%d")
    if "axon_label" in ax:
        lab = ax.pivot_table(index="version", columns="axon_label", values="pre_root_id", aggfunc="count").fillna(0)
        summ = summ.join(lab.astype(int))
    (ROOT / "notebooks" / "08_history_summary.md").write_text(
        f"# Proofreading history of the {len(pre_ids)} fully extended axons (auto-generated by "
        f"scripts/edit_history.py build, v{VERSION})\n\nPer materialization version: axons with > 95% of today's "
        f"cable attached, median cable fraction attached, l2l synapses attached vs today, and axon labels at that "
        f"version where the tables could be queried.\n\n" + summ.to_markdown() + "\n")
    print(summ.to_string())
    print(f"saved {len(versions)} hist sets + index_hist.csv; {time.time() - t0:.0f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["fetch", "build"])
    ap.add_argument("--limit", type=int, default=0, help="fetch: first N axons only (test; writes data/history_test/)")
    ap.add_argument("--versions", default="", help="fetch: comma-separated versions only (test)")
    args = ap.parse_args()
    fetch(args) if args.stage == "fetch" else build(args)


if __name__ == "__main__":
    main()
