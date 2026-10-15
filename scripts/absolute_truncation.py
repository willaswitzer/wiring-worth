"""Session 11: truncation at one absolute radius for every axon (C. Schneider-Mizell's suggestion).

Partial axon extension in MICrONS followed axons out to a fixed distance around the densely proofread column,
not to a per-axon fraction of cable. The Session 7 radial truncation keeps a fraction f of each axon's own
cable, so every axon gets a different cutoff. This script writes two conditions with one cutoff for all axons:

- pathabs : keep axon cable within path length R of the soma (R = 150, 250, 400, 600 um).
- cyl     : keep axon cable inside a vertical cylinder of radius Rc around the column axis
            (Rc = 100, 150, 200, 300 um), and only cable whose whole path back to the soma stays inside
            (once an axon leaves the cylinder, everything beyond is lost, as when a proofreader stops at the
            border). The first 30 um of path is always kept. The axis is the line along y (cortical depth)
            through the median (x, z) of the 193 presynaptic somata; its position is written to the summary.
            Axons whose soma lies more than 100 um from the axis (mostly HVA cells outside the column) are left
            complete, and counted in the summary.

Ld and Nsyn are recomputed exactly as in truncate_like_to_like.py. Run from the repo root:
    .venv/bin/python scripts/absolute_truncation.py
then fit/shuffle/bootstrap with like_to_like_trunc.R using TRUNC_INDEX=index_pathabs.csv TRUNC_TAG=pathabs
(and index_cyl.csv / cyl); see scripts/run_session11.sh.
Writes data/truncation/sets_v1822/{pathabs,cyl}_<level>_0.parquet, index_pathabs.csv, index_cyl.csv,
data/truncation/abs_axons_v1822.parquet and data/truncation/abs_summary_v1822.md.
"""
import argparse
import itertools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from config import VERSION  # noqa: E402
from compute_proximity import AXON, DEND, PROX, load_skeleton, resample  # noqa: E402
from calibrated_truncation import resample_with_edges, vertex_paths  # noqa: E402

DATA = ROOT / "data"
OUT = DATA / "truncation"
SET_DIR = OUT / f"sets_v{VERSION}"
MAX_DIST_UM = 3.0
PATH_R_UM = [150, 250, 400, 600]
CYL_R_UM = [100, 150, 200, 300]
CYL_SOMA_MAX_UM = 100   # axons with somata farther than this from the axis are left complete
PROXIMAL_UM = 30        # always kept (axon initial segment region)


def cylinder_keep_edges(sk, path_nm, axis_xz_nm, r_nm):
    """Kept flag per edge (child side): child vertex and all its ancestors inside the cylinder (or proximal)."""
    v, e = sk["vertices"], sk["edges"]
    child, parent = e[:, 0], e[:, 1]
    lat = np.hypot(v[:, 0] - axis_xz_nm[0], v[:, 2] - axis_xz_nm[1])
    ok_v = (lat <= r_nm) | (path_nm <= PROXIMAL_UM * 1e3)
    par = np.full(len(v), -1)
    par[child] = parent
    kept_v = ok_v.copy()
    for u in np.argsort(path_nm, kind="stable"):  # parents are closer to the root than children
        p = par[u]
        if p >= 0:
            kept_v[u] = kept_v[p] and ok_v[u]
    return kept_v[child]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="first N axons only (timing test); writes nothing")
    args = ap.parse_args()
    t0 = time.time()

    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    usable["root_id"] = usable["root_id"].astype("int64")
    neurons = pd.read_parquet(DATA / f"l2l_neurons_v{VERSION}.parquet")
    neurons["root_id"] = neurons["root_id"].astype("int64")
    pairs = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet")
    pairs["pre_root_id"] = pairs["pre_root_id"].astype("int64")
    pairs["post_root_id"] = pairs["post_root_id"].astype("int64")
    root_idx = pd.Series(np.arange(len(usable)), index=usable["root_id"])
    pre_ids = pairs["pre_root_id"].unique()
    soma = neurons.set_index("root_id").loc[pre_ids, ["soma_x_um", "soma_z_um"]]
    axis_um = soma.median().to_numpy()
    soma_lat = np.hypot(soma["soma_x_um"] - axis_um[0], soma["soma_z_um"] - axis_um[1])
    in_col = (soma_lat <= CYL_SOMA_MAX_UM).to_dict()
    if args.limit:
        pre_ids = pre_ids[: args.limit]

    syn = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    syn["pre_root_id"] = syn["pre_root_id"].astype("int64")
    syn["post_root_id"] = syn["post_root_id"].astype("int64")
    syn = syn[(syn["pre_root_id"] != syn["post_root_id"]) & (syn["pre_axon_dist_um"] <= MAX_DIST_UM)
              & syn["pre_root_id"].isin(pre_ids) & syn["post_root_id"].isin(usable["root_id"])]

    d_pts, d_len, d_owner = [], [], []
    for i, r in enumerate(usable["root_id"]):
        sk = load_skeleton(r)
        if sk is None:
            continue
        pts, ln, _ = resample(sk, DEND)
        d_pts.append(pts), d_len.append(ln), d_owner.append(np.full(len(ln), i))
    d_pts, d_len, d_owner = np.vstack(d_pts), np.concatenate(d_len), np.concatenate(d_owner)
    dtree = cKDTree(d_pts)
    print(f"dendrite tree: {len(d_pts):,} pieces; column axis (x, z) = ({axis_um[0]:.0f}, {axis_um[1]:.0f}) um; "
          f"{time.time() - t0:.0f}s", flush=True)

    sets = [("pathabs", R) for R in PATH_R_UM] + [("cyl", R) for R in CYL_R_UM]
    n_post = len(usable)
    ld = {s: np.zeros((len(pre_ids), n_post), np.float32) for s in sets}
    ns = {s: np.zeros((len(pre_ids), n_post), np.int16) for s in sets}
    axon_rows = []
    max_ld_err = 0.0

    for j, r in enumerate(pre_ids):
        sk = load_skeleton(r)
        a_pts, a_len, a_path = resample(sk, AXON)
        a_edge = resample_with_edges(sk, AXON)
        assert len(a_edge) == len(a_len)
        path_nm = vertex_paths(sk)
        atree = cKDTree(a_pts)
        me = root_idx[r]
        nbrs = atree.query_ball_tree(dtree, PROX)
        lens = np.fromiter(map(len, nbrs), int, len(nbrs))
        a_idx = np.repeat(np.arange(len(nbrs)), lens)
        d_idx = np.fromiter(itertools.chain.from_iterable(nbrs), int, lens.sum())
        other = d_owner[d_idx] != me
        a_idx, d_idx = a_idx[other], d_idx[other]
        uniq_d, inv = np.unique(d_idx, return_inverse=True)
        owner_u, len_u = d_owner[uniq_d], d_len[uniq_d]

        def ld_row(keep_a):
            hit = np.zeros(len(uniq_d), bool)
            hit[inv[keep_a[a_idx]]] = True
            return np.bincount(owner_u[hit], weights=len_u[hit], minlength=n_post) / 1e3

        s = syn[syn["pre_root_id"] == r]
        s_post = root_idx.reindex(s["post_root_id"]).to_numpy().astype(int)
        s_near = atree.query(s[["ctr_x_nm", "ctr_y_nm", "ctr_z_nm"]].to_numpy(float))[1] if len(s) else np.zeros(0, int)
        total = a_len.sum()
        full_ld = ld_row(np.ones(len(a_len), bool))
        ref = pairs.loc[pairs["pre_root_id"] == r].set_index("post_root_id")["Ld_um"]
        max_ld_err = max(max_ld_err, np.abs(full_ld[root_idx[ref.index].to_numpy()] - ref.to_numpy()).max())

        for cond, R in sets:
            if cond == "pathabs":
                keep_a = a_path <= R * 1e3
            elif in_col[r]:
                keep_a = cylinder_keep_edges(sk, path_nm, axis_um * 1e3, R * 1e3)[a_edge]
            else:
                keep_a = np.ones(len(a_len), bool)
            keep_s = keep_a[s_near]
            ld[(cond, R)][j] = ld_row(keep_a)
            ns[(cond, R)][j] = np.bincount(s_post[keep_s], minlength=n_post)
            axon_rows.append((r, cond, R, bool(in_col[r]), a_len[keep_a].sum() / total, int(keep_s.sum()), len(s)))
        if (j + 1) % 25 == 0:
            print(f"  axon {j + 1}/{len(pre_ids)}: {time.time() - t0:.0f}s", flush=True)

    print(f"max |Ld recomputed - l2l_pairs Ld| at full cable: {max_ld_err:.4f} um")
    assert max_ld_err < 0.01
    ax = pd.DataFrame(axon_rows, columns=["pre_root_id", "cond", "level", "in_column", "cable_frac",
                                          "n_syn_kept", "n_syn_full"])
    if args.limit:
        print(ax.groupby(["cond", "level"])[["cable_frac"]].median())
        return

    sel = pairs["pre_root_id"].isin(pre_ids).to_numpy()
    rows = np.flatnonzero(sel)
    pj = pd.Series(np.arange(len(pre_ids)), index=pre_ids)[pairs.loc[sel, "pre_root_id"]].to_numpy()
    pk = root_idx[pairs.loc[sel, "post_root_id"]].to_numpy()
    for cond in ("pathabs", "cyl"):
        index = [("radial_1.00_0", "radial", 1.0, 0)]
        for c, R in sets:
            if c != cond:
                continue
            name = f"{c}_{R:.2f}_0"
            pd.DataFrame({"row": rows.astype(np.int32), "Ld_um": ld[(c, R)][pj, pk],
                          "Nsyn": ns[(c, R)][pj, pk]}).to_parquet(SET_DIR / f"{name}.parquet")
            index.append((name, c, float(R), 0))
        pd.DataFrame(index, columns=["set", "cond", "level", "rep"]).to_csv(SET_DIR / f"index_{cond}.csv", index=False)
    ax.to_parquet(OUT / f"abs_axons_v{VERSION}.parquet")

    meta = neurons.set_index("root_id")
    ax["cell_type"] = meta.loc[ax["pre_root_id"], "cell_type"].to_numpy()
    ax["area"] = meta.loc[ax["pre_root_id"], "area"].to_numpy()
    v1 = ax[ax["area"] == "V1"]
    g = v1.groupby(["cond", "level"]).agg(median_cable_kept=("cable_frac", "median"),
                                          synapses_kept=("n_syn_kept", "sum"), synapses_full=("n_syn_full", "sum"))
    g["synapses_kept_frac"] = g["synapses_kept"] / g["synapses_full"]
    by_ct = v1[v1["cell_type"].isin(["23P", "4P"])].pivot_table(index=["cond", "level"], columns="cell_type",
                                                                 values="cable_frac", aggfunc="median")
    by_ct.columns = [f"median cable kept, {c}" for c in by_ct.columns]
    n_out = int((~soma_lat.le(CYL_SOMA_MAX_UM)).sum())
    (OUT / f"abs_summary_v{VERSION}.md").write_text(
        f"# Absolute-radius truncation sets (auto-generated by scripts/absolute_truncation.py, v{VERSION})\n\n"
        f"Column axis (x, z) = ({axis_um[0]:.0f}, {axis_um[1]:.0f}) um, the median of the {len(soma)} presynaptic somata. "
        f"{n_out} axons have somata more than {CYL_SOMA_MAX_UM} um from it and are left complete under `cyl` "
        f"({int(((~soma_lat.le(CYL_SOMA_MAX_UM)) & (meta.loc[soma.index, 'area'] == 'V1')).sum())} of them in V1).\n\n"
        f"V1 presynaptic axons:\n\n" + g.join(by_ct).round(3).to_markdown() + "\n")
    print(f"saved {len(sets)} sets + index_pathabs.csv, index_cyl.csv; {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
