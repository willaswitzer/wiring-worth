"""Session 7: calibrated truncation -- incompleteness that matches real partial proofreading, per cell type.

Run from the repo root after scripts/truncate_like_to_like.py:
    .venv/bin/python scripts/calibrated_truncation.py [--reps 20] [--seed 11]
then fit with
    TRUNC_INDEX=index_calib.csv TRUNC_TAG=calib .r-env/bin/Rscript scripts/like_to_like_trunc.R fit poisson 4

Model (branch cuts). Each axon edge of a fully extended 23P or 4P axon is cut at random; a cut removes
the whole subtree beyond it (as when an unfinished reconstruction stops at a point on a branch). The cut
probability of an edge from path length d0 to d1 (from the soma) is 1 - S(d1)/S(d0), so a point at
path length d survives with probability S(d) / S(d_axon_start). S is the target survival profile:
the partial/full ratio of mean axon cable per neuron by path-length bin, measured on real partially vs
fully extended axons of the same cell type (Session 5, 05_results.md, section 2), at bin midpoints,
made non-increasing (running minimum), linearly interpolated, S(0) = 1:

    path bin (um)    0-200  200-400  400-600  600-800  800+
    23P              0.92   0.83     0.58     0.55     0.40
    4P               0.94   0.76     0.67     0.75*    0.69*     (* -> 0.67 after running minimum)

Other cell types (26 of 193 axons) have no calibration and are left complete.

Conditions written to data/truncation/sets_v1822/ (one set per rep), listed in index_calib.csv together
with the full-axon reference radial_1.00_0:
- calib     : branch cuts as above; Ld and synapses follow the kept cable.
- calibblk  : null with the same cable kept per axon as that rep's calib, removed as random 25 um cubes.
With --profile-boot (calibration uncertainty), each rep instead draws its own target profile: real
partial and fully extended axons of each cell type are resampled with replacement (within cell type x
axon level) and the partial/full ratio of mean cable per bin is recomputed from their skeletons; that
rep's cuts use the drawn profile. Sets are written as condition `calibunc` (no block null), listed in
index_calibunc.csv; the drawn profiles go to data/truncation/calibunc_profiles_v1822.csv.
Writes data/truncation/calib_axons_v1822.parquet (per axon x rep: cable kept, synapses kept) and
data/truncation/calib_profile_v1822.md (achieved vs. target kept-cable profile).
"""
import argparse
import itertools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from config import VERSION  # noqa: E402
from compute_proximity import AXON, DEND, PROX, STEP, load_skeleton, resample  # noqa: E402
from truncate_like_to_like import MAX_DIST_UM, block_order, blocks_keep  # noqa: E402

DATA = ROOT / "data"
OUT = DATA / "truncation"
SET_DIR = OUT / f"sets_v{VERSION}"
BIN_MID_UM = np.array([100, 300, 500, 700, 900])
PROFILE = {  # Session 5, partial/full ratio of mean axon cable per neuron by path bin
    "23P": np.array([0.92, 0.83, 0.58, 0.55, 0.40]),
    "4P": np.array([0.94, 0.76, 0.67, 0.75, 0.69]),
}
BINS_UM = [0, 200, 400, 600, 800, np.inf]


def survival(ct, profile=None):
    s = np.minimum.accumulate(PROFILE[ct] if profile is None else profile)
    x, y = np.r_[0, BIN_MID_UM], np.r_[1.0, s]
    return lambda d_um: np.interp(d_um, x, y)  # flat beyond the last midpoint


def neuron_bin_cable(usable):
    """Axon cable (um) per path-length bin for every usable 23P/4P neuron with a full or partial axon."""
    u = usable[usable["cell_type"].isin(list(PROFILE)) &
               usable["axon_level"].isin(["axon_fully_extended", "axon_partially_extended"])]
    rows = []
    for r in u.itertuples():
        sk = load_skeleton(r.root_id)
        if sk is None:
            continue
        _, ln, pth = resample(sk, AXON)
        rows.append((r.cell_type, r.axon_level, *(np.histogram(pth / 1e3, bins=BINS_UM, weights=ln)[0] / 1e3)))
    return pd.DataFrame(rows, columns=["cell_type", "axon_level"] + [f"b{i}" for i in range(len(BINS_UM) - 1)])


def draw_profile(cab, ct, rng):
    """Partial/full ratio of mean cable per bin, neurons resampled within cell type x axon level."""
    cols = [c for c in cab.columns if c.startswith("b")]
    means = {}
    for lvl in ("axon_fully_extended", "axon_partially_extended"):
        x = cab[(cab["cell_type"] == ct) & (cab["axon_level"] == lvl)][cols].to_numpy()
        means[lvl] = x[rng.integers(0, len(x), len(x))].mean(axis=0)
    return means["axon_partially_extended"] / means["axon_fully_extended"]


def resample_with_edges(sk, compartment):
    """As compute_proximity.resample, also returning each piece's edge index."""
    v, e, comp = sk["vertices"], sk["edges"], sk["compartment"]
    child, parent = e[:, 0], e[:, 1]
    seg = np.linalg.norm(v[child] - v[parent], axis=1)
    keep = (comp[child] == compartment) & (seg > 0)
    idx = np.flatnonzero(keep)
    n = np.maximum(np.ceil(seg[idx] / STEP), 1).astype(int)
    return np.repeat(idx, n)


def vertex_paths(sk):
    v, e = sk["vertices"], sk["edges"]
    seg = np.linalg.norm(v[e[:, 0]] - v[e[:, 1]], axis=1)
    g = coo_matrix((np.r_[seg, seg], (np.r_[e[:, 0], e[:, 1]], np.r_[e[:, 1], e[:, 0]])), shape=(len(v),) * 2).tocsr()
    return dijkstra(g, directed=False, indices=int(sk["root"]))


def branch_cut_keep(sk, path_nm, S, rng):
    """Kept flag per edge (child side) under random branch cuts on axon edges."""
    e, comp = sk["edges"], sk["compartment"]
    child, parent = e[:, 0], e[:, 1]
    d0, d1 = path_nm[parent] / 1e3, path_nm[child] / 1e3
    s0, s1 = S(np.minimum(d0, d1)), S(np.maximum(d0, d1))
    p_cut = np.where(comp[child] == AXON, 1 - s1 / np.maximum(s0, 1e-12), 0.0)
    cut = rng.random(len(e)) < p_cut
    # propagate from the root: a vertex is kept iff no cut on its path to the root
    n_v = len(sk["vertices"])
    par = np.full(n_v, -1)
    par[child] = parent
    cut_v = np.zeros(n_v, bool)
    cut_v[child] = cut
    order = np.argsort(path_nm, kind="stable")  # parents are closer to the root than children
    kept_v = np.ones(n_v, bool)
    for u in order:
        p = par[u]
        if p >= 0:
            kept_v[u] = kept_v[p] and not cut_v[u]
    return kept_v[child]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--block-um", type=float, default=25.0)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--profile-boot", action="store_true", help="draw the target profile per rep (calibration uncertainty)")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    t0 = time.time()

    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    usable["root_id"] = usable["root_id"].astype("int64")
    pairs = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        pairs[c] = pairs[c].astype("int64")
    root_idx = pd.Series(np.arange(len(usable)), index=usable["root_id"])
    pre_meta = pairs.groupby("pre_root_id")["pre_cell_type"].first()
    pre_ids = pre_meta.index.to_numpy()
    syn = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        syn[c] = syn[c].astype("int64")
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
    print(f"dendrite tree: {len(d_pts):,} pieces; {time.time() - t0:.0f}s", flush=True)

    n_post = len(usable)
    conds = ("calibunc",) if args.profile_boot else ("calib", "calibblk")
    ld = {(c, r): np.zeros((len(pre_ids), n_post), np.float32) for c in conds for r in range(args.reps)}
    ns = {k: np.zeros((len(pre_ids), n_post), np.int16) for k in ld}
    axon_rows, prof_rows = [], []
    if args.profile_boot:
        cab = neuron_bin_cable(usable)
        drawn = {(ct, rep): draw_profile(cab, ct, rng) for rep in range(args.reps) for ct in PROFILE}
        surv_rep = {k: survival(k[0], v) for k, v in drawn.items()}
        pd.DataFrame([(ct, rep, *v) for (ct, rep), v in drawn.items()],
                     columns=["cell_type", "rep"] + [f"bin_{lo}" for lo in BINS_UM[:-1]]).to_csv(
            OUT / f"calibunc_profiles_v{VERSION}.csv", index=False)
        print(f"drew {args.reps} profiles per cell type; {time.time() - t0:.0f}s", flush=True)
    else:
        surv_rep = {(ct, rep): survival(ct) for ct in PROFILE for rep in range(args.reps)}

    for j, r in enumerate(pre_ids):
        ct = pre_meta[r]
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
        full_ns = np.bincount(s_post, minlength=n_post)
        for rep in range(args.reps):
            if ct in PROFILE:
                keep_a = branch_cut_keep(sk, path_nm, surv_rep[(ct, rep)], rng)[a_edge]
                frac = a_len[keep_a].sum() / total
                ld[(conds[0], rep)][j] = ld_row(keep_a)
                ns[(conds[0], rep)][j] = np.bincount(s_post[keep_a[s_near]], minlength=n_post)
                for lo, hi in zip(BINS_UM[:-1], BINS_UM[1:]):
                    m = (a_path / 1e3 >= lo) & (a_path / 1e3 < hi)
                    prof_rows.append((ct, r, rep, lo, a_len[m].sum(), a_len[m & keep_a].sum()))
                keep_b = keep_a
                if not args.profile_boot:
                    bid, perm, blen = block_order(a_pts, a_len, args.block_um * 1e3, rng)
                    keep_b = blocks_keep(bid, perm, blen, frac)
                    ld[("calibblk", rep)][j] = ld_row(keep_b)
                    ns[("calibblk", rep)][j] = np.bincount(s_post[keep_b[s_near]], minlength=n_post)
                axon_rows.append((r, ct, rep, frac, keep_a[s_near].sum(), keep_b[s_near].sum(), len(s),
                                  a_len[keep_b].sum() / total))
            else:  # no calibration available: left complete in both conditions
                for c in conds:
                    ld[(c, rep)][j], ns[(c, rep)][j] = full_ld, full_ns
                axon_rows.append((r, ct, rep, 1.0, len(s), len(s), len(s), 1.0))
        if (j + 1) % 25 == 0:
            print(f"  axon {j + 1}/{len(pre_ids)}: {time.time() - t0:.0f}s", flush=True)

    # write sets aligned with l2l_pairs rows (same format as truncate_like_to_like.py)
    pj = pd.Series(np.arange(len(pre_ids)), index=pre_ids)[pairs["pre_root_id"]].to_numpy()
    pk = root_idx[pairs["post_root_id"]].to_numpy()
    rows = np.arange(len(pairs), dtype=np.int32)
    index = [("radial_1.00_0", "radial", 1.0, 0)]
    for (c, rep) in ld:
        name = f"{c}_0.00_{rep}"
        pd.DataFrame({"row": rows, "Ld_um": ld[(c, rep)][pj, pk], "Nsyn": ns[(c, rep)][pj, pk]}).to_parquet(
            SET_DIR / f"{name}.parquet")
        index.append((name, c, 0.0, rep))
    idx_name = "index_calibunc.csv" if args.profile_boot else "index_calib.csv"
    pd.DataFrame(index, columns=["set", "cond", "level", "rep"]).to_csv(SET_DIR / idx_name, index=False)

    ax = pd.DataFrame(axon_rows, columns=["pre_root_id", "cell_type", "rep", "cable_frac", "n_syn_calib",
                                          "n_syn_calibblk", "n_syn_full", "cable_frac_calibblk"])
    tag = "calibunc" if args.profile_boot else "calib"
    ax.to_parquet(OUT / f"{tag}_axons_v{VERSION}.parquet")
    prof = pd.DataFrame(prof_rows, columns=["cell_type", "pre_root_id", "rep", "bin_lo_um", "full_nm", "kept_nm"])
    # achieved partial/full ratio of mean cable per neuron, per bin (the Session 5 statistic)
    g = prof.groupby(["cell_type", "bin_lo_um"]).agg(full=("full_nm", "sum"), kept=("kept_nm", "sum"))
    g["achieved"] = (g["kept"] / g["full"]).round(2)
    g["target"] = [PROFILE[ct][BINS_UM.index(lo)] for ct, lo in g.index]
    by_ct = ax.groupby("cell_type").agg(axons=("pre_root_id", "nunique"), cable_kept=("cable_frac", "mean"),
                                        cable_kept_blk=("cable_frac_calibblk", "mean"),
                                        syn_kept=("n_syn_calib", "sum"), syn_kept_blk=("n_syn_calibblk", "sum"),
                                        syn_full=("n_syn_full", "sum"))
    by_ct["syn_frac"] = by_ct["syn_kept"] / by_ct["syn_full"]
    by_ct["syn_frac_blk"] = by_ct["syn_kept_blk"] / by_ct["syn_full"]
    by_ct = by_ct.drop(columns=["syn_kept", "syn_kept_blk"]).round(3)
    txt = (f"# Calibrated truncation (auto-generated by scripts/calibrated_truncation.py, v{VERSION})\n\n"
           f"{args.reps} reps; seed {args.seed}; block {args.block_um:g} um. Targets: Session 5 partial/full ratio "
           f"of mean axon cable per neuron (4P made non-increasing for the cut model).\n\n"
           f"## Achieved vs. target kept-cable ratio by path bin (all reps)\n\n{g[['achieved', 'target']].to_markdown()}\n\n"
           f"## Per cell type (means over reps)\n\n{by_ct.to_markdown()}\n")
    (OUT / f"{tag}_profile_v{VERSION}.md").write_text(txt)
    print(txt)
    print(f"saved {len(conds) * args.reps} sets + {idx_name} to {SET_DIR.relative_to(ROOT)}/; {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
