"""Session 7: simulate axon incompleteness on the fully extended reference set.

Run from the repo root after scripts/prepare_like_to_like.py:
    .venv/bin/python scripts/truncate_like_to_like.py [--reps 5] [--block-um 25] [--seed 7]

For each of the 193 fully extended presynaptic axons, and each retention level f (fraction of the
axon's own cable kept), three conditions remove cable/synapses; Ld and Nsyn are then recomputed
for every pair in data/l2l_pairs_v1822.parquet:

- radial  : keep axon cable within path length r_f of the soma, where r_f is chosen per axon so
            that a fraction f of its cable is kept. Synapses are kept if their nearest axon piece
            is kept (equivalent to pre_axon_path_um <= r_f). Ld is recomputed from the kept cable.
            This is the Session 5 picture of how 23P partial proofreading ends an axon.
- blocks  : null with the same cable loss at random positions. The axon is cut into cubes of
            --block-um (random grid offset per rep); cubes are removed in random order until a
            fraction 1 - f of cable is gone. Synapses and Ld follow the kept cable, as above.
- synsub  : null of the project plan. Each axon's synapses are subsampled at random to the count
            that radial truncation keeps at the same f; Ld is left complete.

Levels are nested within a rep (each lower f removes a superset), so curves are smooth.
The like-to-like models are refitted on every (condition, level, rep) set by
scripts/like_to_like_trunc.R.

Ld uses the Session 4/5 definition (compute_proximity.py): dendrite pieces of the post within
5 um of the (kept) pre axon, summed length. The full-cable recomputation is checked against
l2l_pairs Ld_um before anything is written.

Outputs (data/truncation/, not committed):
- sets_v1822/<cond>_<level>_<rep>.parquet: row (0-based index into l2l_pairs), Ld_um, Nsyn;
  sets_v1822/index.csv lists the sets
- trunc_axons_v1822.parquet: per axon x set: cable kept (frac), synapses kept, r_um (radial)
- trunc_summary_v1822.md   : totals per condition x level
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

DATA = ROOT / "data"
OUT = DATA / "truncation"
LEVELS = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
MAX_DIST_UM = 3.0  # synapse exclusion, as in prepare_like_to_like.py


def radial_keep(a_len, a_path, f):
    """Mask of axon pieces within the path length that keeps fraction f of cable; and that length."""
    if f >= 1:
        return np.ones(len(a_len), bool), np.inf
    order = np.argsort(a_path, kind="stable")
    cum = np.cumsum(a_len[order])
    k = np.searchsorted(cum, f * cum[-1])
    r = a_path[order][min(k, len(order) - 1)]
    return a_path <= r, r


def block_order(a_pts, a_len, block_nm, rng):
    """Block id per axon piece and a random removal order of blocks."""
    off = rng.uniform(0, block_nm, 3)
    vox = np.floor((a_pts + off) / block_nm).astype(np.int64)
    _, bid = np.unique(vox, axis=0, return_inverse=True)
    bid = bid.ravel()
    blen = np.bincount(bid, weights=a_len)
    return bid, rng.permutation(len(blen)), blen


def blocks_keep(bid, perm, blen, f):
    if f >= 1:
        return np.ones(len(bid), bool)
    removed_cum = np.cumsum(blen[perm])
    n_rm = np.searchsorted(removed_cum, (1 - f) * blen.sum()) + 1
    gone = np.zeros(len(blen), bool)
    gone[perm[:n_rm]] = True
    return ~gone[bid]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--block-um", type=float, default=25.0)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--limit", type=int, default=0, help="first N axons only (timing test)")
    ap.add_argument("--only-blocks", default="", metavar="NAME",
                    help="sensitivity run: write only the block-null sets, as condition NAME (e.g. blocks10), "
                         "plus radial_1.00_0, listed in index_NAME.csv; the main sets and tables are untouched")
    args = ap.parse_args()
    block_nm = args.block_um * 1e3
    t0 = time.time()
    OUT.mkdir(exist_ok=True)

    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    usable["root_id"] = usable["root_id"].astype("int64")
    pairs = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet")
    pairs["pre_root_id"] = pairs["pre_root_id"].astype("int64")
    pairs["post_root_id"] = pairs["post_root_id"].astype("int64")
    root_idx = pd.Series(np.arange(len(usable)), index=usable["root_id"])
    pre_ids = pairs["pre_root_id"].unique()
    if args.limit:
        pre_ids = pre_ids[: args.limit]

    syn = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    syn["pre_root_id"] = syn["pre_root_id"].astype("int64")
    syn["post_root_id"] = syn["post_root_id"].astype("int64")
    syn = syn[(syn["pre_root_id"] != syn["post_root_id"]) & (syn["pre_axon_dist_um"] <= MAX_DIST_UM)
              & syn["pre_root_id"].isin(pre_ids) & syn["post_root_id"].isin(usable["root_id"])]

    # Dendrites of all usable neurons -> one KD-tree (as compute_proximity.py)
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

    # Set definitions: (cond, level, rep)
    sets = [("radial", f, 0) for f in LEVELS]
    for cond in ("blocks", "synsub"):
        sets += [(cond, f, rep) for rep in range(args.reps) for f in LEVELS if f < 1]
    n_post = len(usable)
    ld = {s: np.zeros((len(pre_ids), n_post), np.float32) for s in sets}
    ns = {s: np.zeros((len(pre_ids), n_post), np.int16) for s in sets}
    axon_rows = []
    rng = np.random.default_rng(args.seed)
    max_ld_err = 0.0

    for j, r in enumerate(pre_ids):
        sk = load_skeleton(r)
        a_pts, a_len, a_path = resample(sk, AXON)
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
            """Ld (um) to every usable neuron from the kept axon pieces."""
            hit = np.zeros(len(uniq_d), bool)
            hit[inv[keep_a[a_idx]]] = True
            return np.bincount(owner_u[hit], weights=len_u[hit], minlength=n_post) / 1e3

        s = syn[syn["pre_root_id"] == r]
        s_post = root_idx.reindex(s["post_root_id"]).to_numpy()
        s_near = np.zeros(0, int)
        if len(s):
            _, s_near = atree.query(s[["ctr_x_nm", "ctr_y_nm", "ctr_z_nm"]].to_numpy(float))
            path_err = np.abs(a_path[s_near] / 1e3 - s["pre_axon_path_um"].to_numpy())
            assert path_err.max() < 1e-3, f"{r}: synapse path mismatch {path_err.max()}"

        def ns_row(keep_syn):
            return np.bincount(s_post[keep_syn].astype(int), minlength=n_post)

        total = a_len.sum()
        full_ld = ld_row(np.ones(len(a_len), bool))
        ref = pairs.loc[pairs["pre_root_id"] == r].set_index("post_root_id")["Ld_um"]
        max_ld_err = max(max_ld_err, np.abs(full_ld[root_idx[ref.index].to_numpy()] - ref.to_numpy()).max())

        # radial
        n_radial = {}
        for f in LEVELS:
            keep_a, r_f = radial_keep(a_len, a_path, f)
            keep_s = keep_a[s_near]
            n_radial[f] = keep_s.sum()
            key = ("radial", f, 0)
            ld[key][j], ns[key][j] = ld_row(keep_a), ns_row(keep_s)
            axon_rows.append((r, *key, a_len[keep_a].sum() / total, keep_s.sum(), len(s), r_f / 1e3))
        for rep in range(args.reps):
            # blocks
            bid, perm, blen = block_order(a_pts, a_len, block_nm, rng)
            for f in LEVELS[1:]:
                keep_a = blocks_keep(bid, perm, blen, f)
                keep_s = keep_a[s_near]
                key = ("blocks", f, rep)
                ld[key][j], ns[key][j] = ld_row(keep_a), ns_row(keep_s)
                axon_rows.append((r, *key, a_len[keep_a].sum() / total, keep_s.sum(), len(s), np.nan))
            # synsub: random synapses, radial count, full Ld
            sperm = rng.permutation(len(s))
            for f in LEVELS[1:]:
                keep_s = np.zeros(len(s), bool)
                keep_s[sperm[: n_radial[f]]] = True
                key = ("synsub", f, rep)
                ld[key][j], ns[key][j] = full_ld, ns_row(keep_s)
                axon_rows.append((r, *key, 1.0, keep_s.sum(), len(s), np.nan))
        if (j + 1) % 25 == 0:
            print(f"  axon {j + 1}/{len(pre_ids)}: {time.time() - t0:.0f}s", flush=True)

    print(f"max |Ld recomputed - l2l_pairs Ld| at full cable: {max_ld_err:.4f} um")
    assert max_ld_err < 0.01

    # One file per set, rows aligned with l2l_pairs (the R side reads one set at a time)
    sel = pairs["pre_root_id"].isin(pre_ids).to_numpy()
    rows = np.flatnonzero(sel)
    pj = pd.Series(np.arange(len(pre_ids)), index=pre_ids)[pairs.loc[sel, "pre_root_id"]].to_numpy()
    pk = root_idx[pairs.loc[sel, "post_root_id"]].to_numpy()
    full_ns = ns[("radial", 1.0, 0)][pj, pk]
    assert np.array_equal(full_ns, pairs.loc[sel, "Nsyn"].to_numpy()), "Nsyn at full cable != l2l_pairs"
    suffix = f"_limit{args.limit}" if args.limit else ""
    set_dir = OUT / f"sets_v{VERSION}{suffix}"
    set_dir.mkdir(exist_ok=True)
    if args.only_blocks:
        keep_sets = [k for k in sets if k == ("radial", 1.0, 0) or k[0] == "blocks"]
        index = []
        for key in keep_sets:
            cond = args.only_blocks if key[0] == "blocks" else key[0]
            name = f"{cond}_{key[1]:.2f}_{key[2]}"
            if cond != "radial":  # radial_1.00_0 already exists from the main run and is identical
                pd.DataFrame({"row": rows.astype(np.int32), "Ld_um": ld[key][pj, pk],
                              "Nsyn": ns[key][pj, pk]}).to_parquet(set_dir / f"{name}.parquet")
            index.append((name, cond, key[1], key[2]))
        pd.DataFrame(index, columns=["set", "cond", "level", "rep"]).to_csv(
            set_dir / f"index_{args.only_blocks}.csv", index=False)
        print(f"saved {len(index) - 1} {args.only_blocks} sets + index_{args.only_blocks}.csv; {time.time() - t0:.0f}s")
        return
    index, parts = [], []
    for key in sets:
        name = f"{key[0]}_{key[1]:.2f}_{key[2]}"
        df = pd.DataFrame({"row": rows.astype(np.int32), "Ld_um": ld[key][pj, pk], "Nsyn": ns[key][pj, pk]})
        df.to_parquet(set_dir / f"{name}.parquet")
        index.append((name, *key))
        parts.append(df.assign(cond=key[0], level=np.float32(key[1]), rep=np.int8(key[2])))
    pd.DataFrame(index, columns=["set", "cond", "level", "rep"]).to_csv(set_dir / "index.csv", index=False)
    out = pd.concat(parts, ignore_index=True)
    ax = pd.DataFrame(axon_rows, columns=["pre_root_id", "cond", "level", "rep", "cable_frac",
                                          "n_syn_kept", "n_syn_full", "r_um"])
    ax.to_parquet(OUT / f"trunc_axons_v{VERSION}{suffix}.parquet")

    # Summary
    per_rep = out.assign(Ld_gt0=out["Ld_um"] > 0).groupby(["cond", "level", "rep"], observed=True).agg(
        Nsyn=("Nsyn", "sum"), Ld_mm=("Ld_um", "sum"), pairs_Ld_gt0=("Ld_gt0", "sum"))
    summ = per_rep.groupby(["cond", "level"], observed=True).mean()
    summ["Ld_mm"] /= 1e3
    ax["level"] = ax["level"].astype(np.float32)
    summ = summ.join(ax.groupby(["cond", "level"])["cable_frac"].median().rename("median_cable_frac"))
    summ["Nsyn_frac"] = summ["Nsyn"] / summ.loc[("radial", np.float32(1.0)), "Nsyn"]
    summ["Ld_frac"] = summ["Ld_mm"] / summ.loc[("radial", np.float32(1.0)), "Ld_mm"]
    print(summ.round(3).to_string())
    (OUT / f"trunc_summary_v{VERSION}{suffix}.md").write_text(
        f"# Truncation sets (auto-generated by scripts/truncate_like_to_like.py, v{VERSION})\n\n"
        f"{len(pre_ids)} fully extended axons; reps {args.reps}; block {args.block_um:g} um; seed {args.seed}. "
        f"Means over reps.\n\n" + summ.round(3).to_markdown() + "\n")
    print(f"saved {len(sets)} sets to {set_dir.relative_to(ROOT)}/; "
          f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
