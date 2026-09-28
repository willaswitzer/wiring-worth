"""Session 4: axon path length of synapses, and axon-dendrite proximity (Ld) for all pairs.

Run from the repo root after scripts/fetch_skeletons.py and scripts/pull_synapses.py:
    .venv/bin/python scripts/compute_proximity.py [--labels synapse|service]

--labels synapse (default, Session 5): axon/dendrite labels from scripts/relabel_skeletons.py;
  writes the main outputs below and notebooks/05_proximity_counts.md.
--labels service (Session 4): the skeleton service's labels; writes *_servicelabels.parquet
  and notebooks/04_proximity_counts.md.

Method (following Ding et al. 2025):
- Skeletons (nm) are resampled so every edge is <= 1 um: each edge (child -> parent) is cut into
  ceil(L / 1 um) equal pieces, and each piece is represented by its midpoint and its length.
  A piece belongs to the child vertex's compartment (2 axon, 3 dendrite).
- Proximity for a pre x post pair: resampled pre-axon points and post-dendrite points within 5 um
  (scipy cKDTree.query_ball_tree). Ld = summed length of the post's dendritic pieces in proximity.
- A synapse is in a proximity if it is within 3 um of any proximity point (axon or dendrite) of
  its own pre x post pair.
- Axon path length of a synapse: soma-to-point path length along the pre skeleton, at the nearest
  resampled pre-axon point.
Pairs: every presynaptic usable neuron (fully/partially extended axon) x every other usable neuron.
"""
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
from config import VERSION  # noqa: E402

DATA = ROOT / "data"
SKEL_DIR = DATA / "skeletons"
SKEL_VERSION = 4
SYN_PATH = DATA / f"synapses_usable_v{VERSION}.parquet"  # input (synapse list); computed columns are replaced
LABELS = "synapse"
OUTPUTS = {
    "synapse": (SYN_PATH, DATA / f"pairs_usable_v{VERSION}.parquet", ROOT / "notebooks" / "05_proximity_counts.md"),
    "service": (DATA / f"synapses_usable_v{VERSION}_servicelabels.parquet",
                DATA / f"pairs_usable_v{VERSION}_servicelabels.parquet", ROOT / "notebooks" / "04_proximity_counts.md"),
}
PRE_AXON = ["axon_fully_extended", "axon_partially_extended"]
AXON, DEND = 2, 3
STEP = 1000.0   # nm, max resampled edge length
PROX = 5000.0   # nm, axon-dendrite proximity radius
SYN_R = 3000.0  # nm, synapse-to-proximity radius


def load_skeleton(root_id):
    p = SKEL_DIR / f"{root_id}_sv{SKEL_VERSION}.npz"
    if not p.exists():
        return None
    z = np.load(p)
    sk = {k: z[k] for k in z.files}
    if LABELS == "synapse":
        sk["compartment"] = np.load(SKEL_DIR / f"{root_id}_sv{SKEL_VERSION}_synlabels.npz")["compartment"]
    return sk


def resample(sk, compartment):
    """Midpoints, lengths (nm) and soma path length (nm) of <= STEP pieces of one compartment."""
    v, e, comp = sk["vertices"], sk["edges"], sk["compartment"]
    child, parent = e[:, 0], e[:, 1]
    seg = np.linalg.norm(v[child] - v[parent], axis=1)
    n_v = len(v)
    g = coo_matrix((np.r_[seg, seg], (np.r_[child, parent], np.r_[parent, child])), shape=(n_v, n_v)).tocsr()
    path = dijkstra(g, directed=False, indices=int(sk["root"]))

    keep = (comp[child] == compartment) & (seg > 0)
    c, p, L = child[keep], parent[keep], seg[keep]
    n = np.maximum(np.ceil(L / STEP), 1).astype(int)
    rep = np.repeat(np.arange(len(L)), n)
    k = np.concatenate([np.arange(m) for m in n]) if len(n) else np.array([], int)
    frac = (k + 0.5) / n[rep]
    pts = v[p[rep]] + frac[:, None] * (v[c[rep]] - v[p[rep]])
    # path length grows from parent to child along the edge
    plen = path[p[rep]] + frac * (path[c[rep]] - path[p[rep]])
    return pts, L[rep] / n[rep], plen


def main():
    global LABELS
    if "--labels" in sys.argv:
        LABELS = sys.argv[sys.argv.index("--labels") + 1]
    syn_out, PAIRS_PATH, COUNTS_MD = OUTPUTS[LABELS]
    print(f"labels: {LABELS}")
    t0 = time.time()
    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    usable["root_id"] = usable["root_id"].astype("int64")
    pre_meta = usable[usable["axon_level"].isin(PRE_AXON)]
    syn = pd.read_parquet(SYN_PATH)
    syn = syn.drop(columns=[c for c in ["pre_axon_path_um", "pre_axon_dist_um", "in_proximity"] if c in syn])
    autapse = syn["pre_root_id"] == syn["post_root_id"]

    # Dendrites of all usable neurons -> one KD-tree
    d_pts, d_len, d_owner, missing = [], [], [], []
    for i, r in enumerate(usable["root_id"]):
        sk = load_skeleton(r)
        if sk is None:
            missing.append(r)
            continue
        pts, ln, _ = resample(sk, DEND)
        d_pts.append(pts), d_len.append(ln), d_owner.append(np.full(len(ln), i))
    d_pts, d_len, d_owner = np.vstack(d_pts), np.concatenate(d_len), np.concatenate(d_owner)
    dtree = cKDTree(d_pts)
    dend_mm = np.bincount(d_owner, weights=d_len, minlength=len(usable)) / 1e6
    print(f"dendrite pieces: {len(d_pts):,} from {len(usable) - len(missing)} skeletons "
          f"({len(missing)} missing); {time.time() - t0:.0f}s", flush=True)

    root_idx = pd.Series(np.arange(len(usable)), index=usable["root_id"])
    syn_cols = {"pre_axon_path_um": np.full(len(syn), np.nan), "pre_axon_dist_um": np.full(len(syn), np.nan),
                "in_proximity": np.zeros(len(syn), bool)}
    ld = np.full((len(pre_meta), len(usable)), np.nan)
    axon_mm, no_axon = np.full(len(pre_meta), np.nan), []

    for j, r in enumerate(pre_meta["root_id"]):
        sk = load_skeleton(r)
        if sk is None:
            continue
        a_pts, a_len, a_path = resample(sk, AXON)
        if len(a_pts) == 0:
            no_axon.append(r)
            continue
        axon_mm[j] = a_len.sum() / 1e6
        atree = cKDTree(a_pts)
        me = root_idx[r]

        # Proximity: axon point i -> dendrite pieces within PROX
        nbrs = atree.query_ball_tree(dtree, PROX)
        lens = np.fromiter(map(len, nbrs), int, len(nbrs))
        a_idx = np.repeat(np.arange(len(nbrs)), lens)
        d_idx = np.fromiter(itertools.chain.from_iterable(nbrs), int, lens.sum())
        other = d_owner[d_idx] != me
        a_idx, d_idx = a_idx[other], d_idx[other]
        uniq_d = np.unique(d_idx)
        row = np.bincount(d_owner[uniq_d], weights=d_len[uniq_d], minlength=len(usable)) / 1e3  # um
        row[me] = np.nan
        row[np.isin(usable["root_id"], missing)] = np.nan
        ld[j] = row

        # Synapses from this pre (non-autapse)
        s_mask = (syn["pre_root_id"] == r).to_numpy() & ~autapse.to_numpy()
        if not s_mask.any():
            continue
        s_rows = np.flatnonzero(s_mask)
        s_pts = syn.loc[s_mask, ["ctr_x_nm", "ctr_y_nm", "ctr_z_nm"]].to_numpy(float)
        dist, near = atree.query(s_pts)
        syn_cols["pre_axon_path_um"][s_rows] = a_path[near] / 1e3
        syn_cols["pre_axon_dist_um"][s_rows] = dist / 1e3
        posts = root_idx.reindex(syn.loc[s_mask, "post_root_id"]).to_numpy()
        pair_post = d_owner[d_idx]
        for pi in np.unique(posts[~np.isnan(posts)]).astype(int):
            sel = pair_post == pi
            if not sel.any():
                continue
            prox_pts = np.vstack([a_pts[np.unique(a_idx[sel])], d_pts[np.unique(d_idx[sel])]])
            k = posts == pi
            dd, _ = cKDTree(prox_pts).query(s_pts[k], distance_upper_bound=SYN_R)
            syn_cols["in_proximity"][s_rows[k]] = np.isfinite(dd)
        if (j + 1) % 50 == 0:
            print(f"  pre {j + 1}/{len(pre_meta)}: {time.time() - t0:.0f}s", flush=True)

    # %% Save synapses with the new columns
    for c, v in syn_cols.items():
        syn[c] = v
    syn["in_proximity"] = syn["in_proximity"].astype("boolean")
    syn.loc[autapse, "in_proximity"] = pd.NA
    syn.to_parquet(syn_out)

    # %% Pairs table
    pre_ids = pre_meta["root_id"].to_numpy()
    post_ids = usable["root_id"].to_numpy()
    pairs = pd.DataFrame({
        "pre_root_id": np.repeat(pre_ids, len(post_ids)),
        "post_root_id": np.tile(post_ids, len(pre_ids)),
        "Ld_um": ld.ravel(),
    })
    pairs = pairs[pairs["pre_root_id"] != pairs["post_root_id"]]
    s = syn[~autapse]
    nsyn = s.groupby(["pre_root_id", "post_root_id"]).agg(Nsyn=("synapse_id", "size"),
                                                         Nsyn_in_prox=("in_proximity", "sum"))
    pairs = pairs.join(nsyn, on=["pre_root_id", "post_root_id"]).fillna({"Nsyn": 0, "Nsyn_in_prox": 0})
    pairs[["Nsyn", "Nsyn_in_prox"]] = pairs[["Nsyn", "Nsyn_in_prox"]].astype(int)
    meta = usable.set_index("root_id")[["cell_type", "axon_level", "dendrite_level"]]
    pairs = pairs.join(meta[["cell_type", "axon_level"]].add_prefix("pre_"), on="pre_root_id")
    pairs = pairs.join(meta[["cell_type", "axon_level", "dendrite_level"]].add_prefix("post_"), on="post_root_id")
    pairs = pairs.reset_index(drop=True)
    pairs.attrs.update(version=str(VERSION), skeleton_version=str(SKEL_VERSION), labels=LABELS,
                       prox_um=str(PROX / 1e3), step_um=str(STEP / 1e3), syn_r_um=str(SYN_R / 1e3))
    pairs.to_parquet(PAIRS_PATH)
    runtime = time.time() - t0

    # %% Report
    ok = pairs["Ld_um"].notna()
    p = pairs[ok]
    ns = s.merge(p[["pre_root_id", "post_root_id"]], on=["pre_root_id", "post_root_id"])
    outside = 1 - ns["in_proximity"].mean()
    by = p.groupby(["pre_cell_type", "pre_axon_level"]).agg(
        pre_neurons=("pre_root_id", "nunique"), pairs=("Ld_um", "size"),
        pairs_Ld_gt0=("Ld_um", lambda x: (x > 0).sum()), Ld_mm=("Ld_um", lambda x: x.sum() / 1e3),
        Nsyn=("Nsyn", "sum"), Nsyn_in_prox=("Nsyn_in_prox", "sum"))
    by.loc[("**All**", ""), :] = by.sum()
    by["syn_per_mm_Ld"] = (by["Nsyn_in_prox"] / by["Ld_mm"]).round(3)
    by["all_syn_per_mm_Ld"] = (by["Nsyn"] / by["Ld_mm"]).round(3)
    by = by.rename(index={"axon_fully_extended": "fully_ext", "axon_partially_extended": "partially_ext"})
    by_ct = p.groupby("pre_cell_type").agg(Ld_mm=("Ld_um", lambda x: x.sum() / 1e3),
                                           Nsyn_in_prox=("Nsyn_in_prox", "sum"), Nsyn=("Nsyn", "sum"))
    by_ct["syn_per_mm_Ld"] = (by_ct["Nsyn_in_prox"] / by_ct["Ld_mm"]).round(3)

    ax = pd.DataFrame({"axon_mm": axon_mm, "axon_level": pre_meta["axon_level"].to_numpy(),
                       "cell_type": pre_meta["cell_type"].to_numpy()})
    ax_tab = ax.groupby(["cell_type", "axon_level"])["axon_mm"].median().unstack().round(2)
    path_q = s["pre_axon_path_um"].describe(percentiles=[.1, .5, .9]).round(1)
    dist_q = s["pre_axon_dist_um"].describe(percentiles=[.5, .9, .99]).round(2)
    no_axon_syn = s["pre_root_id"].isin(no_axon).sum()
    far = s["pre_axon_dist_um"] > SYN_R / 1e3
    far_out = (far & ~s["in_proximity"].fillna(False)).sum()

    fails = pd.read_csv(SKEL_DIR / "failures.csv") if (SKEL_DIR / "failures.csv").exists() else pd.DataFrame()
    lines = [
        f"# Proximity counts, {LABELS} labels (auto-generated by scripts/compute_proximity.py --labels {LABELS}, v{VERSION})\n",
        f"Skeleton version {SKEL_VERSION}; resampled to <= {STEP / 1e3:g} um; proximity radius "
        f"{PROX / 1e3:g} um; synapse-to-proximity radius {SYN_R / 1e3:g} um. Runtime {runtime / 60:.1f} min.\n",
        f"- Skeletons missing: {len(missing)} (fetch failures logged: {len(fails)}); "
        f"presynaptic skeletons with no axon vertices: {len(no_axon)} {no_axon}",
        f"- Pairs: {len(pairs):,} ({ok.sum():,} with Ld computed); **Ld > 0: {(p['Ld_um'] > 0).sum():,}**",
        f"- Non-autapse synapses in computed pairs: {len(ns):,}; **outside any proximity: "
        f"{(~ns['in_proximity'].astype(bool)).sum():,} ({outside:.1%})**",
        f"- Synapses from the {len(no_axon)} no-axon presynaptic skeletons (no path length, pairs Ld = NaN): {no_axon_syn}",
        f"- Synapses > {SYN_R / 1e3:g} um from the nearest labeled pre-axon point: {far.sum()} "
        f"({far.mean():.1%} of non-autapse); {far_out} of the {(~ns['in_proximity'].astype(bool)).sum()} "
        f"outside-proximity synapses are among them",
        f"- Median dendrite length per usable neuron: {np.median(dend_mm[dend_mm > 0]):.2f} mm\n",
        "## Synapse density by pre cell type\n",
        "syn_per_mm_Ld = synapses in proximity / summed Ld (mm); all_syn_per_mm_Ld uses all synapses.\n",
        by_ct.to_markdown(), "",
        "## By pre cell type x axon level\n", by.to_markdown(), "",
        "## Median resampled axon length (mm) per presynaptic neuron\n", ax_tab.to_markdown(), "",
        "## Synapse axon path length from soma (um), non-autapse\n", path_q.to_frame().T.to_markdown(index=False), "",
        "## Synapse distance to nearest resampled axon point (um)\n", dist_q.to_frame().T.to_markdown(index=False), "",
    ]
    COUNTS_MD.write_text("\n".join(lines))
    print("\n".join(lines))
    print(f"saved {PAIRS_PATH.relative_to(ROOT)}, {syn_out.relative_to(ROOT)}; {runtime:.0f}s")


if __name__ == "__main__":
    main()
