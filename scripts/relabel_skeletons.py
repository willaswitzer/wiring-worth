"""Session 5: relabel skeleton axon/dendrite compartments from synapse locations.

Run from the repo root after scripts/fetch_skeletons.py and scripts/fetch_io_synapses.py:
    .venv/bin/python scripts/relabel_skeletons.py

Method: synapse-flow split (Schneider-Mizell et al. 2016, as in pcg_skel/meshparty), repeated.
Every output (presynaptic) and input (postsynaptic) synapse of the neuron, autapses excluded,
is mapped to its nearest skeleton vertex. For each edge, flow = in_below * out_above +
out_below * in_above ("below" = the subtree away from the soma).
- Cut 1: the axon is the subtree below the max-flow edge, accepted if that subtree has a higher
  output fraction than the rest of the neuron.
- Extra cuts (up to MAX_CUTS - 1): the flow split is repeated on the remaining non-axon part; a
  cut is the max-flow edge among subtrees that are output-majority with >= MIN_EXTRA_OUT outputs. This
  recovers axon parts that form a second subtree (e.g. attached separately at the soma) without
  relabeling dendrites that carry scattered false-positive output detections.
Split quality is the segregation index SI = 1 - H(split) / H(unsplit) of the in/out mix;
SI near 1 = outputs and inputs cleanly separated.

Neurons whose split is unreliable (no accepted cut, or SI < MIN_SI) keep the service labels.
Writes data/skeletons/{root}_sv4_synlabels.npz (compartment: 1 soma, 2 axon, 3 dendrite) and
data/skeletons/relabel_summary.csv. The service labels in {root}_sv4.npz are left unchanged.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import PINNED, VERSION, VSFX  # noqa: E402

DATA = ROOT / "data"
SKEL_DIR = DATA / "skeletons"
SKEL_VERSION = 4
IO_PATH = DATA / f"synapses_io_usable_v{VERSION}.parquet"
SUMMARY = SKEL_DIR / f"relabel_summary{VSFX}.csv"
MAX_SYN_DIST = 10_000.0  # nm; synapses farther than this from the skeleton are ignored
MIN_SI = 0.5             # below this the split is flagged as unreliable
MAX_CUTS = 5             # first cut + up to 4 extra axon subtrees
MIN_EXTRA_OUT = 5        # outputs needed to accept an extra axon subtree
XYZ = ["ctr_x_nm", "ctr_y_nm", "ctr_z_nm"]


def labels_path(root_id):
    return SKEL_DIR / f"{root_id}_sv{SKEL_VERSION}_synlabels.npz"


def entropy(p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -(p * np.log(p) + (1 - p) * np.log(1 - p))


def split_axon(sk, out_pts, in_pts):
    v, e = sk["vertices"], sk["edges"]
    root = int(sk["root"])
    n = len(v)
    parent = np.full(n, -1)
    parent[e[:, 0]] = e[:, 1]
    children = [[] for _ in range(n)]
    for c, p in e:
        children[p].append(c)
    order, stack = [], [root]  # soma outward
    while stack:
        x = stack.pop()
        order.append(x)
        stack.extend(children[x])
    order = np.array(order)

    tree = cKDTree(v)
    counts = {}
    for name, pts in [("out", out_pts), ("in", in_pts)]:
        d, i = tree.query(pts) if len(pts) else (np.array([]), np.array([], int))
        counts[name] = np.bincount(i[d <= MAX_SYN_DIST], minlength=n).astype(float)
    t_out, t_in = counts["out"].sum(), counts["in"].sum()

    def subtree(x):
        sub, stack = np.zeros(n, bool), [x]
        while stack:
            y = stack.pop()
            sub[y] = True
            stack.extend(children[y])
        return sub

    def below_counts(mask):
        """Subtree output/input counts using only synapses on vertices in mask."""
        bo, bi = np.where(mask, counts["out"], 0.0), np.where(mask, counts["in"], 0.0)
        for x in order[::-1]:
            if parent[x] >= 0:
                bo[parent[x]] += bo[x]
                bi[parent[x]] += bi[x]
        return bo, bi

    axon = np.zeros(n, bool)
    n_cuts = 0
    for k in range(MAX_CUTS):
        bo, bi = below_counts(~axon)
        r_o, r_i = bo[root], bi[root]
        flow = bi * (r_o - bo) + bo * (r_i - bi)
        if k > 0:  # extra cuts: only subtrees that are output-majority with enough outputs
            flow[(bo < MIN_EXTRA_OUT) | (bo < bi)] = -1
        flow[root] = -1
        flow[axon] = -1
        if flow.max() < 0:
            break
        cut = int(np.argmax(flow))
        s_o, s_i = bo[cut], bi[cut]
        frac_s = s_o / max(s_o + s_i, 1)
        frac_r = (r_o - s_o) / max(r_o - s_o + r_i - s_i, 1)
        accept = (frac_s > frac_r and s_o > 0) if k == 0 else (frac_s >= 0.5 and s_o >= MIN_EXTRA_OUT)
        if not accept:
            break
        axon |= subtree(cut)
        n_cuts += 1
    axon[root] = False

    a_out, a_in = counts["out"][axon].sum(), counts["in"][axon].sum()
    r_out, r_in = t_out - a_out, t_in - a_in
    frac_a = a_out / max(a_out + a_in, 1)
    frac_r = r_out / max(r_out + r_in, 1)
    tot = t_out + t_in
    h_split = ((a_out + a_in) * entropy(frac_a) + (r_out + r_in) * entropy(frac_r)) / max(tot, 1)
    si = 1 - h_split / entropy(t_out / max(tot, 1)) if t_out > 0 and t_in > 0 else np.nan

    comp = np.full(n, 3, dtype=np.int8)
    comp[axon] = 2
    comp[root] = 1
    info = dict(n_out=int(t_out), n_in=int(t_in), n_out_far=int(len(out_pts) - t_out),
                axon_out=int(a_out), axon_in=int(a_in), rest_out=int(r_out), rest_in=int(r_in),
                axon_out_frac=frac_a, rest_out_frac=frac_r, si=si, n_cuts=n_cuts, split_ok=n_cuts > 0)
    return comp, info


def main():
    io = pd.read_parquet(IO_PATH)
    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    groups = dict(tuple(io.groupby("root_id")))
    rows = []
    roots = usable["root_id"].astype("int64")
    if VERSION != PINNED:
        # root ids are immutable: a root already labelled for the pinned release is the same object with the same
        # synapses, so keep its labels (and leave the pinned release's files untouched); label only new roots
        roots = [r for r in roots if not labels_path(r).exists()]
        print(f"v{VERSION}: {len(roots)} new roots to label")
    for r in roots:
        z = np.load(SKEL_DIR / f"{r}_sv{SKEL_VERSION}.npz")
        sk = {k: z[k] for k in z.files}
        g = groups.get(r, io.iloc[:0])
        out_pts = g.loc[g["direction"] == "out", XYZ].to_numpy(float)
        in_pts = g.loc[g["direction"] == "in", XYZ].to_numpy(float)
        comp, info = split_axon(sk, out_pts, in_pts)
        old = sk["compartment"]
        # Unreliable split (no accepted cut, or SI < MIN_SI): keep the service labels.
        reliable = info["split_ok"] and info["si"] >= MIN_SI
        if not reliable:
            comp = old.astype(np.int8)
        info.update(reliable=reliable, label_source="synapse_flow" if reliable else "service")
        np.savez_compressed(labels_path(r), compartment=comp, label_source=info["label_source"])

        e = sk["edges"]
        seg = np.linalg.norm(sk["vertices"][e[:, 0]] - sk["vertices"][e[:, 1]], axis=1)
        changed = old[e[:, 0]] != comp[e[:, 0]]
        info.update(
            root_id=r,
            axon_mm_service=seg[old[e[:, 0]] == 2].sum() / 1e6,
            axon_mm_syn=seg[comp[e[:, 0]] == 2].sum() / 1e6,
            relabeled_mm=seg[changed].sum() / 1e6,
            dend_to_axon_mm=seg[changed & (comp[e[:, 0]] == 2)].sum() / 1e6,
            axon_to_dend_mm=seg[changed & (comp[e[:, 0]] == 3)].sum() / 1e6,
        )
        rows.append(info)
    summ = pd.DataFrame(rows).merge(usable[["root_id", "cell_type", "axon_level"]], on="root_id")
    summ.to_csv(SUMMARY, index=False)
    print(f"relabeled {len(summ)} skeletons -> {SUMMARY.relative_to(ROOT)}")
    print(summ[["si", "n_cuts", "axon_out_frac", "rest_out_frac", "relabeled_mm", "dend_to_axon_mm",
                "axon_to_dend_mm"]].describe(percentiles=[.05, .5, .95]).round(3).to_string())
    print("n_cuts:", summ["n_cuts"].value_counts().sort_index().to_dict())
    print("unreliable (service labels kept):", (~summ["reliable"]).sum(), "by axon level:",
          summ.loc[~summ["reliable"], "axon_level"].value_counts().to_dict())
    return summ


if __name__ == "__main__":
    main()
