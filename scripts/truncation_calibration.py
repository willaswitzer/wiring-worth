"""Session 5, task 2: are partially extended axons radially truncated or thinned across the arbor?

Run from the repo root after scripts/relabel_skeletons.py and scripts/compute_proximity.py:
    .venv/bin/python scripts/truncation_calibration.py

Within 23P and within 4P, fully vs partially extended axons (synapse-flow labels):
1. Synapse level: pre_axon_path_um of non-autapse synapses onto usable neurons (synapses more than
   3 um from the labeled axon excluded): median, IQR, two-sample KS test. Synapses from one neuron
   are not independent, so a neuron-level test (Mann-Whitney U on per-neuron median path length)
   is reported alongside.
2. Supplement: the same using ALL output synapses of each neuron (any postsynaptic partner).
3. Neuron level: axon cable length, number of axon endpoints, max axon path length; and axon
   cable length per path-distance bin. Radial truncation predicts a partial/full cable ratio
   that falls with distance; thinning predicts a flat ratio below 1.
Writes figures/05_truncation_calibration.png and notebooks/05_truncation_counts.md.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import dijkstra  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import VERSION  # noqa: E402

DATA = ROOT / "data"
SKEL_DIR = DATA / "skeletons"
FIG = ROOT / "figures" / "05_truncation_calibration.png"
COUNTS_MD = ROOT / "notebooks" / "05_truncation_counts.md"
CELL_TYPES = ["23P", "4P"]
LEVELS = {"axon_fully_extended": "fully extended", "axon_partially_extended": "partially extended"}
COLORS = {"axon_fully_extended": "#2a78d6", "axon_partially_extended": "#eb6834"}  # validated pair
MAX_DIST_UM = 3.0
BIN_UM = 100
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#898781", "#e1e0d9", "#fcfcfb"


def load(root_id):
    z = np.load(SKEL_DIR / f"{root_id}_sv4.npz")
    sk = {k: z[k] for k in z.files}
    sk["compartment"] = np.load(SKEL_DIR / f"{root_id}_sv4_synlabels.npz")["compartment"]
    return sk


def axon_morphology(sk):
    """Per-edge axon cable with path distance, axon tips, and vertex path lengths (all in um)."""
    v, e, comp = sk["vertices"], sk["edges"], sk["compartment"]
    child, parent = e[:, 0], e[:, 1]
    seg = np.linalg.norm(v[child] - v[parent], axis=1)
    g = coo_matrix((np.r_[seg, seg], (np.r_[child, parent], np.r_[parent, child])), shape=(len(v),) * 2)
    path = dijkstra(g.tocsr(), directed=False, indices=int(sk["root"])) / 1e3
    is_ax = comp[child] == 2
    n_children = np.bincount(parent, minlength=len(v))
    ax_v = np.flatnonzero(comp == 2)
    tips = ax_v[n_children[ax_v] == 0]
    mid = (path[child] + path[parent]) / 2
    return dict(cable_len=seg[is_ax] / 1e3, cable_mid=mid[is_ax], n_tips=len(tips),
                max_path=path[ax_v].max() if len(ax_v) else np.nan, v_path=path, ax_v=ax_v)


def describe(x):
    q1, med, q3 = np.percentile(x, [25, 50, 75]) if len(x) else (np.nan,) * 3
    return med, q1, q3


def compare(df, value, group_col="pre_root_id"):
    """Synapse-level KS + neuron-level Mann-Whitney on per-neuron medians, full vs partial."""
    full = df[df["axon_level"] == "axon_fully_extended"]
    part = df[df["axon_level"] == "axon_partially_extended"]
    ks = stats.ks_2samp(full[value], part[value])
    nf, np_ = full.groupby(group_col)[value].median(), part.groupby(group_col)[value].median()
    mw = stats.mannwhitneyu(nf, np_)
    row = {}
    for lab, d, n in [("full", full, nf), ("partial", part, np_)]:
        med, q1, q3 = describe(d[value])
        row.update({f"{lab}_n_syn": len(d), f"{lab}_n_neurons": d[group_col].nunique(),
                    f"{lab}_median": med, f"{lab}_IQR": f"{q1:.0f}-{q3:.0f}"})
    row.update(KS_D=ks.statistic, KS_p=ks.pvalue, neuron_MW_p=mw.pvalue)
    return row


def main():
    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    pre = usable[usable["cell_type"].isin(CELL_TYPES) & usable["axon_level"].isin(LEVELS)].copy()
    pre["root_id"] = pre["root_id"].astype("int64")

    # %% Neuron-level morphology + all-output synapse path lengths
    io = pd.read_parquet(DATA / f"synapses_io_usable_v{VERSION}.parquet")
    io = io[(io["direction"] == "out") & io["root_id"].isin(pre["root_id"])]
    io_groups = dict(tuple(io.groupby("root_id")))
    morph, cable_rows, all_out = [], [], []
    for r in pre["root_id"]:
        sk = load(r)
        m = axon_morphology(sk)
        if len(m["ax_v"]) == 0:
            continue
        morph.append(dict(root_id=r, axon_cable_mm=m["cable_len"].sum() / 1e3, n_axon_tips=m["n_tips"],
                          max_axon_path_um=m["max_path"]))
        b = (m["cable_mid"] // BIN_UM).astype(int)
        cable_rows.append(pd.DataFrame({"root_id": r, "bin": b, "cable_um": m["cable_len"]})
                          .groupby(["root_id", "bin"], as_index=False)["cable_um"].sum())
        g = io_groups.get(r)
        if g is not None:
            pts = g[["ctr_x_nm", "ctr_y_nm", "ctr_z_nm"]].to_numpy(float)
            d, i = cKDTree(sk["vertices"][m["ax_v"]]).query(pts)
            keep = d / 1e3 <= MAX_DIST_UM
            all_out.append(pd.DataFrame({"pre_root_id": r, "path_um": m["v_path"][m["ax_v"][i[keep]]]}))
    morph = pd.DataFrame(morph).merge(pre[["root_id", "cell_type", "axon_level"]], on="root_id")
    all_out = pd.concat(all_out).merge(pre[["root_id", "cell_type", "axon_level"]],
                                       left_on="pre_root_id", right_on="root_id")
    cable = pd.concat(cable_rows).merge(pre[["root_id", "cell_type", "axon_level"]], on="root_id")

    # %% Synapses onto usable neurons (primary)
    syn = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    syn = syn[(syn["pre_root_id"] != syn["post_root_id"]) & (syn["pre_axon_dist_um"] <= MAX_DIST_UM)]
    syn = syn.merge(pre[["root_id", "cell_type", "axon_level"]], left_on="pre_root_id", right_on="root_id")

    lines = [f"# Session 5 truncation calibration (auto-generated by scripts/truncation_calibration.py, v{VERSION})\n",
             f"Synapse-flow axon labels. Synapses > {MAX_DIST_UM:g} um from the labeled axon excluded. "
             "KS = two-sample Kolmogorov-Smirnov on synapses (not independent within a neuron); "
             "neuron_MW_p = Mann-Whitney U on per-neuron medians.\n"]
    for title, df, val in [("Synapses onto usable neurons: pre_axon_path_um", syn, "pre_axon_path_um"),
                           ("Supplement: all output synapses, axon path length (um)", all_out, "path_um")]:
        rows = {ct: compare(df[df["cell_type"] == ct], val) for ct in CELL_TYPES}
        lines += [f"## {title}\n", pd.DataFrame(rows).T.to_markdown(floatfmt=".3g"), ""]

    mrows = []
    for ct in CELL_TYPES:
        for col in ["axon_cable_mm", "n_axon_tips", "max_axon_path_um"]:
            d = morph[morph["cell_type"] == ct]
            f, p = (d.loc[d["axon_level"] == lv, col] for lv in LEVELS)
            mrows.append({"cell_type": ct, "measure": col, "full_n": len(f),
                          "full_median": f.median(), "full_IQR": "{:.1f}-{:.1f}".format(*np.percentile(f, [25, 75])),
                          "partial_n": len(p), "partial_median": p.median(),
                          "partial_IQR": "{:.1f}-{:.1f}".format(*np.percentile(p, [25, 75])),
                          "partial/full": p.median() / f.median(), "MW_p": stats.mannwhitneyu(f, p).pvalue})
    lines += ["## Neuron-level axon morphology\n", pd.DataFrame(mrows).to_markdown(index=False, floatfmt=".3g"), ""]

    # Mean cable per neuron by path-distance bin (neurons with no cable in a bin count as 0)
    n_neur = morph.groupby(["cell_type", "axon_level"]).size()
    tot =cable.groupby(["cell_type", "axon_level", "bin"])["cable_um"].sum().reset_index()
    tot["mean_cable_mm"] = tot.apply(lambda r: r["cable_um"] / n_neur[(r["cell_type"], r["axon_level"])] / 1e3, axis=1)
    prof = tot.pivot_table(index=["cell_type", "bin"], columns="axon_level", values="mean_cable_mm", fill_value=0)
    prof["partial/full"] = prof["axon_partially_extended"] / prof["axon_fully_extended"]
    coarse = prof.reset_index()
    coarse["range_um"] = pd.cut(coarse["bin"] * BIN_UM, [0, 200, 400, 600, 800, 1e9], right=False,
                                labels=["0-200", "200-400", "400-600", "600-800", "800+"])
    c = coarse.groupby(["cell_type", "range_um"], observed=True)[list(LEVELS)].sum()
    c["partial/full"] = c["axon_partially_extended"] / c["axon_fully_extended"]
    lines += ["## Mean axon cable per neuron (mm) by path distance from soma\n",
              "Ratio falling with distance = radial truncation; flat ratio < 1 = thinning.\n",
              c.round(3).to_markdown(), ""]
    COUNTS_MD.write_text("\n".join(lines))
    print("\n".join(lines))

    # %% Figure: rows = cell type; cols = synapse path histogram (all outputs) | cable profile
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.2), facecolor=SURFACE, constrained_layout=True)
    bins = np.arange(0, 1600 + 50, 50)
    for i, ct in enumerate(CELL_TYPES):
        ax = axes[i, 0]
        for lv, lab in LEVELS.items():
            x = all_out.loc[(all_out["cell_type"] == ct) & (all_out["axon_level"] == lv), "path_um"]
            n = all_out.loc[(all_out["cell_type"] == ct) & (all_out["axon_level"] == lv), "pre_root_id"].nunique()
            ax.hist(x, bins=bins, density=True, histtype="step", lw=2, color=COLORS[lv],
                    label=f"{lab} (n = {n} neurons)")
            ax.axvline(np.median(x), color=COLORS[lv], lw=1, ls=":")
        ax.set_title(f"{ct}: output synapses along the axon", loc="left", color=INK, fontsize=10)
        ax.set_xlabel("Path length from soma (µm)")
        ax.set_ylabel("Density")
        ax.legend(frameon=False, fontsize=8)
        ax = axes[i, 1]
        p = prof.loc[ct]
        for lv, lab in LEVELS.items():
            ax.plot(p.index * BIN_UM + BIN_UM / 2, p[lv], lw=2, color=COLORS[lv], label=lab, marker="o", ms=4)
        ax.set_title(f"{ct}: mean axon cable per neuron", loc="left", color=INK, fontsize=10)
        ax.set_xlabel("Path length from soma (µm, 100 µm bins)")
        ax.set_ylabel("Cable (mm)")
        ax.set_xlim(0, 1600)
        ax.legend(frameon=False, fontsize=8)
        for a in axes[i]:
            a.set_facecolor(SURFACE)
            a.grid(axis="y", color=GRID, lw=0.6)
            a.set_axisbelow(True)
    fig.savefig(FIG, dpi=160, facecolor=SURFACE)
    print(f"saved {FIG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
