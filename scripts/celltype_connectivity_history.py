"""Session 14: does incompleteness bias a second, structural wiring rule -- which cell types an axon connects to?

Run from the repo root after scripts/edit_history.py build:
    .venv/bin/python scripts/celltype_connectivity_history.py [--boot 1000] [--reps 20]

Rule: for presynaptic cell type P, the log odds ratio of connecting (>= 1 synapse) to a target of type A rather than
reference type B, over all V1->V1 candidate pairs of the 193 complete axons (same pairs as the like-to-like
analysis). Two versions:
  raw         pooled 2x2 table
  distance    Mantel-Haenszel odds ratio across 50 um soma-distance bins (the analogue of the distance-matched
              shuffle: compares A and B targets at the same distance)
  depth       Mantel-Haenszel across 50 um bins of signed target depth relative to the presynaptic soma x 50 um
              lateral bins (compares A and B targets at the same position above/below and beside the soma)
States: complete (v1822), Apr 2023, Nov 2025, Feb 2026 (presynaptic rollback, hist sets), and count-matched random
loss (each axon keeps a random subset of its synapses of its Nov 2025 size; 20 reps). Shift CIs: paired bootstrap
over presynaptic axons. Partial extension did not follow axons to layer 5 and below (C. Schneider-Mizell), so
connections onto 5P targets are the expected casualty.
Estimability: a stratified estimate is reported only if at least MIN_STRATA strata contribute to both the MH
numerator (A connected, B unconnected) and denominator (A unconnected, B connected) in the complete state; otherwise
"n/e". Bootstrap samples where a log OR is undefined (0 or infinite odds ratio) are dropped and their share is
reported when non-zero.
Writes notebooks/14_celltype_connectivity.md and data/celltype_history/{shifts,attachment_by_depth}_v<VERSION>.csv.
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
STATES = {"complete": None, "Apr 2023": 202304, "Nov 2025": 202511, "Feb 2026": 202602}
CONTRASTS = [("23P", "5P-ET", "23P"), ("23P", "4P", "23P"),
             ("4P", "5P-ET", "4P"), ("4P", "5P-IT", "4P"), ("4P", "23P", "4P")]
BIN_UM = 50
MAX_DIST_UM = 3.0
MIN_STRATA = 3


def tables(pairs, conn, pre_ct, a, b, strat="soma_dist"):
    """Per presynaptic axon x target type (a/b) x stratum: connected and candidate pair counts.
    strat = "soma_dist" (50 um soma-distance bins) or "depth" (50 um bins of signed depth of the target soma
    relative to the presynaptic soma x 50 um bins of lateral distance)."""
    m = (pairs["pre_cell_type"] == pre_ct) & pairs["post_cell_type"].isin([a, b])
    if strat == "soma_dist":
        bins = (pairs.loc[m, "soma_dist_um"] // BIN_UM).astype(int).to_numpy()
    else:
        bins = ((pairs.loc[m, "dy_um"] // BIN_UM).astype(int) * 1000 + (pairs.loc[m, "dxz_um"] // BIN_UM).astype(int)).to_numpy()
    d = pd.DataFrame({"pre": pairs.loc[m, "pre_root_id"].to_numpy(),
                      "isA": (pairs.loc[m, "post_cell_type"] == a).to_numpy(),
                      "bin": bins,
                      "c": conn[m.to_numpy()].astype(int)})
    g = d.groupby(["pre", "isA", "bin"]).agg(c=("c", "sum"), n=("c", "size")).reset_index()
    return g


def log_or(g, w=None):
    """Raw and Mantel-Haenszel log odds ratio (A vs B); w = per-axon weights (bootstrap counts)."""
    if w is not None:
        g = g.assign(c=g["c"] * g["pre"].map(w).fillna(0), n=g["n"] * g["pre"].map(w).fillna(0))
    s = g.groupby(["isA", "bin"])[["c", "n"]].sum().unstack("isA").fillna(0)
    a, n_a = s[("c", True)], s[("n", True)]
    c, n_b = s[("c", False)], s[("n", False)]
    b_, d_ = n_a - a, n_b - c
    raw = np.log((a.sum() * d_.sum()) / (b_.sum() * c.sum()))
    tot = n_a + n_b
    ok = tot > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        mh = np.log((a * d_ / tot)[ok].sum() / (b_ * c / tot)[ok].sum())
    return raw, mh


def n_informative(g):
    """Strata contributing to the MH numerator and to the denominator (complete-state estimability check)."""
    s = g.groupby(["isA", "bin"])[["c", "n"]].sum().unstack("isA").fillna(0)
    a, n_a, c, n_b = s[("c", True)], s[("n", True)], s[("c", False)], s[("n", False)]
    return int((a * (n_b - c) > 0).sum()), int(((n_a - a) * c > 0).sum())


def ci(x):
    """95% interval over finite bootstrap values; share of undefined samples."""
    f = x[np.isfinite(x)]
    return np.quantile(f, [.025, .975]), 1 - len(f) / len(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--seed", type=int, default=14)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    pairs_all = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        pairs_all[c] = pairs_all[c].astype("int64")
    v1 = (pairs_all["proj"] == "V1->V1").to_numpy()
    pairs = pairs_all[v1].reset_index(drop=True)
    neu = pd.read_parquet(DATA / f"l2l_neurons_v{VERSION}.parquet")
    neu["root_id"] = neu["root_id"].astype("int64")
    ni = neu.set_index("root_id")
    pre_xyz = ni.loc[pairs["pre_root_id"], ["soma_x_um", "soma_y_um", "soma_z_um"]].to_numpy()
    post_xyz = ni.loc[pairs["post_root_id"], ["soma_x_um", "soma_y_um", "soma_z_um"]].to_numpy()
    pairs["dy_um"] = post_xyz[:, 1] - pre_xyz[:, 1]          # y increases with cortical depth
    pairs["dxz_um"] = np.hypot(post_xyz[:, 0] - pre_xyz[:, 0], post_xyz[:, 2] - pre_xyz[:, 2])

    nsyn = {"complete": pairs_all["Nsyn"].to_numpy()[v1]}
    for lab, v in STATES.items():
        if v is not None:
            nsyn[lab] = pd.read_parquet(SET_DIR / f"hist_{v}_0.parquet")["Nsyn"].to_numpy()[v1]
    assert (nsyn["complete"] == pd.read_parquet(SET_DIR / "hist_1822_0.parquet")["Nsyn"].to_numpy()[v1]).all()

    # count-matched random loss (Nov 2025 counts per axon), on the l2l synapse list
    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    syn = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        syn[c] = syn[c].astype("int64")
    syn = syn[(syn["pre_root_id"] != syn["post_root_id"]) & (syn["pre_axon_dist_um"] <= MAX_DIST_UM)
              & syn["pre_root_id"].isin(pairs_all["pre_root_id"].unique())
              & syn["post_root_id"].isin(usable["root_id"].astype("int64"))].reset_index(drop=True)
    pidx = pd.Series(np.arange(len(pairs_all)), index=pd.MultiIndex.from_frame(pairs_all[["pre_root_id", "post_root_id"]]))
    s_row = pidx.reindex(pd.MultiIndex.from_frame(syn[["pre_root_id", "post_root_id"]])).to_numpy().astype(int)
    att = pd.read_parquet(DATA / "history" / "attached_202511.parquet").set_index("synapse_id")["attached"]
    a_now = att.reindex(syn["synapse_id"]).fillna(False).to_numpy(bool)
    groups = syn.groupby("pre_root_id").indices
    k = {r: int(a_now[ix].sum()) for r, ix in groups.items()}
    # synapse attachment by depth relative to the presynaptic soma (V1->V1 synapses of the complete axons)
    v1_syn = v1[s_row]
    dy = syn["ctr_y_nm"].to_numpy() / 1000 - ni.loc[syn["pre_root_id"], "soma_y_um"].to_numpy()
    ct = ni.loc[syn["pre_root_id"], "cell_type"].to_numpy()
    band_lab = ["> 150 above", "75-150 above", "25-75 above", "within 25", "25-75 below", "75-150 below", "> 150 below"]
    band = pd.cut(dy, [-np.inf, -150, -75, -25, 25, 75, 150, np.inf], labels=band_lab)
    att_by = {lab: pd.read_parquet(DATA / "history" / f"attached_{v}.parquet").set_index("synapse_id")["attached"]
              .reindex(syn["synapse_id"]).fillna(False).to_numpy(bool) for lab, v in STATES.items() if v is not None}
    depth_rows = []
    for pre_ct in ("4P", "23P"):
        for bl in band_lab:
            msk = v1_syn & (ct == pre_ct) & (band == bl)
            depth_rows.append({"presynaptic": pre_ct, "synapse depth vs. presynaptic soma (um)": bl,
                               "synapses (complete)": int(msk.sum()),
                               **{f"attached {lab}": f"{att_by[lab][msk].mean():.2f}" for lab in att_by}})

    rand = []
    for _ in range(args.reps):
        keep = np.zeros(len(syn), bool)
        for r, ix in groups.items():
            keep[rng.choice(ix, k[r], replace=False)] = True
        rand.append(np.bincount(s_row[keep], minlength=len(pairs_all))[v1])

    pres = pairs["pre_root_id"].unique()
    boots = [pd.Series(np.bincount(rng.integers(0, len(pres), len(pres)), minlength=len(pres)), index=pres)
             for _ in range(args.boot)]
    rows, num = [], []
    for pre_ct, a, b in CONTRASTS:
        g = {lab: tables(pairs, n > 0, pre_ct, a, b) for lab, n in nsyn.items()}
        gd = {lab: tables(pairs, n > 0, pre_ct, a, b, "depth") for lab, n in nsyn.items()}
        est = {lab: (*log_or(g[lab]), log_or(gd[lab])[1]) for lab in g}
        bs = {lab: np.array([(*log_or(g[lab], w), log_or(gd[lab], w)[1]) for w in boots]) for lab in g}
        rr = np.array([(*log_or(tables(pairs, n > 0, pre_ct, a, b)), log_or(tables(pairs, n > 0, pre_ct, a, b, "depth"))[1])
                       for n in rand])
        inf_strata = {1: n_informative(g["complete"]), 2: n_informative(gd["complete"])}
        for j, kind in enumerate(("raw", "distance-stratified", "depth x lateral-stratified")):
            nn = inf_strata.get(j)
            row = {"presynaptic": pre_ct, "contrast": f"{a} vs {b}", "log OR": kind,
                   "informative strata (num, den)": "" if nn is None else f"{nn[0]}, {nn[1]}"}
            if nn is not None and min(nn) < MIN_STRATA:
                row["complete"] = "n/e"
                for lab in ("Apr 2023", "Nov 2025", "Feb 2026"):
                    row[f"{lab} shift"] = "n/e"
                row["random loss (Nov 2025 counts) shift"] = "n/e"
                rows.append(row)
                continue
            row["complete"] = f"{est['complete'][j]:+.2f}"
            for lab in ("Apr 2023", "Nov 2025", "Feb 2026"):
                q, und = ci(bs[lab][:, j] - bs["complete"][:, j])
                num.append({"presynaptic": pre_ct, "contrast": f"{a} vs {b}", "kind": kind, "state": lab,
                            "complete": est["complete"][j], "shift": est[lab][j] - est["complete"][j],
                            "lo": q[0], "hi": q[1], "undefined_frac": und})
                row[f"{lab} shift"] = (f"{est[lab][j] - est['complete'][j]:+.2f} [{q[0]:+.2f}, {q[1]:+.2f}]"
                                       + (f" ({und:.0%} undef.)" if und > 0 else ""))
            fr = rr[:, j][np.isfinite(rr[:, j])]
            num.append({"presynaptic": pre_ct, "contrast": f"{a} vs {b}", "kind": kind, "state": "random loss (Nov 2025 counts)",
                        "complete": est["complete"][j], "shift": fr.mean() - est["complete"][j],
                        "lo": np.quantile(fr, .025) - est["complete"][j], "hi": np.quantile(fr, .975) - est["complete"][j],
                        "undefined_frac": 1 - len(fr) / len(rr)})
            row["random loss (Nov 2025 counts) shift"] = (f"{fr.mean() - est['complete'][j]:+.2f} "
                                                          f"[{np.quantile(fr, .025) - est['complete'][j]:+.2f}, "
                                                          f"{np.quantile(fr, .975) - est['complete'][j]:+.2f}]"
                                                          + (f" ({1 - len(fr) / len(rr):.0%} undef.)" if len(fr) < len(rr) else ""))
            rows.append(row)
        print(*rows[-3:], sep="\n", flush=True)

    # connections kept per target type, Nov 2025
    kept = []
    for pre_ct in ("23P", "4P"):
        for post_ct in ("23P", "4P", "5P-IT", "5P-ET"):
            m = ((pairs["pre_cell_type"] == pre_ct) & (pairs["post_cell_type"] == post_ct)).to_numpy()
            full = (nsyn["complete"][m] > 0).sum()
            kept.append({"presynaptic": pre_ct, "target": post_ct, "connections (complete)": int(full),
                         **{f"kept {lab}": f"{(nsyn[lab][m] > 0).sum() / full:.2f}" for lab in ("Apr 2023", "Nov 2025", "Feb 2026")},
                         "median soma distance (um)": int(pairs.loc[m & (nsyn["complete"] > 0), "soma_dist_um"].median())})
    t, kt, dt = pd.DataFrame(rows), pd.DataFrame(kept), pd.DataFrame(depth_rows)
    out = DATA / "celltype_history"
    out.mkdir(exist_ok=True)
    pd.DataFrame(num).to_csv(out / f"shifts_v{VERSION}.csv", index=False)
    dt.to_csv(out / f"attachment_by_depth_v{VERSION}.csv", index=False)
    (ROOT / "notebooks" / "14_celltype_connectivity.md").write_text(
        f"# Session 14: cell-type connection preference across the proofreading history (auto-generated, v{VERSION})\n\n"
        "V1->V1 candidate pairs of the 193 complete axons. log OR = log odds of connecting to a target of type A "
        "rather than B; distance-stratified = Mantel-Haenszel over 50 um soma-distance bins; depth x lateral-stratified = "
        "Mantel-Haenszel over 50 um bins of the target's signed depth relative to the presynaptic soma x 50 um lateral bins. Shift = state minus "
        f"complete [95% paired bootstrap over presynaptic axons, {args.boot} samples]. Random loss = each axon keeps a "
        f"random subset of its synapses of its Nov 2025 size ({args.reps} reps; mean [2.5%, 97.5%]). "
        f"n/e = not estimable: fewer than {MIN_STRATA} strata contribute to the Mantel-Haenszel numerator or "
        "denominator in the complete state (the two target types barely overlap in space). "
        "Informative strata (num, den) = strata with A connected & B unconnected, and A unconnected & B connected.\n\n"
        "## Fraction of complete-state connections present, by target type\n\n" + kt.to_markdown(index=False) +
        "\n\n## Synapse attachment by depth relative to the presynaptic soma (V1->V1)\n\n"
        "Depth = synapse y minus presynaptic soma y (y increases with cortical depth). Fraction of the complete axons' "
        "V1->V1 synapses attached to the soma at each state.\n\n" + dt.to_markdown(index=False) +
        "\n\n## Target-type preference (log odds ratio)\n\n" + t.to_markdown(index=False) + "\n")
    print(kt.to_string(index=False))


if __name__ == "__main__":
    main()
