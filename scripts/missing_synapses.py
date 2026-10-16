"""Session 12: which synapses were still missing in the partial states? (Why did the orientation rule look weaker?)

Run from the repo root after scripts/edit_history.py build:
    .venv/bin/python scripts/missing_synapses.py [--perms 2000]

For the V1->V1 synapses of the 193 fully extended axons (the l2l synapse set), at each history state
(Apr 2023, Nov 2025, Feb 2026) a synapse is "missing" if it was not yet attached to its soma. If the missing
synapses were onto orientation-similar targets more often than the attached ones, losing them weakens the
orientation rule directly; if not, the weaker rule in the partial data is loss of power.

Statistic: mean similarity (orientation difference; RF distance as a positive control) of missing minus attached
synapses. Null: attachment permuted among each axon's own synapses, (a) within axon only, (b) within axon and
100 um path-length bin from the soma (so the null keeps where on the arbor synapses were missing), (c) within
axon and 50 um soma-distance bin of the target. p two-sided, centred on the null mean. Also by presynaptic cell
type (4P, 23P), and the cell-type composition of attached synapses per state.
Writes notebooks/12_missing_synapses.md.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from config import VERSION  # noqa: E402
from prepare_like_to_like import GOSI_MIN, MAX_DIST_UM, ori_diff_deg, rf_angle_deg  # noqa: E402

DATA = ROOT / "data"
STATES = {202304: "Apr 2023", 202511: "Nov 2025", 202602: "Feb 2026"}
NULLS = {"within axon": None, "within axon x path bin (100 um)": ("pre_axon_path_um", 100),
         "within axon x target-distance bin (50 um)": ("soma_dist_um", 50)}


def load():
    n = pd.read_parquet(DATA / f"l2l_neurons_v{VERSION}.parquet")
    n["root_id"] = n["root_id"].astype("int64")
    pairs = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet", columns=["pre_root_id"])
    pres = set(pairs["pre_root_id"].astype("int64"))
    s = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        s[c] = s[c].astype("int64")
    s = s[s["pre_root_id"].isin(pres) & s["post_root_id"].isin(n["root_id"]) & (s["pre_root_id"] != s["post_root_id"])
          & (s["pre_axon_dist_um"] <= MAX_DIST_UM)].copy()
    ni = n.set_index("root_id")
    pre, post = ni.loc[s["pre_root_id"]], ni.loc[s["post_root_id"]]
    s = s[(pre["area"].to_numpy() == "V1") & (post["area"].to_numpy() == "V1")]
    pre, post = ni.loc[s["pre_root_id"]], ni.loc[s["post_root_id"]]
    s["pre_cell_type"] = pre["cell_type"].to_numpy()
    s["rf_dist_deg"] = rf_angle_deg(pre["readout_loc_x"].to_numpy(), pre["readout_loc_y"].to_numpy(),
                                    post["readout_loc_x"].to_numpy(), post["readout_loc_y"].to_numpy())
    tuned = (pre["gOSI"].to_numpy(float) > GOSI_MIN) & (post["gOSI"].to_numpy(float) > GOSI_MIN)
    s["ori_diff_deg"] = np.where(tuned, ori_diff_deg(pre["pref_ori"].to_numpy(), post["pref_ori"].to_numpy()), np.nan)
    d = pre[["soma_x_um", "soma_y_um", "soma_z_um"]].to_numpy() - post[["soma_x_um", "soma_y_um", "soma_z_um"]].to_numpy()
    s["soma_dist_um"] = np.linalg.norm(d, axis=1)
    return s.reset_index(drop=True)


def perm_test(x, att, groups, n_perm, rng):
    """Mean x of missing minus attached; null by permuting att within groups."""
    obs = x[~att].mean() - x[att].mean()
    order = np.argsort(groups, kind="stable")
    g_sorted = groups[order]
    bounds = np.flatnonzero(np.r_[True, g_sorted[1:] != g_sorted[:-1], True])
    a_sorted = att[order].copy()
    xs = x[order]
    null = np.empty(n_perm)
    for k in range(n_perm):
        for lo, hi in zip(bounds[:-1], bounds[1:]):
            if hi - lo > 1:
                a_sorted[lo:hi] = rng.permutation(a_sorted[lo:hi])
        null[k] = xs[~a_sorted].mean() - xs[a_sorted].mean()
    mu = null.mean()
    p = (1 + np.sum(np.abs(null - mu) >= abs(obs - mu))) / (1 + n_perm)
    return obs, mu, np.quantile(null, [.025, .975]), p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=12)
    args = ap.parse_args()
    s = load()
    rows = []
    for v, label in STATES.items():
        att_all = pd.read_parquet(DATA / "history" / f"attached_{v}.parquet").set_index("synapse_id")["attached"]
        a = att_all.reindex(s["synapse_id"]).fillna(False).to_numpy(bool)
        for ct in ("all", "4P", "23P"):
            m_ct = np.ones(len(s), bool) if ct == "all" else (s["pre_cell_type"] == ct).to_numpy()
            for sim in ("ori_diff_deg", "rf_dist_deg"):
                m = m_ct & s[sim].notna().to_numpy()
                x = s.loc[m, sim].to_numpy(float)
                att = a[m]
                if att.all() or (~att).all():
                    continue
                for null_name, spec in NULLS.items():
                    g = s.loc[m, "pre_root_id"].astype(str)
                    if spec is not None:
                        col, w = spec
                        g = g + "_" + (s.loc[m, col] // w).astype(int).astype(str)
                    rng = np.random.default_rng(args.seed)
                    obs, mu, (lo, hi), p = perm_test(x, att, pd.factorize(g)[0], args.perms, rng)
                    rows.append(dict(state=label, cell_type=ct, similarity=sim.replace("_deg", "").replace("_", " "),
                                     null=null_name, n_syn=int(m.sum()), missing=int((~att).sum()),
                                     mean_missing=round(x[~att].mean(), 2), mean_attached=round(x[att].mean(), 2),
                                     diff=round(obs, 2), null_mean=round(mu, 2), null_95=f"[{lo:+.2f}, {hi:+.2f}]",
                                     diff_beyond_null=round(obs - mu, 2), p=round(p, 4)))
                    print(rows[-1], flush=True)
    t = pd.DataFrame(rows)

    # cell-type composition of attached V1->V1 synapses per state
    comp = []
    for v, label in {**STATES, VERSION: "complete (v1822)"}.items():
        att_all = pd.read_parquet(DATA / "history" / f"attached_{v}.parquet").set_index("synapse_id")["attached"]
        a = att_all.reindex(s["synapse_id"]).fillna(False).to_numpy(bool)
        for ct in ("4P", "23P"):
            m = (s["pre_cell_type"] == ct).to_numpy()
            comp.append(dict(state=label, cell_type=ct, attached=int(a[m].sum()), total=int(m.sum()),
                             frac_attached=round(a[m].mean(), 3),
                             share_of_attached=round(a[m].sum() / a.sum(), 3)))
    c = pd.DataFrame(comp)

    lines = [f"# Session 12: which synapses were missing in the partial states? (auto-generated, v{VERSION})\n",
             "V1->V1 synapses of the 193 fully extended axons. diff = mean similarity (degrees) of synapses not yet "
             "attached minus attached ones; negative = the missing synapses were onto MORE similar targets (their "
             "loss weakens the rule). Null: attachment permuted among each axon's synapses (within the stated "
             f"bins), {args.perms} permutations; p two-sided, centred on the null mean.\n",
             "## Orientation difference\n", t[t["similarity"] == "ori diff"].drop(columns="similarity").to_markdown(index=False), "",
             "## RF distance (positive control: missing synapses are distal, so expected more RF-dissimilar)\n",
             t[t["similarity"] == "rf dist"].drop(columns="similarity").to_markdown(index=False), "",
             "## Cell-type composition of attached V1->V1 synapses\n", c.to_markdown(index=False), ""]
    (ROOT / "notebooks" / "12_missing_synapses.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
