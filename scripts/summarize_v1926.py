"""v1926 add-on: complete-axon results on the newest release next to the pinned v1822 results.

Run from the repo root after scripts/run_v1926.sh (no MICRONS_VERSION needed):
    .venv/bin/python scripts/summarize_v1926.py [--new 1926]
V1->V1, Nsyn model, slope per 10°, complete (fully extended) axons only. For each release: number of complete axons,
the slope, and the gap beyond the cell-type and distance-matched shuffles (p two-sided, centred on the shuffled
mean), pooled, by presynaptic cell type, and by presynaptic x postsynaptic cell type.
Writes notebooks/15_v1926.md.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import PINNED  # noqa: E402

SIMS = {"rf_dist_deg": "RF", "ori_diff_deg": "orientation"}
GROUPS = [("all", "", ""), ("4P axons", "4P", ""), ("23P axons", "23P", ""), ("23P → 23P", "23P", "23P"),
          ("4P → 4P", "4P", "4P"), ("4P → 23P", "4P", "23P"), ("23P → 4P", "23P", "4P")]


def v1(f, level=1.0):
    if not f.exists():
        return None
    d = pd.read_csv(f)
    d = d[(d["proj"] == "V1->V1") & (d["model"] == "Nsyn")].copy()
    if "level" in d:
        d = d[d["level"].round(2) == level]
    return d


def gap(real, perm, sim):
    if perm is None or real is None:
        return ""
    v = perm[perm["sim"] == sim]["slope"].dropna().to_numpy()
    if not len(v):
        return ""
    mu = v.mean()
    p = (1 + np.sum(np.abs(v - mu) >= abs(real - mu))) / (1 + len(v))
    return f"{real - mu:+.3f} (p = {p:.3f}, n = {len(v)})"


def release(ver):
    sfx = "" if ver == PINNED else f"_v{ver}"
    ll, td = ROOT / "data" / f"like_to_like{sfx}", ROOT / "data" / f"truncation{sfx}"
    pairs_f = ROOT / "data" / f"l2l_pairs_v{ver}.parquet"
    if not pairs_f.exists():
        return []
    pairs = pd.read_parquet(pairs_f)
    pairs = pairs[pairs["proj"] == "V1->V1"]
    rows = []
    for label, pre, post in GROUPS:
        tag = (f"_{pre}" if pre else "") + (f"_to{post}" if post else "")
        if tag:
            fit = v1(td / f"trunc_slopes_poisson{tag}.csv")
            cell = v1(td / f"trunc_perm_poisson{tag}.csv")
        else:  # pooled: the baseline fit and its cell-type shuffle (like_to_like.R)
            fit = ll / "real_slopes.csv"
            fit = pd.read_csv(fit) if fit.exists() else None
            if fit is not None:
                fit = fit[(fit["variant"] == "main") & (fit["model"] == "Nsyn") & (fit["proj"] == "V1->V1")]
            cell = v1(ll / "perm_slopes_main_poisson.csv")
        dist = v1(td / f"trunc_perm_poisson{tag}_nulldistance.csv")
        m = np.ones(len(pairs), bool)
        if pre:
            m &= (pairs["pre_cell_type"] == pre).to_numpy()
        if post:
            m &= (pairs["post_cell_type"] == post).to_numpy()
        for sim, sl in SIMS.items():
            mm = m & pairs[sim].notna().to_numpy() & (pairs["Nsyn"] > 0).to_numpy()
            real = None
            if fit is not None and len(fit[fit["sim"] == sim]):
                real = float(fit[fit["sim"] == sim]["slope"].iloc[0])
            rows.append({"group": label, "similarity": sl, "release": f"v{ver}",
                         "complete axons with synapses": int(pairs.loc[mm, "pre_root_id"].nunique()),
                         "synapses": int(pairs.loc[mm, "Nsyn"].sum()),
                         "slope": "" if real is None else f"{real:+.3f}",
                         "beyond cell-type shuffle": gap(real, cell, sim),
                         "beyond distance-matched shuffle": gap(real, dist, sim)})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", type=int, default=1926)
    args = ap.parse_args()
    t = pd.DataFrame(release(PINNED) + release(args.new))
    order = {g[0]: i for i, g in enumerate(GROUPS)}
    t = t.sort_values(["similarity", "group", "release"], key=lambda c: c.map(order) if c.name == "group" else c,
                      ascending=[False, True, True])
    (ROOT / "notebooks" / "15_v1926.md").write_text(
        f"# v{args.new} add-on: complete-axon results on the newest release (auto-generated)\n\n"
        "V1->V1, Nsyn model, slope per 10°, complete axons only. Axons and synapses counted over pairs where the similarity is "
        "defined (orientation: both gOSI > 0.25). Shuffles: 200 per null unless n says otherwise; p two-sided, centred "
        "on the shuffled mean.\n\n" + t.to_markdown(index=False) + "\n")
    print(t.to_string(index=False))


if __name__ == "__main__":
    main()
