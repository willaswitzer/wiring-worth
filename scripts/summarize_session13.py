"""Session 13: like-to-like by presynaptic x postsynaptic cell type (complete axons, V1->V1, Nsyn).

Closest match to Ding et al. 2025's restriction of the orientation analysis to V1 L2/3 -> L2/3 (their Extended Data
Fig. 10) is 23P -> 23P. Run after scripts/run_session13.sh:
    .venv/bin/python scripts/summarize_session13.py
Writes notebooks/13_celltype_pairs.md.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import VERSION  # noqa: E402

TD = ROOT / "data" / "truncation"
COMBOS = [("23P", "23P"), ("4P", "4P"), ("4P", "23P"), ("23P", "4P")]
SIMS = {"rf_dist_deg": "RF", "ori_diff_deg": "orientation"}


def v1(f):
    if not f.exists():
        return None
    d = pd.read_csv(f)
    d = d[(d["proj"] == "V1->V1") & (d["model"] == "Nsyn")].copy()
    d["level"] = d["level"].round(2)
    return d[d["level"] == 1]


def gap(real, perm, sim):
    if perm is None:
        return ""
    v = perm[perm["sim"] == sim]["slope"].dropna().to_numpy()
    mu = v.mean()
    p = (1 + np.sum(np.abs(v - mu) >= abs(real - mu))) / (1 + len(v))
    return f"{real - mu:+.3f} (p = {p:.3f})"


def main():
    pairs = pd.read_parquet(ROOT / "data" / f"l2l_pairs_v{VERSION}.parquet")
    pairs = pairs[pairs["proj"] == "V1->V1"]
    rows = []
    for pre, post in COMBOS:
        sfx = f"_{pre}_to{post}"
        fit = v1(TD / f"trunc_slopes_poisson{sfx}.csv")
        if fit is None:
            continue
        for sim, lab in SIMS.items():
            r = fit[fit["sim"] == sim]
            if not len(r):
                continue
            m = (pairs["pre_cell_type"] == pre) & (pairs["post_cell_type"] == post) & pairs[sim].notna()
            rows.append({"pre → post": f"{pre} → {post}", "similarity": lab,
                         "synapses": int(pairs.loc[m, "Nsyn"].sum()), "axons": int(pairs.loc[m, "pre_root_id"].nunique()),
                         "slope [Wald 95%]": f"{r['slope'].iloc[0]:+.3f} [{r['slope'].iloc[0] - 1.96 * r['SE'].iloc[0]:+.3f}, "
                                             f"{r['slope'].iloc[0] + 1.96 * r['SE'].iloc[0]:+.3f}]",
                         "beyond distance-matched shuffle": gap(r["slope"].iloc[0], v1(TD / f"trunc_perm_poisson{sfx}_nulldistance.csv"), sim),
                         "beyond cell-type shuffle": gap(r["slope"].iloc[0], v1(TD / f"trunc_perm_poisson{sfx}.csv"), sim)})
    t = pd.DataFrame(rows).sort_values(["similarity", "pre → post"], ascending=[False, True])
    (ROOT / "notebooks" / "13_celltype_pairs.md").write_text(
        f"# Session 13: like-to-like by presynaptic x postsynaptic cell type (auto-generated, v{VERSION})\n\n"
        "Complete axons, V1->V1, Nsyn model, slope per 10°; 200 shuffles per null, p two-sided centred on the shuffled mean. "
        "Synapses and axons counted over pairs where the similarity is defined (orientation: both gOSI > 0.25).\n\n"
        + t.to_markdown(index=False) + "\n")
    print(t.to_string(index=False))


if __name__ == "__main__":
    main()
