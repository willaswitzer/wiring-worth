"""Session 7: pair table for partially extended presynaptic axons (the plan's third curve).

Run from the repo root after scripts/prepare_like_to_like.py:
    .venv/bin/python scripts/prepare_partial.py

Same construction as prepare_like_to_like.py (Nsyn = non-autapse synapses with pre_axon_dist_um <= 3;
RF-centre distance and orientation difference from l2l_neurons; projection type from area), but the
presynaptic neurons are the usable neurons with *partially* extended axons. Rows for the fully extended
axons are appended unchanged from l2l_pairs, so one table holds both groups (column pre_axon_level).
Writes data/l2l_pairs_fullpartial_v1822.parquet.
"""
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


def main():
    neurons = pd.read_parquet(DATA / f"l2l_neurons_v{VERSION}.parquet")
    neurons["root_id"] = neurons["root_id"].astype("int64")
    pairs = pd.read_parquet(DATA / f"pairs_usable_v{VERSION}.parquet")
    pairs = pairs[pairs["pre_axon_level"] == "axon_partially_extended"].copy()
    for c in ("pre_root_id", "post_root_id"):
        pairs[c] = pairs[c].astype("int64")
    syn = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    for c in ("pre_root_id", "post_root_id"):
        syn[c] = syn[c].astype("int64")
    syn = syn[(syn["pre_root_id"] != syn["post_root_id"]) & (syn["pre_axon_dist_um"] <= MAX_DIST_UM)]
    nsyn = syn.groupby(["pre_root_id", "post_root_id"]).size().rename("Nsyn")
    pairs = pairs.drop(columns=["Nsyn", "Nsyn_in_prox"]).join(nsyn, on=["pre_root_id", "post_root_id"])
    pairs["Nsyn"] = pairs["Nsyn"].fillna(0).astype(int)
    pairs = pairs[pairs["Ld_um"].notna()]

    n = neurons.set_index("root_id")
    pre, post = n.loc[pairs["pre_root_id"]], n.loc[pairs["post_root_id"]]
    out = pd.DataFrame({
        "pre_root_id": pairs["pre_root_id"].to_numpy(), "post_root_id": pairs["post_root_id"].to_numpy(),
        "pre_cell_type": pre["cell_type"].to_numpy(), "post_cell_type": post["cell_type"].to_numpy(),
        "pre_area": pre["area"].to_numpy(), "post_area": post["area"].to_numpy(),
        "post_dendrite_level": post["dendrite_level"].to_numpy(),
        "Ld_um": pairs["Ld_um"].to_numpy(), "Nsyn": pairs["Nsyn"].to_numpy(),
    })
    proj = {("V1", "V1"): "V1->V1", ("HVA", "HVA"): "HVA->HVA", ("V1", "HVA"): "V1->HVA (FF)",
            ("HVA", "V1"): "HVA->V1 (FB)"}
    out["proj"] = [proj[k] for k in zip(out["pre_area"], out["post_area"])]
    d = pre[["soma_x_um", "soma_y_um", "soma_z_um"]].to_numpy() - post[["soma_x_um", "soma_y_um", "soma_z_um"]].to_numpy()
    out["soma_dist_um"] = np.linalg.norm(d, axis=1)
    out["rf_dist_deg"] = rf_angle_deg(pre["readout_loc_x"], pre["readout_loc_y"], post["readout_loc_x"], post["readout_loc_y"])
    tuned = (pre["gOSI"].to_numpy(float) > GOSI_MIN) & (post["gOSI"].to_numpy(float) > GOSI_MIN)
    out["ori_diff_deg"] = np.where(tuned, ori_diff_deg(pre["pref_ori"], post["pref_ori"]), np.nan)
    out["pre_axon_level"] = "partial"

    full = pd.read_parquet(DATA / f"l2l_pairs_v{VERSION}.parquet")
    full["pre_axon_level"] = "full"
    allp = pd.concat([full, out[full.columns]], ignore_index=True)
    allp.to_parquet(DATA / f"l2l_pairs_fullpartial_v{VERSION}.parquet")
    s = allp[allp["proj"] == "V1->V1"].groupby(["pre_cell_type", "pre_axon_level"]).agg(
        axons=("pre_root_id", "nunique"), synapses=("Nsyn", "sum"), Ld_mm=("Ld_um", lambda x: x.sum() / 1e3))
    s["syn_per_axon"] = (s["synapses"] / s["axons"]).round(1)
    s["Ld_mm_per_axon"] = (s["Ld_mm"] / s["axons"]).round(2)
    print(s.round(1).to_string())


if __name__ == "__main__":
    main()
