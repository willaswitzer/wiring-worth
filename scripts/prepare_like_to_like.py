"""Session 6: build the pair and neuron tables for the like-to-like models (fit in R).

Run from the repo root after scripts/compute_proximity.py (synapse labels):
    .venv/bin/python scripts/prepare_like_to_like.py

Pairs: presynaptic = usable neurons with fully extended axons (reference set); postsynaptic = all
other usable neurons. Nsyn counts non-autapse synapses with pre_axon_dist_um <= 3 (Session 5).
Similarity (computed here for the real data, and again in R for each permutation):
- RF-centre distance (deg): angle at the eye between the two readout locations on the monitor.
  Monitor 31.8 x 56.5 cm (height x width), 15 cm from the eye (MICrONS Consortium 2025, Methods,
  "Monitor positioning"), assumed centred on and normal to the eye (Ding et al. 2025, Methods).
  readout_loc_x/y span [-1, 1] across the monitor width/height ((-1,-1) top-left).
- Orientation difference (deg): |pref_ori_pre - pref_ori_post| folded to 0-90; defined only
  when both neurons have gOSI > 0.25.
Projection type from nucleus_functional_area_assignment (V1 vs HVA = AL/RL/LM).
Writes data/l2l_pairs_v1822.parquet and data/l2l_neurons_v1822.parquet.
Optional argument TAG (Session 10 match-rule sensitivity): reads usable_set_v1822_s2_TAG.parquet and writes
l2l_pairs_v1822_TAG.parquet / l2l_neurons_v1822_TAG.parquet; pairs are limited to neurons in that set.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import VERSION  # noqa: E402

DATA = ROOT / "data"
MONITOR_W_CM, MONITOR_H_CM, EYE_CM = 56.5, 31.8, 15.0
GOSI_MIN = 0.25
MAX_DIST_UM = 3.0


def rf_vectors(x, y):
    """Unit vectors from the eye to readout locations (x, y in [-1, 1])."""
    v = np.stack([np.asarray(x) * MONITOR_W_CM / 2, np.asarray(y) * MONITOR_H_CM / 2,
                  np.full(len(x), EYE_CM)], axis=1)
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def rf_angle_deg(x1, y1, x2, y2):
    c = np.einsum("ij,ij->i", rf_vectors(x1, y1), rf_vectors(x2, y2))
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def ori_diff_deg(a, b):
    d = np.abs(np.asarray(a) - np.asarray(b)) % 180
    return np.minimum(d, 180 - d)


def main(tag=""):
    sfx = f"_{tag}" if tag else ""
    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2{sfx}.parquet")
    area = pd.read_parquet(DATA / f"nucleus_functional_area_assignment_v{VERSION}.parquet")
    area = area[["target_id", "tag", "value", "pt_position_x", "pt_position_y", "pt_position_z"]]
    neurons = usable.merge(area, left_on="nucleus_id", right_on="target_id", validate="1:1")
    assert len(neurons) == len(usable)
    neurons["area"] = np.where(neurons["tag"] == "V1", "V1", "HVA")
    for c in "xyz":
        neurons[f"soma_{c}_um"] = neurons[f"pt_position_{c}"].astype(float) / 1e3
    neurons = neurons[["root_id", "nucleus_id", "cell_type", "area", "tag", "axon_level", "dendrite_level",
                       "readout_loc_x", "readout_loc_y", "pref_ori", "gOSI", "soma_x_um", "soma_y_um", "soma_z_um"]]
    neurons.to_parquet(DATA / f"l2l_neurons_v{VERSION}{sfx}.parquet")

    pairs = pd.read_parquet(DATA / f"pairs_usable_v{VERSION}.parquet")
    pairs = pairs[pairs["pre_axon_level"] == "axon_fully_extended"]
    pairs = pairs[pairs["pre_root_id"].isin(neurons["root_id"]) & pairs["post_root_id"].isin(neurons["root_id"])]
    syn = pd.read_parquet(DATA / f"synapses_usable_v{VERSION}.parquet")
    syn = syn[(syn["pre_root_id"] != syn["post_root_id"]) & (syn["pre_axon_dist_um"] <= MAX_DIST_UM)]
    nsyn = syn.groupby(["pre_root_id", "post_root_id"]).size().rename("Nsyn")
    pairs = pairs.drop(columns=["Nsyn", "Nsyn_in_prox"]).join(nsyn, on=["pre_root_id", "post_root_id"])
    pairs["Nsyn"] = pairs["Nsyn"].fillna(0).astype(int)

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
    out["rf_dist_deg"] = rf_angle_deg(pre["readout_loc_x"], pre["readout_loc_y"],
                                      post["readout_loc_x"], post["readout_loc_y"])
    both_tuned = (pre["gOSI"].to_numpy(float) > GOSI_MIN) & (post["gOSI"].to_numpy(float) > GOSI_MIN)
    out["ori_diff_deg"] = np.where(both_tuned, ori_diff_deg(pre["pref_ori"], post["pref_ori"]), np.nan)
    assert out["Ld_um"].notna().all()
    out.to_parquet(DATA / f"l2l_pairs_v{VERSION}{sfx}.parquet")

    print(f"{len(out):,} pairs from {out['pre_root_id'].nunique()} fully extended presynaptic neurons; "
          f"{out['Nsyn'].sum():,} synapses; Ld > 0: {(out['Ld_um'] > 0).sum():,}; "
          f"orientation pairs (both gOSI > {GOSI_MIN}): {out['ori_diff_deg'].notna().sum():,}")
    print(out[["rf_dist_deg", "ori_diff_deg", "soma_dist_um"]].describe().round(2).to_string())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "")
