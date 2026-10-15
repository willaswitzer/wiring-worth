"""Session 10: sensitivity of the like-to-like baseline to the co-registration match rules.

Run from the repo root after scripts/build_usable_set.py:
    .venv/bin/python scripts/match_sensitivity.py

The main usable set (Session 2-3) keeps manual co-registration matches with score >= 0 and, for neurons
with several digital-twin (DT) units, the unit with the highest cc_max. This script re-chooses each
neuron's functional unit under two alternative rules and writes a variant usable set for each:

  residual   score >= 0; among passing matches the lowest residual (ties: higher score), then join DT
             (the Session 2 multi-match rule)
  ding       Ding et al. 2025 match rule (residual < 20, score > -10); among passing matches with a DT
             row, the highest cc_max (ties: higher score)

Both keep the Ding quality filter (cc_abs > 0.2, cc_max > 0.4).

Scope: only neurons in the main usable set are re-matched, because proximity (Ld) and synapses were
computed for that set only (pairs_usable). A neuron whose variant match fails the rule or the quality
filter is dropped. Neurons that would enter only under a variant rule are counted, not added.
Writes data/usable_set_v{VERSION}_s2_{tag}.parquet and notebooks/10_match_rules.md.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import VERSION  # noqa: E402

DATA = ROOT / "data"
UNIT = ["session", "scan_idx", "unit_id"]
FUNC_COLS = ["OSI", "DSI", "gOSI", "gDSI", "pref_ori", "pref_dir", "readout_loc_x", "readout_loc_y",
             "cc_abs", "cc_max", "cc_norm"]
CC_ABS_MIN, CC_MAX_MIN = 0.2, 0.4
RULES = {
    "main": lambda d: d["score"] >= 0,
    "residual": lambda d: d["score"] >= 0,
    "ding": lambda d: (d["residual"] < 20) & (d["score"] > -10),
}


def load(name):
    return pd.read_parquet(DATA / f"{name}_v{VERSION}.parquet")


def choose(coreg, dt, rule, select):
    ok = coreg[RULES[rule](coreg)]
    if select == "cc_max":
        j = ok.merge(dt, on=UNIT, how="inner")
        j = j.sort_values(["target_id", "cc_max", "score"], ascending=[True, False, False])
    else:  # lowest residual among passing matches, then join DT (no DT row -> dropped)
        j = ok.sort_values(["target_id", "residual", "score"], ascending=[True, True, False])
        j = j.drop_duplicates("target_id").merge(dt, on=UNIT, how="inner")
    j = j.drop_duplicates("target_id")
    return j[(j["cc_abs"] > CC_ABS_MIN) & (j["cc_max"] > CC_MAX_MIN)].set_index("target_id")


def main():
    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    coreg = load("coregistration_manual_v4")[["target_id", *UNIT, "field", "residual", "score"]]
    dt = load("digital_twin_properties_bcm_coreg_v4")[[*UNIT, *FUNC_COLS]].drop_duplicates(UNIT)
    s1 = load("proofread_coreg_manual")
    s1_exc = set(s1.loc[s1["classification_system"] == "excitatory_neuron", "nucleus_id"])

    in_main = coreg[coreg["target_id"].isin(usable["nucleus_id"])]
    # Check: the main rule reproduces the Session 2-3 choices exactly.
    m = choose(in_main, dt, "main", "cc_max")
    chk = usable.set_index("nucleus_id")[UNIT]
    assert set(m.index) == set(chk.index), "main rule does not reproduce the usable set"
    assert (m.loc[chk.index, UNIT].to_numpy() == chk.to_numpy()).all(), "main rule picks different units"

    rows = []
    for tag, select in [("residual", "residual"), ("ding", "cc_max")]:
        v = choose(in_main, dt, tag, select)
        keep = usable[usable["nucleus_id"].isin(v.index)].copy()
        vv = v.loc[keep["nucleus_id"]]
        changed = (keep[UNIT].to_numpy() != vv[UNIT].to_numpy()).any(axis=1)
        for c in [*UNIT, *FUNC_COLS, "field", "residual", "score"]:
            keep[c] = vv[c].to_numpy()
        keep["chosen_scan"] = "s" + keep["session"].astype(str) + "_scan" + keep["scan_idx"].astype(str)
        keep.to_parquet(DATA / f"usable_set_v{VERSION}_s2_{tag}.parquet")
        # neurons outside the main set that this rule would admit (excitatory, proofread, in coreg v4);
        # Session 1 cell-type corrections are not applied here, so this count is approximate
        out_set = coreg[coreg["target_id"].isin(s1_exc) & ~coreg["target_id"].isin(usable["nucleus_id"])]
        n_new = len(choose(out_set, dt, tag, select))
        fe = keep["axon_level"] == "axon_fully_extended"
        rows.append(dict(rule=tag, neurons=len(keep), dropped=len(usable) - len(keep),
                         unit_changed=int(changed.sum()), fully_extended_pre=int(fe.sum()),
                         fully_extended_unit_changed=int((changed & fe.to_numpy()).sum()),
                         would_enter_not_added=n_new))
        print(rows[-1])

    t = pd.DataFrame(rows)
    md = [f"# Session 10: co-registration match-rule variants (v{VERSION}, auto-generated)\n",
          f"Main usable set: {len(usable)} neurons, {(usable['axon_level'] == 'axon_fully_extended').sum()} "
          "with fully extended axons. Variants re-choose the functional unit for neurons in the main set only.\n",
          t.to_markdown(index=False), ""]
    (ROOT / "notebooks" / "10_match_rules.md").write_text("\n".join(md))


if __name__ == "__main__":
    main()
