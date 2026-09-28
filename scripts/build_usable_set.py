"""Sessions 2-3: finalize the usable set for the like-to-like baseline.

Run from the repo root: `.venv/bin/python scripts/build_usable_set.py`
Needs the Session 1 parquets in data/ (notebooks/01_setup_and_explore.py).

Pipeline (main):
  start       Session 1 set: proofread (dendrite >= clean) + in coregistration_manual_v4 (1,059),
              current (v1822) root ids, v661 labels + aibs_metamodel_celltypes_v661_corrections
  excitatory  classification_system == excitatory_neuron (after corrections)
  match       drop score < 0 matches
  joined      keep matches with a row in digital_twin_properties_bcm_coreg_v4 (session/scan_idx/unit_id);
              per neuron keep the one with the highest cc_max (ties -> higher score) (Session 3 rule);
              record the chosen scan
  quality     cc_abs > 0.2 and cc_max > 0.4 (Ding et al. 2025)

Sensitivity (supplement): lowest-residual multi-match rule (Session 2); Ding match rule
(residual < 20, score > -10); no quality filter.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import DATASTACK, VERSION  # noqa: E402

DATA = ROOT / "data"
OUT = DATA / f"usable_set_v{VERSION}_s2.parquet"
COUNTS_MD = ROOT / "notebooks" / "02_counts.md"

COREG = "coregistration_manual_v4"
FUNC = "digital_twin_properties_bcm_coreg_v4"
CT_CORR = "aibs_metamodel_celltypes_v661_corrections"
UNIT = ["session", "scan_idx", "unit_id"]
AXON_ORDER = ["axon_fully_extended", "axon_partially_extended", "axon_interareal", "none"]
AXON_SHORT = dict(zip(AXON_ORDER, ["fully_ext", "partially_ext", "interareal", "none"]))
CC_ABS_MIN, CC_MAX_MIN = 0.2, 0.4

pd.set_option("display.width", 200)


def load(name):
    return pd.read_parquet(DATA / f"{name}_v{VERSION}.parquet")


# %% Start: the Session 1 set, with root-id and cell-type corrections
from caveclient import CAVEclient  # noqa: E402

client = CAVEclient(DATASTACK, version=VERSION)
usable = load("proofread_coreg_manual").rename(columns={"valid_id": "proofread_root_id"})
assert len(usable) == 1059, len(usable)

# Root ids: proofreading_status_and_strategy records the root that was assessed (valid_id)
# and its v1822 root (pt_root_id). Everything downstream joins on the v1822 root.
proof = load("proofreading_status_and_strategy")
changed = proof[proof["root_changed"]][["valid_id", "pt_root_id", "strategy_axon", "strategy_dendrite"]]
ts = client.materialize.get_timestamp(VERSION)
root_log = []
for r in changed.itertuples():
    latest = client.chunkedgraph.suggest_latest_roots(r.valid_id, timestamp=ts)
    root_log.append(dict(proofread_root_id=r.valid_id, v1822_root_id=r.pt_root_id,
                         lineage_confirms=latest == r.pt_root_id,
                         in_session1_set=r.pt_root_id in set(usable["pt_root_id"])))
root_log = pd.DataFrame(root_log)
print("Root ids changed since proofreading:\n", root_log.to_string(index=False), "\n")
assert root_log["lineage_confirms"].all()

# Cell-type corrections: keyed by nucleus id (target_id). Apply to all rows; log what changes.
corr = client.materialize.query_table(CT_CORR)
corr = corr[["target_id", "classification_system", "cell_type"]].rename(columns={"target_id": "nucleus_id"})
usable = usable.rename(columns={"cell_type": "cell_type_v661"})
usable = usable.merge(corr, on="nucleus_id", how="left", suffixes=("_v661", "_corr"))
usable["cell_type_corrected"] = usable["cell_type"].notna()
usable["cell_type"] = usable["cell_type"].fillna(usable["cell_type_v661"])
usable["classification_system"] = usable["classification_system_corr"].fillna(usable["classification_system_v661"])
proof_ct = load("proofread_celltypes")
corr_log = corr.merge(proof_ct[["nucleus_id", "pt_root_id", "cell_type"]], on="nucleus_id", suffixes=("_corrected", "_v661"))
corr_log["in_session1_set"] = corr_log["nucleus_id"].isin(usable["nucleus_id"])
print(f"Cell-type corrections ({CT_CORR}) hitting proofread cells:\n",
      corr_log.drop(columns="classification_system").to_string(index=False), "\n")


def celltype_by_axon(df):
    t = pd.crosstab(df["cell_type"], df["strategy_axon"]).reindex(columns=AXON_ORDER, fill_value=0)
    t = t.rename(columns=AXON_SHORT)
    t["All"] = t.sum(axis=1)
    t.loc["**All**"] = t.sum()
    return t


# %% Match rules + best-fit selection

coreg = load(COREG)
coreg = coreg[coreg["pt_root_id"].isin(usable["pt_root_id"])]
func = load(FUNC)
FUNC_COLS = ["OSI", "DSI", "gOSI", "gDSI", "pref_ori", "pref_dir", "readout_loc_x", "readout_loc_y",
             "cc_abs", "cc_max", "cc_norm"]


def pick_best(matches, by):
    """One row per neuron, best by `by` = [(column, ascending), ...]; ties fall through the list."""
    cols, asc = zip(*by)
    m = matches.sort_values(["pt_root_id", *cols], ascending=[True, *asc])
    n = matches.groupby("pt_root_id").agg(n_matches=("unit_id", "size"),
                                          n_scans=("scan_idx", lambda s: matches.loc[s.index, ["session", "scan_idx"]].drop_duplicates().shape[0]))
    best = m.drop_duplicates("pt_root_id").set_index("pt_root_id").join(n).reset_index()
    best["chosen_scan"] = "s" + best["session"].astype(str) + "_scan" + best["scan_idx"].astype(str)
    return best


def run(match_rule, select="cc_max", quality=True):
    """Return {stage: neuron table} for one pipeline variant.

    select="cc_max"   (main) among passing matches with a DT row, highest cc_max (ties: higher score)
    select="residual" (S2 rule) among passing matches, lowest residual (ties: higher score), then join DT
    Stage counts are neurons: "match rules" = >= 1 passing match, "joined DT" = neuron has DT data
    for its chosen match (cc_max) / its best-fit match (residual).
    """
    stages = {"start": usable}
    exc = usable[usable["classification_system"] == "excitatory_neuron"]
    stages["excitatory"] = exc
    ok = coreg[coreg["pt_root_id"].isin(exc["pt_root_id"]) & match_rule(coreg)]
    stages["match rules"] = exc[exc["pt_root_id"].isin(ok["pt_root_id"])]
    keep = ["pt_root_id", *UNIT, "field", "residual", "score"]
    if select == "cc_max":
        ok_dt = ok[keep].merge(func[UNIT + FUNC_COLS], on=UNIT)
        best = pick_best(ok_dt, [("cc_max", False), ("score", False)])
    else:
        best = pick_best(ok[keep], [("residual", True), ("score", False)])
        best = best.merge(func[UNIT + FUNC_COLS], on=UNIT)
    joined = exc.merge(best, on="pt_root_id")
    stages["joined DT"] = joined
    if quality:
        stages["quality"] = joined[(joined["cc_abs"] > CC_ABS_MIN) & (joined["cc_max"] > CC_MAX_MIN)]
    return stages


MAIN_RULE = lambda d: d["score"] >= 0                                  # noqa: E731
DING_RULE = lambda d: (d["residual"] < 20) & (d["score"] > -10)        # noqa: E731
MAIN = "main (score >= 0, max cc_max, QC)"
variants = {
    MAIN: run(MAIN_RULE),
    "lowest residual (S2 rule), QC": run(MAIN_RULE, select="residual"),
    "Ding match (residual < 20, score > -10), max cc_max, QC": run(DING_RULE),
    "main, no QC": run(MAIN_RULE, quality=False),
    "lowest residual, no QC": run(MAIN_RULE, select="residual", quality=False),
    "Ding match, no QC": run(DING_RULE, quality=False),
}
main = variants[MAIN]

# How often the two multi-match rules pick different units (neurons in both final sets).
a = main["quality"].set_index("pt_root_id")[UNIT]
b = variants["lowest residual (S2 rule), QC"]["quality"].set_index("pt_root_id")[UNIT]
both = a.index.intersection(b.index)
n_diff = (a.loc[both] != b.loc[both]).any(axis=1).sum()
multi = (main["quality"]["n_matches"] > 1).sum()
rule_note = (f"Final-set neurons with >1 passing DT match: {multi}. Neurons in both the cc_max and "
             f"lowest-residual final sets: {len(both)}; chosen unit differs for {n_diff}.")
print(rule_note, "\n")

# %% Counts tables
lines = [f"# Session 2 counts (auto-generated by scripts/build_usable_set.py, v{VERSION})\n"]
lines.append("## Main pipeline: cell type x axon level at each stage\n")
for stage, df in main.items():
    lines += [f"### {stage} (n = {len(df):,})\n", celltype_by_axon(df).to_markdown(), ""]

lines.append("## Stage totals, all variants\n")
tot = pd.DataFrame({v: {s: len(d) for s, d in st.items()} for v, st in variants.items()}).T
lines += [tot.fillna("").to_markdown(), ""]

lines.append("## Fully extended axons surviving each stage, by cell type\n")
full = pd.DataFrame({
    f"{v} / {s}": d[d["strategy_axon"] == "axon_fully_extended"]["cell_type"].value_counts()
    for v, st in variants.items() for s, d in st.items() if s != "start"
}).fillna(0).astype(int)
full.loc["**All**"] = full.sum()
lines += [full.T.to_markdown(), ""]

lines.append("## Final set by axon level, sensitivity variants (supplement)\n")
for v, st in variants.items():
    if v == MAIN:
        continue
    final = list(st.values())[-1]
    lines += [f"### {v} (n = {len(final):,})\n", celltype_by_axon(final).to_markdown(), ""]

lines.append("## Provenance\n")
lines += ["Root ids changed since proofreading:\n", root_log.to_markdown(index=False), "",
          f"Cell-type corrections from `{CT_CORR}` (v{VERSION}) on proofread cells:\n",
          corr_log.to_markdown(index=False), "",
          rule_note + "\n"]
COUNTS_MD.write_text("\n".join(lines))
print("\n".join(lines))

# %% Save final set
final = main["quality"].rename(columns={"pt_root_id": "root_id", "strategy_axon": "axon_level",
                                        "strategy_dendrite": "dendrite_level"})
final["interareal"] = final["axon_level"] == "axon_interareal"
final["root_changed_since_proofreading"] = final["root_id"] != final["proofread_root_id"]
cols = ["root_id", "nucleus_id", "cell_type", "cell_type_v661", "cell_type_corrected", "axon_level",
        "dendrite_level", "interareal", *UNIT, "OSI", "DSI", "gOSI", "gDSI", "pref_ori", "pref_dir",
        "readout_loc_x", "readout_loc_y", "cc_abs", "cc_max", "cc_norm", "chosen_scan", "n_scans",
        "n_matches", "field", "residual", "score", "proofread_root_id", "root_changed_since_proofreading"]
final = final[cols].reset_index(drop=True)
final.attrs.update(version=str(VERSION), coreg=COREG, func=FUNC)
final.to_parquet(OUT)
print(f"\nsaved {len(final)} neurons -> {OUT.relative_to(ROOT)}")
