"""Session 11: how much does the newest public release change the axons this analysis uses?

Run from the repo root (needs the CAVE token; ~5-10 min):
    .venv/bin/python scripts/check_new_release.py

The analysis stays pinned to v1822 (config.py). This only answers two questions for the preprint:
1. Of the 193 fully extended presynaptic axons (v1822), how many keep that status in the newest release?
2. How many co-registered excitatory neurons are fully extended in the newest release but were not in v1822,
   by cell type (4P in particular, since the newest release reportedly improves layer 4)?
Neurons are matched across versions by nucleus id. Writes notebooks/11_new_release.md.
"""
import sys
from pathlib import Path

import pandas as pd
from caveclient import CAVEclient

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import DATASTACK, VERSION  # noqa: E402

DATA = ROOT / "data"
CT_TABLE = "aibs_metamodel_celltypes_v661"
COREG = "coregistration_manual_v4"


def by_nucleus(client, name, cols):
    """Query a table and key it by nucleus id (target_id for reference tables, else via the nucleus table)."""
    df = client.materialize.query_table(name)
    if "target_id" in df:
        return df.rename(columns={"target_id": "nucleus_id"})[["nucleus_id", *cols]]
    nuc = client.materialize.query_table("nucleus_detection_v0")[["id", "pt_root_id"]]
    nuc = nuc[nuc["pt_root_id"] != 0].drop_duplicates("pt_root_id", keep=False)
    return df.merge(nuc.rename(columns={"id": "nucleus_id"}), on="pt_root_id")[["nucleus_id", *cols]]


def main():
    client = CAVEclient(DATASTACK)
    versions = sorted(client.materialize.get_versions())
    latest = versions[-1]
    lines = [f"# Session 11: newest public release vs. pinned v{VERSION} (auto-generated)\n",
             f"Available versions: {', '.join(map(str, versions))}. Newest: v{latest} "
             f"(materialized {client.materialize.get_timestamp(latest):%Y-%m-%d}).\n"]
    if latest == VERSION:
        lines.append("The pinned version is the newest; nothing to compare.\n")
        (ROOT / "notebooks" / "11_new_release.md").write_text("\n".join(lines))
        print("\n".join(lines))
        return
    client.version = latest
    tables = client.materialize.get_tables()
    for t in ("proofreading_status_and_strategy", CT_TABLE, COREG):
        if t not in tables:
            near = [x for x in tables if t.split("_")[0] in x]
            raise SystemExit(f"{t} not in v{latest}; similar tables: {near}")

    proof = by_nucleus(client, "proofreading_status_and_strategy", ["strategy_axon", "strategy_dendrite"])
    ct = by_nucleus(client, CT_TABLE, ["classification_system", "cell_type"])
    co = by_nucleus(client, COREG, ["score"])
    co = co[co["score"] >= 0].drop_duplicates("nucleus_id")
    new = proof.drop_duplicates("nucleus_id").merge(ct.drop_duplicates("nucleus_id"), on="nucleus_id", how="left")
    new["coreg"] = new["nucleus_id"].isin(co["nucleus_id"])

    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    ours = usable[usable["axon_level"] == "axon_fully_extended"]
    old_proof = pd.read_parquet(DATA / f"proofread_celltypes_v{VERSION}.parquet")
    old_full = set(old_proof.loc[old_proof["strategy_axon"] == "axon_fully_extended", "nucleus_id"].dropna())

    # 1. the 193
    m = ours[["nucleus_id", "cell_type"]].merge(new[["nucleus_id", "strategy_axon"]], on="nucleus_id", how="left")
    m["strategy_axon"] = m["strategy_axon"].fillna("(not in table)")
    lines += [f"## 1. The {len(ours)} fully extended presynaptic axons (v{VERSION}) in v{latest}\n",
              pd.crosstab(m["cell_type"], m["strategy_axon"], margins=True).to_markdown(), ""]

    # 2. newly fully extended, co-registered excitatory neurons
    exc = new[(new["classification_system"] == "excitatory_neuron") & new["coreg"]]
    full_new = exc[(exc["strategy_axon"] == "axon_fully_extended") & ~exc["nucleus_id"].isin(old_full)]
    full_all = exc[exc["strategy_axon"] == "axon_fully_extended"]
    tab = pd.DataFrame({
        f"fully extended, co-registered, v{latest}": full_all["cell_type"].value_counts(),
        f"of which not fully extended in v{VERSION}": full_new["cell_type"].value_counts(),
    }).fillna(0).astype(int)
    tab.loc["All"] = tab.sum()
    lines += [f"## 2. Co-registered excitatory neurons (manual co-registration, score >= 0) with fully extended axons\n",
              "Before the functional-quality filter and the dendrite requirement, so an upper bound on what a rerun "
              "would add.\n", tab.to_markdown(), ""]
    (ROOT / "notebooks" / "11_new_release.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
