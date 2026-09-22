# %% [markdown]
# # Mouse project — Session 1: access + proofreading/completeness inventory
# Run cell by cell in Jupyter or VS Code (the "# %%" markers make cells),
# or end to end from the repo root: `.venv/bin/python notebooks/01_setup_and_explore.py`.
# Goal: CAVE access working, analysis version pinned, proofreading levels
# tabulated by cell type, and a first look at co-registration tables.

# %%
# !pip install caveclient pandas pyarrow

# %% [markdown]
# ## 1. One-time auth (skip after the first run)
# 1. Accept the MICrONS terms of service (Google account):
#    https://global.daf-apis.com/sticky_auth/api/v1/tos/2/accept
# 2. Run `python scripts/save_token.py` from the repo root (prompts for the token).
#    Note: make_new=True invalidates any older token you may have.

# %%
from caveclient import CAVEclient
# client = CAVEclient()
# client.auth.setup_token(make_new=True)
# client.auth.save_token(token="PASTE_TOKEN_HERE")

# %% [markdown]
# ## 2. Connect at the pinned version
# The version lives in `config.py`; every analysis in the project uses the same number.

# %%
import sys
from pathlib import Path
import pandas as pd

if "__file__" in globals():
    ROOT = Path(__file__).resolve().parents[1]
else:
    ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))
from config import DATASTACK, VERSION

DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 30)

client = CAVEclient(DATASTACK)
versions = sorted(client.materialize.get_versions())
for v in versions:
    print(v, client.materialize.get_timestamp(v))
assert VERSION in versions, f"pinned version {VERSION} is no longer available"
if VERSION != versions[-1]:
    print(f"Note: newer version {versions[-1]} exists; staying on pinned {VERSION}")
client.version = VERSION
print("Pinned version:", VERSION)


def save(df, name):
    df.to_parquet(DATA / f"{name}_v{VERSION}.parquet")

# %% [markdown]
# ## 3. Proofreading status = the completeness axis
# `proofreading_status_and_strategy` gives separate axon/dendrite levels:
# - dendrite: none < dendrite_clean < dendrite_extended
# - axon: none < axon_partially_extended < axon_fully_extended
#   (plus axon_column_truncated / axon_interareal, which are biased samplings)

# %%
proof = client.materialize.query_table(
    "proofreading_status_and_strategy",
    desired_resolution=[1, 1, 1],
    split_positions=True,
)
print(len(proof), "proofread cells")
print(proof["strategy_axon"].value_counts(dropna=False), "\n")
print(proof["strategy_dendrite"].value_counts(dropna=False), "\n")
print(pd.crosstab(proof["strategy_axon"], proof["strategy_dendrite"], margins=True))

# Cells whose root id changed since assessment -> flag, handle carefully
proof["root_changed"] = proof["valid_id"] != proof["pt_root_id"]
print("\nroot id changed since assessment:", proof["root_changed"].sum())
print("duplicate pt_root_id rows:", proof["pt_root_id"].duplicated().sum())
save(proof, "proofreading_status_and_strategy")

# %% [markdown]
# ## 4. Join cell types (the confound check for December)
# Reference table: `id` is the cell-type row, `target_id` the nucleus id.
# Root ids with more than one nucleus (multi-soma mergers) are dropped entirely.

# %%
CT_TABLE = "aibs_metamodel_celltypes_v661"
ct_all = client.materialize.query_table(CT_TABLE, desired_resolution=[1, 1, 1], split_positions=True)
ct_all = ct_all[ct_all["pt_root_id"] != 0]
ct = ct_all.drop_duplicates("pt_root_id", keep=False).rename(columns={"target_id": "nucleus_id"})
print(f"{CT_TABLE}: {len(ct_all)} rows with a root id, {len(ct)} after dropping multi-soma roots")
save(ct_all, CT_TABLE)

proof_ct = proof.merge(
    ct[["pt_root_id", "nucleus_id", "classification_system", "cell_type"]],
    on="pt_root_id", how="left",
)
print("proofread cells without a cell type:", proof_ct["cell_type"].isna().sum())
print(pd.crosstab(proof_ct["cell_type"].fillna("(none)"), proof_ct["strategy_axon"], margins=True))

# Manual corrections to the v661 labels exist; report the overlap but keep v661 labels for now.
corr = client.materialize.query_table(f"{CT_TABLE}_corrections")
corr = corr[corr["pt_root_id"].isin(proof_ct["pt_root_id"])]
m = corr.merge(proof_ct[["pt_root_id", "cell_type"]], on="pt_root_id", suffixes=("_corrected", "_v661"))
print(f"\nproofread cells with a correction entry: {len(corr)}; "
      f"label actually differs: {(m['cell_type_corrected'] != m['cell_type_v661']).sum()}")
save(proof_ct, "proofread_celltypes")

# %% [markdown]
# ## 5. Co-registration tables
# Don't assume names: list what this version actually has.

# %%
tables = sorted(client.materialize.get_tables())
coreg_like = [t for t in tables if any(k in t.lower() for k in ["coreg", "functional", "func"])]
print("\n".join(coreg_like))

# %% [markdown]
# A unique functional unit is (session, scan_idx, unit_id); one EM neuron can match
# several units if it was imaged in more than one scan, so counts below are unique neurons.
# "Proofread" = dendrite at least clean; "axon extended" = partially or fully extended.

# %%
COREG_TABLES = [t for t in tables if t.startswith("coregistration_") or t == "apl_functional_coreg_vess_fwd"]
DENDRITE_OK = ["dendrite_clean", "dendrite_extended"]
AXON_OK = ["axon_partially_extended", "axon_fully_extended"]

proofread = proof_ct[proof_ct["strategy_dendrite"].isin(DENDRITE_OK)]
print("proofread (dendrite >= clean):", len(proofread), "\n")

coreg, rows = {}, []
for t in COREG_TABLES:
    df = client.materialize.query_table(t, desired_resolution=[1, 1, 1], split_positions=True)
    df = df[df["pt_root_id"] != 0]
    coreg[t] = df
    save(df, t)
    roots = set(df["pt_root_id"])
    pc = proofread[proofread["pt_root_id"].isin(roots)]
    rows.append(dict(
        table=t, rows=len(df), unique_neurons=len(roots),
        proofread_coreg=len(pc),
        axon_extended=pc["strategy_axon"].isin(AXON_OK).sum(),
        axon_full=(pc["strategy_axon"] == "axon_fully_extended").sum(),
        axon_partial=(pc["strategy_axon"] == "axon_partially_extended").sum(),
    ))
coreg_summary = pd.DataFrame(rows).set_index("table")
print(coreg_summary)

# %% [markdown]
# Primary table: the manual co-registration (highest-confidence matches).
# Compare the automated tables against it.

# %%
PRIMARY = "coregistration_manual_v4"
print("columns:", coreg[PRIMARY].columns.tolist())

usable = proofread[proofread["pt_root_id"].isin(set(coreg[PRIMARY]["pt_root_id"]))]
print(f"\n{len(usable)} proofread + co-registered ({PRIMARY})")
print(usable["strategy_axon"].value_counts())
print(pd.crosstab(usable["cell_type"].fillna("(none)"), usable["strategy_axon"], margins=True))
save(usable, "proofread_coreg_manual")

pr = set(proofread["pt_root_id"])
man = set(coreg[PRIMARY]["pt_root_id"]) & pr
for t in COREG_TABLES:
    if t == PRIMARY:
        continue
    auto = set(coreg[t]["pt_root_id"]) & pr
    print(f"\n{t} vs manual (proofread neurons): both {len(auto & man)}, "
          f"auto only {len(auto - man)}, manual only {len(man - auto)}")
    # For neurons in both, do they map to the same functional unit(s)?
    key = ["pt_root_id", "session", "scan_idx", "unit_id"]
    a = coreg[t][coreg[t]["pt_root_id"].isin(auto & man)][key]
    b = coreg[PRIMARY][coreg[PRIMARY]["pt_root_id"].isin(auto & man)][key]
    agree = a.merge(b, on=key)["pt_root_id"].nunique()
    print(f"  shared neurons with at least one identical functional unit: {agree}/{len(auto & man)}")

# %% [markdown]
# ## 6. Functional-property tables (listed only; analysis comes later)

# %%
FUNC_TABLES = [t for t in tables if "digital_twin_properties" in t or t == "nucleus_functional_area_assignment"]
for t in FUNC_TABLES:
    df = client.materialize.query_table(t, desired_resolution=[1, 1, 1], split_positions=True)
    save(df, t)
    print(f"\n{t}: {len(df)} rows, {df['pt_root_id'].nunique()} unique roots")
    print("  columns:", df.columns.tolist())

# %%
print("\nsaved to", DATA)
for p in sorted(DATA.glob(f"*_v{VERSION}.parquet")):
    print(" ", p.name)
