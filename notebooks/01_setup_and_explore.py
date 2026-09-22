# %% [markdown]
# # Mouse project — Session 1: access + proofreading/completeness inventory
# Run cell by cell in Jupyter or VS Code (the "# %%" markers make cells).
# Goal: CAVE access working, analysis version pinned, proofreading levels
# tabulated by cell type, and a first look at co-registration tables.

# %%
# !pip install caveclient pandas pyarrow

# %% [markdown]
# ## 1. One-time auth (skip after the first run)
# 1. Accept the MICrONS terms of service (Google account):
#    https://global.daf-apis.com/sticky_auth/api/v1/tos/2/accept
# 2. Run the cell below, log in, copy the token, paste it into save_token.
#    Note: make_new=True invalidates any older token you may have.

# %%
from caveclient import CAVEclient
# client = CAVEclient()
# client.auth.setup_token(make_new=True)
# client.auth.save_token(token="PASTE_TOKEN_HERE")

# %% [markdown]
# ## 2. Connect and pin a version

# %%
import pandas as pd

client = CAVEclient("minnie65_public")
versions = sorted(client.materialize.get_versions())
for v in versions:
    print(v, client.materialize.get_timestamp(v))

# Pin the newest public version and write it in the project doc.
# Every analysis in the project should use this same number.
VERSION = versions[-1]
client.version = VERSION
print("Pinned version:", VERSION)

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
print(proof["strategy_axon"].value_counts(), "\n")
print(proof["strategy_dendrite"].value_counts(), "\n")
print(pd.crosstab(proof["strategy_axon"], proof["strategy_dendrite"]))

# Cells whose root id changed since assessment -> flag, handle carefully
proof["root_changed"] = proof["valid_id"] != proof["pt_root_id"]
print("\nroot id changed since assessment:", proof["root_changed"].sum())

# %% [markdown]
# ## 4. Join cell types (the confound check for December)

# %%
ct = client.materialize.tables.aibs_metamodel_celltypes_v661().query(
    select_columns={
        "nucleus_detection_v0": ["pt_root_id", "id"],
        "aibs_metamodel_celltypes_v661": ["classification_system", "cell_type"],
    }
)
ct = ct.drop_duplicates("pt_root_id", keep=False).rename(columns={"id": "nucleus_id"})

proof_ct = proof.merge(
    ct[["pt_root_id", "nucleus_id", "classification_system", "cell_type"]],
    on="pt_root_id", how="left",
)
print(pd.crosstab(proof_ct["cell_type"], proof_ct["strategy_axon"], margins=True))

# %% [markdown]
# ## 5. Find the co-registration / functional tables
# Don't assume names: list what this version actually has.

# %%
tables = client.materialize.get_tables()
for t in sorted(tables):
    if any(k in t.lower() for k in ["coreg", "functional", "func", "manual_central"]):
        print(t)

# %% [markdown]
# Pick the co-registration table from the list above, then:
#   coreg = client.materialize.query_table("<TABLE_NAME>")
#   print(len(coreg), coreg.columns.tolist())
#   usable = proof_ct.merge(coreg, on="pt_root_id", how="inner")
#   print(len(usable), "proofread + co-registered")
#   print(pd.crosstab(usable["cell_type"], usable["strategy_axon"]))
# That last count is the number that sizes the core model.

# %%
proof_ct.to_parquet(f"data/proofread_celltypes_v{VERSION}.parquet")
print("saved")
