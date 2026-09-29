from caveclient import CAVEclient

client = CAVEclient("minnie65_public")
print("available versions:", client.materialize.get_versions()[-5:])
client.version = 1822

dt_tables = [t for t in client.materialize.get_tables() if "digital_twin" in t]
print("digital-twin tables:", dt_tables)

with open("similarity_metric_check.md", "w") as f:
    f.write("# Digital-twin table check (v1822)\n\n")
    for t in dt_tables:
        meta = client.materialize.get_table_metadata(t)
        df = client.materialize.query_table(t, limit=5)
        f.write(f"## {t}\n")
        f.write(f"schema: {meta.get('schema_type')}\n\n")
        f.write(f"description: {meta.get('description')}\n\n")
        f.write(f"columns: {df.columns.tolist()}\n\n")
        f.write(df.head().to_markdown() + "\n\n")
print("wrote similarity_metric_check.md")