"""Session 5: all input and output synapses of every usable neuron (for axon/dendrite relabeling).

Run from the repo root: `.venv/bin/python scripts/fetch_io_synapses.py`
Queries the datastack synapse table (synapses_pni_2 at v1822) for pre_pt_root_id in chunk
(outputs) and post_pt_root_id in chunk (inputs), for all 694 usable roots, any partner.
Autapses (pre == post) are dropped. Resumable: one cached parquet per chunk in data/skeletons/io/.
"""
import sys
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import DATASTACK, VERSION  # noqa: E402
from caveclient import CAVEclient  # noqa: E402

DATA = ROOT / "data"
IO_DIR = DATA / "skeletons" / "io"
OUT = DATA / f"synapses_io_usable_v{VERSION}.parquet"
CHUNK = 10
WORKERS = 4
RETRIES = 3


def fetch_chunk(client, syn_table, k, roots):
    path = IO_DIR / f"chunk_{k:03d}.parquet"
    if path.exists():
        return k, None
    parts = []
    for col, direction in [("pre_pt_root_id", "out"), ("post_pt_root_id", "in")]:
        for attempt in range(RETRIES):
            try:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter("always")
                    df = client.materialize.query_table(
                        syn_table, filter_in_dict={col: roots},
                        select_columns=["id", "pre_pt_root_id", "post_pt_root_id", "size", "ctr_pt_position"],
                        desired_resolution=[1, 1, 1], split_positions=True)
                if any("Limited query" in str(x.message) for x in w):
                    return k, f"{direction}: truncated"
                break
            except Exception as e:  # noqa: BLE001
                if attempt == RETRIES - 1:
                    return k, f"{direction}: {repr(e)[:200]}"
                time.sleep(10 * (attempt + 1))
        df = df[df["pre_pt_root_id"] != df["post_pt_root_id"]]
        df = df.assign(root_id=df[col], direction=direction)
        parts.append(df)
    pd.concat(parts).to_parquet(path)
    return k, None


def main():
    IO_DIR.mkdir(parents=True, exist_ok=True)
    client = CAVEclient(DATASTACK, version=VERSION)
    syn_table = client.info.get_datastack_info()["synapse_table"]
    roots = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")["root_id"].astype("int64").tolist()
    chunks = [roots[i:i + CHUNK] for i in range(0, len(roots), CHUNK)]
    t0 = time.time()
    fails = {}
    with ThreadPoolExecutor(WORKERS) as ex:
        for n, (k, err) in enumerate(ex.map(lambda a: fetch_chunk(client, syn_table, *a), enumerate(chunks)), 1):
            if err:
                fails[k] = err
            if n % 10 == 0:
                print(f"  {n}/{len(chunks)} chunks, {len(fails)} failed, {time.time() - t0:.0f}s", flush=True)
    if fails:
        print("FAILED chunks (rerun to retry):", fails)
        sys.exit(1)
    io = pd.concat([pd.read_parquet(p) for p in sorted(IO_DIR.glob("chunk_*.parquet"))], ignore_index=True)
    io = io.rename(columns={"id": "synapse_id", "ctr_pt_position_x": "ctr_x_nm", "ctr_pt_position_y": "ctr_y_nm",
                            "ctr_pt_position_z": "ctr_z_nm"})
    io = io[["root_id", "direction", "synapse_id", "pre_pt_root_id", "post_pt_root_id", "size",
             "ctr_x_nm", "ctr_y_nm", "ctr_z_nm"]]
    io.to_parquet(OUT)
    print(f"{syn_table}: {len(io):,} rows ({(io.direction == 'out').sum():,} out, {(io.direction == 'in').sum():,} in) "
          f"for {io.root_id.nunique()} roots -> {OUT.relative_to(ROOT)}; {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
