"""Session 4: fetch skeletons for all usable neurons from the CAVE skeleton service.

Run from the repo root after scripts/build_usable_set.py:
    .venv/bin/python scripts/fetch_skeletons.py

One skeleton per root gives both the axon (for the 585 presynaptic neurons) and the dendrite
(for all 694). Skeletons are pcg_skel skeletons served by the CAVE skeleton service
(skeleton_version 4), with SWC compartment labels (1 soma, 2 axon, 3 dendrite); vertices in nm.
Root ids are immutable, so the skeleton of a v1822 root is the v1822 object.

Resumable: roots with a cached file in data/skeletons/ are skipped. Failures are retried
(RETRIES times) and then logged to data/skeletons/failures.csv; rerun to retry them.
"""
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import DATASTACK, VERSION  # noqa: E402
from caveclient import CAVEclient  # noqa: E402

DATA = ROOT / "data"
SKEL_DIR = DATA / "skeletons"
SKEL_VERSION = 4
WORKERS = 6
BATCH = 60
RETRIES = 3


def skel_path(root_id):
    return SKEL_DIR / f"{root_id}_sv{SKEL_VERSION}.npz"


def fetch_one(client, root_id):
    last = None
    for attempt in range(RETRIES):
        try:
            s = client.skeleton.get_skeleton(int(root_id), skeleton_version=SKEL_VERSION, output_format="dict")
            m = s["meta"]
            tmp = skel_path(root_id).with_suffix(".tmp.npz")
            np.savez_compressed(
                tmp,
                vertices=np.asarray(s["vertices"], dtype=np.float64),
                edges=np.asarray(s["edges"], dtype=np.int64),
                compartment=np.asarray(s["compartment"], dtype=np.int8),
                root=np.int64(s["root"]),
                soma_pt=np.array([m["soma_pt_x"], m["soma_pt_y"], m["soma_pt_z"]]),
                skeleton_version=np.int64(m["skeleton_version"]),
            )
            tmp.rename(skel_path(root_id))  # atomic: a partial file never looks cached
            return root_id, None
        except Exception as e:  # noqa: BLE001
            last = repr(e)[:300]
            time.sleep(5 * (attempt + 1))
    return root_id, last


def main():
    SKEL_DIR.mkdir(parents=True, exist_ok=True)
    client = CAVEclient(DATASTACK, version=VERSION)
    usable = pd.read_parquet(DATA / f"usable_set_v{VERSION}_s2.parquet")
    todo = [r for r in usable["root_id"].astype("int64") if not skel_path(r).exists()]
    print(f"{len(usable)} usable roots, {len(usable) - len(todo)} cached, {len(todo)} to fetch")

    t0 = time.time()
    failures = {}
    for b in range(0, len(todo), BATCH):
        batch = todo[b:b + BATCH]
        with ThreadPoolExecutor(WORKERS) as ex:
            for fut in as_completed([ex.submit(fetch_one, client, r) for r in batch]):
                r, err = fut.result()
                if err:
                    failures[r] = err
        done = b + len(batch)
        print(f"  {done}/{len(todo)} fetched, {len(failures)} failed, {time.time() - t0:.0f}s elapsed", flush=True)

    fail_csv = SKEL_DIR / "failures.csv"
    pd.DataFrame({"root_id": list(failures), "error": list(failures.values())}).to_csv(fail_csv, index=False)
    n_cached = sum(skel_path(r).exists() for r in usable["root_id"].astype("int64"))
    print(f"done in {time.time() - t0:.0f}s: {n_cached}/{len(usable)} cached, {len(failures)} failed -> {fail_csv.name}")


if __name__ == "__main__":
    main()
