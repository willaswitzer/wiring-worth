"""Session 11: table for the absolute-radius truncation (pathabs, cyl) vs. complete axons.

Run from the repo root after scripts/run_session11.sh:
    .venv/bin/python scripts/summarize_abs_truncation.py
V1->V1, Nsyn and Nsyn/Ld. For each condition and level:
  - real slope and its shift from complete axons [95% paired bootstrap over presynaptic neurons];
  - real minus shuffled mean against the cell-type and the distance-matched shuffle (p two-sided,
    centred on the shuffled mean), at the levels that have shuffles; f = 1 values come from the main runs.
Writes notebooks/11_absolute_truncation.md.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import VERSION  # noqa: E402

TD = ROOT / "data" / "truncation"
PROJ = "V1->V1"
SIMS = {"rf_dist_deg": "RF distance", "ori_diff_deg": "orientation"}
MODELS = {"Nsyn": "Nsyn", "density": "Nsyn/Ld"}
COND = {"pathabs": "same path radius R from soma (um)", "cyl": "column cylinder radius Rc (um)"}


def v1(df):
    df = df[(df["proj"] == PROJ) & df["model"].isin(MODELS)].copy()
    df["level"] = df["level"].round(2)
    return df


def gap(real, perm_file, level):
    if not perm_file.exists():
        return {}
    pm = v1(pd.read_csv(perm_file))
    pm = pm[pm["level"] == round(level, 2)]
    out = {}
    for (sim, model), g in pm.groupby(["sim", "model"]):
        v = g["slope"].dropna().to_numpy()
        r = real.get((sim, model))
        if r is None or not len(v):
            continue
        mu = v.mean()
        p = (1 + np.sum(np.abs(v - mu) >= abs(r - mu))) / (1 + len(v))
        out[(sim, model)] = f"{r - mu:+.3f} (p = {p:.3f}, n = {len(v)})"
    return out


def main():
    lines = [f"# Session 11: truncation at one absolute radius (auto-generated, v{VERSION})\n",
             f"{PROJ}, slope per 10°. Shift = slope(level) − slope(complete) [95% paired bootstrap]. "
             "Gap = real − shuffled mean (p, number of shuffles).\n"]
    summ = TD / f"abs_summary_v{VERSION}.md"
    if summ.exists():
        lines += ["## Cable and synapses kept\n", summ.read_text().split("\n", 2)[2], ""]
    full = v1(pd.read_csv(TD / "trunc_slopes_poisson.csv"))
    full = full[full["set"] == "radial_1.00_0"].set_index(["sim", "model"])["slope"]
    ref = {"cell": gap(full.to_dict(), TD / "trunc_perm_poisson.csv", 1.0),
           "distance": gap(full.to_dict(), TD / "trunc_perm_poisson_nulldistance.csv", 1.0)}
    for cond, label in COND.items():
        f = TD / f"trunc_slopes_poisson_{cond}.csv"
        if not f.exists():
            lines.append(f"## {cond}: not run yet\n")
            continue
        fit = v1(pd.read_csv(f))
        fit = fit[fit["cond"] == cond]
        bf = TD / f"trunc_boot_poisson_{cond}.csv"
        ci = {}
        if bf.exists():
            b = v1(pd.read_csv(bf))
            b0 = b[b["set"] == "radial_1.00_0"].set_index(["boot", "sim", "model"])["slope"]
            b = b[b["cond"] == cond].copy()
            b["shift"] = b["slope"].to_numpy() - b0.reindex(list(zip(b["boot"], b["sim"], b["model"]))).to_numpy()
            for k, g in b.groupby(["level", "sim", "model"]):
                q = g["shift"].dropna().quantile([.025, .975]).to_numpy()
                ci[k] = f"[{q[0]:+.3f}, {q[1]:+.3f}]"
        rows = []
        for (lv, sim, model), g in fit.groupby(["level", "sim", "model"]):
            r = g["slope"].iloc[0]
            real = {(sim, model): r}
            gc = gap(real, TD / f"trunc_perm_poisson_{cond}.csv", lv).get((sim, model), "")
            gd = gap(real, TD / f"trunc_perm_poisson_{cond}_nulldistance.csv", lv).get((sim, model), "")
            rows.append({"similarity": SIMS[sim], "model": MODELS[model], label: int(lv),
                         "slope": f"{r:+.3f}", "shift": f"{r - full[(sim, model)]:+.3f} {ci.get((lv, sim, model), '')}",
                         "gap, cell-type shuffle": gc, "gap, distance-matched shuffle": gd})
        t = pd.DataFrame(rows).sort_values(["similarity", "model", label], ascending=[False, True, False])
        lines += [f"## {cond}: {label}\n", t.to_markdown(index=False), ""]
    lines += ["## Complete axons (reference)\n",
              pd.DataFrame([{"similarity": SIMS[s], "model": MODELS[m], "slope": f"{full[(s, m)]:+.3f}",
                             "gap, cell-type shuffle": ref["cell"].get((s, m), ""),
                             "gap, distance-matched shuffle": ref["distance"].get((s, m), "")}
                            for s in SIMS for m in MODELS]).to_markdown(index=False), ""]
    (ROOT / "notebooks" / "11_absolute_truncation.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
