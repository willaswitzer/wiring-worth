"""Supplementary Figure S1: the target-layer rule of 4P axons across the proofreading history (Session 14).

Run from the repo root after scripts/celltype_connectivity_history.py:
    .venv/bin/python scripts/figure_celltype_history.py
Inputs: data/celltype_history/{attachment_by_depth,shifts}_v1822.csv. Writes figures/14_celltype_history.{pdf,png}.
(A) Fraction of the complete axons' V1->V1 synapses attached in November 2025, by depth relative to the presynaptic
soma. (B) Shift of 4P axons' target-type log odds ratio from the complete state at November 2025 (filled, 95% paired
bootstrap over axons) and under count-matched random loss (open, 2.5-97.5% of 20 repeats), for the pooled,
distance-stratified and depth x lateral-stratified estimates.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import VERSION  # noqa: E402

D = ROOT / "data" / "celltype_history"
dep = pd.read_csv(D / f"attachment_by_depth_v{VERSION}.csv")
sh = pd.read_csv(D / f"shifts_v{VERSION}.csv")

plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True,
                               gridspec_kw={"width_ratios": [1, 1.5]})

# (A) attachment by depth, depth on the vertical axis (above the soma at the top)
bands = list(dict.fromkeys(dep["synapse depth vs. presynaptic soma (um)"]))
ypos = {b: i for i, b in enumerate(bands)}
for ct, col in (("4P", "tab:blue"), ("23P", "tab:orange")):
    d = dep[dep.presynaptic == ct]
    y = d["synapse depth vs. presynaptic soma (um)"].map(ypos)
    axA.plot(d["attached Nov 2025"], y, "-", color=col, lw=1.2)
    small = d["synapses (complete)"] < 20
    axA.plot(d["attached Nov 2025"][~small], y[~small], "o", color=col, ms=4, label=f"{ct} axons")
    for _, r in d[small].iterrows():  # too few synapses to read: open marker + n
        yy = ypos[r["synapse depth vs. presynaptic soma (um)"]]
        axA.plot(r["attached Nov 2025"], yy, "o", mfc="white", mec=col, ms=4)
        axA.text(r["attached Nov 2025"], yy + 0.3, f"n = {int(r['synapses (complete)'])}", va="top", ha="center",
                 fontsize=6, color="0.35")
axA.set_yticks(range(len(bands)))
axA.set_yticklabels(bands, fontsize=6.5)
axA.invert_yaxis()
axA.set_xlim(0, 1)
axA.set_xlabel("synapses attached, Nov 2025")
axA.set_ylabel("synapse depth vs. presynaptic soma (µm)")
axA.legend(frameon=False, fontsize=6.5, loc="lower left")
axA.set_title("A", loc="left", fontweight="bold")

# (B) 4P target-type preference: shift at Nov 2025 vs. count-matched random loss
contr = ["23P vs 4P", "5P-IT vs 4P", "5P-ET vs 4P"]
kinds = [("raw", "pooled", "0.15"), ("distance-stratified", "matched on soma distance", "tab:red"),
         ("depth x lateral-stratified", "matched on depth × lateral", "tab:purple")]
s4 = sh[sh.presynaptic == "4P"]
for i, c in enumerate(contr):
    for k, (kind, lab, col) in enumerate(kinds):
        y0 = i + (k - 1) * 0.24
        for st, dy, mfc in (("Nov 2025", -0.05, col), ("random loss (Nov 2025 counts)", 0.05, "white")):
            r = s4[(s4.contrast == c) & (s4.kind == kind) & (s4.state == st)].iloc[0]
            axB.errorbar(r["shift"], y0 + dy, xerr=[[r["shift"] - r.lo], [r.hi - r["shift"]]], fmt="o", color=col,
                         mfc=mfc, ms=3.5, capsize=1.5, lw=0.9,
                         label=(lab if st == "Nov 2025" and i == 0 else None))
axB.axvline(0, color="k", lw=0.6)
axB.set_yticks(range(len(contr)))
axB.set_yticklabels([f"{c.split(' vs ')[0]} vs 4P targets" for c in contr], fontsize=7)
axB.invert_yaxis()
axB.set_xlabel("shift of log odds ratio from complete axons")
h, l = axB.get_legend_handles_labels()
h.append(plt.Line2D([], [], marker="o", ls="", mfc="white", mec="0.3", color="0.3", ms=3.5))
l.append("random loss, same counts (open)")
axB.legend(h, l, frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
axB.set_title("B", loc="left", fontweight="bold")
axB.set_title("4P axons, Nov 2025 vs. complete", loc="right", fontsize=7, color="0.3")

out = ROOT / "figures"
for ext in ("pdf", "png"):
    fig.savefig(out / f"14_celltype_history.{ext}", dpi=300)
print("wrote figures/14_celltype_history.{pdf,png}")
