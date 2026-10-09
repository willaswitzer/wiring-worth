"""Session 7: preprint figure -- shift of the like-to-like slope with simulated axon completeness.

Run from the repo root after the truncation pipeline (needs trunc_boot_poisson.csv and trunc_perm_*.csv):
    .venv/bin/python scripts/figure_truncation_shift.py

Writes figures/07_truncation_shift.png (and .pdf) and notebooks/07_shift_table.md.

Each panel (V1->V1; columns = Nsyn, Nsyn/Ld; rows = RF distance, orientation) shows
slope(f) - slope(1), i.e. how far the like-to-like estimate moves when each axon keeps a fraction f
of its cable:
- radial truncation, random cable loss (blocks) and random synapse subsampling (synsub);
  lines = full-data estimates (nulls: mean over 5 reps), bands = 95% paired bootstrap over
  presynaptic neurons (200 samples; nulls use rep 0);
- diamonds: the shift of the real-minus-shuffled gap, (real - null mean)(f) - (real - null mean)(1),
  at the levels with a shuffled control (f = 1, 0.5, 0.3), i.e. the part of the shift that the
  shuffled-wiring control does not also show.
The dotted vertical line marks the cable ratio of real partially extended 23P axons (~0.75, Session 5).
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from config import VERSION  # noqa: E402

TD = ROOT / "data" / "truncation"
PROJ = "V1->V1"
PARTIAL_23P_CABLE = 0.75  # partial/full median axon cable, 23P (Session 5, 05_results.md)
SIMS = {"rf_dist_deg": "RF-centre distance", "ori_diff_deg": "Orientation difference"}
MODELS = {"Nsyn": "Synapses (Nsyn)", "density": "Synapses per contact (Nsyn / Ld)"}
COND = {  # validated categorical slots 1-3 (dataviz palette, light surface); markers/dashes as secondary encoding
    "radial": ("#2a78d6", "o", "-", "radial truncation"),
    "blocks": ("#eb6834", "s", "--", "random cable loss"),
    "synsub": ("#1baf7a", "^", ":", "random synapse subsample"),
}
INK, MUTED, GRID, SURFACE = "#1f1f1e", "#6b6a63", "#e4e3dd", "#fcfcfb"


def main():
    fit = pd.read_csv(TD / "trunc_slopes_poisson.csv")
    boot = pd.read_csv(TD / "trunc_boot_poisson.csv")
    perm = pd.read_csv(TD / "trunc_perm_poisson.csv")
    for d in (fit, boot, perm):
        d["level"] = d["level"].round(2)
    fit, boot, perm = (d[d["proj"] == PROJ] for d in (fit, boot, perm))
    key = ["sim", "model"]

    # point estimates: shift vs radial f = 1 (nulls: mean over reps)
    full = fit[fit["set"] == "radial_1.00_0"].set_index(key)["slope"].rename("full")
    est = fit.groupby(key + ["cond", "level"])["slope"].mean().reset_index().join(full, on=key)
    est["shift"] = est["slope"] - est["full"]
    f1 = est[(est["cond"] == "radial") & (est["level"] == 1)]
    est = pd.concat([est] + [f1.assign(cond=c) for c in ("blocks", "synsub")], ignore_index=True)

    # paired bootstrap CIs on the shift
    bfull = boot[boot["set"] == "radial_1.00_0"].set_index(["boot"] + key)["slope"].rename("bfull")
    boot = boot.join(bfull, on=["boot"] + key)
    boot["shift"] = boot["slope"] - boot["bfull"]
    ci = boot.groupby(key + ["cond", "level"])["shift"].quantile([.025, .975]).unstack()
    ci.columns = ["lo", "hi"]
    ci = ci.reset_index()
    ci1 = ci[(ci["cond"] == "radial") & (ci["level"] == 1)]
    ci = pd.concat([ci] + [ci1.assign(cond=c) for c in ("blocks", "synsub")], ignore_index=True)
    est = est.merge(ci, on=key + ["cond", "level"], how="left")

    # shift of the real-minus-shuffled gap
    null = perm.groupby(key + ["level"])["slope"].mean().rename("null").reset_index()
    real = est[est["cond"] == "radial"][key + ["level", "slope"]]
    gap = null.merge(real, on=key + ["level"])
    gap["gap"] = gap["slope"] - gap["null"]
    gap = gap.join(gap[gap["level"] == 1].set_index(key)["gap"].rename("gap1"), on=key)
    gap["gap_shift"] = gap["gap"] - gap["gap1"]

    # ------------------------------------------------------------ figure
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "font.family": "DejaVu Sans"})
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), squeeze=False)
    for i, sim in enumerate(SIMS):
        for j, model in enumerate(MODELS):
            ax = axes[i, j]
            ax.axhline(0, color=MUTED, lw=.8, zorder=1)
            ax.axvline(PARTIAL_23P_CABLE, color=MUTED, lw=.8, ls=(0, (1, 2)), zorder=1)
            for cond, (col, mk, ls, lab) in COND.items():
                c = est[(est["sim"] == sim) & (est["model"] == model) & (est["cond"] == cond)].sort_values("level")
                ax.fill_between(c["level"], c["lo"], c["hi"], color=col, alpha=.13, lw=0, zorder=2)
                ax.plot(c["level"], c["shift"], color=col, lw=2, ls=ls, marker=mk, ms=4.5, label=lab, zorder=3,
                        markeredgecolor="white", markeredgewidth=.6)
            g = gap[(gap["sim"] == sim) & (gap["model"] == model)].sort_values("level")
            ax.plot(g["level"], g["gap_shift"], ls="none", marker="D", ms=5.5, color=INK, markeredgecolor="white",
                    markeredgewidth=.6, zorder=4, label="radial, beyond shuffled wiring")
            ax.set_xlim(1.03, 0.27)
            ax.grid(axis="y", color=GRID, lw=.6, zorder=0)
            ax.spines[["top", "right"]].set_visible(False)
            if i == 0:
                ax.set_title(MODELS[model], fontsize=9.5, color=INK)
            if i == 1:
                ax.set_xlabel("fraction of each axon's cable kept")
            if j == 0:
                ax.set_ylabel(f"{SIMS[sim]}\nchange in slope per 10°")
            if i == 0 and j == 0:
                ax.annotate("real 23P\npartial axons", xy=(PARTIAL_23P_CABLE, ax.get_ylim()[1]), xytext=(3, -2),
                            textcoords="offset points", fontsize=7.5, color=MUTED, va="top")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=2, frameon=False, fontsize=8.5, bbox_to_anchor=(.5, -.005))
    fig.suptitle("Change in the V1→V1 like-to-like slope as axons are made incomplete", fontsize=10, color=INK,
                 y=.995)
    fig.text(.5, .952, "below 0 = rule looks stronger; bands: 95% paired bootstrap over presynaptic neurons",
             ha="center", fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, .07, 1, .95))
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"07_truncation_shift.{ext}", dpi=300, facecolor=SURFACE)

    # ------------------------------------------------------------ table
    t = est[est["level"].isin([0.7, 0.5, 0.3])].pivot_table(index=key + ["level"], columns="cond",
                                                             values=["shift", "lo", "hi"])
    rows = []
    for (sim, model, level), r in t.iterrows():
        cell = {c: f"{r[('shift', c)]:+.3f} [{r[('lo', c)]:+.3f}, {r[('hi', c)]:+.3f}]" for c in COND}
        gs = gap[(gap["sim"] == sim) & (gap["model"] == model) & (gap["level"] == level)]["gap_shift"]
        rows.append({"similarity": SIMS[sim], "model": MODELS[model], "f": level, **cell,
                     "radial beyond shuffle": f"{gs.iloc[0]:+.3f}" if len(gs) else ""})
    tab = pd.DataFrame(rows).sort_values(["similarity", "model", "f"], ascending=[False, True, False])
    (ROOT / "notebooks" / "07_shift_table.md").write_text(
        f"# Slope shift vs. completeness (auto-generated by scripts/figure_truncation_shift.py, v{VERSION})\n\n"
        f"{PROJ}; slope(f) - slope(1) per 10°, log scale; full-data estimate [95% paired bootstrap, "
        f"{boot['boot'].nunique()} samples]. 'radial beyond shuffle' = change in (real - shuffled mean) since f = 1.\n\n"
        + tab.to_markdown(index=False) + "\n")
    print("saved figures/07_truncation_shift.png/.pdf, notebooks/07_shift_table.md")
    print(tab.to_string(index=False))


if __name__ == "__main__":
    main()
