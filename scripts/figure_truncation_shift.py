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

Open diamonds: the same quantity against the distance-matched shuffle (Session 8; targets permuted only
within the same 50 um soma-distance bin), i.e. the part of the shift that is not spatial sampling.

Third column (Session 9, if the per-cell-type outputs exist), Nsyn only: the like-to-like slope beyond the
distance-matched shuffle (real - shuffled mean; bars = 95% range of the shuffle) for 4P and 23P presynaptic
axons separately, on complete axons and at radial f = 0.3.
(Session 7's third column, the calibrated branch-cut prediction, moved to the supplement: the real
proofreading history (Figure 1) replaced it as the model of realistic incompleteness.)
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
DIST_COL = "#8c5fc7"  # distance-matched shuffle, as in Figure 1
CTS = ("4P", "23P")
COND = {  # validated categorical slots 1-3 (dataviz palette, light surface); markers/dashes as secondary encoding
    "radial": ("#2a78d6", "o", "-", "radial truncation"),
    "blocks": ("#eb6834", "s", "--", "random cable loss"),
    "synsub": ("#1baf7a", "^", ":", "random synapse subsample"),
}
INK, MUTED, GRID, SURFACE = "#1f1f1e", "#6b6a63", "#e4e3dd", "#fcfcfb"


def load_celltype():
    """Real minus distance-matched shuffle, per presynaptic cell type, at f = 1 and 0.3 (Session 9)."""
    files = [f"trunc_slopes_poisson_{ct}.csv" for ct in CTS] + [f"trunc_perm_poisson_{ct}_nulldistance.csv" for ct in CTS]
    if not all((TD / f).exists() for f in files):
        return None
    rows = []
    for ct in CTS:
        f = pd.read_csv(TD / f"trunc_slopes_poisson_{ct}.csv")
        f = f[(f["proj"] == PROJ) & (f["model"] == "Nsyn") & (f["cond"] == "radial")]
        f["level"] = f["level"].round(2)
        pm = pd.read_csv(TD / f"trunc_perm_poisson_{ct}_nulldistance.csv")
        pm = pm[(pm["proj"] == PROJ) & (pm["model"] == "Nsyn")]
        pm["level"] = pm["level"].round(2)
        for sim in SIMS:
            for lv in (1.0, 0.3):
                r = f[(f["sim"] == sim) & (f["level"] == lv)]["slope"].iloc[0]
                v = pm[(pm["sim"] == sim) & (pm["level"] == lv)]["slope"].dropna().to_numpy()
                mu = v.mean()
                p = (1 + (abs(v - mu) >= abs(r - mu)).sum()) / (1 + len(v))
                rows.append(dict(cell_type=ct, sim=sim, level=lv, real=r, null=mu, gap=r - mu,
                                 lo=r - pd.Series(v).quantile(.975), hi=r - pd.Series(v).quantile(.025), p=p, n=len(v)))
    return pd.DataFrame(rows)


def celltype_panel(ax, ctab, sim, top):
    ax.axhline(0, color=MUTED, lw=.8, zorder=1)
    for k, ct in enumerate(CTS):
        for m, (lv, mk, fc) in enumerate(((1.0, "o", DIST_COL), (0.3, "o", "white"))):
            r = ctab[(ctab["cell_type"] == ct) & (ctab["sim"] == sim) & (ctab["level"] == lv)].iloc[0]
            x = k + (m - .5) * .36
            ax.errorbar([x], [r["gap"]], yerr=[[r["gap"] - r["lo"]], [r["hi"] - r["gap"]]], fmt=mk, color=DIST_COL,
                        mfc=fc, ms=5.5, capsize=2.5, lw=1.2, zorder=3,
                        label="complete axons" if lv == 1 else "radial f = 0.3")
    ax.set_xticks(range(len(CTS)), [f"{ct} axons" for ct in CTS])
    ax.set_xlim(-.6, len(CTS) - .4)
    ax.grid(axis="y", color=GRID, lw=.6, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("real − distance-matched shuffle", fontsize=8, color=INK)
    if top:
        ax.set_title("Beyond distance, by cell type\n(Nsyn)", fontsize=9.5, color=INK)


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
    gapd = None
    if (TD / "trunc_perm_poisson_nulldistance.csv").exists():
        pd_ = pd.read_csv(TD / "trunc_perm_poisson_nulldistance.csv")
        pd_["level"] = pd_["level"].round(2)
        pd_ = pd_[pd_["proj"] == PROJ]
        gapd = pd_.groupby(key + ["level"])["slope"].mean().rename("null").reset_index().merge(real, on=key + ["level"])
        gapd["gap"] = gapd["slope"] - gapd["null"]
        gapd = gapd.join(gapd[gapd["level"] == 1].set_index(key)["gap"].rename("gap1"), on=key)
        gapd["gap_shift"] = gapd["gap"] - gapd["gap1"]

    # ------------------------------------------------------------ figure
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "font.family": "DejaVu Sans"})
    ctab = load_celltype()
    ncol = 3 if ctab is not None else 2
    fig, axes = plt.subplots(2, ncol, figsize=(3.6 * ncol, 6.0), squeeze=False,
                             gridspec_kw={"width_ratios": [1, 1, .62][:ncol]})
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
                    markeredgewidth=.6, zorder=4, label="radial, beyond cell-type shuffle")
            if gapd is not None:
                g = gapd[(gapd["sim"] == sim) & (gapd["model"] == model)].sort_values("level")
                ax.plot(g["level"], g["gap_shift"], ls="none", marker="D", ms=5.5, color=DIST_COL, mfc="white",
                        markeredgewidth=1.3, zorder=5, label="radial, beyond distance-matched shuffle")
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
    if ctab is not None:
        for i, sim in enumerate(SIMS):
            celltype_panel(axes[i, 2], ctab, sim, top=i == 0)
    h, lab = [], []
    for a in axes.ravel():
        for hh, ll in zip(*a.get_legend_handles_labels()):
            if ll not in lab:
                h.append(hh), lab.append(ll)
    fig.legend(h, lab, loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(.5, -.005))
    fig.suptitle("Change in the V1→V1 like-to-like slope as axons are made incomplete", fontsize=10, color=INK,
                 y=.995)
    fig.text(.5, .952, "below 0 = rule looks stronger; bands: 95% paired bootstrap over presynaptic neurons",
             ha="center", fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, .1, 1, .95))
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"07_truncation_shift.{ext}", dpi=300, facecolor=SURFACE)

    # ------------------------------------------------------------ table
    t = est[est["level"].isin([0.7, 0.5, 0.3])].pivot_table(index=key + ["level"], columns="cond",
                                                             values=["shift", "lo", "hi"])
    rows = []
    for (sim, model, level), r in t.iterrows():
        cell = {c: f"{r[('shift', c)]:+.3f} [{r[('lo', c)]:+.3f}, {r[('hi', c)]:+.3f}]" for c in COND}
        gs = gap[(gap["sim"] == sim) & (gap["model"] == model) & (gap["level"] == level)]["gap_shift"]
        gd = gapd[(gapd["sim"] == sim) & (gapd["model"] == model) & (gapd["level"] == level)]["gap_shift"] if gapd is not None else []
        rows.append({"similarity": SIMS[sim], "model": MODELS[model], "f": level, **cell,
                     "radial beyond cell-type shuffle": f"{gs.iloc[0]:+.3f}" if len(gs) else "",
                     "radial beyond distance-matched shuffle": f"{gd.iloc[0]:+.3f}" if len(gd) else ""})
    tab = pd.DataFrame(rows).sort_values(["similarity", "model", "f"], ascending=[False, True, False])
    (ROOT / "notebooks" / "07_shift_table.md").write_text(
        f"# Slope shift vs. completeness (auto-generated by scripts/figure_truncation_shift.py, v{VERSION})\n\n"
        f"{PROJ}; slope(f) - slope(1) per 10°, log scale; full-data estimate [95% paired bootstrap, "
        f"{boot['boot'].nunique()} samples]. 'beyond ... shuffle' = change in (real - shuffled mean) since f = 1.\n\n"
        + tab.to_markdown(index=False) + "\n"
        + ("" if ctab is None else "\n## Beyond the distance-matched shuffle, by presynaptic cell type (Nsyn)\n\n"
           + ctab.assign(similarity=ctab["sim"].map(SIMS)).drop(columns="sim").round(3).to_markdown(index=False) + "\n"))
    print("saved figures/07_truncation_shift.png/.pdf, notebooks/07_shift_table.md")
    print(tab.to_string(index=False))


if __name__ == "__main__":
    main()
