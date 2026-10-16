"""Session 8 preprint figure: real proofreading history vs. the like-to-like rule.

A  completeness of the 193 fully extended V1 axons over time (monthly chunkedgraph states)
B  RF like-to-like (Nsyn) shift vs. the complete state over time; paired bootstrap CIs at 3 states
C  same for orientation
D/E  real minus shuffled-wiring slope across states: distance-matched vs. cell-type shuffle (RF, orientation)

Inputs (data/truncation/, data/history/):
  trunc_slopes_poisson_hist.csv, trunc_boot_poisson_hist.csv,
  trunc_perm_poisson_hist{,_nulldistance}.csv, trunc_perm_poisson{,_nulldistance}.csv,
  trunc_perm_poisson_cyl{,_nulldistance}.csv, trunc_slopes_poisson.csv,
  trunc_slopes_poisson_cyl.csv, history_axons_v1822.parquet
Outputs: figures/08_history.png/.pdf, notebooks/08_figure_table.md
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--trunc", default="data/truncation")
ap.add_argument("--hist", default="data/history")
ap.add_argument("--fig", default="figures")
ap.add_argument("--nb", default="notebooks")
a = ap.parse_args()
T, H = Path(a.trunc), Path(a.hist)

PROJ = "V1->V1"
PIN = 1822
PIN_DATE = pd.Timestamp("2026-07-01")  # plotted position of v1822 (after the 2026-06 state)
SIMS = {"rf_dist_deg": "RF distance", "ori_diff_deg": "orientation difference"}


def v1(df, model="Nsyn"):
    return df[(df.proj == PROJ) & (df.model == model)]


def to_date(level):
    level = int(level)
    if level == PIN:
        return PIN_DATE
    return pd.Timestamp(year=level // 100, month=level % 100, day=1)


# ---------------------------------------------------------------- A: completeness
ax_h = pd.read_parquet(H / "history_axons_v1822.parquet")
g = ax_h.groupby("version")
comp = pd.DataFrame({
    "syn_frac": g.l2l_syn_attached.sum() / g.l2l_syn_full.sum(),
    "cable_med": g.cable_attached.median(),
    "n_complete": g.cable_attached.apply(lambda x: (x > 0.95).mean()),
}).reset_index()
comp["date"] = comp.version.map(to_date)
comp = comp.sort_values("date")

# ---------------------------------------------------------------- B/C: shift over time
s = v1(pd.read_csv(T / "trunc_slopes_poisson_hist.csv"))
real_hist = s[s.cond == "hist"].assign(date=lambda d: d.level.map(to_date))
full = s[(s.cond == "radial") & (s.level == 1)].set_index("sim").slope

b = v1(pd.read_csv(T / "trunc_boot_poisson_hist.csv"))
bref = b[b.level == 1].set_index(["boot", "sim"]).slope
bh = b[b.level != 1].copy()
bh["shift"] = bh.slope.values - bref.reindex(list(zip(bh.boot, bh.sim))).values
bh = bh[bh.level != PIN]
boot_ci = (bh.groupby(["level", "sim"])["shift"]
             .quantile([0.025, 0.975]).unstack().reset_index()
             .rename(columns={0.025: "lo", 0.975: "hi"}))
boot_ci["date"] = boot_ci.level.map(to_date)

# ---------------------------------------------------------------- D: real − null by state
def null_summary(path, level):
    d = v1(pd.read_csv(path))
    d = d[d.level == level]
    return d.groupby("sim").slope.agg(
        mean="mean", lo=lambda x: x.quantile(0.025), hi=lambda x: x.quantile(0.975),
        vals=lambda x: list(x))


real_lookup = {}
for lv in (202304, 202511, 202602):
    real_lookup[f"hist_{lv}"] = real_hist[real_hist.level == lv].set_index("sim").slope
real_lookup["complete"] = full
rad = v1(pd.read_csv(T / "trunc_slopes_poisson.csv"))
rad = rad[rad.cond == "radial"]
for f in (0.5, 0.3):
    real_lookup[f"radial_{f}"] = rad[np.isclose(rad.level, f)].groupby("sim").slope.mean()
cyl_f = T / "trunc_slopes_poisson_cyl.csv"  # Session 11: column cylinder, Rc = 100 um
if cyl_f.exists():
    cyl = v1(pd.read_csv(cyl_f))
    real_lookup["cyl_100"] = cyl[(cyl.cond == "cyl") & np.isclose(cyl.level, 100)].set_index("sim").slope

states = [  # key, label, distance-null (file, level), cell-null (file, level)
    ("complete", "complete\nv1822", (T / "trunc_perm_poisson_nulldistance.csv", 1.0), (T / "trunc_perm_poisson.csv", 1.0)),
    ("hist_202304", "Apr\n2023", (T / "trunc_perm_poisson_hist_nulldistance.csv", 202304), (T / "trunc_perm_poisson_hist.csv", 202304)),
    ("hist_202511", "Nov\n2025", (T / "trunc_perm_poisson_hist_nulldistance.csv", 202511), (T / "trunc_perm_poisson_hist.csv", 202511)),
    ("hist_202602", "Feb\n2026", (T / "trunc_perm_poisson_hist_nulldistance.csv", 202602), (T / "trunc_perm_poisson_hist.csv", 202602)),
    ("cyl_100", "column\nRc = 100 µm", (T / "trunc_perm_poisson_cyl_nulldistance.csv", 100), (T / "trunc_perm_poisson_cyl.csv", 100)),
    ("radial_0.3", "radial\nf = 0.3", (T / "trunc_perm_poisson_nulldistance.csv", 0.3), (T / "trunc_perm_poisson.csv", 0.3)),
]
# drop states whose shuffles haven't been run yet (Feb 2026 is added by scripts/run_session10.sh)
states = [st for st in states if st[2][0].exists() and (pd.read_csv(st[2][0], usecols=["level"]).level.round(2) == round(st[2][1], 2)).any()]
rows = []
for key, label, dn, cn in states:
    for nullname, (path, lv) in (("distance", dn), ("cell", cn)):
        if not path.exists():
            continue
        ns = null_summary(path, lv)
        for sim in SIMS:
            if sim not in ns.index or sim not in real_lookup[key].index:
                continue
            r = real_lookup[key][sim]
            vals = np.asarray(ns.loc[sim, "vals"])
            mu = vals.mean()  # two-sided, centred on the shuffled mean (project convention)
            p = (1 + np.sum(np.abs(vals - mu) >= abs(r - mu))) / (1 + len(vals))
            rows.append(dict(state=key, label=label, null=nullname, sim=sim, real=r,
                             null_mean=ns.loc[sim, "mean"], gap=r - ns.loc[sim, "mean"],
                             gap_lo=r - ns.loc[sim, "hi"], gap_hi=r - ns.loc[sim, "lo"],
                             p=p, n=len(vals)))
gaps = pd.DataFrame(rows)

# ---------------------------------------------------------------- plot
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
import matplotlib.dates as mdates
fig = plt.figure(figsize=(7.2, 5.4), constrained_layout=True)
gs = fig.add_gridspec(2, 6)
axA, axB, axC = fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6])
axD, axE = fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])
plateau = (pd.Timestamp("2023-04-01"), pd.Timestamp("2025-11-01"))
XLIM = (pd.Timestamp("2022-06-01"), pd.Timestamp("2026-09-01"))


def time_axis(ax):
    ax.axvspan(*plateau, color="0.92", lw=0)
    ax.set_xlim(*XLIM)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))


ax = axA
time_axis(ax)
ax.plot(comp.date, comp.syn_frac, color="k", label="synapses (all targets)")
ax.plot(comp.date, comp.cable_med, color="tab:blue", label="median axon cable")
ax.plot(comp.date, comp.n_complete, color="tab:green", label="axons complete")
ax.set_ylim(0, 1.03); ax.set_ylabel("fraction attached to soma")
ax.set_title("A  Proofreading history", loc="left")
ax.legend(frameon=False, loc="upper left", fontsize=6.5)
ax.text(plateau[0] + (plateau[1] - plateau[0]) / 2, 0.03, "partial plateau", ha="center", color="0.4", fontsize=6.5)

for ax, sim, ttl in ((axB, "rf_dist_deg", "B  RF rule"),
                     (axC, "ori_diff_deg", "C  Orientation rule")):
    time_axis(ax)
    rh = real_hist[real_hist.sim == sim].sort_values("date")
    rh = rh[rh.date >= XLIM[0]]  # earlier states have < 20% of synapses attached
    ax.plot(rh.date, rh.slope - full[sim], color="tab:red", lw=1.1)
    pt = rh.set_index("level").slope - full[sim]
    for _, r in boot_ci[boot_ci.sim == sim].iterrows():
        y = pt.get(r.level, np.nan)
        ax.errorbar(r.date, y, yerr=[[y - r.lo], [r.hi - y]], fmt="o", color="tab:red", ms=3, capsize=2, lw=0.9)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_title(ttl, loc="left")
    # complete-axon slopes (RF −0.856, orientation −0.040) are given in the caption
axB.set_ylabel("slope shift vs. complete axons\n(< 0 = rule looks stronger)")

order = [k for k, *_ in states]
xs = np.arange(len(order))
for ax, sim, ttl in ((axD, "rf_dist_deg", "D  RF rule beyond shuffled wiring"),
                     (axE, "ori_diff_deg", "E  Orientation rule beyond shuffled wiring")):
    for nullname, col, off, mfc in (("distance", "tab:purple", -0.12, "tab:purple"), ("cell", "0.55", 0.12, "white")):
        d = gaps[(gaps.sim == sim) & (gaps.null == nullname)].set_index("state").reindex(order)
        ok = d.gap.notna().values
        ax.errorbar(xs[ok] + off, d.gap[ok], yerr=[(d.gap - d.gap_lo)[ok], (d.gap_hi - d.gap)[ok]],
                    fmt="o", color=col, mfc=mfc, ms=3.5, capsize=1.5, lw=0.9,
                    label={"distance": "distance-matched shuffle", "cell": "cell-type shuffle"}[nullname])
        for x, (_, r) in zip(xs[ok], d[ok].iterrows()):
            if r.p > 0.05:
                ax.text(x + off, r.gap_hi, "ns", ha="center", va="bottom", fontsize=6, color=col)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(xs); ax.set_xticklabels([l for _, l, *_ in states], fontsize=6.5)
    ax.set_title(ttl, loc="left")
    ax.legend(frameon=False, fontsize=6.5, loc="lower left")
    lo = gaps[gaps.sim == sim].gap_lo.min()
    ax.set_ylim(lo * 1.25, max(0.0, gaps[gaps.sim == sim].gap_hi.max()) + abs(lo) * 0.12)
axD.set_ylabel("real − shuffled slope\n(bars: shuffle 95% range)")
fig.supxlabel("V1→V1, synapse count (Nsyn) model, slope per 10° of RF distance / orientation difference",
              fontsize=6.5, color="0.35")

Path(a.fig).mkdir(exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(Path(a.fig) / f"08_history.{ext}", dpi=300)

# ---------------------------------------------------------------- table
lines = ["# Figure 08 values (V1→V1, Nsyn, slope per 10°)", "",
         "## Shift vs. complete axons (paired bootstrap, 200 samples over presynaptic neurons)", "",
         "| state | sim | real | shift | 95% CI |", "|---|---|---|---|---|"]
for _, r in boot_ci.sort_values(["sim", "level"]).iterrows():
    real = real_hist[(real_hist.level == r.level) & (real_hist.sim == r.sim)].slope
    real = real.iloc[0] if len(real) else full[r.sim]
    lines.append(f"| {int(r.level)} | {SIMS[r.sim]} | {real:.3f} | {real - full[r.sim]:+.3f} | [{r.lo:+.3f}, {r.hi:+.3f}] |")
lines += ["", "## Real vs. shuffled wiring (two-sided permutation p centred on the shuffled mean, floor 1/(n+1))", "",
          "| state | sim | null | real | null mean | real − null | 95% range of gap | p | n |",
          "|---|---|---|---|---|---|---|---|---|"]
for _, r in gaps.iterrows():
    lines.append(f"| {r.state} | {SIMS[r.sim]} | {r.null} | {r.real:.3f} | {r.null_mean:.3f} | {r.gap:+.3f} | "
                 f"[{r.gap_lo:+.3f}, {r.gap_hi:+.3f}] | {r.p:.3f} | {r.n} |")
lines += ["", "## Completeness by version", "", "| version | synapses attached | median cable | axons > 95% |", "|---|---|---|---|"]
for _, r in comp.iterrows():
    lines.append(f"| {int(r.version)} | {r.syn_frac:.2f} | {r.cable_med:.2f} | {r.n_complete:.2f} |")
Path(a.nb, "08_figure_table.md").write_text("\n".join(lines) + "\n")
print(gaps.drop(columns=["label"]).to_string())
