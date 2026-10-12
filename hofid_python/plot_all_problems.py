"""
plot_all_problems.py -- generate PNG figures from the computed benchmark results.

    python plot_all_problems.py          # uses results/tables/all_results.csv + results/data/*.npz
    python plot_all_problems.py --run    # (re)run run_all_problems.py first
    python plot_all_problems.py --run --quick

Every curve is plotted from data computed by this code; nothing is copied from the paper.
Figures are Python re-creations of the paper's Figs. 2-7 (they are NOT the paper's images).
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = "results"
FIG = f"{OUT}/figures"
# MATLAB-like style, as in the figures of the paper: boxed axes, inward ticks on all sides, no grid, no title
plt.rcParams.update({"font.size": 10, "axes.grid": False, "figure.dpi": 100, "font.family": "DejaVu Sans",
                     "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
                     "legend.fontsize": 9, "legend.framealpha": 1.0, "legend.edgecolor": "black", "legend.fancybox": False,
                     "axes.linewidth": 0.8})
MBLUE, MRED, MMAG = "#0072BD", "#D95319", "#FF00FF"
COL = {4: "tab:green", 6: "tab:blue", 8: "tab:red", 10: "tab:purple", "var": "tab:orange"}


def load_runs():
    runs = []
    df = pd.read_csv(f"{OUT}/tables/all_results.csv")
    for _, r in df.iterrows():
        f = r.get("data_file")
        if isinstance(f, str) and os.path.exists(f):
            d = np.load(f, allow_pickle=False)
            runs.append(dict(row=r, t=d["t"], Y=d["Y"], t0=d["t_initial"], err=d["err"]))
    return runs


def sel(runs, **kw):
    out = []
    for r in runs:
        row = r["row"]
        if all((str(row.get(k)) == str(v)) for k, v in kw.items()):
            out.append(r)
    return out


def okrun(r):
    return bool(r["row"]["converged"])


def label(r):
    o = str(r["row"]["orders"])
    return f"order {o}" if "/" not in o else f"orders {o}"


def key(r):
    o = str(r["row"]["orders"])
    return "var" if "/" in o else int(o)


def paper_mesh(ax, r, show_legend=True):
    """Fig. 3 style: red circles = starting (equidistant) mesh (upper), blue circles = adapted mesh (lower)."""
    t0, t = r["t0"], r["t"]
    ax.plot(t, np.zeros_like(t), "o", mfc="none", mec="#0000FF", ms=5, mew=1.6, label="adapted mesh")
    ax.plot(t0, np.ones_like(t0), "o", mfc="none", mec="#FF2200", ms=5, mew=1.6, label="equidistant mesh")
    ax.set_ylim(-0.3, 1.3); ax.set_yticks([]); ax.set_xlim(t0[0] - 0.0, t0[-1] * 1.0)
    if show_legend: ax.legend(loc="upper right", frameon=True)


def fig_problem1(runs):
    rs = sel(runs, problem=1)
    ok = [r for r in rs if okrun(r)]
    if not ok: return
    r = [x for x in ok if key(x) == 8] or ok
    r = r[0]
    # Fig. 2 style (computed solution of Problem 1; exact solution is t^2 ln t)
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.plot(r["t"], r["Y"], "-", color=MBLUE, lw=2.5)
    ax.set_xlim(0, 1); ax.set_ylim(-0.2, 0); ax.set_xticks(np.arange(0, 1.01, 0.1)); ax.set_yticks(np.arange(-0.2, 0.001, 0.02))
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_1_solution.png", dpi=250); plt.close(fig)
    # true error (extra, not a paper figure)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for x in ok:
        ex = x["t"] ** 2 * np.log(np.maximum(x["t"], 1e-300)); ex[0] = 0
        ax.semilogy(x["t"], np.abs(x["Y"] - ex) + 1e-20, lw=1, label=f"order {key(x)}")
    ax.set_xlabel("t"); ax.set_ylabel("|Y - y_exact|"); ax.legend(); ax.set_xlim(0, 1)
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_1_error.png", dpi=250); plt.close(fig)
    # Fig. 3 style
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    paper_mesh(ax, r); ax.set_xticks(np.arange(0, 1.01, 0.1)); ax.set_xlim(0, 1)
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_1_mesh.png", dpi=250); plt.close(fig)


def fig_problem2(runs):
    rs = sel(runs, problem=2, orders=6)
    if not rs: return
    sty = {"y1": (MBLUE, ".", "-"), "y2": ("#FF0000", "s", "-"), "y3": (MMAG, "o", "-")}
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for r in rs:                       # draw y3 first? keep order y1,y2,y3 as in the paper's legend-less Fig. 4
        b = r["row"]["branch"]
        if okrun(r) and b in sty:
            c, m, ls = sty[b]
            ax.plot(r["t"], r["Y"], ls, color=c, marker=m, ms=3.5 if m == "." else 5, mfc="none" if m != "." else c, lw=1)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks(np.arange(0, 1.01, 0.2)); ax.set_yticks(np.arange(0, 1.01, 0.1))
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_2_branches.png", dpi=250); plt.close(fig)
    fig, axs = plt.subplots(3, 1, figsize=(6.4, 9))
    for i, b in enumerate(("y1", "y2", "y3")):
        r = [x for x in rs if x["row"]["branch"] == b and okrun(x)]
        if r:
            paper_mesh(axs[i], r[0]); axs[i].set_xlim(0, 1)
            axs[i].set_xlabel(f"({'abc'[i]}) solution {b}(t), order 6")
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_2_mesh.png", dpi=250); plt.close(fig)


def _L(r): return float(r["row"]["config"].split("=")[1].split(",")[0])


def fig_problem3(runs):
    rs = sel(runs, problem=3, orders=6)
    if not rs: return
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
    caps = {"robin": "(a) $y'(0)=0,\\quad y(L)+y'(L)=0$", "dirichlet": "(b) $y'(0)=0,\\quad y(L)=0$"}
    for ax, bc in zip(axs, ("robin", "dirichlet")):
        for L, c, m in ((20, MBLUE, "."), (5, "#FF0000", "+")):
            r = [x for x in rs if x["row"]["bc_variant"] == bc and _L(x) == L and okrun(x)]
            if r:
                ax.plot(r[0]["t"], r[0]["Y"], "-", color=c, marker=m, ms=4 if m == "." else 6, lw=1.2, label=f"L = {L}")
        ax.set_xlim(0, 20); ax.set_ylim(0, 2.5); ax.set_xticks([0, 5, 10, 15, 20]); ax.set_yticks([0, .5, 1, 1.5, 2, 2.5])
        ax.legend(loc="upper right"); ax.set_xlabel(caps[bc], fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_3_solutions.png", dpi=250); plt.close(fig)
    fig, axs = plt.subplots(2, 1, figsize=(6.6, 6.5))
    for i, L in enumerate((5, 20)):
        r = [x for x in rs if x["row"]["bc_variant"] == "robin" and _L(x) == L and okrun(x)]
        if r:
            paper_mesh(axs[i], r[0]); axs[i].set_xlim(0, L); axs[i].set_xlabel(f"L = {L}, $y(L)+y'(L)=0$, order 6")
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_3_mesh.png", dpi=250); plt.close(fig)


def fig_problem4(runs):
    allp4 = sel(runs, problem=4)
    sols = {}
    for xi in (0.1, 0.6):
        rs = [x for x in allp4 if x["row"]["config"].startswith(f"xi={xi},") and okrun(x) and str(x["row"]["orders"]) == "8"
              and "extra" not in x["row"]["config"]]
        if rs: sols[xi] = rs[0]
    for xi, fn in ((0.1, "problem_4_xi_0_1.png"), (0.6, "problem_4_xi_0_6.png")):   # Fig. 6 style
        if xi not in sols: continue
        r = sols[xi]
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
        ax.plot(r["t"], r["Y"], ".", color="#0000FF", ms=4)
        ax.set_xlim(0, 1); ax.set_xticks([0, .2, .4, .6, .8, 1])
        fig.tight_layout(); fig.savefig(f"{FIG}/{fn}", dpi=250); plt.close(fig)
    if sols:                                                                      # Fig. 7 style
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        for xi, c in ((0.1, MBLUE), (0.6, MRED)):
            if xi in sols: ax.plot(sols[xi]["t"], sols[xi]["Y"], "-", color=c, lw=1.6, label=f"$\\xi$ = {xi}")
        ax.set_xlim(0, 1); ax.set_ylim(-1, 0.6); ax.set_xticks([0, .2, .4, .6, .8, 1]); ax.legend(loc="upper right")
        fig.tight_layout(); fig.savefig(f"{FIG}/problem_4_xi_comparison.png", dpi=250); plt.close(fig)
    fig, axs = plt.subplots(2, 1, figsize=(6.6, 6.5))
    for i, xi in enumerate((0.1, 0.6)):
        if xi in sols:
            paper_mesh(axs[i], sols[xi]); axs[i].set_xlim(0, 1); axs[i].set_xlabel(f"$\\xi$ = {xi}, order 8")
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_4_mesh.png", dpi=250); plt.close(fig)


def fig_summary():
    df = pd.read_csv(f"{OUT}/tables/all_results.csv")
    df = df[df.converged.astype(bool) & ~df.config.fillna("").str.contains("footnote|extra")].copy()
    df["ord"] = df.orders.astype(str).map(lambda s: "var" if "/" in s else s)
    df = df[df["ord"] != "var"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for pr, br, lab in [(1, "", "Problem 1"), (2, "y1", "Problem 2 y1"), (2, "y2", "Problem 2 y2"), (2, "y3", "Problem 2 y3"),
                        (3, "", "Problem 3 L=5 Robin"), (4, "", "Problem 4 xi=0.1")]:
        d = df[df.problem == pr]
        if pr == 2: d = d[d.branch == br]
        if pr == 3: d = d[(d.bc_variant == "robin") & d.config.str.contains("L=5,")]
        if pr == 4: d = d[d.config.str.startswith("xi=0.1,")]
        d = d.sort_values("ord", key=lambda s: s.astype(int))
        if d.empty: continue
        ax[0].semilogy(d["ord"].astype(int), d.final_mesh_points, "o-", label=lab, ms=4)
        ax[1].semilogy(d["ord"].astype(int), d.est_abs_err_max, "o-", label=lab, ms=4)
    d1 = df[df.problem == 1].sort_values("ord", key=lambda s: s.astype(int))
    ax[1].semilogy(d1["ord"].astype(int), d1.actual_abs_err, "s--", color="k", label="Problem 1 actual error", ms=4)
    ax[1].axhline(1e-8, color="gray", ls=":"); ax[0].set_xlabel("order p"); ax[1].set_xlabel("order p")
    ax[0].set_ylabel("final mesh points (TOL = 1e-8)"); ax[1].set_ylabel("estimated abs. error (dotted: TOL)")
    ax[0].legend(); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{FIG}/convergence_summary.png", dpi=250); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="run the benchmark first")
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args(argv)
    if a.run or not os.path.exists(f"{OUT}/tables/all_results.csv"):
        if not a.run:
            print("results/tables/all_results.csv not found -> running run_all_problems.py (use --run to force)")
        import run_all_problems
        run_all_problems.main(["--quick"] if a.quick else [])
    os.makedirs(FIG, exist_ok=True)
    runs = load_runs()
    for fn in (fig_problem1, fig_problem2, fig_problem3, fig_problem4):
        fn(runs)
    fig_summary()
    print("figures written to", FIG, ":", sorted(os.listdir(FIG)))


if __name__ == "__main__":
    main()
