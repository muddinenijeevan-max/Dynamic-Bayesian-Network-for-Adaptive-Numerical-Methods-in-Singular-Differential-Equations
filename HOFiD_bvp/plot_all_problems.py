"""
plot_all_problems.py -- interactive plotter for the HOFiD_bvp benchmarks.

    python plot_all_problems.py

Same two questions as run_all_problems.py (which problem, how many initial
mesh points).  The answers select the same per-run folder, e.g.

    results/problem_1_npts_15/

For every experiment of that selection:
    * HIT   : results/<run>/data/<experiment>.npz + _summary.json already exist
              (and were made with the same start mesh)  ->  they are loaded.
    * MISS  : the problem definition is IMPORTED from run_all_problems.py, the
              solver is run through solve_and_report(), the files are saved in
              the run folder, and the fresh result is plotted.

Figures are written to results/<run>/figures/.
No numerical code and no problem definitions live in this file.
"""
from __future__ import annotations

import json
import os
from types import SimpleNamespace

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# everything about problems / solving / folders comes from the runner
from run_all_problems import (
    get_experiments, ask_problem, ask_npts,
    run_dir_for, write_inputs_json, solve_and_report,
)

# MATLAB-like style
plt.rcParams.update({
    "font.size": 10, "axes.grid": False, "figure.dpi": 100,
    "font.family": "DejaVu Sans",
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "legend.fontsize": 9, "legend.framealpha": 1.0,
    "legend.edgecolor": "black", "legend.fancybox": False,
    "axes.linewidth": 0.8,
})
BLUE, RED, MAG = "#0072BD", "#D95319", "#FF00FF"


# ########################################################################### #
#   PART A.  CACHE LOOKUP  (hit -> load,  miss -> import runner + solve)      #
# ########################################################################### #
def _load(e, run_dir):
    """Try to load one experiment from run_dir.

    Returns  (result_like, diag_like, summary)  on a hit,
             "FAILED_BEFORE"                    if a previous run stored a failure,
             None                               on a miss / stale / unreadable file.
    """
    slug = e["slug"]
    npz_path = os.path.join(run_dir, "data", f"{slug}.npz")
    json_path = os.path.join(run_dir, "data", f"{slug}_summary.json")
    if not os.path.isfile(json_path):
        return None
    try:
        with open(json_path, "r", encoding="utf-8") as fh:
            summary = json.load(fh)
        if int(summary.get("initial_npts", -1)) != int(e["npts"]):
            return None                         # made with another start mesh
        if not os.path.isfile(npz_path):
            return "FAILED_BEFORE" if not summary.get("converged", False) else None
        with np.load(npz_path) as z:
            t = np.asarray(z["t"], float).copy()
            Y = np.asarray(z["Y"], float).copy()
            t_initial = np.asarray(z["t_initial"], float).copy()
            err = np.asarray(z["err"], float).copy()
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return None

    r = SimpleNamespace(t=t, Y=Y, t_initial=t_initial, err_node=err,
                        converged=bool(summary.get("converged", False)),
                        p_final=int(summary.get("fd_order", e["orders"][0])),
                        tol_final=float(summary.get("tolerance_requested", float("nan"))))
    exact = e["prob"].exact
    dg = {"exact": np.asarray(exact(t), float) if exact is not None else None}
    return r, dg, summary


def get_result(e, run_dir):
    """Load from the run folder, or import-and-solve through the runner."""
    got = _load(e, run_dir)
    if got == "FAILED_BEFORE":
        print("  [CACHE HIT ] earlier run did not converge -> nothing to plot.")
        return None
    if got is not None:
        print("  [CACHE HIT ] loaded saved solution.")
        return got
    print("  [CACHE MISS] not found in run folder -> solving with run_all_problems...")
    solve_and_report(e, out_dir=run_dir, show_table=False, print_progress=False)
    got = _load(e, run_dir)
    if got is None or got == "FAILED_BEFORE":
        print("  Solver did not converge -> skipped in plots.")
        return None
    return got


# ########################################################################### #
#   PART B.  PLOT FUNCTIONS                                                   #
# ########################################################################### #
def _mesh_plot(ax, t_init, t_final):
    ax.plot(t_final, np.zeros_like(t_final), "o",
            mfc="none", mec="#0000FF", ms=4, mew=1.4, label="adapted mesh")
    ax.plot(t_init, np.ones_like(t_init), "o",
            mfc="none", mec="#FF2200", ms=4, mew=1.4, label="equidistant start")
    ax.set_ylim(-0.3, 1.3); ax.set_yticks([])
    ax.legend(loc="upper right", frameon=True)


# =========================================================================== #
# PROBLEM 1 plots : solution vs exact, true error, estimated error, mesh
# =========================================================================== #
def plot_problem1(results, FIG):
    ok = [x for x in results if x[1].converged]
    if not ok:
        print("  Problem 1: no converged result to plot."); return
    e, r, dg, s = next((x for x in ok if x[0]["orders"] == [8]), ok[0])

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.plot(r.t, r.Y, "-", color=BLUE, lw=2.5, label="numerical")
    if dg["exact"] is not None:
        ax.plot(r.t, dg["exact"], "--", color=RED, lw=1.2, label="exact $t^2\\ln t$")
    ax.set_xlabel("$t$"); ax.set_ylabel("$y(t)$"); ax.set_xlim(0, 1); ax.legend()
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_1_solution.png", dpi=250); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for (ei, ri, dgi, _) in ok:
        if dgi["exact"] is not None:
            ax.semilogy(ri.t, np.abs(ri.Y - dgi["exact"]) + 1e-20,
                        lw=1, label=f"order {ei['orders'][0]}")
    ax.set_xlabel("$t$"); ax.set_ylabel("$|Y - y_{\\rm exact}|$"); ax.set_xlim(0, 1); ax.legend()
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_1_error.png", dpi=250); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for (ei, ri, _, _) in ok:
        ax.semilogy(ri.t, ri.err_node + 1e-20, lw=1, label=f"order {ei['orders'][0]}")
    ax.set_xlabel("$t$"); ax.set_ylabel("estimated $|e_i|$"); ax.set_xlim(0, 1); ax.legend()
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_1_est_error.png", dpi=250); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    _mesh_plot(ax, r.t_initial, r.t)
    ax.set_xlabel("$t$"); ax.set_xlim(0, 1)
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_1_mesh.png", dpi=250); plt.close(fig)
    print(f"  Problem 1: saved solution, error, est_error, mesh -> {FIG}/")


# =========================================================================== #
# PROBLEM 2 plots : the three solution branches, mesh per branch
# =========================================================================== #
def plot_problem2(results, FIG):
    sty = {"y1": (BLUE, "."), "y2": (RED, "s"), "y3": (MAG, "o")}
    ok6 = {e["branch"]: (r, dg) for (e, r, dg, s) in results
           if r.converged and e["orders"] == [6] and e["branch"] in sty}
    if not ok6:
        print("  Problem 2: no converged order-6 result."); return

    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for br, (r, dg) in ok6.items():
        c, m = sty[br]
        ax.plot(r.t, r.Y, "-", color=c, marker=m, ms=3.5, mfc=c, lw=1, label=f"branch {br}")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.legend()
    ax.set_xlabel("$t$"); ax.set_ylabel("$y(t)$")
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_2_branches.png", dpi=250); plt.close(fig)

    fig, axs = plt.subplots(3, 1, figsize=(6.4, 9))
    for i, br in enumerate(("y1", "y2", "y3")):
        if br in ok6:
            r, _ = ok6[br]
            _mesh_plot(axs[i], r.t_initial, r.t)
            axs[i].set_xlim(0, 1); axs[i].set_xlabel(f"branch {br}, order 6")
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_2_mesh.png", dpi=250); plt.close(fig)
    print(f"  Problem 2: saved branches, mesh -> {FIG}/")


# =========================================================================== #
# PROBLEM 3 plots : L = 5 / 20, Robin and Dirichlet right BC
# =========================================================================== #
def plot_problem3(results, FIG):
    ok = [x for x in results if x[1].converged and x[0]["orders"] == [6]]
    if not ok:
        print("  Problem 3: no converged order-6 result."); return

    def _L(e):  return float(e["cfg"].split("=")[1].split(",")[0])
    def _bc(e): return "robin" if "y(L)+" in e["cfg"] else "dirichlet"

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
    caps = {"robin": "(a) $y'(0)=0,\\quad y(L)+y'(L)=0$",
            "dirichlet": "(b) $y'(0)=0,\\quad y(L)=0$"}
    for ax, bc in zip(axs, ("robin", "dirichlet")):
        for L, c, m in ((20, BLUE, "."), (5, RED, "+")):
            sub = [(e, r) for (e, r, _, s) in ok if _bc(e) == bc and _L(e) == L]
            if sub:
                ax.plot(sub[0][1].t, sub[0][1].Y, "-", color=c, marker=m, ms=4,
                        lw=1.2, label=f"$L={L}$")
        ax.legend(); ax.set_xlabel(caps[bc], fontsize=11)
        ax.set_xlim(0, 20); ax.set_ylim(0, 2.5)
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_3_solutions.png", dpi=250); plt.close(fig)

    fig, axs = plt.subplots(2, 1, figsize=(6.6, 6.5))
    for i, L in enumerate((5, 20)):
        sub = [(e, r) for (e, r, _, s) in ok if _bc(e) == "robin" and _L(e) == L]
        if sub:
            _mesh_plot(axs[i], sub[0][1].t_initial, sub[0][1].t)
            axs[i].set_xlim(0, L)
            axs[i].set_xlabel(f"$L={L}$, $y(L)+y'(L)=0$, order 6")
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_3_mesh.png", dpi=250); plt.close(fig)
    print(f"  Problem 3: saved solutions, mesh -> {FIG}/")


# =========================================================================== #
# PROBLEM 4 plots : xi = 0.1 and xi = 0.6, comparison, mesh
# =========================================================================== #
def plot_problem4(results, FIG):
    ok = {e["cfg"].split(",")[0].replace("xi=", ""): (r, dg)
          for (e, r, dg, s) in results
          if r.converged and e["orders"] == [8] and "extra" not in e["cfg"]}
    if not ok:
        print("  Problem 4: no converged order-8 result."); return

    for xi, fn in (("0.1", "problem_4_xi_0_1.png"), ("0.6", "problem_4_xi_0_6.png")):
        if xi not in ok: continue
        r, _ = ok[xi]
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
        ax.plot(r.t, r.Y, ".", color="#0000FF", ms=4)
        ax.set_xlim(0, 1); ax.set_xlabel("$\\tau$"); ax.set_ylabel("$y(\\tau)$")
        ax.set_title(f"$\\xi = {xi}$")
        fig.tight_layout(); fig.savefig(f"{FIG}/{fn}", dpi=250); plt.close(fig)

    if len(ok) >= 2:
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        for xi, c in (("0.1", BLUE), ("0.6", RED)):
            if xi in ok:
                ax.plot(ok[xi][0].t, ok[xi][0].Y, "-", color=c, lw=1.6, label=f"$\\xi={xi}$")
        ax.set_xlim(0, 1); ax.set_ylim(-1, 0.6); ax.legend()
        ax.set_xlabel("$\\tau$"); ax.set_ylabel("$y(\\tau)$")
        fig.tight_layout(); fig.savefig(f"{FIG}/problem_4_xi_comparison.png", dpi=250); plt.close(fig)

    xi_list = [xi for xi in ("0.1", "0.6") if xi in ok]
    fig, axs = plt.subplots(len(xi_list), 1, figsize=(6.6, 4 * len(xi_list)), squeeze=False)
    for ax, xi in zip(axs[:, 0], xi_list):
        r, _ = ok[xi]
        _mesh_plot(ax, r.t_initial, r.t)
        ax.set_xlim(0, 1); ax.set_xlabel(f"$\\xi={xi}$, order 8")
    fig.tight_layout(); fig.savefig(f"{FIG}/problem_4_mesh.png", dpi=250); plt.close(fig)
    print(f"  Problem 4: saved xi plots, comparison, mesh -> {FIG}/")


PLOTTERS = {1: plot_problem1, 2: plot_problem2, 3: plot_problem3, 4: plot_problem4}


# ########################################################################### #
#   PART C.  MAIN                                                             #
# ########################################################################### #
def main():
    print("=" * 52)
    print(" HOFiD_bvp — Plotter")
    print("=" * 52)

    problem_filter = ask_problem("plot")
    override_npts = ask_npts(problem_filter)

    exps = get_experiments(problem_filter=problem_filter, override_npts=override_npts)
    run_dir = run_dir_for(problem_filter, override_npts)
    FIG = os.path.join(run_dir, "figures")
    write_inputs_json(run_dir, problem_filter, override_npts, len(exps))
    os.makedirs(FIG, exist_ok=True)

    mesh_msg = (f"{override_npts} points (user override)" if override_npts
                else "defaults (per experiment)")
    print(f"\nInitial mesh : {mesh_msg}")
    print(f"Experiments  : {len(exps)} run(s) selected")
    print(f"Run folder   : {run_dir}/\n")

    by_problem: dict[int, list] = {}
    for idx, e in enumerate(exps):
        print(f"Run {idx+1}/{len(exps)} : {e['prob'].name}"
              + (f" branch {e['branch']}" if e['branch'] else "")
              + f" | orders={e['orders']}", flush=True)
        got = get_result(e, run_dir)
        if got is None:
            continue
        r, dg, summary = got
        by_problem.setdefault(e["problem"], []).append((e, r, dg, summary))

    for pnum, results in sorted(by_problem.items()):
        PLOTTERS[pnum](results, FIG)

    print(f"\nFigures written to {FIG}/:")
    for f in sorted(os.listdir(FIG)):
        print(f"  {f}")


if __name__ == "__main__":
    main()
