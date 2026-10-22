"""
run_all_problems.py -- interactive runner for the HOFiD_bvp benchmark problems.

    python run_all_problems.py

WHAT LIVES WHERE
    hofid_bvp.py          the numerical process only (FD operators, Newton,
                          mesh adaptation, error estimate, diagnostics)
    run_all_problems.py   (this file) Problem 1-4 definitions, start profiles,
                          the experiment list, terminal reporting, file output
    plot_all_problems.py  figures; imports everything it needs from this file

EVERY RUN GETS ITS OWN SUB-FOLDER (named after the two answers you give):

    results/problem_1_npts_15/          you chose Problem 1 and 15 points
    results/problem_1_npts_default/     you chose Problem 1 and pressed Enter
    results/problem_all_npts_21/        you chose "all" and 21 points
        inputs.json                     the answers that created this folder
        tables/        <run>_initial.csv, <run>_final.csv
        data/          <run>.npz, <run>_summary.json
        diagnostics/   <run>_coefficients.csv
        figures/       (written by plot_all_problems.py)

Running the same inputs again re-uses (overwrites) the same folder.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

from hofid_bvp import (
    BVP, hofid_solve, dirichlet_right, robin_right,
    compute_diagnostics, summarize_diagnostics, coefficient_table, stencil_labels,
)

ROOT = "results"


# ########################################################################### #
#                                                                             #
#   PART A.  BENCHMARK PROBLEM DEFINITIONS                                    #
#   (each problem = f, df/dy, df/dy', boundary conditions, optional weight)   #
#                                                                             #
# ########################################################################### #

# =========================================================================== #
# PROBLEM 1  --  singular linear BVP with a known exact solution
#
#     y'' - y'/t - 2 = 0 ,   y(0) = y(1) = 0 ,   exact solution  y = t^2 ln t
#
# Singular at t = 0 (the term y'/t).  The only problem with an exact solution,
# so it is the only one where "actual error" columns are filled in.
# =========================================================================== #
def problem1() -> BVP:
    ex = lambda t: np.where(np.asarray(t) > 0,
                            np.asarray(t)**2 * np.log(np.maximum(t, 1e-300)), 0.0)
    return BVP("Problem 1", 1.0,
               f=lambda t, y, yp: yp / t + 2.0,
               fy=lambda t, y, yp: np.zeros_like(y),
               fyp=lambda t, y, yp: 1.0 / t,
               left=("dirichlet", 0.0), right=dirichlet_right(0.0), exact=ex)


# =========================================================================== #
# PROBLEM 2  --  nonlinear spherical reaction-diffusion (three solution branches)
#
#     y'' + (2/t) y' - phi^2 y exp( gamma*beta*(1-y) / (1+beta*(1-y)) ) = 0
#     y'(0) = 0 ,   y(1) = 1 ,   phi = 0.6, gamma = 40, beta = 0.2
#
# Three different solutions exist; which one Newton finds depends on the start
# profile (y1, y2, y3 below).
# =========================================================================== #
def problem2(phi=0.6, gamma=40.0, beta=0.2) -> BVP:
    def E(y):
        return gamma * beta * (1 - y) / (1 + beta * (1 - y))

    def f(t, y, yp):
        return -2.0 / t * yp + phi**2 * y * np.exp(E(y))

    def fy(t, y, yp):
        dE = -gamma * beta / (1 + beta * (1 - y))**2
        return phi**2 * np.exp(E(y)) * (1 + y * dE)

    return BVP("Problem 2", 1.0, f=f, fy=fy,
               fyp=lambda t, y, yp: -2.0 / t * np.ones_like(y),
               left=("neumann", 0.0), right=dirichlet_right(1.0))


# =========================================================================== #
# PROBLEM 3  --  singular nonlinear BVP on [0, L]  (L = 5 or 20)
#
#     y'' + (4/t) y' + (t y - 1) y = 0 ,   y'(0) = 0
#     right BC:  Robin     y(L) + y'(L) = 0      or     Dirichlet  y(L) = 0
# =========================================================================== #
def problem3(L: float, bc: str) -> BVP:
    right = robin_right() if bc == "robin" else dirichlet_right(0.0)
    return BVP(f"Problem 3 (L={L:g}, {bc})", float(L),
               f=lambda t, y, yp: -4.0 / t * yp - (t * y - 1.0) * y,
               fy=lambda t, y, yp: -(2.0 * t * y - 1.0),
               fyp=lambda t, y, yp: -4.0 / t * np.ones_like(y),
               left=("neumann", 0.0), right=right)


# =========================================================================== #
# PROBLEM 4  --  Cahn-Hilliard "bubble", transformed to tau in [0, 1]
#
#     w(t) (y'' - f) = 0 ,   w(t) = (1-t)^6 / 4   (weight keeps the equation
#     regular at t = 1),     y'(0) = 0 ,   y(1) = xi   (xi = 0.1 or 0.6)
# =========================================================================== #
def problem4(xi: float, Nd=3, lam=1.0) -> BVP:
    w = lambda t: (1 - t)**6 / 4.0
    c = lambda t: (Nd - 1) / (2 * t * (2 - t)) - 0.75

    def f(t, y, yp):
        return (-c(t) * (1 - t)**5 * yp + 4 * lam**2 * (y + 1) * y * (y - xi)) / w(t)

    def fy(t, y, yp):
        return 4 * lam**2 * (3 * y**2 + 2 * (1 - xi) * y - xi) / w(t)

    def fyp(t, y, yp):
        return -c(t) * (1 - t)**5 / w(t) * np.ones_like(y)

    return BVP(f"Problem 4 (xi={xi})", 1.0, f=f, fy=fy, fyp=fyp, weight=w,
               left=("neumann", 0.0), right=dirichlet_right(xi))


# =========================================================================== #
# NEWTON START PROFILES  (initial guesses Y0(t) handed to the solver)
# =========================================================================== #
def y0_problem1(t): return np.zeros_like(t)
def y0_p2_y1(t):    return np.zeros_like(t)
def y0_p2_y2(t):    return np.interp(t, [0, .5, 1], [.38, .75, 1.0])
def y0_p2_y3(t):    return np.ones_like(t)
def y0_problem3(t): return np.where(t <= 1.5, 2.0, 2.0 * np.exp(1.5 - t))


def y0_problem4(xi):
    """Tanh bubble profile (centre value yc, jump near tm, width dl)."""
    yc, tm, dl = (-0.3, 0.55, 0.12) if xi < 0.3 else (-1.0, 0.55, 0.03)
    return lambda t: yc + (xi - yc) * 0.5 * (1 + np.tanh((t - tm) / dl))


# =========================================================================== #
# EXPERIMENT LIST  (which problem / order / tolerance / start mesh is run)
# =========================================================================== #
VAR_ORDERS        = [6, 8, 10]
VAR_TOLS_TABLE    = [1e-2, 1e-5, 1e-8]
VAR_TOLS_FOOTNOTE = [1e-4, 1e-6, 1e-8]


def get_experiments(problem_filter=None, override_npts=None):
    """List of experiment dicts for the selected problem(s).

    problem_filter : None = all, or a set/list of ints {1,2,3,4}
    override_npts  : None = default start mesh per experiment, int = use for all
    """
    pf = None if problem_filter is None else set(problem_filter)
    ex = []

    def add(default, src, fallback, **kw):
        if pf and kw["problem"] not in pf:
            return
        npts = default if override_npts is None else override_npts
        # unique, folder-safe name of this experiment (used for all its files)
        parts = [f"problem_{kw['problem']}"]
        if kw["branch"]:
            parts.append(kw["branch"])
        if kw["problem"] == 3:
            parts.append(kw["pid"][2:])                        # e.g. robin5
        if kw["problem"] == 4:
            parts.append(kw["pid"][2:])                        # e.g. xi0.1
        parts.append("orders" + "_".join(map(str, kw["orders"])))
        if len(kw["orders"]) > 1:
            parts.append("tol" + "_".join(f"{x:g}" for x in kw["tols"]))
        if "extra" in kw["cfg"]:
            parts.append("start51")
        ex.append(dict(kw, slug="_".join(parts),
                       default_npts=default, default_source=src,
                       default_is_fallback=fallback,
                       npts=npts, override=(override_npts is not None)))

    # ---- Problem 1 : orders 4, 6, 8, 10 --------------------------------- #
    for p in (4, 6, 8, 10):
        add(21, "default start mesh", False,
            pid="P1", problem=1, branch="", cfg="",
            prob=problem1(), orders=[p], tols=[1e-8], y0=y0_problem1, L=1.0)

    # ---- Problem 2 : three branches, orders 4..10 ----------------------- #
    for br, y0, key in (("y1", y0_p2_y1, "P2y1"),
                        ("y2", y0_p2_y2, "P2y2"),
                        ("y3", y0_p2_y3, "P2y3")):
        for p in (4, 6, 8, 10):
            npts = (51 if p == 10 else 21) if br == "y1" else 15
            add(npts, "default start mesh", False,
                pid=key, problem=2, branch=br, cfg="phi=.6,gamma=40,beta=.2",
                prob=problem2(), orders=[p], tols=[1e-8], y0=y0, L=1.0)

    # ---- Problem 3 : 2 BCs x 2 L values x fixed + variable order -------- #
    for bc in ("robin", "dirichlet"):
        for L, npts in ((5, 32), (20, 48)):
            runs = [([p], [1e-8]) for p in (4, 6, 8, 10)] + [(VAR_ORDERS, VAR_TOLS_TABLE)]
            for orders, tols in runs:
                bc_label = "y(L)+y'(L)=0" if bc == "robin" else "y(L)=0"
                add(npts, "default start mesh", False,
                    pid=f"P3{bc}{L}", problem=3, branch="",
                    cfg=f"L={L}, {bc_label}",
                    prob=problem3(L, bc), orders=orders, tols=tols,
                    y0=y0_problem3, L=float(L))

    # ---- Problem 4 : two xi values -------------------------------------- #
    for xi, npts in ((0.1, 15), (0.6, 51)):
        runs = [([p], [1e-8]) for p in (4, 6, 8, 10)] + \
               [(VAR_ORDERS, VAR_TOLS_TABLE), (VAR_ORDERS, VAR_TOLS_FOOTNOTE)]
        for orders, tols in runs:
            fb = not (xi == 0.1 and orders == [8])
            add(npts, "default start mesh" if not fb else "fallback start mesh", fb,
                pid=f"P4xi{xi}", problem=4, branch="",
                cfg=f"xi={xi}, N=3, lambda=1",
                prob=problem4(xi), orders=orders, tols=tols,
                y0=y0_problem4(xi), L=1.0)
        # extra rows with a 51-point start (only when using default meshes)
        if override_npts is None and xi == 0.1:
            for orders, tols in [([4], [1e-8]), ([6], [1e-8]),
                                 (VAR_ORDERS, VAR_TOLS_TABLE),
                                 (VAR_ORDERS, VAR_TOLS_FOOTNOTE)]:
                add(51, "extra row, 51-point start", True,
                    pid="P4xi0.1", problem=4, branch="",
                    cfg="xi=0.1, N=3, lambda=1 [extra: 51-pt start]",
                    prob=problem4(0.1), orders=orders, tols=tols,
                    y0=y0_problem4(0.1), L=1.0)
    return ex


def validate_npts(npts, min_p):
    """Return an error string if npts is invalid, else None."""
    min_req = min_p + 4
    if npts < min_req:
        return (f"Too small: order p={min_p} needs at least {min_req} initial mesh points "
                f"(blocks need >= p+4 points; the p+2 estimator needs N_intervals >= p+3). "
                f"Minimum for the selected problem: {min_req}.")
    if npts > 50000:
        return f"Too large: limit is 50 000 initial points (got {npts})."
    return None


# ########################################################################### #
#                                                                             #
#   PART B.  RUN FOLDERS  (one sub-folder per distinct set of inputs)         #
#                                                                             #
# ########################################################################### #
def run_dir_for(problem_filter, override_npts, root=ROOT):
    """results/problem_<n|all>_npts_<N|default>  -- same inputs => same folder."""
    prob = "all" if problem_filter is None else "_".join(map(str, sorted(problem_filter)))
    npts = "default" if override_npts is None else str(int(override_npts))
    return os.path.join(root, f"problem_{prob}_npts_{npts}")


def write_inputs_json(run_dir, problem_filter, override_npts, n_experiments):
    os.makedirs(run_dir, exist_ok=True)
    with open(os.path.join(run_dir, "inputs.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(problem=("all" if problem_filter is None else sorted(problem_filter)),
                       initial_mesh_points=("default" if override_npts is None
                                            else int(override_npts)),
                       experiments=n_experiments), fh, indent=2)


# ########################################################################### #
#                                                                             #
#   PART C.  TERMINAL REPORT + FILE OUTPUT                                    #
#                                                                             #
# ########################################################################### #
def _g(v, fmt=".4e"):
    if v is None or (isinstance(v, float) and v != v):
        return "n/a"
    return format(float(v), fmt)


def _rep_indices(N, p, k=20):
    base = set(range(0, min(5, N+1))) | set(range(max(0, N-4), N+1))
    base |= {p//2, p//2+1, N//4, N//2, 3*N//4, N-p//2-1, N-p//2}
    return sorted(i for i in base if 0 <= i <= N)[:k]


def _pointwise_df(dg):
    """Pointwise solution table as a DataFrame."""
    N, t = dg["N"], dg["t"]
    lab = [stencil_labels(i, N, dg["p"], dg["shifts"][i]) for i in range(N+1)]
    rt = np.array(["interior"]*(N+1), dtype=object)
    rt[0] = "left_BC"; rt[-1] = "right_BC"
    hl = np.r_[np.nan, np.diff(t)]
    hr = np.r_[np.diff(t), np.nan]
    return pd.DataFrame(dict(
        i=np.arange(N+1), t=t, h_left=hl, h_right=hr, row_type=rt,
        Y=dg["Y"], Y_hi=dg["Y_hi"] if dg["Y_hi"] is not None else np.nan,
        y_prime=dg["yp"], y_second=dg["ypp"],
        f_rhs=dg["f"], residual_ode=dg["residual"],
        residual_weighted=dg["residual_weighted"],
        residual_bc=dg["bc_residual"],
        est_abs_err=dg["est_abs_err"], scaled_eta=dg["eta"],
        mixed_scaled=dg["mixed_scaled"],
        y_exact=dg["exact"] if dg["exact"] is not None else np.nan,
        actual_abs_err=dg["actual_abs"] if dg["actual_abs"] is not None else np.nan,
        actual_rel_err=dg["actual_rel"] if dg["actual_rel"] is not None else np.nan,
        rel_valid=dg["rel_valid"] if dg["rel_valid"] is not None else False,
        upwind_shift=dg["shifts"],
        d1_stencil=[l[0] for l in lab],
        d2_stencil=[l[2] for l in lab]))


def _print_section(title):
    bar = "=" * 52
    print(f"\n{bar}\n{title}\n{bar}")


def _print_error_summary(sm, has_exact):
    print(f"\nEst. absolute error      : {_g(sm['max_est_abs_err'])}  (|Y_p+2 - Y_p|)")
    print(f"Scaled error eta         : {_g(sm['max_scaled_eta'])}  (|Y_p+2 - Y_p|/(1+|Y|))")
    if has_exact:
        print(f"Actual absolute error     : {_g(sm['max_actual_abs_err'])}")
        print(f"Actual relative error     : {_g(sm['max_actual_rel_err'])}  "
              f"(valid pts: |y_exact| >= 1e-3*max|y_exact|)")
    else:
        print("Actual absolute error     : n/a (no exact solution for this problem)")
        print("Actual relative error     : n/a")
    print(f"Max interior residual    : {_g(sm['max_interior_residual'])}  (y'' - f)")
    print(f"Max boundary residual    : {_g(sm['max_boundary_residual'])}")


def _print_pointwise_table(dg, indices):
    hdr = ["i", "t_i", "Y_i", "y'_i", "y''_i", "f_i",
           "residual", "est.abs.err", "eta", "actual_abs", "actual_rel"]
    print("\n  " + "  ".join(f"{h:>12}" for h in hdr))
    print("  " + "-"*144)
    for i in indices:
        aa = dg["actual_abs"][i] if dg["actual_abs"] is not None else float("nan")
        ar = float("nan")
        if dg["actual_rel"] is not None and dg["rel_valid"] is not None:
            ar = dg["actual_rel"][i] if dg["rel_valid"][i] else float("nan")
        vals = [dg["t"][i], dg["Y"][i], dg["yp"][i], dg["ypp"][i],
                dg["f"][i], dg["residual"][i],
                dg["est_abs_err"][i], dg["eta"][i], aa, ar]
        print(f"  {i:>12d}  " + "  ".join(f"{_g(float(v), '.4e'):>12}" for v in vals))
    if indices:
        print("  (residual: N/A (singular endpoint) at i=0 and i=N; "
              "actual_rel: N/A where |y_exact|~0 or unavailable)")


def _save_summary_json(path, meta):
    os.makedirs(os.path.dirname(path), exist_ok=True)

    def _serial(v):
        if isinstance(v, np.integer):  return int(v)
        if isinstance(v, np.floating): return float(v)
        if isinstance(v, np.ndarray):  return v.tolist()
        return v
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({k: _serial(v) for k, v in meta.items()}, fh, indent=2)


def _classify_p2(y):
    y0 = y[0]
    if y0 < 0.1:       return "y1"
    if 0.2 < y0 < 0.6: return "y2"
    if y0 > 0.8:       return "y3"
    return "other"


def solve_and_report(e, out_dir=ROOT, show_table=True, max_table_rows=20,
                     print_progress=True):
    """Run ONE experiment e (from get_experiments), print the terminal report,
    save CSV/JSON/NPZ under out_dir, return (result, diagnostics, summary).

    out_dir is the per-run folder (see run_dir_for).  All numerics are done by
    hofid_solve / compute_diagnostics in hofid_bvp.py.
    """
    prob = e["prob"]
    t0 = np.linspace(0.0, e["L"], e["npts"])
    orders_str = "/".join(map(str, e["orders"]))
    tols_str = "/".join(f"{x:g}" for x in e["tols"])
    has_exact = prob.exact is not None
    pnum, branch, slug = e["problem"], e["branch"], e["slug"]
    os.makedirs(f"{out_dir}/tables", exist_ok=True)
    os.makedirs(f"{out_dir}/data", exist_ok=True)
    os.makedirs(f"{out_dir}/diagnostics", exist_ok=True)

    # ---------------- INITIAL MESH (one solve, no refinement) ------------- #
    _print_section(f"PROBLEM {pnum}{' | branch ' + branch if branch else ''} — INITIAL MESH RESULTS")
    print(f"  Configuration  : {e['cfg'] or '(see problem definition)'}")
    print(f"  Orders/tols    : {orders_str} / {tols_str}")
    print(f"  Initial mesh   : {e['npts']} points  "
          f"({'user override' if e['override'] else 'default'}"
          f"{' (fallback)' if e['default_is_fallback'] and not e['override'] else ''})")
    print(f"  Initial domain : [0, {e['L']:g}]  ({e['npts']-1} intervals)")

    p0 = e["orders"][0]
    try:
        r0 = hofid_solve(prob, [p0], [e["tols"][0]], t0, e["y0"], maxref=0)
    except Exception as exc:
        print(f"\n  SOLVER ERROR on initial mesh: {exc!r}")
        r0 = None

    csv0 = None
    if r0 is not None and r0.Y_hi is not None:
        dg0 = compute_diagnostics(prob, r0)
        sm0 = summarize_diagnostics(dg0, r0)
        print(f"\n  Finite diff. order : {r0.p_final}  (estimator uses {r0.p_final+2})")
        print(f"  Newton iterations  : {r0.newton_iters}")
        print(f"  Newton converged   : {'YES' if r0.newton_converged else 'NO'}")
        _print_error_summary(sm0, has_exact)
        print(f"\n  Initial numerical solution (showing up to {max_table_rows} points):")
        print(f"  (full table saved to CSV)")
        if show_table:
            _print_pointwise_table(dg0, _rep_indices(dg0["N"], dg0["p"], max_table_rows))
        csv0 = f"{out_dir}/tables/{slug}_initial.csv"
        _pointwise_df(dg0).to_csv(csv0, index=False)
        print(f"\n  Saved: {csv0}")
    else:
        sm0 = None
        print("\n  Initial solve did not converge (no error estimates available).")
        if r0:
            print(f"  Status: {r0.status}")

    # ---------------- ADAPTIVE REFINEMENT --------------------------------- #
    print(f"\nStarting adaptive mesh refinement...")
    try:
        r = hofid_solve(prob, e["orders"], e["tols"], t0, e["y0"], maxref=25,
                        verbose=print_progress)
    except Exception as exc:
        print(f"\n  SOLVER ERROR: {exc!r}")
        summary = dict(converged=False, status=f"exception: {exc!r}",
                       termination="solver exception", problem=pnum, branch=branch,
                       orders=orders_str, tols=tols_str, initial_npts=e["npts"],
                       final_npts=None, est_abs_err_max=float("nan"))
        _save_summary_json(f"{out_dir}/data/{slug}_summary.json", summary)
        return None, None, summary

    seen = set()
    for h in r.history:
        k = (h["stage"], h["iteration"], h["points"])
        if k not in seen:
            seen.add(k)
            if h["iteration"] > 0 or h["stage"] > 0:
                print(f"  Refinement {h['iteration']+h['stage']}:  "
                      f"mesh points = {h['points']:6d}   "
                      f"est. error = {h['est_abs_err']:.3e}   "
                      f"tol = {h['tol']:.1e}")

    # ---------------- FINAL RESULTS --------------------------------------- #
    _print_section(f"PROBLEM {pnum}{' | branch ' + branch if branch else ''} — FINAL ADAPTED MESH RESULTS")
    print(f"  Initial mesh points  : {e['npts']}")

    if r.Y_hi is None:
        print(f"\n  NO CONVERGED SOLUTION.")
        print(f"  Status    : {r.status}")
        print(f"  Termination: {r.termination}")
        print(f"  Newton converged: {'YES' if r.newton_converged else 'NO'}")
        summary = dict(converged=False, status=r.status, termination=r.termination,
                       problem=pnum, branch=branch, orders=orders_str, tols=tols_str,
                       initial_npts=e["npts"], final_npts=None,
                       est_abs_err_max=float("nan"))
        _save_summary_json(f"{out_dir}/data/{slug}_summary.json", summary)
        return r, None, summary

    dg = compute_diagnostics(prob, r)
    sm = summarize_diagnostics(dg, r)

    print(f"  Final mesh points    : {r.n_points}  ({r.n_blocks} equidistant blocks)")
    print(f"  Number of refinements: {r.refinements}")
    print(f"  Finite diff. order   : {r.p_final}  (estimator uses {r.p_final+2})")
    print(f"  Tolerance (last stage): {r.tol_final:g}")
    print(f"  Newton converged     : {'YES' if r.newton_converged else 'NO'}")
    print(f"  Tolerance satisfied  : {'YES' if r.converged else 'NO'}  "
          f"(scaled criterion: {r.scaled_err_max:.3e} <= 1 required)")
    print(f"  Termination reason   : {r.termination}")
    if not r.converged:
        print(f"  Status               : {r.status}")

    if sm0 is not None:
        print(f"\n  Initial est. abs. error : {_g(sm0['max_est_abs_err'])}")
    print(f"  Final   est. abs. error : {_g(sm['max_est_abs_err'])}")
    _print_error_summary(sm, has_exact)

    if pnum == 2:
        found = _classify_p2(r.Y)
        print(f"\n  Solution branch found  : {found}  "
              f"(started from profile '{branch}' -> "
              f"{'match' if found == branch else 'DIFFERENT BRANCH'})")

    print(f"\n  Final numerical solution (showing up to {max_table_rows} points):")
    print(f"  (full table saved to CSV)")
    if show_table:
        _print_pointwise_table(dg, _rep_indices(dg["N"], dg["p"], max_table_rows))

    # ---------------- SAVE ------------------------------------------------ #
    csv_final = f"{out_dir}/tables/{slug}_final.csv"
    _pointwise_df(dg).to_csv(csv_final, index=False)
    coeff_csv = f"{out_dir}/diagnostics/{slug}_coefficients.csv"
    pd.DataFrame(coefficient_table(dg, r, prob)).to_csv(coeff_csv, index=False)
    npz_path = f"{out_dir}/data/{slug}.npz"
    np.savez(npz_path, t=r.t, Y=r.Y, t_initial=r.t_initial, err=r.err_node,
             t_initial_mesh=t0, Y_initial=e["y0"](t0))

    summary = dict(
        problem=pnum, branch=branch, config=e["cfg"], slug=slug,
        orders=orders_str, tols=tols_str,
        initial_npts=e["npts"], default_npts=e["default_npts"],
        default_source=e["default_source"], override=e["override"],
        converged=bool(r.converged), status=r.status,
        termination=r.termination, newton_converged=bool(r.newton_converged),
        tolerance_satisfied=bool(r.converged), tolerance_requested=r.tol_final,
        final_npts=r.n_points, n_blocks=r.n_blocks, refinements=r.refinements,
        newton_iters=r.newton_iters, fd_order=r.p_final, runtime_s=round(r.runtime, 2),
        est_abs_err_max=sm["max_est_abs_err"],
        est_rel_err_max=r.err_rel_est_max,
        scaled_err_max=r.scaled_err_max,
        max_scaled_eta=sm["max_scaled_eta"],
        max_interior_residual=sm["max_interior_residual"],
        max_boundary_residual=sm["max_boundary_residual"],
        actual_abs_err=sm["max_actual_abs_err"],
        actual_rel_err=sm["max_actual_rel_err"],
        csv_initial=csv0, csv_final=csv_final, coeff_csv=coeff_csv, npz=npz_path)
    _save_summary_json(f"{out_dir}/data/{slug}_summary.json", summary)

    print(f"\n  Saved: {csv_final}")
    print(f"         {coeff_csv}")
    print(f"         {out_dir}/data/{slug}_summary.json")
    return r, dg, summary


# ########################################################################### #
#                                                                             #
#   PART D.  INTERACTIVE INPUT  (shared with plot_all_problems.py)            #
#                                                                             #
# ########################################################################### #
def ask_problem(verb="solve"):
    print(f"\nWhich problem do you want to {verb}?\n")
    print("  1. Problem 1")
    print("  2. Problem 2")
    print("  3. Problem 3")
    print("  4. Problem 4")
    print("  5. All problems\n")
    while True:
        try:
            txt = input("Enter your choice (1/2/3/4/5): ").strip().lower()
        except EOFError:
            print("\n(no input -> all problems)")
            return None
        if txt in ("5", "all", ""):
            return None          # None means all
        if txt in ("1", "2", "3", "4"):
            return {int(txt)}
        print(f"  Invalid choice '{txt}'. Please enter 1, 2, 3, 4, or 5 (or 'all').")


def ask_npts(problem_filter):
    """Ask for the initial mesh size. Returns int, or None (= defaults)."""
    exps = get_experiments(problem_filter=problem_filter)
    if not exps:
        print("  No experiments found for that selection.")
        sys.exit(1)
    # The estimator uses p+2 and build_operators needs N_intervals >= (p+2)+2,
    # so points >= p+5 for the highest selected order.
    max_p = max(max(e["orders"]) for e in exps)
    min_points = max_p + 5

    print("\nEnter the number of initial mesh points.")
    print("Press Enter to use the default start mesh of each experiment:\n")
    print(f"  (Points include both end points; intervals = points - 1.")
    print(f"   Minimum for the selected problem(s): {min_points} points.")
    print(f"   Enter a whole number, or press Enter for the defaults.)\n")

    while True:
        try:
            txt = input("Initial mesh points: ").strip()
        except EOFError:
            print("\n(no input -> using defaults)")
            return None
        if txt == "" or txt.lower() == "default":
            return None
        try:
            n = int(txt)
        except ValueError:
            print(f"  Invalid input '{txt}': please enter a whole number.")
            continue
        err = validate_npts(n, max_p)
        if not err and n < min_points:
            err = (f"Need at least {min_points} points for the highest selected "
                   f"order ({max_p}) and its error estimator.")
        if err:
            print(f"  {err}")
            continue
        return n


# ########################################################################### #
#                                                                             #
#   PART E.  MAIN                                                             #
#                                                                             #
# ########################################################################### #
def main():
    print("=" * 52)
    print(" HOFiD_bvp —  Runner")
    print("=" * 52)

    problem_filter = ask_problem("solve")
    override_npts = ask_npts(problem_filter)

    exps = get_experiments(problem_filter=problem_filter, override_npts=override_npts)
    run_dir = run_dir_for(problem_filter, override_npts)
    write_inputs_json(run_dir, problem_filter, override_npts, len(exps))

    mesh_msg = (f"{override_npts} points (user override)" if override_npts
                else "defaults (per experiment)")
    print(f"\nInitial mesh : {mesh_msg}")
    print(f"Experiments  : {len(exps)} run(s) selected")
    print(f"Output folder: {run_dir}/\n")

    all_summaries = []
    for idx, e in enumerate(exps):
        print(f"\n{'─'*52}")
        print(f"Run {idx+1}/{len(exps)} : {e['prob'].name}"
              + (f" branch {e['branch']}" if e['branch'] else "")
              + f" | orders={e['orders']} tols={e['tols']}")
        _, _, summary = solve_and_report(e, out_dir=run_dir, show_table=True,
                                         max_table_rows=20, print_progress=False)
        all_summaries.append(summary)

    print("\n" + "=" * 52)
    print(" SUMMARY OF ALL RUNS")
    print("=" * 52)
    ok = sum(1 for s in all_summaries if s and s.get("converged"))
    print(f"  {ok}/{len(all_summaries)} runs met the requested tolerance.\n")
    for s in all_summaries:
        if not s:
            continue
        status = "OK  " if s.get("converged") else "FAIL"
        br = f" [{s.get('branch')}]" if s.get("branch") else ""
        try:
            est_err = f"{float(s.get('est_abs_err_max', float('nan'))):.2e}"
        except (TypeError, ValueError):
            est_err = "N/A"
        print(f"  [{status}] P{s.get('problem', '?')}{br}  orders={s.get('orders', '?')}  "
              f"N_final={s.get('final_npts', 'N/A')}  est_err={est_err}  "
              f"-> {s.get('termination', s.get('status', 'failed'))}")

    print(f"\nDone. Output files are in {run_dir}/")


if __name__ == "__main__":
    main()
