"""
run_all_problems.py -- run the four benchmark problems of Hohenegger et al. (2024)
with the Python HOFiD_bvp re-implementation and write tables.

    python run_all_problems.py            # all experiments
    python run_all_problems.py --quick    # skip the most expensive runs (p=4, L=20)

Outputs (results/):  tables/problem_{1..4}.csv, tables/all_results.csv,
                     run_summary.txt, data/*.npz (meshes + solutions for the plots).

PUBLISHED values (columns paper_*) were transcribed by hand from Tables 2, 3, 5, 7 of
the paper and are kept in separate columns; nothing in a computed column is copied
from the paper.
"""
from __future__ import annotations

import argparse
import os
import sys
import numpy as np
import pandas as pd

from hofid_bvp import BVP, hofid_solve, dirichlet_right, robin_right

OUT = "results"
NAN = float("nan")


# --------------------------------------------------------------------------- #
# Benchmark problems (Section 4 of the paper)
# --------------------------------------------------------------------------- #
def problem1() -> BVP:
    """y'' - y'/t - 2 = 0, y(0)=y(1)=0, exact y = t^2 ln t."""
    ex = lambda t: np.where(np.asarray(t) > 0, np.asarray(t) ** 2 * np.log(np.maximum(t, 1e-300)), 0.0)
    return BVP("Problem 1", 1.0,
               f=lambda t, y, yp: yp / t + 2.0,
               fy=lambda t, y, yp: np.zeros_like(y),
               fyp=lambda t, y, yp: 1.0 / t,
               left=("dirichlet", 0.0), right=dirichlet_right(0.0), exact=ex)


def problem2(phi=0.6, gamma=40.0, beta=0.2) -> BVP:
    """y'' + 2/t y' - phi^2 y exp(gamma beta (1-y)/(1+beta(1-y))) = 0, y'(0)=0, y(1)=1."""
    def E(y):
        return gamma * beta * (1 - y) / (1 + beta * (1 - y))
    def f(t, y, yp):
        return -2.0 / t * yp + phi ** 2 * y * np.exp(E(y))
    def fy(t, y, yp):
        dE = -gamma * beta / (1 + beta * (1 - y)) ** 2
        return phi ** 2 * np.exp(E(y)) * (1 + y * dE)
    return BVP("Problem 2", 1.0, f=f, fy=fy, fyp=lambda t, y, yp: -2.0 / t * np.ones_like(y),
               left=("neumann", 0.0), right=dirichlet_right(1.0))


def problem3(L: float, bc: str) -> BVP:
    """y'' + 4/t y' + (t y - 1) y = 0, y'(0)=0; bc = 'robin' (y+y'=0) or 'dirichlet' (y=0) at t=L."""
    right = robin_right() if bc == "robin" else dirichlet_right(0.0)
    return BVP(f"Problem 3 (L={L:g}, {bc})", float(L),
               f=lambda t, y, yp: -4.0 / t * yp - (t * y - 1.0) * y,
               fy=lambda t, y, yp: -(2.0 * t * y - 1.0),
               fyp=lambda t, y, yp: -4.0 / t * np.ones_like(y),
               left=("neumann", 0.0), right=right)


def problem4(xi: float, Nd=3, lam=1.0) -> BVP:
    """Cahn-Hilliard bubble, transformed to tau in [0,1]:
       1/4 (1-tau)^6 y'' + ((N-1)/(2 tau (2-tau)) - 3/4)(1-tau)^5 y' - 4 lam^2 (y+1) y (y-xi) = 0,
       y'(0)=0, y(1)=xi.   Written as w(tau)(y'' - f) with w = (1-tau)^6/4."""
    w = lambda t: (1 - t) ** 6 / 4.0
    c = lambda t: (Nd - 1) / (2 * t * (2 - t)) - 0.75
    def f(t, y, yp):
        return (-c(t) * (1 - t) ** 5 * yp + 4 * lam ** 2 * (y + 1) * y * (y - xi)) / w(t)
    def fy(t, y, yp):
        return 4 * lam ** 2 * (3 * y ** 2 + 2 * (1 - xi) * y - xi) / w(t)
    def fyp(t, y, yp):
        return -c(t) * (1 - t) ** 5 / w(t) * np.ones_like(y)
    return BVP(f"Problem 4 (xi={xi})", 1.0, f=f, fy=fy, fyp=fyp, weight=w,
               left=("neumann", 0.0), right=dirichlet_right(xi))


# initial solution profiles -------------------------------------------------- #
def y0_problem1(t): return np.zeros_like(t)           # not specified in the paper (our choice)
def y0_p2_y1(t): return np.zeros_like(t)              # paper: y == 0
def y0_p2_y2(t): return np.interp(t, [0, .5, 1], [.38, .75, 1.0])   # paper: y(0)=.38,y(.5)=.75,y(1)=1
def y0_p2_y3(t): return np.ones_like(t)               # paper: y == 1
def y0_problem3(t): return np.where(t <= 1.5, 2.0, 2.0 * np.exp(1.5 - t))   # paper
def y0_problem4(xi):
    """NOT given in the paper (refers to report [20]).  Our own tanh bubble profile whose
    centre value / layer position / width were read off Figs. 6-7 of the paper."""
    yc, tm, dl = (-0.3, 0.55, 0.12) if xi < 0.3 else (-1.0, 0.55, 0.03)
    return lambda t: yc + (xi - yc) * 0.5 * (1 + np.tanh((t - tm) / dl))


# --------------------------------------------------------------------------- #
# Published values (HOFiD_bvp tables of the paper): key -> (mesh, abs/est, rel/est[, exact abs, exact rel])
# --------------------------------------------------------------------------- #
PAPER = {
    ("P1", 4): (3585, 3.0059e-09, NAN, 1.4246e-08, 1.4243e-08),
    ("P1", 6): (918, 9.3704e-09, NAN, 7.6807e-08, 7.6806e-08),
    ("P1", 8): (474, 9.8446e-09, NAN, 1.0902e-07, 1.0902e-07),
    ("P1", 10): (571, 9.9055e-09, NAN, 1.3958e-07, 1.3958e-07),
    ("P2y1", 4): (157, 9.8068e-09, 8.4706e-09), ("P2y1", 6): (66, 8.2490e-09, 8.2432e-09),
    ("P2y1", 8): (57, 1.3310e-09, 1.3243e-09), ("P2y1", 10): (54, 1.4665e-10, 1.2778e-10),
    ("P2y2", 4): (167, 6.5398e-09, 4.4201e-09), ("P2y2", 6): (47, 6.5869e-09, 4.8281e-09),
    ("P2y2", 8): (27, 1.1974e-08, 8.6217e-09), ("P2y2", 10): (29, 6.1202e-09, 4.4882e-09),
    ("P2y3", 4): (44, 1.7429e-08, 9.1387e-09), ("P2y3", 6): (18, 6.6477e-09, 3.4857e-09),
    ("P2y3", 8): (15, 6.5300e-10, 3.4240e-10), ("P2y3", 10): (15, 1.3426e-10, 7.0396e-11),
    ("P3robin5", 4): (4033, 2.9461e-09, 9.6905e-10), ("P3robin5", 6): (143, 6.8640e-09, 2.2562e-09),
    ("P3robin5", 8): (169, 1.7346e-09, 5.7010e-10), ("P3robin5", 10): (52, 6.0698e-09, 2.0002e-09),
    ("P3robin5", "var"): (49, 1.2210e-08, 4.0235e-09),
    ("P3robin20", 4): (2077, 1.0856e-08, 6.3889e-09), ("P3robin20", 6): (829, 6.9737e-09, 2.2297e-09),
    ("P3robin20", 8): (188, 8.4732e-09, 8.4655e-09), ("P3robin20", 10): (171, 1.5021e-08, 4.8028e-09),
    ("P3robin20", "var"): (133, 1.5724e-08, 5.0274e-09),
    ("P3dirichlet5", 4): (2753, 3.1070e-08, 9.2501e-09), ("P3dirichlet5", 6): (389, 3.8279e-09, 1.1397e-09),
    ("P3dirichlet5", 8): (78, 6.6038e-09, 4.9040e-09), ("P3dirichlet5", 10): (52, 3.0596e-09, 9.1231e-10),
    ("P3dirichlet5", "var"): (50, 2.6145e-08, 7.7963e-09),
    ("P3dirichlet20", 4): (11105, 2.6538e-09, 8.4660e-10), ("P3dirichlet20", 6): (1781, 4.7281e-09, 1.5121e-09),
    ("P3dirichlet20", 8): (387, 3.8053e-09, 1.2167e-09), ("P3dirichlet20", 10): (261, 3.0595e-08, 9.7821e-09),
    ("P3dirichlet20", "var"): (141, 8.0393e-09, 2.5704e-09),
    ("P4xi0.1", 4): (787, 1.1073e-09, 8.4874e-10), ("P4xi0.1", 6): (247, 4.6386e-09, 3.5608e-09),
    ("P4xi0.1", 8): (86, 8.2062e-09, 6.2899e-09), ("P4xi0.1", 10): (48, 1.2494e-08, 9.5763e-09),
    ("P4xi0.1", "var"): (51, 1.1820e-08, 9.0598e-09),
}
VAR_ORDERS = [6, 8, 10]
VAR_TOLS_TABLE = [1e-2, 1e-5, 1e-8]          # as printed in Tables 5 and 7
VAR_TOLS_FOOTNOTE = [1e-4, 1e-6, 1e-8]       # as stated in footnote 8 (Problem 4) -- contradicts Table 7


def experiments(quick: bool):
    ex = []
    for p in (4, 6, 8, 10):
        ex.append(dict(pid="P1", problem=1, branch="", cfg="", paperkey=("P1", p), prob=problem1(),
                       orders=[p], tols=[1e-8], npts=21, y0=y0_problem1, L=1.0))
    for br, y0, key in (("y1", y0_p2_y1, "P2y1"), ("y2", y0_p2_y2, "P2y2"), ("y3", y0_p2_y3, "P2y3")):
        for p in (4, 6, 8, 10):
            npts = (51 if p == 10 else 21) if br == "y1" else 15
            ex.append(dict(pid=key, problem=2, branch=br, cfg="phi=.6,gamma=40,beta=.2",
                           paperkey=(key, p), prob=problem2(), orders=[p], tols=[1e-8],
                           npts=npts, y0=y0, L=1.0))
    for bc in ("robin", "dirichlet"):
        for L, npts in ((5, 32), (20, 48)):
            runs = [([p], [1e-8], p) for p in (4, 6, 8, 10)] + [(VAR_ORDERS, VAR_TOLS_TABLE, "var")]
            for orders, tols, tag in runs:
                if quick and (L == 20 and tag == 4):
                    continue
                ex.append(dict(pid=f"P3{bc}{L}", problem=3, branch="", cfg=f"L={L}, {'y(L)+y\'(L)=0' if bc=='robin' else 'y(L)=0'}",
                               paperkey=(f"P3{bc}{L}", tag), prob=problem3(L, bc), orders=orders,
                               tols=tols, npts=npts, y0=y0_problem3, L=float(L)))
    for xi, npts in ((0.1, 15), (0.6, 51)):
        runs = [([p], [1e-8], p) for p in (4, 6, 8, 10)] + [(VAR_ORDERS, VAR_TOLS_TABLE, "var"),
                                                              (VAR_ORDERS, VAR_TOLS_FOOTNOTE, "var_footnote")]
        for orders, tols, tag in runs:
            ex.append(dict(pid=f"P4xi{xi}", problem=4, branch="", cfg=f"xi={xi}, N=3, lambda=1",
                           paperkey=(f"P4xi{xi}", tag), prob=problem4(xi), orders=orders, tols=tols,
                           npts=npts, y0=y0_problem4(xi), L=1.0))
    # Extra rows (not in the paper): the paper gives the 15-point start only for p=8.  For the orders where
    # Newton fails from 15 points we additionally try a 51-point start (same profile) and report both.
    for orders, tols, tag in [([4], [1e-8], 4), ([6], [1e-8], 6), (VAR_ORDERS, VAR_TOLS_TABLE, "var"),
                              (VAR_ORDERS, VAR_TOLS_FOOTNOTE, "var_footnote")]:
        ex.append(dict(pid="P4xi0.1", problem=4, branch="", cfg="xi=0.1, N=3, lambda=1 [extra: 51-pt start]",
                       paperkey=("P4xi0.1", tag), prob=problem4(0.1), orders=orders, tols=tols,
                       npts=51, y0=y0_problem4(0.1), L=1.0))
    return ex


def classify_p2(y):
    y0 = y[0]
    if y0 < 0.1: return "y1"
    if 0.2 < y0 < 0.6: return "y2"
    if y0 > 0.8: return "y3"
    return "other"


def run_one(e, idx):
    prob = e["prob"]
    t0 = np.linspace(0.0, e["L"], e["npts"])
    row = dict(problem=e["problem"], branch=e["branch"], config=e["cfg"],
               orders="/".join(map(str, e["orders"])), tols="/".join(f"{x:g}" for x in e["tols"]),
               bc_variant="", initial_mesh_points=e["npts"])
    if e["problem"] == 3:
        row["bc_variant"] = "robin" if "y(L)+" in e["cfg"] else "dirichlet"
    tag = e["paperkey"][1]
    try:
        r = hofid_solve(prob, e["orders"], e["tols"], t0, e["y0"], maxref=25)
        ok = True
    except Exception as exc:                                  # reported, never silently dropped
        row.update(converged=False, status=f"EXCEPTION: {exc!r}")
        return row, None
    row.update(converged=r.converged, status=r.status, final_mesh_points=r.n_points,
               n_intervals=r.n_points - 1, n_blocks=r.n_blocks, refinements=r.refinements,
               newton_iterations=r.newton_iters, est_abs_err_max=r.err_est_max,
               est_rel_err_max=r.err_rel_est_max, est_scaled_err_max=r.scaled_err_max,
               y_at_0=r.Y[0], y_max=float(np.max(r.Y)), runtime_s=round(r.runtime, 2),
               actual_abs_err=NAN, actual_rel_err=NAN)
    if not np.all(np.isfinite(r.err_node)):          # Newton failed on the very first mesh: Y is only the initial guess
        row.update(y_at_0=NAN, y_max=NAN, est_abs_err_max=NAN, est_rel_err_max=NAN, est_scaled_err_max=NAN)
    elif prob.exact is not None and np.all(np.isfinite(r.Y)):
        ye = prob.exact(r.t)
        d = np.abs(r.Y - ye)
        row["actual_abs_err"] = float(d.max())
        row["actual_rel_err"] = float((d / (1.0 + np.abs(ye))).max())   # mixed/relative measure, see README
    if e["problem"] == 2:
        row["solution_branch_found"] = classify_p2(r.Y)
        row["branch_matches_start"] = row["solution_branch_found"] == e["branch"]
    pv = PAPER.get(e["paperkey"])
    cols = ["paper_mesh", "paper_est_abs_err", "paper_rel_err", "paper_exact_abs_err", "paper_exact_rel_err"]
    pv = tuple(pv) + (NAN,) * (5 - len(pv)) if pv else (NAN,) * 5
    row.update(dict(zip(cols, pv)))
    if pv[0] == pv[0] and r.converged:
        row["mesh_ratio_python_over_paper"] = r.n_points / pv[0]
    # store data for the plots
    os.makedirs(f"{OUT}/data", exist_ok=True)
    fname = f"{OUT}/data/run_{idx:03d}.npz"
    np.savez(fname, t=r.t, Y=r.Y, t_initial=r.t_initial, err=r.err_node,
             meta=np.array([e["pid"], str(e["problem"]), e["branch"], e["cfg"], row["orders"], row["tols"],
                            str(r.converged), row["bc_variant"]]))
    row["data_file"] = fname
    return row, r


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="skip the most expensive runs (p=4, L=20)")
    args = ap.parse_args(argv)
    os.makedirs(f"{OUT}/tables", exist_ok=True)
    rows = []
    for idx, e in enumerate(experiments(args.quick)):
        print(f"[{idx:02d}] {e['prob'].name:34s} {e['branch']:3s} orders={e['orders']} tols={e['tols']} ...", end=" ", flush=True)
        row, r = run_one(e, idx)
        row["id"] = idx
        rows.append(row)
        print(f"{'OK ' if row.get('converged') else 'FAIL'} N={row.get('final_mesh_points', '-')} "
              f"est={row.get('est_abs_err_max', NAN):.2e} t={row.get('runtime_s', NAN)}s  {'' if row.get('converged') else row['status']}")
    df = pd.DataFrame(rows)
    lead = ["id", "problem", "branch", "config", "bc_variant", "orders", "tols", "initial_mesh_points",
            "converged", "status", "final_mesh_points", "n_intervals", "n_blocks", "refinements",
            "newton_iterations", "est_abs_err_max", "est_rel_err_max", "est_scaled_err_max",
            "actual_abs_err", "actual_rel_err", "y_at_0", "y_max", "solution_branch_found",
            "branch_matches_start", "runtime_s", "paper_mesh", "paper_est_abs_err", "paper_rel_err",
            "paper_exact_abs_err", "paper_exact_rel_err", "mesh_ratio_python_over_paper", "data_file"]
    df = df.reindex(columns=[c for c in lead if c in df.columns] + [c for c in df.columns if c not in lead])
    for k in (1, 2, 3, 4):
        df[df.problem == k].to_csv(f"{OUT}/tables/problem_{k}.csv", index=False)
    df.to_csv(f"{OUT}/tables/all_results.csv", index=False)
    write_summary(df)
    write_paper_tables(df)
    print(f"\nwrote {OUT}/tables/*.csv and {OUT}/run_summary.txt")
    return df


def write_summary(df):
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
    sh = ["id", "branch", "config", "orders", "tols", "converged", "final_mesh_points", "est_abs_err_max",
          "actual_abs_err", "paper_mesh", "paper_est_abs_err"]
    with open(f"{OUT}/run_summary.txt", "w") as fh:
        fh.write("HOFiD_bvp Python re-implementation -- computed results vs published (Hohenegger et al. 2024)\n")
        fh.write("est_* = ESTIMATED error |Y_p - Y_(p+2)| (not the true error); actual_* = true error vs exact\n"
                 "solution (Problem 1 only).  paper_* = values printed in the paper (separate columns).\n"
                 "converged=False rows are failures / unmet tolerances and are reported as such.\n\n")
        n_ok = int(df.converged.fillna(False).sum())
        fh.write(f"{n_ok}/{len(df)} runs met the requested tolerance.\n\n")
        for k in (1, 2, 3, 4):
            fh.write(f"=== Problem {k} ===\n")
            d = df[df.problem == k]
            fh.write(d[[c for c in sh if c in d.columns]].to_string(index=False, float_format=lambda x: f"{x:.4e}") + "\n\n")
        bad = df[~df.converged.fillna(False).astype(bool)]
        fh.write("=== Failed / tolerance-not-met runs ===\n")
        fh.write(bad[["id", "problem", "branch", "config", "orders", "tols", "status"]].to_string(index=False) + "\n" if len(bad) else "none\n")
        p2 = df[df.problem == 2]
        if len(p2):
            fh.write("\n=== Problem 2: branch reached from each starting profile ===\n")
            fh.write(p2[["branch", "orders", "y_at_0", "solution_branch_found", "branch_matches_start"]].to_string(index=False) + "\n")


# --------------------------------------------------------------------------- #
# Tables laid out exactly like Tables 2, 3, 5, 7 of the paper (computed values only)
# --------------------------------------------------------------------------- #
def _tol(x):
    return f"{x:.0e}".replace("e-0", "e-")


def _row_cells(r, kind):
    o = str(r["orders"]).replace("/", " ")
    tl = " ".join(_tol(float(x)) for x in str(r["tols"]).split("/"))
    if not bool(r["converged"]):
        return [o, tl, "failed", "-", "-", "-"] if kind == "p1" else [o, tl, "failed", "-", "-"]
    e = lambda v: f"{v:.4e}"
    if kind == "p1":
        return [o, tl, f"{int(r['final_mesh_points'])}", e(r["est_abs_err_max"]), e(r["actual_abs_err"]), e(r["actual_rel_err"])]
    return [o, tl, f"{int(r['final_mesh_points'])}", e(r["est_abs_err_max"]), e(r["est_rel_err_max"])]


def _fmt(header, blocks, caption):
    """blocks: list of (subheader or None, list of rows). Monospace table like the paper."""
    allrows = [header] + [r for _, rs in blocks for r in rs]
    w = [max(len(str(x[i])) for x in allrows) + 2 for i in range(len(header))]
    line = lambda cells: "".join(str(c).ljust(w[i]) for i, c in enumerate(cells)).rstrip()
    bar = "-" * sum(w)
    out = [caption, bar, line(header), bar]
    for k, (sub, rs) in enumerate(blocks):
        if sub:
            out += [sub, bar]
        out += [line(r) for r in rs]
        out.append(bar)
    return "\n".join(out) + "\n"


def write_paper_tables(df):
    d = df.copy()
    d["config"] = d["config"].fillna("")
    main = d[~d.config.str.contains("extra")]
    txt = ["Tables in the layout of the paper (Hohenegger et al. 2024), values COMPUTED by this Python code.",
           "In Tables 3, 5, 7 'abs. error' / 'rel. error' are the ESTIMATED errors (as in the paper); in Table 2 'error' is the",
           "estimated absolute error and 'exact abs./rel. error' are true errors. 'mesh' = number of mesh points N+1 in the final mesh.\n"]
    # Table 2
    rs = [_row_cells(r, "p1") for _, r in main[main.problem == 1].iterrows()]
    txt.append(_fmt(["order", "TOL", "mesh", "error", "exact abs. error", "exact rel. error"], [(None, rs)],
                    "Table 2\nProblem 1: Results from HOFiD_bvp (Python re-implementation)."))
    # Table 3
    blocks = []
    for b in ("y1", "y2", "y3"):
        blocks.append((f"{b}(t)   [start mesh {int(main[(main.problem==2)&(main.branch==b)].initial_mesh_points.iloc[0])} pts"
                       f"{' (p=10: 51)' if b=='y1' else ''}]",
                       [_row_cells(r, "p2") for _, r in main[(main.problem == 2) & (main.branch == b)].iterrows()]))
    txt.append(_fmt(["order", "TOL", "mesh", "abs. error", "rel. error"], blocks,
                    "Table 3\nProblem 2: results of the HOFiD_bvp Python re-implementation (starting profiles of the paper)."))
    # Table 5
    hdr = ["L", "order", "TOL", "mesh", "abs. error", "rel. error"]
    blocks = []
    for bc, name in (("robin", "y'(0) = 0,  y(L) + y'(L) = 0"), ("dirichlet", "y'(0) = 0,  y(L) = 0")):
        rows = []
        for L in (5, 20):
            sub = main[(main.problem == 3) & (main.bc_variant == bc) & main.config.str.contains(f"L={L},")]
            for _, r in sub.iterrows():
                c = _row_cells(r, "p2")
                rows.append([str(L)] + c)
            rows.append([""] * len(hdr)) if L == 5 else None
        blocks.append((name, rows))
    txt.append(_fmt(hdr, blocks, "Table 5\nProblem 3: Numerical results computed using the Python HOFiD_bvp with an initial mesh of 32 and 48 points\nfor L = 5 and L = 20, respectively."))
    # Table 7 (xi = 0.1, as in the paper) + xi = 0.6 (no published counterpart)
    for xi, num in ((0.1, "Table 7"), (0.6, "Table 7 (extension, xi = 0.6; not in the paper)")):
        sub = main[(main.problem == 4) & main.config.str.startswith(f"xi={xi},")]
        rows = [_row_cells(r, "p2") for _, r in sub.iterrows()]
        txt.append(_fmt(["order", "TOL", "mesh", "abs. error", "rel. error"], [(None, rows)],
                        f"{num}\nProblem 4, xi = {xi}: Results obtained from the Python HOFiD_bvp "
                        f"(start mesh {int(sub.initial_mesh_points.iloc[0])} pts; last row = tolerances of footnote 8)."))
    ex = d[d.config.str.contains("extra")]
    if len(ex):
        rows = [_row_cells(r, "p2") for _, r in ex.iterrows()]
        txt.append(_fmt(["order", "TOL", "mesh", "abs. error", "rel. error"], [(None, rows)],
                        "Table 7 (additional rows, NOT in the paper): xi = 0.1 with a 51-point start mesh."))
    open(f"{OUT}/tables/paper_format_tables.txt", "w").write("\n".join(txt))
    print(f"wrote {OUT}/tables/paper_format_tables.txt")


if __name__ == "__main__":
    main()
