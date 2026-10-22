"""
hofid_bvp.py -- COMPUTE ENGINE of the HOFiD_bvp solver (high-order finite differences).

THIS FILE CONTAINS ONLY THE NUMERICAL PROCESS.  There are no benchmark problems,
no experiment lists, no printing and no file output here:

    * the problems (Problem 1-4) live in      run_all_problems.py
    * printing / CSV / JSON / NPZ output in   run_all_problems.py
    * figures live in                         plot_all_problems.py

Scalar second-order BVP on [0, L] (singularity allowed at t = 0)

    w(t) * ( y'' - f(t, y, y') ) = 0 ,   t in (0, L)
    left  BC : y(0) = v    ("dirichlet")   or   y'(0) = v   ("neumann")
    right BC : b(y(L), y'(L)) = 0

Method
------
* y'' and y' are discretised separately (no first-order system).
* y'' : central formula, p+1 points in the interior; one-sided p+2 point
  formulas near the two ends (order p preserved).
* y'  : p+1 point formula, upwind-biased (window shifted one node according to the
  sign of df/dy'); one-sided formulas at the ends.
* The ODE is never evaluated at t = 0 nor at t = L: the BC rows replace it.
* Damped Newton iteration (affine-invariant monotonicity test), sparse LU.
* Mesh is piecewise equidistant: blocks with >= p+4 points.
* Error estimate: |Y_p - Y_{p+2}| on the same mesh.
* Error equidistribution + grading restriction -> new blocks.
* Variable order: (orders, tols) are processed as successive stages.
"""


from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.interpolate import CubicSpline

__all__ = [
    "BVP", "HOFiDResult", "fd_weights", "build_operators", "newton_solve",
    "hofid_solve", "hofid_bvp", "dirichlet_right", "robin_right",
    "blocks_from_mesh", "build_mesh", "compute_diagnostics",
    "summarize_diagnostics", "coefficient_table", "stencil_labels", "run_tests",
]



###############################################################################
#                                                                             #
#   1. VANDERMONDE-TYPE MATRIX FORMATION  (FINITE-DIFFERENCE WEIGHTS)
#   
#   For a stencil x_0..x_p and an evaluation point z, the weights w_j of the
#   m-th derivative formula   f^(m)(z) ~ sum_j w_j f(x_j)   solve the
#   Vandermonde system
#   
#           sum_j w_j (x_j - z)^k  =  m! * delta_{k,m}      k = 0..p
#   
#   i.e.  V w = m! e_m  with  V_kj = (x_j - z)^k.  Solving V directly is
#   badly conditioned for large p, so the weights are generated with
#   Fornberg's (1988) recurrence, which produces the SAME numbers stably.
#   On a constant block, offsets are (x_j - z)/h, so the weights are the
#   tabulated integer weights scaled by 1/h^m.  (run_tests() checks
#   fd_weights against a direct Vandermonde solve.)
#                                                                             #
###############################################################################


def fd_weights(z: float, x: np.ndarray, m: int) -> np.ndarray:
    """Weights w_j such that f^(m)(z) ~ sum_j w_j f(x_j) (exact for polynomials
    of degree <= len(x)-1)."""
    n = len(x)
    c = np.zeros((n, m + 1))
    c1, c4 = 1.0, x[0] - z
    c[0, 0] = 1.0
    for i in range(1, n):
        mn = min(i, m)
        c2, c5, c4 = 1.0, c4, x[i] - z
        for j in range(i):
            c3 = x[i] - x[j]
            c2 *= c3
            if j == i - 1:
                for k in range(mn, 0, -1):
                    c[i, k] = c1 * (k * c[i-1, k-1] - c5 * c[i-1, k]) / c2
                c[i, 0] = -c1 * c5 * c[i-1, 0] / c2
            for k in range(mn, 0, -1):
                c[j, k] = (c4 * c[j, k] - k * c[j, k-1]) / c3
            c[j, 0] = c4 * c[j, 0] / c3
        c1 = c2
    return c[:, m]


def build_operators(t: np.ndarray, p: int, shifts: np.ndarray):
    """Sparse D1, D2 approximating y' and y'' at all nodes of t for order p."""
    N = len(t) - 1
    if N < p + 2:
        raise ValueError(f"mesh too coarse for order {p}: need N >= {p+2}, got {N}")
    cache: dict = {}

    def weights(i, s, width, m):
        idx = np.arange(s, s + width)
        h = t[i+1] - t[i] if i < N else t[N] - t[N-1]
        off = (t[idx] - t[i]) / h
        key = (m, width, tuple(np.round(off, 9)))
        w = cache.get(key)
        if w is None:
            w = fd_weights(0.0, off, m)
            cache[key] = w
        return idx, w / h**m

    r1, c1, v1, r2, c2, v2 = [], [], [], [], [], []
    h2 = p // 2
    for i in range(N + 1):
        s = min(max(i - h2 + int(shifts[i]), 0), N - p)
        idx, w = weights(i, s, p + 1, 1)
        r1 += [i]*len(idx); c1 += list(idx); v1 += list(w)
        if h2 <= i <= N - h2:
            s, width = i - h2, p + 1
        elif i < h2:
            s, width = 0, p + 2
        else:
            s, width = N - p - 1, p + 2
        idx, w = weights(i, s, width, 2)
        r2 += [i]*len(idx); c2 += list(idx); v2 += list(w)
    D1 = sp.csr_matrix((v1, (r1, c1)), shape=(N+1, N+1))
    D2 = sp.csr_matrix((v2, (r2, c2)), shape=(N+1, N+1))
    return D1, D2




###############################################################################
#                                                                             #
#   2. OPERATOR MATRICES D1, D2  (SPARSE BANDED VANDERMONDE-TYPE ROWS)
#   
#   Row i of D1 holds the (p+1)-point first-derivative weights, window shifted
#   by shifts[i] in {-1,0,+1} (upwind), clipped to stay inside the mesh.
#   Row i of D2 holds the second-derivative weights: central p+1 points
#   in the interior, one-sided p+2 points in the first/last p/2 rows.
#   Then   y' ~ D1 @ Y   and   y'' ~ D2 @ Y.
#                                                                             #
###############################################################################


###############################################################################
#                                                                             #
#   3. PROBLEM CONTAINERS AND BOUNDARY-CONDITION HELPERS
#   
#   BVP        : generic description of a problem (f, df/dy, df/dy', BCs, weight).
#   HOFiDResult: everything the solver returns.
#   No concrete problem is defined in this file; run_all_problems.py builds them.
#                                                                             #
###############################################################################


def dirichlet_right(value: float):
    return lambda y, yp: (y - value, 1.0, 0.0)


def robin_right():
    """y(L) + y'(L) = 0"""
    return lambda y, yp: (y + yp, 1.0, 1.0)


@dataclass
class BVP:
    """w(t)(y'' - f(t,y,y')) = 0 on (0, L].  f, fy, fyp take arrays (t, y, y')."""
    name: str
    L: float
    f: Callable
    fy: Callable
    fyp: Callable
    left: tuple                      # ("dirichlet", v) | ("neumann", v)
    right: Callable                  # b(y, yp) -> (b, db/dy, db/dyp) scalars
    weight: Optional[Callable] = None
    exact: Optional[Callable] = None


@dataclass
class HOFiDResult:
    t: np.ndarray
    Y: np.ndarray
    err_node: np.ndarray
    err_est_max: float
    err_rel_est_max: float
    scaled_err_max: float
    converged: bool
    status: str
    orders: list
    tols: list
    t_initial: np.ndarray
    n_blocks: int
    refinements: int
    newton_iters: int
    runtime: float
    history: list = field(default_factory=list)
    Y_hi: Optional[np.ndarray] = None
    shifts: Optional[np.ndarray] = None
    p_final: int = 0
    tol_final: float = float("nan")
    termination: str = ""
    newton_converged: bool = True

    @property
    def n_points(self):
        return len(self.t)



###############################################################################
#                                                                             #
#   4. DAMPED NEWTON SOLVER  (NONLINEAR SYSTEM  F(Y) = 0)
#   
#   Residual rows : [ BC_left | w_i*(y''_i - f_i) | BC_right ].
#   Jacobian      : w*D2 - diag(fy) - diag(fyp)*D1 + BC rows.
#   Step damping  : lambda halved until |dY_new| <= (1 - lambda/4) |dY|.
#                                                                             #
###############################################################################


def _upwind_shifts(prob: BVP, t, Y):
    yp = np.gradient(Y, t)
    sh = np.zeros(len(t), dtype=int)
    with np.errstate(all="ignore"):
        fyp = prob.fyp(t[1:-1], Y[1:-1], yp[1:-1])
    sh[1:-1] = np.where(fyp < 0, -1, np.where(fyp > 0, 1, 0))
    return sh


def newton_solve(prob: BVP, t, p, Y0, shifts, tol=1e-10, maxit=60):
    """Damped Newton. Returns (Y, converged, iterations, last |dY|)."""
    N = len(t) - 1
    D1, D2 = build_operators(t, p, shifts)
    ti = t[1:-1]
    w = prob.weight(ti) if prob.weight is not None else np.ones_like(ti)
    wf = np.zeros(N + 1); wf[1:-1] = w
    kind, v = prob.left
    E0 = sp.csr_matrix(([1.0], ([0], [0])), shape=(N+1, N+1))
    EN = sp.csr_matrix(([1.0], ([N], [N])), shape=(N+1, N+1))
    S0D1, SND1 = E0 @ D1, EN @ D1
    left_rows = E0 if kind == "dirichlet" else S0D1

    def residual(Y):
        with np.errstate(all="ignore"):
            yp, ypp = D1 @ Y, D2 @ Y
            R = np.empty(N + 1)
            R[0] = (Y[0] - v) if kind == "dirichlet" else (yp[0] - v)
            R[1:-1] = w * (ypp[1:-1] - prob.f(ti, Y[1:-1], yp[1:-1]))
            R[-1] = prob.right(Y[-1], yp[-1])[0]
        return R, yp

    def jacobian(Y, yp):
        with np.errstate(all="ignore"):
            dy = np.zeros(N+1); dy[1:-1] = w * prob.fy(ti, Y[1:-1], yp[1:-1])
            dyp = np.zeros(N+1); dyp[1:-1] = w * prob.fyp(ti, Y[1:-1], yp[1:-1])
        _, by, byp = prob.right(Y[-1], yp[-1])
        J = sp.diags(wf) @ D2 - sp.diags(dy) - sp.diags(dyp) @ D1
        return (J + left_rows + by * EN + byp * SND1).tocsc()

    Y = Y0.astype(float).copy()
    last = np.inf
    for it in range(1, maxit + 1):
        R, yp = residual(Y)
        if not np.all(np.isfinite(R)):
            return Y, False, it, last
        try:
            lu = spla.splu(jacobian(Y, yp))
        except RuntimeError:
            return Y, False, it, last
        dY = -lu.solve(R)
        n0 = np.linalg.norm(dY, np.inf)
        last = n0
        if not np.isfinite(n0):
            return Y, False, it, last
        if n0 <= tol * (1.0 + np.linalg.norm(Y, np.inf)):
            return Y + dY, True, it, n0
        lam = 1.0
        while True:
            Yn = Y + lam * dY
            Rn, _ = residual(Yn)
            if np.all(np.isfinite(Rn)):
                dn = -lu.solve(Rn)
                if np.linalg.norm(dn, np.inf) <= (1.0 - lam / 4.0) * n0:
                    break
            lam *= 0.5
            if lam < 1e-4:
                ok = n0 <= 1e-7 * (1.0 + np.linalg.norm(Y, np.inf))
                return Y + (dY if ok else 0.0), ok, it, n0
        Y = Yn
    return Y, False, maxit, last



###############################################################################
#                                                                             #
#   5. MESH MANAGEMENT  (EQUIDISTANT BLOCKS, GRADING, ADAPTATION)
#   
#   A mesh is stored as blocks (B, n): block k spans [B_k, B_{k+1}] with n_k
#   equal intervals.  adapt_blocks splits blocks, rescales n_k from the
#   error indicator and restricts neighbouring step ratios to q = 2.
#                                                                             #
###############################################################################


def blocks_from_mesh(t, rtol=1e-8):
    h = np.diff(t)
    B, n = [t[0]], []
    cnt = 1
    for i in range(1, len(h)):
        if abs(h[i] - h[i-1]) <= rtol * h[i-1]:
            cnt += 1
        else:
            n.append(cnt); B.append(t[i]); cnt = 1
    n.append(cnt); B.append(t[-1])
    return np.array(B), np.array(n, dtype=int)


def build_mesh(B, n):
    parts = [np.linspace(B[k], B[k+1], n[k]+1)[:-1] for k in range(len(n))]
    return np.concatenate(parts + [[B[-1]]])


def _merge(B, n):
    Bo, no = [B[0]], []
    cur_n, cur_start = n[0], B[0]
    for k in range(1, len(n)):
        h_prev = (B[k] - cur_start) / cur_n
        h_k = (B[k+1] - B[k]) / n[k]
        if abs(h_prev - h_k) <= 1e-9 * h_prev:
            cur_n += n[k]
        else:
            no.append(cur_n); Bo.append(B[k]); cur_n, cur_start = n[k], B[k]
    no.append(cur_n); Bo.append(B[-1])
    return np.array(Bo), np.array(no, dtype=int)


def enforce_block_minimum(B, n, p):
    return B, np.maximum(n, p + 3)


def adapt_blocks(B, n, t, s, p, q=2.0, target=0.5):
    nmin = p + 3
    pieces, lo = [], 0
    for k in range(len(n)):
        hi = lo + n[k]
        split = False
        if n[k] >= 2 * nmin:
            mid = lo + n[k] // 2
            sl, sr = s[lo:mid+1].max(), s[mid:hi+1].max()
            if max(sl, sr) > target and max(sl, sr) > 4.0 * min(sl, sr):
                pieces += [(t[lo], t[mid], mid-lo, lo, mid),
                           (t[mid], t[hi], hi-mid, mid, hi)]
                split = True
        if not split:
            pieces.append((t[lo], t[hi], n[k], lo, hi))
        lo = hi
    qe = max(2.0, p / 2.0)
    a = np.array([pc[0] for pc in pieces]); b = np.array([pc[1] for pc in pieces])
    nn = np.empty(len(pieces), dtype=int)
    for j, (_, _, cnt, l, h) in enumerate(pieces):
        sm = s[l:h+1].max()
        if sm > target:
            f = min((sm / target)**(1.0 / qe), 8.0)
        elif sm < 0.1:
            f = max((max(sm, 1e-300) / target)**(1.0 / qe), 0.5)
        else:
            f = 1.0
        nn[j] = max(nmin, int(np.ceil(cnt * f)))
    wdt = b - a
    for _ in range(200):
        hh = wdt / nn
        changed = False
        for k in range(len(nn)):
            for j in (k-1, k+1):
                if 0 <= j < len(nn) and hh[k] > q * hh[j] * (1 + 1e-12):
                    nn[k] = int(np.ceil(wdt[k] / (q * hh[j]))); hh[k] = wdt[k] / nn[k]
                    changed = True
        if not changed:
            break
    Bn = np.concatenate([a, [b[-1]]])
    return _merge(Bn, nn)



###############################################################################
#                                                                             #
#   6. ADAPTIVE DRIVER  (ERROR ESTIMATE |Y_P - Y_{P+2}|, REFINEMENT LOOP)
#   
#   For each stage (p, tol):  solve with order p and p+2 on the same mesh,
#   estimate the error, stop if  max |Y_p - Y_{p+2}| / (tol (1+|Y|)) <= 1,
#   otherwise adapt the mesh and repeat (at most maxref refinements).
#                                                                             #
###############################################################################


def _local_indicator(prob, t, Yp, p):
    N = len(t) - 1
    sh = _upwind_shifts(prob, t, Yp)
    D1, D2 = build_operators(t, p + 2, sh)
    ti = t[1:-1]
    w = prob.weight(ti) if prob.weight is not None else np.ones_like(ti)
    with np.errstate(all="ignore"):
        yp, ypp = D1 @ Yp, D2 @ Yp
        r = np.zeros(N + 1)
        r[1:-1] = np.abs(w * (ypp[1:-1] - prob.f(ti, Yp[1:-1], yp[1:-1])))
    r = np.where(np.isfinite(r), r, 0.0)
    r[1:-1] = np.maximum(np.maximum(r[:-2], r[1:-1]), r[2:])
    r[0], r[-1] = r[1], r[-2]
    return r


def _solve_pair(prob, t, p, Y0):
    shifts = _upwind_shifts(prob, t, Y0)
    Yp, ok1, it1, _ = newton_solve(prob, t, p, Y0, shifts)
    if not ok1:
        return None, None, it1, f"Newton failed (order {p})", shifts
    Yh, ok2, it2, _ = newton_solve(prob, t, p + 2, Yp, shifts)
    if not ok2:
        return None, None, it1 + it2, f"Newton failed (order {p+2}, error estimator)", shifts
    return Yp, Yh, it1 + it2, "ok", shifts


def hofid_solve(prob: BVP, orders: Sequence[int], tols: Sequence[float], t0,
                y0fun: Callable, q: Optional[float] = None, maxref: int = 25,
                verbose: bool = False, max_points: int = 60000) -> HOFiDResult:
    """Adaptive solve. Returns the best (last successful) mesh/solution."""
    tic = time.time()
    q = 2.0 if q is None else float(q)
    orders, tols = list(orders), list(tols)
    assert len(orders) == len(tols) and all(o % 2 == 0 and o >= 2 for o in orders)
    t_init = np.asarray(t0, float)
    t = t_init.copy()
    Y = np.asarray(y0fun(t), float)
    B, n = blocks_from_mesh(t)
    hist, total_newton, nref = [], 0, 0
    status, termination, converged = "ok", "not started", False
    good = None
    for stage, (p, tol) in enumerate(zip(orders, tols)):
        B, n = enforce_block_minimum(B, n, p)
        B, n = _merge(B, n)
        t_new = build_mesh(B, n)
        if len(t_new) != len(t) or not np.allclose(t_new, t):
            Y = CubicSpline(t, Y)(t_new); t = t_new
        converged = False
        for k in range(maxref + 1):
            Yp, Yh, nit, msg, sh = _solve_pair(prob, t, p, Y)
            total_newton += nit
            if Yp is None:
                status = f"stage {stage} (p={p}, tol={tol:g}): {msg}"
                termination = "Newton failure"
                break
            Y = Yp
            err = np.abs(Yp - Yh)
            s = err / (tol * (1.0 + np.abs(Yp)))
            smax = s.max()
            good = dict(t=t.copy(), Y=Yp.copy(), Yh=Yh.copy(), sh=sh.copy(),
                        err=err, s=s, p=p, tol=tol, nb=len(n))
            hist.append(dict(stage=stage, order=p, tol=tol, iteration=k,
                             points=len(t), blocks=len(n),
                             est_abs_err=err.max(), scaled_err=smax,
                             newton_iters=nit))
            if verbose:
                print(f"  p={p:2d} tol={tol:7.1e} it={k:2d} N+1={len(t):6d} "
                      f"blocks={len(n):3d} max est err={err.max():.3e} scaled={smax:.3e}")
            if smax <= 1.0:
                converged, status, termination = True, "ok", "tolerance met"
                break
            if k == maxref:
                status = f"stage {stage} (p={p}): max refinements ({maxref}) reached"
                termination = "max refinements"
                break
            loc = _local_indicator(prob, t, Yp, p)
            s_eff = smax * loc / max(loc.max(), 1e-300)
            Bn, nn = adapt_blocks(B, n, t, s_eff, p, q)
            if len(Bn) == len(B) and np.allclose(Bn, B) and np.array_equal(nn, n):
                lo, nn = 0, n.copy()
                for j in range(len(n)):
                    if s[lo:lo+n[j]+1].max() > 1.0:
                        nn[j] = int(np.ceil(1.5 * n[j]))
                    lo += n[j]
                Bn = B
            if int(nn.sum()) + 1 > max_points:
                status = f"stage {stage} (p={p}): mesh limit {max_points} exceeded"
                termination = "mesh limit"
                break
            B, n = Bn, nn
            t_new = build_mesh(B, n)
            Y = CubicSpline(t, Y)(t_new); t = t_new
            nref += 1
        if not converged:
            break
    if good is not None:
        gt, gY, gerr, gs = good["t"], good["Y"], good["err"], good["s"]
        rel = float(np.max(gerr / (1.0 + np.abs(gY))))
        return HOFiDResult(
            t=gt, Y=gY, err_node=gerr, err_est_max=float(gerr.max()),
            err_rel_est_max=rel, scaled_err_max=float(gs.max()),
            converged=converged, status=status, orders=orders, tols=tols,
            t_initial=t_init, n_blocks=good["nb"], refinements=nref,
            newton_iters=total_newton, runtime=time.time() - tic, history=hist,
            Y_hi=good["Yh"], shifts=good["sh"], p_final=good["p"],
            tol_final=good["tol"], termination=termination,
            newton_converged=termination != "Newton failure")
    nan = np.full_like(Y, np.nan)
    return HOFiDResult(
        t=t_init, Y=Y, err_node=nan, err_est_max=np.nan, err_rel_est_max=np.nan,
        scaled_err_max=np.nan, converged=False, status=status, orders=orders,
        tols=tols, t_initial=t_init, n_blocks=len(n), refinements=nref,
        newton_iters=total_newton, runtime=time.time() - tic, history=hist,
        termination=termination, newton_converged=False)


def hofid_bvp(prob, orders, tols, t0, y0fun, q=None, maxref=25, verbose=False):
    """MATLAB-like interface. Returns (t, Y, err_max, err_node)."""
    r = hofid_solve(prob, orders, tols, t0, y0fun, q=q, maxref=maxref, verbose=verbose)
    return r.t, r.Y, r.err_est_max, r.err_node



###############################################################################
#                                                                             #
#   7. DIAGNOSTICS  (RESIDUALS, ERROR MEASURES, STENCIL COEFFICIENT TABLE)
#   
#   Pure post-processing of a finished HOFiDResult: y', y'', f, residuals,
#   estimated/actual errors and the FD weights actually used at each node.
#                                                                             #
###############################################################################


REL_EPS = 1e-12
REL_VALID_FRAC = 1e-3


def stencil_labels(i, N, p, shift):
    s_raw = i - p // 2 + int(shift)
    if s_raw < 0:
        l1, r1 = "initial one-sided (forward), p+1 pts", \
                  "window would leave the mesh on the left -> one-sided forward formula"
    elif s_raw > N - p:
        l1, r1 = "final one-sided (backward), p+1 pts", \
                  "window would leave the mesh on the right -> one-sided backward formula"
    elif shift < 0:
        l1, r1 = "upwind backward-biased, p+1 pts", \
                  "df/dy' < 0 (c=-df/dy' > 0) -> window shifted one node backward"
    elif shift > 0:
        l1, r1 = "upwind forward-biased, p+1 pts", \
                  "df/dy' > 0 (c < 0) -> window shifted one node forward"
    else:
        l1, r1 = "central, p+1 pts", "df/dy' = 0 -> no upwinding, central window"
    h2 = p // 2
    l2 = ("central (main scheme), p+1 pts" if h2 <= i <= N - h2 else
          "initial one-sided, p+2 pts" if i < h2 else "final one-sided, p+2 pts")
    return l1, r1, l2


def compute_diagnostics(prob: BVP, res: HOFiDResult) -> dict:
    t, Y, p, sh = res.t, res.Y, res.p_final, res.shifts
    N = len(t) - 1
    D1, D2 = build_operators(t, p, sh)
    yp, ypp = D1 @ Y, D2 @ Y
    ti = t[1:-1]
    w = np.zeros(N + 1)
    w[1:-1] = prob.weight(ti) if prob.weight is not None else 1.0
    f = np.full(N + 1, np.nan)
    R = np.full(N + 1, np.nan)
    Rw = np.full(N + 1, np.nan)
    f[1:-1] = prob.f(ti, Y[1:-1], yp[1:-1])
    R[1:-1] = ypp[1:-1] - f[1:-1]
    Rw[1:-1] = w[1:-1] * R[1:-1]
    kind, v = prob.left
    bc = np.full(N + 1, np.nan)
    bc[0] = (Y[0] - v) if kind == "dirichlet" else (yp[0] - v)
    b, by, byp = prob.right(Y[-1], yp[-1])
    bc[-1] = b
    est = np.abs(res.Y_hi - Y) if res.Y_hi is not None else np.full(N + 1, np.nan)
    eta = est / (1.0 + np.abs(Y))
    d = dict(
        t=t, Y=Y, Y_hi=res.Y_hi, yp=yp, ypp=ypp, f=f,
        residual=R, residual_weighted=Rw, bc_residual=bc,
        est_abs_err=est, eta=eta,
        mixed_scaled=est / (res.tol_final * (1 + np.abs(Y))),
        D1=D1, D2=D2, shifts=sh, p=p, N=N, weight=w,
        right_bc=(b, by, byp),
        bc_text=dict(
            left=(f"Y_0 - ({v:g}) = 0" if kind == "dirichlet"
                  else f"y'(t_0) - ({v:g}) = 0"),
            right=f"b(Y_N, y'_N) = {b:.3e}, db/dy = {by:g}, db/dy' = {byp:g}"),
        exact=None, actual_abs=None, actual_rel=None, rel_valid=None)
    if prob.exact is not None:
        ye = np.asarray(prob.exact(t), float)
        d["exact"] = ye
        d["actual_abs"] = np.abs(Y - ye)
        d["actual_rel"] = d["actual_abs"] / np.maximum(np.abs(ye), REL_EPS)
        d["rel_valid"] = np.abs(ye) >= REL_VALID_FRAC * np.max(np.abs(ye))
    return d


def summarize_diagnostics(d: dict, res: HOFiDResult) -> dict:
    nan = float("nan")
    out = dict(
        max_est_abs_err=float(np.nanmax(d["est_abs_err"])),
        max_scaled_eta=float(np.nanmax(d["eta"])),
        max_actual_abs_err=nan, max_actual_rel_err=nan,
        max_interior_residual=float(np.nanmax(np.abs(d["residual"]))),
        max_interior_residual_weighted=float(np.nanmax(np.abs(d["residual_weighted"]))),
        max_boundary_residual=float(np.nanmax(np.abs(d["bc_residual"]))),
        tolerance_requested=res.tol_final,
        tolerance_satisfied=bool(res.converged),
        newton_converged=bool(res.newton_converged),
        termination=res.termination)
    if d["actual_abs"] is not None:
        out["max_actual_abs_err"] = float(d["actual_abs"].max())
        m = d["rel_valid"]
        out["max_actual_rel_err"] = float(d["actual_rel"][m].max()) if m.any() else nan
    return out


def coefficient_table(d: dict, res: HOFiDResult, prob: BVP):
    t, p, N, sh = d["t"], d["p"], d["N"], d["shifts"]
    rows = []
    kind = prob.left[0]
    byp = d["right_bc"][2]
    for i in range(N + 1):
        l1, r1, l2 = stencil_labels(i, N, p, sh[i])
        hl = t[i] - t[i-1] if i > 0 else np.nan
        hr = t[i+1] - t[i] if i < N else np.nan
        for order, M, lab, reason in ((1, d["D1"], l1, r1), (2, d["D2"], l2, "central/one-sided by position")):
            lo, hi = M.indptr[i], M.indptr[i+1]
            idx, val = M.indices[lo:hi], M.data[lo:hi]
            o = np.argsort(idx); idx, val = idx[o], val[o]
            if 0 < i < N:
                used = "interior ODE"
            elif order == 2:
                used = "not used (BC row replaces ODE)"
            elif i == 0:
                used = "left BC y'(0)=v" if kind == "neumann" else "not used (Dirichlet BC)"
            else:
                used = "right BC b(y,y')" if byp != 0 else "not used (BC has no y' term)"
            for j, c in zip(idx, val):
                rows.append(dict(i=i, t_i=t[i], h_left=hl, h_right=hr,
                                 deriv_order=order, fd_order_p=p,
                                 stencil_type=lab, selection_reason=reason,
                                 used_in=used, upwind_shift=int(sh[i]),
                                 j=int(j), t_j=t[j], coefficient=float(c)))
    return rows



###############################################################################
#                                                                             #
#   8. BUILT-IN VALIDATION TESTS  (SELF-CONTAINED, USE TINY TEST PROBLEMS ONLY)
#                                                                             #
###############################################################################


def run_tests(verbose=True):
    log = print if verbose else (lambda *a, **k: None)
    rng = np.random.default_rng(0)

    # 1. derivative weights exact for polynomials up to degree p
    for p in (4, 6, 8, 10):
        t = np.linspace(0, 1, 40)**1.3
        D1, D2 = build_operators(t, p, np.zeros(len(t), int))
        for deg in range(0, p + 1):
            c = rng.standard_normal(deg + 1)
            P = np.polynomial.Polynomial(c)
            Y = P(t)
            assert np.max(np.abs(D1@Y - P.deriv(1)(t))) < 1e-6*(1+np.max(np.abs(P.deriv(1)(t)))), (p, deg, "D1")
            assert np.max(np.abs(D2@Y - P.deriv(2)(t))) < 1e-5*(1+np.max(np.abs(P.deriv(2)(t)))), (p, deg, "D2")
    log("[ok] derivative weights exact for polynomials of degree <= p (p=4..10)")

    # 1b. Fornberg weights == direct Vandermonde solve (small stencil)
    for m in (1, 2):
        x = np.array([-2., -1., 0., 1., 2., 3.]) * 0.1
        V = np.vander(x, increasing=True).T          # V[k, j] = x_j**k
        rhs = np.zeros(len(x)); rhs[m] = float(math.factorial(m))
        w_v = np.linalg.solve(V, rhs)
        assert np.max(np.abs(w_v - fd_weights(0.0, x, m))) < 1e-6 * np.max(np.abs(w_v))
    log("[ok] Fornberg weights match the Vandermonde system V w = m! e_m")

    # 2. regular linear BVP
    lin = BVP("sin", np.pi/2, f=lambda t,y,yp: -y, fy=lambda t,y,yp: -np.ones_like(y),
              fyp=lambda t,y,yp: np.zeros_like(y),
              left=("dirichlet", 0.0), right=dirichlet_right(1.0))
    r = hofid_solve(lin, [6], [1e-8], np.linspace(0, np.pi/2, 21), lambda t: 0*t)
    assert r.converged and abs(r.Y[0]) < 1e-14 and abs(r.Y[-1]-1) < 1e-14
    assert np.max(np.abs(r.Y - np.sin(r.t))) < 1e-7
    log(f"[ok] regular linear BVP: true error {np.max(np.abs(r.Y-np.sin(r.t))):.2e}")

    # 3. nonlinear Bratu
    th = 1.5
    for _ in range(200): th = np.sqrt(2.0)*np.cosh(th/4.0)
    ex = lambda x: -2*np.log(np.cosh((x-.5)*th/2)/np.cosh(th/4))
    br = BVP("bratu", 1.0, f=lambda t,y,yp: -np.exp(y), fy=lambda t,y,yp: -np.exp(y),
             fyp=lambda t,y,yp: np.zeros_like(y),
             left=("dirichlet", 0.0), right=dirichlet_right(0.0))
    r = hofid_solve(br, [8], [1e-8], np.linspace(0,1,21), lambda t: 0*t)
    e_val = np.max(np.abs(r.Y - ex(r.t)))
    assert r.converged and e_val < 1e-6, e_val
    log(f"[ok] nonlinear Bratu: true error {e_val:.2e}")

    # 4. singular endpoint
    seen = []
    def f1(t, y, yp):
        seen.append(t.min()); return yp/t + 2.0
    p1 = BVP("p1", 1.0, f=f1, fy=lambda t,y,yp: np.zeros_like(y),
             fyp=lambda t,y,yp: 1.0/t, left=("dirichlet", 0.0), right=dirichlet_right(0.0))
    ex1 = lambda t: np.where(t>0, t**2*np.log(np.maximum(t,1e-300)), 0.0)
    r = hofid_solve(p1, [6], [1e-6], np.linspace(0,1,21), lambda t: 0*t)
    assert min(seen) > 0 and r.converged
    log(f"[ok] singular endpoint: ODE never evaluated at t=0; true error {np.max(np.abs(r.Y-ex1(r.t))):.2e}")

    # 5. adaptation
    r = hofid_solve(p1, [4], [1e-8], np.linspace(0,1,21), lambda t: 0*t)
    assert r.converged and r.n_points > 21 and r.refinements > 0
    Bn, nn = blocks_from_mesh(r.t)
    assert nn.min() >= 4+3
    log(f"[ok] mesh adaptation: 21 -> {r.n_points} pts, {r.n_blocks} blocks")

    # 6. diagnostics
    p1.exact = ex1
    d = compute_diagnostics(p1, r)
    rows = coefficient_table(d, r, p1)
    import collections
    acc = collections.defaultdict(float)
    for c in rows:
        acc[(c["i"], c["deriv_order"])] += c["coefficient"] * r.Y[c["j"]]
    for i in [0, 1, r.n_points//2, r.n_points-1]:
        assert abs(acc[(i,1)] - d["yp"][i]) <= 1e-9*(1+abs(d["yp"][i])), ("D1 sum", i)
        assert abs(acc[(i,2)] - d["ypp"][i]) <= 1e-7*(1+abs(d["ypp"][i])), ("D2 sum", i)
    assert np.isnan(d["residual"][0]) and np.isnan(d["residual"][-1])
    sm = summarize_diagnostics(d, r)
    assert sm["max_boundary_residual"] < 1e-12
    assert sm["max_interior_residual_weighted"] < 1e-6
    log(f"[ok] diagnostics: coeff sums == D1@Y, D2@Y; max residual {sm['max_interior_residual']:.2e}; "
        f"actual abs err {sm['max_actual_abs_err']:.2e}")
    return True



if __name__ == "__main__":
    run_tests()
    print("all built-in tests passed")

