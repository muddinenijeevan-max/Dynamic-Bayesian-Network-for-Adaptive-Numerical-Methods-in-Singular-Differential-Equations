"""
hofid_bvp.py -- Python re-implementation of the HOFiD_bvp methodology.

Reference
---------
Hohenegger, Settanni, Weinmueller, Wolde (2024), "Numerical treatment of singular
ODEs using finite difference and collocation methods", Appl. Numer. Math. 205,
184-194, https://doi.org/10.1016/j.apnum.2024.07.002   (Section 2 = HOFiD_bvp).

Scalar second-order BVP on [0, L] (singularity allowed at t = 0)

    w(t) * ( y'' - f(t, y, y') ) = 0 ,   t in (0, L)
    left  BC : y(0) = v    ("dirichlet")   or   y'(0) = v   ("neumann")
    right BC : b(y(L), y'(L)) = 0

(w(t) == 1 except Problem 4, where the paper's transformed equation is
reproduced literally by choosing w(t) = (1-t)^6 / 4.)

What follows the paper
----------------------
* y'' and y' are discretised separately (no first-order system).
* y'' : central formula, p+1 points ("main scheme") in the interior; one-sided
  "initial/final" formulas with p+2 points (order p) near the ends.
* y'  : p+1 point formula, upwind-biased (stencil shifted one node backward/forward
  according to the sign of df/dy'); one-sided formulas at the ends.
* Neumann/Robin conditions use the one-sided first-derivative formula of order p.
* The ODE is never evaluated at t = 0 (nor at t = L): the BC rows replace it.
* Mesh is piecewise constant: blocks of equidistant points with >= p+4 points.
* Error estimate: |Y_p - Y_{p+2}| on the same mesh (two consecutive even orders)
  -> stopping decision.  Mesh PLACEMENT uses a local indicator (residual of the order
  p+2 scheme at Y_p) because the global difference is a smooth plateau near singularities.
* Error equidistribution + grading restriction -> new blocks; previous solution is
  interpolated to the new mesh; damped Newton for the algebraic system.
* Variable order strategy: orders/tolerances are processed as successive stages.

Documented deviations (see README): the weights are Fornberg weights computed on
the *actual* nodes (identical to diag(1/h^nu) M_nu inside a block); the block
refinement rule, grading factor and Newton details are our own because the paper
does not give them; the code of [33] is not available to us.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.interpolate import CubicSpline

__all__ = ["BVP", "HOFiDResult", "fd_weights", "build_operators", "hofid_solve",
           "hofid_bvp", "dirichlet_right", "robin_right", "run_tests"]


# --------------------------------------------------------------------------- #
# Finite-difference weights (Fornberg 1988, "Generation of finite difference
# formulas on arbitrarily spaced grids")
# --------------------------------------------------------------------------- #
def fd_weights(z: float, x: np.ndarray, m: int) -> np.ndarray:
    """Weights w_j such that  f^(m)(z) ~ sum_j w_j f(x_j)  (exact for polynomials
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
                    c[i, k] = c1 * (k * c[i - 1, k - 1] - c5 * c[i - 1, k]) / c2
                c[i, 0] = -c1 * c5 * c[i - 1, 0] / c2
            for k in range(mn, 0, -1):
                c[j, k] = (c4 * c[j, k] - k * c[j, k - 1]) / c3
            c[j, 0] = c4 * c[j, 0] / c3
        c1 = c2
    return c[:, m]


def build_operators(t: np.ndarray, p: int, shifts: np.ndarray):
    """Sparse matrices D1, D2 (already divided by h^nu) approximating y' and y''
    at all nodes of t for order p (p even).

    Row i of D2: central p+1 points if p/2 <= i <= N-p/2, otherwise the one-sided
    "initial/final" scheme with p+2 points (order p).  Row i of D1: p+1 points,
    window centred at i+shifts[i] (clipped to the mesh -> one-sided at the ends).
    Inside an equidistant block the weights are the classical constant stencils
    (= rows of the quasi-Toeplitz matrices M1, M2 of the paper, scaled by 1/h^nu).
    """
    N = len(t) - 1
    if N < p + 2:
        raise ValueError(f"mesh too coarse for order {p}: need N >= {p + 2}, got {N}")
    cache: dict = {}

    def weights(i, s, width, m):
        idx = np.arange(s, s + width)
        h = t[i + 1] - t[i] if i < N else t[N] - t[N - 1]
        off = (t[idx] - t[i]) / h
        key = (m, width, tuple(np.round(off, 9)))
        w = cache.get(key)
        if w is None:
            w = fd_weights(0.0, off, m)
            cache[key] = w
        return idx, w / h ** m

    r1, c1, v1, r2, c2, v2 = [], [], [], [], [], []
    h2 = p // 2
    for i in range(N + 1):
        # first derivative
        s = min(max(i - h2 + int(shifts[i]), 0), N - p)
        idx, w = weights(i, s, p + 1, 1)
        r1 += [i] * len(idx); c1 += list(idx); v1 += list(w)
        # second derivative
        if h2 <= i <= N - h2:
            s, width = i - h2, p + 1
        elif i < h2:
            s, width = 0, p + 2
        else:
            s, width = N - p - 1, p + 2
        idx, w = weights(i, s, width, 2)
        r2 += [i] * len(idx); c2 += list(idx); v2 += list(w)
    D1 = sp.csr_matrix((v1, (r1, c1)), shape=(N + 1, N + 1))
    D2 = sp.csr_matrix((v2, (r2, c2)), shape=(N + 1, N + 1))
    return D1, D2


# --------------------------------------------------------------------------- #
# Problem definition
# --------------------------------------------------------------------------- #
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
    err_node: np.ndarray             # |Y_p - Y_{p+2}| (ESTIMATED absolute error)
    err_est_max: float               # max of err_node
    err_rel_est_max: float           # max err_node/(1+|Y|)
    scaled_err_max: float            # max err_node/(tol(1+|Y|)); <= 1 means tol met
    converged: bool                  # tolerance of the LAST stage met, all Newton OK
    status: str
    orders: list
    tols: list
    t_initial: np.ndarray
    n_blocks: int
    refinements: int
    newton_iters: int
    runtime: float
    history: list = field(default_factory=list)

    @property
    def n_points(self):
        return len(self.t)


# --------------------------------------------------------------------------- #
# Damped Newton on the discrete system
# --------------------------------------------------------------------------- #
def _upwind_shifts(prob: BVP, t, Y):
    yp = np.gradient(Y, t)
    sh = np.zeros(len(t), dtype=int)
    with np.errstate(all="ignore"):
        fyp = prob.fyp(t[1:-1], Y[1:-1], yp[1:-1])
    # y'' + c y' with c = -f_y': c>0 -> backward stencil, c<0 -> forward stencil
    sh[1:-1] = np.where(fyp < 0, -1, np.where(fyp > 0, 1, 0))
    return sh


def newton_solve(prob: BVP, t, p, Y0, shifts, tol=1e-10, maxit=60):
    """Damped Newton (Deuflhard natural-monotonicity test, LU re-used for the
    simplified correction).  Returns (Y, converged, iterations, last |dY|)."""
    N = len(t) - 1
    D1, D2 = build_operators(t, p, shifts)
    ti = t[1:-1]
    w = prob.weight(ti) if prob.weight is not None else np.ones_like(ti)
    wf = np.zeros(N + 1); wf[1:-1] = w
    kind, v = prob.left
    E0 = sp.csr_matrix(([1.0], ([0], [0])), shape=(N + 1, N + 1))
    EN = sp.csr_matrix(([1.0], ([N], [N])), shape=(N + 1, N + 1))
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
            dy = np.zeros(N + 1); dy[1:-1] = w * prob.fy(ti, Y[1:-1], yp[1:-1])
            dyp = np.zeros(N + 1); dyp[1:-1] = w * prob.fyp(ti, Y[1:-1], yp[1:-1])
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
            return Y + dY, True, it, n0          # quadratically small correction
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
                # roundoff-level stagnation is accepted only if already tiny
                ok = n0 <= 1e-7 * (1.0 + np.linalg.norm(Y, np.inf))
                return Y + (dY if ok else 0.0), ok, it, n0
        Y = Yn
    return Y, False, maxit, last


# --------------------------------------------------------------------------- #
# Piecewise-constant mesh handling and error-equidistribution adaptation
# --------------------------------------------------------------------------- #
def blocks_from_mesh(t, rtol=1e-8):
    """Group consecutive equal steps -> (breakpoints B, interval counts n)."""
    h = np.diff(t)
    B, n = [t[0]], []
    cnt = 1
    for i in range(1, len(h)):
        if abs(h[i] - h[i - 1]) <= rtol * h[i - 1]:
            cnt += 1
        else:
            n.append(cnt); B.append(t[i]); cnt = 1
    n.append(cnt); B.append(t[-1])
    return np.array(B), np.array(n, dtype=int)


def build_mesh(B, n):
    parts = [np.linspace(B[k], B[k + 1], n[k] + 1)[:-1] for k in range(len(n))]
    return np.concatenate(parts + [[B[-1]]])


def _merge(B, n):
    Bo, no = [B[0]], []
    cur_n, cur_start = n[0], B[0]
    for k in range(1, len(n)):
        h_prev = (B[k] - cur_start) / cur_n
        h_k = (B[k + 1] - B[k]) / n[k]
        if abs(h_prev - h_k) <= 1e-9 * h_prev:
            cur_n += n[k]
        else:
            no.append(cur_n); Bo.append(B[k]); cur_n, cur_start = n[k], B[k]
    no.append(cur_n); Bo.append(B[-1])
    return np.array(Bo), np.array(no, dtype=int)


def enforce_block_minimum(B, n, p):
    """Every block keeps >= p+4 points (>= p+3 intervals)."""
    return B, np.maximum(n, p + 3)


def adapt_blocks(B, n, t, s, p, q=2.0, target=0.5):
    """One adaptation sweep.
    s : node-wise scaled estimated error (s<=1 <=> tolerance met).
    1. blocks with >= 2(p+3) intervals whose two halves differ strongly in error
       are split at a node (mesh nodes unchanged);
    2. each block gets n_new = n * (s_max/target)^(1/q_eff)  (error equidistribution;
       q_eff = max(2, p/2) is a conservative order, since the paper reports order
       reduction to 2 near the singularity); mild coarsening where s_max < 0.1;
    3. grading restriction: neighbouring step sizes differ at most by factor q;
    4. neighbouring blocks with identical step are merged."""
    nmin = p + 3
    pieces, lo = [], 0
    for k in range(len(n)):
        hi = lo + n[k]
        split = False
        if n[k] >= 2 * nmin:
            mid = lo + n[k] // 2
            sl, sr = s[lo:mid + 1].max(), s[mid:hi + 1].max()
            if max(sl, sr) > target and max(sl, sr) > 4.0 * min(sl, sr):
                pieces += [(t[lo], t[mid], mid - lo, lo, mid),
                           (t[mid], t[hi], hi - mid, mid, hi)]
                split = True
        if not split:
            pieces.append((t[lo], t[hi], n[k], lo, hi))
        lo = hi
    qe = max(2.0, p / 2.0)
    a = np.array([pc[0] for pc in pieces]); b = np.array([pc[1] for pc in pieces])
    nn = np.empty(len(pieces), dtype=int)
    for j, (_, _, cnt, l, h) in enumerate(pieces):
        sm = s[l:h + 1].max()
        if sm > target:
            f = min((sm / target) ** (1.0 / qe), 8.0)
        elif sm < 0.1:
            f = max((max(sm, 1e-300) / target) ** (1.0 / qe), 0.5)
        else:
            f = 1.0
        nn[j] = max(nmin, int(np.ceil(cnt * f)))
    wdt = b - a
    for _ in range(200):                                   # grading restriction
        hh = wdt / nn
        changed = False
        for k in range(len(nn)):
            for j in (k - 1, k + 1):
                if 0 <= j < len(nn) and hh[k] > q * hh[j] * (1 + 1e-12):
                    nn[k] = int(np.ceil(wdt[k] / (q * hh[j]))); hh[k] = wdt[k] / nn[k]
                    changed = True
        if not changed:
            break
    Bn = np.concatenate([a, [b[-1]]])
    return _merge(Bn, nn)


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def _local_indicator(prob, t, Yp, p):
    """Local truncation-error indicator used ONLY to place mesh points: residual of the
    order-(p+2) scheme evaluated at the order-p solution (again a comparison of two
    consecutive orders).  The global estimate |Y_p - Y_(p+2)| is a smooth plateau
    (the singularity at t=0 pollutes the whole interval) and cannot localise refinement."""
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
    r[1:-1] = np.maximum(np.maximum(r[:-2], r[1:-1]), r[2:])      # remove stencil-zero oscillations
    r[0], r[-1] = r[1], r[-2]
    return r


def _solve_pair(prob, t, p, Y0):
    shifts = _upwind_shifts(prob, t, Y0)
    Yp, ok1, it1, _ = newton_solve(prob, t, p, Y0, shifts)
    if not ok1:
        return None, None, it1, f"Newton failed (order {p})"
    Yh, ok2, it2, _ = newton_solve(prob, t, p + 2, Yp, shifts)
    if not ok2:
        return None, None, it1 + it2, f"Newton failed (order {p + 2}, error estimator)"
    return Yp, Yh, it1 + it2, "ok"


def hofid_solve(prob: BVP, orders: Sequence[int], tols: Sequence[float], t0,
                y0fun: Callable, q: Optional[float] = None, maxref: int = 25,
                verbose: bool = False, max_points: int = 60000) -> HOFiDResult:
    """Adaptive solve.  orders/tols are processed as successive stages
    (variable-order strategy); a single order/tol pair is the fixed-order run."""
    tic = time.time()
    q = 2.0 if q is None else float(q)
    orders, tols = list(orders), list(tols)
    assert len(orders) == len(tols) and all(o % 2 == 0 and o >= 2 for o in orders)
    t_init = np.asarray(t0, float)
    t = t_init.copy()
    Y = np.asarray(y0fun(t), float)
    B, n = blocks_from_mesh(t)
    hist, total_newton, nref = [], 0, 0
    status, converged = "ok", False
    Yp = Y; err = np.zeros_like(Y); s = np.full_like(Y, np.inf)
    for stage, (p, tol) in enumerate(zip(orders, tols)):
        B, n = enforce_block_minimum(B, n, p)
        B, n = _merge(B, n)
        t_new = build_mesh(B, n)
        if len(t_new) != len(t) or not np.allclose(t_new, t):
            Y = CubicSpline(t, Y)(t_new); t = t_new
        converged = False
        for k in range(maxref + 1):
            Yp, Yh, nit, msg = _solve_pair(prob, t, p, Y)
            total_newton += nit
            if Yp is None:
                status = f"stage {stage} (p={p}, tol={tol:g}): {msg}"
                if k == 0 and stage == 0:
                    Yp, err, s = Y, np.full_like(Y, np.nan), np.full_like(Y, np.nan)
                break
            Y = Yp
            err = np.abs(Yp - Yh)
            s = err / (tol * (1.0 + np.abs(Yp)))
            smax = s.max()
            hist.append(dict(stage=stage, order=p, tol=tol, iteration=k, points=len(t),
                             blocks=len(n), est_abs_err=err.max(), scaled_err=smax,
                             newton_iters=nit))
            if verbose:
                print(f"  p={p:2d} tol={tol:7.1e} it={k:2d} N+1={len(t):6d} blocks={len(n):3d} "
                      f"max est err={err.max():.3e} scaled={smax:.3e}")
            if smax <= 1.0:
                converged, status = True, "ok"
                break
            if k == maxref:
                status = f"stage {stage} (p={p}): max refinements ({maxref}) reached, tol not met"
                break
            loc = _local_indicator(prob, t, Yp, p)
            s_eff = smax * loc / max(loc.max(), 1e-300)     # equidistribute the local indicator, scaled to the global need
            Bn, nn = adapt_blocks(B, n, t, s_eff, p, q)
            if len(Bn) == len(B) and np.allclose(Bn, B) and np.array_equal(nn, n):
                # stalled: force growth of blocks that violate the tolerance
                lo, nn = 0, n.copy()
                for j in range(len(n)):
                    if s[lo:lo + n[j] + 1].max() > 1.0:
                        nn[j] = int(np.ceil(1.5 * n[j]))
                    lo += n[j]
                Bn = B
            if int(nn.sum()) + 1 > max_points:
                status = f"stage {stage} (p={p}): mesh limit {max_points} points exceeded"
                break
            B, n = Bn, nn
            t_new = build_mesh(B, n)
            Y = CubicSpline(t, Y)(t_new); t = t_new
            nref += 1
        if not converged:
            break
    rel = float(np.nanmax(err / (1.0 + np.abs(Yp)))) if np.all(np.isfinite(err)) else np.nan
    return HOFiDResult(t=t, Y=Yp, err_node=err,
                       err_est_max=float(np.nanmax(err)) if np.any(np.isfinite(err)) else np.nan,
                       err_rel_est_max=rel,
                       scaled_err_max=float(np.nanmax(s)) if np.any(np.isfinite(s)) else np.nan,
                       converged=converged, status=status, orders=orders, tols=tols,
                       t_initial=t_init, n_blocks=len(n), refinements=nref,
                       newton_iters=total_newton, runtime=time.time() - tic, history=hist)


def hofid_bvp(prob, orders, tols, t0, y0fun, q=None, maxref=25, verbose=False):
    """MATLAB-like interface.  Returns (t, Y, err, e):
       t   final mesh,  Y  numerical solution (order of the last stage),
       err max ESTIMATED absolute error  max_i |Y_p - Y_{p+2}|  (NOT the true error),
       e   node-wise estimated absolute error.
    Use hofid_solve for convergence status / scaled error / history; check
    result.converged before trusting the tolerance."""
    r = hofid_solve(prob, orders, tols, t0, y0fun, q=q, maxref=maxref, verbose=verbose)
    return r.t, r.Y, r.err_est_max, r.err_node


# --------------------------------------------------------------------------- #
# Built-in validation
# --------------------------------------------------------------------------- #
def run_tests(verbose=True):
    log = print if verbose else (lambda *a, **k: None)
    rng = np.random.default_rng(0)
    # 1. derivative weights reproduce polynomials up to the expected degree
    for p in (4, 6, 8, 10):
        t = np.linspace(0, 1, 40) ** 1.3                     # smoothly graded mesh
        D1, D2 = build_operators(t, p, np.zeros(len(t), int))
        for deg in range(0, p + 1):
            c = rng.standard_normal(deg + 1)
            P = np.polynomial.Polynomial(c)
            Y = P(t)
            assert np.max(np.abs(D1 @ Y - P.deriv(1)(t))) < 1e-6 * (1 + np.max(np.abs(P.deriv(1)(t)))), (p, deg, "D1")
            if deg >= 2 or True:
                assert np.max(np.abs(D2 @ Y - P.deriv(2)(t))) < 1e-5 * (1 + np.max(np.abs(P.deriv(2)(t)))), (p, deg, "D2")
    log("[ok] derivative weights exact for polynomials of degree <= p (p = 4..10), incl. one-sided rows")
    # 2. regular linear BVP y'' = -y, y(0)=0, y(pi/2)=1  -> sin t  (Newton: 1 step + check)
    lin = BVP("sin", np.pi / 2, f=lambda t, y, yp: -y, fy=lambda t, y, yp: -np.ones_like(y),
              fyp=lambda t, y, yp: np.zeros_like(y), left=("dirichlet", 0.0),
              right=dirichlet_right(1.0))
    r = hofid_solve(lin, [6], [1e-8], np.linspace(0, np.pi / 2, 21), lambda t: 0 * t)
    assert r.converged and abs(r.Y[0]) < 1e-14 and abs(r.Y[-1] - 1) < 1e-14
    assert np.max(np.abs(r.Y - np.sin(r.t))) < 1e-7
    assert r.history[0]["newton_iters"] <= 6          # linear problem: 2-3 per solve
    log(f"[ok] regular linear BVP: true error {np.max(np.abs(r.Y-np.sin(r.t))):.2e}, "
        f"BCs satisfied, Newton iterations (order p + p+2) = {r.history[0]['newton_iters']}")
    # 3. nonlinear Bratu y'' = -exp(y), y(0)=y(1)=0
    th = 1.5
    for _ in range(200):
        th = np.sqrt(2.0) * np.cosh(th / 4.0)
        th = th  # fixed point iteration for the smaller root
    ex = lambda x: -2 * np.log(np.cosh((x - .5) * th / 2) / np.cosh(th / 4))
    br = BVP("bratu", 1.0, f=lambda t, y, yp: -np.exp(y), fy=lambda t, y, yp: -np.exp(y),
             fyp=lambda t, y, yp: np.zeros_like(y), left=("dirichlet", 0.0), right=dirichlet_right(0.0))
    r = hofid_solve(br, [8], [1e-8], np.linspace(0, 1, 21), lambda t: 0 * t)
    e = np.max(np.abs(r.Y - ex(r.t)))
    assert r.converged and e < 1e-6, e
    log(f"[ok] nonlinear (Bratu) damped Newton: converged, true error {e:.2e}")
    # 4. singular endpoint: Problem 1 -- f must never be called at t=0
    seen = []
    def f1(t, y, yp):
        seen.append(t.min()); return yp / t + 2.0
    p1 = BVP("p1", 1.0, f=f1, fy=lambda t, y, yp: np.zeros_like(y),
             fyp=lambda t, y, yp: 1.0 / t, left=("dirichlet", 0.0), right=dirichlet_right(0.0))
    ex1 = lambda t: np.where(t > 0, t ** 2 * np.log(np.maximum(t, 1e-300)), 0.0)
    r = hofid_solve(p1, [6], [1e-6], np.linspace(0, 1, 21), lambda t: 0 * t)
    assert min(seen) > 0 and r.converged and abs(r.Y[0]) < 1e-14 and abs(r.Y[-1]) < 1e-14
    log(f"[ok] singular endpoint: ODE never evaluated at t=0 (min t used {min(seen):.2e}); "
        f"true error {np.max(np.abs(r.Y-ex1(r.t))):.2e}")
    # 5. adaptive refinement changes the mesh when the estimate demands it
    r = hofid_solve(p1, [4], [1e-8], np.linspace(0, 1, 21), lambda t: 0 * t)
    assert r.converged and r.n_points > 21 and r.refinements > 0
    Bn, nn = blocks_from_mesh(r.t)
    assert nn.min() >= 4 + 3, "block minimum p+4 points violated"
    log(f"[ok] mesh adaptation: 21 -> {r.n_points} points in {r.n_blocks} blocks, "
        f"est. error {r.err_est_max:.2e}, all blocks >= p+4 points")
    # 6. error estimator is an estimate, not the true error: both reported
    log(f"     (estimated {r.err_est_max:.2e}  vs  true {np.max(np.abs(r.Y-ex1(r.t))):.2e})")
    return True


if __name__ == "__main__":
    run_tests()
    print("all built-in tests passed")
