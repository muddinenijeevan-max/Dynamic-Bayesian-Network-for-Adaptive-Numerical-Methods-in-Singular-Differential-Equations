function out = dbn_controller(cmd, varargin)
% DBN_CONTROLLER  Dynamic Bayesian Network that decides the next mesh-adaptation
% action of bvpsuite2.0 (pure inference code, no toolbox needed).
%
%   dbn = dbn_controller('init', N, m)
%   [dbn, dec] = dbn_controller('step', dbn, obs)
%
% ---------------------------------------------------------------------------
% Model (two time slices, exact forward filtering)
% ---------------------------------------------------------------------------
% Hidden variables at iteration t
%   S_t  (mesh state)         1 Misplaced     error >> TOL and unevenly distributed
%                                             -> moving points is cheaper than adding points
%                             2 Underresolved error > TOL, already well distributed
%                                             -> more points are needed
%                             3 Adequate      error moderately below TOL
%                             4 Overresolved  error far below TOL (points can be removed)
%   R_t  (estimator regime)   1 Pre-asymptotic  error estimate / order model not yet reliable
%                             2 Asymptotic      error ~ C*h^p holds, extrapolating N is safe
%
% Observations (all computed from the current and the previous iteration)
%   f_q = log10( max scaled error estimate / TOL )                 -> S
%   f_s = log10( max/mean of the per-interval error estimate )     -> S   (solution/error variation)
%   f_d = log10( max/mean of the per-interval defect (residual) )  -> S   (defect)
%   f_t = log10( q_t / q_{t-1} ),  only after a relocation step     -> S   (previous iteration's effect)
%   f_o = log10( p_obs / p_nominal ),  after an N-changing step     -> R   (observed convergence order)
% Each is Gaussian given the hidden state (parameters below, set from domain
% knowledge of collocation/mesh adaptation; they are NOT fitted to the paper's problems).
%
% Dynamics: P(S_t,R_t | S_{t-1},R_{t-1}, a_{t-1}) = P(S_t|S_{t-1},a_{t-1}) * P(R_t|R_{t-1},a_{t-1}),
% i.e. the action taken at the previous iteration moves the mesh state.
%
% Decision: expected utility  EU(a) = sum_s b(s) U(a,s)  over the actions that are
% admissible (refine/relocate only if tolerance is violated, accept/coarsen only if
% it is satisfied -- so the DBN can never accept a solution the error estimator rejects).

switch cmd
    case 'init'
        out = dbn_init(varargin{:});
    case 'step'
        [out.dbn, out.dec] = dbn_step(varargin{:});
    otherwise
        error('dbn_controller: unknown command %s', cmd);
end
end

% =========================================================================
function dbn = dbn_init(N, m)
P = dbn_params();
dbn.P = P;
dbn.m = m;
dbn.pnom = m + 1;                 % nominal order of the max-norm error of an m-point collocation scheme
dbn.t = 0;
dbn.prevAction = 0;
dbn.prev = [];                    % previous iteration's observations
dbn.stagn = 0;                    % consecutive relocation passes that did not pay off
% priors over S and R, joint belief b(s,r) as 4x2 matrix
pS = P.priorS(:);
if N >= 40, pR = [0.35 0.65]; else, pR = [0.65 0.35]; end
dbn.belief = pS * pR;
end

% =========================================================================
function [dbn, dec] = dbn_step(dbn, obs)
% obs fields: q, skewE, skewD, N (intervals), qTol-valid flag, canCoarsen
P = dbn.P;

% ---- prediction: apply the action-conditioned transition model --------
if dbn.t > 0
    a = dbn.prevAction;
    b = P.TS{a}' * dbn.belief * P.TR{a};      % (S' x R')
else
    b = dbn.belief;
end

% ---- correction: multiply by the emission likelihoods ------------------
ll = zeros(4,2);                                % log-likelihood per (s,r)
fq = min(max(log10(max(obs.q, 1e-12)), -4), 4);   % q > 1e4 is treated alike (heavy-tail guard)
ll = ll + P.w.q * repmat(gl(fq, P.q.mu, P.q.sd), 1, 2);
fs = log10(max(obs.skewE, 1));
ll = ll + P.w.s * repmat(gl(fs, P.s.mu, P.s.sd), 1, 2);
fd = log10(max(obs.skewD, 1));
ll = ll + P.w.d * repmat(gl(fd, P.d.mu, P.d.sd), 1, 2);

pobs = NaN;  have_o = false;  have_t = false;
if dbn.t > 0
    prev = dbn.prev;
    % trend after a relocation (same N): informative about S
    if abs(obs.N - prev.N) <= 1 && prev.q > 0 && obs.q > 0
        ft = log10(obs.q / prev.q);
        ll = ll + P.w.t * repmat(gl(ft, P.t.mu, P.t.sd), 1, 2);
        have_t = true;
    end
    % observed order after an N-changing step: informative about R
    if obs.N > 1.05*prev.N || obs.N < prev.N/1.05
        qa = max(prev.q,1e-9);  qb = max(obs.q,1e-9);
        pobs = log(qa/qb) / log(obs.N/prev.N);
        % sign convention: N up & q down -> positive; N down & q up -> positive
        if pobs > 0.2 && min(qa,qb) > 1e-7
            fo = log10(min(pobs, 3*dbn.pnom) / dbn.pnom);
            ll = ll + repmat(gl(fo, P.o.mu, P.o.sd)', 4, 1);   % depends on R only
            have_o = true;
        else
            pobs = NaN;
        end
    end
end
b = b .* exp(ll - max(ll(:)));
b = b / sum(b(:));
dbn.belief = b;

bS = sum(b, 2)';  bR = sum(b, 1);

% ---- decision ----------------------------------------------------------
tolOK = obs.q < 1;
EU = P.U * bS';                                % 4 actions: relocate refine coarsen accept
valid = false(4,1);
if tolOK
    valid([3 4]) = true;
    if ~obs.canCoarsen, valid(3) = false; end
else
    valid([1 2]) = true;
end
EUv = EU;  EUv(~valid) = -Inf;
[~, act] = max(EUv);

% guard: relocation passes that no longer pay off (error reduced by less than 2x) are not repeated
if dbn.t > 0 && dbn.prevAction == 1 && ~tolOK
    if obs.q > P.stagnRatio * dbn.prev.q, dbn.stagn = dbn.stagn + 1; else, dbn.stagn = 0; end
else
    dbn.stagn = 0;
end
if ~tolOK && act == 1 && dbn.stagn >= 1
    act = 2;
end

% order used to extrapolate the needed number of intervals
pnom = dbn.pnom;
if ~isnan(pobs) && bR(2) > 0.5
    peff = min(max(pobs, 1), 2*pnom);            % trust measured order
elseif ~isnan(pobs)
    peff = min(max(pobs, 1), pnom);              % pre-asymptotic: never be more optimistic than nominal
else
    peff = pnom;
end
if bR(1) > 0.5, gmax = P.gmaxPre; else, gmax = P.gmaxAsy; end
if obs.q > 1e3, gmax = 1.75*gmax; end            % error thousands of times the tolerance: doubling cannot be enough

dec.action = act;                                % 1 relocate 2 refine 3 coarsen 4 accept
dec.relocate = (act == 1) || (act == 2 && bS(1) > P.relocThresh);
dec.peff = peff;
dec.gmax = gmax;
dec.bS = bS;  dec.bR = bR;
dec.pobs = pobs;
dec.have_o = have_o;  dec.have_t = have_t;
dec.EU = EU';

% mesh-size proposal (number of intervals)
Mi = obs.N;
switch act
    case 2   % refine: q -> P.qTarget
        g = (obs.q / P.qTarget)^(1/peff);
        g = min(max(g, P.gmin), gmax);
        dec.Mnew = max(Mi + 2, ceil(g * Mi));
    case 3   % coarsen: only as far as the (conservative) model predicts
        qe = max(obs.q, P.qFloor);
        g = (qe / P.qTarget)^(1/pnom);
        g = min(max(g, P.gcoarseMin), 0.9);
        dec.Mnew = max(P.Nmin, floor(g * Mi));
    otherwise
        dec.Mnew = Mi;
end

dbn.t = dbn.t + 1;
dbn.prevAction = act;
dbn.prev = struct('q', obs.q, 'N', obs.N);
end

% =========================================================================
function ll = gl(x, mu, sd)
% Gaussian log-likelihood (up to a constant) for each entry of mu/sd
ll = -0.5*((x - mu(:)')./sd(:)').^2 - log(sd(:)');
ll = ll(:);
end

% =========================================================================
function P = dbn_params()
% ---------------- priors / emissions (S = Misplaced, Underres., Adequate, Overres.) ----
P.priorS = [0.35 0.40 0.17 0.08];
P.q.mu = [ 0.9  1.6 -0.6 -3.0];   P.q.sd = [0.7 1.2 0.45 1.3];
P.s.mu = [ 0.9  0.30 0.45  0.50]; P.s.sd = [0.35 0.22 0.45 0.55];
P.d.mu = [ 0.8  0.30 0.40  0.45]; P.d.sd = [0.40 0.25 0.45 0.55];
P.t.mu = [-0.9 -0.10 -0.30  0.00]; P.t.sd = [0.40 0.25 0.60 0.80];
% evidence weights (tempering): f_q, f_s, f_d are strongly correlated, so each counts less than
% once; the trend f_t (did the last action pay off?) is the most direct evidence and counts more
P.w.q = 0.6;  P.w.s = 0.5;  P.w.d = 0.4;  P.w.t = 1.5;
% R = Pre-asymptotic, Asymptotic
P.o.mu = [-0.45 0.00];            P.o.sd = [0.50 0.22];

% ---------------- transitions, indexed by the previous action ----------
%            to:  M     U     A     O
TSrel = [ 0.25  0.50  0.20  0.05      % from Misplaced     (a = relocate)
          0.08  0.77  0.15  0.00      % from Underresolved
          0.00  0.10  0.80  0.10      % from Adequate
          0.00  0.00  0.20  0.80 ];   % from Overresolved
TSref = [ 0.35  0.35  0.28  0.02      % a = refine (also relocates when Misplaced)
          0.05  0.30  0.60  0.05
          0.00  0.05  0.85  0.10
          0.00  0.00  0.20  0.80 ];
TScoa = [ 0.20  0.40  0.35  0.05      % a = coarsen
          0.05  0.45  0.45  0.05
          0.00  0.30  0.65  0.05
          0.00  0.05  0.45  0.50 ];
TSacc = eye(4);
P.TS = {TSrel, TSref, TScoa, TSacc};
%            to:  pre   asym
TRrel = [ 0.85 0.15; 0.10 0.90 ];            % relocation does not change the regime much
TRref = [ 0.45 0.55; 0.05 0.95 ];            % adding points moves towards the asymptotic regime
TRcoa = [ 0.90 0.10; 0.30 0.70 ];            % coarsening can push back into the pre-asymptotic regime
P.TR = {TRrel, TRref, TRcoa, eye(2)};

% ---------------- utilities U(a,s) ------------------------------------
%                  M     U     A     O
P.U = [  0.60 -0.30 -0.50 -0.80      % relocate
         0.25  0.70 -0.40 -0.90      % refine (+ relocation)
        -1.00 -0.80  0.10  0.80      % coarsen
        -2.00 -2.00  1.00  0.60 ];   % accept

% ---------------- step-size / policy hyper-parameters -------------------
P.qTarget     = 0.35;    % aim at 35 % of the tolerance after a refinement step
P.gmin        = 1.15;    % minimum growth factor of the number of intervals when refining
P.gmaxPre     = 2.0;     % maximum growth per step while the estimator is unreliable
P.gmaxAsy     = 4.0;     % ... and when the asymptotic regime is detected
P.stagnRatio  = 0.5;     % a relocation pass 'pays off' if it reduces q by at least this factor
P.relocThresh = 0.20;    % relocate points together with a refinement if P(Misplaced) > this
P.gcoarseMin  = 0.5;     % never remove more than half of the intervals in one coarsening step
P.qFloor      = 1e-4;    % error estimates below 1e-4*TOL are round-off dominated, do not extrapolate
P.Nmin        = 8;       % never coarsen below this many intervals
end
