function results = run_dbn_comparison(which_problems, opts)
% RUN_DBN_COMPARISON  OPTIONAL add-on: classic bvpsuite2.0 mesh adaptation versus a
% Dynamic-Bayesian-Network (DBN) mesh controller, on the four problems.
%
% FILES (same folder):
%   problems.m   the four problems and all their settings (start meshes, orders, TOL, ...)  <- used here
%   bvp_main.m   the bvpsuite2.0 solver (all its files in one file)                          <- used here
%   this file    contains the DBN mesh controller (Section 2) and the comparison driver.
%   problems.m, bvp_main.m and run_paper_results.m never use this file.
%   The DBN is plugged into the solver through its generic setting 'meshAdaptationFcn'; the solver
%   routines it needs are handed over by bvp_main itself (bvp_main('api')).
%
% Usage (MATLAB or GNU Octave):
%     run_dbn_comparison                 % problems 1-4, all scenarios
%     run_dbn_comparison([1 4])          % only Problems 1 and 4
%     run_dbn_comparison(1:4, struct('orders',6))      % quicker
% Options (struct): orders, plots, show -- defaults are in section 1 below.
%
% Both strategies use the SAME solver, error estimate and acceptance test (scaled error
% estimate < 1); only the rule that chooses "relocate / refine / coarsen / accept" differs.
% Effort measures:  solves (= adaptation iterations), work (sum of 3N), N (final intervals), time.

%% ======================= SECTION 1.  SETTINGS (edit here) =======================
D.orders      = [6 8];        % collocation orders compared
D.plots       = true;         % draw / save the comparison figures
D.show        = ~exist('OCTAVE_VERSION','builtin');   % show them on screen
D.figDir      = 'figures_dbn';
D.logFile     = 'dbn_comparison_output.txt';
D.coarsen     = 0;            % 1 = let the DBN also remove points
% scenarios:  A coarse start mesh | B the start meshes set in problems.m | C tight TOL
D.A.start = 11;                      % scenario A: start mesh points (TOL of the problem)
D.C.start = 21;   D.C.TOL = 1e-11;   % scenario C: start mesh points and tolerance
%% ================================================================================

if nargin < 1 || isempty(which_problems), which_problems = 1:4; end
if nargin < 2, opts = struct(); end
fn = fieldnames(opts);  for i = 1:numel(fn), D.(fn{i}) = opts.(fn{i}); end

here = fileparts(mfilename('fullpath'));
warning('off','all');
global PAPERLOG BVPSTATS DBNTRACE
PAPERLOG = '';
cfg = problems('config');                     % problems + settings (problems.m)
PC  = cfg.P;
figdir = fullfile(here, D.figDir);
if D.plots && ~exist(figdir,'dir'), mkdir(figdir); end
opts = D;  opts.show = D.show;

emit('\nClassic bvpsuite2.0 mesh adaptation  vs.  Dynamic Bayesian Network controller\n');
emit('Platform: %s\n', version_string());
emit('Same solver, same error estimate, same acceptance test; only the adaptation strategy differs.\n');

% ---- list of cases -------------------------------------------------------------------
cases = {};
for o = D.orders
    for pr = which_problems
        cases = [cases, mk(pr, o, D.A.start, PC{pr}.TOL, 'A', 'all')];                %#ok<AGROW>
    end
end
for o = D.orders
    for pr = which_problems
        for s = PC{pr}.startMesh
            cases = [cases, mk(pr, o, s, PC{pr}.TOL, 'B', 'all')];                    %#ok<AGROW>
        end
    end
end
for pr = which_problems       % scenario C: hard variants only (last variant of the problem)
    if pr == 1
        for o = D.orders, cases = [cases, mk(1, o, D.C.start, D.C.TOL, 'C', 'all')]; end    %#ok<AGROW>
    elseif pr >= 3
        cases = [cases, mk(pr, D.orders(1), D.C.start, D.C.TOL, 'C', 'last')];              %#ok<AGROW>
    end
end

% ---- run ----------------------------------------------------------------------------
R = [];  tstart = tic;  curScen = '';
titles = struct('A',sprintf('A: coarse start mesh (%d points), TOL of the problem',D.A.start), ...
                'B','B: start meshes set in problems.m, TOL of the problem', ...
                'C',sprintf('C: tight tolerance TOL = %.0e, %d-point start',D.C.TOL,D.C.start));
for ic = 1:numel(cases)
    c = cases{ic};
    if ~strcmp(c.scen, curScen)
        curScen = c.scen;
        emit('\n==== Scenario %s ====\n', titles.(curScen));
        emit('%-26s %-3s | %6s %7s %5s %7s | %6s %7s %5s %7s | %s\n','case','ord','solves','work','N','time','solves','work','N','time','');
        emit('%-26s %-3s | %-29s | %-29s | %s\n','','','        classic','        DBN','solves  work  time');
    end
    for im = 1:2
        if im == 1, adapt = '';                    extra = struct();
        else,       adapt = @meshadaptation_dbn; extra = struct('dbnCoarsen',D.coarsen,'dbnVerbose',0);
        end
        B = problems('build', c.k, c.variant, c.start, c.order, c.TOL, adapt, extra);
        problem = B.problem;  settings = B.settings;
        BVPSTATS = struct('iters',0,'work',0,'Nhist',[],'errhist',[]);
        DBNTRACE = struct('rows',zeros(0,9));
        t0 = tic;  ok = true;  lastmsg = '';
        try
            evalc('[x,~,s] = bvp_main(problem, settings);');
            err = max(abs(s.errest(:)));
        catch ME
            ok = false;  x = [];  err = NaN;  lastmsg = ME.message;
        end
        r = struct('msg',lastmsg,'ok',ok,'iters',BVPSTATS.iters,'work',BVPSTATS.work,'N',numel(x)-1, ...
                   'time',toc(t0),'err',err,'Nhist',BVPSTATS.Nhist,'errhist',BVPSTATS.errhist,'trace',DBNTRACE.rows);
        res{im} = r;                                               %#ok<AGROW>
    end
    a = res{1};  b = res{2};
    if ~b.ok, emit('  [DBN run failed: %s]\n', strtok(b.msg, char(10))); end
    emit('%-26s %-3d | %6s %7d %5d %6.1fs | %6s %7d %5d %6.1fs | %s\n', c.label, c.order, ...
         itxt(a), a.work, a.N, a.time, itxt(b), b.work, b.N, b.time, ratio(a,b));
    R = [R; struct('case',c,'classic',a,'dbn',b)];                 %#ok<AGROW>
end

% ---- summary ------------------------------------------------------------------------
emit('\n==== Summary (cases where both strategies reached the tolerance) ====\n');
for sc = 'ABC'
    idx = find(arrayfun(@(r) r.case.scen == sc && r.classic.ok && r.dbn.ok, R));
    if isempty(idx), continue; end
    si = sum(arrayfun(@(r) r.classic.iters, R(idx)));  di = sum(arrayfun(@(r) r.dbn.iters, R(idx)));
    sw = sum(arrayfun(@(r) r.classic.work,  R(idx)));  dw = sum(arrayfun(@(r) r.dbn.work,  R(idx)));
    st = sum(arrayfun(@(r) r.classic.time,  R(idx)));  dt = sum(arrayfun(@(r) r.dbn.time,  R(idx)));
    emit('Scenario %s: solves %3d -> %3d (%.0f%% fewer) | work %6d -> %6d (%.0f%% less) | time %6.1fs -> %6.1fs (%.0f%% less)\n', ...
         sc, si, di, 100*(1-di/si), sw, dw, 100*(1-dw/sw), st, dt, 100*(1-dt/st));
end
nfail = sum(arrayfun(@(r) ~r.classic.ok, R));
if nfail > 0
    emit('Classic adaptation failed (iteration limit reached) in %d case(s): marked "fail".\n', nfail);
end
emit('Total run time: %.1f s\n', toc(tstart));
results = R;

% ---- figures -----------------------------------------------------------------------
if D.plots
    save_fig(figdir, opts, 'dbn_fig1_solves_classic_vs_dbn', @() plot_bars(R));
    k = find(arrayfun(@(r) r.case.scen == 'A' && r.case.k == 1 && r.case.order == D.orders(1), R), 1);
    if ~isempty(k)
        save_fig(figdir, opts, 'dbn_fig2_problem1_error_history', @() plot_hist(R(k)));
        save_fig(figdir, opts, 'dbn_fig3_problem1_beliefs',        @() plot_belief(R(k)));
    end
    k = find(arrayfun(@(r) r.case.scen == 'C' && r.case.k == 1, R), 1);
    if ~isempty(k), save_fig(figdir, opts, 'dbn_fig4_problem1_tight_tol', @() plot_hist(R(k))); end
end
try
    fid = fopen(fullfile(here,D.logFile),'w');  fprintf(fid,'%s',PAPERLOG);  fclose(fid);
    emit('Output also written to %s\n', fullfile(here,D.logFile));
catch
end
end

% ==========================================================================================
function c = mk(k, order, start, TOL, scen, which)
% cases of problem k: all variants of the main file (or only the last, "hard", one)
V = problems('variants', k);
if strcmp(which,'last'), V = V(end); end
c = cell(1,numel(V));
for i = 1:numel(V)
    lab = sprintf('P%d %s', k, V(i).label);
    if k == 1, lab = sprintf('P1 start %d', start); end
    if k == 3, lab = sprintf('P3 bc%d %s', V(i).par.bc, V(i).label); end
    c{i} = struct('k',k,'prob',k,'order',order,'start',start,'TOL',TOL,'scen',scen,'label',lab,'variant',V(i));
end
end

function s = itxt(r)
if r.ok, s = sprintf('%d', r.iters); else, s = 'fail'; end
end

function t = ratio(a, b)
if ~a.ok && b.ok, t = 'classic failed';  return; end
if a.iters == b.iters, t = 'same solves'; else, t = sprintf('%d -> %d solves', a.iters, b.iters); end
end

function emit(fmt, varargin)
global PAPERLOG
str = sprintf(fmt, varargin{:});
PAPERLOG = [PAPERLOG str];
fprintf('%s', str);
end

function v = version_string()
if exist('OCTAVE_VERSION','builtin'), v = ['GNU Octave ' OCTAVE_VERSION]; else, v = ['MATLAB ' version]; end
end

function save_fig(figdir, opts, name, plotfun)
try
    if opts.show, vis = 'on'; else, vis = 'off'; end
    f = figure('Visible',vis,'Color','w','Name',name,'NumberTitle','off','Position',[100 100 900 560]);
    plotfun();  drawnow;
    print(f, fullfile(figdir,[name '.png']), '-dpng', '-r150');
    if ~opts.show, close(f); end
    emit('  figure saved: %s\n', fullfile(opts.figDir,[name '.png']));
catch err
    emit('  (could not draw/save %s: %s)\n', name, err.message);
end
end

% ---------------------------------------------------------------- figures
function plot_bars(R)
% number of solves per case, classic vs DBN (only cases where some adaptation happens)
keep = find(arrayfun(@(r) r.classic.iters > 1 || ~r.classic.ok || r.dbn.iters > 1, R));
n = numel(keep);  a = zeros(n,1);  b = zeros(n,1);  lab = cell(n,1);
for k = 1:n
    r = R(keep(k));
    a(k) = r.classic.iters;  b(k) = r.dbn.iters;
    if ~r.classic.ok, a(k) = 30; end                      % iteration limit
    lab{k} = sprintf('%s  o%d  [%s]', r.case.label, r.case.order, r.case.scen);
end
h = barh(flipud([a b]), 1);
try, set(h(1),'FaceColor',[0.85 0.33 0.10]); set(h(2),'FaceColor',[0 0.45 0.74]); catch, end
set(gca,'YTick',1:n,'YTickLabel',flipud(lab),'FontSize',9);
legend({'classic bvpsuite2.0','DBN controller'},'Location','southeast');
xlabel('collocation solves until the tolerance is met  (30 = classic hit its iteration limit)');
title({'Cases that need mesh adaptation', '[A] 11-point start, TOL 1e-8   [C] 21-point start, TOL 1e-11   (o = collocation order)'});
box on; grid on;
end

function plot_hist(r)
% error estimate and mesh size per iteration for one case
c = r.case;  hold on;
semilogy(1:numel(r.classic.errhist), r.classic.errhist, 'o-', 'Color',[0.85 0.33 0.1], 'LineWidth',1.8, 'MarkerSize',6);
semilogy(1:numel(r.dbn.errhist),     r.dbn.errhist,     's-', 'Color',[0 0.45 0.74], 'LineWidth',1.8, 'MarkerSize',6);
set(gca,'YScale','log');
plot([1 max(numel(r.classic.errhist), numel(r.dbn.errhist))], [c.TOL c.TOL], 'k--');
hold off; grid on; box on; xlabel('mesh-adaptation iteration (one solve each)'); ylabel('max error estimate');
legend({sprintf('classic (%d solves, final N = %d)', r.classic.iters, r.classic.N), ...
        sprintf('DBN (%d solves, final N = %d)', r.dbn.iters, r.dbn.N), 'TOL'}, 'Location','northeast');
title(sprintf('%s, order %d, TOL = %.0e', c.label, c.order, c.TOL));
end

function plot_belief(r)
% DBN belief over the mesh state S and chosen action, per iteration
T = r.dbn.trace;  if isempty(T), return; end
nit = size(T,1);
bS = T(:,4:7);  act = T(:,9);
h = bar(bS, 'stacked'); hold on;
cols = [0.80 0.30 0.20; 0.95 0.65 0.15; 0.30 0.70 0.40; 0.30 0.50 0.85];
for kk = 1:4, try, set(h(kk),'FaceColor',cols(kk,:)); catch, end, end
an = {'relocate','refine','coarsen','accept'};
for i = 1:nit
    text(i, 1.04, an{act(i)}, 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
    text(i, -0.06, sprintf('M=%d', T(i,2)), 'HorizontalAlignment','center','FontSize',8);
end
hold off; ylim([-0.1 1.12]); xlim([0.4 nit+0.6]);
legend({'Misplaced','Underresolved','Adequate','Overresolved'}, 'Location','eastoutside');
xlabel('iteration'); ylabel('belief P(S_t | observations)');
title(sprintf('DBN beliefs and decisions: %s, order %d, TOL = %.0e', r.case.label, r.case.order, r.case.TOL));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SECTION 2.   DYNAMIC BAYESIAN NETWORK MESH CONTROLLER  (formerly dbn_controller.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SECTION 3.   MESH ADAPTATION STEERED BY THE DBN  (formerly meshadaptation_dbn.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [finalsol,toleriert,halve] = meshadaptation_dbn(problem,settings,x1,rhoColl,initProfile,SV)
% MESHADAPTATION_DBN  Mesh adaptation for bvpsuite2.0 steered by a Dynamic Bayesian Network.
%
% (SV = struct of handles to the solver routines, passed in by bvp_main.)
%
% Same ingredients as meshadaptation.m (collocation solve, h-h/2 error estimate,
% residual-based density update), same acceptance test (scaled error estimate < 1),
% but the choice of the next action  {relocate, refine(+relocate), coarsen, accept}
% is made by dbn_controller from the observations of the current and previous
% iterations (see dbn_controller.m and README).  Differences to the classic loop:
%   * relocation and refinement can be combined in one step,
%   * the new number of intervals is extrapolated with the *observed* convergence order,
%   * relocation passes stop as soon as they stop paying off,
%   * (optional) coarsening if the error is far below the tolerance; a failed
%     coarsening step is rolled back to the last verified mesh, so the returned
%     solution always satisfies the same tolerance test as in the classic code.

global DBNTRACE
aTOL = feval(settings,'absTolMeshAdaptation');
rTOL = feval(settings,'relTolMeshAdaptation');
K    = feval(settings,'K');
maxiter = 30;
try, allowCoarsen = feval(settings,'dbnCoarsen'); catch, allowCoarsen = 0; end
try, verbose = feval(settings,'dbnVerbose');      catch, verbose = 0; end
halve = 0;  toleriert = 0;

ordnung = SV.feval_problem(problem,'orders');
n   = length(ordnung);
linear = SV.feval_problem(problem,'linear');
m   = length(rhoColl);
if ~exist('initProfile','var'), initProfile = []; end
if length(aTOL) <= n, aTOL = ones(1,n)*aTOL(1); end
if length(rTOL) <= n, rTOL = ones(1,n)*rTOL(1); end
aTOL = aTOL(:); rTOL = rTOL(:);

left_end = x1(1);  right_end = x1(end);
rho = (1./diff(x1)/length(x1))';
LPFilter_k = m*2.75;

dbn = dbn_controller('init', length(x1)-1, m);
best = [];  probing = false;  nCoarsen = 0;  prevMesh = [];
x1_try = x1;

for j = 1:maxiter
    Mi = length(x1_try)-1;

    % ------------------------------------------------------------------ solve
    solved = true;
    try
        if linear
            [tcol,ycol,sol] = SV.solveLinearProblem(problem,x1_try,rhoColl);
        else
            [tcol,ycol,sol,halve] = SV.solveNonLinearProblem(problem,settings,x1_try,rhoColl,initProfile,1-(j>1),[]);
            if halve ~= 0, solved = false; end
        end
        if solved
            [sol_c,~,halve] = SV.errorestimate(sol,problem,settings,[]);
            if halve ~= 0, solved = false; end
        end
    catch err
        if ~probing, rethrow(err); end
        solved = false;
    end

    if ~solved
        if probing && ~isempty(best)          % failed coarsening probe: roll back
            if verbose, fprintf('  [DBN] probe on %d intervals failed -> keep %d intervals\n', Mi, length(best.x1)-1); end
            halve = 0;  finalsol = best;  return
        end
        toleriert = 0; finalsol = 0; return
    end

    % ---------------------------------------------------------- observations
    ecol = abs(sol_c.errest);
    ysc  = aTOL + abs(sol.valx1tau).*rTOL;           % n x npts (aTOL/rTOL are n x 1)
    escl = ecol ./ ysc;
    q = max(escl(:));
    colmax = max(escl,[],1);
    E = max(reshape(colmax(1:Mi*(m+1)), m+1, Mi), [], 1);
    E(Mi) = max(E(Mi), colmax(end));
    skewE = max(E)/max(mean(E), realmin);

    % defect (residual) per interval, as in the classic code ('tcolmid')
    tau = sol.x1;  par = sol.parameters;
    resmesh = tcol(1:end-1) + diff(tcol)/2;
    canCoarsen = logical(allowCoarsen) && ~probing && nCoarsen < 2 && Mi > 12;
    residualfinal = E;  skewD = skewE;
    if q >= 1 || canCoarsen        % the defect is only needed to relocate points / judge coarsening
    try
        residual = abs(SV.computeResidual(problem,sol.coeff,par,tau,rhoColl,resmesh,0));
        residualfinal = zeros(1,Mi);
        for i = 1:size(residual,1)
            radd = residual(i,:).*diff(tcol);
            radd2 = zeros(1,Mi);
            for k = 1:Mi
                radd2(k) = sum(radd(k*(m+1)-m:k*(m+1)))/(x1_try(k+1)-x1_try(k));
            end
            residualfinal = residualfinal + radd2.^2;
        end
        residualfinal = sqrt(residualfinal);
        skewD = max(residualfinal)/max(mean(residualfinal), realmin);
    catch
        residualfinal = E;  skewD = skewE;
    end
    end

    obs = struct('q',q,'skewE',skewE,'skewD',skewD,'N',Mi,'canCoarsen',canCoarsen);
    S = dbn_controller('step', dbn, obs);
    dbn = S.dbn;  dec = S.dec;

    if isstruct(DBNTRACE)
        DBNTRACE.rows(end+1,:) = [j Mi q dec.bS dec.bR(2) dec.action];
    end
    if verbose
        an = {'relocate','refine','coarsen','accept'};
        fprintf('  [DBN] it %2d  M=%3d  q=%9.2e skewE=%5.1f skewD=%5.1f | P(S)=[%.2f %.2f %.2f %.2f] P(asym)=%.2f p_obs=%5.2f -> %s (M_new=%d)\n', ...
            j, Mi, q, skewE, skewD, dec.bS, dec.bR(2), dec.pobs, an{dec.action}, dec.Mnew);
    end

    % remember the solution of this mesh
    cur = struct('x1',sol.x1,'valx1',sol.valx1,'x1tau',sol.x1tau,'valx1tau',sol.valx1tau, ...
                 'coeff',sol.coeff,'parameters',sol.parameters,'errest',sol_c.errest,'q',q);
    if q < 1
        if probing
            probing = false;
        end
        best = cur;
    elseif probing                         % coarsening failed the tolerance test -> roll back
        if verbose, fprintf('  [DBN] coarsened mesh (%d intervals) violates TOL -> keep %d intervals\n', Mi, length(best.x1)-1); end
        finalsol = best;  return
    end

    % ------------------------------------------------------------------ act
    if dec.action == 4                      % accept
        finalsol = best;  return
    end

    if ~linear
        initProfile.initialMesh   = sol.x1;
        initProfile.initialValues = sol.valx1;
        initProfile.initialCoeff  = sol.coeff;
        initProfile.parameters    = sol.parameters;
    end

    Mnew = dec.Mnew;
    if dec.action == 3
        probing = true;  nCoarsen = nCoarsen + 1;
    end

    % new mesh: update the density only for relocation / refine+relocate
    if dec.relocate
        rho = LPfilter(residualfinal', (1./diff(x1_try)/length(x1_try))', LPFilter_k);
    else
        rho = (1./diff(x1_try)/length(x1_try))';
    end
    rho = max(rho, max(rho)/K);
    rho = bcTDFlogV4(rho);
    xn = (x1_try - x1_try(1))/(x1_try(end)-x1_try(1));
    I  = cumtrapz(xn, ([rho(1); rho]+[rho; rho(end)])'/2);
    xn = pchip(I/I(end), xn, 0:1/Mnew:1);
    x1_try = left_end + xn*(right_end-left_end);
end

% iteration limit: return the best verified solution (or the last one)
if ~isempty(best), finalsol = best; else, finalsol = cur; end
toleriert = 0;
end

% ----------------------------------------------------------------------------
function newrho = LPfilter(err,rho,LPFilter_k)
% Generate step density profile update
% For use with 2pBVP solver for 1st order systems
% Euler-Lagrange optimal grid is generated using
% using deadbeat control law. Local error is
% equidistributed over the interval.

% err is already on staggered grid in this version
k = LPFilter_k;       % k = 3 for DBC, k = 4 conv filter (root 1/3)
M = length(rho);

% Process input error
scalederr = M*abs(err).^(1/k);
% scalederr = processV4(bcTDFlogV4(scalederr)); % SWITCHED OFF

ratio = scalederr;
% ratio is suggested update of rho; compute next rho
rhonew = rho.*ratio;
%rhonew = bcTDFlogV4(rhonew);                   % SWITCHED OFF
newrho = rhonew;%*sum(1./rhonew)/M;
end


function out = bcTDFlogV4(in)
% Boundary corrected Topelitz digital filter, multiplicative (logarithmic)
% Applies 2nd order convolution LP filter repeatedly
% to input signal until output signal is sufficiently smooth

smooth = 50;     % Define smoothness requirement
N = length(in);
old = in;
signal = old;    % Mem alloc
snratio = 0;

while snratio < smooth
    % Apply filter until S/N ratio at this stage is at least = smooth
    for i=2:N-1
        signal(i) = (old(i-1)*old(i)^2*old(i+1))^(1/4);
    end
    % Boundary corrections
    signal(1) = (old(1)^3*old(2)^3/old(3)^2/old(4)*old(5))^(1/4);
    % signal(1) = (old(1)^2*old(2)^3/old(3))^(1/4);
    signal(N) = (old(N)^3*old(N-1)^3/old(N-2)^2/old(N-3)*old(N-4))^(1/4);
    % signal(N) = (old(N)^2*old(N-1)^3/old(N-2))^(1/4);
    
    % Compute noise
    noise = old - signal;
    s = norm(signal);
    n = norm(noise);
    snratio = s/(n + 1e-3*s);
    old = signal;
end

out = signal;
end
