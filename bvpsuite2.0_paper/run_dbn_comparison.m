function results = run_dbn_comparison(which_problems, opts)
% RUN_DBN_COMPARISON  Classic bvpsuite2.0 mesh adaptation versus the Dynamic-Bayesian-Network
% controller, on the four problems of Hohenegger/Settanni/Weinmueller/Wolde (ANM 205, 2024).
%
% Usage (MATLAB or GNU Octave; cd into this folder first):
%     run_dbn_comparison                 % problems 1-4, all scenarios
%     run_dbn_comparison([1 4])          % only Problems 1 and 4
%     run_dbn_comparison(1:4, struct('orders',6))      % quicker: collocation order 6 only
%
% Options (struct): orders (default [6 8]), plots (default true), show (MATLAB: true)
%
% Every run uses the SAME solver, the SAME h-h/2 error estimate and the SAME acceptance
% test (scaled error estimate < 1, abs. tol = rel. tol. = TOL); only the strategy that
% decides "relocate / refine / coarsen / accept" differs.  The DBN therefore changes
% the effort, not the accuracy requirement.
%
% Effort measures printed per run
%     solves  number of collocation solves + error estimates (= mesh-adaptation iterations)
%     work    sum over the iterations of (coarse + fine) mesh intervals = 3*N
%     N       number of mesh intervals of the final mesh
%     time    wall-clock seconds
%
% Scenarios
%   A  coarse start mesh (11 points), TOL = 1e-8  -- the user does not know how many points are needed
%   B  paper start meshes (51 points; 21 for Problem 1 as in Table 1), TOL = 1e-8
%   C  tight tolerance TOL = 1e-11 from a 21-point start (needs much more adaptation)

if nargin < 1 || isempty(which_problems), which_problems = 1:4; end
if nargin < 2, opts = struct(); end
if ~isfield(opts,'orders'), opts.orders = [6 8]; end
if ~isfield(opts,'plots'),  opts.plots = true;   end
if ~isfield(opts,'show'),   opts.show = ~exist('OCTAVE_VERSION','builtin'); end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'bvpsuite2.0'), fullfile(here,'PaperProblems'));
warning('off','all');
global PAPER PAPERLOG BVPSTATS DBNTRACE
PAPERLOG = '';
figdir = fullfile(here,'figures_dbn');
if opts.plots && ~exist(figdir,'dir'), mkdir(figdir); end

emit('\nClassic bvpsuite2.0 mesh adaptation  vs.  Dynamic Bayesian Network controller\n');
emit('Platform: %s\n', version_string());
emit('Same solver, same error estimate, same acceptance test; only the adaptation strategy differs.\n');

% ---- list of cases -----------------------------------------------------------------
% each case: struct(prob, order, start, TOL, label, fields...)
cases = {};
for o = opts.orders
    for pr = which_problems
        % scenario A: coarse start
        cases = [cases, mk(pr, o, 11, 1e-8, 'A')];                  %#ok<AGROW>
    end
end
for o = opts.orders
    for pr = which_problems
        % scenario B: paper start meshes
        if pr == 1, st = [51 21]; else, st = 51; end
        for s = st, cases = [cases, mk(pr, o, s, 1e-8, 'B')]; end   %#ok<AGROW>
    end
end
for pr = which_problems         % scenario C: tight tolerance, hard variants only
    switch pr
        case 1, for o = opts.orders, cases = [cases, mk(1, o, 21, 1e-11, 'C')]; end
        case 3, cases = [cases, mk(3, opts.orders(1), 21, 1e-11, 'C', 'hard')];
        case 4, cases = [cases, mk(4, opts.orders(1), 21, 1e-11, 'C', 'hard')];
    end
end

% ---- run ----------------------------------------------------------------------------
R = [];  tstart = tic;  curScen = '';
titles = struct('A','A: coarse start mesh (11 points), TOL = 1e-8', ...
                'B','B: starting meshes of the paper, TOL = 1e-8', ...
                'C','C: tight tolerance TOL = 1e-11, 21-point start');
for ic = 1:numel(cases)
    c = cases{ic};
    if ~strcmp(c.scen, curScen)
        curScen = c.scen;
        emit('\n==== Scenario %s ====\n', titles.(curScen));
        emit('%-26s %-3s | %6s %7s %5s %7s | %6s %7s %5s %7s | %s\n','case','ord','solves','work','N','time','solves','work','N','time','');
        emit('%-26s %-3s | %-29s | %-29s | %s\n','','','        classic','        DBN','solves  work  time');
    end
    for im = 1:2
        meth = {'classic','dbn'};
        setup(c, meth{im});
        BVPSTATS = struct('iters',0,'work',0,'Nhist',[],'errhist',[]);
        DBNTRACE = struct('rows',zeros(0,9));
        t0 = tic;  ok = true;
        try
            evalc('[x,~,s] = bvpsuite2(c.name,''bvps2_paper_settings'');');
            err = max(abs(s.errest(:)));
        catch ME
            ok = false;  x = [];  err = NaN;  lastmsg = ME.message;
        end
        r.msg = ''; if ~ok, r.msg = lastmsg; end
        r.ok = ok; r.iters = BVPSTATS.iters; r.work = BVPSTATS.work; r.N = numel(x)-1;
        r.time = toc(t0); r.err = err; r.Nhist = BVPSTATS.Nhist; r.errhist = BVPSTATS.errhist;
        r.trace = DBNTRACE.rows;
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
    emit('Classic adaptation failed (iteration limit reached) in %d case(s) where the DBN succeeded: marked "fail".\n', nfail);
end
emit('Total run time: %.1f s\n', toc(tstart));
results = R;

% ---- figures -----------------------------------------------------------------------
if opts.plots
    save_fig(figdir, opts, 'dbn_fig1_solves_classic_vs_dbn', @() plot_bars(R));
    k = find(arrayfun(@(r) r.case.scen == 'A' && r.case.prob == 1 && r.case.order == opts.orders(1), R), 1);
    if ~isempty(k)
        save_fig(figdir, opts, 'dbn_fig2_problem1_error_history', @() plot_hist(R(k)));
        save_fig(figdir, opts, 'dbn_fig3_problem1_beliefs',        @() plot_belief(R(k)));
    end
    k = find(arrayfun(@(r) r.case.scen == 'C' && r.case.prob == 1, R), 1);
    if ~isempty(k), save_fig(figdir, opts, 'dbn_fig4_problem1_tight_tol', @() plot_hist(R(k))); end
end
try
    fid = fopen(fullfile(here,'dbn_comparison_output.txt'),'w');  fprintf(fid,'%s',PAPERLOG);  fclose(fid);
    emit('Output also written to %s\n', fullfile(here,'dbn_comparison_output.txt'));
catch
end
end

% ==========================================================================================
function c = mk(pr, order, start, TOL, scen, variant)
% one or two cases per problem (the easy variant and the hard one)
if nargin < 6, variant = 'all'; end
names = {'bvps2_paper_p1','bvps2_paper_p2','bvps2_paper_p3','bvps2_paper_p4'};
base = struct('prob',pr,'name',names{pr},'order',order,'start',start,'TOL',TOL,'scen',scen, ...
              'label','','p2_init',0,'p3_L',20,'p3_bc',1,'p4_xi',0.6,'solTol',1e-12);
if pr == 4, base.solTol = 1e-10; end
switch pr
    case 1
        b = base;  b.label = sprintf('P1 start %d', start);  c = {b};
    case 2
        b1 = base; b1.p2_init = 0; b1.label = 'P2 y1';  b2 = base; b2.p2_init = 1; b2.label = 'P2 y3';  c = {b1, b2};
    case 3
        b1 = base; b1.p3_L = 5;  b1.p3_bc = 1; b1.label = 'P3 L=5, y+y''=0';
        b2 = base; b2.p3_L = 20; b2.p3_bc = 2; b2.label = 'P3 L=20, y=0';
        c = {b1, b2};  if strcmp(variant,'hard'), c = {b2}; end
    case 4
        b1 = base; b1.p4_xi = 0.1; b1.label = 'P4 xi=0.1';  b2 = base; b2.p4_xi = 0.6; b2.label = 'P4 xi=0.6';
        c = {b1, b2};  if strcmp(variant,'hard'), c = {b2}; end
end
if pr == 1 && strcmp(scen,'B'), c{1}.label = sprintf('P1 start %d', start); end
end

function setup(c, meth)
global PAPER
PAPER = struct('TOL',c.TOL,'Nmesh',c.start,'order',c.order,'solTol',c.solTol,'mode',meth,'coarsen',0, ...
               'p2_init',c.p2_init,'p3_L',c.p3_L,'p3_bc',c.p3_bc,'p4_xi',c.p4_xi);
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
    emit('  figure saved: %s\n', fullfile('figures_dbn',[name '.png']));
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
