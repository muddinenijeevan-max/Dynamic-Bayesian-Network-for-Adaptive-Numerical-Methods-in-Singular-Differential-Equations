function results = run_paper_results(which_problems, opts)
% RUN_PAPER_RESULTS  Runs the four test problems of
%
%   M. Hohenegger, G. Settanni, E.B. Weinmueller, M. Wolde,
%   "Numerical treatment of singular ODEs using finite difference and
%    collocation methods", Appl. Numer. Math. 205 (2024) 184-194,
%
%   with bvpsuite2.0 (polynomial collocation, Gauss points, mesh adaptation with the
%   h-h/2 error estimate) -- the runs behind Tables 1, 4, 6 and 8 and Figs. 2-5, 7.
%
% HOW TO USE
%   run_paper_results                  % every problem whose  enable = true  in problems.m
%   run_paper_results(1)               % Problem 1 only        (ignores the enable flags)
%   run_paper_results([2 4])           % Problems 2 and 4
%   run_paper_results(1:4, struct('verbose',true))      % also show the solver's own output
%
%   Everything you may want to change -- which problems run, their initial mesh points,
%   collocation orders, tolerance, problem parameters, and whether their figures are
%   drawn / shown / saved -- is in SECTION 1 of problems.m, one block per problem.
%
% WHAT IS PRINTED
%   For every problem a separate block with one row per run (case x start mesh x order).
%   Only values computed in this run are printed (nothing is compared with the paper).
%       start       number of points of the equidistant START mesh
%       order       number m of Gauss collocation points (the "order" of the paper's tables)
%       mesh        number of points of the FINAL (adapted) mesh
%       error est.  max of bvpsuite2.0's own h-h/2 error estimate on the final mesh
%                   (the "error" / "abs. error" column of the paper's tables)
%       exact abs.  Problem 1 only: max |y_numerical - y_exact| over the mesh points
%       exact rel.  Problem 1 only: exact abs. / max|y_exact|
%
% FILES  (all in the same folder)
%   run_paper_results.m   this file: runs the problems, prints the tables, draws the figures
%   problems.m            the four problems and ALL their settings (start meshes, orders, TOL, plots, ...)
%   bvp_main.m            the bvpsuite2.0 collocation solver (all bvpsuite2.0 files in one file)
%   run_dbn_comparison.m  OPTIONAL: Dynamic-Bayesian-Network mesh controller; uses problems.m and
%                         bvp_main.m.  None of the other three files uses it.
%
% Works in MATLAB and GNU Octave.

cfg   = problems('config');          % problems and all their settings (problems.m)
PROBS = cfg.P;
G     = cfg.G;

if nargin < 1 || isempty(which_problems)
    which_problems = find(cellfun(@(P) P.enable, PROBS));
end
if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts,'verbose'), opts.verbose = false; end

here = fileparts(mfilename('fullpath'));
warning('off','all');          % Octave: harmless "colon arguments" warnings inside the solver
emit('@reset');  coef_text('@reset');
results = struct();
tstart = tic;

emit('\n');
emit('%s\n', repmat('#',1,96));
emit('#  bvpsuite2.0  --  singular BVPs of Hohenegger, Settanni, Weinmueller, Wolde (ANM 205, 2024)\n');
emit('#  Platform: %s\n', version_string());
emit('#  Gauss collocation, mixed error control (abs. tol. = rel. tol. = TOL), mesh adaptation, h-h/2 error estimate\n');
emit('%s\n', repmat('#',1,96));

runners = {@run_problem, @run_problem, @run_problem, @run_problem};
for k = which_problems(:).'
    P = PROBS{k};
    results.(sprintf('problem%d',k)) = runners{k}(k, P, G, here, opts);
end

emit('\nTotal run time: %.1f s\n', toc(tstart));
try
    fid = fopen(fullfile(here,G.logFile),'w');
    fprintf(fid,'%s',emit('@get'));  fclose(fid);
    emit('Output also written to %s\n', fullfile(here,G.logFile));
catch
end
try
    ctext = coef_text('@get');
    if ~isempty(ctext)
        fid = fopen(fullfile(here,G.coeffFile),'w');
        fprintf(fid,'%s',ctext);  fclose(fid);
        emit('Polynomial coefficients written to %s\n', fullfile(here,G.coeffFile));
    end
catch
end
end


%% ========================================================================================
%%  SECTION 2.   RUNNER  (one generic loop for all problems)
%% ========================================================================================

function R = run_problem(k, P, G, here, opts)
% Run every (case, start mesh, order) combination of problem k, print the block, draw figures.
V = problems('variants', k);
hasExact = (k == 1);

emit('\n%s\n', repmat('=',1,96));
emit(' PROBLEM %d    %s\n', k, P.title);
emit('%s\n', repmat('=',1,96));
emit(' start meshes: %s     orders: %s     TOL = %.0e     Newton tol = %.0e\n', ...
     mat2str(P.startMesh), mat2str(P.orders), P.TOL, P.solTol);
if k == 4, emit(' (solved on tau = 1 - 1/sqrt(1+t) in [0,1], as in the paper)\n'); end
if k == 2, emit(' (the paper could not obtain y2(t) with bvpsuite2.0 either)\n'); end

if hasExact
    emit('%s\n', repmat('-',1,96));
    emit(' %-14s %6s %6s %6s %6s %14s %14s %14s\n','case','start','order','TOL','mesh','error est.','exact rel.','exact abs.');
else
    emit('%s\n', repmat('-',1,96));
    emit(' %-14s %6s %6s %6s %6s %14s\n','case','start','order','TOL','mesh','error est.');
end

rows = struct('case',{},'start',{},'order',{},'mesh',{},'errest',{},'exactAbs',{},'exactRel',{});
coefRuns = struct('label',{},'start',{},'order',{},'iters',{});   % coefficients of every iteration
sols = cell(1,numel(V));          % solution of (first start mesh, first order) per case, for the figures
meshFinal = cell(1,numel(P.startMesh));   % final mesh of (start mesh, first order), for Fig. 3
lastGroup = '';

for iv = 1:numel(V)
    v = V(iv);
    if ~strcmp(v.group, lastGroup)
        emit('%s\n', repmat('-',1,96));
        if ~isempty(v.group), emit(' %s\n', v.group); end
        lastGroup = v.group;
    end
    for is = 1:numel(P.startMesh)
        for io = 1:numel(P.orders)
            extra = struct();
            if P.printCoeffs, extra.iterationHook = @coef_hook; coef_hook('@reset'); end
            B = problems('build', k, v, P.startMesh(is), P.orders(io), P.TOL, '', extra);
            problem = B.problem;  settings = B.settings;
            [x, s] = solve_quiet(problem, settings, opts.verbose);
            if P.printCoeffs
                coefRuns(end+1) = struct('label',v.label,'start',P.startMesh(is),'order',P.orders(io),'iters',coef_hook('@take')); %#ok<AGROW>
            end

            r.case = v.label;  r.start = P.startMesh(is);  r.order = P.orders(io);
            r.mesh = numel(x); r.errest = max(abs(s.errest(:)));
            r.exactAbs = NaN;  r.exactRel = NaN;
            if hasExact
                r.exactAbs = max(abs(s.valx1 - P.exact(s.x1)));    % max over the mesh points
                r.exactRel = r.exactAbs / P.exactMax;              % relative to max|y_exact|
            end
            rows(end+1) = r;                                        %#ok<AGROW>
            if is == 1 && io == 1, sols{iv} = s; end
            if io == 1 && iv == 1, meshFinal{is} = x; end

            if hasExact
                emit(' %-14s %6d %6d %6.0e %6d %14.4e %14.4e %14.4e\n', r.case, r.start, r.order, P.TOL, r.mesh, r.errest, r.exactRel, r.exactAbs);
            else
                emit(' %-14s %6d %6d %6.0e %6d %14.4e\n', r.case, r.start, r.order, P.TOL, r.mesh, r.errest);
            end
        end
        if numel(P.startMesh) > 1 && is < numel(P.startMesh), emit('\n'); end
    end
end
emit('%s\n', repmat('-',1,96));
R = rows;

% ---- polynomial coefficients of every iteration ---------------------------------------------------
if P.printCoeffs
    txt = format_coeffs(k, V, coefRuns);
    coef_text(txt);
    if G.coeffsOnScreen, emit('%s', txt); else, emit('\n (polynomial coefficients of problem %d: see %s)\n', k, G.coeffFile); end
end

% ---- figures -----------------------------------------------------------------------------------
if isfield(opts,'plots') && ~isempty(opts.plots), P.plots = opts.plots; end
if P.plots
    o1 = P.orders(1);
    switch k
        case 1
            save_fig(here,G,P,'fig2_problem1_exact', @() plot_exact1());
            for is = 1:numel(P.startMesh)
                save_fig(here,G,P,sprintf('fig3_problem1_meshes_start%d',P.startMesh(is)), ...
                         @() plot_meshes(P.startMesh(is), meshFinal{is}, o1));
            end
        case 2
            save_fig(here,G,P,'fig4_problem2_solutions', @() plot_fig4(sols, P.initNames, P.startMesh(1), o1));
        case 3
            letters = 'abcdefgh';
            for ib = 1:numel(P.bc)
                idx = find(arrayfun(@(q) q.par.bc == P.bc(ib), V));
                save_fig(here,G,P,sprintf('fig5%s_problem3_bc%d',letters(ib),P.bc(ib)), ...
                         @() plot_fig5(sols(idx), V(idx), P.bc(ib), P.startMesh(1), o1));
            end
        case 4
            save_fig(here,G,P,'fig7_problem4_solutions', @() plot_fig7(sols, P.xi, P.startMesh(1), o1));
    end
end
end

function [x, s] = solve_quiet(problem, settings, verbose)
% call bvp_main (bvpsuite2.0) and (unless verbose) swallow its console output
if verbose
    [x,~,s] = bvp_main(problem, settings);
else
    evalc('[x,~,s] = bvp_main(problem, settings);');
end
end

%% ========================================================================================
%%  SECTION 5.   FIGURES
%% ========================================================================================

function save_fig(here, G, P, name, plotfun)
% draw one figure, show it (if P.show) and save it as PNG (if P.save)
try
    if P.show, vis = 'on'; else, vis = 'off'; end
    f = figure('Visible',vis,'Color','w','Name',name,'NumberTitle','off','Position',[100 100 760 560]);
    plotfun();
    drawnow;
    if P.save
        if ~exist(fullfile(here,G.figDir),'dir'), mkdir(fullfile(here,G.figDir)); end
        print(f, fullfile(here,G.figDir,[name '.png']), '-dpng', '-r150');
        emit('  figure saved: %s/%s.png\n', G.figDir, name);
    end
    if ~P.show, close(f); end
catch err
    emit('  (could not draw/save %s: %s)\n', name, err.message);
end
end

function nice_axes()
box on; set(gca,'FontSize',11,'TickDir','in','LineWidth',0.8);
end

function c = palette(j)
cols = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556; 0.929 0.694 0.125];
c = cols(mod(j-1,size(cols,1))+1,:);
end

function plot_exact1()
% Fig. 2: exact solution y(t) = t^2 ln t
t = linspace(0,1,500);
plot(t, t.^2.*log(max(t,realmin)), 'Color',[0 0.447 0.741], 'LineWidth',3);
xlim([0 1]); ylim([-0.2 0]); set(gca,'XTick',0:0.1:1,'YTick',-0.2:0.02:0);
nice_axes(); xlabel('t'); title('Fig. 2.  Problem 1: exact solution y(t) = t^2 ln t');
end

function plot_meshes(Nstart, xfinal, order)
% Fig. 3: equidistant start mesh (Nstart points) and adapted final mesh (any number of points)
hold on;
plot(linspace(0,1,Nstart), ones(1,Nstart), 'o', 'Color',[0.85 0.2 0.1], 'LineWidth',2, 'MarkerSize',7);
plot(xfinal, zeros(1,numel(xfinal)), 'o', 'Color',[0 0 0.9], 'LineWidth',2, 'MarkerSize',7);
hold off;
xlim([-0.01 1.01]); ylim([-1.2 2.2]); set(gca,'XTick',0:0.1:1,'YTick',[]);
legend({sprintf('equidistant start mesh (%d points)',Nstart), ...
        sprintf('adapted final mesh (%d points)',numel(xfinal))},'Location','northeast');
nice_axes(); xlabel('t');
title(sprintf('Fig. 3.  Problem 1: start mesh and adapted mesh (bvpsuite2.0, order %d)', order));
end

function plot_fig4(sols, names, Nstart, order)
% Fig. 4: Problem 2, one curve per starting profile (y1: start y = 0, y3: start y = 1)
hold on;  mk = {'.','o','s','^','d'};  h = zeros(1,numel(sols));
for j = 1:numel(sols)
    h(j) = plot(sols{j}.x1tau, sols{j}.valx1tau(1,:), '-', 'Color',palette(j), 'LineWidth',2);
end
for j = 1:numel(sols)
    plot(sols{j}.x1, sols{j}.valx1(1,:), mk{mod(j-1,5)+1}, 'Color',palette(j), 'MarkerSize',6);
end
hold off;
xlim([0 1]); ylim([0 1]); set(gca,'XTick',0:0.2:1,'YTick',0:0.1:1);
legend(h, names, 'Location','northwest');
nice_axes(); xlabel('t');
title(sprintf('Fig. 4.  Problem 2: solutions (bvpsuite2.0, order %d, %d-point start mesh)', order, Nstart));
end

function plot_fig5(sols, V, bc, Nstart, order)
% Fig. 5: Problem 3 for one boundary condition, one curve per interval length L
hold on;  mk = {'+','.','x','o','s'};  h = zeros(1,numel(sols));  lab = cell(1,numel(sols));  Lmax = 0;
for j = numel(sols):-1:1
    h(j) = plot(sols{j}.x1tau, sols{j}.valx1tau(1,:), '-', 'Color',palette(j+1), 'LineWidth',1.5);
    Lmax = max(Lmax, V(j).par.L);
end
for j = 1:numel(sols)
    plot(sols{j}.x1, sols{j}.valx1(1,:), mk{mod(j-1,5)+1}, 'Color',palette(j+1), 'MarkerSize',6);
    lab{j} = V(j).label;
end
hold off;
xlim([0 Lmax]); ylim([0 2.5]); set(gca,'YTick',0:0.5:2.5);
legend(h, lab, 'Location','northeast');
bcn = {'y''(0) = 0,  y(L) + y''(L) = 0','y''(0) = 0,  y(L) = 0'};
nice_axes(); xlabel('t');
title(sprintf('Fig. 5.  Problem 3, %s   (bvpsuite2.0, order %d, %d-point start mesh)', bcn{bc}, order, Nstart));
end

function plot_fig7(sols, xis, Nstart, order)
% Fig. 7: Problem 4 on tau in [0,1], one curve per xi
hold on;  h = zeros(1,numel(sols));  lab = cell(1,numel(sols));
for j = 1:numel(sols)
    h(j) = plot(sols{j}.x1tau, sols{j}.valx1tau(1,:), 'Color',palette(j), 'LineWidth',2);
    lab{j} = sprintf('\\xi = %g', xis(j));
end
hold off;
xlim([0 1]); ylim([-1 0.7]); set(gca,'XTick',0:0.2:1,'YTick',-1:0.2:0.6);
legend(h, lab, 'Location','northwest');
nice_axes(); xlabel('\tau');
title(sprintf('Fig. 7.  Problem 4: solutions (bvpsuite2.0, order %d, %d-point start mesh)', order, Nstart));
end


%% ========================================================================================
%%  SECTION 6.   SMALL HELPERS
%% ========================================================================================

function out = emit(fmt, varargin)
% printf to the screen and keep a copy for the log file.  emit('@reset') / emit('@get')
persistent LOG
if isempty(LOG), LOG = ''; end
out = [];
if strcmp(fmt,'@reset'), LOG = '';  return; end
if strcmp(fmt,'@get'),   out = LOG; return; end
str = sprintf(fmt, varargin{:});
LOG = [LOG str];
fprintf('%s', str);
end

function v = version_string()
if exist('OCTAVE_VERSION','builtin')
    v = ['GNU Octave ' OCTAVE_VERSION];
else
    v = ['MATLAB ' version];
end
end

%% ========================================================================================
%%  SECTION 7.   POLYNOMIAL COEFFICIENTS OF THE COLLOCATION SOLUTION
%% ========================================================================================
% On every mesh interval [tau_i, tau_{i+1}] the collocation solution is a polynomial of degree
% m+1 (m Gauss points, second-order ODE).  It is reported in the LOCAL variable s = t - tau_i:
%        y(t) = c0 + c1 s + c2 s^2 + ... + c_{m+1} s^{m+1},      s = t - tau_i

function out = coef_hook(sol, problem, settings)
% called by the solver once per mesh-adaptation iteration; keeps a snapshot of the coefficients
persistent SNAP
if isempty(SNAP), SNAP = {}; end
out = [];
if ischar(sol)
    if strcmp(sol,'@reset'), SNAP = {}; end
    if strcmp(sol,'@take'),  out = {SNAP}; end        % cell wrapper: stored as one struct field
    return
end
api = bvp_main('api');
rho = api.getStandardCollocationPoints(feval(settings,'collMethod'), feval(settings,'collPoints'));
l   = feval(problem,'orders');
snap.x1 = sol.x1(:).';
snap.C  = interval_polys(sol.coeff, sol.x1, l(1), rho);
SNAP{end+1} = snap;
end

function C = interval_polys(coeff, x1, l, rho)
% N x (m+l) matrix, row i = [c0 c1 ... c_{m+l-1}] of the polynomial on interval i (scalar ODE of order l)
m = numel(rho);  N = numel(x1)-1;  len = m+l;
C = zeros(N,len);
for i = 1:N
    h = x1(i+1) - x1(i);
    base = (l+m)*(i-1);
    poly = zeros(1,len);                                   % descending powers of s
    for k = 1:l                                            % Taylor part at tau_i
        q = zeros(1,len);  q(len-k+1) = 1/factorial(k-1);
        poly = poly + coeff(base+k)*q;
    end
    for kk = 1:m                                           % collocation part (Lagrange basis integrated l times)
        pk = 1;
        for j = 1:m
            if j ~= kk, pk = conv(pk,[1 -rho(j)]/(rho(kk)-rho(j))); end
        end
        for j = 1:l, pk = polyint(pk); end                 % polynomial in x = s/h, length m+l
        pw = (len-1):-1:0;                                 % powers of x of the entries of pk
        poly = poly + coeff(base+l+kk) * h^l * (pk ./ h.^pw);
    end
    C(i,:) = fliplr(poly);
end
end

function txt = format_coeffs(k, V, runs)
% text block with the coefficients of every run and every iteration of problem k
txt = sprintf('\n%s\n POLYNOMIAL COEFFICIENTS, problem %d:   y(t) = c0 + c1 s + c2 s^2 + ...   on [tau_i, tau_{i+1}],   s = t - tau_i\n%s\n', ...
              repmat('=',1,96), k, repmat('=',1,96));
for r = 1:numel(runs)
    R = runs(r);
    txt = [txt sprintf('\n##### case %s | start mesh %d | order %d | %d iteration(s)\n', ...
                       caselabel(R.label), R.start, R.order, numel(R.iters))];
    for it = 1:numel(R.iters)
        S = R.iters{it};  [N,nc] = size(S.C);
        txt = [txt sprintf('--- iteration %d: %d intervals, %d coefficients each\n', it, N, nc)];
        txt = [txt sprintf('%5s %14s %14s', 'i','tau_i','tau_i+1')];
        for c = 0:nc-1, txt = [txt sprintf(' %20s', sprintf('c%d',c))]; end
        txt = [txt sprintf('\n')];
        for i = 1:N
            txt = [txt sprintf('%5d %14.8g %14.8g', i, S.x1(i), S.x1(i+1))];
            txt = [txt sprintf(' %20.12e', S.C(i,:)) sprintf('\n')];
        end
    end
end
end

function s = caselabel(l)
if isempty(l), s = '-'; else, s = l; end
end

function out = coef_text(str)
% keeps the coefficient text of the whole run for the coefficient file
persistent BUF
if isempty(BUF), BUF = ''; end
out = [];
if strcmp(str,'@reset'), BUF = ''; return; end
if strcmp(str,'@get'),   out = BUF; return; end
BUF = [BUF str];
end
