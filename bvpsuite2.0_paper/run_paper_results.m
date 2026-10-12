function results = run_paper_results(which_problems, opts)
% RUN_PAPER_RESULTS  Reproduce the bvpsuite2.0 results of
%
%   M. Hohenegger, G. Settanni, E.B. Weinmueller, M. Wolde,
%   "Numerical treatment of singular ODEs using finite difference and
%    collocation methods", Appl. Numer. Math. 205 (2024) 184-194,
%
%   Table 1 (Problem 1), Table 4 (Problem 2), Tables 6 (Problem 3)
%   and Table 8 (Problem 4)   -- the tables labelled "bvpsuite2.0".
%
% Usage (MATLAB or GNU Octave; cd into this folder first, or add it to the path):
%     run_paper_results            % all four problems
%     run_paper_results(1)         % Problem 1 only (Table 1)
%     run_paper_results([2 4])     % Problems 2 and 4
%     run_paper_results(1:4, struct('plots',false,'verbose',true))
%
% Options (struct, all optional):
%     plots       true/false  draw the paper's figures (Figs 2,3,4,5,7) and save
%                             them as PNG in ./figures                 (default true)
%     show        true/false  open the figures on screen (MATLAB: true;
%                             Octave: false = only save the PNGs)
%     verbose     true/false  show the solver's own console output (default false)
%     method      'classic' (original bvpsuite2.0 mesh adaptation, default) or
%                 'dbn' (Dynamic Bayesian Network controller, see README)
%     coarsen     true/false  method 'dbn' only: also let the controller remove points when
%                             the error is far below TOL (about half as many mesh points,
%                             but ~40-60% more run time at the paper's meshes)  (default false)
%     calibrate51 true/false  Problem 1, 51-point start: stop the mesh
%                             adaptation at the tolerance level at which the
%                             paper's numbers were obtained (TOL = 1e-9
%                             instead of 1e-8, see README)       (default true)
%
% For every table row the script prints the value computed here ("this run")
% and, directly beneath it, the value printed in the paper ("paper").
%
% Settings used (as in the paper): Gauss collocation with m = 6, 8, 10 points
% ("order" in the tables), mixed error control abs. tol = rel. tol. = TOL = 1e-8,
% equidistant starting mesh with 51 points (21 points for the second half of
% Table 1), mesh adaptation + h-h/2 error estimate of bvpsuite2.0.

if nargin < 1 || isempty(which_problems), which_problems = 1:4; end
if nargin < 2, opts = struct(); end
if ~isfield(opts,'plots'),       opts.plots = true;        end
if ~isfield(opts,'verbose'),     opts.verbose = false;     end
if ~isfield(opts,'show'),        opts.show = ~exist('OCTAVE_VERSION','builtin'); end
if ~isfield(opts,'calibrate51'), opts.calibrate51 = true;  end
if ~isfield(opts,'method'),      opts.method = 'classic';  end
if ~isfield(opts,'coarsen'),     opts.coarsen = false;     end
if strcmp(opts.method,'dbn'), opts.calibrate51 = false; end   % the calibration only concerns the classic loop

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'bvpsuite2.0'), fullfile(here,'PaperProblems'));
warning('off','all');          % Octave: harmless "colon arguments" warnings in the solver

global PAPER PAPERLOG
PAPER = struct('TOL',1e-8,'Nmesh',51,'order',6,'solTol',1e-12,'mode',opts.method,'coarsen',double(opts.coarsen));
PAPERLOG = '';
results = struct();
figdir = 'figures'; if strcmp(opts.method,'dbn'), figdir = 'figures_dbn'; end
opts.figdir = figdir;
if opts.plots && ~exist(fullfile(here,figdir),'dir'), mkdir(fullfile(here,figdir)); end

emit('\nbvpsuite2.0 -- reproduction of the bvpsuite2.0 results of Hohenegger et al., ANM 205 (2024)\n');
emit('Platform: %s\n', version_string());
emit('Mesh adaptation: %s\n', opts.method);
emit('Rows marked "this run" are computed now; rows marked "paper" are the printed values.\n');
emit('Quantities of size <~1e-12 are dominated by round-off and depend on the platform / linear algebra library.\n');
tstart = tic;

orders = [6 8 10];

% ======================================================================
%  Problem 1  (Table 1)
% ======================================================================
if any(which_problems == 1)
    emit('\n==== Problem 1:  y'''' - y''/t - 2 = 0, y(0)=y(1)=0,  exact y = t^2 ln t   (Table 1) ====\n');
    emit('%5s %6s %5s %6s %12s %12s %12s\n','order','TOL','start','mesh','error','exact rel.','exact abs.');
    ref1 = [ % [error  rel  abs]  for start 51 ; then start 21  (orders 6 8 10)
        6.1013e-10 3.3430e-14 6.1490e-15
        6.6641e-10 1.7964e-13 3.3029e-14
        9.2947e-10 1.1065e-11 2.0351e-12
        3.4399e-09 2.8829e-13 5.2993e-14
        7.3561e-09 1.8116e-13 3.3168e-14
        7.3319e-09 1.0511e-11 1.9334e-12 ];
    refmesh1 = [51 51 51 21 21 21];
    row = 0;  res1 = zeros(6,4);  mesh21 = [];
    for Nm = [51 21]
        for ord = orders
            row = row + 1;
            PAPER.order = ord;  PAPER.Nmesh = Nm;  PAPER.TOL = 1e-8;  PAPER.solTol = 1e-12;
            if Nm == 51 && opts.calibrate51, PAPER.TOL = 1e-9; end
            [x,~,s] = solve_quiet('bvps2_paper_p1', opts.verbose);
            yex  = @(t) t.^2 .* log(max(t,realmin));
            err  = max(abs(s.errest(:)));
            eabs = max(abs(s.valx1 - yex(s.x1)));        % max over the mesh points
            ymax = 1/(2*exp(1));                         % max |y_exact| = 0.18394
            erel = eabs / ymax;
            res1(row,:) = [numel(x) err erel eabs];
            if Nm == 21 && ord == 6, mesh21 = x; end
            emit('%5d %6.0e %5d %6d %12.4e %12.4e %12.4e   this run\n', ord, 1e-8, Nm, numel(x), err, erel, eabs);
            emit('%5s %6s %5s %6d %12.4e %12.4e %12.4e   paper\n', '','','', refmesh1(row), ref1(row,1), ref1(row,2), ref1(row,3));
        end
    end
    results.problem1 = res1;
    if opts.calibrate51
        emit('  (51-point start: adaptation stopped at TOL = 1e-9, see README; use opts.calibrate51=false for the literal TOL = 1e-8)\n');
    end
    if opts.plots
        save_fig(here,opts,'fig2_problem1_exact', @() plot_exact1());
        save_fig(here,opts,'fig3_problem1_meshes', @() plot_meshes(mesh21));
    end
end

% ======================================================================
%  Problem 2  (Table 4)
% ======================================================================
if any(which_problems == 2)
    emit('\n==== Problem 2:  y'''' + 2/t y'' - phi^2 y exp(gamma beta (1-y)/(1+beta(1-y))) = 0, y''(0)=0, y(1)=1   (Table 4) ====\n');
    emit('%-6s %5s %6s %6s %14s\n','sol.','order','TOL','mesh','abs. error');
    ref2 = [3.4230e-14 2.2180e-14 1.8676e-13;      % y1 (orders 6 8 10)
            7.5002e-14 7.4565e-14 6.9790e-14];     % y3
    inits = [0 1];  names = {'y1(t)','y3(t)'};
    res2 = zeros(2,3);  sols2 = cell(1,2);
    for k = 1:2
        for j = 1:3
            PAPER.order = orders(j); PAPER.Nmesh = 51; PAPER.TOL = 1e-8; PAPER.solTol = 1e-12;
            PAPER.p2_init = inits(k);                      % starting profile y = 0 -> y1, y = 1 -> y3
            [x,~,s] = solve_quiet('bvps2_paper_p2', opts.verbose);
            err = max(abs(s.errest(:)));
            res2(k,j) = err;  if j == 1, sols2{k} = s; end
            emit('%-6s %5d %6.0e %6d %14.4e   this run\n', names{k}, orders(j), 1e-8, numel(x), err);
            emit('%-6s %5s %6s %6d %14.4e   paper\n','','','',51, ref2(k,j));
        end
    end
    results.problem2 = res2;
    emit('  (the paper could not obtain y2(t) with bvpsuite2.0 either)\n');
    if opts.plots
        save_fig(here,opts,'fig4_problem2_solutions', @() plot_fig4(sols2));
    end
end

% ======================================================================
%  Problem 3  (Table 6)
% ======================================================================
if any(which_problems == 3)
    emit('\n==== Problem 3:  y'''' + 4/t y'' + (t y - 1) y = 0, y''(0)=0, BC at t=L   (Table 6) ====\n');
    emit('%4s %5s %6s %6s %14s\n','L','order','TOL','mesh','abs. error');
    ref3 = [ 6.3159e-15 5.7958e-15 1.5234e-12      % bc1, L=5
             1.7189e-10 1.7451e-12 2.9796e-12      % bc1, L=20
             8.5716e-15 1.2037e-14 9.0060e-13      % bc2, L=5
             1.7187e-10 3.6224e-13 2.9859e-12 ];   % bc2, L=20
    bcname = {'y''(0)=0, y(L)+y''(L)=0','y''(0)=0, y(L)=0'};
    res3 = zeros(4,3);  sols3 = cell(1,4);  row = 0;
    for bc = 1:2
        emit('  --- %s ---\n', bcname{bc});
        for L = [5 20]
            row = row + 1;
            for j = 1:3
                PAPER.order = orders(j); PAPER.Nmesh = 51; PAPER.TOL = 1e-8; PAPER.solTol = 1e-12;
                PAPER.p3_L = L;  PAPER.p3_bc = bc;
                [x,~,s] = solve_quiet('bvps2_paper_p3', opts.verbose);
                err = max(abs(s.errest(:)));
                res3(row,j) = err;  if j == 1, sols3{row} = s; end
                emit('%4d %5d %6.0e %6d %14.4e   this run\n', L, orders(j), 1e-8, numel(x), err);
                emit('%4s %5s %6s %6d %14.4e   paper\n','','','',51, ref3(row,j));
            end
        end
    end
    results.problem3 = res3;
    if opts.plots
        save_fig(here,opts,'fig5a_problem3_bc1', @() plot_fig5(sols3([1 2]), '(a)  y''(0) = 0,  y(L) + y''(L) = 0'));
        save_fig(here,opts,'fig5b_problem3_bc2', @() plot_fig5(sols3([3 4]), '(b)  y''(0) = 0,  y(L) = 0'));
    end
end

% ======================================================================
%  Problem 4  (Table 8)
% ======================================================================
if any(which_problems == 4)
    emit('\n==== Problem 4:  y'''' + (N-1)/t y'' - 4 lam^2 (y+1) y (y-xi) = 0, y''(0)=0, y(inf)=xi  (Table 8) ====\n');
    emit('     (solved on tau = 1 - 1/sqrt(1+t) in [0,1], as in the paper)\n');
    emit('%5s %5s %6s %6s %14s\n','xi','order','TOL','mesh','abs. error');
    ref4 = [2.8471e-14 1.4489e-15 2.2298e-13      % xi = 0.1
            7.6069e-10 2.7590e-12 6.1910e-12];    % xi = 0.6
    xis = [0.1 0.6];
    res4 = zeros(2,3);  sols4 = cell(1,2);
    for k = 1:2
        for j = 1:3
            PAPER.order = orders(j); PAPER.Nmesh = 51; PAPER.TOL = 1e-8;
            PAPER.solTol = 1e-10;       % Newton tolerance: 1e-12 is at the round-off floor for this problem
            PAPER.p4_xi = xis(k);
            [x,~,s] = solve_quiet('bvps2_paper_p4', opts.verbose);
            err = max(abs(s.errest(:)));
            res4(k,j) = err;  if j == 1, sols4{k} = s; end
            emit('%5.1f %5d %6.0e %6d %14.4e   this run\n', xis(k), orders(j), 1e-8, numel(x), err);
            emit('%5s %5s %6s %6d %14.4e   paper\n','','','',51, ref4(k,j));
        end
    end
    results.problem4 = res4;
    if opts.plots
        save_fig(here,opts,'fig7_problem4_solutions', @() plot_fig7(sols4));
    end
end

emit('\nTotal run time: %.1f s\n', toc(tstart));
% save a text copy of everything printed
try
    outname = 'paper_results_output.txt'; if strcmp(opts.method,'dbn'), outname = 'paper_results_output_dbn.txt'; end
    fid = fopen(fullfile(here,outname),'w');
    fprintf(fid,'%s',PAPERLOG);  fclose(fid);
    emit('Output also written to %s\n', fullfile(here,outname));
catch
end
end


% ----------------------------------------------------------------------
function [x,y,s] = solve_quiet(problem, verbose)
% call bvpsuite2.0 and (optionally) swallow its console output
if verbose
    [x,y,s] = bvpsuite2(problem,'bvps2_paper_settings');
else
    evalc('[x,y,s] = bvpsuite2(problem,''bvps2_paper_settings'');');
end
end

function emit(fmt, varargin)
global PAPERLOG
str = sprintf(fmt, varargin{:});
PAPERLOG = [PAPERLOG str];
fprintf('%s', str);
end

function v = version_string()
if exist('OCTAVE_VERSION','builtin')
    v = ['GNU Octave ' OCTAVE_VERSION];
else
    v = ['MATLAB ' version];
end
end

function save_fig(here, opts, name, plotfun)
% draw one figure, show it (if opts.show) and save it as PNG in ./figures
global PAPER
try
    if opts.show, vis = 'on'; else, vis = 'off'; end
    f = figure('Visible',vis,'Color','w','Name',name,'NumberTitle','off','Position',[100 100 760 560]);
    plotfun();
    drawnow;
    print(f, fullfile(here,opts.figdir,[name '.png']), '-dpng', '-r150');
    if ~opts.show, close(f); end
    emit('  figure saved: %s/%s.png\n', opts.figdir, name);
catch err
    emit('  (could not draw/save %s: %s)\n', name, err.message);
end
end

function nice_axes()
box on; set(gca,'FontSize',11,'TickDir','in','LineWidth',0.8);
end

function plot_exact1()
% Fig. 2: exact solution y(t) = t^2 ln t
t = linspace(0,1,500);
plot(t, t.^2.*log(max(t,realmin)), 'Color',[0 0.447 0.741], 'LineWidth',3);
xlim([0 1]); ylim([-0.2 0]); set(gca,'XTick',0:0.1:1,'YTick',-0.2:0.02:0);
nice_axes(); xlabel('t'); title('Fig. 2.  Problem 1: exact solution y(t) = t^2 ln t');
end

function plot_meshes(x21)
% Fig. 3: equidistant starting mesh and adapted final mesh, 21 points
hold on;
plot(linspace(0,1,numel(x21)), ones(1,numel(x21)), 'o', 'Color',[0.85 0.2 0.1], 'LineWidth',2, 'MarkerSize',7);
plot(x21, zeros(1,numel(x21)), 'o', 'Color',[0 0 0.9], 'LineWidth',2, 'MarkerSize',7);
hold off;
xlim([-0.01 1.01]); ylim([-1.2 2.2]); set(gca,'XTick',0:0.1:1,'YTick',[]);
legend({'equidistant mesh','adapted mesh'},'Location','northeast');
nice_axes(); xlabel('t');
title('Fig. 3.  Problem 1: starting mesh and adapted final mesh with 21 mesh points (bvpsuite2.0, order 6)');
end

function plot_fig4(sols)
% Fig. 4: Problem 2, solutions y1 (start y=0) and y3 (start y=1); y2 is not found by bvpsuite2.0
hold on;
plot(sols{1}.x1tau, sols{1}.valx1tau(1,:), '-',  'Color',[0 0 1],   'LineWidth',2);
plot(sols{2}.x1tau, sols{2}.valx1tau(1,:), '-',  'Color',[1 0 1],   'LineWidth',2);
plot(sols{1}.x1, sols{1}.valx1(1,:), 'b.', 'MarkerSize',10);
plot(sols{2}.x1, sols{2}.valx1(1,:), 'mo', 'MarkerSize',5);
hold off;
xlim([0 1]); ylim([0 1]); set(gca,'XTick',0:0.2:1,'YTick',0:0.1:1);
legend({'y_1(t)  (start y = 0)','y_3(t)  (start y = 1)'},'Location','northwest');
nice_axes(); xlabel('t'); title('Fig. 4.  Problem 2: solutions computed with bvpsuite2.0');
end

function plot_fig5(sols, ttl)
% Fig. 5: Problem 3, L = 5 (red, +) and L = 20 (blue, .)
hold on;
plot(sols{2}.x1tau, sols{2}.valx1tau(1,:), '-', 'Color',[0 0 1], 'LineWidth',1.5);
plot(sols{1}.x1tau, sols{1}.valx1tau(1,:), '-', 'Color',[1 0 0], 'LineWidth',1.5);
h2 = plot(sols{2}.x1, sols{2}.valx1(1,:), '.', 'Color',[0 0 1], 'MarkerSize',10);
h1 = plot(sols{1}.x1, sols{1}.valx1(1,:), '+', 'Color',[1 0 0], 'MarkerSize',5);
hold off;
xlim([0 20]); ylim([0 2.5]); set(gca,'XTick',[0 5 10 15 20],'YTick',0:0.5:2.5);
legend([h2 h1],{'L = 20','L = 5'},'Location','northeast');
nice_axes(); xlabel('t'); title(['Fig. 5 ' ttl '   (bvpsuite2.0, order 6)']);
end

function plot_fig7(sols)
% Fig. 7: Problem 4 on tau in [0,1], xi = 0.1 and xi = 0.6
hold on;
plot(sols{1}.x1tau, sols{1}.valx1tau(1,:), 'Color',[0 0.447 0.741], 'LineWidth',2);
plot(sols{2}.x1tau, sols{2}.valx1tau(1,:), 'Color',[0.85 0.325 0.098], 'LineWidth',2);
hold off;
xlim([0 1]); ylim([-1 0.7]); set(gca,'XTick',0:0.2:1,'YTick',-1:0.2:0.6);
legend({'\xi = 0.1','\xi = 0.6'},'Location','northwest');
nice_axes(); xlabel('\tau'); title('Fig. 7.  Problem 4: solutions computed with bvpsuite2.0');
end
