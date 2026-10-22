function out = problems(request, varargin)
% PROBLEMS  The four test problems of
%
%   M. Hohenegger, G. Settanni, E.B. Weinmueller, M. Wolde,
%   "Numerical treatment of singular ODEs using finite difference and collocation methods",
%   Appl. Numer. Math. 205 (2024) 184-194,
%
% together with ALL their settings: initial (start) mesh points, collocation orders, tolerances,
% problem parameters, and the plot / print switches.  Edit SECTION 1 below.
%
% This file only DESCRIBES the problems; it does not solve anything.
%   run_paper_results.m and run_dbn_comparison.m take the problems and settings from here and
%   solve them with bvp_main.m.
%
% Requests (used by those two files)
%   cfg = problems('config')                   % cfg.P = {P1,P2,P3,P4},  cfg.G = general settings
%   V   = problems('variants', k)              % cases of problem k (e.g. L = 5, 20 and bc = 1, 2)
%   B   = problems('build', k, variant, Nmesh, order, TOL, adaptFcn, extra)
%         B.problem, B.settings  are the two arguments of  bvp_main(problem, settings)
%
% Works in MATLAB and GNU Octave.

%% ========================================================================================
%%  SECTION 1.   SETTINGS OF ALL PROBLEMS  (edit here)
%% ========================================================================================

SHOW_DEFAULT = ~exist('OCTAVE_VERSION','builtin');   % open figure windows? (MATLAB yes, Octave no)

% ---- general -------------------------------------------------------------------------------
G.figDir  = 'figures';                       % PNG files go to this folder (next to this file)
G.logFile = 'paper_results_output.txt';      % everything printed is also written to this file
G.coeffFile = 'paper_coefficients_output.txt'; % polynomial coefficients of every iteration (see printCoeffs)
G.coeffsOnScreen = true;                     % true: coefficients are also printed on screen / in the log
                                             % false: only written to G.coeffFile (they are long)

% ---- common meaning of the per-problem fields below ----------------------------------------
%   enable     true/false   run this problem when run_paper_results is called without arguments
%   startMesh  vector       numbers of points of the equidistant start mesh; one run for each entry
%   orders     vector       numbers m of Gauss collocation points;           one run for each entry
%   TOL        scalar       absolute = relative tolerance of the mesh adaptation (mixed error control)
%   solTol     scalar       tolerance of the Newton iteration
%   plots      true/false   draw this problem's figure(s)
%   show       true/false   show them on screen   (false: only save the PNG files)
%   save       true/false   save them as PNG in G.figDir
%   printCoeffs true/false  print the polynomial coefficients of the numerical solution on every
%                           interval, for EVERY mesh-adaptation iteration of every run
%                           (printed after the problem's table; also written to G.coeffFile)
% The figures are made from the FIRST entry of startMesh and the FIRST entry of orders
% (except the mesh figure of Problem 1, which is made for every entry of startMesh).

% ======================= PROBLEM 1  (Table 1, Figs. 2 and 3) ===============================
%   y'' - y'/t - 2 = 0,   y(0) = y(1) = 0,    exact solution  y = t^2 ln t
P1.title     = 'y'''' - y''/t - 2 = 0,   y(0) = y(1) = 0,   exact solution y = t^2 ln t   (Table 1)';
P1.enable    = true;
P1.startMesh = [11 151];
P1.orders    = [6 8 10];
P1.TOL       = 1e-8;
P1.solTol    = 1e-12;
P1.exact     = @(t) t.^2 .* log(max(t,realmin));    % exact solution (used for the "exact" columns)
P1.exactMax  = 1/(2*exp(1));                        % max |y_exact| = 0.18394
P1.plots     = true;   P1.show = SHOW_DEFAULT;   P1.save = true;
P1.printCoeffs = true;

% ======================= PROBLEM 2  (Table 4, Fig. 4) ======================================
%   y'' + 2/t y' - phi^2 y exp( gamma*beta*(1-y)/(1+beta*(1-y)) ) = 0,   y'(0) = 0,  y(1) = 1
%   Three solutions exist; the constant starting profile y = inits(j) selects which one is found
%   (y = 0 -> y1,  y = 1 -> y3;  bvpsuite2.0 does not find y2).
P2.title     = 'y'''' + 2/t y'' - phi^2 y exp(gamma beta (1-y)/(1+beta (1-y))) = 0,   y''(0) = 0, y(1) = 1   (Table 4)';
P2.enable    = true;
P2.startMesh = 151;
P2.orders    = [6 8 10];
P2.TOL       = 1e-8;
P2.solTol    = 1e-12;
P2.phi = 0.6;   P2.gamma = 40;   P2.beta = 0.2;       % problem parameters
P2.inits     = [0 1];                                 % constant starting profiles ...
P2.initNames = {'y1(t)','y3(t)'};                     % ... and the names of the solutions they give
P2.plots     = true;   P2.show = SHOW_DEFAULT;   P2.save = true;
P2.printCoeffs = true;

% ======================= PROBLEM 3  (Table 6, Fig. 5) ======================================
%   y'' + 4/t y' + (t y - 1) y = 0,   y'(0) = 0,   boundary condition at t = L:
%       bc = 1:  y(L) + y'(L) = 0        bc = 2:  y(L) = 0
%   Starting profile: y = 2 for t <= 1.5,  y = 2 exp(1.5 - t) for t > 1.5.
P3.title     = 'y'''' + 4/t y'' + (t y - 1) y = 0,   y''(0) = 0,  condition at t = L   (Table 6)';
P3.enable    = true;
P3.startMesh = 31;
P3.orders    = [6 8 10];
P3.TOL       = 1e-8;
P3.solTol    = 1e-12;
P3.L         = [5 20];                                % interval lengths
P3.bc        = [1 2];                                 % boundary conditions (see above)
P3.plots     = true;   P3.show = SHOW_DEFAULT;   P3.save = true;
P3.printCoeffs = true;

% ======================= PROBLEM 4  (Table 8, Fig. 7) ======================================
%   y'' + (N-1)/t y' - 4 lam^2 (y+1) y (y-xi) = 0,   y'(0) = 0,   y(inf) = xi
%   solved on  tau = 1 - 1/sqrt(1+t) in [0,1]  exactly as in the paper.
%   Starting profiles (accurate reference solutions) are built in for xi = 0.1 and xi = 0.6 only.
P4.title     = 'y'''' + (N-1)/t y'' - 4 lam^2 (y+1) y (y-xi) = 0,   y''(0) = 0, y(inf) = xi   (Table 8)';
P4.enable    = true;
P4.startMesh = 31;
P4.orders    = [6 8 10];
P4.TOL       = 1e-8;
P4.solTol    = 1e-10;               % 1e-12 sits at the round-off floor of this problem
P4.N = 3;   P4.lam = 1;             % problem parameters
P4.xi        = [0.1 0.6];
P4.plots     = true;   P4.show = SHOW_DEFAULT;   P4.save = true;
P4.printCoeffs = true;

%% ========================================================================================
%%  END OF SECTION 1.   Nothing below needs to be edited.
%% ========================================================================================

PROBS = {P1, P2, P3, P4};
switch request
    case 'config'
        out = struct('P', {PROBS}, 'G', G);
    case 'variants'
        k = varargin{1};
        out = variants_of(k, PROBS{k});
    case 'build'
        [k, variant, Nmesh, order, TOL, adaptFcn, extra] = varargin{:};
        [out.problem, out.settings] = build_case(k, PROBS{k}, variant, Nmesh, order, TOL, adaptFcn, extra);
    otherwise
        error('problems: unknown request ''%s''.', request);
end
end


%% ========================================================================================
%%  SECTION 2.   CASES AND HANDLES FOR THE SOLVER
%% ========================================================================================

function V = variants_of(k, P)
% The list of "cases" of problem k (e.g. the two boundary conditions x two interval lengths of
% Problem 3).  Each case has a label, a group heading for the printed table and its parameters.
V = struct('label',{},'group',{},'par',{});
switch k
    case 1
        V(1) = struct('label','','group','','par',struct());
    case 2
        for j = 1:numel(P.inits)
            par = struct('phi',P.phi,'gamma',P.gamma,'beta',P.beta,'init',P.inits(j));
            V(end+1) = struct('label',P.initNames{j},'group','','par',par);     %#ok<AGROW>
        end
    case 3
        bcname = {'y''(0)=0,  y(L)+y''(L)=0','y''(0)=0,  y(L)=0'};
        for b = P.bc
            for L = P.L
                par = struct('L',L,'bc',b);
                V(end+1) = struct('label',sprintf('L = %g',L),'group',['boundary condition ' num2str(b) ':  ' bcname{b}],'par',par); %#ok<AGROW>
            end
        end
    case 4
        for xi = P.xi
            par = struct('N',P.N,'lam',P.lam,'xi',xi);
            V(end+1) = struct('label',sprintf('xi = %g',xi),'group','','par',par);   %#ok<AGROW>
        end
end
end


function [problem, settings] = build_case(k, P, variant, Nmesh, order, TOL, adaptFcn, extra)
% Handles that bvp_main() accepts as "problem" and "settings".
%   adaptFcn  '' = original bvpsuite2.0 mesh adaptation, or a function handle of a plug-in routine
%   extra     struct of additional settings requests (used by run_dbn_comparison.m)
S = struct('Nmesh',Nmesh,'order',order,'TOL',TOL,'solTol',P.solTol,'adaptFcn',adaptFcn,'extra',extra);
settings = @(request) paper_settings(request, S);
par = variant.par;
switch k
    case 1, problem = @(varargin) bvp_problem1(varargin{:});
    case 2, problem = @(varargin) bvp_problem2(par, varargin{:});
    case 3, problem = @(varargin) bvp_problem3(par, varargin{:});
    case 4, problem = @(varargin) bvp_problem4(par, varargin{:});
end
end

%% ========================================================================================
%%  SECTION 3.   SOLVER SETTINGS  (same for all problems; values come from section 1)
%% ========================================================================================

function ret = paper_settings(request, S)
switch request
    case 'mesh',                  ret = linspace(0,1,S.Nmesh);     % equidistant start mesh
    case 'collMethod',            ret = 'gauss';
    case 'collPoints',            ret = S.order;                   % m Gauss points per subinterval
    case 'meshAdaptation',        ret = 1;
    case 'errorEstimate',         ret = 1;
    case 'absTolSolver',          ret = S.solTol;
    case 'relTolSolver',          ret = S.solTol;
    case 'absTolMeshAdaptation',  ret = S.TOL;
    case 'relTolMeshAdaptation',  ret = S.TOL;
    case 'meshAdaptationFcn',     ret = S.adaptFcn;                % '' = original routine
    case 'iterationHook'                                           % called once per adaptation iteration
        if isfield(S.extra,'iterationHook'), ret = S.extra.iterationHook; else, ret = []; end
    case 'minInitialMesh',        ret = 2;
    case 'finemesh',              ret = 0;
    case 'allowTRM',              ret = 1;
    case 'maxFunEvalsTRM',        ret = 90000000;
    case 'maxIterationsTRM',      ret = 90000000;
    case 'lambdaMin',             ret = 0.001;
    case 'maxAdaptations',        ret = 18;
    case 'switchToFFNFactor',     ret = 0.5;
    case 'updateJacFactor',       ret = 0.5;
    case 'K',                     ret = 200;
    % pathfollowing settings (unused here)
    case 'thetaMax',              ret = 0.05;
    case 'maxCorrSteps',          ret = 5;
    case 'maxSteplengthGrowth',   ret = 4;
    case 'angleMin',              ret = 0.75;
    case 'meshFactorMax',         ret = 2;
    case 'PredLengthFactor',      ret = 1;
    case 'CorrLengthGrowth',      ret = 4;
    otherwise
        if isfield(S.extra, request)
            ret = S.extra.(request);                               % extra requests of plug-ins
        else
            error('paper_settings:unknown','unknown settings request ''%s''', request);
        end
end
end


%% ========================================================================================
%%  SECTION 4.   THE FOUR PROBLEMS  (bvpsuite2.0 problem definitions)
%%               Every problem is  F(z) = 0  with z(1,1) = y, z(1,2) = y', z(1,3) = y''.
%% ========================================================================================

% ------------------------------------------------------------------------------ Problem 1
function ret = bvp_problem1(request,z,za,zb,zc,t,p,lambda,lambda_p)  %#ok<INUSD>
% y'' - y'/t - 2 = 0,   y(0) = 0,  y(1) = 0
switch request
    case 'n',        ret = 1;
    case 'orders',   ret = 2;
    case 'problem',  ret = [ z(1,3) - z(1,2)/t - 2 ];
    case 'jacobian'
        ret = zeros(1,1,3);
        ret(1,1,2) = -1/t;
        ret(1,1,3) = 1;
    case 'interval', ret = [0,1];
    case 'linear',   ret = 1;
    case 'parameters', ret = 0;
    case 'c',        ret = [];
    case 'BV',       ret = [ za(1,1);
                             zb(1,1) ];
    case 'dBV'
        ret = zeros(2,1,1,2);
        ret(1,1,1,1) = 1;
        ret(2,2,1,1) = 1;
    case 'dP',       ret = [];
    case 'dP_BV',    ret = [];
    case 'initProfile'
        ret.initialMesh   = linspace(0,1,11);
        ret.initialValues = ones(1,11);
    case 'EVP',      ret = 0;
    case 'dLambda',  ret = 0;
    case 'pathfollowing', ret.activate = 0;
    otherwise,       ret = 0;
end
end

% ------------------------------------------------------------------------------ Problem 2
function ret = bvp_problem2(par,request,z,za,zb,zc,t,p,lambda,lambda_p)  %#ok<INUSD>
% y'' + 2/t y' - phi^2 y exp( gamma*beta*(1-y)/(1+beta*(1-y)) ) = 0,   y'(0) = 0,  y(1) = 1
phi = par.phi;  gamma = par.gamma;  beta = par.beta;
switch request
    case 'n',        ret = 1;
    case 'orders',   ret = 2;
    case 'problem'
        y = z(1,1);
        E = gamma*beta*(1-y)/(1+beta*(1-y));
        ret = [ z(1,3) + 2/t*z(1,2) - phi^2*y*exp(E) ];
    case 'jacobian'
        ret = zeros(1,1,3);
        y  = z(1,1);
        D  = 1+beta*(1-y);
        E  = gamma*beta*(1-y)/D;
        dE = -gamma*beta/D^2;
        ret(1,1,1) = -phi^2*exp(E)*(1 + y*dE);
        ret(1,1,2) = 2/t;
        ret(1,1,3) = 1;
    case 'interval', ret = [0,1];
    case 'linear',   ret = 0;
    case 'parameters', ret = 0;
    case 'c',        ret = [];
    case 'BV',       ret = [ za(1,2);
                             zb(1,1)-1 ];
    case 'dBV'
        ret = zeros(2,1,1,2);
        ret(1,1,1,2) = 1;
        ret(2,2,1,1) = 1;
    case 'dP',       ret = [];
    case 'dP_BV',    ret = [];
    case 'initProfile'                       % constant starting profile y = par.init
        x = linspace(0,1,51);
        ret.initialMesh   = x;
        ret.initialValues = par.init*ones(1,51);
    case 'EVP',      ret = 0;
    case 'dLambda',  ret = 0;
    case 'pathfollowing', ret.activate = 0;
    otherwise,       ret = 0;
end
end

% ------------------------------------------------------------------------------ Problem 3
function ret = bvp_problem3(par,request,z,za,zb,zc,t,p,lambda,lambda_p)  %#ok<INUSD>
% y'' + 4/t y' + (t y - 1) y = 0,   y'(0) = 0,   on [0,L]
%   par.bc = 1:  y(L) + y'(L) = 0       par.bc = 2:  y(L) = 0
switch request
    case 'n',        ret = 1;
    case 'orders',   ret = 2;
    case 'problem',  ret = [ z(1,3) + 4/t*z(1,2) + (t*z(1,1) - 1)*z(1,1) ];
    case 'jacobian'
        ret = zeros(1,1,3);
        ret(1,1,1) = 2*t*z(1,1) - 1;
        ret(1,1,2) = 4/t;
        ret(1,1,3) = 1;
    case 'interval', ret = [0, par.L];
    case 'linear',   ret = 0;
    case 'parameters', ret = 0;
    case 'c',        ret = [];
    case 'BV'
        if par.bc == 1
            ret = [ za(1,2);
                    zb(1,1) + zb(1,2) ];
        else
            ret = [ za(1,2);
                    zb(1,1) ];
        end
    case 'dBV'
        ret = zeros(2,1,1,2);
        ret(1,1,1,2) = 1;
        ret(2,2,1,1) = 1;
        if par.bc == 1
            ret(2,2,1,2) = 1;
        end
    case 'dP',       ret = [];
    case 'dP_BV',    ret = [];
    case 'initProfile'                       % y = 2 (t <= 1.5),  y = 2 exp(1.5 - t) (t > 1.5)
        x = linspace(0,par.L,51);
        ret.initialMesh   = x;
        ret.initialValues = (x<=1.5).*2 + (x>1.5).*2.*exp(1.5-x);
    case 'EVP',      ret = 0;
    case 'dLambda',  ret = 0;
    case 'pathfollowing', ret.activate = 0;
    otherwise,       ret = 0;
end
end

% ------------------------------------------------------------------------------ Problem 4
function ret = bvp_problem4(par,request,z,za,zb,zc,t,p,lambda,lambda_p)  %#ok<INUSD>
% y'' + (N-1)/t y' - 4 lam^2 (y+1) y (y-xi) = 0,   y'(0) = 0,   y(inf) = xi.
% Transformed to the finite domain  tau = 1 - 1/sqrt(1+t),  tau in [0,1]
% (the independent variable of the code is tau, called t here):
%   1/4 (1-tau)^6 y'' + ( (N-1)/(2 tau (2-tau)) - 3/4 ) (1-tau)^5 y'
%       - 4 lam^2 (y+1) y (y-xi) = 0,      y'(0) = 0,   y(1) = xi.
N = par.N;  lam = par.lam;  xi = par.xi;
switch request
    case 'n',        ret = 1;
    case 'orders',   ret = 2;
    case 'problem'
        y = z(1,1);
        ret = [ 0.25*(1-t)^6*z(1,3) ...
                + ((N-1)/(2*t*(2-t)) - 3/4)*(1-t)^5*z(1,2) ...
                - 4*lam^2*(y+1)*y*(y-xi) ];
    case 'jacobian'
        ret = zeros(1,1,3);
        y = z(1,1);
        ret(1,1,1) = -4*lam^2*(3*y^2 + 2*(1-xi)*y - xi);
        ret(1,1,2) = ((N-1)/(2*t*(2-t)) - 3/4)*(1-t)^5;
        ret(1,1,3) = 0.25*(1-t)^6;
    case 'interval', ret = [0,1];
    case 'linear',   ret = 0;
    case 'parameters', ret = 0;
    case 'c',        ret = [];
    case 'BV',       ret = [ za(1,2);
                             zb(1,1) - xi ];
    case 'dBV'
        ret = zeros(2,1,1,2);
        ret(1,1,1,2) = 1;
        ret(2,2,1,1) = 1;
    case 'dP',       ret = [];
    case 'dP_BV',    ret = [];
    case 'initProfile'
        % Starting profile on tau = linspace(0,1,101): accurate reference solution of the original
        % problem (computed once by shooting from y(0)=a, y'(0)=0), mapped to tau.
        % (The paper takes its start profile from its reference [20].)
        tau = linspace(0,1,101);
        if abs(xi-0.1) < 1e-12
            yt = [-0.3046629136, -0.3046393539, -0.3045657623, -0.3044375009, -0.3042495504, -0.30399648, -0.3036724146, -0.303271001, -0.3027853704, -0.3022080998, -0.301531171, -0.3007459265, -0.2998430247, -0.2988123923, -0.2976431756, -0.2963236906, -0.2948413729, -0.2931827274, -0.2913332801, -0.2892775322, -0.2869989189, -0.2844797751, -0.2817013097, -0.2786435943, -0.2752855671, -0.2716050596, -0.2675788513, -0.2631827589, -0.2583917688, -0.2531802237, -0.2475220716, -0.241391193, -0.2347618166, -0.2276090404, -0.219909471, -0.2116419955, -0.2027886974, -0.1933359236, -0.1832755036, -0.1726061107, -0.1613347442, -0.1494782895, -0.1370650951, -0.124136477, -0.1107480357, -0.0969706419, -0.08289092702, -0.06861110504, -0.05424795607, -0.03993083506, -0.02579862682, -0.01199565713, 0.001333317294, 0.01404879113, 0.02602249527, 0.03714322133, 0.0473219791, 0.05649600284, 0.06463121426, 0.07172290342, 0.07779457884, 0.0828951342, 0.08709464865, 0.09047925445, 0.09314555779, 0.09519508518, 0.09672916256, 0.09784453606, 0.09862993219, 0.09916364846, 0.09951217363, 0.09972976355, 0.0998588439, 0.0999310711, 0.099968856, 0.09998714049, 0.09999522022, 0.09999842989, 0.0999995544, 0.09999989378, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
        elseif abs(xi-0.6) < 1e-12
            yt = [-0.9953000352, -0.9952979839, -0.9952915723, -0.9952803823, -0.9952639493, -0.9952417563, -0.9952132254, -0.9951777097, -0.9951344825, -0.9950827249, -0.9950215114, -0.9949497925, -0.9948663738, -0.9947698917, -0.9946587834, -0.9945312509, -0.9943852177, -0.9942182764, -0.9940276234, -0.9938099803, -0.9935614962, -0.9932776263, -0.9929529819, -0.9925811411, -0.9921544117, -0.9916635305, -0.9910972803, -0.9904419994, -0.9896809491, -0.9887934931, -0.9877540248, -0.9865305578, -0.9850828578, -0.983359952, -0.9812967831, -0.9788096839, -0.975790215, -0.9720967199, -0.967542687, -0.9618806392, -0.954779774, -0.9457949333, -0.9343237372, -0.919548057, -0.9003559581, -0.8752421654, -0.8421918811, -0.7985701141, -0.7410758106, -0.6658866855, -0.5692101086, -0.4484994806, -0.3044045866, -0.142854857, 0.02429961122, 0.1817457737, 0.3158415891, 0.4192383202, 0.4919613902, 0.5391040485, 0.5675498658, 0.5836428251, 0.5922121726, 0.5965096548, 0.5985349168, 0.5994278155, 0.5997937495, 0.5999320341, 0.599979753, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6];
        else
            error('Problem 4: a starting profile is built in only for xi = 0.1 and xi = 0.6 (got xi = %g).', xi);
        end
        ret.initialMesh   = tau;
        ret.initialValues = yt;
    case 'EVP',      ret = 0;
    case 'dLambda',  ret = 0;
    case 'pathfollowing', ret.activate = 0;
    otherwise,       ret = 0;
end
end
