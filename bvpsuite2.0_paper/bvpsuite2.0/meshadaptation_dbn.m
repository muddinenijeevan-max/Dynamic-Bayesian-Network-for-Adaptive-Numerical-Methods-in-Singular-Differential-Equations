function [finalsol,toleriert,halve] = meshadaptation_dbn(problem,settings,x1,rhoColl,initProfile)
% MESHADAPTATION_DBN  Mesh adaptation for bvpsuite2.0 steered by a Dynamic Bayesian Network.
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

ordnung = feval_problem(problem,'orders');
n   = length(ordnung);
linear = feval_problem(problem,'linear');
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
            [tcol,ycol,sol] = solveLinearProblem(problem,x1_try,rhoColl);
        else
            [tcol,ycol,sol,halve] = solveNonLinearProblem(problem,settings,x1_try,rhoColl,initProfile,1-(j>1),[]);
            if halve ~= 0, solved = false; end
        end
        if solved
            [sol_c,~,halve] = errorestimate(sol,problem,settings,[]);
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
        residual = abs(computeResidual(problem,sol.coeff,par,tau,rhoColl,resmesh,0));
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

