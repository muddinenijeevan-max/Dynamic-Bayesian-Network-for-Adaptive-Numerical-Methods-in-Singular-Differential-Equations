function [ret] = bvps2_paper_settings(request)
% Solver settings used for all bvpsuite2.0 runs of the paper.
%   order  = global order of the collocation scheme (Gauss: order = 2*m)
%   TOL    = tolerance for the mixed (abs/rel) error control = 1e-8
%   Nmesh  = number of points of the (equidistant) starting mesh (51 or 21)
% These values are set by run_paper_results.m in the global struct PAPER.

global PAPER
if ~isfield(PAPER,'solTol'), PAPER.solTol = 1e-12; end   % Newton tolerance (default of bvpsuite2.0)
if ~isfield(PAPER,'mode'),   PAPER.mode   = 'classic'; end % mesh adaptation: 'classic' or 'dbn'
if ~isfield(PAPER,'coarsen'),PAPER.coarsen = 0;       end % dbn only: allow coarsening probes

switch request
    case 'mesh'
        ret = linspace(0,1,PAPER.Nmesh);
    case 'collMethod'
        ret = 'gauss';
    case 'collPoints'
        ret = PAPER.order;
    case 'meshAdaptation'
        ret = 1;
    case 'errorEstimate'
        ret = 1;
    case 'absTolSolver'
        ret = PAPER.solTol;
    case 'relTolSolver'
        ret = PAPER.solTol;
    case 'absTolMeshAdaptation'
        ret = PAPER.TOL;
    case 'relTolMeshAdaptation'
        ret = PAPER.TOL;
    case 'minInitialMesh'
        ret = 2;
    case 'meshAdaptationMode'
        ret = PAPER.mode;
    case 'dbnCoarsen'
        ret = PAPER.coarsen;
    case 'dbnVerbose'
        ret = isfield(PAPER,'dbnVerbose') && PAPER.dbnVerbose;
    case 'finemesh'
        ret = 0;
    case 'allowTRM'
        ret = 1;
    case 'maxFunEvalsTRM'
        ret = 90000000;
    case 'maxIterationsTRM'
        ret = 90000000;
    case 'lambdaMin'
        ret = 0.001;
    case 'maxAdaptations'
        ret = 18;
    case 'switchToFFNFactor'
        ret = 0.5;
    case 'updateJacFactor'
        ret = 0.5;
    case 'K'
        ret = 200;
    % pathfollowing settings (unused here)
    case 'thetaMax'
        ret = 0.05;
    case 'maxCorrSteps'
        ret = 5;
    case 'maxSteplengthGrowth'
        ret = 4;
    case 'angleMin'
        ret = 0.75;
    case 'meshFactorMax'
        ret = 2;
    case 'PredLengthFactor'
        ret = 1;
    case 'CorrLengthGrowth'
        ret = 4;
end

end
