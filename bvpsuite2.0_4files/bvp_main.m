% BVP_MAIN  The bvpsuite2.0 collocation solver for singular BVPs -- ALL files of the package in ONE file.
%
%   [x, y, sol] = bvp_main(problem, settings)
%
%   problem, settings   function handles (or names) with the bvpsuite2.0 conventions; they are
%                       produced by problems.m:   B = problems('build', ...);  bvp_main(B.problem, B.settings)
%   sol                 solution structure (x1, valx1, x1tau, valx1tau, coeff, errest, ...)
%   api = bvp_main('api')   handles to the internal routines (for add-ons such as the DBN controller)
%
% Contents (in this order, each introduced by a banner):
%   bvpsuite2 (entry, renamed bvp_main), default_settings, feval_problem, meshadaptation,
%   solveLinearProblem, solveNonLinearProblem, solve_nonlinear_sys, functionFDF, errorestimate,
%   initial_coefficients, coeffToValues, computeResidual, getStandardCollocationPoints,
%   higherchainrule, trafo, computeEVPStart, pathfollowing, dispDebug, template_bvp.
% Changes against the original package:
%   (a) optional hooks:  settings 'meshAdaptationFcn' (plug-in mesh adaptation, called with the
%       extra argument solver_api()) and settings 'iterationHook' (called once per iteration),
%   (b) helper functions that had the same name in several files are renamed <file>_<name>,
%   (c) comment characters outside ASCII were transliterated.
%
%   bvpsuite2.0   Copyright (c) 2018  W. Auzinger, M. Fallahpour, G. Kitzhofer, O. Koch,
%                 G. Pulverer, G. Soederlind, E.B. Weinmueller, S. Wurm,  TU Wien
%                 Copyright (c) 2018, Ewa Weinmueller.  All rights reserved.
%                 Redistribution and use in source and binary forms, with or without modification,
%                 are permitted provided that the copyright notice, the conditions and the
%                 disclaimer of the original bvpsuite2.0 license.txt are retained.
%                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%                 ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED (see license.txt).

function [x1,valx1,solution] = bvp_main(problem,settings,initProfile,pathfoll)

% Matlab Code BVPSUITE2.0 (beta-version July 2018).
% Authors: Winfried Auzinger, Merlin Fallahpour, Georg Kitzhofer,
% Othmar Koch, Gernot Pulverer, Gustaf S??derlind, Ewa B. Weinmueller and
% Stefan Wurm.
% Institute for Analysis and Scientific Computing, Vienna University of
% Technology, Vienna, Austria.
%
% The package BVPSUITE2.0 has been developed at the Institute for Analysis
% and Scientific Computing, Vienna University of Technology, and can be
% used for the numerical solution of implicit boundary value problems
% (BVPs) in ordinary differential equations (ODEs) as well as eigenvalue
% problems (EVPs), and Index-1 Differential-Algebraic Equations (DAEs). The
% ODE system can have a general implicit form and be of mixed order subject
% to multi-point boundary conditions. Furthermore the interval can be
% finite or semi-infinite. The BVP may also include unknown parameters to
% be calculated together with the unknown solution z, in which case
% additional boundary conditions are required.
% BVPSUITE2.0 is a product of years of research at the Institute for
% Analysis and Scientific computing on the analysis and numerical solution
% of ODEs with time and space singularities. Through the years, the code
% evolved from SBVP1.0, a collocation code for singular ODEs of the first
% order, to BVPSUITE1.1, which can also handle explicit and implicit
% problems of arbitrary order, to BVPSUITE2.0, which improves on usability
% and readability of the code. The code in the package is now structured
% in a simpler modular form.
%
%   Copyright (c) July 2018 Winfried Auzinger
%                           Merlin Fallahpour
%                           Georg Kitzhofer
%                           Othmar Koch
%                           Gernot Pulverer
%                           Gustaf S??derlind
%                           Ewa B. Weinm??ller
%                           Stefan Wurm
%                           Vienna University of Technology

if ischar(problem)                 % bvp_main('api'): handles to the internal routines
    x1 = solver_api();  valx1 = [];  solution = [];
    return
end

global chainrulecoeffs;
chainrulecoeffs = higherchainrule(max(feval(problem,'orders')));
solution=[];

x1 = feval(settings,'mesh');
n = feval_problem(problem,'n');
collMethod = feval(settings,'collMethod');
collPoints = feval(settings,'collPoints');

try
    tmp_pathfoll=feval(problem,'pathfollowing');
    if exist('pathfoll','var')
        tmp_fn=fieldnames(pathfoll);
        for kk=1:numel(tmp_fn)
            tmp_pathfoll.(tmp_fn{kk})=pathfoll.(tmp_fn{kk});
        end
    end
    pathfoll=tmp_pathfoll;
catch
    pathfoll.activate=0;
end

showPlot=0;
debug=0;

if(nargin<3) || isempty(initProfile)
    try
        initProfile = feval_problem(problem,'initProfile');
    catch MException
        disp('Verwende Standardprofil');
        disp(MException);
        initProfile.initialMesh = x1;
        initProfile.initialValues = ones(n,length(x1));
    end
end

try
    initProfile.initialMesh = initProfile.x1;
    initProfile.initialValues = initProfile.valx1;
catch MException

end

if(isempty(initProfile.initialValues))
    disp('Verwende Standardprofil');
    disp(MException);
    initProfile.initialMesh = x1;
    initProfile.initialValues = ones(n,length(x1));
end

if(feval(problem,'EVP'))
    initProfile.initialValues(end+1,1)=norm(initProfile.initialValues(:,1),2)^2;
    for j=2:length(initProfile.initialMesh)
        initProfile.initialValues(end,j)=norm(initProfile.initialValues(1:end-1,1:j),2)^2;
    end
    if isfield(initProfile,'parameters')
        initProfile.parameters = [initProfile.parameters; initProfile.lambda];
    else
        initProfile.parameters = initProfile.lambda;
    end
end

%Gitter aus Settings wird auf Intervall aus bvp-file transformiert
interval = feval_problem(problem,'interval');
x1 = (interval(end)-interval(1))/(x1(end)-x1(1))*(x1+x1(1))+interval(1);

%Berechne Kollokationspunkte auf [0,1]
m = collPoints;
if (strcmp(collMethod,'user'))
    [rho] = collPoints;
else
    [rho] = getStandardCollocationPoints(collMethod,m);
end

if pathfoll.activate
    if ~strcmp(pathfoll.startat,'start')
        fprintf('\n<strong>Pathfollowing starting at %s</strong>\n',pathfoll.startat);
    else
        fprintf('\n<strong>Pathfollowing initializing ... </strong>\n');
        try
            pathfoll.initProfile=initProfile;
        catch
            fprintf('\nInitial profile not needed:')
        end
    end
    [p_save,p_exact,tur_pts]=pathfollowing(problem,settings,pathfoll,x1,rho);

    speichername=pathfoll.name;
    tmp=1;
    while exist( strcat(speichername,'_',num2str(tmp),'.mat'),'file' ) ==2
        tmp=tmp+1;
    end
    profilname=strcat(speichername,'_',num2str(tmp),'.mat');
    filename=strcat(pathfoll.dir,profilname);
    save(filename,'p_save','p_exact','tur_pts');

    fprintf('\nPathfollowing file saved as:\n%s\n',strcat(pathfoll.dir,profilname))

    fprintf('\nStarted run at lambda_p = %f.\n',p_save{2,1})
    if ~isempty(tur_pts)
        fprintf('Found <strong>TURNING POINT</strong> inbetween parameter values')
        for ii=1:size(tur_pts,2)
            fprintf('\n %i. <strong>***</strong>   %f and %f and %f  <strong>***</strong>',ii,tur_pts(1,ii),tur_pts(2,ii),tur_pts(3,ii))
        end
        fprintf('.\n')
    else
        fprintf('No turning points were encountered in this run.\n')
    end
    fprintf('Finished run after %i steps at lambda_p = %f.\n',size(p_save,2)-2,p_save{2,end-1})

    xfin=p_save{1,end-1}.x1;
    if length(x1)~=length(xfin) || max(abs(x1-xfin))>1e-12
        fprintf('Mesh adaptation was used during the run, ')
        if length(x1)~=length(xfin)
            fprintf('to augment the number of mesh points from %i to %i.\n',length(x1),length(xfin))
        else
            fprintf('to redistribute the %i mesh points fittingly.\n',length(x1))
        end
    elseif (feval(settings,'meshAdaptation'))
        fprintf('Mesh adaptation was not needed in this run.\n')
    end

    solution = p_save{1,end-1};
    solution.filename = filename;
end

%L??se lineares Problem
if (feval_problem(problem,'linear') == 1) && isempty(solution)
    dispDebug('Lineares Problem',debug);


    if(feval(settings,'meshAdaptation'))
        [solution] = meshadapt(problem,settings,x1,rho);
    elseif(feval(settings,'errorEstimate'))
        [~,~,solution] = solveLinearProblem(problem,x1,rho);
        solution = errorestimate(solution,problem,settings);
    else
        [~,~,solution] = solveLinearProblem(problem,x1,rho);
    end



elseif isempty(solution)
    dispDebug('Nicht-Lineares Problem',debug);


    if(feval(settings,'meshAdaptation'))
        [solution] = meshadapt(problem,settings,x1,rho,initProfile);
    elseif(feval(settings,'errorEstimate'))
        [~,~,solution] = solveNonLinearProblem(problem,settings,x1,rho,initProfile);
        solution = errorestimate(solution,problem,settings);
    else
        [~,~,solution] = solveNonLinearProblem(problem,settings,x1,rho,initProfile);
    end

end

interval = feval(problem,'interval');
if(interval(2) == Inf)
    if(trafo('splitting'))
        [x2,sortind] = sort(trafo('xi',solution.x1,0,1));
        solution.x1=[solution.x1(1:end-1),x2];
        solution.valx1 = [solution.valx1(1:end/2,1:(end-1)),solution.valx1((end/2+1):end,sortind)];
        [x2tau,sortind] = sort(trafo('xi',solution.x1tau,0,1));
        solution.x1tau=[solution.x1tau(1:end-1),x2tau];
        solution.valx1tau = [solution.valx1tau(1:end/2,1:(end-1)),solution.valx1tau((end/2+1):end,sortind)];
    else
        [solution.x1,sortind] = sort(trafo('xi',solution.x1,0,interval(1)));
        solution.valx1 = solution.valx1(:,sortind);
        [solution.x1tau,sortind] = sort(trafo('xi',solution.x1tau,0,interval(1)));
        solution.valx1tau = solution.valx1tau(:,sortind);
    end
end
x1 = solution.x1;
valx1 = solution.valx1;
if(feval(problem,'EVP')>0)
    valx1 = valx1(1:end-1,:);
    solution.lambda = solution.parameters(end);
    solution.parameters = solution.parameters(1:end-1);
    if isempty(solution.parameters)
        solution.parameters = [];
    end
end
if(showPlot>0)
    plot(x1,valx1);
end

end


function [solution] = meshadapt(problem,settings,x1,rho,initProfile)
% Choose the mesh adaptation routine.
% Default: the original bvpsuite2.0 routine meshadaptation.m.
% Optional plug-in: if settings(...,'meshAdaptationFcn') returns a function name or
% handle, that routine is called instead.  It must have the same interface as
% meshadaptation.m:  [solution] = fcn(problem,settings,x1,rho,initProfile).
% (Nothing in this folder refers to any particular plug-in.)
fcn = [];
try
    fcn = feval(settings,'meshAdaptationFcn');
catch
end
if nargin < 5, initProfile = []; end
if ~isempty(fcn)
    [solution] = feval(fcn,problem,settings,x1,rho,initProfile,solver_api());
elseif nargin < 5
    [solution] = meshadaptation(problem,settings,x1,rho);
else
    [solution] = meshadaptation(problem,settings,x1,rho,initProfile);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  default_settings.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret] = default_settings(request)

switch request
    case 'mesh'
        ret=0:1/50:1;
    case 'collMethod'
        ret='gauss';
    case 'collPoints'
        ret=3;
    case 'meshAdaptation'
        ret=1;
    case 'errorEstimate'
        ret=0;
    case 'absTolSolver'
        ret=1e-12;
    case 'relTolSolver'
        ret=1e-12;
    case 'absTolMeshAdaptation'
        ret=1e-9;
    case 'relTolMeshAdaptation'
        ret=1e-9;
    case 'minInitialMesh'
        ret=50;

    case 'finemesh'
        ret = 0;

    % Solver settings -- see https://repositum.tuwien.at/handle/20.500.12708/9088
    case 'allowTRM'
        ret=1;
    case 'maxFunEvalsTRM'
        ret=90000000;
    case 'maxIterationsTRM'
        ret=90000000;
    case 'lambdaMin'
        ret = 0.001;
    case 'maxAdaptations'
        ret=18;
    case 'switchToFFNFactor'
        ret=0.5;
    case 'updateJacFactor'
        ret=0.5;
    case 'K'
        ret=200;

    % Pathfollowing settings -- see https://repositum.tuwien.at/handle/20.500.12708/1167
    case 'thetaMax'
        ret=0.05;
    case 'maxCorrSteps'
        ret=5;
    case 'maxSteplengthGrowth'
        ret=4;
    case 'angleMin'
        ret=0.75;
    case 'meshFactorMax'
        ret=2;
    case 'PredLengthFactor'
        ret=1;
    case 'CorrLengthGrowth'
        ret=4;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  feval_problem.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret] = feval_problem(problem,request,varargin) %problem,request,z,za,zb,zc,t,p,lambda,EVPtrafo,lambda_p)

interval=feval(problem,'interval');
% if isempty(interval)
%     switch request
%         case 'interval'
%             ret = 0;
%             return
%         otherwise
%             ret=[0,0];
%     end
% end
EVPtrafo=1;
if(nargin>=10 && ~isempty(varargin{8}))
    EVPtrafo = varargin{8};
end
if(~strcmp(request,'n') && ~strcmp(request,'orders') && ~strcmp(request,'parameters'))
    z=zeros(feval_problem(problem,'n',[],[],[],[],[],[],[],EVPtrafo,[]),max(feval_problem(problem,'orders',[],[],[],[],[],[],[],EVPtrafo,[])));
    za=z;
    zb=z;
    zc=[];
    p=zeros(feval_problem(problem,'parameters',[],[],[],[],[],[],[],EVPtrafo,[]),1);
end
t=0;
lambda=0;
lambda_p=0;

if(nargin>=3 && ~isempty(varargin{1}))
    z = varargin{1};
end
if(nargin>=4 && ~isempty(varargin{2}))
    za = varargin{2};
end
if(nargin>=5 && ~isempty(varargin{3}))
    zb = varargin{3};
end
if(nargin>=6 && ~isempty(varargin{4}))
    zc = varargin{4};
end
if(nargin>=7 && ~isempty(varargin{5}))
    t = varargin{5};
end
if(nargin>=8 && ~isempty(varargin{6}))
    p = varargin{6};
end
if(nargin>=9 && ~isempty(varargin{7}))
    lambda = varargin{7};
end
if(nargin>=11 && ~isempty(varargin{9}))
    lambda_p = varargin{9};
end

if(interval(2)<Inf)
    if(feval(problem,'EVP')>0 && EVPtrafo)
        switch request
            case 'n'
                ret=feval(problem,request)+1;
            case 'orders'
                ret=feval(problem,request);
                ret=[ret 1];
            case 'problem'
                ret=feval(problem,request,z,[],[],[],t,p(1:end-1),p(end),lambda_p);
                ret=[ret; z(end,2)-norm(z(1:end-1,1),2)^2];
            case 'jacobian'
                ret_old=feval(problem,request,z,[],[],[],t,p(1:end-1),p(end),lambda_p);
                [s(1),s(2),s(3)]=size(ret_old);
                ret = zeros(s(1)+1,s(2)+1,max(2,s(3)));
                ret(1:s(1),1:s(2),1:s(3))=ret_old;
                ret(s(1)+1,s(2)+1,2)=1;
                for i=1:s(2)
                    ret(s(1)+1,i,1)=-2*z(i,1);
                end
            case 'interval'
                ret = feval(problem,request);
            case 'linear'
                ret = 0;
            case 'parameters'
                ret=feval(problem,request)+1;
            case 'c'
                ret = feval(problem,request);
                if(~isempty(ret))
                    ret = sort([interval(1) ret interval(2)]);
                end
            case 'BV'
                ret=feval(problem,request,[],za,zb,zc,[],p(1:end-1),p(end),lambda_p);
                ret=[ret; za(end,1); zb(end,1)-1];
            case 'dBV'
                ret_old=feval(problem,request,[],za,zb,zc,[],p(1:end-1),p(end),lambda_p);
                [s(1),s(2),s(3),s(4)]=size(ret_old);
                if(isempty(feval(problem,'c')))
                    ret = zeros(s(1),s(2)+2,s(3)+1,s(4));
                    ret(1:s(1),1:s(2),1:s(3),1:s(4))=ret_old;
                    ret(1,s(2)+1,s(3)+1,1)=1;
                    ret(s(1),s(2)+2,s(3)+1,1)=1;
                else
                    ret = zeros(s(1)+2,s(2)+2,s(3)+1,s(4));
                    ret(2:s(1)+1,1:s(2),1:s(3),1:s(4))=ret_old;
                    ret(1,s(2)+1,s(3)+1,1)=1;
                    ret(s(1)+2,s(2)+2,s(3)+1,1)=1;
                end
            case 'dP'
                ret=feval(problem,request,z,[],[],[],t,p(1:end-1),p(end),lambda_p);
                dLambda=feval(problem,'dLambda',z,[],[],[],t,p(1:end-1),p(end),lambda_p);
                ret = [ret,dLambda];
                ret = [ret;0];
            case 'dP_BV'
                ret=zeros(2+length(feval(problem,'BV',[],za,zb,zc,[],p(1:end-1),p(end),lambda_p)),1+feval(problem,'parameters'));
                ret(1:end-2,1:end-1)=feval(problem,request,[],za,zb,zc,[],p(1:end-1),p(end),lambda_p);
            case 'initProfile'
                ret=feval(problem,request);
            otherwise
                ret=feval(problem,request,z,za,zb,zc,t,p(1:end-1),p(end),lambda_p);
        end
    else
        if(strcmp(request,'n') || strcmp(request,'orders') || strcmp(request,'parameters'))
            ret = feval(problem,request);
        else
            ret = feval(problem,request,z,za,zb,zc,t,p,lambda,lambda_p);
        end
    end
else
    if(feval(problem,'EVP')>0 && EVPtrafo)
        switch request
            case 'n'
                ret=feval(problem,request)+1;
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=2*ret;
                end
            case 'orders'
                ret=feval(problem,request);
                ret=[ret 1];
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=[ret ret];
                end
            case 'problem'
                s=size(z);
                if(interval(1) == 0 && interval(2) == Inf)
                    dimz = s(1)/2-1;
                    ret=[feval(problem,request,z(1:dimz,:),[],[],[],t,p(1:end-1),p(end),lambda_p);...
                        z(dimz+1,2)-norm(z(1:dimz,1),2)^2;...
                        feval(problem,request,z((dimz+2):(s(1)-1),:)*(feval_problem_trafomatrix(1,1/t,s(2)))',[],[],[],1/t,p(1:end-1),p(end),lambda_p);...
                        z(end,2)*(-t^2)-norm(z(dimz+2:(s(1)-1),1),2)^2];
                else
                    zA=z*(feval_problem_trafomatrix(interval(1),interval(1)/t,s(2)))';
                    ret=feval(problem,request,zA,[],[],[],interval(1)/t,p(1:end-1),p(end),lambda_p);
                    ret=[ret; (-t^2/interval(1))*z(end,2)-norm(z(1:end-1,1),2)^2];
                end
            case 'jacobian'
                if(interval(1) == 0 && interval(2) == Inf)
                    s = size(z);
                    dimz = s(1)/2-1;
                    ret1=feval(problem,request,z(1:dimz,:),[],[],[],t,p(1:end-1),p(end),lambda_p);
                    ret2=feval(problem,request,z((dimz+2):(s(1)-1),:)*(feval_problem_trafomatrix(1,1/t,s(2)))',[],[],[],1/t,p(1:end-1),p(end),lambda_p);
                    [s(1),s(2),s(3)]=size(ret2);
                    A=feval_problem_trafomatrix(1,1/t,s(3));
                    for i=1:s(1)
                        ret2(i,:,:)=reshape(ret2(i,:,:),s(2),s(3))*A;
                    end
                    ret = zeros(2*s(1)+2,2*s(2)+2,s(3));

                    ret(1:s(1),1:s(2),1:s(3))=ret1; %Die urspruenglichen Gleichungen

                    ret(s(1)+1,s(2)+1,2)=1; %Die Normierungsbed auf [0,1]
                    for i=1:s(2)
                        ret(s(1)+1,i,1)=-2*z(i,1);
                    end

                    ret(s(1)+1+(1:s(1)),s(2)+1+(1:s(2)),1:s(3))=ret2; %Die transf Glgen

                    ret(2*s(1)+2,2*s(2)+2,2)=-t^2; %Die Normierungsbed transf von [1,Infty)
                    for i=s(2)+1+(1:s(2))
                        ret(2*s(1)+2,i,1)=-2*z(i,1);
                    end
                else
                    s = size(z);
                    A=feval_problem_trafomatrix(interval(1),interval(1)/t,s(2));
                    zA=z*(feval_problem_trafomatrix(interval(1),interval(1)/t,s(2)))';
                    ret_old=feval(problem,request,zA,[],[],[],interval(1)/t,p(1:end-1),p(end),lambda_p);
                    [s(1),s(2),s(3)]=size(ret_old);
                    ret = zeros(s(1)+1,s(2)+1,max(2,s(3)));
                    ret(1:s(1),1:s(2),1:s(3))=ret_old;
                    ret(s(1)+1,s(2)+1,2)=1;
                    for i=1:s(2)
                        ret(s(1)+1,i,1)=-2*zA(i,1);
                    end
                    [s(1),s(2),s(3)]=size(ret);
                    for i=1:s(1)
                        ret(i,:,:)=reshape(ret(i,:,:),s(2),s(3))*A;
                    end
                end
            case 'interval'
                ret = [0,1];
            case 'linear'
                ret = 0;
            case 'parameters'
                ret=feval(problem,request)+1;
            case 'c'
                ret = feval(problem,request);
                if(~isempty(ret))
                    ret = sort([interval(1) ret interval(2)]);
                end
            case 'BV'
                s=size(za);
                if(interval(1) == 0 && interval(2) == Inf)
                    orders = feval_problem(problem,'orders');
                    n=s(1)/2-1;
                    ret=[feval(problem,request,[],za(1:n,:),[za(n+2:(s(1)-1),1),zeros(n,s(2)-1)],[],[],p(1:end-1),p(end),lambda_p);...
                        za(n+1,1); za(end,1)-1];%zc hasn't been worked out yet
                    ret2 = zeros(sum(orders)/2,1);
                    for i=1:n+1
                        A=feval_problem_trafomatrix(1,1,orders(i));
                        ret2(sum(orders(1:i-1))+(1:orders(i)))=(zb(i,1:orders(i))-zb(i+n+1,1:orders(i))*A')';
                    end
                    ret=[ret;ret2];
                else
                    A=feval_problem_trafomatrix(interval(1),interval(1),s(2));
                    zaneu=[zb(:,1),zeros(s(1),s(2)-1)];
                    zbneu=za*A';
                    ret=feval(problem,request,z,zaneu,zbneu,[],[],p(1:end-1),p(end),lambda_p);
                    ret=[ret; zaneu(end,1); zbneu(end,1)-1];
                end
            case 'dBV'
                s=size(za);
                if(interval(1) == 0 && interval(2) == Inf)
                    orders = feval_problem(problem,'orders');
                    n=s(1)/2-1;
                    ret1=feval(problem,request,[],za(1:n,:),[za(n+2:(s(1)-1),1),zeros(n,s(2)-1)],[],[],p(1:end-1),p(end),lambda_p);
                    [s(1),s(2),s(3),s(4)]=size(ret1);
                    ret = zeros(2,s(2)+2+sum(orders)/2,2*n+2,s(4));
                    for i=1:s(2) %transform original BCs
                        ret(1,i,1:n,:) = ret1(1,i,:,:);
                        ret(1,i,(n+2):(2*n+1),1) = ret1(2,i,:,1);
                    end
                    ret(1,s(2)+1,n+1,1)=1;
                    ret(1,s(2)+2,2*(n+1),1)=1;

                    for i=1:n+1 %smoothness conditions
                        A=feval_problem_trafomatrix(1,1,orders(i));
                        for j=1:orders(i)
                            ret(2,s(2)+2+sum(orders(1:i-1))+j,i,j) = 1;
                            ret(2,s(2)+2+sum(orders(1:i-1))+j,n+1+i,:) = -A(j,:);
                        end
                    end
                else
                    A=feval_problem_trafomatrix(interval(1),interval(1),s(2));
                    ret1=feval(problem,request,[],[zb(:,1),zeros(s(1),s(2)-1)],za*A',[],[],p(1:end-1),p(end),lambda_p);
                    [s(1),s(2),s(3),s(4)]=size(ret1);
                    ret2 = zeros(s(1),s(2)+2,s(3)+1,s(4));
                    ret2(1:s(1),1:s(2),1:s(3),1:s(4))=ret1;
                    ret2(1,s(2)+1,s(3)+1,1)=1;
                    ret2(s(1),s(2)+2,s(3)+1,1)=1;
                    [s(1),s(2),s(3),s(4)]=size(ret2);
                    ret=zeros(s);
                    A=feval_problem_trafomatrix(interval(1),interval(1),s(4));
                    for i=1:s(2)
                        ret(1,i,:,:)=reshape(ret2(2,i,:,:),s(3),s(4))*A';
                    end
                    A=feval_problem_trafomatrix(interval(1),interval(2),s(4));
                    for i=1:s(2)
                        ret(2,i,:,:)=reshape(ret2(1,i,:,:),s(3),s(4))*A';
                    end
                end
            case 'dP'
                s = size(z);
                if(interval(1) == 0 && interval(2) == Inf)
                    n=s(1)/2-1;
                    ret1=[feval(problem,request,z(1:n,:),[],[],[],t,p(1:end-1),p(end),lambda_p),...
                        feval(problem,'dLambda',z(1:n,:),[],[],[],t,p(1:end-1),p(end),lambda_p)];
                    A=feval_problem_trafomatrix(1,1/t,s(2));
                    zA=z((n+2):(2*n+1),:)*A';
                    ret2=[feval(problem,request,zA,[],[],[],1/t,p(1:end-1),p(end),lambda_p),...
                        feval(problem,'dLambda',zA,[],[],[],1/t,p(1:end-1),p(end),lambda_p)];
                    ret=[ret1;0;ret2;0];
                else
                    A=feval_problem_trafomatrix(interval(1),interval(1)/t,s(2));
                    zA=z*A';
                    ret_old=feval(problem,request,zA,[],[],[],interval(1)/t,p(1:end-1),p(end),lambda_p);
                    dLambda=feval(problem,'dLambda',zA,[],[],[],interval(1)/t,p(1:end-1),p(end),lambda_p);
                    ret_old = [ret_old,dLambda];
                    ret = [ret_old;0];
                    %                      [s(1),s(2),s(3)]=size(ret);
                    %                      A=trafomatrix(interval(1),interval(1)/t,s(2));
                    %                      for i=1:s(1)
                    %                          ret(i,:,:)=reshape(ret(i,:,:),s(2),s(3))*A;
                    %                      end
                end
            case 'dP_BV'
                s=size(za);
                if(interval(1) == 0 && interval(2) == Inf)
                    n=s(1)/2-1;
                    ret1 = feval(problem,'BV',[],za(1:n,:),[za(n+2:(s(1)-1),1),zeros(n,s(2)-1)],[],[],p(1:end-1),p(end),lambda_p);
                    [sBV(1),sBV(2),sBV(3),sBV(4)]=size(ret1);
                    ret = zeros(sum(feval_problem(problem,'orders'))+1+feval(problem,'parameters'),1+feval(problem,'parameters'));
                    retold = feval(problem,request,[],za(1:n,:),[za(n+2:(s(1)-1),1),zeros(n,s(2)-1)],[],[],p(1:end-1),p(end),lambda_p);
                    [sdP(1),sdP(2)]=size(retold);
                    ret(1:sdP(1),1:sdP(2))=retold;
                else
                    A=feval_problem_trafomatrix(interval(1),interval(1),s(2));
                    ret1=feval(problem,request,[],[zb(:,1),zeros(s(1),s(2)-1)],za*A',[],[],p(1:end-1),p(end),lambda_p);
                    [sdP(1),sdP(2)]=size(ret1);
                    ret = zeros(sdP(1)+2,1+feval(problem,'parameters'));
                    ret(1:sdP(1),1:sdP(2))=ret1;
                end
            case 'initProfile'
                ret=feval(problem,request);
            otherwise
                ret=feval(problem,request,z,za,zb,zc,t,p(1:end-1),p(end),lambda_p);
        end
    else
        switch request
            case 'n'
                ret=feval(problem,request);
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=2*ret;
                end
            case 'orders'
                ret=feval(problem,request);
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=[ret ret];
                end
            case 'problem'
                s=size(z);
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=[feval(problem,request,z(1:s(1)/2,:),[],[],[],t,p,lambda,lambda_p);...
                        feval(problem,request,z((s(1)/2+1):s(1),:)*(feval_problem_trafomatrix(1,1/t,s(2)))',[],[],[],1/t,p,lambda,lambda_p)];
                else
                    ret=feval(problem,request,z*(feval_problem_trafomatrix(interval(1),interval(1)/t,s(2)))',[],[],[],interval(1)/t,p,lambda,lambda_p);
                end
            case 'jacobian'
                s = size(z);
                if(interval(1) == 0 && interval(2) == Inf)
                    ret1=feval(problem,request,z(1:s(1)/2,:),[],[],[],t,p,lambda,lambda_p);
                    ret2=feval(problem,request,z((s(1)/2+1):s(1),:)*(feval_problem_trafomatrix(1,1/t,s(2)))',[],[],[],1/t,p,lambda,lambda_p);
                    [s(1),s(2),s(3)]=size(ret2);
                    A=feval_problem_trafomatrix(1,1/t,s(3));
                    for i=1:s(1)
                        ret2(i,:,:)=reshape(ret2(i,:,:),s(2),s(3))*A;
                    end
                    ret = zeros(2*s(1),2*s(2),s(3));
                    ret(1:s(1),1:s(2),1:s(3))=ret1;
                    ret(s(1)+(1:s(1)),s(2)+(1:s(2)),1:s(3))=ret2;
                else
                    A=feval_problem_trafomatrix(interval(1),interval(1)/t,s(2));
                    ret=feval(problem,request,z*A',[],[],[],interval(1)/t,p,lambda,lambda_p);
                    [s(1),s(2),s(3)]=size(ret);
                    for i=1:s(1)
                        ret(i,:,:)=reshape(ret(i,:,:),s(2),s(3))*A;
                    end
                end
            case 'interval'
                ret = [0,1];
            case 'linear'
                ret = feval(problem,'linear');
            case 'parameters'
                ret= feval(problem,'parameters');%0;
            case 'c'
                ret = [];
            case 'BV'
                s=size(za);
                if(interval(1) == 0 && interval(2) == Inf)
                    orders = feval_problem(problem,'orders',[],[],[],[],[],[],[],EVPtrafo,lambda_p);
                    n=s(1)/2;
                    ret1=feval(problem,request,[],za(1:n,:),[za(n+1:2*n,1),zeros(n,s(2)-1)],[],[],p,lambda,lambda_p); %zc,p hasn't been worked out yet
                    ret2 = zeros(sum(orders)/2,1);
                    for i=1:n
                        A=feval_problem_trafomatrix(1,1,orders(i));
                        ret2(sum(orders(1:i-1))+(1:orders(i)))=(zb(i,1:orders(i))-zb(i+n,1:orders(i))*A')';
                    end
                    ret=[ret1;ret2];
                else
                    A=feval_problem_trafomatrix(interval(1),interval(1),s(2));
                    ret=feval(problem,request,z,[zb(:,1),zeros(s(1),s(2)-1)],za*A',[],[],p,lambda,lambda_p);
                end
            case 'dBV'
                if(isempty(feval(problem,'c')))
                    orders = feval_problem(problem,'orders',[],[],[],[],[],[],[],EVPtrafo,lambda_p);
                    s=size(za);
                    if(interval(1) == 0 && interval(2) == Inf)
                        n=s(1)/2;
                        ret1=feval(problem,request,[],za(1:n,:),[za(n+1:2*n,1),zeros(n,s(2)-1)],[],[],p,lambda,lambda_p);
                        [s(1),s(2),s(3),s(4)]=size(ret1);
                        ret = zeros(2,sum(orders),2*n,s(4));
                        for i=1:s(2) %transform original BCs
                            ret(1,i,1:n,:) = ret1(1,i,:,:);
                            ret(1,i,n+1:2*n,1) = ret1(2,i,:,1);
                        end

                        for i=1:n %smoothness conditions
                            A=feval_problem_trafomatrix(1,1,orders(i));
                            for j=1:orders(i)
                                ret(2,s(2)+sum(orders(1:i-1))+j,i,j) = 1;
                                ret(2,s(2)+sum(orders(1:i-1))+j,n+i,:) = -A(j,:);
                            end
                        end

                    else
                        A=feval_problem_trafomatrix(interval(1),interval(1),s(2));
                        ret1=feval(problem,request,[],[zb(:,1),zeros(s(1),s(2)-1)],za*A',[],[],p,lambda,lambda_p);
                        [s(1),s(2),s(3),s(4)]=size(ret1);
                        ret=zeros(s);
                        A=feval_problem_trafomatrix(interval(1),interval(1),s(4));
                        for i=1:s(2)
                            ret(1,i,:,:)=reshape(ret1(2,i,:,:),s(3),s(4))*A;
                        end
                        A=feval_problem_trafomatrix(interval(1),interval(2),s(4));
                        for i=1:s(2)
                            ret(2,i,:,:)=reshape(ret1(1,i,:,:),s(3),s(4))*A';
                        end
                    end
                end
            case 'initProfile'
                ret = feval(problem,request);
            case 'dP'
                s=size(z);
                if(interval(1) == 0 && interval(2) == Inf)
                    ret=[feval(problem,request,z(1:s(1)/2,:),[],[],[],t,p,lambda,lambda_p);...
                        feval(problem,request,z((s(1)/2+1):s(1),:)*(feval_problem_trafomatrix(1,1/t,s(2)))',[],[],[],1/t,p,lambda,lambda_p)];
                else
                    ret=feval(problem,request,z*(feval_problem_trafomatrix(interval(1),interval(1)/t,s(2)))',[],[],[],interval(1)/t,p,lambda,lambda_p);
                end
            case 'dP_BV'
                ret = zeros( 2*sum(feval(problem,'orders'))+feval(problem,'parameters'),1 );

            otherwise
                ret=feval(problem,request,z,za,zb,zc,t,p,lambda,lambda_p);
        end
    end
end

end

function [A] = feval_problem_trafomatrix(a,t,n)
n=n-1;
A=diag([1,(-1).^(1:n).*a.^(1:n)./t.^(2*(1:n))]);
if(n>1)
    A(1+(2:n),2) = a*((-1).^(2:n).*factorial(2:n)./t.^(1+(2:n)))';
    for i=4:n+1
        A(i,3:i-1) = -((A(i-1,2:i-2)*a)/t^2 + (A(i-1,3:i-1).*((i-1)+(3:i-1)))/t);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  meshadaptation.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [finalsol,toleriert,halve]=meshadaptation(problem,settings,x1,rhoColl,initProfile,predictor)

% meshadaptation ... calculate solution using adaptive meshes
%
% meshadaptation2(aTOL,rTOL,K,bvpfile,plotsol,x1,start,bvpopt,plotres,max
% iter) calculates the solution of a BVP, given in
% bvpfile using an adaptive mesh technique in order to solve the prescribed
% tolerance requirements given in aTOL (absolute tolerance) and rTOL
% (relative tolerance).
%
% K limits the maximal ratio between the shortest and the longest
% subinterval. plotsol (0,1) toggles the plot of the solution after the
% solution process. An initial mesh and an initial solution can be passed
% to the routine using x1 and start. If left empty, the mesh is determined
% from the bvpfile and the initial profile 0 is used. maxiter regulates the
% number of allowed adaptation steps.
%
% See also BVPSUITE, RUN
%

aTOL = feval(settings,'absTolMeshAdaptation');
rTOL = feval(settings,'relTolMeshAdaptation');
K = feval(settings,'K');
maxiter = feval(settings,'maxAdaptations');
finemesh = feval(settings,'finemesh');
halve=0; % Needed in case of pathfollowing, otherwise not

plotres = 0;

ordnung=feval_problem(problem,'orders');
n=length(ordnung);%number of equations
c=feval_problem(problem,'c');%boundary values inside the interval
collMethod=feval(settings,'collMethod');
collPoints=feval(settings,'collPoints');%position of the collocation points transformed to the interval [0 1]
parameter=feval_problem(problem,'parameters');
linear=feval_problem(problem,'linear');

m = length(rhoColl);

% prepare lambda_p and predictor
try
    pathfoll=feval(problem,'pathfollowing');
    lambda_p=pathfoll.start;
catch
    lambda_p=0;
end
if ~exist('predictor','var')
    predictor=[];
elseif ~isfield(predictor,'tangent')
    tmp_disp.display=predictor.display;
    predictor=tmp_disp;
    clear('tmp_disp')
else
    % The step-length might get halved, whenever N gets bigger then
    % meshFactorMax*N_orig
    % Before points are added to the mesh, the step-length will always get
    % halved first as a precautionary measure.
    meshFactorMax = max(1,predictor.meshFactorMax);
    N_orig = length(predictor.x1)-2;
end

% find initial mesh if not transmitted to meshadaptation
if length(x1) < 1
    x1 = feval(settings,'mesh');
    interval = feval_problem(problem,'interval');
    x1 = (interval(2)-interval(1))/(x1(end)-x1(1))*(x1+x1(1))+interval(1);
end
left_end = x1(1);
right_end = x1(end);

% Calculate mesh density rho from the mesh
rho=(1./diff(x1)/length(x1))';
update_mode=1; % 1 .. update rho, 2 .. update N

M=length(x1)-1; % number of intervals

N = M-1;
Ncompare = 1e8; % initialised with a high value
order_method = m;

% Variables
minimprove = 0.1; % 0.1 means at least 10 % points less to try further density adjustment
safetyN = 0; % 0.1 means at least 10 % more points than actually suggested - not accurate anymore - use safety_sigma instead
safety_sigma = 0.9; % safety factor used in the formula of N-preestimation
LPFilter_k = order_method*2.75; % controls the gain of the adjustment of the mesh density - should be coppled to the method order
res_mode = 'tcolmid'; % triggers the choice of points for the evaluation of the residual

numequ=length(ordnung);
if length(aTOL) <= numequ
    aTOL=ones(1,numequ)*aTOL(1);
end
if length(rTOL) <= numequ
    rTOL=ones(1,numequ)*rTOL(1);
end


% DEFINITIONS ==================================
% Number of INTERVALS is M
% Number of internal GRID POINTS is N = M-1
% STEP SIZE on the auxiliary grid is dxi = 1/M
% rho is a vector with M components
% ==============================================

% Main loop: solution of problem + grid refinement

densityUpdates = 0;

for j=1:maxiter
    % ***************************
    % Calculate next mesh from rho
    %*****************************

    % *************
    % SOLVE PROBLEM
    % *************

    if(feval_problem(problem,'linear')) && ( isempty(predictor) || ~isfield(predictor,'tangent') )
        [tcol,ycol,solstruct] = solveLinearProblem(problem,x1,rhoColl);
        tau = solstruct.x1;
        y = solstruct.valx1;
        coeff = solstruct.coeff;
    else
        [tcol,ycol,solstruct,halve] = solveNonLinearProblem(problem,settings,x1,rhoColl,initProfile,1-(j>1),predictor);
        if halve==0
            tau = solstruct.x1;
            y = solstruct.valx1;
            coeff = solstruct.coeff;
            if ( isempty(predictor) || ~isfield(predictor,'tangent') )
                initProfile.initialMesh = tau;
                initProfile.initialValues = y;
                initProfile.initialCoeff = solstruct.coeff;
                initProfile.parameters = solstruct.parameters;
            else
                initProfile = [solstruct.coeff', solstruct.parameters, solstruct.lambda_p ]';
                solstruct.initCoeff = initProfile;
                lambda_p = solstruct.lambda_p ;
            end
        else
            toleriert=0; finalsol=0;
            return
        end
    end

    par = solstruct.parameters;

    % **************
    % ESTIMATE ERROR
    % **************

    [sol_coarse,sol_fine,halve] = errorestimate(solstruct,problem,settings,predictor);
    if halve~=0
        toleriert=0; finalsol=0;
        return
    end
    ecol = sol_coarse.errest;
    ecoeff = sol_coarse.coeff;
    tcol2 = sol_fine.x1tau;
    ycol2 = sol_fine.valx1tau;
    tau2 = sol_fine.x1;
    y2 = sol_fine.valx1;
    ecol2 = sol_fine.errest;

    % Calculate error norms - global error ;
    %ycol_norm = max(abs(ycol),[],1)
    %err_norm = max(abs(ecol),[],1)
    ycol_norm = abs(ycol);
    err_norm = abs(ecol);

    %ycol_norm2 = max(abs(ycol2),[],1);
    %err_norm2 = max(abs(ecol2),[],1);
    ycol_norm2 = abs(ycol2);
    err_norm2 = abs(ecol2);

    qTol = zeros(size(err_norm,1),1);
    jmax = qTol;
    qTol2 = qTol;
    jmax2 = qTol;
    for tolt=1:size(err_norm,1)
        [qTol(tolt),jmax(tolt)] = max(err_norm(tolt,:) ./ (aTOL(tolt) + ycol_norm(tolt,:) .* rTOL(tolt))); % qTol is about tolfactor = errmax/tol
        [qTol2(tolt),jmax2(tolt)] = max(err_norm2(tolt,:) ./ (aTOL(tolt) + ycol_norm2(tolt,:) .* rTOL(tolt))); % qTol is about tolfactor = errmax/tol
    end

    if plotres
        disp(['Tolerance Factor: ' num2str(qTol)])
    end
    if finemesh
        if plotres
            disp(['Tolerance Factor on the fine grid: ' num2str(qTol2)])
        end
    end

    if qTol < 1
        if plotres
            disp('Tolerance satisfied')
        end
        x1=tau;
        valx1=y;
        x1tau=tcol;
        valx1tau=ycol;
        toleriert=0;
        break
    end
    if finemesh == 1
        if qTol2 < 1
            if plotres
                disp('Tolerance satisfied on the fine grid')
            end
            coeff=ecoeff;
            x1=tau2;
            valx1=y2;
            x1tau=tcol2;
            valx1tau=ycol2;
            toleriert=1;
            break
        end
    end

    %%%%%%%%%%%%%%%%%%%%;
    % Update Procedure %;
    %%%%%%%%%%%%%%%%%%%%;

    if update_mode == 2; % update the number of mesh_points
        if densityUpdates <2
            update_mode = 1;
            densityUpdates = 1;
        else
            Nold = N;
            N = ceil((1+safetyN)*N*(max(max(err_norm,[],2)' ./ (safety_sigma*(aTOL))).^(1/(order_method+1))));
            %N = max(N,ceil(1.5^(1/order_method)*Nold)); % max increase
            % of number of points
        end
    elseif update_mode == 1; % update the mesh density
        % Decide about the next step - Update N or rho ;
        Ncompare_old = Ncompare;
        Ncompare = ceil((1+safetyN)*N*(max(max(err_norm,[],2)' ./ (safety_sigma*(aTOL))).^(1/(order_method+1))));
        if plotres
            disp(['N suggested: ' num2str(Ncompare)])
        end

        if Ncompare > (1-minimprove)*Ncompare_old % check if density is finally adjusted
            update_mode = 2;
            % here we have to do the following - it can happen, that
            % the suggested number of points increases by a huge amount
            % so it would be a bad idea to update the density in this
            % last step - we have to use the older (better) density
            % here - this may happen in rare cases
            if Ncompare <= Ncompare_old || (~exist('rhold','var'))
                N = Ncompare;
                Ncompare = 1e8;
            else
                N = Ncompare_old;
                Ncompare = 1e8;
                rho=rhold;
                x1=x1old;
            end
        else
            densityUpdates = densityUpdates + 1;
        end
    end

    % For pathfollowing a maximal number of allowed mesh points is given
    if isfield(predictor,'tangent') && N>=N_orig*meshFactorMax
        fprintf('\n  Suggested number of mesh points: %i\n Number of mesh points allowed in this step: %i\n', N+2, max(N_orig,N_orig*meshFactorMax))
        toleriert=0; finalsol=0;
        halve=1;
        return
    end

    % ******************
    % CALCULATE RESIDUAL
    % ******************
    if update_mode==1
        switch res_mode
            case 'mesh'
                resmesh = x1;
                %residual=abs(equations('residual',coeff,problem,tau,resmesh));
                residual=abs(computeResidual(problem,coeff,par,tau,rhoColl,resmesh,lambda_p));
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it

                residualfinal=zeros(1,size(residual,2)-1);
                for i=1:size(residual,1)
                    residualadd=(residual(i,1:end-1)+residual(i,2:end)).*diff(x1); % this means value on the lhs and rhs times the length
                    residualfinal=residualfinal+residualadd.^2;
                end
                residualfinal=sqrt(residualfinal);
            case 'tcol'
                resmesh = tcol;

                residual=abs(computeResidual(problem,coeff,par,tau,rhoColl,resmesh,lambda_p));
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it

                residualfinal=zeros(1,length(x1)-1);
                for i=1:size(residual,1)
                    residualadd=(residual(i,1:end-1)+residual(i,2:end)).*diff(tcol);
                    for k=1:length(x1)-1
                        residualadd2(k)=sum(residualadd(k*(order_method+1)-order_method:k*(order_method+1)));
                    end
                    residualfinal=residualfinal+residualadd2.^2;
                end
                residualfinal=sqrt(residualfinal);

            case 'midpoints'
                resmesh = x1(1:end-1)+diff(x1)/2;
                residual=abs(computeResidual(problem,coeff,par,tau,rhoColl,resmesh,lambda_p));
                residualfinal=zeros(1,size(residual,2));
                for i=1:size(residual,1)
                    residualadd=(residual(i,:)).*diff(x1); % this is the value at the midpoint times the length
                    residualfinal=residualfinal+residualadd.^2;
                end
                residualfinal=sqrt(residualfinal);
            case 'tcolmid'
                resmesh = tcol(1:end-1)+diff(tcol)/2;
                residual=abs(computeResidual(problem,coeff,par,tau,rhoColl,resmesh,lambda_p));
                % to get the residual for adjusting the mesh density we need a single
                % vector, so we take the 2-norm of the integral of it

                residualfinal=zeros(1,length(x1)-1);
                for i=1:size(residual,1)
                    residualadd=residual(i,:).*diff(tcol);
                    for k=1:length(x1)-1
                        residualadd2(k)=sum(residualadd(k*(order_method+1)-order_method:k*(order_method+1)))/(x1(k+1)-x1(k));
                    end
                    residualfinal=residualfinal+residualadd2.^2;
                end
                residualfinal=sqrt(residualfinal);
        end
    end

    Nnew = N;
    Mnew = Nnew+1;

    % Display stuff
    if plotres
        disp(['N used:      ' num2str(Nnew)])
        disp(['---------------'])
        disp(['Iteration ' num2str(j)])
    end
    if update_mode == 1
        if plotres
            disp(['Updating density'])
        end
        fprintf(1,'%s%d\n','  Density Update N=',N);
    elseif update_mode == 2
        if plotres
            disp('Updating N')
        end
        fprintf(1,'%s%d\n','N Update N=',N);
    end
    rhold = rho;
    x1old = x1;
    if update_mode == 1;
        % Update density function
        rho = LPfilter(residualfinal',rho,LPFilter_k);
    end;

    % check if the grid is within K-range
    maxrho = max(rho);
    rho=max(rho,maxrho/K);

    rho = bcTDFlogV4(rho);           % Smoothing
    rho = rho*sum(1./rho)/Mnew;      % Normalize

    x1 = (x1-x1(1))/(x1(end)-x1(1));
    I = cumtrapz(x1,([rho(1); rho]+[rho; rho(end)])'/2);
    x1 = pchip(I/I(end),x1,0:1/(Mnew):1);
    rho = 1./diff(x1)';
    x1 = left_end+(x1-0)*(right_end-left_end)/(1-0); %umrechnen von [0,1] auf [a,b]

    % Adapt the predictor at the end of each iteration
    fprintf('  Maximal error estimate value: %2.3e \n',max(max(abs(sol_coarse.errest))))
    if ~isempty(predictor) && isfield(predictor,'tangent')
        % update predictor
        initP.initialMesh=x1;

        % N_tmp = length(predictor.x1)-1;
        % n_tmp = N_tmp*(sum(ordnung)+n*m);
        % coeff_a_0 = predictor.a_0(1:n_tmp);
        % param_a_0 = predictor.a_0(n_tmp+1:n_tmp+parameter);
        %
        % initP.initialValues = coeffToValues(coeff_a_0, predictor.x1,ordnung,rhoColl,x1);
        % initP.parameters = param_a_0;
        % initP = initial_coefficients(problem,x1,initP,rhoColl,0);
        % predictor.a_0 = [ initP.initialCoeff ; predictor.lambda_p_0 ]' ;
        % predictor.tangent = update_tangent(problem,predictor.a_0,x1,rhoColl,ordnung).';
        % predictor.steplength = (predictor.lambda_p_p - predictor.a_0(end)) / predictor.tangent(end);

        % update initial Profile
        initP.initialValues = coeffToValues(initProfile, predictor.x1,ordnung,rhoColl,x1);
        initP.parameters = sol_coarse.parameters;
        initP = initial_coefficients(problem,x1,initP,rhoColl,0);
        initProfile = [ initP.initialCoeff ; lambda_p ] ;

        % For the computation on the new mesh keep lambda_p fixed
        predictor.lpfix = 1;

        predictor.x1 = x1;
    end
end

finalsol.x1 = x1;
finalsol.valx1 = valx1;
finalsol.x1tau = x1tau;
finalsol.valx1tau = valx1tau;
finalsol.coeff = coeff;
finalsol.parameters = par;
finalsol.errest = ecol;
if ~isempty(predictor)
    finalsol.lambda_p =  lambda_p;
end

end


% **********************************************************
% **********************************************************
% Functions for processing the residual ********************
% **********************************************************
% **********************************************************


function tosignal = resampleV4(fromgrid,fromsignal,togrid)
% Oversampling that guarantees that tosignal remains positive
% Written by GS, TU Wien, 24 August 2006
% V4 is original version

signal = abs(fromsignal);                   % For robustness only
meansig = norm(signal,1)/length(signal);
mag = max(signal);
signal = signal/mag;                        % Normalize to maximum of 1
signal = signal + 1e-10*exp(-10*signal);      % Limiter lifts values near zero 1e-3
signal = log(signal);                       % Transform to logarithmic scale
hi = max(signal);
lo = min(signal);
%signal = pchip(fromgrid,signal,togrid); % Resample on new grid
signal = max(lo,min(hi,signal));         % Make sure amplitude doesn't grow
signal = exp(signal);                       % Now mapped to new grid and positive
signal = signal + 1e-10*exp(-10*signal);      % Limiter lifts values near zero 1e-2
signal = bcTDFlogV4(signal);                % Smoothing
signal = mag*signal/max(signal);            % Restore signal amplitude
newmean = norm(signal,1)/length(signal);    % Correct the signal's mean value
signal = signal + meansig - newmean;
signal = max(0,signal);                     % Protect against negative values
signal = signal + 1e-10*exp(-10*signal);      % Limiter lifts values near zero 1e-3
tosignal = signal;                          % Output
end

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

function grid = rhogrid(rho)
% Construct the staggered grid for the vector function rho

% Staggered grid
M = length(rho);
offset = 1/M/2;
grid = linspace(offset,1-offset,M)';
end

function Psireturn=meshadaptation_Psi(n,rho,nr)
i=length(rho);
prod=1;
for s=1:i
    if (s~=n)
        prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
    end
end
for s=1:(nr)
    prod=polyint(prod);
end
Psireturn=prod;
%tabulars for collocation points

end

function new_tangent = update_tangent(problem,initProfile,x1_2,rho,ordnung)
m=length(rho);
psi=zeros(m,max(ordnung)+m,max(ordnung));
for ord=1:max(ordnung)
    for ii=1:m
        psi(ii,1+max(ordnung)-ord:m+max(ordnung),ord)=meshadaptation_Psi(ii,rho,ord);
    end
end
psival=zeros(max(ordnung),m,m+2);
for ord=1:max(ordnung)
    for ii=1:m
        %evaluation of psi
        psival(ord,ii,1:m)=polyval(psi(ii,:,ord),rho(1:m));
        psival(ord,ii,m+1)=polyval(psi(ii,:,ord),1);
        psival(ord,ii,m+2)=polyval(psi(ii,:,ord),0);
    end
end
jac_F = functionFDF( 'DF', problem ,initProfile,x1_2,psival,psi,rho,[]);
help1=zeros(length(jac_F(:,1)),1);
help1(end)=1;
new_tangent =jac_F\help1;
new_tangent=new_tangent/sqrt(sum(new_tangent.*new_tangent));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  solveLinearProblem.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1tau,valx1tau,solStruct] = solveLinearProblem(problem,x1,rho)

m=length(rho);
n = feval_problem(problem,'n');
ordnung = feval_problem(problem,'orders');
parameter = feval_problem(problem,'parameters');
interval = feval_problem(problem,'interval');
c=feval_problem(problem,'c');
a=interval(1);
b=interval(2);
N=length(x1)-1;

[ctmp,sortind] = sort(c);
cint=zeros(length(ctmp),1);
for i=1:length(ctmp)
    j=0;
    if (ctmp(i)<a) || (ctmp(i)>b)
        err('equations_err1');
    end
    while ctmp(i)>=x1((j)+1) && j<N
        j=j+1;
    end
    cint(i)=j-1;
end
cint(sortind)=cint;

h(1:N)=x1(2:N+1)-x1(1:N);
tau(1:N,1:m)=x1(1:N).'*ones(1,m)+h(1:N).'*rho(1:m);
faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800;];

p=zeros(parameter);

psi=zeros(m,max(ordnung)+m,max(ordnung));
for ord=1:max(ordnung)
    for i=1:m
        psi(i,1+max(ordnung)-ord:m+max(ordnung),ord)=solveLinearProblem_Psi(i,rho,ord);
    end
end

psival=zeros(max(ordnung),m,m+2);
for ord=1:max(ordnung)
    for i=1:m
        %evaluation of psi
        psival(ord,i,1:m)=polyval(psi(i,:,ord),rho(1:m));
        psival(ord,i,m+1)=polyval(psi(i,:,ord),1);
        psival(ord,i,m+2)=polyval(psi(i,:,ord),0);
    end
end

%For 0-th order special definitions for psi are necessary

%if min(ordnung)==0 || (length(was)==7 && length(strfind(was,'wertabl'))>0) || (length(was)==8 && length(strfind(was,'residual'))>0)
%   for i=1:m
%s       psi0(i,:,1)=Psi(i,rho,0);
%   end
%end

for i=1:m
    psi0(i,:,1)=solveLinearProblem_Psi(i,rho,0);
end

x0=zeros((length(x1)-1)*(n*m+sum(ordnung))+parameter,1);
for j=1:(length(x1)-1)*(n*m+sum(ordnung))+parameter
    x0(j,1)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE RESIDUUM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
residuum=functionFDF( 'F', problem,zeros(length(x0),1),x1,psival,psi,rho,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPUTE SYSTEMMATRIX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
systemmatrix=functionFDF( 'DF', problem,zeros(length(x0),1),x1,psival,psi,rho,[]);

global M1
M1 = systemmatrix;

a=-systemmatrix\residuum;

j=0;
y=zeros(n,max(ordnung),N);
z=zeros(n,m,N);
for i=0:N-1
    for q=1:max(ordnung)
        for ni=1:n
            if q<=ordnung(ni)
                j=j+1;y(ni,q,i+1)=a(j);
            end
        end
    end
    z(1:n,1:m,i+1)=reshape(a(j+1:j+n*m),n,m);
    j=j+n*m;
end
p(1:parameter)=a(j+1:j+parameter).';
output=zeros(n,N+N*m+1);
for i=1:n
    stelle=0;
    for j=0:N-1
        if ordnung(i)~=0
            if j~=0
                help=solveLinearProblem_Poptimiert(0,j,i,x1((j)+1),m,x1,h,y,z,psival,ordnung(i),m+2);
            else
                help=solveLinearProblem_P(0,j,i,x1((j)+1),m,x1,h,y,z,psi,ordnung(i));
            end
        else
            help=solveLinearProblem_P(0,j,i,x1((j)+1),m,x1,h,[],z,psi0,0);
        end
        stelle=stelle+1;output(i,stelle)=help;
        for k=1:m
            if ordnung(i)~=0
                help=solveLinearProblem_Poptimiert(0,j,i,tau((j)+1,k),m,x1,h,y,z,psival,ordnung(i),k);
            else
                help=solveLinearProblem_P(0,j,i,tau((j)+1,k),m,x1,h,[],z,psi0,0);
            end
            stelle=stelle+1;output(i,stelle)=help;
        end
    end
    if ordnung(i)~=0
        help=solveLinearProblem_P(0,N-1,i,x1((N)+1),m,x1,h,y,z,psi,ordnung(i));
    else
        help=solveLinearProblem_P(0,N-1,i,x1((N)+1),m,x1,h,[],z,psi0,0);
    end
    stelle=stelle+1;
    output(i,stelle)=help;
end
valx1tau=output;

stelle=0;
x1tau=zeros(1,N);
for j=0:N-1
    stelle=stelle+1;x1tau(stelle)=x1((j)+1);
    for k=1:m
        stelle=stelle+1;x1tau(stelle)=tau((j)+1,k);
    end
end
stelle=stelle+1;
x1tau(stelle)=x1((N)+1);

solStruct.x1 = x1tau(1:m+1:end);
solStruct.valx1 = valx1tau(:,1:m+1:end);
solStruct.x1tau = x1tau;
solStruct.valx1tau = valx1tau;
solStruct.parameters = p;
solStruct.coeff = a(1:end-length(p));
end




function ret=solveLinearProblem_P(abl,i,komp,x,m,x1,h,y,z,psi,ord) %Polynomial of derivation abl

faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];
sum=0;
for j=1:ord
    if j-abl-1>=0
        sum=sum+((x-x1((i)+1))^(j-abl-1))/(faktorielle(j-abl))*y(komp,j,(i)+1);
    end
end
if ord==0
    for j=1:m
        sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1+1),(x-x1((i)+1))/h((i)+1)));
    end
else
    for j=1:m
        sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1),(x-x1((i)+1))/h((i)+1)));
    end
end
ret=sum;

end

function ret=solveLinearProblem_Poptimiert(abl,i,komp,x,m,x1,h,y,z,psival,ord,k)
%Polynomial opitmized version
faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];

ret=sum(((x-x1((i)+1)).^((abl+1:ord)-abl-1))./(faktorielle((abl+1:ord)-abl)).*y(komp,abl+1:ord,(i)+1))+...
    sum(h((i)+1)^(ord-abl).*(z(komp,1:m,(i)+1).*psival(ord-(abl+1)+1,1:m,k)));

end

function Psireturn=solveLinearProblem_Psi(n,rho,nr)
i=length(rho);
prod=1;
for s=1:i
    if (s~=n)
        prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
    end
end
for s=1:(nr)
    prod=polyint(prod);
end
Psireturn=prod;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  solveNonLinearProblem.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1tau,valx1tau,solStruct,halve] = solveNonLinearProblem(problem,settings,x1,rho,initProfile,transformInitProfile,predictor)

m=length(rho);
n = feval_problem(problem,'n');
N=length(x1)-1;
h(1:N)=x1(2:N+1)-x1(1:N);
tau(1:N,1:m)=x1(1:N).'*ones(1,m)+h(1:N).'*rho(1:m);
ordnung = feval_problem(problem,'orders');
parameter = feval_problem(problem,'parameters');

p=zeros(parameter);
y=zeros(n,max(ordnung),N);
z=zeros(n,m,N);

if(nargin<6)
    transformInitProfile = 1;
end

if ~exist('predictor','var') || isempty(predictor)
    predictor=[];
    % transformInitProfile=0;
elseif ~isfield(predictor,'tangent')
    tmp_disp.display=predictor.display;
    predictor=tmp_disp;
    clear('tmp_disp')
end

if isempty(predictor) || ~isfield(predictor,'tangent')
    initProfile=initial_coefficients(problem,x1,initProfile,rho,transformInitProfile);
    x0 = initProfile.initialCoeff;
else
    x0 = initProfile;
end
solStruct.initProfile = initProfile;

psi=zeros(m,max(ordnung)+m,max(ordnung));
for ord=1:max(ordnung)
    for i=1:m
        psi(i,1+max(ordnung)-ord:m+max(ordnung),ord)=solveNonLinearProblem_Psi(i,rho,ord);
    end
end
psival=zeros(max(ordnung),m,m+2);
for ord=1:max(ordnung)
    for i=1:m
        %evaluation of psi
        psival(ord,i,1:m)=polyval(psi(i,:,ord),rho(1:m));
        psival(ord,i,m+1)=polyval(psi(i,:,ord),1);
        psival(ord,i,m+2)=polyval(psi(i,:,ord),0);
    end
end

%For 0-th order special definitions for psi are necessary

%if min(ordnung)==0 || (length(was)==7 && length(strfind(was,'wertabl'))>0) || (length(was)==8 && length(strfind(was,'residual'))>0)
%   for i=1:m
%s       psi0(i,:,1)=Psi(i,rho,0);
%   end
%end

for i=1:m
       psi0(i,:,1)=solveNonLinearProblem_Psi(i,rho,0);
end

[coeff,halve] = solve_nonlinear_sys(@FDF,x0,problem,x1,settings,psival,psi,rho,predictor);

% Check for pathfollowing whether the tolerances have been satisfied as
% wished or the step-length needs to be halved
if halve~=0
    x1tau=0; valx1tau=0; solStruct=0;
    return
end

%-----------------------------------------------
%entries of "coeff" will be transformed to the names of the documentation
%y,z,p see manual!
j=0;
for i=0:N-1
    for q=1:max(ordnung)
        for ni=1:n
            if q<=ordnung(ni)
                j=j+1;y(ni,q,(i)+1)=coeff(j);
            end
        end
    end
    z(1:n,1:m,(i)+1)=reshape(coeff(j+1:j+n*m),n,m);
    j=j+n*m;
end
p(1:parameter)=coeff(j+1:j+parameter).';

valx1tau=zeros(n,N+N*m+1);
for i=1:n
   stelle=0;
   for j=0:N-1
       if ordnung(i)~=0
           if j~=0
               help=solveNonLinearProblem_Poptimiert(0,j,i,x1((j)+1),m,x1,h,y,z,psival,ordnung(i),m+2);
           else
               help=solveNonLinearProblem_P(0,j,i,x1((j)+1),m,x1,h,y,z,psi,ordnung(i));
           end
       else
           help=solveNonLinearProblem_P(0,j,i,x1((j)+1),m,x1,h,[],z,psi0,0);
       end
       stelle=stelle+1;valx1tau(i,stelle)=help;
       for k=1:m
           if ordnung(i)~=0
               help=solveNonLinearProblem_Poptimiert(0,j,i,tau((j)+1,k),m,x1,h,y,z,psival,ordnung(i),k);
           else
               help=solveNonLinearProblem_P(0,j,i,tau((j)+1,k),m,x1,h,[],z,psi0,0);
           end
           stelle=stelle+1;valx1tau(i,stelle)=help;
       end
   end
   if ordnung(i)~=0
       help=solveNonLinearProblem_P(0,N-1,i,x1((N)+1),m,x1,h,y,z,psi,ordnung(i));
   else
       help=solveNonLinearProblem_P(0,N-1,i,x1((N)+1),m,x1,h,[],z,psi0,0);
   end
   stelle=stelle+1;
   valx1tau(i,stelle)=help;
end

stelle=0;
x1tau=zeros(1,N);
for j=0:N-1
     stelle=stelle+1;x1tau(stelle)=x1((j)+1);
     for k=1:m
         stelle=stelle+1;x1tau(stelle)=tau((j)+1,k);
     end
end
stelle=stelle+1;x1tau(stelle)=x1((N)+1);


solStruct.x1 = x1tau(1:m+1:end);
solStruct.valx1 = valx1tau(:,1:m+1:end);
solStruct.x1tau = x1tau;
solStruct.valx1tau = valx1tau;
solStruct.parameters = p;
if isempty(predictor) || ~isfield(predictor,'tangent')
    solStruct.coeff = coeff(1:end-length(p));
else
    solStruct.coeff = coeff(1:end-length(p)-1);
    solStruct.lambda_p = coeff(end);
end
solStruct.initProfile = initProfile;

end %end of Switch

%--------Local functions

%C: help=P(0,j,i,x1((j)+1),m,x1,h,y,z,psi,ordnung(i));


function ret=solveNonLinearProblem_P(abl,i,komp,x,m,x1,h,y,z,psi,ord) %Polynomial of derivation abl

faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];
sum=0;
for j=1:ord
   if j-abl-1>=0
       sum=sum+((x-x1((i)+1))^(j-abl-1))/(faktorielle(j-abl))*y(komp,j,(i)+1);
   end
end
if ord==0
   for j=1:m
       sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1+1),(x-x1((i)+1))/h((i)+1)));
   end
else
   for j=1:m
       sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1),(x-x1((i)+1))/h((i)+1)));
   end
end
ret=sum;
end

function ret=solveNonLinearProblem_Poptimiert(abl,i,komp,x,m,x1,h,y,z,psival,ord,k)
%Polynomial opitmized version
faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];

ret=sum(((x-x1((i)+1)).^((abl+1:ord)-abl-1))./(faktorielle((abl+1:ord)-abl)).*y(komp,abl+1:ord,(i)+1))+...
   sum(h((i)+1)^(ord-abl).*(z(komp,1:m,(i)+1).*psival(ord-(abl+1)+1,1:m,k)));
end

function Psireturn=solveNonLinearProblem_Psi(n,rho,nr)
 i=length(rho);
 prod=1;
 for s=1:i
     if (s~=n)
         prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
     end
 end
 for s=1:(nr)
     prod=polyint(prod);
 end
 Psireturn=prod;
     %tabulars for collocation points

end

function [outputF,outputDF] = FDF(x0,problem,x1,psival,psi,rho,predictor)
outputF = functionFDF('F',problem,x0,x1,psival,psi,rho,predictor);
outputDF = functionFDF('DF',problem,x0,x1,psival,psi,rho,predictor);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  solve_nonlinear_sys.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [new_x,halve,logstruct,itcount,fcount] = ...
    solve_nonlinear_sys(FDFhandle,x0,problem,tau,settings,psival,psi,rho,predictor)
%   SOLVE_NONLINEAR_SYS  Solve nonlinear systems obtained by SBVPCOL and SBVPERR
%
%   This routine is private and should not be accessed by other callers than SBVPCOL and SBVPERR

% ********** store some fields of bvpopt as variables for quick reference
AbsTol      = feval(settings,'absTolSolver');
RelTol      = feval(settings,'relTolSolver');
max_F_evals = feval(settings,'maxFunEvalsTRM');
max_iter    = feval(settings,'maxIterationsTRM');
if ~isempty(predictor) && isfield(predictor,'display')
    display     =  predictor.display; % Pathfollowing display messages
else
    display     = 1; % BVP solving display messages
end
log         = 0;
halve       = 0; % Needed in case of pathfollowing, otherwise not
TRM         = feval(settings,'allowTRM');

% set 1 to use lsqnonlin for the TRM. If set to 0, fsolve is used instead
lsq = 1;


% ********** some initializations
fcount = 1;   % number of function evaluations
itcount= 1;   % number of iterations
x = x0;
lambda = 1;
lambdamin = feval(settings,'lambdaMin');
nPreviousTRMIterates = 0;
logstruct = [];
updateJacFactor = feval(settings,'updateJacFactor');
switchToFFNFactor = feval(settings,'switchToFFNFactor');

if log
    logstruct.DOC = [];
    logstruct.tol_factor = [];
    logstruct.G0 = [];
    logstruct.G = [];
    logstruct.jac_update = [];
    logstruct.lambda = [];
    logstruct.nCorrections = [];

    logstruct.runtime.DFeval = [];
    logstruct.runtime.Feval = [];
    logstruct.runtime.lu = [];
    logstruct.runtime.resubstitute = [];
end

% Structure of the Algorithm:
% The Zerofinder consists of 3 different algorithms:
%    *) Fast Frozen Newton (cheap, but small domain of convergence (DOC) G1)
%    *) Predictor-Corrector Linesearch (expensive, but large DOC G2)
%    *) Trust Region Method (even more expensive, but even larger DOC G3)
%
% For the choice of the appropriate algorithm to use, it is necessary to determine the
% position of the current approximation w.r.t. the DOCs G1, G2, G3

% ********** Determine position of initial approximation

if display && isempty(predictor) % display information every iteration step
    fprintf('\n Determining zerofinder algorithm ... \n\n');
elseif display
    fprintf('  Performing Newton step...')
end

[TolFactor, DOC, new_x, U, L, G0, G, F, delta_x, simplified_delta_x, fcount, logstruct] = ...
    determine_position2(x,problem,tau,settings,psival,psi,lambdamin,fcount,logstruct,[],[],rho,predictor);


if TolFactor < 1
    if display && isempty(predictor)
        fprintf('\n\n Tolerances satisfied\n');
    elseif display
        fprintf(' Tolerances satisfied.\n')
    end
    return
end

if DOC == 1
    x = new_x; % Accept new approximation
end

% ********** initialize Zerofinder Display
if display  % display information every iteration step
    if ~isempty(predictor)
        fprintf('\n\n')
    end
    switch DOC
        case 1
            fprintf('  Fast Frozen Newton\n');
            fprintf('  STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
        case 2
            fprintf('  Predictor-Corrector Line Search\n');
            fprintf('  STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
        case 3
            fprintf('  Trust Region Method\n');
    end
end

% ********** Choice-of-Algorithm Loop
while 1

    if log
        logstruct.DOC = [logstruct.DOC DOC];
    end

    if itcount > max_iter
        error(' Maximum number of iterations exceeded');
    end

    if fcount > max_F_evals
        error(' Maximum number of function evaluations exceeded');
    end



    switch DOC


        % **************************************************************************
        case 1 % ******************* Fast Frozen Newton ****************************
            % **************************************************************************
            % variables that have to be available: a (=x), G0, G, F (=F(x)), simplified_delta_x

            % ********** Determine if Jacobian has to be updated
            UpdateJac = G > updateJacFactor * G0;

            if log
                logstruct.lambda = [logstruct.lambda 1];
                logstruct.nCorrections = [logstruct.nCorrections 0];
                logstruct.jac_update = [logstruct.jac_update UpdateJac];
            end

            % ********** Determine Jacobian
            if UpdateJac % Update Jacobian
                cpt=cputime;
                DF = functionFDF('DF',problem,x,tau,psival,psi,rho,predictor);
                if log
                    logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
                    logstruct.condestDF = condest(DF);
                end

                cpt=cputime;
                [L,U]=lu(DF);                  % LU-decomposition of DF
                if log
                    logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt];
                end

                lastwarn('');                  % initialize warning state

                cpt=cputime;
                delta_x = U\(L\(- F));         % Newton correction
                if log
                    logstruct.runtime.resubstitute = ...
                        [logstruct.runtime.resubstitute cputime-cpt];
                end

                G0 = norm(delta_x);            % Norm of the Newton correction
            else % Iterate with frozen Jacobian
                delta_x = simplified_delta_x;
                G0 = G; % = norm(simplified_delta_x)
            end

            % ********** Perform Newton step
            new_x(:) = x(:) + delta_x;      % new approximation

            cpt = cputime;
            new_F = functionFDF('F',problem,new_x,tau,psival,psi,rho,predictor);   % new residual
            if log
                logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
            end

            fcount = fcount +1;             % increase number of function evaluations

            % ********** Determine Position of new solution approximation
            cpt = cputime;
            simplified_delta_x = U\(L\(-new_F));            % simplified delta_x
            if log
                logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
            end

            G = norm(simplified_delta_x);

            if log
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G =  [logstruct.G G];
            end

            if G < (1-lambdamin/2) * G0 % Approximation improved

                % ***** Check if Tolerances are satisfied
                % ***** TolFactor: Factor by which the correction had to be
                % ***** reduced in order to satisfy the tolerances
                TolFactor = check_tolerances(new_x,simplified_delta_x,AbsTol,RelTol);

                if log
                    logstruct.tol_factor = [logstruct.tol_factor TolFactor];
                end

                if TolFactor < 1
                    % Perform the last step towards the solution
                    new_x(:) = new_x(:) + simplified_delta_x;

                    if display  % display information every iteration step
                        fprintf('  %3i   %5i    %8i     %13.2e    %.2e   %5.2e\n',...
                            itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);
                    end

                    if display && isempty(predictor) % off
                        fprintf('\n\n Tolerances satisfied\n');
                    elseif display
                        fprintf('\n  Tolerances satisfied.\n')
                    end

                    return
                end

                DOC = 1;           % Keep iterating with the Fast Frozen Newton

                x = new_x;         % Accept new approximation
                F = new_F;

                if display  % display information every iteration step
                    fprintf('  %3i   %5i    %8i     %13.2e    %.2e   %5.2e\n',...
                        itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);
                end
            elseif UpdateJac == 0 % Jacobian has not been updated in the previous step
                DOC = 1;           % -> Try FFN once again,
                %    but update Jacobian this time
                %    (We would have to evaluate the Jacobian
                %    anyway in the Line Search Algorithm)

                if display  % display information every iteration step
                    fprintf('  %3i   %5i    %8i     %13.2e    %.2e   %5.2e   NOT ACCEPTED\n',...
                        itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);
                end
            else
                DOC = 2;           % -> Switch to Predictor Corrector Line Search

                % ********** Prepare for next step with the respective method
                lambda = 1;

                % Use Jacobian computed here
                if log
                    logstruct.jac_update  = [logstruct.jac_update  0];
                end

                if display
                    fprintf('  %3i   %5i    %8i     %13.2e    %.2e   %5.2e   NOT ACCEPTED\n',...
                        itcount,fcount,UpdateJac,G/norm(x(:)),G0/G,TolFactor);

                    fprintf('\n  Switching to Predictor-Corrector Line Search\n');
                    fprintf('  STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                end
            end

            itcount = itcount + 1;


            % **************************************************************************
        case 2 % *************** Predictor-Corrector Line Search *******************
            % **************************************************************************
            % Necessary variables at this point: a (=x), delta_x, simplified_delta_x,
            % lambda (=lambda_pred), G0, G (=G(lambda_pred))

            % Let G(lambda) = ||DF^(-1) * F(x + lambda*delta_x)||
            % At this point, x, delta_x, lambda, G(0), G'(0), G(lambda) are known
            % We check if lambda is accepted and correct lambda otherwise

            Accept_Lambda = G < (1-lambdamin/2)*G0;
            nCorrections = 0;

            while ~Accept_Lambda
                % ********** Repeatedly correct lambda until it is accepted or a termination condition is met

                if nCorrections == 0 % First correction, quadratical method
                    % We know G(0), G'(0) and G=G(lambda) by previous computations.
                    % To obtain a correction, we interpolate a quadratic polynomial through
                    % these values and define the minimum of this polynomial as the corrected lambda
                    % The respective polynomial reads G(x)=G(0) + x * G'(0) + x^2/lambda^2 *
                    % (G(lambda) - G(0) - lambda * G'(0)). Zeroing the first derivative yields:

                    % store lambda (=lambda_pred) for possible further correction
                    l2 = lambda;

                    Gprime0 = -G0;
                    lambda = - lambda^2 * Gprime0 / (2* (G - G0 - lambda * Gprime0));

                    % Allow lambda_cor between lambda_pred / 10 and 1
                    lambda = max(l2/10 , min(lambda, 1));

                    %***** Prepare for possible further correction (cubical, match notation with Num. Recipes)
                    l1 = lambda;

                    nCorrections = 1;
                else % Subsequent corrections, cubical method
                    coeff = 1/(l1-l2) * [1/l1^2 -1/l2^2 ; -l2/l1^2 l1/l2^2] * ...
                        [G1 - Gprime0*l1 - G0 ; G2 - Gprime0*l2 - G0];

                    ca = coeff(1);
                    cb = coeff(2);

                    lambda = (-cb + sqrt(cb^2-3*ca*Gprime0))/(3*ca);

                    % Allow lambda_cor between lambda_cor_old / 10 and lambda_cor_old /2
                    lambda = max(l1/10 , min(lambda, l1/2));

                    % Prepare for possible further correction
                    l2 = l1;
                    l1 = lambda;

                    nCorrections = nCorrections + 1;
                end

                if lambda < lambdamin
                    DOC = 3;  % switch to Trust Region Method
                    break;
                end

                new_x(:) = x(:) + lambda * delta_x;  % get new approximation

                cpt = cputime;
                F = functionFDF('F',problem,new_x,tau,psival,psi,rho,predictor);   % evaluate F at lambda_cor
                if log
                    logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
                end

                fcount = fcount + 1;

                cpt = cputime;
                simplified_delta_x = U\(L\-F);
                if log
                    logstruct.runtime.resubstitute = ...
                        [logstruct.runtime.resubstitute cputime-cpt];
                end

                G_old = G; % store old G for possible further correction
                G = norm(simplified_delta_x);

                Accept_Lambda = G < (1-lambdamin/2) * G0;

                % ***** Prepare for possible further correction
                G1 = G;
                G2 = G_old;
            end

            if Accept_Lambda % a lambda has been accepted
                x = new_x;    % Accept new solution

                if G < switchToFFNFactor * G0
                    DOC = 1;
                else
                    DOC = 2;
                end
            end

            itcount = itcount + 1;

            if log
                logstruct.lambda = [logstruct.lambda lambda];
                logstruct.nCorrections = [logstruct.nCorrections nCorrections];
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G =  [logstruct.G G];
            end

            % At this point, DOC has been determined. Check Termination conditions if DOC ~= 3
            if DOC ~= 3
                TolFactor = check_tolerances(x,simplified_delta_x,AbsTol,RelTol);

                if log
                    logstruct.tol_factor = [logstruct.tol_factor TolFactor];
                end

                if TolFactor < 1
                    new_x(:) = x(:) + simplified_delta_x;

                    if display  % display information every iteration step
                        fprintf('  %3i   %5i        %5.4f    %12.2e    %.2e   %5.2e\n',...
                            itcount,fcount,lambda,G/norm(x(:)),G0/G,TolFactor);
                    end

                    if display && isempty(predictor) % off
                        fprintf('\n\n Tolerances satisfied\n');
                    elseif display
                        fprintf('\n  Tolerances satisfied.\n')
                    end

                    return
                end
            elseif log
                logstruct.tol_factor = [logstruct.tol_factor -1];
            end




            if display  % display information every iteration step
                fprintf('  %3i   %5i        %5.4f    %12.2e    %.2e   %5.2e\n',...
                    itcount,fcount,lambda,G/norm(x(:)),G0/G,TolFactor);

                if DOC == 1
                    fprintf('\n  Switching to Fast Frozen Newton\n');
                    fprintf('  STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                elseif DOC == 3
                    fprintf('\n  Switching to Trust Region Method\n');
                end
            end


            % If we keep up doing the Line Search, we have to predict lambda for the next step
            if DOC == 2

                % ********** Keep old values of delta_x for prediction
                G0_old = G0;

                % ********** Determine delta_x
                cpt = cputime;
                DF = functionFDF('DF',problem,x,tau,psival,psi,rho,predictor);
                if log
                    logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
                    logstruct.jac_update = [logstruct.jac_update 1];
                    logstruct.condestDF = condest(DF);
                end

                cpt=cputime;
                [L,U]=lu(DF);                  % LU-decomposition of DF
                if log
                    logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt];
                end

                lastwarn('');                  % initialize warning state

                cpt=cputime;
                delta_x = U\(L\(- F));         % Newton correction
                if log
                    logstruct.runtime.resubstitute = ...
                        [logstruct.runtime.resubstitute cputime-cpt];
                end

                G0 = norm(delta_x);

                % ********** Predict lambda [Deuflhardt74]
                %dpr%      warning off MATLAB:divideByZero  % If Jacobians are the same, delta_x = simplified_delta_x
                mu = lambda * G0_old / norm(simplified_delta_x - delta_x);
                %dpr%      warning on MATLAB:divideByZero

                if mu > 0.7
                    lambda = 1;
                else
                    lambda = mu;
                end

                % ********** Make sure all necessary variables are available for correction step
                new_x(:) = x(:) + lambda * delta_x;

                cpt = cputime;
                F = functionFDF('F',problem,new_x,tau,psival,psi,rho,predictor);
                if log
                    logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
                end

                fcount = fcount + 1;

                cpt = cputime;
                simplified_delta_x = U\(L\-F);
                if log
                    logstruct.runtime.resubstitute = ...
                        [logstruct.runtime.resubstitute cputime-cpt];
                end

                G = norm(simplified_delta_x);
            end


            % **************************************************************************
        case 3 % ********************* Trust Region Method *************************
            % **************************************************************************

            if isfield(predictor,'tangent')
                fprintf('  Trust Region Method is turned off during pathfollowing.\n')
                halve=3;
                return
            end

            if(~TRM)
                disp('tried to use TRM. Solution may be imprecise or wrong');
                return;
            end


            % Necessary variables at this point: a (=x), nPreviousTRMIterates
            % ********** Determine accuracy for Trust Region Method
            options = optimset('Display','off','MaxIter',max_iter,'MaxFunEvals',max_F_evals);

            options.TolX = sqrt(max(RelTol,AbsTol));
            options.TolFun = 0;        % Assemble necessary options for FSOLVE
            options.Jacobian = 'on';   % (Display, MaxIter, MaxFunEvals are user parameters)
            options.LargeScale = 'on';
            switch nPreviousTRMIterates
                case 0
                    options.TolX = sqrt(max(RelTol,AbsTol));
                    options.TolFun = 0;        % Assemble necessary options for FSOLVE
                    options.Jacobian = 'on';   % (Display, MaxIter, MaxFunEvals are user parameters)
                    options.LargeScale = 'on';
                case 1
                    options.TolX = max(RelTol,AbsTol);
                case 2
                    options.TolX = 1000*eps;
                case 3
                    options.Algorithm =  'levenberg-marquardt';
                case 4
                    fprintf('Zerofinder cannot fulfill the given tolerances');
                    return;
            end


            options.Display = 'iter';
            % *********** Limit the maximum number of iterations and FunEvals for TRM
            options.MaxIter = max_iter - itcount;
            options.MaxFunEvals = max_F_evals - fcount;

            % ********** Invoke Trust Region Method
            if lsq
                [x,~,F,ExitFlag,lsqLog,~,DF] = ...
                    lsqnonlin(FDFhandle,x,[],[],options,problem,tau,psival,psi,rho,predictor);
            else
                [x,F,ExitFlag,lsqLog,DF] = ...
                    fsolve(FDFhandle,x,options,problem,tau,psival,psi,rho,predictor);
            end

            nPreviousTRMIterates = nPreviousTRMIterates +1;
            fcount = fcount + lsqLog.funcCount;

            % ***** Determine whether anything has gone wrong
            if nPreviousTRMIterates==3
                switch ExitFlag
                    case 0
                        fprintf('Number of iterations or number of function evaluations exeed the given borders');
                        return;
                    case -1
                        fprintf('Zerofinding algorithm did not converge');
                        return;
                    case -2
                        fprintf('Zerofinding algorithm converged against a point that is no solution');
                        return;
                    case -3
                        fprintf('Zerofinding algorithm did not converge');
                        return;
                    case -4
                        fprintf('Zerofinding algorithm did not converge');
                        return;
                end
            end



            % ********** Determine position of new approximation
            [TolFactor, DOC, new_x, U, L, G0, G, F, delta_x, simplified_delta_x, fcount, logstruct] = ...
                determine_position2(x,problem,tau,settings,psival,psi,lambdamin,fcount,logstruct,F,DF,rho,predictor);

            if TolFactor < 1
                if display % off
                    fprintf('\n\n Tolerances satisfied\n');
                end
                return
            end

            % ***** We only log those quantities if we continue iterating
            if log
                logstruct.lambda = [logstruct.lambda -1];
                logstruct.nCorrections = [logstruct.nCorrections -1];
                logstruct.jac_update = [logstruct.jac_update -1];
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G =  [logstruct.G G];
                logstruct.tol_factor = [logstruct.tol_factor TolFactor];
            end


            if DOC == 1
                x = new_x; % Accept new approximation
            elseif DOC == 2  % use Jacobian computed here
                if log
                    logstruct.jac_update  = [logstruct.jac_update  0];
                end %dpr%
            end


            if display  % display information every iteration step
                switch DOC
                    case 1
                        fprintf('  Switching to Fast Frozen Newton\n');
                        fprintf('  STEP   F-COUNT   JAC_UPDATE   |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                    case 2
                        fprintf('  Switching to Predictor-Corrector Line Search\n');
                        fprintf('  STEP   F-COUNT     LAMBDA     |DELTA_X|/|X|   IMP_FAC   TOLDIST_FAC \n');
                    case 3
                        fprintf('  Continuing with Trust Region Method\n');
                end
            end
    end
end
end



% *************************************************************************
% ********** DETERMINE POSITION OF INITIAL APPROXIMATION ******************
% *************************************************************************

function [TolFactor, DOC, new_x, U, L, G0, G, F, delta_x, simplified_delta_x, fcount, logstruct] = ...
    determine_position2(x,problem,tau,settings,psival,psi,lambdamin,fcount,logstruct,F,DF,rho,predictor)

log = 0;
AbsTol      = feval(settings,'absTolSolver');
RelTol      = feval(settings,'relTolSolver');
switchToFFNFactor = feval(settings,'switchToFFNFactor');

% ********** Determine Newton correction
% ***** Evaluate Jabobian, if not provided (by fsolve)
if isempty(F)
    cpt = cputime;
    DF=functionFDF('DF',problem,x,tau,psival,psi,rho,predictor);
    if log
        logstruct.runtime.DFeval = [logstruct.runtime.DFeval cputime-cpt];
        logstruct.condestDF = condest(DF);
    end

    cpt = cputime;
    F = functionFDF('F',problem,x,tau,psival,psi,rho,predictor);   % new residual
    if log
        logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
    end

    fcount = fcount +1;            % increase number of function evaluations
end

cpt=cputime;
[L,U]=lu(DF);                  % LU-decomposition of DF
if log
    logstruct.runtime.lu = [logstruct.runtime.lu cputime-cpt];
end

lastwarn('');                  % initialize warning state

cpt=cputime;

delta_x = U\(L\(- F));         % Newton correction
if log
    logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
end


G0 = norm(delta_x);            % Norm of the Newton correction

new_x = zeros(size(x));
new_x(:) = x(:) + delta_x;     % new approximation

% ***** In case we have very good initial approximation (e.g. from a previous mesh),
% ***** tolerances might be satisfied even now and we reduce computational
% ***** costs to the absolute minimum of one F- and one DF-evaluation
TolFactor = check_tolerances(new_x,delta_x,AbsTol,RelTol);
if TolFactor < 1
    DOC = [];  % dummy outputs
    G = [];
    simplified_delta_x = [];

    if log
        logstruct.tol_factor = [logstruct.tol_factor TolFactor];
        logstruct.G0 = [logstruct.G0 G0];
    end

    return
end

% ********** Check Monotonicity Condition for lambda = 1 and lambda = lambda_min
cpt = cputime;
F = functionFDF('F',problem,new_x,tau,psival,psi,rho,predictor);   % new residual
if log
    logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
end

fcount = fcount +1;             % increase number of function evaluations

cpt = cputime;
simplified_delta_x = U\(L\(-F));            % simplified delta_x
if log
    logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
end

G = norm(simplified_delta_x);

if G <(1-lambdamin/2) * G0 % Approximation improved
    if G < switchToFFNFactor * G0         % Approximation improved significantly
        DOC = 1;             % -> Fast Frozen Newton

        % See if tolerances are satisfied
        TolFactor = check_tolerances(new_x,simplified_delta_x,AbsTol,RelTol);
        if TolFactor < 1
            % Perform one last step towards the solution
            new_x(:) = new_x(:) + simplified_delta_x;

            if log
                logstruct.tol_factor = [logstruct.tol_factor TolFactor];
                logstruct.G0 = [logstruct.G0 G0];
                logstruct.G = [logstruct.G G];
            end

            return
        end
    else % Approximation improved slightly
        DOC = 2;             % -> Predictor Corrector Line Search
        if log
            logstruct.jac_update  = [logstruct.jac_update  0];   % Use Jacobian computed here
        end
    end
else % try another Newton step with smallest possible lambda

    new_x(:) = x(:) + lambdamin * delta_x;      % new approximation

    cpt = cputime;
    Fmin = functionFDF('F',problem,new_x,tau,psival,psi,rho,predictor);   % new residual
    if log
        logstruct.runtime.Feval = [logstruct.runtime.Feval cputime-cpt];
    end

    fcount = fcount +1;

    cpt = cputime;
    simplified_delta_x_min = U\(L\(-Fmin));
    if log
        logstruct.runtime.resubstitute = [logstruct.runtime.resubstitute cputime-cpt];
    end

    Gmin = norm(simplified_delta_x_min);

    if Gmin < (1-lambdamin/2) * G0 % Approximation improved
        DOC = 2; % Predictor Corrector Line Search
        if log
            logstruct.jac_update  = [logstruct.jac_update  0];   % Use Jacobian computed here
        end
    else
        DOC = 3; % -> Trust Region Method
    end
end
end


%********************************************************************************
%******************************* CHECK TOLERANCES *******************************
%********************************************************************************

function TolFactor = check_tolerances(x,delta,AbsTol,RelTol)
% Calculates the factor by which delta has to be reduced to satisfy the
% tolerances.
% out <= 1   => Tolerances satisfied
% out >  1   => Tolerances not satisfied

% ***** Determine Dimension of solution array x
[d , ~, N] = size(x);

% ***** Extract solution components at mesh points (last meshpoint (b) not included
InfNorm_y = max(abs(reshape(x(:,1,:),d,N)));

% ***** Extract last correction at mesh points
delta = reshape(delta,size(x));
InfNorm_delta = max(abs(reshape(delta(:,1,:),d,N)));

% ***** Componentwise check if tolerances are satisfied
TolFactor = max(InfNorm_delta ./ (AbsTol + InfNorm_y * RelTol));

if TolFactor<1
    pause(1)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  functionFDF.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ output ] = functionFDF( request, problem,a,x1,psival,psi,rho,predictor)
% request can be either 'F' to compute F(a)
%                    or 'DF' to compute DF(a)

m=length(rho);
n = feval_problem(problem,'n');
ordnung = feval_problem(problem,'orders');
parameter = feval_problem(problem,'parameters');
interval = feval_problem(problem,'interval');
interval_feval = feval(problem,'interval'); % needed for DF in case of pathfollowing
a1=interval(1);
b=interval(2);
c=feval_problem(problem,'c');
N=length(x1)-1;

[ctmp,sortind] = sort(c);
cint=zeros(length(ctmp),1);
for i=1:length(ctmp)
    j=0;
    if (ctmp(i)<a1) || (ctmp(i)>b)
        err('equations_err1');
    end
    while ctmp(i)>=x1((j)+1) && j<N
        j=j+1;
    end
    cint(i)=j-1;
end
cint(sortind)=cint;

h(1:N)=x1(2:N+1)-x1(1:N);
tau(1:N,1:m)=x1(1:N).'*ones(1,m)+h(1:N).'*rho(1:m);
faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800];

%-----------------------------------------------
%entries of "a" will be transformed to the names of the documentation
%y,z,p see Runge-Kutta Basis in the manual!
j=0;
y=zeros(n,max(ordnung),N);
z=zeros(n,m,N);
for i=0:N-1
    for q=1:max(ordnung)
        for ni=1:n
            if q<=ordnung(ni)
                j=j+1; y(ni,q,(i)+1)=a(j);
            end
        end
    end
    z(1:n,1:m,(i)+1)=reshape(a(j+1:j+n*m),n,m);
    j=j+n*m;
end
p(1:parameter)=a(j+1:j+parameter).';

n_tmp = N*(sum(ordnung)+n*m)+parameter;
lambda = [];
if length(a)>n_tmp
    lambda_p=a(end);
    if feval(problem,'EVP')
        lambda = a(end-1);
    end
else
    try
        pathfoll=feval(problem,'pathfollowing');
        lambda_p=pathfoll.start;
    catch
        lambda_p=0;
    end
end

switch request
    % ***********************************************************************************
    case 'F' % *************************** Computation of F(X) ************************************
        % ******************************************************************************
        output=zeros(length(a),1);
        j=0;
        if isempty(c)
            fza=zeros(n,max(ordnung));
            fzb=zeros(n,max(ordnung));
            for oi=1:max(ordnung)
                for ni=1:n
                    if oi<=ordnung(ni)
                        fza(ni,oi)=y(ni,oi,(0)+1);
                        fzb(ni,oi)=functionFDF_Poptimiert(oi-1,N-1,ni,b,m,x1,h,y,z,psival,ordnung(ni),m+1);
                    end
                end
            end
            Rx=feval_problem(problem,'BV',[],fza,fzb,[],[],p,lambda,[],lambda_p);
            if sum(ordnung)+parameter>0
                output(1+j:sum(ordnung)+parameter+j,1)=Rx(1:sum(ordnung)+parameter);
                j=j+sum(ordnung)+parameter;
            end
        else
            fzc=zeros(n,max(ordnung),length(c));
            for oi=1:max(ordnung)
                for ni=1:n
                    if oi<=ordnung(ni)
                        for ci=1:length(c)
                            fzc(ni,oi,ci)=functionFDF_P(oi-1,cint(ci),ni,c(ci),m,x1,h,y,z,psi,ordnung(ni));
                        end
                    end
                end
            end
            Rx=feval_problem(problem,'BV',[],[],[],fzc,[],p,lambda,[],lambda_p);

            if sum(ordnung)+parameter>0
                output(1+j:sum(ordnung)+parameter+j,1)=Rx(1:sum(ordnung)+parameter);
                j=j+sum(ordnung)+parameter;
            end
        end

        %-----------------------------------------------
        %Evaluation of the equations for the solver
        fz=zeros(n,max(ordnung)+1);
        for i=0:N-1
            for k=1:m
                for oi=0:max(ordnung)-1
                    for ni=1:n
                        if oi<ordnung(ni)
                            fz(ni,oi+1)=functionFDF_Poptimiert(oi,i,ni,tau((i)+1,k),m,x1,h,y,z,psival,ordnung(ni),k); % psi wird immer von der Ordnung der entsprechenden Gleichung gewaehlt
                        end
                    end
                end
                for ni=1:n
                    fz(ni,ordnung(ni)+1)=z(ni,k,(i)+1);
                end
                if n>0
                    g_=feval_problem(problem,'problem',fz,[],[],[],tau((i)+1,k),p,lambda,[],lambda_p);
                    output(j+1:j+n,1)=g_(1:n);
                    j=j+n;
                end
            end
            if (i<N-1)
                for oi=0:max(ordnung)-1
                    for ni=1:n
                        if oi<=ordnung(ni)-1
                            j=j+1;u=functionFDF_Poptimiert(oi,i,ni,x1((i+1)+1),m,x1,h,y,z,psival,ordnung(ni),m+1)-y(ni,oi+1,(i+1)+1);
                            output(j,1)=u;
                        end
                    end
                end
            end
        end

        % In case of pathfollowing
        if length(a)>n_tmp && ~isempty(predictor)
            if ~predictor.lpfix
                output(N*(sum(ordnung)+n*m)+parameter+1,1)=(a.'-predictor.a_0)*predictor.tangent.'-predictor.steplength;
            else
                output(N*(sum(ordnung)+n*m)+parameter+1,1)=0;
            end
        end

        % ***********************************************************************************
    case 'DF' % *************************** Computation of DF(X) ************************************
        % ******************************************************************************

        %----Start Boundary conditions
        %R_0
        j=0;
        % Dimensioning of non-zero entries
        % Boundary conditions + Matrix
        nonzeroent=(sum(ordnung)+parameter)*(max(ordnung)*n*max(2,length(c))+m*n+parameter)+...
            N*m*n*max(ordnung)*n+...
            N*m*n*m*n+...
            N*max(ordnung)*n*max(ordnung)*n+...
            N*max(ordnung)*n*m+...
            N*max(ordnung)*n*max(ordnung)*n;
        ROW=zeros(nonzeroent,1);
        COL=zeros(nonzeroent,1);
        OUTPUT=zeros(nonzeroent,1);
        gen_c=0; %General count for sparse matrix efficiency
        if isempty(c)
            fza=zeros(n,max(ordnung));
            fzb=zeros(n,max(ordnung));
            for oi=1:max(ordnung)
                for ni=1:n
                    if oi<=ordnung(ni)
                        fza(ni,oi)=y(ni,oi,(0)+1);
                        fzb(ni,oi)=functionFDF_Poptimiert(oi-1,N-1,ni,b,m,x1,h,y,z,psival,ordnung(ni),m+1);
                    end
                end
            end
            % p1 index for boundary conditions of the equations, oi index for the
            % derivatives
            %p2 the "respect variables" of the differentiation have more than one
            %component
            DpR=feval_problem(problem,'dP_BV',[],fza,fzb,[],[],p,lambda,[],lambda_p);
            dBV=feval_problem(problem,'dBV',[],fza,fzb,[],[],p,lambda,[],lambda_p);
            [size_dBV(1),size_dBV(2),size_dBV(3),size_dBV(4)]=size(dBV);
            for p1=1:sum(ordnung)+parameter
                DRa_=zeros(max(ordnung),size_dBV(3));
                DRb_=zeros(max(ordnung),size_dBV(3));
                for k=1:max(ordnung)
                    temp = reshape(dBV(1,p1,:,k),size_dBV(3),1);
                    for l=1:numel(temp)
                        DRa_(k,l)=temp(l);
                    end
                    temp = reshape(dBV(2,p1,:,k),size_dBV(3),1);
                    for l=1:numel(temp)
                        DRb_(k,l)=temp(l);
                    end
                end

                if (N~=1)
                    j=j+1;
                    zaehler=0;
                    zaehler2=0;
                    %the derivatives with respect to y_0k follow
                    for k=1:max(ordnung)
                        help = dBV(1,p1,:,k);
                        for p2=1:n
                            if k<=ordnung(p2)
                                zaehler=zaehler+1;
                                gen_c=gen_c+1;
                                ROW(gen_c)=j;
                                COL(gen_c)=zaehler;
                                OUTPUT(gen_c)=help(p2);
                            end
                        end
                    end
                    %the derivatives with respect to y_(N-1)k follow
                    for k=1:max(ordnung)
                        for p2=1:n
                            if k<=ordnung(p2)
                                zaehler2=zaehler2+1;
                                help=0;
                                for abl=0:k-1
                                    help=help+DRb_(abl+1,p2)*((b-x1((N-1)+1))^(k-abl-1)/faktorielle(k-abl));
                                end
                                if (k-abl-1)>=0
                                    col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
                                    gen_c=gen_c+1;
                                    ROW(gen_c)=j;
                                    COL(gen_c)=col;
                                    OUTPUT(gen_c)=help;
                                else
                                    col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
                                    gen_c=gen_c+1;
                                    ROW(gen_c)=j;
                                    COL(gen_c)=col;
                                    OUTPUT(gen_c)=0;
                                end
                            end
                        end
                    end
                else
                    j=j+1;
                    zaehler2=0;
                    for k=1:max(ordnung)
                        for p2=1:n
                            if k<=ordnung(p2)
                                zaehler2=zaehler2+1;
                                help=0;
                                for abl=0:k-1
                                    help=help+DRb_(abl+1,p2)*((b-a1)^(k-abl-1)/faktorielle(k-abl));
                                end
                                if (k-abl-1)>=0
                                    gen_c=gen_c+1;
                                    ROW(gen_c)=j;
                                    COL(gen_c)=zaehler2;
                                    OUTPUT(gen_c)=DRa_(k,p2)+help;
                                else
                                    gen_c=gen_c+1;
                                    ROW(gen_c)=j;
                                    COL(gen_c)=zaehler2;
                                    OUTPUT(gen_c)=DRa_(k,p2);
                                end
                            end
                        end
                    end
                end
                %derivatives with respect to z_(N-1)q
                for q=1:m
                    for p2=1:n
                        zaehler2=zaehler2+1;
                        help=0;
                        for abl=0:ordnung(p2)-1
                            help=help+DRb_(abl+1,p2)*h((N-1)+1)^(ordnung(p2)-abl)*polyval(psi(q,:,ordnung(p2)-(abl+1)+1),1);
                        end
                        col=(N-1)*(sum(ordnung)+m*n)+zaehler2;
                        gen_c=gen_c+1;
                        ROW(gen_c)=j;
                        COL(gen_c)=col;
                        OUTPUT(gen_c)=help;
                    end
                end
                %derivatives with respect to the unknown parameters p_1 to
                %p_parameter
                if parameter>0
                    col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
                    gen_c_start=gen_c+1;
                    gen_c=gen_c+parameter+1;%length(col);
                    ROW(gen_c_start:gen_c)=j;
                    COL(gen_c_start:gen_c)=col;
                    OUTPUT(gen_c_start:gen_c)=DpR(p1,1:parameter);
                end
            end
        else   % for additional conditions
            fzc=zeros(n,max(ordnung),length(c));
            for oi=1:max(ordnung)
                for ni=1:n
                    if oi<=ordnung(ni)
                        for ci=1:length(c)
                            fzc(ni,oi,ci)=functionFDF_P(oi-1,cint(ci),ni,c(ci),m,x1,h,y,z,psi,ordnung(ni));
                        end
                    end
                end
            end
            % p1 index for boundary conditions of the equations, oi index for the
            % derivatives
            %p2 the "respect variables" of the differentiation have more than one
            %component
            %derivatives with respect to z_(N-1)q
            DpR=feval_problem(problem,'dP_BV',[],[],[],fzc,[],p,lambda,[],lambda_p);
            dBV=feval_problem(problem,'dBV',[],[],[],fzc,[],p,lambda,[],lambda_p);
            for p1=1:sum(ordnung)+parameter
                DR_=zeros(max(ordnung),n,length(c));
                for k=1:max(ordnung)
                    for ci=1:length(c)
                        temp = dBV(ci,p1,:,k);
                        for l=1:length(temp)
                            DR_(k,l,ci)=temp(l);
                        end
                    end
                end
                % Derivation with respect to y_c1_k, z_c1_k, y_c2_k, z_c2_k, ...
                % -----------------------------------------------------------
                j=j+1;
                for ci=1:length(c)
                    zaehler2=0;
                    for k=1:max(ordnung)
                        for p2=1:n
                            if k<=ordnung(p2)
                                zaehler2=zaehler2+1;
                                help2=0;
                                % Every polynomial containing the variable y_ci has
                                % to be differentiated and added
                                for cii=1:length(c)
                                    if cint(cii)==cint(ci)
                                        help=0;
                                        for abl=0:k-1
                                            help=help+DR_(abl+1,p2,cii)*((c(cii)-x1((  cint(cii)  )+1))^(k-abl-1)/faktorielle(k-abl));
                                        end
                                        help2=help2+help;
                                    end
                                end
                                col=cint(ci)*(sum(ordnung)+m*n)+zaehler2;
                                gen_c=gen_c+1;
                                ROW(gen_c)=j;
                                COL(gen_c)=col;
                                OUTPUT(gen_c)=help2;
                            end
                        end
                    end
                    %Derivation with respect to z_(N-1)q
                    for q=1:m
                        for p2=1:n
                            zaehler2=zaehler2+1;
                            help2=0;
                            for cii=1:length(c)
                                if cint(cii)==cint(ci)
                                    help=0;
                                    for abl=0:ordnung(p2)-1
                                        help=help+DR_(abl+1,p2,cii)*h((cint(cii))+1)^(ordnung(p2)-abl)*polyval(psi(q,:,ordnung(p2)-(abl+1)+1),(c(cii)-x1((  cint(cii)  )+1))/h( cint(cii) +1));
                                    end
                                    help2=help2+help;
                                end
                            end
                            col=cint(ci)*(sum(ordnung)+m*n)+zaehler2;
                            gen_c=gen_c+1;
                            ROW(gen_c)=j;
                            COL(gen_c)=col;
                            OUTPUT(gen_c)=help2;
                        end
                    end
                end

                %-----------------------------------------------------------
                %Derivatives with respect to unknown parameters p_1 bis p_parameter
                if parameter>0
                    col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
                    gen_c_start=gen_c+1;
                    gen_c=gen_c+length(col);
                    ROW(gen_c_start:gen_c)=j;
                    COL(gen_c_start:gen_c)=col;
                    OUTPUT(gen_c_start:gen_c)=DpR(p1,1:parameter);
                end
            end
        end

        %----------End Boundary Conditions

        %----------Start J_i with C_i (notation according to manual)
        fz=zeros(n,max(ordnung)+1);
        for i=0:N-1
            for k=1:m
                for oi=0:max(ordnung)-1
                    for ni=1:n
                        if oi<ordnung(ni)
                            fz(ni,oi+1)=functionFDF_Poptimiert(oi,i,ni,tau((i)+1,k),m,x1,h,y,z,psival,ordnung(ni),k);
                        end
                    end
                end
                for ni=1:n
                    fz(ni,ordnung(ni)+1)=z(ni,k,(i)+1);
                end
                %p1 indicates, which component of g_ik will be differentiated with respect
                %to y und z
                %        for oi=1:max(ordnung)+1
                %            % Dg_(:,:,oi)=feval_problem(bf,'Dg',oi,[],tau((i)+1,k),fz,[],[],[],[],p,lambda,[],lambda_p);
                %            temp = feval_problem(bf,'jacobian',fz,[],[],[],tau((i)+1,k),p,lambda,[],lambda_p);
                %            temp=temp(:,:,oi);
                %            for ind1=1:size(temp,2)
                %                for ind2=1:size(temp,1)
                %                    Dg_(ind1,ind2,oi)=temp(ind2,ind1);
                %                end
                %            end
                %        end %%% Da war glaub ich ein Fehler, deshalb gleich:
                Dg_=feval_problem(problem,'jacobian',fz,[],[],[],tau((i)+1,k),p,lambda,[],lambda_p);
                % Dpg=feval_problem(bf,'Dpg',[],[],tau((i)+1,k),fz,[],[],[],[],p,lambda,[],lambda_p);
                Dpg=feval_problem(problem,'dP',fz,[],[],[],tau((i)+1,k),p,lambda,[],lambda_p);
                for p1=1:n
                    j=j+1;
                    %Differentiation with respect to y_i, sequence: y_i1 to
                    %y_i(max(ordnung))
                    zaehler=0;
                    for oi2=1:max(ordnung)
                        %p2 indicates, which component of y_i(oi2) will be
                        %differentiated
                        for p2=1:n
                            %Differentiation with respect to y_i(oi2) only if
                            %maximum order of the equation equals oi2
                            if oi2<=ordnung(p2)
                                zaehler=zaehler+1;
                                help=0;
                                for r1=1:oi2
                                    help=help+Dg_(p1,p2,r1)*(((tau((i)+1,k)-x1((i)+1))^(oi2-r1))/(faktorielle(oi2-r1+1)));
                                end
                                col=i*m*n+sum(ordnung)*i+zaehler;
                                gen_c=gen_c+1;
                                ROW(gen_c)=j;
                                COL(gen_c)=col;
                                OUTPUT(gen_c)=help;
                            end
                        end
                    end
                    zaehler=zaehler+1;
                    for p2=1:n
                        help(p2,1:m)=0;
                        for r1=1:ordnung(p2)
                            help(p2,1:m)=help(p2,1:m)+Dg_(p1,p2,r1)*h((i)+1)^(ordnung(p2)-r1+1).*psival(ordnung(p2)-r1+1,1:m,k);
                        end
                        help(p2,1:m)=help(p2,1:m)+Dg_(p1,p2,ordnung(p2)+1).*(k==(1:m));
                    end
                    help2=help(:).';
                    col=i*m*n+sum(ordnung)*i+zaehler:i*m*n+sum(ordnung)*i+zaehler+n*m-1;
                    gen_c_start=gen_c+1;
                    gen_c=gen_c+length(col);
                    ROW(gen_c_start:gen_c)=j;
                    COL(gen_c_start:gen_c)=col;
                    OUTPUT(gen_c_start:gen_c)=help2;
                    if parameter>0
                        col=N*(n*m+sum(ordnung))+1:N*(n*m+sum(ordnung))+parameter;
                        gen_c_start=gen_c+1;
                        gen_c=gen_c+length(col);
                        ROW(gen_c_start:gen_c)=j;
                        COL(gen_c_start:gen_c)=col;
                        OUTPUT(gen_c_start:gen_c)=Dpg(p1,1:parameter);
                    end
                end
            end
            if i<N-1
                for abl=0:max(ordnung)-1
                    for p1=1:n
                        if abl<=ordnung(p1)-1
                            j=j+1;
                            %Differentiation with respect to y_i(oi2)
                            zaehler=0;
                            for oi2=1:max(ordnung)
                                for p2=1:n
                                    if oi2<=ordnung(p2)
                                        zaehler=zaehler+1;
                                        if p1==p2
                                            if oi2-abl-1>=0
                                                col=i*(m*n+sum(ordnung))+zaehler;
                                                gen_c=gen_c+1;
                                                ROW(gen_c)=j;
                                                COL(gen_c)=col;
                                                OUTPUT(gen_c)=((x1((i+1)+1)-x1((i)+1))^(oi2-abl-1))/(faktorielle(oi2-abl));
                                            else
                                                col=i*(m*n+sum(ordnung))+zaehler;
                                                gen_c=gen_c+1;
                                                ROW(gen_c)=j;
                                                COL(gen_c)=col;
                                                OUTPUT(gen_c)=0;
                                            end
                                        end
                                    end
                                end
                            end
                            %D with respect to z_iq
                            zaehler=0;
                            for q=1:m
                                zaehler=zaehler+1;
                                col=i*m*n+(i+1)*sum(ordnung)+(zaehler-1)*n+p1;
                                gen_c=gen_c+1;
                                ROW(gen_c)=j;
                                COL(gen_c)=col;
                                OUTPUT(gen_c)=h((i)+1)^(ordnung(p1)-abl)*psival(ordnung(p1)-(abl+1)+1,q,m+1);
                            end
                            zaehler=0;
                            for oi2=1:max(ordnung)
                                for p2=1:n
                                    if oi2<=ordnung(p2)
                                        zaehler=zaehler+1;
                                        if p1==p2
                                            col=(i+1)*m*n+(i+1)*sum(ordnung)+zaehler;
                                            gen_c=gen_c+1;
                                            ROW(gen_c)=j;
                                            COL(gen_c)=col;
                                            OUTPUT(gen_c)=-(oi2==abl+1);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        dimsparse=N*(n*m+sum(ordnung))+parameter;
        ROW=ROW(1:gen_c);
        COL=COL(1:gen_c);
        OUTPUT=OUTPUT(1:gen_c);
        output=sparse(ROW,COL,OUTPUT,dimsparse,dimsparse);

        % In case of pathfollowing
        if length(a)>n_tmp

            if isempty(c)
                fza=zeros(n,max(ordnung));
                fzb=zeros(n,max(ordnung));
                for oi=1:max(ordnung)
                    for ni=1:n
                        if oi<=ordnung(ni)
                            fza(ni,oi)=y(ni,oi,(0)+1);
                            fzb(ni,oi)=functionFDF_Poptimiert(oi-1,N-1,ni,b,m,x1,h,y,z,psival,ordnung(ni),m+1);
                        end
                    end
                end
                % If the problem is posed on [0,Inf)
                if interval_feval(2) == Inf && interval_feval(1) == 0
                    ret1 = feval(problem,'path_dBV',[],fza,fzb,[],[],p,lambda,lambda_p);
                    ret2 = zeros(sum(ordnung)/2,1);
                    p_dB = [ret1 ; ret2];
                else
                    p_dB=feval(problem,'path_dBV',[],fza,fzb,[],[],p,lambda,lambda_p);
                end
                % If the problem is an Eigenvalue problem
                if ~isempty(lambda)
                    p_dB = [p_dB ; 0; 0];
                end
            else   % for additional conditions
                fzc=zeros(n,max(ordnung),length(c));
                for oi=1:max(ordnung)
                    for ni=1:n
                        if oi<=ordnung(ni)
                            for ci=1:length(c)
                                fzc(ni,oi,ci)=functionFDF_P(oi-1,cint(ci),ni,c(ci),m,x1,h,y,z,psi,ordnung(ni));
                            end
                        end
                    end
                end
                % If the problem is posed on [0,Inf)
                if interval_feval(2) == Inf && interval_feval(1) == 0
                    ret1 = feval(problem,'path_dBV',[],[],[],fzc,[],p,lambda,lambda_p);
                    ret2 = zeros(sum(ordnung)/2,1);
                    p_dB = [ret1 ; ret2];
                else
                    p_dB=feval(problem,'path_dBV',[],[],[],fzc,[],p,lambda,lambda_p);
                end
                % If the problem is an Eigenvalue problem
                if ~isempty(lambda)
                    p_dB = [p_dB ; 0; 0];
                end
            end
            output(1:sum(ordnung)+parameter,N*(sum(ordnung)+n*m)+parameter+1)=p_dB;

            j=sum(ordnung)+parameter+1;

            fz=zeros(n,max(ordnung)+1);
            for i=0:N-1
                for k=1:m
                    for oi=0:max(ordnung)-1
                        for ni=1:n
                            if oi<ordnung(ni)
                                fz(ni,oi+1)=functionFDF_Poptimiert(oi,i,ni,tau((i)+1,k),m,x1,h,y,z,psival,ordnung(ni),k);
                            end
                        end
                    end
                    for ni=1:n
                        fz(ni,ordnung(ni)+1)=z(ni,k,(i)+1);
                    end
                    s=size(fz);
                    if interval_feval(2) == Inf && interval_feval(1) == 0
                        ret1 = feval(problem,'path_jac',fz(1:s(1)/2,:),[],[],[],tau((i)+1,k),p,lambda,lambda_p);
                        ret2 = feval(problem,'path_jac',fz((s(1)/2+1):s(1),:)*(functionFDF_trafomatrix(1,1/tau((i)+1,k),s(2)))',[],[],[],1/tau((i)+1,k),p,lambda,lambda_p);
                        stmp=size(ret2,1);
                        A=functionFDF_trafomatrix(1,1/tau((i)+1,k),1);
                        for iter=1:stmp
                            ret2(iter)=reshape(ret2(iter,:,:),1,1)*A;
                        end
                        p_jac = zeros(2*stmp,1);
                        p_jac(1:stmp,1)=ret1;
                        p_jac(stmp+(1:stmp),1)=ret2;
                        %p_jac = [ret1 ; ret2];
                    else
                        p_jac=feval(problem,'path_jac',fz,[],[],[],tau((i)+1,k),p,lambda,lambda_p);
                        % If the problem is an Eigenvalue problem
                        if ~isempty(lambda)
                            p_jac = [p_jac ; 0];
                        end
                    end
                    output(j:j+n-1,N*(sum(ordnung)+n*m)+parameter+1)=p_jac;
                    j=j+n;
                end
                j=j+sum(ordnung);
            end

            output(n_tmp+1,n_tmp+1) = 1;
            if ~isempty(predictor) && ~predictor.lpfix
                output(N*(sum(ordnung)+n*m)+parameter+1,:)=predictor.tangent;
            end
        end

end

end

% Local functions

function ret=functionFDF_P(abl,i,komp,x,m,x1,h,y,z,psi,ord) %Polynomial of derivation abl

faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];
sum=0;
for j=1:ord
    if j-abl-1>=0
        sum=sum+((x-x1((i)+1))^(j-abl-1))/(faktorielle(j-abl))*y(komp,j,(i)+1);
    end
end
if ord==0
    for j=1:m
        sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1+1),(x-x1((i)+1))/h((i)+1)));
    end
else
    for j=1:m
        sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1),(x-x1((i)+1))/h((i)+1)));
    end
end
ret=sum;
end

function ret=functionFDF_Poptimiert(abl,i,komp,x,m,x1,h,y,z,psival,ord,k) %Polynomial opitmized version
faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];

ret=sum(((x-x1((i)+1)).^((abl+1:ord)-abl-1))./(faktorielle((abl+1:ord)-abl)).*y(komp,abl+1:ord,(i)+1))+...
    sum(h((i)+1)^(ord-abl).*(z(komp,1:m,(i)+1).*psival(ord-(abl+1)+1,1:m,k)));
end

function [A] = functionFDF_trafomatrix(a,t,n)
n=n-1;
A=diag([1,(-1).^(1:n).*a.^(1:n)./t.^(2*(1:n))]);
if(n>1)
    A(1+(2:n),2) = a*((-1).^(2:n).*factorial(2:n)./t.^(1+(2:n)))';
    for i=4:n+1
        A(i,3:i-1) = -((A(i-1,2:i-2)*a)/t^2 + (A(i-1,3:i-1).*((i-1)+(3:i-1)))/t);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  errorestimate.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sol_coarse,sol_fine,halve] = errorestimate(sol_coarse,problem,settings,predictor)

%Definiere Werte der zu schaetzenden Loesung
N=length(sol_coarse.x1)-1;

% effort bookkeeping (optional, used by the benchmark scripts): one call = one coarse solve (N intervals)
% plus one fine solve (2N intervals)
global BVPSTATS
if isstruct(BVPSTATS)
    BVPSTATS.iters = BVPSTATS.iters + 1;
    BVPSTATS.work  = BVPSTATS.work + 3*N;
    BVPSTATS.Nhist(end+1) = N;
end
% optional per-iteration hook: settings(...,'iterationHook') may return a function handle
% hook(sol_coarse,problem,settings), called once per mesh-adaptation iteration
iterHook = [];
try
    iterHook = feval(settings,'iterationHook');
catch
end
if ~isempty(iterHook)
    feval(iterHook,sol_coarse,problem,settings);
end
%b=x1((N)+1);
%a1=x1((0)+1);
ordnung=feval_problem(problem,'orders');
n=length(ordnung);%number of equations
collMethod=feval(settings,'collMethod');
collPoints=feval(settings,'collPoints');%position of the collocation points transformed to the interval [0 1]
parameter=feval_problem(problem,'parameters');
linear=feval_problem(problem,'linear');
halve=0; % needed for pathfollowing

% prepare predictor
if ~exist('predictor','var')
    predictor=[];
end
if ~isempty(predictor) && predictor.display
    fprintf('  <strong>Error estimation...</strong>\n')
end
if isempty(predictor)
    predictor = [];
elseif ~isfield(predictor,'tangent')
    tmp_disp.display=predictor.display;
    predictor=tmp_disp;
    clear('tmp_disp')
end

m = collPoints;
if (strcmp(collMethod,'user'))
    [rho] = collPoints;
    m=length(rho);
else
    [rho] = getStandardCollocationPoints(collMethod,m);
end

x1tau_coarse = sol_coarse.x1tau;
valx1tau_coarse = sol_coarse.valx1tau;
coeff_coarse = sol_coarse.coeff;

x1 = x1tau_coarse(1:m+1:end);
valx1 = valx1tau_coarse(:,1:m+1:end);

% Add a mesh point in-between each two mesh points
x1_2(1:2:2*length(x1)) = x1;
x1_2(2:2:2*length(x1)-1) = (x1(2:length(x1))+x1(1:length(x1)-1))/2;

if(linear) && ~isfield(predictor,'tangent')
    [~,~,sol_fine] = solveLinearProblem(problem,x1_2,rho);
else
    try
        sol_coarse.lambda_p;
        predictor.a_0;

        initP.initialMesh=x1_2;

        % n_tmp = N*(sum(ordnung)+n*m);
        % coeff_a_0 = predictor.a_0(1:n_tmp);
        % param_a_0 = predictor.a_0(n_tmp+1:n_tmp+parameter);
        % coeff_tang = predictor.tangent(1:n_tmp);
        % param_tang = predictor.tangent(n_tmp+1:n_tmp+parameter);
        %
        % initP.initialValues = coeffToValues(coeff_a_0, predictor.x1,ordnung,rho,x1_2);
        % initP.parameters = param_a_0;
        % initP = initial_coefficients(problem,x1_2,initP,rho,0);
        % predictor.a_0 = [ initP.initialCoeff ; predictor.lambda_p_0 ]' ;
        %
        % initP.initialValues = coeffToValues(coeff_tang, predictor.x1,ordnung,rho,x1_2);
        % initP.parameters = param_tang;
        % initP = initial_coefficients(problem,x1_2,initP,rho,0);
        % predictor.tangent = [ initP.initialCoeff ; 0 ]';

        initP.initialValues = coeffToValues(sol_coarse.coeff, sol_coarse.x1,ordnung,rho,x1_2);
        initP.parameters = sol_coarse.parameters;
        initP = initial_coefficients(problem,x1_2,initP,rho,0);
        initProfile = [ initP.initialCoeff ; sol_coarse.lambda_p ] ;

        % For the computation on the fine mesh keep lambda_p fixed
        predictor.lpfix = 1;
    catch
        initProfile.initialMesh = x1;
        initProfile.initialValues = valx1;
        initProfile.initialCoeff = coeff_coarse;
        initProfile.parameters=sol_coarse.parameters;
    end
    [~,~,sol_fine,halve] = solveNonLinearProblem(problem,settings,x1_2,rho,initProfile,0,predictor);
    if halve~=0
        sol_coarse=0; sol_fine=0;
        return
    end
end

sol_coarse.x1tau = x1tau_coarse;
sol_coarse.valx1tau = valx1tau_coarse;
valx1tau2_2 = coeffToValues(sol_fine.coeff,sol_fine.x1,ordnung,rho,sol_coarse.x1tau); %fine solution on coarse grid + collpoints
sol_coarse.errest = (2^m / (1-2^m)) * (valx1tau2_2 - sol_coarse.valx1tau);

valx1tau = coeffToValues(sol_coarse.coeff,sol_coarse.x1,ordnung,rho,sol_fine.x1tau); %coarse solution on fine grid + collpoints
sol_fine.errest = (1 / (1-2^m)) * (sol_fine.valx1tau - valx1tau);

if isstruct(BVPSTATS)
    BVPSTATS.errhist(end+1) = max(abs(sol_coarse.errest(:)));
end

end

function Psireturn=errorestimate_Psi(n,rho,nr)
i=length(rho);
prod=1;
for s=1:i
    if (s~=n)
        prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
    end
end
for s=1:(nr)
    prod=polyint(prod);
end
Psireturn=prod;
%tabulars for collocation points
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  initial_coefficients.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ initProfile ] = initial_coefficients( problem,x1,initProfile,rho,transform )
% Compute polynomial coefficients for the polynomial interpolating the
% values in initProfile

m=length(rho);
n=feval_problem(problem,'n'); %number of solution components
ordnung=feval_problem(problem,'orders');
k=n;

if(nargin<5)
    transform = 1;
end

if(transform == 0)
    interval=feval_problem(problem,'interval');
else
    interval=feval(problem,'interval');
end

if(trafo('splitting') && interval(1)==0 && interval(2)==Inf)
    k=n/2;
end

N=length(x1)-1;
initialmesh = initProfile.initialMesh;
initialvalues = initProfile.initialValues;

parameter=feval_problem(problem,'parameters');

faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800;];

if(parameter>0)
    par = initProfile.parameters;
else
    par=[];
end

if max(ordnung)>0
    for ord=1:max(ordnung)
        for i=1:m
            psi(i,1+max(ordnung)-ord:m+max(ordnung),ord)=initial_coefficients_Psi(i,rho,ord);
        end
    end
end
if min(ordnung)==0
    for i=1:m
        psi0(i,:,1)=initial_coefficients_Psi(i,rho,0);
    end
end
for i=0:N-1
    h(i+1)=x1(i+2)-x1(i+1);
end


for ord=1:max(ordnung)
    proj_intervallstellen=0:1/(m+ord-1):1;
    for i=1:m
        psival(ord,i,1:m+ord)=polyval(psi(i,:,ord),proj_intervallstellen);
    end
end

if ~isempty(initialvalues)

    %removal of entries that appear twice in "initialmesh"
    jj=0;
    for ii=2:length(initialmesh)
        if initialmesh(ii)~=initialmesh(ii-1)
            jj=jj+1;
            initialmesh2(jj)=initialmesh(ii-1);
            initialvalues2(1:k,jj)=initialvalues(1:k,jj);
        end
    end

    initialmesh2(jj+1)=initialmesh(end);

    initialvalues2(1:k,jj+1)=initialvalues(1:k,end);
end

for p=1:n
    ordp=ordnung(p);
    interpolationsstellen=zeros(1,(N-1)*(ordp+m-1)+1);

    if(interval(2)==Inf)
        if(trafo('splitting'))
            if p>k %components on [1,infty)

                ind_1=find(initialmesh2>=1);

                if length(ind_1)<=1
                    error('initial2');
                end

                %transformation to the interval [0,1]

                %********************************************************
                %define a computational grid according to the collocation
                %rules
                %********************************************************
                for i=0:N-1
                    %constructs grid (=mesh + coll.points)
                    %inserts coll. points between elements of x1
                    interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
                end

                %****************
                %Interpolation
                %***************

                if ind_1(1)-1 <=0
                    error('initial1');
                end

                if abs(initialmesh2(ind_1(1)-1)-0)< 10^-4
                    initialmesh2(ind_1(1)-1)=10^-4;
                end

                interpolationsstellen=[interpolationsstellen,1./initialmesh2(ind_1(1)-1)];
                interpolationsstellenwerte=interp1([1./initialmesh2(ind_1(1)-1),1./initialmesh2(ind_1)],[initialvalues2(p-k,ind_1(1)-1),initialvalues2(p-k,ind_1)],interpolationsstellen,'linear','extrap');

            else
                %if the number of the solution component is <=n/2
                %then no transformation has to be carried out

                %
                ind_1=find(initialmesh2<=1);

                if length(ind_1)<=1
                    error('initial1');
                end

                for i=0:N-1
                    %construction of the computational grid (=mesh + coll.points)
                    %inserts coll. points between elements of x1
                    interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
                end

                if ind_1(end)+1 > length(initialmesh2)
                    error('initial2');
                end

                interpolationsstellen=[interpolationsstellen,initialmesh2(ind_1(end)+1)];
                interpolationsstellenwerte=interp1([initialmesh2(ind_1),initialmesh2(ind_1(end)+1)],[initialvalues2(p,ind_1),initialvalues2(p,ind_1(end)+1)],interpolationsstellen,'spline');

                %figure
                %plot(interpolationsstellen,interpolationsstellenwerte,'r')

            end
        else
            %the given problem is posed on [a,infty) a ~= 0

            ind_1=find(initialmesh2>=1);

            %transformation to the interval [0,1]

            %********************************************************
            %define a computational grid according to the collocation
            %rules
            %********************************************************
            for i=0:N-1
                %constructs grid (=mesh + coll.points)
                %inserts coll. points between elements of x1
                interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
            end
            %****************
            %Interpolation
            %***************
            interpolationsstellenwerte=interp1(trafo('xi_inv',initialmesh2(ind_1),0,interval(1)),initialvalues2(p,ind_1),interpolationsstellen,'linear','extrap');
        end

    else

        for i=0:N-1

            % x1 ... new mesh
            %constructs grid (=mesh + coll.points)
            %the number of added collocation points varies with the order
            %of derivative of the solution component z_i
            interpolationsstellen(i*(m+ordp-1)+1:(i+1)*(m+ordnung(p)-1)+1)=x1(i+1):h(i+1)/(m+ordp-1):x1(i+2);
        end
        initialvalues2(p,:);
        interpolationsstellenwerte=interp1(initialmesh2,initialvalues2(p,:),interpolationsstellen,'spline');
    end


    %***********************************************
    %calculates coefficients in Runge-Kutta basis!!
    %***********************************************


    for i=0:N-1
        if (m==1 && ordnung(p)==0)
            intervallstellen(1)=x1((i+1));
        else
            intervallstellen=x1(i+1):h(i+1)/(m+ordnung(p)-1):x1(i+2);
        end
        % if length(koeff)==0
        %intervallwerte=interp1(stellen,werte(p,:),intervallstellen,'spline');
        intervallwerte=interpolationsstellenwerte(i*(m+ordnung(p)-1)+1:(i+1)*(m+ordnung(p)-1)+1);
        %  else
        %      intervallwerte=equations('wert',koeff,bvpfile,stellen,intervallstellen);
        %      intervallwerte=intervallwerte(p,:);
        %  end

        A=zeros(m+ordnung(p),ordnung(p)+m);
        for j=1:m+ordnung(p)
            for q=1:ordnung(p)
                A(j,q)=(intervallstellen(j)-x1((i)+1))^(q-1)/faktorielle(q);
            end

            for sp=1:m
                if ordp>0
                    %A(j,sp+ordnung(p))=h((i)+1)^ordnung(p)*polyval(psi(sp,:,ordnung(p)),(intervallstellen(j)-x1(i+1))/h(i+1));
                    A(j,sp+ordnung(p))=h((i)+1)^ordnung(p)*psival(ordnung(p),sp,j);
                else
                    A(j,sp)= polyval(psi0(sp,:),(intervallstellen(j)-x1((i)+1))/h((i)+1));
                end
            end
        end
        help=A\intervallwerte';
        intervallstart((i)+1,1:length(help),p)=help;
    end
end

j=0;
coeff=zeros(N*max(ordnung)*n+N*m*n+parameter,1);
for i=0:N-1
    for k=1:max(ordnung)
        for p=1:n
            if k<=ordnung(p)
                j=j+1;
                coeff(j,1)=intervallstart((i)+1,k,p);
            end
        end
    end
    for k=1:m
        for p=1:n
            j=j+1;
            coeff(j,1)=intervallstart((i)+1,k+ordnung(p),p);
        end
    end
end
for pii=1:parameter
    j=j+1;coeff(j,1)=par(pii);
end
initProfile.initialCoeff=coeff(1:j,1);


end

function Psireturn=initial_coefficients_Psi(n,rho,nr)
i=length(rho);
prod=[1];
for s=1:i
    if (s~=n)
        prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
    end
end
for s=1:(nr)
    prod=polyint(prod);
end
Psireturn=prod;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  coeffToValues.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ z ] = coeffToValues( coeff, tau, l, rho, teval, d)
%COEFFTOVALUES Evalutes from the coefficients computed by bvpsuite values
%at teval
%   INPUT:  coeff       ... coefficient vector returned by bvpsuite
%           tau         ... grid on which the solution was computed
%           l           ... vector containing the orders of the components
%                           (i.e. [2 3], if the highest order of z1 is 2
%                           and of z2 is 3)
%           rho         ... collocation points given on [0,1]
%           teval       ... points at which solution should be evaluated
%                           (default: tau)
%           d           ... computes the d-th derivative of solution
%                           (default: 0)
%   OUTPUT: z           ... the values corresponding to teval
m = length(rho);
N=length(l);
if(nargin < 5)
    teval = tau;
end
if(nargin < 6)
    d = 0;
end

%Compute coefficents of RK-basis-polynomials in monom-basis
PHI = zeros(max(l),max(l));
for i=1:max(l)
    tmp = coeffToValues_phi(i);
    PHI(i,max(l)-length(tmp)+1:max(l))=tmp;
end

PSI = zeros(m,max(l)+1,m+max(l)+1);
for k=1:m
    for i=0:max(l)
        tmp=coeffToValues_psi(rho,k,i);
        PSI(k,i+1,end-length(tmp)+1:end)=tmp;
    end
end

z=zeros(N,length(teval));
for tpind = 1:length(teval)
    tp = teval(tpind);
    i=min(length(tau)-1,find(tp>=tau,1,'last'));
    for k=d+1:max(l)
        for j=1:N
            if(k<=l(j))
                z(j,tpind) = z(j,tpind) + coeff((sum(l)+m*N)*(i-1)+sum(min(l,k-1))+sum(l(1:j)>=k))*polyval(PHI(k-d,:),tp-tau(i));
            end
        end
    end
    for k=1:m
        for j=1:N
            if(l(j)-d+1>0)
                z(j,tpind) = z(j,tpind) + coeff((sum(l)+m*N)*(i-1)+sum(l)+N*(k-1)+j) * polyval(reshape(PSI(k,l(j)-d+1,:),m+max(l)+1,1),(tp-tau(i))/(tau(i+1)-tau(i))) * (tau(i+1)-tau(i))^(l(j)-d);
            else
                z(j,tpind) = z(j,tpind) + coeff((sum(l)+m*N)*(i-1)+sum(l)+N*(k-1)+j) * polyval(coeffToValues_psi(rho,k,l(j)-d+1),(tp-tau(i))/(tau(i+1)-tau(i))) * (tau(i+1)-tau(i))^(l(j)-d);
            end
        end
    end
end

end


function poly = coeffToValues_phi(k)
poly = [1/factorial(k-1),zeros(1,k-1)];
end

function poly = coeffToValues_psi(rho,k,order)
poly = 1;
for j=1:length(rho)
    if(j~=k)
        poly = conv(poly,[1/(rho(k)-rho(j)),-rho(j)/(rho(k)-rho(j))]);
    end
end
for j=1:order
    poly = polyint(poly);
end
for j=1:(-order)
    poly=polyder(poly);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  computeResidual.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ residual ] = computeResidual( problem, coeff, parameters, tau, rho, teval, lambda_p )
% For a given (parameter-dependent) problem with approximate solution given
% by coeff on a mesh defined by tau and rho, the problem is evaluated at
% the points in teval.

orders = feval_problem(problem,'orders');

vals = zeros(length(orders),length(teval),max(orders)+1); %Pre-locate space for the values of the solution, the residual shall be computed of.

for oi = 0:max(orders)
    vals(:,:,oi+1) = coeffToValues(coeff,tau,orders,rho,teval,oi);
end

residual=zeros(length(orders),length(teval));
for k=1:length(teval)
   residual(:,k) =  feval_problem(problem,'problem',reshape(vals(:,k,:),length(orders),max(orders)+1),[],[],[],teval(k),parameters,[],[],lambda_p);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  getStandardCollocationPoints.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rho]=getStandardCollocationPoints(type,k)
if (strcmp(type,'gauss'))
%psi(:,:,i) for n-th order psi(:,:,n-i) equals i-th derivative
%for example for 3-rd order:
%psi(:,:,1) ... second derivative
%psi(:,:,2) ... first derivative
%psi(:,:,3) ... no derivative

switch k
   case 1
  rho=[0.5];
  psi(1,:,1)=[0 1 0];
  psi(1,:,2)=[0.5 0 0];
      case 2
       %C: rho: ci in Gaussian collocation, zeros of Legendre Polynomials shiftet to [0,1]
  %C: transformation from [-1,1] to [0,1]: y= 1/2x+1/2
  rho=[0.211324865405187117745425609748  0.788675134594812882254574390252];
  psi(1,:,1)=[0 -0.866025403784438646763723170755  1.36602540378443864676372317076 0];
  psi(2,:,1)=[0  0.866025403784438646763723170755   -0.366025403784438646763723170754 0];
  %0.th derivative
  %firs coll.point
  psi(1,:,2)=[-0.288675134594812882254574390252  0.683012701892219323381861585378 0 0];
  %second coll.point
  psi(2,:,2)=[0.288675134594812882254574390252  -0.183012701892219323381861585377 0 0];
     case 3
  rho=[1/2-1/10*15^(1/2) 1/2 1/2+1/10*15^(1/2)];
  psi(1,:,1)=[0 10/9 -5/3-1/6*15^(1/2) 5/6+1/6*15^(1/2) 0];
  psi(2,:,1)=[0 -20/9 10/3 -2/3 0];
  psi(3,:,1)=[0 10/9 -5/3+1/6*15^(1/2) 5/6-1/6*15^(1/2) 0];
  psi(1,:,2)=[5/18 -5/9-1/18*15^(1/2) 5/12+1/12*15^(1/2) 0 0];
  psi(2,:,2)=[-5/9 10/9 -1/3 0 0];
  psi(3,:,2)=[5/18 -5/9+1/18*15^(1/2) 5/12-1/12*15^(1/2) 0 0];
     case 4

  rho=[1/2-1/70*(525+70*30^(1/2))^(1/2) 1/2-1/70*(525-70*30^(1/2))^(1/2) 1/2+1/70*(525-70*30^(1/2))^(1/2) 1/2+1/70*(525+70*30^(1/2))^(1/2)];
  psi(1,:,1)=[0 -245/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/36*(105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             -7/24*(45+(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             1/120*(35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0];
  psi(2,:,1)=[0 245/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/36*(105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             -7/24*(-45+30^(1/2)-(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             1/120*(35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0];
  psi(3,:,1)=[0 -245/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/36*(-105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             7/24*(-45+30^(1/2)+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             1/120*(-35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0];
  psi(4,:,1)=[0 245/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/36*(-105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             7/24*(45-(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             1/120*(-35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0];
  psi(1,:,2)=[-49/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/144*(105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             -7/72*(45+(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             1/240*(35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0 0];
  psi(2,:,2)=[49/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/144*(105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             7/72*(45-30^(1/2)+(525-70*30^(1/2))^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) ...
             1/240*(35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0 0];
  psi(3,:,2)=[-49/24/(525-70*30^(1/2))^(1/2)*30^(1/2) -7/144*(-105+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             7/72*(-45+30^(1/2)+(525-70*30^(1/2))^(1/2))/(525-70*30^(1/2))^(1/2)*30^(1/2) ...
             1/240*(-35+(525-70*30^(1/2))^(1/2))*(-10+30^(1/2))*30^(1/2)/(525-70*30^(1/2))^(1/2) 0 0];
  psi(4,:,2)=[49/24/(525+70*30^(1/2))^(1/2)*30^(1/2) 7/144*(-105+(525+70*30^(1/2))^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             7/72*(45-(525+70*30^(1/2))^(1/2)+30^(1/2))/(525+70*30^(1/2))^(1/2)*30^(1/2) ...
             1/240*(-35+(525+70*30^(1/2))^(1/2))*(10+30^(1/2))*30^(1/2)/(525+70*30^(1/2))^(1/2) 0 0];

 case 5

  rho=[1/2-1/42*(245+14*70^(1/2))^(1/2) 1/2-1/42*(245-14*70^(1/2))^(1/2) 1/2 1/2+1/42*(245-14*70^(1/2))^(1/2) 1/2+1/42*(245+14*70^(1/2))^(1/2)];
  psi(1,:,1)=[0 567/25/(35+2*70^(1/2))*70^(1/2) -27/40*(84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/20*(343+9*(245+14*70^(1/2))^(1/2)+2*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             -3/280*(1911+77*(245+14*70^(1/2))^(1/2)+42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/280*(21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0];
  psi(2,:,1)=[0 567/25/(-35+2*70^(1/2))*70^(1/2) -27/40*(84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             -3/20*(-343+2*70^(1/2)-9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             3/280*(-1911+42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             -3/280*(21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0];
  psi(3,:,1)=[0 336/25 -168/5 1232/45 -112/15 8/15 0];
  psi(4,:,1)=[0 567/25/(-35+2*70^(1/2))*70^(1/2) 27/40*(-84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             -3/20*(-343+2*70^(1/2)+9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             -3/280*(1911-42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             3/280*(-21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0];
  psi(5,:,1)=[0 567/25/(35+2*70^(1/2))*70^(1/2) 27/40*(-84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             -3/20*(-343+9*(245+14*70^(1/2))^(1/2)-2*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/280*(-1911+77*(245+14*70^(1/2))^(1/2)-42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             -3/280*(-21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0];
  psi(1,:,2)=[189/50/(35+2*70^(1/2))*70^(1/2) -27/200*(84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/80*(343+9*(245+14*70^(1/2))^(1/2)+2*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             -1/280*(1911+77*(245+14*70^(1/2))^(1/2)+42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/560*(21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0 0];
  psi(2,:,2)=[189/50/(-35+2*70^(1/2))*70^(1/2) -27/200*(84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             3/80*(343-2*70^(1/2)+9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             1/280*(-1911+42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
             -3/560*(21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0 0];
  psi(3,:,2)=[56/25 -168/25 308/45 -112/45 4/15 0 0];
  psi(4,:,2)=[189/50/(-35+2*70^(1/2))*70^(1/2) 27/200*(-84+(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) -3/80*(-343+2*70^(1/2)+9*(245-14*70^(1/2))^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
            -1/280*(1911-42*70^(1/2)-77*(245-14*70^(1/2))^(1/2)+(245-14*70^(1/2))^(1/2)*70^(1/2))/(-35+2*70^(1/2))*70^(1/2) ...
            3/560*(-21+(245-14*70^(1/2))^(1/2))*(-14+70^(1/2))*70^(1/2)/(-35+2*70^(1/2)) 0 0];
  psi(5,:,2)=[189/50/(35+2*70^(1/2))*70^(1/2) 27/200*(-84+(245+14*70^(1/2))^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             3/80*(343-9*(245+14*70^(1/2))^(1/2)+2*70^(1/2))*70^(1/2)/(35+2*70^(1/2)) ...
             1/280*(-1911+77*(245+14*70^(1/2))^(1/2)-42*70^(1/2)+(245+14*70^(1/2))^(1/2)*70^(1/2))/(35+2*70^(1/2))*70^(1/2) ...
             -3/560*(-21+(245+14*70^(1/2))^(1/2))*(14+70^(1/2))*70^(1/2)/(35+2*70^(1/2)) 0 0];
        case 6

  rho=[0.33765242898423986094e-1 0.16939530676686774317 0.38069040695840154568 0.61930959304159845432 0.83060469323313225683 0.96623475710157601391];
  psi(1,:,1)=[0 -8.1412617290086756177 28.978672208695686824 -40.408362140839295291 27.785390555068868256 -9.6944498478780709320 1.5656732001510719330 0];
  psi(2,:,1)=[0 24.533738740695531559 -83.334379226361945114 107.81105264377978693 -64.863346973441684340 16.973778445028729194 -.94046284317634892902 0];
  psi(3,:,1)=[0 -36.168340495472103373 113.68329746902200059 -130.85406519770904404 65.101388388983683720 -12.145253252968680079 .61693005543048870860 0];
  psi(4,:,1)=[0 36.168340495472103373 -103.32674550381061965 104.96268528468059170 -44.848707621074554040 7.6576120141334378913 -.37922770211461375460 0];
  psi(5,:,1)=[0 -24.533738740695531559 63.868053217811244240 -59.145237622403034740 23.711846151968643422 -3.9123422341959200109 .19180001403866795482 0];
  psi(6,:,1)=[0 8.1412617290086756177 -19.868898165356366883 17.633927032490995438 -6.8865705015049570236 1.1206548758805039360 -.54712724329265912892e-1 0];
  psi(1,:,2)=[-1.1630373898583822311 4.8297787014492811373 -8.0816724281678590582 6.9463476387672170640 -3.2314832826260236440 .78283660007553596650 0 0];
  psi(2,:,2)=[3.5048198200993616513 -13.889063204393657519 21.562210528755957386 -16.215836743360421085 5.6579261483429097313 -.47023142158817446451 0 0];
  psi(3,:,2)=[-5.1669057850674433390 18.947216244837000098 -26.170813039541808808 16.275347097245920930 -4.0484177509895600263 .30846502771524435430 0 0];
  psi(4,:,2)=[5.1669057850674433390 -17.221124250635103275 20.992537056936118340 -11.212176905268638510 2.5525373380444792971 -.18961385105730687730 0 0];
  psi(5,:,2)=[-3.5048198200993616513 10.644675536301874040 -11.829047524480606948 5.9279615379921608555 -1.3041140780653066703 .95900007019333977410e-1 0 0];
  psi(6,:,2)=[1.1630373898583822311 -3.3114830275593944805 3.5267854064981990876 -1.7216426253762392559 .37355162529350131200 -.27356362164632956446e-1 0 0];
  case 7
  rho=[.254460438286207377369051579761e-1, .129234407200302780068067613360, .297077424311301416546696793962, .500000000000000000000000000000, .702922575688698583453303206038, .870765592799697219931932386640, .974553956171379262263094842024];
 case 8
  rho=[.198550717512318841582195657153e-1, .101666761293186630204223031762, .237233795041835507091130475405, .408282678752175097530261928820, .591717321247824902469738071180, .762766204958164492908869524595, .898333238706813369795776968238, .980144928248768115841780434285];
 case 9
  rho=[.159198802461869550822118985482e-1, .819844463366821028502851059651e-1, .193314283649704801345648980329, .337873288298095535480730992678, .500000000000000000000000000000, .662126711701904464519269007322, .806685716350295198654351019671, .918015553663317897149714894035, .984080119753813044917788101452];
 case 10
  rho=[.130467357414141399610179939578e-1, .674683166555077446339516557883e-1, .160295215850487796882836317443, .283302302935376404600367028417, .425562830509184394557586999435, .574437169490815605442413000565, .716697697064623595399632971583, .839704784149512203117163682557, .932531683344492255366048344212, .986953264258585860038982006042];
 case 11
  rho=[.108856709269715035980309994386e-1, .564687001159523504624211153480e-1, .134923997212975337953291873984, .240451935396594092037137165271, .365228422023827513834234007300, .500000000000000000000000000000, .634771577976172486165765992700, .759548064603405907962862834729, .865076002787024662046708126016, .943531299884047649537578884652, .989114329073028496401969000561];
 case 12
  rho=[.921968287664037465472545492536e-2, .479413718147625716607670669405e-1, .115048662902847656481553083394, .206341022856691276351648790530, .316084250500909903123654231678, .437383295744265542263779315268, .562616704255734457736220684732, .683915749499090096876345768322, .793658977143308723648351209470, .884951337097152343518446916606, .952058628185237428339232933060, .990780317123359625345274545075];
 case 13
  rho=[.790847264070592526358527559645e-2, .412008003885110173967260817496e-1, .992109546333450436028967552086e-1, .178825330279829889678007696502, .275753624481776573561043573936, .384770842022432602967235939451, .500000000000000000000000000000, .615229157977567397032764060549, .724246375518223426438956426064, .821174669720170110321992303498, .900789045366654956397103244791, .958799199611488982603273918250, .992091527359294074736414724404];
 case 14
  rho=[.685809565159383057920136664797e-2, .357825581682132413318044303111e-1, .863993424651175034051026286748e-1, .156353547594157264925990098490, .242375681820922954017354640724, .340443815536055119782164087916, .445972525646328168966877674890, .554027474353671831033122325110, .659556184463944880217835912084, .757624318179077045982645359276, .843646452405842735074009901510, .913600657534882496594897371325, .964217441831786758668195569689, .993141904348406169420798633352];
 case 15
  rho=[.600374098975728575521714070669e-2, .313633037996470478461205261449e-1, .758967082947863918996758396129e-1, .137791134319914976291906972693, .214513913695730576231386631373, .302924326461218315051396314509, .399402953001282738849685848303, .500000000000000000000000000000, .600597046998717261150314151697, .697075673538781684948603685491, .785486086304269423768613368627, .862208865680085023708093027307, .924103291705213608100324160387, .968636696200352952153879473855, .993996259010242714244782859293];
end
end
if(strcmp(type,'lobatto'))
if k<2 || k>15
   err('equations_err4');
end
switch k
   case 2
       rho=[0,1];
   case 3
       rho=[0., .500000000000000000000000000000, 1.];
   case 4
       rho=[0., .276393202250021030359082633127, .723606797749978969640917366873, 1.];
   case 5
       rho=[0., .172673164646011428100853771876, .500000000000000000000000000000, .827326835353988571899146228124, 1.];
   case 6
       rho=[0., .117472338035267653574498513019, .357384241759677451842924502980, .642615758240322548157075497020, .882527661964732346425501486981, 1.];
   case 7
       rho=[0., .84888051860716535063983893017e-1, .265575603264642893098114059046, .500000000000000000000000000000, .734424396735357106901885940954, .915111948139283464936016106983, 1.];
   case 8
       rho=[0., .641299257451966923312771193897e-1, .204149909283428848927744634301, .395350391048760565615671369827, .604649608951239434384328630173, .795850090716571151072255365699, .935870074254803307668722880610, 1.00000000000000000000000000000];
   case 9
       rho=[0., .501210022942699213438273777908e-1, .161406860244631123277057286454, .318441268086910920644623965646, .500000000000000000000000000000, .681558731913089079355376034354, .838593139755368876722942713546, .949878997705730078656172622209, 1.00000000000000000000000000000];
   case 10
       rho=[0., .402330459167705930855336695888e-1, .130613067447247462498446912570, .261037525094777752169412453634, .417360521166806487686890117021, .582639478833193512313109882979, .738962474905222247830587546366, .869386932552752537501553087430, .959766954083229406914466330411, 1.00000000000000000000000000000];
   case 11
       rho=[0., .329992847959704328338629319503e-1, .107758263168427790688791091946, .217382336501897496764518015261, .352120932206530304284044242220, .500000000000000000000000000000, .647879067793469695715955757780, .782617663498102503235481984739, .892241736831572209311208908054, .967000715204029567166137068050, 1.00000000000000000000000000000];
   case 12
       rho=[0., .275503638885588882962099308484e-1, .903603391779966608256792091415e-1, .183561923484069661168797572778, .300234529517325533867825104217, .431723533572536222567969072130, .568276466427463777432030927870, .699765470482674466132174895783, .816438076515930338831202427222, .909639660822003339174320790858, .972449636111441111703790069152, 1.00000000000000000000000000000];
   case 13
       rho=[0., .233450766789180440515472676223e-1, .768262176740638415670371964506e-1, .156905765459121286963620480217, .258545089454331899126531383182, .375356534946880003715663149813, .500000000000000000000000000000, .624643465053119996284336850187, .741454910545668100873468616818, .843094234540878713036379519783, .923173782325936158432962803549, .976654923321081955948452732378, 1.00000000000000000000000000000];
   case 14
       rho=[0., .200324773663695493224499189923e-1, .660994730848263744998898985459e-1, .135565700454336929707663799740, .224680298535676472341688647070, .328637993328643577478048298179, .441834065558148066170611645132, .558165934441851933829388354868, .671362006671356422521951701821, .775319701464323527658311352930, .864434299545663070292336200260, .933900526915173625500110101454, .979967522633630450677550081008, 1.00000000000000000000000000000];
   case 15
       rho=[0., .173770367480807136020743039652e-1, .574589778885118505872991842589e-1, .118240155024092399647940762012, .196873397265077144438235030682, .289680972643163759539051530631, .392323022318102880887160276864, .500000000000000000000000000000, .607676977681897119112839723136, .710319027356836240460948469369, .803126602734922855561764969318, .881759844975907600352059237988, .942541022111488149412700815741, .982622963251919286397925696035, 1.00000000000000000000000000000];
end
end
if(strcmp(type,'uniform'))
   rho = 1/(k+1):1/(k+1):k/(k+1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  higherchainrule.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ der ] = higherchainrule( n )
%CHAINRULE computes the chainrule up to the n-th derivative
der = cell(n+1,1);

der{1}.A=0; %is the 0-th derivate, i.e. f(tau(t))
der{1}.b=0;
der{1}.c=1;

for j=2:n+1 %compute the (j-1)-th derivative
    [M,N]=size(der{j-1}.A);
    [K,L]=find(der{j-1}.A);
    A = zeros(length(K),j);
    b = zeros(length(K),1);
    c = zeros(length(K),1);
    %first the tau-derivatives from the chain rule
    for ii = 1:length(K)
        A(ii,:)=[der{j-1}.A(K(ii),:),0]; %copy the row to the next A-matrix
        A(ii,L(ii)) = A(ii,L(ii)) - 1;
        A(ii,L(ii)+1) = A(ii,L(ii)+1) + 1;
        b(ii)=der{j-1}.b(K(ii)); %the tau-derivatives corresponds to the same f-derivative
        c(ii)=der{j-1}.c(K(ii))*der{j-1}.A(K(ii),L(ii)); %the new coefficient is the old coefficient times the power of the tau-derivative
    end
    %now the f-derivatives (note: the inner derivative gives one tau')
    c = [c;der{j-1}.c];
    b = [b;der{j-1}.b+1];
    A = [A;[der{j-1}.A,zeros(M,1)]+[zeros(M,1),ones(M,1),zeros(M,N-1)]];

    %Relaxation
    [U,IA,IU]=unique([A,b],'rows');
    der{j}.c = zeros(length(IA),1);
    for k = 1:length(IU)
       der{j}.c(IU(k))=der{j}.c(IU(k))+c(k);
    end
    der{j}.A=U(:,1:end-1);
    der{j}.b=U(:,end);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  trafo.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ret = trafo(varargin)
switch varargin{1}
    case 'xi'
        ret = trafo_xi(varargin{2},varargin{3},varargin{4});
    case 'xi_inv'
        ret = trafo_xi_inv(varargin{2},varargin{3},varargin{4});
    case 'trafomatrix'
        ret = trafo_trafomatrix(varargin{2},varargin{3},varargin{4});
    case 'trafomatrix_inv'
        ret = trafo_trafomatrix_inv(varargin{2},varargin{3},varargin{4});
    case 'splitting'
        ret = 1;
end
end

function ret = trafo_xi(tau,k,a)
%is the k-th derivative of the transformation from [0,1] on [a,infty)
alpha=-1;
if(a>0)
    ret=a./tau.^(k+1).*(-1)^(k);
else
    if(k>0)
        ret=prod(alpha-(0:k-1))*(1-tau).^(alpha-k);
    else
        ret=(1-tau).^alpha-1;
    end
end
end

function ret = trafo_xi_inv(t,k,a)
%is the k-th derivative of the transformation from [a,infty) on [0,1]
alpha=-1;
if(a>0)
    ret=a./t.^(k+1)*(-1)^(k);
else
    if(k>0)
        ret=prod(1/alpha-(0:k-1))*(1+t).^(1/alpha-k);
    else
        ret=1-(1+t).^(1/alpha);
    end
end
end

function [A] = trafo_trafomatrix(tau,n,a)
%trafo from z* defined on [0,1] on z defined on [a,infty)
global chainrulecoeffs;
A=zeros(n,n);
xitmp = zeros(n+1,1);
for i=1:n
    xitmp(i) = trafo_xi(tau,i-1,a);
end
for i=1:n
    b=chainrulecoeffs{i}.b;
    for j=1:length(b)
       A(i,chainrulecoeffs{i}.b(j)+1) = chainrulecoeffs{i}.c(j)*prod(xitmp(1:i)'.^(chainrulecoeffs{i}.A(j,:)));
    end
end
end

function [A] = trafo_trafomatrix_inv(tau,n,a)
%trafo from z defined on [a,infty) on z* defined on [0,1]
B=trafo_trafomatrix(tau,n,a);
if(norm(B)<Inf)
    A=inv(trafo_trafomatrix(tau,n,a));
else
    s=size(B);
    A=zeros(s(1),s(2));
end
A(1,:)=0;
A(:,1)=0;
A(1,1)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  computeEVPStart.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [initProfiles,D,V,L,R] = computeEVPStart(problem,settings,numberofguesses)

evalMesh = feval(settings,'mesh');
collMethod = 'gauss';
collPoints = 4;
try
    pathfoll=feval(problem,'pathfollowing');
    lambda_p=pathfoll.start;
catch
    lambda_p=0;
end
interval = feval_problem(problem,'interval',[],[],[],[],[],[],0,0,lambda_p);
x1=linspace(interval(1),interval(2),ceil(5*numberofguesses/collPoints)+1);
evalMesh = (interval(2)-interval(1))/(evalMesh(end)-evalMesh(1))*(evalMesh+evalMesh(1))+interval(1);
m = collPoints;
if (strcmp(collMethod,'user'))
    [rho] = collPoints;
else
    [rho] = getStandardCollocationPoints(collMethod,m);
end
n = feval_problem(problem,'n',[],[],[],[],[],[],0,0,lambda_p);
ordnung = feval_problem(problem,'orders',[],[],[],[],[],[],0,0,lambda_p);
parameter = feval_problem(problem,'parameters',[],[],[],[],[],[],0,0,lambda_p);
interval = feval_problem(problem,'interval',[],[],[],[],[],[],0,0,lambda_p);
c=feval_problem(problem,'c',[],[],[],[],[],[],0,0,lambda_p);
%a=interval(1);
b=interval(2);
N=length(x1)-1;
h(1:N)=x1(2:N+1)-x1(1:N);
tau(1:N,1:m)=x1(1:N).'*ones(1,m)+h(1:N).'*rho(1:m);
faktorielle=[1;1;2;6;24;120;720;5040;40320;362880;3628800;39916800;479001600;6227020800;];


psi=[];
psival=[];
p=zeros(parameter,1);
y=zeros(n,max(ordnung),N);
z=zeros(n,m,N);

phishift =  [0; sum(tril(ones(max(ordnung),max(ordnung)))*(ones(max(ordnung),1)*ordnung>=(1:max(ordnung))'*ones(1,length(ordnung))),2)];

if isempty(psi)% is psi defined by user?
   for ord=1:max(ordnung)
       for i=1:m
           psi(i,1+max(ordnung)-ord:m+max(ordnung),ord)=computeEVPStart_Psi(i,rho,ord);
       end
   end
end
if isempty(psival)%is psi defined by user?
   for ord=1:max(ordnung)
       for i=1:m
           %evaluation of psi
           psival(ord,i,1:m)=polyval(psi(i,:,ord),rho(1:m));
           psival(ord,i,m+1)=polyval(psi(i,:,ord),1);
           psival(ord,i,m+2)=polyval(psi(i,:,ord),0);
       end
   end
end

%For 0-th order special definitions for psi are necessary

%if min(ordnung)==0 || (length(was)==7 && length(strfind(was,'wertabl'))>0) || (length(was)==8 && length(strfind(was,'residual'))>0)
%   for i=1:m
%s       psi0(i,:,1)=Psi(i,rho,0);
%   end
%end

for j=1:(length(x1)-1)*(n*m+sum(ordnung))+parameter
             x0(j,1)=1;
end

nonzeroent=(sum(ordnung)+parameter)*(max(ordnung)*n*max(2,length(c))+m*n+parameter)+...
   N*m*n*max(ordnung)*n+...
   N*m*n*m*n+...
   N*max(ordnung)*n*max(ordnung)*n+...
   N*max(ordnung)*n*m+...
   N*max(ordnung)*n*max(ordnung)*n;
ROWl=zeros(nonzeroent,1);
COLl=zeros(nonzeroent,1);
OUTPUTl=zeros(nonzeroent,1);
ROWr=zeros(nonzeroent,1);
COLr=zeros(nonzeroent,1);
OUTPUTr=zeros(nonzeroent,1);

rowcount = 1;
countl = 1;
countr = 1;

% Randbedingungen links
% Anmerkung: muessen homogen sein!
dBV = feval_problem(problem,'dBV',[],[],[],[],[],[],0,0,lambda_p);
[size_dBV(1),size_dBV(2),size_dBV(3),size_dBV(4)] = size(dBV);

for i=1:size_dBV(2)
   for j=1:n
       for ll=1:size_dBV(4)
           OUTPUTl(countl) = dBV(1,i,j,ll);
           ROWl(countl) = rowcount;
           COLl(countl) = phishift(ll)+j;
           countl = countl + 1;

           for k=ll:ordnung(j)
               OUTPUTl(countl) = h(end)^(k-ll)/faktorielle(1+k-ll)*dBV(2,i,j,ll);
               ROWl(countl) = rowcount;
               COLl(countl) = (N-1)*(sum(ordnung)+n*m)+phishift(k)+j;
               countl = countl + 1;
           end

           for k=1:m
               OUTPUTl(countl) = polyval(computeEVPStart_Psi(k,rho,ordnung(j)-ll+1),1)*h(end)^(ordnung(j)-ll+1)*dBV(2,i,j,ll);
               ROWl(countl)=rowcount;
               COLl(countl)=(N-1)*(sum(ordnung)+n*m)+sum(ordnung)+(k-1)*n+j;
               countl = countl + 1;
           end
       end
   end
   %fprintf('Randbed in Zeile %d \n',rowcount);
   rowcount = rowcount+1;
end

%Innere Bedingungen: Kollokation und Stetigkeit

for i=1:N % Bedingung in welchem Intervall?
   for k=1:m % Bedingung an welchem Kollokationspunkt?
       tik = tau(i,k);
       dFl = feval_problem(problem,'jacobian',[],[],[],[],tik,p,0,0,lambda_p); % here: lambda = 0
       dFr = dFl-feval_problem(problem,'jacobian',[],[],[],[],tik,p,1,0,lambda_p); % here: lambda = 1
       [size_dF(1),size_dF(2),size_dF(3)] = size(dFl);
       %Linke Seite
       for ii=1:size_dF(1) % Nummer der Gleichung
          for jj=1:n % Welche Komponente
              %PHI-Koeffizienten
              for ll=1:ordnung(jj) % Bis zur welchen Ableitung
                  %Linke Seite
                 temp = 0;
                 for kk=0:ll-1 % Welche Ableitung
                     temp = temp + dFl(ii,jj,kk+1)*(tik-x1(i))^(ll-kk-1)/faktorielle(1+ll-kk-1);
                 end
                 OUTPUTl(countl) = temp;
                 ROWl(countl)=rowcount;
                 COLl(countl)=(i-1)*(sum(ordnung)+n*m)+phishift(ll)+jj;
                 countl = countl + 1;

                  %Rechte Seite
                 temp = 0;
                 for kk=0:ll-1
                     temp = temp + dFr(ii,jj,kk+1)*(tik-x1(i))^(ll-kk-1)/faktorielle(1+ll-kk-1);
                 end
                 OUTPUTr(countr) = temp;
                 ROWr(countr)=rowcount;
                 COLr(countr)=(i-1)*(sum(ordnung)+n*m)+phishift(ll)+jj;
                 countr = countr + 1;
              end
              %PSI-Koeffizienten
              for ll=1:m %PSI-Polynome
                  %Linke Seite
                  temp = 0;
                  for kk=0:ordnung(jj) %Ableitung
                      temp = temp + dFl(ii,jj,kk+1)*h(i)^(ordnung(jj)-kk)*polyval(computeEVPStart_Psi(ll,rho,ordnung(jj)-kk),(tik-x1(i))/h(i));
                  end
                  OUTPUTl(countl) = temp;
                  ROWl(countl)=rowcount;
                  COLl(countl)=(i-1)*(sum(ordnung)+n*m)+sum(ordnung)+(ll-1)*n+jj;
                  countl = countl + 1;
                  %Rechte Seite
                  temp = 0;
                  for kk=0:ordnung(jj) %Ableitung
                      temp = temp + dFr(ii,jj,kk+1)*h(i)^(ordnung(jj)-kk)*polyval(computeEVPStart_Psi(ll,rho,ordnung(jj)-kk),(tik-x1(i))/h(i));
                  end
                  OUTPUTr(countr) = temp;
                  ROWr(countr)=rowcount;
                  COLr(countr)=(i-1)*(sum(ordnung)+n*m)+sum(ordnung)+(ll-1)*n+jj;
                  countr = countr + 1;
              end
          end
          %fprintf('Kollokationsbed in Zeile %d \n',rowcount);
          rowcount = rowcount + 1;
       end
   end

   %Stetigkeit
   if(i<N)
       for j=1:n %Komponente
           for k=1:ordnung(j) %(k-1). Ableitung
               %Von links kommend
               for ll=k:ordnung(j)
                   OUTPUTl(countl)=h(i)^(ll-k)/faktorielle(1+ll-k);
                   ROWl(countl)=rowcount;
                   COLl(countl)=(i-1)*(sum(ordnung)+n*m)+phishift(ll)+j;
                   countl = countl + 1;
               end

               for ll=1:m
                   OUTPUTl(countl) = polyval(computeEVPStart_Psi(ll,rho,ordnung(j)-k+1),1)*h(i)^(ordnung(j)-k+1);
                   ROWl(countl)=rowcount;
                   COLl(countl)=(i-1)*(sum(ordnung)+n*m)+sum(ordnung)+(ll-1)*n+j;
                   countl = countl + 1;
               end

               %Von rechts kommend
               OUTPUTl(countl)=-1;
               ROWl(countl)=rowcount;
               COLl(countl)=(i)*(sum(ordnung)+n*m)+phishift(k)+j;
               countl = countl + 1;

               %fprintf('Stetigkeitsbed in Zeile %d \n',rowcount);
               rowcount = rowcount + 1;
           end
       end
   end
end

countl = countl-1;
countr = countr-1;

ROWl = ROWl(1:countl);
ROWr = ROWr(1:countr);
COLl = COLl(1:countl);
COLr = COLr(1:countr);
OUTPUTl = OUTPUTl(1:countl);
OUTPUTr = OUTPUTr(1:countr);

L = sparse(ROWl,COLl,OUTPUTl,N*(sum(ordnung)+n*m),N*(sum(ordnung)+n*m));
R = sparse(ROWr,COLr,OUTPUTr,N*(sum(ordnung)+n*m),N*(sum(ordnung)+n*m));

[V,D] = eig(full(L),full(R));

D=diag(D);
[~,I] = sort(abs(D));


V=V(:,I(1:N*m));
D=D(I(1:N*m));


%EF = zeros(length(evalMesh),N*m);

for i = 1:N*m
   temp = coeffToValues(V(:,i),x1,ordnung,rho,evalMesh);
   % EF(:,i) = sqrt(length(evalMesh))*EF(:,i)/norm(EF(:,i),2);
   temp =  sqrt(n*length(evalMesh))*temp/norm(temp,2);
   initProfile.initialMesh = evalMesh;
   initProfile.initialValues = temp;
   interval = feval(problem,'interval');
   if(interval(2) == Inf)
       if(interval(1) == 0 && trafo('splitting'))
           initProfile.initialMesh = unique([initProfile.initialMesh(1:end-1),1./initProfile.initialMesh(2:end)]);
           initProfile.initialValues = [initProfile.initialValues(1:end/2,1:(end-1)),initProfile.initialValues((end/2+1):end,end-1:-1:1)];
       else
           [initProfile.initialMesh,sortInd] = sort(trafo('xi',initProfile.initialMesh,0,interval(1)));
           initProfile.initialValues = initProfile.initialValues(:,sortInd);
       end
   end

   initProfile.lambda = D(i);
   initProfile.parameters = [];
   initProfiles{i} = initProfile;
end



end

function ret=computeEVPStart_P(abl,i,komp,x,m,x1,h,y,z,psi,ord) %Polynomial of derivation abl

faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];
sum=0;
for j=1:ord
   if j-abl-1>=0
       sum=sum+((x-x1((i)+1))^(j-abl-1))/(faktorielle(j-abl))*y(komp,j,(i)+1);
   end
end
if ord==0
   for j=1:m
       sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1+1),(x-x1((i)+1))/h((i)+1)));
   end
else
   for j=1:m
       sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1),(x-x1((i)+1))/h((i)+1)));
   end
end
ret=sum;

end


function ret=Pablspezial(abl,i,komp,x,m,x1,h,y,z,psi,ord) %Polynom der Ableitung abl

faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];
sum=0;
for j=1:ord
   if j-abl-1>=0
       sum=sum+((x-x1((i)+1))^(j-abl-1))/(faktorielle(j-abl))*y(komp,j,(i)+1);
   end
end
if ord-abl==0
   for j=1:m
       sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1+1),(x-x1((i)+1))/h((i)+1)));
   end
else
   for j=1:m
       sum=sum+h((i)+1)^(ord-abl)*(z(komp,j,(i)+1)*polyval(psi(j,:,ord-(abl+1)+1),(x-x1((i)+1))/h((i)+1)));
   end
end
ret=sum;
end

function ret=computeEVPStart_Poptimiert(abl,i,komp,x,m,x1,h,y,z,psival,ord,k)
%Polynomial opitmized version
faktorielle=[1 1 2 6 24 120 720 5040 40320 362880 3628800 39916800 479001600 6227020800];

ret=sum(((x-x1((i)+1)).^([abl+1:ord]-abl-1))./(faktorielle([abl+1:ord]-abl)).*y(komp,abl+1:ord,(i)+1))+...
   sum(h((i)+1)^(ord-abl).*(z(komp,1:m,(i)+1).*psival(ord-(abl+1)+1,1:m,k)));

end

function ret=evalP0(i,komp,m,x1,h,y,z,psi,ord)
if (nargin<9) ord=2; end;

     sum=0;
     for j=1:ord

        poly=[1 -x1((i)+1)];
        for k=2:(j-1)
                poly=conv(poly,[1 -x1((i)+1)]);
        end
        if (j-1)==0
            poly=[1];
        end
        zwisch(m+ord-j+1:m+ord)=poly/(factorial(j-1))*y(komp,j,(i)+1);
        sum=sum+zwisch;
    end
for j=1:m
     if ord==0
         help=polyvalpoly(psi(j,:,1),[1 -x1((i)+1)]/h((i)+1));
     else
         help=polyvalpoly(psi(j,:,ord-1+1),[1 -x1((i)+1)]/h((i)+1));
     end
     %Attetion: help has to be of length m+ord - correction for
     %mixed order
     help=help(length(help)-(m+ord-1):length(help));
     zwisch(1:m+ord)=h((i)+1)^ord*(z(komp,j,(i)+1)*help);
     sum=sum+zwisch;
end
ret=sum;
end

function out=polyvalpoly(psi,E)
%polynomial as value in polynomial receiving polynomial

PsiK(1,1)=1;t=1;
for k=2:length(psi)
   t=conv(t,E);
   for m=1:length(t)
     PsiK(k,m)=t(m);
   end
end
PsiE(((length(E)-1))*(length(psi)-1)+1)=PsiK(1,1)*psi(length(psi));
PsiEout=PsiE;
for k=2:length(psi)
   PsiE(((length(E)-1))*(length(psi)-1)+1-(length(E)-1)*(k-1):((length(E)-1))*(length(psi)-1)+1)=...
                      psi(length(psi)-k+1)*PsiK(k,1:(length(E)-1)*(k-1)+1);
   PsiEout=PsiEout+PsiE;
end

out=PsiEout;
end

function Psireturn=computeEVPStart_Psi(n,rho,nr)
 i=length(rho);
 prod=[1];
 for s=1:i
     if (s~=n)
         prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
     end
 end
 for s=1:(nr)
     prod=polyint(prod);
 end
 Psireturn=prod;
end

function [A] = computeEVPStart_trafomatrix(t,n)
n=n-1;
A=diag([1,(-1).^(1:n)./t.^(2*(1:n))]);
if(n>1)
    A(1+(2:n),2) = ((-1).^(2:n).*factorial(2:n)./t.^(1+(2:n)))';
    for i=4:n+1
        A(i,3:i-1) = -(A(i-1,2:i-2)/t^2 + A(i-1,3:i-1)/t);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  pathfollowing.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [speicher,speicher_exact,tur_pts]=pathfollowing(problem,settings,pathfoll,x1,rho)
% Pathfollowing routine

%% Set options
eval_exact  = 1; % evaluate at the given parameter values
save_ws     = 1; % enable taking step(s) back if desired
dispexact   = 1; % plot at the given parameter values
dispres     = 1; % display the solution and parameter evolution plots
disppredcor = 1; % display with predictor and corrector step
log         = 1; % log some values
displog     = 1; % plot the logged values
dispmdep    = 1; % plot mesh density evolution
display     = 1; % display messages
% For secondary options, press Ctrl+F and search for 'Secondary'

%% Initializing
m           = length(rho);
xfin        = x1;
counter_mc  = 0; % Counter for mesh correction - DISABLED FOR NOW
ordnung     = feval_problem(problem,'orders');
n           = length(ordnung);
interval    = feval(problem,'interval');
infsplit    = 0;
se_nn       = 40;
% In case the right hand side of the interval is infinite, then x1 is splitted in [0,1] and [1,infty) when using splitting interval transformation
selim       = 0;
if interval(2) == Inf
    if interval(1) == 0 && trafo('splitting')
        infsplit = 1;
        n = n/2;
    end
    se_nn = 40;
    % Secondary option
    % The solution evolution variable 'selim' sets up to which value on the x-axis the solutions evolution should be drawn when the solution is computed for a problem posed on a semi-infinite interval
    selim = 10;
end
AbsTol     = feval(settings,'absTolSolver');
RelTol     = feval(settings,'relTolSolver');
minmeshpts = max(feval(settings,'minInitialMesh'),length(feval(settings,'mesh'))); % DISABLED for now

% Set the variables from user input or solver settings
tmp_fn = {'thetaMax','maxCorrSteps','maxSteplengthGrowth','angleMin','PredLengthFactor','CorrLengthGrowth','meshFactorMax'};
tmp_val = cell(1,numel(tmp_fn));
for ii=1:numel(tmp_fn)
    if isfield(pathfoll,tmp_fn{ii})
        tmp_val{ii}=pathfoll.(tmp_fn{ii});
    else
        tmp_val{ii}=feval(settings,tmp_fn{ii});
    end
end
[theta_max,maxCorrSteps,maxSteplengthGrowth,cos_min,pred_lf,corr_lf,meshFactorMax] = deal(tmp_val{:});

try
    data = pathfoll.pathdata;
catch
    msg = 'There is no ret.pathdata defined in the problem definition. Please refer to the manual of bvpsuite2.0 for some examples.';
    error(msg)
end
halve       = 0;
halve_nb    = 0; % times the step-length has been halved
corr_old    = Inf;
min_sl      = 1e-8; % Secondary option: minimal steplength
min_sl_flag = 0; % flag for minimal steplength

psi=zeros(m,max(ordnung)+m,max(ordnung));
for ord=1:max(ordnung)
    for ii=1:m
        psi(ii,1+max(ordnung)-ord:m+max(ordnung),ord)=pathfollowing_Psi(ii,rho,ord);
    end
end
psival=zeros(max(ordnung),m,m+2);
for ord=1:max(ordnung)
    for ii=1:m
        %evaluation of psi
        psival(ord,ii,1:m)=polyval(psi(ii,:,ord),rho(1:m));
        psival(ord,ii,m+1)=polyval(psi(ii,:,ord),1);
        psival(ord,ii,m+2)=polyval(psi(ii,:,ord),0);
    end
end

itnum = 0; % Number of completed steps
skip = 0; % Marker in case tagents angle is too steep
% Current number of steps until prompt
cter_flag = 0;
if isfield(pathfoll,'counter')
    counter = pathfoll.counter;
    if counter==Inf || (isfield(pathfoll,'only_counter') && pathfoll.only_counter)
        cter_flag = 1;
    end
else
    counter = 1;
end

prompt='\nType <strong>enter</strong> to carry out the chosen number of steps (ret.counter/1 when not modified),\n<strong>p</strong> and <strong>enter</strong> to change the maximal length of the predictor step,\n<strong>n</strong> and <strong>enter</strong> to change this chosen number,\n<strong>f</strong> and <strong>enter</strong> to plot&display solutions of computed steps, and \n<strong>s</strong> and <strong>enter</strong> to stop&save the computations\n';
subprompt1='Choose a number of steps: ';
subprompt2='Choose a maximal predictor step length: ';
subprompt3='Choose number of steps to go back: ';
subprompt4='Input a row vector [ . . . ] of steps: ';
subprompt5='Enter one of the options above to proceed: ';
f_size=8; % font size for the plots

% keep the previous interpreter settings and set all to latex
dtiii=get(0, 'DefaultTextInterpreter');
dliii=get(0, 'DefaultLegendInterpreter');
datliii=get(groot, 'DefaultAxesTickLabelInterpreter');
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')

%% Prepare what is needed for the options
% load the exact values of the parameter that shall be computed
if isfield(pathfoll,'require_exact')
    exact_val   = sort(unique(pathfoll.require_exact));
    if strcmp(pathfoll.startat,'start')
        speicher_exact = {};
    end
else
    exact_val = [];
    speicher_exact = [];
end
if isfield(pathfoll,'max_pred_length')
    max_pred = pathfoll.max_pred_length;
    adapt=0;
else
    max_pred = [];
    adapt=0;
end
if isfield(pathfoll,'pit_stop')
    pit_stop = pathfoll.pit_stop;
else
    pit_stop = [];
end
% If only_exact is activated, then only speicher_exact is computed
if isfield(pathfoll,'only_exact') && ~isempty(exact_val)
    only_exact=pathfoll.only_exact;
    if only_exact
        eval_exact=1;
    end
else
    only_exact=0;
end
% to save the workspace in order to go back 'counter'-number of steps in the iteration
if save_ws
    prompt_ws='\nDo you want to go back some steps?\nIf yes, press <strong>b</strong> and <strong>enter</strong>.\nIf no just <strong>enter</strong>.\n';
    jump_cell = {};
end
% Which plots should be drawn?
if ~eval_exact
    dispexact = 0;
end
if isempty(data) % data should actually never be empty...
    dispres = 0;
end
if ~dispres
    disppredcor = 0;
end
if ~log && displog
    displog = 0;
end
if~(feval(settings,'meshAdaptation'))
    dispmdep = 0;
end
if dispmdep || dispres || displog
    pathfollplot=figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6],'PaperUnits','normalized');
    ax = gca;
    ax.FontSize=f_size;
    clf ;
    pathfollplot.Color = 'White' ;
end
if dispmdep
    % Mesh density evolution plot
    mdep_nn = 50;
    X_mdep  = linspace(0,1,mdep_nn);
    Y_mdep  = [];
    Z_mdep  = [];

    % Mesh evolution plot
    pts_nb = 20;
    if length(x1)<=pts_nb
        filler = interval(2)*ones(1,length(x1)-pts_nb);
        Pts_Mat=[x1 filler];
    else
        pts_step = ceil(length(x1)/pts_nb);
        filler = interval(2)*ones(1,pts_nb-length(x1(1:pts_step:end-1)));
        Pts_Mat = [x1(1:pts_step:end-1) filler];
    end

    % Call the figures
    if dispres
        path_pos = [0 0.5 1 0.5];
        mdep_pos = [0 0 1 0.5];
    else
        mdep_pos = [0 0 1 1];
    end
elseif dispres
    path_pos = [0 0 1 1];
end
if dispres
    figure(pathfollplot);
end
if log
    logstruct.theta_0         = [] ;
    logstruct.theta_max       = [] ;
    logstruct.steplength      = [] ;
    logstruct.steplength_pred = [pathfoll.steplength] ;
    logstruct.sl_adapt        = [];

    logstruct.norm_F          = [];
    logstruct.cpt             = [] ;
    logstruct.orthogonal      = [] ;
    logstruct.max_error       = [];
    logstruct.mesh_length     = [] ;

    logstruct.c_s             = Inf  ;
    logstruct.cos_ab          = Inf;
    logstruct.corr_dist       = Inf ;
    logstruct.norm_delta_0    = [] ;
    logstruct.norm_delta_1    = [] ;

    logstruct.mesha_halved    = [];
    logstruct.angle_halved    = [];
    logstruct.trm_halved      = [];
    logstruct.pdist_halved    = [];
    logstruct.cdist_halved    = [];

    % logstruct fields that are marks in the plot whenever steplength got halved for one of the reasons
    % 1. mesh adaptation was too long
    % 2. angle between consecutive steps was too big
    % 3. trust region method was activated
    % 4. corrector step was too long compared to predictor step
    % 5. corrector step was too long compared to previous corrector step
    fn_halve = fieldnames(logstruct);
    fn_halve = fn_halve(end-4:end);
    val_halve = [ 0 1 0 1 10 ];

    if displog
        if dispres && dispmdep
            path_pos = [0.5 0.5 0.5 0.5];
            mdep_pos = [0.5 0 0.5 0.5];
            log_pos  = [0 0 0.5 1];
        elseif dispres
            path_pos = [0.5 0 0.5 1];
            log_pos  = [0 0 0.5 1];
        elseif dispmdep
            mdep_pos = [0.5 0 0.5 1];
            log_pos  = [0 0 0.5 1];
        else
            log_pos  = [0 0 1 1];
        end
    end
else
    logstruct = [];
end
exact_lp = Inf;

%% Preprocess the two cases: Start from initial point or loaded data
% Start at the given initial point on the solution path
if strcmp(pathfoll.startat,'start')
    % Set the starting value for lambda_p
    lambda_p=pathfoll.start;

    % Before computing a predictor corrector step, the initial solution needs to be computed first
    itnum=itnum-1;

    % If the predictor and corrector steps should be drawn
    if disppredcor
        val = {};
    else
        val = [];
    end

    % Inital steplength and initalize struct for predictor, used to carry over information from the pathfollowing routine into the bvpsuite2.0 routine
    steplength = pathfoll.steplength;
    predictor  = struct;

    % Save the important values in predictor, that matter for the computation
    predictor.steplength = steplength;
    predictor.infsplit   = infsplit;
    predictor.lpfix      = 0;
    predictor.x1         = x1;
    predictor.lambda_p_0 = lambda_p;
    predictor.skip       = 0;

    % Also some solver settings (needed in PredCorrStrat and meshadaptation.m) and display setting
    predictor.maxCorrSteps        = maxCorrSteps;
    predictor.maxSteplengthGrowth = maxSteplengthGrowth;
    predictor.meshFactorMax       = meshFactorMax;
    predictor.display             = display;

    % Starting point data for the Newton iteration in case of nonlinear problem
    a_p=pathfoll.initProfile;

    % Search for the next given parameter value, for which the value should be computed, if required
    if ~isempty(exact_val)
        if steplength >0
            exact_lp = find(exact_val>=lambda_p,1);
        else
            exact_lp = find(exact_val<=lambda_p,1,'last');
        end
    end

    % Define some values needing to exist for the function call of Postprocess later, they do not have any specific values as of now
    tangent_new=0; theta_0=0; delta_0=0;
else
    try
        filename=strcat(pathfoll.dir,pathfoll.startat);
        load(filename);
    catch MException
        error(MException.message);
    end
    if ~only_exact
        speicher_exact = p_exact;
    else
        speicher_exact = {};
    end
    speicher = p_save(:,1:end-1); % savestruct column is not needed

    % Define x1 and a_p
    sol       = speicher{1,end-1}; % sol is the penultimate solution
    x1        = sol.x1;
    a_0       = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
    % Define sol, predictor, xfin and a_c
    sol       = speicher{1,end}; % sol is the last computed solution
    predictor = sol.predictor;
    a_p       = a_0 + predictor.steplength.*predictor.tangent';
    xfin      = sol.x1;
    a_c       = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;

    savestruct = p_save{1,end};
    if log
        logstruct   = p_save{2,end};
    end

    jump_cell   = p_save{3,end};

    % Reassign all the variables from savestruct
    delta_0         = savestruct.delta_0;
    theta_0         = savestruct.theta_0;
    itnum           = savestruct.itnum;
    val             = savestruct.val;
    exact_lp        = savestruct.exact_lp;
    pred_tmp        = savestruct.pred_tmp;
    corr_old        = savestruct.corr_old;
    if isfield(savestruct,'halve_nb')
        halve_nb    = savestruct.halve_nb;
    else
        halve_nb    = 5;
    end
    if ~exist('pit_stop','var')
        pit_stop    = savestruct.pit_stop;
    end
    if dispmdep
        try
            Pts_Mat = savestruct.Pts_Mat;
            Y_mdep  = savestruct.y_mdep;
            Z_mdep  = savestruct.z_mdep;
        catch
            Y_mdep  = [];
        end
    end

    % If data has been changed, then recompute for the new data
    if max(data(xfin,a_c,ordnung,rho)~=speicher{3,end})
        [speicher,speicher_exact,val,corr_old] = correct_data(speicher,speicher_exact,val,data,ordnung,rho);
    end

    steplength = predictor.steplength;
    tangent_new = predictor.tangent';
    % Compute eval_exact if necessary
    if ~isempty(exact_val) && ~only_exact
        eval_exact = 1;
        % if sol.lambda_p<exact_val(1) and step-length<0 or
        tmp_bool1 = sol.lambda_p<exact_val(1) && steplength<0;
        % if sol.lambda_p>exact_val(end) and step-length>0
        tmp_bool2 = exact_val(end)<sol.lambda_p && steplength>0;
        % then for now, eval_exact is set to 0
        if tmp_bool1 || tmp_bool2
            eval_exact=0;
        end
    end

    % Update solver settings in case it was changed by user when reloading the run
    predictor.maxCorrSteps = maxCorrSteps;
    predictor.maxSteplengthGrowth = maxSteplengthGrowth;
    predictor.meshFactorMax = meshFactorMax;

    if ~only_exact
        % Prepare for the next step
        [a_0,a_p,tangent_new,steplength,predictor,delta_0,nda,nsd,theta_0,logstruct,eval_exact,exact_lp] = PredCorrStrat(problem,settings,predictor,a_c,a_p,x1,xfin,sol,tangent_new,ordnung,rho,psival,psi,steplength,exact_val,exact_lp,delta_0,theta_0,theta_max,log,logstruct,eval_exact,itnum,halve_nb,min_sl);
    else
        eval_exact = 1;
        exact_lp = Inf;
    end

    % Accept the new mesh for the next step
    x1 = xfin;

    lambda_p = speicher{2,1};

    itnum=itnum+1;
end

%% Pathfollowing routine
if itnum==-1 % Initial solution is computed first
    jj = 0;
else % Initial solution was already computed
    jj = 1;
end
% Secondary option:
% After jj_max steps, the routine will stop & save. This number may be changed. This is a safety measure for the case where 'counter' is set to Inf.
jj_max = 1e4;
if counter==Inf
    counter=jj_max;
end

% Main loop: is only exited as soon as jj==counter and then 's' is chosen in the user prompt, either by the user or automatically.
while 1
    % Prechecks of the step-length
    if itnum~=-1 && ~only_exact
        % Check whether the step-length is shorter than the minimal steplength
        if abs(steplength)<min_sl
            min_sl_flag=1; % Set flag to stop later
            steplength=sign(steplength)*min_sl;
            predictor.steplength=steplength;
            a_p = a_0 + steplength.*tangent_new;
            predictor.lambda_p_p=a_p(end);
            if display
                fprintf('\n  The steplength was shorter than the minimal steplength %1.1e!\n  It has been augmented to the minimal steplength.\n',min_sl)
            end
        end
        % Check whether the user-defined maximal step-length is overreached, if it is, adapt the predictor to a closer step
        data_tmp=data(x1,a_p,ordnung,rho);
        pred_tmp=sqrt((a_p(end)-speicher{2,end})^2+(max(abs(speicher{3,end}-data_tmp)))^2);
        if ~isempty(max_pred) && pred_tmp>max_pred % pred_tmp~=max_pred
            div_tmp=max_pred/pred_tmp;
            steplength=steplength*div_tmp;
            predictor.steplength=steplength;
            a_p = a_0 + steplength.*tangent_new;
            predictor.lambda_p_p=a_p(end);

            pred_tmp = max_pred;

            if display
                fprintf('\n  Predictor step was too long. Steplength was reduced to %f.\n',steplength)
            end

            if log
                % save whether the step-length was dimished or not
                adapt = 1;
            end
        end
    end

    if only_exact
        fprintf('\nComputing the approximate solutions for lambda_p=[ require_exact ]:\n')
    elseif itnum~=-1 && display
        fprintf('\n  Starting at predictor value %f:\n\n',a_p(end))
    elseif display
        fprintf('\n  Starting at predictor value %f:\n\n',lambda_p)
    end

    % Corrector Step
    if halve==0
        cpt = cputime;
    end
    if log && itnum~=-1 && ~only_exact
        if skip==0
            for ii=1:length(fn_halve)
                logstruct.(fn_halve{ii}) = [ logstruct.(fn_halve{ii}) val_halve(ii)/(halve==ii) ];
            end
        else
            logstruct.(fn_halve{halve})(end) = val_halve(halve);
        end
    end
    if only_exact
        itnum = itnum-1;
    elseif (feval(settings,'meshAdaptation'))
        [sol,~,halve] = meshadaptation(problem,settings,x1,rho,a_p,predictor);
        % If points have been added in the mesh, then reset the counter, otherwise wait for ? number of iterations in which the mesh does not change, to halve the mesh length if needed
        % ..................... DISABLED FOR NOW .....................
        if length(x1)>length(xfin)
            counter_mc = 0;
        else
            counter_mc = counter_mc + 1;
        end
    elseif(feval(settings,'errorEstimate'))
        if ~isfield(predictor,'tangent') && feval(problem,'linear')
            [~,~,sol]=solveLinearProblem(problem,x1,rho);
        else
            [~,~,sol,halve] = solveNonLinearProblem(problem,settings,x1,rho,a_p,1,predictor);
        end
        if itnum==-1
            sol.lambda_p = lambda_p;
        end
        if halve==0
            initCoeff = [sol.coeff', sol.parameters, sol.lambda_p ]';
            sol.initCoeff = initCoeff;
            [sol,~,halve] = errorestimate(sol,problem,settings,predictor);
        end
    else
        if ~isfield(predictor,'tangent_new') && feval(problem,'linear')
            [~,~,sol]=solveLinearProblem(problem,x1,rho);
        else
            [~,~,sol,halve] = solveNonLinearProblem(problem,settings,x1,rho,a_p,1,predictor);
        end
        if log
            sol.errest = Inf;
        end
    end
    if itnum==-1 && ~only_exact && isstruct(sol)
        sol.lambda_p = lambda_p;
    end

    % If the number of mesh points does not need halving, then compute the angle between the tangents and the lengths of the predictor and corrector steps
    if halve==0 && ~only_exact
        % Save the predictor for possible future use
        sol.predictor=predictor;

        % Save the result of the corrector step on the solution path and the new mesh xfin
        a_c = [ sol.coeff ; sol.parameters ; sol.lambda_p ];
        predictor.lambda_p_0 = a_c(end);
        xfin = sol.x1;

        while itnum~=-1
            % Check whether the corrector step is short enough compared to the predictor step
            speicher_tmp=speicher;
            speicher_tmp{1,end+1}=sol;
            speicher_tmp{2,end}=a_c(end);
            speicher_tmp{3,end}=data(xfin,a_c,ordnung,rho);
            corr_tmp=sqrt((a_c(end)-a_p(end))^2+(max(abs(speicher_tmp{3,end}-data_tmp)))^2);
            factor1_tmp=corr_tmp/pred_tmp;
            factor2_tmp=corr_tmp/corr_old;
            if display
                fprintf('\n  Corrector step is %f-times as long as predictor step.\n',factor1_tmp)
            end
            if itnum~=0
                if display
                    fprintf('  Corrector step is %f-times as long as the previous corrector step.\n',factor2_tmp)
                end
            end
            if factor1_tmp>1/pred_lf
                halve=4;
                if display
                    fprintf('\n  Corrector step is at most allowed to be %f-times as long as predictor step!\n',1/pred_lf)
                end
                if ~min_sl_flag
                    break
                end
            end
            if factor2_tmp>corr_lf
                halve=5;
                if display
                    fprintf('\n  Corrector step is at most allowed to be %f-times as long as previous corrector step!\n',corr_lf)
                end
                if ~min_sl_flag
                    break
                end
            end

            if itnum~=0
                % Compute the cosine between the last computed step and the next tangent
                if (feval(settings,'meshAdaptation')) && (length(x1)~=length(xfin) || max(abs(x1-xfin))>1e-12)
                    initP.initialMesh=xfin;
                    initP.parameters = sol.parameters;

                    initP.initialValues = coeffToValues(tangent_new, x1,ordnung,rho,xfin);
                    initP = initial_coefficients(problem,xfin,initP,rho,0);
                    tang_a = [ initP.initialCoeff ; tangent_new(end) ] ;
                else
                    tang_a = tangent_new;
                end
                jac_F = functionFDF( 'DF', problem ,a_c,xfin,psival,psi,rho,[]);
                tang_b = tangente_berechnen(jac_F);
                a_p_b = a_c + sign(steplength)*sign(tang_a.'*tang_b)*tang_b;
                data_tang=data(xfin,a_p_b,ordnung,rho);
                nn_data = length(data_tang);
                a=[(speicher_tmp{2,end}-speicher_tmp{2,end-1})*ones(nn_data,1) speicher_tmp{3,end}-speicher_tmp{3,end-1}];
                b=[(a_p_b(end)-speicher_tmp{2,end})*ones(nn_data,1) data_tang-speicher_tmp{3,end}];
                tmp_ab=zeros(1,nn_data);
                for ii=1:nn_data
                    n_a=a(ii,:)*a(ii,:)';
                    n_b=b(ii,:)*b(ii,:)';
                    tmp_ab(ii)=a(ii,:)*b(ii,:)'/sqrt(n_a*n_b);
                end
                cos_ab=min(tmp_ab);
                if display
                    fprintf('  Cosine of the angle between current step and next tangent: %1.1e\n',cos_ab)
                end
                % If the cosine smaller than the user set cosine, then the steplength is halved
                if cos_ab<cos_min
                    halve=2;
                    if display
                        fprintf('\n  The cosine is smaller than the minimal allowed cosine %1.1e\n',cos_min)
                    end
                end
            end
            break
        end
    end

    if itnum~=-1 && halve~=0 && ~min_sl_flag
        steplength = steplength/2;
        predictor.steplength=steplength;
        a_p = a_0 + steplength.*tangent_new;
        predictor.lambda_p_p=a_p(end);

        if display
            fprintf('  Steplength was halved from %f to %f!\n',2*steplength,steplength)
        end

        itnum=itnum-1; % itnum is augmented by 1 at the end of each loop
        skip = skip+1;
        predictor.skip = skip;
        counter = counter + 1;
    else % Proceed to finishing the step
        if ~only_exact
            halve_nb = skip;
            halve = 0;
            counter = counter-skip;
            jj = jj-skip;
            skip=0;
            predictor.skip = 0;
            if min_sl_flag && ~isstruct(sol) % pathfollowing is stuck
                msg1='The steplength is reduced too much! The pathfollowing routine has trouble following the path with the current solver settings! The current mininmal steplength (min_sl) is ';
                msg2='. To change this value look for min_sl at the beginning of the file pathfollowing.m in the bvpsuite2.0 package. Otherwise adapting the solver settings, especially using mesh adaptation (MA) and choosing the MA tolerances only slightly looser than the solver tolerances (i.e. solver tol: 10^-6, MA tol: 10^-4), may help.';
                msg=[msg1 num2str(min_sl) msg2];
                error(msg)
            end

            % Save the pathdata of interest in the variable speicher
            if itnum==-1
                speicher{1,1}=sol;
                speicher{2,end}=a_c(end);
                speicher{3,end}=data(xfin,a_c,ordnung,rho);
            else
                speicher=speicher_tmp;
                corr_old=corr_tmp;
            end
        end

        % Check if one of the required values of lambda_p has been passed if yes, compute the solution in this point and save the solution in speicher_exact
        if itnum~=-1 && only_exact==1 % only speicher_exact is computed in this run
            % Restart at the beginning & find the next exact value to compute
            lp_a0 = cell2mat(speicher(2,:));
            % Find the value exact_val(exact_lp) that is closest to the starting value of the parameter lambda_p
            for ii=2:length(lp_a0)
                tmp_sign = sign(lp_a0(ii)-lp_a0(ii-1));
                if tmp_sign==1
                    tmp_line = 1:length(exact_val);
                else
                    tmp_line = length(exact_val):-1:1;
                end
                for kk = tmp_line
                    tmp_bool1=tmp_sign*lp_a0(ii-1)<=tmp_sign*exact_val(kk) && tmp_sign*exact_val(kk)<tmp_sign*lp_a0(ii);
                    % If the required exact value is in-between two values of the path, then save ii and kk and break
                    if tmp_bool1
                        exact_lp=kk; % is the entry number in exact_val
                        cter_a0=ii-1; % is the step number in the path leading to a_0
                        sol = speicher{1,cter_a0+1}; % is the solution in the path leading to a_c exact_val(exact_lp) is in-between a_0(end) and a_c(end)
                        break
                    end
                end
                if exact_lp~=Inf
                    break
                end
            end
            % If the for-loop did not break
            if exact_lp==Inf
                msg='None of the required values in ret.require_exact is within the range of the computed path!';
                error(msg)
            end
        elseif itnum~=-1 && ~isempty(exact_val) % Compute lp_a0 and cter_a0 to check whether a required value has been passed or not
            lp_a0 = a_0(end);
            cter_a0 = 1;
        end
        % Check whether or not one of the values that should be computed at exactly, was already passed, if yes, compute at the value(s) that have been passed and are required
        if itnum~=-1 && ~isempty(exact_val) && eval_exact && abs(lp_a0(cter_a0)-sol.lambda_p)> abs(lp_a0(cter_a0)-exact_val(exact_lp))
            % More than one value can be in the interval that was computed, thus a while statement here
            while eval_exact && abs(lp_a0(cter_a0)-sol.lambda_p)> abs(lp_a0(cter_a0)-exact_val(exact_lp))
                if display
                    fprintf('\n<strong>Exact solution</strong> at %f computation:\n\n',exact_val(exact_lp))
                end

                % If by chance the value of lambda in the last computation was exactly the required point, then no computaiton is needed, otherwise it is
                if abs(lp_a0(cter_a0)-exact_val(exact_lp)) > 0
                    % Define sol_tmp and a_0_tmp whether only_exact is 1 or 0
                    if only_exact
                        sol_tmp = speicher{1,cter_a0};
                        a_0_tmp = [ sol_tmp.coeff ; sol_tmp.parameters ; sol_tmp.lambda_p ] ;
                    else
                        sol_tmp = sol;
                        a_0_tmp = a_0;
                    end
                    predictor = sol_tmp.predictor;
                    tangent = predictor.tangent';
                    sl_tmp = predictor.steplength;
                    predictor.lpfix = 1;
                    predictor.lambda_p_0 = exact_val(exact_lp);
                    a_p_exact = a_0_tmp + (exact_val(exact_lp)-lp_a0(cter_a0))/tangent(end)*tangent;
                    if(feval(settings,'meshAdaptation'))
                        [sol_exact] = meshadaptation(problem,settings,x1,rho,a_p_exact,predictor);
                    elseif(feval(settings,'errorEstimate'))
                        [~,~,sol_exact] = solveNonLinearProblem(problem,settings,x1,rho,a_p_exact,1,predictor);
                        initCoeff = [sol_exact.coeff', sol_exact.parameters, sol_exact.lambda_p ]';
                        sol_exact.initCoeff = initCoeff;
                        sol_exact = errorestimate(sol_exact,problem,settings,predictor);
                    else
                        [~,~,sol_exact] = solveNonLinearProblem(problem,settings,x1,rho,a_p_exact,1,predictor);
                    end
                    predictor.lpfix = 0;
                else % Solution was already computed
                    if ~only_exact
                        sol_exact = speicher{1,end-1};
                    else
                        sol_exact = speicher{1,cter_a0};
                    end
                    % Load the steplength of the last computed step
                    predictor = sol.predictor;
                    sl_tmp = predictor.steplength;
                end
                speicher_exact{1,end+1} = sol_exact;
                tmp_vct = [ sol_exact.coeff ; sol_exact.parameters ; sol_exact.lambda_p ] ;
                speicher_exact{2,end} = tmp_vct(end);
                speicher_exact{3,end} = data(sol_exact.x1,tmp_vct,ordnung,rho);

                if display
                    fprintf('\n<strong>Exact solution</strong> computed.\n')
                end
                predictor.lambda_p_0 = a_c(end); % In case it was changed

                if dispexact
                    plot_step(sol_exact,f_size,infsplit,n,selim)
                end

                % take the next value of the required values for next time
                if ~only_exact
                    exact_lp = exact_lp + sign(sl_tmp);
                    if exact_lp==0 || exact_lp>length(exact_val)
                        exact_lp = exact_lp - sign(sl_tmp);
                        eval_exact=0;
                    end
                else
                    % exact_lp = exact_lp+1?
                    tmp_bool1= sl_tmp>0 && exact_lp<length(exact_val) && (abs(lp_a0(cter_a0)-sol.lambda_p) > abs(lp_a0(cter_a0)-exact_val(exact_lp+1)));
                    if tmp_bool1
                        exact_lp = exact_lp+1;
                    end
                    % or exact_lp = exact_lp-1?
                    tmp_bool2= sl_tmp<0 && exact_lp>1 && (abs(lp_a0(cter_a0)-sol.lambda_p) > abs(lp_a0(cter_a0)-exact_val(exact_lp-1)));
                    if tmp_bool2
                        exact_lp = exact_lp-1;
                    end
                    % if neither one, then search for the next exact_val(exact_lp) along the path
                    if ~(tmp_bool1 || tmp_bool2)
                        exact_lp = Inf;
                        for ii=cter_a0+2:length(lp_a0)
                            for kk=1:length(exact_val)
                                tmp_bool1=exact_val(kk)<lp_a0(ii) && exact_val(kk)>lp_a0(ii-1);
                                tmp_bool2=exact_val(kk)>lp_a0(ii) && exact_val(kk)<lp_a0(ii-1);
                                % If the required exact value is in-between two values of the path, then save ii and kk and break
                                if tmp_bool1 || tmp_bool2
                                    exact_lp=kk;
                                    cter_a0=ii-1;
                                    sol = speicher{1,cter_a0+1};
                                    break
                                end
                            end
                            if exact_lp~=Inf
                                break
                            end
                        end
                        if exact_lp==Inf
                            if display
                                fprintf('\nAll the required exact solutions were computed and saved in p_exact.\n')
                            end
                            % Compute the next required value from here, which is not used in this run, since only_exact would lead straight to save&return
                            if speicher{1,end}.predictor.steplength >0
                                exact_lp = find(exact_val>=sol.lambda_p,1);
                                if isempty(exact_lp)
                                    exact_val=0;
                                    exact_lp=length(exact_val);
                                end
                            else
                                exact_lp = find(exact_val<=sol.lambda_p,1,'last');
                                if isempty(exact_lp)
                                    exact_val=0;
                                    exact_lp=1;
                                end
                            end
                            break
                        end
                    end
                end
            end
        end

        % When only_exact is activated, then stop here
        if only_exact==1
            jj=counter;
        else
            % Compute the predictor and corrector steps, if they are to be plotted
            if itnum~=-1 && isa(val,'cell')
                val{1,end+1} = a_p(end);
                val{2,end} = data(x1,a_p,ordnung,rho);
            end

            if log
                err = max(max(abs(sol.errest)));
                logstruct.max_error = [ logstruct.max_error err];
                logstruct.mesh_length = [logstruct.mesh_length [length(x1) ; length(xfin)]];
                logstruct.cpt = [logstruct.cpt cputime-cpt ] ;

                if itnum~=-1
                    logstruct.theta_0 = [ logstruct.theta_0 theta_0 ] ;
                    logstruct.theta_max = [ logstruct.theta_max theta_max ];
                    logstruct.steplength = [ logstruct.steplength steplength] ;

                    logstruct.norm_delta_0 = [logstruct.norm_delta_0 nda ] ;
                    logstruct.norm_delta_1 = [logstruct.norm_delta_1 nsd ] ;

                    if adapt
                        logstruct.sl_adapt = [ logstruct.sl_adapt 0];
                    else
                        logstruct.sl_adapt = [ logstruct.sl_adapt Inf];
                    end

                    if length(a_p)~=length(a_c) || max(abs(xfin-x1))>1e-12 % mesh got adapted
                        p_ma.lpfix = 1;
                        F_tmp1 = functionFDF( 'F', problem ,a_p ,x1,psival,psi,rho,p_ma) ;
                        F_tmp2 = functionFDF( 'F', problem ,a_c ,xfin,psival,psi,rho,p_ma) ;
                    else
                        F_tmp1 = functionFDF( 'F', problem ,a_p ,x1,psival,psi,rho,predictor) ;
                        F_tmp2 = functionFDF( 'F', problem ,a_c ,xfin,psival,psi,rho,predictor) ;
                    end
                    logstruct.norm_F = [logstruct.norm_F [ norm(F_tmp1) ; norm(F_tmp2) ] ] ;
                    logstruct.orthogonal = [logstruct.orthogonal [abs(F_tmp1(end)) ; abs(F_tmp2(end)) ] ] ;

                    if itnum~=0
                        logstruct.cos_ab = [logstruct.cos_ab cos_ab];
                    end
                end
            end

            % Compute the data for the plot of the mesh density evolution plot
            if dispmdep && itnum~=-1
                x_scaled = (xfin-xfin(1))/(xfin(end)-xfin(1)); % Scaling
                dens =1./(diff(x_scaled)*(length(x_scaled)-1)); % Mesh density
                x_int = x_scaled(1:end-1)+diff(x_scaled)/2;
                new_dens = pchip(x_int,dens,X_mdep);
                Y_mdep = [Y_mdep, itnum+1];
                Z_mdep = [Z_mdep; new_dens ];

                % Mesh evolution
                if length(xfin)<=pts_nb
                    filler = interval(2)*ones(1,length(xfin)-pts_nb);
                    Pts_Mat=[Pts_Mat; xfin filler];
                else
                    pts_step = ceil(length(xfin)/pts_nb);
                    filler = interval(2)*ones(1,pts_nb-length(xfin(1:pts_step:end-1)));
                    Pts_Mat = [Pts_Mat; xfin(1:pts_step:end-1) filler];
                end
            end

            % If the number of mesh points has been changed from the original number, then divide the number of mesh points by 2 if possible, otherwise reduce or keep the number of mesh points to minmeshpts
            % ..................... DISABLED FOR NOW .....................
            if counter_mc == 10^10 % Secondary option: Set a lower number
                if length(xfin)/2 > minmeshpts && mod(length(xfin),2)==1
                    xfin = xfin(1:2:end);
                elseif length(xfin)/2 > minmeshpts && mod(length(xfin),2)==0
                    xfin = [xfin(1:2:end),xfin(end)];
                elseif length(xfin)>minmeshpts
                    help = xfin(1:floor(length(xfin)/minmeshpts):end);
                    if help(end) == xfin(end)
                        xfin = help;
                    else
                        xfin = [help xfin(end)];
                    end
                end
                counter_mc = 0;
            end

            if halve_nb~=0
                fprintf('\nThe step-length was halved %i times in this step.', halve_nb)
            end

            if itnum==-1
                fprintf('\n<strong>Initial solution computed</strong> for value %f.\n', a_c(end) )
            elseif isempty(max_pred)
                fprintf('\n<strong>Step %i complete</strong> for value %f. (Step %i/%i)\n', itnum+1, a_c(end), jj, counter )
            else
                fprintf('\n<strong>Step %i complete</strong> for value %f with pl=%1.3f and max_pl=%1.3f. (Step %i/%i)\n', itnum+1, a_c(end), pred_tmp, max_pred, jj, counter )
            end

            % If min_sl was reached, then the user can choose whether to continue or not
            if min_sl_flag
                min_sl_flag = 0;
                halve = 0;
                jj=counter;
                if display
                    fprintf('\nThe minimal steplength %1.1e was reached!\nChoose whether to continue or stop here.\n',min_sl)
                end
            end

            % Check whether the pit_stop values have been reached
            if itnum~=-1 && ~isempty(pit_stop)
                tmp_bool1= a_0(end)<pit_stop(1) & pit_stop(1)<=a_c(end);
                tmp_bool2= a_0(end)>pit_stop(1) & pit_stop(1)>=a_c(end);
                if tmp_bool1 || tmp_bool2
                    jj=counter;
                    if display
                        fprintf('\nParameter value %f reached!!\nTime for a <strong>pit stop</strong>.\n',pit_stop(1))
                    end
                    pit_stop(1)=Inf;
                end
                tmp_bool1= speicher{3,end-1}<pit_stop(2) & pit_stop(2)<=speicher{3,end};
                tmp_bool2= speicher{3,end-1}>pit_stop(2) & pit_stop(2)>=speicher{3,end};
                if max(tmp_bool1) || max(tmp_bool2)
                    jj=counter;
                    if display
                        fprintf('\nData value %f reached!!\nTime for a <strong>pit stop</strong>.\n',pit_stop(2))
                    end
                    pit_stop(2)=Inf;
                end
            end

            % Save the values to jump back to this state from a later step
            if itnum>0
                savestruct = struct;

                savestruct.delta_0     = delta_0;
                savestruct.theta_0     = theta_0;
                savestruct.itnum       = itnum;
                savestruct.val         = val;
                savestruct.exact_lp    = exact_lp;
                savestruct.pred_tmp    = pred_tmp;
                savestruct.corr_old    = corr_old;
                savestruct.pit_stop    = pit_stop;
                savestruct.halve_nb    = halve_nb;
                if dispmdep
                    savestruct.Pts_Mat = Pts_Mat;
                    savestruct.y_mdep  = Y_mdep;
                    savestruct.z_mdep  = Z_mdep;
                end

                jump_cell{end+1}=savestruct;
            end
        end
    end

    % when the number of steps has been carried out
    if jj==counter
        jj=0;

        % Call the plot
        clf(pathfollplot) ;
        pathfollplot.Color = 'White' ;

        % Draw the path
        if dispres
            % Values of lambda_p
            X_se = cell2mat(speicher(2,:));

            % Mesh points
            if infsplit % Mesh was transformed from [0,infty) to [0,1]
                % In [0,1], only plot every 'DD'th value. This number is a divisor of se_nn.
                KK = 1:30;
                DD = KK(rem(se_nn,KK)==0);
                DD = min(DD(DD>2));
                Y_se = [ linspace(0,1,se_nn/DD+1), fliplr(1./linspace(1/se_nn,1-1/se_nn,se_nn-1)) ];
                selimpt = find(Y_se >= selim,1);
                Y_se = Y_se(1:selimpt);
            else
                Y_se = linspace(interval(1),interval(2),se_nn);
            end

            % Function values at X_se and Y_se
            Z_se = cell(1,n);
            for tt=1:n
                Z_se{tt} = zeros(length(X_se),length(Y_se));
                if infsplit
                    for kk=1:itnum+2
                        x1tmp = speicher{1,kk}.x1;
                        coefftmp = speicher{1,kk}.coeff;
                        help = coeffToValues( coefftmp, x1tmp,ordnung,rho,linspace(0,1,se_nn+1));
                        help1 = help(tt,1:DD:end); % values on [0,1]
                        help2 = help(n+tt,end-1:-1:size(help,2)+se_nn/DD+1-selimpt); % values on (1,selim]
                        Z_se{tt}(kk,:)=[help1,help2];
                    end
                else
                    for kk=1:itnum+2
                        help = coeffToValues( speicher{1,kk}.coeff, speicher{1,kk}.x1,ordnung,rho,Y_se);
                        Z_se{tt}(kk,:) = help(tt,:);
                    end
                end

            end

            % Draw the plots
            if disppredcor && itnum~=-1
                figdataplot(itnum+2,speicher,X_se,Y_se,Z_se,pathfollplot,path_pos,val,problem)
            elseif dispres && itnum~=-1
                figdataplot(itnum+2,speicher,X_se,Y_se,Z_se,pathfollplot,path_pos,zeros(length(speicher{3,1}),1),problem)
            end
        end

        if log && itnum~=-1 && displog % Draw log plot
            if only_exact
                p_ma.lpfix = 1;
                F_tmp2 = functionFDF( 'F', problem ,a_c ,xfin,psival,psi,rho,p_ma) ;
            end
            figlogdataplot(itnum+1,logstruct,pathfollplot,log_pos,AbsTol+max(F_tmp2)*RelTol)
        end

        if dispmdep && size(Z_mdep,1)>1 % Draw mesh and mesh density evolution plot
            figure(pathfollplot);

            % Mesh density evolution plot
            axes( 'Position', [ mdep_pos(1)+mdep_pos(3)/10, mdep_pos(2)+mdep_pos(4)/10, 5*mdep_pos(3)/10, 7*mdep_pos(4)/10] ) ;
            surf(X_mdep,Y_mdep,Z_mdep)
            zl = zlim;
            zlim([0 zl(2)])
            title('mesh density')
            ax = gca;
            ax.FontSize=f_size;
            ax.Title.FontWeight='normal';

            % Mesh evolution plot
            axes( 'Position', [ mdep_pos(1)+7*mdep_pos(3)/10, mdep_pos(2)+mdep_pos(4)/10, 2*mdep_pos(3)/10, 7*mdep_pos(4)/10] ) ;
            hold on
            for ii = Y_mdep(1)-1:itnum+1
                % Rescale in case meshadaptation was not used from the start but just from a loaded point
                oo=ii-Y_mdep(1)+1;
                plot(Pts_Mat(oo+1,:),ii,'k.','LineWidth',1.5);
            end
            title('mesh evolution')
            ax = gca;
            ax.FontSize=f_size;
            ax.Title.FontWeight='normal';
        end

        % If enabled, go back a step if wished
        if itnum>1 && save_ws && ~cter_flag && ~only_exact
            x_inp = input(prompt_ws,'s');
            switch x_inp
                case 'b' % Go back a number of steps
                    % Let the user input a number of steps to go back
                    while 1
                        y_inp=input(subprompt3);
                        if y_inp>=0 && y_inp>=0 && y_inp<=length(jump_cell)-2 && rem(y_inp,floor(y_inp))==0
                            back_cter=y_inp;
                            break
                        end
                        fprintf('Input should be a natural number between 0 and %i.\n',length(jump_cell)-2)
                    end

                    % Adjust jump_cell and define savestruct
                    jump_cell = jump_cell(1:end-back_cter);
                    savestruct = jump_cell{end};

                    % Load from savestruct
                    delta_0         = savestruct.delta_0;
                    theta_0         = savestruct.theta_0;
                    itnum           = savestruct.itnum;
                    val             = savestruct.val;
                    exact_lp        = savestruct.exact_lp;
                    pred_tmp        = savestruct.pred_tmp;
                    corr_old        = savestruct.corr_old;
                    pit_stop        = savestruct.pit_stop;
                    if isfield(savestruct,'halve_nb')
                        halve_nb    = savestruct.halve_nb;
                    else
                        halve_nb    = 5;
                    end
                    if dispmdep
                        try
                            Pts_Mat = savestruct.Pts_Mat;
                            Y_mdep  = savestruct.y_mdep;
                            Z_mdep  = savestruct.z_mdep;
                        catch
                            Y_mdep  = [];
                        end
                    end

                    % Adjust speicher
                    speicher = speicher(:,1:end-back_cter);

                    % Define x1 and a_p
                    sol       = speicher{1,end-1}; % sol is the penultimate solution
                    x1        = sol.x1;
                    a_0       = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
                    % Define sol, predictor, xfin and a_c
                    sol       = speicher{1,end}; % sol is the last computed solution
                    predictor = sol.predictor;
                    a_p       = a_0 + predictor.steplength.*predictor.tangent';
                    xfin      = sol.x1;
                    a_c       = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;

                    % Readjust all the fields in logstruct
                    if log && back_cter~=0
                        tmp_fn = fieldnames(logstruct);
                        for kk=1:numel(tmp_fn)
                            logstruct.(tmp_fn{kk})=logstruct.(tmp_fn{kk})(:,1:end-back_cter);
                        end
                    end

                    predictor.maxCorrSteps        = maxCorrSteps;
                    predictor.maxSteplengthGrowth = maxSteplengthGrowth;
                    predictor.meshFactorMax       = meshFactorMax;

                    % Load steplength and tangent_new from predictor
                    steplength = predictor.steplength;
                    tangent_new = predictor.tangent';
                    % Compute eval_exact if necessary
                    if ~isempty(exact_val)
                        tmp_bool1= exact_lp==1 && exact_val(exact_lp)<sol.lambda_p && steplength>0;
                        tmp_bool2= exact_lp==length(exact_val) && exact_val(exact_lp)>sol.lambda_p && steplength<0;
                        if tmp_bool1 || tmp_bool2
                            eval_exact=1;
                        else
                            eval_exact=0;
                        end
                    end
            end
        end

        % If only_exact is activated or the counter was set to Inf, then speicher and speicher_exact are saved, otherwise it depends on the user input
        if only_exact || cter_flag
            x_inp = 's';
        else
            x_inp = input(prompt,'s');
            while x_inp == 'f'
                while 1
                    y_inp=input(subprompt4);
                    if ~isempty(y_inp) && min(y_inp>=0) && min(y_inp<=itnum+1) && max(rem(y_inp+1,floor(y_inp+1)))==0
                        break
                    end
                    fprintf('Inputs must be natural numbers between 0 and %i!\n',itnum+1)
                end
                for ii=y_inp
                    plot_step(speicher{1,ii+1},f_size,infsplit,n,selim)
                end
                x_inp=input(subprompt5,'s');
            end
        end
        switch x_inp
            case 'p' % Change maximal steplength
                while 1
                    y_inp=input(subprompt2);
                    if ~isempty(y_inp) && y_inp>0
                        break
                    end
                    fprintf('Input must be a positive number!\n')
                end
                max_pred=y_inp;
            case 'n' % Change number of steps to carry out at once
                while 1
                    y_inp=input(subprompt1);
                    if ~isempty(y_inp) && y_inp>0 && rem(y_inp,floor(y_inp))==0
                        break
                    end
                    fprintf('Input must be a positive natural number!\n')
                end
                counter=y_inp;
            case 's' % Save & stop
                if ~only_exact
                    % Save
                    savestruct = struct;

                    savestruct.delta_0     = delta_0;
                    savestruct.theta_0     = theta_0;
                    savestruct.itnum       = itnum;
                    savestruct.val         = val;
                    savestruct.exact_lp    = exact_lp;
                    savestruct.pred_tmp    = pred_tmp;
                    savestruct.corr_old    = corr_old;
                    savestruct.pit_stop    = pit_stop;
                    savestruct.halve_nb    = halve_nb;
                    if dispmdep
                        savestruct.Pts_Mat = Pts_Mat;
                        savestruct.y_mdep  = Y_mdep;
                        savestruct.z_mdep  = Z_mdep;
                    end
                end

                % savestruct is saved at the end of speicher, when the program is started from this point, it will access the last available data in this way
                speicher{1,end+1} = savestruct;
                if log
                    speicher{2,end} = logstruct;
                end
                if save_ws
                    speicher{3,end} = jump_cell;
                end

                lp_a0 = cell2mat(speicher(2,1:end-1));
                % Find the turning points in lp_a0
                tmp1 = find( lp_a0(1:end-1)>lp_a0(2:end) );
                tmp2 = find( lp_a0(1:end-1)<lp_a0(2:end) );
                if isempty(tmp1) || isempty(tmp2)
                    tur_pts=[];
                else
                    tur_pts=zeros(3,length(lp_a0));
                    kk=0;
                    for ii=2:length(lp_a0)
                        if (max(tmp1==ii) && max(tmp2==ii-1)) ||  (max(tmp2==ii) && max(tmp1==ii-1))
                            kk=kk+1;
                            tur_pts(1,kk)=lp_a0(ii-1);
                            tur_pts(2,kk)=lp_a0(ii);
                            tur_pts(3,kk)=lp_a0(ii+1);
                        end
                    end
                    tur_pts=tur_pts(:,1:kk);
                end

                % Reinstate the previous interpreter settings
                set(0, 'DefaultTextInterpreter', dtiii)
                set(0, 'DefaultLegendInterpreter', dliii)
                set(groot, 'DefaultAxesTickLabelInterpreter', datliii)

                return
        end
    end

    if skip==0 % If the step-length needs to be predicted
        % Prepare for the next step
        [a_0,a_p,tangent_new,steplength,predictor,delta_0,nda,nsd,theta_0,logstruct,eval_exact,exact_lp] = PredCorrStrat(problem,settings,predictor,a_c,a_p,x1,xfin,sol,tangent_new,ordnung,rho,psival,psi,steplength,exact_val,exact_lp,delta_0,theta_0,theta_max,log,logstruct,eval_exact,itnum,halve_nb,min_sl);

        % Accept the new mesh for the next step
        x1 = xfin;
    end

    itnum = itnum+1;
    jj = jj+1;

end

end

%% Local functions
function [tangente]=tangente_berechnen(DF)

help1=zeros(length(DF(:,1)),1);
help1(end)=1;

tangente=DF\help1;

tangente=tangente/sqrt(sum(tangente.*tangente));

end

function Psireturn=pathfollowing_Psi(n,rho,nr)
i=length(rho);
prod=1;
for s=1:i
    if (s~=n)
        prod=conv(prod,[1 -rho(s)])/(rho(n)-rho(s));
    end
end
for s=1:(nr)
    prod=polyint(prod);
end
Psireturn=prod;
%tabulars for collocation points
end

function [] = figdataplot(ii,speicher,X_se,Y_se,Z_se,pathfollplot,path_pos,val,problem)

blue_rgb = [0 0.4470 0.7410];
red_rgb = [0.8500 0.3250 0.0980];
green_rgb = [0.4660 0.6740 0.1880];

figure(pathfollplot);
if isa(val,'cell')
    val1=cell2mat(val(1,:));
    val2=cell2mat(val(2,:));
else
    val2=val;
end
nn=size(val2,1);
n=length(Z_se);

f_size=8; % font size for plot

% data evolution plots
for jj=1:nn
    axes( 'Position', [ path_pos(1)+2*path_pos(3)/20, path_pos(2)+(1+(nn-jj)*10)*path_pos(4)/(nn*10), 7*path_pos(3)/20, 8*path_pos(4)/(nn*10)] ) ;
    hold on
    if ~isempty(val)
        for kk=1:ii-1
            if strcmp(problem,'testK_8')
                plot(1./[speicher{2,kk},val1(kk)],[speicher{jj+2,kk},val2(jj,kk)],'Color',green_rgb,'LineStyle','--') ;
                plot(1./[val1(kk),speicher{2,kk+1},],[val2(jj,kk),speicher{jj+2,kk+1}],'Color',red_rgb,'LineStyle',':','Marker','o','MarkerSize',6) ;
            else
                plot([speicher{2,kk},val1(kk)],[speicher{3,kk}(jj),val2(jj,kk)],'Color',green_rgb,'LineStyle','--') ;
                plot([val1(kk),speicher{2,kk+1},],[val2(jj,kk),speicher{3,kk+1}(jj)],'Color',red_rgb,'LineStyle',':','Marker','o','MarkerSize',6) ;
            end
        end
    end
    if strcmp(problem,'testK_8')
        a=speicher(2,1:ii);c=1./cell2mat(a);
    else
        a=speicher(2,1:ii);c=cell2mat(a);
    end
    b=speicher(3,1:ii);d=cell2mat(b);e=d(jj,:);
    plot(c,e,'Color',blue_rgb,'LineStyle','-','Marker','o','MarkerFaceColor',blue_rgb,'MarkerSize',3,'LineWidth',1.5)
    xlabel('pathfollowing parameter')
    if ii>50
        plot(c(50:50:ii),e(50:50:ii),'xw','MarkerSize',5)
    end
    ax = gca;
    ax.FontSize=f_size;
end

% solution evolution plot
if strcmp(problem,'testM_7')
    xstop=8;
    axes( 'Position', [ path_pos(1)+12*path_pos(3)/20, path_pos(2)+11*log_pos(4)/20, 7*path_pos(3)/20, 8*path_pos(4)/20] ) ;
    subplot(n*nn,2,2:2:2*3*nn)
    surf( [Y_se fliplr(1./Y_se(xstop:end-1))] ,X_se, [Z_se{1} fliplr(Z_se{3}(:,xstop:end-1))] )
    title('real part of the solution')
    xlabel('interval on mesh')
    ylabel('pathfollowing parameter')
    zlabel('function value')
    ax = gca;
    ax.FontSize=f_size;

    axes( 'Position', [ path_pos(1)+12*path_pos(3)/20, path_pos(2)+1*log_pos(4)/20, 7*path_pos(3)/20, 8*path_pos(4)/20] ) ;
    subplot(n*nn,2,2*3*nn+2:2:2*6*nn)
    surf( [Y_se fliplr(1./Y_se(xstop:end-1))] ,X_se, [Z_se{2} fliplr(Z_se{4}(:,xstop:end-1))] )
    title('real part of the solution')
    xlabel('interval on mesh')
    ylabel('pathfollowing parameter')
    zlabel('function value')
    ax = gca;
    ax.FontSize=f_size;
else
    for kk=1:n
        axes( 'Position', [ path_pos(1)+12*path_pos(3)/20, path_pos(2)+(1+(n-kk)*10)*path_pos(4)/(n*10), 7*path_pos(3)/20, 8*path_pos(4)/(n*10)] ) ;
        if strcmp(problem,'testK_8')
            surf(Y_se,1./X_se,Z_se{kk})
        else
            surf(Y_se,X_se,Z_se{kk})
        end
        name = strcat('solution function',num2str(kk));
        title(name)
        xlabel('interval on mesh')
        ylabel('pathfollowing parameter')
        zlabel('function value')
        ax = gca;
        ax.FontSize=f_size;
    end
end

end

function [] = figlogdataplot(ii,logstruct,pathfollplot,log_pos,tol)
nrow = 5;
ncol = 2;
f_size=8;
figure(pathfollplot);

axes( 'Position', [ log_pos(1)+2*log_pos(3)/20, log_pos(2)+42*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(0:ii,logstruct.cpt)
hold on ;
plot(1:ii,logstruct.trm_halved(1:ii),'*r')
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)+1]) ;
title('cputime')
sp(1) = gca;
sp(1).FontSize=f_size;
sp(1).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+12*log_pos(3)/20, log_pos(2)+42*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(1:ii,logstruct.steplength)
hold on ;
plot(1:ii,logstruct.steplength_pred,'k')
plot(1:ii,zeros(1,ii),'k:')
plot(1:ii,logstruct.sl_adapt(1:ii),'*g')
set(gca,'Ylim',[-5/4*max(abs(logstruct.steplength)) 5/4*max(abs(logstruct.steplength))]) ;
legend('$s$','predicted $s$')
title('steplength $s$')
sp(2) = gca;
sp(2).FontSize=f_size;
sp(2).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+2*log_pos(3)/20, log_pos(2)+32*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(1:ii,logstruct.theta_0)
hold on
plot(1:ii,ones(1,ii),'k--')
plot(1:ii,logstruct.theta_max,'r--')
%limsy=get(gca,'YLim');
set(gca,'Ylim',[0 3/2*max(abs(logstruct.theta_0))]);%max(1.1,limsy(2))]) ;
title('$\theta_0$')
sp(3) = gca;
sp(3).FontSize=f_size;
sp(3).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+12*log_pos(3)/20, log_pos(2)+32*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
semilogy(1:ii,logstruct.theta_max./logstruct.theta_0)
title('$\Theta_m / \Theta_0$')
sp(4) = gca;
sp(4).FontSize=f_size;
sp(4).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+2*log_pos(3)/20, log_pos(2)+22*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(0:ii,logstruct.mesh_length(1,:))
hold on
plot(0:ii,logstruct.mesh_length(2,:),'k')
plot(1:ii,logstruct.mesha_halved(1:ii),'*r')
legend('at beginning','at end')
set(gca,'Ylim',[0 3/2*max(abs(logstruct.mesh_length(2,:)))]);
title('\# mesh points in step')
sp(5) = gca;
sp(5).FontSize=f_size;
sp(5).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+12*log_pos(3)/20, log_pos(2)+22*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
semilogy(1:ii,logstruct.norm_F(1,:))
hold on
semilogy(1:ii,logstruct.norm_F(2,:))
semilogy(1:ii,tol.*ones(1,ii),'k--')
legend('$||F(a_p)||$','$||F(a_c)||$')
title('$||F||$')
sp(6) = gca;
sp(6).FontSize=f_size;
sp(6).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+2*log_pos(3)/20, log_pos(2)+12*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
semilogy(1:ii,logstruct.orthogonal(1,:))
hold on
semilogy(1:ii,logstruct.orthogonal(2,:))
legend('$|F(a_p)(end)|$','$|F(a_c)(end)|$')
title('$|F(end)|$')
sp(7) = gca;
sp(7).FontSize=f_size;
sp(7).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+12*log_pos(3)/20, log_pos(2)+12*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
semilogy(1:ii,logstruct.norm_delta_0,'--')
hold on
semilogy(1:ii,logstruct.norm_delta_1,'--')
semilogy(1:ii,logstruct.corr_dist,'--')
semilogy(1:ii,logstruct.pdist_halved(1:ii),'*r')
semilogy(1:ii,logstruct.cdist_halved(1:ii),'*g')
legend('$||\Delta y_0||$','$||\Delta y_1||$','$||a_p-a_c||$')
title('$||...||$')
sp(8) = gca;
sp(8).FontSize=f_size;
sp(8).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+2*log_pos(3)/20, log_pos(2)+2*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(1:ii,logstruct.c_s)
hold on
plot(1:ii,logstruct.cos_ab,'k')
plot(1:ii,logstruct.angle_halved(1:ii),'*r')
legend('$c_s$','cos($ab$)')
title('$c_s$ and cos($ab$)')
sp(9) = gca;
sp(9).FontSize=f_size;
sp(9).Title.FontWeight='normal';

axes( 'Position', [ log_pos(1)+12*log_pos(3)/20, log_pos(2)+2*log_pos(4)/50, 7*log_pos(3)/20, 7*log_pos(4)/50] ) ;
plot(0:ii,logstruct.max_error)
title('Maximal Error')
sp(10) = gca;
sp(10).FontSize=f_size;
sp(10).Title.FontWeight='normal';

for jj=1:nrow*ncol
    xlim(sp(jj),[0 ii+1])
end

end

function [a_0,a_p,tangent_new,steplength,predictor,delta_0,nda,nsd,theta_0,logstruct,eval_exact,exact_lp] = PredCorrStrat(problem,settings,predictor,a_c,a_p,x1,xfin,sol,tangent_new,ordnung,rho,psival,psi,steplength,exact_val,exact_lp,delta_0,theta_0,theta_max,log,logstruct,eval_exact,itnum,halve_nb,min_sl)

% In the case of mesh adaptation enabled, take special care to set the new starting point of the predictor step
if itnum~=-1 && (feval(settings,'meshAdaptation'))
    % If the mesh was adapted, then compute the tangent and a_0 for the new number of mesh points
    if length(x1)~=length(xfin) || max(abs(x1-xfin))>1e-12
        initP.initialMesh=xfin;
        initP.parameters = sol.parameters;

        initP.initialValues = coeffToValues(tangent_new, x1,ordnung,rho,xfin);
        initP = initial_coefficients(problem,xfin,initP,rho,0);
        tangent_new = [ initP.initialCoeff ; tangent_new(end) ] ;
        tangent_new(1:end-1) = tangent_new(1:end-1)/sqrt(sum(tangent_new.'*tangent_new)) ;

        initP.initialValues = coeffToValues(a_p, x1,ordnung,rho,xfin);
        initP = initial_coefficients(problem,xfin,initP,rho,0);
        a_p = [ initP.initialCoeff ; a_p(end) ] ;

        a_0 = [ sol.coeff ; sol.parameters ; a_c(end) ] ;

        predictor.x1 = xfin;
    else
        a_0 = a_c;
    end
else
    a_0 = a_c;
end

if itnum~=-1
    tangent_old = tangent_new;
end

% Compute the tangent at the new point
jac_F = functionFDF( 'DF', problem ,a_0,xfin,psival,psi,rho,[]);
tangent_new = tangente_berechnen(jac_F);

% Check for turning point
if itnum~=-1 && sign(tangent_old.'*tangent_new)<0
    fprintf('\n\n			 <strong>***</strong> Found <strong>TURNING POINT</strong> near %f! <strong>***</strong> \n\n',a_0(end));
    steplength=-steplength;
    % take the next value of the required values for next time
    if ~isempty(exact_val) && eval_exact==1
        if (exact_lp>1 && steplength<0) || (exact_lp<length(exact_val) && steplength>0)
            exact_lp = exact_lp + sign(steplength);
        else
            eval_exact=0;
        end
    elseif ~isempty(exact_val)
        eval_exact = 1;
    end
end

if itnum~=-1
    % Compute the new steplength
    c_s = abs(tangent_new.'*tangent_old) ;
    corr_dist = norm(a_p-a_0) ;
    tmp=2*norm(delta_0)/(c_s^2*corr_dist) ;
    beta = sqrt((theta_max/theta_0)*tmp);
    beta_max = max(min(1,predictor.maxSteplengthGrowth),predictor.maxSteplengthGrowth/(2^(0.5*halve_nb)));
    beta=min(beta,beta_max);
    steplength = beta*steplength;

    if predictor.display
        fprintf('\n  Steplength was %3.1e and is predicted to be %3.1e!\n',predictor.steplength,steplength);
    end
else
    c_s = 1;
end

if log && itnum~=-1
    logstruct.c_s = [ logstruct.c_s c_s ] ;

    logstruct.corr_dist = [logstruct.corr_dist corr_dist ] ;

    logstruct.steplength_pred = [logstruct.steplength_pred steplength];
end

% Update predictor
predictor.a_0 = a_0.';
predictor.tangent = tangent_new.';
predictor.steplength = steplength;
predictor.lambda_p_0 = a_0(end);

correct=0;
countercorr=0;
ncorrsteps=predictor.maxCorrSteps;
theta_0_min=Inf; steplength_min=steplength;
while correct~=1
    if abs(countercorr)==ncorrsteps
        steplength=steplength_min;
        predictor.steplength = steplength;
    end

    % Carry out predictor step
    a_p = a_0 + steplength.*tangent_new;

    % Update the value of the pathfollowing variable
    predictor.lambda_p_p = a_p(end);

    % Compute what will be needed to compute the new steplength after the corrector step
    F_a = functionFDF( 'F', problem ,a_p,xfin,psival,psi,rho,predictor);
    jac_F_a = functionFDF( 'DF', problem ,a_p,xfin,psival,psi,rho,predictor);
    delta_a = - jac_F_a\F_a ;
    new_a = a_p + delta_a ;
    F_new = functionFDF( 'F', problem ,new_a ,xfin,psival,psi,rho,predictor);
    simplified_delta = - jac_F_a\F_new ;

    nda = norm(delta_a);
    nsd = norm(simplified_delta);

    theta_0 = nsd / nda ;
    if theta_0<theta_0_min
        steplength_min=steplength;
    end

    delta_0=delta_a;
    beta = sqrt(1/2);%sqrt(theta_max/theta_0)
    predictor.steplength = steplength;

    % No Prediction-correction is needed if steplength is smaller than the minimally allowed step-length
    if abs(steplength)<min_sl
        break
    end

    % Prediction-correction is needed if theta_0 is too big
    if theta_0>theta_max*c_s^2 && abs(countercorr)~=ncorrsteps
        sl_tmp=steplength;
        steplength=beta*steplength;
        if predictor.display && countercorr==0
            fprintf('  <strong>Prediction correction:</strong>\n  Steplength is decreased from %1.1e to %1.1e',sl_tmp,steplength);
        elseif predictor.display
            fprintf(' to %1.1e',steplength)
        end
        countercorr=countercorr-1;
        % Procedure to increase steplength if needed -- disabled for now
        % elseif theta_0<theta_max*c0^2/4 && abs(countercorr)~=ncorrsteps
        %  sl_tmp=steplength;
        %  steplength=steplength/beta*0.99;
        %  countercorr=countercorr+1;
        %  fprintf('\n Steplength is increased from %f to %f!\n',sl_tmp,steplength);
    else
        if predictor.display && abs(countercorr)~=ncorrsteps
            fprintf('\n  Correction procedure should converge with steplength %3.1e!...\n',steplength);
        elseif predictor.display
            fprintf('\n  The steplength %f is used!...\n',steplength);
        end
        correct=1;
    end
end

end

function [speicher,speicher_exact,val,corr_old] = correct_data(speicher,speicher_exact,val,data,ordnung,rho)
sol = speicher{1,1};
a_0 = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
x1 = sol.x1;
speicher{3,1} = data(x1,a_0,ordnung,rho);
for ii=2:size(speicher,2)
    predictor = speicher{1,ii}.predictor;
    a_p = a_0 + predictor.steplength.*predictor.tangent';
    if iscell(val)
        val{2,ii-1} = data(x1,a_p,ordnung,rho);
    end
    sol = speicher{1,ii};
    a_0 = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
    x1 = sol.x1;
    speicher{3,ii} = data(x1,a_0,ordnung,rho);
end

data_tmp = data(speicher{1,end-1}.x1,a_p,ordnung,rho);
corr_old = sqrt((a_0(end)-a_p(end))^2+(max(abs(speicher{3,end}-data_tmp)))^2);

for ii=1:size(speicher_exact,2)
    sol = speicher_exact{1,ii};
    a_0 = [ sol.coeff ; sol.parameters ; sol.lambda_p ] ;
    x1 = sol.x1;
    speicher_exact{3,ii} = data(x1,a_0,ordnung,rho);
end

end

function [] = plot_step(sol,f_size,infsplit,n,selim)
exact_plot=figure('units','normalized','outerposition',[0.3 0.2 0.4 0.6],'PaperUnits','normalized');
ax = gca;
ax.FontSize=f_size;
clf ;
exact_plot.Color = 'White' ;
if infsplit
    selimpt = find((1./sol.x1)>=selim,1,'last');
end
for ii = 1:n
    axes( 'Position', [ 1/10, (1+(n-ii)*10)/(n*10), 8/10, 8/(n*10)] ) ;
    if infsplit
        plot([sol.x1 fliplr(1./sol.x1(selimpt:end-1))],[sol.valx1(ii,:) fliplr(sol.valx1(n+ii,selimpt:end-1))],'LineWidth',1.5)
    else
        plot(sol.x1,sol.valx1(ii,:),'LineWidth',1.5)
    end
    plottitle=['solution ',num2str(ii),' for $\lambda_p=$',num2str(sol.lambda_p)];
    title(plottitle)
    ax = gca;
    ax.FontSize=f_size;
    ax.Title.FontWeight='normal';
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  dispDebug.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = dispDebug(message,show)
if(show>0)
   disp(message);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bvpsuite2.0 file:  template_bvp.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ret] = template_bvp(request,z,za,zb,zc,t,p,lambda,lambda_p)
switch request
    case 'n'
        ret=1;
    case 'orders'
        ret=[3];
    case 'problem'
        %Problemstellung als System der F(z(i,j))=0 (jede Zeile entspricht einer Gleichung) mit z(i,j), wobei
        %i..Loesungskomponente
        %j..Ableitung
        %also z(i,j)=z_i^(j)(t)
        ret=[z(1,1)+z(1,3)-z(1,2)-z(1,4)];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1),0)),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        %ret(i,j,k): Ableitung von
        %i-ter Gleichung nach z(j,k)
         ret(1,1,1)=1;
         ret(1,1,2)=-1;
         ret(1,1,3)=1;
         ret(1,1,4)=-1;
    case 'interval'
        ret = [-1,1];
    case 'linear'
        ret=1;
    case 'parameters'
        ret=0;
    case 'c'
        ret = [];
    case 'BV'
        ret=[za(1,1)-exp(-1);
             za(1,2)-exp(-1);
             za(1,3)-exp(-1)];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1),0)),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        %ret(cInd,i,j,k)
        %cInd..spezifiziert Punkt, an dem Bedingung gestellt wird. (falls c=[] entspricht
        %1..a und 2..b, ansonsten Index von c)
        %i,j,k wie bei jacobian
         ret(1,1,1,1)=1;
         ret(1,2,1,2)=1;
         ret(1,3,1,3)=1;
    case 'dP'
        ret=[0];
    case 'dP_BV'
        ret=[0;0;0];
    case 'initProfile'
        ret.initialMesh = linspace(-1,1,20);
        ret.initialValues = linspace(-1,1,20).^2;
    case 'EVP'
        ret=0;
    case 'dLambda'
        ret = 0;
    case 'pathfollowing'
        ret.activate=0;
    case 'path_jac'
        ret=zeros(feval(mfilename,'n'),1);
    case 'path_dBV'
        ret=zeros(sum(feval(mfilename,'orders'))+feval(mfilename,'parameters'),1);
    otherwise
        ret = 0;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  access to the internal routines (used by plug-ins, e.g. the DBN mesh controller)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = solver_api()
S = struct('solveLinearProblem',@solveLinearProblem, ...
           'solveNonLinearProblem',@solveNonLinearProblem, ...
           'errorestimate',@errorestimate, ...
           'computeResidual',@computeResidual, ...
           'feval_problem',@feval_problem, ...
           'getStandardCollocationPoints',@getStandardCollocationPoints);
end
