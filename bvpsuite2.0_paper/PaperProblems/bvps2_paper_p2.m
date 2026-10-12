function [ret] = bvps2_paper_p2(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Problem 2 (Shampine's problem, refs [34,35] of the paper):
%   y'' + 2/t y' - phi^2 y exp( gamma*beta*(1-y)/(1+beta*(1-y)) ) = 0,
%   y'(0) = 0,  y(1) = 1,    phi = 0.6, gamma = 40, beta = 0.2.
% Three solutions exist; the initial profile y = 0 gives y_1, y = 1 gives y_3.
% The initial profile (constant value) is taken from the global PAPER.p2_init.

global PAPER
phi = 0.6;  gamma = 40;  beta = 0.2;

switch request
    case 'n'
        ret = 1;
    case 'orders'
        ret = 2;
    case 'problem'
        y = z(1,1);
        E = gamma*beta*(1-y)/(1+beta*(1-y));
        ret = [ z(1,3) + 2/t*z(1,2) - phi^2*y*exp(E) ];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        y  = z(1,1);
        D  = 1+beta*(1-y);
        E  = gamma*beta*(1-y)/D;
        dE = -gamma*beta/D^2;
        ret(1,1,1) = -phi^2*exp(E)*(1 + y*dE);
        ret(1,1,2) = 2/t;
        ret(1,1,3) = 1;
    case 'interval'
        ret = [0,1];
    case 'linear'
        ret = 0;
    case 'parameters'
        ret = 0;
    case 'c'
        ret = [];
    case 'BV'
        ret = [ za(1,2);
                zb(1,1)-1 ];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        ret(1,1,1,2) = 1;
        ret(2,2,1,1) = 1;
    case 'dP'
        ret = [];
    case 'dP_BV'
        ret = [];
    case 'initProfile'
        x = linspace(0,1,51);
        ret.initialMesh   = x;
        ret.initialValues = PAPER.p2_init*ones(1,51);
    case 'EVP'
        ret = 0;
    case 'dLambda'
        ret = 0;
    case 'pathfollowing'
        ret.activate = 0;
    otherwise
        ret = 0;
end

end
