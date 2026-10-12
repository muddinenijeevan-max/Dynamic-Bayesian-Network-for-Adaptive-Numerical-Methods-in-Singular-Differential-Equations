function [ret] = bvps2_paper_p3(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Problem 3 (electromagnetic self-interaction theory, ref [7] of the paper):
%   y'' + 4/t y' + (t y - 1) y = 0,   y'(0) = 0,
%   y(L) + y'(L) = 0   (PAPER.p3_bc = 1)   or   y(L) = 0   (PAPER.p3_bc = 2),
% on [0,L], L = 5, 20.  Initial profile:  y = 2 (t<=1.5),  y = 2 exp(1.5-t) (t>1.5).

global PAPER

switch request
    case 'n'
        ret = 1;
    case 'orders'
        ret = 2;
    case 'problem'
        ret = [ z(1,3) + 4/t*z(1,2) + (t*z(1,1) - 1)*z(1,1) ];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        ret(1,1,1) = 2*t*z(1,1) - 1;
        ret(1,1,2) = 4/t;
        ret(1,1,3) = 1;
    case 'interval'
        ret = [0, PAPER.p3_L];
    case 'linear'
        ret = 0;
    case 'parameters'
        ret = 0;
    case 'c'
        ret = [];
    case 'BV'
        if PAPER.p3_bc == 1
            ret = [ za(1,2);
                    zb(1,1) + zb(1,2) ];
        else
            ret = [ za(1,2);
                    zb(1,1) ];
        end
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        ret(1,1,1,2) = 1;
        ret(2,2,1,1) = 1;
        if PAPER.p3_bc == 1
            ret(2,2,1,2) = 1;
        end
    case 'dP'
        ret = [];
    case 'dP_BV'
        ret = [];
    case 'initProfile'
        x = linspace(0,PAPER.p3_L,51);
        ret.initialMesh   = x;
        ret.initialValues = (x<=1.5).*2 + (x>1.5).*2.*exp(1.5-x);
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
