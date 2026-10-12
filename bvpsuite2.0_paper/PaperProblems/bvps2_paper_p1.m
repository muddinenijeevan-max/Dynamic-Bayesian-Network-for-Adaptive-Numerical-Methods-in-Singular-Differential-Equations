function [ret] = bvps2_paper_p1(request,z,za,zb,zc,t,p,lambda,lambda_p)
% Problem 1 of Hohenegger/Settanni/Weinmueller/Wolde, Appl. Numer. Math. 205 (2024):
%     y'' - y'/t - 2 = 0,   y(0) = 0,  y(1) = 0,    exact solution y = t^2 ln t.

switch request
    case 'n'
        ret = 1;
    case 'orders'
        ret = 2;
    case 'problem'
        ret = [ z(1,3) - z(1,2)/t - 2 ];
    case 'jacobian'
        %DON'T CHANGE THIS LINE:
        ret = zeros(length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders'))+1);
        ret(1,1,2) = -1/t;
        ret(1,1,3) = 1;
    case 'interval'
        ret = [0,1];
    case 'linear'
        ret = 1;
    case 'parameters'
        ret = 0;
    case 'c'
        ret = [];
    case 'BV'
        ret = [ za(1,1);
                zb(1,1) ];
    case 'dBV'
        %DON'T CHANGE THIS LINE:
        ret = zeros(max(length(feval(mfilename,'c')),2-length(feval(mfilename,'c'))),length(feval(mfilename,'problem',zeros(feval(mfilename,'n'),max(feval(mfilename,'orders'))+1),[],[],[],0,zeros(feval(mfilename,'parameters'),1))),feval(mfilename,'n'),max(feval(mfilename,'orders')));
        ret(1,1,1,1) = 1;
        ret(2,2,1,1) = 1;
    case 'dP'
        ret = [];
    case 'dP_BV'
        ret = [];
    case 'initProfile'
        ret.initialMesh   = linspace(0,1,11);
        ret.initialValues = ones(1,11);
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
