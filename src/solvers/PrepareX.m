function [x0,x] = PrepareX(mu,n,type,xtransf)
% Prepare X and X0 depending of the number of unknwons and transformation
% required
% A. Farina - CNR-IFN, November 2020
switch lower(type)
    case 'mua'
        x0 = mu(1) * ones(n,1);
    case 'd'
        x0 = mu(2) * ones(n,1);
    case 'muad'
        x0 = [mu(1)*ones(n,1);mu(2)*ones(n,1)];
end

switch lower(xtransf)
    case 'x'
        x = x0;
    case '(x/x0)'
        x = ones(size(x0));
    case 'logx'
        x = log(x0);
    case 'log(x/x0)'
        x = zeros(size(x0));
    otherwise
        disp('Unknown X transform');
end