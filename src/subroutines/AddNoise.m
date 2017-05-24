function [yn,sd,yu] = AddNoise(x,type,param)
% function [y,sd] = AddNoise(x,type,param)
% INPUT
% x: data
% type: 'Gaussian', 'Poisson'
% param: if type 'Gaussian' param is the relative standard deviation
%        if type 'Poisson' param is the average number.
%         - for multiple data the maximum intenisty measurement will be
%           scaled to the maximum average number.
%         - if param is not spedificed it applies poissrnd assuming
%           integer data
% OUTPUT
% yn: noisy data
% sd: standard deviation
% yu: un-noisy data (for Poisson the normalized data before noise)
% Copyright Andrea Farina


switch lower(type)
    case 'gaussian'
        yn = x.*(1 + param.*randn(size(x)));
        if nargout > 1
            sd = x.*param;
            if nargout > 2
                yu = x;
            end
        end
    case 'poisson'
        if varargin > 2
            m = max(sum(x));
            k = param./m;
            yuu = k.*x;
            yn = poissrnd(round(full(yuu)));
            if nargout > 2
                yu = yuu;
            end
                
        else
            yn = poissrnd(round(full(x)));
            if nargout > 2
                yu = x;
            end
        end
        if nargout > 1
            sd = sqrt(y);
        end
        
end
