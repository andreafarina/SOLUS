function [y,Area] = NormalizeTPSF(x)
% Area normalization of x(Ntime,Nmeas)
% optionally it returns the Area
nmeas = size(x,2);
y = x * spdiags(1./sum(x,'omitnan')',0,nmeas,nmeas);
if nargout > 1
    Area = sum(x,1);
end