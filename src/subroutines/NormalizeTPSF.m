function [y,A] = NormalizeTPSF(x)
% Area normalization of x(Ntime,Nmeas)
nmeas = size(x,2);
y = x * spdiags(1./sum(x,'omitnan')',0,nmeas,nmeas);
if nargout>1
    A=spdiags(1./sum(x,'omitnan')',0,nmeas,nmeas);
end
end