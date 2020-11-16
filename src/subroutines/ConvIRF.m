function y = ConvIRF(x,irf)
% convolve temporal data x(ntime, nmeas) with irf and save the valid part
% if numel(irf) == 0 then returns x

nstep = size(x,1);
% Convolution with IRF
if numel(irf)>1
    z = convn(x,irf);
    nmax = max(nstep,numel(irf));
    y = z(1:nmax,:);
    clear nmax z
else
    y = x;
end