function y = ConvIRF(x,irf)
% convolve temporal data x(ntime, nmeas) with irf and save the valid part
% if numel(irf) == 0 then returns x

nstep = size(x,1);
% Convolution with IRF
if numel(irf)>1
    if size(irf,2) ==1
    z = convn(x,irf);
    nmax = max(nstep,numel(irf));
    y = z(1:nmax,:);
    clear nmax z
    elseif size(irf,2) == size(x,2)
        nmax = max(nstep,numel(irf(:,1)));
        y = zeros(nmax,size(irf,2));
        for i = 1:size(irf,2)
            z = convn(x(:,i),irf(:,i));      
            y(:,i) = z(1:nmax,:);           
        end
        clear nmax z
    end
else
    y = x;
end