function [bmua,bmus,bconc,bAB] = XtoMuaMus_spectral(x,mua0,mus0,type,spe,conc0,a0,b0)

switch lower(type)
    case 'mua'
        nV = numel(x)/spe.nCromo;
        x=reshape(x,nV,spe.nCromo);
        bconc = x;
        bmua = (spe.ext_coeff0*bconc')';
        bmus = ones(size(bmua)) .* mus0;
        bAB = ones(nV,2).*[a0 b0];
    case 'd'
        nV = numel(x)/2;
        x=reshape(x,nV,2);
        bAB = x;
        bconc = ones(nV,spe.nCromo).*(conc0');
        bmus = zeros(nV,numel(spe.lambda));
        for in = 1:nV
           bmus(in,:) = bAB(in,1).*(spe.lambda./spe.lambda0).^(-bAB(in,2));
        end
        bmua = ones(size(bmus)) .* mua0;
    case 'muad'
        nV = numel(x)/(spe.nCromo+2);
        x=reshape(x,nV,spe.nCromo+2);
        bconc = x(:,1:spe.nCromo);
        bAB = x(:,spe.nCromo+1:end);
        bmua = (spe.ext_coeff0*bconc')';
        
        bmus = zeros(nV,numel(spe.lambda));
        for in = 1:nV
           bmus(in,:) = bAB(in,1).*(spe.lambda./spe.lambda0).^(-bAB(in,2));
        end
end
