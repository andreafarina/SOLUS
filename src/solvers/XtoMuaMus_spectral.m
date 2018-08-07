function [bmua,bmus,bconc,bAB] = XtoMuaMus_spectral(x,mua0,mus0,type,Vars)

switch lower(type)
    case 'mua'
        nV = numel(x)/Vars.nCromo;
        x=reshape(x,nV,Vars.nCromo);
        bconc = x;
        bmua = (Vars.ext_coeff0*bconc')';
        bmus = ones(size(bmua)) .* mus0;
        bAB = ones(nV,2).*[Vars.DOT.opt.aB Vars.DOT.opt.bB];
    case 'd'
        nV = numel(x)/2;
        x=reshape(x,nV,2);
        bAB = x;
        bconc = ones(nV,Vars.nCromo).*(Vars.DOT.opt.concB');
        for in = 1:nV
           bmus(in,:) = bAB(in,1).*(Vars.lambda/Vars.lambda0).^(-bAB(in,2));
        end
        bmua = ones(size(bmus)) .* mua0;
    case 'muad'
        nV = numel(x)/(Vars.nCromo+2);
        x=reshape(x,nV,Vars.nCromo+2);
        bconc = x(:,1:Vars.nCromo);
        bAB = x(:,Vars.nCromo+1:end);
        bmua = (Vars.ext_coeff0*bconc')';
        for in = 1:nV
           bmus(in,:) = bAB(in,1).*(Vars.lambda/Vars.lambda0).^(-bAB(in,2));
        end
end
