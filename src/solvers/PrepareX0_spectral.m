function x0 = PrepareX0_spectral(spe,n,type)

switch lower(type)
    case 'mua'
        if iscolumn(spe.opt.conc0), spe.opt.conc0 = spe.opt.conc0'; end
        x0 = repmat(spe.opt.conc0,n,1);
        x0 = reshape(x0,spe.nCromo*n,1);
    case 'd'
        x0 = repmat([spe.opt.a0 spe.opt.b0],n,1);
        x0 = reshape(x0,n*2,1);
    case 'muad'
        if iscolumn(spe.opt.conc0), spe.opt.conc0 = spe.opt.conc0'; end
        x0 = repmat(spe.opt.conc0,n,1);
        x0 = [x0 repmat([spe.opt.a0 spe.opt.b0],n,1)];
        x0 = reshape(x0,spe.nCromo*n+n*2,1);
end

