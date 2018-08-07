function x0 = PrepareX0_spectral(Vars,n,type)

switch lower(type)
    case 'mua'
        if iscolumn(Vars.REC.opt.conc0), Vars.REC.opt.conc0 = Vars.REC.opt.conc0'; end
        x0 = repmat(Vars.REC.opt.conc0,n,1);
        x0 = reshape(x0,Vars.nCromo*n,1);
    case 'd'
        x0 = repmat([Vars.REC.opt.a0 Vars.REC.opt.b0],n,1);
        x0 = reshape(x0,n*2,1);
    case 'muad'
        if iscolumn(Vars.REC.opt.conc0), Vars.REC.opt.conc0 = Vars.REC.opt.conc0'; end
        x0 = repmat(Vars.REC.opt.conc0,n,1);
        x0 = [x0 repmat([Vars.REC.opt.a0 Vars.REC.opt.b0],n,1)];
        x0 = reshape(x0,Vars.nCromo*n+n*2,1);
end

