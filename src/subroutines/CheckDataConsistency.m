iserr = false;
if SPECTRA
    if numel(DOT.spe.cromo_label)~=numel(Xd{conc_pos})
        err_message = '----DOT: Number of chromophors must be equal to incput concentations----';
        iserr = true;
    elseif numel(DOT.spe.cromo_label)~=numel(Xr{conc_pos})
        err_message = '----REC: Number of chromophors must be equal to incput concentations----';
        iserr = true;
    elseif isfield(DOT.opt.hete1,'conc')
        if numel(DOT.opt.hete1.conc)~=numel(Xr{conc_pos})
            err_message = '----PERT: Number of chromophors must be equal to incput concentations----';
            iserr = true;
        end
    elseif isfield(Vars,'absorption')
        if numel(absorption)~=DOT.radiometry.nL
            err_message = '----PERT: Number of chromophors must be equal to incput concentations----';
            iserr = true;
        end
    elseif isfield(Vars,'scattering')
        if numel(scattering)~=DOT.radiometry.nL
            err_message = '----PERT: Number of chromophors must be equal to incput concentations----';
            iserr = true;
        end
    end
end
if numel(DOT.opt.muaB)~=numel(DOT.opt.muspB) %rimetti questo
    err_message = '----DOT: Number of Mua and Mus must be the same----';
    iserr = true;
elseif numel(DOT.opt.muaB)~=DOT.radiometry.nL %rimetti questo
    err_message = '----DOT: Check number of Wavelenghts or muaB/muspB----';
    iserr = true;
elseif numel(REC.opt.mua0)~=numel(REC.opt.musp0) %rimetti questo
    err_message = '----REC: Number of Mua and Mus must be the same----';
    iserr = true;
elseif numel(REC.opt.mua0)~=DOT.radiometry.nL
    err_message = '----REC: Check number of Wavelenghts or mua/musp----';
    iserr = true;
elseif numel(Xp)~=numel(DOT.opt.hete1.type)
    err_message = '----DOT: Number of perturbations must be equal to DOT.opt.hete1.type----';
    iserr = true;
elseif SPECTRA == 0
    if numel(Xp) == 2
        if numel(Xp{1})~=numel(Xp{2})
            err_message = '----DOT: Number of wavelenghs in perturbations must be equal----';
            iserr = true;
        end
    elseif numel(Xp{1})~=DOT.radiometry.nL
        err_message = '----DOT: Number of wavelenghs in perturbations must be equal----';
        iserr = true;
    end
end
if iserr
    MEx = MException('TranslateX_2Var:Dimension_mismatch',err_message);
    throwAsCaller(MEx);
end