function Vars = TranslateX_2Var(SPECTRA,Xd,Xp,Xr,lambda,type,lambda0,...
    spectra_file,cromo_label,cromo_units,cromo_factor,ForceConstitSolution)
mua_pos = 1; musp_pos = 2;
a_pos = 1; b_pos = 2; conc_pos = 1; ab_pos = 2;
Vars.Xp = Xp; Vars.Xd = Xd; Vars.Xr = Xr;
Vars.lambda = lambda; Vars.type = type;
Vars.lambda0 = lambda0; Vars.spectra_file = spectra_file;
Vars.cromo_label = cromo_label; Vars.SPECTRA = SPECTRA;
Vars.cromo_factor = cromo_factor;
Vars.conc_pos = conc_pos; Vars.cromo_units = cromo_units;
Vars.nCromo = numel(Vars.cromo_label);
Vars.ext_coeffB = LoadSpectra(Vars.spectra_file,Vars.lambda,ones(1,Vars.nCromo));
Vars.ext_coeff0 = Vars.ext_coeffB;
Vars.ForceConstitSolution = ForceConstitSolution;
switch SPECTRA
    case 0
        Vars.DOT.opt.muaB = Xd{mua_pos}; Vars.DOT.opt.muspB = Xd{musp_pos};
        Vars.REC.opt.mua0 = Xr{mua_pos}; Vars.REC.opt.musp0 = Xr{musp_pos};
        Vars.lambda = lambda;
        Vars.nLambda = numel(lambda);
        Vars.DOT.opt.hete1.val = [];
        for in = 1:numel(Vars.Xp)
            Vars.DOT.opt.hete1.val = [Vars.DOT.opt.hete1.val Xp{in}];
        end
    case 1
        
        ab = Xd{ab_pos};
        Vars.DOT.opt.aB = ab(a_pos); Vars.DOT.opt.bB = ab(b_pos); Vars.DOT.opt.concB = Xd{conc_pos};
        Vars.ext_coeffB = LoadSpectra(Vars.spectra_file,Vars.lambda,Vars.DOT.opt.concB);
        if isrow(Vars.DOT.opt.concB), Vars.DOT.opt.concB = Vars.DOT.opt.concB'; end
        Vars.DOT.opt.muaB = Vars.ext_coeffB*Vars.DOT.opt.concB;
        Vars.DOT.opt.muspB = Vars.DOT.opt.aB.*(Vars.lambda/Vars.lambda0).^(-Vars.DOT.opt.bB);
        if iscolumn(Vars.DOT.opt.muaB), Vars.DOT.opt.muaB = Vars.DOT.opt.muaB'; end
        if iscolumn(Vars.DOT.opt.muspB), Vars.DOT.opt.muspB =Vars.DOT.opt.muspB'; end
        
        ab = Xr{ab_pos};
        Vars.REC.opt.a0 = ab(a_pos); Vars.REC.opt.b0 = ab(b_pos); Vars.REC.opt.conc0 = Xr{conc_pos};
        Vars.ext_coeff0 = LoadSpectra(Vars.spectra_file,Vars.lambda,Vars.REC.opt.conc0);
        if isrow(Vars.REC.opt.conc0), Vars.REC.opt.conc0 = Vars.REC.opt.conc0'; end
        Vars.REC.opt.mua0 = Vars.ext_coeff0*Vars.REC.opt.conc0;
        Vars.REC.opt.musp0 = Vars.REC.opt.a0.*(Vars.lambda/Vars.lambda0).^(-Vars.REC.opt.b0);
        if iscolumn(Vars.REC.opt.mua0), Vars.REC.opt.mua0 = Vars.REC.opt.mua0'; end
        if iscolumn(Vars.REC.opt.musp0), Vars.REC.opt.musp0 =Vars.REC.opt.musp0'; end
        
        if numel(type) == 2
            ab = Xp{ab_pos};
            Vars.DOT.opt.hete1.a = ab(a_pos); Vars.DOT.opt.hete1.b = ab(b_pos); Vars.DOT.opt.hete1.conc = Xp{conc_pos};
            Vars = CheckDataConsistency(Vars);
            if isrow(Vars.DOT.opt.hete1.conc), Vars.DOT.opt.hete1.conc = Vars.DOT.opt.hete1.conc'; end
            absorption = Vars.ext_coeffB*Vars.DOT.opt.hete1.conc;
            if iscolumn(absorption), absorption = absorption'; end
            
            scattering = Vars.DOT.opt.hete1.a .*(Vars.lambda/Vars.lambda0).^(-Vars.DOT.opt.hete1.b);
            if iscolumn(scattering), scattering = scattering'; end
            
            Vars.scattering = scattering;
            Vars.absorption = absorption;
            
            Vars.DOT.opt.hete1.val = [absorption scattering];
        elseif strcmpi(type,'mua')
            Vars.DOT.opt.hete1.conc = Xp{1};
            if isrow(Vars.DOT.opt.hete1.conc), Vars.DOT.opt.hete1.conc = Vars.DOT.opt.hete1.conc'; end
            Vars = CheckDataConsistency(Vars);
            absorption = Vars.ext_coeffB*Vars.DOT.opt.hete1.conc;
            if iscolumn(absorption), absorption = absorption'; end
            Vars.absorption = absorption;
            Vars.DOT.opt.hete1.val = absorption;
        elseif strcmpi(type,'musp')
            ab = Xp{1};
            Vars.DOT.opt.hete1.a = ab(a_pos); Vars.DOT.opt.hete1.b = ab(b_pos);
            scattering = Vars.DOT.opt.hete1.a .*(Vars.lambda/Vars.lambda0).^(-Vars.DOT.opt.hete1.b);
            if iscolumn(scattering), scattering = scattering'; end
            Vars.scattering = scattering;
            Vars.DOT.opt.hete1.val = scattering;
        end
%         DOT.opt.spectral_hete1 = DOT.opt.hete1;
%         if contains(lower(DOT.opt.hete1.type),'mua')
%             DOT.opt.spectral_hete1.type = 'conc';
%         end
%         if contains(lower(DOT.opt.hete1.type),'mus')
%             DOT.opt.spectral_hete1.type = 'ab';
%         end
%         if contains(lower(DOT.opt.hete1.type),'mua')&&contains(lower(DOT.opt.hete1.type),'mus')
%             DOT.opt.spectral_hete1.type = 'concab';
%         end
end

Vars = CheckDataConsistency(Vars);

end

function Vars = CheckDataConsistency(Vars)
iserr = false;
if Vars.SPECTRA
    if numel(Vars.cromo_label)~=numel(Vars.Xd{Vars.conc_pos})
        err_message = '----DOT: Number of chromophors must be equal to incput concentations----';
        iserr = true;
    elseif numel(Vars.cromo_label)~=numel(Vars.Xr{Vars.conc_pos})
        err_message = '----REC: Number of chromophors must be equal to incput concentations----';
        iserr = true;
    elseif isfield(Vars.DOT.opt.hete1,'conc')
        if numel(Vars.DOT.opt.hete1.conc)~=numel(Vars.Xr{Vars.conc_pos})
            err_message = '----PERT: Number of chromophors must be equal to incput concentations----';
            iserr = true;
        end
    elseif isfield(Vars,'absorption')
        if numel(Vars.absorption)~=numel(Vars.lambda)
            err_message = '----PERT: Number of chromophors must be equal to incput concentations----';
            iserr = true;
        end
    elseif isfield(Vars,'scattering')
        if numel(Vars.scattering)~=numel(Vars.lambda)
            err_message = '----PERT: Number of chromophors must be equal to incput concentations----';
            iserr = true;
        end
    end
end
if numel(Vars.DOT.opt.muaB)~=numel(Vars.DOT.opt.muspB)
    err_message = '----DOT: Number of Mua and Mus must be the same----';
    iserr = true;
elseif numel(Vars.DOT.opt.muaB)~=numel(Vars.lambda)
    err_message = '----DOT: Check number of Wavelenghts or muaB/muspB----';
    iserr = true;
elseif numel(Vars.REC.opt.mua0)~=numel(Vars.REC.opt.musp0)
    err_message = '----REC: Number of Mua and Mus must be the same----';
    iserr = true;
elseif numel(Vars.REC.opt.mua0)~=numel(Vars.lambda)
    err_message = '----REC: Check number of Wavelenghts or mua/musp----';
    iserr = true;
elseif numel(Vars.Xp)~=numel(Vars.type)
    err_message = '----DOT: Number of perturbations must be equal to DOT.opt.hete1.type----';
    iserr = true;
elseif Vars.SPECTRA == 0
    if numel(Vars.Xp) == 2
        if numel(Vars.Xp{1})~=numel(Vars.Xp{2})
            err_message = '----DOT: Number of wavelenghs in perturbations must be equal----';
            iserr = true;
        end
    elseif numel(Vars.Xp{1})~=numel(Vars.lambda)
        err_message = '----DOT: Number of wavelenghs in perturbations must be equal----';
        iserr = true;
    end
end
if iserr
    MEx = MException('TranslateX_2Var:Dimension_mismatch',err_message);
    throwAsCaller(MEx);
end
end
