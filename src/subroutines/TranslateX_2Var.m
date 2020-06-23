function [DOT,REC]= TranslateX_2Var(SPECTRA,Xd,Xp,Xr,...
    spectra_file,DOT,REC)
mua_pos = 1; musp_pos = 2;
a_pos = 1; b_pos = 2; conc_pos = 1; ab_pos = 2;
DOT.spe.lambda0 = DOT.radiometry.lambda0;
DOT.spe.lambda = DOT.radiometry.lambda;
DOT.spe.nLambda = numel(DOT.spe.lambda);
DOT.spe.spectra_file = spectra_file;
DOT.spe.SPECTRA = SPECTRA;
DOT.spe.nCromo = numel(DOT.spe.cromo_label);
DOT.spe.ext_coeffB = LoadSpectra(spectra_file,DOT.radiometry.lambda,ones(1,DOT.spe.nCromo));
DOT.spe.ext_coeffB = DOT.spe.ext_coeffB.*DOT.spe.active_cromo;
DOT.spe.ext_coeff0 = DOT.spe.ext_coeffB;
REC.spe.ext_coeff0 = DOT.spe.ext_coeffB;
switch SPECTRA
    case 0
        DOT.opt.muaB = Xd{mua_pos}; DOT.opt.muspB = Xd{musp_pos}; 
        REC.opt.mua0 = Xr{mua_pos}; REC.opt.musp0 = Xr{musp_pos};
        DOT.opt.hete1.val = [];
        for in = 1:numel(Xp)
            DOT.opt.hete1.val = [DOT.opt.hete1.val Xp{in}];
        end
    case 1
        
        ab = Xd{ab_pos};
        DOT.opt.aB = ab(a_pos); DOT.opt.bB = ab(b_pos); DOT.opt.concB = Xd{conc_pos}.*DOT.spe.active_cromo;
        DOT.spe.ext_coeffB = LoadSpectra(spectra_file,DOT.radiometry.lambda,DOT.opt.concB);
        DOT.spe.ext_coeffB = DOT.spe.ext_coeffB.*DOT.spe.active_cromo;
        if isrow(DOT.opt.concB), DOT.opt.concB = DOT.opt.concB'; end
        DOT.opt.muaB = DOT.spe.ext_coeffB*DOT.opt.concB;
        DOT.opt.muspB = DOT.opt.aB.*(DOT.radiometry.lambda/DOT.radiometry.lambda0).^(-DOT.opt.bB);
        if iscolumn(DOT.opt.muaB), DOT.opt.muaB = DOT.opt.muaB'; end
        if iscolumn(DOT.opt.muspB), DOT.opt.muspB =DOT.opt.muspB'; end
        DOT.spe.opt.aB = DOT.opt.aB;DOT.spe.opt.bB = DOT.opt.bB;
        DOT.spe.opt.concB = DOT.opt.concB;
        
        ab = Xr{ab_pos};
        REC.opt.a0 = ab(a_pos); REC.opt.b0 = ab(b_pos); REC.opt.conc0 = Xr{conc_pos}.*DOT.spe.active_cromo;
        REC.spe.ext_coeff0 = LoadSpectra(spectra_file,DOT.radiometry.lambda,REC.opt.conc0);
        REC.spe.ext_coeff0 = REC.spe.ext_coeff0.*DOT.spe.active_cromo;
        if isrow(REC.opt.conc0), REC.opt.conc0 = REC.opt.conc0'; end
        REC.opt.mua0 = REC.spe.ext_coeff0*REC.opt.conc0;
        REC.opt.musp0 = REC.opt.a0.*(DOT.radiometry.lambda/DOT.radiometry.lambda0).^(-REC.opt.b0);
        if iscolumn(REC.opt.mua0), REC.opt.mua0 = REC.opt.mua0'; end
        if iscolumn(REC.opt.musp0), REC.opt.musp0 =REC.opt.musp0'; end
        REC.spe.opt.a0 = REC.opt.a0;REC.spe.opt.b0 = REC.opt.b0;
        REC.spe.opt.conc0 = REC.opt.conc0;
        if numel(DOT.opt.hete1.type) == 2
            ab = Xp{ab_pos};
            DOT.opt.hete1.a = ab(a_pos); DOT.opt.hete1.b = ab(b_pos); DOT.opt.hete1.conc = Xp{conc_pos}.*DOT.spe.active_cromo;
            CheckDataConsistency;
            if isrow(DOT.opt.hete1.conc), DOT.opt.hete1.conc = DOT.opt.hete1.conc'; end
            absorption = DOT.spe.ext_coeffB*DOT.opt.hete1.conc;
            DOT.spe.opt.hete1.conc = DOT.opt.hete1.conc;
            if iscolumn(absorption), absorption = absorption'; end
            
            scattering = DOT.opt.hete1.a .*(DOT.radiometry.lambda/DOT.radiometry.lambda0).^(-DOT.opt.hete1.b);
            DOT.spe.opt.hete1.a = DOT.opt.hete1.a;DOT.spe.opt.hete1.b = DOT.opt.hete1.b;
            if iscolumn(scattering), scattering = scattering'; end
            DOT.opt.hete1.val = [absorption scattering];
        elseif strcmpi(DOT.opt.hete1.type,'mua')
            DOT.opt.hete1.conc = Xp{1}.*DOT.spe.active_cromo;
            if isrow(DOT.opt.hete1.conc), DOT.opt.hete1.conc = DOT.opt.hete1.conc'; end
            CheckDataConsistency;
            absorption = DOT.spe.ext_coeffB*DOT.opt.hete1.conc;
            if iscolumn(absorption), absorption = absorption'; end
            DOT.opt.hete1.val = absorption;
            DOT.spe.opt.hete1.conc = DOT.opt.hete1.conc;
        elseif strcmpi(DOT.opt.hete1.type,'musp')
            ab = Xp{1};
            DOT.opt.hete1.a = ab(a_pos); DOT.opt.hete1.b = ab(b_pos);
            scattering = DOT.opt.hete1.a .*(DOT.radiometry.lambda/DOT.radiometry.lambda0).^(-DOT.opt.hete1.b);
            if iscolumn(scattering), scattering = scattering'; end
            DOT.opt.hete1.val = scattering;
            DOT.spe.opt.hete1.a = DOT.opt.hete1.a;DOT.spe.opt.hete1.b = DOT.opt.hete1.b;
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

CheckDataConsistency

end