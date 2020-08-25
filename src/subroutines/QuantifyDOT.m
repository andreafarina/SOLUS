function Q = QuantifyDOT(REC,ref_true)
% =========================================================================
%%                            Quantify DOT
% =========================================================================
if any(strcmpi(REC.solver.variables,'mua'))
    for inl = 1:numel(REC.opt.muaB)
        fprintf(['<strong>------- Absorbtion: Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
        Q.mua(inl) = QuantifyX(REC.grid,REC.opt.muaB(inl),REC.opt.bmua(:,inl),...
            ref_true,REC.opt.hete1.c,REC.opt.Mua(:,:,:,inl));
    end
end
if any(strcmpi(REC.solver.variables,'mus'))
    for inl = 1:numel(REC.opt.muspB)
        fprintf(['<strong>------- Scattering: Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
        Q.mus(inl) = QuantifyX(REC.grid,REC.opt.muspB(inl),REC.opt.bmusp(:,inl),...
            ref_true,REC.opt.hete1.c,REC.opt.Musp(:,:,:,inl));
    end
end
if contains(lower(REC.solver.type),'spectral')
    if any(strcmpi(REC.solver.variables,'mua'))
        for inc = 1:REC.spe.nCromo
            fprintf(['<strong>------- Chromo ',REC.spe.cromo_label{inc},'-------</strong>\n'])
            Q.cromo(inc) = QuantifyX(REC.grid,REC.opt.concB(inc),REC.opt.bConc(:,inc),...
                ref_true,REC.opt.hete1.c,REC.opt.Conc(:,:,:,inc));
        end
    end
    if any(strcmpi(REC.solver.variables,'mus'))
        fprintf('<strong>------- A coeff -------</strong>\n')
         Q.a = QuantifyX(REC.grid,REC.opt.aB,REC.opt.bA,...
             ref_true,REC.opt.hete1.c,REC.opt.a);
        fprintf('<strong>------- B coeff -------</strong>\n')
         Q.b = QuantifyX(REC.grid,REC.opt.bB,REC.opt.bbB,...
             ref_true,REC.opt.hete1.c,REC.opt.b);
    end
end

%% FIGURE OF MERIT SOLUS
%load /scratch0/NOT_BACKED_UP/gdisciac/SOLUS_spectralfitUCL_bis/SOLUS/example/example_200202_series/dtprior/ex_bulk10_0.1_incl10_0.05/Test_Standard_REC.mat
for nlambda = 1:REC.radiometry.nL
    for n_mu = 1:2
    % extract data
    if n_mu ==  1
        mu = reshape( REC.opt.bmua(:,nlambda), REC.grid.dim);
        mu0 = REC.opt.mua0(nlambda);
        muin0 = REC.opt.hete1.val(nlambda);
        target_mu = REC.opt.Mua(:,:,:,nlambda);
        coeff = 'mua';
    else
        mu = reshape( REC.opt.bmusp(:,nlambda), REC.grid.dim);
        mu0 = REC.opt.musp0(nlambda);
        muin0 = REC.opt.hete1.val(nlambda + REC.radiometry.nL);
        target_mu = REC.opt.Musp(:,:,:,nlambda);
        coeff = 'musp';
    end
    % get region
    if mean(target_mu) >= mu0
        char_target = logical(target_mu > 0.5*(max(target_mu(:)) + min(target_mu(:))) );
    else
        char_target = logical(target_mu < 0.5*(max(target_mu(:)) + min(target_mu(:))));
    end
    idx_incl = identify_inclusion(mu(:));
    % get sensitivity
    [C,  CNR] = recon_sensitivity(mu(:), idx_incl);
    [C_vol,  CNR_vol] = recon_sensitivity(mu(:), idx_incl, 'volume', sum(char_target(:)), mu0);
    % get localisation
    [displ, brd, brdRMS] = recon_localisation(mu, char_target, REC.grid, REC.grid.dim, idx_incl);
    if isrow(displ)
        displ = displ';
    end
    if isrow(brd)
        brd = brd';
    end
    if isrow(brdRMS)
        brdRMS = brdRMS';
    end
    % get accuracy
    acc = recon_accuracy(mu, muin0, idx_incl);
    acc_vol = recon_accuracy(mu, muin0, idx_incl, 'volume', sum(char_target(:)), mu0);

    cmd = sprintf('Q.SOLUS_FigMerit.%s', coeff);
    cmd_end = sprintf('(:,%g)',nlambda);
    for i_field = {'C','C_vol', 'CNR','CNR_vol', 'displ',...
            'brd','brdRMS', 'acc', 'acc_vol'}
        eval([cmd,'.',  char(i_field), cmd_end,'=', char(i_field),';' ])                
    end
    end
end

end



% 
% %% subtract background
% dmua_rec3d = reshape(REC.opt.bmua-REC.opt.muaB,REC.grid.dim); 
% %dmua_rec3d = reshape(REC.opt.bmua,REC.grid.dim); 
% 
% %% Prepare Mask for gaussian fit
% radius_max = 20; %max radius of region of interest;
% if ref_true == 0
%     temp = dmua_rec3d;
%     %temp(:,:,1:10) = 0;
%     [mxm,mind] = max(temp(:));
% %mind = find(dmua_rec3d == mxm);
% [Xrec,Yrec,Zrec] = ind2sub(size(dmua_rec3d),mind);
% Xrec = x(Xrec);
% Yrec = y(Yrec);
% Zrec = z(Zrec);
% else
%     Xrec = REC.opt.hete1.c(1);
%     Yrec = REC.opt.hete1.c(2);
%     Zrec = REC.opt.hete1.c(3);
%     temp = dmua_rec3d;
%     
% end
% V = zeros(size(dmua_rec3d));
% nzm = sqrt((X-Xrec).^2+(Y-Yrec).^2+(Z-Zrec).^2)<radius_max;
% V(nzm) = 1;
% %V = ones(size(dmua_rec3d));
% dmua_rec3d = dmua_rec3d.*V;
% if ~exist('mxm','var')
%     mxm = max(dmua_rec3d(:));
% end
% %% Gaussian 3d estimation
% try
%     [COM_REC,COV_REC,TOT_VOL_REC] = Gaussian3D(x,y,z,dx,dy,dz,dmua_rec3d,V);
% 
% %% compare with reference
% if ref_true == 1
%     dmua_true3d = REC.opt.Mua-REC.opt.muaB; %% subtract background
%     dmua_true3d = dmua_true3d.*V;
%     [COM_TRUE,COV_TRUE,TOT_VOL_TRUE] = Gaussian3D(x,y,z,dx,dy,dz,dmua_true3d,V);
% end
% catch ME
%     if contains(ME.message,'lsqcurvefit')
%         disp('---- Error in fitting procedure ----')
%         Q = [];
%         return
%     else
%         throw(ME);
%     end
% end
% 
% %% Print report 
% 
% for iinc = 1:1
% disp(['------- REPORT for inclusion ',num2str(iinc),' ----------']);
% if ref_true == 1
%     disp(['REF_Tot_Volume: ',num2str(TOT_VOL_TRUE(iinc)),' mm^2']);
% end
% disp(['REC_Tot_Volume: ',num2str(TOT_VOL_REC(iinc)),' mm^2']);
% disp('--------------------------------------------');
% if ref_true == 1
%     disp(['REF_COM(x,y,z): ',num2str(COM_TRUE(iinc,:))]);
% end
% disp(['REC_COM(x,y,z): ',num2str(COM_REC(iinc,:))]);
% if ref_true == 1
%     disp(['X_COM_error: ',num2str(COM_REC(iinc,1)-COM_TRUE(iinc,1))]);
%     disp(['Y_COM_error: ',num2str(COM_REC(iinc,2)-COM_TRUE(iinc,2))]);
%     disp(['Z_COM_error: ',num2str(COM_REC(iinc,3)-COM_TRUE(iinc,3))]);
%     disp(['total COM_error: ',num2str(norm(COM_TRUE(iinc,:)-COM_REC(iinc,:)))]);
% end
% disp('--------------------------------------------');
% if ref_true == 1
%     disp('REF_COVARIANCE: ');
% %disp(num2str(COV_TRUE(:,:,iinc),3));
%     COV_TRUE(:,:,iinc)
% end
% disp('REC_COVARIANCE: ');
% %disp(num2str(COV_REC(:,:,iinc)));
% COV_REC(:,:,iinc)
% %disp('------------------------------------------');
% end
% 
% %% Sigma and CNR
% sigma=std(dmua_rec3d(:));
% disp(['Standard deviation ',num2str(sigma)]);
% 
% av=max(dmua_rec3d(:));
% disp(['Max ',num2str(av)]);
% disp(['CNR ',num2str(av./sigma)]);
% 
% av=max(V(:).*dmua_rec3d(:));
% disp(['Max in region ',num2str(av)]);
% disp(['CNR in region ',num2str(av./sigma)]);
% 
% 
% %% Quantitation
% %factor = REC.opt.hete1.sigma^3/(REC.grid.dx*REC.grid.dy*REC.grid.dz)*(4/3*pi);
% if ref_true == 1
%  DmuaVolTrue = sum(dmua_true3d(:)*dV);%/factor;
% end
%  DmuaVolRec = sum(dmua_rec3d(:)*dV);%/factor;
%  %disp(['DmuaVolTrue ',num2str(DmuaVolTrue)]);
%  %disp(['DmuaVolRec ',num2str(DmuaVolRec)]);
%  %disp(['error ',num2str(DmuaVolRec/DmuaVolTrue-1)]);
% %F(relMUAP).Value(iL)=(DmuaVolRec-DmuaVolTrue)/DmuaVolTrue;
% %% prepare mask for CNR (calculated as Arridge, Ducros..)
% W = zeros(size(dmua_rec3d));
% W(temp>mxm/2) = 1;
% Wback = 1 - W; 
% dmuaROI = W(:).*dmua_rec3d(:);
% dmuaROI(W==0) = [];
% muROI = mean(dmuaROI);
% dmuaBack = Wback(:).*dmua_rec3d(:);
% dmuaBack(Wback==0) = [];
% muBack = mean(dmuaBack);
% 
% CNR2 =(muROI - muBack)./...
%     sqrt(sum(W(:))./numel(W).*(std(dmuaROI(:)).^2) + ...
%      sum(Wback(:))./numel(W).*(std(dmuaBack(:)).^2));
% 
% Q.COM.rec = COM_REC;
% Q.COV.rec = COV_REC;
% Q.max.rec = mxm;
% Q.volumeG.rec = TOT_VOL_REC;
% Q.volume.rec = DmuaVolRec;
% Q.cnr = max(dmua_rec3d(:))./std(dmua_rec3d(:));
% Q.cnr2 = CNR2;
% disp(['CNR2 = ',num2str(CNR2)]);
% if ref_true == 1
%     Q.COM.true = COM_TRUE;
%     Q.COM.error = COM_REC - COM_TRUE;
%     Q.COV.true = COV_TRUE;
%     
%     Q.max.true = max(dmua_true3d(:));
%     Q.max.error = mxm -Q.max.true;
%     Q.max.rel_error = Q.max.error./Q.max.true;
%     
%     disp(['Relative max error = ',num2str(Q.max.rel_error)]);
%     
%     Q.volumeG.true = TOT_VOL_TRUE;
%     Q.volume.true = DmuaVolTrue;
%     Q.volumeG.rel_error = (TOT_VOL_REC - TOT_VOL_TRUE)./TOT_VOL_TRUE;
%     disp(['Relative volume (gaussian fit) error = ',num2str(Q.volumeG.rel_error)]);
%     Q.volume.rel_error = (DmuaVolRec-DmuaVolTrue)/DmuaVolTrue;
%     disp(['Relative volume error = ',num2str(Q.volume.rel_error)]);
%     Q.total.rel_error = norm(dmua_true3d(:)-dmua_rec3d(:))./norm(REC.opt.bmua+dmua_true3d(:));
%     disp(['Integral relative error = ',num2str(Q.total.rel_error)]);
%     
% end