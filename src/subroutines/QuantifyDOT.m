function Q = QuantifyDOT(REC,ref_true)
% =========================================================================
%                            Quantify DOT
% =========================================================================
FullPrint = false;
Q = [];
if any(strcmpi(REC.solver.variables,'mua'))
    for inl = 1:numel(REC.opt.muaB)
        if FullPrint == true
            fprintf(['<strong>------- Absorbtion: Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
        end
        Q.mua(inl) = QuantifyX(REC.grid,REC.opt.mua0(inl),REC.opt.bmua(:,inl),...
            ref_true,REC.opt.hete1.c,REC.opt.Mua(:,:,:,inl),FullPrint);
        
    end
end
if any(strcmpi(REC.solver.variables,'mus'))
    for inl = 1:numel(REC.opt.muspB)
        if FullPrint == true
            fprintf(['<strong>------- Scattering: Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
        end
        Q.mus(inl) = QuantifyX(REC.grid,REC.opt.muspB(inl),REC.opt.bmusp(:,inl),...
            ref_true,REC.opt.hete1.c,REC.opt.Musp(:,:,:,inl),FullPrint);
    end
end
if contains(lower(REC.solver.type),'spectral')
    if any(strcmpi(REC.solver.variables,'mua'))
        for inc = 1:REC.spe.nCromo
            if FullPrint == true
                fprintf(['<strong>------- Chromo ',REC.spe.cromo_label{inc},'-------</strong>\n'])
            end
            Q.cromo(inc) = QuantifyX(REC.grid,REC.opt.concB(inc),REC.opt.bConc(:,inc),...
                ref_true,REC.opt.hete1.c,REC.opt.Conc(:,:,:,inc),FullPrint);
        end
    end
    if any(strcmpi(REC.solver.variables,'mus'))
        if FullPrint == true
            fprintf('<strong>------- A coeff -------</strong>\n')
        end
        Q.a = QuantifyX(REC.grid,REC.opt.aB,REC.opt.bA,...
            ref_true,REC.opt.hete1.c,REC.opt.A,FullPrint);
        if FullPrint == true
            fprintf('<strong>------- B coeff -------</strong>\n')
        end
        Q.b = QuantifyX(REC.grid,REC.opt.bB,REC.opt.bbB,...
            ref_true,REC.opt.hete1.c,REC.opt.B,FullPrint);
    end
end

% Add fields of the Q structure you want to display.
% Fields which have an array of values are NOT supported!
DispQuantifyTable(Q,REC,{'cnr' 'max.rec'});

%% FIGURE OF MERIT SOLUS
n_mu = numel(REC.opt.hete1.type);
for nlambda = 1:REC.radiometry.nL
    for n_mu = 1:n_mu
        if numel(REC.opt.hete1.type)==1
            if strcmpi(REC.opt.hete1.type,'Mua')
                n_mu = 1;
            elseif strcmpi(REC.opt.hete1.type,'Musp')
                n_mu = 2;
            end
        end
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
            if numel(REC.opt.hete1.type) == 2
                muin0 = REC.opt.hete1.val(nlambda + REC.radiometry.nL);
            else
                muin0 = REC.opt.hete1.val(nlambda);
            end
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


if contains(lower(REC.solver.type),'spectral')
    if numel(REC.opt.hete1.type)==1
        if strcmpi(REC.opt.hete1.type,'Mua')
            n_mu = 1;
        elseif strcmpi(REC.opt.hete1.type,'Musp')
            n_mu = 2;
        end
    end
    
    if n_mu == 1 || (numel(REC.opt.hete1.type) == 2)
        for ic = 1:REC.spe.nCromo
            mu = reshape(REC.opt.bConc(:,ic), REC.grid.dim);
            mu0 = REC.opt.conc0(ic);
            muin0 = REC.opt.hete1.conc(ic);
            target_mu = REC.opt.Conc(:,:,:,ic);
            if ic == 1
                coeff = 'Hb';
            elseif ic == 2
                coeff = 'HbO2';
            elseif ic == 3
                coeff = 'Lipid';
            elseif ic == 4
                coeff = 'H2O';
            elseif ic == 5
                coeff = 'Collagen';
            end
            %     get region
            if mean(target_mu(:)) >= mu0
                char_target = logical(target_mu > 0.5*(max(target_mu(:)) + min(target_mu(:))) );
            else
                char_target = logical(target_mu < 0.5*(max(target_mu(:)) + min(target_mu(:))));
            end
            idx_incl = identify_inclusion(mu(:));
            %     get sensitivity
            [C,  CNR] = recon_sensitivity(mu(:), idx_incl);
            [C_vol,  CNR_vol] = recon_sensitivity(mu(:), idx_incl, 'volume', sum(char_target(:)), mu0);
            %     get localisation
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
            %     get accuracy
            acc = recon_accuracy(mu, muin0, idx_incl);
            acc_vol = recon_accuracy(mu, muin0, idx_incl, 'volume', sum(char_target(:)), mu0);
            
            cmd = sprintf('Q.SOLUS_FigMerit.%s', coeff);
            cmd_end = sprintf('(:,%g)',ic);
            for i_field = {'C','C_vol', 'CNR','CNR_vol', 'displ',...
                    'brd','brdRMS', 'acc', 'acc_vol'}
                eval([cmd,'.',  char(i_field), cmd_end,'=', char(i_field),';' ])
            end
            % % %     create mask manually selected by user
            %     mask = manual_mask(REC.spe,ic,REC.grid,REC.opt.bConc(:,ic),coeff);
            %     mask = reshape(mask,size(REC.opt.bConc(:,ic)));
            %     manual = 1;
            %     idx_incl_manual = find(mask == 1);
            %     [C_manual,  CNR_manual, rec(ic), back(ic)] = recon_sensitivity(mu(:), idx_incl_manual);
            %     [C_vol_manual,  CNR_vol_manual] = recon_sensitivity(mu(:), idx_incl_manual, 'volume', sum(char_target(:)), mu0);
            %     [displ_manual, brd_manual, brdRMS_manual] = recon_localisation(mu, char_target, REC.grid, REC.grid.dim, idx_incl_manual, manual);
            %     if isrow(displ_manual)
            %         displ_manual = displ_manual';
            %     end
            %     if isrow(brd_manual)
            %         brd_manual = brd_manual';
            %     end
            %     if isrow(brdRMS_manual)
            %         brdRMS_manual = brdRMS_manual';
            %     end
            %     acc_manual = recon_accuracy(mu, muin0, idx_incl_manual);
            %     acc_vol_manual = recon_accuracy(mu, muin0, idx_incl_manual, 'volume', sum(char_target(:)), mu0);
            %     close all
            %     clear mask; clear idx_incl_manual;
            %     cmd_manual = sprintf('Q.SOLUS_FigMerit_manual.%s', coeff);
            %     cmd_end_manual = sprintf('(:,%g)',ic);
            %     for i_field = {'C_manual','C_vol_manual', 'CNR_manual','CNR_vol_manual', 'displ_manual',...
            %             'brd_manual','brdRMS_manual', 'acc_manual', 'acc_vol_manual'}
            %         eval([cmd_manual,'.',  char(i_field), cmd_end_manual,'=', char(i_field),';' ])
            %     end
            
        end
    end
    
    if n_mu == 2 || (numel(REC.opt.hete1.type) == 2)
        for ic = 1:2
            if ic == 1
                coeff = 'a';
                mu = reshape(REC.opt.bA(:), REC.grid.dim);
                mu0 = REC.opt.a0;
                muin0 = REC.opt.hete1.a;
                target_mu = REC.opt.A(:,:,:);
            else
                coeff = 'b';
                mu = reshape(REC.opt.bbB(:), REC.grid.dim);
                mu0 = REC.opt.b0;
                muin0 = REC.opt.hete1.b;
                target_mu = REC.opt.B(:,:,:);
            end
            % get region
            if mean(target_mu(:)) >= mu0
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
            cmd_end = sprintf('(:,%g)',ic);
            for i_field = {'C','C_vol', 'CNR','CNR_vol', 'displ',...
                    'brd','brdRMS', 'acc', 'acc_vol'}
                eval([cmd,'.',  char(i_field), cmd_end,'=', char(i_field),';' ])
            end
            %create mask manually selected by user
            %     if coeff == 'a'
            %          mask = manual_mask(REC.spe,ic,REC.grid,REC.opt.bA(:,1),coeff);
            %          mask = reshape(mask,size(REC.opt.bA(:,1)));
            %     else
            %         mask = manual_mask(REC.spe,ic,REC.grid,REC.opt.bbB(:,1),coeff);
            %         mask = reshape(mask,size(REC.opt.bbB(:,1)));
            %     end
            %     manual = 1;
            %     idx_incl_manual = find(mask == 1);
            %     [C_manual,  CNR_manual, rec_ab(ic), back_ab(ic)] = recon_sensitivity(mu(:), idx_incl_manual);
            %     [C_vol_manual,  CNR_vol_manual] = recon_sensitivity(mu(:), idx_incl_manual, 'volume', sum(char_target(:)), mu0);
            %     [displ_manual, brd_manual, brdRMS_manual] = recon_localisation(mu, char_target, REC.grid, REC.grid.dim, idx_incl_manual, manual);
            %     if isrow(displ_manual)
            %         displ_manual = displ_manual';
            %     end
            %     if isrow(brd_manual)
            %         brd_manual = brd_manual';
            %     end
            %     if isrow(brdRMS_manual)
            %         brdRMS_manual = brdRMS_manual';
            %     end
            %     acc_manual = recon_accuracy(mu, muin0, idx_incl_manual);
            %     acc_vol_manual = recon_accuracy(mu, muin0, idx_incl_manual, 'volume', sum(char_target(:)), mu0);
            %     close all
            %     cmd_manual = sprintf('Q.SOLUS_FigMerit_manual.%s', coeff);
            %     cmd_end_manual = sprintf('(:,%g)',ic);
            %     for i_field = {'C_manual','C_vol_manual', 'CNR_manual','CNR_vol_manual', 'displ_manual',...
            %             'brd_manual','brdRMS_manual', 'acc_manual', 'acc_vol_manual'}
            %         eval([cmd_manual,'.',  char(i_field), cmd_end_manual,'=', char(i_field),';' ])
            %     end
        end
    end
    
end

end



