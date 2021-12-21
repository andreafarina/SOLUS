function [fh,data] = displayImages(h,REC)
% [fh,data] = displayImages(h,REC)
% display images as the reconstruction core would, depending on whether a spectral solver was used
% h is the input figure handle
% REC is the REC structure after reconstruction
% fh returns the figure handle(s)
% data is a string describing what was plotted
	

    res_fact = [1,1,1];
    if ~contains(REC.solver.type,'fit')
        if ~contains(REC.solver.type,'spectral')
            Nr = ceil(sqrt(REC.radiometry.nL));Nc = round(sqrt(REC.radiometry.nL));
            nfig = 500;
            h
            %fh(1) = figure(nfig);
            for inl = 1:REC.radiometry.nL
                SubPlotMap(reshape(REC.opt.bmua(:,inl),REC.grid.dim),...
                    [num2str(REC.radiometry.lambda(inl)) ' nm'],[],Nr,Nc,inl,res_fact);
            end
            %fh(1).NumberTitle = 'off';fh.Name = 'Reconstructed Absorption';
            %    ------------------------ Reference musp ----------------------------------
            nfig = 600;
            %fh(2) = figure(nfig);
            fh = h;
            h(2)=fh;
            for inl = 1:REC.radiometry.nL
                SubPlotMap(reshape(REC.opt.bmusp(:,inl),REC.grid.dim),...
                    [num2str(REC.radiometry.lambda(inl)) ' nm'],[],Nr,Nc,inl,res_fact);
            end
            %fh(2).NumberTitle = 'off';fh.Name = 'Reconstructed Scattering';

            disp('recon: finished');
            data={'mua','musp'};
        elseif contains(REC.solver.type,'spectral')
            %% plot conc
            Nr = 3; Nc = 3;
            h;
            nfig=[];
            fh=[]
            for ic = 1:REC.spe.nCromo
                SubPlotMap(reshape(REC.opt.bConc(:,ic),REC.grid.dim),...
                    [REC.spe.cromo_label{ic} ' Map'],[],Nr,Nc,ic,res_fact);
            end
            % Hbtot and SO2
            REC.opt.HbTot = REC.opt.bConc(:,strcmpi(REC.spe.cromo_label,'hb'))+...
                REC.opt.bConc(:,strcmpi(REC.spe.cromo_label,'hbo2'));
            REC.opt.So2 = REC.opt.bConc(:,strcmpi(REC.spe.cromo_label,'hbo2'))./REC.opt.HbTot;
            SubPlotMap(reshape(REC.opt.HbTot,REC.grid.dim),'HbTot Map',nfig,Nr,Nc,ic+1,res_fact);
            SubPlotMap(reshape(REC.opt.So2,REC.grid.dim),'So2 Map',nfig,Nr,Nc,ic+2,res_fact);

            % a b scattering
            SubPlotMap(reshape(REC.opt.bA,REC.grid.dim),'a Map',nfig,Nr,Nc,ic+3,res_fact);
            SubPlotMap(reshape(REC.opt.bbB,REC.grid.dim),'b Map',nfig,Nr,Nc,ic+4,res_fact);

            %fh.NumberTitle = 'off';fh.Name = 'Reconstructed';
            
            data={'concs'};

        end
    end



end
