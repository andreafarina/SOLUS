%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOT_core     %
%                                                                       %
% v3:                                                                   %
% Simulation and reconstruction using Green's function and TOAST_FEM.                 %
% Experimental data are not yet implemented.                            %
% Half-space geometry with fwd and Jacobian for the fluence under PCBC. %
% CW and time-domain forward and reconstruction.                        %
% Absorption only reconstruction.                                       %
% Possibility to load an a-priori dictionary matrix (see row 526).      %
% L1 ISTA reconstruction is implemented.
% Laplacian Structured priors regularization is implemented
% Regularization parameter criteria: L-curve, gcv, manual.              %
% Solver: Tickhonov, Truncated SVD, pcg, gmres                          %
%                                                                       %
% A. Farina 15/06/2018                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%clearvars;
if ~any(strcmpi([fileparts(mfilename('fullpath')) filesep 'solvers'],regexp(path, pathsep, 'split')))
    ComingFrom = cd; cd(fileparts([mfilename('fullpath') '.m'])); cd('../');
    disp('Path are not set! DOT_install.m will be run');
    run('DOT_install');
    cd(ComingFrom);
end

close all;
global mesh
mesh = [];
% setPath;
% addpath subroutines solvers
% addpath(genpath('util'))
t_start = cputime;
spacer = ' ------- ';
% =========================================================================
%%                          INITIALIZATION
% =========================================================================
frontpage;
disp('Initializing DOT parameters');
%run(fullfile(folder,'Init_DOT.m'));
Init_DOT;
%disp('Setting paths and filenames');
%run(fullfile(folder,'SetIO_DOT.m'));
SetIO_DOT;
%==========================================================================
%%  SETTING SOURCES (QVEC), DETECTORS (MVEC) AND THEIR PERMUTATIONS (DMASK)
%==========================================================================
disp('Setting sources(Q) and detectors(M)');
%run(fullfile(folder,'SetQM_DOT.m'));
SetQM_DOT;
%%
if RECONSTRUCTION == 1
    disp('Setting reconstruction parameters');
    %run(fullfile(folder,'RecSettings_DOT.m'));
    RecSettings_DOT;
end
% =========================================================================
%%               OVERRIDE param using P values (multiSIM)
% =========================================================================
Override_MultiSim

% =========================================================================
%%                     Create folder for saving data
% =========================================================================
rdir = res_path;
%['..',filesep,'results',filesep,datestr(clock, 'yyyymmdd'),filesep];
if ~exist(rdir,'dir')
    mkdir(rdir)
end
disp(['Results path: ',rdir]);
%% Get True Variables from X
[DOT,REC] = TranslateX_2Var(SPECTRA,Xd,Xp,Xr,spectra_file,DOT,REC);

%==========================================================================
%%                          Experimental data
% Load the structure EXP containing all the data for analyzing experimental
% data. This structure is created with the routine 'CreateExperiment.m'.
%==========================================================================
if (EXPERIMENTAL == 1)
    load([exp_path,'EXP_',exp_file])
    % Select measurements by dmask
    if (EXP_DATA == 1)
        if (size(EXP.data.spc,2)>sum(DOT.dmask))
            EXP.data.spc = EXP.data.spc(:,repmat(DOT.dmask,1,DOT.radiometry.nL));
            EXP.data.ref = EXP.data.ref(:,repmat(DOT.dmask,1,DOT.radiometry.nL));
            if DOT.time.self_norm
                nm = sum(DOT.dmask(:));
                for inl = 1:DOT.radiometry.nL
                    meas_set = (1:nm)+(inl-1)*nm;
                    EXP.data.spc(:,meas_set) = EXP.data.spc(:,meas_set)*spdiags(1./sum(EXP.data.spc(:,meas_set),1)',0,nm,nm);
                    EXP.data.ref(:,meas_set) = EXP.data.ref(:,meas_set)*spdiags(1./sum(EXP.data.ref(:,meas_set),1)',0,nm,nm);
                end
            end
        end
    end
end
disp(['Experiment file: ',exp_path,'EXP_',exp_file]);

%% Background optical properties
DOT.opt.cs = 0.299/DOT.opt.nB;       % speed of light in medium
DOT.opt.kap = 1./(3*DOT.opt.muspB);
DOT.A = A_factor(DOT.opt.nB/DOT.opt.nE); % A factor for boundary conditions
%==========================================================================
%%                              SET FWD FEM MESH/GRID
%==========================================================================
if strcmpi(TYPE_FWD,'fem')
    [vtx,idx,eltp] = mkslab([DOT.grid.x1,DOT.grid.y1,DOT.grid.z1;...
        DOT.grid.x2,DOT.grid.y2,DOT.grid.z2],...
        round(([DOT.grid.x2,DOT.grid.y2,DOT.grid.z2]-...
        [DOT.grid.x1,DOT.grid.y1,DOT.grid.z1])./...
        [DOT.grid.dx,DOT.grid.dy,DOT.grid.dz]));
    mesh.hMesh=toastMesh(vtx,idx,eltp);
    clear vtx idx eltp
end
%==========================================================================
%%                                  SET GRID
%==========================================================================
DOT.grid = setGrid(DOT);

[DOT.opt.Mua, DOT.opt.Musp] = applyGrid(DOT.grid,DOT.opt.muaB,DOT.opt.muspB);

%==========================================================================
%%                      Set Heterogeneities
%==========================================================================
%--------------------------- INCLUSIONS ---------------------------------%
for i = 1:NUM_HETE
    h_string = ['hete',num2str(i)];
    [DOT,DOT.opt.(h_string)] = setHete(DOT,DOT.opt.(h_string));
end

% Map optical properties to mesh and QM
if strcmpi(TYPE_FWD,'fem')
    for inl = 1: DOT.radiometry.nL
        mesh.opt.mua(:,inl) = DOT.grid.hBasis.Map('B->M',DOT.opt.Mua(:,:,:,inl));
        mesh.opt.musp(:,inl) = DOT.grid.hBasis.Map('B->M',DOT.opt.Musp(:,:,:,inl));
        %Qds = 1; % width of Sources
        %Mds = 1; % width of Detectors
    end        
    mesh.hMesh.SetQM(DOT.Source.Pos,DOT.Detector.Pos);
    mesh.qvec = real(mesh.hMesh.Qvec('Neumann','Gaussian',DOT.Source.Area));
    %mesh.hMesh.SetQM(DOT.Source.Pos+[0,0,1./DOT.opt.muspB],DOT.Detector.Pos);
    %mesh.qvec = real(mesh.hMesh.Qvec('Isotropic','Gaussian',DOT.Source.Area));
    mesh.mvec = 1./(2*DOT.A)*real(mesh.hMesh.Mvec('Gaussian',DOT.Detector.Area, 0));    
end

%

%==========================================================================
%%                         Time domain parameters
%==========================================================================
if DOT.TD == 1
    % --------- resample IRF following the parameters in DOT.time -------------
    if ((EXPERIMENTAL == 1) && (EXP_IRF == 1))
        switch lower(EXP_DELTA)
            case 'peak'
                EXP.irf.data = zeros(size(EXP.irf.data));
                peak_pos = (0:DOT.radiometry.nL-1).*(size(EXP.irf.data,1)) + [EXP.irf.peak.pos];
                EXP.irf.data(peak_pos) = 1;
            case 'baric'
                EXP.irf.data = zeros(size(EXP.irf.data));
                baric_pos = (0:DOT.radiometry.nL-1).*(size(EXP.irf.data,1)) + [EXP.irf.peak.pos];
                EXP.irf.data(baric_pos) = 1;
        end
        for inl = 1:DOT.radiometry.nL
            DOT.time.irf.data(:,inl) = resampleSPC(EXP.irf.data(:,inl),EXP.time.axis,DOT.time.dt,'norm');
        end
        % ------ setting channel 0 based on the peak and width of the IRF ---------
        %resc_fact = round(DOT.time.dt./EXP.spc.factor);
        [DOT.time.irf.peak.value,...
            DOT.time.irf.peak.pos] = max(DOT.time.irf.data,[],1);
        DOT.time.irf.ch0 = round(DOT.time.irf.peak.pos - ...
            500./DOT.time.dt);
        
    else
        DOT.time.irf.data = 0;
    end
    
end


% -------------------------------------------------------------------------
% plot DMASK
figure(10), imagesc(DOT.dmask),xlabel('Sources'),ylabel('Meas'), axis image,title('Dmask');

% plot source-detectors and heterogeneities
figure(11)
subplot(1,2,1),PlotHeteQM(DOT,squeeze(DOT.opt.Mua(:,:,:,1)),DOT.opt.muaB(1)),title('Mua');
subplot(1,2,2),PlotHeteQM(DOT,squeeze(DOT.opt.Musp(:,:,:,1)),DOT.opt.muspB(1)),title('Musp');
drawnow;

%==========================================================================
% The structure DOT contains all geometrical parameters needed also for
% the inverse problem
%==========================================================================
%%                             FORWARD PROBLEM
%==========================================================================
%disp('Forwar model computation');
nmeas = sum(DOT.dmask(:));
if FORWARD == 1
    DataCW = ForwardCW_multi_wave(DOT.grid,DOT.Source.Pos, DOT.Detector.Pos, DOT.dmask, ...
        DOT.opt.muaB, DOT.opt.muspB, DOT.opt.Mua, ...
        DOT.opt.Musp, DOT.A, geom, 'born',DOT.radiometry);
    sdCW = zeros(nmeas,DOT.radiometry.nL);
    for inl = 1:DOT.radiometry.nL
        [DataCW(:,inl),sdCW(:,inl)] = AddNoise(DataCW(:,inl),'gaussian',DOT.sigma);
    end
    if REF == 1
        [MuaB,MuspB]=applyGrid(DOT.grid,DOT.opt.muaB,DOT.opt.muspB);
        RefCW  = ForwardCW_multi_wave(DOT.grid,DOT.Source.Pos, DOT.Detector.Pos, DOT.dmask,...
            DOT.opt.muaB, DOT.opt.muspB, MuaB,...
            MuspB, DOT.A, geom, 'homo',DOT.radiometry);
        RefsdCW = zeros(nmeas,DOT.radiometry.nL);
        for inl = 1:DOT.radiometry.nL
            [RefCW(:,inl),RefsdCW(:,inl)] = AddNoise(RefCW(:,inl),'gaussian',DOT.sigma);
        end
        clear MuaB MuspB
    end
    %==========================================================================
    %%                             TOAST2DOT
    %==========================================================================
    if TOAST2DOT == 1
        toastfilename = [getenv('DOTDIR'),filesep,'TOAST_Sim',filesep,filename];
        [time_domain, meshdata, DataTD, RefTD, DataCW, RefCW, nominalcoefficients ,...
            refindex, Muatoast, Mustoast, prior] = toast2dot(toastfilename);
        save([rdir,filename,'_', 'FwdTeo'],'RefTD','DataTD');
        %     DOT.time.dt = time_domain(1);
        %     DOT.time.nstep = time_domain(2);
        %     DOT.time.time = (1:DOT.time.nstep) * DOT.time.dt;
        %     DOT.time.irf.data = resampleSPC(DOT.time.irf.data,DOT.time.time,...
        %         DOT.time.dt,'norm');
    end
    
    %==========================================================================
    %%                             TD Forward
    %==========================================================================
    if DOT.TD == 1
        DOT.time.time = (1:DOT.time.nstep) * DOT.time.dt;
        if LOAD_FWD_TEO == 0
            DataTD = ForwardTD_multi_wave(DOT.grid,DOT.Source.Pos, DOT.Detector.Pos, DOT.dmask,...
                DOT.opt.muaB, DOT.opt.muspB,DOT.opt.nB, DOT.opt.Mua,...
                DOT.opt.Musp, DOT.A, DOT.time.dt,...
                length(DOT.time.time), DOT.time.self_norm, geom, TYPE_FWD,DOT.radiometry);
            save([rdir,filename,'_', 'FwdTeo'],'DataTD');
            if REF == 1
                [MuaB,MuspB]=applyGrid(DOT.grid,DOT.opt.muaB,DOT.opt.muspB);
                RefTD = ForwardTD_multi_wave(DOT.grid,DOT.Source.Pos, DOT.Detector.Pos, DOT.dmask,...
                    DOT.opt.muaB, DOT.opt.muspB,DOT.opt.nB, ...
                    MuaB,MuspB, ...
                    DOT.A, DOT.time.dt,length(DOT.time.time), DOT.time.self_norm,...
                    geom, TYPE_FWD,DOT.radiometry);
                save([rdir,filename,'_', 'FwdTeo'],'RefTD','-append');
                clear MuaB MuspB
            end
        else
            load([rdir,filename,'_', 'FwdTeo']);
            CheckZerosPhi(DataTD,DOT.radiometry,'DataTD');
            CheckZerosPhi(RefTD,DOT.radiometry,'RefTD');
        end
        
        % Convolution with IRF
        % AF: to optimize memory assign directly DataTD
        if ((EXPERIMENTAL == 1) && (EXP_IRF == 1))
            nsub=numSubplots(DOT.radiometry.nL);
            figure(12);
            for inl = 1:DOT.radiometry.nL
                meas_set = (1:nmeas)+(inl-1)*nmeas;
                subplot(nsub(1),nsub(2),inl)
                z = convn(full(DataTD(:,meas_set)),DOT.time.irf.data(:,inl));
                semilogy(z(:,1),'b'),hold,semilogy(DataTD(:,meas_set(1)),'k'),
                semilogy(DOT.time.irf.data(:,inl),'r'),%ylim([max(DataTD(:))/10000 max(DataTD(:))])
                title(['Wavelength: ' num2str(DOT.radiometry.lambda(inl))]);
                nmax = max(DOT.time.nstep,numel(DOT.time.irf.data(:,inl)));
                if inl == 1
                    DataTD(1:nmax,meas_set) = z(1:nmax,:);
                    DataTD(nmax+1:end,:) = [];
                else
                    DataTD(:,meas_set) = z(1:nmax,:);
                end
                
                if REF == 1
                    z = convn(full(RefTD(:,meas_set)),DOT.time.irf.data(:,inl));
                    if inl == 1
                        RefTD(1:nmax,meas_set) = z(1:nmax,:);
                        RefTD(nmax+1:end,:) = [];
                    else
                        RefTD(:,meas_set) = z(1:nmax,:);
                    end
                end
                clear nmax
            end
            drawnow
        end
        clear z inl meas_set
        
        % ---- Radiometry -------
        if RADIOMETRY == 1
            RealFactor = Radiometry(DOT.radiometry);
        else
            RealFactor = ones(1,DOT.radiometry.nL);
        end
        %-------------------- Add noise to TD data ------------------------
        sdTD = zeros(size(DataTD));
        for inl = 1:DOT.radiometry.nL
            meas_set = (1:nmeas)+(inl-1)*nmeas;
            DataTD_single_wave = DataTD(:,meas_set);
            RefTD_single_wave = RefTD(:,meas_set);
            sdTD_single_wave = ones(size(DataTD_single_wave));
            if ~strcmpi(DOT.time.noise,'none')
                [DataTD_single_wave,~] = AddNoise(DataTD_single_wave,'gaussian',DOT.time.sigma);
                if REF == 1
                    RefTD_single_wave = AddNoise(RefTD_single_wave,'gaussian',DOT.time.sigma);
                end
            end
            
            if strcmpi(DOT.time.noise,'poisson')
                
                if (REF == 1) %&& (DOT.time.self_norm == 0))
                    [factor,RefTD_single_wave,sdTD_single_wave] = CutCounts(1:size(RefTD_single_wave,1),DOT.time.dt,...
                        RealFactor(inl)*RefTD_single_wave,DOT.time.TotCounts(inl)*DOT.radiometry.acqtime,RADIOMETRY,CUT_COUNTS,...
                        NumDelays);
                    switch CUT_COUNTS
                        case 0  %
                            DataTD_single_wave = poissrnd(round(...
                                DataTD_single_wave * factor * RealFactor(inl)));
                        case 1
                            %idz = find(factor>0, 1, 'last' );
                            DataTD_single_wave = poissrnd(round(DataTD_single_wave * RealFactor(inl) .* factor));
                            %DataTD = bsxfun(@times,DataTD,factor(idz)./factor');
                    end
                    
                    
                else
                    [factor,DataTD_single_wave,sdTD_single_wave] = CutCounts(1:size(RefTD_single_wave,1),DOT.time.dt,...
                        RealFactor(inl)*DataTD_single_wave,DOT.time.TotCounts(inl),RADIOMETRY,CUT_COUNTS,...
                        NumDelays);
                end
            end
            % self normalized data
            if DOT.time.self_norm == true
                Area = sum(DataTD_single_wave,'omitnan');
                DataTD_single_wave = DataTD_single_wave * spdiags(1./Area',0,nmeas,nmeas);
                sdTD_single_wave = sqrt(DataTD_single_wave * spdiags(1./Area',0,nmeas,nmeas));  % Standard deviation
                if REF == 1
                    Area = sum(RefTD_single_wave,'omitnan');
                    RefTD_single_wave = RefTD_single_wave * spdiags(1./Area',0,nmeas,nmeas);
                    sdTD_single_wave = sqrt(DataTD_single_wave * spdiags(1./Area',0,nmeas,nmeas));  % Standard deviation
                end
            end
            DataTD(:,meas_set) = DataTD_single_wave;
            RefTD(:,meas_set) = RefTD_single_wave;
            sdTD(:,meas_set) = sdTD_single_wave;
        end
    end
    
    
end
clear inl sh meas_set
%==========================================================================
t_end_fwd = cputime - t_start;
disp(['CPU time FWD: ' num2str(t_end_fwd)])
% =========================================================================
%%                            Save forward data
% =========================================================================
if SAVE_FWD == 1
    %rdir = ['..',filesep,'results',filesep,session,filesep];
    if ~exist(rdir,'dir')
        mkdir(rdir)
    end
    disp(['Results will be stored in: ', rdir,filename,'_','Data.m'])
    %-------------------------- Save simulated data ---------------------------
    if FORWARD == 1
        save([rdir,filename,'_', 'Data'],'DataCW','sdCW');
        if REF == 1
            save([rdir,filename,'_', 'Data'],'RefCW','-append');
        end
        if DOT.TD == 1
            save([rdir,filename,'_', 'Data'],'DataTD','sdTD','-append');
            if REF == 1
                save([rdir,filename,'_', 'Data'],'RefTD','-append');
            end
        end
        
    end
end
%==========================================================================
%%                              END FORWARD
%==========================================================================
%==========================================================================
%%                              EXPERIMENTAL
% if EXP_DATA==1 in this section the experimental SPC data are resampled
% accordingly to the time-step defined in DOT.time substructure. If DOT.TD
% == 0 (CW) the resmapling is skipped and the integral of the time-resolved
% curves is calculated.
%==========================================================================
if ((EXPERIMENTAL == 1) && (EXP_DATA == 1))
    if DOT.TD == 1
        
        % % ----- Resample experimental data to match the time-scale of DOT.time ----
        %        z = zeros(DOT.time.nstep,size(EXP.data.spc,2));
        %        f = zeros(DOT.time.nstep,size(EXP.data.spc,2));
        
        for i = 1:size(EXP.data.spc,2)
            z(:,i) = resampleSPC(EXP.data.spc(:,i),EXP.time.axis,DOT.time.dt);
            f(:,i) = resampleSPC(EXP.data.ref(:,i),EXP.time.axis,DOT.time.dt);
        end
        % ------ Comparison between conv(irf,fwd) and the experimental data -------
        if FORWARD == 1
            h = 1;
            nsub=numSubplots(DOT.radiometry.nL);
            for h=1:size(z,2)
                fh=figure(202);
                fh.Name = ['Ref measurement number ',num2str(h)];
                fh=figure(203);
                fh.Name = ['Data measurement number ',num2str(h)];
                for inl = 1:DOT.radiometry.nL
                    meas_set = (1:nmeas)+(inl-1)*nmeas;
                    figure(202);subplot(nsub(1),nsub(2),inl);
                    semilogy(1:size(z,1),[f(:,meas_set(h))./sum(f(:,meas_set(h)),'omitnan'),RefTD(:,meas_set(h))./sum(RefTD(:,meas_set(h)),'omitnan')]),
                    legend('Data','conv(irf,fwd)'),
                    title(['Wavelength :' num2str(DOT.radiometry.lambda(inl))]);
                    figure(203);subplot(nsub(1),nsub(2),inl);
                    semilogy(1:size(z,1),[z(:,meas_set(h))./sum(z(:,meas_set(h)),'omitnan'),DataTD(:,meas_set(h))./sum(DataTD(:,meas_set(h)),'omitnan')]),
                    legend('Data','conv(irf,fwd)'),
                    title(['Wavelength :' num2str(DOT.radiometry.lambda(inl))]);
                end
                ImageMeasure(DOT,h);
                pause(0.01);
            end
            %             for h = 1:size(z,2)
            %             semilogy(1:size(z,1),[f(:,h)./sum(f(:,h)),RefTD(:,h)./sum(RefTD(:,h),'omitnan')]),
            %             legend('Data','conv(irf,fwd)'),
            %             title(['Ref measurement number ',num2str(h)]);
            %             figure(204);semilogy(1:size(z,1),[z(:,h)./sum(z(:,h)),DataTD(:,h)./sum(DataTD(:,h),'omitnan')]),
            %             legend('Data','conv(irf,fwd)'),
            %             title(['Data measurement number ',num2str(h)]);
            %             ImageMeasure(DOT,h);
            %             pause(0.01)
            %             end
        end
        %figure;semilogy(f)
        %figure;semilogy(f./z)
        % -------- Save Time-resolved Data,Reference and standard deviation -------
        DataTD = z;
        RefTD = f;
        sdTD = sqrt(DataTD);    % Poisson
        clear z f sh inl meas_set
    else
        % ------------------------------ CW case ----------------------------------
        DataCW = sum(EXP.data.spc,1);
        RefCW = sum(EXP.data.ref,1);
        sdCW = sqrt(DataCW);    % Poisson
    end
end
%% ========================================================================
if RECONSTRUCTION == 1
    disp('-------------------------------------------------------------------')
    disp('                         RECONSTRUCTION                            ')
    disp('-------------------------------------------------------------------')
    %==========================================================================
    %%                             RECONSTRUCTION
    %==========================================================================
    % By default all the structure DOT is renamed in REC. Is it possible to
    % load a different mesh, create a different grid but, in this case, Qvec
    % and Mvec need to be regenerated.
    %==========================================================================
    %REC = DOT;clear DOT;
    for fn = fieldnames(DOT)'
        if ~isfield(REC,fn{1})
            REC.(fn{1}) = DOT.(fn{1});
        else
            for fn1 = fieldnames(DOT.(fn{1}))'
                if ~isfield(REC.(fn{1}),(fn1{1}))
                    REC.(fn{1}).(fn1{1}) = DOT.(fn{1}).(fn1{1});
                end
            end
        end
    end
    clear DOT
    
    
    %==========================================================================
    %%                               Load data
    %==========================================================================
    if (EXPERIMENTAL == 0)||(EXP_DATA == 0)
        load([rdir,filename,'_', 'Data'])
        if REF == 0
            REC.ref = 0;
        else
            REC.ref = RefCW;
        end
        REC.Data = DataCW;
        REC.sd = sdCW;
        clear DataCW sdCW
    elseif (REC.TD == 0)
        REC.ref = RefCW;
        REC.Data = DataCW;
        REC.sd = sdCW;
        clear DataCW sdCW
    end
    %==========================================================================
    %%                      Create time windows for TD
    %==========================================================================
    if strcmpi(REC.domain,'td')
        %     if NUM_TW > 0
        %         nmax = max(REC.time.nstep,numel(REC.time.irf.data));
        %         twin = CreateTimeWindows(REC.time.nstep,[REC.time.irf.ch0 - 1,nmax],'even',NUM_TW);
        %     else
        %         twin = CreateTimeWindows(REC.time.nstep,REC.time.roi,'even',diff(REC.time.roi)+1);
        %     end
        % select ROI dinamically if not defined in RecSettings_DOT.m
        if REC.solver.prejacobian.load
            warning('off','verbose')
            warning('off','backtrace')
            warning('Pre loaded ROI will be used. Any value will be overwritten');
            warning('Pre loaded NUM_TW will be used. Any value will be overwritten');
            warning('on','verbose')
            warning('on','backtrace')
            if strcmpi(REC.solver.type,'born')||...
                strcmpi(REC.solver.type,'born_spectral_post_proc')||...
                strcmpi(REC.solver.type,'usprior')
                load([REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(1)),'.mat'],'ROI','NW');
            else
                load(REC.solver.prejacobian.path,'ROI','NW');
            end
            REC.time.roi = ROI;
            NUM_TW = NW;
        end
        if isempty(REC.time.roi)
            for inl = 1:REC.radiometry.nL
                meas_set = (1:nmeas)+(inl-1)*nmeas;
                fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                REC.time.roi(inl,:) = SelectROI(DataTD(:,meas_set),REC.time.irf.data(:,inl));
            end
        end
        REC.ref = []; REC.Data = []; REC.sd = [];
        nsub = numSubplots(REC.radiometry.nL);
        figure(20);
        for inl = 1:REC.radiometry.nL
            meas_set = (1:nmeas)+(inl-1)*nmeas;
            twin_set = (1:2)+(inl-1)*2;
            fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
            REC.time.twin(:,twin_set) = CreateTimeWindows(REC.time.nstep,REC.time.roi(inl,:),'even',NUM_TW);
            REC.time.nwin(inl) = size(REC.time.twin(:,twin_set),1);
            subplot(nsub(1),nsub(2),inl)
            ShowTimeWindows(DataTD(:,meas_set),REC.time.twin(:,twin_set),REC.time.dt);
            title(['Wavelength: ' num2str(REC.radiometry.lambda(inl))]);
            %save_figure('time_ROI');
            drawnow;
            if exist('RefTD','var')
                REC.ref(:,meas_set) = WindowTPSF(RefTD(:,meas_set),REC.time.twin(:,twin_set));
                if inl == REC.radiometry.nL,  clear RefTD; end
            end
            REC.Data(:,meas_set) = WindowTPSF(DataTD(:,meas_set),REC.time.twin(:,twin_set));
            if (EXP_DATA == 0)
                if ~strcmpi(REC.time.noise,'none')
                    REC.sd(:,meas_set) = sqrt(WindowTPSF(sdTD(:,meas_set).^2,REC.time.twin(:,twin_set)));     % ATTENTION NOT TO SUM SD
                else
                    REC.sd(:,meas_set) = ones(size(REC.Data(:,meas_set)));
                end
            end
            if (EXP_DATA == 1)
                REC.sd(:,meas_set) = sqrt(WindowTPSF(sdTD(:,meas_set).^2,REC.time.twin(:,twin_set)));
            end
            if inl == REC.radiometry.nL
                clear DataTD sdTD
                tilefigs;
            end
        end
    end
    clear inl shift meas_set
    % % =========================================================================
    % %%                        Initial parameter estimates
    % % =========================================================================
    % ---------------------- Set nodes optical properties ---------------------
    REC.opt.mua = ones(REC.grid.N,1) .* REC.opt.mua0;
    REC.opt.musp = ones(REC.grid.N,1) .* REC.opt.musp0;
    REC.opt.n = ones(REC.grid.N,1) * REC.opt.nB;
    REC.opt.kap = 1./(3*(REC.opt.musp));
    REC.cm = 0.299/REC.opt.nB;
    
    %==========================================================================
    %%              PLOT REFERENCE VALUES IF PRESENT
    %==========================================================================
    % ------------------------ Reference mua ----------------------------------
    if ~SPECTRA
        for inl = 1:REC.radiometry.nL
            PlotMua = squeeze(REC.opt.Mua(:,:,:,inl));
            PlotMusp = squeeze(REC.opt.Musp(:,:,:,inl));
            fh=figure(300+inl);
            fh.Name = ['Wave: ' num2str(REC.radiometry.lambda(inl))];
            if isfield(REC.opt,'Mua')
                ShowRecResults(REC.grid,PlotMua,...
                    REC.grid.z1,REC.grid.z2,REC.grid.dz,1,...
                    'auto',0,max(PlotMua(:)));
                suptitle('Mua');
            end
            % ------------------------ Reference musp ---------------------------------
            fh=figure(400+inl);
            fh.Name = ['Wave: ' num2str(REC.radiometry.lambda(inl))];
            if isfield(REC.opt,'Musp')
                ShowRecResults(REC.grid,PlotMusp,...
                    REC.grid.z1,REC.grid.z2,REC.grid.dz,1,...
                    'auto',0,max(PlotMusp(:)));
                suptitle('Mus');
            end
        end
    end
    if SPECTRA
        mask = REC.opt.Mua(:,:,:,1)./REC.opt.muaB(1);
        for ic = (1:REC.spe.nCromo)
            if REC.opt.conc0(ic)==REC.opt.hete1.conc(ic)
                mask = ones(size(mask));
            else
                mask = REC.opt.Mua(:,:,:,1)./REC.opt.muaB(1);
            end
            REC.opt.Conc(:,:,:,ic) = mask.*REC.opt.concB(ic);
            fh=figure(900+ic);fh.NumberTitle = 'off';fh.Name = [REC.spe.cromo_label{ic} ' Map'];
            ShowRecResults(REC.grid,REC.opt.Conc(:,:,:,ic),...
                REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
            suptitle(REC.spe.cromo_label{ic});
        end
        if REC.spe.active_cromo(strcmpi(REC.spe.cromo_label,'hb'))
            REC.opt.HbTot = REC.opt.Conc(:,:,:,strcmpi(REC.spe.cromo_label,'hb'))+...
                REC.opt.Conc(:,:,:,strcmpi(REC.spe.cromo_label,'hbo2'));
            REC.opt.So2 = REC.opt.Conc(:,:,:,strcmpi(REC.spe.cromo_label,'hbo2'))./REC.opt.HbTot;
            fh=figure(900+ic+1);fh.NumberTitle = 'off';fh.Name = 'HbTot Map';
            ShowRecResults(REC.grid,REC.opt.HbTot,...
                REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
            suptitle('HbTot');
            fh=figure(900+ic+2);fh.NumberTitle = 'off';fh.Name = 'So2 Map';
            ShowRecResults(REC.grid,REC.opt.So2,...
                REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
            suptitle('So2');
        else
            REC.opt.HbTot = zeros(REC.grid.dim);REC.opt.So2 = zeros(REC.grid.dim);
        end
        mask = REC.opt.Musp(:,:,:,1)./REC.opt.muspB(1);
        if REC.opt.aB==REC.opt.hete1.a, mask = ones(size(mask)); end
        REC.opt.a = mask.*REC.opt.aB;
        fh=figure(900+ic+3);fh.NumberTitle = 'off';fh.Name = ('a Map');
        ShowRecResults(REC.grid,REC.opt.a,...
            REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
        suptitle('a');
        mask = REC.opt.Musp(:,:,:,1)./REC.opt.muspB(1);
        if REC.opt.bB==REC.opt.hete1.b, mask = ones(size(mask)); end
        REC.opt.b = mask.*REC.opt.bB;
        fh=figure(900+ic+4);fh.NumberTitle = 'off';fh.Name = ('b Map');
        ShowRecResults(REC.grid,REC.opt.b,...
            REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
        suptitle('b');
    end
    drawnow
    % =========================================================================
    %%                      Reconstruction Solver
    % =========================================================================
    switch lower(REC.domain)
        case 'cw'
            switch lower(REC.solver.type)
                case 'gn'
                    
                case 'lm'
                    
                case 'cg'
                    
                case 'born'
                    %    [REC.opt.bmua,REC.opt.bmusp] = RecSolverBORN_CW(REC.solver,...
                    %        REC.grid,REC.opt.mua0,REC.opt.musp0,REC.A,...
                    %        REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.Data,REC.ref,0);
                    
                case 'fit'
                    %    [REC.opt.bmua,REC.opt.bmusp] = FitMuaMusAmpl_CW(REC.solver,REC.mesh.hMesh,...
                    %        REC.grid.hBasis,REC.opt.notroi,...
                    %        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,...
                    %        REC.Q,REC.M,REC.dmask,REC.Data,REC.ref,0);
                    
                    
            end
        case 'td'
            switch lower(REC.solver.type)
                case 'gn'
                    [REC.opt.bmua,REC.opt.bmusp,REC.erri] = RecSolverGN_TD(REC.solver,...
                        REC.grid,mesh.hMesh,...
                        REC.grid.hBasis,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,...
                        mesh.qvec,mesh.mvec,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm, REC.Data,...
                        REC.time.irf.data,REC.ref,REC.sd,0);
                case 'gn-usprior'
                    if ~isempty(REC.solver.prior.path)
                        REC.solver.prior.refimage = ...
                            priormask3D(REC.solver.prior.path,REC.grid);
                    else
                        disp('No prior is provided in RECSettings_DOT. Reference mua will be used');
                        REC.solver.prior.refimage = REC.opt.Mua;
                    end
                    REC.solver.prior.refimage = REC.opt.Mua;
                    [REC.opt.bmua,REC.opt.bmusp,REC.erri] = RecSolverGN_TD(REC.solver,...
                        REC.grid,mesh.hMesh,...
                        REC.grid.hBasis,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,...
                        mesh.qvec,mesh.mvec,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm, REC.Data,...
                        REC.time.irf.data,REC.ref,REC.sd,0);
                case 'born'
                    REC.solver.prior.refimage = [];
                    original_path = REC.solver.prejacobian.path;
                    for inl = 1:REC.radiometry.nL
                        if REC.solver.prejacobian.load
                            REC.solver.prejacobian.path=strcat(original_path,num2str(REC.radiometry.lambda(inl)),'.mat');
                        end
                        fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                        [REC.opt.bmua(:,inl),REC.opt.bmusp(:,inl)] = RecSolverBORN_TD(REC.solver,...
                            REC.grid,...
                            REC.opt.mua0(inl),REC.opt.musp0(inl),REC.opt.nB,REC.A,...
                            REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                            REC.time.twin(:,(1:2)+(inl-1)*2),REC.time.self_norm,REC.Data(:,(1:nmeas)+(inl-1)*nmeas),...
                            REC.time.irf.data(:,inl),REC.ref(:,(1:nmeas)+(inl-1)*nmeas),REC.sd(:,(1:nmeas)+(inl-1)*nmeas),REC.type_fwd);
                        if REC.solver.prejacobian.load==0
                            movefile([REC.solver.prejacobian.path,'.mat'],strcat(REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'));
                        end
                        drawnow
                    end
                    REC.solver.prejacobian.path = original_path; clear original_path
                case 'born_spectral_post_proc'
                    REC.solver.prior.refimage = [];
                    original_path = REC.solver.prejacobian.path;
                    for inl = 1:REC.radiometry.nL
                        if REC.solver.prejacobian.load
                            REC.solver.prejacobian.path=strcat(original_path,num2str(REC.radiometry.lambda(inl)),'.mat');
                        end
                        fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                        [REC.opt.bmua(:,inl),REC.opt.bmusp(:,inl)] = RecSolverBORN_TD(REC.solver,...
                            REC.grid,...
                            REC.opt.mua0(inl),REC.opt.musp0(inl),REC.opt.nB,REC.A,...
                            REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                            REC.time.twin(:,(1:2)+(inl-1)*2),REC.time.self_norm,REC.Data(:,(1:nmeas)+(inl-1)*nmeas),...
                            REC.time.irf.data(:,inl),REC.ref(:,(1:nmeas)+(inl-1)*nmeas),REC.sd(:,(1:nmeas)+(inl-1)*nmeas),REC.type_fwd);
                        if REC.solver.prejacobian.load==0
                            movefile([REC.solver.prejacobian.path,'.mat'],strcat(REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'));
                        end
                        drawnow
                    end
                    REC.solver.prejacobian.path = original_path; clear original_path
                    [bmua,bmusp,REC.opt.bConc,REC.opt.bA,REC.opt.bbB]=FitVoxel(REC.opt.bmua,REC.opt.bmusp,REC.spe);
                    REC.opt.bmua = []; REC.opt.bmua = bmua;REC.opt.bmusp = []; REC.opt.bmusp = bmusp;
                case 'spectral_born'
                    REC.solver.prior.refimage = [];
                    [REC.opt.bmua,REC.opt.bmusp,REC.opt.bConc,REC.opt.bA,REC.opt.bbB] = RecSolverBORN_TD_spectral(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm,REC.Data,...
                        REC.time.irf.data,REC.ref,REC.sd,REC.type_fwd,REC.radiometry,REC.spe);
                case 'spectral_usprior'
                    if ~isempty(REC.solver.prior.path)
                        REC.solver.prior.refimage = ...
                            priormask3D(REC.solver.prior.path,REC.grid);
                    else
                        disp('No prior is provided in RECSettings_DOT. Reference mua will be used');
                        REC.solver.prior.refimage = REC.opt.Mua;
                    end
                    [REC.opt.bmua,REC.opt.bmusp,REC.opt.bConc,REC.opt.bA,REC.opt.bbB] = RecSolverBORN_TD_USPrior_spectral(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm,REC.Data,...
                        REC.time.irf.data,REC.ref,REC.sd,REC.type_fwd,REC.radiometry,REC.spe);
                case 'born2reg'
                    if ~isempty(REC.solver.prior.path)
                        REC.solver.prior.refimage = ...
                            priormask3D(REC.solver.prior.path,REC.grid);
                    else
                        disp('No prior is provided in RECSettings_DOT. Reference mua will be used');
                        REC.solver.prior.refimage = REC.opt.Mua;
                    end
                    [REC.opt.bmua,REC.opt.bmusp] = RecSolverBORN_TD(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm,REC.Data,...
                        REC.time.irf.data,REC.ref,REC.sd,0);
                    
                    
                    
                    %% @Simon: US prior inversion
                case 'usprior'
                    if ~isempty(REC.solver.prior.path)
                        REC.solver.prior.refimage = ...
                            priormask3D(REC.solver.prior.path,REC.grid);
                    else
                        disp('No prior is provided in RECSettings_DOT. Reference mua will be used');
                        REC.solver.prior.refimage = REC.opt.Mua;
                    end
                    REC.solver.prior.refimage =  double(REC.solver.prior.refimage)*10 + 0.1;
                    original_path = REC.solver.prejacobian.path;
                    for inl = 1:REC.radiometry.nL
                        if REC.solver.prejacobian.load
                            REC.solver.prejacobian.path=strcat(original_path,num2str(REC.radiometry.lambda(inl)),'.mat');
                        end
                        fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                        [REC.opt.bmua(:,inl),REC.opt.bmusp(:,inl)] = RecSolverBORN_TD_USPrior(REC.solver,...
                            REC.grid,...
                            REC.opt.mua0(inl),REC.opt.musp0(inl),REC.opt.nB,REC.A,...
                            REC.Source.Pos,REC.Detector.Pos,REC.dmask,...
                            REC.time.dt,REC.time.nstep,REC.time.twin(:,(1:2)+(inl-1)*2),...
                            REC.time.self_norm,REC.Data(:,(1:nmeas)+(inl-1)*nmeas),...
                            REC.time.irf.data(:,inl),REC.ref(:,(1:nmeas)+(inl-1)*nmeas),REC.sd(:,(1:nmeas)+(inl-1)*nmeas),REC.type_fwd);
                        if REC.solver.prejacobian.load==0
                            movefile([REC.solver.prejacobian.path,'.mat'],strcat(REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'));
                        end
                        drawnow
                    end
                    REC.solver.prejacobian.path = original_path; clear original_path
                    %%
                case 'fit'
                    [REC.opt.bmua,REC.opt.bmusp] = FitMuaMus_TD(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,...
                        REC.time.dt,REC.time.nstep,REC.time.twin,...
                        REC.time.self_norm,REC.Data,...
                        REC.time.irf.data,REC.ref,REC.sd,1);
                case 'spectral_fit'
                    [REC.opt.bmua,REC.opt.bmusp] = SpectralFitMuaMus_TD(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,...
                        REC.time.dt,REC.time.nstep,REC.time.twin,...
                        REC.time.self_norm,REC.Data,...
                        REC.time.irf.data,REC.ref,REC.sd,1,REC.radiometry,REC.spe);
                case 'cg'
                    
                case 'l1'
                    if ~isempty(REC.solver.prior)
                        REC.solver.prior = REC.opt.Mua;
                    end
                    [REC.opt.bmua,REC.opt.bmusp] = RecSolverL1_TD(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm,REC.Data,...
                        REC.time.irf.data,REC.ref,REC.sd,0);
                    
                case 'fit4param' %% you require a TOAST installation
                    % check if TOAST is correctly installed
                    if ~isempty(REC.solver.prior.path)
                        REC.solver.prior.refimage = ...
                            priormask3D(REC.solver.prior.path,REC.grid, REC.solver.type);
                    else
                        disp('No prior is provided in RECSettings_DOT. Reference mua will be used');
                        REC.solver.prior.refimage = REC.opt.Mua;
                        
                    end
                    [REC.opt.bmua,REC.opt.bmusp, REC.opt.fitOUTPUT] = Fit2Mua2Mus_TD(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,[],...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,...
                        REC.time.dt,REC.time.nstep,REC.time.twin,...
                        REC.time.self_norm,REC.Data,...
                        REC.time.irf.data,REC.ref,REC.sd,1);
                    
            end
    end
    ROI = REC.time.roi; %#ok<NASGU>
    NW = NUM_TW; %#ok<NASGU>
    if strcmpi(REC.solver.type,'born')||...
            strcmpi(REC.solver.type,'born_spectral_post_proc')||...
            strcmpi(REC.solver.type,'usprior')
        for inl = 1:REC.radiometry.nL
            if ~exist([REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'],'file')
                save([REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'],'ROI','NW');
            else
                save([REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'],'ROI','NW','-append');
            end
        end
    else
        if ~exist([REC.solver.prejacobian.path '.mat'],'file')
            save(REC.solver.prejacobian.path,'ROI','NW');
        else
            save(REC.solver.prejacobian.path,'ROI','NW','-append');
        end
    end
    clear ROI NW
    % ---------------------------- display mua --------------------------------
    if ~contains(REC.solver.type,'fit')
        drawnow;
        if ~contains(REC.solver.type,'spectral')
            for inl = 1:REC.radiometry.nL
                PlotMua = reshape(REC.opt.bmua,[REC.grid.dim REC.radiometry.nL]);
                PlotMus = reshape(REC.opt.bmusp,[REC.grid.dim REC.radiometry.nL]);
                PlotMua = PlotMua(:,:,:,inl);
                PlotMus = PlotMus(:,:,:,inl);
                fh=figure(500+inl);fh.NumberTitle = 'off';fh.Name = ['Recon Mua. Wave: ' num2str(REC.radiometry.lambda(inl))];
                ShowRecResults(REC.grid,PlotMua,...
                    REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto',0.00,0.05);
                suptitle('Recon Mua');
                % ---------------------------- display musp -------------------------------
                fh=figure(600+inl);fh.NumberTitle = 'off';fh.Name = ['Recon Mus. Wave: ' num2str(REC.radiometry.lambda(inl))];
                ShowRecResults(REC.grid,PlotMus,...
                    REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
                suptitle('Recon Mus');
                
                drawnow;
                tilefigs;
                disp('recon: finished')
                fh = figure(700+inl);fh.NumberTitle = 'off';fh.Name = ['PlotHete. Wave: ' num2str(REC.radiometry.lambda(inl))];
                subplot(1,2,1),PlotHeteQM(REC,PlotMua,REC.opt.mua0(inl)),
                title('Recon Mua');
                subplot(1,2,2),PlotHeteQM(REC,PlotMus,REC.opt.musp0(inl)),
                title('Recon Mus');
                
                drawnow;
            end
        end
        if contains(REC.solver.type,'spectral')
            for ic = 1:REC.spe.nCromo
                Conc=reshape(REC.opt.bConc(:,ic),REC.grid.dim);
                fh=figure(800+ic);fh.NumberTitle = 'off';fh.Name = ['Recon ' REC.spe.cromo_label{ic} ' Map'];
                ShowRecResults(REC.grid,Conc,...
                    REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
                suptitle(['Recon ' REC.spe.cromo_label{ic}]);
            end
            if REC.spe.active_cromo(strcmpi(REC.spe.cromo_label,'hb'))
                REC.opt.HbTot = reshape(REC.opt.bConc(:,strcmpi(REC.spe.cromo_label,'hb')),REC.grid.dim)+...
                    reshape(REC.opt.bConc(:,strcmpi(REC.spe.cromo_label,'hbo2')),REC.grid.dim);
                REC.opt.So2 = reshape(REC.opt.bConc(:,strcmpi(REC.spe.cromo_label,'hbo2')),REC.grid.dim)./REC.opt.HbTot;
                fh=figure(800+ic+1);fh.NumberTitle = 'off';fh.Name = 'HbTot Map';
                ShowRecResults(REC.grid,REC.opt.HbTot,...
                    REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
                suptitle('HbTot');
                fh=figure(800+ic+2);fh.NumberTitle = 'off';fh.Name = 'So2 Map';
                ShowRecResults(REC.grid,REC.opt.So2,...
                    REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
                suptitle('So2');
            else
                REC.opt.HbTot = zeros(REC.grid.dim);REC.opt.So2 = zeros(REC.grid.dim);
            end
            fh=figure(800+ic+3);fh.NumberTitle = 'off';fh.Name = ('Recon a Map');
            ShowRecResults(REC.grid,reshape(REC.opt.bA,REC.grid.dim),...
                REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
            suptitle('Recon a');
            fh=figure(800+ic+4);fh.NumberTitle = 'off';fh.Name = ('Recon b Map');
            ShowRecResults(REC.grid,reshape(REC.opt.bbB,REC.grid.dim),...
                REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
            suptitle('Recon b');
        end
        
    end
    %% Save reconstruction
    if ~contains(REC.solver.type,'fit')
        disp(['Reconstruction will be stored in: ', rdir,filename,'_', 'REC.mat']);
        save([rdir,filename,'_', 'REC'],'REC');
    end
    % =========================================================================
    %%                            Quantify DOT
    % =========================================================================
    if ~contains(REC.solver.type,'fit')
        Q = QuantifyDOT(REC,~EXP_DATA);
    end
    Quantification_MultiSim
    % =========================================================================
    %%                            Save quantification
    % =========================================================================
    if ~contains(REC.solver.type,'fit')
        disp(['Quantification will be stored in: ', rdir,filename,'_', 'REC.mat']);
        save([rdir,filename,'_', 'REC'],'Q','-append')
    end
    clearvars REC
end
