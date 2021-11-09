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
%% correct missing flags for compatibility with older examples
if ~exist('SHOWPLOTS','var'), SHOWPLOTS = 0; end
if ~exist('REMOVE_VOXELS','var'), REMOVE_VOXELS = 0; end
if ~exist('SHOWIMAGEMEAS','var'), SHOWIMAGEMEAS = 0; end
if ~exist('RECFILE_APPEND','var'), RECFILE_APPEND = '';end
if ~exist('RECFILE_APPEND_DATA','var'), RECFILE_APPEND_DATA = '';end
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
%% set suffixes to filename
filename_DATA = [filename, RECFILE_APPEND_DATA];
filename = [filename, RECFILE_APPEND];

%% Get True Variables from X
if 0 == 0
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
end

if contains(REC.solver.type,'spectral') && ~exist('conc_','var')% if concentrations are not given, fit voxels
    
    
    %[~,~,DOT.opt.concB, DOT.opt.aB, DOT.opt.bB] = 
    [~,~,conc_,a_,b_]=FitVoxel(Xd{1},Xd{2},DOT.spe);
    Xd = {conc_,[a_ b_]};
    [~,~,conc_,a_,b_]=FitVoxel(Xp{1},Xp{2},DOT.spe);
    Xp = {conc_,[a_ b_]};
    [~,~,conc_,a_,b_]=FitVoxel(Xr{1},Xr{2},DOT.spe);
    Xr = {conc_,[a_ b_]};


    
end


[DOT,REC] = TranslateX_2Var(SPECTRA,Xd,Xp,Xr,spectra_file,DOT,REC);
%% check dimensions of dmask for compatibily purposes
if ismatrix(DOT.dmask)
    DOT.dmask = repmat(DOT.dmask, [1,1,DOT.radiometry.nL]);
end

%==========================================================================
%%                          Experimental data
% Load the structure EXP containing all the data for analyzing experimental
% data. This structure is created with the routine 'CreateExperiment.m'.
%==========================================================================

if (EXPERIMENTAL == 1)
    load([exp_path,'EXP_',exp_file])

    if (EXP_DATA == 1)
        % Select measurements by dmask
        DOT.dmask = logical(DOT.dmask.*EXP.grid.dmask(:,:,:));
        if (size(EXP.data.spc,2)>sum(DOT.dmask(:)))
            
            % do not select data that are not available experimentally      
            EXP.data.spc = EXP.data.spc(:,DOT.dmask(:));
            EXP.data.ref = EXP.data.ref(:,DOT.dmask(:));
            if DOT.time.self_norm
                indexmeas = findMeasIndex(DOT.dmask);
                for inl = 1:DOT.radiometry.nL
                    meas_set = indexmeas{inl};
                    nm = numel(meas_set);
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
[DOT.opt.Mua, DOT.opt.Musp,DOT.opt.Conc, DOT.opt.A, DOT.opt.B] = applyGrid(DOT.grid,DOT.opt.muaB,DOT.opt.muspB,REC.solver.type,DOT.opt.concB,DOT.opt.aB,DOT.opt.bB);
%==========================================================================
%%                      Set Heterogeneities
%==========================================================================
%--------------------------- INCLUSIONS ---------------------------------%
for i = 1:NUM_HETE
    h_string = ['hete',num2str(i)];
    [DOT,DOT.opt.(h_string)] = setHete(DOT,DOT.opt.(h_string),REC.solver.type);
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
        % check dimensions of IRF for compatibily purposes
        if size(EXP.irf.data,2) == 8
            % push waveleght dimension as last one
            tmp = permute(EXP.irf.data,[1,3,2]);
            tmp = repmat(tmp, [1, size(DOT.dmask,1)*size(DOT.dmask,2) ,1]);
            % save old irf for debug
            EXP.irf.dataMat = EXP.irf.data;
            % update IRF
            EXP.irf.data = reshape(tmp, [size(tmp,1),size(tmp,2)*size(tmp,3)]);
            % update peak positions 
            EXP.irf.peak = repmat(permute(EXP.irf.peak, [1,2]),[size(DOT.dmask,1)*size(DOT.dmask,2) ,1]);
            % update baric
            EXP.irf.baric = repmat(permute(EXP.irf.baric, [1,2]),[size(DOT.dmask,1)*size(DOT.dmask,2) ,1]);
            clear tmp
        end
        switch lower(EXP_DELTA)
            case 'peak'
                EXP.irf.data = zeros(size(EXP.irf.data));
                peak_pos = reshape([EXP.irf.peak.pos], size(EXP.irf.peak));
                [I,J]  = ind2sub(size(peak_pos),[1:numel(peak_pos)]);
                idx = sub2ind(size(EXP.irf.data),peak_pos(:)' ,J,I);
                EXP.irf.data(idx) = 1;
            case 'baric'                
                EXP.irf.data = zeros(size(EXP.irf.data));
                baric_pos = reshape([EXP.irf.baric.pos], size(EXP.irf.baric));
                [I,J]  = ind2sub(size(baric_pos),[1:numel(baric_pos)]);
                idx = sub2ind(size(EXP.irf.data),baric_pos(:)' ,I,J);
                EXP.irf.data(idx) = 1; 
            case 'peak0'
                EXP.irf.data = zeros(1,size(EXP.irf.data,2));
                EXP.irf.data(1,:) = 1;
        end
        
        tmp = EXP.irf.data;
        % set nan to 0
        tmp(isnan(tmp)) = 0;
        tmp = reshape(tmp, [size(tmp,1),size(tmp,2)*size(tmp,3)]);
        disp('experimental irf resampling')
        tic
        for i_sam = 1:size(tmp,2)          
            tmp_(:,i_sam)= resampleSPC(tmp(:,i_sam),EXP.time.axis,DOT.time.dt,'norm');
        end
       toc
        DOT.time.irf.data = reshape(tmp_, [size(tmp_,1),size(EXP.irf.data,2)*size(EXP.irf.data,3)]);
        % get set of IRFs used in convolution
        DOT.time.irf.data_mask = DOT.time.irf.data(:, DOT.dmask(:));
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
if SHOWPLOTS
    % plot DMASK
    figure(10) 
    for inl = 1:DOT.radiometry.nL
        subplot(1,size(DOT.dmask,3), inl);
        imagesc(DOT.dmask(:,:,inl)),xlabel('Sources'),ylabel('Meas'), axis image;
        title(sprintf('Dmask: lam %g',  DOT.radiometry.lambda(inl)));
    end
    % plot source-detectors and heterogeneities
    figure(11)
    subplot(1,2,1),PlotHeteQM(DOT,squeeze(DOT.opt.Mua(:,:,:,1)),DOT.opt.muaB(1)),title('Mua');
    subplot(1,2,2),PlotHeteQM(DOT,squeeze(DOT.opt.Musp(:,:,:,1)),DOT.opt.muspB(1)),title('Musp');
    drawnow;
end
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
        [MuaB,MuspB] = applyGrid(DOT.grid,DOT.opt.muaB,DOT.opt.muspB,REC.solver.type,DOT.opt.concB,DOT.opt.aB,DOT.opt.bB);
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
                    DOT.opt.muaB, DOT.opt.muspB,DOT.opt.nB,...
                    DOT.opt.Mua,DOT.opt.Musp,...
                    DOT.A, DOT.time.dt,...
                    length(DOT.time.time), DOT.time.self_norm, geom, TYPE_FWD,DOT.radiometry);
                save([rdir,filename,'_', 'FwdTeo'],'DataTD');
            if REF == 1
                [MuaB,MuspB] = applyGrid(DOT.grid,DOT.opt.muaB,DOT.opt.muspB,REC.solver.type,DOT.opt.concB,DOT.opt.aB,DOT.opt.bB);
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
            
                % convolve by IRF
                DataTD = ConvIRF(full(DataTD),...
                    DOT.time.irf.data_mask);
       
                
                if REF == 1
                    % convolve by IRF
                    RefTD = ConvIRF(full(RefTD),...
                            DOT.time.irf.data_mask);
                end
                clear nmax
                %% plots
                if SHOWPLOTS
                    idx = findMeasIndex(DOT.dmask);
                    tmpdata = DataTD;
                    for inl = 1:DOT.radiometry.nL
                        nsub=numSubplots(DOT.radiometry.nL);
                        figure(12);    
                        subplot(nsub(1),nsub(2),inl)
                        semilogy(tmpdata(:,idx{inl}(1)),'b'),hold,semilogy(DataTD(:,idx{inl}(1)),'k'),
                        semilogy(DOT.time.irf.data_mask(:,idx{inl}(1)),'r'),%ylim([max(DataTD(:))/10000 max(DataTD(:))])
                        title(['Wavelength: ' num2str(DOT.radiometry.lambda(inl))]);
                    end
                    drawnow
                end

        end

        
        % ---- Radiometry -------
        if RADIOMETRY == 1
            RealFactor = Radiometry(DOT.radiometry);
        else
            RealFactor = ones(1,DOT.radiometry.nL);
        end
    end
   %-------------------- Add noise to TD data
    %------------------------
    sdTD = zeros(size(DataTD));
    idx = findMeasIndex(DOT.dmask);
    for inl = 1:DOT.radiometry.nL
        meas_set = idx{inl};
        nmeas = numel(meas_set);
        DataTD_single_wave = DataTD(:,meas_set);
        if REF == 1
        RefTD_single_wave = RefTD(:,meas_set);
        end
        sdTD_single_wave = ones(size(DataTD_single_wave));
        if ~strcmpi(DOT.time.noise,'none')
            [DataTD_single_wave,~] = AddNoise(DataTD_single_wave,'gaussian',DOT.time.sigma);
            if REF == 1
                RefTD_single_wave = AddNoise(RefTD_single_wave,'gaussian',DOT.time.sigma);
            end
        end

        if strcmpi(DOT.time.noise,'poisson') && FORWARD==1

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
    disp(['Results will be stored in: ', rdir,filename_DATA,'_','Data.m'])
    %-------------------------- Save simulated data ---------------------------
    if FORWARD == 1
        save([rdir,filename_DATA,'_', 'Data'],'DataCW','sdCW');
        if REF == 1
            save([rdir,filename_DATA,'_', 'Data'],'RefCW','-append');
        end
        if DOT.TD == 1
            nstep = DOT.time.nstep;
            dt =  DOT.time.dt;
            save([rdir,filename_DATA,'_', 'Data'],'DataTD','sdTD','dt','nstep','-append');
            if REF == 1
                save([rdir,filename_DATA,'_', 'Data'],'RefTD','-append');
            end
        end
        
    end
end
drawnow
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
            tmpspc = EXP.data.spc(:,i);
            tmpspc(isnan(tmpspc)) = 0;
            tmpref = EXP.data.ref(:,i);
            tmpref(isnan(tmpref)) = 0;
            z(:,i) = resampleSPC(tmpspc,EXP.time.axis,DOT.time.dt);
            f(:,i) = resampleSPC(tmpref,EXP.time.axis,DOT.time.dt);
        end
        % ------ Comparison between conv(irf,fwd) and the experimental data -------
        if FORWARD == 1
            if SHOWPLOTS && SHOWIMAGEMEAS
                h = 1;
                nsub=numSubplots(DOT.radiometry.nL);
                fh=figure(202);
                fh.Name = ['Ref measurement number ',num2str(h)];
                fh=figure(203);
                fh.Name = ['Data measurement number ',num2str(h)];
                idxMeas = findMeasIndex(DOT.dmask);
                for inl = 1:DOT.radiometry.nL
                    meas_set = idxMeas{inl};
                    for h=1:size(meas_set,2)                
                        figure(202);subplot(nsub(1),nsub(2),inl);
                        semilogy(1:size(z,1),[f(:,meas_set(h))./sum(f(:,meas_set(h)),'omitnan'),RefTD(:,meas_set(h))./sum(RefTD(:,meas_set(h)),'omitnan')]),
                        legend('Data','conv(irf,fwd)'),
                        title(['Wavelength :' num2str(DOT.radiometry.lambda(inl))]);
                        figure(203);subplot(nsub(1),nsub(2),inl);
                        semilogy(1:size(z,1),[z(:,meas_set(h))./sum(z(:,meas_set(h)),'omitnan'),DataTD(:,meas_set(h))./sum(DataTD(:,meas_set(h)),'omitnan')]),
                        legend('Data','conv(irf,fwd)'),
                        title(['Wavelength :' num2str(DOT.radiometry.lambda(inl))]);
                        hMeas =findSQcouple(DOT.dmask(:,:,inl), h, 'sub');
                        ImageMeasure(DOT,hMeas(2),hMeas(1) ,DOT.dmask(:,:,inl));
                        pause(0.01);
                    end
                    

                end
            end

        end

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
    %% delete Toast files
    if isfield(DOT.grid,'hBasis') && RECONSTRUCTION == 0
        DOT.grid.hBasis.delete();
        DOT.grid.hBasis = [];
        DOT.grid = rmfield(REC.grid,'hBasis');
    end
    clear DOT
    %% update mask
    lamda_idMask = (1:8).*squeeze(sum(sum(REC.dmask,1),2)~=0)';
    lamda_id = sort(intersect(lamda_id, lamda_idMask), 'ascend');
    REC.dmask = REC.dmask(:,:,lamda_id);
    if numel(REC.radiometry.lambda)>numel(lamda_id)
        REC.radiometry.lambda = REC.radiometry.lambda(lamda_id);
    end
    REC.radiometry.nL = numel(REC.radiometry.lambda);


    %==========================================================================
    %%                               Load data
    %==========================================================================
    if (EXPERIMENTAL == 0)||(EXP_DATA == 0)
        load([rdir,filename_DATA,'_', 'Data'])
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
    
    %% resample DT data if time step or nstep are different
    if exist('nstep','var') && exist('nstep','var') %% This will need to change in a future version
        if strcmpi(REC.domain,'td') && (nstep ~= REC.time.nstep || dt ~= REC.time.dt )
            taxis = 0:dt:(dt*(nstep-1));

            tmpspc = DataTD;
            tmpspc(isnan(tmpspc)) = 0;
            tmpref = RefTD;
            tmpref(isnan(tmpref)) = 0;
            tmpsd = sdTD;
            tmpsd(isnan(sdTD)) = 0;

            disp('Resampling Data')
            tic
            for i = 1:size(DataTD,2)
                tmpspc_(:,i) = resampleSPC(tmpspc(:,i),taxis,REC.time.dt);
                tmpref_(:,i) = resampleSPC(tmpref(:,i),taxis,REC.time.dt);
                tmpsd_(:,i) = resampleSPC(tmpsd(:,i),taxis,REC.time.dt);
            end
            toc
            DataTD = tmpspc_;
            RefTD = tmpref_;
            sdTD = tmpsd_;
            clear tmpspc_ tmpref_ tmpsd_ tmpspc tmpref tmpsd
        end
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
            if strcmpi(REC.solver.type,'born')||...
                    strcmpi(REC.solver.type,'born_spectral_post_proc')||...
                    strcmpi(REC.solver.type,'usprior')||strcmpi(REC.solver.type,'tk1')
                load([REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(1)),'.mat'],'ROI','NW');
            else
                load([REC.solver.prejacobian.path '.mat'],'ROI','NW');
            end
            if ~exist('ROI','var')

                warning('Loading of ROI and NW FAILED. Not present in J');
            end
            warning('on','verbose')
            warning('on','backtrace')
        end
        if isempty(REC.time.roi)
            idxmeas = findMeasIndex(REC.dmask);
            for inl = 1:REC.radiometry.nL
                meas_set = idxmeas{inl};
                fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                REC.time.roi(inl,:) = SelectROI(DataTD(:,meas_set),REC.time.irf.data_mask(:,meas_set));
            end
        elseif strcmpi(REC.time.roi, 'auto')
            % selects the ROI automatically for each curve
            REC.time.roi = zeros(2, size(DataTD,2));
            for i = 1:size(DataTD,2)
                REC.time.roi(:,i) = auto_selectROI(DataTD(:,i),REC.time.roiRange1,REC.time.roiRange2);
            end
        elseif strcmpi(REC.time.roi, 'auto_by_wave')
            % write routing for selecting one ROI automaticaly per wavelength 
            REC.time.roi = zeros(2, size(DataTD,2));
            if ~isfield(REC.time, 'roiRange1') || ~isfield(REC.time, 'roiRange2') 
                REC.time.roiRange1 = [];
                REC.time.roiRange2 = [];
            end
            for i = 1:size(DataTD,2)
                REC.time.roi(:,i) = auto_selectROI(DataTD(:,i),REC.time.roiRange1,REC.time.roiRange2);
            end
            idxmeas = findMeasIndex(REC.dmask);
            for inl = 1:REC.radiometry.nL
                meas_set = idxmeas{inl};
                roi0 = REC.time.roi(1,meas_set);
                roi1 = REC.time.roi(2,meas_set);
                REC.time.roi(1,meas_set) = min(roi0(:));
                REC.time.roi(2,meas_set) = max(roi1(:));
            end
        end
        % foresee one roi per curve by repeating the twin

        tmp_roi = zeros(2,sum(REC.dmask(:)));
        idxmeas = findMeasIndex(REC.dmask);
        if any(size(REC.time.roi) ~= [2,sum(REC.dmask(:))]) % if not the correct dimensions go for a repmat
            for inl = 1:REC.radiometry.nL
                meas_set = idxmeas{inl};
                nmeas = numel(meas_set);
                tmp_roi(:,meas_set) =  repmat(REC.time.roi(inl, :)', [1,nmeas]);
            end
            REC.time.roi = tmp_roi;
            clear tmp_roi
        end

        
        REC.ref = []; REC.Data = []; REC.sd = [];
        
        % set time windows
        REC.time.twin = zeros(NUM_TW,2,size(DataTD,2));
        for i_twin = 1:size(REC.time.roi,2)
            SQ = findSQcouple(REC.dmask, i_twin, 'sub');
            fprintf('lambda: %g , source: %g , detector: %g ->',REC.radiometry.lambda(SQ(3)), SQ(2), SQ(1));
            REC.time.twin(:,:,i_twin) = CreateTimeWindows(REC.time.nstep,REC.time.roi(:,i_twin)','even',NUM_TW);
            REC.time.nwin(i_twin) = size(REC.time.twin(:,:,i_twin),1);
        end
        if exist('RefTD','var')
            REC.ref = WindowTPSF(RefTD,REC.time.twin);
            %clear RefTD
        end
        
        REC.Data = WindowTPSF(DataTD,REC.time.twin);
        
        
        
        if (EXP_DATA == 0)
            if ~strcmpi(REC.time.noise,'none')
                REC.sd = sqrt(WindowTPSF(sdTD.^2,REC.time.twin));     % ATTENTION NOT TO SUM SD
            else
                REC.sd = ones(size(REC.Data));
            end
%         elseif (EXP_DATA == 1)
%             REC.sd = sqrt(WindowTPSF(sdTD.^2,REC.time.twin));
        end        
        REC.sd(REC.sd<0 ) = 0;
        
        
        % normalise data after windowing
        if REC.time.self_norm 
            Area = sum(REC.Data,'omitnan');
            REC.Data = NormalizeTPSF(REC.Data);
            
            REC.sd = sqrt(REC.Data * spdiags(1./Area',0,sum(REC.dmask(:)),sum(REC.dmask(:))));
           if REF
               REC.ref = NormalizeTPSF(REC.ref);
           end
        end
        
        
        
        if SHOWPLOTS
            idxmeas = findMeasIndex(REC.dmask);
            nsub = numSubplots(REC.radiometry.nL);
            for inl = 1:REC.radiometry.nL
                figure(20);
                meas_set = idxmeas{inl};
                if ~isempty(meas_set)
                twin_set = meas_set;
                fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                subplot(nsub(1),nsub(2),inl)
                ShowTimeWindows(DataTD(:,meas_set),REC.time.twin(:,:,twin_set(1)),REC.time.dt);
                title(['Wavelength: ' num2str(REC.radiometry.lambda(inl))]);
                %save_figure('time_ROI');
                if exist('RefTD','var')
                    REC.ref(:,meas_set) = WindowTPSF(RefTD(:,meas_set),REC.time.twin(:,:,twin_set));
                end
                REC.Data(:,meas_set) = WindowTPSF(DataTD(:,meas_set),REC.time.twin(:,:,twin_set));
                if (EXP_DATA == 0)
                    if ~strcmpi(REC.time.noise,'none')
                        REC.sd(:,meas_set) = sqrt(WindowTPSF(sdTD(:,meas_set).^2,REC.time.twin(:,:,twin_set)));     % ATTENTION NOT TO SUM SD
                    else
                        REC.sd(:,meas_set) = ones(size(REC.Data(:,meas_set)));
                    end
                end
                if (EXP_DATA == 1)
                    REC.sd(:,meas_set) = sqrt(WindowTPSF(sdTD(:,meas_set).^2,REC.time.twin(:,:,twin_set)));
                end
                if inl == REC.radiometry.nL
                    %clear DataTD sdTD
                    tilefigs;
                end
                end
            end
        end
    end
    clear inl shift meas_set
    %clear DataTD sdTD
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
    res_fact = [1 1 1];
    % ------------------------ Reference mua ----------------------------------
    if SHOWPLOTS
        if ~SPECTRA
              Nr = ceil(sqrt(REC.radiometry.nL));Nc = round(sqrt(REC.radiometry.nL));
              nfig = 300;
              fh = figure(nfig);
              for inl = 1:REC.radiometry.nL
                  SubPlotMap(REC.opt.Mua(:,:,:,inl),...
                    [num2str(REC.radiometry.lambda(inl)) ' nm'],nfig,Nr,Nc,inl,res_fact);
              end
              fh.NumberTitle = 'off';fh.Name = 'Simulated Absorption';
        % ------------------------ Reference musp ----------------------------------          
             nfig = 400; 
             fh = figure(nfig);
             for inl = 1:REC.radiometry.nL
                  SubPlotMap(REC.opt.Musp(:,:,:,inl),...
                    [num2str(REC.radiometry.lambda(inl)) ' nm'],nfig,Nr,Nc,inl,res_fact);
             end
             fh.NumberTitle = 'off';fh.Name = 'Simulated Scattering'; 
        end
        
        if SPECTRA
            %% plot conc
            Nr = 3; Nc = 3;
            nfig = 900;
            fh = figure(nfig);
            for ic = 1:REC.spe.nCromo
                SubPlotMap(REC.opt.Conc(:,:,:,ic),...
                    [REC.spe.cromo_label{ic} ' Map'],nfig,Nr,Nc,ic,res_fact);
            end
            % Hbtot and SO2
            REC.opt.HbTot = REC.opt.Conc(:,:,:,strcmpi(REC.spe.cromo_label,'hb'))+...
                    REC.opt.Conc(:,:,:,strcmpi(REC.spe.cromo_label,'hbo2'));
            REC.opt.So2 = REC.opt.Conc(:,:,:,strcmpi(REC.spe.cromo_label,'hbo2'))./REC.opt.HbTot;
            SubPlotMap(REC.opt.HbTot,'HbTot Map',nfig,Nr,Nc,ic+1,res_fact);
            SubPlotMap(REC.opt.So2,'So2 Map',nfig,Nr,Nc,ic+2,res_fact);

            % a b scattering
            SubPlotMap(REC.opt.A,'a Map',nfig,Nr,Nc,ic+3,res_fact);
            SubPlotMap(REC.opt.B,'b Map',nfig,Nr,Nc,ic+4,res_fact);

            fh.NumberTitle = 'off';fh.Name = 'Simulated';
                %%
    %         for ic = 1:REC.spe.nCromo
    %             PlotConc = squeeze(REC.opt.Conc(:,:,:,ic));
    %             fh=figure(900+ic);fh.NumberTitle = 'off';fh.Name = [REC.spe.cromo_label{ic} ' Map'];
    %             ShowRecResults(REC.grid,PlotConc,...
    %                 REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto',0,max(PlotConc(:)));%,0.,0.64);
    %             suptitle(REC.spe.cromo_label{ic});
    %         end

            %% detailed visualization
    %         REC.opt.HbTot = REC.opt.Conc(:,:,:,strcmpi(REC.spe.cromo_label,'hb'))+...
    %                 REC.opt.Conc(:,:,:,strcmpi(REC.spe.cromo_label,'hbo2'));
    %         REC.opt.So2 = REC.opt.Conc(:,:,:,strcmpi(REC.spe.cromo_label,'hbo2'))./REC.opt.HbTot;
    %         fh=figure(900+ic+1);fh.NumberTitle = 'off';fh.Name = 'HbTot Map';
    %         ShowRecResults(REC.grid,REC.opt.HbTot,...
    %             REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
    %         suptitle('HbTot');
    %         fh=figure(900+ic+2);fh.NumberTitle = 'off';fh.Name = 'So2 Map';
    %         ShowRecResults(REC.grid,REC.opt.So2,...
    %             REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
    %         suptitle('So2');
    %         fh=figure(900+ic+3);fh.NumberTitle = 'off';fh.Name = ('a Map');
    %         ShowRecResults(REC.grid,REC.opt.A,...
    %             REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
    %         suptitle('a');
    %         fh=figure(900+ic+4);fh.NumberTitle = 'off';fh.Name = ('b Map');
    %         ShowRecResults(REC.grid,REC.opt.B,...
    %             REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
    %         suptitle('b');
        end
        drawnow
    end
%% db
% for i= 1:size(REC.Data,2)
% a = findSQcouple(REC.dmask,i, 'sub');
% d = REC.Data(:,:);
% r = REC.ref(:,:);
% irf = REC.time.irf.data_mask(:,:);
% t = linspace(REC.time.roi(1,i),REC.time.roi(2,i),NUM_TW);
% figure(1); semilogy(t,d(:,i)),hold on,
% %semilogy(irf(:,i)), 
% semilogy(t,r(:,i)),title(sprintf( '(d%d,s%d,l%d) %g', a(1),a(2),a(3),dist(a(1),a(2),a(3))))
% pause
% % if mod(i,8)==0
% %     hold off
% % end
% end
% 
%cleclose all;for i = 1:size(DataTD,2);figure(2000);dsl=findSQcouple(REC.dmask, i, 'sub'),dsl, semilogy(DataTD(:,i),'r'),hold on;semilogy(RefTD(:,i),'k'),hold off;pause();end    
    
    %% fit reference to get starting values for Jacobian based solvers
    if isfield(REC.solver,'fit_reference') 
        if REC.solver.fit_reference == true
            
            if REC.solver.fit_reference_far == true
                fit_dmask = logical(zeros(size(REC.dmask)));                         
                dist = sqrt((REC.Detector.Pos(:,1) - REC.Source.Pos(:,1)' ).^2 + ...
                     (REC.Detector.Pos(:,2) - REC.Source.Pos(:,2)' ).^2 + ...
                  (REC.Detector.Pos(:,3) - REC.Source.Pos(:,3)' ).^2);
                dist = repmat(dist,[1,1,REC.radiometry.nL]);
                Imax = find(dist == max(dist(:)));
                fit_dmask(Imax) = true; 
                fit_dmask = fit_dmask .* REC.dmask;
                fit_dmaskMAT = fit_dmask;
                fit_dmask(~REC.dmask(:)) = [];
            else
                fit_dmask = REC.dmask;
                fit_dmaskMAT = fit_dmask;
            end
            if ~isfield(REC.solver,'fit_type') || strcmpi(REC.solver.fit_type,'spectral')==0
                idxmeas = findMeasIndex(REC.dmask);
                for inl = 1:REC.radiometry.nL
                    meas_set = idxmeas{inl};
                    fitref = REC.ref(:,meas_set);
                    fitsd = REC.sd(:,meas_set);
                    fittwin = REC.time.twin(:,:,meas_set);
                    fitirf = REC.time.irf.data_mask(:,meas_set);
                    if REC.solver.fit_reference_far == true
                        % delete unselected
                        fitref(:,~fit_dmask(meas_set)) = [];
                        fitsd(:,~fit_dmask(meas_set)) = [];
                        fitirf(:,~fit_dmask(meas_set)) = [];
                        fittwin(:,:,~fit_dmask(meas_set)) = [];
                    end
                    if REC.solver.fit_reference_tw == 1 && REC.solver.fit_reference_far == 1
                        fitref = RefTD(:,meas_set);
                        fitsd = sdTD(:,meas_set);
                        fitref(:,~fit_dmask(meas_set)) = [];
                        fitsd(:,~fit_dmask(meas_set)) = [];
                        fittwin = repmat(CreateTimeWindows(REC.time.nstep,[1,REC.time.nstep],'even',REC.time.nstep), [1,1,size(fittwin,3)]); 
                    end
                    fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                       [tmp_mua0(:,inl),tmp_musp0(:,inl)] = FitMuaMus_TD(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0(inl),REC.opt.musp0(inl),REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,fit_dmaskMAT(:,:,inl),...
                        REC.time.dt,REC.time.nstep,fittwin,...
                        REC.time.self_norm,fitref,...
                        fitirf,fitref,fitsd,REC.solver.fit_reference_fwd);
                    if SHOWPLOTS
                        dh=gcf;dh.Name = ['Wavelength ',num2str(REC.radiometry.lambda(inl))];
                        fh(inl)=copyobj(dh,0); delete(dh);
                    end
                end
                [~,~,REC.opt.conc0,REC.opt.a0,REC.opt.b0]=FitVoxel(tmp_mua0(1,:),tmp_musp0(1,:),REC.spe);
            elseif  strcmpi(REC.solver.fit_type, 'spectral')
                    fitref = REC.ref;
                    fitsd = REC.sd;
                    fittwin = REC.time.twin;
                    fitirf = REC.time.irf.data_mask;
                
                if REC.solver.fit_reference_far == true
                    % delete unselected
                    fitref(:,~fit_dmask(:)) = [];
                    fitsd(:,~fit_dmask(:)) = [];
                    fitirf(:,~fit_dmask(:)) = [];
                    fittwin(:,:,~fit_dmask(:)) = [];
                end
                [tmp_mua0,tmp_musp0,tmp_conc,tmp_a,tmp_b] = SpectralFitMuaMus_TD(REC.solver,...
                REC.grid,...
                REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                REC.Source.Pos,REC.Detector.Pos,fit_dmask,...
                REC.time.dt,REC.time.nstep,REC.time.twin,...
                REC.time.self_norm,fitref,...
                fitirf,fitref,fitsd,1,REC.radiometry,REC.spe);
                REC.opt.conc0 = tmp_conc(1,:)';
                REC.opt.a0 = tmp_a(1,:);
                REC.opt.b0 = tmp_b(1,:);
            end
            REC.opt.mua0 = tmp_mua0(1,:);
            REC.opt.musp0 = tmp_musp0(1,:);
            
            REC.opt.muaB = REC.opt.mua0;
            REC.opt.muspB = REC.opt.musp0;
            
            REC.opt.concB = REC.opt.conc0;
            REC.opt.aB = REC.opt.a0;
            REC.opt.bB = REC.opt.b0;
            REC.spe.opt.conc0=REC.opt.conc0;
            REC.spe.opt.a0=REC.opt.a0;
            REC.spe.opt.b0=REC.opt.b0;
            clear tmp_mua0
            clear tmp_musp0
%             if 0==1
%                 [REC.opt.bmua,REC.opt.bmusp] = SpectralFitMuaMus_TD(REC.solver,...
%                 REC.grid,...
%                 REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
%                 REC.Source.Pos,REC.Detector.Pos,REC.dmask,...
%                 REC.time.dt,REC.time.nstep,REC.time.twin,...
%                 REC.time.self_norm,REC.Data,...
%                 REC.time.irf.data,REC.ref,REC.sd,1,REC.radiometry,REC.spe);
%             end
            
            
        end
    end
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
                        REC.grid,mesh.hMesh,REC.grid.hBasis,REC.opt.mua0,REC.opt.musp0,...
                        REC.opt.nB,mesh.qvec,mesh.mvec,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm, REC.Data,REC.time.irf.data,REC.ref,REC.sd,0);
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
                        REC.grid,mesh.hMesh,REC.grid.hBasis,REC.opt.mua0,REC.opt.musp0,...
                        REC.opt.nB,mesh.qvec,mesh.mvec,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm, REC.Data,REC.time.irf.data,REC.ref,REC.sd,0);
                case {'born','tk0'}
                    
                    REC.solver.prior.refimage = [];
                    original_path = REC.solver.prejacobian.path;
                    reg_par = zeros(REC.radiometry.nL,1);
                    idxmeas = findMeasIndex(REC.dmask);
                    for inl = 1:REC.radiometry.nL
                        meas_set = idxmeas{inl};
                        if REC.solver.prejacobian.load
                            REC.solver.prejacobian.path=strcat(original_path,num2str(REC.radiometry.lambda(inl)),'.mat');
                        end
                        fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),' -------</strong>\n'])
                        [REC.opt.bmua(:,inl),REC.opt.bmusp(:,inl),reg_par(inl,1)] = RecSolverTK0_TD(...
                            REC.solver,REC.grid,REC.opt.mua0(inl),REC.opt.musp0(inl),REC.opt.nB,REC.A,REC.Source.Pos,REC.Detector.Pos,REC.dmask(:,:, inl),REC.time.dt,REC.time.nstep,...
                            REC.time.twin(:,:, meas_set),REC.time.self_norm,REC.Data(:,meas_set),...
                            REC.time.irf.data_mask(:, meas_set),REC.ref(:,meas_set),REC.sd(:,meas_set),REC.type_fwd, REC.radiometry.lambda(inl));
                        if REC.solver.prejacobian.load==0
                            movefile([REC.solver.prejacobian.path,'.mat'],strcat(REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'));
                        end
                        drawnow
                    end
                    REC.solver.prejacobian.path = original_path; clear original_path
                case 'born_spectral_post_proc'
                    REC.solver.prior.refimage = [];
                    original_path = REC.solver.prejacobian.path;
                    idxmeas = findMeasIndex(REC.dmask);
                    for inl = 1:REC.radiometry.nL
                        meas_set = idxmeas{inl}; 
                        if REC.solver.prejacobian.load
                            REC.solver.prejacobian.path=strcat(original_path,num2str(REC.radiometry.lambda(inl)),'.mat');
                        end
                        fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),' -------</strong>\n'])
                        [REC.opt.bmua(:,inl),REC.opt.bmusp(:,inl)] = RecSolverTK0_TD(REC.solver,...
                            REC.grid,...
                            REC.opt.mua0(inl),REC.opt.musp0(inl),REC.opt.nB,REC.A,...
                            REC.Source.Pos,REC.Detector.Pos,REC.dmask(:,:,inl),REC.time.dt,REC.time.nstep,...
                            REC.time.twin(:,:,meas_set),REC.time.self_norm,REC.Data(:,meas_set),...
                            REC.time.irf.data_mask(:,meas_set),REC.ref(:,meas_set),REC.sd(:,meas_set),REC.type_fwd);
                        if REC.solver.prejacobian.load==0
                            movefile([REC.solver.prejacobian.path,'.mat'],strcat(REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'));
                        end
                        drawnow
                    end
                    REC.solver.prejacobian.path = original_path; clear original_path
                    [bmua,bmusp,REC.opt.bConc,REC.opt.bA,REC.opt.bbB]=FitVoxel(REC.opt.bmua,REC.opt.bmusp,REC.spe);
                case 'spectral_tk0'
                    REC.solver.prior.refimage = [];
                    [REC.opt.bmua,REC.opt.bmusp,REC.opt.bConc,REC.opt.bA,REC.opt.bbB] = RecSolverTK0_spectra_TD(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm,REC.Data,...
                        REC.time.irf.data,REC.ref,REC.sd,REC.type_fwd,REC.radiometry,REC.spe,REC.opt.conc0,REC.opt.a0,REC.opt.b0);
                case {'spectral_usprior','spectral_tk1'}
                    if ~isempty(REC.solver.prior.path)
                        REC.solver.prior.refimage = ...
                            priormask3D(REC.solver.prior.path,REC.grid);
                    else
                        disp('No prior is provided in RECSettings_DOT. Reference mua will be used');
                        REC.solver.prior.refimage = REC.opt.Mua(:,:,:,1);
                    end
                    if strcmpi(REC.solver.type,'spectral_tk1')
                        REC.solver.prior.refimage = ones(size(REC.opt.Mua(:,:,:,1)));
                    end
                    REC.solver.prior.refimage = 10*double(REC.solver.prior.refimage)+0.1;
                    [REC.opt.bmua,REC.opt.bmusp,REC.opt.bConc,REC.opt.bA,REC.opt.bbB] = RecSolverTK1_spectra_TD(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm,REC.Data,...
                        REC.time.irf.data_mask,REC.ref,REC.sd,REC.type_fwd,REC.radiometry,REC.spe,REC.opt.conc0,REC.opt.a0,REC.opt.b0);
    
                    
                    %% @Simon: US prior inversion
                case {'usprior','tk1'}
                    if ~isempty(REC.solver.prior.path)
                        REC.solver.prior.refimage = ...
                            priormask3D(REC.solver.prior.path,REC.grid);
                    else
                        disp('No prior is provided in RECSettings_DOT. Reference mua will be used');
                        REC.solver.prior.refimage = REC.opt.Mua(:,:,:,1);
                    end
                    REC.solver.prior.refimage =  double(REC.solver.prior.refimage)*10 + 0.1;
                    if strcmpi(REC.solver.type,'tk1')
                        REC.solver.prior.refimage = ones(size(REC.opt.Mua(:,:,:,1)));
                    end
                    original_path = REC.solver.prejacobian.path;
                    idxmeas = findMeasIndex(REC.dmask);
                    for inl = 1:REC.radiometry.nL
                        if REC.solver.prejacobian.load
                            REC.solver.prejacobian.path=strcat(original_path,num2str(REC.radiometry.lambda(inl)),'.mat');
                        end
                        meas_set = idxmeas{inl};
                        fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                        [REC.opt.bmua(:,inl),REC.opt.bmusp(:,inl)] = RecSolverTK1_TD(REC.solver,...
                            REC.grid,...
                            REC.opt.mua0(inl),REC.opt.musp0(inl),REC.opt.nB,REC.A,...
                            REC.Source.Pos,REC.Detector.Pos,REC.dmask(:,:,inl),...
                            REC.time.dt,REC.time.nstep,REC.time.twin(:,:,meas_set),...
                            REC.time.self_norm,REC.Data(:,meas_set),...
                            REC.time.irf.data_mask(:,meas_set),REC.ref(:,meas_set),REC.sd(:,meas_set),REC.type_fwd);
                        if REC.solver.prejacobian.load==0
                            movefile([REC.solver.prejacobian.path,'.mat'],strcat(REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'));
                        end
                        drawnow
                    end
                    REC.solver.prejacobian.path = original_path; clear original_path
                    %%
                case {'combotk0tk1'}
                    if ~isempty(REC.solver.prior.path)
                        REC.solver.prior.refimage = ...
                            priormask3D(REC.solver.prior.path,REC.grid);
                    else
                        disp('No prior is provided in RECSettings_DOT. Reference mua will be used');
                        REC.solver.prior.refimage = REC.opt.Mua(:,:,:,1);
                    end
                    REC.solver.prior.refimage =  double(REC.solver.prior.refimage)*10 + 0.1;
                    original_path = REC.solver.prejacobian.path;
                    idxmeas = findMeasIndex(REC.dmask);
                    for inl = 1:REC.radiometry.nL
                        if REC.solver.prejacobian.load
                            REC.solver.prejacobian.path=strcat(original_path,num2str(REC.radiometry.lambda(inl)),'.mat');
                        end
                        meas_set = idxmeas{inl};
                        fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                        [REC.opt.bmua(:,inl),REC.opt.bmusp(:,inl)] = RecSolverTK0plusTK1_TD(REC.solver,...
                            REC.grid,...
                            REC.opt.mua0(inl),REC.opt.musp0(inl),REC.opt.nB,REC.A,...
                            REC.Source.Pos,REC.Detector.Pos,REC.dmask(:,:,inl),...
                            REC.time.dt,REC.time.nstep,REC.time.twin(:,:,meas_set),...
                            REC.time.self_norm,REC.Data(:,meas_set),...
                            REC.time.irf.data_mask(:,meas_set),REC.ref(:,meas_set),REC.sd(:,meas_set),REC.type_fwd);
                        if REC.solver.prejacobian.load==0
                            movefile([REC.solver.prejacobian.path,'.mat'],strcat(REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'));
                        end
                        drawnow
                    end
                    REC.solver.prejacobian.path = original_path; clear original_path

                case 'fit'
                    idxmeas = findMeasIndex(REC.dmask);
                    for inl = 1:REC.radiometry.nL
                        meas_set = idxmeas{inl};
                        fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                        [REC.opt.bmua(:,inl),REC.opt.bmusp(:,inl)] = FitMuaMus_TD(REC.solver,...
                            REC.grid,...
                            REC.opt.mua0(inl),REC.opt.musp0(inl),REC.opt.nB,REC.A,...
                            REC.Source.Pos,REC.Detector.Pos,REC.dmask(:,:,inl),...
                            REC.time.dt,REC.time.nstep,REC.time.twin(:,:,meas_set),...
                            REC.time.self_norm,REC.Data(:,meas_set),...
                            REC.time.irf.data_mask(:,meas_set),REC.ref(:,meas_set),REC.sd(:,meas_set),REC.type_fwd);
                        dh=gcf;dh.Name = ['Wavelength ',num2str(REC.radiometry.lambda(inl))];
                        fh(inl)=copyobj(dh,0); delete(dh);
                    end
                case 'spectral_fit'
                    [REC.opt.bmua,REC.opt.bmusp, REC.opt.bConc,REC.opt.bA,REC.opt.bbB] = SpectralFitMuaMus_TD(REC.solver,...
                        REC.grid,...
                        REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,...
                        REC.time.dt,REC.time.nstep,REC.time.twin,...
                        REC.time.self_norm,REC.Data,...
                        REC.time.irf.data_mask,REC.ref,REC.sd,1,REC.radiometry,REC.spe);

                case 'components_fit'
                    if ~isempty(REC.solver.prior.path)
                        REC.solver.prior.refimage = ...
                            priormask3D(REC.solver.prior.path,REC.grid);
                    else
                        disp('No prior is provided in RECSettings_DOT. Reference mua will be used');
                        REC.solver.prior.refimage = REC.opt.Mua;
                    end
                     [REC.opt.bmua,REC.opt.bmusp, REC.opt.bConc,REC.opt.bA,REC.opt.bbB] = FitConcAB_2reg_TD( ... RecSolverTK0_spectra_TD(...
                        REC.solver,...
                        REC.grid,...
                        REC.opt.conc0,REC.opt.b0, REC.opt.a0, REC.opt.nB,REC.A,...
                        REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                        REC.time.twin,REC.time.self_norm,REC.Data,...
                        REC.time.irf.data_mask,REC.ref,REC.sd,REC.type_fwd,REC.radiometry,REC.spe, REC.opt.hete1.geometry);
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
                        REC.time.irf.data,REC.ref,REC.sd,REC.type_fwd);
                    
                case 'fit4param' %% you require a TOAST installation
                    % check if TOAST is correctly installed
                    if ~isempty(REC.solver.prior.path)
                        REC.solver.prior.refimage = ...
                            priormask3D(REC.solver.prior.path,REC.grid, REC.solver.type);
                    else
                        disp('No prior is provided in RECSettings_DOT. Reference mua will be used');
                        REC.solver.prior.refimage = REC.opt.Mua;
                        
                    end
                    idxmeas = findMeasIndex(REC.dmask);
                    for inl = 1:REC.radiometry.nL
                        meas_set = idxmeas{inl};
                        fprintf(['<strong>------- Wavelength ',num2str(REC.radiometry.lambda(inl)),'-------</strong>\n'])
                        [REC.opt.bmua(:,inl),REC.opt.bmusp(:,inl), REC.opt.fitOUTPUT] = Fit2Mua2Mus_TD(REC.solver,...
                            REC.grid,...
                            REC.opt.mua0(inl),REC.opt.musp0(inl),REC.opt.nB,REC.A,...
                            REC.Source.Pos,REC.Detector.Pos,REC.dmask(:,:,inl),...
                            REC.time.dt,REC.time.nstep,REC.time.twin(:,:,meas_set),...
                            REC.time.self_norm,REC.Data(:,meas_set),...
                            REC.time.irf.data_mask(:,meas_set),REC.ref(:,meas_set),REC.sd(:,meas_set),1);
                        dh=gcf;dh.Name = ['Wavelength ',num2str(REC.radiometry.lambda(inl))];
                        fh(inl)=copyobj(dh,0); delete(dh);
                    end
                    
            end
    end
    ROI = REC.time.roi; %#ok<NASGU>
    NW = NUM_TW; %#ok<NASGU>
    if strcmpi(REC.solver.type,'born')||...
            strcmpi(REC.solver.type,'born_spectral_post_proc')||...
            strcmpi(REC.solver.type,'usprior')
        if ~isdir(REC.solver.prejacobian.path(1:end-2))
            mkdir(REC.solver.prejacobian.path(1:end-2))
        end
        for inl = 1:REC.radiometry.nL
            if ~exist([REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'],'file')                
                save([REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'],'ROI','NW');
            else
                save([REC.solver.prejacobian.path,num2str(REC.radiometry.lambda(inl)),'.mat'],'ROI','NW','-append');
            end
        end
    else
        if ~exist([REC.solver.prejacobian.path '.mat'],'file')
            if ~exist(REC.solver.prejacobian.path,'dir'),mkdir(REC.solver.prejacobian.path);end
            save(REC.solver.prejacobian.path,'ROI','NW');
        else
            save(REC.solver.prejacobian.path,'ROI','NW','-append');
        end
    end
    clear ROI NW
    
    if REMOVE_VOXELS
        [REC.opt.bmua, REC.opt.bmusp, REC.opt.bConc] = remove_voxels(REC.opt.bmua, REC.opt.bmusp, REC.opt.bConc,...
            REC.grid.dim, REC.opt.mua0, REC.opt.musp0, REC.opt.conc0);
    end
    %% =============================== DISPLAY RESULTS ===========================================
    % ---------------------------- display mua --------------------------------
    if SHOWPLOTS
        if ~contains(REC.solver.type,'fit')
            if ~contains(REC.solver.type,'spectral')
                Nr = ceil(sqrt(REC.radiometry.nL));Nc = round(sqrt(REC.radiometry.nL));
                nfig = 500;
                fh = figure(nfig);
                for inl = 1:REC.radiometry.nL
                    SubPlotMap(reshape(REC.opt.bmua(:,inl),REC.grid.dim),...
                        [num2str(REC.radiometry.lambda(inl)) ' nm'],nfig,Nr,Nc,inl,res_fact);
                end
                fh.NumberTitle = 'off';fh.Name = 'Reconstructed Absorption';
                %    ------------------------ Reference musp ----------------------------------
                nfig = 600;
                fh = figure(nfig);
                for inl = 1:REC.radiometry.nL
                    SubPlotMap(reshape(REC.opt.bmusp(:,inl),REC.grid.dim),...
                        [num2str(REC.radiometry.lambda(inl)) ' nm'],nfig,Nr,Nc,inl,res_fact);
                end
                fh.NumberTitle = 'off';fh.Name = 'Reconstructed Scattering';

                disp('recon: finished');
            elseif contains(REC.solver.type,'spectral')
                %% plot conc
                Nr = 3; Nc = 3;
                nfig = 800;
                fh = figure(nfig);
                for ic = 1:REC.spe.nCromo
                    SubPlotMap(reshape(REC.opt.bConc(:,ic),REC.grid.dim),...
                        [REC.spe.cromo_label{ic} ' Map'],nfig,Nr,Nc,ic,res_fact);
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

                fh.NumberTitle = 'off';fh.Name = 'Reconstructed';

            end
        end
    end
    
    if isfield(REC.grid,'hBasis')
        mesh.hMesh.delete()
        delete(mesh.hMesh);
        mesh.hMesh = [];
        mesh = rmfield(mesh,'hMesh');
        clear mesh
        REC.grid.hBasis.delete();
        REC.grid.hBasis = [];
        REC.grid = rmfield(REC.grid,'hBasis');
    end
    %% Save reconstruction
    %if ~contains(REC.solver.type,'fit')
    disp(['Reconstruction will be stored in: ', rdir,filename,'_', 'REC.mat']);
    save([rdir,filename,'_', 'REC'],'REC');
    %end
    % =========================================================================
    %%                            Quantify DOT
    % =========================================================================

    if ~contains(REC.opt.hete1.geometry,'LOAD_BKG')
        Q = QuantifyDOT(REC,[],0);
    else
        REC.opt.muaBB = REC.opt.muaB; REC.opt.mua0;
        REC.opt.muaB = REC.opt.mua0;
        REC.opt.muspBB = REC.opt.muspB; 
        REC.opt.muspB = REC.opt.musp0;
        Q = QuantifyDOT(REC,[],0);
    end

    Quantification_MultiSim
    % =========================================================================
    %%                            Save quantification
    % =========================================================================
    %if ~contains(REC.solver.type,'fit')
    disp(['Quantification will be stored in: ', rdir,filename,'_', 'REC.mat']);
    save([rdir,filename,'_', 'REC'],'Q','-append')
    %end

end
