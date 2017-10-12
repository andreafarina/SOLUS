%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOT_core     %
%                                                                       %
% v2:                                                                   %
% Simulation and reconstruction using Green's function.                 %
% Experimental data are not yet implemented.                            %
% Half-space geometry with fwd and Jacobian for the fluence under PCBC. %
% CW and time-domain forward and reconstruction.                        %
% Absorption only reconstruction.                                       %
% Possibility to load an a-priori dictionary matrix (see row 526).      % 
% L1 ISTA reconstruction is implemented.
% Structural priors regularization is implemented
% Regularization parameter criteria: L-curve, gcv, manual.              %
% Solver: Tickhonov, Truncated SVD, pcg, gmres                          %
%                                                                       %
% A. Farina 22/09/2017                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%clearvars;
close all;

setPath;
addpath subroutines solvers
addpath(genpath('util'))
t_start = cputime;
spacer = ' ------- ';
% ========================================================================= 
%%                          INITIALIZATION 
% =========================================================================
frontpage;
disp('Initializing DOT parameters');
Init_DOT

disp('Setting paths and filenames');
SetIO_DOT

if RECONSTRUCTION == 1
    disp('Setting reconstruction parameters');
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
%==========================================================================
%%                          Experimental data
% Load the structure EXP containing all the data for analyzing experimental
% data. This structure is created with the routine 'CreateExperiment.m'.
%==========================================================================
if (EXPERIMENTAL == 1)
    load([exp_path,'EXP_',exp_file])
end
disp(['Experiment file: ',exp_path,'EXP_',exp_file]);

%% Background optical properties
DOT.opt.cs = 0.299/DOT.opt.nB;       % speed of light in medium
DOT.opt.kap = 1/(3*DOT.opt.muspB);
DOT.A = A_factor(DOT.opt.nB/DOT.opt.nE); % A factor for boundary conditions
%==========================================================================
%%                                  SET GRID
%==========================================================================
DOT.grid = setGrid(DOT); 

DOT.opt.Mua = ones(DOT.grid.dim) * DOT.opt.muaB;
DOT.opt.Musp = ones(DOT.grid.dim) * DOT.opt.muspB;
%==========================================================================
%%                      Set Heterogeneities
%==========================================================================
%--------------------------- INCLUSIONS ---------------------------------%
for i = 1:NUM_HETE
    h_string = ['hete',num2str(i)];
    [DOT,DOT.opt.(h_string)] = setHete(DOT,DOT.opt.(h_string));
end
%==========================================================================
%%                         Time domain parameters
%==========================================================================
if DOT.TD == 1
% --------- resample IRF following the parameters in DOT.time -------------
if ((EXPERIMENTAL == 1) && (EXP_IRF == 1))
    switch lower(EXP_DELTA)
        case 'peak'
            EXP.irf.data = zeros(size(EXP.irf.data));
            EXP.irf.data(EXP.irf.peak.pos) = 1;
        case 'baric'
            EXP.irf.data = zeros(size(EXP.irf.data));
            EXP.irf.data(EXP.irf.baric.pos) = 1;
    end
    DOT.time.irf.data = resampleSPC(EXP.irf.data,EXP.time.axis,DOT.time.dt,'norm');
% ------ setting channel 0 based on the peak and width of the IRF ---------
    %resc_fact = round(DOT.time.dt./EXP.spc.factor);
    [DOT.time.irf.peak.value,...
        DOT.time.irf.peak.pos] = max(DOT.time.irf.data);
    DOT.time.irf.ch0 = round(DOT.time.irf.peak.pos - ...
                    500./DOT.time.dt);

else
    DOT.time.irf.data = 0;
end

end

%==========================================================================
%%  SETTING SOURCES (QVEC), DETECTORS (MVEC) AND THEIR PERMUTATIONS (DMASK)
%==========================================================================
disp('Setting sources(Q) and detectors(M)');
SetQM_DOT

% -------------------------------------------------------------------------
% plot DMASK
figure, imagesc(DOT.dmask),xlabel('Sources'),ylabel('Meas'),title('Dmask');

% plot source-detectors and heterogeneities
PlotHeteQM(DOT,NUM_HETE)
drawnow;
%==========================================================================
% The structure DOT contains all geometrical parameters needed also for 
% the inverse problem
%==========================================================================
%%                             FORWARD PROBLEM
%==========================================================================
%disp('Forwar model computation');
if FORWARD == 1
    nmeas = sum(DOT.dmask(:));
    DataCW = ForwardCW(DOT.grid,DOT.Source.Pos, DOT.Detector.Pos, DOT.dmask, ...
        DOT.opt.muaB, DOT.opt.muspB, DOT.opt.Mua, ...
        DOT.opt.Musp, DOT.A, geom, type);
    [DataCW,sdCW] = AddNoise(DataCW,'gaussian',DOT.sigma);
    if REF == 1
        RefCW  = ForwardCW(DOT.grid,DOT.Source.Pos, DOT.Detector.Pos, DOT.dmask,...
            DOT.opt.muaB, DOT.opt.muspB, DOT.opt.muaB*ones(DOT.grid.dim),...
            DOT.opt.muspB*ones(DOT.grid.dim), DOT.A, geom, 'homo');
        [RefCW,RefsdCW] = AddNoise(RefCW,'gaussian',DOT.sigma);
    end
    
%==========================================================================
%%                             TD Forward 
%==========================================================================
if DOT.TD == 1
    DOT.time.time = (1:DOT.time.nstep) * DOT.time.dt;
    
    DataTD = ForwardTD(DOT.grid,DOT.Source.Pos, DOT.Detector.Pos, DOT.dmask,...
        DOT.opt.muaB, DOT.opt.muspB,DOT.opt.nB, DOT.opt.Mua,...
                DOT.opt.Musp, DOT.A, DOT.time.dt,...
                length(DOT.time.time), DOT.time.self_norm, geom, 'born');
        if REF == 1
            RefTD = ForwardTD(DOT.grid,DOT.Source.Pos, DOT.Detector.Pos, DOT.dmask,...
                DOT.opt.muaB, DOT.opt.muspB,DOT.opt.nB, ...
                DOT.opt.muaB*ones(DOT.grid.dim),DOT.opt.muspB*ones(DOT.grid.dim), ...
                DOT.A, DOT.time.dt,length(DOT.time.time), DOT.time.self_norm,...
                geom, 'homo');
        end
      
    % Convolution with IRF
    % AF: to optimize memory assign directly DataTD
    if ((EXPERIMENTAL == 1) && (EXP_IRF == 1))
        z = convn(full(DataTD),DOT.time.irf.data);
        figure(101);semilogy(z(:,1),'b'),hold,semilogy(DataTD(:,1),'k'),
        semilogy(DOT.time.irf.data,'r'),%ylim([max(DataTD(:))/10000 max(DataTD(:))])
        nmax = max(DOT.time.nstep,numel(DOT.time.irf.data));
        DataTD = z(1:nmax,:);
                
        if REF == 1
            z = convn(full(RefTD),DOT.time.irf.data);
            RefTD = z(1:nmax,:);
        end
        clear nmax
    end
    clear z
    
    % ---- Radiometry -------
    if RADIOMETRY == 1
        RealFactor = Radiometry(DOT.radiometry);
    else
        RealFactor = 1;
    end
    %-------------------- Add noise to TD data ------------------------
    sdTD = ones(size(DataTD)); 
    if ~strcmpi(DOT.time.noise,'none')
        [DataTD,~] = AddNoise(DataTD,'gaussian',DOT.time.sigma);
        if REF == 1
            RefTD = AddNoise(RefTD,'gaussian',DOT.time.sigma);
        end
    end
       
    if strcmpi(DOT.time.noise,'poisson')
        
        if (REF == 1) %&& (DOT.time.self_norm == 0))
            [factor,RefTD,sdTD] = CutCounts(1:size(RefTD,1),DOT.time.dt,...
                RealFactor*RefTD,DOT.time.TotCounts,RADIOMETRY,CUT_COUNTS,...
                            NumDelays);
            switch CUT_COUNTS
                case 0  % 
                    DataTD = poissrnd(round(...
                        DataTD * factor * RealFactor));
                case 1
                    %idz = find(factor>0, 1, 'last' );
                    DataTD = poissrnd(round(DataTD * RealFactor .* factor));
                    %DataTD = bsxfun(@times,DataTD,factor(idz)./factor');
            end
        
            
        else
            [factor,DataTD,sdTD] = CutCounts(1:size(RefTD,1),DOT.time.dt,...
                RealFactor*DataTD,DOT.time.TotCounts,RADIOMETRY,CUT_COUNTS,...
                            NumDelays);
        end
    end
    % self normalized data
    if DOT.time.self_norm == true
        Area = sum(DataTD);
        DataTD = DataTD * spdiags(1./Area',0,nmeas,nmeas);
        sdTD = sqrt(DataTD * spdiags(1./Area',0,nmeas,nmeas));  % Standard deviation
        if REF == 1
            Area = sum(RefTD);
            RefTD = RefTD * spdiags(1./Area',0,nmeas,nmeas);
            sdTD = sqrt(DataTD * spdiags(1./Area',0,nmeas,nmeas));  % Standard deviation
        end
    end
end

  
end
%==========================================================================
t_end_fwd = cputime - t_start;
disp(['CPU time FWD: ' num2str(t_end_fwd)])
% =========================================================================
%%                            Save forward data
% =========================================================================
if SAVE_FWD == 1
    rdir = ['..',filesep,'results',filesep,session,filesep];
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
    nview = numel(EXP.views);
    nqm = size(EXP.data.spc,2)./nview;
     if DOT.TD == 1
  
% % ----- Resample experimental data to match the time-scale of DOT.time ----
        z = zeros(DOT.time.nstep,size(EXP.data.spc,2));
        f = zeros(DOT.time.nstep,size(EXP.data.spc,2));
        
        for i = 1:size(EXP.data.spc,2)
            z(:,i) = resampleSPC(EXP.data.spc(:,i),EXP.time.axis,DOT.time.dt);
            f(:,i) = resampleSPC(EXP.data.ref(:,i),EXP.time.axis,DOT.time.dt);
        end
% ------ Comparison between conv(irf,fwd) and the experimental data -------        
        if FORWARD == 1
            h = 1;
            figure(202);semilogy(1:size(z,1),[f(:,h)./sum(f(:,h)),RefTD(:,h)./sum(RefTD(:,h))]),
            legend('Data','conv(irf,fwd)'),
            title(['Ref measurement number ',num2str(h)]);
        end
        %figure;semilogy(f)
        %figure;semilogy(f./z)
% -------- Save Time-resolved Data,Reference and standard deviation -------
        DataTD = z;
        RefTD = f;
        sdTD = sqrt(DataTD);    % Poisson
        clear z f
    else
% ------------------------------ CW case ----------------------------------        
        DataCW = sum(EXP.data.spc);
        RefCW = sum(EXP.data.ref);
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
  nmax = max(REC.time.nstep,numel(REC.time.irf.data));   
  twin = CreateTimeWindows(REC.time.nstep,[REC.time.irf.ch0 - 1,nmax],'even',NUM_TW);
  REC.time.twin = twin;% + REC.time.irf.ch0 - 1;%+ Chan0-1; % Chan0 is IRF peak channel, add -1 since twin starts from 1
  REC.time.nwin = size(REC.time.twin,1);

  % plot temporal windows and data
  figure(1000);
  semilogy(DataTD),ylim([max(DataTD(:))/10000 max(DataTD(:))])
  for i = 1:REC.time.nwin
    rectangle('Position',[double(REC.time.twin(i,1)),min(DataTD(DataTD(:)>0)),...
        double(REC.time.twin(i,2)-REC.time.twin(i,1)+1),double(max(DataTD(:,1)))]);
  end
  drawnow;
  if exist('RefTD','var')
      REC.ref = WindowTPSF(RefTD,REC.time.twin);
      clear RefTD
  end
  REC.Data = WindowTPSF(DataTD,REC.time.twin);
  if (EXP_DATA == 0)
      if ~strcmpi(REC.time.noise,'none')
          REC.sd = sqrt(WindowTPSF(sdTD.^2,REC.time.twin));     % ATTENTION NOT TO SUM SD
      else
          REC.sd = ones(size(REC.Data));
      end
  end
  if (EXP_DATA == 1)
      REC.sd = sqrt(WindowTPSF(sdTD.^2,REC.time.twin));
  end
clear DataTD sdTD  
tilefigs;
end

% % =========================================================================
% %%                        Initial parameter estimates 
% % =========================================================================
% ---------------------- Set nodes optical properties --------------------- 
REC.opt.mua = ones(REC.grid.N,1) * REC.opt.mua0;
REC.opt.musp = ones(REC.grid.N,1) * REC.opt.musp0;
REC.opt.n = ones(REC.grid.N,1) * REC.opt.nB;
REC.opt.kap = 1./(3*(REC.opt.musp));
REC.cm = 0.299/REC.opt.nB;

%==========================================================================
%%              PLOT REFERENCE VALUES IF PRESENT              
%==========================================================================        
% ------------------------ Reference mua ----------------------------------
if isfield(REC.opt,'Mua')
    figure(301);
    ShowRecResults(REC.grid,REC.opt.Mua,...
        REC.grid.z1,REC.grid.z2,REC.grid.dz,1,...
          'auto',0,max(REC.opt.Mua(:)));
    suptitle('Mua');
end
% ------------------------ Reference musp ---------------------------------
% if isfield(REC.opt,'Musp')
%     figure(302);
%     ShowRecResults(REC.grid,REC.opt.Musp,...
%         REC.grid.z1,REC.grid.z2,REC.grid.dz,1,...
%         0,max(REC.opt.Musp(:)));
%     suptitle('Mus');
% end
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
                
            case 'born'
                if ~isempty(REC.solver.prior)
                    REC.solver.prior = REC.opt.Mua;
                end
                [REC.opt.bmua,REC.opt.bmusp] = RecSolverBORN_TD(REC.solver,...
                    REC.grid,...
                    REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                    REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                    REC.time.twin,REC.time.self_norm,REC.Data,...
                    REC.time.irf.data,REC.ref,REC.sd,0);
%% @Simon: US prior inversion
            case 'usprior'
                REC.solver.prior = REC.opt.Mua;   % dictionnary matrix, or []

                [REC.opt.bmua,REC.opt.bmusp] = RecSolverBORN_TD_USPrior(REC.solver,...
                    REC.grid,...
                    REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                    REC.Source.Pos,REC.Detector.Pos,REC.dmask,...
                    REC.time.dt,REC.time.nstep,REC.time.twin,...
                    REC.time.self_norm,REC.Data,...
                    REC.time.irf.data,REC.ref,REC.sd,0);
%%                
            case 'fit'
                [REC.opt.bmua,REC.opt.bmusp] = FitMuaMus_TD(REC.solver,...
                    REC.grid,...
                    REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                    REC.Source.Pos,REC.Detector.Pos,REC.dmask,...
                    REC.time.dt,REC.time.nstep,REC.time.twin,...
                    REC.time.self_norm,REC.Data,...
                    REC.time.irf.data,REC.ref,REC.sd,1);
                
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
                
                
        end
end

% ---------------------------- display mua --------------------------------
drawnow;
figure(304);
ShowRecResults(REC.grid,reshape(REC.opt.bmua,REC.grid.dim),...
   REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto',0.00,0.05);
suptitle('Recon Mua');
% ---------------------------- display musp -------------------------------
% figure(305);
% ShowRecResults(REC.grid,reshape(REC.opt.bmusp,REC.grid.dim),...
%    REC.grid.z1,REC.grid.z2,REC.grid.dz,1);%,0.,0.64);
% suptitle('Recon Mus');
drawnow;
tilefigs;
disp('recon: finished')
% =========================================================================
%%                            Quantify DOT 
% =========================================================================
if ~strcmpi(REC.solver.type,'fit')
    Q = QuantifyDOT(REC,~EXP_DATA);
end
Quantification_MultiSim
% =========================================================================
%%                            Save reconstruction 
% =========================================================================
if ~strcmpi(REC.solver.type,'fit')
    disp(['Reconstruction will be stored in: ', rdir,filename,'_', 'REC.m']);
    save([rdir,filename,'_', 'REC'],'REC','Q')
end
clearvars REC
end