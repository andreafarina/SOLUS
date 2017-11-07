%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLUS: Smart OpticaL and UltraSound diagnostics of breast cancer      %
%                                                                       %
% v1:                                                                   %
% Simulation and reconstruction using Green's function.                 %
% Experimental data are not yet implemented.                            %
% Half-space geometry with fwd and Jacobian for the fluence under PCBC. %
% CW and time-domain forward and reconstruction.                        %
% Absorption only reconstruction.                                       %
% Possibility to load an a-priori dictionary matrix (see row 526).      % 
% L1 ISTA reconsrtuction is implemented.                                %
% Regularization parameter criteria: L-curve, gcv, manual.              %
% Solver: Tickhonov, Truncated SVD, pcg, gmres                          %
%                                                                       %
% A. Farina 26/12/2016                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
%clear all;
setPath;
addpath subroutines solvers
addpath(genpath('util'))
t_start = cputime;
spacer = ' ------- ';
% ========================================================================= 
%%                          INITIALIZATION 
% =========================================================================
verbosity = 0;
%toastCatchErrors();                     % redirect toast library errors
%toastSetVerbosity(verbosity);
SHOW_MESH = 0;          % 1 to show fluence and projected pattern on the mesh
filename = 'SOLUS_test';% Filename  prefix 
session = '201612';
exp_path = ['../data/',session,'/'];
res_path = ['../results/',session,'/'];
exp_file = 'SOLUS_test';
%==========================================================================
%%                              OPTIONS 
%==========================================================================
%SET_QM = 1;             % 1: create qvec,mvec. 
                        % 0: load external qvec,mvec from _DOT file.
% ----------------------------- FORWARD -----------------------------------
FORWARD = 1;            % Simulated forward data and save into _Data file
REF = 1;                % 1: create also the homogeneous data
% ------------------------- RECONSTRUCTION --------------------------------
RECONSTRUCTION = 1;     % Enable the reconstruction section.
% ------------------------- EXPERIMENTAL ----------------------------------
EXPERIMENTAL = 1;       % Enable experimental options below
EXP_IRF = 1;            % Use the experimental IRF for forward and reconstruction.
EXP_DELTA = 'all';    % Substitute the IRF with delta function on the 
                        % baricenter ('baric') or peak ('peak') of the IRF.
                        % 'all' to use the experimental IRF.
                    
EXP_DMD = 0;            % Use the experimental registration data for surce and detector
EXP_DATA = 0;           % Load experimental data and use them for reconstruction
% -------------------------------------------------------------------------
DOT.TYPE = 'pointlike';   % 'pointlike','linesources' or 'pattern'
DOT.TD = 1;             % Time-domain: enable the calculation of TPSF
% -------------------------------------------------------------------------
DOT.sigma = 0;%1e-3;%1e-3;       % add gaussian noise to CW data
% -------------------------------------------------------------------------
geom = 'semi-inf';      % geometry
type = 'Born';          % heterogeneous model  
% =========================================================================
%%                     Create folder for saving data
% =========================================================================
rdir = res_path;
%['..',filesep,'results',filesep,datestr(clock, 'yyyymmdd'),filesep];
if ~exist(rdir,'dir')
    mkdir(rdir)
end
%==========================================================================
%%                          Experimental data
% Load the structure EXP containing all the data for analyzing experimental
% data. This structure is created with the routine 'CreateExperiment.m'.
%==========================================================================
if (EXPERIMENTAL == 1)
    load([exp_path,'EXP_',exp_file])
end

%==========================================================================

%% ====================== VOLUME DEFINITION ===============================
%% Background optical properties
DOT.opt.muaB = 0.01;          % mm-1
DOT.opt.muspB = 1;             % mm-1
DOT.opt.nB = 1.4;
DOT.opt.nE = 1.;   % external refractive index
DOT.opt.cs = 0.299/DOT.opt.nB;       % speed of light in medium
DOT.opt.kap = 1/(3*DOT.opt.muspB);
DOT.A = A_factor(DOT.opt.nB/DOT.opt.nE); % A factor for boundary conditions

%==========================================================================
%%                                  SET GRID
%==========================================================================
DOT.grid.x1 = 0;
DOT.grid.x2 = 64;
DOT.grid.dx = 2;

DOT.grid.y1 = 0;
DOT.grid.y2 = 58;           
DOT.grid.dy = DOT.grid.dx;

DOT.grid.z1 = 0;        
DOT.grid.z2 = 32;         
DOT.grid.dz = DOT.grid.dx;

DOT.grid = setGrid(DOT); 

DOT.opt.Mua = ones(DOT.grid.dim) * DOT.opt.muaB;
DOT.opt.Musp = ones(DOT.grid.dim) * DOT.opt.muspB;
%==========================================================================
%%                      Set Heterogeneities
%==========================================================================
%--------------------------- INCLUSION 1 ---------------------------------%
DOT.opt.hete.type  = 'Mua';
DOT.opt.hete.geometry = 'Sphere';
DOT.opt.hete.c     = [35, 25, 10];   % down
% DOT.opt.hete.d     = (M * [0, 0, -1]')';   % down
% DOT.opt.hete.l     = 20;
DOT.opt.hete.sigma = 5;
DOT.opt.hete.distrib = 'OFF';
DOT.opt.hete.profile = 'Step';%'Gaussian';
DOT.opt.hete.val   = 2. * DOT.opt.muaB;
[DOT,DOT.opt.hete] = setHete(DOT,DOT.opt.hete);
%--------------------------- INCLUSION 2 ---------------------------------%
% DOT.opt.hete2.type  = 'Mua';
% DOT.opt.hete2.geometry = 'Sphere';
% DOT.opt.hete2.c     = [40, 20, 5];   % down
% % DOT.opt.hete.d     = (M * [0, 0, -1]')';   % down
% % DOT.opt.hete.l     = 20;
% DOT.opt.hete2.sigma = 3;
% DOT.opt.hete2.distrib = 'OFF';
% DOT.opt.hete2.profile = 'Gaussian';%'Gaussian';
% DOT.opt.hete2.val   = 1.5 * DOT.opt.muaB;
% [DOT,DOT.opt.hete2] = setHete(DOT,DOT.opt.hete2);
%==========================================================================
%%                         Time domain parameters
%==========================================================================
if DOT.TD == 1
    DOT.time.dt = (50e3/4096/6) * 4;    % time step in picoseconds
    DOT.time.nstep = 260;               % number of temporal steps
    DOT.time.noise = 'Poisson';         % 'Poisson','Gaussian','none'
                                        % if 'Poisson' and sigma>0 a
                                        % Gaussian noise is added before
                                        % Poisson noise.
    DOT.time.sigma = 1e-2;              % variance for gaussian noise
    DOT.time.self_norm = false;         % true for self-normalized TPSF   
    DOT.time.TotCounts = 1e6;           % total counts for the maximum-energy
                                        % TPSF. The other are consequently
                                        % rescaled

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
else
    DOT.time.irf.data = 0;
end
end

%==========================================================================
%%  SETTING SOURCES (QVEC), DETECTORS (MVEC) AND THEIR PERMUTATIONS (DMASK)
%==========================================================================
% SOLUS SOURCES - DETECTOR POSITIONS
xs = linspace(5,55,4);
ys = [10,30];
%ys = linspace(5,50,8);
zs = 0;

[xxs,yys,zzs] = ndgrid(xs,ys,zs);

DOT.Source.Pos = [xxs(:),yys(:),zzs(:)];
DOT.Detector.Pos = [xxs(:),yys(:),zzs(:)];
DOT.Source.Ns=size(DOT.Source.Pos,1);
DOT.Detector.Nd=size(DOT.Detector.Pos,1);

%% Define permutation matrix
% ALL COMBINATIONS: null-distances + all the other combinations
%DOT.dmask = logical(ones(DOT.Detector.Nd,DOT.Source.Ns));

% NULL-DISTANCE ONLY
%DOT.dmask = logical(eye(DOT.Detector.Nd,DOT.Source.Ns));

% ALL EXCEPT NULL-DISTANCE
DOT.dmask = logical(ones(DOT.Detector.Nd,DOT.Source.Ns) - ...
    diag(diag(ones(DOT.Detector.Nd,DOT.Source.Ns))));              
% -------------------------------------------------------------------------
% plot DMASK
figure, imagesc(DOT.dmask),xlabel('Sources'),ylabel('Meas'),title('Dmask');

%% plot source-detectors and sphere
figure,plot3(DOT.Source.Pos(:,1),DOT.Source.Pos(:,2),DOT.Source.Pos(:,3),'r*'),grid,
xlabel('x'),ylabel('y'),zlabel('z'),hold on
plot3(DOT.opt.hete.c(1),DOT.opt.hete.c(2),DOT.opt.hete.c(3),'bo'),
[a,b,c]= sphere(100);
surf(a*DOT.opt.hete.sigma + DOT.opt.hete.c(1), ...
    b*DOT.opt.hete.sigma + DOT.opt.hete.c(2),...
    c*DOT.opt.hete.sigma + DOT.opt.hete.c(3))
set(gca,'zdir','reverse'),axis equal,
xlim([DOT.grid.x1 DOT.grid.x2]),...
    ylim([DOT.grid.y1 DOT.grid.y2]),...
    zlim([DOT.grid.z1 DOT.grid.z2])

%==========================================================================
% The structure DOT contains all geometrical parameters needed also for 
% the inverse problem
%==========================================================================
%%                             FORWARD PROBLEM
%==========================================================================
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
        z = zeros(size(DataTD,1) + size(DOT.time.irf.data,1)-1,nmeas);
        for i = 1:nmeas
            z(:,i) = conv(full(DataTD(:,i)),DOT.time.irf.data);
        end
        figure(101);semilogy(z(:,1),'b'),hold,semilogy(DataTD(:,1),'k'),
        semilogy(DOT.time.irf.data,'r'),ylim([max(DataTD(:))/10000 max(DataTD(:))])
        DataTD = z(1:numel(DOT.time.irf.data),:);
                
        if REF == 1
            z = zeros(size(z));
            for i = 1:nmeas
                z(:,i) = conv(full(RefTD(:,i)),DOT.time.irf.data);
            end
            RefTD = z(1:numel(DOT.time.irf.data),:);
        end
    end
    clear z

    %-------------------- Add noise to TD data ------------------------
    sdTD = ones(size(DataTD));
    if ~strcmpi(DOT.time.noise,'none')
        [DataTD,sdTD] = AddNoise(DataTD,'gaussian',DOT.time.sigma);
        if REF == 1
            RefTD = AddNoise(RefTD,'gaussian',DOT.time.sigma);
        end
   end
    if strcmpi(DOT.time.noise,'poisson')
        if (REF == 1) %&& (DOT.time.self_norm == 0))
            [m,j] = max(RefCW);
            factor = 1./sum(RefTD(:,j)) .* DOT.time.TotCounts;
            RefTD = RefTD * factor;
            RefTD = poissrnd(round(full(RefTD)));
            DataTD = DataTD * factor;
            DataTD = poissrnd(round(full(DataTD)));
        else
            [m,j] = max(DataCW);
            factor = 1./sum(DataTD(:,j)) .* DOT.time.TotCounts;
            DataTD = DataTD * factor;
            DataTD = poissrnd(round(full(DataTD)));
        end
    end
    sdTD = sqrt(DataTD);        
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
% =========================================================================
%%                            Save forward data
% =========================================================================
% rdir = ['..',filesep,'Results',filesep];
% if ~exist(rdir,'dir')
%     mkdir(rdir) 
% end
% disp(['Results will be stored in: ', rdir,filename,'_'])
% %-------------------------- Save simulated data ---------------------------
% if FORWARD == 1
%     save([rdir,filename,'_', 'Data'],'DataCW','sdCW');
%     if REF == 1
%         save([rdir,filename,'_', 'Data'],'RefCW','-append');
%     end
%     if DOT.TD == 1
%         save([rdir,filename,'_', 'Data'],'DataTD','sdTD','-append');
%         if REF == 1
%             save([rdir,filename,'_', 'Data'],'RefTD','-append');
%         end
%     end
%     
% end

%==========================================================================
t_end_fwd = cputime - t_start;
disp(['CPU time FWD: ' num2str(t_end_fwd)])
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
% AF: to save memory allocate DataTD and RefTD directly on lines 782/783
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
REC = DOT;clear DOT;
%==========================================================================
%%                      RECONSTRUCTION DOMAIN: CW or TD
%==========================================================================
REC.domain = 'td';          % CW or TD: data type to be inverted
%==========================================================================
%%                               Load data
%==========================================================================
if (EXP_DATA == 0)
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
      
  twin = CreateTimeWindows(REC.time.nstep,[10,REC.time.nstep],'even',20);
  REC.time.twin = twin + 90;
  REC.time.nwin = size(REC.time.twin,1);

  % plot roi on the first measruement
  figure(1000);
  semilogy(DataTD),ylim([max(DataTD(:))/10000 max(DataTD(:))])
  for i = 1:REC.time.nwin
    rectangle('Position',[double(REC.time.twin(i,1)),min(DataTD(DataTD(:)>0)),...
        double(REC.time.twin(i,2)-REC.time.twin(i,1)+1),double(max(DataTD(:,1)))]);
  end
  if exist('RefTD','var')
      REC.ref =WindowTPSF(RefTD,REC.time.twin);
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

% =========================================================================
%%                        Initial parameter estimates 
% =========================================================================
% In this section all the parameter for the inverse solver are setted.
% --------------------------- Optical properties --------------------------
REC.opt.mua0 = 0.01;    % absorption [mm-1]
REC.opt.musp0 = 1;    % reduced scattering [mm-1]
REC.opt.nB = 1.4;
REC.cm = 0.3/REC.opt.nB;
%REC.freq = 0;
% ---------------------- Solver and regularization ------------------------
REC.solver.tau = 1e-1;            % regularisation parameter
REC.solver.type = 'USprior';      % 'born','GN': gauss-newton, 
                                  % 'USprior': Simon's strutural prior
                                  % 'LM': Levenberg-Marquardt,
                                  % 'l1': L1-based minimization
                                  % 'fit': fitting homogeneous data
REC.solver.prior = REC.opt.Mua;   % dictionnary matrix, or []
% ---------------------- Set nodes optical properties --------------------- 
REC.opt.mua = ones(REC.grid.N,1) * REC.opt.mua0;
REC.opt.musp = ones(REC.grid.N,1) * REC.opt.musp0;
REC.opt.n = ones(REC.grid.N,1) * REC.opt.nB;
REC.opt.kap = 1./(3*(REC.opt.musp));
%==========================================================================
%%              PLOT REFERENCE VALUES IF PRESENT              
%==========================================================================        
% ------------------------ Reference mua ----------------------------------
if isfield(REC.opt,'Mua')
    figure(301);
    ShowRecResults(REC.grid,REC.opt.Mua,...
        REC.grid.z1,REC.grid.z2,REC.grid.dz,1,...
          0,max(REC.opt.Mua(:)));
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
%                 [REC.opt.bmua,REC.opt.bmusp] = RecSolverBORN_TD(REC.solver,...
%                     REC.grid,...
%                     REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
%                     REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
%                     REC.time.twin,REC.time.self_norm,REC.Data,...
%                     REC.time.irf.data,REC.ref,REC.sd,0);
%% @Simon: US prior inversion
            case 'usprior'
                [REC.opt.bmua,REC.opt.bmusp] = RecSolverBORN_TD_USPrior(REC.solver,...
                    REC.grid,...
                    REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                    REC.Source.Pos,REC.Detector.Pos,REC.dmask,...
                    REC.time.dt,REC.time.nstep,REC.time.twin,...
                    REC.time.self_norm,REC.Data,...
                    REC.time.irf.data,REC.ref,REC.sd,0);
%%                
            case 'fit'
%                 [REC.opt.bmua,REC.opt.bmusp] = FitMuaMus_TD(REC.solver,REC.mesh.hMesh,...
%                     REC.grid.hBasis,REC.opt.notroi,...
%                     REC.opt.mua0,REC.opt.musp0,REC.opt.nB,...
%                     REC.Q,REC.M,REC.dmask,REC.time.dt,REC.time.nstep,...
%                     REC.time.twin,REC.time.self_norm, REC.ref,REC.sd,...
%                     REC.time.irf.data,1);
                
            case 'cg'
            
            case 'l1'
%                 [REC.opt.bmua,REC.opt.bmusp] = RecSolverL1_TD(REC.solver,...
%                     REC.grid,...
%                     REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
%                     REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
%                     REC.time.twin,REC.time.self_norm,REC.Data,...
%                     REC.time.irf.data,REC.ref,REC.sd,0);
                
                
        end
end

% ---------------------------- display mua --------------------------------
drawnow;
figure(304);
ShowRecResults(REC.grid,reshape(REC.opt.bmua,REC.grid.dim),...
   REC.grid.z1,REC.grid.z2,REC.grid.dz,1,0.00,0.05);
suptitle('Recon Mua');
% ---------------------------- display musp -------------------------------
% figure(305);
% ShowRecResults(REC.grid,reshape(REC.opt.bmusp,REC.grid.dim),...
%    REC.grid.z1,REC.grid.z2,REC.grid.dz,1);%,0.,0.64);
% suptitle('Recon Mus');
drawnow;
tilefigs;
disp('recon: finished')

end