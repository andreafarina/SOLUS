%==========================================================================
%%                              OPTIONS 
%==========================================================================
% ----------------------------- FORWARD -----------------------------------
FORWARD = 1;            % Simulated forward data and save into _Data file
REF = 1;                % 1: create also the homogeneous data
TYPE_FWD = 'fem';       % fwd model computation: 'linear','fem'
geom = 'semi-inf';      % geometry
% ------------------------- RECONSTRUCTION --------------------------------
RECONSTRUCTION = 1;     % Enable the reconstruction section.
% ------------------------- EXPERIMENTAL ----------------------------------
EXPERIMENTAL = 1;       % Enable experimental options below
EXP_IRF = 1;            % Use the experimental IRF for forward and 
                        % reconstruction.
EXP_DELTA = 'all';      % Substitute the IRF with delta function on the 
                        % baricenter ('baric') or peak ('peak') of the IRF.
                        % 'all' to use the experimental IRF.                    
EXP_DATA = 0;           % Load experimental data and use them for 
                        % reconstruction
% -------------------------------------------------------------------------
%DOT.TYPE = 'pointlike'; % 'pointlike','linesources' or 'pattern'
DOT.TD = 1;             % Time-domain: enable the calculation of TPSF
% -------------------------------------------------------------------------
DOT.sigma = 0;          % add gaussian noise to CW data
% -------------------------------------------------------------------------
geom = 'semi-inf';      % geometry
type = 'Born';          % heterogeneous model  
% -------------------------------------------------------------------------
RADIOMETRY = 1;         % apply radiometric inputs to simulated data
% -------------------------------------------------------------------------
SAVE_FWD = 1;           % Save forward data (possibly with noise) 
                        % in a _Data.m file
LOAD_FWD_TEO = 0;       % if 0: save the raw TPSF(un-noisy) in a _FwdTeo.m file.
                        % if 1: load the raw TPSF for speed up
% -------------------------------------------------------------------------
TOAST2DOT = 0;          % if 1 the function toast2dot is used for conversion 
SPECTRA = 1;
DOT.spe.cromo_label = {'Hb','HbO2','Lipid','H2O','Collagen'};
DOT.spe.active_cromo = [1,1,1,1,1];
DOT.spe.cromo_factor = [1,1,1,1,1];
DOT.spe.cromo_units = {'microM','microM','mg/cm^3','mg/cm^3','mg/cm^3'};
DOT.spe.ForceConstitSolution = 1;
if SPECTRA == 0
mua_ = [0.003807,0.001342,0.001387,0.010060,0.007554,0.003942,0.007613,0.004933];
musp_ = [1.22842,1.189904,0.913511,0.803037,0.769676,0.7560782,0.702287,0.675179];
Xd = {mua_,musp_};
else
a_ = 1.2883;	b_ = 1.2641;
conc_ = [1.91, 8.8, 575, 270, 108];
Xd = {conc_,[a_ b_]};
end
% ========================================================================= 
%% ====================== VOLUME DEFINITION ===============================
%% Background optical properties
% DOT.opt.muaB = 0.01;    % mm-1
% DOT.opt.muspB = 1;      % mm-1
DOT.opt.nB = 1.4;       % internal refractive index   
DOT.opt.nE = 1.;        % external refractive index
%==========================================================================
%%                                  SET GRID
%==========================================================================
DOT.grid.x1 = -32; %32
DOT.grid.x2 = 32;
DOT.grid.dx = 2;

DOT.grid.y1 = -29;%29
DOT.grid.y2 = 29;           
DOT.grid.dy = DOT.grid.dx;

DOT.grid.z1 = 0;    %32    
DOT.grid.z2 = 30;         
DOT.grid.dz = DOT.grid.dx;
%==========================================================================
%%                      Set Heterogeneities
%==========================================================================
NUM_HETE = 1;
%--------------------------- INCLUSION 1 ---------------------------------%
DOT.opt.hete1.type  = {'Mua','Musp'};
DOT.opt.hete1.geometry = 'USimage';
%DOT.opt.hete1.randinhom = [10; 0.001]; %set typical correlation length (mm) and intensity of the random inhomogenoities
DOT.opt.hete1.c     = [0, -5, 15];   % down
% DOT.opt.hete1.d     = [0, 0, -1];   % down
% DOT.opt.hete1.l     = 20;
DOT.opt.hete1.sigma = 5;
DOT.opt.hete1.distrib = 'OFF';
DOT.opt.hete1.profile = 'Gaussian';%'Step';%'Gaussian';
% DOT.opt.hete1.val   = 5 * DOT.opt.muaB;
if SPECTRA == 0
muap_ = [0.037321,0.021343,0.010411,0.015747,0.023087,0.051566,0.029411,0.019642];
muspp_ = [0.371159,0.314731,0.214725,0.208025,0.222794,0.297037,0.220381,0.189757];
Xp = {muap_,muspp_};
else
a_ = 0.5 * 1.2883;	b_ = 2 * 1.2641;    
%a_ = (1+10/100)*1.5453;	b_ = 1; % a (mm-1), b(adimensionale)
concp_ = [0.5 0.5 3 2 0.5].*[1.91 8.8 575 270 108];%*2
Xp = {concp_,[a_ b_]};
end
DOT.opt.hete1.path ='ORIGINALkwaved.mat';%'ham_maskYX_flip.mat' ;   % down

%--------------------------- INCLUSION 2 ---------------------------------%
% DOT.opt.hete2.type  = 'Mua';
% DOT.opt.hete2.geometry = 'Sphere';
% DOT.opt.hete2.c     = [5, 20, 5];   % down
% DOT.opt.hete2.d     = [ 1, 0, 0];   % down
% % DOT.opt.hete.d     = (M * [0, 0, -1]')';   % down
% % DOT.opt.hete.l     = 20;
% DOT.opt.hete2.sigma = 5;
% DOT.opt.hete2.distrib = 'OFF';
% DOT.opt.hete2.profile = 'Gaussian';%'Gaussian';
% DOT.opt.hete2.val   = 2 * DOT.opt.muaB;

%==========================================================================
%%                         Time domain parameters
%==========================================================================
DOT.time.dt = 50;%(50e3/4096/6)*4;        % time step in picoseconds
DOT.time.nstep = 90;%4096/4;               % number of temporal steps
DOT.time.noise = 'poisson';         % 'Poisson','Gaussian','none'
                                    % if 'Poisson' and sigma>0 a
                                    % Gaussian noise is added before
                                    % Poisson noise.
DOT.time.sigma = 1e-3;              % variance for gaussian noise
DOT.time.self_norm = false;         % true for self-normalized TPSF
DOT.time.TotCounts = 1e6;           % total counts for the maximum-energy
                                    % TPSF. The other are consequently
                                    % rescaled
%==========================================================================
%%                         Radiometry
%==========================================================================
DOT.radiometry.power = [1];    % (mW) laser input power %AAA
DOT.radiometry.timebin = ...
    DOT.time.dt;                % (ps) width of the time bin
DOT.radiometry.acqtime = 1;     % (s) acquisition time %AAA
DOT.radiometry.opteff = 0.9;    % Typical efficiency of the optical path
DOT.radiometry.lambda = [635,670,830,915,940,980,1030,1065];
DOT.radiometry.lambda0 = 635;
DOT.radiometry.lambda = DOT.radiometry.lambda;
DOT.radiometry.nL = numel(DOT.radiometry.lambda);
                                % (just for calculation of input photons)
DOT.radiometry.area = 4;        % mm^2
DOT.radiometry.qeff = 0.05;     % Quantum efficiency
if size(DOT.time.TotCounts,2)<DOT.radiometry.nL
    DOT.time.TotCounts = repmat(DOT.time.TotCounts(1),1,DOT.radiometry.nL);
    DOT.radiometry.power = repmat(DOT.radiometry.power(1),1,DOT.radiometry.nL);
    warning('backtrace','off'),warning('verbose','off')
    warning(['Power will be set to ' num2str(DOT.radiometry.power(1)) ' for any wavelenght'])
    warning(['AcqTime will be set to ' num2str(DOT.radiometry.acqtime(1)) ' for any wavelenght'])
    warning('backtrace','on'),warning('verbose','on')
end
%==========================================================================
%%                   Cutting of counts
%==========================================================================
% CUT_COUNTS = 0 --> The count-rate is fixed accordingly to the higher
%                    power measurement. The other are consequently
%                    rescaled. RADIOMETRY in this case has no effect.
% CUT_COUNTS = 1 --> A Time-gated measurement is simulated.
%                    Measurements on each dealay are rescaled accordingly 
%                    to DOT.time.TotCounts. In particular, if:
%                    RADIOMETRY==1, the count-rate for each delay is cut to 
%                         DOT.time.TotCounts if photons are available, 
%                         otherwise no;
%                   RADIOMETRY==0, the count-rate for each delay is cut to 
%                         DOT.time.TotCounts in any case.  
CUT_COUNTS = 1;         
NumDelays = 1;      % number of delays