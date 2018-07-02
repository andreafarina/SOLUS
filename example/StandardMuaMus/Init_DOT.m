%==========================================================================
%%                              OPTIONS 
%==========================================================================
% ----------------------------- FORWARD -----------------------------------
FORWARD = 1;            % Simulated forward data and save into _Data file
REF = 1;                % 1: create also the homogeneous data
TYPE_FWD = 'linear';    % fwd model computation: 'linear','fem'
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
RADIOMETRY = 1;         % apply radiometric inputs to simulated data
% -------------------------------------------------------------------------
SAVE_FWD = 1;           % Save forward data (possibly with noise) 
                        % in a _Data.m file
LOAD_FWD_TEO = 0;       % if 0: save the raw TPSF(un-noisy) in a _FwdTeo.m file.
                        % if 1: load the raw TPSF for speed up
% -------------------------------------------------------------------------
TOAST2DOT = 0;          % if 1 the function toast2dot is used for conversion 
% ========================================================================= 
%% ====================== VOLUME DEFINITION ===============================
%% Background optical properties
DOT.opt.muaB = 0.01;    % mm-1
DOT.opt.muspB = 1;      % mm-1
DOT.opt.nB = 1.4;       % internal refractive index   
DOT.opt.nE = 1.;        % external refractive index
%==========================================================================
%%                                  SET GRID
%==========================================================================
DOT.grid.x1 = -32;
DOT.grid.x2 = 32;
DOT.grid.dx = 4;

DOT.grid.y1 = -29;
DOT.grid.y2 = 29;           
DOT.grid.dy = DOT.grid.dx;

DOT.grid.z1 = 0;        
DOT.grid.z2 = 32;         
DOT.grid.dz = DOT.grid.dx;
%==========================================================================
%%                      Set Heterogeneities
%==========================================================================
NUM_HETE = 1;
%--------------------------- INCLUSION 1 ---------------------------------%
DOT.opt.hete1.type  = {'Mua','Musp'};
DOT.opt.hete1.geometry = 'sphere';
DOT.opt.hete1.c     = [0, 0, 05];   % down
% DOT.opt.hete1.d     = [0, 0, -1];   % down
% DOT.opt.hete1.l     = 20;
DOT.opt.hete1.sigma = 5;
DOT.opt.hete1.distrib = 'OFF';
DOT.opt.hete1.profile = 'Gaussian';%'Step';%'Gaussian';
DOT.opt.hete1.val   = [2 * DOT.opt.muaB,2 * DOT.opt.muspB];%2 * DOT.opt.muspB;
DOT.opt.hete1.path ='../3DMasks/Mask3d_mus1_05-1_56_malignant_3.mat' ;   % down

%--------------------------- INCLUSION 2 ---------------------------------%
% DOT.opt.hete2.type  = 'Mua';
% DOT.opt.hete2.geometry = 'sphere';
% DOT.opt.hete2.c     = [0, 0, 5];   % down
% % DOT.opt.hete1.d     = [0, 0, -1];   % down
% % DOT.opt.hete1.l     = 20;
% DOT.opt.hete2.sigma = 5;
% DOT.opt.hete2.distrib = 'OFF';
% DOT.opt.hete2.profile = 'Gaussian';%'Step';%'Gaussian';
% DOT.opt.hete2.val   = 2 * DOT.opt.muaB;
% DOT.opt.hete2.path ='../3DMasks/Mask3d_mus1_05-1_56_malignant_3.mat' ;   % down


%==========================================================================
%%                         Time domain parameters
%==========================================================================
DOT.time.dt = (50e3/1024/4);        % time step in picoseconds
DOT.time.nstep = 400;               % number of temporal steps
DOT.time.noise = 'Poisson';         % 'Poisson','Gaussian','none'
                                    % if 'Poisson' and sigma>0 a
                                    % Gaussian noise is added before
                                    % Poisson noise.
DOT.time.sigma = 1e-3;              % variance for gaussian noise
DOT.time.self_norm = false;         % true for self-normalized TPSF
DOT.time.TotCounts = 1e10;           % total counts for the maximum-energy
                                    % TPSF. The other are consequently
                                    % rescaled
%==========================================================================
%%                         Radiometry
%==========================================================================
DOT.radiometry.power = 1;    % (mW) laser input power %AAA
DOT.radiometry.timebin = ...
    DOT.time.dt;                % (ps) width of the time bin
DOT.radiometry.acqtime = 1;     % (s) acquisition time %AAA
DOT.radiometry.opteff = 0.9;    % Typical efficiency of the optical path
DOT.radiometry.lambda = 800;    % Wavelength (nm) 
                                % (just for calculation of input photons)
DOT.radiometry.area = 4;        % mm^2
DOT.radiometry.qeff = 0.05;     % Quantum efficiency
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