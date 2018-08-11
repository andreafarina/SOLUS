%==========================================================================
%%                              OPTIONS 
%==========================================================================
% ----------------------------- FORWARD -----------------------------------
FORWARD = 0;            % Simulated forward data and save into _Data file
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
EXP_DATA = 1;           % Load experimental data and use them for 
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
DOT.spe.cromo_factor = [1,1,10*100*0.91,10*100, 10*100*0.196];
DOT.spe.cromo_units = {'microM','microM','mg/cm^3','mg/cm^3','mg/cm^3'};
DOT.spe.ForceConstitSolution = 0;
lamda_id = 1:8;
if SPECTRA == 0
mua_ = [0.00380793160000000,0.00134266120000000,0.00138762010000000,0.0100602733000000,0.00755439960000000,0.00394251390000000,0.00761375700000000,0.00493305200000000];
musp_ = [1.22842091900000,1.18990461900000,0.913511025000000,0.803037154800000,0.769676676200000,0.756078261900000,0.702287059500000,0.675179692900000];
mua_ = mua_(lamda_id); musp_ = musp_(lamda_id);
Xd = {mua_,musp_};
else
a_ = 1.5221;	b_ = 1.1415;
conc_ = [1.4492 0.48527 0.97134 0.037411 0.2021];
Xd = {conc_,[a_ b_]};
end
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
DOT.grid.dx = 2;

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
DOT.opt.hete1.c     = [0, -5, 10];   % down
% DOT.opt.hete1.d     = [0, 0, -1];   % down
% DOT.opt.hete1.l     = 20;
DOT.opt.hete1.sigma = 5;
DOT.opt.hete1.distrib = 'OFF';
DOT.opt.hete1.profile = 'Gaussian';%'Step';%'Gaussian';
DOT.opt.hete1.val   = 5 * DOT.opt.muaB;
if SPECTRA == 0
muap_ = [0.0373210250000000,0.0213434250000000,0.0104115095000000,0.0157479800000000,0.0230874250000000,0.0515664600000000,0.0294118300000000,0.0196425550000000];
muap_ = muap_(lamda_id);
muspp_ = [0.371159400000000,0.314731450000000,0.214725500000000,0.208025850000000,0.222794350000000,0.297037400000000,0.220381200000000,0.189757350000000];
Xp = {muap_,muspp_};
else
a_ = (1+10/100)*1.5453;	b_ = 1; % a (mm-1), b(adimensionale)
concp_ = 2.*[0.40855	0.71424 0.48844	0.5868	0.37051]; %microMolare
Xp = {concp_,[a_ b_]};
end
DOT.opt.hete1.path ='../3DMasks/Mask3D_Mask_malignant_4.mat' ;   % down

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
DOT.time.dt = (50e3/4096/6)*4;        % time step in picoseconds
DOT.time.nstep = 4096/4;               % number of temporal steps
DOT.time.noise = 'none';         % 'Poisson','Gaussian','none'
                                    % if 'Poisson' and sigma>0 a
                                    % Gaussian noise is added before
                                    % Poisson noise.
DOT.time.sigma = 1e-3;              % variance for gaussian noise
DOT.time.self_norm = true;         % true for self-normalized TPSF
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
DOT.radiometry.lambda = DOT.radiometry.lambda(lamda_id);
DOT.radiometry.nL = numel(DOT.radiometry.lambda);
                                % (just for calculation of input photons)
DOT.radiometry.area = 4;        % mm^2
DOT.radiometry.qeff = 0.05;     % Quantum efficiency
DOT.radiometry.lambda0 = 635;
if size(DOT.time.TotCounts,2)<DOT.radiometry.nL
    DOT.time.TotCounts = repmat(DOT.time.TotCounts,1,DOT.radiometry.nL);
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