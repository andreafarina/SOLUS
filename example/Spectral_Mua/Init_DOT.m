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
EXP_DATA = 0;           % Load experimental data and use them for  ERA 0
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
REMOVE_VOXELS = 0; %if 1 the function Voxels_constant removes the most superficial voxels 
                   %(to remove the source-detector configuration pattern in the maps)
DOT.spe.cromo_label = {'Hb','HbO2','Lipid','H2O','Collagen'};
DOT.spe.active_cromo = [1,1,1,1,1];
DOT.spe.cromo_factor = [1,1,10*100*0.91,10*100,10*100*0.196];
% DOT.spe.cromo_factor = [1,1,1,1,1];
DOT.spe.cromo_units = {'microM','microM','mg/cm^3','mg/cm^3','mg/cm^3'};
DOT.spe.ForceConstitSolution = 0;
if SPECTRA == 0
mua_ = [0.00335394973467767,0.00208272121572369,0.00328333490577450,0.0101311816365722,0.0113533356313992,0.00950065021456835,0.00907688587133433,0.00617935497869766];
musp_ = [1.23400000000000,1.20133635788488,1.07935169073345,1.02799611957338,1.01423383723663,0.993319553161961,0.968909937001070,0.952855880308506];
Xd = {mua_,musp_};
else
a_ = 1.234; b_ = 0.50;
conc_ = [1.45 9.38 0.8186*10*100*0.91 0.1172*10*100 0.0453*10*100*0.196];
Xd = {conc_,[a_ b_]};
end
% ========================================================================= 
%% ====================== VOLUME DEFINITION ===============================
%% Background optical properties
% DOT.opt.muaB = 0.01;    % mm-1 %era commentato
% DOT.opt.muspB = 1;      % mm-1 %era commentato
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
DOT.opt.hete1.type  = {'Mua'};
DOT.opt.hete1.geometry = 'usimage';
DOT.opt.hete1.c     = [3, 4, 15];   % down 
% DOT.opt.hete1.d     = [0, 0, -1];   % down
% DOT.opt.hete1.l     = 20;
DOT.opt.hete1.sigma = 5; 
DOT.opt.hete1.distrib = 'OFF';
DOT.opt.hete1.profile = 'Gaussian';%'Step';%'Gaussian';
% DOT.opt.hete1.val   = 5 * DOT.opt.muaB;
if SPECTRA == 0
muap_=[0.0131844085972107,0.00811400967510129,0.00662420840605220,0.0128323383331261,0.0155414766883930,0.0187324032316907,0.0127792247467340,0.00824432384889349];
muspp_=[1.23400000000000,1.20133635788488,1.07935169073345,1.02799611957338,1.01423383723663,0.993319553161961,0.968909937001070,0.952855880308506];
Xp = {muap_,muspp_};
else
a_ = 1.634; b_ = 0.80; % a (mm-1), b(adimensionale)
concp_ = [9.36 16.22 0.5371*10*100*0.91 0.2539*10*100 0.1291*10*100*0.196]; 
Xp = {concp_};
end
% DOT.opt.hete1.path ='../3DMasks/Mask3D_Mask_malignant_4.mat' ;   % down
DOT.opt.hete1.path ='../../3DMasks/malignant_4.mat' ;   % down

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
DOT.time.dt = (50e3/4096/6)*4;   % time step in picoseconds
DOT.time.nstep = 1024;           % number of temporal steps
DOT.time.noise = 'Poisson';         % 'Poisson','Gaussian','none'
                                 % if 'Poisson' and sigma>0 a
                                 % Gaussian noise is added before
                                 % Poisson noise.
DOT.time.sigma = 1e-3;              % variance for gaussian noise
DOT.time.self_norm = false;         % true for self-normalized TPSF
DOT.time.TotCounts = 1e6*ones(1,8);  % total counts for the maximum-energy
                                      % TPSF. The other are consequently
                                      % rescaled
%==========================================================================
%%                         Radiometry
%==========================================================================
[P] = xlsread([getenv('DOTSRC'),filesep,'experimental',filesep,'LaserPower.xlsx'],'B2:B9'); 
DOT.radiometry.power = P;
DOT.radiometry.timebin = ...
    DOT.time.dt;                % (ps) width of the time bin
DOT.radiometry.acqtime = 1;    % (s) acquisition time
DOT.radiometry.opteff = 0.9;    % Typical efficiency of the optical path
DOT.radiometry.lambda = [635,670,830,915,940,980,1030,1065]; 
DOT.radiometry.lambda0 = 635;
DOT.radiometry.lambda = DOT.radiometry.lambda;
DOT.radiometry.nL = numel(DOT.radiometry.lambda);
                                % (just for calculation of input photons)
DOT.radiometry.area = 5.1083;      % mm^2
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
%                    to 2. In particular, if:
%                    RADIOMETRY==1, the count-rate for each delay is cut to 
%                         DOT.time.TotCounts if photons are available, 
%                         otherwise no;
%                   RADIOMETRY==0, the count-rate for each delay is cut to 
%                         DOT.time.TotCounts in any case.  
CUT_COUNTS = 1;         
NumDelays = 1;      % number of delays  
