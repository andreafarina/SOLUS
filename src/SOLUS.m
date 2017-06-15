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
% Solver: Tickhonov, Truhncated SVD, pcg, gmres                          %
%                                                                       %
% A. Farina 26/12/2016                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
%clear all;
setPath;
%addpath toastUtil subroutines solvers
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
% exp_path = ['../data/',session,'/'];
% res_path = ['../results/',session,'/'];

exp_path = ['D:/Beta/Simulations/Solus/data/',session,'/'];
res_path = ['D:/Beta/Simulations/Solus/results/',session,'/'];
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
EXP_DELTA = 'peak';    % Substitute the IRF with delta function on the 
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


%% ====================== VOLUME DEFINITION ===============================
%% Background optical properties
DOT.opt.muaB = 0.01;          % mm-1
DOT.opt.muspB = 1.0;             % mm-1
%vi = 1000;          % (mm3) volume of the inclusion
DOT.opt.nB = 1.4;
DOT.opt.nE = 1.;   % external refractive index

%==========================================================================
%%                                  SET GRID
%==========================================================================
DOT.grid.x1 = 0;%0
DOT.grid.x2 = 64;%64
DOT.grid.dx = 2;%2

DOT.grid.y1 = 0;%0
DOT.grid.y2 = 58;%58           

DOT.grid.z1 = 0;%0        
DOT.grid.z2 = 32;%32         
%==========================================================================
%%                      Set Heterogeneities
%==========================================================================
%--------------------------- INCLUSION 1 ---------------------------------%
DOT.opt.hete.type  = 'Mua';
DOT.opt.hete.geometry = 'Sphere';
DOT.opt.hete.c     = [35, 25, 16];   % down
DOT.opt.hete.sigma = 5;
DOT.opt.hete.distrib = 'OFF';
DOT.opt.hete.profile = 'Step';%'Gaussian';
DOT.opt.hete.val   = 2. * DOT.opt.muaB;
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
%%                         Radiometry
%==========================================================================
%InputParm;
timebin=10; % (ps) width of the time bin %AAA later is DOT.time.dt
power=1; % (mW) laser input power %AAA
acqtime=1; % (s) acquisition time %AAA

% unitary parameters
DetAreaUnitary=1; %(mm2) area of the detector unitary (will be adjusted for actual area after)
SourcePowerUnitary=1; % mW Power of the Unitary Source (will be adjusted for power later on)
TaqUnitary=1; %s Acquisition Time unitary
OpticalEfficiency=0.9; % Typical efficiency of the optical path
LambdaUnitary=800; %Wavelength (nm) (just for calculation of input photons)
RepUnitary=1E6; % (MHz) Reference Repetition rate Unitary
dtUnitary=1;  % (ps) Unitary bin for the temporal axis
Area=4; %mm
QuantumEfficiency=0.05;

%==========================================================================
%%                         Time domain parameters
%==========================================================================
DOT.time.dt = 10;    % time step in picoseconds
DOT.time.maxtime = 4000; % max time in picosecond %AAA
DOT.time.noise = 'Poisson';         % 'Poisson','Gaussian','none'
                                    % if 'Poisson' and sigma>0 a
                                    % Gaussian noise is added before
                                    % Poisson noise.
DOT.time.sigma = 1e-2;              % variance for gaussian noise
DOT.time.self_norm = false;         % true for self-normalized TPSF   
DOT.time.TotCounts = 1e6;           % total counts for the maximum-energy
                                    % TPSF. The other are consequently
                                    % rescaled

%==========================================================================
%%  SETTING SOURCES (QVEC), DETECTORS (MVEC) AND THEIR PERMUTATIONS (DMASK)
%==========================================================================
% SOLUS SOURCES - DETECTOR POSITIONS
xs = linspace(-25+6,25-6,4);
ys = [-5-6-3,5+6+3];
zs = 0;

xd = linspace(-25+6,25-6,4);
yd = [-5-6+3,5+6-3];
zd = 0;

%% Define permutation matrix
% Source Detector arrangement (1=all; 2=only rhozero; 3=no rhozero)
SourceDetMatrix=1; % ALL COMBINATIONS

%% Type of Cutting for statistics
% Cutting of counts (1=gated; 2=no gated; 3=each gate=countrate/numgate)
CutCounts = 1; % all gated
NumDelays = 3; % number of delays
NumGates = 10; % number of gates

%% Regularization Parameters
% ---------------------- Solver and regularization ------------------------
SolverTau = 0.01;  % regularisation parameter
SolverBeta = 0.01;           % TV regularisation parameter
SolverTol = 1e-6;            % Gauss-Newton convergence criterion
SolverTolKrylov = 1e-6;      % Krylov convergence criterion
SolverItrmax = 1;            % Gauss-Newton max. iterations
SolverType = 'born';           % 'born','GN': gauss-newton, 



%% END OF INPUT FOR SOLUS

%
%
%
%
%

% ========================================================================= 
%%                          OVERRIDE param using P values (multiSIM) 
% =========================================================================

% Background optical properties
if exist('P(MUAB).Value','var'), DOT.opt.muaB = eval('P(MUAB).Value'); end            % mm-1
if exist('P(MUSB).Value','var'), DOT.opt.muspB = eval('P(MUSB).Value'); end            % mm-1

% Grid
if exist('P(X1).Value','var'), DOT.grid.x1 = eval('P(X1).Value'); end %0
if exist('P(X2).Value','var'), DOT.grid.x2 = eval('P(X2).Value'); end %64
if exist('P(DX).Value','var'), DOT.grid.dx = eval('P(DX).Value'); end %2
if exist('P(Y1).Value','var'), DOT.grid.y1 = eval('P(Y1).Value'); end %0
if exist('P(Y2).Value','var'), DOT.grid.y2 = eval('P(Y2).Value'); end %58           
if exist('P(Z1).Value','var'), DOT.grid.z1 = eval('P(Z1).Value'); end %0        
if exist('P(Z2).Value','var'), DOT.grid.z2 = eval('P(Z2).Value'); end %32         

% Inclusion 1
if exist('P(XP).Value','var'), DOT.opt.hete.c = eval('[P(XP).Value, P(YP).Value, P(ZP).Value]'); end % down
if exist('P(MUAP).Value','var'), DOT.opt.hete.val = eval('P(MUAP).Value')+DOT.opt.muaB; end

% Radiometry
if exist('P(DT).Value','var'), timebin=eval('P(DT).Value'); end % (ps) time bin
if exist('P(PW).Value','var'), power=eval('P(DT).Value'); end % (mW) injected laser power
if exist('P(TA).Value','var'), acqtime=eval('P(TA).Value'); end % (s) acquisition time

% time domain parameters
if exist('P(DT).Value','var'), DOT.time.dt = eval('P(DT).Value'); end % time step in picoseconds
if exist('P(MT).Value','var'), DOT.time.maxtime = eval('P(MT).Value'); end % max time in picosecond
if exist('P(CR).Value','var'), DOT.time.TotCounts = eval('P(CR).Value'); end % total counts for the maximum-energy
    
% permutation matrix
if exist('P(SD).Value','var'), SourceDetMatrix = eval('P(SD).Value'); end

% Type of Cutting for statistics
if exist('P(CUT).Value','var'), CutCounts = eval('P(CUT).Value'); end; % all gated
if exist('P(ND).Value','var'), NumDelays = eval('P(ND).Value'); end; % number of delays
if exist('P(NG).Value','var'), NumGates = eval('P(NG).Value'); end; % number of gates

% Regularization Parameters
if exist('P(TAU).Value','var'), SolverTau = eval('P(TAU).Value'); end;  % regularisation parameter

%% END OF OVERRIDE P FOR MULTISIM

%
%
%
%
%

%% START SOLUS PROGRAM

% =========================================================================
%%                     Create folder for saving data
% =========================================================================
rdir = res_path;
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
%     if (EXP_DMD == 1)
%      %   DOT.DMDpath = [EXP.path.data_folder,EXP.path.day,filesep,...
%      %       EXP.path.calib_dir,filesep];
%         DOT.DMDpath = [exp_path,'Calibration',filesep]; 
%         %DOT.DMDpath = './experimental/data/201510_Calibration/';
%         DOT.DMDfile = EXP.path.calib_file;
%     end
end

%==========================================================================


%% ====================== VOLUME DEFINITION ===============================
%% Background optical properties
DOT.opt.cs = 0.299/DOT.opt.nB;       % speed of light in medium
DOT.opt.kap = 1/(3*DOT.opt.muspB);
DOT.A = A_factor(DOT.opt.nB/DOT.opt.nE); % A factor for boundary conditions

%==========================================================================
%%                                  SET GRID
%==========================================================================
DOT.grid.dy = DOT.grid.dx;
DOT.grid.dz = DOT.grid.dx;
DOT.grid = setGrid(DOT); 
DOT.opt.Mua = ones(DOT.grid.dim) * DOT.opt.muaB;
DOT.opt.Musp = ones(DOT.grid.dim) * DOT.opt.muspB;

%==========================================================================
%%                      Set Heterogeneities
%==========================================================================
%--------------------------- INCLUSION 1 ---------------------------------%
[DOT,DOT.opt.hete] = setHete(DOT,DOT.opt.hete);
%--------------------------- INCLUSION 2 ---------------------------------%
% [DOT,DOT.opt.hete2] = setHete(DOT,DOT.opt.hete2);

%==========================================================================
%%                         Radiometry
%==========================================================================
% constants
h=6.63E-34; % (m2.kg/s) Plank Constant
c=3E8; % m/s
mW_2_W=1E-3;
mm2_2_m2=1E-6;
nm_2_m=1E-9;

% unitary factor
phE=h*c/(LambdaUnitary*nm_2_m); % (J) photon energy
phN=SourcePowerUnitary*mW_2_W*TaqUnitary/phE; % Number of Injected Photons
ResponsivityUnitary=DetAreaUnitary*OpticalEfficiency; % mm2
RealFactor=phN*ResponsivityUnitary*dtUnitary; % Note: you need to divide for PI to transform Reflectance into Radiance

% attualize to effective parameters
RealFactor=RealFactor*Area*QuantumEfficiency*timebin*power*acqtime;
%==========================================================================
%%                         Time domain parameters
%==========================================================================
if DOT.TD == 1
    DOT.time.nstep = floor(DOT.time.maxtime/DOT.time.dt)+1;    % number of temporal steps        260

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
    % PATCH IRF
    temp=load('SIPMNew.txt');
    EXP.irf.data=temp(:,2);
    EXP.time.axis=temp(:,1);
    % END PATCH IRF
    %DOT.time.irf.data = resampleSPC(EXP.irf.data,EXP.time.axis,DOT.time.dt,'norm');
    DOT.time.irf.data = resampleSPCinterp(EXP.irf.data,EXP.time.axis,DOT.time.dt,'norm');
    %DOT.time.irf.data = interp1(EXP.time.axis,EXP.irf.data,EXP.time.axis(1):DOT.time.dt:EXP.time.axis(end));
    [MaxIRF,Chan0]=max(DOT.time.irf.data);
else
    DOT.time.irf.data = 0;
    Chan0=1;
end
end

%==========================================================================
%%  SETTING SOURCES (QVEC), DETECTORS (MVEC) AND THEIR PERMUTATIONS (DMASK)
%==========================================================================
% SOLUS SOURCES - DETECTOR POSITIONS
xs = linspace(-25+6,25-6,4);
ys = [-5-6-3,5+6+3];
zs = 0;

xd = linspace(-25+6,25-6,4);
yd = [-5-6+3,5+6-3];
zd = 0;

[xxs,yys,zzs] = ndgrid(xs,ys,zs);
[xxd,yyd,zzd] = ndgrid(xd,yd,zd);

DOT.Source.Pos = [xxs(:),yys(:),zzs(:)];
DOT.Detector.Pos = [xxd(:),yyd(:),zzd(:)];
DOT.Source.Ns=size(DOT.Source.Pos,1);
DOT.Detector.Nd=size(DOT.Detector.Pos,1);

%% Define permutation matrix
% Source Detector arrangement (1=all; 2=only rhozero; 3=no rhozero)
% ALL COMBINATIONS
if SourceDetMatrix==1, DOT.dmask = logical(ones(DOT.Detector.Nd,DOT.Source.Ns)); end
% NULL-DISTANCE ONLY
if SourceDetMatrix==2, DOT.dmask = logical(eye(DOT.Detector.Nd,DOT.Source.Ns)); end
% ALL EXCEPT NULL-DISTANCE
if SourceDetMatrix==3, DOT.dmask = logical(ones(DOT.Detector.Nd,DOT.Source.Ns) - ...
    diag(diag(ones(DOT.Detector.Nd,DOT.Source.Ns))));                       end              
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
    figure,suptitle('Data CW')
    ShowData(DOT.dmask,(DataCW(:)'));
    if REF == 1
        RefCW  = ForwardCW(DOT.grid,DOT.Source.Pos, DOT.Detector.Pos, DOT.dmask,...
            DOT.opt.muaB, DOT.opt.muspB, DOT.opt.muaB*ones(DOT.grid.dim),...
                DOT.opt.muspB*ones(DOT.grid.dim), DOT.A, geom, 'homo');
        [RefCW,RefsdCW] = AddNoise(RefCW,'gaussian',DOT.sigma);
        figure, suptitle('CW Born contrast')
        ShowData(DOT.dmask,(DataCW(:)' - RefCW(:)')./RefCW(:)');
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

    %-------------------- Add  noise to TD data ------------------------
    sdTD = ones(size(DataTD));
    if ~strcmpi(DOT.time.noise,'none')
        [DataTD,sdTD] = AddNoise(DataTD,'gaussian',DOT.time.sigma);
        if REF == 1
            RefTD = AddNoise(RefTD,'gaussian',DOT.time.sigma);
        end
   end
    if strcmpi(DOT.time.noise,'poisson')
        if (REF == 1) %&& (DOT.time.self_norm == 0))
            % Cutting of counts (1=gated; 2=no gated; 3=each gate=countrate/numgate)
            RefTD = RefTD*RealFactor;
            DataTD = DataTD*RealFactor;
            
            % 1=GATED
            if CutCounts==1 
                for id=1:NumDelays     
                    deltachan=floor((DOT.time.nstep-1)/NumDelays);
                    startdelay=(id-1)*deltachan;
                    Chan1=Chan0+startdelay;
                    Chan2=Chan0+DOT.time.nstep-1;
                    Chan3=Chan0+id*deltachan-1;
                    area=sum(RefTD(Chan1:Chan2,:),1);
                    %maxrate=P(CR).Value*P(TA).Value/P(NG).Value*ones(size(area)); % divide by NG to be fair with other approaches
                    maxrate=DOT.time.TotCounts*acqtime*ones(size(area));
                    cutfactor=min(area,maxrate)./area./NumDelays;
                    RefTD(Chan1:Chan3,:)=bsxfun(@times,RefTD(Chan1:Chan3,:),cutfactor);
                    DataTD(Chan1:Chan3,:)=bsxfun(@times,DataTD(Chan1:Chan3,:),cutfactor);
                end
            end
            
            % 2=NOT GATED
            if CutCounts==2 
                area=sum(RefTD,1);
                cutfactor=min(area,DOT.time.TotCounts*acqtime*ones(size(area)))./area;
                RefTD=bsxfun(@times,RefTD,cutfactor);
                DataTD=bsxfun(@times,DataTD,cutfactor);
            end
            
            % 3=EVENLY DISTRIBUTED OVER ALL GATES
            if CutCounts==3,
                cutfactor=(DOT.time.TotCounts*acqtime)/DOT.time.nstep./RefTD;
                cutfactor(RefTD==0)=0;
                RefTD=bsxfun(@times,RefTD,cutfactor);
                RefTD=bsxfun(@times,RefTD,cutfactor);
            end
            
            % Add Poisson Noise
            RefTD = poissrnd(round(full(RefTD)));
            DataTD = poissrnd(round(full(DataTD)));

        % NOTE: THIS SECTION IS NOT ADJUSTED FOR THE DIFFERENT CONDITIONS
        else
            [m,j] = max(DataCW);
            %factor = 1./sum(DataTD(:,j)) .* DOT.time.TotCounts;
            factor = RealFactor;
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
disp(['Results will be stored in: ', rdir,filename,'_'])
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
%         figure(200);
%         for k=1:nview
%             subplot(ceil(sqrt(nview/1.5)),ceil(1.5*sqrt(nview/1.5)),k);
%             y = EXP.data.spc(:,(1:nqm)+(k-1)*nqm)-EXP.data.ref(:,(1:nqm)+(k-1)*nqm);
%             %y = y./EXP.data.ref(:,(1:nqm)+(k-1)*nqm);
%             y(isnan(y)) = [];
%             plot(y)%,ylim([min(y(:)) max(y(:))*1.1]),xlabel('time steps');
%             title(['Exp data - TD K-space view ',num2str(k)]);
%             clear y;
%         end
%         
% % ----- Resample experimental data to match the time-scale of DOT.time ----
        for i = 1:size(EXP.data.spc,2)
            z(:,i) = resampleSPC(EXP.data.spc(:,i),EXP.time.axis,DOT.time.dt);
            f(:,i) = resampleSPC(EXP.data.ref(:,i),EXP.time.axis,DOT.time.dt);
        end
%         figure(201);
%         for k=1:nview
%             subplot(ceil(sqrt(nview/1.5)),ceil(1.5*sqrt(nview/1.5)),k);
%             %y = EXP.data.spc(:,(1:nqm)+(k-1)*nqm)-EXP.data.ref(:,(1:nqm)+(k-1)*nqm);
%             y = z(:,(1:nqm)+(k-1)*nqm)-f(:,(1:nqm)+(k-1)*nqm);%./f(:,(1:nqm)+(k-1)*nqm);
%             y(isnan(y)) = [];
%             plot(y),%ylim([min(y(:)) max(y(:))]),xlabel('time steps');
%             title(['Exp resampled-data TD K-space view ',num2str(k)]);
%             clear y;
%         end
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

if RECONSTRUCTION == 1
%disp('-------------------------------------------------------------------')
%disp('                         RECONSTRUCTION                            ')
%disp('-------------------------------------------------------------------')
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
%%                               OPTIONS
%==========================================================================
NEW_GRID = false;           % true to load a new grid (typically coarser)
%==========================================================================
%%                               Load data
%==========================================================================
if (EXP_DATA == 0)
    %if isfield(REC,'Datafile')
    %    load([REC.Datapath,REC.Datafile]);
    %else
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
      
  %twin = CreateTimeWindows(REC.time.nstep,[10,REC.time.nstep],'even',20);
  twin = CreateTimeWindows(REC.time.nstep,[1,REC.time.nstep],'even',NumGates);
  REC.time.twin = twin + Chan0-1; % Chan0 is IRF peak channel, add -1 since twin starts from 1
  REC.time.nwin = size(REC.time.twin,1);

  % plot roi on the first measruement
  figure(1000);
  %semilogy(DataTD),ylim([max(DataTD(:))/1E9 max(DataTD(:))])
  semilogy(DataTD),ylim([1 max(DataTD(:))])
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
%==========================================================================
%%                  NEW GRID (OPTIONAL)
%==========================================================================
if NEW_GRID == true
    disp('----- REGULAR GRID MAPPING FOR RECONSTRUCTION-----'); tic;
    %-- Regular mesh --%
    REC.grid.x1 = 0;           
    REC.grid.x2 = 64;            
    REC.grid.dx = 4;
    
    %
    REC.grid.y1 = 0;           
    REC.grid.y2 = 58;            
    REC.grid.dy = REC.grid.dx;
    
    REC.grid.z1 = 0;          
    REC.grid.z2 = 32;          
    REC.grid.dz = REC.grid.dx;
    
    REC.grid = setGrid(REC); toc;
    REC.opt = rmfield(REC.opt,'Mua');
    REC.opt = rmfield(REC.opt,'Musp');
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
REC.solver.tau = SolverTau;       % regularisation parameter
REC.solver.beta = SolverBeta;     % TV regularisation parameter
REC.solver.tol = SolverTol;       % Gauss-Newton convergence criterion
REC.solver.tolKrylov = SolverTolKrylov;      % Krylov convergence criterion
REC.solver.itrmax = SolverItrmax;            % Gauss-Newton max. iterations
REC.solver.type = SolverType;           % 'born','GN': gauss-newton, 
                                  % 'LM': Levenberg-Marquardt,
                                  % 'fit': fitting homogeneous data
REC.solver.prior = [];%REC.opt.Mua;   % dictionnary matrix, or []
%REC.solver.Himplicit = true;                  
REC.solver.filename = filename;
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
%    set(gcf,'Position',get(0,'ScreenSize'));
max(REC.opt.Mua(:));
min(REC.opt.Mua(:));
    ShowRecResults(REC.grid,REC.opt.Mua,...
        REC.grid.z1,REC.grid.z2,REC.grid.dz,1,...
          0,max(REC.opt.Mua(:)));
    suptitle('Mua');
   % export_fig '../Results/20151116/ref_mua_1incl.pdf' -transparent
   %save_figure([rdir,filename,'_', 'mua_ref']);
end
% ------------------------ Reference musp ---------------------------------
if isfield(REC.opt,'Musp')
    figure(302);
%    set(gcf,'Position',get(0,'ScreenSize'));

    ShowRecResults(REC.grid,REC.opt.Musp,...
        REC.grid.z1,REC.grid.z2,REC.grid.dz,1,...
        0,max(REC.opt.Musp(:)));
    suptitle('Mus');
   % export_fig '../Results/20151116/ref_mus_1incl.pdf' -transparent
   %save_figure([rdir,filename,'_', 'mus_ref']);
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
                [REC.opt.bmua,REC.opt.bmusp] = RecSolverBORN_CW(REC.solver,...
                    REC.grid,REC.opt.mua0,REC.opt.musp0,REC.A,...
                    REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.Data,REC.ref,0);
                        
            case 'fit'
                [REC.opt.bmua,REC.opt.bmusp] = FitMuaMusAmpl_CW(REC.solver,REC.mesh.hMesh,...
                    REC.grid.hBasis,REC.opt.notroi,...
                    REC.opt.mua0,REC.opt.musp0,REC.opt.nB,...
                    REC.Q,REC.M,REC.dmask,REC.Data,REC.ref,0);
                
                
        end
    case 'td'
        switch lower(REC.solver.type)
            case 'gn'
                
            case 'born'
                [REC.opt.bmua,REC.opt.bmusp] = RecSolverBORN_TD(REC.solver,...
                    REC.grid,...
                    REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                    REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                    REC.time.twin,REC.time.self_norm,REC.Data,...
                    REC.time.irf.data,REC.ref,REC.sd,0);
                
            case 'fit'
                [REC.opt.bmua,REC.opt.bmusp] = FitMuaMus_TD(REC.solver,REC.mesh.hMesh,...
                    REC.grid.hBasis,REC.opt.notroi,...
                    REC.opt.mua0,REC.opt.musp0,REC.opt.nB,...
                    REC.Q,REC.M,REC.dmask,REC.time.dt,REC.time.nstep,...
                    REC.time.twin,REC.time.self_norm, REC.ref,REC.sd,...
                    REC.time.irf.data,1);
                
            case 'cg'
            
            case 'l1'
                [REC.opt.bmua,REC.opt.bmusp] = RecSolverL1_TD(REC.solver,...
                    REC.grid,...
                    REC.opt.mua0,REC.opt.musp0,REC.opt.nB,REC.A,...
                    REC.Source.Pos,REC.Detector.Pos,REC.dmask,REC.time.dt,REC.time.nstep,...
                    REC.time.twin,REC.time.self_norm,REC.Data,...
                    REC.time.irf.data,REC.ref,REC.sd,0);
                
                
        end
end
        %% display the reconstructions
% if ~isfield(REC,'Datafile') % if simulated I have reference values
%     clim_mua = [min(REC.opt.Mua(:)),max(REC.opt.Mua(:))];
%     clim_musp = [min(REC.opt.Musp(:)),max(REC.opt.Musp(:))];
% else
%    clim_mua = [min(REC.opt.bmua),max(REC.opt.bmua)];
%    clim_musp = [min(REC.opt.bmusp),max(REC.opt.bmusp)];
% end
% ---------------------------- display mua --------------------------------
drawnow;
figure(304);
%set(gcf,'Position',get(0,'ScreenSize'));
ShowRecResults(REC.grid,reshape(REC.opt.bmua,REC.grid.dim),...
   REC.grid.z1,REC.grid.z2,REC.grid.dz,1,0.00,0.05);
suptitle('Recon Mua');
pause(5);
%export_fig '../Results/20151130/recGN_mua_2incl_late.pdf' -transparent
%SSsave_figure([rdir,filename,'_', 'mua_rec']);
% ---------------------------- display musp -------------------------------
figure(305);
%set(gcf,'Position',get(0,'ScreenSize'));
ShowRecResults(REC.grid,reshape(REC.opt.bmusp,REC.grid.dim),...
   REC.grid.z1,REC.grid.z2,REC.grid.dz,1);%,0.,0.64);
suptitle('Recon Mus');
%export_fig '../Results/20151130/recGN_mus_2incl_late.pdf' -transparent   
%save_figure([rdir,filename,'_', 'mus_rec']);
drawnow;
tilefigs;
disp('recon: finished')
toc;
% =========================================================================
%%                            Save reconstruction 
% =========================================================================
save([rdir,filename,'_', 'REC'],'REC')
%save('CW','dxx','erri','stepi')

% =========================================================================
%%                            Fit Gauss reconstruction 
% =========================================================================
% xData=REC.grid.x;
% yData=REC.grid.y;
% zData=REC.grid.z;
% zslice=REC.grid.z(5);
% vData=reshape(REC.opt.bmua,REC.grid.dim);
% FitGauss3D(xData,yData,zData,vData);

% =========================================================================
%%                            Quantify DOT 
% =========================================================================

%% Prepare variables
%[X_in,Y_in,Z_in] = meshgrid(REC.grid.x,REC.grid.y,REC.grid.z);
[X_in,Y_in,Z_in] = ndgrid(REC.grid.x,REC.grid.y,REC.grid.z);
dmua_rec3d=reshape(REC.opt.bmua-REC.opt.muaB,REC.grid.dim); %% subtract background
dmua_true3d=REC.opt.Mua-REC.opt.muaB; %% subtract background

%% Prepare Mask
radius_max=20; %max radius of region of interest;
temp=dmua_rec3d;
temp(:,:,1:2)=0;
[mxm,mind] = max(temp(:));
%mind = find(dmua_rec3d == mxm);
[Xrec,Yrec,Zrec]=ind2sub(size(dmua_rec3d),mind);

V=zeros(size(dmua_rec3d));
nzm = find(sqrt((X_in-Xrec).^2+(Y_in-Yrec).^2+(Z_in-Zrec).^2)<radius_max);
V(nzm) = 1;

%% Centre of Mass
COM_TRUE = zeros(3,1);
vv = sum(dmua_true3d(:));
COM_TRUE(1) = sum(dmua_true3d(:).*X_in(:))/vv;
COM_TRUE(2) = sum(dmua_true3d(:).*Y_in(:))/vv;
COM_TRUE(3) = sum(dmua_true3d(:).*Z_in(:))/vv;
disp(['TRUE centre of mass [',num2str(COM_TRUE(1)),',',num2str(COM_TRUE(2)),',',num2str(COM_TRUE(3)),']']);
 
COM_REC = zeros(3,1);
vr = sum(dmua_rec3d(:));
COM_REC(1) = sum(dmua_rec3d(:).*X_in(:))/vr;
COM_REC(2) = sum(dmua_rec3d(:).*Y_in(:))/vr;
COM_REC(3) = sum(dmua_rec3d(:).*Z_in(:))/vr;
disp(['REC centre of mass [',num2str(COM_REC(1)),',',num2str(COM_REC(2)),',',num2str(COM_REC(3)),']']);
 
COM_err = norm(COM_REC-COM_TRUE);
disp(['Difference in centre of mass [',num2str(COM_REC(1)-COM_TRUE(1)),',',num2str(COM_REC(2)-COM_TRUE(2)),',',num2str(COM_REC(3)-COM_TRUE(3)),']']);
disp(['Error in centre of mass ',num2str(COM_err)]);
F(errXP).Value(iL)=COM_REC(1)-COM_TRUE(1);
F(errYP).Value(iL)=COM_REC(2)-COM_TRUE(2);
F(errZP).Value(iL)=COM_REC(3)-COM_TRUE(3);

COM_REC = zeros(3,1);
 vr = sum(V(:).*dmua_rec3d(:));
 COM_REC(1) = sum(V(:).*dmua_rec3d(:).*X_in(:))/vr;
 COM_REC(2) = sum(V(:).*dmua_rec3d(:).*Y_in(:))/vr;
 COM_REC(3) = sum(V(:).*dmua_rec3d(:).*Z_in(:))/vr;
COM_err = norm(COM_REC-COM_TRUE);
disp(['Difference in centre of mass in region [',num2str(COM_REC(1)-COM_TRUE(1)),',',num2str(COM_REC(2)-COM_TRUE(2)),',',num2str(COM_REC(3)-COM_TRUE(3)),']']);
disp(['Error in centre of mass in region ',num2str(COM_err)]);

%% covariance of mass
VAR_TRUE = zeros(3,1);
VAR_TRUE(1) = sum(dmua_true3d(:).*((X_in(:)-COM_TRUE(1)).^2 ))/vv;
VAR_TRUE(2) = sum(dmua_true3d(:).*((Y_in(:)-COM_TRUE(2)).^2 ))/vv;
VAR_TRUE(3) = sum(dmua_true3d(:).*((Z_in(:)-COM_TRUE(3)).^2 ))/vv;

VAR_REC = zeros(3,1);
vr = sum(dmua_rec3d(:));
VAR_REC(1) = sum(dmua_rec3d(:).*((X_in(:)-COM_REC(1)).^2 ))/vr;
VAR_REC(2) = sum(dmua_rec3d(:).*((Y_in(:)-COM_REC(2)).^2 ))/vr;
VAR_REC(3) = sum(dmua_rec3d(:).*((Z_in(:)-COM_REC(3)).^2 ))/vr;

disp(['Sqrt of covariance of mass [',num2str(sqrt(VAR_REC(1))),',',num2str(sqrt(VAR_REC(2))),',',num2str(sqrt(VAR_REC(3))),']']);

vr = sum(V(:).*dmua_rec3d(:));
VAR_REC(1) = sum(V(:).*dmua_rec3d(:).*((X_in(:)-COM_REC(1)).^2 ))/vr;
VAR_REC(2) = sum(V(:).*dmua_rec3d(:).*((Y_in(:)-COM_REC(2)).^2 ))/vr;
VAR_REC(3) = sum(V(:).*dmua_rec3d(:).*((Z_in(:)-COM_REC(3)).^2 ))/vr;
disp(['Sqrt of covariance of mass in region [',num2str(sqrt(VAR_REC(1))),',',num2str(sqrt(VAR_REC(2))),',',num2str(sqrt(VAR_REC(3))),']']);

%% Sigma and CNR
sigma=std(dmua_rec3d(:));
disp(['Standard deviation ',num2str(sigma)]);

av=max(dmua_rec3d(:));
disp(['Max ',num2str(av)]);
disp(['CNR ',num2str(av./sigma)]);

av=max(V(:).*dmua_rec3d(:));
disp(['Max in region ',num2str(av)]);
disp(['CNR in region ',num2str(av./sigma)]);


%% Quantitation
factor=REC.opt.hete.sigma^3/(REC.grid.dx*REC.grid.dy*REC.grid.dz)*(4/3*pi);
DmuaVolTrue=sum(dmua_true3d(:))/factor;
DmuaVolRec=sum(dmua_rec3d(:))/factor;
disp(['DmuaVolTrue ',num2str(DmuaVolTrue)]);
disp(['DmuaVolRec ',num2str(DmuaVolRec)]);
disp(['error ',num2str(DmuaVolRec/DmuaVolTrue-1)]);
F(relMUAP).Value(iL)=(DmuaVolRec-DmuaVolTrue)/DmuaVolTrue;


end
