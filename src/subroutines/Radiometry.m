function RealFactor = Radiometry(radiometry)
% Radiometry
% Converts simulated TPSF data to photon counts tacking into account
% radiometry is a structure with the following fields
%
% lambda:   reference wavelength
% power:    injected power (mW)
% acqtime:  total acquisition time (s)
% area:     detector area (mm2)
% opteff:   optical efficiency
% qeff:     quantum efficiency
%
% RealFactor: is the multiplication factor for converting 
%             TPSF units to photon-counts unit

% unitary parameters
DetAreaUnitary = 1; %(mm2) area of the detector unitary (will be adjusted for actual area after)
SourcePowerUnitary = 1; % mW Power of the Unitary Source (will be adjusted for power later on)
TaqUnitary = 1; %s Acquisition Time unitary
OpticalEfficiency = radiometry.opteff; % Typical efficiency of the optical path
LambdaUnitary = radiometry.lambda; %Wavelength (nm) (just for calculation of input photons)
%RepUnitary = 1E6; % (MHz) Reference Repetition rate Unitary
dtUnitary = 1;  % (ps) Unitary bin for the temporal axis
Area = radiometry.area; %mm
QuantumEfficiency = radiometry.qeff;

% constants
h=6.63E-34; % (m2.kg/s) Plank Constant
c=3E8; % m/s
mW_2_W=1E-3;
%mm2_2_m2=1E-6;
nm_2_m=1E-9;

% unitary factor
phE = h*c./(LambdaUnitary*nm_2_m); % (J) photon energy
phN = SourcePowerUnitary*mW_2_W*TaqUnitary./phE; % Number of Injected Photons
ResponsivityUnitary = DetAreaUnitary*OpticalEfficiency; % mm2
RealFactor = phN.*ResponsivityUnitary*dtUnitary; % Note: you need to divide for PI to transform Reflectance into Radiance

% attualize to effective parameters
RealFactor = RealFactor.*Area*QuantumEfficiency*...
    radiometry.timebin.*...
    radiometry.power.*...
    radiometry.acqtime;
