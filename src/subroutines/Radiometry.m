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


% % unitary parameters
% 
% DetAreaUnitary = 1; %(mm2) area of the detector unitary (will be adjusted for actual area after)
% SourcePowerUnitary = 1; % mW Power of the Unitary Source (will be adjusted for power later on)
% TaqUnitary = 1; %s Acquisition Time unitary
% OpticalEfficiency = radiometry.opteff; % Typical efficiency of the optical path
% LambdaUnitary = radiometry.lambda; %Wavelength (nm) (just for calculation of input photons)
% RepUnitary = 1E6; % (MHz) Reference Repetition rate Unitary %TOGLI QUESTA
% dtUnitary = 1;  % (ps) Unitary bin for the temporal axis
% QuantumEfficiency = radiometry.qeff;
% 
% % constants
% h=6.63E-34; % (m2.kg/s) Plank Constant
% c=3E8; % m/s
% mW_2_W=1E-3;
% %mm2_2_m2=1E-6;
% nm_2_m=1E-9;
% 
% % unitary factor
% phE = h*c./(LambdaUnitary*nm_2_m); % (J) photon energy 
% phN = SourcePowerUnitary*mW_2_W*TaqUnitary./phE; % Number of Injected Photons 
% ResponsivityUnitary = DetAreaUnitary*OpticalEfficiency; % mm2 
% RealFactor = phN.*ResponsivityUnitary*dtUnitary; % Note: you need to divide for PI to transform Reflectance into Radiance 
% 
% %attualize to effective parameters
% RealFactor = RealFactor.*radiometry.area*QuantumEfficiency*...
%     radiometry.timebin.*...
%     radiometry.power.*...
%     radiometry.acqtime;












lambda = [635; 670; 830; 915; 940; 980; 1030; 1065];
responsivity670 = 1.8e-6;
PDE670 = 0.08351;
const = responsivity670/PDE670;
[A] = xlsread('D:\programs\DOT\SOLUS\src\experimental\PDE.xlsx');
responsivity = zeros(size(lambda));   
for j = 1 : max(size(lambda))
    for i = 1 : max(size(A))
        if A(i,1) == lambda(j)
            responsivity(j,1) = const * A(i,2);
        end
    end
end



% unitary parameters
Area = radiometry.area; %mm2
power = radiometry.power'; 
TaqUnitary = 1; %s Acquisition Time unitary
LambdaUnitary = radiometry.lambda; %Wavelength (nm) (just for calculation of input photons)
dtUnitary = 1;  % (ps) Unitary bin for the temporal axis

% constants 
h=6.63E-34; % (m2.kg/s) Plank Constant
c=3E8; % m/s
mW_2_W=1E-3;
%mm2_2_m2=1E-6;
nm_2_m=1E-9;

phE = h*c./(LambdaUnitary*nm_2_m); % (J) photon energy 
phN = power*mW_2_W*TaqUnitary./phE; % Number of Injected Photons
phN = phN';
RealFactor = phN.*responsivity.*dtUnitary; % Note: you need to divide for
%PI to transform Reflectance into Radiance 

RealFactor = Area*RealFactor.*... 
    radiometry.timebin.*...
    radiometry.acqtime;

RealFactor = RealFactor'; 
