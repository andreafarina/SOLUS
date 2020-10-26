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

% lambda = [635; 670; 830; 915; 940; 980; 1030; 1065];
lambda = radiometry.lambda';
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
% power = radiometry.power(1);
TaqUnitary = 1; %s Acquisition Time unitary
LambdaUnitary = radiometry.lambda; %Wavelength (nm) (just for calculation of input photons)
dtUnitary = 1;  % (ps) Unitary bin for the temporal axis

% constants 
h=6.63E-34; % (m2.kg/s) Plank Constant
c=3E8; % m/s
mW_2_W=1E-3;
%mm2_2_m2=1E-6;
nm_2_m=1E-9;
m2_2_mm2=1e6;

phE = h*c./(LambdaUnitary*nm_2_m); % (J) photon energy 
phN = power*mW_2_W*TaqUnitary./phE; % Number of Injected Photons
phN = phN';
RealFactor = phN.*responsivity.*dtUnitary; % Note: you need to divide for
%PI to transform Reflectance into Radiance 

RealFactor = Area*RealFactor.*... 
    radiometry.timebin.*...
    radiometry.acqtime;
RealFactor = RealFactor.*m2_2_mm2;%need to get the same unit of measurement of RealFactor before the hardware SOLUS
RealFactor = RealFactor'; 
