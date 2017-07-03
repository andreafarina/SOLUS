%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MultiECBO                                                               %
%                                                                       %
% v1:                                                                   %
% Try to produce time course
%                                                                       %
% A. Pifferi 24/06/2017                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

IS_DIN=1;
FileNameDyn='DYN.mat';

nT=150; %s max macro time
dT=50; %s macro time delta over which to integrate

setPath;

for iDyn=1:(nT/dT)
    ExpScan;
    RhoZero_Scan;
end

save([PathServer,PathDataOut,FileNameDyn],'DYN');
