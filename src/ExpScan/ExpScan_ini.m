%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ini File for ExpScan                                                  %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CHOICES
IS_BKG(1)=1;
IS_BKG(2)=1;
iL=2;
BIN_X=2;



%% FILES
PathServer='D:\Beta\';
if BIN_X==2
    PathDataIn='Data\ScanHead\for_polimi_dat_1x\'; % use 'Data\ScanHead\for_polimi_dat_1x\' for 1x, 'Data\ScanHead\for_polimi_dat\' for 4x
else
    PathDataIn='Data\ScanHead\for_polimi_dat\'; % use 'Data\ScanHead\for_polimi_dat_1x\' for 1x, 'Data\ScanHead\for_polimi_dat\' for 4x
end    
PathDataOut='Simulations\Solus\data\201612\';
FileNameIn='S04_2_DTOFarray_6D_single_binxy_nobs_blocksum.mat';
if BIN_X==2
    FileNameOut='EXP_DATA_EXP_2x.mat';
else
    FileNameOut='EXP_DATA_EXP_4x.mat';
end   
FileNamePeak2='peak_001_2_1.sdt';

%% TIMING
DTOF_length=1024;
Gain=4;
Factor=50000/Gain/DTOF_length; % ps; range 50 ns, gain 4; ok for campaign 2013

%% LAMBDA
l(1) = 760; % first line in scan
l(2) = 860;

%% ROI
BkgChan{1}=[300:400];
BkgChan{2}=[300:400];
RoiChan{1}=[30:100];
RoiChan{2}=[90:160];

%% Time t0
DTOF_reduced = 301:700;

% REGION OF INTEREST
%EXP.time.roi1 = [1:64,65:128]; %bkg 256:end, good up to 200
%EXP.time.roi2 = [65:128,193:256]; %bkg 300:330, good up to 400 for bkg
%But it is better to start from 30, as from Heidrun (the previous zone is zer) and to avoid the secon part of the windows since already into noise. Therefore, just 2 windows of 64 bins:
% REGION OF INTEREST
%EXP.time.roi1 = [31:94]; %bkg 256:end, good up to 200
%EXP.time.roi2 = [95:158]; %bkg 300:330, good up to 400 for bkg
%Please remember that the two channels (early and late) are not aligned since there is around 700 pss time difference among the two. Therefore the two time axis must be calculated separately. (for rigorous 



%% ACTIVATION
rangeAct=20+[20:40,80:100];
rangeRest=20+[-10:10];





