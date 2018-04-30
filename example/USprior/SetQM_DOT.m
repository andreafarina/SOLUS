%==========================================================================
%%  SETTING SOURCES (QVEC), DETECTORS (MVEC) AND THEIR PERMUTATIONS (DMASK)
%==========================================================================
% SOLUS SOURCES POSITIONS
relSource = 0;% 0.5 % 1% only 1 out of 7 wavelegths is aligned to the detectors
xs = [-14.95 , 14.95];%[-15 , 15];%
ys = linspace(-19.5,19.5,4) + relSource;%[-20, -7, 7, 20];%
%ys = linspace(5,50,8);
zs = 0;

[yys,xxs,zzs] = ndgrid(ys,xs,zs);

DOT.Source.Pos = [xxs(:),yys(:),zzs(:)];

% SOLUS DETECTORS POSITIONS
xd = [-10.95 , 10.95];%[-11 , 11];%
yd = linspace(-19.5, 19.5 , 4);%[-20, -7, 7, 20];%
zd = 0;

[yyd,xxd,zzd] = ndgrid(yd,xd,zd);

DOT.Detector.Pos = [xxd(:),yyd(:),zzd(:)];

% non-contact PTB setup 40x40 mm2 scan with 8x8 and s-d 5mm
% rhozero
% rhosd = 5;
% % 8x8
% % DOT.Source.Pos = RasterScan(-20-rhosd/2,20-rhosd/2,-20,20,8,8,0);
% % DOT.Detector.Pos = RasterScan(-20+rhosd/2,20+rhosd/2,-20,20,8,8,0);
% % 16x16
% DOT.Source.Pos = RasterScan(-20-rhosd/2,20-rhosd/2,-20,20,16,16,0);
% DOT.Detector.Pos = RasterScan(-20+rhosd/2,20+rhosd/2,-20,20,16,16,0);

DOT.Source.Ns=size(DOT.Source.Pos,1);
DOT.Detector.Nd=size(DOT.Detector.Pos,1);
%% Define permutation matrix
% ALL COMBINATIONS: null-distances + all the other combinations
%DOT.dmask = logical(ones(DOT.Detector.Nd,DOT.Source.Ns));

% NULL-DISTANCE ONLY
%DOT.dmask = logical(eye(DOT.Detector.Nd,DOT.Source.Ns));

% ALL EXCEPT NULL-DISTANCE
DOT.dmask = logical(ones(DOT.Detector.Nd,DOT.Source.Ns) - ...
    diag(diag(ones(DOT.Detector.Nd,DOT.Source.Ns))));              
% -------------------------------------------------------------------------
