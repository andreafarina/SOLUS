%==========================================================================
%%  SETTING SOURCES (QVEC), DETECTORS (MVEC) AND THEIR PERMUTATIONS (DMASK)
%==========================================================================
% SOLUS SOURCES - DETECTOR POSITIONS
ys = [-10 10];
xs = linspace(-15,15,4);
%ys = linspace(5,50,8);
zs = 0;

[xxs,yys,zzs] = ndgrid(xs,ys,zs);

DOT.Source.Pos = [xxs(:),yys(:),zzs(:)];
DOT.Detector.Pos = [xxs(:),yys(:),zzs(:)];

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
