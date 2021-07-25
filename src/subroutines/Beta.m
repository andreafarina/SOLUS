function kap3D = Beta(refimage,th,dt_th)
%% Generate a diffusivity coefficient for the anisotropic laplacian
% th: threshold for the gradient
% 
dt_th = 1;
VERBOSITY = 1;
nvox = numel(refimage);
refimage = (refimage - min(refimage(:)))/(max(refimage(:)) - min(refimage(:)));

if sum(refimage(:) == 1) > sum(refimage(:) == 0) % inclusion should be smaller than the bulk
    refimage = 1 - refimage;
end

refimage = double(bwdist(logical(refimage)) <=dt_th*sqrt(3)); %pad image with one voxel
%refimage = smooth3(refimage,'gaussian',[9 9 9]);
gg = refimage;
Mgg = - 1e-2 + max(gg(:)); % 
if Mgg==0
    Mgg = 1; % to manage a flat mask
end
kappa = exp(- th*gg/Mgg); % factor 5 makes quite an extreme setting..
kap3D = spdiags(reshape(kappa,[],1),0:0,nvox,nvox);
%% plot
if VERBOSITY == 1
    figure(180),
    SubPlotMap(kappa,'Beta Map',181,1,1,1,[1 1 1]);
    drawnow
end