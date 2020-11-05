function kap3D = Kappa(refimage,th)
%% Generate a diffusivity coefficient for the anisotropic laplacian
% th: threshold for the gradient
% 
VERBOSITY = 1;
nvox = numel(refimage);
[gx,gy,gz] = gradient(refimage);
gg = sqrt(gx.^2+gy.^2+gz.^2);
Mgg = max(gg(:));
if Mgg==0
    Mgg = 1; % to manage a flat mask
end
kappa = exp(- th*gg/Mgg); % factor 5 makes quite an extreme setting..
kap3D = spdiags(reshape(kappa,[],1),0:0,nvox,nvox);
%% plot
if VERBOSITY == 1
    figure(180),
    SubPlotMap(kappa,'Diffusivity Map',180,1,1,1,[4 4 1]);
    drawnow
end