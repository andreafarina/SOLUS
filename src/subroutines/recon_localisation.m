% function [displ, brd, brdRMS] = recon_localisation(mu_recon, bin_target, grid,dim, idx)
% 
% % displ, brd, brdRMS] = recon_localisation(mu_recon, bin_target, grid,dim, idx)
% % Geometrical localisation of inclusion.
% % mu_recon: reconstruction
% % bin_target: binary characteristic function which is 1 where the inclusion
% % ground truth is and zeros otherwise
% % grid: grid as defined ins the SOLUS software
% % dim: dimensions of the domain
% % idx: reconstructed characteristic function
% % displ: displacement vector of the centre of the inclusion in x,y,z
% % brd: broadening of the inclusion along x,y,z considering the maximum extesion
% % of the considered inclusion along each axis as its characteristic dimension
% % brdRMS: broadening of the inclusion along x,y,z considering the RMS extesion
% % of the considered inclusion along each axis as its characteristic dimension
% 
%   
%     if nargin < 6
%         idx = identify_inclusion(mu_recon);
%     end
% 
%     if nargin < 5
%         dim = size(mu_recon);
%     end
%     
%     if numel(dim) ~= 3
%         disp('Wrong dimensions');
%         return;
%     end
%     
%     [xx, yy, zz] = ndgrid(grid.x,grid.y,grid.z);
%     
%     % displacement
%     xc = mean(xx(idx));
%     yc = mean(yy(idx));
%     zc = mean(zz(idx));
%     
%     xc_target = mean(xx(bin_target(:)));
%     yc_target = mean(yy(bin_target(:)));
%     zc_target = mean(zz(bin_target(:)));
%     
%     displ = [xc - xc_target;yc - yc_target;zc - zc_target;];
%     
%     % broadening
%     brd = zeros(3,1);
%     brdRMS = zeros(3,1);
%     bin_recon = zeros(size(mu_recon));
%     bin_recon(idx) = 1;
%    
%     dL = [grid.dx;grid.dy;grid.dz];
%     for dd = 1:numel(dim)
%         dl = dL(dd);
%         cumul_sq = squeeze(dl*sum(bin_recon, dd)).^2;
%         RMS = sqrt( mean( cumul_sq(cumul_sq > 0 ) ));
%         
%         cumul_sq_target = squeeze(dl*sum(bin_target, dd)).^2;
%         RMS_target = sqrt( mean( cumul_sq_target(cumul_sq_target > 0 ) ));
%         
%         brdRMS(dd) = RMS - RMS_target;
%         brd(dd) = max(sqrt(cumul_sq(:))) - max(sqrt(cumul_sq_target(:)));
%     end
%     
%     return;
% end

function [displ, brd, brdRMS] = recon_localisation(mu_recon, bin_target, grid,dim, idx, manual, coeff)

% [displ, brd, brdRMS] = recon_localisation(mu_recon, bin_target, grid,dim, idx)
% Geometrical localisation of inclusion.
% mu_recon: reconstruction
% bin_target: binary characteristic function which is 1 where the inclusion
% ground truth is and zeros otherwise
% grid: grid as defined ins the SOLUS software
% dim: dimensions of the domain
% idx: reconstructed characteristic function
% displ: displacement vector of the centre of the inclusion in x,y,z
% brd: broadening of the inclusion along x,y,z considering the maximum extesion
% of the considered inclusion along each axis as its characteristic dimension
% brdRMS: broadening of the inclusion along x,y,z considering the RMS extesion
% of the considered inclusion along each axis as its characteristic dimension

  
    if nargin < 6
        idx = identify_inclusion(mu_recon);
    end

    if nargin < 5
        dim = size(mu_recon);
    end
    
    if numel(dim) ~= 3
        disp('Wrong dimensions');
        return;
    end
    
    [xx, yy, zz] = ndgrid(grid.x,grid.y,grid.z);
    
    % displacement
    xc = mean(xx(idx));
    yc = mean(yy(idx));
    zc = mean(zz(idx));
    
    xc_target = mean(xx(bin_target(:))); 
    yc_target = mean(yy(bin_target(:)));
    zc_target = mean(zz(bin_target(:)));
    
    displ = [xc - xc_target;yc - yc_target;zc - zc_target;];
    
    % broadening
    brd = zeros(3,1);
    brdRMS = zeros(3,1);
    bin_recon = zeros(size(mu_recon));
    bin_recon(idx) = 1;
   
    dL = [grid.dx;grid.dy;grid.dz];
    for dd = 1:numel(dim)
        dl = dL(dd);
        cumul_sq = squeeze(dl*sum(bin_recon, dd)).^2;
        RMS = sqrt( mean( cumul_sq(cumul_sq > 0 ) ));
        
        cumul_sq_target = squeeze(dl*sum(bin_target, dd)).^2;
        RMS_target = sqrt( mean( cumul_sq_target(cumul_sq_target > 0 ) )); %mettilo a 0
        if exist ('coeff')
            if coeff == 'a' || coeff == 'b'
                RMS_target = 0;
            end
        end
        brdRMS(dd) = RMS - RMS_target;
        brd(dd) = max(sqrt(cumul_sq(:))) - max(sqrt(cumul_sq_target(:)));
    end
    
    return;
end

