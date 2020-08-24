function [displacement, broadening, accuracy] = Quantify_SOLUS(recon, true, vol_true, grid)
% computes figure of merit described by SOLUS for single lambda and single coefficient


    idx = find_recon_inclusion(recon);
    
    % displacement
    vol_recon = zeros(grid.dim);
    vol_recon(idx) = 1;
    [xx,yy,zz] = meshgrid(grid.x,grid.y,grid.z);
    centre_true = zeros(3,1);
    centre_recon = zeros(3,1);
    centre_true = [sum(vol_true(:).*xx(:));sum(vol_true(:).*yy(:));sum(vol_true(:).*zz(:))];
    centre_recon = [sum(vol_recon(:).*xx(:));sum(vol_recon(:).*yy(:));sum(vol_recon(:).*zz(:))];
    displacement = centre_recon - centre_true;
    
    % broadening
    sigma_true = zeros(3,1);
    sigma_recon = zeros(3,1);
    sigma_true = [max(squeeze(sum(vol_true,1)),[1,2]),max(squeeze(sum(vol_true,2)),[1,2]),max(squeeze(sum(vol_true,3)),[1,2])];
    sigma_recon = [max(squeeze(sum(vol_recon,1)),[1,2]),max(squeeze(sum(vol_recon,2)),[1,2]),max(squeeze(sum(vol_recon,3)),[1,2])];
    broadening = sigma_recon - sigma_true .*[grid.dx;grid.dy;grid.dz];
    
    % accuracy
    pert_recon = mean(recon(idx));
    pert_true = mean (true(vol_true(:)));
    accuracy = (pert_recon - pert_true) / pert_true; 


end