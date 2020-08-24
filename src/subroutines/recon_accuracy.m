function err = recon_accuracy(mu_recon, target, idx, str_norm, vol_target)

% err = recon_accuracy(mu_recon, target, idx)
% retrieves the relative error in the quantification of the inclusion,
% where the reconstructed value for the inclusion is assumed to be
% the average value inside the identified region of the inclusion
% mu_recon: reconstruction value
% target: 3D matrix with target of the reconstruction
% idx: 3D matrix which is one where the inclusion is recognised to be
% str_norm: string, if set to "volume", computes the results taking into
% account the volume mismatch
% vol_target: volume of the target in number of voxel

    if nargin < 3
        idx = identify_inclusion(mu_recon(:));
    end
    
    if strcompi(str_norm,'volume') && nargin > 4
        vol_norm = sum(idx(:))/vol_target;
        
    else
        vol_norm = 1;
    end
    
    err = mean(mu_recon(idx))*vol_norm - target;
    err = err/target;
    return
end