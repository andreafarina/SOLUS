function err = recon_accuracy(mu_recon, target, idx, str_norm, vol_target, mu0)

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
% mu0: reference value for which the jacobian was calculated
    if nargin < 3
        idx = identify_inclusion(mu_recon(:));
    end
    if nargin < 4
       str_norm = 'else'; 
       mu0 = 0;
    end

    
    if strcmpi(str_norm,'volume') && nargin > 5
        vol_norm = sum(idx(:))/vol_target;
        
    else
        vol_norm = 1;
    end
    
    err = mu0 + ( mean(mu_recon(idx)-mu0)*vol_norm) - target;
    err = err/target;
    return
end