function [contrast,  CNR] = recon_sensitivity(mu_recon, incl_idx,str_norm, vol_target)

% [contrast,  CNR] = recon_sensitivity(mu_recon, incl_idx,str_norm, vol_target)
% retrieves the sensisitivity figure of merit in the quantification of the inclusion,
% where the reconstructed value for the inclusion is assumed to be
% the average value inside the identified region of the inclusion
% mu_recon: reconstruction value
% incl_idx: 3D matrix which is one where the inclusion is recognised to be
% str_norm: string, if set to "volume", computes the results taking into
% account the volume mismatch
% vol_target: volume of the target in number of voxel




    mu_recon = mu_recon(:);
    if nargin < 2
        incl_idx = identify_inclusion(mu_recon(:));
    end
        bin_incl_idx = false(numel(mu_recon),1);
        bin_incl_idx(incl_idx) = true;
        bin_bulk_idx = logical(true(numel(mu_recon),1) - bin_incl_idx);
        
        incl = mu_recon(bin_incl_idx);
        bulk = mu_recon(bin_bulk_idx);
        
        if strcompi(str_norm,'volume') && nargin > 3
            vol_norm = sum(incl_idx(:))/vol_target;
        else
            vol_norm = 1;
        end

        diff_mean = mean(incl)*vol_norm - mean(bulk);
        
        contrast = diff_mean/ mean(bulk);
        CNR = diff_mean/std(bulk);

end