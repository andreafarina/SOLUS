function [contrast,  CNR] = recon_sensitivity(mu_recon, incl_idx,str_norm, vol_target, mu0)

% [contrast,  CNR] = recon_sensitivity(mu_recon, incl_idx,str_norm, vol_target)
% retrieves the sensisitivity figure of merit in the quantification of the inclusion,
% where the reconstructed value for the inclusion is assumed to be
% the average value inside the identified region of the inclusion
% mu_recon: reconstruction value
% incl_idx: 3D matrix which is one where the inclusion is recognised to be
% str_norm: string, if set to "volume", computes the results taking into
% account the volume mismatch
% vol_target: volume of the target in number of voxel
% mu0: reference value for which the jacobian was calculated





    mu_recon = mu_recon(:);
    if nargin < 2
        incl_idx = identify_inclusion(mu_recon(:));
    end
    if nargin < 3
       str_norm = 'else'; 
       mu0 = 0;
    end    
        bin_incl_idx = false(numel(mu_recon),1); 
        bin_incl_idx(incl_idx) = true; 
        bin_bulk_idx = logical(true(numel(mu_recon),1) - bin_incl_idx); 
        
        incl = mu_recon(bin_incl_idx);
        bulk = mu_recon(bin_bulk_idx);
        
        if strcmpi(str_norm,'volume') && nargin > 4
%             vol_norm = sum(incl_idx(:))/vol_target;
              vol_norm = sum(bin_incl_idx(:))/vol_target;
        else
            vol_norm = 1;
        end

%         diff_mean = mu0 + ( (mean(incl)- mu0)*vol_norm) - mean(bulk);
%         CNR = diff_mean/std(bulk); 
            diff_mean = (mu0 + (mean(incl)- mu0)*vol_norm) - mean(bulk);       
            contrast = diff_mean/ mean(bulk);
            rec = (mu0 + (mean(incl)- mu0)*vol_norm);
            back = mean(bulk);
            CNR = diff_mean/std(bulk);
end






