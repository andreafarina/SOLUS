function idx = identify_inclusion(mu_recon)
%
%  idx = identify_inclusion(mu_recon)
%
% Retrieves the indexes idx where the inclusion is reconstructed. 
% It operates a normalisation by median and standard deviation
% Selects the indexes of the input values that are above 4 or below -4 standard deviations
% If no inclusion is found with this condition, it returns all the indexes
% which correspond to a nonzero value.

    mu_recon = mu_recon(:);
    prcmax = 98; 
    prcmin = 2;
    mu_recon(mu_recon > prctile(mu_recon, prcmax)) = prctile(mu_recon, prcmax);
    mu_recon(mu_recon < prctile(mu_recon, prcmin)) = prctile(mu_recon, prcmin);
    mu_recon = mu_recon - median(mu_recon(:));
    mu_recon = mu_recon/(6*std(mu_recon(:)));
    idx = find(mu_recon > 1 | mu_recon < - 1);
    if isempty(idx)
        idx= find(mu_recon > 0 |mu_recon < 0);
    end

end

