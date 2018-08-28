function [data_o,sd_o,ref_o,sdj,proj,mask] = PrepDataLH(data,sd,ref,proj,kind,...
    opt_calib, opt_sd)
% Prepare data for Likelihood
% for TD, data, ref, sd, proj have the dimension of size(twin,1)
% opt_data:
%   'none':
%
nwin = size(proj,1);
switch lower(kind)
    case 'nonlin'
        
        %% 1 - factor
        switch lower(opt_calib)
            case 'factor_td'
                factor = proj./ref;
            case 'factor_area'
                factor = repmat(sum(proj),nwin,1)./repmat(sum(ref),nwin,1);
            case 'factor_norm'
                factor = norm(proj(:))./norm(ref(:));
            case 'factor_first'
                factor = norm(proj(:,1))./norm(ref(:,1));
            otherwise
                factor = 1;
        end
        %
        data_o = data.*factor;
        ref_o = ref.*factor;
        %% 2 - data
        switch lower(opt_sd)
            case 'poisson'
                sd_o = sd.*factor;
            case 'proj'
                sd_o = proj;
            case 'none'
                sd_o = ones(size(proj));
        end
        sdj = ones(size(sd_o));
        %ref_o = ref;
        obj = (data_o-proj)./sd_o;
        mask = isinf(obj)|isnan(obj);
             
    case 'born'
        data_o = data;
        ref_o = ref;
        sd_o = ref;
        sdj = proj;
        obj = (data-ref)./ref;
        mask = isinf(obj)|isnan(obj)|(proj==0);
    case 'rytov'
        data_o = data./ref;
        proj = 0;
        ref_o = ref;
        sdj = proj;
        mask = isinf(data_o)||isnan(data_0);
end
%% apply mask
% data_o(mask) = [];
% ref_o(mask) = [];
% sd_o(mask) = [];
% sdj(mask) = [];
% proj(mask) = [];
%% vectorize
data_o = data_o(:);
ref_o = ref_o(:);
sd_o = sd_o(:);
sdj = sdj(:);
mask = mask(:);
proj = proj(:);


