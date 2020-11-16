function [dphi,sd] = PrepareDataFitting(data,ref,sd,type_ratio,type_ref,proj)
% Prepare data fitting term
% data, ref, sd derived by measurements
% proj is the theoretical model at the bkg optical properties
%
% type_ratio:   'Born' :  (data-ref)/ref
%               'Rytov':  log(data/ref)
%               'gauss':  (data - ref)/sd    
%
% type_ref:     'meas' :  use measured ref
%               'theor':  rescale with model
%               'area' :  normalize to ref area

switch lower(type_ref)
    case 'meas'
        ref = ref(:);
        factor = 1;
    case 'theor'
        factor = proj(:)./ref(:);
    case 'area'
        nwin = size(ref,1);
        Aref = sum(ref);
        Aproj = sum(proj);
        factor = repmat(Aproj./Aref,[nwin,1]);
        factor = factor(:);
        
    otherwise
        warning ('Unknown reference type');
        return;
end

%% apply factors
data = data(:) .* factor;
ref = ref(:) .* factor;
sd = sd(:) .* factor;
        
switch lower(type_ratio)
    case 'born'
        dphi = (data - ref)./ref;
        sd = ones(size(dphi));
    case 'rytov'
        dphi = log(data) - log(ref);
        sd = ones(size(dphi));
    case 'gauss'
        dphi = (data - ref)./sd;
        if strcmpi(type_ref,'meas')
            warning('Attention! Data are not normalized to the model!');
        end
        
end
        
        