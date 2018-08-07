%==========================================================================
%%                      RECONSTRUCTION DOMAIN: CW or TD
%==========================================================================
REC.domain = 'td';          % CW or TD: data type to be inverted
REC.type_fwd = 'linear';    % 'linear' or 'fem'.
% -------------------------------------------------------------------------
REC.time.roi = [];% ROI in time-step unit. If omitted, the ROI will be % REC.time.roi = [103 501;100 499;77 508;69 507;67 504;80 496;70 497;65 506];
                        % selected dinamically by the user.
NUM_TW = 20;            % Number of Time Windows within ROI
% =========================================================================
%%                        Initial parameter estimates 
% =========================================================================
% In this section all the parameter for the inverse solver are setted.
% --------------------------- Optical properties --------------------------
REC.solver.variables = {'mua','mus'}; % variables mua,mus.
REC.opt.mua0 = 0.01;    % absorption [mm-1]
REC.opt.musp0 = 1.0;      % reduced scattering [mm-1]
REC.opt.nB = 1.4;
if SPECTRA == 0
mua_ = [0.00380740474000000 0.00286718196800000 0.00353571764000000 0.0109674880900000 0.0170446766000000 0.0314137863400000 0.0185177919450000 0.0111878307750000];
musp_ = [1.54530000000000 1.46457537313433 1.18224759036145 1.07242131147541 1.04389946808511 1.00129132653061 0.952684951456311 0.921376056338028];
mua_ = mua_(lamda_id); musp_ = musp_(lamda_id);
Xr = {mua_,musp_};
else
a_ = 1.5453;	b_ = 1; % a (mm-1), b(adimensionale)
conc_ = [0.40855 0.71424 0.48844 0.5868	0.37051]; 
Xr = {conc_,[a_ b_]};
end

% ---------------------- Solver and regularization ------------------------
REC.solver.tau = 1e-4;            % regularisation parameter
REC.solver.type = 'spectral_born';         % 'born','GN': gauss-newton, 
                                  % 'USprior': Simon's strutural prior
                                  % 'LM': Levenberg-Marquardt,
                                  % 'l1': L1-based minimization
                                  % 'fit': fitting homogeneous data
if strcmpi(REC.solver.type,'spectral_born')&&SPECTRA == 0
MEx = MException('spectral_born:SpectralDataInput','Set SPECTRA = 1 to use spectra_born');
throwAsCaller(MEx);    
end
% =========================================================================
%%                            US prior 
% =========================================================================
REC.solver.prior.path = [];
% =========================================================================
%%                     load a precomputed jacobian 
% =========================================================================
% Pay attention! The jacobian depends on source-detectors configuration,
% optical properties of the background and number of time-windows.
REC.solver.prejacobian.load = false;
REC.solver.prejacobian.path = '../results/precomputed_jacobians/J';
