%==========================================================================
%%                      RECONSTRUCTION DOMAIN: CW or TD
%==========================================================================
REC.domain = 'td';          % CW or TD: data type to be inverted
REC.type_fwd = 'linear';    % 'linear' or 'fem'.
% -------------------------------------------------------------------------
REC.time.roi = [227 937;211 939;183 932;169 914;169 891;167 870;169 785;160 914];%[250 980;241 971;222 976;221 969;212 980;211 967;215 978;220 971];% ROI in time-step unit. If omitted, the ROI will be % REC.time.roi = [103 501;100 499;77 508;69 507;67 504;80 496;70 497;65 506];
%[257 600;234 653;209 649;199 639;190 584;185 469;179 582;188 676];%[250 980;241 971;222 976;221 969;212 980;211 967;215 978;220 971];% ROI in time-step unit. If omitted, the ROI will be % REC.time.roi = [103 501;100 499;77 508;69 507;67 504;80 496;70 497;65 506];
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
mua_ = [0.00380793160000000,0.00134266120000000,0.00138762010000000,0.0100602733000000,0.00755439960000000,0.00394251390000000,0.00761375700000000,0.00493305200000000];
musp_ = [1.22842091900000,1.18990461900000,0.913511025000000,0.803037154800000,0.769676676200000,0.756078261900000,0.702287059500000,0.675179692900000];
mua_ = mua_(lamda_id); musp_ = musp_(lamda_id);
Xr = {mua_,musp_};
else
a_ = 1.5221;	b_ = 1.1415;
conc_ = [1.4492 0.48527 0.97134 0.037411 0.2021];
Xr = {conc_,[a_ b_]};
end

% ---------------------- Solver and regularization ------------------------
REC.solver.tau = 1e-2;            % regularisation parameter
REC.solver.type = 'spectral_born';         % 'born','GN': gauss-newton, 
                                  % 'USprior': Simon's strutural prior
                                  % 'LM': Levenberg-Marquardt,
                                  % 'l1': L1-based minimization
                                  % 'fit': fitting homogeneous data
                                  % 'spectral_fit': fitting homogeneous
                                  % data with spectral model
                                  % 'spectral_usprior': spectral approach
                                  % to USprior
                                  % 'spectral_born': recon with spectral
                                  % Jacobian
                                  % 'born_spectral_post_proc': multi_wave
                                  % classical born with post proc
                                  % chromophores estimation
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
