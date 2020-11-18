%==========================================================================
%%                      RECONSTRUCTION DOMAIN: CW or TD
%==========================================================================
REC.domain = 'td';          % CW or TD: data type to be inverted
REC.type_fwd = 'linear';    % 'linear' or 'fem'.
% -------------------------------------------------------------------------
% REC.time.roi = [];
REC.time.roi = [259 919];%;243 790;220 661;225 446;224 420;225 462;226 468;236 538];
                        % ROI in time-step unit. 
                        % If omitted, the ROI will be
                        % selected dinamically by the user.
NUM_TW = 20;            % Number of Time Windows within ROI
% =========================================================================
%%                        Initial parameter estimates 
% =========================================================================
% In this section all the parameter for the inverse solver are setted.
% --------------------------- Optical properties --------------------------
REC.solver.variables = {'mua'}; % variables mua,mus.
REC.opt.mua0 = 0.01;    % absorption [mm-1]
REC.opt.musp0 = 1.0;      % reduced scattering [mm-1]
REC.opt.nB = 1.4;
%if SPECTRA == 0
    mua_ = 0.01;
    musp_ = 1.0;
    Xr = {mua_,musp_};
% else
%     a_ = 1.234; b_ = 0.50;
%     conc_ = [1.45 9.38 0.8186*10*100*0.91 0.1172*10*100 0.0453*10*100*0.196];
%     Xr = {conc_,[a_ b_]};
% end

% ---------------------- Solver and regularization ------------------------
REC.solver.tau = 0.005;%logspace(-5,0,30);            % regularisation parameter
REC.solver.type = 'tk1';          % 'tk0','GN': gauss-newton, 
                                  % 'tk1': first order Tichonov regul
                                  % 'USprior': tk1+structural prior
                                  % 'LM': Levenberg-Marquardt,
                                  % 'l1': L1-based minimization
                                  % 'fit': fitting homogeneous data
                                  % 'spectral_fit': fitting homogeneous
                                  % data with spectral model
                                  % 'spectral_usprior': spectral approach
                                  % to USprior
                                  % 'spectral_born': recon with spectral
                                  % Jacobian
                                  % 'spectral_tk1': first order Tichonov
                                  % regul for spectral Jacobian
                                  % 'born_spectral_post_proc': multi_wave
                                  % classical born with post proc
                                  % chromophores estimation
if (strcmpi(REC.solver.type,'spectral_born')||strcmpi(REC.solver.type,'spectral_usprior'))&&SPECTRA == 0
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
