%==========================================================================
%%                      RECONSTRUCTION DOMAIN: CW or TD
%==========================================================================
REC.domain = 'td';          % CW or TD: data type to be inverted
REC.type_fwd = 'linear';    % 'linear' or 'fem'.
% -------------------------------------------------------------------------
REC.time.roi = [282 946;264 941;241 942;232 946;222 946;218 942;213 941;218 944];

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
mua_ = [0.09118862 0.08252885 0.07457543 0.114975675 0.07321214 0.07269354 0.09152519 0.076381]/10;
musp_ = [14.267035 13.06632 9.491283 7.8722685 7.707295 7.060895 6.352237 5.989931]/10;
mua_ = mua_(lamda_id); musp_ = musp_(lamda_id);
Xr = {mua_,musp_};
else
a_ = 0.789550;	b_ =1.05756;
conc_ = [9.50498	0.92704	0.15076	0.63606	1.11337];
Xr = {conc_,[a_ b_]};
end

% ---------------------- Solver and regularization ------------------------
REC.solver.tau = 0.1;            % regularisation parameter
REC.solver.type = 'born';         % 'born','GN': gauss-newton, 
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
REC.solver.prior.path = []%['DTsilicon_mask.mat'];%;['ham_maskYX_flip.mat'];
% =========================================================================
%%                     load a precomputed jacobian 
% =========================================================================
% Pay attention! The jacobian depends on source-detectors configuration,
% optical properties of the background and number of time-windows.
REC.solver.prejacobian.load = false;
REC.solver.prejacobian.path = 'precomputed_jacobians/J';
