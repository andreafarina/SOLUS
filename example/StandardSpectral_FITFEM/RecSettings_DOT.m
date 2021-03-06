%==========================================================================
%%                      RECONSTRUCTION DOMAIN: CW or TD
%==========================================================================
REC.domain = 'td';          % CW or TD: data type to be inverted
REC.type_fwd = 'fem';    % 'linear' or 'fem'.
% -------------------------------------------------------------------------
REC.time.roi =[];%[32   105;    32   106;    25    98;    24    99;    25   103;    23   103;    26   102;   19   111;];
    %[227 937;211 939;183 932;169 914;169 891;167 870;169 785;160 914];
                        % ROI in time-step unit. 
                        % If omitted, the ROI will be
                        % selected dinamically by the user.
NUM_TW = 50;            % Number of Time Windows within ROI
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
mua_ = [0.003807,0.001342,0.001387,0.010060,0.007554,0.003942,0.007613,0.004933];
musp_ = [1.22842,1.189904,0.913511,0.803037,0.769676,0.7560782,0.702287,0.675179];
Xr = {mua_,musp_};
else
a_ = 0.5*1.2883;	b_ = 0.5*1.2641;
conc_ =  [0.5 0.5 0.5 0.5 0.5].*[1.91 8.8 575 270 108];
Xr = {conc_,[a_ b_]};
end

% ---------------------- Solver and regularization ------------------------
REC.solver.tau = 0.01;            % regularisation parameter
REC.solver.type = 'components_fit';         % 'born','GN': gauss-newton, 
                                  % 'USprior': Simon's strutural prior
                                  % 'LM': Levenberg-Marquardt,
                                  % 'l1': L1-based minimization
                                  % 'fit': fitting homogeneous data
                                  % 'spectral_fit': fitting homogeneous
                                  % data with spectral model
                                  % 'spectral_usprior': spectral approach
                                  % to USprior
                                  % 'spectral_fitFEM':  
                                  % 'spectral_born': recon with spectral
                                  % Jacobian
                                  % 'components_fit': sperctral fit 
                                  % 'born_spectral_post_proc': multi_wave
                                  % classical born with post 
                                  % chromophores estimation
if (strcmpi(REC.solver.type,'spectral_born')||strcmpi(REC.solver.type,'spectral_usprior'))&&SPECTRA == 0
MEx = MException('spectral_born:SpectralDataInput','Set SPECTRA = 1 to use spectra_born');
throwAsCaller(MEx);    
end
% =========================================================================
%%                            US prior 
% =========================================================================
REC.solver.prior.path = 'DTkwaved.mat';%'ham_maskYX_flip.mat';
% =========================================================================
%%                     load a precomputed jacobian 
% =========================================================================
% Pay attention! The jacobian depends on source-detectors configuration,
% optical properties of the background and number of time-windows.
REC.solver.prejacobian.load = false;
REC.solver.prejacobian.path = '../results/precomputed_jacobians/J';
