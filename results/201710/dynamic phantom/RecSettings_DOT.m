%==========================================================================
%%                      RECONSTRUCTION DOMAIN: CW or TD
%==========================================================================
REC.domain = 'td';          % CW or TD: data type to be inverted

% -------------------------------------------------------------------------
REC.time.roi = round([1021,2500]/4);
NUM_TW = 10;%diff(REC.time.roi)+1;%10;%diff(REC.time.roi)+1;             % Number of Time Windows within REC.time.roi
% =========================================================================
%%                        Initial parameter estimates 
% =========================================================================
% In this section all the parameter for the inverse solver are setted.
% --------------------------- Optical properties --------------------------
REC.opt.mua0 = 0.012;    % absorption [mm-1]
REC.opt.musp0 = 1.2;      % reduced scattering [mm-1]
REC.opt.nB = 1.5;
% ---------------------- Solver and regularization ------------------------
REC.solver.tau = 1e-3;            % regularisation parameter
REC.solver.type = 'born';         % 'born','GN': gauss-newton, 
                                  % 'USprior': Simon's strutural prior
                                  % 'LM': Levenberg-Marquardt,
                                  % 'l1': L1-based minimization
                                  % 'fit': fitting homogeneous data
% =========================================================================
%%                            US prior 
% =========================================================================
REC.solver.prior = [];
% =========================================================================
%%                     load a precomputed jacobian 
% =========================================================================
% Pay attention! The jacobian depends on source-detectors configuration,
% optical properties of the background and number of time-windows.
REC.solver.prejacobian.load = true;
REC.solver.prejacobian.path = [getenv('DOTDIR'),filesep,'/precomputed_jacobians/J'];
