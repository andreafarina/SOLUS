%==========================================================================
%%                      RECONSTRUCTION DOMAIN: CW or TD
%==========================================================================
REC.domain = 'td';          % CW or TD: data type to be inverted

% -------------------------------------------------------------------------
NUM_TW = 6;             % Number of Time Windows within ROI
% =========================================================================
%%                        Initial parameter estimates 
% =========================================================================
% In this section all the parameter for the inverse solver are setted.
% --------------------------- Optical properties --------------------------
REC.opt.mua0 = 0.01;    % absorption [mm-1]
REC.opt.musp0 = 1;    % reduced scattering [mm-1]
REC.opt.nB = 1.4;
% ---------------------- Solver and regularization ------------------------
REC.solver.tau = 1;            % regularisation parameter
REC.solver.type = 'Usprior';         % 'born','GN': gauss-newton, 
                                  % 'USprior': Simon's strutural prior
                                  % 'LM': Levenberg-Marquardt,
                                  % 'l1': L1-based minimization
                                  % 'fit': fitting homogeneous data

REC.solver.prior = [];