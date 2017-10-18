if exist('isComingFromInterface','var')
    warning('off','backtrace');
    % Inverse problem: parameters for the inverse solver
    if exist('TW','var'), NUM_TW = eval('P(TW).Value'); else, warning('Parameter P(TW).Value not in override. Default value will be used'); end
    if exist('MUA0','var'), REC.opt.mua0 = eval('P(MUA0).Value'); else, warning('Parameter P(MUA0).Value not in override. Default value will be used'); end
    if exist('MUS0.Value','var'), REC.opt.musp0 = eval('P(MUS0).Value'); else, warning('Parameter P(MUS0).Value not in override. Default value will be used'); end
    if exist('NB','var'), REC.opt.nB = eval('P(NB).Value'); else, warning('Parameter P(NB).Value not in override. Default value will be used'); end
    
    % Regularization Parameters
    if exist('TAU','var'), REC.solver.tau = eval('P(TAU).Value'); else, warning('Parameter P(TAU).Value not in override. Default value will be used'); end
    if exist('SOLVTYPE','var'), REC.solver.type = eval('P(SOLVTYPE).Value'); else, warning('Parameter P(SOLVTYPE).Value not in override. Default value will be used'); end                                                                                             % 'born','GN': gauss-newton,
    % 'USprior': Simon's strutural prior
    % 'LM': Levenberg-Marquardt,
    % 'l1': L1-based minimization
    % 'fit': fitting homogeneous data
    % Forward
    if exist('EXPDELTA','var'), EXP_DELTA = eval('P(EXPDELTA).Value'); else, warning('Parameter P(EXPDELTA).Value not in override. Default value will be used'); end
    
    % Background optical properties
    if exist('MUAB','var'), DOT.opt.muaB = eval('P(MUAB).Value'); else, warning('Parameter P(MUAB).Value not in override. Default value will be used'); end
    if exist('MUSB','var'), DOT.opt.muspB = eval('P(MUSB).Value'); else, warning('Parameter P(MUSB).Value not in override. Default value will be used'); end                      % mm-1
    if exist('INRI','var'), DOT.opt.nB = eval('P(INRI).Value'); else, warning('Parameter P(INRI).Value not in override. Default value will be used'); end                     % mm-1
    if exist('EXRI','var'), DOT.opt.nE = eval('P(EXRI).Value'); else, warning('Parameter P(EXRI).Value not in override. Default value will be used'); end                      % mm-1
    
    % Grid
    if exist('X1','var'), DOT.grid.x1 = eval('P(X1).Value'); else, warning('Parameter P(X1).Value not in override. Default value will be used'); end           %0
    if exist('X2','var'), DOT.grid.x2 = eval('P(X2).Value'); else, warning('Parameter P(X2).Value not in override. Default value will be used'); end           %64
    if exist('DX','var'), DOT.grid.dx = eval('P(DX).Value'); else, warning('Parameter P(DX).Value not in override. Default value will be used'); end           %2
    if exist('Y1','var'), DOT.grid.y1 = eval('P(Y1).Value'); else, warning('Parameter P(Y1).Value not in override. Default value will be used'); end           %0
    if exist('Y2','var'), DOT.grid.y2 = eval('P(Y2).Value'); else, warning('Parameter P(Y2).Value not in override. Default value will be used'); end           %58
    if exist('Z1','var'), DOT.grid.z1 = eval('P(Z1).Value'); else, warning('Parameter P(Z1).Value not in override. Default value will be used'); end           %0
    if exist('Z2','var'), DOT.grid.z2 = eval('P(Z2).Value'); else, warning('Parameter P(Z2).Value not in override. Default value will be used'); end           %32
    
    % Inclusion 1
    if exist('XP','var'), DOT.opt.hete1.c = eval('[P(XP).Value, P(YP).Value, P(ZP).Value]'); else, warning('Parameter P(XP).Value not in override. Default value will be used'); end           % down
    if exist('MUAP','var'), DOT.opt.hete1.val = eval('P(MUAP).Value')+DOT.opt.muaB; else, warning('Parameter P(MUAP).Value not in override. Default value will be used'); end
    if exist('SIG','var'), DOT.opt.hete1.sigma = eval('P(SIG).Value'); else, warning('Parameter P(SIG).Value not in override. Default value will be used'); end
    
    % time domain parameters
    if exist('DT','var'), DOT.time.dt = eval('P(DT).Value'); else, warning('Parameter P(DT).Value not in override. Default value will be used'); end           % time step in picoseconds
    if exist('NT','var'), DOT.time.nstep = eval('P(NT).Value'); else, warning('Parameter P(NT).Value not in override. Default value will be used'); end           % number of temporal steps
    if exist('CR','var'), DOT.time.TotCounts = eval('P(CR).Value'); else, warning('Parameter P(CR).Value not in override. Default value will be used'); end           % total counts for the maximum-energy
    
    % Radiometry
    if exist('RAD','var'), RADIOMETRY = eval('P(RAD).Value'); else, warning('Parameter P(RAD).Value not in override. Default value will be used'); end           % Apply radiometry
    if exist('DT','var'), DOT.radiometry.timebin = eval('P(DT).Value'); else, warning('Parameter P(DT).Value not in override. Default value will be used'); end           % time step in picoseconds
    if exist('PW','var'), DOT.radiometry.power = eval('P(PW).Value'); else, warning('Parameter P(PW).Value not in override. Default value will be used'); end           % (mW) injected laser power
    if exist('TA','var'), DOT.radiometry.acqtime = eval('P(TA).Value'); else, warning('Parameter P(TA).Value not in override. Default value will be used'); end           % (s) acquisition time
    if exist('OE','var'), DOT.radiometry.opteff = eval('P(OE).Value'); else, warning('Parameter P(OE).Value not in override. Default value will be used'); end           % typical efficiency of the optical path
    if exist('EA','var'), DOT.radiometry.area = eval('P(EA).Value'); else, warning('Parameter P(EA).Value not in override. Default value will be used'); end           % effective area of the detector
    if exist('QE','var'), DOT.radiometry.qeff = eval('P(QE).Value'); else, warning('Parameter P(QE).Value not in override. Default value will be used'); end           % quantum efficiency
    
    % permutation matrix
    %if exist('SD','var'), DOT.dmask = eval('P(SD).Value'); else, warning('Parameter P(SD).Value not in override. Default value will be used'); end
    
    % Type of Cutting for statistics
    if exist('CUT','var'), CUT_COUNTS = eval('P(CUT).Value'); else, warning('Parameter P(CUT).Value not in override. Default value will be used'); end % 0 non gated, 1 gated
    if exist('ND','var'), NumDelays = eval('P(ND).Value'); else, warning('Parameter P(ND).Value not in override. Default value will be used'); end % number of delays
    %if exist('NG','var'), NumGates = eval('P(NG).Value'); else, warning('Parameter P(NG).Value not in override. Default value will be used'); end % number of gates
    
    % Load Jacobian
    if exist('LJ','var')
        if iL == 1
            REC.solver.prejacobian.load = false;
            BuffMuaB = REC.opt.mua0;
            BuffMusB = REC.opt.musp0;
        else
            if (BuffMuaB == REC.opt.mua0 && BuffMusB == REC.opt.musp0)==false
                BuffMuaB = REC.opt.mua0;
                BuffMusB = REC.opt.musp0;
                REC.solver.prejacobian.load = false;
            else
                REC.solver.prejacobian.load = eval('P(LJ).Value');
            end
        end
    else
        warning('Parameter P(LJ).Value not in override. Default value will be used');
    end
    %if exist('LJ','var'), REC.solver.loadjacobian = eval('P(LJ).Value'); else, warning('Parameter P(LJ).Value not in override. Default value will be used'); end
    warning('on','backtrace');
end