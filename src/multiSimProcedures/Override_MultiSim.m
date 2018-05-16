if exist('isComingFromInterface','var')
    warning('off','backtrace');
    % Inverse problem: parameters for the inverse solver
    if exist('ROISTART','var'), REC.time.roi(1)= eval('P(ROISTART).Value'); else, warning('Parameter P(ROISTART).Value not in override. Default value will be used'); end
    if exist('ROISTOP','var'), REC.time.roi(2)=eval('P(ROISTOP).Value'); else, warning('Parameter P(ROISTOP).Value not in override. Default value will be used'); end
    if exist('TW','var'), NUM_TW = eval('P(TW).Value'); else, warning('Parameter P(TW).Value not in override. Default value will be used'); end
    if exist('MUA0','var'), REC.opt.mua0 = eval('P(MUA0).Value'); else, warning('Parameter P(MUA0).Value not in override. Default value will be used'); end
    if exist('MUS0','var'), REC.opt.musp0 = eval('P(MUS0).Value'); else, warning('Parameter P(MUS0).Value not in override. Default value will be used'); end
    if exist('NB','var'), REC.opt.nB = eval('P(NB).Value'); else, warning('Parameter P(NB).Value not in override. Default value will be used'); end
    
    % Regularization Parameters
    if exist('VARS','var'), REC.solver.variables = eval('P(VARS).Value'); else, warning('Parameter P(VARS).Value not in override. Default value will be used'); end
    if exist('TAU','var'), REC.solver.tau = eval('P(TAU).Value'); else, warning('Parameter P(TAU).Value not in override. Default value will be used'); end
    if exist('SOLVTYPE','var'), REC.solver.type = eval('P(SOLVTYPE).Value'); else, warning('Parameter P(SOLVTYPE).Value not in override. Default value will be used'); end                                                                                             % 'born','GN': gauss-newton,
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
    
    if exist('NUMBER_HETE','var'), NUM_HETE = eval('P(NUMBER_HETE).Value'); else, warning('Parameter P(NUMBER_HETE).Value not in override. Default value will be used'); end           % NUM_HETE
    % Inclusion 1
    if exist('INCTYPE1','var'), DOT.opt.hete1.type = eval('P(INCTYPE1).Value'); else, warning('Parameter P(INCTYPE1).Value not in override. Default value will be used'); end           % type of the inclusion
    if exist('XP1','var'), DOT.opt.hete1.c = eval('[P(XP1).Value, P(YP1).Value, P(ZP1).Value]'); else, warning('Parameter P(XP1).Value not in override. Default value will be used'); end           % down
    if exist('SIG1','var'), DOT.opt.hete1.sigma = eval('P(SIG1).Value'); else, warning('Parameter P(SIG1).Value not in override. Default value will be used'); end
    if exist('PROFILE1','var'), DOT.opt.hete1.profile = eval('P(PROFILE1).Value'); else, warning('Parameter P(PROFILE1).Value not in override. Default value will be used'); end %profile of the inclusion
    if exist('INCPEAKVAL1','var'), DOT.opt.hete1.val = eval('P(INCPEAKVAL1).Value'); else, warning('Parameter P(INCPEAKVAL1).Value not in override. Default value will be used'); end
    % Inclusion 2
    if NUM_HETE > 1
        if exist('INCTYPE2','var'), DOT.opt.hete2.type = eval('P(INCTYPE2).Value'); else, warning('Parameter P(INCTYPE2).Value not in override. Default value will be used'); end           % type of the inclusion
        if exist('XP2','var'), DOT.opt.hete2.c = eval('[P(XP1).Value, P(YP1).Value, P(ZP1).Value]'); else, warning('Parameter P(XP2).Value not in override. Default value will be used'); end           % down
        if exist('SIG2','var'), DOT.opt.hete2.sigma = eval('P(SIG2).Value'); else, warning('Parameter P(SIG).Value not in override. Default value will be used'); end
        if exist('PROFILE2','var'), DOT.opt.hete2.profile = eval('P(PROFILE2).Value'); else, warning('Parameter P(PROFILE2).Value not in override. Default value will be used'); end %profile of the inclusion
        if exist('INCPEAKVAL2','var'), DOT.opt.hete2.val = eval('P(INCPEAKVAL2).Value'); else, warning('Parameter P(INCPEAKVAL2).Value not in override. Default value will be used'); end
    end
    
    % time domain parameters
    if exist('DT','var'), DOT.time.dt = eval('P(DT).Value'); else, warning('Parameter P(DT).Value not in override. Default value will be used'); end           % time step in picoseconds
    if exist('NT','var'), DOT.time.nstep = eval('P(NT).Value'); else, warning('Parameter P(NT).Value not in override. Default value will be used'); end           % number of temporal steps
    if exist('NOISETYPE','var'), DOT.time.noise = eval('P(NOISETYPE).Value'); else, warning('Parameter P(NOISETYPE).Value not in override. Default value will be used'); end
    if exist('SELF_NORM','var'), DOT.time.self_norm = eval('P(SELF_NORM).Value'); else, warning('Parameter P(SELF_NORM).Value not in override. Default value will be used'); end           % total counts for the maximum-energy
    if exist('CR','var'), DOT.time.TotCounts = eval('P(CR).Value'); else, warning('Parameter P(CR).Value not in override. Default value will be used'); end           % total counts for the maximum-energy
    
    % Radiometry
    if exist('RAD','var'), RADIOMETRY = eval('P(RAD).Value'); else, warning('Parameter P(RAD).Value not in override. Default value will be used'); end           % Apply radiometry
    if exist('PW','var'), DOT.radiometry.power = eval('P(PW).Value'); else, warning('Parameter P(PW).Value not in override. Default value will be used'); end           % (mW) injected laser power
    if exist('TA','var'), DOT.radiometry.acqtime = eval('P(TA).Value'); else, warning('Parameter P(TA).Value not in override. Default value will be used'); end           % (s) acquisition time
    if exist('OE','var'), DOT.radiometry.opteff = eval('P(OE).Value'); else, warning('Parameter P(OE).Value not in override. Default value will be used'); end           % typical efficiency of the optical path
    if exist('LAMBDA','var'), DOT.radiometry.lambda = eval('P(LAMBDA).Value'); else, warning('Parameter P(LAMBDA).Value not in override. Default value will be used'); end           % wavelenghts
    if exist('EA','var'), DOT.radiometry.area = eval('P(EA).Value'); else, warning('Parameter P(EA).Value not in override. Default value will be used'); end           % effective area of the detector
    if exist('QE','var'), DOT.radiometry.qeff = eval('P(QE).Value'); else, warning('Parameter P(QE).Value not in override. Default value will be used'); end           % quantum efficiency
    
    % permutation matrix
    %if exist('SD','var'), DOT.dmask = eval('P(SD).Value'); else, warning('Parameter P(SD).Value not in override. Default value will be used'); end
    
    % Type of Cutting for statistics
    if exist('CUT','var'), CUT_COUNTS = eval('P(CUT).Value'); else, warning('Parameter P(CUT).Value not in override. Default value will be used'); end % 0 non gated, 1 gated
    if exist('ND','var'), NumDelays = eval('P(ND).Value'); else, warning('Parameter P(ND).Value not in override. Default value will be used'); end % number of delays
    
    % Load Jacobian
    if exist('LJ','var')
        REC.solver.prejacobian.load = eval('P(LJ).Value');
        if REC.solver.prejacobian.load
            if iL == 1
                REC.solver.prejacobian.load = false;
                BuffOverride.NumTW = NUM_TW;
                BuffOverride.JSourcePos = DOT.Source.Pos;
                BuffOverride.JDetPos = DOT.Detector.Pos;
                BuffOverride.Jdmask = DOT.dmask;
                BuffOverride.Mua0 = REC.opt.mua0;
                BuffOverride.Mus0 = REC.opt.musp0;
                BuffOverride.NumDelays = NumDelays;
                BuffOverride.Roi = REC.time.roi;
                BuffOverride.JDOTgrid = DOT.grid;
                BuffOverride.JDOTtime = DOT.time;
            else
                if (... 
                        BuffOverride.NumTW == NUM_TW...
                        && isequaln(BuffOverride.JSourcePos,DOT.Source.Pos)...
                        && isequaln(BuffOverride.JDetPos,DOT.Detector.Pos)...
                        && isequaln(BuffOverride.Jdmask,DOT.dmask)...
                        && BuffOverride.Mua0 == REC.opt.mua0 ...
                        && BuffOverride.Mus0 == REC.opt.musp0 ...
                        && BuffOverride.NumDelays == NumDelays ...
                        && all(BuffOverride.Roi == REC.time.roi)...
                        && isequaln(BuffOverride.JDOTgrid,DOT.grid)...
                        && isequaln(BuffOverride.JDOTtime,DOT.time)...
                        )==false
                    BuffOverride.NumTW = NUM_TW;
                    BuffOverride.JSourcePos = DOT.Source.Pos;
                    BuffOverride.JDetPos = DOT.Detector.Pos;
                    BuffOverride.Jdmask = DOT.dmask;
                    BuffOverride.Mua0 = REC.opt.mua0;
                    BuffOverride.Mus0 = REC.opt.musp0;
                    BuffOverride.NumDelays = NumDelays;
                    BuffOverride.Roi = REC.time.roi;
                    BuffOverride.JDOTgrid = DOT.grid;
                    BuffOverride.JDOTtime = DOT.time;
                    REC.solver.prejacobian.load = false;
                else
                    REC.solver.prejacobian.load = true;
                end
            end
        end
    else
        warning('Parameter P(LJ).Value not in override. Default value will be used');
    end
    
    % Load Forward
    if exist('LF','var')
        LOAD_FWD_TEO =  eval('P(LF).Value');
        if LOAD_FWD_TEO
            if iL == 1
                LOAD_FWD_TEO = false;
                for ih = 1:NUM_HETE
                    BuffOverride.(['Hete' num2str(ih)]) = DOT.opt.(['hete' num2str(ih)]);
                end
                BuffOverride.FSourcePos = DOT.Source.Pos;
                BuffOverride.FDetPos = DOT.Detector.Pos;
                BuffOverride.Fdmask = DOT.dmask;
                BuffOverride.MuaB = DOT.opt.muaB;
                BuffOverride.MusB = DOT.opt.muspB;
                BuffOverride.FDOTgrid = DOT.grid;
                BuffOverride.FDOTtime = DOT.time;
            else
                HeteStructCond = true;
                for ih = 1:NUM_HETE
                    HeteStructCond = HeteStructCond && isequaln(BuffOverride.(['Hete' num2str(ih)]),DOT.opt.(['hete' num2str(ih)]));
                end
                if (...
                        HeteStructCond...
                        && isequaln(BuffOverride.FSourcePos,DOT.Source.Pos)...
                        && isequaln(BuffOverride.FDetPos,DOT.Detector.Pos)...
                        && isequaln(BuffOverride.Fdmask,DOT.dmask)...
                        && BuffOverride.MuaB == DOT.opt.muaB...
                        && BuffOverride.MusB == DOT.opt.muspB...
                        && isequaln(BuffOverride.FDOTgrid,DOT.grid)...
                        && isequaln(BuffOverride.FDOTtime,DOT.time)...
                        )==false
                    for ih = 1:NUM_HETE
                        BuffOverride.(['Hete' num2str(ih)]) = DOT.opt.(['hete' num2str(ih)]);
                    end
                    BuffOverride.FSourcePos = DOT.Source.Pos;
                    BuffOverride.FDetPos = DOT.Detector.Pos;
                    BuffOverride.Fdmask = DOT.dmask;
                    BuffOverride.MuaB = DOT.opt.muaB;
                    BuffOverride.MusB = DOT.opt.muspB;
                    BuffOverride.FDOTgrid = DOT.grid;
                    BuffOverride.FDOTtime = DOT.time;
                    LOAD_FWD_TEO = false;
                else
                    LOAD_FWD_TEO = true;
                end
            end
        end
    else
        warning('Parameter P(LF).Value not in override. Default value will be used');
    end
    warning('on','backtrace');
end