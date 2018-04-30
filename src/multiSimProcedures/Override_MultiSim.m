if exist('isComingFromInterface','var')
    warning('off','backtrace');
    % Inverse problem: parameters for the inverse solver
    if exist('TW','var'), NUM_TW = eval('P(TW).Value'); else, warning('Parameter P(TW).Value not in override. Default value will be used'); end
    if exist('MUA0','var'), REC.opt.mua0 = eval('P(MUA0).Value'); else, warning('Parameter P(MUA0).Value not in override. Default value will be used'); end
    if exist('MUS0','var'), REC.opt.musp0 = eval('P(MUS0).Value'); else, warning('Parameter P(MUS0).Value not in override. Default value will be used'); end
    if exist('NB','var'), REC.opt.nB = eval('P(NB).Value'); else, warning('Parameter P(NB).Value not in override. Default value will be used'); end
    
    % Regularization Parameters
    if exist('TAU','var'), REC.solver.tau = eval('P(TAU).Value'); else, warning('Parameter P(TAU).Value not in override. Default value will be used'); end
    if exist('SOLVTYPE','var'), REC.solver.type = eval('P(SOLVTYPE).Value'); else, warning('Parameter P(SOLVTYPE).Value not in override. Default value will be used'); end                                                                                             % 'born','GN': gauss-newton,
    if exist('EXPDELTA','var'), EXP_DELTA = eval('P(EXPDELTA).Value'); else, warning('Parameter P(EXPDELTA).Value not in override. Default value will be used'); end
    if exist('NOISETYPE','var'), DOT.time.noise = eval('P(NOISETYPE).Value'); else, warning('Parameter P(NOISETYPE).Value not in override. Default value will be used'); end
    
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
    if exist('MUAP','var'), DOT.opt.hete1.val = eval('P(MUAP).Value'); else, warning('Parameter P(MUAP).Value not in override. Default value will be used'); end
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
    if exist('ROISTART','var')
        REC.time.roi(1)= eval('P(ROISTART).Value');
        if REC.time.roi(1)<0
            if iL > 1
                nREC = REC;
                load([res_path,filename,'_', 'REC'],'REC')
                nREC.time.roi = REC.time.roi;
                clearvars REC; REC = nREC; clearvars nREC
                P(ROISTART).Default = REC.time.roi(1);
                if HIDE_FIG == 1
                    HideFig(true);
                end
            else
                REC.time=rmfield(REC.time,'roi');
                if HIDE_FIG == 1
                    HideFig(false);
                end
            end
        end
    else
        warning('Parameter P(ROISTART).Value not in override. Default value will be used');
    end % 0 non gated, 1 gated
    if exist('ROISTOP','var')
        REC.time.roi(2)=eval('P(ROISTOP).Value');
        if REC.time.roi(2)<0
            if iL > 1
                nREC = REC;
                load([res_path,filename,'_', 'REC'],'REC')
                nREC.time.roi = REC.time.roi;
                clearvars REC; REC = nREC; clearvars nREC
                P(ROISTOP).Default = REC.time.roi(2);
                if HIDE_FIG == 1
                    HideFig(true);
                end
            else
                REC.time=rmfield(REC.time,'roi');
                if HIDE_FIG == 1
                    HideFig(false);
                end
            end
        end
    else
        warning('Parameter P(ROISTOP).Value not in override. Default value will be used');
    end % 0 non gated, 1 gated
    
    if exist('ND','var')
        NumDelays = eval('P(ND).Value');
        %         if isfield(REC.time,'roi')
        %             if iL == 1
        %                 BuffOverride.Roi=REC.time.roi;
        %                 BuffOverride.ND=NumDelays;
        %                 %             if isfield(REC.time,'roi')
        %                 %                 button = questdlg({['Present ROI: ',num2str(REC.time.roi)],'Accept?'});
        %                 %                 if strcmpi(button,'no')
        %                 %                     REC.time=rmfield(REC.time,'roi');
        %                 %                 end
        %                 %             end
        %                 if HIDE_FIG == 1
        %                     HideFig(false);
        %                 end
        %             else
        %                 if (BuffOverride.ND == NumDelays)==false
        %                     if isfield(REC.time,'roi'), REC.time=rmfield(REC.time,'roi'); end
        %                     BuffOverride.ND = NumDelays;
        %                     if HIDE_FIG == 1
        %                         HideFig(false);
        %                     end
        %                 else
        %                     nREC = REC;
        %                     load([res_path,filename,'_', 'REC'],'REC')
        %                     nREC.time.roi = REC.time.roi;
        %                     clearvars REC; REC = nREC; clearvars nREC
        %                     if HIDE_FIG == 1
        %                         HideFig(true);
        %                     end
        %
        %                 end
        %             end
        %         end
    else
        warning('Parameter P(ND).Value not in override. Default value will be used');
    end % number of delays
    
    % Load Jacobian
    if exist('LJ','var')
        if eval('P(LJ).Value') == 0
            REC.solver.prejacobian.load = false;
        else
            if iL == 1
                REC.solver.prejacobian.load = false;
                BuffOverride.MuaB = REC.opt.mua0;
                BuffOverride.NumDelays = NumDelays;
                if isfield(REC.time,'Roi'), BuffOverride.Roi = REC.time.roi; end
                BuffOverride.MusB = REC.opt.musp0;
            else
                if ~isfield(BuffOverride,'Roi'), BuffOverride.Roi = [0 0]; end
                if (BuffOverride.MuaB == REC.opt.mua0 && BuffOverride.MusB == REC.opt.musp0 && BuffOverride.NumDelays == NumDelays && all(BuffOverride.Roi == REC.time.roi))==false
                    BuffOverride.NumDelays = NumDelays;
                    BuffOverride.MuaB = REC.opt.mua0;
                    BuffOverride.MusB = REC.opt.musp0;
                    BuffOverride.Roi = REC.time.roi;
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
        if eval('P(LF).Value') == 0
            LOAD_FWD_TEO = false;
        else
            if iL == 1
                LOAD_FWD_TEO = false;
                BuffOverride.IncPos = DOT.opt.hete1.c;
                BuffOverride.IncVal = DOT.opt.hete1.val;
                BuffOverride.IncSigma = DOT.opt.hete1.sigma;
            else
                if (all(BuffOverride.IncPos == DOT.opt.hete1.c) && BuffOverride.IncVal == DOT.opt.hete1.val && BuffOverride.IncSigma == DOT.opt.hete1.sigma)==false
                    BuffOverride.IncPos = DOT.opt.hete1.c;
                    BuffOverride.IncVal = DOT.opt.hete1.val;
                    BuffOverride.IncSigma = DOT.opt.hete1.sigma;
                    LOAD_FWD_TEO = false;
                else
                    LOAD_FWD_TEO = true;
                end
            end
        end
    else
        warning('Parameter P(LF).Value not in override. Default value will be used');
    end
    %if exist('LJ','var'), REC.solver.loadjacobian = eval('P(LJ).Value'); else, warning('Parameter P(LJ).Value not in override. Default value will be used'); end
    warning('on','backtrace');
end