function [phi,Area]=ForwardTD_multi_wave(grid,SourcePos,DetectorPos, dmask,...
        muaB,muspB,nB,Mua,Musp,A, dt, nstep, self_norm, geom, TYPE_FWD, radiometry)
<<<<<<< HEAD
nQM = sum(dmask(:));   
phi = zeros(nstep,nQM*radiometry.nL);
Area = zeros(nQM,radiometry.nL);

    
    for inl = 1:radiometry.nL
        fprintf(['<strong>------- Wavelength ',num2str(radiometry.lambda(inl)),'-------</strong>\n'])
        meas_set = (1:nQM)+(inl-1)*nQM;
        if isempty(Mua)||isempty(Musp)
            MuaData = []; MusData = [];
        else
            MuaData = squeeze(Mua(:,:,:,inl));
            MusData = squeeze(Musp(:,:,:,inl));
        end
        [phi(:,meas_set),Area(1:nQM,inl)] = ForwardTD(grid,SourcePos, DetectorPos, dmask,...
            muaB(inl), muspB(inl),nB, MuaData,...
            MusData, A, dt,...
            nstep, self_norm, geom, TYPE_FWD);
        dummy_phi = phi(:,meas_set);
        if sum(dummy_phi(:)<0) > 0
            warning('off','verbose')
            warning('off','backtrace')
            warning([num2str(sum(dummy_phi(:)<0)),' elements of DataTD < 0!']);
            warning('on','verbose')
            warning('on','backtrace')
        end
    end

end
=======

nQM = sum(dmask(:));   
phi = zeros(nstep,nQM*radiometry.nL);
Area = zeros(nQM,radiometry.nL);
for inl = 1:radiometry.nL
    fprintf(['<strong>------- Wavelength ',num2str(radiometry.lambda(inl)),'-------</strong>\n'])
    meas_set = (1:nQM)+(inl-1)*nQM;
    if isempty(Mua)||isempty(Musp)
        MuaData = []; MusData = [];
    else
        MuaData = squeeze(Mua(:,:,:,inl));
        MusData = squeeze(Musp(:,:,:,inl));
    end
    [phi(:,meas_set),Area(1:nQM,inl)] = ForwardTD(grid,SourcePos, DetectorPos, dmask,...
        muaB(inl), muspB(inl),nB, MuaData,...
        MusData, A, dt,...
        nstep, self_norm, geom, TYPE_FWD);
    dummy_phi = phi(:,meas_set);
    if sum(dummy_phi(:)<0) > 0
        warning('off','verbose')
        warning('off','backtrace')
        warning([num2str(sum(dummy_phi(:)<0)),' elements of DataTD < 0!']);
        warning('on','verbose')
        warning('on','backtrace')
    end
end
end
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
