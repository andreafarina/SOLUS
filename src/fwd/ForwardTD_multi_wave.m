function [phi,Area]=ForwardTD_multi_wave(grid,SourcePos,DetectorPos, dmask,...
        muaB,muspB,nB,Mua,Musp,A, dt, nstep, self_norm, geom, TYPE_FWD, radiometry)
nQM = sum(dmask(:));   
phi = zeros(nstep,nQM,radiometry.nL);
Area = zeros(nQM,radiometry.nL);

    lambda = radiometry.lambda;
    for inl = 1:radiometry.nL
        fprintf(['<strong>------- Wavelength ',num2str(lambda(inl)),'-------</strong>\n'])
        %meas_set = (1:nQM)+(inl-1)*nQM;
        if isempty(Mua)||isempty(Musp)
            MuaData = []; MusData = [];
        else
            MuaData = squeeze(Mua(:,:,:,inl));
            MusData = squeeze(Musp(:,:,:,inl));
        end
        %[phi(:,meas_set),Area(1:nQM,inl)] 
        [phi_,Area_] = ForwardTD(grid,SourcePos, DetectorPos, dmask,...
            muaB(inl), muspB(inl),nB, MuaData,...
            MusData, A, dt,...
            nstep, self_norm, geom, TYPE_FWD);
            phi(:,:,inl) = phi_;
            Area(:,inl) = Area_;
        dummy_phi = phi_;
        if sum(dummy_phi(:)<0) > 0
            warning('off','verbose')
            warning('off','backtrace')
            warning([num2str(sum(dummy_phi(:)<0)),' elements of DataTD < 0!']);
            warning('on','verbose')
            warning('on','backtrace')
        end
    end
    phi = reshape(phi,[nstep,nQM*radiometry.nL]);
    

end
