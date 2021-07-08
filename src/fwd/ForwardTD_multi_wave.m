function [phi,Area]=ForwardTD_multi_wave(grid,SourcePos,DetectorPos, dmask,...
        muaB,muspB,nB,Mua,Musp,A, dt, nstep, self_norm, geom, TYPE_FWD, radiometry,irf)
%dmask= reshape(dmask, [size(dmask,1)*size(dmask,2),size(dmask,3)] );
nQM = sum(dmask(:));   
phi = zeros(nstep,nQM);
Area = zeros(nQM,1);
nQM_lam = findMeasIndex(dmask);
if nargin<17
    irf = [];
end
    lambda = radiometry.lambda;
    for inl = 1:radiometry.nL
        fprintf(['<strong>------- Wavelength ',num2str(lambda(inl)),'-------</strong>\n'])
        idxQM = nQM_lam{inl};
        if isempty(Mua)||isempty(Musp)
            MuaData = []; MusData = [];
        else
            MuaData = squeeze(Mua(:,:,:,inl));
            MusData = squeeze(Musp(:,:,:,inl));
        end
        if ~isempty(irf)
            irf_lambda = irf(:,idxQM);
        else
            irf_lambda = [];
        end
        [phi_,Area_] = ForwardTD(grid,SourcePos, DetectorPos, dmask(:,:,inl),...
            muaB(inl), muspB(inl),nB, MuaData,...
            MusData, A, dt,...
            nstep, self_norm, geom, TYPE_FWD,irf_lambda);
            phi(:,idxQM) = phi_;
            Area(idxQM) = Area_;
        dummy_phi = phi_;
        if sum(dummy_phi(:)<0) > 0
            warning('off','verbose')
            warning('off','backtrace')
            warning([num2str(sum(dummy_phi(:)<0)),' elements of DataTD < 0!']);
            warning('on','verbose')
            warning('on','backtrace')
        end
    end
    phi = reshape(phi,[nstep,nQM]);
    

end
