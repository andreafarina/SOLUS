function J = JacobianTD_multiwave_spectral(grid,Spos, Dpos, dmask, muaB, muspB, n, ...
    A, dt, nstep, twin, irf, geom,type,type_fwd,radiometry,spe) %#ok<INUSL>

nQM = sum(dmask(:));
nV = grid.N;
nTW = size(twin,1);
ishift = 0;
switch lower(type)
    case 'mua'
        J = zeros(nTW*nQM*radiometry.nL,nV*spe.nCromo);
    case 'd'
        J = zeros(nTW*nQM*radiometry.nL,nV*2);
    case 'muad'
        J = zeros(nTW*nQM*radiometry.nL,nV*spe.nCromo+nV*2);
        ishift = nV*spe.nCromo;
end
if contains(lower(type),'mua')
    for inl = 1:radiometry.nL
        twin_set = (1:2)+(inl-1)*2;
        meas_set = (1:nTW*nQM)+(inl-1)*nTW*nQM;
        disp('-------');
        fprintf('<strong>Absorption operations</strong>\n');
        fprintf(['<strong>------- Wavelength ',num2str(radiometry.lambda(inl)),'-------</strong>\n'])
        ext_coeff_mat = [];
        for ic = 1:spe.nCromo
            ext_coeff_mat = [ext_coeff_mat spdiags(repmat(spe.ext_coeff0(inl,ic),nV,1),0,nV,nV)];
        end
        J(meas_set,(1:nV*spe.nCromo)) = JacobianTD(grid,Spos, Dpos, dmask, muaB(inl), muspB(inl), n, ...
            A, dt, nstep, twin(:,twin_set), irf(:,inl), geom,'mua','linear')*ext_coeff_mat;
        clear ext_coeff_mat
    end
end
if contains(lower(type),'d')
    for inl = 1:radiometry.nL
        twin_set = (1:2)+(inl-1)*2;
        meas_set = (1:nTW*nQM)+(inl-1)*nTW*nQM;
        disp('-------');
        fprintf('<strong>Scattering operations</strong>\n');
        fprintf(['<strong>------- Wavelength ',num2str(radiometry.lambda(inl)),'-------</strong>\n'])
        dD = -1/(3*(muspB(inl))^2);
        J(meas_set,ishift+(1:nV*2)) = JacobianTD(grid,Spos, Dpos, dmask, muaB(inl), muspB(inl), n, ...
            A, dt, nstep, twin(:,twin_set), irf(:,inl), geom,'d','linear')*...
            [spdiags(repmat(dD.*(radiometry.lambda(inl)./radiometry.lambda0).^(-spe.opt.b0),nV,1),0,nV,nV) spdiags(...
            repmat(-dD.*(spe.opt.a0.*(radiometry.lambda(inl)./radiometry.lambda0).^(-spe.opt.b0))*log(radiometry.lambda(inl)./radiometry.lambda0),nV,1),0,nV,nV)];
    end
end
end


