function J = JacobianTD_multiwave(grid,Spos, Dpos, dmask, muaB, muspB, n, ...
    A, dt, nstep, twin, irf, geom,type,type_fwd,radiometry) %#ok<INUSL>

nQM = sum(dmask(:));
nV = grid.N;
nTW = size(twin,1);
ishift = 0;
switch lower(type)
    case 'mua'
        J = zeros(nTW*nQM*radiometry.nL,nV);
    case 'd'
        J = zeros(nTW*nQM*radiometry.nL,nV);
    case 'muad'
        J = zeros(nTW*nQM*radiometry.nL,nV*2);
        ishift = nV;
end
if contains(lower(type),'mua')
    for inl = 1:radiometry.nL
        twin_set = (1:2)+(inl-1)*2;
        meas_set = (1:nTW*nQM)+(inl-1)*nTW*nQM;
        disp('-------');
        fprintf('<strong>Absorption operations</strong>\n');
        fprintf(['<strong>------- Wavelength ',num2str(radiometry.lambda(inl)),'-------</strong>\n'])
        J(meas_set,(1:nV)) = JacobianTD(grid,Spos, Dpos, dmask, muaB(inl), muspB(inl), n, ...
            A, dt, nstep, twin(:,twin_set), irf(:,inl), geom,'mua','linear');
    end
end
if contains(lower(type),'d')
    for inl = 1:radiometry.nL
        twin_set = (1:2)+(inl-1)*2;
        meas_set = (1:nTW*nQM)+(inl-1)*nTW*nQM;
        disp('-------');
        fprintf('<strong>Scattering operations</strong>\n');
        fprintf(['<strong>------- Wavelength ',num2str(radiometry.lambda(inl)),'-------</strong>\n'])
        J(meas_set,ishift+(1:nV)) = JacobianTD(grid,Spos, Dpos, dmask, muaB(inl), muspB(inl), n, ...
            A, dt, nstep, twin(:,twin_set), irf(:,inl), geom,'d','linear');
    end
end
end


