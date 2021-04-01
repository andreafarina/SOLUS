function J = JacobianTD_multiwave_spectral(grid,Spos, Dpos, dmask, muaB, muspB, n, ...
    A, dt, nstep, twin, irf, geom,type,type_fwd,radiometry,spe,self_norm,logdata) %#ok<INUSL>

nQM = sum(dmask(:));
nV = grid.N;
nTW = size(twin,1);
ishift = 0;
switch lower(type)
    case 'mua'
        J = zeros(nTW*nQM,nV*spe.nCromo);
    case 'd'
        J = zeros(nTW*nQM,nV*2);
    case 'muad'
        J = zeros(nTW*nQM,nV*spe.nCromo+nV*2);
        ishift = nV*spe.nCromo;
end
lambda = radiometry.lambda;
if contains(lower(type),'mua')
    ext_coeff0 = spe.ext_coeff0;
    %Ja = zeros(nTW*nQM,nV*spe.nCromo);
    idxmeas = findMeasIndex(dmask);
    parfor inl = 1:radiometry.nL
        meas_set = idxmeas{inl};
        disp('-------');
        fprintf('<strong>Absorption operations</strong>\n');
        fprintf(['<strong>------- Wavelength ',num2str(lambda(inl)),' -------</strong>\n'])
        
        %J(meas_set,(1:nV*spe.nCromo))
        Ja{inl} = JacobianTD(grid,Spos, Dpos, dmask(:,:,inl), muaB(inl), muspB(inl), n, ...
            A, dt, nstep, twin(:,:,meas_set), irf(:,meas_set), geom,'mua','linear',self_norm,logdata)*kron(ext_coeff0(inl,:),speye(nV));
        %clear ext_coeff_mat
    end
    J(:,(1:nV*spe.nCromo)) = cell2mat(Ja');%reshape(permute(Ja,[1,3,2]),nTW*nQM,[]);
    %J(:,(1:nV*spe.nCromo)) = reshape(permute(Ja,[1,3,2]),nTW*nQM,[]);
    clear Ja
end
if contains(lower(type),'d')
    b0 = spe.opt.b0;
    a0 = spe.opt.a0;
    lambda0 = radiometry.lambda0;
    %Js = (zeros(nTW*nQM,nV*2,radiometry.nL));
    idxmeas = findMeasIndex(dmask);
    for inl = 1:radiometry.nL
        meas_set = idxmeas{inl};
        disp('-------');
        fprintf('<strong>Scattering operations</strong>\n');
        fprintf(['<strong>------- Wavelength ',num2str(lambda(inl)),' -------</strong>\n'])
        dD = -1/(3*(muspB(inl))^2);
        %Js(meas_set,ishift+(1:nV*2))
        Js{inl}= JacobianTD(grid,Spos, Dpos, dmask(:,:,inl), muaB(inl), muspB(inl), n, ...
            A, dt, nstep, twin(:,:,meas_set), irf(:,meas_set), geom,'d','linear',self_norm,logdata)*dD*...
            kron([(lambda(inl)./lambda0).^(-b0),(-a0.*(lambda(inl)./lambda0).^(-b0))*...
            log(lambda(inl)./lambda0)],speye(nV));
    end
    J(:,ishift+(1:nV*2)) = cell2mat(Js');
    %J(:,ishift+(1:nV*2)) = reshape(permute(Js,[1,3,2]),nTW*nQM*radiometry.nL,[]);
end
end


