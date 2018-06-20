function J = JacobianTD(grid,Spos, Dpos, dmask, muaB, muspB, n, ...
    A, dt, nstep, twin, irf, geom,type,type_fwd)
global mesh

v = 0.2999/n;
switch lower(type_fwd)
    case 'linear'
        Ns = size(Spos,1);
        %Nd = size(Dpos,1);
        Nopt = sum(dmask(:));
        nwin = size(twin,1);
        if strcmpi(type,'muad')
            n = 2*grid.N;
        else
            n = grid.N;
        end
        J = zeros(nwin*Nopt,n);
        t = (1:nstep) * dt;
        t = t';
        
        x = grid.x;%(1:end-1) + grid.dx/2;
        y = grid.y;%(1:end-1) + grid.dy/2;
        z = grid.z;%(1:end-1) + grid.dz/2;
        [XX,YY,ZZ]= ndgrid(x,y,z);
        clear x y z
        
        if irf == 0
            irf = 1;
        end
        disp(['calculating Jacobian (',geom,')']);
        %textprogressbar('progress: ');
        switch lower(geom)
            case 'infinite'
                
            case 'semi-inf'
                row_off = 0;
                for i = 1:Ns
                    ind_d = find(dmask(:,i));
                    %textprogressbar(i/Ns*100);
                    pause(0.01);
                    for j=1:numel(ind_d)
                        m = ind_d(j);
                        J(row_off + (1:nwin),:) = WindowTPSF(...
                            convn(...
                            J_Semi_Inf_PCBC_TR_3D_muaD(t, muaB,muspB,v,A,Spos(i,:),Dpos(m,:),XX,YY,ZZ,...
                            grid.dV,type), irf), twin);
                        row_off = row_off + nwin;
                    end
                end
                
            case 'slab'
        end
        %textprogressbar('done');
    case 'fem'
        N = mesh.hMesh.NodeCount;
        if ((abs(max(muaB(:))-min(muaB(:)))<eps)&&((abs(max(muspB(:))-min(muspB(:)))<eps)))
                mua = ones(N,1)*muaB;
                musp = ones(N,1) * muspB;
        end
         if strcmpi(type,'muad')
             GRADIENT = 'operator';
         else
             GRADIENT = 'none';
         end
        J = toastJacobianTimedomain_INPROGRESS(mesh.hMesh,grid.hBasis,...
            mesh.qvec, mesh.mvec, dmask, mua, musp, n*ones(N,1),dt,...
            nstep,twin,irf, GRADIENT);
        
end


