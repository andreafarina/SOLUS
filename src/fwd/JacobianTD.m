function J = JacobianTD(grid,Spos, Dpos, dmask, muaB, muspB, refind, ...
    A, dt, nstep, twin, irf, geom,type,type_fwd,selfnorm,logdata)
global mesh
DEBUG = 0;
v = 0.2999/refind;
nQM = sum(dmask(:));
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
        J = (zeros(nwin*Nopt,n));
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
                    %pause(0.01);
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
                mua = ones(N,1) .* muaB;
                musp = ones(N,1) .* muspB;
        else
            mua = muaB;
            musp = muspB;
        end
         if strcmpi(type,'muad')
             GRADIENT = 'operator';%'matlab';%'toast';%'operator';
         else
             GRADIENT = 'none';
         end
         BB = mesh.hMesh.BoundingBox();
        bmin = BB(1,:); bmax = BB(2,:);
        bdim = grid.hBasis.Dims();bdim = bdim';
        vscale = prod((bmax(:)-bmin(:))./(bdim(:)-1));%1;%2.3;


        J = toastJacobianTimedomain_INPROGRESS(mesh.hMesh,grid.hBasis,...
            mesh.qvec, mesh.mvec, dmask, mua, musp, refind*ones(N,1),dt,...
            nstep,twin,irf, GRADIENT) * vscale * dt;
        if DEBUG == 1
            Ns = size(Spos,1);
            %Nd = size(Dpos,1);
            Nopt = sum(dmask(:));
            nwin = size(twin,1);
            if strcmpi(type,'muad')
                n = 2*grid.N;
            else
                n = grid.N;
            end
            J2 = zeros(nwin*Nopt,n);
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
            
            row_off = 0;
            for i = 1:Ns
                ind_d = find(dmask(:,i));
                %textprogressbar(i/Ns*100);
                %pause(0.01);
                for j=1:numel(ind_d)
                    m = ind_d(j);
                    J2(row_off + (1:nwin),:) = WindowTPSF(...
                        convn(...
                        J_Semi_Inf_PCBC_TR_3D_muaD(t, muaB,muspB,v,A,Spos(i,:),Dpos(m,:),XX,YY,ZZ,...
                        grid.dV,type), irf), twin);
                    row_off = row_off + nwin;
                end
            end
            ss = 0;
            while ss~=-1
                ss = input('Jac row?');
                if ss>0
            figure,subplot(1,2,1),plot(squeeze(J2(ss,:))),subplot(1,2,2),
            plot(J(ss,:));%,subplot(2,2,3),plot(squeeze(J2(67,:))),subplot(2,2,4),
            %plot(J(67,:))
                end
            end
        end
        
end

%% self-norm and logdata
if selfnorm||logdata
    [proj, Aproj] = ForwardTD(grid,Spos, Dpos, dmask, muaB, muspB, refind, ...
    [],[], A, dt, nstep, selfnorm, geom, type_fwd, irf);
    proj = WindowTPSF(proj,twin);
end

%% case self-normalized
if selfnorm == true
    for i=1:nQM
        sJ = sum(J((1:nwin)+(i-1)*nwin,:),'omitnan');
        sJ = repmat(sJ,nwin,1);
        sJ = spdiags(proj((1:nwin)+(i-1)*nwin)',0,nwin,nwin) * sJ;
        J((1:nwin)+(i-1)*nwin,:) = (J((1:nwin)+(i-1)*nwin,:) - sJ)./Aproj(i);
    end
end
%% case logData
if logdata
    for i = 1:nQM*nwin
        J(i,:) = J(i,:)./proj(i);
    end
end



