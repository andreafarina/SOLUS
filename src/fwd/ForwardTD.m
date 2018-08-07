function [phi, Area] = ForwardTD(grid,Spos, Dpos, dmask, muaB, muspB, n, ...
    Mua, Musp, A, dt, nstep, self_norm, geom, FWD_TYPE)
global mesh
DEBUG = 0;
v = 0.2999/n;

Ns = size(Spos,1);
%Nd = size(Dpos,1);
Nopt = sum(dmask(:));
phi = zeros(nstep,Nopt);
t = (1:nstep) * dt;
t = t';
DiffB = 1./(3*muspB);
Diff = 1./(3*Musp);
switch lower(FWD_TYPE)
    case 'linear'
        % understand homo/hetero
        if (~isempty(Mua))&&(~isempty(Musp))
            
            if ((abs(max(Mua(:))-min(Mua(:)))<eps)&&((abs(max(Musp(:))-min(Musp(:)))<eps)))
                type = 'homo';
            else
                type = 'born';
            end
        else
            type = 'homo';
        end
        disp(['calculating Forward(',type,', ',geom,')']);
        
        switch lower(type)
            case 'homo'
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
                                phi(row_off + (1:nstep)) = SemiInfinite_TR(t, Spos(i,:),...
                                    Dpos(m,:),muaB,muspB,v,A);
                                row_off = row_off + nstep;
                            end
                        end
                        
                    case 'slab'
                        
                        
                end
                
            case 'born'
                x = grid.x;%(1:end-1) + grid.dx/2;
                y = grid.y;%(1:end-1) + grid.dy/2;
                z = grid.z;%(1:end-1) + grid.dz/2;
                [XX,YY,ZZ]= ndgrid(x,y,z);
                clear x y z
                
                % dimension permutation so to have f(x,y,z)
                %XX=ipermute(XX,[2 1 3]);
                %YY=ipermute(YY,[2 1 3]);
                %ZZ=ipermute(ZZ,[2 1 3]);
                %Mua=ipermute(Mua,[2,1,3]);
                % grid to voxel
                %Mua = Mua(1:end-1,1:end-1,1:end-1);
                switch lower(geom)
                    case 'infinite'
                        
                    case 'semi-inf'
                        row_off = 0;
                        for i = 1:Ns
                            ind_d = find(dmask(:,i));
                            %textprogressbar(i/Ns*100);
                            %pause(0.01);
                            for jj=1:numel(ind_d)
                                m = ind_d(jj);
                                phi(row_off + (1:nstep)) = SemiInfinite_TR(t, Spos(i,:), Dpos(m,:),muaB,muspB,v,A) + ...
                                    J_Semi_Inf_PCBC_TR_3D_muaD(t, muaB,muspB,v,A,Spos(i,:),Dpos(m,:),XX,YY,ZZ,...
                                    grid.dV) * [(Mua(:) - muaB);(Diff(:) - DiffB)];
                                row_off = row_off + nstep;
                            end
                        end
                        
                    case 'slab'
                end
                
        end
        
        %textprogressbar('done');
    case 'fem'
        %% PLACE here FEM forward solver
        if (~isempty(Mua))&&(~isempty(Musp))
            mua = grid.hBasis.Map('B->M',Mua);
            musp = grid.hBasis.Map('B->M',Musp);
        else
            mua = ones(mesh.hMesh.NodeCount,1).*muaB;
            musp = ones(mesh.hMesh.NodeCount,1).*muspB;
        end
        [phi,~] = ProjectFieldTD(mesh.hMesh,mesh.qvec,mesh.mvec,...
            dmask, mua,musp,0,0,n.*ones(size(mesh.opt.mua)),dt,nstep,0,0,'diff',0);
        if DEBUG == 1
            row_off = 0;
            for i = 1:Ns
                ind_d = find(dmask(:,i));
                %textprogressbar(i/Ns*100);
                %pause(0.01);
                for j=1:numel(ind_d)
                    m = ind_d(j);
                    phi2(row_off + (1:nstep)) = SemiInfinite_TR(t, Spos(i,:),...
                        Dpos(m,:),muaB,muspB,v,A);
                    row_off = row_off + nstep;
                end
            end
            figure,subplot(1,2,1),semilogy(1:numel(phi),phi(:),1:numel(phi2),phi2(:)),...
                legend('fem','analytic'),ylim(max(phi(:))*[1e-4 10]),
            subplot(1,2,2),plot((phi(:)-phi2(:)))
        end
        
end
if nargout > 1
    Area = sum(phi,'omitnan');
end
if self_norm == true
    Area = sum(phi,'omitnan');
    phi = bsxfun(@times,phi,1./Area);
end


