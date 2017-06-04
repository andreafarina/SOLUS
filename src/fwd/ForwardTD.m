function [phi, Area] = ForwardTD(grid,Spos, Dpos, dmask, muaB, muspB, n, ...
    Mua, Musp, A, dt, nstep, self_norm, geom, type)

v = 0.2999/n;

Ns = size(Spos,1);
%Nd = size(Dpos,1);
Nopt = sum(dmask(:));
phi = zeros(nstep,Nopt);
t = (1:nstep) * dt;
t = t';
switch lower(type)
    case 'homo'
        switch lower(geom)
            case 'infinite'
                
            case 'semi-inf'
                row_off = 0;
                for i = 1:Ns
                    ind_d = find(dmask(:,i));
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
                    for j=1:numel(ind_d)
                        m = ind_d(j);
                        phi(row_off + (1:nstep)) = SemiInfinite_TR(t, Spos(i,:), Dpos(m,:),muaB,muspB,v,A) + ... 
                        J_Semi_Inf_PCBC_TR_3D(t, muaB,muspB,v,A,Spos(i,:),Dpos(m,:),XX,YY,ZZ,...
                            grid.dV) * (Mua(:) - muaB);
                        row_off = row_off + nstep;
                    end
                end
                
            case 'slab'
        end

end
if nargout > 1        
    Area = sum(phi);
end


