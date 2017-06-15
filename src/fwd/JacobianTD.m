function J = JacobianTD(grid,Spos, Dpos, dmask, muaB, muspB, n, ...
    A, dt, nstep, twin, irf, geom)

v = 0.2999/n;

Ns = size(Spos,1);
%Nd = size(Dpos,1);
Nopt = sum(dmask(:));
nwin = size(twin,1);
J = zeros(nwin*Nopt,grid.N);
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
switch lower(geom)
    case 'infinite'
        
    case 'semi-inf'
        row_off = 0;
        for i = 1:Ns
            ind_d = find(dmask(:,i));
            for j=1:numel(ind_d)
                m = ind_d(j);
                J(row_off + (1:nwin),:) = WindowTPSF(...
                  convn(...
                  J_Semi_Inf_PCBC_TR_3D(t, muaB,muspB,v,A,Spos(i,:),Dpos(m,:),XX,YY,ZZ,...
                    grid.dV), irf), twin);
                row_off = row_off + nwin;
            end
        end
        
    case 'slab'
end
end


