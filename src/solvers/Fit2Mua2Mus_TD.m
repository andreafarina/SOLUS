%==========================================================================
% This function contains a solver for fitting optical properties of
% 2 regions mesh using TOAST as forward and Matlab Optimization Toolbox
%
% Andrea Farina 02/18
%==========================================================================

function [bmua,bmus] = Fit2Mua2Mus_TD(solver,grid,mua0,mus0, n, ~,...
    Qpos,Mpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd,verbosity)
verbosity = 1;
self_norm = true;
INCL_ONLY = false;

%% initial setting the FEM problem
% create the mesh
mdim = [grid.Nx,grid.Ny,grid.Nz];
[vtx,idx,eltp] = mkslab([grid.x1,grid.y1,grid.z1;...
                    grid.x2,grid.y2,grid.z2],mdim);
hmesh = toastMesh(vtx,idx,eltp);
refind = n * ones(hmesh.NodeCount,1);
% create basis
bdim = mdim;% + 1;
hbasis = toastBasis(hmesh,bdim);
% map prior to mesh
priorM = hbasis.Map('B->M',double(solver.prior.refimage));
% create Q/M
ds = 1;
hmesh.SetQM(Qpos,Mpos);
qvec = hmesh.Qvec('Neumann','Gaussian',ds);
mvec = hmesh.Mvec('Gaussian',ds,0);
nQM = sum(dmask(:));
%% normalize data
if self_norm == true
        data = data * spdiags(1./sum(data)',0,nQM,nQM);
        ref = ref * spdiags(1./sum(ref)',0,nQM,nQM);
end
%% mask for excluding zeros
mask = (data(:) == 0) | (isnan(data(:)));

sd = ones(size(data));%sqrt(data);%%ones(size(proj));%proj(:);
%sd = ones(size(data));
data = data./sd;
data(mask) = [];
sd(mask) = [];

%% fitting procedure
if INCL_ONLY
    x0 = [mua0,mus0];
%    fitfun = @forward2;
else
    x0 = [mua0,mua0,mus0,mus0]; %start from homogeneous combination
%    fitfun = @forward;
end

% setting optimization
opts = optimoptions('lsqcurvefit',...
     'Jacobian','off',...
  ...'Algorithm','levenberg-marquardt',...
     'DerivativeCheck','off',...
     'MaxIter',100,'Display','iter-detailed',...%'FinDiffRelStep',[1e-3,1e-1],...;%,
     'TolFun',1e-10,'TolX',1e-10);

 x = lsqcurvefit(@forward,x0,[],data(:),[],[],opts);


%% display fit result
display(['mua_IN = ',num2str(x(1))]);
display(['musp_IN = ',num2str(x(2))]);
if ~INCL_ONLY
    display(['mua_BK = ',num2str(x(3))]);
    display(['musp_BK = ',num2str(x(4))]);
end

%% Map parameters back to basis
optmua = x(1) * priorM;
optmus = x(2) * priorM;
if INCL_ONLY
    optmua = optmua + (~priorM)*x0(1);
    optmus = optmus + (~priorM)*x0(2);
else
    optmua = optmua + (~priorM)*x(3);
    optmus = optmus + (~priorM)*x(4);
end
%% Map parameter back to basis
bmua = hbasis.Map('M->B',optmua(:));
bmus = hbasis.Map('M->B',optmus(:));




%% forward solvers
function [proj] = forward(x, ~)
        
        if INCL_ONLY
            mua = x(1) * priorM + mua0 * ( ~priorM);
            mus = x(2) * priorM + mus0 * ( ~priorM);
        else
            mua = x(1) * priorM + x(3) * ( ~priorM);
            mus = x(2) * priorM + x(4) * ( ~priorM);
        end
            
        
        
%         Mus = basis.Map('B->M',mus);
%         Mua = basis.Map('B->M',mua);
        
        [proj,~] = ProjectFieldTD(hmesh,qvec,mvec,dmask, mua,mus,0,0,refind,dt,nstep,0,0,'diff',0);
        proj = proj * spdiags(1./sum(proj)',0,nQM,nQM);
        
        if numel(irf)>1
            z = convn(proj,irf);
            nmax = max(nstep,numel(irf));
            proj = z(1:nmax,:);
            clear nmax
            if self_norm == true
                proj = proj * spdiags(1./sum(proj)',0,nQM,nQM);
            end
            clear z
        end
        %proj = circshift(proj,round(t0/dt));
        proj = WindowTPSF(proj,twin);
        if self_norm == true
            proj = proj * spdiags(1./sum(proj)',0,nQM,nQM);
        end
        proj(mask) = [];
        if verbosity
            % plot forward
            t = (1:numel(data)) * dt;
            figure(1003);
            semilogy(t,proj(:),'-',t,data(:),'.'),ylim([1e-3 1])
            drawnow;
        end
        proj = proj(:)./sd(:);
        
end


end