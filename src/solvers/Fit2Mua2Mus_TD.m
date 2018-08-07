%==========================================================================
% This function contains a solver for fitting optical properties of
% 2 regions mesh using TOAST as forward and Matlab Optimization Toolbox
%
% Andrea Farina 02/18
%==========================================================================

function [bmua,bmus, OUTPUT] = Fit2Mua2Mus_TD(solver,grid,mua0,mus0, n, ~,...
    Qpos,Mpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd,verbosity)
verbosity = 1;
self_norm = true;
INCL_ONLY = false;
MUA_ONLY = true; musIN = 1.5;%real value of scattering to be used for MUA_ONLY

%% initial setting the FEM problem
% create the mesh
mdim = [grid.Nx,grid.Ny,grid.Nz] ;
[vtx,idx,eltp] = mkslab([grid.x1,grid.y1,grid.z1;...
                    grid.x2,grid.y2,grid.z2],mdim);
hmesh = toastMesh(vtx,idx,eltp);
refind = n * ones(hmesh.NodeCount,1);
% create basis
bdim = mdim;

% if priormask does not retunr back the same number of elemets as the DOT grid
%bdim = size(solver.prior.refimage);

hbasis = toastBasis(hmesh,bdim, 'LINEAR');

% map prior to mesh
priorM = hbasis.Map('B->M',double(solver.prior.refimage));

% set intermidiate values to 1;
% inter_UP = find(priorM >= 0.01);
% inter_DW = find( priorM < 0.5);
% priorM(inter_UP) = 1;
% priorM(inter_DW) = 0;
% create Q/M

Qds = 1; % width of Sources 
Mds = 1.5; % width of Detectors
hmesh.SetQM(Qpos,Mpos);
qvec = hmesh.Qvec('Neumann','Gaussian',Qds);
mvec = hmesh.Mvec('Gaussian',Mds, n);

%     mtot = mvec(:,1) + mvec(:,2) + mvec(:,3) + mvec(:,4) + mvec(:,5) + mvec(:,6) + mvec(:,7) + mvec(:,8); % FOR DISPLAY
%     qtot = qvec(:,1) + qvec(:,2) + qvec(:,3) + qvec(:,4) + qvec(:,5) + qvec(:,6) + qvec(:,7) + qvec(:,8);
%     tot = (max(qtot) / max(mtot)) * mtot + qtot;
%     hmesh.Display(qtot);

nQM = sum(dmask(:));
%% normalize data
if self_norm == true
        data = data * spdiags(1./sum(data,'omitnan')',0,nQM,nQM);
        ref = ref * spdiags(1./sum(ref,'omitnan')',0,nQM,nQM);
        sd = sqrt(data) * spdiags(1./sum(data,'omitnan')',0,nQM,nQM); 
end
%% mask for excluding zeros
mask = (data(:) == 0) | (isnan(data(:)));

%sd = ones(size(data));%;%%ones(size(proj));%proj(:);
sd = ones(size(data));
data = data./sd;
data(mask) = [];
sd(mask) = [];
%data2 (mask) = [];

%% fitting procedure
if INCL_ONLY
    x0 = [mua0,mus0]; lb = [0,0]; ub = [1, 10];
%    fitfun = @forward2;
elseif MUA_ONLY
    x0 = [mua0, mua0];lb = [0,0]; ub = [1, 10];
else
   x0 =[mua0 + 0.001,mus0 + 0.1,mua0,mus0];  % [muaIN, musIN, muaOUT, musOUT]
   %x0 = [0.001,1,0.001,1]; %start from homogeneous combination
   lb = [0,0,0,0]; ub = [1, 10, 1, 10];
%    fitfun = @forward;
end

% setting optimization
opts = optimoptions('lsqcurvefit',...
     'Jacobian','off',...
     ...'Algorithm','trust-region-reflective',...
     'DerivativeCheck','off',...
     'MaxIter',100,'Display','iter-detailed',...%'FinDiffRelStep',[1e-4,1e-2],...%,
     'TolFun',1e-10,'TolX',1e-10);

[x,~,~,~,OUTPUT] = lsqcurvefit(@forward,x0,[],data(:),lb,ub,opts);x_lsqr = x;
%[x] = fmincon(@forwardfmincon,x0); OUTPUT = 0; 
%x_fmin = x 
x_lsqr
%x = x0;
%x = bicg(@forward, data(:));
%x_pcg = x;

%% display fit result

if INCL_ONLY
    display(['mua_IN = ',num2str(x(1))]);
    display(['musp_IN = ',num2str(x(2))]);
elseif MUA_ONLY
    display(['mua_IN = ',num2str(x(1))]);
    display(['mua_BK = ',num2str(x(2))]);
else
    display(['mua_BK = ',num2str(x(3))]);
    display(['musp_BK = ',num2str(x(4))]);
    display(['mua_IN = ',num2str(x(1))]);
    display(['musp_IN = ',num2str(x(2))]);
end


%% Map parameters back to basis
if INCL_ONLY
    optmua = x(1) * priorM;
    optmus = x(2) * priorM;
    optmua = optmua + (1-priorM)*x0(1);
    optmus = optmus + (1-priorM)*x0(2);
elseif MUA_ONLY
    optmua = x(1) * priorM;
    optmua = optmua +(1-priorM)*x(2);
    optmus = priorM*musIN + (1-priorM)*mus0;
else
    optmua = x(1) * priorM;
    optmus = x(2) * priorM;
    optmua = optmua + (1-priorM)*x(3);
    optmus = optmus + (1-priorM)*x(4);
end

% if priorMask does not return back the same number of elements 
% hbasis_out = toastBasis(hmesh,mdim);
% bmua = hbasis_out.Map('M->B',optmua(:));
% bmus = hbasis_out.Map('M->B',optmus(:));

bmua = hbasis.Map('M->B',optmua(:));
bmus = hbasis.Map('M->B',optmus(:));

%bmua = optmua(:);
%bmus = optmus(:);

%% Delete Mesh and Basis
% hbasis_out.delete;
hbasis.delete;
hmesh.delete;
clearvars -except bmua bmus OUTPUT
return;

%% forward solvers
function [proj] = forward(x, ~)
        
        if INCL_ONLY
            mua = x(1) * priorM + mua0 * ( 1-priorM);
            mus = x(2) * priorM + mus0 * ( 1-priorM);
        elseif MUA_ONLY
            mua = x(1) * priorM + x(2) * ( 1-priorM);
            mus = musIN * priorM + mus0 * ( 1-priorM);
        else
            mua = x(1) * priorM + x(3) * ( 1-priorM);
            mus = x(2) * priorM + x(4) * ( 1-priorM);
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
                proj = proj * spdiags(1./sum(proj,'omitnan')',0,nQM,nQM);
            end
            clear z
        end
        %proj = circshift(proj,round(t0/dt));
        proj = WindowTPSF(proj,twin);
        if self_norm == true
            proj = proj * spdiags(1./sum(proj,'omitnan')',0,nQM,nQM);
        end
        proj(mask) = [];
        proj = proj(:)./sd(:);
        if verbosity
            % plot forward
            t = (1:numel(data)) * dt;
            figure(1003);
            semilogy(t,proj(:),'-',t,data(:),'.'),ylim([1e-3 1])
            title(['||proj-data||=',num2str(norm(proj-data(:)))])
            drawnow,
            x
        end
        
        
end

    function [out ] = forwardfmincon(x)
    
        out = norm(forward(x) - data(:));
    end
end