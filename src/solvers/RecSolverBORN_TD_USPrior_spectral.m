%==========================================================================
% This function contains solvers for DOT or fDOT.
% Andrea Farina 04/17
%==========================================================================

function [bmua,bmus,bconc,bA,bB] = RecSolverBORN_TD_USPrior_spectral(solver,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, fwd_type,radiometry,spe,conc0,a0,b0)
%% Jacobain options
%% Jacobain options
USEGPU = gpuDeviceCount;
LOAD_JACOBIAN = solver.prejacobian.load;      % Load a precomputed Jacobian
geom = 'semi-inf';
%% REGULARIZATION PARAMETER CRITERION
NORMDIFF = 'sd';   % 'ref', 'sd'
% -------------------------------------------------------------------------
nQM = sum(dmask(:));
nwin = size(twin,1);
% -------------------------------------------------------------------------
[p,type_jac] = ExtractVariables_spectral(solver.variables,spe);
Jacobian = @(mua, mus) JacobianTD_multiwave_spectral (grid, Spos, Dpos, dmask, mua, mus, n, A, ...
    dt, nstep, twin, irf, geom,type_jac,fwd_type,radiometry,spe);
%% Inverse solver
% homogeneous forward model
[proj, Aproj] = ForwardTD_multi_wave(grid,Spos, Dpos, dmask, mua0, mus0, n, ...
    [],[], A, dt, nstep, self_norm,...
    geom, fwd_type,radiometry);

% Convolution with IRF
if numel(irf)>1
    for inl = 1:radiometry.nL
        meas_set =(1:nQM)+(inl-1)*nQM;
        z = convn(proj(:,meas_set),irf(:,inl));
        nmax = max(nstep,numel(irf(:,inl)));
        if inl == 1
            proj(1:nmax,meas_set) = z(1:nmax,:);
            proj(nmax+1:end,:) = [];
        else
            proj(:,meas_set) = z(1:nmax,:);
        end
        clear nmax z
    end
end
if self_norm == true
    for inl = 1:radiometry.nL
        meas_set = (1:nQM)+(inl-1)*nQM;
        proj(:,meas_set) = proj(:,meas_set)*spdiags(1./sum(proj(:,meas_set),1,'omitnan')',0,nQM,nQM);
    end
end
dummy_proj = zeros(size(twin,1),nQM*radiometry.nL);
for inl = 1:radiometry.nL
    meas_set = (1:nQM)+(inl-1)*nQM; twin_set = (1:2)+(inl-1)*2;
    proj_single = proj(:,meas_set);
    proj_single = WindowTPSF(proj_single,twin(:,twin_set));
    dummy_proj(:,meas_set) = proj_single;
end
proj = dummy_proj;
proj = proj(:);
ref = ref(:);
data = data(:);

factor = proj./ref;
% load('factor_ref.mat')
% factor = repmat(factor,[nwin 1]);

factor = factor(:);

data = data .* factor;
ref = proj(:);%ref .* factor;
%% data scaling
sd = sd(:).*factor;%sqrt(factor);   % Because of the Poisson noise
%sd = proj(:);
%sd = ones(size(proj(:)));
%% mask for excluding zeros,nan,inf
mask = ((ref(:).*data(:)) == 0) | ...
    (isnan(ref(:))) | (isnan(data(:))) | ...
    (isinf(data(:))) | (isinf(ref(:)));
%mask = false(size(mask));

if ref == 0 %#ok<*BDSCI>
   ref = proj(:);
end

ref(mask) = [];
data(mask) = [];

%sd(mask) = [];
% solution vector
x0 = PrepareX0_spectral(spe,grid.N,type_jac);
x = ones(size(x0));

if strcmpi(NORMDIFF,'sd'), dphi = (data(:)-ref(:))./sd(~mask); end
if strcmpi(NORMDIFF,'ref'), dphi = (data(:)-ref(:))./ref(:); end
%sd = proj(:);
%dphi = log(data(:)) - log(ref(:));
%save('dphi','dphi');
% ---------------------- Construct the Jacobian ---------------------------
if LOAD_JACOBIAN == true
    fprintf (1,'Loading Jacobian\n');
    tic;
    %load([jacdir,jacfile])
    load(solver.prejacobian.path);
    toc;
else
    %fprintf (1,'Calculating Jacobian\n');
    tic;
    J = Jacobian ( mua0, mus0);
    [jpath,jname,jext] = fileparts(solver.prejacobian.path);
    if ~exist(jpath,'dir')
        mkdir(jpath)
    end
    save(solver.prejacobian.path,'J','-v7.3');
    toc;
end

if self_norm == true
    for inl = 1:radiometry.nL
        for i=1:nQM
            sJ = sum(J((1:nwin)+(i-1)*nwin+(inl-1)*nQM*nwin,:),1,'omitnan');
            sJ = repmat(sJ,nwin,1);
            sJ = spdiags(proj((1:nwin)+(i-1)*nwin+(inl-1)*nQM*nwin),0,nwin,nwin) * sJ;
            J((1:nwin)+(i-1)*nwin+(inl-1)*nQM*nwin,:) =...
                (J((1:nwin)+(i-1)*nwin+(inl-1)*nQM*nwin,:) - sJ)./Aproj(i,inl);
        end
    end
end

if strcmpi(NORMDIFF,'sd'), J = spdiags(1./sd(:),0,numel(sd),numel(sd)) * J;  end % data normalisation
if strcmpi(NORMDIFF,'ref'), J = spdiags(1./proj(:),0,numel(proj),numel(proj)) * J;  end % data normalisation

nsol = size(J,2);
%   parameter normalisation (scale x0)
J = J * spdiags(x0,0,length(x0),length(x0));
   
J(mask,:) = [];


%% Structured laplacian prior

%siz_prior = size(solver.prior.refimage);
%solver.prior(solver.prior == max(solver.prior(:))) = 1.1*min(solver.prior(:)); 
%solver.prior = solver.prior .* (1 + 0.01*randn(size(solver.prior)));
%[L,~] = StructuredLaplacianPrior(solver.prior.refimage,siz_prior(1),siz_prior(2),siz_prior(3));
%% new gradient
k3d = Kappa(solver.prior.refimage,5);
[Dy,Dx,Dz] = gradientOperator(permute(grid.dim,[2,1,3])',...
   [grid.dy,grid.dx,grid.dz],[],'none');
L = [sqrt(k3d)*Dx;sqrt(k3d)*Dy;sqrt(k3d)*Dz];

%% Solver
disp('Calculating singular values');
s = svds(J,1);

L1 = [];
for ip = 1:p
     L1 = blkdiag(L1,L);
end
%% case Lcurve or direct solution
if USEGPU
    gpu = gpuDevice;
    disp('Using GPU');
    J = gpuArray(J);
    b = gpuArray([full(dphi);zeros(p*3*nsol/p,1)]);
    %dphiG = gpuArray(dphi);
    L1 = gpuArray(full(L1));
end

if numel(solver.tau)>1
    for i = 1:numel(solver.tau)
        alpha = solver.tau(i)*s(1); 
        disp(['Solving for tau = ',num2str(solver.tau(i))]);
        tic;
        dx = lsqr([J;alpha*L1],[dphi;zeros(p*3*nsol/p,1)],1e-6,1000);
        toc;
        %dx = [J;(alpha)*L]\[dphi;zeros(3*nsol,1)];
        res(i) = gather(norm(J*dx-dphi)); %#ok<AGROW>
        prior(i) = gather(norm(L1*dx)); %#ok<AGROW>
        figure(144),loglog(res,prior,'-o'),title('L-curve');
        text(res,prior,num2cell(solver.tau(1:i)));xlabel('residual');ylabel('prior');
    end
           
    tau = solver.tau;
    save('LcurveData','res','prior','tau');
%     tau_suggested = l_corner(flip(res)',flip(prior)',flip(tau));
%     disp(['Suggested tau = ',num2str(tau_suggested)]);
    pause;
    tau_sel = inputdlg('Choose tau');
    solver.tau = str2double(tau_sel{1});
end
%% final solution
alpha = solver.tau * s(1);
disp(['Solving for tau = ',num2str(solver.tau)]);
dx = lsqr([J;alpha*L1],[dphi;zeros(p*3*nsol/p,1)],1e-6,1000);
if USEGPU>0
    dx = gather(dx);
end
%==========================================================================
%%                        Add update to solution
%==========================================================================

x = x + dx;
%logx = logx + dx;
%x = exp(logx);
x = x.*x0;
[bmua,bmus,bconc,bAB] = XtoMuaMus_spectral(x,mua0,mus0,type_jac,spe,conc0,a0,b0);
bA = bAB(:,1);bB = bAB(:,2);



end
