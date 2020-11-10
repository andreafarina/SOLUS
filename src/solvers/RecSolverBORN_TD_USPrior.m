%==========================================================================
% This function contains solvers for DOT or fDOT.
% Andrea Farina 04/17
%==========================================================================

function [bmua,bmus] = RecSolverBORN_TD_USPrior(solver,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, type_fwd)
%% Jacobain options
USEGPU = 0;%gpuDeviceCount;
LOAD_JACOBIAN = solver.prejacobian.load;      % Load a precomputed Jacobian
geom = 'semi-inf';
% -------------------------------------------------------------------------
nQM = sum(dmask(:));
nwin = size(twin,1);
% -------------------------------------------------------------------------
[p,type] = ExtractVariables(solver.variables);
Jacobian = @(mua, mus) JacobianTD (grid, Spos, Dpos, dmask, mua, mus, n, A, ...
    dt, nstep, twin, irf, geom,type,type_fwd);
%% self normalise (probably useless because input data are normalize)
if self_norm == true
        data = data * spdiags(1./sum(data)',0,nQM,nQM);
        ref = ref * spdiags(1./sum(ref)',0,nQM,nQM);
end
%% Inverse solver
% homogeneous forward model
[proj, Aproj] = ForwardTD(grid,Spos, Dpos, dmask, mua0, mus0, n, ...
    [],[], A, dt, nstep, self_norm,...
    geom, 'linear');

% Convolution with IRF
if numel(irf)>1
    z = convn(proj,irf);
    nmax = max(nstep,numel(irf));
    proj = z(1:nmax,:);
    clear nmax z
end
    if self_norm == true
        proj = proj * spdiags(1./sum(proj,'omitnan')',0,nQM,nQM);
    end
    

proj = WindowTPSF(proj,twin);
proj = proj(:);
ref = ref(:);
data = data(:);

factor = proj./ref;
% load('factor_ref.mat')
% factor = repmat(factor,[nwin 1]);

factor = factor(:);
if self_norm == true
    factor = 1;
end
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
x0 = PrepareX0([mua0,1./(3*mus0)],grid.N,type);
x = ones(size(x0));

dphi = (data(:)-ref(:))./sd(~mask);%./ref(:);
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
    save(solver.prejacobian.path,'J');
    toc;
end

if self_norm == true
    for i=1:nQM
        sJ = sum(J((1:nwin)+(i-1)*nwin,:),'omitnan');
        sJ = repmat(sJ,nwin,1);
        sJ = spdiags(proj((1:nwin)+(i-1)*nwin),0,nwin,nwin) * sJ;
        J((1:nwin)+(i-1)*nwin,:) = (J((1:nwin)+(i-1)*nwin,:) - sJ)./Aproj(i);
    end
end

J = spdiags(1./sd(:),0,numel(sd),numel(sd)) * J;  % data normalisation
nsol = size(J,2);
%   parameter normalisation (scale x0)
J = J * spdiags(x0,0,length(x0),length(x0));
   
J(mask,:) = [];


%% Structured laplacian prior

%siz_prior = size(solver.prior.refimage);
%solver.prior(solver.prior == max(solver.prior(:))) = 1.1*min(solver.prior(:)); 
%solver.prior = solver.prior .* (1 + 0.01*randn(size(solver.prior)));
%[L,~] = StructuredLaplacianPrior(solver.prior.refimage,siz_prior(1),siz_prior(2),siz_prior(3));
k3d = Kappa(solver.prior.refimage,5);
[Dx,Dy,Dz] = gradientOperator(grid.dim,[1,1,1],[],'none');
L = [sqrt(k3d)*Dx;sqrt(k3d)*Dy;sqrt(k3d)*Dz];
%% Solver
disp('Calculating the larger singular value');
s = svds(J,1);
% structured Laplacian
L1 = [];
for ip = 1:p
     L1 = blkdiag(L1,L);
end

%% case Lcurve or direct solution
b = [dphi;zeros(p*3*nsol/p,1)];
if USEGPU
    gpu = gpuDevice; %#ok<UNRCH>
    disp('Using GPU');
    %J = gpuArray(J);
    b = gpuArray(b);
    %dphiG = gpuArray(dphi);
    %L1 = gpuArray(L1);
end

if numel(solver.tau)>1
    for i = 1:numel(solver.tau)
        alpha = solver.tau(i)*s(1); 
        disp(['Solving for tau = ',num2str(solver.tau(i))]);
        tic;
        if USEGPU
            A = [J;alpha*L1]; %#ok<UNRCH>
            A = gpuArray(A);
            dx = lsqr(A*1e10,b*1e10,1e-6,1000);
        else
            dx = lsqr([J;alpha*L1],b,1e-6,1000);
        end
        toc;
        %dx = [J;(alpha)*L]\[dphi;zeros(3*nsol,1)];
        res(i) = gather(norm(J*dx-dphi)); %#ok<AGROW>
        prior(i) = gather(norm(L1'*L1*dx)); %#ok<AGROW>
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
if USEGPU > 0
    A = [J;alpha*L1];
    A = gpuArray((A));
    dx = lsqr(A*1e10,b*1e10,1e-6,1000);
    dx = gather(dx);
else
    dx = lsqr([J;alpha*L1],b,1e-6,1000);
end

%==========================================================================
%%                        Add update to solution
%==========================================================================

x = x + dx;
%logx = logx + dx;
%x = exp(logx);
x = x.*x0;
[bmua,bmus] = XtoMuaMus(x,mua0,mus0,type);


end
