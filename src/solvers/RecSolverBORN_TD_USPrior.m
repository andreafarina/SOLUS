%==========================================================================
% This function contains solvers for DOT or fDOT.
% Andrea Farina 04/17
%==========================================================================

function [bmua,bmus] = RecSolverBORN_TD_USPrior(solver,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, type_fwd)
%% Jacobain options
LOAD_JACOBIAN = solver.prejacobian.load;      % Load a precomputed Jacobian
geom = 'semi-inf';
% -------------------------------------------------------------------------
nQM = sum(dmask(:));
nwin = size(twin,1);
% -------------------------------------------------------------------------
[p,type] = ExtractVariables(solver.variables);
Jacobian = @(mua, mus) JacobianTD (grid, Spos, Dpos, dmask, mua, mus, n, A, ...
    dt, nstep, twin, irf, geom,type,type_fwd);
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
    clear nmax
    if self_norm == true
        proj = proj * spdiags(1./sum(proj,'omitnan')',0,nQM,nQM);
    end
    clear z
end

proj = WindowTPSF(proj,twin);
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

siz_prior = size(solver.prior.refimage);
%solver.prior(solver.prior == max(solver.prior(:))) = 1.1*min(solver.prior(:)); 
%solver.prior = solver.prior .* (1 + 0.01*randn(size(solver.prior)));
[L,~] = StructuredLaplacianPrior(solver.prior.refimage,siz_prior(1),siz_prior(2),siz_prior(3));
%% Solver
s = svd(J);
alpha = solver.tau*s(1) %#ok<NOPRT>
%dx = [J;(alpha)*speye(nsol)]\[dphi;zeros(nsol,1)];
%dx = [J;(alpha)*L]\[dphi;zeros(3*nsol,1)];
dx = lsqr([J;repmat(alpha*L,1,p)],[dphi;zeros(3*nsol/p,1)],1e-6,1000);
%dx = lsqr([J;alpha*speye(nsol)],[dphi;zeros(nsol,1)],1e-6,100);
%==========================================================================
%%                        Add update to solution
%==========================================================================

x = x + dx;
%logx = logx + dx;
%x = exp(logx);
x = x.*x0;
[bmua,bmus] = XtoMuaMus(x,mua0,mus0,type);


end
