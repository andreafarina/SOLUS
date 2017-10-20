%==========================================================================
% This function contains solvers for DOT or fDOT.
% Andrea Farina 04/17
%==========================================================================

function [bmua,bmus] = RecSolverBORN_TD_USPrior(solver,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, ~)
%% Jacobain options
LOAD_JACOBIAN = solver.prejacobian.load;      % Load a precomputed Jacobian
geom = 'semi-inf';
% -------------------------------------------------------------------------
nQM = sum(dmask(:));
nwin = size(twin,1);
Jacobian = @(mua, mus) JacobianTD (grid, Spos, Dpos, dmask, mua, mus, n, A, ...
    dt, nstep, twin, irf, geom);
%% Inverse solver
% homogeneous forward model
[proj, Aproj] = ForwardTD(grid,Spos, Dpos, dmask, mua0, mus0, n, ...
    [],[], A, dt, nstep, self_norm,...
    geom, 'homo');

% Convolution with IRF
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

proj = WindowTPSF(proj,twin);
proj = proj(:);
ref = ref(:);
data = data(:);
factor = proj./ref;

data = data .* factor;
ref = ref .* factor;
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
x = ones(grid.N,1) * mua0;

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

% translate the spatial structure into 2 binary complementary masks
% if ~isempty(solver.prior)
%     d1 = (solver.prior(:) > 0 )&(solver.prior(:)==min(solver.prior(:)));
%     d2 = (solver.prior(:) > min(solver.prior(:)))&(solver.prior(:)==max(solver.prior(:)));
%     D = [d1(:),d2(:)];
%     J = J * D;
% end
if self_norm == true
    for i=1:nQM
        sJ = sum(J((1:nwin)+(i-1)*nwin,:));
        sJ = repmat(sJ,nwin,1);
        sJ = spdiags(proj((1:nwin)+(i-1)*nwin),0,nwin,nwin) * sJ;
        J((1:nwin)+(i-1)*nwin,:) = (J((1:nwin)+(i-1)*nwin,:) - sJ)./Aproj(i);
    end
end

J = spdiags(1./sd(:),0,numel(sd),numel(sd)) * J;  % data normalisation
nsol = size(J,2);
%   parameter normalisation (map to log)
%     for i = 1:p
%         J(:,i) = J(:,i) * x(i);
%     end
%proj(mask) = [];
J(mask,:) = [];


%% Structured laplacian prior
siz_prior = size(solver.prior);
solver.prior(solver.prior == max(solver.prior(:))) = 1.1*min(solver.prior(:)); 
%solver.prior = solver.prior .* (1 + 0.01*randn(size(solver.prior)));
[L,~] = StructuredLaplacianPrior(solver.prior,siz_prior(1),siz_prior(2),siz_prior(3));
%% Solver
s = svd(J);
alpha = solver.tau*s(1) %#ok<NOPRT>
%dx = [J;(alpha)*speye(nsol)]\[dphi;zeros(nsol,1)];
%dx = [J;(alpha)*L]\[dphi;zeros(3*nsol,1)];
dx = lsqr([J;alpha*L],[dphi;zeros(3*nsol,1)],1e-6,1000);
%dx = lsqr([J;alpha*speye(nsol)],[dphi;zeros(nsol,1)],1e-6,100);
%==========================================================================
%%                        Add update to solution
%==========================================================================

x = x + dx;
%logx = logx + dx;
%x = exp(logx);

bmua = x;
bmus = ones(size(bmua)) * mus0;


end
