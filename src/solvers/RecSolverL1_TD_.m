%==========================================================================
% This function contains solvers for DOT or fDOT.
% To have available all the functionalty install REGU toolbox
% Andrea Farina 12/16
%==========================================================================

function [bmua,bmus] = RecSolverL1_TD(solver,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, ~)
%% Jacobain options
LOAD_JACOBIAN = false;      % Load a precomputed Jacobian
geom = 'semi-inf';
%% SOLVER PARAMETER CRITERION
LSQRit = 40; % solve LSQR with small number of iterations
LSQRtol = 1e-6; % tolerence of LSQR
% ISTA FLAGS
ISTA_FLAGS.FISTA= true;             
ISTA_FLAGS.pos = true;
ISTA_FLAGS.Wavelets =true;
ISTA_FLAGS.Iterates = true;


%% path
%rdir = ['../results/test/precomputed_jacobians/'];
jacdir = ['../results/test/precomputed_jacobians/'];
jacfile = 'J';
%mkdir(rdir);
%disp(['Intermediate results will be stored in: ' rdir])
%save([rdir 'REC'], '-v7.3','REC'); %save REC structure which describes the experiment
% -------------------------------------------------------------------------

bdim = (grid.dim);
nQM = sum(dmask(:));
nwin = size(twin,1);
Jacobian = @(mua, mus) JacobianTD (grid, Spos, Dpos, dmask, mua, mus, n, A, ...
    dt, nstep, twin, irf, geom);
%% Inverse solver
[proj, Aproj] = ForwardTD(grid,Spos, Dpos, dmask, mua0, mus0, n, ...
                [],[], A, dt, nstep, self_norm,...
                geom, 'homo');
if numel(irf)>1
    for i = 1:nQM
        z(:,i) = conv(proj(:,i),irf);
    end
    proj = z(1:numel(irf),:);
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
sd = sd(:).*sqrt(factor);
%sd = proj(:);
%sd = ones(size(proj(:)));
%% mask for excluding zeros
mask = ((ref(:).*data(:)) == 0) | ...
    (isnan(ref(:))) | (isnan(data(:)));
%mask = false(size(mask));


if ref == 0
    ref = proj(:);
end

ref(mask) = [];
data(mask) = [];


%sd(mask) = [];
% solution vector
x = ones(grid.N,1) * mua0; 
x0 = x;
p = length(x);
dphi = (data(:)-ref(:))./sd(~mask);%./ref(:);
%dphi = log(data(:)) - log(ref(:));
%save('dphi','dphi');
% ---------------------- Construct the Jacobian ---------------------------

if LOAD_JACOBIAN == true
    fprintf (1,'Loading Jacobian\n');
    tic;
    load([jacdir,jacfile])
    toc;
else
    fprintf (1,'Calculating Jacobian\n');
    tic;
    J = Jacobian ( mua0, mus0);
    save([jacdir,jacfile],'J');
    toc;
end

if ~isempty(solver.prior)
    d1 = (solver.prior(:) > 0 )&(solver.prior(:)==min(solver.prior(:)));
    d2 = (solver.prior(:) > min(solver.prior(:)))&(solver.prior(:)==max(solver.prior(:)));
    D = [d1(:),d2(:)];
    J = J * D;
end
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
% ------- to solve only for mua uncomment the following sentence ----------
%J(:,nsol+(1:nsol)) = 0;
proj(mask) = [];
J(mask,:) = [];

%% now ISTA (slow!)
lambda = sum(sum(J.^2));
xista = zeros(size(x0));
Niter = 5000/1;
tau = 2/lambda;
T = 1e-5;
tic
h = waitbar(0,'ISTA iterations');
for k = 1:Niter
    xista = SoftThresh(xista + tau*J'*(dphi - J*xista),T);
    waitbar(k/Niter);
end
toc;
disp('solved using ISTA');
dx = xista;
%% Let's try shrinkage-Newton
% xsn = zeros(size(x0));
% NSNit = ceil(5000/4);
% T = 1e-5;%0.02;
% alpha = 0;
% tic;
% h = waitbar(0,'iteration');
% for k = 1:NSNit
% %    xsn = SoftThresh(xsn + lsqr([J; alpha*lambda*speye(nsol)],[dphi-J*xsn;zeros(nsol,1)],1e-1,5),T );    
% %   xsn = SoftThresh(xsn + (J'*a + alpha*lambda*speye(nsol))\J'*(dphi-J*xsn),T );
% %xsn = SoftThresh(xsn + pcg(J'*J + alpha*lambda*speye(nsol),J'*(dphi-J*xsn),1e-6,100),T);
%     waitbar(k/NSNit);
% end
%toc;
%disp('solved using shrinkage-Newton');
%% update solution
x = x + dx;
%logx = logx + dx;
%x = exp(logx);

bmua = x;
bmus = ones(size(bmua)) * mus0;

end
