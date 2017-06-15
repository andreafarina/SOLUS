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
SOLVER = 'ISTA';
%ISTA: "normal" (i.e. pixel) ISTA
%FISTA: "normal" (i.e. pixel) FISTA
%ADMM: "normal" (i.e. pixel) ADMM
NISTAit  = 1000;
NFISTAit = 200;
NADMMit = 10;
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
alpha = 1e-6 * lambda;
tau = 2/lambda;
T = alpha * tau;
%% set up operators for Shrinkage - pixel case

pops.W    = @(x) x;
pops.Winv = @(c) c;
pops.WTr  = @(c) c;
pops.B    = @(x)x;
pops.Binv = @(x)x;
pops.A    = @(x)J*x;
pops.ATr  = @(y) J'*y;
pops.ntc  = nsol;

pops.ndat = numel(data(:));
pops.nsol = nsol;
pops.Dims = grid.dim;

%% now ISTA (slow!)
if strcmpi(SOLVER,'ista')
    xista = zeros(nsol,1);
    Niter = NISTAit;
 %   T = alpha*tau;
    ISTA_FLAGS.Rayleighstep = true;
    ISTA_FLAGS.pos = false;
    ISTA_FLAGS.FISTA = false;
    ISTA_FLAGS.Wavelets = false;
    ISTA_FLAGS.Iterates = true;
    tic;
%    [x,c,postista] = DotWavDeconvByISTA(hBasis,Ja,y,alpha,Niter,tau,ISTA_FLAGS);
    [xx,c,postista] = LinearShrinkage(pops,dphi,alpha,NISTAit,tau,ISTA_FLAGS);
    toc;
    disp('--------------  solved using ISTA  --------------');
    xista = xx(:,end);
%     for k = 1:NISTAit
%         lerrista(k) = norm(y - Ja*x(:,k));
%         xistim = reshape(hBasis.Map('S->B',x(:,k)),bx,by);
%         xerrista(k) = norm(xistim-tgtmuaim );
%         perrista(k) = norm(c(:,k),1);        
%     end
dx = xista;
end
%% now FISTA (faster?)
if strcmpi(SOLVER,'fista')
    Niter = NFISTAit;
    ISTA_FLAGS.Rayleighstep = true;
    ISTA_FLAGS.pos = false;
    ISTA_FLAGS.FISTA = true;
    ISTA_FLAGS.Wavelets = false;
    ISTA_FLAGS.Iterates = true;
    tic;
%    [x,c,postfista] = DotWavDeconvByISTA(hBasis,Ja,y,alpha,Niter,tau,ISTA_FLAGS);
    [xx,c,postfista] = LinearShrinkage(pops,dphi,alpha,NFISTAit,tau,ISTA_FLAGS);
    toc;
    disp('--------------  solved using FISTA  --------------');
    xfistay = xx(:,end);
%     for k = 1:NFISTAit
%         lerrfista(k) = norm(y - Ja*x(:,k));
%         xfistim = reshape(hBasis.Map('S->B',x(:,k)),bx,by);
%         xerrfista(k) = norm(xfistim-tgtmuaim );
%         perrfista(k) = norm(c(:,k),1);        
%     end
    dx = xfistay;
end

%% ADMM - pixels
if strcmpi(SOLVER,'admm')
ISTA_FLAGS.Rayleighstep = true;
    ISTA_FLAGS.pos = false;
    ISTA_FLAGS.FISTA = false;
    ISTA_FLAGS.Wavelets = false;
    ISTA_FLAGS.Iterates = true;
rho = 1e-2*lambda; % 1e-2*alpha/T; % what is best way to set rho ?
sqrho = sqrt(rho);
tic; 
[xpadm,padm,postpadm] = LinearADMM(pops,dphi,alpha,NADMMit,tau,rho,LSQRtol,LSQRit,ISTA_FLAGS);
toc;
xpadmin = xpadm(:,end);
% for k = 1:NADMMit
%     lerrpadm(k) = norm(y - Ja*xpadm(:,k));
%     xpadmmim = reshape(hBasis.Map('S->B',xpadm(:,k)),bx,by);
%     xerrpadm(k) = norm(xpadmmim-tgtmuaim );
%     perrpadm(k) = norm(padm(:,k),1);
% end
disp('------------- solved using ADMM  -------------- ');
dx = xpadmin;
end



%% update solution
x = x + dx;
%logx = logx + dx;
%x = exp(logx);

bmua = x;
bmus = ones(size(bmua)) * mus0;

end
