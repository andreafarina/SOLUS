%==========================================================================
% This function contains solvers for DOT or fDOT.
% To have available all the functionalty install REGU toolbox
% Andrea Farina 12/16
%==========================================================================

function [bmua,bmus] = RecSolverL1_TD(solver,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, fwd_type)
%% Jacobain options
LOAD_JACOBIAN = solver.prejacobian.load;      % Load a precomputed Jacobian
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
jacdir = ['../results/precomputed_jacobians/'];
jacfile = 'J';
% -------------------------------------------------------------------------

bdim = (grid.dim);
nQM = sum(dmask(:));
nwin = size(twin,1);
% -------------------------------------------------------------------------
[p,type_jac] = ExtractVariables(solver.variables);
Jacobian = @(mua, mus) JacobianTD (grid, Spos, Dpos, dmask, mua, mus, n, A, ...
    dt, nstep, twin, irf, geom,type_jac,fwd_type);
%% self normalise (probably useless because input data are normalize)
% if self_norm == true
%         data = data * spdiags(1./sum(data,'omitnan')',0,nQM,nQM);
%         ref = ref * spdiags(1./sum(ref,'omitnan')',0,nQM,nQM);
% end%% Inverse solver
[proj, Aproj] = ForwardTD(grid,Spos, Dpos, dmask, mua0, mus0, n, ...
                [],[], A, dt, nstep, self_norm,...
                geom, fwd_type);
if numel(irf)>1
    z = convn(proj,irf);
    nmax = max(nstep,numel(irf));
    proj = z(1:nmax,:);
    clear nmax z
end
    if self_norm == true
        proj = proj * spdiags(1./sum(proj)',0,nQM,nQM);
    end
    
proj = WindowTPSF(proj,twin);
proj = proj(:);
ref = ref(:);
data = data(:);
factor = proj./ref;

data = data .* factor;
ref = ref .* factor;
%% data scaling
sd = sd(:).*factor;
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
x0 = PrepareX0([mua0,1./(3*mus0)],grid.N,type_jac);
x = ones(size(x0));

dphi = (data(:)-ref(:))./sd(~mask);%./ref(:);
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


% if ~isempty(solver.prior.refimage)
%     d1 = (solver.prior.refimage(:) > 0 )&(solver.prior.refimage(:)==min(solver.prior.refimage(:)));
%     d2 = (solver.prior.refimage(:) > min(solver.prior.refimage(:)))&(solver.prior.refimage(:)==max(solver.prior.refimage(:)));
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
%   parameter normalisation (scale x0)
J = J * spdiags(x0,0,length(x0),length(x0));


%proj(mask) = [];
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
%     disp('--------------  solved using ISTA AF --------------');
%     %lambda = solver.tau;
%     %alpha = max(svd(J));
%     [dx,J] = ista(dphi,J,lambda,alpha,NISTAit);


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
x = x.*x0;
[bmua,bmus] = XtoMuaMus(x,mua0,mus0,type_jac);

end
