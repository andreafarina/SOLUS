%==========================================================================
% This function contains solvers for DOT or fDOT.
% To have available all the functionalty install REGU toolbox
% Andrea Farina 12/16
%==========================================================================

function [bmua,bmus] = RecSolverBORN_TD(solver,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, ~)
%% Jacobain options
LOAD_JACOBIAN = false;      % Load a precomputed Jacobian
geom = 'semi-inf';
%% REGULARIZATION PARAMETER CRITERION
NORMDIFF = 'ref';   % 'ref', 'sd'
REGU = 'external';        % 'lcurve', 'gcv', 'external'
BACKSOLVER = 'tikh';    % 'tikh', 'tsvd', 'simon', 'gmres', 'pcg'
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
if strcmpi(NORMDIFF,'sd'), dphi = (data(:)-ref(:))./sd(~mask); end
if strcmpi(NORMDIFF,'ref'), dphi = (data(:)-ref(:))./ref(:); end
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

if strcmpi(NORMDIFF,'sd'), J = spdiags(1./sd(:),0,numel(sd),numel(sd)) * J;  end % data normalisation
if strcmpi(NORMDIFF,'ref'), J = spdiags(1./proj(:),0,numel(proj),numel(proj)) * J;  end % data normalisation
nsol = size(J,2);
%   parameter normalisation (map to log)
%     for i = 1:p
%         J(:,i) = J(:,i) * x(i);
%     end
% ------- to solve only for mua uncomment the following sentence ----------
%J(:,nsol+(1:nsol)) = 0;
proj(mask) = [];
J(mask,:) = [];

%% Solver 
if ~strcmpi((BACKSOLVER),'simon')
[U,s,V]=csvd(J);     % compact SVD (Regu toolbox)
        figure(402);
        picard(U,s,dphi);    % Picard plot (Regu toolbox)    
end
if (~strcmpi(REGU,'lcurve')&&(~strcmpi(REGU,'gcv')))
    alpha = solver.tau * s(1);
end
if ~exist('alpha','var')
    figure(403);
    if strcmpi(REGU,'lcurve')
       alpha = l_curve(U,s,dphi);%,BACKSOLVER);  % L-curve (Regu toolbox)
    elseif strcmpi(REGU,'gcv')
        alpha = gcv(U,s,dphi);%,BACKSOLVER)
    end
    disp(['alpha=' num2str(alpha)]);
end

switch lower(BACKSOLVER)
    case 'tikh'
        disp('Tikhonov');
        
        [dx,rho] = tikhonov(U,s,V,dphi,alpha);
    case 'tsvd'
        disp('TSVD');
        
        [dx,rho] = tsvd(U,s,V,dphi,alpha);
    case 'simon'
        disp('Simon');
        tic;
        s1 = svds(J,1);
        toc;
        alpha = solver.tau * s1;
        dx = [J;sqrt(alpha)*speye(nsol)]\[dphi;zeros(nsol,1)];
        %dx = [dx;zeros(nsol,1)];
   
    
%rho
%cond(J)
%dx = [dx;zeros(nsol,1)];
%dx = [J;alpha*eye(2*nsol)]\[dphi;zeros(2*nsol,1)];
%dx = lsqr(J,dphi,1e-4,100);
%dx = pcg(@(x)JTJH(x,J,eye(2*nsol),alpha),J'*dphi,1e-4,100);
    case 'gmres'
        disp('gmres')
        if strcmpi(REGU,'external')
            lambda = sum(sum(J.^2));
            alpha = lambda * alpha;
        end
        dx = gmres(@(x)JTJH(x,J,eye(nsol),alpha),J'*dphi,[],1e-6,100);

    case 'pcg'
        disp('pcg')
        tic;
        if strcmpi(REGU,'external')
            lambda = sum(sum(J.^2));
            alpha = alpha * lambda;
            disp(['alpha=' num2str(alpha)]);
        end
        %dx = pcg(J'*J + alpha*speye(nsol),J'*dphi,1e-6,100);   
        dx = pcg(@(x)JTJH(x, J, speye(nsol),alpha),J'*dphi,1e-6,100);   
        
        toc;

%dx = J' * ((J*J' + alpha*eye(length(data)))\dphi);
end  
%==========================================================================
%%                        Add update to solution
%==========================================================================
if ~isempty(solver.prior)
    dx = D * dx;
end
x = x + dx;
%logx = logx + dx;
%x = exp(logx);

bmua = x;
bmus = ones(size(bmua)) * mus0;

%% display the reconstructions
% figure(304);
% %figure('Position',get(0,'ScreenSize'));
% ShowRecResults(grid,reshape(bmua,bdim),grid.z1,grid.z2,grid.dz,1,...
%     min(bmua),max(bmua));
% suptitle('Recon Mua');
% %export_fig '../Results/20151111/rec_mua_1incl.pdf' -transparent
% 
% figure(305);
% %figure('Position',get(0,'ScreenSize'));
% ShowRecResults(grid,reshape(bmus,bdim),grid.z1,grid.z2,grid.dz,1,...
%     min(bmus),max(bmus));
% suptitle('Recon Mus');
% %export_fig '../Results/20151111/rec_mus_1incl.pdf' -transparent
% drawnow;
% tilefigs;
% %pause
% x = x0;
% %save([rdir 'sol'], 'xiter');
% %save([rdir 'err_pattern'],'erri');
end





