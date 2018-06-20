%==========================================================================
% This function contains solvers for DOT or fDOT.
% To have available all the functionalty install REGU toolbox
% Andrea Farina 12/16
%==========================================================================

function [bmua,bmus] = RecSolverBORN_TD(solver,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, fwd_type)
%global factor
%ref = 0;
%% Jacobain options
LOAD_JACOBIAN = solver.prejacobian.load;      % Load a precomputed Jacobian
geom = 'semi-inf';
%% REGULARIZATION PARAMETER CRITERION
NORMDIFF = 'sd';   % 'ref', 'sd'
REGU = 'external'; % 'lcurve', 'gcv', 'external'
BACKSOLVER = 'tikh'; % 'tikh', 'tsvd', 'discrep','simon', 'gmres', 'pcg', 'lsqr'
% -------------------------------------------------------------------------
nQM = sum(dmask(:));
nwin = size(twin,1);
% -------------------------------------------------------------------------
[p,type_jac] = ExtractVariables(solver.variables);
Jacobian = @(mua, mus) JacobianTD (grid, Spos, Dpos, dmask, mua, mus, n, A, ...
    dt, nstep, twin, irf, geom,type_jac,fwd_type);
%% self normalise
if self_norm == true
        data = data * spdiags(1./sum(data)',0,nQM,nQM);
end
%% Inverse solver
[proj, Aproj] = ForwardTD(grid,Spos, Dpos, dmask, mua0, mus0, n, ...
                [],[], A, dt, nstep, self_norm,...
                geom, fwd_type);
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
factor = factor(:);


data = data .* factor;
ref = proj(:);%ref .* factor;
%% data scaling
sd = sd(:).*(factor);
%sd = proj(:);%ref(:);%proj(:);
%sd = ones(size(proj(:)));
if ref == 0 %#ok<BDSCI>
    ref = proj(:);
end
%% mask for excluding zeros
mask = ((ref(:).*data(:)) == 0) | ...
    (isnan(ref(:))) | (isnan(data(:)));
%mask = false(size(mask));




ref(mask) = [];
data(mask) = [];

%sd(mask) = [];
% solution vector
x0 = PrepareX0([mua0,1./(3*mus0)],grid.N,type_jac);
x = ones(size(x0));

if strcmpi(NORMDIFF,'sd'), dphi = (data(:)-ref(:))./sd(~mask); end
if strcmpi(NORMDIFF,'ref'), dphi = (data(:)-ref(:))./ref(:); end
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

if ~isempty(solver.prior.refimage)
    %d1 = (solver.prior(:) > 0 )&(solver.prior(:)==min(solver.prior(:)));
    d1 = solver.prior.refimage(:) > mean(solver.prior.refimage(:));
    d2 = ~d1; 
    if mean(solver.prior.refimage(d1))<mean(solver.prior.refimage(d2))
        d1 = ~d1;
        d2 = ~d2;
    end
    %(solver.prior(:) > min(solver.prior(:)))&(solver.prior(:)==max(solver.prior(:)));
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
%   parameter normalisation (scale x0)
J = J * spdiags(x0,0,length(x0),length(x0));
    
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
    
end
disp(['alpha = ' num2str(alpha)]);
switch lower(BACKSOLVER)
    case 'tikh'
        disp('Tikhonov');
        
        [dx,~] = tikhonov(U,s,V,dphi,alpha);
    case 'tsvd'
        disp('TSVD');
        
        [dx,~] = tsvd(U,s,V,dphi,alpha);
    case 'discrep'
        disc_value = norm(sd(~mask))*10;
        disp(['Discrepancy principle with value=' num2str(disc_value)]);
            dx = discrep(U,s,V,dphi,disc_value);
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
        dx = pcg(@(x)JTJH(x, J, speye(nsol),alpha),J'*dphi,1e-6,1000);   
    case 'lsqr'
        disp('lsqr');
        tic;
        if strcmpi(REGU,'external')
            lambda = sum(sum(J.^2));
            alpha = alpha * lambda;
            disp(['alpha=' num2str(alpha)]);
        end
        dx = lsqr([J;alpha*speye(nsol)],[dphi;zeros(nsol,1)],1e-6,100);

        toc;

%dx = J' * ((J*J' + alpha*eye(length(data)))\dphi);
end  
%==========================================================================
%%                        Add update to solution
%==========================================================================
if ~isempty(solver.prior.refimage)
    dx = D * dx;
end
x = x + dx;
%logx = logx + dx;
%x = exp(logx);
x = x.*x0;
[bmua,bmus] = XtoMuaMus(x,mua0,mus0,type_jac);


end




