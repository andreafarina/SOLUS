%==========================================================================
% This function contains solvers for DOT or fDOT.
% To have available all the functionalty install REGU toolbox
% Andrea Farina 12/16
%==========================================================================

function [bmua,bmus,bconc,bA,bB] = RecSolverTK0_spectral_TD(solver,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, fwd_type,radiometry,spe,conc0,a0,b0)

%global factor
%ref = 0;
%% Jacobain options
LOAD_JACOBIAN = solver.prejacobian.load;      % Load a precomputed Jacobian
geom = 'semi-inf';
xtransf = '(x/x0)'; %log(x),x,log(x/x0)
type_ratio = 'gauss';   % 'gauss', 'born', 'rytov';
type_ref = 'theor';      % 'theor', 'meas', 'area'

%% REGULARIZATION PARAMETER CRITERION
NORMDIFF = 'sd';   % 'ref', 'sd'
REGU = 'lcurve'; % 'lcurve', 'gcv', 'external'
BACKSOLVER = 'tikh'; % 'tikh', 'tsvd', 'discrep','simon', 'gmres', 'pcg', 'lsqr'
% -------------------------------------------------------------------------
nQM = sum(dmask(:));
nwin = size(twin,1);
% -------------------------------------------------------------------------
[p,type_jac] = ExtractVariables_spectral(solver.variables,spe);

% if rytov and born jacobian is normalized to proj
if strcmpi(type_ratio,'rytov')||strcmpi(type_ratio,'born')
    logdata = true;
else
    logdata = false;
end
Jacobian = @(mua, mus) JacobianTD_multiwave_spectral (grid, Spos, Dpos, dmask, mua, mus, n, A, ...
    dt, nstep, twin, irf, geom,type_jac,fwd_type,radiometry,spe,self_norm,logdata);

% homogeneous forward model
[proj, ~] = ForwardTD_multi_wave(grid,Spos, Dpos, dmask, mua0, mus0, n, ...
    [],[], A, dt, nstep, self_norm,...
    geom, fwd_type,radiometry,irf);

%% Window TPSF for each wavelength
dummy_proj = zeros(size(twin,1),nQM*radiometry.nL);
for inl = 1:radiometry.nL
    meas_set = (1:nQM)+(inl-1)*nQM; twin_set = (1:2)+(inl-1)*2;
    proj_single = proj(:,meas_set);
    proj_single = WindowTPSF(proj_single,twin(:,twin_set));
    dummy_proj(:,meas_set) = proj_single;
end
proj = dummy_proj;
clear dummy_proj

[dphi,sd] = PrepareDataFitting(data,ref,sd,type_ratio,type_ref,proj);

% creat mask for nan, isinf
mask = (isnan(dphi(:))) | (isinf(dphi(:)));
dphi(mask) = [];

% solution vector
[x0,x] = PrepareX_spectral(spe,grid.N,type_jac,xtransf);

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


% sd jacobian normalization
J = spdiags(1./sd(:),0,numel(sd),numel(sd)) * J;
nsol = size(J,2);

% parameter normalisation (scale x0)
if ~strcmpi(xtransf,'x')
    J = J * spdiags(x0,0,length(x0),length(x0));
end

J(mask,:) = [];

%% Solver
%if ~strcmpi((BACKSOLVER),'simon')
if (~strcmpi(REGU,'lcurve')&&(~strcmpi(REGU,'gcv')))
    disp('Calculating larger singular value');
    s = svds(J,1)
    alpha = solver.tau * s;
end
if (strcmpi(BACKSOLVER,'tikh'))
    disp('Calculating compact SVD');
    [U,s,V]=csvd(J);     % compact SVD (Regu toolbox)
    figure(402);
    picard(U,s,dphi);    % Picard plot (Regu toolbox)
end

if ~exist('alpha','var')
    fh = figure(403);
    fh.NumberTitle = 'off'; fh.Name = 'L-curve';
    if strcmpi(REGU,'lcurve')
        alpha = l_curve(U,s,dphi);%,BACKSOLVER);  % L-curve (Regu toolbox)
    elseif strcmpi(REGU,'gcv')
        alpha = gcv(U,s,dphi);%,BACKSOLVER)
    end
    savefig(fh,[fh.Name '.fig'])
end
disp(['alpha = ' num2str(alpha), 'tau = ',num2str(alpha/s(1))]);
disp('Solving...')
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
        %         tic;
        %         s1 = svds(J,1);
        %         toc;
        %         alpha = solver.tau * s1;
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
            %lambda = sum(sum(J.^2));
            %alpha = alpha * lambda;
            disp(['alpha=' num2str(alpha)]);
        end
        
        
        dx = lsqr([J;alpha*speye(nsol)],[dphi;zeros(nsol,1)],1e-6,1000);
        
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

x = BackTransfX(x,x0,xtransf);

[bmua,bmus,bconc,bAB] = XtoMuaMus_spectral(x,mua0,mus0,type_jac,spe,conc0,a0,b0);
bA = bAB(:,1);bB = bAB(:,2);

end




