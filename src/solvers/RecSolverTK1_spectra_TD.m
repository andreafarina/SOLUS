%==========================================================================
% This function contains solvers for DOT or fDOT.
% Andrea Farina 04/17
%==========================================================================

function [bmua,bmus,bconc,bA,bB] = RecSolverTK1_spectra_TD(solver,grid,mua0,mus0, n, A,...
    Spos,Dpos,dmask, dt, nstep, twin, self_norm, data, irf, ref, sd, fwd_type,radiometry,spe,conc0,a0,b0)
%% Jacobain options
USEGPU = 0;%gpuDeviceCount;
LOAD_JACOBIAN = solver.prejacobian.load;      % Load a precomputed Jacobian
geom = 'semi-inf';
xtransf = '(x/x0)'; %log(x),x,log(x/x0)
type_ratio = 'gauss';   % 'gauss', 'born', 'rytov';
type_ref = 'theor';      % 'theor', 'meas', 'area'

% -------------------------------------------------------------------------
nQM = sum(dmask(:));
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
dummy_proj = zeros(size(twin,1),nQM);
idxmeas = findMeasIndex(dmask);
for inl = 1:radiometry.nL
    meas_set = idxmeas{inl};
    proj_single = proj(:,meas_set);
    proj_single = WindowTPSF(proj_single,twin(:,:,meas_set));
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
    load(solver.prejacobian.path);
    toc;
else
    tic;
    J = Jacobian ( mua0, mus0);
    [jpath,~, ~] = fileparts(solver.prejacobian.path);
    if ~exist(jpath,'dir')
        mkdir(jpath)
    end
    save(solver.prejacobian.path,'J','-v7.3');
    toc;
end

% sd normalisation
J = spdiags(1./sd(:),0,numel(sd),numel(sd)) * J;
nsol = size(J,2);

% parameter normalisation (scale x0)
if ~strcmpi(xtransf,'x')
    J = J * spdiags(x0,0,length(x0),length(x0));
end 
J(mask,:) = [];


%% Structured laplacian prior

siz_prior = size(solver.prior.refimage);
% %solver.prior.refimage(solver.prior.refimage == max(solver.prior(:))) = 1.1*min(solver.prior.refimage(:)); 
% %solver.prior.refimage = solver.prior.refimage .* (1 + 0.01*randn(size(solver.prior)));
%[L,~] = StructuredLaplacianPrior(solver.prior.refimage,siz_prior(1),siz_prior(2),siz_prior(3));
%% new gradient
k3d = Kappa(solver.prior.refimage,5);
[Dx,Dy,Dz] = gradientOperator(grid.dim,[1,1,1],[],'1st');
L = [sqrt(k3d)*Dx;sqrt(k3d)*Dy;sqrt(k3d)*Dz];

%% Solver
disp('Calculating singular values');
s = svds(J,1);

L1 = [];
for ip = 1:p
     L1 = blkdiag(L1,L);
end
%% case Lcurve or direct solution
b = [dphi;zeros(p*3*nsol/p,1)];
if USEGPU
    gpu = gpuDevice; %#ok<UNRCH>
    disp('Using GPU');
    b = gpuArray([full(dphi);zeros(p*3*nsol/p,1)]);
end

if numel(solver.tau)>1
    for i = 1:numel(solver.tau)
        alpha = solver.tau(i)*s(1); 
        disp(['Solving for tau = ',num2str(solver.tau(i))]);
        tic;
        if USEGPU
            A = [J;alpha*L1]; %#ok<UNRCH>
            A = gpuArray(sparse(A));
            dx = lsqr(A*1e10,b*1e10,1e-6,1000);
        else
            dx = lsqr([J;alpha*L1],b,1e-6,1000);
        end
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
if USEGPU > 0
    %J(abs(J)<1e-3) = 0;
    A = [J;alpha*L1];
    A = gpuArray(sparse(A));
    dx = lsqr(A*1e10,b*1e10,1e-6,1000);
    dx = gather(dx);
else
    dx = lsqr([J;alpha*L1],b,1e-6,1000);
end
%==========================================================================
%%                        Add update to solution
%==========================================================================

x = x + dx;

x = BackTransfX(x,x0,xtransf);

[bmua,bmus,bconc,bAB] = XtoMuaMus_spectral(x,mua0,mus0,type_jac,spe,conc0,a0,b0);
bA = bAB(:,1);bB = bAB(:,2);



end
