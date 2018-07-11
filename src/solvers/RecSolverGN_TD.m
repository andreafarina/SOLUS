%==========================================================================
% This function contains solvers for DOT or fDOT.
% To have available all the functionalty install REGU toolbox
% Andrea Farina 05/15
%==========================================================================

function [bmua,bmus,erri] = RecSolverGN_TD(solver,grid,hMesh,hBasis,mua0,mus0,refind,...
    qvec,mvec,dmask,dt,nstep,twin,self_norm,data,irf,ref,sd,verbosity)

%global xd
step = 0.1;
LOAD_JACOBIAN = solver.prejacobian.load;
LOGX = false;
solver.Himplicit = true;
solver.itrmax = 3;
solver.tol = 1e-6;
solver.tolKrylov = 1e-6;
prm.method = 'TV';


% -------------------------------------------------------------------------
[~,type_jac] = ExtractVariables(solver.variables);
Jacobian = @(mua, mus) JacobianTD (grid, [], [], dmask, mua, mus, refind, [], ...
    dt, nstep, twin, irf, [],type_jac,'fem');
%% self normalise (probably useless because input data are normalize)
nQM = sum(dmask(:));
% if self_norm == true
%         data = data * spdiags(1./sum(data)',0,nQM,nQM);
%         ref = ref * spdiags(1./sum(ref)',0,nQM,nQM);
% end
N = hMesh.NodeCount;

% Set up homogeneous initial parameter estimates
mua = ones(N,1) * mua0;                           % initial mua estimate
mus = ones(N,1) * mus0;                         % initial mus estimate
n = ones(N,1) * refind;                           % refractive index estimate
kap = 1./(3*mus);                          % diffusion coefficient
solmask = find(hBasis.GridElref>0);
slen = numel(solmask);
%% proj[x0] %%
[proj, Aproj] = ForwardTD([],[], [], dmask, mua0, mus0, n, ...
                [],[], [], dt, nstep, self_norm,...
                [], 'fem');
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

%% prepareData
if self_norm == true
    strfactor = 'factor_one';
    strsd = 'poisson';
else
    strfactor = 'factor_td';
    strsd = 'Poisson';
end
[data,sd,ref,sdj,proj,mask] = PrepDataLH(data,sd,ref,proj,...
    'nonlin',strfactor,strsd);

   
% map initial parameter estimates to solution basis
bdim = (hBasis.Dims())';
bmua = hBasis.Map('M->B', mua);           % mua mapped to full grid
bmus = hBasis.Map('M->B', mus);           % mus mapped to full grid
bkap = hBasis.Map('M->B', kap);           % kap mapped to full grid
% mask non valid data
%ref(mask) = [];
%data(mask) = [];

% solution vector
x0 = PrepareX0([mua0,1./(3*mus0)],slen,type_jac);
x = x0./x0;%ones(size(x0)); 
if LOGX == true
    logx =log(x);%x./x0;%log(x);            % transform to log
else
    logx = x;
end                                      % transform to log
nsol = length(x);

% REGULARIZATION STRUCTURE
prm.basis = hBasis;
prm.x0 = logx;
prm.tau = solver.tau;  
%prm.tv.beta = lprm.beta;
if isfield(solver.prior,'refimage')
if ~isempty(solver.prior.path)
    ref_img = solver.prior.refimage;
    prm.prior.refimg = [ref_img(:);ref_img(:)];
    prm.prior.threshold = 0.25;
    prm.prior.smooth = 0.1;
end
end
prm.tau = 1; %no tau is set
hreg = toastRegul(prm,logx);
% -------------------- Initial data error (=2 due to data scaling) --------
err0 = toastObjective (proj(~mask), data(~mask), sd(~mask), hreg, logx);  %initial error
err = err0;                                         % current error
errp = inf;%1e10;                                         % previous error
erri(1) = err0;
itr = 1; % iteration counter
fprintf (1, '\n**** INITIAL ERROR %e\n\n', err);
tic;

% set figures for visualization
habs = figure;
hsca = figure;
% --------------------------- Exit condition ------------------------------
while (itr <= solver.itrmax) && (err > solver.tol*err0)...
        && (errp-err > solver.tol)
    
    % -------------------------------------------------------------------------
    errp = err;
    % ---------------------- Construct the Jacobian ---------------------------
    nwin = size(twin,1);
    
    
    if LOAD_JACOBIAN == true
        if (~exist('J','var'))
            fprintf (1,'Loading Jacobian\n');
            tic;
            load(solver.prejacobian.path);
            toc;
        end
        LOAD_JACOBIAN = false;
        [jpath,~,~] = fileparts(solver.prejacobian.path);
        if ~exist(jpath,'dir')
            mkdir(jpath)
        end
        save(solver.prejacobian.path,'J');
        toc;
    else
        fprintf (1,'Calculating Jacobian\n');
        tj = tic;
        J = Jacobian (mua, mus);
        disp(['Jacobian computation time: ',num2str(toc(tj))]);
        %save([jacdir,jacfile,'_it',num2str(itr)],'J');
        %load Jpoint
        if itr == 1
            [jpath,~,~] = fileparts(solver.prejacobian.path);
            if ~exist(jpath,'dir')
                mkdir(jpath)
            end
            save(solver.prejacobian.path,'J');
        end
        %J(:,nsol+(1:nsol)) = 0;
    end
    
       
    %% Normalized measurements
    if self_norm == true
        for i=1:nQM
            sJ = sum(J((1:nwin)+(i-1)*nwin,:));
            sJ = repmat(sJ,nwin,1);
            sJ = spdiags(proj((1:nwin)+(i-1)*nwin),0,nwin,nwin) * sJ;
            J((1:nwin)+(i-1)*nwin,:) = (J((1:nwin)+(i-1)*nwin,:) - sJ)./Aproj(i);
        end
    end
   
    %   parameter normalisation (scale x0)
       J = J * spdiags(x0,0,length(x0),length(x0));
 
    % % Transform to   df/ d logx
    if LOGX == true
       J =  J * spdiags(x, 0,length(x0),length(x0)) ;
    end
    
    % ----------------------- Data normalization ------------------------------
    %
    for i = 1:numel(proj)
        J(i,:) = J(i,:) / sdj(i);
    end  
    
    J = spdiags(1./sd,0,numel(data),numel(data)) * J;
    J(mask,:) = [];
   
    
   % Normalisation of Hessian (map to diagonal 1)
    %psiHdiag = hreg.HDiag(logx);
    %M = zeros(nsol,1);
%     for i = 1:p
%         M(i) = sum(J(:,i) .* J(:,i));
%         M(i) = M(i) + psiHdiag(i);
%         M(i) = 1 ./ sqrt(M(i));
%     end
    M = ones(nsol,1);
    for i = 1:nsol
        J(:,i) = J(:,i) * M(i);
    end

 %M = ones(size(M));
    % Gradient of cost function
    tau = solver.tau * max(svd(J))
    %tau = solver.tau * sum(sum(J.^2))
    r = J' * (2*(data(~mask)-proj(~mask))./sd(~mask));
    r = r - tau * hreg.Gradient (logx) .* M;
    
    if solver.Himplicit == true
        % Update with implicit Krylov solver
        fprintf (1, 'Entering Krylov solver\n');
        %dx = toastKrylov (x, J, r, M, 0, hreg, lprm.tolKrylov);
        if exist('hreg','var')
            HessReg = hreg.Hess (logx);
        end
        dx = krylov(r);
        disp(['|dx| = ',num2str(norm(dx))])
    else
        % Update with explicit Hessian
        H = J' * J;
        lambda = 0.01;
        H = H + eye(size(H)).* lambda;
        dx = H \ r;
        clear H;
    end
    
    %clear J;
%             %% AF   TBU is it needed??
%     for i = 1:p
%         dx(i) = dx(i) ./ M(i);
%     end   
    % Line search
    fprintf (1, 'Entering line search\n');
    step0 = step;
    %step0 = .1;
    [step, err] = toastLineSearch (logx, dx, step0, err, @objective, 'verbose', 1);
    if errp-err <= solver.tol
        dx = r; % try steepest descent
        fprintf (1, 'Try steepest descent\n');
        step = toastLineSearch (logx, dx, step0, err, @objective, 'verbose', 1);
    end

    % Add update to solution
    logx = logx + dx*step;
    if LOGX == true
        x = exp(logx);
    else
        x = logx;
    end
    x = x.*x0;
    % Map parameters back to mesh
    [smua,smus] = XtoMuaMus(x,mua0,mus0,type_jac);
    
    mua = hBasis.Map('S->M',smua);
    bmua = hBasis.Map('S->B',smua);
    mus = hBasis.Map('S->M', smus);
    bmus = hBasis.Map('S->B',smus);
        
    % display the reconstructions
    
    figure(habs);
    ShowRecResults(grid,reshape(bmua,bdim),...
        grid.z1,grid.z2,grid.dz,1,'auto',0.00,0.05);
    %suptitle('Recon Mua');
    %save_figure([rdir, 'mua_rec_it_',num2str(itr)]);
    %saveas(gcf,[rdir,lprm.filename,'_mua_rec_it_',num2str(itr)],'tif');
    
    figure(hsca);
    ShowRecResults(grid,reshape(bmus,bdim),...
        grid.z1,grid.z2,grid.dz,1,'auto');
    %suptitle('Recon Mus');
    %save_figure([rdir, 'mus_rec_it_',num2str(itr)]);
    %saveas(gcf,[rdir,lprm.filename,'_mus_rec_it_',num2str(itr)],'tif');
    drawnow;
    tilefigs;
    
    %==========================================================================
    %%                        Update field projections
    %==========================================================================
    disp('----- Update field projections -------');tic;
    
    [proj, Aproj] = ForwardTD([],[], [], dmask, mua, mus, n, ...
                [],[], [], dt, nstep, self_norm,...
                [], 'fem');
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
    proj = proj(:);%-proj0(:);% * factor;
    %proj(mask) = [];
    %proj = log(proj(:)) - log(proj0(:));
    %sdj = proj(:);  
    %==========================================================================
    %%                        Update objective function
    %==========================================================================
    err = toastObjective (proj(~mask), data(~mask), sd(~mask), hreg, logx);
    fprintf (1, '**** GN ITERATION %d, ERROR %e\n\n', itr, err);
    
    erri(itr) = err;
        bmua_itr(itr,:) = bmua;
        bmus_itr(itr,:) = bmus;
    %   REC.Data = REC.Data./REC.Data(1)*REC.proj(1);
    % REC.Data = REC.Data./max(REC.Data)*max(REC.proj);
%     it_dir = [rdir,filesep,solver.filename,filesep];
%     mkdir(it_dir);
%     save([it_dir,'res_',num2str(itr)],'bmua_itr','bmus_itr');
    itr = itr+1;
end
%save([rdir,'sol'],'bmua_itr','bmus_itr');

% =====================================================================
    % Callback function for objective evaluation (called by toastLineSearch)
    function p = objective(x)
    if LOGX == true
        xx = exp(x);
    else
        xx = x;
    end
    xx = x.*x0;
    if strcmpi(type_jac,'mua')
        mua = hBasis.Map('S->M',xx);
    else
        mua = hBasis.Map('S->M',xx(1:slen));
        mus = hBasis.Map('S->M',1./(3*xx(slen+1:end)));
    end
    %mua(mua<1e-4) = 1e-4;
    %mus(mus<0.2) = 0.2;
    %% if fixed scattering
   % mus(:) = mus0;
%     for j = 1:length(mua) % ensure positivity
%         mua(j) = max(1e-4,mua(j));
%         mus(j) = max(0.2,mus(j));
%      end
%  
    [proj, Aproj] = ForwardTD([],[], [], dmask, mua, mus, n, ...
                [],[], [], dt, nstep, self_norm,...
                [], 'fem');
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
    proj = proj(:);%-proj0(:);% * factor;
    %proj (mask) = [];
    %proj = log(proj(:)) - log(proj0(:));
    %sdj = proj(:);
    %proj = privProject (hMesh, hBasis, x, ref, freq, qvec, mvec);
    [p, p_data, p_prior] = toastObjective (proj(~mask), data(~mask), sd(~mask), hreg, x);
   % if verbosity > 0
        fprintf (1, '    [LH: %e, PR: %e]\n', p_data, p_prior);
   % end
    end

%save([rdir 'sol'], 'xiter');
%save([rdir 'err_pattern'],'erri');

% =====================================================================
    % Krylov solver subroutine
    function dx = krylov(r)
        k_t = cputime;
        %switch prm.solver.krylov.method
        %    case 'bicgstab'
        %        [dx, k_flag, k_res, k_iter] = bicgstab(@jtjx, r, lprm.tolKrylov,100);
       %     otherwise
       %         [dx, k_flag, k_res, k_iter] = gmres (@jtjx, r, 30,lprm.tolKrylov, 100);
                [dx, k_flag, k_res, k_iter] = pcg (@jtjx, r, solver.tolKrylov, 100);
       % end
        k_dt = cputime-k_t;
        % fprintf (1, '--> iter=%0.0f(%0.0f), time=%0.1fs, res=%g\n', ...
       %     k_iter(1), k_iter(2), k_dt, k_res);
         fprintf (1, '--> iter=%0.0f, time=%0.1fs, res=%g\n', ...
            k_iter(1), k_dt, k_res);
       clear k_t k_dt k_flag k_res k_iter
    end
    
    
    % =====================================================================
    % Callback function for matrix-vector product (called by toastKrylov)
    function b = jtjx(x)
        b = J' * (J*x);
        if exist('hreg','var')
            b = b + M .* (tau * HessReg * (M .* x));
        end
    end


end