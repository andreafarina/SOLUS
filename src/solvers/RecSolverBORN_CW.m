%==========================================================================
% This function contains solvers for DOT or fDOT.
% To have available all the functionalty install REGU toolbox
% Andrea Farina 12/16
%==========================================================================

function [bmua,bmus] = RecSolverBORN_CW(~,grid,mua0,mus0, A,...
    Spos,Dpos,dmask,data,ref,~)
%% Jacobain options
SAVE_JACOBIAN = 1;      % Save the Jacobian into a file
LOAD_JACOBIAN = 0;      % Load a precomputed Jacobian
COMPUTE_JACOBIAN = 1;  % Compute the Jacobian
geom = 'semi-inf';
%% REGULARIZATION PARAMETER CRITERION
REGU = 'lcurve';        % 'lcurve' or 'gcv'
BACKSOLVER = 'Tikh';    % 'Tikh' or 'tsvd' or 'simon'

rdir = ['../results/test/precomputed_jacobians/'];
mkdir(rdir);
disp(['Intermediate results will be stored in: ' rdir])
%save([rdir 'REC'], '-v7.3','REC'); %save REC structure which describes the experiment
% -------------------------------------------------------------------------

bdim = (grid.dim);
nQM = sum(dmask(:));

Jacobian = @(mua, mus) JacobianCW (grid, Spos, Dpos, dmask, mua, mus, A, geom);
%% Inverse solver
proj = ForwardCW(grid, Spos, Dpos, dmask, ...
        mua0, mus0, [], [], A, geom, 'homo');
%% data scaling
%sd = proj(:);
sd = ones(size(proj));

if ref == 0
    ref = proj(:);
end
% figure(400);
% subplot(2,2,1),
% plot([data(:) ref(:) proj(:)./sum(proj(:))]),legend('data','ref','proj'),grid;
% subplot(2,2,2),
% plot(ref(:)./proj(:)),legend('ratio ref/proj'),
% subplot(2,2,3),grid;
% plot([ref(:)./sum(ref(:)) proj(:)./sum(proj(:))]),
% legend('ref','proj'),grid;
% subplot(2,2,4),
% plot([ref(:)./sum(ref(:))-proj(:)./sum(proj(:))]),
% legend('ref - proj'),grid;
% drawnow;

% solution vector
x = ones(grid.N,1) * mua0; 
x0 = x;
p = length(x);
dphi = (data(:)-ref(:));%./ref(:);
%dphi = log(data(:)) - log(ref(:));
%save('dphi','dphi');
% ---------------------- Construct the Jacobian ---------------------------
if COMPUTE_JACOBIAN == 1
    fprintf (1,'Calculating Jacobian\n');
    tic;
    J = Jacobian ( mua0, mus0);
    toc;
end
if LOAD_JACOBIAN == 1
    fprintf (1,'Loading Jacobian\n');
    tic;
    load([rdir,'Jacobian'])
    toc;
end
if SAVE_JACOBIAN == 1
    fprintf (1,'Saving Jacobian\n');
    tic;
    save([rdir,'Jacobian'],'J');
    toc;
end

J = spdiags(1./sd,0,nQM,nQM) * J;  % data normalisation
nsol = size(J,2);
%   parameter normalisation (map to log)
%     for i = 1:p
%         J(:,i) = J(:,i) * x(i);
%     end
% ------- to solve only for mua uncomment the following sentence ----------
%J(:,nsol+(1:nsol)) = 0;


%% Solver 
%alpha_vec = logspace(-6,-3,10);        % Regularization parameter
%close all
alpha_vec = 0;%1e-10;        % Regularization parameter

if ~strcmpi((BACKSOLVER),'simon')
[U,s,V]=csvd(J);     % compact SVD (Regu toolbox)
        figure(402);
        picard(U,s,dphi);    % Picard plot (Regu toolbox)    
end
for i = 1:numel(alpha_vec)
alpha = alpha_vec(i)
if ~exist('alpha','var')
    figure(403);
    if strcmpi(REGU,'lcurve')
       alpha = l_curve(U,s,dphi)%,BACKSOLVER);  % L-curve (Regu toolbox)
    elseif strcmpi(REGU,'gcv')
        alpha = gcv(U,s,dphi,BACKSOLVER);
    end
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
        alpha = 1e0*s1
        dx = [J;sqrt(alpha)*speye(nsol)]\[dphi;zeros(nsol,1)];
        %dx = [dx;zeros(nsol,1)];
end   
    
%rho
%cond(J)
%dx = [dx;zeros(nsol,1)];
%dx = [J;alpha*eye(2*nsol)]\[dphi;zeros(2*nsol,1)];
%dx = lsqr(J,dphi,1e-4,100);
%dx = pcg(@(x)JTJH(x,J,eye(2*nsol),alpha),J'*dphi,1e-4,100);
%dx = gmres(@(x)JTJH(x,J,eye(2*nsol),alpha),J'*dphi,[],1e-4,100);

%dx = J' * ((J*J' + alpha*eye(length(data)))\dphi);
  
%==========================================================================
%%                        Add update to solution
%==========================================================================
x = x + dx;
%logx = logx + dx;
%x = exp(logx);



bmua = x;
bmus = ones(size(bmua)) * mus0;

% display the reconstructions
figure(304);
%figure('Position',get(0,'ScreenSize'));
ShowRecResults(grid,reshape(bmua,bdim),grid.z1,grid.z2,grid.dz,1,...
    min(bmua),max(bmua));
suptitle('Recon Mua');
%export_fig '../Results/20151111/rec_mua_1incl.pdf' -transparent

figure(305);
%figure('Position',get(0,'ScreenSize'));
ShowRecResults(grid,reshape(bmus,bdim),grid.z1,grid.z2,grid.dz,1,...
    min(bmus),max(bmus));
suptitle('Recon Mus');
%export_fig '../Results/20151111/rec_mus_1incl.pdf' -transparent
drawnow;
tilefigs;
%pause
x = x0;
%save([rdir 'sol'], 'xiter');
%save([rdir 'err_pattern'],'erri');
end
end