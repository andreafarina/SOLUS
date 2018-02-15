function [proj, Area, proj_lifetime] = ProjectFieldTD(hMesh,qvec,mvec,dmask,...
    mua,musp,conc,tau,n,dt,nstep,freq,self_norm,fwd,verbosity)
% PROJECTFIELDTD Time Domain forward projection
% function [proj, Area, proj_lifetime] = ProjectFieldTD(hMesh,qvec,mvec,dmask,...
%            mua,musp,conc,tau,n,dt,nstep,freq,self_norm,fwd,verbosity)
% March 2017: implemented GPU iterative solver
USEGPU = 1;
tol = 1e-8;

spacer = ' ------- ';
disp([spacer mfilename ': setup' spacer]);
tic;

% Setting parameters
N = hMesh.NodeCount;
nQ = size(dmask,2);       % effective sources
nM = size(dmask,1);       % effective detectors
nQM = sum(dmask(:));      % effective measurements
mvecT = mvec';
proj = zeros(nQM,nstep);
% initial condition
%phi = qvec;

% Time stepping scheme: 0 - Forward Euler, 1 - Backward Euler, 1/2 - Crank Nicholson
theta = 1;

%[smat,bmat] = toastSysmat (hMesh, mua, mus, ref, 0);
Smat = real(dotSysmat(hMesh, mua, musp, n, freq));
mmat = hMesh.Massmat; % the bmat is only the boundary part!

K0 = -(Smat * (1-theta) - mmat * 1/dt);            % backward difference matrix
K1 = Smat * theta + mmat * 1/dt;                   % forward difference matrix

%dK1 = sparse(diag(diag(K1)));

%Lc = ichol(K1);
%hLc = @(x) Lc'\(Lc\x);
%hPK1 = @(x) Lc\(K1\(Lc'\x));

dstep = ceil(nstep/4);
    
%% CPU-based version
if ~USEGPU
    tic
    disp('---computing LU decomp---');
    [L,U,P,Q,R] = lu(K1);
    toc
    disp([spacer mfilename ': time stepping' spacer]);
    tstart = cputime;
    % AF initial conditions see toast examples fwd_tpsf.m
    % initial condition
    q = qvec/dt;
    %q = zeros(size(phi));
    phi = Q*(U\(L\(P*(R\q))));
    tmp = mvecT * phi;
    proj(:,1) = tmp(dmask(:));
    for i = 2:nstep
        q = K0 * phi;
        %prec - standard lu
        %phi = U\(L\q);
        %prec umfpack: P*(R\A)*Q = L*U
        phi = Q*(U\(L\(P*(R\q))));
        tmp = mvecT * phi;
        proj(:,i) = tmp(dmask(:));
        ddisp(['Step: ' num2str(i) '. Elapsed time: ' num2str(cputime-tstart)], mod(i,dstep)==0, mfilename);
    end
    proj = permute(proj,[2 1]);
%% GPU-based version with iterative solver
else
    %Select GPU device (this shoud be done on higher level)
    gpu = gpuDevice;%(USEGPU);
    
    K0 = gpuArray((K0));
    K1 = gpuArray((K1));
    %dK1 = gpuArray(dK1);
    mvecT = gpuArray((mvecT));
    proj = gpuArray(proj);
    
    tstart = cputime;
    % AF initial conditions see toast examples fwd_tpsf.m
    % initial condition
    q = qvec/dt;
    phi = gpuArray(zeros(size(q)));
    for iq = 1:nQ
        qg = gpuArray((q(:,iq)));
        %[phi(:,iq),flag] = bicgstab(K1,q(:,iq),1e-12,100);%,U);%,PP',PP);
        [phi(:,iq),~] = pcg(K1,qg,tol,100);
        %[phi(:,iq),flag] = gmres(K1,q(:,iq),30,1e-12,100);
    end
    tmp = mvecT * phi;
    %tmp = gather(tmp);
    proj(:,1) = tmp(dmask(:));
    
    
    for i = 2:nstep
        q = K0 * phi;
        for iq = 1:nQ
            qg = gpuArray(q(:,iq));
            %[phi(:,iq),flag] = bicgstab(K1,q(:,iq),1e-12,100);%,L);%,U);%,PP',PP);
            [phi(:,iq),~] = pcg(K1,qg,tol,100);
            %[phi(:,iq),flag] = gmres(K1,q(:,iq),30,1e-12,100);
        end
        
        tmp = mvecT * phi;
        proj(:,i) = tmp(dmask(:));
        ddisp(['Step: ' num2str(i) '. Elapsed time: ' num2str(cputime-tstart)], mod(i,dstep)==0, mfilename);
    end
    proj = permute(proj,[2 1]);
    
    proj = gather(proj);
    % Reset memory on GPU
    gpu.delete;
    %reset(gpu);
end

%% Normalized measurement
Area = sum(proj);
if self_norm == true
    proj = proj * spdiags(1./Area',0,nQM,nQM);
end
%proj = proj * diag(1./Area);

%% Fluorescence
% to be completed.


