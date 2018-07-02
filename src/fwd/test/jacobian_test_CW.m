% MATLAB-TOAST sample script:
% 2D image reconstruction with nonlinear conjugate gradient solver

clear all
close all
% ======================================================================
% User-defined parameters
% ======================================================================
meshname  = 'circle25_32.msh';              % mesh file
qmname    = 'circle25_1x1.qm';            % source-detector file
refind = 1.4;                           % refractive index
dimRect = [64 32];                          % solution basis: grid dimension
grd = round(dimRect);
freq = 0;                             % modulation frequency [MHz]
test_meshbasis = false;
% ======================================================================
% End user-defined parameters
% ======================================================================

% Initialisations
tol = 1e-12;

disp('jacobian_test');
disp('Test the integrity of toastJacobian by comparing to explicit');
disp('calculation obtained by pixel-wise parameter perturbation.');
% Set up some variables

c0 = 0.299;
cm = c0/refind;

% time domain parameters
dt = 1;                             % time step in picoseconds
nstep = 1000;
a = 1:1:nstep;
b = a + 99;
twin = [a', a'];
%twin = [ a',b'];

% Read a TOAST mesh definition from file.
%% Make rectangle
[vtx,idx,eltp] = mkrect ([[0 0];dimRect], dimRect);
hMesh = toastMesh(vtx,idx,eltp);
%hMesh = toastMesh(meshname);

hMesh.SetQM([36,0],[28,32]);
%hMesh.ReadQM (qmname);
n = hMesh.NodeCount();
dmask = hMesh.DataLinkList();

% Set up the mapper between FEM and solution bases
BB = hMesh.BoundingBox();
hBasis = toastBasis (hMesh,grd);
blen = hBasis.blen;
solmask = find(hBasis.GridElref>0);
slen = numel(solmask);
% element size

bmin = BB(1,:); bmax = BB(2,:);
bdim = hBasis.Dims();bdim = bdim';
vscale = prod((bmax(:)-bmin(:))./(bdim(:)));%1;%2.3;

dd = (bmax-bmin)./(bdim);    % the basis can be smaller than the mesh! attention to z
ElSize = hMesh.ElementSize;

% Set up homogeneous initial parameter estimates
blen = hBasis.blen;
bmua0 = ones(blen,1) * 0.01;
bmus0 = ones(blen,1) * 1;
mua0 = hBasis.Map('B->M',bmua0);
mus0 = hBasis.Map('B->M',bmus0);
ref = ones(n,1) * refind;

% Generate source vectors
qvec = real(hMesh.Qvec('Neumann', 'Gaussian', 2));

% Generate measurement vectors
A = toastDotBndterm(ref,'Contini');
mvec = real(hMesh.Mvec ('Gaussian', 2,0).*ones(n,1)./(2*A));
% SINUSIODAL Q / M
    Area = dimRect(1);
    x = 1:Area;
    N = 16;
    p = sin(2*pi*x./N);
    qvec = zeros(1,n);
    qbasis = zeros(bdim);
    qbasis(:,1) = 1/2 + 1/2*p;
    qvec = hBasis.Map('B->M',qbasis);
    mvec = zeros(1,n);
    mbasis = zeros(bdim);
    mbasis(:,end) = 1/2 + 1/2*p;
    mvec = hBasis.Map('B->M',mbasis);
    
    
% Initial data set f[x0]
proj0 = ProjectFieldCW(hMesh,qvec,mvec,dmask,...
            mua0,mus0,[],ref,freq,'diff',0);

dmua = max(mua0)*1e-3;
bkap0 = 1./(3*bmus0);
k0 = hBasis.Map('B->M',bkap0);
dkap = max(k0)*1e-3;
irf = 0;

%dmask = true;
GRADIENT = 'toast';
%% Test 1: Jacobian in mesh basis
if test_meshbasis == true
tic;  
mesh_el_size = ElSize;
%load Je_mesh

J = toastJacobianTimedomain_(hMesh,0,qvec,mvec, dmask, mua0, mus0, ref,dt,nstep,twin,irf);
%J = -toastJacobianTimedomain_conv_time(hMesh,0,qvec,mvec, dmask, mua0, mus0, ref,dt,nstep,twin);
J = toastJacobianTimedomain_INPROGRESS(hMesh,hBasis,qvec,mvec, ...
    dmask, mua0, mus0, ref,dt,nstep,twin,irf, GRADIENT);
toc;
% h = waitbar(0,'Calculating explicit Jacobian');
% 
% for i=1:n
%     mua = mua0;
%     mua(i) = mua0(i)+dmua;
%     proj = ProjectFieldTD(hMesh,qvec,mvec,dmask,...
%             mua,mus0,[],[],ref,dt,nstep,freq,0,'diff',0);
%     dy = (proj-proj0)/dmua;
%     Je(:,i) = dy / cm; % parameter is c*mua
%     waitbar(i/n);
% end
% delete(h);
%load Je
%%
for i = 1:5:nstep
mua_d = J(i,1:n) * dt;
mua_e = Je(i,1:n);
mn = min([mua_d mua_e]);
mx = max([mua_d mua_e]);

figure(1);
img = reshape(toastMapMeshToBasis(hBasis,mua_d),grd);
subplot(1,3,1);imagesc(img); axis equal tight; title('direct'); colorbar
%subplot(1,3,1);imagesc(img); axis equal tight; title('direct'); colorbar

img = reshape(toastMapMeshToBasis(hBasis,mua_e),grd);
subplot(1,3,2);imagesc(img); axis equal tight; title('explicit'); colorbar
%subplot(1,3,2);imagesc(img); axis equal tight; title('explicit'); colorbar

ratio = mua_d./mua_e;
img = reshape(toastMapMeshToBasis(hBasis,ratio),grd);
subplot(1,3,3);imagesc(img);axis equal tight;colorbar;title('ratio')
suptitle('J_e_MESH')
mean(mean(img(8:24,8:24)));
pause(0.01)
end
save('Je_mesh','J','Je');
clear J Je
end

%% Test 2: Jacobian in grid basis
tic;
disp('Jacobian in grid basis')
tic;
J = toastJacobianCW2(hMesh,hBasis,qvec,mvec, ...
    dmask, mua0, mus0, ref, GRADIENT);
toc;
%save('J_CW2','J');
%load J_CW
h = waitbar(0,'Calculating explicit Jacobian');
tic;
for i=1:slen
    idx = solmask(i);
    bkap = bkap0;
    bmua = bmua0;
    bmua(idx) = bmua0(idx)+dmua;
    mua = hBasis.Map('B->M',bmua);
    proj = ProjectFieldCW(hMesh,qvec,mvec,dmask,...
            mua,mus0,[],ref,freq,'diff',0);
    dy = (proj-proj0)/dmua;
    Je(:,i) = dy;
    bkap(idx) = bkap0(idx) + dkap;
    kap = hBasis.Map('B->M',bkap);
    mus = 1./(3*kap);
    proj = ProjectFieldCW(hMesh,qvec,mvec,dmask,...
            mua0,mus,[],ref,freq,'diff',0);
    dy = (proj-proj0)/dkap;
    Je(:,i+slen) = dy;% / cm;
    waitbar(i/slen);
end
toc
delete(h);
%save('Je_CW2','Je');
%load Je_CW2
%Je = Je;
%% plot
for i = 1:1
mua_d = -J(i,1:slen) * vscale;% * dt;% / cm;
mua_e = -Je(i,1:slen);
mus_d = -J(i,slen+1:end) * vscale;
mus_e = -Je(i,slen+1:end);
%mua_c = Jconv(i,1:slen);% * vscale  * dt / cm;
mn = min([mua_d mua_e]);% mua_c]);
mx = max([mua_d mua_e]);% mua_c]);
mn_e = min(mua_e(:));
mx_e = max(mua_e(:));
% [bbmin bbmax] = toastMeshBB(hMesh);
% elsize = prod ((bbmax-bbmin) ./ grd');
% mua_d = mua_d * elsize; % scale with element size
% plot
figure(2);
img_d = reshape(hBasis.Map('S->B',mua_d),grd);
%subplot(1,3,1);imagesc(img,[mn mx]); axis equal tight; title('direct'); colorbar
subplot(2,2,1);imagesc(img_d,[mn mx]); axis equal tight; title('direct Fourier'); colorbar

img_e = reshape(hBasis.Map('S->B',mua_e),grd);
%subplot(1,3,2);imagesc(img,[mn mx]); axis equal tight; title('explicit'); colorbar
subplot(2,2,2);imagesc(img_e,[mn mx]); axis equal tight; title('explicit'); colorbar

ratio = img_d./img_e;
%img = reshape(hBasis.Map,ratio),grd);
subplot(2,2,3);imagesc(ratio);axis equal tight;colorbar;title('ratio Fourier/explicit')

figure(3);
mn = min([mus_d mus_e]);% mua_c]);
mx = max([mus_d mus_e]);% mua_c]);
img_d = reshape(hBasis.Map('S->B',mus_d),grd);
%subplot(1,3,1);imagesc(img,[mn mx]); axis equal tight; title('direct'); colorbar
subplot(2,2,1);imagesc(img_d,[mn mx]); axis equal tight; title('direct Fourier'); colorbar

img_e = reshape(hBasis.Map('S->B',mus_e),grd);
%subplot(1,3,2);imagesc(img,[mn mx]); axis equal tight; title('explicit'); colorbar
subplot(2,2,2);imagesc(img_e,[mn mx]); axis equal tight; title('explicit'); colorbar

ratio = img_d./img_e;
%img = reshape(hBasis.Map,ratio),grd);
subplot(2,2,3);imagesc(ratio);axis equal tight;colorbar;title('ratio Fourier/explicit')


% img_c = reshape(hBasis.Map('S->M',mua_c),grd);
% subplot(2,3,4);imagesc(img_c); axis equal tight; title('direct Convolution'); colorbar
% 
% ratio2 = img_d./img_c;
% %img = reshape(toastMapSolToBasis(hBasis,ratio2),grd);
% subplot(2,3,5);imagesc(ratio2);axis equal tight;colorbar;title('ratio Fourier/Conv')
% 
% ratio3 = img_c./img_e;
% img = reshape(toastMapSolToBasis(hBasis,ratio3),grd);
% subplot(2,3,6);imagesc(img);axis equal tight;colorbar;title('ratio Conv/explicit')
% suptitle('J_e_BASIS')
% mean (ratio);
pause(0.1);
end
% for i=1:5:nstep
% figure(3); % scattering
% mus_d = J(i,slen+(1:slen));% * vscale  * dt / cm;
% mus_e = Je(i,1:slen);
% mus_c = Jconv(i,slen+(1:slen));% * vscale  * dt / cm;
% mn = min([mus_d mua_e mus_c]);
% mx = max([mus_d mua_e mus_c]);
% 
% img_d = reshape(toastMapSolToBasis(hBasis,mus_d),grd);
% %subplot(1,3,1);imagesc(img,[mn mx]); axis equal tight; title('direct'); colorbar
% subplot(1,3,1);imagesc(img_d); axis equal tight; title('direct Fourier'); colorbar
% 
% img_c = reshape(toastMapSolToBasis(hBasis,mus_c),grd);
% %subplot(1,3,2);imagesc(img,[mn mx]); axis equal tight; title('explicit'); colorbar
% subplot(1,3,2);imagesc(img_c); axis equal tight; title('direct Conv'); colorbar
% 
% ratio = img_d./img_c;
% %img = reshape(toastMapSolToBasis(hBasis,ratio),grd);
% subplot(1,3,3);imagesc(ratio,[0.9 1.3]);axis equal tight;colorbar;title('ratio Fourier/Conv')
% suptitle('J_basis')
% pause(0.001)
% end
%save('Je_basis','J','Je');
%toastDeleteBasis(hBasis);
%toastDeleteMesh(hMesh);
