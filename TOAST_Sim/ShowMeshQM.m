close all
%file = 'SIMULATIONMask3D_Mask_benign_1.mat';
%file = 'dt5_muao0,0036_muai0,0091_muso1,26_musi1,50_benign_3.mat';
file = 'mus1_05-1_56_malignant_3_FwdTeo.mat';
%% save path
res_path = './';
out_folder = 'figure/';
SAVE = 0;
TD = 0;
%% Choose grid
GRID = 'native';%'native'; 'defined';
%% load
load(file);
%% create TOAST mesh
hmesh = toastMesh('slab_test.msh');
hmesh.Make(test.prm.mesh.toastmeshdata.vtx,...
    test.prm.mesh.toastmeshdata.idx,...
    test.prm.mesh.toastmeshdata.eltp);
%% create Q/M
ds = test.prm.mqwidth;
ds = 2;
hmesh.SetQM(test.prm.mesh.QM.Qpos,test.prm.mesh.QM.Mpos);
Q = hmesh.Qvec('Neumann','Gaussian',ds);
M = hmesh.Mvec('Gaussian',ds,0);
%% plot
figure,hmesh.Display(sum(Q,2));
switch lower(GRID)
    case 'native'
        nx = test.prm.mesh.nx+1;
        ny = test.prm.mesh.ny+1;
        nz = test.prm.mesh.nz+1;
        %% axis
        xx = linspace(test.prm.mesh.x_dim(1),test.prm.mesh.x_dim(2),...
        nx);
        yy = linspace(test.prm.mesh.y_dim(1),test.prm.mesh.y_dim(2),...
        ny);
        zz = linspace(test.prm.mesh.z_dim(1),test.prm.mesh.z_dim(2),...
        nz);
        hbasis = toastBasis(hmesh,[nx,ny,nz]);
    case 'defined'
        Init_DOT;
        xx = DOT.grid.x1:DOT.grid.dx:DOT.grid.x2;
        yy = DOT.grid.y1:DOT.grid.dy:DOT.grid.y2;
        zz = DOT.grid.z1:-DOT.grid.dz:-DOT.grid.z2;
        nx = numel(xx);ny = numel(yy);nz = numel(zz);
        
        hbasis = toastBasis(hmesh,[nx,ny,nz],[nx,ny,nz],...
                            [DOT.grid.x1, DOT.grid.x2; ...
                             DOT.grid.y1,DOT.grid.y2; ...
                             zz(1),zz(end)]);
end
        
%% set basis
Mua = hbasis.Map('M->B',test.prm.toastmeshdata.Mua);
Mua = reshape(Mua,hbasis.Dims');
Qb = zeros(numel(Mua),1);
for i = 1:size(Q,2)
    Qb = Qb + hbasis.Map('M->B',Q(:,i));
end
Qb = reshape(Qb,hbasis.Dims');
figure(2),imagesc(yy,xx,Qb(:,:,1)),
set(gca,'ydir','normal'),xlabel('y(mm)'),ylabel('x(mm)')

for i = 1:nz
    figure(3),imagesc(yy,xx,Mua(:,:,i)),
    set(gca,'ydir','normal'),xlabel('y(mm)'),ylabel('x(mm)')
    pause(0.1)
end
%% plot Q and M
figure,set(gcf,'Position',get(0,'ScreenSize'))
plot3(test.prm.mesh.QM.Qpos(:,1),test.prm.mesh.QM.Qpos(:,2),test.prm.mesh.QM.Qpos(:,3),'r*','MarkerSize',20),grid,
xlabel('x'),ylabel('y'),zlabel('z'),hold on
%% plot detectors
%plot3(test.prm.mesh.QM.Mpos(:,1),test.prm.mesh.QM.Mpos(:,2),test.prm.mesh.QM.Mpos(:,3),'bx','MarkerSize',20),
%xlabel('x(mm)'),ylabel('y(mm)'),zlabel('z(mm)')

%% plot heterogeneity
RES_FACT = 5;
Muaplot = imresizen(Mua,RES_FACT);
xxx = imresizen(xx',RES_FACT);
yyy = imresizen(yy',RES_FACT);
zzz = imresizen(zz',RES_FACT);
[x,y,z] = ind2sub(size(Muaplot),find(Muaplot>1/2*max(Muaplot(:))));
plot3(xxx(x), yyy(y), zzz(z), 'k-','MarkerSize',15);
%set(gca,'zdir','reverse'),
axis equal,
xlim([xx(1) xx(end)]),...
    ylim([yy(1) yy(end)]),...
    zlim([min(zz) max(zz)]),
% %legend('source','detector','pert');
%vol3d('Cdata',Muaplot,'XData',[xxx(1) xxx(end)],...)
%'YData',[yyy(1) yyy(end)],'ZData',[zzz(1) zzz(end)])
set(gca,'FontSize',20);
if SAVE == 1
    save_path = [res_path,out_folder,file(1:end-11)];
    save_figure(save_path);
end
save(['Mask3d_',file(1:end-11)],'Mua');

%hmesh.delete;
if TD == 1
%% test time-domain
dmask = true(size(M,2),size(Q,2));
%dmask(1,6) = true;
mua0 = test.prm.mua.out;
mus0 = test.prm.mus.out;
n_node = hmesh.NodeCount;
refind = test.prm.ref.out;
dt = test.prm.dt;
nstep = test.prm.n_step;
proj = ProjectFieldTD(hmesh,Q,M,dmask,...
    mua0*ones(n_node,1),mus0*ones(n_node,1),...
    0,0,refind*ones(n_node,1),dt,nstep,0,0,'diff',1);
figure,semilogy(proj)
%% compare diffusion
dmask_single = false(size(dmask));
dmask_single(1,5) = true;
proj0 = proj(:,dmask_single(:));
n_in = refind;
n_ext = 1.0; % air
c0 = 0.299;  % mm/ps
Tmax = dt * nstep; % ps
pos_Q = test.prm.mesh.QM.Qpos(1,:);
pos_M = test.prm.mesh.QM.Mpos(5,:);
% 
t = 1:dt:dt*nstep;
A = A_factor(n_in./n_ext);
y = SemiInfinite_TR(t,pos_Q,pos_M,mua0,mus0,c0./n_in,A);
figure,semilogy(t,[y'./max(y),proj0./max(proj0)]),ylim([1e-3 1]),
xlabel('time (ps)'),ylabel('mm^{-2}ps^{-1}')
end