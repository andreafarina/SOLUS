%% SHOW RESULTS

close all
%clear all
%% Data file
%data_file='TRallNormal';
%data_file='TRall';
%data_file='TRNull';
%data_file='CW';
%% may 2014
%data_file='TRNull_GMRES_alpha1_threshold40mm'
%data_file='TRAll_GMRES_alpha0540mm'
%data_file='CW_GMRES_alpha140mm'
%data_file='CW_GMRES_noGauss40mm'
diary('Tomo_Comparison.txt')
%% FINAL DATA
%data_file='TRNull_GMRES_alpha1_threshold40mm'; disp('Time resolved'); % TR
data_file='CW_GMRES_alpha140mm';   disp('Continuous Wave');        % CW

SAVE=0;
load ([data_file '.mat']);
%% Geometry
nx=length(x);ny=length(y);nz=length(z);%nt=length(T);
Dx=x(2)-x(1);Dy=y(2)-y(1);Dz=z(2)-z(1);%Dt=T(2)-T(1);
%% Color scale
%clim=[-22e-4 5e-4];
%% inclusion position
r_inc=[20 20 40];
[X_in Y_in Z_in]=meshgrid(x,y,z);

%% mask for the ROI
V1=AddInclusion(X_in,Y_in,Z_in,r_inc,15);
% nzv = length(find(V1==1)); % number of points in V
nzv = (find(V1==1)); % number of points in V
 
% different version of mask :
nzm = find(DMUA3d > 0);
V = zeros(size(DMUA3d));
V(nzm) = 1;

% dmua_rec3d=dmua_rec3d.*V;

%% Error norms
RMSE=norm(dmua_rec3d(:)-DMUA3d(:));
RMSE_V=norm(V(:).*(dmua_rec3d(:)-DMUA3d(:)));
disp(['RMSE: ',num2str(RMSE),' in region: ',num2str(RMSE_V)]);
ABSE=norm(dmua_rec3d(:)-DMUA3d(:),1);
ABSE_V=norm(V(:).*(dmua_rec3d(:)-DMUA3d(:)),1);
disp(['1-norm: ',num2str(ABSE),' in region: ',num2str(ABSE_V)]);


idx=r_inc./[Dx Dy Dz];

%% Centre of Mass
 mxm = max(dmua_rec3d(:));
 mind = find(dmua_rec3d == mxm);
 
 COM_TRUE = zeros(3,1);
 vv = sum(DMUA3d(:));
 COM_TRUE(1) = sum(DMUA3d(:).*X_in(:))/vv;
 COM_TRUE(2) = sum(DMUA3d(:).*Y_in(:))/vv;
 COM_TRUE(3) = sum(DMUA3d(:).*Z_in(:))/vv;
 
 COM_REC = zeros(3,1);
 vr = sum(dmua_rec3d(:));
 COM_REC(1) = sum(dmua_rec3d(:).*X_in(:))/vr;
 COM_REC(2) = sum(dmua_rec3d(:).*Y_in(:))/vr;
 COM_REC(3) = sum(dmua_rec3d(:).*Z_in(:))/vr;
 
COM_err = norm(COM_REC-COM_TRUE);
disp(['Difference in centre of mass [',num2str(COM_REC(1)-COM_TRUE(1)),',',num2str(COM_REC(2)-COM_TRUE(2)),',',num2str(COM_REC(3)-COM_TRUE(3)),']']);
disp(['Error in centre of mass ',num2str(COM_err)]);

COM_REC = zeros(3,1);
 vr = sum(V(:).*dmua_rec3d(:));
 COM_REC(1) = sum(V(:).*dmua_rec3d(:).*X_in(:))/vr;
 COM_REC(2) = sum(V(:).*dmua_rec3d(:).*Y_in(:))/vr;
 COM_REC(3) = sum(V(:).*dmua_rec3d(:).*Z_in(:))/vr;
 
COM_err = norm(COM_REC-COM_TRUE);
disp(['Difference in centre of mass in region [',num2str(COM_REC(1)-COM_TRUE(1)),',',num2str(COM_REC(2)-COM_TRUE(2)),',',num2str(COM_REC(3)-COM_TRUE(3)),']']);
disp(['Error in centre of mass in region ',num2str(COM_err)]);

%% covariance of mass
VAR_TRUE = zeros(3,1);
VAR_TRUE(1) = sum(DMUA3d(:).*((X_in(:)-COM_TRUE(1)).^2 ))/vv;
VAR_TRUE(2) = sum(DMUA3d(:).*((Y_in(:)-COM_TRUE(2)).^2 ))/vv;
VAR_TRUE(3) = sum(DMUA3d(:).*((Z_in(:)-COM_TRUE(3)).^2 ))/vv;

VAR_REC = zeros(3,1);
vr = sum(dmua_rec3d(:));
VAR_REC(1) = sum(dmua_rec3d(:).*((X_in(:)-COM_REC(1)).^2 ))/vr;
VAR_REC(2) = sum(dmua_rec3d(:).*((Y_in(:)-COM_REC(2)).^2 ))/vr;
VAR_REC(3) = sum(dmua_rec3d(:).*((Z_in(:)-COM_REC(3)).^2 ))/vr;

disp(['Sqrt of covariance of mass [',num2str(sqrt(VAR_REC(1))),',',num2str(sqrt(VAR_REC(2))),',',num2str(sqrt(VAR_REC(3))),']']);

vr = sum(V(:).*dmua_rec3d(:));
VAR_REC(1) = sum(V(:).*dmua_rec3d(:).*((X_in(:)-COM_REC(1)).^2 ))/vr;
VAR_REC(2) = sum(V(:).*dmua_rec3d(:).*((Y_in(:)-COM_REC(2)).^2 ))/vr;
VAR_REC(3) = sum(V(:).*dmua_rec3d(:).*((Z_in(:)-COM_REC(3)).^2 ))/vr;

disp(['Sqrt of covariance of mass in region [',num2str(sqrt(VAR_REC(1))),',',num2str(sqrt(VAR_REC(2))),',',num2str(sqrt(VAR_REC(3))),']']);

%% Sigma and CNR
sigma=std(dmua_rec3d(:));
disp(['Standard deviation ',num2str(sigma)]);

av=max(V(:).*dmua_rec3d(:));
disp(['Max ',num2str(av)]);

disp(['CNR ',num2str(av./sigma)]);

%% Sections Data
figure(1); clf; hold on;
nnn=sqrt(nz/1.5);
nrow=ceil(nnn);
ncol=ceil(1.5*nnn);
for iz = 1:nz
    subplot(nrow,ncol,iz); imagesc(x,y,-DMUA3d(:,:,iz)'),colorbar;
    axis square; shading interp; colormap pink;
    xlabel('x (mm)'); ylabel('y (mm)');
    title(['z = ' num2str(z(iz)) ' (mm)']);
end
hold off;

%% Sections Rec
figure(2); clf; hold on;
nnn=sqrt(nz/1.5);
nrow=ceil(nnn);
ncol=ceil(1.5*nnn);

%% Normalize 
%dmua_rec3d=dmua_rec3d./max(dmua_rec3d(:))*100;
for iz = 1:nz
    subplot(nrow,ncol,iz); imagesc(x,y,-dmua_rec3d(:,:,iz)'),colorbar;
    axis square; shading interp; colormap pink;
    xlabel('x (mm)'); ylabel('y (mm)');
    title(['z = ' num2str(z(iz)) ' (mm)']);
end
hold off;
%% Slice at inclusion
figure(3); clf; hold on;
cm=flipud(pink);
%clim=[-50 10];
subplot(2,2,1); imagesc(z,x,squeeze(dmua_rec3d(:,idx(2),:))); 
shading interp; axis square; colormap(cm); colorbar;%title(['xz image at y = ' num2str(y(idx(2)))]);
ylabel('x (mm)'); xlabel('z (mm)');
subplot(2,2,2); imagesc(y,x,squeeze(dmua_rec3d(:,:,idx(3))));
shading interp; axis square; colormap(cm);colorbar; %title(['xy image at z = ' num2str(z(idx(3)))]);
ylabel('x (mm)'); xlabel('y (mm)');
subplot(2,2,3); plot(z,squeeze(dmua_rec3d(idx(1),idx(2),:)),'linewidth',2); axis square; grid off;%ylim(clim);
xlabel('z (mm)'); %title(['y = ' num2str(y(idx(2))) ' (mm), x = ' num2str(x(idx(1)))]);
ylabel('\Delta\mu_a (mm^{-1})');xlim([0 60]);
subplot(2,2,4); plot(x,squeeze(dmua_rec3d(:,idx(2),idx(3))),'linewidth',2); axis square; grid off;%ylim(clim);
xlabel('x (mm)'); %title(['y = ' num2str(y(idx(2))) ' (mm), z = ' num2str(z(idx(3)))]); 
ylabel('\Delta\mu_a (mm^{-1})');xlim([0 60]);
hold off;
if SAVE
    saveas(gcf,[data_file num2str(r_inc(3)) 'mm_TRall'],'eps2c');
end