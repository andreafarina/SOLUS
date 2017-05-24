function clim = ShowRecResults(grid,Data,zmin,zmax,dz,step,cmin,cmax)
%% preset the subplot grid
%BB = hMesh.BoundingBox();
%[bmin bmax] = toastMeshBB(hMesh);
bmin = [grid.x1;grid.y1;grid.z1]; 
bmax = [grid.x2;grid.y2;grid.z2];
%bdim  = toastGetBasisSize(hBasis);
bdim = grid.dim;
%solmask = toastSolutionMask(hBasis);

%solmask = hBasis.Map('B->S',ones(bdim')); % *** this might need looking at ***
%solmask = find(hBasis.GridElref>0);
%solmask = ones(size(Data(:)));
dd = (bmax(:)-bmin(:))./bdim(:);%-1);    % the basis can be smaller than the mesh! attention!

x = grid.x;%bmin(1):dd(1):bmax(1);
y = grid.y;%bmin(2):dd(2):bmax(2);
z = (zmin:dz:(zmax-grid.dz)) + grid.dz/2;
N = numel(z);
Nsub = ceil(N./step);
%%
%figure;clf;
l = 1;
if nargin < 9
  cmin = min(Data(:));
  cmax = max(Data(:));
end
clim = [cmin-eps cmax+eps];
for i=1:step:N
            subplot(ceil(sqrt(Nsub/1.5)),ceil(1.5*sqrt(Nsub/1.5)),l),
                     imagesc(Data(:,:,i)',clim),
              axis image%,axis xy
            
            %end
            %title(['mu_a at z = ' num2str(REC.grid.z(i))]);axis square;colorbar
            %title(['z = ' num2str(z(i))]);%axis square;
            
            %xlabel('x'),ylabel('y')
           l = l+1;
end
% subplot(ceil(sqrt(Nsub/1.5)),ceil(1.5*sqrt(Nsub/1.5)),N+2);
% imagesc((mean(clim,2))*ones(1,1),clim),colorbar,
% axis image,axis xy;
           




% if ~isfield(REC,'Datafile')
%     subplot(2,3,1),imagesc(REC.grid.x,REC.grid.y,REC.ref.Mua(:,:,offset+zsec)'),
%     title(['mu_a at z = ' num2str(offset+zsec)]);axis square;colorbar
%     xlabel('x (mm)'),ylabel('y (mm)')
%     subplot(2,3,2),imagesc(REC.grid.x,REC.grid.y,REC.ref.Musp(:,:,offset+zsec)'),
%     title(['mu_{sp} at z = ' num2str(offset+zsec)]);axis square;colorbar
%     xlabel('x (mm)'),ylabel('y (mm)')
%     subplot(2,3,3),imagesc(REC.grid.x,REC.grid.y,REC.ref.Conc(:,:,offset+zsec)'),
%     title(['[c] at z = ' num2str(offset+zsec)]);axis square;colorbar
%     xlabel('x (mm)'),ylabel('y (mm)')
% end
% subplot(2,3,4),imagesc(REC.grid.x,REC.grid.y,REC.opt.bmua(:,:,zsec)'),
% title(['Rec mu_a at z = ' num2str(zsec)]);axis square;colorbar
% xlabel('x (mm)'),ylabel('y (mm)')
% subplot(2,3,5),imagesc(REC.grid.x,REC.grid.y,REC.opt.bmusp(:,:,zsec)'),
% title(['Rec mu_{sp} z = ' num2str(zsec)]);axis square;colorbar
% xlabel('x (mm)'),ylabel('y (mm)')
% subplot(2,3,6),imagesc(REC.grid.x,REC.grid.y,REC.opt.bconc(:,:,zsec)'),
% title(['Rec [c] at z = ' num2str(zsec)]);axis square;colorbar
% xlabel('x (mm)'),ylabel('y (mm)')
% drawnow;