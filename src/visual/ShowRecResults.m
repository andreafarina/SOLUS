function clim = ShowRecResults(grid,Data,zmin,zmax,dz,step,scale,cmin,cmax)
% scale:    'adaptive': each slice with his own colorbar, subsequent
%            arguments will be ignored
%           'auto':     fixed scale between min(Data) and max(Data)
%           (default if no argument is provided)
%           'user':     cmin, cmax defined by the user
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
if ((nargin < 7)||strcmpi(scale,'auto'))
  cmin = min(Data(:));
  cmax = max(Data(:));
%  cmin = prctile(Data(:),3);
%  cmax = prctile(Data(:),97);
<<<<<<< HEAD
end
Data_dummy = Data(:,:,z>10);
if ((nargin < 7)||strcmpi(scale,'auto'))
  cmin = min(Data_dummy(:));
  cmax = max(Data_dummy(:));
end
clim = double([cmin * (1-eps) cmax * (1+eps)]);
=======
end
Data_dummy = Data(:,:,z>10);
if ((nargin < 7)||strcmpi(scale,'auto'))
  cmin = min(Data_dummy(:));
  cmax = max(Data_dummy(:));
end
clim = [cmin-eps cmax+eps];
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
H=subplot1(ceil(sqrt(Nsub/1.5)),ceil(1.5*sqrt(Nsub/1.5)),'Gap',[0.001 0.001],'Min',[0.05 0.07],'Max',[1.01 0.98],'XTickL', 'All', 'YTickL', 'Margin','YScale','linear','FontS',18);
for i=1:step:N
           % Data(:,:,i) = Data(:,:,i)';
%            subplot(ceil(sqrt(Nsub/1.5)),ceil(1.5*sqrt(Nsub/1.5)),l),
                    subplot1(i);
                    if strcmpi(scale,'adaptive')
                     imagesc(x,y,Data(:,:,i)')
                     set(gca,'ydir','normal'),
                     if rem(i,ceil(1.5*sqrt(Nsub/1.5)))==1
                         ylabel('y(mm)'),
                     end
                     if rem(i/ceil(1.5*sqrt(Nsub/1.5)),ceil(sqrt(Nsub/1.5)))==0
                         xlabel('x(mm)'),
                     end
                     axis('image')
                    else
                     imagesc(x,y,Data(:,:,i)',clim)
                     set(gca,'ydir','normal'),
                     if rem(i,ceil(1.5*sqrt(Nsub/1.5)))==1
                         ylabel('y(mm)'),
                     end
                     if i/ceil(1.5*sqrt(Nsub/1.5))>(ceil(sqrt(Nsub/1.5))-1)
                         xlabel('x(mm)'),
                     end
                     axis('image')
                    end
                     %pcolor(Data(:,:,i)'), shading interp
              %axis image%,axis xy
            
            %end
            %title(['mu_a at z = ' num2str(REC.grid.z(i))]);axis square;colorbar
            title(['z = ' num2str(z(i)) ' mm']);%axis square;
            
            %xlabel('x'),ylabel('y')
           l = l+1;
           cm = colorbar;
           %set(cm, 'YTickLabel', cellstr(num2str(reshape(get(cm, 'YTick'),[],1),'%0.4f')) )
end
ExSub = numel(H)-N;
for ie = 1:(N-ExSub-1)
    set(H(ie),'XTickLabel',[]);
end
for ie = (N-ExSub):(N-1)
    axes(H(ie))
    xlabel('x(mm)')
end
delete(H(i+1:end))
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