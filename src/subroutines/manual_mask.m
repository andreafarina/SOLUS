function [mask] = manual_mask(spe, ic, grid, conc, coeff)
Data = reshape(conc,grid.dim); zmin = grid.z1; zmax = grid.z2; 
dz = grid.dz; step = 1; 
x = grid.x;
y = grid.y;
z = (zmin:dz:(zmax-grid.dz)) + grid.dz/2;
N = numel(z);
l = 1;
cmin = min(Data(:));
cmax = max(Data(:));
Data_dummy = Data(:,:,z>10);
if (isempty(Data_dummy) == 0)
    cmin = min(Data_dummy(:));
    cmax = max(Data_dummy(:));
end
clim = double([cmin * (1-eps) cmax * (1+eps)]);
for i=1:step:N
    figure(i);
    imagesc(x,y,Data(:,:,i)',clim)
    set(gca,'ydir','normal'),
    ylabel('y(mm)'),
    xlabel('x(mm)'),
    axis('image')
    cm = colorbar;
    title(['Recon ' coeff ' at z = ' num2str(z(i)) ' mm']);
    l = l+1;
    roi = drawcircle('Color','k','FaceAlpha',0.4); %press ESC if no ROI has to be selected
    if isempty(roi.Center) == 1
        mask(:,:,i) = zeros(grid.dim(1),grid.dim(2));
    else
        a = createMask(roi);
        mask(:,:,i) = a';
    end
end
end

% 
% for ic = 1:REC.spe.nCromo
%     Conc=reshape(REC.opt.bConc(:,ic),REC.grid.dim);
%     fh=figure(800+ic);fh.NumberTitle = 'off';fh.Name = ['Recon ' REC.spe.cromo_label{ic} ' Map'];
%     ShowRecResults(REC.grid,Conc,...
%         REC.grid.z1,REC.grid.z2,REC.grid.dz,1,'auto');%,0.,0.64);
%     suptitle(['Recon ' REC.spe.cromo_label{ic}]);
%     savefig(fh,['./figures/' fh.Name '.fig'])
%     savefig(fh,[fh.Name '.fig'])
% end











