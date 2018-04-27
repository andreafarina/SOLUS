function [ellipsoid] = retrieve_ellipsoid(im_mask)
% find 3D shape from 2D mask.
% str upload the name of the mask
% saveme = 1 saves the result


close all;
sm_fact = 1; % set to 1 for 'circular' extrusion

% 
% if (exist('str','var') == 1 && exist('direc','var') == 1)
%     CaseName = str;
%     im_mask = imread( [direc CaseName]);
% end
% if (im_mask == 0)
%     [CaseName, direc, ~] = uigetfile('Documents/Segmentation/Milan_Data/B_splined/*.jpg', 'Choose your mask');
%     im_mask = imread( [direc CaseName]);
%     
% % elseif(length(mask) == 1);
% %     saveme = 0;
% %     [CaseName, direc, ~] = uigetfile('Documents/Segmentation/Milan_Data/B_splined/*.jpg', 'Choose your mask');
% %     im_mask = imread( [direc CaseName]);
% end

%im_mask = im_mask/255;
%im_mask = imresize(im_mask, [size(im_mask,1), size(im_mask,2)], 'method', 'nearest');
neg_mask = ~im_mask;
[D, IDX] = bwdist(neg_mask, 'euclidean');


MAX_third = max(D(:));
MAX_D = MAX_third; 
MAX_x = size(D, 1);
MAX_y = size(D, 2);

%MAX_third = 0.5 * ceil(max(MAX_x, MAX_y)/2);
%MAX_third = ceil(MAX_third);

%% Proportionality of smoothing factor
Area = sum(D(:) ~= 0);
R = MAX_D;
sm_fact = sqrt(Area/ pi)/R;

MAX_third = ceil(2 * sm_fact * MAX_D);
new_mask = zeros(MAX_x, MAX_y, MAX_third);
%% for procedure

tic;

for i = 1:1:MAX_x   
    for j = 1:1:MAX_y
        new_mask(i,j,MAX_third) = im_mask(i,j);
        for k = 1:1: MAX_third - 1
                
                if(k <= sm_fact * sqrt( 2 * D(i,j) * MAX_D -  D(i,j)^2) )
                    new_mask(i, j, MAX_third - k) = 1; 
                    new_mask(i, j, MAX_third + k) = 1; 
                else
                    new_mask(i,j,  MAX_third - k) = 0;
                    new_mask(i,j,  MAX_third + k) = 0;
                end
                
        end
    end
end
 
ellipsoid = permute(new_mask,[1,3,2]);
ellipsoid = permute(ellipsoid,[3,2,1]);
ellipsoid = flip(ellipsoid, 3);
ellipsoid = permute(ellipsoid,[2,1,3]);
ellipsoid = flip(ellipsoid,3);
disp('3D mask has been generated')
toc

return;
end
%new_mask = 255*new_mask;
% figure;
% for i = 1:1:25
% subplot(5,5,i);image(ellipsoid(:,:,round((i/25)*size(ellipsoid,3))));
% end
% 
% figure;
% for i = 1:1:25
% subplot(5,5,i);image(new_mask(:,:,round((i/25)*size(new_mask,3))));
% end
% 
% 
% [x, y, z] = ind2sub(size(new_mask), find(new_mask));
% figure;plot3(x, y, z, 'k.');
% 
% 
% [x, y, z] = ind2sub(size(ellipsoid), find(ellipsoid));
% figure;plot3(x, y, z, 'k.');



% if (saveme == 1)
%    
%    CaseName = strrep(CaseName,'.jpg', '');
%    direc = strrep(direc, 'B_splined', 'Masks_3D' );
%    maskName = [direc, 'Mask3D_', CaseName];
%    save(maskName, 'ellipsoid');
%     
% end


        


