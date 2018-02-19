%% SEGMENT AND RETRIEVE THE 3D SHAPE FROM A 2D IMAGE  
clear;
close all;

CROP = 1;   % enables autocropping of the 2D image to get a better approximation 
            % of the position of the inclusion with respect to the US probe. TO BE CHANGED WITH DICOM
SAVE_3D = 0;   % saves 3D image
SAVE_2D = 0;     % save segmented image

loadname = 'benign_3.jpg';
savename3D = '../../3DMasks/myshape.mat';
savename2D = [loadname, '_SEGMENTED'];


%% Segmentation 

[segmented, delta] = pointSplineSegs(loadname, CROP); % get segmented shape % delta is set to be 0.1 mm/px (to be fixed when dicom protocol is present)

%% Extrusion

mask3D = retrieve_ellipsoid(segmented);

%% show shape
figure;
for i = 1:1:25
subplot(5,5,i);imagesc(mask3D(:,:,round((i/25)*size(mask3D,3)))');axis image;
end

figure;
[x, y, z] = ind2sub(size(mask3D), find(mask3D));
plot3(x, y, z, 'k.'); axis image;

%% SAVE
if SAVE_3D == 1
    save(savename3D, 'mask3D', 'delta');
end

if SAVE_2D == 1
    save(savename2D, 'segmented', 'delta');
end



