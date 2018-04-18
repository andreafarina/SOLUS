%% SEGMENT AND RETRIEVE THE 3D SHAPE FROM A 2D IMAGE  
clear;
close all;

CROP = 1;   % enables autocropping of the 2D image to get a better approximation 
            % of the position of the inclusion with respect to the US
            % probe. Preferably set to 1 with dicoms
SAVE_3D = 1;   % saves 3D image
SAVE_2D = 1;     % save segmented image

FORMAT = 'DICOM'; % or ELSE for other formats (jpg, png)


loadname = 'Mdicom4.dcm';
savename3D = 'Mdicom1_3D.mat';
savename2D = [loadname, '_SEGMENTED'];
%% segmentation
if strcmpi(FORMAT, 'DICOM') == 1 
    
    im = dicomread(loadname);
    info = dicominfo(loadname);
    
    if CROP == 1
        xmin = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0;
        xmax = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxX1;
        ymin = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0;
        ymax = info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxY1;
        im = im(ymin:ymax, xmin:xmax, :);
    end
    
    % find delta
    delta_x = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX;
    delta_y = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaY;
    unx = info.SequenceOfUltrasoundRegions.Item_1.PhysicalUnitsXDirection; % gives in which units physical delta x is expressed 
    uny = info.SequenceOfUltrasoundRegions.Item_1.PhysicalUnitsYDirection;
    
    
    if unx == 3
        delta_x = delta_x * 10; % transform deltax from cm/px to mm/px
    else
        disp('Error: see dicom dictionary to know what units are assigned to delta_x')
    end
    
    if uny == 3
        delta_y = delta_y * 10;
    else
        disp('Error: see dicom dictionary to know what units are assigned to delta_x')
    end
    
    if delta_x ~= delta_y
       
        disp('DICOM image has problems with calibrations: pixels are not squares')
        delta = 0.5 * (delta_x + delta_y);
    else
        delta = delta_x;
    
    end
    
else
    im = imread(loadname);
    
    if CROP == 1
        disp('Cropping...')
        import = autocropper(import);
    end         
    % delta is set to be 0.1 mm/px just as an example 
    delta = 0.1; %mm/px
end

segmented = pointSplineSegs(im);

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



