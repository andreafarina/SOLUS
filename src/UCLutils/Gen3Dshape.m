%% SEGMENT AND RETRIEVE THE 3D SHAPE FROM A 2D IMAGE  
clear;
close all;

CROP = 1;   % enables autocropping of the 2D image to get a better approximation 
            % of the position of the inclusion with respect to the US
            % probe. Preferably set to 1 with dicoms
SAVE_3D = 0;   % saves 3D image
SAVE_2D = 0;     % save segmented image
SAVE_JPG = 0; % save images for display
FORMAT = 'DICOM'; % or ELSE for other formats (jpg, png)


loadname = 'DICOMimages/Mdicom3.dcm';
savename3D = [loadname(1: end - 4), '_3D.mat'];
savename2D = [loadname(1: end - 4), '_SEGMENTED.mat'];
savenameJPG = [loadname(1: end - 4), '.jpg'];
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
        im = autocropper(im);
    end         
    % delta is set to be 0.1 mm/px just as an example 
    %delta = 0.1; %mm/px
    delta = 0.07; %mm/px FIELDII
end

[segmented, cor] = pointSplineSegs(im);
close all
sgm = snake_fitting(im, cor);
%% Extrusion
disp('Extrusion...');
mask3D = retrieve_ellipsoid(segmented);


%% SAVE
mask3D = logical(mask3D(:,:,1:min(500, size(mask3D,3))));
if SAVE_3D == 1
    save(savename3D, 'mask3D', 'delta');
end

if SAVE_2D == 1
    save(savename2D, 'segmented', 'delta');
end

%% show and save shape
if SAVE_JPG == 1;
    
    h_2D = figure();
    imagesc(segmented(:,:)), axis image;
    saveas(h_2D, ['MASK',savenameJPG]);
    
    h_2Ds = figure();
    
    if size(im,3) ~=3
        todisp(:,:,:) = uint8(cat(3, im,im,im));
    else
        todisp(:,:,:) = im(:,:,:);
    end
    
    todisp(:,:,2) = todisp(:,:,2) + 100 * uint8(segmented(:,:));
    image(todisp), axis image
    saveas(h_2Ds, ['MASKoverDCM_',savenameJPG]);
    
    h_dcm = figure();
    imagesc(im), axis image, colormap(gray(127))
    saveas(h_dcm, ['readDCM_',savenameJPG]);
    
    
    
    h_dt2D = figure();
    for i = 1:1:25
    subplot(5,5,i);image(1:round(size(mask3D,1) * delta), 1:round(size(mask3D,2) * delta), 255* mask3D(:,:,round((i/25)*size(mask3D,3)))');axis image;
    %xlabel('x(mm)'),ylabel('y(mm)')%title(sprintf('%d of %d', round((i/25)*size(mask3D,3)), size(mask3D,3))); 
    end
    saveas(h_dt2D, ['3Dsliced', savenameJPG]);
    
    h_dt3D = figure();
    [x, y, z] = ind2sub(size(mask3D), find(mask3D));
    plot3(x*delta, y*delta, z*delta, 'k.'); axis image;xlabel('x(mm)'),ylabel('y(mm)'),zlabel('z(mm)')
    saveas(h_dt3D, ['3Dwhole', savenameJPG]);

end


