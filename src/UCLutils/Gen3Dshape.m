%% SEGMENT AND RETRIEVE THE 3D SHAPE FROM A 2D IMAGE  
clear;
close all;

<<<<<<< HEAD
CROP = 1;   % enables autocropping of the 2D image to get a better approximation 
=======
CROP = 0;   % enables autocropping of the 2D image to get a better approximation 
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
            % of the position of the inclusion with respect to the US
            % probe. Preferably set to 1 with dicoms
SAVE_3D = 1;   % saves 3D image
SAVE_2D = 1;     % save segmented image
SAVE_JPG = 1; % save images for display
<<<<<<< HEAD
FORMAT = 'DICOM'; % or ELSE for other formats (jpg, png)

[filename, pathname] = uigetfile('./*.dcm');
loadname = [pathname,filename];pathname = './';
%loadname = 'DICOMkwaved-BALL=xpos_0.03__rad_0.015__deltaScat_-50.dcm';%'Bmode_FieldII.jpg';
savename3D = ['DT',filename(1: end - 4), '_3D.mat'];
savename2D = [filename(1: end - 4), '_SEGMENTED.mat'];
savenameJPG = [filename(1: end - 4), '.jpg'];
=======
FORMAT = 'ELSE'; % or ELSE for other formats (jpg, png)


loadname = 'Bmode_FieldII.jpg';
savename3D = [loadname(1: end - 4), '_3D.mat'];
savename2D = [loadname(1: end - 4), '_SEGMENTED.mat'];
savenameJPG = [loadname(1: end - 4), '.jpg'];
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
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
<<<<<<< HEAD
        disp('Error: see dicom dictionary to know what units are assigned to delta_y')
=======
        disp('Error: see dicom dictionary to know what units are assigned to delta_x')
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
    end
    
    if delta_x ~= delta_y
       
<<<<<<< HEAD
        disp('DICOM image has problems with calibrations: pixels are not squares. Setting delta_y as standard')
        im = uint8(imresizen(double(im), [delta_x/delta_y,1,1] ));
        delta_x = delta_y;
        %delta_x = 0.5 * (delta_x + delta_y);
        delta = delta_x;
        
=======
        disp('DICOM image has problems with calibrations: pixels are not squares')
        delta = 0.5 * (delta_x + delta_y);
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
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

segmented = pointSplineSegs(im);

%% Extrusion
disp('Extrusion...');
mask3D = retrieve_ellipsoid(segmented);


%% SAVE
<<<<<<< HEAD
siz3 = floor(min([ 70/delta , size(mask3D,3) ]));
mask3D = logical(mask3D(:,:,1:siz3));
if SAVE_3D == 1
    save([pathname,savename3D], 'mask3D', 'delta');
end

if SAVE_2D == 1
    save([pathname,savename2D], 'segmented', 'delta');
end

%% show and save shape
if SAVE_JPG == 1
    
    h_2D = figure();
    imagesc(segmented(:,:)), axis image;
    saveas(h_2D, [pathname,'MASK',savenameJPG]);
=======
mask3D = logical(mask3D(:,:,1:500));
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
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
    
    h_2Ds = figure();
    
    if size(im,3) ~=3
        todisp(:,:,:) = uint8(cat(3, im,im,im));
    else
        todisp(:,:,:) = im(:,:,:);
    end
    
    todisp(:,:,2) = todisp(:,:,2) + 100 * uint8(segmented(:,:));
    image(todisp), axis image
<<<<<<< HEAD
    saveas(h_2Ds, [pathname,'MASKoverDCM_',savenameJPG]);
    
    h_dcm = figure();
    imagesc(im), axis image, colormap(gray(127))
    saveas(h_dcm, [pathname,'readDCM_',savenameJPG]);
=======
    saveas(h_2Ds, ['MASKoverDCM_',savenameJPG]);
    
    h_dcm = figure();
    imagesc(im), axis image, colormap(gray(127))
    saveas(h_dcm, ['readDCM_',savenameJPG]);
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
    
    
    
    h_dt2D = figure();
    for i = 1:1:25
    subplot(5,5,i);image(1:round(size(mask3D,1) * delta), 1:round(size(mask3D,2) * delta), 255* mask3D(:,:,round((i/25)*size(mask3D,3)))');axis image;
    %xlabel('x(mm)'),ylabel('y(mm)')%title(sprintf('%d of %d', round((i/25)*size(mask3D,3)), size(mask3D,3))); 
    end
<<<<<<< HEAD
    saveas(h_dt2D, [pathname,'3Dsliced', savenameJPG]);
=======
    saveas(h_dt2D, ['3Dsliced', savenameJPG]);
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98
    
    h_dt3D = figure();
    [x, y, z] = ind2sub(size(mask3D), find(mask3D));
    plot3(x*delta, y*delta, z*delta, 'k.'); axis image;xlabel('x(mm)'),ylabel('y(mm)'),zlabel('z(mm)')
<<<<<<< HEAD
    saveas(h_dt3D, [pathname,'3Dwhole', savenameJPG]);
=======
    saveas(h_dt3D, ['3Dwhole', savenameJPG]);
>>>>>>> 601f76d0bf4879d27e7741946f367a2a80e5bb98

end


