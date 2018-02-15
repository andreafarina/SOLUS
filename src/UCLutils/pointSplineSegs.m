function maskIm = pointSplineSegs(str)

close all
varargin

[caseName, direc,~] = uigetfile('3D_simulations/Segmentation/Milan_Data/*.jpg', 'Choose your Image');
imName = [direc caseName];
import = imread(imName); % think of using dicomread(imName) with dicom images...to be checked

im = autocropper(import);

maskIm = roispline(im);


if (exist('str','var') == 1)
    if (strcmp(str, 'save') == 1)
        maskName = [direc,'B_splined/', 'Mask_',caseName];
        imwrite(maskIm,maskName);
    end
end

figure;
image(maskIm*255); axis image;  
return;
end
