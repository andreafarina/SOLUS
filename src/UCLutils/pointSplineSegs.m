function [maskIm, delta] = pointSplineSegs(imName, CROP)


if (strcmpi(imName,'choose') == 1);
    
    [caseName, direc,~] = uigetfile('*.jpg', 'Choose your Image');
    imName = [direc caseName];
end

import = imread(imName); % think of using dicomread(imName) with dicom images...to be checked

if (CROP == 1)
    disp('Cropping...')
    import = autocropper(import);

end

disp('Segmentation...');
maskIm = roispline(import);


% if (exist('str','var') == 1)
%     if (strcmp(str, 'save') == 1)
%         maskName = [direc,'B_splined/', 'Mask_',caseName];
%         imwrite(maskIm,maskName);
%     end
% end

figure;
image(maskIm*255); axis image; 

%% set delta
delta = 0.1; %mm/px % TO BE FIXED WITH DICOM

return;
end
