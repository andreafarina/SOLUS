function [maskIm, cor] = pointSplineSegs(import)
%-
% 
% Takes an image as input and starts the segmentation routine


% if (strcmpi(imName,'choose') == 1);
%     
%     [caseName, direc,~] = uigetfile('*.jpg', 'Choose your Image');
%     imName = [direc caseName];
% end

 % think of using dicomread(imName) with dicom images...to be checked
DISP = 0;
 
%disp('Segmentation...');
[maskIm, cor] = roispline(uint8(import));

if DISP 
    figure;
    image(maskIm*255); axis image; 
end

return;
end
