function [maskIm] = pointSplineSegs(import)
%-
% 
% Takes an image as input and starts the segmentation routine


% if (strcmpi(imName,'choose') == 1);
%     
%     [caseName, direc,~] = uigetfile('*.jpg', 'Choose your Image');
%     imName = [direc caseName];
% end

 % think of using dicomread(imName) with dicom images...to be checked

disp('Segmentation...');
maskIm = roispline(uint8(import));


figure;
image(maskIm*255); axis image; 


return;
end
