 %A = imread('fromDeliverable.png');
 A = imread('INCL_2cm.png');
 info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0 = 0;
 dicomwrite(A, 'US_image_dicom_LARGE.dcm', info,'CreateMode','Create');

 A = dicomread('US_image_dicom_LARGE.dcm');
 imshow(A);res = ginput;
 A = dicomread('US_image_dicom_LARGE.dcm'); info = dicominfo('US_image_dicom_LARGE.dcm');
 info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0 = res(1,1);
 info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxX1 = res(2,1);
info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxY1 = res(3,2);
info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0 = res(2,2);
info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX = 1.5 / (res(6,2) - res(5,2));
info.SequenceOfUltrasoundRegions.Item_1.PhysicalUnitsYDirection = 3;
info.SequenceOfUltrasoundRegions.Item_1.PhysicalUnitsXDirection = 3;
info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaY = 1.5 / (res(6,2) - res(5,2));
dicomwrite(A, 'US_image_dicom_LARGE.dcm', info, 'CreateMode','copy');


%% create
filename = 'TEST_VICRE7VICTRE.mat';
load(filename);
%imshow(out_us);
% out_us = zeros(500,500);
% lx = 20*10;
% lz = 8*10;
% ld = 100;
% out_us(ld:ld+lz,250-lx/2:250+lx/2)=1;
out_us = uint8(repmat(out_us*255,[1 1 3]));
info = dicominfo('../DICOMimages/Mdicom1.dcm');

info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinX0 = 1;
info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxX1 = size(out_us,2);
info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMaxY1 = size(out_us,1);
info.SequenceOfUltrasoundRegions.Item_1.RegionLocationMinY0 = 1;
info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX = 0.01;
info.SequenceOfUltrasoundRegions.Item_1.PhysicalUnitsYDirection = 3;
info.SequenceOfUltrasoundRegions.Item_1.PhysicalUnitsXDirection = 3;
info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaY = 0.01;
dicomwrite(out_us,[filename(1:end-4),'.dcm'], info, 'CreateMode','Copy');
