function mask3D = dicom_to_DT(im,flag_snake)
%  dicom_to_DT(im,flag_snake)
% given an image starts the manual segmentation and extraxtion procedure.
% if flag_snake is set to 1, then the snake fitting is also applied to the
% segmentation
% returns a 3D binary matrix as a result of the segmentation and DT
% extrapolation


[sgm, cor] = pointSplineSegs(im);
if nargin > 1
    if flag_snake ==1
        [sgm, ~ ] = snake_fitting(im, cor);
    end
end
mask3D = logical(retrieve_ellipsoid(sgm));

end