function [mu,cov_matrix,area] = Gaussian3D(x,y,z,dx,dy,dz,data)
% =========================================================================
%                            Gaussian 3D
% Find center-of-mass (mu), covariance matrix(cov_matrix), and area of 
% a 3-fold variate gaussian distribution.
% size(data) = [numel(x),numel(y),numel(z)];
% =========================================================================
%% check dimensions consistency
size_space = [numel(x),numel(y),numel(z)];
if sum(size_space-size(data))~=0
    disp('Attention! Check the dimensions of x,y,z and data!');
    return
end
%% Prepare space variables
[X,Y,Z] = ndgrid(x,y,z);
space = [X(:), Y(:), Z(:)];
dV = dx*dy*dz;
%% Center of mass
data = data(:)';
mu = (data * space )./ sum(data);
%% Covariance matrix
data3dd = repmat(data(:)',[3 1]);
off_space = bsxfun(@minus,space,mu);
cov_matrix = ((data3dd.*off_space')*off_space)./sum(data(:));
%% Area
area = sum(data)*dV;
