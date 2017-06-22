function [y,H] = PreProcessECBO(data,imeas,twin,counts)
% process RhoZero forward data to simulate the experimental conditions
% of nwin detectors
% data:        nt x nm    td data matrix
% imeas:       index of the measruement to take as reference for the
%              scaling factor
% twin:     [nwin,2]    twin matrx
% counts:   [1,nwin]    vector containing the area to impose in each
% window
% y:           scaled data
% H:           stepwise function used for the correction

% Andrea Farina    22/06/2017
%

% create a stepwise function
[nt,nm] = size(data);
nwin = size(twin,1);
H = zeros(1,nt);
area = zeros(nwin,nm);
for i = 1:nwin
    H(twin(i,1):twin(i,2)) = counts(i) / ...
        sum(data(twin(i,1):twin(i,2),imeas));
end
%y = H' .* data;    % for Matlab 2015 and earlier versions
y = bsxfun(@times,H',data);



