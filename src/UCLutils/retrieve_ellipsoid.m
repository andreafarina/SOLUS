function ellipsoid = retrieve_ellipsoid(im_mask)
% find 3D shape from 2D mask using an extrusion approach based on the
% distance transform.

im_mask = logical(im_mask ~= 0);
sm_fact = 1; % set to 1 for 'circular' extrusion
[D, ~] = bwdist(~im_mask, 'euclidean');
clear im_mask

MAX_third = max(D(:));
MAX_D = MAX_third; 
MAX_x = size(D, 1);
MAX_y = size(D, 2);

%MAX_third = 0.5 * ceil(max(MAX_x, MAX_y)/2);
%MAX_third = ceil(MAX_third);

%% Smoothing factor  
Area = sum(D(:) ~= 0);
R = MAX_D;
sm_fact = sqrt(Area/ pi)/R;

MAX_third = ceil(sm_fact * MAX_D);
half_mask = false(MAX_x, MAX_y, ceil(MAX_third));
%% for procedure 
tic;
idxk = (find(D >= 1));
[idxi, idxj] = ( ind2sub(size(D),find( D >= 1)));  
%idxi = (idxi);
%idxj = (idxj);
idxk = (ceil( sm_fact * sqrt( 2 * D(idxk) * MAX_D -  D(idxk).^2)) );
disp('Generating 3D Mask')

 for j = 1: numel(idxk)
    half_mask(idxi(j), idxj(j), 1:idxk(j) ) = true;
 end
 
 
new_mask = cat(3, half_mask(:,:, end: -1:2), half_mask);
ellipsoid = flip(new_mask, 2);
ellipsoid = flip(new_mask, 1);
ellipsoid = permute(new_mask,[2,3,1]); % x oriented towards line of optodes, x across it, z is depth

disp('3D mask has been generated')
toc

return;


%% Uncomment to display thwe results of the procedure
% %new_mask = 255*new_mask;
% figure;
% for i = 1:1:25
% subplot(5,5,i);imagesc(ellipsoid(:,:,round((i/25)*size(ellipsoid,3)))); axis image
% end
% 
% 
% [x, y, z] = ind2sub(size(ellipsoid), find(ellipsoid));
% figure;plot3(x, y, z, 'k.'); axis image



end

% if (saveme == 1)
%    
%    CaseName = strrep(CaseName,'.jpg', '');
%    direc = strrep(direc, 'B_splined', 'Masks_3D' );
%    maskName = [direc, 'Mask3D_', CaseName];
%    save(maskName, 'ellipsoid');
%     
% end


        


