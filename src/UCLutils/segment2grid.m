function out_grid =  segment2grid(mask, delta, final_dims)
%-
%-
%-
% out_grid = segment2grid(MASK, DELTA, FINAL_DIMS)
%
% Takes a slab (MASK 3d volume), a factor of conversion DELTA (mm/pixel) and 
% the desired dimensions (FINAL_DIMS) and gives back a 3D matrix with the same real 
% dimensions given in input and same factor of conversion. 
% The new matrix will have the input matrix inserted in its centre 
% along x and y, and at index = 1 along z.
% final_dims will be a matrix 2x3, with first row representing the lower
% extremes x, y and z and second row representing the upper ones.
% If out dimensions are smaller than the mask then this will be cut.
%
% MASK: 3d array
% DELTA: number mm/pixel
% FINAL_DIMS: [x1,y1,z1; x2,y2,z2] (mm)
% OUT GRID:  3d array

DISPLAY = 0;


% rewriting the dimensions of the out grid
infx = final_dims(1,1);
infy = final_dims(1,2);
infz = final_dims(1,3);
supx = final_dims(2,1);
supy = final_dims(2,2);
supz = final_dims(2,3);

% No of voxels of the in mask
[smaskx, smasky, smaskz] = size(mask);

% N0 of voxels of the out grid
outdimx = numel( infx: delta : supx);
outdimy = numel( infy: delta : supy);
outdimz = numel( infz: delta : supz);

out_grid = false(outdimx, outdimy, outdimz);

% Defining middle indexes
mx = round(outdimx/2);
my = round(outdimy/2);

midmaskx = ceil(smaskx/2);
midmasky = ceil(smasky/2);

smaskx = min( smaskx, outdimx);
smasky = min( smasky, outdimy);
smaskz = min( smaskz, outdimz);

% insert slab in out grid
out_grid(mx - ceil(smaskx/2) + (1:smaskx),...
         my - ceil(smasky/2) + (1:smasky),... 
                                   1:smaskz) = ...
                mask(midmaskx - ceil(smaskx/2) + (1:smaskx),...
                        midmasky - ceil(smasky/2) + (1:smasky),...
                        1:smaskz);

if (DISPLAY == 1)

    for i = 1:1:size(out_grid,3)
        
        figure(1); imagesc(1:size(out_grid,1) * delta , 1:size(out_grid,2) * delta, out_grid(:,:,i)); title(sprintf('%d', i));axis image, pause
    end
    
end

 
end