
function mask = priormask3D(path,grid, type)

Nx = grid.Nx;
Ny = grid.Ny;
Nz = grid.Nz;

DOT_GRID = 1;FACT = [1,1,1];



%% here!
smask = load(path);
fn = fieldnames(smask);
mask = smask.(fn{1});
delta = smask.(fn{2});


%% swap fields... in case
if ~isvector(delta)
    dd = mask;
    mask = delta;
    delta = dd;
end

final_dims = [grid.x1,grid.y1,grid.z1;...
    grid.x2,grid.y2,grid.z2];
mask_oversampled =  segment2grid(mask, delta, final_dims);

 
if DOT_GRID == 1 %if manual mapping of prior of DOT grid 
   
    mask = logical(imresizen(single(mask_oversampled),...
            [Nx,Ny,Nz]./size(mask_oversampled),'nearest'));
    return;

else
    if exist('type', 'var') == 1
         if strcmpi(type,'fit4param') == 1 % when considering TOAST is better not to resize it, proportions are already handled  

        mask = logical(imresizen(single(mask_oversampled),...
            FACT,'nearest'));
        return;

         end
     else

        mask = logical(imresizen(single(mask_oversampled),...
            [Nx,Ny,Nz]./size(mask_oversampled),'nearest'));
        return;
    end

% end
end



