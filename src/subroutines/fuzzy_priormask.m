
function mask = fuzzy_priormask(mu,delta,grid)

Nx = grid.Nx;
Ny = grid.Ny;
Nz = grid.Nz;

DOT_GRID = 1;FACT = [1,1,1];


final_dims = [grid.x1,grid.y1,grid.z1;...
    grid.x2,grid.y2,grid.z2];
mask_oversampled =  segment2grid(mu, delta, final_dims, 'fuzzy');

 
if DOT_GRID == 1 %if manual mapping of prior of DOT grid 
   
    mask = (imresizen(single(mask_oversampled),...
            [Nx,Ny,Nz]./size(mask_oversampled),'linear'));
    return;

else
    if exist('type', 'var') == 1
         if strcmpi(type,'fit4param') == 1 % when considering TOAST is better not to resize it, proportions are already handled  

        mask = (imresizen(single(mask_oversampled),...
            FACT,'linear'));
        return;

         end
     else

        mask = (imresizen(single(mask_oversampled),...
            [Nx,Ny,Nz]./size(mask_oversampled),'linear'));
        return;
    end

% end
end



