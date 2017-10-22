%-------------------------------------------------------------------------%
% Add a generic inhomogeneity
%
% var       -- string   -- the field you wanna mofify
% back      -- [1x1]    -- background value of 'var'
% INTENSITY -- [1x1]    -- maximum value of the inhomogeneous 'var'
% 
% No profiles are implemented. We get a binary tridimensional mask.
%-------------------------------------------------------------------------%





function [DOT, hete] = priormask3D(DOT, hete)

var  = hete.type;
back = eval(['DOT.opt.',lower(hete.type),'B']);
intensity = hete.val;
%distrib   = hete.distrib;


if isfield(DOT,'grid'),
    nx = DOT.grid.Nx;
    ny = DOT.grid.Ny;
    nz = DOT.grid.Nz;
    add = zeros(DOT.grid.N,1);
    Nx = DOT.grid.Nx;
    Ny = DOT.grid.Ny;
    Nz = DOT.grid.Nz;
    
%-- Regular mesh --%   
else
    M = [DOT.mesh.pos(:,1)-c(1) DOT.mesh.pos(:,2)-c(2) DOT.mesh.pos(:,3)-c(3)];
    add = zeros(DOT.mesh.N,1);
    Nx = DOT.mesh.Nx;
    Ny = DOT.mesh.Ny;
    Nz = DOT.mesh.Nz;
end

mask = load(hete.path);
mask = mask.ellipsoid;
mask = mask/255;
%mask = permute(mask, [1, 3, 2]);
resized = imresizen(mask,[ (nx/size(mask,1)) (ny/size(mask,2))  (nz/size(mask,3))], 'nearest');
 

idx = find(resized > 0);
% resized is 1 where the inclusion is
param = getfield(DOT.opt, var); 


% add(:) = resized(:) * (intensity - back) + back; % what should this be in
% for neat boundaries?
% 
param(idx) = param(idx) + (intensity - back);



hete    = setfield(hete, ['d',hete.type], add);
DOT.opt = setfield(DOT.opt, var, param);

% Uncomment to plot the shape of the inclusion
[x, y, z] = ind2sub(size(mask), find(mask));
figure; xlabel('plot3D mask'); plot3(x, y, z, 'k.');

[x, y, z] = ind2sub(size(resized), find(resized));
figure; xlabel('plot3D resized'); plot3(x, y, z, 'k.'),
set(gca,'zdir','reverse'),axis equal,
xlim([DOT.grid.x1 DOT.grid.x2]),...
    ylim([DOT.grid.y1 DOT.grid.y2]),...
    zlim([DOT.grid.z1 DOT.grid.z2])
return;
end



