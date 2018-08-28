%-------------------------------------------------------------------------%
% Add a generic inhomogeneity
%
% var       -- string   -- the field you wanna mofify
% back      -- [1x1]    -- background value of 'var'
% INTENSITY -- [1x1]    -- maximum value of the inhomogeneous 'var'
% 
% No profiles are implemented. We get a binary tridimensional mask.
%-------------------------------------------------------------------------%





function [DOT, hete] = prior3D(DOT, hete)
for itype = 1:size(hete.type,2)
var  = hete.type{itype};
back = eval(['DOT.opt.',lower(hete.type{itype}),'B']);
intensity = hete.val((1:DOT.radiometry.nL)+(itype-1)*DOT.radiometry.nL);
%distrib   = hete.distrib;


if isfield(DOT,'grid')
    nx = DOT.grid.Nx;
    ny = DOT.grid.Ny;
    nz = DOT.grid.Nz;
    add = zeros(DOT.grid.N,DOT.radiometry.nL);
    Nx = DOT.grid.Nx;
    Ny = DOT.grid.Ny;
    Nz = DOT.grid.Nz;
    
%-- Regular mesh --%   
else
    M = [DOT.mesh.pos(:,1)-c(1) DOT.mesh.pos(:,2)-c(2) DOT.mesh.pos(:,3)-c(3)];
    add = zeros(DOT.mesh.N,DOT.radiometry.nL);
    Nx = DOT.mesh.Nx;
    Ny = DOT.mesh.Ny;
    Nz = DOT.mesh.Nz;
end

%% here!
smask = load(hete.path);
fn = fieldnames(smask);
mask = smask.(fn{1});
delta = smask.(fn{2});
%% swap fields... in case
if ~isvector(delta)
    dd = mask;
    mask = delta;
    delta = dd;
end

prior = priormask3D(hete.path,DOT.grid);
param = getfield(DOT.opt, var);
for inl = 1:DOT.radiometry.nL
    param(:,:,:,inl) = param(:,:,:,inl) + double(prior).*(intensity(inl)-back(inl));
end
DOT.opt = setfield(DOT.opt, var, param);
end
return;
end



