%-------------------------------------------------------------------------%
% Add a spherical inhomogeneity
%
% c         -- [1x3]    -- position of the center of the sphere
% var       -- string   -- the field you wanna mofify
% back      -- [1x1]    -- background value of 'var'
% DISTRIB   -- string   -- indicates wheter you wanna distribute the
%                          intensity within the inhomogeneity ('ON' or
%                          'OFF')
% INTENSITY -- [1x1]    -- maximum value of the inhomogeneous 'var'
% SHAPE     -- string   -- 'gaussian' or 'step'
% SIGMA     -- [1x1]    -- Width of the inhomogeneity
%
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 29/11/10
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 16/12/09
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 22/12/09
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 14/01/09
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 02/02/09
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 17/05/09
%-------------------------------------------------------------------------%


function [DOT,hete] = sphere3D(DOT, hete)
%%--
c    = hete.c;
var  = hete.type;
back = eval(['DOT.opt.',lower(hete.type),'B']);
sigma= hete.sigma;
intensity = hete.val;
shape     = hete.profile;
distrib   = hete.distrib;

%% ----------------------- distances au carree mesh centre ---------------%
%-- Unstructured mesh --%
if isfield(DOT,'grid'),
    [X,Y,Z] = ndgrid(DOT.grid.x, DOT.grid.y, DOT.grid.z);
    X = reshape(X,DOT.grid.N,[]);
    Y = reshape(Y,DOT.grid.N,[]);
    Z = reshape(Z,DOT.grid.N,[]);
    M = [X-c(1) Y-c(2) Z-c(3)];
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

dist2 = sum(M.^2,2);

switch upper(shape)
%% ---------------------- profil de concentration gaussien ---------------%       
case 'GAUSSIAN'        
    %-- selection des indices --%
    indx = find(dist2 < 9*sigma*sigma);
    %-- update concentration --%
    param = getfield(DOT.opt, var);    
    switch upper(distrib)
    case 'ON'    
        add(indx,1) = add(indx,1) + exp(-dist2(indx,1)/sigma/sigma);
        add = add*intensity/sum(add);   
    case 'OFF'
        add(indx,1) = add(indx,1) + (intensity-back)*exp(-dist2(indx,1)/sigma/sigma);
    end   
    
%% --------------------- profil de concentration creneau -----------------%       
case 'STEP'    
    %-- selection des indices --%
    indx = find(dist2 <= sigma*sigma);  
    %-- update concentration --%
    param = getfield(DOT.opt, var);  
    switch upper(distrib)
    case 'ON',      add = intensity/length(indx);    
    case 'OFF',     param(indx) = param(indx) + (intensity-back);
    end
end

%% ---------------------------- Update -----------------------------------%
add = reshape(add,Nx,Ny,Nz);
param(indx) = param(indx) + add(indx);

hete    = setfield(hete, ['d',hete.type], add);
DOT.opt = setfield(DOT.opt, var, param);