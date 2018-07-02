%-------------------------------------------------------------------------%
% Add a cylindrical inhomogeneity
%
% c         -- [1x3]    -- a point on the axis of the cylinder
% d         -- [1x3]    -- the direction of the axis of the cylinder
% var       -- string   -- the field you wanna mofify
% back      -- [1x1]    -- background value of 'var'
% INTENSITY -- [1x1]    -- maximum value of the inhomogeneous 'var'
% SHAPE     -- string   -- 'gaussian' or 'step'
% SIGMA     -- [1x1]    -- radius of the cylinder
%
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 08/02/10
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 06/07/10
% A. Farina - CNR-IFN - Dip. Fisica  - Politecnico di Milano - 10/04/15
%-------------------------------------------------------------------------%


function [DOT,hete] = cylinder3D(DOT,hete,varargin)

%--
d = hete.d;     % direction vector
d = d./norm(d); % unitary direction vector
c = hete.c;    % first point on the axis
l = hete.l;    % length of the cylinder
for itype = 1:size(hete.type,2)
back = eval(['DOT.opt.',lower(hete.type{itype}),'B']);
sigma = hete.sigma;
intensity = hete.val(itype);
shape = hete.profile;
distrib = hete.distrib;

%-- --%
if isfield(DOT,'grid'),
    %-- distance to the axis of the cylinder --%
    D = [d(1)*ones(DOT.grid.N,1),d(2)*ones(DOT.grid.N,1),d(3)*ones(DOT.grid.N,1)];
    [X,Y,Z] = ndgrid(DOT.grid.x,DOT.grid.y,DOT.grid.z);
    X = reshape(X,DOT.grid.N,[]);
    Y = reshape(Y,DOT.grid.N,[]);
    Z = reshape(Z,DOT.grid.N,[]);
    M = [X-c(1) Y-c(2) Z-c(3)];
    dist = sum(cross(D,M).^2,2).^0.5./norm(d);
    add = zeros(DOT.grid.N,1);
else
    %-- distance to the axis of the cylinder --%
    D = [d(1)*ones(DOT.mesh.N,1),d(2)*ones(DOT.mesh.N,1),d(3)*ones(DOT.mesh.N,1)];
    M = [DOT.mesh.pos(:,1)-c(1) DOT.mesh.pos(:,2)-c(2) DOT.mesh.pos(:,3)-c(3)];
    dist = sum(cross(D,M).^2,2).^0.5./norm(d);
    add = zeros(DOT.mesh.N,1);
end

%-- axial distance --%
L = M*d';
ind1 = find(L>=0);
ind2 = find(L<=l);
ind = intersect(ind1,ind2);

switch upper(shape)
    
%-- profil de concentration gaussien --%       
case 'GAUSSIAN'
    
    %-- selection des indices --%
    indx = find(dist < 3*sigma);
    indx = intersect(indx,ind);
    
    %-- update concentration --%
    if isfield(DOT,'grid'),
        param = getfield(DOT.opt, hete.type{itype});
    else
        a=lower(hete.type{itype});
        param=getfield(DOT.opt,a);
    end
    switch upper(distrib)
    case 'ON'
        add(indx,1) = add(indx,1) + exp(-dist(indx,1).^2/sigma/sigma);
        add = add*intensity/sum(add);
        param(indx) = param(indx) + add(indx);    
    
    case 'OFF'  
        add(indx,1) = add(indx,1) + (intensity-back)*exp(-dist(indx,1).^2/sigma/sigma);
        param(indx) = param(indx) + add(indx);

    end
    
    hete    = setfield(hete, ['d',hete.type{itype}], add);
    if isfield(DOT,'grid'),
       DOT.opt = setfield(DOT.opt, hete.type{itype}, param);
    else
        DOT.opt = setfield(DOT.opt, lower(hete.type{itype}), param);
    end
    
%-- profil de concentration creneau --%       
case 'STEP'
    
    %add = zeros(MOL.mesh.N,1);
    
    %-- selection des indices --%
    indx = find(dist <= sigma);
    indx = intersect(indx,ind);
    
    %-- update concentration --%
    if isfield(DOT,'grid'),
        param = getfield(DOT.opt, hete.type{itype});
    else
        a=lower(hete.type{itype});
        param=getfield(DOT.opt,a);
    end
    switch upper(distrib)
    case 'ON'
        add(indx,1) = intensity/length(indx);
        param(indx) = param(indx) + add(indx);
    
    case 'OFF'
        add(indx,1) = (intensity-back);
        param(indx) = param(indx) + add(indx);
    end
    hete    = setfield(hete, ['d',hete.type{itype}], add);
    if isfield(DOT,'grid'),
       DOT.opt = setfield(DOT.opt, hete.type{itype}, param);
    else
        DOT.opt = setfield(DOT.opt, lower(hete.type{itype}), param);
    end
end
end   