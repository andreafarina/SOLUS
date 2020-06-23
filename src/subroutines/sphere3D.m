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
for itype = 1:size(hete.type,2)
var  = hete.type{itype}; %modifico il tipo di inclusione Mua in assorbimento e Musp in scattering
back = eval(['DOT.opt.',lower(hete.type{itype}),'B']); %il muaB (e nel ciclo successivo il muspB) per il background calcolato dalle concentrazioni in TranslateXVar
sigma= hete.sigma;
intensity = hete.val((1:DOT.radiometry.nL)+(itype-1)*DOT.radiometry.nL); %val ha gli 8 valori di assorbimento e gli 8 di scattering per l'inclusione calcolati in TranslateXVar
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
    add = zeros(DOT.grid.N,DOT.radiometry.nL); %posso creare questa con 5 colonne per le conc, 1 colonna per a e una per b
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

dist2 = sum(M.^2,2);

switch upper(shape)
%% ---------------------- profil de concentration gaussien ---------------%       
case 'GAUSSIAN'        
    %-- selection des indices --%
    indx = find(dist2 < 9*sigma*sigma); %in valore assoluto tutti gli indici dei numeri di dist2 diversi da zero minori di 9*sigma*sigma
    %-- update concentration --%
    param = getfield(DOT.opt, var);    
    switch upper(distrib)
    case 'ON'    
        add(indx,:) = add(indx,:) + exp(-dist2(indx,1)/sigma/sigma);
        add = add.*intensity./sum(add,1);   
    case 'OFF'
        add(indx,:) = add(indx,:) + (intensity-back).*exp(-dist2(indx,1)/sigma/sigma); 
        %qui si crea l'inclusione togliendo per ogni valore di assorbimento
        %e scattering i valori quelli del background
    end   
    
%% --------------------- profil de concentration creneau -----------------%       
case 'STEP'    
    %-- selection des indices --%
    indx = find(dist2 <= sigma*sigma);  
    %-- update concentration --%
    param = getfield(DOT.opt, var);  
    switch upper(distrib)
    case 'ON',      add = intensity./length(indx);    
    case 'OFF',     
        param = reshape(param,DOT.grid.N,DOT.radiometry.nL);
        param(indx,:) = param(indx,:) + (intensity-back);
        param = reshape(param,Nx,Ny,Nz,DOT.radiometry.nL);
    end
end

%% ---------------------------- Update -----------------------------------%
add = squeeze(reshape(add,Nx,Ny,Nz,DOT.radiometry.nL));
indx_dummy = indx;
for inl = 1:DOT.radiometry.nL-1
    indx = [indx;indx_dummy+inl*DOT.grid.N];
end
param(indx) = param(indx) + add(indx);

hete    = setfield(hete, ['d',hete.type{itype}], add);
DOT.opt = setfield(DOT.opt, var, param);
end