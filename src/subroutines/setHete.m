%==========================================================================
% This version of steHete requires the geometry of the honomogeneity
%
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 15/01/09
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 02/02/09
% N. Ducros - Departamento di Fisica - Politecnico di Milano - 06/07/09
%==========================================================================

function [DOT,hete] = setHete(DOT,hete)

%==========================================================================
%%                           OPTIONS
%==========================================================================
%-- default --%
if not(isfield(hete,'profile')),  hete.profile  = 'gaussian'; end
if not(isfield(hete,'distrib')),  hete.distrib  = 'OFF';      end 
if not(isfield(hete,'geometry')), hete.geometry = 'sphere';   end

%==========================================================================
%%                           MAIN
%==========================================================================
switch upper(hete.geometry)
    
    case 'SPHERE'
    disp('+++ Sphere')
    [DOT, hete] = sphere3D(DOT, hete);
    
    case 'CYLINDER'
    disp('+++ Cylinder')
    [DOT, hete] = cylinder3D(DOT, hete);
    
    case 'USIMAGE'
    disp('+++ Distance Transform');
    hete = hete;
    load(hete.path);
    Mua_r = imresizen(Mua,DOT.grid.dim./size(Mua));
    DOT.opt.Mua = (Mua_r);% - DOT.opt.muaB);%*1e3;
    %[DOT, hete] = priormask3D(DOT, hete);
end