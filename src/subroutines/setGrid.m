% =========================================================================
% This function sets the reconstruction grid
%
% N. Ducros - Departimento di Fisica - Politecnico di Milano - 24/09/10
% A. Farina - CNR-IFN - Dip. di Fisica - Politecnico di Milano 14/04/15
% A. Farina - CNR-IFN  - Dip di Fisica - Politecnico di Milano 20/12/16
% if hMesh doesn't exist it jump toastBasis
% =========================================================================

function grid = setGrid(DOT)
global mesh
grid = DOT.grid;
%
grid.dV = DOT.grid.dx*DOT.grid.dy*DOT.grid.dz;
%
if (sum(size(mesh)) == 0)
grid.x = (DOT.grid.x1:DOT.grid.dx:(DOT.grid.x2-DOT.grid.dx)) + DOT.grid.dx/2;
grid.y = (DOT.grid.y1:DOT.grid.dy:(DOT.grid.y2-DOT.grid.dy)) + DOT.grid.dy/2;
grid.z = (DOT.grid.z1:DOT.grid.dz:(DOT.grid.z2-DOT.grid.dz)) + DOT.grid.dz/2;
%
grid.Nx = length(grid.x);
grid.Ny = length(grid.y);
grid.Nz = length(grid.z);
grid.N  = grid.Nx*grid.Ny*grid.Nz;
grid.dim = [grid.Nx, grid.Ny, grid.Nz];
%
else %isfield(DOT,'mesh')
    bbox =  mesh.hMesh.BoundingBox;
    minX = min(bbox(:,1));
    minY = min(bbox(:,2));
    minZ = min(bbox(:,3));
    maxX = max(bbox(:,1));
    maxY = max(bbox(:,2));
    maxZ = max(bbox(:,3));
    grid.x = DOT.grid.x1:DOT.grid.dx:DOT.grid.x2;
    grid.y = DOT.grid.y1:DOT.grid.dy:DOT.grid.y2;
    grid.z = DOT.grid.z1:DOT.grid.dz:DOT.grid.z2;
    %
    grid.Nx = length(grid.x);
    grid.Ny = length(grid.y);
    grid.Nz = length(grid.z);
    grid.N  = grid.Nx*grid.Ny*grid.Nz;
    grid.dim=[grid.Nx, grid.Ny, grid.Nz];
    
    
    %------- no bounding box, i.e. the whole mesh volume is considered -------%
    if (DOT.grid.x1 == minX) && (DOT.grid.y1 == minY) && (DOT.grid.z1 == minZ) && ...
            (DOT.grid.x2 == maxX) && (DOT.grid.y2 == maxY) && (DOT.grid.z2 == maxZ)
        
        %    grid.hBasis = toastSetBasis('LINEAR', DOT.mesh.hMesh, ...
        grid.hBasis = toastBasis(mesh.hMesh, ...
            [grid.Nx, grid.Ny, grid.Nz],'LINEAR' );
        %-------- with bounding box, i.e. the mesh is truncated ------------------%
    else
        %    grid.hBasis = toastSetBasis('LINEAR', DOT.mesh.hMesh, ...
        grid.hBasis = toastBasis(mesh.hMesh, ...
            [grid.Nx, grid.Ny, grid.Nz], ...
            [grid.Nx, grid.Ny, grid.Nz], ...
            [DOT.grid.x1, DOT.grid.x2; ...
            DOT.grid.y1,DOT.grid.y2; ...
            DOT.grid.z1,DOT.grid.z2],'LINEAR' );
    end
    %
    grid.mask = find(grid.hBasis.GridElref>0);
    % following is a bit of a guess..
    %grid.mask = grid.hBasis.Map('B->S',ones(grid.Nx, grid.Ny, grid.Nz));
    %grid.NN =length(grid.mask);
end