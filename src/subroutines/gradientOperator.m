function [Dy, Dx, Dz, Lap] = gradientOperator(volume, h, c, BC, mask)
% GRADIENT OPERATOR Computes gradient matrices using forward finite difference
%  [Dy, Dx, Dz] = gradientOperator(volume)
%
%  INPUTS:
%   volume   - dimensions of the volume, (1,2 or 3 dimensional)
%   c        - vector of weights for Laplace operator 
%   BC       - boundary conditions {'none', 'Neumann' or 'Dirichlet'}
%    'none'      : do not impose BC. Takes the forward difference 
%                  at each point in the volume. Results in a square block 
%                  and hence matrix and Dy_block(ny,ny) = -1; 
%    'neumann'   : Neumann 0 BC. Removes last row of Dy
%    'dirichlet' : Dirichlet 0 BC. Preappends an extra row [1 0 0 ...]  
%                  Dy_block = [[1 0 0 ...]; Dy_block]
%   mask     - mask for the domain 
%  OUTPUTS:
%   Dy,Dx,Dz - gradients along y,x,z dimension (dimensions as in meshgrid)
%  OPTIONAL OUTPUTS:
%   Lap      - Weighted Laplace operator: 
%               - Dy'*C*Dy - Dx'*C*Dx - Dz'*C*Dz, 
%              where C is diagonal weight matrix C = diag(c).
%              The such constructed Laplacian will has  
%
% Default: does not impose BC. Takes the forward difference 
%          at each point in the volume. Results in a square block 
%          and hence matrix and Dy(ny,ny) = -1; 
%
% Reference: Arridge, Betcke, Harhanen'2014 
%
% Copyright C Marta M. Betcke 2013
% substituted LinIndex with the buit-in sub2ind A. Farina 2020
% Corrected bug on nargin

if nargin < 5
  mask = ones(volume);
end
if nargin < 1
  BC = 'none';
end

% Determine dimension
if volume(end) == 1, volume = volume(1:end-1); end
dim = length(volume);

n = prod(volume);
ny = volume(1);

if dim > 1
  nx = volume(2);
else
  nx = 1; 
end
if dim > 2
  nz = volume(3);
else
  nz = 1;
end

% Compute Dy
% Compute block of Dy
% Dy_block =    -spdiags(ones(ny,       1),   0,     ny,    ny) ...
%              + spdiags(ones(ny,       1),   1,     ny,    ny);
Dy_block =    -1/(2*h(1))*spdiags(ones(ny,       1),   -1,     ny,    ny) ...
             + 1/(2*h(1))*spdiags(ones(ny,       1),    1,     ny,    ny);

% Boundary conditions
switch lower(BC)
  case 'neumann'
    % Remove last row
    %Dy_block = Dy_block(1:ny-1,:);
    Dy_block(ny,:) = 0;
%   case 'dirichlet'
%     % Preappend an extra row
%     extra_row = sparse(1,ny);
%     extra_row(1) = 1;
%     Dy_block = [extra_row; Dy_block];
end
  
% Construct Dy repeating Dy_block nx*nz times 
Dy = kron(speye(nx*nz), Dy_block);


% Compute Dx
if dim > 1
  % Compute block of Dx
%   Dx_block =  -spdiags(ones(ny*nx,    1),   0,     ny*nx, ny*nx) ...
%              + spdiags(ones(ny*nx   , 1),   ny,    ny*nx, ny*nx);  
  Dx_block =  -1/(2*h(2))*spdiags(ones(ny*nx,    1),   -ny,    ny*nx, ny*nx) ...
             + 1/(2*h(2))*spdiags(ones(ny*nx   , 1),    ny,    ny*nx, ny*nx);  
           
  % Boundary conditions
  switch lower(BC)
  case 'neumann'
    % Remove last row ny rows
    %Dx_block = Dx_block(1:(ny*nx-ny),:);
    Dx_block(ny*nx-ny+(1:ny),:) = 0;
%   case 'dirichlet'
%     % Preappend extra ny rows
%     extra_row = spdiags(ones(ny,1), 0, ny, ny*nx);
%     Dx_block = [extra_row; Dx_block];
  end
                    
  % Construct Dx repeating Dx_block nz times
  Dx = kron(speye(nz), Dx_block);
else
  Dx = sparse(n,n); 
end

if dim > 2
  % Compute Dz
%   Dz       =  -spdiags(ones(n,        1),   0,     n,     n) ...
%              + spdiags(ones(n,        1),   ny*nx, n,     n);
  Dz       =  -1/(2*h(3))*spdiags(ones(n,        1),   -ny*nx,     n,     n) ...
             + 1/(2*h(3))*spdiags(ones(n,        1),    ny*nx,     n,     n);
           
  % Boundary conditions
  switch lower(BC)
  case 'neumann'
    % Remove last row ny*nx rows
    %Dz = Dz(1:(n-ny*nx),:);
    Dz(n-ny*nx+(1:ny*nx),:) = 0;
%  case 'dirichlet'
%    % Preappend extra ny*nx rows
%    extra_row = spdiags(ones(ny*nx,1), 0, ny*nx, n);
%    Dz = [extra_row; Dz];
  end

else
  Dz = sparse(n,n); 
end


% Constructing Weighted Laplacian
if nargout > 3
  m = size(Dy,1);
  % Default standard homogenous Laplacian
  if nargin < 3 || isempty(c)
    c = ones(m,1);
  end
  C = spdiags(c,0,m,m);
  Lapy = Dy'*C*Dy;
  Lapx = Dx'*C*Dx;
  Lapz = Dz'*C*Dz;
  
  switch lower(BC)
    case 'dirichlet'
      % Correct the BC in the Laplacian
      for j = 1:ny:n
        % (1,1) entry in each block is 1 instead of 2
        Lapy(j,j) = 2*Lapy(j,j);
      end
      if dim > 1
        for j = kron((0:nz-1)*ny*nx, ones(1,ny)) + kron(ones(1,nz), 1:ny)
          % (j,j) where j = 1,..,ny in each block is 1 istead of 2
          Lapx(j,j) = 2*Lapx(j,j);
        end
      end
      if dim > 2
        for j = 1:ny*nx
          % (j,j) where j = 1,..,ny*nx in each block is 1 istead of 2
          Lapz(j,j) = 2*Lapz(j,j);
        end
      end
  end
  
  Lap = Lapx + Lapy + Lapz;
end


% Find the linear indecies of points at the boundary of the mask (boundary here is the outmost pixel 1)
if dim > 1
  nn = prod(volume(2:dim));
  jy1 = zeros([volume(2:dim),1]);
  jy2 = zeros([volume(2:dim),1]);
  for j = 1:prod(volume(2:dim))
    indy1 =  find(squeeze(mask(:,j)),1,'first');
    if ~isempty(indy1)
      jy1(j) = indy1 + (j-1)*ny;
    else
      jy1(j) = 0; 
    end
    indy2 = find(squeeze(mask(:,j)),1,'last');
    if ~isempty(indy2)
      jy2(j) = indy2 + (j-1)*ny;
    else
      jy2(j) = 0;
    end
  end
else
  jy1 = find(mask(:),1,'first');
  jy2 = find(mask(:),1,'last');
end

if dim == 2
  jx1 = zeros([volume(1),1]);
  jx2 = zeros([volume(1),1]);
  for j = 1:volume(1)
    indx1 = find(reshape(mask(j,:),[],1),1,'first');
    if ~isempty(indx1)
      jx1(j) = sub2ind(volume,j,indx1);
    else
      jx1(j) = 0;
    end
    indx2 = find(reshape(mask(j,:),[],1),1,'last');
    if ~isempty(indx2)
      jx2(j) = sub2ind(volume,j,indx2);
    else
      jx2(j) = 0;
    end
  end
end

if dim > 2
  jx1 = zeros([volume(1),volume(3)]);
  jx2 = zeros([volume(1),volume(3)]);
  for jy = 1:volume(1)
    for jz = 1:volume(3)
      indx1 = find(reshape(mask(jy,:,jz),[],1),1,'first');
      if ~isempty(indx1)
        jx1(jy,jz) = sub2ind(volume,jy,indx1,jz); %pronta da usare
      else
        jx1(jy,jz) = 0;
      end
      indx2 = find(reshape(mask(jy,:,jz),[],1),1,'last');
      if ~isempty(indx2)
        jx2(jy,jz) = sub2ind(volume,jy,indx2,jz);
      else
        jx2(jy,jz) = 0;
      end
    end
  end
end

if dim == 3
  jz1 = zeros(volume(1:2));
  jz2 = zeros(volume(1:2));
  for jy = 1:volume(1)
    for jx = 1:volume(2)
      indz1 = find(reshape(mask(jy,jx,:),[],1),1,'first');
      if ~isempty(indz1)
        jz1(jy,jx) = sub2ind(volume,jy,jx,indz1);
      else
        jz1(jy,jx) = 0;
      end
      indz2 = find(reshape(mask(jy,jx,:),[],1),1,'last');
      if ~isempty(indz2)
        jz2(jy,jx) = sub2ind(volume,jy,jx,indz2);
      else
        jz2(jy,jx) = 0;
      end
    end
  end
end

% Replace the central difference at 
% - the entry to the mask domain with forward difference
% - at the exit of the mask domain with backward difference
% This can be dome just by replacing the corresponding rows
for j = 1:length(jy1(:))
  if jy1(j) ~= jy2(j) && jy1(j) && jy2(j) %if only one pixel in mask domain, leave centrel differences.
    %central -> forward difference
    if jy1(j)-1 > 0
      Dy(jy1(j),jy1(j)-1) = 0; 
    end
    Dy(jy1(j),jy1(j)) = -1/h(1);   %new entry
    Dy(jy1(j),jy1(j)+1) = 1/h(1);  %rescale
    %central -> backward difference
    if jy2(j)+1 < size(Dy,2)
      Dy(jy2(j),jy2(j)+1) = 0; 
    end
    Dy(jy2(j),jy2(j)-1) = -1/h(1); %rescale
    Dy(jy2(j),jy2(j)) = 1/h(1);    %new entry
  end
end
if dim > 1
for j = 1:length(jx1(:))
  if jx1(j) ~= jx2(j) && jx1(j) && jx2(j) %if only one pixel in mask domain, leave centrel differences.
    if jx1(j)-ny > 0 
      Dx(jx1(j),jx1(j)-ny) = 0; %central -> forward difference
    end
    Dx(jx1(j),jx1(j)) = -1/h(2);    %new entry 
    Dx(jx1(j),jx1(j)+ny) = 1/h(2);  %rescale
    
    if jx2(j)+ny <= size(Dx,2)
      Dx(jx2(j),jx2(j)+ny) = 0; %central -> backward difference
    end
    Dx(jx2(j),jx2(j)-ny) = -1/h(2); %rescale
    Dx(jx2(j),jx2(j)) = 1/h(2);     %new entry
   
  end
end
end
if dim > 2
for j = 1:length(jz1(:))
  if jz1(j) ~= jz2(j) && jz1(j) && jz2(j) %if only one pixel in mask domain, leave centrel differences.
    if jz1(j)-ny*nx > 0
      Dz(jz1(j),jz1(j)-ny*nx) = 0; %central -> forward difference
    end
    Dz(jz1(j),jz1(j)) = -1/h(3);       %new entry
    Dz(jz1(j),jz1(j)+ny*nx) = 1/h(3);  %rescale
    
    if jz2(j)+ny*nx < size(Dx,2)
      Dz(jz2(j),jz2(j)+ny*nx) = 0; %central -> backward difference
    end
    Dz(jz2(j),jz2(j)-ny*nx) = -1/h(3); %rescale
    Dz(jz2(j),jz2(j)) = 1/h(3);        %new entry
    
  end
end
end

% if dim >= 2
%   mask = permute(mask, [2 1 3:dim]);
%   for j = 1:volume(1)*prod(volume(3:dim))
%     jx1(j) = find(squeeze(mask(:,j)),1,'first');
%     jx2(j) = find(squeeze(mask(:,j)),1,'last');
%   end
%   mask = permute(mask, [2 1 3:dim]);
%   sub2ind(volume, )
% end
% 
% if dim == 3
%   mask = permute(mask, [3 1 2]);
%   for j = 1:prod(volume(1:2))
%     jz1(j) = find(squeeze(mask(:,j)),1,'first');
%     jz2(j) = find(squeeze(mask(:,j)),1,'last');
%   end
%   mask = permute(mask, [2 3 1]);
% end




switch lower(BC)
  case 'neumann'
    % Eliminate 0-rows, they were only used to compute the Laplacian
    Dy(sum(abs(Dy),2)==0,:) = [];
    Dx(sum(abs(Dx),2)==0,:) = [];
    Dz(sum(abs(Dz),2)==0,:) = [];
end



  
  
