function [L,C3D] = StructuredLaplacianPrior(refimage,bx,by,bz)

% create 3D Laplacian, and scale by edge weights.
% refimage is of size bx X by X bz
% at the moment we have to assume bx = by = bz (could be changed)
% return is the decompositon of the Laplacian, 
% i.e. C3D = L'*L
% L is a 3*nvox X nvox sparse matrix where nvox = bx*by*bz. 
% To use in e.g. LSQR solve for A x = b call use
%  x = lsqr([A;lamda*L],[b;nzeros(3*nvox,1],tol,maxit);  
% Or use a symettric solver e.g.
%  x = pcg(A'*A + lambda*C3D,A*b,tol,maxit);
%
nvox = bx*by*bz;
[gx,gy,gz] = gradient(refimage);
gg = sqrt(gx.^2+gy.^2+gz.^2);
kappa = exp(- 5*gg/max(gg(:))); % factor 5 makes quite an extreme setting..
kap3D = spdiags(reshape(kappa,[],1),0:0,nvox,nvox);
kapsqrt = spdiags(reshape(sqrt(kappa),[],1),0:0,nvox,nvox);
% form sparse Laplacian 
% D1d = sparse(bx,bx);
% for i = 1:bx
%     D1d(i,i) = 1;
%     if(i > 1)
%         D1d(i,i-1) = -1;
%     end
% end

D1dz = spdiags( [-ones(bz,1),ones(bz,1)],-1:0,bz,bz);
Dz3D = kron(kron(D1dz,speye(by)),speye(bx));
D1dy = spdiags( [-ones(by,1),ones(by,1)],-1:0,by,by);
Dy3D = kron(kron(speye(bz),D1dy),speye(bx));
D1dx = spdiags( [-ones(bx,1),ones(bx,1)],-1:0,bx,bx);
Dx3D = kron(speye(bz),kron(speye(by),D1dx));

%kapsqrt = 1;
L = [kapsqrt*Dx3D; kapsqrt*Dy3D; kapsqrt*Dz3D];
C3D = L'*L;
