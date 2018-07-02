function j = toastJacobianCW2 (hMesh, hBasis, qvec, mvec,dmask,mua, mus, n,GRADIENT)
% Absolute Jacobian for CW absorption and k
% A. Farina CNR - Polimi	07/2016
% M. Betcke UCL
% possibility to calculate the gradient in 3 ways
% GRADIENT = 
%   'toast'     : use toastImageGradient
%   'operator'  : use gradientOperator function by UCL
%   'matlab'    : use gradient function in Matlab (only 3D)
%   'none'      : avoid calcualation of pdmf_K 

%N=toastMeshNodeCount(hMesh);
if nargin < 9
    GRADIENT = 'operator';
end
%GRADIENT = 'operator';
display(['Jacobian option GRADIENT = ',GRADIENT]);
if strcmpi(GRADIENT,'none')
%    display('J_K part is set to zero');
end
%% Direct and adjoint
%c = 0.299/n(1);
Smat = dotSysmat2_noC(hMesh, mua, mus, n);
%mvec = mvec;
dphi=Smat\qvec;
aphi=Smat\mvec;

BB = hMesh.BoundingBox();
bmin = BB(1,:); bmax = BB(2,:);
bdim = hBasis.Dims;%bdim = bdim';
%hBasis.Map('B->S',ones(bdim'));
mask = (find(hBasis.GridElref>0));

dd = (bmax-bmin)./(bdim-1)';
%se = strel('disk',1);
%bmask = zeros(bdim);
%bmask(mask) = 1;
% for i = 1:bdim(3)
%     bmask2(:,:,i) = imerode(squeeze(bmask(:,:,i)),se);
% end
%bound = bmask-bmask2;
%mask2 = uint32(find(bound(:)==1));
%% gradient operator
if strcmpi(GRADIENT,'operator')
    tmp = zeros(bdim'); tmp(mask) = 1;
    [Dy, Dx, Dz] = gradientOperator(bdim',dd,ones(prod(bdim'),1),'none',tmp);
    clear tmp
end

nsol = hBasis.slen;
nqm = sum(dmask(:));
k=1;
if strcmpi(GRADIENT,'none')
    j = zeros(nqm,nsol);
else
j = zeros(nqm,2*nsol);
end
for q = 1:size(qvec,2)
    ind_m = find(dmask(:,q));	% A. Farina: consistent with dmask order
%   dphi_1=toastMapBasis(hBasis,'M->B',dphi(:,i));
    dphi_1=hBasis.Map('M->B',dphi(:,q));
    switch GRADIENT
        case lower('toast')
            Gdphi = hBasis.ImageGradient (dphi_1);
        case lower('matlab')
            [Gdx,Gdy,Gdz] = gradient(reshape(dphi_1,bdim'),dd(1),dd(2),dd(3));
    end
  for im = 1:length(ind_m)
    tic
    m = ind_m(im);
    %aphi_1=toastMapBasis(hBasis,'M->B',aphi(:,j));
    aphi_1 = hBasis.Map('M->B', aphi(:,m));
    pmdf_mua = -hBasis.Map('B->S',dphi_1.*aphi_1);%./proj(m,q);
    switch GRADIENT
        case lower('toast')
            Gaphi = hBasis.ImageGradient(aphi_1);
            %pmdf_K = - hBasis.Map('B->S',sum(Gdphi .* Gaphi, 1)); %scattering
            pmdf_K = - sum(Gdphi .* Gaphi, 1);
        case lower('matlab')
            [Gadx,Gady,Gadz] = gradient(reshape(aphi_1,bdim'),dd(1),dd(2),dd(3));
            pmdf_K = -(Gdx.*Gadx + Gdy.*Gady + Gdz.*Gadz);
        case lower('operator')
            pmdf_K = - ( (Dy*dphi_1) .* (Dy*aphi_1) ...
            +(Dx*dphi_1) .* (Dx*aphi_1) ...
            +(Dz*dphi_1) .* (Dz*aphi_1) );            
    end
    
    j(k,1:nsol)=pmdf_mua;
    if ~strcmpi(GRADIENT,'none')
        j(k,nsol+(1:nsol))=pmdf_K(mask);%./proj(m,q);
    end
%     mask3 = ismember(mask,mask2);
%     mask3 = [mask3 mask3];
%     j(k,mask3) = 0;
    k=k+1;
  end
end
