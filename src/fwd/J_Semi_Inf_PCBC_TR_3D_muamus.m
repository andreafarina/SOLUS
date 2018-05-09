function J = J_Semi_Inf_PCBC_TR_3D_muamus(t,mua,mus,cs,A,ri,rj,XX,YY,ZZ,vi,type)
if nargin < 12
    type = 'muaD';
end
switch lower(type)
    case 'mua'
        J = J_Semi_Inf_PCBC_TR_3D_muaD(t,mua,mus,cs,A,ri,rj,XX,YY,ZZ,vi,type);
    case 'mus'
        J = -1./(3*mus^2)*J_Semi_Inf_PCBC_TR_3D_muaD(t,mua,mus,cs,A,ri,rj,XX,YY,ZZ,vi,type);
    case 'muad'
        J = J_Semi_Inf_PCBC_TR_3D_muaD(t,mua,mus,cs,A,ri,rj,XX,YY,ZZ,vi,type);
        nvox = size(J,2)/2;
        J(:,nvox+1:end) = J(:,nvox+1:end) .* (-1./(3*mus^2));
end