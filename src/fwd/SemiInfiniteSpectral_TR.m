% TIME-RESOLVED FLUENCE INSIDE A Semi-INFINITE MEDIUM using PCBC 
% function phi = SemiInfinite_TR(time,rs,rd,mua,mus,v,A);
function phi = SemiInfiniteSpectral_TR(time,rs,rd,v,A,S)
%%
a=S.a; b=S.b; lambda=S.lambda; lambda0=S.lambda0; conc=S.conc; ext_c=S.ext_c;
nL = S.nL;
[dim1,dim2]=size(lambda);
if dim1<dim2, lambda=lambda'; S.lambda = lambda; end
% rs source position
% rd detector position
mus = a.*(lambda./lambda0).^(-b);
rs=repmat(rs,[nL 1]);
rd=repmat(rd,[nL 1]);
if any(rs(:,3)) > 0
    z0=rs(:,3);     
 elseif any(rs(:,3)==0)
    z0=1./mus;
    rs(:,3)=z0;
end
D_0 = 1./(3*a.*(lambda./lambda0).^(-b));
ze=2*A*D_0;
rs_min=rs; rs_min(:,3)=-2*ze-z0;

%phi = zeros(size(time));

phi=(InfiniteSpectral_TR(time,rs,rd,v,S)-InfiniteSpectral_TR(time,rs_min,rd,v,S))./(2*A);
phi(isnan(phi))=0;
return