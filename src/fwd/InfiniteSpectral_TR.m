
% SPECTRAL TIME-RESOLVED FLUENCE INSIDE AN INFINITE MEDIUM
% function phi = InfiniteSpectral_TR(time,rs,rd,mua,mus, v,~);
function phi = InfiniteSpectral_TR(time,rs,rd,v,S)
%%
a=S.a; b=S.b; lambda=S.lambda; lambda0=S.lambda0; conc=S.conc; ext_c=S.ext_c;
nL = S.nL;
[dim1,dim2]=size(lambda);
if dim1>dim2, lambda=lambda'; S.lambda = lambda; end

% rs source position
% rd detector position

delta_r=rs-rd;
%rhosq=delta_r*delta_r';
rhosq=dot(delta_r,delta_r,2);
[dim1,dim2]=size(rhosq);
if dim1>dim2, rhosq=rhosq'; end

%cm/ps velocita' della luce
D_0 = 1./(3*a.*(lambda./lambda0).^(-b));
mu_0 = 1./(4*D_0*v.*time);

%phi = zeros(size(time));

phi=v./(4*pi*D_0*v.*time).^(1.5).*exp(-v*time*ext_c*conc').*exp(-mu_0.*rhosq);
phi(isnan(phi))=0;
return