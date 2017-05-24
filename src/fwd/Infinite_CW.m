
% CW FLUENCE INSIDE AN INFINITE MEDIUM 
function phi = Infinite_CW(rs,rd,mua,mus)
%%
% rs source position
% rd detector position

delta_r=rs-rd;
rho=sqrt(delta_r*delta_r');
mueff=sqrt(mua*mus*3);


D = 1/(3*mus);




phi=1./(4*pi*D*rho).*exp(-mueff*rho);
%phi(isnan(phi))=0;
return