
% TIME-RESOLVED FLUENCE INSIDE AN INFINITE MEDIUM 
% function phi = Infinite_TR(time,rs,rd,mua,mus, v,~);
function phi = Infinite_TR(time,rs,rd,mua,mus, v,~,~)
%%
% rs source position
% rd detector position

delta_r=rs-rd;
rhosq=delta_r*delta_r';


%cm/ps velocita' della luce
D = 1/(3*mus);
mu = 1./(4*D*v*time);

%phi = zeros(size(time));

phi=v./(4*pi*D*v*time).^(1.5).*exp(-mua*v*time-mu*rhosq);
phi(isnan(phi))=0;
return