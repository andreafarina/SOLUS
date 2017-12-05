% TIME-RESOLVED FLUENCE INSIDE A Semi-INFINITE MEDIUM using PCBC 
% function phi = SemiInfinite_TR(time,rs,rd,mua,mus,v,A);
function phi = SemiInfinite_TR(time,rs,rd,mua,mus, v,A)
%%
% rs source position
% rd detector position

if rs(3) > 0
    z0=rs(3);     
 elseif rs(3)==0
    z0=1/mus;
    rs(3)=z0;
end
D = 1/(3*mus);
ze=2*A*D;
rs_min=rs; rs_min(3)=-2*ze-z0;

%phi = zeros(size(time));

phi=(Infinite_TR(time,rs,rd,mua,mus,v)-Infinite_TR(time,rs_min,rd,mua,mus,v))./(2*A);
phi(isnan(phi))=0;
return